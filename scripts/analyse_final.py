#!/usr/bin/env python3
"""
analyze_nanoindent.py

Reads LAMMPS-derived analysis files produced in your working directory,
computes the requested metrics, detects pop-in events, and produces
annotated plots and a CSV summary.

Expected files (from your listing) - adjust paths if needed:
 - force_evolution_300K.dat
 - force_evolution_1000K.dat
 - dump.indent_300K_with_indenter (if needed for OVITO; not parsed here)
 - temp_profile_300K.dat, temp_profile_1000K.dat (optional)
 - temp_profile_*.dat (any)
 - nanoindentation_summary.csv (optional summary)
 - nanoindentation_analysis.txt (optional)
 - crack_analysis_300K_t500000.png, crack_analysis_1000K_t500000.png (used as figures)
 - files in analysis_plots/ such as structure_evolution.png
Output:
 - annotated PNG files saved to ./figures/
 - summary_results.csv with computed metrics
"""

import os, re, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit

plt.rcParams.update({'font.size': 12})

# -------------------------
# User-editable settings
# -------------------------
R_INDENTER = 5.0            # Angstrom
EARLY_FRACTION = 0.15       # fraction of maximum depth used for Hertz fit
DERIV_SMOOTH = 5            # points for smoothing derivative (moving average)
EVENT_DF_THRESHOLD = -0.5   # fallback threshold for dF/dz detection (adjust as needed)
OUTDIR = "figures"
os.makedirs(OUTDIR, exist_ok=True)

# -------------------------
# Utility / physics funcs
# -------------------------
EVA3_TO_GPA = 160.21766208  # 1 eV/Å^3 -> GPa

def eV_per_A3_to_GPa(x):
    return x * EVA3_TO_GPA

def contact_radius(R, d):
    """Contact radius a = sqrt(2 R d - d^2), guard for d>2R"""
    arg = 2.0 * R * d - d * d
    arg = np.where(arg > 0, arg, 0.0)
    return np.sqrt(arg)

def contact_area(R, d):
    a = contact_radius(R, d)
    return np.pi * a**2

def hertz_fit(F, d, R):
    """
    Fit F = (4/3) * E_r * sqrt(R) * d^(3/2) for small d.
    Returns E_r and R^2 for the fit.
    """
    x = d**1.5
    # linear fit F = C * x -> C = (4/3) * E_r * sqrt(R)
    slope, intercept, r_value, p_value, stderr = stats.linregress(x, F)
    C = slope
    E_r = C / ((4.0/3.0) * math.sqrt(R))
    return E_r, r_value**2, slope, intercept

# -------------------------
# Parsing helpers
# -------------------------
def try_read_space_table(fname, colnames=None, skiprows=0):
    """Try to read a whitespace separated file; return DataFrame or None"""
    if not os.path.exists(fname):
        return None
    try:
        df = pd.read_csv(fname, delim_whitespace=True, comment='#', header=None, skiprows=skiprows)
        if colnames is not None and df.shape[1] >= len(colnames):
            df = df.iloc[:, :len(colnames)]
            df.columns = colnames
        return df
    except Exception as e:
        print(f"[warn] couldn't parse {fname}: {e}")
        return None

def parse_force_evolution(fname):
    """
    Expect file format with columns: time_ps, fz, true_stress_zz, contact_area (or similar).
    We will be flexible and infer columns.
    """
    df = try_read_space_table(fname)
    if df is None:
        return None
    # heuristics to name columns if unknown
    ncol = df.shape[1]
    if ncol >= 4:
        df.columns = ['time_ps', 'fz', 'true_stress_zz', 'contact_area'] + list(df.columns[4:])
    elif ncol == 3:
        df.columns = ['time_ps', 'fz', 'contact_area']
    elif ncol == 2:
        df.columns = ['time_ps', 'fz']
    else:
        df.columns = [f'c{i}' for i in range(ncol)]
    # make sure numeric
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    return df

# -------------------------
# Core analysis
# -------------------------
def detect_popin_by_derivative(z, F, time=None):
    """
    Compute derivative dF/dz and detect sharp negative spikes (pop-in).
    Return event index and derivative array.
    """
    # smooth F to reduce noise
    kernel = np.ones(DERIV_SMOOTH)/DERIV_SMOOTH
    Fs = np.convolve(F, kernel, mode='same')
    # derivative w.r.t z (avoid zero division)
    dz = np.gradient(z)
    dF = np.gradient(Fs)
    with np.errstate(divide='ignore', invalid='ignore'):
        dFdz = np.where(dz != 0, dF / dz, np.nan)
    # pick candidate index as most negative dFdz within top percent of |dFdz|
    # fallback threshold
    idx = np.nanargmin(dFdz)
    minval = dFdz[idx]
    # find first index where derivative drops below some threshold relative to distribution
    med = np.nanmedian(dFdz)
    std = np.nanstd(dFdz)
    thr = med - max(3*std, abs(EVENT_DF_THRESHOLD))
    candidates = np.where(dFdz < thr)[0]
    if candidates.size:
        idx = int(candidates[0])
    return idx, dFdz

def integrate_work(z, F):
    """Simple trapezoidal integration for work (area under F(z))."""
    # Ensure monotonic z; otherwise sort by z
    order = np.argsort(z)
    zsorted = z[order]
    Fsorted = F[order]
    work = np.trapz(Fsorted, zsorted)
    return work

def parse_summary_txt(fname="nanoindentation_analysis.txt"):
    """Try to extract some key numbers from the text file you generated."""
    out = {}
    if not os.path.exists(fname):
        return out
    txt = open(fname).read()
    # regex patterns:
    patterns = {
        'Pmax_300K': r"300K: Pmax=([0-9Ee.+-]+)",
        'Pmax_1000K': r"1000K: Pmax=([0-9Ee.+-]+)",
        'E_r_300K': r"300K: .*E_r\([^\)]*\)=([0-9Ee.+-]+)",
        'E_r_1000K': r"1000K: .*E_r\([^\)]*\)=([0-9Ee.+-]+)"
    }
    for k,p in patterns.items():
        m = re.search(p, txt)
        if m:
            out[k] = float(m.group(1))
    return out

def parse_crack_analysis_text(fname="nanoindentation_analysis.txt"):
    """If a fuller txt output contains crack numbers, extract them (fallback to patterns found earlier)."""
    # We also look at crack_analysis.py printed outputs saved in working dir if you have a textual log.
    if not os.path.exists("crack_analysis.py") and not os.path.exists("crack_analysis_300K_t500000.png"):
        pass
    # fallback: the script cannot reliably parse graphical PNGs. Use the summary CSV if present.
    return {}

# -------------------------
# Main pipeline
# -------------------------
def main():
    # Targets / files
    files = {
        'force_300': 'force_evolution_300K.dat',
        'force_1000': 'force_evolution_1000K.dat',
        'temp_300': 'temp_profile_300K.dat',
        'temp_1000': 'temp_profile_1000K.dat',
        'summary_txt': 'nanoindentation_analysis.txt',
        'summary_csv': 'nanoindentation_summary.csv'
    }

    # read summary TXT for quick numbers
    quick = parse_summary_txt(files['summary_txt'])

    results = []

    for tag, fname in [('300K', files['force_300']), ('1000K', files['force_1000'])]:
        df = parse_force_evolution(fname)
        if df is None:
            print(f"[warn] {fname} not found or not parsed. Skipping {tag}.")
            continue

        # decide which columns are which
        # prefer columns named 'time_ps','fz','contact_area'
        time = df.columns[0]
        time_arr = df[time].to_numpy()
        # heuristics for F and z and contact area / height columns
        # you used thermo to print v_zpos_indent (v_zpos_indent) and v_current_pillar_height
        # common names: v_zpos_indent, v_fz, v_current_pillar_height
        col_candidates = df.columns.tolist()
        # guess F
        Fcol = [c for c in col_candidates if 'fz' in c.lower() or 'f_' in c.lower() or c.lower()=='fz']
        if not Fcol:
            # try second column
            Fcol = [col_candidates[1]]
        Fcol = Fcol[0]
        # guess zpos indent
        zcol = [c for c in col_candidates if 'zpos' in c.lower() or 'zpos_indent' in c.lower() or 'v_zpos_indent' in c.lower()]
        if not zcol:
            zcol = [c for c in col_candidates if 'zpos' in c.lower()]
        if not zcol:
            # fallback to a column with numbers in expected range (10-100)
            zcol = [col_candidates[2]] if len(col_candidates) > 2 else col_candidates[0]
        zcol = zcol[0]

        # Also try to find pillar height column
        hcol = [c for c in col_candidates if 'height' in c.lower() or 'current_pillar_height' in c.lower() or 'v_current_pillar_height' in c.lower()]
        hcol = hcol[0] if hcol else None

        t = df[time].to_numpy()
        F = pd.to_numeric(df[Fcol], errors='coerce').to_numpy()
        z = pd.to_numeric(df[zcol], errors='coerce').to_numpy()
        if hcol:
            height = pd.to_numeric(df[hcol], errors='coerce').to_numpy()
        else:
            height = None

        # compute simple metrics
        valid = (~np.isnan(F)) & (~np.isnan(z))
        if valid.sum() < 10:
            print(f"[warn] not enough valid F vs z data for {tag}.")
            continue
        t = t[valid]; F = F[valid]; z = z[valid]

        # peak force
        idx_pmax = np.nanargmax(F)
        Pmax = float(F[idx_pmax])
        depth_at_Pmax = float(z[idx_pmax])

        # detect pop-in via derivative dF/dz
        event_idx, dFdz = detect_popin_by_derivative(z, F, time=t)
        event_step = int(t[event_idx]) if len(t)>event_idx else None

        # Hertz fit on early region (d < EARLY_FRACTION * max depth)
        maxd = np.nanmax(z)
        early_mask = (z >= 0) & (z <= max(EARLY_FRACTION * maxd, 0.5))
        try:
            E_r, R2, slope, intercept = hertz_fit(F[early_mask], z[early_mask], R_INDENTER)
        except Exception as e:
            E_r, R2, slope, intercept = np.nan, np.nan, np.nan, np.nan

        # dF/dz early slope (elastic stiffness) via linear fit F vs z in early region
        try:
            slope_elastic, intercept_e, r_e, p_e, stderr_e = stats.linregress(z[early_mask], F[early_mask])
        except:
            slope_elastic = np.nan

        # contact area at Pmax (compute from R and depth)
        a = contact_area(R_INDENTER, depth_at_Pmax)
        hardness = Pmax / a if a>0 else np.nan

        # integrate work (area under loading curve)
        work_total = integrate_work(z, F)

        # prepare results
        res = {
            'tag': tag,
            'Pmax_eV_per_A': Pmax,
            'depth_at_Pmax_A': depth_at_Pmax,
            'event_time_ps': event_step,
            'event_index': int(event_idx),
            'E_r_eV_per_A3': E_r,
            'E_r_R2': R2,
            'early_stiffness_dFdz_eV_per_A2': slope_elastic,
            'contact_area_A2': a,
            'hardness_eV_per_A3': hardness,
            'work_total_eV': work_total,
        }

        # if height available record initial and final
        if height is not None:
            res['initial_height_A'] = float(np.nanmedian(height[:10]))
            res['final_height_A'] = float(np.nanmedian(height[-20:]))
            res['height_change_A'] = res['final_height_A'] - res['initial_height_A']

        # save arrays for plotting
        res['t'] = t; res['F'] = F; res['z'] = z; res['dFdz'] = dFdz
        results.append(res)

        # produce annotated plots for this run
        # Force vs depth with pop-in vertical line and derivative inset
        fig, ax = plt.subplots(figsize=(7,5))
        ax.plot(z, F, '-', lw=1.0)
        ax.set_xlabel("Indentation depth z (Å)")
        ax.set_ylabel("Force F (eV/Å)")
        ax.set_title(f"Force vs Depth ({tag})")
        # mark Pmax
        ax.plot(depth_at_Pmax, Pmax, 'ro', label=f"P_max={Pmax:.2f}")
        ax.axvline(depth_at_Pmax, color='r', alpha=0.3)
        # mark pop-in by mapping event index z[event_idx]
        if event_idx is not None and 0 <= event_idx < len(z):
            z_ev = z[event_idx]; t_ev = t[event_idx]
            ax.axvline(z_ev, color='k', linestyle='--', lw=1.2)
            ax.text(z_ev, ax.get_ylim()[1]*0.8, f"event@t={t_ev:.0f} ps", rotation=90, va='top', ha='right', bbox=dict(alpha=0.6))
        ax.legend(loc='best')
        # inset derivative
        left, bottom, width, height = 0.55, 0.55, 0.35, 0.33
        axins = fig.add_axes([left, bottom, width, height])
        # plot dF/dz vs z (smoothed)
        # compute dF/dz
        axins.plot(z, res['dFdz'], '-', lw=0.8)
        axins.axhline(0, color='k', lw=0.5)
        axins.set_xlabel("z (Å)", fontsize=9)
        axins.set_ylabel("dF/dz", fontsize=9)
        axins.tick_params(axis='both', which='major', labelsize=8)
        plt.tight_layout()
        outpng = os.path.join(OUTDIR, f"force_depth_{tag}.png")
        plt.savefig(outpng, dpi=200)
        plt.close(fig)
        print(f"[info] saved {outpng}")

        # Indenter position vs time (z vs t) with event marker
        fig, ax = plt.subplots(figsize=(7,4))
        ax.plot(t, z, '-', lw=0.8)
        ax.set_xlabel("Time (ps)")
        ax.set_ylabel("Indenter z (Å)")
        ax.set_title(f"Indenter z vs time ({tag})")
        if event_idx is not None and 0 <= event_idx < len(t):
            ax.axvline(t[event_idx], color='r', linestyle='--')
            ax.text(t[event_idx], np.nanmedian(z), f" event@{int(t[event_idx])} ps", color='r')
        plt.tight_layout()
        outpng = os.path.join(OUTDIR, f"indenter_position_{tag}.png")
        plt.savefig(outpng, dpi=200)
        plt.close(fig)
        print(f"[info] saved {outpng}")

        # Save force/time/depth arrays to CSV for quick viewing
        arrdf = pd.DataFrame({'time_ps': t, 'z_A': z, 'F_eV_per_A': F, 'dFdz': res['dFdz']})
        arrdf.to_csv(os.path.join(OUTDIR, f"timeseries_{tag}.csv"), index=False)

    # Summarize results to CSV
    if results:
        rows = []
        for r in results:
            row = {k: v for k,v in r.items() if not isinstance(v, np.ndarray) and k not in ('t','F','z','dFdz')}
            rows.append(row)
        ssum = pd.DataFrame(rows)
        ssum.to_csv("summary_results.csv", index=False)
        print("[info] wrote summary_results.csv")
    else:
        print("[warn] no results to summarize.")

    # Try to parse static precomputed summary CSV if present and merge (optional)
    if os.path.exists(files['summary_csv']):
        try:
            pre = pd.read_csv(files['summary_csv'])
            pre.to_csv(os.path.join(OUTDIR, "precomputed_summary.csv"), index=False)
            print("[info] copied existing precomputed summary to figures/")
        except:
            pass

if __name__ == "__main__":
    main()
