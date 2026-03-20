#!/usr/bin/env python3
"""
Recompute nanoindentation metrics from LAMMPS log files (thermo output),
recalculate indentation depth = zstart - zpos_indent, detect events,
fit Hertz (early region), compute hardness, work, and save annotated plots.

Place this script in the same directory as:
 - log_300K.lammps
 - log_1000K.lammps
 - dump.indent_300K  (for mapping; optional here)
 - dump.indent_1000K (optional)
Outputs:
 - summary_results_corrected.csv
 - figures/force_depth_300K_fromlog.png, ...1000K...
 - timeseries_fromlog_300K.csv etc.
"""
import re, os, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

plt.rcParams.update({'font.size':11})
OUTDIR = "figures"
os.makedirs(OUTDIR, exist_ok=True)
R_IND=5.0  # Å

# ---------- utility ----------
def parse_lammps_log_for_thermo(logfile):
    """
    Parse a LAMMPS log file for the thermo table printed repeatedly.
    Returns:
      header_cols: list of column names (strings)
      data: numpy array of shape (N, len(header_cols))
      preamble: dict of variables read from header region (zstart_indent etc)
    """
    header_cols = None
    rows = []
    preamble = {}
    with open(logfile,'r') as f:
        lines = f.readlines()

    # First, search for debug print lines with zstart, initial_pillar_height if present
    # simple regex for DEBUG: z_top_atoms = X or "DEBUG: zstart_indent ="
    for L in lines[:200]:
        if 'DEBUG: zstart_indent' in L or 'zstart_indent' in L:
            try:
                v = float(re.findall(r"[-+]?\d*\.\d+|\d+", L)[0])
                preamble['zstart_indent'] = v
            except: pass
        if 'DEBUG: initial_pillar_height' in L or 'initial_pillar_height' in L:
            try:
                v = float(re.findall(r"[-+]?\d*\.\d+|\d+", L)[0]); preamble['initial_pillar_height'] = v
            except: pass
        # also try reading lines like "DEBUG: z_top_atoms = 34.9"
        m = re.search(r"z_top_atoms\s*=\s*([0-9eE\+\-\.]+)", L)
        if m:
            preamble['z_top_atoms']=float(m.group(1))

    # find last thermo header pattern: "Step Temp ..." or lines with column names
    # We'll collect the most frequent header line and use it.
    maybe_headers = []
    for i,L in enumerate(lines):
        if re.match(r'^\s*Step\b', L):
            hdr = L.strip()
            maybe_headers.append((i,hdr))
    if not maybe_headers:
        # sometimes header formatting differs; search for a line with many words and next lines numeric
        for i,L in enumerate(lines):
            tokens=L.strip().split()
            if len(tokens)>=5 and tokens[0].isdigit()==False and all(not t.isdigit() for t in tokens[:2]):
                # check next line if numeric
                if i+1 < len(lines):
                    if re.match(r'^\s*[-+0-9eE\.]', lines[i+1]):
                        maybe_headers.append((i,L.strip()))
    if not maybe_headers:
        raise RuntimeError("No thermo header (Step ...) found in log: " + logfile)

    # pick the last header found (to match final runs)
    hdr_idx, hdr_line = maybe_headers[-1]
    header_cols = hdr_line.split()
    # now read following numeric lines until a non-numeric appears
    j = hdr_idx + 1
    while j < len(lines):
        L = lines[j].strip()
        if L == "" or re.match(r'^[^0-9\-\+\.]', L):
            break
        # check if line is numeric (starts with digit or -)
        if re.match(r'^[\-\d]', L):
            # split and take numeric values, but preserve scientific notation
            toks = L.split()
            # make sure tokens count matches header length; if not, try to skip
            if len(toks) >= len(header_cols):
                rows.append(toks[:len(header_cols)])
            else:
                # if mismatch, just break
                break
        else:
            break
        j += 1

    data = np.array(rows, dtype=float) if rows else np.empty((0,len(header_cols)))
    return header_cols, data, preamble

# ---------- compute ----------
def contact_area(R, d):
    arg = 2.0*R*d - d*d
    arg = np.where(arg>0, arg, 0.0)
    return np.pi * arg  # since arg = a^2

def hertz_fit(F, d, R):
    mask = (d>0) & (~np.isnan(F))
    if mask.sum() < 6:
        return np.nan, np.nan, None
    x = d[mask]**1.5
    y = F[mask]
    slope, intercept, r, p, stderr = stats.linregress(x, y)
    E_r = slope / ((4.0/3.0)*math.sqrt(R))
    return E_r, r**2, (slope, intercept)

def detect_popin(d, F, smooth_window=9):
    if len(F) < smooth_window+2:
        return None, None
    Fs = np.convolve(F, np.ones(smooth_window)/smooth_window, mode='same')
    dz = np.gradient(d)
    dF = np.gradient(Fs)
    with np.errstate(divide='ignore', invalid='ignore'):
        dFdz = np.where(np.abs(dz)>1e-12, dF/dz, np.nan)
    med = np.nanmedian(dFdz); std = np.nanstd(dFdz)
    thr = med - 3*std
    idxs = np.where(dFdz < thr)[0]
    idx = int(idxs[0]) if len(idxs)>0 else int(np.nanargmin(dFdz))
    return idx, dFdz

def analyze_log(prefix):
    logf = f'log_{prefix}.lammps' if os.path.exists(f'log_{prefix}.lammps') else f'log{"" if prefix=="" else "_"+prefix}.lammps'
    if not os.path.exists(logf):
        raise RuntimeError("Log file not found: " + logf)
    hdr, data, pre = parse_lammps_log_for_thermo(logf)
    print("Parsed log:", logf, "columns:", hdr, "rows:", data.shape[0])
    # known variable names to find
    # try to find v_zpos_indent or zpos or zpos_indent
    col_map = {c: i for i,c in enumerate(hdr)}
    # candidate names for indenter zpos and force, pillar height
    zpos_names = [c for c in hdr if 'zpos' in c.lower() or 'v_zpos' in c.lower() or 'z_pos' in c.lower() or 'zpos_indent' in c.lower()]
    fz_names = [c for c in hdr if 'fz' in c.lower() or 'f_' in c.lower() or 'force' in c.lower()]
    height_names = [c for c in hdr if 'pillar' in c.lower() or 'current_pillar_height' in c.lower() or 'height' in c.lower()]
    time_names = [c for c in hdr if 'time' in c.lower() or 'step' in c.lower()]

    # fallback choices
    time_col = time_names[0] if time_names else hdr[0]
    zpos_col = zpos_names[0] if zpos_names else (hdr[hdr.index('v_zpos_indent')] if 'v_zpos_indent' in hdr else None)
    fz_col = fz_names[0] if fz_names else (hdr[hdr.index('v_fz')] if 'v_fz' in hdr else None)
    if zpos_col is None or fz_col is None:
        # try to find columns by substring 'v_zpos' etc robustly
        for c in hdr:
            if 'zpos' in c.lower(): zpos_col=c
            if 'v_zpos' in c.lower(): zpos_col=c
            if 'v_fz' in c.lower(): fz_col=c
            if 'fz' in c.lower() and fz_col is None: fz_col=c

    # convert arrays
    t = data[:, col_map[time_col]]
    zpos = data[:, col_map[zpos_col]] if zpos_col in col_map else None
    fz = data[:, col_map[fz_col]] if fz_col in col_map else None
    # compute zstart from preamble or initial zpos
    if 'zstart_indent' in pre:
        zstart = pre['zstart_indent']
    else:
        zstart = float(zpos[0]) if zpos is not None else None

    # compute indentation depth d = zstart - zpos
    d = zstart - zpos if (zstart is not None and zpos is not None) else None

    # detect Pmax (peak force)
    idx_p = int(np.nanargmax(fz))
    Pmax = float(fz[idx_p])
    depth_at_pmax = float(d[idx_p])
    time_at_pmax = float(t[idx_p])

    # detect pop-in
    event_idx, dFdz = detect_popin(d, fz)
    event_time = float(t[event_idx])
    event_depth = float(d[event_idx])

    # Hertz
    E_r, R2, fitparams = hertz_fit(fz, d, R_IND)
    a_p = contact_area(R_IND, depth_at_pmax)
    hardness = Pmax / a_p if a_p>0 else np.nan
    work = np.trapz(fz[np.argsort(d)], np.sort(d))

    # prepare dataframe for timeseries
    df_ts = pd.DataFrame({'time_ps': t, 'zpos_A': zpos, 'depth_A': d, 'F_eV_per_A': fz})
    outcsv = os.path.join(OUTDIR, f"timeseries_fromlog_{prefix}.csv")
    df_ts.to_csv(outcsv, index=False)
    # plot force vs depth and mark event
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7,5))
    ax.plot(d, fz, '-', lw=1)
    ax.set_xlabel("Indentation depth d (Å)")
    ax.set_ylabel("Force F (eV/Å)")
    ax.set_title(f"Force vs depth ({prefix})")
    ax.plot(depth_at_pmax, Pmax, 'ro', label=f"Pmax={Pmax:.2f}")
    ax.axvline(event_depth, color='k', linestyle='--', label=f"event @ {int(event_time)} ps")
    ax.legend()
    # inset dF/dz
    left, bottom, width, height = 0.55, 0.55, 0.35, 0.33
    axins = fig.add_axes([left, bottom, width, height])
    axins.plot(d, dFdz, '-', lw=0.8)
    axins.axhline(0, color='k', lw=0.6)
    axins.set_xlabel('d (Å)', fontsize=9); axins.set_ylabel('dF/dz', fontsize=9)
    plt.tight_layout()
    outpng = os.path.join(OUTDIR, f"force_depth_fromlog_{prefix}.png")
    fig.savefig(outpng, dpi=200)
    plt.close(fig)

    return {
        'tag': prefix,
        'Pmax_eV_per_A': Pmax,
        'depth_at_Pmax_A': depth_at_pmax,
        'time_at_Pmax_ps': time_at_pmax,
        'event_idx': int(event_idx),
        'event_time_ps': event_time,
        'event_depth_A': event_depth,
        'E_r_eV_per_A3': E_r,
        'E_r_R2': R2,
        'hardness_eV_per_A3': hardness,
        'work_eV': work,
        'timeseries_csv': outcsv,
        'figure_png': outpng
    }

# ---------- MAIN ----------
if __name__ == "__main__":
    reslist = []
    for prefix in ['300K','1000K']:
        try:
            res = analyze_log(prefix)
            reslist.append(res)
            print("Done", prefix, "Pmax:", res['Pmax_eV_per_A'], "depth@Pmax:", res['depth_at_Pmax_A'], "event_time:", res['event_time_ps'])
        except Exception as e:
            print("Error processing", prefix, e)

    if reslist:
        df = pd.DataFrame(reslist)
        df.to_csv("summary_results_corrected.csv", index=False)
        print("Wrote summary_results_corrected.csv and plots in", OUTDIR)
    else:
        print("No results computed.")
