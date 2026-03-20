# file: analyze_nanoindent.py
import os, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress
from dataclasses import dataclass
from datetime import datetime

# ----------------------------- CONFIG ---------------------------------
INDENTER_RADIUS_A = 5.0          # Å
INDENTER_SPEED_A_PS = 0.1        # Å/ps
EARLY_FRACTION_FOR_ELASTIC = 0.10  # first 10% for elastic/Hertz fit
THERMO_COLUMNS = [
    # must match your thermo_style custom (see in.indent.*)
    # Step Temp Press PE KE Etotal posz fz pillar_height pillar_volume
    # stress_zz time_ps dz strain temp_avg
    "step","temp","press","pe","ke","etotal",
    "posz","fz","pillar_height","pillar_volume",
    "stress_zz","time_ps","dz","strain","temp_avg"
]
# ----------------------------------------------------------------------

@dataclass
class CaseResult:
    temp_label: str
    max_force: float = np.nan
    force_std: float = np.nan
    initial_pe: float = np.nan
    final_pe: float = np.nan
    energy_change: float = np.nan
    avg_ke: float = np.nan
    etotal_std: float = np.nan
    max_depth: float = np.nan
    initial_height: float = np.nan
    final_height: float = np.nan
    height_change: float = np.nan
    max_stress: float = np.nan
    avg_pressure: float = np.nan
    pressure_std: float = np.nan
    avg_temp: float = np.nan
    temp_std: float = np.nan
    total_work: float = np.nan
    elastic_modulus_lin: float = np.nan
    strain_rate_eff: float = np.nan
    hertz_Er: float = np.nan
    hertz_fit_R2: float = np.nan
    stiffness_max: float = np.nan
    hardness_at_Pmax: float = np.nan
    popin_indices: list = None

def parse_lammps_log(filepath: str) -> pd.DataFrame | None:
    if not os.path.exists(filepath):
        print(f"Error: {filepath} not found")
        return None

    data, collecting = [], False
    with open(filepath, "r") as f:
        for raw in f:
            line = raw.strip()
            if line.startswith("Step"):
                collecting = True
                continue
            if collecting:
                parts = line.split()
                if len(parts) == len(THERMO_COLUMNS):
                    # verify numeric
                    try:
                        [float(p) for p in parts]
                        data.append(parts)
                    except ValueError:
                        pass
                else:
                    # end of this thermo block
                    if data:
                        break
    if not data:
        print(f"Warning: no thermo rows parsed from {filepath}")
        return None

    df = pd.DataFrame(data, columns=THERMO_COLUMNS)
    df = df.apply(pd.to_numeric, errors="coerce").dropna()
    if df.empty:
        print(f"Warning: thermo became empty after numeric cleaning for {filepath}")
        return None
    return df

def trapezoid_work(df):
    try:
        return float(np.trapz(df["fz"].values, df["dz"].values))
    except Exception:
        return np.nan

def elastic_modulus_linear(df):
    # linear slope stress vs strain in early regime
    n = len(df)
    if n < 10:
        return np.nan
    k = max(5, int(n * EARLY_FRACTION_FOR_ELASTIC))
    x = df["strain"].values[:k]
    y = df["stress_zz"].values[:k]
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return np.nan
    m, _, r, _, _ = linregress(x[mask], y[mask])
    return float(m)

def strain_rate(df):
    t = df["time_ps"].values
    e = df["strain"].values
    dt = t.max() - t.min()
    if dt <= 0:
        return np.nan
    return float((e.max() - e.min()) / dt)

def stiffness_curve(df):
    # numerical derivative dF/dz
    z = df["dz"].values
    f = df["fz"].values
    if len(z) < 3:
        return np.nan, None
    # guard monotonicity for dz (should increase)
    order = np.argsort(z)
    z, f = z[order], f[order]
    df_dz = np.gradient(f, z, edge_order=2)
    return float(np.nanmax(df_dz)), pd.Series(df_dz, index=df.index[order])

def detect_popins(df, threshold_sigma=4.0):
    # sudden displacement bursts: look at incremental dz/step spikes
    if len(df) < 5:
        return []
    dz = df["dz"].values
    ddz = np.diff(dz, prepend=dz[0])
    mu, sd = np.mean(ddz), np.std(ddz)
    if sd == 0:
        return []
    spikes = np.where((ddz - mu) > threshold_sigma * sd)[0].tolist()
    return spikes

def hertz_fit_spherical(df, R_angstrom=INDENTER_RADIUS_A, fraction=EARLY_FRACTION_FOR_ELASTIC):
    """
    Fit F = (4/3) * E_r * sqrt(R) * delta^(3/2)  (pure elastic, spherical indenter)
    Use early portion to avoid plasticity.
    Returns (E_r, R^2)
    Units: F in eV/Å, delta in Å, R in Å -> E_r in eV/Å^3 (same dimension as stress)
    """
    n = len(df)
    if n < 12:
        return np.nan, np.nan
    k = max(8, int(n * fraction))
    delta = df["dz"].values[:k]
    F = df["fz"].values[:k]
    # Only fit where both are positive
    mask = (delta > 0) & (F > 0) & np.isfinite(delta) & np.isfinite(F)
    if mask.sum() < 6:
        return np.nan, np.nan
    x = delta[mask] ** 1.5
    y = F[mask]
    # linear fit y = a * x, where a = (4/3) E_r sqrt(R)
    a, _, r, _, _ = linregress(x, y)
    if not np.isfinite(a) or a <= 0:
        return np.nan, np.nan
    Er = (3.0 / 4.0) * a / np.sqrt(R_angstrom)
    return float(Er), float(r**2)

def hardness_at_pmax(df):
    # H = Pmax / A_contact at Pmax (your true_stress_zz uses same idea)
    if "contact_area" in df.columns:
        # If your thermo had it, great. (Current thermo prints it as v_contact_area)
        # But your posted thermo_style doesn't include it. We'll reconstruct it to be safe.
        pass
    # reconstruct Hertz contact area for sphere: a = sqrt(2 R delta - delta^2), A = pi a^2
    R = INDENTER_RADIUS_A
    delta = df["dz"].values
    a_sq = np.maximum(0.0, 2.0*R*delta - delta*delta)
    A = np.pi * a_sq + 1e-12
    F = df["fz"].values
    imax = int(np.nanargmax(F))
    H = F[imax] / A[imax]  # units: (eV/Å) / Å^2 = eV/Å^3 (stress-like)
    return float(H)

def analyze_case(df: pd.DataFrame, label: str) -> CaseResult:
    res = CaseResult(temp_label=label)

    res.max_force = float(df["fz"].max())
    res.force_std = float(df["fz"].std())
    res.initial_pe = float(df["pe"].iloc[0])
    res.final_pe = float(df["pe"].iloc[-1])
    res.energy_change = res.final_pe - res.initial_pe
    res.avg_ke = float(df["ke"].mean())
    res.etotal_std = float(df["etotal"].std())
    res.max_depth = float(df["dz"].max())
    res.initial_height = float(df["pillar_height"].iloc[0])
    res.final_height = float(df["pillar_height"].iloc[-1])
    res.height_change = res.initial_height - res.final_height
    res.max_stress = float(df["stress_zz"].max())
    res.avg_pressure = float(df["press"].mean())
    res.pressure_std = float(df["press"].std())
    res.avg_temp = float(df["temp"].mean())
    res.temp_std = float(df["temp"].std())
    res.total_work = trapezoid_work(df)
    res.elastic_modulus_lin = elastic_modulus_linear(df)
    res.strain_rate_eff = strain_rate(df)
    res.hertz_Er, res.hertz_fit_R2 = hertz_fit_spherical(df)
    res.stiffness_max, dfdz = stiffness_curve(df)
    res.hardness_at_Pmax = hardness_at_pmax(df)
    res.popin_indices = detect_popins(df)

    return res

def load_and_analyze(log_300k, log_1000k):
    out = {}
    df300 = parse_lammps_log(log_300k)
    df1000 = parse_lammps_log(log_1000k)

    if df300 is not None:
        out["300K"] = analyze_case(df300, "300K")
    if df1000 is not None:
        out["1000K"] = analyze_case(df1000, "1000K")
    return out, df300, df1000

def plots(df300, df1000, outdir="analysis_plots"):
    os.makedirs(outdir, exist_ok=True)
    pairs = [
        ("dz","fz","Depth (Å)","Force (eV/Å)","force_depth"),
        ("strain","stress_zz","Strain","Stress (eV/Å³)","stress_strain"),
        ("time_ps","temp","Time (ps)","Temperature (K)","temperature"),
        ("time_ps","pe","Time (ps)","Potential Energy (eV)","pe"),
        ("time_ps","pillar_height","Time (ps)","Pillar Height (Å)","height"),
        ("time_ps","press","Time (ps)","Pressure (bars)","pressure"),
        ("time_ps","ke","Time (ps)","Kinetic Energy (eV)","ke"),
        ("time_ps","etotal","Time (ps)","Total Energy (eV)","etotal"),
    ]
    for x,y,xlab,ylab,name in pairs:
        plt.figure(figsize=(10,6))
        plotted=False
        if df300 is not None and x in df300 and y in df300:
            plt.plot(df300[x], df300[y], label="300 K", linewidth=2)
            plotted=True
        if df1000 is not None and x in df1000 and y in df1000:
            plt.plot(df1000[x], df1000[y], label="1000 K", linewidth=2)
            plotted=True
        if plotted:
            plt.xlabel(xlab); plt.ylabel(ylab); plt.legend(); plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(outdir, f"{name}.png"), dpi=300)
        plt.close()

def write_report(results: dict, df300, df1000, path_txt="nanoindentation_analysis.txt", path_csv="nanoindentation_summary.csv"):
    with open(path_txt,"w") as f:
        f.write("=== Nanoindentation Analysis Report ===\n")
        f.write(f"Generated: {datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')} UTC\n")
        f.write("Indenter: spherical, R = 5.0 Å; speed = 0.1 Å/ps; units = metal\n\n")
        for label,res in results.items():
            f.write(f"--- {label} ---\n")
            f.write(f"Max Force (eV/Å): {res.max_force:.4f}\n")
            f.write(f"Hardness @ Pmax (eV/Å³): {res.hardness_at_Pmax:.4e}\n")
            f.write(f"Max Depth (Å): {res.max_depth:.4f}\n")
            f.write(f"Height Change (Å): {res.height_change:.4f}\n")
            f.write(f"Total Work (eV): {res.total_work:.4f}\n")
            f.write(f"Elastic slope (stress/strain, early) (eV/Å³): {res.elastic_modulus_lin:.4e}\n")
            f.write(f"Hertz reduced modulus E_r (eV/Å³) [R²={res.hertz_fit_R2:.3f}]: {res.hertz_Er:.4e}\n")
            f.write(f"Max stiffness dF/dz (eV/Å²): {res.stiffness_max:.4e}\n")
            f.write(f"Max stress_zz (eV/Å³): {res.max_stress:.4e}\n")
            f.write(f"Avg Temp ± std (K): {res.avg_temp:.2f} ± {res.temp_std:.2f}\n")
            if res.popin_indices:
                f.write(f"Pop-in candidates (indices): {res.popin_indices[:10]}")
                if len(res.popin_indices) > 10: f.write(" ...")
                f.write("\n")
            f.write("\n")
    # CSV summary
    rows=[]
    for label,res in results.items():
        rows.append({
            "case": label,
            "max_force_eV_per_A": res.max_force,
            "hardness_at_Pmax_eV_per_A3": res.hardness_at_Pmax,
            "max_depth_A": res.max_depth,
            "height_change_A": res.height_change,
            "total_work_eV": res.total_work,
            "elastic_modulus_lin_eV_per_A3": res.elastic_modulus_lin,
            "hertz_Er_eV_per_A3": res.hertz_Er,
            "hertz_R2": res.hertz_fit_R2,
            "stiffness_max_eV_per_A2": res.stiffness_max,
            "max_stress_eV_per_A3": res.max_stress,
            "avg_temp_K": res.avg_temp,
            "temp_std_K": res.temp_std,
            "strain_rate_ps_inv": res.strain_rate_eff,
            "avg_pressure_bar": res.avg_pressure,
            "pressure_std_bar": res.pressure_std,
            "etotal_std_eV": res.etotal_std
        })
    pd.DataFrame(rows).to_csv(path_csv, index=False)
    print(f"✓ wrote {path_txt} and {path_csv}")

def main():
    results, df300, df1000 = load_and_analyze("log_300K.lammps","log_1000K.lammps")
    if not results:
        print("No results.")
        return
    plots(df300, df1000)
    write_report(results, df300, df1000)
    # quick console summary
    for k,v in results.items():
        print(f"{k}: Pmax={v.max_force:.2f} eV/Å, depth={v.max_depth:.2f} Å, E_r(hertz)={v.hertz_Er:.2e} eV/Å³")

if __name__ == "__main__":
    main()
