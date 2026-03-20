import pandas as pd
import matplotlib.pyplot as plt

def parse_logfile(filename):
    """
    Reads a LAMMPS log file and returns a DataFrame of thermo output.
    Works even if log has multiple thermo blocks and summary lines.
    """
    blocks = []
    with open(filename, "r") as f:
        lines = f.readlines()

    header = None
    data_block = []

    for line in lines:
        parts = line.strip().split()
        if not parts:
            continue

        if parts[0].lower() == "step":  # thermo header
            if header is not None and data_block:
                try:
                    df_block = pd.DataFrame(data_block, columns=header, dtype=float)
                    blocks.append(df_block)
                except Exception:
                    pass
            header = parts
            data_block = []
            continue

        if header is not None:
            if parts[0].isdigit():
                data_block.append(parts)
            else:  # end of thermo block
                if data_block:
                    try:
                        df_block = pd.DataFrame(data_block, columns=header, dtype=float)
                        blocks.append(df_block)
                    except Exception:
                        pass
                header = None
                data_block = []

    if header is not None and data_block:
        df_block = pd.DataFrame(data_block, columns=header, dtype=float)
        blocks.append(df_block)

    if not blocks:
        raise ValueError(f"No thermo data found in {filename}")

    df = pd.concat(blocks, ignore_index=True)
    df.columns = [c.lower() for c in df.columns]
    return df

def remove_jump_region(df, start, end):
    """Remove rows where step is between start and end (inclusive)."""
    return df[~df["step"].between(start, end)].copy()

def compare_stats(df_raw, df_clean, label):
    """Make side-by-side statistics for raw vs cleaned data."""
    stats_raw = df_raw.describe().T
    stats_clean = df_clean.describe().T

    comparison = pd.DataFrame({
        "raw_mean": stats_raw["mean"],
        "raw_std": stats_raw["std"],
        "clean_mean": stats_clean["mean"],
        "clean_std": stats_clean["std"],
        "raw_min": stats_raw["min"],
        "clean_min": stats_clean["min"],
        "raw_max": stats_raw["max"],
        "clean_max": stats_clean["max"],
    })
    comparison.to_csv(f"stats_comparison_{label}.csv")
    print(f"\n=== {label} statistics comparison ===")
    print(comparison[["raw_mean", "clean_mean", "raw_std", "clean_std"]])
    return comparison

def plot_relevant(df_raw, df_clean, label):
    """Plot only key quantities: indenter position and load (if present)."""
    if "v_zpos_indent" in df_raw.columns:
        plt.figure(figsize=(8,5))
        plt.plot(df_raw["step"], df_raw["v_zpos_indent"], label="Raw", alpha=0.5)
        plt.plot(df_clean["step"], df_clean["v_zpos_indent"], label="Clean", lw=1.5)
        plt.xlabel("Step")
        plt.ylabel("Indenter Z Position (Å)")
        plt.title(f"Indenter position vs step ({label})")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"indenter_position_{label}.png")
        plt.close()

    if "v_fz" in df_raw.columns:
        plt.figure(figsize=(8,5))
        plt.plot(df_raw["step"], df_raw["v_fz"], label="Raw", alpha=0.5)
        plt.plot(df_clean["step"], df_clean["v_fz"], label="Clean", lw=1.5)
        plt.xlabel("Step")
        plt.ylabel("Force on Indenter (eV/Å)")
        plt.title(f"Indentation load vs step ({label})")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"indentation_force_{label}.png")
        plt.close()

def main():
    print("Reading logs...")

    df300 = parse_logfile("log_300K.lammps")
    df1000 = parse_logfile("log_1000K.lammps")

    df300_clean = remove_jump_region(df300, 100000, 106000)
    df1000_clean = remove_jump_region(df1000, 300000, 311000)

    print("Row counts (raw → clean):")
    print(f"300K: {len(df300)} → {len(df300_clean)}")
    print(f"1000K: {len(df1000)} → {len(df1000_clean)}")

    df300_clean.to_csv("cleaned_300K.csv", index=False)
    df1000_clean.to_csv("cleaned_1000K.csv", index=False)

    compare_stats(df300, df300_clean, "300K")
    compare_stats(df1000, df1000_clean, "1000K")

    plot_relevant(df300, df300_clean, "300K")
    plot_relevant(df1000, df1000_clean, "1000K")

    print("\nOutputs: cleaned CSVs, stats_comparison CSVs, relevant plots.")

if __name__ == "__main__":
    main()
