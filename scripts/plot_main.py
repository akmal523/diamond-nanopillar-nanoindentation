import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from datetime import datetime

def parse_lammps_log(filepath):
    """
    Parse LAMMPS log file with all energy terms and geometric parameters
    """
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: File {filepath} not found")
        return None

    # Updated columns based on your thermo_style
    columns = [
        "step", "temp", "press", "pe", "ke", "etotal",
        "posz", "fz", "pillar_height", "pillar_volume",
        "stress_zz", "time_ps", "dz", "strain", "temp_avg"
    ]

    data = []
    collecting = False
    
    for line in lines:
        line = line.strip()
        if line.startswith("Step"):
            collecting = True
            continue
        if collecting:
            parts = line.split()
            # Check if we have the right number of columns
            if len(parts) == len(columns):
                try:
                    # Convert all parts to float to validate
                    float_parts = [float(p) for p in parts]
                    data.append(parts)
                except ValueError:
                    # Skip lines that can't be converted to numbers
                    continue
            else:
                # Stop collecting if format changes
                if len(parts) < len(columns) and len(data) > 0:
                    break

    if not data:
        print(f"Warning: No valid data found in {filepath}")
        return None

    df = pd.DataFrame(data, columns=columns)
    
    # Convert to numeric, handling any remaining non-numeric values
    for col in columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    
    # Drop rows with any NaN values
    initial_rows = len(df)
    df = df.dropna()
    final_rows = len(df)
    
    if final_rows < initial_rows:
        print(f"Warning: Dropped {initial_rows - final_rows} rows with invalid data from {filepath}")
    
    if len(df) == 0:
        print(f"Error: No valid data remaining after cleaning {filepath}")
        return None
        
    return df

def setup_plot_style():
    """Set consistent plot style"""
    plt.rcParams.update({
        'figure.figsize': (12, 7),
        'lines.linewidth': 2,
        'font.size': 12,
        'axes.labelsize': 14,
        'axes.titlesize': 16,
        'legend.fontsize': 12,
        'axes.grid': True,
        'grid.alpha': 0.5,
        'grid.linestyle': '--',
        'axes.axisbelow': True,
        'axes.labelpad': 10,
        'figure.dpi': 100
    })

def save_plot_dual_x(df, y_column, ylabel, title, filename, format_y=None):
    """
    Create plots with dual x-axes (timestep and time in ps)
    """
    if df is None or len(df) == 0:
        print(f"Warning: Cannot create plot {filename} - no data available")
        return
        
    fig, ax1 = plt.subplots()
    
    # Main plot
    ax1.plot(df["step"], df[y_column], color='#1f77b4', label=ylabel)
    ax1.set_xlabel("Timestep")
    ax1.set_ylabel(ylabel)
    ax1.set_title(title)
    ax1.grid(True, linestyle='--', alpha=0.7)
    
    if format_y:
        ax1.yaxis.set_major_formatter(plt.FuncFormatter(format_y))

    # Secondary x-axis for time - only if we have time data
    if "time_ps" in df.columns and not df["time_ps"].isna().all():
        ax2 = ax1.twiny()
        ax2.set_xlim(ax1.get_xlim())
        
        # Get current x-ticks from primary axis
        xticks = ax1.get_xticks()
        
        # Only proceed if we have valid data for interpolation
        if len(df) > 1 and not df["step"].isna().all() and not df["time_ps"].isna().all():
            try:
                # Filter out any NaN values for interpolation
                mask = ~(df["step"].isna() | df["time_ps"].isna())
                if mask.sum() > 1:  # Need at least 2 points for interpolation
                    step_vals = df.loc[mask, "step"].values
                    time_vals = df.loc[mask, "time_ps"].values
                    
                    # Only interpolate for xticks within our data range
                    valid_xticks = xticks[(xticks >= step_vals.min()) & (xticks <= step_vals.max())]
                    
                    if len(valid_xticks) > 0:
                        times_ps = np.interp(valid_xticks, step_vals, time_vals)
                        ax2.set_xticks(valid_xticks)
                        ax2.set_xticklabels([f"{t:.2f}" for t in times_ps])
                        ax2.set_xlabel("Time (ps)")
            except Exception as e:
                print(f"Warning: Could not create time axis for {filename}: {e}")

    plt.tight_layout()
    plt.savefig(str(filename), dpi=300, bbox_inches='tight')
    plt.close()

def cumulative_work(df):
    """Calculate cumulative work from force and displacement"""
    if df is None or len(df) < 2:
        return np.array([])
        
    depth = df["dz"].to_numpy()
    force = df["fz"].to_numpy()
    work = np.zeros_like(force)
    
    for i in range(1, len(force)):
        dx = depth[i] - depth[i-1]
        work[i] = work[i-1] + 0.5 * (force[i] + force[i-1]) * dx
    return work

def save_comparison_plot(df300, df1000, x_col, y_col, xlabel, ylabel, title, filename, format_y=None):
    """
    Create comparison plots between 300K and 1000K data
    """
    fig, ax = plt.subplots()
    
    # Plot both temperatures if data exists
    if df300 is not None and len(df300) > 0 and x_col in df300.columns and y_col in df300.columns:
        ax.plot(df300[x_col], df300[y_col], label="300K", color='#1f77b4')
    
    if df1000 is not None and len(df1000) > 0 and x_col in df1000.columns and y_col in df1000.columns:
        ax.plot(df1000[x_col], df1000[y_col], label="1000K", color='#d62728')
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.7)
    
    if format_y:
        ax.yaxis.set_major_formatter(plt.FuncFormatter(format_y))
    
    plt.tight_layout()
    plt.savefig(str(filename), dpi=300, bbox_inches='tight')
    plt.close()

def plot_energy_components(df, temp, filename):
    """Plot all energy components together"""
    if df is None or len(df) == 0:
        print(f"Warning: Cannot create energy plot for {temp} - no data available")
        return
        
    fig, ax = plt.subplots()
    
    # Check which energy columns exist
    if "pe" in df.columns:
        ax.plot(df["time_ps"], df["pe"], label="Potential", color='#1f77b4')
    if "ke" in df.columns:
        ax.plot(df["time_ps"], df["ke"], label="Kinetic", color='#d62728')
    if "etotal" in df.columns:
        ax.plot(df["time_ps"], df["etotal"], label="Total", color='#2ca02c')
    
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("Energy (eV)")
    ax.set_title(f"Energy Components @{temp}")
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(str(filename), dpi=300, bbox_inches='tight')
    plt.close()

def format_scientific(x, p):
    """Format numbers in scientific notation"""
    return f"{x:.2e}"

def print_stats(df, temp):
    """Print statistics for a dataset"""
    if df is None or len(df) == 0:
        print(f"\nNo data available for {temp}")
        return
        
    print(f"\nStatistics for {temp}:")
    print(f"Data points: {len(df)}")
    
    stats = {}
    if "fz" in df.columns:
        stats["Maximum Force"] = f"{df['fz'].max():.2f} eV/Å"
    if "stress_zz" in df.columns:
        stats["Maximum Stress"] = f"{df['stress_zz'].max():.2e} eV/Å³"
    if "dz" in df.columns and len(df) > 0:
        stats["Final Indentation Depth"] = f"{df['dz'].iloc[-1]:.2f} Å"
    if "temp" in df.columns:
        stats["Average Temperature"] = f"{df['temp'].mean():.2f} ± {df['temp'].std():.2f} K"
    if "etotal" in df.columns:
        stats["Energy Conservation"] = f"{df['etotal'].max() - df['etotal'].min():.2e} eV"
    if "pillar_height" in df.columns and len(df) > 0:
        stats["Final Pillar Height"] = f"{df['pillar_height'].iloc[-1]:.2f} Å"
    if "press" in df.columns:
        stats["Maximum Pressure"] = f"{df['press'].max():.2e} bars"
    
    for key, value in stats.items():
        print(f"{key:.<30} {value}")

def main():
    print(f"Nanoindentation Analysis Started: {datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')} UTC")
    print("User: akmal523\n")
    
    # Set plot style
    setup_plot_style()
    
    # Load logs
    print("Loading log files...")
    df300 = parse_lammps_log("log_300K.lammps")
    df1000 = parse_lammps_log("log_1000K.lammps")
    
    if df300 is None and df1000 is None:
        print("Error: No valid data files found. Please check your log files.")
        return

    # Calculate work if possible
    if df300 is not None and "dz" in df300.columns and "fz" in df300.columns:
        df300["work"] = cumulative_work(df300)
    if df1000 is not None and "dz" in df1000.columns and "fz" in df1000.columns:
        df1000["work"] = cumulative_work(df1000)

    # Generate individual temperature plots
    print("Generating individual plots...")
    for df, temp in [(df300, "300K"), (df1000, "1000K")]:
        if df is None:
            print(f"Skipping {temp} - no data available")
            continue
            
        # Force vs depth
        if "fz" in df.columns:
            save_plot_dual_x(
                df, "fz", "Force (eV/Å)", 
                f"Force vs Depth @{temp}", 
                f"force_depth_{temp}.png"
            )
        
        # Indenter position
        if "posz" in df.columns:
            save_plot_dual_x(
                df, "posz", "Indenter Position (Å)", 
                f"Indenter Position @{temp}", 
                f"indenter_position_{temp}.png"
            )
        
        # Stress evolution
        if "stress_zz" in df.columns:
            save_plot_dual_x(
                df, "stress_zz", "Stress (eV/Å³)", 
                f"Stress Evolution @{temp}", 
                f"stress_{temp}.png",
                format_scientific
            )
        
        # Temperature evolution
        if "temp" in df.columns:
            save_plot_dual_x(
                df, "temp", "Temperature (K)", 
                f"Temperature Evolution @{temp}", 
                f"temperature_{temp}.png"
            )
        
        # Energy components
        plot_energy_components(df, temp, f"energy_components_{temp}.png")

    # Generate comparison plots
    print("Generating comparison plots...")
    comparisons = [
        ("dz", "fz", "Depth (Å)", "Force (eV/Å)", "Force vs Depth Comparison", "force_depth_compare.png"),
        ("strain", "stress_zz", "Strain", "Stress (eV/Å³)", "Stress vs Strain Comparison", "stress_strain_compare.png", format_scientific),
        ("time_ps", "pe", "Time (ps)", "Potential Energy (eV)", "Potential Energy Evolution", "pe_time_compare.png"),
        ("time_ps", "ke", "Time (ps)", "Kinetic Energy (eV)", "Kinetic Energy Evolution", "ke_time_compare.png"),
        ("time_ps", "etotal", "Time (ps)", "Total Energy (eV)", "Total Energy Evolution", "etotal_time_compare.png"),
        ("time_ps", "press", "Time (ps)", "Pressure (bars)", "Pressure Evolution", "pressure_time_compare.png", format_scientific),
        ("time_ps", "pillar_height", "Time (ps)", "Pillar Height (Å)", "Pillar Height Evolution", "height_time_compare.png"),
        ("time_ps", "posz", "Time (ps)", "Indenter Position (Å)", "Indenter Position Comparison", "indenter_position_compare.png")
    ]

    # Add work comparison if work was calculated
    if (df300 is not None and "work" in df300.columns) or (df1000 is not None and "work" in df1000.columns):
        comparisons.append(("dz", "work", "Depth (Å)", "Cumulative Work (eV)", "Work vs Depth Comparison", "work_depth_compare.png"))

    for comparison in comparisons:
        try:
            save_comparison_plot(df300, df1000, *comparison)
        except Exception as e:
            print(f"Warning: Could not create comparison plot {comparison[4]}: {e}")

    # Print statistics
    print("\n=== Nanoindentation Analysis Results ===")
    print_stats(df300, "300K")
    print_stats(df1000, "1000K")
    
    print(f"\n✅ Analysis completed at {datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')} UTC")
    print("Generated plots saved to current directory.")

if __name__ == "__main__":
    main()