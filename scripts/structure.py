import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

def parse_lammps_log(filepath):
    """Parse LAMMPS log file to get structural data"""
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: File {filepath} not found")
        return None

    # Updated columns based on your actual thermo_style
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
            if len(parts) == len(columns):
                try:
                    # Validate that all parts can be converted to float
                    [float(p) for p in parts]
                    data.append(parts)
                except ValueError:
                    continue
            else:
                # Stop if format changes
                if len(data) > 0:
                    break

    if not data:
        print(f"Warning: No valid data found in {filepath}")
        return None

    df = pd.DataFrame(data, columns=columns)
    
    # Convert to numeric
    for col in columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    
    # Remove rows with NaN values
    df = df.dropna()
    
    if len(df) == 0:
        print(f"Error: No valid data after cleaning {filepath}")
        return None
        
    return df

def create_structure_evolution_plot(log_300k, log_1000k):
    """Create comprehensive structure evolution plot"""
    print(f"Creating structure evolution plot at {datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')} UTC")
    print(f"User: akmal523")

    # Parse log files
    df300 = parse_lammps_log(log_300k)
    df1000 = parse_lammps_log(log_1000k)

    if df300 is None and df1000 is None:
        print("Error: No valid data available from either log file")
        return

    # Create figure with subplots
    n_plots = 2 if ('pillar_height' in (df300.columns if df300 is not None else []) and 
                    'pillar_volume' in (df300.columns if df300 is not None else [])) else 1
    
    fig, axes = plt.subplots(n_plots, 1, figsize=(12, 6*n_plots))
    if n_plots == 1:
        axes = [axes]
    
    fig.suptitle('Structural Evolution During Nanoindentation', fontsize=14)

    plot_idx = 0

    # Plot height evolution if data available
    if ((df300 is not None and 'pillar_height' in df300.columns) or 
        (df1000 is not None and 'pillar_height' in df1000.columns)):
        
        ax1 = axes[plot_idx]
        
        if df300 is not None and 'pillar_height' in df300.columns and len(df300) > 0:
            ax1.plot(df300['time_ps'], df300['pillar_height'], 
                     label='300K', color='blue', linewidth=2)
            
            # Add final height annotation
            final_height_300 = df300['pillar_height'].iloc[-1]
            final_time_300 = df300['time_ps'].iloc[-1]
            ax1.annotate(f'Final: {final_height_300:.2f} Å', 
                        xy=(final_time_300, final_height_300),
                        xytext=(10, 10), textcoords='offset points',
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="blue", alpha=0.3))
        
        if df1000 is not None and 'pillar_height' in df1000.columns and len(df1000) > 0:
            ax1.plot(df1000['time_ps'], df1000['pillar_height'], 
                     label='1000K', color='red', linewidth=2)
            
            # Add final height annotation
            final_height_1000 = df1000['pillar_height'].iloc[-1]
            final_time_1000 = df1000['time_ps'].iloc[-1]
            ax1.annotate(f'Final: {final_height_1000:.2f} Å',
                        xy=(final_time_1000, final_height_1000),
                        xytext=(10, -10), textcoords='offset points',
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="red", alpha=0.3))
        
        ax1.set_xlabel('Time (ps)')
        ax1.set_ylabel('Pillar Height (Å)')
        ax1.set_title('Height Evolution')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        plot_idx += 1

    # Plot volume evolution if data available
    if ((df300 is not None and 'pillar_volume' in df300.columns) or 
        (df1000 is not None and 'pillar_volume' in df1000.columns)) and plot_idx < len(axes):
        
        ax2 = axes[plot_idx]
        
        volume_change_300 = None
        volume_change_1000 = None
        
        if df300 is not None and 'pillar_volume' in df300.columns and len(df300) > 0:
            ax2.plot(df300['time_ps'], df300['pillar_volume'], 
                     label='300K', color='blue', linewidth=2)
            
            # Calculate volume change
            volume_change_300 = df300['pillar_volume'].iloc[-1] - df300['pillar_volume'].iloc[0]
            final_time_300 = df300['time_ps'].iloc[-1]
            final_volume_300 = df300['pillar_volume'].iloc[-1]
            
            ax2.annotate(f'ΔV: {volume_change_300:.2f} Å³', 
                        xy=(final_time_300, final_volume_300),
                        xytext=(10, 10), textcoords='offset points',
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="blue", alpha=0.3))
        
        if df1000 is not None and 'pillar_volume' in df1000.columns and len(df1000) > 0:
            ax2.plot(df1000['time_ps'], df1000['pillar_volume'], 
                     label='1000K', color='red', linewidth=2)
            
            # Calculate volume change
            volume_change_1000 = df1000['pillar_volume'].iloc[-1] - df1000['pillar_volume'].iloc[0]
            final_time_1000 = df1000['time_ps'].iloc[-1]
            final_volume_1000 = df1000['pillar_volume'].iloc[-1]
            
            ax2.annotate(f'ΔV: {volume_change_1000:.2f} Å³',
                        xy=(final_time_1000, final_volume_1000),
                        xytext=(10, -10), textcoords='offset points',
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="red", alpha=0.3))
        
        ax2.set_xlabel('Time (ps)')
        ax2.set_ylabel('Pillar Volume (Å³)')
        ax2.set_title('Volume Evolution')
        ax2.grid(True, alpha=0.3)
        ax2.legend()

    # Adjust layout and save
    plt.tight_layout()
    plt.savefig('structure_evolution.png', dpi=300, bbox_inches='tight')
    plt.close()

    print("✓ Created structure_evolution.png")
    
    # Print summary
    if df300 is not None and len(df300) > 0:
        print(f"300K Summary:")
        if 'pillar_height' in df300.columns:
            print(f"  - Initial height: {df300['pillar_height'].iloc[0]:.2f} Å")
            print(f"  - Final height: {df300['pillar_height'].iloc[-1]:.2f} Å")
            print(f"  - Height change: {df300['pillar_height'].iloc[0] - df300['pillar_height'].iloc[-1]:.2f} Å")
        if 'pillar_volume' in df300.columns:
            volume_change = df300['pillar_volume'].iloc[-1] - df300['pillar_volume'].iloc[0]
            print(f"  - Volume change: {volume_change:.2f} Å³")
    
    if df1000 is not None and len(df1000) > 0:
        print(f"1000K Summary:")
        if 'pillar_height' in df1000.columns:
            print(f"  - Initial height: {df1000['pillar_height'].iloc[0]:.2f} Å")
            print(f"  - Final height: {df1000['pillar_height'].iloc[-1]:.2f} Å")
            print(f"  - Height change: {df1000['pillar_height'].iloc[0] - df1000['pillar_height'].iloc[-1]:.2f} Å")
        if 'pillar_volume' in df1000.columns:
            volume_change = df1000['pillar_volume'].iloc[-1] - df1000['pillar_volume'].iloc[0]
            print(f"  - Volume change: {volume_change:.2f} Å³")

def main():
    try:
        create_structure_evolution_plot("log_300K.lammps", "log_1000K.lammps")
    except FileNotFoundError as e:
        print(f"Error: {e}")
        print("Please ensure log files exist in the current directory")
    except Exception as e:
        print(f"An error occurred: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()