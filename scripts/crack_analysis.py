import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from datetime import datetime

class DumpReader:
    def __init__(self, filename):
        self.filename = filename
        self.timesteps = {}

    def read_dump(self):
        """Read LAMMPS dump file"""
        current_timestep = None
        atoms = []

        with open(self.filename, 'r') as f:
            lines = f.readlines()
            i = 0
            while i < len(lines):
                line = lines[i].strip()

                if line.startswith("ITEM: TIMESTEP"):
                    if current_timestep is not None and atoms:
                        self.timesteps[current_timestep] = np.array(atoms)
                        atoms = []
                    current_timestep = int(lines[i+1])
                    i += 2

                elif line.startswith("ITEM: NUMBER OF ATOMS"):
                    self.num_atoms = int(lines[i+1])
                    i += 2

                elif line.startswith("ITEM: BOX BOUNDS"):
                    self.box_bounds = [
                        [float(x) for x in lines[i+1].split()],
                        [float(x) for x in lines[i+2].split()],
                        [float(x) for x in lines[i+3].split()]
                    ]
                    i += 4

                elif line.startswith("ITEM: ATOMS"):
                    self.columns = line.split()[2:]
                    i += 1

                else:
                    try:
                        values = [float(x) for x in line.split()]
                        if len(values) == len(self.columns):
                            atoms.append(values)
                    except:
                        pass
                    i += 1

            if current_timestep is not None and atoms:
                self.timesteps[current_timestep] = np.array(atoms)

        return self.timesteps

def calculate_3d_coordination(positions, cutoff=2.0):
    """Calculate coordination numbers in 3D (simple distance-based cutoff)"""
    n_atoms = len(positions)
    coordination = np.zeros(n_atoms)

    for i in range(n_atoms):
        distances = np.sqrt(np.sum((positions - positions[i])**2, axis=1))
        coordination[i] = np.sum((distances > 0) & (distances < cutoff))

    return coordination

def create_visualization(dump_file, temperature):
    """Create 3D and 2D visualizations for the last timestep"""
    reader = DumpReader(dump_file)
    timesteps = reader.read_dump()

    # Get last timestep data
    last_timestep = max(timesteps.keys())
    positions = timesteps[last_timestep][:, 2:5]  # assuming x,y,z are columns 3,4,5
    coordination_3d = calculate_3d_coordination(positions)

    # Create figure with both 3D and 2D views
    fig = plt.figure(figsize=(20, 10))

    # 3D plot
    ax1 = fig.add_subplot(121, projection='3d')
    scatter1 = ax1.scatter(positions[:,0], positions[:,1], positions[:,2],
                          c=np.clip(coordination_3d, 0, 5), cmap='viridis',
                          s=20, alpha=0.6, vmin=0, vmax=5)
    ax1.set_title(f'{temperature}K - Final Structure (t={last_timestep}ps)')
    ax1.set_xlabel('X (Å)')
    ax1.set_ylabel('Y (Å)')
    ax1.set_zlabel('Z (Å)')

    # 2D projection
    ax2 = fig.add_subplot(122)
    scatter2 = ax2.scatter(positions[:,0], positions[:,1],
                          c=np.clip(coordination_3d, 0, 5), cmap='viridis',
                          s=20, alpha=0.6, vmin=0, vmax=5)
    ax2.set_title(f'{temperature}K - 2D Projection (t={last_timestep}ps)')
    ax2.set_xlabel('X (Å)')
    ax2.set_ylabel('Y (Å)')

    # Add colorbar manually on the right side
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
    cbar = fig.colorbar(scatter2, cax=cbar_ax)
    cbar.set_label('Coordination Number')

    # Adjust layout (avoid warning with 3D plots)
    fig.subplots_adjust(wspace=0.3, right=0.88)
    return fig, positions, coordination_3d, last_timestep

def analyze_crack_regions(positions, coordination):
    """Analyze regions with potential cracks (CN < 3)"""
    crack_mask = coordination < 3
    crack_positions = positions[crack_mask]
    crack_coordination = coordination[crack_mask]

    return {
        'total_atoms': len(positions),
        'crack_atoms': len(crack_positions),
        'crack_percentage': 100 * len(crack_positions) / len(positions),
        'mean_crack_cn': np.mean(crack_coordination) if len(crack_coordination) > 0 else 0,
        'crack_volume': np.prod(np.max(crack_positions, axis=0) - np.min(crack_positions, axis=0))
        if len(crack_positions) > 0 else 0
    }

def print_coordination_histogram(coordination, temp):
    """Print histogram of coordination numbers, clipped at 5"""
    max_cn = 5
    clipped = np.clip(coordination, 0, max_cn).astype(int)
    print(f"\nCoordination histogram for {temp}K (clipped to 0..{max_cn}):")
    for cn in range(1, max_cn+1):
        count = np.sum(clipped == cn)
        print(f"  CN={cn:2d}: {count:6d} atoms")

def main():
    print(f"\nCrack Analysis - Final Timestep")
    print(f"-------------------------------")
    print(f"Current Date and Time (UTC): {datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Current User's Login: akmal523\n")

    # Analyze both temperatures separately
    for temp, dump_file in [(300, "dump.indent_300K"), (1000, "dump.indent_1000K")]:
        fig, positions, coordination, last_timestep = create_visualization(dump_file, temp)

        # Save visualization
        fig.savefig(f'crack_analysis_{temp}K_t{last_timestep}.png', dpi=300, bbox_inches='tight')
        plt.close(fig)

        # Print histogram for this system
        print_coordination_histogram(coordination, temp)

        # Analyze crack regions
        crack_analysis = analyze_crack_regions(positions, coordination)

        print(f"\n{temp}K System Analysis (t = {last_timestep} ps):")
        print("-" * 50)
        print(f"Total atoms analyzed: {crack_analysis['total_atoms']}")
        print(f"Atoms in crack regions (CN < 3): {crack_analysis['crack_atoms']}")
        print(f"Percentage of structure damaged: {crack_analysis['crack_percentage']:.2f}%")
        print(f"Average coordination in crack regions: {crack_analysis['mean_crack_cn']:.2f}")
        print(f"Approximate crack volume: {crack_analysis['crack_volume']:.2f} Å³")

        # Analyze coordination by height
        z_positions = positions[:,2]
        z_min, z_max = np.min(z_positions), np.max(z_positions)
        n_layers = 5
        layer_thickness = (z_max - z_min) / n_layers

        print("\nCoordination analysis by height:")
        for i in range(n_layers):
            z_start = z_min + i * layer_thickness
            z_end = z_start + layer_thickness
            layer_mask = (z_positions >= z_start) & (z_positions < z_end)
            layer_cn = coordination[layer_mask]

            print(f"Layer {i+1} ({z_start:.1f}-{z_end:.1f} Å):")
            if len(layer_cn) > 0:
                print(f"  Mean CN: {np.mean(layer_cn):.2f} ± {np.std(layer_cn):.2f}")
                print(f"  Damaged atoms: {np.sum(layer_cn < 3)}/{len(layer_cn)}")
            else:
                print("  No atoms in this layer")

if __name__ == "__main__":
    main()
