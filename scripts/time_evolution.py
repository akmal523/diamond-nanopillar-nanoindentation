import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

class TimeEvolutionAnalyzer:
    def __init__(self, dump_file):
        self.dump_file = dump_file
        self.timesteps = {}
        self.data = {}
        
    def read_dump(self):
        """Read LAMMPS dump file timestep by timestep"""
        current_timestep = None
        atoms = []
        
        with open(self.dump_file, 'r') as f:
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

    def calculate_coordination(self, positions, cutoff=2.0):
        """Calculate coordination numbers"""
        n_atoms = len(positions)
        coordination = np.zeros(n_atoms)
        
        for i in range(n_atoms):
            distances = np.sqrt(np.sum((positions - positions[i])**2, axis=1))
            coordination[i] = np.sum((distances > 0) & (distances < cutoff))
            
        return coordination

    def analyze_timestep(self, positions, coordination):
        """Analyze structural properties for a timestep"""
        crack_mask = coordination < 3
        z_positions = positions[:,2]
        
        return {
            'mean_cn': np.mean(coordination),
            'std_cn': np.std(coordination),
            'crack_percentage': 100 * np.sum(crack_mask) / len(coordination),
            'height': np.max(z_positions) - np.min(z_positions),
            'damaged_atoms': np.sum(crack_mask)
        }

    def analyze_time_evolution(self):
        """Analyze evolution of structural properties over time"""
        time_data = {
            'timesteps': [],
            'mean_cn': [],
            'std_cn': [],
            'crack_percentage': [],
            'height': [],
            'damaged_atoms': []
        }
        
        for timestep, atoms in sorted(self.timesteps.items()):
            positions = atoms[:, 2:5]  # assuming x,y,z are columns 3,4,5
            coordination = self.calculate_coordination(positions)
            analysis = self.analyze_timestep(positions, coordination)
            
            time_data['timesteps'].append(timestep)
            time_data['mean_cn'].append(analysis['mean_cn'])
            time_data['std_cn'].append(analysis['std_cn'])
            time_data['crack_percentage'].append(analysis['crack_percentage'])
            time_data['height'].append(analysis['height'])
            time_data['damaged_atoms'].append(analysis['damaged_atoms'])
            
        return time_data

def create_time_evolution_plots(temp):
    """Create time evolution plots for a given temperature"""
    analyzer = TimeEvolutionAnalyzer(f"dump.indent_{temp}K")
    analyzer.read_dump()
    time_data = analyzer.analyze_time_evolution()
    
    # Create figure with multiple subplots
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle(f'Time Evolution Analysis - {temp}K System', fontsize=16)
    
    # Plot 1: Mean Coordination Number
    ax1 = axes[0,0]
    ax1.plot(time_data['timesteps'], time_data['mean_cn'], 'b-')
    ax1.fill_between(time_data['timesteps'],
                     np.array(time_data['mean_cn']) - np.array(time_data['std_cn']),
                     np.array(time_data['mean_cn']) + np.array(time_data['std_cn']),
                     alpha=0.2)
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('Mean Coordination Number')
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Crack Percentage
    ax2 = axes[0,1]
    ax2.plot(time_data['timesteps'], time_data['crack_percentage'], 'r-')
    ax2.set_xlabel('Time (ps)')
    ax2.set_ylabel('Damaged Structure (%)')
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: System Height
    ax3 = axes[1,0]
    ax3.plot(time_data['timesteps'], time_data['height'], 'g-')
    ax3.set_xlabel('Time (ps)')
    ax3.set_ylabel('System Height (Å)')
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Number of Damaged Atoms
    ax4 = axes[1,1]
    ax4.plot(time_data['timesteps'], time_data['damaged_atoms'], 'm-')
    ax4.set_xlabel('Time (ps)')
    ax4.set_ylabel('Number of Damaged Atoms')
    ax4.grid(True, alpha=0.3)
    
    # Add timestamp and user info
    plt.figtext(0.02, 0.02, 
                f"Generated: {datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')} UTC\n"
                f"User: akmal523",
                fontsize=8)
    
    plt.tight_layout()
    
    # Save the figure
    plt.savefig(f'time_events_{temp}K.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return time_data

def main():
    print(f"Current Date and Time (UTC - YYYY-MM-DD HH:MM:SS formatted): "
          f"{datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Current User's Login: akmal523")
    
    # Process both temperatures
    for temp in [300, 1000]:
        print(f"\nProcessing {temp}K system...")
        time_data = create_time_evolution_plots(temp)
        
        # Print summary statistics
        print(f"\nSummary for {temp}K:")
        print(f"Initial mean CN: {time_data['mean_cn'][0]:.2f}")
        print(f"Final mean CN: {time_data['mean_cn'][-1]:.2f}")
        print(f"Maximum damage percentage: {max(time_data['crack_percentage']):.2f}%")
        print(f"Final height: {time_data['height'][-1]:.2f} Å")
        print(f"Total damaged atoms at end: {time_data['damaged_atoms'][-1]}")

if __name__ == "__main__":
    main()
