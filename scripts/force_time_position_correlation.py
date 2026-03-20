import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

class ForcePositionAnalyzer:
    def __init__(self, dump_file):
        self.dump_file = dump_file
        self.timesteps = {}
        
    def read_dump(self):
        """Read LAMMPS dump file with force and position data"""
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

    def analyze_force_position(self):
        """Extract force and position data over time"""
        time_data = {
            'timesteps': [],
            'mean_force': [],
            'max_force': [],
            'min_z': [],
            'max_z': [],
            'indentor_position': []
        }
        
        for timestep, atoms in sorted(self.timesteps.items()):
            positions = atoms[:, 2:5]  # x,y,z positions
            forces = atoms[:, 5:8]  # fx,fy,fz forces
            
            time_data['timesteps'].append(timestep)
            time_data['mean_force'].append(np.mean(np.abs(forces[:,2])))  # mean |Fz|
            time_data['max_force'].append(np.max(np.abs(forces[:,2])))    # max |Fz|
            time_data['min_z'].append(np.min(positions[:,2]))
            time_data['max_z'].append(np.max(positions[:,2]))
            # Assuming the highest z-position corresponds to the indentor
            time_data['indentor_position'].append(np.max(positions[:,2]))
            
        return time_data

def create_correlation_plot(temp):
    """Create force-time-position correlation plot"""
    analyzer = ForcePositionAnalyzer(f"dump.indent_{temp}K")
    analyzer.read_dump()
    time_data = analyzer.analyze_force_position()
    
    # Create figure with three subplots sharing x-axis
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 15), sharex=True)
    fig.suptitle(f'Force-Time-Position Correlation - {temp}K System', fontsize=16)
    
    # Plot 1: Force Evolution
    ax1.plot(time_data['timesteps'], time_data['mean_force'], 'b-', label='Mean Force')
    ax1.plot(time_data['timesteps'], time_data['max_force'], 'r--', label='Max Force')
    ax1.set_ylabel('Force (eV/Å)')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Plot 2: Position Evolution
    ax2.plot(time_data['timesteps'], time_data['indentor_position'], 'g-', label='Indentor Position')
    ax2.set_ylabel('Z Position (Å)')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    # Plot 3: System Height
    height = np.array(time_data['max_z']) - np.array(time_data['min_z'])
    ax3.plot(time_data['timesteps'], height, 'm-', label='System Height')
    ax3.set_xlabel('Time (ps)')
    ax3.set_ylabel('Height (Å)')
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    
    # Add timestamp and user info
    plt.figtext(0.02, 0.02, 
                f"Generated: {datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')} UTC\n"
                f"User: akmal523",
                fontsize=8)
    
    # Identify and mark critical events
    critical_events = identify_critical_events(time_data)
    
    for event in critical_events:
        t = event['time']
        for ax in [ax1, ax2, ax3]:
            ax.axvline(x=t, color='k', linestyle=':', alpha=0.3)
            ax.annotate(event['description'], 
                       xy=(t, ax.get_ylim()[1]),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=8, rotation=90)
    
    plt.tight_layout()
    
    # Save the figure
    plt.savefig(f'force_time_position_correlation_{temp}K.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return time_data, critical_events

def identify_critical_events(time_data):
    """Identify critical events in the data"""
    events = []
    
    # Find maximum force point
    max_force_idx = np.argmax(time_data['max_force'])
    events.append({
        'time': time_data['timesteps'][max_force_idx],
        'description': f'Max Force: {time_data["max_force"][max_force_idx]:.2f} eV/Å'
    })
    
    # Find rapid height change points
    height = np.array(time_data['max_z']) - np.array(time_data['min_z'])
    height_gradient = np.gradient(height)
    rapid_changes = np.where(np.abs(height_gradient) > np.std(height_gradient) * 2)[0]
    
    if len(rapid_changes) > 0:
        for idx in rapid_changes[:2]:  # Only mark first two major changes
            events.append({
                'time': time_data['timesteps'][idx],
                'description': f'Rapid height change: {height_gradient[idx]:.2f} Å/ps'
            })
    
    return events

def main():
    print(f"Current Date and Time (UTC): {datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Current User's Login: akmal523")
    
    # Process both temperatures
    for temp in [300, 1000]:
        print(f"\nProcessing {temp}K system...")
        time_data, events = create_correlation_plot(temp)
        
        print(f"\nCritical Events for {temp}K:")
        for event in events:
            print(f"Time: {event['time']} ps - {event['description']}")
        
        print("\nSummary Statistics:")
        print(f"Maximum force: {max(time_data['max_force']):.2f} eV/Å")
        print(f"Final height: {time_data['max_z'][-1] - time_data['min_z'][-1]:.2f} Å")
        print(f"Total simulation time: {time_data['timesteps'][-1]} ps")

if __name__ == "__main__":
    main()
