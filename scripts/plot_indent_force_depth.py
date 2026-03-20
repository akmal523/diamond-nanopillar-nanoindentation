import matplotlib.pyplot as plt

# Use Agg backend to avoid opening any window
plt.switch_backend('Agg')

# Example: load indentation data from a CSV/log file
# Replace with your actual data loading logic
# Here we simulate with sample data:
depth = [0, 0.5, 1.0, 1.5, 2.0]  # indentation depth (Angstrom)
force = [0, -1, -2.3, -3.8, -5.2]  # indentation force (eV/Angstrom)

plt.figure(figsize=(6,4))
plt.plot(depth, force, marker='o', linestyle='-', color='blue')
plt.xlabel('Indentation Depth (Å)')
plt.ylabel('Indentation Force (eV/Å)')
plt.title('Nanoindentation Force-Depth Curve')
plt.grid(True)

# Save figure to file (PDF, PNG, etc)
plt.tight_layout()
plt.savefig('indentation_force_depth.png', dpi=300)
print("Plot saved as indentation_force_depth.png")
