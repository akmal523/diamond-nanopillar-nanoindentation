# Nanoindentation of a Diamond Nanopillar — MD Simulation

Molecular dynamics study of nanoindentation on a carbon nanopillar at **300 K** and **1000 K** using LAMMPS.  
Potential: AIREBO. Indenter: spherical, R = 5.0 Å, k = 2.0 eV/Å³, speed = 0.1 Å/ps.

---

## Repository structure

```
diamond-nanopillar-nanoindentation/
├── README.md
├── environment.yml              # Conda environment (Python + LAMMPS)
├── lammps/                      # LAMMPS input scripts (run order below)
│   ├── in.relax                 # Step 1 — geometry minimisation + short NVT
│   ├── in.nvt                   # Step 2 — NVT equilibration (alternate protocol)
│   ├── in.indent.300K           # Step 3a — nanoindentation at 300 K
│   └── in.indent.1000K          # Step 3b — nanoindentation at 1000 K
├── potentials/
│   └── CH.airebo                # AIREBO potential file (not tracked by git — see below)
├── scripts/                     # Python post-processing scripts
│   ├── correct.py               # Log/data correction utilities
│   ├── no_jump_data.py          # Remove artefact jump regions from log data
│   ├── compr_analysis.py        # Comprehensive analysis: Hertz fit, hardness, pop-in, stiffness
│   ├── analyse_final.py         # Final analysis pipeline from force_evolution_*.dat files
│   ├── recompute_from_log_and_plot.py  # Recompute depth & metrics directly from LAMMPS logs
│   ├── plot_main.py             # Main comparison plots (300K vs 1000K)
│   ├── plot_indent_force_depth.py  # Quick force-depth curve plot
│   ├── structure.py             # Pillar height and volume evolution
│   ├── time_evolution.py        # Coordination number and damage evolution over time
│   ├── force_time_position_correlation.py  # Force/position/height correlation plots
│   ├── crack_analysis.py        # 3D/2D crack and coordination analysis from dump files
│   └── add_indenter_to_dump.py  # Post-process dump: add spherical indenter ghost atoms
└── results/
    └── RESULTS.md               # Full description of all output figures and data files
```

---

## Simulation workflow

### 1. Environment setup

```bash
conda env create -f environment.yml
conda activate diamond-nanopillar-md
```

### 2. Prepare initial structure

Use `atomsk` to generate the diamond nanopillar and save as `nanopillar_bottom.lmp` and `bulk_diamond.lmp`.  
These `.lmp` data files are not tracked in git (they can be regenerated with atomsk).

### 3. Run LAMMPS simulations

```bash
# Step 1 — relax structure
lammps -in lammps/in.relax

# Step 2 — NVT equilibration (if using alternate protocol)
lammps -in lammps/in.nvt

# Step 3 — nanoindentation (run both in parallel or sequentially)
lammps -in lammps/in.indent.300K   > nohup_300K.out &
lammps -in lammps/in.indent.1000K  > nohup_1000K.out &
```

Each indentation run produces:
- `log_300K.lammps` / `log_1000K.lammps` — thermo output
- `dump.indent_300K` / `dump.indent_1000K` — atom trajectories
- `force_evolution_300K.dat` / `force_evolution_1000K.dat` — force data
- `temp_profile_300K.dat` / `temp_profile_1000K.dat` — temperature profiles

### 4. Post-processing

Run scripts in this order:

```bash
cd <working_dir_with_logs>

# 1. Clean jump artefacts from log data
python scripts/no_jump_data.py

# 2. Full comprehensive analysis (Hertz, hardness, pop-in)
python scripts/compr_analysis.py

# 3. Recompute metrics directly from logs
python scripts/recompute_from_log_and_plot.py

# 4. Main comparison plots
python scripts/plot_main.py

# 5. Structural evolution
python scripts/structure.py

# 6. Time evolution of coordination and damage
python scripts/time_evolution.py

# 7. Force-time-position correlation
python scripts/force_time_position_correlation.py

# 8. Crack analysis (requires dump files)
python scripts/crack_analysis.py

# 9. Add indenter to dump for OVITO visualisation
python scripts/add_indenter_to_dump.py
```

All figures are saved to `./figures/` and `./analysis_plots/`.

---

## Key simulation parameters

| Parameter | Value |
|---|---|
| Potential | AIREBO (Stuart et al., 2000) |
| Indenter geometry | Spherical |
| Indenter radius | 5.0 Å |
| Indenter stiffness | 2.0 eV/Å³ |
| Indentation speed | 0.1 Å/ps |
| Timestep | 1 fs (metal units) |
| Total run | 500,000 steps (500 ps) |
| Thermostat | NVT (Nosé-Hoover, τ = 0.1 ps) |
| Temperatures | 300 K and 1000 K |
| Boundary conditions | Periodic x, y; Free surface z |

---

## Key results summary

| Quantity | 300 K | 1000 K |
|---|---|---|
| Max force | 58.44 eV/Å | 113.81 eV/Å |
| Max force fluctuation | 18.77 eV/Å | 39.79 eV/Å |
| Elastic modulus (slope) | 1.16 × 10⁷ eV/Å³ | 5.19 × 10⁵ eV/Å³ |
| Initial pillar height | 29.91 Å | 29.91 Å |
| Final pillar height | 55.96 Å | 33.68 Å |
| Total work | 875.33 eV | 845.57 eV |
| Avg temperature | 260.54 ± 29.44 K | 859.87 ± 49.27 K |
| Deformation mode | Brittle (structural failure) | Ductile (plastic deformation) |

Full description of all output figures and data files: see [`results/RESULTS.md`](results/RESULTS.md).

---

## Files not tracked by git

The following files are excluded via `.gitignore` because they are either too large, binary, or regenerable:

- `*.lmp` structure files (regenerate with atomsk)
- `dump.*` trajectory files (can be hundreds of MB)
- `potentials/CH.airebo` (download from LAMMPS distribution)
- `nohup*.out`, `*.log` runtime output
- `__pycache__/`, `.ipynb_checkpoints/`

---

## Reference

Razikulov, A. (2025). *Nanoindentation Analysis of a Carbon Nanopillar at Different Temperatures*. Molecular Dynamics Simulation Report.

Stuart, S. J., Tutein, A. B., & Harrison, J. A. (2000). A reactive potential for hydrocarbons with intermolecular interactions. *Journal of Chemical Physics*, 112(14), 6472–6486.
