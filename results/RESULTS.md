# Results Description

This document describes all output figures, data files, and their physical meaning.  
All results were produced from LAMMPS MD simulations of nanoindentation on a diamond nanopillar at 300 K and 1000 K using the AIREBO potential.

---

## Simulation overview

A spherical indenter (R = 5.0 Å, k = 2.0 eV/Å³) was driven into a carbon nanopillar at 0.1 Å/ps for 500 ps (500,000 steps × 1 fs timestep). Two independent runs were performed — one at 300 K and one at 1000 K — to compare temperature-dependent mechanical response.

---

## Output figures

### Force–depth curves

| File | Description |
|---|---|
| `force_depth_300K.png` | Force (eV/Å) vs indentation depth (Å) at 300 K. Shows initial elastic rise followed by abrupt force loss consistent with catastrophic brittle failure. |
| `force_depth_1000K.png` | Same at 1000 K. Smoother profile with higher peak force, indicating ductile plastic deformation. |
| `force_depth_compare.png` | Direct overlay of both temperatures. The 1000 K curve reaches ~113.8 eV/Å vs ~58.4 eV/Å at 300 K. |
| `indentation_force_depth.png` | Quick diagnostic plot (from `plot_indent_force_depth.py`). |

### Indenter position

| File | Description |
|---|---|
| `indenter_position_300K.png` | z-position of the indenter centre vs time at 300 K. |
| `indenter_position_1000K.png` | Same at 1000 K. |
| `indenter_position_compare.png` | Side-by-side comparison. Both start at z_top + 5 Å and descend at 0.1 Å/ps. |

### Force evolution over time

| File | Description |
|---|---|
| `indentation_force_300K.png` | Indentation load vs simulation step at 300 K (after jump-region removal). |
| `indentation_force_1000K.png` | Same at 1000 K. |

### Force–time–position correlation

| File | Description |
|---|---|
| `force_time_position_correlation_300K.png` | Three-panel plot: force magnitude (mean and max), indenter z-position, and system height vs time at 300 K. Vertical dashed lines mark critical events detected automatically. |
| `force_time_position_correlation_1000K.png` | Same at 1000 K. The smoother height evolution confirms stable plastic deformation. |

### Stress and strain

| File | Description |
|---|---|
| `stress_300K.png` | True stress (eV/Å³) vs simulation step at 300 K. Max stress: 3.65 × 10⁶ eV/Å³. |
| `stress_1000K.png` | Same at 1000 K. Max stress: 3.11 × 10⁶ eV/Å³. |
| `stress_strain_compare.png` | Stress vs strain overlay for both temperatures. 300 K shows steeper initial slope (higher elastic modulus: 1.16 × 10⁷ vs 5.19 × 10⁵ eV/Å³). |

### Energy components

| File | Description |
|---|---|
| `energy_components_300K.png` | Potential, kinetic, and total energy vs time at 300 K. |
| `energy_components_1000K.png` | Same at 1000 K. Average kinetic energy: 64.96 eV (300 K) vs 214.40 eV (1000 K), confirming effective temperature control. |
| `pe_time_compare.png` | Potential energy evolution comparison. |
| `ke_time_compare.png` | Kinetic energy evolution comparison. |
| `etotal_time_compare.png` | Total energy evolution. Total work: 875.33 eV (300 K), 845.57 eV (1000 K). |

### Pressure and thermodynamics

| File | Description |
|---|---|
| `pressure_time_compare.png` | System pressure (bars) vs time. 300 K average pressure: −69.4 bars (tensile); 1000 K: +136 bars (compressive). Different signs indicate fundamentally different internal stress states. |

### Pillar geometry

| File | Description |
|---|---|
| `height_time_compare.png` | Pillar height vs time for both temperatures. 300 K: 29.91 Å → 55.96 Å (+26 Å, structural failure artefact). 1000 K: 29.91 Å → 33.68 Å (+3.78 Å, controlled plastic deformation). |
| `structure_evolution.png` | Two-panel plot of height and volume evolution over time for both systems. |
| `work_depth_compare.png` | Cumulative work (eV) vs indentation depth for both temperatures. Area under the curve represents total energy deposited into the system. |

### Temperature profiles

| File | Description |
|---|---|
| `temperature_300K.png` | Temperature vs time at 300 K. Mean: 260.54 ± 29.44 K (11.3% fluctuation). |
| `temperature_1000K.png` | Temperature vs time at 1000 K. Mean: 859.87 ± 49.27 K (5.7% fluctuation — better stability at higher T). |

### Time evolution (coordination and damage)

| File | Description |
|---|---|
| `time_events_300K.png` | Four-panel: mean coordination number, damaged structure percentage, system height, and damaged atom count vs time at 300 K. Coordination cutoff: 2.0 Å. Atoms with CN < 3 classified as damaged. |
| `time_events_1000K.png` | Same at 1000 K. |

### Crack and structural analysis

| File | Description |
|---|---|
| `crack_analysis_300K_t500000.png` | 3D scatter and 2D projection of atom positions at the final timestep (t = 500 ps) for the 300 K system. Atoms coloured by coordination number (0–5). Atoms with CN < 3 indicate bond-breaking and crack-like regions. |
| `crack_analysis_1000K_t500000.png` | Same for 1000 K. Fewer low-coordination atoms, confirming better structural integrity. |

### Diagnostic plots

| File | Description |
|---|---|
| `diagnostic_300K_temp.png` | Temperature diagnostic at 300 K. |
| `diagnostic_300K_v_fz.png` | Raw force diagnostic at 300 K. |
| `diagnostic_300K_v_true_stress_zz.png` | True stress diagnostic at 300 K. |
| `diagnostic_1000K_temp.png` | Temperature diagnostic at 1000 K. |
| `diagnostic_1000K_v_fz.png` | Raw force diagnostic at 1000 K. |
| `diagnostic_1000K_v_true_stress_zz.png` | True stress diagnostic at 1000 K. |

---

## Output data files

| File | Generated by | Description |
|---|---|---|
| `force_evolution_300K.dat` | `in.indent.300K` (LAMMPS) | Per-timestep: time (ps), force fz (eV/Å), true stress (eV/Å³), contact area (Å²). Columns separated by spaces. |
| `force_evolution_1000K.dat` | `in.indent.1000K` (LAMMPS) | Same for 1000 K. |
| `temp_profile_300K.dat` | `in.indent.300K` (LAMMPS) | Temperature profile along z-axis averaged every 1000 steps. Chunk-based output. |
| `temp_profile_1000K.dat` | `in.indent.1000K` (LAMMPS) | Same for 1000 K. |
| `cleaned_300K.csv` | `no_jump_data.py` | Log thermo data with artefact jump region removed (steps 100,000–106,000). |
| `cleaned_1000K.csv` | `no_jump_data.py` | Same for 1000 K (steps 300,000–311,000 removed). |
| `stats_300K.csv` | `no_jump_data.py` | Descriptive statistics (mean, std, min, max) for all thermo columns at 300 K. |
| `stats_1000K.csv` | `no_jump_data.py` | Same for 1000 K. |
| `stats_comparison_300K.csv` | `no_jump_data.py` | Side-by-side comparison of raw vs cleaned statistics at 300 K. |
| `stats_comparison_1000K.csv` | `no_jump_data.py` | Same for 1000 K. |
| `nanoindentation_summary.csv` | `compr_analysis.py` | One row per temperature case. Columns: max force, hardness at Pmax, max depth, height change, total work, elastic modulus (linear fit), Hertz reduced modulus (E_r), Hertz R², max stiffness (dF/dz), max stress, avg temperature, strain rate, avg pressure. |
| `nanoindentation_analysis.txt` | `compr_analysis.py` | Human-readable text report with all metrics from `nanoindentation_summary.csv`, plus pop-in candidate indices. |
| `summary_results.csv` | `analyse_final.py` | Per-case: Pmax, depth at Pmax, event time, Hertz E_r, early stiffness, contact area, hardness, total work, height change. |
| `summary_results_corrected.csv` | `recompute_from_log_and_plot.py` | Recalculated metrics using depth = zstart − zpos, with Hertz fit, hardness, and detected pop-in times. More accurate than raw log values. |
| `corrected_auto_refined_results.csv` | Intermediate analysis | Refined merged result table. |

---

## Physical interpretation

### 300 K — Brittle failure

At 300 K the nanopillar has insufficient thermal energy for atomic rearrangement. The elastic modulus is ~22× higher than at 1000 K. Once the indenter exceeds a critical depth, the structure fails catastrophically: atoms are displaced violently upward (final height 55.96 Å vs initial 29.91 Å). This is not plastic deformation — it is structural collapse. The negative average pressure (−69.4 bars) indicates a tensile-dominated internal stress state that accelerates fracture.

### 1000 K — Ductile response

At 1000 K the enhanced atomic mobility allows the structure to absorb indenter work through controlled plastic deformation. The much lower elastic modulus (5.19 × 10⁵ vs 1.16 × 10⁷ eV/Å³) reflects a more compliant material. The pillar height increases only 3.78 Å. Higher peak force (113.81 vs 58.44 eV/Å) combined with larger fluctuations (39.79 vs 18.77 eV/Å) is consistent with dynamic atomic rearrangements rather than elastic loading followed by fracture.

### Energy efficiency

Energy conversion efficiency (ΔE / W_total): 61.9% at 300 K vs 72.9% at 1000 K. At 300 K, the remaining ~33% is dissipated in the fracture event. At 1000 K, more energy goes into permanent structural modification.
