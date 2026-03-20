# Potential files

This directory is a placeholder for the AIREBO potential file required by the simulations.

## Required file

`CH.airebo` — AIREBO potential for carbon/hydrogen interactions.

## How to obtain

The file is distributed with LAMMPS and is available at:

- Official LAMMPS potentials repository: https://github.com/lammps/lammps/blob/develop/potentials/CH.airebo
- Your local LAMMPS installation: `<lammps_root>/potentials/CH.airebo`

Copy or symlink it into this directory before running any simulation:

```bash
cp /path/to/lammps/potentials/CH.airebo potentials/CH.airebo
```

## Reference

Stuart, S. J., Tutein, A. B., & Harrison, J. A. (2000).  
A reactive potential for hydrocarbons with intermolecular interactions.  
*Journal of Chemical Physics*, 112(14), 6472–6486.  
https://doi.org/10.1063/1.481208
