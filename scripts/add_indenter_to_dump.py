# file: add_indenter_to_dump.py
import re
import numpy as np
import pandas as pd
from math import sqrt, pi
from datetime import datetime

DT_PS = 0.001                 # ps per step (from timestep 0.001 in metal)
INDENTER_RADIUS_A = 5.0
INDENTER_TYPE = 2
POINT_SPACING_A = 1.0         # ~ surface point spacing (larger -> fewer points)

THERMO_COLS = [
    "step","temp","press","pe","ke","etotal",
    "posz","fz","pillar_height","pillar_volume",
    "stress_zz","time_ps","dz","strain","temp_avg"
]

def parse_log_for_posz(log_path):
    # returns DataFrame with step, time_ps, posz
    df = []
    collecting=False
    with open(log_path,'r') as f:
        for raw in f:
            line = raw.strip()
            if line.startswith("Step"):
                collecting=True; continue
            if collecting:
                parts=line.split()
                if len(parts)==len(THERMO_COLS):
                    try:
                        vals = list(map(float, parts))
                        df.append(vals)
                    except ValueError:
                        pass
                elif df:
                    break
    if not df:
        raise RuntimeError(f"No thermo block parsed from {log_path}")
    df = pd.DataFrame(df, columns=THERMO_COLS)
    return df[["step","time_ps","posz"]].copy()

def sphere_surface_points(center, radius, spacing=1.0):
    """
    Even-ish point set via latitude-longitude grid sized by spacing.
    Returns (N,3) array.
    """
    # circumference ~ 2*pi*R -> number of longitudes ~ 2*pi*R/spacing
    n_theta = max(8, int(2*np.pi*radius/spacing))
    # latitudes from 0..pi
    n_phi = max(6, int(np.pi*radius/spacing))
    phi = np.linspace(0, np.pi, n_phi)
    theta = np.linspace(0, 2*np.pi, n_theta)
    PH, TH = np.meshgrid(phi, theta, indexing="ij")
    x = center[0] + radius*np.sin(PH)*np.cos(TH)
    y = center[1] + radius*np.sin(PH)*np.sin(TH)
    z = center[2] + radius*np.cos(PH)
    pts = np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1)
    # remove near-duplicates and ensure within box numerically later if needed
    return pts

def map_timestep_to_posz(df_log):
    # straightforward mapping by integer step
    d = dict(zip(df_log["step"].astype(int).tolist(), df_log["posz"].tolist()))
    return d

def process_dump(dump_in, log_path, dump_out):
    print(f"[add-indenter] start {datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')} UTC")
    df_log = parse_log_for_posz(log_path)
    step2z = map_timestep_to_posz(df_log)

    with open(dump_in,'r') as fin, open(dump_out,'w') as fout:
        current_step = None
        n_atoms_original = None
        atoms_written_this_frame = 0
        in_atoms = False
        extra_cols = 0
        header_atoms_line = ""

        # per-frame box bounds for center
        xlo=xhi=ylo=yhi=zlo=zhi=None

        for line in fin:
            if line.startswith("ITEM: TIMESTEP"):
                # flush previous frame state
                current_step = int(next(fin).strip())
                # write both lines out
                fout.write("ITEM: TIMESTEP\n")
                fout.write(f"{current_step}\n")
                # reset frame state
                in_atoms = False
                atoms_written_this_frame = 0
                n_atoms_original = None
                xlo=xhi=ylo=yhi=zlo=zhi=None
                continue

            if line.startswith("ITEM: NUMBER OF ATOMS"):
                n_atoms_original = int(next(fin).strip())
                # we don't yet know how many points; compute later using posz and box center
                # BUT we must immediately write NUMBER OF ATOMS. So compute now.
                # We can estimate points now: need center z => requires posz for this timestep.
                posz = step2z.get(current_step, None)
                if posz is None:
                    raise RuntimeError(f"No posz for timestep {current_step} from {log_path}")
                # temporary center; will get x,y after BOX BOUNDS, but count is independent of center
                pts = sphere_surface_points((0,0,posz), INDENTER_RADIUS_A, POINT_SPACING_A)
                total_atoms = n_atoms_original + len(pts)
                fout.write("ITEM: NUMBER OF ATOMS\n")
                fout.write(f"{total_atoms}\n")
                continue

            if line.startswith("ITEM: BOX BOUNDS"):
                # write header and three bounds lines, but keep values to compute center
                fout.write(line)
                # LAMMPS can have "pp f" flags; we just read numbers
                bounds=[]
                for _ in range(3):
                    bline = next(fin)
                    fout.write(bline)
                    vals = list(map(float, bline.strip().split()[:2]))
                    bounds.append(vals)
                (xlo,xhi),(ylo,yhi),(zlo,zhi) = bounds
                continue

            if line.startswith("ITEM: ATOMS"):
                header_atoms_line = line.rstrip("\n")
                fout.write(header_atoms_line + "\n")
                in_atoms = True
                atoms_written_this_frame = 0
                # figure column count from header
                # example header: ITEM: ATOMS id type x y z c_pe_atom c_disp_vec[1] ...
                header_cols = header_atoms_line.split()[2:]
                # ensure standard layout: id, type, x, y, z are first 5 variables for us
                # count total columns in lines by reading one atom line lookahead
                # BUT we cannot consume the next line prematurely; instead, we’ll infer from each atom line length.
                continue

            if in_atoms:
                if atoms_written_this_frame < n_atoms_original:
                    # copy original atom
                    fout.write(line)
                    atoms_written_this_frame += 1
                    if atoms_written_this_frame == n_atoms_original:
                        # append the indenter particles now
                        if any(v is None for v in [xlo,xhi,ylo,yhi,zlo,zhi]):
                            # in rare dump orders, BOX BOUNDS always precedes ATOMS; safety check
                            raise RuntimeError("Missing BOX BOUNDS before ATOMS block.")
                        xc = 0.5*(xlo+xhi); yc = 0.5*(ylo+yhi)
                        posz = step2z.get(current_step, None)
                        if posz is None:
                            raise RuntimeError(f"No posz for timestep {current_step}")
                        pts = sphere_surface_points((xc,yc,posz), INDENTER_RADIUS_A, POINT_SPACING_A)

                        # detect number of columns from the last real atom line we just wrote
                        cols = len(line.strip().split())
                        next_id = n_atoms_original + 1
                        zeros_tail = ["0.000000"] * (cols - 5) if cols > 5 else []

                        for i, (x,y,z) in enumerate(pts, start=0):
                            new_fields = [
                                str(next_id + i),                 # id
                                str(INDENTER_TYPE),               # type
                                f"{x:.6f}", f"{y:.6f}", f"{z:.6f}"  # x y z
                            ] + zeros_tail
                            fout.write(" ".join(new_fields) + "\n")
                    continue
                else:
                    # after we've written added particles, any stray lines (shouldn't happen) are passed
                    fout.write(line)
                    continue

            # default: copy through
            fout.write(line)

    print(f"✓ created {dump_out}")
    print(f"  - type {INDENTER_TYPE} sphere added (R={INDENTER_RADIUS_A} Å)")
    print(f"  - per-frame center from BOX BOUNDS; z from log posz")
    print(f"  - spacing ~{POINT_SPACING_A} Å → particles per sphere scales with R/spacing")

if __name__ == "__main__":
    for T in ["300K","1000K"]:
        dump_in  = f"dump.indent_{T}"
        log_in   = f"log_{T}.lammps"
        dump_out = f"dump.indent_{T}_with_indenter"
        process_dump(dump_in, log_in, dump_out)
