#!/usr/bin/env python3

import numpy as np
import MDAnalysis as mda
from tqdm import tqdm
import sys, os, argparse, shutil
from glob import glob
try:
    import waterlib as wl
except ImportError:
    sys.exit(
        "\nERROR: Could not import waterlib.\n"
        "Please compile it first (from the 'water_triplets' directory):\n\n"
        "    python setup.py build_ext --inplace\n\n"
        "Make sure waterlib.c exists and that you are running Python in the same directory.\n"
    )


parser = argparse.ArgumentParser(
    description=("Compile water triplet angles around the residue/"
                 "custom group in your protein.\n"
                 "Example: python triplet.py protein.gro traj.dcd -res 42"))
parser.add_argument('protein')
parser.add_argument('trajectory')
parser.add_argument('-res','--resid', type=int)
parser.add_argument('-ch','--chain')
parser.add_argument('-t','--time',   type=float, default=5.0,
                    help='Last x ns to analyse (default 5)')
parser.add_argument('--groupsFile')
parser.add_argument('--groupNum', type=int)
parser.add_argument('--selection')
parser.add_argument('--hydrationCutoff', type=float, default=4.25)
parser.add_argument('--hydrationLowCutoff', type=float, default=None)
args = parser.parse_args()

# some checks
protein_processed = args.protein if args.protein.endswith('.gro') else args.protein+'.gro'
if not os.path.exists(protein_processed):
    sys.exit(f"Error: cannot find {protein_processed}")
if not os.path.exists(args.trajectory):
    sys.exit(f"Error: cannot find {args.trajectory}")

protein_name = os.path.splitext(os.path.basename(protein_processed))[0].replace('_processed','')

if args.selection:
    group_string = args.selection.replace(" ","_")
    out_base = f'{protein_name}_{group_string}'
elif args.groupsFile and args.groupNum:
    out_base = f'{protein_name}_group{args.groupNum}'
elif args.chain:
    out_base = f'{protein_name}_res{args.resid}_chain{args.chain}'
else:
    out_base = f'{protein_name}_res{args.resid}'

if args.hydrationLowCutoff is not None:
    out_base += f'_lowC_{args.hydrationLowCutoff}_highC_{args.hydrationCutoff}'

angles_dir   = 'angles'
os.makedirs(angles_dir, exist_ok=True)
output_path  = os.path.join(angles_dir, out_base+'_angles.txt')

# reading trajectory
u = mda.Universe(protein_processed, args.trajectory) 
total_frames   = len(u.trajectory)
frame_dt_ps    = u.trajectory.dt
frames_to_load = int((args.time*1e3)/frame_dt_ps)
first_idx      = max(0, total_frames - frames_to_load)

# MDAnalysis atom selection
if args.selection:
    sel = args.selection
elif args.groupsFile and args.groupNum is not None:
    with open(args.groupsFile) as fh:
        sel = fh.readlines()[args.groupNum-1].strip()
elif args.resid is not None:
    sel = (f'resid {args.resid} and segid {args.chain} and not name H*'
           if args.chain else f'resid {args.resid} and not name H*')
else:
    sys.exit("Error: must specify -res/-ch, --groupsFile/--groupNum, or --selection")

# test selection
if len(u.select_atoms(sel)) == 0:
    sys.exit("Your atom selection picked zero atoms – please check.")

waters_sel = 'resname SOL and name O*'
BoxDims    = u.dimensions[:3]
lowCut     = 0.0
highCut    = 3.5

# checkpointing
done_lines = 0
if os.path.exists(output_path):
    with open(output_path) as fh:
        for done_lines,_ in enumerate(fh,1):
            pass
    print(f"Resuming – {done_lines} frames already present in {output_path}")

traj_slice = u.trajectory[first_idx+done_lines:]
with open(output_path, 'a') as out_f:
    for i, ts in tqdm(enumerate(traj_slice, start=done_lines),
                      total=frames_to_load-done_lines, unit='frame'):

        # hydration‑water selection
        if args.hydrationLowCutoff:
            shell = u.select_atoms(
                f'({waters_sel}) and around {args.hydrationCutoff} ({sel}) '
                f'and not around {args.hydrationLowCutoff} ({sel})')
        else:
            shell = u.select_atoms(
                f'({waters_sel}) and around {args.hydrationCutoff} ({sel})')

        # compute triplet angles (using C-implementation from waterlib)
        ang_np, _ = wl.triplet_angles(shell.positions,
                                      u.select_atoms(waters_sel).positions,
                                      BoxDims, lowCut, highCut)

        # write immediately
        out_f.write(" ".join(f"{a:.1f}" for a in ang_np) + "\n")

        # flush every 500 frames
        if i % 500 == 0:
            out_f.flush()

print(f"Done! Angles in {output_path}")
