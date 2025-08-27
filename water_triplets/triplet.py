#!/usr/bin/env python3

import numpy as np
import MDAnalysis as mda
from tqdm import tqdm
import sys, os, argparse, shutil
from glob import glob
import re
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module=r".*MDAnalysis\.topology\.PDBParser")
warnings.filterwarnings("ignore", category=UserWarning, module=r".*MDAnalysis\.topology\.guessers")
warnings.filterwarnings("ignore", category=DeprecationWarning, module=r".*MDAnalysis\.coordinates\.DCD")


try:
    import waterlib as wl
except ImportError:
    sys.exit(
        "\nERROR: Could not import waterlib.\n"
        "Please compile it first (from the 'water_triplets' directory):\n\n"
        "    python setup.py build_ext --inplace\n\n"
        "Make sure waterlib.c exists in the same directory as triplet.py.\n"
    )

# -------- helpers for insertion codes --------
def parse_res_token(tok: str):
    """
    Parse residue tokens like '112', '112A', '7c' -> (112, 'A'/'C' or '')
    """
    m = re.fullmatch(r'\s*(\d+)\s*([A-Za-z]?)\s*', str(tok))
    if not m:
        raise ValueError(f"Invalid residue token: {tok!r} (expected e.g. 112 or 112B)")
    resseq = int(m.group(1))
    icode = m.group(2).upper()
    return resseq, icode

def residue_selection(segid: str, res_token: str, heavy_only=True):
    """
    Build an MDAnalysis selection string for a residue, including insertion code if present.
    If segid is None/empty, omit it from the selection.
    """
    resseq, icode = parse_res_token(res_token)
    resid_clause = f"{resseq}{icode}"  # e.g., '112B' or '112'
    if segid:
        sel = f"segid {segid} and resid {resid_clause}"
    else:
        sel = f"resid {resid_clause}"
    if heavy_only:
        sel += " and not name H*"
    return sel, resseq, icode

# -------- CLI --------
parser = argparse.ArgumentParser(
    description=("Compile water triplet angles around the residue/"
                 "custom group in your protein.\n"
                 "Example: python triplet.py protein_processed.pdb traj.dcd -res 42"))
parser.add_argument('protein')
parser.add_argument('trajectory')
parser.add_argument('-res','--resid', type=str,  # allow tokens like 112B
                    help='Residue number (optionally with insertion code), e.g., 112 or 112B')
parser.add_argument('-ch','--chain')
parser.add_argument('-t','--time',   type=float, default=5.0,
                    help='Last x ns to analyse (default 5)')
parser.add_argument('--groupsFile')
parser.add_argument('--groupNum', type=int)
parser.add_argument('--selection')
parser.add_argument('--hydrationCutoff', type=float, default=4.25)
parser.add_argument('--hydrationLowCutoff', type=float, default=None)
parser.add_argument('--outdir', default='angles',
                    help='Output directory (default: angles)')
parser.add_argument('--skip', type=int, default=1,
                    help='Process every Nth frame (default 1 = every frame)')
args = parser.parse_args()

# some checks
protein_processed = args.protein if args.protein.endswith('.pdb') else args.protein+'.pdb'
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
elif args.resid is not None and args.chain:
    # include insertion code in filename if present (args.resid is a string token, e.g. '112B')
    out_base = f'{protein_name}_res{args.resid}_chain{args.chain}'
elif args.resid is not None:
    out_base = f'{protein_name}_res{args.resid}'
else:
    out_base = f'{protein_name}_res{args.resid}'

if args.hydrationLowCutoff is not None:
    out_base += f'_lowC_{args.hydrationLowCutoff}_highC_{args.hydrationCutoff}'

angles_dir   = args.outdir
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
    # Build selection supporting insertion codes; fallback to chainID if needed
    sel, resseq, icode = residue_selection(args.chain, args.resid, heavy_only=True)
    if len(u.select_atoms(sel)) == 0:
        # Optional fallback if topology exposes chainID but not segid
        resid_clause = f"{resseq}{icode}"
        if args.chain:
            sel_try = f'chainID {args.chain} and resid {resid_clause} and not name H*'
            if len(u.select_atoms(sel_try)) > 0:
                sel = sel_try
            else:
                sys.exit(f"Your atom selection picked zero atoms – tried '{sel}' and '{sel_try}'.")
        else:
            sys.exit(f"Your atom selection picked zero atoms – tried '{sel}'.")
else:
    sys.exit("Error: must specify -res/-ch, --groupsFile/--groupNum, or --selection")

# test selection
if len(u.select_atoms(sel)) == 0:
    sys.exit("Your atom selection picked zero atoms – please check.")

waters_sel = 'resname SOL and name O*'
BoxDims    = u.dimensions[:3]
lowCut     = 0.0
highCut    = 3.5

# how many frames back from the end we care about
frames_to_load = int((args.time*1e3)/u.trajectory.dt)
first_idx      = max(0, total_frames - frames_to_load)

# --- checkpointing on number of lines already in the output file ---
done_lines = 0
if os.path.exists(output_path):
    with open(output_path) as fh:
        for done_lines,_ in enumerate(fh,1):
            pass
    print(f"Resuming – {done_lines} frames already present in {output_path}")

# compute the list of frame‐indices we actually want to step through
start_frame   = first_idx + done_lines * args.skip
end_frame     = total_frames
stride        = args.skip
frame_indices = range(start_frame, end_frame, stride)

with open(output_path, 'a') as out_f:
    # total iterations = number of indices left
    total_iters = len(frame_indices)
    for write_idx, frame_idx in tqdm(
          enumerate(frame_indices, start=done_lines),
          total=total_iters, unit='frame'
    ):
        # explicitly jump to that frame
        u.trajectory[frame_idx]

        # hydration-water selection (unchanged)…  
        if args.hydrationLowCutoff:
            shell = u.select_atoms(
                f'({waters_sel}) and around {args.hydrationCutoff} ({sel}) '
                f'and not around {args.hydrationLowCutoff} ({sel})')
        else:
            shell = u.select_atoms(
                f'({waters_sel}) and around {args.hydrationCutoff} ({sel})')

        # compute triplet angles as before
        ang_np, _ = wl.triplet_angles(
            shell.positions,
            u.select_atoms(waters_sel).positions,
            BoxDims, lowCut, highCut
        )

        # write one line per processed frame
        out_f.write(" ".join(f"{a:.1f}" for a in ang_np) + "\n")

        # flush periodically
        if write_idx % 500 == 0:
            out_f.flush()

print(f"Done! Angles in {output_path}")
