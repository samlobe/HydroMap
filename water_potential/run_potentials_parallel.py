#!/usr/bin/env python3
"""
Run many potential.py jobs in parallel on the local machine.
Ideal for CPU parallelization if CPUs are cheap and available.

Example
-------
python run_potentials_parallel.py ../protein.pdb ../traj.dcd --top ../topol.top -t 5 --skip 50 --nprocs 8
"""

import argparse, os, sys, subprocess, itertools
from time import time
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import MDAnalysis as mda

SCRIPT_DIR     = Path(__file__).resolve().parent
POTENTIAL_PATH = SCRIPT_DIR / "potential.py"      # => /full/path/to/potential.py

# throw error if potential.py is not found in the same directory as this script
if not POTENTIAL_PATH.exists():
    sys.exit(f"ERROR: cannot find potential.py in {SCRIPT_DIR}. "
             "Please put potential.py and run_potentials_parallel.py are in the same directory.")

def run_cmd(cmd):
    """Run a shellless subprocess and return (returncode, cmd, stderr string)."""
    ret = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ret.returncode != 0:
        return False, cmd, ret.stderr.decode()
    return True, cmd, ""

def residue_key(res):
    """
    Deterministic key/token for residues with insertion codes.
    Returns (sort_tuple, res_token, segid_str)
    """
    rid   = int(getattr(res, "resid", res.resid))
    icode = (getattr(res, "icode", "") or "").upper()
    segid = (getattr(res, "segid", "") or "").strip()
    if not segid:
        try:
            chain_ids = np.unique(res.atoms.chainIDs)
            if len(chain_ids) == 1:
                segid = str(chain_ids[0]).strip()
        except Exception:
            segid = ""
    sort_key = (rid, "" if icode == "" else icode)   # plain number first, then A<B<C...
    token    = f"{rid}{icode}"                        # e.g., '112', '112A'
    return sort_key, token, segid

if __name__ == "__main__":
    multiprocessing.set_start_method("spawn", force=True)

    parser = argparse.ArgumentParser(
        description="Launch one potential.py job per residue (or custom group) "
                    "in parallel on local CPUs."
    )
    parser.add_argument("protein",  help="input PDB (unprocessed)")
    parser.add_argument("trajectory", help="trajectory file (e.g. traj.dcd)")
    parser.add_argument("--top", default="topol.top",
                        help="Gromacs topology file (default topol.top)")
    parser.add_argument("--multiChain", action="store_true",
                        help="Protein has multiple chains")
    parser.add_argument("--groupsFile",
                        help="File with MDAnalysis selection strings, one per line")
    parser.add_argument("-t","--time", type=int, default=5,
                        help="Last X ns to analyse in each job (default 5)")
    parser.add_argument("--nprocs", type=int, default=os.cpu_count(),
                        help="Number of parallel processes (default = all CPUs)")
    parser.add_argument("--nogpu", action="store_true", default=True,  # default True: this runner prefers CPU
                        help="Force CPU platform (pass --nogpu to potential.py)")
    parser.add_argument("--skip", type=int, default=100,
                        help="Frame stride (default 100; reduce for better estimate of mean)")
    parser.add_argument("--cutoff", type=float, default=4.25,
                        help="Cutoff distance in Angstroms (default 4.25)")
    parser.add_argument("--outdir", type=str, default="potentials",
                        help="Output directory (default: potentials)")
    
    args = parser.parse_args()
    if args.multiChain and args.groupsFile:
        sys.exit("ERROR: use EITHER --multiChain OR --groupsFile, not both")

    pdb_path  = args.protein if args.protein.endswith(".pdb") else args.protein + ".pdb"
    if not os.path.exists(pdb_path):
        sys.exit(f"ERROR: cannot find {pdb_path}")

    processed_pdb_path  = pdb_path[:-4] + "_processed.pdb"
    if not os.path.exists(processed_pdb_path):
        sys.exit(f"ERROR: cannot find processed PDB {processed_pdb_path}")

    if not os.path.exists(args.trajectory):
        sys.exit(f"ERROR: cannot find trajectory {args.trajectory}")

    protein_name = os.path.splitext(os.path.basename(pdb_path))[0]

    #  parse the groups / residues
    u = mda.Universe(pdb_path)

    if args.groupsFile:
        with open(args.groupsFile) as fh:
            groups = [line.strip() for line in fh if line.strip() and not line.startswith("#")]
        work_items = [
            dict(kind="group",
                 cmd=["python", str(POTENTIAL_PATH), processed_pdb_path, args.trajectory,
                      "--top", args.top,
                      "--groupsFile", args.groupsFile, "--groupNum", str(i+1),
                      "-t", str(args.time),
                      "--skip", str(args.skip),
                      "--cutoff", str(args.cutoff),
                      "--outdir", str(args.outdir)]
                 + (["--nogpu"] if args.nogpu else []))
            for i in range(len(groups))
        ]

    elif args.multiChain:
        residues = sorted(u.residues, key=lambda r: residue_key(r)[0])
        work_items = []
        for res in residues:
            _, res_token, segid = residue_key(res)
            cmd = ["python", str(POTENTIAL_PATH), processed_pdb_path, args.trajectory,
                   "--top", args.top,
                   "-res", str(res_token), "-ch", str(segid),
                   "-t", str(args.time),
                   "--skip", str(args.skip),
                   "--cutoff", str(args.cutoff),
                   "--outdir", str(args.outdir)]
            if args.nogpu:
                cmd.append("--nogpu")
            work_items.append(dict(kind="residue", cmd=cmd))

    else:  # single-chain residue mode (still include insertion codes)
        residues = sorted(u.residues, key=lambda r: residue_key(r)[0])
        work_items = []
        for res in residues:
            _, res_token, _ = residue_key(res)
            cmd = ["python", str(POTENTIAL_PATH), processed_pdb_path, args.trajectory,
                   "--top", args.top,
                   "-res", str(res_token),
                   "-t", str(args.time),
                   "--skip", str(args.skip),
                   "--cutoff", str(args.cutoff),
                   "--outdir", str(args.outdir)]
            if args.nogpu:
                cmd.append("--nogpu")
            work_items.append(dict(kind="residue", cmd=cmd))

    total_jobs = len(work_items)
    print(f"Preparing {total_jobs} jobs to measure potentials "
          f"({args.nprocs} parallel processes)…\n")

    #  parallel execution with progress
    start = time()
    fails = []

    with ProcessPoolExecutor(max_workers=args.nprocs) as pool:
        futures = {pool.submit(run_cmd, w["cmd"]): w for w in work_items}
        done_so_far = 0
        for fut in as_completed(futures):
            ok, cmd, err = fut.result()
            done_so_far += 1
            if not ok:
                fails.append((cmd, err.strip()))
            print(f"\rCompleted {done_so_far}/{total_jobs}", end="", flush=True)

    elapsed = time() - start
    print("\n")

    #  summarize
    if fails:
        print(f"Finished with {len(fails)} failures out of {total_jobs} jobs.\n")
        for cmd, err in fails:
            print("Command:", " ".join(cmd))
            print("Error  :", err or "(no stderr)")
            print("-"*60)
    else:
        print(f"All {total_jobs} jobs finished successfully in "
              f"{elapsed/60:.1f} min.")
