#!/usr/bin/env python3
"""
Run many triplet.py jobs in parallel on the local machine.

Example
-------
python run_triplets_parallel.py protein.pdb traj.dcd -t 5 --nprocs 8
"""

import argparse, os, sys, subprocess, itertools, time
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import MDAnalysis as mda

#  Worker helper (must be top-level function; enables it to work on both Mac and Linux)
def run_cmd(cmd):
    """Run a shellless subprocess and return (returncode, cmd, stderr string)."""
    ret = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ret.returncode != 0:
        return False, cmd, ret.stderr.decode()
    return True, cmd, ""

if __name__ == "__main__":
    multiprocessing.set_start_method("spawn", force=True)

    parser = argparse.ArgumentParser(
        description="Launch one triplet.py job per residue (or custom group) "
                    "in parallel on local CPUs."
    )
    parser.add_argument("protein",  help="input PDB (unprocessed)")
    parser.add_argument("trajectory", help="trajectory file (e.g. traj.dcd)")
    parser.add_argument("--multiChain", action="store_true",
                        help="Protein has multiple chains")
    parser.add_argument("--groupsFile",
                        help="File with MDAnalysis selection strings, one per line")
    parser.add_argument("-t","--time", type=int, default=5,
                        help="Last X ns to analyse in each job (default 5)")
    parser.add_argument("--hydrationCutoff", type=float, default=4.25,
                        help="Hydration cutoff in Angstroms (default 4.25)")
    parser.add_argument("--outdir", type=str, default="angles",
                        help="Output directory for angles files (default 'angles')")
    parser.add_argument("--nprocs", type=int, default=os.cpu_count(),
                        help="Number of parallel processes (default = all CPUs)")
    
    args = parser.parse_args()
    if args.multiChain and args.groupsFile:
        sys.exit("ERROR: use EITHER --multiChain OR --groupsFile, not both")

    pdb_path  = args.protein if args.protein.endswith(".pdb") else args.protein + ".pdb"
    # check if "processed" is in pdb_path and throw a warning asking for the unprocessed pdb as input
    if "processed" in pdb_path:
        sys.exit("ERROR: I see 'processed' in the protein file name. Please provide the unprocessed PDB file as input (i.e. no waters/ions).")

    if not os.path.exists(pdb_path):
        sys.exit(f"ERROR: cannot find {pdb_path}")

    processed_pdb_path  = pdb_path[:-4] + "_processed.pdb"
    if not os.path.exists(processed_pdb_path):
        sys.exit(f"ERROR: cannot find processed PDB {processed_pdb_path}")

    if not os.path.exists(args.trajectory):
        sys.exit(f"ERROR: cannot find trajectory {args.trajectory}")

    protein_name = os.path.splitext(os.path.basename(pdb_path))[0]

    #  enumerate residues or custom groups
    u = mda.Universe(pdb_path)

    if args.groupsFile:
        with open(args.groupsFile) as fh:
            groups = [line.strip() for line in fh if line.strip() and not line.startswith("#")]
        work_items = [
            dict(kind="group",
                 cmd=["python", "triplet.py", processed_pdb_path, args.trajectory,
                      "--groupsFile", args.groupsFile, "--groupNum", str(i+1),
                      "-t", str(args.time),
                      "--hydrationCutoff", str(args.hydrationCutoff),
                      "--outdir", args.outdir])
            for i in range(len(groups))
        ]

    elif args.multiChain:
        resids = u.residues.resids
        segids = u.residues.segids
        work_items = [
            dict(kind="residue",
                 cmd=["python", "triplet.py", processed_pdb_path, args.trajectory,
                      "-res", str(rid), "-ch", str(sid),
                      "-t", str(args.time),
                      "--hydrationCutoff", str(args.hydrationCutoff),
                      "--outdir", args.outdir])
            for rid, sid in zip(resids, segids)
        ]

    else:  # single‑chain residue mode
        resids = u.residues.resids
        work_items = [
            dict(kind="residue",
                 cmd=["python", "triplet.py", processed_pdb_path, args.trajectory,
                      "-res", str(rid),
                      "-t", str(args.time),
                      "--hydrationCutoff", str(args.hydrationCutoff),
                      "--outdir", args.outdir])
            for rid in resids
        ]

    total_jobs = len(work_items)
    print(f"Preparing {total_jobs} triplet jobs "
          f"({args.nprocs} parallel processes)…\n")

    # parallel execution with progress
    start = time.time()
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

    elapsed = time.time() - start
    print("\n")

    # summarize
    if fails:
        print(f"Finished with {len(fails)} failures out of {total_jobs} jobs.\n")
        for cmd, err in fails:
            print("Command:", " ".join(cmd))
            print("Error  :", err or "(no stderr)")
            print("-"*60)
    else:
        print(f"All {total_jobs} jobs finished successfully in "
              f"{elapsed/60:.1f} min.")
