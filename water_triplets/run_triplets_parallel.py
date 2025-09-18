#!/usr/bin/env python3
"""
Run many triplet.py jobs in parallel on the local machine.

Example
-------
python run_triplets_parallel.py protein.pdb traj.dcd -t 5 --nprocs 8
"""

import argparse, os, sys, subprocess, itertools, time
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import MDAnalysis as mda

SCRIPT_DIR   = Path(__file__).resolve().parent
TRIPLET_PATH = SCRIPT_DIR / "triplet.py"      # => /full/path/to/triplet.py

# throw error if triplet.py is not found in the same directory as this script
if not TRIPLET_PATH.exists():
    sys.exit(f"ERROR: cannot find triplet.py in {SCRIPT_DIR}. "
             "Please put triplet.py, run_triplets_parallel.py, and the waterlib executable in the same directory.")

#  Worker helper (must be top-level function; enables it to work on both Mac and Linux)
def run_cmd(cmd):
    """
    Run a subprocess without buffering its output in RAM.
    Writes combined stdout/stderr to logs/<short-id>.log
    Returns (ok, cmd, err_hint).
    """
    from hashlib import md5
    log_dir = Path("logs")
    log_dir.mkdir(exist_ok=True)

    short = md5(" ".join(cmd).encode()).hexdigest()[:8]
    log_path = log_dir / f"{short}.log"

    # Clamp threading to avoid oversubscription crashes / stalls
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("MKL_NUM_THREADS", "1")
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    env.setdefault("NUMEXPR_NUM_THREADS", "1")

    with open(log_path, "w") as lf:
        ret = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, env=env)

    if ret.returncode != 0:
        return False, cmd, f"exit={ret.returncode} (see {log_path})"
    return True, cmd, ""


def residue_key(res):
    """
    Deterministic sort key and token for residues with insertion codes.
    Returns (sort_tuple, res_token, segid_str)
    sort_tuple ensures plain number < 'A' < 'B' < ...
    """
    rid   = int(getattr(res, "resid", res.resid))
    icode = (getattr(res, "icode", "") or "").upper()
    # segid usually mirrors chainID for PDB; derive from chainIDs if empty
    segid = (getattr(res, "segid", "") or "").strip()
    if not segid:
        try:
            chain_ids = np.unique(res.atoms.chainIDs)
            if len(chain_ids) == 1:
                segid = str(chain_ids[0]).strip()
        except Exception:
            segid = ""
    sort_key = (rid, "" if icode == "" else icode)
    token    = f"{rid}{icode}"  # e.g. '112', '112A'
    return sort_key, token, segid

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
    parser.add_argument("--skip", type=int, default=1,
                        help="Frame stride (default 1)")
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
                 cmd=["python", str(TRIPLET_PATH), processed_pdb_path, args.trajectory,
                      "--groupsFile", args.groupsFile, "--groupNum", str(i+1),
                      "-t", str(args.time),
                      "--hydrationCutoff", str(args.hydrationCutoff),
                      "--skip", str(args.skip),
                      "--outdir", args.outdir])
            for i in range(len(groups))
        ]

    elif args.multiChain:
        # sort residues deterministically by (resid, icode)
        residues = sorted(u.residues, key=lambda r: residue_key(r)[0])
        work_items = []
        for res in residues:
            _, res_token, segid = residue_key(res)
            cmd = ["python", str(TRIPLET_PATH), processed_pdb_path, args.trajectory,
                   "-res", str(res_token),
                   "-ch", str(segid),
                   "-t", str(args.time),
                   "--hydrationCutoff", str(args.hydrationCutoff),
                   "--skip", str(args.skip),
                   "--outdir", args.outdir]
            work_items.append(dict(kind="residue", cmd=cmd))

    else:  # single-chain residue mode (still include insertion codes)
        residues = sorted(u.residues, key=lambda r: residue_key(r)[0])
        work_items = []
        for res in residues:
            _, res_token, segid = residue_key(res)
            cmd = ["python", str(TRIPLET_PATH), processed_pdb_path, args.trajectory,
                   "-res", str(res_token),
                   "-t", str(args.time),
                   "--hydrationCutoff", str(args.hydrationCutoff),
                   "--skip", str(args.skip),
                   "--outdir", args.outdir]
            work_items.append(dict(kind="residue", cmd=cmd))

    total_jobs = len(work_items)
    print(f"Preparing {total_jobs} triplet jobs "
          f"({args.nprocs} parallel processes)…\n")

    # parallel execution with progress
    start = time.time()
    fails = []

    with ThreadPoolExecutor(max_workers=args.nprocs) as pool:
        futures = {}
        for w in work_items:
            futures[pool.submit(run_cmd, w["cmd"])] = w
            time.sleep(0.1)   # small ramp to prevent a thundering herd

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
