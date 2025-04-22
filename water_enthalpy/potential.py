#!/usr/bin/env python3
"""
Compute non‑bonded potential energy between a selected atom *group*
(e.g. one residue) and the solvent.

Outputs a CSV with three columns:
    coulomb, lj, total  (kJ/mol)

Example
-------
python potential.py ../protein_processed.gro ../traj.dcd --top ../topol.top -res 10 -t 5 --skip 50
"""

import os, sys, argparse, time, warnings
import numpy as np
import pandas as pd
import MDAnalysis as mda
from tqdm import tqdm

from openmm import unit, Platform, Context, VerletIntegrator
from openmm.app import PDBFile, GromacsGroFile, GromacsTopFile
from openmm.unit import nanometer, kilojoule_per_mole, picosecond
from openmm.app.forcefield import CutoffPeriodic

# ----------------------------------------------------------------------
#  Helper : build OpenMM system & tag forces
# ----------------------------------------------------------------------
def prepare_system(topfile, pdb, group_idx, solvent_idx,
                   cutoff_nm=1.0):
    """
    Returns (system, residue_set, solvent_set)
    """

    system = topfile.createSystem(nonbondedMethod = CutoffPeriodic,
                                  nonbondedCutoff = cutoff_nm * nanometer,
                                  constraints      = None)

    residue = set(group_idx)
    solvent = set(solvent_idx)

    for force in system.getForces():
        # NonbondedForce: separate Coulomb term via parameter offsets
        if force.__class__.__name__ == "NonbondedForce":
            force.setForceGroup(0)
            force.addGlobalParameter("group_scale",   1.0)
            force.addGlobalParameter("solvent_scale", 1.0)

            for i in range(force.getNumParticles()):
                q, s, eps = force.getParticleParameters(i)
                label = "group_scale" if i in residue else "solvent_scale"
                # zero out parameters
                force.setParticleParameters(i, 0, 0, 0)
                force.addParticleParameterOffset(label, i, q, s, eps)

            for i in range(force.getNumExceptions()):
                p1, p2, qprod, s, eps = force.getExceptionParameters(i)
                force.setExceptionParameters(i, p1, p2, 0, 0, 0)

        # CustomNonbondedForce: Amber FF stores LJ here → restrict to cross term
        elif force.__class__.__name__ == "CustomNonbondedForce":
            force.setForceGroup(1)
            force.addInteractionGroup(residue, solvent)
        else:
            force.setForceGroup(2)   # everything else ignored

    return system

#  energy helper
def coulomb_energy(ctx, group_scale, solvent_scale):
    ctx.setParameter("group_scale",   group_scale)
    ctx.setParameter("solvent_scale", solvent_scale)
    state = ctx.getState(getEnergy=True, groups={0})
    return state.getPotentialEnergy()

def main():
    parser = argparse.ArgumentParser(
        description="Measure potential energy between a selected atom "
                    "group and solvent.")
    parser.add_argument("protein", help="processed GRO or PDB (matches topology)")
    parser.add_argument("trajectory", help="trajectory file (e.g. DCD)")
    parser.add_argument("--top", default="topol.top",
                        help="Gromacs topology file (default topol.top)")
    parser.add_argument("-res","--resid", type=int)
    parser.add_argument("-ch","--chain")
    parser.add_argument("--selection", help="Raw MDAnalysis selection string")
    parser.add_argument("--groupsFile")
    parser.add_argument("--groupNum", type=int)
    parser.add_argument("-t","--time", type=float, default=5.0,
                        help="Last X ns to analyse (default 5)")
    parser.add_argument("--skip", type=int, default=50,
                        help="Frame stride (default 50; reduce for better confidence bars)")
    parser.add_argument("--nogpu", action="store_true",
                        help="Force CPU platform")
    parser.add_argument("--cutoff", type=float, default=0.55,
                        help="Cutoff distance in nm (default 0.55)")
    parser.add_argument("--norm_per_water", action="store_true",
                        help="Normalize per water residues in cutoff")
    parser.add_argument("--outdir", type=str, default="energies",
                        help="Output directory (default: energies)")
    args = parser.parse_args()

    # checks
    prot_path = args.protein
    if not os.path.exists(prot_path):
        sys.exit(f"ERROR: cannot find {prot_path}")
    if not os.path.exists(args.trajectory):
        sys.exit(f"ERROR: cannot find {args.trajectory}")
    if not os.path.exists(args.top):
        sys.exit(f"ERROR: cannot find topology {args.top}")

    # load trajectory and select atoms from your group
    u = mda.Universe(prot_path, args.trajectory)

    if args.selection:
        sel_str = args.selection
    elif args.groupsFile and args.groupNum:
        with open(args.groupsFile) as fh:
            sel_str = fh.readlines()[args.groupNum-1].strip()
    elif args.resid is not None:
        sel_str = (f"resid {args.resid} and segid {args.chain}"
                   if args.chain else f"resid {args.resid}")
    else:
        sys.exit("ERROR: need --selection OR --groupsFile/--groupNum OR -res")

    group_atoms = u.select_atoms(sel_str)
    if len(group_atoms) == 0:
        sys.exit("ERROR: selection returned zero atoms")

    group_idx = group_atoms.indices.tolist()
    solvent_atoms = u.select_atoms("resname SOL")
    solvent_idx = solvent_atoms.indices.tolist()

    # setup openmm system
    gro = GromacsGroFile(prot_path) # for later pdb support: pdb = PDBFile(prot_path)
    top = GromacsTopFile(args.top,
                         periodicBoxVectors=gro.getPeriodicBoxVectors())
    system = prepare_system(top, gro, group_idx, solvent_idx, cutoff_nm=args.cutoff) # pdb instead of gro for pdb support laters

    # pick platform
    if args.nogpu:
        plat = Platform.getPlatformByName("CPU")
    else:
        try:
            plat = Platform.getPlatformByName("CUDA")
        except Exception:
            warnings.warn("Using CPU platform (CUDA not requested / not available)")
            plat = Platform.getPlatformByName("CPU")

    integrator = VerletIntegrator(0.001*picosecond)
    context = Context(system, integrator, plat)

    # trajectory details
    nframes = len(u.trajectory)
    frames_to_analyse = int(args.time*1000 / u.trajectory.dt)
    first = max(0, nframes - frames_to_analyse)

    coul_list, lj_list = [], []

    for ts in tqdm(u.trajectory[first::args.skip]):
        positions_nm = ts.positions / 10.0  # Å → nm
        context.setPositions(positions_nm * nanometer)

        if args.norm_per_water:
            water_oxygens = 'resname SOL and name O*'
            group = sel_str
            shell_waters = u.select_atoms(f'{water_oxygens} and around {args.cutoff * 10} ({group})')
            n_within_cutoff = len(shell_waters)
            if n_within_cutoff == 0:
                n_within_cutoff = 1 # avoid division by zero

        total    = coulomb_energy(context, 1, 1)
        group_ce = coulomb_energy(context, 1, 0)
        solv_ce  = coulomb_energy(context, 0, 1)
        coul     = (total - group_ce - solv_ce).value_in_unit(kilojoule_per_mole)
        lj = context.getState(getEnergy=True, groups={1}).getPotentialEnergy()
        
        if args.norm_per_water:
            lj = lj / n_within_cutoff
            coul = coul / n_within_cutoff

        coul_list.append(coul)
        lj_list.append(lj.value_in_unit(kilojoule_per_mole))

    total_arr = np.array(coul_list) + np.array(lj_list)

    # save energies to CSV
    out_dir = args.outdir
    os.makedirs(out_dir, exist_ok=True)

    if args.selection:
        tag = sel_str.replace(" ","_")
    elif args.groupsFile:
        tag = f"group{args.groupNum}"
    elif args.chain:
        tag = f"res{args.resid}_chain{args.chain}"
    else:
        tag = f"res{args.resid}"

    base = os.path.splitext(os.path.basename(prot_path))[0]
    # exclude _processed from filename
    if base.endswith("_processed"):
        base = base[:-10]
    out_csv = os.path.join(out_dir, f"{base}_{tag}_energies.csv")

    df = pd.DataFrame({
        "coulomb": coul_list,
        "lj":      lj_list,
        "total":   total_arr
    })
    df.to_csv(out_csv, index=False)
    print(f"\nSaved energies from {len(df)} frames of the MD simulation to {out_csv}")

# ----------------------------------------------------------------------
if __name__ == "__main__":
    main()

