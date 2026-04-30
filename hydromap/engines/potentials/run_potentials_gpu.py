#!/usr/bin/env python3
"""
Run many potential jobs on one process.

Execution engine:
- inproc: in-process OpenMM runner with optional GPU shell counting for n_waters.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import MDAnalysis as mda
import numpy as np
import pandas as pd
from tqdm import tqdm

from openmm import Platform, Context, VerletIntegrator, Vec3
from openmm.unit import nanometer, kilojoule_per_mole, picosecond
from hydromap.analysis_inputs import NORMALIZED_WATER_OXYGEN_SELECTION, NORMALIZED_WATER_SELECTION
from hydromap.forcefield import load_openmm_system
from hydromap.selection import residue_sort_token_chain

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
POTENTIAL_PATH = REPO_ROOT / "hydromap" / "potentials" / "potential.py"
DEFAULT_GPU_GROUP_BATCH_SIZE = 10
MAX_GPU_GROUP_BATCH_SIZE = 64

if not POTENTIAL_PATH.exists():
    sys.exit(
        f"ERROR: cannot find potential.py in {SCRIPT_DIR}. "
        "Please ensure hydromap/potentials/potential.py is available."
    )


_GPU_SHELL_KERNEL = None
_CUPY = None
_GPU_SHELL_CUDA_SRC = r"""
extern "C" __global__
void shell_mask_kernel(
    const double* group_atoms,
    const int* group_offsets,
    int n_groups,
    const double* waters,
    int n_waters,
    const double* box,
    double high2,
    unsigned char* out_mask
) {
    long long tid = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    long long total = (long long)n_groups * (long long)n_waters;
    if (tid >= total) return;

    int g = (int)(tid / n_waters);
    int w = (int)(tid - (long long)g * n_waters);
    int a0 = group_offsets[g];
    int a1 = group_offsets[g + 1];

    if (a0 >= a1) {
        out_mask[tid] = 0;
        return;
    }

    double wx = waters[3 * w + 0];
    double wy = waters[3 * w + 1];
    double wz = waters[3 * w + 2];
    double min_r2 = 1.0e300;

    for (int a = a0; a < a1; ++a) {
        double dx = wx - group_atoms[3 * a + 0];
        double dy = wy - group_atoms[3 * a + 1];
        double dz = wz - group_atoms[3 * a + 2];
        dx -= box[0] * nearbyint(dx / box[0]);
        dy -= box[1] * nearbyint(dy / box[1]);
        dz -= box[2] * nearbyint(dz / box[2]);
        double r2 = dx * dx + dy * dy + dz * dz;
        if (r2 < min_r2) min_r2 = r2;
    }
    out_mask[tid] = (min_r2 <= high2) ? 1 : 0;
}
"""


def _configure_cupy_runtime_env() -> None:
    # Default to non-CUB reductions unless user explicitly overrides.
    os.environ.setdefault("CUPY_ACCELERATORS", "")

    # Opportunistically prepend common CUDA library locations if present.
    cuda_lib_candidates = (
        "/usr/local/cuda/lib64",
        "/usr/local/cuda/targets/x86_64-linux/lib",
        "/usr/local/cuda-12.4-temp/targets/x86_64-linux/lib",
    )
    ld_parts = [p for p in os.environ.get("LD_LIBRARY_PATH", "").split(":") if p]
    for libdir in cuda_lib_candidates:
        if os.path.isdir(libdir) and libdir not in ld_parts:
            ld_parts.insert(0, libdir)
    if ld_parts:
        os.environ["LD_LIBRARY_PATH"] = ":".join(ld_parts)


def _get_cupy():
    global _CUPY
    if _CUPY is None:
        _configure_cupy_runtime_env()
        try:
            import cupy as cp
        except ImportError as exc:
            raise RuntimeError("CuPy is required for GPU shell counting.") from exc
        _CUPY = cp
    return _CUPY


def _get_gpu_shell_kernel(cp):
    global _GPU_SHELL_KERNEL
    if _GPU_SHELL_KERNEL is None:
        mod = cp.RawModule(
            code=_GPU_SHELL_CUDA_SRC,
            options=("-std=c++11",),
            name_expressions=["shell_mask_kernel"],
        )
        _GPU_SHELL_KERNEL = mod.get_function("shell_mask_kernel")
    return _GPU_SHELL_KERNEL


@dataclass
class GroupSpec:
    out_base: str
    sel: str
    atom_indices: np.ndarray


@dataclass
class MutableSystemTemplate:
    system: Any
    nb_force: Any
    cnb_force: Any | None
    orig_nb_params: list[tuple[Any, Any, Any]]
    group_offset_idx: list[int]
    solvent_set: set[int]
    active_group: np.ndarray
    cnb_interaction_slot: int | None


def residue_key(res):
    return residue_sort_token_chain(res)


def _frame_indices(total_frames: int, dt_ps: float, tail_ns: float, skip: int) -> Iterable[int]:
    frames_to_load = int((tail_ns * 1e3) / dt_ps)
    first_idx = max(0, total_frames - frames_to_load)
    return range(first_idx, total_frames, skip)


def _parse_box_vectors_from_pdb(pdb_path: str):
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("CRYST1"):
                continue
            fields = line.split()
            a = float(fields[1]) * 0.1
            b = float(fields[2]) * 0.1
            c = float(fields[3]) * 0.1
            alpha, beta, gamma = float(fields[4]), float(fields[5]), float(fields[6])
            if not (abs(alpha - 90.0) < 1e-3 and abs(beta - 90.0) < 1e-3 and abs(gamma - 90.0) < 1e-3):
                raise ValueError(
                    "Box vectors are not orthogonal. Orthogonal vectors are required."
                )
            return (Vec3(a, 0, 0), Vec3(0, b, 0), Vec3(0, 0, c))
    raise ValueError("Cannot find CRYST1 in PDB file.")


def _resolve_residue_selection(u_analysis, potmod, chain: str | None, token: str) -> str:
    _, resseq, icode = potmod.residue_selection(None, token)
    resid_clause = f"{resseq}{icode}"
    tried = []
    if chain:
        sel_try = f"chainID {chain} and resid {resid_clause}"
        tried.append(sel_try)
        if len(u_analysis.select_atoms(sel_try)) > 0:
            return sel_try
    else:
        sel_any_chain = f"resid {resid_clause}"
        tried.append(sel_any_chain)
        if len(u_analysis.select_atoms(sel_any_chain)) > 0:
            return sel_any_chain
    raise RuntimeError(f"Your atom selection picked zero atoms – tried {tried}.")


def build_group_specs(args, protein_name: str, u_analysis, potmod) -> list[GroupSpec]:
    groups: list[GroupSpec] = []

    if args.selection:
        sel = args.selection
        out_base = f"{protein_name}_{sel.replace(' ', '_')}"
        atoms = u_analysis.select_atoms(sel)
        if atoms.n_atoms == 0:
            raise RuntimeError(f"Your atom selection picked zero atoms: '{sel}'.")
        groups.append(GroupSpec(out_base=out_base, sel=sel, atom_indices=np.asarray(atoms.indices, dtype=np.int32)))
        return groups

    if args.groupsFile:
        with open(args.groupsFile, encoding="utf-8") as fh:
            sel_lines = [line.strip() for line in fh if line.strip() and not line.startswith("#")]
        for i, sel in enumerate(sel_lines, start=1):
            out_base = f"{protein_name}_group{i}"
            atoms = u_analysis.select_atoms(sel)
            if atoms.n_atoms == 0:
                raise RuntimeError(f"Your atom selection picked zero atoms: '{sel}'.")
            groups.append(GroupSpec(out_base=out_base, sel=sel, atom_indices=np.asarray(atoms.indices, dtype=np.int32)))
        return groups

    protein_res = u_analysis.select_atoms("protein").residues
    residues = sorted(protein_res if len(protein_res) > 0 else u_analysis.residues, key=lambda r: residue_key(r)[0])
    if args.maxResidues is not None:
        residues = residues[: args.maxResidues]

    for res in residues:
        _, token, chain, _ = residue_key(res)
        if args.multiChain:
            sel = _resolve_residue_selection(u_analysis, potmod, chain, token)
            if chain:
                out_base = f"{protein_name}_res{token}_chain{chain}"
            else:
                out_base = f"{protein_name}_res{token}"
        else:
            sel = _resolve_residue_selection(u_analysis, potmod, None, token)
            out_base = f"{protein_name}_res{token}"
        atoms = u_analysis.select_atoms(sel)
        if atoms.n_atoms == 0:
            raise RuntimeError(f"Your atom selection picked zero atoms: '{sel}'.")
        groups.append(GroupSpec(out_base=out_base, sel=sel, atom_indices=np.asarray(atoms.indices, dtype=np.int32)))

    return groups


def _build_group_batches(groups: list[GroupSpec], batch_size: int):
    batches = []
    for start in range(0, len(groups), batch_size):
        end = min(start + batch_size, len(groups))
        sub = groups[start:end]
        offsets = np.zeros((len(sub) + 1,), dtype=np.int32)
        total_atoms = int(sum(int(g.atom_indices.size) for g in sub))
        atom_indices = np.empty((total_atoms,), dtype=np.int32)
        off = 0
        for i, g in enumerate(sub):
            offsets[i] = off
            n = int(g.atom_indices.size)
            atom_indices[off : off + n] = g.atom_indices
            off += n
        offsets[len(sub)] = off
        batches.append((start, end, atom_indices, offsets))
    return batches


def _select_platform(use_cpu: bool):
    if use_cpu:
        return Platform.getPlatformByName("CPU")
    try:
        return Platform.getPlatformByName("CUDA")
    except Exception:
        warnings.warn("CUDA unavailable; using CPU platform.")
        return Platform.getPlatformByName("CPU")


def _build_mutable_system_template(system, cutoff_nm: float, solvent_idx: list[int]) -> MutableSystemTemplate:
    # Build one base System and keep it mutable between context creations.
    nb_force = None
    cnb_force = None
    for force in system.getForces():
        name = force.__class__.__name__
        if name == "NonbondedForce":
            force.setNonbondedMethod(force.CutoffPeriodic)
            force.setCutoffDistance(cutoff_nm * nanometer)
            nb_force = force
        elif name == "CustomNonbondedForce":
            force.setNonbondedMethod(force.CutoffPeriodic)
            force.setCutoffDistance(cutoff_nm * nanometer)
            cnb_force = force

    if nb_force is None:
        raise RuntimeError("Failed to locate NonbondedForce in the OpenMM System.")

    solvent_set = {int(i) for i in solvent_idx}

    nb_force.setForceGroup(0)
    nb_force.addGlobalParameter("group_scale", 1.0)
    nb_force.addGlobalParameter("solvent_scale", 1.0)

    n_particles = nb_force.getNumParticles()
    orig_nb_params: list[tuple[Any, Any, Any]] = [nb_force.getParticleParameters(i) for i in range(n_particles)]
    group_offset_idx: list[int] = [0] * n_particles

    # Add one mutable group offset entry for every particle; most remain zero.
    for i in range(n_particles):
        group_offset_idx[i] = int(nb_force.addParticleParameterOffset("group_scale", i, 0, 0, 0))

    # Solvent offsets are fixed across all groups.
    for i in solvent_set:
        q, sigma, eps = orig_nb_params[i]
        nb_force.addParticleParameterOffset("solvent_scale", int(i), q, sigma, eps)

    # Base parameters are zero; active params come only through offsets.
    for i in range(n_particles):
        nb_force.setParticleParameters(i, 0, 0, 0)

    for i in range(nb_force.getNumExceptions()):
        p1, p2, qprod, sigma, eps = nb_force.getExceptionParameters(i)
        nb_force.setExceptionParameters(i, p1, p2, 0, 0, 0)

    cnb_interaction_slot = None
    for force in system.getForces():
        name = force.__class__.__name__
        if name == "CustomNonbondedForce":
            force.setForceGroup(1)
            if force.getNumInteractionGroups() == 0:
                # Placeholder that will be overwritten before each context creation.
                force.addInteractionGroup({0}, {0})
            cnb_interaction_slot = 0
        elif name != "NonbondedForce":
            force.setForceGroup(2)

    return MutableSystemTemplate(
        system=system,
        nb_force=nb_force,
        cnb_force=cnb_force,
        orig_nb_params=orig_nb_params,
        group_offset_idx=group_offset_idx,
        solvent_set=solvent_set,
        active_group=np.empty((0,), dtype=np.int32),
        cnb_interaction_slot=cnb_interaction_slot,
    )


def _activate_group_on_template(template: MutableSystemTemplate, atom_indices: np.ndarray) -> None:
    # Clear previous active residue offsets.
    for idx in template.active_group:
        i = int(idx)
        off = template.group_offset_idx[i]
        template.nb_force.setParticleParameterOffset(off, "group_scale", i, 0, 0, 0)

    # Set new residue offsets.
    for idx in atom_indices:
        i = int(idx)
        q, sigma, eps = template.orig_nb_params[i]
        off = template.group_offset_idx[i]
        template.nb_force.setParticleParameterOffset(off, "group_scale", i, q, sigma, eps)

    if template.cnb_force is not None and template.cnb_interaction_slot is not None:
        template.cnb_force.setInteractionGroupParameters(
            template.cnb_interaction_slot,
            {int(i) for i in atom_indices},
            template.solvent_set,
        )

    template.active_group = np.asarray(atom_indices, dtype=np.int32)


def run_inproc(args, processed_pdb_path: str):
    from hydromap.potentials import potential as potmod

    u = mda.Universe(processed_pdb_path, args.trajectory)
    protein_name = os.path.splitext(os.path.basename(processed_pdb_path))[0].replace("_processed", "")
    groups = build_group_specs(args, protein_name, u, potmod)
    if not groups:
        raise RuntimeError("No groups/residues were built for analysis.")

    os.makedirs(args.outdir, exist_ok=True)
    output_paths = [Path(args.outdir) / f"{g.out_base}_potentials.csv" for g in groups]

    all_waters = u.select_atoms(NORMALIZED_WATER_SELECTION)
    water_ox = u.select_atoms(NORMALIZED_WATER_OXYGEN_SELECTION)
    solvent_idx = all_waters.indices.tolist()
    frame_indices = list(_frame_indices(len(u.trajectory), u.trajectory.dt, args.time, args.skip))

    _parse_box_vectors_from_pdb(processed_pdb_path)
    base_system = load_openmm_system(args.top, processed_pdb_path)

    platform = _select_platform(args.nogpu)
    use_gpu_shell = False
    cp = None
    shell_kernel = None
    mask_buf = None
    pool_peak_bytes = 0
    shell_peak_mask_bytes = 0
    if platform.getName() == "CUDA":
        try:
            cp = _get_cupy()
            shell_kernel = _get_gpu_shell_kernel(cp)
            use_gpu_shell = True
        except Exception:
            warnings.warn("CuPy unavailable for shell counting; falling back to MDAnalysis around selections.")

    if platform.getName() == "CUDA":
        if args.groupBatchSize <= 0:
            batch_size = min(DEFAULT_GPU_GROUP_BATCH_SIZE, len(groups))
        else:
            requested = max(1, args.groupBatchSize)
            if requested > MAX_GPU_GROUP_BATCH_SIZE:
                warnings.warn(
                    f"--groupBatchSize={requested} is large for CUDA and may be slower/unstable; "
                    f"clamping to {MAX_GPU_GROUP_BATCH_SIZE}."
                )
                requested = MAX_GPU_GROUP_BATCH_SIZE
            batch_size = min(requested, len(groups))
    else:
        batch_size = len(groups) if args.groupBatchSize <= 0 else max(1, min(args.groupBatchSize, len(groups)))
    batches = _build_group_batches(groups, batch_size)

    t_ctx_build = 0.0
    t_setpos = 0.0
    t_energy = 0.0
    t_nwaters = 0.0
    kernel_calls = 0
    used_mutable_template = False

    mutable_template = None
    try:
        mutable_template = _build_mutable_system_template(base_system, cutoff_nm=args.cutoff / 10.0, solvent_idx=solvent_idx)
        used_mutable_template = True
    except Exception as exc:
        warnings.warn(f"Mutable system template unavailable; falling back to legacy per-group build. ({exc})")
        mutable_template = None

    for b_start, b_end, batch_atom_idx, batch_offsets in batches:
        batch_groups = groups[b_start:b_end]

        contexts = []
        integrators = []
        t0 = time.perf_counter()
        for g in batch_groups:
            if mutable_template is not None:
                _activate_group_on_template(mutable_template, g.atom_indices)
                system = mutable_template.system
            else:
                system = XmlSerializer.deserialize(XmlSerializer.serialize(base_system))
                system = potmod.prepare_system(system, g.atom_indices.tolist(), solvent_idx, cutoff_nm=args.cutoff / 10.0)
            integ = VerletIntegrator(0.001 * picosecond)
            ctx = Context(system, integ, platform)
            contexts.append(ctx)
            integrators.append(integ)
        t_ctx_build += time.perf_counter() - t0

        coul_data = [[] for _ in batch_groups]
        lj_data = [[] for _ in batch_groups]
        nwat_data = [[] for _ in batch_groups]

        if use_gpu_shell:
            offsets_cp = cp.asarray(batch_offsets, dtype=cp.int32)
            if shell_peak_mask_bytes < int(offsets_cp.nbytes):
                shell_peak_mask_bytes = int(offsets_cp.nbytes)

        for frame_idx in tqdm(frame_indices, total=len(frame_indices), unit="frame", leave=False):
            u.trajectory[frame_idx]
            pos_nm = u.atoms.positions / 10.0
            t1 = time.perf_counter()
            pos_quantity = pos_nm * nanometer
            for ctx in contexts:
                ctx.setPositions(pos_quantity)
            t_setpos += time.perf_counter() - t1

            if use_gpu_shell:
                t2 = time.perf_counter()
                n_w = int(water_ox.n_atoms)
                if n_w > 0:
                    grp_pos_gpu = cp.asarray(u.atoms.positions[batch_atom_idx], dtype=cp.float64)
                    wat_pos_gpu = cp.asarray(water_ox.positions, dtype=cp.float64)
                    box_gpu = cp.asarray(u.dimensions[:3].copy(), dtype=cp.float64)

                    total_mask = len(batch_groups) * n_w
                    if mask_buf is None or int(mask_buf.size) < total_mask:
                        mask_buf = cp.empty((total_mask,), dtype=cp.uint8)
                    mask_view = mask_buf[:total_mask]
                    blocks = (total_mask + 255) // 256
                    shell_kernel(
                        (blocks,),
                        (256,),
                        (
                            grp_pos_gpu,
                            offsets_cp,
                            np.int32(len(batch_groups)),
                            wat_pos_gpu,
                            np.int32(n_w),
                            box_gpu,
                            np.float64(args.cutoff * args.cutoff),
                            mask_view,
                        ),
                    )
                    counts = cp.asnumpy(cp.sum(mask_view.reshape((len(batch_groups), n_w)), axis=1, dtype=cp.int64)).astype(int)
                    kernel_calls += 1
                    shell_peak_mask_bytes = max(shell_peak_mask_bytes, int(mask_view.nbytes))
                    pool_peak_bytes = max(pool_peak_bytes, int(cp.get_default_memory_pool().used_bytes()))
                else:
                    counts = np.zeros((len(batch_groups),), dtype=int)
                t_nwaters += time.perf_counter() - t2
            else:
                t2 = time.perf_counter()
                counts = np.zeros((len(batch_groups),), dtype=int)
                for i, g in enumerate(batch_groups):
                    shell = u.select_atoms(f"{NORMALIZED_WATER_OXYGEN_SELECTION} and around {args.cutoff} ({g.sel})")
                    counts[i] = int(shell.n_atoms)
                t_nwaters += time.perf_counter() - t2

            solv_ce = None
            if contexts:
                t3 = time.perf_counter()
                solv_ce = potmod.coulomb_energy(contexts[0], 0, 1)
                t_energy += time.perf_counter() - t3

            for i, ctx in enumerate(contexts):
                t3 = time.perf_counter()
                total = potmod.coulomb_energy(ctx, 1, 1)
                group_ce = potmod.coulomb_energy(ctx, 1, 0)
                coul = (total - group_ce - solv_ce).value_in_unit(kilojoule_per_mole)
                lj = ctx.getState(getEnergy=True, groups={1}).getPotentialEnergy().value_in_unit(kilojoule_per_mole)
                t_energy += time.perf_counter() - t3

                coul_data[i].append(float(coul))
                lj_data[i].append(float(lj))
                nwat_data[i].append(int(counts[i]))

        for i, g in enumerate(batch_groups):
            c = np.asarray(coul_data[i], dtype=float)
            l = np.asarray(lj_data[i], dtype=float)
            n = np.asarray(nwat_data[i], dtype=int)
            df = pd.DataFrame({"coulomb": c, "lj": l, "total": c + l, "n_waters": n})
            df.to_csv(output_paths[b_start + i], index=False)

        for ctx in contexts:
            del ctx
        for integ in integrators:
            del integ

    print(f"Done! Wrote {len(groups)} potential files to {args.outdir}")
    print(f"Frames processed: {len(frame_indices)}")
    print(f"Groups processed: {len(groups)}")
    print(f"Group batch size: {batch_size} ({len(batches)} batch(es))")
    print(f"Context build total: {t_ctx_build:.6f}s")
    print(f"SetPositions total: {t_setpos:.6f}s")
    print(f"Energy query total: {t_energy:.6f}s")
    print(f"n_waters total: {t_nwaters:.6f}s")
    print(f"Mutable system template: {'enabled' if used_mutable_template else 'disabled'}")
    if use_gpu_shell:
        print(f"GPU shell kernel calls: {kernel_calls}")
        print(f"GPU shell mask peak bytes: {shell_peak_mask_bytes}")
        print(f"CuPy pool peak used bytes: {pool_peak_bytes}")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Run production potential analysis serially in-process (GPU/CPU)."
    )
    parser.add_argument("protein", help="input PDB (unprocessed)")
    parser.add_argument("trajectory", help="trajectory file (e.g. traj.dcd)")
    parser.add_argument(
        "--top",
        required=True,
        help="Topology for the processed PDB: OpenMM System XML or GROMACS .top.",
    )
    parser.add_argument("--multiChain", action="store_true", help="Protein has multiple chains")
    parser.add_argument("--groupsFile", help="File with MDAnalysis selection strings, one per line")
    parser.add_argument("--selection", help="One custom MDAnalysis selection string")
    parser.add_argument("-t", "--time", type=float, default=5.0, help="Last X ns to analyse (default 5)")
    parser.add_argument("--skip", type=int, default=100, help="Frame stride (default 100)")
    parser.add_argument("--cutoff", type=float, default=4.25, help="Cutoff distance in Angstroms (default 4.25)")
    parser.add_argument("--outdir", type=str, default="potentials", help="Output directory for CSVs (default potentials)")
    parser.add_argument("--nogpu", action="store_true", default=False, help="Force CPU platform")
    parser.add_argument(
        "--engine",
        default="auto",
        choices=["auto", "inproc"],
        help="Execution engine (default auto -> inproc).",
    )
    parser.add_argument(
        "--groupBatchSize",
        type=int,
        default=0,
        help=(
            "Groups per in-process batch. 0 means auto: CUDA uses a safe default "
            f"({DEFAULT_GPU_GROUP_BATCH_SIZE}), CPU uses all groups."
        ),
    )
    parser.add_argument(
        "--maxResidues",
        type=int,
        default=None,
        help="Optional cap on number of residues (for benchmarking/debugging)",
    )

    args = parser.parse_args(argv)

    if args.multiChain and args.groupsFile:
        sys.exit("ERROR: use EITHER --multiChain OR --groupsFile, not both")
    if args.selection and args.groupsFile:
        sys.exit("ERROR: use EITHER --selection OR --groupsFile, not both")
    if args.selection and args.multiChain:
        sys.exit("ERROR: --selection is mutually exclusive with --multiChain")
    if args.groupBatchSize < 0:
        sys.exit("ERROR: --groupBatchSize must be >= 0.")

    pdb_path = args.protein if args.protein.endswith(".pdb") else args.protein + ".pdb"
    if "processed" in pdb_path:
        sys.exit(
            "ERROR: input protein appears processed. Please provide unprocessed PDB path "
            "(the script will use *_processed.pdb for analysis)."
        )
    if not os.path.exists(pdb_path):
        sys.exit(f"ERROR: cannot find {pdb_path}")

    processed_pdb_path = pdb_path[:-4] + "_processed.pdb"
    if not os.path.exists(processed_pdb_path):
        sys.exit(f"ERROR: cannot find processed PDB {processed_pdb_path}")
    if not os.path.exists(args.trajectory):
        sys.exit(f"ERROR: cannot find trajectory {args.trajectory}")
    if not os.path.exists(args.top):
        sys.exit(f"ERROR: cannot find topology {args.top}")

    engine = "inproc" if args.engine == "auto" else args.engine
    start = time.time()
    print(f"Running in-process potential analysis with engine={engine}, nogpu={args.nogpu}")
    run_inproc(args, processed_pdb_path)
    print(f"Elapsed wall time: {time.time() - start:.2f}s")


if __name__ == "__main__":
    main()
