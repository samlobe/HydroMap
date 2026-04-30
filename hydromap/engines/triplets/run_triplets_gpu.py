#!/usr/bin/env python3
"""
Run many triplet jobs serially on one process.

Production engine:
- batch_groups: batches all groups' shell centers per frame into one
  triplet kernel call (CPU or GPU), then demultiplexes outputs.

Shell construction modes:
- mda: MDAnalysis updating selections.
- gpu: CUDA kernel builds all group shell masks per frame (batch_groups + GPU only).
"""

import argparse
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import MDAnalysis as mda
import numpy as np
from tqdm import tqdm
from hydromap.analysis_inputs import NORMALIZED_WATER_OXYGEN_SELECTION
from hydromap.selection import residue_sort_token_chain

_GPU_SHELL_KERNEL = None
_GPU_SHELL_CUDA_SRC = r"""
extern "C" __global__
void shell_mask_kernel(
    const double* group_atoms,
    const int* group_offsets,
    int n_groups,
    const double* waters,
    int n_waters,
    const double* box,
    double low2,
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
        if (r2 < min_r2) {
            min_r2 = r2;
            if (low2 > 0.0 && min_r2 <= low2) {
                break;
            }
        }
    }

    out_mask[tid] = (min_r2 <= high2 && min_r2 > low2) ? 1 : 0;
}
"""


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
    shell_sel_expr: str
    shell_ag: Optional[object] = None
    target_atom_indices: Optional[np.ndarray] = None


def residue_key(res):
    """
    Deterministic key/token for residues with insertion codes.
    Returns (sort_tuple, res_token, chain_token, alternate_chain_tokens)
    """
    return residue_sort_token_chain(res)


def _resolve_residue_selection(u_analysis, tri, chain: Optional[str], token: str) -> str:
    _, resseq, icode = tri.residue_selection(None, token, heavy_only=True)
    resid_clause = f"{resseq}{icode}"
    tried = []
    if chain:
        sel_try = f"chainID {chain} and resid {resid_clause} and not name H*"
        tried.append(sel_try)
        if len(u_analysis.select_atoms(sel_try)) > 0:
            return sel_try
    else:
        sel_any_chain = f"resid {resid_clause} and not name H*"
        tried.append(sel_any_chain)
        if len(u_analysis.select_atoms(sel_any_chain)) > 0:
            return sel_any_chain
    raise RuntimeError(f"Your atom selection picked zero atoms – tried {tried}.")


def build_group_specs(args, protein_name: str, u_analysis, tri) -> list[GroupSpec]:
    waters_sel = NORMALIZED_WATER_OXYGEN_SELECTION
    groups: list[GroupSpec] = []

    def with_cutoff_suffix(base: str) -> str:
        if args.hydrationLowCutoff is None:
            return base
        return f"{base}_lowC_{args.hydrationLowCutoff}_highC_{args.hydrationCutoff}"

    if args.selection:
        out_base = with_cutoff_suffix(f"{protein_name}_{args.selection.replace(' ', '_')}")
        sel = args.selection
        if args.hydrationLowCutoff is not None:
            shell_expr = (
                f"({waters_sel}) and around {args.hydrationCutoff} ({sel}) "
                f"and not around {args.hydrationLowCutoff} ({sel})"
            )
        else:
            shell_expr = f"({waters_sel}) and around {args.hydrationCutoff} ({sel})"
        groups.append(GroupSpec(out_base=out_base, sel=sel, shell_sel_expr=shell_expr))
        return groups

    if args.groupsFile:
        with open(args.groupsFile) as fh:
            sel_lines = [line.strip() for line in fh if line.strip() and not line.startswith("#")]
        for i, sel in enumerate(sel_lines, start=1):
            out_base = with_cutoff_suffix(f"{protein_name}_group{i}")
            if args.hydrationLowCutoff is not None:
                shell_expr = (
                    f"({waters_sel}) and around {args.hydrationCutoff} ({sel}) "
                    f"and not around {args.hydrationLowCutoff} ({sel})"
                )
            else:
                shell_expr = f"({waters_sel}) and around {args.hydrationCutoff} ({sel})"
            groups.append(GroupSpec(out_base=out_base, sel=sel, shell_sel_expr=shell_expr))
        return groups

    # default: per-residue groups
    residues = sorted(u_analysis.select_atoms("protein").residues, key=lambda r: residue_key(r)[0])
    if args.maxResidues is not None:
        residues = residues[: args.maxResidues]

    for res in residues:
        _, token, chain, _ = residue_key(res)
        if args.multiChain:
            sel = _resolve_residue_selection(u_analysis, tri, chain, token)
            if chain:
                out_base = with_cutoff_suffix(f"{protein_name}_res{token}_chain{chain}")
            else:
                out_base = with_cutoff_suffix(f"{protein_name}_res{token}")
        else:
            sel = _resolve_residue_selection(u_analysis, tri, None, token)
            out_base = with_cutoff_suffix(f"{protein_name}_res{token}")

        if args.hydrationLowCutoff is not None:
            shell_expr = (
                f"({waters_sel}) and around {args.hydrationCutoff} ({sel}) "
                f"and not around {args.hydrationLowCutoff} ({sel})"
            )
        else:
            shell_expr = f"({waters_sel}) and around {args.hydrationCutoff} ({sel})"
        groups.append(GroupSpec(out_base=out_base, sel=sel, shell_sel_expr=shell_expr))

    return groups


def _frame_indices(total_frames: int, dt_ps: float, tail_ns: float, skip: int) -> Iterable[int]:
    frames_to_load = int((tail_ns * 1e3) / dt_ps)
    first_idx = max(0, total_frames - frames_to_load)
    return range(first_idx, total_frames, skip)


def run_inproc(args, processed_pdb_path: str):
    from hydromap.triplets import triplet as tri

    u = mda.Universe(processed_pdb_path, args.trajectory)
    total_frames = len(u.trajectory)

    protein_name = os.path.splitext(os.path.basename(processed_pdb_path))[0].replace("_processed", "")
    groups = build_group_specs(args, protein_name, u, tri)
    if not groups:
        raise RuntimeError("No groups/residues were built for analysis.")

    os.makedirs(args.outdir, exist_ok=True)

    output_paths = [Path(args.outdir) / f"{g.out_base}_angles.txt" for g in groups]
    out_handles = [open(p, "w", encoding="utf-8") for p in output_paths]

    all_waters = u.select_atoms(NORMALIZED_WATER_OXYGEN_SELECTION)
    frame_indices = list(_frame_indices(total_frames, u.trajectory.dt, args.time, args.skip))

    low_cut = 0.0
    high_cut = 3.5
    shell_low_cut = 0.0 if args.hydrationLowCutoff is None else float(args.hydrationLowCutoff)
    shell_high_cut = float(args.hydrationCutoff)
    shell_low2 = shell_low_cut * shell_low_cut
    shell_high2 = shell_high_cut * shell_high_cut

    t_shell = 0.0
    t_kernel = 0.0
    kernel_calls = 0
    group_batch_size = len(groups) if args.groupBatchSize is None or args.groupBatchSize <= 0 else int(args.groupBatchSize)
    group_batch_size = max(1, min(group_batch_size, len(groups)))
    batch_ranges = [(i, min(i + group_batch_size, len(groups))) for i in range(0, len(groups), group_batch_size)]

    use_gpu_shell = args.shellMode == "gpu"
    cp = None
    shell_kernel = None
    gpu_shell_batches = []
    shell_mask_cp = None
    gpu_peak_pool_bytes = 0
    shell_peak_mask_bytes = 0

    if use_gpu_shell:
        cp = tri._get_cupy()
        shell_kernel = _get_gpu_shell_kernel(cp)
        mempool = cp.get_default_memory_pool()

        group_atom_arrays = []
        for g in groups:
            ag = u.select_atoms(g.sel)
            idx = np.asarray(ag.indices, dtype=np.int32)
            if idx.size == 0:
                raise RuntimeError(f"Your atom selection picked zero atoms: '{g.sel}'.")
            g.target_atom_indices = idx
            group_atom_arrays.append(idx)

        for start, end in batch_ranges:
            batch_atom_arrays = group_atom_arrays[start:end]
            batch_n_groups = end - start
            batch_offsets = np.zeros((batch_n_groups + 1,), dtype=np.int32)
            total_batch_atoms = int(sum(int(x.size) for x in batch_atom_arrays))
            batch_atom_indices = np.empty((total_batch_atoms,), dtype=np.int32)
            off = 0
            for i, idx in enumerate(batch_atom_arrays):
                batch_offsets[i] = off
                n = int(idx.size)
                batch_atom_indices[off : off + n] = idx
                off += n
            batch_offsets[batch_n_groups] = off
            gpu_shell_batches.append(
                {
                    "start": start,
                    "end": end,
                    "n_groups": batch_n_groups,
                    "atom_indices": batch_atom_indices,
                    "offsets_cp": cp.asarray(batch_offsets, dtype=cp.int32),
                }
            )
    else:
        for g in groups:
            g.shell_ag = u.select_atoms(g.shell_sel_expr, updating=True)

    try:
        for write_idx, frame_idx in tqdm(enumerate(frame_indices), total=len(frame_indices), unit="frame"):
            u.trajectory[frame_idx]
            box = u.dimensions[:3].copy()
            all_water_pos = all_waters.positions

            center_arrays = None
            if not use_gpu_shell:
                t0 = time.perf_counter()
                center_arrays = [g.shell_ag.positions for g in groups]
                t_shell += time.perf_counter() - t0

            if use_gpu_shell:
                n_waters = int(all_water_pos.shape[0])
                all_water_pos_gpu = cp.asarray(all_water_pos, dtype=cp.float64)
                box_gpu = cp.asarray(box, dtype=cp.float64)
                if int(all_water_pos_gpu.nbytes) > shell_peak_mask_bytes:
                    shell_peak_mask_bytes = int(all_water_pos_gpu.nbytes)

                for batch in gpu_shell_batches:
                    b_start = int(batch["start"])
                    b_end = int(batch["end"])
                    b_ngroups = int(batch["n_groups"])

                    if n_waters <= 0:
                        for gi in range(b_start, b_end):
                            out_handles[gi].write("\n")
                        continue

                    t_shell_batch = time.perf_counter()
                    group_atom_pos = u.atoms.positions[batch["atom_indices"]]
                    group_atom_pos_gpu = cp.asarray(group_atom_pos, dtype=cp.float64)
                    total_mask = b_ngroups * n_waters
                    if shell_mask_cp is None or int(shell_mask_cp.size) < total_mask:
                        shell_mask_cp = cp.empty((total_mask,), dtype=cp.uint8)
                    shell_mask_view = shell_mask_cp[:total_mask]
                    threads = 256
                    blocks = (total_mask + threads - 1) // threads
                    shell_kernel(
                        (blocks,),
                        (threads,),
                        (
                            group_atom_pos_gpu,
                            batch["offsets_cp"],
                            np.int32(b_ngroups),
                            all_water_pos_gpu,
                            np.int32(n_waters),
                            box_gpu,
                            np.float64(shell_low2),
                            np.float64(shell_high2),
                            shell_mask_view,
                        ),
                    )
                    center_counts = cp.asnumpy(
                        cp.sum(shell_mask_view.reshape((b_ngroups, n_waters)), axis=1, dtype=cp.int64)
                    ).astype(int)
                    t_shell += time.perf_counter() - t_shell_batch
                    shell_peak_mask_bytes = max(shell_peak_mask_bytes, int(shell_mask_view.nbytes))

                    total_centers = int(np.sum(center_counts, dtype=np.int64))
                    if total_centers > 0:
                        selected_linear = cp.where(shell_mask_view > 0)[0]
                        selected_water_idx = selected_linear % n_waters
                        sub_all = all_water_pos_gpu[selected_water_idx]

                        t1 = time.perf_counter()
                        ang_all, counts_all = tri.triplet_angles_gpu(
                            sub_all,
                            all_water_pos_gpu,
                            box_gpu,
                            low_cut,
                            high_cut,
                        )
                        t_kernel += time.perf_counter() - t1
                        kernel_calls += 1
                    else:
                        ang_all = np.array([], dtype=float)
                        counts_all = np.array([], dtype=int)

                    center_off = 0
                    angle_off = 0
                    for i_rel in range(b_ngroups):
                        gi = b_start + i_rel
                        n_centers = int(center_counts[i_rel])
                        if n_centers <= 0:
                            out_handles[gi].write("\n")
                            continue
                        group_counts = counts_all[center_off : center_off + n_centers]
                        n_angles = int(np.sum(group_counts, dtype=np.int64))
                        group_angles = ang_all[angle_off : angle_off + n_angles]
                        out_handles[gi].write(" ".join(f"{a:.1f}" for a in group_angles) + "\n")
                        center_off += n_centers
                        angle_off += n_angles

                    gpu_peak_pool_bytes = max(gpu_peak_pool_bytes, int(mempool.used_bytes()))
            else:
                for b_start, b_end in batch_ranges:
                    batch_centers = center_arrays[b_start:b_end]
                    non_empty = [a for a in batch_centers if len(a) > 0]
                    if non_empty:
                        sub_all = np.concatenate(non_empty, axis=0)
                        t1 = time.perf_counter()
                        if args.gpu:
                            ang_all, counts_all = tri.triplet_angles_gpu(
                                sub_all,
                                all_water_pos,
                                box,
                                low_cut,
                                high_cut,
                            )
                        else:
                            ang_all, counts_all = tri.wl.triplet_angles(
                                sub_all,
                                all_water_pos,
                                box,
                                low_cut,
                                high_cut,
                            )
                        t_kernel += time.perf_counter() - t1
                        kernel_calls += 1
                    else:
                        ang_all = np.array([], dtype=float)
                        counts_all = np.array([], dtype=int)

                    center_off = 0
                    angle_off = 0
                    for i_rel, centers in enumerate(batch_centers):
                        gi = b_start + i_rel
                        n_centers = int(centers.shape[0])
                        if n_centers <= 0:
                            out_handles[gi].write("\n")
                            continue
                        group_counts = counts_all[center_off : center_off + n_centers]
                        n_angles = int(np.sum(group_counts, dtype=np.int64))
                        group_angles = ang_all[angle_off : angle_off + n_angles]
                        out_handles[gi].write(" ".join(f"{a:.1f}" for a in group_angles) + "\n")
                        center_off += n_centers
                        angle_off += n_angles

            if write_idx % 200 == 0:
                for fh in out_handles:
                    fh.flush()
    finally:
        for fh in out_handles:
            fh.close()

    print(f"Done! Wrote {len(groups)} angle files to {args.outdir}")
    print(f"Frames processed: {len(frame_indices)}")
    print(f"Groups processed: {len(groups)}")
    print(f"Group batch size: {group_batch_size} ({len(batch_ranges)} batch(es)/frame)")
    print(f"Shell extraction total: {t_shell:.6f}s ({t_shell / max(1, len(frame_indices)):.6e} s/frame)")
    print(f"Triplet kernel total: {t_kernel:.6f}s")
    print(f"Triplet kernel calls: {kernel_calls}")
    if use_gpu_shell:
        print(f"GPU shell mask peak bytes: {shell_peak_mask_bytes}")
        print(f"CuPy pool peak used bytes: {gpu_peak_pool_bytes}")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Run production triplet analysis serially in-process with batched groups."
        )
    )
    parser.add_argument("protein", help="input PDB (unprocessed)")
    parser.add_argument("trajectory", help="trajectory file (e.g. traj.dcd)")
    parser.add_argument("--multiChain", action="store_true", help="Protein has multiple chains")
    parser.add_argument("--groupsFile", help="File with MDAnalysis selection strings, one per line")
    parser.add_argument("--selection", help="One custom MDAnalysis selection string")
    parser.add_argument("-t", "--time", type=float, default=5.0, help="Last X ns to analyse (default 5)")
    parser.add_argument("--hydrationCutoff", type=float, default=4.25, help="Hydration cutoff in Angstroms (default 4.25)")
    parser.add_argument("--hydrationLowCutoff", type=float, default=None, help="Optional lower shell cutoff in Angstroms")
    parser.add_argument("--skip", type=int, default=1, help="Frame stride (default 1)")
    parser.add_argument("--outdir", type=str, default="angles", help="Output directory for angles files (default angles)")
    parser.add_argument("--gpu", action="store_true", default=False, help="Use GPU triplet kernel")
    parser.add_argument(
        "--shellMode",
        default="auto",
        choices=["auto", "mda", "gpu"],
        help=(
            "Shell construction backend. "
            "'mda': MDAnalysis updating selections; "
            "'gpu': CUDA shell mask kernel (batch_groups + --gpu only); "
            "'auto': picks gpu when using batch_groups+gpu."
        ),
    )
    parser.add_argument(
        "--engine",
        default="auto",
        choices=["auto", "batch_groups"],
        help="Execution engine (default auto -> batch_groups).",
    )
    parser.add_argument(
        "--groupBatchSize",
        type=int,
        default=0,
        help=(
            "Number of groups per batched triplet call in batch_groups mode. "
            "0 means all groups in one batch (fastest, highest memory)."
        ),
    )
    parser.add_argument(
        "--maxResidues",
        type=int,
        default=None,
        help="Optional cap on number of residues (for quick benchmarking/debugging)",
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
            "ERROR: input protein appears processed. Please provide the unprocessed PDB path "
            "(the script will use *_processed.pdb for analysis)."
        )
    if not os.path.exists(pdb_path):
        sys.exit(f"ERROR: cannot find {pdb_path}")

    processed_pdb_path = pdb_path[:-4] + "_processed.pdb"
    if not os.path.exists(processed_pdb_path):
        sys.exit(f"ERROR: cannot find processed PDB {processed_pdb_path}")
    if not os.path.exists(args.trajectory):
        sys.exit(f"ERROR: cannot find trajectory {args.trajectory}")

    engine = args.engine
    if engine == "auto":
        engine = "batch_groups"

    shell_mode = args.shellMode
    if shell_mode == "auto":
        shell_mode = "gpu" if (engine == "batch_groups" and args.gpu) else "mda"
    args.shellMode = shell_mode

    if args.shellMode == "gpu" and not args.gpu:
        sys.exit("ERROR: --shellMode gpu requires --gpu.")

    start = time.time()
    print(
        f"Running in-process triplet analysis with engine={engine}, "
        f"gpu={args.gpu}, shellMode={args.shellMode}"
    )
    run_inproc(args, processed_pdb_path)
    print(f"Elapsed wall time: {time.time() - start:.2f}s")


if __name__ == "__main__":
    main()
