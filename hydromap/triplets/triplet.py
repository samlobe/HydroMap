#!/usr/bin/env python3

import numpy as np
import MDAnalysis as mda
from tqdm import tqdm
import sys, os, argparse, shutil
from glob import glob
import re
import time
import warnings
from hydromap.analysis_inputs import NORMALIZED_WATER_OXYGEN_SELECTION
warnings.filterwarnings("ignore", category=UserWarning, module=r".*MDAnalysis\.topology\.PDBParser")
warnings.filterwarnings("ignore", category=UserWarning, module=r".*MDAnalysis\.topology\.guessers")
warnings.filterwarnings("ignore", category=DeprecationWarning, module=r".*MDAnalysis\.coordinates\.DCD")


try:
    import waterlib as wl
except ImportError:
    try:
        from hydromap.triplets import waterlib as wl
    except ImportError:
        sys.exit(
            "\nERROR: Could not import waterlib.\n"
            "Please compile it first (from the 'hydromap/triplets' directory):\n\n"
            "    python setup.py build_ext --inplace\n\n"
            "Make sure waterlib.c exists in the same directory as triplet.py.\n"
        )


def _configure_cupy_runtime_env() -> None:
    """
    Make CuPy runtime setup more robust on systems where CUDA libs are installed
    in non-default locations and where CUB reduction kernels are incompatible.
    """
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


_CUPY = None
_FUSED_KERNELS = None
_FUSED_MAX_NEIGH = 1024

_FUSED_CUDA_SRC = r'''
#define FUSED_MAX_NEIGH 1024

extern "C" __global__
void count_neighbors(const double* sub, const double* pos, const double* box,
                     int Ns, int Np, double low2, double high2,
                     int* neigh_counts) {
    int i = blockIdx.x;
    if (i >= Ns) return;
    __shared__ int sh[256];
    int t = threadIdx.x;
    int local = 0;

    double cx = sub[3*i+0];
    double cy = sub[3*i+1];
    double cz = sub[3*i+2];

    for (int j = t; j < Np; j += blockDim.x) {
        double dx = pos[3*j+0] - cx;
        double dy = pos[3*j+1] - cy;
        double dz = pos[3*j+2] - cz;
        dx -= box[0] * nearbyint(dx / box[0]);
        dy -= box[1] * nearbyint(dy / box[1]);
        dz -= box[2] * nearbyint(dz / box[2]);
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 > low2 && r2 <= high2) local++;
    }

    sh[t] = local;
    __syncthreads();
    for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
        if (t < stride) sh[t] += sh[t + stride];
        __syncthreads();
    }
    if (t == 0) neigh_counts[i] = sh[0];
}

extern "C" __global__
void gather_neighbors_seq(const double* sub, const double* pos, const double* box,
                          int Ns, int Np, double low2, double high2,
                          const long long* neigh_offsets,
                          double* neigh_vecs) {
    int i = blockIdx.x;
    if (i >= Ns) return;
    if (threadIdx.x != 0) return;

    double cx = sub[3*i+0];
    double cy = sub[3*i+1];
    double cz = sub[3*i+2];
    long long off = neigh_offsets[i];
    int n = 0;

    for (int j = 0; j < Np; ++j) {
        double dx = pos[3*j+0] - cx;
        double dy = pos[3*j+1] - cy;
        double dz = pos[3*j+2] - cz;
        dx -= box[0] * nearbyint(dx / box[0]);
        dy -= box[1] * nearbyint(dy / box[1]);
        dz -= box[2] * nearbyint(dz / box[2]);
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 > low2 && r2 <= high2) {
            long long idx = off + n;
            neigh_vecs[3*idx+0] = dx;
            neigh_vecs[3*idx+1] = dy;
            neigh_vecs[3*idx+2] = dz;
            n++;
        }
    }
}

extern "C" __global__
void angle_pairs_from_neighbors(const double* neigh_vecs,
                                const int* neigh_counts,
                                const long long* neigh_offsets,
                                const long long* pair_offsets,
                                int Ns,
                                double* out_angles) {
    int i = blockIdx.x;
    if (i >= Ns) return;
    int t = threadIdx.x;

    int n = neigh_counts[i];
    if (n < 2) return;
    long long total = ((long long)n * (n - 1)) / 2;
    long long voff = neigh_offsets[i];
    long long poff = pair_offsets[i];

    for (long long p = t; p < total; p += blockDim.x) {
        // triangular index p -> lexicographic pair (a,b), 0<=a<b<n
        double dn = (double)n;
        double disc = (2.0*dn - 1.0)*(2.0*dn - 1.0) - 8.0*(double)p;
        int a = (int)floor((2.0*dn - 1.0 - sqrt(disc)) * 0.5);
        long long before = ((long long)a * (2LL*n - a - 1LL)) / 2LL;
        int b = (int)(p - before + a + 1);

        long long ia = voff + a;
        long long ib = voff + b;

        double ax = neigh_vecs[3*ia+0];
        double ay = neigh_vecs[3*ia+1];
        double az = neigh_vecs[3*ia+2];
        double bx = neigh_vecs[3*ib+0];
        double by = neigh_vecs[3*ib+1];
        double bz = neigh_vecs[3*ib+2];

        double na = sqrt(ax*ax + ay*ay + az*az);
        double nb = sqrt(bx*bx + by*by + bz*bz);
        double ang = 0.0;
        if (na > 0.0 && nb > 0.0) {
            double c = (ax*bx + ay*by + az*bz) / (na*nb);
            c = c > 1.0 ? 1.0 : (c < -1.0 ? -1.0 : c);
            ang = acos(c) * 57.29577951308232;
        }
        out_angles[poff + p] = ang;
    }
}

extern "C" __global__
void direct_angles_sorted(const double* sub, const double* pos, const double* box,
                          int Ns, int Np, double low2, double high2,
                          const long long* pair_offsets,
                          int* pair_counts,
                          int* overflow_flags,
                          double* out_angles) {
    int i = blockIdx.x;
    if (i >= Ns) return;
    int t = threadIdx.x;

    __shared__ int n;
    __shared__ int ids[FUSED_MAX_NEIGH];
    __shared__ double vx[FUSED_MAX_NEIGH];
    __shared__ double vy[FUSED_MAX_NEIGH];
    __shared__ double vz[FUSED_MAX_NEIGH];
    __shared__ double vn[FUSED_MAX_NEIGH];

    if (t == 0) {
        n = 0;
        overflow_flags[i] = 0;
    }
    __syncthreads();

    double cx = sub[3*i+0];
    double cy = sub[3*i+1];
    double cz = sub[3*i+2];

    for (int j = t; j < Np; j += blockDim.x) {
        double dx = pos[3*j+0] - cx;
        double dy = pos[3*j+1] - cy;
        double dz = pos[3*j+2] - cz;
        dx -= box[0] * nearbyint(dx / box[0]);
        dy -= box[1] * nearbyint(dy / box[1]);
        dz -= box[2] * nearbyint(dz / box[2]);
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 > low2 && r2 <= high2) {
            int idx = atomicAdd(&n, 1);
            if (idx < FUSED_MAX_NEIGH) {
                ids[idx] = j;
                vx[idx] = dx;
                vy[idx] = dy;
                vz[idx] = dz;
            }
        }
    }
    __syncthreads();

    if (t == 0) {
        if (n > FUSED_MAX_NEIGH) {
            overflow_flags[i] = 1;
            pair_counts[i] = -1;
        } else {
            // Sort by water index to match CPU's deterministic neighbor order.
            for (int a = 1; a < n; ++a) {
                int key = ids[a];
                double kx = vx[a], ky = vy[a], kz = vz[a];
                int b = a - 1;
                while (b >= 0 && ids[b] > key) {
                    ids[b+1] = ids[b];
                    vx[b+1] = vx[b];
                    vy[b+1] = vy[b];
                    vz[b+1] = vz[b];
                    --b;
                }
                ids[b+1] = key;
                vx[b+1] = kx;
                vy[b+1] = ky;
                vz[b+1] = kz;
            }
            for (int a = 0; a < n; ++a) {
                vn[a] = sqrt(vx[a]*vx[a] + vy[a]*vy[a] + vz[a]*vz[a]);
            }
            pair_counts[i] = n * (n - 1) / 2;
        }
    }
    __syncthreads();

    if (pair_counts[i] <= 0) return;

    long long total = (long long)pair_counts[i];
    long long base = pair_offsets[i];

    for (long long p = t; p < total; p += blockDim.x) {
        // triangular index p -> lexicographic pair (a,b), 0<=a<b<n
        double dn = (double)n;
        double disc = (2.0*dn - 1.0)*(2.0*dn - 1.0) - 8.0*(double)p;
        int a = (int)floor((2.0*dn - 1.0 - sqrt(disc)) * 0.5);
        long long before = ((long long)a * (2LL*n - a - 1LL)) / 2LL;
        int b = (int)(p - before + a + 1);

        double na = vn[a], nb = vn[b];
        double ang = 0.0;
        if (na > 0.0 && nb > 0.0) {
            double dot = vx[a]*vx[b] + vy[a]*vy[b] + vz[a]*vz[b];
            double c = dot / (na*nb);
            c = c > 1.0 ? 1.0 : (c < -1.0 ? -1.0 : c);
            ang = acos(c) * 57.29577951308232;
        }
        out_angles[base + p] = ang;
    }
}
'''


def _get_cupy():
    global _CUPY
    _configure_cupy_runtime_env()
    if _CUPY is None:
        try:
            import cupy as cp
        except ImportError as exc:
            raise RuntimeError("GPU triplet mode requires CuPy. Install cupy for your CUDA version.") from exc
        _CUPY = cp
    return _CUPY


def _get_fused_kernels(cp):
    global _FUSED_KERNELS
    if _FUSED_KERNELS is None:
        mod = cp.RawModule(
            code=_FUSED_CUDA_SRC,
            options=("-std=c++11",),
            name_expressions=[
                "count_neighbors",
                "gather_neighbors_seq",
                "angle_pairs_from_neighbors",
                "direct_angles_sorted",
            ],
        )
        _FUSED_KERNELS = {
            "count_neighbors": mod.get_function("count_neighbors"),
            "gather_neighbors_seq": mod.get_function("gather_neighbors_seq"),
            "angle_pairs_from_neighbors": mod.get_function("angle_pairs_from_neighbors"),
            "direct_angles_sorted": mod.get_function("direct_angles_sorted"),
        }
    return _FUSED_KERNELS


def _gpu_triplets_fused_cuda(cp, sub, pos, box, low2, high2, profile):
    profile_times = {"h2d": 0.0, "compute": 0.0, "d2h": 0.0}
    kernels = _get_fused_kernels(cp)
    k_count = kernels["count_neighbors"]
    k_gather = kernels["gather_neighbors_seq"]
    k_angle = kernels["angle_pairs_from_neighbors"]

    ns = int(sub.shape[0])
    npw = int(pos.shape[0])
    neigh_counts = cp.zeros((ns,), dtype=cp.int32)

    t_compute = time.perf_counter() if profile else None

    # Pass 1: count neighbors per center.
    k_count(
        (ns,),
        (256,),
        (sub, pos, box, np.int32(ns), np.int32(npw), np.float64(low2), np.float64(high2), neigh_counts),
    )

    pair_counts = (neigh_counts.astype(cp.int64) * (neigh_counts.astype(cp.int64) - 1)) // 2
    neigh_offsets = cp.cumsum(neigh_counts.astype(cp.int64), dtype=cp.int64) - neigh_counts.astype(cp.int64)
    pair_offsets = cp.cumsum(pair_counts, dtype=cp.int64) - pair_counts

    total_neigh = int(cp.asnumpy(cp.sum(neigh_counts, dtype=cp.int64)))
    total_pairs = int(cp.asnumpy(cp.sum(pair_counts, dtype=cp.int64)))
    out_counts = cp.asnumpy(pair_counts).astype(int)

    if total_pairs == 0:
        if profile:
            cp.cuda.Stream.null.synchronize()
            profile_times["compute"] += time.perf_counter() - t_compute
        return np.array([], dtype=float), out_counts, profile_times

    neigh_vecs = cp.empty((total_neigh, 3), dtype=cp.float64)
    out_angles = cp.empty((total_pairs,), dtype=cp.float64)

    # Pass 2: deterministic gather in water-index order.
    k_gather(
        (ns,),
        (1,),
        (sub, pos, box, np.int32(ns), np.int32(npw), np.float64(low2), np.float64(high2), neigh_offsets, neigh_vecs),
    )

    # Pass 3: angle computation from compact neighbor vectors.
    k_angle(
        (ns,),
        (256,),
        (neigh_vecs, neigh_counts, neigh_offsets, pair_offsets, np.int32(ns), out_angles),
    )

    if profile:
        cp.cuda.Stream.null.synchronize()
        profile_times["compute"] += time.perf_counter() - t_compute
        t_d2h = time.perf_counter()
        out_angles_np = cp.asnumpy(out_angles)
        profile_times["d2h"] += time.perf_counter() - t_d2h
    else:
        cp.cuda.Stream.null.synchronize()
        out_angles_np = cp.asnumpy(out_angles)

    return out_angles_np, out_counts, profile_times


def _gpu_triplets_fused_cuda_v2(cp, sub, pos, box, low2, high2, profile):
    profile_times = {"h2d": 0.0, "compute": 0.0, "d2h": 0.0}
    kernels = _get_fused_kernels(cp)
    k_count = kernels["count_neighbors"]
    k_direct = kernels["direct_angles_sorted"]

    ns = int(sub.shape[0])
    npw = int(pos.shape[0])
    neigh_counts = cp.zeros((ns,), dtype=cp.int32)

    t_compute = time.perf_counter() if profile else None

    # Pass 1: count neighbors and precompute pair offsets.
    k_count(
        (ns,),
        (256,),
        (sub, pos, box, np.int32(ns), np.int32(npw), np.float64(low2), np.float64(high2), neigh_counts),
    )
    pair_counts = (neigh_counts.astype(cp.int64) * (neigh_counts.astype(cp.int64) - 1)) // 2
    pair_offsets = cp.cumsum(pair_counts, dtype=cp.int64) - pair_counts

    total_pairs = int(cp.asnumpy(cp.sum(pair_counts, dtype=cp.int64)))
    if total_pairs == 0:
        if profile:
            cp.cuda.Stream.null.synchronize()
            profile_times["compute"] += time.perf_counter() - t_compute
        return np.array([], dtype=float), cp.asnumpy(pair_counts).astype(int), profile_times

    # Pass 2: gather+sort+angle directly in shared memory per center.
    out_angles = cp.empty((total_pairs,), dtype=cp.float64)
    pair_counts_out = cp.zeros((ns,), dtype=cp.int32)
    overflow_flags = cp.zeros((ns,), dtype=cp.int32)
    k_direct(
        (ns,),
        (256,),
        (
            sub,
            pos,
            box,
            np.int32(ns),
            np.int32(npw),
            np.float64(low2),
            np.float64(high2),
            pair_offsets,
            pair_counts_out,
            overflow_flags,
            out_angles,
        ),
    )
    cp.cuda.Stream.null.synchronize()

    overflow_any = bool(int(cp.asnumpy(cp.sum(overflow_flags, dtype=cp.int64))))
    if overflow_any:
        raise RuntimeError(
            f"fused_cuda_v2 overflow: center exceeded {_FUSED_MAX_NEIGH} neighbors within triplet cutoff; "
            "falling back to fused_cuda is required for this frame."
        )

    if profile:
        profile_times["compute"] += time.perf_counter() - t_compute
        t_d2h = time.perf_counter()
        out_angles_np = cp.asnumpy(out_angles)
        out_counts_np = cp.asnumpy(pair_counts_out).astype(int)
        profile_times["d2h"] += time.perf_counter() - t_d2h
    else:
        out_angles_np = cp.asnumpy(out_angles)
        out_counts_np = cp.asnumpy(pair_counts_out).astype(int)

    return out_angles_np, out_counts_np, profile_times


def triplet_angles_gpu(
    sub_positions: np.ndarray,
    all_water_positions: np.ndarray,
    box_dims: np.ndarray,
    low_cut: float,
    high_cut: float,
    profile: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """
    GPU implementation of triplet angle computation using CuPy.
    This path is locked to the fast fused CUDA kernel (fused_cuda_v2), with a
    fused_cuda fallback only when fused_cuda_v2 overflows shared-memory neighbor
    capacity for a center.
    Returns (angles_deg, per-center-pair-counts), optionally plus timing details in profile mode.
    """
    if sub_positions.size == 0:
        if profile:
            return np.array([], dtype=float), np.array([], dtype=int), {"h2d": 0.0, "compute": 0.0, "d2h": 0.0, "mode": "empty"}
        return np.array([], dtype=float), np.array([], dtype=int)

    cp = _get_cupy()

    chosen_mode = "fused_cuda_v2"
    if chosen_mode != "fused_cuda_v2":
        raise ValueError("GPU mode is locked to fused_cuda_v2.")

    h2d_time = 0.0
    t_h2d = time.perf_counter() if profile else None
    sub = cp.asarray(sub_positions, dtype=cp.float64)
    pos = cp.asarray(all_water_positions, dtype=cp.float64)
    box = cp.asarray(box_dims, dtype=cp.float64)
    if profile:
        cp.cuda.Stream.null.synchronize()
        h2d_time = time.perf_counter() - t_h2d
    low2 = float(low_cut) * float(low_cut)
    high2 = float(high_cut) * float(high_cut)

    try:
        angles_np, counts_np, profile_times = _gpu_triplets_fused_cuda_v2(cp, sub, pos, box, low2, high2, profile)
    except RuntimeError as exc:
        if "overflow" not in str(exc):
            raise
        # Rare safety fallback for centers with very large neighbor counts.
        angles_np, counts_np, profile_times = _gpu_triplets_fused_cuda(cp, sub, pos, box, low2, high2, profile)
        chosen_mode = "fused_cuda_fallback"

    if profile:
        profile_times["h2d"] += h2d_time
        profile_times["mode"] = chosen_mode
        return angles_np, counts_np, profile_times
    return angles_np, counts_np

# -------- helpers for insertion codes --------
def parse_res_token(tok: str):
    """
    Parse residue tokens like '112', '-2', '112A', '7c' -> (number, insertion code).
    """
    m = re.fullmatch(r'\s*(-?\d+)\s*([A-Za-z]?)\s*', str(tok))
    if not m:
        raise ValueError(f"Invalid residue token: {tok!r} (expected e.g. -2, 112, or 112B)")
    resseq = int(m.group(1))
    icode = m.group(2).upper()
    return resseq, icode

def residue_selection(chain_id: str, res_token: str, heavy_only=True):
    """
    Build an MDAnalysis selection string for a residue, including insertion code if present.
    If chain_id is None/empty, omit it from the selection.
    """
    resseq, icode = parse_res_token(res_token)
    resid_clause = f"{resseq}{icode}"  # e.g., '112B' or '112'
    if chain_id:
        sel = f"chainID {chain_id} and resid {resid_clause}"
    else:
        sel = f"resid {resid_clause}"
    if heavy_only:
        sel += " and not name H*"
    return sel, resseq, icode

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compile water triplet angles around the residue/"
            "custom group in your protein.\n"
            "Example: python triplet.py protein_processed.pdb traj.dcd -res 42"
        )
    )
    parser.add_argument("protein")
    parser.add_argument("trajectory")
    parser.add_argument(
        "-res",
        "--resid",
        type=str,  # allow tokens like 112B
        help="Residue number (optionally with insertion code), e.g., 112 or 112B",
    )
    parser.add_argument("-ch", "--chain")
    parser.add_argument("-t", "--time", type=float, default=5.0, help="Last x ns to analyse (default 5)")
    parser.add_argument("--groupsFile")
    parser.add_argument("--groupNum", type=int)
    parser.add_argument("--selection")
    parser.add_argument("--hydrationCutoff", type=float, default=4.25)
    parser.add_argument("--hydrationLowCutoff", type=float, default=None)
    parser.add_argument("--outdir", default="angles", help="Output directory (default: angles)")
    parser.add_argument("--skip", type=int, default=1, help="Process every Nth frame (default 1 = every frame)")
    parser.add_argument(
        "--gpu",
        action="store_true",
        help="Use CuPy GPU kernel for triplet-angle math (selection remains MDAnalysis-based).",
    )
    parser.add_argument("--profile", action="store_true", help="Print runtime breakdown for selections and triplet kernel.")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    # some checks
    protein_processed = args.protein if args.protein.endswith(".pdb") else args.protein + ".pdb"
    if not os.path.exists(protein_processed):
        sys.exit(f"Error: cannot find {protein_processed}")
    if not os.path.exists(args.trajectory):
        sys.exit(f"Error: cannot find {args.trajectory}")

    protein_name = os.path.splitext(os.path.basename(protein_processed))[0].replace("_processed", "")

    if args.selection:
        group_string = args.selection.replace(" ", "_")
        out_base = f"{protein_name}_{group_string}"
    elif args.groupsFile and args.groupNum:
        out_base = f"{protein_name}_group{args.groupNum}"
    elif args.resid is not None and args.chain:
        # include insertion code in filename if present (args.resid is a string token, e.g. '112B')
        out_base = f"{protein_name}_res{args.resid}_chain{args.chain}"
    elif args.resid is not None:
        out_base = f"{protein_name}_res{args.resid}"
    else:
        out_base = f"{protein_name}_res{args.resid}"

    if args.hydrationLowCutoff is not None:
        out_base += f"_lowC_{args.hydrationLowCutoff}_highC_{args.hydrationCutoff}"

    angles_dir = args.outdir
    os.makedirs(angles_dir, exist_ok=True)
    output_path = os.path.join(angles_dir, out_base + "_angles.txt")

    # reading trajectory
    u = mda.Universe(protein_processed, args.trajectory)
    total_frames = len(u.trajectory)
    frame_dt_ps = u.trajectory.dt
    frames_to_load = int((args.time * 1e3) / frame_dt_ps)
    first_idx = max(0, total_frames - frames_to_load)

    # MDAnalysis atom selection
    if args.selection:
        sel = args.selection
    elif args.groupsFile and args.groupNum is not None:
        with open(args.groupsFile) as fh:
            sel = fh.readlines()[args.groupNum - 1].strip()
    elif args.resid is not None:
        # Build selection supporting insertion codes and canonical PDB chainID.
        sel, _, _ = residue_selection(args.chain, args.resid, heavy_only=True)
        if len(u.select_atoms(sel)) == 0:
            sys.exit(f"Your atom selection picked zero atoms – tried '{sel}'.")
    else:
        sys.exit("Error: must specify -res/-ch, --groupsFile/--groupNum, or --selection")

    # test selection
    if len(u.select_atoms(sel)) == 0:
        sys.exit("Your atom selection picked zero atoms – please check.")

    waters_sel = NORMALIZED_WATER_OXYGEN_SELECTION
    all_waters = u.select_atoms(waters_sel)
    lowCut = 0.0
    highCut = 3.5

    if args.hydrationLowCutoff:
        shell_sel_expr = (
            f"({waters_sel}) and around {args.hydrationCutoff} ({sel}) "
            f"and not around {args.hydrationLowCutoff} ({sel})"
        )
    else:
        shell_sel_expr = f"({waters_sel}) and around {args.hydrationCutoff} ({sel})"

    # how many frames back from the end we care about
    frames_to_load = int((args.time * 1e3) / u.trajectory.dt)
    first_idx = max(0, total_frames - frames_to_load)

    # --- checkpointing on number of lines already in the output file ---
    done_lines = 0
    if os.path.exists(output_path):
        with open(output_path) as fh:
            for done_lines, _ in enumerate(fh, 1):
                pass
        print(f"Resuming – {done_lines} frames already present in {output_path}")

    # compute the list of frame-indices we actually want to step through
    start_frame = first_idx + done_lines * args.skip
    end_frame = total_frames
    stride = args.skip
    frame_indices = range(start_frame, end_frame, stride)

    profile_totals = {
        "sel_shell": 0.0,
        "all_water_pos": 0.0,
        "cpu_kernel": 0.0,
        "gpu_kernel_total": 0.0,
        "gpu_h2d": 0.0,
        "gpu_compute": 0.0,
        "gpu_d2h": 0.0,
    }
    gpu_mode_counts = {}

    with open(output_path, "a") as out_f:
        # total iterations = number of indices left
        total_iters = len(frame_indices)
        for write_idx, frame_idx in tqdm(
            enumerate(frame_indices, start=done_lines),
            total=total_iters,
            unit="frame",
        ):
            # explicitly jump to that frame
            u.trajectory[frame_idx]

            box_dims = u.dimensions[:3].copy()

            if args.profile:
                t_sel = time.perf_counter()
            shell = u.select_atoms(shell_sel_expr)
            if args.profile:
                profile_totals["sel_shell"] += time.perf_counter() - t_sel

            if args.profile:
                t_all = time.perf_counter()
            all_water_pos = all_waters.positions
            if args.profile:
                profile_totals["all_water_pos"] += time.perf_counter() - t_all

            # compute triplet angles (CPU C-extension by default; optional GPU kernel)
            if args.gpu:
                if args.profile:
                    t_gpu = time.perf_counter()
                    ang_np, _, gpu_times = triplet_angles_gpu(
                        shell.positions,
                        all_water_pos,
                        box_dims,
                        lowCut,
                        highCut,
                        profile=True,
                    )
                    profile_totals["gpu_kernel_total"] += time.perf_counter() - t_gpu
                    profile_totals["gpu_h2d"] += gpu_times["h2d"]
                    profile_totals["gpu_compute"] += gpu_times["compute"]
                    profile_totals["gpu_d2h"] += gpu_times["d2h"]
                    mode_key = gpu_times.get("mode", "fused_cuda_v2")
                    gpu_mode_counts[mode_key] = gpu_mode_counts.get(mode_key, 0) + 1
                else:
                    ang_np, _ = triplet_angles_gpu(
                        shell.positions,
                        all_water_pos,
                        box_dims,
                        lowCut,
                        highCut,
                    )
            else:
                if args.profile:
                    t_cpu = time.perf_counter()
                    ang_np, _ = wl.triplet_angles(
                        shell.positions,
                        all_water_pos,
                        box_dims,
                        lowCut,
                        highCut,
                    )
                    profile_totals["cpu_kernel"] += time.perf_counter() - t_cpu
                else:
                    ang_np, _ = wl.triplet_angles(
                        shell.positions,
                        all_water_pos,
                        box_dims,
                        lowCut,
                        highCut,
                    )

            # write one line per processed frame
            out_f.write(" ".join(f"{a:.1f}" for a in ang_np) + "\n")

            # flush periodically
            if write_idx % 500 == 0:
                out_f.flush()

    print(f"Done! Angles in {output_path}")
    if args.profile and total_iters > 0:
        print("\nTiming summary (--profile)")
        print(f"Frames processed: {total_iters}")
        print(
            f"Shell selection total: {profile_totals['sel_shell']:.6f}s "
            f"({profile_totals['sel_shell'] / total_iters:.6e} s/frame)"
        )
        print(
            f"All-water positions total: {profile_totals['all_water_pos']:.6f}s "
            f"({profile_totals['all_water_pos'] / total_iters:.6e} s/frame)"
        )
        if args.gpu:
            print(
                f"GPU triplet total: {profile_totals['gpu_kernel_total']:.6f}s "
                f"({profile_totals['gpu_kernel_total'] / total_iters:.6e} s/frame)"
            )
            if gpu_mode_counts:
                modes = ", ".join(f"{k}:{v}" for k, v in sorted(gpu_mode_counts.items()))
                print(f"  GPU mode usage: {modes}")
            print(f"  GPU H2D: {profile_totals['gpu_h2d']:.6f}s")
            print(f"  GPU compute: {profile_totals['gpu_compute']:.6f}s")
            print(f"  GPU D2H: {profile_totals['gpu_d2h']:.6f}s")
        else:
            print(
                f"CPU triplet total: {profile_totals['cpu_kernel']:.6f}s "
                f"({profile_totals['cpu_kernel'] / total_iters:.6e} s/frame)"
            )


if __name__ == "__main__":
    main()
