from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
import copy
import hashlib
import json
from typing import Any

from .errors import ConfigError


@dataclass
class MDConfig:
    device: str = "gpu"  # gpu|cpu|auto
    nanoseconds: float = 0.5
    restrain: str | None = None
    restrain_selection: str | None = None
    restraint_k: float = 1000.0
    equilibration_ns: float = 0.1
    timestep_ps: float = 0.003
    report_interval_ps: float = 1.0
    strip_non_protein: bool = True
    random_seed: int | None = None
    preprocess_seed: int | None = None
    velocity_seed: int | None = None
    barostat_seed: int | None = None
    checkpoint_policy: str = "error"  # error|resume|overwrite
    no_cuda: bool = False
    allow_cpu_md: bool = False
    deterministic: bool = True
    cuda_precision: str = "mixed"  # single|mixed|double
    repair_missing_atoms: str = "none"  # none|pdbfixer
    capping_mode: str = "none"  # none|termini|breaks|termini_and_breaks
    prep_policy: str = "permissive"  # permissive|strict
    neutralize: bool = True
    ionic_strength_molar: float = 0.0
    positive_ion: str = "Na+"
    negative_ion: str = "Cl-"
    histidine_mode: str = "auto"  # auto|hid|hie|hip
    histidine_overrides: dict[str, str] = field(default_factory=dict)


@dataclass
class AnalysisConfig:
    compute_potentials: bool = True
    triplet_histogram_bin_width_deg: float = 10.0
    triplets_frame_stride: int = 1
    triplets_sample_ps: float | None = None
    triplets_hydration_cutoff: float = 4.25
    triplets_device: str = "gpu"  # gpu|cpu|auto
    potentials_frame_stride: int = 100
    potentials_sample_ps: float | None = None
    potentials_cutoff: float = 4.25
    discard_initial_ns: float = 0.0
    tail_ns: float | None = None
    triplets_tail_ns: float | None = None
    potentials_tail_ns: float | None = None
    potentials_device: str = "gpu"  # gpu|cpu|auto
    existing_processed_pdb: str | None = None
    existing_trajectory: str | None = None
    existing_topology: str | None = None
    min_waters: float = 5.0
    color_properties: list[str] = field(default_factory=list)


@dataclass
class ResourceConfig:
    max_cpu_workers: int = 8
    reserve_cpus: int = 2
    max_gpu_jobs: int = 1


@dataclass
class ExecutionConfig:
    profile: str = "balanced"  # gpu_fast|balanced


@dataclass
class ParityConfig:
    pooled_pearson_threshold: float = 0.96


@dataclass
class HydroMapConfig:
    config_path: Path
    repo_root: Path
    input_dir: Path
    artifacts_root: Path
    run_id: str | None
    proteins: list[str]
    seeds: list[int]
    model_path: Path
    forcefield: str
    groups_file: Path | None
    md: MDConfig = field(default_factory=MDConfig)
    analysis: AnalysisConfig = field(default_factory=AnalysisConfig)
    resources: ResourceConfig = field(default_factory=ResourceConfig)
    execution: ExecutionConfig = field(default_factory=ExecutionConfig)
    parity: ParityConfig = field(default_factory=ParityConfig)

    def clone(self) -> "HydroMapConfig":
        return copy.deepcopy(self)

    def config_digest(self) -> str:
        payload = {
            "input_dir": str(self.input_dir),
            "artifacts_root": str(self.artifacts_root),
            "proteins": self.proteins,
            "seeds": self.seeds,
            "model_path": str(self.model_path),
            "forcefield": self.forcefield,
            "groups_file": str(self.groups_file) if self.groups_file else None,
            "md": asdict(self.md),
            "analysis": asdict(self.analysis),
            "resources": asdict(self.resources),
            "execution": asdict(self.execution),
            "parity": asdict(self.parity),
        }
        raw = json.dumps(payload, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def _as_list(raw: Any) -> list[Any]:
    if raw is None:
        return []
    if isinstance(raw, list):
        return raw
    return [raw]


def _resolve_path(base: Path, maybe_path: str | None) -> Path | None:
    if maybe_path is None:
        return None
    p = Path(maybe_path)
    if not p.is_absolute():
        p = (base / p).resolve()
    return p


def _optional_float(value: Any) -> float | None:
    if value is None:
        return None
    return float(value)


def _optional_path_template(value: Any) -> str | None:
    if value is None:
        return None
    return str(value).strip() or None


def _normalize_analysis_device(raw_device: Any, raw_backend: Any, default_device: str) -> str:
    if raw_device is not None:
        return str(raw_device).strip().lower()

    if raw_backend is None:
        return default_device

    backend = str(raw_backend).strip().lower()
    backend_map = {
        "gpu_serial": "gpu",
        "cpu_serial": "cpu",
        "cpu_parallel": "cpu",
    }
    return backend_map.get(backend, backend)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ConfigError(message)


def package_repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def load_config(config_path: str | Path, repo_root: str | Path | None = None) -> HydroMapConfig:
    try:
        import yaml
    except ImportError as exc:
        raise ImportError("PyYAML is required. Install with `pip install pyyaml`.") from exc

    config_path = Path(config_path).resolve()
    _require(config_path.exists(), f"Config file not found: {config_path}")

    cfg_dir = config_path.parent
    with config_path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}

    if "multi_chain" in data:
        raise ConfigError(
            "`multi_chain` mode was removed in v2. Residue mode is chain-aware by default; "
            "use `groups_file` for custom grouping."
        )

    repo_root_path = Path(repo_root).resolve() if repo_root else package_repo_root()

    proteins = [str(p) for p in _as_list(data.get("proteins", data.get("protein")))]
    seeds = [int(s) for s in _as_list(data.get("seeds", data.get("seed", [42])))]

    input_dir = _resolve_path(cfg_dir, data.get("input_dir", "."))
    if data.get("artifacts_root") is None:
        artifacts_root = (repo_root_path / "artifacts").resolve()
    else:
        artifacts_root = _resolve_path(cfg_dir, data.get("artifacts_root"))
    if data.get("model_path") is None:
        model_path = (repo_root_path / "hydromap" / "models" / "Fdewet.joblib").resolve()
    else:
        model_path = _resolve_path(cfg_dir, data.get("model_path"))
    groups_file = _resolve_path(cfg_dir, data.get("groups_file"))

    md_raw = data.get("md", {}) or {}
    analysis_raw = data.get("analysis", {}) or {}
    resources_raw = data.get("resources", {}) or {}
    execution_raw = data.get("execution", {}) or {}
    parity_raw = data.get("parity", {}) or {}

    prepare_backend_raw = md_raw.get("prepare_backend")
    if prepare_backend_raw is not None:
        prepare_backend = str(prepare_backend_raw).strip().lower()
        if prepare_backend != "openmm":
            raise ConfigError(
                "md.prepare_backend was removed. HydroMap now uses the OpenMM preparation path only."
            )

    restrain = md_raw.get("restrain", MDConfig.restrain)
    restrain_selection = md_raw.get("restrain_selection", MDConfig.restrain_selection)
    if restrain is not None and restrain_selection is not None:
        raise ConfigError("Use only one of md.restrain or md.restrain_selection, not both.")
    resolved_restrain = restrain_selection if restrain_selection is not None else restrain

    legacy_no_cuda = bool(md_raw.get("no_cuda", MDConfig.no_cuda))
    md_device = str(md_raw.get("device", "cpu" if legacy_no_cuda else MDConfig.device)).strip().lower()

    md = MDConfig(
        device=md_device,
        nanoseconds=float(md_raw.get("nanoseconds", MDConfig.nanoseconds)),
        restrain=resolved_restrain,
        restrain_selection=resolved_restrain,
        restraint_k=float(md_raw.get("restraint_k", MDConfig.restraint_k)),
        equilibration_ns=float(md_raw.get("equilibration_ns", MDConfig.equilibration_ns)),
        timestep_ps=float(md_raw.get("timestep_ps", MDConfig.timestep_ps)),
        report_interval_ps=float(md_raw.get("report_interval_ps", MDConfig.report_interval_ps)),
        strip_non_protein=bool(md_raw.get("strip_non_protein", MDConfig.strip_non_protein)),
        random_seed=(None if md_raw.get("random_seed") is None else int(md_raw.get("random_seed"))),
        preprocess_seed=(None if md_raw.get("preprocess_seed") is None else int(md_raw.get("preprocess_seed"))),
        velocity_seed=(None if md_raw.get("velocity_seed") is None else int(md_raw.get("velocity_seed"))),
        barostat_seed=(None if md_raw.get("barostat_seed") is None else int(md_raw.get("barostat_seed"))),
        checkpoint_policy=str(md_raw.get("checkpoint_policy", MDConfig.checkpoint_policy)),
        no_cuda=(md_device == "cpu"),
        allow_cpu_md=bool(md_raw.get("allow_cpu_md", MDConfig.allow_cpu_md)),
        deterministic=bool(md_raw.get("deterministic", MDConfig.deterministic)),
        cuda_precision=str(md_raw.get("cuda_precision", MDConfig.cuda_precision)),
        repair_missing_atoms=str(md_raw.get("repair_missing_atoms", MDConfig.repair_missing_atoms)).strip().lower(),
        capping_mode=str(md_raw.get("capping_mode", MDConfig.capping_mode)).strip().lower(),
        prep_policy=str(md_raw.get("prep_policy", MDConfig.prep_policy)).strip().lower(),
        neutralize=bool(md_raw.get("neutralize", MDConfig.neutralize)),
        ionic_strength_molar=float(md_raw.get("ionic_strength_molar", MDConfig.ionic_strength_molar)),
        positive_ion=str(md_raw.get("positive_ion", MDConfig.positive_ion)),
        negative_ion=str(md_raw.get("negative_ion", MDConfig.negative_ion)),
        histidine_mode=str(md_raw.get("histidine_mode", MDConfig.histidine_mode)).strip().lower(),
        histidine_overrides={
            str(k).strip(): str(v).strip().upper()
            for k, v in (md_raw.get("histidine_overrides", {}) or {}).items()
        },
    )

    triplets_stride = int(
        analysis_raw.get("triplets_frame_stride", analysis_raw.get("triplets_skip", AnalysisConfig.triplets_frame_stride))
    )
    potentials_stride = int(
        analysis_raw.get(
            "potentials_frame_stride",
            analysis_raw.get("potentials_skip", AnalysisConfig.potentials_frame_stride),
        )
    )

    analysis = AnalysisConfig(
        compute_potentials=bool(
            analysis_raw.get("compute_potentials", AnalysisConfig.compute_potentials)
        ),
        triplet_histogram_bin_width_deg=float(
            analysis_raw.get(
                "triplet_histogram_bin_width_deg",
                AnalysisConfig.triplet_histogram_bin_width_deg,
            )
        ),
        triplets_frame_stride=triplets_stride,
        triplets_sample_ps=_optional_float(analysis_raw.get("triplets_sample_ps")),
        triplets_hydration_cutoff=float(
            analysis_raw.get("triplets_hydration_cutoff", AnalysisConfig.triplets_hydration_cutoff)
        ),
        triplets_device=_normalize_analysis_device(
            analysis_raw.get("triplets_device"),
            analysis_raw.get("triplets_backend"),
            AnalysisConfig.triplets_device,
        ),
        potentials_frame_stride=potentials_stride,
        potentials_sample_ps=_optional_float(analysis_raw.get("potentials_sample_ps")),
        potentials_cutoff=float(analysis_raw.get("potentials_cutoff", AnalysisConfig.potentials_cutoff)),
        discard_initial_ns=float(analysis_raw.get("discard_initial_ns", AnalysisConfig.discard_initial_ns)),
        tail_ns=_optional_float(analysis_raw.get("tail_ns")),
        # Backward-compatible aliases: *_time_ns means "last X ns".
        triplets_tail_ns=_optional_float(
            analysis_raw.get("triplets_tail_ns", analysis_raw.get("triplets_time_ns"))
        ),
        potentials_tail_ns=_optional_float(
            analysis_raw.get("potentials_tail_ns", analysis_raw.get("potentials_time_ns"))
        ),
        potentials_device=_normalize_analysis_device(
            analysis_raw.get("potentials_device"),
            analysis_raw.get("potentials_backend"),
            AnalysisConfig.potentials_device,
        ),
        existing_processed_pdb=_optional_path_template(analysis_raw.get("existing_processed_pdb")),
        existing_trajectory=_optional_path_template(analysis_raw.get("existing_trajectory")),
        existing_topology=_optional_path_template(analysis_raw.get("existing_topology")),
        min_waters=float(analysis_raw.get("min_waters", AnalysisConfig.min_waters)),
        color_properties=[str(x) for x in _as_list(analysis_raw.get("color_properties", []))],
    )

    resources = ResourceConfig(
        max_cpu_workers=int(resources_raw.get("max_cpu_workers", ResourceConfig.max_cpu_workers)),
        reserve_cpus=int(resources_raw.get("reserve_cpus", ResourceConfig.reserve_cpus)),
        max_gpu_jobs=int(resources_raw.get("max_gpu_jobs", ResourceConfig.max_gpu_jobs)),
    )

    execution = ExecutionConfig(profile=str(execution_raw.get("profile", ExecutionConfig.profile)))

    parity = ParityConfig(
        pooled_pearson_threshold=float(
            parity_raw.get("pooled_pearson_threshold", ParityConfig.pooled_pearson_threshold)
        )
    )

    cfg = HydroMapConfig(
        config_path=config_path,
        repo_root=repo_root_path,
        input_dir=input_dir,
        artifacts_root=artifacts_root,
        run_id=data.get("run_id"),
        proteins=proteins,
        seeds=seeds,
        model_path=model_path,
        forcefield=str(data.get("forcefield", "a99SBdisp")),
        groups_file=groups_file,
        md=md,
        analysis=analysis,
        resources=resources,
        execution=execution,
        parity=parity,
    )

    validate_config(cfg)
    return cfg


def validate_config(cfg: HydroMapConfig) -> None:
    _require(bool(cfg.proteins), "Config must define at least one protein via `protein` or `proteins`.")
    _require(bool(cfg.seeds), "Config must define at least one seed via `seed` or `seeds`.")
    _require(cfg.model_path is not None and cfg.model_path.exists(), f"model_path not found: {cfg.model_path}")

    has_external_inputs = (
        cfg.analysis.existing_processed_pdb is not None
        and cfg.analysis.existing_trajectory is not None
    )
    if has_external_inputs:
        pass
    else:
        _require(cfg.input_dir.exists(), f"input_dir does not exist: {cfg.input_dir}")
        for protein in cfg.proteins:
            pdb_path = cfg.input_dir / f"{protein}.pdb"
            _require(pdb_path.exists(), f"Missing input PDB: {pdb_path}")

    if cfg.groups_file is not None:
        _require(cfg.groups_file.exists(), f"groups_file not found: {cfg.groups_file}")

    _require(
        cfg.md.checkpoint_policy in {"error", "resume", "overwrite"},
        "md.checkpoint_policy must be one of: error, resume, overwrite.",
    )
    _require(cfg.md.equilibration_ns >= 0.0, "md.equilibration_ns must be >= 0.")
    _require(cfg.md.device in {"gpu", "cpu", "auto"}, "md.device must be one of: gpu, cpu, auto.")
    _require(cfg.md.timestep_ps > 0.0, "md.timestep_ps must be > 0.")
    _require(cfg.md.report_interval_ps > 0.0, "md.report_interval_ps must be > 0.")
    _require(cfg.md.ionic_strength_molar >= 0.0, "md.ionic_strength_molar must be >= 0.")
    _require(
        cfg.md.report_interval_ps >= cfg.md.timestep_ps,
        "md.report_interval_ps must be >= md.timestep_ps.",
    )
    _require(
        cfg.md.cuda_precision in {"single", "mixed", "double"},
        "md.cuda_precision must be one of: single, mixed, double.",
    )
    _require(
        cfg.md.repair_missing_atoms in {"none", "pdbfixer"},
        "md.repair_missing_atoms must be one of: none, pdbfixer.",
    )
    _require(
        cfg.md.capping_mode in {"none", "termini", "breaks", "termini_and_breaks"},
        "md.capping_mode must be one of: none, termini, breaks, termini_and_breaks.",
    )
    _require(
        cfg.md.prep_policy in {"permissive", "strict"},
        "md.prep_policy must be one of: permissive, strict.",
    )
    _require(
        cfg.md.positive_ion in {"Na+", "K+", "Li+", "Rb+", "Cs+"},
        "md.positive_ion must be one of: Na+, K+, Li+, Rb+, Cs+.",
    )
    _require(
        cfg.md.negative_ion in {"Cl-", "Br-", "F-", "I-"},
        "md.negative_ion must be one of: Cl-, Br-, F-, I-.",
    )
    _require(
        cfg.md.histidine_mode in {"auto", "hid", "hie", "hip"},
        "md.histidine_mode must be one of: auto, hid, hie, hip.",
    )
    _require(
        all(value in {"HID", "HIE", "HIP"} for value in cfg.md.histidine_overrides.values()),
        "md.histidine_overrides values must be one of: HID, HIE, HIP.",
    )
    _require(
        cfg.execution.profile in {"gpu_fast", "balanced"},
        "execution.profile must be one of: gpu_fast, balanced.",
    )

    _require(cfg.resources.max_cpu_workers >= 1, "resources.max_cpu_workers must be >= 1")
    _require(cfg.resources.reserve_cpus >= 0, "resources.reserve_cpus must be >= 0")
    _require(cfg.resources.max_gpu_jobs >= 1, "resources.max_gpu_jobs must be >= 1")

    _require(cfg.analysis.triplets_frame_stride >= 1, "analysis.triplets_frame_stride must be >= 1")
    _require(
        cfg.analysis.triplet_histogram_bin_width_deg > 0.0
        and 140.0 / cfg.analysis.triplet_histogram_bin_width_deg
        == round(140.0 / cfg.analysis.triplet_histogram_bin_width_deg),
        "analysis.triplet_histogram_bin_width_deg must be a positive divisor of 140 degrees.",
    )
    _require(cfg.analysis.potentials_frame_stride >= 1, "analysis.potentials_frame_stride must be >= 1")
    if cfg.analysis.triplets_sample_ps is not None:
        _require(cfg.analysis.triplets_sample_ps > 0.0, "analysis.triplets_sample_ps must be > 0.")
    if cfg.analysis.potentials_sample_ps is not None:
        _require(cfg.analysis.potentials_sample_ps > 0.0, "analysis.potentials_sample_ps must be > 0.")
    _require(
        cfg.analysis.triplets_device in {"gpu", "cpu", "auto"},
        "analysis.triplets_device must be one of: gpu, cpu, auto.",
    )
    _require(cfg.analysis.discard_initial_ns >= 0.0, "analysis.discard_initial_ns must be >= 0.")
    if cfg.analysis.tail_ns is not None:
        _require(cfg.analysis.tail_ns > 0.0, "analysis.tail_ns must be > 0 when set.")
    if cfg.analysis.triplets_tail_ns is not None:
        _require(cfg.analysis.triplets_tail_ns > 0.0, "analysis.triplets_tail_ns must be > 0 when set.")
    if cfg.analysis.potentials_tail_ns is not None:
        _require(cfg.analysis.potentials_tail_ns > 0.0, "analysis.potentials_tail_ns must be > 0 when set.")
    _require(
        cfg.analysis.potentials_device in {"gpu", "cpu", "auto"},
        "analysis.potentials_device must be one of: gpu, cpu, auto.",
    )
    existing_inputs = [
        cfg.analysis.existing_processed_pdb,
        cfg.analysis.existing_trajectory,
    ]
    if any(existing_inputs):
        _require(
            all(existing_inputs),
            "analysis.existing_processed_pdb and analysis.existing_trajectory must either both be set or both be omitted.",
        )
        if cfg.analysis.existing_topology:
            existing_topology = Path(cfg.analysis.existing_topology)
            _require(
                existing_topology.suffix.lower() in {".xml", ".top"},
                "analysis.existing_topology must point to a serialized OpenMM System XML file or a GROMACS .top file.",
            )
        _require(
            len(cfg.proteins) == 1,
            "analysis.existing_* inputs support exactly one protein in the config.",
        )
    _require(
        0.0 <= cfg.parity.pooled_pearson_threshold <= 1.0,
        "parity.pooled_pearson_threshold must be in [0, 1].",
    )


def config_to_dict(cfg: HydroMapConfig) -> dict[str, Any]:
    return {
        "config_path": str(cfg.config_path),
        "repo_root": str(cfg.repo_root),
        "input_dir": str(cfg.input_dir),
        "artifacts_root": str(cfg.artifacts_root),
        "run_id": cfg.run_id,
        "proteins": cfg.proteins,
        "seeds": cfg.seeds,
        "model_path": str(cfg.model_path),
        "forcefield": cfg.forcefield,
        "groups_file": str(cfg.groups_file) if cfg.groups_file else None,
        "md": asdict(cfg.md),
        "analysis": asdict(cfg.analysis),
        "resources": asdict(cfg.resources),
        "execution": asdict(cfg.execution),
        "parity": asdict(cfg.parity),
        "config_digest": cfg.config_digest(),
    }
