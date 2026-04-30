from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np

from .config import HydroMapConfig
from .manifest import load_json, write_json
from .selection import canonical_group_key
from .types import ParityCaseMetrics
from .workflow import STAGES, WorkflowRunner


@dataclass(frozen=True)
class ParityFixture:
    name: str
    input_dir: Path
    protein: str
    seed: int

    @property
    def case_key(self) -> str:
        return f"{self.name}:{self.protein}:seed{self.seed}"


def fixed_smoke_fixtures(repo_root: Path) -> list[ParityFixture]:
    fixture_root = repo_root / "tests" / "fixtures" / "benchmark_smoke"
    return [
        ParityFixture(name="single_chain", input_dir=fixture_root, protein="alpha", seed=7),
        ParityFixture(name="multi_chain", input_dir=fixture_root, protein="SLS1_dimer", seed=7),
        ParityFixture(
            name="insertion_code",
            input_dir=fixture_root,
            protein="F_insertion_chainP_start10",
            seed=9,
        ),
        ParityFixture(
            name="dirty_stripped",
            input_dir=fixture_root,
            protein="F_dirty_with_hetatm",
            seed=11,
        ),
    ]


def _rows_from_results_csv(results_csv: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with results_csv.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            raw_key = row.get("normalized_group_key")
            selector = row.get("MDAnalysis_selection_strings", "")
            key = raw_key or canonical_group_key(selector)
            pred = float(row["Fdewet_pred"])
            rows.append(
                {
                    "normalized_group_key": key,
                    "MDAnalysis_selection_strings": selector,
                    "Fdewet_pred": pred,
                }
            )
    rows.sort(key=lambda r: str(r["normalized_group_key"]))
    return rows


def _fixture_config(cfg: HydroMapConfig, fixture: ParityFixture) -> HydroMapConfig:
    c = cfg.clone()
    c.input_dir = fixture.input_dir
    c.proteins = [fixture.protein]
    c.seeds = [fixture.seed]
    c.groups_file = None
    c.md.strip_non_protein = True
    return c


def run_parity_baseline(cfg: HydroMapConfig, label: str) -> Path:
    fixtures = fixed_smoke_fixtures(cfg.repo_root)
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")

    cases: dict[str, dict[str, object]] = {}
    for fx in fixtures:
        case_cfg = _fixture_config(cfg, fx)
        run_id = f"parity_{label}_baseline_{fx.name}_{timestamp}"
        runner = WorkflowRunner(case_cfg, run_id=run_id)
        manifest_path = runner.run()[0]
        manifest = load_json(manifest_path)
        results_csv = Path(manifest["files"]["results_csv"])
        cases[fx.case_key] = {
            "fixture": {
                "name": fx.name,
                "protein": fx.protein,
                "seed": fx.seed,
                "input_dir": str(fx.input_dir),
            },
            "manifest": str(manifest_path),
            "rows": _rows_from_results_csv(results_csv),
        }

    payload = {
        "version": 1,
        "label": label,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "threshold": cfg.parity.pooled_pearson_threshold,
        "cases": cases,
    }

    out_path = cfg.artifacts_root / "baselines" / "parity" / f"{label}.json"
    write_json(out_path, payload)
    return out_path


def _pearson(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 2 or len(ys) < 2:
        return None
    x = np.asarray(xs, dtype=float)
    y = np.asarray(ys, dtype=float)
    if float(np.std(x)) == 0.0 or float(np.std(y)) == 0.0:
        return None
    return float(np.corrcoef(x, y)[0, 1])


def _mae(xs: list[float], ys: list[float]) -> float | None:
    if not xs:
        return None
    x = np.asarray(xs, dtype=float)
    y = np.asarray(ys, dtype=float)
    return float(np.mean(np.abs(x - y)))


def _rmse(xs: list[float], ys: list[float]) -> float | None:
    if not xs:
        return None
    x = np.asarray(xs, dtype=float)
    y = np.asarray(ys, dtype=float)
    return float(np.sqrt(np.mean((x - y) ** 2)))


def _compare_case(case_key: str, baseline_rows: list[dict[str, object]], candidate_rows: list[dict[str, object]]) -> ParityCaseMetrics:
    bmap = {str(r["normalized_group_key"]): float(r["Fdewet_pred"]) for r in baseline_rows}
    cmap = {str(r["normalized_group_key"]): float(r["Fdewet_pred"]) for r in candidate_rows}

    bkeys = set(bmap)
    ckeys = set(cmap)
    common = sorted(bkeys & ckeys)

    bx = [bmap[k] for k in common]
    cx = [cmap[k] for k in common]

    return ParityCaseMetrics(
        case_key=case_key,
        matched_rows=len(common),
        baseline_rows=len(bkeys),
        candidate_rows=len(ckeys),
        pearson=_pearson(bx, cx),
        mae=_mae(bx, cx),
        rmse=_rmse(bx, cx),
        missing_keys=sorted(list(bkeys - ckeys)),
        extra_keys=sorted(list(ckeys - bkeys)),
    )


def run_parity_check(cfg: HydroMapConfig, baseline_path: Path | str) -> tuple[Path, bool]:
    baseline_path = Path(baseline_path)
    baseline = load_json(baseline_path)

    fixtures = fixed_smoke_fixtures(cfg.repo_root)
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")

    case_metrics: list[ParityCaseMetrics] = []
    pooled_baseline: list[float] = []
    pooled_candidate: list[float] = []

    for fx in fixtures:
        case_key = fx.case_key
        baseline_case = baseline["cases"].get(case_key)
        if baseline_case is None:
            raise ValueError(f"Baseline missing case: {case_key}")

        case_cfg = _fixture_config(cfg, fx)
        run_id = f"parity_check_{fx.name}_{timestamp}"
        runner = WorkflowRunner(case_cfg, run_id=run_id)
        manifest_path = runner.run()[0]
        manifest = load_json(manifest_path)
        results_csv = Path(manifest["files"]["results_csv"])
        candidate_rows = _rows_from_results_csv(results_csv)

        metrics = _compare_case(case_key, baseline_case["rows"], candidate_rows)
        case_metrics.append(metrics)

        bmap = {str(r["normalized_group_key"]): float(r["Fdewet_pred"]) for r in baseline_case["rows"]}
        cmap = {str(r["normalized_group_key"]): float(r["Fdewet_pred"]) for r in candidate_rows}
        for key in sorted(set(bmap) & set(cmap)):
            pooled_baseline.append(bmap[key])
            pooled_candidate.append(cmap[key])

    pooled_pearson = _pearson(pooled_baseline, pooled_candidate)
    threshold = float(cfg.parity.pooled_pearson_threshold)
    success = pooled_pearson is not None and pooled_pearson >= threshold

    report = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "baseline_path": str(baseline_path),
        "threshold": threshold,
        "pooled_pearson": pooled_pearson,
        "pooled_matched_rows": len(pooled_baseline),
        "success": success,
        "cases": [m.__dict__ for m in case_metrics],
        "diagnostics": {
            "pooled_mae": _mae(pooled_baseline, pooled_candidate),
            "pooled_rmse": _rmse(pooled_baseline, pooled_candidate),
        },
    }

    report_path = cfg.artifacts_root / "reports" / f"parity_check_{timestamp}.json"
    write_json(report_path, report)
    return report_path, success


def run_benchmark_smoke(
    cfg: HydroMapConfig,
    stages: Iterable[str] | None = None,
) -> tuple[Path, bool]:
    fixtures = fixed_smoke_fixtures(cfg.repo_root)
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    stage_list = list(stages) if stages is not None else ["prepare"]
    invalid = [stage for stage in stage_list if stage not in STAGES]
    if invalid:
        raise ValueError(f"Unknown benchmark stages: {invalid}. Supported: {STAGES}")

    cases: list[dict[str, object]] = []
    success = True
    for fx in fixtures:
        case_cfg = _fixture_config(cfg, fx)
        run_id = f"benchmark_smoke_{fx.name}_{timestamp}"
        runner = WorkflowRunner(case_cfg, run_id=run_id)
        try:
            manifest_path = runner.run(stages=stage_list)[0]
            manifest = load_json(manifest_path)
            prepare_report_path = Path(manifest["files"]["prepare_report"])
            prepare_report = load_json(prepare_report_path) if prepare_report_path.exists() else None
            cases.append(
                {
                    "case_key": fx.case_key,
                    "fixture": {
                        "name": fx.name,
                        "protein": fx.protein,
                        "seed": fx.seed,
                        "input_dir": str(fx.input_dir),
                    },
                    "success": True,
                    "manifest": str(manifest_path),
                    "stages": stage_list,
                    "results_summary": manifest.get("results_summary"),
                    "prepare_status": prepare_report.get("status") if prepare_report else None,
                    "prepare_warnings": prepare_report.get("warnings") if prepare_report else [],
                }
            )
        except Exception as exc:
            success = False
            cases.append(
                {
                    "case_key": fx.case_key,
                    "fixture": {
                        "name": fx.name,
                        "protein": fx.protein,
                        "seed": fx.seed,
                        "input_dir": str(fx.input_dir),
                    },
                    "success": False,
                    "stages": stage_list,
                    "error": str(exc),
                }
            )

    report = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "suite": "smoke",
        "stages": stage_list,
        "success": success,
        "fixture_count": len(fixtures),
        "passed_count": sum(1 for case in cases if case["success"]),
        "failed_count": sum(1 for case in cases if not case["success"]),
        "cases": cases,
    }
    report_path = cfg.artifacts_root / "reports" / f"benchmark_smoke_{timestamp}.json"
    write_json(report_path, report)
    return report_path, success
