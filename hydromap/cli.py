from __future__ import annotations

import argparse
from pathlib import Path
import sys

from .config import load_config, package_repo_root
from .parity import run_benchmark_smoke, run_parity_baseline, run_parity_check
from .workflow import STAGES, WorkflowRunner


def _add_common_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--config", required=True, help="Path to HydroMap YAML config.")
    parser.add_argument("--run-id", default=None, help="Optional run identifier; defaults to UTC timestamp.")
    parser.add_argument("--proteins", nargs="+", default=None, help="Optional override list of proteins.")
    parser.add_argument("--seeds", nargs="+", type=int, default=None, help="Optional override list of seeds.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="hydromap", description="HydroMap v2 workflow CLI")
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="Run prepare+simulate+analyze+predict+color")
    _add_common_options(run_parser)
    run_parser.add_argument("--stages", nargs="+", choices=STAGES, default=None, help="Optional subset of stages.")

    for stage in STAGES:
        stage_parser = subparsers.add_parser(stage, help=f"Run only the {stage} stage")
        _add_common_options(stage_parser)

    parity_parser = subparsers.add_parser("parity", help="Prediction parity utilities")
    parity_sub = parity_parser.add_subparsers(dest="parity_command", required=True)

    baseline = parity_sub.add_parser("baseline", help="Generate a parity baseline from the smoke fixtures")
    baseline.add_argument("--config", required=True, help="Path to HydroMap YAML config.")
    baseline.add_argument("--label", required=True, help="Baseline label stored under artifacts/baselines/parity.")

    check = parity_sub.add_parser("check", help="Compare current predictions against a saved baseline")
    check.add_argument("--config", required=True, help="Path to HydroMap YAML config.")
    check.add_argument("--baseline", required=True, help="Baseline JSON path from parity baseline command.")

    benchmark_parser = subparsers.add_parser("benchmark", help="Run fixed-fixture validation suites")
    benchmark_sub = benchmark_parser.add_subparsers(dest="benchmark_command", required=True)

    smoke = benchmark_sub.add_parser("smoke", help="Run the smoke fixture suite and summarize outcomes")
    smoke.add_argument("--config", required=True, help="Path to HydroMap YAML config.")
    smoke.add_argument(
        "--stages",
        nargs="+",
        choices=STAGES,
        default=["prepare"],
        help="Workflow stages to run for each smoke fixture (default: prepare).",
    )

    return parser


def _run_workflow(command: str, args: argparse.Namespace) -> int:
    cfg = load_config(args.config, repo_root=package_repo_root())
    runner = WorkflowRunner(cfg, run_id=args.run_id)

    stages = args.stages if command == "run" else [command]
    manifests = runner.run(stages=stages, proteins=args.proteins, seeds=args.seeds)
    print(f"Run complete. Wrote {len(manifests)} manifest(s):")
    for path in manifests:
        print(path)
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        if args.command in STAGES or args.command == "run":
            return _run_workflow(args.command, args)

        if args.command == "parity" and args.parity_command == "baseline":
            cfg = load_config(args.config, repo_root=package_repo_root())
            path = run_parity_baseline(cfg, label=args.label)
            print(f"Wrote parity baseline: {path}")
            return 0

        if args.command == "parity" and args.parity_command == "check":
            cfg = load_config(args.config, repo_root=package_repo_root())
            report_path, success = run_parity_check(cfg, baseline_path=args.baseline)
            print(f"Wrote parity report: {report_path}")
            return 0 if success else 2

        if args.command == "benchmark" and args.benchmark_command == "smoke":
            cfg = load_config(args.config, repo_root=package_repo_root())
            report_path, success = run_benchmark_smoke(cfg, stages=args.stages)
            print(f"Wrote benchmark report: {report_path}")
            return 0 if success else 2

        parser.error("Unsupported command")
        return 2
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
