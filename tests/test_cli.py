from hydromap.cli import build_parser


def test_cli_has_expected_commands() -> None:
    parser = build_parser()
    args = parser.parse_args(["run", "--config", "configs/example_gpu.yaml"])
    assert args.command == "run"

    args = parser.parse_args(["parity", "baseline", "--config", "configs/tiny_repro.yaml", "--label", "v2"])
    assert args.command == "parity"
    assert args.parity_command == "baseline"

    args = parser.parse_args(["parity", "check", "--config", "configs/tiny_repro.yaml", "--baseline", "x.json"])
    assert args.command == "parity"
    assert args.parity_command == "check"

    args = parser.parse_args(["benchmark", "smoke", "--config", "configs/tiny_repro.yaml"])
    assert args.command == "benchmark"
    assert args.benchmark_command == "smoke"
    assert args.stages == ["prepare"]
