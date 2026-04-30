from pathlib import Path

from hydromap.parity import _compare_case, fixed_smoke_fixtures, run_benchmark_smoke


def test_compare_case_metrics() -> None:
    baseline = [
        {"normalized_group_key": "k1", "Fdewet_pred": 1.0},
        {"normalized_group_key": "k2", "Fdewet_pred": 2.0},
        {"normalized_group_key": "k3", "Fdewet_pred": 3.0},
    ]
    candidate = [
        {"normalized_group_key": "k1", "Fdewet_pred": 1.1},
        {"normalized_group_key": "k2", "Fdewet_pred": 1.9},
        {"normalized_group_key": "k4", "Fdewet_pred": 4.0},
    ]

    m = _compare_case("case", baseline, candidate)
    assert m.matched_rows == 2
    assert m.baseline_rows == 3
    assert m.candidate_rows == 3
    assert "k3" in m.missing_keys
    assert "k4" in m.extra_keys


def test_fixed_smoke_fixtures_point_to_checked_in_pdbs() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    fixtures = fixed_smoke_fixtures(repo_root)
    assert fixtures
    for fixture in fixtures:
        assert (fixture.input_dir / f"{fixture.protein}.pdb").exists()


def test_run_benchmark_smoke_writes_report(tmp_path: Path, monkeypatch) -> None:
    class DummyCfg:
        def __init__(self) -> None:
            self.repo_root = tmp_path
            self.artifacts_root = tmp_path / "artifacts"
            self.md = type("MD", (), {"strip_non_protein": False})()
            self.input_dir = tmp_path
            self.proteins = []
            self.seeds = []
            self.groups_file = None

        def clone(self):
            clone = DummyCfg()
            clone.repo_root = self.repo_root
            clone.artifacts_root = self.artifacts_root
            clone.md.strip_non_protein = self.md.strip_non_protein
            return clone

    fixture = type(
        "Fixture",
        (),
        {
            "name": "tiny",
            "protein": "alpha",
            "seed": 7,
            "input_dir": tmp_path,
            "case_key": "tiny:alpha:seed7",
        },
    )()

    class DummyRunner:
        def __init__(self, cfg, run_id=None):
            self.cfg = cfg
            self.run_id = run_id

        def run(self, stages=None):
            case_root = tmp_path / "artifacts" / "runs" / "dummy" / "alpha" / "seed_7"
            case_root.mkdir(parents=True, exist_ok=True)
            prepare_report = case_root / "prepare_report.json"
            prepare_report.write_text('{"status":"ok","warnings":[]}\n', encoding="utf-8")
            manifest = case_root / "manifest.json"
            manifest.write_text(
                (
                    "{"
                    '"files":{"prepare_report":"%s"},'
                    '"results_summary":{"rows":0}'
                    "}\n"
                )
                % str(prepare_report),
                encoding="utf-8",
            )
            return [manifest]

    monkeypatch.setattr("hydromap.parity.fixed_smoke_fixtures", lambda repo_root: [fixture])
    monkeypatch.setattr("hydromap.parity.WorkflowRunner", DummyRunner)

    report_path, success = run_benchmark_smoke(DummyCfg(), stages=["prepare"])
    assert success is True
    payload = report_path.read_text(encoding="utf-8")
    assert "benchmark_smoke" in report_path.name
    assert '"success": true' in payload
    assert '"prepare_status": "ok"' in payload
