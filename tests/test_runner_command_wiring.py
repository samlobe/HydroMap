from pathlib import Path

import pytest

from hydromap.config import load_config
from hydromap.errors import GPUNotAvailableError
from hydromap.runner import CaseSpec, HydroMapRunner


REPO_ROOT = Path(__file__).resolve().parents[1]


def make_cfg(
    tmp_path: Path,
    allow_cpu_md: bool = True,
    analysis_block: str = "",
) -> Path:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture_dir = REPO_ROOT / "tests" / "fixtures"
    cfg = tmp_path / "config.yaml"
    analysis_yaml = (
        analysis_block.strip()
        if analysis_block.strip()
        else "triplets_device: cpu\npotentials_device: cpu"
    )
    analysis_yaml = "\n".join(f"  {line}" for line in analysis_yaml.splitlines())
    cfg.write_text(
        f"""
input_dir: {fixture_dir}
protein: tiny_protein
seed: 123
model_path: {model}
artifacts_root: {tmp_path / 'artifacts'}
md:
  device: gpu
  nanoseconds: 0.005
  random_seed: 77
  velocity_seed: 88
  barostat_seed: 99
  checkpoint_policy: overwrite
  deterministic: true
  allow_cpu_md: {'true' if allow_cpu_md else 'false'}
execution:
  profile: gpu_fast
analysis:
{analysis_yaml}
""",
        encoding="utf-8",
    )
    return cfg


def test_simulate_falls_back_to_cpu_only_if_allowed(tmp_path: Path) -> None:
    cfg = load_config(make_cfg(tmp_path, allow_cpu_md=True), repo_root=REPO_ROOT)
    runner = HydroMapRunner(cfg, run_id="test")
    case = CaseSpec("tiny_protein", 123)
    paths = runner._case_paths(case)

    paths.simulation.mkdir(parents=True, exist_ok=True)
    paths.processed_pdb.parent.mkdir(parents=True, exist_ok=True)
    paths.processed_pdb.write_text("dummy", encoding="utf-8")
    paths.topology.write_text("dummy", encoding="utf-8")

    captured: list[list[str]] = []

    def fake_run_command(case_arg, stage, cmd, cwd):
        captured.append(cmd)
        paths.trajectory.write_text("dummy", encoding="utf-8")

    runner._run_command = fake_run_command  # type: ignore[method-assign]
    runner._cuda_available = lambda: False  # type: ignore[method-assign]
    meta = runner._run_simulate(case, paths)

    cmd = captured[0]
    assert cmd[3] == paths.topology.name
    assert "--noCUDA" in cmd
    assert "--equilibration_ps" in cmd
    assert meta["cpu_md_fallback"] is True


def test_prepare_passes_ion_controls_to_openmm_prep(tmp_path: Path) -> None:
    cfg_path = make_cfg(tmp_path, allow_cpu_md=True)
    cfg = load_config(cfg_path, repo_root=REPO_ROOT)
    cfg.md.neutralize = False
    cfg.md.ionic_strength_molar = 0.15
    cfg.md.positive_ion = "K+"
    cfg.md.negative_ion = "Br-"

    runner = HydroMapRunner(cfg, run_id="test")
    case = CaseSpec("tiny_protein", 123)
    paths = runner._case_paths(case)

    captured: list[list[str]] = []

    def fake_run_command(case_arg, stage, cmd, cwd):
        captured.append(cmd)
        paths.processed_pdb.parent.mkdir(parents=True, exist_ok=True)
        paths.processed_pdb.write_text("dummy", encoding="utf-8")
        paths.topology.parent.mkdir(parents=True, exist_ok=True)
        paths.topology.write_text("dummy", encoding="utf-8")

    def fake_sanitize(src, dst):
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text("dummy", encoding="utf-8")
        return {"ok": True}

    runner._run_command = fake_run_command  # type: ignore[method-assign]
    runner._sanitize_input_pdb = fake_sanitize  # type: ignore[method-assign]
    runner._audit_and_restore_processed_chains = lambda raw, processed: {"ok": True}  # type: ignore[method-assign]

    meta = runner._run_prepare(case, paths)

    cmd = captured[0]
    assert cmd[cmd.index("--audit_json") + 1] == str(paths.root / "prep_audit.json")
    assert cmd[cmd.index("--capping_mode") + 1] == "none"
    assert cmd[cmd.index("--neutralize") + 1] == "false"
    assert cmd[cmd.index("--ionic_strength_molar") + 1] == "0.15"
    assert cmd[cmd.index("--positive_ion") + 1] == "K+"
    assert cmd[cmd.index("--negative_ion") + 1] == "Br-"
    assert meta["prepare_policy"] == "permissive"
    assert meta["prepare_status"] == "ok"
    assert (paths.root / "prepare_report.json").exists()


def test_prepare_passes_capping_mode_to_openmm_prep(tmp_path: Path) -> None:
    cfg_path = make_cfg(tmp_path, allow_cpu_md=True)
    cfg = load_config(cfg_path, repo_root=REPO_ROOT)
    cfg.md.capping_mode = "termini_and_breaks"

    runner = HydroMapRunner(cfg, run_id="test")
    case = CaseSpec("tiny_protein", 123)
    paths = runner._case_paths(case)

    captured: list[list[str]] = []

    def fake_run_command(case_arg, stage, cmd, cwd):
        captured.append(cmd)
        paths.processed_pdb.parent.mkdir(parents=True, exist_ok=True)
        paths.processed_pdb.write_text("dummy", encoding="utf-8")
        paths.topology.parent.mkdir(parents=True, exist_ok=True)
        paths.topology.write_text("dummy", encoding="utf-8")
        (paths.root / "prep_audit.json").write_text("{\"warnings\": []}\n", encoding="utf-8")

    def fake_sanitize(src, dst):
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text("dummy", encoding="utf-8")
        return {"ok": True}

    runner._run_command = fake_run_command  # type: ignore[method-assign]
    runner._sanitize_input_pdb = fake_sanitize  # type: ignore[method-assign]
    runner._audit_and_restore_processed_chains = lambda raw, processed: {"ok": True}  # type: ignore[method-assign]

    meta = runner._run_prepare(case, paths)

    cmd = captured[0]
    assert cmd[cmd.index("--capping_mode") + 1] == "termini_and_breaks"
    assert meta["prepare_status"] == "ok"


def test_prepare_passes_histidine_controls_to_openmm_prep(tmp_path: Path) -> None:
    cfg_path = make_cfg(tmp_path, allow_cpu_md=True)
    cfg = load_config(cfg_path, repo_root=REPO_ROOT)
    cfg.md.histidine_mode = "hip"
    cfg.md.histidine_overrides = {"B:417": "HID", "505": "HIP"}

    runner = HydroMapRunner(cfg, run_id="test")
    case = CaseSpec("tiny_protein", 123)
    paths = runner._case_paths(case)

    captured: list[list[str]] = []

    def fake_run_command(case_arg, stage, cmd, cwd):
        captured.append(cmd)
        paths.processed_pdb.parent.mkdir(parents=True, exist_ok=True)
        paths.processed_pdb.write_text("dummy", encoding="utf-8")
        paths.topology.parent.mkdir(parents=True, exist_ok=True)
        paths.topology.write_text("dummy", encoding="utf-8")

    def fake_sanitize(src, dst):
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text("dummy", encoding="utf-8")
        return {"ok": True}

    runner._run_command = fake_run_command  # type: ignore[method-assign]
    runner._sanitize_input_pdb = fake_sanitize  # type: ignore[method-assign]
    runner._audit_and_restore_processed_chains = lambda raw, processed: {"ok": True}  # type: ignore[method-assign]

    runner._run_prepare(case, paths)

    cmd = captured[0]
    assert cmd[cmd.index("--histidine_mode") + 1] == "hip"
    override_pairs = [
        cmd[i + 1]
        for i, token in enumerate(cmd[:-1])
        if token == "--histidine_override"
    ]
    assert "B:417=HID" in override_pairs
    assert "505=HIP" in override_pairs


def test_prepare_uses_pdbfixer_repaired_temp_input_when_enabled(tmp_path: Path) -> None:
    cfg_path = make_cfg(tmp_path, allow_cpu_md=True)
    cfg = load_config(cfg_path, repo_root=REPO_ROOT)
    cfg.md.repair_missing_atoms = "pdbfixer"

    runner = HydroMapRunner(cfg, run_id="test")
    case = CaseSpec("tiny_protein", 123)
    paths = runner._case_paths(case)

    captured: list[list[str]] = []

    def fake_run_command(case_arg, stage, cmd, cwd):
        captured.append(cmd)
        paths.processed_pdb.parent.mkdir(parents=True, exist_ok=True)
        paths.processed_pdb.write_text("dummy", encoding="utf-8")
        paths.topology.parent.mkdir(parents=True, exist_ok=True)
        paths.topology.write_text("dummy", encoding="utf-8")

    def fake_sanitize(src, dst):
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text("raw", encoding="utf-8")
        return {
            "warnings": [],
            "remaining_incomplete_residues": ["LEU A:6"],
        }

    def fake_repair(src, dst):
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text("repaired", encoding="utf-8")
        return {"warnings": ["Applied pdbfixer"]}

    runner._run_command = fake_run_command  # type: ignore[method-assign]
    runner._sanitize_input_pdb = fake_sanitize  # type: ignore[method-assign]
    runner._repair_missing_atoms_with_pdbfixer = fake_repair  # type: ignore[method-assign]
    runner._audit_and_restore_processed_chains = lambda raw, processed: {"ok": True}  # type: ignore[method-assign]

    runner._run_prepare(case, paths)

    cmd = captured[0]
    assert cmd[1].endswith("prepare_with_openmm.py")
    assert cmd[2] == str(paths.simulation / "tiny_protein_prep_input.pdb")


def test_prepare_report_collects_warnings_in_permissive_mode(tmp_path: Path) -> None:
    cfg_path = make_cfg(tmp_path, allow_cpu_md=True)
    cfg = load_config(cfg_path, repo_root=REPO_ROOT)

    runner = HydroMapRunner(cfg, run_id="test")
    case = CaseSpec("tiny_protein", 123)
    paths = runner._case_paths(case)

    def fake_run_command(case_arg, stage, cmd, cwd):
        paths.processed_pdb.parent.mkdir(parents=True, exist_ok=True)
        paths.processed_pdb.write_text("dummy", encoding="utf-8")
        paths.topology.parent.mkdir(parents=True, exist_ok=True)
        paths.topology.write_text("dummy", encoding="utf-8")
        (paths.root / "prep_audit.json").write_text(
            '{"warnings":["Skipped one requested cap"],"skipped_caps":[{"residue":"A:7"}]}\n',
            encoding="utf-8",
        )

    def fake_sanitize(src, dst):
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text("dummy", encoding="utf-8")
        return {
            "warnings": ["Dropped incomplete terminal residues before preparation: GLY A:1"],
            "dropped_incomplete_terminal_residues": ["GLY A:1"],
            "dropped_incomplete_terminal_residue_count": 1,
            "remaining_incomplete_residues": [],
        }

    runner._run_command = fake_run_command  # type: ignore[method-assign]
    runner._sanitize_input_pdb = fake_sanitize  # type: ignore[method-assign]
    runner._audit_and_restore_processed_chains = lambda raw, processed: {"chain_restore_applied": False}  # type: ignore[method-assign]

    meta = runner._run_prepare(case, paths)
    report = (paths.root / "prepare_report.json").read_text(encoding="utf-8")

    assert meta["prepare_status"] == "warning"
    assert '"mode": "permissive"' in report
    assert '"would_block": false' in report
    assert "dropped_incomplete_terminal_residues" in report
    assert "skipped_caps_due_to_clash" in report


def test_prepare_policy_strict_blocks_warning_cases(tmp_path: Path) -> None:
    cfg_path = make_cfg(tmp_path, allow_cpu_md=True)
    cfg = load_config(cfg_path, repo_root=REPO_ROOT)
    cfg.md.prep_policy = "strict"

    runner = HydroMapRunner(cfg, run_id="test")
    case = CaseSpec("tiny_protein", 123)
    paths = runner._case_paths(case)

    def fake_run_command(case_arg, stage, cmd, cwd):
        paths.processed_pdb.parent.mkdir(parents=True, exist_ok=True)
        paths.processed_pdb.write_text("dummy", encoding="utf-8")
        paths.topology.parent.mkdir(parents=True, exist_ok=True)
        paths.topology.write_text("dummy", encoding="utf-8")
        (paths.root / "prep_audit.json").write_text('{"warnings":[]}\n', encoding="utf-8")

    def fake_sanitize(src, dst):
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text("dummy", encoding="utf-8")
        return {
            "warnings": ["Dropped incomplete terminal residues before preparation: GLY A:1"],
            "dropped_incomplete_terminal_residues": ["GLY A:1"],
            "dropped_incomplete_terminal_residue_count": 1,
            "remaining_incomplete_residues": [],
        }

    runner._run_command = fake_run_command  # type: ignore[method-assign]
    runner._sanitize_input_pdb = fake_sanitize  # type: ignore[method-assign]
    runner._audit_and_restore_processed_chains = lambda raw, processed: {"chain_restore_applied": False}  # type: ignore[method-assign]

    with pytest.raises(RuntimeError, match="Strict prepare policy blocked"):
        runner._run_prepare(case, paths)

    report = (paths.root / "prepare_report.json").read_text(encoding="utf-8")
    assert '"mode": "strict"' in report
    assert '"would_block": true' in report


def test_simulate_errors_when_gpu_missing_and_cpu_disallowed(tmp_path: Path) -> None:
    cfg = load_config(make_cfg(tmp_path, allow_cpu_md=False), repo_root=REPO_ROOT)
    runner = HydroMapRunner(cfg, run_id="test")
    case = CaseSpec("tiny_protein", 123)
    paths = runner._case_paths(case)

    paths.simulation.mkdir(parents=True, exist_ok=True)
    paths.processed_pdb.parent.mkdir(parents=True, exist_ok=True)
    paths.processed_pdb.write_text("dummy", encoding="utf-8")
    paths.topology.write_text("dummy", encoding="utf-8")

    runner._cuda_available = lambda: False  # type: ignore[method-assign]

    with pytest.raises(GPUNotAvailableError):
        runner._run_simulate(case, paths)


def test_analyze_is_chain_aware_by_default(tmp_path: Path) -> None:
    cfg = load_config(make_cfg(tmp_path, allow_cpu_md=True), repo_root=REPO_ROOT)
    runner = HydroMapRunner(cfg, run_id="test")
    case = CaseSpec("tiny_protein", 123)
    paths = runner._case_paths(case)

    paths.trajectory.parent.mkdir(parents=True, exist_ok=True)
    paths.trajectory.write_text("dummy", encoding="utf-8")
    paths.topology.parent.mkdir(parents=True, exist_ok=True)
    paths.topology.write_text("dummy", encoding="utf-8")
    paths.processed_pdb.parent.mkdir(parents=True, exist_ok=True)
    paths.processed_pdb.write_text("ATOM      1  N   ALA A   1      1.0 1.0 1.0 1.00 0.0           N\n", encoding="utf-8")
    paths.raw_pdb.parent.mkdir(parents=True, exist_ok=True)
    paths.raw_pdb.write_text("ATOM      1  N   ALA A   1      1.0 1.0 1.0 1.00 0.0           N\n", encoding="utf-8")

    captured: list[list[str]] = []

    def fake_run_command(case_arg, stage, cmd, cwd):
        captured.append(cmd)

    runner._run_command = fake_run_command  # type: ignore[method-assign]
    runner._trajectory_timing = lambda *_: (0.005, 1.0, 5)  # type: ignore[method-assign]
    meta = runner._run_analyze(case, paths)

    assert any("--multiChain" in cmd for cmd in captured)
    assert meta["potentials_device"] == "cpu"
    assert any("run_potentials_cpu.py" in " ".join(cmd) for cmd in captured)
    assert any("--nogpu" in cmd for cmd in captured if "run_potentials_cpu.py" in " ".join(cmd))


def test_analyze_can_use_gpu_devices(tmp_path: Path) -> None:
    cfg = load_config(
        make_cfg(
            tmp_path,
            allow_cpu_md=True,
            analysis_block=(
                "triplets_device: gpu\n"
                "potentials_device: gpu\n"
                "discard_initial_ns: 0.002\n"
                "triplets_frame_stride: 3\n"
                "potentials_frame_stride: 9"
            ),
        ),
        repo_root=REPO_ROOT,
    )
    runner = HydroMapRunner(cfg, run_id="test")
    case = CaseSpec("tiny_protein", 123)
    paths = runner._case_paths(case)

    paths.trajectory.parent.mkdir(parents=True, exist_ok=True)
    paths.trajectory.write_text("dummy", encoding="utf-8")
    paths.topology.parent.mkdir(parents=True, exist_ok=True)
    paths.topology.write_text("dummy", encoding="utf-8")
    paths.processed_pdb.parent.mkdir(parents=True, exist_ok=True)
    paths.processed_pdb.write_text("ATOM      1  N   ALA A   1      1.0 1.0 1.0 1.00 0.0           N\n", encoding="utf-8")
    paths.raw_pdb.parent.mkdir(parents=True, exist_ok=True)
    paths.raw_pdb.write_text("ATOM      1  N   ALA A   1      1.0 1.0 1.0 1.00 0.0           N\n", encoding="utf-8")

    captured: list[list[str]] = []

    def fake_run_command(case_arg, stage, cmd, cwd):
        captured.append(cmd)

    runner._run_command = fake_run_command  # type: ignore[method-assign]
    runner._triplet_gpu_available = lambda: True  # type: ignore[method-assign]
    runner._cuda_available = lambda: True  # type: ignore[method-assign]
    runner._trajectory_timing = lambda *_: (0.005, 1.0, 5)  # type: ignore[method-assign]
    meta = runner._run_analyze(case, paths)

    pot_cmds = [cmd for cmd in captured if "run_potentials_gpu.py" in " ".join(cmd)]
    tri_cmds = [cmd for cmd in captured if "run_triplets_gpu.py" in " ".join(cmd)]
    assert len(pot_cmds) == 1
    assert len(tri_cmds) == 1
    assert "--nogpu" not in pot_cmds[0]
    assert "--gpu" in tri_cmds[0]
    assert pot_cmds[0][pot_cmds[0].index("--skip") + 1] == "9"
    assert tri_cmds[0][tri_cmds[0].index("--skip") + 1] == "3"
    assert tri_cmds[0][tri_cmds[0].index("--frame-interval-ps") + 1] == "1.0"
    assert meta["triplets_device"] == "gpu"
    assert meta["potentials_device"] == "gpu"
    assert meta["discard_initial_ns"] == pytest.approx(0.002)


def test_analyze_can_skip_potentials_and_write_histograms(tmp_path: Path) -> None:
    cfg = load_config(
        make_cfg(
            tmp_path,
            allow_cpu_md=True,
            analysis_block=(
                "triplets_device: cpu\n"
                "potentials_device: gpu\n"
                "compute_potentials: false\n"
                "triplet_histogram_bin_width_deg: 20"
            ),
        ),
        repo_root=REPO_ROOT,
    )
    runner = HydroMapRunner(cfg, run_id="test")
    case = CaseSpec("tiny_protein", 123)
    paths = runner._case_paths(case)
    for path in (paths.trajectory, paths.processed_pdb, paths.raw_pdb):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("dummy", encoding="utf-8")

    captured: list[tuple[str, list[str]]] = []

    def fake_run_command(case_arg, stage, cmd, cwd):
        captured.append((stage, cmd))

    runner._run_command = fake_run_command  # type: ignore[method-assign]
    runner._cuda_available = lambda: False  # type: ignore[method-assign]
    runner._trajectory_timing = lambda *_: (0.005, 1.0, 5)  # type: ignore[method-assign]
    meta = runner._run_analyze(case, paths)

    stages = [stage for stage, _ in captured]
    assert "analyze_potentials" not in stages
    summary = next(cmd for stage, cmd in captured if stage == "summarize_triplets")
    assert summary[summary.index("--bin-width-deg") + 1] == "20.0"
    assert meta["compute_potentials"] is False
    assert meta["potentials_device"] is None

    with pytest.raises(RuntimeError, match="requires water-protein potentials"):
        runner._run_predict(case, paths)


def test_analyze_errors_if_triplet_gpu_device_requested_without_cupy(tmp_path: Path) -> None:
    cfg = load_config(
        make_cfg(
            tmp_path,
            allow_cpu_md=True,
            analysis_block="triplets_device: gpu",
        ),
        repo_root=REPO_ROOT,
    )
    runner = HydroMapRunner(cfg, run_id="test")
    case = CaseSpec("tiny_protein", 123)
    paths = runner._case_paths(case)

    paths.trajectory.parent.mkdir(parents=True, exist_ok=True)
    paths.trajectory.write_text("dummy", encoding="utf-8")
    paths.topology.parent.mkdir(parents=True, exist_ok=True)
    paths.topology.write_text("dummy", encoding="utf-8")
    paths.processed_pdb.parent.mkdir(parents=True, exist_ok=True)
    paths.processed_pdb.write_text("ATOM      1  N   ALA A   1      1.0 1.0 1.0 1.00 0.0           N\n", encoding="utf-8")
    paths.raw_pdb.parent.mkdir(parents=True, exist_ok=True)
    paths.raw_pdb.write_text("ATOM      1  N   ALA A   1      1.0 1.0 1.0 1.00 0.0           N\n", encoding="utf-8")

    runner._triplet_gpu_available = lambda: False  # type: ignore[method-assign]
    runner._trajectory_timing = lambda *_: (0.005, 1.0, 5)  # type: ignore[method-assign]
    with pytest.raises(RuntimeError, match="CuPy/CUDA"):
        runner._run_analyze(case, paths)


def test_analyze_errors_if_equilibration_consumes_all_trajectory(tmp_path: Path) -> None:
    cfg = load_config(
        make_cfg(
            tmp_path,
            allow_cpu_md=True,
            analysis_block="discard_initial_ns: 0.005",
        ),
        repo_root=REPO_ROOT,
    )
    runner = HydroMapRunner(cfg, run_id="test")
    case = CaseSpec("tiny_protein", 123)
    paths = runner._case_paths(case)

    paths.trajectory.parent.mkdir(parents=True, exist_ok=True)
    paths.trajectory.write_text("dummy", encoding="utf-8")
    paths.topology.parent.mkdir(parents=True, exist_ok=True)
    paths.topology.write_text("dummy", encoding="utf-8")
    paths.processed_pdb.parent.mkdir(parents=True, exist_ok=True)
    paths.processed_pdb.write_text("ATOM      1  N   ALA A   1      1.0 1.0 1.0 1.00 0.0           N\n", encoding="utf-8")
    paths.raw_pdb.parent.mkdir(parents=True, exist_ok=True)
    paths.raw_pdb.write_text("ATOM      1  N   ALA A   1      1.0 1.0 1.0 1.00 0.0           N\n", encoding="utf-8")

    runner._trajectory_timing = lambda *_: (0.005, 1.0, 5)  # type: ignore[method-assign]
    with pytest.raises(RuntimeError, match="discard_initial_ns"):
        runner._run_analyze(case, paths)
