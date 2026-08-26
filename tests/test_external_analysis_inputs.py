from __future__ import annotations

from pathlib import Path
from textwrap import dedent

import pytest

from hydromap.analysis_inputs import normalize_external_processed_pdb
from hydromap.config import load_config
from hydromap.runner import CaseSpec, HydroMapRunner


REPO_ROOT = Path(__file__).resolve().parents[1]


def _write_config(tmp_path: Path, analysis_block: str) -> Path:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture_dir = REPO_ROOT / "tests" / "fixtures"
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        f"""
input_dir: {fixture_dir}
protein: tiny_protein
seed: 123
model_path: {model}
artifacts_root: {tmp_path / 'artifacts'}
analysis:
{analysis_block}
""",
        encoding="utf-8",
    )
    return cfg


def _write_external_processed_pdb(path: Path, water_resname: str = "HOH") -> None:
    path.write_text(
        dedent(
            f"""
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1           1
ATOM      1  N   ALA A   1       8.000   8.000   8.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       9.300   8.400   8.200  1.00  0.00           C
ATOM      3  C   ALA A   1      10.100   7.300   8.900  1.00  0.00           C
ATOM      4  O   ALA A   1       9.900   6.100   8.700  1.00  0.00           O
ATOM      5  N   GLY A   2      11.100   7.700   9.700  1.00  0.00           N
ATOM      6  CA  GLY A   2      12.000   6.700  10.300  1.00  0.00           C
ATOM      7  C   GLY A   2      11.400   5.300  10.600  1.00  0.00           C
ATOM      8  O   GLY A   2      12.100   4.300  10.600  1.00  0.00           O
TER
HETATM    9  O   {water_resname:>3} A 101      14.000  14.000  14.000  1.00  0.00           O
HETATM   10  H1  {water_resname:>3} A 101      14.700  14.000  14.000  1.00  0.00           H
HETATM   11  H2  {water_resname:>3} A 101      13.500  14.700  14.000  1.00  0.00           H
TER
END
"""
        ).lstrip(),
        encoding="utf-8",
    )


def test_normalize_external_processed_pdb_rewrites_water_and_derives_raw(tmp_path: Path) -> None:
    src = tmp_path / "uploaded_processed.pdb"
    dst_processed = tmp_path / "normalized_processed.pdb"
    dst_raw = tmp_path / "normalized_raw.pdb"
    _write_external_processed_pdb(src, water_resname="HOH")

    summary = normalize_external_processed_pdb(src, dst_processed, dst_raw)

    normalized = dst_processed.read_text(encoding="utf-8")
    raw = dst_raw.read_text(encoding="utf-8")
    assert "SOL" in normalized
    assert " OW " in normalized
    assert "HOH" not in normalized
    assert "SOL" not in raw
    assert summary["water_residue_count"] == 1
    assert summary["water_oxygen_count"] == 1


def test_normalize_external_processed_pdb_requires_explicit_solvent(tmp_path: Path) -> None:
    src = tmp_path / "no_water.pdb"
    src.write_text((REPO_ROOT / "tests" / "fixtures" / "tiny_protein.pdb").read_text(encoding="utf-8"), encoding="utf-8")

    with pytest.raises(ValueError, match="explicit-solvent"):
        normalize_external_processed_pdb(src, tmp_path / "processed.pdb", tmp_path / "raw.pdb")


def test_analysis_inputs_accept_user_gromacs_topology(tmp_path: Path) -> None:
    uploaded_pdb = tmp_path / "uploaded_processed.pdb"
    uploaded_traj = tmp_path / "traj.dcd"
    uploaded_top = tmp_path / "topol.top"
    _write_external_processed_pdb(uploaded_pdb)
    uploaded_traj.write_bytes(b"dummy")
    uploaded_top.write_text("#include \"foo.itp\"\n", encoding="utf-8")

    cfg = load_config(
        _write_config(
            tmp_path,
            """
  triplets_device: cpu
  potentials_device: cpu
  existing_processed_pdb: {pdb}
  existing_trajectory: {traj}
  existing_topology: {top}
""".format(pdb=uploaded_pdb, traj=uploaded_traj, top=uploaded_top),
        ),
        repo_root=REPO_ROOT,
    )
    runner = HydroMapRunner(cfg, run_id="external-top")
    case = CaseSpec("tiny_protein", 123)
    paths = runner._case_paths(case)

    runner._ensure_analysis_inputs(case, paths)

    assert paths.topology == uploaded_top
    assert paths.processed_pdb.exists()
    assert paths.raw_pdb.exists()
    summary = (paths.root / "analysis_input_summary.json").read_text(encoding="utf-8")
    assert "user_gromacs_top" in summary


def test_analysis_inputs_build_xml_when_topology_omitted(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    uploaded_pdb = tmp_path / "uploaded_processed.pdb"
    uploaded_traj = tmp_path / "traj.dcd"
    _write_external_processed_pdb(uploaded_pdb)
    uploaded_traj.write_bytes(b"dummy")

    cfg = load_config(
        _write_config(
            tmp_path,
            f"""
  triplets_device: cpu
  potentials_device: cpu
  existing_processed_pdb: {uploaded_pdb}
  existing_trajectory: {uploaded_traj}
""",
        ),
        repo_root=REPO_ROOT,
    )
    runner = HydroMapRunner(cfg, run_id="external-fallback")
    case = CaseSpec("tiny_protein", 123)
    paths = runner._case_paths(case)

    called = {}

    def fake_build(processed_pdb: Path, topology_path: Path):
        called["processed_pdb"] = processed_pdb
        called["topology_path"] = topology_path
        topology_path.parent.mkdir(parents=True, exist_ok=True)
        topology_path.write_text("<System/>", encoding="utf-8")
        return {"topology_mode": "generated_xml", "topology_source": str(topology_path), "prediction_support": "a99SBdisp_supported"}

    monkeypatch.setattr(runner, "_build_analysis_topology_fallback", fake_build)

    runner._ensure_analysis_inputs(case, paths)

    assert called["processed_pdb"] == paths.processed_pdb
    assert called["topology_path"] == paths.topology
    assert paths.topology.read_text(encoding="utf-8") == "<System/>"


def test_triplet_only_external_analysis_does_not_build_topology(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    uploaded_pdb = tmp_path / "uploaded_processed.pdb"
    uploaded_traj = tmp_path / "traj.dcd"
    _write_external_processed_pdb(uploaded_pdb)
    uploaded_traj.write_bytes(b"dummy")

    cfg = load_config(
        _write_config(
            tmp_path,
            f"""
  triplets_device: cpu
  compute_potentials: false
  existing_processed_pdb: {uploaded_pdb}
  existing_trajectory: {uploaded_traj}
""",
        ),
        repo_root=REPO_ROOT,
    )
    runner = HydroMapRunner(cfg, run_id="external-triplets-only")
    case = CaseSpec("tiny_protein", 123)
    paths = runner._case_paths(case)

    def fail_if_called(*_args):
        raise AssertionError("triplet-only analysis must not build a potential topology")

    monkeypatch.setattr(runner, "_build_analysis_topology_fallback", fail_if_called)
    runner._ensure_analysis_inputs(case, paths)

    assert not paths.topology.exists()
    summary = (paths.root / "analysis_input_summary.json").read_text(encoding="utf-8")
    assert "not_required_triplet_only" in summary


def test_analysis_inputs_external_mode_does_not_require_source_pdb(tmp_path: Path) -> None:
    uploaded_pdb = tmp_path / "uploaded_processed.pdb"
    uploaded_traj = tmp_path / "traj.dcd"
    uploaded_top = tmp_path / "topol.top"
    _write_external_processed_pdb(uploaded_pdb)
    uploaded_traj.write_bytes(b"dummy")
    uploaded_top.write_text("#include \"foo.itp\"\n", encoding="utf-8")

    empty_input_dir = tmp_path / "empty_inputs"
    empty_input_dir.mkdir()
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(
        f"""
input_dir: {empty_input_dir}
protein: missing_source_name
seed: 123
model_path: {model}
artifacts_root: {tmp_path / 'artifacts'}
analysis:
  triplets_device: cpu
  potentials_device: cpu
  existing_processed_pdb: {uploaded_pdb}
  existing_trajectory: {uploaded_traj}
  existing_topology: {uploaded_top}
""",
        encoding="utf-8",
    )

    cfg = load_config(cfg_path, repo_root=REPO_ROOT)
    runner = HydroMapRunner(cfg, run_id="external-no-source")
    case = CaseSpec("missing_source_name", 123)
    paths = runner._case_paths(case)

    runner._ensure_analysis_inputs(case, paths)

    assert paths.raw_pdb.exists()
    assert paths.processed_pdb.exists()
    assert paths.trajectory.exists()
    summary = (paths.root / "analysis_input_summary.json").read_text(encoding="utf-8")
    assert "source_processed_pdb" in summary


def test_external_mode_does_not_require_existing_input_dir(tmp_path: Path) -> None:
    uploaded_pdb = tmp_path / "uploaded_processed.pdb"
    uploaded_traj = tmp_path / "traj.dcd"
    _write_external_processed_pdb(uploaded_pdb)
    uploaded_traj.write_bytes(b"dummy")

    missing_input_dir = tmp_path / "does_not_exist"
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(
        f"""
input_dir: {missing_input_dir}
protein: external_case
seed: 123
model_path: {model}
artifacts_root: {tmp_path / 'artifacts'}
analysis:
  triplets_device: cpu
  potentials_device: cpu
  existing_processed_pdb: {uploaded_pdb}
  existing_trajectory: {uploaded_traj}
""",
        encoding="utf-8",
    )

    cfg = load_config(cfg_path, repo_root=REPO_ROOT)
    assert cfg.input_dir == missing_input_dir
