from pathlib import Path

from hydromap.config import load_config
from hydromap.runner import HydroMapRunner


REPO_ROOT = Path(__file__).resolve().parents[1]


def _write_cfg(tmp_path: Path, strip_non_protein: bool = True) -> Path:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture = tmp_path / "in"
    fixture.mkdir(parents=True, exist_ok=True)
    (fixture / "tiny_protein.pdb").write_text(
        "CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1           1\n"
        "ATOM      1  N   ALA P  10      8.000   8.000   8.000  1.00  0.00           N\n"
        "HETATM    2  O   HOH P 501      1.000   1.000   1.000  1.00  0.00           O\n"
        "END\n",
        encoding="utf-8",
    )
    cfg = tmp_path / "cfg.yaml"
    cfg.write_text(
        f"""
input_dir: {fixture}
protein: tiny_protein
seed: 1
model_path: {model}
md:
  strip_non_protein: {'true' if strip_non_protein else 'false'}
""",
        encoding="utf-8",
    )
    return cfg


def test_strip_non_protein_default(tmp_path: Path) -> None:
    cfg = load_config(_write_cfg(tmp_path, strip_non_protein=True), repo_root=REPO_ROOT)
    runner = HydroMapRunner(cfg, run_id="t")
    src = cfg.input_dir / "tiny_protein.pdb"
    dst = tmp_path / "sanitized.pdb"
    summary = runner._sanitize_input_pdb(src, dst)
    text = dst.read_text(encoding="utf-8")
    assert "HETATM" not in text
    assert summary["removed_hetatm_lines"] == 1


def test_keep_non_protein_when_disabled(tmp_path: Path) -> None:
    cfg = load_config(_write_cfg(tmp_path, strip_non_protein=False), repo_root=REPO_ROOT)
    runner = HydroMapRunner(cfg, run_id="t")
    src = cfg.input_dir / "tiny_protein.pdb"
    dst = tmp_path / "kept.pdb"
    summary = runner._sanitize_input_pdb(src, dst)
    text = dst.read_text(encoding="utf-8")
    assert "HETATM" in text
    assert summary["strip_non_protein"] is False


def test_strip_non_protein_preserves_protein_like_caps_from_hetatm(tmp_path: Path) -> None:
    cfg = load_config(_write_cfg(tmp_path, strip_non_protein=True), repo_root=REPO_ROOT)
    runner = HydroMapRunner(cfg, run_id="t")

    src = tmp_path / "capped_input.pdb"
    src.write_text(
        "HETATM    1  C   ACE A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  O   ACE A   1       1.220   0.000   0.000  1.00  0.00           O\n"
        "HETATM    3  CH3 ACE A   1      -0.750   1.200   0.000  1.00  0.00           C\n"
        "ATOM      4  N   GLY A   2      -0.700  -1.150   0.000  1.00  0.00           N\n"
        "ATOM      5  CA  GLY A   2      -0.050  -2.440   0.000  1.00  0.00           C\n"
        "ATOM      6  C   GLY A   2       1.430  -2.250   0.000  1.00  0.00           C\n"
        "ATOM      7  O   GLY A   2       2.150  -3.220   0.000  1.00  0.00           O\n"
        "HETATM    8  O   HOH A 501       5.000   5.000   5.000  1.00  0.00           O\n"
        "END\n",
        encoding="utf-8",
    )

    dst = tmp_path / "sanitized_caps.pdb"
    summary = runner._sanitize_input_pdb(src, dst)
    text = dst.read_text(encoding="utf-8")

    assert "ATOM      1  C   ACE A   1" in text
    assert " HOH " not in text
    assert summary["preserved_protein_like_hetatm_lines"] == 3
    assert summary["removed_hetatm_lines"] == 1


def test_drop_incomplete_terminal_residues(tmp_path: Path) -> None:
    cfg = load_config(_write_cfg(tmp_path, strip_non_protein=True), repo_root=REPO_ROOT)
    runner = HydroMapRunner(cfg, run_id="t")

    src = tmp_path / "terminal_incomplete.pdb"
    src.write_text(
        "CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1           1\n"
        "ATOM      1  N   LEU A   4      0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  LEU A   4      1.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C   LEU A   4      2.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      4  O   LEU A   4      3.000   0.000   0.000  1.00  0.00           O\n"
        "ATOM      5  N   ALA A   5      4.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      6  CA  ALA A   5      5.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      7  C   ALA A   5      6.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      8  O   ALA A   5      7.000   0.000   0.000  1.00  0.00           O\n"
        "ATOM      9  CB  ALA A   5      5.000   1.000   0.000  1.00  0.00           C\n"
        "ATOM     10  N   VAL A   6      8.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM     11  CA  VAL A   6      9.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM     12  C   VAL A   6     10.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM     13  O   VAL A   6     11.000   0.000   0.000  1.00  0.00           O\n"
        "END\n",
        encoding="utf-8",
    )

    dst = tmp_path / "sanitized_terminal_incomplete.pdb"
    summary = runner._sanitize_input_pdb(src, dst)
    text = dst.read_text(encoding="utf-8")

    assert " LEU A   4" not in text
    assert " VAL A   6" not in text
    assert " ALA A   5" in text
    assert summary["dropped_incomplete_terminal_residue_count"] == 2
    assert summary["dropped_incomplete_terminal_residues"] == ["LEU A:4", "VAL A:6"]


def test_report_incomplete_internal_residue_warning(tmp_path: Path) -> None:
    cfg = load_config(_write_cfg(tmp_path, strip_non_protein=True), repo_root=REPO_ROOT)
    runner = HydroMapRunner(cfg, run_id="t")

    src = tmp_path / "internal_incomplete.pdb"
    src.write_text(
        "CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1           1\n"
        "ATOM      1  N   ALA A   5      4.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   5      5.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   5      6.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      4  O   ALA A   5      7.000   0.000   0.000  1.00  0.00           O\n"
        "ATOM      5  CB  ALA A   5      5.000   1.000   0.000  1.00  0.00           C\n"
        "ATOM      6  N   LEU A   6      8.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      7  CA  LEU A   6      9.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      8  C   LEU A   6     10.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      9  O   LEU A   6     11.000   0.000   0.000  1.00  0.00           O\n"
        "ATOM     10  N   GLY A   7     12.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM     11  CA  GLY A   7     13.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM     12  C   GLY A   7     14.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM     13  O   GLY A   7     15.000   0.000   0.000  1.00  0.00           O\n"
        "END\n",
        encoding="utf-8",
    )

    dst = tmp_path / "sanitized_internal_incomplete.pdb"
    summary = runner._sanitize_input_pdb(src, dst)

    assert summary["remaining_incomplete_residues"] == ["LEU A:6"]
    assert any("incomplete non-terminal residues" in msg for msg in summary["warnings"])


def test_restore_single_chain_id_in_processed_pdb(tmp_path: Path) -> None:
    cfg = load_config(_write_cfg(tmp_path, strip_non_protein=True), repo_root=REPO_ROOT)
    runner = HydroMapRunner(cfg, run_id="t")

    raw = tmp_path / "raw.pdb"
    raw.write_text(
        "ATOM      1  N   ALA P  10      0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  C   ALA P  10      1.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  N   GLY P  11A     2.000   0.000   0.000  1.00  0.00           N\n"
        "END\n",
        encoding="utf-8",
    )
    processed = tmp_path / "processed.pdb"
    processed.write_text(
        "ATOM      1  N   ALA    10      0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  C   ALA    10      1.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  N   GLY    11A     2.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      4  O   SOL    20      3.000   0.000   0.000  1.00  0.00           O\n"
        "END\n",
        encoding="utf-8",
    )

    audit = runner._audit_and_restore_processed_chains(raw, processed)
    text = processed.read_text(encoding="utf-8")
    assert audit["chain_restore_applied"] is True
    assert " ALA P  10" in text
    assert " GLY P  11A" in text


def test_restore_residue_ids_and_chain_from_raw_segid(tmp_path: Path) -> None:
    cfg = load_config(_write_cfg(tmp_path, strip_non_protein=True), repo_root=REPO_ROOT)
    runner = HydroMapRunner(cfg, run_id="t")

    raw = tmp_path / "raw_segid.pdb"
    raw.write_text(
        "ATOM      1  N   ACE   332      73.510  70.008  63.097  1.00  0.00      B    N  \n"
        "ATOM      2  N   THR   333      72.843  68.973  62.711  1.00  0.00      B    N  \n"
        "ATOM      3  N   ASN   334A     71.935  68.020  60.360  1.00  0.00      B    N  \n"
        "END\n",
        encoding="utf-8",
    )
    processed = tmp_path / "processed_segid.pdb"
    processed.write_text(
        (
            "HETATM    1  N   ACE A   1      0.000   0.000   0.000  1.00  0.00           N\n"
            "ATOM      2  N   THR A   2      1.000   0.000   0.000  1.00  0.00           N\n"
            "ATOM      3  N   ASN A   3      2.000   0.000   0.000  1.00  0.00           N\n"
            "HETATM    4  O   SOL C  10      3.000   0.000   0.000  1.00  0.00           O\n"
            "END\n"
        ),
        encoding="utf-8",
    )

    audit = runner._audit_and_restore_processed_chains(raw, processed)
    text = processed.read_text(encoding="utf-8")
    assert audit["chain_restore_applied"] is True
    assert audit["chain_restore_reason"] == "restored_protein_chain_and_residue_ids"
    assert " ACE B 332" in text
    assert " THR B 333" in text
    assert " ASN B 334A" in text


def test_keep_existing_multichain_ids(tmp_path: Path) -> None:
    cfg = load_config(_write_cfg(tmp_path, strip_non_protein=True), repo_root=REPO_ROOT)
    runner = HydroMapRunner(cfg, run_id="t")

    raw = tmp_path / "raw_multi.pdb"
    raw.write_text(
        "ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  N   GLY B   1      1.000   0.000   0.000  1.00  0.00           N\n"
        "END\n",
        encoding="utf-8",
    )
    processed = tmp_path / "processed_multi.pdb"
    original = (
        "ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  N   GLY B   1      1.000   0.000   0.000  1.00  0.00           N\n"
        "END\n"
    )
    processed.write_text(original, encoding="utf-8")

    audit = runner._audit_and_restore_processed_chains(raw, processed)
    text = processed.read_text(encoding="utf-8")
    assert audit["chain_restore_applied"] is False
    assert audit["chain_restore_reason"] == "no_atom_lines_changed"
    assert text == original
