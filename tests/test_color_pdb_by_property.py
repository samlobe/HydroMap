from __future__ import annotations

import sys
from pathlib import Path

from hydromap.utils import color_pdb_by_property


def test_color_accepts_fdewet_alias_and_chainid_selection(tmp_path: Path, monkeypatch) -> None:
    pdb_path = tmp_path / "protein.pdb"
    csv_path = tmp_path / "results.csv"
    outdir = tmp_path / "colored"

    pdb_path.write_text(
        (
            "ATOM      1  N   ALA A  10      0.000   0.000   0.000  1.00  0.00           N\n"
            "ATOM      2  CA  ALA A  10      1.000   0.000   0.000  1.00  0.00           C\n"
            "ATOM      3  C   GLY A  11      2.000   0.000   0.000  1.00  0.00           C\n"
            "ATOM      4  O   GLY A  11      3.000   0.000   0.000  1.00  0.00           O\n"
            "END\n"
        ),
        encoding="utf-8",
    )
    csv_path.write_text(
        (
            "MDAnalysis_selection_strings,Fdewet_pred,avg_n_waters\n"
            "\"resid 10 and chainID A\",4.25,10\n"
            "\"resid 11 and chainID A\",5.50,12\n"
        ),
        encoding="utf-8",
    )

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "color_pdb_by_property.py",
            str(pdb_path),
            str(csv_path),
            "--outdir",
            str(outdir),
            "--properties",
            "Fdewet",
        ],
    )
    color_pdb_by_property.main()

    out_pdb = outdir / "protein_Fdewet_colored.pdb"
    assert out_pdb.exists()
    text = out_pdb.read_text(encoding="utf-8")
    assert "  4.25" in text
    assert "  5.50" in text
