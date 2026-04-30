from pathlib import Path
from types import SimpleNamespace

import pytest

from hydromap.selection import add_normalized_group_key_column, canonical_group_key, canonical_residue_key
from hydromap.selection import residue_chain_token, residue_sort_token_chain
from hydromap.utils.process_and_predict import build_residue_records, residue_scoped_path


def test_canonical_residue_key_distinguishes_icode() -> None:
    assert canonical_residue_key("A", 11, "") != canonical_residue_key("A", 11, "A")


def test_canonical_group_key_parses_resid_chain() -> None:
    sel = "resid 11A and segid P"
    assert canonical_group_key(sel) == "chain:P|resid:11|icode:A"


def test_canonical_group_key_chainid_precedence_over_segid() -> None:
    sel = "resid 11 and segid Z and chainID P"
    assert canonical_group_key(sel) == "chain:P|resid:11|icode:-"


def test_residue_chain_token_prefers_chainid_over_segid() -> None:
    residue = SimpleNamespace(
        resid=11,
        icode="A",
        segid="B",
        atoms=SimpleNamespace(chainIDs=["E", "E"]),
    )

    with pytest.warns(RuntimeWarning, match="conflicting PDB chainID and segid"):
        chain, alternates = residue_chain_token(residue)

    assert chain == "E"
    assert alternates == ["B"]
    sort_key, token, chain, alternates = residue_sort_token_chain(residue)
    assert sort_key == (11, "A")
    assert token == "11A"
    assert chain == "E"
    assert alternates == ["B"]


def test_residue_scoped_path_uses_canonical_chainid(tmp_path: Path) -> None:
    angles = tmp_path / "SARS2_res333_chainB_angles.txt"
    angles.write_text("90 100\n", encoding="utf-8")

    assert (
        residue_scoped_path(tmp_path, "SARS2", "333", "E", "angles.txt")
        == str(tmp_path / "SARS2_res333_chainE_angles.txt")
    )


def test_residue_chain_token_ignores_placeholder_system_segid() -> None:
    residue = SimpleNamespace(
        resid=332,
        icode="",
        segid="SYSTEM",
        atoms=SimpleNamespace(chainIDs=[]),
    )

    chain, alternates = residue_chain_token(residue)

    assert chain == ""
    assert alternates == []
    sort_key, token, chain, alternates = residue_sort_token_chain(residue)
    assert sort_key == (332, "")
    assert token == "332"
    assert chain == ""
    assert alternates == []


def test_residue_scoped_path_omits_empty_chain_suffix(tmp_path: Path) -> None:
    assert (
        residue_scoped_path(tmp_path, "alpha", "332", "", "angles.txt")
        == str(tmp_path / "alpha_res332_angles.txt")
    )


def test_process_and_predict_records_use_chainid_first() -> None:
    residue = SimpleNamespace(
        resid=333,
        icode="",
        segid="B",
        atoms=SimpleNamespace(chainIDs=["E"]),
    )
    universe = SimpleNamespace(residues=[residue])

    records = build_residue_records(universe, multi_chain=True)

    assert records == [
        {
            "sel": "resid 333 and chainID E",
            "token": "333",
            "chain": "E",
            "alternate_chains": ["B"],
        }
    ]


def test_canonical_group_key_parses_explicit_icode_clause() -> None:
    sel = "resid 11 and icode A and chainID P"
    assert canonical_group_key(sel) == "chain:P|resid:11|icode:A"


def test_add_normalized_group_key_column(tmp_path: Path) -> None:
    p = tmp_path / "results.csv"
    p.write_text(
        "MDAnalysis_selection_strings,Fdewet_pred\n"
        "resid 11 and segid P,1.0\n"
        "resid 11A and segid P,2.0\n",
        encoding="utf-8",
    )
    changed = add_normalized_group_key_column(p)
    text = p.read_text(encoding="utf-8")
    assert changed == 2
    assert "normalized_group_key" in text
    assert "chain:P|resid:11|icode:-" in text
    assert "chain:P|resid:11|icode:A" in text
