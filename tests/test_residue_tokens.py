import pytest

from hydromap.potentials.potential import (
    parse_res_token as parse_potential_res_token,
    residue_selection as potential_residue_selection,
)
from hydromap.triplets.triplet import (
    parse_res_token as parse_triplet_res_token,
    residue_selection as triplet_residue_selection,
)


@pytest.mark.parametrize(
    ("token", "expected"),
    [("112", (112, "")), ("112b", (112, "B")), ("-2", (-2, "")), (" -2A ", (-2, "A"))],
)
def test_residue_token_parsers_support_pdb_residue_numbers(token, expected):
    assert parse_triplet_res_token(token) == expected
    assert parse_potential_res_token(token) == expected


def test_negative_residue_selections():
    assert triplet_residue_selection("A", "-2", heavy_only=True) == (
        "chainID A and resid -2 and not name H*",
        -2,
        "",
    )
    assert potential_residue_selection("A", "-2") == (
        "chainID A and resid -2",
        -2,
        "",
    )


def test_insertion_code_selections_are_preserved():
    assert triplet_residue_selection("A", "112b", heavy_only=True) == (
        "chainID A and resid 112B and not name H*",
        112,
        "B",
    )
    assert potential_residue_selection("A", "-2a") == (
        "chainID A and resid -2A",
        -2,
        "A",
    )


@pytest.mark.parametrize("token", ["", "A12", "1AB", "--2"])
def test_residue_token_parsers_reject_malformed_values(token):
    with pytest.raises(ValueError, match="Invalid residue token"):
        parse_triplet_res_token(token)
    with pytest.raises(ValueError, match="Invalid residue token"):
        parse_potential_res_token(token)
