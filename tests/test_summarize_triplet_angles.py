import csv
from pathlib import Path

import pytest

from hydromap.utils.summarize_triplet_angles import summarize


def test_complete_histogram_preserves_group_mapping_and_outliers(tmp_path: Path) -> None:
    angles = tmp_path / "angles"
    angles.mkdir()
    (angles / "sample_group1_angles.txt").write_text(
        "39.9 40 49.9 50 179.9 180 180.1\n\n", encoding="utf-8"
    )
    groups = tmp_path / "groups.txt"
    groups.write_text(
        "# audited groups\nchainID B and resid 417 and not name H*\n",
        encoding="utf-8",
    )
    output = tmp_path / "histograms.csv"

    summarize(angles, output, groups_file=groups, bin_width_deg=10.0)

    with output.open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle))
    assert row["selection"] == "chainID B and resid 417 and not name H*"
    assert row["frame_count"] == "2"
    assert row["angle_count"] == "7"
    assert row["in_range_angle_count"] == "5"
    assert row["below_40_count"] == "1"
    assert row["above_180_count"] == "1"
    assert row["40-50_count"] == "2"
    assert row["50-60_count"] == "1"
    assert row["170-180_count"] == "2"
    fractions = [float(value) for key, value in row.items() if key.endswith("_fraction")]
    assert sum(fractions) == pytest.approx(1.0)


def test_histogram_rejects_bin_width_that_cannot_cover_range(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="positive divisor of 140"):
        summarize(tmp_path, tmp_path / "out.csv", bin_width_deg=12.0)
