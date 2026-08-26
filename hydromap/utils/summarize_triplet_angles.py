#!/usr/bin/env python3
"""Write complete, auditable angle histograms from HydroMap raw triplet files."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


def _selections(groups_file: Path | None, count: int) -> list[str]:
    if groups_file is None:
        return [""] * count
    rows = [
        line.strip()
        for line in groups_file.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if len(rows) != count:
        raise ValueError(
            f"Found {count} group angle files but {len(rows)} selections in {groups_file}"
        )
    return rows


def summarize(
    angles_dir: Path,
    output_csv: Path,
    groups_file: Path | None = None,
    bin_width_deg: float = 10.0,
) -> None:
    if bin_width_deg <= 0.0 or not np.isclose(
        140.0 / bin_width_deg, round(140.0 / bin_width_deg)
    ):
        raise ValueError("bin_width_deg must be a positive divisor of 140 degrees")
    files = sorted(
        angles_dir.glob("*_angles.txt"),
        key=lambda path: (
            int(path.stem.rsplit("group", 1)[1].split("_", 1)[0])
            if "group" in path.stem
            else 10**9,
            path.name,
        ),
    )
    if not files:
        raise FileNotFoundError(f"No *_angles.txt files found in {angles_dir}")
    selections = _selections(groups_file, len(files))
    edges = np.arange(40.0, 180.0 + bin_width_deg / 2.0, bin_width_deg)
    bin_labels = [f"{edge:g}-{edges[i + 1]:g}" for i, edge in enumerate(edges[:-1])]
    fieldnames = [
        "group_index",
        "selection",
        "angle_file",
        "frame_count",
        "angle_count",
        "in_range_angle_count",
        "below_40_count",
        "above_180_count",
        *[f"{label}_count" for label in bin_labels],
        *[f"{label}_fraction" for label in bin_labels],
    ]
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for index, (path, selection) in enumerate(zip(files, selections), start=1):
            values: list[float] = []
            frame_count = 0
            with path.open(encoding="utf-8") as source:
                for line in source:
                    frame_count += 1
                    values.extend(float(token) for token in line.split())
            counts, _ = np.histogram(np.asarray(values), bins=edges)
            total = int(counts.sum())
            row: dict[str, object] = {
                "group_index": index,
                "selection": selection,
                "angle_file": str(path.resolve()),
                "frame_count": frame_count,
                "angle_count": len(values),
                "in_range_angle_count": total,
                "below_40_count": sum(value < 40.0 for value in values),
                "above_180_count": sum(value > 180.0 for value in values),
            }
            for label, count in zip(bin_labels, counts):
                row[f"{label}_count"] = int(count)
                row[f"{label}_fraction"] = (float(count) / total if total else 0.0)
            writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--angles-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--groups-file", type=Path)
    parser.add_argument("--bin-width-deg", type=float, default=10.0)
    args = parser.parse_args()
    summarize(
        args.angles_dir,
        args.output,
        groups_file=args.groups_file,
        bin_width_deg=args.bin_width_deg,
    )


if __name__ == "__main__":
    main()
