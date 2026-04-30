from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def maybe_hash(path: Path) -> str | None:
    if not path.exists() or not path.is_file():
        return None
    return sha256_file(path)


def summarize_results_csv(path: Path) -> dict[str, Any]:
    summary: dict[str, Any] = {
        "rows": 0,
        "fdewet_mean": None,
        "fdewet_min": None,
        "fdewet_max": None,
    }
    if not path.exists():
        return summary

    values: list[float] = []
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            summary["rows"] += 1
            val = row.get("Fdewet_pred")
            if val is None:
                continue
            try:
                values.append(float(val))
            except ValueError:
                continue

    if values:
        summary["fdewet_mean"] = sum(values) / len(values)
        summary["fdewet_min"] = min(values)
        summary["fdewet_max"] = max(values)
    return summary


def stable_payload_hash(payload: dict[str, Any]) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, sort_keys=True, indent=2)
        handle.write("\n")


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)
