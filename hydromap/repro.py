from __future__ import annotations

from pathlib import Path

from .config import HydroMapConfig
from .parity import run_parity_check


def run_repro_check(cfg: HydroMapConfig, run_id: str | None = None) -> tuple[Path, bool]:
    raise RuntimeError(
        "`hydromap repro check` was removed in v2. Use `hydromap parity baseline` and `hydromap parity check`."
    )
