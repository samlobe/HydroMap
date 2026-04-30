from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class CaseSpec:
    protein: str
    seed: int

    @property
    def case_key(self) -> str:
        return f"{self.protein}_seed{self.seed}"


@dataclass
class CasePaths:
    root: Path
    simulation: Path
    angles: Path
    potentials: Path
    results: Path
    colored: Path
    logs: Path
    raw_pdb: Path
    processed_pdb: Path
    topology: Path
    trajectory: Path
    energies_log: Path
    groups_file: Path | None


@dataclass
class CommandRecord:
    stage: str
    cwd: str
    command: str
    log: str
    exit_code: int


@dataclass
class StageResult:
    stage: str
    success: bool
    outputs: dict[str, str] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class RunRequest:
    run_id: str | None = None
    stages: list[str] | None = None
    proteins: list[str] | None = None
    seeds: list[int] | None = None


@dataclass
class RunResult:
    run_id: str
    manifest_paths: list[Path]


@dataclass
class ParityCaseMetrics:
    case_key: str
    matched_rows: int
    baseline_rows: int
    candidate_rows: int
    pearson: float | None
    mae: float | None
    rmse: float | None
    missing_keys: list[str] = field(default_factory=list)
    extra_keys: list[str] = field(default_factory=list)


@dataclass
class ParityReport:
    success: bool
    threshold: float
    pooled_pearson: float | None
    pooled_matched_rows: int
    cases: list[ParityCaseMetrics] = field(default_factory=list)
