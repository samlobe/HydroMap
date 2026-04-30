from __future__ import annotations

"""Workflow runner aliases."""

from .types import CasePaths, CaseSpec, StageResult
from .workflow import STAGES, WorkflowRunner


HydroMapRunner = WorkflowRunner

__all__ = [
    "HydroMapRunner",
    "WorkflowRunner",
    "STAGES",
    "CaseSpec",
    "CasePaths",
    "StageResult",
]
