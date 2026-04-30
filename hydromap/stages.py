from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from .types import StageResult

if TYPE_CHECKING:
    from .types import CasePaths, CaseSpec
    from .workflow import WorkflowRunner


@dataclass(frozen=True)
class Stage:
    name: str

    def run(self, runner: "WorkflowRunner", case: "CaseSpec", paths: "CasePaths") -> StageResult:
        raise NotImplementedError


class PrepareStage(Stage):
    def __init__(self) -> None:
        super().__init__(name="prepare")

    def run(self, runner: "WorkflowRunner", case: "CaseSpec", paths: "CasePaths") -> StageResult:
        meta = runner._run_prepare(case, paths)
        return StageResult(
            stage=self.name,
            success=True,
            outputs={"processed_pdb": str(paths.processed_pdb), "topology": str(paths.topology)},
            metadata=meta,
        )


class SimulateStage(Stage):
    def __init__(self) -> None:
        super().__init__(name="simulate")

    def run(self, runner: "WorkflowRunner", case: "CaseSpec", paths: "CasePaths") -> StageResult:
        meta = runner._run_simulate(case, paths)
        return StageResult(
            stage=self.name,
            success=True,
            outputs={"trajectory": str(paths.trajectory), "energies_log": str(paths.energies_log)},
            metadata=meta,
        )


class AnalyzeStage(Stage):
    def __init__(self) -> None:
        super().__init__(name="analyze")

    def run(self, runner: "WorkflowRunner", case: "CaseSpec", paths: "CasePaths") -> StageResult:
        meta = runner._run_analyze(case, paths)
        return StageResult(
            stage=self.name,
            success=True,
            outputs={"angles": str(paths.angles), "potentials": str(paths.potentials)},
            metadata=meta,
        )


class PredictStage(Stage):
    def __init__(self) -> None:
        super().__init__(name="predict")

    def run(self, runner: "WorkflowRunner", case: "CaseSpec", paths: "CasePaths") -> StageResult:
        meta = runner._run_predict(case, paths)
        return StageResult(
            stage=self.name,
            success=True,
            outputs={"results": str(paths.results / f"{case.protein}_results.csv")},
            metadata=meta,
        )


class ColorStage(Stage):
    def __init__(self) -> None:
        super().__init__(name="color")

    def run(self, runner: "WorkflowRunner", case: "CaseSpec", paths: "CasePaths") -> StageResult:
        runner._run_color(case, paths)
        return StageResult(
            stage=self.name,
            success=True,
            outputs={"colored": str(paths.colored)},
        )


def default_stage_registry() -> dict[str, Stage]:
    stages: list[Stage] = [PrepareStage(), SimulateStage(), AnalyzeStage(), PredictStage(), ColorStage()]
    return {s.name: s for s in stages}
