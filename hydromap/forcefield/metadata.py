from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import json
from pathlib import Path


PACKAGE_DIR = Path(__file__).resolve().parent
A99SBDISP_METADATA_JSON = PACKAGE_DIR / "a99SBdisp_metadata.json"


@dataclass(frozen=True)
class ExplicitImproperDef:
    atom1: str
    atom2: str
    atom3: str
    atom4: str
    phase_degrees: float
    k_kj_per_mol: float
    periodicity: int


@dataclass(frozen=True)
class ResidueTemplateDef:
    name: str
    explicit_impropers: tuple[ExplicitImproperDef, ...]


@lru_cache(maxsize=1)
def load_residue_template_metadata() -> dict[str, ResidueTemplateDef]:
    raw = json.loads(A99SBDISP_METADATA_JSON.read_text(encoding="utf-8"))
    residues = raw.get("residues", {})
    return {
        name: ResidueTemplateDef(
            name=name,
            explicit_impropers=tuple(
                ExplicitImproperDef(
                    atom1=item["atom1"],
                    atom2=item["atom2"],
                    atom3=item["atom3"],
                    atom4=item["atom4"],
                    phase_degrees=float(item["phase_degrees"]),
                    k_kj_per_mol=float(item["k_kj_per_mol"]),
                    periodicity=int(item["periodicity"]),
                )
                for item in payload.get("explicit_impropers", [])
            ),
        )
        for name, payload in residues.items()
    }
