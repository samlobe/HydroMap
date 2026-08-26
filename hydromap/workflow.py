from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import json
import math
import platform
import shutil
import subprocess
import sys
from typing import Iterable

from .analysis_inputs import normalize_external_processed_pdb
from .config import HydroMapConfig, config_to_dict
from .errors import GPUNotAvailableError, StageExecutionError
from .forcefield.openmm_utils import canonical_base_residues, canonicalize_residue_name
from .manifest import load_json, maybe_hash, stable_payload_hash, summarize_results_csv, write_json
from .resources import compute_cpu_workers
from .selection import add_normalized_group_key_column
from .stages import default_stage_registry
from .types import CasePaths, CaseSpec, CommandRecord, StageResult


STAGES = ["prepare", "simulate", "analyze", "predict", "color"]
GPU_POTENTIAL_GROUP_BATCH_SIZE = 10
NON_PROTEIN_RESNAMES = {
    "SOL",
    "HOH",
    "WAT",
    "TIP3",
    "TIP4",
    "TIP4P",
    "NA",
    "CL",
    "K",
    "CA",
    "MG",
}
HEAVY_ATOM_NAME_ALIASES = {
    ("ILE", "CD1"): "CD",
}


class WorkflowRunner:
    def __init__(self, cfg: HydroMapConfig, run_id: str | None = None) -> None:
        self.cfg = cfg
        self.repo_root = cfg.repo_root
        self.python = sys.executable

        timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        self.run_id = run_id or cfg.run_id or timestamp

        self.runs_root = cfg.artifacts_root / "runs" / self.run_id
        self.logs_root = cfg.artifacts_root / "logs" / self.run_id
        self.reports_root = cfg.artifacts_root / "reports"
        self.baselines_root = cfg.artifacts_root / "baselines"

        self.runs_root.mkdir(parents=True, exist_ok=True)
        self.logs_root.mkdir(parents=True, exist_ok=True)
        self.reports_root.mkdir(parents=True, exist_ok=True)
        self.baselines_root.mkdir(parents=True, exist_ok=True)

        self.command_history: dict[str, list[dict[str, object]]] = {}
        self.stage_results: dict[str, list[StageResult]] = {}
        self.stage_registry = default_stage_registry()

    def iter_cases(self, proteins: Iterable[str] | None = None, seeds: Iterable[int] | None = None) -> list[CaseSpec]:
        proteins = list(proteins) if proteins is not None else list(self.cfg.proteins)
        seeds = [int(x) for x in (list(seeds) if seeds is not None else self.cfg.seeds)]
        return [CaseSpec(protein=p, seed=s) for p in proteins for s in seeds]

    def run(
        self,
        stages: list[str] | None = None,
        proteins: Iterable[str] | None = None,
        seeds: Iterable[int] | None = None,
    ) -> list[Path]:
        stage_list = STAGES if stages is None else stages
        self._validate_stages(stage_list)

        manifests: list[Path] = []
        for case in self.iter_cases(proteins=proteins, seeds=seeds):
            manifest = self.run_case(case, stage_list)
            manifests.append(manifest)
        return manifests

    def run_case(self, case: CaseSpec, stages: list[str]) -> Path:
        paths = self._case_paths(case)
        self.command_history.setdefault(case.case_key, [])
        self.stage_results.setdefault(case.case_key, [])

        for stage_name in stages:
            stage = self.stage_registry[stage_name]
            result = stage.run(self, case, paths)
            self.stage_results[case.case_key].append(result)

        manifest_path = paths.root / "manifest.json"
        write_json(manifest_path, self._build_manifest(case, paths, stages))
        return manifest_path

    def _validate_stages(self, stages: list[str]) -> None:
        invalid = [s for s in stages if s not in STAGES]
        if invalid:
            raise ValueError(f"Unknown stages: {invalid}. Supported: {STAGES}")

    def _case_paths(self, case: CaseSpec) -> CasePaths:
        root = self.runs_root / case.protein / f"seed_{case.seed}"
        simulation = root / "simulation"
        angles = root / "angles"
        potentials = root / "potentials"
        results = root / "results"
        colored = root / "colored"
        logs = root / "logs"

        raw_pdb = root / f"{case.protein}.pdb"
        processed = root / f"{case.protein}_processed.pdb"
        topology = simulation / f"{case.protein}.xml"
        traj = simulation / f"{case.protein}_traj.dcd"
        energies = simulation / f"{case.protein}_energies.log"

        groups_file = None
        if self.cfg.groups_file is not None:
            groups_file = root / self.cfg.groups_file.name

        return CasePaths(
            root=root,
            simulation=simulation,
            angles=angles,
            potentials=potentials,
            results=results,
            colored=colored,
            logs=logs,
            raw_pdb=raw_pdb,
            processed_pdb=processed,
            topology=topology,
            trajectory=traj,
            energies_log=energies,
            groups_file=groups_file,
        )

    def _stage_log(self, case: CaseSpec, stage: str) -> Path:
        return self.logs_root / f"{case.case_key}_{stage}.log"

    def _source_pdb(self, protein: str) -> Path:
        return self.cfg.input_dir / f"{protein}.pdb"

    def _ensure_simulation_assets(self, simulation_dir: Path) -> None:
        simulation_dir.mkdir(parents=True, exist_ok=True)

    def _run_command(self, case: CaseSpec, stage: str, cmd: list[str], cwd: Path) -> None:
        cwd.mkdir(parents=True, exist_ok=True)
        log_path = self._stage_log(case, stage)
        log_path.parent.mkdir(parents=True, exist_ok=True)

        with log_path.open("w", encoding="utf-8") as handle:
            handle.write(f"cwd: {cwd}\n")
            handle.write(f"cmd: {' '.join(cmd)}\n\n")
            handle.flush()
            proc = subprocess.run(cmd, cwd=cwd, stdout=handle, stderr=subprocess.STDOUT)

        record = CommandRecord(
            stage=stage,
            cwd=str(cwd),
            command=" ".join(cmd),
            log=str(log_path),
            exit_code=int(proc.returncode),
        )
        self.command_history[case.case_key].append(record.__dict__)

        if proc.returncode != 0:
            raise StageExecutionError(f"Stage '{stage}' failed for {case.case_key}. See log: {log_path}")

    def _run_prepare(self, case: CaseSpec, paths: CasePaths) -> dict[str, object]:
        paths.root.mkdir(parents=True, exist_ok=True)
        paths.logs.mkdir(parents=True, exist_ok=True)
        self._ensure_simulation_assets(paths.simulation)

        src_pdb = self._source_pdb(case.protein)
        if not src_pdb.exists():
            raise FileNotFoundError(f"Missing source PDB: {src_pdb}")

        sanitize_summary = self._sanitize_input_pdb(src_pdb, paths.raw_pdb)
        prep_audit_path = paths.root / "prep_audit.json"
        topology_audit_path = paths.root / "topology_audit.json"
        write_json(paths.root / "input_sanitization.json", sanitize_summary)
        for warning in sanitize_summary.get("warnings", []):
            print(f"WARNING: {warning}", file=sys.stderr)

        prep_input_pdb = paths.raw_pdb
        prep_audit: dict[str, object] = {}
        topology_audit: dict[str, object] = {}

        try:
            remaining_incomplete = list(sanitize_summary.get("remaining_incomplete_residues", []))
            if remaining_incomplete:
                if self.cfg.md.repair_missing_atoms == "pdbfixer":
                    prep_input_pdb = paths.simulation / f"{case.protein}_prep_input.pdb"
                    repair_summary = self._repair_missing_atoms_with_pdbfixer(paths.raw_pdb, prep_input_pdb)
                    sanitize_summary["repair_summary"] = repair_summary
                    for warning in repair_summary.get("warnings", []):
                        print(f"WARNING: {warning}", file=sys.stderr)
                    write_json(paths.root / "input_sanitization.json", sanitize_summary)
                else:
                    residue_list = ", ".join(remaining_incomplete[:8])
                    if len(remaining_incomplete) > 8:
                        residue_list += ", ..."
                    raise RuntimeError(
                        "Detected incomplete non-terminal protein residues after sanitization: "
                        f"{residue_list}. HydroMap can only drop incomplete terminal residues by default. "
                        "To attempt repair of missing heavy atoms, set md.repair_missing_atoms: pdbfixer."
                    )

            shutil.copy2(prep_input_pdb, paths.simulation / f"{case.protein}.pdb")

            if self.cfg.groups_file is not None and paths.groups_file is not None:
                shutil.copy2(self.cfg.groups_file, paths.groups_file)

            preprocess_seed = self.cfg.md.preprocess_seed if self.cfg.md.preprocess_seed is not None else case.seed
            cmd = [
                self.python,
                str(self.repo_root / "hydromap" / "engines" / "simulation" / "prepare_with_openmm.py"),
                str(prep_input_pdb),
                "--output_pdb",
                str(paths.processed_pdb),
                "--output_system",
                str(paths.topology),
                "--audit_json",
                str(prep_audit_path),
                "--seed",
                str(preprocess_seed),
                "--capping_mode",
                self.cfg.md.capping_mode,
                "--neutralize",
                str(self.cfg.md.neutralize).lower(),
                "--ionic_strength_molar",
                str(self.cfg.md.ionic_strength_molar),
                "--positive_ion",
                self.cfg.md.positive_ion,
                "--negative_ion",
                self.cfg.md.negative_ion,
                "--histidine_mode",
                self.cfg.md.histidine_mode,
            ]
            for key, value in self.cfg.md.histidine_overrides.items():
                cmd.extend(["--histidine_override", f"{key}={value}"])
            self._run_command(case, "prepare", cmd, cwd=paths.simulation)

            if not paths.processed_pdb.exists():
                raise FileNotFoundError(f"Expected processed PDB missing: {paths.processed_pdb}")
            if not paths.topology.exists():
                raise FileNotFoundError(f"Expected topology missing: {paths.topology}")
            if prep_audit_path.exists():
                prep_audit = load_json(prep_audit_path)
                for warning in prep_audit.get("warnings", []):
                    print(f"WARNING: {warning}", file=sys.stderr)

            topology_audit = self._audit_and_restore_processed_chains(paths.raw_pdb, paths.processed_pdb)
            write_json(topology_audit_path, topology_audit)
            prepare_report = self._write_prepare_report(
                case=case,
                paths=paths,
                sanitize_summary=sanitize_summary,
                prep_audit=prep_audit,
                topology_audit=topology_audit,
            )
            self._enforce_prepare_policy(prepare_report, paths)

            sim_processed = paths.simulation / f"{case.protein}_processed.pdb"
            if sim_processed.exists() and sim_processed.resolve() != paths.processed_pdb.resolve():
                shutil.copy2(paths.processed_pdb, sim_processed)

            return {
                "prepare_report": str(paths.root / "prepare_report.json"),
                "prepare_status": str(prepare_report.get("status", "ok")),
                "prepare_warning_count": int(prepare_report.get("summary", {}).get("warning_count", 0)),
                "prepare_policy": self.cfg.md.prep_policy,
            }
        except Exception as exc:
            if not prep_audit and prep_audit_path.exists():
                prep_audit = load_json(prep_audit_path)
            if not topology_audit and topology_audit_path.exists():
                topology_audit = load_json(topology_audit_path)
            self._write_prepare_report(
                case=case,
                paths=paths,
                sanitize_summary=sanitize_summary,
                prep_audit=prep_audit,
                topology_audit=topology_audit,
                error_message=str(exc),
            )
            raise

    @staticmethod
    def _unique_messages(messages: list[str]) -> list[str]:
        ordered: list[str] = []
        seen: set[str] = set()
        for message in messages:
            if not message or message in seen:
                continue
            seen.add(message)
            ordered.append(message)
        return ordered

    def _build_prepare_report(
        self,
        case: CaseSpec,
        paths: CasePaths,
        sanitize_summary: dict[str, object],
        prep_audit: dict[str, object] | None = None,
        topology_audit: dict[str, object] | None = None,
        error_message: str | None = None,
    ) -> dict[str, object]:
        prep_audit = prep_audit or {}
        topology_audit = topology_audit or {}
        repair_summary = sanitize_summary.get("repair_summary") or {}
        repair_summary = repair_summary if isinstance(repair_summary, dict) else {}

        warnings = self._unique_messages(
            list(sanitize_summary.get("warnings", []))
            + list(repair_summary.get("warnings", []))
            + list(prep_audit.get("warnings", []))
        )
        if topology_audit.get("chain_restore_applied"):
            warnings.append(
                "Restored chain IDs and residue numbers onto the processed PDB after OpenMM preparation."
            )
        warnings = self._unique_messages(warnings)

        policy_triggers: list[str] = []
        if int(sanitize_summary.get("dropped_incomplete_terminal_residue_count", 0)) > 0:
            policy_triggers.append("dropped_incomplete_terminal_residues")
        if list(sanitize_summary.get("remaining_incomplete_residues", [])):
            policy_triggers.append("remaining_incomplete_residues")
        if repair_summary:
            policy_triggers.append("applied_missing_atom_repair")
        if list(prep_audit.get("skipped_caps", [])):
            policy_triggers.append("skipped_caps_due_to_clash")
        if topology_audit.get("chain_restore_applied"):
            policy_triggers.append("restored_chain_and_residue_ids")
        if warnings:
            policy_triggers.append("prepare_warnings")
        policy_triggers = self._unique_messages(policy_triggers)

        notes: list[str] = []
        removed_hetatm = int(sanitize_summary.get("removed_hetatm_lines", 0) or 0)
        if removed_hetatm > 0:
            notes.append(f"Removed {removed_hetatm} non-protein HETATM records during sanitization.")
        preserved_hetatm = int(sanitize_summary.get("preserved_protein_like_hetatm_lines", 0) or 0)
        if preserved_hetatm > 0:
            notes.append(f"Preserved {preserved_hetatm} protein-like HETATM records as ATOM records.")
        dropped_residues = list(sanitize_summary.get("dropped_incomplete_terminal_residues", []))
        if dropped_residues:
            notes.append(
                "Dropped incomplete terminal residues before preparation: " + ", ".join(dropped_residues)
            )
        if repair_summary:
            notes.append("Applied pdbfixer missing-atom repair to a temporary preparation input.")
        detected_breaks = int(prep_audit.get("detected_peptide_break_count", 0) or 0)
        if detected_breaks > 0:
            notes.append(f"Detected {detected_breaks} peptide break boundary/boundaries before peptide bonding.")
        added_ace = int(prep_audit.get("added_ace_cap_count", 0) or 0)
        added_nme = int(prep_audit.get("added_nme_cap_count", 0) or 0)
        if added_ace or added_nme:
            notes.append(f"Added caps during prep: ACE={added_ace}, NME={added_nme}.")
        skipped_caps = list(prep_audit.get("skipped_caps", []))
        if skipped_caps:
            notes.append(f"Skipped {len(skipped_caps)} requested cap placement(s) after clash screening.")
        if topology_audit.get("chain_restore_applied"):
            notes.append("Restored original chain IDs/residue numbers onto the processed protein records.")
        if error_message:
            notes.append(f"Preparation ended with an error: {error_message}")

        summary = {
            "removed_hetatm_lines": removed_hetatm,
            "preserved_protein_like_hetatm_lines": preserved_hetatm,
            "dropped_incomplete_terminal_residue_count": int(
                sanitize_summary.get("dropped_incomplete_terminal_residue_count", 0) or 0
            ),
            "remaining_incomplete_residue_count": len(list(sanitize_summary.get("remaining_incomplete_residues", []))),
            "detected_peptide_break_count": detected_breaks,
            "added_ace_cap_count": added_ace,
            "added_nme_cap_count": added_nme,
            "preserved_input_cap_count": len(list(prep_audit.get("preserved_input_caps", []))),
            "skipped_cap_count": len(skipped_caps),
            "chain_restore_applied": bool(topology_audit.get("chain_restore_applied", False)),
            "warning_count": len(warnings),
        }

        status = "error" if error_message else ("warning" if policy_triggers else "ok")
        return {
            "case": {
                "protein": case.protein,
                "seed": case.seed,
                "case_key": case.case_key,
            },
            "status": status,
            "policy": {
                "mode": self.cfg.md.prep_policy,
                "triggers": policy_triggers,
                "would_block": bool(policy_triggers) and self.cfg.md.prep_policy == "strict",
            },
            "inputs": {
                "source_pdb": str(paths.raw_pdb),
                "processed_pdb": str(paths.processed_pdb),
                "topology": str(paths.topology),
                "sanitization_report": str(paths.root / "input_sanitization.json"),
                "prep_audit": str(paths.root / "prep_audit.json"),
                "topology_audit": str(paths.root / "topology_audit.json"),
            },
            "summary": summary,
            "warnings": warnings,
            "notes": notes,
            "sanitization": sanitize_summary,
            "prep_audit": prep_audit,
            "topology_audit": topology_audit,
            "error": error_message,
        }

    def _write_prepare_report(
        self,
        case: CaseSpec,
        paths: CasePaths,
        sanitize_summary: dict[str, object],
        prep_audit: dict[str, object] | None = None,
        topology_audit: dict[str, object] | None = None,
        error_message: str | None = None,
    ) -> dict[str, object]:
        report = self._build_prepare_report(
            case=case,
            paths=paths,
            sanitize_summary=sanitize_summary,
            prep_audit=prep_audit,
            topology_audit=topology_audit,
            error_message=error_message,
        )
        write_json(paths.root / "prepare_report.json", report)
        return report

    def _enforce_prepare_policy(self, prepare_report: dict[str, object], paths: CasePaths) -> None:
        if self.cfg.md.prep_policy != "strict":
            return
        triggers = list(prepare_report.get("policy", {}).get("triggers", []))
        if not triggers:
            return
        trigger_text = ", ".join(triggers)
        raise RuntimeError(
            "Strict prepare policy blocked this structure after intake/prep review: "
            f"{trigger_text}. See {paths.root / 'prepare_report.json'} for details."
        )

    def _sanitize_input_pdb(self, src_pdb: Path, dst_pdb: Path) -> dict[str, object]:
        if not self.cfg.md.strip_non_protein:
            shutil.copy2(src_pdb, dst_pdb)
            return {
                "strip_non_protein": False,
                "source_pdb": str(src_pdb),
                "sanitized_pdb": str(dst_pdb),
                "removed_hetatm_lines": 0,
                "kept_atom_lines": None,
                "dropped_incomplete_terminal_residues": [],
                "dropped_incomplete_terminal_residue_count": 0,
                "remaining_incomplete_residues": [],
                "warnings": [],
            }

        lines = src_pdb.read_text(encoding="utf-8", errors="ignore").splitlines(keepends=True)
        out: list[str] = []
        removed_hetatm = 0
        preserved_hetatm = 0
        kept_atom = 0
        for line in lines:
            if line.startswith("HETATM"):
                if self._should_preserve_hetatm_line(line):
                    preserved_hetatm += 1
                    kept_atom += 1
                    out.append("ATOM  " + line[6:])
                    continue
                removed_hetatm += 1
                continue
            if line.startswith("ATOM"):
                kept_atom += 1
            out.append(line)

        if kept_atom == 0:
            raise RuntimeError(
                f"After stripping HETATM lines, no ATOM records remain in {src_pdb}. "
                "Set md.strip_non_protein=false if you intentionally need non-standard records."
            )

        out, dropped_residues = self._drop_incomplete_terminal_residues(out)
        remaining_incomplete = self._find_incomplete_protein_residues(out)
        kept_atom_after_drop = sum(1 for line in out if line.startswith("ATOM"))
        if kept_atom_after_drop == 0:
            raise RuntimeError(
                f"After removing incomplete terminal residues, no ATOM records remain in {src_pdb}."
            )

        dst_pdb.write_text("".join(out), encoding="utf-8")
        warnings: list[str] = []
        if dropped_residues:
            warnings.append(
                "Dropped incomplete terminal residues before preparation: "
                + ", ".join(dropped_residues)
            )
        if remaining_incomplete:
            residue_list = ", ".join(remaining_incomplete[:8])
            if len(remaining_incomplete) > 8:
                residue_list += ", ..."
            warnings.append(
                "Detected incomplete non-terminal residues that may require repair: "
                + residue_list
            )
        return {
            "strip_non_protein": True,
            "source_pdb": str(src_pdb),
            "sanitized_pdb": str(dst_pdb),
            "removed_hetatm_lines": removed_hetatm,
            "preserved_protein_like_hetatm_lines": preserved_hetatm,
            "kept_atom_lines": kept_atom_after_drop,
            "dropped_incomplete_terminal_residues": dropped_residues,
            "dropped_incomplete_terminal_residue_count": len(dropped_residues),
            "remaining_incomplete_residues": remaining_incomplete,
            "warnings": warnings,
        }

    @staticmethod
    def _should_preserve_hetatm_line(line: str) -> bool:
        resname = line[17:20].strip().upper()
        if resname in {"ACE", "NME", "NHE", "NH2"}:
            return True
        return canonicalize_residue_name(resname) in canonical_base_residues()

    @staticmethod
    def _residue_key_from_pdb_line(line: str) -> tuple[str, str, str, str]:
        chain = line[21].strip() or line[72:76].strip()
        resid = line[22:26].strip()
        icode = line[26].strip()
        resname = line[17:20].strip().upper()
        return chain, resid, icode, resname

    @staticmethod
    def _format_residue_key(key: tuple[str, str, str, str]) -> str:
        chain, resid, icode, resname = key
        label = f"{resname} {chain or '?'}:{resid}"
        if icode:
            label += icode
        return label

    @staticmethod
    def _expected_heavy_atoms_by_residue() -> dict[str, set[str]]:
        import xml.etree.ElementTree as ET

        from .forcefield import A99SBDISP_OPENMM_XML

        tree = ET.parse(A99SBDISP_OPENMM_XML)
        root = tree.getroot()

        type_elements: dict[str, str | None] = {}
        atom_types = root.find("AtomTypes")
        if atom_types is not None:
            for node in atom_types.findall("Type"):
                type_elements[node.attrib["name"]] = node.attrib.get("element")

        residues_node = root.find("Residues")
        expected: dict[str, set[str]] = {}
        if residues_node is None:
            return expected

        base_names = canonical_base_residues()
        for residue in residues_node.findall("Residue"):
            resname = residue.attrib["name"].upper()
            if resname not in base_names:
                continue
            heavy_atoms: set[str] = set()
            for atom in residue.findall("Atom"):
                atom_name = atom.attrib["name"].upper()
                atom_type = atom.attrib.get("type", "")
                element = type_elements.get(atom_type)
                if element == "H" or atom_name.startswith("H"):
                    continue
                heavy_atoms.add(atom_name)
            expected[resname] = heavy_atoms
        return expected

    @staticmethod
    def _normalize_heavy_atom_name(resname: str, atom_name: str) -> str:
        return HEAVY_ATOM_NAME_ALIASES.get((resname.upper(), atom_name.upper()), atom_name.upper())

    @classmethod
    def _drop_incomplete_terminal_residues(cls, lines: list[str]) -> tuple[list[str], list[str]]:
        expected_by_residue = cls._expected_heavy_atoms_by_residue()

        residue_atoms: dict[tuple[str, str, str, str], set[str]] = {}
        chain_residues: dict[str, list[tuple[str, str, str, str]]] = {}
        for line in lines:
            if not line.startswith("ATOM"):
                continue
            key = cls._residue_key_from_pdb_line(line)
            chain = key[0]
            if key not in residue_atoms:
                residue_atoms[key] = set()
                chain_residues.setdefault(chain, []).append(key)
            residue_atoms[key].add(cls._normalize_heavy_atom_name(key[3], line[12:16].strip()))

        def is_incomplete_terminal(key: tuple[str, str, str, str]) -> bool:
            resname = canonicalize_residue_name(key[3])
            expected = expected_by_residue.get(resname)
            if not expected:
                return False
            present = residue_atoms.get(key, set())
            return not expected.issubset(present)

        dropped: list[tuple[str, str, str, str]] = []
        for chain, residues in chain_residues.items():
            active = list(residues)
            while len(active) > 1 and is_incomplete_terminal(active[0]):
                dropped.append(active.pop(0))
            while len(active) > 1 and is_incomplete_terminal(active[-1]):
                dropped.append(active.pop())

        if not dropped:
            return lines, []

        dropped_keys = set(dropped)
        kept_lines = [
            line
            for line in lines
            if not (line.startswith("ATOM") and cls._residue_key_from_pdb_line(line) in dropped_keys)
        ]
        dropped_labels = [cls._format_residue_key(key) for key in dropped]
        return kept_lines, dropped_labels

    @classmethod
    def _find_incomplete_protein_residues(cls, lines: list[str]) -> list[str]:
        expected_by_residue = cls._expected_heavy_atoms_by_residue()
        residue_atoms: dict[tuple[str, str, str, str], set[str]] = {}
        residue_order: list[tuple[str, str, str, str]] = []
        seen: set[tuple[str, str, str, str]] = set()
        for line in lines:
            if not line.startswith("ATOM"):
                continue
            key = cls._residue_key_from_pdb_line(line)
            if key not in seen:
                seen.add(key)
                residue_order.append(key)
                residue_atoms[key] = set()
            residue_atoms[key].add(cls._normalize_heavy_atom_name(key[3], line[12:16].strip()))

        incomplete: list[str] = []
        for key in residue_order:
            resname = canonicalize_residue_name(key[3])
            expected = expected_by_residue.get(resname)
            if not expected:
                continue
            present = residue_atoms.get(key, set())
            if not expected.issubset(present):
                incomplete.append(cls._format_residue_key(key))
        return incomplete

    def _repair_missing_atoms_with_pdbfixer(self, src_pdb: Path, dst_pdb: Path) -> dict[str, object]:
        try:
            from pdbfixer import PDBFixer
            from openmm.app import PDBFile
        except ImportError as exc:
            raise RuntimeError(
                "md.repair_missing_atoms is set to pdbfixer, but pdbfixer is not installed. "
                "Install pdbfixer in the local/Docker HydroMap environment or set md.repair_missing_atoms: none."
            ) from exc

        fixer = PDBFixer(filename=str(src_pdb))
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

        dst_pdb.parent.mkdir(parents=True, exist_ok=True)
        with dst_pdb.open("w", encoding="utf-8") as handle:
            PDBFile.writeFile(fixer.topology, fixer.positions, handle)

        repaired_incomplete = self._find_incomplete_protein_residues(
            dst_pdb.read_text(encoding="utf-8", errors="ignore").splitlines(keepends=True)
        )
        warnings = [
            "Applied pdbfixer missing-atom repair to the preparation input. "
            "HydroMap will still restore chain IDs/residue numbers onto the processed output."
        ]
        if repaired_incomplete:
            residue_list = ", ".join(repaired_incomplete[:8])
            if len(repaired_incomplete) > 8:
                residue_list += ", ..."
            raise RuntimeError(
                "pdbfixer was unable to fully repair incomplete residues in the preparation input: "
                f"{residue_list}"
            )
        return {
            "repaired_pdb": str(dst_pdb),
            "missing_residue_ranges_detected": len(getattr(fixer, "missingResidues", {})),
            "warnings": warnings,
        }

    @staticmethod
    def _iter_residues(
        pdb_path: Path,
        atom_only: bool,
        exclude_non_protein: bool,
    ) -> list[tuple[str, str, str, str]]:
        residues: list[tuple[str, str, str, str]] = []
        seen: set[tuple[str, str, str, str]] = set()

        with pdb_path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                record = line[:6]
                if atom_only:
                    if record != "ATOM  ":
                        continue
                elif record not in {"ATOM  ", "HETATM"}:
                    continue

                resname = line[17:20].strip()
                if exclude_non_protein and resname in NON_PROTEIN_RESNAMES:
                    continue

                chain = line[21].strip()
                if not chain:
                    chain = line[72:76].strip()
                resid = line[22:26].strip()
                icode = line[26].strip()
                key = (chain, resid, icode, resname)
                if key not in seen:
                    seen.add(key)
                    residues.append(key)

        return residues

    def _audit_and_restore_processed_chains(self, raw_pdb: Path, processed_pdb: Path) -> dict[str, object]:
        raw_residues = self._iter_residues(raw_pdb, atom_only=True, exclude_non_protein=False)
        processed_before = self._iter_residues(processed_pdb, atom_only=False, exclude_non_protein=True)

        raw_chain_ids = sorted({r[0] for r in raw_residues if r[0].strip()})
        processed_chain_ids_before = sorted({r[0] for r in processed_before if r[0].strip()})
        same_residue_names_in_order = [r[3] for r in raw_residues] == [r[3] for r in processed_before]
        same_residues_ignoring_chain = [(r[1], r[2], r[3]) for r in raw_residues] == [
            (r[1], r[2], r[3]) for r in processed_before
        ]

        chain_restore_applied = False
        chain_restore_reason = "not_needed"
        changed_atom_lines = 0

        can_restore = (
            len(raw_residues) == len(processed_before)
            and len(processed_before) > 0
            and same_residue_names_in_order
        )

        if can_restore:
            residue_index: dict[tuple[str, str, str, str], int] = {}
            rewritten: list[str] = []
            next_idx = 0

            with processed_pdb.open("r", encoding="utf-8", errors="ignore") as handle:
                for line in handle:
                    record = line[:6]
                    if record in {"ATOM  ", "HETATM"}:
                        resname = line[17:20].strip()
                        if resname not in NON_PROTEIN_RESNAMES:
                            key = (line[21].strip(), line[22:26].strip(), line[26].strip(), resname)
                            if key not in residue_index:
                                residue_index[key] = next_idx
                                next_idx += 1

                            idx = residue_index[key]
                            if idx < len(raw_residues):
                                target_chain, target_resid, target_icode, _ = raw_residues[idx]
                                current_chain = line[21]
                                current_resid = line[22:26]
                                current_icode = line[26]
                                desired_chain = target_chain[:1] if target_chain else " "
                                desired_resid = target_resid.rjust(4)[:4]
                                desired_icode = target_icode[:1] if target_icode else " "
                                if (
                                    current_chain != desired_chain
                                    or current_resid != desired_resid
                                    or current_icode != desired_icode
                                ):
                                    line = f"{line[:21]}{desired_chain}{desired_resid}{desired_icode}{line[27:]}"
                                    changed_atom_lines += 1
                    rewritten.append(line)

            if changed_atom_lines > 0:
                processed_pdb.write_text("".join(rewritten), encoding="utf-8")
                chain_restore_applied = True
                chain_restore_reason = "restored_protein_chain_and_residue_ids"
            else:
                chain_restore_reason = "no_atom_lines_changed"
        else:
            if len(processed_before) == 0:
                chain_restore_reason = "no_protein_residues_in_processed"
            elif len(raw_residues) != len(processed_before):
                chain_restore_reason = "protein_residue_count_differs"
            elif not same_residue_names_in_order:
                chain_restore_reason = "protein_residue_names_differ"
            else:
                chain_restore_reason = "conditions_not_met"

        processed_after = self._iter_residues(processed_pdb, atom_only=False, exclude_non_protein=True)
        processed_chain_ids_after = sorted({r[0] for r in processed_after if r[0].strip()})

        return {
            "raw_pdb": str(raw_pdb),
            "processed_pdb": str(processed_pdb),
            "raw_protein_residue_count": len(raw_residues),
            "processed_protein_residue_count_before": len(processed_before),
            "processed_protein_residue_count_after": len(processed_after),
            "raw_chain_ids": raw_chain_ids,
            "processed_chain_ids_before": processed_chain_ids_before,
            "processed_chain_ids_after": processed_chain_ids_after,
            "same_residues_ignoring_chain": same_residues_ignoring_chain,
            "same_residue_names_in_order": same_residue_names_in_order,
            "chain_restore_applied": chain_restore_applied,
            "chain_restore_reason": chain_restore_reason,
            "changed_atom_lines": changed_atom_lines,
        }

    def _cuda_available(self) -> bool:
        try:
            from openmm import Platform

            Platform.getPlatformByName("CUDA")
            return True
        except Exception:
            return False

    def _triplet_gpu_available(self) -> bool:
        try:
            import cupy as cp

            return int(cp.cuda.runtime.getDeviceCount()) > 0
        except Exception:
            return False

    def _resolve_case_template_path(self, raw: str, case: CaseSpec) -> Path:
        try:
            rendered = raw.format(protein=case.protein, seed=case.seed)
        except Exception as exc:
            raise RuntimeError(
                f"Invalid case path template '{raw}'. "
                "Supported placeholders: {protein}, {seed}."
            ) from exc

        path = Path(rendered).expanduser()
        if not path.is_absolute():
            path = (self.cfg.config_path.parent / path).resolve()
        return path

    @staticmethod
    def _materialize_file(src: Path, dst: Path) -> None:
        if dst.exists():
            return
        dst.parent.mkdir(parents=True, exist_ok=True)
        try:
            dst.symlink_to(src)
        except OSError:
            shutil.copy2(src, dst)

    def _write_analysis_input_summary(self, paths: CasePaths, summary: dict[str, object]) -> None:
        summary_path = paths.root / "analysis_input_summary.json"
        with summary_path.open("w", encoding="utf-8") as handle:
            json.dump(summary, handle, indent=2, sort_keys=True)
            handle.write("\n")

    def _build_analysis_topology_fallback(self, processed_pdb: Path, topology_path: Path) -> dict[str, object]:
        if str(self.cfg.forcefield).strip().lower() != "a99sbdisp":
            raise RuntimeError(
                "analysis.existing_topology is required for external trajectories unless forcefield=a99SBdisp, "
                "because HydroMap can only rebuild a packaged OpenMM system for a99SBdisp."
            )
        from .forcefield import serialize_a99sbdisp_system_for_processed_pdb

        serialize_a99sbdisp_system_for_processed_pdb(processed_pdb, topology_path)
        return {
            "topology_mode": "generated_xml",
            "topology_source": str(topology_path),
            "prediction_support": "a99SBdisp_supported",
        }

    def _ensure_analysis_inputs(self, case: CaseSpec, paths: CasePaths) -> None:
        paths.root.mkdir(parents=True, exist_ok=True)
        paths.logs.mkdir(parents=True, exist_ok=True)
        paths.simulation.mkdir(parents=True, exist_ok=True)

        if self.cfg.groups_file is not None and paths.groups_file is not None and not paths.groups_file.exists():
            shutil.copy2(self.cfg.groups_file, paths.groups_file)

        ext_processed = self.cfg.analysis.existing_processed_pdb
        ext_traj = self.cfg.analysis.existing_trajectory
        ext_top = self.cfg.analysis.existing_topology
        if ext_processed and ext_traj:
            src_processed = self._resolve_case_template_path(ext_processed, case)
            src_traj = self._resolve_case_template_path(ext_traj, case)
            for src in (src_processed, src_traj):
                if not src.exists():
                    raise FileNotFoundError(f"Configured analysis existing input not found: {src}")

            normalization_summary = normalize_external_processed_pdb(
                src_processed,
                paths.processed_pdb,
                paths.raw_pdb,
            )
            self._materialize_file(src_traj, paths.trajectory)

            if ext_top:
                src_top = self._resolve_case_template_path(ext_top, case)
                if not src_top.exists():
                    raise FileNotFoundError(f"Configured analysis existing input not found: {src_top}")
                if src_top.suffix.lower() == ".top":
                    paths.topology = src_top
                    topology_summary = {
                        "topology_mode": "user_gromacs_top",
                        "topology_source": str(src_top),
                        "prediction_support": "a99SBdisp_only",
                    }
                else:
                    self._materialize_file(src_top, paths.topology)
                    topology_summary = {
                        "topology_mode": "user_openmm_xml",
                        "topology_source": str(src_top),
                        "prediction_support": "a99SBdisp_only",
                    }
            elif self.cfg.analysis.compute_potentials:
                topology_summary = self._build_analysis_topology_fallback(paths.processed_pdb, paths.topology)
            else:
                topology_summary = {
                    "topology_mode": "not_required_triplet_only",
                    "topology_source": None,
                    "prediction_support": "disabled_without_potentials",
                }

            self._write_analysis_input_summary(paths, {**normalization_summary, **topology_summary})
            return

        if not paths.raw_pdb.exists():
            src = self._source_pdb(case.protein)
            if not src.exists():
                raise FileNotFoundError(f"Missing source PDB: {src}")
            shutil.copy2(src, paths.raw_pdb)

    def _resolve_analysis_devices(self) -> tuple[str, str]:
        triplets_device = self.cfg.analysis.triplets_device
        potentials_device = self.cfg.analysis.potentials_device

        if triplets_device == "auto":
            triplets_device = "gpu" if self._triplet_gpu_available() else "cpu"
        if self.cfg.analysis.compute_potentials and potentials_device == "auto":
            potentials_device = "gpu" if self._cuda_available() else "cpu"
        return triplets_device, potentials_device

    @staticmethod
    def _resolve_stride(frame_dt_ps: float, frame_stride: int, sample_ps: float | None) -> int:
        if sample_ps is None:
            return frame_stride
        if frame_dt_ps <= 0.0:
            raise RuntimeError("Trajectory dt must be > 0 to use analysis.*_sample_ps controls.")
        return max(1, int(math.ceil(sample_ps / frame_dt_ps)))

    def _trajectory_timing(self, processed_pdb: Path, trajectory: Path) -> tuple[float, float, int]:
        try:
            import MDAnalysis as mda
        except Exception as exc:
            raise RuntimeError("MDAnalysis is required to inspect trajectory timing.") from exc

        u = mda.Universe(str(processed_pdb), str(trajectory))
        n_frames = len(u.trajectory)
        if n_frames <= 0:
            raise RuntimeError("Trajectory has zero frames.")
        dt_ps = float(u.trajectory.dt)
        if dt_ps <= 0.0:
            raise RuntimeError("Trajectory dt must be > 0 ps.")
        total_ns = (n_frames * dt_ps) / 1000.0
        return total_ns, dt_ps, n_frames

    def _run_simulate(self, case: CaseSpec, paths: CasePaths) -> dict[str, object]:
        if not paths.processed_pdb.exists() or not paths.topology.exists():
            raise RuntimeError("Missing prepared inputs. Run `prepare` first.")

        sim_seed = self.cfg.md.random_seed if self.cfg.md.random_seed is not None else case.seed
        vel_seed = self.cfg.md.velocity_seed if self.cfg.md.velocity_seed is not None else sim_seed
        bar_seed = self.cfg.md.barostat_seed if self.cfg.md.barostat_seed is not None else sim_seed

        md_device = self.cfg.md.device
        if md_device == "auto":
            md_device = "gpu" if self._cuda_available() else "cpu"

        use_no_cuda = (md_device == "cpu")
        cpu_fallback = False

        if use_no_cuda and not self.cfg.md.allow_cpu_md:
            raise GPUNotAvailableError(
                "CPU MD requested (md.device=cpu/auto) but disallowed by policy (md.allow_cpu_md=false)."
            )

        if not use_no_cuda and not self._cuda_available():
            if not self.cfg.md.allow_cpu_md:
                raise GPUNotAvailableError(
                    "CUDA platform unavailable and CPU MD is disallowed (md.allow_cpu_md=false)."
                )
            use_no_cuda = True
            cpu_fallback = True

        cmd = [
            self.python,
            str(self.repo_root / "hydromap" / "engines" / "simulation" / "simulate_with_openmm.py"),
            f"../{case.protein}_processed.pdb",
            paths.topology.name,
            "-ns",
            str(self.cfg.md.nanoseconds),
            "--random_seed",
            str(sim_seed),
            "--velocity_seed",
            str(vel_seed),
            "--barostat_seed",
            str(bar_seed),
            "--checkpoint_policy",
            self.cfg.md.checkpoint_policy,
            "--cuda_precision",
            self.cfg.md.cuda_precision,
            "--equilibration_ps",
            str(self.cfg.md.equilibration_ns * 1000.0),
            "--equilibration_protocol",
            self.cfg.md.equilibration_protocol,
            "--timestep_ps",
            str(self.cfg.md.timestep_ps),
            "--report_interval_ps",
            str(self.cfg.md.report_interval_ps),
            "--checkpoint_interval_ps",
            str(self.cfg.md.checkpoint_interval_ps),
            "-o",
            f"{case.protein}_traj.dcd",
        ]

        if self.cfg.md.restrain is not None:
            if self.cfg.md.restrain == "":
                cmd.append("-r")
            else:
                cmd.extend(["-r", self.cfg.md.restrain])
            cmd.extend(["--restraint_k", str(self.cfg.md.restraint_k)])

        if use_no_cuda:
            cmd.append("--noCUDA")

        if self.cfg.md.deterministic:
            cmd.append("--deterministic")
        if self.cfg.md.constant_volume:
            cmd.append("--nvt")
        if self.cfg.md.initial_state is not None:
            initial_state = self._resolve_case_template_path(self.cfg.md.initial_state, case)
            if not initial_state.is_file():
                raise FileNotFoundError(f"Initial OpenMM state not found: {initial_state}")
            cmd.extend(["--initial_state", str(initial_state)])

        self._run_command(case, "simulate", cmd, cwd=paths.simulation)

        if not paths.trajectory.exists():
            raise FileNotFoundError(f"Expected trajectory missing: {paths.trajectory}")

        return {
            "cpu_md_fallback": cpu_fallback,
            "used_no_cuda": use_no_cuda,
            "md_device": md_device,
            "execution_profile": self.cfg.execution.profile,
        }

    def _resolve_analysis_windows(self, total_ns: float) -> dict[str, float]:
        discard_ns = float(self.cfg.analysis.discard_initial_ns)
        available_ns = total_ns - discard_ns
        if available_ns <= 0:
            raise RuntimeError(
                f"analysis.discard_initial_ns ({discard_ns}) must be smaller than available trajectory "
                f"time ({total_ns:.4f} ns)."
            )

        base_tail_ns = self.cfg.analysis.tail_ns if self.cfg.analysis.tail_ns is not None else available_ns
        triplets_tail = (
            self.cfg.analysis.triplets_tail_ns if self.cfg.analysis.triplets_tail_ns is not None else base_tail_ns
        )
        potentials_tail = (
            self.cfg.analysis.potentials_tail_ns if self.cfg.analysis.potentials_tail_ns is not None else base_tail_ns
        )

        triplets_tail = max(0.0, min(float(triplets_tail), available_ns))
        potentials_tail = max(0.0, min(float(potentials_tail), available_ns))
        if triplets_tail <= 0.0:
            raise RuntimeError("Resolved triplets analysis window is <= 0 ns. Check analysis time-window settings.")
        if self.cfg.analysis.compute_potentials and potentials_tail <= 0.0:
            raise RuntimeError("Resolved potentials analysis window is <= 0 ns. Check analysis time-window settings.")

        return {
            "discard_initial_ns": discard_ns,
            "available_ns": available_ns,
            "triplets_tail_ns": triplets_tail,
            "potentials_tail_ns": potentials_tail,
        }

    def _run_analyze(self, case: CaseSpec, paths: CasePaths) -> dict[str, object]:
        self._ensure_analysis_inputs(case, paths)
        if not paths.trajectory.exists():
            raise RuntimeError(
                "Missing trajectory. Run `simulate` first or set analysis.existing_processed_pdb and analysis.existing_trajectory."
            )
        if self.cfg.analysis.compute_potentials and not paths.topology.exists():
            raise RuntimeError(
                "Missing topology. Run `prepare`/`simulate` first, provide analysis.existing_topology, "
                "or use forcefield=a99SBdisp so HydroMap can rebuild a fallback XML."
            )
        if not paths.processed_pdb.exists():
            raise RuntimeError(
                "Missing processed PDB. Run `prepare` first or set analysis.existing_* inputs."
            )

        traj_total_ns, frame_dt_ps, n_frames = self._trajectory_timing(paths.processed_pdb, paths.trajectory)
        windows = self._resolve_analysis_windows(traj_total_ns)
        workers = compute_cpu_workers(
            max_cpu_workers=self.cfg.resources.max_cpu_workers,
            reserve_cpus=self.cfg.resources.reserve_cpus,
            profile=self.cfg.execution.profile,
        )
        triplet_time_ns = float(windows["triplets_tail_ns"])
        potential_time_ns = float(windows["potentials_tail_ns"])
        triplets_skip = self._resolve_stride(
            frame_dt_ps=frame_dt_ps,
            frame_stride=self.cfg.analysis.triplets_frame_stride,
            sample_ps=self.cfg.analysis.triplets_sample_ps,
        )
        potentials_skip = None
        if self.cfg.analysis.compute_potentials:
            potentials_skip = self._resolve_stride(
                frame_dt_ps=frame_dt_ps,
                frame_stride=self.cfg.analysis.potentials_frame_stride,
                sample_ps=self.cfg.analysis.potentials_sample_ps,
            )

        triplets_device, potentials_device = self._resolve_analysis_devices()
        if triplets_device == "gpu" and not self._triplet_gpu_available():
            raise RuntimeError(
                "analysis.triplets_device=gpu requested but CuPy/CUDA is unavailable. "
                "Install CuPy for your CUDA version or use triplets_device=cpu."
            )
        if (
            self.cfg.analysis.compute_potentials
            and potentials_device == "gpu"
            and not self._cuda_available()
        ):
            raise RuntimeError(
                "analysis.potentials_device=gpu requested but OpenMM CUDA platform is unavailable. "
                "Use potentials_device=cpu on this machine."
            )

        if triplets_device == "cpu":
            triplets_cmd = [
                self.python,
                str(self.repo_root / "hydromap" / "engines" / "triplets" / "run_triplets_cpu.py"),
                str(paths.raw_pdb),
                str(paths.trajectory),
                "--nprocs",
                str(workers),
                "-t",
                str(triplet_time_ns),
                "--skip",
                str(triplets_skip),
                "--hydrationCutoff",
                str(self.cfg.analysis.triplets_hydration_cutoff),
                "--outdir",
                str(paths.angles),
            ]
        else:
            triplets_cmd = [
                self.python,
                str(self.repo_root / "hydromap" / "engines" / "triplets" / "run_triplets_gpu.py"),
                str(paths.raw_pdb),
                str(paths.trajectory),
                "-t",
                str(triplet_time_ns),
                "--skip",
                str(triplets_skip),
                "--frame-interval-ps",
                str(frame_dt_ps),
                "--hydrationCutoff",
                str(self.cfg.analysis.triplets_hydration_cutoff),
                "--outdir",
                str(paths.angles),
                "--gpu",
            ]

        potentials_cmd = None
        if self.cfg.analysis.compute_potentials and potentials_device == "cpu":
            potentials_cmd = [
                self.python,
                str(self.repo_root / "hydromap" / "engines" / "potentials" / "run_potentials_cpu.py"),
                str(paths.raw_pdb),
                str(paths.trajectory),
                "--top",
                str(paths.topology),
                "--nprocs",
                str(workers),
                "-t",
                str(int(max(1, math.ceil(potential_time_ns)))),
                "--skip",
                str(potentials_skip),
                "--cutoff",
                str(self.cfg.analysis.potentials_cutoff),
                "--outdir",
                str(paths.potentials),
                "--nogpu",
            ]
        elif self.cfg.analysis.compute_potentials:
            potentials_cmd = [
                self.python,
                str(self.repo_root / "hydromap" / "engines" / "potentials" / "run_potentials_gpu.py"),
                str(paths.raw_pdb),
                str(paths.trajectory),
                "--top",
                str(paths.topology),
                "-t",
                str(potential_time_ns),
                "--skip",
                str(potentials_skip),
                "--cutoff",
                str(self.cfg.analysis.potentials_cutoff),
                "--outdir",
                str(paths.potentials),
                "--groupBatchSize",
                str(GPU_POTENTIAL_GROUP_BATCH_SIZE),
            ]

        if paths.groups_file is not None:
            triplets_cmd.extend(["--groupsFile", str(paths.groups_file)])
            if potentials_cmd is not None:
                potentials_cmd.extend(["--groupsFile", str(paths.groups_file)])
        else:
            # Default v2 behavior: always chain-aware residue mode.
            triplets_cmd.append("--multiChain")
            if potentials_cmd is not None:
                potentials_cmd.append("--multiChain")

        self._run_command(case, "analyze_triplets", triplets_cmd, cwd=paths.root)
        histogram_cmd = [
            self.python,
            str(self.repo_root / "hydromap" / "utils" / "summarize_triplet_angles.py"),
            "--angles-dir",
            str(paths.angles),
            "--output",
            str(paths.results / f"{case.protein}_triplet_histograms.csv"),
            "--bin-width-deg",
            str(self.cfg.analysis.triplet_histogram_bin_width_deg),
        ]
        if paths.groups_file is not None:
            histogram_cmd.extend(["--groups-file", str(paths.groups_file)])
        self._run_command(case, "summarize_triplets", histogram_cmd, cwd=paths.root)
        if potentials_cmd is not None:
            self._run_command(case, "analyze_potentials", potentials_cmd, cwd=paths.root)
        return {
            "workers": workers,
            "triplets_device": triplets_device,
            "potentials_device": (
                potentials_device if self.cfg.analysis.compute_potentials else None
            ),
            "compute_potentials": self.cfg.analysis.compute_potentials,
            "triplet_histogram_csv": str(
                paths.results / f"{case.protein}_triplet_histograms.csv"
            ),
            "discard_initial_ns": windows["discard_initial_ns"],
            "available_ns": windows["available_ns"],
            "triplets_tail_ns": windows["triplets_tail_ns"],
            "potentials_tail_ns": windows["potentials_tail_ns"],
            "frame_dt_ps": frame_dt_ps,
            "trajectory_total_ns": traj_total_ns,
            "trajectory_frames": n_frames,
            "triplets_skip": triplets_skip,
            "potentials_skip": potentials_skip,
        }

    def _run_predict(self, case: CaseSpec, paths: CasePaths) -> dict[str, object]:
        if not self.cfg.analysis.compute_potentials:
            raise RuntimeError(
                "The bundled Fdewet predictor requires water-protein potentials. "
                "Set analysis.compute_potentials=true or run only through the analyze stage."
            )
        paths.results.mkdir(parents=True, exist_ok=True)
        predict_cmd = [
            self.python,
            str(self.repo_root / "hydromap" / "utils" / "process_and_predict.py"),
            str(paths.raw_pdb),
            "--anglesDir",
            str(paths.angles),
            "--potentialsDir",
            str(paths.potentials),
            "--model",
            str(self.cfg.model_path),
            "--forcefield",
            self.cfg.forcefield,
            "--outdir",
            str(paths.results),
        ]

        if paths.groups_file is not None:
            predict_cmd.extend(["--groupsFile", str(paths.groups_file)])
        else:
            predict_cmd.append("--multiChain")

        self._run_command(case, "predict", predict_cmd, cwd=self.repo_root)

        results_csv = paths.results / f"{case.protein}_results.csv"
        if not results_csv.exists():
            raise FileNotFoundError(f"Missing results CSV: {results_csv}")

        changed_rows = add_normalized_group_key_column(results_csv)
        return {"normalized_group_key_rows_changed": changed_rows}

    def _run_color(self, case: CaseSpec, paths: CasePaths) -> None:
        paths.colored.mkdir(parents=True, exist_ok=True)
        results_csv = paths.results / f"{case.protein}_results.csv"
        if not results_csv.exists():
            raise RuntimeError(f"Missing results CSV: {results_csv}. Run `predict` first.")

        color_cmd = [
            self.python,
            str(self.repo_root / "hydromap" / "utils" / "color_pdb_by_property.py"),
            str(paths.raw_pdb),
            str(results_csv),
            "--outdir",
            str(paths.colored),
            "--minWaters",
            str(self.cfg.analysis.min_waters),
        ]

        if self.cfg.analysis.color_properties:
            color_cmd.extend(["--properties", *self.cfg.analysis.color_properties])

        self._run_command(case, "color", color_cmd, cwd=paths.root)

    def _build_manifest(self, case: CaseSpec, paths: CasePaths, stages: list[str]) -> dict[str, object]:
        results_csv = paths.results / f"{case.protein}_results.csv"
        colored_files = sorted(paths.colored.glob("*_colored.pdb"))

        files = {
            "raw_pdb": str(paths.raw_pdb),
            "processed_pdb": str(paths.processed_pdb),
            "topology": str(paths.topology),
            "trajectory": str(paths.trajectory),
            "energies_log": str(paths.energies_log),
            "results_csv": str(results_csv),
            "groups_file": (str(paths.groups_file) if paths.groups_file else None),
            "input_sanitization": str(paths.root / "input_sanitization.json"),
            "prep_audit": str(paths.root / "prep_audit.json"),
            "topology_audit": str(paths.root / "topology_audit.json"),
            "prepare_report": str(paths.root / "prepare_report.json"),
            "analysis_input_summary": str(paths.root / "analysis_input_summary.json"),
            "colored_pdbs": [str(p) for p in colored_files],
        }

        hashes = {
            "raw_pdb_sha256": maybe_hash(paths.raw_pdb),
            "processed_pdb_sha256": maybe_hash(paths.processed_pdb),
            "topology_sha256": maybe_hash(paths.topology),
            "trajectory_sha256": maybe_hash(paths.trajectory),
            "energies_log_sha256": maybe_hash(paths.energies_log),
            "results_csv_sha256": maybe_hash(results_csv),
            "model_sha256": maybe_hash(self.cfg.model_path),
            "input_sanitization_sha256": maybe_hash(paths.root / "input_sanitization.json"),
            "prep_audit_sha256": maybe_hash(paths.root / "prep_audit.json"),
            "topology_audit_sha256": maybe_hash(paths.root / "topology_audit.json"),
            "prepare_report_sha256": maybe_hash(paths.root / "prepare_report.json"),
            "analysis_input_summary_sha256": maybe_hash(paths.root / "analysis_input_summary.json"),
            "colored_pdbs_sha256": {p.name: maybe_hash(p) for p in colored_files},
        }

        payload: dict[str, object] = {
            "run_id": self.run_id,
            "case": {
                "protein": case.protein,
                "seed": case.seed,
                "case_key": case.case_key,
            },
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "host": {
                "python": sys.version,
                "platform": platform.platform(),
            },
            "stages_executed": stages,
            "stage_results": [s.__dict__ for s in self.stage_results.get(case.case_key, [])],
            "config": config_to_dict(self.cfg),
            "files": files,
            "hashes": hashes,
            "results_summary": summarize_results_csv(results_csv),
            "commands": self.command_history.get(case.case_key, []),
        }

        payload["manifest_sha256"] = stable_payload_hash(payload)
        return payload
