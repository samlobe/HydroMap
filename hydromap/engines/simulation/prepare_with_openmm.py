#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from openmm import XmlSerializer
from openmm import Vec3
from openmm.app import HBonds, Modeller, PDBFile, PME
from openmm.unit import molar, nanometer

from hydromap.forcefield import (
    add_adjacent_peptide_bonds,
    add_disulfide_bonds,
    add_explicit_rtp_impropers,
    assign_residue_templates,
    canonicalize_residue_names_for_hydrogens,
    compute_orthorhombic_box,
    ensure_c_terminal_oxygens,
    load_a99sbdisp_forcefield,
    prepare_peptide_topology,
    rename_solvent_for_forcefield_input,
    rename_solvent_for_output,
    rename_terminal_oxygen_atoms,
    restore_pdb_atom_names,
)
from hydromap.forcefield.openmm_utils import infer_variant_name


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare a HydroMap simulation system using OpenMM only.")
    parser.add_argument("input_pdb", help="Input protein PDB")
    parser.add_argument("--output_pdb", required=True, help="Prepared, solvated PDB output path")
    parser.add_argument("--output_system", required=True, help="Serialized OpenMM System XML output path")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for ion placement")
    parser.add_argument("--padding_nm", type=float, default=1.0, help="Solute-box padding in nm (default: 1.0)")
    parser.add_argument("--neutralize", choices=["true", "false"], default="true", help="Neutralize the solute charge with counterions (default: true)")
    parser.add_argument("--ionic_strength_molar", type=float, default=0.0, help="Additional monovalent salt concentration in molar after neutralization (default: 0.0)")
    parser.add_argument("--positive_ion", choices=["Na+", "K+", "Li+", "Rb+", "Cs+"], default="Na+", help="Positive ion species to use during solvation (default: Na+)")
    parser.add_argument(
        "--negative_ion",
        choices=["Cl-", "Br-", "F-", "I-"],
        default="Cl-",
        help="Negative ion species to use during solvation (default: Cl-)",
    )
    parser.add_argument("--histidine_mode", choices=["auto", "hid", "hie", "hip"], default="auto", help="Global histidine protonation mode (default: auto)")
    parser.add_argument("--histidine_override", action="append", default=[], help="Per-residue histidine override in the form CHAIN:RESID=HIP or RESID=HID. May be repeated.")
    parser.add_argument("--capping_mode", choices=["none", "termini", "breaks", "termini_and_breaks"], default="none", help="Optional peptide capping policy for true termini and/or internal breaks (default: none)")
    parser.add_argument("--audit_json", help="Optional JSON file capturing preparation audit details.")
    return parser.parse_args()


def _parse_histidine_overrides(raw_overrides: list[str]) -> dict[str, str]:
    overrides: dict[str, str] = {}
    for raw in raw_overrides:
        if "=" not in raw:
            raise ValueError(
                f"Invalid --histidine_override {raw!r}. Expected CHAIN:RESID=HIP or RESID=HID."
            )
        key, value = raw.split("=", 1)
        key = key.strip()
        value = value.strip().upper()
        if not key:
            raise ValueError(f"Invalid --histidine_override {raw!r}. Residue key is empty.")
        if value not in {"HID", "HIE", "HIP"}:
            raise ValueError(
                f"Invalid --histidine_override {raw!r}. Variant must be one of HID, HIE, HIP."
            )
        overrides[key.upper()] = value
    return overrides


def _histidine_override_for_residue(residue, overrides: dict[str, str]) -> str | None:
    chain_id = residue.chain.id.strip().upper()
    residue_id = residue.id.strip().upper()
    candidates = []
    if chain_id:
        candidates.append(f"{chain_id}:{residue_id}")
    candidates.append(residue_id)
    for candidate in candidates:
        if candidate in overrides:
            return overrides[candidate]
    return None


def choose_hydrogen_variants(
    topology,
    histidine_mode: str = "auto",
    histidine_overrides: dict[str, str] | None = None,
) -> list[str | None]:
    histidine_mode = histidine_mode.lower()
    histidine_overrides = histidine_overrides or {}
    variants: list[str | None] = []
    for residue in topology.residues():
        resname = residue.name.upper().strip()
        if resname == "HIS":
            override = _histidine_override_for_residue(residue, histidine_overrides)
            if override is not None:
                variants.append(override)
            elif histidine_mode == "auto":
                variants.append(infer_variant_name(residue))
            else:
                variants.append(histidine_mode.upper())
        elif resname == "ASP":
            variants.append(infer_variant_name(residue) if infer_variant_name(residue) == "ASH" else "ASP")
        elif resname == "GLU":
            variants.append(infer_variant_name(residue) if infer_variant_name(residue) == "GLH" else "GLU")
        elif resname == "LYS":
            variants.append(infer_variant_name(residue) if infer_variant_name(residue) == "LYN" else "LYS")
        elif resname == "CYS":
            variants.append(None)
        else:
            variants.append(None)
    return variants


def residue_templates_for_add_hydrogens(topology):
    templates = assign_residue_templates(topology)
    return {
        residue: template
        for residue, template in templates.items()
        if template not in {"ACE", "NME", "NHE", "NH2"}
        and template == residue.name.upper().strip()
    }


def rename_cap_atoms_for_openmm_hydrogens(topology) -> None:
    for residue in topology.residues():
        if residue.name.upper().strip() != "NME":
            continue
        for atom in residue.atoms():
            if atom.name == "CH3":
                atom.name = "C"


def restore_cap_atom_names_after_hydrogens(topology) -> None:
    for residue in topology.residues():
        resname = residue.name.upper().strip()
        if resname == "ACE":
            for atom in residue.atoms():
                if atom.name == "H1":
                    atom.name = "HH31"
                elif atom.name == "H2":
                    atom.name = "HH32"
                elif atom.name == "H3":
                    atom.name = "HH33"
        elif resname == "NME":
            for atom in residue.atoms():
                if atom.name == "C":
                    atom.name = "CH3"
                elif atom.name == "H1":
                    atom.name = "HH31"
                elif atom.name == "H2":
                    atom.name = "HH32"
                elif atom.name == "H3":
                    atom.name = "HH33"


def _jitter_positions(positions, amplitude_nm: float, seed: int):
    rng = random.Random(seed)
    jittered = []
    for pos in positions:
        dx = (rng.random() - 0.5) * 2.0 * amplitude_nm
        dy = (rng.random() - 0.5) * 2.0 * amplitude_nm
        dz = (rng.random() - 0.5) * 2.0 * amplitude_nm
        jittered.append(Vec3(pos[0] + dx, pos[1] + dy, pos[2] + dz))
    return jittered


def add_hydrogens_with_fallback(
    ff,
    topology,
    positions,
    residue_templates,
    seed: int,
    histidine_mode: str = "auto",
    histidine_overrides: dict[str, str] | None = None,
):
    variants = choose_hydrogen_variants(
        topology,
        histidine_mode=histidine_mode,
        histidine_overrides=histidine_overrides,
    )
    rename_cap_atoms_for_openmm_hydrogens(topology)
    modeller = Modeller(topology, positions)
    try:
        modeller.addHydrogens(ff, variants=variants, residueTemplates=residue_templates)
        restore_cap_atom_names_after_hydrogens(modeller.topology)
        return modeller
    except Exception as exc:
        if "Particle coordinate is NaN" not in str(exc):
            raise
        original_exc = exc

    seeds_to_try = [seed, 0, seed + 1, seed + 2]
    amplitudes_nm = [1e-3, 2e-3, 5e-3]
    for amplitude_nm in amplitudes_nm:
        for jitter_seed in seeds_to_try:
            jittered_positions = _jitter_positions(
                positions.value_in_unit(nanometer),
                amplitude_nm=amplitude_nm,
                seed=jitter_seed,
            )
            modeller = Modeller(topology, jittered_positions * nanometer)
            try:
                modeller.addHydrogens(ff, variants=variants, residueTemplates=residue_templates)
                restore_cap_atom_names_after_hydrogens(modeller.topology)
                return modeller
            except Exception as exc:
                if "Particle coordinate is NaN" not in str(exc):
                    raise
    raise original_exc


def main() -> None:
    args = parse_args()
    input_pdb = Path(args.input_pdb)
    output_pdb = Path(args.output_pdb)
    output_system = Path(args.output_system)
    neutralize = (args.neutralize == "true")
    if args.ionic_strength_molar < 0.0:
        raise ValueError("--ionic_strength_molar must be >= 0.")
    histidine_overrides = _parse_histidine_overrides(args.histidine_override)

    ff = load_a99sbdisp_forcefield()
    pdb = PDBFile(str(input_pdb))
    topology = pdb.topology
    restore_pdb_atom_names(topology, input_pdb)
    positions = [Vec3(*pos.value_in_unit(nanometer)) for pos in pdb.positions]

    add_disulfide_bonds(topology, positions)
    topology, positions, prep_audit = prepare_peptide_topology(topology, positions, capping_mode=args.capping_mode)
    add_adjacent_peptide_bonds(topology)
    topology, positions, _ = ensure_c_terminal_oxygens(topology, positions)
    rename_terminal_oxygen_atoms(topology)
    canonicalize_residue_names_for_hydrogens(topology)

    initial_templates = residue_templates_for_add_hydrogens(topology)
    random.seed(int(args.seed))
    modeller = add_hydrogens_with_fallback(
        ff,
        topology,
        positions * nanometer,
        initial_templates,
        seed=int(args.seed),
        histidine_mode=args.histidine_mode,
        histidine_overrides=histidine_overrides,
    )

    solute_templates = assign_residue_templates(modeller.topology)
    box_vectors = compute_orthorhombic_box(modeller.topology, modeller.positions, padding_nm=float(args.padding_nm))

    random.seed(int(args.seed))
    modeller.addSolvent(
        ff,
        model="tip4pew",
        boxVectors=box_vectors,
        positiveIon=args.positive_ion,
        negativeIon=args.negative_ion,
        ionicStrength=float(args.ionic_strength_molar) * molar,
        neutralize=neutralize,
        residueTemplates=solute_templates,
    )

    rename_solvent_for_forcefield_input(modeller.topology)
    final_templates = assign_residue_templates(modeller.topology)
    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * nanometer,
        constraints=HBonds,
        residueTemplates=final_templates,
    )
    add_explicit_rtp_impropers(system, modeller.topology, final_templates)
    rename_solvent_for_output(modeller.topology)

    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    output_system.parent.mkdir(parents=True, exist_ok=True)
    with output_pdb.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(modeller.topology, modeller.positions, handle)
    output_system.write_text(XmlSerializer.serialize(system), encoding="utf-8")
    if args.audit_json:
        audit_path = Path(args.audit_json)
        audit_path.parent.mkdir(parents=True, exist_ok=True)
        audit_path.write_text(json.dumps(prep_audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
