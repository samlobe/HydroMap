from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Iterable

import numpy as np
from openmm import Vec3
from openmm import XmlSerializer
from openmm import unit
from openmm.app import ForceField, GromacsTopFile, HBonds, PDBFile, PME, Topology, element as elem
from openmm.unit import nanometer

from .metadata import load_residue_template_metadata


PACKAGE_DIR = Path(__file__).resolve().parent
A99SBDISP_OPENMM_XML = PACKAGE_DIR / "a99SBdisp_openmm.xml"

CAP_NTERM_BLOCKERS = {"ACE"}
CAP_CTERM_BLOCKERS = {"NH2", "NHE", "NME"}
SPECIAL_WATER_NAMES = {"HOH", "WAT", "SOL"}
CANONICAL_RESNAME_MAP = {
    "HID": "HIS",
    "HIE": "HIS",
    "HIP": "HIS",
    "ASH": "ASP",
    "GLH": "GLU",
    "LYN": "LYS",
    "CYM": "CYS",
    "CYX": "CYS",
    "NHID": "HIS",
    "NHIE": "HIS",
    "NHIP": "HIS",
    "CHID": "HIS",
    "CHIE": "HIS",
    "CHIP": "HIS",
    "NASP": "ASP",
    "NGLU": "GLU",
    "NLYS": "LYS",
    "NCYM": "CYS",
    "NCYX": "CYS",
    "CASP": "ASP",
    "CGLU": "GLU",
    "CLYS": "LYS",
    "CCYX": "CYS",
}
DEGREES_TO_RADIANS = (1.0 * unit.degrees).value_in_unit(unit.radian)
PEPTIDE_BOND_MAX_DISTANCE_NM = 0.20
CAP_CLASH_DISTANCE_NM = 0.16
CAP_REFERENCE_ANCHORS_NM = np.array(
    [
        [1.550, 1.171, 1.199],
        [1.603, 1.307, 1.202],
        [1.759, 1.316, 1.233],
    ],
    dtype=float,
)
ACE_REFERENCE_ATOMS_NM: dict[str, np.ndarray] = {
    "C": np.array([1.461, 1.128, 1.114], dtype=float),
    "O": np.array([1.428, 1.191, 1.012], dtype=float),
    "CH3": np.array([1.426, 0.977, 1.124], dtype=float),
}
NME_REFERENCE_ATOMS_NM: dict[str, np.ndarray] = {
    "N": np.array([3.038, 2.346, 6.090], dtype=float) / 10.0,
    "CH3": np.array([3.145, 2.256, 6.130], dtype=float) / 10.0,
}
CAPPING_MODES = {"none", "termini", "breaks", "termini_and_breaks"}


def load_a99sbdisp_forcefield() -> ForceField:
    if not A99SBDISP_OPENMM_XML.exists():
        raise FileNotFoundError(
            f"Missing packaged OpenMM force field XML: {A99SBDISP_OPENMM_XML}."
        )
    return ForceField(str(A99SBDISP_OPENMM_XML))


def load_openmm_system(topology_path: str | Path, processed_pdb_path: str | Path):
    topology_path = Path(topology_path)
    processed_pdb_path = Path(processed_pdb_path)
    suffix = topology_path.suffix.lower()
    if suffix == ".xml":
        return XmlSerializer.deserialize(topology_path.read_text(encoding="utf-8"))
    if suffix != ".top":
        raise ValueError(
            f"Unsupported topology file {topology_path}. HydroMap analysis accepts .xml and .top files."
        )

    pdb = PDBFile(str(processed_pdb_path))
    periodic_box_vectors = pdb.topology.getPeriodicBoxVectors()
    top = GromacsTopFile(
        str(topology_path),
        periodicBoxVectors=periodic_box_vectors,
        includeDir=str(topology_path.parent),
    )
    return top.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * nanometer,
        constraints=HBonds,
    )


def serialize_a99sbdisp_system_for_processed_pdb(input_pdb: str | Path, output_system: str | Path) -> None:
    input_pdb = Path(input_pdb)
    output_system = Path(output_system)

    ff = load_a99sbdisp_forcefield()
    pdb = PDBFile(str(input_pdb))
    topology = pdb.topology
    restore_pdb_atom_names(topology, input_pdb)
    positions = [Vec3(*pos.value_in_unit(nanometer)) for pos in pdb.positions]

    add_disulfide_bonds(topology, positions)
    topology, positions, _ = prepare_peptide_topology(topology, positions, capping_mode="none")
    add_adjacent_peptide_bonds(topology)
    topology, positions, _ = ensure_c_terminal_oxygens(topology, positions)
    rename_terminal_oxygen_atoms(topology)
    canonicalize_residue_names_for_hydrogens(topology)
    rename_solvent_for_forcefield_input(topology)

    residue_templates = assign_residue_templates(topology)
    system = ff.createSystem(
        topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * nanometer,
        constraints=HBonds,
        residueTemplates=residue_templates,
    )
    add_explicit_rtp_impropers(system, topology, residue_templates)

    output_system.parent.mkdir(parents=True, exist_ok=True)
    output_system.write_text(XmlSerializer.serialize(system), encoding="utf-8")


def restore_pdb_atom_names(topology: Topology, pdb_path: str | Path) -> None:
    atom_names: list[str] = []
    previous_atom: tuple[str, str, str, str, str] | None = None
    pdb_path = Path(pdb_path)
    with pdb_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith(("ATOM  ", "HETATM")):
                # OpenMM collapses alternate locations onto a single atom entry.
                # Mirror that here so atom-name restoration stays aligned with
                # the parsed topology.
                key = (
                    line[21],
                    line[22:26],
                    line[26],
                    line[17:20],
                    line[12:16],
                )
                # Alternate locations for one atom are consecutive. Do not
                # deduplicate globally: residue identifiers legitimately wrap
                # in large explicit-solvent PDB files.
                if key == previous_atom:
                    continue
                previous_atom = key
                atom_names.append(line[12:16].strip())

    topology_atoms = list(topology.atoms())
    if len(atom_names) != len(topology_atoms):
        raise ValueError(
            f"PDB atom-name restore failed for {pdb_path}: file has {len(atom_names)} atoms, "
            f"but topology has {len(topology_atoms)} atoms."
        )

    for atom, name in zip(topology_atoms, atom_names):
        atom.name = name


@lru_cache(maxsize=1)
def known_residue_templates() -> dict[str, object]:
    return load_residue_template_metadata()


@lru_cache(maxsize=1)
def protein_base_residues() -> set[str]:
    residues = known_residue_templates()
    out = {
        name
        for name in residues
        if name not in SPECIAL_WATER_NAMES and name not in CAP_NTERM_BLOCKERS and name not in CAP_CTERM_BLOCKERS
        and not name.startswith(("N", "C"))
    }
    # Internal monatomic ions and water are not proteins.
    out -= {"NA", "CL", "K", "CA", "MG", "CS", "LI", "RB", "ZN", "IB+", "URE", "HOH", "HO4"}
    return out


@lru_cache(maxsize=1)
def canonical_base_residues() -> set[str]:
    return {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}


def canonicalize_residue_name(name: str) -> str:
    key = name.upper().strip()
    return CANONICAL_RESNAME_MAP.get(key, key)


def add_disulfide_bonds(topology: Topology, positions, cutoff_nm: float = 0.25) -> int:
    residues = list(topology.residues())
    sg_atoms = []
    for residue in residues:
        resname = residue.name.upper().strip()
        if canonicalize_residue_name(resname) != "CYS":
            continue
        for atom in residue.atoms():
            if atom.name == "SG":
                sg_atoms.append(atom)
                break

    added = 0
    for i, atom1 in enumerate(sg_atoms):
        pos1 = _pos_nm(positions[atom1.index])
        for atom2 in sg_atoms[i + 1 :]:
            if atom2.residue is atom1.residue:
                continue
            if any(
                (bond[0] is atom1 and bond[1] is atom2) or (bond[0] is atom2 and bond[1] is atom1)
                for bond in topology.bonds()
            ):
                continue
            pos2 = _pos_nm(positions[atom2.index])
            dx = pos1[0] - pos2[0]
            dy = pos1[1] - pos2[1]
            dz = pos1[2] - pos2[2]
            if dx * dx + dy * dy + dz * dz <= cutoff_nm * cutoff_nm:
                topology.addBond(atom1, atom2)
                added += 1
    return added


def add_explicit_rtp_impropers(system, topology: Topology, residue_templates: dict[object, str]) -> int:
    torsion_force = None
    for force in system.getForces():
        if force.__class__.__name__ == "PeriodicTorsionForce":
            torsion_force = force
            break
    if torsion_force is None:
        raise RuntimeError("Expected PeriodicTorsionForce in generated OpenMM System.")

    adjacency = _topology_adjacency(topology)
    existing: set[tuple[int, int, int, int, int, int, int]] = set()
    for idx in range(torsion_force.getNumTorsions()):
        a1, a2, a3, a4, periodicity, phase, k = torsion_force.getTorsionParameters(idx)
        existing.add(
            (
                a1,
                a2,
                a3,
                a4,
                int(periodicity),
                _phase_key(phase.value_in_unit(unit.radian)),
                _energy_key(k.value_in_unit(unit.kilojoule_per_mole)),
            )
        )

    added = 0
    residue_defs = known_residue_templates()
    for residue in topology.residues():
        template_name = residue_templates.get(residue)
        if template_name is None:
            continue
        residue_def = residue_defs.get(template_name)
        if residue_def is None or not residue_def.explicit_impropers:
            continue
        for improper in residue_def.explicit_impropers:
            atoms = []
            for token in (improper.atom1, improper.atom2, improper.atom3, improper.atom4):
                atom = _resolve_improper_atom(residue, token, adjacency)
                if atom is None:
                    atoms = []
                    break
                atoms.append(atom)
            if len(atoms) != 4:
                continue
            key = (
                atoms[0].index,
                atoms[1].index,
                atoms[2].index,
                atoms[3].index,
                improper.periodicity,
                _phase_key(improper.phase_degrees * DEGREES_TO_RADIANS),
                _energy_key(improper.k_kj_per_mol),
            )
            if key in existing:
                continue
            torsion_force.addTorsion(
                atoms[0].index,
                atoms[1].index,
                atoms[2].index,
                atoms[3].index,
                improper.periodicity,
                improper.phase_degrees * unit.degrees,
                improper.k_kj_per_mol * unit.kilojoule_per_mole,
            )
            existing.add(key)
            added += 1
    return added


def prepare_peptide_topology(
    topology: Topology,
    positions,
    capping_mode: str = "none",
) -> tuple[Topology, list[Vec3], dict[str, object]]:
    mode = str(capping_mode).strip().lower()
    if mode not in CAPPING_MODES:
        raise ValueError(f"Unsupported capping mode {capping_mode!r}. Expected one of {sorted(CAPPING_MODES)}.")

    segmented_topology, segmented_positions, segment_infos, segment_audit = _segment_topology_by_peptide_continuity(
        topology,
        positions,
    )
    capped_topology, capped_positions, cap_audit = _apply_capping_mode(
        segmented_topology,
        segmented_positions,
        segment_infos,
        capping_mode=mode,
    )
    audit = {
        "capping_mode": mode,
        **segment_audit,
        **cap_audit,
    }
    return capped_topology, capped_positions, audit


def _segment_topology_by_peptide_continuity(
    topology: Topology,
    positions,
) -> tuple[Topology, list[Vec3], list[dict[str, object]], dict[str, object]]:
    original_positions = [_pos_nm(pos) for pos in positions]
    rebuilt = Topology()
    rebuilt.setPeriodicBoxVectors(topology.getPeriodicBoxVectors())
    atom_map: dict[object, object] = {}
    rebuilt_positions: list[Vec3] = []
    segment_infos: list[dict[str, object]] = []
    break_boundaries: list[dict[str, str]] = []

    for original_chain_index, chain in enumerate(topology.chains()):
        current_new_chain = None
        current_segment: dict[str, object] | None = None
        previous_residue = None
        for residue in chain.residues():
            start_new_chain = current_new_chain is None
            if previous_residue is not None:
                prev_peptide = _is_peptide_residue(previous_residue)
                curr_peptide = _is_peptide_residue(residue)
                if prev_peptide or curr_peptide:
                    if not _should_link_adjacent_residues(previous_residue, residue, original_positions):
                        start_new_chain = True
                        if _is_protein_like_residue(previous_residue) and _is_protein_like_residue(residue):
                            break_boundaries.append(
                                {
                                    "left_residue": _format_residue_label(previous_residue),
                                    "right_residue": _format_residue_label(residue),
                                    "chain_id": chain.id.strip() or "",
                                }
                            )

            if start_new_chain:
                current_new_chain = rebuilt.addChain(chain.id)
                current_segment = {
                    "chain": current_new_chain,
                    "chain_id": chain.id.strip() or "",
                    "group_key": (chain.id.strip() or f"__orig_chain_{original_chain_index}"),
                    "original_chain_index": original_chain_index,
                    "contains_protein": False,
                    "n_boundary_type": None,
                    "c_boundary_type": None,
                }
                segment_infos.append(current_segment)

            assert current_new_chain is not None
            assert current_segment is not None
            new_residue = rebuilt.addResidue(residue.name, current_new_chain, residue.id, residue.insertionCode)
            if _is_protein_like_residue(residue):
                current_segment["contains_protein"] = True
            for atom in residue.atoms():
                new_atom = rebuilt.addAtom(atom.name, atom.element, new_residue, atom.id, atom.formalCharge)
                atom_map[atom] = new_atom
                rebuilt_positions.append(_pos_nm(original_positions[atom.index]))
            previous_residue = residue

    for atom1, atom2 in topology.bonds():
        if _should_skip_copied_bond(atom1, atom2, original_positions):
            continue
        rebuilt.addBond(atom_map[atom1], atom_map[atom2])

    protein_segment_indices_by_group: dict[str, list[int]] = {}
    for idx, info in enumerate(segment_infos):
        if info["contains_protein"]:
            protein_segment_indices_by_group.setdefault(str(info["group_key"]), []).append(idx)

    for group_indices in protein_segment_indices_by_group.values():
        for offset, segment_index in enumerate(group_indices):
            info = segment_infos[segment_index]
            info["n_boundary_type"] = "break" if offset > 0 else "terminus"
            info["c_boundary_type"] = "break" if offset < len(group_indices) - 1 else "terminus"

    return (
        rebuilt,
        rebuilt_positions,
        segment_infos,
        {
            "detected_peptide_break_count": len(break_boundaries),
            "detected_peptide_break_boundaries": break_boundaries,
            "protein_segment_count": sum(1 for info in segment_infos if info["contains_protein"]),
        },
    )


def _apply_capping_mode(
    topology: Topology,
    positions,
    segment_infos: list[dict[str, object]],
    capping_mode: str,
) -> tuple[Topology, list[Vec3], dict[str, object]]:
    original_positions = [_pos_nm(pos) for pos in positions]
    rebuilt = Topology()
    rebuilt.setPeriodicBoxVectors(topology.getPeriodicBoxVectors())
    atom_map: dict[object, object] = {}
    rebuilt_positions: list[Vec3] = []
    added_ace_caps: list[dict[str, str]] = []
    added_nme_caps: list[dict[str, str]] = []
    skipped_caps: list[dict[str, str]] = []
    warnings: list[str] = []
    preserved_input_ace_caps = 0
    preserved_input_c_caps = 0

    chains = list(topology.chains())
    chain_to_info = {id(info["chain"]): info for info in segment_infos}

    for chain in chains:
        info = chain_to_info.get(id(chain), {})
        residues = list(chain.residues())
        new_chain = rebuilt.addChain(chain.id)
        protein_residues = [residue for residue in residues if _is_protein_like_residue(residue)]
        first_protein = protein_residues[0] if protein_residues else None
        last_protein = protein_residues[-1] if protein_residues else None
        existing_n_cap = residues and residues[0].name.upper().strip() in CAP_NTERM_BLOCKERS
        existing_c_cap = residues and residues[-1].name.upper().strip() in CAP_CTERM_BLOCKERS
        if existing_n_cap:
            preserved_input_ace_caps += 1
        if existing_c_cap:
            preserved_input_c_caps += 1

        insert_ace = (
            first_protein is not None
            and not existing_n_cap
            and _boundary_matches_mode(str(info.get("n_boundary_type") or ""), capping_mode)
        )
        insert_nme = (
            last_protein is not None
            and not existing_c_cap
            and _boundary_matches_mode(str(info.get("c_boundary_type") or ""), capping_mode)
        )

        if insert_ace and first_protein is not None:
            ace_record = _build_cap_record(
                cap_name="ACE",
                anchor_residue=first_protein,
                topology=topology,
                positions=original_positions,
            )
            if ace_record is None:
                skipped_caps.append(
                    {
                        "cap": "ACE",
                        "site": _format_residue_label(first_protein),
                        "boundary_type": str(info.get("n_boundary_type") or ""),
                        "reason": "missing_anchor_atoms_or_clash",
                    }
                )
                warnings.append(
                    f"Skipped ACE capping at {_format_residue_label(first_protein)} due to missing anchor atoms or clashes."
                )
            else:
                _add_cap_residue(rebuilt, new_chain, anchor_residue=first_protein, cap_record=ace_record)
                rebuilt_positions.extend(ace_record["positions"])
                added_ace_caps.append(
                    {
                        "site": _format_residue_label(first_protein),
                        "boundary_type": str(info.get("n_boundary_type") or ""),
                    }
                )

        for residue in residues:
            new_residue = rebuilt.addResidue(residue.name, new_chain, residue.id, residue.insertionCode)
            for atom in residue.atoms():
                new_atom = rebuilt.addAtom(atom.name, atom.element, new_residue, atom.id, atom.formalCharge)
                atom_map[atom] = new_atom
                rebuilt_positions.append(_pos_nm(original_positions[atom.index]))

        if insert_nme and last_protein is not None:
            nme_record = _build_cap_record(
                cap_name="NME",
                anchor_residue=last_protein,
                topology=topology,
                positions=original_positions,
            )
            if nme_record is None:
                skipped_caps.append(
                    {
                        "cap": "NME",
                        "site": _format_residue_label(last_protein),
                        "boundary_type": str(info.get("c_boundary_type") or ""),
                        "reason": "missing_anchor_atoms_or_clash",
                    }
                )
                warnings.append(
                    f"Skipped NME capping at {_format_residue_label(last_protein)} due to missing anchor atoms or clashes."
                )
            else:
                _add_cap_residue(rebuilt, new_chain, anchor_residue=last_protein, cap_record=nme_record)
                rebuilt_positions.extend(nme_record["positions"])
                added_nme_caps.append(
                    {
                        "site": _format_residue_label(last_protein),
                        "boundary_type": str(info.get("c_boundary_type") or ""),
                    }
                )

    for atom1, atom2 in topology.bonds():
        rebuilt.addBond(atom_map[atom1], atom_map[atom2])

    return (
        rebuilt,
        rebuilt_positions,
        {
            "added_ace_caps": added_ace_caps,
            "added_ace_cap_count": len(added_ace_caps),
            "added_nme_caps": added_nme_caps,
            "added_nme_cap_count": len(added_nme_caps),
            "preserved_input_nterm_cap_count": preserved_input_ace_caps,
            "preserved_input_cterm_cap_count": preserved_input_c_caps,
            "skipped_caps": skipped_caps,
            "warnings": warnings,
        },
    )


def _build_cap_record(cap_name: str, anchor_residue, topology: Topology, positions: list[Vec3]) -> dict[str, object] | None:
    anchor_positions = _anchor_positions(anchor_residue, positions)
    if anchor_positions is None:
        return None

    if cap_name == "ACE":
        atom_positions = _transformed_cap_positions(anchor_positions, ACE_REFERENCE_ATOMS_NM)
        atom_specs = [("CH3", elem.carbon), ("C", elem.carbon), ("O", elem.oxygen)]
    elif cap_name == "NME":
        atom_positions = _transformed_cap_positions(anchor_positions, NME_REFERENCE_ATOMS_NM)
        atom_specs = [("N", elem.nitrogen), ("CH3", elem.carbon)]
    else:
        raise ValueError(f"Unsupported cap name {cap_name!r}")

    if _cap_positions_clash(atom_positions, anchor_residue, topology, positions):
        return None

    return {
        "cap_name": cap_name,
        "atom_specs": atom_specs,
        "atom_positions": atom_positions,
        "positions": [Vec3(*atom_positions[name]) for name, _ in atom_specs],
    }


def _add_cap_residue(rebuilt: Topology, chain, anchor_residue, cap_record: dict[str, object]) -> None:
    cap_name = str(cap_record["cap_name"])
    new_residue = rebuilt.addResidue(cap_name, chain, anchor_residue.id, anchor_residue.insertionCode)
    created_atoms: dict[str, object] = {}
    for atom_name, element in cap_record["atom_specs"]:
        created_atoms[atom_name] = rebuilt.addAtom(atom_name, element, new_residue)
    if cap_name == "ACE":
        rebuilt.addBond(created_atoms["C"], created_atoms["O"])
        rebuilt.addBond(created_atoms["C"], created_atoms["CH3"])
    elif cap_name == "NME":
        rebuilt.addBond(created_atoms["N"], created_atoms["CH3"])


def _anchor_positions(residue, positions: list[Vec3]) -> np.ndarray | None:
    anchor_atoms = [_find_atom(residue, name) for name in ("N", "CA", "C")]
    if any(atom is None for atom in anchor_atoms):
        return None
    return np.vstack([np.asarray(_pos_nm(positions[atom.index]), dtype=float) for atom in anchor_atoms])


def _transformed_cap_positions(anchor_positions: np.ndarray, reference_atoms: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    rotation, translation = _rigid_transform(CAP_REFERENCE_ANCHORS_NM, anchor_positions)
    return {name: (rotation @ coords) + translation for name, coords in reference_atoms.items()}


def _rigid_transform(reference_points: np.ndarray, target_points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    ref_centroid = reference_points.mean(axis=0)
    target_centroid = target_points.mean(axis=0)
    ref_centered = reference_points - ref_centroid
    target_centered = target_points - target_centroid
    covariance = ref_centered.T @ target_centered
    u, _, vt = np.linalg.svd(covariance)
    rotation = vt.T @ u.T
    if np.linalg.det(rotation) < 0:
        vt[-1, :] *= -1
        rotation = vt.T @ u.T
    translation = target_centroid - (rotation @ ref_centroid)
    return rotation, translation


def _cap_positions_clash(
    atom_positions: dict[str, np.ndarray],
    anchor_residue,
    topology: Topology,
    positions: list[Vec3],
) -> bool:
    excluded_indices = {atom.index for atom in anchor_residue.atoms()}
    threshold_sq = CAP_CLASH_DISTANCE_NM * CAP_CLASH_DISTANCE_NM
    for atom in topology.atoms():
        if atom.index in excluded_indices:
            continue
        if atom.element is None or atom.element == elem.hydrogen:
            continue
        existing = np.asarray(_pos_nm(positions[atom.index]), dtype=float)
        for coords in atom_positions.values():
            delta = coords - existing
            if float(delta @ delta) < threshold_sq:
                return True
    return False


def _boundary_matches_mode(boundary_type: str, capping_mode: str) -> bool:
    if capping_mode == "none" or not boundary_type:
        return False
    if capping_mode == "termini":
        return boundary_type == "terminus"
    if capping_mode == "breaks":
        return boundary_type == "break"
    if capping_mode == "termini_and_breaks":
        return boundary_type in {"terminus", "break"}
    return False


def add_adjacent_peptide_bonds(topology: Topology) -> int:
    added = 0
    for chain in topology.chains():
        residues = list(chain.residues())
        for current, nxt in zip(residues, residues[1:]):
            if not _can_form_peptide_bond(current, nxt):
                continue
            c_atom = _find_atom(current, "C")
            n_atom = _find_atom(nxt, "N")
            if c_atom is None or n_atom is None:
                continue
            if any(
                (bond[0] is c_atom and bond[1] is n_atom) or (bond[0] is n_atom and bond[1] is c_atom)
                for bond in topology.bonds()
            ):
                continue
            topology.addBond(c_atom, n_atom)
            added += 1
    return added


def ensure_c_terminal_oxygens(topology: Topology, positions) -> tuple[Topology, list[Vec3], int]:
    original_positions = [_pos_nm(pos) for pos in positions]
    residues_needing_oxt: set[int] = set()
    new_oxt_positions: dict[int, Vec3] = {}
    for chain in topology.chains():
        chain_residues = list(chain.residues())
        for idx, residue in enumerate(chain_residues):
            if not _is_protein_like_residue(residue):
                continue
            next_residue = chain_residues[idx + 1] if idx + 1 < len(chain_residues) else None
            if next_residue is not None and next_residue.name.upper().strip() in CAP_CTERM_BLOCKERS:
                continue
            if next_residue is not None and _is_protein_like_residue(next_residue):
                continue
            oxygen_names = {atom.name for atom in residue.atoms() if atom.element == elem.oxygen}
            if "OC2" in oxygen_names or "OXT" in oxygen_names:
                continue
            c_atom = _find_atom(residue, "C")
            o_atom = _find_first_existing_atom(residue, ["O", "OC1"])
            ca_atom = _find_atom(residue, "CA")
            if c_atom is None or o_atom is None or ca_atom is None:
                continue
            c_pos = _pos_nm(original_positions[c_atom.index])
            o_pos = _pos_nm(original_positions[o_atom.index])
            ca_pos = _pos_nm(original_positions[ca_atom.index])
            vec_co = _normalize(c_pos - o_pos)
            vec_ca = _normalize(c_pos - ca_pos)
            direction = _normalize(vec_co + vec_ca)
            new_oxt_positions[residue.index] = Vec3(
                c_pos[0] + direction[0] * 0.123,
                c_pos[1] + direction[1] * 0.123,
                c_pos[2] + direction[2] * 0.123,
            )
            residues_needing_oxt.add(residue.index)

    if not residues_needing_oxt:
        return topology, original_positions, 0

    rebuilt = Topology()
    rebuilt.setPeriodicBoxVectors(topology.getPeriodicBoxVectors())
    atom_map: dict[object, object] = {}
    rebuilt_positions: list[Vec3] = []

    for chain in topology.chains():
        new_chain = rebuilt.addChain(chain.id)
        for residue in chain.residues():
            new_residue = rebuilt.addResidue(residue.name, new_chain, residue.id, residue.insertionCode)
            for atom in residue.atoms():
                new_atom = rebuilt.addAtom(atom.name, atom.element, new_residue, atom.id, atom.formalCharge)
                atom_map[atom] = new_atom
                rebuilt_positions.append(_pos_nm(original_positions[atom.index]))
            if residue.index in residues_needing_oxt:
                new_atom = rebuilt.addAtom("OXT", elem.oxygen, new_residue)
                rebuilt_positions.append(new_oxt_positions[residue.index])
                c_atom = _find_atom(residue, "C")
                if c_atom is not None:
                    atom_map[(residue.index, "OXT")] = new_atom
                    rebuilt.addBond(atom_map[c_atom], new_atom)

    for atom1, atom2 in topology.bonds():
        rebuilt.addBond(atom_map[atom1], atom_map[atom2])

    return rebuilt, rebuilt_positions, len(residues_needing_oxt)


def rename_terminal_oxygen_atoms(topology: Topology) -> None:
    for chain in topology.chains():
        chain_residues = list(chain.residues())
        for idx, residue in enumerate(chain_residues):
            if not _is_protein_like_residue(residue):
                continue
            next_residue = chain_residues[idx + 1] if idx + 1 < len(chain_residues) else None
            if next_residue is not None and next_residue.name.upper().strip() in CAP_CTERM_BLOCKERS:
                continue
            if next_residue is not None and _is_protein_like_residue(next_residue):
                continue
            o_atom = _find_first_existing_atom(residue, ["O", "OC1"])
            oxt_atom = _find_first_existing_atom(residue, ["OXT", "OC2"])
            if o_atom is not None:
                o_atom.name = "OC1"
            if oxt_atom is not None:
                oxt_atom.name = "OC2"


def canonicalize_residue_names_for_hydrogens(topology: Topology) -> None:
    for residue in topology.residues():
        resname = residue.name.upper().strip()
        residue.name = canonicalize_residue_name(resname)


def infer_variant_name(residue) -> str:
    atom_names = {atom.name for atom in residue.atoms()}
    base = residue.name.upper().strip()
    has_explicit_hydrogens = any(name.startswith("H") for name in atom_names)
    if base == "HIS":
        has_hd1 = "HD1" in atom_names
        has_he2 = "HE2" in atom_names
        if has_hd1 and has_he2:
            return "HIP"
        if has_hd1:
            return "HID"
        return "HIE"
    if base == "ASP" and "HD2" in atom_names:
        return "ASH"
    if base == "GLU" and "HE2" in atom_names:
        return "GLH"
    if base == "LYS" and has_explicit_hydrogens and not {"HZ1", "HZ2", "HZ3"}.issubset(atom_names):
        return "LYN"
    if base == "CYS":
        sg = _find_atom(residue, "SG")
        if sg is not None:
            for bond in residue.external_bonds():
                if bond[0] is sg or bond[1] is sg:
                    return "CYX"
        if has_explicit_hydrogens and "HG" not in atom_names:
            return "CYM"
    return base


def assign_residue_templates(topology: Topology) -> dict[object, str]:
    templates: dict[object, str] = {}
    for chain in topology.chains():
        chain_residues = list(chain.residues())
        for idx, residue in enumerate(chain_residues):
            resname = residue.name.upper().strip()
            if resname in SPECIAL_WATER_NAMES:
                templates[residue] = "SOL" if resname == "SOL" else "HOH"
                continue
            if resname in {"NA", "NA+", "NA+1", "NA+1.0"}:
                templates[residue] = "NA" if resname == "NA" else "Na+"
                continue
            if resname in {"CL", "CL-", "CL-1", "CL-1.0"}:
                templates[residue] = "CL" if resname == "CL" else "Cl-"
                continue
            if resname in {"BR", "BR-", "BR-1", "BR-1.0"}:
                templates[residue] = "BR" if resname == "BR" else "Br-"
                continue
            if resname in {"F", "F-", "F-1", "F-1.0"}:
                templates[residue] = "F" if resname == "F" else "F-"
                continue
            if resname in {"I", "I-", "I-1", "I-1.0"}:
                templates[residue] = "I" if resname == "I" else "I-"
                continue
            if resname in {"K", "K+"}:
                templates[residue] = "K+" if resname.endswith("+") else "K"
                continue
            if resname in {"ACE", "NH2", "NHE", "NME"}:
                templates[residue] = resname
                continue
            if not _is_protein_like_residue(residue):
                continue

            variant = infer_variant_name(residue)
            prev_residue = chain_residues[idx - 1] if idx > 0 else None
            next_residue = chain_residues[idx + 1] if idx + 1 < len(chain_residues) else None
            has_prev_link = False
            if prev_residue is not None:
                prev_name = prev_residue.name.upper().strip()
                has_prev_link = prev_name in CAP_NTERM_BLOCKERS or _is_protein_like_residue(prev_residue)
            has_next_link = False
            if next_residue is not None:
                next_name = next_residue.name.upper().strip()
                has_next_link = next_name in CAP_CTERM_BLOCKERS or _is_protein_like_residue(next_residue)
            is_nterm = not has_prev_link
            is_cterm = not has_next_link
            if is_nterm and is_cterm:
                raise ValueError(
                    f"Standalone protein residue {residue.name} {residue.id} is unsupported by a99SB-disp. "
                    "Use a capped residue or a longer peptide/protein."
                )

            template = variant
            if is_nterm:
                template = f"N{variant}"
            elif is_cterm:
                template = f"C{variant}"
            if template not in known_residue_templates():
                raise ValueError(f"No HydroMap/OpenMM template available for residue {residue.name} -> {template}")
            templates[residue] = template
    return templates


def rename_solvent_for_output(topology: Topology) -> None:
    for residue in topology.residues():
        resname = residue.name.upper().strip()
        if resname == "HOH":
            residue.name = "SOL"
            for atom in residue.atoms():
                if atom.name == "O":
                    atom.name = "OW"
                elif atom.name == "H1":
                    atom.name = "HW1"
                elif atom.name == "H2":
                    atom.name = "HW2"
                elif atom.name == "M":
                    atom.name = "MW"
        elif resname == "NA+":
            residue.name = "NA"
            atom = next(residue.atoms())
            atom.name = "NA"
        elif resname == "NA":
            atom = next(residue.atoms())
            atom.name = "NA"
        elif resname == "CL-":
            residue.name = "CL"
            atom = next(residue.atoms())
            atom.name = "CL"
        elif resname == "CL":
            atom = next(residue.atoms())
            atom.name = "CL"
        elif resname == "BR-":
            residue.name = "BR"
            atom = next(residue.atoms())
            atom.name = "BR"
        elif resname == "BR":
            atom = next(residue.atoms())
            atom.name = "BR"
        elif resname == "F-":
            residue.name = "F"
            atom = next(residue.atoms())
            atom.name = "F"
        elif resname == "F":
            atom = next(residue.atoms())
            atom.name = "F"
        elif resname == "I-":
            residue.name = "I"
            atom = next(residue.atoms())
            atom.name = "I"
        elif resname == "I":
            atom = next(residue.atoms())
            atom.name = "I"


def rename_solvent_for_forcefield_input(topology: Topology) -> None:
    for residue in topology.residues():
        resname = residue.name.upper().strip()
        if resname == "SOL":
            residue.name = "HOH"
            for atom in residue.atoms():
                if atom.name == "OW":
                    atom.name = "O"
                elif atom.name == "HW1":
                    atom.name = "H1"
                elif atom.name == "HW2":
                    atom.name = "H2"
                elif atom.name == "MW":
                    atom.name = "M"
        elif resname == "BR":
            residue.name = "BR-"
        elif resname == "F":
            residue.name = "F-"
        elif resname == "I":
            residue.name = "I-"


def compute_orthorhombic_box(topology: Topology, positions, padding_nm: float = 1.0) -> tuple[Vec3, Vec3, Vec3]:
    coords = [_pos_nm(pos) for pos in positions]
    min_xyz = [min(p[i] for p in coords) for i in range(3)]
    max_xyz = [max(p[i] for p in coords) for i in range(3)]
    sizes = [max_xyz[i] - min_xyz[i] + 2.0 * padding_nm for i in range(3)]
    return (
        Vec3(sizes[0], 0.0, 0.0),
        Vec3(0.0, sizes[1], 0.0),
        Vec3(0.0, 0.0, sizes[2]),
    )


def _is_cap_residue(residue) -> bool:
    return residue.name.upper().strip() in CAP_NTERM_BLOCKERS | CAP_CTERM_BLOCKERS


def _is_protein_like_residue(residue) -> bool:
    resname = residue.name.upper().strip()
    if resname in SPECIAL_WATER_NAMES or resname in CAP_NTERM_BLOCKERS or resname in CAP_CTERM_BLOCKERS:
        return False
    if resname in {"NA", "NA+", "CL", "CL-", "BR", "BR-", "F", "F-", "I", "I-", "K", "K+", "CA", "MG", "CS", "RB", "LI", "ZN"}:
        return False
    return canonicalize_residue_name(resname) in canonical_base_residues() or resname in protein_base_residues()


def _is_peptide_residue(residue) -> bool:
    return _is_protein_like_residue(residue) or _is_cap_residue(residue)


def _can_form_peptide_bond(current, nxt) -> bool:
    current_name = current.name.upper().strip()
    next_name = nxt.name.upper().strip()
    current_ok = _is_protein_like_residue(current) or current_name in CAP_NTERM_BLOCKERS
    next_ok = _is_protein_like_residue(nxt) or next_name in CAP_CTERM_BLOCKERS
    return current_ok and next_ok


def _should_link_adjacent_residues(current, nxt, positions: list[Vec3]) -> bool:
    if not _can_form_peptide_bond(current, nxt):
        return False
    c_atom = _find_atom(current, "C")
    n_atom = _find_atom(nxt, "N")
    if c_atom is None or n_atom is None:
        return False
    c_pos = np.asarray(_pos_nm(positions[c_atom.index]), dtype=float)
    n_pos = np.asarray(_pos_nm(positions[n_atom.index]), dtype=float)
    delta = c_pos - n_pos
    return float(delta @ delta) <= PEPTIDE_BOND_MAX_DISTANCE_NM * PEPTIDE_BOND_MAX_DISTANCE_NM


def _should_skip_copied_bond(atom1, atom2, positions: list[Vec3]) -> bool:
    if atom1.residue is atom2.residue:
        return False
    if {atom1.name, atom2.name} != {"C", "N"}:
        return False
    res1 = atom1.residue
    res2 = atom2.residue
    if not (_is_peptide_residue(res1) and _is_peptide_residue(res2)):
        return False
    current = res1 if atom1.name == "C" else res2
    nxt = res2 if atom2.name == "N" else res1
    return not _should_link_adjacent_residues(current, nxt, positions)


def _format_residue_label(residue) -> str:
    chain = residue.chain.id.strip() or "?"
    icode = getattr(residue, "insertionCode", "") or ""
    label = f"{residue.name.upper().strip()} {chain}:{residue.id}"
    if icode:
        label += icode
    return label


def _find_atom(residue, name: str):
    for atom in residue.atoms():
        if atom.name == name:
            return atom
    return None


def _topology_adjacency(topology: Topology) -> dict[int, list[object]]:
    adjacency: dict[int, list[object]] = {}
    for atom1, atom2 in topology.bonds():
        adjacency.setdefault(atom1.index, []).append(atom2)
        adjacency.setdefault(atom2.index, []).append(atom1)
    return adjacency


def _resolve_improper_atom(residue, token: str, adjacency: dict[int, list[object]]):
    if not token.startswith(("+", "-")):
        return _find_atom(residue, token)
    target_name = token[1:]
    forward = token.startswith("+")
    candidates = []
    seen = set()
    for atom in residue.atoms():
        for other in adjacency.get(atom.index, []):
            if other.residue is residue or other.name != target_name:
                continue
            if other.index in seen:
                continue
            if forward and other.residue.index <= residue.index:
                continue
            if not forward and other.residue.index >= residue.index:
                continue
            candidates.append(other)
            seen.add(other.index)
    if len(candidates) == 1:
        return candidates[0]
    return None


def _phase_key(value: float) -> int:
    return int(round(value * 1_000_000))


def _energy_key(value: float) -> int:
    return int(round(value * 1_000_000))


def _find_first_existing_atom(residue, names: Iterable[str]):
    name_set = set(names)
    for atom in residue.atoms():
        if atom.name in name_set:
            return atom
    return None


def _normalize(vec) -> Vec3:
    if isinstance(vec, Vec3):
        x, y, z = vec[0], vec[1], vec[2]
    else:
        x, y, z = vec
    norm = (x * x + y * y + z * z) ** 0.5
    if norm == 0.0:
        return Vec3(1.0, 0.0, 0.0)
    return Vec3(x / norm, y / norm, z / norm)


def _pos_nm(pos) -> Vec3:
    if hasattr(pos, "value_in_unit"):
        val = pos.value_in_unit(nanometer)
        return Vec3(val[0], val[1], val[2])
    return Vec3(pos[0], pos[1], pos[2])
