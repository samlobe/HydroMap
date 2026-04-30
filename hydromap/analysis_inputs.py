from __future__ import annotations

from pathlib import Path

WATER_RESIDUE_NAMES = {
    "SOL",
    "HOH",
    "WAT",
    "TIP3",
    "TIP3P",
    "TIP4",
    "TIP4P",
    "TIP4PEW",
    "SPC",
    "SPCE",
}
WATER_OXYGEN_NAMES = {"O", "OW", "OH2"}
WATER_HYDROGEN_MAP = {
    "H1": "HW1",
    "H2": "HW2",
    "HW1": "HW1",
    "HW2": "HW2",
}
WATER_M_SITE_NAMES = {"M", "MW", "EPW"}
ION_RESNAME_MAP = {
    "NA+": "NA",
    "NA": "NA",
    "CL-": "CL",
    "CL": "CL",
    "BR-": "BR",
    "BR": "BR",
    "F-": "F",
    "F": "F",
    "I-": "I",
    "I": "I",
    "K+": "K",
    "K": "K",
    "CA": "CA",
    "MG": "MG",
    "CS": "CS",
    "RB": "RB",
    "LI": "LI",
    "ZN": "ZN",
}
STRIP_FROM_RAW_RESNAMES = set(ION_RESNAME_MAP.values()) | WATER_RESIDUE_NAMES | {"SOL"}
GENERIC_WATER_SELECTION = "resname SOL HOH WAT TIP3 TIP3P TIP4 TIP4P TIP4PEW SPC SPCE"
GENERIC_WATER_OXYGEN_SELECTION = f"({GENERIC_WATER_SELECTION}) and name O OW OH2"
NORMALIZED_WATER_SELECTION = "resname SOL"
NORMALIZED_WATER_OXYGEN_SELECTION = "resname SOL and name OW"


def is_water_residue_name(name: str) -> bool:
    return name.upper().strip() in WATER_RESIDUE_NAMES or name.upper().strip() == "SOL"


def normalize_external_processed_pdb(src_pdb: str | Path, dst_processed_pdb: str | Path, dst_raw_pdb: str | Path) -> dict[str, int | bool | str]:
    src_pdb = Path(src_pdb)
    dst_processed_pdb = Path(dst_processed_pdb)
    dst_raw_pdb = Path(dst_raw_pdb)

    normalized_lines: list[str] = []
    raw_lines: list[str] = []
    solvent_residue_keys: set[tuple[str, str, str, str]] = set()
    water_oxygen_count = 0
    protein_atom_count = 0

    with src_pdb.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            record = line[:6]
            if record.startswith(("ATOM  ", "HETATM")) and len(line) >= 26:
                atom_name = line[12:16].strip().upper()
                resname = line[17:20].strip().upper()
                chain_id = line[21:22]
                resid = line[22:26]
                icode = line[26:27] if len(line) > 26 else " "

                out_resname = resname
                out_atom_name = atom_name
                is_water = is_water_residue_name(resname)
                if is_water:
                    out_resname = "SOL"
                    if atom_name in WATER_OXYGEN_NAMES:
                        out_atom_name = "OW"
                        solvent_residue_keys.add((chain_id, resid, icode, out_resname))
                        water_oxygen_count += 1
                    elif atom_name in WATER_HYDROGEN_MAP:
                        out_atom_name = WATER_HYDROGEN_MAP[atom_name]
                    elif atom_name in WATER_M_SITE_NAMES:
                        out_atom_name = "MW"
                elif resname in ION_RESNAME_MAP:
                    out_resname = ION_RESNAME_MAP[resname]
                    out_atom_name = out_resname
                else:
                    protein_atom_count += 1

                rewritten = _rewrite_pdb_fields(line, atom_name=out_atom_name, resname=out_resname)
                normalized_lines.append(rewritten)
                if out_resname not in STRIP_FROM_RAW_RESNAMES:
                    raw_lines.append(rewritten)
                continue

            normalized_lines.append(line)
            if line.startswith("CRYST1"):
                raw_lines.append(line)

    if water_oxygen_count <= 0 or not solvent_residue_keys:
        raise ValueError(
            "HydroMap water analysis requires an explicit-solvent PDB/trajectory. "
            "No water molecules were detected in the provided processed PDB."
        )
    if protein_atom_count <= 0:
        raise ValueError("No non-solvent atoms were found in the provided processed PDB.")

    if raw_lines and not any(line.startswith("TER") for line in raw_lines if line.startswith(("ATOM  ", "HETATM"))):
        raw_lines.append("TER\n")
    raw_lines.append("END\n")
    if not normalized_lines or not normalized_lines[-1].rstrip().endswith("END"):
        normalized_lines.append("END\n")

    dst_processed_pdb.parent.mkdir(parents=True, exist_ok=True)
    dst_raw_pdb.parent.mkdir(parents=True, exist_ok=True)
    dst_processed_pdb.write_text("".join(normalized_lines), encoding="utf-8")
    dst_raw_pdb.write_text("".join(raw_lines), encoding="utf-8")

    return {
        "source_processed_pdb": str(src_pdb),
        "normalized_processed_pdb": str(dst_processed_pdb),
        "derived_raw_pdb": str(dst_raw_pdb),
        "has_explicit_solvent": True,
        "water_residue_count": len(solvent_residue_keys),
        "water_oxygen_count": water_oxygen_count,
        "protein_atom_count": protein_atom_count,
    }


def _rewrite_pdb_fields(line: str, atom_name: str, resname: str) -> str:
    atom_field = atom_name.rjust(4)[:4]
    resname_field = resname.rjust(3)[:3]
    suffix = line[20:] if len(line) > 20 else "\n"
    return f"{line[:12]}{atom_field}{line[16:17]}{resname_field}{suffix}"
