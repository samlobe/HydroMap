from __future__ import annotations

import csv
import re
import warnings
from pathlib import Path
from typing import Any


_RESID_TOKEN_RE = re.compile(r"\bresid\s+(-?\d+)([A-Za-z]?)\b", re.IGNORECASE)
_ICODE_RE = re.compile(r"\bicode\s+([A-Za-z])\b", re.IGNORECASE)
_CHAIN_ID_RE = re.compile(r"\bchainID\s+([A-Za-z0-9_]+)\b", re.IGNORECASE)
_SEGID_RE = re.compile(r"\bsegid\s+([A-Za-z0-9_]+)\b", re.IGNORECASE)
_WHITESPACE_RE = re.compile(r"\s+")


CANONICAL_CHAIN_SENTINEL = "_"
_WARNED_CHAIN_SEGID_MISMATCHES: set[tuple[str, str]] = set()
_PLACEHOLDER_SEGIDS = {"SYSTEM"}


def canonical_residue_key(chain_id: str | None, resid: int, icode: str | None = None) -> str:
    chain = (chain_id or "").strip() or CANONICAL_CHAIN_SENTINEL
    ic = (icode or "").strip().upper() or "-"
    return f"chain:{chain}|resid:{int(resid)}|icode:{ic}"


def _unique_nonempty(values: Any) -> list[str]:
    seen: set[str] = set()
    unique: list[str] = []
    for value in values if values is not None else []:
        text = str(value).strip()
        if text and text not in seen:
            seen.add(text)
            unique.append(text)
    return unique


def residue_chain_token(residue: Any) -> tuple[str, list[str]]:
    """Return the chain token HydroMap should use for residue-scoped outputs.

    PDB ``chainID`` is the canonical chain identifier when present. ``segid``
    is reported as diagnostic metadata when it disagrees, but it is not used
    for residue-scoped output names or generated selections.
    """
    segid = str(getattr(residue, "segid", "") or "").strip()
    if segid.upper() in _PLACEHOLDER_SEGIDS:
        segid = ""
    chain_ids: list[str] = []
    try:
        chain_ids = _unique_nonempty(getattr(residue.atoms, "chainIDs", []))
    except Exception:
        chain_ids = []

    preferred = chain_ids[0] if len(chain_ids) == 1 else (segid or "")
    alternates: list[str] = []
    for token in [segid, *chain_ids]:
        if token and token != preferred and token not in alternates:
            alternates.append(token)

    if len(chain_ids) == 1 and segid and segid != chain_ids[0]:
        mismatch = (chain_ids[0], segid)
        if mismatch not in _WARNED_CHAIN_SEGID_MISMATCHES:
            _WARNED_CHAIN_SEGID_MISMATCHES.add(mismatch)
            warnings.warn(
                "Residue has conflicting PDB chainID and segid values; "
                f"using chainID '{chain_ids[0]}' and ignoring segid '{segid}' "
                "for HydroMap residue-scoped outputs.",
                RuntimeWarning,
                stacklevel=2,
            )

    return preferred, alternates


def residue_sort_token_chain(residue: Any) -> tuple[tuple[int, str], str, str, list[str]]:
    """Return sort key, residue token, preferred chain token, and alternates."""
    resid = int(getattr(residue, "resid", residue.resid))
    icode = (getattr(residue, "icode", "") or "").upper()
    sort_key = (resid, "" if icode == "" else icode)
    token = f"{resid}{icode}"
    chain, alternates = residue_chain_token(residue)
    return sort_key, token, chain, alternates


def canonical_group_key(selection: str) -> str:
    text = (selection or "").strip()
    m_res = _RESID_TOKEN_RE.search(text)
    if m_res is None:
        normalized = _WHITESPACE_RE.sub(" ", text)
        return f"selection:{normalized}"

    resid = int(m_res.group(1))
    icode = (m_res.group(2) or "").upper()
    if not icode:
        # Allow explicit insertion-code selectors like: "resid 11 and icode A".
        m_icode = _ICODE_RE.search(text)
        icode = (m_icode.group(1) if m_icode else "").upper()

    # Resolver precedence: chainID > segid > sentinel.
    m_chain = _CHAIN_ID_RE.search(text)
    if m_chain is not None:
        chain = m_chain.group(1)
    else:
        m_segid = _SEGID_RE.search(text)
        chain = m_segid.group(1) if m_segid else CANONICAL_CHAIN_SENTINEL

    return canonical_residue_key(chain, resid, icode)


def add_normalized_group_key_column(csv_path: Path) -> int:
    if not csv_path.exists():
        return 0

    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            return 0
        rows = list(reader)
        fields = list(reader.fieldnames)

    if "MDAnalysis_selection_strings" not in fields:
        return 0

    if "normalized_group_key" not in fields:
        fields.append("normalized_group_key")

    changed = 0
    for row in rows:
        selector = row.get("MDAnalysis_selection_strings", "")
        value = canonical_group_key(selector)
        if row.get("normalized_group_key") != value:
            changed += 1
        row["normalized_group_key"] = value

    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    return changed
