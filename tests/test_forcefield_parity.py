from __future__ import annotations

from pathlib import Path
import hashlib
import json
import subprocess
import sys

import pytest
import numpy as np

from openmm import Context, Platform, VerletIntegrator, XmlSerializer, unit
from openmm.app import HBonds, PDBFile, PME

from hydromap.forcefield import (
    add_adjacent_peptide_bonds,
    add_disulfide_bonds,
    add_explicit_rtp_impropers,
    assign_residue_templates,
    load_a99sbdisp_forcefield,
    rename_solvent_for_forcefield_input,
    restore_pdb_atom_names,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURE_ROOT = REPO_ROOT / "tests" / "fixtures" / "forcefield_parity"
PARITY_CASES = ("F", "SLS1_monomer", "alpha")


def _value(x, unit_obj=None) -> float:
    return x.value_in_unit(unit_obj) if hasattr(x, "value_in_unit") else float(x)


def _force(system, cls_name: str):
    for force in system.getForces():
        if force.__class__.__name__ == cls_name:
            return force
    raise RuntimeError(f"Missing force: {cls_name}")


def _digest_rows(rows: list[tuple]) -> str:
    payload = json.dumps(sorted(rows), separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _optional_force(system, cls_name: str):
    for force in system.getForces():
        if force.__class__.__name__ == cls_name:
            return force
    return None


def _nonbonded_particle_arrays(system):
    nb = _force(system, "NonbondedForce")
    charges = np.empty(system.getNumParticles(), dtype=np.float64)
    sigma = np.empty(system.getNumParticles(), dtype=np.float64)
    epsilon = np.empty(system.getNumParticles(), dtype=np.float64)
    for idx in range(system.getNumParticles()):
        charge, s, e = nb.getParticleParameters(idx)
        charges[idx] = _value(charge, unit.elementary_charge)
        sigma[idx] = _value(s, unit.nanometer)
        epsilon[idx] = _value(e, unit.kilojoule_per_mole)
    return charges, sigma, epsilon


def _custom_nonbonded_exclusions_hash(system) -> str:
    exclusions: list[tuple] = []
    custom_nonbonded = _optional_force(system, "CustomNonbondedForce")
    if custom_nonbonded is not None:
        for idx in range(custom_nonbonded.getNumExclusions()):
            a1, a2 = custom_nonbonded.getExclusionParticles(idx)
            exclusions.append((min(int(a1), int(a2)), max(int(a1), int(a2))))
    return _digest_rows(exclusions)


def _build_openmm_candidate(name: str):
    pdb_path = FIXTURE_ROOT / f"{name}_processed.pdb"
    pdb = PDBFile(str(pdb_path))
    topology = pdb.topology
    restore_pdb_atom_names(topology, pdb_path)
    positions = list(pdb.positions)
    add_adjacent_peptide_bonds(topology)
    add_disulfide_bonds(topology, positions)
    rename_solvent_for_forcefield_input(topology)

    ff = load_a99sbdisp_forcefield()
    residue_templates = assign_residue_templates(topology)
    system = ff.createSystem(
        topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=HBonds,
        residueTemplates=residue_templates,
    )
    added_impropers = add_explicit_rtp_impropers(system, topology, residue_templates)
    return pdb, system, added_impropers


@pytest.mark.parametrize("name", PARITY_CASES)
def test_a99sbdisp_openmm_matches_checked_in_legacy_baselines(name: str) -> None:
    baseline = json.loads((FIXTURE_ROOT / f"{name}_baseline.json").read_text(encoding="utf-8"))
    pdb, candidate, added_impropers = _build_openmm_candidate(name)

    assert candidate.getNumParticles() == baseline["num_particles"]
    assert candidate.getNumConstraints() == baseline["num_constraints"]
    assert added_impropers == baseline["added_impropers"]
    assert _custom_nonbonded_exclusions_hash(candidate) == baseline["custom_nonbonded_exclusions_hash"]

    nonbonded = np.load(FIXTURE_ROOT / f"{name}_nonbonded_particles.npz")
    q_ref = nonbonded["q"]
    sigma_ref = nonbonded["sigma"]
    epsilon_ref = nonbonded["epsilon"]
    q_candidate, sigma_candidate, epsilon_candidate = _nonbonded_particle_arrays(candidate)

    assert q_candidate.shape == q_ref.shape
    assert sigma_candidate.shape == sigma_ref.shape
    assert epsilon_candidate.shape == epsilon_ref.shape
    assert np.max(np.abs(q_candidate - q_ref)) <= 1e-5
    assert np.max(np.abs(sigma_candidate - sigma_ref)) == 0.0
    assert np.max(np.abs(epsilon_candidate - epsilon_ref)) == 0.0

    platform = Platform.getPlatformByName("CPU")
    ctx_candidate = Context(candidate, VerletIntegrator(0.001 * unit.picoseconds), platform)
    ctx_candidate.setPositions(pdb.positions)
    energy_candidate = _value(ctx_candidate.getState(getEnergy=True).getPotentialEnergy(), unit.kilojoule_per_mole)
    energy_baseline = float(baseline["potential_energy_kj_per_mol"])
    relative_delta = abs(energy_candidate - energy_baseline) / max(1.0, abs(energy_baseline))
    assert relative_delta <= 1e-3


def test_prepare_with_openmm_tiny_smoke(tmp_path: Path) -> None:
    input_pdb = REPO_ROOT / "tests" / "fixtures" / "tiny_protein.pdb"
    output_pdb = tmp_path / "tiny_processed.pdb"
    output_xml = tmp_path / "tiny.xml"

    subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "hydromap" / "engines" / "simulation" / "prepare_with_openmm.py"),
            str(input_pdb),
            "--output_pdb",
            str(output_pdb),
            "--output_system",
            str(output_xml),
            "--seed",
            "7",
        ],
        check=True,
        cwd=REPO_ROOT,
    )

    assert output_pdb.exists()
    assert output_xml.exists()

    pdb_lines = output_pdb.read_text(encoding="utf-8").splitlines()
    solvent_lines = [line for line in pdb_lines if line.startswith(("ATOM  ", "HETATM")) and line[17:20].strip() == "SOL"]
    assert solvent_lines
    assert [line[12:16].strip() for line in solvent_lines[:4]] == ["OW", "HW1", "HW2", "MW"]

    with output_xml.open("r", encoding="utf-8") as handle:
        system = XmlSerializer.deserialize(handle.read())
    atom_count = sum(1 for line in pdb_lines if line.startswith(("ATOM  ", "HETATM")))
    assert system.getNumParticles() == atom_count
