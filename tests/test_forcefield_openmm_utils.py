import io
import random
from pathlib import Path

from openmm import Vec3
from openmm.unit import molar, nanometer
from openmm.app import PDBFile, Topology, element as elem

from hydromap.engines.simulation.prepare_with_openmm import (
    add_hydrogens_with_fallback,
    choose_hydrogen_variants,
    remove_pdb_conect_records,
    residue_templates_for_add_hydrogens,
)
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
)
from hydromap.forcefield.openmm_utils import infer_variant_name, restore_pdb_atom_names


def _build_residue(name: str, atom_names: list[str]) -> Topology:
    topology = Topology()
    chain = topology.addChain("A")
    residue = topology.addResidue(name, chain, "1")
    for atom_name in atom_names:
        element = elem.hydrogen if atom_name.startswith("H") else elem.carbon
        if atom_name.startswith("N"):
            element = elem.nitrogen
        elif atom_name.startswith("O"):
            element = elem.oxygen
        elif atom_name.startswith("S"):
            element = elem.sulfur
        topology.addAtom(atom_name, element, residue)
    return topology


def test_remove_pdb_conect_records_preserves_coordinate_records(tmp_path: Path) -> None:
    pdb_path = tmp_path / "large_system.pdb"
    original = (
        "ATOM  A0000  O   SOL A   1       0.000   0.000   0.000  1.00  0.00           O\n"
        "CONECTA0000A0001\n"
        "CONECT    1    2\n"
        "END\n"
    )
    pdb_path.write_text(original, encoding="utf-8")

    assert remove_pdb_conect_records(pdb_path) == 2
    assert pdb_path.read_text(encoding="utf-8") == original.split("CONECT", 1)[0] + "END\n"
    assert remove_pdb_conect_records(pdb_path) == 0


def test_infer_variant_name_defaults_lys_without_hydrogens_to_lys() -> None:
    topology = _build_residue("LYS", ["N", "CA", "CB", "CG", "CD", "CE", "NZ"])
    residue = next(topology.residues())
    assert infer_variant_name(residue) == "LYS"


def test_infer_variant_name_defaults_cys_without_hydrogens_to_cys() -> None:
    topology = _build_residue("CYS", ["N", "CA", "CB", "SG"])
    residue = next(topology.residues())
    assert infer_variant_name(residue) == "CYS"


def test_infer_variant_name_keeps_cyx_for_disulfide() -> None:
    topology = Topology()
    chain = topology.addChain("A")
    res1 = topology.addResidue("CYS", chain, "1")
    res2 = topology.addResidue("CYS", chain, "2")
    sg1 = topology.addAtom("SG", elem.sulfur, res1)
    sg2 = topology.addAtom("SG", elem.sulfur, res2)
    topology.addBond(sg1, sg2)
    assert infer_variant_name(res1) == "CYX"
    assert infer_variant_name(res2) == "CYX"


def test_choose_hydrogen_variants_does_not_force_cys_variant() -> None:
    topology = _build_residue("CYS", ["N", "CA", "CB", "SG"])
    assert choose_hydrogen_variants(topology) == [None]


def test_choose_hydrogen_variants_can_globally_force_hip() -> None:
    topology = _build_residue("HIS", ["N", "CA", "CB", "CG", "ND1", "CD2", "CE1", "NE2"])
    assert choose_hydrogen_variants(topology, histidine_mode="hip") == ["HIP"]


def test_choose_hydrogen_variants_histidine_override_beats_global_mode() -> None:
    topology = Topology()
    chain = topology.addChain("B")
    residue = topology.addResidue("HIS", chain, "417")
    for atom_name, element in [
        ("N", elem.nitrogen),
        ("CA", elem.carbon),
        ("CB", elem.carbon),
        ("CG", elem.carbon),
        ("ND1", elem.nitrogen),
        ("CD2", elem.carbon),
        ("CE1", elem.carbon),
        ("NE2", elem.nitrogen),
    ]:
        topology.addAtom(atom_name, element, residue)

    assert choose_hydrogen_variants(
        topology,
        histidine_mode="hie",
        histidine_overrides={"B:417": "HIP"},
    ) == ["HIP"]


def test_restore_pdb_atom_names_ignores_altloc_duplicates(tmp_path: Path) -> None:
    pdb_text = """\
ATOM      1  N   ARG A   1      17.260  26.959  16.899  1.00 27.66           N
ATOM      2  CA  ARG A   1      18.612  26.675  17.351  1.00 27.90           C
ATOM      3  C   ARG A   1      19.326  28.000  17.050  1.00 28.53           C
ATOM      4  O   ARG A   1      18.867  29.082  17.461  1.00 30.10           O
ATOM      5  CB AARG A   1      18.670  26.364  18.891  0.50 27.40           C
ATOM      6  CB BARG A   1      18.527  26.310  18.839  0.50 25.76           C
TER
END
"""
    pdb_path = tmp_path / "altloc.pdb"
    pdb_path.write_text(pdb_text, encoding="utf-8")
    pdb = PDBFile(str(pdb_path))
    restore_pdb_atom_names(pdb.topology, pdb_path)
    assert [atom.name for atom in pdb.topology.atoms()] == ["N", "CA", "C", "O", "CB"]


def test_add_hydrogens_with_fallback_retries_nan(monkeypatch) -> None:
    topology = _build_residue("ALA", ["N", "CA", "C", "O", "CB"])
    positions = [Vec3(0.0, 0.0, 0.0) for _ in topology.atoms()]

    class FakeModeller:
        calls = 0

        def __init__(self, topology, positions):
            self.topology = topology
            self.positions = positions

        def addHydrogens(self, ff, variants=None, residueTemplates=None):
            FakeModeller.calls += 1
            if FakeModeller.calls == 1:
                raise Exception("Particle coordinate is NaN.  For more information, see https://github.com/openmm/openmm/wiki/Frequently-Asked-Questions#nan")

    monkeypatch.setattr("hydromap.engines.simulation.prepare_with_openmm.Modeller", FakeModeller)

    modeller = add_hydrogens_with_fallback(
        ff=object(),
        topology=topology,
        positions=positions * nanometer,
        residue_templates={next(topology.residues()): "ALA"},
        seed=123,
    )

    assert isinstance(modeller, FakeModeller)
    assert FakeModeller.calls == 2


def test_add_hydrogens_with_fallback_supports_capped_peptide(tmp_path: Path) -> None:
    pdb_text = """\
ATOM      1  C   ACE A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  O   ACE A   1       1.220   0.000   0.000  1.00  0.00           O
ATOM      3  CH3 ACE A   1      -0.750   1.200   0.000  1.00  0.00           C
ATOM      4  N   GLY A   2      -0.700  -1.150   0.000  1.00  0.00           N
ATOM      5  CA  GLY A   2      -0.050  -2.440   0.000  1.00  0.00           C
ATOM      6  C   GLY A   2       1.430  -2.250   0.000  1.00  0.00           C
ATOM      7  O   GLY A   2       2.150  -3.220   0.000  1.00  0.00           O
ATOM      8  N   NME A   3       1.950  -1.020   0.000  1.00  0.00           N
ATOM      9  CH3 NME A   3       3.360  -0.760   0.000  1.00  0.00           C
TER
END
"""
    pdb_path = tmp_path / "capped_gly.pdb"
    pdb_path.write_text(pdb_text, encoding="utf-8")

    ff = load_a99sbdisp_forcefield()
    pdb = PDBFile(str(pdb_path))
    topology = pdb.topology
    restore_pdb_atom_names(topology, pdb_path)
    positions = [Vec3(*pos.value_in_unit(nanometer)) for pos in pdb.positions]

    add_adjacent_peptide_bonds(topology)
    add_disulfide_bonds(topology, positions)
    topology, positions, added = ensure_c_terminal_oxygens(topology, positions)
    assert added == 0
    rename_terminal_oxygen_atoms(topology)
    canonicalize_residue_names_for_hydrogens(topology)

    modeller = add_hydrogens_with_fallback(
        ff=ff,
        topology=topology,
        positions=positions * nanometer,
        residue_templates=residue_templates_for_add_hydrogens(topology),
        seed=1,
    )

    residues = list(modeller.topology.residues())
    assert [atom.name for atom in residues[0].atoms()] == ["C", "O", "CH3", "HH31", "HH32", "HH33"]
    assert [atom.name for atom in residues[-1].atoms()] == ["N", "H", "CH3", "HH31", "HH32", "HH33"]

    system = ff.createSystem(modeller.topology, residueTemplates=assign_residue_templates(modeller.topology))
    assert system.getNumParticles() == modeller.topology.getNumAtoms()


def test_prepare_pipeline_is_reproducible_for_same_seed() -> None:
    fixture = Path(__file__).resolve().parent / "fixtures" / "tiny_protein.pdb"
    ff = load_a99sbdisp_forcefield()

    def prepare_once(seed: int) -> str:
        pdb = PDBFile(str(fixture))
        topology = pdb.topology
        restore_pdb_atom_names(topology, fixture)
        positions = [Vec3(*pos.value_in_unit(nanometer)) for pos in pdb.positions]

        add_disulfide_bonds(topology, positions)
        topology, positions, _ = prepare_peptide_topology(topology, positions, capping_mode="none")
        add_adjacent_peptide_bonds(topology)
        topology, positions, _ = ensure_c_terminal_oxygens(topology, positions)
        rename_terminal_oxygen_atoms(topology)
        canonicalize_residue_names_for_hydrogens(topology)

        random.seed(seed)
        modeller = add_hydrogens_with_fallback(
            ff=ff,
            topology=topology,
            positions=positions * nanometer,
            residue_templates=residue_templates_for_add_hydrogens(topology),
            seed=seed,
        )

        solute_templates = assign_residue_templates(modeller.topology)
        box_vectors = compute_orthorhombic_box(modeller.topology, modeller.positions, padding_nm=1.0)

        random.seed(seed)
        modeller.addSolvent(
            ff,
            model="tip4pew",
            boxVectors=box_vectors,
            positiveIon="Na+",
            negativeIon="Cl-",
            ionicStrength=0.0 * molar,
            neutralize=True,
            residueTemplates=solute_templates,
        )

        rename_solvent_for_forcefield_input(modeller.topology)
        final_templates = assign_residue_templates(modeller.topology)
        system = ff.createSystem(modeller.topology, residueTemplates=final_templates)
        add_explicit_rtp_impropers(system, modeller.topology, final_templates)
        rename_solvent_for_output(modeller.topology)

        handle = io.StringIO()
        PDBFile.writeFile(modeller.topology, modeller.positions, handle)
        return handle.getvalue()

    assert prepare_once(123) == prepare_once(123)


def test_ensure_c_terminal_oxygens_handles_multichain_topology() -> None:
    pdb_text = """\
ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  GLY A   1       1.450   0.000   0.000  1.00  0.00           C
ATOM      3  C   GLY A   1       2.100   1.300   0.000  1.00  0.00           C
ATOM      4  O   GLY A   1       1.600   2.350   0.000  1.00  0.00           O
TER
ATOM      5  N   ALA B   1       5.000   0.000   0.000  1.00  0.00           N
ATOM      6  CA  ALA B   1       6.450   0.000   0.000  1.00  0.00           C
ATOM      7  C   ALA B   1       7.100   1.300   0.000  1.00  0.00           C
ATOM      8  O   ALA B   1       6.600   2.350   0.000  1.00  0.00           O
ATOM      9  CB  ALA B   1       6.950  -0.850  -1.200  1.00  0.00           C
END
"""
    handle = io.StringIO(pdb_text)
    pdb = PDBFile(handle)
    topology = pdb.topology
    positions = [Vec3(*pos.value_in_unit(nanometer)) for pos in pdb.positions]

    topology, positions, added = ensure_c_terminal_oxygens(topology, positions)

    assert added == 2
    residues = list(topology.residues())
    assert [atom.name for atom in residues[0].atoms()] == ["N", "CA", "C", "O", "OXT"]
    assert [atom.name for atom in residues[1].atoms()] == ["N", "CA", "C", "O", "CB", "OXT"]


def test_prepare_peptide_topology_caps_true_termini(tmp_path: Path) -> None:
    fixture = Path(__file__).resolve().parent / "fixtures" / "tiny_protein.pdb"
    ff = load_a99sbdisp_forcefield()
    pdb = PDBFile(str(fixture))
    topology = pdb.topology
    restore_pdb_atom_names(topology, fixture)
    positions = [Vec3(*pos.value_in_unit(nanometer)) for pos in pdb.positions]

    add_disulfide_bonds(topology, positions)
    topology, positions, audit = prepare_peptide_topology(topology, positions, capping_mode="termini")
    add_adjacent_peptide_bonds(topology)
    topology, positions, added_oxt = ensure_c_terminal_oxygens(topology, positions)
    rename_terminal_oxygen_atoms(topology)
    canonicalize_residue_names_for_hydrogens(topology)

    residues = list(topology.residues())
    assert [res.name for res in residues] == ["ACE", "ALA", "GLY", "NME"]
    assert audit["added_ace_cap_count"] == 1
    assert audit["added_nme_cap_count"] == 1
    assert added_oxt == 0

    modeller = add_hydrogens_with_fallback(
        ff=ff,
        topology=topology,
        positions=positions * nanometer,
        residue_templates=residue_templates_for_add_hydrogens(topology),
        seed=7,
    )
    system = ff.createSystem(modeller.topology, residueTemplates=assign_residue_templates(modeller.topology))
    assert system.getNumParticles() == modeller.topology.getNumAtoms()


def test_prepare_peptide_topology_splits_and_caps_internal_breaks() -> None:
    pdb_text = """\
ATOM      1  N   ALA A   1       8.000   8.000   8.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       9.300   8.400   8.200  1.00  0.00           C
ATOM      3  C   ALA A   1      10.100   7.300   8.900  1.00  0.00           C
ATOM      4  O   ALA A   1       9.900   6.100   8.700  1.00  0.00           O
ATOM      5  CB  ALA A   1       9.800   9.700   7.500  1.00  0.00           C
ATOM      6  N   GLY A   2      11.100   7.700   9.700  1.00  0.00           N
ATOM      7  CA  GLY A   2      12.000   6.700  10.300  1.00  0.00           C
ATOM      8  C   GLY A   2      11.400   5.300  10.600  1.00  0.00           C
ATOM      9  O   GLY A   2      12.100   4.300  10.600  1.00  0.00           O
ATOM     10  N   ALA A   5      18.000   8.000   8.000  1.00  0.00           N
ATOM     11  CA  ALA A   5      19.300   8.400   8.200  1.00  0.00           C
ATOM     12  C   ALA A   5      20.100   7.300   8.900  1.00  0.00           C
ATOM     13  O   ALA A   5      19.900   6.100   8.700  1.00  0.00           O
ATOM     14  CB  ALA A   5      19.800   9.700   7.500  1.00  0.00           C
ATOM     15  N   GLY A   6      21.100   7.700   9.700  1.00  0.00           N
ATOM     16  CA  GLY A   6      22.000   6.700  10.300  1.00  0.00           C
ATOM     17  C   GLY A   6      21.400   5.300  10.600  1.00  0.00           C
ATOM     18  O   GLY A   6      22.100   4.300  10.600  1.00  0.00           O
END
"""
    pdb = PDBFile(io.StringIO(pdb_text))
    topology = pdb.topology
    positions = [Vec3(*pos.value_in_unit(nanometer)) for pos in pdb.positions]

    topology, positions, audit = prepare_peptide_topology(topology, positions, capping_mode="breaks")
    add_adjacent_peptide_bonds(topology)
    ff = load_a99sbdisp_forcefield()
    topology, positions, added_oxt = ensure_c_terminal_oxygens(topology, positions)
    rename_terminal_oxygen_atoms(topology)
    canonicalize_residue_names_for_hydrogens(topology)

    chains = list(topology.chains())
    assert len(chains) == 2
    assert [res.name for res in chains[0].residues()] == ["ALA", "GLY", "NME"]
    assert [res.name for res in chains[1].residues()] == ["ACE", "ALA", "GLY"]
    assert audit["detected_peptide_break_count"] == 1
    assert audit["added_ace_cap_count"] == 1
    assert audit["added_nme_cap_count"] == 1
    assert added_oxt == 1

    modeller = add_hydrogens_with_fallback(
        ff=ff,
        topology=topology,
        positions=positions * nanometer,
        residue_templates=residue_templates_for_add_hydrogens(topology),
        seed=9,
    )
    system = ff.createSystem(modeller.topology, residueTemplates=assign_residue_templates(modeller.topology))
    assert system.getNumParticles() == modeller.topology.getNumAtoms()


def test_prepare_peptide_topology_splits_internal_breaks_without_capping() -> None:
    pdb_text = """\
ATOM      1  N   ALA A   1       8.000   8.000   8.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       9.300   8.400   8.200  1.00  0.00           C
ATOM      3  C   ALA A   1      10.100   7.300   8.900  1.00  0.00           C
ATOM      4  O   ALA A   1       9.900   6.100   8.700  1.00  0.00           O
ATOM      5  CB  ALA A   1       9.800   9.700   7.500  1.00  0.00           C
ATOM      6  N   GLY A   2      11.100   7.700   9.700  1.00  0.00           N
ATOM      7  CA  GLY A   2      12.000   6.700  10.300  1.00  0.00           C
ATOM      8  C   GLY A   2      11.400   5.300  10.600  1.00  0.00           C
ATOM      9  O   GLY A   2      12.100   4.300  10.600  1.00  0.00           O
ATOM     10  N   ALA A   5      18.000   8.000   8.000  1.00  0.00           N
ATOM     11  CA  ALA A   5      19.300   8.400   8.200  1.00  0.00           C
ATOM     12  C   ALA A   5      20.100   7.300   8.900  1.00  0.00           C
ATOM     13  O   ALA A   5      19.900   6.100   8.700  1.00  0.00           O
ATOM     14  CB  ALA A   5      19.800   9.700   7.500  1.00  0.00           C
ATOM     15  N   GLY A   6      21.100   7.700   9.700  1.00  0.00           N
ATOM     16  CA  GLY A   6      22.000   6.700  10.300  1.00  0.00           C
ATOM     17  C   GLY A   6      21.400   5.300  10.600  1.00  0.00           C
ATOM     18  O   GLY A   6      22.100   4.300  10.600  1.00  0.00           O
END
"""
    pdb = PDBFile(io.StringIO(pdb_text))
    topology = pdb.topology
    positions = [Vec3(*pos.value_in_unit(nanometer)) for pos in pdb.positions]

    topology, positions, audit = prepare_peptide_topology(topology, positions, capping_mode="none")
    added = add_adjacent_peptide_bonds(topology)

    assert len(list(topology.chains())) == 2
    assert added == 0
    assert audit["detected_peptide_break_count"] == 1
