from pathlib import Path

import pytest

from hydromap.config import load_config, package_repo_root
from hydromap.errors import ConfigError


REPO_ROOT = Path(__file__).resolve().parents[1]


def write_config(tmp_path: Path, body: str) -> Path:
    cfg = tmp_path / "config.yaml"
    cfg.write_text(body, encoding="utf-8")
    return cfg


def test_load_config_defaults_and_lists(tmp_path: Path) -> None:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture_dir = REPO_ROOT / "tests" / "fixtures"

    cfg_path = write_config(
        tmp_path,
        f"""
input_dir: {fixture_dir}
protein: tiny_protein
seed: 11
model_path: {model}
""",
    )

    cfg = load_config(cfg_path, repo_root=REPO_ROOT)
    assert cfg.proteins == ["tiny_protein"]
    assert cfg.seeds == [11]
    assert cfg.execution.profile == "balanced"
    assert cfg.md.allow_cpu_md is False
    assert cfg.md.nanoseconds == pytest.approx(0.5)
    assert cfg.md.equilibration_ns == pytest.approx(0.1)
    assert cfg.md.device == "gpu"
    assert cfg.md.timestep_ps == pytest.approx(0.003)
    assert cfg.md.report_interval_ps == pytest.approx(1.0)
    assert cfg.md.neutralize is True
    assert cfg.md.ionic_strength_molar == pytest.approx(0.0)
    assert cfg.md.positive_ion == "Na+"
    assert cfg.md.negative_ion == "Cl-"
    assert cfg.md.histidine_mode == "auto"
    assert cfg.md.histidine_overrides == {}
    assert cfg.md.repair_missing_atoms == "none"
    assert cfg.md.capping_mode == "none"
    assert cfg.md.prep_policy == "permissive"
    assert cfg.analysis.potentials_device == "gpu"
    assert cfg.analysis.tail_ns is None
    assert cfg.analysis.triplets_frame_stride == 1
    assert cfg.analysis.potentials_frame_stride == 100


def test_multi_chain_removed(tmp_path: Path) -> None:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture_dir = REPO_ROOT / "tests" / "fixtures"

    cfg_path = write_config(
        tmp_path,
        f"""
input_dir: {fixture_dir}
protein: tiny_protein
seed: 11
model_path: {model}
multi_chain: true
""",
    )

    with pytest.raises(ConfigError, match="multi_chain"):
        load_config(cfg_path, repo_root=REPO_ROOT)


def test_invalid_execution_profile(tmp_path: Path) -> None:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture_dir = REPO_ROOT / "tests" / "fixtures"

    cfg_path = write_config(
        tmp_path,
        f"""
input_dir: {fixture_dir}
protein: tiny_protein
seed: 11
model_path: {model}
execution:
  profile: bad_profile
""",
    )

    with pytest.raises(ConfigError, match="execution.profile"):
        load_config(cfg_path, repo_root=REPO_ROOT)


def test_invalid_potentials_device(tmp_path: Path) -> None:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture_dir = REPO_ROOT / "tests" / "fixtures"

    cfg_path = write_config(
        tmp_path,
        f"""
input_dir: {fixture_dir}
protein: tiny_protein
seed: 11
model_path: {model}
analysis:
  potentials_device: bad_device
""",
    )

    with pytest.raises(ConfigError, match="analysis.potentials_device"):
        load_config(cfg_path, repo_root=REPO_ROOT)


def test_invalid_triplets_device(tmp_path: Path) -> None:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture_dir = REPO_ROOT / "tests" / "fixtures"

    cfg_path = write_config(
        tmp_path,
        f"""
input_dir: {fixture_dir}
protein: tiny_protein
seed: 11
model_path: {model}
analysis:
  triplets_device: bad_device
""",
    )

    with pytest.raises(ConfigError, match="analysis.triplets_device"):
        load_config(cfg_path, repo_root=REPO_ROOT)


def test_legacy_backend_alias_maps_to_device(tmp_path: Path) -> None:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture_dir = REPO_ROOT / "tests" / "fixtures"

    cfg_path = write_config(
        tmp_path,
        f"""
input_dir: {fixture_dir}
protein: tiny_protein
seed: 11
model_path: {model}
analysis:
  triplets_backend: cpu_parallel
  potentials_backend: gpu_serial
""",
    )

    cfg = load_config(cfg_path, repo_root=REPO_ROOT)
    assert cfg.analysis.triplets_device == "cpu"
    assert cfg.analysis.potentials_device == "gpu"


def test_analysis_time_aliases_map_to_tail_fields(tmp_path: Path) -> None:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture_dir = REPO_ROOT / "tests" / "fixtures"

    cfg_path = write_config(
        tmp_path,
        f"""
input_dir: {fixture_dir}
protein: tiny_protein
seed: 11
model_path: {model}
analysis:
  triplets_time_ns: 3
  potentials_time_ns: 2
""",
    )

    cfg = load_config(cfg_path, repo_root=REPO_ROOT)
    assert cfg.analysis.triplets_tail_ns == 3
    assert cfg.analysis.potentials_tail_ns == 2


def test_analysis_stride_aliases_map_to_frame_stride(tmp_path: Path) -> None:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture_dir = REPO_ROOT / "tests" / "fixtures"

    cfg_path = write_config(
        tmp_path,
        f"""
input_dir: {fixture_dir}
protein: tiny_protein
seed: 11
model_path: {model}
analysis:
  triplets_skip: 7
  potentials_skip: 55
""",
    )

    cfg = load_config(cfg_path, repo_root=REPO_ROOT)
    assert cfg.analysis.triplets_frame_stride == 7
    assert cfg.analysis.potentials_frame_stride == 55


def test_restrain_selection_alias(tmp_path: Path) -> None:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture_dir = REPO_ROOT / "tests" / "fixtures"

    cfg_path = write_config(
        tmp_path,
        f"""
input_dir: {fixture_dir}
protein: tiny_protein
seed: 11
model_path: {model}
md:
  restrain_selection: "backbone and name CA"
""",
    )

    cfg = load_config(cfg_path, repo_root=REPO_ROOT)
    assert cfg.md.restrain == "backbone and name CA"


def test_restrain_and_restrain_selection_conflict(tmp_path: Path) -> None:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture_dir = REPO_ROOT / "tests" / "fixtures"

    cfg_path = write_config(
        tmp_path,
        f"""
input_dir: {fixture_dir}
protein: tiny_protein
seed: 11
model_path: {model}
md:
  restrain: "protein"
  restrain_selection: "backbone"
""",
    )

    with pytest.raises(ConfigError, match="restrain"):
        load_config(cfg_path, repo_root=REPO_ROOT)


def test_md_ion_controls_are_loaded(tmp_path: Path) -> None:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture_dir = REPO_ROOT / "tests" / "fixtures"

    cfg_path = write_config(
        tmp_path,
        f"""
input_dir: {fixture_dir}
protein: tiny_protein
seed: 11
model_path: {model}
md:
  neutralize: false
  ionic_strength_molar: 0.15
  positive_ion: K+
  negative_ion: Br-
""",
    )

    cfg = load_config(cfg_path, repo_root=REPO_ROOT)
    assert cfg.md.neutralize is False
    assert cfg.md.ionic_strength_molar == pytest.approx(0.15)
    assert cfg.md.positive_ion == "K+"
    assert cfg.md.negative_ion == "Br-"


def test_histidine_controls_are_loaded(tmp_path: Path) -> None:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture_dir = REPO_ROOT / "tests" / "fixtures"

    cfg_path = write_config(
        tmp_path,
        f"""
input_dir: {fixture_dir}
protein: tiny_protein
seed: 11
model_path: {model}
md:
  histidine_mode: hip
  histidine_overrides:
    "B:417": HID
    "505": HIP
""",
    )

    cfg = load_config(cfg_path, repo_root=REPO_ROOT)
    assert cfg.md.histidine_mode == "hip"
    assert cfg.md.histidine_overrides == {"B:417": "HID", "505": "HIP"}


def test_repair_missing_atoms_mode_is_loaded(tmp_path: Path) -> None:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture_dir = REPO_ROOT / "tests" / "fixtures"

    cfg_path = write_config(
        tmp_path,
        f"""
input_dir: {fixture_dir}
protein: tiny_protein
seed: 11
model_path: {model}
md:
  repair_missing_atoms: pdbfixer
""",
    )

    cfg = load_config(cfg_path, repo_root=REPO_ROOT)
    assert cfg.md.repair_missing_atoms == "pdbfixer"


def test_capping_mode_is_loaded(tmp_path: Path) -> None:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture_dir = REPO_ROOT / "tests" / "fixtures"

    cfg_path = write_config(
        tmp_path,
        f"""
input_dir: {fixture_dir}
protein: tiny_protein
seed: 11
model_path: {model}
md:
  capping_mode: termini_and_breaks
""",
    )

    cfg = load_config(cfg_path, repo_root=REPO_ROOT)
    assert cfg.md.capping_mode == "termini_and_breaks"


def test_prep_policy_is_loaded(tmp_path: Path) -> None:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture_dir = REPO_ROOT / "tests" / "fixtures"

    cfg_path = write_config(
        tmp_path,
        f"""
input_dir: {fixture_dir}
protein: tiny_protein
seed: 11
model_path: {model}
md:
  prep_policy: strict
""",
    )

    cfg = load_config(cfg_path, repo_root=REPO_ROOT)
    assert cfg.md.prep_policy == "strict"


def test_invalid_positive_ion_is_rejected(tmp_path: Path) -> None:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture_dir = REPO_ROOT / "tests" / "fixtures"

    cfg_path = write_config(
        tmp_path,
        f"""
input_dir: {fixture_dir}
protein: tiny_protein
seed: 11
model_path: {model}
md:
  positive_ion: Ca2+
""",
    )

    with pytest.raises(ConfigError, match="md.positive_ion"):
        load_config(cfg_path, repo_root=REPO_ROOT)


def test_invalid_negative_ion_is_rejected(tmp_path: Path) -> None:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture_dir = REPO_ROOT / "tests" / "fixtures"

    cfg_path = write_config(
        tmp_path,
        f"""
input_dir: {fixture_dir}
protein: tiny_protein
seed: 11
model_path: {model}
md:
  negative_ion: SO4--
""",
    )

    with pytest.raises(ConfigError, match="md.negative_ion"):
        load_config(cfg_path, repo_root=REPO_ROOT)


def test_invalid_histidine_mode_is_rejected(tmp_path: Path) -> None:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture_dir = REPO_ROOT / "tests" / "fixtures"

    cfg_path = write_config(
        tmp_path,
        f"""
input_dir: {fixture_dir}
protein: tiny_protein
seed: 11
model_path: {model}
md:
  histidine_mode: charged
""",
    )

    with pytest.raises(ConfigError, match="md.histidine_mode"):
        load_config(cfg_path, repo_root=REPO_ROOT)


def test_invalid_repair_missing_atoms_mode_is_rejected(tmp_path: Path) -> None:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture_dir = REPO_ROOT / "tests" / "fixtures"

    cfg_path = write_config(
        tmp_path,
        f"""
input_dir: {fixture_dir}
protein: tiny_protein
seed: 11
model_path: {model}
md:
  repair_missing_atoms: auto
""",
    )

    with pytest.raises(ConfigError, match="md.repair_missing_atoms"):
        load_config(cfg_path, repo_root=REPO_ROOT)


def test_invalid_capping_mode_is_rejected(tmp_path: Path) -> None:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture_dir = REPO_ROOT / "tests" / "fixtures"

    cfg_path = write_config(
        tmp_path,
        f"""
input_dir: {fixture_dir}
protein: tiny_protein
seed: 11
model_path: {model}
md:
  capping_mode: auto
""",
    )

    with pytest.raises(ConfigError, match="md.capping_mode"):
        load_config(cfg_path, repo_root=REPO_ROOT)


def test_invalid_prep_policy_is_rejected(tmp_path: Path) -> None:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture_dir = REPO_ROOT / "tests" / "fixtures"

    cfg_path = write_config(
        tmp_path,
        f"""
input_dir: {fixture_dir}
protein: tiny_protein
seed: 11
model_path: {model}
md:
  prep_policy: auto
""",
    )

    with pytest.raises(ConfigError, match="md.prep_policy"):
        load_config(cfg_path, repo_root=REPO_ROOT)


def test_invalid_histidine_override_is_rejected(tmp_path: Path) -> None:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture_dir = REPO_ROOT / "tests" / "fixtures"

    cfg_path = write_config(
        tmp_path,
        f"""
input_dir: {fixture_dir}
protein: tiny_protein
seed: 11
model_path: {model}
md:
  histidine_overrides:
    "B:417": HSP
""",
    )

    with pytest.raises(ConfigError, match="histidine_overrides"):
        load_config(cfg_path, repo_root=REPO_ROOT)


def test_negative_ionic_strength_is_rejected(tmp_path: Path) -> None:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture_dir = REPO_ROOT / "tests" / "fixtures"

    cfg_path = write_config(
        tmp_path,
        f"""
input_dir: {fixture_dir}
protein: tiny_protein
seed: 11
model_path: {model}
md:
  ionic_strength_molar: -0.1
""",
    )

    with pytest.raises(ConfigError, match="md.ionic_strength_molar"):
        load_config(cfg_path, repo_root=REPO_ROOT)


def test_existing_analysis_inputs_require_all_fields(tmp_path: Path) -> None:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture_dir = REPO_ROOT / "tests" / "fixtures"

    cfg_path = write_config(
        tmp_path,
        f"""
input_dir: {fixture_dir}
protein: tiny_protein
seed: 11
model_path: {model}
analysis:
  existing_trajectory: /tmp/traj.dcd
""",
    )

    with pytest.raises(ConfigError, match="existing_processed_pdb and analysis.existing_trajectory"):
        load_config(cfg_path, repo_root=REPO_ROOT)


def test_prepare_backend_gromacs_is_rejected(tmp_path: Path) -> None:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture_dir = REPO_ROOT / "tests" / "fixtures"

    cfg_path = write_config(
        tmp_path,
        f"""
input_dir: {fixture_dir}
protein: tiny_protein
seed: 11
model_path: {model}
md:
  prepare_backend: gromacs
""",
    )

    with pytest.raises(ConfigError, match="OpenMM preparation path only"):
        load_config(cfg_path, repo_root=REPO_ROOT)


def test_existing_topology_may_be_gromacs_top(tmp_path: Path) -> None:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture_dir = REPO_ROOT / "tests" / "fixtures"

    cfg_path = write_config(
        tmp_path,
        f"""
input_dir: {fixture_dir}
protein: tiny_protein
seed: 11
model_path: {model}
analysis:
  existing_processed_pdb: /tmp/protein_processed.pdb
  existing_trajectory: /tmp/traj.dcd
  existing_topology: /tmp/topol.top
""",
    )

    cfg = load_config(cfg_path, repo_root=REPO_ROOT)
    assert cfg.analysis.existing_topology == "/tmp/topol.top"


def test_existing_topology_rejects_unknown_suffix(tmp_path: Path) -> None:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture_dir = REPO_ROOT / "tests" / "fixtures"

    cfg_path = write_config(
        tmp_path,
        f"""
input_dir: {fixture_dir}
protein: tiny_protein
seed: 11
model_path: {model}
analysis:
  existing_processed_pdb: /tmp/protein_processed.pdb
  existing_trajectory: /tmp/traj.dcd
  existing_topology: /tmp/topology.psf
""",
    )

    with pytest.raises(ConfigError, match="OpenMM System XML file or a GROMACS .top file"):
        load_config(cfg_path, repo_root=REPO_ROOT)


def test_load_config_defaults_repo_root_to_package_root(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    model = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
    fixture_dir = REPO_ROOT / "tests" / "fixtures"

    cfg_path = write_config(
        tmp_path,
        f"""
input_dir: {fixture_dir}
protein: tiny_protein
seed: 11
model_path: {model}
""",
    )

    monkeypatch.chdir(tmp_path)
    cfg = load_config(cfg_path)
    assert cfg.repo_root == package_repo_root()
