# HydroMap

Modeling water-protein interactions to map local water structure, water-protein potential energy, and dewetting free energy.

<p align="center">
  <img
    src="https://raw.githubusercontent.com/samlobe/HydroMap/main/images/HydroMap_image.png"
    alt="HydroMap overview"
    width="720"
  />
</p>

---

## What It Does

HydroMap maps local hydration around proteins from explicit-solvent molecular dynamics. For each residue (or user-defined group), it computes:

- water triplet-angle distributions (projected into PCs)
- water-protein interaction energies
- predicted dewetting free energy (`Fdewet`)

## Input / Output

**Input**
- protein structure file (pdb)

**Output**
- CSV with per-group features and predictions
- colored PDB files (properties in B-factor/tempfactor column)
- run artifacts under `artifacts/`

---

## Installation

### Recommended: Docker

For the most reproducible GPU installation, use Docker.

```bash
# build Docker image for GPU runs
docker build -f Dockerfile.gpu -t hydromap:gpu .

# verify the container and GPU runtime:
docker run --rm -it --gpus all --user "$(id -u):$(id -g)" -v "$(pwd)":/workspace/HydroMap hydromap:gpu hydromap --help
```

If your cluster blocks Docker daemon access, run the same images with Apptainer/Singularity: [Apptainer user docs](https://apptainer.org/docs/user/latest/).

Docker quick start: [Docker Get Started](https://docs.docker.com/get-started/).

### Local GPU installation

This is the main non-Docker alternative. If Docker works for you, prefer the Docker path above.

First clone the repository and enter it, since the editable install step below (`pip install -e .`) installs HydroMap from your local checkout:

```bash
git clone https://github.com/samlobe/HydroMap.git
cd HydroMap
```

Then create the pinned local GPU environment. `environment.gpu.yml` defines the tested `Python 3.12` GPU stack:

```bash
conda install -y -n base -c conda-forge mamba
mamba env create -f environment.gpu.yml -n hydromap
conda activate hydromap
```

Then install HydroMap:
```bash
pip install -e .
```

Then build the water-triplet C extension:
```bash
cd hydromap/triplets
python setup.py build_ext --inplace
python test_waterlib_compilation.py
cd ..
```

GPU verification:
```bash
python -m openmm.testInstallation
python tools/verify_gpu_install.py

python - <<'PY'
import openmm, MDAnalysis, yaml
import importlib.util
print('openmm', openmm.__version__)
print('MDAnalysis', MDAnalysis.__version__)
print('yaml', yaml.__version__)
print('pdbfixer installed:', importlib.util.find_spec('pdbfixer') is not None)
PY

hydromap --help
```

If CUDA checks fail, fix the environment or driver mismatch rather than falling back silently.

If you must deviate from `environment.gpu.yml`, make sure your OpenMM build, CUDA target, and NVIDIA driver are compatible. See the [OpenMM installation guide](https://docs.openmm.org/development/userguide/application/01_getting_started.html) and NVIDIA's [CUDA compatibility documentation](https://docs.nvidia.com/deploy/cuda-compatibility/).

---

## Quick Start

Run full pipeline:

```bash
hydromap run --config configs/example_gpu.yaml
```

Docker equivalent (GPU):

```bash
docker run --rm -it --gpus all \
  --user "$(id -u):$(id -g)" \
  -v "$(pwd)":/workspace/HydroMap \
  -w /workspace/HydroMap \
  hydromap:gpu hydromap run --config configs/example_gpu.yaml
```

Run by stage:

```bash
hydromap prepare  --config configs/example_gpu.yaml
hydromap simulate --config configs/example_gpu.yaml
hydromap predict  --config configs/example_gpu.yaml
hydromap color    --config configs/example_gpu.yaml
hydromap analyze  --config configs/example_external_trajectory.yaml
```

## Config Essentials

For most users, the settings that matter are:

- run setup:
  `input_dir`, `proteins`/`protein`, `seeds`/`seed`
- simulation length:
  `md.nanoseconds`
- restraints:
  `md.restrain_selection`, `md.restraint_k`
- grouping:
  default is residue-by-residue within each chain; use `groups_file` for custom MDAnalysis selections
- coloring and filtering:
  `analysis.min_waters`, `analysis.color_properties`
- external-trajectory mode:
  `analysis.existing_processed_pdb`, `analysis.existing_trajectory`, optional `analysis.existing_topology`
  HydroMap uses the processed PDB and trajectory for analysis. In this mode, `protein` is used for output naming and `input_dir` is not used for structure loading.
- prep safeguards:
  `md.repair_missing_atoms`, `md.capping_mode`, `md.prep_policy`

Common advanced controls:

- sampling windows:
  `analysis.discard_initial_ns`, `analysis.tail_ns`
- frame sampling:
  `analysis.triplets_frame_stride`, `analysis.potentials_frame_stride`
- auditable/triplet-only analysis:
  `analysis.triplet_histogram_bin_width_deg`, `analysis.compute_potentials`
- ion controls:
  `md.neutralize`, `md.ionic_strength_molar`, `md.positive_ion`, `md.negative_ion`

For the full schema and advanced options, see [configs/README.md](configs/README.md).

---

## Analyze Existing Trajectories

If you already have an explicit-solvent PDB and matching trajectory, HydroMap can analyze them directly without running MD.

Set:

- `analysis.existing_processed_pdb`
- `analysis.existing_trajectory`
- optional `analysis.existing_topology`

Then run:

```bash
hydromap analyze --config <your_config>.yaml
```

If you also want predictions and colored PDBs:

```bash
hydromap run --stages analyze predict color --config <your_config>.yaml
```

Note:
- `Fdewet` predictions and PC values are only meaningful for simulations run with the `a99SBdisp` force field.
- `existing_topology` can be an OpenMM `System` `.xml` or a GROMACS `.top`
- if `existing_topology` is omitted, HydroMap will try rebuilding an OpenMM system if `forcefield: a99SBdisp`

## Notes
- Visualize your prepared system or trajectory in [ChimeraX](https://www.cgl.ucsf.edu/chimerax/download.html) or [PyMOL](https://pymol.org/) to confirm the solvation, ions, and overall geometry.
- If a protein fails during preparation, please leave a GitHub issue and include the input structure plus the HydroMap error output.
- In colored PDB outputs, gray atoms are usually ones HydroMap left at the pad B-factor (`-999` by default), meaning that group was not colored, most commonly because it had fewer hydration waters than the `min_waters` cutoff (typically 7, i.e. the low end of the training data) or never received a property assignment/selection match.
- General PTMs are not supported in the OpenMM preparation path.
- For prep behavior such as missing-atom repair, capping, strict/permissive prep policy, and `prepare_report.json`, see [configs/README.md](configs/README.md).
- For PC sign conventions and interpretation, see [configs/README.md](configs/README.md).

---

## Coloring Presets

Common value ranges:

- `Fdewet`: `4.0` to `6.5`
- `PC1`: `-8` to `8`
- `PC2`: `-2` to `8`
- `PC3`: `2` to `-2`

ChimeraX coloring commands:

- `Fdewet`: `color bfactor range 4.0,6.5 palette ^lipophilicity; color @@bfactor<-998 gray`
- `PC1`: `color bfactor range -8,8 palette red-white-blue; color @@bfactor<-998 gray`
- `PC2`: `color bfactor range -2,8 palette cyanmaroon; color @@bfactor<-998 gray`
- `PC3`: `color bfactor range -2,2 palette ^lipophilicity; color @@bfactor<-998 gray`

PyMOL coloring commands:

- `Fdewet`: `spectrum b, red_white_blue, minimum=4.0, maximum=6.5; color gray, b<-998`
- `PC1`: `spectrum b, red_white_blue, minimum=-8, maximum=8; color gray, b<-998`
- `PC2`: `spectrum b, cyan_maroon, minimum=-2, maximum=8; color gray, b<-998`
- `PC3`: `spectrum b, red_white_blue, minimum=-2, maximum=2; color gray, b<-998`

---

## Acknowledgements

[Shell Lab](https://theshelllab.org) and [Shea Group](https://labs.chem.ucsb.edu/shea/joan-emma/)  
[UCSB CNSI](https://www.cnsi.ucsb.edu/) for computing resources  
[Patel Group](https://patelgroup.seas.upenn.edu/) for [INDUS](https://github.com/patellab511/indus) used in model fitting  
[DE Shaw Group](https://www.deshaw.com/) for [a99SB-disp](https://github.com/paulrobustelli/Force-Fields/tree/master/Gromacs_FFs) force field
