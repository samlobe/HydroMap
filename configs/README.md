# HydroMap Configs

Files in this directory:

- `example_gpu.yaml`: recommended starting point for a standard GPU workflow
- `example_external_trajectory.yaml`: analyze an existing explicit-solvent trajectory without running MD
- `example_advanced.yaml`: copy-paste source for common advanced options

Suggested starting points:

- start from `example_gpu.yaml` for ordinary runs
- use `example_external_trajectory.yaml` when you already have a solvated PDB + trajectory
- copy blocks from `example_advanced.yaml` when you want to customize restraints, salt, histidines, or analysis windows

---

# HydroMap v2 Config Schema

Paths you set are resolved relative to the config file location. Omitted built-in
paths, such as the default model and output directory, resolve from the HydroMap
repository root.

## Top-Level

- `input_dir`: directory containing `<protein>.pdb` input files.
- `artifacts_root`: output root (default: repository `artifacts/`).
- `protein` or `proteins`: one protein or a list (without `.pdb`).
- `seed` or `seeds`: one seed or a list.
- `model_path`: prediction model path (default: packaged `Fdewet.joblib` model).
- `forcefield`: forcefield label for prediction processing (`a99SBdisp` default).
- `groups_file`: optional custom groups file (one MDAnalysis selection per line).

## MD

- `device`: use `gpu` for HydroMap runs.
- `nanoseconds` (default `0.5`)
- `restrain_selection` (alias: `restrain`), `restraint_k`
- `equilibration_ns`
- `equilibration_protocol`: `constant` (default) or staged `gradual` heating
- `timestep_ps`: integrator timestep (default `0.003`)
- `report_interval_ps`: trajectory/report interval (default `1.0`)
- `checkpoint_interval_ps`: equilibration/production checkpoint interval (default `10.0`)
- `constant_volume`: set `true` for NVT instead of the default NPT ensemble
- `initial_state`: optional OpenMM State XML; supports `{protein}` and `{seed}` placeholders and skips minimization/equilibration when no production checkpoint exists
- `strip_non_protein` (default `true`)
- `random_seed`, `preprocess_seed`, `velocity_seed`, `barostat_seed`
- `checkpoint_policy`: `error`, `resume`, or `overwrite`
- `allow_cpu_md`: advanced control; leave `false` for standard HydroMap runs
- `deterministic`, `cuda_precision`
- `repair_missing_atoms`: `none` (default) or `pdbfixer`
- `capping_mode`: `none` (default), `termini`, `breaks`, or `termini_and_breaks`
- `prep_policy`: `permissive` (default) or `strict`
- ion controls:
  - `neutralize` (default `true`)
  - `ionic_strength_molar` (default `0.0`)
  - `positive_ion` (`Na+`, `K+`, `Li+`, `Rb+`, `Cs+`)
  - `negative_ion` (`Cl-`, `Br-`, `F-`, `I-`)
- histidine controls:
  - `histidine_mode`: `auto`, `hid`, `hie`, or `hip`
  - `histidine_overrides`: per-residue overrides such as `"B:417": HIP` or `"505": HID`

## Analysis

- `triplets_device`: use `gpu`
- `triplets_frame_stride` (alias: `triplets_skip`)
- optional `triplets_sample_ps` (if set, overrides `triplets_frame_stride`)
- `triplets_hydration_cutoff`
- `triplet_histogram_bin_width_deg`: width for the auditable 40–180 degree histogram CSV (default `10` and must divide 140)
- `compute_potentials`: set `false` for triplet-only analysis; no topology is
  required, and prediction then remains unavailable
- `potentials_device`: use `gpu`
- `potentials_frame_stride` (alias: `potentials_skip`)
- optional `potentials_sample_ps` (if set, overrides `potentials_frame_stride`)
- `potentials_cutoff`
- time windows:
  - `discard_initial_ns`
  - `tail_ns` (default: all post-discard frames)
  - per-analysis overrides: `triplets_tail_ns`, `potentials_tail_ns`
  - aliases: `triplets_time_ns`, `potentials_time_ns`
- optional pre-existing simulation inputs (analyze without MD):
  - `existing_processed_pdb`
  - `existing_trajectory`
  - optional `existing_topology` (`.xml` OpenMM `System` or GROMACS `.top`)
  - a GROMACS `.top` is loaded through OpenMM; HydroMap does not call `gmx`
  - `existing_processed_pdb` must contain explicit solvent; HydroMap will reject solvent-free uploads because water structure/potential analysis would be undefined
  - if `existing_topology` is omitted, HydroMap will try to rebuild an OpenMM XML from the uploaded solvated PDB only when `forcefield: a99SBdisp`
  - external trajectory mode requires exactly one `protein` in the config
  - optional placeholders: `{protein}`, `{seed}`
- `min_waters`, `color_properties`

Every analysis run writes raw per-group angle files plus a complete histogram CSV in
the case `results` directory. The CSV records frame/angle totals and counts outside
the standard 40–180 degree histogram range so dropped values are visible.

## Resources

- `max_cpu_workers`
- `reserve_cpus`
- `max_gpu_jobs`

## Execution

- `profile`: `gpu_fast` or `balanced`

## Prediction Support

- HydroMap can compute water triplets and water-protein potentials from external explicit-solvent trajectories.
- `Fdewet`/PC prediction support is only validated for `a99SBdisp`.

## Histidine Example

```yaml
md:
  histidine_mode: auto
  histidine_overrides:
    "B:417": HIP
    "505": HID
```

## Histidine Behavior

- `histidine_mode: auto` uses the input structure to infer `HID`/`HIE`/`HIP` when possible.
- `histidine_mode: hip` forces all histidines positive unless a per-residue override says otherwise.
- Per-residue overrides take precedence over the global mode.
- Residue-only keys such as `"505"` are useful when chain IDs are unreliable in the raw PDB.

## Restart behavior

- `checkpoint_policy: resume` resumes production from the production checkpoint, or an interrupted equilibration from its separate equilibration checkpoint.
- `checkpoint_policy: overwrite` starts over and removes prior checkpoint/output files for that case.
- Production time and step count start at zero after equilibration, so resumed production is not offset by the equilibration duration.
- OpenMM checkpoints are platform- and version-sensitive. Use `initial_state` for a portable conditioned State XML; it is ignored when a production checkpoint is present.

## Preparation Behavior

- HydroMap drops incomplete terminal protein residues automatically during sanitization.
- If incomplete non-terminal residues are detected, HydroMap will raise a clear error by default; set `md.repair_missing_atoms: pdbfixer` to attempt broader missing-heavy-atom repair on a prep-only temporary copy.
- HydroMap detects peptide breaks before adding peptide bonds. Set `md.capping_mode` to cap true termini, break boundaries, or both with `ACE`/`NME`.
- The prepare stage writes `prepare_report.json`, a consolidated intake/prep report that references `input_sanitization.json`, `prep_audit.json`, and `topology_audit.json`.
- `md.prep_policy: permissive` records prep interventions in `prepare_report.json` and continues. `md.prep_policy: strict` writes the same report and then exits if preparation had to modify or repair the structure.

## PC Interpretation

- `+PC1`: icosahedral (~60 degrees) ordering shift (hydrophilic signature) at the expense of tetrahedral waters
- `+PC2`: enrichment near ~90 degrees (hydrophobic signature)
- `+PC3`: enrichment near ~50 degrees waters (hydrophilic signature)
