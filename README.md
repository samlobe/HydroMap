# hydrophobicityPod
This is a tutorial explaining how to color a pdb file based on predicted dewetting free energy using Sam Lobo & Saeed Najafi's hydrophobicity model.
Using Sam Lobo &amp; Saeed Najafi's hydrophobicity model on UCSB's pod cluster.

Input = pdb file of just protein (no water or ions)
Output = pdb file with dewetting free energies of each residue listed in the beta column

## Procedure summary:

1. **Use GROMACS to build your system**: protein + water + ions
2. **Use OpenMM to simulate your system for 10ns**.
3. **Measure water angles around each residue**.
4. **Predict dewetting free energy from water angles and output new pdb**.
5. **Color pdb by dewetting free energy**.

## Installation requirements:

- **GROMACS (Step 1)**:
  - Load on Pod cluster:
  - Install on personal computer (not required):

- **OpenMM (Step 2)**:
  - Install on Pod cluster:
  - Install on personal computer (not required):

- **Fortran compiled code (Step 3)**:
  - Compile fortran library (waterlib.f90) on Pod cluster:
  - Compile fortran library (waterlib.f90) on personal computer (not required):

- **MDAnalysis for atom selection (Step 3)**:
  - Install on Pod cluster:
  - Install on personal computer (recommended):

- **ChimeraX (Step 5)**:
  - Don't install on Pod cluster.
  - Install on your personal computer:

> The `hydrophobicity_tutorial` folder is on Pod at `/home/lobo/hydrophobicity_tutorial`

## Contents of hydrophobicity_tutorial:
[List of contents you provided, formatted as a list or sub-sections]

## Full Procedure:

### STEP 1) Use GROMACS to build your system:
[Detailed steps you provided, formatted with bullet points or ordered lists]

### STEP 2) Use OpenMM to simulate your system for 10ns:
[Detailed steps...]

### STEP 3) Measure water angles around each residue:
[Detailed steps...]

### STEP 4) Predict dewetting free energy from water angles and output new pdb:
[Detailed steps...]

### STEP 5) Color pdb by dewetting free energy:
[Detailed steps...]


