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

## [Click here for the complete tutorial](https://roamresearch.com/#/app/SamLobo/page/P2_MRPX_6) which includes detailed install instructions, the commands necessary for each step in the procedure, and some background info.
