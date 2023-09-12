# hydrophobicityPod
This is a tutorial explaining how to color a pdb file based on predicted dewetting free energy using Sam Lobo & Saeed Najafi's hydrophobicity model.  
Using Sam Lobo &amp; Saeed Najafi's hydrophobicity model on UCSB's pod cluster.  

Input = pdb file of just protein (no water or ions)  
Output = pdb file with dewetting free energies of each residue listed in the beta column  

## Procedure summary:

1. **Use GROMACS to build your system**: protein + water + ions
2. **Use OpenMM to run a short simulation**.
3. **Measure water angles around each residue/group**.
4. **Process water angles, and output pdbs with stored dewetting free energy predictions and principal component contributions**.
5. **Color pdb by dewetting free energy and principal component contributions.**.

## [Click here for the complete tutorial](https://roamresearch.com/#/app/SamLobo/page/P2_MRPX_6) which includes detailed install instructions, the commands necessary for each step in the procedure, and some background info.
