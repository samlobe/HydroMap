# protein_WaterStructure_HydrophobicityModelling
This is a tutorial explaining how to color a protein based on predicted dewetting free energy (using Sam Lobo & Saeed Najafi's hydrophobicity model) and water triplet distribution principal component contributions (from Dennis Robinson and Sally Jiao's work modelling small hydrophobes). This analysis takes <25 min of computation for a 100 residue protein, including ~20 min on a GPU (NVIDIA V100) and <3min on parallel CPUs (Intel Xeon Gold 6148 with 1 thread per residue / custom group). This can be reduced depending on your error bar tolerance.

Input : pdb file of just protein (no water or ions)
Output : 
1. pdb file with dewetting free energy predictions of each residue (or custom group) listed in the tempfactor column
2. 3 pdb files with the 3 principal component contributions of each residue (or custom group) listed in the tempfactor column
3. csv file with triplet distribution for each residue (or custom group), and many txt files of raw data triplet angles

## Procedure summary:

1. **Use GROMACS to build your system**: protein + water + ions
2. **Use OpenMM to run a short simulation**
3. **Measure water angles around each residue/group**
4. **Process water angles, and output pdbs with stored dewetting free energy predictions and principal component contributions**
5. **Color pdb by dewetting free energy and principal component contributions**

## [Click here for the complete tutorial](https://roamresearch.com/#/app/SamLobo/page/P2_MRPX_6) which includes detailed install instructions, the commands necessary for each step in the procedure, and some background info.
