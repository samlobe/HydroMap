# protein_WaterStructure_HydrophobicityModel

Use water structure analysis to color your protein based on predicted dewetting free energy (using Sam Lobo & Saeed Najafi's hydrophobicity model) and water triplet distribution principal component contributions (from Dennis Robinson and Sally Jiao's work modelling small hydrophobes).  
  
This analysis takes <25 min of computation for a 100 residue protein, including ~20 min on a GPU (NVIDIA V100) and <3min on parallel CPUs (Intel Xeon Gold 6148 with 1 thread per solvated residue / custom group). This can be reduced depending on your error bar tolerance.  

---

### **Input:**
- pdb file of just protein (no water or ions)  

### **Output:**
- pdb file with dewetting free energy predictions of each residue (or custom group) listed in the tempfactor column
- 3 pdb files with the 3 principal component contributions of each residue (or custom group) listed in the tempfactor column
- csv file with triplet distribution for each residue (or custom group), and many txt files of raw data triplet angles  

### **Tutorial:**
[Click here](https://roamresearch.com/#/app/SamLobo/page/P2_MRPX_6) for the complete tutorial on Sam Lobo's Roam page which includes detailed install instructions, the commands necessary for each step in the procedure, some background info, and FAQs (e.g. "How can I get this working on my cluster?").

---

## Procedure summary:

1. Use GROMACS to ***build your system***: protein + water + ions
2. Use OpenMM to ***run a short MD simulation*** (GPU recommended)
3. ***Measure water angles*** around each residue/group (parallel CPUs recommended)
4. ***Process water angles***, and ***output pdbs*** with stored dewetting free energy predictions and principal component contributions
5. **Color pdbs** by dewetting free energy and principal component contributions

You can skip steps 1 & 2 if you already have a simulation to analyze.  

---

## Install requirements:

1. [Anaconda](https://www.anaconda.com/download) or miniconda for python tools
2. [GROMACS](https://manual.gromacs.org/documentation/current/install-guide/index.html) (Step 1)
3. [OpenMM](http://docs.openmm.org/latest/userguide/application/01_getting_started.html#installing-openmm) (Step 2) - you need to fix for TIP4P water constraints and you may need to install CUDA drivers depending on your architecture
4. Fortran compiler (Step 3) - anaconda's f2py should work (see tutorial's FAQs)
5. [MDAnalysis](https://www.mdanalysis.org/pages/installation_quick_start/) (Step 3 & 4)
6. [ChimeraX](https://www.cgl.ucsf.edu/chimerax/download.html) or [Pymol](https://pymol.org/2/) (recommended for Step 5)
See [tutorial](https://roamresearch.com/#/app/SamLobo/page/P2_MRPX_6) for more install instructions.

---

## Main Code:  
-
-
---

## Code specifically for UCSB's Pod Cluster (SLURM scheduler):  
- **submit_simulation.py**
  - Usage: `python submit_simulation.py <protein_processed[.gro]>`
  - Uses *simulate_with_openmm.py* to submit a GPU job (Step 2)
- **submit_triplets.py**
  - Usage: `python submit_triplets.py <protein[.pdb]> <trajectory[.dcd]`
  - Uses *triplet.py* and srun to parallelize and submit many CPU jobs, 1 per solvated residue/group (Step 3)
- **full_hydrophobicity_procedure.sh**
  - Usage: `sbatch full_hydrophobicity_procedure <protein_name>`
  - Manages and submits SLURM jobs (Steps 1-4)
  - Designed so that one command can fully process your protein & get colored pdbs
  - Calls *process_with_gromacs.sh*, *submit_simulation.sh*, *submit_triplets.py*, *process_angles.py*, and *analyze_groups.py*

