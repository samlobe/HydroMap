#!/bin/bash
set -e
protein="example_protein"  
bash process_with_gromacs.sh ${protein}.pdb
python simulate_with_openmm.py ${protein}_processed.pdb topol.top -ns 5 -r
python water_triplets/run_triplets_parallel.py ${protein}.pdb traj.dcd --nprocs 30 --outdir angles --skip 5
python water_potential/run_potentials_parallel.py ${protein}.pdb traj.dcd --top topol.top --nprocs 30 --outdir potentials --skip 100 
python process_and_predict.py ${protein}.pdb --anglesDir angles --potentialsDir potentials --model ../models/Fdewet.joblib --outdir results
python color_pdb_by_property.py ${protein}.pdb results/${protein}_results.csv --outdir results # outputs pdbs with Fdewet_pred, PC1, PC2, PC3, and per-water potential in bfactor columns
