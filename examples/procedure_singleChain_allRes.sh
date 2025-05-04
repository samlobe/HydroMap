#!/bin/bash
set -e
protein="example_protein"  # name of the protein (without .pdb)

# preprocess with gromcas
cp ${protein}.pdb simulation; cd simulation
bash process_with_gromacs.sh ${protein}.pdb
# run MD simulation
python simulate_with_openmm.py ${protein}_processed.pdb ${protein}.top -ns 5 ; cd ..
# analyze water triplets and water-protein potentials on 30 CPUs
python water_triplets/run_triplets_parallel.py ${protein}.pdb simulation/${protein}_traj.dcd --nprocs 30 --outdir water_triplets/angles_${protein} --skip 5
python water_potential/run_potentials_parallel.py ${protein}.pdb simulation/${protein}_traj.dcd --top simulation/${protein}.top --nprocs 30 --outdir water_potential/potentials_${protein} --skip 100
# output results for Fdewet, PC1, PC2, PC3, and per-water potential (processing angles/potentials, applying model, and outputting pdbs with bfactor columns)
python utils/process_and_predict.py ${protein}.pdb --anglesDir water_triplets/angles_${protein} --potentialsDir water_potential/potentials_${protein} --model models/Fdewet.joblib --outdir results/${protein} 
python utils/color_pdb_by_property.py ${protein}.pdb results/${protein}/${protein}_results.csv --outdir results/${protein}
