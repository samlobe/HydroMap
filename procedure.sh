#!/bin/bash
set -e
protein="example_protein"  
bash process_with_gromacs.sh ${protein}.pdb
python simulate_with_openmm.py ${protein}_processed.pdb topol.top -ns 5 -r
python run_triplets_parallel.py ${protein}.pdb traj.dcd --hydrationCutoff 5.5 --nprocs 30 --outdir triplets_55
python run_triplets_parallel.py ${protein}.pdb traj.dcd --hydrationCutoff 4.25 --nprocs 30 --outdir triplets_425
python run_potentials_parallel.py ${protein}.pdb traj.dcd --top topol.top --nprocs 30 --cutoff 5.5 --outdir potentials_55 --skip 100
python run_potentials_parallel.py ${protein}.pdb traj.dcd --top topol.top --nprocs 30 --cutoff 4.25 --outdir potentials_425 --skip 100 

python process_angles.py ${protein}.pdb --anglesDir triplets_55 -o triplet_data_${protein}_55.csv
python process_angles.py ${protein}.pdb --anglesDir triplets_425 -o triplet_data_${protein}_425.csv
python process_potentials.py ${protein}.pdb --energiesDir potentials_55 -o potentials_${protein}_55.csv
python process_potentials.py ${protein}.pdb --energiesDir potentials_425 -o potentials_${protein}_425.csv
python process_PCs.py ${protein}.pdb  --tripletFile triplet_data_${protein}_55.csv --PCs_outputFile PCs_${protein}_55.csv
python process_PCs.py ${protein}.pdb  --tripletFile triplet_data_${protein}_425.csv --PCs_outputFile PCs_${protein}_425.csv

