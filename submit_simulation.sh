#!/bin/bash

#SBATCH --nodes=1 --ntasks-per-node=1 --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=00:40:00
#SBATCH --job-name=humanPDZ3
#SBATCH --mail-user=lobo@ucsb.edu
#SBATCH --mail-type=all    # Send email at begin and end of job

cd $SLURM_SUBMIT_DIR
module load cuda/11.2
conda activate openmm

srun --gres=gpu:1 python simulate_with_openmm.py HSPDZ3_capped_processed

