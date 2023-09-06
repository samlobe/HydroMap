import argparse
import subprocess

def minutes_to_hms(minutes):
    """Converts minutes to hh:mm:ss format."""
    h = minutes // 60
    m = minutes % 60
    return f"{h:02d}:{m:02d}:00"

def create_slurm_script(protein_name, time_limit):
    time_limit_str = minutes_to_hms(time_limit)
    slurm_script = f"""#!/bin/bash

#SBATCH --nodes=1 --ntasks-per-node=1 --partition=gpu # Pod cluster's GPU queue
#SBATCH --gres=gpu:1
#SBATCH --time={time_limit_str}
#SBATCH --job-name={protein_name}
# #SBATCH --mail-user=<yourEmail> # uncomment these two lines and include email if desired
# #SBATCH --mail-type=END,FAIL    # Send email at begin and end of job

cd $SLURM_SUBMIT_DIR
module load cuda/11.2
conda activate hydrophobicity

srun --gres=gpu:1 python simulate_with_openmm.py {protein_name}
"""

    with open("submit_simulation.sh", "w") as file:
        file.write(slurm_script)

def main():
    parser = argparse.ArgumentParser(description='Create and submit a SLURM job for protein simulation.')
    parser.add_argument('-p','--protein_name', required=True, help='Name of the protein for the simulation job.')
    parser.add_argument('-t','--time_limit', type=int, default=40, help='Time limit for the SLURM job in minutes. Default is 40 min.')
    
    args = parser.parse_args()

    create_slurm_script(args.protein_name, args.time_limit)
    
    subprocess.run(["sbatch", "submit_simulation.sh"])

if __name__ == "__main__":
    main()

