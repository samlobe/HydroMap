#!/bin/bash
#SBATCH --job-name=ProteinHydrophobicity
#SBATCH --output=ProteinHydrophobicity_%j.out
#SBATCH --error=ProteinHydrophobicity_%j.err
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=05:00:00
# #SBATCH --mail-user=<yourEmail> # uncomment these two lines and include email if desired
# #SBATCH --mail-type=END,FAIL    

# Check if a protein name was provided
if [ -z "$1" ]; then
    echo "Usage: $0 <protein_name>"
    exit 1
fi

PROTEIN_NAME="$1"

# Step 1: Process with GROMACS
bash process_with_gromacs.sh "${PROTEIN_NAME}.pdb"
if [ $? -ne 0 ]; then
    echo "Error processing with GROMACS. Exiting."
    exit 1
fi

# Step 2: Submit the simulation and capture the job ID
JOB_ID_1=$(python submit_simulation.py "${PROTEIN_NAME}_processed.gro")
if [ $? -ne 0 ]; then
    echo "Error submitting simulation. Exiting."
    exit 1
fi

# Step 3: Submit the triplets (with dependency on the previous job) and capture job IDs
JOB_ID_2_ARRAY=($(cd water_triplets && python submit_triplets.py "../${PROTEIN_NAME}.pdb" ../traj.dcd --dependency=afterany:$JOB_ID_1))
if [ $? -ne 0 ]; then
    echo "Error submitting triplet analysis. Exiting."
    exit 1
fi

# Wait until EACH job from Step 3 is finished
for job in "${JOB_ID_2_ARRAY[@]}"; do
    while squeue -j "$job" | grep -q "$job"; do
        sleep 60  # Check every 60 seconds
    done
done

# Step 4: Directly run the Python scripts without sbatch
python process_angles.py "../${PROTEIN_NAME}.pdb"
if [ $? -ne 0 ]; then
    echo "Error processing triplet angles. Exiting."
    exit 1
fi

python analyze_groups.py "../${PROTEIN_NAME}.pdb"
if [ $? -ne 0 ]; then
    echo "Error analyzing triplet angles or outputting pdbs. Exiting."
    exit 1
fi
