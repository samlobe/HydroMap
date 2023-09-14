#!/bin/bash

#SBATCH --job-name=Protein_Hydrophobicity
#SBATCH --output=Hydrophobicity_%j.out
#SBATCH --error=Hydrophobicity_%j.err
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=05:00:00
# #SBATCH --mail-user=<your_email> # uncomment these two lines and include email if desired
# #SBATCH --mail-type=ALL

# Function to display section headers
section() {
    echo -e "\n\n################################"
    echo "# $1"
    echo "################################\n"
}

# Check if a protein name was provided
if [ -z "$1" ]; then
    echo "Usage: $0 <protein_name>"
    exit 1
fi

PROTEIN_NAME="$1"
echo "ANALYZING ${PROTEIN_NAME}\n\n"
# Change to the submit directory
cd $SLURM_SUBMIT_DIR
source ~/.bashrc # Pod cluster uses bash shell
echo "Activating conda environment: 'hydrophobicity'"
conda activate hydrophobicity

# Step 1: Process with GROMACS
section "STEP 1: Build System with GROMACS"
bash process_with_gromacs.sh "${PROTEIN_NAME}.pdb"
if [ $? -ne 0 ]; then
    echo "Error in STEP 1. Exiting."
    exit 1
fi

# Step 2: Submit the simulation and capture the job ID
section "STEP 2: Run Short MD Simulation"
OUTPUT_STRING=$(python submit_simulation.py "${PROTEIN_NAME}_processed.gro")
echo "Job info from STEP 2: $OUTPUT_STRING"
if [ $? -ne 0 ]; then
    echo "Error in STEP 2. Exiting."
    exit 1
fi
JOB_ID_1=$(echo $OUTPUT_STRING | grep -oP '\d+')

# Step 3: Submit the triplets
section "STEP 3: Measure Water Triplets"
cd water_triplets
OUTPUT_STRING_2=$(python submit_triplets.py "../${PROTEIN_NAME}.pdb" ../traj.dcd --dependency=afterany:$JOB_ID_1)
echo "Output from STEP 3: $OUTPUT_STRING_2"
if [ $? -ne 0 ]; then
    echo "Error in STEP 3. Exiting."
    exit 1
fi
JOB_ID_2_ARRAY=($(echo $OUTPUT_STRING_2 | grep -oP '\d+'))
echo "Triplet submission job IDs from STEP 3: ${JOB_ID_2_ARRAY[*]}"

# Validate job IDs
echo "Validating job IDs from STEP 3..."
for job in "${JOB_ID_2_ARRAY[@]}"; do
    if ! [[ $job =~ ^[0-9]+$ ]]; then
        echo "Error: Invalid job ID '$job'. Exiting."
        exit 1
    fi
done

# Wait for jobs
echo "Waiting for jobs from STEP 3 to complete..."
for job in "${JOB_ID_2_ARRAY[@]}"; do
    while squeue -j "$job" | grep -q "$job"; do
        sleep 60
    done
done
echo "Jobs from STEP 3 completed!"

# Step 4: Run Python scripts
section "STEP 4: Process Triplet Angles & Output PDBs"
python process_angles.py "../${PROTEIN_NAME}.pdb"
if [ $? -ne 0 ]; then
    echo "Error in STEP 4a (processing angles). Exiting."
    exit 1
fi

python analyze_groups.py "../${PROTEIN_NAME}.pdb"
if [ $? -ne 0 ]; then
    echo "Error in STEP 4b (analyzing angles & outputting PDBs). Exiting."
    exit 1
fi

# Completion message
section "COMPLETION"
echo "Completed all 4 steps of the water structure-hydrophobicity procedure!"
echo "Open the pdbs and color them in ChimeraX or Pymol. Check out the plots."
