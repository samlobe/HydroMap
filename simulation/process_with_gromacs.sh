#!/bin/bash

# Ensure the protein name is passed to the script
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 protein_name"
    exit 1
fi

protein=$1

# If the protein name ends with .pdb, trim it
if [[ $protein == *.pdb ]]; then
    protein=${protein%.pdb}
fi

# Check for the existence of the .pdb file
if [ ! -f "$protein.pdb" ]; then
    echo "$protein.pdb not found!"
    exit 1
fi

# create the log file and delete the old one if it exists
log_file="${protein}_gromacs_processing.log"
if [ -f "$log_file" ]; then
    rm "$log_file"
fi
touch "$log_file"

topol_file="${protein}.top"

# Function to run gmx command and display messages
run_gmx_command() {
    command=$1
    message=$2

    echo -e "\n$message"
    echo "Executing: $command"

    {
        echo -e "\n=== [$protein] $message ==="
        echo "Command: $command"
    } >> "$log_file"

    if [[ "$command" == *"|"* ]]; then
        eval "$command" >> "$log_file" 2>&1
    else
        $command >> "$log_file" 2>&1
    fi
}

echo "Processing with GROMACS..."
echo "(writing progress/errors to $log_file)"

# pdb2gmx
run_gmx_command "echo -e \"1\n1\" | gmx pdb2gmx -f ${protein}.pdb -o ${protein}_noSolvent.gro -ignh -p ${topol_file}" \
    "Creating topology file (topol.top) and GROMACS coordinate file (.gro) from the pdb file using a force field. If using a custom force field file, place the force field folder in the working directory."

# editconf
run_gmx_command "gmx editconf -f ${protein}_noSolvent.gro -o ${protein}_newbox.gro -c -d 1.0 -bt triclinic" \
    "Setting the size of the simulation box..."

# solvate
run_gmx_command "gmx solvate -cp ${protein}_newbox.gro -cs tip4p -o ${protein}_solv.gro -p ${topol_file}" \
    "Adding water molecules to the simulation box..."

# add_ions (two commands)
run_gmx_command "gmx grompp -f ions.mdp -c ${protein}_solv.gro -p ${topol_file} -o ions.tpr" \
    "Adding ions to the simulation box. First reading the instructions from ions.mdp and processing all the interactions..."
run_gmx_command "echo 13 | gmx genion -s ions.tpr -o ${protein}_processed.gro -p ${topol_file} -pname NA -nname CL -neutral" \
    "Now we add in the ions... Selecting SOL to replace the SOL (i.e. water) atoms with ions."

# convert from .gro back to .pdb
run_gmx_command "gmx grompp -f ions.mdp -c ${protein}_processed.gro -p ${topol_file} -o ${protein}_processed.tpr" \
    "Getting a .tpr file so we can convert the processed file from .gro back to .pdb..."
run_gmx_command "echo 0 | gmx trjconv -f ${protein}_processed.gro -s ${protein}_processed.tpr -o ${protein}_processed.pdb" \
    "Converting the processed .gro file back to .pdb format..."

# copy the {protein}_processed.pdb to the above directory
cp ${protein}_processed.pdb ../

echo "Cleaning up files..."
rm ${protein}_newbox.gro ${protein}_noSolvent.gro ${protein}_solv.gro posre*.itp mdout.mdp ions.tpr ${protein}_processed.tpr ${protein}_processed.gro
find . -maxdepth 1 -type f -name "#*" -exec rm -f {} \;

echo -e "\nDone processing with GROMACS. Use ${protein}_processed.pdb and topol.top for the simulation.\n"
echo "To run the simulation, use 'python simulate_with_openmm.py ${protein}_processed.pdb topol.top'"