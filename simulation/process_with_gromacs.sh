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

# Function to run gmx command and display messages
# Function to run gmx command and display messages
run_gmx_command() {
    command=$1
    message=$2

    echo -e "\n$message"
    echo "Executing: $command"

    # Check if the command has a pipe
    if [[ "$command" == *"|"* ]]; then
        # Using eval to execute piped commands
        eval "$command >> gromacs_processing.log 2>&1"
    else
        # Execute the command normally
        $command >> gromacs_processing.log 2>&1
    fi

    if [ $? -ne 0 ]; then
        echo "Error executing: $command"
        exit 1
    fi  
}

echo "Processing with GROMACS..."
echo "(writing progress/errors to gromacs_processing.log)"

# pdb2gmx
run_gmx_command "echo -e \"1\n1\" | gmx pdb2gmx -f ${protein}.pdb -o ${protein}_noSolvent.gro -ignh" \
    "Creating topology file (topol.top) and GROMACS coordinate file (.gro) from the pdb file using a force field. If using a custom force field file, place the force field folder in the working directory."

# Remove the old log file if it exists and create a new empty one
[ -f gromacs_processing.log ] && rm gromacs_processing.log
touch gromacs_processing.log

# editconf
run_gmx_command "gmx editconf -f ${protein}_noSolvent.gro -o ${protein}_newbox.gro -c -d 1.0 -bt triclinic" \
    "Setting the size of the simulation box..."

# solvate
run_gmx_command "gmx solvate -cp ${protein}_newbox.gro -cs tip4p -o ${protein}_solv.gro -p topol.top" \
    "Adding water molecules to the simulation box..."

# add_ions (two commands)
run_gmx_command "gmx grompp -f ions.mdp -c ${protein}_solv.gro -p topol.top -o ions.tpr" \
    "Adding ions to the simulation box. First reading the instructions from ions.mdp and processing all the interactions..."
run_gmx_command "echo 13 | gmx genion -s ions.tpr -o ${protein}_processed.gro -p topol.top -pname NA -nname CL -neutral" \
    "Now we add in the ions... Selecting SOL to replace the SOL (i.e. water) atoms with ions."

# convert from .gro back to .pdb
run_gmx_command "gmx grompp -f ions.mdp -c ${protein}_processed.gro -p topol.top -o ${protein}_processed.tpr" \
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