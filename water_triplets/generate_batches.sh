#!/bin/bash

#SBATCH -J residue_analysis       # Job name
#SBATCH -o outLog                 # Name of stdout output file
#SBATCH -e errLog                 # Name of stderr error file
#SBATCH -p batch                  # Queue (partition) name
#SBATCH -N 1                      # Total # of nodes (must be 1 for serial)
#SBATCH --cpus-per-task=1
#SBATCH -t 00:04:00               # Run time (hh:mm:ss)
#SBATCH --mail-user=lobo@ucsb.edu
#SBATCH --mail-type=all           # Send email at begin and end of job

# Other commands must follow all #SBATCH directives...

module list
pwd
date

start=323
end=502

batch_count=5
total_residues=$((end - start + 1))
batch_size=$((total_residues / batch_count))  # Calculate the size of each batch
remainder=$((total_residues % batch_count))  # Calculate the remainder residues

for ((i = 0; i < batch_count; i++))
do
  start_residue=$((start + i * batch_size))
  end_residue=$((start_residue + batch_size - 1))

  # Add one additional residue to the last batch to accommodate the remainder
  if ((i == batch_count - 1)); then
    end_residue=$((end_residue + remainder))
  fi

  cat > submit_batch$i.sh <<EOL
#!/bin/bash

#SBATCH -J covid_triplet           # Job name
#SBATCH -o outLog       # Name of stdout output file
#SBATCH -e errLog       # Name of stderr error file
#SBATCH -p batch          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 40              # Total # of mpi tasks (should be 1 for serial)
#SBATCH --cpus-per-task=1
#SBATCH -t 00:40:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=wwebb@ucsb.edu
#SBATCH --mail-type=all    # Send email at begin and end of job

# Other commands must follow all #SBATCH directives...

module list
pwd
date

for residue in {$start_residue..$end_residue}; do
    srun --ntasks 1 --exclusive -c 1 python residue_triplet.py \$residue &
done

wait

EOL
done

