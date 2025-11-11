#!/bin/bash

################## RESOURCES REQUEST (time, memory and cores)
#SBATCH -J tree_trimming
#SBATCH -t 240:00:00
#SBATCH --mem-per-cpu=9G
#SBATCH --cpus-per-task 20


################## OUTPUT LOG FILES
#SBATCH -o /work/%u/%x-%j.out   #  here they will dumped on your /work/<your_user> folder
#SBATCH -e /work/%u/%x-%j.err

################## MODULES

# Set the requested number of cores to the number of Threads your app should use

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}


################## COMMAND

cd /gpfs1/schlecker/home/kabiruna/Treemer

echo "start trimming"
date

singularity exec treemmer_sb/ python3 /Treemmer_v0.3.py gtdbtk_unrooted_bacteria.tree -X 2470 -lm CBB.csv -mc 1 -c ${SLURM_CPUS_PER_TASK:-1}


echo "end trimming"
date


