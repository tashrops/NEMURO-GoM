#!/bin/bash
#SBATCH -p coaps14_q
#SBATCH -n 80
#SBATCH -t 120:00:00

# Load correct modules
module unload intel-openmpi
module unload intel-mvapich2
module load intel-openmpi

# Directory where the executable is located 
cd /gpfs/research/coaps/home/tas14j/expt/run

# Copy configuration to the experiment folder

# Submit job
srun -n 80 ./mitgcmuv >& log.txt


