#!/bin/bash
#SBATCH -p eoas_q
#SBATCH -n 80
#SBATCH -t 00:10:00

# Load correct modules
module unload intel-openmpi
module unload intel-mvapich2
module load intel-openmpi

# Directory where the executable is located 
cd /gpfs/home/tas14j/MITgcm_scp5/verification/sim_default_011221/run

# Submit job
srun -n 80 ./mitgcmuv >& log.txt


