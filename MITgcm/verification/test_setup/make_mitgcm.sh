#!/bin/tcsh
# ---------- #
# Description: 
# Creates directories, files, and executable to run MITgcm 
# Written by: Taylor Shropshire
# Date: 6/6/16 
# --------- #

# Remove and Create Working Directories
/bin/rm -rf ./build ./run 
/bin/mkdir ./build ./run 

#load correct environment
module load intel-openmpi

# Compile MITgcm
cd ./build
../../../tools/genmake2 -mods=../code -enable=mnc -of=../../../tools/build_options/linux_amd64_ifort+mpi_generic >& build.log
make depend
make 

# Link files to run directory 
cd /gpfs/home/tas14j/MITgcm_scp5/verification/sim_default_011221/run
ln -s /gpfs/home/tas14j/MITgcm_scp5/verification/sim_default_011221/input/* ./
ln -s /gpfs/home/tas14j/MITgcm_scp5/verification/sim_default_011221/build/mitgcmuv ./ 
ln -s /gpfs/home/tas14j/MITgcm_scp5/verification/sim_default_011221/fields/grid/depth_gom.bin ./
ln -s /gpfs/home/tas14j/MITgcm_scp5/verification/sim_default_011221/fields/grid/grid.face001.bin ./
ln -s /gpfs/home/tas14j/MITgcm_scp5/verification/sim_default_011221/fields/ini_cond/*field* ./
ln -s /gpfs/home/tas14j/MITgcm_scp5/verification/sim_default_011221/fields/force_fields/*SurfRad* ./
ln -s /gpfs/home/tas14j/MITgcm_scp5/verification/sim_default_011221/fields/force_fields/*Riv* ./
ln -s /gpfs/home/tas14j/MITgcm_scp5/verification/sim_default_011221/fields/OB_cond/*OB* ./
ln -s /gpfs/home/tas14j/MITgcm_scp5/verification/sim_default_011221/fields/off_fields_link ./


