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

# Organize input files and the MITgcm exe. created above
cd /gpfs/research/coaps/home/tas14j/final_w_biomass_082619/run
ln -s /gpfs/home/tas14j/MITgcm_scp4/verification/NEMURO_GoM/input/* ./
ln -s /gpfs/home/tas14j/MITgcm_scp4/verification/NEMURO_GoM/build/mitgcmuv ./ 
ln -s /gpfs/research/coaps/home/tas14j/input_version/one_grid_point/depth_gom.bin ./
ln -s /gpfs/research/coaps/home/tas14j/input_version/one_grid_point/grid.face001.bin ./

ln -s /gpfs/research/coaps/home/tas14j/input_version/one_grid_point/flt_positions_071619/* ./
ln -s /gpfs/research/coaps/home/tas14j/input_version/one_grid_point/inital_and_OB_conditons_071819/*field* ./
ln -s /gpfs/research/coaps/home/tas14j/input_version/one_grid_point/inital_and_OB_conditons_022619/*light* ./
ln -s /gpfs/research/coaps/home/tas14j/input_version/one_grid_point/inital_and_OB_conditons_071819/*riv* ./
ln -s /gpfs/research/coaps/home/tas14j/input_version/one_grid_point/inital_and_OB_conditons_120618/*OB* ./
ln -s /gpfs/research/coaps/home/tas14j/input_version/one_grid_point/input_off_link ./

# To Run MITgcm, 
#> sbatch mitgcm_mpi_script.sh


