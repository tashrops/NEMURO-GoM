#!/bin/bash
# ------------------------------------------------ #
# Purpose: Creates soft links and renames flow field files for running MITgcm 
# Written by: Taylor Shropshire
# Date: 06/29/17
# ------------------------------------------------- #

# Parameters
dir1="/gpfs/home/tas14j/MITgcm_scp5/verification/sim_default_011221/fields/off_fields"
dir2="/gpfs/home/tas14j/MITgcm_scp5/verification/sim_default_011221/fields/off_fields_link"
field_freq=86400
physical_dt=1800

# Make directory for soft links
/bin/rm -rf $dir2
/bin/mkdir $dir2
cd $dir2
/bin/mkdir tmp sal uvel vvel wvel

mit_count1=0000000000

# Create soft links
ln -s ${dir1}/tmp/* ${dir2}/tmp
ln -s ${dir1}/sal/* ${dir2}/sal
ln -s ${dir1}/uvel/* ${dir2}/uvel
ln -s ${dir1}/vvel/* ${dir2}/vvel
ln -s ${dir1}/wvel/* ${dir2}/wvel 

# Rename soft links for MITgcm file nameing
num_files=`find "${dir1}/sal" -name "*.data" | wc -l`
file_count1=0001
file_inc=1
let mit_inc=$field_freq/$physical_dt

for ii in $(seq 1 $num_files); do
file_count2=`echo $file_count1 | awk '{printf("%04d", $1)}'`
mit_count2=`echo $mit_count1 | awk '{printf("%010d", $1)}'`

in_u="${dir2}/uvel/u.${file_count2}.data"
in_v="${dir2}/vvel/v.${file_count2}.data"
in_w="${dir2}/wvel/w.${file_count2}.data"
in_t="${dir2}/tmp/t.${file_count2}.data"
in_s="${dir2}/sal/s.${file_count2}.data"

out_u="${dir2}/uvel/U.${mit_count2}.data"
out_v="${dir2}/vvel/V.${mit_count2}.data"
out_w="${dir2}/wvel/W.${mit_count2}.data"
out_t="${dir2}/tmp/T.${mit_count2}.data"
out_s="${dir2}/sal/S.${mit_count2}.data"

/bin/mv $in_u $out_u
/bin/mv $in_v $out_v
/bin/mv $in_w $out_w
/bin/mv $in_t $out_t
/bin/mv $in_s $out_s

let file_count1=$file_count1+$file_inc
let mit_count1=$mit_count1+$mit_inc
echo "File Number = "$file_count2", MITgcm Number = "$mit_count2

done 	# for ii in $(seq 1 $num_files); do

/bin/ls ${dir2}/uvel/ | wc -l
echo "Done ${dir1}" 
echo "Finished Renaming Flow Field File Names"






