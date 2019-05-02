#!/bin/sh

# Determine the number folders/Directory
base_dir="/storage/htc/bdm/tosin/GSDB_Scripts/Build_Scripts/LorDG_Build/Script_Batch2"


for local_file in $base_dir/Build*
do				
	[ -f "${local_file}" ] && ((files++))
	
	sbatch $local_file
	
done
				