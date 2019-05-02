#!/bin/sh

#!/bin/bash -l
#SBATCH -J Run_Program_run
#SBATCH -o Run_Program_run-%j.out
#SBATCH -p Interactive
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem 2G
#SBATCH --mail-type=end

# Determine the number folders/Directory
base_dir="/storage/htc/bdm/tosin/GSDB_Scripts/scripts"


for local_file in $base_dir/create_ChromSDE_*
do				
	[ -f "${local_file}" ] && ((files++))
	
	sbatch $local_file
	
done
				