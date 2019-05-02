#!/bin/sh

#!/bin/bash -l
#SBATCH -J Run_Program_run
#SBATCH -o Run_Program_run-%j.out
#SBATCH -p Interactive
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 1G
#SBATCH --mail-type=end
#SBATCH -t 0-03:00:00

# Determine the number folders/Directory
base_dir="/storage/htc/bdm/tosin/GSDB_Scripts/scripts/SIMBA3D"


for local_file in $base_dir/create*
do				
	[ -f "${local_file}" ] && ((files++))
	
	chmod +x $local_file
	$local_file
	
done
				