#!/bin/sh

#!/bin/bash -l
#SBATCH -J Run_Program_run
#SBATCH -o Run_Program_run-%j.out
#SBATCH -p Interactive
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem 2G
#SBATCH --mail-type=end
#SBATCH -t 0-02:00:00

# Determine the number folders/Directory
base_dir="/storage/htc/bdm/tosin/GSDB_Scripts/Build_Scripts/3DMax_Build/Scripts_Batch1"


for local_file in $base_dir/Build*
do				
	[ -f "${local_file}" ] && ((files++))
	
	chmod +x $local_file
	$local_file
	
done
				