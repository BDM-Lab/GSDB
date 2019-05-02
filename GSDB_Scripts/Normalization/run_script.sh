#!/bin/sh

#!/bin/bash -l
#SBATCH -J Run_Normalizer_run
#SBATCH -o Run_Normalizer_run-%j.out
#SBATCH -p Interactive
#SBATCH -N 1
#SBATCH -n 3
#SBATCH --mem 10G
#SBATCH --mail-type=end
#SBATCH -t 0-02:00:00

# Determine the number folders/Directory
base_dir="/storage/htc/bdm/tosin/GSDB_Scripts/Normalization"


for local_file in $base_dir/40kb_Normalizer_*
do				
	[ -f "${local_file}" ] && ((files++))
	
	sbatch $local_file
	
done
				