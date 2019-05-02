#!/bin/sh

#!/bin/bash -l
#SBATCH -J Extract_run
#SBATCH -o Extract_run-%j.out
#SBATCH -p Lewis
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 5G
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=oeow39@mail.missouri.edu
#SBATCH -t 2-00:00:00


input_dir="/storage/htc/bdm/tosin/GSDB/Data"


for dir in "${input_dir}"/GSE106*; do 
	[ -d "${dir}" ] && ((dir++))
	if [ -d "$dir" ] 
	then
		rootdir=$dir	
		rootname="${rootdir##*/}"    #Get the name
		echo "Root Directory Name: $rootname "	
		newdir="$rootdir/Processed_Data"
		mkdir $newdir
		gse="$rootdir/GSE*"
		mv $gse $newdir		
	fi
done