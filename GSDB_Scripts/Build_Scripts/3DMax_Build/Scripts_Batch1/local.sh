#!/bin/bash -l
#SBATCH -J 3DMax_run
#SBATCH -o 3DMax_run-%j.out
#SBATCH -p Lewis
#SBATCH -N 1
#SBATCH -n 6
#SBATCH --mem 80G
#SBATCH -t 2-00:00:00
#Generate a random parameter
n=6
rand=$(tr -cd '[:alnum:]' < /dev/urandom | head -c$n)
param="parameters_3DMax_${rand}.txt"
#Create algorithm parameters
{
echo "NUM = 1"
echo "OUTPUT_FOLDER = /storage/htc/bdm/tosin/GSDB/Structures/ENCSR662QKG_VH9178CA/GSE92819_ENCFF876LAW/VC_5kb/3DMax/"
 echo "INPUT_FILE = /storage/htc/bdm/tosin/GSDB/Data/ENCSR662QKG_VH9178CA/Extracted_Data/GSE92819_ENCFF876LAW/VC_5kb/chrY_list.txt"	
echo "VERBOSE = false"
echo "LEARNING_RATE = 1.0"
echo "MAX_ITERATION = 20000"
} > $param
# Call the 3DMax alg_parameters
java -jar /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/3DMax.jar $param
rm $param	
#Change log_Name
NAME="chrY"
oldlogname="/storage/htc/bdm/tosin/GSDB/Structures/ENCSR662QKG_VH9178CA/GSE92819_ENCFF876LAW/VC_5kb/3DMax/${NAME}_list_log.txt"
newlogname="/storage/htc/bdm/tosin/GSDB/Structures/ENCSR662QKG_VH9178CA/GSE92819_ENCFF876LAW/VC_5kb/3DMax/${NAME}.log"
mv $oldlogname $newlogname
#Read all the pdb files in directory and print chromosome name only
algorithm_dir="/storage/htc/bdm/tosin/GSDB/Structures/ENCSR662QKG_VH9178CA/GSE92819_ENCFF876LAW/VC_5kb/3DMax/chr*_*.pdb"
echo "Renaming Chromosomes..............."
for local_file in $algorithm_dir
do
	[ -f "${local_file}" ] && ((files++))
	echo ""
	echo "Processing $local_file file..."	
	echo ""
	NAME=`basename "$local_file"`
	tmp=$(echo "$NAME" | awk -F '_' '{print $1}' )
	newfname="/storage/htc/bdm/tosin/GSDB/Structures/ENCSR662QKG_VH9178CA/GSE92819_ENCFF876LAW/VC_5kb/3DMax/${tmp}.pdb"
	echo " New file name = ${newfname}"
	#delete filename if exist
	if [ -f $newfname ] ; then
		echo "Deleting old file = $newfname "
		rm $newfname 
	fi
	mv $local_file $newfname
done
