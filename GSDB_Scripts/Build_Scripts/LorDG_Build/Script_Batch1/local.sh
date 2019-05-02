#!/bin/sh
#!/bin/bash -l
#SBATCH -J LorDG_run
#SBATCH -o LorDG_run-%j.out
#SBATCH -p Lewis
#SBATCH -N 1
#SBATCH -n 6
#SBATCH --mem 80G
#SBATCH -t 2-00:00:00
#Generate a random parameter
n=6
rand=$(tr -cd '[:alnum:]' < /dev/urandom | head -c$n)
param="parameters_LorDG_${rand}.txt"
#Create algorithm parameters
{
echo "NUM = 1"
echo "OUTPUT_FOLDER = /storage/htc/bdm/tosin/GSDB/Structures/WT9059TG/GSM455138_30305AAXX.2/VC_1mb/LorDG/"
 echo "INPUT_FILE = /storage/htc/bdm/tosin/GSDB/Data/WT9059TG/Extracted_Data/GSM455138_30305AAXX.2/VC_1mb/chr9_list.txt"	
echo "VERBOSE = false"
echo "LEARNING_RATE = 0.07"
echo "MAX_ITERATION = 20000"
} > $param
# Call the LorDG alg_parameters
java -jar /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/LorDG.jar $param
rm $param	
#Change log_Name
NAME="chr9"
oldlogname="/storage/htc/bdm/tosin/GSDB/Structures/WT9059TG/GSM455138_30305AAXX.2/VC_1mb/LorDG/${NAME}_list_log.txt"
newlogname="/storage/htc/bdm/tosin/GSDB/Structures/WT9059TG/GSM455138_30305AAXX.2/VC_1mb/LorDG/${NAME}.log"
mv $oldlogname $newlogname
#Read all the pdb files in directory and print chromosome name only
algorithm_dir="/storage/htc/bdm/tosin/GSDB/Structures/WT9059TG/GSM455138_30305AAXX.2/VC_1mb/LorDG/chr*_*.pdb"
echo "Renaming Chromosomes..............."
for local_file in $algorithm_dir
do
	[ -f "${local_file}" ] && ((files++))
	echo ""
	echo "Processing $local_file file..."	
	echo ""
	NAME=`basename "$local_file"`
	tmp=$(echo "$NAME" | awk -F '_' '{print $1}' )
	newfname="/storage/htc/bdm/tosin/GSDB/Structures/WT9059TG/GSM455138_30305AAXX.2/VC_1mb/LorDG/${tmp}.pdb"
	echo " New file name = ${newfname}"
	#delete filename if exist
	if [ -f $newfname ] ; then
		echo "Deleting old file = $newfname "
		rm $newfname 
	fi
	mv $local_file $newfname
done
