#!/bin/bash -l
#!/bin/bash -l
#SBATCH -J MOGEN_run
#SBATCH -o MOGEN_run-%j.out
#SBATCH -p Lewis
#SBATCH -N 1
#SBATCH -n 6
#SBATCH --mem 80G
#SBATCH -t 2-00:00:00
#Generate a random parameter
n=6
rand=$(tr -cd '[:alnum:]' < /dev/urandom | head -c$n)
param="parameters_MOGEN_${rand}.txt"
#Create algorithm parameters
{
 echo "NUM = 1"
 echo "NBR_OF_CHR = 1"
 echo "INTRA_IF_THRESHOLD = 80%"
 echo "ADJACENT_DIST = 1.5"
 echo "CONTACT_DIST = 6.0"
 echo "POS_MIN_DIST = 0.2"
 echo "NEG_MAX_DIST_INTRA = 30"
 echo "POS_MAX_DIST_WEIGHT_FILE = /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/parameters/pos_max_dist_weight_1_normal.txt"
 echo "POS_MIN_DIST_WEIGHT_FILE = /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/parameters/pos_min_dist_weight_1_normal.txt"
 echo "NEG_MIN_DIST_WEIGHT_FILE = /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/parameters/neg_min_dist_weight_1_normal.txt"
 echo "NEG_MAX_DIST_WEIGHT_FILE = /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/parameters/neg_max_dist_weight_1_normal.txt"
 echo ""
 echo "OUTPUT_FOLDER = /storage/htc/bdm/tosin/GSDB/Structures/ENCSR662QKG_VH9178CA/GSE92819_ENCFF876LAW/VC_5kb/MOGEN/"
 echo "INPUT_FILE = /storage/htc/bdm/tosin/GSDB/Data/ENCSR662QKG_VH9178CA/Extracted_Data/GSE92819_ENCFF876LAW/VC_5kb/chrY_list.txt"
 echo "VERBOSE = false"	
 echo "LEARNING_RATE = 0.001"
 echo "MAX_ITERATION = 200000"
} > $param
# Call the LorDG alg_parameters
java -jar /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/MOGEN.jar $param
echo "construction completed... "
rm $param	
#Read all the pdb files in directory and print chromosome name only
algorithm_dir="/storage/htc/bdm/tosin/GSDB/Structures/ENCSR662QKG_VH9178CA/GSE92819_ENCFF876LAW/VC_5kb/MOGEN/chr*_*.pdb"
echo ".........................................."
echo ".......Renaming Chromosomes..............."
echo ".........................................."
for local_file in $algorithm_dir
do
	[ -f "${local_file}" ] && ((files++))
	echo ""
	echo "Processing $local_file file..."	
	echo ""
	NAME=`basename "$local_file"`
	tmp=$(echo "$NAME" | awk -F '_' '{print $1}' )
	newfname="/storage/htc/bdm/tosin/GSDB/Structures/ENCSR662QKG_VH9178CA/GSE92819_ENCFF876LAW/VC_5kb/MOGEN/${tmp}.pdb"
	echo " New file name = ${newfname}"
	#delete filename if exist
	if [ -f $newfname ] ; then
		echo "Deleting old file = $newfname "
		rm $newfname 
	fi
	mv $local_file $newfname
 MatrixNAME="/storage/htc/bdm/tosin/GSDB/Data/ENCSR662QKG_VH9178CA/Extracted_Data/GSE92819_ENCFF876LAW/VC_5kb/${tmp}_matrix.txt"
 echo 
 # Analysis	
 java -jar /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/pdbAnalyze.jar $MatrixNAME $newfname /storage/htc/bdm/tosin/GSDB/Structures/ENCSR662QKG_VH9178CA/GSE92819_ENCFF876LAW/VC_5kb/MOGEN
	#Change log_Name 
	oldlogname="/storage/htc/bdm/tosin/GSDB/Structures/ENCSR662QKG_VH9178CA/GSE92819_ENCFF876LAW/VC_5kb/MOGEN/${tmp}_matrix_log.txt"
	newlogname="/storage/htc/bdm/tosin/GSDB/Structures/ENCSR662QKG_VH9178CA/GSE92819_ENCFF876LAW/VC_5kb/MOGEN/${tmp}.log"
	mv $oldlogname $newlogname
done
