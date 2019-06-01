#!/bin/bash

#==============================================================================================================
#Input:
# (1) Input Matrix textfile (textfile): Specify the path to  input Matrix
# (2) Pdb file (.pdb): Specify the path to  the pdb
# (3) Destination directory (dest_dir): Specify the path to the output directory to place the log file
#  ./IF_Structure_evaluate.sh [matrix] [pdb] [output]
#==============================================================================================================


echo "******************************************************************************"
echo "Assessing the similarity between Chromosome IF Matrix and Structure "
echo "******************************************************************************"
echo ""

if [ $# -lt 4 ]
  then
    echo "Not enough input!. function needs 2 input."
	echo "Input 1: IF Matrix "
	echo "Input 2: Structure in Pdb format "
	exit
	
else

	matrix_file="$1"
	pdb_file="$2"
	dest_dir="$3"
	jobid="$4";

	# Loop to copy

	/usr/bin/java -jar /var/www/html/3dgenome/GSDB/evaluate/IF_Structure_Evaluator.jar $matrix_file $pdb_file $dest_dir

	#Change output name to include the job id

	#change output name
	old_name="$dest_dir/Result.log"
	new_name="$dest_dir/${jobid}_Result.log"

	mv $old_name $new_name

fi



echo ""

echo "******************************************************************************"
echo "STATUS: Task  Completed"
echo "******************************************************************************"
