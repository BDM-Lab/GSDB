#!/bin/bash

#==============================================================================================================
#Input:
# (1) Pdb file1 (.pdb): Specify the path to  the pdb
# (2) Pdb file2 (.pdb): Specify the path to  the pdb
# (3) Destination directory (dest_dir): Specify the path to the output directory to place the log file
#  ./Structure_Structure_evaluate.sh [pdb1] [pdb2]  [output]
#==============================================================================================================


# Specify the pawprint textfile, base_dir and dest_dir

echo "==================================="
echo "STATUS: Evaluation started"
echo "==================================="
echo ""

if [ $# -lt 4 ]
  then
    echo "Not enough input!. Script needs 3 input. Use the format below to execute the evaluation script."
	echo "To execute: ./Structure_Structure_evaluate.sh [pdb1] [pdb2] [Dest_dir] [jobid]"
	exit
	
else

	pdb_file1="$1"
	pdb_file2="$2"
	dest_dir="$3"
	jobid="$4";

	# Loop to copy

	/usr/bin/java -jar /var/www/html/3dgenome/GSDB/evaluate/Structure_Structure_Evaluator.jar $pdb_file1 $pdb_file2 $dest_dir

	#Change output name to include the job id

	#change output name
	old_name="$dest_dir/Result.log"
	new_name="$dest_dir/${jobid}_Result.log"

	mv $old_name $new_name

fi



echo ""

echo "==================================="
echo "STATUS: Evaluation Completed"
echo "==================================="


