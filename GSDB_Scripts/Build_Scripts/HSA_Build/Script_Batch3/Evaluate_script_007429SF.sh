#!/bin/sh


# Determine the number folders/Directory
base_dir="/storage/htc/bdm/tosin/GSDB/Data/OO7429SF/primary/GM12878_250kb/GM12878/KR_250kb/"
output_dir="/storage/htc/bdm/tosin/GSDB/Structures/OO7429SF/GM12878/KR_250kb/HSA/"



for local_file in $output_dir/chr*.pdb
do				
	[ -f "${local_file}" ] && ((files++))
	fnme=`basename "$local_file" .pdb`
	#echo "filename = ${fnme}"
	matrixfile="$base_dir${fnme}_matrix.txt"
	#inputfile="$output_dir${fnme}_HSA_out.txt"
	echo "matrix file = $matrixfile"
	echo "inputfile = $local_file"
	#rm "$output_dir${fnme}_HSA.pdb"
	java -jar /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/pdbAnalyze.jar $matrixfile $local_file $output_dir
	#java -jar /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/HSA/HSA2PDB.jar $matrixfile $inputfile $output_dir
	oldname="$output_dir/${fnme}_matrix_log.txt"
	newname="$output_dir/${fnme}.log"
	mv $oldname $newname
done


