#!/bin/sh

directory="/storage/htc/bdm/tosin/GSDB_Scripts/example/GSE105710_ENCFF063WNA/VC"

alg_dir="/storage/htc/bdm/tosin/GSDB_Scripts/Outputs/GSE105710_ENCFF063WNA/VC/MOGEN"

local_file="/storage/htc/bdm/tosin/GSDB_Scripts/example/GSE105710_ENCFF063WNA/VC/chr5_list.txt"
#copy the model out to the outside directory
NAME=`basename "$local_file" _list.txt`
echo $NAME
pdbfile="${alg_dir}/${NAME}.pdb"												
MatrixNAME="${directory}/${NAME}_matrix.txt"

echo $MatrixNAME

# Analysis									
java -jar /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/pdbAnalyze.jar $MatrixNAME $pdbfile $alg_dir

oldlogname="${alg_dir}/${NAME}_matrix_log.txt"	
newlogname="${alg_dir}/${NAME}.log"	

mv $oldlogname $newlogname