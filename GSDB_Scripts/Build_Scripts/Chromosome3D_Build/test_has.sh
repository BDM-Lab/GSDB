#!/bin/sh
local_file="/storage/htc/bdm/tosin/GSDB_Scripts/example/GSE105710_ENCFF063WNA/VC/chr5_matrix.txt"
NAME="chr5_matrix"
alg_dir="/storage/htc/bdm/tosin/GSDB_Scripts/Outputs/GSE105710_ENCFF063WNA/VC/Chromosome3D"
#copy the model out
pdbfile="${alg_dir}/${NAME}/${NAME}_model1.pdb"
NAME=`basename "$pdbfile" _matrix_model1.pdb`						
pdbNAME="${NAME}.pdb"
newpdbfile="${alg_dir}/${pdbNAME}"
echo $newpdbfile
cp $pdbfile $newpdbfile
# Analysis									
java -jar /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/pdbAnalyze.jar $local_file $newpdbfile $alg_dir	