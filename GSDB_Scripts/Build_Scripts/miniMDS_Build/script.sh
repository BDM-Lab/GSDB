#!/bin/sh


n=6
rand=$(tr -cd '[:alnum:]' < /dev/urandom | head -c$n)
matlab_name="extract_pdb_scc_${rand}.m"
call_script="/storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/miniMDS/Access/${matlab_name}";

directory="/storage/htc/bdm/tosin/GSDB_Scripts/example/GSE105710_ENCFF063WNA/VC"

alg_dir="/storage/htc/bdm/tosin/GSDB_Scripts/Outputs/GSE105710_ENCFF063WNA/VC/miniMDS"

local_file="/storage/htc/bdm/tosin/GSDB_Scripts/example/GSE105710_ENCFF063WNA/VC/chr5_miniMDS.bed"

 #Call the miniMDS alg_parameters	
NAME=`basename "$local_file" .bed`

outname="${alg_dir}/${NAME}.tsv"
#Call Assessment code here
Name=`basename "$local_file" _miniMDS.bed`	
listName="${directory}/${Name}_list.txt"
matrixName="${directory}/${Name}_matrix.txt"						
Res=`sed -n '2{p;q}' ${outname}`
echo "Resolution = ${Res}"

outlog="${alg_dir}/${Name}.log"
newfname="${alg_dir}/${Name}.pdb"

echo "matlab script = $call_script"