#!/bin/bash -l
#SBATCH -J ShRec3D_run
#SBATCH -o ShRec3D_run-%j.out
#SBATCH -p Lewis
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 90G
#SBATCH -t 2-00:00:00
#SBATCH --licenses=matlab:1 
n=6
rand=$(tr -cd '[:alnum:]' < /dev/urandom | head -c$n)
matlab_name="call_ShRec3D_${rand}.m"
call_script="/storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/ShRec3D/${matlab_name}";
#Use Normalized local file before prediction
echo 
echo "Processing /storage/htc/bdm/tosin/GSDB/Data/ENCSR662QKG_VH9178CA/Extracted_Data/GSE92819_ENCFF876LAW/VC_5kb/chrY_matrix.txt file ........."
echo 
NAME="chrY"
newfname="/storage/htc/bdm/tosin/GSDB/Structures/ENCSR662QKG_VH9178CA/GSE92819_ENCFF876LAW/VC_5kb/ShRec3D/${NAME}.pdb"
outlog="/storage/htc/bdm/tosin/GSDB/Structures/ENCSR662QKG_VH9178CA/GSE92819_ENCFF876LAW/VC_5kb/ShRec3D/${NAME}.log"
#Create algorithm parameters
{
echo "filename='/storage/htc/bdm/tosin/GSDB/Data/ENCSR662QKG_VH9178CA/Extracted_Data/GSE92819_ENCFF876LAW/VC_5kb/chrY_matrix.txt';"
echo "outname='${newfname}';"
echo "logout='${outlog}';" 
echo "ShRec3D"
} > $call_script
## Module Commands
module load matlab/matlab-R2018a
module list
# Run matlab non-interactively
SCRIPT=$call_script
srun matlab -nodesktop -nosplash -nodisplay -r "run('${SCRIPT}');exit"
rm $call_script
