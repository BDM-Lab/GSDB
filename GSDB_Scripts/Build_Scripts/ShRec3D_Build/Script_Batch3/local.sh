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
echo "Processing /storage/htc/bdm/tosin/GSDB/Data/OO7429SF/primary/GM_100kb/GM12878/KR_100kb/chr9_matrix.txt file ........."
echo 
NAME="chr9"
newfname="/storage/htc/bdm/tosin/GSDB/Structures/OO7429SF/GM12878/KR_100kb/ShRec3D/${NAME}.pdb"
outlog="/storage/htc/bdm/tosin/GSDB/Structures/OO7429SF/GM12878/KR_100kb/ShRec3D/${NAME}.log"
#Create algorithm parameters
{
echo "filename='/storage/htc/bdm/tosin/GSDB/Data/OO7429SF/primary/GM_100kb/GM12878/KR_100kb/chr9_matrix.txt';"
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
