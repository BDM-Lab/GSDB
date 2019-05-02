#!/bin/sh
#!/bin/bash -l
#SBATCH -J miniMDS_run
#SBATCH -o miniMDS_run-%j.out
#SBATCH -p Lewis
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 80G
#SBATCH -t 2-00:00:00
#SBATCH --licenses=matlab:1 
n=6
rand=$(tr -cd '[:alnum:]' < /dev/urandom | head -c$n)
matlab_name="call_MDS_${rand}.m"
call_script="/storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/miniMDS/Access/${matlab_name}";
 #Call the miniMDS alg_parameters	
NAME="chr9_miniMDS"
outname="/storage/htc/bdm/tosin/GSDB/Structures/OO7429SF/GM12878/KR_500kb/miniMDS/${NAME}.tsv"
python  /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/miniMDS/minimds.py -m 0.01 -p 0.01 -o $outname  /storage/htc/bdm/tosin/GSDB/Data/OO7429SF/primary/GM12878_normalized/GM12878/KR_500kb/chr9_miniMDS.bed
#Call Assessment code here
Name="chr9"
listName="/storage/htc/bdm/tosin/GSDB/Data/OO7429SF/primary/GM12878_normalized/GM12878/KR_500kb/${Name}_list.txt"
matrixName="/storage/htc/bdm/tosin/GSDB/Data/OO7429SF/primary/GM12878_normalized/GM12878/KR_500kb/${Name}_matrix.txt"
Res=`sed -n '2{p;q}' ${outname}`
echo "Resolution = ${Res}"
outlog="/storage/htc/bdm/tosin/GSDB/Structures/OO7429SF/GM12878/KR_500kb/miniMDS/${Name}.log"
newfname="/storage/htc/bdm/tosin/GSDB/Structures/OO7429SF/GM12878/KR_500kb/miniMDS/${Name}.pdb"
#Create algorithm parameters
{
echo "filename='${outname}';"
echo "fname='${listName}';"
echo "filenm='${matrixName}';"
echo "outname='${newfname}';"
echo "logout='${outlog}';"
echo "Res=${Res};"
echo "extract_pdb_scc"
}> $call_script
## Module Commands
module load matlab/matlab-R2018a
module list
# Run matlab non-interactively
SCRIPT=$call_script
srun matlab -nodesktop -nosplash -nodisplay -r "run('${SCRIPT}');exit"
