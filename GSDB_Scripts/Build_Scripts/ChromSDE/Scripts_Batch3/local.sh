#!/bin/sh
#!/bin/bash -l
#SBATCH -J ChromSDE_run
#SBATCH -o ChromSDE_run-%j.out
#SBATCH -p Lewis
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 80G
#SBATCH -t 2-00:00:00
#SBATCH --licenses=matlab:1 
echo "structure construction started....."
echo "Processing /storage/htc/bdm/tosin/GSDB/Data/OO7429SF/primary/GM12878_250kb/GM12878/KR_500kb/chr9.n_contact file ........."
alg_dir="/storage/htc/bdm/tosin/GSDB/Structures/OO7429SF/GM12878/KR_500kb/ChromSDE"
NAME="chr9"
n=6
rand=$(tr -cd '[:alnum:]' < /dev/urandom | head -c$n)
matlab_name="call_ChromSDE_${rand}.m"
call_script="/storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/ChromSDE/${matlab_name}";
filename="/storage/htc/bdm/tosin/GSDB/Data/OO7429SF/primary/GM12878_250kb/GM12878/KR_500kb/${NAME}"
outname="${alg_dir}/${NAME}"
#Create algorithm parameters
{
echo "filename='${filename}';"
echo "out='${outname}';"
echo "[binAnno,normFreqMat]=readpipeline_output(filename);"
echo "ChromSDE_new(binAnno,normFreqMat,0,out);"
} > $call_script
## Module Commands
module load matlab/matlab-R2018a
module list
# Run matlab non-interactively
SCRIPT=$call_script
srun matlab -nodesktop -nosplash -nodisplay -r "run('${SCRIPT}');exit"
rm $call_script
