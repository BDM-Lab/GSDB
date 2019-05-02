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
n=6
rand=$(tr -cd '[:alnum:]' < /dev/urandom | head -c$n)
matlab_name="call_ChromSDE_${rand}.m"
call_script="${matlab_name}";
#Create algorithm parameters
{
echo "filepath='/storage/htc/bdm/tosin/GSDB/Data/OO7429SF/primary/GM12878_250kb/GM12878/KR_500kb';"
echo "output_path='/storage/htc/bdm/tosin/GSDB/Data/OO7429SF/primary/GM12878_250kb/GM12878/KR_500kb';"
echo "filename = [filepath,'/chr',num2str(23),'_HSA.txt'];"
echo "disp (['Running Job for filename =', filename]);"
echo "contact = dlmread(filename);"
echo "nn=23;"
echo "ChromSDE_input;"
} > $call_script
## Module Commands
module load matlab/matlab-R2018a
module list
# Run matlab non-interactively
SCRIPT=$call_script
srun matlab -nodesktop -nosplash -nodisplay -r "run('${SCRIPT}');exit"
rm $call_script
