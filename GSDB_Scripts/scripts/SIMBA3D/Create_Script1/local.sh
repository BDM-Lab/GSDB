#!/bin/sh
#!/bin/bash -l
#SBATCH -J SIMBA3D_run
#SBATCH -o SIMBA3D_run-%j.out
#SBATCH -p Interactive
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 4G
#SBATCH -t 0-04:00:00
n=6
rand=$(tr -cd '[:alnum:]' < /dev/urandom | head -c$n)
python_name="call_SIMBA3D_${rand}.py"
call_script="${python_name}";
#Create algorithm parameters
filename="/storage/htc/bdm/tosin/GSDB/Data/OO7429SF/primary/GM12878_normalized/GM12878/KR_50kb/chr23_matrix.txt"
echo "Processing file ${filename}.............."
new_filename="/storage/htc/bdm/tosin/GSDB/Data/OO7429SF/primary/GM12878_normalized/GM12878/KR_50kb/chr23"
{
echo "import numpy as np"
echo "data=np.loadtxt('$filename')"
echo "np.save('$new_filename',data)"
echo""
} > $call_script
python $call_script
echo "Python Script Execution completed............."
rm $call_script
