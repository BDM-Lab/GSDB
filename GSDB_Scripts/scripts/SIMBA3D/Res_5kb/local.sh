#!/bin/sh
#!/bin/bash -l
#SBATCH -J SIMBA3D_run
#SBATCH -o SIMBA3D_run-%j.out
#SBATCH -p Lewis
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem 80G
#SBATCH -t 2-00:00:00
n=6
rand=$(tr -cd '[:alnum:]' < /dev/urandom | head -c$n)
python_name="call_SIMBA3D_${rand}.py"
call_script="${python_name}";
#Create algorithm parameters
filename="/storage/htc/bdm/tosin/GSDB/Data/ENCSR662QKG_VH9178CA/Extracted_Data/GSE92819_ENCFF876LAW/VC_5kb/chrM_matrix.txt"
echo "Processing file ${filename}.............."
new_filename="/storage/htc/bdm/tosin/GSDB/Data/ENCSR662QKG_VH9178CA/Extracted_Data/GSE92819_ENCFF876LAW/VC_5kb/chrM"
{
echo "import numpy as np"
echo "data=np.loadtxt('$filename')"
echo "np.save('$new_filename',data)"
echo""
} > $call_script
python $call_script
echo "Python Script Execution completed............."
rm $call_script
