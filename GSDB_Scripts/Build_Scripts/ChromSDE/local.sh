#!/bin/sh
echo "structure construction started....."
echo "Processing /storage/htc/bdm/tosin/GSDB_Scripts/example/GSE105710_ENCFF278TZH/VC/chrY.n_contact file ........."
alg_dir="/storage/htc/bdm/tosin/GSDB_Scripts/Outputs/GSE105710_ENCFF278TZH/VC/ChromSDE"
NAME="chrY"
n=6
rand=$(tr -cd '[:alnum:]' < /dev/urandom | head -c$n)
matlab_name="call_ChromSDE_${rand}.m"
call_script="/storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/ChromSDE/${matlab_name}";
filename="/storage/htc/bdm/tosin/GSDB_Scripts/example/GSE105710_ENCFF278TZH/VC/${NAME}"
outname="${alg_dir}/${NAME}"
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
