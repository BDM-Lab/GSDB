#!/bin/bash -l
#SBATCH -J SIMBA3D_run
#SBATCH -o SIMBA3D_run-%j.out
#SBATCH -p Lewis
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem 90G
#SBATCH -t 2-00:00:00

#Generate a random parameter
n=6
rand=$(tr -cd '[:alnum:]' < /dev/urandom | head -c$n)
param="parameters_SIMBA3D_${rand}.json"
call_script="get_SIMBA3D_XYZ_${rand}.py"
#Create algorithm parameters
{
echo "["
echo "	{"
echo "	\"taskname\":\"SIMBA3D_chrY\","
echo "	\"description\":\"json script for each chromosome construction\","
echo "	\"uuid\":\"oVLxGn\","
echo "	\"file_names\": {"
echo "		\"inputdir\": \"/storage/htc/bdm/tosin/GSDB/Data/ENCSR998ZSP_XF4844ZB/Extracted_Data/GSE106015_ENCFF777JOF/VC/\","
echo "		\"outputdir\": \"/storage/htc/bdm/tosin/GSDB/Structures/XF4844ZB/GSE106015_ENCFF777JOF/VC/SIMBA3D/\","
echo "		\"pairwise_contact_matrix\":\"chrY.npy\","
echo "		\"output_filename\" : \"chrY.npz\""
echo "		}"
echo "	}"
echo " ]"
} > $param
outname="/storage/htc/bdm/tosin/GSDB/Structures/XF4844ZB/GSE106015_ENCFF777JOF/VC/SIMBA3D/chrY.npz"
# Call the SIMBA3D alg_parameters
matname="/storage/htc/bdm/tosin/GSDB/Structures/XF4844ZB/GSE106015_ENCFF777JOF/VC/SIMBA3D/chrY.mat"
xyzname="/storage/htc/bdm/tosin/GSDB/Structures/XF4844ZB/GSE106015_ENCFF777JOF/VC/SIMBA3D/chrY.xyz"
simba3d -r $param
rm $param
# convert npz to mat
simba3d-convertion ${outname} --ext_out .mat
#Evaluate Result
 echo " Java Evaluation Started............"
 MatrixNAME="/storage/htc/bdm/tosin/GSDB/Data/ENCSR998ZSP_XF4844ZB/Extracted_Data/GSE106015_ENCFF777JOF/VC/chrY_matrix.txt"
 echo $MatrixNAME
newfname="/storage/htc/bdm/tosin/GSDB/Structures/XF4844ZB/GSE106015_ENCFF777JOF/VC/SIMBA3D/chrY.pdb"
outlog="/storage/htc/bdm/tosin/GSDB/Structures/XF4844ZB/GSE106015_ENCFF777JOF/VC/SIMBA3D/chrY.log"
matfile="/storage/htc/bdm/tosin/GSDB/Structures/XF4844ZB/GSE106015_ENCFF777JOF/VC/SIMBA3D/chrY.mat"
{
echo "import scipy.io as spio"
echo "import numpy as np"
echo "mat = spio.loadmat('${matname}', squeeze_me=True)"
echo "M = mat['X_evol'] "
echo "G = np.transpose(M)"
echo "np.savetxt('${xyzname}', G,fmt='%5.5f' ,delimiter='\t') "
} > $call_script
## Module Commands
python $call_script
#Evaluate algorithm parameters
newpdbname="/storage/htc/bdm/tosin/GSDB/Structures/XF4844ZB/GSE106015_ENCFF777JOF/VC/SIMBA3D/chrY.pdb"
outlog="/storage/htc/bdm/tosin/GSDB/Structures/XF4844ZB/GSE106015_ENCFF777JOF/VC/SIMBA3D/chrY.log"
# Analysis
java -jar /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/XYZ2PDB.jar $MatrixNAME ${xyzname} /storage/htc/bdm/tosin/GSDB/Structures/XF4844ZB/GSE106015_ENCFF777JOF/VC/SIMBA3D
oldlogname="/storage/htc/bdm/tosin/GSDB/Structures/XF4844ZB/GSE106015_ENCFF777JOF/VC/SIMBA3D/chrY_matrix_log.txt"
oldpdbname="/storage/htc/bdm/tosin/GSDB/Structures/XF4844ZB/GSE106015_ENCFF777JOF/VC/SIMBA3D/chrY_matrix.pdb"	
mv $oldlogname $outlog
mv $oldpdbname $newpdbname
rm $call_script
