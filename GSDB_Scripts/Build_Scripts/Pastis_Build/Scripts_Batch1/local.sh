#!/bin/sh
#SBATCH -J pastis_run
#SBATCH -o pastis_run-%j.out
#SBATCH -p Interactive
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 2G
#SBATCH -t 0-02:00:00
#SBATCH --licenses=matlab:1 
echo "structure construction started....."
echo "Processing /storage/htc/bdm/tosin/GSDB_Scripts/example/GSE105710_ENCFF278TZH/VC/chrY.n_contact file ........."
alg_dir="/storage/htc/bdm/tosin/GSDB_Scripts/Outputs/GSE105710_ENCFF278TZH/VC/Pastis"
NAME="chrY"
n=6
rand=$(tr -cd '[:alnum:]' < /dev/urandom | head -c$n)
matlab_name="call_Pastis_${rand}.m"
call_script="/storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/pastis/Evaluate/${matlab_name}";
dir_name="${NAME}_${rand}"
mkdir $dir_name
cbins="/storage/htc/bdm/tosin/GSDB_Scripts/example/GSE105710_ENCFF278TZH/VC/${NAME}.cbins "
counts="$dir_name/${NAME}.matrix"
lenght="$dir_name/${NAME}.bed"
cp /storage/htc/bdm/tosin/GSDB_Scripts/example/GSE105710_ENCFF278TZH/VC/chrY.n_contact $counts
cp $cbins $lenght
cd $dir_name
pwd
{
echo "[all]"
echo "output_name: structure"
echo "verbose: 1"
echo "max_iter: 100"
echo "counts:${NAME}.matrix"
echo "lengths: ${NAME}.bed"
echo "normalize: True"
} > config.ini
## Module Commands
 pastis-pm2 .
echo "Pastis structure generated Succesffully..............."
cord="/storage/htc/bdm/tosin/GSDB_Scripts/Build_Scripts/Pastis_Build/Scripts_Batch1/${dir_name}/PM2.structure"
#Evaluate the result
 MatrixNAME="/storage/htc/bdm/tosin/GSDB_Scripts/example/GSE105710_ENCFF278TZH/VC/${NAME}_matrix.txt"
 echo $MatrixNAME
#Create algorithm parameters
newfname="/storage/htc/bdm/tosin/GSDB_Scripts/Outputs/GSE105710_ENCFF278TZH/VC/Pastis/${NAME}.pdb"
outlog="/storage/htc/bdm/tosin/GSDB_Scripts/Outputs/GSE105710_ENCFF278TZH/VC/Pastis/${NAME}.log"
{
echo "pdbout='${newfname}';"
echo "pdblog='${outlog}';" 

echo "XYZ=dlmread('${cord}');"
echo "Data2=dlmread('$MatrixNAME');"
echo "fprintf('Previous Length = %d\n',length(XYZ)); "
echo "distance;"
echo "%output pdb and image  "
echo "output_structure;  "
echo " % correlate structure"
echo "[dist] = StructureDistance(XYZ,n); % distance from structure "
echo "n=length(XYZ); "
echo "fprintf('New Length = %d\n',n); "
echo "Data2(ntfd_index,:)=[];"
echo "Data2(:,ntfd_index)=[];"
echo "if (n~=length(Data2))"
echo "	disp('The Matrices are not of the same length');"
echo "end"
echo "%convert IF to distance "
echo " for j = 1:length(Data2)  "
echo "	  for k = j:n "
echo "	 	if(Data2(j,k)~=0  ) "
echo "	 		Data(j,k)=1./(Data2(j,k));"
echo "	 		Data(k,j)= Data(j,k);"
echo "	 	end"
echo "	  end"
echo "end"
echo " [RHO,RHO2] = Spearman_corr(Data,dist,Data2,n);"
echo "fileID=fopen(pdblog,'w');  "
echo "fprintf(fileID, 'Spearman correlation Expected Dist vs. Reconstructed Dist: %f\n',RHO);;"
echo "fclose(fileID);"
} > $call_script
## Module Commands
module load matlab/matlab-R2018a
module list
# Run matlab non-interactively
SCRIPT=$call_script
srun matlab -nodesktop -nosplash -nodisplay -r "run('${SCRIPT}');exit"
rm $call_script
rm -r $(pwd)
