#!/bin/sh
#!/bin/bash -l
#SBATCH -J InfMod3DGen_run
#SBATCH -o InfMod3DGen_run-%j.out
#SBATCH -p Lewis
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 90G
#SBATCH -t 2-00:00:00
#SBATCH --licenses=matlab:1

#Generate a random parameter
n=6
rand=$(tr -cd '[:alnum:]' < /dev/urandom | head -c$n)
matlab_name="call_InfMod3DGen_${rand}.m"
call_script="/storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/InfMod3DGen/${matlab_name}";
#Create algorithm parameters
MatrixNAME="/storage/htc/bdm/tosin/GSDB/Data/WT9059TG/Cheng_Prop_WT/GSM455133_30E0LAAXX.1/VC_1mb/chr21_matrix.txt"
fname="/storage/htc/bdm/tosin/GSDB/Data/WT9059TG/Cheng_Prop_WT/GSM455133_30E0LAAXX.1/VC_1mb/chr21_list.txt"
new_fname="/storage/htc/bdm/tosin/GSDB/Data/WT9059TG/Cheng_Prop_WT/GSM455133_30E0LAAXX.1/VC_1mb/chr21_InfMod3DGen.txt"
outname="/storage/htc/bdm/tosin/GSDB/Structures/WT9059TG/GSM455133_30E0LAAXX.1/VC_1mb/InfMod3DGen/chr21"
#Matlab for Execution and Evaluation
 echo " MATLAB Evaluation Started............"
newfname="/storage/htc/bdm/tosin/GSDB/Structures/WT9059TG/GSM455133_30E0LAAXX.1/VC_1mb/InfMod3DGen/chr21.pdb"
outlog="/storage/htc/bdm/tosin/GSDB/Structures/WT9059TG/GSM455133_30E0LAAXX.1/VC_1mb/InfMod3DGen/chr21.log"
matfile="/storage/htc/bdm/tosin/GSDB/Structures/WT9059TG/GSM455133_30E0LAAXX.1/VC_1mb/InfMod3DGen/chr21.mat"
echo "Processing file ${fname}.............."
{
echo "ChrMod_main('${fname}','${new_fname}','${matfile}',21,1)"
echo "disp('InfMod3DGen Execution Completed........');"
echo "data=load('${matfile}');"
echo "B=data.Ensemble;"
echo "[n,m,p]=size(B)"
echo "B = reshape(B,[n*m p]);"
echo "pdbout='${newfname}';"
echo "pdblog='${outlog}';" 

echo "XYZ=B"
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
