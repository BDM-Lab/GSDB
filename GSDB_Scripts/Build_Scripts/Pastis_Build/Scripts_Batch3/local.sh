#!/bin/sh
#SBATCH -J pastis_run
#SBATCH -o pastis_run-%j.out
#SBATCH -p Lewis
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem 80G
#SBATCH -t 2-00:00:00
echo "structure construction started....."
echo "Processing /storage/htc/bdm/tosin/GSDB/Data/OO7429SF/primary/GM12878_250kb/GM12878/KR_500kb/chr9.n_contact file ........."
alg_dir="/storage/htc/bdm/tosin/GSDB/Structures/OO7429SF/GM12878/KR_500kb/Pastis"
NAME="chr9"
n=6
rand=$(tr -cd '[:alnum:]' < /dev/urandom | head -c$n)
dir_name="${NAME}_${rand}"
mkdir $dir_name
cbins="/storage/htc/bdm/tosin/GSDB/Data/OO7429SF/primary/GM12878_250kb/GM12878/KR_500kb/${NAME}.cbins "
counts="$dir_name/${NAME}.matrix"
lenght="$dir_name/${NAME}.bed"
cp /storage/htc/bdm/tosin/GSDB/Data/OO7429SF/primary/GM12878_250kb/GM12878/KR_500kb/chr9.n_contact $counts
cp $cbins $lenght
cd $dir_name
pwd
{
echo "[all]"
echo "output_name: structure"
echo "verbose: 1"
echo "counts:${NAME}.matrix"
echo "lengths: ${NAME}.bed"
echo "normalize: False"
} > config.ini
## Module Commands
 pastis-pm2 .
echo "Pastis structure generated Succesffully..............."
pdbfile="/storage/htc/bdm/tosin/GSDB_Scripts/Build_Scripts/Pastis_Build/Scripts_Batch3/${dir_name}/PM2.structure.pdb"
#Evaluate the result
 MatrixNAME="/storage/htc/bdm/tosin/GSDB/Data/OO7429SF/primary/GM12878_250kb/GM12878/KR_500kb/${NAME}_matrix.txt"
 echo $MatrixNAME
#Create algorithm parameters
newfname="/storage/htc/bdm/tosin/GSDB/Structures/OO7429SF/GM12878/KR_500kb/Pastis/${NAME}.pdb"
outlog="/storage/htc/bdm/tosin/GSDB/Structures/OO7429SF/GM12878/KR_500kb/Pastis/${NAME}.log"
mv $pdbfile $newfname
# Analysis
java -jar /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/pdbAnalyze.jar $MatrixNAME $newfname $alg_dir
oldlogname="${alg_dir}/${NAME}_matrix_log.txt"
newlogname="${alg_dir}/${NAME}.log"	
mv $oldlogname $newlogname
rm -r $(pwd)
