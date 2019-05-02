#!/bin/bash -l
#SBATCH -J HSA_run
#SBATCH -o HSA_run-%j.out
#SBATCH -p hpc4
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem 50G
#SBATCH -t 2-00:00:00
alg_dir="/storage/htc/bdm/tosin/GSDB/Structures/OO7429SF/GM12878/KR_500kb/HSA"
NAME="chr9_HSA"
outname="/storage/htc/bdm/tosin/GSDB/Structures/OO7429SF/GM12878/KR_500kb/HSA/${NAME}_out"
#Generate a random parameter
n=6
rand=$(tr -cd '[:alnum:]' < /dev/urandom | head -c$n)
param="parameter_HSA_${rand}.txt"
echo "output file name PREFIX = ${outname}"
{
echo "options('expressions'=50000)"
echo "source('/storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/HSA/cstruct1.R')"
echo "mat=read.table(\"/storage/htc/bdm/tosin/GSDB/Data/OO7429SF/primary/GM12878_250kb/GM12878/KR_500kb/chr9_HSA.txt\",header=F,sep="")"
echo "lsmap0=as.matrix(mat)"
echo "outfile=\"${outname}\""
echo "lscov0=0"
echo "mak=0"
echo "lsmap0=vector(\"list\",1)"
echo "lsmap0[[1]]=as.matrix(mat)"
echo "out=fmain(lsmap0,lscov0,outfile,300,150,50,50,0.001,0,0,mak)"
} > $param
echo ""
#Call the HSA alg	
module load  R/R-3.3.1
/usr/bin/time --verbose R --no-save -f $param 
#convert to pdb and assessment
echo "pdb file conversion....."
java -jar /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/HSA/HSA2PDB.jar /storage/htc/bdm/tosin/GSDB/Data/OO7429SF/primary/GM12878_250kb/GM12878/KR_500kb/chr9_HSA.txt "${outname}.txt" /storage/htc/bdm/tosin/GSDB/Structures/OO7429SF/GM12878/KR_500kb/HSA
rm $param 
#Change log_Name
NAME="chr9"
oldlogname="/storage/htc/bdm/tosin/GSDB/Structures/OO7429SF/GM12878/KR_500kb/HSA/${NAME}_HSA_log.txt"
newlogname="/storage/htc/bdm/tosin/GSDB/Structures/OO7429SF/GM12878/KR_500kb/HSA/${NAME}.log"
mv $oldlogname $newlogname	
#Change pdb_Name
oldlogname="/storage/htc/bdm/tosin/GSDB/Structures/OO7429SF/GM12878/KR_500kb/HSA/${NAME}_HSA.pdb"
newlogname="/storage/htc/bdm/tosin/GSDB/Structures/OO7429SF/GM12878/KR_500kb/HSA/${NAME}.pdb"
mv $oldlogname $newlogname
