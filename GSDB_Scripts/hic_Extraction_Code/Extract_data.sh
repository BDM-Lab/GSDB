#!/bin/sh

#!/bin/bash -l
#SBATCH -J Extract_run
#SBATCH -o Extract_run-%j.out
#SBATCH -p Lewis
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 5G
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=oeow39@mail.missouri.edu
#SBATCH -t 2-00:00:00

# Determine the number folders/Directory
input_dir="/storage/htc/bdm/tosin/GSDB/Data/GSE92825/ProcessedData"
tool="/storage/htc/bdm/tosin/GSDB/hic_Extraction/juicer_tools.jar"
cur_dir="/storage/htc/bdm/tosin/GSDB/Data/GSE92825/Extracted_Data"
echo "input directory = ${input_dir}"

## declare an array variable
declare -a arr=("2500000" "1000000" "500000" "250000" "100000" "50000" "25000" "10000" "5000")

#Unzip the input files

for file in "${input_dir}"/*.hic.gz
do
	echo "Processing $file"	
	gunzip $file	
done

#finds files with same extention and loop through

for file in "${input_dir}"/*.hic
do
	echo "Processing $file"		
	tmp=$(echo "$file" | awk -F '_' '{print $2}' )	
	work_dir="${cur_dir}/${tmp}"
	echo "HiC Data = ${tmp} ......"
	mkdir $work_dir
	# do something with file
	for i in `seq 1 22`; 
	do
		echo "Processing chromosome $i..........."
		## now loop through the above array
		for a in "${arr[@]}"
		do
			#echo $a
			# or do whatever with individual element of the array
			fname="${work_dir}/chr${i}_${a}.txt"				   
			java -jar $tool dump observed VC $file $i $i BP $a $fname
			 
		done		

	done
		
	echo "Processing chromosome M............."
	## Chromosome M
	for a in "${arr[@]}"
	do
		#echo "resolution = $a"
		# or do whatever with individual element of the array
		fM="${work_dir}/chrM_${a}.txt"			   
		java -jar $tool dump observed VC $file M M BP $a $fM			   
	  
	done	
	
	echo "Processing chromosome X............."
	# Chromosome X
	for a in "${arr[@]}"
		do
		#echo "resolution = $a"
		# or do whatever with individual element of the array
		fX="${work_dir}/chrX_${a}.txt"			   
		java -jar $tool dump observed VC $file X X BP $a $fX			   
	  
	done
	
	echo "Processing chromosome Y............."
	# Chromosome Y
	for a in "${arr[@]}"
		do
		#echo "resolution = $a"
		# or do whatever with individual element of the array
		fY="${work_dir}/chrY_${a}.txt"			   
		java -jar $tool dump observed VC $file Y Y BP $a $fY			   
	  
	done
	  
done 



  
  		