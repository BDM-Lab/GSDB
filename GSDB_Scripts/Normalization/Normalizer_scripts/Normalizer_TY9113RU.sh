#!/bin/sh

#!/bin/bash -l
#SBATCH -J Run_Normalizer_run
#SBATCH -o Run_Normalizer_run-%j.out
#SBATCH -p Lewis
#SBATCH -N 1
#SBATCH -n 3
#SBATCH --mem 40G
#SBATCH --mail-type=end
#SBATCH -t 2-00:00:00

# Determine the number folders/Directory
base_dir="/storage/htc/bdm/tosin/GSDB/Data/ENCSR244BBG_TY9113RU/Extracted_Data/*"
norm="VC"

echo "base directory = ${base_dir}"



for directory in $base_dir; do
        [ -d "${directory}" ] && ((directories++))
		if [ -d "$directory" ] 
		then
			echo "Directory Name: $directory "	
			chr_path=$directory			
			# make directory in outputs
			namedir="${chr_path##*/}"    #Get the name	
			path_to_output=$output_dir$namedir  #Create path to output
						
			
			# Read files with only chr -prefix
			normalized_path="${chr_path}/${norm}" 
			mkdir $normalized_path
			
			chr_path+="/chr*"
			genome_file="${directory}/${namedir}_All_Genome.txt"
			genome_length_file="${directory}/${namedir}_chrom_sequence_length.txt"
			genome_mapping_file="${directory}/${namedir}_mapping.txt"
					
			#copy the file content inside new directory
			
			chromseq_normalized_outpath="${normalized_path}/${namedir}_chrom_sequence_length.txt"
			mapping_normalized_outpath="${normalized_path}/${namedir}_mapping.txt"
			
			
			if [ -f $chromseq_normalized_outpath ] ; then
				echo "Deleting file = $chromseq_normalized_outpath"
				rm $chromseq_normalized_outpath
			fi
			
			if [ -f $mapping_normalized_outpath ] ; then
				echo "Deleting file = $mapping_normalized_outpath"
				rm $mapping_normalized_outpath
			fi
			
			cat $genome_length_file >> $chromseq_normalized_outpath
			cat $genome_mapping_file >> $mapping_normalized_outpath
			
			
			genome_file_normalized_outpath="${normalized_path}/${namedir}_All_Genome"
			chromosome_number="All"
			java -jar Normalizer_v1.jar  -i $genome_file -o $genome_file_normalized_outpath -c  $chromosome_number
			
			: 'Multiline comment:				
			files=0	
			
			for local_file in $chr_path
			do				
				[ -f "${local_file}" ] && ((files++))
				#Normalize local file before prediction
				NAME=`basename "$local_file" .txt`
				normalized_outpath="${normalized_path}/${NAME}"
				echo "Normalized Data filepath = $normalized_outpath"
				chromosome_number=`echo $NAME | cut -c 4-`
				java -jar Normalizer_v1.jar -i $local_file -o $normalized_outpath -c  $chromosome_number
							
			done
			
			#construct genome structure
			echo ""
			echo "Processing ${genome_file} file..."	
			echo ""
						
			echo "Number of files: ${files-0}"			
			'

		fi
		
done

echo "Number of directories: ${directories-0}"