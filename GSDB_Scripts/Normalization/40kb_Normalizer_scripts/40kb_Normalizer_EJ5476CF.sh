#!/bin/sh

#!/bin/bash -l
#SBATCH -J Normalizer_run_40kb
#SBATCH -o Normalizer_run_40kb-%j.out
#SBATCH -p Interactive
#SBATCH -N 1
#SBATCH -n 3
#SBATCH --mem 250G
#SBATCH --mail-type=end
#SBATCH -t 0-04:00:00

# Determine the number folders/Directory
norm="VC"
directory="/storage/htc/bdm/tosin/GSDB/Data/ENCSR982KWR_EJ5476CF/Extracted_Data/"
## declare an array variable
declare -a arr=("GSE105988_ENCFF925ITS")
# Read files with only chr -prefix

	for namedir in "${arr[@]}"
	do
		echo "Processing directory $namedir..........."
		normalized_path="${directory}/${namedir}/${norm}" 
		mkdir $normalized_path
		genome_file="${directory}/${namedir}/${namedir}_All_Genome.txt"
		genome_length_file="${directory}/${namedir}/${namedir}_chrom_sequence_length.txt"
		genome_mapping_file="${directory}/${namedir}/${namedir}_mapping.txt"
				
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
		java -Xms10g -Xmx200g -jar Normalizer_v1.jar  -i $genome_file -o $genome_file_normalized_outpath -c  $chromosome_number

	done