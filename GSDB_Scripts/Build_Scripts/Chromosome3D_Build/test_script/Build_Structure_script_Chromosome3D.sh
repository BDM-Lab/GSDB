#!/bin/sh

#!/bin/bash -l
#SBATCH -J Chromosome3D_run
#SBATCH -o Chromosome3D_run-%j.out
#SBATCH -p Interactive
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu 5G
#SBATCH --mail-type=end
#SBATCH -t 0-04:00:00


# Determine the number folders/Directory
base_dir="/storage/htc/bdm/tosin/GSDB_Scripts/example/*"
output_dir="/storage/htc/bdm/tosin/GSDB_Scripts/Outputs/"
algorithm="Chromosome3D"


echo "base directory = ${base_dir}"


for dir in $base_dir; do 
    [ -d "${dir}" ] && ((dir++))
	if [ -d "$dir" ] 
	then
		rootdir=$dir	
	    rootname="${rootdir##*/}"    #Get the name
		echo "Root Directory Name: $rootname "	
		rootdir+="/*"
		path_to_out=$output_dir$rootname  #Create path to output
		mkdir $path_to_out 
		
		echo "Directory created successfully: $path_to_output"
		#nested for loop for new directory
		for directory in $rootdir; do
				[ -d "${directory}" ] && ((directories++))
				if [ -d "$directory" ] 
				then			
					echo "Nest Directory Path: $directory "	
					chr_path=$directory			
					# make directory in outputs
					namedir="${chr_path##*/}"    #Get the name	
					echo "Nest directory Name:  ${namedir}"					
					path_to_output="${path_to_out}/${namedir}"  #Create path to output
					mkdir $path_to_output	
					
					echo "Directory created successfully........"
					
					#Crate directory for algorithm in each folder
					alg_dir="${path_to_output}/${algorithm}"
					mkdir $alg_dir

					# Read files with only chr -prefix
							
					chr_path+="/chr*_matrix.txt"
					genome_file="${directory}/${rootname}_All_Genome_matrix.txt"
					genome_length_file="${directory}/${rootname}_chrom_sequence_length.txt"
					value=$(<$genome_length_file) #Load length file
					echo "Length: ${value}"
					echo "Chromosome Directory Name: $chr_path"	
					
					# Read the Genome sequesce file
					echo ""
					echo "Construct Chromosome Structure"				
					files=0	
					
					for local_file in $chr_path
					do				
						[ -f "${local_file}" ] && ((files++))
						
						#Use Normalized local file before prediction
						echo ""
						echo "Running job for $local_file... "	
						echo ""
						
						 #Call the miniMDS alg_parameters	
						NAME=`basename "$local_file" .txt`
						outname="${alg_dir}/${NAME}"
						echo ""	
						/storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/Chromosome3D/chromosome3D.pl -i $local_file -o $outname 
						
						#copy the model out to the outside directory
						pdbfile="${alg_dir}/${NAME}/${NAME}_model1.pdb"
						NAME=`basename "$pdbfile" _matrix_model1.pdb`						
						pdbNAME="${NAME}.pdb"
						newpdbfile="${alg_dir}/${pdbNAME}"
						echo $newpdbfile
						cp $pdbfile $newpdbfile
						
						# Analysis									
						java -jar /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/pdbAnalyze.jar $local_file $newpdbfile $alg_dir		
						#Change log_Name
						oldlogname="${alg_dir}/${NAME}_matrix_log.txt"	
						newlogname="${alg_dir}/${NAME}.log"	
						mv $oldlogname $newlogname
					done
		
					: 'Multiline comment:
					#construct genome structure
					echo "Construct Genome Structure"
					echo "Processing ${genome_file} file..."	
					echo ""
					
					#Create algorithm genome_parameters	
					genomepath="${alg_dir}/genome" 
					mkdir $genomepath	
					
					NAME=`basename "${genome_file}" .txt`
					outname="${genomepath}/${NAME}"
					/storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/Chromosome3D/chromosome3D.pl -i $genome_file -o $outname 
					
					#copy the model out to the outside directory
					pdbfile="${genomepath}/${NAME}/${NAME}_model1.pdb"
					NAME=`basename "$pdbfile" _matrix_model1.pdb`						
					pdbNAME="${NAME}.pdb"
					newpdbfile="${genomepath}/${pdbNAME}"
					echo $newpdbfile
					cp $pdbfile $newpdbfile
					# Analysis									
					java -jar /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/pdbAnalyze.jar $local_file $newpdbfile $genomepath
					
					#Change log_Name
					oldlogname="${genomepath}/${NAME}_matrix_log.txt"	
					newlogname="${genomepath}/${NAME}.log"	
					mv $oldlogname $newlogname
					
					echo "Number of files: ${files-0}"	
					'

				fi
				
		done
	
	fi
done
	
echo ""

echo "miniMDS Run Successfully!!!"

# For each directory 
	# Get the name of the folders/Directory
	# Create a new directory/folder with name under the Parent Directory Structure
	# Get the chromosome data in the directory
	# Construct Structure
	# Structure name equals parent directory_chromosome name
	# For genome modelling - A solution will be to 
			#(1) Load the length file
			#(2) Create a new parameter file 
			#(3) Naming stucture folows the above