#!/bin/sh

#!/bin/bash -l
#SBATCH -J MOGEN_run
#SBATCH -o MOGEN_run-%j.out
#SBATCH -p Lewis
#SBATCH -N 1
#SBATCH -n 5
#SBATCH --mem 80G
#SBATCH --mail-type=end
#SBATCH -t 2-00:00:00

# Determine the number folders/Directory
base_dir="/storage/htc/bdm/tosin/GSDB/Data/VH2561BL/*"
output_dir="/storage/htc/bdm/tosin/GSDB/Structures/VH2561BL/"

algorithm="MOGEN"
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
					rm -r $alg_dir		
					mkdir $alg_dir

					# Read files with only chr -prefix
							
					chr_path+="/chr*_list.txt"
					genome_file="${directory}/${rootname}_All_Genome_list.txt"
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
						echo "Processing $local_file file... "	
						echo ""
						#Generate a random parameter
						n=6
						rand=$(tr -cd '[:alnum:]' < /dev/urandom | head -c$n)
						param="parameters_MOGEN_${rand}.txt"	
						#Create algorithm parameters
						{
						 echo "NUM = 1"
						 echo "NBR_OF_CHR = 1"
						 echo "INTRA_IF_THRESHOLD = 80%"
						 echo "ADJACENT_DIST = 1.5"
						 echo "CONTACT_DIST = 6.0"
						 echo "POS_MIN_DIST = 0.2"
						 echo "NEG_MAX_DIST_INTRA = 30"
						 echo "POS_MAX_DIST_WEIGHT_FILE = /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/parameters/pos_max_dist_weight_1_normal.txt"
						 echo "POS_MIN_DIST_WEIGHT_FILE = /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/parameters/pos_min_dist_weight_1_normal.txt"
						 echo "NEG_MIN_DIST_WEIGHT_FILE = /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/parameters/neg_min_dist_weight_1_normal.txt"
						 echo "NEG_MAX_DIST_WEIGHT_FILE = /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/parameters/neg_max_dist_weight_1_normal.txt"
						 echo ""
						 echo "OUTPUT_FOLDER = ${alg_dir}/"
						 echo "INPUT_FILE = ${local_file}"
						 echo "VERBOSE = false"				 
						 echo "LEARNING_RATE = 0.001"
						 echo "MAX_ITERATION = 200000"
						 
						} > $param
						
						# Call the LorDG alg_parameters	
						
						java -jar /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/MOGEN.jar $param
						echo "construction completed... "				
						rm $param
		
						
					done

					#Read all the pdb files in directory and print chromosome name only
					algorithm_dir="${alg_dir}/chr*_*.pdb"
					echo ".........................................."
					echo ".......Renaming Chromosomes..............."
					echo ".........................................."
					for local_file in $algorithm_dir
					do				
						[ -f "${local_file}" ] && ((files++))
						echo ""
						echo "Processing $local_file file..."	
						echo ""
						 NAME=`basename "$local_file"`
						 tmp=$(echo "$NAME" | awk -F '_' '{print $1}' )
						 newfname="${alg_dir}/${tmp}.pdb"
						 echo " New file name = ${newfname}"
						 #delete filename if exist
						if [ -f $newfname ] ; then
							echo "Deleting file = $newfname "
							rm $newfname 
						fi						 
						 mv $local_file $newfname
						 
						 						
						#copy the model out to the outside directory																					
						MatrixNAME="${directory}/${tmp}_matrix.txt"						
						echo $MatrixNAME						
						# Analysis									
						java -jar /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/pdbAnalyze_Batch2.jar $MatrixNAME $newfname $alg_dir
						#Change log_Name
						oldlogname="${alg_dir}/${tmp}_matrix_log.txt"	
						newlogname="${alg_dir}/${tmp}.log"	
						mv $oldlogname $newlogname	
						
					done
					
					: 'Multiline comment:
					#construct genome structure
					echo "Construct Genome Structure"
					echo "Processing ${genome_file} file..."	
					echo ""
					
					#Create algorithm genome_parameters	
					
					genomepath="${alg_dir}/genome" 
					rm -r $genomepath
					mkdir $genomepath
					#Generate a random parameter
					n=6
					rand=$(tr -cd '[:alnum:]' < /dev/urandom | head -c$n)
					param="genome_parameters_MOGEN_${rand}.txt"					
					
					{
					 echo "NUM = 1"
					 echo "NBR_OF_CHR = 25"
					 echo "INTRA_IF_THRESHOLD = 80%"
					 ehco "INTER_IF_THRESHOLD = 18%"
					 echo "ADJACENT_DIST = 1.5"
					 echo "CONTACT_DIST = 6.0"
					 echo "POS_MIN_DIST = 0.2"
					 echo "NEG_MAX_DIST_INTRA = 30"
					 echo "POS_MAX_DIST_WEIGHT_FILE = /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/parameters/pos_max_dist_weight_1_normal.txt"
					 echo "POS_MIN_DIST_WEIGHT_FILE = /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/parameters/pos_min_dist_weight_1_normal.txt"
					 echo "NEG_MIN_DIST_WEIGHT_FILE = /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/parameters/neg_min_dist_weight_1_normal.txt"
					 echo "NEG_MAX_DIST_WEIGHT_FILE = /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/parameters/neg_max_dist_weight_1_normal.txt"
					 echo ""
					 echo "OUTPUT_FOLDER = ${genomepath}/"
					 echo "INPUT_FILE = ${genome_file}"	
					 echo "CHR_UPPER_BOUND_ID_FILE = ${value}"
					 echo "LEARNING_RATE = 0.0001"
					 echo "VERBOSE = true"
					 echo "MAX_ITERATION = 200000"
					} > $param
					
					java -jar /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/MOGEN.jar $param
					rm $param			
					
					#Name file
					local_file="${genomepath}/*.pdb"
					echo "Processing $local_file file..."	
					echo ""
					 NAME=`basename "$local_file"`
					 tmp=$(echo "$NAME" | awk -F '_' '{print $1}' )
					 newfname="${genomepath}/${tmp}.pdb"
					 echo " New file name = ${newfname}"
					 #delete filename if exist
					if [ -f $newfname ] ; then
						echo "Deleting file = $newfname "
						rm $newfname 
					fi						 
					 mv $local_file $newfname
					
					#copy the model out to the outside directory																					
					MatrixNAME="${directory}/${tmp}_matrix.txt"						
					echo $MatrixNAME						
					# Analysis									
					java -jar /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/pdbAnalyze.jar $MatrixNAME $newfname $genomepath
					#Change log_Name
					oldlogname="${genomepath}/${tmp}_matrix_log.txt"	
					newlogname="${genomepath}/${tmp}.log"	
					mv $oldlogname $newlogname	
					
					echo "Number of files: ${files-0}"		
					'					

				fi
				
		done
	
	fi
done
	
echo ""

echo "LorDG Run Successfully!!!"

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