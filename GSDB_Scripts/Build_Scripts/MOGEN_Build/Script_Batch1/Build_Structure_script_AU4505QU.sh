#!/bin/sh



# Determine the number folders/Directory
base_dir="/storage/htc/bdm/tosin/GSDB/Data/ENCSR757IVO_AU4505QU/Extracted_Data/*"
output_dir="/storage/htc/bdm/tosin/GSDB/Structures/ENCSR757IVO_AU4505QU/"

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
					#remove output directory if exist
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
						
						n=6
						rand=$(tr -cd '[:alnum:]' < /dev/urandom | head -c$n)
						fname=`basename "$local_file" _list.txt`
						param="local_MOGEN_${rootname}_${namedir}_${fname}_${rand}.sh"	
						param="local.sh"
						
						{
							echo "#!/bin/bash -l"
							echo "#!/bin/bash -l"
							echo "#SBATCH -J MOGEN_run"
							echo "#SBATCH -o MOGEN_run-%j.out"
							echo "#SBATCH -p Lewis"
							echo "#SBATCH -N 1"
							echo "#SBATCH -n 6"
							echo "#SBATCH --mem 80G"
							echo "#SBATCH -t 2-00:00:00"
							
							echo "#Generate a random parameter"
							echo "n=6"
							echo "rand=\$(tr -cd '[:alnum:]' < /dev/urandom | head -c\$n)"
							echo "param=\"parameters_MOGEN_\${rand}.txt\""
							echo "#Create algorithm parameters"
							echo "{"
							echo " echo \"NUM = 1\""
							echo " echo \"NBR_OF_CHR = 1\""
							echo " echo \"INTRA_IF_THRESHOLD = 80%\""
							echo " echo \"ADJACENT_DIST = 1.5\""
							echo " echo \"CONTACT_DIST = 6.0\""
							echo " echo \"POS_MIN_DIST = 0.2\""
							echo " echo \"NEG_MAX_DIST_INTRA = 30\""
							echo " echo \"POS_MAX_DIST_WEIGHT_FILE = /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/parameters/pos_max_dist_weight_1_normal.txt\""
							echo " echo \"POS_MIN_DIST_WEIGHT_FILE = /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/parameters/pos_min_dist_weight_1_normal.txt\""
							echo " echo \"NEG_MIN_DIST_WEIGHT_FILE = /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/parameters/neg_min_dist_weight_1_normal.txt\""
							echo " echo \"NEG_MAX_DIST_WEIGHT_FILE = /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/parameters/neg_max_dist_weight_1_normal.txt\""
							echo " echo \"\""
							echo " echo \"OUTPUT_FOLDER = ${alg_dir}/\""
							echo " echo \"INPUT_FILE = ${local_file}\""
							echo " echo \"VERBOSE = false\"	"			 
							echo " echo \"LEARNING_RATE = 0.001\""
							echo " echo \"MAX_ITERATION = 200000\""
							 
							echo "} > \$param"
							
							echo "# Call the LorDG alg_parameters"	
							
							echo "java -jar /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/MOGEN.jar \$param"
							echo "echo \"construction completed... \""				
							echo "rm \$param	"
							
							echo "#Read all the pdb files in directory and print chromosome name only"
							echo "algorithm_dir=\"${alg_dir}/chr*_*.pdb\""
							echo "echo \"..........................................\""
							echo "echo \".......Renaming Chromosomes...............\""
							echo "echo \"..........................................\""
							echo "for local_file in \$algorithm_dir"
							echo "do"				
							echo "	[ -f \"\${local_file}\" ] && ((files++))"
							echo "	echo \"\""
							echo "	echo \"Processing \$local_file file...\"	"
							echo "	echo \"\""
							echo "	NAME=\`basename \"\$local_file\"\`"
							echo "	tmp=\$(echo \"\$NAME\" | awk -F '_' '{print \$1}' )"
							echo "	newfname=\"${alg_dir}/\${tmp}.pdb\""
							echo "	echo \" New file name = \${newfname}\""
							echo "	#delete filename if exist"
							echo "	if [ -f \$newfname ] ; then"
							echo "		echo \"Deleting old file = \$newfname \""
							echo "		rm \$newfname "
							echo "	fi"
							echo "	mv \$local_file \$newfname"
								 
														
								#copy the model out to the outside directory																					
							echo " MatrixNAME=\"${directory}/\${tmp}_matrix.txt\""						
							echo " echo $MatrixNAME"						
							echo " # Analysis	"							
							echo " java -jar /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/pdbAnalyze.jar \$MatrixNAME \$newfname $alg_dir"
							echo "	#Change log_Name "
							echo "	oldlogname=\"${alg_dir}/\${tmp}_matrix_log.txt\""
							echo "	newlogname=\"${alg_dir}/\${tmp}.log\""
							echo "	mv \$oldlogname \$newlogname"	
								
							echo "done"
					
		
						}>$param
					
						sbatch $param
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