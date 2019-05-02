#!/bin/sh

# Determine the number folders/Directory
base_dir="/storage/htc/bdm/tosin/GSDB/Data/ENCSR499RVD_NR6150KF/Extracted_Data/*"
output_dir="/storage/htc/bdm/tosin/GSDB/Structures/ENCSR499RVD_NR6150KF/"
algorithm="3DMax"

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
							
					chr_path+="/chr*_list.txt"
					genome_file="${directory}/${rootname}_All_Genome_list.txt"
					genome_length_file="${directory}/${rootname}_chrom_sequence_length.txt"
					value=$(<$genome_length_file) #Load length file
					echo "Length: ${value}"
					echo "Chromosome Directory Name: $chr_path"	
					
					# Read the Genome sequesce file
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
						param="local_3DMax_${rootname}_${namedir}_${fname}_${rand}.sh"	
						#param="local.sh"
						
						{
							echo "#!/bin/bash -l"
							echo "#SBATCH -J 3DMax_run"
							echo "#SBATCH -o 3DMax_run-%j.out"
							echo "#SBATCH -p Lewis"
							echo "#SBATCH -N 1"
							echo "#SBATCH -n 6"
							echo "#SBATCH --mem 80G"
							echo "#SBATCH -t 2-00:00:00"
							
							
							echo "#Generate a random parameter"
							echo "n=6"
							echo "rand=\$(tr -cd '[:alnum:]' < /dev/urandom | head -c\$n)"
							echo "param=\"parameters_3DMax_\${rand}.txt\""
							echo "#Create algorithm parameters"
							echo "{"
							echo "echo \"NUM = 1\""
							echo "echo \"OUTPUT_FOLDER = ${alg_dir}/\""
							echo " echo \"INPUT_FILE = ${local_file}\"	"					 
							echo "echo \"VERBOSE = false\""				 
							echo "echo \"LEARNING_RATE = 1.0\""
							echo "echo \"MAX_ITERATION = 20000\""
							echo "} > \$param"
							
							echo "# Call the 3DMax alg_parameters"							
							echo "java -jar /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/3DMax.jar \$param"
							echo "rm \$param	"
							
							echo "#Change log_Name"
							echo "NAME=\"`basename "$local_file" _list.txt`\""
							echo "oldlogname=\"${alg_dir}/\${NAME}_list_log.txt\""	
							echo "newlogname=\"${alg_dir}/\${NAME}.log\""	
							echo "mv \$oldlogname \$newlogname"
							
							
							echo "#Read all the pdb files in directory and print chromosome name only"
							echo "algorithm_dir=\"${alg_dir}/chr*_*.pdb\""
							
							echo "echo \"Renaming Chromosomes...............\""
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
					mkdir $genomepath
					limit=1
					for i in $(seq 0.1 0.2 $limit );do
						#Generate a random parameter
						n=6
						rand=$(tr -cd '[:alnum:]' < /dev/urandom | head -c$n)
						param="genome_parameters_3DMax_${rand}.txt"
						echo "Processing.... conversion factor = ${i}"
						{
						 echo "NUM = 1"
						 echo "OUTPUT_FOLDER = ${genomepath}/"
						 echo "INPUT_FILE = ${genome_file}"
						 echo "VERBOSE = true"
						 echo "CONVERT_FACTOR = ${i}"
						 echo "CHROMOSOME_LENGTH = ${value}"
						 echo "LEARNING_RATE = 1.0"
						 echo "MAX_ITERATION = 20000"
						} > $param		
						java -jar /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/3DMax.jar $param
						rm $param			
						
						#Change log_Name
						NAME=`basename "$genome_file" _list.txt`
						oldlogname="${genomepath}/${NAME}_list_log.txt"	
						newlogname="${genomepath}/${NAME}.log"	
						mv $oldlogname $newlogname						
					done
					
					'
					
					echo "Number of files: ${files-0}"			

				fi
				
		done
	
	fi
done
	
echo ""

echo "3DMax Run Successfully!!!"

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