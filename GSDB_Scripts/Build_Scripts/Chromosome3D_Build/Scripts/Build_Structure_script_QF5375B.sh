


# Determine the number folders/Directory


base_dir="/storage/htc/bdm/tosin/GSDB/Data/ENCSR323QIP_QF5375B/Extracted_Data/*"
output_dir="/storage/htc/bdm/tosin/GSDB/Structures/ENCSR323QIP_QF5375B/"

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
						
						n=6
						rand=$(tr -cd '[:alnum:]' < /dev/urandom | head -c$n)
						fname=`basename "$local_file" _matrix.txt`
						param="local_Chromosome3D_${rootname}_${namedir}_${fname}_${rand}.sh"	
						
						{		
							echo "#!/bin/sh"

							echo "#!/bin/bash -l"
							echo "#SBATCH -J Chromosome3D_run"
							echo "#SBATCH -o Chromosome3D_run-%j.out"
							echo "#SBATCH -p Lewis"
							echo "#SBATCH -N 1"
							echo "#SBATCH -n 5"
							echo "#SBATCH --mem 20G"
							echo "#SBATCH -t 2-00:00:00"
							
							#Call the miniMDS alg_parameters
							echo "echo \"structure construction started.....\""	
							echo "alg_dir=\"${path_to_output}/${algorithm}\""
							echo "NAME=\"`basename "$local_file" .txt`\""
							echo "outname=\"${alg_dir}/\${NAME}\""
							echo "echo \"\"	"
							echo "/storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/Chromosome3D/chromosome3D.pl -i $local_file -o \$outname "
							
							#copy the model out to the outside directory
							echo "pdbfile=\"\${alg_dir}/\${NAME}/\${NAME}_model1.pdb\""
							echo "NAME=\`basename \"\$pdbfile\" _matrix_model1.pdb\`"						
							echo "pdbNAME=\"\${NAME}.pdb\""
							echo "newpdbfile=\"${alg_dir}/\${pdbNAME}\""
							echo "echo \$newpdbfile"
							echo "cp \$pdbfile \$newpdbfile"
							
							# Analysis									
							echo "java -jar /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/pdbAnalyze.jar $local_file \$newpdbfile $alg_dir"		
							#Change log_Name
							echo "oldlogname=\"${alg_dir}/\${NAME}_matrix_log.txt\""	
							echo "newlogname=\"${alg_dir}/\${NAME}.log\""	
							echo "mv \$oldlogname \$newlogname"
							echo "echo \"structure construction and analysis completed.....\""	
						}> $param
						
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

echo "chromosome3D Run Successfully!!!"

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