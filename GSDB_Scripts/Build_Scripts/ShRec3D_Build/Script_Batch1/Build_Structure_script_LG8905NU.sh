#!/bin/sh



# Determine the number folders/Directory
base_dir="/storage/htc/bdm/tosin/GSDB/Data/ENCSR393LOP_LG8905NU/Extracted_Data/*"
output_dir="/storage/htc/bdm/tosin/GSDB/Structures/ENCSR393LOP_LG8905NU/"
algorithm="ShRec3D"
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
						echo "Processing ${local_file} file ........."	
						echo ""
						
						
						n=6
						rand=$(tr -cd '[:alnum:]' < /dev/urandom | head -c$n)
						fname=`basename "$local_file" .n_contact`
						param="local_ChromSDE_${rootname}_${namedir}_${fname}_${rand}.sh"	
						param="local.sh"
						
					
					{
						echo "#!/bin/bash -l"
						echo "#SBATCH -J ShRec3D_run"
						echo "#SBATCH -o ShRec3D_run-%j.out"
						echo "#SBATCH -p Lewis"
						echo "#SBATCH -N 1"
						echo "#SBATCH -n 1"
						echo "#SBATCH --mem 90G"
						echo "#SBATCH -t 2-00:00:00"
						echo "#SBATCH --licenses=matlab:1 "
						
						echo "n=6"
						echo "rand=\$(tr -cd '[:alnum:]' < /dev/urandom | head -c\$n)"
						echo "matlab_name=\"call_ShRec3D_\${rand}.m\""
						echo "call_script=\"/storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/ShRec3D/\${matlab_name}\";"
					
						echo "#Use Normalized local file before prediction"
						echo "echo """
						echo "echo \"Processing ${local_file} file .........\""	
						echo "echo """
						echo "NAME=\"`basename "$local_file" _matrix.txt`\""						
						echo "newfname=\"${alg_dir}/\${NAME}.pdb\""
						echo "outlog=\"${alg_dir}/\${NAME}.log\""
						
						echo "#Create algorithm parameters"
						echo "{"
						
						echo "echo \"filename='${local_file}';\""
						echo "echo \"outname='\${newfname}';\""
						echo "echo \"logout='\${outlog}';\" " 
						echo "echo \"ShRec3D\""
						
						echo "} > \$call_script"
						
						echo "## Module Commands"
						echo "module load matlab/matlab-R2018a"
						echo "module list"
						echo "# Run matlab non-interactively"
						echo "SCRIPT=\$call_script"

						echo "srun matlab -nodesktop -nosplash -nodisplay -r \"run('\${SCRIPT}');exit\""
										
						echo "rm \$call_script"					
						
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
					
					NAME=`basename "$genome_file" .txt`						
					newfname="${genomepath}/${NAME}.pdb"
					outlog="${genomepath}/${NAME}.log"
					{
					
					echo "filename='${genome_file}';"
					echo "outname='${newfname}';"
					echo "logout='${outlog}';"  
					echo "ShRec3D";
					
					} > $call_script
						
					## Module Commands
					module load matlab/matlab-R2018a
					module list
					# Run matlab non-interactively
					SCRIPT=$call_script

					srun matlab -nodesktop -nosplash -nodisplay -r "run('${SCRIPT}');exit"
					'
					
					echo "Number of files: ${files-0}"	
					

				fi
				
		done
	
	fi
done
	
echo ""

echo "ShRec3D Run Successfully!!!"

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