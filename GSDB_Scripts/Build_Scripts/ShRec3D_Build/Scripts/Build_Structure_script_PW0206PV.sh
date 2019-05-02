#!/bin/sh

#!/bin/bash -l
#SBATCH -J ShRec3D_run
#SBATCH -o ShRec3D_run-%j.out
#SBATCH -p Lewis
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 80G
#SBATCH --mail-type=end
#SBATCH -t 2-00:00:00
#SBATCH --licenses=matlab:1 

# Determine the number folders/Directory
base_dir="/storage/htc/bdm/tosin/GSDB/Data/ENCSR079VIJ_PW0206PV/Extracted_Data/*"
output_dir="/storage/htc/bdm/tosin/GSDB/Structures/ENCSR079VIJ_PW0206PV/"

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
					
					
					n=6
					rand=$(tr -cd '[:alnum:]' < /dev/urandom | head -c$n)
					matlab_name="call_ShRec3D_${rand}.m"
					call_script="/storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/ShRec3D/${matlab_name}";
					
					for local_file in $chr_path
					do				
						[ -f "${local_file}" ] && ((files++))
						
						#Use Normalized local file before prediction
						echo ""
						echo "Processing ${local_file} file ........."	
						echo ""
						NAME=`basename "$local_file" _matrix.txt`						
						newfname="${alg_dir}/${NAME}.pdb"
						outlog="${alg_dir}/${NAME}.log"
						
						#Create algorithm parameters
						{
						
						echo "filename='${local_file}';"
						echo "outname='${newfname}';"
						echo "logout='${outlog}';"  
						echo "ShRec3D"
						
						} > $call_script
						
						## Module Commands
						module load matlab/matlab-R2018a
						module list
						# Run matlab non-interactively
						SCRIPT=$call_script

						srun matlab -nodesktop -nosplash -nodisplay -r "run('${SCRIPT}');exit"
					
						
									
					done
					
					rm $call_script
					

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