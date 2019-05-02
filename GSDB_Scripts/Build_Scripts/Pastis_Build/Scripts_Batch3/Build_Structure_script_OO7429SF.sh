#!/bin/bash

# Determine the number folders/Directory
base_dir="/storage/htc/bdm/tosin/GSDB/Data/OO7429SF/primary/GM12878_250kb/*"
output_dir="/storage/htc/bdm/tosin/GSDB/Structures/OO7429SF/"
	

algorithm="Pastis"
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
							
					chr_path+="/chr*.n_contact"
					
					
					# Read the chromosome sequesce file
					echo ""
					echo "Construct Chromosome Structure"				
					files=0	
										
					n=6
					rand=$(tr -cd '[:alnum:]' < /dev/urandom | head -c$n)
				
					
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
						logname="${fname}_${rootname}_${namedir}.log"
						param="${fname}_${rootname}_${namedir}.sh"						
						param="local.sh"
						
					{
					
						echo "#!/bin/sh"						
						echo "#SBATCH -J pastis_run"
						echo "#SBATCH -o pastis_run-%j.out"
						echo "#SBATCH -p Lewis"
						echo "#SBATCH -N 1"
						echo "#SBATCH -n 8"
						echo "#SBATCH --mem 80G"
						echo "#SBATCH -t 2-00:00:00"
											
												
						echo "echo \"structure construction started.....\""	
						echo "echo \"Processing ${local_file} file .........\""	
						echo "alg_dir=\"${path_to_output}/${algorithm}\""
						echo "NAME=\"`basename "$local_file" .n_contact`\""	
						
						echo "n=6"
						echo "rand=\$(tr -cd '[:alnum:]' < /dev/urandom | head -c\$n)"
					
																
						echo "dir_name=\"\${NAME}_\${rand}\""
						echo "mkdir \$dir_name"
						echo "cbins=\"${directory}/\${NAME}.cbins \""
						echo "counts=\"\$dir_name/\${NAME}.matrix\""
						echo "lenght=\"\$dir_name/\${NAME}.bed\""
						
						echo "cp $local_file \$counts"
						echo "cp \$cbins \$lenght"
						
						echo "cd \$dir_name"	
						echo "pwd"
						#Create algorithm parameters
						echo "{"
						echo "echo \"[all]\""
						echo "echo \"output_name: structure\""
						echo "echo \"verbose: 1\""						
						echo "echo \"counts:\${NAME}.matrix\""
						echo "echo \"lengths: \${NAME}.bed\""
						echo "echo \"normalize: False\""
						echo "} > config.ini"
						
						echo "## Module Commands"
						#echo "pastis-mds ."
						echo " pastis-pm2 ."
						echo "echo \"Pastis structure generated Succesffully...............\""
						echo "pdbfile=\"$(pwd)/\${dir_name}/PM2.structure.pdb\""
						echo "#Evaluate the result"
						#copy the model out to the outside directory																					
						echo " MatrixNAME=\"${directory}/\${NAME}_matrix.txt\""						
						echo " echo \$MatrixNAME"		
						
						echo "#Create algorithm parameters"
						echo "newfname=\"${alg_dir}/\${NAME}.pdb\""
						echo "outlog=\"${alg_dir}/\${NAME}.log\""
						
						echo "mv \$pdbfile \$newfname"
						echo "# Analysis"									
						echo "java -jar /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/pdbAnalyze.jar \$MatrixNAME \$newfname \$alg_dir"
						echo "oldlogname=\"\${alg_dir}/\${NAME}_matrix_log.txt\""	
						echo "newlogname=\"\${alg_dir}/\${NAME}.log\"	"
						echo "mv \$oldlogname \$newlogname"
						
						echo "rm -r \$(pwd)"		
					}> $param					
								
						sbatch $param
					done
							
					echo "Number of files: ${files-0}"	
					

				fi
				
		done
	
	fi
done
	
echo ""

echo "Pastis Run Successfully!!!"

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