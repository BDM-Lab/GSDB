

# Determine the number folders/Directory
base_dir="/storage/htc/bdm/tosin/GSDB/Data/TE1402WS/*"
output_dir="/storage/htc/bdm/tosin/GSDB/Structures/TE1402WS/"	

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
					#matlab_name="call_Pastis_${rand}.m"
					#call_script="/storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/Pastis/${matlab_name}";
					
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
						#param="local_ChromSDE_${rootname}_${namedir}_${fname}_${rand}.sh"						
						param="local.sh"
						
					{
					
						echo "#!/bin/sh"						
						echo "#SBATCH -J pastis_run"
						echo "#SBATCH -o pastis_run-%j.out"
						echo "#SBATCH -p Lewis"
						echo "#SBATCH -N 1"
						echo "#SBATCH -n 1"
						echo "#SBATCH --mem 90G"
						echo "#SBATCH -t 2-00:00:00"
						echo "#SBATCH --licenses=matlab:1 "
												
						echo "echo \"structure construction started.....\""	
						echo "echo \"Processing ${local_file} file .........\""	
						echo "alg_dir=\"${path_to_output}/${algorithm}\""
						echo "NAME=\"`basename "$local_file" .n_contact`\""	
						
						echo "n=6"
						echo "rand=\$(tr -cd '[:alnum:]' < /dev/urandom | head -c\$n)"
						echo "matlab_name=\"call_Pastis_\${rand}.m\""
						echo "call_script=\"/storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/pastis/Evaluate/\${matlab_name}\";"
											
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
						echo "echo \"max_iter: 100\""
						echo "echo \"counts:\${NAME}.matrix\""
						echo "echo \"lengths: \${NAME}.bed\""
						echo "echo \"normalize: True\""
						echo "} > config.ini"
						
						echo "## Module Commands"
						#echo "pastis-mds ."
						echo " pastis-pm2 ."
						echo "echo \"Pastis structure generated Succesffully...............\""
						echo "cord=\"$(pwd)/\${dir_name}/PM2.structure\""
						echo "#Evaluate the result"
						#copy the model out to the outside directory																					
						echo " MatrixNAME=\"${directory}/\${NAME}_matrix.txt\""						
						echo " echo \$MatrixNAME"		
						
						echo "#Create algorithm parameters"
						echo "newfname=\"${alg_dir}/\${NAME}.pdb\""
						echo "outlog=\"${alg_dir}/\${NAME}.log\""
						
						echo "{"						
												
						echo "echo \"pdbout='\${newfname}';\""
						echo "echo \"pdblog='\${outlog}';\" " 
						echo ""
						echo "echo \"XYZ=dlmread('\${cord}');\""
						echo "echo \"Data2=dlmread('\$MatrixNAME');\""
						echo "echo \"fprintf('Previous Length = %d\\n',length(XYZ)); \""
						echo "echo \"distance;\""						
						echo "echo \"%output pdb and image  \""
						echo "echo \"output_structure;  \""
						echo "echo \" % correlate structure\""
						echo "echo \"[dist] = StructureDistance(XYZ,n); % distance from structure \""   
						echo "echo \"n=length(XYZ); \""
						echo "echo \"fprintf('New Length = %d\\n',n); \""
						echo "echo \"Data2(ntfd_index,:)=[];\""
						echo "echo \"Data2(:,ntfd_index)=[];\""
						echo "echo \"if (n~=length(Data2))\""
						echo "echo \"	disp('The Matrices are not of the same length');\""
						echo "echo \"end\""
						echo "echo \"%convert IF to distance \"" 
						echo "echo \" for j = 1:length(Data2)  \"" 
						echo "echo \"	  for k = j:n \"" 
						echo "echo \"	 	if(Data2(j,k)~=0  ) \""  
						echo "echo \"	 		Data(j,k)=1./(Data2(j,k));\""  
						echo "echo \"	 		Data(k,j)= Data(j,k);\""
						echo "echo \"	 	end\""
						echo "echo \"	  end\""				
						echo "echo \"end\""
						
						echo "echo \" [RHO,RHO2] = Spearman_corr(Data,dist,Data2,n);\""
						echo "echo \"fileID=fopen(pdblog,'w');  \""
						echo "echo \"fprintf(fileID, 'Spearman correlation Expected Dist vs. Reconstructed Dist: %f\n',RHO);;\""
						echo "echo \"fclose(fileID);\""
						echo "} > \$call_script"

						echo "## Module Commands"
						echo "module load matlab/matlab-R2018a"
						echo "module list"
						echo "# Run matlab non-interactively"
						echo "SCRIPT=\$call_script"

						echo "srun matlab -nodesktop -nosplash -nodisplay -r \"run('\${SCRIPT}');exit\""
										
						echo "rm \$call_script"		
									
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