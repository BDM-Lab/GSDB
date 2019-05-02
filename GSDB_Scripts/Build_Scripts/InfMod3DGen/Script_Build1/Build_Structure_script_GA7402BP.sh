#!/bin/sh



# Determine the number folders/Directory
base_dir="/storage/htc/bdm/tosin/GSDB/Data/ENCSR489OCU_GA7402BP/Extracted_Data/*"
output_dir="/storage/htc/bdm/tosin/GSDB/Structures/ENCSR489OCU_GA7402BP/"

algorithm="InfMod3DGen"
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

					
					# Read the Genome sequesce file
					echo ""
					echo "Construct Chromosome Structure"				
					files=0	
					
					
						
				for i in $(seq 1 25);
				do					
						#Create path from ChromSDE data	
						n=6
						rand=$(tr -cd '[:alnum:]' < /dev/urandom | head -c$n)						
						param="local_InfMod3DGen_${rootname}_${i}_${rand}.sh"
						param="local.sh"	
					{	
						
						echo "#!/bin/sh"
						echo "#!/bin/bash -l"
						echo "#SBATCH -J InfMod3DGen_run"
						echo "#SBATCH -o InfMod3DGen_run-%j.out"
						echo "#SBATCH -p Lewis"
						echo "#SBATCH -N 1"
						echo "#SBATCH -n 1"
						echo "#SBATCH --mem 90G"
						echo "#SBATCH -t 2-00:00:00"
						echo "#SBATCH --licenses=matlab:1"
						
						echo ""
						echo "#Generate a random parameter"					
						echo "n=6"
						echo "rand=\$(tr -cd '[:alnum:]' < /dev/urandom | head -c\$n)"					
						echo "matlab_name=\"call_InfMod3DGen_\${rand}.m\""
						echo "call_script=\"/storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/InfMod3DGen/\${matlab_name}\";"
						
						if [ "$i" -lt 23 ]
						then
							echo "#Create algorithm parameters"
							
							echo "MatrixNAME=\"${chr_path}/chr${i}_matrix.txt\""
							
							echo "fname=\"${chr_path}/chr${i}_list.txt\""
							echo "new_fname=\"${chr_path}/chr${i}_InfMod3DGen.txt\""
							echo "outname=\"${alg_dir}/chr${i}\""
							
							echo "#Matlab for Execution and Evaluation"
							
							echo " echo \" MATLAB Evaluation Started............\""
							
							echo "newfname=\"${alg_dir}/chr${i}.pdb\""
							echo "outlog=\"${alg_dir}/chr${i}.log\""
							echo "matfile=\"${alg_dir}/chr${i}.mat\""
							
							echo "echo \"Processing file \${fname}..............\""
							
						elif [ "$i" = 23 ]
						then
							
							echo "#Create algorithm parameters"
							echo "MatrixNAME=\"${chr_path}/chrX_matrix.txt\""
													
							
							echo "fname=\"${chr_path}/chrX_list.txt\""
							echo "new_fname=\"${chr_path}/chrX_InfMod3DGen.txt\""
							echo "outname=\"${alg_dir}/chrX\""
							
							echo "#Matlab for Execution and Evaluation"
							
							echo " echo \" MATLAB Execution Started............\""
							
							echo "newfname=\"${alg_dir}/chrX.pdb\""
							echo "outlog=\"${alg_dir}/chrX.log\""
							echo "matfile=\"${alg_dir}/chrX.mat\""
							
							echo "echo \"Processing file \${fname}..............\""
							
						elif [ "$i" = 24 ]
						then
							echo "#Create algorithm parameters"
							echo "MatrixNAME=\"${chr_path}/chrY_matrix.txt\""
						
							echo "fname=\"${chr_path}/chrY_list.txt\""
							echo "new_fname=\"${chr_path}/chrY_InfMod3DGen.txt\""
							echo "outname=\"${alg_dir}/chrY\""
						
							
							echo "#Matlab for Execution and Evaluation"
							
							echo " echo \" MATLAB Execution Started............\""
							
							echo "newfname=\"${alg_dir}/chrY.pdb\""
							echo "outlog=\"${alg_dir}/chrY.log\""
							echo "matfile=\"${alg_dir}/chrY.mat\""
							
							echo "echo \"Processing file \${fname}..............\""
						
						else
							echo "#Create algorithm parameters"
							echo "MatrixNAME=\"${chr_path}/chrM_matrix.txt\""
							
							echo "fname=\"${chr_path}/chrM_list.txt\""
						
							echo "new_fname=\"${chr_path}/chrM_InfMod3DGen.txt\""
							echo "outname=\"${alg_dir}/chrM\""
							
							echo "#Matlab for Execution and Evaluation"
							
							echo " echo \" MATLAB Execution Started............\""
							
							echo "newfname=\"${alg_dir}/chrM.pdb\""
							echo "outlog=\"${alg_dir}/chrM.log\""
							echo "matfile=\"${alg_dir}/chrM.mat\""
							
							echo "echo \"Processing file \${fname}..............\""
							
						fi
						
							echo "{"	
							
							echo "echo \"ChrMod_main('\${fname}','\${new_fname}','\${matfile}',${i},1)\""
							echo "echo \"disp('InfMod3DGen Execution Completed........');\""	
							echo "echo \"data=load('\${matfile}');\""	
							echo "echo \"B=data.Ensemble;\""
							echo "echo \"[n,m,p]=size(B)\""
							echo "echo \"B = reshape(B,[n*m p]);\""
							
							echo "echo \"pdbout='\${newfname}';\""
							echo "echo \"pdblog='\${outlog}';\" " 
							echo ""
							
							echo "echo \"XYZ=B\""
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
										
							
					}> $param
					
					sbatch $param
					
				done
				
					
					
					

				
					

				fi
				
		done
	
	fi
done
	
echo ""

echo "InfMod3DGen Run Successfully!!!"

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