#!/bin/sh

# Determine the number folders/Directory
base_dir="/storage/htc/bdm/tosin/GSDB/Data/ENCSR549MGQ_SU3340NA/Extracted_Data/*"
output_dir="/storage/htc/bdm/tosin/GSDB/Structures/SU3340NA/"

algorithm="SIMBA3D"

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
							
					chr_path+="/chr*.npy"
		
					
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
						fname=`basename "$local_file" .npy`
						param="local_SIMBA3D_${rootname}_${namedir}_${fname}.sh"	
						param="local.sh"
						
						{
							echo "#!/bin/bash -l"
							echo "#SBATCH -J SIMBA3D_run"
							echo "#SBATCH -o SIMBA3D_run-%j.out"
							echo "#SBATCH -p Lewis"
							echo "#SBATCH -N 1"
							echo "#SBATCH -n 8"
							echo "#SBATCH --mem 90G"
							echo "#SBATCH -t 2-00:00:00"							
						
							
							echo ""
							echo "#Generate a random parameter"
							echo "n=6"
							echo "rand=\$(tr -cd '[:alnum:]' < /dev/urandom | head -c\$n)"
							echo "param=\"parameters_SIMBA3D_\${rand}.json\""
							echo "call_script=\"get_SIMBA3D_XYZ_\${rand}.py\""
							
							echo "#Create algorithm parameters"
						
							echo "{"
							echo "echo \"[\""
							echo "echo \"	{\""
							echo "echo \"	\\\"taskname\\\":\\\"SIMBA3D_${fname}\\\",\""
							echo "echo \"	\\\"description\\\":\\\"json script for each chromosome construction\\\",\""
							echo "echo \"	\\\"uuid\\\":\\\"$rand\\\",\""
							echo "echo \"	\\\"file_names\\\": {\""
							echo "echo \"		\\\"inputdir\\\": \\\"${directory}/\\\",\""
							echo "echo \"		\\\"outputdir\\\": \\\"${alg_dir}/\\\",\""
							echo "echo \"		\\\"pairwise_contact_matrix\\\":\\\"${fname}.npy\\\",\""
							echo "echo \"		\\\"output_filename\\\" : \\\"${fname}.npz\\\"\""
							echo "echo \"		}\""							
							echo "echo \"	}\""
							echo "echo \" ]\""
							
							echo "} > \$param"
							
							echo "outname=\"${alg_dir}/${fname}.npz\""
							echo "# Call the SIMBA3D alg_parameters"	
							echo "matname=\"${alg_dir}/${fname}.mat\""
							echo "xyzname=\"${alg_dir}/${fname}.xyz\""
							echo "simba3d -r \$param"
							echo "rm \$param"
							
							echo "# convert npz to mat"	
							echo "simba3d-convertion \${outname} --ext_out .mat"
							
							echo "#Evaluate Result"							
							echo " echo \" Java Evaluation Started............\""
							echo " MatrixNAME=\"${directory}/${fname}_matrix.txt\""						
							echo " echo \$MatrixNAME"
							echo "newfname=\"${alg_dir}/${fname}.pdb\""
							echo "outlog=\"${alg_dir}/${fname}.log\""
							echo "matfile=\"${alg_dir}/${fname}.mat\""
							
							echo "{"						
							echo "echo \"import scipy.io as spio\""	
							echo "echo \"import numpy as np\""	
							echo "echo \"mat = spio.loadmat('\${matname}', squeeze_me=True)\""
							echo "echo \"M = mat['X_evol'] \""
							echo "echo \"G = np.transpose(M)\""							
							echo "echo \"np.savetxt('\${xyzname}', G,fmt='%5.5f' ,delimiter='\\t') \""
						
							echo "} > \$call_script"	
													
							echo "## Module Commands"
							echo "python \$call_script"
							
							echo "#Evaluate algorithm parameters"
							echo "newpdbname=\"${alg_dir}/${fname}.pdb\""
							echo "outlog=\"${alg_dir}/${fname}.log\""
							
							echo "# Analysis"									
						    echo "java -jar /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/XYZ2PDB.jar \$MatrixNAME \${xyzname} $alg_dir"
							
							echo "oldlogname=\"${alg_dir}/${fname}_matrix_log.txt\""	
							echo "oldpdbname=\"${alg_dir}/${fname}_matrix.pdb\"	"
							
							echo "mv \$oldlogname \$outlog"
							echo "mv \$oldpdbname \$newpdbname"	
							
							echo "rm \$call_script"		
						
								
						}>$param
					
						sbatch $param
				
					done

					
					echo "Number of files: ${files-0}"			

				fi
				
		done
	
	fi
done
	
echo ""

echo "SIMBA3D Run Successfully!!!"

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