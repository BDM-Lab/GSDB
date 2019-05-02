#!/bin/bash -l

# Determine the number folders/Directory
chr_path="/storage/htc/bdm/tosin/GSDB/Data/ENCSR499RVD_NR6150KF/Extracted_Data/GSE92811_ENCFF604YDD/VC_5kb"

				
for i in $(seq 1 25);
do					
		#Create path from ChromSDE data	
		n=6
		rand=$(tr -cd '[:alnum:]' < /dev/urandom | head -c$n)	
		param="local.sh"	
	{	
		
		echo "#!/bin/sh"
		echo "#!/bin/bash -l"
		echo "#SBATCH -J SIMBA3D_run"
		echo "#SBATCH -o SIMBA3D_run-%j.out"
		echo "#SBATCH -p Lewis"
		echo "#SBATCH -N 1"
		echo "#SBATCH -n 4"
		echo "#SBATCH --mem 80G"
		echo "#SBATCH -t 2-00:00:00"
		
							
		echo "n=6"
		echo "rand=\$(tr -cd '[:alnum:]' < /dev/urandom | head -c\$n)"
		echo "python_name=\"call_SIMBA3D_\${rand}.py\""
		echo "call_script=\"\${python_name}\";"
		
		if [ "$i" -lt 23 ]
		then
			echo "#Create algorithm parameters"
			
			echo "filename=\"${chr_path}/chr${i}_matrix.txt\""
			echo "echo \"Processing file \${filename}..............\""
			echo "new_filename=\"${chr_path}/chr${i}\""
			echo "{"				
			echo "echo \"import numpy as np\""
			echo "echo \"data=np.loadtxt('\$filename')\""
			echo "echo \"np.save('\$new_filename',data)\""
			echo "echo\"\""														
			echo "} > \$call_script"
			
		elif [ "$i" = 23 ]
		then
			
			echo "#Create algorithm parameters"
			echo "filename=\"${chr_path}/chrX_matrix.txt\""
			echo "echo \"Processing file \${filename}..............\""
			echo "new_filename=\"${chr_path}/chrX\""
			echo "{"				
			echo "echo \"import numpy as np\""
			echo "echo \"data=np.loadtxt('\$filename')\""
			echo "echo \"np.save('\$new_filename',data)\""
			echo "echo\"\""														
			echo "} > \$call_script"
		elif [ "$i" = 24 ]
		then
			echo "#Create algorithm parameters"
			echo "filename=\"${chr_path}/chrY_matrix.txt\""
			echo "echo \"Processing file \${filename}..............\""
			echo "new_filename=\"${chr_path}/chrY\""
			echo "{"				
			echo "echo \"import numpy as np\""
			echo "echo \"data=np.loadtxt('\$filename')\""
			echo "echo \"np.save('\$new_filename',data)\""
			echo "echo\"\""														
			echo "} > \$call_script"
		else
			echo "#Create algorithm parameters"
			echo "filename=\"${chr_path}/chrM_matrix.txt\""
			echo "echo \"Processing file \${filename}..............\""
			echo "new_filename=\"${chr_path}/chrM\""
			echo "{"				
			echo "echo \"import numpy as np\""
			echo "echo \"data=np.loadtxt('\$filename')\""
			echo "echo \"np.save('\$new_filename',data)\""
			echo "echo\"\""														
			echo "} > \$call_script"
		fi
		
			echo "python \$call_script"
			echo "echo \"Python Script Execution completed.............\""
							
			echo "rm \$call_script"					
			
	}> $param
	
	sbatch $param
	
done




echo "SIMBA3D Data created Run Successfully!!!"

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