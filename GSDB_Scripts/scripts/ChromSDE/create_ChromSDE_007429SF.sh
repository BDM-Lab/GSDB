#!/bin/bash -l

# Determine the number folders/Directory
base_dir="/storage/htc/bdm/tosin/GSDB/Data/OO7429SF/primary/GM12878_250kb/*"

for dir in $base_dir; do 
    [ -d "${dir}" ] && ((dir++))
	if [ -d "$dir" ] 
	then
		rootdir=$dir	
	    rootname="${rootdir##*/}"    #Get the name
		echo "Root Directory Name: $rootname "	
		rootdir+="/*"
				
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
					
			
					
				for i in $(seq 1 23);
				do					
						#Create path from ChromSDE data	
						n=6
						rand=$(tr -cd '[:alnum:]' < /dev/urandom | head -c$n)						
						param="local_ChromSDE_${rootname}_${i}_${rand}.sh"							
						param="local.sh"	
					{	
						
						echo "#!/bin/sh"
						
						echo "#!/bin/bash -l"
						echo "#SBATCH -J ChromSDE_run"
						echo "#SBATCH -o ChromSDE_run-%j.out"
						echo "#SBATCH -p Lewis"
						echo "#SBATCH -N 1"
						echo "#SBATCH -n 1"
						echo "#SBATCH --mem 80G"
						echo "#SBATCH -t 2-00:00:00"
						echo "#SBATCH --licenses=matlab:1 "
						
						echo "n=6"
						echo "rand=\$(tr -cd '[:alnum:]' < /dev/urandom | head -c\$n)"
						echo "matlab_name=\"call_ChromSDE_\${rand}.m\""
						echo "call_script=\"\${matlab_name}\";"
						
						if [ "$i" -lt 24 ]
						then
							echo "#Create algorithm parameters"
							echo "{"
							
							echo "echo \"filepath='${chr_path}';\""
							echo "echo \"output_path='${directory}';\""
							echo "echo \"filename = [filepath,'/chr',num2str($i),'_HSA.txt'];\""
							echo "echo \"disp (['Running Job for filename =', filename]);\""
							echo "echo \"contact = dlmread(filename);\""
							echo "echo \"nn=$i;\""
							echo "echo \"ChromSDE_input;\""							
							echo "} > \$call_script"							

						fi
						
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

echo "ChromSDE Data created Run Successfully!!!"

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