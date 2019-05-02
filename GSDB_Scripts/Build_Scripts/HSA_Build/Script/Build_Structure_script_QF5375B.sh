#!/bin/sh


# Determine the number folders/Directory
base_dir="/storage/htc/bdm/tosin/GSDB/Data/N_Data/QF5375B/*"
output_dir="/storage/htc/bdm/tosin/GSDB/Structures/QF5375B/"


algorithm="HSA"


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
							
					chr_path+="/chr*_HSA.txt"
					genome_file="${directory}/${rootname}_All_Genome_HSA.txt"
										
					echo "Chromosome Directory Name: $chr_path"	
					echo ""
					echo "Construct Chromosome Structure"								
					files=0	
						
					for local_file in $chr_path
					do	
						echo "Processing ${local_file}............."						
						[ -f "${local_file}" ] && ((files++))
												
						
						n=6
						rnd=$(tr -cd '[:alnum:]' < /dev/urandom | head -c$n)
						fnme=`basename "$local_file" _matrix.txt`
						#sriptname="local_HSA_${rootname}_${namedir}_${fnme}_${rnd}.sh"	
						
						sriptname="local_HSA.sh"
						{
							
							
							echo "#!/bin/bash -l"
							echo "#SBATCH -J HSA_run"
							echo "#SBATCH -o HSA_run-%j.out"
							echo "#SBATCH -p Lewis"
							echo "#SBATCH -N 1"
							echo "#SBATCH -n 8"
							echo "#SBATCH --mem 15G"
							echo "#SBATCH -t 2-00:00:00"
							
							echo "alg_dir=\"${path_to_output}/${algorithm}\""
							echo "NAME=\"`basename "$local_file" .txt`\""
							echo "outname=\"${alg_dir}/\${NAME}_out\""
							echo "#Generate a random parameter"
							echo "n=6"
							echo "rand=\$(tr -cd '[:alnum:]' < /dev/urandom | head -c\$n)"
							echo "param=\"parameter_HSA_\${rand}.txt\""
							
							echo "echo \"output file name PREFIX = \${outname}\""
							echo "{"
							echo "echo \"options('expressions'=50000)\""
							echo "echo \"source('/storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/HSA/cstruct1.R')\""
							echo "echo \"mat=read.table(\\\"${local_file}\\\",header=F,sep=\"\")\""
							echo "echo \"lsmap0=as.matrix(mat)\""
							echo "echo \"outfile=\\\"\${outname}\\\"\""
							echo "echo \"lscov0=0\""
							echo "echo \"mak=0\""
							echo "echo \"lsmap0=vector(\\\"list\\\",1)\""
							echo "echo \"lsmap0[[1]]=as.matrix(mat)\""
							echo "echo \"out=fmain(lsmap0,lscov0,outfile,300,150,50,50,0.001,0,0,mak)\""
							echo "} > \$param"
							echo "echo \"\""
							
							echo "#Call the HSA alg	"
							echo "module load  R/R-3.3.1"
							echo "/usr/bin/time --verbose R --no-save -f \$param " 
							
							echo "#convert to pdb and assessment"
							echo "echo \"pdb file conversion.....\""
							echo "java -jar /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/HSA/HSA2PDB.jar $local_file \"\${outname}.txt\" $alg_dir"
							
							echo "rm \$param "
							
							echo "#Change log_Name"
							echo "NAME=\"`basename "$local_file" _HSA.txt`\""
							echo "oldlogname=\"${alg_dir}/\${NAME}_HSA_log.txt\""	
							echo "newlogname=\"${alg_dir}/\${NAME}.log\""	
							echo "mv \$oldlogname \$newlogname	"	
							
							echo "#Change pdb_Name"
							echo "oldlogname=\"${alg_dir}/\${NAME}_HSA.pdb\""	
							echo "newlogname=\"${alg_dir}/\${NAME}.pdb\""	
							echo "mv \$oldlogname \$newlogname"	
							
						}	> $sriptname
						
						sbatch $sriptname
						
					done
					
					
					#construct genome structure
					: 'Multiline comment:
					echo "Construct Genome Structure"
					echo "Processing ${genome_file} file..."	
					echo ""
					
					#Create algorithm genome_parameters	
					
					genomepath="${alg_dir}/genome" 
					mkdir $genomepath
					
					NAME=`basename "$genome_file" .txt`
					outname="${genomepath}/${NAME}_out"
					param="genome_parameter_HSA.txt"
					echo "output file name PREFIX = ${outname}"
					{
					echo "options('expressions'=50000)"
					echo "source('/storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/HSA/cstruct1.R')"
					echo "mat=read.table(\"${genome_file}\",header=F,sep="")"
					echo "lsmap0=as.matrix(mat)"
					echo "outfile=\"${outname}\""
					echo "lscov0=0"
					echo "mak=0"
					echo "lsmap0=vector(\"list\",1)"
					echo "lsmap0[[1]]=as.matrix(mat)"
					echo "out=fmain(lsmap0,lscov0,outfile,300,150,50,50,0.001,0,0,mak)"
					} > $param
					echo ""
					
					#Call the HSA alg	
					module load  R/R-3.3.1
					/usr/bin/time --verbose R --no-save -f $param  
					
					echo "Genome Structure Completed successfully!!!"
					'
					
					echo "Number of files: ${files-0}"			

				fi
				
		done
	
	fi
done
	
echo ""

echo "HSA Run Successfully!!!"

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