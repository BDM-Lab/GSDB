#!/bin/sh

#!/bin/bash -l
#SBATCH -J Chrom3D_Remove_run
#SBATCH -o Chrom3D_Remove_run-%j.out
#SBATCH -p Interactive
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu 1G
#SBATCH -t 0-04:00:00


# Determine the number folders/Directory
base_dir="/storage/htc/bdm/tosin/GSDB/Structures/*"

algorithm="Chromosome3D"

echo "base directory = ${base_dir}"


#GSDB ID Level
for dir in $base_dir; do 
    [ -d "${dir}" ] && ((dir++))
	if [ -d "$dir" ] 
	then		
		rootdir=$dir
		rootname="${base_dir##*/}"    #Get the name
		echo ""
		echo ""
		echo "Root Directory Name: ${rootname} "
		rootdir+="/*"				
		
		#nested for loop for new directory : Filename Level
		for directory in $rootdir; do
			[ -d "${directory}" ] && ((directories++))
			if [ -d "$directory" ] 
			then
				nextdir=$directory	
				namedir="${nextdir##*/}"    #Get the name					
				echo "Nest directory Name:  ${namedir}"		
				#nested for loop for new directory : Normalization Level
				for direct in $nextdir/*; do
					[ -d "${direct}" ] && ((directs++))
					if [ -d "$direct" ] 
					then
						directo=$direct
						namedir="${directo##*/}"   						
						echo "Nest directory Name:  ${namedir}"		
												
						#Algorithm Directory
						cur_alg_dir="$direct/$algorithm"					
						echo $cur_alg_dir 
						echo "....Removal Started................."
						rm -r  $cur_alg_dir 
						echo "....Removal Ended................."
						echo ""
						echo ""
					
					fi	
				
				done
			fi	
		
		done
		
	fi
done
