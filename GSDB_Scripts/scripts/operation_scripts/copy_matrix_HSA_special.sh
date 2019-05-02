#!/bin/sh



# Determine the number folders/Directory
base_dir="/storage/htc/bdm/tosin/GSDB/Data/TE1402WS"
output_dir="/storage/htc/bdm/tosin/GSDB/HSA_Data/"

echo "base directory = ${base_dir}"


#GSDB ID Level
for dir in $base_dir; do 
    [ -d "${dir}" ] && ((dir++))
	if [ -d "$dir" ] 
	then		
		rootdir=$dir
		rootname="${dir##*/}"    #Get the name
		echo ""
		echo ""
		echo "Root Directory Name: ${rootname} "				
		path_to_out=$output_dir$rootname  #Create path to output
		echo "Root Directory path_to_out: ${path_to_out} "
		mkdir $path_to_out 
		
		subdir="${rootdir}/*"
		#nested for loop for new directory : Filename Level
		for directory in $subdir; do
			[ -d "${directory}" ] && ((directories++))
			if [ -d "$directory" ] 
			then
				nextdir=$directory	
				namedir="${nextdir##*/}"    #Get the name					
				echo "Nest directory Name:  ${namedir}"		
				path_to_output="${path_to_out}/${namedir}"  #Create path to output
				echo "path to output1: $path_to_output"
				mkdir $path_to_output	
				
				#nested for loop for new directory : Normalization Level
				for direct in $nextdir/*; do
					[ -d "${direct}" ] && ((directs++))
					if [ -d "$direct" ] 
					then
						directo=$direct
						namedir="${directo##*/}"   						
						echo "Nest directory Name:  ${namedir}"		
						path_to_output2="${path_to_output}/${namedir}"  #Create path to output
						echo "path to output2: $path_to_output2"
						mkdir $path_to_output2		
						
						#Matrix Directory
						
						echo "....Copy Started................."
						cur_matrix="$direct/chr*_HSA.txt"						
						cp $cur_matrix $path_to_output2						
						echo "....Copy Ended................."
						echo ""
						echo ""
					
					fi	
				
				done
			fi	
		
		done
		
	fi
done
