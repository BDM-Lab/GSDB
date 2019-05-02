#!/bin/bash



# Determine the number folders/Directory
base_dir="/storage/htc/bdm/tosin/GSDB/Data/ENCSR393LOP_LG8905NU"

output_dir="/storage/htc/bdm/tosin/GSDB/Unused_Data/"


echo "base directory = ${base_dir}"


#GSDB ID Level
for dir in $base_dir/Extracted_Data; do 
    [ -d "${dir}" ] && ((dir++))
	if [ -d "$dir" ] 
	then		
		rootdir=$dir
		rootname="${base_dir##*/}"    #Get the name
		echo ""
		echo ""
		echo "Root Directory Name: ${rootname} "
		rootdir+="/*"
		path_to_out=$output_dir$rootname  #Create path to output
		echo "Root Directory path_to_out: ${path_to_out} "
		mkdir $path_to_out 
		
		#nested for loop for new directory : Filename Level
		for directory in $rootdir; do
			[ -d "${directory}" ] && ((directories++))
			if [ -d "$directory" ] 
			then
				nextdir=$directory	
				namedir="${nextdir##*/}"    #Get the name					
				echo "Nest directory Name:  ${namedir}"					
				path_to_output="${path_to_out}/${namedir}"  #Create path to output
				
				echo "Root Directory path_to_output: ${path_to_output} "
				mkdir $path_to_output	
				
				#VC 10kb
				oldvc10="$directory/VC_10kb"
				newvc10="${path_to_output}/VC_10kb"
				mv $oldvc10 $newvc10
				
				#VC 5kb
				oldvc5="$directory/VC_5kb"
				newvc5="${path_to_output}/VC_5kb"
				mv $oldvc5 $newvc5

			
				


			fi	
		
		done
		
		
		
		
	fi
done
