#!/bin/sh


# Determine the number folders/Directory
base_dir="/storage/htc/bdm/tosin/GSDB/SIMBA3D_Data"


for dir in $base_dir/*; do     
	if [ -d "$dir" ] 
	then
	namedir="${dir##*/}"    #Get the name	
	#echo $namedir
	NAME=`basename "$namedir"`
	tmp=$(echo "$NAME" | awk -F '_' '{print $2}' )	
	oldname="${base_dir}/${namedir}"	
	newname="${base_dir}/${tmp}"
	
	if [ -z "$tmp"  ];then
		echo "moving not necessary"
		
	else	
		echo $oldname
		echo $newname
		mv $oldname $newname
	fi
	
	fi
	
done
