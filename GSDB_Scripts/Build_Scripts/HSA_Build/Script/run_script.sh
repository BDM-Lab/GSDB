#!/bin/sh

# Determine the number folders/Directory
base_dir="/storage/htc/bdm/tosin/GSDB_Scripts/Build_Scripts/HSA_Build/Script"


for local_file in $base_dir/Build*
do				
	[ -f "${local_file}" ] && ((files++))
	
	dos2unix $local_file
	
	#chmod +x $local_file
	
	#echo "$local_file"
	$local_file
	
done
				