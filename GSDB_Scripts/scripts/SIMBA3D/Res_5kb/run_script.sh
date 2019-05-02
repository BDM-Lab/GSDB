#!/bin/sh


# Determine the number folders/Directory
base_dir="/storage/htc/bdm/tosin/GSDB_Scripts/scripts/SIMBA3D/Res_5kb"


for local_file in $base_dir/create*
do				
	[ -f "${local_file}" ] && ((files++))
	
	chmod +x $local_file
	$local_file
	
done
				