#!/bin/sh


# Determine the number folders/Directory
base_dir="/storage/htc/bdm/tosin/GSDB_Scripts/Build_Scripts/ChromSDE/Scripts_Batch1"


for local_file in $base_dir/Build*
do				
	[ -f "${local_file}" ] && ((files++))
	
	chmod +x $local_file
	$local_file
	
done
				