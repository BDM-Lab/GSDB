#!/bin/sh

# Determine the number folders/Directory
base_dir="/storage/htc/bdm/tosin/GSDB_Scripts/scripts/ChromSDE"


for local_file in $base_dir/create_ChromSDE_*
do				
	[ -f "${local_file}" ] && ((files++))
	
	chmod +x $local_file
	echo $local_file
	$local_file
done
				