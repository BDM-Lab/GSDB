#!/bin/sh


# Determine the number folders/Directory

path="/storage/htc/bdm/tosin/GSDB/Data/OO7429SF/primary/GM12878_normalized/GM12878/KR_10kb"

chr_path="${path}/chr*_10kb.KRnorm.txt"

for local_file in $chr_path
do				
	[ -f "${local_file}" ] && ((files++))
	
	#Use Normalized local file before prediction
	echo ""
	echo "Processing $local_file file... "	
	echo ""
	
	fname=`basename "$local_file" _10kb.KRnorm.txt`
	newname="${path}/${fname}_list.txt"
	echo "New file name... ${newname}"	
	mv $local_file $newname
	
done

	
chr_path="${path}/chr*_10kb.KRnormSQUARE.txt"

for local_file in $chr_path
do				
	[ -f "${local_file}" ] && ((files++))
	
	#Use Normalized local file before prediction
	echo ""
	echo "Processing $local_file file... "	
	echo ""
	
	fname=`basename "$local_file" _10kb.KRnormSQUARE.txt`
	newname="${path}/${fname}_matrix.txt"
	echo "New file name... ${newname}"	
	mv $local_file $newname
	
done