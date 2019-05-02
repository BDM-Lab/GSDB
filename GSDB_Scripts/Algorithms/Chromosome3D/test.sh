#!/bin/bash -l


for i in `seq 18 21`; do
	echo "Running job for Chromosome ${i} at 1MB.."
	./chromosome3D.pl -i "input/chr${i}_1mb_matrix.txt" -o "output/chr${i}_1mb" > "output/chr${i}_1mb.log" &
done

