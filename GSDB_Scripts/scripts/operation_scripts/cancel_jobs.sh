#!/bin/bash
#-------------------------------------------------------------------------------
#  SBATCH CONFIG
#-------------------------------------------------------------------------------

for i in $(seq 7646132 7646216);
do

	scancel $i

done
