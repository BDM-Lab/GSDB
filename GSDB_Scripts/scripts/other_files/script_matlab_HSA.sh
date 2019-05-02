#!/bin/bash
#-------------------------------------------------------------------------------
#  SBATCH CONFIG
#-------------------------------------------------------------------------------

#SBATCH --job-name=matlab_interactive
#SBATCH --mem=5G
#SBATCH --output="matlab_%j.out"
#SBATCH -p Interactive
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 0-04:00:00
#SBATCH --licenses=matlab:1     
#-------------------------------------------------------------------------------

echo "### Starting at: $(date) ###"

## Module Commands
module load matlab/matlab-R2018a
module list


## Run matlab non-interactively
SCRIPT='create_HSA_data.m'

srun matlab -nodesktop -nosplash -nodisplay -r "run('${SCRIPT}');exit"

echo "### Ending at: $(date) ###"