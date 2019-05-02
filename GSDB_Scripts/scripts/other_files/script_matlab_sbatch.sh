#!/bin/bash
#-------------------------------------------------------------------------------
#  SBATCH CONFIG
#-------------------------------------------------------------------------------

#SBATCH --job-name=matlab_interactive
#SBATCH --mem=40G
#SBATCH --output="matlab_%j.out"
#SBATCH -p hpc4
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=oeow39@mail.missouri.edu
#SBATCH -t 2-00:00:00
#SBATCH --licenses=matlab:1     
#-------------------------------------------------------------------------------

echo "### Starting at: $(date) ###"

## Module Commands
module load matlab/matlab-R2018a
module list


## Run matlab non-interactively
SCRIPT='create_miniMDS_data.m'

srun matlab -nodesktop -nosplash -nodisplay -r "run('${SCRIPT}');exit"

echo "### Ending at: $(date) ###"