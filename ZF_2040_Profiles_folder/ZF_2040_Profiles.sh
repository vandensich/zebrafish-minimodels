#!/bin/bash

# Job name
#SBATCH --job-name=ZF_2040_Profiles
# Define format of output, deactivated
#SBATCH --output=ZF_2040_Profiles_%j-%a.out
# Define format of errorfile, deactivated
#SBATCH --error=ZF_2040_Profiles_%j-%a.err
# Define partition
#SBATCH --partition=single
# Define number of nodes per task
#SBATCH --nodes=1
# Define number of cores per node
#SBATCH --ntasks-per-node=64
# Define walltime
#SBATCH --time=12:00:00
# Define of repetition
#SBATCH -a 0-0


# Load R modules
module load math/R
export OPENBLAS_NUM_THREADS=64

# Run R script
Rscript ZF_2040_Profiles.R