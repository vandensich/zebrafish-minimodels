#!/bin/bash

# Job name
#SBATCH --job-name=ZebraFish2040_TG_all_Hill_more
# Define format of output, deactivated
#SBATCH --output=ZebraFish2040_TG_all_Hill_more_%j-%a.out
# Define format of errorfile, deactivated
#SBATCH --error=ZebraFish2040_TG_all_Hill_more_%j-%a.err
# Define partition
#SBATCH --partition=single
# Define number of nodes per task
#SBATCH --nodes=1
# Define number of cores per node
#SBATCH --ntasks-per-node=64
# Define walltime
#SBATCH --time=120:00:00
# Define of repetition
#SBATCH -a 0-56


# Load R modules
module load math/R
export OPENBLAS_NUM_THREADS=64

# Run R script
Rscript ZebraFish2040_TG_all_Hill_more.R