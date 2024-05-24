#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --time=05:00:00
#SBATCH --mem=128GB
#SBATCH --job-name=OdorResponses
#SBATCH --mail-type=END
#SBATCH --mail-user=karadm01@nyumc.org
module purge
module load matlab/R2023a
cd /gpfs/data/rinberglab/Mursel/ImagingPreProcess
matlab -nodisplay -r "RUN_stim_prepocess" -logfile BarrelCortex_200520.log
