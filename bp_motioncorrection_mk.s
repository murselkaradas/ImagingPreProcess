#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=05:00:00
#SBATCH --mem=60GB
#SBATCH --job-name=motionCorrectionTest
#SBATCH --mail-type=END
#SBATCH --mail-user=karadm01@nyumc.org
#SBATCH --array=0-250
module purge
module load matlab/R2018a
cd /gpfs/scratch/karadm01/2Pdata/SC18/220908/stim2
tifs=(/gpfs/scratch/karadm01/2Pdata/SC18/220908/stim2/SC18_220908_field1stim*.tif)

tif=${tifs[$SLURM_ARRAY_TASK_ID]}
{
    echo $tif
    matlab -nodisplay -r "normcorremotioncorrection_single('$tif','Ref.tif'); exit"
} > $tif.log 2>&1

exit