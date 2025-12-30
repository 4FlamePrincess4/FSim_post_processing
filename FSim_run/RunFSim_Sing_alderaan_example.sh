#!/bin/bash

## This is code to run the script that launches FSim within the firemodels.sif container

#SBATCH --job-name=F1_20_p1
#SBATCH --partition=math-alderaan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=laurel.sindewald@ucdenver.edu
#SBATCH --output=FOA1c_r16_pt1_LF2020_RecOff_slurmlog_%A.out

singularity exec firemodels.sif /data001/projects/sindewal/RunFSim_foa1c_r16_pt1_LF2020.sh

