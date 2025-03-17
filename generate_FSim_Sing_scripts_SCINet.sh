#!/bin/bash

# Usage: ./generate_FSim_Sing_scripts.sh <foa_run> <scenario> <run_timepoint> <num_parts> <study_area> <job_name>

if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <foa_run> <scenario> <run_timepoint> <num_parts> <study_area> <job_name>"
    exit 1
fi

foa_run=$1
scenario=$2
run_timepoint=$3
num_parts=$4
study_area=$5
job_name=$6

output_dir="/project/wildland_fire_smoke_tradeoff/${study_area}_${foa_run}_${scenario}_${run_timepoint}/_error"

for part in $(seq 1 $num_parts); do
    script_name="RunFSim_${foa_run}_pt${part}_${scenario}_${run_timepoint}.sh"
    slurm_script="RunFSim_Sing_${foa_run}_pt${part}_${scenario}_${run_timepoint}.sh"
    
    cat > "$slurm_script" <<EOL
#!/bin/bash

#SBATCH --job-name=${job_name}
#SBATCH --partition=ceres
#SBATCH --time=07-00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
#SBATCH --mail-user=laurel.sindewald@usda.gov
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=${output_dir}/${foa_run}_pt${part}_${scenario}_${run_timepoint}_slurmlog_%A.out
#SBATCH --error=${output_dir}/${foa_run}_pt${part}_${scenario}_${run_timepoint}_error_%A.log

module load apptainer

singularity exec /project/wildland_fire_smoke_tradeoff/firemodels.sif ./$script_name
EOL

done

echo "Generated $num_parts SLURM scripts for ${study_area} ${foa_run}, scenario ${scenario}, run timepoint ${run_timepoint}."
