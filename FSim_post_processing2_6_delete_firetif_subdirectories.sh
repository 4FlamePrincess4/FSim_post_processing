#!/bin/bash
# Note: Before running this on Linux, you'll need to run this line to remove stupid Windows characters:
# sed -i -e 's/\r$//' rename_tifs.sh

#SBATCH --job-name=120BOcln
#SBATCH --account=wildland_fire_smoke_tradeoff
#SBATCH --partition=ceres
#SBATCH --time=02-00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
#SBATCH --output=/project/wildland_fire_smoke_tradeoff/okawen_foa1c_r17_LF2020_TM_baseline_time0/_error/foa1c_r17_LF2020_TM_baseline_time0_cleanup_slurmlog_%A.out
#SBATCH --error=/project/wildland_fire_smoke_tradeoff/okawen_foa1c_r17_LF2020_TM_baseline_time0/_error/foa1c_r17_LF2020_TM_baseline_time0_cleanup_error_%A.log
#SBATCH --mail-user=laurel.sindewald@usda.gov
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

BASE_DIR="/project/wildland_fire_smoke_tradeoff/okawen_foa1c_r17_LF2020_TM_baseline_time0"   # Replace with the base directory 
ERROR_FILE="${SLURM_SUBMIT_DIR}/slurm-${SLURM_JOB_ID}.err"  # SLURM error file path

delete_dir() {
  local d="$1"
  local max_tries=5
  local try=1

  while [[ -d "$d" && $try -le $max_tries ]]; do
    rm -rf "$d" 2>/dev/null || true
    sleep 2
    ((try++))
  done

  if [[ -d "$d" ]]; then
    echo "ERROR: Could not fully delete $d after $max_tries attempts" >&2
    return 1
  fi
}

set -euo pipefail

for dir in "$BASE_DIR"/*{ArrivalDays,FlameLengths,ArrivalTimes,gdb}; do
  # Skip if glob didn't match anything
  [[ -d "$dir" ]] || continue

  zipfile="${dir}.zip"

  if [[ -f "$zipfile" ]]; then
    echo "Deleting directory: $dir (zip found: $zipfile)"
    delete_dir "$dir" || echo "WARNING: Incomplete delete: $dir" >&2
  else
    echo "No zip found for: $dir"
  fi

done





