#!/bin/bash
# Note: Before running this on Linux, you'll need to run this line to remove stupid Windows characters:
# sed -i -e 's/\r$//' rename_tifs.sh

#SBATCH --job-name=120BOzip
#SBATCH --account=wildland_fire_smoke_tradeoff
#SBATCH --partition=ceres
#SBATCH --time=02-00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
#SBATCH --output=/project/wildland_fire_smoke_tradeoff/okawen_foa1c_r17_LF2020_TM_baseline_time0/_error/foa1c_r17_LF2020_TM_baseline_time0_zip_tifs_slurmlog_%A.out
#SBATCH --error=/project/wildland_fire_smoke_tradeoff/okawen_foa1c_r17_LF2020_TM_baseline_time0/_error/foa1c_r17_LF2020_TM_baseline_time0_zip_tifs_error_%A.log
#SBATCH --mail-user=laurel.sindewald@usda.gov
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

set -euo pipefail

for dir in *_FlameLengths *_ArrivalDays *_ArrivalTimes; do
  # Skip if glob didn't match anything
  [[ -d "$dir" ]] || continue

  zipfile="${dir}.zip"

  if [[ -f "$zipfile" ]]; then
    echo "[DRY RUN] Would delete directory: $dir (zip found: $zipfile)"
  else
    echo "[SKIP] No zip found for: $dir"
  fi
done