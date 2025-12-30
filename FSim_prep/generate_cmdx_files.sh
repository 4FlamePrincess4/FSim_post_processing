#!/bin/bash

# Usage: ./generate_cmdx_files.sh <foa> <foa_run> <scenario> <run_timepoint> <num_parts> <study_area> <lcp> <fdist> <erc>

if [ "$#" -ne 9 ]; then
    echo "Usage: $0 <foa> <foa_run> <scenario> <run_timepoint> <num_parts> <study_area> <lcp> <fdist> <erc>"
    exit 1
fi

foa=$1
foa_run=$2
scenario=$3
run_timepoint=$4
num_parts=$5
study_area=$6
lcp=$7
fdist=$8
erc=$9

output_dir="/project/wildland_fire_smoke_tradeoff/${study_area}_${foa_run}_${scenario}_${run_timepoint}"

for part in $(seq 1 $num_parts); do
    cmdx_file="${foa_run}_pt${part}_${scenario}_${run_timepoint}.cmdx"
    
    cat > "$cmdx_file" <<EOL
IgnitionProbabilityGrid:      _inputs/idg/${foa}_idg.tif
landscape:                    _inputs/lcp/${lcp}
FireDayDistributionFile:      _inputs/fdist/${fdist}.fdist
ROSAdjust:                    _inputs/adj/${foa_run}.adj
FMS80:                        _inputs/fms/fms80.fms
FMS90:                        _inputs/fms/fms90.fms
FMS97:                        _inputs/fms/fms97.fms
SeasonErcFile:                _inputs/erc/${erc}_SeasonERC${part}.csv
FriskFile:                    _inputs/frisk/Colville_FOA4d_r1.frisk
#Customfmd:                    _inputs/fmd/XXX.fmd (if present)
#IgnitionMask:                XXX_mask.tif (if present)
#PolicyGrid:                  XXX_policy.tif (if present)
#PolicyCSV:                   XXX_policy.csv (if present)
FireSizeLimit:                -1
DisableConvergenceAngleCheck: 1
PowerLawMaxFires:             0
CrownFireMethod:              1
JulianStart:                  1
SortFires:                    1
FireListOnly:                 0
#
OutputFirePerims:             1
OutputIgnitions:              1
OutputFlameLengths:           1
OutputArrivalDays:            1
OutputArrivalTimes:           1
EmberOutputs:                 1
#
Suppression:                  1
SuppressionFactor:            2.0
#
NumSimulations:               5000
#
OutputsName:                  ${foa_run}_pt${part}_${scenario}_${run_timepoint}
ProgressFilePathname:	      ./${foa_run}_pt${part}_${scenario}_${run_timepoint}_progress.txt

EOL

done

echo "Generated $num_parts cmdx files for ${study_area} ${foa_run}, scenario ${scenario}, run timepoint ${run_timepoint}."
