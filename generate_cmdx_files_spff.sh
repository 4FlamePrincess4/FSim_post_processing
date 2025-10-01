#!/bin/bash

# Usage: ./generate_cmdx_files.sh <pyrome> <GCM> <timepoint> <num_parts> <study_area> <lcp> <fdist> <erc> <idg> <frisk> <adj>

if [ "$#" -ne 9 ]; then
    echo "Usage: $0 <pyrome> <GCM> <timepoint> <num_parts> <study_area> <lcp> <fdist> <erc> <idg> <frisk> <adj>"
    exit 1
fi

pyrome=$1
GCM=$2
timepoint=$3
num_parts=$4
study_area=$5
lcp=$6
fdist=$7
erc=$8
idg=$9
frisk=$10
adj=$11

output_dir="/project/spff/${study_area}_${pyrome}_${GCM}_${timepoint}"

for part in $(seq 1 $num_parts); do
    cmdx_file="${pyrome}_pt${part}_${GCM}_${timepoint}.cmdx"
    
    cat > "$cmdx_file" <<EOL
IgnitionProbabilityGrid:      _inputs/idg/${idg}
landscape:                    _inputs/lcp/${lcp}
FireDayDistributionFile:      _inputs/fdist/${fdist}
ROSAdjust:                    _inputs/adj/${adj}
FMS80:                        _inputs/fms/fms80.fms
FMS90:                        _inputs/fms/fms90.fms
FMS97:                        _inputs/fms/fms97.fms
SeasonErcFile:                _inputs/erc/${erc}_SeasonERC${part}.csv
FriskFile:                    _inputs/frisk/{frisk}
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
OutputsName:                  ${pyrome}_pt${part}_${GCM}_${timepoint}
ProgressFilePathname:	      ./${pyrome}_pt${part}_${GCM}_${timepoint}_progress.txt

EOL

done

echo "Generated $num_parts cmdx files for ${study_area} ${pyrome}, scenario ${GCM}, run timepoint ${timepoint}."
