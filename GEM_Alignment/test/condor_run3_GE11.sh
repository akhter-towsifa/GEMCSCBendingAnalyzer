#!/bin/bash

set -euo pipefail

scratch_dir="$PWD"
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/t/toakhter/private/mual_tamu/CMSSW_16_0_4/src
eval "$(scramv1 runtime -sh)"
cd /afs/cern.ch/work/t/toakhter/private/mual_tamu/CMSSW_16_0_4/src/GEMCSCBendingAnalyzer/GEM_Alignment/test/

job_index="${1:?usage: $0 <condor-process-index>}"
line_number=$((job_index + 1))
input_file=$(sed -n "${line_number}p" file_names.list)

if [[ -z "${input_file}" ]]; then
    echo "No input file found for job index ${job_index}"
    exit 1
fi

echo "Processing file: ${input_file} with output number: ${line_number}"
cmsRun condor_run3_GE11.py inputFile="${input_file}" outputFileNumber="${line_number}"