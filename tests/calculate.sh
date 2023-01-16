#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

module load jq/1.6

# enter the workflow's final output directory ($1)
cd $1

# bamqc is a one-line JSON, making diff a horror. This pretty prints it
jq '.' bamQC.bamQC_results.json
# These are the merged fastp outputs from the WDL scatter. Sort in case scatter order is random
# Then drop the command field, as it will be different for each run
jq '. |= sort_by(.command)' output.fastp.json | jq '.[] |= del(.command)'
grep -v ^# output.lambda.controlstats
zcat output.methyldackel_CpG.bedGraph.gz
grep -v ^# output.puc19.controlstats
