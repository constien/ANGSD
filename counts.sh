#!/bin/bash

set -e

# Courtesy of StackOverflow
SCRIPT_DIR=$(cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd)

pushd Reads

featureCounts\
	-a ${SCRIPT_DIR}/T2TGenome/genomic.gtf\
	-o counts.txt\
	-T 4\
	*trim*.bam

popd
