#!/bin/bash

set -e

pushd Reads

for FILE in *.bam; do
	samtools index ${FILE} &
done
wait

popd
