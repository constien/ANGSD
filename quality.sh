#!/bin/bash

set -e

pushd Reads

for FILE in *.fastq.gz; do
	conda run -n angsd fastqc ${FILE} &
done
wait

conda run -n multiqc multiqc .

mv multiqc_data untrimmed_multiqc_data
mv multiqc_report.html untrimmed_multiqc_report.html

for FILE in *.fastq.gz; do
	conda run -n trim-galore trim_galore --fastqc --stringency=3 ${FILE} &
done
wait

conda run -n multiqc multiqc *trim*

mv multiqc_data trimmed_multiqc_data
mv multiqc_report.html trimmed_multiqc_report.html

popd
