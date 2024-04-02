#!/bin/bash

set -e 

curl -X GET -o assembly.zip https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_009914755.1/download?include_annotation_type=GENOME_FASTA,GENOME_GTF

yes | unzip assembly.zip

mv $(find ncbi_dataset -name "*.fna") .
mv $(find ncbi_dataset -name "*.gtf") .

rm assembly.zip
rm -r ncbi_dataset

STAR\
	--runMode genomeGenerate\
	--runThreadN 8\
	--genomeDir T2TGenome\
	--genomeFastaFiles *.fna\
	--sjdbGTFfile *.gtf\
	--sjdbOverhang 99

mv *.fna T2TGenome
mv *.gtf T2TGenome
