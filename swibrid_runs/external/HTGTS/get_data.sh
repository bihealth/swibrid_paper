#!/bin/bash -e

mkdir -p input tmp
for acc in SRR21047{31..47} SRR62934{56..63}; do 
        prefetch ${acc} --output-directory tmp
	fasterq-dump -p --split-3 -O tmp tmp/${acc}/${acc}.sra
	bbmerge.sh in1=tmp/${acc}_1.fastq in2=tmp/${acc}_2.fastq out=input/${acc}.fastq.gz
	python make_info.py -i input/${acc}.fastq.gz -o input/${acc}_info.csv
done
