#!/bin/bash 

mkdir -p demux
cat index/barcodes_primers.fa | makeblastdb -in - -title "barcodes_primers" -out demux/barcodes_primers -input_type fasta -dbtype nucl

for sample in $(grep SAMPLES config.yaml | tr ',' '\n' | sed 's/SAMPLES: //' | tr -d '"][")'); do

	date=$(echo ${sample} | cut -c 1-8)
	grep ${sample} sample_sheets/${date}_sample_sheet.tsv > demux/sample_sheet.tsv

	gunzip -c raw_data/${sample}.fastq.gz | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' | blastn -db demux/barcodes_primers -query - -task blastn-short -max_target_seqs 50 -outfmt "6 saccver qaccver slen qlen pident length qstart qend evalue" -gapopen 5 -gapextend 2 -reward 3 -penalty -4 -evalue 1 -num_threads 1 -perc_identity 70 > demux/blast_output.txt

	swibrid demultiplex -i raw_data/${sample}.fastq.gz -b demux/blast_output.txt --collapse -o input/ -s demux/sample_sheet.tsv -c 70
done

