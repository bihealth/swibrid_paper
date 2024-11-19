#!/bin/bash 

cat index/barcodes_primers.fa | makeblastdb -in - -title "barcodes_plasmid_primers" -out demux/barcodes_plasmid_primers -input_type fasta -dbtype nucl

gunzip -c raw_data/${dataset}.fastq.gz | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' | blastn -db demux/barcodes_plasmid_primers -query - -task blastn-short -max_target_seqs 50 -outfmt "6 saccver qaccver slen qlen pident length qstart qend evalue" -gapopen 5 -gapextend 2 -reward 3 -penalty -4 -evalue 1 -num_threads 1 -perc_identity 70 > demux/blast_output_${dataset}.txt

mkdir -p demux/${dataset}
swibrid demultiplex -i raw_data/${dataset}.fastq.gz -b demux/blast_output_${dataset}.txt -f demux/${dataset}_demultiplexing.png --collapse -o demux/${dataset} -r demux/${dataset}_demultiplexing_stats.csv --split-reads -s  sample_sheets/${dataset}_sample_sheet.tsv -c 70

mkdir -p input
ln -rfs demux/${dataset}/* input
