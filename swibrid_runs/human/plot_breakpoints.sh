#!/bin/bash -e

for sample in 20231211_rFP_01_downsampled 20250113_SKA1_downsampled 20211122_85_n500_1_downsampled; do
	swibrid get_breakpoint_stats -g pipeline/${sample}/${sample}_gaps.npz -c pipeline/${sample}/${sample}_clustering.csv -a pipeline/${sample}/${sample}_cluster_analysis.csv -b 50 --max_gap 75 --switch_coords chr14:105583000-105872000:- --switch_annotation hg38_switch_regions.bed -o ${sample}_breakpoints.csv -p ${sample}_breakpoints.pdf --weights cluster
done
