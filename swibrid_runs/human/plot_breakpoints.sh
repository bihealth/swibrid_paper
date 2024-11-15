#!/bin/bash -e

#sample=$1
mkdir -p paper_figures

#for sample in 20210422_plasmid_downsampled 20210618_FakePoly_downsampled 20211122_85_n500_1_downsampled 20220329_86_100000_1_downsampled 20231020_HD_K_downsampled 20231020_IRF4_downsampled 20231026_HD_K_downsampled 20231026_IRF4_downsampled 20231211_rFP_01_downsampled 20231211_rFP_04_downsampled 20231211_rFP_15_downsampled 20231211_rFP_22_downsampled 20231211_rFP_29_downsampled; do 
for sample in 20231211_rFP_01_downsampled 20231211_rFP_29_downsampled 20211122_85_n500_1_downsampled; do
	swibrid get_breakpoint_stats -g pipeline/${sample}/${sample}_gaps.npz -c pipeline/${sample}/${sample}_clustering.csv -a pipeline/${sample}/${sample}_cluster_analysis.csv -b 50 --max_gap 75 --switch_coords chr14:105583000-105872000:- --switch_annotation hg38_switch_regions.bed -o paper_figures/${sample}_breakpoints.csv -p paper_figures/${sample}_breakpoints.pdf --weights cluster
done
