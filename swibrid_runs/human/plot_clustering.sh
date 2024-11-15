#!/bin/bash -e

#sample=$1
mkdir -p paper_figures

#for sample in 20210422_plasmid_downsampled 20210618_FakePoly_downsampled 20211122_85_n500_1_downsampled 20220329_86_100000_1_downsampled 20231020_HD_K_downsampled 20231020_IRF4_downsampled 20231026_HD_K_downsampled 20231026_IRF4_downsampled 20231211_rFP_01_downsampled 20231211_rFP_04_downsampled 20231211_rFP_15_downsampled 20231211_rFP_22_downsampled 20231211_rFP_29_downsampled; do 
for sample in 20231211_rFP_01_downsampled 20231211_rFP_29_downsampled 20211122_85_n500_1_downsampled 20210422_plasmid_downsampled; do
	swibrid plot_clustering --msa pipeline/${sample}/${sample}_msa.npz --figure paper_figures/${sample}.png --linkage pipeline/${sample}/${sample}_linkage.npz --info input/${sample}_info.csv  --clustering_results pipeline/${sample}/${sample}_clustering.csv --switch_annotation hg38_switch_regions.bed --switch_coords "chr14:105583000-105872000:-" --color_by cluster --clustering_stats pipeline/${sample}/${sample}_cluster_stats.csv --sample ${sample} --dpi 500 --fig_width 4 --fig_height 1.25 --cmax 0.5 --linkage_border .2 --cutoff 0.01 --cmax 0.1 --no_x_ticks --omit_scale_bar --dpi 600
done

for sample in 20211122_85_n500_1_downsampled; do
	swibrid plot_clustering --msa pipeline/${sample}/${sample}_msa.npz --figure paper_figures/${sample}_cov.png --linkage pipeline/${sample}/${sample}_linkage.npz --info input/${sample}_info.csv  --clustering_results pipeline/${sample}/${sample}_clustering.csv --switch_annotation hg38_switch_regions.bed --switch_coords "chr14:105583000-105872000:-" --color_by coverage --clustering_stats pipeline/${sample}/${sample}_cluster_stats.csv --sample ${sample} --dpi 500 --fig_width 4 --fig_height 1.25 --cmax 0.5 --linkage_border .2 --cutoff 0.01 --cmax 0.1 --no_x_ticks --omit_scale_bar --dpi 600 --variants_matrix pipeline/${sample}/${sample}_variants.npz 
done
