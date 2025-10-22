#!/bin/bash -e

for sample in 20231211_rFP_01_downsampled 20250113_SKA1_downsampled 20211122_85_n500_1_downsampled 20210422_plasmid_downsampled; do
	swibrid plot_clustering --msa pipeline/${sample}/${sample}_msa.npz --figure ${sample}_reads_by_cluster.png --linkage pipeline/${sample}/${sample}_linkage.npz --info input/${sample}_info.csv  --clustering_results pipeline/${sample}/${sample}_clustering.csv --switch_annotation hg38_switch_regions.bed --switch_coords "chr14:105583000-105872000:-" --color_by cluster --clustering_stats pipeline/${sample}/${sample}_cluster_stats.csv --sample ${sample} --dpi 500 --fig_width 4 --fig_height 1.25 --cmax 0.5 --linkage_border .2 --cutoff 0.01 --cmax 0.1 --no_x_ticks --omit_scale_bar --dpi 600
done

for sample in 20211122_85_n500_1_downsampled; do
	swibrid plot_clustering --msa pipeline/${sample}/${sample}_msa.npz --figure ${sample}_reads_by_coverage.png --linkage pipeline/${sample}/${sample}_linkage.npz --info input/${sample}_info.csv  --clustering_results pipeline/${sample}/${sample}_clustering.csv --switch_annotation hg38_switch_regions.bed --switch_coords "chr14:105583000-105872000:-" --color_by coverage --clustering_stats pipeline/${sample}/${sample}_cluster_stats.csv --sample ${sample} --dpi 500 --fig_width 4 --fig_height 1.25 --cmax 0.5 --linkage_border .2 --cutoff 0.01 --cmax 0.1 --no_x_ticks --omit_scale_bar --dpi 600 --variants_matrix pipeline/${sample}/${sample}_variants.npz 
done
