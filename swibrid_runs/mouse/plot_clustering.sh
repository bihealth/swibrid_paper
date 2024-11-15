#!/bin/bash -e

#sample=$1
mkdir -p paper_figures

for sample in 20240817_S.A3_downsampled 20220425_WT_3_downsampled; do
	swibrid plot_clustering --msa pipeline/${sample}/${sample}_msa.npz --figure paper_figures/${sample}.png --linkage pipeline/${sample}/${sample}_linkage.npz --info input/${sample}_info.csv  --clustering_results pipeline/${sample}/${sample}_clustering.csv --switch_annotation mm10_switch_regions.bed --switch_coords "chr12:113200000-113500000:-" --color_by cluster --clustering_stats pipeline/${sample}/${sample}_cluster_stats.csv --sample ${sample} --dpi 500 --fig_width 4 --fig_height 1.25 --cmax 0.5 --linkage_border .2 --cutoff 0.01 --cmax 0.1 --no_x_ticks --omit_scale_bar --dpi 600
done

