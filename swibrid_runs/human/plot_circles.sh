#!/bin/bash -e

python plot_circles.py --samples "20231211_rFP_01_downsampled,20231211_rFP_29_downsampled,20211122_85_n500_1_downsampled" --ncols 3 --figure paper_figures/20240320_bubble_plots_downsampled_cluster.pdf --color cluster
#python plot_circles.py --samples "20231211_rFP_01_downsampled,20231211_rFP_29_downsampled,20211122_85_n500_1_downsampled" --ncols 3 --figure paper_figures/20240320_bubble_plots_downsampled_isotype.pdf --color isotype

