#!/bin/bash -e

# use this to cluster two random samples of donor ${donor} with ${ncells} cells
#python meta_clustering.py -s ${seed} -d ${donor} -n ${ncells} -o meta_clustering/${donor}_${seed}_${ncells}_meta_clustering.csv
# use this to cluster all replicates of donor ${donro} with ${ncells} cells
#python meta_clustering.py -d ${donor} -n ${ncells} -o meta_clustering/${donor}_${ncells}_meta_clustering.csv
# use this to cluster minION and pacBio samples of HD19005
#python meta_clustering.py --samples 20230313_HD19005,20230504_HD19005 -o meta_clustering/HD19005_meta_clustering.csv
# use this to cluster replicates with different PCR cycle numbers
#python meta_clustering.py --samples 20250703_${donor}x15,20250703_${donor}x20,20250703_${donor}x25 -o meta_clustering/${donor}_meta_clustering.csv


