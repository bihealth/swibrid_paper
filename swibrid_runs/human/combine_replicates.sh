#!/bin/bash -e

donors=$(cut -f 1 ../../sodar/2023_FB/a_2023_FB.txt | grep -v Sample | sort -u )

for donor in ${donors}; do
	samples=$(cut -f 1,17 ../../sodar/2023_FB/a_2023_FB.txt | grep -v Extract | grep ${donor} | cut -f 2 | sort -u | tr '\n' ',' | sed 's/,$//')
	[[ ! -s input/pooled_${donor}_info.csv ]] && echo ${donor} ${samples} && swibrid combine_replicates -s ${samples} -c pooled_${donor} --nmax 50000
done
