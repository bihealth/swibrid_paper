import numpy as np
import pandas as pd
import fastcluster
import glob
import scipy.sparse
import os
from swibrid.utils import remove_gaps
from argparse import ArgumentParser

parser=ArgumentParser()
parser.add_argument('-d','--donor')
parser.add_argument('-r','--replicate', default=3)
parser.add_argument('-s','--seed', type=int)
parser.add_argument('-n','--ncells',type=int)
parser.add_argument('-o','--out')
parser.add_argument('--samples')
args=parser.parse_args()

study=pd.read_csv('../../sodar/2020_CellSwitch/s_2020_CellSwitch.txt', sep='\t',header=0)
assay=pd.read_csv('../../sodar/2020_CellSwitch/a_2020_CellSwitch.txt', sep='\t',header=0)
QC=assay.set_index('Sample Name')['Parameter Value[QC]']

if args.samples:
    samples=args.samples.split(',')
else:
    ncells=args.ncells
    if args.donor and args.donor!='mixed':
        donor=args.donor
        if args.seed is not None:
            seed=args.seed
            np.random.seed(seed)
            samples=study['Sample Name'][(study['Source Name'].astype(str)==donor).values & (QC[study['Sample Name']] == 'PASS').values & (study['Parameter Value[Cell Number]']==ncells).values].sample(2).values
        else:
            samples=study['Sample Name'][(study['Source Name'].astype(str)==donor).values & (QC[study['Sample Name']] == 'PASS').values & (study['Parameter Value[Cell Number]']==ncells).values].values
    else:
        seed=args.seed
        np.random.seed(seed)
        samples=[study['Sample Name'][(study['Source Name'].astype(str)==donor).values & (QC[study['Sample Name']] == 'PASS').values & (study['Parameter Value[Cell Number]']==ncells).values].sample(1).values[0] for donor in np.random.choice(['21084','21085','21086'], size=2)]
    
nsamples=len(samples)
print(f"meta-clustering for samples " + ",".join(samples))

msa=dict((sample,scipy.sparse.load_npz(f'pipeline/{sample}/{sample}_msa.npz')) for sample in samples)
clusterings=dict((sample,pd.read_csv(f'pipeline/{sample}/{sample}_clustering.csv',header=0,index_col=0)) for sample in samples)
gaps=dict((sample,np.load(f'pipeline/{sample}/{sample}_gaps.npz')) for sample in samples)

msa_cleaned=dict((s,remove_gaps(msa[s], gaps=gaps[s], max_gap=75)) for s in samples)

mm={}
clusters={}
for s in samples:
    clusters[s], cinv, csize = np.unique(
        clusterings[s]["cluster"].dropna(), return_inverse=True, return_counts=True
    )
    mm[s] = scipy.sparse.csr_matrix(
        (1.0 / csize[cinv], (cinv, np.arange(len(cinv)))),
        shape=(len(clusters[s]), len(cinv)),
    )

tot_avg_msa = np.vstack([np.asarray(mm[s].dot(msa_cleaned[s]).todense()) for s in samples])
nreads=tot_avg_msa.shape[0]

tot_clusters=[f'{s}_{c}' for s in samples for c in clusters[s]]

print(f'clustering {nsamples} samples with {nreads} consensus reads')

Z = fastcluster.linkage(tot_avg_msa, metric='cosine', method='average')

cc = scipy.cluster.hierarchy.cut_tree(Z, height=.01)
meta_cluster=pd.DataFrame({'meta_cluster': cc[:,0]}, index=tot_clusters)

for s in samples:
    clusterings[s]['sample']=s

meta=pd.concat([clusterings[s] for s in samples],axis=0)[['isotype','orientation','cluster','filtered_cluster','sample']]
meta['cluster']=meta['sample']+'_'+meta['cluster'].astype(str)
meta['meta_cluster']=meta_cluster.loc[meta['cluster']]['meta_cluster'].values

meta.to_csv(args.out, index=True)
