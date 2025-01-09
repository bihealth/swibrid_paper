import os
import sys
from Bio import SeqIO
import pandas as pd
import gzip
from argparse import ArgumentParser

parser=ArgumentParser()
parser.add_argument('-i','--input',dest='input',help="""Input fastq.gz file""")
parser.add_argument('-o','--output',dest='output',help="""output info.csv file""")

args=parser.parse_args()

info={}
read_id = ""
with gzip.open(args.input,'rt') as inf:
    reads=SeqIO.parse(inf,'fastq')
    while True:
        try:
            read = next(reads)
        except StopIteration:
            break
        if read.id == read_id:
            info[read.id]['length'] += len(read)
            info[read.id]['num_mates'] += 1
        else:
            info[read.id] = {
                "length": len(read),
                "num_mates": 1,
                "barcodes": "BC01@1-30:+",
                "primers": "primer_fw@0-20:+;primer_rv@{0}-{1}:+".format(
                    max(0,len(read) - 20), len(read)
                ),
                "internal": "",
            }
        read_id = read.id

pd.DataFrame.from_dict(info, orient="index").to_csv(args.output)

