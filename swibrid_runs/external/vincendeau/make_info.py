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
with gzip.open(args.input,'rt') as inf:
    for read in SeqIO.parse(inf,'fastq'):
        info[read.id] = {
            "length": len(read),
            "barcodes": "BC01@1-30:+",
            "primers": "primer_fw@0-50:+;primer_rv@{0}-{1}:+".format(
                len(read) - 50, len(read)
            ),
            "internal": "",
        }

pd.DataFrame.from_dict(info, orient="index").to_csv(args.output)

