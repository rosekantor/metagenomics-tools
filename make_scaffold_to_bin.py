#!/usr/bin/env python3

### make scaffold-to-bin lookup file

from Bio import SeqIO
import argparse
from glob import glob
import pandas as pd
import re
import os

def get_headers(fasta):
    headers_list = []
    for record in SeqIO.parse(fasta, "fasta"):
        headers_list.append(record.id)
    return headers_list

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--indir', type=str,
                        help='input directory containing fasta files of bins')
    parser.add_argument('-o', '--output_tsv', type=str,
                        help='output file')
    parser.add_argument('-x', '--suffix', type=str,
                        help='suffix of fasta', default='fasta')


    args = parser.parse_args()
    indir = args.indir
    output_tsv = args.output_tsv
    suffix = args.suffix

    out_df = []
    for fasta in glob(f'{indir}/*.{suffix}'):
        headers_list = get_headers(fasta)
        df = pd.DataFrame(headers_list, columns=['contig'])
        # get the bin name from the fasta file name
        # and make it pretty for anvio (rsplit to cut off suffix)
        bin = re.sub('[\.-]', '_', os.path.basename(fasta).rsplit('.', 1)[0])
        # add it to the dataframe
        df['bin'] = bin
        out_df.append(df)

    out_df = pd.concat(out_df)
    out_df.to_csv(output_tsv, sep='\t', index=False)

if __name__ == '__main__':
    main()
