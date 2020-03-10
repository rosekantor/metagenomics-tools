#!/usr/bin/env python
import pandas as pd
import numpy as np
import argparse

def import_table(anvio_file, sample_name):
    # import coverage means table (we requested it in --contigs-mode
    # so even though we are working with splits,
    # all splits will have the same values across a contig)
    df = pd.read_csv(anvio_file, sep='\t')
    df = df.rename(columns = {'contig': 'split', '__parent__' : 'contig', sample_name : 's'})
    df = df.drop(columns = 'contig')
    df = df.set_index('split')

    return(df)

def score_contigs(cov_means, cov_stds):
    # calculate relative coverage
    cov_means_norm = cov_means.div(cov_means.sum(axis=0), axis=1)
    # compare each column to s and make a matrix that is binary control < sample --> noncontam contig
    cov_means_norm_lt = cov_means_norm.lt(cov_means_norm.s, axis=0)

    # get coefficient of variance (mean-normalized standard deviation of coverage for each contig)
    cov_stds_norm = cov_stds / cov_means
    cov_stds_norm = cov_stds_norm.replace([np.nan, np.inf], 0)
    # compare each column to s and make a matrix that is binary control > sample --> noncontam contig
    cov_stds_norm_gt = cov_stds_norm.gt(cov_stds_norm.s, axis=0)

    ## Rules
    # identify contig x control combinations where:
    # iscontam = 'contam': (cov_means_norm_lt==False and cov_stds_norm_gt==False) score=0 or cov_means_norm.s==0
    # iscontam = 'real': (cov_means_norm_lt==True and cov_stds_norm_gt==True) score=2 or cov_means_norm.c==0
    # iscontam = 'unknown': ((cov_means_norm_lt==False and cov_stds_norm_gt==True) or
    #                        (cov_means_norm_lt==True and cov_stds_norm_gt==False)) score=1

    ## First look at just combinations of cov_means_norm_lt and cov_stds_norm_gt
    # convert boolean to number and then add, so that True&True = 2, True&False = 1, False&False = 0
    iscontam = cov_means_norm_lt*1 + cov_stds_norm_gt*1

    ## resolve cases of False&True and True&False (both have value of 1)

    ## for any contig where (cov_means_norm.s == 0), change iscontam to 0 for all controls
    s_no_cov = (cov_means_norm.s == 0) * 1
    iscontam = iscontam.subtract(s_no_cov, axis=0) # subtract 1 in cases where sample has no cov.

    # for any contig x control where (cov_means_norm.c == 0), change iscontam to 2 for that control
    c_no_cov = (cov_means_norm==0 * 1)
    iscontam = iscontam.add(c_no_cov) # add 1 in cases where control has no cov.

    iscontam = iscontam.drop(columns = 's')
    iscontam = iscontam.add_prefix('score_')

    ## convert score [0,1,2] to ['contam', 'unknown', 'noncontam']
    iscontam = iscontam.replace({0: 'contam', 1: 'unknown', 2: 'noncontam'})

    return(iscontam)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--sample_name', type=str,
                        help='sample name used in anvi-gen-contigs-database')
    parser.add_argument('-m', '--cov_means_file', type=str,
                        help='file produced by anvi-export-table -t mean_coverage_Q2Q3_contigs PROFILE.db')
    parser.add_argument('-s', '--cov_stds_file', type=str,
                        help='file produced by anvi-export-table -t std_coverage_contigs PROFILE.db')
    parser.add_argument('-o', '--out_prefix', type=str,
                        default='contam',
                        help='file name for collection of contaminant contigs, default: contam_collection.txt')

    args = parser.parse_args()
    sample_name = args.sample_name
    cov_means_file = args.cov_means_file
    cov_stds_file = args.cov_stds_file
    out_prefix = args.out_prefix

    cov_means = import_table(cov_means_file, sample_name)
    cov_stds = import_table(cov_stds_file, sample_name)
    iscontam = score_contigs(cov_means, cov_stds)

    # collect only scaffolds that are possible contam by excluding any row that contains all 'noncontam'
    possible_contams = iscontam[~((iscontam == 'noncontam').all(axis = 'columns'))]
    # output table of scores to be imported as items (e.g. anvi-import-misc-data -p PROFILE.db -t items scores.txt)
    possible_contams.to_csv(f'{out_prefix}_scores.txt', index=True, sep='\t')
    # output table of contig names to be used as a collection and split (anvi-import-collection and anvi-split)
    collection = pd.DataFrame(possible_contams.index)
    collection['bin'] = 'contam'
    collection.to_csv(f'{out_prefix}_collection.txt', header=False, index=False, sep='\t')

if __name__ == '__main__':
    main()
