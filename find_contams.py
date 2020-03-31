#!/usr/bin/env python
import pandas as pd
import numpy as np
import argparse
import subprocess
import shlex

def get_tables(profileDB, contigsDB, sample_name):
    cov_means_file = f'{sample_name}_coverage_Q2Q3_splits.txt'
    cov_stds_file = f'{sample_name}_std_coverage_splits.txt'
    detect_file = f'{sample_name}_detection_splits.txt'
    taxonomy_file = f'{sample_name}_taxonomy_splits.txt'
    mcov_cmd = f'anvi-export-table --table mean_coverage_Q2Q3_splits {profileDB} -o {cov_means_file}'
    scov_cmd = f'anvi-export-table --table std_coverage_splits {profileDB} -o {cov_stds_file}'
    dcov_cmd = f'anvi-export-table --table detection_splits {profileDB} -o {detect_file}'
    tax_cmd = f'anvi-export-splits-taxonomy -c {contigsDB} -o {taxonomy_file}'
    all_cmd = [mcov_cmd, scov_cmd, dcov_cmd, tax_cmd]

    for i in all_cmd:
        subprocess.run(shlex.split(i))

    return(cov_means_file, cov_stds_file, detect_file, taxonomy_file)

def import_table(anvio_file, sample_name):
    df = pd.read_csv(anvio_file, sep='\t')
    df = df.rename(columns = {'contig': 'split', '__parent__' : 'contig', sample_name : 's'})
    df = df.drop(columns = 'contig')
    df = df.set_index('split')

    return(df)

def score_contigs_cov(cov_means, cov_stds):
    '''Check whether the contig mean coverage is greater in sample than control and
    whether contig normalized standard deviation of coverage is less in sample than in control.
    If both are true, then the contig is not a contaminant.
    If one is false, contig will be in the returned dataframe for further analysis.'''

    # calculate relative coverage
    cov_means_norm = cov_means.div(cov_means.sum(axis=0), axis=1)
    # compare each column to s and make a matrix that is boolean control < sample (True is noncontam contig)
    cov_means_norm_lt = cov_means_norm.lt(cov_means_norm.s, axis=0)

    # get coefficient of variance (mean-normalized standard deviation of coverage for each contig)
    cov_stds_norm = cov_stds / cov_means
    cov_stds_norm = cov_stds_norm.replace([np.nan, np.inf], 0)
    # compare each column to s and make a matrix that is boolean control > sample (True is noncontam contig)
    cov_stds_norm_gt = cov_stds_norm.gt(cov_stds_norm.s, axis=0)

    ## combinations of cov_means_norm_lt and cov_stds_norm_gt
    # convert boolean to number and then add, so that:
    # True & True = 2 (noncontam)
    # True & False, False & True = 1 (unknown)
    # False & False = 0 (contam)

    iscontam = cov_means_norm_lt*1 + cov_stds_norm_gt*1

    ## special cases of False&True and True&False (both have value of 1)

    ## for any contig where (cov_means_norm.s == 0), change iscontam to 0 for all controls
    s_no_cov = (cov_means_norm.s == 0) * 1
    iscontam = iscontam.subtract(s_no_cov, axis=0) # subtract 1 in cases where sample has no cov.

    # for any contig x control where (cov_means_norm.c == 0), change iscontam to 2 for that control
    c_no_cov = (cov_means_norm==0 * 1)
    iscontam = iscontam.add(c_no_cov) # add 1 in cases where control has no cov.

    iscontam = iscontam.drop(columns = 's')

    # Remove all contigs where every control comparison yielded a 2 (noncontam)
    cov_contams = iscontam[~((iscontam == 2).all(axis = 'columns'))]

    return(cov_contams)

def score_contigs_detection(detection, cov_contams):
    '''Check detection on suspected contaminant contigs in all controls:
    if detection in control is less than detection in sample - 0.05 then remove from contaminant list
    because coverage is probably due to highly conserved sequence or mismapping, not true presence in sample'''

    # subset detection consider to just possible contaminant contigs
    detect_posscontam = detection[detection.index.isin(cov_contams.index)]

    # make masking df where:
    # True means detection in control is not within 0.05 of the detection in sample, so it is likely to be real
    # False means detection in control is similar to detection in sample, so it is likely to be a contaminant
    detect_posscontam_lt = detect_posscontam.lt(detect_posscontam.s - 0.05, axis=0)
    detect_posscontam_lt = detect_posscontam_lt.drop(columns = 's')

    # Use detect_posscontam_lt as a mask
    # remove any contigs that are 'True' for all controls
    # leaving scores for the remaining contigs
    final_contams = cov_contams[~detect_posscontam_lt]
    final_contams = final_contams[~((final_contams.isna()).all(axis = 'columns'))]

    return(final_contams)

def taxonomy_remove(tax_names, taxonomy_file):
    tax = pd.read_csv(taxonomy_file, sep = '\t', names=['split', 'genus'], index_col=0)
    tax_contam_splits = pd.DataFrame(tax[tax.genus.isin(tax_names)].index)
    return(tax_contam_splits)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--sample_name', type=str,
                        help='sample name used in anvi-gen-contigs-database')
    parser.add_argument('-p', '--profileDB', type=str,
                        help='anvio merged PROFILE.db')
    parser.add_argument('-c', '--contigsDB', type=str,
                        help='anvio CONTIGS.db')
    parser.add_argument('-o', '--out_prefix', type=str,
                        default='contam',
                        help='file name for collection of contaminant contigs, default: contam')
    parser.add_argument('-g', '--tax_names', type=str,
                        help='comma-separated list of known contaminant genera,\
                        e.g. Lysobacter,Asaccharospora')

    args = parser.parse_args()
    sample_name = args.sample_name
    out_prefix = '_'.join([sample_name, args.out_prefix])
    profileDB = args.profileDB
    contigsDB = args.contigsDB

    cov_means_file, cov_stds_file, detect_file, taxonomy_file = get_tables(profileDB, contigsDB, sample_name)

    cov_means = import_table(cov_means_file, sample_name)
    cov_stds = import_table(cov_stds_file, sample_name)
    detection = import_table(detect_file, sample_name)
    cov_contams = score_contigs_cov(cov_means, cov_stds)
    detect_contams = score_contigs_detection(detection, cov_contams)

    detect_contams = detect_contams.replace({0: 'contam', 1: 'unknown', 2: 'noncontam'})
    detect_contams = detect_contams.add_prefix('score_')

    # output table of scores to be imported as items for layers in anvio_file
    #(e.g. anvi-import-misc-data -p PROFILE.db -t items scores.txt)
    detect_contams.to_csv(f'{out_prefix}_scores.txt', index=True, sep='\t')
    collection = pd.DataFrame(detect_contams.index)

    if args.tax_names:
        tax_names = args.tax_names.split(',')
        tax_contam_splits = taxonomy_remove(tax_names, taxonomy_file)
        #combine contam_splits and collection, drop duplicate rows
        collection = pd.concat([collection, tax_contam_splits]).drop_duplicates('split')

    # output table of contig names to be used as a collection and split the DB
    #(with anvi-import-collection and anvi-split)
    collection['bin'] = 'contam'
    collection.to_csv(f'{out_prefix}_collection.txt', header=False, index=False, sep='\t')

if __name__ == '__main__':
    main()
