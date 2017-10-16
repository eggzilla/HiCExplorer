from __future__ import division

import argparse
import sys
import os
import numpy as np
from builtins import range
from past.builtins import map

from intervaltree import IntervalTree, Interval

# from scipy.sparse import triu
# from scipy.stats import pearsonr, spearmanr

import hicexplorer.HiCMatrix as hm
from hicexplorer._version import __version__

# for plotting
from matplotlib import use as mplt_use
mplt_use('Agg')

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FixedLocator


def parse_arguments(args=None):

    # heatmap_parser = heatmap_options()

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Computes differentally analysis based on the comparison of the found tad domains.')

    # define the arguments
    parser.add_argument('--tadFiles', '-m',
                        help='Tad score files to compute the differntial expression of the chromatine.',
                        nargs=2,
                        type=str,
                        required=True,
                        metavar='BED file')

    parser.add_argument('--chromosomes',
                        help='List of chromosomes to be included in the '
                        'correlation.',
                        default=None,
                        nargs='+')
    parser.add_argument('--outFileName', '-o',
                        help='File name to save the resulting matrix. '
                        'The output is also a .h5 file. But don\'t add '
                        'the suffix',
                        required=True)
    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser

def main(args=None):
    """This differentally analysis is based on the comparison of the found tad domains."""
    args = parse_arguments().parse_args(args)

    # check for same bin size
    # accept only if match


    tad_domains = []
    print args.tadFiles
    for tadFile in args.tadFiles:
        print tadFile
        # tad_domains[tadFile] = []
        intervals = []
        with open(tadFile, 'r') as tadFile_:
            chrom_old = None
            for line in tadFile_.readlines():
                chrom, start, end  = line.split('\t')[:3]
                intervals.append((chrom, start, end))
        tad_domains.append(intervals)
        intervals = None
    
    hicMatrix = hm.hiCMatrix()
    interval_tree = [None] * len(tad_domains)
    i = 0
    while i < len(interval_tree):
        interval_tree[i], _ = hicMatrix.intervalListToIntervalTree(tad_domains[i])
        i += 1
    i = None

    result = {}
    # print()
    print("interval_tree", len(interval_tree))
    print("interval_tree", len(interval_tree[0]))
    print("interval_tree", len(interval_tree[1]))
    
    # compute intersection of interval trees per chromosome
    for chrom in interval_tree[0]:
        i = 1
        print("chrom", chrom)
        while i < len(interval_tree) and chrom in interval_tree[i]:
            print("chrom_inner")
            
            if chrom in result:
                print('chrom seen, create intersection: ', chrom)
                result[chrom] &= interval_tree[i][chrom]
            else:
                print('chrom new, create intersection: ', chrom)
                print('interval_tree[0][chrom]', interval_tree[0][chrom])
                print('interval_tree[i][chrom]', interval_tree[i][chrom])
                
                result[chrom] = interval_tree[0][chrom] & interval_tree[i][chrom] 
            i += 1
        if i != len(interval_tree):
            result[chrom] = None
        

        # for tree in interval_tree[1:]:
        #     if chrom in tree:
        #         # print('tree', tree)
        #         if chrom in result:
                    
        #             result[chrom] &= tree[chrom]
        #         else:
        #             result[chrom] = interval_tree[0][chrom] & tree[chrom] 
        #     else:
        #         result[chrom] = None
        #         break

    print result
    # return tads which occure in all files

    # with open(args.outFileName, 'w') as outFile:
    #     for chrom in result:
    #         if result[chrom] is not None:
    #             for interval in result[chrom]:
    #                 outFile.write("{}\t{}\t{}\t{}".format(chrom, interval[0], interval[1]))
