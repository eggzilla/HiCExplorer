from __future__ import division
import argparse

from scipy.sparse import csr_matrix, lil_matrix
import numpy as np

from hicexplorer import HiCMatrix as hm
from hicexplorer._version import __version__
from hicexplorer.utilities import exp_obs_matrix_lieberman
from hicexplorer.utilities import convertNansToZeros, convertInfsToZeros


import logging
log = logging.getLogger(__name__)



def parse_arguments(args=None):
    """
    get command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Converts the (interaction) matrix to a observed/expected matrix or a pearson_correlated.')

    # define the arguments
    parser.add_argument('--matrix', '-m',
                        help='input file. The computation is done per chromosome.',
                        required=True)

    parser.add_argument('--outFileName', '-o',
                        help='File name to save the exported matrix.',
                        required=True)
    parser.add_argument('--threads', '-t',
                        help='Number of threads for pearson correlation.',
                        required=False,
                        default=4,
                        type=int)

    parser.add_argument('--method', '-me',
                        help='Transformation method to use. If the option all is used, all three matrices in '
                        'consecutively way (input -> obs_exp -> pearson -> covariance) are created.',
                        choices=['obs_exp', 'pearson', 'covariance', 'all'],
                        default='obs_exp')

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))
    parser.add_argument('--chromosomes',
                        help='List of chromosomes to be included in the '
                        'correlation.',
                        default=None,
                        nargs='+')
    return parser


def __obs_exp(pSubmatrix, pLengthChromosome, pChromosomeCount):

    exp_obs_matrix_ = exp_obs_matrix_lieberman(pSubmatrix, pLengthChromosome, pChromosomeCount)
    exp_obs_matrix_ = convertNansToZeros(csr_matrix(exp_obs_matrix_))
    exp_obs_matrix_ = convertInfsToZeros(csr_matrix(exp_obs_matrix_)).todense()
    return exp_obs_matrix_


def __pearson(pSubmatrix):
    pearson_correlation_matrix = np.corrcoef(pSubmatrix)
    pearson_correlation_matrix = convertNansToZeros(csr_matrix(pearson_correlation_matrix))
    pearson_correlation_matrix = convertInfsToZeros(csr_matrix(pearson_correlation_matrix)).todense()
    return pearson_correlation_matrix


def main(args=None):
    
    args = parse_arguments().parse_args(args)

    hic_ma = hm.hiCMatrix(matrixFile=args.matrix)

    if args.chromosomes:
        hic_ma.keepOnlyTheseChr(args.chromosomes)

    length_chromosome = 0
    chromosome_count = len(hic_ma.getChrNames())
    for chrname in hic_ma.getChrNames():
        chr_range = hic_ma.getChrBinRange(chrname)
        length_chromosome += chr_range[1] - chr_range[0]
    trasf_matrix = lil_matrix(hic_ma.matrix.shape)

    if args.method == 'obs_exp':
        for chrname in hic_ma.getChrNames():
            chr_range = hic_ma.getChrBinRange(chrname)
            submatrix = hic_ma.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]]
            submatrix.astype(float)
            trasf_matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = lil_matrix(__obs_exp(submatrix, length_chromosome, chromosome_count))

    elif args.method == 'pearson':
        for chrname in hic_ma.getChrNames():
            chr_range = hic_ma.getChrBinRange(chrname)
            submatrix = hic_ma.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]]
            submatrix.astype(float)
            trasf_matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = lil_matrix(__pearson(submatrix.todense()))

    elif args.method == 'covariance':
        for chrname in hic_ma.getChrNames():
            chr_range = hic_ma.getChrBinRange(chrname)
            submatrix = hic_ma.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]]
            submatrix.astype(float)
            corrmatrix = np.cov(submatrix.todense())
            trasf_matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = lil_matrix(corrmatrix)

    elif args.method == 'all':
        trasf_matrix_obs_exp = lil_matrix(hic_ma.matrix.shape)
        trasf_matrix_pearson = lil_matrix(hic_ma.matrix.shape)
        trasf_matrix_corr = lil_matrix(hic_ma.matrix.shape)

        for chrname in hic_ma.getChrNames():
            chr_range = hic_ma.getChrBinRange(chrname)
            submatrix = hic_ma.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]]
            submatrix.astype(float)
            submatrix = __obs_exp(submatrix, length_chromosome, chromosome_count)
            trasf_matrix_obs_exp[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = lil_matrix(submatrix)
            submatrix = __pearson(submatrix)
            trasf_matrix_pearson[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = lil_matrix(submatrix)
            corrmatrix = np.cov(submatrix)
            trasf_matrix_corr[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = lil_matrix(submatrix)

        hic_ma.setMatrix(trasf_matrix_obs_exp.tocsr(), cut_intervals=hic_ma.cut_intervals)
        hic_ma.save("obs_exp_" + args.outFileName, pSymmetric=False)

        hic_ma.setMatrix(trasf_matrix_pearson.tocsr(), cut_intervals=hic_ma.cut_intervals)
        hic_ma.save("pearson_" + args.outFileName, pSymmetric=False)

        hic_ma.setMatrix(trasf_matrix_corr.tocsr(), cut_intervals=hic_ma.cut_intervals)
        hic_ma.save("covariance_" + args.outFileName, pSymmetric=False)

    if not args.method == 'all':
        hic_ma.setMatrix(trasf_matrix.tocsr(), cut_intervals=hic_ma.cut_intervals)
        hic_ma.save(args.outFileName, pSymmetric=False)
