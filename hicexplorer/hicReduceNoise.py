import argparse

from builtins import range

import numpy as np
from scipy.sparse import coo_matrix

from multiprocessing import Process, Queue
from multiprocessing.sharedctypes import Array, RawArray

from hicexplorer import HiCMatrix as hm

# import matplotlib.pyplot as plt

# from sklearn.preprocessing import normalize
from scipy.sparse import csr_matrix

import _c_noise_reduction


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=('Corrects a given hicMatrix based on z-scores and fdr.'))

    # define the arguments
    parser.add_argument('--matrix', '-m',
                        help='Path of the Hi-C matrix to plot',
                        required=True)
    parser.add_argument('--threads',
                        help='Number of threads. Using the python multiprocessing module.',
                        required=False,
                        default=4,
                        type=int
                        )
    parser.add_argument('--window_size',
                        help='',
                        required=False,
                        default=5,
                        type=int
                        )
    parser.add_argument('--threshold_variance',
                        help='',
                        required=False,
                        default=0.1,
                        type=float
                        )
    parser.add_argument('--threshold_abs_mean',
                        help='',
                        required=False,
                        default=100,
                        type=float
                        )
    parser.add_argument('--power',
                        help='',
                        required=False,
                        default=-3,
                        type=float
                        )
    parser.add_argument('--iterations',
                        help='',
                        required=False,
                        default=5,
                        type=int
                        )
    return parser


class ReduceNoise():

    def __init__(self):
        pass

    def calculateDistributionMatrix(self, args=None):
        args = parse_arguments().parse_args(args)
        hicMatrix = hm.hiCMatrix(args.matrix)
        # print(hicMatrix.matrix.shape[0], hicMatrix.matrix.shape[1])
        # distribution = np.zeros(hicMatrix.matrix.shape[0])
        instances, features = hicMatrix.matrix.nonzero()
        # for i, j in zip(instances.tolist(), features.tolist()):
        #     if i - j < 0:
        #         continue
        #     distribution[i - j] += hicMatrix.matrix[i, j]
        # plt.yscale('log')
        # plt.xscale('log')
        # plt.xlabel("Genomic distances")
        # plt.ylabel("Number of interactions")
        # plt.plot(list(range(len(distribution))), distribution)
        # plt.savefig("interactions per genomic distance")
        # sliding window over genomic distances
        window_size = args.window_size
        threshold_variance = args.threshold_variance
        threshold_abs_mean = args.threshold_abs_mean
        power = args.power
        threads = args.threads
        iterations = args.iterations
        # calculate changes in c++
        instances_new, features_new, data_new = _c_noise_reduction.c_powerLawNoiseReduction(instances.tolist(), features.tolist(), hicMatrix.matrix.data.tolist(),
                                                                                            window_size, threshold_variance, threshold_abs_mean,
                                                                                            len(instances.tolist()), hicMatrix.matrix.shape[0], power, threads,
                                                                                            iterations)

        hicMatrix.matrix = csr_matrix((data_new, (instances_new, features_new)), shape=(hicMatrix.matrix.shape[0], hicMatrix.matrix.shape[0]))
        distribution = np.zeros(hicMatrix.matrix.shape[0])
        
        # for i, j in zip(instances_new, features_new):
        #     if i - j < 0:
        #         continue
        #     distribution[abs(i - j)] += hicMatrix.matrix[i, j]
        # plt.yscale('log')
        # plt.xscale('log')
        # plt.xlabel("Genomic distances")
        # plt.ylabel("Number of interactions")
        # plt.plot(list(range(len(distribution))), distribution)
        # plt.savefig("interactions per genomic distance_corrected")

        hicMatrix.save("corrected.h5")


if __name__ == '__main__':
    reduceNoise = ReduceNoise()
    reduceNoise.calculateDistributionMatrix()
