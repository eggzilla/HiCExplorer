import argparse

from builtins import range
from six import iteritems

import numpy as np
from scipy.sparse import coo_matrix

from multiprocessing import Process, Queue
from multiprocessing.sharedctypes import Array, RawArray

from hicexplorer import HiCMatrix as hm
import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt

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
                        default=100,
                        type=int
                        )
    parser.add_argument('--threshold_variance',
                        help='',
                        required=False,
                        default=0.1,
                        type=float
                        )
    parser.add_argument('--removeLowInteractionCount',
                        help='',
                        required=False,
                        default=1,
                        type=int
                        )
    parser.add_argument('--power',
                        help='',
                        required=False,
                        default=-2,
                        type=float
                        )
    parser.add_argument('--output',
                        help='',
                        required=False,
                        default='output',
                        type=str
                        )
    return parser


class ReduceNoise():

    def __init__(self):
        pass

    def zero_to_nan(self, values):
        """Replace every 0 with 'nan' and return a copy."""
        return [float('nan') if x == 0 else x for x in values]

    def calculateDistributionMatrix(self, args=None):
        # print("Calcuate genomic distance for plot...")
        args = parse_arguments().parse_args(args)
        # hicMatrix = hm.hiCMatrix(args.matrix)

        # print(hicMatrix.matrix.shape[0], hicMatrix.matrix.shape[1])
        # distribution = np.zeros(hicMatrix.matrix.shape[0])
        # instances, features = hicMatrix.matrix.nonzero()
        # for i, j in zip(instances.tolist(), features.tolist()):
        #     if i - j < 0:
        #         continue
        #     distribution[i - j] += hicMatrix.matrix[i, j]
        # distribution = self.zero_to_nan(distribution)
        # plt.yscale('log')
        # plt.xscale('log')
        # plt.xlabel("Genomic distances")
        # plt.ylabel("Number of interactions")
        # plt.plot(list(range(len(distribution))), distribution)
        # plt.savefig("interactions per genomic distance")
        # print("Calcuate genomic distance for plot...Done!")

        # # print(distribution)
        # sliding window over genomic distances
        window_size = args.window_size
        threshold_variance = args.threshold_variance
        power = args.power
        threads = args.threads
        # iterations = args.iterations
        # calculate changes in c++
        print("Start correction in C++...")
        # instances_new, features_new, data_new = _c_noise_reduction.c_powerLawNoiseReduction(hicMatrix.matrix.nonzero()[0].tolist(), hicMatrix.matrix.nonzero()[1].tolist(), hicMatrix.matrix.data.tolist(),
        #                                                                                     window_size, threshold_variance, hicMatrix.matrix.count_nonzero(),
        #                                                                                     hicMatrix.matrix.shape[0], power, threads,
        #                                                                                     args.removeLowInteractionCount)
        _c_noise_reduction.c_powerLawNoiseReduction_h5(args.matrix,
                                                       window_size, threshold_variance,
                                                       power, threads,
                                                       args.removeLowInteractionCount,
                                                       args.output)

        print("Start correction in C++...Done!")

        # hicMatrix.matrix = csr_matrix((data_new, (instances_new, features_new)), shape=(hicMatrix.matrix.shape[0], hicMatrix.matrix.shape[0]))
        # print(len(hicMatrix.matrix.nonzero()[0]), len(hicMatrix.matrix.nonzero()[1]))
        # distribution = np.zeros(hicMatrix.matrix.shape[0])
        # print("Calcuate genomic distance for plot...")

        # for i, j in zip(instances_new, features_new):
        #     if i - j < 0:
        #         continue
        #     distribution[abs(i - j)] += hicMatrix.matrix[i, j]
        # distribution = self.zero_to_nan(distribution)

        # plt.yscale('log')
        # plt.xscale('log')
        # plt.xlabel("Genomic distances")
        # plt.ylabel("Number of interactions")
        # plt.plot(list(range(len(distribution))), distribution)
        # plt.savefig("interactions per genomic distance_corrected")
        # print("Calcuate genomic distance for plot...Done!")

        # hicMatrix.save(args.output)
        # print(distribution)


if __name__ == '__main__':
    reduceNoise = ReduceNoise()
    reduceNoise.calculateDistributionMatrix()
