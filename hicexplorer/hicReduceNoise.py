import argparse

from builtins import range

import numpy as np
from scipy.sparse import coo_matrix

from multiprocessing import Process, Queue
from multiprocessing.sharedctypes import Array, RawArray

from hicexplorer import HiCMatrix as hm

import matplotlib.pyplot as plt

from sklearn.preprocessing import normalize

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
    return parser


class ReduceNoise():

    def __init__(self):
        pass

    def calculateDistributionMatrix(self, args=None):
        args = parse_arguments().parse_args(args)
        hicMatrix = hm.hiCMatrix(args.matrix)
        print(hicMatrix.matrix.shape[0], hicMatrix.matrix.shape[1])
        distribution = np.zeros(hicMatrix.matrix.shape[0])
        instances, features = hicMatrix.matrix.nonzero()
        for i, j in zip(instances.tolist(), features.tolist()):
            if i - j < 0:
                continue
            distribution[i - j] += hicMatrix.matrix[i, j]
        plt.yscale('log')
        plt.xscale('log')
        plt.xlabel("Genomic distances")
        plt.ylabel("Number of interactions")
        plt.plot(list(range(len(distribution))), distribution)
        plt.savefig("interactions per genomic distance")
        # sliding window over genomic distances
        window_size = 5
        threshold_variance = 0.0008
        threshold_abs_mean = 80
        power = -3
        for i in xrange(0, len(distribution) - window_size, window_size):
            # compute mean
            # compute variance
            # if variance is too big, skip
            normalized_values = normalize(distribution[i:i + window_size].reshape(1, -1))
            mean_absolute = np.floor(np.mean(distribution[i:i + window_size]))
            variance_normalized = np.var(normalized_values)
            if variance_normalized > threshold_variance and mean_absolute > threshold_abs_mean:
                continue

            # set all genomic distances within window size to absolute mean value
            values = []
            # index_values_row = []
            # index_values_column = []
            print("sliding window iteration: ", i)
            for j in range(i, i + window_size):
                for k, l in zip(range(0, len(distribution)), range(j, len(distribution))):
                    # TODO get correct values
                    if hicMatrix.matrix[k, l] != 0:
                        values.append(hicMatrix.matrix[k, l])

                mean_interaction_count_difference = abs(mean_absolute - distribution[i])
                change_values = [0] * len(values)
                if mean_absolute < distribution[i]:
                    # rich-getting-richer
                    for k in range(len(values)):
                        change_values[k] = ((1 - (values[k] / distribution[i])) ** power) * mean_interaction_count_difference
                else:
                    # poor-getting-poorer
                    for k in range(len(values)):
                        change_values[k] = (((values[k] / distribution[i])) ** power) * mean_interaction_count_difference
                # normalize to range [0, mean_interaction_count_difference]
                min_value = min(change_values)
                max_value = max(change_values)
                for k in range(len(change_values)):
                    change_values[k] = mean_interaction_count_difference * ((change_values[k] - min_value) / (max_value - min_value))
                    if np.isnan(change_values[k]):
                        change_values[k] = 0
                if mean_absolute < distribution[i]:
                    for k in range(len(values)):
                        values[k] += int(change_values[k])
                else:
                    for k in range(len(values)):
                        values[k] -= int(change_values[k])
                        if values[k] < 0:
                            values[k] = 0
                # TODO: write data back to matrix

                values.reverse()
                for k, l in zip(range(0, len(distribution)), range(j, len(distribution))):
                    # TODO get correct values
                    if hicMatrix.matrix[k, l] != 0:
                        hicMatrix.matrix[k, l] = values.pop()
                values = []
                change_values = []
        distribution = np.zeros(hicMatrix.matrix.shape[0])
        instances, features = hicMatrix.matrix.nonzero()
        for i, j in zip(instances.tolist(), features.tolist()):
            if i - j < 0:
                continue
            distribution[abs(i - j)] += hicMatrix.matrix[i, j]
        plt.yscale('log')
        plt.xscale('log')
        plt.xlabel("Genomic distances")
        plt.ylabel("Number of interactions")
        plt.plot(list(range(len(distribution))), distribution)
        plt.savefig("interactions per genomic distance_corrected")

        hicMatrix.save("corrected.h5")

if __name__ == '__main__':
    reduceNoise = ReduceNoise()
    reduceNoise.calculateDistributionMatrix()
