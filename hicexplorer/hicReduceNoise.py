import argparse

from builtins import range

import numpy as np
from scipy.sparse import coo_matrix

from multiprocessing import Process, Queue
from multiprocessing.sharedctypes import Array, RawArray

from hicexplorer import HiCMatrix as hm

import matplotlib.pyplot as plt

from sklearn.preprocessing import normalize

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
            distribution[abs(i - j)] += hicMatrix.matrix[i, j]
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
        for i in range(len(distribution) - window_size):
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
            index_values_row = []
            index_values_column = []
            
            
            for j in range(i, i+window_size):
                for k in range():
                    # TODO get correct values
                    value = 9
                    values.append(value)
                
                mean_interaction_count_difference = abs(mean_absolute - distribution[i])
                change_values = [0] * len(values)
                if mean_absolute < distribution[i]:
                    # rich-getting-richer
                    for j in range(len(values)):
                        change_values[j] = ((1-(values[j] / distribution[i] )) ** power) * mean_interaction_count_difference
                    
                else:
                    # poor-getting-poorer
                    for j in range(len(values)):
                        change_values[j] = (((values[j] / distribution[i] )) ** power) * mean_interaction_count_difference
                # normalize to range [0, mean_interaction_count_difference]
                min_value = min(change_values)
                max_value = max(change_values)
                for j in range(len(change_values)):
                    change_values[j] = mean_interaction_count_difference * ((change_values[j] - min_value) / (max_value - min_value))
                
                if mean_absolute >= distribution[i]:
                    change_values = change_values * -1
                for j in range(len(values)):
                    values[j] += int(np.around(change_values))
                # TODO: write data back to matrix



            # print("mean_absolute", mean_absolute)
            # print("variance_normalized", variance_normalized)
            # print("values:",  distribution[i:i + window_size])
            # print("normalized:", normalized_values)


if __name__ == '__main__':
    reduceNoise = ReduceNoise()
    reduceNoise.calculateDistributionMatrix()
