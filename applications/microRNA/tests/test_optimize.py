__author__ = 'asish'


import sys
import numpy as np
from Satyanveshi.optimize import init_parameters,\
    fold_parameters, flatten_parameters


def test_parameter_folding_unfolding():
    params = init_parameters(M=100, R=20, N=100, KM=10, KR=5)
    flat_params = flatten_parameters(params)
    unflattened_params = fold_parameters(flat_params)

    print "Flattened Params:", flat_params.shape
    print "Same number of params:", params.keys() == unflattened_params.keys()
    for k in params.keys():
        print "Param %s match? %s" % (k, np.all(params[k] == unflattened_params[k]))


if __name__ == '__main__':
    sys.exit(test_parameter_folding_unfolding())