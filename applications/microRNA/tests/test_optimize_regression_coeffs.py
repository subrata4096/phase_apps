__author__ = 'asish'


import sys
import numpy as np
from Satyanveshi.optimize_regression_coeffs import compute_regression_coeffs
from Satyanveshi.simulation import generate_simulated_data
import matplotlib.pyplot as plt
import time


def hinton(matrix, max_weight=None, ax=None):
    """Draw Hinton diagram for visualizing a weight matrix."""
    ax = ax if ax is not None else plt.gca()

    if not max_weight:
        max_weight = 2**np.ceil(np.log(np.abs(matrix).max())/np.log(2))

    ax.patch.set_facecolor('gray')
    ax.set_aspect('equal', 'box')
    ax.xaxis.set_major_locator(plt.NullLocator())
    ax.yaxis.set_major_locator(plt.NullLocator())

    for (x, y), w in np.ndenumerate(matrix):
        color = 'white' if w > 0 else 'black'
        size = np.sqrt(np.abs(w))
        rect = plt.Rectangle([x - size / 2, y - size / 2], size, size,
                             facecolor=color, edgecolor=color)
        ax.add_patch(rect)

    ax.autoscale_view()
    ax.invert_yaxis()


def test_compute_regression_coeffs():
    sim_data = generate_simulated_data(20, 10, 200, 5, 2, 0.01, 1.0)
    C, D = sim_data["group_interactions"], sim_data["reg_coeffs"]
    Z, Y = sim_data["interactions"], sim_data["mrna_expression"]
    X = sim_data["mirna_expression"]

    start_time = time.time()
    D_est = compute_regression_coeffs(Z, Y, X, 0.01, 1.0, method="LBFGS")
    end_time = time.time()

    print "Computed in %f sec" % (end_time - start_time)

    diff = np.abs(D_est - D)
    max_diff = np.max(diff)
    argmax_diff = np.unravel_index(np.argmax(diff), D.shape)

    min_diff = np.min(diff)
    argmin_diff = np.unravel_index(np.argmin(diff), D.shape)

    print "Maximum difference %.3f. Actual %.3f. Estimated %.3f" % (
        max_diff, D[argmax_diff], D_est[argmax_diff])

    print "Minimum difference %.3f. Actual %.3f. Estimated %.3f" % (
        min_diff, D[argmin_diff], D_est[argmin_diff])

    f, axarr = plt.subplots(1, 2)

    hinton(D, ax=axarr[0])
    axarr[0].set_title(r'Original matrix ($\mathbf{D}$)')

    hinton(D_est, ax=axarr[1])
    axarr[1].set_title('Estimated matrix ($\hat{\mathbf{D}}$)')
    # plt.savefig('plot_regression_coeffs.pdf', bbox_inches='tight')
    plt.show()


if __name__ == '__main__':
    sys.exit(test_compute_regression_coeffs())