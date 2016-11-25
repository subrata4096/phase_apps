#!/usr/bin/env python

from Satyanveshi.variational_inference import VariationalInference
from Satyanveshi.bayesian_posterior_check import BayesianPosteriorCheck
from Satyanveshi.simulation import generate_simulated_data
import sys
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import numpy as np
import argparse


# CMAP = plt.cm.PRGn
# CMAP = plt.cm.PuOr
CMAP = plt.cm.RdBu_r


M, R, N, KM, KR = 6, 4, 10, 2, 2
#M, R, N, KM, KR = 6, 4, 100, 2, 2


def parse_args():
    parser = argparse.ArgumentParser(description='Run Variational Inference.')
    parser.add_argument(
        "--save",
        help="Save variational inference object to the specified file",
        dest="save_fname", default=None
    )
    parser.add_argument(
        "--load", help="Load variational inference object from specified file",
        dest="load_fname", default=None
    )
    parser.add_argument("--plot", action="store_true",
                        help="plot results.")
    parser.add_argument("--save-plot", dest="save_plot", default=None,
                        help="File to save the plot to")

    args = parser.parse_args()
    return args


def plot_mrna_mirna_interaction(axarr, mrna_mem, mirna_mem, group_inter,
                                reg_coeffs, norm):
    """
    Plot mRNA-miRNA interactions in different samples.
    :param mrna_mem: per-sample mRNA memberships
    :param mirna_mem: per-sample miRNA memberships
    :param group_inter: matrix of group interactions
    :param reg_coeffs: Regression coefficients
    :param norm
    :return:
    """
    N, KM, M = mrna_mem.shape
    _, KR, R = mirna_mem.shape
    # D = DM
    # D -= np.mean(DM)
    # D /= (np.max(D) - np.min(D))
    # D *= -1.0
    X = np.einsum('...ji,jk,...kl', mrna_mem, group_inter, mirna_mem)
    for n in range(N):
        im =axarr[n].matshow(X[n] * reg_coeffs, cmap=CMAP, norm=norm)

    return im


def plot_results(var_inf, sim_data, args):
    D = sim_data["reg_coeffs"]
    fig, axarr = plt.subplots(2, N)
    vmax = max(np.max(D), np.max(var_inf.params.D))
    vmin = min(np.min(D), np.min(var_inf.params.D))
    norm = Normalize(vmin=vmin, vmax=vmax)
    mappable = ScalarMappable(norm=norm, cmap=CMAP)
    mappable.set_array(np.linspace(vmin, vmax, 100))

    plot_mrna_mirna_interaction(axarr[0],
                                sim_data["mrna_memberships"],
                                sim_data["mirna_memberships"],
                                sim_data["group_interactions"],
                                sim_data["reg_coeffs"],
                                norm)

    im = plot_mrna_mirna_interaction(axarr[1],
                                     var_inf.params.GammaM,
                                     var_inf.params.GammaR,
                                     var_inf.C,
                                     var_inf.params.D, norm)

    fig.subplots_adjust(bottom=0.3, right=0.8, top=0.7)
    cbar_ax = fig.add_axes([0.85, 0.3, 0.01, 0.4])

    fig.colorbar(mappable, cax=cbar_ax)
    if args.save_plot is not None:
        plt.savefig(args.save_plot, bbox_inches='tight')
    plt.show()


def test(args):
    sim_data = generate_simulated_data(M, R, N, KM, KR, 0.01, 1.0, seed=13)
    if args.load_fname is not None:
        with open(args.load_fname, 'rU') as fp:
            var_inf = VariationalInference.load(fp)
    else:
        var_inf = VariationalInference(M, R, N, KM, KR,
                                       (sim_data["mirna_expression"],
                                        sim_data["mrna_expression"]),
                                       max_iterations=1000)
        var_inf()

    if args.save_fname is not None:
        with open(args.save_fname, 'w') as fp:
            VariationalInference.save(var_inf, fp)

    bayesian_post_check = BayesianPosteriorCheck(sim_data["mrna_expression"],
                                                 sim_data["mirna_expression"],
                                                 var_inf.C,
                                                 var_inf.params)
    bayesian_post_check.fit()
    per_sample_pvals, pval = bayesian_post_check.pvalue_chi_sq()
    print "Per-sample pvalues:", per_sample_pvals.round(3)
    print "Overall pvalue:", round(pval, 3)

    if args.plot:
        plot_results(var_inf, sim_data, args)

    return var_inf


def main():
    args = parse_args()
    test(args)
    return 0

if __name__ == '__main__':
    sys.exit(main())
