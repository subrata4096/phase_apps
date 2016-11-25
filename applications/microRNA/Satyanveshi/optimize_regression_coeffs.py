__author__ = 'asish'
"""
Given prior interaction data Z (M x R), mRNA
expression data Y (N x M) mRNAs and
expression data X (N x R) for miRNAs.
Estimate the regression coefficient matrix D (M x R)
"""

import numpy as np
import numpy.linalg as la
from scipy.optimize import fmin_l_bfgs_b, fmin_ncg


def __objective_func(x, z_mr, y_m, x_r, mu_m, s, alpha, beta):
    d_mr = x[0]
    N = y_m.size
    p = 1.0/s
    P = np.diag(p)
    Sigma = np.diag(s)

    mu_mr = mu_m - (d_mr * x_r)
    diff = y_m - mu_mr
    ss = (- 0.5) * (N * np.log(2 * np.pi) + np.log(la.det(Sigma)) +
          diff.T.dot(P).dot(diff))

    ss += z_mr * np.log(beta * d_mr**2 + alpha) - np.log(1 + alpha + beta * d_mr**2)
    return np.array([-ss])


def __grad(x, z_mr, y_m, x_r, mu_m, s, alpha, beta):
    d_mr = x[0]
    p = 1.0/s
    P = np.diag(p)

    ss = - x_r.T.dot(P.dot(y_m) - mu_m * p + d_mr * (P.dot(x_r)))
    ss += (2 * z_mr * beta * d_mr)/(alpha + beta * d_mr**2) -\
         (2 * beta * d_mr)/(1 + alpha + beta * d_mr**2)

    return np.array([-ss])


def __hessian(x, z_mr, y_m, x_r, mu_m, s, alpha, beta):
    d_mr = x[0]
    p = 1.0/s
    P = np.diag(p)

    ss = x_r.T.dot(P).dot(x_r) + \
         ((alpha - beta * d_mr**2) * (2 * z_mr * beta))/((alpha + beta * d_mr**2)**2) -\
         ((1 + alpha - beta * d_mr**2) * (2 * z_mr * beta))/((1 + alpha + beta * d_mr**2)**2)

    return np.array([[-ss]])


def __maximize(z_mr, y_m, x_r, mu_m, s, alpha, beta, method="LBFGS"):
    if method == "LBFGS":
        d_mr, f_min, info = fmin_l_bfgs_b(__objective_func, np.zeros(1),
                                          fprime=__grad, approx_grad=False,
                                          args=(z_mr, y_m, x_r, mu_m, s, alpha, beta))
    elif method == "NEWTON":
        d_mr = fmin_ncg(__objective_func, np.zeros(1),
                     fprime=__grad, fhess=__hessian,
                     args=(z_mr, y_m, x_r, mu_m, s, alpha, beta), disp=0)
        # d_mr = r[0]
        f_min = 0
        info = {"warnflag": 0}

    return d_mr, f_min, info


def compute_regression_coeffs(Z, Y, X, alpha, beta, method="LBFGS"):
    # mRNA base expression
    mu = np.mean(Y, axis=0)

    # mRNA per sample variance
    s = np.var(Y, axis=1)

    M = Y.shape[1]
    R = X.shape[1]
    D = np.ones((M, R))

    for m in range(M):
        for r in range(R):

            d_mr, f_min, info = __maximize(Z[m,r], Y[:,m], X[:,r], mu[m], s, alpha, beta, method=method)
            if info['warnflag'] == 0:
                D[m, r] = d_mr
            else:
                D[m, r] = np.nan
                print "Error for mRNA %d, miRNA %d" % (m, r)
                print info
    return D
