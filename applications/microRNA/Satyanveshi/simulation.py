__author__ = 'asish'
"""
    Generate simulated data based on the model
"""
import numpy as np


def generate_simulated_data(M, R, N, KM, KR, alpha, beta, noise=True,
                            seed=17):
    """
    :param M: # mRNAs
    :param R: # miRNAs
    :param N: # samples
    :param KM: # mRNA groups
    :param KR: # miRNA groups
    :param alpha:
    :param beta:
    :return:
    """
    # Fix the seed so that we get deterministic results each time we run.
    np.random.seed(seed=seed)
    # Generate base expression values for mRNAs.
    mu = np.random.uniform(0, 5.0, M)

    # Generate per sample variance
    s = np.repeat(0.1, N)

    # Gaussian error
    e = np.zeros(N)
    if noise:
        e = np.random.randn(N)
        e = e * np.sqrt(s)

    # Generate miRNA expression values.
    X = np.random.randn(N, R) * 2

    # Generate regression coefficients.
    D = np.random.randn(M, R)

    # Generate the group interaction matrix
    C = np.random.randint(2, size=(KM, KR))

    alpha_M = np.random.uniform(1, 100, KM)
    alpha_R = np.random.uniform(1, 100, KR)

    Y = np.zeros((N, M))

    # Assign miRNA to groups in each sample
    V = np.zeros((N, KR , R))
    for r in range(R):
        pi_r = np.random.dirichlet(alpha_R)
        V[:, :, r] = np.random.multinomial(1, pi_r, N)

    # Assign mRNAs to groups in each sample
    U = np.zeros((N, KM , M))
    for m in range(M):
        pi_m = np.random.dirichlet(alpha_M)
        U[:, :, m] = np.random.multinomial(1, pi_m, N)

    for m in range(M):
        for n in range(N):
            u_mn = U[n, :, m]
            w_mn = D[m, :] * X[n, :]
            Y[n,m] = mu[m] - u_mn.T.dot(C).dot(V[n]).dot(w_mn) + e[n]

    pz = (beta * D**2 + alpha)/(1 + (beta * D**2 + alpha))
    Z = np.random.binomial(1, pz)
    return dict(group_interactions=C,
                reg_coeffs=D,
                interactions=Z,
                mrna_expression=Y,
                mirna_expression=X,
                base_mrna_expression=mu,
                mrna_memberships=U,
                mirna_memberships=V)
