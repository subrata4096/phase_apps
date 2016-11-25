"""
Use Bayesian Posterior Check to compute the
p-value of the observed data under the model.
"""
import numpy as np


def compute_mean(m, n, X, U, V, C, params):
    """
    Compute mean for the m-th mRNA in the n-th sample.
    :param m: mRNA
    :param n: sample
    :param X: miRNA expression
    :param U: mRNA membership vectors
    :param V: miRNA membership vectors
    :param C: Group interaction matrix
    :param params: instance of Params class.
    :return: mean.
    """
    u_mn = U[n, :, m]
    w_mn = params.D[m, :] * X[n, :]
    mu = params.mu[m] - u_mn.T.dot(C).dot(V[n]).dot(w_mn)
    return mu


def chi_sq_statistic(obs_Y, replicates, params):
    """
    The chi square statistics if (Y - E[Y])^2/Var[Y].
    The observed data Y (N x M) is the expression of M mRNAs across
    N samples.
    :param obs_Y: The observed mRNA expression.
    :param replicates: Samples and parameters drawn from posterior.
    :param params: Inferred params (Class Params).
    :return: per-sample p-values and overall pvalue
    """
    N, M = obs_Y.shape
    # counts the number of replicated samples
    # whose statistic exceed that of the
    # the observed data
    bad_count = np.zeros(N)
    overall_badcount = 0.0 # This is across all N samples.
    for (repl_data, repl_param) in replicates:
        d_repl_all = 0.0
        d_obs_all = 0.0
        for n in range(N):
            # compute statistic for replicate
            # repl_data is of dim N x M
            diff = repl_data[n, :] - repl_param['Mu'][n]
            d_repl = np.dot(diff, diff)/params.var[n]
            d_repl_all += d_repl
            # compute statistic for observed data
            # Note that the mean remains the same as that for
            # replicated data becuase we are trying to explain
            # the observed data in terms of the replicated params.
            diff = obs_Y[n, :] - repl_param['Mu'][n]
            d_obs = np.dot(diff, diff)/params.var[n]
            d_obs_all += d_obs
            if d_repl > d_obs:
                bad_count[n] += 1.0

        if d_repl_all > d_obs_all:
            overall_badcount += 1.0

    pvals = bad_count/replicates.num_replicates
    pval = overall_badcount/replicates.num_replicates
    return pvals, pval


def generate_replicate(N, M, R, KM, KR, X, C, params):
    Y = np.zeros((N, M))
    Mu = np.zeros((N, M)) # mean expression for each mRNA across N samples.

    # Assign miRNA to groups in each sample
    V = np.zeros((N, KR , R))
    for r in range(R):
        for n in range(N):
            V[n, :, r] = np.random.multinomial(1, params.GammaR[n, :, r], 1)

    # Assign mRNAs to groups in each sample
    U = np.zeros((N, KM , M))
    for m in range(M):
        for n in range(N):
            U[n, :, m] = np.random.multinomial(1, params.GammaM[n, :, m], 1)

    for m in range(M):
        for n in range(N):
            mu = compute_mean(m, n, X, U, V, C, params)
            Mu[n,m] = mu
            v = params.var[n]
            Y[n,m] = np.random.randn() * v + mu

    return Y, dict(Mu=Mu)


class Replicates:
    def __init__(self, num_replicates):
        self.num_replicates = num_replicates
        self.data_replicates = [None] * self.num_replicates
        self.param_replicates = [None] * self.num_replicates
        self._count = 0

    def __iter__(self):
        return zip(self.data_replicates, self.param_replicates).__iter__()

    def _check_count(self):
        if self._count >= self.num_replicates:
            raise Exception("Maximum number of replicates exceeded.")

    def add_replicate(self, data, param):
        self._check_count()
        self.data_replicates[self._count] = data
        self.param_replicates[self._count] = param
        self._count += 1


class BayesianPosteriorCheck:
    def __init__(self, Y, X, C, posterior_params,
                 num_replicates=500):
        self.X = X
        self.Y = Y
        self.C = C
        self.params = posterior_params
        self.N, self.M = Y.shape
        _, self.R = X.shape
        self.KM = posterior_params.GammaM.shape[1]
        self.KR = posterior_params.GammaR.shape[1]
        self.replicates = Replicates(num_replicates)
        self.initialized = False

    def fit(self):
        if self.initialized:
            return

        # Generate replicates
        print "Generating replicates from posterior distribution."
        for s in range(self.replicates.num_replicates):
            repl = generate_replicate(self.N, self.M, self.R,
                                      self.KM, self.KR, self.X,
                                      self.C, self.params)
            self.replicates.add_replicate(repl[0], repl[1])
        print "Done."
        self.initialized = True

    def pvalue_chi_sq(self):
        if not self.initialized:
            raise Exception("Call the method fit() before calling this function.")
        # Compute the pvalue for the chi squared test statistic
        return chi_sq_statistic(self.Y, self.replicates, self.params)




