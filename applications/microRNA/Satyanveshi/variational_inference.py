import numpy as np
import scipy as sp
import copy
from functools import partial
from scipy.optimize import fmin_l_bfgs_b
import traceback
import warnings
import pickle


def trigamma(x):
    return sp.special.polygamma(1, x)


def log_beta_fun(alpha):
    """
        Beta functions: prod(gamma(alpha))/gamma(sum(alpha))
    """
    return np.sum(np.log(sp.special.gamma(alpha))) - \
           np.log(sp.special.gamma(sum(alpha)))


def rho_func(x):
    return sp.special.digamma(x) - sp.special.digamma(np.sum(x))


def construct_rho_matrix(Lambda):
    # Lambda is a matrix of Dirichlet parameters.
    lambda0 = np.sum(Lambda, axis=0)
    return rho_func(Lambda) - rho_func(lambda0)


def __wrapper_func(func, args, x):
    return - func(x, *args)


def __maximize(func, grad_func, init_value, bounds=None, args=[]):
    f_obj = partial(__wrapper_func, func, args)
    f_grad = partial(__wrapper_func, grad_func, args)

    x_opt, _, info = fmin_l_bfgs_b(f_obj, init_value,
                                   fprime=f_grad, approx_grad=False, bounds=bounds,
                                   pgtol=0.001)
    if info['warnflag'] != 0:
        print info
        raise Exception("Failed to maximize function.")

    return x_opt


def obj_Kmn(gamma_mn, GammaR_n, w_mn, C):
    S = C.T.dot(np.diag(gamma_mn)).dot(C)
    K_mn =  GammaR_n.T.dot(S).dot(GammaR_n)
    np.fill_diagonal(K_mn, GammaR_n.T.dot(np.diag(S)))
    return w_mn.dot(K_mn).dot(w_mn)


def grad_Kmn(gamma_mn, GammaR_n, w_mn, C):
    def create_k(c_k):
        K = GammaR_n.T.dot(np.outer(c_k, c_k)).dot(GammaR_n)
        np.fill_diagonal(K, GammaR_n.T.dot(c_k))
        return K

    T = np.array([create_k(C[k, :]) for k in range(C.shape[0])])
    k_m_mn = np.einsum('i,...ij,j', w_mn, T, w_mn)
    return k_m_mn


def update_gamma_mn(y_mn, mu_m, C, GammaR_n, w_mn, lambda_m, rho_m, KM, var_n):
    t = 0
    gamma_mn = np.zeros(KM)
    try:
        k_m_mn = grad_Kmn(0, GammaR_n, w_mn, C)
        f_mn = C.dot(GammaR_n).dot(w_mn)
        t = (((-2) * (y_mn - mu_m)) * f_mn - k_m_mn)/(2 * var_n) + rho_m - 1
        gamma_mn = np.exp(t)
        gamma_mn /= np.sum(gamma_mn)
    except:
        print "Value of t", t
        print "Value of gamma_mn", gamma_mn
        print "Vaue of variance", var_n
        print "rho_m", rho_m
        print "f_mn", f_mn
        print "k_m_mn", k_m_mn
        print "mu_m", mu_m
        print "y_mn", y_mn
        raise Exception("Couldnt update gamma_mn")
    return gamma_mn


def objective_func_gamma_rn(gamma_rn, r, GammaR_n, y_n, mu, GammaM_n, C, W_n, rho_r, KR, var_n):
    GammaR_n[:, r] = gamma_rn
    ss = 0
    M = mu.size
    for m in range(M):
        S = C.T.dot(np.diag(GammaM_n[:, m])).dot(C)
        K_mn = GammaR_n.T.dot(S).dot(GammaR_n)
        np.fill_diagonal(K_mn, GammaR_n.T.dot(np.diag(S)))
        f = (2 * (y_n[m] - mu[m]) * (GammaM_n[:, m].T.dot(C).dot(gamma_rn)) * W_n[m, r])
        f += W_n[m].T.dot(K_mn).dot(W_n[m])
        f = - f/(2 * var_n)
        ss += f

    ss += gamma_rn.dot(rho_r)
    ss -= gamma_rn.dot(np.where(gamma_rn > 0.0, np.log(gamma_rn), 0.0))
    return ss


def gradient_gamma_rn(gamma_rn, r, GammaR_n, y_n, mu, GammaM_n, C, W_n, rho_r, KR, var_n):
    GammaR_n[:, r] = gamma_rn
    M = mu.size
    ss = np.zeros(gamma_rn.size)
    for m in range(M):
        S = C.T.dot(np.diag(GammaM_n[:, m])).dot(C)
        k_r_mn = 2 * W_n[m, r] * S.dot(GammaR_n).dot(W_n[m, :]) -\
                 (W_n[m, r]**2) * (2 * S.dot(GammaR_n[:, r]) - np.diag(S))

        f = 2 * ((y_n[m] - mu[m]) * (GammaM_n[:, m].T.dot(C)) * W_n[m, r])
        f += k_r_mn
        f = - f/(2 * var_n)
        ss += f

    ss += rho_r - np.log(gamma_rn) - 1 # This is vector 1 but because of numpy broadcasting
                                       # rules we can get away with a scalar.
    return ss


def update_gamma_rn(r, oGammaR_n, y_n, mu, GammaM_n, C, W_n, rho_r, KR, var_n):
    GammaR_n = np.copy(oGammaR_n)
    gamma_rn = GammaR_n[:, r]
    bounds = [(0.001, 1)] * gamma_rn.size
    gamma_rn = __maximize(objective_func_gamma_rn, gradient_gamma_rn,
                          gamma_rn, args=[r, GammaR_n, y_n, mu, GammaM_n, C, W_n, rho_r, KR, var_n],
                          bounds=bounds)
    gamma_rn /= np.sum(gamma_rn)
    return gamma_rn


def _jacobian_rho(lambd):
    J = np.diag(trigamma(lambd)) - np.ones((lambd.size, lambd.size)) * trigamma(np.sum(lambd))
    return J


def objective_func_lambda(lambda_x, Gamma_x, alphaX):
    rho_x = rho_func(lambda_x)
    f = (np.sum(Gamma_x, axis=0) + alphaX - lambda_x).T.dot(rho_x) + log_beta_fun(lambda_x)
    return f


def gradient_lambda(lambda_x, Gamma_x, alphaX):
    J = _jacobian_rho(lambda_x)
    f = J.dot(np.sum(Gamma_x, axis=0) + alphaX - lambda_x)
    return f


def update_lambda(lambda_x, Gamma_x, alphaX):
    bounds = [(0, None)] * lambda_x.size
    lambda_x = __maximize(objective_func_lambda, gradient_lambda, lambda_x,
                          args=[Gamma_x, alphaX], bounds=bounds)
    return lambda_x


def objective_func_alpha(alpha, Rho):
    n = Rho.shape[1]
    f = np.sum((alpha - 1).T.dot(Rho)) - n * log_beta_fun(alpha)
    return f


def gradient_func_alpha(alpha, Rho):
    n = Rho.shape[1]
    f = np.sum(Rho, axis=1) - n * rho_func(alpha)
    return f


def update_alpha(alpha, Rho):
    bounds = [(0, None)] * alpha.size
    alpha = __maximize(objective_func_alpha, gradient_func_alpha, alpha,
                       args=[Rho], bounds=bounds)
    return alpha


def update_mu_m(y_m, GammaM_m, C, GammaR, W_m, var):
    p = 1.0/var
    mu_m = 0
    for n in range(y_m.size):
        mu_m += (y_m[n] + GammaM_m[n, :].dot(C).dot(GammaR[n]).dot(W_m[n])) * p[n]
    mu_m /= np.sum(p)
    return mu_m


def update_var_n(y_n, mu, C, GammaM_n, GammaR_n, W_n):
    S_mn = np.array([C.T.dot(np.diag(GammaM_n[:, m])).dot(C) for m in range(y_n.size)])
    K_mn = np.einsum('ij,...jk,kl', GammaR_n.T, S_mn, GammaR_n)

    t = np.einsum('...i,ij,jk,k...', GammaM_n.T, C, GammaR_n, W_n.T)
    phi = (y_n - mu)**2 + 2 * (y_n - mu) * t
    phi += np.einsum('...i,...ij,...j', W_n, K_mn, W_n)
    var_n = np.mean(phi)
    var_n = max(var_n, 0.01)
    return var_n


def obj_func_dmr(d_mr, r, y_m, mu_m, GammaM_m, GammaR, K_mr, C, X, D_m, var):
    D_m[r] = d_mr
    W_m = D_m * X # N x R
    X_r = X[:, r]
    L = (y_m - mu_m)[:, np.newaxis] * GammaM_m
    R = (GammaR[:, :, r] * X_r[:, np.newaxis])
    f = np.einsum('i...,ij,...j', L.T, C, R)  * d_mr # This is a vector of size N
    f += np.einsum('...i, ...i', K_mr, W_m) * X_r * d_mr
    f -= (W_m[:, r]**2 * K_mr[:, r])/2.0
    f /= -var
    return np.sum(f)


def grad_dmr(d_mr, r, y_m, mu_m, GammaM_m, GammaR, K_mr, C, X, D_m, var):
    D_m[r] = d_mr
    W_m = D_m * X
    X_r = X[:, r]
    L = (y_m - mu_m)[:, np.newaxis] * GammaM_m
    R = (GammaR[:, :, r] * X_r[:, np.newaxis])
    f = np.einsum('i...,ij,...j', L.T, C, R) # This is a vector of size N
    f += np.einsum('...i, ...i', K_mr, W_m) * X_r
    f /= -var
    return np.array([np.sum(f)])


def update_dmr(d_mr, r, y_m, mu_m, GammaM_m, GammaR, K_mr, C, X, oD_m, var):
    D_m = np.copy(oD_m)
    d_mr = __maximize(obj_func_dmr, grad_dmr, d_mr,
                       args=[r, y_m, mu_m, GammaM_m, GammaR, K_mr, C, X, D_m, var])
    return d_mr


def generate_gamma(N, *shape):
    Gamma = [None] * N
    for n in range(N):
        Gamma[n] = np.ones(shape)
        Gamma[n] /= np.sum(Gamma[n], axis=0)
    return np.array(Gamma)


def generate_K_matrices(M, N, C, GammaM, GammaR):
    K = {}
    for m in range(M):
        for n in range(N):
            K[(m,n)] = GammaR[n].T.dot(C.T).dot(np.diag(GammaM[n,:,m])).dot(C).dot(GammaR[n])
    return K


def generate_lambda_matrix(*shape):
    return np.random.rand(*shape) + 1.0


def has_converged(x, prev_x, tolerance=0.1):
    max_diff = np.max(np.abs(x - prev_x))
    return max_diff <= tolerance


class VariationalInference:
    def __init__(self, M, R, N, KM, KR, data,
                 max_iterations=20, tolerance=0.0001,
                 damping_factor=0.7, maintain_history=False):

        """
        Initialize variational inference.
        :param M: # mRNAs
        :param R: # miRNAs
        :param N: # samples
        :param KM: # mRNA groups
        :param KR: # miRNA groups
        :param data: (X, Y) where X is N x R matrix of expression values of miRNAs
            and Y is N x M matrix of expression values of mRNAs.
        :param D: M x R matrix of regression coefficients.
        :param max_iterations: Maximum number of iterations of parameter updates to perform.
        :param tolerance: Tolerance for checking convergence of parameters. If the maximum absolute difference
            between previous and current value of a parameter is less than the tolerance value then the
            parameter has converged.
        """
        warnings.filterwarnings('error')
        self.M = M
        self.R = R
        self.N = N
        self.KM = KM
        self.KR = KR + 1
        self.C = np.ones((self.KM, self.KR))
        self.C[:, -1] = 0
        self.params = Params(self.M, self.R, self.N, self.KM, self.KR)
        self.param_hist = None
        if maintain_history:
            self.param_hist = ParamHistory(20)
            self.param_hist.add(self.params)
        self.damping_factor = damping_factor

        self.X, self.Y = data
        # W is (M x N x R)
        # self.W = np.einsum('i...,k...->ik...', self.D, self.X)
        self.converged = False
        self.max_iters = max_iterations

    def update_gamma_mn(self, m, n, RhoM):
        W_mn = self.X[n] * self.params.D[m]
        return update_gamma_mn(self.Y[n, m], self.params.mu[m], self.C, self.params.GammaR[n],
                               W_mn, self.params.LambdaM[:, m], RhoM[:, m],
                               self.KM, self.params.var[n])

    def update_gamma_rn(self, r, n, RhoR):
        W_n = self.X[n] * self.params.D
        return update_gamma_rn(r, self.params.GammaR[n], self.Y[n], self.params.mu, self.params.GammaM[n],
                                   self.C, W_n, RhoR[:, r], self.KR, self.params.var[n])

    def update_lambda_m(self, m):
        lambda_m = self.params.LambdaM[:, m]
        return update_lambda(lambda_m, self.params.GammaM[:, :, m], self.params.alpha_M)

    def update_lambda_r(self, r):
        lambda_r = self.params.LambdaR[:, r]
        return update_lambda(lambda_r, self.params.GammaR[:, :, r], self.params.alpha_R)

    def update_alphaR(self, RhoR):
        return update_alpha(self.params.alpha_R, RhoR)

    def update_alphaM(self, RhoM):
        return update_alpha(self.params.alpha_M, RhoM)

    def update_mu_m(self, m):
        W_m = self.params.D[m] * self.X
        return update_mu_m(self.Y[:, m], self.params.GammaM[:, :, m], self.C, self.params.GammaR,
                           W_m, self.params.var)

    def update_var_n(self, n):
        W_n = self.X[n] * self.params.D
        return update_var_n(self.Y[n], self.params.mu, self.C,
                            self.params.GammaM[n], self.params.GammaR[n], W_n)

    def update_dmr(self, m, r):
        S_mn = np.array([self.C.T.dot(np.diag(self.params.GammaM[n, :, m])).dot(self.C)
                         for n in range(self.N)])
        K_mr = np.einsum('...i,...ij,...jk',
                         self.params.GammaR[:, :, r], S_mn, self.params.GammaR) # N x R

        return update_dmr(self.params.D[m,r], r,
                          self.Y[:, m], self.params.mu[m],
                          self.params.GammaM[:, :, m],
                          self.params.GammaR, K_mr, self.C, self.X,
                          self.params.D[m], self.params.var)

    def update_param(self, new_param, old_param, damping_factor):
        return new_param * damping_factor + (1 - damping_factor) * old_param

    def iterate(self):
        converged = True
        new_params = copy.deepcopy(self.params)

        # Update the miRNA gamma parameters
        RhoR = construct_rho_matrix(self.params.LambdaR)
        for n in range(self.N):
	    #Manish : loop
	    #if(n%3==0):
		#continue
            for r in range(self.R):
                new_params.GammaR[n, :, r] = self.update_param(
                    self.update_gamma_rn(r, n, RhoR), new_params.GammaR[n, :, r],
                    self.damping_factor
                )
                conv = has_converged(new_params.GammaR[n, :, r],
                                     self.params.GammaR[n, :, r])
                converged = converged and conv

        # Update the mRNA gamma parameters
        RhoM = construct_rho_matrix(self.params.LambdaM)
        '''
        for n in range(self.N):
            for m in range(self.M):
                new_params.GammaM[n, :, m] = self.update_param(
                    self.update_gamma_mn(m, n, RhoM), new_params.GammaM[n, :, m],
                    self.damping_factor
                )
                conv = has_converged(new_params.GammaM[n, :, m],
                                    self.params.GammaM[n, :, m])
                converged = converged and conv
        '''
        # Update the mRNA lambda parameters
        for m in range(self.M):
            new_params.LambdaM[:, m] = self.update_param(
                self.update_lambda_m(m), new_params.LambdaM[:, m],
                self.damping_factor
            )
            conv = has_converged(new_params.LambdaM[:, m],
                                 self.params.LambdaM[:, m])
            converged = converged and conv

        # Update the miRNA lambda parameters
        for r in range(self.R):
            new_params.LambdaR[:, r] = self.update_param(
                self.update_lambda_r(r), new_params.LambdaR[:, r],
                self.damping_factor
            )
            conv = has_converged(new_params.LambdaR[:, r],
                                 self.params.LambdaR[:, r])
            converged = converged and conv

        # Update the mRNA dirichlet parameters
        new_params.alpha_M = self.update_param(
            update_alpha(self.params.alpha_M, RhoM), new_params.alpha_M,
            self.damping_factor
        )
        conv = has_converged(new_params.alpha_M, self.params.alpha_M)
        converged = converged and conv

        # Update the miRNA dirichlet parameters
        new_params.alpha_R = self.update_param(
                update_alpha(self.params.alpha_R, RhoR), new_params.alpha_R,
                self.damping_factor
        )
        conv = has_converged(new_params.alpha_R,
                             self.params.alpha_R)
        converged = converged and conv

        # Update regression coeffs
        for m in range(self.M):
            for r in range(self.R):
                new_params.D[m, r] = self.update_param(
                        self.update_dmr(m, r), new_params.D[m, r],
                        self.damping_factor
                )
                conv = has_converged(new_params.D[m, r], self.params.D[m, r])
                converged = converged and conv

        # Update mean mRNA parameter
        for m in range(self.M):
            new_params.mu[m] = self.update_param(
                    self.update_mu_m(m), new_params.mu[m],
                    self.damping_factor
            )
        conv = has_converged(new_params.mu, self.params.mu)
        converged = converged and conv

        # Update per-sample variance parameters
        for n in range(self.N):
            new_params.var[n] = self.update_param(
                self.update_var_n(n), new_params.var[n],
                self.damping_factor
            )
        conv = has_converged(new_params.var, self.params.var)
        converged = converged and conv
        
        self.params = new_params
        if self.param_hist is not None:
            self.param_hist.add(self.params)
        return converged

    def __call__(self, *args, **kwargs):
        if self.converged:
            print "Already converged."
            return

        for i in range(self.max_iters):
            try:
                self.converged = self.iterate()
            except Exception as e:
                print "Caught exception."
                traceback.print_exc()
                print
                print self.params
                return

            if self.converged:
                print "Converged in %d iterations" % i
                return
            print "Iteration %d done" % i
            # print self.params
        print "Failed to converged after %d iterations" % self.max_iters
        print self.params

    @staticmethod
    def save(obj, fname):
        pickle.dump(obj, fname)

    @staticmethod
    def load(fname):
        return pickle.load(fname)


class Params:
    def __init__(self, M, R, N, KM, KR):
        """
            Initializes parameters to default values.
        """
        self.GammaR = generate_gamma(N, KR, R)
        self.GammaM = generate_gamma(N, KM, M)
        self.LambdaM = generate_lambda_matrix(KM, M)
        self.LambdaR = generate_lambda_matrix(KR, R)
        self.alpha_M = np.ones(KM) * 0.1
        self.alpha_R = np.ones(KR) * 0.1
        self.mu = np.zeros(M)
        self.var = np.ones(N)
        self.D = np.zeros((M, R))

    def __str__(self):
        s = ""
        for k, v in self.__dict__.iteritems():
            s += "%s\n%s\n" % (k, v.round(2))
        return s

class ParamHistory:
    def __init__(self, length):
        self.params = [None] * length
        self.count = 0

    def add(self, param):
        if self.count < len(self.params):
            self.params[self.count] = param
            self.count += 1
        else:
            self.params[0:(self.count -  1)] = self.params[1:]
            self.params[-1] = param

    def get_history(self, param_key,
                    index_func=lambda x: x):
        hist = map(lambda param: index_func(param.__dict__[param_key]), self.params)
        return np.array(hist)
