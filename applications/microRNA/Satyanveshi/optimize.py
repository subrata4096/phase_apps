"""
    Maximizes the Evidence Lower Bound (ELBO) using BFGS

N: Number of samples
M: Number of mRNAs
R: Number of miRNAs
KM: Number of mRNA modules
KR: Number of miRNA modules

Data:
    X: (N x R) Matrix of miRNA expression values.
    Y: (N x M) Matrix of mRNA expression values.
    Z: (M x R) Matrix of prior interaction probabilities.

Parameters:
    mu: (M x 1) vector of base expression values for mRNAs.
    sigmasq: (N x 1) vector of per sample variance
    C: (KM x KR) matrix of module interactions
    D: (M x R) matrix of regression coefficients
    alphaM: (KM x 1) Vector of mRNA Dirichlet parameters
    alphaR: (KR x 1) Vector of miRNA Dirichlet parameters

Variational Parameters:
    GammaM: Tensor of dimension (N x KM x M), where GammaM[n] is the matrix for the
        n-th sample. Each matrix has size (KM x M) and represent the KM dimensional
        multinomial mean parameters for the M mRNA.
    LambdaM: (KM x M) Matrix of Dirichlet parameters
    GammaR: Tensor of dimension (N x KR x R), where GammaR[n] is the matrix for the
        n-th sample. Each matrix has size (KR x R) and represents the KR
        dimensional multinomial parameters for the R miRNA.
    LambdaR: (KR x R) Matrix of Dirichlet parameters for the miRNAs.
"""
import numpy as np
from scipy.special import digamma, gamma


# Global variables
M = -1  # mRNAs
N = -1  # Samples
R = -1   # miRNAs
KM = -1   # mRNA modules
KR = -1   # miRNA modules


def construct_rho_matrix(Lambda):
    # Lambda is a matrix of Dirichlet parameters.
    lambda0 = np.sum(Lambda, axis=0)
    return digamma(Lambda) - digamma(lambda0)


def beta_fun(alpha):
    """
        Beta functions: prod(gamma(alpha))/gamma(sum(alpha))
    """
    return np.prod(gamma(alpha))/gamma(sum(alpha))


def log_beta_fun(alpha):
    """
        Beta functions: prod(gamma(alpha))/gamma(sum(alpha))
    """
    return np.sum(np.log(gamma(alpha))) - np.log(gamma(sum(alpha)))


def elbo_1(Y, X, mu, sigmasq, GammaM, GammaR, C, D):
    """
    First term of the ELBO.
    """
    W = np.einsum('ij,kj->ikj', D, X)
    s = 0.0    
    P = np.dot(C, GammaR)
    for m in range(M):
        for n in range(N):
            Temp = np.dot(np.dot(C.T, np.diag(GammaM[n][:, m])), C)
            K_mn = np.dot(np.dot(GammaR[n].T, Temp), GammaR[n]) 
            s -= np.dot(np.dot(W[m,n].T, K_mn), W[m,n])
            
            s -= (2*(Y[n,m] - mu[m])**2)*(GammaM[n][:,m].T.dot(P).dot(W[m,n]))

            s -= ((Y[n,m] - mu[m])**2)/(2*(sigmasq[n])**2)

            s -= np.log((2*(sigmasq[n])**2)*np.pi)/2.0
    return s


def elbo_3(GammaM, PM):
    return np.sum([np.sum(np.einsum('...i,...i', GammaM[i], PM)) for i in range(N)])


def elbo_4(alphaM, PM):
    return np.dot((alphaM - 1).T, PM) - M * log_beta_fun(alphaM)


def elbo_5(GammaR, PR):
    return np.sum([np.sum(np.einsum('...i,...i', GammaR[i], PR)) for i in range(N)])


def elbo_6(alphaR, PR):
    return np.dot((alphaR - 1).T, PR) - R * log_beta_fun(alphaR)


def elbo_7(GammaM):
    return (- np.sum([np.sum(np.einsum('i...,i...', GammaM[i], np.log(GammaM[i]))) 
            for i in range(N)]))


def elbo_8(LambdaM, PM):
    return (np.sum(np.apply_along_axis(log_beta_fun, 0, LambdaM)) - 
            np.sum(np.einsum('i...,i...', (LambdaM - 1), PM)))


def elbo_9(GammaR):
    return (- np.sum([np.sum(np.einsum('i...,i...', GammaR[i], np.log(GammaR[i]))) 
        for i in range(N)]))


def elbo_10(LambdaR, PR):
    return (np.sum(np.apply_along_axis(log_beta_fun, 0, LambdaR)) - 
            np.sum(np.einsum('i...,i...', (LambdaR - 1), PR)))

    
def flatten_parameters(param_map):    
    """
        Flatten a map of the following parameters into an array
        Parameters:
            mu: (M x 1) vector of base expression values for mRNAs.
            sigmasq: (N x 1) vector of per sample variance
            C: (KM x KR) matrix of module interactions
            alphaM: (KM x 1) Vector of mRNA Dirichlet parameters
            alphaR: (KR x 1) Vector of miRNA Dirichlet parameters

        Variational Parameters:
            GammaM: Tensor of dimension (N x KM x M), where GammaM[n] is the matrix for the
                n-th sample. Each matrix has size (KM x M) and represent the KM dimensional
                multinomial mean parameters for the M mRNA.
            LambdaM: (KM x M) Matrix of Dirichlet parameters
            GammaR: Tensor of dimension (N x KR x R), where GammaR[n] is the matrix for the
                n-th sample. Each matrix has size (KR x R) and represents the KR
                dimensional multinomial parameters for the R miRNA.
            LambdaR: (KR x R) Matrix of Dirichlet parameters for the miRNAs.
    """
    params = []
    bounds = []
    params.append(param_map['mu'])
    bounds.append()
    params.append(param_map['sigmasq'])
    params.append(param_map['C'].flatten())
    params.append(param_map['alphaM'])
    params.append(param_map['alphaR'])
    params.append(param_map['GammaM'].flatten())
    params.append(param_map['LambdaM'].flatten())
    params.append(param_map['GammaR'].flatten())
    params.append(param_map['LambdaR'].flatten())

    params = np.hstack(params)
    print "Created parameter array of size:", params.shape
    return params

    
def fold_parameters(param_arr):
    """
        Create Map of Params from array. Inverse of flatten_parameters.
    """
    offset = 0
    mu = param_arr[offset: (offset + M)]
    offset += mu.size
    
    sigmasq = param_arr[offset: (offset + N)]
    offset += sigmasq.size
    
    C = param_arr[offset: (offset + KM * KR)].reshape((KM, KR))
    offset += C.size

    alphaM = param_arr[offset: (offset + KM)]
    offset += alphaM.size

    alphaR = param_arr[offset: (offset + KR)]
    offset += alphaR.size

    GammaM = param_arr[offset: (offset + N * KM * M)].reshape((N, KM, M))
    offset = GammaM.size

    LambdaM = param_arr[offset: (offset + KM * M)].reshape((KM, M))
    offset += LambdaM.size

    GammaR = param_arr[offset: (offset + N * KR * R)].reshape((N, KR, R))
    offset += GammaR.size

    LambdaR = param_arr[offset: (offset + KR * R)].reshape((KR, R))

    return dict(mu=mu, sigmasq=sigmasq, C=C, alphaM=alphaM, alphaR=alphaR,
                GammaM=GammaM, LambdaM=LambdaM, GammaR=GammaR, LambdaR=LambdaR)


def setup_model_parameters(MODEL_PARAMS):
    global M
    M = MODEL_PARAMS['M']
    global R
    R = MODEL_PARAMS['R']
    global N
    N = MODEL_PARAMS['N']
    global KM
    KM = MODEL_PARAMS['KM']
    global KR
    KR = MODEL_PARAMS['KR']


def init_parameters(**MODEL_PARAMS):
    setup_model_parameters(MODEL_PARAMS)

    mu = np.zeros(M)
    sigmasq = np.ones(N)
    C = np.ones((KM, KR))
    alphaM = np.ones(KM)
    alphaR = np.ones(KR)
    GammaM = np.repeat(1.0, (N * KM * M)).reshape((N, KM, M))
    LambdaM = np.ones((KM, M))
    GammaR = np.repeat(1.0, (N * KR * R)).reshape((N, KR, R))
    LambdaR = np.ones((KR, R))
    return dict(mu=mu, sigmasq=sigmasq, C=C, alphaM=alphaM, alphaR=alphaR,
                GammaM=GammaM, LambdaM=LambdaM, GammaR=GammaR, LambdaR=LambdaR)



