import unittest
from Satyanveshi.simulation import generate_simulated_data
from Satyanveshi.variational_inference import *
from check_gradient import check_gradient
import copy
import scipy.stats


M = 20  # Number of mRNAs.
R = 9  # Number of miRNAs.
N = 10  # Number of samples.
KM = 5  # Number of mRNA groups.
KR = 2  # Number of miRNA groups.


class GradientCheckTestCase(unittest.TestCase):
    def setUp(self):
        sim_data = generate_simulated_data(M, R, N, KM, KR, 0.1, 1.0, seed=17)
        (self.C, self.D, self.Y, self.X, self.mu, self.U, self.V) = (
            sim_data["group_interactions"], sim_data["reg_coeffs"],
            sim_data["mrna_expression"], sim_data["mirna_expression"],
            sim_data["mrna_base_expression"], sim_data["mrna_memberships"],
            sim_data["mirna_memberships"]
        )

        self.GammaR = generate_gamma(N, KR, R)
        self.GammaM = generate_gamma(N, KM, M)
        self.LambdaM = generate_lambda_matrix(KM, M)
        self.LambdaR = generate_lambda_matrix(KR, R)
        self.alpha_M = np.random.uniform(1, 5, KM)
        self.alpha_R = np.random.uniform(1, 5, KR)

    def test_grad_Kmn(self):
        m, n = np.random.randint(0, M), np.random.randint(0, N)
        w_mn = self.D[m, :] * self.X[n, :]
        GammaR_n = self.GammaR[n]
        gamma_mn = self.GammaM[n, :, m]
        match = check_gradient(obj_Kmn, grad_Kmn, gamma_mn, args=[GammaR_n, w_mn, self.C])
        self.assertTrue(match, "Gradient doesnt match for K_mn with respect to gamma_mn")

    def test_gradient_gamma_rn(self):
        r, n = np.random.randint(0, R), np.random.randint(0, N)
        w_n = self.D * self.X[n, :]
        rho_r = rho_func(self.LambdaR[:,r])
        y_n = self.Y[n, :]
        mu = np.ones(M)
        match =  check_gradient(objective_func_gamma_rn, gradient_gamma_rn, self.GammaR[n, :, r],
                               args=(r, self.GammaR[n], y_n, mu, self.GammaM[n,:,:],
                                     self.C, w_n, rho_r, KR, 1.0))

        self.assertTrue(match, "Gradient doesnt match for Gamma_rn")

    def test_gradient_lambda_m(self):
        m = np.random.randint(0, M)
        Gamma_m = self.GammaM[:, :, m]
        match = check_gradient(objective_func_lambda, gradient_lambda,
                               self.LambdaM[:, KM],
                               args=[Gamma_m, self.alpha_M])

        self.assertTrue(match, "Gradient doesnt match for lambda_m")

    def test_gradient_lambda_r(self):
        r = np.random.randint(0, R)
        Gamma_r = self.GammaR[:, :, r]
        match = check_gradient(objective_func_lambda, gradient_lambda,
                               self.LambdaR[:, KR],
                               args=[Gamma_r, self.alpha_R])

        self.assertTrue(match, "Gradient doesnt match for lambda_r")

    def test_gradient_alphaM(self):
        RhoM = construct_rho_matrix(self.LambdaM)
        match = check_gradient(objective_func_alpha, gradient_func_alpha,
                               self.alpha_M, args=[RhoM])
        self.assertTrue(match, "Gradient doesnt match for alphaM")

    def test_gradient_alphaR(self):
        RhoR = construct_rho_matrix(self.LambdaR)
        match = check_gradient(objective_func_alpha, gradient_func_alpha,
                               self.alpha_R, args=[RhoR])
        self.assertTrue(match, "Gradient doesnt match for alphaR")

    def test_gradient_log_beta_func(self):
        match = check_gradient(log_beta_fun, rho_func, np.random.rand(5))
        self.assertTrue(match, "Gradient doesnt match log of Beta function")

    def test_gradient_d_mr(self):
        m = np.random.randint(0, M)
        r = np.random.randint(0, R)
        S_mn = np.array([self.C.T.dot(np.diag(self.GammaM[n, :, m])).dot(self.C) for n in range(N)])
        K_mr = np.einsum('...i,...ij,...jk', self.GammaR[:, :, r], S_mn, self.GammaR) # N x R
        match = check_gradient(obj_func_dmr, grad_dmr, np.random.randn(1),
                               args=[r, self.Y[:, m], self.mu[m], self.GammaM[:, :, m],
                                     self.GammaR, K_mr, self.C,  self.X, self.D[m, :],
                                     np.ones(N)])
        self.assertTrue(match, "Gradient doesnt match for d_mr")


class VariationalUpdateTestCase(unittest.TestCase):
    def setUp(self):
        sim_data = generate_simulated_data(M, R, N, KM, KR, 0.1, 1.0)
        (self.C, self.D, self.Y, self.X, self.mu) = (
            sim_data["group_interactions"], sim_data["reg_coeffs"],
            sim_data["mrna_expression"], sim_data["mirna_expression"],
            sim_data["mrna_base_expression"]
        )
        self.variational_inference = VariationalInference(M, R, N, KM, KR, (X, Y))

    def test_update_gamma_mn(self):
        m = np.random.randint(0, M)
        n = np.random.randint(0, N)
        RhoM = construct_rho_matrix(self.variational_inference.params.LambdaM)
        gamma_mn = self.variational_inference.update_gamma_mn(m, n, RhoM)

        print "Updated gamma_mn:", gamma_mn.round(2)
        self.assertTrue(np.all(gamma_mn <= 1) and np.all(gamma_mn >= 0), "gamma_mn is not between 0 and 1.")
        self.assertAlmostEquals(np.sum(gamma_mn), 1.0)

    def test_update_gamma_rn(self):
        r = 0
        n = 0
        RhoR = construct_rho_matrix(self.variational_inference.params.LambdaR)
        gamma_rn = self.variational_inference.update_gamma_rn(r, n, RhoR)
        # Check that gamma_rn is a valid multinomial mean parameter i.e. the elements sum to 1.
        print "Updated gamma_rn:", gamma_rn
        self.assertTrue(np.all(gamma_rn <= 1) and np.all(gamma_rn >= 0), "gamma_rn is not between 0 and 1.")
        self.assertAlmostEquals(np.sum(gamma_rn), 1.0)

    def test_update_lambda_m(self):
        m = np.random.randint(0, M)
        updated_lambda_m = self.variational_inference.update_lambda_m(m)
        print "Lambda_m after updating:", updated_lambda_m
        grad_lambda_m = gradient_lambda(updated_lambda_m, self.variational_inference.params.GammaM[:, :, m],
                                        self.variational_inference.params.alpha_M)
        print "Magnitude of grad of updated Lambda_m:", grad_lambda_m.dot(grad_lambda_m).round(6)
        self.assertTrue(updated_lambda_m.size == KM)
        self.assertAlmostEquals(grad_lambda_m.dot(grad_lambda_m), 0.0, 6)

    def test_update_lambda_r(self):
        r = np.random.randint(0, R)
        updated_lambda_r = self.variational_inference.update_lambda_r(r)
        print "Lambda_r after updating:", updated_lambda_r
        grad_lambda_r = gradient_lambda(updated_lambda_r, self.variational_inference.params.GammaR[:, :, r],
                                        self.variational_inference.params.alpha_R)
        print "Magnitude of grad of updated Lambda_r:", grad_lambda_r.dot(grad_lambda_r).round(6)
        self.assertAlmostEquals(grad_lambda_r.dot(grad_lambda_r), 0.0, 6)

    def test_update_alpha(self):
        RhoR = construct_rho_matrix(self.variational_inference.params.LambdaR)
        RhoM = construct_rho_matrix(self.variational_inference.params.LambdaM)
        updated_alphaR = self.variational_inference.update_alphaR(RhoR)
        updated_alphaM = self.variational_inference.update_alphaM(RhoM)
        print "Updated alphaR:", updated_alphaR.round(2)
        print "Updated alphaM:", updated_alphaM.round(2)
        grad_alphaM = gradient_func_alpha(updated_alphaM, RhoM)
        grad_alphaR = gradient_func_alpha(updated_alphaR, RhoR)
        print "Magnitude of gradient at updated_alphaM:", grad_alphaM.dot(grad_alphaM).round(6)
        print "Magnitude of gradient at updated_alphaR:", grad_alphaR.dot(grad_alphaR).round(6)
        self.assertAlmostEquals(grad_alphaM.dot(grad_alphaM), 0.0, 5)
        self.assertAlmostEquals(grad_alphaR.dot(grad_alphaR), 0.0, 5)

    def test_update_mu_m(self):
        m = np.random.randint(0, M)
        mu_m = self.variational_inference.update_mu_m(m)
        # Can't do much here apart from testing that updating mu didn't fail.
        print "Updated mu_m:", mu_m.round(2)

    def test_update_var_n(self):
        n = np.random.randint(0, N)
        var_n = self.variational_inference.update_var_n(n)
        # We can't do much here apart from testing that the variance is positive.
        print "Updated var_n:", var_n
        self.assertTrue(var_n > 0, "Variance is <= 0.")


class ParamTestCase(unittest.TestCase):
    def test_copy(self):
        param1 = Params(M, R, N, KM, KR)
        param2 = copy.deepcopy(param1)
        max_diff = np.max(np.abs(param2.GammaR[0, :, 0] - param1.GammaR[0, :, 0]))
        self.assertAlmostEquals(max_diff, 0)

        # Modify a parameter
        param1.GammaR[0, :, 0] = np.random.randn(KR)
        max_diff = np.max(np.abs(param2.GammaR[0, :, 0] - param1.GammaR[0, :, 0]))
        self.assertNotAlmostEquals(max_diff, 0)


def empirical_expectation(func_x, distribution_x, args_func, args_dist,
                          nsamples = 1000):
    """
    Compute the empirical expectation of a scalar valued
    function of random variable X by drawing samples from the given distribution.
    :param func_x: Function of random variable X.
    :param distribution_x: Distribution over random variable X.
    :param args_func: Arguments for @func_x
    :param args_dist: Arguments for @distribution_x
    :param nsamples: Number of samples to draw.
    :return: Expected value.
    """
    f_dist = partial(distribution_x, *args_dist)
    samples = f_dist(nsamples)
    func_val = np.apply_along_axis(lambda x: func_x(x, *args_func), 1, samples)
    if len(func_val.shape) != 1:
        raise Exception("Function is not scalar valued.")

    return np.mean(func_val)


def entropy_7(u_mn, gamma_mn):
    return u_mn.dot(np.log(gamma_mn))


def exp_entropy_7(gamma_mn):
    return gamma_mn.dot(np.log(gamma_mn))


def entropy_8(pi_m, lambda_m):
    return - np.log(scipy.stats.dirichlet.pdf(pi_m, lambda_m))


def exp_entropy_8(lambda_m):
    rho_m = rho_func(lambda_m)
    return log_beta_fun(lambda_m) - (lambda_m - 1).dot(rho_m)


def term_L1(x, y_mn, mu_m, C, w_mn):
    u_mn = x[0: KM]
    V_n = x[KM:].reshape((R, KR)).T
    mu_p_mn = mu_m - u_mn.dot(C).dot(V_n).dot(w_mn) # mu prime
    return (y_mn - mu_p_mn)**2


def distr_L1(gamma_mn, GammaR_n, size):
    us = np.random.multinomial(1, gamma_mn, size)
    Vs = [np.random.multinomial(1, GammaR_n[:, i], size)
          for i in range(GammaR_n.shape[1])]
    Vs = np.hstack(Vs)
    return np.hstack((us, Vs))


def cov(x, S):
    return x.dot(np.diag(S))


def exp_term_L1(gamma_mn, GammaR_n, y_mn, mu_m, C, w_mn):
    f = (y_mn - mu_m)**2
    f += 2 * (y_mn - mu_m) * gamma_mn.dot(C).dot(GammaR_n).dot(w_mn)
    S_mn = (C.T).dot(np.diag(gamma_mn)).dot(C)
    K_mn = GammaR_n.T.dot(S_mn).dot(GammaR_n)
    np.fill_diagonal(K_mn, GammaR_n.T.dot(np.diag(S_mn)))
    # K_mn_diag = np.diag(np.apply_along_axis(lambda x: cov(x, S_mn), 0, GammaR_n))
    # K_mn += K_mn_diag
    f += w_mn.dot(K_mn).dot(w_mn)
    return f


class ExpectationTestCase(unittest.TestCase):
    def setUp(self):
        np.random.seed(13)

    def test_exp_entropy_7(self):
        gamma_mn = np.random.rand(5)
        gamma_mn /= np.sum(gamma_mn)
        em_exp = empirical_expectation(entropy_7, np.random.multinomial,
                                       args_func=[gamma_mn],
                                       args_dist=[1, gamma_mn],
                                       nsamples=1000)

        true_exp = exp_entropy_7(gamma_mn)
        self.assertAlmostEquals(em_exp, true_exp, places=1,
                                msg="Empirical expectation of entropy term 7 didnt match.")

    def test_exp_entropy_8(self):
        lambda_m = np.random.rand(5) + 1.0
        em_exp = empirical_expectation(entropy_8, np.random.dirichlet,
                                       args_func=[lambda_m], args_dist=[lambda_m])
        true_exp = exp_entropy_8(lambda_m)
        msg = "Empirical expectation of entropy term 8 didnt match."
        msg += " True: %.1f Empirical: %.1f" %(em_exp, true_exp)
        self.assertAlmostEquals(em_exp, true_exp, places=1,
                                msg=msg)

    def test_exp_term1(self):
        gamma_mn = np.random.rand(KM)
        gamma_mn /= np.sum(gamma_mn)
        GammaR_n = generate_gamma(1, KR, R)[0]

        y_mn = 2.0
        mu_m = 1.0
        C = np.random.randint(0, 2, size=KM * KR).reshape((KM, KR))
        w_mn = np.random.randn(R) * 2

        em_exp = empirical_expectation(term_L1, distr_L1,
                                       args_func=[y_mn, mu_m, C, w_mn],
                                       args_dist=[gamma_mn, GammaR_n],
                                       nsamples=2000)
        true_exp = exp_term_L1(gamma_mn, GammaR_n, y_mn, mu_m, C, w_mn)
        msg = "Empirical expectation of term 1 didnt match."
        msg += " True: %.1f Empirical: %.1f" %(em_exp, true_exp)
        self.assertAlmostEquals(em_exp, true_exp, places=1,
                                msg=msg)

if __name__ == '__main__':
    unittest.main()