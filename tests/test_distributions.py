"""
Unit tests for distributions module
"""
import unittest
import pytest
import random
import numpy as np
import parameterized as ptzd
from scipy.stats import invweibull, lognorm, kstest
from psrpoppy import distributions as dists

class test_ppdot_cohen27_posterior_draws(unittest.TestCase):

    def setUp(self):
        self.seed = 12345
        self.nsamp = 10000
        
    @ptzd.parameterized.expand([("k_P", 2.142086, 0.136349),
                                ("theta_nu_ms", 0.325553, 0.011591),
                                ("mu_logPdot", -46.344261, 0.098300),
                                ("sigma_logPdot", 1.269051, 0.074506),
                                ("rho", 0.290428, 0.072079),],
                               name_func=lambda fxn, n, par : "{}".format(ptzd.parameterized.to_safe_name(str(par.args[0]))).join(fxn.__name__.split("PARAM")))
    def test_PARAM_normalapprox_marginal_mean_std(self, param, mean, sigma):
        """
        Test marginal distributions of each parameter match sample means and std
        from transformed variables
        """
        N = self.nsamp
        joint_samp, names = dists.ppdot_cohen27_posterior_draws(
            fpath="correlated_invWeibull-P_lognormal-Pdot_multivarnorm_approx.npz",
            size=self.nsamp,
            rng=np.random.RandomState(self.seed))
        param_draws = np.squeeze(joint_samp[np.argwhere(names == param)])
        
        mean_sterr = sigma / np.sqrt(N)
        mu4 = 3 * sigma ** 4.
        std_sterr = np.sqrt((mu4 - (sigma ** 4.) * (N - 3) / (N - 1)) / N) / \
            (2 * sigma)
        np.testing.assert_allclose(
            np.mean(param_draws), mean,
            rtol=0, atol=5 * mean_sterr,
            err_msg="x==E({0}), y==np.mean({0})".format(param))
        np.testing.assert_allclose(
            np.std(param_draws, ddof=1), sigma,
            rtol=0, atol=5 * std_sterr,
            err_msg="x==s_{0}, y==np.std({0}, ddof=1)".format(param))
