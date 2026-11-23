"""
Unit tests for distributions module
"""
import unittest
try:
    from unittest.mock import patch
except ImportError:
    from mock import patch
import pytest
import random
from os import path
import numpy as np
import parameterized as ptzd
from scipy.stats import invweibull, lognorm, kstest
from scipy.special import ndtr, ndtri
from psrpoppy import distributions as dists
from psrpoppy.pulsar import Pulsar

class test_ppdot_cohen27_posterior_draws(unittest.TestCase):

    def setUp(self):
        self.seed = 12345
        self.nsamp = 10000

    def test_npz_parameter_names_match_expected(self):
        npz = "correlated_invWeibull-P_lognormal-Pdot_multivarnorm_approx.npz"
        dat = np.load(path.join(dists.__dir__, npz))
        np.testing.assert_array_equal(
            dat["parameter_names"],
            np.array(['k_P', 'theta_nu_ms', 'mu_logPdot', 'sigma_logPdot', 'rho'])
            )
        
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

@ptzd.parameterized_class(
    ("k_p", "theta_nu", "mu_logpdot", "sigma_logpdot","rho"),
    [(2.142, 0.326, -46.344, 1.269, 0.290),]
                          )
class test_ppdot_cohen27(unittest.TestCase):

    def setUp(self):
        self.seed = 12345
        self.nsamp = 100000
        self.p, self.pdot = dists.ppdot_cohen27(
            self.k_p, self.theta_nu,
            self.mu_logpdot, self.sigma_logpdot,
            self.rho,
            n_samp=self.nsamp,
            rng=np.random.RandomState(self.seed)
        )
        # KS test statistic should be within DVK inequality bounds
        dvk_alpha = 1e-6 # false failure probability
        self.ks_limit = np.sqrt(np.log(2.0 / dvk_alpha) / (2.0 * self.nsamp))
        
    def test_P_marginal_vs_scipy_kstest(self):
        
        p_dist = invweibull(c=self.k_p, scale=1.0 / self.theta_nu)

        ks_result = kstest(self.p, p_dist.cdf)
        self.assertLess(ks_result.statistic, self.ks_limit)

    def test_Pdot_marginal_vs_scipy_kstest(self):
        pdot_dist = lognorm(s=self.sigma_logpdot, scale=np.exp(self.mu_logpdot))

        ks_result = kstest(self.pdot, pdot_dist.cdf)
        self.assertLess(ks_result.statistic, self.ks_limit)

    def test_P_marginal_vs_scipy_draws(self):
        # Reproduce the Gaussian inputs used by the sampler.
        z_p, z_2 = np.random.RandomState(self.seed).standard_normal((2, self.nsamp))

        p_reference = invweibull.ppf(ndtr(z_p),
                                     c=self.k_p,
                                     scale=1.0 / self.theta_nu)

        np.testing.assert_allclose(self.p, p_reference, rtol=1e-9, atol=0)

    def test_Pdot_correlation_recoveres_zscores(self):
        """
        test if correlation recovers expected Gaussian zscores
        """
        p_dist = invweibull(c=self.k_p, scale=1.0 / self.theta_nu)
        pdot_dist = lognorm(s=self.sigma_logpdot, scale=np.exp(self.mu_logpdot))

        u = p_dist.cdf(self.p) 
        v = pdot_dist.cdf(self.pdot)
        # Protect against numerical rounding to 0 or 1.
        eps = np.finfo(float).eps
        # recover normal z scores
        z_p = ndtri(np.clip(u, eps, 1.0 - eps))
        z_pdot = ndtri(np.clip(v, eps, 1.0 - eps))
        z_2_recovered = (z_pdot - self.rho * z_p) / np.sqrt(1.0 - self.rho**2)
        z_1_expected, z_2_expected = \
            np.random.RandomState(self.seed).standard_normal((2, self.nsamp))
        np.testing.assert_allclose(z_2_recovered, z_2_expected,
                                   rtol=1e-7, atol=1e-7)

class test_luminosity_fk06(unittest.TestCase):

    @ptzd.parameterized.expand([
        (1000., 1e-15, 0.18),
        (250., 4e-15, 0.18 * 16.),
    ])
    def test_fk06_table6_parameters_no_dither(self, period, pdot, lum1400):
        pulsar = Pulsar(period=period, pdot=pdot)

        with patch('psrpoppy.distributions.random.gauss', return_value=0.0):
            dists.luminosity_fk06(pulsar,
                                  alpha=-1.5,
                                  beta=0.5,
                                  gamma=0.18)

        self.assertAlmostEqual(pulsar.lum_1400, lum1400)
