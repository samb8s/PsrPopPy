"""
Integration tests for populate module
"""
import unittest
import pytest
import random
import numpy as np
import parameterized as ptzd
from psrpoppy import populate

class test_populate_generate_distributions(unittest.TestCase):

    def test_population_size_equals_ngen(self):
        pop = populate.generate(100,
                                nostdout=True)
        self.assertEqual(100, len(pop.population))

    @ptzd.parameterized.expand([(100,), (1000,), (10000,)],
                               name_func=lambda fxn, n, par : "_{}".format(ptzd.parameterized.to_safe_name(str(par.args[0]))).join(fxn.__name__.split("_npsrs")))
    def test_pDistType_lorimer12_peak_ngen_npsrs(self, npsrs):
        """
        Lorimer (2012, IAU) dist supposed to have peak at 3ms
        As implemented, its actually slightly more than 3
        """
        logpmin = 0.
        logpmax = 1.5
        dist = [1., 3., 5., 16., 9., 5., 5., 3., 2.]
        bins = np.logspace(logpmin, logpmax, len(dist) + 1)
        pop = populate.generate(npsrs,
                                pDistType='lorimer12',
                                nostdout=True)
        periods = [p.period for p in pop.population]
        expected_mode_bin = np.argmax(dist)
        histdata, _ = np.histogram(periods, bins=bins)
        self.assertEqual(np.argmax(histdata), expected_mode_bin)

    @ptzd.parameterized.expand([(1000,),],
                               name_func=lambda fxn, n, par : "_{}".format(ptzd.parameterized.to_safe_name(str(par.args[0]))).join(fxn.__name__.split("_npsrs")))
    def test_pDistType_lnorm_mean_std_ngen_npsrs(self, npsrs):
        """
        Test base-10 lnorm period sample mean and standard deviation match 
        pop.pmean and pop.psigma to 5 x standard error
        """
        random.seed(12345)
        mu_log = 2.7
        sig_log = -0.34
        mu_ln = np.log(10) * mu_log
        sig_ln = np.log(10) * sig_log
        pop = populate.generate(npsrs,
                                pDistType='lnorm',
                                pDistPars=[mu_log, sig_log],
                                nostdout=True)
        periods = np.array([p.period for p in pop.population])

        sampmean = (10 ** mu_log) * np.exp(0.5 * sig_ln ** 2.)
        sampstd = sampmean * np.sqrt(np.exp(sig_ln ** 2.) - 1)
        sampmean_sterr = sampstd / np.sqrt(npsrs)
        beta2 = np.exp(4 * sig_ln**2) +\
            2 * np.exp(3 * sig_ln**2) +\
            3 * np.exp(2 * sig_ln**2) - 3
        sampstd_sterr = sampstd / (2 * np.sqrt(npsrs)) * \
            np.sqrt(beta2 - (npsrs - 3.) / (npsrs - 1.))
        
        np.testing.assert_allclose(np.mean(periods), sampmean,
                                   rtol=0, atol=5 * sampmean_sterr,
                                   err_msg="x==E(P), y==np.mean(periods)")
        np.testing.assert_allclose(np.std(periods, ddof=1), sampstd,
                                   rtol=0, atol=5 * sampstd_sterr,
                                   err_msg="x==sP, y==np.std(periods, ddof=1)")

    @ptzd.parameterized.expand([("norm", [0., 1.]),
                                ("lnorm", [1e-3, 2.]),
                                ("lorimer12", [0., 0.]),
                                ("cc97", [0., 0.])],
                               name_func=lambda fxn, n, par : "{}".format(ptzd.parameterized.to_safe_name(str(par.args[0]))).join(fxn.__name__.split("_disttype_")))
    def test_pDistType__disttype__disallows_negative_periods(self,
                                                             pdisttype,
                                                             pdistpars):
        npsrs = 1000
        pop = populate.generate(npsrs,
                                pDistType=pdisttype,
                                pDistPars=pdistpars,
                                nostdout=True)
        periods = np.array([p.period for p in pop.population])
        np.testing.assert_array_less(np.zeros(npsrs), periods)
        
    @ptzd.parameterized.expand([(1000,),],
                               name_func=lambda fxn, n, par : "_{}".format(ptzd.parameterized.to_safe_name(str(par.args[0]))).join(fxn.__name__.split("_npsrs")))
    def test_pDistType_lnorm_e_mean_std_ngen_npsrs(self, npsrs):
        """
        Test base-e lnorm_e period sample mean and standard deviation match 
        pop.pmean and pop.psigma to 5 x standard error
        """
        self.skipTest("This should be for a base-e log-normal distribution")
        random.seed(12345)
        mu_ln = 1.5
        sig_ln = 0.58
        pop = populate.generate(npsrs,
                                pDistType='lnorm_e',
                                pDistPars=[mu_ln, sig_ln],
                                nostdout=True)
        periods = np.array([p.period for p in pop.population])
        sampmean = np.exp(mu_ln + 0.5 * (sig_ln ** 2.))
        sampstd = np.sqrt((np.exp(sig_ln ** 2.) -1) * \
                          np.exp(2 * mu_ln + sig_ln ** 2.))
        sampmean_sterr = sampstd / np.sqrt(npsrs)
        beta2 = np.exp(4 * sig_ln**2) + 2 * np.exp(3 * sig_ln**2) +\
             3 * np.exp(2 * sig_ln**2) - 3
        sampstd_sterr = sampstd / (2 * np.sqrt(npsrs)) *\
            np.sqrt(beta2 - (npsrs - 3.) / (npsrs - 1.))
        np.testing.assert_allclose(np.mean(periods), sampmean,
                                   rtol=0, atol=5 * sampmean_sterr,
                                   err_msg="x==E(P), y==np.mean(periods)")
        np.testing.assert_allclose(np.std(periods, ddof=1), sampstd,
                                   rtol=0, atol=5 * sampstd_sterr,
                                   err_msg="x==sP, y==np.std(periods, ddof=1)")
        
    @ptzd.parameterized.expand([(1000,),],
                               name_func=lambda fxn, n, par : "_{}".format(ptzd.parameterized.to_safe_name(str(par.args[0]))).join(fxn.__name__.split("_npsrs")))
    def test_pDistType__cc97_mean_std_ngen_npsrs(self, npsrs):
        """
        Cordes & Chernoff (1998) :math:`p(P) \propto P^{-\\alpha}, 1 < P < 20 ms`
        """
        random.seed(12345)
        N = npsrs
        alpha = 2.0
        a = 1.
        b = 20.
        pop = populate.generate(npsrs,
                                pDistType='cc97',
                                nostdout=True)
        periods = np.array([p.period for p in pop.population])
        m1 = powerlaw_moment(1, alpha, a, b)
        m2 = powerlaw_moment(2, alpha, a, b)
        m3 = powerlaw_moment(3, alpha, a, b)
        m4 = powerlaw_moment(4, alpha, a, b)
        mean = m1
        std = np.sqrt(m2 - m1 ** 2)
        mu4 = m4 - 4 * m1 * m3 + 6 * (m1 ** 2) * m2 - 3 * (m1 ** 4)
        beta2 = mu4 / (std ** 4)
        mean_sterr = std / np.sqrt(N)
        std_sterr = std / (2 * np.sqrt(N)) * np.sqrt(beta2 - (N - 3.) / (N - 1.))
        np.testing.assert_allclose(np.mean(periods), mean,
                                   rtol=0, atol=5 * mean_sterr,
                                   err_msg="x==E(P), y==np.mean(periods)")
        np.testing.assert_allclose(np.std(periods, ddof=1), std,
                                   rtol=0, atol=5 * std_sterr,
                                   err_msg="x==sP, y==np.std(periods, ddof=1)")

    def test_radialDistType_disk_have_zero_scale_height(self):
        random.seed(12345)
        npsrs = 1000
        pop = populate.generate(npsrs,
                                radialDistType='disk',
                                nostdout=True)
        zs = np.array([p.galCoords[2] for p in pop.population])
        np.testing.assert_array_equal(zs, np.zeros(len(zs)))

    @ptzd.parameterized.expand([("X", 0, -15., 15.),
                                ("Y", 1, -15., 15.),],
                               name_func=lambda fxn, n, par : "{}".format(ptzd.parameterized.to_safe_name(str(par.args[0]))).join(fxn.__name__.split("COORD")))
    def test_radialDistType_disk_galCOORD_mean_std(self,
                                                   coordname,
                                                   coordidx,
                                                   a, b):
        random.seed(12345)
        npsrs = 1000
        pop = populate.generate(npsrs,
                                radialDistType='disk',
                                nostdout=True)
        galcoord = np.array([p.galCoords[coordidx] for p in pop.population])
        mean = 0
        std = (b - a) / np.sqrt(12)
        mean_sterr = std / np.sqrt(npsrs)
        std_sterr = std * np.sqrt(9 / 5. - ((npsrs - 3) / (npsrs - 2))) /\
            (2 * np.sqrt(npsrs))
        np.testing.assert_allclose(
            np.mean(galcoord), mean,
            rtol=0, atol=5 * mean_sterr,
            err_msg="x==E(gal{0}), y==np.mean(gal{0})".format(coordname))
        np.testing.assert_allclose(
            np.std(galcoord, ddof=1), std,
            rtol=0, atol=5 * std_sterr,
            err_msg="x==sgal{0}, y==np.std(gal{0}, ddof=1)".format(coordname))

    @ptzd.parameterized.expand([("X", 0, -15., 15.),
                                ("Y", 1, -15., 15.),
                                ("Z", 2, -5., 5.),],
                               name_func=lambda fxn, n, par : "{}".format(ptzd.parameterized.to_safe_name(str(par.args[0]))).join(fxn.__name__.split("COORD")))
    def test_radialDistType_slab_galCOORD_mean_std(self,
                                                   coordname,
                                                   coordidx,
                                                   a, b):
        """
        Test sample mean and sample std for galXYZ consistent with uniform dist
        """
        random.seed(12345)
        npsrs = 1000
        pop = populate.generate(npsrs,
                                radialDistType='slab',
                                nostdout=True)
        galcoord = np.array([p.galCoords[coordidx] for p in pop.population])
        mean = 0
        std = (b - a) / np.sqrt(12)
        mean_sterr = std / np.sqrt(npsrs)
        momentfac = np.sqrt(9 / 5. - ((npsrs - 3) / (npsrs - 2))) /\
            (2 * np.sqrt(npsrs))
        std_sterr = std * momentfac
        np.testing.assert_allclose(
            np.mean(galcoord), mean,
            rtol=0, atol=5 * mean_sterr,
            err_msg="x==E(gal{0}), y==np.mean(gal{0})".format(coordname))
        np.testing.assert_allclose(
            np.std(galcoord, ddof=1), std,
            rtol=0, atol=5 * std_sterr,
            err_msg="x==sgal{0}, y==np.std(gal{0}, ddof=1)".format(coordname))
        
def powerlaw_moment(k, alpha, a, b):
    if alpha == 1.:
        c = 1. / np.log(b / a)
    else:
        c = (1. - alpha) / (b ** (1. - alpha) - a ** (1. - alpha))

    q = k + 1. - alpha

    if q == 0.:
        return c * np.log(b / a)
    else:
        return c * (b ** q - a ** q) / q
