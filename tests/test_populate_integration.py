"""
Integration tests for populate module
"""
import unittest
import pytest
import random
import numpy as np
try:
    from unittest.mock import patch
except ImportError:
    from mock import patch
from scipy import stats
from scipy.special import ndtri
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
        random.seed(12345)
        mu_ln = 1.5
        sig_ln = 0.58
        pop = populate.generate(npsrs,
                                pDistType='lnorm-e',
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

    @ptzd.parameterized.expand([("x",),
                                ("y",),
                                ("z",),],
                               name_func=lambda fxn, n, par : "{}".format(ptzd.parameterized.to_safe_name(str(par.args[0]))).join(fxn.__name__.split("COORD")))
    def test_velDistType_gauss_vCOORD_mean_std(self, coord):
        random.seed(12345)
        npsrs = 1000
        N = npsrs
        sigma = 196.3
        mean = 0
        pop = populate.generate(npsrs,
                                velDistType='gauss',
                                velDistPars=[0, sigma],
                                nostdout=True)
        v = np.array([getattr(p, "v" + coord) for p in pop.population])
        mean_sterr = sigma / np.sqrt(npsrs)
        mu4 = 3 * sigma ** 4.
        std_sterr = np.sqrt((mu4 - (sigma ** 4.) * (N - 3) / (N - 1)) / N) / \
            (2 * sigma)
        np.testing.assert_allclose(
            np.mean(v), mean,
            rtol=0, atol=5 * mean_sterr,
            err_msg="x==E(v{0}), y==np.mean(v{0})".format(coord))
        np.testing.assert_allclose(
            np.std(v, ddof=1), sigma,
            rtol=0, atol=5 * std_sterr,
            err_msg="x==s_v{0}, y==np.std(v{0}, ddof=1)".format(coord))

    def test_pdotDistType_lnorm_e_mean_std(self):
        """
        Test lognormal pdot sample mean and standard deviation match 
        pop.pdotmean and pop.pdotsigma to 5 x standard error
        """
        npsrs = 1000
        random.seed(12345)
        mu_ln = -33.75
        sig_ln = 2.59
        pop = populate.generate(npsrs,
                                pdotDistType='lnorm-e',
                                pdotDistPars=[mu_ln, sig_ln],
                                nostdout=True)
        pdots = np.array([p.pdot for p in pop.population])
        sampmean = np.exp(mu_ln + 0.5 * (sig_ln ** 2.))
        sampstd = np.sqrt((np.exp(sig_ln ** 2.) -1) * \
                          np.exp(2 * mu_ln + sig_ln ** 2.))
        sampmean_sterr = sampstd / np.sqrt(npsrs)
        beta2 = np.exp(4 * sig_ln**2) + 2 * np.exp(3 * sig_ln**2) +\
             3 * np.exp(2 * sig_ln**2) - 3
        sampstd_sterr = sampstd / (2 * np.sqrt(npsrs)) *\
            np.sqrt(beta2 - (npsrs - 3.) / (npsrs - 1.))
        np.testing.assert_allclose(np.mean(pdots), sampmean,
                                   rtol=0, atol=5 * sampmean_sterr,
                                   err_msg="x==E(Pdot), y==np.mean(pdots)")
        np.testing.assert_allclose(np.std(pdots, ddof=1), sampstd,
                                   rtol=0, atol=5 * sampstd_sterr,
                                   err_msg="x==sPdot, y==np.std(pdots, ddof=1)")

    def test_lumDistType_fk06_mean_std(self):
        """
        reconstruct lum_14000 scatter from fk06 params and check if
        delta_L ~ N(0, 0.8)
        """
        npsrs = 1000
        alpha = -1.5
        beta = 0.5
        gamma = 0.18

        pop = populate.generate(npsrs,
                                lumDistType='fk06',
                                lumDistPars=[alpha, beta, gamma],
                                nostdout=True)

        # reconstruct lum_1400 scatter from fk06 params
        delta_l = np.array([np.log10(p.lum_1400)
                            - np.log10(gamma)
                            - alpha * np.log10(p.period / 1000.)
                            - beta * np.log10(p.pdot * 1.e15)
                            for p in pop.population])

        n = len(delta_l)
        expected_mean = 0.0
        expected_std = 0.8

        mean_se = expected_std / np.sqrt(n)
        std_se = expected_std / np.sqrt(2.0 * (n - 1))

        # reconstructed residual should be < 5 * MC standard err
        self.assertLess(abs(np.mean(delta_l) - expected_mean),
                        5.0 * mean_se)

        self.assertLess(abs(np.std(delta_l, ddof=1) - expected_std),
                        5.0 * std_se)

        

class test_pdotDistType_logt(unittest.TestCase):

    def setUp(self):
        self.dvk_alpha = 1e-6
        self.dof = 2.67
        self.mu = -33.65
        self.sig = 1.48
    
    def test_generated_pdots_match_scipy_kstest(self):
        npsrs = 5000
        logpdot_dist = stats.t(df=self.dof, loc=self.mu, scale=self.sig)
        pop = populate.generate(npsrs,
                                pdotDistType='logt',
                                pdotDistPars=[self.dof, self.mu, self.sig],
                                nostdout=True)
            
        logpdots = np.log(np.array([p.pdot for p in pop.population]))
        
        # DVK inequality bound
        ks_limit = np.sqrt(np.log(2.0 / self.dvk_alpha) / (2.0 * len(logpdots)))
        ks_result = stats.kstest(logpdots, logpdot_dist.cdf)
        self.assertLess(ks_result.statistic, ks_limit)

    def test_pdot_population_attributes_set(self):
        npsrs = 10
        pop = populate.generate(npsrs,
                                pdotDistType='logt',
                                pdotDistPars=[self.dof, self.mu, self.sig],
                                nostdout=True)
        self.assertEqual(pop.pdotdof, self.dof)
        self.assertEqual(pop.pdotmu, self.mu)
        self.assertEqual(pop.pdotsigma, self.sig)        
        
class test_ppdotDistType_cohen27(unittest.TestCase):

    def setUp(self):
        self.dvk_alpha = 1e-6
        # expected distribution
        self.k_p = 2.142086
        self.theta_nu = 0.325553
        self.mu_logpdot = -46.344261
        self.sig_logpdot = 1.269051
        self.rho = 0.290428        
        self.draws = np.array([[self.k_p],   
                               [self.theta_nu],
                               [self.mu_logpdot],
                               [self.sig_logpdot],
                               [self.rho]])
        self.names = np.array(['k_P',
                               'theta_nu_ms',
                               'mu_logPdot',
                               'sigma_logPdot',
                               'rho'])
    
    def test_population_ppdotDistParams_match_expected(self):
        npsrs = 10
        pop = populate.generate(npsrs,
                                ppdotDistType='cohen27',
                                nostdout=True)
        expected = ["k_P", "theta_nu_ms", "mu_logPdot", "sigma_logPdot", "rho"]
        result = pop.ppdotDistParams.keys()
        np.testing.assert_array_equal(sorted(expected), sorted(result))

    def test_population_pDistType_is_None(self):
        npsrs = 10
        pop = populate.generate(npsrs,
                                ppdotDistType='cohen27',
                                nostdout=True)
        self.assertEqual(pop.pDistType, None)

    def test_period_invweibull_kstest(self):
        """
        Tests whether the CDF of npsrs period samples has a Kolmogorov-Smirnov 
        test statistic within the Dvoretsky-Kieffer-Wolfowitz inequality with 
        false failure rate given by dvk_alpha, 
        mocks ppdot_cohen27_posterior_draws, 
        can still fail by (small) random chance
        """
        npsrs = 5000
        p_dist = stats.invweibull(c=self.k_p, scale=1.0 / self.theta_nu)        
        
        with patch(
                'psrpoppy.populate.dists.ppdot_cohen27_posterior_draws',
                return_value=(self.draws, self.names)):
            pop = populate.generate(npsrs,
                                    ppdotDistType='cohen27',
                                    nostdout=True)
            
        periods = np.array([p.period for p in pop.population])
        
        # DVK inequality bound
        ks_limit = np.sqrt(np.log(2.0 / self.dvk_alpha) / (2.0 * len(periods)))
        ks_result = stats.kstest(periods, p_dist.cdf)
        self.assertLess(ks_result.statistic, ks_limit)

    def test_pdot_mean_std(self):
        npsrs = 5000
        
        with patch(
                'psrpoppy.populate.dists.ppdot_cohen27_posterior_draws',
                return_value=(self.draws, self.names)):
            pop = populate.generate(npsrs,
                                    ppdotDistType='cohen27',
                                    nostdout=True)
            
        pdots = np.array([p.pdot for p in pop.population])
        sampmean = np.exp(self.mu_logpdot + 0.5 * (self.sig_logpdot ** 2.))
        sampstd = np.sqrt((np.exp(self.sig_logpdot ** 2.) -1) * \
                          np.exp(2 * self.mu_logpdot + self.sig_logpdot ** 2.))
        sampmean_sterr = sampstd / np.sqrt(npsrs)
        beta2 = np.exp(4 * self.sig_logpdot**2) + 2 \
            * np.exp(3 * self.sig_logpdot**2) \
            + 3 * np.exp(2 * self.sig_logpdot**2) - 3
        sampstd_sterr = sampstd / (2 * np.sqrt(npsrs)) \
            * np.sqrt(beta2 - (npsrs - 3.) / (npsrs - 1.))
        np.testing.assert_allclose(np.mean(pdots), sampmean,
                                   rtol=0, atol=5 * sampmean_sterr,
                                   err_msg="x==E(Pdot), y==np.mean(pdots)")
        np.testing.assert_allclose(np.std(pdots, ddof=1), sampstd,
                                   rtol=0, atol=5 * sampstd_sterr,
                                   err_msg="x==sPdot, y==np.std(pdots, ddof=1)")
        
    def test_p_pdot_sample_correlation_matches_input(self):
        """
        test if correlation recovers expected Gaussian zscores
        """
        npsrs = 5000
        seed = 12345
        p_dist = stats.invweibull(c=self.k_p,
                                  scale=1.0 / self.theta_nu)
        pdot_dist = stats.lognorm(s=self.sig_logpdot,
                                  scale=np.exp(self.mu_logpdot))

        with patch(
                'psrpoppy.populate.dists.ppdot_cohen27_posterior_draws',
                return_value=(self.draws, self.names)):
            pop = populate.generate(npsrs,
                                    ppdotDistType='cohen27',
                                    nostdout=True)
        periods = np.array([p.period for p in pop.population])
        pdots = np.array([p.pdot for p in pop.population])        
        u = p_dist.cdf(periods) 
        v = pdot_dist.cdf(pdots)
        # Protect against numerical rounding to 0 or 1.
        eps = np.finfo(float).eps
        # recover normal z scores
        z_p = ndtri(np.clip(u, eps, 1.0 - eps))
        z_pdot = ndtri(np.clip(v, eps, 1.0 - eps))
        rho_recovered = np.corrcoef(z_p, z_pdot)[0, 1]

        # fisher transform to ~normal dist
        fisher_recovered = np.arctanh(rho_recovered)
        fisher_expected = np.arctanh(self.rho)
        fisher_se = 1.0 / np.sqrt(npsrs - 3.0)

        # difference must be within 5 fisher std of expected
        self.assertLess(abs(fisher_recovered - fisher_expected),
                        5.0 * fisher_se)

class test_populate_generate_error_handling(unittest.TestCase):
    
        @ptzd.parameterized.expand(
            [
                ("lumDistPars", [1,2,3,4], "lumDistType", "pow"),
                ("lumDistPars", [1,2,3,4], "lumDistType", "lnorm"),                
                ("velDistPars", [1,2,3,4], "velDistType", "gauss"),
                ("lumDistPars", [1,], "lumDistType", "pow"),
                ("lumDistPars", [1,], "lumDistType", "lnorm"),                
                ("velDistPars", [1,], "velDistType", "gauss"),
            ],
            name_func=lambda fxn, n, par: fxn.__name__.replace(
                "DISTVAR",
                str(par.args[2])
            ).replace(
                "DISTTYPE",
                str(par.args[3])
            ).replace(
                "DISTPARS",
                ptzd.parameterized.to_safe_name(str(par.args[0]))
            ).replace(
                "PARLEN",
                str(len(par.args[1]))
            )
        )
        def test_assert_DISTVAR_DISTTYPE_with_DISTPARS_len_PARLEN_raises_PopulateException(
                self,
                pars_varname,
                params,
                dist_varname,
                disttype
        ):
            with self.assertRaises(populate.PopulateException):
                 pop = populate.generate(10,
                                         nostdout=True,
                                         **{pars_varname: params,
                                            dist_varname: disttype})

        
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
