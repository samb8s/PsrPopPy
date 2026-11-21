"""
Integration tests for evolve module
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
from psrpoppy import evolve

class test_lumDistType_fk06(unittest.TestCase):

    def setUp(self):
        self.lumdisttype = 'fk06'
        self.alpha = -1.5
        self.beta = 0.5
        self.gamma = 0.18
        self.lumdistpars = [self.alpha, self.beta, self.gamma]

    def test_population_parameter_attrs_set(self):
        pop = evolve.generate(2,
                              lumDistType=self.lumdisttype,
                              lumDistPars=self.lumdistpars,
                              nostdout=True)
        self.assertEqual(pop.lumPar1, self.alpha)
        self.assertEqual(pop.lumPar2, self.beta)
        self.assertEqual(pop.lumPar3, self.gamma)        

    def test_population_disttype_attr_set(self):
        pop = evolve.generate(2,
                              lumDistType=self.lumdisttype,
                              lumDistPars=self.lumdistpars,
                              nostdout=True)
        self.assertEqual(pop.lumDistType, self.lumdisttype)
        
    def test_lum_1400_mean_std(self):
        """
        reconstruct lum_14000 scatter from fk06 params and check if
        delta_L ~ N(0, 0.8)
        """
        npsrs = 1000

        pop = evolve.generate(npsrs,
                              lumDistType=self.lumdisttype,
                              lumDistPars=self.lumdistpars,
                              beamModel='none',
                              nodeathline=True,
                              nostdout=True)

        # reconstruct lum_1400 scatter from fk06 params
        delta_l = np.array([np.log10(p.lum_1400)
                            - np.log10(self.gamma)
                            - self.alpha * np.log10(p.period / 1000.)
                            - self.beta * np.log10(p.pdot * 1.e15)
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

class test_lumDistType_lnorm(unittest.TestCase):

    def setUp(self):
        self.lumdisttype = 'lnorm'
        self.lumdistpars = [-1.1, 0.9]

    def test_population_parameter_attrs_set(self):
        pop = evolve.generate(2,
                              lumDistType=self.lumdisttype,
                              lumDistPars=self.lumdistpars,
                              nostdout=True)
        self.assertEqual(pop.lumPar1, self.lumdistpars[0])
        self.assertEqual(pop.lumPar2, self.lumdistpars[1])        

    def test_population_disttype_attr_set(self):
        pop = evolve.generate(2,
                              lumDistType=self.lumdisttype,
                              lumDistPars=self.lumdistpars,
                              nostdout=True)
        self.assertEqual(pop.lumDistType, self.lumdisttype)        
        
    def test_lum_1400_mean_std(self):
        """
        Test base-10 lnorm period sample mean and standard deviation match 
        pop.pmean and pop.psigma to 5 x standard error
        """
        npsrs = 1000
        mu_log = self.lumdistpars[0]
        sig_log = self.lumdistpars[1]
        mu_ln = np.log(10) * mu_log
        sig_ln = np.log(10) * sig_log
        pop = evolve.generate(npsrs,
                              lumDistType=self.lumdisttype,
                              lumDistPars=self.lumdistpars,
                              beamModel='none',
                              nodeathline=True,
                              nostdout=True)
        lum1400 = np.array([p.lum_1400 for p in pop.population])

        sampmean = (10 ** mu_log) * np.exp(0.5 * sig_ln ** 2.)
        sampstd = sampmean * np.sqrt(np.exp(sig_ln ** 2.) - 1)
        sampmean_sterr = sampstd / np.sqrt(npsrs)
        beta2 = np.exp(4 * sig_ln**2) +\
            2 * np.exp(3 * sig_ln**2) +\
            3 * np.exp(2 * sig_ln**2) - 3
        sampstd_sterr = sampstd / (2 * np.sqrt(npsrs)) * \
            np.sqrt(beta2 - (npsrs - 3.) / (npsrs - 1.))
        
        np.testing.assert_allclose(np.mean(lum1400), sampmean,
                                   rtol=0, atol=5 * sampmean_sterr,
                                   err_msg="x==E(L), y==np.mean(lum_1400)")
        np.testing.assert_allclose(np.std(lum1400, ddof=1), sampstd,
                                   rtol=0, atol=5 * sampstd_sterr,
                                   err_msg="x==sL, y==np.std(lum_1400, ddof=1)")
        
    
class test_evolve_generate_error_handling(unittest.TestCase):
    
        @ptzd.parameterized.expand(
            [
                ("lumDistPars", [1,2,3,4], "lumDistType", "pow"),
                ("lumDistPars", [1,2,3,4], "lumDistType", "lnorm"),
                ("lumDistPars", [1,2,3,4], "lumDistType", "fk06"),
                ("lumDistPars", [1,], "lumDistType", "pow"),
                ("lumDistPars", [1,], "lumDistType", "lnorm"),
                ("lumDistPars", [1,], "lumDistType", "fk06"),
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
        def test_DISTVAR_DISTTYPE_with_DISTPARS_len_PARLEN_raises_EvolveException(
                self,
                pars_varname,
                params,
                dist_varname,
                disttype
        ):
            with self.assertRaises(evolve.EvolveException):
                 pop = evolve.generate(10,
                                       nostdout=True,
                                       **{pars_varname: params,
                                          dist_varname: disttype})

