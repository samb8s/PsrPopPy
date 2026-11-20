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

    def test_lum_1400_mean_std(self):
        """
        reconstruct lum_14000 scatter from fk06 params and check if
        delta_L ~ N(0, 0.8)
        """
        npsrs = 1000
        alpha = -1.5
        beta = 0.5
        gamma = 0.18

        pop = evolve.generate(npsrs,
                              lumDistType='fk06',
                              lumDistPars=[alpha, beta, gamma],
                              beamModel='none',
                              nodeathline=True,
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
