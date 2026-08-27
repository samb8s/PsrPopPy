import sys
import os

import math
import unittest
import numpy as np
import psrpoppy.galacticops as go

class test_double_sided_exp(unittest.TestCase):

    def test_double_sided_exp_zero(self):
        result = go._double_sided_exp(0.0, origin=0.0)
        self.assertEqual(result, 0.0)

    def test_double_sided_exp_origin(self):
        result = go._double_sided_exp(0.0, origin=10.0)
        self.assertEqual(result, 10.0)

    def test_double_sided_exp_positive_sign(self):
        """
        Mock random.random and choice.choice to make the test deterministic
        """
        old_random = go.random.random
        old_choice = go.random.choice
        try:
            go.random.random = lambda: math.exp(-1.0)
            go.random.choice = lambda x: 1.0

            result = go._double_sided_exp(10.0, origin=0.0)
            self.assertAlmostEqual(result, -10.0)
        finally:
            go.random.random = old_random
            go.random.choice = old_choice

    def test_double_sided_exp_negative_sign(self):
        """
        Mock random.random and choice.choice to make the test deterministic
        """
        old_random = go.random.random
        old_choice = go.random.choice
        try:
            go.random.random = lambda: math.exp(-1.0)
            go.random.choice = lambda x: -1.0

            result = go._double_sided_exp(10.0, origin=0.0)
            self.assertAlmostEqual(result, 10.0)
        finally:
            go.random.random = old_random
            go.random.choice = old_choice

class test_vcirc(unittest.TestCase):
    """
    Test implementation of Kuijken & Gilmore (1989) galactic potential
    for computing circular velocity 
    """

    def test_vcirc_zero_at_GC(self):
        self.assertAlmostEqual(go.vcirc(0), 0.0)

    def test_vcirc_at_R0(self):
        """See K&J (89) Table 2"""
        np.testing.assert_allclose(go.vcirc(7.8), 222., rtol=0.01)

class test_scale_bhat(unittest.TestCase):

    def test_default_tscatter_increases_as_frequency_increases(self):
        tscatter_1400 = go.scale_bhat(1., 1400.)
        tscatter_2000 = go.scale_bhat(1., 2000.)
        np.testing.assert_array_less(tscatter_2000, tscatter_1400)
