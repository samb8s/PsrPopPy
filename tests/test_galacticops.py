import sys
import os

import math
import unittest

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
