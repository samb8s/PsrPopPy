"""
Integration tests NE2001
"""
import unittest
import pytest
import random
import numpy as np
import parameterized as ptzd
from psrpoppy.galacticops import ne2001_dist_to_dm
import psrpoppy

class test_NE2001_DM_clumps(unittest.TestCase):
    """
    DM should never be negative, especially around DM 'clumps' it should increase
    Remember to rebuild fortran SOs after checking out a new branch
    """
    
    def test_positive_dm_around_NGC6334N(self):
        lng = np.linspace(-179.9, 179.9, 1000)
        lat = np.linspace(-89.9, 89.9, 1000)
        dtrue = 1.63555909561 #kpc
        gl = -8.82877643574 # deg
        gb = 0.53676598385 # deg
        dm_bslice = np.array([ne2001_dist_to_dm(dtrue, l, gb) for l in lng])
        dm_lslice = np.array([ne2001_dist_to_dm(dtrue, gl, b) for b in lat])
        stack = np.concatenate([dm_bslice, dm_lslice])
        np.testing.assert_array_less(np.zeros(len(stack)), stack)
