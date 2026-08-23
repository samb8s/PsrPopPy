"""
Unit and integration tests for pulsar module
"""
import unittest
import numpy as np
from mock import patch
from psrpoppy.pulsar import Pulsar
from psrpoppy import pulsar

class test_proper_motion(unittest.TestCase):
    """
    test integration of proper motion property of Pulsar class
    (includes galacticops)
    """

    def test_transverse_velocity_2xVsun_on_opposite_side_GC(self):
        """cancel out prop motion constant by setting dtrue = 0.211"""
        v_expected = 2 * pulsar.VSUN_CIRC
        testpsr = Pulsar(galCoords=(0, -8.5, 0), dtrue=0.211) # opposite side GC
        testpsr.vx = 0.
        testpsr.vy = 0.
        testpsr.vz = 0.
        v_calc = testpsr.pm
        self.assertAlmostEqual(v_expected, v_calc)

    def test_pulsar_close_to_GC_has_higher_transverse_velocity(self):
        """cancel out prop motion constant by setting dtrue = 0.211"""
        testpsr_far = Pulsar(galCoords=(0, 6., 0), dtrue=0.211)
        testpsr_near = Pulsar(galCoords=(0, 0.5, 0), dtrue=0.211)
        testpsr_far.vx = 0.
        testpsr_far.vy = 0.
        testpsr_far.vz = 0.
        testpsr_near.vx = 0.
        testpsr_near.vy = 0.
        testpsr_near.vz = 0.
        vfar = testpsr_far.pm
        vnear = testpsr_near.pm
        np.testing.assert_array_less(vfar, vnear)

    def test_pulsar_with_zero_transverse_velocity_has_zero_proper_motion(self):
        """true unit test, monkeypatch pulsar.go.vcirc and pulsar.VSUN_CIRC"""
        def patched_vcirc(r):
            if r == 10.0:
                return 200.0
            elif r == 8.5:
                return 220.0
        with patch('psrpoppy.pulsar.go.vcirc', side_effect=patched_vcirc):
            with patch('psrpoppy.pulsar.VSUN_CIRC', 220.0):
                testpsr = Pulsar(galCoords=(0, 10., 0), dtrue=1.5)
                testpsr.vx = 20.
                testpsr.vy = 100.
                testpsr.vz = 0.
                self.assertAlmostEqual(testpsr.pm, 0.)
