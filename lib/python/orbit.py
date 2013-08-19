#!/usr/bin/python

import math
import random

class OrbitException(Exception):
    pass

class Orbit(object):
    """ Class to store an individual pulsar"""
    def __init__(self,
                 orbital_period_days=None,
                 pulsar_mass_msolar = None,
                 companion_mass_msolar=None,
                 long_peri_degrees=None,
                 inclination_degrees=None,
                 eccentricity=None):
        """___init___ function for the Orbit class"""

        self.orbital_period_days = orbital_period_days

        self.pulsar_mass_msolar = pulsar_mass_msolar
        self.companion_mass_msolar = companion_mass_msolar
        self.long_peri_degrees = long_peri_degrees
        self.inclination_degrees = inclination_degrees
        self.ecc = eccentricity
