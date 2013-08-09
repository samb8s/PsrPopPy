#!/usr/bin/python

import math
import random

class OrbitException(Exception):
    pass

class Orbit(object):
    """ Class to store an individual pulsar"""
    def __init__(self,
                 orbital_period_seconds=None,
                 pulsar_mass_msolar = None,
                 companion_mass_msolar=None,
                 long_peri_degrees=None,
                 inclination_degrees=None,
                 eccentricity=None):
        """___init___ function for the Orbit class"""

        self.orbital_period_seconds = orbital_period_seconds

        self.pulsar_mass_msolar = pulsar_mass_msolar
        self.companion_mass_msolar = companion_mass_msolar
        self.long_peri_degrees = long_peri_degrees
        self.inclination_degrees = inclination_degrees
        self.ecc = eccentricity

        # get/set properties
        self.__orbital_period_days = None
        self.__omega_orbit = None

    @property
    def orbital_period_days(self):
        return self.__orbital_period_days

    @orbital_period_days.setter
    def orbital_period_days(self, orbital_period_days):
        return self.orbital_period_seconds / 86400.

    @property
    def omega_orbit(self):
        return self.__omega_orbit

    @omega_orbit.setter
    def omega_orbit(self, omega_orbit):
        return math.pi * 2.0 / self.orbital_period_seconds
