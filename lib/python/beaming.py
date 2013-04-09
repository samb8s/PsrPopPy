#!/usr/bin/python

import math
import random

def tm98_fraction(pulsar):
    """
    Tauris and Manchester 1998 beaming fraction
    """
    periodterm = math.log10(pulsar.period / 1000.) - 1.0

    return 0.03 + 0.09 * periodterm**2.0

def wj08_fraction(pulsar):
    """ 
    Weltevrede & Johnston 2008 beaming model
    """
    if pulsar.chi < 1.0E-5:
        pulsar.chi = 0.

    rho = 5.4 * (pulsar.period/1000.)**(-0.5)
    beta_temp = math.degrees(math.acos(random.random()))
    if beta_temp <= rho:
        fraction = 1.
    else:
        fraction = 0.

    # get rho and chi in radians for simplicity
    chi_rad = math.radians(pulsar.chi)
    rho_rad = math.radians(rho)

    if pulsar.chi > rho and (pulsar.chi + rho)<90.:
        fraction = 2.0 * math.sin(chi_rad) * math.sin(rho_rad)
    elif pulsar.chi>rho and (pulsar.chi + rho) > 90.:
        fraction= math.cos(chi_rad - rho_rad)
    elif pulsar.chi <= rho and (pulsar.chi + rho) < 90.:
        fraction = 1.0 - math.cos(chi_rad + rho_rad)
    else:
        fraction = 1.


    return fraction
