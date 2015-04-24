#!/usr/bin/python

import os

import ctypes as C

# get the FORTRAN libraries
__dir__ = os.path.dirname(os.path.abspath(__file__))
__libdir__ = os.path.dirname(__dir__)
fortranpath = os.path.join(__libdir__, 'fortran')

gammalib = C.CDLL(os.path.join(fortranpath, 'libgamma.so'))
gammalib.gamma1_.restype = C.c_double
gammalib.gamma2_.restype = C.c_double
gammalib.gamma3_.restype = C.c_double


# BEGIN FUNCTION DEFINITIONS
def gamma1(pulsar, tobs, harmonic):

    """S/N ratio decrease in a standard pulsar search"""
    gamma = C.c_double(0.)

    harm = float(harmonic)
    m1 = pulsar.pulsar_mass_msolar
    m2 = pulsar.companion_mass_msolar
    period_s = pulsar.period / 1000.
    long_peri = pulsar.long_peri_degrees
    inc_deg = pulsar.inclination_degrees
    ecc = pulsar.ecc
    orb_p = pulsar.orbital_period_days

    gammalib.gamma1_(
                     C.byref(C.c_double(harm)),
                     C.byref(C.c_double(tobs)),
                     C.byref(C.c_double(m1)),
                     C.byref(C.c_double(m2)),
                     C.byref(C.c_double(period_s)),
                     C.byref(C.c_double(long_peri)),
                     C.byref(C.c_double(inc_deg)),
                     C.byref(C.c_double(ecc)),
                     C.byref(C.c_double(orb_p)),
                     C.byref(gamma)
                     )

    return gamma.value


def gamma2(pulsar, tobs, harmonic):
    """S/N ratio decrease in an accel pulsar search"""

    gamma = C.c_double(0.)

    harm = float(harmonic)
    m1 = pulsar.pulsar_mass_msolar
    m2 = pulsar.companion_mass_msolar
    period_s = pulsar.period / 1000.
    long_peri = pulsar.long_peri_degrees
    inc_deg = pulsar.inclination_degrees
    ecc = pulsar.ecc
    orb_p = pulsar.orbital_period_days

    gammalib.gamma2_(
                     C.byref(C.c_double(harm)),
                     C.byref(C.c_double(tobs)),
                     C.byref(C.c_double(m1)),
                     C.byref(C.c_double(m2)),
                     C.byref(C.c_double(period_s)),
                     C.byref(C.c_double(long_peri)),
                     C.byref(C.c_double(inc_deg)),
                     C.byref(C.c_double(ecc)),
                     C.byref(C.c_double(orb_p)),
                     C.byref(gamma)
                     )

    return gamma.value


def gamma3(pulsar, tobs, harmonic):
    """S/N ratio decrease in an accel & jerk pulsar search"""

    gamma = C.c_double(0.)

    harm = float(harmonic)
    m1 = pulsar.pulsar_mass_msolar
    m2 = pulsar.companion_mass_msolar
    period_s = pulsar.period / 1000.
    long_peri = pulsar.long_peri_degrees
    inc_deg = pulsar.inclination_degrees
    ecc = pulsar.ecc
    orb_p = pulsar.orbital_period_days

    gammalib.gamma3_(
                     C.byref(C.c_double(harm)),
                     C.byref(C.c_double(tobs)),
                     C.byref(C.c_double(m1)),
                     C.byref(C.c_double(m2)),
                     C.byref(C.c_double(period_s)),
                     C.byref(C.c_double(long_peri)),
                     C.byref(C.c_double(inc_deg)),
                     C.byref(C.c_double(ecc)),
                     C.byref(C.c_double(orb_p)),
                     C.byref(gamma)
                     )

    return gamma.value
