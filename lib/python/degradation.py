#!/usr/bin/python

import os

import ctypes as C

# get the FORTRAN libraries
__dir__ = os.path.dirname(os.path.abspath(__file__))
__libdir__ = os.path.dirname(__dir__)
fortranpath = os.path.join(__libdir__, 'fortran')

gammalib = C.CDLL(os.path.join(fortranpath,'libgamma.so'))
gammalib.gamma1_.restype = C.c_double
gammalib.gamma2_.restype = C.c_double
gammalib.gamma3_.restype = C.c_double


# BEGIN FUNCTION DEFINITIONS
def gamma1(pulsar):

    gamma = C.c_double(0.)
    pulsar.pulsar_mass_msolar

    gammalib.gamma1_(
                     C.byref(C.c_double(4.)),
                     C.byref(C.c_double(1000.)), 
                     C.byref(C.c_double(1.4)), 
                     C.byref(C.c_double(0.3)), 
                     C.byref(C.c_double(0.08)), 
                     C.byref(C.c_double(30.)), 
                     C.byref(C.c_double(60.)), 
                     C.byref(C.c_double(0.1)), 
                     C.byref(C.c_double(0.1)), 
                     C.byref(gamma)
                     )

    return gamma.value

def gamma2(pulsar):

    gamma = C.c_double(0.)
    pulsar.pulsar_mass_msolar

    gammalib.gamma2_(
                     C.byref(C.c_double(4.)),
                     C.byref(C.c_double(1000.)), 
                     C.byref(C.c_double(1.4)), 
                     C.byref(C.c_double(0.3)), 
                     C.byref(C.c_double(0.08)), 
                     C.byref(C.c_double(30.)), 
                     C.byref(C.c_double(60.)), 
                     C.byref(C.c_double(0.1)), 
                     C.byref(C.c_double(0.1)), 
                     C.byref(gamma)
                     )

    return gamma.value

def gamma3(pulsar):

    gamma = C.c_double(0.)
    pulsar.pulsar_mass_msolar

    gammalib.gamma3_(
                     C.byref(C.c_double(4.)),
                     C.byref(C.c_double(1000.)), 
                     C.byref(C.c_double(1.4)), 
                     C.byref(C.c_double(0.3)), 
                     C.byref(C.c_double(0.08)), 
                     C.byref(C.c_double(30.)), 
                     C.byref(C.c_double(60.)), 
                     C.byref(C.c_double(0.1)), 
                     C.byref(C.c_double(0.1)), 
                     C.byref(gamma)
                     )

    return gamma.value
