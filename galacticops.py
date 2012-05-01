#!/usr/bin/python

import os

import math
import random

import ctypes as C


# get the c libraries
__dir__ = os.path.dirname(os.path.abspath(__file__))
filepath = os.path.join(__dir__, 'fortran')

ne2001lib = C.CDLL(filepath+'/libne2001.so')
ne2001lib.dm_.restype = C.c_float
tskylib = C.CDLL(filepath+'/libtsky.so')
tskylib.psr_tsky_.restype = C.c_float
slalib =  C.CDLL(filepath+'/libsla.so')



# Class definition
class GalacticOps:
    """Class for keeping all the galactic methods together"""
    def __init__(self):
        pass
    
    def calc_dtrue(self, (x, y, z)):
        """Calculate true distance to pulsar from the sun"""
        rsun = 8.5 # kpc
        return math.sqrt( x*x + (y-rsun)*(y-rsun) + z*z)

    def calcXY(self, r0):
        """Calculate the X, Y, Z alactic coords for the pulsar"""
        # calculate a random theta in a unifrom distribution
        theta = 2.0 * math.pi * random.random()

        # calc x and y
        x = r0 * math.cos(theta)
        y = r0 * math.sin(theta)

        # calculate z - for now use exponential distribution of expenential tail zscale
        return x, y

    def ne2001_dist_to_dm(self, dist, gl, gb):
        """Use NE2001 distance model"""
        dist = C.c_float(dist)
        gl = C.c_float(gl)
        gb = C.c_float(gb)
        return ne2001lib.dm_(C.byref(dist),
                             C.byref(gl),
                             C.byref(gb),
                             C.byref(C.c_int(4)),
                             C.byref(C.c_float(0.0)))

    def lm98_dist_to_dm(self, dist, gl, gb):
        """ Use Lyne & Manchester distance model to compute DM"""
        dist = C.c_float(dist)
        gl = C.c_float(gl)
        gb = C.c_float(gb)
        return ne2001lib.dm_(C.byref(dist),
                             C.byref(gl),
                             C.byref(gb),
                             C.byref(C.c_int(0)),
                             C.byref(C.c_float(0.0)))

    def lb_to_radec(self, l, b):
        """Convert l, b to RA, Dec using SLA fortran (should be faster)"""
        ra = C.c_float(0.)
        dec = C.c_float(0.)
        l = C.c_float(l)
        b = C.c_float(b)
        # call with final arg 1 to do conversion in right direction!
        slalib.galtfeq_(C.byref(l), 
                        C.byref(b),
                        C.byref(ra),
                        C.byref(dec),
                        C.byref(C.c_int(1)))
        return ra.value, dec.value

    def tsky(self, gl, gb, freq):
        """ Calculate Galactic sky temperature for a given observing frequency (MHz), 
            l and b"""
        gl =  C.c_float(gl)
        gb = C.c_float(gb)
        freq = C.c_float(freq)
        return tskylib.psr_tsky_(C.byref(gl),
                                 C.byref(gb),
                                 C.byref(freq))
        

    def xyz_to_lb(self, (x, y, z)):
        """ Convert galactic xyz in kpc to l and b in degrees"""
        rsun = 8.5 # kpc

        # distance to pulsar
        d = math.sqrt( x*x + (rsun-y)*(rsun-y) + z*z)
        # radial distance
        b = math.asin(z/d)

        # take cosine
        dcb = d * math.cos(b)

        if y<=rsun:
            if math.fabs(x/dcb) > 1.0:
                l = 1.57079632679
            else:
                l = math.asin(x/dcb)
        else:
            if math.fabs(x/dcb) > 1.0:
                l = 0.0
            else:
                l = math.acos(x/dcb)

            l += 1.57079632679
            if x < 0.:
                l -= 6.28318530718

        # convert back to degrees
        l = math.degrees(l)
        b = math.degrees(b)

        # convert to -180 < l < 180
        if l > 180.0:
            l -= 360.0

        return l, b


    def scatter_bhat(self, dm, scatterindex, freq_mhz):
        """Calculate bhat et al 2004 scattering timescale for freq in MHz"""
        logtau = -6.46 + 0.154*math.log10(dm) + \
                    1.07*math.log10(dm)*math.log10(dm) + \
                    scatterindex*math.log10(freq_mhz/1000.0)
        
        # seems like scattering timescale is tooooooo big? Definitely my weffs are too big
        # return tau with power scattered with a gaussian, width 0.8
        return math.pow(10.0, random.gauss(logtau, 0.8)) 
