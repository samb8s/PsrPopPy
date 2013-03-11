#!/usr/bin/python

import os

import math
import random

import ctypes as C


# get the c libraries
__dir__ = os.path.dirname(os.path.abspath(__file__))
filepath = os.path.join(__dir__, 'fortran')

ne2001lib = C.CDLL(os.path.join(filepath,'libne2001.so'))
ne2001lib.dm_.restype = C.c_float

tskylib = C.CDLL(os.path.join(filepath,'libtsky.so'))
tskylib.psr_tsky_.restype = C.c_float

slalib =  C.CDLL(os.path.join(filepath,'libsla.so'))
slalib.galtfeq_.restype= C.c_float

vxyzlib = C.CDLL(os.path.join(filepath, 'libvxyz.so'))
vxyzlib.vxyz_.restype = C.c_float

# Class definition
class GalacticOps:
    """
    Class for keeping all the galactic methods together
    
    """

    def __init__(self):
        pass
    
    def vxyz(self, pulsar):
        """Evolve a pulsar through galactic potential"""
        x,y,z = pulsar.galCoords
        vx,vy,vz = pulsar.vx, pulsar.vy, pulsar.vz
        age_Myr = pulsar.age/1.0E6

        x,y,z = C.c_float(x), C.c_float(y),C.c_float(z)
        vx,vy,vz = C.c_float(vx), C.c_float(vy),C.c_float(vz)
        age_Myr = C.c_float(age_Myr)
        bound = C.c_long(0)

        # run the evolution code
        vxyzlib.vxyz_(C.byref(C.c_float(0.005)),
                      C.byref(x),
                      C.byref(y),
                      C.byref(z),
                      C.byref(age_Myr),
                      C.byref(vx),
                      C.byref(vy),
                      C.byref(vz),
                      C.byref(x),
                      C.byref(y),
                      C.byref(z),
                      C.byref(vx),
                      C.byref(vy),
                      C.byref(vz),
                      C.byref(bound)
                      )

        # convert the output C types to python numbers
        pulsar.galCoords = x.value, y.value, z.value
        pulsar.vx = vx.value
        pulsar.vy = vy.value
        pulsar.vz = vz.value



    def calc_dtrue(self, (x, y, z)):
        """Calculate true distance to pulsar from the sun."""
        rsun = 8.5 # kpc
        return math.sqrt( x*x + (y-rsun)*(y-rsun) + z*z)

    def calcXY(self, r0):
        """Calculate the X, Y, Z alactic coords for the pulsar."""
        # calculate a random theta in a unifrom distribution
        theta = 2.0 * math.pi * random.random()

        # calc x and y
        x = r0 * math.cos(theta)
        y = r0 * math.sin(theta)

        return x, y

    def ne2001_dist_to_dm(self, dist, gl, gb):
        """Use NE2001 distance model."""
        dist = C.c_float(dist)
        gl = C.c_float(gl)
        gb = C.c_float(gb)
        return ne2001lib.dm_(C.byref(dist),
                             C.byref(gl),
                             C.byref(gb),
                             C.byref(C.c_int(4)),
                             C.byref(C.c_float(0.0)))

    def lmt85_dist_to_dm(self, dist, gl, gb):
        """ Use Lyne, Manchester & Taylor distance model to compute DM."""
        dist = C.c_float(dist)
        gl = C.c_float(gl)
        gb = C.c_float(gb)
        return ne2001lib.dm_(C.byref(dist),
                             C.byref(gl),
                             C.byref(gb),
                             C.byref(C.c_int(0)),
                             C.byref(C.c_float(0.0)))

    def lb_to_radec(self, l, b):
        """Convert l, b to RA, Dec using SLA fortran (should be faster)."""
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

    def radec_to_lb(self, ra, dec):
        """Convert RA, Dec to l, b using SLA fortran.
        Be sure to return l in range -180 -> +180"""
        l = C.c_float(0.)
        b = C.c_float(0.)
        ra = C.c_float(ra)
        dec = C.c_float(dec)
        # call with arg = -1 to convert in reverse!
        slalib.galtfeq_(C.byref(l), 
                        C.byref(b),
                        C.byref(ra),
                        C.byref(dec),
                        C.byref(C.c_int(-1)))
        if l.value>180.:
            l.value -= 360.
        return l.value, b.value

    def tsky(self, gl, gb, freq):
        """ Calculate Galactic sky temperature for a given observing frequency (MHz), 
            l and b."""
        gl =  C.c_float(gl)
        gb = C.c_float(gb)
        freq = C.c_float(freq)
        return tskylib.psr_tsky_(C.byref(gl),
                                 C.byref(gb),
                                 C.byref(freq))


    def xyz_to_lb(self, (x, y, z)):
        """ Convert galactic xyz in kpc to l and b in degrees."""
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

    def lb_to_xyz(self, gl, gb, dist):
        """ Convert galactic coords to Galactic XYZ."""
        rsun = 8.5 # kpc

        l = math.radians(gl)
        b = math.radians(gb)

        x = dist * math.cos(b) * math.sin(l)
        y = rsun - dist * math.cos(b) * math.cos(l)
        z = dist * math.sin(b)

        return (x, y, z)

    def scatter_bhat(self, dm, scatterindex, freq_mhz):
        """Calculate bhat et al 2004 scattering timescale for freq in MHz."""
        logtau = -6.46 + 0.154*math.log10(dm) + \
                    1.07*math.log10(dm)*math.log10(dm) + \
                    scatterindex*math.log10(freq_mhz/1000.0)
        
        # seems like scattering timescale is tooooooo big? Definitely my weffs are too big
        # return tau with power scattered with a gaussian, width 0.8
        return math.pow(10.0, random.gauss(logtau, 0.8)) 

    def _glgboffset(self, gl1, gb1, gl2, gb2):
        """Calculate the angular distance (deg) between two
        points in galactic coordinates"""
        # Angular offset in polar coordinates
        # taken brazenly from
        #http://www.atnf.csiro.au/people/Tobias.Westmeier/tools_hihelpers.php

        # requires gb conversion from +90 -> -90 to 0 -> 180
        gb1 = 90.0 - gb1
        gb2 = 90.0 - gb2

        term1 = math.cos(math.radians(gb1)) * math.cos(math.radians(gb2))
        term2 = math.sin(math.radians(gb1)) * math.sin(math.radians(gb2)) 
        term3 = math.cos(math.radians(gl1) - math.radians(gl2))
        cosalpha = term1 + term2*term3

        return math.degrees(math.acos(cosalpha))
