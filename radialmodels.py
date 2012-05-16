#!/usr/bin/python

import os
import sys
import random
import math

import ctypes as C

# get the c libraries
__dir__ = os.path.dirname(os.path.abspath(__file__))
filepath = os.path.join(__dir__, 'fortran')
yklib = C.CDLL(filepath+'/libykarea.so')
yklib.ykr_.restype = C.c_float
yklib.llfr_.restype = C.c_float

seedlib = C.CDLL(filepath+'/libgetseed.so')
seedlib.getseed_.restype = C.c_int

# I'd like to implement this without resorting to globals. But for now, it should work
global FIRST_CALL
FIRST_CALL=True

class RadialModels:
    """ Class to store radial model methods."""
    def __init__(self):
        """Initialize the class."""
        pass

    def seed(self):
        return C.c_int(seedlib.getseed_(C.byref(C.c_int(-1))))
        #return C.c_int(5)
    
    # class methods

    def slabdist(self):
        x = -15.0 + random.random()*30.0
        y = -15.0 + random.random()*30.0
        z = -5.0 + random.random() * 10.0 

        return (x, y, z)

    def diskdist(self):
        x = -15.0 + random.random()*30.0
        y = -15.0 + random.random()*30.0
        return (x, y, 0.0)
    
    def lfl06(self):
        """lfl06 model, using Y&K"""
        #print self.seed
        return yklib.llfr_(C.byref(self.seed()))
        #rand = C.c_float(random.random())
        #return yklib.llfr_(C.byref(rand))

    def ykr(self):
        """ Y&K Model"""
        #rand = C.c_float(random.random())
        return yklib.ykr_(C.byref(self.seed()))


    def llfr(self):
        """Python implementation of the lfl06 model. Not in use."""

        return self.ykr0(3.51, 7.89, 0.0)

    def ykr0(self, a, b, r1):
        """Python implementation of the Y&K model. Not in use."""
        if FIRST_CALL:
            # set the global to false so we don't end up here again
            global FIRST_CALL
            FIRST_CALL = False

            self.amax = self.ykarea(500.0, 0.0, a, b, r1)

        # now get the result
        area = random.random() * self.amax
        return self.ykarea(500.0, area, a, b, r1)
            
    def ykarea(self, r, amax, a, b, r1):
        """Python implementation of Y+K, not in use."""

        # two "constants"
        DX = 0.01
        RSUN = 8.5 #kpc

        # do the code
        r0 = RSUN + r1
        ykarea = 0.0

        # create list to loop over
        nloops = int(math.ceil(r/DX))
        looplist = [x * DX for x in range(nloops)]

        for x in looplist:
            ykarea = ykarea + (((x+r1)/r0)**a) * math.exp(-1.0*b*(x-RSUN)/r0)*x*DX
            if amax > 0.0 and ykarea > amax:
                return x

        return ykarea
