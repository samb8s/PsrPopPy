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
seedlib = C.CDLL(filepath+'/libgetseed.so')
yklib.ykr_.restype = C.c_float
yklib.llfr_.restype = C.c_float


# I'd like to implement this without resorting to globals. But for now, it should work
global FIRST_CALL
FIRST_CALL=True

class RadialModels:
    """ Class to store radial model methods"""
    def __init__(self):

        self._seed = None

    # what if I make seed a property, so we run it with same seed every time?
    
    def seed(self):
        #return C.c_int(seedlib.getseed_(C.byref(C.c_int(-1))))
        return C.c_int(5)
    """

    # get/set methods for seed
    @property
    def seed(self):
        return self._seed

    @seed.setter
    def seed(self):
        self._seed = C.c_int(seedlib.getseed_(C.byref(C.c_int(-1))))

    @seed.deleter
    def seed(self):
        del self._seed
    """
    # class methods
    
    def lfl06(self):
        """lfl06 model, using Y&K"""
        #print self.seed
        #return yklib.llfr_(C.byref(self.seed()))
        rand = C.c_float(random.random())
        return yklib.llfr_(C.byref(rand))

    def ykr(self):
        """ Y&K Model"""
        return yklib.ykr_(C.byref(self.seed()))


    def llfr(self):
        """Python implementation of the lfl06 model"""

        return self.ykr0(3.51, 7.89, 0.0)

    def ykr0(self, a, b, r1):
        """Python implementation of the Yusifov + Kucuk raidal model"""
        if FIRST_CALL:
            # set the global to false so we don't end up here again
            global FIRST_CALL
            FIRST_CALL = False

            self.amax = self.ykarea(500.0, 0.0, a, b, r1)

        # now get the result
        area = random.random() * self.amax
        return self.ykarea(500.0, area, a, b, r1)
            
    def ykarea(self, r, amax, a, b, r1):
        """Python implementation of Y+K"""

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
