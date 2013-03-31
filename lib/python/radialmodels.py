#!/usr/bin/python

import os
import sys
import random
import math

import ctypes as C

# get the FORTRAN libraries
__dir__ = os.path.dirname(os.path.abspath(__file__))
__libdir__ = os.path.dirname(__dir__)
filepath = os.path.join(__libdir__, 'fortran')

yklib = C.CDLL(filepath+'/libykarea.so')
yklib.ykr_.restype = C.c_float
yklib.llfr_.restype = C.c_float

seedlib = C.CDLL(filepath+'/libgetseed.so')
seedlib.getseed_.restype = C.c_int

def seed():
    return C.c_int(seedlib.getseed_(C.byref(C.c_int(-1))))

# class methods

def slabdist():
    x = -15.0 + random.random()*30.0
    y = -15.0 + random.random()*30.0
    z = -5.0 + random.random() * 10.0 

    return (x, y, z)

def diskdist():
    x = -15.0 + random.random()*30.0
    y = -15.0 + random.random()*30.0
    return (x, y, 0.0)

def lfl06():
    """lfl06 model, using Y&K"""
    #print self.seed
    return yklib.llfr_(C.byref(seed()))

def ykr():
    """ Y&K Model"""
    return yklib.ykr_(C.byref(seed()))

def spiralize( r):
    """ Make spiral arms, as seen in Fuacher-Giguere & Kaspi 2006"""

    # definitions
    k_list = [4.25,4.25,4.89,4.89]
    r0_list = [3.48,3.48,4.9,4.9]
    theta0_list = [1.57,4.71,4.09,0.95]
    
    # select a spiral arm ( 1->4)
    arm = random.choice([0,1,2,3]) 
    k = k_list[arm]
    r0 = r0_list[arm]
    theta0 = theta0_list[arm]

    # pick an angle
    theta = k * math.log(r/r0) + theta0

    # blurring angle
    angle = 2*math.pi * random.random() * math.exp(-0.35* r)

    if random.random()<0.5:
        angle = 0 - angle

    #modify theta
    theta += angle

    # blur in radial direction a little
    dr = math.fabs(random.gauss(0.0, 0.5 * r))
    angle = random.random() * 2.0 * math.pi
    dx=dr * math.cos(angle)
    dy = dr * math.sin(angle)

    x = r * math.cos(theta) + dx
    y = r * math.cos(theta) + dy

    return x, y

def _double_sided_exp( scale, origin=0.0):
    """Exponential distribution around origin, with scale height scale."""
    if scale == 0.0:
        return origin

    rn = random.random()
    sign = random.choice([-1.0, 1.0])

    return origin + sign * scale * math.log(rn)
