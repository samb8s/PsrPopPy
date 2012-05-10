#!/usr/bin/python

class Pulsar:
    """ Class to store pulsar information"""
    def __init__(self=None,
                 period=None,
                 dm=None,
                 gl=None,
                 gb=None,
                 galCoords=None,
                 dtrue=None,
                 lum1400=None,
                 spindex=None,
                 scindex=None,
                 alpha=None,
                 rho=None,
                 width_degree=None):
        """___init___ function for the Pulsar class"""
        self.period = period
        self.dm = dm
        
        self.gl = gl
        self.gb = gb
        self.galCoords = galCoords
        self.dtrue = dtrue

        self.lum1400 = lum1400
        self.spindex = spindex
        self.scindex = scindex
        self.alpha = alpha
        self.rho = rho
        self.width_degree = width_degree

