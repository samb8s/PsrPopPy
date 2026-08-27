#!/usr/bin/python

import math
import random

from orbit import Orbit
import galacticops as go

RSUN = 8.5 #kpc
VSUN_CIRC = go.vcirc(RSUN)

class PulsarException(Exception):
    pass


class Pulsar(Orbit):
    """ Class to store an individual pulsar"""
    def __init__(self,
                 period=None,
                 pdot=None,
                 dm=None,
                 gl=None,
                 gb=None,
                 galCoords=None,
                 r0=None,
                 dtrue=None,
                 lum_1400=None,
                 spindex=None,
                 alpha=None,
                 rho=None,
                 width_degree=None,
                 snr=None,
                 beaming=None,
                 scindex=-3.86,
                 t_scatter=None,
                 gpsFlag=0,
                 gpsA=None,
                 brokenFlag=0,
                 brokenSI=None,
                 *args,
                 **kwargs):
        """___init___ function for the Pulsar class"""

        # initialise the inherited orbit class
        super(Pulsar, self).__init__(*args, **kwargs)

        self.period = period
        self.pdot = pdot
        self.dm = dm

        # convert to -180->+180 range
        if gl > 180.:
            gl -= 360.

        self.gl = gl
        self.gb = gb
        self.galCoords = galCoords
        self.r0 = r0
        self.dtrue = dtrue
        self.vx = None
        self.vy = None
        self.vz = None
        
        self.lum_1400 = lum_1400
        self.spindex = spindex
        self.scindex = scindex

        # set the scattering timescale
        # for 1.4 GHz (will scale to obs freq)
        self.t_scatter = None

        self.alpha = alpha
        self.rho = rho
        self.width_degree = width_degree

        self.beaming = beaming

        self.gpsFlag = gpsFlag
        self.gpsA = gpsA

        self.brokenFlag = brokenFlag
        self.brokenSI = brokenSI

        self.snr = snr

        # add this little flag which can be
        # switched if the pulsar is detected in any survey
        self.detected = False

        # need to add pulsar dead/alive for evolution code
        self.dead = False

    # methods to calculate derived properties
    def s_1400(self):
        """Calculate the flux of the pulsar"""
        if self.lum_1400 is None:
            raise PulsarException(
                   'Luminosity is not defined')
        elif self.dtrue is None:
            raise PulsarException(
                   'Distance not defined')

        return self.lum_1400 / self.dtrue / self.dtrue

    def efficiency(self):
        """Calculate pulsar efficiency at 1400 MHz, L1400 / Edot"""

        edot = self.edot()
        if edot is None:
            return None
        else:
            # return self.lum_1400 * 1.0e-26 * 3.086e21 * 3.086E21 / edot
            return self.lum_1400 * 7.4E27 / edot

    def edot(self):
        """Return the Edot of the pulsar in erg / s"""

        if self.pdot is None:
            return None

        pdot_15 = self.pdot * 1.0E15
        return 3.95E31 * pdot_15 / (self.period/1000.)**3

    @property
    def pm(self):
        """
        Compute proper motion in mas/yr
        """
        required_attr = ['galCoords', 'vx', 'vy', 'vz', 'dtrue']
        for attr in required_attr:
            if getattr(self, attr) is None:
                raise PulsarException(
                    '{} is not defined.'.format(attr))
            
        x, y, z = self.galCoords
        r = math.sqrt(x ** 2 + y ** 2) # cyclindrical radius
        t = math.atan2(y, x)

        v_r = go.vcirc(r)
        
        # velocity components
        vtot_x = self.vx + v_r * math.sin(t) - VSUN_CIRC
        vtot_y = self.vy - v_r * math.cos(t)
        vtot_z = self.vz
        vtot = math.sqrt(vtot_x ** 2 + vtot_y ** 2 + vtot_z ** 2)
        norm = math.sqrt(x ** 2 + (y - RSUN) ** 2 + z ** 2)
        v_para = (vtot_x * x + vtot_y * (y - RSUN) + vtot_z * z) / norm
        v_perp = math.sqrt(vtot ** 2 - v_para ** 2)

        # Compute proper motion from 2D v component perp to LOS
        pm = v_perp * 0.211 / self.dtrue
        
        return pm

    @property
    def width_ms(self):
        """W50 in ms"""
        return self.width_degree * self.period / 360.
    
