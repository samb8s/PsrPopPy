#!/usr/bin/python

import math
import random

from orbit import Orbit
import galacticops as go


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
