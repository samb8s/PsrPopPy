#!/usr/bin/python

class Pulsar:
    """ Class to store pulsar information"""
    def __init__(self):
        """___init___ function for the Pulsar class"""
        self._period = None
        self._dm = None
        
        self._gl = None
        self._gb = None
        self._galCoords = None
        self._dtrue = None    

        self._lum1400 = None
        self._spindex = None
        self._scindex = None
        self._alpha = None
        self._rho = None
        self._width_degree = None

    # get/set methods for the period
    @property
    def period(self):
        return self._period

    @period.setter
    def period(self, value):
        self._period = value

    @period.deleter
    def period(self):
        del self._period

    # get/set methods for DM
    @property
    def dm(self):
        return self._dm

    @dm.setter
    def dm(self, value):
        self._dm = value

    @dm.deleter
    def dm(self):
        del self._dm

    # get/set methods for GL
    @property
    def gl(self):
        return self._gl

    @gl.setter
    def gl(self, value):
        self._gl = value

    @gl.deleter
    def gl(self):
        del self._gl

    # get/set methods for GB
    @property
    def gb(self):
        return self._gb

    @gb.setter
    def gb(self, value):
        self._gb = value

    @gb.deleter
    def gb(self):
        del self._gb

    # get/set methods for lum_1400
    @property
    def lum_1400(self):
        return self._lum_1400

    @lum_1400.setter
    def lum_1400(self, value):
        self._lum_1400 = value

    @lum_1400.deleter
    def lum_1400(self):
        del self._lum_1400

    # get/set methods for alpha
    @property
    def alpha(self):
        return self._alpha

    @alpha.setter
    def alpha(self, value):
        self._alpha = value

    @alpha.deleter
    def alpha(self):
        del self._alpha

    # get/set methods for rho
    @property
    def rho(self):
        return self._rho

    @rho.setter
    def rho(self, value):
        self._rho = value

    @rho.deleter
    def rho(self):
        del self._rho

    # get/set methods for the width (measured in phase)
    @property
    def width_degree(self):
        return self._width_degree

    @width_degree.setter
    def width_degree(self, value):
        self._width_degree = value

    @width_degree.deleter
    def width_degree(self):
        del self._width_degree

    # get/set methods for the spectral index
    @property
    def spindex(self):
        return self._spindex

    @spindex.setter
    def spindex(self, value):
        self._spindex = value

    @spindex.deleter
    def spindex(self):
        del self._spindex
  
    # get/set methods for spectral index of scattering 
    @property
    def scindex(self):
        return self._scindex

    @scindex.setter
    def scindex(self, value):
        self._scindex = value

    @scindex.deleter
    def scindex(self):
        del self._scindex

    # get /set methods for galactic XYZ
    @property
    def galCoords(self):
        return self._galCoords

    @galCoords.setter
    def galCoords(self, value):
        self._galCoords = value

    @galCoords.deleter
    def galCoords(self):
        del self._galCoords


    # get/set methods for true pulsar distance
    @property
    def dtrue(self):
        return self._dtrue

    @dtrue.setter
    def dtrue(self, value):
        self._dtrue = value

    @dtrue.deleter
    def dtrue(self):
        del self._dtrue

    # methods to calculate derived properties
    def s_1400(self):
        return self.lum_1400 / self.dtrue / self.dtrue

    def width_ms(self):
        return self.width_degree * self.period / 360.0
