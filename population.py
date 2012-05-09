#!/usr/bin/python

class Population:
    def __init__(self):
        self.population = []

        self._pDistType = None
        self._pmean = None
        self._psigma = None
        
        self._simean = None
        self._sisigma = None

        self._lummean = None
        self._lumsigma = None

        self._zscale = None

        self.ref_freq = 1400.0

        # non-properties
        self.ndet = 0

    def __str__(self):
        """Define how we print the population to screen"""
        s = "Population model:"
        s = '\n\t'.join([s, "Population size = {0}".format(self.size())])
        s = '\n\t'.join([s, "Reference Frequency = {0} MHz".format(self.ref_freq)])
        s = '\n\n\t'.join([s, "Mean period = {0}".format(self.pmean)])
        s = '\n\t'.join([s, "Period std dev = {0}".format(self.psigma)])
        s = '\n\t'.join([s, "Spectral index mean = {0}".format(self.simean)])
        s = '\n\t'.join([s, "Spectral index std dev = {0}".format(self.sisigma)])
        s = '\n\t'.join([s, "Luminosity mean = {0}".format(self.lummean)])
        s = '\n\t'.join([s, "Luminosity std dev = {0}".format(self.lumsigma)])
        return s

    def size(self):
        return len(self.population)
    
    # get/get methods for p distn type
    @property
    def pDistType(self):
        return self._pDistType

    @pDistType.setter
    def pDistType(self, value):
        self._pDistType = value

    @pDistType.deleter
    def pDistType(self):
        del self._pDistType


    # get/set methods for the mean period
    @property
    def pmean(self):
        return self._pmean

    @pmean.setter
    def pmean(self, value):
        self._pmean = value

    @pmean.deleter
    def pmean(self):
        del self._pmean

    # get/set methods for the period sigma
    @property
    def psigma(self):
        return self._psigma

    @psigma.setter
    def psigma(self, value):
        self._psigma = value

    @psigma.deleter
    def psigma(self):
        del self._psigma

    # get/set methods for the mean spindex
    @property
    def simean(self):
        return self._simean

    @simean.setter
    def simean(self, value):
        self._simean = value

    @simean.deleter
    def simean(self):
        del self._simean

    # get/set methods for the spindex sigma
    @property
    def sisigma(self):
        return self._sisigma

    @sisigma.setter
    def sisigma(self, value):
        self._sisigma = value

    @sisigma.deleter
    def sisigma(self):
        del self._sisigma

    # get/set methods for the mean luminosity
    @property
    def lummean(self):
        return self._lummean

    @lummean.setter
    def lummean(self, value):
        self._lummean = value

    @lummean.deleter
    def lummean(self):
        del self._lummean

    # get/set methods for the spindex luminosity
    @property
    def lumsigma(self):
        return self._lumsigma

    @lumsigma.setter
    def lumsigma(self, value):
        self._lumsigma = value

    @lumsigma.deleter
    def lumsigma(self):
        del self._lumsigma

    # get/set methods for the zscale
    @property
    def zscale(self):
        return self._zscale

    @zscale.setter
    def zscale(self, value):
        self._zscale = value

    @zscale.deleter
    def zscale(self):
        del self._zscale
