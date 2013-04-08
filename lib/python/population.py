#!/usr/bin/python

import copy 
import cPickle

class Population:

    """
    Class to describe a pulsar population
    
    """

    def __init__(self,
                 pDistType=None,
                 radialDistType=None,
                 lumDistType=None,
                 pmean=None,
                 psigma=None,
                 simean=None,
                 sisigma=None,
                 lummean=None,
                 lumsigma=None,
                 zscale=None,
                 electronModel=None,
                 gpsFrac=None,
                 gpsA=None,
                 brokenFrac=None,
                 brokenSI=None,
                 ref_freq=1400.0):

        """Initialise the population object."""

        self.population = []

        # distribution types
        self.pDistType = pDistType
        self.radialDistType = radialDistType
        self.lumDistType = lumDistType
        self.electronModel=electronModel

        # distribution values
        self.pmean = pmean
        self.psigma = psigma
        
        self.simean = simean
        self.sisigma = sisigma

        self.lummean = lummean
        self.lumsigma = lumsigma

        self.zscale = zscale

        self.ref_freq = ref_freq

        # GPS and double SI values
        self.gpsFrac = gpsFrac
        self.gpsA = gpsA
        self.brokenFrac = brokenFrac
        self.brokenSI = brokenSI

        # non-"properties"
        self.ndet = 0

    def __str__(self):
        """Define how we print the population to screen."""
        s = "Population model:"
        s = '\n\t'.join([s, "Population size = {0}".format(self.size())])
        s = '\n\t'.join([s, "Reference Frequency = {0} MHz".format(self.ref_freq)])
        
        s = '\n\n\t'.join([s, "Mean period = {0}".format(self.pmean)])
        s = '\n\t'.join([s, "Period std dev = {0}".format(self.psigma)])
        s = '\n\t'.join([s, "Spectral index mean = {0}".format(self.simean)])
        s = '\n\t'.join([s, "Spectral index std dev = {0}".format(self.sisigma)])
        s = '\n\t'.join([s, "Luminosity mean = {0}".format(self.lummean)])
        s = '\n\t'.join([s, "Luminosity std dev = {0}".format(self.lumsigma)])

        s = '\n\n\t'.join([s, "Period Distribution = {0}".format(self.pDistType)])
        s = '\n\t'.join([s, "Radial distribution type = {0}".format(self.radialDistType)])
        s = '\n\t'.join([s, "Electron Model = {0}".format(self.electronModel)])

        return s

    def size(self):
        """Returns the number of pulsars in the object."""
        return len(self.population)

    def join(self,poplist):
        """
        Join pops in poplist to this population object. Returns a new object
        
        Be aware that population properties such as distribution type etc will be
        retained ONLY for the original population model.
        """

        # want to return a new object, so need to do a copy. 
        newpop = copy.deepcopy(self)

        for pop in poplist:
            # append each population to the original
            newpop.population = newpop.population + pop.population

        # return the new object
        return newpop

    def write(self, outf):
        """Write the population object to a file"""
        output = open(outf, 'wb')
        cPickle.dump(self, output)
        output.close()

    def write_asc(self, outf):
        """Write population to an ascii file"""
        with open(outf, 'w') as f:
            titlestr = "Period_ms DM Width_ms GL GB S1400"
            titlestr = " ".join([titlestr, "L1400 SPINDEX SNR DTRUE X Y Z\n"])
            f.write(titlestr)
            for psr in self.population:
            
                s = "{0}".format(self.period)
                s = "\t".join([s, "{0}".format(psr.dm)])
                s = "\t".join([s, "{0}".format(psr.width_ms())])
                s = "\t".join([s, "{0}".format(psr.gl)])
                s = "\t".join([s, "{0}".format(psr.gb)])
                s = "\t".join([s, "{0}".format(psr.s_1400())])
                s = "\t".join([s, "{0}".format(psr.lum_1400)])
                s = "\t".join([s, "{0}".format(psr.spindex)])
                s = "\t".join([s, "{0}".format(psr.snr)])
                s = "\t".join([s, "{0}".format(psr.dtrue)])
                s = "\t".join([s, "{0}".format(psr.galCoords[0])])
                s = "\t".join([s, "{0}".format(psr.galCoords[1])])
                s = "\t".join([s, "{0}".format(psr.galCoords[2])])
                s = "".join([s, "\n"])

                f.write(s)
