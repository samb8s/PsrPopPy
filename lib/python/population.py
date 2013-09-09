#!/usr/bin/python

import copy 
import cPickle

import numpy as np

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

        # list to store pulsar objects
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
        self.arguments = None

    def __str__(self):
        """Define how we print the population to screen."""
        s = "Population model:"

        for key, value in self.arguments.iteritems():
            s = '\n\t'.join([s, "{0} = {1}".format(key, value)])

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
        cPickle.dump(self, output, 2)
        output.close()

    def write_asc(self, outf):
        """Write population to an ascii file"""
        with open(outf, 'w') as f:
            titlestr = "Period_ms DM Width_ms GL GB S1400"
            titlestr = " ".join([titlestr, "L1400 SPINDEX SNR DTRUE X Y Z\n"])
            f.write(titlestr)
            for psr in self.population:
            
                s = "{0}".format(psr.period)
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

    def make_plotting_dicts(self):
        """
        Make an object that can be plotted using wxView, based on the attributes of 
        the pulsars in the population
        """

        labelDict = {
                    'Period'   : 'Period (ms)',
                    'Pdot'  : 'Period Derivative',
                    'DM'    : r'DM (cm$^{-3}$ pc)',
                    'Gal X' : 'X (kpc)',
                    'Gal Y' : 'Y (kpc)',
                    'Gal Z' : 'Z (kpc)',
                    'W'     : 'Width (degrees)' ,
                    'alpha' : 'alpha (deg)',
                    'rho'   : 'rho (deg)',
                    'SI'    : 'Spectral Index',
                    'S1400' : 'S1400 (mJy)',
                    'gl'    : 'Galactic Longitude (degrees)',
                    'gb'    : 'Galactic Latitude (degrees)',
                    'D'     : 'Distance (kpc)',
                    'r0'    : 'GalacticRadius (kpc)',
                    'n'     : 'Array Index'
                    }

        dataDict = {
            'Period': np.array([psr.period for psr in self.population]),
            'Pdot'  : np.array([psr.pdot for psr in self.population]),
            'DM'    : np.array([psr.dm for psr in self.population]),
            'Gal X' : np.array([psr.galCoords[0] for psr in self.population]),
            'Gal Y' : np.array([psr.galCoords[1] for psr in self.population]),
            'Gal Z' : np.array([psr.galCoords[2] for psr in self.population]),
            'W'     : np.array([psr.width_degree for psr in self.population]),
            'alpha' : np.array([psr.alpha for psr in self.population]),
            'rho'   : np.array([psr.rho for psr in self.population]),
            'SI'    : np.array([psr.spindex for psr in self.population]),
            'S1400' : np.array([psr.s_1400() for psr in self.population]),
            'gl'    : np.array([psr.gl for psr in self.population]),
            'gb'    : np.array([psr.gb for psr in self.population]),
            'D'     : np.array([psr.dtrue for psr in self.population]),
            'r0'    : np.array([psr.r0 for psr in self.population]),
            'n'     : np.array([x for x in range(len(self.population))]),
                    }

        return labelDict, dataDict
