#!/usr/bin/python

import sys
import argparse
import math

import cPickle

from galacticops import GalacticOps
from population import Population 
from pulsar import Pulsar
from radialmodels import RadialModels
from survey import Survey

from progressbar import ProgressBar

import matplotlib.pyplot as plt



class View:
    """Generate basic histograms of pulsar population data."""
    def __init__(self, popfile):
        """Reads in the population file"""
        try:
            f = open(popfile, 'rb')
        except IOError:
            print "File {0} not found".format(popfile)
            sys.exit()

        self.pop = cPickle.load(f)
        f.close()

    def histogram(self, prop, log=False):
        """Create list and make histogram for the selected parameter"""

        if prop == 'period' :
            proplist = [p.period for p in self.pop.population]
        elif prop == 'dm':
            proplist = [p.dm for p in self.pop.population]
        elif prop == 'gl':
            proplist = [p.gl for p in self.pop.population]
        elif prop == 'gb':
            proplist = [p.gb for p in self.pop.population]
        elif prop == 'lum' :
            proplist = [p.lum_1400 for p in self.pop.population]
        elif prop == 'alpha':
            proplist = [p.alpha for p in self.pop.population]
        elif prop == 'r0':
             proplist = [math.fabs(p.r0) for p in self.pop.population]
        elif prop == 'rho':
            proplist = [p.rho for p in self.pop.population]
        elif prop == 'width':
            proplist = [p.width_degree for p in self.pop.population]
        elif prop == 'spindex':
            proplist = [p.spindex for p in self.pop.population]
        elif prop == 'scindex':
            proplist = [p.scindex for p in self.pop.population]
        elif prop == 'dist':
            proplist = [p.dtrue for p in self.pop.population]
        else:
            print "Property '{0}' not recognised".format(prop)
            sys.exit()

        # if asked for a log plot, take logs
        print "Plotting property '{0}'".format(prop)
        if log:
            proplist = [math.log10(p) for p in proplist]

        # show histogram
        plt.hist(proplist)
        plt.show()



if __name__ == '__main__':
    """ 'Main' function; read in options, then generate population"""
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='View a histogram of your population model')
    parser.add_argument('-f',
                        metavar = 'fname',
                        default='populate.model',
                        help='file containing population model')
    parser.add_argument('-p',
                        metavar = 'param',
                        required = True,
                        help = 'parameter to plot',
                        choices = ['period','dm','gl','gb','lum',
                                    'alpha','r0','rho','width',
                                    'spindex', 'scindex', 'dist'])
    parser.add_argument('--log', nargs='?', const=True, default=False,
                         help = 'flag to indicate a log scale plot')

    args = parser.parse_args()
    v = View(args.f)
    v.histogram(args.p, log=args.log)
