#!/usr/bin/python

import sys
import argparse
import math
import random

import cPickle

from population import Population 
from pulsar import Pulsar
from survey import Survey

class DoSurvey:

    """
    Run a population model through a (set of) survey(s)

    """

    def __init__(self, popfile='populate.model'):
        f = open(popfile, 'rb')
        self.pop = cPickle.load(f)
        f.close()

        self.surveyPops = []

    def write(self, extension='.results'):
        """Write a survey results population to a binary file."""

        for surv, survpop in self.surveyPops:
            # create an output file
            s = surv + extension

            # write the survpop to the file
            output = open(s,'wb')
            cPickle.dump(survpop, output)
            output.close()
    
    def run(self, surveyList):
        """ Run the surveys and detect the pulsars."""

        # print the population
        print "Running doSurvey on population..."
        print self.pop

        # loop over the surveys we want to run on the pop file
        for surv in surveyList:
            s = Survey(surv)
            print "\nRunning survey {0}".format(surv)

            # create a new population object to store discovered pulsars in 
            survpop = Population()
            # HERE SHOULD INCLUDE THE PROPERTIES OF THE ORIGINAL POPULATION
            
            # counters 
            nsmear = 0
            nout = 0
            ntf = 0
            ndet = 0
            # loop over the pulsars in the population list
            for psr in self.pop.population:
                # is the pulsar over the detection threshold?
                if s.SNRcalc(psr, self.pop) > s.SNRlimit:
                    ndet += 1
                    survpop.population.append(psr)

                elif s.SNRcalc(psr, self.pop) == -1.0:
                    nsmear += 1
                elif s.SNRcalc(psr, self.pop) == -2.0:
                    nout += 1
                else:
                    ntf += 1

            # report the results
            print "Number of pulsars detected by survey {0} = {1}".format(surv,ndet)
            print "Number too faint = {0}".format(ntf)
            print "Number smeared = {0}".format(nsmear)
            print "Number out = {0}".format(nout)
            print "\n"

            self.surveyPops.append([surv,survpop])


if __name__ == '__main__':
    """ 'Main' function; read in options, then survey the population"""
    # Parse command line arguments
     
    parser = argparse.ArgumentParser(description='Run a survey on your population model')
    parser.add_argument('-f', metavar = 'fname', default='populate.model', 
                         help='file containing population model')
    parser.add_argument('-surveys', metavar='S', nargs='+', required=True,
                         help='surveys to use to detect pulsars')
    
    args = parser.parse_args()

    ds = DoSurvey(popfile=args.f)

    # run the survey and write the results to file
    ds.run(args.surveys)
    ds.write()
