#!/usr/bin/python

import sys
import argparse
import math
import random

import cPickle

from population import Population 
from pulsar import Pulsar
from survey import Survey

class Detections:
    """Just a simple object to store survey detection summary"""
    def __init__(self,
                 ndet=None,
                 nsmear=None,
                 nout=None,
                 ntf=None):
        self.ndet = ndet
        self.nsmear=nsmear
        self.nout = nout
        self.nfaint = ntf

class DoSurvey:

    """
    Run a population model through a (set of) survey(s)

    """

    def __init__(self, popfile='populate.model'):
        with open(popfile, 'rb') as f:
            self.pop = cPickle.load(f)

        self.surveyPops = []

    def write(self, extension='.results', asc=False, summary=False):
        """Write a survey results population to a binary file."""

        for surv, survpop, detected in self.surveyPops:
            # create an output file
            s = ''.join([surv,'.results'])

            # write the survpop to the file
            with open(s,'wb') as output:
                cPickle.dump(survpop, output)

            # Write ascii file if required
            if asc:
                surv.pop.write_asc(surv+ '.det')

            if summary:
                # Write a summary file for the survey (if true)
                filename = ''.join([surv,'.summary'])
                s = 'Detected {0}'.format(detected.ndet)
                s = '\n'.join([s, 'Nsmear {0}'.format(detected.nsmear)])
                s = '\n'.join([s, 'Nfaint {0}'.format(detected.nfaint)])
                s = '\n'.join([s, 'Nout {0}'.format(detected.nout)])
                s = ''.join([s, '\n'])

                with open(filename, 'w') as output:
                    output.write(s)


    
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
                snr = s.SNRcalc(psr, self.pop)
                if snr > s.SNRlimit:
                    ndet += 1
                    psr.snr = snr
                    survpop.population.append(psr)

                elif snr == -1.0:
                    nsmear += 1
                elif snr == -2.0:
                    nout += 1
                else:
                    ntf += 1

            # report the results
            print "Number of pulsars detected by survey {0} = {1}".format(surv,ndet)
            print "Number too faint = {0}".format(ntf)
            print "Number smeared = {0}".format(nsmear)
            print "Number out = {0}".format(nout)
            print "\n"

            d = Detections(ndet=ndet, ntf=ntf, nsmear=nsmear, nout=nout)
            self.surveyPops.append([surv,survpop,d])


if __name__ == '__main__':
    """ 'Main' function; read in options, then survey the population"""
    # Parse command line arguments
     
    parser = argparse.ArgumentParser(description='Run a survey on your population model')
    parser.add_argument('-f', metavar = 'fname', default='populate.model', 
                         help='file containing population model')
    parser.add_argument('-surveys', metavar='S', nargs='+', required=True,
                         help='surveys to use to detect pulsars')
    parser.add_argument('--asc', nargs='?', const=True, default=False,
                         help='flag for ascii output')
    parser.add_argument('--summary', nargs='?', const=True, default=False,
                         help='flag for ascii summary file output')
    
    args = parser.parse_args()

    ds = DoSurvey(popfile=args.f)

    # run the survey and write the results to file
    ds.run(args.surveys)
    ds.write(asc=args.asc, summary=args.summary)
