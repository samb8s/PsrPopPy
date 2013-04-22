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


def loadModel(popfile='populate.model', popmodel=None):
    """Loads in either a model from disk (popfile, cPickle), 
       or pass in a model from memory (popmodel)"""
    if popmodel is None:
        with open(popfile, 'rb') as f:
            pop = cPickle.load(f)
    else:
        pop = popmodel

    return pop

def write(surveyPops,
          extension='.results', 
          nores=True, 
          asc=False, 
          summary=False):
    """Write a survey results population to a binary file."""

    for surv, survpop, detected in surveyPops:
        # create an output file, if required
        if not nores:
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
    
def run(pop, surveyList, nostdout=False, pattern='gaussian'):
    """ Run the surveys and detect the pulsars."""

    # print the population
    if not nostdout:
        print "Running doSurvey on population..."
        print pop

    # loop over the surveys we want to run on the pop file
    surveyPops = []
    for surv in surveyList:
        s = Survey(surv,pattern)
        if not nostdout:
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
        for psr in pop.population:
            # pulsar could be dead (evolve!) - continue if so
            if psr.dead:
                continue

            # is the pulsar over the detection threshold?
            snr = s.SNRcalc(psr, pop)
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
        if not nostdout:
            print "Number detected by survey {0} = {1}".format(surv,ndet)
            print "Number too faint = {0}".format(ntf)
            print "Number smeared = {0}".format(nsmear)
            print "Number out = {0}".format(nout)
            print "\n"

        d = Detections(ndet=ndet, ntf=ntf, nsmear=nsmear, nout=nout)
        surveyPops.append([surv,survpop,d])

    return surveyPops


if __name__ == '__main__':
    """ 'Main' function; read in options, then survey the population"""
    # Parse command line arguments
     
    parser = argparse.ArgumentParser(description='Run a survey on your population model')
    parser.add_argument('-f', metavar = 'fname', default='populate.model', 
                         help='file containing population model (def=populate.model')

    parser.add_argument('-surveys', metavar='S', nargs='+', required=True,
                         help='surveys to use to detect pulsars (required)')

    parser.add_argument('--noresults', nargs='?', const=True, default=False,
                         help='flag to switch off pickled .results file (def=False)')

    parser.add_argument('--asc', nargs='?', const=True, default=False,
                         help='flag to create ascii population file (def=False)')
    
    parser.add_argument('--summary', nargs='?', const=True, default=False,
                         help='flag to create ascii summary file (def=False)')
    
    parser.add_argument('--nostdout', nargs='?', const=True, default=False,
                         help='flag to switch off std output (def=False)')

    args = parser.parse_args()

    # Load a model population
    population = loadModel(popfile=args.f)

    # run the population through the surveys
    surveyPopulations = run(population,
                            args.surveys,
                            nostdout=args.nostdout)

    # write the output files
    write(surveyPopulations,
          nores=args.noresults,
          asc=args.asc, 
          summary=args.summary)
