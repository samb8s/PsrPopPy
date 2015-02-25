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
                 ndisc=None,
                 nsmear=None,
                 nout=None,
                 ntf=None):
        self.ndet = ndet
        self.ndisc = ndisc
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
          nores=False, 
          asc=False, 
          summary=False):
    """Write a survey results population to a binary file."""

    for surv, survpop, detected in surveyPops:
        # create an output file, if required
        if not nores:
            if surv is not None:
                s = ''.join([surv,extension])

                survpop.write(s)
            else:
                s = 'allsurveys.results'
                survpop.write(s)

        # Write ascii file if required
        if asc and surv is not None:
            survpop.write_asc(surv+ '.det')

        if summary and surv is not None:
            # Write a summary file for the survey (if true)
            filename = ''.join([surv,'.summary'])
            s = 'Detected {0}'.format(detected.ndet)
            s = '\n'.join([s, 'Ndiscovered {0}'.format(detected.ndisc)])
            s = '\n'.join([s, 'Nsmear {0}'.format(detected.nsmear)])
            s = '\n'.join([s, 'Nfaint {0}'.format(detected.nfaint)])
            s = '\n'.join([s, 'Nout {0}'.format(detected.nout)])
            s = ''.join([s, '\n'])

            with open(filename, 'w') as output:
                output.write(s)
    
def run(pop,
        surveyList, 
        nostdout = False, 
        allsurveyfile = False,
        scint = False,
        accelsearch = False,
        jerksearch = False):

    """ Run the surveys and detect the pulsars."""

    # print the population
    if not nostdout:
        print "Running doSurvey on population..."
        print pop

    # loop over the surveys we want to run on the pop file
    surveyPops = []
    for surv in surveyList:
        s = Survey(surv)
        s.discoveries = 0
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

            # add scintillation, if required
            # modifying S/N rather than flux is sensible because then
            # a pulsar can have same flux but change S/N in repeated surveys
            if scint:
                snr = s.scint(psr, snr)

            if snr > s.SNRlimit:
                ndet += 1
                psr.snr = snr
                survpop.population.append(psr)
                
                # check if the pulsar has been detected in other 
                # surveys
                if not psr.detected:
                    # if not, set to detected and increment
                    # number of discoveries by the survey
                    psr.detected = True
                    s.discoveries += 1

            elif snr == -1.0:
                nsmear += 1
            elif snr == -2.0:
                nout += 1
            else:
                ntf += 1

        # report the results
        if not nostdout:
            print "Total pulsars in model = {0}".format(len(pop.population))
            print "Number detected by survey {0} = {1}".format(surv,ndet)
            print "Of which are discoveries = {0}".format(s.discoveries)
            print "Number too faint = {0}".format(ntf)
            print "Number smeared = {0}".format(nsmear)
            print "Number out = {0}".format(nout)
            print "\n"

        d = Detections(ndet=ndet, 
                       ntf=ntf, 
                       nsmear=nsmear, 
                       nout=nout,
                       ndisc=s.discoveries)
        surveyPops.append([surv,survpop,d])

    if allsurveyfile:
        allsurvpop = Population()
        allsurvpop.population = [psr for psr in pop.population if psr.detected]
        surveyPops.append([None, allsurvpop, None])

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

    parser.add_argument('--allsurveys', nargs='?', const=True, default=False,
                        help = 'write additional allsurv.results file (def=False)')

    parser.add_argument('--scint', nargs='?', const=True, default=False,
                        help = 'include model scintillation effects (def=False)')

    parser.add_argument('--accel', nargs='?', const=True, default=False,
                        help = 'use accel search for MSPs (def=False)')

    parser.add_argument('--jerk', nargs='?', const=True, default=False,
                        help = 'use accel & jerk search for MSPs (def=False)')

    args = parser.parse_args()

    # Load a model population
    population = loadModel(popfile=args.f)

    # run the population through the surveys
    surveyPopulations = run(population,
                            args.surveys,
                            nostdout=args.nostdout,
                            allsurveyfile=args.allsurveys,
                            scint=args.scint,
                            accelsearch=args.accel,
                            jerksearch=args.jerk)

    # write the output files
    write(surveyPopulations,
          nores=args.noresults,
          asc=args.asc, 
          summary=args.summary)
