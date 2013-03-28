#!/usr/bin/python

import sys
import argparse
import math
import random

import cPickle
import scipy.integrate
import numpy as np

import galacticops as go
import beaming as beammodels

from population import Population 
from pulsar import Pulsar
from survey import Survey

from progressbar import ProgressBar

class evolveException(Exception):
    pass

def write(pop, outf="evolve.model"):
    """Writes the population object to a binary file using cPickle."""
    output = open(outf,'wb')
    cPickle.dump(pop, output)
    output.close()

def generate(ngen,
             surveyList=None,
             age_max=1.0e9,
             pDistPars=[.3, .15],
             bFieldPars=[12.65, 0.55],
             birthvPars=[0.0, 180.],
             siDistPars= [-1.6,0.35],
             alignModel='orthogonal',
             spinModel = 'fk06',
             beamModel = 'tm98',
             spatialModel = 'spiralarms',
             birthVModel = 'gaussian',
             electronModel='ne2001',
             braking_index=0,
             zscale=0.33,
             deathline=True,
             nostdout=False):


    pop = Population()

    # set the parameters in the population object
    pop.pmean, pop.psigma = pDistPars
    pop.bmean, pop.bsigma = bFieldPars
    pop.simean, pop.sisigma = siDistPars
    pop.birthvmean, pop.birthvsigma = birthvPars
    
    pop.alignModel = alignModel
    pop.spinModel = spinModel
    pop.beamModel = beamModel
    pop.spatialModel = spatialModel
    pop.birthVModel = birthVModel
    pop.electronModel = electronModel

    pop.braking_index = braking_index
    pop.deathline = deathline

    pop.zscale = zscale

    if not nostdout:
        print "\tGenerating evolved pulsars with parameters:"
        print "\t\tngen = {0}".format(ngen)
        print "\t\tUsing electron distn model {0}".format(
                                        pop.electronModel)
        print "\n\t\tPeriod mean, sigma = {0}, {1}".format(
                                                    pop.pmean,
                                                    pop.psigma)
        print "\t\tLuminosity mean, sigma = {0}, {1}".format(
                                                    pop.lummean,
                                                    pop.lumsigma)
        print "\t\tSpectral index mean, sigma = {0}, {1}".format(
                                                    pop.simean,
                                                    pop.sisigma)
        print "\t\tGalactic z scale height = {0} kpc".format(
                                                    pop.zscale)

        # set up progress bar for fun :)
        prog = ProgressBar(min_value = 0,
                           max_value=ngen,
                           width=65,
                           mode='dynamic')

    # create survey objects here and put them in a list
    if surveyList is not None:
        surveys = [Survey(s) for s in surveyList]
    else:
        # make an empty list here - makes some code just a little
        # simpler - can still loop over an empty list (ie zero times)
        surveys=[]

    # initialise these counters to zero 
    for surv in surveys:
        surv.ndet =0 # number detected
        surv.nout=0 # number outside survey region
        surv.nsmear=0 # number smeared out
        surv.ntf=0 # number too faint

    # this is the nitty-gritty loop for generating the pulsars
    while pop.ndet < ngen:
        pulsar = Pulsar()

        # initial age for pulsar
        pulsar.age = random.random() * age_max

        # initial period
        pulsar.p0 = -1.
        while pulsar.p0 <= 0.:
            pulsar.p0 = random.gauss(pop.pmean, pop.psigma)

        # initial magnetic field (in Gauss)
        pulsar.bfield_init = 10**random.gauss(pop.bmean, pop.bsigma)

        # aligment angle
        alignpulsar(pulsar, pop.alignModel)

        # braking index
        if pop.braking_index == 0:
            pulsar.braking_index = 2.5 + 0.5 * random.random()
        else:
            pulsar.braking_index = float(pop.braking_index)

        #apply relevant spin down model
        pulsar.dead = False # pulsar should start alive! 
        if pop.spinModel == 'fk06':
            spindown_fk06(pulsar)

            # apply deathline if relevant
            if pop.deathline:
                bhattacharya_deathperiod_92(pulsar)

        elif pop.spinModel == 'cs06':
            # contopoulos and spitkovsky
            spindown_cs06(pulsar)

        # define pulse width as 5% (18 degrees)
        pulsar.width_degree = 18.

        if not pulsar.dead:
            print pulsar.dead
        # plough on - only if the pulsar isn't dead!
        if not pulsar.dead:
            # is the pulsar beaming? 
            pulsar_beaming(pulsar, pop.beamModel)
            # position of the pulsar       
            galacticDistribute(pulsar, pop)
            # birthvelocity
            birthVelocity(pulsar, pop)

            # model the xyz velocity
            go.vxyz(pulsar)

            # luminosity
            luminosity_fk06(pulsar)

            # spectral index
            pulsar.spindex = random.gauss(pop.simean, pop.sisigma)

            # calculate galactic coords and distance
            pulsar.gl, pulsar.gb = go.xyz_to_lb(pulsar.galCoords)
            pulsar.dtrue = go.calc_dtrue(pulsar.galCoords)

            # then calc DM  using fortran libs
            if pop.electronModel == 'ne2001':
                pulsar.dm = go.ne2001_dist_to_dm(pulsar.dtrue,
                                                   pulsar.gl, 
                                                   pulsar.gb)
            elif pop.electronModel == 'lmt85':
                pulsar.dm = go.lmt85_dist_to_dm(pulsar.dtrue,
                                                  pulsar.gl, 
                                                  pulsar.gb)

            # if no surveys, just generate ngen pulsars
            if surveyList is None:
                pop.population.append(pulsar)
                pop.ndet += 1

                if not nostdout:
                    prog.increment_amount()
                    print prog, '\r',
                    sys.stdout.flush()
            # if surveys are given, check if pulsar detected or not
            # in ANY of the surveys
            else:
                detect_int = 0 # just a flag if pulsar is detected
                for surv in surveys:
                    SNR = surv.SNRcalc(pulsar, pop)

                    if SNR > surv.SNRlimit:
                        # SNR is over threshold
                        # increment the flag 
                        # and survey ndetected
                        detect_int += 1
                        surv.ndet += 1
                        continue

                    elif SNR == -1:
                        # pulse is smeared out
                        surv.nsmear += 1
                        continue

                    elif SNR == -2:
                        # pulsar is outside survey region
                        surv.nout += 1
                        continue

                    else:
                        #pulsar is just too faint
                        surv.ntf += 1
                        continue 
            
                # add the pulsar to the population
                pop.population.append(pulsar)

                # if detected, increment ndet (for whole population)
                # and redraw the progress bar
                if detect_int:
                    pop.ndet += 1
                    if not nostdout:
                        prog.increment_amount()
                        print prog, '\r',
                        sys.stdout.flush()
    return pop

def luminosity_fk06( pulsar):
    """ Equation 14 from  Ridley & Lorimer """
    # variables to use in the equation
    alpha = -1.5
    beta = 0.5
    delta_l = random.gauss(0.0, 0.8)

    # the equation
    logL = math.log10(0.18) + alpha*math.log10(pulsar.period/1000.) + \
            beta*math.log10(pulsar.pdot * 1.0e15) + delta_l

    # set L
    pulsar.lum_1400 = 10.0**logL

def birthVelocity( pulsar, pop):
    """ Get a birth veolocity for the pulsar"""

    bvM = pop.birthVModel
    mean = pop.birthvmean
    sigma = pop.birthvsigma
    if bvM == 'gaussian':
        pulsar.vx = random.gauss(mean, sigma)
        pulsar.vy = random.gauss(mean, sigma)
        pulsar.vz = random.gauss(mean, sigma)
    elif bvM == 'exp':
        pulsar.vx = go._double_sided_exp(sigma, origin = mean)
        pulsar.vy = go._double_sided_exp(sigma, origin = mean)
        pulsar.vz = go._double_sided_exp(sigma, origin = mean)
    else:
        raise evolveException('Invalid velocity model selected')


def galacticDistribute(pulsar, pop):
    """ select a galactic position """
    if pop.spatialModel == 'spiralarms':
        r0 = go.ykr()
        x, y = go.spiralize(r0)
        z = go._double_sided_exp(pop.zscale)

        pulsar.galCoords = (x, y, z)
        pulsar.r0 = math.sqrt(x**2 + y**2)

def alignpulsar(pulsar, aM):
    """
    Pick an alignment angle for pulsar, depending on model chosen
    """
    if aM == 'orthogonal':
        pulsar.chi = 90.0
        pulsar.sinchi_init = 1.0
        pulsar.sinchi = pulsar.sinchi_init
        pulsar.coschi = 0.

    elif aM == 'random':
        chi = math.acos(random.random()) # in radians
        
        pulsar.chi = math.degrees(chi) # -> degrees
        pulsar.sinchi_init = math.sin(chi)
        pulsar.sinchi = pulsar.sinchi_init
        pulsar.coschi = math.cos(chi)
    
    else:
        raise evolveException('Invalid alignment model selected')

    # more models to add here, but it'll do for now

def pulsar_beaming(pulsar, bM):
    """
    Work out if the pulsar is beaming --- model-dependent
    """

    # Tauris & Manchester beaming model (default)
    if bM == 'tm98':
        fraction = beammodels.tm98_fraction(pulsar.period/1000.)
    ## more models to add here!
    elif bM == 'none':
        fraction = 1.0
    elif bM == 'const':
        fraction = 0.2

    if random.random()<fraction:
        pulsar.beaming = True
    else:
        pulsar.beaming = False

def spindown_fk06(pulsar):
    """
    Spindown model from Faucher-Giguere & Kaspi 2006
    """
    
    # k is a constant from the FK06 paper. 
    # defined in cgs as
    # k = (8 * pi**2 * R**6)/(3 * I * c**3)
    k = 9.768E-40 # in cgs
    kprime = k * pulsar.bfield_init**2

    # calculate p(t) - convert to milliseconds
    pulsar.period = _p_of_t_fk06(pulsar, kprime) * 1000.0

    # calculate pdot
    pulsar.pdot = _pdot_fk06(pulsar, kprime)

def _p_of_t_fk06(pulsar, kprime):
    """
    Equation 6 - Ridley & Lorimer 2010
    """
    index_term = pulsar.braking_index - 1.0
    age_s = pulsar.age * 3.15569e7 # convert to seconds

    brackets = pulsar.p0**(index_term) + index_term * kprime * age_s \
                    * pulsar.sinchi_init**2

    return brackets**(1.0/index_term)

def _pdot_fk06(pulsar, kprime):
    """
    Equation 3 - Ridley & Lorimer 2010
    """

    index_term = 2.0 - pulsar.braking_index
    period_s = pulsar.period / 1000.
    return kprime * pulsar.sinchi_init**2 * period_s**(index_term)

def spindown_cs06(pulsar):
    """
    Equations 9 and 10 in Ridley & Lorimer - 
    due to Contopoulos & Spitkovsky 2006
    """
    # this is equation 10
    index = (2.0 / (pulsar.braking_index + 1.0))
    pdeath = (0.81 * pulsar.bfield_init / 1.0E12 / pulsar.p0) ** index
    # convert to milliseconds
    pdeath *= 1000. 

    # this method needs integrals.
    #converting the qsimp (fortran) to using scipy.integrate package
    lower_limit = pulsar.p0*1000.
    upper_limit = 1.0E5
    const_of_integration = pulsar.coschi**2.0 / pdeath
    
    min_value = 1.0E24
    min_p = pulsar.p0*1000.

    index = pulsar.braking_index -3.0
    temp_const = 3.3E-40*pulsar.bfield_init**2 * (pulsar.p0*1000.)**index \
                        * pulsar.age * 3.15569e7

    # do the integral
    result = scipy.integrate.quad(_cs06_poft,
                                  lower_limit, 
                                  pdeath,
                                  args = (const_of_integration, 
                                          pulsar.braking_index)
                                  )[0]
    
    if result < temp_const or result>1.0E14:
        #print result, temp_const
        #print "step one"
        pulsar.period = 1.0E6
    else:
        looparray = np.arange(lower_limit, upper_limit+1)
        for m in looparray:
            result = scipy.integrate.quad(_cs06_poft,
                                          lower_limit, 
                                          pdeath,
                                          args = (const_of_integration, 
                                                  pulsar.braking_index)
                                         )[0]
            if result > 1.0E14:
                #print "step two"
                pulsar.period = 1.0E6
                break

            tempmin = math.fabs(result - temp_const)
            print tempmin, "= tempmin;", min_value, " = minval"
            if tempmin <= min_value:
                min_value = tempmin
                min_p = m
            else:
                pulsar.period = min_p
                break

            # end of loop
            if m == looparray[-1]:
                #print "step three"
                pulsar.period = 1.0E6

    #print pulsar.period, pdeath
    if pulsar.period>pdeath and pop.deathline:
        pulsar.dead = True
    else:
        pulsar.dead = False
        index = 2.0-pulsar.braking_index
        pulsar.pdot = 3.3E-40 * pulsar.bfield_init**2.0 * (1./pulsar.p0)\
                      *(pulsar.period / (pulsar.p0*1000.))**index * \
                      (1. - const_of_integration * pulsar.period)
                             
def _cs06_poft(x, a, n):
    return x**(n-2.0) / (1.0-a*x)


def bhattacharya_deathperiod_92(pulsar):
    """
    Eq 8 in Ridley & Lorimer - but it's the deathline from 
    Bhattacharya et al 1992
    """
    # deathline described by 
    # B/P^2 = 0.17E12 G s^-2

    # magnetic field = 3.2e19 * sqrt( P Pdot)
    B = 3.2E19 * math.sqrt(pulsar.period* pulsar.pdot/1000.)
    
    if B/(pulsar.period/1000.)**2 < 0.17E12:
        pulsar.dead = True
    

            
if __name__ == '__main__':
    """ 'Main' function; read in options, then generate population"""

    # set defaults here
    parser = argparse.ArgumentParser(description='Generate a population of pulsars')
    # number of pulsars to detect!
    parser.add_argument('-n', type=int, required=True,
                         help='number of pulsars to generate/detect')

    # list of surveys to use (if any)
    parser.add_argument('-surveys', metavar='S', nargs='+', default=None,
                         help='surveys to use to check if pulsars are detected'
                         )

    parser.add_argument('-tmax', type=float, required=False,
                         help = 'maximum initial age of pulsars')

    parser.add_argument('-spinmodel', type=str, required=False,
                        nargs=1, default=['fk06'],
                        choices=['fk06', 'cs06'],
                        help = 'spin-down model to employ')

    # output file name
    parser.add_argument('-o', type=str, metavar='outfile', required=False,
                        default='evolve.model',
                        help='Output filename for population model')

    parser.add_argument('--nostdout', nargs='?', const=True, default=False,
                         help='flag to switch off std output (def=False)')

    args = parser.parse_args()

    # write command line to populate.cmd file
    with open('evolve.cmd', 'a') as f:
        f.write(' '.join(sys.argv))
        f.write('\n')


    # run the code!
    pop = generate(args.n,
                    surveyList=args.surveys,
                    spinModel=args.spinmodel[0],
                    nostdout=args.nostdout)

    write(pop, outf=args.o)
