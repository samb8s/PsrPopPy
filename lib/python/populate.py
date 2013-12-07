#!/usr/bin/python

import sys
import argparse
import math
import random

import inspect
import cPickle

import galacticops as go

from population import Population 
from pulsar import Pulsar
from survey import Survey

from progressbar import ProgressBar

class PopulateException(Exception):
    pass

def generate(ngen,
             surveyList=None,
             pDistType='lnorm',
             radialDistType='lfl06',
             radialDistPars=7.5,
             electronModel='ne2001',
             pDistPars=[2.7, -0.34],
             siDistPars= [-1.6,0.35], 
             lumDistType='lnorm',
             lumDistPars=[-1.1, 0.9],
             zscaleType='exp',
             zscale=0.33, 
             duty=6.,
             scindex=-3.86,
             gpsArgs=[None, None],
             doubleSpec=[None, None],
             nostdout=False,
             pattern='gaussian'):

    """
    Generate a population of pulsars.

    Keyword args:
    ngen -- the number of pulsars to generate (or detect)
    surveyList -- a list of surveys you want to use to try and detect the pulsars
    pDistType -- the pulsar period distribution model to use (def=lnorm)
    radialDistType -- radial distribution model
    electronModel -- mode to use for Galactic electron distribution
    pDistPars -- parameters to use for period distribution
    siDistPars -- parameters to use for spectral index distribution
    lumDistPars -- parameters to use for luminosity distribution
    radialDistPars -- parameters for radial distribution
    zscale -- if using exponential z height, set it here (in kpc)
    scindex -- spectral index of the scattering model
    gpsArgs -- add GPS-type spectrum sources
    doubleSpec -- add double-spectrum type sources
    nostdout -- (bool) switch off stdout
    """

    pop = Population()

    # check that the distribution types are supported....
    if lumDistType not in ['lnorm', 'pow']:
        print "Unsupported luminosity distribution: {0}".format(lumDistType)

    if pDistType not in ['lnorm', 'norm', 'cc97', 'lorimer12']:
        print "Unsupported period distribution: {0}".format(pDistType)

    if radialDistType not in ['lfl06', 'yk04', 'isotropic',
                                 'slab', 'disk', 'gauss']:
        print "Unsupported radial distribution: {0}".format(radialDistType)
    
    if electronModel not in ['ne2001', 'lmt85']:
        print "Unsupported electron model: {0}".format(electronModel)

    if pattern not in ['gaussian','airy']:
        print "Unsupported gain pattern: {0}".format(pattern)

    if duty<0.:
        print "Unsupported value of duty cycle: {0}".format(duty)

    # need to use properties in this class so they're get/set-type props
    pop.pDistType = pDistType
    pop.radialDistType = radialDistType
    pop.electronModel = electronModel
    pop.lumDistType = lumDistType

    pop.rsigma = radialDistPars
    pop.pmean, pop.psigma = pDistPars
    pop.simean, pop.sisigma = siDistPars

    pop.gpsFrac, pop.gpsA = gpsArgs
    pop.brokenFrac, pop.brokenSI = doubleSpec

    if pop.lumDistType == 'lnorm':
        pop.lummean, pop.lumsigma = \
                lumDistPars[0], lumDistPars[1]
    else:
        try:
            pop.lummin, pop.lummax, pop.lumpow = \
                lumDistPars[0], lumDistPars[1], lumDistPars[2]
        except ValueError:
            raise PopulateException('Not enough lum distn parameters')

    pop.zscaleType = zscaleType
    pop.zscale = zscale

    # store the dict of arguments inside the model. Could be useful.
    argspec = inspect.getargspec(generate)
    key_values = [(arg, locals()[arg]) for arg in argspec.args]
    pop.arguments = {key: value for (key,value) in key_values}

    if not nostdout:
        print "\tGenerating pulsars with parameters:"
        param_string_list = []
        for key, value in key_values:
            s = ": ".join([key, str(value)])
            param_string_list.append(s)

        # join this list of strings, and print it
        s = "\n\t\t".join(param_string_list)
        print "\t\t{0}".format(s)
            
        # set up progress bar for fun :)
        prog = ProgressBar(min_value = 0,
                           max_value=ngen,
                           width=65,
                           mode='dynamic')


    # create survey objects here and put them in a list
    if surveyList is not None:
        surveys = [Survey(s,pattern) for s in surveyList]
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

        #surv.gainpat=pattern      

    while pop.ndet < ngen:
        # Declare new pulsar object
        p = Pulsar()

        # period, alpha, rho, width distribution calls
        # Start creating the pulsar!
        if pop.pDistType == 'lnorm':
            p.period = _drawlnorm(pop.pmean, pop.psigma)
        elif pop.pDistType == 'norm':
            p.period = random.gauss(pop.pmean, pop.psigma)
        elif pop.pDistType == 'cc97':
            p.period = _cc97()
        elif pop.pDistType == 'gamma':
            print "Gamma function not yet supported"
            sys.exit()
        elif pop.pDistType == 'lorimer12':
            p.period = _lorimer2012_msp_periods()

        if duty>0.:
            # use a simple duty cycle for each pulsar
            # with a log-normal scatter
            width = (float(duty)/100.) * p.period**0.9
            width = math.log10(width)
            width = _drawlnorm(width, 0.3)

            p.width_degree = width*360./p.period
        else:
            # use the model to caculate if beaming
            p.alpha = _genAlpha()
            
            p.rho, p.width_degree = _genRhoWidth(p)

            if p.width_degree == 0.0 and p.rho ==0.0:
                continue
            # is pulsar beaming at us? If not, move on!
            p.beaming = _beaming(p)
            if not p.beaming:
                continue

        # Spectral index stuff here

        # suppose it might be nice to be able to have GPS sources 
        # AND double spectra. But for now I assume only have one or
        # none of these types.
        if random.random() > pop.gpsFrac:
            # This will evaluate true when gpsArgs[0] is NoneType
            # might have to change in future
            p.gpsFlag = 0
        else:
            p.gpsFlag = 1
            p.gpsA = pop.gpsA
            
        if random.random() > pop.brokenFrac:
            p.brokenFlag=0
        else:
            p.brokenFlag=1
            p.brokenSI = pop.brokenSI

        p.spindex = random.gauss(pop.simean, pop.sisigma)

        # get galactic position 
        # first, Galactic distribution models
        if pop.radialDistType == 'isotropic': 
            # calculate gl and gb randomly
            p.gb = math.degrees(math.asin(random.random()))
            if random.random() < 0.5:
                p.gb = 0.0 - p.gb
            p.gl = random.random() * 360.0

            # use gl and gb to compute galactic coordinates
            # pretend the pulsar is at distance of 1kpc
            # not sure why, ask Dunc!
            p.galCoords = go.lb_to_xyz(1.0, p.gl, p.gb)
        
        elif pop.radialDistType == 'slab':
            p.galCoords= go.slabDist()
            p.gl, p.gb = go.xyz_to_lb(p.galCoords)

        elif pop.radialDistType == 'disk':
            p.galCoords = go.diskDist()
            p.gl, p.gb = go.xyz_to_lb(p.galCoords)              

        else: # we want to use exponential z and a radial dist
            if pop.radialDistType == 'lfl06':
                p.r0 = go.lfl06()
            elif pop.radialDistType == 'yk04':
                p.r0 = go.ykr()
            elif pop.radialDistType == 'gauss':
                # guassian of mean 0
                # and stdDev given by parameter (kpc)
                p.r0 = random.gauss(0., pop.rsigma) 

            # then calc xyz,distance, l and b
            if pop.zscaleType == 'exp':
                zheight  = go._double_sided_exp(zscale)
            else:
                zheight = random.gauss(0., zscale)
            gx,gy  = go.calcXY(p.r0)
            p.galCoords = gx, gy, zheight
            p.gl, p.gb = go.xyz_to_lb(p.galCoords)

        p.dtrue = go.calc_dtrue(p.galCoords)

        # then calc DM  using fortran libs
        if pop.electronModel == 'ne2001':
            p.dm = go.ne2001_dist_to_dm(p.dtrue, p.gl, p.gb)
        elif pop.electronModel == 'lmt85':
            p.dm = go.lmt85_dist_to_dm(p.dtrue, p.gl, p.gb)

        p.scindex = scindex
        # then calc scatter time
        
        if pop.lumDistType == 'lnorm':
            p.lum_1400 = _drawlnorm(pop.lummean,
                                         pop.lumsigma)
        else:
            p.lum_1400 = _powerlaw(pop.lummin,
                                        pop.lummax,
                                        pop.lumpow)

        # if no surveys, just generate ngen pulsars
        if surveyList is None:
            pop.population.append(p)
            pop.ndet += 1
            if not nostdout:
                prog.increment_amount()
                print prog, '\r',
                sys.stdout.flush()
        # if surveys are given, check if pulsar detected or not
        # in ANY of the surveys
        else:
            detect_int = 0 # just a flag to increment if pulsar is detected
            for surv in surveys:
                # do SNR calculation
                SNR = surv.SNRcalc(p, pop)

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
            pop.population.append(p)

            # if detected, increment ndet (for whole population)
            # and redraw the progress bar
            if detect_int:
                pop.ndet += 1
                if not nostdout:
                    prog.increment_amount()
                    print prog, '\r',
                    sys.stdout.flush()

    # print info to stdout 
    if not nostdout:
        print "\n"
        print "  Total pulsars = {0}".format(len(pop.population))
        print "  Total detected = {0}".format(pop.ndet)
        #print "  Number not beaming = {0}".format(surv.nnb)

        for surv in surveys:
            print "\n  Results for survey '{0}'".format(surv.surveyName)
            print "    Number detected = {0}".format(surv.ndet)
            print "    Number too faint = {0}".format(surv.ntf)
            print "    Number smeared = {0}".format(surv.nsmear)
            print "    Number outside survey area = {0}".format(surv.nout)


    return pop


def _lorimer2012_msp_periods():
    """Picks a period at random from Dunc's
       distribution as mentioned in IAU 
       (China 2012) proceedings
    """
    # min max and n in distribution
    logpmin = 0.
    logpmax = 1.5
    dist = [1.,3.,5.,16.,9.,5.,5.,3.,2.]

    # calculate which bin to take value of
    #
    bin_num = _draw1d(dist)
    
    # assume linear distn inside the bins
    logp = logpmin + (logpmax-logpmin)*(bin_num+random.random())/len(dist)

    return 10.**logp

def _draw1d(dist):
    """Draw a bin number form a home-made distribution
        (dist is a list of numbers per bin)
    """
    # sum of distribution
    total = sum(dist)
    #cumulative distn
    cumulative = [sum(dist[:x+1])/total for x in range(len(dist))]

    rand_num = random.random()
    for i, c in enumerate(cumulative):
        if rand_num <= c:
            return i


def _drawlnorm(mean, sigma):
    """Get a random log-normal number."""
    return 10.0**random.gauss(mean, sigma)

def _powerlaw( minval, maxval, power):
    """Draw a value randomly from the specified power law"""

    logmin = math.log10(minval)
    logmax = math.log10(maxval)

    c = -1.0 * logmax * power
    nmax = 10.0**(power*logmin + c)

    # Dunc's code uses a goto statement
    # slightly worried about inf loops here...
    while True:
        log = logmin + (logmax-logmin)*random.random()
        n = 10.0**(power*log + c)

        if nmax*random.random() <= n:
            break

    return 10.0**log

def _cc97():
    """A model for MSP period distribution."""
    p = 0.0

    # pick p from the distribution, but cut off at 1 and 30 ms
    while p < 1.0 or p > 30.0:
        p = 0.65 / (1 - random.random())

    return p

def _beaming(psr):
    """Returns a boolean indicating if the pulsar is beaming at Earth."""
    # Emmering & Chevalier 89
    # find fraction of 4pi steradians that the pulsars beams to
    # for alpha and rho in degrees
    
    thetal = math.radians(max([0.0, psr.alpha-psr.rho]))
    thetau = math.radians(min([90.0, psr.alpha+psr.rho]))

    beamfrac = math.cos(thetal) - math.cos(thetau)

    # compare beamfrac vs a random number
    return random.random() < beamfrac

def _genRhoWidth(psr):
    """Calculate the opening angle of pulsar, and the beamwidth.
        Based on model outlined in Smits et al. 2009"""
    # cut off period for model
    perCut = 30.0

    # calclate rho
    randfactor = random.uniform(-.15, .15)
    if psr.period > perCut:
        rho = _rhoLaw(psr.period)
    else:
        rho = _rhoLaw(perCut)
    
    logrho = math.log10(rho) + randfactor
    rho = 10. ** logrho

    # generate beta and pulse width
    beta = random.uniform(-1,1) * rho
    width = _sindegree(0.5 * rho) * _sindegree(0.5 * rho)
    width = width - (_sindegree(0.5 * beta) * _sindegree(0.5 * beta))
    width = width /(_sindegree(psr.alpha) * _sindegree(psr.alpha + beta))

    if width < 0.0 or width > 1.0:
        width = 0.0
        rho = 0.0
    else:
        width = math.sqrt(width)

        # convert the width into degrees 0 -> 360 (ie. 90*4)
        width = math.degrees(math.asin(width))*4.0

    return rho, width

def _genAlpha():
    """Pick an inclination angle from 0-> 90 degrees.
        Based on model outlined in Smits et al. 2009"""
    angle = math.degrees(math.acos(random.uniform(0,1)))
    return math.fabs(angle)

def _rhoLaw(p_ms):
    """Calculate rho based on Rankin law of rho(period_ms)."""
    return 5.4 / math.sqrt(0.001 * p_ms)
    
def _sindegree(angle):
    """Return the sine of an angle in degrees."""
    return math.sin(math.radians(angle))


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
    # galactic-Z distn
    parser.add_argument('-zdist', nargs=1, required=False, default=['exp'],
                        help = 'type of distribution for z-scale',
                        choices=['exp', 'gauss'])
    parser.add_argument('-z', type=float, required=False, default=0.33,
                         help='exponential z-scale to use (def=0.33kpc)')

    # pulse width model
    parser.add_argument('-w', type=float, required=False, default=6,
                         help='pulse width, percent (def=6%) ')

    # spectral index distribution
    parser.add_argument('-si', nargs=2, type=float,
                         required=False, default=[-1.6, 0.35],
                         help='mean and std dev of spectral index distribution\
                                 (def = -1.6, 0.35)')
    # scattering index
    parser.add_argument('-sc', type=float, required=False, default=-3.86,
                        help='modify the frequency-dependance of Bhat et al\
                                scattering formula (def = -3.86)')

    # period distribution parameters
    parser.add_argument('-pdist', nargs=1, required=False, default=['lnorm'],
                        help='type of distribution to use for pulse periods',
                        choices=['lnorm', 'norm', 'cc97', 'lorimer12'])

    parser.add_argument('-p', nargs=2, required=False, type=float,
                         default=[2.7, -0.34],
                         help='period distribution mean and std dev \
                                 (def= [2.7, -0.34], Lorimer et al. 2006)')

    # luminosity distribution parameters
    parser.add_argument('-ldist', nargs=1, required=False, default=['lnorm'],
                        help='distribution to use for luminosities',
                        choices=['lnorm', 'pow'])
    parser.add_argument('-l', nargs='+', required=False, type=float,
                        default=[-1.1, 0.9],
                        help='luminosity distribution mean and std dev \
                             (def = [-1.1, 0.9], Faucher-Giguere&Kaspi, 2006)')

    # radial distribution type
    parser.add_argument('-rdist', type=str, nargs=1, required=False,
                        default=['lfl06'],
                        help='type of distrbution to use for Galactic radius',
                        choices=['lfl06', 'yk04', 'isotropic', 'slab',
                                   'disk', 'gauss'])

    parser.add_argument('-r', required = False,
                        default = 7.5, type=float,
                        help = 'radial distribution parameter \
                                (required for "-rdist gauss")')


    # electron/dm model
    parser.add_argument('-dm', type=str, nargs=1, required=False,
                        default=['ne2001'],
                        help='Galactic electron distribution model to use',
                        choices=['ne2001', 'lmt85'])
    
    # GPS sources
    parser.add_argument('-gps', type=float, nargs=2, required=False,
                        default=[None,None],
                        help='GPS fraction and "a" value')

    # double-spectral-index sources
    parser.add_argument('-doublespec', type=float, nargs=2, required=False,
                        default=[None,None],
                        help='Dbl spec fraction and alpha value')

    # output file name
    parser.add_argument('-o', type=str, metavar='outfile', required=False,
                        default='populate.model',
                        help='Output filename for population model')

    parser.add_argument('--nostdout', nargs='?', const=True, default=False,
                         help='flag to switch off std output (def=False)')

    args = parser.parse_args()


    # write command line to populate.cmd file
    with open('populate.cmd', 'a') as f:
        f.write(' '.join(sys.argv))
        f.write('\n')

    # run the code and write out a cPickle population class
    
    pop = generate(args.n,
                 surveyList=args.surveys,
                 pDistType=args.pdist[0],
                 lumDistType = args.ldist[0],
                 radialDistType=args.rdist[0],
                 radialDistPars=args.r,
                 pDistPars=args.p,
                 lumDistPars=args.l,
                 siDistPars=args.si,
                 zscaleType=args.zdist[0],
                 zscale=args.z,
                 duty=args.w,
                 scindex=args.sc,
                 electronModel=args.dm[0],
                 gpsArgs=args.gps,
                 doubleSpec = args.doublespec,
                 nostdout=args.nostdout
                 )

    pop.write(outf=args.o)
