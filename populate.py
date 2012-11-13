#!/usr/bin/python

import sys
import argparse
import math
import random

import cPickle

from galacticops import GalacticOps
from population import Population 
from pulsar import Pulsar
from radialmodels import RadialModels
from survey import Survey

from progressbar import ProgressBar



class Populate(RadialModels, GalacticOps):

    """
    A class to generate a population of pulsars.

    """

    def __init__(self):
        """Initialise the populate class"""
        # create population object
        self.pop = Population()

    def write(self, outf="populate.model"):
        """Writes the population object to a binary file using cPickle."""
        output = open(outf,'wb')
        cPickle.dump(self.pop, output)
        output.close()

    def generate(self, 
                 ngen,
                 surveyList=None,
                 pDistType='lnorm',
                 radialDistType='lfl06',
                 electronModel='ne2001',
                 pDistPars=[2.7, -0.34],
                 siDistPars= [-1.6,0.35], 
                 lumDistType='lnorm',
                 lumDistPars=[-1.1, 0.9],
                 zscale=0.33, 
                 duty=0,
                 scindex=-3.86,
                 gpsArgs=[None, None],
                 doubleSpec=[None, None],
                 nostdout=False):

        """
        Generate a population of pulsars.

        Keyword args:
        ngen -- the number of pulsars to generate (or detect)
        surveyList -- a list of surveys you want to use to try and detect the pulsars
        pDistType -- the pulsar period distribution model to use (def=lnorm)
        electronModel -- mode to use for Galactic electron distribution
        pDistPars -- parameters to use for period distribution
        siDistPars -- parameters to use for spectral index distribution
        lumDistPars -- parameters to use for luminosity distribution
        zscale -- if using exponential z height, set it here (in kpc)
        scindex -- spectral index of the scattering model
        """

        # check that the distribution types are supported....
        if lumDistType not in ['lnorm', 'pow']:
            print "Unsupported luminosity distribution: {0}".format(lumDistType)

        if pDistType not in ['lnorm', 'norm', 'cc97']:
            print "Unsupported period distribution: {0}".format(pDistType)

        if radialDistType not in ['lfl06', 'yk04', 'isotropic',
                                     'slab', 'disk']:
            print "Unsupported radial distribution: {0}".format(radialDistType)
        
        if electronModel not in ['ne2001', 'lm98']:
            print "Unsupported electron model: {0}".format(electronModel)

        if duty<0.:
            print "Unsupported value of duty cycle: {0}".format(duty)

        # need to use properties in this class so they're get/set-type props
        self.pop.pDistType = pDistType
        self.pop.radialDistType = radialDistType
        self.pop.electronModel = electronModel
        self.pop.lumDistType = lumDistType
    
        self.pop.pmean, self.pop.psigma = pDistPars
        self.pop.simean, self.pop.sisigma = siDistPars

        self.pop.gpsFrac, self.pop.gpsA = gpsArgs
        self.pop.brokenFrac, self.pop.brokenSI = doubleSpec

        if self.pop.lumDistType == 'lnorm':
            self.pop.lummean, self.pop.lumsigma = lumDistPars
        else:
            self.pop.lummin, self.pop.lummax, self.pop.lumpow = lumDistPars

        self.pop.zscale = zscale

        if not nostdout:
            print "\tGenerating pulsars with parameters:"
            print "\t\tngen = {0}".format(ngen)
            print "\t\tUsing electron distn model {0}".format(
                                            self.pop.electronModel)
            print "\n\t\tPeriod mean, sigma = {0}, {1}".format(
                                                        self.pop.pmean,
                                                        self.pop.psigma)
            print "\t\tLuminosity mean, sigma = {0}, {1}".format(
                                                        self.pop.lummean,
                                                        self.pop.lumsigma)
            print "\t\tSpectral index mean, sigma = {0}, {1}".format(
                                                        self.pop.simean,
                                                        self.pop.sisigma)
            print "\t\tGalactic z scale height = {0} kpc".format(
                                                        self.pop.zscale)

            print "\t\tWidth {0}% -- (0 == model)".format(duty)
        
            if self.pop.gpsFrac:
                print "\n\t\tGPS Fraction = {0}, a = {1}".format(
                                                        self.pop.gpsFrac,
                                                        self.pop.gpsA)
            if self.pop.brokenFrac:
                print "\n\t\tDbl Spectrum Fraction = {0}, a = {1}".format(
                                                        self.pop.brokenFrac,
                                                        self.pop.brokenSI)

            # set up progress bar for fun :)
            prog = ProgressBar(min_value = 0,
                               max_value=ngen,
                               width=65,
                               mode='dynamic')

        nnb = 0
        ntf = 0
        nsmear = 0
        nout=0

        # create survey objects here and put them in a list
        if surveyList is not None:
            surveys = [Survey(surv) for surv in surveyList]

        while self.pop.ndet < ngen:
            # Declare new pulsar object
            p = Pulsar()

            # period, alpha, rho, width distribution calls
            # Start creating the pulsar!
            if self.pop.pDistType == 'lnorm':
                p.period = self._drawlnorm(self.pop.pmean, self.pop.psigma)
            elif self.pop.pDistType == 'norm':
                p.period = random.gauss(self.pop.pmean, self.pop.psigma)
            elif self.pop.pDistType == 'cc97':
                p.period = self._cc97()
            elif self.pop.pDistType == 'gamma':
                print "Gamma function not yet supported"
                sys.exit()

            if duty>0.:
                # use a simple duty cycle for each pulsar
                # with a log-normal scatter
                width = (float(duty)/100.) * p.period**0.9
                width = math.log10(width)
                width = self._drawlnorm(width, 0.3)

                p.width_degree = width*360./p.period
            else:
                # use the model to caculate if beaming
                p.alpha = self._genAlpha()
                
                p.rho, p.width_degree = self._genRhoWidth(p)

                if p.width_degree == 0.0 and p.rho ==0.0:
                    continue
                # is pulsar beaming at us? If not, move on!
                p.beaming = self._beaming(p)
                if not p.beaming:
                    nnb += 1
                    continue

            # Spectral index stuff here

            # suppose it might be nice to be able to have GPS sources 
            # AND double spectra. But for now I assume only have one or
            # none of these types.
            if random.random() > self.pop.gpsFrac:
                # This will evaluate true when gpsArgs[0] is NoneType
                # might have to change in future
                p.gpsFlag = 0
            else:
                p.gpsFlag = 1
                p.gpsA = self.pop.gpsA
                
            if random.random() > self.pop.brokenFrac:
                p.brokenFlag=0
            else:
                p.brokenFlag=1
                p.brokenSI = self.pop.brokenSI

            p.spindex = random.gauss(self.pop.simean, self.pop.sisigma)

            # get galactic position (method in radialmodels)
            # first, Galactic distribution models
            if self.pop.radialDistType == 'isotropic': 
                # calculate gl and gb randomly
                p.gb = math.degrees(math.asin(random.random()))
                if random.random() < 0.5:
                    p.gb = 0.0 - p.gb
                p.gl = random.random() * 360.0

                # use gl and gb to compute galactic coordinates
                # pretend the pulsar is at distance of 1kpc
                # not sure why, ask Dunc!
                p.galCoords = self.lb_to_xyz(1.0, p.gl, p.gb)
            
            elif self.pop.radialDistType == 'slab':
                p.galCoords= self.slabDist()
                p.gl, p.gb = self.xyz_to_lb(p.galCoords)

            elif self.pop.radialDistType == 'disk':
                p.galCoords = self.diskDist()
                p.gl, p.gb = self.xyz_to_lb(p.galCoords)              

            else: # we want to use exponential z and a radial dist
                if self.pop.radialDistType == 'lfl06':
                    p.r0 = self.lfl06()
                elif self.pop.radialDistType == 'yk04':
                    p.r0 = self.ykr()

                # then calc xyz,distance, l and b
                zheight = self._double_sided_exp(zscale)
                gx,gy  = self.calcXY(p.r0)
                p.galCoords = gx, gy, zheight
                p.gl, p.gb = self.xyz_to_lb(p.galCoords)

            p.dtrue = self.calc_dtrue(p.galCoords)

            # then calc DM  using fortran libs
            if self.pop.electronModel == 'ne2001':
                p.dm = self.ne2001_dist_to_dm(p.dtrue, p.gl, p.gb)
            elif self.pop.electronModel == 'lm98':
                p.dm = self.lm98_dist_to_dm(p.dtrue, p.gl, p.gb)

            p.scindex = scindex
            # then calc scatter time
            
            if self.pop.lumDistType == 'lnorm':
                p.lum_1400 = self._drawlnorm(self.pop.lummean,
                                             self.pop.lumsigma)
            else:
                p.lum_1400 = self._powerlaw(self.pop.lummin,
                                            self.pop.lummax,
                                            self.pop.lumpow)

            # if no surveys, just generate ngen pulsars
            if surveyList is None:
                self.pop.population.append(p)
                self.pop.ndet += 1
            # if surveys are given, check if pulsar detected or not
            # in ANY of the surveys
            else:
                for surv in surveys:
                    # pulsar in survey region 
                    # is pulsar detectable in the survey?
                    SNR = surv.SNRcalc(p, self.pop)

                    if SNR > surv.SNRlimit:
                        self.pop.population.append(p)
                        self.pop.ndet += 1
                        if not nostdout:
                            prog.increment_amount()
                            print prog, '\r',
                            sys.stdout.flush()
                        # the pulsar was detected in one of the surveys,
                        # so we can break out of the surveys loop now
                        break
                    elif SNR == -1:
                        nsmear += 1
                        self.pop.population.append(p)
                        break
                    elif SNR == -2:
                        nout += 1
                        self.pop.population.append(p)
                        break
                    else:
                        ntf += 1
                        self.pop.population.append(p)
                        break
               # print p.lum_1400

        if not nostdout:
            print "\n\n"
            print "  Total pulsars = {0}".format(len(self.pop.population))
            print "  Number detected = {0}".format(self.pop.ndet)
            print "  Number not beaming = {0}".format(nnb)
            print "  Number too faint = {0}".format(ntf)
            print "  Number smeared = {0}".format(nsmear)
            print "  Number outside survey area = {0}".format(nout)


    def _double_sided_exp(self, scale, origin=0.0):
        """Exponential distribution around origin, with scale height scale."""
        if scale == 0.0:
            return origin

        rn = random.random()
        sign = random.choice([-1.0, 1.0])

        return origin + sign * scale * math.log(rn)

    def _drawlnorm(self, mean, sigma):
        """Get a random log-normal number."""
        #print mean, sigma
        return 10.0**random.gauss(mean, sigma)

    def _powerlaw(self, minval, maxval, power):
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
    
    def _cc97(self):
        """A model for MSP period distribution."""
        p = 0.0

        # pick p from the distribution, but cut off at 1 and 30 ms
        while p < 1.0 or p > 30.0:
            p = 0.65 / (1 - random.random())

        return p

    def _beaming(self, psr):
        """Returns a boolean indicating if the pulsar is beaming at Earth."""
        # Emmering & Chevalier 89
        # find fraction of 4pi steradians that the pulsars beams to
        # for alpha and rho in degrees
        
        thetal = math.radians(max([0.0, psr.alpha-psr.rho]))
        thetau = math.radians(min([90.0, psr.alpha+psr.rho]))

        beamfrac = math.cos(thetal) - math.cos(thetau)

        # compare beamfrac vs a random number
        return random.random() < beamfrac

    def _genRhoWidth(self, psr):
        """Calculate the opening angle of pulsar, and the beamwidth."""
        # cut off period for model
        perCut = 30.0
        drho = 0.3

        # calclate rho
        randfactor = random.uniform(-.5, .5) * drho
        if psr.period > perCut:
            rho = self._rhoLaw(psr.period)
        else:
            rho = self._rhoLaw(perCut)
        
        logrho = math.log10(rho) + randfactor
        rho = 10. ** logrho

        # generate beta and pulse width
        beta = random.uniform(-1,1) * rho
        width = self._sindegree(0.5 * rho) * self._sindegree(0.5 * rho)
        width = width - (self._sindegree(0.5 * beta) * self._sindegree(0.5 * beta))
        width = width /(self._sindegree(psr.alpha) * self._sindegree(psr.alpha + beta))

        if width < 0.0 or width > 1.0:
            width = 0.0
            rho = 0.0
        else:
            width = math.sqrt(width)

            # convert the width into degrees 0 -> 360 (ie. 90*4)
            width = math.degrees(math.asin(width))*4.0

        return rho, width

    def _genAlpha(self):
        """Pick an inclination angle from 0-> 90 degrees."""
        angle = math.degrees(math.asin(random.uniform(-1,1)))
        return math.fabs(angle)

    def _rhoLaw(self, p_ms):
        """Calculate rho based on Rankin law of rho(period_ms)."""
        return 5.4 / math.sqrt(0.001 * p_ms)
        
    def _sindegree(self, angle):
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
    parser.add_argument('-z', type=float, required=False, default=0.33,
                         help='exponential z-scale to use (def=0.33kpc)')

    # galactic-Z distn
    parser.add_argument('-w', type=float, required=False, default=0,
                         help='pulse width % (def=0=rankin model)')

    # spectral index distribution
    parser.add_argument('-si', nargs=2, type=float,
                         required=False, default=[-1.6, 0.35],
                         help='mean and std dev of spectral index distribution \
                                 (def = -1.6, 0.35)')
    # scattering index
    parser.add_argument('-sc', type=float, required=False, default=-3.86,
                        help='modify the frequency-dependance of Bhat et al \
                                scattering formula (def = -3.86)')

    # period distribution parameters
    parser.add_argument('-pdist', nargs=1, required=False, default=['lnorm'],
                        help='type of distribution to use for pulse periods',
                        choices=['lnorm', 'norm', 'cc97'])#, 'gamma', '1d'])
    parser.add_argument('-p', nargs=2, required=False,
                         default=[2.7, -0.34],
                         help='period distribution mean and std dev \
                                 (def= [2.7, -0.34], Lorimer et al. 2006)')

    # luminosity distribution parameters
    parser.add_argument('-ldist', nargs=1, required=False, default=['lnorm'],
                        help='distribution to use for luminosities',
                        choices=['lnorm', 'pow'])
    parser.add_argument('-l', nargs=2, required=False,
                        default=[-1.1, 0.9],
                        help='luminosity distribution mean and std dev \
                                (def = [-1.1, 0.9], Faucher-Giguere&Kaspi, 2006)'

    # radial distribution type
    parser.add_argument('-rdist', type=str, nargs=1, required=False, default=['lfl06'],
                        help='type of distrbution to use for Galactic radius',
                        choices=['lfl06', 'yk04', 'isotropic', 'slab', 'disk'])

    # electron/dm model
    parser.add_argument('-dm', type=str, nargs=1, required=False,
                        default=['ne2001'],
                        help='Galactic electron distribution model to use',
                        choices=['ne2001', 'lm98'])
    
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
    pop = Populate()
    
    pop.generate(args.n,
                 surveyList=args.surveys,
                 pDistType=args.pdist[0],
                 radialDistType=args.rdist[0],
                 pDistPars=args.p,
                 siDistPars=args.si,
                 zscale=args.z,
                 duty=args.w,
                 scindex=args.sc,
                 electronModel=args.dm[0],
                 gpsArgs=args.gps,
                 doubleSpec = args.doublespec,
                 nostdout=args.nostdout
                 )

    pop.write(outf=args.o)
