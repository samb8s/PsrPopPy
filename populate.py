#!/usr/bin/python

import sys
import argparse
import getopt #consider replacing with argparse
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
    def __init__(self):
        # creat epopulation object
        #self.population = Population()
        # generate pulsars
        self.pop = Population()
        self.pop.population=[]
        self.lums = []

    def write(self, outf="populate.model"):
        """Writes the population object to a binary file using cPickle"""
        output = open(outf,'wb')
        cPickle.dump(self.pop, output)
        output.close()

    def generate(self, 
                 ngen,
                 surveyList=None,
                 pDistType='lnorm',
                 pDistPars=[2.7, -0.34],
                 siDistPars= [-1.6,0.35], 
                 lumDistPars=[-1.1, 0.9],
                 zscale=0.33, 
                 scindex=-3.86):

        """Generate the pulsars"""

        # need to use properties in this class so they're get/set-type props
        self.pop.pDistType = pDistType
        self.pop.pmean, self.pop.psigma = pDistPars
        self.pop.simean, self.pop.sisigma = siDistPars
        self.pop.lummean, self.pop.lumsigma = lumDistPars
        self.pop.zscale = zscale

        print "\tGenerating pulsars with parameters:"
        print "\t\tngen = {0}".format(ngen)
        print "\t\tPeriod mean, sigma = {0}, {1}".format(self.pop.pmean,
                                                         self.pop.psigma)
        print "\t\tLuminosity mean, sigma = {0}, {1}".format(self.pop.lummean,
                                                             self.pop.lumsigma)
        print "\t\tSpectral index mean, sigma = {0}, {1}".format(self.pop.simean,
                                                                 self.pop.sisigma)
        print "\t\tGalactic z scale height = {0} kpc".format(self.pop.zscale)
        
        # set up progress bar for fun :)
        prog = ProgressBar(min_value = 0,max_value=ngen, width=65, mode='dynamic')

        nnb = 0
        ntf = 0
        nsmear = 0
        nout=0
        while self.pop.ndet < ngen:
            # Declare new pulsar object
            p = Pulsar()


            # period, alpha, rho, width distribution calls

            # check that the distribution type is supported....
            if self.pop.pDistType not in ['lnorm', 'norm']:
                #, 'gamma', '1d']:
                print "Unsupported period distribution type {0}".format(
                                                                self.pop.pDistType)
            if self.pop.pDistType == 'lnorm':
                p.period = self._drawlnorm(self.pop.pmean, self.pop.psigma)
            elif self.pop.pDistType == 'norm':
                p.period = random.gauss(self.pop.pmean, self.pop.psigma)
            elif self.pop.pDistType == 'gamma':
                print "Gamma function not yet supported"
                sys.exit()
                """
                p.period = self._drawgamma(0.0,
                                           10.0,
                                           self.pop.pmean,
                                           self.pop.psigma)
            elif self.pop.pDistType == '1d':
                p.period = self._draw1d
                """

            p.alpha = self._genAlpha()
            
            p.rho, p.width_degree = self._genRhoWidth(p)
            if p.width_degree == 0.0 and p.rho ==0.0:
                continue
            # is pulsar beaming at us? If not, move on!
            p.beaming = self._beaming(p)
            if not p.beaming:
                nnb += 1
                continue

            p.spindex = random.gauss(self.pop.simean, self.pop.sisigma)
            # get galactic radial position (method in radialmodels)
            p.r0 = self.lfl06()
            #p.r0 =  self.llfr() # python implementation (testing!)

            # then calc xyz,distance, l and b
            zheight = self._double_sided_exp(zscale)
            gx,gy  = self.calcXY(p.r0)
            p.galCoords = gx, gy, zheight

            p.dtrue = self.calc_dtrue(p.galCoords)
            p.gl, p.gb = self.xyz_to_lb(p.galCoords)

            # then calc DM  using fortran libs
            p.dm = self.ne2001_dist_to_dm(p.dtrue, p.gl, p.gb)
            p.scindex = scindex
            # then calc scatter time
            
            p.lum_1400 = self._drawlnorm(self.pop.lummean, self.pop.lumsigma)

            self.lums.append(p.s_1400())

            # if no surveys, just generate ngen pulsars
            if surveyList is None:
                self.pop.population.append(p)
                self.pop.ndet += 1
            # if surveys are given, check if pulsar detected or not
            # in ANY of the surveys
            else:
                for surv in surveyList:
                    # pulsar in survey region 
                    # is pulsar detectable in the survey
                    s = Survey(surv)
                    #print s.SNRcalc(p, self.pop), s.SNRlimit, s.inRegion(p)
                    if s.SNRcalc(p, self.pop) > s.SNRlimit:
                        self.pop.population.append(p)
                        self.pop.ndet += 1
                        prog.increment_amount()
                        print prog, '\r',
                        sys.stdout.flush()
                        # ok, the pulsar was detected in one of the surveys,
                        # so we can break out of the surveys loop now
                        break
                    elif s.SNRcalc(p, self.pop) == -1.0:
                        nsmear += 1
                        self.pop.population.append(p)
                        break
                    elif s.SNRcalc(p, self.pop) == -2.0:
                        nout += 1
                        self.pop.population.append(p)
                        break
                    else:
                        ntf += 1
                        self.pop.population.append(p)
                        break
               # print p.lum_1400

        print "\n\n"
        print "  Total pulsars = {0}".format(len(self.pop.population))
        print "  Number detected = {0}".format(self.pop.ndet)
        print "  Number not beaming = {0}".format(nnb)
        print "  Number too faint = {0}".format(ntf)
        print "  Number smeared = {0}".format(nsmear)
        print "  Number outside survey area = {0}".format(nout)


    def _double_sided_exp(self, scale, origin=0.0):
        """Exponential distribution around origin, with scale height scale"""
        rn1 = random.random()
        rn2 = random.random()
        sgn = 1.0
        if rn1 < 0.5:
            sgn = -1.0

        return origin + sgn * scale * math.log(rn2)

    def _drawlnorm(self, mean, sigma):
        """Get a random log-normal number to represent period in milliseconds"""
        #print mean, sigma
        return 10.0**random.gauss(mean, sigma)
    
    def _genAlpha(self):
        """Pick an inclination angle from 0-> 90 degree uniform omega"""
        angle = math.degrees(math.asin(random.uniform(-1,1)))
        return math.fabs(angle)

    def _rhoLaw(self, p_ms):
        return 5.4 / math.sqrt(0.001 * p_ms)

    def _sindegree(self, angle):
        return math.sin(math.radians(angle))

    def _beaming(self, psr):
        """Is the pulsar beaming to us? Boolean"""
        # Emmering & Chevalier 89
        # find fraction of 4pi steradians that the pulsars beams to
        # for alpha and rho in degrees
        
        thetal = math.radians(max([0.0, psr.alpha - psr.rho]))
        thetau = math.radians(min([90.0, psr.alpha + psr.rho]))

        beamfrac = math.cos(thetal) - math.cos(thetau)

        # compare beamfrac vs a random number
        return random.random() < beamfrac
        #if random.uniform(0,1) < beamfrac:
        #    return True
        #else:
        #    return False

    def _genRhoWidth(self, psr):
        """Get a pulsar opening angle"""
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
            # convert the width into degrees
            width = math.degrees(math.asin(width)*4.0) 

        return rho, width


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
    parser.add_argument('-pdist', nargs=1, required=False, default='lnorm',
                        help='type of distribution to use for pulse periods',
                        choices=['lnorm', 'norm'])#, 'gamma', '1d'])
    parser.add_argument('-p', nargs=2, required=False,
                         default=[2.7, -0.34],
                         help='period distribution mean and std dev \
                                 (def= [-2.7, -0.34])')

    args = parser.parse_args()

    # run the code and write out a cPickle population class
    pop = Populate()
    
    pop.generate(args.n,
                 surveyList=args.surveys,
                 pDistPars=args.p,
                 siDistPars=args.si,
                 zscale=args.z,
                 scindex=args.sc
                 )

    pop.write()




    # some stuff for testing
    """
    l = [p.lum_1400 for p in populate.pop.population]
    l2 = [math.log10(x) for x in populate.lums]
    plt.hist(l2, alpha=0.7, normed=True)
    #plt.hist(l3, alpha=0.3, normed=True)


    plt.show()
    """
