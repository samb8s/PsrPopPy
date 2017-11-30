#!/usr/bin/python

import sys
import argparse
import math
import random

import inspect
import cPickle
import scipy.integrate
import numpy as np

import distributions as dists
import galacticops as go
import beaming as beammodels
import populate

from population import Population
from pulsar import Pulsar
from survey import Survey

from progressbar import ProgressBar


class EvolveException(Exception):
    pass


def generate(ngen,
             surveyList=None,
             age_max=1.0E9,
             pDistPars=[.3, .15],
             bFieldPars=[12.65, 0.55],
             birthVPars=[0.0, 180.],
             siDistPars=[-1.6, 0.35],
             alignModel='orthogonal',
             lumDistType='fk06',
             lumDistPars=[-1.5, 0.5],
             alignTime=None,
             spinModel='fk06',
             beamModel='tm98',
             birthVModel='gaussian',
             electronModel='ne2001',
             braking_index=0,
             zscale=0.05,
             duty=5.,
             scindex=-3.86,
             widthModel=None,
             nodeathline=False,
             efficiencycut=None,
             nostdout=False,
             nospiralarms=False,
             keepdead=False):

    pop = Population()

    # set the parameters in the population object
    pop.pmean, pop.psigma = pDistPars
    pop.bmean, pop.bsigma = bFieldPars

    if lumDistType == 'pow':
        try:
            pop.lummin, pop.lummax, pop.lumpow = \
                lumDistPars[0], lumDistPars[1], lumDistPars[2]
        except ValueError:
            raise EvolveException('Not enough lum distn parameters for "pow"')

    elif lumDistType == 'fk06':
        pop.lumPar1, pop.lumPar2 = lumDistPars[0], lumDistPars[1]
        if len(lumDistPars) == 3:
            pop.lumPar3 = lumDistPars[2]
        else:
            pop.lumPar3 = 0.18

    else:
        pop.lumPar1, pop.lumPar2 = lumDistPars

    pop.simean, pop.sisigma = siDistPars
    pop.birthvmean, pop.birthvsigma = birthVPars

    pop.alignModel = alignModel
    pop.alignTime = alignTime
    pop.spinModel = spinModel
    pop.beamModel = beamModel
    pop.birthVModel = birthVModel
    pop.electronModel = electronModel

    pop.braking_index = braking_index
    pop.deathline = not nodeathline
    pop.nospiralarms = nospiralarms

    pop.zscale = zscale

    if widthModel == 'kj07':
        print "\tLoading KJ07 models...."
        kj_p_vals, kj_pdot_vals, kj_dists = beammodels.load_kj2007_models()
        print "\tDone\n"

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
        if widthModel is None:
            print "\t\tWidth {0}% ".format(duty)
        else:
            print "\t\tUsing Karastergiou & Johnston beam width model"

        # set up progress bar for fun :)
        prog = ProgressBar(min_value=0,
                           max_value=ngen,
                           width=65,
                           mode='dynamic')

    # create survey objects here and put them in a list
    if surveyList is not None:
        surveys = [Survey(s) for s in surveyList]
    else:
        # make an empty list here - makes some code just a little
        # simpler - can still loop over an empty list (ie zero times)
        surveys = []

    # initialise these counters to zero
    for surv in surveys:
        surv.ndet = 0  # number detected
        surv.nout = 0  # number outside survey region
        surv.nsmear = 0  # number smeared out
        surv.ntf = 0  # number too faint

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
        alignpulsar(pulsar, pop)

        # braking index
        if pop.braking_index == 0:
            pulsar.braking_index = 2.5 + 0.5 * random.random()
        else:
            pulsar.braking_index = float(pop.braking_index)

        # apply relevant spin down model
        pulsar.dead = False  # pulsar should start alive!

        if pop.spinModel == 'fk06':
            spindown_fk06(pulsar)

            # apply deathline if relevant
            if pop.deathline:
                bhattacharya_deathperiod_92(pulsar)

        elif pop.spinModel == 'cs06':
            # contopoulos and spitkovsky
            spindown_cs06(pulsar, pop)

        # if period > 10 seconds, just try a new one
        if pulsar.period > 10000.0 or pulsar.period < 10.:
            continue
        # cut on pdot too - this doesn't help
        if pulsar.pdot > 1.e-11 or pulsar.pdot < 1.e-18:
            continue

        # define pulse width (default = 6% = 18 degrees)
        if widthModel is None:
            width = (float(duty)/100.) * pulsar.period**0.9
            width = math.log10(width)
            width = dists.drawlnorm(width, 0.3)
        elif widthModel == 'kj07':
            # Karastergiou & Johnston beam model

            # find closest p, pdot in the kj lists
            logp = math.log10(pulsar.period)
            logpdot = math.log10(pulsar.pdot)

            p_idx = (np.abs(kj_p_vals - logp)).argmin()
            pd_idx = (np.abs(kj_pdot_vals - logpdot)).argmin()

            # pick a width from relevant model
            dist = kj_dists[p_idx][pd_idx][2]
            width = np.random.choice(dist)

            """
            width = beammodels.kj2007_width(pulsar)
            no_width = 0
            while width  == 0.0:
                no_width+=1
                width = beammodels.kj2007_width(pulsar)

                # if we get to five, then try a new pulsar
                if no_width == 5:
                    no_width = 0
                    continue
            """

        else:
            print "Undefined width model!"
            sys.exit()
        # print width
        # print pulsar.period, width, pulsar.pdot
        if width == 0.0:
            # some kj2007 models make many zero-width sources. Skip!
            continue
        pulsar.width_degree = 360. * width / pulsar.period

        # plough on - only if the pulsar isn't dead!
        if not pulsar.dead or keepdead:
            # is the pulsar beaming?
            pulsar_beaming(pulsar, pop.beamModel)
            # if not, then skip onto next pulsar
            if not pulsar.beaming:
                continue

            # position of the pulsar
            galacticDistribute(pulsar, pop)
            # birthvelocity
            birthVelocity(pulsar, pop)

            # model the xyz velocity
            go.vxyz(pulsar)

            # luminosity
            if lumDistType == 'fk06':
                luminosity_fk06(pulsar,
                                alpha=pop.lumPar1,
                                beta=pop.lumPar2,
                                gamma=pop.lumPar3)

            elif lumDistType == 'lnorm':
                pulsar.lum_1400 = dists.drawlnorm(pop.lumPar1, pop.lumPar2)

            elif lumDistType == 'pow':
                pulsar.lum_1400 = dists.powerlaw(pop.lummin,
                                                 pop.lummax,
                                                 pop.lumpow)
            else:
                # something's wrong!
                raise EvolveException('Invalid luminosity distn selected')

            # apply efficiency cutoff
            if efficiencycut is not None:
                if pulsar.efficiency() > efficiencycut:
                    pulsar.dead = True
                    if not keepdead:
                        continue

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
            else:
                raise EvolveException('Invalid electron dist model selected')

            pulsar.scindex = scindex
            pulsar.t_scatter = go.scatter_bhat(pulsar.dm,
                                               pulsar.scindex)

            # if surveys are given, check if pulsar detected or not
            # in ANY of the surveys
            if surveyList is not None:
                detect_int = 0  # just a flag if pulsar is detected
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
                        # pulsar is just too faint
                        surv.ntf += 1
                        continue

                # if detected, increment ndet (for whole population)
                # and redraw the progress bar
                if detect_int:
                    pop.ndet += 1
                    # update the counter
                    if not nostdout:
                        prog.increment_amount()
                        print prog, '\r',
                        sys.stdout.flush()

            else:
                # no survey list, just add the pulsar to population,
                # and increment number of pulsars
                pop.ndet += 1
                # update the counter
                if not nostdout:
                    prog.increment_amount()
                    print prog, '\r',
                    sys.stdout.flush()

            # pulsar isn't dead, add to population!
            pop.population.append(pulsar)

        else:
            # pulsar is dead. If no survey list,
            # just increment number of pulsars
            if surveyList is None:
                pop.ndet += 1
                # update the counter
                if not nostdout:
                    prog.increment_amount()
                    print prog, '\r',
                    sys.stdout.flush()

    if not nostdout:
        print "\n\n"
        print "  Total pulsars = {0}".format(len(pop.population))
        print "  Total detected = {0}".format(pop.ndet)

        for surv in surveys:
            print "\n  Results for survey '{0}'".format(surv.surveyName)
            print "    Number detected = {0}".format(surv.ndet)
            print "    Number too faint = {0}".format(surv.ntf)
            print "    Number smeared = {0}".format(surv.nsmear)
            print "    Number outside survey area = {0}".format(surv.nout)

    # save list of arguments into the pop
    try:
        argspec = inspect.getargspec(generate)
        key_values = [(arg, locals()[arg]) for arg in argspec.args]
        pop.arguments = {key: value for (key, value) in key_values}
    except SyntaxError:
        pass

    return pop


def birthVelocity(pulsar, pop):
    """ Get a birth veolocity for the pulsar"""

    bvM = pop.birthVModel
    mean = pop.birthvmean
    sigma = pop.birthvsigma
    if bvM == 'gaussian':
        pulsar.vx = random.gauss(mean, sigma)
        pulsar.vy = random.gauss(mean, sigma)
        pulsar.vz = random.gauss(mean, sigma)
    elif bvM == 'exp':
        pulsar.vx = dists.draw_double_sided_exp(sigma, origin=mean)
        pulsar.vy = dists.draw_double_sided_exp(sigma, origin=mean)
        pulsar.vz = dists.draw_double_sided_exp(sigma, origin=mean)
    else:
        raise EvolveException('Invalid velocity model selected')


def galacticDistribute(pulsar, pop):
    """ select a galactic position - spiral arms on or off?"""
    if not pop.nospiralarms:
        # using spiral arms
        r0 = go.ykr()
        x, y = go.spiralize(r0)
    else:
        # distribute randomly in x-y plane
        x = -20. + random.random()*40.
        y = -20. + random.random()*40.

    # calculate z and r0
    z = dists.draw_double_sided_exp(pop.zscale)
    pulsar.galCoords = (x, y, z)
    pulsar.r0 = math.sqrt(x**2 + y**2)


def alignpulsar(pulsar, pop):
    """
    Pick an alignment angle for pulsar, depending on model chosen

    """

    # need to be careful to use chi in degrees, but
    # consistently remember to take the sins/cosines of radians
    if pop.alignModel == 'orthogonal':
        pulsar.chi = 90.0
        pulsar.sinchi_init = 1.0
        pulsar.sinchi = pulsar.sinchi_init
        pulsar.coschi = 0.

    elif pop.alignModel == 'random':
        chi = math.acos(random.random())  # in radians

        pulsar.chi = math.degrees(chi)  # -> degrees
        pulsar.sinchi_init = math.sin(math.radians(chi))
        pulsar.sinchi = pulsar.sinchi_init
        pulsar.coschi = math.cos(math.radians(chi))

    elif pop.alignModel == 'rand45':
        pulsar.coschi = random.random() * (1.0 - math.sqrt(0.5)) \
                            + math.sqrt(0.5)
        pulsar.chi = math.degrees(math.acos(pulsar.coschi))
        pulsar.sinchi_init = math.sin(math.radians(pulsar.chi))
        pulsar.sinchi = pulsar.sinchi_init

    elif pop.alignModel == 'wj08':
        # check an alignment timescale has been provided
        if pop.alignTime is None:
            raise EvolveException('Align timescale needed for WJ08 model')

        # weltevrede & johnston alignment model
        pulsar.coschi = random.random()
        chi = math.degrees(math.acos(pulsar.coschi))
        pulsar.sinchi_init = math.sin(math.radians(chi))
        pulsar.sinchi = pulsar.sinchi_init * math.exp(
            -pulsar.age/pop.alignTime)
        pulsar.chi = math.degrees(math.asin(pulsar.sinchi))

    else:
        raise EvolveException('Invalid alignment model: {0}'.format(aM))

    # more models to add here, but it'll do for now


def pulsar_beaming(pulsar, bM):
    """
    Work out if the pulsar is beaming --- model-dependent
    """

    # Tauris & Manchester beaming model (default)
    if bM == 'tm98':
        fraction = beammodels.tm98_fraction(pulsar)
    # more models to add here!
    elif bM == 'none':
        fraction = 1.0
    elif bM == 'const':
        fraction = 0.2
    elif bM == 'wj08':
        fraction = beammodels.wj08_fraction(pulsar)
    else:
        raise EvolveException('Invalid beaming model: {0}'.format(bM))

    if random.random() < fraction:
        pulsar.beaming = True
    else:
        pulsar.beaming = False


def luminosity_fk06(pulsar, alpha=-1.5, beta=0.5, gamma=0.18):
    """ Equation 14 from  Ridley & Lorimer """
    # variables to use in the equation
    delta_l = random.gauss(0.0, 0.8)

    # the equation
    logL = math.log10(gamma) + alpha*math.log10(pulsar.period/1000.) + \
        beta * math.log10(pulsar.pdot * 1.0e15) + delta_l

    # set L
    pulsar.lum_1400 = 10.0**logL


def spindown_fk06(pulsar):
    """
    Spindown model from Faucher-Giguere & Kaspi 2006
    """

    # k is a constant from the FK06 paper.
    # defined in cgs as
    # k = (8 * pi**2 * R**6)/(3 * I * c**3)
    k = 9.768E-40  # in cgs
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
    age_s = pulsar.age * 3.15569e7  # convert to seconds

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


def spindown_cs06(pulsar, pop):
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
    # converting the qsimp (fortran) to using scipy.integrate package
    # tested and they should give same output
    lower_limit = pulsar.p0*1000.
    upper_limit = 1.0E5
    const_of_integration = pulsar.coschi**2.0 / pdeath

    min_value = 1.0E24
    min_p = pulsar.p0*1000.

    index = pulsar.braking_index - 3.0
    temp_const = 3.3E-40*pulsar.bfield_init**2 * (pulsar.p0*1000.)**index \
                        * pulsar.age * 365.25 * 24. * 3.6E9  # 3.15569e7

    # do the integral
    result = scipy.integrate.quad(_cs06_poft,
                                  lower_limit,
                                  pdeath,
                                  args=(const_of_integration,
                                        pulsar.braking_index)
                                  )[0]

    if result < temp_const or result > 1.0E14:
        pulsar.period = 1.0E6
    else:
        looparray = np.arange(lower_limit, upper_limit+1)
        count = 0
        for m in looparray:
            count += 1
            result = scipy.integrate.quad(_cs06_poft,
                                          lower_limit,
                                          pdeath,
                                          args=(const_of_integration,
                                                pulsar.braking_index)
                                          )[0]
            if result > 1.0E14:
                pulsar.period = 1.0E6
                break

            tempmin = math.fabs(result - temp_const)

            if tempmin <= min_value:
                min_value = tempmin
                min_p = m
            else:
                pulsar.period = min_p
                break

            # end of loop
        else:
            pulsar.period = 1.0E6

    # see whether the pulsar should be dead
    # (are we using deathline? Is pulsar above it?)
    if pulsar.period > pdeath and pop.deathline:
        pulsar.dead = True
    else:
        pulsar.dead = False
        index = 2.0-pulsar.braking_index
        pulsar.pdot = 3.3E-40 * pulsar.bfield_init**2.0 * (1./pulsar.p0)\
            * (pulsar.period / (pulsar.p0*1000.))**index * \
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
    B = 3.2E19 * math.sqrt(pulsar.period * pulsar.pdot / 1000.)

    if B/(pulsar.period/1000.)**2 < 0.17E12:
        pulsar.dead = True


if __name__ == '__main__':
    """ 'Main' function; read in options, then generate population"""

    # set defaults here
    parser = argparse.ArgumentParser(
        description='Generate a population of pulsars')
    # number of pulsars to detect!
    parser.add_argument('-n', type=int, required=True,
                        help='number of pulsars to generate/detect')

    # list of surveys to use (if any)
    parser.add_argument('-surveys', metavar='S', nargs='+', default=None,
                        help='surveys to use to check if pulsars are detected'
                        )

    # maximum initial age of pulsars
    parser.add_argument('-tmax', type=float, required=False,
                        default=1.0E9,
                        help='maximum initial age of pulsars, in years (def=1.0E9)')

    # period distribution
    parser.add_argument('-p', nargs=2, type=float,
                        required=False, default=[0.3, 0.15],
                        help='mean and std dev of period distribution, secs \
                                 (def = 0.3, 0.15)')

    # luminosity distribution parameters
    parser.add_argument('-ldist', nargs=1, required=False, default=['fk06'],
                        help='distribution to use for luminosities',
                        choices=['fk06', 'lnorm', 'pow'])

    parser.add_argument('-l', nargs='+', required=False, type=float,
                        default=[-1.5, 0.5],
                        help='luminosity distribution parameters \
                             (def = [-1.5, 0.5], Faucher-Giguere&Kaspi, 2006)')

    # B field distribution
    parser.add_argument('-b', nargs=2, type=float,
                        required=False, default=[12.65, 0.55],
                        help='mean and std dev of log normal B field distn \
                         (Gauss, def = 12.65, 0.55)')

    # velocity distn model
    parser.add_argument('-vmodel', type=str, required=False,
                        nargs=1, default=['gaussian'],
                        choices=['gaussian', 'exp'],
                        help='velocity model to use')

    # velocity disttribution values
    parser.add_argument('-v', nargs=2, type=float,
                        required=False, default=[0.0, 180.],
                        help='velocity distn values\
                        (def = 0, 180 km/s)')

    # spectral index distribution
    parser.add_argument('-si', nargs=2, type=float,
                        required=False, default=[-1.4, 0.96],
                        help='mean and std dev of spectral index distribution \
                                 (def = -1.4, 0.96)')

    # spindown model
    parser.add_argument('-spinmodel', type=str, required=False,
                        nargs=1, default=['fk06'],
                        choices=['fk06', 'cs06'],
                        help='spin-down model to employ')

    # alignment model
    parser.add_argument('-alignmodel', type=str, required=False,
                        nargs=1, default=['orthogonal'],
                        choices=['orthogonal', 'random', 'rand45', 'wj08'],
                        help='pulsar alignment model to use')

    # alignment timescale
    parser.add_argument('-aligntime', type=float, required=False,
                        default=None,
                        help='alignment timescale')

    # beaming model
    parser.add_argument('-beammodel', type=str, required=False,
                        nargs=1, default=['tm98'],
                        choices=['tm98', 'none', 'const', 'wj08'],
                        help='beaming model to use (def=tm98)')

    # pulse width
    parser.add_argument('-w', type=float, required=False, default=5.,
                        help='pulse width, percent (def=5.0)')

    # pulse width model
    parser.add_argument('-wmod', type=str, required=False, default=None,
                        choices=[None, 'kj07'],
                        help='pulse width model to use (def=using -w opt')

    # scattering index
    parser.add_argument('-sc', type=float, required=False, default=-3.86,
                        help='modify the frequency-dependance of Bhat et al\
                                scattering formula (def = -3.86)')

    # efficiency cut off, if wanted
    parser.add_argument('-eff', type=float, required=False,
                        default=None,
                        help='efficiency cutoff value (def=None)')

    # galactic-Z distn
    parser.add_argument('-z', type=float, required=False, default=0.05,
                        help='exponential z-scale to use (def=0.05kpc)')

    # braking index
    parser.add_argument('-bi', type=float, required=False, default=0,
                        help='braking index value to use (def=0 = model)')

    # electron/dm model
    parser.add_argument('-dm', type=str, nargs=1, required=False,
                        default=['ne2001'],
                        help='Galactic electron distribution model to use',
                        choices=['ne2001', 'lmt85'])

    # output file name
    parser.add_argument('-o', type=str, metavar='outfile', required=False,
                        default='evolve.model',
                        help='Output filename for population model \
                               (def=evolve.model)')

    # turn off printing to stdout
    parser.add_argument('--nostdout', nargs='?', const=True, default=False,
                        help='switch off std output')

    # Flag for NOT using deathline
    parser.add_argument('--nodeathline', nargs='?', const=True, default=False,
                        help='turn OFF the deathline ')

    # flag to turn off spiral arms
    parser.add_argument('--nospiralarms', nargs='?', const=True, default=False,
                        help='turn off spiral arms galactic distn')

    parser.add_argument('--keepdead', nargs='?', const=True, default=False,
                        help='keep dead pulsars in the population model')

    args = parser.parse_args()

    # write command line to populate.cmd file
    with open('evolve.cmd', 'a') as f:
        f.write(' '.join(sys.argv))
        f.write('\n')

    # run the code!
    pop = generate(args.n,
                   surveyList=args.surveys,
                   age_max=args.tmax,
                   lumDistType=args.ldist[0],
                   lumDistPars=args.l,
                   pDistPars=args.p,
                   bFieldPars=args.b,
                   birthVModel=args.vmodel[0],
                   birthVPars=args.v,
                   alignModel=args.alignmodel[0],
                   alignTime=args.aligntime,
                   spinModel=args.spinmodel[0],
                   beamModel=args.beammodel[0],
                   siDistPars=args.si,
                   zscale=args.z,
                   nostdout=args.nostdout,
                   duty=args.w,
                   widthModel=args.wmod,
                   scindex=args.sc,
                   braking_index=args.bi,
                   electronModel=args.dm[0],
                   nodeathline=args.nodeathline,
                   efficiencycut=args.eff,
                   nospiralarms=args.nospiralarms,
                   keepdead=args.keepdead)

    pop.write(outf=args.o)
