#!/usr/bin/python

import os
import sys
import numpy as np

import math
import random
from scipy.special import j1

import galacticops as go
from population import Population

class SurveyException(Exception):
    pass

class CoordinateException(Exception):
    pass

class Pointing:
    """Simple class -- pointing has a gl and gb position (***NEWLY ADDED GAIN/TOBS***)"""
    """This class works with pointing list files, which have GL/GB/GAIN/TOBS columns. """

    def __init__(self, coord1, coord2, coordtype, gain, tobs):
        """Set the coordinates of the pointing.
            Convert to gl and gb if in RA and Dec"""

        if coordtype not in ['eq', 'gal']:
            raise CoordinateException('Wrong coordtype passed to Pointing')

        if coordtype == 'eq':
            # assume pointings in decimal degrees
            ra = coord1
            dec = coord2

            # convert to l and b :)
            gl,gb = go.radec_to_lb(ra, dec)

            if gl>180.:
                gl -= 360.

            self.gl = gl
            self.gb = gb
        else:
            if coord1>180.:
                coord1 -= 360.

            self.gl = coord1
            self.gb = coord2
            self.tobs = tobs 
            self.gain = gain

class Survey:
    """Class to store survey parameters and methods"""
    def __init__(self, surveyName, pattern='gaussian'):
        """Read in a survey file and obtain the survey parameters"""
        
        # try to open the survey file locally first
        if os.path.isfile(surveyName):
            f = open(surveyName, 'r')
        else:
            try:
                # try to open file in lib
                # get path to surveys directory
                __dir__ = os.path.dirname(os.path.abspath(__file__))
                __libdir__ = os.path.dirname(__dir__)
                filepath = os.path.join(__libdir__, 'surveys', surveyName)
                f = open(filepath, 'r')

            except IOError:
                # couldn't find the file
                s = 'File {0} does not exist!!!'.format(surveyName)
                raise SurveyException(s)

        self.surveyName = surveyName
        # initialise the pointings list to None
        # only change this is there is a list of pointings to be used
        self.pointingslist = None
        self.gainpat       = pattern 

        # adding AA parameter, so can scale s/n if the survey is
        # an aperture array
        self.AA = False

        # Parse the file line by line
        for line in f:
            # ignore any lines starting '#'
            if line.strip()[0] == '#':
                continue
            # otherwise, parse!
            a = line.split('!')
            
            # new feature - possible to have a list of positions in survey,
            # rather than a range of l,b or ra,dec
            if a[1].count('pointing list'):
                pointfname = a[0].strip()
                pointfpath = os.path.join(__libdir__, 'surveys', pointfname)

                # try to open the pointing list
                if os.path.isfile(pointfpath):
                    pointfptr = open(pointfpath, 'r')
                else:
                    try:
                        # try to open a local file
                        pointfptr = open(pointfname, 'r')
                    except:
                        s = 'File {0} does not exist!!!'.format(pointfpath)
                        raise CoordinateException(s)

                # read in the pointing list
                self.pointingslist = []
                if a[1].count('galactic'):
                    # already galactic coordinates
                    for line in pointfptr:
                        a = line.split()
                        # Expected form of pointing list files is hard-coded.
                        if len(a) != 4:
                            s = 'File {0} should have columns: gl/gb/gain/tobs!'.format(pointfpath)
                            raise CoordinateException(s)
                        p = Pointing(float(a[0]), float(a[1]), 'gal', float(a[2]), float(a[3]))
                        self.pointingslist.append(p)

                elif a[1].count('equatorial'):
                    # need to be converted from ra dec to gl gb
                    for line in pointfptr:
                        a = line.split()
                        if len(a) != 4:
                            s = 'File {0} should have columns: gl/gb/gain/tobs!'.format(pointfpath)
                            raise CoordinateException(s)
                        p = Pointing(float(a[0]), float(a[1]), 'eq', float(a[2]), float(a[3]))
                        self.pointingslist.append(p)

                else:
                    s = 'Coordinate type unspecified in {0}.'.format(surveyName)
                    raise CoordinateException(s)

                pointfptr.close()

            elif a[1].count('survey degradation'):
                # beta
                self.beta = float(a[0].strip())
            elif a[1].count('antenna gain'):
                # gain
                self.gain = float(a[0].strip())
            elif a[1].count('integration time'):
                # tobs
                self.tobs = float(a[0].strip())
            elif a[1].count('sampling'):
                # tsamp
                self.tsamp = float(a[0].strip())
            elif a[1].count('system temperature'):
                # tsys
                self.tsys = float(a[0].strip())
            elif a[1].count('centre frequency'):
                # centre frequency
                self.freq = float(a[0].strip())
            elif a[1].strip().startswith('bandwidth'):
                # bandwidth
                self.bw = float(a[0].strip())
            elif a[1].count('channel bandwidth'):
                # bw_chan
                self.bw_chan = float(a[0].strip())
            elif a[1].count('polarizations'):
                # num polns
                self.npol = float(a[0].strip())
            elif a[1].count('half maximum'):
                # FWHM
                self.fwhm = float(a[0].strip())
            elif a[1].count('minimum RA'):
                # min RA
                self.RAmin = float(a[0].strip())
            elif a[1].count('maximum RA'):
                # max RA
                self.RAmax = float(a[0].strip())
            elif a[1].count('minimum DEC'):
                # min dec
                self.DECmin = float(a[0].strip())
            elif a[1].count('maximum DEC'):
                # mac dec
                self.DECmax = float(a[0].strip())
            elif a[1].count('minimum Galactic'):
                # min longitude
                self.GLmin = float(a[0].strip())
            elif a[1].count('maximum Galactic'):
                # max longitude
                self.GLmax = float(a[0].strip())
            elif a[1].count('minimum abs'):
                # min latitude
                self.GBmin = float(a[0].strip())
            elif a[1].count('maximum abs'):
                # max latitude
                self.GBmax = float(a[0].strip())
            elif a[1].count('coverage'):
                # coverage fraction
                self.coverage = float(a[0].strip())
                if self.coverage > 1.0:
                    self.coverage = 1.0
            elif a[1].count('signal-to-noise'):
                # SNR limit
                self.SNRlimit = float(a[0].strip())
            elif a[1].count('gain pattern'):
                # Gain pattern (can specify airy, default = gaussian)
                self.gainpat = a[0].strip()
            elif a[1].count('Aperture Array'):
                # turn on AA
                self.AA = True
            else:
                print "Parameter '",a[1].strip(),"' not recognized!"

        f.close()

        # get tsky array from file
        self.tskylist = go.readtskyfile()

    def __str__(self):
        """Method to define how to print the class"""
        s = "Survey class for {0}:".format(self.surveyName)
        s = '\n\t'.join([s, "beta = {0}".format(self.beta)])
        s = '\n\t'.join([s, "gain = {0}".format(self.gain)])
        s = '\n\t'.join([s, "tobs = {0} s".format(self.tobs)])
        s = '\n\t'.join([s, "tsamp = {0} ms".format(self.tsamp)])
        s = '\n\t'.join([s, "Tsys = {0} K".format(self.tsys)])
        s = '\n\t'.join([s, "Centre frequency = {0} MHz".format(self.freq)])
        s = '\n\t'.join([s, "Bandwidth = {0} MHz".format(self.bw)])
        s = '\n\t'.join([s, "Chan BW = {0} MHz".format(self.bw_chan)])
        s = '\n\t'.join([s, "Num polarisations = {0}".format(self.npol)])
        s = '\n\t'.join([s, "FWHM = {0} arcmin".format(self.fwhm)])
        s = '\n\t'.join([s, "SNR limit = {0}".format(self.SNRlimit)])

        return s

    def nchans(self):
        """ Returns the number of channels in the survey backend."""
        return self.bw / self.bw_chan


    def inRegion(self, pulsar):
        """Test if pulsar is inside region bounded by survey."""
        # check if l, b are outside region first of all
        if pulsar.gl>180.:
            pulsar.gl -= 360.
        if pulsar.gl > self.GLmax or pulsar.gl < self.GLmin:
            return False
        if math.fabs(pulsar.gb) > self.GBmax or math.fabs(pulsar.gb) < self.GBmin:
            return False

        # need to compute ra/dec of pulsar from the l and b (galtfeq)
        ra, dec = go.lb_to_radec(pulsar.gl, pulsar.gb)

        # are ra, dec outside region?
        if ra > self.RAmax or ra < self.RAmin:
            return False
        if dec > self.DECmax or dec < self.DECmin:
            return False
        
        # randomly decide if pulsar is in completed area of survey
        if random.random() > self.coverage:
            return False
        
        return True 

    def inPointing(self, pulsar):
        """Calculate whether pulsar is inside FWHM/2 of pointing position.
        Currently breaks as soon as it finds a match. !!!Could be a closer
        position further down the list!!!"""
        # initialise offset_deg to be a big old number
        # FWHM is in arcmin so always multiply by 60
        offset_deg = 5.

        # loop over pointings
        for point in self.pointingslist:
            # do a really basic check first

            glterm = (pulsar.gl - point.gl)**2
            gbterm = (pulsar.gb - point.gb)**2
            offset_new = math.sqrt(glterm + gbterm)

            # if the beam is close enough, break out of the loop
            if offset_new < offset_deg:
                offset_deg = offset_new
                self.gain  = point.gain
                self.tobs  = point.tobs
                
        return offset_deg

    def SNRcalc(self, pulsar, pop):
        """Calculate the S/N ratio of a given pulsar in the survey"""
        # if not in region, S/N = 0

        # if we have a list of pointings, use this bit of code
        # haven't tested yet, but presumably a lot slower
        # (loops over the list of pointings....)
 
        if pulsar.dead:
            return 0.
        # otherwise check if pulsar is in entire region
        if self.inRegion(pulsar):
            # If pointing list is provided, check how close nearest 
            # pointing is
            if self.pointingslist is not None:
                # convert offset from degree to arcmin
                offset = self.inPointing(pulsar) * 60.0

            else:
                # calculate offset as a random offset within FWHM/2
                offset = self.fwhm * math.sqrt(random.random()) / 2.0
        else:
            return -2

        # Get degfac depending on self.gainpat
        if self.gainpat == 'airy':
            conv    = math.pi/(60*180.)         # Conversion arcmins -> radians
            eff_diam = 3.0e8/(self.freq*self.fwhm*conv*1.0e6)  # Also MHz -> Hz
            a       = eff_diam/2.               # Effective radius of telescope
            lamda   = 3.0e8/(self.freq*1.0e6)   # Obs. wavelength
            kasin   = (2*math.pi*a/lamda)*np.sin(offset*conv)
            degfac = 4*(j1(kasin)/kasin)**2
        else:
            degfac = math.exp(-2.7726 * offset * offset / (self.fwhm *self.fwhm))

        # calculate total temperature
        Ttot = self.tsys + self.tskypy(pulsar)

        # calc dispersion smearing across single channel
        tdm = self._dmsmear(pulsar)

        # calculate bhat et al scattering time (inherited from GalacticOps)
        # in units of ms
        tscat = go.scatter_bhat(pulsar.dm, pulsar.scindex, self.freq)

        # Calculate the effective width
        width_ms = pulsar.width_degree * pulsar.period / 360.0
        weff_ms = math.sqrt(width_ms**2 + self.tsamp**2 + tdm**2 + tscat**2)

        # calculate duty cycle (period is in ms)
        delt = weff_ms / pulsar.period
        #print weff_ms, pulsar.period

        # if pulse is smeared out, return -1.0
        if delt > 1.0:
            #print weff_ms, tscat, pulsar.dm, pulsar.gl, pulsar.gb, pulsar.dtrue
            return -1
        else:
            sig_to_noise = self._SNfac(pulsar, pop.ref_freq, degfac, Ttot) \
                                  * math.sqrt((1.0 -delt)/delt)

            if self.AA:
                sig_to_noise *= self._AA_factor(pulsar)
                
            return sig_to_noise

    def _AA_factor(self, pulsar):
        """ Aperture array factor """

        # need to compute ra/dec of pulsar from the l and b (galtfeq)
        ra, dec = go.lb_to_radec(pulsar.gl, pulsar.gb)

        offset_from_zenith = dec - (self.DECmax + self.DECmin)/2.0

        return math.cos(math.radians(offset_from_zenith))

    def _SNfac(self, psr, ref_freq, degfac, Ttot):
        """The S/N factor from system parameters"""
        # scale flux to survey frequency
        if psr.gpsFlag == 1:
            # do crazy flux calculation for GPS sources
            flux = self._gpsFlux(psr, ref_freq)

        elif psr.brokenFlag == 1 and self.freq < ref_freq:
            # assuming the alpha_1 value is for freq<ref_freq (which is ~1GHz)
            flux = psr.s_1400() * (self.freq / ref_freq)**psr.brokenSI
        
        else:
            flux = psr.s_1400() * (self.freq / ref_freq)**psr.spindex

        return flux * degfac * self.gain * \
                  math.sqrt(self.npol * self.bw * self.tobs) \
                  / self.beta / Ttot

    def _gpsFlux(self, psr, ref_freq):
        """Calculate the flux assuming GPS spectrum shape, spindex===b"""
        log_nu_1 = math.log10(ref_freq/1000.)
        log_nu_2 = math.log10(self.freq/1000.)
        gpsC = math.log10(psr.s_1400()) - (psr.gpsA * log_nu_1**2) \
                                            - psr.spindex * log_nu_1
        return 10.**(psr.gpsA * log_nu_2**2 + psr.spindex * log_nu_2 + gpsC)

    def _dmsmear(self, psr):
        """Calculate the smearing across a channel due to the pulsar DM"""
        return 8.3E6 * psr.dm * self.bw_chan / math.pow(self.freq, 3.0)

    def tskypy(self, psr):
        """ Calculate tsky from Haslam table, scale to survey frequency"""
        # ensure l is in range 0 -> 360
        b = psr.gb
        if psr.gl < 0.:
            l = 360 + psr.gl
        else:
            l = psr.gl

        # convert from l and b to list indices
        j = b + 90.5
        if j > 179:
            j = 179

        nl = l - 0.5
        if l < 0.5:
            nl = 359
        i = float(nl) / 4.
        
        tsky_haslam = self.tskylist[180*int(i) + int(j)]
        # scale temperature before returning
        return tsky_haslam * (self.freq/408.0)**(-2.6)

    def scint(self, psr, snr):
        """ Add scintillation effects and modify the pulsar's S/N"""

        # calculate the scintillation strength (commonly "u")
        # first, calculate scint BW, assume Kolmogorov, C=1.16
        tscat = go.scatter_bhat(psr.dm, psr.scindex, self.freq) # ms
        tscat /= 1000. # convert to seconds

        scint_bandwidth = 1.16 / 2.0 / math.pi / tscat # BW in Hz
        scint_bandwidth /= 1.0E6 # convert to MHz (self.freq is in MHz)

        scint_strength = math.sqrt(self.freq / scint_bandwidth)

        if scint_strength < 1.0:
            # weak scintillation
            # modulation index
            u_term = math.pow(scint_strength, 1.666666)
            mod_indx = math.sqrt(u_term)

        else:
            # strong scintillation

            # m^2 = m_riss^2 + m_diss^2 + m_riss * m_diss
            # e.g. Lorimer and Kramer ~eq 4.44 
            m_riss = math.pow(scint_strength, -0.33333)

            # lorimer & kramer eq 4.44
            kappa = 0.15 # taking this as avrg for now
            
            # calculate scintillation timescale
            scint_ts, scint_bw = go.ne2001_scint_time_bw(psr.dtrue, 
                                                         psr.gl, 
                                                         psr.gb, 
                                                         self.freq)
            
            # calc n_t and n_f
            if scint_ts is None:
                n_t = 1.
            else:
                n_t = self._calc_n_t(kappa, scint_ts)
                
            if scint_bw is None:
                n_f = 1.
            else:
                n_f = self._calc_n_f(kappa, scint_bw)
            
            # finally calc m_diss
            m_diss = 1. / math.sqrt(n_t * n_f)

            m_tot_sq = m_diss * m_diss + m_riss * m_riss + m_riss * m_diss
            
            # modulation index for strong scintillation
            mod_indx = math.sqrt(m_tot_sq)
      
        return self._modulate_flux_scint(snr, mod_indx)

    def _calc_n_t(self, kappa, delt_t):
        """Number of scintles sampled in time"""
        return 1. + kappa * self.tobs / delt_t
        
    def _calc_n_f(self, kappa, delt_f):
        """Number of scintles sampled in frequency"""
        return 1. + kappa * self.bw / delt_f

    def _modulate_flux_scint(self, snr, mod_indx):
        """Modify pulsar flux (actually S/N) according to the modulation index"""
        # flux and S/N are obviously proportional so it's simple to do this 
        # sigma of scintillation
        sig_scint = mod_indx * snr
        return random.gauss(snr, sig_scint)
