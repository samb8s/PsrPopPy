"""
Integration tests for dosurvey module
"""
import unittest
import pytest
import copy
from os import path
import shutil
import tempfile
import numpy as np
from psrpoppy.population import Population
from psrpoppy.pulsar import Pulsar
from psrpoppy import dosurvey
from psrpoppy import survey

class test_dosurvey_detection_scintillation(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

        self.pointing_file = path.join(self.tmpdir, "test_pointing.txt")
        self.survey_file = path.join(self.tmpdir, "test_survey.txt")
        self.gl = 0.0
        self.gb = 0.0

        # put beam at pulsar location so no random off-axis degradation
        with open(self.pointing_file, "w") as f:
            f.write("{0} {1}\n".format(self.gl, self.gb))

        survey_contents = """\
1.2  ! survey degradation factor
0.735  ! antenna gain (K/Jy)
4200  ! integration time (s)
0.125  ! sampling time (ms)
21.  ! system temperature (K)
1400  ! centre frequency (MHz)
256  ! bandwidth (MHz)
0.5  ! channel bandwidth (MHz)
2  ! number of polarizations
14.  ! full-width half maximum (arcmin)
0.  ! minimum RA (deg)
360.  ! maximum RA (deg)
-90.  ! minimum DEC (deg)
90.  ! maximum DEC (deg)
-180.  ! minimum Galactic longitude (deg)
180.  ! maximum Galactic longitude (deg)
0.  ! minimum abs(Galactic latitude) (deg)
90.  ! maximum abs(Galactic latitude) (deg)
1.0  ! fractional survey coverage (0-1)
8.0  ! signal-to-noise ratio
{0}  ! galactic pointing list
""".format(self.pointing_file)

        with open(self.survey_file, "w") as f:
            f.write(survey_contents)

        self.psr = Pulsar(
            period=10.0,
            dm=10.0,
            gl=self.gl,
            gb=self.gb,
            dtrue=1.0,
            lum_1400=1.0,
            spindex=-1.4,
            width_degree=18.0,
            scindex=-3.86,
        )

        # Pulsar.__init__ currently ignores the t_scatter argument,
        # so set it explicitly.  Use a very small non-zero value so
        # scintillation is in the weak regime and does not invoke NE2001.
        self.psr.t_scatter = 1.0e-8

        self.pop = Population(ref_freq=1400.0)
        self.pop.population.append(self.psr)

        # adjust luminosity so that the real survey calculation gives
        # an unscintillated S/N of exactly 7, just below SNRlimit=8.
        surv = survey.Survey(self.survey_file)
        snr = surv.SNRcalc(self.psr, self.pop)

        self.psr.lum_1400 *= 7.0 / snr

        # save the real random.gauss so tearDown can restore it
        self.real_gauss = survey.random.gauss

    def tearDown(self):
        survey.random.gauss = self.real_gauss
        shutil.rmtree(self.tmpdir)

    def test_unscintillated_snr(self):
        surv = survey.Survey(self.survey_file)
        snr_noscint = surv.SNRcalc(self.psr, self.pop)
        np.testing.assert_allclose(snr_noscint, 7.)
        
    def test_scintillation_changes_detection(self):
        # use independent pulsar/population objects bc dosurvey.run()
        # modifies pulsar.detected
        pop_noscint = copy.deepcopy(self.pop)
        psr_noscint = pop_noscint.population[0]

        result_noscint = dosurvey.run(pop_noscint,
                                      [self.survey_file],
                                      nostdout=True,
                                      scint=False)

        self.assertFalse(psr_noscint.detected)
        self.assertEqual(result_noscint[0][2].ndet, 0)

        # force the stochastic scintillation draw to be +1 sigma.
        survey.random.gauss = lambda mean, sigma: mean + sigma

        pop_scint = copy.deepcopy(self.pop)
        psr_scint = pop_scint.population[0]

        result_scint = dosurvey.run(pop_scint,
                                    [self.survey_file],
                                    nostdout=True,
                                    scint=True)

        self.assertTrue(psr_scint.detected)
        self.assertEqual(result_scint[0][2].ndet, 1)
