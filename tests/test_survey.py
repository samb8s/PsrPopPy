"""
Unit and integration tests for survey module
"""
import unittest
import tempfile
import os
from os import path
from psrpoppy.survey import Survey
import psrpoppy

class test_survey_pointing_list_IO(unittest.TestCase):
    """
    integration test for Survey reads pointing list correctly
    """
    def test_survey_reads_pointing_file_in_surveys_dir(self):
        surveydir = path.join(path.dirname(psrpoppy.__file__), "surveys")
        surv_fname = "testsurv"
        surv_ptfname = "testsurv.glgb"
        surv_path = path.join(surveydir, surv_fname)
        point_path = path.join(surveydir, surv_ptfname)
        self.assertFalse(path.exists(surv_path))
        self.assertFalse(path.exists(point_path))
        self.addCleanup(os.remove, surv_path)
        self.addCleanup(os.remove, point_path)
        surv_contents = (
"""1.2  ! survey degradation factor
0.6  ! antenna gain (K/Jy)
2100 ! integration time (s)
0.25 ! sampling time (ms)
25.  ! system temperature (K)
1374 ! centre frequency (MHz)
288  ! bandwidth (MHz)
3.0  ! channel bandwidth (MHz)
2    ! number of polarizations
14.  ! full-width half maximum (arcmin)
0.   ! minimum RA (deg)
360. ! maximum RA (deg)
-90.  ! minimum DEC (deg)
90.  ! maximum DEC (deg)
-100  ! minimum Galactic longitude (deg)
50.  ! maximum Galactic longitude (deg)
0.   ! minimum abs(Galactic latitude) (deg)
5.   ! maximum abs(Galactic latitude) (deg)
testsurv.glgb ! pointing list galactic
1.0  ! fractional survey coverage (0-1)
9.0  ! signal-to-noise ratio
""")
        with open(surv_path, "w") as f:
            f.write(surv_contents)
        with open(point_path, "w") as f:
            f.write("0.\t0.")
        surv = Survey("testsurv")
