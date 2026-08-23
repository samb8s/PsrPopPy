"""
Unit tests for dosurvey module
"""
import unittest
from os import path
import numpy as np
from psrpoppy.population import Population
from psrpoppy.pulsar import Pulsar
from psrpoppy import dosurvey
from psrpoppy import survey

class test_loadModel(unittest.TestCase):

    def test_unpickle_succeeds(self):
        testdir = path.join(path.dirname(path.dirname(dosurvey.__file__)), "tests")
        pop = dosurvey.loadModel(popfile=path.join(testdir, "populate.model"))
