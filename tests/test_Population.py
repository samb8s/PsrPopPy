import sys
import os
import unittest

from psrpoppy.population import Population


# not sure yet that population needs much testing... it's so simple
class testPopulation(unittest.TestCase):
    def setUp(self):
        self.population = Population()
