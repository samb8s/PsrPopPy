#!/usr/bin/python

import math
import random

class Beaming:
    """
    Various beaming models
    """

    def __init__(self):
        pass

    def tm98_fraction(self, period_s):
        """
        Tauris and Manchester 1998 beaming fraction
        """
        periodterm = math.log10(period_s) - 1.0

        return 0.03 + 0.09 * periodterm**2.0
