#!/usr/bin/python

import sys
import math
import random
from os import path
import numpy as np
from scipy.special import expit

__dir__ = path.dirname(path.abspath(__file__))

def drawlnorm(mean, sigma):
    """Draw a random number from a log-normal distribution"""

    return 10.0**random.gauss(mean, sigma)


def powerlaw(minval, maxval, power):
    """Draw a value randomly from the specified power law"""

    logmin = math.log10(minval)
    logmax = math.log10(maxval)

    c = -1.0 * logmax * power
    nmax = 10.0**(power*logmin + c)

    # slightly worried about inf loops here...
    while True:
        log = random.uniform(logmin, logmax)
        n = 10.0**(power*log + c)

        if nmax*random.random() <= n:
            break

    return 10.0**log


def draw1d(dist):
    """Draw a bin number form a home-made distribution
        (dist is a list of numbers per bin)
    """
    # sum of distribution
    total = sum(dist)
    # cumulative distn
    cumulative = [sum(dist[:x+1])/total for x in range(len(dist))]

    rand_num = random.random()
    for i, c in enumerate(cumulative):
        if rand_num <= c:
            return i


def draw_double_sided_exp(scale, origin=0.0):
    """Exponential distribution around origin, with scale height scale."""
    if scale == 0.0:
        return origin

    rn = random.random()
    sign = random.choice([-1.0, 1.0])

    return origin + sign * scale * math.log(rn)

def ppdot_cohen27_posterior_draws(
        fpath="correlated_invWeibull-P_lognormal-Pdot_multivarnorm_approx.npz",
        size=1,
        rng=None
):
    """
    Return inverse-transformed joint population parameter posterior draws
    for multivariate normal approximation to correlated inverse-Weibull-period,
    lognormal-Pdot joined by a Gaussian copula from Cohen et al. (2027, in prep)

    Parameters:
    ----------
    fpath : (str, optional) numpy zip file containing numpy array with
    transformed parameter means, covariance matrix, and names. Local file
    with same name overrides default.
    size : (int, optional) number of samples for each parameter
    rng : (numpy.random.RandomState, optional) random state

    Returns:
    -------
    draws : (numpy.ndarray) (5,size) array of param draws in order of names
    names : (numpy.ndarray) (5,) array of param names
    """
    if path.isfile(fpath):
        dat = np.load(fpath)
    elif path.isfile(path.join(__dir__, fpath)):
        dat = np.load(path.join(__dir__,
                                fpath))
    else:
        raise IOError(2, "'fpath' does not exist. ",
                      fpath)
    
    if rng is None:
        rng = np.random.RandomState()

    mean = dat["mean"]
    covariance = dat["covariance"]
    names = dat["parameter_names"]
        
    transformed = rng.multivariate_normal(mean,
                                          covariance,
                                          size=size,
                                          check_valid="raise")
    samples = np.column_stack([
        0.1 + (3.0 - 0.1) * expit(transformed[:, 0]),
        0.1 + (8.0 - 0.1) * expit(transformed[:, 1]),
        transformed[:, 2],
        np.exp(transformed[:, 3]),
        np.tanh(transformed[:, 4] / 2.0),
    ])
    draws = samples.T

    return draws, names
