#!/usr/bin/python

import math

def calc_gamma1(pulsar, survey):
    """Calculate the degradation of SNR For a pulsar in a survey with no
        acceleration searching (Bagchi, Lorimer & Wolfe 2013"""


    # angles from 0 -> 2pi
    full_circle = [x*10.*3.142/180. for x in range(36)]
    
    # calculate probabilities, and normalise
    probabilities = [_calc_prob(pulsar, x) for x in full_circle]
    probabilities = _normalise(probabilities)

    # calculate keplerian parameters
    ecc_anomaly = [_calc_ecc_anomaly(pulsar, x) for x in full_circle]
    mean_anomaly = [_calc_mean_anomaly(pulsar, x) for x in ecc_anomaly]
    peri_epoch = [-x/pulsar.omega_orbit for x in mean_anomaly]


    # calculate alpha and gamma


def find_alpha1(pulsar, survey):
    ### THIS BIT NEEDS FILLING IN

    

    return False


def kepler_solve_1(ecc, mean_anom):

    ecc_anom = mean_anom + ecc * math.sin(mean_anom) * \ 
                     (1. + ecc * math.cos(mean_anom))
    e_fit = ecc_anom + 10.

    # newton-raphson method
    n_inter = 0 
    while math.fabs(e_fit - ecc_anom)>1.E-10:
        
        e_fit = _newton_raphson_kepler1(mean_anom, ecc_anom, ecc)
        ecc_anom = e_fit
        n_iter += 1
        if n_iter > 100:
            break

    ecc_anom = _newton_raphson_kepler1(mean_anom, ecc_anom, acc)

    return ecc_anom


def _newton_raphson_kepler1(ma, ea, ecc):
    return ma+ecc * math.sin(ea) - ea*ecc*math.cos(ea)/(1. - ecc*math.cos(ea))

def kepler_solve_2(eyy, Etyy):

    a1 = math.sqrt((1.0 + eyy)/(1.0 - eyy)) * math.tan(Etyy/2.0)
    return 2.0 * math.atan(a1)

def _normalise(l):
    """Normalise a list l"""
    m = max(l)
    return [x / m for x in l]

def _calc_prob(pulsar, angle_rad):
    
    denom = (1. + pulsar.ecc * math.cos(angle))**2.
    return 1./denom

def _calc_ecc_anomaly(pulsar, angle_rad):

    sqrt_term = math.sqrt((1. - pulsar.ecc) / (1. + pulsar.ecc))
    e0ar = 2.0 * math.atan(sqrt_term * math.tan(angle_rad/2.0) )
    return e0ar

def _calc_mean_anomaly(pulsar, ecc_anom):

    return ecc_anom - pulsar.ecc * math.sin(ecc_anom)
