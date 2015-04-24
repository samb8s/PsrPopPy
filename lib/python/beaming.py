#!/usr/bin/python

import os
import sys
import math
import random

from operator import itemgetter

import numpy as np


def tm98_fraction(pulsar):
    """
    Tauris and Manchester 1998 beaming fraction
    """
    periodterm = math.log10(pulsar.period / 1000.) - 1.0

    return 0.03 + 0.09 * periodterm**2.0


def wj08_fraction(pulsar):
    """
    Weltevrede & Johnston 2008 beaming model
    """
    if pulsar.chi < 1.0E-5:
        pulsar.chi = 0.

    rho = 5.4 * (pulsar.period/1000.)**(-0.5)
    beta_temp = math.degrees(math.acos(random.random()))
    if beta_temp <= rho:
        fraction = 1.
    else:
        fraction = 0.

    # get rho and chi in radians for simplicity
    chi_rad = math.radians(pulsar.chi)
    rho_rad = math.radians(rho)

    if pulsar.chi > rho and (pulsar.chi + rho) < 90.:
        fraction = 2.0 * math.sin(chi_rad) * math.sin(rho_rad)
    elif pulsar.chi > rho and (pulsar.chi + rho) > 90.:
        fraction = math.cos(chi_rad - rho_rad)
    elif pulsar.chi <= rho and (pulsar.chi + rho) < 90.:
        fraction = 1.0 - math.cos(chi_rad + rho_rad)
    else:
        fraction = 1.

    return fraction


def load_kj2007_models():
    """
    Karastergiou & Johnston 2007 beam model -
    Load the pre-run reesults
    """

    # get current directory
    __dir__ = os.path.dirname(os.path.abspath(__file__))
    __models__ = os.path.join(__dir__, 'models')

    # lists to save p, pdot values in
    pd_vals = []
    p_vals = []
    all_dists = []

    # get dist_XXX_YYY.npy file names
    for model in os.listdir(__models__):
        # get full path and root of filename
        filepath = os.path.join(__models__, model)
        filename = os.path.splitext(model)[0]

        # get p, pdot vals
        vals = filename.split("_")
        p, pd = float(vals[1]), float(vals[2])
        p_vals.append(p)
        pd_vals.append(pd)

        dist = np.load(filepath)
        all_dists.append((p, pd, dist))

    # now get unique values of p, pdot
    p_vals = sorted(set(p_vals))
    pd_vals = sorted(set(pd_vals), reverse=True)

    # and sort the dists accordingly
    sorted_dists = []
    for p in p_vals:
        dists = []
        # get dists with correct p value
        for dist in all_dists:
            if dist[0] == p:
                dists.append(dist)

        # sort according to pdot value, i.e. dist[1]
        sorted(dists, key=itemgetter(1))

        # add dists to sorted_dists
        sorted_dists.append(dists)

    return np.array(p_vals), np.array(pd_vals), sorted_dists


def kj2007_width(pulsar):
    """
    Karastergiou & Johnston 2007 beam model - get width from pulsar
    parameters
    """

    # get period and edot in correct format
    period = pulsar.period / 1000.
    edot = math.log10(pulsar.edot())

    if edot < 35.:
        hmin = 20.
        hmax = 1000.
        ncomp = 4
        npatch = 5
    else:
        hmin = 950.
        hmax = 1000.
        ncomp = 1
        npatch = 10

    rhomax = np.degrees(3. * np.sqrt(np.pi * hmax /
                        (2. * 3.e5 * period)))
    beta = 1000.

    while beta > rhomax:
        alpha = math.degrees(math.acos(random.random()))
        zeta = math.degrees(math.acos(random.random()))
        beta = zeta - alpha

    # if we got here, everything's good, I guess
    alpha = alpha
    beta = beta
    zeta = beta + alpha
    zeta = zeta

    hrange = hmax - hmin
    hcomp = hmin + hrange * np.random.random(ncomp)**2
    # NEW HCOMP??? (only freq dependent, so din't include for now)

    rhocomp = np.degrees(3. * np.sqrt(np.pi * hcomp /
                                      (2. * 3.e5 * period)))
    pwcomp = 0.49 * np.sqrt(hcomp / 10. / period)

    delta = 360. / 1000.
    xarray = np.arange(-180., +180., delta)
    yarray = np.arange(-180., +180., delta)

    stokes_i = np.zeros((len(xarray), len(yarray)))

    for rho, pw in zip(rhocomp, pwcomp):
        stokes_i = patchbeam(stokes_i,
                             xarray,
                             yarray,
                             rho,
                             pw,
                             npatch)

    # compute max profile width
    rhomax = rhocomp.max()
    pwmax = pwcomp.max()
    wtmp = sin_deg(rhomax/2. + pwmax/2.)**2
    wtmp -= sin_deg(beta)**2

    """
    if wtmp < 0 or alpha + beta == 0.:
        widthmax = 360.
    else:
        albe = alpha + beta
        widthmax = math.sqrt(wtmp /
                               sin_deg(alpha) /
                               math.fabs(sin_deg(albe)))
        if widthmax > 1.:
            widthmax = 360.
        else:
            widthmax = 4. * math.asin(widthmax)
            widthmax = math.degrees(widthmax)
    """
    # compute line of sight angle
    theta_los, xlos, ylos = get_lineofsight(alpha, beta)

    # compute rvm swing
    # pa = rvm()

    # get the line-of-sight values for the pulse profile
    # I'll just consider stokes-I for now
    xinds = np.rint((xlos + 180.) / delta).astype(int)
    yinds = np.rint((ylos + 180.) / delta).astype(int)

    # get profile
    prof = stokes_i[xinds, yinds]

    # need to make sure this function returns width in
    # whatever units are consistent with my other code
    return calcwidth(prof) * pulsar.period


def get_stokes_index(val, x1, dx):
    n = round((val - x1) / dx)
    return int(n)


def get_lineofsight(alpha, beta):
    delta = 360./1024.
    phiarray = np.arange(-180., +180., delta)

    xlos, ylos, thetalos = [], [], []
    for phi in phiarray:
        xp, yp = mapphi(alpha, beta, phi)
        th = math.degrees(math.atan2(xp, yp)) - 90.
        thetalos.append(th)
        xlos.append(xp)
        ylos.append(yp)

    return thetalos, np.array(xlos), np.array(ylos)


def mapphi(alpha, beta, phi):
    cosR = cos_deg(alpha + beta) \
            * cos_deg(alpha) \
            + sin_deg(alpha + beta)\
            * sin_deg(alpha) \
            * cos_deg(phi)

    cosR = correct(cosR)

    R = math.degrees(math.acos(cosR))
    if int(R * 100.0) == 18000:
        R = int(R*100.0)/100.0

    if (R != 0.0) and (R != 180.) and (alpha > 0.):
        cosgamma = (cos_deg(alpha + beta) - cos_deg(alpha) * cosR)
        cosgamma /= (sin_deg(alpha) * sin_deg(R))
    else:
        cosgamma = 0.0

    cosgamma = correct(cosgamma)
    gamma = math.degrees(math.acos(cosgamma))
    xp = R * sin_deg(gamma)
    if phi > 0:
        xp = 0. - xp

    yp = -R * cosgamma

    return xp, yp


def rvm(alpha, zeta):

    delta = 360./1024.
    xprof = np.arange(-180., +180., delta)
    phi0, psi0 = 0.0, 0.0

    pa = [rvm_calc(x, phi0, psi0, alpha, zeta) for x in xprof]

    return pa


def rvm_calc(x, phi0, psi0, alpha, zeta):
    """Do the rvm calculation"""
    phi = x * math.pi / 180.

    num = sin_deg(alpha) * \
        sin_deg(phi - phi0)
    denom = sin_deg(zeta) * cos_deg(alpha) - \
        cos_deg(zeta) * sin_deg(alpha) * \
        cos_deg(phi - phi0)

    result = psi0 + math.atan(num/denom)

    if result < math.pi/2.0:
        result += math.pi
    if result > math.pi/2.0:
        result -= math.pi

    return math.degrees(result)


def correct(x):
    """Correct value x for rounding errors"""
    tol = 1.0e-7

    if x >= 0.:
        sign = 1.0
    else:
        sign = -1.0

    y = math.fabs(x)
    if math.fabs(y - 1.0) < tol:
        y = 1.0
    if math.fabs(y - 0.0) < tol:
        y = 0.0

    return y * sign


def sin_deg(angle):
    a = math.radians(angle)
    return math.sin(a)


def cos_deg(angle):
    a = math.radians(angle)
    return math.cos(a)


def patchbeam(stokes_i,
              xarray,
              yarray,
              rho,
              pw,
              npatch):

    # patch centre, radians
    pcr = 2.0 * math.pi * np.random.random(npatch)
    pcx = rho * np.sin(pcr)
    pcy = rho * np.cos(pcr)

    # radii = np.zeros((len(xarray), len(yarray)))
    # do numpy stuff
    for px, py in zip(pcx, pcy):
        xtmp = xarray - px
        ytmp = yarray - py

        # x2 = xtmp * xtmp
        # y2 = ytmp * ytmp
        radii = np.sqrt(np.add.outer(xtmp*xtmp, ytmp*ytmp))
        # for i,x in enumerate(xtmp):
        #    radii[i] = np.sqrt(x*x + ytmp*ytmp)

        # Then masks here and operate for the if statements
        index1 = radii/3.0 > pw
        index2 = radii/3.0 <= pw

        radii[index1] = 0.0
        radii[index2] = np.exp(-((radii[index2]/pw)**2)) * 10.

        stokes_i += radii

    return stokes_i


def calcwidth(prof):

    # "area under" profile
    area = np.sum(prof)

    if area == 0.0:
        return 0.0
    # cumulative sum of profile
    sumprof = np.cumsum(prof)

    # find all position > 25%
    indx1 = sumprof > 0.25 * area
    # find values > 75%
    indx2 = sumprof > 0.75 * area

    # first values in each array are the quartile points, basically
    low_x = np.where(indx1)[0][0]
    high_x = np.where(indx2)[0][0]

    # so then compute width of profile
    width = float(high_x) - float(low_x)
    return width / len(prof)
