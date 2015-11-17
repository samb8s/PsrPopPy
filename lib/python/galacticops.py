#!/usr/bin/python

import os

import math
import random

import ctypes as C


# get the FORTRAN libraries
__dir__ = os.path.dirname(os.path.abspath(__file__))
__libdir__ = os.path.dirname(__dir__)
fortranpath = os.path.join(__libdir__, 'fortran')

ne2001lib = C.CDLL(os.path.join(fortranpath, 'libne2001.so'))
ne2001lib.dm_.restype = C.c_float

slalib = C.CDLL(os.path.join(fortranpath, 'libsla.so'))
vxyzlib = C.CDLL(os.path.join(fortranpath, 'libvxyz.so'))

yklib = C.CDLL(os.path.join(fortranpath, 'libykarea.so'))
yklib.ykr_.restype = C.c_float
yklib.llfr_.restype = C.c_float

# BEGIN FUNCTION DEFINITIONS


def vxyz(pulsar):
    """Evolve a pulsar through galactic potential"""
    x, y, z = pulsar.galCoords
    vx, vy, vz = pulsar.vx, pulsar.vy, pulsar.vz
    age_Myr = pulsar.age/1.0E6

    x, y, z = C.c_float(x), C.c_float(y), C.c_float(z)
    vx, vy, vz = C.c_float(vx), C.c_float(vy), C.c_float(vz)
    age_Myr = C.c_float(age_Myr)
    bound = C.c_long(0)

    # run the evolution code
    vxyzlib.vxyz_(C.byref(C.c_float(0.005)),
                  C.byref(x),
                  C.byref(y),
                  C.byref(z),
                  C.byref(age_Myr),
                  C.byref(vx),
                  C.byref(vy),
                  C.byref(vz),
                  C.byref(x),
                  C.byref(y),
                  C.byref(z),
                  C.byref(vx),
                  C.byref(vy),
                  C.byref(vz),
                  C.byref(bound)
                  )

    # convert the output C types to python numbers
    pulsar.galCoords = x.value, y.value, z.value
    pulsar.vx = vx.value
    pulsar.vy = vy.value
    pulsar.vz = vz.value


def calc_dtrue((x, y, z)):
    """Calculate true distance to pulsar from the sun."""
    rsun = 8.5  # kpc
    return math.sqrt(x*x + (y-rsun)*(y-rsun) + z*z)


def calcXY(r0):
    """Calculate the X, Y, Z alactic coords for the pulsar."""
    # calculate a random theta in a unifrom distribution
    theta = 2.0 * math.pi * random.random()

    # calc x and y
    x = r0 * math.cos(theta)
    y = r0 * math.sin(theta)

    return x, y


def ne2001_dist_to_dm(dist, gl, gb):
    """Use NE2001 distance model."""
    # expects -180 < l < 180
    dist = C.c_float(dist)
    gl = C.c_float(gl)
    gb = C.c_float(gb)
    inpath = C.create_string_buffer(fortranpath)
    linpath = C.c_int(len(fortranpath))
    return ne2001lib.dm_(C.byref(dist),
                         C.byref(gl),
                         C.byref(gb),
                         C.byref(C.c_int(4)),
                         C.byref(C.c_float(0.0)),
                         C.byref(inpath),
                         C.byref(linpath))


def lmt85_dist_to_dm(dist, gl, gb):
    """ Use Lyne, Manchester & Taylor distance model to compute DM."""
    dist = C.c_float(dist)
    gl = C.c_float(gl)
    gb = C.c_float(gb)
    # passing path to fortran dir and the length of
    # this path --- removes need to edit getpath.f
    # during installation
    inpath = C.create_string_buffer(fortranpath)
    linpath = C.c_int(len(fortranpath))

    return ne2001lib.dm_(C.byref(dist),
                         C.byref(gl),
                         C.byref(gb),
                         C.byref(C.c_int(0)),
                         C.byref(C.c_float(0.0)),
                         C.byref(inpath),
                         C.byref(linpath))


def ne2001_get_smtau(dist, gl, gb):
    """Use NE2001 model to get the DISS scattering timescale"""
    dist = C.c_float(dist)

    # gl gb need to be in radians
    gl = C.c_float(math.radians(gl))
    gb = C.c_float(math.radians(gb))

    # call dmdsm and get the value out of smtau
    ndir = C.c_int(-1)
    sm = C.c_float(0.)
    smtau = C.c_float(0.)
    inpath = C.create_string_buffer(fortranpath)
    linpath = C.c_int(len(fortranpath))
    ne2001lib.dmdsm_(C.byref(gl),
                     C.byref(gb),
                     C.byref(ndir),
                     C.byref(C.c_float(0.0)),
                     C.byref(dist),
                     C.byref(C.create_string_buffer(' ')),
                     C.byref(sm),
                     C.byref(smtau),
                     C.byref(C.c_float(0.0)),
                     C.byref(C.c_float(0.0)),
                     C.byref(inpath),
                     C.byref(linpath)
                     )
    return sm.value, smtau.value


def ne2001_scint_time_bw(dist, gl, gb, freq):
    sm, smtau = ne2001_get_smtau(dist, gl, gb)
    if smtau <= 0.:
        scint_time = None
    else:
        # reference: eqn (46) of Cordes & Lazio 1991, ApJ, 376, 123
        # uses coefficient 3.3 instead of 2.3. They do this in the code
        # and mention it explicitly, so I trust it!
        scint_time = 3.3 * (freq/1000.)**1.2 * smtau**(-0.6)
    if sm <= 0.:
        scint_bw = None
    else:
        # and eq 48
        scint_bw = 223. * (freq/1000.)**4.4 * sm**(-1.2) / dist

    return scint_time, scint_bw


def lb_to_radec(l, b):
    """Convert l, b to RA, Dec using SLA fortran (should be faster)."""
    ra = C.c_float(0.)
    dec = C.c_float(0.)
    l = C.c_float(l)
    b = C.c_float(b)
    # call with final arg 1 to do conversion in right direction!
    slalib.galtfeq_(C.byref(l),
                    C.byref(b),
                    C.byref(ra),
                    C.byref(dec),
                    C.byref(C.c_int(1)))
    return ra.value, dec.value


def radec_to_lb(ra, dec):
    """Convert RA, Dec to l, b using SLA fortran.
    Be sure to return l in range -180 -> +180"""
    l = C.c_float(0.)
    b = C.c_float(0.)
    ra = C.c_float(ra)
    dec = C.c_float(dec)
    # call with arg = -1 to convert in reverse!
    slalib.galtfeq_(C.byref(l),
                    C.byref(b),
                    C.byref(ra),
                    C.byref(dec),
                    C.byref(C.c_int(-1)))
    if l.value > 180.:
        l.value -= 360.
    return l.value, b.value


def xyz_to_lb((x, y, z)):
    """ Convert galactic xyz in kpc to l and b in degrees."""
    rsun = 8.5  # kpc

    # distance to pulsar
    d = math.sqrt(x*x + (rsun-y)*(rsun-y) + z*z)
    # radial distance
    b = math.asin(z/d)

    # take cosine
    dcb = d * math.cos(b)

    if y <= rsun:
        if math.fabs(x/dcb) > 1.0:
            l = 1.57079632679
        else:
            l = math.asin(x/dcb)
    else:
        if math.fabs(x/dcb) > 1.0:
            l = 0.0
        else:
            l = math.acos(x/dcb)

        l += 1.57079632679
        if x < 0.:
            l -= 6.28318530718

    # convert back to degrees
    l = math.degrees(l)
    b = math.degrees(b)

    # convert to -180 < l < 180
    if l > 180.0:
        l -= 360.0

    return l, b


def lb_to_xyz(gl, gb, dist):
    """ Convert galactic coords to Galactic XYZ."""
    rsun = 8.5  # kpc

    l = math.radians(gl)
    b = math.radians(gb)

    x = dist * math.cos(b) * math.sin(l)
    y = rsun - dist * math.cos(b) * math.cos(l)
    z = dist * math.sin(b)

    return (x, y, z)


def scatter_bhat(dm,
                 scatterindex=-3.86,
                 freq_mhz=1400.0):
    """Calculate bhat et al 2004 scattering timescale for freq in MHz."""
    logtau = -6.46 + 0.154 * math.log10(dm)
    logtau += 1.07 * math.log10(dm)*math.log10(dm)
    logtau += scatterindex * math.log10(freq_mhz/1000.0)

    # return tau with power scattered with a gaussian, width 0.8
    return math.pow(10.0, random.gauss(logtau, 0.8))


def scale_bhat(timescale,
               frequency,
               scaling_power=3.86):
    """Scale the scattering timescale from 1.4 GHz to frequency"""

    return timescale * (frequency/1400.0)**scaling_power


def _glgboffset(gl1, gb1, gl2, gb2):
    """
    Calculate the angular distance (deg) between two
    points in galactic coordinates
    """
    # Angular offset in polar coordinates
    # taken brazenly from
    # http://www.atnf.csiro.au/people/Tobias.Westmeier/tools_hihelpers.php

    # requires gb conversion from +90 -> -90 to 0 -> 180
    gb1 = 90.0 - gb1
    gb2 = 90.0 - gb2

    term1 = math.cos(math.radians(gb1)) * math.cos(math.radians(gb2))
    term2 = math.sin(math.radians(gb1)) * math.sin(math.radians(gb2))
    term3 = math.cos(math.radians(gl1) - math.radians(gl2))
    cosalpha = term1 + term2*term3

    return math.degrees(math.acos(cosalpha))


def seed():
    return C.c_int(random.randint(1, 9999))


def slabdist():
    x = -15.0 + random.random()*30.0
    y = -15.0 + random.random()*30.0
    z = -5.0 + random.random() * 10.0

    return (x, y, z)


def diskdist():
    x = -15.0 + random.random()*30.0
    y = -15.0 + random.random()*30.0
    return (x, y, 0.0)


def lfl06():
    """lfl06 model, using Y&K"""
    return yklib.llfr_(C.byref(seed()))


def ykr():
    """ Y&K Model"""
    return yklib.ykr_(C.byref(seed()))


def spiralize(r):
    """ Make spiral arms, as seen in Fuacher-Giguere & Kaspi 2006"""

    # definitions
    k_list = [4.25, 4.25, 4.89, 4.89]
    r0_list = [3.48, 3.48, 4.9, 4.9]
    theta0_list = [1.57, 4.71, 4.09, 0.95]

    # select a spiral arm ( 1->4)
    arm = random.choice([0, 1, 2, 3])
    k = k_list[arm]
    r0 = r0_list[arm]
    theta0 = theta0_list[arm]

    # pick an angle
    theta = k * math.log(r/r0) + theta0

    # blurring angle
    angle = 2*math.pi * random.random() * math.exp(-0.35 * r)

    if random.random() < 0.5:
        angle = 0 - angle

    # modify theta
    theta += angle

    # blur in radial direction a little
    dr = math.fabs(random.gauss(0.0, 0.5 * r))
    angle = random.random() * 2.0 * math.pi
    dx = dr * math.cos(angle)
    dy = dr * math.sin(angle)

    x = r * math.cos(theta) + dx
    y = r * math.cos(theta) + dy

    return x, y


def _double_sided_exp(scale, origin=0.0):
    """Exponential distribution around origin, with scale height scale."""
    if scale == 0.0:
        return origin

    rn = random.random()
    sign = random.choice([-1.0, 1.0])

    return origin + sign * scale * math.log(rn)


def readtskyfile():
    """Read in tsky.ascii into a list from which temps can be retrieved"""

    tskypath = os.path.join(fortranpath, 'lookuptables/tsky.ascii')
    tskylist = []
    with open(tskypath) as f:
        for line in f:
            str_idx = 0
            while str_idx < len(line):
                # each temperature occupies space of 5 chars
                temp_string = line[str_idx:str_idx+5]
                try:
                    tskylist.append(float(temp_string))
                except:
                    pass
                str_idx += 5

    return tskylist
