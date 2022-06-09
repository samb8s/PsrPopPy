#!/usr/bin/python


def test_1802_2124(pulsar):

    """For testing, give the pulsar the orbital parameters of
        PSR J1802-2142"""

    pulsar.period = 12.648
    pulsar.dm = 149.63

    pulsar.is_binary = True
    pulsar.orbital_period_days = 0.6989
    pulsar.ecc = 2.47e-06
    pulsar.companion_mass_msolar = 0.8
    pulsar.long_peri_degrees = 20.34
    pulsar.inclination_degrees = 78.52
    pulsar.pulsar_mass_msolar = 1.24

    print(pulsar.gb, pulsar.gl)
    pulsar.gb = 0.61
    pulsar.gl = 4.38
    print(pulsar.gb, pulsar.gl)
    pulsar.galcoords = (0.49, 5.2, 0.04)
    pulsar.lum_1400 = 8.54
    pulsar.t_scatter = 0.0
    #return pulsar
