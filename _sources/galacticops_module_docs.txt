:mod:`galacticops` -- Container for functions relating to the Galaxy
====================================================================
.. module:: galactic ops
   :synopsis: Calculate Galactic coords, electron distn, etc
.. moduleauthor:: Sam Bates <sam.d.bates@gmail.com>

.. class:: GalacticOps

.. method:: calc_dtrue((x, y, z))

   Calculate the distance from the Sun to Galactic coords (x, y, z) (NB. tuple)

.. method:: calcXY(r0)

   Given a Galactic radius r0, choose an (x, y) position at random :math:`\theta`

.. method:: ne2001_dist_to_dm(dist, gl, gb)

   Given a distance and Galactic coordinates, calculate DM according to NE2001

.. method:: lm98_dist_to_dm(dist, gl, gb)

   Given a distance and Galactic coordinates, calculate DM according to lm98

.. method:: lb_to_radec(gl, gb)

   Convert Galactic coordinates to equatorial

.. method:: ra_dec_to_lb(ra, dec)

   Convert equatorial coordinates to Galactic

.. method:: tsky( gl, gb, freq)

   Calculate sky temperature at observing frequency freq and at Galactic coordinates gl, gb according to Haslam et al

.. method:: xyz_to_lb((x, y, z))

   Convert the tuple (x, y, z) to Galactic sky coordinates.

   Returns l, b in degrees

.. method:: lb_to_xyz(l, b, dist)

   Convert Galactic sky coordinates at a distance dist to x,y,z coordinates.

   Returns position as a tuple

.. method:: scatter_bhat(dm, scatterindex, freq_mhz)

   Calculate the scatter time according to Bhat et al at. Frequency in MHz, pulsar with dispersion measure dm, and using a scattering spectral index of scatterindex.

   Calculated as

   :math:`\tau = -6.46 + 0.154 \log_{10}({\rm{dm}}) + 1.07\log_{10}({\rm{dm}})^2 + {\rm{scatterindex}} \times \log_{10}(\frac{\rm{freq\_mhz}}{1000})`

   and typically :math:`{\rm{scatterindex}} = -3.86` (but there is an option to vary it!)
