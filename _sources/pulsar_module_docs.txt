:mod:`pulsar` -- Creates/stores a pulsar object
===============================================
.. module:: pulsar
   :synopsis: The pulsar object
.. moduleauthor:: Sam Bates <sam.d.bates@gmail.com>

.. class:: Pulsar

.. method:: Pulsar.__init__([period, dm, gl, gb, galCoords, r0, dtrue, lum_1400, spindex, alpha, rho, width_degree, snr, beaming, scindex, gpsFlag, gpsA, brokenFlag, brokenSI])

   Initialise the pulsar object

.. method:: Pulsar.s_1400()
   
   Returns the flux at 1400 MHz, calculated as
   
   :math:`S_{\rm{1400}} = \frac{L_{\rm{1400}}}{D_{\rm{true}}^2}`

.. method:: Pulsar.width_ms()

   Returns the pulse width in milliseconds, calculated as

   :math:`W_{\rm{ms}} = P_{\rm{ms}} \times W_{\rm{degree}} / 360`
