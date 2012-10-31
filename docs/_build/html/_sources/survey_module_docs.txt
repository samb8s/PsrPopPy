:mod:`survey` -- Read a survey file into a survey object
========================================================
.. module:: survey
   :synopsis: Store survey information
.. moduleauthor:: Sam Bates <sam.d.bates@gmail.com>

.. class:: Survey

.. method:: Survey.__init__(surveyName)

   Read in a (correctly formatted!) survey file

.. method:: Survey.__str__()

   Define how to perform ``print Survey``

.. method:: Survey.nchans()

   Returns the number of channels, calculated as

   :math:`n_{\rm{chans}} = \frac{ {\rm{BW}}_{\rm{total}} }{ {\rm{BM}}_{\rm{chan}} }`

.. method:: Survey.inRegion(pulsar)

   Determines if :class:`~pulsar.Pulsar` is inside survey region. Returns True or False accordingly

.. method:: Survey.inPointing(pulsar)

   Determines if :class:`~pulsar.Pulsar` is inside one of the survey's :class:`pointings<survey.Pointing>`. Returns the offset from beam centre to the pulsar.

.. method:: Survey.SNRcalc(pulsar, pop)

   Calculates the SNR of a :class:`~pulsar.Pulsar` from :class:`~population.Population` ``pop`` in the survey. Returns -1 if pulse is smeared, and -2 if pulsar is outside survey region. SNR is calculated (with familiar terms) as

   :math:`{\rm{SNR}} = \frac{S_{1400} G \sqrt{n_{\rm{pol}} BW \tau}}{\beta T_{\rm{tot}}} \sqrt{\frac{1-\delta}{\delta}} \eta`

   where

   :math:`\eta = \exp(-2.7727 \times {\rm{offset}}^2 / {\rm{fwhm}}^2)`

.. class:: Pointing

.. method:: Pointing.__init__(coord1, coord2, coordtype)

   Converts (coord1, coord2) into correctly formatted (l, b). Coordtype must be
   either ``eq`` or ``gal``. If ``eq``, the RA and Dec are converted internally
