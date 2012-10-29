.. _command_line_docs:

******************************************
Documentation for the command-line scripts
******************************************

Populate.py
===========
.. program:: populate.py

.. cmdoption:: -n <number of pulsars>
  
   Required: Number of pulsars to generate; or to detect in a survey

.. cmdoption:: -surveys <SURVEY NAME(S)>

   List of surveys to use when trying to detect pulsars (default=None)

.. cmdoption:: -z <scale height>

   Scale height of pulsars about Galactic plane, in kpc (default=0.33)

.. cmdoption:: -w <width>

   Pulse width to use when generating pulsars (default=0, use beaming model)

.. cmdoption:: -si <SImean SIsigma>

   Spectral index mean and standard deviation (default=-1.6, 0.35)

.. cmdoption:: -sc <scatter index>

   Spectral index of scattering law to use (default=-3.86, Bhat et al model)
