
populate.py
===========
.. program:: populate.py

.. cmdoption:: -n <number of pulsars>

   Required: Number of pulsars to generate; or to detect in a survey

.. cmdoption:: -o <output>

   Output file name for population model (def=populate.model)

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

.. cmdoption:: -pdist <distribution name>

   Distribution type for pulse periods (default=lnorm)

   Supported: 'lnorm', 'norm', 'cc97'

.. cmdoption:: -p <mean stddev>
  
   Mean and standard deviation to use in period dist 'lnorm' or 'norm' (def=-2.7, 0.34)

.. cmdoption:: -rdist <radial model>

   Model to use for Galactic radial distribution of pulsars

   Supported: 'lfl06', 'yk04', 'isotropic', 'slab', 'disk'

.. cmdoption:: -dm <Electron model>

   Model to describe the Galactic electron distribution
   
   Supported: 'ne2001', 'lm98'

.. cmdoption:: -gps <fraction 'a'>
  
   Add <fraction> pulsars with GHz-frequency turnovers with index a

.. cmdoption:: -doublespec <fraction alpha1>

   Add <fraction> pulsars with low-frequency (below 1GHz) spectral index of alpha1

.. cmdoption:: --nostdout

   Turn off writing to stdout. Useful for many iterations eg. in large simulations
