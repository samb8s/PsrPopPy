:mod:`population` -- Creates/stores a population
=================================================
.. module:: population
   :synopsis: The population object
.. moduleauthor:: Sam Bates <sam.d.bates@gmail.com>

.. class:: Population

.. method:: Population.__init__([pDistType, radialDistType, lumDistType, pmean, psigma, simean, sisigma, lummean, lumsigma, zscale, electronModel, gpsFrac, gpsA, brokenFrac, brokenSI, ref_freq])

   Initialise the population object

.. method:: Population.__str__()
   
   Defines how the operation ``print Population`` is performed
   
.. method:: Population.size()

   Returns the number of pulsars in the population object

.. method:: Population.join(poplist)

   Joins each of the populations in list ``poplist`` to the current population

.. method:: Population.write(outf)

   Uses cPickle to dump the population to file ``outf``

.. method:: Population.write_asc(outf)

   Writes the population to an ascii file in the old psrpop way

