:mod:`populate` -- Create a population object
=============================================
.. module:: populate
   :synopsis: Create a population object
.. moduleauthor:: Sam Bates <sam.d.bates@gmail.com>

.. class:: Populate

.. method:: Populate.generate(ngen [, surveyList, pDistType, radialDistType, electronModel, pDistPars, siDistPars, lumDistType, lumDistPars, zscale, duty, scindex, gpsArgs, doubleSpec, nostdout])

   The method called by the ``populate.py`` command-line-script

.. method:: Populate.write(outf=populate.model)

   Writes the :class:`~population.Population` model into file outf as a cPickle dump


