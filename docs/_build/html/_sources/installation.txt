.. _installation:

************
Installation
************

To get started with PsrPopPy there are a few steps you'll need to 
go through.

PsrPopPy is currently supported on Linux and Mac OS X, and for
full feature support, it is recommended to install the 
`Matplotlib <http://matplotlib.org/>`_ package and either use
Python versions >2.6, or install the argparse module.

.. _download_package:

Download the package
====================
The source for PsrPopPy can be downloaded from `GitHub <http://github.com>`_ 
from the `PsrPopPy page <https://github.com/samb8s/PsrPopPy>`_.
The source will contain the Python modules and scripts needed both for
basic and advanced use.

.. _compiling_fortran:

Compile the FORTRAN
===================
Although PsrPopPy is a Python-based package, some of the algorithms
have been kept in their native FORTRAN for speed and ease of 
programming (e.g. the NE2001 electron model, coordinate conversion...).
Therefore, it is necessary to compile the FORTRAN.

Two scripts are provided for just this purpose (though I hope to find
someone willing to contribute configure scripts). From the base 
directory::

  > cd fortran

then edit either ``make_mac.csh`` or ``make_linux.csh``, depending upon 
your system, so that the ``gf`` variable points to your local gfortran/f77
compiler. Running the script::

  > tcsh make_<os>.csh

`should` then compile a series of ``.so`` files in the ``fortran`` directory.
Assuming this went to plan, the installation process is completed.

