PsrPopPy
========

(For full documentation, see http://samb8s.github.com/PsrPopPy/ or manual.pdf)

Python implementation of PSRPOP (which was written by D Lorimer).
Several of the old models from that (e.g. NE2001) are still included in their native fortran, since re-writing those is beyond the scope of this work. Currently, only a rudimentary makefile is included. This is something that needs work from a willing volunteer!

Note that when compiling, be sure to edit fortran/getpath.f to the actual directory path! If anyone knows a better way to implement this (in FORTRAN, I don't think there is...) please let me know.

The other main requirement is [matplotlib](matplotlib.sourceforge.net), which is used for the visualization stuff. It has a very useful API for making simple GUIs, as well as making beautiful plots.

Compiling
---------

Inside the lib/fortran directory, edit getpath.f and set the path variable to point to the current (ie. the fortran) directory. Still in the fortran directory, edit make_mac.csh or make_linux.csh (as appropriate) -- change the variable `gf' to point to your local gfortran compiler, then run the script. Fingers crossed, it should all work.

Note mac users should be sure to use a suitable version of gfotran - available from http://gcc.gnu.org/wiki/GFortranBinaries
It may also be necessary to do `setenv CFLAGS -m32' before running the compile script on a mac.


Usage
=====

I'd recommend adding the `lib/python` directory to your PYTHONPATH and adding the `bin` directory to your PATH. This should then leave you set up to run the code from wherever you like.


A brief description of the "executables" follows.

populate
--------

Create a population mode using user-defined parameters.

dosurvey 
--------

Run a population model through a survey. Pre-defined surveys are given, but a user may also create their own.

view.py
-------

Program for making quick histograms of population model

wxView.py
---------

More detailed population model viewer. Make histograms, scatter plots, etc. All based off the wx backend for matplotlib.
