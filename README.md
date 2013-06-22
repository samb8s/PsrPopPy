PsrPopPy
========

(For full documentation, see http://samb8s.github.com/PsrPopPy/ or manual.pdf)

Python implementation of PSRPOP (which was written by D Lorimer).
Several of the old models from that (e.g. NE2001) are still included in their native fortran, since re-writing those is beyond the scope of this work. Currently, only a rudimentary makefile is included. This is something that needs work from a willing volunteer!

The main external requirement is [matplotlib](matplotlib.sourceforge.net), which is used for the visualization stuff. It has a very useful API for making simple GUIs, as well as making beautiful plots.

Compiling
---------

I've just incorporated a real makefile for the first time - be sure to edit `makefile.linux` or `makefile.darwin` as appropriate to point to the correct location of your gfortran compiler. After that, it should be as simple as typing `make`.

If this fails, try using the scripts. Inside the lib/fortran directory, edit `make_mac.csh` or `make_linux.csh` (as appropriate) -- change the variable `gf' to point to your local gfortran compiler, then run the script. Fingers crossed, it should all work.

Note mac users should be sure to use a suitable version of gfortran - available from http://gcc.gnu.org/wiki/GFortranBinaries


Usage
=====

I'd recommend adding the `lib/python` directory to your PYTHONPATH and adding the `bin` directory to your PATH. This should then leave you set up to run the code from wherever you like.


A brief description of the "executables" follows.

populate
--------

Create a population mode using user-defined parameters using the snapshot method

evolve
------

Create a populate model using the Ridley & Lorimer evolution method

dosurvey 
--------

Run a population model through a survey. Pre-defined surveys are given, but a user may also create their own.

viewHist
--------

Program for making quick histograms of population models

wxView
------

More detailed population model viewer. Make histograms, scatter plots, etc. All based off the wx backend for matplotlib.
