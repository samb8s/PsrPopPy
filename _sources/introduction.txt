.. _installation:

************
Introduction
************

What is PsrPopPy?
=================
PsrPopPy is a Python (for the most part) implementation of Duncan Lorimer's `PSRPOP 
code <http://psrpop.sourceforge.net/>`_. All of the user-facing (in normal 
circumstances) code is written in Python, with some remaining FORTRAN functions 
(e.g. NE2001, coordinate conversion) for speed.

Why is PsrPopPy?
================
For the development of a research project, I was modifying the PSRPOP code, but 
found it somewhat tricky. I decided to re-write the code with much heavier reliance on 
functions and with added OO design. This makes modifying the code and addition of new
features much more simple, with little speed loss since the difficult number crunching
is still done in FORTRAN. The added bonus of re-writing the code is the detection
of a number of bugs in the original code, which have hopefully been eliminated.

Who can contribute?
===================
In short - anybody. The code is up on github and I'll be happy to accept suggestions for 
future modifications and improvements. 

Acknowledgements
================
Many thanks to Duncan Lorimer for giving his blessing to this project and of course
for generating the codebase which has inspired this project. 
