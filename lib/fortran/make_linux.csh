#!/bin/tcsh

set gf = /usr/bin/gfortran

$gf -O2 -fPIC -fno-second-underscore -c -I. *.f

$gf -shared -o libne2001.so -fno-second-underscore ne2001.o dm.o psr_ne.o dist.o calc_xyz.o density.o glun.o
$gf -shared -o libtsky.so -fno-second-underscore psr_tsky.o glun.o
$gf -shared -o libykarea.so -fno-second-underscore ykarea.o psrran.o
$gf -shared -o libgetseed.so -fno-second-underscore getseed.o psrran.o clock.o
$gf -shared -o libsla.so -fno-second-underscore galtfeq.o sla.o

rm *.o
