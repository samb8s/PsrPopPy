#!/bin/tcsh

set gf = /usr/bin/gfortran

$gf -O2 -fPIC -fno-second-underscore -c -I. -std=legacy *.f

$gf -shared -o libne2001.so -fno-second-underscore ne2001.o dm.o psr_ne.o dist.o calc_xyz.o density.o glun.o
$gf -shared -o libykarea.so -fno-second-underscore ykarea.o psrran.o
$gf -shared -o libsla.so -fno-second-underscore galtfeq.o sla.o
$gf -shared -o libvxyz.so -fno-second-underscore vxyz.o rkqc.o rk4.o
$gf -shared -o libgamma.so -fno-second-underscore gamma.o

rm *.o
