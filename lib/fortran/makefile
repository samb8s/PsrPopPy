include makefile.$(OSTYPE)

all: libne2001.so libykarea.so libsla.so libvxyz.so libgamma.so

library : *.f
	$(GF) -O2 -fPIC -fno-second-underscore -c -I. -std=legacy *.f

libne2001.so : library
	$(LINK.f) -o libne2001.so ne2001.o dm.o psr_ne.o dist.o calc_xyz.o density.o glun.o

libykarea.so : library
	$(LINK.f) -o libykarea.so ykarea.o psrran.o

libsla.so : library
	$(LINK.f) -o libsla.so galtfeq.o sla.o

libvxyz.so : library
	$(LINK.f) -o libvxyz.so vxyz.o rkqc.o rk4.o

libgamma.so : library
	$(LINK.f) -o libgamma.so gamma.o

clean : 
	rm -f *.o *.so *~
