###############################################################################
# "makefile" for PSRPOP - a pulsar population synthesis package
# Report any problems compiling to drl@jb.man.ac.uk
# YOU SHOULD NOT NORMALLY HAVE TO CHANGE ANYTHING IN THIS FILE !!!!!!!!!!!!!
###############################################################################
include makefile.darwin
###############################################################################
LIB     =  libpsrpop_$(OSTYPE).a
LIBOBJS =         \
$(LIB)(asind.o) \
$(LIB)(beamfrac.o) \
$(LIB)(beaming.o) \
$(LIB)(calc_xyz.o) \
$(LIB)(clock.o) \
$(LIB)(density.o) \
$(LIB)(deviate.o) \
$(LIB)(dist.o) \
$(LIB)(dmdist.o) \
$(LIB)(dm.o) \
$(LIB)(dmidx.o) \
$(LIB)(draw1d.o) \
$(LIB)(drawgamma.o) \
$(LIB)(galactic.o) \
$(LIB)(galtfeq.o) \
$(LIB)(gamma.o) \
$(LIB)(getpath.o) \
$(LIB)(get_pkspos.o) \
$(LIB)(getrhowidth.o) \
$(LIB)(getseed.o) \
$(LIB)(glun.o) \
$(LIB)(lfunc2.o) \
$(LIB)(lfunc.o) \
$(LIB)(ne2001.o) \
$(LIB)(pksoff.o) \
$(LIB)(pmfft.o) \
$(LIB)(psr_exp.o) \
$(LIB)(psr_ne.o) \
$(LIB)(psr_normal.o) \
$(LIB)(psrran.o) \
$(LIB)(psr_tsky.o) \
$(LIB)(readbinned.o) \
$(LIB)(readfn.o) \
$(LIB)(readsurv.o) \
$(LIB)(rholaw.o) \
$(LIB)(RFMrhowidth.o) \
$(LIB)(sind.o) \
$(LIB)(sla.o) \
$(LIB)(spiralize.o) \
$(LIB)(stopjob.o) \
$(LIB)(surveysn.o) \
$(LIB)(whittle.o) \
$(LIB)(ykarea.o) 

all: survey populate visualize realpsrs pkspsrs underlying popjoin


underlying : library underlying.o $(LIB)
	$(LINK.f) -o $(BIN)/underlying underlying.o  $(LIB)

realpsrs : library realpsrs.o $(LIB)
	$(LINK.f) -o $(BIN)/realpsrs realpsrs.o  $(LIB)

pkspsrs : library pkspsrs.o $(LIB)
	$(LINK.f) -o $(BIN)/pkspsrs pkspsrs.o  $(LIB)

getpsrs : library getpsrs.o $(LIB)
	$(LINK.f) -o $(BIN)/getpsrs getpsrs.o  $(LIB)

namedpsrs : library namedpsrs.o $(LIB)
	$(LINK.f) -o $(BIN)/namedpsrs namedpsrs.o  $(LIB)

popjoin : popjoin.o $(LIB)
	$(LINK.f) -o $(BIN)/popjoin popjoin.o $(LIB)

populate : library populate.o $(LIB)
	$(LINK.f) -o $(BIN)/populate populate.o $(LIB)

survey : library survey.o $(LIB)
	$(LINK.f) -o $(BIN)/survey survey.o  $(LIB)

visualize: visualize.tk
	cp visualize.tk $(BIN)/visualize

radial: library radial.o $(LIB)
	$(LINK.f) -o $(BIN)/radial radial.o $(LIB) $(LPGPLOT)

rlzp: library rlzp.o $(LIB)
	$(LINK.f) -o $(BIN)/rlzp rlzp.o $(LIB)

ner: library ner.o $(LIB)
	$(LINK.f) -o $(BIN)/ner ner.o $(LIB)

nez: library nez.o $(LIB)
	$(LINK.f) -o $(BIN)/nez nez.o $(LIB)

library : $(LIBOBJS)

help:
	./help.csh populate
	./help.csh survey

doc:
	make help
	./document.csh
	latex psrpop
	latex psrpop
	dvipdf psrpop

export:
	./exporter.csh

pulsarcurrent :
	echo "#!`which csh`" >               $(BIN)/pulsarcurrent
	cat pca.csh         >>              $(BIN)/pulsarcurrent
	chmod +x $(BIN)/pulsarcurrent

#
###############################################################################
