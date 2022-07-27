#--------1---------2---------3---------4---------5---------6---------7---------8
#
#  File     : Makefile.am
#  Authors  : IY 
#
#  (c) Copyright 2014 RIKEN. All rights reserved.
#
#--------1---------2---------3---------4---------5---------6---------7---------8

if INTERNAL_LAPACK
  INT_LAPACK_LIB=$(LAPACK_LIB3)
else
  INT_LAPACK_LIB=
endif

bin_PROGRAMS = sasa_analysis

sasa_analysis_SOURCES = \
        sasa_option_str.fpp \
        sasa_option.fpp \
        sasa_control.fpp \
        sasa_analyze.fpp\
	sasa_main.fpp

EXTRA_DIST = \
	Makefile.depends 

sasa_analysis_ADDLIB = ../../../lib/lib.a $(INT_LAPACK_LIB)
sasa_analysis_ADDLIBANA = ../../libana/libana.a
sasa_analysis_ADDLIBSPANA = ../libspana/libspana.a

SUFFIXES = .fpp
DEPENDS = Makefile.depends


.fpp.o:
if USEKCOMP
	cp $*.fpp $*.f90
	$(FPP) $(PPFLAGS) $(DEFS) $*.f90
	mv $*.cpp.f90 $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) \
		-I../../../lib -I../../libana -I../libspana $(FCFLAGS) -c $*.f90
else
	$(FPP) $(PPFLAGS) $(DEFS) $*.fpp $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) \
		-I../../../lib -I../../libana -I../libspana $(FCFLAGS) -c $*.f90
endif

sasa_analysis$(EXEEXT) : makelib $(sasa_analysis_OBJECTS)
	$(FC) -o $(bin_PROGRAMS) $(sasa_analysis_OBJECTS) \
		 $(sasa_analysis_ADDLIBSPANA) $(sasa_analysis_ADDLIBANA) $(sasa_analysis_ADDLIB) $(LDFLAGS)

clean:
	rm -f $(bin_PROGRAMS) $(sasa_analysis_OBJECTS) *.f90 *~ *.mod

makelib:
	cd ../../../lib; make
	cd ../../libana; make
	cd ../libspana; make

depend: clean_depend
	python ../../../../fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
