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

bin_PROGRAMS = rdf_analysis

rdf_analysis_SOURCES = \
        rdf_option_str.fpp \
        rdf_option.fpp \
        rdf_control.fpp \
        rdf_analyze.fpp\
	rdf_main.fpp

EXTRA_DIST = \
	Makefile.depends 

rdf_analysis_ADDLIB = ../../../lib/lib.a $(INT_LAPACK_LIB)
rdf_analysis_ADDLIBANA = ../../libana/libana.a
rdf_analysis_ADDLIBSPANA = ../libspana/libspana.a

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

rdf_analysis$(EXEEXT) : makelib $(rdf_analysis_OBJECTS)
	$(FC) -o $(bin_PROGRAMS) $(rdf_analysis_OBJECTS) \
		 $(rdf_analysis_ADDLIBSPANA) $(rdf_analysis_ADDLIBANA) $(rdf_analysis_ADDLIB) $(LDFLAGS)

clean:
	rm -f $(bin_PROGRAMS) $(rdf_analysis_OBJECTS) *.f90 *~ *.mod

makelib:
	cd ../../../lib; make
	cd ../../libana; make
	cd ../libspana; make

depend: clean_depend
	python ../../../../fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
