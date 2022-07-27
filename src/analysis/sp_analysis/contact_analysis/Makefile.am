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

bin_PROGRAMS = contact_analysis

contact_analysis_SOURCES = \
        contact_option_str.fpp \
        contact_option.fpp \
        contact_control.fpp \
        contact_analyze.fpp\
	contact_main.fpp

EXTRA_DIST = \
	Makefile.depends 

contact_analysis_ADDLIB = ../../../lib/lib.a $(INT_LAPACK_LIB)
contact_analysis_ADDLIBANA = ../../libana/libana.a
contact_analysis_ADDLIBSPANA = ../libspana/libspana.a

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

contact_analysis$(EXEEXT) : makelib $(contact_analysis_OBJECTS)
	$(FC) -o $(bin_PROGRAMS) $(contact_analysis_OBJECTS) \
		 $(contact_analysis_ADDLIBSPANA) $(contact_analysis_ADDLIBANA) $(contact_analysis_ADDLIB) $(LDFLAGS)

clean:
	rm -f $(bin_PROGRAMS) $(contact_analysis_OBJECTS) *.f90 *~ *.mod

makelib:
	cd ../../../lib; make
	cd ../../libana; make
	cd ../libspana; make

depend: clean_depend
	python ../../../../fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
