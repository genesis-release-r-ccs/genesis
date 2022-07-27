#--------1---------2---------3---------4---------5---------6---------7---------8
#
#  File     : Makefile.am
#  Authors  : NT
#
#  (c) Copyright 2014 RIKEN. All rights reserved.
#
#--------1---------2---------3---------4---------5---------6---------7---------8
if INTERNAL_LAPACK
  INT_LAPACK_LIB=$(LAPACK_LIB3)
else
  INT_LAPACK_LIB=
endif

SUBDIRS := ../../libana .

bin_PROGRAMS = emmap_generator

emmap_generator_SOURCES = \
        eg_option_str.fpp \
        eg_option.fpp \
        eg_control.fpp \
        eg_setup.fpp \
        eg_generate.fpp \
        eg_main.fpp

EXTRA_DIST = \
	Makefile.depends 

emmap_generator_ADDLIB = ../../../lib/lib.a $(INT_LAPACK_LIB)
emmap_generator_ADDLIBANA = ../../libana/libana.a

SUFFIXES = .fpp
DEPENDS = Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = $(bin_PROGRAMS) $(emmap_generator_OBJECTS) *.f90 *~ *.mod

.fpp.o:
if USEKCOMP
	cp $*.fpp $*.f90
	$(FPP) $(PPFLAGS) $(DEFS) $*.f90
	mv $*.cpp.f90 $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) \
		-I../../../lib -I../../libana $(FCFLAGS) -c $*.f90
else
	$(FPP) $(PPFLAGS) $(DEFS) $*.fpp $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) \
		-I../../../lib -I../../libana $(FCFLAGS) -c $*.f90
endif

emmap_generator$(EXEEXT) : $(emmap_generator_OBJECTS) \
		$(emmap_generator_ADDLIBANA) $(emmap_generator_ADDLIB)
	$(FC) -o $(bin_PROGRAMS) $(emmap_generator_OBJECTS) \
		$(emmap_generator_ADDLIBANA) $(emmap_generator_ADDLIB) $(LDFLAGS)

depend: clean_depend
	python ../../../../fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
