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

bin_PROGRAMS = rpath_generator

rpath_generator_SOURCES = \
        rg_option_str.fpp \
        rg_option.fpp \
        rg_control.fpp \
        rg_setup.fpp \
        rg_generate.fpp \
        rg_main.fpp

EXTRA_DIST = \
	Makefile.depends 

rpath_generator_ADDLIB = ../../../lib/lib.a $(INT_LAPACK_LIB)
rpath_generator_ADDLIBANA = ../../libana/libana.a

SUFFIXES = .fpp
DEPENDS = Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = $(bin_PROGRAMS) $(rpath_generator_OBJECTS) *.f90 *~ *.mod

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

rpath_generator$(EXEEXT) : $(rpath_generator_OBJECTS) \
		$(rpath_generator_ADDLIBANA) $(rpath_generator_ADDLIB)
	$(FC) -o $(bin_PROGRAMS) $(rpath_generator_OBJECTS) \
		$(rpath_generator_ADDLIBANA) $(rpath_generator_ADDLIB) $(LDFLAGS)

depend: clean_depend
	python ../../../../fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
