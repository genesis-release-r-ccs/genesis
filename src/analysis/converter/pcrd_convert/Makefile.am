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

bin_PROGRAMS = pcrd_convert

pcrd_convert_SOURCES = \
	parallel_trj.fpp \
	pc_option_str.fpp \
	pc_option.fpp \
	pc_control.fpp \
	pc_setup.fpp \
	pc_convert.fpp \
	pc_main.fpp

EXTRA_DIST = \
	Makefile.depends 

pcrd_convert_ADDLIB = ../../../lib/lib.a $(INT_LAPACK_LIB)
pcrd_convert_ADDLIBANA = ../../libana/libana.a

SUFFIXES = .fpp
DEPENDS = Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = $(bin_PROGRAMS) $(pcrd_convert_OBJECTS) *.f90 *~ *.mod

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

pcrd_convert$(EXEEXT) :	$(pcrd_convert_ADDLIBANA) $(pcrd_convert_ADDLIB) \
						$(pcrd_convert_OBJECTS)
	$(FC) -o $(bin_PROGRAMS) $(pcrd_convert_OBJECTS) \
		$(pcrd_convert_ADDLIBANA) $(pcrd_convert_ADDLIB) $(LDFLAGS)

depend: clean_depend
	python ../../../../fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
