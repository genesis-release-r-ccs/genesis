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

bin_PROGRAMS = crd_convert

crd_convert_SOURCES = \
	cc_option_str.fpp \
	cc_option.fpp \
	cc_control.fpp \
	cc_setup.fpp \
	cc_convert.fpp \
	cc_main.fpp

EXTRA_DIST = \
	Makefile.depends 

crd_convert_ADDLIB = ../../../lib/lib.a $(INT_LAPACK_LIB)
crd_convert_ADDLIBANA = ../../libana/libana.a

SUFFIXES = .fpp
DEPENDS = Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = $(bin_PROGRAMS) $(crd_convert_OBJECTS) *.f90 *~ *.mod

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

crd_convert$(EXEEXT) : $(crd_convert_ADDLIBANA) $(crd_convert_ADDLIB) \
						$(crd_convert_OBJECTS)
	$(FC) -o $(bin_PROGRAMS) $(crd_convert_OBJECTS) \
		$(crd_convert_ADDLIBANA) $(crd_convert_ADDLIB) $(LDFLAGS)

depend: clean_depend
	python ../../../../fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
