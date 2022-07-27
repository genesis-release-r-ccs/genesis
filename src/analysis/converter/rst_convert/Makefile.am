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

bin_PROGRAMS = rst_convert

rst_convert_SOURCES = \
	fileio_namd_xyz.fpp \
	fileio_namd_xsc.fpp \
	rs_option_str.fpp \
	rs_option.fpp \
	rs_control.fpp \
	rs_setup.fpp \
	rs_convert.fpp \
	rs_main.fpp

EXTRA_DIST = \
	Makefile.depends 

rst_convert_ADDLIB = ../../../lib/lib.a $(INT_LAPACK_LIB)
rst_convert_ADDLIBANA = ../../libana/libana.a

SUFFIXES = .fpp
DEPENDS = Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = $(bin_PROGRAMS) $(rst_convert_OBJECTS) *.f90 *~ *.mod

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

rst_convert$(EXEEXT) : $(rst_convert_ADDLIBANA) $(rst_convert_ADDLIB) \
						$(rst_convert_OBJECTS)
	$(FC) -o $(bin_PROGRAMS) $(rst_convert_OBJECTS) \
		$(rst_convert_ADDLIBANA) $(rst_convert_ADDLIB) $(LDFLAGS)

depend: clean_depend
	python ../../../../fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
