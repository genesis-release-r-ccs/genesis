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

bin_PROGRAMS = qmmm_generator

qmmm_generator_SOURCES = \
        qg_option_str.fpp \
        qg_option.fpp \
        qg_control.fpp \
        qg_setup.fpp \
        qg_extract.fpp \
        qg_main.fpp

EXTRA_DIST = \
	Makefile.depends 

qmmm_generator_ADDLIB = ../../../lib/lib.a $(INT_LAPACK_LIB)
qmmm_generator_ADDLIBANA = ../../libana/libana.a

SUFFIXES = .fpp
DEPENDS = Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = $(bin_PROGRAMS) $(qmmm_generator_OBJECTS) *.f90 *~ *.mod

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

qmmm_generator$(EXEEXT) : $(qmmm_generator_OBJECTS) \
				$(qmmm_generator_ADDLIBANA) $(qmmm_generator_ADDLIB) 
	$(FC) -o $(bin_PROGRAMS) $(qmmm_generator_OBJECTS) \
		$(qmmm_generator_ADDLIBANA) $(qmmm_generator_ADDLIB) $(LDFLAGS)

depend: clean_depend
	python ../../../../fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
