#--------1---------2---------3---------4---------5---------6---------7---------8
#
#  File     : Makefile.am
#  Authors  : NT
#
#  (c) Copyright 2014 RIKEN. All rights reserved.
#
#--------1---------2---------3---------4---------5---------6---------7---------8

if QSIMULATE
  QSIMULATE_DIR=qsimulate
  JSON_FORTRAN_DIR=json-fortran
else
  QSIMULATE_DIR=
  JSON_FORTRAN_DIR=
endif

SUBDIRS = $(LAPACK_DIR) $(SCALAPACK_DIR) $(QSIMULATE_DIR) $(JSON_FORTRAN_DIR)

noinst_LIBRARIES = lib.a

lib_a_SOURCES = \
	constants.fpp \
	mpi_parallel.fpp \
	messages.fpp \
	random.fpp \
	atom_libs.fpp \
	math_libs.fpp \
	string.fpp \
	table_libs.fpp \
	timers.fpp \
	nbond_list.fpp \
	ffte_fft235.fpp \
	ffte_kernel.fpp \
	ffte_zfft1d.fpp \
	ffte_pzfft3dv.fpp \
	ffte_pdzfft3d.fpp \
	ffte_pzdfft3d.fpp \
	ffte_dzfft3d.fpp \
	ffte_zdfft3d.fpp \
	fft3d.fpp \
	fft3d_1dalltoall.fpp \
	fft3d_2dalltoall.fpp \
	fft3d_slab.fpp \
	fileio.fpp \
	fileio_control.fpp \
	fileio_data.fpp \
	fileio_data_.c \
	fileio_pdb.fpp \
	fileio_crd.fpp \
	fileio_top.fpp \
	fileio_par.fpp \
	fileio_str.fpp \
	fileio_gpr.fpp \
	fileio_psf.fpp \
	fileio_rst.fpp \
	fileio_prmtop.fpp \
	fileio_ambcrd.fpp \
	fileio_gropp.fpp \
	fileio_grotop.fpp \
	fileio_grocrd.fpp \
	fileio_mode.fpp \
	fileio_eef1.fpp \
	fileio_sit.fpp \
	fileio_mrc.fpp \
	fileio_localres.fpp \
	fileio_rstmep.fpp \
	fileio_morph.fpp \
	fileio_table.fpp \
	fileio_spot.fpp \
	hardwareinfo.fpp \
	dSFMT/dSFMT.c \
	dSFMT/dSFMT_if.c \
	dSFMT/dSFMT-common.h \
	dSFMT/dSFMT-params.h \
	dSFMT/dSFMT-params11213.h \
	dSFMT/dSFMT-params1279.h \
	dSFMT/dSFMT-params132049.h \
	dSFMT/dSFMT-params19937.h \
	dSFMT/dSFMT-params216091.h \
	dSFMT/dSFMT-params2203.h \
	dSFMT/dSFMT-params4253.h \
	dSFMT/dSFMT-params44497.h \
	dSFMT/dSFMT-params521.h \
	dSFMT/dSFMT-params86243.h \
	dSFMT/dSFMT.h \
	Lbfgsb.3.0/lbfgsb.f \
	molecules_str.fpp \
	select_atoms_str.fpp \
	fitting_str.fpp \
	molecules.fpp \
	dihedral_libs.fpp \
	select_lexer.fpp \
	select_parser.fpp \
	select_contacts.fpp \
	select_atoms.fpp \
	select.fpp \
	getopt.fpp \
	fitting.fpp \
	fileio_minfo.fpp \
	structure_check.fpp \
	ffte_param.h

EXTRA_DIST = \
	Makefile.depends 

SUFFIXES = .fpp .c .f
DEPENDS	= Makefile.depends

CLEANFILES = $(noinst_LIBRARIES) $(lib_a_OBJECTS) *.f90 *~ *.mod

ffte_kernel.o : ffte_kernel.fpp
if USEKCOMP
	cp $*.fpp $*.f90
	$(FPP) $(PPFLAGS) $(DEFS) $*.f90
	mv $*.cpp.f90 $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) -I. $(FCFLAGS_FFTE_KERNEL_FPP__PGI__SP) -c $*.f90
else
	$(FPP) $(PPFLAGS) $(DEFS) $*.fpp $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) -I. $(FCFLAGS_FFTE_KERNEL_FPP__PGI__SP) -c $*.f90
endif

.fpp.o:
if USEKCOMP
	cp $*.fpp $*.f90
	$(FPP) $(PPFLAGS) $(DEFS) $*.f90
	mv $*.cpp.f90 $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) -I. $(FCFLAGS) -c $*.f90
else
	$(FPP) $(PPFLAGS) $(DEFS) $*.fpp $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) -I. $(FCFLAGS) -c $*.f90
endif

.f.o:
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) -I. $(FCFLAGS) -c -o $@ $*.f 

lib.a : $(lib_a_OBJECTS)
	$(AR) $(ARFLAGS) cr $(noinst_LIBRARIES) $(lib_a_OBJECTS)

depend: clean_depend
	python ../../fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
