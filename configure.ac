#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(genesis, [2.1.1], genesis@riken.jp)
AC_CONFIG_SRCDIR([src/spdyn/spdyn.fpp])
AC_CONFIG_HEADER([src/config.h])

AM_INIT_AUTOMAKE([foreign subdir-objects])

AC_PREFIX_DEFAULT([$PWD/])

m4_ifdef([AM_EXTRA_RECURSIVE_TARGETS], [AM_EXTRA_RECURSIVE_TARGETS([depend])], [])

## make set
AC_PROG_MAKE_SET


## ranlib
AC_PROG_RANLIB
AC_LANG_PUSH(Fortran)

## real precision
#
DEFINED_VARIABLES=""

AC_ARG_ENABLE([double],
    [AS_HELP_STRING(
      [--enable-double],
      [enable double precision (default = yes)])],
    [],[enable_double=yes])

AC_ARG_ENABLE([mixed],
    [AS_HELP_STRING(
      [--enable-mixed],
      [enable mixed precision (default = no)])],
    [],[enable_mixed=no])

AC_ARG_ENABLE([single],
    [AS_HELP_STRING(
      [--enable-single],
      [enable single precision (default = no)])],
    [],[enable_single=no])

if test x"${enable_mixed}" = x"yes"; then
  enable_double=no
fi

if test x"${enable_single}" = x"yes"; then
  enable_double=no
  enable_mixed=no
fi

mainly_single=no
if test x"${enable_single}" = x"yes"; then
  AC_DEFINE(_SINGLE, 1, [defined if _SINGLE is used.])
  DEFINED_VARIABLES+=" -D_SINGLE"
  mainly_single=yes
elif test x"${enable_mixed}" = x"yes"; then
  AC_DEFINE(_MIXED, 1, [defined if _MIXED is used.])
  DEFINED_VARIABLES+=" -D_MIXED"
  mainly_single=yes
fi

## integer precision for degree of freedom
#
DEFINED_VARIABLES=""

AC_ARG_ENABLE([large_int],
    [AS_HELP_STRING(
      [--enable-large_int],
      [enable 8byte integer precision (default = no)])],
    [],[enable_large_int=no])

mainly_large_int=no
if test x"${enable_large_int}" = x"yes"; then
  AC_DEFINE(_LARGE_INT, 1, [defined if _LARGE_INT is used.])
  DEFINED_VARIABLES+=" -D_LARGE_INT"
  mainly_large_int=yes
fi

DEFINED_VARIABLES=""

## Debug
#

AC_ARG_ENABLE([debug],
    [AS_HELP_STRING(
      [--enable-debug@<:@=ARG@:>@],
      [set Debug level (from 1 to 4, default = 1)])],
    [
      if test x"${enableval}" = x"yes"; then
	    enable_debug=1
      elif test x"${enableval}" = x"1" -o \
                x"${enableval}" = x"2" -o \
                x"${enableval}" = x"3" -o \
                x"${enableval}" = x"4"; then
        enable_debug=${enableval}
      else
        AC_MSG_ERROR([Not supported Debug level ${enableval}.])

      fi
    ],
	[
      enable_debug=0
    ])


## MPI
#

AC_ARG_ENABLE([mpi],
        [AS_HELP_STRING([--disable-mpi], [disable MPI parallel])],
        [], [enable_mpi=yes])

if test x"${enable_mpi}" = x"yes"; then
  AC_DEFINE(HAVE_MPI_GENESIS, 1, [defined if MPI is used.])
  DEFINED_VARIABLES+=" -DHAVE_MPI_GENESIS"
fi
AM_CONDITIONAL(SPDYN, test x"${enable_mpi}" = x"yes")
AM_CONDITIONAL(CGDYN, test x"${enable_mpi}" = x"yes")

AC_ARG_WITH([msmpi],
	[AS_HELP_STRING([--with-msmpi],
	[use MSMPI (default: does not use MSMPI)])],
	[], [with_msmpi=no])


## fortran compiler
#

if test x"${enable_mpi}" = x"yes"; then
  genesis_default_fc="mpiifort mpifrtpx mpifrt mpif90"
else
  genesis_default_fc="ifort frtpx frt ifort gfortran"
fi

# check usablility
AC_PROG_FC($genesis_default_fc)
AC_PROG_F77($genesis_default_fc)

## C++ compiler
qsimulate_default_cxx="mpicxx mpic++"
AC_PROG_CXX($qsimulate_default_cxx)

# check actual compiler
AC_MSG_CHECKING([for actual compiler])

if test x"${FC}" = x"mpifrtpx"; then
  VER_OPT="-V"
else
  VER_OPT="--version"
fi


fc_ver_info=$(${FC} -E ${VER_OPT})
for vtok in ${fc_ver_info} ; do

  if test x"${vtok}" = x"ifort"; then
    FC_ACT="ifort"
    break

  elif test x"${vtok}" = x"GNU"; then
    FC_ACT="gfortran"
    break

  elif test x"${vtok}" = x"pgfortran"; then
    FC_ACT="pgf90"
    break

  elif test x"${vtok}" = x"pgf90"; then
    FC_ACT="pgf90"
    break

  elif test x"${vtok}" = x"pgf95"; then
    FC_ACT="pgf90"
    break

  fi

done

if test -z "${FC_ACT}"; then
  FC_ACT=${FC}
fi
F77_ACT=FC_ACT

AC_MSG_RESULT([$FC_ACT])

AC_MSG_CHECKING([for MPI type])
mpitype="no"
if test x"${with_msmpi}" = x"yes"; then
  mpitype="MS"
elif test x"${enable_mpi}"=x"yes" ; then
  mpi_ver_info=$(${FC} --showme:version)
  prev=""
  for vtok in ${mpi_ver_info} ; do
    if echo ${vtok} | egrep -q "MPI"; then
      mpitype=${prev}
    fi
    prev=${vtok}
  done
fi

version_warning=yes
# check gfortran & intel version
if test x"${FC_ACT}" = x"frtpx" ; then
  if test x"${enable_mpi}" != x"yes"; then
    AC_MSG_ERROR([disable_mpi is not available on Fugaku])
  fi
elif test x"${FC_ACT}" = x"ifort" ; then

  AC_MSG_CHECKING([for fortran version])

  version=""
  build=""

  fc_ver_info=$(${FC} -E ${VER_OPT})
  for vtok in ${fc_ver_info} ; do
    if echo ${vtok} | egrep -q "^@<:@0-9@:>@+\.@<:@0-9@:>@+\.@<:@0-9@:>@+"; then
      version=${vtok}
    fi
    if echo ${vtok} | egrep -q "^@<:@0-9@:>@{8}"; then
      build=${vtok}
    fi
  done

  if test x"${mpitype}" = x"no" -a x"${with_msmpi}" = x"no" ; then
    fc_ver_info=$(${FC} -E -v)
	prev=""
    for vtok in ${fc_ver_info} ; do
      if echo ${vtok} | egrep -q "MPI"; then
        mpitype=${prev}
      fi
	  prev=${vtok}
    done
  fi

  AC_MSG_RESULT([$version])

  FC_GNU_VER=`echo ${version}\
  | awk -F'[.]' '{printf "%02d%02d%02d",$1,$2,$3}'`


  if test ${FC_GNU_VER} -eq 120100 -a ${build} -eq 20110811 ; then
    AC_MSG_ERROR([ifort version 12.1.0 (20110811) has a serious problem in GENESIS])

  elif test ${FC_GNU_VER} -eq 150001 -a ${build} -eq 20141023 ; then
    AC_MSG_ERROR([ifort version 15.0.1 (20141023) has a serious problem in GENESIS])
  fi

elif test x"${FC_ACT}" = x"gfortran"; then

  FC_GNU_VER=`echo "__GNUC__ __GNUC_MINOR__ __GNUC_PATCHLEVEL__"\
   | ${FC_ACT} -E -P - | xargs printf "%02d%02d%02d"`

  if test ${FC_GNU_VER} -lt 040407; then
    AC_MSG_ERROR([Gfortran version must be >= 4.4.7])
  fi

elif test x"${FC_ACT}" = x"pgf90"; then
  version=""
  build=""

  fc_ver_info=$(${FC} -E ${VER_OPT} )
  for vtok in ${fc_ver_info} ; do
    if echo ${vtok} | egrep -q "^@<:@0-9@:>@+\.@<:@0-9@:>@+\-@<:@0-9@:>@+$"; then
      version=${vtok}
    fi
  done
  FC_GNU_VER=`echo ${version}\
   | sed -e 's/\(@<:@0-9@:>@*\).\(@<:@0-9@:>@*\)-/\1 \2/'\
   | xargs printf "%02d%02d"`

  if test ${FC_GNU_VER} -lt 1610 ; then
    AC_MSG_ERROR([pgi version must be >= 16.1])
  fi

fi
AC_MSG_RESULT([${mpitype} MPI])


## C compiler
#

if test x"${enable_mpi}" = x"yes"; then
  genesis_default_cc="mpiicc mpifccpx mpifcc mpicc"
else
  genesis_default_cc="icc fccpx fcc gcc"
fi

# check usablility
AC_PROG_CC($genesis_default_cc)

# check actual compiler
AC_MSG_CHECKING([for actual compiler])

cc_ver_info=$(${CC} -E ${VER_OPT})
for vtok in ${cc_ver_info} ; do

  if test x"${vtok}" = x"icc"; then
    CC_ACT="icc"
    break

  elif test x"${vtok}" = x"fccpx"; then
    CC_ACT="fccpx"
    if test x"${enable_mpi}" = x"yes"; then
      CC_ACT="mpifccpx"
    fi
    break

  elif test x"${vtok}" = x"gcc"; then
    CC_ACT="gcc"
    break

  elif test x"${vtok}" = x"pgcc"; then
    CC_ACT="pgcc"
    break

  fi

done

if test -z "${CC_ACT}"; then
  CC_ACT=${CC}
fi

AC_MSG_RESULT([$CC_ACT])


## preprocessor
#

if test -n "${FPP}"; then

  if test ` echo ${FPP} | grep / ` ; then
    AC_CHECK_FILES([${FPP}], check_fpp=${FPP})

  else
    AC_CHECK_PROGS(check_fpp, [${FPP}])

  fi

  if test -z "${check_fpp}"; then
    AC_MSG_ERROR([Preprocessor ${FPP} is not found.])

  fi

else

  if test x"${FC_ACT}" = x"frtpx" -o x"${FC_ACT}" = x"mpifrtpx"; then
    AC_CHECK_PROGS(FPP, frtpx fpp cpp)

  elif test x"${FC_ACT}" = x"frt" -o x"${FC_ACT}" = x"mpifrt"; then
    AC_CHECK_PROGS(FPP, frt fpp cpp)

  else
    AC_CHECK_PROGS(FPP, fpp cpp)

  fi

  if test -z "${FPP}"; then
    AC_MSG_ERROR([There is no any Preprocessor for Fortran])
  fi
fi


## OpenMP
#

AC_ARG_ENABLE([openmp],
	[AS_HELP_STRING([--disable-openmp], [disable OpenMP])],
	[], [enable_openmp=yes])

if test x"${enable_openmp}" = x"yes"; then
  AC_DEFINE(OMP, 1, [defined if OpenMP is used.])
  DEFINED_VARIABLES+=" -DOMP"

  if   test x"${FC_ACT}" = x"ifort"; then
    OPT_OPENMP="-openmp"
    if test ${FC_GNU_VER} -gt 160000; then
       OPT_OPENMP="-qopenmp"
    fi
    if test ${FC_GNU_VER} -gt 190000; then
       AC_DEFINE(IFORT19, 1, [defined if ifort version (>19) is used.])
       DEFINED_VARIABLES+=" -DIFORT19"
    fi

  elif test x"${FC_ACT}" = x"gfortran"; then
    OPT_OPENMP="-fopenmp"

  elif test x"${FC_ACT}" = x"pgf90"; then
    OPT_OPENMP="-mp"

  elif test x"${FC_ACT}" = x"frtpx"; then
    OPT_OPENMP="-Kopenmp"

  elif test x"${FC_ACT}" = x"mpifrtpx"; then
    OPT_OPENMP="-Kopenmp"

  elif test x"${FC_ACT}" = x"frt"; then
    OPT_OPENMP="-Kopenmp"

  elif test x"${FC_ACT}" = x"mpifrt"; then
    OPT_OPENMP="-Kopenmp"

  fi
  if test x"${with_msmpi}" = x"yes"; then
     AC_MSG_ERROR([OpenMP is not available with MSMPI])
  fi
fi

## ffte
#

AC_DEFINE(FFTE, 1, [defined if FFTE is used.])
DEFINED_VARIABLES+=" -DFFTE"


## lapack
#

AC_ARG_WITH([lapack],
	[AS_HELP_STRING([--without-lapack],
	[do not use LAPACK (default: uses LAPACK)])],
	[], [with_lapack=yes])

AC_ARG_VAR([LAPACK_PATH], 
    [LAPACK library path (linker flags, e.g. -L<lapack dir>)])

AC_ARG_VAR([LAPACK_LIBS],
    [LAPACK library files (libraries, e.g. -l<lapack library>)])


if test x"${with_lapack}" == x"yes"; then

  AC_DEFINE(LAPACK, 1, [defined if LAPACK is used.])
  DEFINED_VARIABLES+=" -DLAPACK"


  blaslink="-lblas"
  lapacklib=""

  if test -n "${LAPACK_PATH}"; then

    lapack_path_string=${LAPACK_PATH}

    if test x"${FC_ACT}" = x"gfortran" -o \
  	      x"${FC_ACT}" = x"pgf90" ; then
      lapack_lib_path=""
      blas_lib_path=""
      for vtok in ${lapack_path_string} ; do

        if test -z "${lapack_lib_path}"; then
          AC_CHECK_FILE([${vtok}/liblapack.a],[lapack_lib_path=${vtok} ; lapacklib="lapack"],[])
          if test -z "${lapack_lib_path}"; then
            AC_CHECK_FILE([${vtok}/liblapack.so],[lapack_lib_path=${vtok} ; lapacklib="lapack"],[])
          fi
        fi

        if test -z "${blas_lib_path}"; then
          AC_CHECK_FILE([${vtok}/libblas.a],[blas_lib_path=${vtok} ; blaslink="-lblas"],[])
          if test -z "${blas_lib_path}"; then
            AC_CHECK_FILE([${vtok}/libblas.so],[blas_lib_path=${vtok} ; blaslink="-lblas"],[])
          fi
          if test -z "${blas_lib_path}"; then
            AC_CHECK_FILE([${vtok}/librefblas.a],[blas_lib_path=${vtok} ; blaslink="-lrefblas"],[])
          fi
          if test -z "${blas_lib_path}"; then
            AC_CHECK_FILE([${vtok}/librefblas.so],[blas_lib_path=${vtok} ; blaslink="-lrefblas" ],[])
          fi
        fi
      done

      if test -z "${lapack_lib_path}"; then
         AC_MSG_ERROR([Lapack is not found.])
      fi
      if test -z "${blas_lib_path}"; then
         AC_MSG_ERROR([BLAS is not found])
      fi
      if test x"${lapack_lib_path}" == x"${blas_lib_path}"; then
        OPT_LD_LAPACK+=" -L${lapack_lib_path}"
      else
        OPT_LD_LAPACK+=" -L${lapack_lib_path} -L${blas_lib_path}"
      fi

    else

      for vtok in ${lapack_path_string} ; do
        OPT_LD_LAPACK+=" -L${vtok}"
  	  done

    fi

  fi

  if test -n "${LAPACK_LIBS}"; then
    OPT_LD_LAPACK+="${LAPACK_LIBS}"

  elif test x"${FC_ACT}" = x"ifort"; then
    OPT_LD_LAPACK+=" -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_lapack95_lp64"
    lapacklib="mkl_lapack95_lp64"

    LDFLAGS=${OPT_LD_LAPACK}
    AC_CHECK_LIB([${lapacklib}], [dsyev],
          [], [AC_MSG_ERROR([lapack is not correct])],[-lmkl_intel_lp64 -lmkl_sequential -lmkl_core])
    LDFLAGS=""

  elif test x"${FC_ACT}" = x"pgf90"; then
#    OPT_LD_LAPACK+=" -llapack ${blaslink}"
    lapacklib="lapack"
    LDFLAGS=${OPT_LD_LAPACK}
    AC_CHECK_LIB([${lapacklib}], [dsyev],
          [OPT_LD_LAPACK+=" -llapack ${blaslink}"], [AC_MSG_ERROR([lapack is not correct])],["${blaslink} ${OPT_LD_LAPACK}"])
    LDFLAGS=""

  elif test x"${FC_ACT}" = x"gfortran"; then
#    OPT_LD_LAPACK+=" -llapack ${blaslink}"
    lapacklib="lapack"
    LDFLAGS=${OPT_LD_LAPACK}
    AC_CHECK_LIB([${lapacklib}], [dsyev],
          [OPT_LD_LAPACK+=" -llapack ${blaslink}"], [AC_MSG_ERROR([lapack is not correct])],["${blaslink}"])
    LDFLAGS=""
  fi


fi
with_internallapack=no
AM_CONDITIONAL(INTERNAL_LAPACK, test x"${with_internallapack}" = x"yes")
AM_CONDITIONAL(INTERNAL_SCALAPACK, test x"${with_internallapack}" = x"yes" -a x"${enable_mpi}" = x"yes")


## GPU precision
#
AC_ARG_ENABLE([gpu],
    [AS_HELP_STRING(
      [--enable-gpu],
      [enable gpu (default = no)])],
    [],[enable_gpu=no])

if test x"${enable_gpu}" = x"yes"; then
  AC_DEFINE(USE_GPU, 1, [defined if gpu is used.])
  DEFINED_VARIABLES+=" -DUSE_GPU"
fi


## parallel I/O
#
AC_ARG_ENABLE([parallel_IO],
	[AS_HELP_STRING([--enable-parallel_IO], [enable parallel_IO (defalt: no)])],
	[], [enable_parallel_IO=no])
if test x"${FC_ACT}" = x"frtpx" -o x"${FC_ACT}" = x"mpifrtpx"; then
  enable_parallel_IO=no
fi
if test x"${FC_ACT}" = x"frt" -o x"${FC_ACT}" = x"mpifrt"; then
  enable_parallel_IO=no
fi
if test x"${enable_single}" = x"yes"; then
  enable_parallel_IO=no
fi
if test x"${enable_mixed}" = x"yes"; then
  enable_parallel_IO=no
fi
if test x"${enable_gpu}" = x"yes"; then
  enable_parallel_IO=no
fi

AM_CONDITIONAL(PARALLEL_IO, test x"${enable_parallel_IO}" = x"yes")


## huge molecule on hard disk
#

AC_ARG_ENABLE([hmdisk],
    [AS_HELP_STRING(
      [--enable-hmdisk],
      [enable huge molecule on hard disk (default = no)])],
    [],[enable_hmdisk=no])

if test x"${enable_hmdisk}" = x"yes"; then
  AC_DEFINE(HM_DISK, 1, [defined if HM_DISK is used.])
fi


## CUDA
#
OPT_FC_CUDA=""

AC_ARG_WITH([cuda],
	[AS_HELP_STRING([--with-cuda@<:@=PATH@:>@],
	[use cuda (default: no)])],
	[], [with_cuda=no])

#AC_ARG_VAR([CUDA_LIB_PATH],
#    [CUDA library path (linker flags, e.g. -L<cuda lib>)])

if test x"${enable_gpu}" = x"yes"; then
  if test x"${with_cuda}" = x"no"; then
    with_cuda=yes
  fi
fi

if test x"${with_cuda}" != x"no"; then
  cuda_lib_path=""
  cuda_path=""
  if test x"${with_msmpi}" = x"yes"; then
    AC_MSG_ERROR([GPU is not available with MSMPI])
  fi
  if test x"${with_cuda}" != x"yes"; then
    cuda_path=${with_cuda}
    AC_CHECK_FILE([${cuda_path}/libcudart.so],
	     [cuda_lib_path=${cuda_path}],[])

    AC_CHECK_FILE([${cuda_path}/lib64/libcudart.so],
	     [cuda_lib_path=${cuda_path}/lib64],[])
  fi

  NVCC=no
  if test x"${with_cuda}" != x"yes"; then
    AC_CHECK_FILE([${cuda_path}/bin/nvcc],
    	   [NVCC=${cuda_path}/bin/nvcc],[])
  fi

  AC_PATH_PROG(NVCC, nvcc ,no)
  if test x"${NVCC}" = x"no"; then
    AC_MSG_ERROR([nvcc is not found])
  fi

  if test -z "${cuda_lib_path}"; then
    lib_dirs=$(echo ${LIBRARY_PATH} ${LD_LIBRARY_PATH} |tr ':' ' ')
    for dir in ${lib_dirs}; do
    if echo ${dir} | egrep -q "cuda"; then
      AC_CHECK_FILE([${dir}/libcudart.so], [cuda_lib_path=${dir}], [])
      if test x"${cuda_lib_path}" != x""; then
        break
      fi
    fi
    done
  fi

  if test -z "${cuda_lib_path}"; then
    # for debian & ubuntsu (non nvidia-package)
    AC_CHECK_FILE([/usr/lib/x86_64-linux-gnu/libcudart.so],
		  [cuda_lib_path=/usr/lib/x86_64-linux-gnu], [])
  fi

  if test -z "${cuda_lib_path}"; then
    # for fedora
    AC_CHECK_FILE([/usr/local/cuda/lib64/libcudart.so],
		  [cuda_lib_path=/usr/local/cuda/lib64], [])
  fi

  if test -z "${cuda_lib_path}"; then
    AC_CHECK_FILE([/usr/lib64/libcudart.so],
		  [cuda_lib_path=/usr/lib64], [])
  fi

  if test -z ${cuda_lib_path}; then
    AC_MSG_ERROR([CUDA is not found, set --with-cuda or CUDA_LIB_PATH])
  fi

#    AC_CHECK_PROG(NVCC, nvcc)
  version=""
  nvcc_ver_info=$(${NVCC} -V)
  for vtok in ${nvcc_ver_info} ; do
    if test x"${vtok}" = x"(NVCC)"; then
      check_ver=yes
    elif echo ${vtok} | egrep -q "V@<:@0-9@:>@*\.@<:@0-9@:>@*\.@<:@0-9@:>@*"; then
      version=${vtok}
    fi
  done

  AC_MSG_RESULT([$version])

  NVCC_VER=`echo ${version}\
	   | sed -e 's/\V\(@<:@0-9@:>@*\).\(@<:@0-9@:>@*\).\(@<:@0-9@:>@*\)/\1 \2 \3/'\
	   | xargs printf "%02d%02d%02d"`


  AC_DEFINE(CUDAGPU, 1, [defined if cuda_gpu is used.])
  DEFINED_VARIABLES+=" -DCUDAGPU"

  if test x"${FC_ACT}" = x"pgf90"; then
    OPT_LD_CUDA+="-Mcuda -L${cuda_lib_path} -lcudart -lstdc++"
    OPT_FC_CUDA+="-Mcuda "
    AC_DEFINE(PGICUDA, 1, [defined if pgi and cuda are used.])
    DEFINED_VARIABLES+=" -DPGICUDA"

  elif test x"${FC_ACT}" = x"ifort"; then
    OPT_LD_CUDA+="-L${cuda_lib_path} -lcudart -lstdc++"

  elif test x"${FC_ACT}" = x"gfortran"; then
    OPT_LD_CUDA+="-L${cuda_lib_path} -lcudart -lstdc++"

  else
    AC_MSG_ERROR([CUDA is not allowed])
  fi

# CUDA VERSION and compile option
  if test ${NVCC_VER} -lt 080000; then
    AC_MSG_ERROR([CUDA version must be >= 8.0.0 for GENESIS])
  fi
  if test ${NVCC_VER} -lt 120000; then
    GENCODEFLAG='--generate-code=arch=compute_35,code=\"sm_35,compute_35\"'
    if test ${NVCC_VER} -gt 065000; then
      GENCODEFLAG+=' --generate-code=arch=compute_37,code=\"sm_37,compute_37\"' #K80
    fi
  fi
  if test ${NVCC_VER} -gt 080000; then
    GENCODEFLAG+=' --generate-code=arch=compute_60,code=\"sm_60,compute_60\"' #P100
    GENCODEFLAG+=' --generate-code=arch=compute_61,code=\"sm_61,compute_61\"' #GTX1080
    GENCODEFLAG+=' --generate-code=arch=compute_62,code=\"sm_62,compute_62\"'
  fi

# NVCCFLAG
  NVCCFLAG="-c -g -O3 ${GENCODEFLAG} --ftz=false --prec-div=true --prec-sqrt=true --fmad=true"
  AC_SUBST(NVCCFLAG)
fi

AM_CONDITIONAL(CUDAGPU, test x"${with_cuda}" != x"no")

AC_ARG_ENABLE([package],
    [AS_HELP_STRING(
      [--enable-package],
      [enable package (default = no)])],
    [],[enable_package=no])


## ATDYN calculation
#

enable_atdyn=yes
if test x"${FC_ACT}" = x"pgif90" -o \
        x"${FC_ACT}" = x"pgf90" -o \
        x"${FC_ACT}" = x"pgfortran"; then
  enable_atdyn=no
fi
if test x"${enable_single}" = x"yes"; then
  enable_atdyn=no
fi
if test x"${enable_mixed}" = x"yes"; then
  enable_atdyn=no
fi
if test x"${enable_gpu}" = x"yes"; then
  enable_atdyn=no
fi
AM_CONDITIONAL(ATDYN, test x"${enable_atdyn}" = x"yes")


## ANALSIS calculation
#

enable_analysis=yes
if test x"${enable_single}" = x"yes"; then
  enable_analysis=no
fi
if test x"${enable_mixed}" = x"yes"; then
  enable_analysis=no
fi
if test x"${enable_gpu}" = x"yes"; then
  enable_analysis=no
fi
if test x"${FC_ACT}" = x"pgif90" -o \
        x"${FC_ACT}" = x"pgf90" -o \
        x"${FC_ACT}" = x"pgfortran"; then
  enable_analysis=no
fi
AM_CONDITIONAL(ANALYSIS, test x"${enable_analysis}" = x"yes")


## gprof
#
AC_ARG_ENABLE([gprof],
    [AS_HELP_STRING(
      [--enable-gprof],
      [enable gprof profiler (default = no)])],
    [],[enable_gprof=no])

OPT_GPROF=""
if test x"${enable_gprof}" = x"yes"; then
  OPT_GPROF="-pg -g"
fi


## dSFMT settings
#

AC_DEFINE(DSFMT_MEXP, 19937, [defined always.])
DEFINED_VARIABLES+=" -DDSFMT_MEXP=19937"


# architecture check
#
cpuinfo=`grep "^flags.*" /proc/cpuinfo | head -1 `
have_avx512f=no
have_avx512er=no
have_avx512cd=no
have_avx512pf=no
have_avx512dq=no
have_avx512bw=no
have_avx512vl=no
have_avx2=no
have_avx=no
have_sse4_2=no
have_rtm=no
have_rdseed=no
have_mmxext=no
have_fxsr_opt=no
have_3dnowprefetch=no
for vtok in ${cpuinfo} ; do
  if test x"${vtok}" = x"avx512f"; then
	  have_avx512f=yes
  elif test x"${vtok}" = x"avx512er"; then
	  have_avx512er=yes
  elif test x"${vtok}" = x"avx512cd"; then
	  have_avx512cd=yes
  elif test x"${vtok}" = x"avx512pf"; then
	  have_avx512pf=yes
  elif test x"${vtok}" = x"avx512dq"; then
	  have_avx512dq=yes
  elif test x"${vtok}" = x"avx512bw"; then
	  have_avx512bw=yes
  elif test x"${vtok}" = x"avx512vl"; then
	  have_avx512vl=yes
  elif test x"${vtok}" = x"avx2"; then
	  have_avx2=yes
  elif test x"${vtok}" = x"avx"; then
	  have_avx=yes
  elif test x"${vtok}" = x"sse4_2"; then
	  have_sse4_2=yes
  elif test x"${vtok}" = x"rtm"; then
	  have_rtm=yes
  elif test x"${vtok}" = x"rdseed"; then
	  have_rdseed=yes
  elif test x"${vtok}" = x"mmxext"; then
      have_mmxext=yes
  elif test x"${vtok}" = x"3dnowprefetch"; then
      have_3dnowprefetch=yes
  elif test x"${vtok}" = x"fxsr_opt"; then
      have_fxsr_opt=yes
  fi
done


## compile flags
#

AC_ARG_WITH([simd],
	[AS_HELP_STRING([--with-simd=@<:@auto/MIC-AVX512/CORE-AVX512/CORE-AVX2@:>@],
	[SIMD instruction (default: auto)])],
	[], [with_simd=auto])

if test x"${FC_ACT}" != x"ifort" -a x"${with_simd}" != x"auto" ; then
  AC_MSG_ERROR([with-simd is only used in intel fortran])
fi

if test x"${FC_ACT}" = x"ifort" ; then
  if test ${FC_GNU_VER} -lt 152000; then

    if test x"${with_simd}" = x"CORE-AVX512"; then
      AC_MSG_ERROR([ifort > 15.2 is required for CORE-AVX512])
    fi

    if test x"${with_simd}" = x"MIC-AVX512"; then
      AC_MSG_ERROR([ifort > 15.2 is required for MIC-AVX512])
    fi
  fi
fi

FARCFLAGS=""
FARCFLAGS_NOOPT=""
FARCFLAGS_KERNELS=""
FARCFLAGS_KERNELS_PAIRLIST=""
FARCFLAGS_DIHEDRAL=""

# FFLAGS
if test x"${ac_test_FFLAGS}" = x"set"; then
  FCFLAGS=${ac_save_FFLAGS}

else

  if test x"${enable_debug}" == x"0"; then

# no-debug no report
    if test x"${FC_ACT}" = x"ifort"; then
      if test x"${have_fxsr_opt}"  = x"yes" -o \
              x"${have_mmxext}" = x"yes" -o \
              x"${have_3dnowprefetch}" = x"yes"; then
#AMD
        FARCFLAGS=""
        FARCFLAGS_NOOPT=""
        FCFLAGS+=" -O3 "
      else
        FARCFLAGS=" -xHost "
        FARCFLAGS_NOOPT=" -axCORE-AVX2 -xHost "
        FCFLAGS+=" -O3 -ip "
      fi
      if test x"${with_simd}" != x"auto"; then
        FARCFLAGS+=" -ax${with_simd}"
      elif test ${FC_GNU_VER} -ge 152000 ; then
        if test x"${have_avx512f}"  = x"yes" -a \
                x"${have_avx512er}" = x"yes" -a \
                x"${have_avx512cd}" = x"yes" -a \
                x"${have_avx512pf}" = x"yes"; then
          with_simd="MIC-AVX512"
          FARCFLAGS+=" -ax${with_simd}"
        elif test x"${have_avx512f}"  = x"yes" -a \
                  x"${have_avx512dq}" = x"yes" -a \
                  x"${have_avx512cd}" = x"yes" -a \
                  x"${have_avx512vl}" = x"yes" -a \
                  x"${have_avx512bw}" = x"yes"; then
          with_simd="CORE-AVX512"
          FARCFLAGS+=" -ax${with_simd}"
        elif test x"${have_avx2}"  = x"yes" -a \
                  x"${have_rdseed}" = x"yes" -a \
                  x"${have_rtm}" = x"yes" ; then
          with_simd="CORE-AVX512"
          FARCFLAGS+=" -ax${with_simd}"
        elif test x"${have_avx2}"  = x"yes" ; then
          with_simd="CORE-AVX2"
          FARCFLAGS+=" -ax${with_simd}"
        fi
      fi

	  if test ${FC_GNU_VER} -ge 20210000; then
        FCFLAGS+=" -qmkl=parallel "
      else
        FCFLAGS+=" -mkl=parallel "
      fi

    elif test x"${FC_ACT}" = x"gfortran"; then
      FCFLAGS="-O3 -ffast-math -march=native -ffree-line-length-none"
      echo ${FC_GNU_VER}
      if test ${FC_GNU_VER} -gt 100000; then
        if test x"${with_msmpi}" = x"no"; then
          FCFLAGS+=" -fallow-argument-mismatch"
        else
          FCFLAGS+=" -Wno-argument-mismatch"
        fi
      fi

    elif test x"${FC_ACT}" = x"pgf90"; then
      FCFLAGS="-O2 -fastsse -Mipa=fast,inline -m64 -Mflushz -pc 64 -mcmodel=medium -Msmartalloc=huge"

    elif test x"${FC_ACT}" = x"frtpx"; then
      FCFLAGS="-m -Kocl -Kfast,openmp,parallel,simd=2,noalias -Nlst=t -Koptmsg=2"

    elif test x"${FC_ACT}" = x"mpifrtpx"; then
      FCFLAGS="-m -Kocl -Kfast,openmp,parallel,simd=2,noalias -Nlst=t -Koptmsg=2"
      FCFLAGS_FFTE_KERNEL_FPP__PGI__SP=" -m -Kocl  -Kfast,openmp,parallel,simd=2,noalias -Nlst=t -Koptmsg=2"
      FCFLAGS_FPP__PGI_NOOPT="-m -Kocl -Kfast,openmp,parallel,simd=2,noalias -Nlst=t -Koptmsg=2"
      FCFLAGS_KERNELS_PAIRLIST="-m -Kocl  -Kfast,openmp,parallel,simd=2,nounroll,swp_strong,noalias=s -Nlst=t -Koptmsg=2"
      FCFLAGS_DIHEDRAL="-m  -xsp_energy_dihedrals_mod.calculate_dihedral_2 -Kocl -Kfast,openmp,parallel,simd=2,noalias -Nlst=t -Koptmsg=2"
      FCFLAGS_KERNELS=" -m -Kocl  -Kfast,openmp,parallel,simd=2,swp_strong,noalias=s -Nlst=t -Koptmsg=2"

    elif test x"${FC_ACT}" = x"frt"; then
      FCFLAGS="-m -Am -Kfast"

    elif test x"${FC_ACT}" = x"mpifrt"; then
      FCFLAGS="-m -Am -Kfast"

    fi

  elif test x"${enable_debug}" == x"1"; then

# debug level 1 (without higher optmization ) (default)
    if   test x"${FC_ACT}" = x"ifort"; then
      FCFLAGS=" -axCORE-AVX2 -vec-report"
	  if test ${FC_GNU_VER} -ge 20210000; then
        FCFLAGS+=" -qmkl=parallel "
      else
        FCFLAGS+=" -mkl=parallel "
      fi

    elif test x"${FC_ACT}" = x"gfortran"; then
      FCFLAGS="-Wall -ffree-line-length-none"
      if test ${FC_GNU_VER} -gt 100000; then
        if test x"${with_msmpi}" = x"no"; then
          FCFLAGS+=" -fallow-argument-mismatch"
        else
          FCFLAGS+=" -O2 -march=x86-64 -Wno-argument-mismatch"
        fi
      fi

    elif test x"${FC_ACT}" = x"pgf90"; then
      FCFLAGS="-m64 -Mflushz -pc 64 -mcmodel=medium -Msmartalloc=huge -Mbackslash"

    elif test x"${FC_ACT}" = x"frtpx"; then
      FCFLAGS="-m -Kocl -Kopenmp -Nlst=t  -Koptmsg=2"

    elif test x"${FC_ACT}" = x"mpifrtpx"; then
      FCFLAGS="-m -Kocl -Kopenmp -Nlst=t  -Koptmsg=2"

    elif test x"${FC_ACT}" = x"frt"; then
      FCFLAGS="-m "

    elif test x"${FC_ACT}" = x"mpifrt"; then
      FCFLAGS="-m "

    fi

  elif test x"${enable_debug}" == x"2"; then

# debug level 2 (with g)
    AC_DEFINE(DEBUG, 1, [defined if Debug is used.])
    DEFINED_VARIABLES+=" -DDEBUG"

    if   test x"${FC_ACT}" = x"ifort"; then
      FCFLAGS="-g -axCORE-AVX2 -vec-report"
	  if test ${FC_GNU_VER} -ge 20210000; then
        FCFLAGS+=" -qmkl=parallel "
      else
        FCFLAGS+=" -mkl=parallel "
      fi

    elif test x"${FC_ACT}" = x"gfortran"; then
      FCFLAGS="-g -Wall -ffree-line-length-none"
      if test ${FC_GNU_VER} -gt 100000; then
        if test x"${with_msmpi}" = x"no"; then
          FCFLAGS+=" -fallow-argument-mismatch"
        else
          FCFLAGS+=" -Wno-argument-mismatch"
        fi
      fi

    elif test x"${FC_ACT}" = x"pgf90"; then
      FCFLAGS="-g  -m64 -Mflushz -pc 64 -mcmodel=medium -Msmartalloc=huge -Mbackslash"

    elif test x"${FC_ACT}" = x"frtpx"; then
      FCFLAGS="-g -m -Kocl -Kopenmp"

    elif test x"${FC_ACT}" = x"mpifrtpx"; then
      FCFLAGS=" -g -m -Kocl -Kopenmp  -Nlst=t -Koptmsg=2"

    elif test x"${FC_ACT}" = x"frt"; then
      FCFLAGS="-g -m "

    elif test x"${FC_ACT}" = x"mpifrt"; then
      FCFLAGS="-g -m "

    fi

  elif test x"${enable_debug}" == x"3"; then

# debug level 3 (with g and traceback and without optimization)
    AC_DEFINE(DEBUG, 1, [defined if Debug is used.])
    DEFINED_VARIABLES+=" -DDEBUG"

    if   test x"${FC_ACT}" = x"ifort"; then
      FCFLAGS="-O0 -g -traceback -CB"

    elif test x"${FC_ACT}" = x"gfortran"; then
      FCFLAGS="-O0 -g -Wall -fbounds-check -ffree-line-length-none -fbacktrace"
      if test ${FC_GNU_VER} -gt 100000; then
        if test x"${with_msmpi}" = x"no"; then
          FCFLAGS+=" -fallow-argument-mismatch"
        else
          FCFLAGS+=" -Wno-argument-mismatch"
        fi
      fi

    elif test x"${FC_ACT}" = x"pgf90"; then
      FCFLAGS="-O0 -g -m64 -Msmartalloc=huge"

    elif test x"${FC_ACT}" = x"frtpx"; then
      FCFLAGS=" -O0 -g -m -Kocl -Kopenmp -Hs -Hu -Ha -Nlst=t -Koptmsg=2"

    elif test x"${FC_ACT}" = x"mpifrtpx"; then
      FCFLAGS="-O0 -g -Hs -Hu -Ha"
      FCFLAGS_KERNELS="-O0 -g -Hs -Hu -Ha -Kopenmp"
      FCFLAGS_KERNELS_PAIRLIST="-O0 -g -Hs -Hu -Ha -Kopenmp"
      FCFLAGS_FFTE_KERNEL_FPP__PGI__SP="-O0 -g -Hs -Hu -Ha -Kopenmp"
      FCFLAGS_FPP__PGI_NOOPT="-O0 -g -Hs -Hu -Ha -Kopenmp"
      FCFLAGS_DIHEDRAL="-O0 -g -Hs -Hu -Ha -Kopenmp"

    elif test x"${FC_ACT}" = x"frt"; then
      FCFLAGS="-O0 -g"

    elif test x"${FC_ACT}" = x"mpifrt"; then
      FCFLAGS="-O0 -g"

    fi
  else

# debug level 4 (with g and traceback and fully check (only for intel))
    AC_DEFINE(DEBUG, 1, [defined if Debug is used.])
    DEFINED_VARIABLES+=" -DDEBUG"

    if   test x"${FC_ACT}" = x"ifort"; then
      FCFLAGS="-O0 -g -traceback -check all,noarg_temp_created"

    elif   test x"${FC_ACT}" = x"gfortran"; then
      FCFLAGS="-O0 -g -fbacktrace -fcheck=all,no-array-temps  -ffree-line-length-none -ffpe-trap=invalid,zero,overflow,underflow"
      if test ${FC_GNU_VER} -gt 100000; then
        if test x"${with_msmpi}" = x"no"; then
          FCFLAGS+=" -fallow-argument-mismatch"
        else
          FCFLAGS+=" -Wno-argument-mismatch"
        fi
      fi

    else
      AC_MSG_ERROR([DEBUG=4 is only allowed in IFORT/GFORTRAN])

    fi

  fi

fi

if test x"${FC_ACT}" = x"ifort"; then
  FCFLAGS+=" -assume byterecl"
fi

FCFLAGS+=" ${OPT_OPENMP}"
FCFLAGS+=" ${OPT_FC_CUDA}"
FCFLAGS+=" ${OPT_GPROF}"
TMP="${FARCFLAGS} ${FCFLAGS}"
FCFLAGS_NOOPT="${FARCFLAGS_NOOPT} ${FCFLAGS}"
FCFLAGS_DISABLEFARC=" ${FCFLAGS}"
FCFLAGS=${TMP}

AC_SUBST(FCFLAGS_NOOPT)
AC_SUBST(FCFLAGS_KERNELS)
AC_SUBST(FCFLAGS_KERNELS_PAIRLIST)
AC_SUBST(FCFLAGS_DIHEDRAL)
AC_SUBST(FCFLAGS_DISABLEFARC)

# CFLAGS
if test x"${ac_test_CFLAGS}" = x"set"; then
  CFLAGS=${ac_save_CFLAGS}

else

  if test x"${enable_debug}" == x"0"; then

# no-debug no report
    if   test x"${CC_ACT}" = x"icc"; then
      CFLAGS="-O3 -ip -axCORE-AVX2 "

    elif test x"${CC_ACT}" = x"gcc"; then
      CFLAGS="-O3 -ffast-math -march=native"
      if test ${FC_GNU_VER} -gt 100000; then
        if test x"${with_msmpi}" = x"no"; then
          CFLAGS+=" -fallow-argument-mismatch"
        else
          CFLAGS+=" -Wno-argument-mismatch"
        fi
      fi

    elif test x"${CC_ACT}" = x"pgcc"; then
      CFLAGS="-O3 -fastsse -m64 -Mflushz -pc 64 -mcmodel=medium -Msmartalloc=huge"

    elif test x"${CC_ACT}" = x"fccpx"; then
      CFLAGS="-m -Kvisimpact -Kocl -Kswp"

    elif test x"${CC_ACT}" = x"mpifccpx"; then
      CFLAGS="-m -Kvisimpact -Kocl -Kswp"

    elif test x"${CC_ACT}" = x"fcc"; then
      CFLAGS="-m -Kfast"

    elif test x"${CC_ACT}" = x"mpifcc"; then
      CFLAGS="-m -Kfast"

    fi

  elif test x"${enable_debug}" == x"1"; then

# debug level 1 (without O3) (default)
    if   test x"${CC_ACT}" = x"icc"; then
      CFLAGS="-ip -axCORE-AVX2 -vec-report"

    elif test x"${CC_ACT}" = x"gcc"; then
      CFLAGS="-Wall"
      if test ${FC_GNU_VER} -gt 100000; then
        if test x"${with_msmpi}" = x"no"; then
          CFLAGS+=" -fallow-argument-mismatch"
        else
          CFLAGS+=" -Wno-argument-mismatch"
        fi
      fi

    elif test x"${CC_ACT}" = x"pgcc"; then
      CFLAGS="-fastsse -m64 -Mflushz -pc 64 -mcmodel=medium -Msmartalloc=huge"

    elif test x"${CC_ACT}" = x"fccpx"; then
      CFLAGS="-m -Kocl"

    elif test x"${CC_ACT}" = x"mpifccpx"; then
      CFLAGS="-m -Kocl"

    elif test x"${CC_ACT}" = x"fcc"; then
      CFLAGS="-m -Kfast"

    elif test x"${CC_ACT}" = x"mpifcc"; then
      CFLAGS="-m -Kfast"

    fi

  elif test x"${enable_debug}" == x"2"; then

# debug level 2 (with g)
    AC_DEFINE(DEBUG, 1, [defined if Debug is used.])
    DEFINED_VARIABLES+=" -DDEBUG"

    if   test x"${CC_ACT}" = x"icc"; then
      CFLAGS="-g -ip -axCORE-AVX2 -vec-report"

    elif test x"${CC_ACT}" = x"gcc"; then
      CFLAGS="-g -Wall"
      if test ${FC_GNU_VER} -gt 100000; then
        if test x"${with_msmpi}" = x"no"; then
          CFLAGS+=" -fallow-argument-mismatch"
        else
          CFLAGS+=" -Wno-argument-mismatch"
        fi
      fi

    elif test x"${CC_ACT}" = x"pgcc"; then
      CFLAGS="-g -fastsse -m64 -Mflushz -pc 64 -mcmodel=medium -Msmartalloc=huge"

    elif test x"${CC_ACT}" = x"fccpx"; then
      CFLAGS="-g -m -Kocl"

    elif test x"${CC_ACT}" = x"mpifccpx"; then
      CFLAGS="-g -m -Kocl"

    elif test x"${CC_ACT}" = x"fcc"; then
      CFLAGS="-g -m -Kfast"

    elif test x"${CC_ACT}" = x"mpifcc"; then
      CFLAGS="-g -m -Kfast"

    fi

  elif test x"${enable_debug}" == x"3"; then

# debug level 3 (with g and traceback and without optimization)
    AC_DEFINE(DEBUG, 1, [defined if Debug is used.])
    DEFINED_VARIABLES+=" -DDEBUG"

    if   test x"${CC_ACT}" = x"icc"; then
      CFLAGS="-O0 -g -traceback -CB"

    elif test x"${CC_ACT}" = x"gcc"; then
      CFLAGS="-O0 -g -Wall -fbounds-check"
      if test ${FC_GNU_VER} -gt 100000; then
        if test x"${with_msmpi}" = x"no"; then
          CFLAGS+=" -fallow-argument-mismatch"
        else
          CFLAGS+=" -Wno-argument-mismatch"
        fi
      fi

    elif test x"${CC_ACT}" = x"pgcc"; then
      CFLAGS="-O0 -g -m64 -Msmartalloc=huge"

    elif test x"${CC_ACT}" = x"fccpx"; then
      CFLAGS="-O0 -g -Kocl"

    elif test x"${CC_ACT}" = x"mpifccpx"; then
      CFLAGS="-O0 -g -Kocl"

    elif test x"${CC_ACT}" = x"fcc"; then
      CFLAGS="-O0 -g"

    elif test x"${CC_ACT}" = x"mpifcc"; then
      CFLAGS="-O0 -g"

    fi
  else

# debug level 4 (with g and traceback and fully check (only for intel))
    AC_DEFINE(DEBUG, 1, [defined if Debug is used.])
    DEFINED_VARIABLES+=" -DDEBUG"

    if test x"${FC_ACT}" = x"ifort"; then
      CFLAGS="-O0 -g -traceback -check all,noarg_temp_created"

    elif test x"${FC_ACT}" = x"gfortran"; then
      CFLAGS="-O0 -g -fbacktrace -fcheck=all,no-array-temps  -ffree-line-length-none -ffpe-trap=invalid,zero,overflow,underflow"
      if test ${FC_GNU_VER} -gt 100000; then
        if test x"${with_msmpi}" = x"no"; then
          CFLAGS+=" -fallow-argument-mismatch"
        else
          CFLAGS+=" -Wno-argument-mismatch"
        fi
      fi

    else
      AC_MSG_ERROR([DEBUG=4 is only allowed in IFORT/GFORTRAN])

    fi

  fi

fi

CFLAGS+=" ${OPT_OPENMP}"
CFLAGS+=" ${OPT_GPROF}"

# LDFLAGS
if test -z "${LDFLAGS}"; then
  if   test x"${FC_ACT}" = x"ifort"; then
    LDFLAGS=""

  elif test x"${FC_ACT}" = x"pgf90"; then
    LDFLAGS="-Msmartalloc=huge"

  elif test x"${FC_ACT}" = x"gfortran"; then
    LDFLAGS=""
    #LDFLAGS="-march=native -fuse-linker-plugin"

  elif test x"${FC_ACT}" = x"frtpx"; then
    LDFLAGS="-SSL2 -Kparallel -Nfjomplib"

  elif test x"${FC_ACT}" = x"mpifrtpx"; then
    LDFLAGS="-SSL2BLAMP -Kparallel -Nfjomplib"

  elif test x"${FC_ACT}" = x"frt"; then
    LDFLAGS="-SSL2 -Kparallel -Nfjomplib"

  elif test x"${FC_ACT}" = x"mpifrt"; then
    LDFLAGS="-SSL2BLAMP -Kparallel -Nfjomplib"

  fi
fi

if test x"${FC_ACT}" = x"ifort"; then
  LDFLAGS+=" -assume byterecl -lgfortran"
fi

LDFLAGS+=" ${OPT_OPENMP}"
if test x"${with_msmpi}" = x"yes"; then
  LDFLAGS+=" -lmsmpi -static"
fi
LDFLAGS+=" ${OPT_LD_LAPACK}"
LDFLAGS+=" ${OPT_LD_CUDA}"
LDFLAGS+=" ${OPT_GPROF}"

# PPFLAGS
if test -z "${PPFLAGS}"; then

  if   test x"${FPP}" = x"frtpx"; then
    PPFLAGS="-P -Cpp -Cfpp"

  elif test x"${FPP}" = x"frt"; then
    PPFLAGS="-P -Cpp -Cfpp"

  elif test x"${FPP}" = x"cpp"; then
    PPFLAGS="-traditional-cpp -traditional"

  fi
fi
AC_SUBST(PPFLAGS)


## K computer settings
#

# USEKCOMP
AM_CONDITIONAL(USEKCOMP, test x"${FC_ACT}" = x"frtpx" -o \
			      x"${FC_ACT}" = x"mpifrtpx" -o \
                              x"${FC_ACT}" = x"frt" -o \
			      x"${FC_ACT}" = x"mpifrt")

# KCOMP preprocess macro
AC_ARG_ENABLE([pktimer],
    [AS_HELP_STRING(
      [--enable-pktimer],
      [enable enhanced pktimer (default = no)])],
    [],[enable_pktimer=no])

AC_ARG_WITH([fapp],
	[AS_HELP_STRING([--with-fapp],
	[use fapp (default: no)])],
	[], [with_fapp=no])

AC_ARG_WITH([fipp],
	[AS_HELP_STRING([--with-fipp],
	[use fipp (default: no)])],
	[], [with_fipp=no])

AC_ARG_WITH([fj_timer_detail],
	[AS_HELP_STRING([--with-fj_timer_detail],
	[use FJ-timer in detail (default: not use)])],
	[], [with_fj_timer_detail=no])

AC_ARG_WITH([fj_timer_2],
	[AS_HELP_STRING([--with-fj_timer_2],
	[use FJ-timer in main loop (default: not use)])],
	[], [with_fj_timer_2=no])

if test x"${FC_ACT}" = x"frtpx" -o x"${FC_ACT}" = x"mpifrtpx"; then
  AC_DEFINE(KCOMP, 1, [defined if K-computer compiler is used.])
  DEFINED_VARIABLES+=" -DKCOMP"

  if test x"${host_alias}" = x"Fugaku"; then
    AC_DEFINE(FUGAKU, 1, [defined if --host=Fugaku is set.])
  fi
## Fujitsu timer
#
  if test x"${with_fapp}" = x"yes" -o \
          x"${with_fipp}" = x"yes" -o \
          x"${with_fj_timer_detail}" = x"yes" -o \
          x"${with_fj_timer_2}" = x"yes" ; then
     enable_pktimer=yes
  fi
  mainly_pktimer=no
  if test x"${enable_pktimer}" = x"yes"; then

    if test x"${FC_ACT}" = x"mpifrtpx"; then
      AC_DEFINE(PKTIMER, 1, [defined if PKTIMER is used.])
      DEFINED_VARIABLES+=" -DPKTIMER"
      mainly_pktimer=yes
    fi

    if test x"${with_fapp}" = x"yes"; then
      AC_DEFINE(FJ_PROF_FAPP, 1, [defined if fapp is used.])
      DEFINED_VARIABLES+=" -DFJ_PROF_FAPP"
    fi

    if test x"${with_fipp}" = x"yes"; then
      AC_DEFINE(FJ_PROF_FIPP, 1, [defined if fipp is used.])
      DEFINED_VARIABLES+=" -DFJ_PROF_FIPP"
    fi

    if test x"${with_fj_timer_detail}" = x"yes"; then
      AC_DEFINE(FJ_TIMER_DETAIL, 1, [defined if fipp is used.])
      DEFINED_VARIABLES+=" -DFJ_TIMER_DETAIL"
    fi

    if test x"${with_fj_timer_2}" = x"yes"; then
      AC_DEFINE(FJ_TIMER_2, 1, [defined if fipp is used.])
      DEFINED_VARIABLES+=" -DFJ_TIMER_2"
    fi

  fi

fi

if test x"${FC_ACT}" = x"frt" -o x"${FC_ACT}" = x"mpifrt"; then
  AC_DEFINE(KCOMP, 1, [defined if K-computer compiler is used.])
  DEFINED_VARIABLES+=" -DKCOMP"

  if test x"${host_alias}" = x"Fugaku"; then
    AC_DEFINE(FUGAKU, 1, [defined if --host=Fugaku is set.])
  fi
## Fujitsu timer
#
  if test x"${with_fapp}" = x"yes" -o \
          x"${with_fipp}" = x"yes" -o \
          x"${with_fj_timer_detail}" = x"yes" -o \
          x"${with_fj_timer_2}" = x"yes" ; then
     enable_pktimer=yes
  fi
  mainly_pktimer=no
  if test x"${enable_pktimer}" = x"yes"; then

    if test x"${FC_ACT}" = x"mpifrt"; then
      AC_DEFINE(PKTIMER, 1, [defined if PKTIMER is used.])
      DEFINED_VARIABLES+=" -DPKTIMER"
      mainly_pktimer=yes
    fi

    if test x"${with_fapp}" = x"yes"; then
      AC_DEFINE(FJ_PROF_FAPP, 1, [defined if fapp is used.])
      DEFINED_VARIABLES+=" -DFJ_PROF_FAPP"
    fi

    if test x"${with_fipp}" = x"yes"; then
      AC_DEFINE(FJ_PROF_FIPP, 1, [defined if fipp is used.])
      DEFINED_VARIABLES+=" -DFJ_PROF_FIPP"
    fi

    if test x"${with_fj_timer_detail}" = x"yes"; then
      AC_DEFINE(FJ_TIMER_DETAIL, 1, [defined if fipp is used.])
      DEFINED_VARIABLES+=" -DFJ_TIMER_DETAIL"
    fi

    if test x"${with_fj_timer_2}" = x"yes"; then
      AC_DEFINE(FJ_TIMER_2, 1, [defined if fipp is used.])
      DEFINED_VARIABLES+=" -DFJ_TIMER_2"
    fi

  fi

fi

## Intel compiler settings
#

# INTEL preprocess macro
if test x"${FC_ACT}" = x"ifort"; then
  AC_DEFINE(INTEL, 1, [defined if Intel compiler is used.])
  DEFINED_VARIABLES+=" -DINTEL"
fi

# GCC gfortran settings
if test x"${FC_ACT}" = x"gfortran"; then
  AC_DEFINE(__GFORTRAN__, 1, [defined if GCC gfortran compiler is used.])
  DEFINED_VARIABLES+=" -D__GFORTRAN__"
fi

ARFLAGS=""
AC_SUBST(ARFLAGS)

# PGI preprocess macro
if test x"${FC_ACT}" = x"pgf90"; then
  AC_DEFINE(__PGI, 1, [defined if pgi compiler is used.])
  DEFINED_VARIABLES+=" -D__PGI"
fi
AM_CONDITIONAL(PGICOMP, test x"${FC_ACT}" = x"pgf90")

# MS-MPI setting

if test x"${with_msmpi}" = x"yes"; then
  AC_DEFINE(MSMPI, 1, [defined if MS-MPI compiler is used.])
  DEFINED_VARIABLES+=" -DMSMPI"
fi

## PGI compiler settings
#

if test x"${mainly_single}" = x"yes" -a \
        x"${FC_ACT}"        = x"pgf90"; then

  if test x"${ac_test_FFLAGS}" = x"set"; then
    FCFLAGS_FFTE_KERNEL_FPP__PGI__SP=${ac_save_FFLAGS}

  else

    if test x"${enable_debug}" = x"0"; then
      FCFLAGS_FFTE_KERNEL_FPP__PGI__SP=" -m64 -Mflushz -pc 64 -mcmodel=medium -Msmartalloc=huge  -Mbackslash"

    elif test x"${enable_debug}" = x"1"; then
      FCFLAGS_FFTE_KERNEL_FPP__PGI__SP="-m64 -Mflushz -pc 64 -mcmodel=medium -Msmartalloc=huge -Mbackslash"

    elif test x"${enable_debug}" = x"2"; then
      FCFLAGS_FFTE_KERNEL_FPP__PGI__SP="-g -m64 -Mflushz -pc 64 -mcmodel=medium -Msmartalloc=huge -Mbackslash"

    elif test x"${enable_debug}" = x"3"; then
      FCFLAGS_FFTE_KERNEL_FPP__PGI__SP="-O0 -g -m64 -Msmartalloc=huge -Mbackslash"

    fi

  fi

else

  FCFLAGS_FFTE_KERNEL_FPP__PGI__SP=${FCFLAGS}

fi

AC_SUBST(FCFLAGS_FFTE_KERNEL_FPP__PGI__SP)

if test x"${FC_ACT}"        = x"pgf90"; then
  if test x"${ac_test_FFLAGS}" = x"set"; then
    FCFLAGS_FPP__PGI_NOOPT=${ac_save_FFLAGS}

  else

    if test x"${enable_debug}" = x"0"; then
      FCFLAGS_FPP__PGI_NOOPT="-m64 -Mflushz -pc 64 -mcmodel=medium -Msmartalloc=huge -Mbackslash"

    elif test x"${enable_debug}" = x"1"; then
      FCFLAGS_FPP__PGI_NOOPT="-m64 -Mflushz -pc 64 -mcmodel=medium -Msmartalloc=huge -Mbackslash"

    elif test x"${enable_debug}" = x"2"; then
      FCFLAGS_FPP__PGI_NOOPT="-g -m64 -Mflushz -pc 64 -mcmodel=medium -Msmartalloc=huge -Mbackslash"

    elif test x"${enable_debug}" = x"3"; then
      FCFLAGS_FPP__PGI_NOOPT="-O0 -g -m64 -Msmartalloc=huge -Mbackslash"

    fi

  fi
else
  FCFLAGS_FPP__PGI_NOOPT=${FCFLAGS}
fi
AC_SUBST(FCFLAGS_FPP__PGI_NOOPT)


## check compile information
#
test -d .git && rm -f src/GENESIS_VERSION
test -d .git && git checkout src/GENESIS_VERSION

AC_DEFINE_UNQUOTED([COMPILE_HOST],["`hostname`"],[hostname])
AC_DEFINE_UNQUOTED([COMPILE_USER],["`whoami`"],[username])
AC_DEFINE_UNQUOTED([COMPILE_GENESIS_VERSION],["`cat src/GENESIS_VERSION | sed -e 's/\$GENESIS_VERSION\$: //'`"],[genesis version])
AC_DEFINE_UNQUOTED([COMPILE_FC_VER],["`${FC_ACT} ${VER_OPT} 2>&1 | sed -e '/^$/d'| head -1`"],[fortran compiler version])
AC_DEFINE_UNQUOTED([COMPILE_CC_VER],["`${CC_ACT} ${VER_OPT} 2>&1 | sed -e '/^$/d'| head -1`"],[c compiler version])
AC_DEFINE_UNQUOTED([COMPILE_DEFINED_VARIABLES],["${DEFINED_VARIABLES}"],[defined variables])
AC_DEFINE_UNQUOTED([COMPILE_FCFLAGS],["${FCFLAGS}"],[fortran flags])
AC_DEFINE_UNQUOTED([COMPILE_CFLAGS],["${CFLAGS}"],[c flags])
AC_DEFINE_UNQUOTED([COMPILE_LDFLAGS],["${LDFLAGS}"],[ld flags])
AC_DEFINE_UNQUOTED([MPITYPE],["${mpitype}"],[MPI WRAPPER TYPE])
if test x"${with_cuda}" != x"no"; then
  AC_DEFINE_UNQUOTED([COMPILE_NVCC_VER],["`${NVCC} -V 2>&1 | tail -1`"],[cuda version])
fi


## QSimulate-QM/MM calculation
#
AC_ARG_ENABLE([qsimulate],
    [AS_HELP_STRING([--enable-qsimulate],[enable QSimulate integration.])]
	[], [enable_qsimulate=yes])
AM_CONDITIONAL(QSIMULATE, test x"${enable_qsimulate}" = x"yes")
if  test x"${enable_qsimulate}" = x"yes"; then
  AC_DEFINE(QSIMULATE, 1, [defined if QSimulate is used.])
  AC_LANG_PUSH(C++)
  AC_ARG_WITH([include], [AS_HELP_STRING([--with-include],[include flags])], [with_include=$withval], [with_include=no])
  if test "x${with_include}" != xno; then
    CXXFLAGS="${with_include} ${CXXFLAGS}"
    CPPFLAGS="${with_include} ${CPPFLAGS}"
  fi
  AC_DEFINE(HAVE_MPI_H, 1, [defined if MPI is used.]) # there has to be
  AC_CHECK_HEADERS([bagel.h], [], [AC_MSG_ERROR([bagel.h not found or not working, set the include path])], [])
  AC_CHECK_LIB(bagel, main,  [LDFLAGS+=" -lbagel"], [AC_MSG_ERROR("Linking against BAGEL library failed; perhaps add -L/path/to/bagel/lib in LDFLAGS")])
  ## because we changed the linked to CXX
  LDFLAGS+=" -lmpifort -lgfortran -lstdc++"

 if test x"${with_msmpi}" = x"yes"; then
    AC_MSG_ERROR([qsimulate is not available with MSMPI])
 fi
fi


##

AC_CONFIG_FILES([Makefile
                src/Makefile
                src/lib/Makefile
                src/lib/qsimulate/Makefile
                src/lib/json-fortran/Makefile])
AM_COND_IF([SPDYN],[AC_CONFIG_FILES([src/spdyn/Makefile])])
AM_COND_IF([ATDYN],[AC_CONFIG_FILES([src/atdyn/Makefile])])
AM_COND_IF([ANALYSIS],[AC_CONFIG_FILES([src/analysis/Makefile
                src/analysis/libana/Makefile
                src/analysis/trj_analysis/Makefile
                src/analysis/trj_analysis/comcrd_analysis/Makefile
                src/analysis/trj_analysis/msd_analysis/Makefile
                src/analysis/trj_analysis/diffusion_analysis/Makefile
                src/analysis/trj_analysis/distmat_analysis/Makefile
                src/analysis/trj_analysis/drms_analysis/Makefile
                src/analysis/trj_analysis/fret_analysis/Makefile
                src/analysis/trj_analysis/hb_analysis/Makefile
                src/analysis/trj_analysis/lipidthick_analysis/Makefile
                src/analysis/trj_analysis/qval_analysis/Makefile
                src/analysis/trj_analysis/qval_residcg_analysis/Makefile
                src/analysis/trj_analysis/rg_analysis/Makefile
                src/analysis/trj_analysis/rmsd_analysis/Makefile
                src/analysis/trj_analysis/tilt_analysis/Makefile
                src/analysis/trj_analysis/trj_analysis/Makefile
                src/analysis/trj_analysis/ring_analysis/Makefile
                src/analysis/trj_analysis/energy_analysis/Makefile
                src/analysis/sp_analysis/Makefile
                src/analysis/sp_analysis/libspana/Makefile
                src/analysis/sp_analysis/rdf_analysis/Makefile
                src/analysis/sp_analysis/density_analysis/Makefile
                src/analysis/sp_analysis/hbond_analysis/Makefile
                src/analysis/sp_analysis/contact_analysis/Makefile
                src/analysis/sp_analysis/sasa_analysis/Makefile
                src/analysis/mode_analysis/Makefile
                src/analysis/mode_analysis/avecrd_analysis/Makefile
                src/analysis/mode_analysis/eigmat_analysis/Makefile
                src/analysis/mode_analysis/flccrd_analysis/Makefile
                src/analysis/mode_analysis/pcavec_drawer/Makefile
                src/analysis/mode_analysis/prjcrd_analysis/Makefile
                src/analysis/converter/Makefile
                src/analysis/converter/crd_convert/Makefile
                src/analysis/converter/pcrd_convert/Makefile
                src/analysis/converter/remd_convert/Makefile
                src/analysis/converter/rst_convert/Makefile
                src/analysis/converter/rst_upgrade/Makefile
                src/analysis/converter/cg_convert/Makefile
                src/analysis/free_energy/Makefile
                src/analysis/free_energy/mbar_analysis/Makefile
                src/analysis/free_energy/meanforce_analysis/Makefile
                src/analysis/free_energy/pmf_analysis/Makefile
                src/analysis/free_energy/wham_analysis/Makefile
                src/analysis/interface/Makefile
                src/analysis/interface/dssp_interface/Makefile
                src/analysis/clustering/Makefile
                src/analysis/clustering/kmeans_clustering/Makefile
                src/analysis/utilities/Makefile
                src/analysis/utilities/pathcv_analysis/Makefile
                src/analysis/utilities/qmmm_generator/Makefile
                src/analysis/utilities/rpath_generator/Makefile
                src/analysis/utilities/emmap_generator/Makefile
                src/analysis/utilities/morph_generator/Makefile])])
AC_OUTPUT
