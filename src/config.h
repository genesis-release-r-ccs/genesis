/* src/config.h.  Generated from config.h.in by configure.  */
/* src/config.h.in.  Generated from configure.ac by autoheader.  */

/* c compiler version */
#define COMPILE_CC_VER "icc (ICC) 18.0.2 20180210"

/* c flags */
#define COMPILE_CFLAGS "-O3 -ip -axCORE-AVX2  -qopenmp"

/* defined variables */
#define COMPILE_DEFINED_VARIABLES " -DMPI -DOMP -DFFTE -DLAPACK -DDSFMT_MEXP=19937 -DINTEL"

/* fortran flags */
#define COMPILE_FCFLAGS "-xHost -axCORE-AVX512 -g -O3 -ip -mkl=parallel  -assume byterecl -qopenmp "

/* fortran compiler version */
#define COMPILE_FC_VER "ifort (IFORT) 18.0.2 20180210"

/* genesis version */
#define COMPILE_GENESIS_VERSION "$GENESIS_VERSION$"

/* hostname */
#define COMPILE_HOST "dolphin1"

/* ld flags */
#define COMPILE_LDFLAGS " -assume byterecl -qopenmp  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_lapack95_lp64 "

/* cuda version */
/* #undef COMPILE_NVCC_VER */

/* username */
#define COMPILE_USER "jung"

/* defined if cuda_gpu is used. */
/* #undef CUDAGPU */

/* defined if Debug is used. */
/* #undef DEBUG */

/* defined always. */
#define DSFMT_MEXP 19937

/* defined if FFTE is used. */
#define FFTE 1

/* defined if fapp is used. */
/* #undef FJ_PROF_FAPP */

/* defined if fipp is used. */
/* #undef FJ_PROF_FIPP */

/* defined if fipp is used. */
/* #undef FJ_TIMER_2 */

/* defined if fipp is used. */
/* #undef FJ_TIMER_DETAIL */

/* defined if --host=Fugaku is set. */
/* #undef FUGAKU */

/* defined if HM_DISK is used. */
/* #undef HM_DISK */

/* defined if Intel compiler is used. */
#define INTEL 1

/* defined if K-computer compiler is used. */
/* #undef KCOMP */

/* defined if LAPACK is used. */
#define LAPACK 1

/* defined if MPI is used. */
#define MPI 1

/* defined if OpenMP is used. */
#define OMP 1

/* Name of package */
#define PACKAGE "genesis"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "genesis@riken.jp"

/* Define to the full name of this package. */
#define PACKAGE_NAME "genesis"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "genesis 2.0.0"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "genesis"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "2.0.0"

/* defined if pgi and cuda are used. */
/* #undef PGICUDA */

/* defined if PKTIMER is used. */
/* #undef PKTIMER */

/* defined if gpu is used. */
/* #undef USE_GPU */

/* Version number of package */
#define VERSION "2.0.0"

/* defined if _LARGE_INT is used. */
/* #undef _LARGE_INT */

/* defined if _MIXED is used. */
/* #undef _MIXED */

/* defined if _SINGLE is used. */
#define _SINGLE 1

/* defined if GCC gfortran compiler is used. */
/* #undef __GFORTRAN__ */

/* defined if pgi compiler is used. */
/* #undef __PGI */
