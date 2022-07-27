.. highlight:: bash
.. _getting_started:

==============================================================
Getting Started
==============================================================

Installation of GENESIS
=======================

Requirements
------------

Compilers
^^^^^^^^^

**GENESIS** works on various systems: laptop PCs, workstations, cluster machines, and supercomputers.
Since the source code of **GENESIS** is mainly written in
Fortran language, Fortran compiler is the first requirement for installation.
In addition, "preprocessor" is required, because the source code is
"processed" according to the user's computer environment before the compilation.
One of the commonly used Fortran compilers is ``gfortran``, which
is freely available as part of the GNU Compiler Collection (GCC).
In this case, ``cpp`` is selected as a preprocessor, which is also available freely.
Another recommended Fortran compiler is ``ifort`` provided
by Intel Corporation which enables us to run the program much faster on Intel CPU.
In the Intel compiler package, ``fpp`` is provided as a preprocessor.
Fujitsu compiler ``frtpx``, which also functions as a preprocessor,
is suitable for Fujitsu machines like FX100.

MPI and OpenMP
^^^^^^^^^^^^^^

Both **ATDYN** and **SPDYN** work on multiple CPU cores using MPI
(Message Passing Interface) and OpenMP protocols (hybrid MPI+OpenMP).
MPI and OpenMP are commonly used for parallel computing.
In general, MPI is employed for communication between different machines,
nodes, or processors, where the memory is not shared among them (distributed-memory).
On the other hand, OpenMP is employed in a single processor, and thus,
memory is shared in the parallel computation.

OpenMP is natively supported in most modern Fortran compilers.
As for MPI, however, the users may have to install MPI libraries by
themselves, especially, in the case of laptop PCs and workstations.
One of the commonly used MPI software is OpenMPI (https://www.open-mpi.org/).
When the users install the OpenMPI libraries in the computer, the users must
specify Fortran and C compilers (e.g., ``gfortran`` and ``gcc``) to be used with MPI.
After installing the libraries, the users can use ``mpif90``, ``mpicc``, and ``mpirun``,
which are necessary to compile and run the program that is parallelized with MPI. 
OpenMPI is available freely, and the example installation scheme is shown in :ref:`Appendix`.
Intel and Fujitsu Corporations are also providing their own MPI libraries for parallel computation.

Mathematical libraries
^^^^^^^^^^^^^^^^^^^^^^

**GENESIS** utilizes mathematical libraries such as LAPACK and BLAS (http://www.netlib.org/lapack/).
These libraries enable us to efficiently solve complicated mathematical
equations such as eigenvalue problems and singular value decomposition.
The users have to install these libraries by themselves,
if they are not installed in the computer (see :ref:`Appendix`).
In the case of the Intel and Fujitsu compilers, Intel MKL and Fujitsu SSL II are automatically selected, respectively.

GPGPU
^^^^^

**SPDYN** works not only with CPU but also with CPU+GPU.
Some of the source code in **SPDYN** are written in CUDA,
which enables us to effectively run the program on NVIDIA GPU cards.
If the users want to run **SPDYN** with GPGPU calculations, 
the CUDA toolkit (https://developer.nvidia.com/cuda-toolkit) must be also installed in the computer.
Note that OpenACC is not employed in **GENESIS** currently.

-------------------------------------------------------------------------

The recommended compilers, preprocessors, and libraries for **GENESIS** are listed below.
Please make sure that at least one of them in each section is installed on your system (GPU is optional).
If the users do not use the Intel or Fujitsu compilers,
the combination of GCC compiler, GCC preprocessor, and OpenMPI is recommended.

* Operating systems (see :ref:`Appendix`)

  * Linux
  * macOS
  * Windows 10/11

* Fortran and C compilers

  * GCC compiler ``gfortran``, ``gcc`` (version 4.4.7 or higher is required)
  * Intel compiler ``ifort``, ``icc``
  * Fujitsu compiler ``frtpx``, ``fccpx``

* Preprocessors

  * GCC preprocessor ``cpp``
  * Intel preprocessor ``fpp``
  * Fujitsu compiler ``frtpx`` 

* MPI libraries for parallel computing

  * OpenMPI ``mpirun``, ``mpif90``, ``mpicc``
  * Intel MPI
  * Fujitsu MPI

* Numerical libraries for mathematical algorithms

  * LAPACK/BLAS
  * Intel Math Kernel Library (MKL)
  * Fujitsu Scientific Subroutine Library (SSL II)

* GPU **(Optional)**

  * NVIDIA GPU cards which support Compute Capability (CC) 3.5 or higher
  * The following GPU cards and CUDA versions have been tested by the GENESIS developers

    * NVIDIA K20, K40, P100, TITAN V, GTX 1080, GTX 1080Ti, RTX 2080, RTX 2080Ti, RTX 3090
    * CUDA ver. 8.0, 9.0, 9.1, 9.2, 10.0, 11.4, 11.5

.. note::
   If you are using a supercomputer at a university or research institute,
   the above requirements are likely already provided, so you do not need to install them yourself.
   Please refer to the users' guide of the supercomputer, or consult the system administrator.

   In general, the latest version of CUDA does not support the latest version of GCC compiler.
   If you cannot compile GENESIS with new CUDA and new GCC compiler,
   please first make an attempt to install the new CUDA using an older GCC compiler,
   and then install GENESIS with those CUDA and GCC compilers.


.. raw:: latex

    \clearpage


General scheme for installation
-------------------------------

Step1. Download the source code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The source code of **GENESIS** is available in the GENESIS website (https://www.r-ccs.riken.jp/labs/cbrt/download/).
The users have to first uncompress the download file in an appropriate directory.
Here, we assume that the users install **GENESIS** in "$HOME/genesis".
The "src" directory contains the source code, and "COPYING" is the software license.
:: 
  $ mkdir $HOME/genesis
  $ cd $HOME/genesis
  $ mv ~/Downloads/genesis-2.0.0.tar.bz2 ./
  $ tar xvfj genesis-2.0.0.tar.bz2
  $ cd genesis-2.0.0
  $ ls
  AUTHORS         Makefile        bootstrap      depcomp
  COPYING         Makefile.am     cleanup        fortdep.py
  COPYING.LESSER  Makefile.in     compile        install-sh
  ChangeLog       NEWS            config.log     missing
  HOW_TO          README          config.status  src
  INSTALL         aclocal.m4      configure      version_scripts
  LICENSE         autom4te.cache  configure.ac


Step2. Configure
^^^^^^^^^^^^^^^^

In order to compile the source code, the users execute the "configure" script in the directory.
This script automatically detects appropriate compilers, preprocessors, and libraries
in the users' computer, and create "``Makefile``".
:: 

   $ ./configure

If you encountered a failure in the configure command, please check the error message carefully.
You may have to add appropriate options in this command according to your computer environment (see :ref:`Advanced installation`).
The followings are possible suggestions to solve frequent problems.
Other solutions might be found in the online page (https://www.r-ccs.riken.jp/labs/cbrt/installation/).

  * First of all, please check whether the Fortran and C compilers are installed in your computer.
    If you are going to run GENESIS with multiple CPUs, you should additionally install
    MPI libraries such as OpenMPI before compiling GENESIS (see :ref:`appendix`).

  * If you see the error message "configure: error: lapack is not correct",
    make sure that the BLAS or LAPACK libraries are installed in the computer (see also :ref:`Appendix`).
    The users may also have to set the path to the libraries in the "configure" command
    with the "``LAPACK_LIBS``" or "``LAPACK_PATH``" option (see :ref:`Advanced installation`).

  * If you see the error message "configure: error: Fortran compiler cannot create executables",
    it may imply that the path to the installed compilers or MPI libraries might not be correctly
    set in "~/.bashrc" or "~/.bash_profile" (see :ref:`appendix`).
    This configure script automatically detects "mpiifort", "mpif90", "mpifrtpx", or "mpifrt" for Fortran compiler,
    and "mpicc", "mpifccpx", or "mpifrt" for C compiler.
    The error message may indicate that the detection was failed due to some reasons.
    For example, if you installed OpenMPI in your computer, both "mpif90" and "mpicc" should be detected.
    Please check the path to these executables by typing the "which" command (e.g., which mpif90) in the terminal window.
    If you cannot see any paths, setting of the path in "~/.bashrc" or "~/.bash_profile" might have a mistake (see :ref:`appendix`).
    You should also check typing mistakes of the path.

  * If the recommended software are not used in compilation, warning messages might be
    displayed in the terminal when the configure command is executed.
    Those messages are just a warning (not an error), and you may continue the compilation.
    However, we strongly recommended you to verity the installation in such cases
    (see :ref:`Verify installation`).

  * In some supercomputer systems, "module load [module]" command is required to use
    compliers, and need to be set before the configure. See the user guide of the system.

  * Try "autoreconf" or "./bootstrap" before the configure command,
    if your computer environment is significantly different from what we assume
    and/or if you modify "configure.ac" or "Makefile.am" by yourself.


Step3. Make install
^^^^^^^^^^^^^^^^^^^

After the "configure" command is successful, type the following command to compile and install **GENESIS**.
All programs in **GENESIS** are compiled and installed into the "./bin" directory by default.
:: 

   $ make install

If you encountered a failure, please check the error message carefully.
In many cases, errors are caused by invalid path of compilers and libraries.
The followings are possible suggestions to solve frequent problems. 
Other solutions might be found in the online page
(https://www.r-ccs.riken.jp/labs/cbrt/installation/).

  * If the error message is like "/usr/bin/ld: cannot find -lblas" or "/usr/bin/ld: cannot find -llapack",
    make sure that the BLAS or LAPACK libraries are installed in the computer (see also :ref:`Appendix`).
    The users may also have to set the path to the libraries in the "configure" command
    with the "``LAPACK_LIBS``" or "``LAPACK_PATH``" option (see :ref:`Advanced installation`).

  * If you have installed additional software or libraries to solve a make error,
    please execute "make clean", and try Step2 and "make install" again.


Step4. Confirmation
^^^^^^^^^^^^^^^^^^^

After the installation is successfully finished, the following binary files are found in the "bin" directory.
There are 44 programs in total. Brief description of each program is shown in :ref:`available_programs`.
:: 

  $ ls ./bin
  atdyn               hb_analysis          qval_residcg_analysis
  avecrd_analysis     hbond_analysis       rdf_analysis
  cg_convert          kmeans_clustering    remd_convert
  comcrd_analysis     lipidthick_analysis  rg_analysis
  contact_analysis    mbar_analysis        ring_analysis
  crd_convert         meanforce_analysis   rmsd_analysis
  density_analysis    morph_generator      rpath_generator
  diffusion_analysis  msd_analysis         rst_convert
  distmat_analysis    pathcv_analysis      rst_upgrade
  drms_analysis       pcavec_drawer        sasa_analysis
  dssp_interface      pcrd_convert         spdyn
  eigmat_analysis     pmf_analysis         tilt_analysis
  emmap_generator     prjcrd_analysis      trj_analysis
  flccrd_analysis     qmmm_generator       wham_analysis
  fret_analysis       qval_analysis


.. raw:: latex

    \clearpage

.. _Advanced installation:

Advanced installation
---------------------

In the above scheme, **GENESIS** is installed with default options,
and all installed programs run on CPU with double precision calculation.
The users can specify additional options in the configure command
according to the users' computer environment or desired conditions.
The full lists of the available options are obtained by "``./configure --help``".
The representative options are as follows.

``--enable-single``

  Turn on single precision calculation. In this case, only **SPDYN** is installed.

``--enable-mixed``

  Turn on mixed precision calculation. In this case, only **SPDYN** is installed.

``--enable-double`` (default)

  Turn on double precision calculation. In this case, all binaries are installed.

``--enable-gpu``

  Turn on GPGPU calculation. In this case, only **SPDYN** is installed.

``--with-cuda=PATH``

  Define path to the CUDA libraries manually.

``--disable-mpi``

  Turn off MPI parallelization. In this case, **SPDYN** is not installed.

``--disable-parallel_IO`` (default)

  Do not install the parallel I/O tool (**prst_setup**)

``--enable-debug``

  Turn on program debugging (see below)

``--prefix=PREFIX``

  Install the programs in the directory designated by ``PREFIX``

``--with-msmpi``

  Turn on use of MSMPI. Compilation and execution must be done on Windows10/11.

-------------------------------------------------------------------------

Configuration with specified compilers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The users can explicitly specify the compiler in the configure command.
Fortran compiler is specified with ``FC`` and ``F77``, and C compiler with ``CC``.
For example, in the case of "mpiifort" and "mpiicc", the following options are added:
:: 
   $ ./configure CC=mpiicc FC=mpiifort F77=mpiifort


Configuration with an explicit path to LAPACK/BLAS libraries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The following is an example command to set the path to LAPACK and BLAS libraries
that are installed in $HOME/Software/lapack-3.10.1/ (see also :ref:`Appendix`).
Please be careful about the filename of the installed libraries.
If the BLAS libraries are installed as "librefblas.a", the option "-lrefblas" must be used.
If "librefblas.a" is renamed to "libblas.a", the following command can be used.
Linking with the reverse order of "-llapack" and "-lblas" might
also cause a failure of installation of **GENESIS**. :: 
   $ ./configure LAPACK_LIBS="-L$HOME/Software/lapack-3.10.1 -llapack -lblas"

or use the "``LAPACK_PATH``" option:
:: 
   $ ./configure LAPACK_PATH=$HOME/Software/lapack-3.10.1


Configuration for single-precision calculation on CPU
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The following command is used to turn on single-precision calculation in **SPDYN**.
In this case, force calculations are carried out with single precision,
while integration of the equations of motion as well as
accumulation of the force and energy are still done with double-precision.
:: 
   $ ./configure --enable-single

The mixed precision calculation is also prepared.
:: 
   $ ./configure --enable-mixed

Only **SPDYN** that works on CPU will be installed with these options.
If the user additionally needs analysis tools as well as **ATDYN**,
one must prepare another GENESIS directory, and install without "``--enable-single``" or "``--enable-mixed``" option.


Configuration for GPGPU calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In the following command, the users install **SPDYN** that works on CPU+GPU with single-precision calculation.
If "``--enable-single``" or "``--enable-mixed``" is omitted in the command, **SPDYN** works on CPU+GPU with double-precision calculation.
:: 
   $ ./configure --enable-single --enable-gpu

Here, if the users encountered an error message like "nvcc: command not found",
make sure that the CUDA Toolkit is installed in the computer.
In typical Linux workstations or cluster machines, 
CUDA is installed in "/usr/local/cuda-x.x/" or "/usr/lib/x86_64-linux-gnu/",
and "nvcc" should be in a "bin" directory of the install directory.
The path to "nvcc" and CUDA libraries should be set in a startup file
such as "~/.bashrc".  For example, add the following information to "~/.bashrc" 
in the case of CUDA 11.4, 
:: 
  CUDAROOT=/usr/local/cuda-11.4
  export PATH=$CUDAROOT/bin:$PATH
  export LD_LIBRARY_PATH=$CUDAROOT/lib64:/lib:$LD_LIBRARY_PATH

then reload "~/.bashrc" and try the configure command again.
If there still remain some troubles, explicitly specify a path to CUDA 
libraries in the configure command by:
:: 
   $ ./configure --enable-single --enable-gpu --with-cuda=/usr/local/cuda-11.4


Configuration for supercomputer systems
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The configuration for supercomputer systems may require non-standard setups. 
In the online usage page, we describe recommended configure options for some supercomputers
(https://www.r-ccs.riken.jp/labs/cbrt/usage/).

For example, the following commands are used to compile **GENESIS** on Fugaku in RIKEN.
Note that the parallel I/O tool (**prst_setup**) is not compiled in this configuration,
because Fujitsu compiler has a trouble in compiling **prst_setup** (see aldo :ref:`available_programs`).
:: 
   $ ./configure --enable-mixed --host=Fugaku


Configuration for single CPU calculations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
By specifying the "``--disable-mpi``" option, the users can install **GENESIS** that can work on one CPU.
The configure script automatically looks for "gfortran", "ifort", "frt", or "frtpx" for Fortran compiler, 
and "gcc", "icc", "fcc", or "fccpx" for C compiler. Therefore, in this case MPI libraries are not required
for the installation and execution of **GENESIS**. **ATDYN** and analysis tools are installed.
:: 
   $ ./configure --disable-mpi


Configuration for program debugging
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If the users encountered memory leak errors during the simulation using **GENESIS**,
the origin of the error might be tracked by using a program compiled with a debug option.
Note that the debug option makes the calculation much slow.
In this case, the runtime check is activated only for CPU codes,
even if the "--enable-gpu" option is added to the command.
:: 
   $ ./configure --enable-debug=3

Note that ``--enable-debug`` corresponds to ``--enable-debug=1``.

  * 0 = no debugging (default)

  * 1 = debugging without intensive optimization

  * 2 = LEVEL1 + debug information (``-g`` and ``-DDEBUG``)

  * 3 = LEVEL2 + memory check (if possible)

  * 4 = LEVEL3 + full check (intel compiler only)

----------------------------------------------------------------------------

Installation using multiple CPU cores (parallel compile)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In Step3, ``-j`` option is available, which enables quick compilation 
of the program using multiple CPU cores. The following command uses 4 CPU cores.
:: 
  $ make -j 4 install

If you encountered an error message like "Fatal Error: Can't delete temporary module
file '...': No such file or directory", please try "make install"
without the "``-j``" option.


.. raw:: latex

    \clearpage


.. _Verify installation:

Verify the installation
-----------------------

The users can verify the installation of **GENESIS** by using test sets 
which are available in the **GENESIS** website
(https://www.r-ccs.riken.jp/labs/cbrt/download/).
Please uncompress the downloaded file in an appropriate directory,
and move to the "regression_test" directory.
Note that the file name of the tar.bz2 file contains the date (year, month, and day),
so please change the following execution commands accordingly.
:: 
  $ cd $HOME/genesis
  $ mv ~/Downloads/tests-2.0.0_YYMMDD.tar.bz2 ./
  $ tar xvfj tests-2.0.0_YYMMDD.tar.bz2
  $ cd tests-2.0.0_YYMMDD/regression_test
  $ ls
  build       test_analysis    test_gamd_spdyn    test_rpath_atdyn
  charmm.py   test_atdyn       test_nonstrict.py  test_rpath_spdyn
  cleanup.sh  test_common      test_parallel_IO   test_spana
  fep.py      test_fep         test_remd.py       test_spdyn
  genesis.py  test_fep.py      test_remd_common   test_vib
  param       test_gamd.py     test_remd_spdyn    test_vib.py
  test.py     test_gamd_atdyn  test_rpath.py

In the sub-directories in "regression_test", the users can find a lot of input files ("``inp``"),
in which various combinations of simulation parameters are specified.
In addition, each sub-directory contains output file ("``ref``") obtained by the developers.
The users run "test.py", "test_remd.py", "test_rpath.py", and so on,
which enable automatic comparison between the users' and developers'
results for each MD algorithm.


Run the basic tests
^^^^^^^^^^^^^^^^^^^
The following is an example command to verify the two simulators
**atdyn** and **spdyn** for basic MD and energy minimization.
Here, the programs are executed using 1 CPU core with the "mpirun" command.
The users can increase the number of MPI processors according to
the users' computer environment, but only 1, 2, 4, or 8 are allowed in these tests.
Other MPI launchers such as "mpiexec" are also available in the command.
There are about 50 test sets, and each test should finish in a few seconds.
:: 
  $ export OMP_NUM_THREADS=1
  $ ./test.py "mpirun -np 1 ~/genesis/genesis-2.0.0/bin/atdyn"
  $ ./test.py "mpirun -np 1 ~/genesis/genesis-2.0.0/bin/spdyn"

If any tests cannot run, please check the following points:

* Number of OpenMP threads should be specified before running the tests (one is recommended).

* Original executable file name (e.g., **spdyn** and **atdyn**) must not be changed.

* Python ver. 3 is used.

.. note::
For **spdyn** on Fugaku in RIKEN, the number of MPI processors be greater than or equal to 8. We recommend to use mixed precision on Fugaku.

The "test.py" script compares energy in log file between the developer’s and user’s ones.
If the energy differences are less than the tolerance (default = 1.00e-08),
"Passed" is displayed for each test. Among the physical quantities in the log file,
virial is the most sensitive to numerical factors, and thus, the tolerance
for virial is set to a larger value (1.00e-06).
After all tests are finished, the total number of succeeded, 
failed, and aborted runs will be displayed at the end.
:: 
  Passed  46 / 46
  Failed  0 / 46
  Aborted 0 / 46

If all tests were passed, it means that your **GENESIS** can generate
identical results to the developer's **GENESIS**.
Note that the developer’s **GENESIS** was compiled with Intel compilers,
Intel MKL, OpenMPI library, and the double precision option on Intel CPUs.
If your computer system is significantly different from the developer’s one,
unexpected numerical errors may happen, which can cause failures in some tests.
If there were any aborted tests, the users had better check their log or error files carefully,
which exist in the tested sub-directory, and figure out why the error happened.
The followings are suggestions to solve typical problems:

* If some tests were aborted due to "memory allocation error", 
  the reason might come from limitation of the memory size.
  Namely, those tested systems were too large for your computer.
  The problem should not be so serious.

* Available number of MPI slots in your computer might be actually
  smaller than the given number of MPI processors.
  Try to use less number of MPI processors.

* Try to specify the "absolute path" to the program instead of using "relative path".

* Make sure that the MPI environment is properly set.

* Detailed solutions in specific supercomputer systems might be found in the GENESIS website
  (https://www.r-ccs.riken.jp/labs/cbrt/usage/).

Run the additional tests
^^^^^^^^^^^^^^^^^^^^^^^^
By using a similar way, the users can check other functions in **atdyn** and **spdyn**,
such as GaMD, REMD, path sampling, vibrational analysis, parallel I/O, and GPGPU calculation.
Available number of MPI processors depends on each test (test_gamd: 1, 2, 4, 8;
test_remd: 4, 8, 16, 32; test_rpath: 8; test_vib: 8; parallel_io: 8; gpu: 1, 2, 4, 8).
As for the GPGPU tests, the users must use **spdyn** that was installed with the "--enable-gpu" option.
The parallel_io tests require both **spdyn** and **prst_setup**.
Note that **prst_setup** is not installed in some cases according to the configure options or compilers (see :ref:`Advanced installation`).
In order to run the analysis tool tests, the users first move to "test_analysis",
and then execute "./test_analysis.py". Note that MPI is not used in the analysis tool tests.
In a similar way, the users can test the SPANA (spatial decomposition analysis) tool sets. 
SPANA tool sets are tested with 1, 2, 4, and 8 MPI processes.
:: 
  $ export OMP_NUM_THREADS=1
  $ ./test_gamd.py "mpirun -np 1 ~/genesis/genesis-2.0.0/bin/atdyn"
  $ ./test_gamd.py "mpirun -np 1 ~/genesis/genesis-2.0.0/bin/spdyn"
  $ ./test_remd.py "mpirun -np 4 ~/genesis/genesis-2.0.0/bin/atdyn"
  $ ./test_remd.py "mpirun -np 4 ~/genesis/genesis-2.0.0/bin/spdyn"
  $ ./test_fep.py "mpirun -np 8 ~/genesis/genesis-2.0.0/bin/spdyn"
  $ ./test_rpath.py "mpirun -np 8 ~/genesis/genesis-2.0.0/bin/atdyn"
  $ ./test_rpath.py "mpirun -np 8 ~/genesis/genesis-2.0.0/bin/spdyn"
  $ ./test_vib.py   "mpirun -np 8 ~/genesis/genesis-2.0.0/bin/atdyn"
  $ ./test.py "mpirun -np 8 ~/genesis/genesis-2.0.0/bin/spdyn" parallel_io
  $ ./test.py "mpirun -np 8 ~/genesis/genesis-2.0.0/bin/spdyn" gpu
  
  $ cd test_analysis
  $ ./cleanup.sh
  $ export OMP_NUM_THREADS=1
  $ ./test_analysis.py ~/genesis/genesis-2.0.0/bin/

  $ cd test_spana
  $ ./cleanup.sh
  $ export OMP_NUM_THREADS=1
  $ ./test_spana.py ~/genesis/genesis-2.0.0/bin/

.. note::
   Some tests might be using "abnormal" parameters or conditions in the input files
   for the sake of simple tests. Do not use such parameters in your research.
   "Normal" parameters are mainly introduced in this user manual or online tutorials.


.. raw:: latex

    \clearpage


Clean install and re-compilation
--------------------------------

The following commands are used to fully recompile **GENESIS**.
Note that the direct "make clean" command may not work
in the case where ``Makefiles`` were created in another machine.
In this case, the users must run the "./configure" command before "make clean".
:: 
  $ make clean
  $ ./configure [option]
  $ make install


Uninstall
---------

The user can uninstall **GENESIS** by just removing the program directory.
If the user changed the install directory by specifying "``--prefix=PREFIX``" in the configure command,
please remove the programs (**atdyn**, **spdyn**, and so on) in the "``PREFIX``" directory.
:: 
  $ rm -rf $HOME/genesis/genesis-2.0.0


.. raw:: latex

    \clearpage


Basic usage of GENESIS
======================

Running GENESIS on a command line
---------------------------------

The **GENESIS** programs are executed on a command line. 
The first argument is basically interpreted as an input file of the program.
The input file, which we call *control file* hereafter, contains parameters for simulations.
The following examples show typical usage of the **GENESIS** programs.
In the case of serial execution,
:: 
  $ [program_name] [control_file]

In the case of parallel execution with "mpirun",
:: 
  $ mpirun -np n [program_name] [control_file]

For example, **SPDYN** is executed in the following way using 8 MPI processors:
:: 
  $ mpirun -np 8 ~/genesis/genesis-2.0.0/bin/spdyn INP

The users should specify an OpenMP thread number explicitly before running the program.
Appropriate number of CPU cores must be used according to the user's computer environment (see also :ref:`available_programs`).
For example, if the users want to use 32 CPU cores in the calculation,
the following command might be executed.
:: 
  $ export OMP_NUM_THREADS=4
  $ mpirun -np 8 ~/genesis/genesis-2.0.0/bin/spdyn INP

As for the analysis tools, the usage is almost same, but mpirun is not used.
Note that some analysis tools (e.g., mbar_analysis, wham_analysis, msd_analysis, and drms_analysis) are 
parallelized with OpenMP.
:: 
  # RMSD analysis tool 
  $ ~/genesis/genesis-2.0.0/bin/rmsd_analysis INP

  # MBAR analysis
  $ export OMP_NUM_THREADS=4
  $ ~/genesis/genesis-2.0.0/bin/mbar_analysis INP


Automatic generation of a template control file
-----------------------------------------------

Basic usage of each program is shown by executing the program with the ``-h`` option. In addition, sample control file of each program can be obtained with the ``-h ctrl`` option:
:: 
  # Show the usage of the program
  $ [program_name] -h

  # Display a template control file
  $ [program_name] -h ctrl [module_name]

For example, in the case of **SPDYN**, the following messages are displayed:
:: 

  $ spdyn -h 

  # normal usage
    % mpirun -np XX ./spdyn INP

  # check control parameters of md
    % ./spdyn -h ctrl md

  # check control parameters of min
    % ./spdyn -h ctrl min

  # check control parameters of remd
    % ./spdyn -h ctrl remd

  # check control parameters of rpath
    % ./spdyn -h ctrl rpath

  # check all control parameters of md
    % ./spdyn -h ctrl_all md

  (skipped...)

This message tells the users that **SPDYN** can be executed with mpirun.
A template control file for molecular dynamics simulation (md) can be
generated by executing **SPDYN** with the ``-h ctrl md`` option.
The same way is applicable for energy minimization (min),
replica exchange simulation (remd), and replica path sampling simulation (rpath).
The template control file for energy minimization is shown below.
If the users want to show all available options, please specify ``ctrl_all`` instead of ``ctrl``.
The users can edit this template control file to perform the simulation that the users want to do. 
:: 

  $ ~/genesis/genesis-2.0.0/bin/spdyn -h ctrl min > INP

  $ less INP

  [INPUT]
  topfile = sample.top      # topology file
  parfile = sample.par      # parameter file
  psffile = sample.psf      # protein structure file
  pdbfile = sample.pdb      # PDB file

  [ENERGY]
  forcefield    = CHARMM    # [CHARMM,AMBER,GROAMBER,GROMARTINI]
  electrostatic = PME       # [CUTOFF,PME]
  switchdist    = 10.0      # switch distance
  cutoffdist    = 12.0      # cutoff distance
  pairlistdist  = 13.5      # pair-list distance

  [MINIMIZE]
  method        = SD        # [SD]
  nsteps        = 100       # number of minimization steps

  [BOUNDARY]
  type          = PBC       # [PBC, NOBC]

.. raw:: latex

    \clearpage


Control file
============

In the control file, detailed simulation conditions are specified.
The control file consists of several sections (e.g., **[INPUT]**,
**[OUTPUT]**, **[ENERGY]**, **[ENSEMBLE]**, and so on), 
each of which contains closely-related keywords.
For example, in the **[ENERGY]** section, parameters are specified for the potential energy calculation
such as a force field type and cut-off distance.
In the **[ENSEMBLE]** section, there are parameters to select the algorithm to control the temperature 
and pressure in addition to the target temperature and pressure of the system.
Here, we show example control files for the energy minimization and normal molecular dynamics simulations.

Example control file for the energy minimization
------------------------------------------------

The control file for the energy minimization must include a **[MINIMIZE]** section (see :ref:`Minimize`).
By using the following control file, the users carry out 2,000-step
energy minimization with the steepest descent algorithm (SD).
The CHARMM36m force field is used, and the particle mesh Ewald (PME) method
is employed for the calculation of long-range interaction.
:: 

  [INPUT]
  topfile = top_all36_prot.rtf     # topology file
  parfile = par_all36m_prot.prm    # parameter file
  strfile = toppar_water_ions.str  # stream file
  psffile = build.psf              # protein structure file
  pdbfile = build.pdb              # PDB file

  [OUTPUT]
  dcdfile = min.dcd                # coordinates trajectory file
  rstfile = min.rst                # restart file

  [ENERGY]
  forcefield        = CHARMM       # CHARMM force field
  electrostatic     = PME          # Particl mesh Ewald method
  switchdist        = 10.0         # switch distance (Ang)
  cutoffdist        = 12.0         # cutoff distance (Ang)
  pairlistdist      = 13.5         # pair-list cutoff distance (Ang)
  vdw_force_switch  = YES          # turn on van der Waals force switch
  contact_check     = YES          # turn on clash checker
  
  [MINIMIZE]
  method            = SD           # Steepest descent method
  nsteps            = 2000         # number of steps
  eneout_period     =  100         # energy output freq
  crdout_period     =  100         # coordinates output frequency
  rstout_period     = 2000         # restart output frequency
  nbpdate_period    =   10         # pairlist update frequency

  [BOUNDARY]
  type              = PBC          # periodic boundary condition
  box_size_x        = 64.0         # Box size in X dimension (Ang)
  box_size_y        = 64.0         # Box size in Y dimension (Ang)
  box_size_z        = 64.0         # Box size in Z dimension (Ang)

.. raw:: latex

    \clearpage

Example control file for normal MD simulations
----------------------------------------------

The control file for normal MD simulations must include a **[DYNAMICS]** section (see :ref:`Dynamics`).
By using the following control file, the users carry out a
100-ps MD simulation at :math:`T` = 298.15 K and :math:`P` = 1 atm in the NPT ensemble.
The equations of motion are integrated by the RESPA algorithm
with a time step of 2.5 fs, and the bonds of light atoms (hydrogen atoms)  
are constrained using the SHAKE/RATTLE and SETTLE algorithms.
The temperature and pressure are controlled with the Bussi thermostat and barostat.
:: 

  [INPUT]
  topfile = top_all36_prot.rtf     # topology file
  parfile = par_all36m_prot.prm    # parameter file
  strfile = toppar_water_ions.str  # stream file
  psffile = build.psf              # protein structure file
  pdbfile = build.pdb              # PDB file
  rstfile = min.rst                # restart file

  [OUTPUT]
  dcdfile = md.dcd                 # coordinates trajectory file
  rstfile = md.rst                 # restart file

  [ENERGY]
  forcefield        = CHARMM       # CHARMM force field
  electrostatic     = PME          # Particl mesh Ewald method
  switchdist        = 10.0         # switch distance (Ang)
  cutoffdist        = 12.0         # cutoff distance (Ang)
  pairlistdist      = 13.5         # pair-list cutoff distance (Ang)
  vdw_force_switch  = YES          # turn on van der Waals force switch

  [DYNAMICS]
  integrator        =   VRES       # RESPA integrator
  timestep          = 0.0025       # timestep (2.5fs)
  nsteps            =  40000       # number of MD steps (100ps)
  eneout_period     =    400       # energy output period (1ps)
  crdout_period     =    400       # coordinates output period (1ps)
  rstout_period     =  40000       # restart output period
  nbupdate_period   =     10       # nonbond update period
  elec_long_period  =      2       # period of reciprocal space calculation
  thermostat_period =     10       # period of thermostat update
  barostat_period   =     10       # period of barostat update

  [CONSTRAINTS]
  rigid_bond        = YES          # constraint all bonds involving hydrogen

  [ENSEMBLE]
  ensemble          = NPT          # NPT ensemble
  tpcontrol         = BUSSI        # BUSSI thermostat and barostat
  temperature       = 300          # target temperature (K)
  pressure          = 1.0          # target pressure (atm)
  group_tp          = YES          # usage of group tempeature and pressure

  [BOUNDARY]
  type              = PBC          # periodic boundary condition

