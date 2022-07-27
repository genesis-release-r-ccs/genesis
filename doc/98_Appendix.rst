.. highlight:: bash
.. _appendix:

=======================================================================
Appendix
=======================================================================

Install the requirements in Linux
=================================

In the first sub-section, we explain how to install the requirements
using the package manager "apt" in Ubuntu/Debian.
If you want to install them from the source codes,
or if you want to use other Linux systems like CentOS,
please see the second sub-section.


For Ubuntu/Debian users
-----------------------

GNU compilers and build tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, we install compilers and build tools.
:: 
  $ sudo apt update
  $ sudo apt install build-essential
  $ sudo apt install gfortran autoconf automake


OpenMPI
^^^^^^^

Then, we install OpenMPI. Note that the development version (XXX-dev) should be installed. 
:: 
  $ sudo apt install openmpi-bin libopenmpi-dev

  $ which mpirun mpif90 mpicc
  /usr/bin/mpirun
  /usr/bin/mpif90
  /usr/bin/mpicc


LAPACK/BLAS libraries
^^^^^^^^^^^^^^^^^^^^^

Finally, we install LAPACK/BLAS libraries. Again, development version (XXX-dev) is installed.
:: 
  $ sudo apt install liblapack-dev

  $ ls /usr/lib/x86_64-linux-gnu/liblapack.*
  $ ls /usr/lib/x86_64-linux-gnu/libblas.*


.. raw:: latex

    \clearpage


For CentOS/RedHat users
-----------------------

Here, we explain how to install OpenMPI and LAPACK/BLAS libraries from the source codes.
We assume that the users already installed GNU compilers.
The following schemes are commonly applicable to typical Linux systems including CentOS and Red Hat.

OpenMPI
^^^^^^^

The source code of OpenMPI is availabe in https://www.open-mpi.org/.
The following commands install OpenMPI 3.1.5 in the user's local directory
"$HOME/Software/mpi" as an example.
Here, we use GNU compilers (gcc, g++, and gfortran).
:: 
  $ cd $HOME
  $ mkdir Software
  $ cd Software

  $ mkdir build
  $ cd build

  $ wget https://download.open-mpi.org/release/open-mpi/v3.1/openmpi-3.1.5.tar.gz

  $ tar -xvf openmpi-3.1.5.tar.gz
  $ cd openmpi-3.1.5

  $ ./configure --prefix=$HOME/Software/mpi CC=gcc CXX=g++ F77=gfortran FC=gfortran

  $ make all
  $ make install

The following information is added in "~/.bash_profile" (or "~/.bashrc").
:: 
  MPIROOT=$HOME/Software/mpi
  export PATH=$MPIROOT/bin:$PATH
  export LD_LIBRARY_PATH=$MPIROOT/lib:$LD_LIBRARY_PATH
  export MANPATH=$MPIROOT/share/man:$MANPATH

Launch another terminal window or reload "~/.bash_profile" (or "~/.bashrc"):
:: 
  $ source ~/.bash_profile

The OpenMPI tools should be installed in "$HOME/Software/mpi/bin".
:: 
  $ which mpirun mpif90 mpicc
  ~/Software/mpi/bin/mpirun
  ~/Software/mpi/bin/mpif90
  ~/Software/mpi/bin/mpicc

If you want to uninstall OpenMPI, just remove the directory "mpi" in "Software".


LAPACK/BLAS libraries
^^^^^^^^^^^^^^^^^^^^^

The source code of LAPACK/BLAS is availabe in http://www.netlib.org/lapack/.
The following commands install LAPACK 3.10.1
in the user's local directory "$HOME/Software/lapack-3.10.1" as an example.
The BLAS library is also installed. We use GNU compilers (gcc and gfortran).
:: 
  $ cd $HOME/Software
  $ wget https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.10.1.tar.gz
  $ tar -xvf lapack-3.10.1.tar.gz
  $ cd lapack-3.10.1

  $ cp make.inc.example make.inc
  $ make blaslib
  $ make lapacklib

  $ ls lib*
  liblapack.a  librefblas.a

  $ ln -s librefblas.a ./libblas.a

The following information is added in “~/.bash_profile” (or “~/.bashrc”).
:: 
  export LAPACK_PATH=$HOME/Software/lapack-3.10.1

Launch another terminal window or reload "~/.bash_profile" (or "~/.bashrc"):
:: 
  $ source ~/.bash_profile

If you want to uninstall LAPACK/BLAS, just remove the directory "lapack-3.10.1" in "Software".

.. raw:: latex

    \clearpage



Install the requirements in Mac
===============================

We recommend the Mac users to utilize "Xcode" for the installation of GENESIS,
and also to install "OpenMPI" from the source code to avoid a "clang" problem (see below).

Install general tools
---------------------

Xcode and Homebrew
^^^^^^^^^^^^^^^^^^
 
"Xcode" is available in the Mac App Store (https://developer.apple.com/xcode/),
and it is free of charge.
After the installation of Xcode, all tasks described below will be done on "Terminal".
The "Terminal app" is in the "Utilities" folder in Applications.
Please launch the Terminal. This terminal is almost same with that in Linux.

We recommend you to further install "Homebrew", which enables easy installation
of various tools such as compilers. If you have already installed "MacPorts",
you do not need to install "Homebrew" to avoid a conflict between "Homebrew" and "MacPorts".
In the Homebrew website (https://brew.sh/), you can find a long command
like "``/usr/bin/ruby -e "$(curl -fsSL https://...``".
To install homebrew, execute that command in the Terminal prompt.


GNU compilers and build tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, we install "gcc", "autoconf", "automake", and other tools via homebrew:
:: 
  $ brew install gcc
  $ brew install autoconf
  $ brew install automake
  $ brew install wget

To confirm the installation of "gcc", let us type the following commands:
:: 
  $ which gcc
  /usr/bin/gcc

  $ gcc --version
  ...
  Apple LLVM version 10.0.1 (clang-1001.0.46.4)

These messages tell us that "gcc" is installed in the "/usr/bin" directory.
However, this gcc is not a "real" GNU compiler, and it is linked to another compiler "clang".
If you use this gcc for the installation of OpenMPI, 
it can cause a trouble in compiling **GENESIS** with a certain option.
Therefore, you have to use a "real" GNU compiler, which is actually installed in "/usr/local/bin".
For example, if you have installed gcc ver. 9, you can find it as "gcc-9" in "/usr/local/bin".
:: 
  $ ls /usr/local/bin/gcc*
  /usr/local/bin/gcc-9      /usr/local/bin/gcc-ar-9      ...

  $ gcc-9 --version
  gcc-9 (Homebrew GCC 9.2.0) 9.2.0


Install libraries
-----------------

OpenMPI
^^^^^^^

We then install "OpenMPI". We specify "real" GNU compilers explicitly in the configure command.
The following commands install OpenMPI in the user's local directory "$HOME/Software/mpi".
:: 
  $ cd $HOME
  $ mkdir Software
  $ cd Software
  $ mkdir build
  $ cd build

  $ wget https://download.open-mpi.org/release/open-mpi/v3.1/openmpi-3.1.5.tar.gz

  $ tar -xvf openmpi-3.1.5.tar.gz
  $ cd openmpi-3.1.5

  $ ./configure --prefix=$HOME/Software/mpi CC=gcc-9 CXX=g++-9 F77=gfortran-9 FC=gfortran-9

  $ make all
  $ make install

The following information is added in "~/.bash_profile" (or "~/.bashrc").
:: 
  MPIROOT=$HOME/Software/mpi
  export PATH=$MPIROOT/bin:$PATH
  export LD_LIBRARY_PATH=$MPIROOT/lib:$LD_LIBRARY_PATH
  export MANPATH=$MPIROOT/share/man:$MANPATH

Launch another terminal window or reload "~/.bash_profile" (or "~/.bashrc"):
:: 
  $ source ~/.bash_profile

The OpenMPI tools should be installed in "$HOME/Software/mpi/bin".
:: 
  $ which mpirun mpif90 mpicc
  /Users/[username]/Software/mpi/bin/mpirun
  /Users/[username]/Software/mpi/bin/mpif90
  /Users/[username]/Software/mpi/bin/mpicc

Make sure that "mpicc" and "mpif90" are linked to "gcc-9" and "gfortran-9", respectively.
:: 
  $ mpicc --version
  gcc-9 (Homebrew GCC 9.2.0) 9.2.0

  $ mpif90 --version
  GNU Fortran (Homebrew GCC 9.2.0) 9.2.0

If you want to uninstall OpenMPI, just remove the directory "mpi" in "Software".


LAPACK/BLAS libraries
^^^^^^^^^^^^^^^^^^^^^

Finally, we install LAPACK and BLAS libraries.
Again, "real" GNU compilers are used for the install.
:: 
  $ cd $HOME/Software
  $ wget https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.10.1.tar.gz
  $ tar -xvf v3.10.1.tar.gz
  $ cd lapack-3.10.1

  $ cp make.inc.example make.inc

In the "make.inc" file, there are three lines to be modified:
:: 
  #  CC is the C compiler, normally invoked with options CFLAGS.
  #
  CC     = gcc-9
  ...

  #  should not compile LAPACK with flags such as -ffpe-trap=overflow.
  #
  FORTRAN = gfortran-9
  ...

  #  load options for your machine.
  #
  LOADER   = gfortran-9
  ...

After the modification, we install BLAS and LAPACK libraries:
:: 
  $ make blaslib
  $ make lapacklib

  $ ls lib*
  liblapack.a  librefblas.a

  $ ln -s librefblas.a ./libblas.a

The following information is added in “~/.bash_profile” (or “~/.bashrc”).
:: 
  export LAPACK_PATH=$HOME/Software/lapack-3.10.1

Launch another terminal window or reload "~/.bash_profile" (or "~/.bashrc"):
:: 
  $ source ~/.bash_profile

If you want to uninstall LAPACK/BLAS, just remove the directory "lapack-3.10.1" in "Software".

