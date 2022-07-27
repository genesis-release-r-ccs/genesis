.. highlight:: bash
.. _introduction:

==============================================================
Introduction
==============================================================

**GENESIS** (*Gen*\ eralized-\ *E*\ nsemble *Si*\ mulation *S*\ ystem)
is a suite of 
computer programs for carrying out molecular dynamics (MD) simulations
of biomolecular systems. MD simulations of biomolecules such as proteins,
nucleic acids, lipid bilayers, N-glycans, are used as
important research tools in structural and molecular biology. Many
useful MD simulation packages :cite:`AMBER18:Program`
:cite:`CHARMM:Program` :cite:`GROMACS:Program` :cite:`DESMOND:Program`
:cite:`NAMD:Program` are now available together with accurate
molecular force field parameter sets :cite:`Case:2005`
:cite:`MacKerell:1998` :cite:`MacKerell:2004` :cite:`Jorgensen:1996`
:cite:`Oostenbrink:2004`. Most of the MD software have been optimized
and parallelized for distributed-memory parallel supercomputers or
PC-clusters. Therefore hundreds of CPUs or CPU cores can be used
efficiently for a single MD simulation of a relatively large
biomolecular system, typically composed of several hundred
thousands of atoms. In recent years, the number of available CPUs
or CPU cores is rapidly increasing. The implmentation of highly
efficient parallel schemes is therefore required in modern MD simulation
programs. Accelerators such as GPGPU (General-Purpose
computing on Graphics Processing Units) also become popular,
and thus their utilization is also desired.
Actually, many MD program packages support various accelerators.

Our major motivation is to develop MD simulation software
with a scalable performance on such modern supercomputers. For this
purpose, we have developed the software from scratch, introducing the
hybrid (MPI + OpenMP) parallelism, several new parallel
algorithms :cite:`Jung:2013` :cite:`Jung:2014`, and GPGPU calculation.
Another motivation is to develop a simple MD program, which can be
easily understood and modified for methodological developments. These two
policies (high parallel performance and simplicity) usually conflict
each other in computer software. 
To avoid the conflict, we have developed two MD programs in **GENESIS**, 
namely **SPDYN** (*Sp*\ atial decomposition *dyn*\ amics)
and **ATDYN** (*At*\ omic decomposition *dyn*\ amics).

**SPDYN** and **ATDYN** share almost the same data structures, subroutines,
and modules, but differ in their parallelization schemes.
In **SPDYN**, the spatial decomposition scheme is implemented with new
parallel algorithms :cite:`Jung:2013` :cite:`Jung:2014` and GPGPU calculation.
In **ATDYN**, the atomic decomposition scheme is introduced for simplicity.
The performance of **ATDYN** is not comparable to **SPDYN** due to the
simple parallelization scheme.
However, **ATDYN** is easier to modify for development of new algorithms
or novel molecular models. We hope that users develop new methodologies
in **ATDYN** at first and, eventually, port them to **SPDYN** for the
better performance.  As we maintain consistency between the source codes
of **ATDYN** and **SPDYN**, switching from **ATDYN** to **SPDYN** is
not quite hard.

Other features in **GENESIS** are listed below:

* Not only atomistic molecular force field (CHARMM, AMBER) but
  also some coarse-grained models are available in **ATDYN**.

* For extremely large biomolecular systems (more than 10 million
  atoms), parallel input/output (I/O) scheme is implemented.

* **GENESIS** is optimized for K and Fugaku computer (developed by RIKEN and
  Fujitsu company), but it is also available on Intel-based supercomputers
  and PC-clusters.

* **GENESIS** is written in modern Fortran language (90/95/2003)
  using modules and dynamic memory allocation. No common blocks are used.

* **GENESIS** is free software under the GNU Lesser General Public License (LGPL)
  version 3 or later. We allow users to use/modify **GENESIS** and
  redistribute the modified version under the same license.

This user manual mainly provides detailed description of 
keywords used in the control file. Tutorials for standard MD
simulations, REMD simulations, and some analyses are
`available online <https://www.r-ccs.riken.jp/labs/cbrt/>`_.
We recommend new users of **GENESIS** to start from 
the next chapter to learn a basic idea, installation, and work flow of the program.

Comparing to other MD software, e.g. AMBER, CHARMM, or NAMD,
**GENESIS** is a very young MD simulation program. Before releasing the
program, the developers and contributors in **GENESIS** development team
worked hard to fix all bugs in the program, and performed a bunch of
test simulations. Still, there might be defects or bugs in
**GENESIS**. Since we cannot bear any responsibility for the
simulation results produced by **GENESIS**, we strongly recommend 
the users to check the results carefully.

The **GENESIS** development team has a rich plan for future development
of methodology and molecular models. We would like to
make **GENESIS** one of the most powerful and feasible MD software packages,
contributing to computational chemistry and biophysics. Computational
studies in life science is still at a very early stage (like
‘GENESIS’) compared to established experimental researches. We hope
that **GENESIS** pushes forward the computational science and
contribute to bio-tech and medical applications in the future.

