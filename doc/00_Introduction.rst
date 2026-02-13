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

Our primary motivation is to develop MD simulation software that delivers
scalable performance on modern supercomputers. To achieve this, we developed
GENESIS from scratch, incorporating hybrid (MPI + OpenMP) parallelism, several
original parallel algorithms :cite:`Jung:2013` :cite:`Jung:2014`, and GPGPU
acceleration. Another key objective is to maintain a codebase that is simple
enough to be easily understood and modified for methodological developments.
Since high parallel performance and code simplicity often conflict in software
engineering, we address this challenge by providing three distinct MD engines in
**GENESIS**: **SPDYN** (*Sp*\ atial decomposition *dyn*\ amics)
**ATDYN** (*At*\ omic decomposition *dyn*\ amics), and 
**CGDYN** (*Coarse-Grained* molecular *dyn*\ amics).

**SPDYN**, **CGDYN**, and **ATDYN** share many common data structures,
subroutines, and modules, but diverge significantly in their parallelization
strategies to target specific needs. **SPDYN** and **CGDYN** are dedicated to
high-performance simulations using spatial decomposition techniques: **SPDYN**
employs a static decomposition scheme optimized for all-atom systems
:cite:`Jung:2013` :cite:`Jung:2014`, while **CGDYN** utilizes an advanced
dynamic domain decomposition scheme to ensure efficient load balancing in
coarse-grained (CG) simulations :cite:`Jung:2024`. In contrast, **ATDYN** adopts
an atomic decomposition scheme, prioritizing algorithmic simplicity over raw
performance. Although **ATDYN** cannot match the scalability of the domain
decomposition engines, its straightforward architecture makes it the ideal
platform for developing new algorithms or novel molecular models. We encourage
users to prototype new methodologies in **ATDYN** and subsequently port them to
**SPDYN** or **CGDYN** for production runs. Thanks to the strict consistency
maintained across the source codes of all three engines, transitioning from the
development phase (**ATDYN**) to high-performance execution (**SPDYN** or
**CGDYN**) is seamless.

Other features in **GENESIS** are listed below:

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
`available online <https://mdgenesis.org/>`_.
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

