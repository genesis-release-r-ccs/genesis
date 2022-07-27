.. highlight:: bash
.. _qmmm:

=======================================================================
QMMM section
=======================================================================

Quantum mechanics/Molecular mechanics method (QM/MM)
======================================================

*QM/MM is available only in* **ATDYN**.

The QM/MM method, first proposed in seminal papers by Warshel, Levitt, 
and Karplus :cite:`Warshel:1972` :cite:`Warshel:1976`, is a multi-scale
approach, which treats a partial region of interest (QM region) by 
quantum chemistry, and the surrounding environment (MM region) by 
force field.  The method is valid, particularly when the QM region 
involves an event that cannot be described by the standard force field; 
for example, chemical reactions, spectroscopy, etc.

In the QM/MM method, the potential energy of the system is written as,

   .. math::
      V (\mathbf{R}_a, \mathbf{R}_m) 
      = V^\mathrm{QM} (\mathbf{R}_a, \mathbf{R}_m) 
      + V^\mathrm{QM-MM}_\mathrm{LJ} (\mathbf{R}_a, \mathbf{R}_m)
      + V^\mathrm{MM} (\mathbf{R}_m), 

where :math:`\mathbf{R}_a` and :math:`\mathbf{R}_m` denote the 
positions of atoms in QM and MM regions, respectively.
:math:`V^\mathrm{QM-MM}_\mathrm{LJ}` and :math:`V^{\mathrm{MM}}` 
are the Lennard-Jones interactions between QM-MM atoms and the
force field for MM atoms, respectively.  The QM energy, 
:math:`V^\mathrm{QM}`, is written in 
terms of the electronic energy and the Coulomb interaction 
between nucleus-nucleus and nucleus-MM atoms, 

   .. math::
      V^\mathrm{QM} (\mathbf{R}_a, \mathbf{R}_m) = E_e (\mathbf{R}_a, \mathbf{R}_m) 
      + \sum_{a > a^\prime} \frac{Z_a Z_{a^\prime}}{r_{a{a^\prime}}}
      + \sum_{a, m} \frac{Z_a q_m}{r_{am}},

where :math:`Z_a` and :math:`q_m` are the charge of nucleus and MM atoms, 
respectively, and :math:`r_{a{a^\prime}}` and :math:`r_{am}` denote 
the nucleus-nucleus and nucleus-MM distantces, respectively.
The electronic energy is given by solving the Schrödinger 
equation for electrons,

  .. math::  
     \left[ - \frac{1}{2} \sum_i \nabla_i^2 
     + \sum_{i>j} \frac{1}  {r_{ij}} 
     - \sum_{i,a} \frac{Z_a}{r_{ia}} 
     - \sum_{i,m} \frac{q_m}{r_{im}} \right] \ket{\Psi_e} = E_e \ket{\Psi_e},

where *i*, *a*, and *m* are indices for electrons, nucleus, and MM 
atoms, respectively, and :math:`r_{XY}` denotes the distance between 
particle *X* and *Y*.  

**GENESIS** does not have a function to solve the electronic Schrödinger 
equation but relies on external QM programs, which provide 
the QM energy, its derivatives, and other information.
The interface is currently avaliable for Gaussian, Q-Chem, TeraChem, 
DFTB+, and QSimulate. 

   * `Gaussian09/Gaussian16 <http://gaussian.com>`_

   * `Q-Chem <http://www.q-chem.com>`_

   * `TeraChem <http://www.petachem.com>`_

   * `DFTB+ <https://www.dftbplus.org>`_

   * `QSimulate <https://qsimulate.com/academic>`_

GENESIS/QSimulate is interfaced via shared libraries and seamlessly uses 
the MPI parallelization, thereby facilitating high performace QM/MM-MD 
simulations :cite:`Yagi:2021`. A ready-to-use Singularity image is 
provided by QSimualte Inc.  See `QSimulate <https://qsimulate.com/academic>`_ 
for further information.

Other QM programs are invoked via a system call function of Fortran. In this 
scheme, the input file of a QM calculation is first generated, followed by 
executing a script to run a QM program and reading the information from 
QM output files.  For more information on the method and implementation, 
see Ref.\ :cite:`Yagi:2019`.  Samples of a QM input file
(**qmcnt**) and a script (**qmexe**) are available in our 
`github <https://github.com/yagikiyoshi/QMMMscripts>`_

In order to run QM/MM calculations, users add the **[QMMM]** section
in the control file. Avaliable options are listed in the following.  

-----------------------------------------------------------------------


**qmtyp** *DFTB+ / GAUSSIAN / QCHEM / TERACHEM / QSIMULATE*

  The QM program to use in QM/MM calculations.

**qmatm_select_index** *Integer*

  Index of a group of atoms which is treated as QM atoms. 
  Link hydrogen atoms are automatically added based on a 
  bond connectivity (e.g., given by a PSF file). The index 
  must be defined in **[SELECTION]** (see :ref:`selection`).

**qmcnt**  *Character*

  A template input file of QM calculations.

**qmexe** *Character*

  A script file to execute a QM program.


**workdir** *Character*

  **Default : qmmm**

  The name of a directory where QM input/output files are generated.
  The replica ID is added after this name, e.g., qmmm.0, qmmm.1, etc.

**basename** *Character*

  **Default : N/A**

  The basename of input / output files of QM calculations.

**qmsave_period** *Integer*

  **Default : 1**

  Frequency to save input / output files for QM calculations. 

**savedir** *Character*

  **Default : N/A**

  If present, QM files are copied from **workdir** to this directory. 
  It is typically the case that QM calculations are carried out 
  within a node, and the whole simulation (such as REMD) accross 
  nodes. Then, it is useful for a better performance to set 
  **workdir** to a local disk of each node with fast access (e.g., 
  /dev/shm), and copy the QM files to **savedir** with a frequency
  specified by qmsave_period.

.. Other options:

**qmmaxtrial** *Integer*

  **Default : 1**

  The maximum number of trial run for QM calculations.
  When a QM calculation fails, **GENESIS** repeats the calculation until 
  the iteration reaches this number. The SCF threshold is lowered, if 
  the SCF threshold option is present in the QM control file.

**exclude_charge** *ATOM / GROUP / AMBER* 

  **Default: GROUP**

  This option specifies how to exclude the MM charge in the vicinity of a 
  QM-MM boundary to avoid overpolarization of QM electron density.  When 
  the CHARMM force field is used, *ATOM* excludes only the charge of MM 
  link atom, while *GROUP* excludes the charges of all MM atoms that 
  belongs to the same group as a MM atom at the boundary. When the AMBER
  force field is used, *AMBER* excludes the charge of MM link atom and 
  distributes it to rest of the system evenly.

.. note::
  Since version 1.6.1, QM/MM supports both CHARMM and AMBER force field for 
  the MM. Please keep in mind that **qmmm_generator** is available to only CHARMM.

.. note::
  The QM/MM calculation must be carried out in non-PBC. A non-PBC system 
  can be created from MD trajectory (pdb, dcd) using **qmmm_generator** in 
  the analysis tool. 
  See `the tutorial of QM/MM <https://www.r-ccs.riken.jp/labs/cbrt/tutorials2022/tutorial-15-3>`_
  for more details.


Examples
========

In the following example, the atoms from # 1 to 14 are selected as QM 
atoms by **[SELECTION]** section. The QM program is Gaussian. A 
directory ``qmmm_min`` is created, where input and output files for 
Gaussian (jobXXXX.inp and jobXXXX.log) are saved every 10 steps.
:: 

  [SELECTION]
  group1  = atomno:1-14

  [QMMM]
  qmatm_select_index  = 1
  qmtyp               = gaussian
  qmcnt               = gaussian.com
  qmexe               = runGau.sh
  workdir             = qmmm_min
  basename            = job
  qmsave_period       = 10

The following example is for DFTB+. Because DFTB calculations 
typically take < 1 sec per one snapshot, I/O to generate the 
input and read output could be non-negligible.  It is 
recommended to set **workdir** to a fast disk such as 
``/dev/shm``.  The input and output files of DFTB+ will be 
copied to ``qmmm_min`` every 100 steps.
::

  [QMMM]
  qmatm_select_index = 1
  qmtyp              = dftb+
  qmcnt              = dftb.hsd
  qmexe              = runDFTB.sh
  workdir            = /dev/shm/qmmm_min
  savedir            = qmmm_min
  basename           = job
  qmsave_period      = 100

