.. highlight:: bash
.. _output:

=======================================================================
Output section
=======================================================================

**GENESIS** yields trajectory data (coordinates and velocities)
in the *DCD* file format regardless of the force field or MD algorithm. 
**GENESIS** can also generate a restart file (*rstfile*) 
during or at the end of the simulation,
which can be used to restart and extend the simulation continuously.
Output frequency of each file (e.g., crdout_period and velout_period)
is specified in the **[DYNAMICS]** section in the case of the MD, 
REMD, and RPATH simulations, or **[MINIMIZE]** section in the case of
the energy minimization.

-----------------------------------------------------------------------

General output files
====================

**dcdfile**

  Filename for the coordinates trajectory data.
  Coordinates are written in the DCD format,
  which is commonly used in various MD software such as CHARMM and NAMD.
  The filename must be given in the case of ``crdout_period > 0``.
  However, if ``crdout_period = 0`` is specified in the control file, 
  no **dcdfile** is generated, even if the filename is specified
  in the **[OUTPUT]** section.

**dcdvelfile**

  Filename for the velocity trajectory data.
  Velocities are written in the DCD format.
  The filename must be given in the case of ``velout_period > 0``.
  However, if ``velout_period = 0`` is specified in the control file, 
  no **dcdvelfile** is generated, even if the filename is specified
  in the **[OUTPUT]** section.

**rstfile**

  Filename for the restart data.
  The rstfile contains coordinates, velocities, simulation box size, and so on. 
  This file can be used to extend the simulation continuously.
  In addition, it can be used to switch the simulation algorithms
  (e.g., from minimization to MD, from MD to REMD, from REMD to minimization, etc).
  The filename must be given in the case of ``rstout_period > 0``.
  However, if ``rstout_period = 0`` is specified in the control file,
  no **rstfile** is generated, even if the filename is specified
  in the **[OUTPUT]** section.

**pdbfile** (for ATDYN only)

  Filename for the restart PDB file. This file is updated every ``rstout_period`` steps.


Output files in REMD and RPATH simulations
=====================================================================
When the user performs REMD or RPATH simulations,
the user must include '{}' in the output filename.
This {} is automatically replaced with the replica index.

**remfile** (only for REMD simulations)

  This file contains parameter index data from the REMD simulation, which is
  written for each replica every ``exchange_period`` steps.
  This is used as an input file for the *remd_convert* tool
  to sort the coordinates trajectory data by parameters. 
  The filename must contain '{}', which is automatically replaced with the replica index.
  Note that the information about the parameter index as well as replica index
  in the entire REMD simulation is written in the standard (single) output file
  (see online Tutorials).

**logfile** (only for REMD and RPATH simulations)

  This file contains the energy trajectory data from the REMD or RPATH simulations, which is
  written for each replica every ``exchange_period`` steps.
  This is used as an input file for the *remd_convert* tool
  to sort the coordinates trajectory data by parameters.
  The filename must contain '{}', which is automatically replaced with the replica index.

**rpathfile** (only for RPATH simulations)

  This file contains the trajectory of image coordinates in the string method,
  which are reference values used in the restraint functions.
  Columns correspond to the collective variables, and rows are time steps.
  This data is written with the same timing as the *dcdfile*.
  For details, see :ref:`RPATH`.


Output file in GaMD simulations
===================================

**gamdfile**

  This file provides GaMD parameters determined during the GaMD simulation.
  The filename must be given in the case of ``update_period > 0`` in **[GAMD]** section.
  The GaMD simulation updates its parameters every ``update_period`` steps, and then the updated parameters
  are output to **gamdfile**.
  This file includes the maximum, minimum, average, and deviation of the total potential or dihedral potential,
  which are calculated within the interval ``update_period``.


Output file in Vibrational analysis
===================================

**minfofile**

  This file provides the coordinates and normal mode vectors of the molecules
  specified for vibrational analysis. It is used in **SINDO** for visualizing
  the vibrational motion. It is also an input file to start anharmonic
  vibrational calculations. See :ref:`Vibration` for the vibrational 
  analysis.


Output file in FEP simulations
===================================

**fepfile**

  This file provides energy differences between adjacent windows in FEP simulations.
  The filename must be given in the case of ``fepout_period > 0`` in **[ALCHEMY]** section.
  However, if ``fepout_period = 0`` is specified in the control file, 
  no **fepfile** is generated, even if the filename is specified
  in the **[OUTPUT]** section.
  If FEP/:math:`\lambda`-REMD simulations are performed
  (i.e., both **[REMD]** and **[ALCHEMY]** sections are specified in the control file),
  the filename must contain '{}', which is automatically replaced with the replica index.


Examples
========

For normal MD simulations:
::
  
  [OUTPUT]
  dcdfile = run.dcd
  rstfile = run.rst

For REMD simulations:
::

  [OUTPUT]
  logfile = run_rep{}.log
  dcdfile = run_rep{}.dcd
  remfile = run_rep{}.rem
  rstfile = run_rep{}.rst
