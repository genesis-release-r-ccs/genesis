.. highlight:: bash
.. _fitting:

=======================================================================
Fitting section
=======================================================================

Structure fitting
=================

*(In GENESIS 1.1.5 or later only)*
Keywords in the **[FITTING]** section define a structure superimposition scheme,
which is often employed in targeted MD, steered MD, or the String method
(see :ref:`RPATH`) with positional restraint.
In the String method, the reference coordinate for fitting is given by **fitfile**
in the **[INPUT]** section. Otherwise (MD, MIN, REMD), the reference coordinate
is given by **reffile**, **ambreffile**, or **groreffile** in the
**[INPUT]** section. Note that this section is not related to
cryo-EM flexible fitting (see :ref:`Experiments`).

-----------------------------------------------------------------------

**fitting_method** *NO / TR+ROT / XYTR+ZROT*

  **Default: TR+ROT**

  Type of fitting method.

  * NO: No fitting routine is applied
  * TR+ROT: Remove both of translation and rotation
  * XYTR+ZROT: Remove translation in the XY-plane and rotation along the Z-axis 

**fitting_atom**: *Integer*

  **Default: N/A**

  Index of an atom group which is to be fitted to the reference structure.
  In RMSD restraints, Steered MD, or Targeted MD, this should be identical
  to the group where the restraint potential is applied. 
  The index must be defined in **[SELECTION]** (see :ref:`selection`).
  For example, if you specify ``fitting_atom = 1``, the reference atoms 
  are members of ``group1`` in the **[SELECTION]** section.

**mass_weight**: *YES / NO*

  **Default: NO**

  If the parameter is set to *YES*, mass-weighted fitting is employed.
  This parameter should be *YES* for RMSDCOM/PCCOM restraints
  and *NO* for RMSD/PC restraints. Please make sure that this
  parameter is correctly specified when you perform RMSD/RMSDCOM/PC/PCOM calculations.
  In the String method, mass-weighted superimposition is not supported. 

**force_no_fitting**: *YES / NO*

  **Default: NO**

  *This parameter must not be changed for standard MD runs.*
  If the parameter is set to YES and the fitting_method is set to NO, 
  the fitting routine is turned off.
  Translational and rotational fittings are
  usually required to calculate correct RMSD values.
  Therefore, GENESIS simulators (ATDYN and SPDYN) do not allow
  *fitting_method = NO* for simulations involving RMSD calculations
  (targeted/steered MD, for example).
  However, such fitting is not desirable when generating an initial
  structure set for the String method using Cartesian coordinates as CV (see :ref:`RPATH`).
  For this only reason, *fitting_method = NO* was implemented.
  If you want to turn off fittings of RMSD calculations for the preparation of
  an initial structure set for the String method,
  please specify *fitting_method = NO* and *force_no_fitting = YES*.


Examples
========

Example of **[FITTING]** section
::
    
  [FITTING]
  fitting_method = TR+ROT
  fitting_atom   = 1
  mass_weight    = NO

