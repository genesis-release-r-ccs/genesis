.. highlight:: bash
.. _reaccoord:

=======================================================================
ReacCoord section
=======================================================================

Reaction coordinates for ABMD
=============================

The **[ReacCoord]** section contains keywords to define reaction
coordinates for ABMD and its extensions. The bias potential are 
applied to the selected atom groups in the **[SELECTION]** section 
to enchance the motions of those atoms.

-----------------------------------------------------------------------

General keywords
----------------

**nreac_coord** *Integer*

  **Default: 0**

  Number of bias potentials.
  
**function**:math:`\textbf{\textit{N}}` *DIST / DISTMASS / ANGLE / ANGLEMASS / DIHED / DIHEDMASS*

  **Default: N/A**

  Type of bias potential.

  * **DIST or DISTMASS**: distance restraint. 

  * **ANGLE or ANGLEMASS**: angle restraint. 

  * **DIHED or DIHEDMASS**: dihedral angle restraint. 

  *DIST*, *ANGLE*, *DIHED* impose bias potential on a distance/angle/dihedral
  defined by groups in the **[SELECTION]** section specified by 
  ``select_indexN`` (see examples below).
  *MASS* indicates that the restraint force is applied to the center of mass 
  of the selected group.
  When *MASS* is omitted, the force is applied to the geometric center of
  the coordinates.
  The *MASS* keyword does nothing for groups consisting of a single particle.

**monitor_period** *Integer*

  **Default: 0**

  Output frequency for the reaction coordinates. The reaction coordinates 
  are written in the reaccoord file every **monitor_period** steps during the simulation.
  For example, if "timestep=0.001" and "monitor_period=1000" are specified, 
  the reaction coordinattes are written every 1 ps.

**nxi** *Integer*

  **Default : 0**

  Number of bins along reaction coordinates.

**ximin** *Real*

  **Default : 0.0**

  The minimum velues of bins along reaction coordinates.

**ximax** *Real*

  **Default : 0.0**

  The macimum velues of bins along reaction coordinates.  

**periodicity**  *YES / NO*

  **Default : NO**

  Select Yes if the CV is periodic.

**select_index**:math:`\textbf{\textit{N}}` *Integer*

  **Default: N/A**

  Index of a group, to which bias potentials are applied.
  The index must be defined in **[SELECTION]** (see :ref:`selection`).
  For example, if you specify ``select_index1 = 1``, this restraint function
  is applied for ``group1`` in the **[SELECTION]** section.

  The number of groups required depends on the type of the restraint function.

  * *DIST[MASS]*: 2
  * *ANGLE[MASS]*: 3
  * *DIHED[MASS]*: 4

  A group can contain more than a single atom.
  Suppose we have the following input.
  
  ::

    [SELECTION]
    group1        = ai:1-10
    group2        = ai:11-20

    [ReacCoord]
    function1        = DIST
    select_index1    = 1 2
    nxi1             = 51
    ximin1           =  5.0
    ximax1           = 20.0
    periodicity1     = NO

  In this case, the distance restraint is applied to the distance between
  the geometric centers of group1 and group2. The calculated force is then
  distributed to each atom. If *DISTMASS* is given instead of *DIST*,
  mass centers (mass-weighted average positions) are used instead of
  geometric centers (non-weighted average position).

  *REPUL[MASS]/FB[MASS]* restraint can be simultaneuosly applied to more than 2 groups.
  In the case, the restraint forces are applied to all pairs in the groups.
  See example inputs in the end of this section.

Restraints of a linear combination of distances
-----------------------------------------------

A linear combination of :math:`m` distances, :math:`r_i (i = 1, 2, ..., m)`, can 
be restrained by specifying :math:`2m` groups in *DIST[MASS]*. The combination of distances
is expressed as,

   .. math::
      r_{\rm sum} = \sum_{i=1}^m{w_i|r_i|^{n_i}}.

where :math:`n_i` and :math:`w_i` are specified by **exponent_dist**:math:`\textbf{\textit{N}}`
and **weight_dist**:math:`\textbf{\textit{N}}`, respectively.

**weight_dist**:math:`\textbf{\textit{N}}` *Real* (**DIST[MASS]** only)

  **Default : 1.0**

  The coefficients, :math:`w_i`, of :math:`r_{\rm sum}`. Note that there must be :math:`m` entries
  of real numbers.

:: 
  
  [ReacCoord]

  # |r12| - |r13|
  function1        = DIST
  select_index1    = 1 2 1 3
  nxi1             = 51
  ximin1           = -1.2
  ximax1           =  1.2
  weight_dist1     = 1.0 -1.0


Examples
========

Example of **[ReacCoord]** section:
:: 
  
  [ReacCoord]
  nreac_coord      = 1
  monitor_period   = 100000

  function1        = DIST
  select_index1    = 1 2
  nxi1             = 51
  ximin1           =  5.0
  ximax1           = 20.0
  periodicity1     = NO


Example of multiple restraints:
:: 

  [ReacCoord]
  nreac_coord      = 2
  monitor_period   = 100000

  function1        = DIST
  select_index1    = 1 2
  nxi1             = 51
  ximin1           =  5.0
  ximax1           = 20.0
  periodicity1     = NO

  function2        = ANGLE
  select_index2    = 3 4
  nxi2             = 51
  ximin2           = 20.0
  ximax2           = 60.0
  periodicity2     = NO

