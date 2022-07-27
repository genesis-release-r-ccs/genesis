.. highlight:: bash
.. _restraints:

=======================================================================
Restraints section
=======================================================================

Restraint potential
===================

The **[RESTRAINTS]** section contains keywords to define external
restraint functions. The restraint functions are applied to the
selected atom groups in the **[SELECTION]** section to restrict the
motions of those atoms.

The potential energy of a restaint can be written as:

  .. raw:: latex
     
     \vspace{-5mm}

  .. math::  
     U(x) = k\;(x-x_{0})^{n}

  .. raw:: latex
     
     \vspace{-3mm}

where :math:`x` is a variable (see bellow), :math:`x_0` is a reference value,
:math:`k` is a force constant, and :math:`n` is an exponent factor.

If a repulsive distance restraint is applied,
the above restraint function is modified as following:

  .. raw:: latex
     
     \vspace{-5mm}

  .. math::  
     U(x) = 
     \begin{cases}
     k\;(x-x_{0})^{n} & (x < x_{0})\\
     0 & (x \ge x_{0})
     \end {cases}.

  .. raw:: latex
     
     \vspace{-3mm}

Or if a flat-bottom distance restraint is applied,
the function is modified as following:

  .. raw:: latex
     
     \vspace{-5mm}

  .. math::  
     U(x) = 
     \begin{cases}
     k\;(x-x_{0})^{n} & (x > x_{0})\\
     0 & (x \le x_{0})
     \end {cases}.

  .. raw:: latex
     
     \vspace{-3mm}

-----------------------------------------------------------------------

General keywords
----------------

**nfunctions** *Integer*

  **Default: 0**

  Number of restraint functions.
  
**function**:math:`\textbf{\textit{N}}` *POSI / DIST / DISTMASS / ANGLE / ANGLEMASS / DIHED / DIHEDMASS / RMSD /RMSDMASS / PC / PCMASS / REPUL / REPULMASS / FB / FBMASS / EM*

  **Default: N/A**

  Type of restraint.

  * **POSI**: positional restraint. The reference coordinates are given
    by **reffile**, **ambreffile**, or **groreffile** in **[INPUT]**.
    (see :ref:`input`)

  * **DIST or DISTMASS**: distance restraint. 

  * **ANGLE or ANGLEMASS**: angle restraint. 

  * **DIHED or DIHEDMASS**: dihedral angle restraint. 

  * **RMSD or RMSDMASS**: RMSD restraint. *MASS* means mass-weighted RMSD.
    Translational and rotational fittings to the reference coordinate
    are done before calculating RMSD. The reference coodinate is specified
    in the same manner as *POSI*.

  * **PC or PCMASS**: principal component constraint. This option
    requires *modefile* in the :ref:`input`.

  * **REPUL or REPULMASS**: repulsive distance restraint. 

  * **FB or FBMASS**: flat-bottom distance restraint. 

  * **EM**: cryo-EM flexible fitting (see :ref:`Experiments`)

  *DIST*, *ANGLE*, *DIHED* impose restraints on a distance/angle/dihedral
  defined by groups in the **[SELECTION]** section specified by 
  ``select_indexN`` (see examples below).
  *REPUL/FB* also impose restraint on the distance defined by the selected
  groups, but the restraint is imposed to only one side of the harmonic
  potential centered on the given distance.
  *MASS* indicates that the restraint force is applied to the center of mass 
  of the selected group.
  When *MASS* is omitted, the force is applied to the geometric center of
  the coordinates.
  The *MASS* keyword does nothing for groups consisting of a single particle.

  In **SPDYN**, *POSI* and *RMSD[MASS]* restraints are mutually exclusive;
  you can use either one or none of them. Additionally, two different *POSI*
  restraints might not be applied simultaneously, either.

  *Note: POSI, PC, and RMSD restraints can be influenced by the removal
  of translational/rotational momentum. See also the stoptr_period parameter
  in the [DYNAMICS] section for more information.*

**constant**:math:`\textbf{\textit{N}}` *Real*

  **Default: 0.0** (unit: depend on the restraint type)

  Force constant of a restraint function.
  The unit depends on the type of restraint.
  Namely, :math:`{\rm{kcal/mol/{\text{\AA}}^{\textit{n}}}}` is used in the case of *DIST*, *RMSD*, *REPUL*, and *FB*,
  while :math:`{\rm{kcal/mol/{rad}^{\textit{n}}}}` in the case of *ANGLE* and *DIHED*,
  where :math:`n` is **exponent**:math:`\textbf{\textit{N}}` specified in this section.


**reference**:math:`\textbf{\textit{N}}` *Real*

  **Default: 0.0** (unit: depend on the restraint type)
  
  Reference value of a restraint function. For the positional restraint, the value is ignored.
  The unit depends the type of restraint.
  Namely, angstroms are used in the case of *DIST*, *REPUL*, and *FB*,
  while degrees (NOT radians) are used in the case of *ANGLE* and *DIHED*.

**select_index**:math:`\textbf{\textit{N}}` *Integer*

  **Default: N/A**

  Index of a group, to which restraint potentials are applied.
  The index must be defined in **[SELECTION]** (see :ref:`selection`).
  For example, if you specify ``select_index1 = 1``, this restraint function
  is applied for ``group1`` in the **[SELECTION]** section.

  The number of groups required depends on the type of the restraint function.

  * *POSI/RMSD[MASS]*: 1
  * *DIST[MASS]*: 2
  * *ANGLE[MASS]*: 3
  * *DIHED[MASS]*: 4
  * *PC[MASS]*: :math:`\ge 1`
  * *REPUL[MASS]*: :math:`\ge 2`
  * *FB[MASS]*: :math:`\ge 2`

  A group can contain more than a single atom.
  Suppose we have the following input.
  
  ::

    [SELECTION]
    group1        = ai:1-10
    group2        = ai:11-20

    [RESTRAINTS]
    nfunctions    = 1
    function1     = DIST
    constant1     = 3.0
    reference1    = 10.0
    select_index1 = 1 2

  In this case, the distance restraint is applied to the distance between
  the geometric centers of group1 and group2. The calculated force is then
  distributed to each atom. If *DISTMASS* is given instead of *DIST*,
  mass centers (mass-weighted average positions) are used instead of
  geometric centers (non-weighted average position).

  *REPUL[MASS]/FB[MASS]* restraint can be simultaneuosly applied to more than 2 groups.
  In the case, the restraint forces are applied to all pairs in the groups.
  See example inputs in the end of this section.

**direction**:math:`\textbf{\textit{N}}` *ALL / X / Y / Z*

  **Default : ALL**

  Direction of the *POSI* restraint. If *X* or *Y* or *Z* is specified, restraints
  along the other two axes are not applied.

**exponent**:math:`\textbf{\textit{N}}` *Integer*

  **Default : 2**

  Exponent factor of the restraint function. The default is harmonic (2).
  This parameter does not work for *POSI* and *RMSD[MASS]* restraints
  in **SPDYN**, where the default value, 2, is always used.

**mode**:math:`\textbf{\textit{N}}` *Integer*

  Specifies the mode index which is used for the PC (principal component)
  restraint. For example, the 1st PC mode can be restrained by specifying
  ``mode1=1``. 

Restraints of a linear combination of distances
-----------------------------------------------

A linear combination of :math:`m` distances, :math:`r_i (i = 1, 2, ..., m)`, can 
be restrained by specifying :math:`2m` groups in *DIST[MASS]*. The combination of distances
is expressed as,

   .. math::
      r_{\rm sum} = \sum_{i=1}^m{w_i|r_i|^{n_i}}.

where :math:`n_i` and :math:`w_i` are specified by **exponent_dist**:math:`\textbf{\textit{N}}`
and **weight_dist**:math:`\textbf{\textit{N}}`, respectively.

**exponent_dist**:math:`\textbf{\textit{N}}` *Integer* (**DIST[MASS]** only)

  **Default : 1**

  The exponent factors, :math:`n_i`, of :math:`r_{\rm sum}`. Note that there must be :math:`m` entries
  of integers.

**weight_dist**:math:`\textbf{\textit{N}}` *Real* (**DIST[MASS]** only)

  **Default : 1.0**

  The coefficients, :math:`w_i`, of :math:`r_{\rm sum}`. Note that there must be :math:`m` entries
  of real numbers.

:: 
  
  [RESTRAINTS]

  # |r12| - |r13|
  nfunctions       = 1
  function1        =  DIST
  constant1        =  10.0 
  reference1       =  -1.0  
  exponent_dist1   = 1    1
  weight_dist1     = 1.0 -1.0
  select_index1    = 1 2  1 3   


Pressure derived from restraints
--------------------------------

The users can select whether or not to include the pressure calculated from a restraint potential 
in the total internal pressure, which is the target pressure of a simulation with the NPT ensemble,
using the keywords **pressure_position** and **pressure_rmsd**.
By default, the pressure derived from positional and RMSD restraints are treated as an external pressure.
However, if simulations with POSI or RMSD restraint show a strange behaviour (unphysical box deformation), 
especially, when strong force constants are applied, it is recommended to turn on these options.

**pressure_position** *YES / NO*

  **Default : NO**

  If YES, the virial terms from positional restraints are included in pressure evaluation.

**pressure_rmsd** *YES / NO*

  **Default : NO**

  If YES, the virial terms from RMSD restraints are included in pressure evaluation.


Advanced definition of restraints
---------------------------------

Restraints can be also defined in an external input file (**localresfile**).
In this case, the number of local restraints must NOT be included in **nfunctions**.
This option is availabe in **SPDYN** only.
For details, see :ref:`input`.


Restraints in REUS simulations
------------------------------

If you employ a certain restraint term for REUS runs, *nreplica*
force constants and reference values must be given as a space-separated list.
The above keywords, except for **nfunctions**, **pressure_position**,
and **pressure_rmsd**, must have a serial number,
'N', of the function (:math:`N \geq 1`).
This serial number is referred when selecting restraints in REUS runs.
For details, see :ref:`REMD`.


Examples
========

Example of **[RESTRAINTS]** section:
:: 
  
  [RESTRAINTS]
  nfunctions    = 1
  function1     = DIST
  reference1    = 10.0
  constant1     = 2.0
  select_index1 = 1 2        # group1 and group2 in [SELECTION]


Example of multiple restraints:
:: 

  [RESTRAINTS]
  nfunctions    = 2
  
  function1     = DIST
  constant1     = 2.0
  reference1    = 10.0       # in angstroms
  select_index1 = 1 2

  function2     = DIHED
  constant2     = 3.0
  reference2    = 120.0      # in degrees
  select_index2 = 3 4 5 6


Example of repulsive restraints. The repulsive restraint is simultaneously
applied to the centers of mass (COMs) of groups 1, 2, and 3. In the case, 
restraint forces are applied to three distances of the COMs between 1 and 2, 
those between 1 and 3, and those between 2 and 3.
Groups 1, 2, and 3 cannot approach each other.
:: 

  [RESTRAINTS]
  nfunctions    = 1

  function1     = REPULMASS
  reference1    = 10.0
  constant1     = 1.0
  select_index1 = 1 2 3

