.. highlight:: bash
.. _constraints:

=======================================================================
Constraints section
=======================================================================

SHAKE/RATTLE algorithms
=======================

In the **[CONSTRAINTS]** section, keywords related to bond constraints are specified.
In the leapfrog integrator, the SHAKE algorithm is applied for
covalent bonds involving hydrogen :cite:`Rycaert:1997`.
In the velocity Verlet and multiple time-step integrators,
not only SHAKE but also RATTLE is used :cite:`Andersen:1983`.
Note that bond constraint between heavy atoms is currently not available.

-----------------------------------------------------------------------

**rigid_bond** *YES / NO*

  **Default : NO**

  Turn on or off the SHAKE/RATTLE algorithms for covalent bonds involving hydrogen.

**shake_iteration** *Integer*

  **Default : 500**

  Maximum number of iterations for SHAKE/RATTLE constraint.
  If SHAKE/RATTLE does not converge within the given number of iterations,
  the program terminates with an error message.

**shake_tolerance** *Real*

  **Default : 1.0e-10** (unit : :math:`\text{\AA}`)

  Tolerance of SHAKE/RATTLE convergence.

**hydrogen_type** *NAME / MASS*

  **Default : NAME**

  This parameter defines how hydrogen atoms are detected. This parameter
  is ignored when *rigid_bond = NO*. Usually, the user does not need to 
  take care about this parameter.

  * **MASS** : detect hydrogens based on atomic mass. If the mass
    of an atom is less than *hydrogen_mass_upper_bound* and greater than 0, 
    that atom is considered as a hydrogen.

  * **NAME** : detect hydrogens based on atom name, type, and mass. If the mass
    of an atom is less than *hydrogen_mass_upper_bound* and the name
    or type begins with 'h', 'H', 'd', or 'D', that atom is considered as
    a hydrogen.  
    
  +------------------+------+------------+--------+
  | atom name (type) | mass | *NAME*     | *MASS* |
  +==================+======+============+========+
  | HX               | 1.0  | o          | o      |
  +------------------+------+------------+--------+
  | XX               | 1.0  | x          | o      |
  +------------------+------+------------+--------+
  | HY               | 3.0  | x          | x      |
  +------------------+------+------------+--------+
  | YY               | 3.0  | x          | x      |
  +------------------+------+------------+--------+

  o: treated as hydrogen, x: not treated as hydrogen. Here, we assumed *hydrogen_mass_upper_bound* 2.1.

**hydrogen_mass_upper_bound** *Real*

  **Default : 2.1**

  This parameter defines the upper limit of atomic mass to define a hydrogen atom. 
  For exmaple, if 3.0 is used, atoms with an atomic mass of less than 3.0 are treated as hydrogens.
  You should write it in the case of hydrogen mass repartitioning scheme.
  This option must be used in **GENESIS 1.2** or later.

**noshake_index** *Integer*

  **Default : N/A**

  This is the index of the atom where the constraint is NOT applied. 
  The index must be defined in **[SELECTION]** (see :ref: `selection`).
  For example, if you specify ``select_index1 = 1``, this constraint function
  s NOT applied for ``group1`` in the **[SELECTION]** section.

SETTLE algorithm
================

**fast_water** *YES / NO*

  **Default : YES**

  Turn on or off the SETTLE algorithm for the constraints of the water molecules :cite:`Miyamoto:1992`.
  Although the default is "fast_water=YES", the user must specify "rigid_bond=YES" to use the SETTLE algorithm.
  If "rigid_bond=YES" and "fast_water=NO" are specified, the SHAKE/RATTLE algorithm is
  applied to water molecules, which is computationally inefficient.

**water_model** *expression or NONE*

  **Default : TIP3**

  Residue name of the water molecule to be rigidified in the SETTLE algorithm.
  In the case of the AMBER force field, "water_model = WAT" must be specified.


.. note::
  TIP4P water model is availabe in **GENESIS 1.2** or later.
  In the case of using TIP4P water model, it is regarded as rigid. In
  molecular dynamics simulations, please use **rigid_bond** and **fast_water**
  "YES". In minimization, **[Constraints]** has not been used before, but now
  "**fast_water** = YES" can be used when TIP4P water model is used.
  However, please keep in mind that
  other parameters cannot be defined in minimizations, and constraints are not
  applied except for water molecules. TIP4P water model can be used only in SPDYN.


LINCS algorithm
===============

**fast_bond** *YES / NO* (**LEAP** integrator in **ATDYN** only)

  **Default : NO**

  Turn on or off the LINCS algorithm.
  To use the LINCS algorithm, "rigid_bond=YES" should be also specified.

**lincs_iteration** *Integer* (**ATDYN** only)

  **Default : 1**

  Number of iterations in the LINCS algorithm.

**lincs_order** *Integer* (**ATDYN** only)

  **Default : 4**

  Matrix expansion order in the LINCS algorithm.


Examples
========

In the case of the CHARMM force field:
:: 
  [CONSTRAINTS]
  rigid_bond  = YES   # Turn on SHAKE/RATTLE
  fast_water  = YES   # Turn on SETTLE

In the case of the AMBER force field:
:: 
  [CONSTRAINTS]
  rigid_bond  = YES   # Turn on SHAKE/RATTLE
  fast_water  = YES   # Turn on SETTLE
  water_model = WAT   # residue name of the rigid water

Turn off all constraints in the system
:: 
  [CONSTRAINTS]
  rigid_bond  = NO
  fast_water  = NO
