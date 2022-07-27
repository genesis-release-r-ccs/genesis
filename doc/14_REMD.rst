.. highlight:: bash
.. _remd:

=======================================================================
REMD section
=======================================================================

Replica-exchange molecular-dynamics simulation (REMD)
=====================================================

In the **[REMD]** section, the users can specify keywords for Replica-Exchange
Molecular Dynamics (REMD) simulation.
REMD method is one of the enhanced conformational sampling methods 
used for systems with rugged free-energy landscapes.
The original temperature-exchange method (T-REMD) is one of the most
widely used methods in biomolecules' simulations
:cite:`Sugita:1999` :cite:`Mitsutake:2001`.
Here, replicas (or copies) of the original system are prepared,
and different temperatures are assigned to each replica. 
Each replica runs in a canonical (NVT) or isobaric-isothermal (NPT)
ensemble, and the temperatures are periodically exchanged between
the neighboring replicas during a simulation. 
Exchanging temperature enforces a random walk in temperature space,
allowing the system to overcome energy barriers and sample a much wider
conformational space.

In REMD methods, the transition probability of the replica exchange process is given by the usual Metropolis criterion,

  .. raw:: latex

     \vspace{-5mm}

  .. math::  
     w(X \to X') & = \mathrm{min}(1, \frac{P(X')}{P(X)}) = \mathrm{min}(1,\mathrm{exp}(-\Delta)).


  .. raw:: latex

     \vspace{-3mm}

In the T-REMD method, we have

  .. raw:: latex

     \vspace{-5mm}

  .. math:: 
     \Delta = (\beta_{m} - \beta_{n})\left\{ {E(q^{[j]}) -E(q^{[i]})} \right\},

  .. raw:: latex

     \vspace{-3mm}

where :math:`E` is the potential energy, :math:`q` is the position of atoms,
:math:`\beta` is the inverse temperature defined by :math:`\beta = 1/k_B T`,
:math:`i` and :math:`j` are the replica indexes, and :math:`m` and :math:`n`
are the parameter indexes.
After the replica exchange, atomic momenta are rescaled as follows:

  .. raw:: latex

     \vspace{-5mm}

  .. math:: 
     p^{[i]'} = \sqrt{\frac{T_{n}}{T_{m}}} p^{[i]}, \qquad p^{[j]'} = \sqrt{\frac{T_{m}}{T_{n}}} p^{[j]},

  .. raw:: latex

     \vspace{-3mm}

where :math:`T` is the temperature and :math:`p` is the momenta of atoms.

The transition probability should be independent of the algorithms used,
i.e. constant temperature and constant pressure algorithms.
On the other hand, the momenta-rescaling scheme depends on the algorithm
used in the simulation. 
If thermostat and barostat momenta are included in the equations of motion,
these variables should be also rescaled after replica exchange
:cite:`Mori:2010-1` :cite:`Mori:2010-2`.
In GENESIS, barostat momentum is rescaled in the case of T-REMD with
Langevin or Bussi method in NPT, NPAT, and NPgT ensembles.
In other cases, only atomic momenta are rescaled.

In GENESIS, not only Temperature REMD but also pressure REMD :cite:`Okabe:2001`,
surface-tension REMD :cite:`Mori:2013`, REUS (or Hamiltonian REMD) :cite:`Sugita:2000` :cite:`Fukunishi:2002`, replica exchange with solute tempering (REST) :cite:`Terakawa:REST2` :cite:`Kamiya:gREST`,
and their multi-dimensional combinations are available in both **ATDYN** and **SPDYN**.
Basically, these methods can be employed in the NVT, NPT, NPAT, NPgT ensembles,
except for the surface-tension REMD, which is only used in the NPgT ensemble.
REMD simulations in GENESIS require an MPI environment. 
At least one MPI process must be assigned to one replica.
For example, when the user wants to employ 32 replicas,
:math:`32n` MPI processes are required.

In the following parameters excluding *dimension*, *exchange_period*,
and *iseed*, the last character *N* must be replaced with a positive
integer number (i.e. :math:`N \geq 1`), which defines the index of
the replica dimension. For example, *type1*, *nreplica1* are the replica type
and number of replicas for the first dimension, respectively.
For details, see the examples below.

-----------------------------------------------------------------------

**dimension** *Integer*

  **Default: 1**

  Number of dimensions (i.e. number of parameter types to be exchanged) 

**type**:math:`\textbf{\textit{N}}` *TEMPERATURE / PRESSURE / GAMMA / RESTRAINT / REST*

  **Default: TEMPERATURE**

  Type of parameter to be exchanged in the \ :math:`N`\-th dimension

  * **TEMPERATURE**: Temperature REMD :cite:`Sugita:1999`
  * **PRESSURE**: Pressure REMD :cite:`Okabe:2001`
  * **GAMMA**: Surface-tension REMD :cite:`Mori:2013`
  * **RESTRAINT**: REUS (or Hamiltonian REMD) :cite:`Sugita:2000` :cite:`Fukunishi:2002`
  * **REST**: replica exchange with solute tempering (REST2 or gREST) :cite:`Terakawa:REST2` :cite:`Kamiya:gREST`, which is totally different from the original version of REST :cite:`Berne:REST`. Currently, only AMBER and CHARMM force fields are supported.
  * **ALCHEMY**: FEP/:math:`\lambda`-REMD :cite:`Jiang:2009` 

**nreplica**:math:`\textbf{\textit{N}}` *Integer*

  **Default: 0**

  Number of replicas (or parameters) in the \ :math:`N`\-th dimension

**parameters**:math:`\textbf{\textit{N}}` *Real*

  **Default: N/A**

  List of parameters for each replica in the \ :math:`N`\-th dimension.
  Parameters must be given as a space-separated list,
  and the total number of parameters must be equal to **nreplicaN**:math:`\textbf{\textit{N}}`.
  In case of REUS (type = RESTRAINT), parameters must be 
  specified in **[RESTRAINTS]** section (see the sample below).
  In case of gREST (type = REST), these parameters are considered 
  as temperature of solute region.
  Note that the order of the parameters in this list must NOT be changed
  before and after the restart run, even if the parameters are exchanged
  during the REMD simulation.

**exchange_period** *Integer*

  **Default: 100**

  Frequency of the parameter exchange attempt. 
  If "exchange_period = 0" is specified, 
  REMD simulation is carried out without parameter exchange,
  which is useful to equilibrate the system in a condition
  assinged to each replica before performing the production run.

**cyclic_params**:math:`\textbf{\textit{N}}` *YES / NO*

  **Default: NO**

  Turn on or off the periodicity of the parameters in the \ :math:`N`\-th dimension.
  If "cyclic_paramsN = YES" is specified, the first and last parameters are
  considered as neighbouring parameters.
  This option can be applicabe to all parameter types.
  Basically, this is useful in the case of REUS in dihedral angle space,
  since the dihedral angle is a periodic variable.

**iseed** *Integer*

  **Default: 3141592**

  Random number seed in the replica exchange scheme.
  If this is not specified explicitly, iseed is taken over from the restart file.

.. note::
   In the [ENSEMBLE] section, there is also a parameter "temperature".
   In the T-REMD simulation, this temperature is ignored, even if it is specified explicitly.
   Similarly, pressure and gamma in the [ENSEMBLE] section are ignored in the P-REMD 
   and surface-tension REMD simulations, respectively.

.. note::
   When multi-dimensional REMD is carried out, parameters are exchanged alternatively.
   For example, in TP-REMD (type1 = TEMPERATURE and type2 = PRESSURE),
   there is a temperature exchange first, followed by a pressure exchange.
   This is repeated during the simulations.

.. note::
  **ALCHEMY** is not implemented in GENESIS 2.0.0. Please use GENESIS 1.7.0 or later version for this purpose.

Replica-exchange umbrella-sampling (REUS)
=========================================

**rest_function**:math:`\textbf{\textit{N}}` (for **REUS** only)

  Index of the restraint function to be used in the REUS simulation.
  The detailed parameters in the restraint function (e.g., force constant and reference)
  are defined in the **[RESTRAINTS]** section (see :ref:`restraints`).
  Note that the order of the parameters in the **[RESTRAINTS]** section 
  must NOT be changed before and after the restart run, 
  even if the parameters are exchanged during the REUS simulation.

  **GENESIS** supports not only on-grid but also off-grid schemes.
  In the off-grid REUS, multiple restraints are merged into a single reaction coordinate (see example below).
  Those restraints are defined in **[RESTRAINTS]** section,
  where the number of parameters (const and reference) must be equal to *nreplicaN*.
  Note that this kinds of combined axis can be used only for restraints,
  other types (such as combined tempreature-pressure or temperature-restraint
  coordinate) are not available currently.

.. note::
   Positional restraint is not available for REUS.
   In **SPDYN**, PCA restraint is not available for REUS.
   The control file format was completely changed after verion 1.1.0,
   since the off-grid REUS scheme was introduced.
   When the users use the old control file, please be careful.


Replica-exchange with solute-tempering (gREST)
==============================================

**select_index**:math:`\textbf{\textit{N}}` *Integer*

  **Default: N/A**

  Index of an atom group. The selected atoms are considered as "solute" in gREST. 
  The index must be defined in **[SELECTION]** (see :ref:`selection`).

**param_type**:math:`\textbf{\textit{N}}` *ALL / BOND / ANGLE / UREY / DIHEDRAL / IMPROPER / CMAP / CHARGE / LJ*

  **Default: ALL**

  Solute energy terms for gREST :cite:`Kamiya:gREST` simulations.
  Energy terms selected by this parameter in the solute atom group
  (defined by *select_indexN*) are considered as "solute"
  (scaled according to solute temperature) in
  gREST. Other terms are considered as "solvent" (kept intact).
  Solute-solvent terms are automatically determined from the solute
  selection. You can specify multiple terms (see examples).
  The parameter names are case-insensitive as follows:

  * **ALL**: all the available energy terms.
  * **BOND**: (aliases: **B**, **BONDS**): 1-2 bonding terms.
  * **ANGLE**: (aliases: **A**, **ANGLES**): 1-2-3 angle terms.
  * **UREY**: (aliases: **U**, **UREYS**): Urey-Bradley terms.
  * **DIHEDRAL**: (aliases: **D**, **DIHEDRALS**): 1-2-3-4 dihedral terms.
  * **IMPROPER**: (aliases: **I**, **IMPROPERS**): improper torsion terms.
  * **CMAP**: (aliases: **CM**, **CMAPS**): CMAP terms.
  * **CHARGE**: (aliases: **C**, **CHARGES**): coulombic interaction terms.
  * **LJ**: (aliases: **L**, **LJS**): Lennard-Jones interaction terms.

**analysis_grest** *Yes /No*

  **Default: No**
  
  Logical flag to write energy output for free energy calculations with REST. 
  If it is assigned to be *Yes*, the energy output file for each replica should be written in **[OUTPUT]** section.

.. note::
   Note that restraint energy terms defined in **[RESTRAINTS]** 
   cannot be treated as solute terms. 
   They never be affected by gREST solute temperatures.
   In SPDYN, water atoms cannot be specified as "solute" now.
   This limitation will be removed in the future version.

.. note::
   When the coulombic interaction terms are considered as the solute, 
   the solute region should have a net charge of 0 for an adequate 
   PME calculation.


Examples
========

Basically, REMD simulations in **GENESIS** can be carried out by 
just adding the **[REMD]** section 
in the control file of a normal MD simulation.
For details, see the online Tutorial (https://www.r-ccs.riken.jp/labs/cbrt/).

T-REMD
------

If the users want to carry out T-REMD simulations with 4 replicas in the NVT ensemble,
where each replica has the temperature 298.15, 311.79, 321.18, or 330.82 K,
and replica exchange is attempted every 1000 steps,
the following section is added to the control file of a normal MD simulation in the NVT ensemble:
:: 
  [REMD]
  dimension       = 1
  exchange_period = 1000
  type1           = TEMPERATURE
  nreplica1       = 4
  parameters1     = 298.15 311.79 321.18 330.82

As for the T-REMD simulation in the NPT ensemble, the users add this section
to the control file of a normal MD simulation in the NPT ensemble.
The REMD temperature generator (http://folding.bmc.uu.se/remd/) is
a useful tool to set the target temperature of each replica.


Two-dimensional REMD (T-REMD/REUS)
----------------------------------

The following is an example of two-dimensinal REMD, where
temperature and restraint are exchanged alternatively,
The 1st dimension is T-REMD with 8 parameters, 
and 2nd dimension is REUS in distance space with 4 parameters.
In total, 8 x 4 = 32 replicas are used:
:: 
  [REMD]
  dimension       = 2
  exchange_period = 1000
  type1           = TEMPERATURE
  nreplica1       = 8
  parameters1     = 298.15 311.79 321.18 330.82 340.70 350.83 361.23 371.89
  type2           = RESTRAINT
  nreplica2       = 4
  rest_function2  = 1

  [SELECTION]
  group1          = ai:25
  group2          = ai:392

  [RESTRAINTS]
  nfunctions      = 1
  function1       = DIST
  constant1       =  2.0   2.0   2.0   2.0
  reference1      = 10.0  10.5  11.0  11.5
  select_index1   = 1 2

These sections are added to the control file of a normal MD simulation.


Off-grid REUS
-------------

Example of off-grid REUS (merge two restraints into single reaction
coordinate), where distance and dihedral restraints are merged into
single reaction coordinate. First values of restraints ((2.0,10.0) for
distance, (10,-40) for dihedral) will be used for the first replica,
the fourth parameters ((2.0,11.5) for distance, (10,-10) for dihedral)
will be used for the fourth replica:
:: 
  [REMD]
  dimension       = 1
  exchange_period = 1000
  type1           = RESTRAINT           # REUS
  nreplica1       = 4
  rest_function1  = 1 2                 # off-grid REUS

  [SELECTION]
  group1          = ai:25
  group2          = ai:392
  group3          = ai:72
  group4          = ai:73
  group5          = ai:74
  group6          = ai:75

  [RESTRAINTS]
  nfunctions      = 2

  function1       = DIST
  constant1       =  2.0  2.0  2.0  2.0  # num of values must be nreplica1
  reference1      = 10.0 10.5 11.0 11.5
  select_index1   = 1 2

  function2       = DIHED
  constant2       =   10   10   10   10  # num of values must be nreplica1
  reference2      =  -40  -30  -20  -10
  select_index2   = 3 4 5 6


gREST
-----

In this example, the dihedral, CMAP, and LJ energy terms
in the selected atom groups are treated as "solute".
:: 
  [REMD]
  dimension       = 1
  exchange_period = 1000
  type1           = REST
  nreplica1       = 4
  parameters1     = 300.0 310.0 320.0 330.0  # solute temperatures
  param_type1     = D CM L                   # dihedral, CMAP, and LJ
  select_index1   = 1

  [SELECTION]
  group1          = ai:1-313

T-REMD in the two-dimensional REMD (T-REMD/REUS) may 
be replaced with gREST (gREST/REUS :cite:`Re:2019`) to reduce the required number of replicas. 
:: 
  [REMD]
  dimension       = 2
  exchange_period = 1000
  type1           = REST
  nreplica1       = 4
  parameters1     = 300.0 310.0 320.0 330.0  # solute temperatures
  param_type1     = D CM L                   # dihedral, CMAP, and LJ
  select_index1   = 3
  type2           = RESTRAINT
  nreplica2       = 4
  rest_function2  = 1

  [SELECTION]
  group1          = ai:25
  group2          = ai:392
  group3          = ai:1-313

  [RESTRAINTS]
  nfunctions      = 1
  function1       = DIST
  constant1       =  2.0   2.0   2.0   2.0
  reference1      = 10.0  10.5  11.0  11.5
  select_index1   = 1 2


