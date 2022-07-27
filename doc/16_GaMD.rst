.. highlight:: bash
.. _gamd:

=======================================================================
GAMD section
=======================================================================

Gaussian accelerated Molecular Dynamics
=======================================

In the **[GAMD]** section, the users can specify keywords for Gaussian 
accelerated Molecular Dynamics (GaMD) simulation.
The GaMD method
:cite:`Miao:2015,Pang:2017`
accelerates the conformational sampling of biomolecules by adding
a harmonic boost potential to smooth their potential energy surface.
GaMD has the advantage that reaction coordinates do not need to be predefined,
thus setting up the system for the simulation is rather easy.
The use of the harmonic boost potential allows to recover the
unbiased free-energy changes through cumulant expansion to the second order,
which resolves the practical reweighting problem in the original accelerated MD method.

GaMD was developed as a potential-biasing method for enhanced sampling.
It accelerates the conformational sampling of a biomolecule by adding a 
non-negative boost potential to the system potential energy :math:`U(\vec{x})`: 

  .. raw:: latex

     \vspace{-5mm}

  .. math::  
     U'(\vec{x}) = U(\vec{x}) + \Delta U^{\mathrm{GaMD}}\left(U(\vec{x})\right),


  .. raw:: latex

     \vspace{-3mm}

where :math:`\vec{x}` is the configuration of the system,
:math:`U'(\vec{x})` is the modified potential energy,
and :math:`\Delta U^{\mathrm{GaMD}}` is the boost potential 
depending only on :math:`U(\vec{x})`.

In conventional accelerated MD :cite:`Hamelberg:2004,Hamelberg:2007,Shen:2008`,
the average of the Boltzmann factors of the boost potential terms appears in the 
reweighting equation of the probability along the selected reaction coordinates,
causing a large statistical error.
In order to reduce the noise, GaMD uses a harmonic boost potential, which adopts
a positive value only when the system potential is lower than an energy threshold :math:`E`:

  .. raw:: latex

     \vspace{-5mm}

  .. math::  
     \Delta U^{\mathrm{GaMD}}\left(U(\vec{x})\right) = 
     \begin{cases}
       \frac{1}{2}k\{E-U(\vec{x})\}^2 & (U(\vec{x}) < E) \\
       0 & (U(\vec{x}) \ge E)
     \end{cases},

  .. raw:: latex

     \vspace{-3mm}

where :math:`k` is a harmonic force constant.
:math:`U'(\vec{x})` should satisfy the following relationships :cite:`Miao:2015,Pang:2017`:
:math:`U'(\vec{x}_1) < U'(\vec{x}_2)` and :math:`U'(\vec{x}_2) - U'(\vec{x}_1) < U(\vec{x}_2) - U(\vec{x}_1)`
if :math:`U(\vec{x}_1) < U(\vec{x}_2)`.
To keep the relationships, the threshold energy needs to be set as:

  .. raw:: latex

     \vspace{-5mm}

  .. math::  
     U_{\mathrm{max}} \le E \le U_{\mathrm{min}} + \frac{1}{k},

  .. raw:: latex

     \vspace{-3mm}

where :math:`U_{\mathrm{max}}` and :math:`U_{\mathrm{min}}` are
maximum and minimum energies of the system, respectively.
To ensure accurate reweighting, the deviation of the potential
must also satisfy the relation:

  .. raw:: latex

     \vspace{-5mm}

  .. math::  
     k(E-U_{\mathrm{ave}}) \sigma_{U} \le \sigma_0,

  .. raw:: latex

     \vspace{-3mm}

where :math:`U_{\mathrm{ave}}` and :math:`\sigma_{U}` are the average
and standard deviation of :math:`U(\vec{x})`, respectively.
:math:`\sigma_0` is a user-specified upper limit.
:math:`k_0` is defined as :math:`k_0 \equiv k(U_{\mathrm{max}}-U_{\mathrm{min}})`,
then :math:`0 < k_0 \le 1`.

When :math:`E` is set to the lower bound :math:`U_{\mathrm{max}}`,
:math:`k_0` is determined by

  .. raw:: latex

     \vspace{-5mm}

  .. math::  
     k_0 = \min \left(1,\frac{\sigma_0}{\sigma_U}\frac{U_{\mathrm{max}}-U_{\mathrm{min}}}{U_{\mathrm{max}}-U_{\mathrm{ave}}}\right)

  .. raw:: latex

     \vspace{-3mm}

When :math:`E` is set to the upper bound :math:`U_{\mathrm{min}}+1/k`,
:math:`k_0` is set to

  .. raw:: latex

     \vspace{-5mm}

  .. math::  
     k''_0 \equiv \left(1-\frac{\sigma_0}{\sigma_U}\right)\frac{U_{\mathrm{max}}-U_{\mathrm{min}}}{U_{\mathrm{ave}}-U_{\mathrm{min}}}

  .. raw:: latex

     \vspace{-3mm}

if :math:`0 < k''_0 < 1`, and :math:`k_0` is set to 1 otherwise.

The above parameters
(:math:`U_{\mathrm{max}}`, :math:`U_{\mathrm{min}}`, :math:`U_{\mathrm{ave}}`, and :math:`\sigma_{U}`)
are determined from short-time simulations a priori.
When the distribution of the boost potential approaches Gaussian distribution,
the cumulant expansion of the average of :math:`\exp[\beta \Delta U^\mathrm{GaMD}]`
to the second order provides a good approximation for the free energy :cite:`Miao:2014`.

GaMD can be combined with REUS in such a way that each replica in REUS is accelerated by the GaMD boost potential:

  .. raw:: latex

     \vspace{-5mm}

  .. math::  
     U''_i(\vec{x}) &= U'(\vec{x}) + \Delta U_i^{\mathrm{REUS}}\left(\xi(\vec{x})\right) \\
     &= U(\vec{x}) + \Delta U^{\mathrm{GaMD}}\left(U(\vec{x})\right) + \Delta U_i^{\mathrm{REUS}}\left(\xi(\vec{x})\right),

  .. raw:: latex

     \vspace{-3mm}

where :math:`U''_i(\vec{x})` is the modified potential energy of replica :math:`i`,
:math:`\Delta U_i^{\mathrm{REUS}}` is the bias potential of REUS for replica :math:`i`, and
:math:`\xi(\vec{x})` is the collective variable of REUS.
This method is referred to as Gaussian accelerated replica exchange umbrella sampling (GaREUS) :cite:`Oshima:2019`.
The parameters in the GaMD boost potential are used in all replicas of GaREUS simulations.
By using this combination, the simulated system in each replica becomes more flexible, or the energy
barrier irrelevant to the collective variable is lowered, enhancing the sampling efficiency.
When performing GaREUS simulations, the user must specify [REMD] section to use REUS
and define a collective variable in the [SELECTION] and [RESTRAINTS] sections.
Please check the example below.

-----------------------------------------------------------------------

**gamd** *YES / NO*

  **Default : NO**

  Enable the GaMD method.

**boost** *YES / NO*

  **Default : YES**

  Flag to apply GaMD boost to the system.
  If *boost = NO*, boost is not applied but GaMD parameters
  are updated from the trajectory.

**boost_type** *DUAL / DIHEDRAL / POTENTIAL*

  **Default: DUAL**

  Type of boost.

  * **DUAL**: Boost is applied on both the dihedral and total potential energies.
  * **DIHEDRAL**: Boost is applied on only the dihedral energy.
  * **POTENTIAL**: Boost is applied on only the total potential energy.

**thresh_type** *LOWER / HIGHER*

  **Default: LOWER**

  Type of threshold.

  * **LOWER**: :math:`E` is set to the lower bound :math:`E = U_{\text{max}}`.
  * **HIGHER**: :math:`E` is set to its upper bound :math:`E = U_{\text{min}} + 1/k`.

**update_period** *Integer*

  **Default: 0**

  Period of updating parameters in units of time step.
  When **update_period** = 0, GaMD parameters are not updated during the GaMD simulation.
  When **update_period** > 0, the GaMD simulation updates its parameters every **update_period** steps, and then the updated parameters
  are output to **gamdfile** specified in the **[OUTPUT]** section.
  The file includes the maximum, minimum, average, and deviation of the total potential or dihedral potential (**pot_max**, **pot_min**, **pot_ave**, **pot_dev**, **dih_max**, **dih_min**, **dih_ave**, **dih_dev**), which are calculated within the interval **update_period**.

**sigma0_pot** *Real*

  **Default: 6.0** (unit: kcal/mol)

  Upper limit of the standard deviation of the total potential
  boost (:math:`\sigma_0^{\mathrm{pot}}`) that allows for
  accurate reweighting.

**pot_max** *Real*

  **Default: -99999999.0** (unit: kcal/mol)

  Maximum of the total potential energy of the system, :math:`U_{\mathrm{max}}^{\mathrm{pot}}`.
  When **update_period** is not zero, :math:`U_{\mathrm{max}}^{\mathrm{pot}}` is updated every time step.
  If :math:`U^{\mathrm{pot}}` becomes larger than :math:`U_{\mathrm{max}}^{\mathrm{pot}}`, 
  :math:`U_{\mathrm{max}}^{\mathrm{pot}}` is set to the larger value.
  To determine the initial value of :math:`U_{\mathrm{max}}^{\mathrm{pot}}`, the default value is set to a large negative number.

**pot_min** *Real*

  **Default: 99999999.0** (unit: kcal/mol)

  Minimum of the total potential energy of the system, :math:`U_{\mathrm{min}}^{\mathrm{pot}}`.
  When **update_period** is not zero, :math:`U_{\mathrm{min}}^{\mathrm{pot}}` is updated every time step.
  If :math:`U^{\mathrm{pot}}` becomes smaller than :math:`U_{\mathrm{min}}^{\mathrm{pot}}`, 
  :math:`U_{\mathrm{min}}^{\mathrm{pot}}` is set to the smaller value.
  To determine the initial value of :math:`U_{\mathrm{min}}^{\mathrm{pot}}`, the default value is set to a large positive number.

**pot_ave** *Real*

  **Default: 0.0** (unit: kcal/mol)

  Average of the total potential energy of the system, :math:`U_{\mathrm{ave}}^{\mathrm{pot}}`.

**pot_dev** *Real*

  **Default: 0.0** (unit: kcal/mol)

  Standard deviation of the total potential energy of the system, :math:`\sigma_{U}^{\mathrm{pot}}`.

**sigma0_dih** *Real*

  **Default: 6.0** (unit: kcal/mol)

  Upper limit of the standard deviation of the dihedral
  boost (:math:`\sigma_0^{\mathrm{dih}}`) that allows 
  for accurate reweighting.

**dih_max** *Real*

  **Default: -99999999.0** (unit: kcal/mol)

  Maximum of the dihedral energy of the system, :math:`U_{\mathrm{max}}^{\mathrm{dih}}`.
  When **update_period** is not zero, :math:`U_{\mathrm{max}}^{\mathrm{dih}}` is updated every time step.
  If :math:`U^{\mathrm{dih}}` becomes larger than :math:`U_{\mathrm{max}}^{\mathrm{dih}}`, 
  :math:`U_{\mathrm{max}}^{\mathrm{dih}}` is set to the larger value.
  To determine the initial value of :math:`U_{\mathrm{max}}^{\mathrm{dih}}`, the default value is set to a large negative number.

**dih_min** *Real*

  **Default: 99999999.0** (unit: kcal/mol)

  Minimum of the dihedral energy of the system, :math:`U_{\mathrm{min}}^{\mathrm{dih}}`.
  When **update_period** is not zero, :math:`U_{\mathrm{min}}^{\mathrm{dih}}` is updated every time step.
  If :math:`U^{\mathrm{dih}}` becomes smaller than :math:`U_{\mathrm{min}}^{\mathrm{dih}}`, 
  :math:`U_{\mathrm{min}}^{\mathrm{dih}}` is set to the smaller value.
  To determine the initial value of :math:`U_{\mathrm{min}}^{\mathrm{dih}}`, the default value is set to a large positive number.

**dih_ave** *Real*

  **Default: 0.0** (unit: kcal/mol)

  Average of the dihedral energy of the system, :math:`U_{\mathrm{ave}}^{\mathrm{dih}}`.

**dih_dev** *Real*

  **Default: 0.0** (unit: kcal/mol)

  Standard deviation of the dihedral energy of the system, :math:`\sigma_{U}^{\mathrm{dih}}`.


Examples
========

Example of a GaMD simulation to determine initial parameters.
To obtain the initial guess of the boost potential, 
(pot_max, pot_min, pot_ave, pot_dev, dih_max, dih_min, dih_ave, dih_dev)
are calculated from a short simulation without boosting.
::
  
  [GAMD]
  gamd          = yes
  boost         = no
  boost_type    = DUAL
  thresh_type   = LOWER
  sigma0_pot    = 6.0
  sigma0_dih    = 6.0
  update_period = 50000

Example of a GaMD simulation updating parameters.
The boost potential is updated every *update_period* during the simulation.
::

  [GAMD]
  gamd          = yes
  boost         = yes
  boost_type    = DUAL
  thresh_type   = LOWER
  sigma0_pot    = 6.0
  sigma0_dih    = 6.0
  update_period = 500
  pot_max       = -20935.8104
  pot_min       = -21452.3778
  pot_ave       = -21183.9911
  pot_dev       = 78.1207
  dih_max       = 16.4039
  dih_min       = 8.5882
  dih_ave       = 11.0343
  dih_dev       = 1.0699

Example of a GaMD simulation for production.
In order to fix the parameters (pot_max, pot_min, pot_ave, pot_dev, 
dih_max, dih_min, dih_ave, dih_dev), *update_period* is set to 0.
::

  [GAMD]
  gamd          = yes
  boost         = yes
  boost_type    = DUAL
  thresh_type   = LOWER
  sigma0_pot    = 6.0
  sigma0_dih    = 6.0
  update_period = 0
  pot_max       = -20669.2404
  pot_min       = -21452.3778
  pot_ave       = -20861.5224
  pot_dev       = 48.9241
  dih_max       = 23.2783
  dih_min       = 8.5882
  dih_ave       = 13.3806
  dih_dev       = 1.7287

Example of a GaREUS simulation.
The same GaMD parameters are applied in each replica of REUS.
After the simulation, the two-step reweighting procedure using the
multistate Bennett acceptance ratio method and the cumulant expansion 
for the exponential average is required to obtain the unbiased free-energy landscapes.
::

  [REMD]
  dimension        = 1
  exchange_period  = 5000
  type1            = RESTRAINT
  nreplica1        = 4
  rest_function1   = 1

  [GAMD]
  gamd          = yes
  boost         = yes
  boost_type    = DUAL
  thresh_type   = LOWER
  sigma0_pot    = 6.0
  sigma0_dih    = 6.0
  update_period = 0
  pot_max       = -26491.7344
  pot_min       = -27447.4316
  pot_ave       = -26744.5742
  pot_dev       = 52.5674
  dih_max       = 135.8921
  dih_min       = 91.2309
  dih_ave       = 116.8572
  dih_dev       = 3.6465

  [SELECTION]
  group1 = rno:1  and an:CA
  group2 = rno:10 and an:CA

  [RESTRAINTS]
  nfunctions    = 1
  function1     = DISTMASS
  constant1     = 1.0 1.0 1.0 1.0
  reference1    = 5.0 6.0 7.0 8.0
  select_index1 = 1 2
