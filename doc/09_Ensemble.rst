.. highlight:: bash
.. _ensemble:

=======================================================================
Ensemble section
=======================================================================

Thermostat and barostat
=======================

In the **[ENSEMBLE]** section, the type of ensemble, temperature and
pressure control algorithm, and parameters used in these algorithms
(such as temperature and pressure) can be specified.

In the Berendsen thermostat with weak coupling ("tpcontrol=berendsen"), velocities are rescaled by

  .. raw:: latex

     \vspace{-5mm}

  .. math::
    \lambda = \left [ 1 + \frac{\Delta t}{\tau} \left ( \frac{K_{\textrm{target}}}{K} - 1 \right) \right]^{\frac{1}{2}}

  .. raw:: latex

     \vspace{-3mm}

at every step where :math:`{\tau}` is the dampling time :cite:`Berendsen:1984`. It almost conserves dynamics
of NVE and one of the most popular thermostats. However, it fails to generate canonical distributions and
anomalous effects can happern :cite:`Braun:2018`. Bussi suggested a modification of Berendsen thermostat to
generate canonical distribution :cite:`Bussi:2007`. According to the scheme ("tpcontrol=bussi"),
            the target temperature is not a fixed value but stochastically driven from the distributions:

  .. raw:: latex

     \vspace{-5mm}

  .. math::

   P \left( {K_{\textrm{target}}} \right ) \propto {K_{\textrm{target}}}^{\frac{N_{f}}{2}-1} e ^ {-\beta K_{\textrm{target}}} .
  .. raw:: latex

     \vspace{-3mm}

In addition, random force :math:`{\Delta W}` is applied to update temperature based on the instantaneousand target temperatures:

  .. raw:: latex

     \vspace{-5mm}

  .. math::


   \Delta K = \left ( K_{\textrm{target}} - K \right ) \frac{\Delta t}{\tau} + 2 \sqrt{\frac{K K_{\textrm{target}}}{N_{f}}} \frac{\Delta W}{\tau} .
  .. raw:: latex

     \vspace{-3mm}

Another famous thermostat is suggested by Martyna et al., which is named as Nose-Hoover chain (NHC) :cite:`Martyna:1992`.
This thermostat was designed to solve non-ergodic progem in Nose-Hoover thermostat.
In this scheme, we solve the euqation of motion with additional \textit{m} chains:

  .. raw:: latex

     \vspace{-5mm}

  .. math::

   \frac{d\mathbf{r}_{k}}{dt} & = \frac{\mathbf{p}_{k}}{m_{k}}

   \frac{d\mathbf{p}_{k}}{dt} & = \mathbf{F}_{k} - \frac{p_{\eta_{1}}}{Q_{1}} \mathbf{p}_{k}

   \frac{d \eta_{\alpha}}{dt} & = \frac{p_{\eta_{\alpha}}}{Q_{\alpha}}, \  \alpha=1,2,...,m

   \frac{d p_{\eta_{1}}}{dt} & = K - K_{\textrm{target}} - \frac{p_{\eta_{2}}}{Q_{2}} p_{\eta_{1}}

   \frac{d p_{\eta_{\alpha}}}{dt} & = \frac{p_{\eta_{\alpha-1}}}{Q_{\alpha-1}} - k_{B}T - \frac{p_{\eta_{\alpha + 1}}}{Q_{\alpha + 1}} p_{\eta_{\alpha}}, \ \alpha = 2,...,m-1

   \frac{d p_{\eta_{m}}}{dt} & = \frac{p_{\eta_{m-1}}}{Q_{m-1}} - k_{B}T .

  .. raw:: latex

     \vspace{-3mm}


In **SPDYN**, we evaluate temperature in more accurate way than existing ones. With velocity Verlet
integration, two kinetic energy types can be used in temperature evaluations:

  .. raw:: latex

     \vspace{-5mm}

  .. math::

   K_{\textrm{full}} & = \sum_{k=1}^{N} \frac{\mathbf{p}_{k}\left( t \right) ^2}{2m_{k}}

   K_{\textrm{half}} & = \frac{1}{2} \sum_{k=1}^{N} \left ( \frac{\mathbf{p}_{k} \left(t - \frac{\Delta t}{2} \right)^2}{2m_{k}} + \frac{\mathbf{p}_{k} \left( t + \frac{\Delta t}{2} \right)^2}{2m_{k}} \right)

  .. raw:: latex

     \vspace{-3mm}

These kinetic energies are accurate up to the first order of a time step and not recommended for a large time step.
In **SPDYN**, temperature is evaluted by combining the two kinetic energies:

  .. raw:: latex

     \vspace{-5mm}

  .. math::

   N_{f} k_{B} T = \frac{4}{3} \left < K_{\textrm{half}} (t) \right > + \frac{2}{3} \left < K_{\textrm{full}} (t) \right >

  .. raw:: latex

     \vspace{-3mm}

Temperature evaluted in this way is accurate up to the third order of the time step.

Barostat in **SPDYN** is based on MTK barostat type suggested by :cite:`Martyna:1996`
with three thermostats described above. To make use of kinetic energy at :math:`t + \frac{\Delta t}{2}`,
we recommend group temperature/pressure evaluations where kinetic energy and virial are evaluted by
considering XHn group one particle.

In the Langevin thermostat algorithm ("ensemble=NVT" with
"tpcontrol=LANGEVIN"), every particles are coupled with a viscous background
and a stochastic heat bath :cite:`Adelman:1976`:

  .. raw:: latex
     
     \vspace{-5mm}

  .. math:: 
    \frac{d\mathbf{v}(t)}{dt} = \frac{\mathbf{F}(t)+\mathbf{R}(t)}{m}-\gamma \mathbf{v}(t)

  .. raw:: latex
     
     \vspace{-3mm}

where :math:`\gamma` is the thermostat friction parameter (*gamma_t* keyword) 
and :math:`\mathbf{R}(t)` is the stochastic force. 
In the Langevin thermostat and barostat method ("ensemble=NPT" with
"tpcontrol=LANGEVIN"),
the equation of motion is given by :cite:`Quigley:2004`:

  .. raw:: latex
     
     \vspace{-5mm}

  .. math::
    \frac{d\mathbf{r}(t)}{dt} & = \mathbf{v}(t)+v_{\epsilon}\mathbf{r}(t) \\ \frac{d\mathbf{v}(t)}{dt} & = \frac{\mathbf{F}(t)+\mathbf{R}(t)}{m}-[\gamma_p+(1+\frac{3}{f})v_{\epsilon}]\mathbf{v}(t) \\ \frac{dv_{\epsilon}(t)}{dt} & = [3V(P(t)-P_0(t))+\frac{3K}{f}-\gamma_p v_{\epsilon} + R_p ] / p_{mass}

  .. raw:: latex
     
     \vspace{-3mm}

where :math:`K` is the kinetic energy, :math:`\gamma_p` is the barostat 
friction parameter (*gamma_p* keyword), :math:`R_p` is the stochastic 
pressure variable. 

-----------------------------------------------------------------------

**ensemble** *NVE / NVT / NPT / NPAT / NPgT*

  **Default : NVE**

  Type of ensemble.

  * **NVE**: Microcanonical ensemble.

  * **NVT**: Canonical ensemble.

  * **NPT**: Isothermal-isobaric ensemble.

  * **NPAT**: Constant area A (XY), pressure along the normal (Z), temperature
    :cite:`Zhang:1995`. In this case, *isotropy* must be set to 'XY-FIXED' (see below).

  * **NPgT**: Constant surface-tension :math:`\gamma` (XY), pressure along the normal (Z),
    temperature :cite:`Zhang:1995`. In this case, *isotropy* must be set to 'SEMI-ISO' (see below).

**temperature** *Real*

  **Default : 298.15** (unit : Kelvin)

  Initial and target temperature.

**pressure** *Real*

  **Default : 1.0** (unit : atm)

  Target pressure in the NPT ensemble.
  In the case of the NPAT and NPgT ensembles,
  this is the pressure along the 'Z' axis.

**group_tp** *Yes / No*

  **Default : Yes**

    Use of group temperature/pressure evaluation in thermostat/barostat. :cite:`Jung:2020`

**gamma** *Real*

  **Default : 0.0** (unit : dyn/cm)

  Target surface tension in NPgT ensemble.

**tpcontrol** *NO / BERENDSEN / BUSSI / NHC / LANGEVIN*

  **Default : NO**

  Type of thermostat and barostat. The availabe algorithm depends on the integrator.

  * **NO**: Do not use temperature/pressure control algorithm (for NVE only)

  * **BERENDSEN**: Berendsen thermostat/barostat :cite:`Berendsen:1984`

  * **BUSSI**: Bussi's thermostat/barostat :cite:`Bussi:2007` :cite:`Bussi:2009`

  * **NHC**: Nose-Hoover chain thermostat with MTK barostat :cite:`Berendsen:1984` :cite:`Martyna:1996`

  * **LANGEVIN**: Langevin thermostat/barostat :cite:`Quigley:2004` (for ATDYN only)

  +------------+-------------+------------------------------+
  | integrator |  ensemble   |  tpcontrol                   |
  +============+=============+==============================+
  | LEAP       |  NVT        |  BERENDSEN, LANGEVIN         |
  |            +-------------+------------------------------+
  | (ATDYN)    |  NPT        |  BERENDSEN, LANGEVIN         |
  |            +-------------+------------------------------+
  |            |  NPAT/NPgT  |  BERENDSEN, LANGEVIN         |
  +------------+-------------+------------------------------+
  | VVER       |  NVT        |  BERENDSEN, LANGEVIN, BUSSI  |
  |            +-------------+------------------------------+
  | (ATDYN)    |  NPT        |  LANGEVIN, BUSSI             |
  |            +-------------+------------------------------+
  |            |  NPAT/NPgT  |  LANGEVIN                    |
  +------------+-------------+------------------------------+
  | VVER       |  NVT        |  BUSSI, BERENDSEN, NHC       |
  |            +-------------+------------------------------+
  | (SPDYN)    |  NPT        |  BUSSI, BERENDSEN, NHC       |
  |            +-------------+------------------------------+
  |            |  NPAT/NPgT  |  BUSSI, BERENDSEN, NHC       |
  +------------+-------------+------------------------------+
  | VRES       |  NVT        |  BUSSI, BERENDSEN, NHC       |
  |            +-------------+------------------------------+
  | (SPDYN)    |  NPT        |  BUSSI, BERENDSEN, NHC       |
  |            +-------------+------------------------------+
  |            |  NPAT/NPgT  |  BUSSI, BERENDSEN, NHC       |
  +------------+-------------+------------------------------+

**tau_t** *Real*

  **Default : 5.0** (unit : ps)

  Temperature coupling time in the Berendsen and Bussi thermostats.

**tau_p** *Real*

  **Default : 5.0** (unit : ps)

  Pressure coupling time in the Berendsen and Bussi barostats.

**compressibility** *Real*

  **Default : 0.0000463** (unit : atm\ :sup:`-1`) 

  Compressibility parameter in the Berendsen barostat.

**gamma_t** *Real*

  **Default : 1.0** (unit : ps\ :sup:`-1`) 

  Friction parameter of the Langevin thermostat.

**gamma_p** *Real*

  **Default : 0.1** (unit : ps\ :sup:`-1`)

  Friction parameter of the Langevin barostat.

**isotropy** *ISO / ANISO / SEMI-ISO / XY-FIXED*

  **Default : ISO**

  Isotropy of the simulation system.
  This parameter specifies how X, Y, Z dimensions of the simulation 
  box change in NPT, NPgT, and NPAT ensembles.

  * **ISO**: X, Y, and Z dimensions are coupled together.

  * **ANISO**: X, Y, and Z dimensions fluctuate independently.

  * **SEMI-ISO**: X, Y, and Z dimensions fluctuate, where the ratio of
    X and Y dimensions are kept constant, and Z dimension can change
    independently :cite:`Kandt:2007`. This setting with NPT or NPAT or NPgT
    ensemble is expected to be useful for bio-membrane systems.

  * **XY-FIXED**: X and Y dimensions are fixed, while Z dimension can change
    (NPAT only).

Examples
========

NVT ensemble with Bussi thermostat:
:: 
  [ENSEMBLE] 
  ensemble    = NVT       # Canonical ensemble 
  tpcontrol   = BUSSI     # Bussi thermostat 
  temperature = 300.0     # target temperature (K)

NPT ensemble with isotropic pressure coupling:
:: 
  [ENSEMBLE] 
  ensemble    = NPT       # Isothermal-isobaric ensemble
  tpcontrol   = BUSSI     # Bussi thermostat and barostat
  temperature = 300.0     # target temperature (K)
  pressure    = 1.0       # target pressure (atm)

NPT ensemble with semi-isotropic pressure coupling, which is usually used for lipid bilayer systems:
:: 
  [ENSEMBLE] 
  ensemble    = NPT       # Isothermal-isobaric ensemble
  tpcontrol   = BUSSI     # Bussi thermostat and barostat
  temperature = 300.0     # target temperature (K)
  pressure    = 1.0       # target pressure (atm)
  isotropy    = SEMI-ISO  # Ratio of X to Y is kept constant

NPAT ensemble:
:: 
  [ENSEMBLE] 
  ensemble    = NPAT      # Constant area ensemble
  tpcontrol   = BUSSI     # Bussi thermostat and barostat
  temperature = 300.0     # target temperature (K)
  pressure    = 1.0       # target normal pressure (atm)
  isotropy    = XY-FIXED  # the system area is kept constant

NP\ :math:`\gamma`\T ensemble:
:: 
  [ENSEMBLE] 
  ensemble    = NPgT      # Constant surface-tension ensemble
  tpcontrol   = BUSSI     # Bussi thermostat and barostat
  temperature = 300.0     # target temperature (K)
  pressure    = 1.0       # target normal pressure (atm)
  gamma       = 200.0     # target surface tension (dyn/cm)
  isotropy    = SEMI-ISO  # Ratio of X to Y is kept constant

