.. highlight:: bash
.. _dynamics:

=======================================================================
Dynamics section
=======================================================================

Molecular dynamics simulations
=======================================================================

In MD simulations, Newton's equation of motion (*F* = *ma*) is integrated numerically,
where the force *F* is derived from the first derivative of the potential energy function
with respect to the atomic position.
To date, various integrators have been proposed.
In the leap-frog algorithm, velocities are updated with

  .. raw:: latex
          
     \vspace{-5mm}

  .. math::  

     {{\bf{v}}_i}(t + \frac{{\Delta t}}{2}) = {{\bf{v}}_i}(t - \frac{{\Delta t}}{2}) + \frac{{\Delta t}}{{{m_i}}}{{\bf{F}}_i}(t),

  .. raw:: latex
          
     \vspace{-3mm}

and coordinates are updated with

  .. raw:: latex

     \vspace{-5mm}

  .. math::  

     {{\bf{r}}_i}(t + \Delta t) = {{\bf{r}}_i}(t) + \Delta t{{\bf{v}}_i}(t + \frac{{\Delta t}}{2}).

  .. raw:: latex
               
     \vspace{-3mm}

In the velocity Verlet algorithm, coordinates and velocities are obtained at the same time.
The velocities are updated with

  .. raw:: latex

     \vspace{-5mm}

  .. math::  

     {{\bf{v}}_i}\left( t + \Delta t \right) = {{\bf{v}}_i}\left( t \right) + \frac{{\Delta {t}}}{{{2m_i}}} 
     \left( {{\bf{F}}_i}\left( t \right) + {{\bf{F}}_i}\left( t + {\Delta t} \right) \right).

  .. raw:: latex

     \vspace{-3mm}

and the coordinates are updated with

  .. raw:: latex

     \vspace{-5mm}

  .. math::  

     {{\bf{r}}_i}\left( t + \Delta t \right) = {{\bf{r}}_i}\left( t \right) + \Delta t{{\bf{v}}_i}\left( t \right) + \frac{{\Delta {t^2}}}{{{2m_i}}}{{\bf{F}}_i}\left( t \right).

  .. raw:: latex

     \vspace{-3mm}

These velocity/coordinate updates are performed consecutively using the following procedure:

  .. raw:: latex

     \vspace{-5mm}

  .. math::

    {{\bf{v}}_i}\left( t + \frac{\Delta t}{2} \right) &= {{\bf{v}}_i}(t) &+& \frac{\Delta t}{2{m_i}}{{\bf{F}}_i}(t)

    {{\bf{r}}_i}\left( t + \Delta t \right) &= {{\bf{r}}_i}\left( t \right) &+& \Delta t{{\bf{v}}_i}\left( t + \frac{\Delta t}{2} \right)

    {{\bf{v}}_i}\left( t + \Delta t \right) &= {{\bf{v}}_i}\left( t + \frac{\Delta t}{2} \right) &+& \frac{\Delta t}{2{m_i}} {{\bf{F}}_i}\left( t + \Delta t \right)

  .. raw:: latex

     \vspace{-3mm}

In the case of multiple time step integration with r-RESPA, the force is split into slow and fast motion forces.
Slow motion force is less frequently evaluated to increase the performance while keeping the accuracy.

In **ATDYN**, both leap-frog and velocity Verlet integrators are available. 
In **SPDYN**, velocity Verlet and multiple time step integrator with veclotiy Verlet type are available.
The users must pay attention to the **[ENSEMBLE]** section as well,
because the algorithms that control the temperature and pressure are involved in the integrator.
For details, see :ref:`ensemble`.

-----------------------------------------------------------------------

**integrator** *LEAP / VVER / VRES / VVER_CG*

  **Default : VVER**
  
  * **LEAP**: leap-frog integrator (**ATDYN** only).
    
  * **VVER**: velocity Verlet integrator.
    
  * **VRES**: RESPA integrator (**SPDYN** only).

  * **VVER_CG**: velocity Verlet integrator for the coarse-grained simulations (**ATDYN** and Langevin thermostat only).


**timestep** *Real*

  **Default : 0.001** (unit : ps)

  Time step in the MD run.
  In general, timestep can be extended to 2 fs or longer,
  when the SHAKE, RATTLE, or SETTLE algorithms are employed.
  (see :ref:`constraints`).

**nsteps** *Integer*

  **Default : 100**

  Total number of steps in one MD run.
  If "timestep=0.001" and "nsteps=1000000" are specified,
  the users can carry out 1-ns MD simulation.

**eneout_period** *Integer*

  **Default : 10**

  Output frequency for the energy data. The trajectories are written
  in the log file every **eneout_period** steps during the simulation.
  For example, if "timestep=0.001" and "eneout_period=1000" are specified, 
  the energy is written every 1 ps.

**crdout_period** *Integer*

  **Default : 0**

  Output frequency for the coordinates data. The trajectories are written
  in the "dcdfile" specified in the **[OUTPUT]** section
  every **crdout_period** steps during the simulation.

**velout_period** *Integer*

  **Default : 0**

  Output frequency for the velocities data. The trajectories are written
  in the "dcdvelfile" specified in the **[OUTPUT]** section
  every **velout_period** steps during the simulation.

**rstout_period** *Integer*

  **Default : 0**

  Output frequency for the restart file. The restart information is written
  in the "rstfile" specified in the **[OUTPUT]** section
  every **rstout_period** steps during the simulation.

.. note::
  In the REMD or RPATH simulations, the value of **rstout_period** must be a multiple of *exchange_period* (REMD) or *rpath_period* (RPATH). 

**stoptr_period** *Integer* 

  **Default : 10**

  Frequency of removing translational and rotational motions of the whole system.
  Note that the rotational motion is not removed when the periodic boundary condition is employed.
  When you use positional restraints or RMSD restraints in the simulation,
  you may have to take care about removal of those motions.
  In some cases, such restraints can generate translational or rotational momentum in the system.
  If the momentum is frequently removed, the dynamics can be significantly disturbed.

**nbupdate_period** *Integer* 

  **Default : 10**

  Update frequency of the non-bonded pairlist.

**elec_long_period** *Integer* (**VRES** in **SPDYN** only)

  **Default : 1**

  Frequency of long-range interaction calculation.

**thermostat_period** *Integer* (**VRES** in **SPDYN** only)

  **Default : 1**

  Frequency of thermostat integration. It must be multiple of **elec_long_period**.

**barostat_period** *Integer* (**VRES** in **SPDYN** only)

  **Default : 1**

  Frequency of barostat integration. It must be multiple of **thermostat_period**.

**initial_time** *Real*

  **Default : 0.0** (unit : ps)

  Initial time of the MD run. Basically, you do not need to specify a certain value.
  This option is useful in the case of the restart MD run, because the initial time is reset to 0 ps.

**iseed** *Integer*

  **Default : automatically generated according to the current date and time**

  Seed for the pseudo-random number generator.
  This random number seed is used in the Langevin and Bussi thermostats (see `ensemble`).
  If **iseed** is not specified in the control file,
  it is automatically generated according to the current date and time.
  In the restart MD run, the random number seed is taken over from rstfile.
  However, if **iseed** value is specified in the control file in the restart run,
  it is alternatively used, and the seed in rstfile is neglected.

**verbose** *YES / NO*

  **Default : NO**

  Turn on or off the verbose output of the log information.
  For example, if "verbose=YES" is specified,
  virial and pressure of the system are written in the log file
  even in the case of the NVE or NVT ensemble.


Hydrogen mass repartitioning (HMR)
==================================

In **GENESIS2.0.0**, we can increase the time step in NVT/NPT conditions by evaluating temperature
and pressure in more accurate ways than conventional schemes. If we want to increase the time step, however,
we should be careful because there could be constraint errors. To avoid the error, we recommend the user to
make use of the HMR scheme with a large time step. In HMR, the mass of hydrogen atoms is increased by two to four
fold whereas the bonded heavy becomes lighter to conservere the total mass :cite:`Jung:2021_2`.

---------------------------------

**hydrogen_mr** *YES / NO*

  **Default : NO**

  Turn on or off the usage of HMR

**hmr_ratio** *Real*

  **Default : 3.0**

  Mass scaling factor of hydrogen atoms

**hmr_ratio_xh1** *Real*

  **Default : 3.0**

  Mass scaling facto of hydrogen atoms in XH1 group. If it is not written, scaling factor is decided by **hmr_ratio**.

**hmr_target** *All / Solute*

  **Default : All** (only when **hydrogen_mr** = *YES*)

  Target of HMR application. If **hmr_target** is *Solute*, HMR is not applied to the water molecules.


Simulated annealing and heating
=======================================================================

In simulated annealing or heating protocol, the following keywords
are additionally specified in the conventional MD simulation.
In the protocol used in **GENESIS**, the target temperature is changed linearly.
Note that the protocol is available only in the *LEAP* integrator.

-----------------------------------------------------------------------

**annealing** *YES / NO*

  **Default : NO**

  Turn on or off the simulated annealing or heating protocol.

**anneal_period** *Integer*

  **Default : 0**

  The target temperature is changed every **anneal_period** steps during the simulation.

**dtemperature** *Real*

  **Default : 0.0** (unit : Kelvin)

  Magnitude of changes of the target temperature.
  If **dtemperature** > 0, the temperature is increased by **dtemperature** every **anneal_period** steps.
  If **dtemperature** < 0, the temperature is decreased.


Targeted MD and Steered MD simulations
=======================================================================

In GENESIS, targeted MD (TMD) and steered MD (SMD) methods are available.
These methods are useful to guide a protein structure towards a target.
In SMD, restraint forces (or steering forces) are applied on the selected atoms,
where the RMSD with respect to the target is changed during the MD simulation. 
The restraint force is calculated from the derivative of the RMSD restraint potential:

  .. raw:: latex

     \vspace{-5mm}

  .. math::  
     U = \frac{1}{2}k\left(RMSD(t)-RMSD_{0}(t)\right)^{2}

  .. raw:: latex

     \vspace{-3mm}


where :math:`RMSD(t)` is the instantaneous RMSD of the current coordinates
from the target coordinates, and :math:`RMSD_{0}` is the target RMSD value.
The target RMSD value is changed linearly from the initial to final RMSD values:

  .. raw:: latex

     \vspace{-5mm}

  .. math::  
     RMSD_0(t) = RMSD_{\text{initial}}+\frac{t}{T}\left(RMSD_{\text{final}}-RMSD_{\text{initial}}\right)

  .. raw:: latex

     \vspace{-3mm}

where :math:`T` is the total MD simulation time. Targeted MD (TMD), 
originally suggested by J. Schlitter et al. :cite:`Schlitter:1994`,
is different from SMD in that force constants are changed during MD
simulations. If the users perform SMD, there is a possibility observing the
large difference between the instantaneous RMSD and target RMSD.
In TMD, force constants are given by 
Lagrangian multipliers to overcome the energy barrier between the
instantaneous and target RMSDs. Therefore, the users could find 
trajectories where RMSD is almost identical to the target RMSD at 
each time. In **[SELECTION]** section, the users select atoms involved
in RMSD calculations for SMD or TMD. Users should specify either
*RMSD* or *RMSDMASS* (mass-weighted RMSD) in **[RESTRAINTS]** section to 
run TMD or SMD. In SMD, force constants defined in **[RESTRAINTS]** 
section are used, but force constants are automatically determined using
Lagrangian multipliers during simulation in TMD.

-----------------------------------------------------------------------

**target_md** *YES / NO* 

  **Default : NO**

  Turn on or off the targeted MD simulation. 

**steered_md** *YES / NO* 

  **Default : NO**

  Turn on or off the steered MD simulation.

**initial_rmsd** *Real* 

  **Default : 0.0** (unit : :math:`\text{\AA}`)

  Initial value of the reference rmsd.
  If not specified explicitly, it is calculated from the initial and reference structures.

**final_rmsd** *Real* 

  **Default : 0.0** (unit : :math:`\text{\AA}`)

  Final value of the reference rmsd.

.. note::
   In the RMSD restraint, structure fitting scheme is specified in the **[FITTING]**
   section (see :ref:`fitting`). Since the default behavior was significantly
   changed in ver. 1.1.5 (no fitting applied on the default setting),
   the users of 1.1.4 or before must pay special attention on the fitting scheme.
   In versions of 1.1.4 or before, structure fitting is automatically
   applied for the atoms concerning restraint potential.


Examples
========

100-ps MD simulation with the velocity Verlet integrator with the timestep of 2 fs:
:: 
  [DYNAMICS]
  integrator        =   VVER  # velocity Verlet
  nsteps            =  50000  # number of MD steps (100ps)
  timestep          =  0.002  # timestep (2fs)
  eneout_period     =    500  # energy output period (1ps)
  crdout_period     =    500  # coordinates output period (1ps)
  rstout_period     =  50000  # restart output period
  nbupdate_period   =     10  # nonbond pair list update period

100-ps MD simulation with the RESPA integrator with the timestep of 2.5 fs:
:: 
  [DYNAMICS]
  integrator        =   VRES  # RESPA integrator
  nsteps            =  40000  # number of MD steps (100ps)
  timestep          = 0.0025  # timestep (2.5fs)
  eneout_period     =    400  # energy output period (1ps)
  crdout_period     =    400  # coordinates output period (1ps)
  rstout_period     =  40000  # restart output period
  nbupdate_period   =     10  # nonbond pair list update period
  elec_long_period  =      2  # period of reciprocal space calculation
  thermostat_period =     10  # period of thermostat update
  barostat_period   =     10  # period of barostat update

The following is an example for simulated annealing in the NVT ensemble (see :ref:`Ensemble`),
where the temperature is decreased from 500 K by 2 K
every 250 steps in the 250,000-steps MD simulation (1 step = 2 fs).
Thus, the temperature eventually reaches to 300 K during 50 ps.
:: 
  [DYNAMICS]
  integrator       =   VVER   # leap-frog integrator
  nsteps           =  25000   # number of MD steps
  timestep         =  0.002   # timestep (ps)
  nbupdate_period  =     10   # nonbond pair list update period
  annealing        =    YES   # simulated annealing
  dtemperature     =   -2.0   # delta temperature
  anneal_period    =    250   # temperature change period

  [ENSEMBLE]
  ensemble         = NVT      # [NVT,NPT,NPAT,NPgT]
  tpcontrol        = BUSSI    # [BERENDSEN,NHC,BUSSI]
  temperature      = 500.0    # initial temperature (K)
