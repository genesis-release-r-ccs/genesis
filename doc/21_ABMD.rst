.. highlight:: bash
.. _abmd:

=======================================================================
ABMD section
=======================================================================

Adaptively Biassed Molecular Dynamics
=======================================

In the **[ABMD]** section, the users can specify keywords for Adaptively 
Biassed Molecular Dynamics (ABMD) simulation :cite:`Babin:2008,Babin:2009`.
ABMD is an enhanced sampling method that applies a time-dependent bias 
potential along predefined collective variables (CVs) to accelerate 
sampling of rare events and explore the free energy landscape.
The ABMD method requires the **[ReacCoord]** section to define the 
reaction coordinates. Please refer to the **[ReacCoord]** section for details 
on how to set up the reaction coordinates.

The time-dependent bias potential, :math:`U(\xi \mid t)`, at time *t* along 
predefined CV, :math:`\xi`, is added to the system potential energy :math:`U(\vec{x})`:

  .. raw:: latex

     \vspace{-5mm}

  .. math::  
     U'(\vec{x}) = \sum_{m=-1}^{n} U(\vec{x}) + U(\xi \mid t),


  .. raw:: latex

     \vspace{-3mm}

where :math:`\vec{x}` is the configuration of the system,
:math:`U'(\vec{x})` is the modified potential energy.
The :math:`U(\xi \mid t)` is formulated in the following equation:


  .. raw:: latex

     \vspace{-5mm}

  .. math::
     U(\xi \mid t) = \sum_{m=-1}^{n} U_m(t)\, B\left(\frac{\xi - m}{\Delta \xi}\right),


  .. raw:: latex

     \vspace{-3mm}

where :math:`B(\vec{x})` is a cubic B-spline, :math:`U_m(t)` is B-spline
coeffcients which are updated at every time step :math:`\Delta t`. The :math:`U_m(t)`
are updated using Euler-like discretization scheme with the biweight kernel, *G*, 
given by,


  .. raw:: latex

     \vspace{-5mm}

  .. math::
     U_m(t+\Delta t) = U_m(t) + \Delta t \, \frac{k_B T}{\tau_F} \, G\!\left(\frac{\sigma}{\Delta \xi} - m\right)


  .. raw:: latex

     \vspace{-3mm}

where :math:`\tau_F` is a flooding time scale, :math:`k_B` is a Boltzmann constant,
and :math:`\sigma` is a CV at time *t*. As :math:`t \rightarrow \infty`, the bias potential 
converges and allows the free-energy landscape :math:`F(\xi)` to be estimated as


  .. raw:: latex

     \vspace{-5mm}

  .. math::
     F(\xi) = -U(\xi) + \mathrm{Const}.


  .. raw:: latex

     \vspace{-3mm}

ABMD is commonly combined with the well-tempered (WT) scheme :cite:`Barducci:2008`. to 
obtain convergence of the bias potential within a finite simulation time.
In this approach, an additional factor


  .. raw:: latex

     \vspace{-5mm}

  .. math::
     \omega = \exp\!\left(-\frac{U(t \mid \xi)}{k_B T'}\right)


  .. raw:: latex

     \vspace{-3mm}

is introduced, where :math:`T'` denotes the well-tempered temperature.
This factor modifies the update rate of the bias potential, resulting in


  .. raw:: latex

     \vspace{-5mm}

  .. math::
     U_m(t+\Delta t) = U_m(t) + \omega\,\Delta t\,\frac{k_B T}{\tau_F}\,G\!\left(\frac{\sigma}{\Delta\xi}-m\right).


  .. raw:: latex

     \vspace{-3mm}


ABMD has four commonly used extensions: multiple-walker ABMD (MW-ABMD) 
:cite:`Babin:2008,Babin:2009`, replica-exchange ABMD (RE-ABMD),
:cite:`Babin:2008,Babin:2009`,
generalized replica-exchange with soulte tempering ABMD (gREST-ABMD) 
[in preparation], and bias-exchange ABMD (BE-ABMD) :cite:`Ito:2024`.
In MW-ABMD, multiple simulations are performed simultaneously while 
sharing a common bias potential :math:`U(\xi \mid t)`.
The shared bias potential is updated at every MD step according to


  .. raw:: latex

     \vspace{-5mm}

  .. math::
  
     \frac{\partial}{\partial t} U(\xi \mid t)
     =
     \frac{k_B T}{\tau_F}
     \sum_{\alpha}
     G\!\left[\xi-\sigma\left(\mathbf{r}_1^{\alpha},\ldots,\mathbf{r}_N^{\alpha}\right)\right].


  .. raw:: latex

     \vspace{-3mm}


Here, :math:`N` denotes the total number of atoms and :math:`\mathbf{r}`
represents the atomic positions in each MD trajectory labeled by :math:`\alpha`.
RE-ABMD combines temperature replica-exchange MD (T-REMD) with ABMD in 
order to improve conformational sampling using multiple replicas simulated 
at different target temperatures. In BE-ABMD, replicas exchange different 
bias potentials :math:`U(\xi \mid t)` between replicas :math:`a` and :math:`b`
based on the Metropolis criterion,


  .. raw:: latex

     \vspace{-5mm}

  .. math::

     P^{ab} = \min\left(1,\exp(-\Delta)\right),


  .. raw:: latex

     \vspace{-3mm}


where


  .. raw:: latex

     \vspace{-5mm}


  .. math::

     \Delta =
     (\beta_b-\beta_a)(E^a-E^b)
     +
     \beta_a\left[U^a(\xi^b)-U^a(\xi^a)\right]
     -
     \beta_b\left[U^b(\xi^b)-U^b(\xi^a)\right].

  .. raw:: latex

     \vspace{-3mm}


Here, :math:`P^{ab}` is the exchange probability and :math:`E` denotes 
the potential energy. These methods can be applied 
to MD simulations with all-atom (AA), coarse-grained (CG), and QM/MM 
models using multiple collective variables available in ATDYN and SPDYN.

-----------------------------------------------------------------------

**tau** *Real*

  **Default : 0.0** (unit: ps)

  Characteristic time controlling the magnitude of the bias added at each step.

**wt_temperature** *Real*

  **Default : 0.0** (unit: K)

  Well-tempered temperature controlling the decay of the bias update factor.

**monitor_period** *Integer*

  **Default: 0**

  Period of output of bias potential file.

**monitor_period_spline** *Integer*

  **Default: 0**

  Period of output of b-spline coeffcients file.


Examples
========

Example of a ABMD simulation for production.
::

  [ABMD]
  tau                    = 100.0
  wt_temperature         = 3000.0
  monitor_period         = 1000
  monitor_period_spline  = 10000

Example of a MW-ABMD simulation for production.
::

  [ABMD]
  tau                    = 100.0
  wt_temperature         = 3000.0
  monitor_period         = 1000
  monitor_period_spline  = 10000

  [REMD]
  dimension        = 1
  exchange_period  = 0
  type1            = TEMPERATURE
  nreplica1        = 4
  parameters1      = 300.00 300.00 300.00 300.00

Example of a RE-ABMD simulation for production.
::

  [ABMD]
  tau                    = 100.0
  wt_temperature         = 3000.0
  monitor_period         = 1000
  monitor_period_spline  = 10000

  [REMD]
  dimension        = 1
  exchange_period  = 1000
  type1            = TEMPERATURE
  nreplica1        = 4
  parameters1      = 300.00 350.00 400.00 450.00

Example of a gREST-ABMD simulation for production.
::

  [ABMD]
  tau                    = 100.0
  wt_temperature         = 3000.0
  monitor_period         = 1000
  monitor_period_spline  = 10000

  [REMD]
  dimension        = 1
  exchange_period  = 1000
  type1            = REST
  nreplica1        = 4
  parameters1      = 300.00 350.00 400.00 450.00
  select_index1    = 1
  param_type1      = C LJ D CM
  analysis_grest   = no

  [SELECTION]
  group1           = ((segid:PROA or segid:PROB) and (rno:79 or rno:83)) or (rnam:DMPC and (an:N or an:P))