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
