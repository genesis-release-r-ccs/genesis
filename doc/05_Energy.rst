.. highlight:: bash
.. _energy:

=======================================================================
Energy section
=======================================================================


Force fields
=======================================================================

In general, potential energy function is given by:

  .. raw:: latex
     
     \vspace{-5mm}

  .. math::  
     E(r) & = \sum_{\mathrm{bond}}K_b(b-b_0)^2+\sum_{\mathrm{angle}}K_{\theta}(\theta-\theta_0)^2  \\
          & + \sum_{\mathrm{dihedral}}K_{\phi}(1+\cos\left(n\phi-\delta)\right)+\sum_{\mathrm{improper}}K_{\phi_i}(\phi_i-\phi_{i,0})^2 \\
          & + \sum_{\mathrm{nonbond}}\epsilon\left[\left(\frac{R_{min,ij}}{r_{ij}}\right)^{12}-2\left(\frac{R_{min,ij}}{r_{ij}}\right)^6\right]+\sum_{\mathrm{nonbond}}\frac{q_iq_j}{\epsilon_1r_{ij}}

  .. raw:: latex
     
     \vspace{-3mm}

where :math:`K_b`, :math:`K_{\theta}`, :math:`K_{\phi}`,
and :math:`K_{\phi_i}` are the force constant of the bond, angle, 
dihedral angle, and improper dihedral angle term, respectively, and
:math:`b_0`, :math:`\theta_0`, :math:`\phi_0`, and :math:`\phi_{i,0}` 
are corresponding equilibrium values.
:math:`\delta` is a phase shift of the dihedral angle potential,
:math:`\epsilon` is a Lennard-Jones potential well depth,
:math:`R_{min,ij}` is a distance of the Lennard-Jones potential minimum,
:math:`q_i` is an atomic charge, :math:`\epsilon_1` is an effective
dielectric constant, and :math:`r_{ij}` is a distance between two atoms. 
The detailed formula and parameters in the potential energy function 
depend on the force field and molecular model.

-----------------------------------------------------------------------

**forcefield** *CHARMM / CHARMM19 / AMBER / GROAMBER / GROMARTINI / KBGO / CAGO / AAGO / RESIDCG*

  **Default : CHARMM**

  Type of the force field used for energy and force calculation.
  For the AMBER force field, the scheme used in the GROMACS program package
  is availble in addition to that used in the AMBER package.
  In this case, calculation for the dispersion correction term and
  truncation of the non-bonded energy term are different between AMBER and GROMACS.

  * **CHARMM**: CHARMM force field with the all-atom model (CHARMM22, 27, 36, 36m) :cite:`MacKerell:1998` :cite:`Klauda:2010` :cite:`Best:2012` :cite:`Huang:2013`
  * **CHARMM19**: CHARMM force field with the united-atom model (**ATDYN** only)
  * **AMBER**: AMBER force field with the original AMBER scheme :cite:`Case:2005`
  * **GROAMBER**: AMBER force field with the GROMACS scheme
  * **GROMARTINI**: MARTINI model :cite:`Marrink:2007` :cite:`Monticelli:2008`
  * **KBGO**: model by Karanicolas and Brooks :cite:`Karanicolas:2002hl` :cite:`karanicolas2003improved` (**ATDYN** only)
  * **CAGO**: C\ :math:`\alpha` Go-model :cite:`Clementi:2000` (**ATDYN** only)
  * **AAGO**: All-atom Go-model :cite:`Whitford:2009`
  * **RESIDCG**: Residue-level coarse-grained models (**ATDYN** only)

**Note:** Recently, ff19SB is provided in AMBER force field :cite:`Tian:2020`.
ff19SB is recommended to use with OPC four-point water model, which is available only in **SPDYN**.
Therefore, please use **SPDYN** to make use of ff19SB. 

Non-bonded interactions
=======================================================================

Calculation of the non-bonded interaction is the most time consuming part in MD simulations.
Computational time for the non-bonded interaction terms without
any approximation is proportional to :math:`O(N^2)`. 
To reduce the computational cost, a cut-off approximation is introduced,
where the energy and force calculation is truncated at a given cut-off value (keyword *cutoffdist*).
Simple truncation at the cut-off distance leads to discontinuous
energy and forces. So it is necessary to introduce a polynomial function
(so called *switching function*) that smoothly turn off the interaction from
another given value (so called *switch cut-off*), which is generally applied
to the van der Waals interactions (keyword *switchdist*). 
There are two kinds of switching: "potential switch" and "force switch".
In **GENESIS**, potential switching is turned on as the default.
However, in the case of the AMBER force field, potential switching is still turned off,
since the original AMBER program package is not using the potential switching.
To turn on the "force switching", ``vdw_force_switch=YES`` must be specified.
Note that the cut-off scheme for the electrostatic energy term is different from
that for the van der Waals energy term, where the former uses a shift function.
Such shift is turned on when ``Electrostatic=Cutoff`` is specified.

-----------------------------------------------------------------------

**electrostatic** *CUTOFF / PME*

  **Default : PME**

  * **CUTOFF**: Non-bonded interactions including the van der Waals interaction are just truncated at *cutoffdist*.

  * **PME**: Particle mesh Ewald (PME) method is employed for long-range interactions. This option is only availabe in the periodic boundary condition.

**switchdist** *Real*

  **Default : 10.0** (unit :  :math:`\text{\AA}`)

  Switch-on distance for nonbonded interaction energy/force quenching.
  If *switchdist* is set to be equal to *cutoffdist*, switching can be turned off.
  Switching scheme depends on the selected force field, *vdw_shift*, and *vdw_force_switch* parameters. 
  In the case of AMBER force field, this switching must be disabled, because the switching function is not available.
  In the case of "forcefield = GROMARTINI" and "electrostatic = CUTOFF", *switchdist* is used only in the van der Waals potential energy. 
  The switching-on distance for the electrostatic energy is automatically defined as 0.0.

**cutoffdist** *Real*

  **Default : 12.0** (unit :  :math:`\text{\AA}`)

  Cut-off distance for the non-bonded interactions.
  This distance must be larger than *switchdist*, while smaller than *pairlistdist*. 
  In the case of the AMBER force field, this value must be equal to *switchdist*.

**pairlistdist** *Real*

  **Default : 13.5** (unit :  :math:`\text{\AA}`)

  Distance used to make a Verlet pair list for non-bonded interactions :cite:`Verlet:1967`.
  This distance must be larger than *cutoffdist*.

**dielec_const** *Real*

  **Default : 1.0**

  Dielectric constant of the system.
  In GENESIS, the distance dependent dielectric constant is not available
  except for a specific case like the implicit membrane/micelle models (IMM1/IMIC).
  Note that in the IMM1/IMIC models this parameter is neglected.

**vdw_force_switch** *YES / NO*

  **Default : NO**

  This paramter determines whether the force switch function for van der Waals
  interactions is employed or not. :cite:`Steinbach:1994`
  The users must take care about this parameter, when the CHARMM force field is used.
  Typically, "vdw_force_switch=YES" should be specified in the case of CHARMM36.

**vdw_shift** *YES / NO*

  **Default : NO**

  This parameter determines whether the energy shift for the van der Waals interactions
  is employed or not. If it is turned on, potential energy at the cut-off distance is
  shifted by a constant value so as to nullify the energy at that distance,
  instead of the default smooth quenching function.
  This parameter is available only when "forcefield = GROAMBER" or "forcefield = GROMARTINI".

**dispersion_corr** *NONE / ENERGY / EPRESS*

  **Default : NONE** (automatically set to **EPRESS** in the case of AMBER)

  This parameter determines how to deal with the long-range correction 
  about the cut-off for the van der Waals interactions.
  Note that the formula used for the correction is different
  between the GROMACS and AMBER schemes.
  In the case of the CHARMM force filed, "dispersion_corr=NONE" is always used.

  * **NONE**: No correction is carried out.
  
  * **ENERGY**: Only energy correction is carried out.
    
  * **EPRESS**: Both energy and internal pressure corrections are carried out.

**implicit_solvent** *NONE / GBSA / EEF1 / IMM1 / IMIC* (**ATDYN** only)

  **Default : NONE**

  Use implicit solvent or not.

  * **NONE**: Do not use implicit solvent model

  * **GBSA**: Use the GB/SA implicit water model (Only available with the CHARMM or AMBER all-atom force fields in non-boundary condition ("type=NOBC" in the **[BOUNDARY]** section). :cite:`Onufriev:2004` :cite:`Weiser:1999`

  * **EEF1**: Use the EEF1 implicit water model (Only available with the CHARMM force fields in NOBC) :cite:`Lazaridis:1999`

  * **IMM1**: Use the IMM1 implicit membrane model (Only available with the CHARMM force fields in NOBC) :cite:`Lazaridis:2003`

  * **IMIC**: Use the IMIC implicit micelle model (Only available with the CHARMM force fields in NOBC) :cite:`Mori:2020`

**contact_check** *YES / NO*

  **Default : NO**

  If this parameter is set to *YES*, length of all covalent bonds 
  as well as distance between non-bonded atom pairs are checked at the begining of the simulation.
  If long covalent bonds or clashing atoms are detected, those atom indexes are displayed in the log file.
  If *contact_check* is turned on, *nonb_limiter* is also automatically enabled.
  If the users want to turn on only "contact_check", 
  please specify "contact_check = YES" and "nonb_limiter = NO" explicitly.
  Note that this contact_check does not work in the parallel-io scheme.
  If you are using **SPDYN**, please see also *structure_check*.

**structure_check** *NONE / FIRST / DOMAIN* (**SPDYN** only)

  **Default : NONE**

  If this parameter is set to FIRST or DOMAIN, length of all covalent bonds
  as well as distance between non-bonded atom pairs are checked at the begininig or during the simulation.
  This option is similar to *contact_check*, 
  but has an improved capability when the parallel-io scheme is employed.
  In **SPDYN**, we recommend the users to use this option instead of *contact_check*.
  Since the structure check spends additional computational time,
  the users had better turn off this option in the production run.

  * **NONE**: Do not check the structure

  * **FIRST**: Check the structure only at the beginning of the simulation

  * **DOMAIN**: Check the structure whenever the pairlist is updated

**nonb_limiter** *YES / NO*

  **Default : NO** (automatically set to be equal to **contact_check**)

  If this parameter is set to *YES*,
  large force caused by the atomic clash is suppressed during the simulation.
  Here, the atomic clash can be defined by *minimum_contact* (see below).
  If "contact_check = YES" is specified, this parameter is automatically set to "YES".
  If the users want to turn on only "contact_check",
  please specify "contact_check = YES" and "nonb_limiter = NO" explicitly.
  This option is basically useful for the energy minimization or equilibration of the system.
  However, we strongly recommend the users to turn off this option in the production run,
  because suppression of large forces is an "unphysical" manipulation to avoid unstable simulations.

**minimum_contact** *Real*

  **Default : 0.5** (unit :  :math:`\text{\AA}`)

  This parameter defines the clash distance, when ``contact_check = YES`` is specified.
  If the distance between the non-bonded atoms is less than this value,
  energy and force are computed using this distance instead of the actual distance.

**nonbond_kernel**  *AUTOSELECT / GENERIC / FUGAKU / INTEL / GPU*

  **Default: AUTOSELECT**

  If this parameter is set to AUTOSELECT, the program automatically decide the kernel
  for the real-space non-bonded interaction. Please remember that you should compile with GPU
  option to define GPU as nonboned_kernel.

Particle mesh Ewald method
=======================================================================

Electrostatic energy in the conventional Ewald sum method is expressed as:

  .. raw:: latex
     
     \vspace{-5mm}

  .. math::
     E_{elec}=\sum_{i<j}\frac{q_i q_j}{\epsilon_1} \frac{\text{erfc}(\alpha r_{ij})}{r_{ij}}+\frac{2\pi}{V}\sum_{{\left|\mathbf{G}\right|}^2\ne 0}\frac{\exp(-{\left|\mathbf{G}\right|}^2/4\alpha^2)}{\left|\mathbf{G}\right|^2}\sum_{ij}\frac{q_i q_j}{\epsilon_1}\exp(i{\mathbf{G}} \cdot \mathbf{r}_{ij})-\sum_{ij}\frac{q_i q_j}{\epsilon_1}\frac{\alpha}{\sqrt{\pi}}

  .. raw:: latex
     
     \vspace{-3mm}

where :math:`\mathbf{G}` is the three-dimensional grid vectors in reciprocal space. Here, the cut-off scheme can be used for the first term, because it decreases
rapidly as distance between atoms increases. The third term is so called
*self-energy*, and is calculated only once. The second term can be rewritten as:

  .. raw:: latex
     
     \vspace{-5mm}

  .. math::
     \sum_{{\left|\mathbf{G}\right|}^2\ne 0}\frac{\exp(-{\left|\mathbf{G}\right|}^2/4\alpha^2)}{\left|\mathbf{G}\right|^2} {\left|\mathbf{S}(\mathbf{G})\right|}^2

  .. raw:: latex
     
     \vspace{-3mm}

where the structure factor :math:`\mathbf{S}(\mathbf{G})` is defined as:

  .. raw:: latex
     
     \vspace{-5mm}

  .. math::
     \mathbf{S}(\mathbf{G})=\sum_i q_i \exp(i\mathbf{G} \cdot \mathbf{r}_i)

We cannot employ fast Fourier transformation (FFT) for the calculation of
:math:`\mathbf{S}(\mathbf{G})` since atomic positions are usually not
equally spaced. In the smooth particle mesh Ewald (PME) method
:cite:`Darden:1993in` :cite:`Essmann:1995vj`, this structure factor is
approximated by using cardinal B-spline interpolation as:

  .. raw:: latex
     
     \vspace{-5mm}

  .. math::
     \mathbf{S}(\mathbf{G})=\sum_i q_i \exp(i\mathbf{G} \cdot \mathbf{r}_i) \approx b_1(G_1)b_2(G_2)b_3(G_3)\mathbf{F}(\mathbf{Q})(G_1,G_2,G_3)

  .. raw:: latex
     
     \vspace{-3mm}

where :math:`b_1(G_1)`, :math:`b_2(G_2)`, and :math:`b_3(G_3)` are
the coefficients brought by the cardinal B-spline interpolation of order
:math:`n` and :math:`\mathbf{Q}` is a 3D tensor obtained by interpolating
atomic charges on the grids.
Since this :math:`\mathbf{Q}` has equally spaced structure, its Fourier
transformation, :math:`\mathbf{F}(\mathbf{Q})`, can be calculated by using
FFT in the PME method.

-----------------------------------------------------------------------

**pme_alpha** *Real or auto*

  **Default : auto**

  Exponent of complementary error function. 
  If ``pme_alpha=auto`` is specified, 
  the value is automatically determined from *cutoffdist* and *pme_alpha_tol*.

  *Note: The default of pme_alpha was 0.34 in GENESIS ver. 1.1.0 or former.*

**pme_alpha_tol** *Real*

  **Default : 1.0e-5**

  Tolerance to be used for determining *pme_alpha*, when ``pme_alpha=auto`` is specified.

**pme_nspline** *Integer*

  **Default : 4**

  B-spline interpolation order used for the evaluation of :math:`b_1 (G_1)`,
  :math:`b_2 (G_2)`, :math:`b_3 (G_3)`, and :math:`\mathbf{Q}`.
  The order must be :math:`>= 3`.

**pme_max_spacing** *Real*

  **Default : 1.2** (unit :  :math:`\text{\AA}`)

  Max PME grid size used in the automatic grid number determination.
  This parameter is used only when *pme_ngrid_x*, *pme_ngrid_y*,
  and *pme_ngrid_z* are not given in the control file.

**pme_ngrid_x** *Integer*

  **Default : N/A (Optional)**

  Number of FFT grid points along x dimension.
  If not specified, program will determine an appropriate number of grids
  using ``pme_max_spacing``.

**pme_ngrid_y** *Integer*

  **Default : N/A (Optional)**

  Number of FFT grid points along y dimension.
  If not specified, program will determine an appropriate number of grids
  using ``pme_max_spacing``.

**pme_ngrid_z** *Integer*

  **Default : N/A (Optional)**

  Number of FFT grid points along z dimension.
  If not specified, program will determine an appropriate number of grids
  using ``pme_max_spacing``.

**pme_multiple** *YES/NO* (**ATDYN** only)

  **Default : NO**

  IF pme_multiple is set to YES, MPI processes are divided into two groups
  to compute the PME real and reciprocal parts individually.

**pme_mul_ratio** *Integer* (**ATDYN** only)

  **Default : 1**

  Ratio of the MPI processors for real and reciprocal PME term computations
  (only used when "PME_multiple=YES" is specified).

**PME_scheme** *AUTOSELECT / OPT_1DALLTOAlL / NOOPT_1DALLTOALL / OPT_2DALLTOALL / NOOPT_2DALLTOALL* (**SPDYN** only)

  **Default : AUTOSELECTL**

  *AUTOSELECT* chooses the best scheme by running in setup procedure.
  Other schemes can be chosen manually according to your preference..
  For *1DALLTOALL* and *2DALLTOALL*, see ref :cite:`Jung:2016` for details.

**SPDYN** use OpenMP/MPI hybrid parallel fast Fourier transformation library, FFTE :cite:`FFTE:Online`.
The number of PME grid points must be multiples of 2, 3, and 5 due to the restriction of this library.
Moreover, in **SPDYN**, there are several additional rules, which depends on the number of processes, in PME grid numbers.
In **SPDYN**, we first define domain numbers in each dimension such that product of them equals to the total number of
MPI processors. Let us assume that the domain numbers in each dimension are
``domain_x``, ``domain_y``, and ``domain_z``. The restriction condition of the grid numbers are as follows:

1) OPT_1DALLTOALL or NOOPT_1DALLTOALL:

  ``pme_ngrid_x`` should be multiple of 2 :math:`\times` ``domain_x``.

  ``pme_ngrid_y`` should be multiple of ``domain_y`` :math:`\times` ``domain_z`` if ``domain_z`` is an even number and multiple of ``domain_y`` :math:`\times` ``domain_z`` :math:`\times` 2 otherwise.

  ``pme_ngrid_z`` should be multiple of ``domain_x`` :math:`\times` ``domain_z`` and multiple of ``domain_y`` :math:`\times` ``domain_z``.

2) OPT_2DALLTOALL or NOOPT_2DALLTOAlL:

  ``pme_ngrid_x`` should be multiple of 2 :math:`\times` ``domain_x``.

  ``pme_ngrid_y`` should be multiple of ``domain_y`` :math:`\times` ``domain_z`` if ``domain_z`` is an even number and multiple of ``domain_y`` :math:`\times` ``domain_z`` :math:`\times` 2 otherwise.

  ``pme_ngrid_z`` should be multiple of ``domain_x`` :math:`\times` ``domain_z``.

3) NOOPT_1DALLTOALL or NOOPT_2DALLTOALL:

  Cell size in each dimension divided by grid spacing should be greater than ``pme_nsplie``.

4) OPT_1DALLTOALL or OPT_2DALLTOALL:

  ``pme_ngrid_[x,y,z]`` divided by ``domain_[x,y,z]`` should be greater than :math:`\times` ``pme_nspline``

External electric field
=======================================================================

In GENESIS, we can apply constant external electric field on x, y, and z directions. For biological simulations, it can be used for applying electric potential across a membrane.   

-----------------------------------------------------------------------

**efield_x** *Real*

  **Default : N/A**

  Usage of external electric field in x direction (unit: V/:math:`\text{\AA}`)

**efield_y** *Real*

  **Default : N/A**

  Usage of external electric field in y direction (unit: V/:math:`\text{\AA}`)

**efield_z** *Real*

  **Default : N/A**

  Usage of external electric field in z direction (unit: V/:math:`\text{\AA}`)

**efield_virial** *Yes / No*

  **Default : No**

  Logical flag to assign efield contribution to virial term

**efield_normal** *Yes / No*

  **Default : No**

  Logical flag to adjust electric field strength according to system size in NPT simulations. It is used to apply constant electric potential difference across a memebrane.


Lookup table
=======================================================================

The following keywords are relevant when CHARMM/AMBER/GROAMBER force fields are used.
For a linearly-interpolating lookup table with :math:`r_v` (cutoff),
force at pairwisde distance :math:`r` is evaluated according to the unit interval of :math:`r_v^2/r^2`:
:cite:`Jung:2013` 

  .. raw:: latex
     
     \vspace{-5mm}

  .. math::
     F(r^2) \approx F_{\text{tab}}(L)+t(F_{\text{tab}}(L+1)-F_{\text{tab}}(L))

  .. raw:: latex
     
     \vspace{-3mm}

where

  .. raw:: latex
     
     \vspace{-5mm}

  .. math::
     L=\text{INT}(\text{Density} \times r_v^2/r^2)

  .. raw:: latex
     
     \vspace{-3mm}

and

  .. raw:: latex
     
     \vspace{-5mm}

  .. math::
     t=\text{Density} \times r_v^2/r^2-L

  .. raw:: latex
     
     \vspace{-3mm}

Linear interpolation is used if "Electrostatic=PME".

Density is the number of points per unit interval. Lookup table using cubic interpolation is different from that of linear interpolation. In the case of cubic interpolation, monotonic cubic Hermite polynomial  interpolation is used to impose the monotonicity of the energy value. Energy/gradients are evaluated as a function of :math:`r^2` :cite:`Nilsson:2009` using four basis functions for the cubic Hermite spline : :math:`h_{00}(t)`, :math:`h_{10}(t)`, :math:`h_{01}(t)`, :math:`h_{11}(t)`

  .. raw:: latex
     
     \vspace{-5mm}

  .. math::
     F(r^2) & \approx F_{\text{tab}}(L-1)h_{00}(t)+\frac{F_{\text{tab}}(L-2) + F_{\text{tab}}(L-1)}{2} h_{10} \\ & + F_{\text{tab}}(L)h_{10}(t) + \frac{F_{\text{tab}}(L-1)+F_{\text{tab}}(L)}{2} h_{11}(t)

  .. raw:: latex
     
     \vspace{-3mm}

where

  .. raw:: latex
     
     \vspace{-5mm}

  .. math::
     L=\text{INT}(\text{Density} \times r^2)

  .. raw:: latex
     
     \vspace{-3mm}

and

  .. raw:: latex
     
     \vspace{-5mm}

  .. math::
     t=\text{Density} \times r^2-L

  .. raw:: latex
     
     \vspace{-3mm}

Cubic iterpolation is used if "Electrostatic=Cutoff".


Generalized Born/Solvent-Accessible Surface-Area model
======================================================

Implicit solvent model is useful to reduce computational cost 
for the simulations of biomolecules :cite:`Mori:2016`. The GB/SA (Generalized
Born/Solvent accessible surface area) model is one of the popular
implicit solvent models, where the electrostatic contribution
to the solvation free energy (:math:`\Delta {G_{{\rm{elec}}}}`) 
is computed with the GB theory :cite:`Still:1990`, and the non-polar contribution (:math:`\Delta {G_{{\rm{np}}}}`)
is calculated from the solvent accessible surface area :cite:`Eisenberg:1986`.
In the GB theory, solvent molecules
surrounding the solute are approximated as a continuum that has
the dielectric constant of ~80. To date, various GB models have
been developed. In GENESIS, the OBC model :cite:`Onufriev:2004` and LCPO method :cite:`Weiser:1999` are
available in the calculations of the GB and SA energy terms, respectively.
Note that the GB/SA model is implemented in **ATDYN** but NOT **SPDYN**.
The solvation free energy is incorporated into the molecular
mechanics potential energy function as an effective energy term,
namely, :math:`U = U_{{\rm{FF}}} + \Delta {G_{{\rm{elec}}}} + \Delta {G_{{\rm{np}}}}`.

GB energy term
--------------

In the GB theory, the solvation free energy of solute is given by

  .. raw:: latex

     \vspace{-5mm}

  .. math::

    \Delta {G_{{\rm{elec}}}} =  - \frac{1}{2}\left\{ {\frac{1}{{{\varepsilon _{\rm{p}}}}} - \frac{{\exp ( - \kappa {f_{ij}})}}{{{\varepsilon _{\rm{w}}}}}} \right\}\sum\limits_{i,j} {\frac{{{q_i}{q_j}}}{{{f_{ij}}}}},

  .. raw:: latex

     \vspace{-3mm}

where :math:`\varepsilon_{{\rm{p}}}` and :math:`\varepsilon_{{\rm{w}}}` are the dielectric constants of
solute and solvent, respectively, :math:`q_i` and :math:`q_j` are the partial charges 
on the *i*-th and *j*-th atoms, respectively. :math:`\kappa` is the inverse of Debye length.
:math:`f_{ij}` is the effective distance between the *i*- and *j*-th atoms,
which depends on the degree of burial of the atoms, and is given by

  .. raw:: latex

     \vspace{-5mm}
  .. math::

    {f_{ij}} = \sqrt {r_{ij}^2 + {R_i}{R_j}\exp \left( {\frac{{ - r_{ij}^2}}{{4{R_i}{R_j}}}} \right)}.

  .. raw:: latex

     \vspace{-3mm}

Here, :math:`r_{ij}` is the actual distance between the *i*- and *j*-th atoms, 
and :math:`R_i` is the effective Born radius of the *i*-th atom,
which is typically estimated in the Coulomb field approximation by

  .. raw:: latex

     \vspace{-5mm}

  .. math::

     {\frac{1}{R_i}} = {\frac{1}{{{\rho _i}}} - \frac{1}{{4\pi }}\int_{{\rm{solute}},{\rm{ }}r > {\rho _i}} {\frac{1}{{{r^4}}}} dV}.

  .. raw:: latex

     \vspace{-3mm}

:math:`\rho_i` is the radius of the *i*-th atom (mostly set to the atom’s van der Waals radius), 
and the integral is carried out over the volume inside the solute but outside the *i*-th atom.
In the case of an isolated ion, :math:`R_i` is equal to its van der Waals radius.
On the other hand, if the atom is buried inside a solute, :math:`R_i` becomes larger,
resulting in larger :math:`f_{ij}`.
In the OBC model, the effective Born radius is approximated as

  .. raw:: latex

     \vspace{-5mm}

  .. math::

     {\frac{1}{R_i}} = {\frac{1}{{\tilde \rho _i}} - \frac{1}{{{\rho _i}}}\tanh ({\alpha}{\Psi _i} - {\beta}\Psi _i^2 + {\gamma}\Psi _i^3)},

  .. raw:: latex

     \vspace{-3mm}

where :math:`{\tilde \rho _i}` is defied as :math:`\rho_i - \rho_0` (intrinsic offset), and
:math:`\Psi_{i}` describes the dgree of burial of the solute atom,
which is calculated from the pairwise descreening
function: :math:`{\Psi _i} = {{\tilde \rho _i}}\sum\limits_j {{H(r_{ij})}}` :cite:`Schaefer:1990`.

SA energy term
--------------

In general, the non-polar contribution to the solvation free energy is calculated by

  .. raw:: latex

     \vspace{-5mm}

  .. math::

    \Delta {G_{{\rm{np}}}} = \sum\limits_i {{\gamma _i}{A_i}},

  .. raw:: latex

     \vspace{-3mm}

where :math:`\gamma` is the surface tension coefficient, and :math:`A_i` is the surface area of the *i*-th atom.
In the LCPO method, :math:`A_i` is calculated from a linear combination of the overlaps between the
neighboring atoms, given by

  .. raw:: latex

     \vspace{-5mm}

  .. math::

    {A_i} = {P_{1i}}4\pi R_i^2 + {P_{2i}}\sum\limits_{j = 1}^n {{A_{ij}}}  + {P_{3i}}\sum\limits_{j = 1}^n {\sum\limits_{k = 1}^m {{A_{jk}}} }  + {P_{4i}}\sum\limits_{j = 1}^n {\left[ {{A_{ij}}\sum\limits_{j = 1}^n {\sum\limits_{k = 1}^m {{A_{jk}}} } } \right]}.

  .. raw:: latex

     \vspace{-3mm}

:math:`P_{1-4}` are the empirical parameters determined for each atom type,
:math:`R_i` is the radius of the *i*-th atom + probe radius (typically 1.4 :math:`\text{\AA}`),
and :math:`A_{ij}` is the area of the *i*-th atom buried inside the *j*-th atom, given by

  .. raw:: latex

     \vspace{-5mm}

  .. math::

    {A_{ij}} = 2\pi {R_i}\left( {{R_i} - \frac{{{r_{ij}}}}{2} - \frac{{R_i^2 - R_j^2}}{{2{r_{ij}}}}} \right)

  .. raw:: latex

     \vspace{-3mm}

where :math:`r_{ij}` is the distance between the *i*- and *j*-th atoms.

------------------------------------------------------------------

**gbsa_eps_solvent** *Real*

  **Default : 78.5**
  
  Dielectric constant of solvent :math:`\varepsilon_{{\rm{w}}}`.

**gbsa_eps_solute** *Real*

  **Default : 1.0**
  
  Dielectric constant of solute :math:`\varepsilon_{{\rm{p}}}`.

**gbsa_alpha** *Real*

  **Default : 1.0**

  The empirical parameter :math:`\alpha` in the equation for the effective Born radius calculation.
  "gbsa_alpha=0.8" for OBC1, and "gbsa_alpha=1.0" for OBC2.

**gbsa_beta** *Real*

  **Default : 0.8**

  The empirical parameter :math:`\beta` in the equation for the effective Born radius calculation.
  "gbsa_beta=0.0" for OBC1, and "gbsa_beta=0.8" for OBC2.

**gbsa_gamma** *Real*

  **Default : 4.85**

  The empirical parameter :math:`\gamma` in the equation for the effective Born radius calculation.
  "gbsa_gamma=2.91" for OBC1, and "gbsa_gamma=4.85" for OBC2.

**gbsa_salt_cons** *Real*

  **Default : 0.2** (unit : mol/L)

  Concentration of the monovalent salt solution.

**gbsa_vdw_offset** *Real*

  **Default : 0.09** (unit :  :math:`\text{\AA}`)

  Intrinsic offset :math:`\rho_0` for the van der Waals radius.

**gbsa_surf_tens** *Real*

  **Default : 0.005** (unit : kcal/mol/:math:`\text{\AA}^2`)

  Surface tension coefficient :math:`\gamma` in the SA energy term.


.. note::
  Debye length is calculated by
  :math:`{\kappa ^{ - 1}} = \sqrt {{\varepsilon _0}{\varepsilon _w}{k_B}T/2{N_A}{e^2}I}`, where
  *T* is automatically set to the target temperature specified in the **[DYNAMICS]** section.
  In the case of the energy minimization, *T* = 298.15 K is used. In the T-REMD simulations with
  the GB/SA model, each replica has an individual Debye length depending on the assigned temperature.


EEF1, IMM1, and IMIC implicit solvent models
============================================

In the EEF1 implicit solvent model :cite:`Lazaridis:1999`, the effective energy *W* of a solute molecule is defined as the sum of
the molecular mechanics potential energy :math:`E_{\rm{MM}}` and solvation free energy :math:`\Delta {G_{{\rm{solv}}}}`, given by

  .. raw:: latex

     \vspace{-5mm}

  .. math::

     W = E_{\rm{MM}} + \Delta {G_{{\rm{solv}}}},

  .. raw:: latex

     \vspace{-3mm}

where

  .. raw:: latex

     \vspace{-5mm}

  .. math::

     \Delta {G_{{\rm{solv}}}} = \sum\limits_i {\Delta {G_i^{{\rm{ref}}}}} - \sum\limits_i {\sum\limits_{j \neq i} {g_i(r_{ij})V_j}},

     g_i(r_{ij}) = {\frac{\Delta {G_i^{{\rm{free}}}}}{2\pi \sqrt{\pi} \lambda_i r_{ij}^2}} \exp \left\{ - \left( {\frac{r_{ij} - R_i}{\lambda_i}} \right)^2 \right\}.

  .. raw:: latex

     \vspace{-3mm}

:math:`r_{ij}` is the distance between atoms *i* and *j*, and :math:`V_{j}` is the volume
of the *j*-th atom. The function :math:`g_{i}` is the density of the solvation free energy of the *i*-th atom,
defined with the van der Waals radius :math:`R_{i}` and thickness of the first hydration shell :math:`\lambda_{i}`.
:math:`\Delta {G_i^{{\rm{ref}}}}` is the solvation free energy of the atom when it is fully exposed to solvent.
:math:`\Delta {G_i^{{\rm{free}}}}` is similar to :math:`\Delta {G_i^{{\rm{ref}}}}`,
but is determined to satisfy the zero solvation energy of deeply buried atoms.

In the IMM1 implicit membrane :cite:`Lazaridis:2003` and IMIC implicit micelle models :cite:`Mori:2020`,
:math:`\Delta {G_i^{{\rm{free}}}}` as well as :math:`\Delta {G_i^{{\rm{ref}}}}` are defined as a combination
of the solvation free energies of the *i*-th solute atom in water and cyclohexane:

  .. raw:: latex

     \vspace{-5mm}

  .. math::

     \Delta {G_i^{{\rm{ref}}}} &= f_i \Delta {G_i^{{\rm{ref,water}}}} + (1 - f_i) \Delta {G_i^{{\rm{ref,cyclohexane}}}},

  .. raw:: latex

     \vspace{-3mm}

where *f* is a function that describes the transition between water and cyclohexane phases.

In the IMM1 model, :math:`f_i` is given by the sigmoidal function:

  .. raw:: latex

     \vspace{-5mm}

  .. math::

     f(z_i^{\prime}) = {\frac{z_i^{\prime n}}{1 + z_i^{\prime n}}},

  .. raw:: latex

     \vspace{-3mm}

where :math:`z_i^{\prime} = |z_i|/(T/2)`, :math:`z_i` is the *z*-coordinate of the *i*-th atom,
and *T* is the membrane thickness. *n* controls the steepness of the membrane-water interface.
In the IMM1 model, the membrane is centered at :math:`z = 0`.

In the IMIC model, the following function is used for :math:`f_i`:

  .. raw:: latex

     \vspace{-5mm}

  .. math::

     f(d_i) = {\frac{1}{2}} \left\{ {\rm{tanh}}(sd_i) + 1 \right\},

  .. raw:: latex

     \vspace{-3mm}

where :math:`d_i` is the depth of the solute atom *i* from the micelle surface, and
*s* controls the steepness of the micelle-water interface.
The shape of the micelle is defined using a super-ellipsoid function:

  .. raw:: latex

     \vspace{-5mm}

  .. math::

     \left\{ \left( {\frac{|x|}{a}} \right)^{\frac{2}{m_2}}  + \left( {\frac{|y|}{b}} \right)^{\frac{2}{m_2}} \right\}^{\frac{m_2}{m_1}} + \left( {\frac{|z|}{c}} \right)^{\frac{2}{m_1}} = 1.

  .. raw:: latex

     \vspace{-3mm}

*a*, *b*, and *c* are the semi-axes of the super-ellipsoid, and
:math:`m_1` and :math:`m_2` determine the shape of the cross section in the super-ellipsoid.
In the case of :math:`m_1 = m_2 = 1`, the equation gives an ordinary ellipsoid.
If :math:`0 < m_1 < 1` and :math:`m_2 = 1`, the cross section in a plane perpendicular to
the *XY*-plane is expanded, keeping the semi-axes at the given lengths,
and the shape also resembles a bicelle or nanodisc.
If :math:`m_1 = 1` and :math:`0 < m_2 < 1`, the cross section in a plane parallel to the *XY*-plane is expanded.
The shape becomes close to rectangle as both :math:`m_1` and :math:`m_2` decrease.
Note that :math:`m < 0` or :math:`m > 1` is not allowed, because it produces a non-micelle-like shape resembling an octahedron.
In the IMIC model, the micelle is centered at the origin of the system :math:`(x, y, z) = (0, 0, 0)`.

In the IMM1 and IMIC models, a distance-dependent dielectric constant is used for the electrostatic interactions.
The dielectric constant depends on the positions of interacting atoms with respect to the membrane/micelle surface, defined as


  .. raw:: latex

     \vspace{-5mm}

  .. math::

     \epsilon = r^{ p + (1 - p) {\sqrt {f_i f_j}}},

  .. raw:: latex

     \vspace{-3mm}

where *r* is the distance between the *i*-th and *j*-th atoms, and *p* is an empirical parameter to adjust strength
of the interactions (*p* = 0.85 for CHARMM19 and 0.91 for CHARMM36).
Far from the membrane/micelle surface, the dielectric constant :math:`\epsilon` is close to *r*, corresponding to the EEF1 model,
while in the membrane/micelle center, it provides strengthened interactions.
The IMIC model is nearly equivalent to the IMM1 model when *a* and :math:`b \rightarrow \infty`  and *c* is half membrane thickness. 

In the control file of GENESIS, the following parameters are specified,
and the other parameters such as *V*, :math:`\Delta {G_i^{{\rm{ref}}}}`,
:math:`\Delta {G_i^{{\rm{free}}}}`, and :math:`\lambda` in the above equations
are read from the **eef1file**, which is set in the **[INPUT]** section (see :ref:`input`).
*R* is read from the **parfile**.

------------------------------------------------------------------

**imm1_memb_thick** *Real*

  **Default : 27.0** (unit :  :math:`\text{\AA}`)

  Membrane thickness *T* in IMM1

**imm1_exponent_n** *Real*

  **Default : 10**

  Steepness parameter *n* in IMM1

**imm1_factor_a** *Real*

  **Default : 0.91**

  Adjustable empirical parameter *p* in IMM1 and IMIC.
  *p* = 0.85 and 0.91 are recommended for CHARMM19 and CHARMM36, respectively.

**imm1_make_pore**

  **Default : NO**

  Use IMM1-pore model :cite:`Lazaridis:2005`

**imm1_pore_radius** *Real*

  **Default : 5.0** (unit :  :math:`\text{\AA}`)

  Aqueous pore radius in the IMM1-pore model

**imic_axis_a** *Real*

  **Default : 18.0** (unit :  :math:`\text{\AA}`)

  Semi-axis *a* of the super-ellipsoid in IMIC

**imic_axis_b** *Real*

  **Default : 18.0** (unit :  :math:`\text{\AA}`)

  Semi-axis *b* of the super-ellipsoid in IMIC

**imic_axis_c** *Real*

  **Default : 18.0** (unit :  :math:`\text{\AA}`)

  Semi-axis *c* of the super-ellipsoid in IMIC

**imic_exponent_m1** *Real*

  **Default : 1.0**

  Expansion parameter :math:`m_1` in IMIC

**imic_exponent_m2** *Real*

  **Default : 1.0**

  Expansion parameter :math:`m_2` in IMIC

**imic_steepness** *Real*

  **Default : 0.5**

  Steepness parameter *s* in IMIC


Residue-level coarse-grained models
===================================

**Note:** here we only briefly describe the potential energy functions and
corresponding options of the available CG models in
GENESIS :cite:`Tan:PLOSCB2022`.  To prepare the topology and coordinate files,
please make use of the "GENESIS-CG-tool" :cite:`Tan:PLOSCB2022`, which is
included in GENESIS as a submodule and can be found from
"src/analysis/cg_tool".  It is also deposited as a Github repository, with
which a wiki page of the GENESIS-CG-tool can be found:
https://github.com/noinil/genesis_cg_tool/wiki

AICG2+ protein model
--------------------

In the AICG2+ model for proteins, the energy function is given by: :cite:`Li:PNAS2014`

  .. raw:: latex

     \vspace{-5mm}

  .. math::

     E_{AICG2+}(\mathbf{r}) & = \sum_{\mathrm{bond}} K_b\left(b - b_{0}\right)^2 + V_{loc}^{flp} \\
          & + \sum_{\mathrm{1,3\ pair}} \varepsilon_{1,3} \exp\left(-\frac{(d-d_0)^2}{2w_{1,3}^2}\right) + \sum_{\mathrm{dihedral}} \varepsilon_{\phi}\exp \left(-\frac{(\phi - \phi_0)^2}{2w_{\phi}^2}\right) \\
          & + \sum_{\mathrm{native\ contact}} \varepsilon_{Go} \left[ 5\left(\frac{r_0}{r}\right)^{12} - 6 \left(\frac{r_0}{r}\right)^{10} \right] \\
          & + \sum_{\mathrm{non-native\ contact}} \varepsilon_{exv} \left[ \left(\frac{\sigma}{r}\right)^{12} - \frac{1}{2^{12}} \right].


  .. raw:: latex

     \vspace{-3mm}

where :math:`K_{b}`, :math:`\varepsilon_{1,3}`, :math:`\varepsilon_{\phi}`,
:math:`\varepsilon_{Go}`, and :math:`\varepsilon_{exv}` are the force constants
of the bond, 1,3 (the next neighbor) distance, dihedral angle, native contact,
and non-native contact terms, respectively; and :math:`b_{0}`, :math:`d_{0}`,
:math:`\phi_0`, and :math:`r_0` are the native values of the corresponding
terms.  :math:`\sigma` in the non-native contact term is the residue-type
dependent excluded volume radius. :math:`w_{1,3}` and :math:`w_{\phi}` are the
widths of the Gaussian-type local potentials.  The statistical local potential
:math:`V_{loc}^{flp}` is used to describe the intrinsical flexibility of
peptides and includes the following terms: :cite:`Terakawa:BJ2011`

  .. raw:: latex

     \vspace{-5mm}

  .. math::

     V_{loc}^{flp} = \sum_{\mathrm{angle}} -k_B T \ln \frac{P(\theta)}{\sin(\theta)}  + \sum_{\mathrm{dihedral}} - k_B T \ln P(\phi)

  .. raw:: latex

     \vspace{-3mm}

where :math:`k_B` is the Boltzmann constant, :math:`T` is the temperature,
and :math:`P(\theta)` and :math:`P(\phi)` are the probability distributions of
the angles and dihedral angles, respectively.


3SPN.2C DNA model
-----------------

Among the series of the 3SPN DNA models developed by de Pablo's group
:cite:`Freeman:2011,Hinckley:2013,Freeman:2014`, 3SPN.2C is the
one for reproducing sequence-dependent curvature of double-stranded DNA
(dsDNA). :cite:`Freeman:2014` In this model, each nucleotide is represented by
three CG particles, P (phosphate), S (sugar), and B (base), respectively.  The
potential energy is given by::cite:`Freeman:2014`

  .. raw:: latex

     \vspace{-5mm}

  .. math::

      E_{3SPN.2C}(\mathbf{r}) & = \sum_{\mathrm{bond}} k_b (r - r_{0})^2 + 100 k_b (r - r_{0})^4 \\
          & + \sum_{\mathrm{angle}} k_\theta (\theta - \theta_{0})^2 \\
          & + \sum_{\mathrm{backbone\ dihedral}} -k_{\phi, Gaussian} \exp\big( \frac{-(\phi - \phi_0)^2}{2\sigma_{\phi}^2} \big) + \sum_{\mathrm{dihedral}} k_{\phi, periodic} \big[ 1+\cos(\phi - \phi_0) \big] \\
          & + \sum_{\mathrm{bstk}} U_m^{rep}(\epsilon_{bs}, \alpha_{bs}, r) + f(K_{bs}, \Delta\theta_{bs}) U_m^{attr} (\epsilon_{bs}, \alpha_{bs}, r) \\
          & + \sum_{\mathrm{bp}} U_m^{rep}(\epsilon_{bp}, \alpha_{bp}, r) + \frac{1}{2} \big( 1+\cos(\Delta \phi_1) \big) f(K_{bp}, \Delta\theta_{1}) f(K_{bp}, \Delta\theta_{2}) U_m^{attr} (\epsilon_{bp}, \alpha_{bp}, r) \\
          & + \sum_{\mathrm{cstk}} f(K_{BP}, \Delta\theta_{3}) f(K_{cs}, \Delta\theta_{cs}) U_m^{attr} (\epsilon_{cs}, \alpha_{cs}, r) \\
          & + \sum_{exv} \epsilon_r\bigg[ \Big(\frac{\sigma}{r} \Big)^{12} - 2 \Big(\frac{\sigma}{r} \Big)^6 \bigg]+ \epsilon_r \\
          & + \sum_{ele} \frac{q_i q_j e^{-r/\lambda_D}}{4\pi \epsilon_0 \epsilon(T, C) r}.

  .. raw:: latex

     \vspace{-3mm}

where the first four lines are the local terms, for bond, angle, dihedral angle,
and the base-stacking between neighboring bases, respectively.  The non-local
terms include the base-pairing, cross-stacking, excluded volume, and
electrostatic interactions.  :math:`U_m^{rep}` and :math:`U_m^{attr}` are the
splited Morse-type potentials to describe the repulsive and attractive
interactions between DNA bases, while the :math:`f(K, \Delta\theta)` are
bell-shaped modulating functions to smooth out the potentials.  The excluded
volume interactions are considered with a cutoff at :math:`\sigma`.  The
electrostatic interactions are modeled with the Debye-Hückel theory.  For more
details, please refer to reference :cite:`Freeman:2014`.


Protein-DNA interaction models
------------------------------

We usually consider excluded volume and electrostatic terms as the
sequence-non-specific interactions between protein and DNA.  The excluded volume
term has the same form as the non-native contact in the AICG2+ model, and the
electrostatic term uses the Debye-Hückel potential as in the 3SPN.2C model. As
for the protein-DNA sequence-specific recognition, we used the PWMcos model for
the protein-DNA base interactions: :cite:`Tan:PWMcos`

  .. raw:: latex

     \vspace{-5mm}

  .. math::

      E_{PWMcos}(\mathbf{r}) & = \sum_{i} \sum_{j} \sum_{m}\left[ V_{m,j}(b_i, \mathbf{r}) +  V_{m',j}(b_i, \mathbf{r}) \right],

  .. raw:: latex

     \vspace{-3mm}

where :math:`i` is the index of DNA base, :math:`j` is the index of amino acid
residue (:math:`C_\alpha`) that is forming contact with DNA in the reference
(native) complex structure, and :math:`m` is the index of native contacts
between protein and DNA.  Each amino acid residue can form multiple native
contacts with DNA bases.  We consider contributions from both bases (indices
:math:`m` and :math:`m'`) in every base-pair.

The sequence dependent potential :math:`V_{m, j} (b_i, \mathbf{r})` is defined as:

  .. raw:: latex

     \vspace{-5mm}

  .. math::

     V_{m, j} (b_i, \mathbf{r}) = \varepsilon_{PWM}(b_i) \cdot U_{mj} (i, \mathbf{r}),

  .. raw:: latex

     \vspace{-3mm}

where :math:`N_m` is the total number of contacts formed with the base pair
:math:`m-m'`, and :math:`\varepsilon_{PWM} (b_i)` is based on the
position-weight-matrix and dependent on the base type of the :math:`i`-th base:

  .. raw:: latex

     \vspace{-5mm}

  .. math::

     \varepsilon_{PWM} = \gamma \left[  \frac{-k_B T}{N_m} \left( \log P_m(b) - \frac{1}{4} \sum_{b\in \{A, C, G, T\}} \log P_m(b) \right) + \varepsilon'  \right]

  .. raw:: latex

     \vspace{-3mm}

where :math:`\gamma` and :math:`\varepsilon'` are two hyperparameters, which can
be calibrated by comparing simulated quantities with experimental results.

The second term in the formula of :math:`V_{m, j} (b_i, \mathbf{r})` is a
modulating function:

  .. raw:: latex

     \vspace{-5mm}

  .. math::

     U_{mj}(i, \mathbf{r}) = f_{mj}(r)g_{1,mj} (\theta_1) g_{2, mj} (\theta_2) g_{3,mj} (\theta_{3})  

  .. raw:: latex

     \vspace{-3mm}

where :math:`r` is the distance between the :math:`j`-th :math:`C_\alpha` and
the :math:`i`-th DNA base.  :math:`\theta_1`, :math:`\theta_2`, and
:math:`\theta_3` are the angle of sugar-base-:math:`C_\alpha`, the one between
vectors :math:`m-1`-base-:math:`m+1`-base and base-:math:`C_\alpha`, and the one
between vectors :math:`j-1`-:math:`C_\alpha`-:math:`j+1`-:math:`C_\alpha` and
base-:math:`C_\alpha`.  :math:`f` is a Gaussian centered at the native value of
:math:`r`, and :math:`g` is a bell-shaped function centered the native values of
the :math:`\theta` s.  For more details, please refer to reference
:cite:`Tan:PWMcos`.

A similar potential is useful in the modeling of hydrogen-bonds formed between
proteins and DNA backbone phosphate groups, which has been successfully applied
in the studies of nucleosome dynamics :cite:`Niina:nuc2017,Brandani:nuc2018`. We
also provide this model in GENESIS.


HPS/KH IDR models
-----------------

GENESIS provides two models for intrinsically disordered proteins, namely the
HPS and KH IDR models :cite:`Dignon:plos2018`.  These two models share the same
potential function formulas and are different in the parameters.  The local
interactions include the bond terms.  The nonlocal interactions include the
electrostatic term (Debye-Hückel form, the same as described in the 3SPN.2C
section) and the Ashbaugh-Hatch type hydrophobicity-scale (HPS) term:
:cite:`Dignon:plos2018`

  .. raw:: latex

     \vspace{-5mm}

  .. math::

     \Phi(r) = \begin{cases} \Phi_{LJ} + (1+\lambda)\epsilon, & \mathrm{if}\ r \le 2^{1/6}\sigma \\ \lambda \Phi_{LJ}, & \mathrm{otherwise} \end{cases}

  .. raw:: latex

     \vspace{-3mm}

where :math:`\lambda` is the hydrophobicity value and :math:`\Phi_{LJ}` is the standard Lennard-Jones potential:

  .. raw:: latex

     \vspace{-5mm}

  .. math::

     \Phi_{LJ} =  4\epsilon \left[  \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6    \right].

  .. raw:: latex

     \vspace{-3mm}

For more details such as the parameter values in the HPS and KH models, please refer to reference :cite:`Dignon:plos2018`.

Notably, user can easily change the hydrophobicity parameters in the topology
files (see more descriptions in the wiki-page of the GENESIS-CG-tool).  


Structure and context-based RNA model
-------------------------------------

GENESIS also provides a structure and context based RNA model :cite:`Hori:JCTC2012`, which has the potential functions:

  .. raw:: latex

     \vspace{-5mm}

  .. math::

     E_{RNA}(\mathbf{r}) & = \sum_{\mathrm{bond}} K_b\left(b - b_{0}\right)^2 + V_{loc}^{flp} \\
          & + \sum_{\mathrm{angle}} k_a \left(\theta - \theta_0\right)^2 \\
          & + \sum_{\mathrm{dihedral}} k_{\phi} \big[ 1-\cos(\phi - \phi_0) \big] + \frac{1}{2}k_{\phi} \big[ 1-\cos\left(3(\phi - \phi_0)\right) \big] \\
          & + \sum_{\mathrm{native\ contact}} \varepsilon_{Go} \left[ 5\left(\frac{r_0}{r}\right)^{12} - 6 \left(\frac{r_0}{r}\right)^{10} \right] \\
          & + \sum_{\mathrm{non-native\ contact}} \varepsilon_{exv} \left[ \left(\frac{\sigma}{r}\right)^{12} - \frac{1}{2^{12}} \right] \\
          & + \sum_{ele} \frac{q_i q_j e^{-r/\lambda_D}}{4\pi \epsilon_0 \epsilon(T, C) r}.

  .. raw:: latex

     \vspace{-3mm}

where the first three lines are for the local terms (bond, angle, and dihedral
angles), the forth term is the structure-based Go-like potential for native
contacts, the fifth term is for the non-native contacts, and the sixth term is
the Debye-Hückel type electrostatic interaction.


------------------------------------------------------------------

**cg_cutoffdist_ele** *Real* (unit:  :math:`\text{\AA}`)

  **Default : 52.0**

  Cutoff for Debye-Hückel type electrostatic interactions.

**cg_cutoffdist_126** *Real* (unit:  :math:`\text{\AA}`)

  **Default : 39.0**

  Cutoff for 12-6 type Lennard-Jones interactions.

**cg_cutoffdist_DNAbp** *Real* (unit:  :math:`\text{\AA}`)

  **Default : 18.0**

  Cutoff for DNA base-pairing interactions.

**cg_pairlistdist_ele** *Real* (unit:  :math:`\text{\AA}`)

  **Default : 57.0**

  Distance for the pair list of Debye-Hückel type electrostatic interactions.

**cg_pairlistdist_126** *Real* (unit:  :math:`\text{\AA}`)

  **Default : 44.0**

  Distance for the pair list of 12-6 type Lennard-Jones interactions.

**cg_pairlistdist_PWMcos** *Real* (unit:  :math:`\text{\AA}`)

  **Default : 23.0**

  Distance for the pair list of PWMcos type protein-DNA interactions.

**cg_pairlistdist_DNAbp** *Real* (unit:  :math:`\text{\AA}`)

  **Default : 23.0**

  Distance for the pair list of DNA base-pairing interactions.

**cg_pairlistdist_exv** *Real* (unit:  :math:`\text{\AA}`)

  **Default : 15.0**

  Distance for the pair list of all the excluded volume interactions.

**cg_sol_temperature** *Real* (unit: K)

  **Default : 300.0**

  Solution temperature used for the calculation of Debye length and dielectric
  constant. **Note:** this option is only used when "temperature" in the "[
  ENSEMBLE ]" section is set to 0, otherwise the later is used.
  
**cg_sol_ionic_strength** *Real* (unit: M)

  **Default : 0.15**

  Ionic strength used for the calculation of Debye length and dielectric constant.

**cg_pro_DNA_ele_scale_Q** *Real* (unit: :math:`e^{-}`)

  **Default : -1.0**

  Charge of the CG phosphate in the 3SPN.2C DNA is set to this value when
  protein-DNA electrostatic interactions are calculated.  **Note:** this option
  does not change the DNA-DNA interaction, in which the phosphate charge is
  :math:`-0.6e^{-}`.

**cg_exv_sigma_scaling** *Real* 

  **Default : 1.0**

  The radii of CG particles for the excluded volume interactions (including
  non-native contacts) are scaled by this value.  **Note:** this option does not
  affect the DNA-DNA excluded volume interactions.

**cg_IDR_HPS_epsilon** *Real* (unit: kcal/mol)

  **Default : 0.2**

  The force constant :math:`\epsilon` in the HPS-IDR model.

**cg_infinite_DNA** *Boolean*

  **Default : NO**

  Use the infinite DNA model or not.


Examples
========

Simulation with the CHARMM36 force field in the periodic boundary condition:
:: 
  [ENERGY]
  forcefield       = CHARMM  # CHARMM force field
  electrostatic    = PME     # use Particle mesh Ewald method
  switchdist       = 10.0    # switch distance
  cutoffdist       = 12.0    # cutoff distance
  pairlistdist     = 13.5    # pair-list distance
  vdw_force_switch = YES     # force switch option for van der Waals
  pme_nspline      = 4       # order of B-spline in [PME]
  pme_max_spacing  = 1.2     # max grid spacing allowed 

Simulation with the AMBER force field in the periodic boundary condition:
:: 
  [ENERGY]
  forcefield       = AMBER   # AMBER force field
  electrostatic    = PME     # use Particle mesh Ewald method
  switchdist       = 8.0     # switch distance
  cutoffdist       = 8.0     # cutoff distance
  pairlistdist     = 9.5     # pair-list distance
  pme_nspline      = 4       # order of B-spline in [PME]
  pme_max_spacing  = 1.2     # max grid spacing allowed 
 
Recommended options in the case of energy minimization (see :ref:`minimize`) for the initial structure with the CHARMM36 force field:
:: 
  [ENERGY]
  forcefield       = CHARMM  # CHARMM force field
  electrostatic    = PME     # use Particle mesh Ewald method
  switchdist       = 10.0    # switch distance
  cutoffdist       = 12.0    # cutoff distance
  pairlistdist     = 13.5    # pair-list distance
  vdw_force_switch = YES     # force switch option for van der Waals
  pme_nspline      = 4       # order of B-spline in [PME]
  pme_max_spacing  = 1.2     # max grid spacing allowed
  contact_check    = YES     # check atomic clash
  nonb_limiter     = YES     # avoid failure due to atomic clash
  minimum_contact  = 0.5     # definition of atomic clash distance
 
Simulations with the GB/SA implicit solvent model:
:: 
  [ENERGY]
  forcefield       = CHARMM  # CHARMM force field
  electrostatic    = CUTOFF  # use cutoff scheme
  switchdist       = 23.0    # switch distance
  cutoffdist       = 25.0    # cutoff distance
  pairlistdist     = 27.0    # pair-list distance
  implicit_solvent = GBSA    # Turn on GBSA calculation
  gbsa_eps_solvent = 78.5    # solvent dielectric constant in GB
  gbsa_eps_solute  = 1.0     # solute dielectric constant in GB
  gbsa_salt_cons   = 0.2     # salt concentration (mol/L) in GB
  gbsa_surf_tens   = 0.005   # surface tension (kcal/mol/A^2) in SA

Simulations with the IMM1 implicit membrane model:
:: 
  [ENERGY]
  forcefield       = CHARMM  # CHARMM force field
  electrostatic    = CUTOFF  # use cutoff scheme
  switchdist       = 16.0    # switch distance
  cutoffdist       = 18.0    # cutoff distance
  pairlistdist     = 20.0    # pair-list distance
  implicit_solvent = IMM1    # Turn on IMM1 calculation
  imm1_memb_thick  = 27.0    # membrane thickness in IMM1

Simulations with the residue-level CG models
:: 
  [ENERGY]
  forcefield             = RESIDCG
  electrostatic          = CUTOFF
  cg_cutoffdist_ele      = 52.0
  cg_cutoffdist_126      = 39.0
  cg_cutoffdist_DNAbp    = 18.0
  cg_pairlistdist_ele    = 57.0
  cg_pairlistdist_126    = 44.0
  cg_pairlistdist_PWMcos = 23.0
  cg_pairlistdist_DNAbp  = 23.0
  cg_pairlistdist_exv    = 15.0
  cg_sol_ionic_strength  = 0.15
  cg_pro_DNA_ele_scale_Q = -1.0
  cg_IDR_HPS_epsilon     = 0.2
