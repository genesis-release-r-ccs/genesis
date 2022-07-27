.. highlight:: bash
.. _boundary:

=======================================================================
Boundary section
=======================================================================

Boundary condition
===========================

**type** *PBC / NOBC*

  **Default : PBC**

  Type of boundary condition.

  * **PBC**: Periodic boundary condition (rectangular or cubic box).
  * **NOBC**: Non-boundary condition (vacuum system).
    In **ATDYN**, **NOBC** is applicable to various force fields and models.
    However, in **SPDYN**, it cannot be used.

**box_size_x** *Real*

  **Default : N/A** (unit : :math:`\text{\AA}`)

  Box size along the x dimension.

**box_size_y** *Real*

  **Default : N/A** (unit : :math:`\text{\AA}`)

  Box size along the y dimension.

**box_size_z** *Real*

  **Default : N/A** (unit : :math:`\text{\AA}`)

  Box size along the z dimension.

**local_pbc** *Yes / NO*

  **Default : NO**

  Logical flag to use PBC coordinates in force/energy calculations. (only available for coarse-grained models and PBC in **ATDYN**, In **SPDYN**, **local_pbc** is automatically defined in CHARMM and GROAMBER force fields)


.. note::
  If the simulation system has a periodic boundary condition (**PBC**),
  the user must specify the box size in the control file
  (in the energy minimization stage in most cases).
  During the simulations, box size is saved in a restart file. 
  If the restart file is used as an input for the subsequent simulation, 
  the box size is overwritten with the restart information. 
  Note that in this case the box size given in the control file is ignored.


Domain decomposition
====================

With domain decomposition, the total domain number is the same as the number of MPIs. Each domain has at least 
two cells in each dimension. In GENESIS 1.0 to 1.7.1, the minimum cell size is decided 
according to the working conditions. In the case of minimization or NVE/NVT MD without constraints, it is 

  .. math::

     {\frac{\text{pairlistlidst}}{2}}

If we apply constraints, it becomes

  .. math::

     {\frac{\text{pairlistlidst}+2.0}{2}}

It is further increased when we perform NPT simulation with constraints as follows:

  .. math::

     {\frac{\text{pairlistlidst}+2.6}{2}}

From GENESIS 2.0, the minimum cell size does not change by changing the working condition, and calculated as

  .. math::

     {\frac{\text{pairlistlidst}+2.0+\text{cell}\_\text{size}\_\text{buffer}}{2}}

-----------------------------------------------------------------------

**domain_x** *Integer*

  **Default : N/A (Optional)** (**SPDYN** only)

  Number of domains along the x dimension.

**domain_y** *Integer*

  **Default : N/A (Optional)** (**SPDYN** only)

  Number of domains along the y dimension.

**domain_z** *Integer*

  **Default : N/A (Optional)** (**SPDYN** only)

  Number of domains along the z dimension.

**cell_size_buffer** *Real*

  **Default : 0.6**
  Additional buffer size for determining the cell size.
  
.. note::
  If the number of domains (domain_x, domain_y, and domain_z) are not specified
  in the control file, they are automatically determined based on the
  number of MPI processes. 
  When the user specifies the number of domains explicitly,
  please make sure that the product of the domain numbers in each dimension
  (i.e., domain_x * domain_y * domain_z) is equal to the total number of
  MPI processes.


Spherical potential
===================

In MD simulations with **NOBC**, molecules may evaporate from a system, 
and once such an event happens, the molecule escapes in the vacuum with 
constant velocity to infinity. Therefore, it is useful to set a 
potential which pulls the molecule back to the system.

In **ATDYN**, the user can set a spherical potential,

  .. math::
     V & = k (r_i - r_b)^n \hspace{10pt} ( r_i > r_b) \\
       & = 0,              \hspace{10pt} ( r_i \leq r_b)

where :math:`k, n` and :math:`r_b` represent the force constant, an 
exponent, and the radius of the sphere, respectively, and :math:`r_i` 
is the distance between the :math:`i`-th atom and the center of the sphere,

  .. math::
     r_i = | \mathbf{x}_i - \mathbf{x}_0 | .


  .. figure:: _figures/fig10_1.pdf
     :width: 30 %
     :align: center
     :name: fig10_1
     :alt: 

     An illustration of a combination of two spherical potentials (black thin circles), 
     which pulls back atoms that are out of the range towards the center of 
     sphere (1 and 2).

Multiple spheres with different centers and radii can be combined to construct 
the potential; for example, two spheres are combined in :numref:`fig10_1`. The 
atoms that went out of the sphere (thin line) are pulled back to the nearest 
center; the red atom to center 1 and the blue atoms to center 2.

The coordinates of the center can be specified in two ways. The first is
to set the center to the position of the atoms in the initial structure (pdbfile) 
using **[SELECTOR]**.  The other way is to directly specify the coordinates of the 
center. See the description of options and the examples below for details.

The following options are available when using the spherical potential:

-----------------------------------------------------------------------

**spherical_pot** *YES / NO*

  **Default : NO** 

  If YES (with type=NOBC), use the spherical boundary potential.

**constant** *Real*

  **Default : 10.0** (unit : :math:`\mathrm{kcal mol}^{-1}`)

  The force constant of the potential.

**exponent** *Integer*

  **Default : 2** 

  The exponent of the potential.

**nindex** *Integer*

  **Default : 0** 

  The number of indices, used with center_select_index *N*.

**center_select_index** `N` *Integer*

  **Default : N/A** 

  The index of center in **[SELECTION]**.

**nfunction** *Integer*

  **Default : 0** 

  The number of function, used with center *N*.

**center** `N` *Real* :math:`\times 3`

  **Default : N/A** 

  The xyz coordinates of the center.

**radius** `N`

  **Default : 0.0** (unit : :math:`\text{\AA}`)

  The radius of the sphere.

**fixatom** *YES / NO*

  **Default : YES**

  Atoms out of the sphere in the input structure are fixed.

**fix_layer** *Real*

  **Default : 1.0** (unit : :math:`\text{\AA}`)

  If fixatom = YES, atoms within this distance from the potential in the input structure are also fixed.

**restart** *YES / NO*

  **Default : YES**

  Use the information in the restart file.

.. note::
   The information of the sphere and fixed atoms are saved in a restart file. If the 
   information exists in rstfile, the options for the spherical potential in **[BOUNDARY]** will be 
   ignored. If you want to re-set the potential, you need to specify **restart** = NO. 
  

Examples
========

* Simulations in the gas-phase:
  :: 
   [BOUNDARY]
   type       = NOBC      # non-periodic system

* Simulations with the periodic boundary condition, where the
  box size is set to 64 x 64 x 64.
  In this case, the user should not use a restart file as an input,
  because the box size in the control is overwritten with that in the restart file.
  :: 
   [BOUNDARY]
   type       = PBC       # periodic boundary condition
   box_size_x = 64.0      # Box size in the x dimension (Ang)
   box_size_y = 64.0      # Box size in the y dimension (Ang)
   box_size_z = 64.0      # Box size in the z dimension (Ang)

* Simulations with two spherical potentials around atom number 1 and 100 with a radius of 22.0 :math:`\text{\AA}`. 
  ::
   [BOUNDARY]
   type          = NOBC
   spherical_pot = yes
   constant      = 2.0
   exponent      = 2
   nindex        = 1
   center_select_index1 = 2
   radius1       = 22.0
   fix_layer     = 0.0
   fixatom       = no

   [SELECTION]
   ...
   group2         = ano:1 or ano:100

  .. note::
     Be careful not to set too many spheres because it may slow down 
     the performance. If you want to set the spheres around a protein, instead
     of specifying all atoms in a protein, select part of the atoms, for example, by
     :: 
       group2 = segid:PROA and an:CA

* For simulations with two spherical potentials, the center 
  coordinates are explicitly set by **center1** and **center2**.
  With **fixatom** =yes and **fix_layer** =1.0 :math:`\text{\AA}`, atoms
  farther than 34 :math:`\text{\AA}` from the centers are fixed. 
  ::
    [BOUNDARY]
    type          = NOBC      # [PBC,NOBC]
    spherical_pot = yes
    constant      = 10.0
    exponent      = 2
    nfunctions    = 2
    center1       =  17.0, 0.0, 0.0    # [x,y,z]
    radius1       = 35.0
    center2       = -17.0, 0.0, 0.0    # [x,y,z]
    radius2       = 35.0
    fixatom       = YES
    fix_layer     = 1.0
