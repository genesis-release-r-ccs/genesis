.. highlight:: bash
.. _rpath:

=======================================================================
RPATH section
=======================================================================

Reaction Path Search
====================

In the **[RPATH]** section, users can specify keywords for finding the 
reaction path. The path search is carried out in two modes: the 
minimum energy path (MEP) and the minimum free-energy path (MFEP). 
The former searches for the energetically most favorable pathway on the 
potential energy surface (PES), while the latter does the same on the 
free-energy surface (FES). 
The MEP search is used to find relatively fast 
processes such as chemical reactions (very likely along with QM/MM), in
which the environment can be regarded more or less rigid. 
On the other hand, the MFEP reveals large-scale conformational changes 
of biomolecules by searching the path on an FES, where fast molecular 
motions are averaged out.
In both cases, the path is 
represented by a chain-of-replica, which evolves on the energy surface
so as to minimize the forces in the transverse direction.

-----------------------------------------------------------------------

**rpathmode** *MFEP/MEP*

  **Default: MFEP**

  Specify **MFEP** or **MEP** to invoke the MFEP or MEP search.

-----------------------------------------------------------------------

Minimum Free-Energy Path (MFEP) Search
======================================

The MFEP search is invoked by specifying **rpathmode = MFEP**.
The path search is carried out by the string method, which is a powerful 
sampling technique to find a path connecting two stable conformational 
states.  This method is widely used for investigating large-scale
conformational changes of biomolecules where time-scale of the
transitions is not reachable in brute-force simulations.

There are three major algorithms in the string method:
the mean forces string method :cite:`string_mean-forces`,
the on-the-fly string method :cite:`string_on-the-fly`,
and the string method of swarms of trajectories :cite:`string_swarms`.
Among these algorithms, the mean forces
string method is available in **ATDYN** and **SPDYN** :cite:`string_genesis`.

In the mean-forces string method, the pathway is represented by discretized 
points (called images) in the collective variable (CV) space. 
The current GENESIS supports distances, angles, dihedrals, Cartesian
coordinates, and principal components for CVs
(note that different types of CVs cannot be mixed in GENESIS. For example, 
users cannot mix distance and angle). 

In the calculation, each image is assigned to each replica, and 
a replica samples mean forces and an average metric
tensor around its own image by short MD simulation (ps to ns length) with
restraints.
The restraints are imposed using the image coordinates as their reference values. 
After the short simulation, each image is
evolved according to the mean force and metric tensor.
Then, smoothing and re-parametrization of the images are performed
and we go to the next cycle.

Image coordinates are written in rpath files (**rpathfile** keyword)
which user can specify in **[OUTPUT]** section.
This file provides the trajectory of the image coordinates.
Columns correspond to CVs and rows are time steps.
These values are written at the same timing as **dcdfile**
(specified by **crdout_period** in **[DYNAMICS]** section).

For the string method calculation, an initial pathway in the CV space
and atomistic coordinates around the pathway are required. For preparing these, 
targeted, or steered MD methods are recommended. 

-----------------------------------------------------------------------

**nreplica** *Integer*

  **Default: 1**

  Number of replicas (images) for representing the pathway.

**rpath_period** *Integer*

  **Default: 0**

  Time-step period during which the mean-forces acting on the
  images are evaluated. After evaluating the mean-forces,
  we update images according to the mean-forces, and go to the
  next cycle. If **rpath_period = 0**, images are not updated. This
  option is used for equilibration or umbrella sampling around the
  pathway. 

**delta** *Real*

  **Default: 0.0**

  Step-size for steepest descent update of images.

**smooth** *Real*

  **Default: 0.0**

  Smoothing parameter which controls the aggressiveness of the smoothing.
  Values from 0.0 to 0.1 are recommended,
  where "smooth = 0.0" means no-smoothing

**rest_function** *List of Integers*

  **Default: N/A**

  List of restraint function indices
  defined in **[RESTRAINTS]** section (see :ref:`restraints`).
  Specified restraints are defined as CVs, and *nreplica* images
  (replicas) are created, where a set of corresponding restraint
  reference values is assigned to each image.
  Force constants in **[RESTRAINTS]** are also used for evaluation
  of mean-forces. 

**fix_terminal** *YES / NO*

  **Default: NO**

  If **fix_terminal = YES** is specified, the two terminal images are
  always fixed and not updated. This is useful when the terminal images
  correspond to crystal structures and users do not want to move them. 

**use_restart** *YES / NO*

  **Default : YES**

  Restart file generated by the string method calculation includes the
  last snapshot of images. If **use_restart = YES** is specified,
  the reference values in **[RESTRAINTS]** will be overwritten by
  the values in the restart file. Note that force constants are not
  overwritten.


.. note::
   The following options are needed in the **[FITTING]** section when
   Cartesian coordinates are used for CVs.

   **fitting_method** *TR+ROT / XYTR+ZROT / NO*

     If this keyword is specified, rotational/translational elements are
     removed from the mean-force estimation by fitting instaneous
     structures to the reference coordinates given by **fitfile**. 

   **fitting_atom** *List of Integers*

     The user can specify the index of an atom group which are fitted  
     to the reference structure.
     Usually, the same atoms as CVs are selected.

Minimum Energy Path (MEP) Search
================================
The MEP search is available only in **ATDYN**. 

The calculation is invoked by specifying **rpathmode = MEP** :cite:`Yagi:2021`.
Cartesian coordinates of atoms selected via **mepatm_select_index** 
(denoted MEP atoms)  are employed for the path search. 
Currently, other coordinates can not be used as CVs.
Starting from a set of images along an initial path, e.g., a linear 
interpolation between the reactant and product,
the coordinates of the surrounding atoms are first energy minimized 
with MEP atoms held fixed. (The [MINIMIZE] section is thus
required in the input as well.) Then, the forces 
acting on MEP atoms are evaluated for each image, and the 
images are evolved to minimize the forces in the transverse 
direction by the string method :cite:`string_min_ene_path`.
The process is repeated either until the convergence 
threshold is met [variation in the energy (**tol_energy**) and the 
path length (**tol_path**)], or 
until the number of iterations reaches the maximum number (**ncycle**). 

Another search algorithm, the nudged elastic band (NEB) method 
:cite:`nudged-elastic-band` is also implemented, which differs from 
the string method in how the images evolve.  Note, however, that 
the NEB is still experimental at this moment.

-----------------------------------------------------------------------

**mepatm_select_index** *Integer*

  Index of a group of atoms which is treated as MEP atoms. 
  The index must be defined in **[SELECTION]** (see :ref:`selection`).

**ncycle** *Integer*

  **Default: 1000**

  Maximum number of cycle.

**nreplica** *Integer*

  **Default: 1**

  Number of replicas (images) for representing the pathway. 

.. note::
  If MPI processes are larger than **nreplica**, the MPI processes must 
  be a multiple of **nreplica**. For example, if **nreplica = 16**, 
  MPI processes must be 16, 32, 48, etc. If MPI processes are smaller 
  than **nreplica**, the MPI processes must be a divisor of **nreplica**.
  For example, the calculation with **nreplica = 16**  can be performed 
  using 1, 2, 4, and 8 MPI processes.

**eneout_period** *Integer*

  **Default : 10**

  Frequency of the output of the energy profile and path length to the 
  standard output.

**crdout_period** *Integer*

  **Default : 0**

  Frequency of coordinates outputs. Note that coordinate outputs are turned
  off for minimization (**crdout_period** in the [MINIMIZE] section).

**rstout_period** *Integer*
 
  **Default : 0**

  Frequency of restart file updates. Note that the updates are turned 
  off for minimization (**rstout_period** in the [MINIMIZE] section).

**tol_energy** *Real*

  **Default : 0.01** (unit : kcal/mol)

  Tolerence of convergence for the energy.

**tol_path** *Real*

  **Default : 0.01** (unit : :math:`\text{\AA}`)

  Tolerence of convergence for the path length.

**massweightcoord** *YES / NO*

  **Default : NO** 

  Use mass weighted Cartesian.

**method** *STRING/NEB*

  **Default: STRING**

  Choose the algorithm of a MEP search.

-----------------------------------------------------------------------

Options for String.

**delta** *Real*

  **Default: 0.001**

  Step-size for steepest descent update of images.

-----------------------------------------------------------------------

Options for NEB.

**k_spring** *Real*

  **Default: 10.0** :math:`{\rm{kcal/mol/{\text{\AA}}^{2}}}`

  Spring constant of the force that connects the images

**ncorrection** *Integer*

  **Default : 10**

  Number of corrections to build the inverse Hessian.

**lbfgs_bnd** *YES / NO*

  **Default : YES**

  Set a boundary to move atoms in each step of image update.

**lbfgs_bnd_qmonly** *YES / NO*

  **Default : NO**

  Set the boundary only to QM atoms.

**lbfgs_bnd_maxmove** *Real*

  **Default : 0.1** (unit : :math:`\text{\AA}`)

  The maximum size of move in each step.

-----------------------------------------------------------------------

Examples
========

Example of alanine-tripeptide with 16 replicas (images). Two dihedral
angles are specified as the collective variables. 
::
  
  [RPATH]
  nreplica          = 16
  rpath_period      = 1000
  delta             = 0.02
  smooth            = 0.0
  rest_function     = 1 2
  
  [SELECTION]
  group1        = atomindex:15
  group2        = atomindex:17
  group3        = atomindex:19
  group4        = atomindex:25
  group5        = atomindex:27
  
  [RESTRAINTS]
  nfunctions    = 2
  
  function1     = DIHED
  constant1     = 100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0 \
                  100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0
  reference1    = -40.0 -40.0 -40.0 -40.0 -40.0 -40.0 -40.0 -40.0 \
                  -40.0 -40.0 -40.0 -40.0 -40.0 -40.0 -40.0 -40.0
  select_index1 = 1 2 3 4  # PHI
  
  function2     = DIHED
  constant2     = 100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0 \
                  100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0
  reference2    = -45.0 -33.0 -21.0 -9.0 3.0 15.0 27.0 39.0 \
                   51.0  63.0  75.0 87.0 99.0 111.0 123.0 135.0
  select_index2 = 2 3 4 5  # PSI

Here is another example of Cartesian coordiante CVs for the same
alanine-tripeptide.
::
  
  [INPUT]
  ... skip ...
  rstfile = ../eq/{}.rst
  reffile = {}.pdb
  fitfile = fit.pdb
  
  [RPATH]
  nreplica          = 16
  rpath_period      = 1000
  delta             = 0.001
  smooth            = 0.00
  rest_function     = 1
  fix_terminal      = NO
  
  [FITTING]
  fitting_method = TR+ROT
  fitting_atom   = 1
  
  [SELECTION]
  group1        = ai:15 or ai:17 or ai:19 or ai:25 or ai:27
  
  [RESTRAINTS]
  nfunctions    = 1
  
  function1     = POSI
  constant1     = 10.0 10.0 10.0 10.0 \
                  10.0 10.0 10.0 10.0 \
                  10.0 10.0 10.0 10.0 \
                  10.0 10.0 10.0 10.0
  select_index1 = 1

Here is an example of a MEP search along with QM/MM
::

  [INPUT]
  topfile = toppar/top_all36_prot.rtf, ...
  parfile = toppar/par_all36_prot.prm, ...
  psffile = prot.psf               # protein structure file
  reffile = prot.pdb               # PDB file
  pdbfile = initial/initial{}.pdb  # initial path
  
  [OUTPUT]
  dcdfile = mep_{}.dcd        # coordinates
  logfile = mep_{}.log        # log files
  rstfile = mep_{}.rst        # restart file
  rpathfile  = mep_{}.rpath   # rpath file
  
  [ENERGY]
  forcefield       = CHARMM    # CHARMM force field
  ... skip ...
  
  [MINIMIZE]
  method              = LBFGS  # MIN using L-BFGS
  nsteps              = 100    # max. number of steps
  eneout_period       = 5      # energy output period
  fixatm_select_index = 2      # fix the outer layer
  macro               = yes    # macro/micro iteration

  [RPATH]
  rpathmode           = MEP    # MEP search
  method              = STRING # String method
  delta               = 0.0005 # step size
  ncycle              = 200    # max. number of cycle
  nreplica            = 16     # number of replica
  eneout_period       = 1      # frequency of the energy output
  crdout_period       = 1      # frequency of the coordinate output
  rstout_period       = 1      # frequency of the restart update
  fix_terminal        = no     # fix the terminal
  massWeightCoord     = no     # mass-weighted Cartesian
  mepatm_select_index = 1      # selection of the MEP atoms
  
  [BOUNDARY]
  type          = NOBC

  [QMMM]
  qmtyp              = gaussian
  qmatm_select_index = 1
  ... skip ...

  [SELECTION]
  group1 = sid:DHA or (sid:TIMA and (rno:95 or rno:165) and \
           not (an:CA | an:C | an:O | an:N | an:HN | an:HA))
  group2 = not (sid:DHA  or sid:DHA  around_res:6.0)


