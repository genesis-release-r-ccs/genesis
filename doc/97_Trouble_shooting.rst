.. highlight:: bash
.. _trouble:

=======================================================================
Trouble Shooting
=======================================================================

The followings are representative error messages that the users frequently encounter during the simulations.
We describe possible reasons for each error message, and provide suggestions to solve the problem.
 
**Compute_Shake> SHAKE algorithm failed to converge**

  This message indicates that the constraint of rigid bonds using the SHAKE algorithm (see :ref:`Constraints`) failed with error.
  In most cases, SHAKE errors originate from insufficient equilibration, bad initial structures, and/or inappropriate input parameters.
  We recommend the users to check the following points:

  * Reconsider the equilibration scheme. More moderate equilibration might be needed.
    For example, heating the system from 0 K, using a shorter timestep (e.g., 1.0 fs),
    and/or performing long energy minimization could be a possible solution.

  * Check the initial structure very carefully.
    One of the frequent mistakes in the initial structure modeling is "ring penetration" of covalent bonds, where one covalent bond is accidentally inserted into an aromatic ring.
    Check and solve such erroneous structures first, and then try the simulation again.

  * Some force field parameters are missing or wrong, which can easily cause unstable simulations.


**Check_Atom_Coord> Some atoms have large clashes**

  This message indicates that there is an atom pair whose distance is zero or close to zero.
  The indices of those atoms and the distance are displayed in a warning message: "WARNING: too short distance:".
  This situation is not allowed, especially in **SPDYN**, since it can cause a numerical error in the lookup table method.
  Check the initial structure first. Even if you cannot see such atomic clashes,
  there may be a clash between the atoms in the unit cell and image cells in the case of the periodic boundary condition.
  One of the automatic solutions is to specify "contact_check = YES" in the control file (see :ref:`Energy`).
  However, this cannot work well, if the distance is exactly zero.
  In such cases, the problem should be solved by the users themselves.
  For example, the users may have to slightly move the clashing atoms manually, or specify larger or smaller box size, or
  rebuild the initial structure more carefully.


**Setup_Processor_Number> Cannot define domains and cells. Smaller MPI processors, or shorter pairlistdist, or larger boxsize should be used**

  This message indicates that the total number of MPI processors used in your calculation is not appropriate for your system.
  The users had better understand relations between the system size and number of MPI processors.
  In **SPDYN**, the system is divided into several domains for parallel computation,
  where the number of domains must be equal to the number of MPI processors (see :ref:`available_programs`).
  In most cases, this message tells you that the system could not be divided into the specified number of domains.
  Although there are mainly three solutions to this problem, the first one is the most recommended way:

  * Use smaller number of MPI processors. If it can work, the previous number was too large to handle the system.

  * Use shorter pairlistdist. This treatment can make a domain size smaller, allowing to use a larger number of MPI processors.
    However, this is not recommended, if you are already using a recommended parameter set for
    switchdist, cutoffdist, and pairlistdist (e.g., 10, 12, and 13.5 :math:`\text{\AA}` in the CHARMM force field)

  * Build a larger initial structure by adding solvent molecules in the system,
    which may allow the users to divide the system into the desired number of domains.


**Update_Boundary_Pbc> too small boxsize/pairdist. larger boxsize or shorter pairdist should be used.**

  This message indicates that your system is too small to handle in the periodic boundary condition.
  In **ATDYN**, cell-linked list method is used to make non-bonded pairlists,
  where the cell size is determined to be close to and larger than the pairlist distance given in the control file.
  In addition, the total number of cells in x, y, and z dimensions must be at least three.
  **SPDYN** has a similar lower limitation in the available box size.
  Therefore, in order to solve this problem, the users may have to set a shorter pairlistdist,
  or build a larger system by adding much solvent molecules.


**Compute_Energy_Experimental_Restraint_Emfit> Gaussian kernel is extending outside the map box**

  This message indicates that the simulated densities were generated outside of the target density map.
  This error can frequently happen, when then atoms to be fitted are located near the edge of the target density map.

  * Create a larger density map by adding a sufficient margin to the map,
    which can be easily accomplished with the "voledit" tool in SITUS (https://situs.biomachina.org/).

  * Examine a normal MD simulation by turning off the EM biasing potential (emfit = NO).
    If the simulation is still unstable, there is likely an issue in the molecular mechanics calculation
    rather than the biasing potential calculation.
    In such cases, please check the initial structure carefully.
    There might be large clashes between some atoms, which can cause explosion of
    the target molecule, and push some atoms out of the density map.
    The problems to be solved are almost the same as those in the SHAKE errors (see above).


**Compute_Energy_Restraints_Pos> Positional restraint energy is too big**

  This message indicates that some atoms to be restrained are significantly deviated
  from the reference position, indicating that the restraint is not properly applied to such atoms.
  This situation is not allowed in **SPDYN**.

  * Use a larger force constant to keep their position near the reference.

  * Turn off the positional restraint for such atoms if it is not essential.
