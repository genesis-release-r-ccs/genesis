.. highlight:: bash
.. _selection:

=======================================================================
Selection section
=======================================================================


Atom selection
==============

This section is used to select atoms, and define them as a group.
The user can select atoms accorging to their name, index, residue number, segment name, and so on.
The selected group index is used in other sections.
For example, restraint potential can be applied on the group selected in this section,
and the force constant of the potential is specified in the **[RESTRAINTS]** section.
**[SELECTION]** section is also used in the GENESIS analysis tools to specify the atoms to be analyzed.

-----------------------------------------------------------------------

**group**:math:`\textbf{\textit{N}}` *expression*

  The user defines selected atoms as "group1", "group2", ..., and :math:`\text{group}N`.
  Here, *N* must be a positive integer (:math:`N \geq 1`).
  The user selects atoms by using keywords and operators with a certain syntax (see table below).
  Note that in the table mname (or moleculename, molname) is a molecule name defined by **mol_name**.

**mole_name**:math:`\textbf{\textit{N}}` *molecule* *starting-residue* *ending-residue*

  The user defines a molecule by specifying its segment name, first and last residue numbers,
  and residue name. :math:`N` must be a positive integer (:math:`N \geq 1`).
  The syntax for the residue selection is as follows:

  ``[segment name]:[residue number]:[residue name]``

  For details, see the example below.


Table. Available keywords and operators in ``group``.

======================== ============================== =============== ==========================
expression               meaning                        example         other available expression  
======================== ============================== =============== ==========================
an:*name*                atom name                      an:CA           atomname, atom_name 
ai:*number[-[number]]*   atom index                     ai:1-5          atomindex, atomidx 
atno:*number[-[number]]* atom number                    atno:6          atomno 
rnam:*name*              residue name                   rnam:GLY        residuename, resname 
rno:*number[-[number]]*  residue number                 rno:1-5         residueno, resno
mname:*name*             molecule name                  mname:molA      moleculename, molname 
segid:*ID*               segment index                  segid:PROA      segmentid, sid 
hydrogen                 hydrogen atoms                                 hydrogenatom
heavy                    heavy atoms                                    heavyatom
all                      all atoms                                      ``*``
and                      conjunction                                    &
or                       logical add                                    \|
not                      negation                                       !
 ()                      assemble                        
======================== ============================== =============== ==========================

Table. Available keywords and operators in ``group`` (continued).

======================== ====================================  =============== ==========================
expression               meaning                               example         other available expression  
======================== ====================================  =============== ==========================
*X* around: *r*          atoms around *r* Angstrom of *X*      see below       around_atoms
*X* around_res: *r*      residues around *r* Angstrom of *X*   see below       around_residues
*X* around_mol: *r*      molecules around *r* Angstrom of *X*  see below       around_molecules
======================== ====================================  =============== ==========================

.. Note::

  ``ai`` and ``atno`` are slightly different. ``ai`` indicates the atom index 
  which is sequentially re-numbered over all atoms in the system.
  On the other hand, ``atno`` is the index of atoms in the PDB file.
  Atom index in PDB file (column 2) does not always start from 1, nor is numbered sequentially.
  In such cases, ``atno`` is useful to select atoms, although it is a very rare case.

.. Note::

  Atoms that are within a distance of a given atom (X) can be selected by ``around``.
  Note that the coordinates in ``reffile`` is used to judge the distance. If ``reffile``
  is not present, those in input files (``pdbfile``, ``crdfile``, etc.) are used instead.
  Coordinates in ``rstfile`` are never used.


Examples
==============

Select atoms based on their atom name, residue name, or residue number:
:: 
  [SELECTION]
  group1     = resno:1-60 and an:CA
  group2     = (segid:PROA and not hydrogen) | an:CA
  mole_name1 = molA PROA:1:TYR PROA:5:MET
  group3     = mname:molA and (an:CA or an:C or an:O or an:N)

Select atoms around an atom X. In the following examples, X = atom number 100.
:: 
  [SELECTION]
  group1     = atno:100 around:10.0
  group2     = atno:100 around_res:10.0
  group3     = atno:100 around_mol:10.0
  group4     = atno:100 around_mol:10.0 or atno:100

In group1, atoms around 10.0  :math:`\text{\AA}` of X are selected. Group 2 selects residues around
10.0 :math:`\text{\AA}` of X, i.e., if the distance between X and any one of atoms in a residue is
less than 10.0 :math:`\text{\AA}`, all atoms of the residue are selected. Group 3 is the same as
group 2, but for a molecule. Note that these commands do NOT select X itself.
In order to include X in the selection, add "or atno:100", as in group 4.

------------------------------------------------------------------

Select atoms around multiple atoms.
:: 
    
  [SELECTION]
  group1     = atno:100-101 around:10.0
  group2     = (sid:PROT around_res:10.0) and rnam:TIP3
  group3     = (rno:1 around:10.0) or rno:1

Group 1 selects atoms around 10.0 :math:`\text{\AA}` of atom 100 *or* 101. Note that it is NOT
"100 *and* 101" nor a center of 100 and 101. Group 2 is an example to select water
molecules around a protein (segname PROT). Group 3 selects not only the atoms around
residue1 but also the atoms of residue1.
