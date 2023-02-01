!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_enefunc_charmm_mod
!> @brief   define potential energy functions
!! @authors Jaewoon Jung (JJ), Takeshi Imai (TI), Chigusa Kobayashi (CK),
!!          Takaharu Mori (TM), Yuji Sugita (YS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_enefunc_charmm_mod

  use sp_enefunc_localres_mod
  use sp_enefunc_restraints_mod
  use sp_enefunc_fit_mod
  use sp_enefunc_table_mod
  use sp_energy_mod
  use sp_restraints_str_mod
  use sp_constraints_str_mod
  use sp_enefunc_str_mod
  use sp_energy_str_mod
  use sp_domain_str_mod
  use dihedral_libs_mod
  use molecules_str_mod
  use fileio_par_mod
  use fileio_localres_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: define_enefunc_charmm
  private :: setup_enefunc_bond
  private :: setup_enefunc_bond_constraint
  private :: setup_enefunc_angl
  private :: setup_enefunc_angl_constraint
  private :: setup_enefunc_dihe
  private :: setup_enefunc_impr
  private :: setup_enefunc_cmap
  private :: setup_enefunc_nonb
  public  :: count_nonb_excl
  public  :: check_pbc
  ! FEP
  public  :: count_nonb_excl_fep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc_charmm
  !> @brief        a driver subroutine for defining potential energy functions
  !! @authors      YS, TI, JJ, CK, HO
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    par         : CHARMM PAR information
  !! @param[in]    localres    : local restraint information
  !! @param[in]    molecule    : molecule information
  !! @param[inout] constraints : constraints information
  !! @param[in]    restraints  : restraints information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc_charmm(ene_info, par, localres, molecule, &
                                   constraints, restraints, domain, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_par),             intent(in)    :: par
    type(s_localres),        intent(in)    :: localres
    type(s_molecule),        intent(in)    :: molecule
    type(s_constraints),     intent(inout) :: constraints
    type(s_restraints),      intent(in)    :: restraints
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: ncel, ncelb


    ! base
    !
    ncel  = domain%num_cell_local
    ncelb = domain%num_cell_local + domain%num_cell_boundary

    call alloc_enefunc(enefunc, EneFuncBase, ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncBond, ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncAngl, ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncDihe, ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncImpr, ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncBondCell, ncel, ncelb)
    if (domain%fep_use) then
      ! FEP
      call alloc_enefunc(enefunc, EneFuncBondCell_FEP, ncel, ncelb)
      call alloc_enefunc(enefunc, EneFuncFEPBonded, ncel, ncel)
    end if

    if (.not. constraints%rigid_bond) then

      ! bond
      !
      call setup_enefunc_bond(par, molecule, constraints, domain, enefunc)

      ! angle
      !
      call setup_enefunc_angl(par, molecule, domain, enefunc)

    else

      ! bond
      !
      call setup_enefunc_bond_constraint( &
                              par, molecule, domain, constraints, enefunc)

      ! angle
      !
      call setup_enefunc_angl_constraint( &
                              par, molecule, domain, constraints, enefunc)

    end if

    ! dihedral
    !
    call setup_enefunc_dihe(par, molecule, domain, enefunc)

    ! improper
    !
    call setup_enefunc_impr(par, molecule, domain, enefunc)

    ! cmap
    !
    call setup_enefunc_cmap(ene_info, par, molecule, domain, enefunc)

    ! nonbonded
    !
    call setup_enefunc_nonb(par, molecule, constraints, domain, enefunc)

    ! lookup table
    !
    if (ene_info%table) &
    call setup_enefunc_table(ene_info, enefunc)

    ! restraint
    !
    call setup_enefunc_restraints(molecule, restraints, domain, enefunc)

    call setup_enefunc_localres(localres, domain, enefunc)

    ! write summary of energy function
    !
    if (main_rank) then
      write(MsgOut,'(A)') &
           'Define_Enefunc_Charmm> Number of Interactions in Each Term'
      write(MsgOut,'(A20,I10,A20,I10)')                  &
           '  bond_ene        = ', enefunc%num_bond_all, &
           '  angle_ene       = ', enefunc%num_angl_all
      write(MsgOut,'(A20,I10,A20,I10)')                  &
           '  torsion_ene     = ', enefunc%num_dihe_all, &
           '  improper_ene    = ', enefunc%num_impr_all
      write(MsgOut,'(A20,I10)')                          &
           '  cmap_ene        = ', enefunc%num_cmap_all
      write(MsgOut,'(A20,I10,A20,I10)')                  &
           '  nb_exclusions   = ', enefunc%num_excl_all, &
           '  nb14_calc       = ', enefunc%num_nb14_all
      if (domain%fep_use) then
        write(MsgOut,'(A20,I10)')                        &
           '  nb14_calc_fep   = ', enefunc%num_nb14_all_fep
      end if
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine define_enefunc_charmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_bond
  !> @brief        define BOND term for each cell in potential energy function
  !! @authors      YS, JJ, TM, HO
  !! @param[in]    par      : CHARMM PAR information
  !! @param[in]    molecule : molecule information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond(par, molecule, constraints, domain, enefunc)

    ! formal arguments
    type(s_par),             intent(in)    :: par
    type(s_molecule),target, intent(in)    :: molecule
    type(s_constraints),     intent(in)    :: constraints
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    real(wp)                 :: m1, m2
    integer                  :: dupl, ioffset
    integer                  :: i, j, icel_local
    integer                  :: icel1, icel2
    integer                  :: nbond, nbond_p, found, found1
    integer                  :: i1, i2, i3, ia, ib, pbc_int
    integer                  :: wat_bonds
    character(6)             :: ci1, ci2, ci3, ri1, ri2
    real(wp)                 :: cwork(3,2), dij(3)

    real(wp),        pointer :: force(:,:), dist(:,:), coord(:,:)
    real(wp),        pointer :: box_size(:)
    integer,         pointer :: nwater(:), bond(:), list(:,:,:)
    integer,         pointer :: ncel
    integer(int2),   pointer :: cell_pair(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: mol_bond_list(:,:)
    character(6),    pointer :: mol_cls_name(:), mol_res_name(:)
    integer,         pointer :: bond_pbc(:,:)

    ! FEP
    integer                  :: nbond_fep, k, fg1, fg2
    logical                  :: flag_singleB
    integer,         pointer :: bond_singleB(:,:)

    mol_bond_list => molecule%bond_list
    mol_cls_name  => molecule%atom_cls_name
    mol_res_name  => molecule%residue_name
    coord         => molecule%atom_coord

    ncel          => domain%num_cell_local
    cell_pair     => domain%cell_pair
    id_g2l        => domain%id_g2l
    nwater        => domain%num_water
    box_size      => domain%system_size

    bond          => enefunc%num_bond
    list          => enefunc%bond_list
    force         => enefunc%bond_force_const
    dist          => enefunc%bond_dist_min
    bond_pbc      => enefunc%bond_pbc

    nbond         = molecule%num_bonds
    nbond_p       = par%num_bonds

    ! FEP
    if (domain%fep_use) then
      bond_singleB => enefunc%bond_singleB
    end if

    do dupl = 1, domain%num_duplicate

      ioffset = (dupl-1) * enefunc%table%num_all

      do i = 1, nbond

        i1  = molecule%bond_list(1,i)
        i2  = molecule%bond_list(2,i)
        ci1 = mol_cls_name(i1)
        ci2 = mol_cls_name(i2)
        ri1 = molecule%residue_name(i1)
        ri2 = molecule%residue_name(i2)

        if (ri1(1:3) .ne. 'TIP' .and. ri1(1:3) .ne. 'WAT' .and. &
            ri1(1:3) .ne. 'SOL' .and. ri2(1:3) .ne. 'TIP' .and. &
            ri2(1:3) .ne. 'WAT' .and. ri2(1:3) .ne. 'SOL') then

          ! FEP
          if (domain%fep_use) then
            flag_singleB = .false.
            fg1 = molecule%fepgrp(i1)
            fg2 = molecule%fepgrp(i2)
            if (molecule%fepgrp_bond(fg1,fg2) == 0) then
              ! FEP: If the bond are not set to any group of FEP, exclude this bond.
              cycle
            else if (molecule%fepgrp_bond(fg1,fg2) == 2) then
              ! FEP: If 2-2 or 2-5, flag is set to true.
              flag_singleB = .true.
            end if
            ! FEP: If i1 is singleB atom, its ID is replaced with singleA's ID.
            if (fg1==2) then
              do k = 1, molecule%num_atoms_fep(2)
                if (molecule%id_singleB(k) == i1) then
                  exit
                end if
              end do
              i1 = molecule%id_singleA(k)
            end if
            ! FEP: If i2 is singleB atom, its ID is replaced with singleA's ID.
            if (fg2==2) then
              do k = 1, molecule%num_atoms_fep(2)
                if (molecule%id_singleB(k) == i2) then
                  exit
                end if
              end do
              i2 = molecule%id_singleA(k)
            end if
          end if

          ia  = i1 + ioffset
          ib  = i2 + ioffset

          icel1 = id_g2l(1,ia)
          icel2 = id_g2l(1,ib)

          ! Check if it is in my domain
          !
          if (icel1 /= 0 .and. icel2 /= 0) then

            icel_local = cell_pair(icel1,icel2)

            if (icel_local > 0 .and. icel_local <= ncel) then

              do j = 1, nbond_p
                if ((ci1 == par%bond_atom_cls(1,j) .and.  &
                     ci2 == par%bond_atom_cls(2,j)) .or.  &
                    (ci1 == par%bond_atom_cls(2,j) .and.  &
                     ci2 == par%bond_atom_cls(1,j))) then
  
                  bond (icel_local) = bond(icel_local) + 1
                  list (1,bond(icel_local),icel_local) = ia
                  list (2,bond(icel_local),icel_local) = ib
                  force(bond(icel_local),icel_local)   = par%bond_force_const(j)
                  dist (bond(icel_local),icel_local)   = par%bond_dist_min(j)
                  cwork(1:3,1) = coord(1:3,i1)
                  cwork(1:3,2) = coord(1:3,i2)
                  dij(1:3) = cwork(1:3,1) - cwork(1:3,2)
                  call check_pbc(box_size, dij, pbc_int)
                  bond_pbc(bond(icel_local),icel_local) = pbc_int

                  if (domain%fep_use) then
                    ! FEP: flag for singleB bond
                    if (flag_singleB) &
                      bond_singleB (bond(icel_local),icel_local) = 1
                  end if

                  exit
  
                end if
              end do

              if (j == nbond_p + 1) &
                write(MsgOut,*) &
                  'Setup_Enefunc_Bond> not found BOND: [', &
                  ci1, ']-[', ci2, '] in parameter file. (ERROR)'
  
            end if
  
          end if

        end if
 
      end do
    end do

    ! bond/angel from water
    !
    wat_bonds = 0
    do i = 1, nbond
      i1 = mol_bond_list(1,i)
      i2 = mol_bond_list(2,i)
      ri1 = mol_res_name(i1)
      ri2 = mol_res_name(i2)
      if (ri1(1:3) .eq. 'TIP' .or. ri1(1:3) .eq. 'WAT' .or. &
          ri1(1:3) .eq. 'SOL') then
        wat_bonds = wat_bonds+1
      end if
    end do

    if (.not.constraints%fast_water .and. enefunc%table%num_water > 0) then

      i1 = enefunc%table%water_list(1,1)
      i2 = enefunc%table%water_list(2,1)
      i3 = enefunc%table%water_list(3,1)
      ci1 = molecule%atom_cls_name(i1)
      ci2 = molecule%atom_cls_name(i2)
      ci3 = molecule%atom_cls_name(i3)
      do j = 1, nbond_p
        if ((ci1 == par%bond_atom_cls(1, j) .and.  &
             ci2 == par%bond_atom_cls(2, j)) .or.  &
            (ci1 == par%bond_atom_cls(2, j) .and.  &
             ci2 == par%bond_atom_cls(1, j))) then
          enefunc%table%OH_bond = par%bond_dist_min(j)
          enefunc%table%OH_force = par%bond_force_const(j)
          enefunc%table%water_bond_calc_OH = .true.
          exit
        end if
      end do
      do j = 1, nbond_p
        if ((ci2 == par%bond_atom_cls(1, j) .and.  &
             ci3 == par%bond_atom_cls(2, j)) .or.  &
            (ci2 == par%bond_atom_cls(2, j) .and.  &
             ci3 == par%bond_atom_cls(1, j))) then
          enefunc%table%HH_bond = par%bond_dist_min(j)
          enefunc%table%HH_force = par%bond_force_const(j)
          enefunc%table%water_bond_calc_HH = .true.
          exit
        end if
      end do
      if (enefunc%table%water_bond_calc_OH .or. &
          enefunc%table%water_bond_calc_HH)     &
          enefunc%table%water_bond_calc = .true.

    end if

    found  = 0
    do i = 1, ncel
      found = found + bond(i)
      if (bond(i) > MaxBond) &
        call error_msg('Setup_Enefunc_Bond> Too many bonds.')
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_bond_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_bond_all = found
#endif
    enefunc%num_bond_all = enefunc%num_bond_all + wat_bonds*domain%num_duplicate

    if (domain%fep_use) then

      nbond_fep = 0
      do i = 1, 5
        nbond_fep = nbond_fep + molecule%num_bonds_fep(i)
      end do
      if (enefunc%num_bond_all /= nbond_fep*domain%num_duplicate) &
        call error_msg('Setup_Enefunc_Bond> Some bond paremeters are missing.')

    else

      if (enefunc%num_bond_all /= nbond*domain%num_duplicate) &
        call error_msg('Setup_Enefunc_Bond> Some bond paremeters are missing.')

    end if

    return

  end subroutine setup_enefunc_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_bond_constraint
  !> @brief        define BOND term between heavy atoms
  !! @authors      JJ, HO
  !! @param[in]    par         : CHARMM PAR information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond_constraint(par, molecule, domain, &
                                           constraints, enefunc)

    ! formal arguments
    type(s_par),                 intent(in)    :: par
    type(s_molecule),    target, intent(in)    :: molecule
    type(s_domain),      target, intent(in)    :: domain
    type(s_constraints), target, intent(inout) :: constraints
    type(s_enefunc),     target, intent(inout) :: enefunc

    ! local variable
    integer                      :: dupl, ioffset
    integer                      :: i, j, k, m, ih, icel_local, connect
    integer                      :: ih1, ih2, icel1, icel2, icel
    integer                      :: i1, i2, ia, ib, pbc_int
    integer                      :: nbond, nbond_p, nbond_a, nbond_c
    integer                      :: wat_bonds, tmp_mole_no, mole_no, wat_found
    character(6)                 :: ci1, ci2
    character(6)                 :: ri1, ri2
    logical                      :: mi1, mi2
    logical                      :: cl1, cl2
    real(wp)                     :: cwork(3,2), dij(3)

    real(wp),            pointer :: force(:,:), dist(:,:), coord(:,:)
    real(wp),            pointer :: box_size(:)
    real(wip),           pointer :: HGr_bond_dist(:,:,:,:)
    integer,             pointer :: bond(:), list(:,:,:), num_water, ncel
    integer,             pointer :: id_l2g(:,:), id_l2g_sol(:,:)
    integer(int2),       pointer :: cell_pair(:,:)
    integer(int2),       pointer :: id_g2l(:,:)
    integer,             pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)
    integer,             pointer :: mol_bond_list(:,:), mol_no(:)
    character(6),        pointer :: mol_cls_name(:), mol_res_name(:)
    logical,             pointer :: mol_light_name(:), mol_light_mass(:)
    integer,             pointer :: bond_pbc(:,:)

    ! FEP
    integer                      :: l, nbond_fep, fg1, fg2
    logical                      :: flag_singleB
    integer,             pointer :: bond_singleB(:,:)

    mol_bond_list  => molecule%bond_list
    mol_no         => molecule%molecule_no
    mol_cls_name   => molecule%atom_cls_name
    mol_res_name   => molecule%residue_name
    mol_light_name => molecule%light_atom_name
    mol_light_mass => molecule%light_atom_mass
    coord          => molecule%atom_coord

    ncel           => domain%num_cell_local
    cell_pair      => domain%cell_pair
    id_g2l         => domain%id_g2l
    id_l2g         => domain%id_l2g
    id_l2g_sol     => domain%id_l2g_solute
    box_size       => domain%system_size

    HGr_local      => constraints%HGr_local
    HGr_bond_list  => constraints%HGr_bond_list
    HGr_bond_dist  => constraints%HGr_bond_dist

    bond           => enefunc%num_bond
    list           => enefunc%bond_list
    force          => enefunc%bond_force_const
    dist           => enefunc%bond_dist_min
    num_water      => enefunc%table%num_water
    bond_pbc       => enefunc%bond_pbc

    nbond          = molecule%num_bonds
    nbond_p        = par%num_bonds

    connect        =  constraints%connect

    nbond_a        = 0
    nbond_c        = 0

    ! FEP
    if (domain%fep_use) then
      bond_singleB => enefunc%bond_singleB
    end if

    do dupl = 1, domain%num_duplicate

      ioffset = (dupl-1) * enefunc%table%num_all

      do i = 1, nbond

        i1  = mol_bond_list(1,i)
        i2  = mol_bond_list(2,i)

        ri1 = mol_res_name(i1)
        ri2 = mol_res_name(i2)

        if (ri1(1:3) .ne. 'TIP' .and. ri1(1:3) .ne. 'WAT' .and. &
            ri1(1:3) .ne. 'SOL' .and. ri2(1:3) .ne. 'TIP' .and. &
            ri2(1:3) .ne. 'WAT' .and. ri2(1:3) .ne. 'SOL') then

          ci1 = mol_cls_name(i1)
          ci2 = mol_cls_name(i2)
          mi1 = mol_light_mass(i1)
          mi2 = mol_light_mass(i2)
          cl1 = mol_light_name(i1)
          cl2 = mol_light_name(i2)
        
          if (constraints%hydrogen_type == ConstraintAtomMass) then
            cl1 = mi1 
            cl2 = mi2 
          else if (constraints%hydrogen_type == ConstraintAtomBoth) then
            cl1 = (cl1 .or. mi1) 
            cl2 = (cl2 .or. mi2) 
          end if
  
          ! FEP
          if (domain%fep_use) then
            flag_singleB = .false.
            fg1 = molecule%fepgrp(i1)
            fg2 = molecule%fepgrp(i2)
            if (molecule%fepgrp_bond(fg1,fg2) == 0) then
              ! FEP: If the bond are not set to any group of FEP, exclude this bond.
              cycle
            else if (molecule%fepgrp_bond(fg1,fg2) == 2) then
              ! FEP: If 2-2 or 2-5, flag is set to true.
              flag_singleB = .true.
            end if
            ! FEP: If i1 is singleB atom, its ID is replaced with the ID corresponding to singleA.
            if (fg1 == 2) then
              do k = 1, molecule%num_atoms_fep(2)
                if (molecule%id_singleB(k) == i1) then
                  exit
                end if
              end do
              i1 = molecule%id_singleA(k)
            end if
            ! FEP: If i2 is singleB atom, its ID is replaced with the ID corresponding to singleA.
            if (fg2 == 2) then
              do k = 1, molecule%num_atoms_fep(2)
                if (molecule%id_singleB(k) == i2) then
                  exit
                end if
              end do
              i2 = molecule%id_singleA(k)
            end if
          end if

          ia  = i1 + ioffset
          ib  = i2 + ioffset
  
          if (.not. (cl1 .or.  cl2)) then
  
            icel1 = id_g2l(1,ia)
            icel2 = id_g2l(1,ib)
  
            ! Check if it is in my domain
            !
            if (icel1 /= 0 .and. icel2 /= 0) then
    
              icel_local = cell_pair(icel1,icel2)
    
              if (icel_local > 0 .and. icel_local <= ncel) then
  
                do j = 1, nbond_p
                  if ((ci1 == par%bond_atom_cls(1, j) .and.  &
                       ci2 == par%bond_atom_cls(2, j)) .or.  &
                      (ci1 == par%bond_atom_cls(2, j) .and.  &
                       ci2 == par%bond_atom_cls(1, j))) then
    
                    nbond_a = nbond_a + 1
                    bond (icel_local) = bond(icel_local) + 1
                    list (1,bond(icel_local),icel_local) = ia
                    list (2,bond(icel_local),icel_local) = ib
                    force(bond(icel_local),icel_local) = par%bond_force_const(j)
                    dist (bond(icel_local),icel_local) = par%bond_dist_min(j)
                    cwork(1:3,1) = coord(1:3,i1)
                    cwork(1:3,2) = coord(1:3,i2)
                    dij(1:3) = cwork(1:3,1) - cwork(1:3,2)
                    call check_pbc(box_size, dij, pbc_int)
                    bond_pbc(bond(icel_local),icel_local) = pbc_int

                    if (domain%fep_use) then
                      ! FEP: flag for singleB bond
                      if (flag_singleB) &
                        bond_singleB (bond(icel_local),icel_local) = 1
                    end if

                    exit
    
                  end if
                end do
    
                if (j == nbond_p + 1) &
                  write(MsgOut,*) &
                    'Setup_Enefunc_Bond_Constraint> not found BOND: [', &
                    ci1, ']-[', ci2, '] in parameter file. (ERROR)'
    
              end if
    
            end if
    
          else

            if (domain%fep_use) then
              ! FEP: skip singleB atoms if the bond connects with hydrogen.
              if (flag_singleB) cycle
            end if
   
            icel1 = id_g2l(1,ia)
            icel2 = id_g2l(1,ib)
  
            if (icel1 /= 0 .and. icel2 /= 0) then
  
              icel = cell_pair(icel1,icel2)
              if (icel > 0 .and. icel <= ncel) then
  
                do j = 1, connect
  
                  do k = 1, HGr_local(j,icel)
                    ih1 = id_l2g_sol(HGr_bond_list(1,k,j,icel),icel)
                    do ih = 1, j
                      ih2 = id_l2g_sol(HGr_bond_list(ih+1,k,j,icel),icel)
    
                      if (ih1 == ia .and. ih2 == ib .or. &
                          ih2 == ia .and. ih1 == ib) then
                      
                        do m = 1, nbond_p
                          if ((ci1 == par%bond_atom_cls(1, m) .and.  &
                               ci2 == par%bond_atom_cls(2, m)) .or.  &
                              (ci1 == par%bond_atom_cls(2, m) .and.  &
                               ci2 == par%bond_atom_cls(1, m))) then
      
                            nbond_c = nbond_c + 1
                            HGr_bond_dist(ih+1,k,j,icel) = par%bond_dist_min(m)
                            exit
                          end if
                        end do
    
                    end if
    
                    end do
                  end do
    
                end do
              end if
            end if
    
          end if
 
        end if
 
      end do
    end do

    ! for water molecule
    !
    wat_bonds = 0
    wat_found = 0

    if (constraints%fast_water) then

      tmp_mole_no = -1

      do i = 1, nbond

        i1 = mol_bond_list(1,i)
        i2 = mol_bond_list(2,i)
        ri1 = mol_res_name(i1)
        ri2 = mol_res_name(i2)

        if (ri1(1:3) .eq. 'TIP' .or. ri1(1:3) .eq. 'WAT' .or. &
            ri1(1:3) .eq. 'SOL') then

          wat_bonds = wat_bonds+1
          mole_no = mol_no(i1)

          if (mole_no /= tmp_mole_no) then
            wat_found = wat_found +1
            tmp_mole_no = mole_no
          end if

        end if

      end do

      if (wat_found /= num_water) &
        call error_msg( &
          'Setup_Enefunc_Bond_Constraint> # of water is incorrect')

    end if

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(nbond_a, enefunc%num_bond_all,  1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)

    call mpi_allreduce(nbond_c, constraints%num_bonds, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_bond_all  = nbond_a
    constraints%num_bonds = nbond_c
#endif

    if (domain%fep_use) then

      nbond_fep = 0
      do i = 1, 5
        nbond_fep = nbond_fep + molecule%num_bonds_fep(i)
      end do
!        if (enefunc%num_bond_all /= (nbond_fep*domain%num_duplicate &
!                                    -2*constraints%num_bonds &
!                                    -wat_bonds*domain%num_duplicate)) then
!          call error_msg( &
!            'Setup_Enefunc_Bond_Constraint> Some bond paremeters are missing.')
!        end if

    else

      if (enefunc%num_bond_all /= (nbond*domain%num_duplicate &
                                     -constraints%num_bonds     &
                                     -wat_bonds*domain%num_duplicate)) then
        call error_msg( &
          'Setup_Enefunc_Bond_Constraint> Some bond paremeters are missing.')
      end if



    end if

    return

  end subroutine setup_enefunc_bond_constraint

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_angl
  !> @brief        define ANGLE term for each cell in potential energy function
  !! @authors      YS, JJ, TM, HO
  !! @param[in]    par      : CHARMM PAR information
  !! @param[in]    molecule : molecule information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_angl(par, molecule, domain, enefunc)

    ! formal arguments
    type(s_par),     target, intent(in)    :: par
    type(s_molecule),target, intent(in)    :: molecule
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    logical                  :: angle_calc
    integer                  :: dupl, ioffset
    integer                  :: i, j, icel_local
    integer                  :: icel1, icel2
    integer                  :: nangl, nangl_p, found
    integer                  :: list(3), lists(3), pbc_int
    character(6)             :: ci1, ci2, ci3, ri1, ri2, ri3
    real(wp)                 :: cwork(3,3), dij(3)

    real(wp),        pointer :: force(:,:), theta(:,:), coord(:,:)
    real(wp),        pointer :: box_size(:)
    real(wp),        pointer :: ubforce(:,:), ubrmin(:,:)
    integer,         pointer :: angle(:), alist(:,:,:)
    integer,         pointer :: ncel
    integer(int2),   pointer :: cell_pair(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: mol_angl_list(:,:)
    integer,         pointer :: nwater(:)
    character(6),    pointer :: mol_cls_name(:), mol_res_name(:)
    integer,         pointer :: angl_pbc(:,:,:)

    ! FEP
    integer                  :: nangl_fep, k, fg(3)
    logical                  :: flag_singleB
    integer,         pointer :: angl_singleB(:,:)


    mol_angl_list => molecule%angl_list
    mol_cls_name  => molecule%atom_cls_name
    mol_res_name  => molecule%residue_name
    coord         => molecule%atom_coord

    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l
    nwater    => domain%num_water
    box_size  => domain%system_size

    angle     => enefunc%num_angle
    alist     => enefunc%angle_list
    force     => enefunc%angle_force_const
    theta     => enefunc%angle_theta_min
    ubforce   => enefunc%urey_force_const
    ubrmin    => enefunc%urey_rmin
    angl_pbc  => enefunc%angle_pbc

    nangl     = molecule%num_angles
    nangl_p   = par%num_angles

    ! FEP
    if (domain%fep_use) then
      angl_singleB => enefunc%angl_singleB
    end if

    do dupl = 1, domain%num_duplicate

      ioffset = (dupl-1)*enefunc%table%num_all

      do i = 1, nangl

        list(1:3) = mol_angl_list(1:3,i)
        ci1 = mol_cls_name(list(1))
        ci2 = mol_cls_name(list(2))
        ci3 = mol_cls_name(list(3))
        ri1 = mol_res_name(list(1))
        ri2 = mol_res_name(list(2))
        ri3 = mol_res_name(list(3))

        if (ri1(1:3) .ne. 'TIP' .and. ri1(1:3) .ne. 'WAT' .and. &
            ri1(1:3) .ne. 'SOL' .and. ri3(1:3) .ne. 'TIP' .and. &
            ri3(1:3) .ne. 'WAT' .and. ri3(1:3) .ne. 'SOL') then

          ! FEP
          if (domain%fep_use) then
            flag_singleB = .false.
            do j = 1, 3
              fg(j) = molecule%fepgrp(list(j))
            end do
            if (molecule%fepgrp_angl(fg(1),fg(2),fg(3)) == 0) then
              ! FEP: If the angle are not set to any group of FEP, exclude this angle.
              cycle
            else if (molecule%fepgrp_angl(fg(1),fg(2),fg(3)) == 2) then
              ! FEP: If 2-2 and 2-5
              flag_singleB = .true.
            end if
            do j = 1, 3
              ! FEP: If list(j) is singleB atom, its ID is replaced with the ID corresonding to singleA.
              if (fg(j) == 2) then
                do k = 1, molecule%num_atoms_fep(2)
                  if (molecule%id_singleB(k) == list(j)) then
                    exit
                  end if
                end do
                list(j) = molecule%id_singleA(k)
              end if
            end do
          end if

          lists(1:3) = list(1:3) + ioffset

          icel1 = id_g2l(1,lists(1))
          icel2 = id_g2l(1,lists(3))

          if (icel1 /= 0 .and. icel2 /= 0) then

            icel_local = cell_pair(icel1,icel2)

            if (icel_local >= 1 .and. icel_local <= ncel) then

              do j = 1, nangl_p
                if ((ci1 == par%angl_atom_cls(1,j) .and. &
                     ci2 == par%angl_atom_cls(2,j) .and. &
                     ci3 == par%angl_atom_cls(3,j)) .or. &
                    (ci1 == par%angl_atom_cls(3,j) .and. &
                     ci2 == par%angl_atom_cls(2,j) .and. &
                     ci3 == par%angl_atom_cls(1,j))) then
  
                  angle(icel_local) = angle(icel_local) + 1
                  alist(1:3,angle(icel_local),icel_local) = lists(1:3)
    
                  force(angle(icel_local),icel_local)   = &
                       par%angl_force_const(j)
                  theta(angle(icel_local),icel_local)   = &
                       par%angl_theta_min(j)*RAD
                  ubforce(angle(icel_local),icel_local) = &
                       par%angl_ub_force_const(j)
                  ubrmin(angle(icel_local),icel_local)  = &
                     par%angl_ub_rmin(j)
                  cwork(1:3,1) = coord(1:3,list(1))
                  cwork(1:3,2) = coord(1:3,list(2))
                  cwork(1:3,3) = coord(1:3,list(3))
                  dij(1:3) = cwork(1:3,1) - cwork(1:3,2)
                  call check_pbc(box_size, dij, pbc_int)
                  angl_pbc(1,angle(icel_local),icel_local) = pbc_int
                  dij(1:3) = cwork(1:3,3) - cwork(1:3,2)
                  call check_pbc(box_size, dij, pbc_int)
                  angl_pbc(2,angle(icel_local),icel_local) = pbc_int
                  dij(1:3) = cwork(1:3,1) - cwork(1:3,3)
                  call check_pbc(box_size, dij, pbc_int)
                  angl_pbc(3,angle(icel_local),icel_local) = pbc_int

                  if (domain%fep_use) then
                    ! FEP: flag for singleB angle
                    if (flag_singleB) &
                      angl_singleB (angle(icel_local),icel_local) = 1
                  end if

                  exit
  
                end if
              end do

              if (j == nangl_p + 1) &
                write(MsgOut,*) &
                  'Setup_Enefunc_Angl> not found ANGL: [',&
                  ci1, ']-[', ci2, ']-[', ci3, '] in parameter file. (ERROR)'
    
            end if
    
          end if

        end if

      end do

    end do

    ! angle from water
    !
    if (enefunc%table%num_water > 0) then

      list(1:3) = enefunc%table%water_list(1:3,1)
      ci1 = molecule%atom_cls_name(list(2))
      ci2 = molecule%atom_cls_name(list(1))
      ci3 = molecule%atom_cls_name(list(3))
      do j = 1, nangl_p
        if ((ci1 .eq. par%angl_atom_cls(1, j) .and.  &
             ci2 .eq. par%angl_atom_cls(2, j) .and.  &
             ci3 .eq. par%angl_atom_cls(3, j)) .or.  &
            (ci1 .eq. par%angl_atom_cls(3, j) .and.  &
             ci2 .eq. par%angl_atom_cls(2, j) .and.  &
             ci1 .eq. par%angl_atom_cls(1, j))) then
          angle_calc = .true.
          enefunc%table%HOH_angle = par%angl_theta_min(j)*RAD
          enefunc%table%HOH_force = par%angl_force_const(j)
          if (enefunc%table%HOH_force > EPS) &
            enefunc%table%water_angle_calc = .true.
          exit
        end if
      end do

    end if

    found = 0
    do i = 1, ncel
      found = found + angle(i)
      if (angle_calc) found = found + nwater(i)
      if (angle(i) > MaxAngle) &
        call error_msg('Setup_Enefunc_Angl> Too many angles.')
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_angl_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_angl_all = found
#endif

    if (domain%fep_use) then

      nangl_fep = 0
      do i = 1, 5
        nangl_fep = nangl_fep + molecule%num_angles_fep(i)
      end do
      if (enefunc%num_angl_all /= nangl_fep*domain%num_duplicate) &
        call error_msg('Setup_Enefunc_Angl> Some angle paremeters are missing.')

    else

      if (enefunc%num_angl_all /= nangl*domain%num_duplicate) &
        call error_msg('Setup_Enefunc_Angl> Some angle paremeters are missing.')

    end if

    return

  end subroutine setup_enefunc_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_angl_constraint
  !> @brief        define ANGLE term for each cell in potential energy function
  !                with SETTLE constraint
  !! @authors      JJ, HO
  !! @param[in]    par         : CHARMM PAR information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    domain      : domain information
  !! @param[in]    constraints : constraints information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_angl_constraint(par, molecule, domain, &
                                           constraints, enefunc)

    ! formal arguments
    type(s_par),                 intent(in)    :: par
    type(s_molecule),    target, intent(in)    :: molecule
    type(s_domain),      target, intent(in)    :: domain
    type(s_constraints), target, intent(in)    :: constraints
    type(s_enefunc),     target, intent(inout) :: enefunc

    ! local variables
    integer                  :: dupl, ioffset
    integer                  :: i, j, icel_local
    integer                  :: icel1, icel2
    integer                  :: nangl, nangl_p, found
    integer                  :: list(3), lists(3), pbc_int
    character(6)             :: ci1, ci2, ci3
    character(6)             :: ri1, ri2, ri3
    integer                  :: nangl_per_water
    real(wp)                 :: cwork(3,3), dij(3)

    real(wp),        pointer :: force(:,:), theta(:,:), coord(:,:)
    real(wp),        pointer :: box_size(:)
    real(wp),        pointer :: ubforce(:,:), ubrmin(:,:)
    integer,         pointer :: angle(:), alist(:,:,:), num_water
    integer,         pointer :: ncel
    integer(int2),   pointer :: cell_pair(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: mol_angl_list(:,:), mol_no(:)
    integer,         pointer :: nwater(:)
    character(6),    pointer :: mol_cls_name(:), mol_res_name(:)
    integer,         pointer :: angl_pbc(:,:,:)

    ! FEP
    integer                  :: nangl_fep, k, fg(3)
    logical                  :: flag_singleB
    integer,         pointer :: angl_singleB(:,:)

    mol_angl_list  => molecule%angl_list
    mol_no         => molecule%molecule_no
    mol_cls_name   => molecule%atom_cls_name
    mol_res_name   => molecule%residue_name
    coord         => molecule%atom_coord

    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l
    nwater    => domain%num_water
    box_size  => domain%system_size

    angle     => enefunc%num_angle
    alist     => enefunc%angle_list
    force     => enefunc%angle_force_const
    theta     => enefunc%angle_theta_min
    ubforce   => enefunc%urey_force_const
    ubrmin    => enefunc%urey_rmin
    num_water => enefunc%table%num_water
    angl_pbc  => enefunc%angle_pbc

    nangl     = molecule%num_angles
    nangl_p   = par%num_angles

    nangl_per_water = 0

    ! FEP
    if (domain%fep_use) then
      angl_singleB => enefunc%angl_singleB
    end if

    if (num_water > 0) then

      do i = 1, nangl
        list(1:3) = mol_angl_list(1:3,i)
        ri1 = mol_res_name(list(1))
        ri2 = mol_res_name(list(2))
        ri3 = mol_res_name(list(3))
        if (ri1(1:4) == constraints%water_model .and. &
            ri2(1:4) == constraints%water_model .and. &
            ri3(1:4) == constraints%water_model) then
          nangl_per_water = nangl_per_water + 1
        end if
      end do

      if (mod(nangl_per_water,num_water) /= 0) then
        write(MsgOut,*) &
             'Setup_Enefunc_Angl_Constraint> invalid ANGL count: ', &
             'number of angle terms in a water molecule is not integer.'
        call error_msg('Setup_Enefunc_Angl_Constraint> abort')
      end if
      nangl_per_water = nangl_per_water / num_water

    end if

    do dupl = 1, domain%num_duplicate

      ioffset = (dupl-1) * enefunc%table%num_all

      do i = 1, nangl

        list(1:3) = mol_angl_list(1:3,i)
        ci1 = mol_cls_name(list(1))
        ci2 = mol_cls_name(list(2))
        ci3 = mol_cls_name(list(3))
        ri1 = mol_res_name(list(1))
        ri2 = mol_res_name(list(2))
        ri3 = mol_res_name(list(3))

        if (ri1(1:4) /= constraints%water_model .and. &
            ri2(1:4) /= constraints%water_model .and. &
            ri3(1:4) /= constraints%water_model) then

          ! FEP
          if (domain%fep_use) then
            flag_singleB = .false.
            do j = 1, 3
              fg(j) = molecule%fepgrp(list(j))
            end do
            if (molecule%fepgrp_angl(fg(1),fg(2),fg(3)) == 0) then
              ! FEP: If the angle are not set to any group of FEP, exclude this angle.
              cycle
            else if (molecule%fepgrp_angl(fg(1),fg(2),fg(3)) == 2) then
              ! FEP: If 2-2 and 2-5
              flag_singleB = .true.
            end if
            do j = 1, 3
              ! FEP: If list(j) is singleB atom, its ID is replaced with singleA's ID.
              if (fg(j) == 2) then
                do k = 1, molecule%num_atoms_fep(2)
                  if (molecule%id_singleB(k) == list(j)) then
                    exit
                  end if
                end do
                list(j) = molecule%id_singleA(k)
              end if
            end do
          end if

          lists(1:3) = list(1:3) + ioffset

          icel1 = id_g2l(1,lists(1))
          icel2 = id_g2l(1,lists(3))

          if (icel1 /= 0 .and. icel2 /= 0) then

            icel_local = cell_pair(icel1,icel2)

            if (icel_local >= 1 .and. icel_local <= ncel) then

              do j = 1, nangl_p
                if ((ci1 == par%angl_atom_cls(1,j) .and. &
                     ci2 == par%angl_atom_cls(2,j) .and. &
                     ci3 == par%angl_atom_cls(3,j)) .or. &
                    (ci1 == par%angl_atom_cls(3,j) .and. &
                     ci2 == par%angl_atom_cls(2,j) .and. &
                     ci3 == par%angl_atom_cls(1,j))) then
  
                  angle(icel_local) = angle(icel_local) + 1
                  alist(1:3,angle(icel_local),icel_local) = lists(1:3)
  
                  force(angle(icel_local),icel_local)   = &
                       par%angl_force_const(j)
                  theta(angle(icel_local),icel_local)   = &
                       par%angl_theta_min(j)*RAD
                  ubforce(angle(icel_local),icel_local) = &
                       par%angl_ub_force_const(j)
                  ubrmin(angle(icel_local),icel_local)  = &
                       par%angl_ub_rmin(j)
                  cwork(1:3,1) = coord(1:3,list(1))
                  cwork(1:3,2) = coord(1:3,list(2))
                  cwork(1:3,3) = coord(1:3,list(3))
                  dij(1:3) = cwork(1:3,1) - cwork(1:3,2)
                  call check_pbc(box_size, dij, pbc_int)
                  angl_pbc(1,angle(icel_local),icel_local) = pbc_int
                  dij(1:3) = cwork(1:3,3) - cwork(1:3,2)
                  call check_pbc(box_size, dij, pbc_int)
                  angl_pbc(2,angle(icel_local),icel_local) = pbc_int
                  dij(1:3) = cwork(1:3,1) - cwork(1:3,3)
                  call check_pbc(box_size, dij, pbc_int)
                  angl_pbc(3,angle(icel_local),icel_local) = pbc_int

                  if (domain%fep_use) then
                    ! FEP: flag for singleB angle
                    if (flag_singleB) &
                      angl_singleB (angle(icel_local),icel_local) = 1
                  end if

                  exit
  
                end if
              end do
  
              if (j == nangl_p + 1) &
                write(MsgOut,*) &
                  'Setup_Enefunc_Angl_Constraint> not found ANGL: [', &
                  ci1, ']-[', ci2, ']-[', ci3, '] in parameter file. (ERROR)'
  
            end if
  
          end if
  
        end if

      end do
    end do

    found = 0
    do i = 1, ncel
      found = found + angle(i)
      if (angle(i) > MaxAngle) &
        call error_msg('Setup_Enefunc_Angl_Constraint> Too many angles.')
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_angl_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_angl_all = found
#endif

    if (domain%fep_use) then

      nangl_fep = 0
      do i = 1, 5
        nangl_fep = nangl_fep + molecule%num_angles_fep(i)
      end do
      if (enefunc%num_angl_all /= &
          domain%num_duplicate*(nangl_fep - nangl_per_water*num_water)) &
        call error_msg( &
        'Setup_Enefunc_Angl> Some angle paremeters are missing.')

    else

      if (enefunc%num_angl_all /= &
          domain%num_duplicate*(nangl - nangl_per_water*num_water)) &
        call error_msg( &
        'Setup_Enefunc_Angl_Constraint> Some angle paremeters are missing.')

    end if

    return

  end subroutine setup_enefunc_angl_constraint

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_dihe
  !> @brief        define DIHEDRAL term in potential energy function
  !! @authors      YS, JJ, TM, HO
  !! @param[in]    par      : CHARMM PAR information
  !! @param[in]    molecule : molecule information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_dihe(par, molecule, domain, enefunc)

    ! formal arguments
    type(s_par),      target, intent(in)    :: par
    type(s_molecule), target, intent(in)    :: molecule
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variables
    integer                   :: dupl, ioffset
    integer                   :: ndihe, ndihe_p
    integer                   :: i, j, icel_local
    integer                   :: icel1, icel2
    integer                   :: found, nw_found
    integer                   :: list(4), lists(4), pbc_int
    character(6)              :: ci1, ci2, ci3, ci4
    real(wp)                  :: cwork(3,4), dij(3)

    real(wp),         pointer :: force(:,:), phase(:,:), coord(:,:)
    real(wp),         pointer :: box_size(:)
    integer,          pointer :: dihedral(:), dlist(:,:,:), period(:,:)
    integer,          pointer :: ncel
    integer(int2),    pointer :: cell_pair(:,:)
    integer(int2),    pointer :: id_g2l(:,:)
    logical,      allocatable :: no_wild(:)
    integer,          pointer :: notation
    integer,          pointer :: mol_dihe_list(:,:)
    character(6),     pointer :: mol_cls_name(:)
    integer,          pointer :: dihe_pbc(:,:,:)

    ! FEP
    integer                   :: ndihe_fep, k, fg(4)
    logical                   :: flag_singleB
    integer,          pointer :: dihe_singleB(:,:)

    mol_dihe_list  => molecule%dihe_list
    mol_cls_name   => molecule%atom_cls_name
    coord          => molecule%atom_coord

    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l
    box_size  => domain%system_size

    dihedral  => enefunc%num_dihedral
    dlist     => enefunc%dihe_list
    force     => enefunc%dihe_force_const
    period    => enefunc%dihe_periodicity
    phase     => enefunc%dihe_phase
    notation  => enefunc%notation_14types
    dihe_pbc  => enefunc%dihe_pbc
    notation = 100

    ndihe     = molecule%num_dihedrals
    ndihe_p   = par%num_dihedrals

    ! FEP
    if (domain%fep_use) then
      dihe_singleB => enefunc%dihe_singleB
    end if

    ! check usage of wild card
    !
    allocate(no_wild(ndihe_p))

    do i = 1, ndihe_p

      if ((par%dihe_atom_cls(1,i) .ne. WildCard) .and. &
          (par%dihe_atom_cls(4,i) .ne. WildCard)) then
        ! A-B-C-D type
        no_wild(i) = .true.
      else
        ! X-B-C-D type
        no_wild(i) = .false.
      end if

    end do

    ! find number of interactions
    !
    do dupl = 1, domain%num_duplicate

      ioffset = (dupl-1) * enefunc%table%num_all

      do i = 1, ndihe

        list(1:4) = mol_dihe_list(1:4,i)
        ci1 = mol_cls_name(list(1))
        ci2 = mol_cls_name(list(2))
        ci3 = mol_cls_name(list(3))
        ci4 = mol_cls_name(list(4))

        ! FEP
        if (domain%fep_use) then
          flag_singleB = .false.
          do j = 1, 4
            fg(j) = molecule%fepgrp(list(j))
          end do
          if (molecule%fepgrp_dihe(fg(1),fg(2),fg(3),fg(4)) == 0) then
            ! FEP: If the dihedral are not set to any group of FEP, exclude this dihedral.
            cycle
          else if (molecule%fepgrp_dihe(fg(1),fg(2),fg(3),fg(4)) == 2) then
            ! FEP: 2-2 and 2-5
            flag_singleB = .true.
          end if
          do j = 1, 4
            ! FEP: If list(j) is singleB atom, its ID is replaced with singleA's ID.
            if (fg(j) == 2) then
              do k = 1, molecule%num_atoms_fep(2)
                if (molecule%id_singleB(k) == list(j)) then
                  exit
                end if
              end do
              list(j) = molecule%id_singleA(k)
            end if
          end do
        end if

        lists(1:4) = list(1:4) + ioffset

        icel1 = id_g2l(1,lists(1))
        icel2 = id_g2l(1,lists(4))

        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1,icel2)

          if (icel_local >= 1 .and. icel_local <= ncel) then

            nw_found = 0
            do j = 1, ndihe_p
              if (no_wild(j)) then
                if (((ci1 .eq. par%dihe_atom_cls(1,j)) .and. &
                     (ci2 .eq. par%dihe_atom_cls(2,j)) .and. &
                     (ci3 .eq. par%dihe_atom_cls(3,j)) .and. &
                     (ci4 .eq. par%dihe_atom_cls(4,j))) .or. &
                    ((ci1 .eq. par%dihe_atom_cls(4,j)) .and. &
                     (ci2 .eq. par%dihe_atom_cls(3,j)) .and. &
                     (ci3 .eq. par%dihe_atom_cls(2,j)) .and. &
                     (ci4 .eq. par%dihe_atom_cls(1,j)))) then
  
                  nw_found = nw_found + 1
                  dihedral(icel_local) = dihedral(icel_local) + 1
                  dlist(1:4,dihedral(icel_local),icel_local) = lists(1:4)
  
                  force (dihedral(icel_local),icel_local) = &
                       par%dihe_force_const(j)
                  period(dihedral(icel_local),icel_local) = &
                       par%dihe_periodicity(j)
                  phase (dihedral(icel_local),icel_local) = &
                       par%dihe_phase(j) * RAD
 
                  cwork(1:3,1) = coord(1:3,list(1))
                  cwork(1:3,2) = coord(1:3,list(2))
                  cwork(1:3,3) = coord(1:3,list(3))
                  cwork(1:3,4) = coord(1:3,list(4))
                  dij(1:3) = cwork(1:3,1) - cwork(1:3,2)
                  call check_pbc(box_size, dij, pbc_int)
                  dihe_pbc(1,dihedral(icel_local),icel_local) = pbc_int
                  dij(1:3) = cwork(1:3,2) - cwork(1:3,3)
                  call check_pbc(box_size, dij, pbc_int)
                  dihe_pbc(2,dihedral(icel_local),icel_local) = pbc_int
                  dij(1:3) = cwork(1:3,4) - cwork(1:3,3)
                  call check_pbc(box_size, dij, pbc_int)
                  dihe_pbc(3,dihedral(icel_local),icel_local) = pbc_int

                  if (domain%fep_use) then
                    ! FEP: flag for singleB dihedral
                    if (flag_singleB) &
                      dihe_singleB (dihedral(icel_local),icel_local) = 1
                  end if
  
                  if (period(dihedral(icel_local),icel_local) >  &
                  enefunc%notation_14types) &
                  call error_msg('Setup_Enefunc_Dihe> Too many periodicity.')
  
                end if
              end if
            end do
  
            if (nw_found == 0) then
              do j = 1, ndihe_p
                if (.not.no_wild(j)) then
                  if (((ci2 .eq. par%dihe_atom_cls(2,j)) .and. &
                       (ci3 .eq. par%dihe_atom_cls(3,j))) .or. &
                      ((ci2 .eq. par%dihe_atom_cls(3,j)) .and. &
                       (ci3 .eq. par%dihe_atom_cls(2,j)))) then
  
                    dihedral(icel_local) = dihedral(icel_local) + 1
                    dlist(1:4,dihedral(icel_local),icel_local) = lists(1:4)
  
                    force (dihedral(icel_local),icel_local) = &
                         par%dihe_force_const(j)
                    period(dihedral(icel_local),icel_local) = &
                         par%dihe_periodicity(j)
                    phase (dihedral(icel_local),icel_local) = &
                         par%dihe_phase(j) * RAD
                    cwork(1:3,1) = coord(1:3,list(1))
                    cwork(1:3,2) = coord(1:3,list(2))
                    cwork(1:3,3) = coord(1:3,list(3))
                    cwork(1:3,4) = coord(1:3,list(4))
                    dij(1:3) = cwork(1:3,1) - cwork(1:3,2)
                    call check_pbc(box_size, dij, pbc_int)
                    dihe_pbc(1,dihedral(icel_local),icel_local) = pbc_int
                    dij(1:3) = cwork(1:3,2) - cwork(1:3,3)
                    call check_pbc(box_size, dij, pbc_int)
                    dihe_pbc(2,dihedral(icel_local),icel_local) = pbc_int
                    dij(1:3) = cwork(1:3,4) - cwork(1:3,3)
                    call check_pbc(box_size, dij, pbc_int)
                    dihe_pbc(3,dihedral(icel_local),icel_local) = pbc_int

                    if (domain%fep_use) then
                      ! FEP: flag for singleB dihedral
                      if (flag_singleB) &
                        dihe_singleB (dihedral(icel_local),icel_local) = 1
                    end if
 
                    if (period(dihedral(icel_local),icel_local) >  &
                    enefunc%notation_14types) &
                    call error_msg('Setup_Enefunc_Dihe> Too many periodicity.')
  
                  end if
                end if
              end do
            end if
  
          end if

        end if

      end do

    end do

    deallocate(no_wild)

    found = 0
    do i = 1, ncel
      found = found + dihedral(i)
      if (dihedral(i) > MaxDihe) &
        call error_msg('Setup_Enefunc_Dihe> Too many dihedral angles.')
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_dihe_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_dihe_all = found
#endif

    if (domain%fep_use) then

      ndihe_fep = 0
      do i = 1, 5
        ndihe_fep = ndihe_fep + molecule%num_dihedrals_fep(i)
      end do
      if (enefunc%num_dihe_all < ndihe_fep*domain%num_duplicate) &
        call error_msg( &
           'Setup_Enefunc_Dihe> Some dihedral paremeters are missing.')

    else

      if (enefunc%num_dihe_all < ndihe*domain%num_duplicate) &
        call error_msg( &
           'Setup_Enefunc_Dihe> Some dihedral paremeters are missing.')

    end if

    return

  end subroutine setup_enefunc_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_impr
  !> @brief        define IMPROPER term in potential energy function
  !! @authors      YS,JJ,HO
  !! @param[in]    par      : CHARMM PAR information
  !! @param[in]    molecule : molecule information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_impr(par, molecule, domain, enefunc)

    ! formal variables
    type(s_par),      target, intent(in)    :: par
    type(s_molecule), target, intent(in)    :: molecule
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variables
    integer                   :: dupl, ioffset
    integer                   :: nimpr, nimpr_p
    integer                   :: i, j, icel_local
    integer                   :: icel1, icel2
    integer                   :: found
    integer                   :: list(4), lists(4), pbc_int
    character(6)              :: ci1, ci2, ci3, ci4
    real(wp)                  :: cwork(3,4), dij(3)

    real(wp),         pointer :: force(:,:), phase(:,:), coord(:,:)
    real(wp),         pointer :: box_size(:)
    integer,          pointer :: improper(:), ilist(:,:,:)
    integer,          pointer :: ncel
    integer(int2),    pointer :: cell_pair(:,:)
    integer(int2),    pointer :: id_g2l(:,:)
    integer,      allocatable :: wc_type(:)
    logical,      allocatable :: no_wild(:)
    integer,          pointer :: mol_impr_list(:,:)
    character(6),     pointer :: mol_cls_name(:)
    integer,          pointer :: impr_pbc(:,:,:)

    ! FEP
    integer                   :: nimpr_fep, k, fg(4)
    logical                   :: flag_singleB
    integer,          pointer :: impr_singleB(:,:)

    mol_impr_list  => molecule%impr_list
    mol_cls_name   => molecule%atom_cls_name
    coord          => molecule%atom_coord

    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l
    box_size  => domain%system_size

    improper  => enefunc%num_improper
    ilist     => enefunc%impr_list
    force     => enefunc%impr_force_const
    phase     => enefunc%impr_phase
    impr_pbc  => enefunc%impr_pbc

    nimpr     = molecule%num_impropers
    nimpr_p   = par%num_impropers

    ! FEP
    if (domain%fep_use) then
      impr_singleB => enefunc%impr_singleB
    end if

    ! check usage of wild card
    !
    allocate(wc_type(nimpr_p), no_wild(nimpr))

    do i = 1, nimpr_p

      if ((par%impr_atom_cls(1,i) .ne. WildCard) .and. &
          (par%impr_atom_cls(2,i) .ne. WildCard) .and. &
          (par%impr_atom_cls(3,i) .ne. WildCard) .and. &
          (par%impr_atom_cls(4,i) .ne. WildCard)) then

        ! A-B-C-D type
        wc_type(i) = 0

      else if ((par%impr_atom_cls(1,i) .eq. WildCard) .and. &
               (par%impr_atom_cls(2,i) .ne. WildCard) .and. &
               (par%impr_atom_cls(3,i) .ne. WildCard) .and. &
               (par%impr_atom_cls(4,i) .ne. WildCard)) then

        ! X-B-C-D type
        wc_type(i) = 1

      else if ((par%impr_atom_cls(1,i) .eq. WildCard) .and. &
               (par%impr_atom_cls(2,i) .eq. WildCard) .and. &
               (par%impr_atom_cls(3,i) .ne. WildCard) .and. &
               (par%impr_atom_cls(4,i) .ne. WildCard)) then

        ! X-X-C-D type
        wc_type(i) = 2

      else if ((par%impr_atom_cls(1,i) .ne. WildCard) .and. &
               (par%impr_atom_cls(2,i) .eq. WildCard) .and. &
               (par%impr_atom_cls(3,i) .eq. WildCard) .and. &
               (par%impr_atom_cls(4,i) .ne. WildCard)) then

        ! A-X-X-D type
        wc_type(i) = 3

      else
        call error_msg('Setup_Enefunc_Impr> Undefined Wild Card')

      end if

    end do

    ! setup parameters
    !
    do dupl = 1, domain%num_duplicate

      ioffset = (dupl-1) * enefunc%table%num_all

      do i = 1, nimpr

        no_wild(i) = .false.
  
        list(1:4) = mol_impr_list(1:4,i)
        ci1 = mol_cls_name(list(1))
        ci2 = mol_cls_name(list(2))
        ci3 = mol_cls_name(list(3))
        ci4 = mol_cls_name(list(4))

        ! FEP
        if (domain%fep_use) then
          flag_singleB = .false.
          do j = 1, 4
            fg(j) = molecule%fepgrp(list(j))
          end do
          if (molecule%fepgrp_dihe(fg(1),fg(2),fg(3),fg(4)) == 0) then
            ! FEP: If the dihedral are not set to any group of FEP, exclude this dihedral.
            cycle
          else if (molecule%fepgrp_dihe(fg(1),fg(2),fg(3),fg(4)) == 2) then
            ! FEP: 2-2 and 2-5
            flag_singleB = .true.
          end if
          do j = 1, 4
            ! FEP: If list(j) is singleB atom, its ID is replaced with singleA's ID.
            if (fg(j) == 2) then
              do k = 1, molecule%num_atoms_fep(2)
                if (molecule%id_singleB(k) == list(j)) then
                  exit
                end if
              end do
              list(j) = molecule%id_singleA(k)
            end if
          end do
        end if

        lists(1:4) = list(1:4) + ioffset

        icel1 = id_g2l(1,lists(1))
        icel2 = id_g2l(1,lists(4))

        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1,icel2)

          if (icel_local >= 1 .and. icel_local <= ncel) then

            ! A-B-C-D type
            !
            do j = 1, nimpr_p
              if (wc_type(j) == 0) then
                if (((ci1 .eq. par%impr_atom_cls(1,j)) .and. &
                     (ci2 .eq. par%impr_atom_cls(2,j)) .and. &
                     (ci3 .eq. par%impr_atom_cls(3,j)) .and. &
                     (ci4 .eq. par%impr_atom_cls(4,j))) .or. &
                    ((ci1 .eq. par%impr_atom_cls(4,j)) .and. &
                     (ci2 .eq. par%impr_atom_cls(3,j)) .and. &
                     (ci3 .eq. par%impr_atom_cls(2,j)) .and. &
                     (ci4 .eq. par%impr_atom_cls(1,j)))) then
  
                  improper(icel_local) = improper(icel_local) + 1
                  ilist(1:4,improper(icel_local),icel_local) = lists(1:4)
  
                  force(improper(icel_local),icel_local) = &
                       par%impr_force_const(j)
                  phase(improper(icel_local),icel_local) = &
                       par%impr_phase(j) * RAD
                  no_wild(i) = .true.
                  cwork(1:3,1) = coord(1:3,list(1))
                  cwork(1:3,2) = coord(1:3,list(2))
                  cwork(1:3,3) = coord(1:3,list(3))
                  cwork(1:3,4) = coord(1:3,list(4))
                  dij(1:3) = cwork(1:3,1) - cwork(1:3,2)
                  call check_pbc(box_size, dij, pbc_int)
                  impr_pbc(1,improper(icel_local),icel_local) = pbc_int
                  dij(1:3) = cwork(1:3,2) - cwork(1:3,3)
                  call check_pbc(box_size, dij, pbc_int)
                  impr_pbc(2,improper(icel_local),icel_local) = pbc_int
                  dij(1:3) = cwork(1:3,4) - cwork(1:3,3)
                  call check_pbc(box_size, dij, pbc_int)
                  impr_pbc(3,improper(icel_local),icel_local) = pbc_int

                  if (domain%fep_use) then
                    ! FEP: flag for singleB improper
                    if (flag_singleB) &
                      impr_singleB (improper(icel_local),icel_local) = 1
                  end if

                  exit
  
                end if
              end if
            end do
  
            ! X-B-C-D type
            !
              if (.not.no_wild(i)) then
              do j = 1, nimpr_p
                if (wc_type(j) == 1) then
                  if (((ci2 .eq. par%impr_atom_cls(2,j)) .and. &
                       (ci3 .eq. par%impr_atom_cls(3,j)) .and. &
                       (ci4 .eq. par%impr_atom_cls(4,j))) .or. &
                      ((ci2 .eq. par%impr_atom_cls(4,j)) .and. &
                       (ci3 .eq. par%impr_atom_cls(3,j)) .and. &
                       (ci4 .eq. par%impr_atom_cls(2,j)))) then
  
                    improper(icel_local) = improper(icel_local) + 1
                    ilist(1:4,improper(icel_local),icel_local) = lists(1:4)
  
                    force(improper(icel_local),icel_local) = &
                         par%impr_force_const(j)
                    phase(improper(icel_local),icel_local) = &
                         par%impr_phase(j) * RAD
                    no_wild(i) = .true.
                    cwork(1:3,1) = coord(1:3,list(1))
                    cwork(1:3,2) = coord(1:3,list(2))
                    cwork(1:3,3) = coord(1:3,list(3))
                    cwork(1:3,4) = coord(1:3,list(4))
                    dij(1:3) = cwork(1:3,1) - cwork(1:3,2)
                    call check_pbc(box_size, dij, pbc_int)
                    impr_pbc(1,improper(icel_local),icel_local) = pbc_int
                    dij(1:3) = cwork(1:3,2) - cwork(1:3,3)
                    call check_pbc(box_size, dij, pbc_int)
                    impr_pbc(2,improper(icel_local),icel_local) = pbc_int
                    dij(1:3) = cwork(1:3,4) - cwork(1:3,3)
                    call check_pbc(box_size, dij, pbc_int)
                    impr_pbc(3,improper(icel_local),icel_local) = pbc_int

                    if (domain%fep_use) then
                      ! FEP: flag for singleB improper
                      if (flag_singleB) &
                        impr_singleB (improper(icel_local),icel_local) = 1
                    end if

                    exit
  
                  end if
                end if
              end do
            end if
    
            ! X-X-C-D type
            !
            if (.not.no_wild(i)) then
              do j = 1, nimpr_p
                if (wc_type(j) == 2) then
                  if (((ci3 .eq. par%impr_atom_cls(3,j)) .and. &
                       (ci4 .eq. par%impr_atom_cls(4,j))) .or. &
                      ((ci3 .eq. par%impr_atom_cls(4,j)) .and. &
                       (ci4 .eq. par%impr_atom_cls(3,j)))) then
  
                    improper(icel_local) = improper(icel_local) + 1
                    ilist(1:4,improper(icel_local),icel_local) = lists(1:4)
  
                    force(improper(icel_local),icel_local) = &
                           par%impr_force_const(j)
                    phase(improper(icel_local),icel_local) = &
                         par%impr_phase(j) * RAD
                    no_wild(i) = .true.
                    cwork(1:3,1) = coord(1:3,list(1))
                    cwork(1:3,2) = coord(1:3,list(2))
                    cwork(1:3,3) = coord(1:3,list(3))
                    cwork(1:3,4) = coord(1:3,list(4))
                    dij(1:3) = cwork(1:3,1) - cwork(1:3,2)
                    call check_pbc(box_size, dij, pbc_int)
                    impr_pbc(1,improper(icel_local),icel_local) = pbc_int
                    dij(1:3) = cwork(1:3,2) - cwork(1:3,3)
                    call check_pbc(box_size, dij, pbc_int)
                    impr_pbc(2,improper(icel_local),icel_local) = pbc_int
                    dij(1:3) = cwork(1:3,4) - cwork(1:3,3)
                    call check_pbc(box_size, dij, pbc_int)
                    impr_pbc(3,improper(icel_local),icel_local) = pbc_int

                    if (domain%fep_use) then
                      ! FEP: flag for singleB improper
                      if (flag_singleB) &
                        impr_singleB (improper(icel_local),icel_local) = 1
                    end if

                    exit
  
                  end if
                end if
              end do
            end if
  
            ! A-X-X-D type
            !
            if (.not.no_wild(i)) then
              do j = 1, nimpr_p
                if (wc_type(j) == 3) then
                  if (((ci1 .eq. par%impr_atom_cls(1,j)) .and. &
                       (ci4 .eq. par%impr_atom_cls(4,j))) .or. &
                      ((ci1 .eq. par%impr_atom_cls(4,j)) .and. &
                       (ci4 .eq. par%impr_atom_cls(1,j)))) then
      
                    improper(icel_local) = improper(icel_local) + 1
                    ilist(1:4,improper(icel_local),icel_local) = lists(1:4)
  
                    force(improper(icel_local),icel_local) = &
                         par%impr_force_const(j)
                    phase(improper(icel_local),icel_local) = &
                         par%impr_phase(j) * RAD
                    no_wild(i) = .true.
                    cwork(1:3,1) = coord(1:3,list(1))
                    cwork(1:3,2) = coord(1:3,list(2))
                    cwork(1:3,3) = coord(1:3,list(3))
                    cwork(1:3,4) = coord(1:3,list(4))
                    dij(1:3) = cwork(1:3,1) - cwork(1:3,2)
                    call check_pbc(box_size, dij, pbc_int)
                    impr_pbc(1,improper(icel_local),icel_local) = pbc_int
                    dij(1:3) = cwork(1:3,2) - cwork(1:3,3)
                    call check_pbc(box_size, dij, pbc_int)
                    impr_pbc(2,improper(icel_local),icel_local) = pbc_int
                    dij(1:3) = cwork(1:3,4) - cwork(1:3,3)
                    call check_pbc(box_size, dij, pbc_int)
                    impr_pbc(3,improper(icel_local),icel_local) = pbc_int

                    if (domain%fep_use) then
                      ! FEP: flag for singleB improper
                      if (flag_singleB) &
                        impr_singleB (improper(icel_local),icel_local) = 1
                    end if

                    exit
  
                  end if
                end if
              end do
            end if

            if (.not.no_wild(i)) &
              write(MsgOut,*) &
                'Setup_Enefunc_Impr> Unknown IMPR type. [', &
                ci1, ']-[', ci2, ']-[', ci3, ']-[', ci4, '] (ERROR)'
  
          end if
        end if

      end do

    end do

    deallocate(wc_type, no_wild)

    found = 0
    do i = 1, ncel
      found = found + improper(i)
      if (improper(i) > MaxImpr) &
        call error_msg('Setup_Enefunc_Impr> Too many improper dihedral angles')
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_impr_all, 1, mpi_integer, mpi_sum, &
                       mpi_comm_country, ierror)
#else
    enefunc%num_impr_all = found
#endif

    if (domain%fep_use) then

      nimpr_fep = 0
      do i = 1, 5
        nimpr_fep = nimpr_fep + molecule%num_impropers_fep(i)
      end do
      if (enefunc%num_impr_all < nimpr_fep*domain%num_duplicate) &
        call error_msg( &
          'Setup_Enefunc_Impr> Some improper paremeters are missing.')

    else

      if (enefunc%num_impr_all < nimpr*domain%num_duplicate) &
        call error_msg( &
          'Setup_Enefunc_Impr> Some improper paremeters are missing.')

    end if

    return

  end subroutine setup_enefunc_impr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cmap
  !> @brief        define cmap term in potential energy function with DD
  !! @authors      TY, TM, HO
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    par         : CHARMM PAR information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    domain      : domain information
  !! @param[inout] enefunc     : energy potential functions informationn
  !!
  !! @note       In str "par", following variables have been defined:
  !!   cmap_atom_cls(8,imap) (char) : 8 atom classes (4 for phi and 4 for psi)
  !!   cmap_resolution(imap) (int ) : = 24 (for all imap) = 360degree/15degree
  !!   cmap_data(i,j,imap)   (real) : Ecmap for grid points. i and j are the
  !!                                  1st (psi?) and the 2nd (phi?) grid IDs.
  !!                                  1<=i,j<=24.
  !!   Where imap is ID of cmap type (1 <= imap <= 6).
  !!
  !! @note       TY determined to use natural spline (periodic = .false.)
  !!             because it is similar to the CHARMM way.
  !!             However, I notice that force will be more accurately continuous
  !!             at phi (or psi) = +-180 degree if periodic spline was used.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cmap(ene_info, par, molecule, domain, enefunc)

    ! formal variables
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_par),             intent(in)    :: par
    type(s_molecule),target, intent(in)    :: molecule
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: dupl, ioffset
    integer                  :: i, j, k, l, ityp, pbc_int
    integer                  :: ncmap_p, found, ngrid0
    integer                  :: list(8), lists(8), icel1, icel2, icel_local
    integer                  :: flag_cmap_type, alloc_stat, dealloc_stat
    character(6)             :: ci1, ci2, ci3, ci4, ci5, ci6, ci7, ci8
    logical                  :: periodic
    real(wp)                 :: cwork(3,8), dij(3)

    integer,         pointer :: ncel
    integer(int2),   pointer :: cell_pair(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: mol_cmap_list(:,:)
    character(6),    pointer :: mol_cls_name(:)
    real(wp),        pointer :: coord(:,:), box_size(:)

    real(wp),    allocatable :: c_ij(:,:,:,:) ! cmap coeffs

    ! FEP
    integer                  :: ncmap_fep, idx, fg(8)
    logical                  :: flag_singleB

    mol_cmap_list => molecule%cmap_list
    mol_cls_name  => molecule%atom_cls_name
    coord         => molecule%atom_coord

    ncel          => domain%num_cell_local
    cell_pair     => domain%cell_pair
    id_g2l        => domain%id_g2l
    box_size      => domain%system_size

    ! If 'periodic' is .true.,
    !   then cubic spline with periodic (in dihedral-angle space) boundary
    !   will be applied.
    ! If 'periodic' is .false.,
    !   then natural cubic spline without periodic boudnary
    !   will be applied to expanded cross-term maps.
    !   This is similar with CHARMM's source code.
    !
    periodic = ene_info%cmap_pspline
    ncmap_p  = par%num_cmaps

    ! begin
    !
    ngrid0 = 0
    do i = 1, ncmap_p
      ngrid0 = max(ngrid0, par%cmap_resolution(i))
    end do

    call alloc_enefunc(enefunc, EneFuncCmap, ncel, ngrid0, ncmap_p)

    alloc_stat = 0
    allocate(c_ij(4,4,ngrid0,ngrid0), stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    do i = 1, ncmap_p
      enefunc%cmap_resolution(i) = par%cmap_resolution(i)
    end do

    ! derive cmap coefficients by bicubic interpolation
    !
    do ityp = 1, ncmap_p

      if (periodic) then
        call derive_cmap_coefficients_p(ityp, par, c_ij)
      else
        call derive_cmap_coefficients_np(ityp, par, c_ij)
      end if

      do l = 1, ngrid0
        do k = 1, ngrid0
          do j = 1, 4
            do i = 1, 4
              enefunc%cmap_coef(i,j,k,l,ityp) = c_ij(i,j,k,l)
            end do
          end do
        end do
      end do
    end do

    enefunc%num_cmap(1:ncel) = 0

    do dupl = 1, domain%num_duplicate

      ioffset = (dupl-1) * enefunc%table%num_all

      do i = 1, molecule%num_cmaps

        list(1:8) = mol_cmap_list(1:8,i) 

        ! ci* will be atom-type strings
        !
        ci1 = mol_cls_name(list(1))
        ci2 = mol_cls_name(list(2))
        ci3 = mol_cls_name(list(3))
        ci4 = mol_cls_name(list(4))
        ci5 = mol_cls_name(list(5))
        ci6 = mol_cls_name(list(6))
        ci7 = mol_cls_name(list(7))
        ci8 = mol_cls_name(list(8))

        ! FEP
        if (domain%fep_use) then
          flag_singleB = .false.
          do j = 1, 8
            fg(j) = molecule%fepgrp(list(j))
          end do
          idx = fg(1) + 5*(fg(2)-1 + 5*(fg(3)-1 + 5*(fg(4)-1 + &
            5*(fg(5)-1 + 5*(fg(6)-1 + 5*(fg(7)-1 + 5*(fg(8)-1)))))))
          if (molecule%fepgrp_cmap(idx) == 0) then
            ! FEP: If the cmap are not set to any group of FEP, exclude this cmap.
            cycle
          else if (molecule%fepgrp_cmap(idx) == 2) then
            ! FEP: 2-2 and 2-5
            flag_singleB = .true.
          end if
          do j = 1, 8
            ! FEP: If list(j) is singleB atom, its ID is replaced with singleA's ID.
            if (fg(j) == 2) then
              do k = 1, molecule%num_atoms_fep(2)
                if (molecule%id_singleB(k) == list(j)) then
                  exit
                end if
              end do
              list(j) = molecule%id_singleA(k)
            end if
          end do
        end if

        lists(1:8) = list(1:8)+ioffset

        icel1 = id_g2l(1,lists(1))
        icel2 = id_g2l(1,lists(8))

        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1,icel2)

          if (icel_local >= 1 .and. icel_local <= ncel) then

            flag_cmap_type = -1
  
            ! assign cmap type ID to each (psi,phi) pair
            !
            do j = 1, ncmap_p
              if (ci1 .eq. par%cmap_atom_cls(1, j) .and. &
                  ci2 .eq. par%cmap_atom_cls(2, j) .and. &
                  ci3 .eq. par%cmap_atom_cls(3, j) .and. &
                  ci4 .eq. par%cmap_atom_cls(4, j) .and. &
                  ci5 .eq. par%cmap_atom_cls(5, j) .and. &
                  ci6 .eq. par%cmap_atom_cls(6, j) .and. &
                  ci7 .eq. par%cmap_atom_cls(7, j) .and. &
                  ci8 .eq. par%cmap_atom_cls(8, j)) then
  
                enefunc%num_cmap(icel_local) = enefunc%num_cmap(icel_local) + 1
                k = enefunc%num_cmap(icel_local)
                enefunc%cmap_list(1:8,k,icel_local) = lists(1:8)
                enefunc%cmap_type(k,icel_local) = j
                flag_cmap_type = j
                cwork(1:3,1) = coord(1:3,list(1))
                cwork(1:3,2) = coord(1:3,list(2))
                cwork(1:3,3) = coord(1:3,list(3))
                cwork(1:3,4) = coord(1:3,list(4))
                cwork(1:3,5) = coord(1:3,list(5))
                cwork(1:3,6) = coord(1:3,list(6))
                cwork(1:3,7) = coord(1:3,list(7))
                cwork(1:3,8) = coord(1:3,list(8))
                dij(1:3) = cwork(1:3,1) - cwork(1:3,2)
                call check_pbc(box_size, dij, pbc_int)
                enefunc%cmap_pbc(1,k,icel_local) = pbc_int
                dij(1:3) = cwork(1:3,2) - cwork(1:3,3)
                call check_pbc(box_size, dij, pbc_int)
                enefunc%cmap_pbc(2,k,icel_local) = pbc_int
                dij(1:3) = cwork(1:3,4) - cwork(1:3,3)
                call check_pbc(box_size, dij, pbc_int)
                enefunc%cmap_pbc(3,k,icel_local) = pbc_int
                dij(1:3) = cwork(1:3,5) - cwork(1:3,6)
                call check_pbc(box_size, dij, pbc_int)
                enefunc%cmap_pbc(4,k,icel_local) = pbc_int
                dij(1:3) = cwork(1:3,6) - cwork(1:3,7)
                call check_pbc(box_size, dij, pbc_int)
                enefunc%cmap_pbc(5,k,icel_local) = pbc_int
                dij(1:3) = cwork(1:3,8) - cwork(1:3,7)
                call check_pbc(box_size, dij, pbc_int)
                enefunc%cmap_pbc(6,k,icel_local) = pbc_int

                if (domain%fep_use) then
                  ! FEP: flag for singleB cmap
                  if (flag_singleB) &
                   enefunc%cmap_singleB(k,icel_local) = 1
                end if

                exit
  
              end if
            end do

            ! if not found, print detail about the error.
            !
  
            if (flag_cmap_type <= 0) then
              write(MsgOut,*) 'Setup_Enefunc_Cmap> not found CMAP: '
              write(MsgOut,*) ' [',ci1,']-[',ci2,']-[',ci3,']-[',ci4,']-'
              write(MsgOut,*) ' [',ci5,']-[',ci6,']-[',ci7,']-[',ci8,'] '
              write(MsgOut,*) '  in parameter file. (ERROR)'
            end if
  
          end if
        end if

      end do

    end do

    deallocate(c_ij, stat=dealloc_stat)
    if (dealloc_stat /= 0) &
      call error_msg_dealloc

    ! write summary
    !
    if (main_rank) then
      if (periodic) then
        write(MsgOut,'(A)') &
    'Setup_Enefunc_Cmap> Periodic-boundary spline is used to derive cmap coefs.'
        write(MsgOut,'(A)') ''
      else
        write(MsgOut,'(A)') &
            'Setup_Enefunc_Cmap> Natural spline is used to derive cmap coefs.'
        write(MsgOut,'(A)') ''
      end if
    end if

    ! stop if parameter is missing
    !
    found = 0
    do i = 1, ncel
      found = found + enefunc%num_cmap(i)
      if (enefunc%num_cmap(i) > MaxCmap) &
        call error_msg('Setup_Enefunc_Cmap> Too many cmaps.')
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_cmap_all, 1, mpi_integer, mpi_sum, &
                       mpi_comm_country, ierror)
#else
    enefunc%num_cmap_all = found
#endif

    if (domain%fep_use) then

      ncmap_fep = 0
      do i = 1, 5
        ncmap_fep = ncmap_fep + molecule%num_cmaps_fep(i)
      end do
      if (enefunc%num_cmap_all /= ncmap_fep*domain%num_duplicate) &
        call error_msg('Setup_Enefunc_Cmap> Some cmap parameters are missing.')

    else

      if (enefunc%num_cmap_all /=  molecule%num_cmaps*domain%num_duplicate) &
        call error_msg('Setup_Enefunc_Cmap> Some cmap parameters are missing.')

    end if

    return

  end subroutine setup_enefunc_cmap

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_nonb
  !> @brief        define NON-BOND term in potential energy function
  !! @authors      YS, JJ, TI
  !! @param[in]    par         : CHARMM PAR information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    constraints : constraints information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_nonb(par, molecule, constraints, domain, enefunc)

    ! formal arguments
    type(s_par),             intent(in)    :: par
    type(s_molecule),        intent(in)    :: molecule
    type(s_constraints),     intent(in)    :: constraints
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: eps14, rmin14, eps, rmin, lamda_i, lamda_j
    real(dp)                 :: vdw_self1, vdw_self2
    integer                  :: i, j, k, ix, jx, kk
    integer                  :: nonb_p, nbfx_p, cls_local, ncel
    character(6)             :: ci1, ci2

    integer,  allocatable    :: nonb_atom_cls(:), check_cls(:)
    integer,  allocatable    :: atmcls_map_g2l(:), atmcls_map_l2g(:)
    real(wp), allocatable    :: nb14_lj6(:,:), nb14_lj12(:,:)
    real(wp), allocatable    :: nonb_lj6(:,:), nonb_lj12(:,:)
    real(wp), allocatable    :: lj_coef(:,:)


    enefunc%num_atom_cls = par%num_atom_cls

    ELECOEF          = ELECOEF_CHARMM

    ! set lennard-jones parameters
    !
    nonb_p = enefunc%num_atom_cls

    allocate(nonb_atom_cls(nonb_p),     &
             check_cls(nonb_p),         &
             atmcls_map_g2l(nonb_p),    &
             atmcls_map_l2g(nonb_p),    &
             nb14_lj6 (nonb_p, nonb_p), &
             nb14_lj12(nonb_p, nonb_p), &
             nonb_lj6 (nonb_p, nonb_p), &
             nonb_lj12(nonb_p, nonb_p), &
             lj_coef(2,nonb_p))

    nonb_atom_cls(1:nonb_p)           = 0
    check_cls(1:nonb_p)               = 0
    nb14_lj6 (1:nonb_p, 1:nonb_p) = 0.0_wp
    nb14_lj12(1:nonb_p, 1:nonb_p) = 0.0_wp
    nonb_lj6 (1:nonb_p, 1:nonb_p) = 0.0_wp
    nonb_lj12(1:nonb_p, 1:nonb_p) = 0.0_wp
    lj_coef  (1:2,1:nonb_p)       = 0.0_wp 

    do i = 1, nonb_p

      lj_coef(1,i) = abs(par%nonb_eps(i))
      lj_coef(2,i) = par%nonb_rmin(i)

      do j = 1, nonb_p

        ! combination rule
        eps14  = sqrt(par%nonb_eps_14(i) * par%nonb_eps_14(j))
        rmin14 = par%nonb_rmin_14(i) + par%nonb_rmin_14(j)
        eps    = sqrt(par%nonb_eps(i) * par%nonb_eps(j))
        rmin   = par%nonb_rmin(i) + par%nonb_rmin(j)

        ! set parameters
        if (eps14 > 1.0e-15) then
          nb14_lj12(i,j) = eps14 * (rmin14 ** 12)
          nb14_lj6(i,j)  = 2.0_wp * eps14 * (rmin14 ** 6)
        else
          nb14_lj12(i,j) = 0.0_wp
          nb14_lj6(i,j)  = 0.0_wp
        end if
        if (eps > 1.0e-15) then
          nonb_lj12(i,j) = eps * (rmin ** 12)
          nonb_lj6(i,j)  = 2.0_wp * eps * (rmin ** 6)
        else
          nonb_lj12(i,j) = 0.0_wp
          nonb_lj6(i,j)  = 0.0_wp
        end if
      end do
    end do

    ! overwrite lennard-jones parameters by NBFIX parameters
    !
    nbfx_p = par%num_nbfix

    do k = 1, nbfx_p
      ci1 = par%nbfi_atom_cls(1,k)
      ci2 = par%nbfi_atom_cls(2,k)
      do i = 1, nonb_p
        do j = 1, nonb_p
          if ((ci1 .eq. par%nonb_atom_cls(i)  .and. &
               ci2 .eq. par%nonb_atom_cls(j)) .or.  &
              (ci2 .eq. par%nonb_atom_cls(i)  .and. &
               ci1 .eq. par%nonb_atom_cls(j))) then

            ! combination rule
            !
            eps14  = abs(par%nbfi_eps_14 (k)) !TODO CHECK
            rmin14 = par%nbfi_rmin_14(k)
            eps    = abs(par%nbfi_eps    (k))
            rmin   = par%nbfi_rmin   (k)

            ! set parameters
            !
            nb14_lj12(i,j) = eps14 * (rmin14 ** 12)
            nb14_lj6 (i,j) = 2.0_wp * eps14 * (rmin14 ** 6)
            nonb_lj12(i,j) = eps * (rmin ** 12)
            nonb_lj6 (i,j) = 2.0_wp * eps * (rmin ** 6)
          end if
        end do
      end do
    end do

    ! check the usage of atom class
    !
    do i = 1, molecule%num_atoms
      k = molecule%atom_cls_no(i)
      if (k < 1) then
        call error_msg( &
        'Setup_Enefunc_Nonb> atom class is not defined: "'&
        //trim(molecule%atom_cls_name(i))//'"')
      end if
      check_cls(k) = check_cls(k) + 1
    end do

    k = 0
    do i = 1, nonb_p
      if (check_cls(i) > 0) then
        k = k + 1
        atmcls_map_g2l(i) = k
        atmcls_map_l2g(k) = i
      end if
    end do
    cls_local = k
    max_class = cls_local

    call alloc_enefunc(enefunc, EneFuncNbon, cls_local)

    do i = 1, cls_local
      ix = atmcls_map_l2g(i)
      enefunc%nonb_lj6_factor(i)  = sqrt(2.0_wp)*8.0_wp*sqrt(lj_coef(1,ix))  &
                                   *lj_coef(2,ix)**3
      enefunc%nonb_atom_cls_no(i) = check_cls(ix)
    end do

    do i = 1, cls_local
      ix = atmcls_map_l2g(i)
      do j = 1, cls_local
        jx = atmcls_map_l2g(j)
        enefunc%nb14_lj12(i,j) = nb14_lj12(ix,jx)
        enefunc%nb14_lj6 (i,j) = nb14_lj6 (ix,jx)
        enefunc%nonb_lj12(i,j) = nonb_lj12(ix,jx)
        enefunc%nonb_lj6 (i,j) = nonb_lj6 (ix,jx)
      end do
    end do

    ! Sum of lamda_i*lamda_j and lamda_i*lamda_i
    !
    vdw_self1 = 0.0_dp
    vdw_self2 = 0.0_dp
    k = 0
    do i = 1, cls_local
      ix = enefunc%nonb_atom_cls_no(i)
      lamda_i = enefunc%nonb_lj6_factor(i)
      vdw_self2 = vdw_self2 + ix*lamda_i*lamda_i
      do j = 1, cls_local 
        jx = enefunc%nonb_atom_cls_no(j)
        lamda_j = enefunc%nonb_lj6_factor(j)
        vdw_self1 = vdw_self1 + ix*jx*lamda_i*lamda_j
      end do
    end do
    enefunc%pme_dispersion_self1 = real(vdw_self1,wp)
    enefunc%pme_dispersion_self2 = real(vdw_self2,wp)

    ! update domain information
    !
    do i = 1, domain%num_cell_local+domain%num_cell_boundary
      do ix = 1, domain%num_atom(i)
        domain%atom_cls_no(ix,i) = atmcls_map_g2l(domain%atom_cls_no(ix,i))
      end do
    end do

    if (enefunc%table%num_water > 0) then

      domain%water%atom_cls_no(1:3)  &
                             = atmcls_map_g2l(domain%water%atom_cls_no(1:3))
      enefunc%table%atom_cls_no_O    &
                             = atmcls_map_g2l(enefunc%table%atom_cls_no_O)
      enefunc%table%atom_cls_no_H    &
                             = atmcls_map_g2l(enefunc%table%atom_cls_no_H)
      if (constraints%water_type == TIP4) then
        domain%water%atom_cls_no(4)  &
                             = atmcls_map_g2l(domain%water%atom_cls_no(4))
        enefunc%table%atom_cls_no_D  &
                             = atmcls_map_g2l(enefunc%table%atom_cls_no_D)
      end if

    end if

    ! FEP
    if (domain%fep_use) then
      do i = 1, domain%num_cell_local+domain%num_cell_boundary
        do ix = 1, domain%num_atom(i)
          if (domain%fepgrp(ix,i) == 1) then
            domain%fep_atmcls_singleB(ix,i) = atmcls_map_g2l(domain%fep_atmcls_singleB(ix,i))
          end if
        end do
      end do
    end if

    deallocate(nonb_atom_cls,  &
               check_cls,      &
               atmcls_map_g2l, &
               atmcls_map_l2g, &
               nb14_lj6,       &
               nb14_lj12,      &
               nonb_lj6,       &
               nonb_lj12,      & 
               lj_coef)

    enefunc%num_atom_cls = cls_local

    ! treatment for 1-2, 1-3, 1-4 interactions
    !
    ncel   = domain%num_cell_local

    call alloc_enefunc(enefunc, EneFuncNonb,     ncel, maxcell_near)
    call alloc_enefunc(enefunc, EneFuncNonbList, ncel, maxcell_near)

    if (domain%fep_use) then
      ! FEP
      call alloc_enefunc(enefunc, EneFuncNonbFEPList, ncel, maxcell_near)
      if (constraints%rigid_bond) then
        call count_nonb_excl_fep(.true., .true., constraints, domain, enefunc)
      else
        call count_nonb_excl_fep(.true., .false., constraints, domain, enefunc)
      end if

    else

      if (constraints%rigid_bond) then
        call count_nonb_excl(.true., .true., constraints, domain, enefunc)
      else
        call count_nonb_excl(.true., .false., constraints, domain, enefunc)
      end if

    end if

    return

  end subroutine setup_enefunc_nonb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    count_nonb_excl
  !> @brief        exclude 1-2, 1-3 interactions and constraints
  !! @authors      JJ
  !! @param[in]    first       : flag for first call or not
  !! @param[in]    constraint  : flag for constraint usage
  !! @param[in]    constraints : constraints information   
  !! @param[inout] domain      : structure of domain
  !! @param[inout] enefunc     : structure of enefunc
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine count_nonb_excl(first, constraint, constraints, domain, enefunc)

    ! formal arguments
    logical,                     intent(in)    :: first
    logical,                     intent(in)    :: constraint
    type(s_constraints), target, intent(in)    :: constraints
    type(s_domain),      target, intent(inout) :: domain
    type(s_enefunc),     target, intent(inout) :: enefunc

    ! local variables
    integer                  :: ncell, ncell_local, i, ii, ix, k, i1, i2, i3
    integer                  :: icel, icel1, icel2
    integer                  :: ic, j, ih, ij, index(4)
    integer                  :: num_excl, num_nb14, id, omp_get_thread_num
    integer                  :: found1, found2
    integer                  :: fkind
    integer                  :: list1, list2

    integer,         pointer :: natom(:), nwater(:), start_atom(:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: water_list(:,:,:)
    integer,         pointer :: nbond(:), bond_list(:,:,:)
    integer,         pointer :: nangle(:), angl_list(:,:,:)
    integer(1),      pointer :: bond_kind(:,:), angl_kind(:,:), dihe_kind(:,:)
    integer,         pointer :: ndihedral(:), dihe_list(:,:,:)
    integer,         pointer :: cell_pairlist2(:,:)
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: nonb_excl_list(:,:,:)
    integer,         pointer :: nb14_calc_list(:,:,:)
    integer(1),      pointer :: exclusion_mask(:,:,:), exclusion_mask1(:,:,:)
    integer,         pointer :: sc_calc_list(:,:)
    integer,         pointer :: num_nonb_excl(:), num_nb14_calc(:)
    integer,         pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)
    real(wp),        pointer :: nb14_qq_scale(:,:), nb14_lj_scale(:,:)
    real(wp),        pointer :: dihe_scnb(:), dihe_scee(:)
    real(wp),        pointer :: charge(:,:)


    natom           => domain%num_atom
    nwater          => domain%num_water
    start_atom      => domain%start_atom
    water_list      => domain%water_list
    id_g2l          => domain%id_g2l
    cell_pair       => domain%cell_pairlist1
    cell_pairlist2  => domain%cell_pairlist2
    charge          => domain%charge

    nbond           => enefunc%num_bond
    nangle          => enefunc%num_angle
    ndihedral       => enefunc%num_dihedral
    bond_list       => enefunc%bond_list
    bond_kind       => enefunc%bond_kind
    angl_list       => enefunc%angle_list
    angl_kind       => enefunc%angle_kind
    dihe_list       => enefunc%dihe_list
    dihe_kind       => enefunc%dihe_kind
    nonb_excl_list  => enefunc%nonb_excl_list
    nb14_calc_list  => enefunc%nb14_calc_list
    sc_calc_list    => enefunc%sc_calc_list
    num_nonb_excl   => enefunc%num_nonb_excl
    num_nb14_calc   => enefunc%num_nb14_calc
    exclusion_mask  => enefunc%exclusion_mask
    exclusion_mask1 => enefunc%exclusion_mask1
    nb14_qq_scale   => enefunc%nb14_qq_scale
    nb14_lj_scale   => enefunc%nb14_lj_scale
    dihe_scnb       => enefunc%dihe_scnb
    dihe_scee       => enefunc%dihe_scee

    ncell_local = domain%num_cell_local
    ncell       = domain%num_cell_local + domain%num_cell_boundary

    domain%max_num_atom = 0
    do i = 1, ncell
      domain%max_num_atom = max(domain%max_num_atom,domain%num_atom(i))
    end do

    ! initialization
    !
    num_nonb_excl(1:ncell_local)  = 0
    num_nb14_calc(1:ncell_local)  = 0

    domain%start_atom(1:ncell) = 0
    ij = natom(1)
    !ocl nosimd
    !dir$ novector
    do i = 2, ncell
      domain%start_atom(i) = domain%start_atom(i-1) + natom(i-1)
      ij = ij + natom(i)
    end do
    domain%num_atom_domain = ij

    ! exclude 1-2 interaction
    !
    !$omp parallel default(shared)                                     &
    !$omp private(id, i, ix, icel1, icel2, icel, i1, i2, i3, num_excl, &
    !$omp         k, num_nb14, ic, j, ih, list1, list2, ii, ij, fkind, &
    !$omp         index)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell_local, nthread
      k = natom(i)
      exclusion_mask1(1:k,1:k,i) = 1
      do i1 = 1, k
        exclusion_mask1(i1,i1,i) = 0
      end do
    end do
    do ij = id+1, maxcell_near, nthread
      i = cell_pair(1,ij)
      j = cell_pair(2,ij)
      i1 = max(natom(i),natom(j))
      exclusion_mask(1:i1,1:i1,ij) = 1
    end do
    !$omp barrier

    if (enefunc%excl_level > 0) then

      do i = id+1, ncell_local, nthread
        do ix = 1, nbond(i)

          fkind = bond_kind(ix,i)

          if (fkind == 0) then

            list1 = bond_list(1,ix,i)
            list2 = bond_list(2,ix,i)
            icel1 = id_g2l(1,list1)
            icel2 = id_g2l(1,list2)
            i1    = id_g2l(2,list1)
            i2    = id_g2l(2,list2)
            num_excl = num_nonb_excl(i) + 1
            num_nonb_excl(i) = num_excl
            if (domain%nonbond_kernel == NBK_Intel  .or. &
                domain%nonbond_kernel == NBK_Fugaku) then
              nonb_excl_list(1,num_excl,i) = start_atom(icel1) + i1
              nonb_excl_list(2,num_excl,i) = start_atom(icel2) + i2
            else
              nonb_excl_list(1,num_excl,i) = icel1
              nonb_excl_list(2,num_excl,i) = icel2
              nonb_excl_list(3,num_excl,i) = i1
              nonb_excl_list(4,num_excl,i) = i2
            end if

            if (icel1 /= icel2) then

              icel  = cell_pairlist2(icel1,icel2)
              if (icel1 < icel2) then
                exclusion_mask(i2,i1,icel) = 0
              else if (icel1 > icel2) then
                exclusion_mask(i1,i2,icel) = 0
              end if
  
            else

              if (i1 < i2) then
                exclusion_mask1(i2,i1,i) = 0
              else if (i1 > i2) then
                exclusion_mask1(i1,i2,i) = 0
              end if

            end if

          end if

        end do
      end do
   
      ! exclude constraint
      !
      if (constraint) then

        HGr_local       => constraints%HGr_local
        HGr_bond_list   => constraints%HGr_bond_list
        do icel = id+1, ncell_local, nthread
          do ic = 1, constraints%connect
            do j = 1, HGr_local(ic,icel)

              i1 = HGr_bond_list(1,j,ic,icel)
!ocl nosimd
              do ih = 1, ic
                i2 = HGr_bond_list(ih+1,j,ic,icel)
                num_excl = num_nonb_excl(icel) + 1
                num_nonb_excl(icel) = num_excl
                if (domain%nonbond_kernel == NBK_Intel  .or. &
                    domain%nonbond_kernel == NBK_Fugaku) then
                  nonb_excl_list(1,num_excl,icel) = start_atom(icel) + i1
                  nonb_excl_list(2,num_excl,icel) = start_atom(icel) + i2
                else
                  nonb_excl_list(1,num_excl,icel) = icel
                  nonb_excl_list(2,num_excl,icel) = icel
                  nonb_excl_list(3,num_excl,icel) = i1
                  nonb_excl_list(4,num_excl,icel) = i2
                end if
                exclusion_mask1(i2,i1,icel) = 0
              end do

            end do
          end do
        end do

      end if

    end if

    !$omp barrier

    ! exclude water
    !
    if (enefunc%excl_level > 1) then

      if (constraints%water_type == TIP4) then

        do icel = id+1, ncell_local, nthread
          do ic = 1, nwater(icel)

            index(1) = water_list(1,ic,icel)
            index(2) = water_list(2,ic,icel)
            index(3) = water_list(3,ic,icel)
            index(4) = water_list(4,ic,icel)

!ocl nosimd
            do i1 = 2, 3
              do i2 = i1+1, 4
                num_excl = num_nonb_excl(icel) + 1
                num_nonb_excl(icel) = num_excl
                if (domain%nonbond_kernel == NBK_Intel  .or. &
                    domain%nonbond_kernel == NBK_Fugaku) then
                  nonb_excl_list(1,num_excl,icel) = start_atom(icel) + index(i1)
                  nonb_excl_list(2,num_excl,icel) = start_atom(icel) + index(i2)
                else
                  nonb_excl_list(1,num_excl,icel) = icel
                  nonb_excl_list(2,num_excl,icel) = icel
                  nonb_excl_list(3,num_excl,icel) = index(i1)
                  nonb_excl_list(4,num_excl,icel) = index(i2)
                end if
                exclusion_mask1(index(i2),index(i1),icel) = 0
              end do
            end do
            i1 = 1
            do i2 = 2, 4
              exclusion_mask1(index(i2),index(i1),icel) = 0
            end do

          end do
        end do

      else

        do icel = id+1, ncell_local, nthread
          do ic = 1, nwater(icel)

            index(1) = water_list(1,ic,icel)
            index(2) = water_list(2,ic,icel)
            index(3) = water_list(3,ic,icel)
          
            num_excl = num_nonb_excl(icel) + 1
            num_nonb_excl(icel) = num_excl
            if (domain%nonbond_kernel == NBK_Intel  .or. &
                domain%nonbond_kernel == NBK_Fugaku) then
              nonb_excl_list(1,num_excl,icel) = start_atom(icel) + index(1)
              nonb_excl_list(2,num_excl,icel) = start_atom(icel) + index(2)
            else
              nonb_excl_list(1,num_excl,icel) = icel
              nonb_excl_list(2,num_excl,icel) = icel
              nonb_excl_list(3,num_excl,icel) = index(1)
              nonb_excl_list(4,num_excl,icel) = index(2)
            end if

            num_excl = num_nonb_excl(icel) + 1
            num_nonb_excl(icel) = num_excl
            if (domain%nonbond_kernel == NBK_Intel  .or. &
                domain%nonbond_kernel == NBK_Fugaku) then
              nonb_excl_list(1,num_excl,icel) = start_atom(icel) + index(1)
              nonb_excl_list(2,num_excl,icel) = start_atom(icel) + index(3)
            else
              nonb_excl_list(1,num_excl,icel) = icel
              nonb_excl_list(2,num_excl,icel) = icel
              nonb_excl_list(3,num_excl,icel) = index(1)
              nonb_excl_list(4,num_excl,icel) = index(3)
            end if

            num_excl = num_nonb_excl(icel) + 1
            num_nonb_excl(icel) = num_excl
            if (domain%nonbond_kernel == NBK_Intel  .or. &
                domain%nonbond_kernel == NBK_Fugaku) then
              nonb_excl_list(1,num_excl,icel) = start_atom(icel) + index(2)
              nonb_excl_list(2,num_excl,icel) = start_atom(icel) + index(3)
            else
              nonb_excl_list(1,num_excl,icel) = icel
              nonb_excl_list(2,num_excl,icel) = icel
              nonb_excl_list(3,num_excl,icel) = index(2)
              nonb_excl_list(4,num_excl,icel) = index(3)
            end if

            exclusion_mask1(index(2),index(1),icel) = 0
            exclusion_mask1(index(3),index(1),icel) = 0
            exclusion_mask1(index(3),index(2),icel) = 0

          end do
        end do

      end if

      ! exclude 1-3 interaction
      !
      do i = id+1, ncell_local, nthread
        do ix = 1, nangle(i)

          fkind = angl_kind(ix,i)

          if (fkind == 0) then

            list1 = angl_list(1,ix,i)
            list2 = angl_list(3,ix,i)
            icel1 = id_g2l(1,list1)
            icel2 = id_g2l(1,list2)
            i1    = id_g2l(2,list1)
            i2    = id_g2l(2,list2)

            if (icel1 /= icel2) then

              icel  = cell_pairlist2(icel1,icel2)

              if (icel1 < icel2) then

                if (exclusion_mask(i2,i1,icel)==1) then
                  if (abs(charge(i1,icel1)) > EPS .and. &
                      abs(charge(i2,icel2)) > EPS) then
                    num_excl = num_nonb_excl(i) + 1
                    num_nonb_excl(i) = num_excl
                    if (domain%nonbond_kernel == NBK_Intel  .or. &
                        domain%nonbond_kernel == NBK_Fugaku) then
                      nonb_excl_list(1,num_excl,i) = start_atom(icel1) + i1
                      nonb_excl_list(2,num_excl,i) = start_atom(icel2) + i2
                    else
                      nonb_excl_list(1,num_excl,i) = icel1
                      nonb_excl_list(2,num_excl,i) = icel2
                      nonb_excl_list(3,num_excl,i) = i1
                      nonb_excl_list(4,num_excl,i) = i2
                    end if
                  end if
                  exclusion_mask(i2,i1,icel) = 0
                end if

              else if (icel1 > icel2) then

                if (exclusion_mask(i1,i2,icel)==1) then
                  if (abs(charge(i1,icel1)) > EPS .and. &
                      abs(charge(i2,icel2)) > EPS) then
                    num_excl = num_nonb_excl(i) + 1
                    num_nonb_excl(i) = num_excl
                    if (domain%nonbond_kernel == NBK_Intel  .or. &
                        domain%nonbond_kernel == NBK_Fugaku) then
                      nonb_excl_list(1,num_excl,i) = start_atom(icel1) + i1
                      nonb_excl_list(2,num_excl,i) = start_atom(icel2) + i2
                    else
                      nonb_excl_list(1,num_excl,i) = icel1
                      nonb_excl_list(2,num_excl,i) = icel2
                      nonb_excl_list(3,num_excl,i) = i1
                      nonb_excl_list(4,num_excl,i) = i2
                    end if
                  end if
                  exclusion_mask(i1,i2,icel) = 0
                end if

              end if

            else

              if (i1 < i2) then

                if (exclusion_mask1(i2,i1,i)==1) then
                  if (abs(charge(i1,icel1)) > EPS .and. &
                      abs(charge(i2,icel2)) > EPS) then
                    num_excl = num_nonb_excl(i) + 1
                    num_nonb_excl(i) = num_excl
                    if (domain%nonbond_kernel == NBK_Intel  .or. &
                        domain%nonbond_kernel == NBK_Fugaku) then
                      nonb_excl_list(1,num_excl,i) = start_atom(icel1) + i1
                      nonb_excl_list(2,num_excl,i) = start_atom(icel2) + i2
                    else
                      nonb_excl_list(1,num_excl,i) = icel1
                      nonb_excl_list(2,num_excl,i) = icel2
                      nonb_excl_list(3,num_excl,i) = i1
                      nonb_excl_list(4,num_excl,i) = i2
                    end if
                  end if
                  exclusion_mask1(i2,i1,i) = 0
                end if

              else if (i1 > i2) then

                if (exclusion_mask1(i1,i2,i)==1) then
                  if (abs(charge(i1,icel1)) > EPS .and. &
                      abs(charge(i2,icel2)) > EPS) then
                    num_excl = num_nonb_excl(i) + 1
                    num_nonb_excl(i) = num_excl
                    if (domain%nonbond_kernel == NBK_Intel  .or. &
                        domain%nonbond_kernel == NBK_Fugaku) then
                      nonb_excl_list(1,num_excl,i) = start_atom(icel1) + i1
                      nonb_excl_list(2,num_excl,i) = start_atom(icel2) + i2
                    else
                      nonb_excl_list(1,num_excl,i) = icel1
                      nonb_excl_list(2,num_excl,i) = icel2
                      nonb_excl_list(3,num_excl,i) = i1
                      nonb_excl_list(4,num_excl,i) = i2
                    end if
                  end if
                  exclusion_mask1(i1,i2,i) = 0
                end if

              end if

            end if

          end if
        end do
      end do

    end if

    !$omp end parallel

    ! count 1-4 interaction
    !
    if (enefunc%excl_level > 2) then

      do ii = 1, 2

        if (ii == 2) then
          ndihedral => enefunc%num_rb_dihedral
          dihe_list => enefunc%rb_dihe_list
        end if

        !$omp parallel default(shared)                                     &
        !$omp private(id, i, ix, icel1, icel2, icel, i1, i2, list1, list2, &
        !$omp         num_nb14, fkind)
        !
#ifdef OMP
        id = omp_get_thread_num()
#else
        id = 0
#endif
        do i = id+1, ncell_local, nthread

          do ix = 1, ndihedral(i)

            fkind = dihe_kind(ix,i)

            if (fkind == 0) then

              list1 = dihe_list(1,ix,i)
              list2 = dihe_list(4,ix,i)
              icel1 = id_g2l(1,list1)
              icel2 = id_g2l(1,list2)
              i1    = id_g2l(2,list1)
              i2    = id_g2l(2,list2)

              if (icel1 /= icel2) then

                icel  = cell_pairlist2(icel1,icel2)

                if (icel1 < icel2) then

                  if (exclusion_mask(i2,i1,icel)==1) then
                    num_nb14 = num_nb14_calc(i) + 1
                    num_nb14_calc(i) = num_nb14
                    if (domain%nonbond_kernel == NBK_Intel  .or. &
                        domain%nonbond_kernel == NBK_Fugaku) then
                      nb14_calc_list(1,num_nb14,i) = start_atom(icel1) + i1
                      nb14_calc_list(2,num_nb14,i) = start_atom(icel2) + i2
                    else
                      nb14_calc_list(1,num_nb14,i) = icel1
                      nb14_calc_list(2,num_nb14,i) = icel2
                      nb14_calc_list(3,num_nb14,i) = i1
                      nb14_calc_list(4,num_nb14,i) = i2
                    end if
                    sc_calc_list(num_nb14,i)     = &
                      int(enefunc%dihe_periodicity(ix,i) &
                          / enefunc%notation_14types)
                    exclusion_mask(i2,i1,icel) = 0
                  end if

                else if (icel1 > icel2) then

                  if (exclusion_mask(i1,i2,icel)==1) then
                    num_nb14 = num_nb14_calc(i) + 1
                    num_nb14_calc(i) = num_nb14
                    if (domain%nonbond_kernel == NBK_Intel  .or. &
                        domain%nonbond_kernel == NBK_Fugaku) then
                      nb14_calc_list(1,num_nb14,i) = start_atom(icel1) + i1
                      nb14_calc_list(2,num_nb14,i) = start_atom(icel2) + i2
                    else
                      nb14_calc_list(1,num_nb14,i) = icel1
                      nb14_calc_list(2,num_nb14,i) = icel2
                      nb14_calc_list(3,num_nb14,i) = i1
                      nb14_calc_list(4,num_nb14,i) = i2
                    end if
                    sc_calc_list(num_nb14,i)     = &
                      int(enefunc%dihe_periodicity(ix,i) &
                          / enefunc%notation_14types)
                    exclusion_mask(i1,i2,icel) = 0
                  end if

                end if

              else

                if (i1 < i2) then

                  if (exclusion_mask1(i2,i1,i)==1) then
                    num_nb14 = num_nb14_calc(i) + 1
                    num_nb14_calc(i) = num_nb14
                    if (domain%nonbond_kernel == NBK_Intel  .or. &
                        domain%nonbond_kernel == NBK_Fugaku) then
                      nb14_calc_list(1,num_nb14,i) = start_atom(icel1) + i1
                      nb14_calc_list(2,num_nb14,i) = start_atom(icel2) + i2
                    else
                      nb14_calc_list(1,num_nb14,i) = icel1
                      nb14_calc_list(2,num_nb14,i) = icel2
                      nb14_calc_list(3,num_nb14,i) = i1
                      nb14_calc_list(4,num_nb14,i) = i2
                    end if
                    sc_calc_list(num_nb14,i)     = &
                      int(enefunc%dihe_periodicity(ix,i) &
                          / enefunc%notation_14types)
                    exclusion_mask1(i2,i1,i) = 0
                  end if

                else if (i1 > i2) then

                  if (exclusion_mask1(i1,i2,i)==1) then
                    num_nb14 = num_nb14_calc(i) + 1
                    num_nb14_calc(i) = num_nb14
                    if (domain%nonbond_kernel == NBK_Intel  .or. &
                        domain%nonbond_kernel == NBK_Fugaku) then
                      nb14_calc_list(1,num_nb14,i) = start_atom(icel1) + i1
                      nb14_calc_list(2,num_nb14,i) = start_atom(icel2) + i2
                    else
                      nb14_calc_list(1,num_nb14,i) = icel1
                      nb14_calc_list(2,num_nb14,i) = icel2
                      nb14_calc_list(3,num_nb14,i) = i1
                      nb14_calc_list(4,num_nb14,i) = i2
                    end if
                    sc_calc_list(num_nb14,i)     = &
                      int(enefunc%dihe_periodicity(ix,i) &
                          / enefunc%notation_14types)
                    exclusion_mask1(i1,i2,i) = 0
                  end if

                end if

              end if

            end if
          end do
        end do
        !$omp end parallel

      end do

    end if

    ! scnb/fudge_lj & scee/fudge_qq
    !
    !$omp parallel default(shared)            &
    !$omp private(id, i, ix, list1, list2)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    if (enefunc%forcefield == ForcefieldAMBER) then
      do i = id+1, ncell_local, nthread
        do ix = 1, num_nb14_calc(i)
          list1 = sc_calc_list(ix,i)
          nb14_lj_scale(ix,i) = dihe_scnb(list1)
          nb14_qq_scale(ix,i) = dihe_scee(list1)
        end do
      end do
    end if
    if (enefunc%forcefield == ForcefieldGROAMBER .or. &
        enefunc%forcefield == ForcefieldGROMARTINI) then
      do i = id+1, ncell_local, nthread
        do ix = 1, num_nb14_calc(i)
          list1 = sc_calc_list(ix,i)
          nb14_lj_scale(ix,i) = enefunc%fudge_lj
          nb14_qq_scale(ix,i) = enefunc%fudge_qq
        end do
      end do
    end if
    !$omp end parallel

    ! Check the total number of exclusion list
    !
    if (first) then

      found1 = 0
      found2 = 0

      do icel = 1, ncell_local
        found1 = found1 + num_nonb_excl(icel)
        found2 = found2 + num_nb14_calc(icel)
      end do

#ifdef HAVE_MPI_GENESIS
      call mpi_reduce(found1, enefunc%num_excl_all, 1, mpi_integer, mpi_sum, &
                      0, mpi_comm_country, ierror)
      call mpi_reduce(found2, enefunc%num_nb14_all, 1, mpi_integer, mpi_sum, &
                      0, mpi_comm_country, ierror)
#else
    enefunc%num_excl_all = found1
    enefunc%num_nb14_all = found2
#endif
    end if

    return

  end subroutine count_nonb_excl

  subroutine check_pbc(box_size, dij, pbc_int)

    real(wp),         intent(in)    :: box_size(:)
    real(wp),         intent(inout) :: dij(:)
    integer,          intent(inout) :: pbc_int

    integer                  :: i, j, k

    if (dij(1) > box_size(1)/2.0_dp) then
      i = 0
      dij(1) = dij(1) - box_size(1)
    else if (dij(1) < -box_size(1)/2.0_dp) then
      i = 2
      dij(1) = dij(1) + box_size(1)
    else
      i = 1
    end if

    if (dij(2) > box_size(2)/2.0_dp) then
      j = 0
      dij(2) = dij(2) - box_size(2)
    else if (dij(2) < -box_size(2)/2.0_dp) then
      j = 2
      dij(2) = dij(2) + box_size(2)
    else
      j = 1
    end if

    if (dij(3) > box_size(3)/2.0_dp) then
      k = 0
      dij(3) = dij(3) - box_size(3)
    else if (dij(3) < -box_size(3)/2.0_dp) then
      k = 2
      dij(3) = dij(3) + box_size(3)
    else
      k = 1
    end if

    pbc_int = i + j*3 + k*9

    return

  end subroutine check_pbc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    count_nonb_excl_fep
  !> @brief        exclude 1-2, 1-3 interactions and constraints for FEP
  !! @authors      HO
  !! @param[in]    first       : flag for first call or not
  !! @param[in]    constraint  : flag for constraint usage
  !! @param[in]    constraints : constraints information   
  !! @param[inout] domain      : structure of domain
  !! @param[inout] enefunc     : structure of enefunc
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine count_nonb_excl_fep(first, constraint, constraints, domain, enefunc)

    ! formal arguments
    logical,                     intent(in)    :: first
    logical,                     intent(in)    :: constraint
    type(s_constraints), target, intent(in)    :: constraints
    type(s_domain),      target, intent(inout) :: domain
    type(s_enefunc),     target, intent(inout) :: enefunc

    ! local variables
    integer                  :: ncell, ncell_local, i, ii, ix, k, i1, i2, i3
    integer                  :: icel, icel1, icel2
    integer                  :: ic, j, ih, ij, index(4)
    integer                  :: num_excl, num_nb14, id, omp_get_thread_num
    integer                  :: found1, found2
    integer                  :: fkind
    integer                  :: list1, list2, list3, step

    integer,         pointer :: natom(:), nwater(:), start_atom(:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: water_list(:,:,:)
    integer,         pointer :: nbond(:), bond_list(:,:,:)
    integer,         pointer :: nangle(:), angl_list(:,:,:)
    integer(1),      pointer :: bond_kind(:,:), angl_kind(:,:), dihe_kind(:,:)
    integer,         pointer :: ndihedral(:), dihe_list(:,:,:)
    integer,         pointer :: cell_pairlist2(:,:)
    integer(int2),   pointer :: cell_pairlist1(:,:), cell_pair(:,:)
    integer,         pointer :: nonb_excl_list(:,:,:)
    integer,         pointer :: nb14_calc_list(:,:,:)
    integer(1),      pointer :: exclusion_mask(:,:,:), exclusion_mask1(:,:,:)
    integer,         pointer :: sc_calc_list(:,:)
    integer,         pointer :: num_nonb_excl(:), num_nb14_calc(:)
    integer,         pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)
    real(wp),        pointer :: nb14_qq_scale(:,:), nb14_lj_scale(:,:)
    real(wp),        pointer :: dihe_scnb(:), dihe_scee(:)
    real(wp),        pointer :: charge(:,:)

    ! FEP
    integer                  :: fg1, fg2
    integer                  :: num_nb14_fep
    integer                  :: found_fep1, found_fep2
    integer(1)               :: flag_bond, flag_angl, flag_dihe
    integer,         pointer :: nb14_calc_list_fep(:,:,:)
    integer,         pointer :: sc_calc_list_fep(:,:)
    integer,         pointer :: num_nb14_calc_fep(:)
    real(wp),        pointer :: nb14_qq_scale_fep(:,:), nb14_lj_scale_fep(:,:)

    natom           => domain%num_atom
    nwater          => domain%num_water
    start_atom      => domain%start_atom
    water_list      => domain%water_list
    id_g2l          => domain%id_g2l
    cell_pair       => domain%cell_pair
    cell_pairlist1  => domain%cell_pairlist1
    cell_pairlist2  => domain%cell_pairlist2
    charge          => domain%charge

    nbond           => enefunc%num_bond
    nangle          => enefunc%num_angle
    ndihedral       => enefunc%num_dihedral
    bond_list       => enefunc%bond_list
    bond_kind       => enefunc%bond_kind
    angl_list       => enefunc%angle_list
    angl_kind       => enefunc%angle_kind
    dihe_list       => enefunc%dihe_list
    dihe_kind       => enefunc%dihe_kind
    nonb_excl_list  => enefunc%nonb_excl_list
    nb14_calc_list  => enefunc%nb14_calc_list
    sc_calc_list    => enefunc%sc_calc_list
    num_nonb_excl   => enefunc%num_nonb_excl
    num_nb14_calc   => enefunc%num_nb14_calc
    exclusion_mask  => enefunc%exclusion_mask
    exclusion_mask1 => enefunc%exclusion_mask1
    nb14_qq_scale   => enefunc%nb14_qq_scale
    nb14_lj_scale   => enefunc%nb14_lj_scale
    dihe_scnb       => enefunc%dihe_scnb
    dihe_scee       => enefunc%dihe_scee

    ! FEP
    nb14_calc_list_fep  => enefunc%nb14_calc_list_fep
    sc_calc_list_fep    => enefunc%sc_calc_list_fep
    num_nb14_calc_fep   => enefunc%num_nb14_calc_fep
    nb14_qq_scale_fep   => enefunc%nb14_qq_scale_fep
    nb14_lj_scale_fep   => enefunc%nb14_lj_scale_fep

    ncell_local = domain%num_cell_local
    ncell       = domain%num_cell_local + domain%num_cell_boundary

    domain%max_num_atom = 0
    do i = 1, ncell
      domain%max_num_atom = max(domain%max_num_atom,domain%num_atom(i))
    end do

    ! initialization
    !
    num_nonb_excl(1:ncell_local)  = 0
    num_nb14_calc(1:ncell_local)  = 0
    ! FEP
    num_nb14_calc_fep(1:ncell_local)  = 0

    domain%start_atom(1:ncell) = 0
    ij = natom(1)
    !ocl nosimd
    !dir$ novector
    do i = 2, ncell
      domain%start_atom(i) = domain%start_atom(i-1) + natom(i-1)
      ij = ij + natom(i)
    end do
    domain%num_atom_domain = ij

    ! exclude 1-2 interaction
    !
    !$omp parallel default(shared)                                     &
    !$omp private(id, i, ix, icel1, icel2, icel, i1, i2, i3, num_excl, &
    !$omp         k, num_nb14, ic, j, ih, list1, list2, ii, ij, fkind, &
    !$omp         index, flag_bond, flag_angl)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell_local, nthread
      k = natom(i)
      exclusion_mask1(1:k,1:k,i) = 1
      do i1 = 1, k
        exclusion_mask1(i1,i1,i) = 0
      end do
    end do
    do ij = id+1, maxcell_near, nthread
      i = cell_pairlist1(1,ij)
      j = cell_pairlist1(2,ij)
      i1 = max(natom(i),natom(j))
      exclusion_mask(1:i1,1:i1,ij) = 1
    end do
    !$omp barrier

    if (enefunc%excl_level > 0) then

      do i = id+1, ncell_local, nthread
        do ix = 1, nbond(i)

          fkind = bond_kind(ix,i)

          ! FEP
          flag_bond = enefunc%bond_singleB(ix,i)

          if ((fkind == 0) .and. (flag_bond == 0)) then

            list1 = bond_list(1,ix,i)
            list2 = bond_list(2,ix,i)
            icel1 = id_g2l(1,list1)
            icel2 = id_g2l(1,list2)
            i1    = id_g2l(2,list1)
            i2    = id_g2l(2,list2)
            num_excl = num_nonb_excl(i) + 1
            num_nonb_excl(i) = num_excl
            if (domain%nonbond_kernel == NBK_Intel  .or. &
                domain%nonbond_kernel == NBK_Fugaku) then
              nonb_excl_list(1,num_excl,i) = start_atom(icel1) + i1
              nonb_excl_list(2,num_excl,i) = start_atom(icel2) + i2
            else
              nonb_excl_list(1,num_excl,i) = icel1
              nonb_excl_list(2,num_excl,i) = icel2
              nonb_excl_list(3,num_excl,i) = i1
              nonb_excl_list(4,num_excl,i) = i2
            end if

            if (icel1 /= icel2) then

              icel  = cell_pairlist2(icel1,icel2)
              if (icel1 < icel2) then
                exclusion_mask(i2,i1,icel) = 0
              else if (icel1 > icel2) then
                exclusion_mask(i1,i2,icel) = 0
              end if
  
            else

              if (i1 < i2) then
                exclusion_mask1(i2,i1,i) = 0
              else if (i1 > i2) then
                exclusion_mask1(i1,i2,i) = 0
              end if

            end if

          end if

        end do
      end do
   
      ! exclude constraint
      !
      if (constraint) then

        HGr_local       => constraints%HGr_local
        HGr_bond_list   => constraints%HGr_bond_list
        do icel = id+1, ncell_local, nthread
          do ic = 1, constraints%connect
            do j = 1, HGr_local(ic,icel)

              i1 = HGr_bond_list(1,j,ic,icel)
!ocl nosimd
              do ih = 1, ic
                i2 = HGr_bond_list(ih+1,j,ic,icel)
                num_excl = num_nonb_excl(icel) + 1
                num_nonb_excl(icel) = num_excl
                if (domain%nonbond_kernel == NBK_Intel  .or. &
                    domain%nonbond_kernel == NBK_Fugaku) then
                  nonb_excl_list(1,num_excl,icel) = start_atom(icel) + i1
                  nonb_excl_list(2,num_excl,icel) = start_atom(icel) + i2
                else
                  nonb_excl_list(1,num_excl,icel) = icel
                  nonb_excl_list(2,num_excl,icel) = icel
                  nonb_excl_list(3,num_excl,icel) = i1
                  nonb_excl_list(4,num_excl,icel) = i2
                end if
                exclusion_mask1(i2,i1,icel) = 0
              end do

            end do
          end do
        end do

      end if

    end if

    !$omp barrier

    ! exclude water
    !
    if (enefunc%excl_level > 1) then

      if (constraints%water_type == TIP4) then

        do icel = id+1, ncell_local, nthread
          do ic = 1, nwater(icel)

            index(1) = water_list(1,ic,icel)
            index(2) = water_list(2,ic,icel)
            index(3) = water_list(3,ic,icel)
            index(4) = water_list(4,ic,icel)

!ocl nosimd
            do i1 = 2, 3
              do i2 = i1+1, 4
                num_excl = num_nonb_excl(icel) + 1
                num_nonb_excl(icel) = num_excl
                if (domain%nonbond_kernel == NBK_Intel  .or. &
                    domain%nonbond_kernel == NBK_Fugaku) then
                  nonb_excl_list(1,num_excl,icel) = start_atom(icel) + index(i1)
                  nonb_excl_list(2,num_excl,icel) = start_atom(icel) + index(i2)
                else
                  nonb_excl_list(1,num_excl,icel) = icel
                  nonb_excl_list(2,num_excl,icel) = icel
                  nonb_excl_list(3,num_excl,icel) = index(i1)
                  nonb_excl_list(4,num_excl,icel) = index(i2)
                end if
                exclusion_mask1(index(i2),index(i1),icel) = 0
              end do
            end do
            i1 = 1
            do i2 = 2, 4
              exclusion_mask1(index(i2),index(i1),icel) = 0
            end do

          end do
        end do

      else

        do icel = id+1, ncell_local, nthread
          do ic = 1, nwater(icel)

            index(1) = water_list(1,ic,icel)
            index(2) = water_list(2,ic,icel)
            index(3) = water_list(3,ic,icel)
          
            num_excl = num_nonb_excl(icel) + 1
            num_nonb_excl(icel) = num_excl
            if (domain%nonbond_kernel == NBK_Intel  .or. &
                domain%nonbond_kernel == NBK_Fugaku) then
              nonb_excl_list(1,num_excl,icel) = start_atom(icel) + index(1)
              nonb_excl_list(2,num_excl,icel) = start_atom(icel) + index(2)
            else
              nonb_excl_list(1,num_excl,icel) = icel
              nonb_excl_list(2,num_excl,icel) = icel
              nonb_excl_list(3,num_excl,icel) = index(1)
              nonb_excl_list(4,num_excl,icel) = index(2)
            end if

            num_excl = num_nonb_excl(icel) + 1
            num_nonb_excl(icel) = num_excl
            if (domain%nonbond_kernel == NBK_Intel  .or. &
                domain%nonbond_kernel == NBK_Fugaku) then
              nonb_excl_list(1,num_excl,icel) = start_atom(icel) + index(1)
              nonb_excl_list(2,num_excl,icel) = start_atom(icel) + index(3)
            else
              nonb_excl_list(1,num_excl,icel) = icel
              nonb_excl_list(2,num_excl,icel) = icel
              nonb_excl_list(3,num_excl,icel) = index(1)
              nonb_excl_list(4,num_excl,icel) = index(3)
            end if

            num_excl = num_nonb_excl(icel) + 1
            num_nonb_excl(icel) = num_excl
            if (domain%nonbond_kernel == NBK_Intel  .or. &
                domain%nonbond_kernel == NBK_Fugaku) then
              nonb_excl_list(1,num_excl,icel) = start_atom(icel) + index(2)
              nonb_excl_list(2,num_excl,icel) = start_atom(icel) + index(3)
            else
              nonb_excl_list(1,num_excl,icel) = icel
              nonb_excl_list(2,num_excl,icel) = icel
              nonb_excl_list(3,num_excl,icel) = index(2)
              nonb_excl_list(4,num_excl,icel) = index(3)
            end if

            exclusion_mask1(index(2),index(1),icel) = 0
            exclusion_mask1(index(3),index(1),icel) = 0
            exclusion_mask1(index(3),index(2),icel) = 0

          end do
        end do

      end if

      ! exclude 1-3 interaction
      !
      do i = id+1, ncell_local, nthread
        do ix = 1, nangle(i)

          fkind = angl_kind(ix,i)

          ! FEP
          flag_angl = enefunc%angl_singleB(ix,i)

          if ((fkind == 0) .and. (flag_angl == 0)) then

            list1 = angl_list(1,ix,i)
            list2 = angl_list(3,ix,i)
            icel1 = id_g2l(1,list1)
            icel2 = id_g2l(1,list2)
            i1    = id_g2l(2,list1)
            i2    = id_g2l(2,list2)

            if (icel1 /= icel2) then

              icel  = cell_pairlist2(icel1,icel2)

              if (icel1 < icel2) then

                if (exclusion_mask(i2,i1,icel)==1) then
!                  if (abs(charge(i1,icel1)) > EPS .and. &
!                      abs(charge(i2,icel2)) > EPS) then
                    num_excl = num_nonb_excl(i) + 1
                    num_nonb_excl(i) = num_excl
                    if (domain%nonbond_kernel == NBK_Intel  .or. &
                        domain%nonbond_kernel == NBK_Fugaku) then
                      nonb_excl_list(1,num_excl,i) = start_atom(icel1) + i1
                      nonb_excl_list(2,num_excl,i) = start_atom(icel2) + i2
                    else
                      nonb_excl_list(1,num_excl,i) = icel1
                      nonb_excl_list(2,num_excl,i) = icel2
                      nonb_excl_list(3,num_excl,i) = i1
                      nonb_excl_list(4,num_excl,i) = i2
                    end if
!                  end if
                  exclusion_mask(i2,i1,icel) = 0
                end if

              else if (icel1 > icel2) then

                if (exclusion_mask(i1,i2,icel)==1) then
!                  if (abs(charge(i1,icel1)) > EPS .and. &
!                      abs(charge(i2,icel2)) > EPS) then
                    num_excl = num_nonb_excl(i) + 1
                    num_nonb_excl(i) = num_excl
                    if (domain%nonbond_kernel == NBK_Intel  .or. &
                        domain%nonbond_kernel == NBK_Fugaku) then
                      nonb_excl_list(1,num_excl,i) = start_atom(icel1) + i1
                      nonb_excl_list(2,num_excl,i) = start_atom(icel2) + i2
                    else
                      nonb_excl_list(1,num_excl,i) = icel1
                      nonb_excl_list(2,num_excl,i) = icel2
                      nonb_excl_list(3,num_excl,i) = i1
                      nonb_excl_list(4,num_excl,i) = i2
                    end if
!                  end if
                  exclusion_mask(i1,i2,icel) = 0
                end if

              end if

            else

              if (i1 < i2) then

                if (exclusion_mask1(i2,i1,i)==1) then
!                  if (abs(charge(i1,icel1)) > EPS .and. &
!                      abs(charge(i2,icel2)) > EPS) then
                    num_excl = num_nonb_excl(i) + 1
                    num_nonb_excl(i) = num_excl
                    if (domain%nonbond_kernel == NBK_Intel  .or. &
                        domain%nonbond_kernel == NBK_Fugaku) then
                      nonb_excl_list(1,num_excl,i) = start_atom(icel1) + i1
                      nonb_excl_list(2,num_excl,i) = start_atom(icel2) + i2
                    else
                      nonb_excl_list(1,num_excl,i) = icel1
                      nonb_excl_list(2,num_excl,i) = icel2
                      nonb_excl_list(3,num_excl,i) = i1
                      nonb_excl_list(4,num_excl,i) = i2
                    end if
!                  end if
                  exclusion_mask1(i2,i1,i) = 0
                end if

              else if (i1 > i2) then

                if (exclusion_mask1(i1,i2,i)==1) then
!                  if (abs(charge(i1,icel1)) > EPS .and. &
!                      abs(charge(i2,icel2)) > EPS) then
                    num_excl = num_nonb_excl(i) + 1
                    num_nonb_excl(i) = num_excl
                    if (domain%nonbond_kernel == NBK_Intel  .or. &
                        domain%nonbond_kernel == NBK_Fugaku) then
                      nonb_excl_list(1,num_excl,i) = start_atom(icel1) + i1
                      nonb_excl_list(2,num_excl,i) = start_atom(icel2) + i2
                    else
                      nonb_excl_list(1,num_excl,i) = icel1
                      nonb_excl_list(2,num_excl,i) = icel2
                      nonb_excl_list(3,num_excl,i) = i1
                      nonb_excl_list(4,num_excl,i) = i2
                    end if
!                  end if
                  exclusion_mask1(i1,i2,i) = 0
                end if

              end if

            end if

          end if
        end do
      end do

    end if
    !$omp end parallel

    ! count 1-4 interaction
    !
    if (enefunc%excl_level > 2) then

      do ii = 1, 2

        if (ii == 2) then
          ndihedral => enefunc%num_rb_dihedral
          dihe_list => enefunc%rb_dihe_list
        end if

        !$omp parallel default(shared)                                     &
        !$omp private(id, i, ix, icel1, icel2, icel, i1, i2, list1, list2, &
        !$omp         num_nb14, fkind, fg1, fg2, num_nb14_fep, flag_dihe)
        !
#ifdef OMP
        id = omp_get_thread_num()
#else
        id = 0
#endif
        do i = id+1, ncell_local, nthread

          do ix = 1, ndihedral(i)

            fkind = dihe_kind(ix,i)

            ! FEP
            if (ii == 2) then
              flag_dihe = enefunc%rb_dihe_singleB(ix,i)
            else
              flag_dihe = enefunc%dihe_singleB(ix,i)
            end if


            if ((fkind == 0) .and. (flag_dihe==0)) then

              list1 = dihe_list(1,ix,i)
              list2 = dihe_list(4,ix,i)
              icel1 = id_g2l(1,list1)
              icel2 = id_g2l(1,list2)
              i1    = id_g2l(2,list1)
              i2    = id_g2l(2,list2)

              ! FEP
              fg1 = domain%fepgrp(i1,icel1)
              fg2 = domain%fepgrp(i2,icel2)

              if (icel1 /= icel2) then

                icel  = cell_pairlist2(icel1,icel2)

                if (icel1 < icel2) then

                  if (exclusion_mask(i2,i1,icel)==1) then

                    if (enefunc%fepgrp_nonb(fg1,fg2) == 5) then
                      ! FEP: non-perturbed
                      num_nb14 = num_nb14_calc(i) + 1
                      num_nb14_calc(i) = num_nb14
                      if (domain%nonbond_kernel == NBK_Intel  .or. &
                          domain%nonbond_kernel == NBK_Fugaku) then
                        nb14_calc_list(1,num_nb14,i) = start_atom(icel1) + i1
                        nb14_calc_list(2,num_nb14,i) = start_atom(icel2) + i2
                      else
                        nb14_calc_list(1,num_nb14,i) = icel1
                        nb14_calc_list(2,num_nb14,i) = icel2
                        nb14_calc_list(3,num_nb14,i) = i1
                        nb14_calc_list(4,num_nb14,i) = i2
                      end if
                      sc_calc_list(num_nb14,i)     = &
                        int(enefunc%dihe_periodicity(ix,i) &
                            / enefunc%notation_14types)
                      exclusion_mask(i2,i1,icel) = 0
                    else if (enefunc%fepgrp_nonb(fg1,fg2) /= 0) then
                      ! FEP: perturbed
                      num_nb14_fep = num_nb14_calc_fep(i) + 1
                      num_nb14_calc_fep(i) = num_nb14_fep
                      nb14_calc_list_fep(1,num_nb14_fep,i) = icel1
                      nb14_calc_list_fep(2,num_nb14_fep,i) = icel2
                      nb14_calc_list_fep(3,num_nb14_fep,i) = i1
                      nb14_calc_list_fep(4,num_nb14_fep,i) = i2
                      sc_calc_list_fep(num_nb14_fep,i)     = int(enefunc%dihe_periodicity(ix,i)/enefunc%notation_14types)
                      exclusion_mask(i2,i1,icel) = 0
                    else
                      exclusion_mask(i2,i1,icel) = 0
                    end if

                  end if

                else if (icel1 > icel2) then

                  if (exclusion_mask(i1,i2,icel)==1) then

                    if (enefunc%fepgrp_nonb(fg1,fg2) == 5) then
                      ! FEP: non-perturbed
                      num_nb14 = num_nb14_calc(i) + 1
                      num_nb14_calc(i) = num_nb14
                      if (domain%nonbond_kernel == NBK_Intel  .or. &
                          domain%nonbond_kernel == NBK_Fugaku) then
                        nb14_calc_list(1,num_nb14,i) = start_atom(icel1) + i1
                        nb14_calc_list(2,num_nb14,i) = start_atom(icel2) + i2
                      else
                        nb14_calc_list(1,num_nb14,i) = icel1
                        nb14_calc_list(2,num_nb14,i) = icel2
                        nb14_calc_list(3,num_nb14,i) = i1
                        nb14_calc_list(4,num_nb14,i) = i2
                      end if
                      sc_calc_list(num_nb14,i)     = &
                        int(enefunc%dihe_periodicity(ix,i) &
                            / enefunc%notation_14types)
                      exclusion_mask(i1,i2,icel) = 0
                    else if (enefunc%fepgrp_nonb(fg1,fg2) /= 0) then
                      ! FEP: perturbed
                      num_nb14_fep = num_nb14_calc_fep(i) + 1
                      num_nb14_calc_fep(i) = num_nb14_fep
                      nb14_calc_list_fep(1,num_nb14_fep,i) = icel1
                      nb14_calc_list_fep(2,num_nb14_fep,i) = icel2
                      nb14_calc_list_fep(3,num_nb14_fep,i) = i1
                      nb14_calc_list_fep(4,num_nb14_fep,i) = i2
                      sc_calc_list_fep(num_nb14_fep,i)     = int(enefunc%dihe_periodicity(ix,i)/enefunc%notation_14types)
                      exclusion_mask(i1,i2,icel) = 0
                    else
                      exclusion_mask(i1,i2,icel) = 0
                    end if

                  end if

                end if

              else

                if (i1 < i2) then

                  if (exclusion_mask1(i2,i1,i)==1) then

                    if (enefunc%fepgrp_nonb(fg1,fg2) == 5) then
                      ! FEP: non-perturbed
                      num_nb14 = num_nb14_calc(i) + 1
                      num_nb14_calc(i) = num_nb14
                      if (domain%nonbond_kernel == NBK_Intel  .or. &
                          domain%nonbond_kernel == NBK_Fugaku) then
                        nb14_calc_list(1,num_nb14,i) = start_atom(icel1) + i1
                        nb14_calc_list(2,num_nb14,i) = start_atom(icel2) + i2
                      else
                        nb14_calc_list(1,num_nb14,i) = icel1
                        nb14_calc_list(2,num_nb14,i) = icel2
                        nb14_calc_list(3,num_nb14,i) = i1
                        nb14_calc_list(4,num_nb14,i) = i2
                      end if
                      sc_calc_list(num_nb14,i)     = &
                        int(enefunc%dihe_periodicity(ix,i) &
                            / enefunc%notation_14types)
                      exclusion_mask1(i2,i1,i) = 0
                    else if (enefunc%fepgrp_nonb(fg1,fg2) /= 0) then
                      ! FEP: perturbed
                      num_nb14_fep = num_nb14_calc_fep(i) + 1
                      num_nb14_calc_fep(i) = num_nb14_fep
                      nb14_calc_list_fep(1,num_nb14_fep,i) = icel1
                      nb14_calc_list_fep(2,num_nb14_fep,i) = icel2
                      nb14_calc_list_fep(3,num_nb14_fep,i) = i1
                      nb14_calc_list_fep(4,num_nb14_fep,i) = i2
                      sc_calc_list_fep(num_nb14_fep,i)     = int(enefunc%dihe_periodicity(ix,i)/enefunc%notation_14types)
                      exclusion_mask1(i2,i1,i) = 0
                    else
                      exclusion_mask1(i2,i1,i) = 0
                    end if

                  end if

                else if (i1 > i2) then

                  if (exclusion_mask1(i1,i2,i)==1) then

                    if (enefunc%fepgrp_nonb(fg1,fg2) == 5) then
                      ! FEP: non-perturbed
                      num_nb14 = num_nb14_calc(i) + 1
                      num_nb14_calc(i) = num_nb14
                      if (domain%nonbond_kernel == NBK_Intel  .or. &
                          domain%nonbond_kernel == NBK_Fugaku) then
                        nb14_calc_list(1,num_nb14,i) = start_atom(icel1) + i1
                        nb14_calc_list(2,num_nb14,i) = start_atom(icel2) + i2
                      else
                        nb14_calc_list(1,num_nb14,i) = icel1
                        nb14_calc_list(2,num_nb14,i) = icel2
                        nb14_calc_list(3,num_nb14,i) = i1
                        nb14_calc_list(4,num_nb14,i) = i2
                      end if
                      sc_calc_list(num_nb14,i)     = &
                        int(enefunc%dihe_periodicity(ix,i) &
                            / enefunc%notation_14types)
                      exclusion_mask1(i1,i2,i) = 0
                    else if (enefunc%fepgrp_nonb(fg1,fg2) /= 0) then
                      ! FEP: perturbed
                      num_nb14_fep = num_nb14_calc_fep(i) + 1
                      num_nb14_calc_fep(i) = num_nb14_fep
                      nb14_calc_list_fep(1,num_nb14_fep,i) = icel1
                      nb14_calc_list_fep(2,num_nb14_fep,i) = icel2
                      nb14_calc_list_fep(3,num_nb14_fep,i) = i1
                      nb14_calc_list_fep(4,num_nb14_fep,i) = i2
                      sc_calc_list_fep(num_nb14_fep,i)     = int(enefunc%dihe_periodicity(ix,i)/enefunc%notation_14types)
                      exclusion_mask1(i1,i2,i) = 0
                    else
                      exclusion_mask1(i1,i2,i) = 0
                    end if

                  end if

                end if

              end if

            end if
          end do
        end do
        !$omp end parallel

      end do

    end if

    ! scnb/fudge_lj & scee/fudge_qq
    !
    !$omp parallel default(shared)            &
    !$omp private(id, i, ix, list1, list2)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    if (enefunc%forcefield == ForcefieldAMBER) then
      do i = id+1, ncell_local, nthread
        do ix = 1, num_nb14_calc(i)
          list1 = sc_calc_list(ix,i)
          nb14_lj_scale(ix,i) = dihe_scnb(list1)
          nb14_qq_scale(ix,i) = dihe_scee(list1)
        end do
        ! FEP
        do ix = 1, num_nb14_calc_fep(i)
          list1 = sc_calc_list_fep(ix,i)
          nb14_lj_scale_fep(ix,i) = dihe_scnb(list1)
          nb14_qq_scale_fep(ix,i) = dihe_scee(list1)
        end do
      end do
    end if
    if (enefunc%forcefield == ForcefieldGROAMBER .or. &
        enefunc%forcefield == ForcefieldGROMARTINI) then
      do i = id+1, ncell_local, nthread
        do ix = 1, num_nb14_calc(i)
          list1 = sc_calc_list(ix,i)
          nb14_lj_scale(ix,i) = enefunc%fudge_lj
          nb14_qq_scale(ix,i) = enefunc%fudge_qq
        end do
        ! FEP
        do ix = 1, num_nb14_calc_fep(i)
          list1 = sc_calc_list_fep(ix,i)
          nb14_lj_scale_fep(ix,i) = enefunc%fudge_lj
          nb14_qq_scale_fep(ix,i) = enefunc%fudge_qq
        end do
      end do
    end if
    !$omp end parallel

    ! Check the total number of exclusion list
    !
    if (first) then

      found1 = 0
      found2 = 0
      ! FEP
      found_fep2 = 0

      do icel = 1, ncell_local
        found1 = found1 + num_nonb_excl(icel)
        found2 = found2 + num_nb14_calc(icel)
        ! FEP
        found_fep2 = found_fep2 + num_nb14_calc_fep(icel)
      end do

#ifdef HAVE_MPI_GENESIS
      call mpi_reduce(found1, enefunc%num_excl_all, 1, mpi_integer, mpi_sum, &
                      0, mpi_comm_country, ierror)
      call mpi_reduce(found2, enefunc%num_nb14_all, 1, mpi_integer, mpi_sum, &
                      0, mpi_comm_country, ierror)
      ! FEP
      call mpi_reduce(found_fep2, enefunc%num_nb14_all_fep, 1, mpi_integer, mpi_sum, &
                      0, mpi_comm_country, ierror)
#else
    enefunc%num_excl_all = found1
    enefunc%num_nb14_all = found2
    ! FEP
    enefunc%num_nb14_all_fep = found_fep2
#endif
    end if

    return

  end subroutine count_nonb_excl_fep

end module sp_enefunc_charmm_mod
