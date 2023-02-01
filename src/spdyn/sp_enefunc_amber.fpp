!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_enefunc_amber_mod
!> @brief   define potential energy functions
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_enefunc_amber_mod

  use sp_enefunc_charmm_mod
  use sp_enefunc_restraints_mod
  use sp_enefunc_table_mod
  use sp_energy_mod
  use sp_restraints_str_mod
  use sp_constraints_str_mod
  use sp_enefunc_str_mod
  use sp_energy_str_mod
  use sp_domain_str_mod
  use dihedral_libs_mod
  use molecules_str_mod
  use fileio_prmtop_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: define_enefunc_amber
  private :: setup_enefunc_bond
  private :: setup_enefunc_bond_constraint
  private :: setup_enefunc_angl
  private :: setup_enefunc_angl_constraint
  private :: setup_enefunc_dihe
  private :: setup_enefunc_impr
  private :: setup_enefunc_cmap
  private :: setup_enefunc_nonb

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc_amber
  !> @brief        a driver subroutine for defining potential energy
  !! @authors      NT, HO
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    prmtop      : AMBER parameter topology information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    constraints : constraints information
  !! @param[in]    restraints  : restraints information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc_amber(ene_info, prmtop, molecule, &
                                  constraints, restraints, domain, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info 
    type(s_prmtop),          intent(in)    :: prmtop
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

    call alloc_enefunc(enefunc, EneFuncBase,     ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncBond,     ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncAngl,     ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncDihe,     ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncRBDihe,   ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncImpr,     ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncBondCell, ncel, ncelb)
    if (domain%fep_use) then
      ! FEP
      call alloc_enefunc(enefunc, EneFuncBondCell_FEP, ncel, ncelb)
      call alloc_enefunc(enefunc, EneFuncFEPBonded, ncel, ncel)
    end if

    if (.not. constraints%rigid_bond) then

      ! bond
      !
      call setup_enefunc_bond(molecule, prmtop, domain, constraints, enefunc)

      ! angle
      !
      call setup_enefunc_angl(molecule, prmtop, domain, enefunc)

    else

      ! bond
      !
      call setup_enefunc_bond_constraint(prmtop, molecule, &
                                         domain, constraints, enefunc)

      ! angle
      !
      call setup_enefunc_angl_constraint(prmtop, molecule, &
                                         domain, constraints, enefunc)

    end if

    ! dihedral
    !
    call setup_enefunc_dihe(molecule, prmtop, domain, enefunc)

    ! improper
    !
    call setup_enefunc_impr(molecule, prmtop, domain, enefunc)

    ! cmap
    !
    call setup_enefunc_cmap(ene_info, molecule, prmtop, domain, enefunc)

    ! nonbonded
    !
    call setup_enefunc_nonb(prmtop, molecule, constraints, domain, enefunc)

    ! lookup table
    !
    if (ene_info%table) &
    call setup_enefunc_table(ene_info, enefunc)

    ! restraints
    !
    call setup_enefunc_restraints(molecule, restraints, domain, enefunc)

    ! write summary of energy function
    !
    if (main_rank) then

      write(MsgOut,'(A)') &
           'Define_Enefunc_Amber> Number of Interactions in Each Term'
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           '  bond_ene        = ', enefunc%num_bond_all,        &
           '  angle_ene       = ', enefunc%num_angl_all
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           '  torsion_ene     = ', enefunc%num_dihe_all,        &
           '  improper_ene    = ', enefunc%num_impr_all
      write(MsgOut,'(A20,I10)')                                 &
           '  cmap_ene        = ', enefunc%num_cmap_all
      if (.not. ene_info%table) then
        write(MsgOut,'(A20,I10,A20,I10)')                       &
           '  nb_exclusions   = ', enefunc%num_excl_all,        &
           '  nb14_calc       = ', enefunc%num_nb14_all
      end if
      if (domain%fep_use) then
        write(MsgOut,'(A20,I10)')                        &
           '  nb14_calc_fep   = ', enefunc%num_nb14_all_fep
      end if
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           ' restraint_groups = ', enefunc%num_restraintgroups, &
           ' restraint_funcs  = ', enefunc%num_restraintfuncs
      write(MsgOut,'(A)') ' '

    end if

    return

  end subroutine define_enefunc_amber

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_bond
  !> @brief        define BOND term in potential energy function
  !! @authors      NT, HO
  !! @param[in]    molecule : molecule information
  !! @param[in]    prmtop   : AMBER parameter topology information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond(molecule, prmtop, domain, constraints, enefunc)

    ! formal arguments
    type(s_molecule),target, intent(in)    :: molecule
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_domain),  target, intent(in)    :: domain
    type(s_constraints),     intent(in)    :: constraints
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: m1, m2
    integer                  :: dupl, ioffset
    integer                  :: i, found, i1, i2, icel1, icel2
    integer                  :: ia, ib, pbc_int
    integer                  :: icel_local, wat_bonds
    character(6)             :: ri1, ri2
    real(wp)                 :: cwork(3,2), dij(3)

    real(wp),        pointer :: force(:,:), dist(:,:), coord(:,:)
    real(wp),        pointer :: box_size(:)
    integer,         pointer :: nwater(:), bond(:), list(:,:,:)
    integer,         pointer :: ncel
    integer(int2),   pointer :: cell_pair(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: nbond
    integer,         pointer :: bond_pbc(:,:)

    ! FEP
    integer                  :: nbond_fep, k, fg1, fg2
    logical                  :: flag_singleB
    integer,         pointer :: bond_singleB(:,:)

    coord     => molecule%atom_coord

    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l
    nwater    => domain%num_water
    box_size  => domain%system_size

    bond      => enefunc%num_bond
    list      => enefunc%bond_list
    force     => enefunc%bond_force_const
    dist      => enefunc%bond_dist_min
    bond_pbc  => enefunc%bond_pbc

    ! FEP
    if (domain%fep_use) then
      bond_singleB => enefunc%bond_singleB
    end if

    do dupl = 1, domain%num_duplicate

      ioffset = (dupl-1) * enefunc%table%num_all

      do i = 1, prmtop%num_bondh

        i1 = prmtop%bond_inc_hy(1,i) / 3 + 1
        i2 = prmtop%bond_inc_hy(2,i) / 3 + 1
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
              ! FEP: If the bond includes singleB atoms
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

          i1 = i1 + ioffset
          i2 = i2 + ioffset

          icel1 = id_g2l(1,i1)
          icel2 = id_g2l(1,i2)

          if (icel1 /= 0 .and. icel2 /= 0) then

            icel_local = cell_pair(icel1,icel2)

            if (icel_local > 0 .and. icel_local <= ncel) then

              nbond => bond(icel_local)
              nbond = nbond + 1

              if (nbond > MaxBond) &
                call error_msg('Setup_Enefunc_Bond> Too many bonds.') 

              list (1,nbond,icel_local) = i1
              list (2,nbond,icel_local) = i2
              force(  nbond,icel_local) = &
                             prmtop%bond_fcons_uniq(prmtop%bond_inc_hy(3,i))
              dist (  nbond,icel_local) = &
                             prmtop%bond_equil_uniq(prmtop%bond_inc_hy(3,i))
              bond_pbc(bond(icel_local),icel_local) = 13

              if (domain%fep_use) then
                ! FEP: flag for singleB bond
                if (flag_singleB) &
                  bond_singleB (nbond,icel_local) = 1
              end if

            end if

          end if

        end if

      end do

      do i = 1, prmtop%num_mbonda

        i1 = prmtop%bond_wo_hy(1,i) / 3 + 1
        i2 = prmtop%bond_wo_hy(2,i) / 3 + 1

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
              ! FEP: If the bond includes singleB atoms
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
 
          i1 = i1 + ioffset
          i2 = i2 + ioffset

          icel1 = id_g2l(1,i1)
          icel2 = id_g2l(1,i2)

          if (icel1 /= 0 .and. icel2 /= 0) then

            icel_local = cell_pair(icel1,icel2)

            if (icel_local > 0 .and. icel_local <= ncel) then

              nbond => bond(icel_local)
              nbond = nbond + 1

              if (nbond > MaxBond) &
                call error_msg('Setup_Enefunc_Bond> Too many bonds.') 

              list (1,nbond,icel_local) = i1
              list (2,nbond,icel_local) = i2
              force(  nbond,icel_local) = &
                               prmtop%bond_fcons_uniq(prmtop%bond_wo_hy(3,i))
              dist (  nbond,icel_local) = &
                               prmtop%bond_equil_uniq(prmtop%bond_wo_hy(3,i))
              bond_pbc(bond(icel_local),icel_local) = 13
  
              if (domain%fep_use) then
                ! FEP: flag for singleB bond
                if (flag_singleB) &
                  bond_singleB (nbond,icel_local) = 1
              end if
 
            end if

          end if

        end if

      end do

    end do

    ! water molecules
    !
    wat_bonds = 0

    do i = 1, prmtop%num_mbonda

      i1 = prmtop%bond_wo_hy(1,i) / 3 + 1
      i2 = prmtop%bond_wo_hy(2,i) / 3 + 1
      ri1 = molecule%residue_name(i1)
      ri2 = molecule%residue_name(i2)
      if (ri1(1:3) .eq. 'TIP' .or. ri1(1:3) .eq. 'WAT' .or. &
          ri1(1:3) .eq. 'SOL') then
        wat_bonds = wat_bonds+1
      end if

    end do

    if (.not.constraints%fast_water .and.  enefunc%table%num_water > 0) then

      do i = 1, prmtop%num_bondh

        i1 = prmtop%bond_inc_hy(1,i) / 3 + 1
        i2 = prmtop%bond_inc_hy(2,i) / 3 + 1
        ri1 = molecule%residue_name(i1)
        ri2 = molecule%residue_name(i2)
        m1 = molecule%mass(i1)
        m2 = molecule%mass(i2)

        if (ri1(1:3) .eq. 'TIP' .or. ri1(1:3) .eq. 'WAT' .or. &
            ri1(1:3) .eq. 'SOL') then

          wat_bonds = wat_bonds + 1
          if (abs(m1-m2) > EPS) then
            enefunc%table%water_bond_calc = .true.
            enefunc%table%water_bond_calc_OH = .true.
            enefunc%table%OH_bond = &
               prmtop%bond_equil_uniq(prmtop%bond_inc_hy(3,i))
            enefunc%table%OH_force = &
               prmtop%bond_fcons_uniq(prmtop%bond_inc_hy(3,i))
          else if (m1 == m2 .and. m1 < LIGHT_ATOM_MASS_LIMIT) then
            enefunc%table%water_bond_calc = .true.
            enefunc%table%water_bond_calc_HH = .true.
            enefunc%table%HH_bond = &
               prmtop%bond_equil_uniq(prmtop%bond_inc_hy(3,i))
            enefunc%table%HH_force = &
               prmtop%bond_fcons_uniq(prmtop%bond_inc_hy(3,i))
          end if
        end if

      end do
    end if

    found = 0 
    do i = 1, ncel
      found = found + bond(i)
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_bond_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_bond_all = found
#endif
    enefunc%num_bond_all = enefunc%num_bond_all + wat_bonds
    enefunc%num_bond_all = enefunc%num_bond_all * domain%num_duplicate

    return

  end subroutine setup_enefunc_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_bond_constraint
  !> @brief        define BOND term in potential energy function
  !! @authors      NT, HO
  !! @param[in]    prmtop      : AMBER parameter topology information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    domain      : domain information
  !! @param[in]    constraints : constraints information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond_constraint(prmtop, molecule, domain, &
                                           constraints, enefunc)

    ! formal arguments
    type(s_prmtop),              intent(in)    :: prmtop
    type(s_molecule),    target, intent(in)    :: molecule
    type(s_domain),      target, intent(in)    :: domain
    type(s_constraints), target, intent(inout) :: constraints
    type(s_enefunc),     target, intent(inout) :: enefunc

    ! local variables
    integer                      :: dupl, ioffset
    integer                      :: i, j, k, ih, icel_local, connect
    integer                      :: i1, i2, ih1, ih2, icel1, icel2, icel
    integer                      :: nbond_a, nbond_c
    integer                      :: wat_bonds, wat_found
    integer                      :: tmp_mole_no, mole_no
    character(6)                 :: ri1, ri2

    real(wp),            pointer :: force(:,:), dist(:,:), coord(:,:)
    real(wip),           pointer :: HGr_bond_dist(:,:,:,:)
    integer,             pointer :: bond(:), list(:,:,:), ncel
    integer,             pointer :: id_l2g(:,:)
    integer(int2),       pointer :: cell_pair(:,:)
    integer(int2),       pointer :: id_g2l(:,:)
    integer,             pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)
    integer,             pointer :: bond_pbc(:,:)

    ! FEP
    integer                      :: l, fg1, fg2
    integer                      :: nbond_c_singleB, nbond_c_singleB_all
    logical                      :: flag_singleB
    integer,             pointer :: bond_singleB(:,:)

    coord         => molecule%atom_coord

    ncel          => domain%num_cell_local
    cell_pair     => domain%cell_pair
    id_g2l        => domain%id_g2l
    id_l2g        => domain%id_l2g_solute

    HGr_local     => constraints%HGr_local
    HGr_bond_list => constraints%HGr_bond_list
    HGr_bond_dist => constraints%HGr_bond_dist

    bond          => enefunc%num_bond
    list          => enefunc%bond_list
    force         => enefunc%bond_force_const
    dist          => enefunc%bond_dist_min
    bond_pbc      => enefunc%bond_pbc

    connect       =  constraints%connect

    nbond_a       =  0
    nbond_c       =  0

    ! FEP
    if (domain%fep_use) then
      bond_singleB => enefunc%bond_singleB
      nbond_c_singleB = 0
    end if

    do dupl = 1, domain%num_duplicate

      ioffset = (dupl-1) * enefunc%table%num_all

      do i = 1, prmtop%num_mbonda

        i1 = prmtop%bond_wo_hy(1,i) / 3 + 1
        i2 = prmtop%bond_wo_hy(2,i) / 3 + 1

        ri1 = molecule%residue_name(i1)
        ri2 = molecule%residue_name(i2)

        ! FEP
        if (domain%fep_use) then
          flag_singleB = .false.
          fg1 = molecule%fepgrp(i1)
          fg2 = molecule%fepgrp(i2)
          if (molecule%fepgrp_bond(fg1,fg2) == 0) then
            ! FEP: If the bond are not set to any group of FEP, exclude this bond.
            cycle
          else if (molecule%fepgrp_bond(fg1,fg2) == 2) then
            ! FEP: If the bond includes singleB atoms
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

        i1  = i1 + ioffset
        i2  = i2 + ioffset

        if (ri1 /= constraints%water_model .and. &
            ri2 /= constraints%water_model) then

          icel1 = id_g2l(1,i1)
          icel2 = id_g2l(1,i2)

          if (icel1 /= 0 .and. icel2 /= 0) then

            icel_local = cell_pair(icel1,icel2)

            if (icel_local > 0 .and. icel_local <= ncel) then

              nbond_a = nbond_a + 1
              bond(icel_local) = bond(icel_local) + 1

              if (bond(icel_local) > MaxBond) &
                call error_msg('Setup_Enefunc_Bond_Constraint> Too many bonds.') 

              list (1,bond(icel_local),icel_local) = i1
              list (2,bond(icel_local),icel_local) = i2
              force(  bond(icel_local),icel_local) = &
                               prmtop%bond_fcons_uniq(prmtop%bond_wo_hy(3,i))
              dist (  bond(icel_local),icel_local) = &
                               prmtop%bond_equil_uniq(prmtop%bond_wo_hy(3,i))
              bond_pbc(bond(icel_local),icel_local) = 13

              if (domain%fep_use) then
                ! FEP: flag for singleB bond
                if (flag_singleB) then
                  bond_singleB (bond(icel_local),icel_local) = 1
                end if
              end if

            end if
          end if

        end if

      end do

      do i = 1, prmtop%num_bondh

        i1 = prmtop%bond_inc_hy(1,i) / 3 + 1
        i2 = prmtop%bond_inc_hy(2,i) / 3 + 1

        ri1 = molecule%residue_name(i1)
        ri2 = molecule%residue_name(i2)

        ! FEP
        if (domain%fep_use) then
          flag_singleB = .false.
          fg1 = molecule%fepgrp(i1)
          fg2 = molecule%fepgrp(i2)
          if (molecule%fepgrp_bond(fg1,fg2) == 0) then
            ! FEP: If the bond are not set to any group of FEP, exclude this bond.
            cycle
          else if (molecule%fepgrp_bond(fg1,fg2) == 2) then
            ! FEP: If the bond includes singleB atoms
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
 
        i1 = i1 + ioffset
        i2 = i2 + ioffset

        if (ri1 /= constraints%water_model .and. &
            ri2 /= constraints%water_model) then

          icel1 = id_g2l(1,i1)
          icel2 = id_g2l(1,i2)

          if (icel1 /= 0 .and. icel2 /= 0) then

            icel = cell_pair(icel1,icel2)

            if (icel > 0 .and. icel <= ncel) then

              do j = 1, connect

                do k = 1, HGr_local(j,icel)
                  ih1 = id_l2g(HGr_bond_list(1,k,j,icel),icel)
                  do ih = 1, j
                    ih2 = id_l2g(HGr_bond_list(ih+1,k,j,icel),icel)
  
                    if (ih1 == i1 .and. ih2 == i2 .or. &
                        ih2 == i1 .and. ih1 == i2) then
  
                      ! FEP
                      if (domain%fep_use) then
                        if (flag_singleB) then
                          ! FEP: skip singleB atoms if the bond includes hydrogens.
                          nbond_c_singleB = nbond_c_singleB + 1
                          cycle
                        end if
                      end if

                      nbond_c = nbond_c + 1
                      HGr_bond_dist(ih+1,k,j,icel) = &
                              prmtop%bond_equil_uniq(prmtop%bond_inc_hy(3,i))
  
                    end if
  
                  end do
                end do
              end do
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

      do i = 1, prmtop%num_mbonda

        i1 = prmtop%bond_wo_hy(1,i) / 3 + 1
        i2 = prmtop%bond_wo_hy(2,i) / 3 + 1

        ri1 = molecule%residue_name(i1)
        ri2 = molecule%residue_name(i2)

        if (ri1 == constraints%water_model .and. &
            ri2 == constraints%water_model) then

          wat_bonds = wat_bonds+1

        end if

      end do

      do i = 1, prmtop%num_bondh
  
        i1 = prmtop%bond_inc_hy(1,i) / 3 + 1
        i2 = prmtop%bond_inc_hy(2,i) / 3 + 1
  
        ri1 = molecule%residue_name(i1)
        ri2 = molecule%residue_name(i2)

        if (ri1 == constraints%water_model .and. &
            ri2 == constraints%water_model) then

          wat_bonds = wat_bonds+1
          mole_no = molecule%molecule_no(molecule%bond_list(1,i))

          if (mole_no /= tmp_mole_no) then
            wat_found = wat_found +1
            tmp_mole_no = mole_no
          end if

        end if

      end do

      if (wat_found /= enefunc%table%num_water) &
        call error_msg( &
          'Setup_Enefunc_Bond_Constraint> # of water is incorrect')

    end if

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(nbond_a, enefunc%num_bond_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
    call mpi_allreduce(nbond_c, constraints%num_bonds, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_bond_all  = nbond_a
    constraints%num_bonds = nbond_c
#endif

    if (domain%fep_use) then

      ! FEP
#ifdef HAVE_MPI_GENESIS
      call mpi_allreduce(nbond_c_singleB, nbond_c_singleB_all, 1, mpi_integer, &
                         mpi_sum, mpi_comm_country, ierror)
#else
      nbond_c_singleB_all = nbond_c_singleB
#endif
      if (constraints%fast_water) then
        if (enefunc%num_bond_all /= prmtop%num_bondh *domain%num_duplicate &
                                    + prmtop%num_mbonda*domain%num_duplicate &
                                    - (constraints%num_bonds + &
                                       nbond_c_singleB_all) &
                                    - wat_bonds*domain%num_duplicate) then
          call error_msg( &
            'Setup_Enefunc_Bond_Constraint> Some bond paremeters are missing.')
        end if
      else
        if (enefunc%num_bond_all /=   prmtop%num_bondh *domain%num_duplicate &
                                      + prmtop%num_mbonda*domain%num_duplicate &
                                      - (constraints%num_bonds + &
                                         nbond_c_singleB_all) &
                                      -3*enefunc%table%num_water*domain%num_duplicate) then
          call error_msg( &
            'Setup_Enefunc_Bond_Constraint> Some bond paremeters are missing.')
        end if
      end if

    else

      if (constraints%fast_water) then

        if (enefunc%num_bond_all /= prmtop%num_bondh *domain%num_duplicate &
                                    + prmtop%num_mbonda*domain%num_duplicate &
                                    - constraints%num_bonds &
                                    - wat_bonds*domain%num_duplicate) then
          call error_msg( &
            'Setup_Enefunc_Bond_Constraint> Some bond paremeters are missing.')
        end if

      else

        if (enefunc%num_bond_all /=   prmtop%num_bondh *domain%num_duplicate &
                                      + prmtop%num_mbonda*domain%num_duplicate &
                                      - constraints%num_bonds &
                                      -3*enefunc%table%num_water*domain%num_duplicate) then
          call error_msg( &
            'Setup_Enefunc_Bond_Constraint> Some bond paremeters are missing.')
        end if

      end if

    end if

    return

  end subroutine setup_enefunc_bond_constraint

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_angl
  !> @brief        define ANGLE term in potential energy function
  !! @authors      NT, HO
  !! @param[in]    molecule : molecule information
  !! @param[in]    prmtop   : AMBER parameter topology information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_angl(molecule, prmtop, domain, enefunc)

    ! formal arguments
    type(s_molecule),target, intent(in)    :: molecule
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: dupl, ioffset
    integer                  :: i, found, i1, i2, i3, icel1, icel2
    integer                  :: icel_local
    character(6)             :: ri1, ri2, ri3

    real(wp),        pointer :: force(:,:), theta(:,:)
    integer,         pointer :: angle(:), alist(:,:,:)
    integer,         pointer :: ncel
    integer(int2),   pointer :: cell_pair(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: nangl, nwater(:)
    character(6),    pointer :: mol_res_name(:)
    integer,         pointer :: angl_pbc(:,:,:)

    ! FEP
    integer                  :: nangl_fep, k, fg(3)
    logical                  :: flag_singleB
    integer,         pointer :: angl_singleB(:,:)

    mol_res_name => molecule%residue_name

    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l
    nwater    => domain%num_water

    angle     => enefunc%num_angle
    alist     => enefunc%angle_list
    force     => enefunc%angle_force_const
    theta     => enefunc%angle_theta_min
    angl_pbc  => enefunc%angle_pbc

    ! FEP
    if (domain%fep_use) then
      angl_singleB => enefunc%angl_singleB
    end if

    do dupl = 1, domain%num_duplicate

      ioffset = (dupl-1)*enefunc%table%num_all

      do i = 1, prmtop%num_anglh

        i1 = prmtop%angl_inc_hy(1,i) / 3 + 1
        i2 = prmtop%angl_inc_hy(2,i) / 3 + 1
        i3 = prmtop%angl_inc_hy(3,i) / 3 + 1

        ri1 = mol_res_name(i1)
        ri2 = mol_res_name(i2)
        ri3 = mol_res_name(i3)

        if (ri1(1:3) .ne. 'TIP' .and. ri1(1:3) .ne. 'WAT' .and. &
            ri1(1:3) .ne. 'SOL' .and. ri3(1:3) .ne. 'TIP' .and. &
            ri3(1:3) .ne. 'WAT' .and. ri3(1:3) .ne. 'SOL') then

          ! FEP
          if (domain%fep_use) then
            flag_singleB = .false.
            fg(1) = molecule%fepgrp(i1)
            fg(2) = molecule%fepgrp(i2)
            fg(3) = molecule%fepgrp(i3)
            if (molecule%fepgrp_angl(fg(1),fg(2),fg(3)) == 0) then
              ! FEP: If the angle are not set to any group of FEP, exclude this angle.
              cycle
            else if (molecule%fepgrp_angl(fg(1),fg(2),fg(3)) == 2) then
              ! FEP: If the angle includes singleB atoms
              flag_singleB = .true.
            end if
            ! FEP: If i1 is singleB atom, its ID is replaced with the ID corresponding to singleA.
            if (fg(1) == 2) then
              do k = 1, molecule%num_atoms_fep(2)
                if (molecule%id_singleB(k) == i1) then
                  exit
                end if
              end do
              i1 = molecule%id_singleA(k)
            end if
            ! FEP: If i2 is singleB atom, its ID is replaced with the ID corresponding to singleA.
            if (fg(2) == 2) then
              do k = 1, molecule%num_atoms_fep(2)
                if (molecule%id_singleB(k) == i2) then
                  exit
                end if
              end do
              i2 = molecule%id_singleA(k)
            end if
            ! FEP: If i3 is singleB atom, its ID is replaced with the ID corresponding to singleA.
            if (fg(3) == 2) then
              do k = 1, molecule%num_atoms_fep(2)
                if (molecule%id_singleB(k) == i3) then
                  exit
                end if
              end do
              i3 = molecule%id_singleA(k)
            end if
          end if

          i1 = i1 + ioffset
          i2 = i2 + ioffset
          i3 = i3 + ioffset

          icel1 = id_g2l(1,i1)
          icel2 = id_g2l(1,i3)

          if (icel1 /= 0 .and. icel2 /= 0) then

            icel_local = cell_pair(icel1,icel2)

            if (icel_local > 0 .and. icel_local <= ncel) then

              nangl => angle(icel_local)
              nangl = nangl + 1

              if (nangl > MaxAngle) &
                call error_msg('Setup_Enefunc_Angl> Too many angles.') 

              alist(1:3,nangl,icel_local) = (/i1, i2, i3/)
              force(    nangl,icel_local) = &
                               prmtop%angl_fcons_uniq(prmtop%angl_inc_hy(4,i))
              theta(    nangl,icel_local) = &
                               prmtop%angl_equil_uniq(prmtop%angl_inc_hy(4,i))
              angl_pbc(1:3,nangl,icel_local) = 13

              if (domain%fep_use) then
                ! FEP: flag for singleB angle
                if (flag_singleB) &
                  angl_singleB (nangl,icel_local) = 1
              end if

            end if

          end if

        end if

      end do

      do i = 1, prmtop%num_mangla

        i1 = prmtop%angl_wo_hy(1,i) / 3 + 1
        i2 = prmtop%angl_wo_hy(2,i) / 3 + 1
        i3 = prmtop%angl_wo_hy(3,i) / 3 + 1

        ! FEP
        if (domain%fep_use) then
          flag_singleB = .false.
          fg(1) = molecule%fepgrp(i1)
          fg(2) = molecule%fepgrp(i2)
          fg(3) = molecule%fepgrp(i3)
          if (molecule%fepgrp_angl(fg(1),fg(2),fg(3)) == 0) then
            ! FEP: If the angle are not set to any group of FEP, exclude this angle.
            cycle
          else if (molecule%fepgrp_angl(fg(1),fg(2),fg(3)) == 2) then
            ! FEP: If the angle includes singleB atoms
            flag_singleB = .true.
          end if
          ! FEP: If i1 is singleB atom, its ID is replaced with the ID corresponding to singleA.
          if (fg(1) == 2) then
            do k = 1, molecule%num_atoms_fep(2)
              if (molecule%id_singleB(k) == i1) then
                exit
              end if
            end do
            i1 = molecule%id_singleA(k)
          end if
          ! FEP: If i2 is singleB atom, its ID is replaced with the ID corresponding to singleA.
          if (fg(2) == 2) then
            do k = 1, molecule%num_atoms_fep(2)
              if (molecule%id_singleB(k) == i2) then
                exit
              end if
            end do
            i2 = molecule%id_singleA(k)
          end if
          ! FEP: If i3 is singleB atom, its ID is replaced with the ID corresponding to singleA.
          if (fg(3) == 2) then
            do k = 1, molecule%num_atoms_fep(2)
              if (molecule%id_singleB(k) == i3) then
                exit
              end if
            end do
            i3 = molecule%id_singleA(k)
          end if
        end if

        i1 = i1 + ioffset
        i2 = i2 + ioffset
        i3 = i3 + ioffset

        icel1 = id_g2l(1,i1)
        icel2 = id_g2l(1,i3)

        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1,icel2)

          if (icel_local > 0 .and. icel_local <= ncel) then

            nangl => angle(icel_local)
            nangl = nangl + 1

            if (nangl > MaxAngle) &
              call error_msg('Setup_Enefunc_Angl> Too many angles.') 

            alist(1:3,nangl,icel_local) = (/i1, i2, i3/)
            force(    nangl,icel_local) = &
                             prmtop%angl_fcons_uniq(prmtop%angl_wo_hy(4,i))
            theta(    nangl,icel_local) = &
                             prmtop%angl_equil_uniq(prmtop%angl_wo_hy(4,i))
            angl_pbc(1:3,nangl,icel_local) = 13

            if (domain%fep_use) then
              ! FEP: flag for singleB angle
              if (flag_singleB) &
                angl_singleB (nangl,icel_local) = 1
            end if

          end if

        end if

      end do

    end do

    ! angle energy for water
    !
    do i = 1, prmtop%num_anglh

      i1 = prmtop%angl_inc_hy(1,i) / 3 + 1
      i2 = prmtop%angl_inc_hy(2,i) / 3 + 1
      i3 = prmtop%angl_inc_hy(3,i) / 3 + 1

      ri1 = mol_res_name(i1)
      ri2 = mol_res_name(i2)
      ri3 = mol_res_name(i3)

      if (ri1(1:3) .eq. 'TIP' .or. ri1(1:3) .eq. 'WAT' .or. &
          ri1(1:3) .eq. 'SOL') then

        enefunc%table%HOH_angle = &
            prmtop%angl_equil_uniq(prmtop%angl_inc_hy(4,i))
        enefunc%table%HOH_force = &
            prmtop%angl_fcons_uniq(prmtop%angl_inc_hy(4,i))
        enefunc%table%water_angle_calc = .true.
        exit
      end if
    end do

    found = 0 
    do i = 1, ncel
      found = found + angle(i)
      if (enefunc%table%water_bond_calc_OH) found = found + nwater(i)
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_angl_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_angl_all = found
#endif

    return

  end subroutine setup_enefunc_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_angl_constraint
  !> @brief        define ANGLE term in potential energy function
  !! @authors      NT, HO
  !! @param[in]    prmtop      : AMBER parameter topology information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    domain      : domain information
  !! @param[in]    constraints : constraints information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_angl_constraint(prmtop, molecule, domain, &
                                           constraints, enefunc)

    ! formal arguments
    type(s_prmtop),              intent(in)    :: prmtop
    type(s_molecule),            intent(in)    :: molecule
    type(s_domain),      target, intent(in)    :: domain
    type(s_constraints), target, intent(in)    :: constraints
    type(s_enefunc),     target, intent(inout) :: enefunc

    ! local variables
    integer                  :: dupl, ioffset
    integer                  :: i, found, i1, i2, i3, icel1, icel2
    integer                  :: icel_local
    character(6)             :: res1, res2, res3

    real(wp),        pointer :: force(:,:), theta(:,:)
    integer,         pointer :: angle(:), alist(:,:,:)
    integer,         pointer :: ncel
    integer(int2),   pointer :: cell_pair(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: nangl
    integer,         pointer :: angl_pbc(:,:,:)

    ! FEP
    integer                  :: k, fg(3), nangle_c_singleB
    integer                  :: nangle_c_singleB_all
    logical                  :: flag_singleB
    integer,         pointer :: angl_singleB(:,:)

    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l

    angle     => enefunc%num_angle
    alist     => enefunc%angle_list
    force     => enefunc%angle_force_const
    theta     => enefunc%angle_theta_min
    angl_pbc  => enefunc%angle_pbc

    ! FEP
    if (domain%fep_use) then
      angl_singleB => enefunc%angl_singleB
      nangle_c_singleB = 0
    end if

    do dupl = 1, domain%num_duplicate

      ioffset = (dupl-1) * enefunc%table%num_all

      do i = 1, prmtop%num_anglh

        i1 = prmtop%angl_inc_hy(1,i) / 3 + 1
        i2 = prmtop%angl_inc_hy(2,i) / 3 + 1
        i3 = prmtop%angl_inc_hy(3,i) / 3 + 1

        res1 = molecule%residue_name(i1)
        res2 = molecule%residue_name(i2)
        res3 = molecule%residue_name(i3)

        ! FEP
        if (domain%fep_use) then
          flag_singleB = .false.
          fg(1) = molecule%fepgrp(i1)
          fg(2) = molecule%fepgrp(i2)
          fg(3) = molecule%fepgrp(i3)
          if (molecule%fepgrp_angl(fg(1),fg(2),fg(3)) == 0) then
            ! FEP: If the angle are not set to any group of FEP, exclude this angle.
            cycle
          else if (molecule%fepgrp_angl(fg(1),fg(2),fg(3)) == 2) then
            ! FEP: If the angle includes singleB atoms
            flag_singleB = .true.
          end if
          ! FEP: If i1 is singleB atom, its ID is replaced with the ID corresponding to singleA.
          if (fg(1) == 2) then
            do k = 1, molecule%num_atoms_fep(2)
              if (molecule%id_singleB(k) == i1) then
                exit
              end if
            end do
            i1 = molecule%id_singleA(k)
          end if
          ! FEP: If i2 is singleB atom, its ID is replaced with the ID corresponding to singleA.
          if (fg(2) == 2) then
            do k = 1, molecule%num_atoms_fep(2)
              if (molecule%id_singleB(k) == i2) then
                exit
              end if
            end do
            i2 = molecule%id_singleA(k)
          end if
          ! FEP: If i3 is singleB atom, its ID is replaced with the ID corresponding to singleA.
          if (fg(3) == 2) then
            do k = 1, molecule%num_atoms_fep(2)
              if (molecule%id_singleB(k) == i3) then
                exit
              end if
            end do
            i3 = molecule%id_singleA(k)
          end if
        end if

        i1 = i1 + ioffset
        i2 = i2 + ioffset
        i3 = i3 + ioffset

        if (res1 /= constraints%water_model .and. &
            res2 /= constraints%water_model .and. &
            res3 /= constraints%water_model) then

          icel1 = id_g2l(1,i1)
          icel2 = id_g2l(1,i3)
  
          if (icel1 /= 0 .and. icel2 /= 0) then
  
            icel_local = cell_pair(icel1,icel2)
  
            if (icel_local > 0 .and. icel_local <= ncel) then
  
              nangl => angle(icel_local)
              nangl = nangl + 1
  
              if (nangl > MaxAngle) &
                call error_msg('Setup_Enefunc_Angl_Constraint> Too many angles.') 
  
              alist(1:3,nangl,icel_local) = (/i1, i2, i3/)
              force(    nangl,icel_local) = &
                               prmtop%angl_fcons_uniq(prmtop%angl_inc_hy(4,i))
              theta(    nangl,icel_local) = &
                               prmtop%angl_equil_uniq(prmtop%angl_inc_hy(4,i))
              angl_pbc(1:3,nangl,icel_local) = 13

              if (domain%fep_use) then
                ! FEP: flag for singleB angle
                if (flag_singleB) &
                  angl_singleB (nangl,icel_local) = 1
              end if

            end if
  
          end if

        end if

      end do

      do i = 1, prmtop%num_mangla

        i1 = prmtop%angl_wo_hy(1,i) / 3 + 1
        i2 = prmtop%angl_wo_hy(2,i) / 3 + 1
        i3 = prmtop%angl_wo_hy(3,i) / 3 + 1

        ! FEP
        if (domain%fep_use) then
          flag_singleB = .false.
          fg(1) = molecule%fepgrp(i1)
          fg(2) = molecule%fepgrp(i2)
          fg(3) = molecule%fepgrp(i3)
          if (molecule%fepgrp_angl(fg(1),fg(2),fg(3)) == 0) then
            ! FEP: If the angle are not set to any group of FEP, exclude this angle.
            cycle
          else if (molecule%fepgrp_angl(fg(1),fg(2),fg(3)) == 2) then
            ! FEP: If the angle includes singleB atoms
            flag_singleB = .true.
          end if
          ! FEP: If i1 is singleB atom, its ID is replaced with the ID corresponding to singleA.
          if (fg(1) == 2) then
            do k = 1, molecule%num_atoms_fep(2)
              if (molecule%id_singleB(k) == i1) then
                exit
              end if
            end do
            i1 = molecule%id_singleA(k)
          end if
          ! FEP: If i2 is singleB atom, its ID is replaced with the ID corresponding to singleA.
          if (fg(2) == 2) then
            do k = 1, molecule%num_atoms_fep(2)
              if (molecule%id_singleB(k) == i2) then
                exit
              end if
            end do
            i2 = molecule%id_singleA(k)
          end if
          ! FEP: If i3 is singleB atom, its ID is replaced with the ID corresponding to singleA.
          if (fg(3) == 2) then
            do k = 1, molecule%num_atoms_fep(2)
              if (molecule%id_singleB(k) == i3) then
                exit
              end if
            end do
            i3 = molecule%id_singleA(k)
          end if
        end if

        i1 = i1 + ioffset
        i2 = i2 + ioffset
        i3 = i3 + ioffset

        icel1 = id_g2l(1,i1)
        icel2 = id_g2l(1,i3)

        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1,icel2)

          if (icel_local > 0 .and. icel_local <= ncel) then

            nangl => angle(icel_local)
            nangl = nangl + 1

            if (nangl > MaxAngle) &
              call error_msg('Setup_Enefunc_Angl_Constraint> Too many angles.') 

            alist(1:3,nangl,icel_local) = (/i1, i2, i3/)
            force(    nangl,icel_local) = &
                             prmtop%angl_fcons_uniq(prmtop%angl_wo_hy(4,i))
            theta(    nangl,icel_local) = &
                             prmtop%angl_equil_uniq(prmtop%angl_wo_hy(4,i))
            angl_pbc(1:3,nangl,icel_local) = 13

            if (domain%fep_use) then
              ! FEP: flag for singleB angle
              if (flag_singleB) &
                angl_singleB (nangl,icel_local) = 1
            end if

          end if

        end if

      end do

    end do

    found = 0 
    do i = 1, ncel
      found = found + angle(i)
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_angl_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_angl_all = found
#endif

    if (enefunc%num_angl_all /= (prmtop%num_anglh + prmtop%num_mangla)       &
                                 *domain%num_duplicate                 .and. &
        enefunc%num_angl_all /= (prmtop%num_anglh + prmtop%num_mangla - enefunc%table%num_water)  &
                                 *domain%num_duplicate ) &
      call error_msg( &
        'Setup_Enefunc_Angl_Constraint> Some angle paremeters are missing.')

    return

  end subroutine setup_enefunc_angl_constraint

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_dihe
  !> @brief        define DIHEDRAL term in potential energy function
  !! @authors      NT, HO
  !! @param[in]    molecule : molecule information
  !! @param[in]    prmtop   : AMBER parameter topology information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_dihe(molecule, prmtop, domain, enefunc)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: dupl, ioffset
    integer                  :: i, found
    integer                  :: i1, i2, i3, i4, icel1, icel2, icel_local

    real(wp),           pointer :: force(:,:), phase(:,:)
    integer,            pointer :: dihe(:), list(:,:,:), period(:,:)
    integer,            pointer :: ncel
    integer(int2),      pointer :: cell_pair(:,:)
    integer(int2),      pointer :: id_g2l(:,:)
    integer,            pointer :: ndihe
    integer,            pointer :: notation
    integer,            pointer :: dihe_pbc(:,:,:)

    ! FEP
    integer                     :: j, k, fg(4), atom_list(4)
    logical                     :: flag_singleB
    integer,            pointer :: dihe_singleB(:,:)

    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l

    dihe      => enefunc%num_dihedral
    list      => enefunc%dihe_list
    force     => enefunc%dihe_force_const
    phase     => enefunc%dihe_phase
    period    => enefunc%dihe_periodicity
    notation  => enefunc%notation_14types
    dihe_pbc  => enefunc%dihe_pbc

    ! FEP
    if (domain%fep_use) then
      dihe_singleB => enefunc%dihe_singleB
    end if

    notation = 100
    if (prmtop%num_uniqdihe > 100) then
      notation = 1000
      if (prmtop%num_uniqdihe > 1000) then
        call error_msg('Setup_Enefunc_Dihe> Too many dihedrals.') 
      end if
    end if

    do dupl = 1, domain%num_duplicate

      ioffset = (dupl-1)*molecule%num_atoms

      do i = 1, prmtop%num_diheh
      
        if (prmtop%dihe_inc_hy(4,i) < 0) &
          cycle

        i1 =      prmtop%dihe_inc_hy(1,i)  / 3 + 1
        i2 =      prmtop%dihe_inc_hy(2,i)  / 3 + 1
        i3 = iabs(prmtop%dihe_inc_hy(3,i)) / 3 + 1
        i4 =      prmtop%dihe_inc_hy(4,i)  / 3 + 1

        ! FEP
        if (domain%fep_use) then
          atom_list = (/i1, i2, i3, i4/)
          flag_singleB = .false.
          do j = 1, 4
            fg(j) = molecule%fepgrp(atom_list(j))
          end do
          if (molecule%fepgrp_dihe(fg(1),fg(2),fg(3),fg(4)) == 0) then
            ! FEP: If the dihedral is not set to any group of FEP, exclude this angle.
            cycle
          else if (molecule%fepgrp_dihe(fg(1),fg(2),fg(3),fg(4)) == 2) then
            ! FEP: If the dihedral includes singleB atoms
            flag_singleB = .true.
          end if
          ! FEP: If atom_list(j) is singleB atom, its ID is replaced with the ID corresponding to singleA.
          do j = 1, 4
            if (fg(j) == 2) then
              do k = 1, molecule%num_atoms_fep(2)
                if (molecule%id_singleB(k) == atom_list(j)) then
                  exit
                end if
              end do
              atom_list(j) = molecule%id_singleA(k)
            end if
          end do
          i1 = atom_list(1)
          i2 = atom_list(2)
          i3 = atom_list(3)
          i4 = atom_list(4)
        end if

        i1 = i1 + ioffset
        i2 = i2 + ioffset
        i3 = i3 + ioffset
        i4 = i4 + ioffset

        icel1 = id_g2l(1,i1)
        icel2 = id_g2l(1,i4)

        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1,icel2)

          if (icel_local > 0 .and. icel_local <= ncel) then

            ndihe => dihe(icel_local)
            ndihe = ndihe + 1

            if (ndihe > MaxDihe) &
              call error_msg('Setup_Enefunc_Dihe> Too many dihedrals.') 

            list(1:4,ndihe,icel_local) = (/i1,i2,i3,i4/)
            force (  ndihe,icel_local) = &
                              prmtop%dihe_fcons_uniq(prmtop%dihe_inc_hy(5,i))
            phase (  ndihe,icel_local) = &
                              prmtop%dihe_phase_uniq(prmtop%dihe_inc_hy(5,i))
            period(  ndihe,icel_local) = &
                          int(prmtop%dihe_perio_uniq(prmtop%dihe_inc_hy(5,i)))
            if (prmtop%lscee_scale_factor .or. prmtop%lscnb_scale_factor ) then
              period(ndihe,icel_local) = period(  ndihe,icel_local) &
                              + prmtop%dihe_inc_hy(5,i)*notation
            end if
            dihe_pbc(1:3,ndihe,icel_local) = 13

            if (domain%fep_use) then
              ! FEP: flag for singleB dihedral
              if (flag_singleB) &
                dihe_singleB (ndihe,icel_local) = 1
            end if

          end if

        end if

      end do

      do i = 1, prmtop%num_mdihea

        if (prmtop%dihe_wo_hy(4,i) < 0) &
          cycle

        i1 =      prmtop%dihe_wo_hy(1,i)  / 3 + 1
        i2 =      prmtop%dihe_wo_hy(2,i)  / 3 + 1
        i3 = iabs(prmtop%dihe_wo_hy(3,i)) / 3 + 1
        i4 =      prmtop%dihe_wo_hy(4,i)  / 3 + 1

        ! FEP
        if (domain%fep_use) then
          atom_list = (/i1, i2, i3, i4/)
          flag_singleB = .false.
          do j = 1, 4
            fg(j) = molecule%fepgrp(atom_list(j))
          end do
          if (molecule%fepgrp_dihe(fg(1),fg(2),fg(3),fg(4)) == 0) then
            ! FEP: If the dihedral is not set to any group of FEP, exclude this angle.
            cycle
          else if (molecule%fepgrp_dihe(fg(1),fg(2),fg(3),fg(4)) == 2) then
            ! FEP: If the dihedral includes singleB atoms
            flag_singleB = .true.
          end if
          ! FEP: If atom_list(j) is singleB atom, its ID is replaced with the ID corresponding to singleA.
          do j = 1, 4
            if (fg(j) == 2) then
              do k = 1, molecule%num_atoms_fep(2)
                if (molecule%id_singleB(k) == atom_list(j)) then
                  exit
                end if
              end do
              atom_list(j) = molecule%id_singleA(k)
            end if
          end do
          i1 = atom_list(1)
          i2 = atom_list(2)
          i3 = atom_list(3)
          i4 = atom_list(4)
        end if

        i1 = i1 + ioffset
        i2 = i2 + ioffset
        i3 = i3 + ioffset
        i4 = i4 + ioffset

        icel1 = id_g2l(1,i1)
        icel2 = id_g2l(1,i4)

        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1,icel2)

          if (icel_local > 0 .and. icel_local <= ncel) then

            ndihe => dihe(icel_local)
            ndihe = ndihe + 1

            if (ndihe > MaxDihe) &
              call error_msg('Setup_Enefunc_Dihe> Too many dihedrals.') 

            list(1:4,ndihe,icel_local) = (/i1,i2,i3,i4/)
            force (  ndihe,icel_local) = &
                              prmtop%dihe_fcons_uniq(prmtop%dihe_wo_hy(5,i))
            phase (  ndihe,icel_local) = &
                              prmtop%dihe_phase_uniq(prmtop%dihe_wo_hy(5,i))
            period(  ndihe,icel_local) = &
                          int(prmtop%dihe_perio_uniq(prmtop%dihe_wo_hy(5,i)))
            if (prmtop%lscee_scale_factor .or. prmtop%lscnb_scale_factor ) then
              period(ndihe,icel_local) = period(  ndihe,icel_local) &
                              + prmtop%dihe_wo_hy(5,i)*notation
            end if
            dihe_pbc(1:3,ndihe,icel_local) = 13

            if (domain%fep_use) then
              ! FEP: flag for singleB dihedral
              if (flag_singleB) &
                dihe_singleB (ndihe,icel_local) = 1
            end if

          end if

        end if

      end do

    end do

    found = 0 
    do i = 1, ncel
      found = found + dihe(i)
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_dihe_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_dihe_all = found
#endif

    return

  end subroutine setup_enefunc_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_impr
  !> @brief        define IMPROPER dihedral term in potential energy function
  !! @authors      NT, HO
  !! @param[in]    molecule : molecule information
  !! @param[in]    prmtop   : AMBER parameter topology information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_impr(molecule, prmtop, domain, enefunc)

    ! formal arguments
    type(s_molecule),target, intent(in)    :: molecule
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: dupl, ioffset, pbc_int
    integer                  :: i, found
    integer                  :: i1, i2, i3, i4, icel1, icel2, icel_local
    integer                  :: ia, ib, ic, id
    real(wp)                 :: cwork(3,4), dij(3)

    real(wp),           pointer :: force(:,:), phase(:,:), coord(:,:)
    real(wp),           pointer :: box_size(:)
    integer,            pointer :: impr(:), list(:,:,:), period(:,:)
    integer,            pointer :: ncel
    integer(int2),      pointer :: cell_pair(:,:)
    integer(int2),      pointer :: id_g2l(:,:)
    integer,            pointer :: nimpr
    integer,            pointer :: notation
    integer,            pointer :: impr_pbc(:,:,:)

    ! FEP
    integer                     :: j, k, fg(4), atom_list(4)
    logical                     :: flag_singleB
    integer,            pointer :: impr_singleB(:,:)

    coord     => molecule%atom_coord

    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l
    box_size  => domain%system_size

    impr      => enefunc%num_improper
    list      => enefunc%impr_list
    force     => enefunc%impr_force_const
    phase     => enefunc%impr_phase
    period    => enefunc%impr_periodicity
    notation  => enefunc%notation_14types
    impr_pbc  => enefunc%impr_pbc

    ! FEP
    if (domain%fep_use) then
      impr_singleB => enefunc%impr_singleB
    end if

    do dupl = 1, domain%num_duplicate

      ioffset = molecule%num_atoms*(dupl-1)

      do i = 1, prmtop%num_diheh
      
        if (prmtop%dihe_inc_hy(4,i) >= 0) &
          cycle

        i1 =      prmtop%dihe_inc_hy(1,i)  / 3 + 1
        i2 =      prmtop%dihe_inc_hy(2,i)  / 3 + 1
        i3 = iabs(prmtop%dihe_inc_hy(3,i)) / 3 + 1
        i4 = iabs(prmtop%dihe_inc_hy(4,i)) / 3 + 1

        ! FEP
        if (domain%fep_use) then
          atom_list = (/i1, i2, i3, i4/)
          flag_singleB = .false.
          do j = 1, 4
            fg(j) = molecule%fepgrp(atom_list(j))
          end do
          if (molecule%fepgrp_dihe(fg(1),fg(2),fg(3),fg(4)) == 0) then
            ! FEP: If the dihedral is not set to any group of FEP, exclude this angle.
            cycle
          else if (molecule%fepgrp_dihe(fg(1),fg(2),fg(3),fg(4)) == 2) then
            ! FEP: If the dihedral includes singleB atoms
            flag_singleB = .true.
          end if
          ! FEP: If atom_list(j) is singleB atom, its ID is replaced with the ID corresponding to singleA.
          do j = 1, 4
            if (fg(j) == 2) then
              do k = 1, molecule%num_atoms_fep(2)
                if (molecule%id_singleB(k) == atom_list(j)) then
                  exit
                end if
              end do
              atom_list(j) = molecule%id_singleA(k)
            end if
          end do
          i1 = atom_list(1)
          i2 = atom_list(2)
          i3 = atom_list(3)
          i4 = atom_list(4)
        end if

        ia = i1 + ioffset
        ib = i2 + ioffset
        ic = i3 + ioffset
        id = i4 + ioffset

        icel1 = id_g2l(1,ia)
        icel2 = id_g2l(1,id)

        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1,icel2)

          if (icel_local > 0 .and. icel_local <= ncel) then

            nimpr => impr(icel_local)
            nimpr = nimpr + 1

            if (nimpr > MaxImpr) &
              call error_msg('Setup_Enefunc_Impr> Too many impropers.') 

            list(1:4,nimpr,icel_local) = (/ia,ib,ic,id/)
            force (  nimpr,icel_local) = &
                              prmtop%dihe_fcons_uniq(prmtop%dihe_inc_hy(5,i))
            phase (  nimpr,icel_local) = &
                              prmtop%dihe_phase_uniq(prmtop%dihe_inc_hy(5,i))
            period(  nimpr,icel_local) = &
                         int(prmtop%dihe_perio_uniq(prmtop%dihe_inc_hy(5,i)))
            impr_pbc(1,nimpr,icel_local) = 13
            impr_pbc(2,nimpr,icel_local) = 13
            impr_pbc(3,nimpr,icel_local) = 13
            if (prmtop%lscee_scale_factor .or. prmtop%lscnb_scale_factor ) then
              period(nimpr,icel_local) = period(  nimpr,icel_local) &
                              + prmtop%dihe_inc_hy(5,i)*notation
            end if

            if (domain%fep_use) then
              ! FEP: flag for singleB improper
              if (flag_singleB) &
                impr_singleB (nimpr,icel_local) = 1
            end if

          end if

        end if

      end do

      do i = 1, prmtop%num_mdihea

        if (prmtop%dihe_wo_hy(4,i) >= 0) &
          cycle

        i1 =      prmtop%dihe_wo_hy(1,i)  / 3 + 1
        i2 =      prmtop%dihe_wo_hy(2,i)  / 3 + 1
        i3 = iabs(prmtop%dihe_wo_hy(3,i)) / 3 + 1
        i4 = iabs(prmtop%dihe_wo_hy(4,i)) / 3 + 1

        ! FEP
        if (domain%fep_use) then
          atom_list = (/i1, i2, i3, i4/)
          flag_singleB = .false.
          do j = 1, 4
            fg(j) = molecule%fepgrp(atom_list(j))
          end do
          if (molecule%fepgrp_dihe(fg(1),fg(2),fg(3),fg(4)) == 0) then
            ! FEP: If the dihedral is not set to any group of FEP, exclude this angle.
            cycle
          else if (molecule%fepgrp_dihe(fg(1),fg(2),fg(3),fg(4)) == 2) then
            ! FEP: If the dihedral includes singleB atoms
            flag_singleB = .true.
          end if
          ! FEP: If atom_list(j) is singleB atom, its ID is replaced with the ID corresponding to singleA.
          do j = 1, 4
            if (fg(j) == 2) then
              do k = 1, molecule%num_atoms_fep(2)
                if (molecule%id_singleB(k) == atom_list(j)) then
                  exit
                end if
              end do
              atom_list(j) = molecule%id_singleA(k)
            end if
          end do
          i1 = atom_list(1)
          i2 = atom_list(2)
          i3 = atom_list(3)
          i4 = atom_list(4)
        end if

        ia = i1 + ioffset
        ib = i2 + ioffset
        ic = i3 + ioffset
        id = i4 + ioffset

        icel1 = id_g2l(1,ia)
        icel2 = id_g2l(1,id)

        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1,icel2)

          if (icel_local > 0 .and. icel_local <= ncel) then

            nimpr => impr(icel_local)
            nimpr = nimpr + 1

            if (nimpr > MaxImpr) &
              call error_msg('Setup_Enefunc_Impr> Too many impropers.') 

            list(1:4,nimpr,icel_local) = (/ia,ib,ic,id/)
            force (  nimpr,icel_local) = &
                              prmtop%dihe_fcons_uniq(prmtop%dihe_wo_hy(5,i))
            phase (  nimpr,icel_local) = &
                              prmtop%dihe_phase_uniq(prmtop%dihe_wo_hy(5,i))
            period(  nimpr,icel_local) = &
                            int(prmtop%dihe_perio_uniq(prmtop%dihe_wo_hy(5,i)))
            impr_pbc(1,nimpr,icel_local) = 13
            impr_pbc(2,nimpr,icel_local) = 13
            impr_pbc(3,nimpr,icel_local) = 13

            if (prmtop%lscee_scale_factor .or. prmtop%lscnb_scale_factor ) then
              period(nimpr,icel_local) = period(  nimpr,icel_local) &
                              + prmtop%dihe_wo_hy(5,i)*notation
            end if

            if (domain%fep_use) then
              ! FEP: flag for singleB improper
              if (flag_singleB) &
                impr_singleB (nimpr,icel_local) = 1
            end if

          end if

        end if

      end do

    end do

    found = 0 
    do i = 1, ncel
      found = found + impr(i)
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_impr_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_impr_all = found
#endif

    return

  end subroutine setup_enefunc_impr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cmap
  !> @brief        define CMAP term in potential energy function
  !! @authors      NT, HO
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[in]    molecule : molecule information
  !! @param[in]    prmtop   : AMBER parameter topology information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cmap(ene_info, molecule, prmtop, domain, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_molecule),        intent(in)    :: molecule
    type(s_prmtop),          intent(in)    :: prmtop
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

    integer,            pointer :: ncel
    integer(int2),      pointer :: cell_pair(:,:)
    integer(int2),      pointer :: id_g2l(:,:)
    integer,            pointer :: notation
    integer,            pointer :: dihe_pbc(:,:,:)

    real(wp),       allocatable :: c_ij(:,:,:,:)

    ! FEP
    integer                     :: idx, fg(8)
    logical                     :: flag_singleB

    ncel           => domain%num_cell_local
    cell_pair      => domain%cell_pair
    id_g2l         => domain%id_g2l

    dihe_pbc       => enefunc%dihe_pbc

    periodic = ene_info%cmap_pspline
    ncmap_p  = prmtop%num_cmaptype

    ngrid0 = 0
    do i = 1, ncmap_p
      ngrid0 = max(ngrid0, prmtop%cmap_resolution(i))
    end do

    call alloc_enefunc(enefunc, EneFuncCmap, ncel, ngrid0, ncmap_p)

    alloc_stat = 0
    allocate(c_ij(4,4,ngrid0,ngrid0), stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    do i = 1, ncmap_p
      enefunc%cmap_resolution(i) = prmtop%cmap_resolution(i)
    end do

    ! derive cmap coefficients by bicubic interpolation
    !
    do ityp = 1, ncmap_p
      if (periodic) then
        call derive_cmap_amber_coefficients_p(ityp, prmtop, c_ij)
      else
        call derive_cmap_amber_coefficients_np(ityp, prmtop, c_ij)
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

      ioffset = (dupl-1)*molecule%num_atoms

      do i = 1, prmtop%num_cmap

        list(1:4) = prmtop%cmap_list(1:4,i)        
        list(5:8) = prmtop%cmap_list(2:5,i)        

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
            ! FEP: If the cmap includes singleB atoms
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

        lists(1:8) = list(1:8) + ioffset

        icel1 = id_g2l(1,lists(1))
        icel2 = id_g2l(1,lists(8))

        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1,icel2)

          if (icel_local > 0 .and. icel_local <= ncel) then

            enefunc%num_cmap(icel_local) = enefunc%num_cmap(icel_local) + 1
            k = enefunc%num_cmap(icel_local)
            enefunc%cmap_list(1:8,k,icel_local) = lists(1:8)
            enefunc%cmap_type(k,icel_local) = prmtop%cmap_list(6,i)
            enefunc%cmap_pbc(1:6,k,icel_local) = 13

            if (domain%fep_use) then
              ! FEP: flag for singleB cmap
              if (flag_singleB) &
                enefunc%cmap_singleB(k,icel_local) = 1
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

    if (enefunc%num_cmap_all /=  prmtop%num_cmap*domain%num_duplicate) &
      call error_msg('Setup_Enefunc_Cmap> Some cmap parameters are missing.')

    return

  end subroutine setup_enefunc_cmap

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_nonb
  !> @brief        define NON-BOND term in potential energy function
  !! @authors      NT
  !! @param[in]    prmtop      : AMBER parameter topology information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    constraints : constraints information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_nonb(prmtop, molecule, constraints, domain, enefunc)

    ! formal arguments
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_molecule),        intent(in)    :: molecule
    type(s_constraints),     intent(in)    :: constraints
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, nnonb, ncel, ndihetype, ndihe
    integer                  :: nb_par_idx
    integer                  :: k, ix, jx, cls_local

    integer,  allocatable    :: check_cls(:)
    integer,  allocatable    :: atmcls_map_g2l(:), atmcls_map_l2g(:)
    real(wp), allocatable    :: nb14_lj6(:,:), nb14_lj12(:,:)
    real(wp), allocatable    :: nonb_lj6(:,:), nonb_lj12(:,:)


    ndihetype = size(prmtop%scee_scale_fact)
    call alloc_enefunc(enefunc, EneFuncAMBERScale, ndihetype, 0)

    ELECOEF          = ELECOEF_AMBER

    ndihe = 0

    do i = 1, prmtop%num_diheh
      if (prmtop%dihe_inc_hy(4,i) < 0) &
        cycle

      ndihe = ndihe + 1

      if (prmtop%lscee_scale_factor) then
        enefunc%dihe_scee(prmtop%dihe_inc_hy(5,i)) = &
         1.0_wp/prmtop%scee_scale_fact(prmtop%dihe_inc_hy(5,i))
      else
        enefunc%dihe_scee(0) = 1.0_wp/1.2_wp
      end if

      if (prmtop%lscnb_scale_factor) then
        enefunc%dihe_scnb(prmtop%dihe_inc_hy(5,i)) = &
         1.0_wp/prmtop%scnb_scale_fact(prmtop%dihe_inc_hy(5,i))
      else
        enefunc%dihe_scnb(0) = 1.0_wp/2.0_wp
      end if

    end do

    do i = 1, prmtop%num_mdihea
      if (prmtop%dihe_wo_hy(4,i) < 0) &
        cycle

      ndihe = ndihe + 1

      if (prmtop%lscee_scale_factor) then
        enefunc%dihe_scee(prmtop%dihe_wo_hy(5,i)) =  &
          1.0_wp/prmtop%scee_scale_fact(prmtop%dihe_wo_hy(5,i))
      else
        enefunc%dihe_scee(0) = 1.0_wp/1.2_wp
      end if

      if (prmtop%lscnb_scale_factor) then
        enefunc%dihe_scnb(prmtop%dihe_wo_hy(5,i)) = &
          1.0_wp/prmtop%scnb_scale_fact(prmtop%dihe_wo_hy(5,i))
      else
        enefunc%dihe_scnb(0) = 1.0_wp/2.0_wp
      end if

    end do

    nnonb = prmtop%num_types

    allocate(check_cls(nnonb),             &
             atmcls_map_g2l(nnonb),        &
             atmcls_map_l2g(nnonb),        &
             nb14_lj6      (nnonb, nnonb), &
             nb14_lj12     (nnonb, nnonb), &
             nonb_lj6      (nnonb, nnonb), &
             nonb_lj12     (nnonb, nnonb))

    check_cls(1:nnonb)          = 0
    nb14_lj6 (1:nnonb, 1:nnonb) = 0.0_wp
    nb14_lj12(1:nnonb, 1:nnonb) = 0.0_wp
    nonb_lj6 (1:nnonb, 1:nnonb) = 0.0_wp
    nonb_lj12(1:nnonb, 1:nnonb) = 0.0_wp

    do i = 1, nnonb
      do j = 1, nnonb

        nb_par_idx = prmtop%nb_par_idx(nnonb*(i-1)+j)

        if (nb_par_idx < 0) then
          nb14_lj12(i,j) = 0.0_wp
          nb14_lj6 (i,j) = 0.0_wp
          nonb_lj12(i,j) = 0.0_wp
          nonb_lj6 (i,j) = 0.0_wp
        else
          nb14_lj12(i,j) = prmtop%lennarda(nb_par_idx)
          nb14_lj6 (i,j) = prmtop%lennardb(nb_par_idx)
          nonb_lj12(i,j) = prmtop%lennarda(nb_par_idx)
          nonb_lj6 (i,j) = prmtop%lennardb(nb_par_idx)
        end if

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
      check_cls(k) = 1
    end do

    k = 0
    do i = 1, nnonb
      if (check_cls(i) == 1) then
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
      do j = 1, cls_local
        jx = atmcls_map_l2g(j)
        enefunc%nb14_lj12(i,j) = nb14_lj12(ix,jx)
        enefunc%nb14_lj6 (i,j) = nb14_lj6 (ix,jx)
        enefunc%nonb_lj12(i,j) = nonb_lj12(ix,jx)
        enefunc%nonb_lj6 (i,j) = nonb_lj6 (ix,jx)
      end do
    end do

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
      enefunc%table%atom_cls_no_O = atmcls_map_g2l(enefunc%table%atom_cls_no_O)
      enefunc%table%atom_cls_no_H = atmcls_map_g2l(enefunc%table%atom_cls_no_H)
      if (constraints%water_type == TIP4) then
        enefunc%table%atom_cls_no_D = atmcls_map_g2l(enefunc%table%atom_cls_no_D)
        domain%water%atom_cls_no(4) = atmcls_map_g2l(domain%water%atom_cls_no(4))
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

    deallocate(check_cls,     &
               atmcls_map_g2l,&
               atmcls_map_l2g,&
               nb14_lj6,      &
               nb14_lj12,     &
               nonb_lj6,      &
               nonb_lj12)

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

end module sp_enefunc_amber_mod
