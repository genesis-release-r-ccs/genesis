!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_fep_topology_mod
!> @brief   define topology information for FEP
!! @authors Nobuhiko Kato (NK), Hiraku Oshima (HO)
!
!  (c) Copyright 2019 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_fep_topology_mod

  use fileio_par_mod
  use fileio_prmtop_mod
  use select_mod
  use select_atoms_mod
  use select_atoms_str_mod
  use molecules_str_mod
  use fileio_control_mod
  use messages_mod
  use mpi_parallel_mod
  use sp_alchemy_mod
  use sp_alchemy_str_mod

  implicit none
  private

  public  :: define_fep_topology
  private :: delete_bond
  private :: delete_angle
  private :: delete_dihedral
  private :: delete_improper
  private :: delete_cmap
  private :: make_bondlist
  private :: make_anglelist
  private :: make_dihedrallist
  private :: make_improperlist
  private :: make_cmaplist
  private :: check_parameter
  private :: count_hbonds_single_fep
  private :: correspondence_single_fep
  private :: make_bond_group_fep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_fep_topology
  !> @brief        a driver subroutine for defining topology in FEP
  !! @authors      NK, HO
  !! @param[inout] molecule : structure of molecule
  !! @param[in]    par      : PAR information (optional)
  !! @param[in]    prmtop   : PRMTOP information (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_fep_topology(molecule, par, prmtop, sel_info, alch_info)

    ! formal arguments
    type(s_molecule),         intent(inout) :: molecule
    type(s_par),              intent(in)    :: par
    type(s_prmtop),           intent(in)    :: prmtop
    type(s_sel_info),         intent(in)    :: sel_info
    type(s_alch_info),        intent(in)    :: alch_info


    ! local variables
    integer                  :: nfunc, ngroup, natom
    integer                  :: alloc_stat
    integer                  :: dealloc_stat
    integer                  :: i, j, k
    integer                  :: nerr, nerr_cls, nerr_chg
    integer                  :: idx_select

    type(s_selatoms),        allocatable :: selatoms(:)

    ! Setup molecule
    !
    molecule%fep_topology = alch_info%fep_topology

    if (FEPTopologyTypes(alch_info%fep_topology) == "Dual") then

      ngroup = 2

      call setup_selection(sel_info, molecule)
      allocate(selatoms(ngroup), stat=alloc_stat)
      if (alloc_stat /= 0) &
        call error_msg('Define_FEP_Topology> Memory allocation error.')

      ! preserved atoms are set to 5
      allocate(molecule%fepgrp(molecule%num_atoms))
      molecule%fepgrp(:) = 5

      molecule%num_atoms_fep(1) = 0
      if (alch_info%singleA /= "NONE") then
        if (main_rank) then
          write(MsgOut,'(A)') 'Define_FEP_Topology> singleA is ignored &
            &in Dual topology.'
        end if
      end if

      molecule%num_atoms_fep(2) = 0
      if (alch_info%singleB /= "NONE") then
        if (main_rank) then
          write(MsgOut,'(A)') 'Define_FEP_Topology> singleB is ignored &
            &in Dual topology.'
        end if
      end if

      ! dualA atoms are set to 3
      molecule%num_atoms_fep(3) = 0
      if (alch_info%dualA /= "NONE") then
        read(alch_info%dualA, *) idx_select
        call select_atom(molecule, sel_info%groups(idx_select), &
          selatoms(1))
        natom = size(selatoms(1)%idx)
        do j = 1, natom
          do k = 1, molecule%num_atoms
            if(selatoms(1)%idx(j) == k) then
              molecule%fepgrp(k) = 3
            end if
          end do
        end do
        molecule%num_atoms_fep(3) = size(selatoms(1)%idx)
      end if

      ! dualB atoms are set to 4
      molecule%num_atoms_fep(4) = 0
      if (alch_info%dualB /= "NONE") then
        read(alch_info%dualB, *) idx_select
        call select_atom(molecule, sel_info%groups(idx_select), &
          selatoms(2))
        natom = size(selatoms(2)%idx)
        do j = 1, natom
          do k = 1, molecule%num_atoms
            if(selatoms(2)%idx(j) == k) then
              molecule%fepgrp(k) = 4
            end if
          end do
        end do
        molecule%num_atoms_fep(4) = size(selatoms(2)%idx)
      end if

      if ((alch_info%dualA == "NONE") .and. (alch_info%dualB == "NONE")) then
        call error_msg('Define_FEP_Topology> Atoms with dual &
          &topology are not defined.')
      end if

      molecule%num_atoms_fep(5) = molecule%num_atoms &
        - molecule%num_atoms_fep(3) &
        - molecule%num_atoms_fep(4)

      if (main_rank) then
        write(MsgOut,'(A)')'Define_FEP_Topology> Dual topology is &
          &assigned for FEP.'
        write(MsgOut,'(A)') ' '
      end if

    else

      ngroup = 4

      call setup_selection(sel_info, molecule)
      allocate(selatoms(ngroup), stat=alloc_stat)
      if (alloc_stat /= 0) &
        call error_msg('Setup_Restraints> Memory allocation error.')

      ! preserved atoms are set to 5
      allocate(molecule%fepgrp(molecule%num_atoms))
      molecule%fepgrp(:) = 5

      ! singleA atoms are set to 1
      molecule%num_atoms_fep(1) = 0
      if (alch_info%singleA /= "NONE") then
        read(alch_info%singleA, *) idx_select
        call select_atom(molecule, sel_info%groups(idx_select), &
          selatoms(1))
        natom = size(selatoms(1)%idx)
        do j = 1, natom
          do k = 1, molecule%num_atoms
            if(selatoms(1)%idx(j) == k) then
              molecule%fepgrp(k) = 1
            end if
          end do
        end do
        molecule%num_atoms_fep(1) = size(selatoms(1)%idx)
      else
        call error_msg('Define_FEP_Topology> singleA must be assigned &
          &in Hybrid topology.')
      end if

      ! singleB atoms are set to 2
      molecule%num_atoms_fep(2) = 0
      if (alch_info%singleB /= "NONE") then
        read(alch_info%singleB, *) idx_select
        call select_atom(molecule, sel_info%groups(idx_select), &
          selatoms(2))
        natom = size(selatoms(2)%idx)
        do j = 1, natom
          do k = 1, molecule%num_atoms
            if(selatoms(2)%idx(j) == k) then
              molecule%fepgrp(k) = 2
            end if
          end do
        end do
        molecule%num_atoms_fep(2) = size(selatoms(2)%idx)
      else
!        call error_msg('Define_FEP_Topology> singleB must be assigned &
!          &in Hybrid topology.')
      end if

      ! dualA atoms are set to 3
      molecule%num_atoms_fep(3) = 0
      if (alch_info%dualA /= "NONE") then
        read(alch_info%dualA, *) idx_select
        call select_atom(molecule, sel_info%groups(idx_select), &
          selatoms(3))
        natom = size(selatoms(3)%idx)
        do j = 1, natom
          do k = 1, molecule%num_atoms
            if(selatoms(3)%idx(j) == k) then
              molecule%fepgrp(k) = 3
            end if
          end do
        end do
        molecule%num_atoms_fep(3) = size(selatoms(3)%idx)
      end if

      ! dualB atoms are set to 4
      molecule%num_atoms_fep(4) = 0
      if (alch_info%dualB /= "NONE") then
        read(alch_info%dualB, *) idx_select
        call select_atom(molecule, sel_info%groups(idx_select), &
          selatoms(4))
        natom = size(selatoms(4)%idx)
        do j = 1, natom
          do k = 1, molecule%num_atoms
            if(selatoms(4)%idx(j) == k) then
              molecule%fepgrp(k) = 4
            end if
          end do
        end do
        molecule%num_atoms_fep(4) = size(selatoms(4)%idx)
      end if

      molecule%num_atoms_fep(5) = molecule%num_atoms &
        - molecule%num_atoms_fep(1) &
        - molecule%num_atoms_fep(2) &
        - molecule%num_atoms_fep(3) &
        - molecule%num_atoms_fep(4)
                                  
      if (main_rank) then
        write(MsgOut,'(A)')'Define_FEP_Topology> Hybrid topology &
          &is assigned for FEP.'
        write(MsgOut,'(A)') ' '
      end if

    end if

    ! make correspondence list of singleA and singleB
    call correspondence_single_fep(molecule)

    ! make fepgrp for bond, angle, dihedral, and cmap
    call make_bond_group_fep(molecule)

    ! delete A-B covalent bond
    call delete_bond(molecule)
    call delete_angle(molecule)
    call delete_dihedral(molecule)
    call delete_improper(molecule)
    call delete_cmap(molecule)

    ! make list for covalent bond
    call make_bondlist(molecule)
    call make_anglelist(molecule)
    call make_dihedrallist(molecule)
    call make_improperlist(molecule)
    call make_cmaplist(molecule)

    ! count hbonds in single topology
    call count_hbonds_single_fep(molecule)

    call check_parameter(molecule)

    ! write summary of molecule information
    !
    if (main_rank) then

      write(MsgOut,'(A)')'Define_FEP_Topology> Summary of molecules for FEP'

      write(MsgOut,'(A25,I10,A25,I10)')                  &
           '  num_atoms_singleA      = ', molecule%num_atoms_fep(1),     &
           '  num_atoms_singleB      = ', molecule%num_atoms_fep(2),     &
           '  num_atoms_dualA        = ', molecule%num_atoms_fep(3),     &
           '  num_atoms_dualB        = ', molecule%num_atoms_fep(4),     &
           '  num_atoms_preserve     = ', molecule%num_atoms_fep(5)
      write(MsgOut,'(A25,I10,A25,I10)')                  &
           '  num_bonds_singleA      = ', molecule%num_bonds_fep(1),     &
           '  num_bonds_singleB      = ', molecule%num_bonds_fep(2),     &
           '  num_bonds_dualA        = ', molecule%num_bonds_fep(3),     &
           '  num_bonds_dualB        = ', molecule%num_bonds_fep(4),     &
           '  num_bonds_preserve     = ', molecule%num_bonds_fep(5)
      write(MsgOut,'(A25,I10,A25,I10)')                  &
           '  num_angles_singleA     = ', molecule%num_angles_fep(1),    &
           '  num_angles_singleB     = ', molecule%num_angles_fep(2),    &
           '  num_angles_dualA       = ', molecule%num_angles_fep(3),    &
           '  num_angles_dualB       = ', molecule%num_angles_fep(4),    &
           '  num_angles_preserve    = ', molecule%num_angles_fep(5)
      write(MsgOut,'(A25,I10,A25,I10)')                  &
           '  num_dihedrals_singleA  = ', molecule%num_dihedrals_fep(1),  &
           '  num_dihedrals_singleB  = ', molecule%num_dihedrals_fep(2),  &
           '  num_dihedrals_dualA    = ', molecule%num_dihedrals_fep(3),  &
           '  num_dihedrals_dualB    = ', molecule%num_dihedrals_fep(4),  &
           '  num_dihedrals_preserve = ', molecule%num_dihedrals_fep(5)
      write(MsgOut,'(A25,I10,A25,I10)')                  &
           '  num_impropers_singleA  = ', molecule%num_impropers_fep(1),  &
           '  num_impropers_singleB  = ', molecule%num_impropers_fep(2),  &
           '  num_impropers_dualA    = ', molecule%num_impropers_fep(3),  &
           '  num_impropers_dualB    = ', molecule%num_impropers_fep(4),  &
           '  num_impropers_preserve = ', molecule%num_impropers_fep(5)
      write(MsgOut,'(A25,I10,A25,I10)')                  &
           '  num_cmaps_singleA      = ', molecule%num_cmaps_fep(1), &
           '  num_cmaps_singleB      = ', molecule%num_cmaps_fep(2), &
           '  num_cmaps_dualA        = ', molecule%num_cmaps_fep(3), &
           '  num_cmaps_dualB        = ', molecule%num_cmaps_fep(4), &
           '  num_cmaps_preserve     = ', molecule%num_cmaps_fep(5)
      write(MsgOut,'(A25,I10,A25,I10)')                    &
           '  num_residues           = ', molecule%num_residues,           &
           '  num_molecules          = ', molecule%num_molecules
      write(MsgOut,'(A25,I10,A25,I10)')                    &
           '  num_segments           = ', molecule%num_segments,           &
           '  num_deg_freedom        = ', molecule%num_deg_freedom
      write(MsgOut,'(A25,F10.3)')                          &
           '  total_charge           = ', molecule%total_charge
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine define_fep_topology

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine delete_bond(molecule)

    ! formal arguments
    type(s_molecule),         intent(inout) :: molecule

    ! local variables
    integer                  :: i
    integer                  :: nbnd
    integer                  :: nbnd_fep
    integer                  :: nbnd_del
    integer                  :: idx1
    integer                  :: idx2

    integer,          allocatable :: bond_list_del(:,:)
    integer,          allocatable :: bond_list_fep(:,:)

    nbnd = molecule%num_bonds

    allocate(bond_list_del(2,nbnd))
    allocate(bond_list_fep(2,nbnd))

    nbnd_fep = 0
    nbnd_del = 0
    do i = 1, nbnd
      idx1 = int(molecule%fepgrp(molecule%bond_list(1, i)))
      idx2 = int(molecule%fepgrp(molecule%bond_list(2, i)))

      if (molecule%fepgrp_bond(idx1,idx2) == 0) then
        nbnd_del = nbnd_del + 1
        bond_list_del(:, nbnd_del) = molecule%bond_list(:, i)
      else
        nbnd_fep = nbnd_fep + 1
        bond_list_fep(:, nbnd_fep) = molecule%bond_list(:, i)
      end if

    end do

    deallocate(molecule%bond_list)
    allocate(molecule%bond_list(2, nbnd_fep))

    do i = 1, nbnd_fep
      molecule%bond_list(:, i) = bond_list_fep(:, i)
    end do

    ! Write delete Bonds
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Delete Bonds> Summary of delete bond lists'
      write(MsgOut,'(A,I10,A)') 'Delete ', nbnd_del, ' bonds'
      write(MsgOut,'(8I8)') (bond_list_del(1, i), &
                             bond_list_del(2, i), i = 1, nbnd_del)
      write(MsgOut,'(A)') ' '
    end if

    deallocate(bond_list_fep)
    deallocate(bond_list_del)

    molecule%num_bonds = nbnd_fep

    return

  end subroutine delete_bond

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine delete_angle(molecule)

    ! formal arguments
    type(s_molecule),         intent(inout) :: molecule

    ! local variables
    integer                  :: i
    integer                  :: nangl
    integer                  :: nangl_fep
    integer                  :: nangl_del
    integer                  :: idx1
    integer                  :: idx2
    integer                  :: idx3

    integer,          allocatable :: angl_list_del(:,:)
    integer,          allocatable :: angl_list_fep(:,:)


    ! Delete Angles contain vanishing atom and appearing atom
    !
    nangl = molecule%num_angles

    allocate(angl_list_del(3,nangl))
    allocate(angl_list_fep(3,nangl))

    nangl_fep = 0
    nangl_del = 0
    do i = 1, nangl
      idx1 = int(molecule%fepgrp(molecule%angl_list(1, i)))
      idx2 = int(molecule%fepgrp(molecule%angl_list(2, i)))
      idx3 = int(molecule%fepgrp(molecule%angl_list(3, i)))

      if (molecule%fepgrp_angl(idx1,idx2,idx3) == 0) then
        nangl_del = nangl_del + 1
        angl_list_del(:, nangl_del) = molecule%angl_list(:, i)
      else
        nangl_fep = nangl_fep + 1
        angl_list_fep(:, nangl_fep) = molecule%angl_list(:, i)
      end if

    end do

    deallocate(molecule%angl_list)
    allocate(molecule%angl_list(3, nangl_fep))

    do i = 1, nangl_fep
      molecule%angl_list(:, i) = angl_list_fep(:, i)
    end do

    ! Write delete Angles
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Delete Angles> Summary of delete angle lists'
      write(MsgOut,'(A,I10,A)') 'Delete ', nangl_del, ' angles'
      write(MsgOut,'(9I8)') (angl_list_del(1, i), &
                             angl_list_del(2, i), &
                             angl_list_del(3, i), i = 1, nangl_del)
      write(MsgOut,'(A)') ' '
    end if

    deallocate(angl_list_fep)
    deallocate(angl_list_del)

    molecule%num_angles = nangl_fep

    return

  end subroutine delete_angle

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine delete_dihedral(molecule)

    ! formal arguments
    type(s_molecule),         intent(inout) :: molecule

    ! local variables
    integer                  :: i
    integer                  :: ndihe
    integer                  :: ndihe_fep
    integer                  :: ndihe_del
    integer                  :: idx1
    integer                  :: idx2
    integer                  :: idx3
    integer                  :: idx4

    integer,          allocatable :: dihe_list_del(:,:)
    integer,          allocatable :: dihe_list_fep(:,:)


    ! Delete Dihedrals contain vanishing atom and appearing atom
    !
    ndihe = molecule%num_dihedrals

    allocate(dihe_list_del(4,ndihe))
    allocate(dihe_list_fep(4,ndihe))

    ndihe_fep = 0
    ndihe_del = 0
    do i = 1, ndihe
      idx1 = int(molecule%fepgrp(molecule%dihe_list(1, i)))
      idx2 = int(molecule%fepgrp(molecule%dihe_list(2, i)))
      idx3 = int(molecule%fepgrp(molecule%dihe_list(3, i)))
      idx4 = int(molecule%fepgrp(molecule%dihe_list(4, i)))

      if (molecule%fepgrp_dihe(idx1,idx2,idx3,idx4) == 0) then
        ndihe_del = ndihe_del + 1
        dihe_list_del(:, ndihe_del) = molecule%dihe_list(:, i)
      else
        ndihe_fep = ndihe_fep + 1
        dihe_list_fep(:, ndihe_fep) = molecule%dihe_list(:, i)
      end if

    end do

    deallocate(molecule%dihe_list)
    allocate(molecule%dihe_list(4, ndihe_fep))

    do i = 1, ndihe_fep
      molecule%dihe_list(:, i) = dihe_list_fep(:, i)
    end do

    ! Write delete Dihedrals
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Delete Dihedrals> Summary of delete dihedral lists'
      write(MsgOut,'(A,I10,A)') 'Delete ', ndihe_del, ' dihedrals'
      write(MsgOut,'(8I8)') (dihe_list_del(1, i), &
                             dihe_list_del(2, i), &
                             dihe_list_del(3, i), &
                             dihe_list_del(4, i), i = 1, ndihe_del)
      write(MsgOut,'(A)') ' '
    end if

    deallocate(dihe_list_fep)
    deallocate(dihe_list_del)

    molecule%num_dihedrals = ndihe_fep

    return

  end subroutine delete_dihedral

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine delete_improper(molecule)

    ! formal arguments
    type(s_molecule),         intent(inout) :: molecule

    ! local variables
    integer                  :: i
    integer                  :: nimpr
    integer                  :: nimpr_fep
    integer                  :: nimpr_del
    integer                  :: idx1
    integer                  :: idx2
    integer                  :: idx3
    integer                  :: idx4

    integer,          allocatable :: impr_list_del(:,:)
    integer,          allocatable :: impr_list_fep(:,:)


    ! Delete Impropers contain vanishing atom and appearing atom
    !
    nimpr = molecule%num_impropers

    allocate(impr_list_del(4,nimpr))
    allocate(impr_list_fep(4,nimpr))

    nimpr_fep = 0
    nimpr_del = 0
    do i = 1, nimpr
      idx1 = int(molecule%fepgrp(molecule%impr_list(1, i)))
      idx2 = int(molecule%fepgrp(molecule%impr_list(2, i)))
      idx3 = int(molecule%fepgrp(molecule%impr_list(3, i)))
      idx4 = int(molecule%fepgrp(molecule%impr_list(4, i)))

      if (molecule%fepgrp_dihe(idx1,idx2,idx3,idx4) == 0) then
        nimpr_del = nimpr_del + 1
        impr_list_del(:, nimpr_del) = molecule%impr_list(:, i)
      else
        nimpr_fep = nimpr_fep + 1
        impr_list_fep(:, nimpr_fep) = molecule%impr_list(:, i)
      end if

    end do

    deallocate(molecule%impr_list)
    allocate(molecule%impr_list(4, nimpr_fep))

    do i = 1, nimpr_fep
      molecule%impr_list(:, i) = impr_list_fep(:, i)
    end do

    ! Write delete Impropers
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Delete Impropers> Summary of delete improper lists'
      write(MsgOut,'(A,I10,A)') 'Delete ', nimpr_del, ' impropers'
      write(MsgOut,'(8I8)') (impr_list_del(1, i), &
                             impr_list_del(2, i), &
                             impr_list_del(3, i), &
                             impr_list_del(4, i), i = 1, nimpr_del)
      write(MsgOut,'(A)') ' '
    end if

    deallocate(impr_list_fep)
    deallocate(impr_list_del)

    molecule%num_impropers = nimpr_fep

    return

  end subroutine delete_improper

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine delete_cmap(molecule)

    ! formal arguments
    type(s_molecule),         intent(inout) :: molecule

    ! local variables
    integer                  :: i
    integer                  :: ncmap
    integer                  :: ncmap_fep
    integer                  :: ncmap_del
    integer                  :: i1, i2, i3, i4
    integer                  :: i5, i6, i7, i8
    integer                  :: idx

    integer,          allocatable :: cmap_list_del(:,:)
    integer,          allocatable :: cmap_list_fep(:,:)

    if (.not. allocated(molecule%cmap_list)) return

    ! Delete Cmaps contain vanishing atom and appearing atom
    !
    ncmap = molecule%num_cmaps

    allocate(cmap_list_del(8,ncmap))
    allocate(cmap_list_fep(8,ncmap))

    ncmap_del = 0
    ncmap_fep = 0
    do i = 1, ncmap
      i1 = int(molecule%fepgrp(molecule%cmap_list(1, i)))
      i2 = int(molecule%fepgrp(molecule%cmap_list(2, i)))
      i3 = int(molecule%fepgrp(molecule%cmap_list(3, i)))
      i4 = int(molecule%fepgrp(molecule%cmap_list(4, i)))
      i5 = int(molecule%fepgrp(molecule%cmap_list(5, i)))
      i6 = int(molecule%fepgrp(molecule%cmap_list(6, i)))
      i7 = int(molecule%fepgrp(molecule%cmap_list(7, i)))
      i8 = int(molecule%fepgrp(molecule%cmap_list(8, i)))

      idx = i1 + 5*(i2-1 + 5*(i3-1 + 5*(i4-1 + 5*(i5-1 + 5*(i6-1 + 5*(i7-1 + 5*(i8-1)))))))

      if (molecule%fepgrp_cmap(idx) == 0) then
        ncmap_del = ncmap_del + 1
        cmap_list_del(:, ncmap_del) = molecule%cmap_list(:, i)
      else
        ncmap_fep = ncmap_fep + 1
        cmap_list_fep(:, ncmap_fep) = molecule%cmap_list(:, i)
      end if

    end do

    deallocate(molecule%cmap_list)
    allocate(molecule%cmap_list(8, ncmap_fep))

    do i = 1, ncmap_fep
      molecule%cmap_list(:, i) = cmap_list_fep(:, i)
    end do

    ! Write delete Cmaps
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Delete Cmaps> Summary of delete cmap lists'
      write(MsgOut,'(A,I10,A)') 'Delete ', ncmap_del, ' cmaps'
      write(MsgOut,'(8I8)') (cmap_list_del(1, i), &
                             cmap_list_del(2, i), &
                             cmap_list_del(3, i), &
                             cmap_list_del(4, i), &
                             cmap_list_del(5, i), &
                             cmap_list_del(6, i), &
                             cmap_list_del(7, i), &
                             cmap_list_del(8, i), i = 1, ncmap_del)
      write(MsgOut,'(A)') ' '
    end if

    deallocate(cmap_list_fep)
    deallocate(cmap_list_del)

    molecule%num_cmaps = ncmap_fep

    return

  end subroutine delete_cmap

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine make_bondlist(molecule)

    ! formal arguments
    type(s_molecule),         intent(inout) :: molecule

    ! local variables
    integer                  :: i, j, k
    integer                  :: nbnd
    integer                  :: nbnd_fep(5)
    integer                  :: nbnd_fep_max
    integer                  :: idx1, idx2
    integer,     allocatable :: bond_list_fep(:,:,:)

    ! Make bondlists for FEP calculation
    !
    nbnd = molecule%num_bonds

    allocate(bond_list_fep(2,nbnd,5))

    do k = 1, 5
      nbnd_fep(k) = 0
    end do

    ! Assign each dihedral to hybrid topology
    do i = 1, nbnd
      idx1 = int(molecule%fepgrp(molecule%bond_list(1, i)))
      idx2 = int(molecule%fepgrp(molecule%bond_list(2, i)))
      do k = 1, 5
        if (molecule%fepgrp_bond(idx1,idx2) == k) then
          nbnd_fep(k) = nbnd_fep(k) + 1
          bond_list_fep(:, nbnd_fep(k), k) = molecule%bond_list(:, i)
        end if
      end do
    end do

    nbnd_fep_max = 0
    do k = 1, 5
      if (nbnd_fep(k) > nbnd_fep_max) &
        nbnd_fep_max = nbnd_fep(k)
    end do
    allocate(molecule%bond_list_fep(2, nbnd_fep_max, 5))

    do k = 1, 5
      do i = 1, nbnd_fep(k)
        molecule%bond_list_fep(:, i, k) = bond_list_fep(:, i, k)
      end do
    end do

    deallocate(bond_list_fep)

    do k = 1, 5
      molecule%num_bonds_fep(k) = nbnd_fep(k)
    end do

    return

  end subroutine make_bondlist

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine make_anglelist(molecule)

    ! formal arguments
    type(s_molecule),         intent(inout) :: molecule

    ! local variables
    integer                  :: i, j, k
    integer                  :: nangl
    integer                  :: nangl_fep(5)
    integer                  :: nangl_fep_max
    integer                  :: idx1, idx2, idx3
    integer,     allocatable :: angl_list_fep(:,:,:)

    ! Make anglelists for FEP calculation
    !
    nangl = molecule%num_angles

    allocate(angl_list_fep(3,nangl,5))

    do k = 1, 5
      nangl_fep(k) = 0
    end do

    ! Assign each dihedral to hybrid topology
    do i = 1, nangl
      idx1 = int(molecule%fepgrp(molecule%angl_list(1, i)))
      idx2 = int(molecule%fepgrp(molecule%angl_list(2, i)))
      idx3 = int(molecule%fepgrp(molecule%angl_list(3, i)))
      do k = 1, 5
        if (molecule%fepgrp_angl(idx1,idx2,idx3) == k) then
          nangl_fep(k) = nangl_fep(k) + 1
          angl_list_fep(:, nangl_fep(k), k) = molecule%angl_list(:, i)
        end if
      end do
    end do

    nangl_fep_max = 0
    do k = 1, 5
      if (nangl_fep(k) > nangl_fep_max) &
        nangl_fep_max = nangl_fep(k)
    end do
    allocate(molecule%angl_list_fep(3, nangl_fep_max, k))

    do k = 1, 5
      do i = 1, nangl_fep(k)
        molecule%angl_list_fep(:, i, k) = angl_list_fep(:, i, k)
      end do
    end do

    deallocate(angl_list_fep)

    do k = 1, 5
      molecule%num_angles_fep(k) = nangl_fep(k)
    end do

    return

  end subroutine make_anglelist

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine make_dihedrallist(molecule)

    ! formal arguments
    type(s_molecule),         intent(inout) :: molecule

    ! local variables
    integer                  :: i, j, k, l
    integer                  :: ndihe
    integer                  :: ndihe_fep(5)
    integer                  :: ndihe_fep_max
    integer                  :: idx1, idx2, idx3, idx4
    integer,     allocatable :: dihe_list_fep(:,:,:)

    ! Make dihedrallists for FEP calculation
    !
    ndihe = molecule%num_dihedrals

    allocate(dihe_list_fep(4,ndihe,5))

    do k = 1, 5
      ndihe_fep(k) = 0
    end do

    ! Assign each dihedral to hybrid topology
    do i = 1, ndihe
      idx1 = int(molecule%fepgrp(molecule%dihe_list(1, i)))
      idx2 = int(molecule%fepgrp(molecule%dihe_list(2, i)))
      idx3 = int(molecule%fepgrp(molecule%dihe_list(3, i)))
      idx4 = int(molecule%fepgrp(molecule%dihe_list(4, i)))
      do k = 1, 5
        if (molecule%fepgrp_dihe(idx1,idx2,idx3,idx4) == k) then
          ndihe_fep(k) = ndihe_fep(k) + 1
          dihe_list_fep(:, ndihe_fep(k), k) = molecule%dihe_list(:, i)
        end if
      end do
    end do

    ndihe_fep_max = 0
    do k = 1, 5
      if (ndihe_fep(k) > ndihe_fep_max) &
        ndihe_fep_max = ndihe_fep(k)
    end do
    allocate(molecule%dihe_list_fep(4, ndihe_fep_max, 5))


    do k = 1, 5
      do i = 1, ndihe_fep(k)
        molecule%dihe_list_fep(:, i, k) = dihe_list_fep(:, i, k)
      end do
    end do

    deallocate(dihe_list_fep)

    do k = 1, 5
      molecule%num_dihedrals_fep(k) = ndihe_fep(k)
    end do

    return

  end subroutine make_dihedrallist

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine make_improperlist(molecule)

    ! formal arguments
    type(s_molecule),         intent(inout) :: molecule

    ! local variables
    integer                  :: i, j, k, l
    integer                  :: nimpr
    integer                  :: nimpr_fep(5)
    integer                  :: nimpr_fep_max
    integer                  :: idx1, idx2, idx3, idx4
    integer,          allocatable :: impr_list_fep(:,:,:)

    ! Make improperlists for FEP calculation
    !
    nimpr = molecule%num_impropers

    allocate(impr_list_fep(4,nimpr,5))

    do k = 1, 5
      nimpr_fep(k) = 0
    end do

    ! Assign each dihedral to hybrid topology
    do i = 1, nimpr
      idx1 = int(molecule%fepgrp(molecule%impr_list(1, i)))
      idx2 = int(molecule%fepgrp(molecule%impr_list(2, i)))
      idx3 = int(molecule%fepgrp(molecule%impr_list(3, i)))
      idx4 = int(molecule%fepgrp(molecule%impr_list(4, i)))
      do k = 1, 5
        if (molecule%fepgrp_dihe(idx1,idx2,idx3,idx4) == k) then
          nimpr_fep(k) = nimpr_fep(k) + 1
          impr_list_fep(:, nimpr_fep(k), k) = molecule%impr_list(:, i)
        end if
      end do
    end do

    nimpr_fep_max = 0
    do k = 1, 5
      if (nimpr_fep(k) > nimpr_fep_max) &
        nimpr_fep_max = nimpr_fep(k)
    end do
    allocate(molecule%impr_list_fep(4, nimpr_fep_max, 5))

    do k = 1, 5
      do i = 1, nimpr_fep(k)
        molecule%impr_list_fep(:, i, k) = impr_list_fep(:, i, k)
      end do
    end do

    deallocate(impr_list_fep)

    do k = 1, 5
      molecule%num_impropers_fep(k) = nimpr_fep(k)
    end do

    return

  end subroutine make_improperlist

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine make_cmaplist(molecule)

    ! formal arguments
    type(s_molecule),         intent(inout) :: molecule

    ! local variables
    integer                  :: i, k
    integer                  :: i1, i2, i3, i4
    integer                  :: i5, i6, i7, i8
    integer                  :: idx
    integer                  :: ncmap
    integer                  :: ncmap_fep(5)
    integer                  :: ncmap_fep_max
    integer,          allocatable :: cmap_list_fep(:,:,:)

    if (.not. allocated(molecule%cmap_list)) return

    ! Make cmaplists for FEP calculation
    !
    ncmap = molecule%num_cmaps

    allocate(cmap_list_fep(8,ncmap,5))

    do k = 1, 5
      ncmap_fep(k) = 0
    end do

    ! Assign each dihedral to hybrid topology
    do i = 1, ncmap
      i1 = int(molecule%fepgrp(molecule%cmap_list(1, i)))
      i2 = int(molecule%fepgrp(molecule%cmap_list(2, i)))
      i3 = int(molecule%fepgrp(molecule%cmap_list(3, i)))
      i4 = int(molecule%fepgrp(molecule%cmap_list(4, i)))
      i5 = int(molecule%fepgrp(molecule%cmap_list(5, i)))
      i6 = int(molecule%fepgrp(molecule%cmap_list(6, i)))
      i7 = int(molecule%fepgrp(molecule%cmap_list(7, i)))
      i8 = int(molecule%fepgrp(molecule%cmap_list(8, i)))
      idx = i1 + 5*(i2-1 + 5*(i3-1 + 5*(i4-1 + 5*(i5-1 + 5*(i6-1 + 5*(i7-1 + 5*(i8-1)))))))
      do k = 1, 5
        if (molecule%fepgrp_cmap(idx) == k) then
          ncmap_fep(k) = ncmap_fep(k) + 1
          cmap_list_fep(:, ncmap_fep(k), k) = molecule%cmap_list(:, i)
        end if
      end do
    end do

    ncmap_fep_max = 0
    do k = 1, 5
      if (ncmap_fep(k) > ncmap_fep_max) &
        ncmap_fep_max = ncmap_fep(k)
    end do
    allocate(molecule%cmap_list_fep(8, ncmap_fep_max, 5))

    do k = 1, 5
      do i = 1, ncmap_fep(k)
        molecule%cmap_list_fep(:, i, k) = cmap_list_fep(:, i, k)
      end do
    end do

    deallocate(cmap_list_fep)

    do k = 1, 5
      molecule%num_cmaps_fep(k) = ncmap_fep(k)
    end do

    return

  end subroutine make_cmaplist

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_parameter(molecule)

    ! formal arguments
    type(s_molecule),         intent(inout) :: molecule

    ! local variables
    integer                  :: i
    integer                  :: nbnd
    integer                  :: nangl
    integer                  :: ndihe
    integer                  :: nimpr
    integer                  :: nerr
    integer                  :: nerr_bond
    integer                  :: nerr_angl
    integer                  :: nerr_dihe
    integer                  :: nerr_impr
    integer                  :: idx1
    integer                  :: idx2
    integer                  :: idx3
    integer                  :: idx4

    integer,          allocatable :: bond_list_del(:,:)
    integer,          allocatable :: bond_list_fep(:,:)

    nbnd  = molecule%num_bonds
    nangl = molecule%num_angles
    ndihe = molecule%num_dihedrals
    nimpr = molecule%num_impropers

    nerr = 0
    nerr_bond = 0
    nerr_angl = 0
    nerr_dihe = 0
    nerr_impr = 0
    do i = 1, nbnd
      idx1 = int(molecule%fepgrp(molecule%bond_list(1, i)))
      idx2 = int(molecule%fepgrp(molecule%bond_list(2, i)))

     select case(idx1)

     case(1)
       if((idx2 == 2) .or. (idx2 == 4)) then
          nerr = nerr + 1
          nerr_bond = nerr_bond + 1
          cycle
       end if

     case(2)
       if((idx2 == 1) .or. (idx2 == 3)) then
          nerr = nerr + 1
          nerr_bond = nerr_bond + 1
          cycle
       end if

     case(3)
       if(idx2 == 2) then
          nerr = nerr + 1
          nerr_bond = nerr_bond + 1
          cycle
       end if

     case(4)
       if(idx2 == 1) then
          nerr = nerr + 1
          nerr_bond = nerr_bond + 1
          cycle
       end if

     end select

    end do

    do i = 1, nangl
      idx1 = int(molecule%fepgrp(molecule%angl_list(1, i)))
      idx2 = int(molecule%fepgrp(molecule%angl_list(2, i)))
      idx3 = int(molecule%fepgrp(molecule%angl_list(3, i)))

     select case(idx1)

     case(1)
       if((idx2 == 2) .or. (idx2 == 4) .or. &
          (idx3 == 2) .or. (idx3 == 4)) then
            nerr = nerr + 1
            nerr_angl = nerr_angl + 1
            cycle
       end if

     case(2)
       if((idx2 == 1) .or. (idx2 == 3) .or. &
          (idx3 == 1) .or. (idx3 == 3)) then
            nerr = nerr + 1
            nerr_angl = nerr_angl + 1
            cycle
       end if

     case(3)
       if((idx2 == 2) .or. (idx3 == 2)) then
            nerr = nerr + 1
            nerr_angl = nerr_angl + 1
            cycle
       end if

     case(4)
       if((idx2 == 1) .or. (idx3 == 1)) then
            nerr = nerr + 1
            nerr_angl = nerr_angl + 1
            cycle
       end if

     end select


     select case(idx2)

     case(1)
       if((idx3 == 2) .or. (idx3 == 4)) then
         nerr = nerr + 1
         nerr_angl = nerr_angl + 1
         cycle
       end if

     case(2)
       if((idx3 == 1) .or. (idx3 == 3)) then
         nerr = nerr + 1
         nerr_angl = nerr_angl + 1
         cycle
       end if

     case(3)
       if(idx3 == 2) then
         nerr = nerr + 1
         nerr_angl = nerr_angl + 1
         cycle
       end if

     case(4)
       if(idx3 == 1) then
         nerr = nerr + 1
         nerr_angl = nerr_angl + 1
         cycle
       end if

     end select

    end do

    do i = 1, ndihe
      idx1 = int(molecule%fepgrp(molecule%dihe_list(1, i)))
      idx2 = int(molecule%fepgrp(molecule%dihe_list(2, i)))
      idx3 = int(molecule%fepgrp(molecule%dihe_list(3, i)))
      idx4 = int(molecule%fepgrp(molecule%dihe_list(4, i)))

     select case(idx1)

     case(1)
       if((idx2 == 2) .or. (idx2 == 4) .or. &
          (idx3 == 2) .or. (idx3 == 4) .or. &
          (idx4 == 2) .or. (idx4 == 4)) then
            nerr = nerr + 1
            nerr_dihe = nerr_dihe + 1
            cycle
       end if

     case(2)
       if((idx2 == 1) .or. (idx2 == 3) .or. &
          (idx3 == 1) .or. (idx3 == 3) .or. &
          (idx4 == 1) .or. (idx4 == 3)) then
            nerr = nerr + 1
            nerr_dihe = nerr_dihe + 1
            cycle
       end if

     case(3)
       if((idx2 == 2) .or. (idx3 == 2) .or. (idx4 == 2)) then
            nerr = nerr + 1
            nerr_dihe = nerr_dihe + 1
            cycle
       end if

     case(4)
       if((idx2 == 1) .or. (idx3 == 1) .or. (idx4 == 1)) then
            nerr = nerr + 1
            nerr_dihe = nerr_dihe + 1
            cycle
       end if

     end select


     select case(idx2)

     case(1)
       if((idx3 == 2) .or. (idx3 == 4) .or. &
          (idx4 == 2) .or. (idx4 == 4)) then
            nerr = nerr + 1
            nerr_dihe = nerr_dihe + 1
            cycle
       end if

     case(2)
       if((idx3 == 1) .or. (idx3 == 3) .or. &
          (idx4 == 1) .or. (idx4 == 3)) then
            nerr = nerr + 1
            nerr_dihe = nerr_dihe + 1
            cycle
       end if

     case(3)
       if((idx3 == 2) .or. (idx4 == 2)) then
            nerr = nerr + 1
            nerr_dihe = nerr_dihe + 1
            cycle
       end if

     case(4)
       if((idx3 == 1) .or. (idx4 == 1)) then
            nerr = nerr + 1
            nerr_dihe = nerr_dihe + 1
            cycle
       end if

     end select


     select case(idx3)

     case(1)
       if((idx4 == 2) .or. (idx4 == 4)) then
            nerr = nerr + 1
            nerr_dihe = nerr_dihe + 1
            cycle
       end if

     case(2)
       if((idx4 == 1) .or. (idx4 == 3)) then
            nerr = nerr + 1
            nerr_dihe = nerr_dihe + 1
            cycle
       end if

     case(3)
       if(idx4 == 2) then
            nerr = nerr + 1
            nerr_dihe = nerr_dihe + 1
            cycle
       end if

     case(4)
       if(idx4 == 1) then
            nerr = nerr + 1
            nerr_dihe = nerr_dihe + 1
            cycle
       end if

     end select

    end do

    do i = 1, nimpr
      idx1 = int(molecule%fepgrp(molecule%impr_list(1, i)))
      idx2 = int(molecule%fepgrp(molecule%impr_list(2, i)))
      idx3 = int(molecule%fepgrp(molecule%impr_list(3, i)))
      idx4 = int(molecule%fepgrp(molecule%impr_list(4, i)))

     select case(idx1)

     case(1)
       if((idx2 == 2) .or. (idx2 == 4) .or. &
          (idx3 == 2) .or. (idx3 == 4) .or. &
          (idx4 == 2) .or. (idx4 == 4)) then
            nerr = nerr + 1
            nerr_impr = nerr_impr + 1
            cycle
       end if

     case(2)
       if((idx2 == 1) .or. (idx2 == 3) .or. &
          (idx3 == 1) .or. (idx3 == 3) .or. &
          (idx4 == 1) .or. (idx4 == 3)) then
            nerr = nerr + 1
            nerr_impr = nerr_impr + 1
            cycle
       end if

     case(3)
       if((idx2 == 2) .or. (idx3 == 2) .or. (idx4 == 2)) then
            nerr = nerr + 1
            nerr_impr = nerr_impr + 1
            cycle
       end if

     case(4)
       if((idx2 == 1) .or. (idx3 == 1) .or. (idx4 == 1)) then
            nerr = nerr + 1
            nerr_impr = nerr_impr + 1
            cycle
       end if

     end select


     select case(idx2)

     case(1)
       if((idx3 == 2) .or. (idx3 == 4) .or. &
          (idx4 == 2) .or. (idx4 == 4)) then
            nerr = nerr + 1
            nerr_impr = nerr_impr + 1
            cycle
       end if

     case(2)
       if((idx3 == 1) .or. (idx3 == 3) .or. &
          (idx4 == 1) .or. (idx4 == 3)) then
            nerr = nerr + 1
            nerr_impr = nerr_impr + 1
            cycle
       end if

     case(3)
       if((idx3 == 2) .or. (idx4 == 2)) then
            nerr = nerr + 1
            nerr_impr = nerr_impr + 1
            cycle
       end if

     case(4)
       if((idx3 == 1) .or. (idx4 == 1)) then
            nerr = nerr + 1
            nerr_impr = nerr_impr + 1
            cycle
       end if

     end select


     select case(idx3)

     case(1)
       if((idx4 == 2) .or. (idx4 == 4)) then
            nerr = nerr + 1
            nerr_impr = nerr_impr + 1
            cycle
       end if

     case(2)
       if((idx4 == 1) .or. (idx4 == 3)) then
            nerr = nerr + 1
            nerr_impr = nerr_impr + 1
            cycle
       end if

     case(3)
       if(idx4 == 2) then
            nerr = nerr + 1
            nerr_impr = nerr_impr + 1
            cycle
       end if

     case(4)
       if(idx4 == 1) then
            nerr = nerr + 1
            nerr_impr = nerr_impr + 1
            cycle
       end if

     end select

    end do

    if(nerr /= 0) then
      if (main_rank) write(MsgOut,*) "Check Hybrid Topology> Inappropriate lists"
      if(nerr_bond /= 0) then
        if (main_rank) write(MsgOut,*) "Inappropriate bond lists>"
        do i = 1, nbnd
          idx1 = int(molecule%fepgrp(molecule%bond_list(1, i)))
          idx2 = int(molecule%fepgrp(molecule%bond_list(2, i)))


         select case(idx1)

         case(1)
            if((idx2 == 2) .or. (idx2 == 4)) then
              if (main_rank) write(MsgOut,*) molecule%bond_list(:, i)
              cycle
            end if

          case(2)
            if((idx2 == 1) .or. (idx2 == 3)) then
               if (main_rank) write(MsgOut,*) molecule%bond_list(:, i)
               cycle
            end if

          case(3)
            if(idx2 == 2) then
              if (main_rank) write(MsgOut,*) molecule%bond_list(:, i)
              cycle
            end if

          case(4)
            if(idx2 == 1) then
               if (main_rank) write(MsgOut,*) molecule%bond_list(:, i)
               cycle
            end if

          end select

        end do
      end if

      if(nerr_angl /= 0) then
        if (main_rank) write(MsgOut,*) "Inappropriate angle lists>"
        do i = 1, nangl
          idx1 = int(molecule%fepgrp(molecule%angl_list(1, i)))
          idx2 = int(molecule%fepgrp(molecule%angl_list(2, i)))
          idx3 = int(molecule%fepgrp(molecule%angl_list(3, i)))


          select case(idx1)

          case(1)
            if((idx2 == 2) .or. (idx2 == 4) .or. &
               (idx3 == 2) .or. (idx3 == 4)) then
                 if (main_rank) write(MsgOut,*) molecule%angl_list(:, i)
                 cycle
            end if

          case(2)
            if((idx2 == 1) .or. (idx2 == 3) .or. &
               (idx3 == 1) .or. (idx3 == 3)) then
                 if (main_rank) write(MsgOut,*) molecule%angl_list(:, i)
                 cycle
            end if

          case(3)
            if((idx2 == 2) .or. (idx3 == 2)) then
                 if (main_rank) write(MsgOut,*) molecule%angl_list(:, i)
                 cycle
            end if

          case(4)
            if((idx2 == 1) .or. (idx3 == 1)) then
                 if (main_rank) write(MsgOut,*) molecule%angl_list(:, i)
                 cycle
            end if

          end select


          select case(idx2)

          case(1)
            if((idx3 == 2) .or. (idx3 == 4)) then
              if (main_rank) write(MsgOut,*) molecule%angl_list(:, i)
              cycle
            end if

          case(2)
            if((idx3 == 1) .or. (idx3 == 3)) then
              if (main_rank) write(MsgOut,*) molecule%angl_list(:, i)
              cycle
            end if

          case(3)
            if(idx3 == 2) then
              if (main_rank) write(MsgOut,*) molecule%angl_list(:, i)
              cycle
            end if

          case(4)
            if(idx3 == 1) then
              if (main_rank) write(MsgOut,*) molecule%angl_list(:, i)
              cycle
            end if

          end select

        end do

      end if

      if(nerr_dihe /= 0) then
        if (main_rank) write(MsgOut,*) "Inappropriate dihedral lists>"
        do i = 1, ndihe
          idx1 = int(molecule%fepgrp(molecule%dihe_list(1, i)))
          idx2 = int(molecule%fepgrp(molecule%dihe_list(2, i)))
          idx3 = int(molecule%fepgrp(molecule%dihe_list(3, i)))
          idx4 = int(molecule%fepgrp(molecule%dihe_list(4, i)))


          select case(idx1)

          case(1)
            if((idx2 == 2) .or. (idx2 == 4) .or. &
               (idx3 == 2) .or. (idx3 == 4) .or. &
               (idx4 == 2) .or. (idx4 == 4)) then
                 if (main_rank) write(MsgOut,*) molecule%dihe_list(:, i)
                 cycle
            end if

          case(2)
            if((idx2 == 1) .or. (idx2 == 3) .or. &
               (idx3 == 1) .or. (idx3 == 3) .or. &
               (idx4 == 1) .or. (idx4 == 3)) then
                 if (main_rank) write(MsgOut,*) molecule%dihe_list(:, i)
                 cycle
            end if

          case(3)
            if((idx2 == 2) .or. (idx3 == 2) .or. (idx4 == 2)) then
                 if (main_rank) write(MsgOut,*) molecule%dihe_list(:, i)
                 cycle
            end if

          case(4)
            if((idx2 == 1) .or. (idx3 == 1) .or. (idx4 == 1)) then
                 if (main_rank) write(MsgOut,*) molecule%dihe_list(:, i)
                 cycle
            end if

          end select


          select case(idx2)

          case(1)
            if((idx3 == 2) .or. (idx3 == 4) .or. &
               (idx4 == 2) .or. (idx4 == 4)) then
                 if (main_rank) write(MsgOut,*) molecule%dihe_list(:, i)
                 cycle
            end if

          case(2)
            if((idx3 == 1) .or. (idx3 == 3) .or. &
               (idx4 == 1) .or. (idx4 == 3)) then
                 if (main_rank) write(MsgOut,*) molecule%dihe_list(:, i)
                 cycle
            end if

          case(3)
            if((idx3 == 2) .or. (idx4 == 2)) then
                 if (main_rank) write(MsgOut,*) molecule%dihe_list(:, i)
                 cycle
            end if

          case(4)
            if((idx3 == 1) .or. (idx4 == 1)) then
                 if (main_rank) write(MsgOut,*) molecule%dihe_list(:, i)
                 cycle
            end if

          end select


          select case(idx3)

          case(1)
            if((idx4 == 2) .or. (idx4 == 4)) then
                 if (main_rank) write(MsgOut,*) molecule%dihe_list(:, i)
                 cycle
            end if

          case(2)
            if((idx4 == 1) .or. (idx4 == 3)) then
                 if (main_rank) write(MsgOut,*) molecule%dihe_list(:, i)
                 cycle
            end if

          case(3)
            if(idx4 == 2) then
                 if (main_rank) write(MsgOut,*) molecule%dihe_list(:, i)
                 cycle
            end if

          case(4)
            if(idx4 == 1) then
                 if (main_rank) write(MsgOut,*) molecule%dihe_list(:, i)
                 cycle
            end if

          end select

        end do
      end if

      if(nerr_impr /= 0) then
        if (main_rank) write(MsgOut,*) "Inappropriate improper lists>"
        do i = 1, nimpr
          idx1 = int(molecule%fepgrp(molecule%impr_list(1, i)))
          idx2 = int(molecule%fepgrp(molecule%impr_list(2, i)))
          idx3 = int(molecule%fepgrp(molecule%impr_list(3, i)))
          idx4 = int(molecule%fepgrp(molecule%impr_list(4, i)))


          select case(idx1)

          case(1)
            if((idx2 == 2) .or. (idx2 == 4) .or. &
               (idx3 == 2) .or. (idx3 == 4) .or. &
               (idx4 == 2) .or. (idx4 == 4)) then
                 if (main_rank) write(MsgOut,*) molecule%impr_list(:, i)
                 cycle
            end if

          case(2)
            if((idx2 == 1) .or. (idx2 == 3) .or. &
               (idx3 == 1) .or. (idx3 == 3) .or. &
               (idx4 == 1) .or. (idx4 == 3)) then
                 if (main_rank) write(MsgOut,*) molecule%impr_list(:, i)
                 cycle
            end if

          case(3)
            if((idx2 == 2) .or. (idx3 == 2) .or. (idx4 == 2)) then
                 if (main_rank) write(MsgOut,*) molecule%impr_list(:, i)
                 cycle
            end if

          case(4)
            if((idx2 == 1) .or. (idx3 == 1) .or. (idx4 == 1)) then
                 if (main_rank) write(MsgOut,*) molecule%impr_list(:, i)
                 cycle
            end if

          end select


          select case(idx2)

          case(1)
            if((idx3 == 2) .or. (idx3 == 4) .or. &
               (idx4 == 2) .or. (idx4 == 4)) then
                 if (main_rank) write(MsgOut,*) molecule%impr_list(:, i)
                 cycle
            end if

          case(2)
            if((idx3 == 1) .or. (idx3 == 3) .or. &
               (idx4 == 1) .or. (idx4 == 3)) then
                 if (main_rank) write(MsgOut,*) molecule%impr_list(:, i)
                 cycle
            end if

          case(3)
            if((idx3 == 2) .or. (idx4 == 2)) then
                 if (main_rank) write(MsgOut,*) molecule%impr_list(:, i)
                 cycle
            end if

          case(4)
            if((idx3 == 1) .or. (idx4 == 1)) then
                 if (main_rank) write(MsgOut,*) molecule%impr_list(:, i)
                 cycle
            end if

          end select


          select case(idx3)

          case(1)
            if((idx4 == 2) .or. (idx4 == 4)) then
                 if (main_rank) write(MsgOut,*) molecule%impr_list(:, i)
                 cycle
            end if

          case(2)
            if((idx4 == 1) .or. (idx4 == 3)) then
                 if (main_rank) write(MsgOut,*) molecule%impr_list(:, i)
                 cycle
            end if

          case(3)
            if(idx4 == 2) then
                 if (main_rank) write(MsgOut,*) molecule%impr_list(:, i)
                 cycle
            end if

          case(4)
            if(idx4 == 1) then
                 if (main_rank) write(MsgOut,*) molecule%impr_list(:, i)
                 cycle
            end if

          end select

        end do
      end if

      if (main_rank) call error_msg(' ')
    end if

    return

  end subroutine check_parameter

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    correspondence_single_fep
  !> @brief        make correspondence list of atoms with single toplogy
  !! @authors      HO
  !! @param[inout] molecule : structure of molecule
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine correspondence_single_fep(molecule)

    ! formal arguments
    type(s_molecule), intent(inout) :: molecule

    ! local variables
    integer                         :: iA, iB, k

    allocate(molecule%id_singleA(molecule%num_atoms_fep(1)))
    allocate(molecule%id_singleB(molecule%num_atoms_fep(2)))

    iA = 1
    iB = 1
    do k = 1, molecule%num_atoms
      if (molecule%fepgrp(k) == 1) then
        molecule%id_singleA(iA) = k
        iA = iA + 1
      else if (molecule%fepgrp(k) == 2) then
        molecule%id_singleB(iB) = k
        iB = iB + 1
      end if
    end do

    return

  end subroutine correspondence_single_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    count_hbonds_single_fep
  !> @brief        count hydrogen atoms in single topology in FEP
  !!               This subroutine should be called after check_light_atom_name
  !!               is called.
  !! @authors      HO
  !! @param[inout] molecule : structure of molecule
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine count_hbonds_single_fep(molecule)

    ! formal arguments
    type(s_molecule),         intent(inout) :: molecule

    ! local variables
    integer                  :: k

    molecule%num_hbonds_singleA = 0
    molecule%num_hbonds_singleB = 0

    do k = 1, molecule%num_atoms
      if (molecule%fepgrp(k) == 1) then
        if (molecule%light_atom_name(k)) then
          molecule%num_hbonds_singleA = molecule%num_hbonds_singleA + 1
        end if
      else if (molecule%fepgrp(k) == 2) then
        if (molecule%light_atom_name(k)) then
          molecule%num_hbonds_singleB = molecule%num_hbonds_singleB + 1
        end if
      end if
    end do

    return

  end subroutine count_hbonds_single_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_bond_group_fep
  !> @brief        make fep group for bond, angle, dihedral, and cmap
  !! @authors      HO
  !! @param[inout] molecule : structure of molecule
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine make_bond_group_fep(molecule)

    ! formal arguments
    type(s_molecule), target, intent(inout) :: molecule

    ! local variables
    integer                  :: i1, i2, i3, i4
    integer                  :: i5, i6, i7, i8
    integer                  :: idx
    integer,         pointer :: fepgrp_bond(:,:)
    integer,         pointer :: fepgrp_angl(:,:,:)
    integer,         pointer :: fepgrp_dihe(:,:,:,:)
    integer,         pointer :: fepgrp_cmap(:)

    fepgrp_bond => molecule%fepgrp_bond
    fepgrp_angl => molecule%fepgrp_angl
    fepgrp_dihe => molecule%fepgrp_dihe
    fepgrp_cmap => molecule%fepgrp_cmap

    fepgrp_bond(:,:) = 0
    fepgrp_angl(:,:,:) = 0
    fepgrp_dihe(:,:,:,:) = 0
    fepgrp_cmap(:) = 0

    ! Bond table
    !
    do i2 = 1, 5
      do i1 = 1, 5
        ! preserve: all 5
        if ((i1 == 5) .and. (i2 == 5)) then
          fepgrp_bond(i1,i2) = 5
        end if
        ! singleA: all 1, 1 and 5
        if (((i1 == 1) .or. (i1 == 5)) .and. ((i2 == 1) .or. (i2 == 5))) then
          if ((i1 == 1) .or. (i2 == 1)) then
            fepgrp_bond(i1,i2) = 1
          end if
        end if
        ! singleB: all 2, 2 and 5
        if (((i1 == 2) .or. (i1 == 5)) .and. ((i2 == 2) .or. (i2 == 5))) then
          if ((i1 == 2) .or. (i2 == 2)) then
            fepgrp_bond(i1,i2) = 2
          end if
        end if
        ! dualA: all 3, 3 and 5, 3 and 1
        if (((i1 == 3) .or. (i1 == 5) .or. (i1 == 1)) .and. &
          ((i2 == 3) .or. (i2 == 5) .or. (i2 == 1))) then
          if ((i1 == 3) .or. (i2 == 3)) then
            fepgrp_bond(i1,i2) = 3
          end if
        end if
        ! dualB: all 4, 4 and 5, 4 and 2
        if (((i1 == 4) .or. (i1 == 5) .or. (i1 == 2)) .and. &
          ((i2 == 4) .or. (i2 == 5) .or. (i2 == 2))) then
          if ((i1 == 4) .or. (i2 == 4)) then
            fepgrp_bond(i1,i2) = 4
          end if
        end if
      end do
    end do

    ! Angle table
    !
    do i3 = 1, 5
      do i2 = 1, 5
        do i1 = 1, 5
          ! preserve: all 5
          if ((i1 == 5) .and. (i2 == 5) .and. (i3 == 5)) then
            fepgrp_angl(i1,i2,i3) = 5
          end if
          ! singleA: all 1, 1 and 5
          if (((i1 == 1) .or. (i1 == 5)) .and. &
            ((i2 == 1) .or. (i2 == 5)) .and. &
            ((i3 == 1) .or. (i3 == 5))) then
            if ((i1 == 1) .or. (i2 == 1) .or. (i3 == 1)) then
              fepgrp_angl(i1,i2,i3) = 1
            end if
          end if
          ! singleB: all 2, 2 and 5
          if (((i1 == 2) .or. (i1 == 5)) .and. &
            ((i2 == 2) .or. (i2 == 5)) .and. &
            ((i3 == 2) .or. (i3 == 5))) then
            if ((i1 == 2) .or. (i2 == 2) .or. (i3 == 2)) then
              fepgrp_angl(i1,i2,i3) = 2
            end if
          end if
          ! dualA: all 3, 3 and 5, 3 and 1
          if (((i1 == 3) .or. (i1 == 5) .or. (i1 == 1)) .and. &
            ((i2 == 3) .or. (i2 == 5) .or. (i2 == 1)) .and. &
            ((i3 == 3) .or. (i3 == 5) .or. (i3 == 1))) then
            if ((i1 == 3) .or. (i2 == 3) .or. (i3 == 3)) then
              fepgrp_angl(i1,i2,i3) = 3
            end if
          end if
          ! dualB: all 4, 4 and 5, 4 and 2
          if (((i1 == 4) .or. (i1 == 5) .or. (i1 == 2)) .and. &
            ((i2 == 4) .or. (i2 == 5) .or. (i2 == 2)) .and. &
            ((i3 == 4) .or. (i3 == 5) .or. (i3 == 2))) then
            if ((i1 == 4) .or. (i2 == 4) .or. (i3 == 4)) then
              fepgrp_angl(i1,i2,i3) = 4
            end if
          end if
        end do
      end do
    end do

    ! Dihedral table
    !
    do i4 = 1, 5
      do i3 = 1, 5
        do i2 = 1, 5
          do i1 = 1, 5
            ! preserve: all 5
            if ((i1 == 5) .and. (i2 == 5) .and. (i3 == 5) .and. (i4 == 5)) then
              fepgrp_dihe(i1,i2,i3,i4) = 5
            end if
            ! singleA: all 1, 1 and 5
            if (((i1 == 1) .or. (i1 == 5)) .and. &
              ((i2 == 1) .or. (i2 == 5)) .and. &
              ((i3 == 1) .or. (i3 == 5)) .and. &
              ((i4 == 1) .or. (i4 == 5))) then
              if ((i1 == 1) .or. (i2 == 1) .or. (i3 == 1) .or. (i4 == 1)) then
                fepgrp_dihe(i1,i2,i3,i4) = 1
              end if
            end if
            ! singleB: all 2, 2 and 5
            if (((i1 == 2) .or. (i1 == 5)) .and. &
              ((i2 == 2) .or. (i2 == 5)) .and. &
              ((i3 == 2) .or. (i3 == 5)) .and. &
              ((i4 == 2) .or. (i4 == 5))) then
              if ((i1 == 2) .or. (i2 == 2) .or. (i3 == 2) .or. (i4 == 2)) then
                fepgrp_dihe(i1,i2,i3,i4) = 2
              end if
            end if
            ! dualA: all 3, 3 and 5, 3 and 1
            if (((i1 == 3) .or. (i1 == 5) .or. (i1 == 1)) .and. &
              ((i2 == 3) .or. (i2 == 5) .or. (i2 == 1)) .and. &
              ((i3 == 3) .or. (i3 == 5) .or. (i3 == 1)) .and. &
              ((i4 == 3) .or. (i4 == 5) .or. (i4 == 1))) then
              if ((i1 == 3) .or. (i2 == 3) .or. (i3 == 3) .or. (i4 == 3)) then
                fepgrp_dihe(i1,i2,i3,i4) = 3
              end if
            end if
            ! dualB: all 4, 4 and 5, 4 and 2
            if (((i1 == 4) .or. (i1 == 5) .or. (i1 == 2)) .and. &
              ((i2 == 4) .or. (i2 == 5) .or. (i2 == 2)) .and. &
              ((i3 == 4) .or. (i3 == 5) .or. (i3 == 2)) .and. &
              ((i4 == 4) .or. (i4 == 5) .or. (i4 == 2))) then
              if ((i1 == 4) .or. (i2 == 4) .or. (i3 == 4) .or. (i4 == 4)) then
                fepgrp_dihe(i1,i2,i3,i4) = 4
              end if
            end if
          end do
        end do
      end do
    end do

    ! CMAP table
    !
    do i8 = 1, 5
      do i7 = 1, 5
        do i6 = 1, 5
          do i5 = 1, 5
            do i4 = 1, 5
              do i3 = 1, 5
                do i2 = 1, 5
                  do i1 = 1, 5
                    idx = i1 + 5*(i2-1 + 5*(i3-1 + 5*(i4-1 + 5*(i5-1 + 5*(i6-1 + 5*(i7-1 + 5*(i8-1)))))))
                    ! preserve: all 5
                    if ((i1 == 5) .and. (i2 == 5) .and. (i3 == 5) .and. (i4 == 5) .and. &
                      (i5 == 5) .and. (i6 == 5) .and. (i7 == 5) .and. (i8 == 5)) then
                      fepgrp_cmap(idx) = 5
                    end if
                    ! singleA: all 1, 1 and 5
                    if (((i1 == 1) .or. (i1 == 5)) .and. &
                      ((i2 == 1) .or. (i2 == 5)) .and. &
                      ((i3 == 1) .or. (i3 == 5)) .and. &
                      ((i4 == 1) .or. (i4 == 5)) .and. &
                      ((i5 == 1) .or. (i5 == 5)) .and. &
                      ((i6 == 1) .or. (i6 == 5)) .and. &
                      ((i7 == 1) .or. (i7 == 5)) .and. &
                      ((i8 == 1) .or. (i8 == 5))) then
                      if ((i1 == 1) .or. (i2 == 1) .or. (i3 == 1) .or. (i4 == 1) .or. &
                        (i5 == 1) .or. (i6 == 1) .or. (i7 == 1) .or. (i8 == 1)) then
                        fepgrp_cmap(idx) = 1
                      end if
                    end if
                    ! singleB: all 2
                    if (((i1 == 2) .or. (i1 == 5)) .and. &
                      ((i2 == 2) .or. (i2 == 5)) .and. &
                      ((i3 == 2) .or. (i3 == 5)) .and. &
                      ((i4 == 2) .or. (i4 == 5)) .and. &
                      ((i5 == 2) .or. (i5 == 5)) .and. &
                      ((i6 == 2) .or. (i6 == 5)) .and. &
                      ((i7 == 2) .or. (i7 == 5)) .and. &
                      ((i8 == 2) .or. (i8 == 5))) then
                      if ((i1 == 2) .or. (i2 == 2) .or. (i3 == 2) .or. (i4 == 2) .or. &
                        (i5 == 2) .or. (i6 == 2) .or. (i7 == 2) .or. (i8 == 2)) then
                        fepgrp_cmap(idx) = 2
                      end if
                    end if
                    ! dualA: all 3, 3 and 5, 3 and 1
                    if (((i1 == 3) .or. (i1 == 5) .or. (i1 == 1)) .and. &
                      ((i2 == 3) .or. (i2 == 5) .or. (i2 == 1)) .and. &
                      ((i3 == 3) .or. (i3 == 5) .or. (i3 == 1)) .and. &
                      ((i4 == 3) .or. (i4 == 5) .or. (i4 == 1)) .and. &
                      ((i5 == 3) .or. (i5 == 5) .or. (i5 == 1)) .and. &
                      ((i6 == 3) .or. (i6 == 5) .or. (i6 == 1)) .and. &
                      ((i7 == 3) .or. (i7 == 5) .or. (i7 == 1)) .and. &
                      ((i8 == 3) .or. (i8 == 5) .or. (i8 == 1))) then
                      if ((i1 == 3) .or. (i2 == 3) .or. (i3 == 3) .or. (i4 == 3) .or. &
                        (i5 == 3) .or. (i6 == 3) .or. (i7 == 3) .or. (i8 == 3)) then
                        fepgrp_cmap(idx) = 3
                      end if
                    end if
                    ! dualB: all 4, 4 and 5, 4 and 2
                    if (((i1 == 4) .or. (i1 == 5) .or. (i1 == 2)) .and. &
                      ((i2 == 4) .or. (i2 == 5) .or. (i2 == 2)) .and. &
                      ((i3 == 4) .or. (i3 == 5) .or. (i3 == 2)) .and. &
                      ((i4 == 4) .or. (i4 == 5) .or. (i4 == 2)) .and. &
                      ((i5 == 4) .or. (i5 == 5) .or. (i5 == 2)) .and. &
                      ((i6 == 4) .or. (i6 == 5) .or. (i6 == 2)) .and. &
                      ((i7 == 4) .or. (i7 == 5) .or. (i7 == 2)) .and. &
                      ((i8 == 4) .or. (i8 == 5) .or. (i8 == 2))) then
                      if ((i1 == 4) .or. (i2 == 4) .or. (i3 == 4) .or. (i4 == 4) .or. &
                        (i5 == 4) .or. (i6 == 4) .or. (i7 == 4) .or. (i8 == 4)) then
                        fepgrp_cmap(idx) = 4
                      end if
                    end if
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do

    return

  end subroutine make_bond_group_fep

end module sp_fep_topology_mod
