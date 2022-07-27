!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   qg_option_str_mod
!> @brief   structure of option information
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module qg_option_str_mod

  use fileio_psf_mod
  use constants_mod
  use select_atoms_str_mod
  use string_mod
  use messages_mod

  implicit none
  private

  ! structures
  type, public :: s_option

    logical              :: check_only
    integer              :: coord_format
    character(MaxLine)   :: qmmm_atom_exp = ''
    type(s_selatoms)     :: qmmm_atom
    type(s_selatoms)     :: qmmm_atom_trj
    character(MaxLine)   :: qm_atom_exp   = ''
    type(s_selatoms)     :: qm_atom
    logical              :: reconcile_num_atoms 
    character(MaxLine)   :: origin_atom_exp  = ''
    type(s_selatoms)     :: origin_atom
    integer, allocatable :: frame_no(:)
    type(s_psf)          :: dup_psf

  end type s_option

  ! parameters
  integer,      public, parameter :: CrdFormatPDB     = 1
  integer,      public, parameter :: CrdFormatAmber   = 2
  integer,      public, parameter :: CrdFormatGromacs = 3
  integer,      public, parameter :: CrdFormatCharmm  = 4

  character(*), public, parameter :: CrdFormatTypes(4) = (/'PDB     ', &
                                                           'AMBER   ', &
                                                           'GROMACS ', &
                                                           'CHARMM  '/)

  ! subroutines
  public  :: dealloc_option
  public  :: alloc_option
  public  :: duplicate_psf

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_option
  !> @brief        deallocate option
  !! @authors      NT
  !! @param[inout] option : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine dealloc_option(option)

    ! formal arguments
    type(s_option),          intent(inout) :: option

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0


    call dealloc_selatoms(option%qmmm_atom)
    call dealloc_selatoms(option%qmmm_atom_trj)
    call dealloc_selatoms(option%qm_atom)
    call dealloc_selatoms(option%origin_atom)

    if (allocated(option%frame_no)) then
      deallocate(option%frame_no, stat = dealloc_stat)
    end if
    if (dealloc_stat /= 0) call error_msg_dealloc

    call dealloc_psf_all(option%dup_psf)

    return

  end subroutine dealloc_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_option
  !> @brief        allocate option
  !! @authors      NT
  !! @param[inout] option : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine alloc_option(option, nframes)

    ! formal arguments
    type(s_option),          intent(inout) :: option
    integer,                 intent(in)    :: nframes

    ! local variables
    integer                  :: alloc_stat, dealloc_stat

    alloc_stat   = 0
    dealloc_stat = 0

    if (allocated(option%frame_no)) then
      deallocate(option%frame_no, stat = dealloc_stat)
    end if
    allocate(option%frame_no(1: nframes), stat = alloc_stat)
    
    if (alloc_stat /= 0) call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    option%frame_no = 0
    return

  end subroutine alloc_option

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine duplicate_psf(org_psf, dup_psf)

    ! formal arguments 
    type(s_psf),             intent(in)    :: org_psf
    type(s_psf),             intent(out)   :: dup_psf

    ! local variables
    integer                  :: i
    integer                  :: natom, nbond, ntheta, nphi, nimphi, &
                                ndon, nacc, nnb, ngrp, nst2, molnt, &
                                numlp, numlph, ncrterm
    real(wp)                 :: total_charge


    natom   = org_psf%num_atoms
    nbond   = org_psf%num_bonds
    ntheta  = org_psf%num_angles
    nphi    = org_psf%num_dihedrals
    nimphi  = org_psf%num_impropers
    ngrp    = org_psf%num_groups
    ncrterm = org_psf%num_cross_terms
    ndon    = org_psf%num_HB_donors
    nacc    = org_psf%num_HB_acceptors
    nnb     = org_psf%num_NB_exclusions

    ! init pdb
    call init_psf(dup_psf)

    ! setup atom data
    call alloc_psf(dup_psf, PsfAtom, natom)
    call alloc_psf(dup_psf, PsfBond, nbond)
    call alloc_psf(dup_psf, PsfAngl, ntheta)
    call alloc_psf(dup_psf, PsfDihe, nphi)
    call alloc_psf(dup_psf, PsfImpr, nimphi)
    call alloc_psf(dup_psf, PsfDonr, ndon)
    call alloc_psf(dup_psf, PsfAcce, nacc)
    call alloc_psf(dup_psf, PsfNb, nnb)
    call alloc_psf(dup_psf, PsfGrp, ngrp)
    call alloc_psf(dup_psf, PsfCmap, ncrterm)

    dup_psf%type              = org_psf%type
    dup_psf%num_atoms         = natom
    dup_psf%num_bonds         = nbond
    dup_psf%num_angles        = ntheta
    dup_psf%num_dihedrals     = nphi
    dup_psf%num_impropers     = nimphi
    dup_psf%num_groups        = ngrp
    dup_psf%num_cross_terms   = ncrterm
    dup_psf%num_HB_donors     = ndon
    dup_psf%num_HB_acceptors  = nacc
    dup_psf%num_NB_exclusions = nnb
    dup_psf%total_charge      = org_psf%total_charge

    do i = 1, natom

      dup_psf%atom_no(i)        = i ! org_psf%atom_no(i)
      dup_psf%segment_name(i)   = org_psf%segment_name(i)
      dup_psf%residue_no(i)     = org_psf%residue_no(i) 
      dup_psf%residue_name(i)   = org_psf%residue_name(i)
      dup_psf%atom_name(i)      = org_psf%atom_name(i)
      dup_psf%atom_cls_name(i)  = org_psf%atom_cls_name(i)
      dup_psf%atom_cls_no(i)    = org_psf%atom_cls_no(i)
      dup_psf%charge(i)         = org_psf%charge(i)
      dup_psf%mass(i)           = org_psf%mass(i)
      dup_psf%imove(i)          = org_psf%imove(i)
      dup_psf%molecule_no(i)    = org_psf%molecule_no(i) 

    end do

    do i = 1, nbond
      dup_psf%bond_list(1: 2, i) = org_psf%bond_list(1: 2, i)
    end do

    do i = 1, ntheta
      dup_psf%angl_list(1: 3, i) = org_psf%angl_list(1: 3, i)
    end do

    do i = 1, nphi
      dup_psf%dihe_list(1: 4, i) = org_psf%dihe_list(1: 4, i)
    end do

    do i = 1, nimphi
      dup_psf%impr_list(1: 4, i) = org_psf%impr_list(1: 4, i)
    end do

    do i = 1, ndon
      dup_psf%donr_list(1: 2, i) = org_psf%donr_list(1: 2, i)
    end do

    do i = 1, nacc
      dup_psf%acce_list(1: 2, i) = org_psf%acce_list(1: 2, i)
    end do

    do i = 1, nnb
      dup_psf%nb_list(i) = org_psf%nb_list(i)
    end do

    do i = 1, ngrp
      dup_psf%grp_list(1: 3, i) = org_psf%grp_list(1: 3, i)
    end do

    do i = 1, ncrterm
      dup_psf%cmap_list(1: 8, i) = org_psf%cmap_list(1: 8, i)
    end do

    return
  end subroutine duplicate_psf

end module qg_option_str_mod
