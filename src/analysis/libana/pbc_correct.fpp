!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pbc_correct_mod
!> @brief   module for PBC-correction
!! @authors Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../config.h"
#endif

module pbc_correct_mod

  use trajectory_str_mod
  use molecules_str_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! parameters
  integer,      public, parameter :: PBCCModeNo       = 1
  integer,      public, parameter :: PBCCModeMolecule = 2

  character(*), public, parameter :: PBCCModeTypes (2) = (/'NO      ', &
                                                           'MOLECULE'/)

  ! subroutines
  public  :: setup_pbc_correct
  public  :: run_pbc_correct
  private :: wrap_molecules

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_pbc_correct
  !> @brief        setup pbc-correction
  !! @authors      NT
  !! @param[in]    pbcc_mode : PBC-correct mode
  !! @param[inout] molecule  : molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine setup_pbc_correct(pbcc_mode, molecule)

    ! formal argments
    integer,                 intent(in)    :: pbcc_mode
    type(s_molecule),        intent(inout) :: molecule

    ! local variables
    integer                  :: i, j, natom, nmole, first, last
    logical                  :: set_mass


    ! setup for PBC-correct mode "molecule"
    if (pbcc_mode /= PBCCModeMolecule) &
      return

    ! check molecules information
    if (size(molecule%molecule_atom_no) == 0) &
      call error_msg('Setup_Pbc_Correct> Molecule information is not defined.')

    set_mass = .true.

    natom = size(molecule%atom_coord(1,:))
    nmole = size(molecule%molecule_atom_no)

    do i = 1, natom
      if (molecule%mass(i) > 0) then
        set_mass = .false.
        exit
      end if
    end do

    if (set_mass) then

      write(MsgOut,'(A)') &
           'Setup_PBC_Correct> Setup molecule mass information.'

      molecule%mass(:) = 1.0_wp
    end if

    do i = 1, nmole

      first = molecule%molecule_atom_no(i)
      if (i /= nmole) then
        last = molecule%molecule_atom_no(i + 1) - 1
      else
        last = natom
      end if

      molecule%molecule_mass(i) = 0.0_wp

      do j = first, last
        molecule%molecule_mass(i) = molecule%molecule_mass(i) + &
             molecule%mass(j)
      end do
    end do


    return

  end subroutine setup_pbc_correct

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_pbc_correct
  !> @brief        run pbc correct
  !! @authors      NT
  !! @param[in]    pbcc_mode  : PBC-correct mode
  !! @param[in]    molecule   : molecule information
  !! @param[inout] trajectory : trajectory information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine run_pbc_correct(pbcc_mode, molecule, trajectory)

    ! formal argments
    integer,                 intent(in)    :: pbcc_mode
    type(s_molecule),        intent(in)    :: molecule
    type(s_trajectory),      intent(inout) :: trajectory


    select case(pbcc_mode)

    case (PBCCModeNo)
      return

    case (PBCCModeMolecule)

      call wrap_molecules(molecule, trajectory)

    end select

    return

  end subroutine run_pbc_correct

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    wrap_molecules
  !> @brief        wrap molecules
  !! @authors      NT
  !! @param[in]    molecule   : molecule information
  !! @param[inout] trajectory : trajectory information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine wrap_molecules(molecule, trajectory)

    ! formal arguments
    type(s_molecule),           intent(in)    :: molecule
    type(s_trajectory), target, intent(inout) :: trajectory

    ! local variables
    real(wp)                    :: com(3)
    real(wp)                    :: boxx, boxy, boxz
    real(wp)                    :: iboxx, iboxy, iboxz
    integer                     :: i, j, initial, final

    real(wp),           pointer :: coord(:,:)


    coord => trajectory%coord

    boxx  = trajectory%pbc_box(1,1)
    boxy  = trajectory%pbc_box(2,2)
    boxz  = trajectory%pbc_box(3,3)
    iboxx = 1.0_wp/boxx
    iboxy = 1.0_wp/boxy
    iboxz = 1.0_wp/boxz

    do i = 1, molecule%num_molecules

      initial  = molecule%molecule_atom_no(i)
      if (i /= molecule%num_molecules) then
        final = molecule%molecule_atom_no(i+1) - 1
      else
        final = molecule%num_atoms
      end if

      ! compute center of mass of a molecule
      !
      com(1) = 0.0_wp
      com(2) = 0.0_wp
      com(3) = 0.0_wp

      do j = initial, final
        com(1) = com(1) + coord(1,j) * molecule%mass(j)
        com(2) = com(2) + coord(2,j) * molecule%mass(j)
        com(3) = com(3) + coord(3,j) * molecule%mass(j)
      end do

      com(1) = com(1)/molecule%molecule_mass(i)
      com(2) = com(2)/molecule%molecule_mass(i)
      com(3) = com(3)/molecule%molecule_mass(i)

      ! move molecule into the unit cell
      !
      do j = initial, final
        coord(1,j) = coord(1,j) - boxx * nint(com(1)*iboxx)
        coord(2,j) = coord(2,j) - boxy * nint(com(2)*iboxy)
        coord(3,j) = coord(3,j) - boxz * nint(com(3)*iboxz)
      end do

    end do

    return

  end subroutine wrap_molecules
  
end module pbc_correct_mod
