!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   qa_cg_setup_mod
!> @brief   setup variables and structures in QVAL_ANALYSIS
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module qa_cg_setup_mod

  use qa_cg_control_mod
  use qa_cg_option_mod
  use qa_cg_option_str_mod
  use trajectory_mod
  use output_mod
  use input_mod
  use trajectory_str_mod
  use output_str_mod
  use input_str_mod
  use select_mod
  use molecules_mod
  use molecules_str_mod
  use fileio_prmtop_mod
  use fileio_ambcrd_mod
  use fileio_grotop_mod
  use fileio_grocrd_mod
  use fileio_psf_mod
  use fileio_pdb_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public  :: setup
  private :: setup_contact

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup
  !> @brief        setup variables and structures in QVAL_ANALYSIS
  !! @authors      NT
  !! @param[in]    ctrl_data  : information of control parameters
  !! @param[inout] molecule   : molecule information
  !! @param[inout] trj_list   : trajectory file list information
  !! @param[inout] output     : output information
  !! @param[inout] option     : option information
  !! @param[inout] contact    : contact information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup(ctrl_data, molecule, trj_list, trajectory, output, option, contact)

    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_molecule),        intent(inout) :: molecule
    type(s_trj_list),        intent(inout) :: trj_list
    type(s_trajectory),      intent(inout) :: trajectory
    type(s_output),          intent(inout) :: output
    type(s_option),          intent(inout) :: option
    type(s_contact),         intent(inout) :: contact

    ! local variables
    type(s_psf)              :: psf
    type(s_pdb)              :: ref
    type(s_prmtop)           :: prmtop
    type(s_ambcrd)           :: ambref
    type(s_grotop)           :: grotop
    type(s_grocrd)           :: grocrd


    ! input files
    !
    call input_files(ctrl_data%inp_info, &
        psf=psf,            &
        ref=ref,            &
        prmtop=prmtop,      &
        ambref=ambref,      &
        grotop=grotop,      &
        grocrd=grocrd)


    ! define molecules
    !
    call define_molecules(molecule, pdb    = ref,    &
        psf    = psf,    &
        prmtop = prmtop, &
        ambcrd = ambref, &
        grotop = grotop, &
        grocrd = grocrd)

    call setup_contact(grotop, contact)

    call dealloc_pdb_all(ref)
    call dealloc_psf_all(psf)
    call dealloc_prmtop_all(prmtop)
    call dealloc_ambcrd_all(ambref)
    call dealloc_grotop_all(grotop)
    call dealloc_grocrd_all(grocrd)


    ! setup trajectory
    !
    call setup_trajectory(ctrl_data%trj_info, molecule, trj_list, trajectory)

    ! setup selection
    !
    call setup_selection(ctrl_data%sel_info, molecule)

    ! setup option
    !
    call setup_option(ctrl_data%opt_info, ctrl_data%sel_info, &
        molecule, option)

    ! setup output
    !
    call setup_output(ctrl_data%out_info, output)

    return

  end subroutine setup


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_contact
  !> @brief        define contact in CG energy functions
  !! @authors      CT
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[inout] contact  : contact information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_contact(grotop, contact)

    implicit none
    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    type(s_contact),         intent(inout) :: contact

    ! local variables
    real(wp)                               :: distance
    integer                                :: i, j, k
    integer                                :: natom, npairs, ioffset

    type(s_grotop_mol), pointer            :: gromol


    ! ---------------------
    ! count number of pairs
    ! ---------------------
    npairs = 0
    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        npairs  = npairs + gromol%num_pairs
      end do
    end do

    contact%n_pair = npairs
    allocate(contact%r0_ij(npairs), contact%cnt_pair(2, npairs))

    ! -----------------------
    ! determine contact pairs
    ! -----------------------
    natom  = 0
    npairs = 0
    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        ioffset = natom
        natom   = natom + gromol%num_atoms
        do k = 1, gromol%num_pairs
          npairs = npairs + 1
          contact%cnt_pair(1,npairs) = gromol%pairs(k)%atom_idx1 + ioffset
          contact%cnt_pair(2,npairs) = gromol%pairs(k)%atom_idx2 + ioffset
          contact%r0_ij(npairs) = 10.0_wp * gromol%pairs(k)%v
        end do                ! k
      end do                  ! j
    end do                    ! i

    return

  end subroutine setup_contact


end module qa_cg_setup_mod
