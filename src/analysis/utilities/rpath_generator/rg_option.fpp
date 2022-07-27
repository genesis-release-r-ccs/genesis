!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   rg_option_mod
!> @brief   module for analysis options
!! @authors Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module rg_option_mod

  use rg_option_str_mod
  use input_mod
  use input_str_mod
  use select_mod
  use select_atoms_mod
  use molecules_str_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_opt_info

    integer                  :: nreplica
    integer                  :: cv_atom
    integer                  :: iseed = 777
    integer                  :: iter_reparam = 1

  end type s_opt_info

  ! subroutines
  public  :: show_ctrl_option
  public  :: read_ctrl_option
  public  :: setup_option

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_option
  !> @brief        show control parameters in OPTION section
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_option

    write(MsgOut,'(A)') '[RPATH]'
    write(MsgOut,'(A)') 'nreplica       = 8'
    write(MsgOut,'(A)') 'cv_atom        = 1'
    write(MsgOut,'(A)') 'iseed          = 777'
    write(MsgOut,'(A)') 'iter_reparam   = 1'
    write(MsgOut,'(A)') ' '

    return

  end subroutine show_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_option
  !> @brief        read control parameters in OPTION section
  !! @authors      NT
  !! @param[in]    handle   : unit number of control file
  !! @param[out]   opt_info : OPTION section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine read_ctrl_option(handle, opt_info)
  
    ! parameters
    character(*),            parameter     :: Section1 = 'RPATH'

    ! formal argments
    integer,                 intent(in)    :: handle
    type(s_opt_info),        intent(inout) :: opt_info

    ! local variables
    character(100)           :: line


    ! read parameters
    !

    ! 'RPATH' section
    call begin_ctrlfile_section(handle, Section1)

    call read_ctrlfile_integer(handle, Section1, 'nreplica', &
                               opt_info%nreplica)

    call read_ctrlfile_integer(handle, Section1, 'cv_atom', &
                               opt_info%cv_atom)

    call read_ctrlfile_integer(handle, Section1, 'iseed', &
                               opt_info%iseed)

    call read_ctrlfile_integer(handle, Section1, 'iter_reparam', &
                               opt_info%iter_reparam)

    call end_ctrlfile_section(handle)


    ! write parameters to MsgOut
    !
    write(MsgOut,'(A)') 'Read_Ctrl_Option> Parameters of Rpath'
    write(MsgOut,'(A20,I10)') &
         ' nreplica         = ', opt_info%nreplica
    write(MsgOut,'(A20,I10)') &
         ' cv_atom          = ', opt_info%cv_atom
    write(MsgOut,'(A20,I10)') &
         ' iseed            = ', opt_info%iseed
    write(MsgOut,'(A20,I10)') &
         ' iter_reparam     = ', opt_info%iter_reparam
    write(MsgOut,'(A)') ' '

    return

900 call error_msg('Read_Ctrl_Option> format error : "md_baro_moment"')

  end subroutine read_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_option
  !> @brief        setup option information
  !! @authors      NT
  !! @param[in]    opt_info : OPTION section control parameters information
  !! @param[in]    sel_info : SELECTION section control parameters information
  !! @param[in]    molecule : molecule information
  !! @param[inout] option   : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine setup_option(opt_info, sel_info, molecule, input, option)

    ! formal arguments
    type(s_opt_info),        intent(in)    :: opt_info
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_molecule),        intent(in)    :: molecule
    type(s_input),           intent(in)    :: input
    type(s_option),  target, intent(inout) :: option

    ! local variables
    integer                  :: igrp
    integer                  :: i

    write(MsgOut,'(a)') 'Setup_Option> :'
    write(MsgOut,'(a)') ''

    option%nreplica = opt_info%nreplica
    option%cv_atom  = opt_info%cv_atom
    option%iseed    = opt_info%iseed
    option%iter_reparam = opt_info%iter_reparam

    if (input%cvfile /= '') then
      option%is_cartesian = .false.
    else
      option%is_cartesian = .true.
    end if

    ! selection settings
    !

    ! num_atoms
    option%num_atoms = molecule%num_atoms

    ! selatoms
    allocate(option%selatoms(1:size(sel_info%groups)))

    do i = 1, size(sel_info%groups)
      call select_atom(molecule, sel_info%groups(i), option%selatoms(i))

      write(MsgOut,'(a,i12,a,a)') &
           ' group: ', i, ' :', trim(sel_info%groups(i))
      write(MsgOut,'(a,i12)') &
           '    # of selected atom = ', size(option%selatoms(i)%idx)
    end do

    allocate(option%groups(size(sel_info%groups)))
    do i = 1, size(sel_info%groups)
      option%groups(i) = sel_info%groups(i)
    end do

    write(MsgOut,'(a)') ''


    return

  end subroutine setup_option

end module rg_option_mod
