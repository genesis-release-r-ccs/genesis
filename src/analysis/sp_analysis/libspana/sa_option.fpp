!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sa_option_mod
!> @brief   module for analysis options
!! @authors Isseki Yu (IY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module sa_option_mod

  use sa_option_str_mod
  use select_mod
  use select_atoms_mod
  use select_atoms_str_mod
  use molecules_str_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use constants_mod
  use mpi_parallel_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! structures
  type, public :: s_opt_info
    logical    :: wrap          = .false.
    real(wp)   :: buffer        = 0.0_wp
    integer    :: determine_box = DetermineBoxTrajectory
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
  !! @authors      IY
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine show_ctrl_option

    write(MsgOut,'(A)') '[SPANA_OPTION]'
    write(MsgOut,'(A)') '# wrap        = no         # wrap trajectories'
    write(MsgOut,'(A)') 'buffer        = 10'
    write(MsgOut,'(A)') 'box_size      = TRAJECTORY # (TRAJECTORY / MAX / MANUAL)'
    write(MsgOut,'(A)') ''

    return

  end subroutine show_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_option
  !> @brief        read control parameters in OPTION section
  !! @authors      IY
  !! @param[in]    handle   : unit number of control file
  !! @param[out]   opt_info : OPTION section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_option(handle, opt_info)

    ! parameters
    character(*),            parameter     :: Section = 'Spana_option'

    ! formal argments
    integer,                 intent(in)    :: handle
    type(s_opt_info),        intent(inout) :: opt_info


    ! read parameters
    !
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_logical(handle, Section, &
                              'wrap', opt_info%wrap)

    call read_ctrlfile_real(handle, Section, &
                              'buffer', opt_info%buffer)

    call read_ctrlfile_type(handle, Section, &
                              'box_size', opt_info%determine_box, &
                              DetermineBoxTypes)

    call end_ctrlfile_section(handle)


    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Ctrl_Option> Parameters of SPANA Options'

      if (opt_info%wrap) then
        write(MsgOut,'(A20,A3)') '  wrap molecule   = ', 'yes'
      else
        write(MsgOut,'(A20,A2)') '  wrap molecule   = ', 'no'
      end if

      write(MsgOut,'(A20,F8.3)') '  buffer distance = ', opt_info%buffer
      write(MsgOut,'(A20,A)')    '  box size        = ', &
           DetermineBoxTypes(opt_info%determine_box)

      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine read_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_option
  !> @brief        setup option information
  !! @authors      IY
  !! @param[in]    opt_info : OPTION section control parameters information
  !! @param[in]    sel_info : SELECTION section control parameters information
  !! @param[in]    molecule : molecule information
  !! @param[inout] option   : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_option(opt_info, sel_info, molecule, option)

    ! formal arguments
    type(s_opt_info),        intent(in)    :: opt_info
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_molecule),        intent(in)    :: molecule
    type(s_option),          intent(inout) :: option
    integer                                :: igroup, alloc_stat


    option%buffer         = opt_info%buffer
    option%wrap           = opt_info%wrap
    option%determine_box  = opt_info%determine_box

    alloc_stat = 0

    allocate(option%analysis_atom_exps(size(sel_info%groups)), stat = alloc_stat)
    if (alloc_stat /= 0)   call error_msg_alloc

    allocate(option%analysis_atoms(size(sel_info%groups)), stat = alloc_stat)
    if (alloc_stat /= 0)   call error_msg_alloc

    do igroup = 1, size(sel_info%groups)
      option%analysis_atom_exps(igroup) = sel_info%groups(igroup)
    end do


    do igroup = 1, size(sel_info%groups)
      call select_atom(molecule, &
                       option%analysis_atom_exps(igroup), &
                       option%analysis_atoms(igroup))
    end do

    ! if (main_rank) then
    !   do igroup = 1, size(sel_info%groups)
    !     write(MsgOut,'(A,A7,I3,A2,I0)') 'Setup_Option> analysis atom count :', &
    !                                     '(group ',igroup,') ',                 &
    !                                     size(option%analysis_atoms(igroup)%idx)
    !   end do
    !   write(MsgOut,'(A)') ' '
    ! end if

    return

  end subroutine setup_option

end module sa_option_mod
