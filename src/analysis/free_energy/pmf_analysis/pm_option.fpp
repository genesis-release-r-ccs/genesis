!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pm_option_mod
!> @brief   module for analysis options
!! @authors Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module pm_option_mod

  use pm_option_str_mod
  use fileio_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_opt_info
    logical                         :: check_only      = .false.
    logical                         :: allow_backup    = .false.
    integer                         :: nreplicas       = 1
    integer                         :: dimension       = 1
    real(wp)                        :: temperature     = 300.0_wp
    real(wp)                        :: cutoff          = 0.0_wp
    character(MaxLine), allocatable :: grids(:)
    real(wp),           allocatable :: band_width(:)
    logical,            allocatable :: is_periodic(:)
    real(wp),           allocatable :: box_size(:)
    integer                         :: output_type      = OutputTypeGNUPLOT
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

    write(MsgOut,'(A)') '[OPTION]'
    write(MsgOut,'(A)') 'check_only     = NO              # only checking input files (YES/NO)'
    write(MsgOut,'(A)') 'allow_backup   = NO              # backup existing output files (YES/NO)'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'nreplica       = 8'
    write(MsgOut,'(A)') 'dimension      = 1'
    write(MsgOut,'(A)') 'temperature    = 300'
    write(MsgOut,'(A)') 'cutoff         = 100             # cutoff distance measured by pathcv_analysis'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'grids1         = 1.0 10.0 100    # (min max num_of_bins)'
    write(MsgOut,'(A)') 'band_width1    = 0.1             # sigma of Gaussian kernel'
    write(MsgOut,'(A)') 'is_periodic1   = YES'
    write(MsgOut,'(A)') 'box_size1      = 360.0'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'output_type    = GNUPLOT  # GNUPLOT MATLAB'
    write(MsgOut,'(A)') ''

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
    character(*),            parameter     :: Section = 'OPTION'

    ! formal argments
    integer,                 intent(in)    :: handle
    type(s_opt_info),        intent(inout) :: opt_info

    integer                                :: i
    character(30)                          :: cdim


    ! read control parameters
    !

    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_logical(handle, Section, &
                               'check_only', opt_info%check_only)

    call read_ctrlfile_logical(handle, Section, &
                               'allow_backup', opt_info%allow_backup)

    call read_ctrlfile_integer(handle, Section, &
                               'nreplica', opt_info%nreplicas)

    call read_ctrlfile_integer(handle, Section, &
                               'dimension', opt_info%dimension)

    call read_ctrlfile_real   (handle, Section, &
                               'temperature', opt_info%temperature)

    call read_ctrlfile_real   (handle, Section, &
                               'cutoff', opt_info%cutoff)

    call read_ctrlfile_type   (handle, Section, &
                              'output_type', opt_info%output_type, OutputType)

    ! call read_ctrlfile_string (handle, Section, &
    !                            'center', opt_info%center)

    allocate(opt_info%grids(opt_info%dimension), &
             opt_info%band_width(opt_info%dimension), &
             opt_info%is_periodic(opt_info%dimension), &
             opt_info%box_size(opt_info%dimension))

    opt_info%is_periodic(1:opt_info%dimension) = .false.
    opt_info%box_size(1:opt_info%dimension)    = 0.0_wp
    opt_info%band_width(1:opt_info%dimension)  = 0.0_wp

    do i = 1, opt_info%dimension
      write(cdim,'(i0)') i
      call read_ctrlfile_string (handle, Section, &
                                 'grids'//cdim, opt_info%grids(i))

      write(cdim,'(i0)') i
      call read_ctrlfile_real (handle, Section, &
                                 'band_width'//cdim, opt_info%band_width(i))

      call read_ctrlfile_logical(handle, Section, &
                               'is_periodic'//cdim, opt_info%is_periodic(i))

      call read_ctrlfile_real(handle, Section, &
                               'box_size'//cdim, opt_info%box_size(i))

    end do

    call end_ctrlfile_section(handle)


    ! write parameter to MsgOut
    !

    write(MsgOut,'(A)') 'Read_Ctrl_Option> Parameters of OPTION'
    write(MsgOut,'(A)') ''

    if (opt_info%check_only) then
      write(MsgOut,'(A20,A3)') '  check only      = ', 'yes'
    else
      write(MsgOut,'(A20,A2)') '  check only      = ', 'no'
    end if

    if (opt_info%allow_backup) then
      write(MsgOut,'(A20,A3)') '  allow backup    = ', 'yes'
    else
      write(MsgOut,'(A20,A2)') '  allow backup    = ', 'no'
    end if

    write(MsgOut,'(A20,I10)') &
         '     dimension    = ', opt_info%dimension

    write(MsgOut,'(A20,I10)') &
         '     nreplica     = ', opt_info%nreplicas

    write(MsgOut,'(A20,F10.2)') &
         '     temperature  = ', opt_info%temperature

    write(MsgOut,'(A20,F10.2)') &
         '     cutoff       = ', opt_info%cutoff

    do i = 1, opt_info%dimension

      write(MsgOut,'(A20,I10)') &
        '   dimension :      ', i
      write(MsgOut,'(A20,A)') &
        '     grids        = ', trim(opt_info%grids(i))
      write(MsgOut,'(A20,F10.2)') &
        '     band_width   = ', opt_info%band_width(i)
      if (opt_info%is_periodic(i)) then
        write(MsgOut,'(A20,A3)') '     is_periodic  = ', 'yes'
        write(MsgOut,'(A20,F10.3)') &
          '     box_size     = ', opt_info%box_size(i)
      else
        write(MsgOut,'(A20,A2)') '     is_periodic  = ', 'no'
      end if
      write(MsgOut,'(A)') ''

    end do

    write(MsgOut,'(A)') ''

    return

  end subroutine read_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_option
  !> @brief        setup option information
  !! @authors      NT
  !! @param[in]    opt_info : OPTION section control parameters information
  !! @param[inout] option   : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine setup_option(opt_info, option)

    ! formal arguments
    type(s_opt_info),        intent(in)    :: opt_info
    type(s_option),          intent(inout) :: option

    ! local variables
    integer                  :: i, j, k
    real(wp)                 :: delta
    integer                  :: num_bins_max = 0
    character(100)           :: str
    real(wp),  allocatable   :: grid(:)


    option%check_only  = opt_info%check_only

    ! backup output file
    !
    call setup_backup(opt_info%allow_backup)

    option%nreplicas   = opt_info%nreplicas
    option%dimension   = opt_info%dimension
    option%temperature = opt_info%temperature
    option%cutoff      = opt_info%cutoff
    option%output_type = opt_info%output_type

    ! !
    allocate(option%grid_min (option%dimension), &
             option%grid_max (option%dimension), &
             option%num_grids(option%dimension), &
             option%band_width(option%dimension), &
             option%is_periodic(option%dimension), &
             option%box_size(option%dimension), &
             option%delta_grid(option%dimension))
    
    do i = 1, option%dimension
      str = opt_info%grids(i)
      do j = 1, len(str)
        if (str(j:j) == '(' .or. &
            str(j:j) == ')' .or. &
            str(j:j) == ',') &
            str(j:j) = ' '
      end do
      if (str /= '') then
        if (split_num(str) /= 3) &
           call error_msg('Setup_Option> grids must be (min max num_grids)')
        read(str,*) option%grid_min(i), &
                    option%grid_max(i), &
                    option%num_grids(i)
      else
        option%grid_min(i)  = 0.0_wp
        option%grid_max(i)  = 0.0_wp
        option%num_grids(i) = 0
      end if
    end do

    if (option%dimension == 1) then
      num_bins_max = option%num_grids(1) - 1
    else if (option%dimension == 2) then
      if (option%num_grids(1) > option%num_grids(2)) then
        num_bins_max = option%num_grids(1) - 1
      else
        num_bins_max = option%num_grids(2) - 1
      end if
    end if

    allocate(option%center(option%dimension, num_bins_max))

    do i = 1, option%dimension
      option%delta_grid(i) = (option%grid_max(i) - option%grid_min(i)) &
                             / (option%num_grids(i) - 1)   
      do j = 1, option%num_grids(i) - 1
        option%center(i, j) = option%grid_min(i) + 0.5_wp*option%delta_grid(i) + &
                               real((j-1),wp)*option%delta_grid(i)
      end do
    end do

    do i = 1, option%dimension
      write(MsgOut, '(a,i3)') 'Setup_Option> centers of grids in dimension ',i
      do j = 1, option%num_grids(i) - 1
        write(MsgOut,'(F10.4,$)') option%center(i, j)
      end do
      write(MsgOut, *)
      write(MsgOut, *)
    end do
    write(MsgOut, *)

    do i = 1, option%dimension
      option%delta_grid(i) = 0.5_wp * option%delta_grid(i)    
      option%band_width(i) = opt_info%band_width(i)
    end do

    option%is_periodic(1:option%dimension) = opt_info%is_periodic(1:option%dimension)
    option%box_size(1:option%dimension) = opt_info%box_size(1:option%dimension)

  end subroutine setup_option
  
end module pm_option_mod
