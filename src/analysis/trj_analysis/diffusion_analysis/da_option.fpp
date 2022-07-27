!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   da_option_mod
!> @brief   module for analysis options
!! @authors Donatas Surblys (DS), Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module da_option_mod

  use da_option_str_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use constants_mod
  use output_str_mod
  use fileio_mod

  implicit none
  private

  ! structures
  type, public :: s_opt_info
    logical                         :: check_only    = .false.
    logical                         :: allow_backup  = .false.
    real(wp)                        :: time_step     =  1.0d0
    real(wp)                        :: distance_unit =  1.0d0
    character(MaxLine)              :: start         =  "20 %"
    character(MaxLine)              :: stop          = "100 %"
    character(MaxLine)              :: ndofs         =     "3"
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
  !! @authors      DS
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine show_ctrl_option

    write(MsgOut,'(A)') '[OPTION]'
    write(MsgOut,'(A)') 'check_only     = NO      # only checking input files (YES/NO)'
    write(MsgOut,'(A)') 'allow_backup   = NO      # backup existing output files (YES/NO)'
    write(MsgOut,'(A)') 'time_step      = 1.00    # time interval (ps) of one step in MSD data'
    write(MsgOut,'(A)') 'distance_unit  = 1.00    # distance unit of MSD data files in angstroms'
    write(MsgOut,'(A)') 'start          = 20 %    # data point from which to start fitting'
    write(MsgOut,'(A)') '                         # n:    start at n-th step'
    write(MsgOut,'(A)') '                         # n %:  start at n percent'
    write(MsgOut,'(A)') '                         # n ps: start at n picoseconds'
    write(MsgOut,'(A)') 'stop           = 100 %   # data point from which to stop fitting'
    write(MsgOut,'(A)') '                         # n:    stop at n-th step'
    write(MsgOut,'(A)') '                         # n %:  stop at n percent'
    write(MsgOut,'(A)') '                         # n ps: stop at n picoseconds'
    write(MsgOut,'(A)') 'ndofs          = 3       # degrees of freedom (i.e. number of axes used) for MSD'
    write(MsgOut,'(A)') '                         # when a single number, sets for all MSD data sets'
    write(MsgOut,'(A)') '                         # when space seperated numbers, sets for each'

    return

  end subroutine show_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_option
  !> @brief        read control parameters in OPTION section
  !! @authors      DS
  !! @param[in]    handle   : unit number of control file
  !! @param[out]   opt_info : OPTION section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_option(handle, opt_info)

    ! formal argments
    integer, intent(in)                :: handle
    type(s_opt_info),  intent(inout)   :: opt_info

    character(*), parameter            :: section = "Option"


    ! read parameters
    !
    call begin_ctrlfile_section(handle, section)

    call read_ctrlfile_logical(handle, Section, 'check_only',   &
                              opt_info%check_only)
    call read_ctrlfile_logical(handle, Section, 'allow_backup', &
                              opt_info%allow_backup)
    call read_ctrlfile_real  (handle, section,  'time_step',    &
                              opt_info%time_step)
    call read_ctrlfile_real  (handle, section,  'distance_unit',&
                              opt_info%distance_unit)
    call read_ctrlfile_string(handle, section,  'start',        &
                              opt_info%start)
    call read_ctrlfile_string(handle, section,  'stop',         &
                              opt_info%stop)
    call read_ctrlfile_string(handle, section,  'ndofs',        &
                              opt_info%ndofs)

    call end_ctrlfile_section(handle)


    ! write parameters to MsgOut
    !
    write(MsgOut,'(A)') 'Read_Ctrl_Option> Parameters of ' // section

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

    write(MsgOut,'(A20,f5.2)') '  time_step       = ', opt_info%time_step
    write(MsgOut,'(A20,f5.2)') '  distance_unit   = ', opt_info%distance_unit
    write(MsgOut,'(A20,A)')    '  start           = ', trim(opt_info%start)
    write(MsgOut,'(A20,A)')    '  stop            = ', trim(opt_info%stop)
    write(MsgOut,'(A20,A)')    '  ndofs           = ', trim(opt_info%ndofs)
    write(MsgOut,'(A)')        ''

    return

  end subroutine read_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_option
  !> @brief        setup option information
  !! @authors      DS
  !! @param[in]    opt_info : OPTION section control parameters information
  !! @param[inout] option   : option information
  !! @param[in]    output   : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_option(opt_info, option, output)

    ! formal arguments
    type(s_opt_info), intent(in)                 :: opt_info
    type(s_option),   intent(inout)              :: option
    type(s_output),   intent(in)                 :: output

    logical :: output_exists
    integer :: i, msd_data, strlen_start, n, strlen_stop

    ! check only
    !
    option%check_only   = opt_info%check_only

    ! backup output file
    !
    call setup_backup(opt_info%allow_backup)

    ! time step
    !
    option%time_step     = opt_info%time_step
    option%distance_unit = opt_info%distance_unit

    ! number of degrees of freedom
    !
    n = split_num(trim(opt_info%ndofs))
    allocate(option%ndofs(n))
    call split(n, n, opt_info%ndofs, option%ndofs)

    ! starting point
    !
    strlen_start = len_trim(opt_info%start)
    if (opt_info%start(max(1, strlen_start):strlen_start) .eq. "%") then
      read(opt_info%start(1:strlen_start-1), *) option%start_percent
    else if (opt_info%start(max(1, strlen_start-1):strlen_start) .eq. "ps") then
      read(opt_info%start(1:strlen_start-2), *) option%start_time
    else
      read(opt_info%start, *) option%start_step
    end if

    ! end point
    !
    strlen_stop = len_trim(opt_info%stop)
    if (opt_info%stop(max(1, strlen_stop):strlen_stop) .eq. "%" ) then
      read(opt_info%stop(1:strlen_stop-1), *) option%stop_percent
    else if (opt_info%stop(max(1, strlen_stop-1):strlen_stop) .eq. "ps") then
      read(opt_info%stop(1:strlen_stop-2), *) option%stop_time
    else
      read(opt_info%stop, *) option%stop_step
    end if


    return

  end subroutine setup_option

end module da_option_mod
