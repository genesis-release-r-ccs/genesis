!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ma_option_mod
!> @brief   module for analysis options
!! @authors Donatas Surblys (DS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ma_option_mod

  use ma_option_str_mod
  use select_mod
  use select_atoms_mod
  use molecules_str_mod
  use fileio_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use select_molecules_mod
  use constants_mod

  implicit none
  private

  character(*), parameter :: default_section = "OPTION"
  character(*), parameter :: default_axes    = "x y z"

  ! structures
  type, public :: s_opt_info
    logical                                       :: check_only   = .false.
    logical                                       :: allow_backup = .false.
    integer                                       :: delta        = 10
    logical                                       :: oversample   = .true.
    character(MaxLine), dimension(:), allocatable :: axes
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
  !! @param[in]    section  : name of section
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine show_ctrl_option

    write(MsgOut,'(A)') '[SELECTION]'
    write(MsgOut,'(A)') 'group1         = rnam:TIP3       # selection group 1'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '[MOLECULE_SELECTION]'
    write(MsgOut,'(A)') 'selection1     = 1               # group number'
    write(MsgOut,'(A)') 'mode1          = ALL             # [ALL,SET,ANY] molecule selection mode for "selection1"'
    write(MsgOut,'(A)') '                                 # ALL: select all molecules in "selection1"'
    write(MsgOut,'(A)') '                                 # SET: define "selection1" as a single molecule'
    write(MsgOut,'(A)') '                                 # ANY: select molecules containing at least one atom in "selection1"'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '[OPTION]'
    write(MsgOut,'(A)') 'check_only     = NO              # only checking input files (YES/NO)'
    write(MsgOut,'(A)') 'allow_backup   = NO              # backup existing output files (YES/NO)'
    write(MsgOut,'(A)') 'delta          = 10              # span of MSD in steps'
    write(MsgOut,'(A)') 'oversample     = yes             # smoother graph'
    write(MsgOut,'(A)') 'axes1          = x y z           # axes for "selection1"'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '# Note1: Molecules in the trajectories should be unwrapped.'
    write(MsgOut,'(A)') '# Note2: OpenMP parallelization is availabe in this program.'
    write(MsgOut,'(A)') '#        For example, if 2 CPU cores are used, please set'
    write(MsgOut,'(A)') '#        $ export OMP_NUM_THREADS=2'

    return

  end subroutine show_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_option
  !> @brief        read control parameters in OPTION section
  !! @authors      DS
  !! @param[in]    handle   : unit number of control file
  !! @param[out]   opt_info : OPTION section control parameters information
  !! @param[in]    section  : name of section
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_option(handle, opt_info, molsel_info, section)

    ! formal argments
    integer,             intent(in)             :: handle
    type(s_molsel_info), intent(in)             :: molsel_info
    type(s_opt_info),    intent(inout)          :: opt_info
    character(*),        intent(in),   optional :: section

    character(MaxLine) :: section_name
    character(MaxLine) :: field_name
    integer            :: i

    if (present(section)) then
      section_name = section
    else
      section_name = default_section
    end if

    ! read parameters
    !
    call begin_ctrlfile_section(handle, trim(section_name))

    call read_ctrlfile_logical(handle, trim(section_name), &
                              'check_only', opt_info%check_only)

    call read_ctrlfile_logical(handle, trim(section_name), &
                              'allow_backup', opt_info%allow_backup)

    call read_ctrlfile_integer(handle, trim(section_name), &
                              'delta', opt_info%delta)

    call read_ctrlfile_logical(handle, trim(section_name), &
                              'oversample', opt_info%oversample)


    allocate(opt_info%axes(size(molsel_info%selections)))
    do i = 1, size(opt_info%axes)
      write(field_name, '("axes", i0)') i
      opt_info%axes(i) = default_axes
      call read_ctrlfile_string(handle, trim(section_name), trim(field_name), &
        opt_info%axes(i))
    end do

    call end_ctrlfile_section(handle)


    ! write parameters to MsgOut
    !
    write(MsgOut,'(A)') 'Read_Ctrl_Option> Parameters of ' // trim(section_name)

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

    write(MsgOut,'(A20,I0)')   '  delta           = ', opt_info%delta
    write(MsgOut,'(A20,L)')    '  oversample      = ', opt_info%oversample
    do i = 1, size(opt_info%axes)
      write(MsgOut,'(A16,i0,A,A)')  'axes', i, ' = ', trim(opt_info%axes(i))
    end do
    write(MsgOut,'(A)') ' '

    return

  end subroutine read_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_option
  !> @brief        setup option information
  !! @authors      DS
  !! @param[in]    opt_info : OPTION section control parameters information
  !! @param[in]    sel_info : SELECTION section control parameters information
  !! @param[in]    molecule : molecule information
  !! @param[inout] option   : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_option(opt_info, selmols, allmols, molecule, option)

    ! formal arguments
    type(s_opt_info),                                intent(in)    :: opt_info
    type(s_selmols),      dimension(:), allocatable, intent(inout) :: selmols
    type(s_one_molecule), dimension(:), allocatable, intent(inout) :: allmols
    type(s_molecule),                                intent(in)    :: molecule
    type(s_option),                                  intent(inout) :: option

    ! local variables
    integer               :: i, j
    logical, dimension(3) :: do_axes


    ! check only
    !
    option%check_only   = opt_info%check_only

    ! backup output file
    !
    call setup_backup(opt_info%allow_backup)

    ! selection
    !
    if (size(selmols) < 1) &
      call error_msg('Setup_Option> There must be at least one selection')

    call move_alloc(selmols, option%analysis_mols)
    call move_alloc(allmols, option%all_mols)

    write(MsgOut, '(A, i0)') &
      "Setup_Option> Total number of molecules in system: ", size(allmols)

    write(MsgOut,'(A,i0)') 'Setup_Option> Analysis molecule set count: ', &
         size(option%analysis_mols)

    do i = 1, size(option%analysis_mols)
      write(MsgOut,'(A,i0,A,i0)') &
        'Setup_Option> Analysis molecule count in set ', i, ': ', &
           size(option%analysis_mols(i)%mols)
    end do


    allocate(option%axes(size(opt_info%axes)))
    do i = 1, size(option%axes)
      do_axes = .false.

      do j = 1, len_trim(opt_info%axes(i))

        select case(opt_info%axes(i)(j:j))

        case("x","X")
          do_axes(1) = .true.

        case("y","Y")
          do_axes(2) = .true.

        case("z","Z")
          do_axes(3) = .true.

        case(" ")
          cycle

        case default
          call error_msg("Setup_Option> Illegal axis value")

        end select

      end do

      allocate(option%axes(i)%i(count(do_axes)))
      option%axes(i)%i(:) = pack([(j,j=1,3)], do_axes)

    end do

    write(MsgOut,'("")')

    option%delta = opt_info%delta
    option%oversample = opt_info%oversample

    return

  end subroutine setup_option

end module ma_option_mod
