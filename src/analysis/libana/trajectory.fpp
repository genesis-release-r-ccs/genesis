!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   trajectory_mod
!> @brief   module for trajectory files
!! @authors Norio Takase (NT), Donatas Surblys (DS)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../config.h"
#endif

module trajectory_mod

  use trajectory_str_mod
  use molecules_str_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use constants_mod
  use mpi_parallel_mod
  use fileio_trj_mod

  implicit none
  private

  ! structures
  type, public :: s_trj_info
    character(MaxFilename), allocatable :: trj_files(:)
    integer,                allocatable :: md_steps(:)
    integer,                allocatable :: mdout_periods(:)
    integer,                allocatable :: ana_periods(:)
    integer,                allocatable :: start_steps(:)
    integer                             :: trj_format = TrjFormatDCD
    integer                             :: trj_type   = TrjTypeCoorBox
    integer                             :: trj_natom  = 0
  end type s_trj_info

  ! subroutines
  public  :: show_ctrl_trajectory
  public  :: read_ctrl_trajectory
  public  :: setup_trajectory

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_trajectory
  !> @brief        show control parameters in TRAJECTORY section
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_trajectory(tool_name)

    ! formal arguments
    character(len=*), optional, intent(in) :: tool_name

    if (.not.present(tool_name)) then
      write(MsgOut,'(A)') '[TRAJECTORY]'
      write(MsgOut,'(A)') '# trjfile1       = sample.dcd      # trajectory file'
      write(MsgOut,'(A)') '# md_step1       = 0               # number of MD steps'
      write(MsgOut,'(A)') '# mdout_period1  = 0               # MD output period'
      write(MsgOut,'(A)') '# ana_period1    = 1               # analysis period'
      write(MsgOut,'(A)') '# repeat1        = 1'
      write(MsgOut,'(A)') 'trj_format     = DCD             # (PDB/DCD)'
      write(MsgOut,'(A)') 'trj_type       = COOR+BOX        # (COOR/COOR+BOX)'
      write(MsgOut,'(A)') 'trj_natom      = 0               # (0:uses reference PDB atom count)'
      write(MsgOut,'(A)') ' '

    else
      select case (tool_name)

      case ('qg')
        write(MsgOut,'(A)') '[TRAJECTORY]'
        write(MsgOut,'(A)') '# trjfile1       = sample1.dcd      # trajectory file'
        write(MsgOut,'(A)') '# trjfile2       = sample2.dcd      # trajectory file'
        write(MsgOut,'(A)') '# trjfile3       = sample3.dcd      # trajectory file'
        write(MsgOut,'(A)') '# md_step1       = 4000000          # number of MD steps (set 1)'
        write(MsgOut,'(A)') '# mdout_period1  =    5000          # MD output period   (set 1)'
        write(MsgOut,'(A)') '# ana_period1    =    5000          # analysis period    (set 1)'
        write(MsgOut,'(A)') '# repeat1        = 3                # repeat count of set 1'
        write(MsgOut,'(A)') '# trjfile4       = sample4.dcd      # trajectory file'
        write(MsgOut,'(A)') '# md_step2       = 2000000          # number of MD steps (set 2)'
        write(MsgOut,'(A)') '# mdout_period2  =    5000          # MD output period   (set 2)'
        write(MsgOut,'(A)') '# ana_period2    =    5000          # analysis period    (set 2)'
        write(MsgOut,'(A)') '# repeat2        = 1                # repeat count of set 2'
        !write(MsgOut,'(A)') '# trjfile5       = sample5.dcd      # trajectory file'
        !write(MsgOut,'(A)') '# md_step3       = 4000000          # number of MD steps (set 3)'
        !write(MsgOut,'(A)') '# mdout_period3  =    5000          # MD output period   (set 3)'
        !write(MsgOut,'(A)') '# ana_period3    =    5000          # analysis period    (set 3)'
        !write(MsgOut,'(A)') '# trjfile6       = sample6.dcd      # trajectory file'
        !write(MsgOut,'(A)') '# repeat3        = 2                # repeat count of set 3'
        !write(MsgOut,'(A)') 'trj_format     = DCD             # (PDB/DCD)'
        write(MsgOut,'(A)') 'trj_type       = COOR+BOX        # (COOR/COOR+BOX)'
        write(MsgOut,'(A)') 'trj_natom      = 0               # (0:uses reference PDB atom count)'
        write(MsgOut,'(A)') '# Please see https://www.r-ccs.riken.jp/labs/cbrt/tutorial/faq/faq2/'
        write(MsgOut,'(A)') ' '

      case default
        write(MsgOut,'(A)') '[TRAJECTORY]'
        write(MsgOut,'(A)') '# trjfile1       = sample.dcd      # trajectory file'
        write(MsgOut,'(A)') '# md_step1       = 0               # number of MD steps'
        write(MsgOut,'(A)') '# mdout_period1  = 0               # MD output period'
        write(MsgOut,'(A)') '# ana_period1    = 1               # analysis period'
        write(MsgOut,'(A)') '# repeat1        = 1'
        write(MsgOut,'(A)') 'trj_format     = DCD             # (PDB/DCD)'
        write(MsgOut,'(A)') 'trj_type       = COOR+BOX        # (COOR/COOR+BOX)'
        write(MsgOut,'(A)') 'trj_natom      = 0               # (0:uses reference PDB atom count)'
        write(MsgOut,'(A)') ' '

      end select
    end if

    return

  end subroutine show_ctrl_trajectory

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_trajectory
  !> @brief        read TRAJECTORY section in the control file
  !! @authors      NT, DS
  !! @param[in]    handle   : unit number of control file
  !! @param[inout] trj_info : TRAJECTORY section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine read_ctrl_trajectory(handle, trj_info)
  
    ! parameters
    character(*),            parameter :: Section = 'Trajectory'

    ! formal argments
    integer,                 intent(in)    :: handle
    type(s_trj_info),        intent(inout) :: trj_info

    integer                  :: i, j, k, ntrj
    integer                  :: md_step, mdout_period, ana_period, repeat
    integer                  :: start_step
    character(MaxFilename)   :: value
    character(20)            :: trjname
    character(10)            :: numstr


    ! read parameters 
    !

    call begin_ctrlfile_section(handle, Section)

    ! check readable trajectory count
    ntrj = 0
    do while (.true.)

      value = ''
      write(trjname,'(a7,i0)') 'trjfile', ntrj + 1
      call read_ctrlfile_string(handle, Section, trjname, value)
      if (value == '') &
        exit

      ntrj = ntrj + 1

    end do

    ! allocate for trajectory file infos.
    allocate(trj_info%trj_files    (ntrj), &
             trj_info%md_steps     (ntrj), &
             trj_info%mdout_periods(ntrj), &
             trj_info%ana_periods  (ntrj), &
             trj_info%start_steps  (ntrj))
 
    ! read trajectory file name
    do i = 1, ntrj
      write(trjname,'(a7,i0)') "trjfile", i
      call read_ctrlfile_string(handle, Section, trjname, &
                                trj_info%trj_files(i))
    end do


    i = 1
    j = 1
    do while(ntrj > 0)
      
      write(numstr, '(I0)') j
      j = j + 1

      ! set default to 0 and read md steps
      md_step = 0
      call read_ctrlfile_integer(handle, Section, &
                                 'md_step'//numstr, md_step)

      ! set default to 1 and read mdout period
      mdout_period = 1
      call read_ctrlfile_integer(handle, Section, &
                                 'mdout_period'//numstr, mdout_period)
      if (mdout_period  <=  0) &
        call error_msg('Read_Ctrl_Trajectory> # of out-period is zero.')

      ! md_step <= 0 triggers automatic trajectory step detection,
      ! so mdout_period must be 1
      if (md_step <= 0 .and. mdout_period /= 1) &
        call error_msg('Read_Ctrl_Trajectory> # of out-period must be one &
        &when md-steps is less than one.')

      ana_period = 1
      repeat     = 1
      start_step = 1

      ! read analysis period
      call read_ctrlfile_integer(handle, Section, &
                                 'ana_period'//numstr, ana_period)

      ! read start steps
      call read_ctrlfile_integer(handle, Section, &
                                 'start_step'//numstr, start_step)

      ! read repeats
      call read_ctrlfile_integer(handle, Section, &
                                 'repeat'//numstr, repeat)

      md_step    = md_step    / mdout_period
      ana_period = ana_period / mdout_period

      if (ana_period <= 0)   ana_period   = 1
      if (repeat <= 0)       repeat       = 1

      do k = 1, repeat
        if (i > ntrj) &
          call error_msg( &
            'Read_Ctrl_Trajectory> trj-files and md-steps is mismatch.')

        trj_info%md_steps(i)      = md_step
        trj_info%mdout_periods(i) = mdout_period
        trj_info%ana_periods(i)   = ana_period
        trj_info%start_steps(i)   = start_step
        i = i + 1

      end do

      if (i > ntrj) &
        exit

    end do

    ! read trajectory format 
    call read_ctrlfile_type  (handle, Section, &
                              'trj_format', trj_info%trj_format, TrjFormatTypes)

    ! read trajectory type
    call read_ctrlfile_type  (handle, Section, &
                              'trj_type', trj_info%trj_type, TrjTypeTypes)

    ! read trajectory # of atoms
    call read_ctrlfile_integer(handle, Section, &
                              'trj_natom', trj_info%trj_natom)

    call end_ctrlfile_section(handle)

    do i = 1, ntrj

      ! automatically detect number of trajectory steps when md_steps(i) < 1
      if (trj_info%md_steps(i) < 1) then

        md_step = get_num_steps_trj(trj_info%trj_files(i), &
          trj_info%trj_format, trj_info%trj_type)

        ! when md_steps(i) < 0, reduce step number by their value
        trj_info%md_steps(i) = md_step + trj_info%md_steps(i)

      end if

    end do

    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Ctrl_Trajectory> Parameters of Trajectory'
      write(MsgOut,'(A20,I10)') '  # of trj files  = ', ntrj
      do i = 1, ntrj
        write(MsgOut,'(A20,I0,A3,A)') &
             '    trajectory file ', i, ' = ', trim(trj_info%trj_files(i))
        write(MsgOut,'(A20,I10)') &
         '         md steps : ', trj_info%md_steps(i)*trj_info%mdout_periods(i)
        write(MsgOut,'(A20,I10)') &
         '     mdout period : ', trj_info%mdout_periods(i)
        write(MsgOut,'(A20,I10)') &
         '       ana period : ', trj_info%ana_periods(i)*trj_info%mdout_periods(i)
        write(MsgOut,'(A20,I10)') &
         '     start step   : ', trj_info%start_steps(i)
      end do

      write(MsgOut,'(A20,A10)') &
           '  trj format      = ', TrjFormatTypes(trj_info%trj_format)
      write(MsgOut,'(A20,A10)') &
           '  trj type        = ', TrjTypeTypes(trj_info%trj_type)
      if (trj_info%trj_natom /= 0) &
        write(MsgOut,'(A20,I10)')  &
           '  # of trj atom   = ', trj_info%trj_natom

      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine read_ctrl_trajectory

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_trajectory
  !> @brief        setup trajectory information
  !! @authors      NT
  !! @param[in]    trj_info   :TRAJECTORY section control parameters information
  !! @param[in]    molecule   : molecule information
  !! @param[inout] trj_list   : trajectory list information
  !! @param[inout] trajectory : trajectory information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine setup_trajectory(trj_info, &
                              molecule, &
                              trj_list, &
                              trajectory)
  
    ! formal argments
    type(s_trj_info),        intent(in)    :: trj_info
    type(s_molecule),        intent(in)    :: molecule
    type(s_trj_list),        intent(inout) :: trj_list
    type(s_trajectory),      intent(inout) :: trajectory

    ! local variables
    integer                  :: ntrj, natom


    ! setup trajectory list
    !
    ntrj = size(trj_info%trj_files(:))

    call alloc_trj_list(trj_list, ntrj)

    trj_list%filenames(:)     = trj_info%trj_files(:)
    trj_list%md_steps(:)      = trj_info%md_steps(:)
    trj_list%mdout_periods(:) = trj_info%mdout_periods(:)
    trj_list%ana_periods(:)   = trj_info%ana_periods(:)
    trj_list%start_steps(:)   = trj_info%start_steps(:)
    trj_list%trj_format       = trj_info%trj_format
    trj_list%trj_type         = trj_info%trj_type

    ! setup trajectory
    !
    natom = trj_info%trj_natom
    if (natom == 0) then
      natom = molecule%num_atoms
    end if

    call alloc_trajectory(trajectory, natom)

    trajectory%coord(:,:)   = 0.0_wp
    trajectory%pbc_box(:,:) = 0.0_wp

    return

  end subroutine setup_trajectory

end module trajectory_mod
