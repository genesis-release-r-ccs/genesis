!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_localres_mod
!> @brief   read local restraint file
!! @authors Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fileio_localres_mod

  use fileio_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_localres

    integer           :: num_funcs    = 0

    integer,          allocatable :: index_atoms(:,:)
    integer,          allocatable :: func(:)
    real(wp),         allocatable :: ref(:)
    real(wp),         allocatable :: const(:)

  end type s_localres

  ! parameters for allocatable variables
  integer,      public, parameter :: LocalRestraint    = 1

  ! parameters
  integer,      public, parameter :: LocalResBonds     = 1
  integer,      public, parameter :: LocalResAngles    = 2
  integer,      public, parameter :: LocalResDihedrals = 3

  ! subroutines
  public  :: input_localres
  public  :: init_localres
  public  :: alloc_localres
  public  :: dealloc_localres
  private :: read_localres
  private :: check_localres

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_locares
  !> @brief        a driver subroutine for reading local restraint
  !! @authors      CK
  !! @param[in]    localres_filename : filename of local restaraint
  !! @param[out]   localres          : structure of local restraint
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_localres(localres_filename, localres)

    ! formal arguments
    character(*),            intent(in)    :: localres_filename
    type(s_localres),        intent(inout) :: localres
   
    ! local variables
    integer                  :: file
    

    ! open local restraint file
    !
    call open_file(file, localres_filename, IOFileInput)

    ! read restraint data from local restraint file
    !
    call read_localres(file, localres)

    ! close local restraint file
    !
    call close_file(file)

    return

  end subroutine input_localres

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_localres
  !> @brief        initialize local restraint information
  !! @authors      CK
  !! @param[out]   localres : structure of local restraint information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_localres(localres)

    ! formal arguments
    type(s_localres),        intent(inout) :: localres


    localres%num_funcs    = 0

    return

  end subroutine init_localres

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_localres
  !> @brief        allocate local restraint information
  !! @authors      CK
  !! @param[inout] localres : structure of local restraint information
  !! @param[in]    variable : selected variable
  !! @param[in]    var_size : size of the selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_localres(localres, variable, var_size)

    ! formal arguments
    type(s_localres),        intent(inout) :: localres
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size
  
    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat

  
    alloc_stat   = 0
    dealloc_stat = 0


    ! allocate selected variables
    !
    select case (variable)

    case(LocalRestraint)

      if (allocated(localres%func)) then
        if (size(localres%func) == var_size) return
        deallocate(localres%func,             &
                   localres%const,            &
                   localres%ref,              &
                   localres%index_atoms,      &
                   stat = dealloc_stat)
      end if

      allocate(localres%func(var_size),                 &
               localres%const(var_size),                &
               localres%ref(var_size),                  &
               localres%index_atoms(1:4,var_size),      &
               stat = alloc_stat)

      localres%func           (1:var_size) = 0
      localres%const          (1:var_size) = 0.0_wp
      localres%ref            (1:var_size) = 0.0_wp
      localres%index_atoms(1:4,1:var_size) = 0

    case default

      call error_msg('Alloc_localres> bad variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_localres

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_localres
  !> @brief        deallocate local restraint information
  !! @authors      CK
  !! @param[inout] localres : structure of local restraint information
  !! @param[in]    variable : selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_localres(localres, variable)

    ! formal arguments
    type(s_localres),        intent(inout) :: localres
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0


    ! deallocate selected variable
    !  
    select case (variable)

    case(LocalRestraint)

      if (allocated(localres%func)) then
        deallocate(localres%func,             &
                   localres%const,            &
                   localres%ref,              &
                   localres%index_atoms,      &
                   stat = dealloc_stat)
      end if

    case default

      call error_msg('Dealloc_Localres> bad variable')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_localres

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_localres
  !> @brief        read data from local restraint file
  !! @authors      CK
  !! @param[in]    file      : unit number of local restraint file
  !! @param[out]   localres  : localres data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_localres(file, localres)

    ! formal arugments
    integer,                 intent(in)    :: file
    type(s_localres),        intent(inout) :: localres

    ! local variables
    integer                  :: i, ifuncs
    character(80)            :: line
    character(10)            :: chara


    ! deallocate old data
    !
    call dealloc_localres(localres, LocalRestraint)
    call init_localres(localres)


    ! check line count, atom number, segment
    !
    localres%num_funcs = 0

    call check_localres(file, localres%num_funcs)

    ! allocate buffer
    !
    call alloc_localres(localres, LocalRestraint, localres%num_funcs)

    ! read data from local restraint file 
    !
    ifuncs = 0

    do while(.true.)

      read(file, fmt='(a80)', end=200, err=910) line

      line = adjustl(line)

      if (line(1:4) == 'BOND' .or. line(1:4) == 'bond') then
        ifuncs = ifuncs + 1
        read(line,*) chara,(localres%index_atoms(i,ifuncs),i=1,2), &
                     localres%const(ifuncs), localres%ref(ifuncs)
        localres%func(ifuncs) = LocalResBonds

      else if (line(1:4) == 'ANGL' .or. line(1:4) == 'angl') then
        ifuncs = ifuncs + 1
        read(line,*) chara,(localres%index_atoms(i,ifuncs),i=1,3), &
                     localres%const(ifuncs), localres%ref(ifuncs)
        localres%func(ifuncs) = LocalResAngles

      else if (line(1:4) == 'DIHE' .or. line(1:4) == 'dihe') then
        ifuncs = ifuncs + 1
        read(line,*) chara,(localres%index_atoms(i,ifuncs),i=1,4), &
                     localres%const(ifuncs), localres%ref(ifuncs)
        localres%func(ifuncs) = LocalResDihedrals

      end if

    end do

200 continue

    ! write summary of Local restraint information
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Localres> Summary of Data in local restraint'
      write(MsgOut,'(A20,I10)') '  num_funcs       = ', localres%num_funcs
      write(MsgOut,'(A)') ' '
    end if

    return

910 call error_msg('Read_Localres> read error (read line)')

  end subroutine read_localres

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_localres
  !> @brief        check format of local restraint file
  !! @authors      CK
  !! @param[in]    file      : unit number of local restraint file
  !! @param[out]   num_funcs : number of read restraint functions
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_localres(file, num_funcs)

    ! formal arugments
    integer,                 intent(in)    :: file
    integer,                 intent(out)   :: num_funcs
 
    ! local variables
    character(80)            :: line
    integer                  :: num_char
    logical                  :: error_flag


    num_funcs  = 0
    error_flag = .false.

    ! check format of local restraint data
    !
    do while(.true.)

      read(file, fmt='(a80)', end=100, err=910) line

      line = adjustl(line)

      if (line(1:4) == 'BOND' .or. line(1:4) == 'bond') then
        num_funcs = num_funcs + 1
        num_char=split_num(line)
        if (num_char /= 5) error_flag = .true.

      else if (line(1:4) == 'ANGL' .or. line(1:4) == 'angl') then
        num_funcs = num_funcs + 1
        num_char=split_num(line)
        if (num_char /= 6) error_flag = .true.

      else if (line(1:4) == 'DIHE' .or. line(1:4) == 'dihe') then
        num_funcs = num_funcs + 1
        num_char=split_num(line)
        if (num_char /= 7) error_flag = .true.

      end if

    end do

100 rewind(file)

    if (error_flag) &
      call error_msg('Check_Locares> incorrect setting')

    return

910 call error_msg('Check_Locares> read error (read line)')

  end subroutine check_localres

end module fileio_localres_mod
