!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_mode_mod
!> @brief   read principal component vector from MODE file
!! @authors Yuji Sugita (YS), Yasuaki Komuro (YK) 
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fileio_mode_mod

  use fileio_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_mode

    integer               :: num_pc_modes   = 0
    real(wp), allocatable :: pc_mode(:)

  end type s_mode

  ! subroutines
  public  :: input_mode
  public  :: output_mode
  public  :: init_mode
  public  :: alloc_mode
  public  :: dealloc_mode
  private :: read_mode
  private :: write_mode
  private :: check_mode

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_mode
  !> @brief        a driver subroutine for reading MODE file
  !! @authors      YS, YK
  !! @param[in]    mode_filename : filename of MODE file
  !! @param[out]   mode          : structure of MODE information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_mode(mode_filename, mode)

    ! formal arguments
    character(*),            intent(in)    :: mode_filename
    type(s_mode),            intent(inout) :: mode
   
    ! local variables
    integer                  :: file
    

    ! open MODE file
    !
    call open_file(file, mode_filename, IOFileInput)

    ! read principal component vector from MODE file
    !
    call read_mode(file, mode)

    ! close MODE file
    !
    call close_file(file)

    return

  end subroutine input_mode

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_mode
  !> @brief        a driver subroutine for writing MODE file
  !! @authors      YS, YK
  !! @param[in]    mode_filename : filename of MODE file
  !! @param[in]    mode          : structure of MODE information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_mode(mode_filename, mode)

    ! formal arguments
    character(*),            intent(in) :: mode_filename
    type(s_mode),            intent(in) :: mode
   
    ! local variables
    integer                  :: file
    

    ! open MODE file
    !
    call open_file(file, mode_filename, IOFileOutputNew)

    ! write coordinate data from MODE file
    !
    call write_mode(file, mode)

    ! close MODE file
    !
    call close_file(file)

    return

  end subroutine output_mode

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_mode
  !> @brief        initialize MODE information
  !! @authors      YS, YK
  !! @param[out]   mode : structure of MODE information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_mode(mode)

    ! formal arguments
    type(s_mode),             intent(inout) :: mode

    mode%num_pc_modes         = 0

    return

  end subroutine init_mode

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_mode
  !> @brief        allocate MODE information
  !! @authors      YS, YK
  !! @param[inout] mode     : structure of MODE information
  !! @param[in]    variable : selected variable
  !! @param[in]    var_size : size of the selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_mode(mode, var_size)

    ! formal arguments
    type(s_mode),            intent(inout) :: mode
    integer,                     intent(in)    :: var_size
  
    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat
  
    ! initialize
    !
    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate variable
    !
    if (allocated(mode%pc_mode)) then
      if (size(mode%pc_mode) == var_size) return
      deallocate(mode%pc_mode, stat = dealloc_stat)
    end if

    allocate(mode%pc_mode(1:var_size), stat = alloc_stat)

    mode%pc_mode(1:var_size) = 0.0_wp


    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_mode

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_mode
  !> @brief        deallocate MODE information
  !! @authors      YS
  !! @param[inout] mode     : structure of MODE information
  !! @param[in]    variable : selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_mode(mode)

    ! formal arguments
    type(s_mode),            intent(inout) :: mode

    ! local variables
    integer                  :: dealloc_stat

    ! initialize
    !
    dealloc_stat = 0

    ! deallocate selected variable
    !  
    if (allocated(mode%pc_mode)) then
      deallocate (mode%pc_mode, stat = dealloc_stat)
    end if

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_mode

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_mode
  !> @brief        read data from MODE file
  !! @authors      YS, YK
  !! @param[in]    file : unit number of MODE file
  !! @param[out]   mode : mode data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_mode(file, mode)

    ! formal arugments
    integer,                 intent(in)    :: file
    type(s_mode),            intent(inout) :: mode

    ! local variables
    integer                  :: i, ios, ipc, num_read_pcs
    real(wp)                 :: line


    ! deallocate old data
    !
    call dealloc_mode(mode)
    call init_mode(mode)


    ! check line count
    !
    num_read_pcs = 0

    call check_mode(file, num_read_pcs)

    ! allocate
    !
    call alloc_mode(mode, num_read_pcs)

    ! read data from MODE file 
    !
    ios  = 1 
    ipc  = 0

    do while(.true.)
      read(file, *, end=1000, err=915) line
      ipc = ipc + 1
      mode%pc_mode(ipc) = line
    end do
1000 rewind(file)

  
    ! check 
    !
    if (mod(ipc,3) /= 0) &
      call error_msg('Read_Mode> read error (number of data)')

    mode%num_pc_modes = ipc

    ! write summary of MODE information
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Mode> Summary of Data in MODE file'
      write(MsgOut,'(A20,I10)') '  num_pc_modes    = ', mode%num_pc_modes
      write(MsgOut,'(A)') ' '
      do i = 1, mode%num_pc_modes
        !write(MsgOut,'(A10,E15.7E3)') 'mode = ', mode%pc_mode(i)
        write(MsgOut,*) 'mode = ', mode%pc_mode(i)
      end do
      write(MsgOut,'(A)') ' '
    end if

    return

915 call error_msg('Read_Pdb> read error (read line)')
916 call error_msg('Read_Pdb> read error (real mode)')

  end subroutine read_mode

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_mode
  !> @brief        write data to MODE file
  !! @authors      YS, YK
  !! @param[in]    file : unit number of MODE file
  !! @param[in]    mode : structure of MODE information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_mode(file, mode)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_mode),            intent(in)    :: mode

    ! local variables
    integer                  :: i, j, len
    character(80)            :: fmt_a, fmt_t
    character(6)             :: crec
    character(4)             :: catm, cres, cstr, cseg
    logical                  :: use_cid
    

    if (.not.allocated(mode%pc_mode)) &
      call error_msg('Out_Mode> not allocated: mode%pc_mode')

    do i = 1, mode%num_pc_modes
      write(file, *) mode%pc_mode(i)
    enddo

    return 

  end subroutine write_mode

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_mode
  !> @brief        check format of MODE file
  !! @authors      YS, YK
  !! @param[in]    file           : unit number of MODE file
  !! @param[out]   num_read_pcs   : number of read principal components
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_mode(file, num_read_pcs)

    ! formal arugments
    integer,                 intent(in)    :: file
    integer,                 intent(out)   :: num_read_pcs
 
    ! local variables
    character(100)            :: line

    ! initialize
    !
    num_read_pcs = 0

    ! count of number of pc
    !
    do while(.true.)
      read(file, '(a100)', end=999, err=917) line
      num_read_pcs = num_read_pcs + 1
    end do

999 rewind(file) 

    return

917 call error_msg('Check_Mode> read error (read line)')

  end subroutine check_mode

end module fileio_mode_mod
