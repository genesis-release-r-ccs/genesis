!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_mrc_mod
!> @brief   read mrc density map file
!! @authors Takaharu Mori (TM), Daisuke Matsuoka (DM)
!
!  (c) Copyright 2016 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fileio_mrc_mod

  use fileio_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_mrc
    integer(kind=4)       :: nx, ny, nz, mode
    integer(kind=4)       :: nxstart, nystart, nzstart
    integer(kind=4)       :: mx, my, mz
    integer(kind=4)       :: mapc, mapr, maps
    integer(kind=4)       :: ispg, nsymbt, extra, nlabl
    integer(kind=1)       :: machst(4)
    real(sp)              :: cella(1:3), cellb(1:3), origin(1:3)
    real(sp)              :: dmin, dmax, dmean, rms
    character(4)          :: map
    character(1)          :: label
    real(sp), allocatable :: map_value(:)
  end type s_mrc

  ! subroutines
  public  :: input_mrc
  public  :: output_mrc
  public  :: alloc_mrc
  public  :: dealloc_mrc
  public  :: read_mrc
  private :: write_mrc

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_mrc
  !> @brief        open, read, and close mrcfile 
  !! @authors      TM
  !! @param[in]    mrc_filename : filename of mrcfile
  !! @param[out]   mrc          : EM data information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_mrc(mrc_filename, mrc)

    ! formal arguments
    character(*),            intent(in)    :: mrc_filename
    type(s_mrc),             intent(inout) :: mrc

    ! local variables
    integer                  :: unit_no
    integer(kind=1)          :: machst


    ! open mrcfile
    !
    call open_binary_file(unit_no, mrc_filename, IOFileInput,     &
                          IOFileLittleEndian, IOFileStreamAccess)

    ! Check endianess in mrcfile
    !
    read(unit_no,pos=213) machst
    if (machst == 68) then
      ! little endian
      rewind(unit_no)
    else if (machst == 17) then
      ! big endian
      call close_file(unit_no)
      call open_binary_file(unit_no, mrc_filename, IOFileInput,     &
                            IOFileBigEndian,  IOFileStreamAccess)
    end if

    ! read mrcfile
    !
    call read_mrc(unit_no, mrc)

    ! close mrcfile
    !
    call close_file(unit_no)

    ! write sumary
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Input_Emap> Summary of mrcfile'
      write(MsgOut,'(A20,3I10  )') '  NX NY NZ        = ', mrc%nx, mrc%ny, mrc%nz
      write(MsgOut,'(A20, I10  )') '  MODE            = ', mrc%mode
      write(MsgOut,'(A20,3I10  )') '  NSTART X Y Z    = ', mrc%nxstart, mrc%nystart, mrc%nzstart
      write(MsgOut,'(A20,3I10  )') '  MX MY MZ        = ', mrc%mx, mrc%my, mrc%mz
      write(MsgOut,'(A20,3F10.3)') '  CELLA           = ', mrc%cella(1:3)
      write(MsgOut,'(A20,3F10.3)') '  CELLB           = ', mrc%cellb(1:3)
      write(MsgOut,'(A20,3I10  )') '  MAPC MAPR MAPS  = ', mrc%mapc, mrc%mapr, mrc%maps
      write(MsgOut,'(A20,3F10.3)') '  DMIN DMAX DMEAN = ', mrc%dmin, mrc%dmax, mrc%dmean
      write(MsgOut,'(A20,3F10.3)') '  ORIGIN X Y Z    = ', mrc%origin(1:3)
      write(MsgOut,'(A20, I10  )') '  ISPG            = ', mrc%ispg
      write(MsgOut,'(A20, F10.3)') '  RMS             = ', mrc%rms
      write(MsgOut,'(A)') ''
    end if

    return

  end subroutine input_mrc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_mrc
  !> @brief        open, write, and close mrcfile
  !! @authors      DM
  !! @param[in]    mrc_filename : filename of mrcfile
  !! @param[out]   mrc          : EM data information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_mrc(mrc_filename, mrc)

    ! formal arguments
    character(*),            intent(in)    :: mrc_filename
    type(s_mrc),             intent(inout) :: mrc

    ! local variables
    integer                  :: unit_no


    ! open mrcfile
    !
    call open_binary_file(unit_no, mrc_filename, IOFileOutputNew,   &
                          IOFileLittleEndian, IOFileStreamAccess)

    ! rwrite mrcfile
    !
    call write_mrc(unit_no, mrc)

    ! close mrcfile
    !
    call close_file(unit_no)

    ! write sumary
    !
    write(MsgOut,'(A)') 'Input_Emap> Summary of mrcfile'
    write(MsgOut,'(A20,3I10  )') '  NX NY NZ        = ', mrc%nx, mrc%ny, mrc%nz
    write(MsgOut,'(A20, I10  )') '  MODE            = ', mrc%mode
    write(MsgOut,'(A20,3I10  )') '  NSTART X Y Z    = ', mrc%nxstart, mrc%nystart, mrc%nzstart
    write(MsgOut,'(A20,3I10  )') '  MX MY MZ        = ', mrc%mx, mrc%my, mrc%mz
    write(MsgOut,'(A20,3F10.3)') '  CELLA           = ', mrc%cella(1:3)
    write(MsgOut,'(A20,3F10.3)') '  CELLB           = ', mrc%cellb(1:3)
    write(MsgOut,'(A20,3I10  )') '  MAPC MAPR MAPS  = ', mrc%mapc, mrc%mapr, mrc%maps
    write(MsgOut,'(A20,3F10.3)') '  DMIN DMAX DMEAN = ', mrc%dmin, mrc%dmax, mrc%dmean
    write(MsgOut,'(A20,3F10.3)') '  ORIGIN X Y Z    = ', mrc%origin(1:3)
    write(MsgOut,'(A20, I10  )') '  ISPG            = ', mrc%ispg
    write(MsgOut,'(A20, F10.3)') '  RMS             = ', mrc%rms
    write(MsgOut,'(A)') ''

    return

  end subroutine output_mrc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_mrc
  !> @brief        allocate EM data information
  !! @authors      TM
  !! @param[inout] mrc      : EM data information
  !! @param[in]    var_size : allocation size
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_mrc(mrc, var_size)

    ! formal arguments
    type(s_mrc),             intent(inout) :: mrc
    integer,                 intent(in)    :: var_size
    
    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat

    
    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    if (allocated(mrc%map_value)) then
      if (size(mrc%map_value(:)) == var_size) &
        return
      deallocate(mrc%map_value,      &
                 stat = dealloc_stat)
    end if

    allocate(mrc%map_value(var_size), &
             stat = alloc_stat)

    mrc%map_value(:) = 0.0d0

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_mrc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_mrc
  !> @brief        deallocate EM data information
  !! @authors      TM
  !! @param[inout] mrc      : EM data information
  !! @param[in]    variable : an variable to be allocated 
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_mrc(mrc)

    ! formal arguments
    type(s_mrc),             intent(inout) :: mrc

    ! local variables
    integer                  :: dealloc_stat
    

    dealloc_stat = 0

    if (allocated(mrc%map_value)) then
      deallocate(mrc%map_value,       &
                 stat = dealloc_stat)
    end if

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_mrc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_mrc
  !> @brief        read data from mrcfile
  !! @authors      TM
  !! @param[in]    unit_no : unit number of mrcfile
  !! @param[out]   mrc     : EM data inforation
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_mrc(unit_no, mrc)

    ! formal arguments
    integer,                 intent(in)    :: unit_no
    type(s_mrc),             intent(inout) :: mrc

    ! local variables
    integer                  :: i, ix, iy, iz, idx


    read(unit_no) mrc%nx, mrc%ny, mrc%nz
    read(unit_no) mrc%mode
    read(unit_no) mrc%nxstart, mrc%nystart, mrc%nzstart
    read(unit_no) mrc%mx, mrc%my, mrc%mz
    read(unit_no) mrc%cella(1:3), mrc%cellb(1:3)
    read(unit_no) mrc%mapc, mrc%mapr, mrc%maps
    read(unit_no) mrc%dmin, mrc%dmax, mrc%dmean
    read(unit_no) mrc%ispg, mrc%nsymbt
    do i = 25, 49
      read(unit_no) mrc%extra
    end do
    read(unit_no) mrc%origin(1:3)
    read(unit_no) mrc%map
    read(unit_no) mrc%machst
    read(unit_no) mrc%rms
    read(unit_no) mrc%nlabl
    do i = 1, 800
      read(unit_no) mrc%label
    end do

    call alloc_mrc(mrc, mrc%nx*mrc%ny*mrc%nz)

    do iz = 0, mrc%nz-1
      do iy = 0, mrc%ny-1
        do ix = 0, mrc%nx-1
          idx = 1 + ix + iy*mrc%nx + iz*mrc%nx*mrc%ny
          read(unit_no) mrc%map_value(idx)
        end do
      end do
    end do

    return

  end subroutine read_mrc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_mrc
  !> @brief        write data from mrcfile
  !! @authors      DM
  !! @param[in]    unit_no : unit number of mrcfile
  !! @param[out]   mrc     : EM data inforation
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_mrc(unit_no, mrc)

    ! formal arguments
    integer,                 intent(in)    :: unit_no
    type(s_mrc),             intent(inout) :: mrc

    ! local variables
    integer                  :: i, ix, iy, iz, idx

    integer(kind=1)       :: machst(4)
    data machst / z'44', z'41', z'00', z'00' /

    ! constants
    mrc%mode  = 2
    mrc%map   = 'MAP '
    mrc%ispg  = 1
    mrc%mapc  = 1
    mrc%mapr  = 2
    mrc%maps  = 3
    mrc%extra = 0
    mrc%nlabl = 0
    mrc%label = char(0)

    ! specify endianness used in an output binary file
    ! (stamp for little endian)
    !
    mrc%machst = machst

    ! write map data
    !
    write(unit_no) mrc%nx, mrc%ny, mrc%nz
    write(unit_no) mrc%mode
    write(unit_no) mrc%nxstart, mrc%nystart, mrc%nzstart
    write(unit_no) mrc%mx, mrc%my, mrc%mz
    write(unit_no) mrc%cella(1:3), mrc%cellb(1:3)
    write(unit_no) mrc%mapc, mrc%mapr, mrc%maps
    write(unit_no) mrc%dmin, mrc%dmax, mrc%dmean
    write(unit_no) mrc%ispg, mrc%nsymbt

    do i = 25, 49
      write(unit_no) mrc%extra
    end do

    write(unit_no) mrc%origin(1:3)
    write(unit_no) mrc%map
    write(unit_no) mrc%machst
    write(unit_no) mrc%rms

    write(unit_no) mrc%nlabl
    do i = 1, 800
      write(unit_no) mrc%label
    end do

    ! density value
    do iz = 0, mrc%nz-1
      do iy = 0, mrc%ny-1
        do ix = 0, mrc%nx-1
          idx = 1 + ix + iy*mrc%nx + iz*mrc%nx*mrc%ny
          write(unit_no) mrc%map_value(idx)
        end do
      end do
    end do

    return

  end subroutine write_mrc

end module fileio_mrc_mod
