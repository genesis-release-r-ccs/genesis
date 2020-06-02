!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_sit_mod
!> @brief   read situs density map file
!! @authors Takaharu Mori (TM)
!
!  (c) Copyright 2016 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fileio_sit_mod

  use fileio_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_sit
    real(wp), allocatable  :: map_data(:,:,:)
    real(wp)               :: dx
    real(wp)               :: dy
    real(wp)               :: dz
    real(wp)               :: x0
    real(wp)               :: y0
    real(wp)               :: z0
    integer                :: nx
    integer                :: ny
    integer                :: nz
  end type s_sit

  ! subroutines
  public  :: input_sit
  public  :: alloc_sit
  public  :: dealloc_sit
  public  :: read_sit

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_sit
  !> @brief        open, read, and close sitfile 
  !! @authors      TM
  !! @param[in]    sit_filename : filename of sitfile
  !! @param[out]   sit          : EM data information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_sit(sit_filename, sit)

    ! formal arguments
    character(*),            intent(in)    :: sit_filename
    type(s_sit),             intent(inout) :: sit

    ! local variables
    integer                  :: unit_no

    ! open sitfile
    !
    call open_file(unit_no, sit_filename, IOFileInput)

    ! read sitfile
    !
    call read_sit(unit_no, sit)

    ! close sitfile
    !
    call close_file(unit_no)

    ! write sumary
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Input_Emap> Summary of sitfile'
      write(MsgOut,'(A20,F10.3)') '  voxel size x    = ', sit%dx
      write(MsgOut,'(A20,F10.3)') '  voxel size y    = ', sit%dy
      write(MsgOut,'(A20,F10.3)') '  voxel size z    = ', sit%dz
      write(MsgOut,'(A20,I10  )') '  num x increments= ', sit%nx
      write(MsgOut,'(A20,I10  )') '  num y increments= ', sit%ny
      write(MsgOut,'(A20,I10  )') '  num z increments= ', sit%nz
      write(MsgOut,'(A20,F10.3)') '  first voxel xcrd= ', sit%x0
      write(MsgOut,'(A20,F10.3)') '  first voxel ycrd= ', sit%y0
      write(MsgOut,'(A20,F10.3)') '  first voxel zcrd= ', sit%z0
      write(MsgOut,'(A)') ''
    end if

    return

  end subroutine input_sit

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_sit
  !> @brief        allocate EM data information
  !! @authors      TM
  !! @param[inout] sit     : EM data information
  !! @param[in]    var_size : allocation size
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_sit(sit, var_size1, var_size2, var_size3)

    ! formal arguments
    type(s_sit),             intent(inout) :: sit
    integer,                 intent(in)    :: var_size1
    integer,                 intent(in)    :: var_size2
    integer,                 intent(in)    :: var_size3
    
    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat

    
    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    if (allocated(sit%map_data)) then
      if (size(sit%map_data(:,0,0)) == var_size1) &
        return
      deallocate(sit%map_data,      &
                 stat = dealloc_stat)
    end if

    allocate(sit%map_data(0:var_size1-1, 0:var_size2-1, 0:var_size3-1), &
             stat = alloc_stat)

    sit%map_data(0:var_size1-1,0:var_size2-1,0:var_size3-1) = 0.0_wp

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_sit

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_sit
  !> @brief        deallocate EM data information
  !! @authors      TM
  !! @param[inout] sit     : EM data information
  !! @param[in]    variable : an variable to be allocated 
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_sit(sit)

    ! formal arguments
    type(s_sit),             intent(inout) :: sit

    ! local variables
    integer                  :: dealloc_stat
    

    dealloc_stat = 0

    if (allocated(sit%map_data)) then
      deallocate(sit%map_data,       &
                 stat = dealloc_stat)
    end if

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_sit

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_sit
  !> @brief        read data from sitfile
  !! @authors      TM
  !! @param[in]    unit_no : unit number of sitfile
  !! @param[out]   sit    : EM data inforation
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_sit(unit_no, sit)

    ! formal arguments
    integer,                 intent(in)    :: unit_no
    type(s_sit),             intent(inout) :: sit

    ! local variables
    real(wp)       :: dx, x0, y0, z0
    integer        :: nx, ny, nz

    read(unit_no, *) dx, x0, y0, z0, nx, ny, nz

    sit%dx = dx
    sit%dy = dx
    sit%dz = dx
    sit%x0 = x0
    sit%y0 = y0
    sit%z0 = z0
    sit%nx = nx
    sit%ny = ny
    sit%nz = nz

    call alloc_sit(sit, nx, ny, nz)

    read(unit_no, *) sit%map_data

    return

  end subroutine read_sit

end module fileio_sit_mod
