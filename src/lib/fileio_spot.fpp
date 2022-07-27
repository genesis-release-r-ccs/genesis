!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_spot_mod
!> @brief   read coordinates and radius of spherical potential functions
!! @authors Kiyoshi Yagi (KY)
!
!  (c) Copyright 2016 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fileio_spot_mod

  use fileio_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_spot

    integer                       :: ncenter = 0
    real(wp),         allocatable :: center(:,:)
    real(wp),         allocatable :: radius(:)

  end type s_spot

  ! subroutines
  public  :: input_spot
!  public  :: alloc_spot
  public  :: dealloc_spot
  public  :: read_spot

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_spot
  !> @brief        open, read, and close spotfile 
  !! @authors      KY
  !! @param[in]    spot_filename : filename of spotfile
  !! @param[out]   spot          : spherical potential information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_spot(spot_filename, spot)

    ! formal arguments
    character(*),            intent(in)    :: spot_filename
    type(s_spot),            intent(inout) :: spot

    ! local variables
    integer                  :: unit_no


    ! open spotfile
    !
    call open_file(unit_no, spot_filename, IOFileInput)

    ! read spotfile
    !
    call read_spot(unit_no, spot)

    ! close spotfile
    !
    call close_file(unit_no)


    if (main_rank) then
      write(MsgOut,'(A)') 'Input_spot> Summary of spotfile'
      write(MsgOut,'(A20,I10)') '  num_center      = ', spot%ncenter
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine input_spot

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_spot
  !> @brief        deallocate spherical potential information
  !! @authors      KY
  !! @param[inout] spot     : spherical potential information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_spot(spot)

    ! formal arguments
    type(s_spot),            intent(inout) :: spot

    ! local variables
    integer                  :: dealloc_stat
    

    dealloc_stat = 0

    if (allocated(spot%center)) then
      deallocate(spot%center, &
                 spot%radius, &
                 stat = dealloc_stat)
    end if

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_spot

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_spot
  !> @brief        read data from spotfile
  !! @authors      KY
  !! @param[in]    unit_no : unit number of spotfile
  !! @param[out]   spot    : spherical potential information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_spot(unit_no, spot)

    ! formal arguments
    integer,                 intent(in)    :: unit_no
    type(s_spot),            intent(inout) :: spot

    ! local variables
    integer                  :: i, nc
    character(6)             :: an


    read(unit_no, *) spot%ncenter
    read(unit_no, *) 

    nc = spot%ncenter
    allocate(spot%center(3,nc), spot%radius(nc))
    do i = 1, spot%ncenter
      read(unit_no, *) an, spot%center(:,i), spot%radius(i)
    end do

    return

  end subroutine read_spot

end module fileio_spot_mod
