!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_namd_xyz_mod
!> @brief   NAMD xyz file I/O
!! @authors Norio Takase (NT), Chigusa Kobayashi (CK)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module fileio_namd_xyz_mod

  use fileio_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_namd_xyz

    real(dp), allocatable :: v(:,:)

  end type s_namd_xyz

  ! subroutines
  public  :: input_namd_xyz
  public  :: output_namd_xyz
  public  :: alloc_namd_xyz
  public  :: dealloc_namd_xyz

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_namd_xyz
  !> @brief        input NAMD binary coordinate or velocity
  !! @authors      NT
  !! @param[in]    filename : file name
  !! @param[inout] xyz      : coordinates or velocities information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_namd_xyz(filename, xyz)

    ! formal arguments
    character(*),            intent(in)    :: filename
    type(s_namd_xyz),        intent(inout) :: xyz

    ! local variables
    real(wp)                 :: v(3)
    integer                  :: unit_no, alloc_stat, dealloc_stat
    integer                  :: i, j, natom


    unit_no = get_unit_no()
    if (unit_no == InvalidUnitNo) &
      goto 900
    open(unit_no, file=filename, status='old', form='unformatted',  &
         access='stream', err=900)

    read(unit_no,err=900) natom
    call alloc_namd_xyz(xyz, natom)

    do i = 1, natom
      read(unit_no,err=900,end=900) (v(j), j = 1,3)
      xyz%v(:,i) = v(:)
    end do

    close(unit_no)
    call free_unit_no(unit_no)

    return

900 write(ErrOut,'(A)') 'Input_Namd_XYZ> ', trim(filename)
    call error_msg('Input_Namd_Xyz> I/O Error.')

  end subroutine input_namd_xyz

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_namd_xyz
  !> @brief        output NAMD binary coordinate or velocity
  !! @authors      CK
  !! @param[in]    filename : file name
  !! @param[in]    xyz      : coordinates or velocities information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_namd_xyz(filename, xyz)

    ! formal arguments
    character(*),            intent(in)    :: filename
    type(s_namd_xyz),        intent(in)    :: xyz

    ! local variables
    real(wp)                 :: v(3)
    integer                  :: unit_no
    integer                  :: i, j, natom


    unit_no = get_unit_no()
    if (unit_no == InvalidUnitNo) &
      goto 900
    open(unit_no, file=filename, status='new', form='unformatted',  &
         access='stream', err=900)

    natom = size(xyz%v(1,:))
    write(unit_no,err=900) natom

    do i = 1, natom
      v(:) = xyz%v(:,i)
      write(unit_no,err=900) (v(j), j = 1,3)
    end do

    close(unit_no)
    call free_unit_no(unit_no)

    return

900 write(ErrOut,'(A)') 'Output_Namd_XYZ> ', trim(filename)
    call error_msg('Output_Namd_Xyz> I/O Error.')

  end subroutine output_namd_xyz

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_namd_xyz
  !> @brief        allocate memory
  !! @authors      CK
  !! @param[inout] xyz : coordinates or velocities information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_namd_xyz(xyz, var_size)

    ! formal arguments
    type(s_namd_xyz),        intent(inout) :: xyz
    integer,                 intent(in)    :: var_size

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat


    ! initialize
    !
    dealloc_stat = 0
      
    if (allocated(xyz%v)) &
      deallocate(xyz%v, stat = dealloc_stat)

    allocate(xyz%v(1:3,1:var_size),  stat = alloc_stat)

    if (dealloc_stat /= 0) &
      call error_msg_dealloc
    if (alloc_stat /= 0) &
      call error_msg_alloc

    return

  end subroutine alloc_namd_xyz

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_namd_xyz
  !> @brief        deallocate memory
  !! @authors      NT
  !! @param[inout] xyz : coordinates or velocities information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_namd_xyz(xyz)

    ! formal arguments
    type(s_namd_xyz),        intent(inout) :: xyz

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0
    if (allocated(xyz%v)) &
      deallocate(xyz%v, stat = dealloc_stat)

    if (dealloc_stat /= 0) &
      call error_msg_dealloc

    return

  end subroutine dealloc_namd_xyz

end module fileio_namd_xyz_mod
