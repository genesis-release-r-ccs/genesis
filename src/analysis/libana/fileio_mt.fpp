!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_mt_mod
!> @brief   read MotionTree information file
!! @authors Chigusa Kobayashi (CK)
!! @reference  Koike et al. J. Mol. Biol, (2014), vol. 426, pp 752-762
!! @reference  Kobayashi et al. Proteins (2015) 
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../config.h"
#endif

module fileio_mt_mod

  use fileio_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  type, public :: s_mt

    integer                    :: num_rigid_domains
    real(wp)                   :: rigid_domain_crit
    real(wp)                   :: rigid_domain_flexible = 20.0_wp
    integer,      allocatable  :: rigid_domain_list(:)
    real(wp),     allocatable  :: rigid_domain_height(:,:)

  end type s_mt

  ! parameters for allocatable variables
  integer, parameter :: MTRDType = 1

  ! subroutines
  public  :: input_mt
  public  :: alloc_mt
  public  :: dealloc_mt
  private :: read_mt

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_mt
  !> @brief        a driver subroutine for reading MT file
  !! @authors      NT
  !! @param[in]    mt_filename : filename of MT file
  !! @param[inout] mt          : structure of MT data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_mt(mt_filename, mt)

    ! formal arguments
    character(*),            intent(in)    :: mt_filename
    type(s_mt),              intent(inout) :: mt

    ! local variables
    integer                  :: file

    ! open MT file
    !
    call open_file(file, mt_filename, IOFileInput)


    ! read coordinate data from ATP file
    !
    call read_mt(file, mt)


    ! close ATP file
    !
    call close_file(file)

    return

  end subroutine input_mt

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_mt
  !> @brief        allocate Motion Tree data
  !! @authors      CK
  !! @param[inout] mt      : structure of Motion Tree data
  !! @param[in]    variable : selected variable
  !! @param[in]    var_size : size of the selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_mt(mt, variable, var_size1, var_size2)

    ! formal arguments
    type(s_mt),              intent(inout) :: mt
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size1
    integer,                 intent(in)    :: var_size2

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat


    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case(MTRDType)

      if (allocated(mt%rigid_domain_list)) then
        if (size(mt%rigid_domain_list) == var_size2) return
        deallocate(mt%rigid_domain_list,    &
                   mt%rigid_domain_height,  &
                   stat = dealloc_stat)
      end if

      allocate(mt%rigid_domain_list(var_size2),&
               mt%rigid_domain_height(var_size1, var_size1), &
               stat = alloc_stat)

      mt%rigid_domain_list(1:var_size2)               = 0
      mt%rigid_domain_height(1:var_size1,1:var_size1) = 0.0_wp

    case default

      call error_msg('Alloc_MT> bad variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return 

  end subroutine alloc_mt

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_mt
  !> @brief        deallocate MT data
  !! @authors      CK
  !! @param[inout] mt      : structure of Motion Tree data
  !! @param[in]    variable : selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_mt(mt, variable)

    ! formal arguments
    type(s_mt),              intent(inout) :: mt
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0

    ! deallocate selected variables
    !
    select case (variable)

    case(MTRDType)

      if (allocated(mt%rigid_domain_list)) then
        deallocate(mt%rigid_domain_list,    &
                   mt%rigid_domain_height,  &
                   stat = dealloc_stat)
      end if

    case default

      call error_msg('Dealloc_MT> bad variable')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return 

  end subroutine dealloc_mt

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_mt
  !> @brief        read data from Motion Tree file
  !! @authors      CK
  !! @param[in]    file : file unit number
  !! @param[inout] mt   : structure of MT data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_mt(file, mt)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_mt),              intent(inout) :: mt

    ! local variables
    integer                  :: i, j, num_residues
    real(wp)                 :: height


    read(file,*) mt%num_rigid_domains, num_residues, mt%rigid_domain_crit

    call alloc_mt(mt, MTRDType, mt%num_rigid_domains, num_residues)

    read(file,*) (mt%rigid_domain_list(j), j = 1, num_residues)

    do i = 1, mt%num_rigid_domains
      read(file,*) (mt%rigid_domain_height(j,i), j=1,mt%num_rigid_domains)
    end do

    do i = 1, mt%num_rigid_domains
      do j = i+1, mt%num_rigid_domains
        height = mt%rigid_domain_crit/mt%rigid_domain_height(j,i)
        mt%rigid_domain_height(i,j)=height
        mt%rigid_domain_height(j,i)=height
      end do
      mt%rigid_domain_height(i,i)=1.0_wp
    end do
    mt%rigid_domain_flexible = mt%rigid_domain_crit/mt%rigid_domain_flexible

    ! write summary of MT information
    !
    if (main_rank) then

      write(MsgOut,'(a)') 'Read_MT> Summary of MT file'
      write(MsgOut,'(a20,i10)') &
           '  num_rigid_domains  = ', mt%num_rigid_domains
      write(MsgOut,'(a)') ''

    end if

    return 

  end subroutine read_mt

end module fileio_mt_mod
