!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sa_setup_mpi_mod
!> @brief   setup MPI
!! @authors Jaewoon Jung (JJ), Isseki Yu (IY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module sa_setup_mpi_mod

  use fileio_control_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: setup_mpi_sa

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_mpi_sa
  !> @brief        setup mpi in GA
  !> @authors      JJ, TM, IY
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_mpi_sa()

    ! local variables
    integer                  :: alloc_stat, dealloc_stat

    integer,     allocatable :: globrank(:), locrank(:)
    logical,     allocatable :: nonbreal(:), nonbrecip(:)


#ifdef HAVE_MPI_GENESIS
    ! Equalize communicators (comm_country & comm_city) between MD
    !
    mpi_comm_country = mpi_comm_world
    nproc_country    = nproc_world
    my_country_rank  = my_world_rank
#endif

#ifdef HAVE_MPI_GENESIS
    reciprocal_calc = .true.
    mpi_comm_city = mpi_comm_country
    nproc_city = nproc_country
    my_city_rank = my_country_rank
#else
    real_calc       = .true.
    reciprocal_calc = .true.
    nproc_city      = nproc_country
    my_city_rank    = my_country_rank
#endif

#ifdef HAVE_MPI_GENESIS

    ! Write the summary of setup MP
    !
    allocate(globrank (nproc_world), &
             locrank  (nproc_world), &
             nonbreal (nproc_world), &
             nonbrecip(nproc_world), stat=alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    call mpi_gather(my_world_rank,   1, mpi_integer,                     &
                    globrank,  1, mpi_integer, 0, mpi_comm_world, ierror)
    call mpi_gather(my_city_rank, 1, mpi_integer,                        &
                    locrank,   1, mpi_integer, 0, mpi_comm_world, ierror)
    call mpi_gather(real_calc, 1, mpi_logical,                           &
                    nonbreal,  1, mpi_logical, 0, mpi_comm_world, ierror)
    call mpi_gather(reciprocal_calc, 1, mpi_logical,                     &
                    nonbrecip, 1, mpi_logical, 0, mpi_comm_world, ierror)

    if (main_rank) then

      write(MsgOut,'(a)') 'Setup_Mpi_Md> Summary of Setup MPI'
      write(MsgOut,'(a,i10)') '  number of MPI processes   = ', nproc_world

#ifdef OMP
      write(MsgOut,'(a,i10)') '  number of OpenMP threads  = ', nthread
#endif

      write(MsgOut,'(a,i10)') '  total number of CPU cores = ', &
                              nproc_world*nthread
      write(MsgOut,'(a)') ''

    end if


    ! deallocate local array
    !
    deallocate(globrank, locrank, nonbreal, nonbrecip, stat=dealloc_stat)
    if (dealloc_stat .ne. 0) &
      call error_msg_dealloc

#else

    write(MsgOut,'(a)') 'Setup_Mpi_Md> Single CPU is used.'
    write(MsgOut,'(a)') ''

#endif

    return

  end subroutine setup_mpi_sa

end module sa_setup_mpi_mod
