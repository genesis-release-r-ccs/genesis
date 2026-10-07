!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_setup_mpi_mod
!> @brief   Sets up MPI
!! @authors Takaharu Mori (TM), Jaewoon Jung (JJ)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_setup_mpi_mod

  use cg_energy_mod
  use cg_energy_str_mod
  use fileio_control_mod
  use math_libs_mod
  use messages_mod
  use mpi_parallel_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif
#ifdef OMP
    use omp_lib, only: omp_get_num_threads
#endif

  implicit none
  private

  ! subroutines
  public  :: setup_mpi_md

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_mpi_md
  !> @brief        Sets up MPI in MD
  !> @authors      JJ, TM
  !! @param[in]    ene_info : ENERGY section parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_mpi_md(ene_info)
  
    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info

    ! local variables
    integer                  :: alloc_stat, dealloc_stat
    integer                  :: nthread_actual
    logical                  :: prime_flag

    integer,     allocatable :: globrank(:), locrank(:)
    logical,     allocatable :: nonbreal(:), nonbrecip(:)


#ifdef HAVE_MPI_GENESIS
      ! Equalize communicators (comm_country & comm_city) between MD
      !
      mpi_comm_country = mpi_comm_world
      nproc_country    = nproc_world
      my_country_rank  = my_world_rank

      ! check factorization of 2,3,5
      if (nproc_country /= 1) then
        prime_flag=factorization_235(nproc_country)
        if (.not. prime_flag)  &
          call error_msg('Setup_Mpi_Md> number of MPI processes should be'//&
                         ' multiples of 2,3,5')
      end if
#endif

    select case (ene_info%electrostatic)

    case (ElectrostaticPME)

#ifdef HAVE_MPI_GENESIS
      mpi_comm_city = mpi_comm_country
#endif

      real_calc       = .true.
      reciprocal_calc = .true.
      nproc_city      = nproc_country
      my_city_rank    = my_country_rank

    case (ElectrostaticCUTOFF)

#ifdef HAVE_MPI_GENESIS
      mpi_comm_city = mpi_comm_country
#endif

      real_calc       = .true.
      reciprocal_calc = .false.
      nproc_city      = nproc_country
      my_city_rank    = my_country_rank

    case (ElectrostaticRF)

#ifdef HAVE_MPI_GENESIS
      mpi_comm_city = mpi_comm_country
#endif

      real_calc       = .true.
      reciprocal_calc = .false.
      nproc_city      = nproc_country
      my_city_rank    = my_country_rank

    end select

#ifdef HAVE_MPI_GENESIS
    ! Write the summary of setup MPI
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
      !$omp parallel shared(nthread_actual)
      nthread_actual = omp_get_num_threads()
      !$omp end parallel
      if (nthread_actual /= nthread) then
        write(MsgOut,'(a,i10)') '  expected OpenMP threads = ', nthread
        write(MsgOut,'(a,i10)') '  actual OpenMP threads   = ', nthread_actual
        call error_msg('Setup_Mpi_Md> OpenMP thread count mismatch')
      end if

      write(MsgOut,'(a,i10)') '  number of OpenMP threads  = ', nthread
#else
      nthread = 1
#endif

      write(MsgOut,'(a,i10)') '  total number of CPU cores = ', &
                              nproc_world*nthread
      write(MsgOut,'(a)') ''

    end if


    ! deallocate local array
    !
    deallocate(globrank, locrank, nonbreal, nonbrecip, stat=dealloc_stat)
    if (dealloc_stat /= 0) &
      call error_msg_dealloc

#else

    write(MsgOut,'(a)') 'Setup_Mpi_Md> Single CPU is used.'
    write(MsgOut,'(a)') ''

#endif

    return

  end subroutine setup_mpi_md

end module cg_setup_mpi_mod
