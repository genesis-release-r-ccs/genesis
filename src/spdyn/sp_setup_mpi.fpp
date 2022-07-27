!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_setup_mpi_mod
!> @brief   setup MPI
!! @authors Takaharu Mori (TM), Jaewoon Jung (JJ)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_setup_mpi_mod

  use sp_energy_mod
  use sp_energy_str_mod
  use sp_remd_mod
  use sp_rpath_mod
  use sp_boundary_mod
  use fileio_control_mod
  use math_libs_mod
  use messages_mod
  use mpi_parallel_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: setup_mpi_md
  public  :: setup_mpi_remd
  public  :: setup_mpi_rpath

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_mpi_md
  !> @brief        setup mpi in MD
  !> @authors      JJ, TM
  !! @param[in]    ene_info : ENERGY section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_mpi_md(ene_info)
  
    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info

    ! local variables
    integer                  :: i, alloc_stat, dealloc_stat
    integer                  :: omp_get_num_threads
    integer                  :: nthread
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

    end select

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
      !$omp parallel shared(nthread)
      nthread = omp_get_num_threads()
      !$omp end parallel
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

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_mpi_remd
  !> @brief        setup mpi in REMD
  !> @authors      TM, JJ
  !! @param[in]    ene_info   : ENERGY section control parameters information
  !! @param[in]    rep_info   : REMD section control parameters information
  !! @param[in]    bound_info : BOUNDARY section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_mpi_remd(ene_info, rep_info, bound_info)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_rep_info),        intent(in)    :: rep_info
    type(s_boundary_info),   intent(in)    :: bound_info

    ! local variables
    integer                  :: i, alloc_stat, dealloc_stat
    integer                  :: Nx, Ny, Nz, Bx, By, Bz
    integer                  :: rep_x, rep_y, rep_z
    integer                  :: domain_x, domain_y, domain_z
    integer                  :: proc_x, proc_y, proc_z
    integer                  :: omp_get_num_threads, nthread
    integer                  :: nreplicas, nproc_per_rep
    integer                  :: color, key, temp
    character(6)             :: mast, nonb1, nonb2
    logical                  :: prime_flag

    integer,     allocatable :: globrank(:), repno(:), reprank(:)
    logical,     allocatable :: repmast(:), nonbreal(:), nonbrecip(:)


#ifdef HAVE_MPI_GENESIS

    nreplicas  = product(rep_info%nreplicas(1:rep_info%dimension))

    ! Check relationship between nreplicas and num_procs
    !
    if (mod(nproc_world,nreplicas) /= 0) then
      call error_msg('Setup_Mpi_Remd> mod(total_MPI_procs,nreplicas)'//&
                     ' is not ZERO')
    else
      nproc_per_rep = int(nproc_world/nreplicas)
    end if

    ! check factorization of 2,3,5
    if (nproc_per_rep /= 1) then
      prime_flag=factorization_235(nproc_per_rep)
      if (.not. prime_flag)  &
        call error_msg('Setup_Mpi_Remd> number of MPI processes per replica'//&
                       ' should be multiples of 2,3,5')
    end if

    ! replica number in each dimension
    !
    rep_x = rep_info%nreplicas_x(1)
    rep_y = rep_info%nreplicas_y(1)
    rep_z = rep_info%nreplicas_z(1)

    ! domain number in each dimension
    !
    domain_x = bound_info%domain_x
    domain_y = bound_info%domain_y
    domain_z = bound_info%domain_z

    ! MPI number
    !
    if (rep_x > 0 .and. rep_y > 0 .and. rep_z > 0) then
      Nx = rep_x * domain_x
      Ny = rep_y * domain_y
      Nz = rep_z * domain_z
      proc_x = mod(my_world_rank,Nx)
      proc_y = mod(my_world_rank/Nx,Ny)
      proc_z = my_world_rank/(Nx*Ny)
      Bx = proc_x / domain_x 
      By = proc_y / domain_y 
      Bz = proc_z / domain_z 
    end if

    ! Communicator within country
    ! we make blocks of each replica only when three-dimensional indices exist
    !
    if (rep_x > 0 .and. rep_y > 0 .and. rep_z > 0) then
      my_country_no = Bx + By*rep_x + Bz*rep_x*rep_y
    else
      my_country_no = int(my_world_rank / nproc_per_rep)
    end if

    call mpi_comm_split(mpi_comm_world, my_country_no, my_world_rank, &
                        mpi_comm_country, ierror)
    call mpi_comm_size (mpi_comm_country, nproc_country,   ierror)
    call mpi_comm_rank (mpi_comm_country, my_country_rank, ierror)

    if (my_country_rank == 0) then
      replica_main_rank = .true.
    else
      replica_main_rank = .false.
    end if


    ! Communicator between replicas (mpi_comm_airplane)
    !
    temp  = my_country_no
    color = my_country_rank
    key   = temp

    call mpi_comm_split(mpi_comm_world, color, key, mpi_comm_airplane, ierror)
    call mpi_comm_size (mpi_comm_airplane, nproc_airplane,   ierror)
    call mpi_comm_rank (mpi_comm_airplane, my_airplane_rank, ierror)


    ! Equalize communicator (mpi_comm_city) between MD and REMD
    !
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

    end select

    ! Write the summary of setup MPI
    !
    allocate(globrank (nproc_world), &
             repno    (nproc_world), &
             reprank  (nproc_world), &
             repmast  (nproc_world), &
             nonbreal (nproc_world), &
             nonbrecip(nproc_world), &
             stat=alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    call mpi_gather(my_world_rank, 1, mpi_integer,                      &
                    globrank, 1, mpi_integer, 0, mpi_comm_world, ierror)
    call mpi_gather(my_country_no, 1, mpi_integer,                      &
                    repno, 1, mpi_integer, 0, mpi_comm_world, ierror)
    call mpi_gather(my_city_rank, 1, mpi_integer,                       &
                    reprank, 1, mpi_integer, 0, mpi_comm_world, ierror)
    call mpi_gather(replica_main_rank, 1, mpi_logical,                  &
                    repmast, 1, mpi_logical, 0, mpi_comm_world, ierror)
    call mpi_gather(real_calc, 1, mpi_logical,                          &
                    nonbreal, 1, mpi_logical, 0, mpi_comm_world, ierror)
    call mpi_gather(reciprocal_calc, 1, mpi_logical,                    &
                    nonbrecip, 1, mpi_logical, 0, mpi_comm_world, ierror)

    if (main_rank) then

      write(MsgOut,'(a)') 'Setup_Mpi_Remd> Summary of Setup MPI'
      write(MsgOut,'(a,i10)') '  number of MPI processes                = ', &
                              nproc_world
      write(MsgOut,'(a,i10)') '  number of MPI processes in one replica = ', &
                              nproc_country
#ifdef OMP

      !$omp parallel shared(nthread)
      nthread = omp_get_num_threads()
      !$omp end parallel
      write(MsgOut,'(a,i10)') '  number of OpenMP threads               = ', &
                              nthread
#else

      nthread = 1

#endif
      write(MsgOut,'(a,i10)') '  total number of CPU cores              = ', &
                              nproc_world*nthread
      write(MsgOut,'(a)') '     world_rank     country_no   country_rank'

      do i = 1, nproc_world

        mast  = '      '
        nonb1 = '      '
        nonb2 = '      '

        if (repmast(i))   mast  = 'MASTER'
        if (nonbreal(i))  nonb1 = 'REAL  '
        if (nonbrecip(i)) nonb2 = 'RECP  '

        write(MsgOut,'(3I15,2X,3A6)') &
              globrank(i), repno(i), reprank(i), NONB1, NONB2, mast
      end do

      write(MsgOut,'(a)') ''

    end if


    ! deallocate local array
    !
    deallocate(globrank, repno, reprank, repmast,      &
               nonbreal, nonbrecip, stat=dealloc_stat)
    if (dealloc_stat /= 0) &
      call error_msg_dealloc

#else

    call error_msg('Setup_Mpi_Remd> REMD is available with MPI')

#endif

    return

  end subroutine setup_mpi_remd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_mpi_rpath
  !> @brief        setup mpi in RPATH
  !> @authors      TM, JJ
  !! @param[in]    ene_info   : ENERGY section control parameters information
  !! @param[in]    rpath_info : RPATH section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_mpi_rpath(ene_info, rpath_info)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_rpath_info),      intent(in)    :: rpath_info

    ! local variables
    integer                  :: i, alloc_stat, dealloc_stat
    integer                  :: omp_get_num_threads, nthread
    integer                  :: nreplica, nproc_per_rep
    integer                  :: color, key, temp
    character(6)             :: mast, nonb1, nonb2
    logical                  :: prime_flag

    integer,     allocatable :: globrank(:), repno(:), reprank(:)
    logical,     allocatable :: repmast(:), nonbreal(:), nonbrecip(:)


#ifdef HAVE_MPI_GENESIS

    nreplica  = rpath_info%nreplica

    ! Check relationship between nreplicas and num_procs
    !
    if (mod(nproc_world,nreplica) /= 0) then
      call error_msg('Setup_Mpi_Rpath> mod(total_MPI_procs,nreplica)'//&
                     ' is not ZERO')
    else
      nproc_per_rep = int(nproc_world/nreplica)
    end if

    ! check factorization of 2,3,5
    if (nproc_per_rep /= 1) then
      prime_flag=factorization_235(nproc_per_rep)
      if (.not. prime_flag)  &
        call error_msg('Setup_Mpi_Rpath> number of MPI processes per replica'//&
                       ' should be multiples of 2,3,5')
    end if


    ! Communicator within country
    !
    my_country_no = int(my_world_rank / nproc_per_rep)

    call mpi_comm_split(mpi_comm_world, my_country_no, my_world_rank, &
                        mpi_comm_country, ierror)
    call mpi_comm_size (mpi_comm_country, nproc_country,   ierror)
    call mpi_comm_rank (mpi_comm_country, my_country_rank, ierror)

    if (my_country_rank == 0) then
      replica_main_rank = .true.
    else
      replica_main_rank = .false.
    end if


    ! Communicator between replicas (mpi_comm_airplane)
    !
    temp  = my_country_no
    color = my_country_rank
    key   = temp

    call mpi_comm_split(mpi_comm_world, color, key, mpi_comm_airplane, ierror)
    call mpi_comm_size (mpi_comm_airplane, nproc_airplane,   ierror)
    call mpi_comm_rank (mpi_comm_airplane, my_airplane_rank, ierror)


    ! Equalize communicator (mpi_comm_city) between MD and RPATH
    !
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

    end select

    ! Write the summary of setup MPI
    !
    allocate(globrank (nproc_world), &
             repno    (nproc_world), &
             reprank  (nproc_world), &
             repmast  (nproc_world), &
             nonbreal (nproc_world), &
             nonbrecip(nproc_world), &
             stat=alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    call mpi_gather(my_world_rank, 1, mpi_integer,                      &
                    globrank, 1, mpi_integer, 0, mpi_comm_world, ierror)
    call mpi_gather(my_country_no, 1, mpi_integer,                      &
                    repno, 1, mpi_integer, 0, mpi_comm_world, ierror)
    call mpi_gather(my_city_rank, 1, mpi_integer,                       &
                    reprank, 1, mpi_integer, 0, mpi_comm_world, ierror)
    call mpi_gather(replica_main_rank, 1, mpi_logical,                  &
                    repmast, 1, mpi_logical, 0, mpi_comm_world, ierror)
    call mpi_gather(real_calc, 1, mpi_logical,                          &
                    nonbreal, 1, mpi_logical, 0, mpi_comm_world, ierror)
    call mpi_gather(reciprocal_calc, 1, mpi_logical,                    &
                    nonbrecip, 1, mpi_logical, 0, mpi_comm_world, ierror)

    if (main_rank) then

      write(MsgOut,'(a)') 'Setup_Mpi_Rpath> Summary of Setup MPI'
      write(MsgOut,'(a,i10)') '  number of MPI processes                = ', &
                              nproc_world
      write(MsgOut,'(a,i10)') '  number of MPI processes in one replica = ', &
                              nproc_country
#ifdef OMP

      !$omp parallel shared(nthread)
      nthread = omp_get_num_threads()
      !$omp end parallel
      write(MsgOut,'(a,i10)') '  number of OpenMP threads               = ', &
                              nthread
#else

      nthread = 1

#endif
      write(MsgOut,'(a,i10)') '  total number of CPU cores              = ', &
                              nproc_world*nthread
      write(MsgOut,'(a)') '     world_rank     country_no   country_rank'

      do i = 1, nproc_world

        mast  = '      '
        nonb1 = '      '
        nonb2 = '      '

        if (repmast(i))   mast  = 'MASTER'
        if (nonbreal(i))  nonb1 = 'REAL  '
        if (nonbrecip(i)) nonb2 = 'RECP  '

        write(MsgOut,'(3I15,2X,3A6)') &
              globrank(i), repno(i), reprank(i), NONB1, NONB2, mast
      end do

      write(MsgOut,'(a)') ''

    end if


    ! deallocate local array
    !
    deallocate(globrank, repno, reprank, repmast,      &
               nonbreal, nonbrecip, stat=dealloc_stat)
    if (dealloc_stat /= 0) &
      call error_msg_dealloc

#else

    call error_msg('Setup_Mpi_Rpath> RPATH is available with MPI')

#endif

    return

  end subroutine setup_mpi_rpath

end module sp_setup_mpi_mod
