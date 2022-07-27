!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_setup_mpi_mod
!> @brief   setup MPI
!! @authors Jaewoon Jung (JJ), Takaharu Mori (TM)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_setup_mpi_mod

  use at_remd_mod
  use at_rpath_mod
  use at_vibration_mod
  use at_energy_mod
  use at_energy_str_mod
  use fileio_control_mod
  use math_libs_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: setup_mpi_md
  public  :: setup_mpi_remd
  public  :: setup_mpi_rpath
  public  :: setup_mpi_vib

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
    real(wp)                 :: ratio, recip
    integer                  :: i, alloc_stat, dealloc_stat
    integer                  :: omp_get_num_threads
    integer                  :: color, key, nthread
    integer                  :: node_recip, node_real, expon
    logical                  :: prime_flag

    integer,     allocatable :: globrank(:), locrank(:)
    logical,     allocatable :: nonbreal(:), nonbrecip(:)


#ifdef HAVE_MPI_GENESIS
    ! Equalize communicators (comm_country & comm_city) between MD
    !
    mpi_comm_country = mpi_comm_world
    nproc_country    = nproc_world
    my_country_rank  = my_world_rank
    replica_main_rank = main_rank

    ! check factorization of 2,3,5
    if (nproc_country /= 1) then
      prime_flag=factorization_235(nproc_country)
      if (.not. prime_flag)  &
        call error_msg('Setup_Mpi_Md> number of MPI processes should be'//&
                       ' multiples of 2,3,5')
    end if
#else
    nproc_country    = 1
    my_country_rank  = 0
    replica_main_rank = .true.
#endif

    select case (ene_info%electrostatic)

    case (ElectrostaticPME)

#ifdef HAVE_MPI_GENESIS

      if (ene_info%pme_multiple) then

        if (nproc_country == 1) &
          call error_msg('Setup_Mpi_Md> "PME_multiple" must be "NO", '//&
                         'if only 1 MPI is used for one replica.')

        ratio = ene_info%pme_mul_ratio
        recip = real(nproc_country,wp) / (1.0_wp + ratio)
        expon = int(log(recip)/log(real(2,wp)))
        node_recip = 2**expon
        node_real  = nproc_country - node_recip

        if (my_country_rank < node_real) then
          real_calc = .true.
          reciprocal_calc = .false.
          color = 1
          key   = my_country_rank
        else
          real_calc = .false.
          reciprocal_calc = .true.
          color = 2
          key   = my_country_rank
        end if

        call mpi_comm_split(mpi_comm_country, color, key, mpi_comm_city, ierror)
        call mpi_comm_size (mpi_comm_city, nproc_city,   ierror)
        call mpi_comm_rank (mpi_comm_city, my_city_rank, ierror)

      else

        real_calc = .true.
        reciprocal_calc = .true.
        mpi_comm_city = mpi_comm_country
        nproc_city    = nproc_country
        my_city_rank  = my_country_rank

      end if
#else
        real_calc = .true.
        reciprocal_calc = .true.
        nproc_city   = nproc_country
        my_city_rank = my_country_rank
#endif

    case (ElectrostaticCUTOFF)

#ifdef HAVE_MPI_GENESIS
      real_calc = .true.
      reciprocal_calc = .false.
      mpi_comm_city = mpi_comm_country
      nproc_city    = nproc_country
      my_city_rank  = my_country_rank
#else
      real_calc = .true.
      reciprocal_calc = .false.
      nproc_city   = nproc_country
      my_city_rank = my_country_rank
#endif

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

#ifdef OMP
    !$omp parallel shared(nthread)
    nthread = omp_get_num_threads()
    !$omp end parallel
    write(MsgOut,'(a,i10)') '  number of OpenMP threads  = ', nthread
#else
    nthread = 1
#endif

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
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[in]    rep_info : REMD section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_mpi_remd(ene_info, rep_info)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_rep_info),        intent(in)    :: rep_info

    ! local variables
    real(wp)                 :: ratio, recip
    integer                  :: i, alloc_stat, dealloc_stat
    integer                  :: omp_get_num_threads, nthread
    integer                  :: nreplicas, nproc_per_rep
    integer                  :: color, key, temp
    integer                  :: node_recip, node_real, expon
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


    ! Equalize communicator (mpi_comm_city) between MD and REMD
    !
    select case (ene_info%electrostatic)

    case (ElectrostaticPME)

      if (ene_info%pme_multiple) then

        if (nproc_country == 1) &
          call error_msg('Setup_Mpi_Remd> "PME_multiple" must be "NO", '//&
                         'if only 1 MPI is used for one replica.')

        ratio = ene_info%pme_mul_ratio
        recip = real(nproc_country,wp) / (1.0_wp + ratio)
        expon = int(log(recip)/log(real(2,wp)))
        node_recip = 2**expon
        node_real  = nproc_country - node_recip

        if (my_country_rank < node_real) then
          real_calc = .true.
          reciprocal_calc = .false.
          color = 1
          key   = my_country_rank
        else
          real_calc = .false.
          reciprocal_calc = .true.
          color = 2
          key   = my_country_rank
        end if

        call mpi_comm_split(mpi_comm_country, color, key, mpi_comm_city, ierror)
        call mpi_comm_size (mpi_comm_city, nproc_city,   ierror)
        call mpi_comm_rank (mpi_comm_city, my_city_rank, ierror)

      else

        real_calc = .true.
        reciprocal_calc = .true.
        mpi_comm_city = mpi_comm_country
        nproc_city    = nproc_country
        my_city_rank  = my_country_rank

      end if

    case (ElectrostaticCUTOFF)

      real_calc = .true.
      reciprocal_calc = .false.
      mpi_comm_city = mpi_comm_country
      nproc_city    = nproc_country
      my_city_rank  = my_country_rank

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
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[in]    rpath_info : RPATH section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_mpi_rpath(ene_info, rpath_info)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_rpath_info),      intent(in)    :: rpath_info

    ! local variables
    real(wp)                 :: ratio, recip
    integer                  :: i, alloc_stat, dealloc_stat
    integer                  :: omp_get_num_threads, nthread
    integer                  :: nreplica, nproc_per_rep
    integer                  :: color, key, temp
    integer                  :: node_recip, node_real, expon
    character(6)             :: mast, nonb1, nonb2
    logical                  :: prime_flag

    integer,     allocatable :: globrank(:), repno(:), reprank(:)
    logical,     allocatable :: repmast(:), nonbreal(:), nonbrecip(:)

#ifdef HAVE_MPI_GENESIS

    nreplica  = rpath_info%nreplica

    ! Check relationship between nreplicas and num_procs
    !
    if (nreplica <= nproc_world) then
      if (mod(nproc_world,nreplica) /= 0)  &
        call error_msg('Setup_Mpi_Rpath> mod(total_MPI_procs,nreplica)'//&
                       ' is not ZERO')
      nproc_per_rep = int(nproc_world/nreplica)
    else
      if ((mod(nreplica,nproc_world) /= 0) .or. &
        (rpath_info%rpathmode == 1) .or. (rpath_info%rpathmode == 3)) &
        call error_msg('Setup_Mpi_Rpath> mod(nreplica,total_MPI_procs)'//&
                       ' is not ZERO')
      nrep_per_proc = int(nreplica/nproc_world)

      nproc_per_rep = 1
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

    ! Equalize communicator (mpi_comm_city) between MD and REMD/RPATH
    !
    select case (ene_info%electrostatic)

    case (ElectrostaticPME)

      if (ene_info%pme_multiple) then

        if (nproc_country == 1) &
          call error_msg('Setup_Mpi_Rpath> "PME_multiple" must be "NO", '//&
                         'if only 1 MPI is used for one replica.')

        ratio = ene_info%pme_mul_ratio
        recip = real(nproc_country,wp) / (1.0_wp + ratio)
        expon = int(log(recip)/log(real(2,wp)))
        node_recip = 2**expon
        node_real  = nproc_country - node_recip

        if (my_country_rank < node_real) then
          real_calc = .true.
          reciprocal_calc = .false.
          color = 1
          key   = my_country_rank
        else
          real_calc = .false.
          reciprocal_calc = .true.
          color = 2
          key   = my_country_rank
        end if

        call mpi_comm_split(mpi_comm_country, color, key, mpi_comm_city, ierror)
        call mpi_comm_size (mpi_comm_city, nproc_city,   ierror)
        call mpi_comm_rank (mpi_comm_city, my_city_rank, ierror)

      else

        real_calc = .true.
        reciprocal_calc = .true.
        mpi_comm_city = mpi_comm_country
        nproc_city    = nproc_country
        my_city_rank  = my_country_rank

      end if

    case (ElectrostaticCUTOFF)

      real_calc = .true.
      reciprocal_calc = .false.
      mpi_comm_city = mpi_comm_country
      nproc_city    = nproc_country
      my_city_rank  = my_country_rank

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

    replica_main_rank = .true.

    select case (ene_info%electrostatic)

    case (ElectrostaticPME)

      real_calc = .true.
      reciprocal_calc = .true.
      nproc_city   = nproc_country
      my_city_rank = my_country_rank

    case (ElectrostaticCUTOFF)

      real_calc = .true.
      reciprocal_calc = .false.
      nproc_city   = nproc_country
      my_city_rank = my_country_rank

    end select

    replica_main_rank = .true.

    select case (ene_info%electrostatic)

    case (ElectrostaticPME)

      real_calc = .true.
      reciprocal_calc = .true.
      nproc_city   = nproc_country
      my_city_rank = my_country_rank

    case (ElectrostaticCUTOFF)

      real_calc = .true.
      reciprocal_calc = .false.
      nproc_city   = nproc_country
      my_city_rank = my_country_rank

    end select

#ifdef OMP

    !$omp parallel shared(nthread)
    nthread = omp_get_num_threads()
    !$omp end parallel
    write(MsgOut,'(a,i10)') '  number of OpenMP threads  = ', &
                            nthread
#else

    nthread = 1

#endif

    write(MsgOut,'(a)') 'Setup_Mpi_Rpath> Single CPU is used.'
    write(MsgOut,'(a)') ''

#endif

    return

  end subroutine setup_mpi_rpath

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_mpi_vib
  !> @brief        setup mpi in VIB
  !> @authors      KY
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[in]    vib_info : VIBRATION section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_mpi_vib(ene_info, vib_info)

    ! formal arguments
    type(s_ene_info),      intent(in)    :: ene_info
    type(s_vib_info),      intent(in)    :: vib_info

    ! local variables
    real(wp)                 :: ratio, recip
    integer                  :: i, alloc_stat, dealloc_stat
    integer                  :: omp_get_num_threads, nthread
    integer                  :: nreplica, nproc_per_rep
    integer                  :: color, key, temp
    integer                  :: node_recip, node_real, expon
    character(6)             :: mast, nonb1, nonb2
    logical                  :: prime_flag

    integer,     allocatable :: globrank(:), repno(:), reprank(:)
    logical,     allocatable :: repmast(:), nonbreal(:), nonbrecip(:)

#ifdef HAVE_MPI_GENESIS

    nreplica  = vib_info%nreplica

    ! Check relationship between nreplicas and num_procs
    !
    if (mod(nproc_world,nreplica) /= 0)  then
      write(MsgOut,'("Setup_Mpi_Vib> nproc_world = ",i5)') nproc_world
      write(MsgOut,'("Setup_Mpi_Vib> nreplica    = ",i5)') nreplica
      call error_msg('Setup_Mpi_Vib> mod(total_MPI_procs,nreplica)'//&
                     ' is not ZERO')
    end if
    nproc_per_rep = int(nproc_world/nreplica)

    ! check factorization of 2,3,5
    if (nproc_per_rep /= 1) then
      prime_flag=factorization_235(nproc_per_rep)
      if (.not. prime_flag)  &
        call error_msg('Setup_Mpi_Vib> number of MPI processes per replica'//&
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

    ! Equalize communicator (mpi_comm_city) between MD and REMD/RPATH
    !
    select case (ene_info%electrostatic)

    case (ElectrostaticPME)

      if (ene_info%pme_multiple) then

        if (nproc_country == 1) &
          call error_msg('Setup_Mpi_Rpath> "PME_multiple" must be "NO", '//&
                         'if only 1 MPI is used for one replica.')

        ratio = ene_info%pme_mul_ratio
        recip = real(nproc_country,wp) / (1.0_wp + ratio)
        expon = int(log(recip)/log(real(2,wp)))
        node_recip = 2**expon
        node_real  = nproc_country - node_recip

        if (my_country_rank < node_real) then
          real_calc = .true.
          reciprocal_calc = .false.
          color = 1
          key   = my_country_rank
        else
          real_calc = .false.
          reciprocal_calc = .true.
          color = 2
          key   = my_country_rank
        end if

        call mpi_comm_split(mpi_comm_country, color, key, mpi_comm_city, ierror)
        call mpi_comm_size (mpi_comm_city, nproc_city,   ierror)
        call mpi_comm_rank (mpi_comm_city, my_city_rank, ierror)

      else

        real_calc = .true.
        reciprocal_calc = .true.
        mpi_comm_city = mpi_comm_country
        nproc_city    = nproc_country
        my_city_rank  = my_country_rank

      end if

    case (ElectrostaticCUTOFF)

      real_calc = .true.
      reciprocal_calc = .false.
      mpi_comm_city = mpi_comm_country
      nproc_city    = nproc_country
      my_city_rank  = my_country_rank

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

      write(MsgOut,'(a)') 'Setup_Mpi_Vib> Summary of Setup MPI'
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

    replica_main_rank = .true.

    select case (ene_info%electrostatic)

    case (ElectrostaticPME)

      real_calc = .true.
      reciprocal_calc = .true.
      nproc_city   = nproc_country
      my_city_rank = my_country_rank

    case (ElectrostaticCUTOFF)

      real_calc = .true.
      reciprocal_calc = .false.
      nproc_city   = nproc_country
      my_city_rank = my_country_rank

    end select

    replica_main_rank = .true.

    select case (ene_info%electrostatic)

    case (ElectrostaticPME)

      real_calc = .true.
      reciprocal_calc = .true.
      nproc_city   = nproc_country
      my_city_rank = my_country_rank

    case (ElectrostaticCUTOFF)

      real_calc = .true.
      reciprocal_calc = .false.
      nproc_city   = nproc_country
      my_city_rank = my_country_rank

    end select

#ifdef OMP

    !$omp parallel shared(nthread)
    nthread = omp_get_num_threads()
    !$omp end parallel
    write(MsgOut,'(a,i10)') '  number of OpenMP threads  = ', &
                            nthread
#else

    nthread = 1

#endif

    write(MsgOut,'(a)') 'Setup_Mpi_Rpath> Single CPU is used.'
    write(MsgOut,'(a)') ''

#endif

    return

  end subroutine setup_mpi_vib

end module at_setup_mpi_mod
