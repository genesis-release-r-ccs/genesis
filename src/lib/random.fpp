!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   random_mod
!> @brief   random number generator
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module random_mod

  use messages_mod
  use constants_mod
  use mpi_parallel_mod
#ifdef MPI
  use mpi
#endif

  implicit none
  private

  ! structures
  logical               :: g_in_use       = .false.
  logical               :: g_stock_pushed = .true.
  character,    pointer :: g_dsfmt(:)
  character,    pointer :: g_dsfmt_stock(:)
  logical               :: legacy_in_use  = .false.

  ! subroutines
  public  :: random_init
  public  :: random_term
!  public  :: random_seed
  public  :: random_get
  public  :: random_get_gauss
  public  :: random_get_legacy
  public  :: random_push_stock
  public  :: random_pull_stock
  public  :: random_stock_tobyte
  public  :: random_stock_frombyte
  public  :: random_stock_mpi_tobyte
  public  :: random_stock_mpi_frombyte
  public  :: random_seed_initial_time

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    random_init
  !> @brief        initialize random number system
  !! @authors      NT
  !! @param[in]    iseed  : seed number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine random_init(iseed)

    ! formal arguments
    integer,                 intent(in)    :: iseed

    ! local variables
    integer                  :: dsfmt_size


    if (g_in_use) &
      return

    call get_size_of_dsfmt_t(dsfmt_size)
    allocate(g_dsfmt(dsfmt_size), g_dsfmt_stock(dsfmt_size))

    call dsfmt_init_gen_rand(g_dsfmt, iseed)

    g_in_use = .true.
    g_stock_pushed = .false.

    if (main_rank) then
      write(MsgOut,'(a)')    'Random_Init> Initialize the random number'
      write(MsgOut,'(a,i0,i0)') '  seed            = ', iseed
      write(MsgOut,'(a)')    ''
    end if

    return

  end subroutine random_init

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    random_term
  !> @brief        terminate random number system
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine random_term


    if (.not. g_in_use) &
      return

    deallocate(g_dsfmt_stock, g_dsfmt)

    g_in_use = .false.

    return

  end subroutine random_term

!  !======1=========2=========3=========4=========5=========6=========7=========8
!  !
!  !  Subroutine    random_seed
!  !> @brief        setup seed number
!  !! @authors      NT
!  !! @param[in]    iseed  : seed number
!  !
!  !======1=========2=========3=========4=========5=========6=========7=========8
!
!  subroutine random_seed(iseed)
!
!    ! formal arguments
!    integer,                 intent(in)    :: iseed
!
!
!    if (.not. g_in_use) &
!      call error_msg('Random_Seed> random system is not initialized.')
!
!    call dsfmt_init_gen_rand(g_dsfmt, iseed)
!
!    return
!
!  end subroutine random_seed

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      random_get
  !> @brief        get random number
  !! @authors      NT
  !! @return       random number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function random_get()

    ! return
    real(dp)                 :: random_get


    call dsfmt_genrand_close0_open1(g_dsfmt, random_get)

    return

  end function random_get

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      random_get_gauss
  !> @brief        get random number
  !! @authors      NT
  !! @return       random number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function random_get_gauss()

    ! return
    real(dp)                 :: random_get_gauss

    ! local variables
    real(dp)                 :: v1, v2, rsq


    rsq = 2.0_dp
    do while (rsq >= 1.0_dp)
      call dsfmt_genrand_close0_open1(g_dsfmt, v1)
      call dsfmt_genrand_close0_open1(g_dsfmt, v2)
      v1  = 2.0_dp*v1 - 1.0_dp
      v2  = 2.0_dp*v2 - 1.0_dp
      rsq = v1*v1 + v2*v2
    end do
    rsq   = sqrt(-2.0_dp * log(rsq) / rsq)
    random_get_gauss = rsq * v1

    return

  end function random_get_gauss

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      random_get_legacy
  !> @brief        random number generator (uniform distribution [0,1])
  !! @authors      Takaharu Mori (TM)
  !! @param[inout] iseed : random number seed
  !! @note         P.A.W.Lewis et al., BM Systems Journal, 8, 136 (1969).
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function random_get_legacy(iseed)

    ! formal arguments
    integer,                 intent(inout) :: iseed

    ! local variables
    real(dp)                 :: dseed
    real(dp)                 :: random_get_legacy


    if (iseed <= 1) iseed = 314159

    if (main_rank .and. .not. legacy_in_use) then
      write(MsgOut,'(a)') 'Random_Init_Velocity> Initialize the random number'
      write(MsgOut,'(a,i0,i0)') '  seed            = ', iseed
      write(MsgOut,'(a)')    ''
    end if

    dseed = 16807.0_dp * iseed
    dseed = mod(dseed, 2147483647.0_dp)

    random_get_legacy = dseed / 2147483711.0_dp
    iseed = dseed

    legacy_in_use = .true.

    return

  end function random_get_legacy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    random_push_stock
  !> @brief        push the random internal stats to stocked internal stats
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine random_push_stock

    ! local variables
    integer                  :: dsfmt_size


    call get_size_of_dsfmt_t(dsfmt_size)

    g_dsfmt_stock(1:dsfmt_size) = g_dsfmt(1:dsfmt_size)
    g_stock_pushed = .true.

    return

  end subroutine random_push_stock

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    random_pull_stock
  !> @brief        pull the random internal stats from stocked internal stats
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine random_pull_stock

    ! local variables
    integer                  :: dsfmt_size


    call get_size_of_dsfmt_t(dsfmt_size)

    g_dsfmt(1:dsfmt_size) = g_dsfmt_stock(1:dsfmt_size)

    return

  end subroutine random_pull_stock

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    random_stock_tobyte
  !> @brief        get the bytes from random stocked internal stats
  !! @authors      NT
  !! @param[inout] bytes : serialized bytes
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine random_stock_tobyte(bytes)

    ! formal arguments
    character, allocatable,  intent(inout) :: bytes(:)

    ! local variables
    integer                  :: dsfmt_size


    call get_size_of_dsfmt_t(dsfmt_size)

    if (allocated(bytes)) then
      if (size(bytes) /= dsfmt_size) then
        deallocate(bytes)
        allocate(bytes(dsfmt_size))
      end if
    else
      allocate(bytes(dsfmt_size))
    end if

    if (g_stock_pushed) then
      bytes(1:dsfmt_size) = g_dsfmt_stock(1:dsfmt_size)
    else
      bytes(1:dsfmt_size) = g_dsfmt(1:dsfmt_size)
    end if

    return

  end subroutine random_stock_tobyte

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    random_stock_frombyte
  !> @brief        put the bytes to random stocked internal stats
  !! @authors      NT
  !! @param[in]    bytes : serialized bytes
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine random_stock_frombyte(bytes)

    ! formal arguments
    character, allocatable,  intent(in)    :: bytes(:)

    ! local variables
    integer                  :: dsfmt_size


    call get_size_of_dsfmt_t(dsfmt_size)

    if (size(bytes) /= dsfmt_size) &
      call error_msg( &
      'Random_Frombyte> different from bytes size and random internal size.')
    
    g_dsfmt_stock(1:dsfmt_size) = bytes(1:dsfmt_size)

    return

  end subroutine random_stock_frombyte

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    random_stock_mpi_tobyte
  !> @brief        get the bytes from random stocked internal stats 
  !!               per MPI communicator
  !! @authors      NT
  !! @param[in]    mpi_comm : MPI communicator
  !! @param[inout] bytes    : serialized bytes
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine random_stock_mpi_tobyte(mpi_comm, bytes)

    ! formal arguments
    integer,                 intent(in)    :: mpi_comm
    character, allocatable,  intent(inout) :: bytes(:)

    ! local variables
    integer                  :: comm_size, comm_rank, ierr
    integer                  :: dsfmt_size, byte_size
    integer                  :: unit_no


    if (.not. g_in_use) &
      call error_msg('Random_Mpi_Tobyte> random system is not initialized.')


    ! gather all rank random internal stats
    !

    call get_size_of_dsfmt_t(dsfmt_size)

#ifdef MPI

    call MPI_Comm_size(mpi_comm, comm_size, ierr)
    call MPI_Comm_rank(mpi_comm, comm_rank, ierr)

    byte_size = dsfmt_size * comm_size

    if (allocated(bytes)) then
      if (size(bytes) /= byte_size) then
        deallocate(bytes)
        allocate(bytes(byte_size))
      end if
    else
      allocate(bytes(byte_size))
    end if

    if (g_stock_pushed) then
      call MPI_Gather(g_dsfmt_stock, dsfmt_size, MPI_BYTE, &
                      bytes,         dsfmt_size, MPI_BYTE, &
                      0, mpi_comm, ierr)
    else
      call MPI_Gather(g_dsfmt      , dsfmt_size, MPI_BYTE, &
                      bytes,         dsfmt_size, MPI_BYTE, &
                      0, mpi_comm, ierr)
    end if

    if (comm_rank /= 0) &
      deallocate(bytes)

#else

    comm_size = 1
    comm_rank = 0

    byte_size = dsfmt_size * comm_size

    if (allocated(bytes)) then
      if (size(bytes) /= byte_size) then
        deallocate(bytes)
        allocate(bytes(byte_size))
      end if
    else
      allocate(bytes(byte_size))
    end if

    if (g_stock_pushed) then
      bytes(1:byte_size) = g_dsfmt_stock(1:dsfmt_size)
    else
      bytes(1:byte_size) = g_dsfmt(1:dsfmt_size)
    end if

#endif

    return

  end subroutine random_stock_mpi_tobyte

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    random_stock_mpi_frombyte
  !> @brief        put the bytes to random stocked internal stats 
  !!               per MPI communicator
  !! @authors      NT
  !! @param[in]    mpi_comm : MPI communicator
  !! @param[in]    bytes    : bytes array
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine random_stock_mpi_frombyte(mpi_comm, bytes)

    ! formal arguments
    integer,                 intent(in)    :: mpi_comm
    character,               intent(in)    :: bytes(:)

    ! local variables
    integer                  :: comm_size, comm_size_f, comm_rank, ierr
    integer                  :: dsfmt_size, byte_size
    integer                  :: unit_no


    if (.not. g_in_use) &
      call error_msg('Random_Mpi_Frombyte> random system is not initialized.')


    call get_size_of_dsfmt_t(dsfmt_size)

#ifdef MPI

    call MPI_Comm_size(mpi_comm, comm_size, ierr)
    call MPI_Comm_rank(mpi_comm, comm_rank, ierr)

#else

    comm_size = 1
    comm_rank = 0

#endif

    byte_size = dsfmt_size * comm_size
    if (size(bytes) /= byte_size) then

      if (main_rank) then
        write(MsgOut,'(a)') &
     'Random_Mpi_Frombyte> MPI communicator size was changed from previous run.'
        write(MsgOut,'(a)') &
     '                     skip read random internal state of restart file.'
      end if

      g_dsfmt_stock(:) = g_dsfmt(:)
      return
    end if

    ! scatt all rank random stats
    !

#ifdef MPI

    call MPI_Scatter(bytes,         dsfmt_size, MPI_BYTE, &
                     g_dsfmt_stock, dsfmt_size, MPI_BYTE, &
                     0, mpi_comm, ierr)

#else

    g_dsfmt_stock(1:dsfmt_size) = bytes(1:byte_size)

#endif

    return

  end subroutine random_stock_mpi_frombyte

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    random_seed_initial_time
  !> @brief        put iseed by time function
  !! @authors      CK
  !! @param[out]   gseed    : out seed from time function
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine random_seed_initial_time(gseed)
    integer,                 intent(out)    :: gseed
    integer                  :: tval(8)
    integer(8)               :: yy, mm, dd, hh, nn, ss, ms

    call date_and_time(values=tval)
    yy  = int(tval(1),kind=8)
    mm  = int(tval(2),kind=8)
    dd  = int(tval(3),kind=8)
    hh  = int(tval(5),kind=8)
    nn  = int(tval(6),kind=8)
    ss  = int(tval(7),kind=8)
    ms  = int(tval(8),kind=8)
    gseed = ms*1000+ss*100+nn*10+hh

    return

  end subroutine random_seed_initial_time

end module random_mod
