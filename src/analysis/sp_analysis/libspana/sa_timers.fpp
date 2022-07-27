!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   timers_mod
!> @brief   process timer (msec)
!! @authors Yuji Sugita (YS), Takaharu Mori (TM), Naoyuki Miyashita (NM), 
!!          Norio Takase (NT), Kiyoshi Yagi (KY)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module sa_timers_mod

  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! parameters
  integer, public, parameter :: TimerOff        = 0
  integer, public, parameter :: TimerOn         = 1

  integer, public, parameter :: TimerTotal      = 1
  integer, public, parameter :: TimerAnalysis   = 2 

  integer,         parameter :: InvalidID       = -1
  integer,         parameter :: NumTimers = 2  !< total number of timers
  integer,         parameter :: MaxProc   = 100 !< maximum number of processes

  ! variables
  real(dp),        public    :: total_time(NumTimers)
  real(dp),        save      :: proc_time_st(NumTimers)

  ! functions
  public  :: timer
  public  :: sa_output_time
  private :: get_unix_time
  private :: days

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    timer
  !> @brief        integrate process time
  !! @authors      TM
  !! @param[in]    mode : timer index
  !! @param[in]    set  : 0 = turn off timer; 1 = turn on timer
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine timer(mode, set)

    ! formal arguments
    integer,                 intent(in) :: mode
    integer,                 intent(in) :: set


    select case (set)

    case (TimerOn)

      proc_time_st(mode) = get_unix_time()

    case (TimerOff)

      total_time(mode) = total_time(mode) &
                       + (get_unix_time()-proc_time_st(mode))*1.0d-6

    end select

    return

  end subroutine timer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_time
  !> @brief        output timer profile
  !! @authors      TM
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine sa_output_time

    ! local variables
    real(dp)                      :: tottime_ebond, tottime_enbond
    real(dp)                      :: tottime_dynamc, avetime_setup
    integer                       :: i, j, alloc_stat, dealloc_stat
    integer                       :: num_realnodes, num_recipnodes

    real(dp)        , allocatable :: tottime(:,:), avetime(:), sumtime(:)
    real(dp)        , allocatable :: maxtime(:), mintime(:)
    logical,          allocatable :: node_real(:), node_recip(:)


    ! allocate local array
    !
    allocate(node_real (nproc_world),            &
             node_recip(nproc_world),            &
             tottime   (NumTimers, nproc_world), &
             sumtime   (NumTimers),              &
             avetime   (NumTimers),              &
             mintime   (NumTimers),              &
             maxtime   (NumTimers),              &
             stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

#ifdef HAVE_MPI_GENESIS
    call mpi_allgather(total_time, NumTimers, mpi_double_precision, &
                       tottime,    NumTimers, mpi_double_precision, &
                       mpi_comm_world, ierror)
    call mpi_allgather(real_calc, 1, mpi_logical,                   &
                       node_real, 1, mpi_logical,                   &
                       mpi_comm_world, ierror)
    call mpi_allgather(reciprocal_calc, 1, mpi_logical,             &
                       node_recip,      1, mpi_logical,             &
                       mpi_comm_world, ierror)
#else
    tottime(1:NumTimers,1) = total_time(1:NumTimers)
    node_real(1)  = real_calc
    node_recip(1) = reciprocal_calc
#endif

    ! count the number of real and resiprocal part nodes
    !
    if (main_rank) then
      num_realnodes  = 0
      num_recipnodes = 0
      do i = 1, nproc_world
        if (node_real(i))  num_realnodes  = num_realnodes  + 1
        if (node_recip(i)) num_recipnodes = num_recipnodes + 1
      end do

    ! calculate total time over all nodes
    !
      sumtime(1:NumTimers) = 0.0_dp
      do i = 1, NumTimers
        maxtime(i) = tottime(i,1)
        mintime(i) = tottime(i,1)
        do j = 1, nproc_world
          sumtime(i) = sumtime(i) + tottime(i,j)
          if (maxtime(i) < tottime(i,j)) then
            maxtime(i) = tottime(i,j)
          endif
          if (mintime(i) > tottime(i,j)) then
            mintime(i) = tottime(i,j)
          endif
        end do
      end do

      ! average over all ranks
      !
      avetime(TimerTotal)      = sumtime(TimerTotal)      / nproc_world
      avetime(TimerAnalysis)   = sumtime(TimerAnalysis)   / nproc_world

      ! write detailed timer profile of the main_rank
      !
!    if (main_rank) then
      avetime_setup   = avetime(TimerTotal) - avetime(TimerAnalysis)

      write(MsgOut,'(a)') 'spana_Output_Time> Averaged timer profile'
      write(MsgOut,'(a,f12.3)') '  total time      =', avetime(TimerTotal)
      write(MsgOut,'(a,f12.3)') '    setup         =', avetime_setup
      write(MsgOut,'(a,f12.3)') '    analysis      =', avetime(TimerAnalysis)
      write(MsgOut,'(a)') ''

    end if


    ! deallocate local array
    !
    deallocate(node_real, node_recip, tottime, sumtime, avetime, &
               mintime, maxtime, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine sa_output_time

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      get_unix_time
  !> @brief        get the local unix time (not GMT)
  !! @authors      JJ, ST(FJ), NT
  !! @return       milliseconds
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_unix_time()

    ! return value
    real(dp)               :: get_unix_time

    ! parameters
    real(dp)               :: usec

#ifdef KCOMP
    call gettod(usec)
#elif __INTEL_COMPILER
    call clockx(usec)
#else
    integer,   save :: init_flg = 0
    integer*8 :: CountPerSec, CountMax
    real(dp), save :: ratio_usec
    integer*8 time_c

    if( init_flg == 0 ) then
      init_flg = 1
      call system_clock(time_c, CountPerSec, CountMax)
      ratio_usec = 1.0d+6 / dble(CountPerSec)
    else
      call system_clock(time_c)
    endif

    usec = dble(time_c)*ratio_usec
#endif

     get_unix_time = usec

    return

  end function get_unix_time

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      days
  !> @brief        get number of days from AD.1
  !! @authors      NT
  !! @param[in]    year : A.D. year 
  !! @return       number of days
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function days(year)

    ! return values
    integer(8)               :: days

    ! formal arguments
    integer(8),              intent(in)    :: year


    days = year * 365
    days = days + year / 4
    days = days - year / 100
    days = days + year / 400

    return

  end function days

end module sa_timers_mod
