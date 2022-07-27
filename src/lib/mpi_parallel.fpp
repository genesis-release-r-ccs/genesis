!--------1---------2---------3---------4---------5---------6---------7---------8
! 
!  Module   mpi_parallel
!> @brief   utilities for mpi related programs                 
!! @authors Jaewoon Jung (JJ), Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
! 
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module mpi_parallel_mod

#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! variables for MPI in MD
  integer, public :: mpi_comm_country  = 0
  integer, public :: mpi_comm_city     = 0
  integer, public :: nproc_world       = 1
  integer, public :: nproc_country     = 1
  integer, public :: nproc_city        = 1
  integer, public :: nprocx            = 1
  integer, public :: nprocy            = 1
  integer, public :: nprocz            = 1
  integer, public :: my_world_rank     = 0
  integer, public :: my_country_rank   = 0
  integer, public :: my_rank_pio       = 0
  integer, public :: my_city_rank      = 0
  integer, public :: my_x_rank         = 0
  integer, public :: my_y_rank         = 0
  integer, public :: my_z_rank         = 0
  integer, public :: mpi_drain         = 50000
  integer, public :: grid_commx        = 0               
  integer, public :: grid_commy        = 0               
  integer, public :: grid_commz        = 0               
  integer, public :: my_xy_rank        = 0
  integer, public :: nprocxy           = 1
  integer, public :: grid_commxy       = 0               
  integer, public :: nthread           = 0
  integer, public :: ierror            = 0
  integer, public :: scalapack_context = 0
  integer, public :: my_prow           = 0
  integer, public :: my_pcol           = 0
  integer, public :: nprow             = 1
  integer, public :: npcol             = 1

  logical, public :: main_rank         = .false.
  logical, public :: real_calc         = .true.
  logical, public :: reciprocal_calc   = .true.

  ! variables for MPI in REMD
  integer, public :: mpi_comm_airplane = 0
  integer, public :: nproc_airplane    = 1
  integer, public :: my_airplane_rank  = 0
  integer, public :: my_replica_no     = 0
  integer, public :: my_country_no     = 0
  logical, public :: replica_main_rank = .false.

  ! variables for MPI in REPLICA
  integer, public :: nrep_per_proc     = 1

#ifdef USE_GPU
  integer, public :: my_node_local_rank = 0
#endif /* USE_GPU */

  ! subroutines
  public   :: get_loop_index
  public   :: get_para_range

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_loop_index
  !> @brief        divide loop indexes for MPI
  !! @authors      JJ, TM
  !! @param[in]    num_indexes : number of indexes
  !! @param[out]   istart      : start index
  !! @param[out]   iend        : end index
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine get_loop_index(num_indexes, istart, iend)

    ! formal arguments
    integer,                 intent(in)    :: num_indexes
    integer,                 intent(out)   :: istart
    integer,                 intent(out)   :: iend

    ! local variables
    integer                  :: quotient, remainder


    quotient  = num_indexes / nproc_city
    remainder = mod(num_indexes, nproc_city)

    if (my_city_rank <= (remainder - 1)) then
      quotient = quotient + 1
      istart   = quotient * my_city_rank + 1
      iend     = istart + quotient - 1
    else
      istart   = (quotient + 1)*remainder + &
                 quotient*(my_city_rank - remainder) + 1
      iend     = istart + quotient - 1
    end if

    return

  end subroutine get_loop_index

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_para_range
  !> @brief        Compute start and end numbers for parallel computation
  !! @authors      TA
  !! @param[in]    istart    : start index
  !! @param[in]    iend      : end index
  !! @param[in]    nthread   : number of threads
  !! @param[in]    id        : core ID
  !! @param[out]   id_istart : start index in ID core
  !! @param[out]   id_iend   : end index in ID core
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine get_para_range(istart, iend, nthread, id, id_istart, id_iend)

    ! formal arguments
    integer,                 intent(in)    :: istart
    integer,                 intent(in)    :: iend
    integer,                 intent(in)    :: nthread
    integer,                 intent(in)    :: id
    integer,                 intent(out)   :: id_istart
    integer,                 intent(out)   :: id_iend

    ! local vairbles
    integer                  :: iwork1, iwork2


    iwork1 =     (iend - istart + 1) / nthread
    iwork2 = mod((iend - istart + 1),  nthread)

    id_istart = id*iwork1 + istart + min(id, iwork2)
    id_iend = id_istart + iwork1 - 1

    if (iwork2 > id) &
      id_iend = id_iend + 1

    return

  end subroutine get_para_range

end module mpi_parallel_mod
