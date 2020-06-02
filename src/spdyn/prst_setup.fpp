!--------1---------2---------3---------4---------5---------6---------7---------8
! 
!> Program  prst_setup
!! @brief   setup restart file for domain decomposition MD
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

program prst_setup

  use pr_run_prst_setup_mod
  use pr_control_mod
  use pr_huge_molecule_mod
  use messages_mod
  use mpi_parallel_mod
  use timers_mod
  use mpi

  implicit none

  ! local variables
  character(100)    :: ctrl_filename
  type(s_ctrl_data) :: ctrl_data
  integer           :: ierr


  ! mpi settings
  ! 
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_prst_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc_prst,   ierr)

  main_rank = (my_prst_rank == 0)


  ! show usage
  !
  call usage(ctrl_filename)

  call timer(TimerDynamics, TimerOn)

  ! [Step1] Read control file
  !
  if (main_rank) then
    write(MsgOut,'(A)') '[STEP1] Read Control Parameters for output restart'
    write(MsgOut,'(A)') ' '
  end if

  call control(ctrl_filename, ctrl_data)


  ! [Step2] Output restart file
  !
  if (main_rank) then
    write(MsgOut,'(A)') '[STEP2] Output parallel I/O restart file'
    write(MsgOut,'(A)') ' '
  end if

  call run_prst_setup(ctrl_data)

  call timer(TimerDynamics, TimerOff)
  call output_time_prst(MPI_COMM_WORLD, nproc_prst)

  ! mpi settings
  !
  call MPI_Finalize(ierr)

  stop

end program prst_setup
