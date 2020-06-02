!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_restart_mod
!> @brief   setup variables and structures in MD simulaton
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_restart_mod

  use sp_dynvars_str_mod
  use sp_dynamics_str_mod
  use molecules_str_mod
  use fileio_rst_mod
  use random_mod
  use messages_mod
  use mpi_parallel_mod
 
  implicit none
  private

  ! subroutines
  public  :: setup_restart_pre
  public  :: setup_restart_post

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_restart_pre
  !> @brief        setup restart for molecule information
  !! @authors      NT
  !! @param[in]    rst        : restart file information
  !! @param[inout] molecule   : molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_restart_pre(rst, molecule)

    ! formal arguments
    type(s_rst),             intent(in)    :: rst
    type(s_molecule),        intent(inout) :: molecule

    ! local variables
    integer                  :: natom


    if (rst%rstfile_type == RstfileTypeMin) then

      natom = molecule%num_atoms

      molecule%atom_coord(1:3,1:natom) = rst%coord(1:3,1:natom)

      if (main_rank) then
        write(MsgOut,'(A)') 'Setup_Restart_Pre> Coordinates were replaced'
        write(MsgOut,'(A)') ''
      end if

    else if (rst%rstfile_type == RstfileTypeMd .or. &
             rst%rstfile_type == RstfileTypeRemd .or. &
             rst%rstfile_type == RstfileTypeRpath) then

      natom = molecule%num_atoms
      molecule%atom_coord(1:3,1:natom)    = rst%coord(1:3,1:natom)
      molecule%atom_velocity(1:3,1:natom) = rst%velocity(1:3,1:natom)

      if (main_rank) then
        write(MsgOut,'(A)') 'Setup_Restart_Pre> Coordinates and velocities were replaced'
        write(MsgOut,'(A)') ''
      end if

    end if

    return

  end subroutine setup_restart_pre

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_restart_post
  !> @brief        setup restart for other variables
  !! @authors      NT
  !! @param[in]    rst      : restart file information
  !! @param[inout] dynamics : dynamics information
  !! @param[inout] dynvars  : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_restart_post(rst, dynamics, dynvars)

    ! formal arguments
    type(s_rst),             intent(in)    :: rst
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_dynvars),         intent(inout) :: dynvars

    if (rst%rstfile_type == RstfileTypeMin) then

      dynamics%restart = .false.

    else if (rst%rstfile_type == RstfileTypeMd .or. &
             rst%rstfile_type == RstfileTypeRemd .or.  &
             rst%rstfile_type == RstfileTypeRpath) then

      dynamics%restart = .true.
      dynvars%thermostat_momentum    = rst%thermostat_momentum
      dynvars%barostat_momentum(1:3) = rst%barostat_momentum(1:3)

      ! overwrite dynamics%iseed by rst%iseed when dyn_info%iseed is not set
      !
      if (dynamics%iseed_read) then
        dynamics%iseed_init_velocity = rst%iseed
      endif

      if (allocated(rst%random) .and. dynamics%iseed_read)  then
        if (main_rank) then
          write(MsgOut,'(A)') 'Setup_Restart_Post> Read random seed '//&
                              ' from RST file and Overwrite it'
          write(MsgOut,'(A)') ''
        endif
        call random_stock_mpi_frombyte(mpi_comm_country, rst%random)
        call random_pull_stock
      else
        call random_term
        call random_init(dynamics%iseed)
      endif

      if (main_rank) then
        write(MsgOut,'(A)') 'Setup_Restart_Post> Parameters were replaced'
        write(MsgOut,'(A)') ''
      end if

    end if

    return

  end subroutine setup_restart_post

end module sp_restart_mod
