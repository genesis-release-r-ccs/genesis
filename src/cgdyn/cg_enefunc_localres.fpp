!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_enefunc_localres_mod
!> @brief   Provides restraint energy functions.
!! @authors Chigusa Kobayashi (CK), Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8
  
#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_enefunc_localres_mod

  use cg_enefunc_str_mod
  use cg_domain_str_mod
  use fileio_localres_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: setup_enefunc_localres
! private :: setup_enefunc_localres_bond
! private :: setup_enefunc_localres_angle
! private :: setup_enefunc_localres_dihed

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_localres
  !> @brief        Defines local restraints for each cell in potential energy 
  !                function.  
  !! @authors      CK
  !! @param[in]    localres : local restraint information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_localres(localres, domain, enefunc)

    ! formal arguments
    type(s_localres), target, intent(in)    :: localres
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(inout) :: enefunc


    ! bond
    !
!   call setup_enefunc_localres_bond(localres, domain, enefunc)

    ! angle
    !
!   call setup_enefunc_localres_angle(localres, domain, enefunc)

    ! dihedral
    !
!   call setup_enefunc_localres_dihed(localres, domain, enefunc)


    return

  end subroutine setup_enefunc_localres

end module cg_enefunc_localres_mod
