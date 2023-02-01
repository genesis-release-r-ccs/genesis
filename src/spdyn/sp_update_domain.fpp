!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_update_domain_mod
!> @brief   library subroutine used for integrator
!! @authors Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_update_domain_mod

  use sp_domain_mod
  use sp_pairlist_mod
  use sp_enefunc_mod
  use sp_communicate_mod
  use sp_migration_mod
  use sp_minimize_str_mod
  use sp_dynamics_str_mod
  use sp_constraints_str_mod
  use sp_boundary_str_mod
  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public  :: domain_interaction_update
  public  :: domain_update_pio
  ! FEP
  public  :: domain_interaction_update_fep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    domain_interaction_update
  !> @brief        update interactions with MD
  !! @authors      JJ
  !! @param[in]    istep         : step number
  !! @param[in]    update_period : update period for pairlist
  !! @param[inout] domain        : domain information
  !! @param[inout] enefunc       : potential energy functions information
  !! @param[inout] pairlist      : pairlist information
  !! @param[inout] boundary      : boundary condition information
  !! @param[inout] constraints   : constraints information
  !! @param[inout] comm          : communication information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine domain_interaction_update(istep, update_period, domain, enefunc, &
                                       pairlist, boundary, constraints, comm)

    ! formal arguments
    integer,                 intent(in)    :: istep
    integer,                 intent(in)    :: update_period
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_boundary),        intent(inout) :: boundary
    type(s_constraints),     intent(inout) :: constraints
    type(s_comm),            intent(inout) :: comm

    ! local variable
    integer                  :: ncell, nb, water_atom


    ncell = domain%num_cell_local
    nb    = domain%num_cell_boundary
    pairlist%univ_update = 0

    if (constraints%water_type == TIP4) then
      water_atom = 4
    else if (constraints%water_type == TIP3) then
      water_atom = 3
    else if (constraints%water_type == TIP2) then
      water_atom = 2
    else if (constraints%water_type == TIP1) then
      water_atom = 1
    end if

    if (mod(istep,update_period) == 0) then

#ifdef DEBUG
      if (main_rank) then
        write(MsgOut,'(a,3f15.8)') &
             "Debuging > bsize_[x,y,z]", &
             boundary%box_size_x,boundary%box_size_y,boundary%box_size_z
        write(MsgOut,'(a,2i8)')    &
             "Debuging > ncell, nb",ncell,nb
      end if
#endif

      pairlist%univ_update = 1

      if (istep == update_period) then

        call alloc_domain(domain, DomainPtlMove, ncell+nb, ncell, 1)
        call alloc_domain(domain, DomainWaterMove, ncell+nb, ncell, &
                          water_atom)

        call alloc_constraints(constraints, ConstraintsHGroupMove, &
                               constraints%connect, ncell+nb)

        call alloc_enefunc(enefunc, EneFuncBondCell, ncell, ncell+nb)

      end if

      call timer(TimerUpdate, TimerOn)
      call update_outgoing_atom (boundary, constraints, domain)
      call update_outgoing_HGr  (boundary, constraints, domain)
      call update_outgoing_water(water_atom, boundary, domain)

      call timer(TimerComm3, TimerOn)
      call communicate_constraints(domain, comm, constraints)
      call timer(TimerComm3, TimerOff)

      call update_incoming_atom (constraints, domain)
      call update_incoming_HGr  (constraints, domain)
      call update_incoming_water(water_atom, domain)

      call update_cell_size_constraints    (domain, comm, constraints)
      call update_cell_boundary(domain, comm, boundary, constraints)

      call update_enefunc(domain, comm, enefunc, constraints)

      call timer(TimerUpdate, TimerOff)

      call timer(TimerPairList, TimerOn)

      if (enefunc%pairlist_check) then
        call update_pairlist_pbc_check(enefunc, domain, pairlist)
      else
        call update_pairlist_pbc(enefunc, domain, pairlist)
      end if

      call timer(TimerPairList, TimerOff)

    end if

    return

  end subroutine domain_interaction_update
 
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    domain_update_pio
  !> @brief        update domain after reading parallel I/O
  !! @authors      JJ
  !! @param[inout] domain        : domain information
  !! @param[inout] enefunc       : potential energy functions information
  !! @param[inout] boundary      : boundary condition information
  !! @param[inout] constraints   : constraints information
  !! @param[inout] comm          : communication information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine domain_update_pio(domain, enefunc, boundary, constraints, comm)

    ! formal arguments
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_boundary),        intent(inout) :: boundary
    type(s_constraints),     intent(inout) :: constraints
    type(s_comm),            intent(inout) :: comm

    ! local variable
    integer                  :: ncell, nb, water_atom


    ncell = domain%num_cell_local
    nb    = domain%num_cell_boundary

    if (constraints%water_type == TIP4) then
      water_atom = 4
    else if (constraints%water_type == TIP3) then
      water_atom = 3
    else if (constraints%water_type == TIP2) then
      water_atom = 2
    else if (constraints%water_type == TIP1) then
      water_atom = 1
    end if

    call alloc_domain(domain, DomainPtlMove, ncell+nb, ncell, 1)
    call alloc_domain(domain, DomainWaterMove, ncell+nb, ncell, &
                      water_atom)

    call alloc_constraints(constraints, ConstraintsHGroupMove, &
                           constraints%connect, ncell+nb)

    call update_outgoing_atom (boundary, constraints, domain)
    call update_outgoing_HGr  (boundary, constraints, domain)
    call update_outgoing_water(water_atom, boundary, domain)

    call communicate_constraints(domain, comm, constraints)

    call update_incoming_atom (constraints, domain)
    call update_incoming_HGr  (constraints, domain)
    call update_incoming_water(water_atom, domain)

    call update_cell_size_constraints    (domain, comm, constraints)
    call update_cell_boundary(domain, comm, boundary, constraints)

    call dealloc_domain(domain, DomainPtlMove)
    call dealloc_domain(domain, DomainWaterMove)
    call dealloc_constraints(constraints, ConstraintsHGroupMove)

    return

  end subroutine domain_update_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    domain_interaction_update_fep
  !> @brief        update interactions with MD for FEP
  !! @authors      HO
  !! @param[in]    istep         : step number
  !! @param[in]    update_period : update period for pairlist
  !! @param[inout] domain        : domain information
  !! @param[inout] enefunc       : potential energy functions information
  !! @param[inout] pairlist      : pairlist information
  !! @param[inout] boundary      : boundary condition information
  !! @param[inout] constraints   : constraints information
  !! @param[inout] comm          : communication information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine domain_interaction_update_fep(istep, update_period, domain, &
                                           enefunc, pairlist, boundary,  &
                                           constraints, comm)

    ! formal arguments
    integer,                 intent(in)    :: istep
    integer,                 intent(in)    :: update_period
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_boundary),        intent(inout) :: boundary
    type(s_constraints),     intent(inout) :: constraints
    type(s_comm),            intent(inout) :: comm

    ! local variable
    integer                  :: ncell, nb, water_atom


    ncell = domain%num_cell_local
    nb    = domain%num_cell_boundary
    pairlist%univ_update = 0

    if (constraints%water_type == TIP4) then
      water_atom = 4
    else if (constraints%water_type == TIP3) then
      water_atom = 3
    else if (constraints%water_type == TIP2) then
      water_atom = 2
    else if (constraints%water_type == TIP1) then
      water_atom = 1
    end if

    if (mod(istep,update_period) == 0) then

#ifdef DEBUG
      if (main_rank) then
        write(MsgOut,'(a,3f15.8)') &
             "Debuging > bsize_[x,y,z]", &
             boundary%box_size_x,boundary%box_size_y,boundary%box_size_z
        write(MsgOut,'(a,2i8)')    &
             "Debuging > ncell, nb",ncell,nb
      end if
#endif

      pairlist%univ_update = 1

      if (istep == update_period) then

        ! FEP
        call alloc_domain(domain, DomainPtlMove_FEP, ncell+nb, ncell, 1)
        call alloc_domain(domain, DomainWaterMove_FEP, ncell+nb, ncell, &
                          water_atom)

        call alloc_constraints(constraints, ConstraintsHGroupMove_FEP, &
                               constraints%connect, ncell+nb)

        call alloc_enefunc(enefunc, EneFuncBondCell_FEP, ncell, ncell+nb)

      end if

      call timer(TimerUpdate, TimerOn)

      ! FEP
      call update_outgoing_atom_fep (boundary, constraints, domain)
      call update_outgoing_HGr_fep  (boundary, constraints, domain)
      call update_outgoing_water_fep(water_atom, boundary, domain)

      call timer(TimerComm3, TimerOn)

      ! FEP
      call communicate_constraints_fep(domain, comm, constraints)

      call timer(TimerComm3, TimerOff)

      ! FEP
      call update_incoming_atom_fep (constraints, domain)
      call update_incoming_HGr_fep  (constraints, domain)
      call update_incoming_water_fep(water_atom, domain)

      call update_cell_size_constraints    (domain, comm, constraints)

      ! FEP
      call update_cell_boundary_fep(domain, comm, boundary, constraints)

      call update_enefunc_fep(domain, comm, enefunc, constraints)

      call timer(TimerUpdate, TimerOff)

      call timer(TimerPairList, TimerOn)

      if (enefunc%pairlist_check) then
        call update_pairlist_pbc_check(enefunc, domain, pairlist)
      else
        call update_pairlist_pbc(enefunc, domain, pairlist)
      endif

      call timer(TimerPairList, TimerOff)

    end if

    return

  end subroutine domain_interaction_update_fep 

end module sp_update_domain_mod
