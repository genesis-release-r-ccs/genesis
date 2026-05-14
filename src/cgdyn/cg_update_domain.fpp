!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_update_domain_mod
!> @brief   Library subroutine used for integrators
!! @authors Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_update_domain_mod

  use cg_domain_mod
  use cg_pairlist_mod
  use cg_enefunc_mod
  use cg_communicate_str_mod
  use cg_communicate_mod
  use cg_migration_mod
  use cg_minimize_str_mod
  use cg_dynamics_str_mod
  use cg_boundary_str_mod
  use cg_pairlist_str_mod
  use cg_enefunc_str_mod
  use cg_domain_str_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use fileio_grotop_mod
  use molecules_str_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: domain_interaction_update_md
  public  :: domain_interaction_update_min
  public  :: do_migration

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    domain_interaction_update_md
  !> @brief        Updates interactions with MD
  !! @authors      JJ
  !! @param[in]    istep       : step number
  !! @param[in]    dynamics    : Dynamics information
  !! @param[inout] domain      : Domain information
  !! @param[inout] enefunc     : Potential energy functions information
  !! @param[inout] pairlist    : Pairlist information
  !! @param[inout] boundary    : Boundary condition information
  !! @param[inout] comm        : Communication information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine domain_interaction_update_md(istep, grotop, dynamics, molecule, &
                                          domain, enefunc, pairlist,         &
                                          boundary, comm)

    ! formal arguments
    integer,                 intent(in)    :: istep
    type(s_grotop),          intent(in)    :: grotop
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_molecule),        intent(inout) :: molecule
    type(s_domain),  target, intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_boundary),        intent(inout) :: boundary
    type(s_comm),            intent(inout) :: comm

    ! local variable
    integer                  :: ncell, nb, num_proc, i, j, k, l
    real(wp)                 :: dx, dy, dz, disp

    if (dynamics%lbupdate_period > 0) then
      if (mod(istep, dynamics%lbupdate_period) == 0) then
        call molecule_accum(boundary, domain, molecule)
        call setup_domain_lb(boundary, molecule, enefunc, domain)

        call define_enefunc_lb(grotop, molecule, domain, enefunc, &
                               comm)
        call setup_communicate(boundary, domain, comm)
        call update_communicate_size(domain, comm)
        call update_enefunc_contact(domain, enefunc)
        call communicate_contact(domain, comm, enefunc)
        call do_migration(domain, enefunc, boundary, comm, pairlist)
        call do_pairlist (domain, enefunc, boundary, pairlist)
        call copy_coord(domain)
      else if (mod(istep,dynamics%nbupdate_period) == 0) then
        call check_displacement(domain, disp)
        disp = sqrt(disp)
        if (disp > enefunc%buffer_min) then
          call do_migration(domain, enefunc, boundary, comm, pairlist)
          call do_pairlist (domain, enefunc, boundary, pairlist)
          call copy_coord(domain)
        end if
      end if
    else if (mod(istep,dynamics%nbupdate_period) == 0) then
      call check_displacement(domain, disp)
      disp = sqrt(disp)
      if (disp > enefunc%buffer_min) then
        call do_migration(domain, enefunc, boundary, comm, pairlist)
        call do_pairlist (domain, enefunc, boundary, pairlist)
        call copy_coord(domain)
      end if
    end if

    return

  end subroutine domain_interaction_update_md

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_displacement
  !> @brief        Checks coordinate displacement to decide migration
  !! @authors      JJ
  !! @param[in   ] domain      : Domain information
  !! @param[inout] disp        : maximum displacement length
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_displacement(domain, disp)

    ! formal arguments
    type(s_domain),  target, intent(in   ) :: domain
    real(wp),                intent(inout) :: disp

    ! local variable
    integer                  :: i
    real(wp)                 :: dx, dy, dz
    real(wip),       pointer :: coord(:,:), coord_prev(:,:)

    coord      => domain%coord
    coord_prev => domain%coord_prev

    disp = 0.0_wp
    do i = 1, domain%num_atom_domain
      dx   = coord_prev(i,1) - coord(i,1)
      dy   = coord_prev(i,2) - coord(i,2)
      dz   = coord_prev(i,3) - coord(i,3)
      disp = max(disp, dx*dx+dy*dy+dz*dz)
    end do
    call mpi_allreduce(mpi_in_place, disp, 1, mpi_wp_real, mpi_max, &
                       mpi_comm_country, ierror)
    disp = sqrt(disp)

    return

  end subroutine check_displacement

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    copy_coord
  !> @brief        Copies coordinates for next displacement check
  !! @authors      JJ
  !! @param[in   ] domain      : Domain information
  !! @param[inout] disp        : maximum displacement length
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine copy_coord(domain)

    ! formal arguments
    type(s_domain),  target, intent(in   ) :: domain

    ! local variable
    integer                  :: i
    real(wip),       pointer :: coord(:,:), coord_prev(:,:)

    coord      => domain%coord
    coord_prev => domain%coord_prev

    do i = 1, domain%num_atom_domain
      coord_prev(i,1) = coord(i,1)
      coord_prev(i,2) = coord(i,2)
      coord_prev(i,3) = coord(i,3)
    end do

    return

  end subroutine copy_coord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    domain_interaction_update_min
  !> @brief        Updates interactions with minimization
  !! @authors      JJ
  !! @param[in]    istep       : step number
  !! @param[in]    minimize    : Minimize information
  !! @param[inout] domain      : Domain information
  !! @param[inout] enefunc     : Potential energy functions information
  !! @param[inout] pairlist    : Pairlist information
  !! @param[inout] boundary    : Boundary condition information
  !! @param[inout] comm        : Communication information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine domain_interaction_update_min(istep, minimize, domain, enefunc, &
                                           pairlist, boundary, comm)

    ! formal arguments
    integer,                 intent(in)    :: istep
    type(s_minimize),        intent(in)    :: minimize
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_boundary),        intent(inout) :: boundary
    type(s_comm),            intent(inout) :: comm

    ! local variable
    integer                  :: ncell, nb, num_proc

    ncell = domain%num_cell_local
    nb    = domain%num_cell_boundary
    pairlist%univ_update = 0
    num_proc = domain%num_comm_proc

    if (mod(istep,minimize%nbupdate_period) == 0) then

      call do_migration(domain, enefunc, boundary, comm, pairlist)
      call do_pairlist (domain, enefunc, boundary, pairlist)

    end if

    return

  end subroutine domain_interaction_update_min

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    do_migration
  !> @brief        Migrates particles and energy functions
  !! @authors      JJ
  !! @param[inout] domain      : Domain information
  !! @param[inout] enefunc     : Potential energy functions information
  !! @param[inout] pairlist    : Pairlist information
  !! @param[inout] boundary    : Boundary condition information
  !! @param[inout] comm        : Communication information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine do_migration(domain, enefunc, boundary, comm, pairlist)

    ! formal arguments
    type(s_domain),   target, intent(inout) :: domain
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_boundary),         intent(inout) :: boundary
    type(s_comm),             intent(inout) :: comm
    type(s_pairlist),         intent(inout) :: pairlist

    call timer(TimerMigration, TimerOn)
    call timer(TimerMptl, TimerOn)
    if (enefunc%forcefield /= ForcefieldGroMartini) then
      call update_outgoing_charge(boundary, domain)
      call update_outgoing_nocharge(boundary, domain)
    else
      call update_outgoing_charge_martini(boundary, domain)
      call update_outgoing_nocharge_martini(boundary, domain)
    end if
    call timer(TimerMptl, TimerOff)
    call timer(TimerComm3, TimerOn)
    if (enefunc%forcefield /= ForcefieldGroMartini) then
      call communicate_ptl(domain, comm)
    else
      call communicate_ptl_martini(domain, comm)
    end if
    call timer(TimerComm3, TimerOff)
    call timer(TimerMptl, TimerOn)
    if (enefunc%forcefield /= ForcefieldGroMartini) then
      call update_incoming_ptl(domain)
    else
      call update_incoming_ptl_martini(domain)
    end if
    call timer(TimerMptl, TimerOff)

    call update_cell_size(domain, comm)
    call update_alloc_size_domain(domain, enefunc, pairlist)
    call update_cell_boundary(domain, enefunc, comm, boundary)

    call update_enefunc_aicg(domain, comm, enefunc)
    call update_alloc_size_enefunc(domain, enefunc, pairlist)
    call timer(TimerMigration, TimerOff)

    return

  end subroutine do_migration

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    do_pairlist
  !> @brief        Generates pairlist after migration
  !! @authors      JJ
  !! @param[inout] domain      : Domain information
  !! @param[inout] enefunc     : Potential energy functions information
  !! @param[inout] pairlist    : Pairlist information
  !! @param[inout] boundary    : Boundary condition information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine do_pairlist(domain, enefunc, boundary, pairlist)

    ! formal arguments
    type(s_domain),   target, intent(inout) :: domain
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_boundary),         intent(inout) :: boundary
    type(s_pairlist),         intent(inout) :: pairlist

    integer                   :: i
    real(wip),        pointer :: coord(:,:)
    real(wp),         pointer :: trans(:,:), coord_pbc(:,:)

    coord                => domain%coord
    trans                => domain%trans_vec
    coord_pbc            => domain%translated

    call timer(TimerPairList, TimerOn)

    if (enefunc%forcefield == ForcefieldRESIDCG) then

      if (boundary%type == BoundaryTypeNOBC) then
        !$omp parallel do private(i)
        do i = 1, domain%num_atom_domain + domain%num_atom_boundary
          coord_pbc(i,1) = real(coord(i,1),wp)
          coord_pbc(i,2) = real(coord(i,2),wp)
          coord_pbc(i,3) = real(coord(i,3),wp)
        end do
        !$omp end parallel do
      else if (boundary%type == BoundaryTypePBC) then
        !$omp parallel do private(i)
        do i = 1, domain%num_atom_domain + domain%num_atom_boundary
          coord_pbc(i,1) = real(coord(i,1),wp) + trans(i,1)
          coord_pbc(i,2) = real(coord(i,2),wp) + trans(i,2)
          coord_pbc(i,3) = real(coord(i,3),wp) + trans(i,3)
        end do
        !$omp end parallel do
      end if

      call update_alloc_size_pairlist(enefunc, pairlist)
      call update_pairlist_cg(coord_pbc, enefunc, domain, pairlist)
      if (enefunc%cg_pwmcos_calc) &
        call update_pairlist_cg_pwmcos(coord_pbc, enefunc, domain, pairlist)
      if (enefunc%cg_pwmcosns_calc) &
        call update_pairlist_cg_pwmcosns(coord_pbc, enefunc, domain, pairlist)

    else if (enefunc%forcefield == ForcefieldGroMartini) then

      !$omp parallel do private(i)
      do i = 1, domain%num_atom_domain + domain%num_atom_boundary
        coord_pbc(i,1) = real(coord(i,1),wp) + trans(i,1)
        coord_pbc(i,2) = real(coord(i,2),wp) + trans(i,2)
        coord_pbc(i,3) = real(coord(i,3),wp) + trans(i,3)
      end do
      !$omp end parallel do
      call update_alloc_size_pairlist(enefunc, pairlist)
      call update_pairlist_martini(coord_pbc, enefunc, domain, pairlist)
   
    end if

    call timer(TimerPairList, TimerOff)

    return

  end subroutine do_pairlist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_alloc_size_domain
  !> @brief        Checks the array size and reallocates if necessary
  !! @authors      JJ  
  !! @param[inout] domain      : Domain information
  !! @param[inout] enefunc     : Potential energy functions information
  !! @param[inout] pairlist    : Pairlist information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_alloc_size_domain(domain, enefunc, pairlist)

    ! formal arguments
    type(s_domain),           intent(inout) :: domain
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_pairlist),         intent(inout) :: pairlist
    integer                   :: num_atom

    num_atom = domain%num_atom_domain + domain%num_atom_boundary

    if (num_atom > MaxAtom_domain*3/4) then
      call copy_domain_information(1, domain)
      MaxAtom_domain = num_atom * 2
      num_atom = MaxAtom_domain
      call alloc_domain(domain, DomainDynvar_Atom, num_atom, 1, 1)
      call copy_domain_information(2, domain)
      call alloc_enefunc(enefunc, EneFuncCGElecInvList, num_atom)
      call alloc_enefunc(enefunc, EneFuncNonb, num_atom)
      call alloc_pairlist(pairlist, PairListCGExvNum, num_atom)
      call alloc_pairlist(pairlist, PairListCGExvList, &
                          Max_exv_nb15, num_atom)

      if (enefunc%cg_DNA_exv_calc) &
        call alloc_enefunc(enefunc, EneFuncCGDNAInvList, num_atom)
      if (enefunc%cg_DNA_base_pair_calc) &
        call alloc_enefunc(enefunc, EneFuncCGBaseInvList, num_atom)
      if (enefunc%cg_KH_calc) &
        call alloc_enefunc(enefunc, EneFuncCGKHInvList, num_atom)
      if (enefunc%cg_IDR_KH_calc) &
        call alloc_enefunc(enefunc, EneFuncCGIDRKHInvList, num_atom)
      if (enefunc%cg_IDR_HPS_calc) &
        call alloc_enefunc(enefunc, EneFuncCGIDRHPSInvList, num_atom)
      if (enefunc%cg_IDR_MPIPI_calc) &
        call alloc_enefunc(enefunc, EneFuncCGIDRMPIPIInvList, num_atom)
    end if

    if (enefunc%cg_ele_calc) then
      if (enefunc%num_cg_elec > Max_cg_elec*3/4) then
        Max_cg_elec = max(enefunc%num_cg_elec*2,min_atom)
        num_atom = Max_cg_elec
        call alloc_enefunc(enefunc, EneFuncCGElecList, num_atom)
        call alloc_pairlist(pairlist, PairListCGEleNum, num_atom)
        call alloc_pairlist(pairlist, PairListCGEleList, &
                            Max_elec_nb15, num_atom)
      end if
    end if
      
    if (enefunc%cg_DNA_exv_calc) then
      if (enefunc%num_cg_DNA > Max_cg_dna*3/4) then
        Max_cg_dna = enefunc%num_cg_DNA * 2
        num_atom = Max_cg_dna
        call alloc_enefunc(enefunc, EneFuncCGDNAList,  num_atom)
        call alloc_pairlist(pairlist, PairListCGDNAExvNum, num_atom)
        call alloc_pairlist(pairlist, PairListCGDNAExvList, &
                            Max_dna_exv_nb15, num_atom)
      end if
    end if

    if (enefunc%cg_DNA_base_pair_calc) then
      if (enefunc%num_cg_base > Max_cg_base*3/4) then
        Max_cg_base = enefunc%num_cg_base * 2
        num_atom = Max_CG_base
        call alloc_enefunc(enefunc, EneFuncCGBaseList, num_atom)
        call alloc_pairlist(pairlist, PairListCGDNABaseNum, num_atom)
        call alloc_pairlist(pairlist, PairListCGDNABaseList, &
                            Max_dna_base_nb15, num_atom)
      end if
    end if

    if (enefunc%cg_KH_calc) then
      if (enefunc%num_cg_KH > Max_cg_KH*3/4) then
        Max_cg_KH = enefunc%num_cg_KH * 2
        num_atom = Max_CG_KH
        call alloc_enefunc(enefunc, EneFuncCGKHList, num_atom)
        call alloc_pairlist(pairlist, PairListCGKHNum, num_atom)
        call alloc_pairlist(pairlist, PairListCGKHList, &
                            Max_KH_nb15, num_atom)
      end if
    end if

    if (enefunc%cg_IDR_KH_calc) then
      if (enefunc%num_cg_IDR_KH > Max_cg_IDR_KH*3/4) then
        Max_cg_IDR_KH = enefunc%num_cg_IDR_KH * 2
        num_atom = Max_CG_IDR_KH
        call alloc_enefunc(enefunc, EneFuncCGIDRKHList, num_atom)
        call alloc_pairlist(pairlist, PairListCGIDRKHNum, num_atom)
        call alloc_pairlist(pairlist, PairListCGIDRKHList, &
                            Max_IDR_KH_nb15, num_atom)
      end if
    end if

    if (enefunc%cg_IDR_HPS_calc) then
      if (enefunc%num_cg_IDR_HPS > Max_cg_IDR_HPS*3/4) then
        Max_cg_IDR_HPS = enefunc%num_cg_IDR_HPS * 2
        num_atom = Max_CG_IDR_HPS
        call alloc_enefunc(enefunc, EneFuncCGIDRHPSList, num_atom)
        call alloc_pairlist(pairlist, PairListCGIDRHPSNum, num_atom)
        call alloc_pairlist(pairlist, PairListCGIDRHPSList, &
                            Max_IDR_HPS_nb15, num_atom)
      end if
    end if

    if (enefunc%cg_IDR_MPIPI_calc) then
      if (enefunc%num_cg_IDR_MPIPI > Max_cg_IDR_MPIPI*3/4) then
        Max_cg_IDR_MPIPI = enefunc%num_cg_IDR_MPIPI * 2
        num_atom = Max_CG_IDR_MPIPI
        call alloc_enefunc(enefunc, EneFuncCGIDRMPIPIList, num_atom)
        call alloc_pairlist(pairlist, PairListCGIDRMPIPINum, num_atom)
        call alloc_pairlist(pairlist, PairListCGIDRMPIPIList, &
            Max_IDR_MPIPI_nb15, num_atom)
      end if
    end if

    return

  end subroutine update_alloc_size_domain

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_alloc_size_enefunc
  !> @brief        Checks the array size and reallocates if necessary
  !! @authors      JJ  
  !! @param[inout] domain      : Domain information
  !! @param[inout] enefunc     : Potential energy functions information
  !! @param[inout] pairlist    : Pairlist information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_alloc_size_enefunc(domain, enefunc, pairlist)

    ! formal arguments
    type(s_domain),           intent(inout) :: domain
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_pairlist),         intent(inout) :: pairlist
    integer                   :: num, num1

    num = enefunc%num_bondsq_domain + enefunc%num_bond_domain
    if (num > MaxBond*3/4) then
      call copy_bond_information(1, domain, enefunc)
      MaxBond = num * 2
      num = MaxBond
      call alloc_enefunc(enefunc, EneFuncBond, num, 1)
      call copy_bond_information(2, domain, enefunc)
    end if

    num = enefunc%num_angle_flexible_domain + enefunc%num_angle_local_domain &
        + enefunc%num_angle_domain
    if (num > MaxAngl*3/4) then
      call copy_angl_information(1, domain, enefunc)
      MaxAngl = num * 2
      num = MaxAngl
      call alloc_enefunc(enefunc, EneFuncAngl, num, 1)
      call copy_angl_information(2, domain, enefunc)
    end if

    num = enefunc%num_dihe_flexible_domain + enefunc%num_dihe_local_domain &
        + enefunc%num_dihe_domain
    if (num > MaxDihe*3/4) then
      call copy_dihe_information(1, domain, enefunc)
      MaxDihe = num * 2
      num = MaxDihe
      call alloc_enefunc(enefunc, EneFuncDihe, num, 1)
      call copy_dihe_information(2, domain, enefunc)
    end if

    if (enefunc%forcefield == ForcefieldGroMartini) then
      if (MaxExcl < MaxBond+MaxAngl+MaxDihe) then
        MaxExcl = MaxBond + MaxAngl + MaxDihe
        call alloc_enefunc(enefunc, EneFuncExcl, MaxExcl, 1)
      end if
    end if

    num = enefunc%num_stack_domain
    if (num > MaxStack*3/4) then
      call copy_stack_information(1, domain, enefunc)
      MaxStack = num * 2
      num = MaxStack
      call alloc_enefunc(enefunc, EneFuncBaseStack, num, 1)
      call copy_stack_information(2, domain, enefunc)
    end if

    num = enefunc%num_contact_domain + enefunc%num_contact_boundary
    if (num > MaxContact*3/4) then
      call copy_contact_information(1, domain, enefunc)
      MaxContact = num * 2
      num = MaxContact
      call alloc_enefunc(enefunc, EneFuncContact, num, 1)
      call copy_contact_information(2, domain, enefunc)
    end if

    num = enefunc%num_pwmcos_domain
    if (num > MaxPwmCos*3/4) then
      call copy_pwmcos_information(1, domain, enefunc)
      MaxPwmCos = num * 2
      num = MaxPwmCos
      call alloc_enefunc(enefunc, EneFuncPWMcos, num)
      call copy_pwmcos_information(2, domain, enefunc)
      call alloc_pairlist(pairlist, PairListPWMCOSNum, num)
      call alloc_pairlist(pairlist, PairListPwmCosList, &
                          Max_pwmcos_nb15, MaxPwmCos)
    end if

    num = enefunc%num_pwmcosns_domain
    if (num > MaxPwmCosns*3/4) then
      call copy_pwmcosns_information(1, domain, enefunc)
      MaxPwmCosns = num * 2
      num = MaxPwmCosns * 2
      call alloc_enefunc(enefunc, EneFuncPWMcosns, num)
      call copy_pwmcosns_information(2, domain, enefunc)
      call alloc_pairlist(pairlist, PairListPWMCOSnsNum, num)
      call alloc_pairlist(pairlist, PairListPwmCosnsList, &
                          Max_pwmcosns_nb15, MaxPwmCosns)
    end if

    num = enefunc%num_rest_domain
    if (num > MaxRest*3/4) then
      call copy_rest_information(1, domain, enefunc)
      MaxRest = num * 2
      num = MaxRest
      call alloc_enefunc(enefunc, EneFuncRestDomain, num)
      call copy_rest_information(2, domain, enefunc)
    end if

    num  = Max_alloc_size1
    num1 = Max_alloc_size2
    Max_alloc_size1 = max(8*MaxAtom_Domain, 2*MaxBond, 4*MaxAngl, &
                          3*MaxDihe, 3*MaxStack, 60*MaxPwmCos,    &
                          24*MaxPwmCosns, 3*MaxContact)
    Max_alloc_size2 = max(6*MaxAtom_Domain, 3*MaxBond, 4*MaxAngl, &
                          6*MaxDihe, 3*MaxStack, 5*MaxPwmCos,     &
                          10*MaxPwmCosns, 4*MaxContact)
    if (num /= Max_alloc_size1 .or. num1 /= Max_alloc_size2)   &
    call alloc_domain(domain, DomainPtlArray, Max_alloc_size1, &
                      Max_alloc_size2, 1)

    return

  end subroutine update_alloc_size_enefunc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_alloc_size_pairlist
  !> @brief        Reallocates pairlist if necessary
  !! @authors      JJ
  !! @param[inout] enefunc     : Potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_alloc_size_pairlist(enefunc, pairlist)

    ! formal arguments
    type(s_enefunc),          intent(in   ) :: enefunc 
    type(s_pairlist),         intent(inout) :: pairlist
    integer                   :: num
  
    num = pairlist%realloc_exv
    call mpi_allreduce(mpi_in_place, num, 1, mpi_integer, mpi_max, &
                       mpi_comm_country, ierror) 
    if (num > 0) then
      Max_exv_nb15 = Max_exv_nb15 * 3 / 2
      call alloc_pairlist(pairlist, PairListCGExvList, &
                          Max_exv_nb15, MaxAtom_domain)
    end if

    if (enefunc%cg_DNA_exv_calc) then
      num = pairlist%realloc_DNA_exv
      call mpi_allreduce(mpi_in_place, num, 1, mpi_integer, mpi_max, &
                         mpi_comm_country, ierror) 
      if (num > 0) then
        Max_dna_exv_nb15 = Max_dna_exv_nb15 * 3 / 2
        call alloc_pairlist(pairlist, PairListCGDNAExvList, &
                            Max_dna_exv_nb15, Max_cg_DNA)
      end if
    end if

    if (enefunc%cg_DNA_base_pair_calc) then
      num = pairlist%realloc_DNA_base
      call mpi_allreduce(mpi_in_place, num, 1, mpi_integer, mpi_max, &
                         mpi_comm_country, ierror) 
      if (num > 0) then
        Max_dna_base_nb15 = Max_dna_base_nb15 * 3 / 2
        call alloc_pairlist(pairlist, PairListCGDNABaseList, &
                            Max_dna_base_nb15, Max_cg_Base)
      end if
    end if

    if (enefunc%cg_ele_calc) then
      num = pairlist%realloc_elec
      call mpi_allreduce(mpi_in_place, num, 1, mpi_integer, mpi_max, &
                         mpi_comm_country, ierror) 
      if (num > 0) then
        Max_elec_nb15 = Max_elec_nb15 * 3 / 2
        call alloc_pairlist(pairlist, PairListCGEleList, &
                            Max_elec_nb15, Max_cg_Elec)
      end if
    end if

    if (enefunc%cg_IDR_KH_calc) then
      num = pairlist%realloc_IDR_KH
      call mpi_allreduce(mpi_in_place, num, 1, mpi_integer, mpi_max, &
                         mpi_comm_country, ierror) 
      if (num > 0) then
        Max_IDR_KH_nb15 = Max_IDR_KH_nb15 * 3 / 2
        call alloc_pairlist(pairlist, PairListCGIDRKHList, &
                            Max_IDR_KH_nb15, Max_cg_IDR_KH)
      end if
    end if

    if (enefunc%cg_KH_calc) then
      num = pairlist%realloc_KH
      call mpi_allreduce(mpi_in_place, num, 1, mpi_integer, mpi_max, &
                         mpi_comm_country, ierror) 
      if (num > 0) then
        Max_KH_nb15 = Max_KH_nb15 * 3 / 2
        call alloc_pairlist(pairlist, PairListCGKHList, &
                          Max_KH_nb15, Max_cg_KH)
      end if
    end if

    if (enefunc%cg_IDR_HPS_calc) then
      num = pairlist%realloc_IDR_HPS
      call mpi_allreduce(mpi_in_place, num, 1, mpi_integer, mpi_max, &
                         mpi_comm_country, ierror) 
      if (num > 0) then
        Max_IDR_HPS_nb15 = Max_IDR_HPS_nb15 * 3 / 2
        call alloc_pairlist(pairlist, PairListCGIDRHPSList, &
                            Max_IDR_HPS_nb15, Max_cg_IDR_HPS)
      end if
    end if

    if (enefunc%cg_IDR_MPIPI_calc) then
      num = pairlist%realloc_IDR_MPIPI
      call mpi_allreduce(mpi_in_place, num, 1, mpi_integer, mpi_max, &
          mpi_comm_country, ierror) 
      if (num > 0) then
        Max_IDR_MPIPI_nb15 = Max_IDR_MPIPI_nb15 * 3 / 2
        call alloc_pairlist(pairlist, PairListCGIDRMPIPIList, &
            Max_IDR_MPIPI_nb15, Max_cg_IDR_MPIPI)
      end if
    end if

    if (enefunc%cg_pwmcos_calc) then
      num = pairlist%realloc_PWMcos
      call mpi_allreduce(mpi_in_place, num, 1, mpi_integer, mpi_max, &
                         mpi_comm_country, ierror) 
      if (num > 0) then
        Max_pwmcos_nb15 = Max_pwmcos_nb15 * 3 / 2
        call alloc_pairlist(pairlist, PairListPwmCosList, &
                            Max_pwmcos_nb15, MaxPwmCos)
      end if
    end if

    if (enefunc%cg_pwmcosns_calc) then
      num = pairlist%realloc_PWMcosns
      call mpi_allreduce(mpi_in_place, num, 1, mpi_integer, mpi_max, &
                         mpi_comm_country, ierror) 
      if (num > 0) then
        Max_pwmcosns_nb15 = Max_pwmcosns_nb15 * 3 / 2
        call alloc_pairlist(pairlist, PairListPwmCosnsList, &
                            Max_pwmcosns_nb15, MaxPwmCosns)
      end if
    end if

    return

  end subroutine update_alloc_size_pairlist

end module cg_update_domain_mod
