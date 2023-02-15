!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_fep_energy_mod
!> @brief   calculate energy differnces between adjacent states in FEP 
!! @authors Hiraku Oshima (HO)
!
!  (c) Copyright 2019 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_fep_energy_mod

  use sp_energy_pme_mod
  use sp_energy_mod          
  use sp_dynvars_str_mod
  use sp_boundary_str_mod
  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_energy_str_mod
  use sp_ensemble_str_mod
  use sp_remd_str_mod
  use sp_alchemy_str_mod
  use sp_ensemble_str_mod
  use sp_domain_str_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use math_libs_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: compute_fep_energy
  public  :: assign_condition_feprest
  public  :: assign_lambda

  private :: compute_energy_ref
  private :: compute_energy_tgt
  private :: assign_condition_fep
  private :: assign_condition_feprest_angle_ub
  private :: assign_condition_feprest_lj
  private :: assign_condition_feprest_internal
  private :: set_lambda_table_fep
  private :: set_fepgrp_nonb

contains

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_fep_energy(domain, enefunc, dynvars, &
                               pairlist, ensemble, boundary,      &
                               alchemy, remd)

    ! formal arguments
    type(s_domain),    target, intent(inout) :: domain
    type(s_enefunc),           intent(inout) :: enefunc
    type(s_dynvars),   target, intent(inout) :: dynvars
    type(s_pairlist),          intent(in)    :: pairlist
    type(s_ensemble),          intent(in)    :: ensemble
    type(s_boundary),          intent(inout) :: boundary
    type(s_alchemy),   target, intent(inout) :: alchemy
    type(s_remd),              intent(inout) :: remd

    ! local variables
    integer,         pointer :: lambid
    integer                  :: lambid_forward, lambid_backward
    integer                  :: replica_id
    integer                  :: param_id, param_id_forward, param_id_backward
    integer                  :: i, ix, j, k
    type(s_energy)           :: energy_ref
    type(s_energy)           :: energy_tgt

    lambid    => alchemy%lambid

    if (alchemy%fep_md_type == FEP_Parallel) then
      ! FEP-REMD case

      ! energy calculation by changing lambdas
      replica_id = my_country_no + 1
      param_id = remd%repid2parmsetid(replica_id)

      call compute_energy_ref(domain, enefunc, boundary, dynvars%energy, energy_ref)
      dynvars%energy%deltU_fep(1) = energy_ref%total

      ! evalutate total energy for target state
      ! and calculate the difference from reference state

      ! shift to the backward state
      param_id_backward = param_id - 1
      if (param_id_backward > 0) then
        call assign_condition_fep(param_id_backward, alchemy, domain,&
                                  ensemble, enefunc, remd)
        call assign_lambda(alchemy, domain, enefunc)
        call compute_energy_tgt(domain, enefunc, pairlist, boundary, energy_tgt)
        dynvars%energy%deltU_fep(2)=energy_tgt%total-energy_ref%total
      end if

      ! shift to the forward state
      param_id_forward  = param_id + 1
      if (param_id_forward <= remd%nreplicas(1)) then
        call assign_condition_fep(param_id_forward, alchemy, domain, &
                                  ensemble, enefunc, remd)
        call assign_lambda(alchemy, domain, enefunc)
        call compute_energy_tgt(domain, enefunc, pairlist, boundary, energy_tgt)
        dynvars%energy%deltU_fep(3)=energy_tgt%total-energy_ref%total
      end if

      ! return to the original state
      call assign_condition_fep(param_id, alchemy, domain, ensemble, &
                                enefunc, remd)
      call assign_lambda(alchemy, domain, enefunc)
    else
      ! FEP-MD case

      ! evaluate total energy for reference state
      call compute_energy_ref(domain, enefunc, boundary, dynvars%energy, energy_ref)
      dynvars%energy%deltU_fep(1) = energy_ref%total

      ! evalutate total energy for target state
      ! and calculate the difference from reference state
      if (alchemy%fep_direction == FEP_Bothsides) then

        if (lambid == 1) then
          lambid_backward = lambid
          lambid_forward  = lambid + 1
        else if (lambid == alchemy%num_fep_windows) then
          lambid_backward = lambid - 1
          lambid_forward  = lambid
        else
          lambid_backward = lambid - 1
          lambid_forward  = lambid + 1
        end if

        enefunc%lambljA      = alchemy%lambljA(lambid_backward)
        enefunc%lambljB      = alchemy%lambljB(lambid_backward)
        enefunc%lambelA      = alchemy%lambelA(lambid_backward)
        enefunc%lambelB      = alchemy%lambelB(lambid_backward)
        enefunc%lambbondA    = alchemy%lambbondA(lambid_backward)
        enefunc%lambbondB    = alchemy%lambbondB(lambid_backward)
        enefunc%lambrest     = alchemy%lambrest(lambid_backward)
        call assign_lambda(alchemy, domain, enefunc)
        call compute_energy_tgt(domain, enefunc, pairlist, boundary, energy_tgt)
        dynvars%energy%deltU_fep(2)=energy_tgt%total-energy_ref%total

        enefunc%lambljA      = alchemy%lambljA(lambid_forward)
        enefunc%lambljB      = alchemy%lambljB(lambid_forward)
        enefunc%lambelA      = alchemy%lambelA(lambid_forward)
        enefunc%lambelB      = alchemy%lambelB(lambid_forward)
        enefunc%lambbondA    = alchemy%lambbondA(lambid_forward)
        enefunc%lambbondB    = alchemy%lambbondB(lambid_forward)
        enefunc%lambrest     = alchemy%lambrest(lambid_forward)
        call assign_lambda(alchemy, domain, enefunc)
        call compute_energy_tgt(domain, enefunc, pairlist, boundary, energy_tgt)
        dynvars%energy%deltU_fep(3)=energy_tgt%total-energy_ref%total

      else if (alchemy%fep_direction == FEP_Forward) then

        lambid_forward       = lambid + 1
        enefunc%lambljA      = alchemy%lambljA(lambid_forward)
        enefunc%lambljB      = alchemy%lambljB(lambid_forward)
        enefunc%lambelA      = alchemy%lambelA(lambid_forward)
        enefunc%lambelB      = alchemy%lambelB(lambid_forward)
        enefunc%lambbondA    = alchemy%lambbondA(lambid_forward)
        enefunc%lambbondB    = alchemy%lambbondB(lambid_forward)
        enefunc%lambrest     = alchemy%lambrest(lambid_forward)
        call assign_lambda(alchemy, domain, enefunc)
        call compute_energy_tgt(domain, enefunc, pairlist, boundary, energy_tgt)
        dynvars%energy%deltU_fep(2)=energy_tgt%total-energy_ref%total

      else if (alchemy%fep_direction == FEP_Reverse) then

        lambid_backward      = lambid - 1
        enefunc%lambljA      = alchemy%lambljA(lambid_backward)
        enefunc%lambljB      = alchemy%lambljB(lambid_backward)
        enefunc%lambelA      = alchemy%lambelA(lambid_backward)
        enefunc%lambelB      = alchemy%lambelB(lambid_backward)
        enefunc%lambbondA    = alchemy%lambbondA(lambid_backward)
        enefunc%lambbondB    = alchemy%lambbondB(lambid_backward)
        enefunc%lambrest     = alchemy%lambrest(lambid_backward)
        call assign_lambda(alchemy, domain, enefunc)
        call compute_energy_tgt(domain, enefunc, pairlist, boundary, energy_tgt)
        dynvars%energy%deltU_fep(2)=energy_tgt%total-energy_ref%total

      end if

      ! return to the original state
      enefunc%lambljA      = alchemy%lambljA(lambid)
      enefunc%lambljB      = alchemy%lambljB(lambid)
      enefunc%lambelA      = alchemy%lambelA(lambid)
      enefunc%lambelB      = alchemy%lambelB(lambid)
      enefunc%lambbondA    = alchemy%lambbondA(lambid)
      enefunc%lambbondB    = alchemy%lambbondB(lambid)
      enefunc%lambrest     = alchemy%lambrest(lambid)
      call assign_lambda(alchemy, domain, enefunc)

    end if

  end subroutine compute_fep_energy

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_ref(domain, enefunc, boundary, energy, energy_ref)

    ! formal arguments
    type(s_domain),   target, intent(inout) :: domain
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_boundary),         intent(in)    :: boundary
    type(s_energy),           intent(inout) :: energy
    type(s_energy),           intent(inout) :: energy_ref

    ! local variables
    real(wip)                 :: volume
    real(dp)                  :: total

    total = energy%bond               + &
            energy%angle              + &
            energy%urey_bradley       + &
            energy%dihedral           + &
            energy%improper           + &
            energy%cmap               + &
            energy%electrostatic      + &
            energy%van_der_waals      + &
            energy%electric_field     + &
            energy%restraint_position
#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(total, energy_ref%total, 1, mpi_real8,  &
                       mpi_sum, mpi_comm_country, ierror)
#endif

    energy_ref%total = energy_ref%total          + &
                       energy%restraint_distance + &
                       energy%restraint_rmsd     + &
                       energy%restraint_emfit    + &
                       energy%contact            + &
                       energy%noncontact

    ! Dispersion correction
    if (enefunc%dispersion_corr /= Disp_corr_NONE) then
      volume =  boundary%box_size_x_ref * &
                boundary%box_size_y_ref * &
                boundary%box_size_z_ref
      energy_ref%disp_corr_energy = &
        enefunc%dispersion_energy_CC + &
        enefunc%dispersion_energy_AC * enefunc%lambljA + &
        enefunc%dispersion_energy_BC * enefunc%lambljB + &
        enefunc%dispersion_energy_AA * enefunc%lambljA * enefunc%lambljA + &
        enefunc%dispersion_energy_BB * enefunc%lambljB * enefunc%lambljB + &
        enefunc%dispersion_energy_AB * enefunc%lambljA * enefunc%lambljB
      energy_ref%disp_corr_energy = energy_ref%disp_corr_energy / volume
      energy_ref%total = energy_ref%total + energy_ref%disp_corr_energy
    end if

    return

  end subroutine compute_energy_ref

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_tgt(domain, enefunc, pairlist, boundary, energy)

    ! formal arguments
    type(s_domain),   target, intent(inout) :: domain
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_pairlist),         intent(in)    :: pairlist
    type(s_boundary),         intent(in)    :: boundary
    type(s_energy),           intent(inout) :: energy

    integer                           :: nc
    real(wip),      allocatable, save :: force(:,:,:)
    real(wip),      allocatable, save :: force_long(:,:,:)
    real(dp)                          :: virial(3,3)
    real(dp)                          :: virial_long(3,3)
    real(dp)                          :: virial_extern(3,3)


    if ( .not. allocated(force) ) then
      nc = size(domain%force(3,MaxAtom,:))
      allocate(force(3,MaxAtom,nc), force_long(3,MaxAtom,nc))
    end if

    call compute_energy(domain, enefunc, pairlist, boundary, &
                        domain%coord, .false., &
                        .true., .true., .true., .false., &
                        energy, domain%atmcls_pbc, domain%translated, &
                        force, force_long, &
                        domain%force_omp, domain%force_pbc, &
                        domain%virial_cellpair, virial, &
                        virial_long, virial_extern)

    energy%total = energy%total + &
               energy%restraint_distance + &
               energy%restraint_rmsd     + &
               energy%restraint_emfit    + &
               energy%contact            + &
               energy%noncontact

    if (enefunc%dispersion_corr /= Disp_Corr_NONE) then
      energy%total = energy%total + energy%disp_corr_energy
    end if

    deallocate(force, force_long)

    return

  end subroutine compute_energy_tgt

  !=============================================================================

  subroutine assign_condition_fep(parmsetid, alchemy, domain, ensemble, &
                                  enefunc, remd)

    ! formal arguments
    integer,                 intent(in)    :: parmsetid
    type(s_alchemy),         intent(in)    :: alchemy
    type(s_domain),          intent(inout) :: domain
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_remd),  optional, intent(inout) :: remd

    ! local variables
    integer                  :: i, j, k
    integer                  :: umbrid, funcid, atmid
    integer                  :: lambid
    real(wp)                 :: rest_param

    do i = 1, remd%dimension
      lambid = remd%iparameters(i,remd%parmidsets(parmsetid,i))
      enefunc%lambljA   = remd%dlambljA(lambid)
      enefunc%lambljB   = remd%dlambljB(lambid)
      enefunc%lambelA   = remd%dlambelA(lambid)
      enefunc%lambelB   = remd%dlambelB(lambid)
      enefunc%lambbondA = remd%dlambbondA(lambid)
      enefunc%lambbondB = remd%dlambbondB(lambid)
      enefunc%lambrest  = remd%dlambrest(lambid)

      call assign_lambda(alchemy, domain, enefunc)

      if (remd%types(i) == RemdAlchemyRest) then

        if (ensemble%ensemble == EnsembleNVE ) then
          call error_msg('Fep_Rest> Fep_Rest is not available for NVE!')
        end if

        rest_param = remd%dparameters(i,remd%parmidsets(parmsetid,i))
        call assign_condition_feprest( remd%fep_rest(i), &
                                                rest_param, domain, &
                                                ensemble, enefunc )

      end if

    end do

    return

  end subroutine assign_condition_fep

  !=============================================================================

  subroutine assign_condition_feprest( soltemp, rest_param, &
                                                domain, ensemble, enefunc )

    ! formal arguments
    type(s_soltemp),         intent(inout) :: soltemp
    real(wp),                intent(in)    :: rest_param
    type(s_domain),          intent(inout) :: domain
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_enefunc),         intent(inout) :: enefunc

    integer  :: alist(8), i, j, ix, num, tgt, n, ig, org, ncell, counter
    real(wp) :: coeff_full, coeff_half, wgt, tmp1, tmp2
    real(wp) :: el_fact, alpha

    coeff_full = ensemble%temperature / rest_param
    coeff_half = sqrt( coeff_full )

    if ( soltemp%done_setup ) then
      tmp1 = coeff_full
      tmp2 = coeff_half

      coeff_full = coeff_full / soltemp%rest_param_full
      coeff_half = coeff_half / soltemp%rest_param_half

      ! remember coeffs for next exchange
      !
      soltemp%rest_param_full = tmp1
      soltemp%rest_param_half = tmp2
    else
      soltemp%rest_param_full = coeff_full
      soltemp%rest_param_half = coeff_half
      soltemp%done_setup = .true.
    end if

    ncell = domain%num_cell_local + domain%num_cell_boundary

    if ( soltemp%sw_charge ) then
      do i = 1, ncell
        ! charge
        do ix = 1, domain%num_atom(i)
          if ( soltemp%is_solute(domain%id_l2g(ix,i)) > 0 ) then
            domain%charge(ix,i) = coeff_half * domain%charge(ix,i)
          end if
        end do
      end do
#ifdef USE_GPU
    if (domain%nonbond_kernel == NBK_GPU) then
      n = domain%num_atom_domain
      do i = 1, ncell
        ig = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)
          domain%translated(    ig+ix,1,1) = domain%coord(1,ix,i) + domain%trans_vec(1,ix,i)
          domain%translated(  n+ig+ix,1,1) = domain%coord(2,ix,i) + domain%trans_vec(2,ix,i)
          domain%translated(2*n+ig+ix,1,1) = domain%coord(3,ix,i) + domain%trans_vec(3,ix,i)
          domain%translated(3*n+ig+ix,1,1) = domain%charge(ix,i)
        end do
      end do
      call gpu_upload_charge( domain%translated );
    end if
#endif /* USE_GPU */

      ! need to update PME self energy (not necessary for NPT?)
      u_self = 0.0_wip
      el_fact = ELECOEF / enefunc%dielec_const
      alpha = enefunc%pme_alpha
      do i = 1, domain%num_cell_local
        do ix = 1, domain%num_atom(i)
          u_self = u_self + domain%charge(ix,i)**2
        end do
      end do
      u_self = - u_self * el_fact * alpha / sqrt(PI)

    end if

    do i = 1, domain%num_cell_local

      ! bond
      if ( soltemp%sw_bonds ) then
        call assign_condition_feprest_internal( soltemp, &
          2, coeff_full, &
          enefunc%num_bond(i), &
          enefunc%bond_list(:,:,i), &
          enefunc%bond_force_const(:,i) )
      end if

      ! do not add sw_angles flags here!
      call assign_condition_feprest_angle_ub( soltemp, coeff_full, &
        enefunc%num_angle(i), &
        enefunc%angle_list(:,:,i), &
        enefunc%angle_force_const(:,i), &
        enefunc%urey_force_const(:,i) )

      ! dihedral
      if ( soltemp%sw_dihedrals ) then
        call assign_condition_feprest_internal( soltemp, &
          4, coeff_full, &
          enefunc%num_dihedral(i), &
          enefunc%dihe_list(:,:,i), &
          enefunc%dihe_force_const(:,i) )
      end if

      ! impropers
      if ( soltemp%sw_impropers ) then
        call assign_condition_feprest_internal( soltemp, &
          4, coeff_full, &
          enefunc%num_improper(i), &
          enefunc%impr_list(:,:,i), &
          enefunc%impr_force_const(:,i) )
      end if

    end do

    ! lj
    if ( soltemp%sw_lj ) then
      n = enefunc%num_atom_cls
      if ( allocated(enefunc%nb14_lj6) ) then
        call assign_condition_feprest_lj( &
                 n, soltemp%rest_param_half, soltemp%rest_param_full, &
                 soltemp, enefunc%nb14_lj6 )
      end if
      if ( allocated(enefunc%nb14_lj12) ) then
        call assign_condition_feprest_lj( &
                 n, soltemp%rest_param_half, soltemp%rest_param_full, &
                 soltemp, enefunc%nb14_lj12 )
      end if
      if ( allocated(enefunc%nonb_lj6) ) then
        call assign_condition_feprest_lj( &
                 n, soltemp%rest_param_half, soltemp%rest_param_full, &
                 soltemp, enefunc%nonb_lj6 )
      end if
      if ( allocated(enefunc%nonb_lj12) ) then
        call assign_condition_feprest_lj( &
                 n, soltemp%rest_param_half, soltemp%rest_param_full, &
                 soltemp, enefunc%nonb_lj12 )
      end if
#ifdef USE_GPU
      if ( allocated(enefunc%nonb_lj12) .or. allocated(enefunc%nonb_lj6) ) then
        call gpu_upload_lj_coeffs( n, enefunc%nonb_lj12, enefunc%nonb_lj6 )
      end if
#endif /* USE_GPU */
    end if

    ! cmap
    if ( soltemp%num_cmap_type > 0 .and. soltemp%sw_cmaps ) then
      do i = 1, soltemp%num_cmap_type
        tgt = i + soltemp%istart_cmap_type - 1
        org = soltemp%cmap_type_org(i)
        wgt = soltemp%rest_param_full ** soltemp%cmap_weight(i)
        enefunc%cmap_coef(1:4,1:4,1:24,1:24,tgt) = &
          wgt * enefunc%cmap_coef(1:4,1:4,1:24,1:24,org)
      end do
    end if

    return

  end subroutine assign_condition_feprest

  !=============================================================================

  subroutine assign_condition_feprest_angle_ub( soltemp, &
                                  coeff_full, nangles, aindex, fc, fc_ub )

    ! formal arguments
    type(s_soltemp),         intent(in)    :: soltemp
    real(wp),                intent(in)    :: coeff_full
    integer,                 intent(in)    :: nangles
    integer,                 intent(in)    :: aindex(:,:)
    real(wp),                intent(inout) :: fc(:)
    real(wp),                intent(inout) :: fc_ub(:)

    ! local variables
    integer  :: alist(3), ix, j, num
    real(wp) :: wgt

    ! angle and urey
    do ix = 1, nangles
      alist(1:3) = aindex(1:3,ix)
      num = 0
      do j = 1, 3
        if (soltemp%is_solute(alist(j)) > 0) num = num + 1
      end do
      if ( num > 0 .and. soltemp%sw_angles ) then
        wgt = real(num,wp) / 3.0_wp
        wgt = coeff_full ** wgt
        fc(ix) = wgt * fc(ix)
      end if
      if ( fc_ub(ix) > EPS .and. soltemp%sw_ureys ) then
        num = 0
        if (soltemp%is_solute(alist(1)) > 0) num = num + 1
        if (soltemp%is_solute(alist(3)) > 0) num = num + 1
        if ( num > 0 ) then
          wgt = real(num,wp) / 2.0_wp
          wgt = coeff_full ** wgt
          fc_ub(ix) = wgt * fc_ub(ix)
        end if
      end if
    end do

    return

  end subroutine assign_condition_feprest_angle_ub

  !=============================================================================

  subroutine assign_condition_feprest_lj( n, coeff_half, coeff_full, &
                                                   soltemp, nbcoeff )

    ! formal arguments
    integer,                 intent(in)    :: n
    real(wp),                intent(in)    :: coeff_half
    real(wp),                intent(in)    :: coeff_full
    type(s_soltemp),         intent(in)    :: soltemp
    real(wp),                intent(inout) :: nbcoeff(n,n)

    ! local variables
    integer :: i, j, oldcount, newcount, org, orgj

    oldcount = soltemp%istart_atom_cls - 1
    newcount = oldcount + soltemp%num_atom_cls

    do i = 1, soltemp%num_atom_cls
      org = soltemp%atom_cls_no_org(i)
      do j = 1, oldcount
        nbcoeff(j,i+oldcount) = coeff_half * nbcoeff(j,org)
        nbcoeff(i+oldcount,j) = coeff_half * nbcoeff(j,org)
      end do
      do j = oldcount + 1, newcount
        orgj = soltemp%atom_cls_no_org(j-oldcount)
        nbcoeff(j,i+oldcount) = coeff_full * nbcoeff(orgj,org)
        nbcoeff(i+oldcount,j) = coeff_full * nbcoeff(orgj,org)
      end do
    end do

    return

  end subroutine assign_condition_feprest_lj

  !=============================================================================

  subroutine assign_condition_feprest_internal( soltemp, n, &
                                  coeff_full, n_internal, aindex, fc )

    ! formal arguments
    type(s_soltemp),         intent(in)    :: soltemp
    integer,                 intent(in)    :: n
    real(wp),                intent(in)    :: coeff_full
    integer,                 intent(in)    :: n_internal
    integer,                 intent(in)    :: aindex(:,:)
    real(wp),                intent(inout) :: fc(:)

    ! local variables
    integer  :: alist(1:n), ix, j, num
    real(wp) :: wgt

    do ix = 1, n_internal
      alist(1:n) = aindex(1:n,ix)
      num = 0
      do j = 1, n
        if (soltemp%is_solute(alist(j)) > 0) num = num + 1
      end do
      if ( num > 0 ) then
        wgt = real(num,wp) / real(n,wp)
        wgt = coeff_full ** wgt
        fc(ix) = wgt * fc(ix)
      end if
    end do

    return

  end subroutine assign_condition_feprest_internal

  !=============================================================================

  subroutine assign_lambda(alchemy, domain, enefunc)

    ! formal arguments
    type(s_alchemy),         intent(in)    :: alchemy
    type(s_domain),  target, intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: ix, icel, n, ig
    integer                  :: umbrid, funcid, atmid
    integer                  :: lambid
    integer                  :: iA, iB, ixA, ixB
    integer                  :: ncell_local, ncell
    real(wp)                 :: lambelA, lambelB, lambljA, lambljB
    real(wp)                 :: el_fact, alpha
    integer,         pointer :: fepgrp(:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: fep_chargeA(:,:), fep_chargeB(:,:)


    fepgrp        => domain%fepgrp
    charge        => domain%charge
    fep_chargeA   => domain%fep_chargeA
    fep_chargeB   => domain%fep_chargeB

    lambelA = enefunc%lambelA
    lambelB = enefunc%lambelB
    lambljA = enefunc%lambljA
    lambljB = enefunc%lambljB

    ncell_local = domain%num_cell_local
    ncell       = ncell_local + domain%num_cell_boundary
 
    ! elec
    do icel = 1, ncell
      do ix = 1, domain%num_atom(icel)
        if(enefunc%fep_topology == 2) then
          if (fepgrp(ix,icel) == 3) then
            charge(ix,icel) = lambelA * fep_chargeA(ix,icel)
          else if (fepgrp(ix,icel) == 4) then
            charge(ix,icel) = lambelB * fep_chargeB(ix,icel)
          end if
        else
          ! If hybrid topology is used, charges in singleA and singleB are merged
          ! only into singleA part. Charges in singleB part is set to 0.
          ! This modification can remove electrostatic interacitons between atoms 
          ! whose distances are 0.
          if (fepgrp(ix,icel) == 1) then
            charge(ix,icel) = lambelA*fep_chargeA(ix,icel) + lambelB*fep_chargeB(ix,icel)
          else if (fepgrp(ix,icel) == 3) then
            charge(ix,icel) = lambelA * fep_chargeA(ix,icel)
          else if (fepgrp(ix,icel) == 4) then
            charge(ix,icel) = lambelB * fep_chargeB(ix,icel)
          end if
        end if
      end do
    end do

#ifdef USE_GPU
    if (domain%nonbond_kernel == NBK_GPU) then
      n = domain%num_atom_domain
      do icel = 1, ncell
        ig = domain%start_atom(icel)
        do ix = 1, domain%num_atom(icel)
          domain%translated(    ig+ix,1,1) = domain%coord(1,ix,icel) + domain%trans_vec(1,ix,icel)
          domain%translated(  n+ig+ix,1,1) = domain%coord(2,ix,icel) + domain%trans_vec(2,ix,icel)
          domain%translated(2*n+ig+ix,1,1) = domain%coord(3,ix,icel) + domain%trans_vec(3,ix,icel)
          domain%translated(3*n+ig+ix,1,1) = charge(ix,icel)
        end do
      end do
      call gpu_upload_charge( domain%translated );
    end if
#endif /* USE_GPU */

    ! elec, self
    u_self = 0.0_wip
    el_fact = ELECOEF / enefunc%dielec_const
    alpha = enefunc%pme_alpha
    do icel = 1, ncell_local
      do ix = 1, domain%num_atom(icel)
        u_self = u_self + domain%charge(ix,icel)**2
      end do
    end do
    u_self = - u_self * el_fact * alpha / sqrt(PI)

    ! lj
    if ( allocated(enefunc%nb14_lj6) ) then
      call assign_lambda_lj(enefunc%num_atom_cls, lambljA, lambljB, &
        alchemy, enefunc%nb14_lj6)
    end if
    if ( allocated(enefunc%nb14_lj12) ) then
      call assign_lambda_lj(enefunc%num_atom_cls, lambljA, lambljB, &
        alchemy, enefunc%nb14_lj12)
    end if
    if ( allocated(enefunc%nonb_lj6) ) then
      call assign_lambda_lj(enefunc%num_atom_cls, lambljA, lambljB, &
        alchemy, enefunc%nonb_lj6)
    end if
    if ( allocated(enefunc%nonb_lj12) ) then
      call assign_lambda_lj(enefunc%num_atom_cls, lambljA, lambljB, &
        alchemy, enefunc%nonb_lj12)
    end if
#ifdef USE_GPU
    if ( allocated(enefunc%nonb_lj12) .or. allocated(enefunc%nonb_lj6) ) then
      call gpu_upload_lj_coeffs( enefunc%num_atom_cls, enefunc%nonb_lj12, enefunc%nonb_lj6 )
    end if
#endif /* USE_GPU */

    ! Make table of lambda and softcore in FEP
    call set_lambda_table_fep(enefunc)

    return

  end subroutine assign_lambda

  !=============================================================================

  subroutine assign_lambda_lj(n, lambljA, lambljB, alchemy, nbcoeff)

    ! formal arguments
    integer,                 intent(in)    :: n
    real(wp),                intent(in)    :: lambljA, lambljB
    type(s_alchemy),         intent(in)    :: alchemy
    real(wp),                intent(inout) :: nbcoeff(n,n)

    ! local variables
    integer :: i, j, orgi, orgiA, orgiB, orgj, orgjA, orgjB
    integer :: localcount_single, localcount_dualA, localcount_dualB
    integer :: oldcount, tmpcount, tmpcount2, newcount

    localcount_single = alchemy%num_atom_cls_single
    localcount_dualA  = alchemy%num_atom_cls_dualA
    localcount_dualB  = alchemy%num_atom_cls_dualB
    oldcount  = alchemy%istart_atom_cls_single - 1
    tmpcount  = oldcount + localcount_single
    tmpcount2 = oldcount + localcount_single + localcount_dualA
    newcount  = oldcount + localcount_single + localcount_dualA + localcount_dualB

    do i = 1, localcount_single
      orgiA = alchemy%atom_cls_no_org_singleA(i)
      orgiB = alchemy%atom_cls_no_org_singleB(i)
      ! single and Common
      do j = 1, oldcount
        nbcoeff(j,i+oldcount) = lambljA*nbcoeff(j,orgiA) + lambljB*nbcoeff(j,orgiB)
        nbcoeff(i+oldcount,j) = lambljA*nbcoeff(j,orgiA) + lambljB*nbcoeff(j,orgiB)
      end do
      ! single and single
      do j = oldcount + 1, tmpcount
        orgjA = alchemy%atom_cls_no_org_singleA(j-oldcount)
        orgjB = alchemy%atom_cls_no_org_singleB(j-oldcount)
        nbcoeff(j,i+oldcount) =   lambljA*lambljA*nbcoeff(orgjA,orgiA) &
                                + lambljA*lambljB*nbcoeff(orgjA,orgiB) &
                                + lambljB*lambljA*nbcoeff(orgjB,orgiA) &
                                + lambljB*lambljB*nbcoeff(orgjB,orgiB)
        nbcoeff(i+oldcount,j) =   lambljA*lambljA*nbcoeff(orgjA,orgiA) &
                                + lambljA*lambljB*nbcoeff(orgjA,orgiB) &
                                + lambljB*lambljA*nbcoeff(orgjB,orgiA) &
                                + lambljB*lambljB*nbcoeff(orgjB,orgiB)
      end do
      ! single and dualA
      do j = tmpcount + 1, tmpcount2
        orgj = alchemy%atom_cls_no_org_dualA(j-tmpcount)
        nbcoeff(j,i+oldcount) = lambljA*(lambljA*nbcoeff(orgj,orgiA) + lambljB*nbcoeff(orgj,orgiB))
        nbcoeff(i+oldcount,j) = lambljA*(lambljA*nbcoeff(orgj,orgiA) + lambljB*nbcoeff(orgj,orgiB))
      end do
      ! single and dualB
      do j = tmpcount2 + 1, newcount
        orgj = alchemy%atom_cls_no_org_dualB(j-tmpcount2)
        nbcoeff(j,i+oldcount) = lambljB*(lambljA*nbcoeff(orgj,orgiA) + lambljB*nbcoeff(orgj,orgiB))
        nbcoeff(i+oldcount,j) = lambljB*(lambljA*nbcoeff(orgj,orgiA) + lambljB*nbcoeff(orgj,orgiB))
      end do
    end do

    do i = 1, localcount_dualA
      orgi = alchemy%atom_cls_no_org_dualA(i)
      ! dualA and Common
      do j = 1, oldcount
        nbcoeff(j,i+tmpcount) = lambljA*nbcoeff(j,orgi)
        nbcoeff(i+tmpcount,j) = lambljA*nbcoeff(j,orgi)
      end do
      ! dualA and dualA
      do j = tmpcount + 1, tmpcount2
        orgj = alchemy%atom_cls_no_org_dualA(j-tmpcount)
        nbcoeff(j,i+tmpcount) = lambljA*lambljA*nbcoeff(orgj,orgi)
        nbcoeff(i+tmpcount,j) = lambljA*lambljA*nbcoeff(orgj,orgi)
      end do
      ! dualA and dualB
      do j = tmpcount2 + 1, newcount
        orgj = alchemy%atom_cls_no_org_dualB(j-tmpcount2)
        nbcoeff(j,i+tmpcount) = lambljA*lambljB*nbcoeff(orgj,orgi)
        nbcoeff(i+tmpcount,j) = lambljA*lambljB*nbcoeff(orgj,orgi)
      end do
    end do

    do i = 1, localcount_dualB
      orgi = alchemy%atom_cls_no_org_dualB(i)
      ! dualB and Common
      do j = 1, oldcount
        nbcoeff(j,i+tmpcount2) = lambljB*nbcoeff(j,orgi)
        nbcoeff(i+tmpcount2,j) = lambljB*nbcoeff(j,orgi)
      end do
      ! dualB and dualB
      do j = tmpcount2 + 1, newcount
        orgj = alchemy%atom_cls_no_org_dualB(j-tmpcount2)
        nbcoeff(j,i+tmpcount2) = lambljB*lambljB*nbcoeff(orgj,orgi)
        nbcoeff(i+tmpcount2,j) = lambljB*lambljB*nbcoeff(orgj,orgi)
      end do
    end do

    return

  end subroutine assign_lambda_lj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    set_lambda_table_fep
  !> @brief        make lambda table for calculations of bonded and nonbonded
  !                energies in FEP
  !! @authors      HO
  !! @param[in]    enefunc  : enefunc information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine set_lambda_table_fep(enefunc)

    ! formal arguments
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: fg1, fg2
    integer                  :: i1, i2, i3, i4
    integer                  :: i5, i6, i7, i8
    integer                  :: idx
    integer,         pointer :: fepgrp_nonb(:,:)
    integer(1),      pointer :: fep_mask(:,:)
    real(wp),        pointer :: table_bond_lambda(:,:)
    real(wp),        pointer :: table_angl_lambda(:,:,:)
    real(wp),        pointer :: table_dihe_lambda(:,:,:,:)

    real(wp),        pointer :: table_sclj(:,:)
    real(wp),        pointer :: table_scel(:,:)

    fepgrp_nonb       => enefunc%fepgrp_nonb
    fep_mask          => enefunc%fep_mask
    table_bond_lambda => enefunc%table_bond_lambda
    table_angl_lambda => enefunc%table_angl_lambda
    table_dihe_lambda => enefunc%table_dihe_lambda

    table_sclj => enefunc%table_sclj
    table_scel => enefunc%table_scel

    ! Make lambda table for calculation of nonbonded energy
    !
    fepgrp_nonb(:,:) = 0
    fep_mask(:,:)    = 0

    table_sclj(:,:) = 0.0_wp
    table_scel(:,:) = 0.0_wp
    do fg1 = 1, 5
      do fg2 = 1, 5

        ! Make table of lambda and softcore in FEP
        !
        if (((fg1==1).or.(fg1==2).or.(fg1==5)) .and. &
          ((fg2==1).or.(fg2==2).or.(fg2==5))) then
          fepgrp_nonb(fg1,fg2) = 5
        else if ((fg1==3).or.(fg2==3)) then
          fepgrp_nonb(fg1,fg2) = 3
          table_sclj(fg1,fg2) = enefunc%sc_alpha*(1.0_wp-enefunc%lambljA)
          table_scel(fg1,fg2) = enefunc%sc_beta*(1.0_wp-enefunc%lambelA)
        else if ((fg1==4).or.(fg2==4)) then
          fepgrp_nonb(fg1,fg2) = 4
          table_sclj(fg1,fg2) = enefunc%sc_alpha*(1.0_wp-enefunc%lambljB)
          table_scel(fg1,fg2) = enefunc%sc_beta*(1.0_wp-enefunc%lambelB)
        end if

        ! Interactons between partA and partB are set to 0,
        ! while others are set to 1.
        !
        if (((fg1==3).and.(fg2==4)) .or. ((fg1==4).and.(fg2==3))) then
          fepgrp_nonb(fg1,fg2) = 0
          table_sclj(fg1,fg2) = 0.0_wp
          table_scel(fg1,fg2) = 0.0_wp
          fep_mask(fg1,fg2) = 0
        else
          fep_mask(fg1,fg2) = 1
        end if

      end do
    end do
#ifdef USE_GPU
    call gpu_upload_table_fep( fep_mask, table_sclj, table_scel )
#endif /* USE_GPU */


    ! Make lambda table for calculation of bonded energy
    !
    table_bond_lambda(:,:) = 0.0_wp
    table_angl_lambda(:,:,:) = 0.0_wp
    table_dihe_lambda(:,:,:,:) = 0.0_wp

    ! Bond table
    !
    do i2 = 1, 5
      do i1 = 1, 5
        ! preserve: all 5
        if ((i1 == 5) .and. (i2 == 5)) then
          table_bond_lambda(i1,i2) = 1.0_wp
        end if
        ! single: all 1, 1 and 5
        if (((i1 == 1) .or. (i1 == 5)) .and. &
          ((i2 == 1) .or. (i2 == 5))) then
          if ((i1 == 1) .or. (i2 == 1)) then
            table_bond_lambda(i1,i2) = enefunc%lambbondA
          end if
        end if
        ! dualA: all 3, 3 and 5, 3 and 1
        if (((i1 == 3) .or. (i1 == 5) .or. (i1 == 1)) .and. &
          ((i2 == 3) .or. (i2 == 5) .or. (i2 == 1))) then
          if ((i1 == 3) .or. (i2 == 3)) then
            table_bond_lambda(i1,i2) = 1.0_wp
          end if
        end if
        ! dualB: all 4, 4 and 5, 4 and 1
        if (((i1 == 4) .or. (i1 == 5) .or. (i1 == 1)) .and. &
          ((i2 == 4) .or. (i2 == 5) .or. (i2 == 1))) then
          if ((i1 == 4) .or. (i2 == 4)) then
            table_bond_lambda(i1,i2) = 1.0_wp
          end if
        end if
      end do
    end do

    ! Angle table
    !
    do i3 = 1, 5
      do i2 = 1, 5
        do i1 = 1, 5
          ! preserve: all 5
          if ((i1 == 5) .and. (i2 == 5) .and. (i3 == 5)) then
            table_angl_lambda(i1,i2,i3) = 1.0_wp
          end if
          ! single: all 1, 1 and 5
          if (((i1 == 1) .or. (i1 == 5)) .and. &
            ((i2 == 1) .or. (i2 == 5)) .and. &
            ((i3 == 1) .or. (i3 == 5))) then
            if ((i1 == 1) .or. (i2 == 1) .or. (i3 == 1)) then
              table_angl_lambda(i1,i2,i3) = enefunc%lambbondA
            end if
          end if
          ! dualA: all 3, 3 and 5, 3 and 1
          if (((i1 == 3) .or. (i1 == 5) .or. (i1 == 1)) .and. &
            ((i2 == 3) .or. (i2 == 5) .or. (i2 == 1)) .and. &
            ((i3 == 3) .or. (i3 == 5) .or. (i3 == 1))) then
            if ((i1 == 3) .or. (i2 == 3) .or. (i3 == 3)) then
              table_angl_lambda(i1,i2,i3) = 1.0_wp
            end if
          end if
          ! dualB: all 4, 4 and 5, 4 and 1
          if (((i1 == 4) .or. (i1 == 5) .or. (i1 == 1)) .and. &
            ((i2 == 4) .or. (i2 == 5) .or. (i2 == 1)) .and. &
            ((i3 == 4) .or. (i3 == 5) .or. (i3 == 1))) then
            if ((i1 == 4) .or. (i2 == 4) .or. (i3 == 4)) then
              table_angl_lambda(i1,i2,i3) = 1.0_wp
            end if
          end if
        end do
      end do
    end do

    ! Dihedral table
    !
    do i4 = 1, 5
      do i3 = 1, 5
        do i2 = 1, 5
          do i1 = 1, 5
            ! preserve: all 5
            if ((i1 == 5) .and. (i2 == 5) .and. (i3 == 5) .and. (i4 == 5)) then
              table_dihe_lambda(i1,i2,i3,i4) = 1.0_wp
            end if
            ! single: all 1, 1 and 5
            if (((i1 == 1) .or. (i1 == 5)) .and. &
              ((i2 == 1) .or. (i2 == 5)) .and. &
              ((i3 == 1) .or. (i3 == 5)) .and. &
              ((i4 == 1) .or. (i4 == 5))) then
              if ((i1 == 1) .or. (i2 == 1) .or. (i3 == 1) .or. (i4 == 1)) then
                table_dihe_lambda(i1,i2,i3,i4) = enefunc%lambbondA
              end if
            end if
            ! dualA: all 3, 3 and 5, 3 and 1
            if (((i1 == 3) .or. (i1 == 5) .or. (i1 == 1)) .and. &
              ((i2 == 3) .or. (i2 == 5) .or. (i2 == 1)) .and. &
              ((i3 == 3) .or. (i3 == 5) .or. (i3 == 1)) .and. &
              ((i4 == 3) .or. (i4 == 5) .or. (i4 == 1))) then
              if ((i1 == 3) .or. (i2 == 3) .or. (i3 == 3) .or. (i4 == 3)) then
                table_dihe_lambda(i1,i2,i3,i4) = 1.0_wp
              end if
            end if
            ! dualB: all 4, 4 and 5, 4 and 1
            if (((i1 == 4) .or. (i1 == 5) .or. (i1 == 1)) .and. &
              ((i2 == 4) .or. (i2 == 5) .or. (i2 == 1)) .and. &
              ((i3 == 4) .or. (i3 == 5) .or. (i3 == 1)) .and. &
              ((i4 == 4) .or. (i4 == 5) .or. (i4 == 1))) then
              if ((i1 == 4) .or. (i2 == 4) .or. (i3 == 4) .or. (i4 == 4)) then
                table_dihe_lambda(i1,i2,i3,i4) = 1.0_wp
              end if
            end if

          end do
        end do
      end do
    end do

    return

  end subroutine set_lambda_table_fep

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine set_fepgrp_nonb(enefunc)

    ! formal arguments
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: fg1, fg2
    integer                  :: i1, i2, i3, i4
    integer                  :: i5, i6, i7, i8
    integer                  :: idx
    integer,         pointer :: fepgrp_nonb(:,:)

    fepgrp_nonb       => enefunc%fepgrp_nonb

    ! Make lambda table for calculation of nonbonded energy
    !
    fepgrp_nonb(:,:) = 0

    do fg1 = 1, 5
      do fg2 = 1, 5

        ! Make table of lambda and softcore in FEP
        !
        if (((fg1==1).or.(fg1==2).or.(fg1==5)) .and. &
          ((fg2==1).or.(fg2==2).or.(fg2==5))) then
          fepgrp_nonb(fg1,fg2) = 5
        else if ((fg1==3).or.(fg2==3)) then
          fepgrp_nonb(fg1,fg2) = 3
        else if ((fg1==4).or.(fg2==4)) then
          fepgrp_nonb(fg1,fg2) = 4
        end if

        ! Interactons between partA and partB are set to 0,
        ! while others are set to 1.
        !
        if (((fg1==3).and.(fg2==4)) .or. ((fg1==4).and.(fg2==3))) then
          fepgrp_nonb(fg1,fg2) = 0
        end if

      end do
    end do

    return

  end subroutine set_fepgrp_nonb

end module sp_fep_energy_mod
