!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_nonbond_fep_mod
!> @brief   calculate nonbonded energy of perturbed regions in FEP
!! @authors Hiraku Oshima (HO)
!  
!  (c) Copyright 2022 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_nonbond_fep_mod

  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use mpi_parallel_mod
  use timers_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  !
  ! for CPU
  public  :: compute_energy_nonbond_tbl_lnr_generic_fep
  public  :: compute_energy_nonbond_tbl_lnr_chk_generic_fep
  public  :: compute_force_nonbond_tbl_lnr_generic_fep
  public  :: compute_energy_nonbond_notbl_generic_fep
  public  :: compute_energy_nonbond_notbl_chk_generic_fep
  public  :: compute_force_nonbond_notbl_generic_fep
  ! for GPU
  public  :: compute_energy_nonbond_tbl_lnr_gpu_fep
  public  :: compute_force_nonbond_tbl_lnr_gpu_fep
  public  :: compute_energy_nonbond_notbl_gpu_fep
  public  :: compute_force_nonbond_notbl_gpu_fep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_tbl_lnr_generic_fep
  !> @brief        calculate nonbonded energy with lookup table for FEP
  !! @authors      HO
  !! @param[in]    domain    : domain information
  !! @param[in]    enefunc   : potential energy functions
  !! @param[in]    pairlist  : interaction list in each domain
  !! @param[inout] coord_pbc : pbc oriented coordinates
  !! @param[inout] force     : forces for each cell
  !! @param[inout] virial    : virial term of target systems
  !! @param[inout] eelec     : electrostatic energy of target systems
  !! @param[inout] evdw      : van der Waals energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_tbl_lnr_generic_fep( &
                                                 domain, enefunc, pairlist, &
                                                 virial,  &
                                                 eelec, evdw)

    ! formal arguments
    type(s_domain),   target, intent(inout) :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(dp),                 intent(inout) :: virial(:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: lj6, lj12
    real(wp)                  :: R
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: trans_x, trans_y, trans_z
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: elec_temp, evdw_temp
    real(wp)                  :: force_local(3)
    real(wp)                  :: density_tmp, qxy
    integer                   :: i, ix, iy, j, k, ij, iix, L, L1
    integer                   :: num_nb15, ini_nb15, fin_nb15
    integer                   :: id, omp_get_thread_num
    integer                   :: ncell, ncell_local
    integer                   :: check_virial
    integer                   :: iatmcls, jatmcls

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: atmcls(:,:)
    integer(int2),    pointer :: cell_pairlist(:,:)
    integer,          pointer :: natom(:)
    integer,          pointer :: nb15_list(:,:)
    integer,          pointer :: nb15_cell(:)
    integer,          pointer :: num_nb15_calc1(:,:), num_nb15_calc(:,:)
    integer,          pointer :: nb15_calc_list1(:,:), nb15_calc_list(:,:)
    integer(int1),    pointer :: virial_check(:,:)

    ! FEP
    real(wp)                  :: rij2_sclj, rij2_scel
    integer                   :: fg1, fg2
    integer,          pointer :: fepgrp(:,:)
    real(wp),         pointer :: table_sclj(:,:)
    real(wp),         pointer :: table_scel(:,:)
    real(wp),         pointer :: coord_pbc(:,:,:)
    real(wp),         pointer :: force(:,:,:,:)

    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff

    ! FEP
    nb15_list       => pairlist%nb15_list_fep
    nb15_cell       => pairlist%nb15_cell_fep
    num_nb15_calc1  => pairlist%num_nb15_calc1_fep
    num_nb15_calc   => pairlist%num_nb15_calc_fep
    nb15_calc_list1 => pairlist%nb15_calc_list1_fep
    nb15_calc_list  => pairlist%nb15_calc_list_fep
    fepgrp          => domain%fepgrp
    table_sclj      => enefunc%table_sclj
    table_scel      => enefunc%table_scel
    coord_pbc       => domain%translated_fep
    force           => domain%force_pbc_fep
    
    density_tmp     = cutoff2 * density

    ! calculate energy and gradient
    !

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, ix, num_nb15, rtmp, qtmp, ini_nb15, fin_nb15, k, iy,  &
    !$omp         rij2, L, R, term_lj12, term_lj6, term_elec, grad_coef, work, &
    !$omp         ij, j, trans_x, trans_y, trans_z, iix, force_local, lj6,     &
    !$omp         lj12, L1, elec_temp, evdw_temp, dij, check_virial,           &
    !$omp         fg1, fg2, rij2_sclj, rij2_scel, iatmcls, jatmcls, qxy)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread
      do ix = 1, natom(i)
        coord_pbc(ix,1,i) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(ix,2,i) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(ix,3,i) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    !$omp barrier

    ! energy within a cell
    !
    do i = id+1, ncell_local, nthread

      num_nb15 = 0
      do ix = 1, natom(i) - 1

        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp = charge(ix,i)
        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc1(ix,i)
        num_nb15  = fin_nb15
        force_local(1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp

        iatmcls = atmcls(ix,i)

        ! FEP: flag for atom (ix,i)
        fg1 = fepgrp(ix,i)

!dir$ simd
        do k = ini_nb15, fin_nb15

          iy = nb15_calc_list1(k,i)
          jatmcls = atmcls(iy,i)
          lj6  = nonb_lj6(iatmcls,jatmcls)
          lj12 = nonb_lj12(iatmcls,jatmcls)
          qxy = qtmp*charge(iy,i)

          ! compute distance
          !
          dij(1) = rtmp(1) - coord_pbc(iy,1,i)
          dij(2) = rtmp(2) - coord_pbc(iy,2,i)
          dij(3) = rtmp(3) - coord_pbc(iy,3,i)
          rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          ! FEP: flag for atom (iy,i)
          fg2 = fepgrp(iy,i)

          ! FEP: LJ with soft core
          rij2_sclj  = density_tmp/(rij2 + table_sclj(fg1,fg2))
          L    = int(rij2_sclj)
          R    = rij2_sclj - L
          L1   = 3*L - 2
          term_lj12 = table_ene(L1)   + R*(table_ene(L1+3)-table_ene(L1))
          term_lj6  = table_ene(L1+1) + R*(table_ene(L1+4)-table_ene(L1+1))
          evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
          term_lj12 = table_grad(L1)   + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))

          ! FEP: elec with soft core
          rij2_scel  = density_tmp/(rij2 + table_scel(fg1,fg2))
          L    = int(rij2_scel)
          R    = rij2_scel - L
          L1   = 3*L
          term_elec = table_ene(L1) + R*(table_ene(L1+3)-table_ene(L1))
          elec_temp = elec_temp + qxy*term_elec
          term_elec = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))

          grad_coef = term_lj12*lj12 - term_lj6*lj6 + qxy*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force(iy,1,i,id+1) = force(iy,1,i,id+1) + work(1)
          force(iy,2,i,id+1) = force(iy,2,i,id+1) + work(2)
          force(iy,3,i,id+1) = force(iy,3,i,id+1) + work(3)

        end do

        force(ix,1:3,i,id+1) = force(ix,1:3,i,id+1) + force_local(1:3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp

      end do

    end do


    ! interaction between different cells
    !
    do ij = id+1, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      num_nb15 = 0

      ! FEP: shift for generic & gpu kernels
      if (domain%nonbond_kernel /= NBK_Fugaku .and. &
          domain%nonbond_kernel /= NBK_Intel) then
        trans_x = real(cell_move(1,j,i),wp) * system_size(1)
        trans_y = real(cell_move(2,j,i),wp) * system_size(2)
        trans_z = real(cell_move(3,j,i),wp) * system_size(3)
      else
        trans_x = 0.0_wp
        trans_y = 0.0_wp
        trans_z = 0.0_wp
      end if

      check_virial = virial_check(j,i)

      do iix = 1, nb15_cell(ij)

        ix   = nb15_list(iix,ij)
        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp = charge(ix,i)

        ini_nb15 = num_nb15 + 1
        fin_nb15 = num_nb15 + num_nb15_calc(ix,ij)
        num_nb15 = fin_nb15

        force_local(1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp

        iatmcls = atmcls(ix,i)

        ! FEP: flag for atom (ix,i)
        fg1 = fepgrp(ix,i)

!dir$ simd
        do k = ini_nb15, fin_nb15
          iy   = nb15_calc_list(k,ij)
          jatmcls = atmcls(iy,j)
          lj6  = nonb_lj6 (iatmcls,jatmcls)
          lj12 = nonb_lj12(iatmcls,jatmcls)
          qxy = qtmp*charge(iy,j)

          ! compute distance
          !
          dij(1) = rtmp(1) - coord_pbc(iy,1,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(iy,2,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(iy,3,j) + trans_z
          rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          ! FEP: flag for atom (iy,j)
          fg2 = fepgrp(iy,j)

          ! FEP: LJ with soft core
          rij2_sclj = density_tmp/(rij2 + table_sclj(fg1,fg2))
          L    = int(rij2_sclj)
          R    = rij2_sclj - L
          L1   = 3*L - 2
          term_lj12 = table_ene(L1)   + R*(table_ene(L1+3)-table_ene(L1))
          term_lj6  = table_ene(L1+1) + R*(table_ene(L1+4)-table_ene(L1+1))
          evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
          term_lj12 = table_grad(L1)   + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))

          ! FEP: elec with soft core
          rij2_scel = density_tmp/(rij2 + table_scel(fg1,fg2))
          L    = int(rij2_scel)
          R    = rij2_scel - L
          L1   = 3*L
          term_elec = table_ene(L1) + R*(table_ene(L1+3)-table_ene(L1))
          elec_temp = elec_temp + qxy*term_elec
          term_elec = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))

          grad_coef = term_lj12*lj12 - term_lj6*lj6 + qxy*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)

          force(iy,1,j,id+1) = force(iy,1,j,id+1) + work(1)
          force(iy,2,j,id+1) = force(iy,2,j,id+1) + work(2)
          force(iy,3,j,id+1) = force(iy,3,j,id+1) + work(3)

        end do

        if (check_virial == 1) &
          virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
        force(ix,1:3,i,id+1) = force(ix,1:3,i,id+1) + force_local(1:3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond_tbl_lnr_generic_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_tbl_lnr_chk_generic_fep
  !> @brief        calculate nonbonded energy with lookup table for FEP
  !! @authors      HO
  !! @param[in]    domain    : domain information
  !! @param[in]    enefunc   : potential energy functions
  !! @param[in]    pairlist  : interaction list in each domain
  !! @param[inout] coord_pbc : pbc oriented coordinates
  !! @param[inout] force     : forces for each cell
  !! @param[inout] virial    : virial term of target systems
  !! @param[inout] eelec     : electrostatic energy of target systems
  !! @param[inout] evdw      : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_tbl_lnr_chk_generic_fep( &
                                                 domain, enefunc, pairlist, &
                                                 virial,  &
                                                 eelec, evdw)

    ! formal arguments
    type(s_domain),   target, intent(inout) :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(dp),                 intent(inout) :: virial(:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: lj6, lj12
    real(wp)                  :: R
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, term_lj, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: trans_x, trans_y, trans_z
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: elec_temp, evdw_temp
    real(wp)                  :: dij_list(3,MaxAtom)
    real(wp)                  :: force_local(3), rij2_list(MaxAtom)
    real(wp)                  :: force_localj(3,MaxAtom)
    real(wp)                  :: minimum_contact
    real(wp)                  :: density_tmp, qxy
    integer                   :: i, ix, iy, j, k, ij, iix, iwater, list, L, L1
    integer                   :: num_nb15, ini_nb15, fin_nb15, num_count
    integer                   :: id, omp_get_thread_num, j_list(MaxAtom)
    integer                   :: ncell, ncell_local
    integer                   :: check_virial
    integer                   :: iatmcls, jatmcls

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: atmcls(:,:)
    integer(int2),    pointer :: cell_pairlist(:,:)
    integer,          pointer :: natom(:)
    integer,          pointer :: nb15_list(:,:)
    integer,          pointer :: nb15_cell(:)
    integer,          pointer :: num_nb15_calc1(:,:), num_nb15_calc(:,:)
    integer,          pointer :: nb15_calc_list1(:,:), nb15_calc_list(:,:)
    integer(int1),    pointer :: virial_check(:,:)

    ! FEP
    real(wp)                  :: rij2_sclj, rij2_scel
    integer                   :: fg1, fg2
    integer,          pointer :: fepgrp(:,:)
    real(wp),         pointer :: table_sclj(:,:)
    real(wp),         pointer :: table_scel(:,:)
    real(wp),         pointer :: coord_pbc(:,:,:)
    real(wp),         pointer :: force(:,:,:,:)

    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff
    minimum_contact =  enefunc%minimum_contact

    ! FEP
    nb15_list       => pairlist%nb15_list_fep
    nb15_cell       => pairlist%nb15_cell_fep
    num_nb15_calc1  => pairlist%num_nb15_calc1_fep
    num_nb15_calc   => pairlist%num_nb15_calc_fep
    nb15_calc_list1 => pairlist%nb15_calc_list1_fep
    nb15_calc_list  => pairlist%nb15_calc_list_fep
    fepgrp          => domain%fepgrp
    table_sclj      => enefunc%table_sclj
    table_scel      => enefunc%table_scel
    coord_pbc       => domain%translated_fep
    force           => domain%force_pbc_fep

    density_tmp     = cutoff2 * density

    ! calculate energy and gradient
    !

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, ix, num_nb15, rtmp, qtmp, ini_nb15, fin_nb15, k, iy,  &
    !$omp         rij2, L, R, term_lj12, term_lj6, term_lj, grad_coef, work,   &
    !$omp         term_elec, iwater, ij, j, trans_x, trans_y,                  &
    !$omp         trans_z, iix, list, num_count, j_list, dij_list, rij2_list,  &
    !$omp         force_local, force_localj, lj6, lj12, L1, elec_temp,         &
    !$omp         evdw_temp, dij, check_virial, iatmcls, jatmcls, qxy,         &
    !$omp         fg1, fg2, rij2_sclj, rij2_scel)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread
      do ix = 1, natom(i)
        coord_pbc(ix,1,i) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(ix,2,i) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(ix,3,i) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    !$omp barrier

    ! energy within a cell
    !
    do i = id+1, ncell_local, nthread

      num_nb15 = 0
      do ix = 1, natom(i) - 1

        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp = charge(ix,i)
        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc1(ix,i)
        num_nb15  = fin_nb15
        num_count = 0

        do k = ini_nb15, fin_nb15

          iy = nb15_calc_list1(k,i)

          ! compute distance
          !
          dij(1) = rtmp(1) - coord_pbc(iy,1,i)
          dij(2) = rtmp(2) - coord_pbc(iy,2,i)
          dij(3) = rtmp(3) - coord_pbc(iy,3,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1:3,num_count)  = dij(1:3)
            rij2_list(num_count) = max(rij2, minimum_contact)
            j_list(num_count)   = iy
          end if
        end do

        force_local(1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp

        iatmcls = atmcls(ix,i)

        ! FEP: flag for atom (ix,i)
        fg1 = fepgrp(ix,i)

!ocl norecurrence
!ocl swp
        do k = 1, num_count
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          iy   = j_list(k)
          jatmcls = atmcls(iy,i)
          lj6  = nonb_lj6(iatmcls,jatmcls)
          lj12 = nonb_lj12(iatmcls,jatmcls)
          qxy = qtmp*charge(iy,i)

          ! FEP: flag for atom (iy,i)
          fg2 = fepgrp(iy,i)

          ! FEP: LJ with soft core
          rij2_sclj = density_tmp/(rij2_list(k) + table_sclj(fg1,fg2))
          L    = int(rij2_sclj)
          R    = rij2_sclj - L
          L1   = 3*L - 2
          term_lj12 = table_ene(L1) + R*(table_ene(L1+3)-table_ene(L1))
          term_lj6  = table_ene(L1+1) + R*(table_ene(L1+4)-table_ene(L1+1))
          evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
          term_lj12 = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))

          ! FEP: elec with soft core
          rij2_scel = density_tmp/(rij2_list(k) + table_scel(fg1,fg2))
          L    = int(rij2_scel)
          R    = rij2_scel - L
          L1   = 3*L
          term_elec = table_ene(L1) + R*(table_ene(L1+3)-table_ene(L1))
          elec_temp = elec_temp + qxy*term_elec
          term_elec = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))

          grad_coef = term_lj12*lj12 - term_lj6*lj6 + qxy*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force_localj(1,k)  = work(1)
          force_localj(2,k)  = work(2)
          force_localj(3,k)  = work(3)

        end do

        do k = 1, num_count
          iy = j_list(k)
          force(iy,1,i,id+1) = force(iy,1,i,id+1) + force_localj(1,k)
          force(iy,2,i,id+1) = force(iy,2,i,id+1) + force_localj(2,k)
          force(iy,3,i,id+1) = force(iy,3,i,id+1) + force_localj(3,k)
        end do
        force(ix,1:3,i,id+1) = force(ix,1:3,i,id+1) + force_local(1:3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp

      end do

    end do


    ! interaction between different cells
    !
    do ij = id+1, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      num_nb15 = 0

      ! FEP: shift for generic & gpu kernels
      if (domain%nonbond_kernel /= NBK_Fugaku .and. &
          domain%nonbond_kernel /= NBK_Intel) then
        trans_x = real(cell_move(1,j,i),wp) * system_size(1)
        trans_y = real(cell_move(2,j,i),wp) * system_size(2)
        trans_z = real(cell_move(3,j,i),wp) * system_size(3)
      else
        trans_x = 0.0_wp
        trans_y = 0.0_wp
        trans_z = 0.0_wp
      end if

      check_virial = virial_check(j,i)

      do iix = 1, nb15_cell(ij)

        ix   = nb15_list(iix,ij)
        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp = charge(ix,i)

        ini_nb15 = num_nb15 + 1
        fin_nb15 = num_nb15 + num_nb15_calc(ix,ij)
        num_nb15 = fin_nb15

        num_count = 0
        do k = ini_nb15, fin_nb15
          iy   = nb15_calc_list(k,ij)
          dij(1) = rtmp(1) - coord_pbc(iy,1,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(iy,2,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(iy,3,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1:3,num_count)  = dij(1:3)
            rij2_list(num_count) = max(rij2, minimum_contact)
            j_list(num_count)   = iy
          end if
        end do

        force_local(1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp

        iatmcls = atmcls(ix,i)

        ! FEP: flag for atom (ix,i)
        fg1 = fepgrp(ix,i)

!ocl norecurrence
!ocl swp
        do k = 1, num_count
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          iy   = j_list(k)
          jatmcls = atmcls(iy,j)
          lj6  = nonb_lj6 (iatmcls,jatmcls)
          lj12 = nonb_lj12(iatmcls,jatmcls)
          qxy = qtmp*charge(iy,j)

          ! FEP: flag for atom (iy,j)
          fg2 = fepgrp(iy,j)

          ! FEP: LJ with soft core
          rij2_sclj = density_tmp/(rij2_list(k) + table_sclj(fg1,fg2))
          L    = int(rij2_sclj)
          R    = rij2_sclj - L
          L1   = 3*L - 2
          term_lj12 = table_ene(L1) + R*(table_ene(L1+3)-table_ene(L1))
          term_lj6  = table_ene(L1+1) + R*(table_ene(L1+4)-table_ene(L1+1))
          evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
          term_lj12 = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))

          ! FEP: elec with soft core
          rij2_scel = density_tmp/(rij2_list(k) + table_scel(fg1,fg2))
          L    = int(rij2_scel)
          R    = rij2_scel - L
          L1   = 3*L
          term_elec = table_ene(L1) + R*(table_ene(L1+3)-table_ene(L1))
          elec_temp = elec_temp + qxy*term_elec
          term_elec = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))

          grad_coef = term_lj12*lj12 - term_lj6*lj6 + qxy*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)

          force_localj(1,k) = work(1)
          force_localj(2,k) = work(2)
          force_localj(3,k) = work(3)

        end do

        do k = 1, num_count
          iy = j_list(k)
          force(iy,1,j,id+1) = force(iy,1,j,id+1) + force_localj(1,k)
          force(iy,2,j,id+1) = force(iy,2,j,id+1) + force_localj(2,k)
          force(iy,3,j,id+1) = force(iy,3,j,id+1) + force_localj(3,k)
        end do
        if (check_virial == 1) &
          virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
        force(ix,1:3,i,id+1) = force(ix,1:3,i,id+1) + force_local(1:3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond_tbl_lnr_chk_generic_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_tbl_lnr_generic_fep
  !> @brief        calculate nonbonded force without solvents with lookup table
  !                for FEP
  !! @authors      HO
  !! @param[in]    domain    : domain information
  !! @param[in]    enefunc   : potential energy functions
  !! @param[in]    pairlist  : interaction list in each domain
  !! @param[inout] coord_pbc : pbc oriented coordinates
  !! @param[inout] force     : forces for each cell
  !! @param[inout] virial    : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_tbl_lnr_generic_fep( &
                                                 domain, enefunc, pairlist, &
                                                 virial)

    ! formal arguments
    type(s_domain),   target, intent(inout) :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(dp),                 intent(inout) :: virial(:,:)

    ! local variables
    real(wp)                  :: dij(1:3), rij2, inv_r2
    real(wp)                  :: R, lj12, lj6
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: trans_x, trans_y, trans_z
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: force_local(3)
    integer                   :: i, ix, iy, j, k, ij, iix, iwater, list, L, L1
    integer                   :: num_nb15, ini_nb15, fin_nb15, num_count
    integer                   :: id, omp_get_thread_num
    integer                   :: iatmcls,jatmcls
    integer                   :: ncell, ncell_local
    integer                   :: check_virial
    real(dp)                  :: sas,eae
    real(dp)                  :: density_tmp

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_grad(:)
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: atmcls(:,:)
    integer(int2),    pointer :: cell_pairlist(:,:)
    integer,          pointer :: natom(:)
    integer,          pointer :: nb15_list(:,:)
    integer,          pointer :: nb15_cell(:)
    integer,          pointer :: num_nb15_calc1(:,:), num_nb15_calc(:,:)
    integer,          pointer :: nb15_calc_list1(:,:), nb15_calc_list(:,:)
    integer(int1),    pointer :: virial_check(:,:)

    ! FEP
    real(wp)                  :: rij2_sclj, rij2_scel
    integer                   :: fg1, fg2
    integer,          pointer :: fepgrp(:,:)
    real(wp),         pointer :: table_sclj(:,:)
    real(wp),         pointer :: table_scel(:,:)
    real(wp),         pointer :: coord_pbc(:,:,:)
    real(wp),         pointer :: force(:,:,:,:)

    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff

    ! FEP
    nb15_list       => pairlist%nb15_list_fep
    nb15_cell       => pairlist%nb15_cell_fep
    num_nb15_calc1  => pairlist%num_nb15_calc1_fep
    num_nb15_calc   => pairlist%num_nb15_calc_fep
    nb15_calc_list1 => pairlist%nb15_calc_list1_fep
    nb15_calc_list  => pairlist%nb15_calc_list_fep
    fepgrp          => domain%fepgrp
    table_sclj      => enefunc%table_sclj
    table_scel      => enefunc%table_scel
    coord_pbc       => domain%translated_fep
    force           => domain%force_pbc_fep

    density_tmp     = cutoff2 * density

    ! calculate energy and gradient
    !
    !$omp parallel default(shared)                                      &
    !$omp private(id, i, ix, num_nb15, rtmp, qtmp, ini_nb15, fin_nb15,  & 
    !$omp         k, iy, iwater, ij, j, iix, list, num_count, L, L1,    &
    !$omp         iatmcls, jatmcls, rij2, R, term_lj12, term_lj6,       &
    !$omp         dij, lj12, lj6, grad_coef, work, term_elec, inv_r2,   &
    !$omp         trans_x, trans_y, trans_z, force_local, check_virial, &
    !$omp         fg1, fg2, rij2_sclj, rij2_scel)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell, nthread
      do ix = 1, natom(i)
        coord_pbc(ix,1,i) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(ix,2,i) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(ix,3,i) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    !$omp barrier

    ! energy within a cell(solute-solute)
    !
    do i = id+1, ncell_local, nthread

      num_nb15 = 0

      do ix = 1, natom(i) - 1

        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp      = charge(ix,i)
        iatmcls   = atmcls(ix,i)

        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc1(ix,i)
        num_nb15  = fin_nb15

        force_local(1:3) = 0.0_wp

        ! FEP: flag for atom (ix,i)
        fg1 = fepgrp(ix,i)

!dir$ simd
        do k = ini_nb15, fin_nb15

          iy = nb15_calc_list1(k,i)
          jatmcls = atmcls(iy,i)
          lj12 = nonb_lj12(jatmcls,iatmcls)
          lj6  = nonb_lj6(jatmcls,iatmcls)

          ! Compute distance
          dij(1) = rtmp(1) - coord_pbc(iy,1,i)
          dij(2) = rtmp(2) - coord_pbc(iy,2,i)
          dij(3) = rtmp(3) - coord_pbc(iy,3,i)
          rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          ! FEP: flag for atom (iy,i)
          fg2 = fepgrp(iy,i)

          ! FEP: LJ with soft core
          rij2_sclj  = density_tmp / (rij2 + table_sclj(fg1,fg2))
          L  = int(rij2_sclj)
          R  = rij2_sclj - L
          L1 = 3*L - 2
          term_lj12 = table_grad(L1  ) + R*(table_grad(L1+3)-table_grad(L1  ))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))

          ! FEP: elec with soft core
          rij2_scel  = density_tmp / (rij2 + table_scel(fg1,fg2))
          L  = int(rij2_scel)
          R  = rij2_scel - L
          L1 = 3*L
          term_elec = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))

          grad_coef = term_lj12*lj12 - term_lj6*lj6   &
                    + qtmp*charge(iy,i)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force(iy,1,i,id+1) = force(iy,1,i,id+1) + work(1)
          force(iy,2,i,id+1) = force(iy,2,i,id+1) + work(2)
          force(iy,3,i,id+1) = force(iy,3,i,id+1) + work(3)

        end do

        force(ix,1:3,i,id+1) = force(ix,1:3,i,id+1) + force_local(1:3)
      end do

    end do

    ! interaction between different cells
    !
    do ij = id+1, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      num_nb15 = 0

      ! FEP: shift for generic & gpu kernels
      if (domain%nonbond_kernel /= NBK_Fugaku .and. &
          domain%nonbond_kernel /= NBK_Intel) then
        trans_x = real(cell_move(1,j,i),wp) * system_size(1)
        trans_y = real(cell_move(2,j,i),wp) * system_size(2)
        trans_z = real(cell_move(3,j,i),wp) * system_size(3)
      else
        trans_x = 0.0_wp
        trans_y = 0.0_wp
        trans_z = 0.0_wp
      end if

      check_virial = virial_check(j,i)

      do iix = 1, nb15_cell(ij)

        ix = nb15_list(iix,ij)
        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp = charge(ix,i)
        iatmcls = atmcls(ix,i)

        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc(ix,ij)
        num_nb15  = fin_nb15

        force_local(1:3) = 0.0_wp

        ! FEP: flag for atom (ix,i)
        fg1 = fepgrp(ix,i)

!dir$ simd
        do k = ini_nb15, fin_nb15

          iy = nb15_calc_list(k,ij)
          jatmcls = atmcls(iy,j)
          lj12 = nonb_lj12(jatmcls,iatmcls)
          lj6  = nonb_lj6(jatmcls,iatmcls)

          ! Compute distance
          dij(1) = rtmp(1) - coord_pbc(iy,1,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(iy,2,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(iy,3,j) + trans_z
          rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          ! FEP: flag for atom (iy,j)
          fg2 = fepgrp(iy,j)

          ! FEP: LJ with soft core
          rij2_sclj = density_tmp / (rij2 + table_sclj(fg1,fg2))
          L  = int(rij2_sclj)
          R  = rij2_sclj - L
          L1 = 3*L - 2
          term_lj12 = table_grad(L1  ) + R*(table_grad(L1+3)-table_grad(L1  ))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))

          ! FEP: elec with soft core
          rij2_scel = density_tmp / (rij2 + table_scel(fg1,fg2))
          L  = int(rij2_scel)
          R  = rij2_scel - L
          L1 = 3*L
          term_elec = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))

          grad_coef = term_lj12*lj12 - term_lj6*lj6   &
                    + qtmp*charge(iy,j)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force(iy,1,j,id+1) = force(iy,1,j,id+1) + work(1)
          force(iy,2,j,id+1) = force(iy,2,j,id+1) + work(2)
          force(iy,3,j,id+1) = force(iy,3,j,id+1) + work(3)

        end do

        if (check_virial == 1) &
          virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
        force(ix,1:3,i,id+1) = force(ix,1:3,i,id+1) + force_local(1:3)
      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_force_nonbond_tbl_lnr_generic_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_notbl_generic_fep
  !> @brief        calculate nonbonded energy without lookup table of vdw
  !                for FEP
  !! @authors      HO
  !! @param[in]    domain    : domain information
  !! @param[in]    enefunc   : potential energy functions
  !! @param[in]    pairlist  : interaction list in each domain
  !! @param[inout] coord_pbc : pbc oriented coordinates
  !! @param[inout] force     : forces for each cell
  !! @param[inout] virial    : virial term of target systems
  !! @param[inout] eelec     : electrostatic energy of target systems
  !! @param[inout] evdw      : van der Waals energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_notbl_generic_fep( &
                                                 domain, enefunc, pairlist, &
                                                 virial,  &
                                                 eelec, evdw)

    ! formal arguments
    type(s_domain),   target, intent(inout) :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(dp),                 intent(inout) :: virial(:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                  :: dij(1:3), rij2, rij2_inv
    real(wp)                  :: lj6, lj12
    real(wp)                  :: R
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, grad_coef, within_cutoff
    real(wp)                  :: work(1:3)
    real(wp)                  :: trans_x, trans_y, trans_z
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: elec_temp, evdw_temp
    real(wp)                  :: force_local(3)
    real(wp)                  :: density_tmp, qxy
    integer                   :: i, ix, iy, j, k, ij, iix, L, L1
    integer                   :: num_nb15, ini_nb15, fin_nb15
    integer                   :: id, omp_get_thread_num
    integer                   :: ncell, ncell_local
    integer                   :: check_virial
    integer                   :: iatmcls, jatmcls

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: atmcls(:,:)
    integer(int2),    pointer :: cell_pairlist(:,:)
    integer,          pointer :: natom(:)
    integer,          pointer :: nb15_list(:,:)
    integer,          pointer :: nb15_cell(:)
    integer,          pointer :: num_nb15_calc1(:,:), num_nb15_calc(:,:)
    integer,          pointer :: nb15_calc_list1(:,:), nb15_calc_list(:,:)
    integer(int1),    pointer :: virial_check(:,:)

    ! FEP
    real(wp)                  :: rij2_sclj, rij2_scel
    integer                   :: fg1, fg2
    integer,          pointer :: fepgrp(:,:)
    real(wp),         pointer :: table_sclj(:,:)
    real(wp),         pointer :: table_scel(:,:)
    real(wp),         pointer :: coord_pbc(:,:,:)
    real(wp),         pointer :: force(:,:,:,:)

    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff

    ! FEP
    nb15_list       => pairlist%nb15_list_fep
    nb15_cell       => pairlist%nb15_cell_fep
    num_nb15_calc1  => pairlist%num_nb15_calc1_fep
    num_nb15_calc   => pairlist%num_nb15_calc_fep
    nb15_calc_list1 => pairlist%nb15_calc_list1_fep
    nb15_calc_list  => pairlist%nb15_calc_list_fep
    fepgrp          => domain%fepgrp
    table_sclj      => enefunc%table_sclj
    table_scel      => enefunc%table_scel
    coord_pbc       => domain%translated_fep
    force           => domain%force_pbc_fep

    density_tmp     = cutoff2 * density

    ! calculate energy and gradient
    !

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, ix, num_nb15, rtmp, qtmp, ini_nb15, fin_nb15, k, iy,  &
    !$omp         rij2, L, R, term_lj12, term_lj6, term_elec, grad_coef, work, &
    !$omp         ij, j, trans_x, trans_y, trans_z, iix, force_local, lj6,     &
    !$omp         lj12, L1, elec_temp, evdw_temp, dij, check_virial, rij2_inv, &
    !$omp         within_cutoff, iatmcls, jatmcls, qxy,                        &
    !$omp         fg1, fg2, rij2_sclj, rij2_scel)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread
      do ix = 1, natom(i)
        coord_pbc(ix,1,i) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(ix,2,i) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(ix,3,i) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    !$omp barrier

    ! energy within a cell
    !
    do i = id+1, ncell_local, nthread

      num_nb15 = 0
      do ix = 1, natom(i) - 1

        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp = charge(ix,i)
        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc1(ix,i)
        num_nb15  = fin_nb15
        force_local(1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp

        iatmcls = atmcls(ix,i)

        ! FEP: flag for atom (ix,i)
        fg1 = fepgrp(ix,i)

!dir$ simd
        do k = ini_nb15, fin_nb15

          iy = nb15_calc_list1(k,i)
          jatmcls = atmcls(iy,i)
          lj6  = nonb_lj6 (iatmcls,jatmcls)
          lj12 = nonb_lj12(iatmcls,jatmcls)
          qxy = qtmp*charge(iy,i)

          ! compute distance
          !
          dij(1) = rtmp(1) - coord_pbc(iy,1,i)
          dij(2) = rtmp(2) - coord_pbc(iy,2,i)
          dij(3) = rtmp(3) - coord_pbc(iy,3,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)

          ! FEP: flag for atom (iy,i)
          fg2 = fepgrp(iy,i)

          ! FEP: LJ with soft core
          rij2_inv = 1.0_wp / (rij2 + table_sclj(fg1,fg2))
          term_lj6  = rij2_inv * rij2_inv * rij2_inv * within_cutoff
          term_lj12 = term_lj6 * term_lj6
          evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
          term_lj12 = -12.0_wp * term_lj12 * rij2_inv
          term_lj6  = - 6.0_wp * term_lj6  * rij2_inv

          ! FEP: elec with soft core
          rij2_scel  = density_tmp / (rij2 + table_scel(fg1,fg2))
          L    = int(rij2_scel)
          R    = rij2_scel - L
          term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
          term_elec = term_elec * within_cutoff
          elec_temp = elec_temp + qxy*term_elec
          term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
          term_elec = term_elec * within_cutoff

          grad_coef = term_lj12*lj12 - term_lj6*lj6 + qxy*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force(iy,1,i,id+1) = force(iy,1,i,id+1) + work(1)
          force(iy,2,i,id+1) = force(iy,2,i,id+1) + work(2)
          force(iy,3,i,id+1) = force(iy,3,i,id+1) + work(3)

        end do

        force(ix,1:3,i,id+1) = force(ix,1:3,i,id+1) + force_local(1:3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp

      end do

    end do


    ! interaction between different cells
    !
    do ij = id+1, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      num_nb15 = 0

      ! FEP: shift for generic & gpu kernels
      if (domain%nonbond_kernel /= NBK_Fugaku .and. &
          domain%nonbond_kernel /= NBK_Intel) then
        trans_x = real(cell_move(1,j,i),wp) * system_size(1)
        trans_y = real(cell_move(2,j,i),wp) * system_size(2)
        trans_z = real(cell_move(3,j,i),wp) * system_size(3)
      else
        trans_x = 0.0_wp
        trans_y = 0.0_wp
        trans_z = 0.0_wp
      end if

      check_virial = virial_check(j,i)

      do iix = 1, nb15_cell(ij)

        ix   = nb15_list(iix,ij)
        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp = charge(ix,i)

        ini_nb15 = num_nb15 + 1
        fin_nb15 = num_nb15 + num_nb15_calc(ix,ij)
        num_nb15 = fin_nb15

        force_local(1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp

        iatmcls = atmcls(ix,i)

        ! FEP: flag for atom (ix,i)
        fg1 = fepgrp(ix,i)

!dir$ simd
        do k = ini_nb15, fin_nb15
          iy   = nb15_calc_list(k,ij)
          jatmcls = atmcls(iy,j)
          lj6  = nonb_lj6 (iatmcls,jatmcls)
          lj12 = nonb_lj12(iatmcls,jatmcls)
          qxy = qtmp*charge(iy,j)

          ! Compute distance
          dij(1) = rtmp(1) - coord_pbc(iy,1,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(iy,2,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(iy,3,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)

          ! FEP: flag for atom (iy,j)
          fg2 = fepgrp(iy,j)

          ! FEP: LJ with soft core
          rij2_inv = 1.0_wp / (rij2 + table_sclj(fg1,fg2))
          term_lj6  = rij2_inv * rij2_inv * rij2_inv * within_cutoff
          term_lj12 = term_lj6 * term_lj6
          evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
          term_lj12 = -12.0_wp * term_lj12 * rij2_inv
          term_lj6  = - 6.0_wp * term_lj6  * rij2_inv

          ! FEP: elec with soft core
          rij2_scel  = density_tmp / (rij2 + table_scel(fg1,fg2))
          L    = int(rij2_scel)
          R    = rij2_scel - L
          term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
          term_elec = term_elec * within_cutoff
          elec_temp = elec_temp + qxy*term_elec
          term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
          term_elec = term_elec * within_cutoff

          grad_coef = term_lj12*lj12 - term_lj6*lj6 + qxy*term_elec
  
          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)
  
          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
  
          force(iy,1,j,id+1) = force(iy,1,j,id+1) + work(1)
          force(iy,2,j,id+1) = force(iy,2,j,id+1) + work(2)
          force(iy,3,j,id+1) = force(iy,3,j,id+1) + work(3)
  
        end do

        if (check_virial == 1) &
          virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
        force(ix,1:3,i,id+1) = force(ix,1:3,i,id+1) + force_local(1:3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond_notbl_generic_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
 !
  !  Subroutine    compute_energy_nonbond_notbl_chk_generic_fep
  !> @brief        calculate nonbonded energy without lookup table of vdw
  !                for FEP
  !! @authors      HO
  !! @param[in]    domain    : domain information
  !! @param[in]    enefunc   : potential energy functions
  !! @param[in]    pairlist  : interaction list in each domain
  !! @param[inout] coord_pbc : pbc oriented coordinates
  !! @param[inout] force     : forces for each cell
  !! @param[inout] virial    : virial term of target systems
  !! @param[inout] eelec     : electrostatic energy of target systems
  !! @param[inout] evdw      : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_notbl_chk_generic_fep( &
                                                 domain, enefunc, pairlist, &
                                                 virial,  &
                                                 eelec, evdw)

    ! formal arguments
    type(s_domain),   target, intent(inout) :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(dp),                 intent(inout) :: virial(:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                  :: dij(1:3), rij2, rij2_inv
    real(wp)                  :: lj6, lj12
    real(wp)                  :: R
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, term_lj, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: trans_x, trans_y, trans_z
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: elec_temp, evdw_temp
    real(wp)                  :: dij_list(3,MaxAtom)
    real(wp)                  :: force_local(3), rij2_list(MaxAtom)
    real(wp)                  :: force_localj(3,MaxAtom)
    real(wp)                  :: minimum_contact
    real(wp)                  :: density_tmp, qxy
    integer                   :: i, ix, iy, j, k, ij, iix, iwater, list, L, L1
    integer                   :: num_nb15, ini_nb15, fin_nb15, num_count
    integer                   :: id, omp_get_thread_num, j_list(MaxAtom)
    integer                   :: ncell, ncell_local
    integer                   :: check_virial
    integer                   :: iatmcls, jatmcls

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: atmcls(:,:)
    integer(int2),    pointer :: cell_pairlist(:,:)
    integer,          pointer :: natom(:)
    integer,          pointer :: nb15_list(:,:)
    integer,          pointer :: nb15_cell(:)
    integer,          pointer :: num_nb15_calc1(:,:), num_nb15_calc(:,:)
    integer,          pointer :: nb15_calc_list1(:,:), nb15_calc_list(:,:)
    integer(int1),    pointer :: virial_check(:,:)

    ! FEP
    real(wp)                  :: rij2_sclj, rij2_scel
    integer                   :: fg1, fg2
    integer,          pointer :: fepgrp(:,:)
    real(wp),         pointer :: table_sclj(:,:)
    real(wp),         pointer :: table_scel(:,:)
    real(wp),         pointer :: coord_pbc(:,:,:)
    real(wp),         pointer :: force(:,:,:,:)

    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff
    minimum_contact =  enefunc%minimum_contact

    ! FEP
    nb15_list       => pairlist%nb15_list_fep
    nb15_cell       => pairlist%nb15_cell_fep
    num_nb15_calc1  => pairlist%num_nb15_calc1_fep
    num_nb15_calc   => pairlist%num_nb15_calc_fep
    nb15_calc_list1 => pairlist%nb15_calc_list1_fep
    nb15_calc_list  => pairlist%nb15_calc_list_fep
    fepgrp          => domain%fepgrp
    table_sclj      => enefunc%table_sclj
    table_scel      => enefunc%table_scel
    coord_pbc       => domain%translated_fep
    force           => domain%force_pbc_fep

    density_tmp     = cutoff2 * density

    ! calculate energy and gradient
    !

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, ix, num_nb15, rtmp, qtmp, ini_nb15, fin_nb15, k, iy,  &
    !$omp         rij2, L, R, term_lj12, term_lj6, term_lj, grad_coef, work,   &
    !$omp         term_elec, iwater, ij, j, trans_x, trans_y,                  &
    !$omp         trans_z, iix, list, num_count, j_list, dij_list, rij2_list,  &
    !$omp         force_local, force_localj, lj6, lj12, L1, elec_temp,         &
    !$omp         evdw_temp, dij, check_virial, rij2_inv, iatmcls, jatmcls,    &
    !$omp         fg1, fg2, rij2_sclj, rij2_scel, qxy)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread
      do ix = 1, natom(i)
        coord_pbc(ix,1,i) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(ix,2,i) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(ix,3,i) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    !$omp barrier

    ! energy within a cell
    !
    do i = id+1, ncell_local, nthread

      num_nb15 = 0
      do ix = 1, natom(i) - 1

        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp = charge(ix,i)
        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc1(ix,i)
        num_nb15  = fin_nb15
        num_count = 0

        do k = ini_nb15, fin_nb15

          iy = nb15_calc_list1(k,i)

          ! compute distance
          !
          dij(1) = rtmp(1) - coord_pbc(iy,1,i)
          dij(2) = rtmp(2) - coord_pbc(iy,2,i)
          dij(3) = rtmp(3) - coord_pbc(iy,3,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1:3,num_count)  = dij(1:3)
            rij2_list(num_count) = max(rij2, minimum_contact)
            j_list(num_count)   = iy
          end if
        end do

        force_local(1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp

        iatmcls = atmcls(ix,i)

        ! FEP: flag for atom (ix,i)
        fg1 = fepgrp(ix,i)

        do k = 1, num_count
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          iy   = j_list(k)
          jatmcls = atmcls(iy,i)
          lj6  = nonb_lj6 (iatmcls,jatmcls)
          lj12 = nonb_lj12(iatmcls,jatmcls)
          qxy = qtmp*charge(iy,i)

          ! FEP: flag for atom (iy,i)
          fg2 = fepgrp(iy,i)

          ! FEP: LJ with soft core
          rij2_inv = 1.0_wp / (rij2_list(k) + table_sclj(fg1,fg2))
          term_lj6  = rij2_inv * rij2_inv * rij2_inv
          term_lj12 = term_lj6 * term_lj6
          evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
          term_lj12 = -12.0_wp*term_lj12*rij2_inv
          term_lj6  = - 6.0_wp*term_lj6*rij2_inv

          ! FEP: elec with soft core
          rij2_scel = density_tmp / (rij2_list(k) + table_scel(fg1,fg2))
          L    = int(rij2_scel)
          R    = rij2_scel - L
          term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
          elec_temp = elec_temp + qxy*term_elec
          term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))

          grad_coef = term_lj12*lj12 - term_lj6*lj6 + qxy*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force_localj(1,k)  = work(1)
          force_localj(2,k)  = work(2)
          force_localj(3,k)  = work(3)

        end do

        do k = 1, num_count
          iy = j_list(k)
          force(iy,1,i,id+1) = force(iy,1,i,id+1) + force_localj(1,k)
          force(iy,2,i,id+1) = force(iy,2,i,id+1) + force_localj(2,k)
          force(iy,3,i,id+1) = force(iy,3,i,id+1) + force_localj(3,k)
        end do
        force(ix,1:3,i,id+1) = force(ix,1:3,i,id+1) + force_local(1:3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp

      end do

    end do


    ! interaction between different cells
    !
    do ij = id+1, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      num_nb15 = 0

      ! FEP: shift for generic & gpu kernels
      if (domain%nonbond_kernel /= NBK_Fugaku .and. &
          domain%nonbond_kernel /= NBK_Intel) then
        trans_x = real(cell_move(1,j,i),wp) * system_size(1)
        trans_y = real(cell_move(2,j,i),wp) * system_size(2)
        trans_z = real(cell_move(3,j,i),wp) * system_size(3)
      else
        trans_x = 0.0_wp
        trans_y = 0.0_wp
        trans_z = 0.0_wp
      end if

      check_virial = virial_check(j,i)

      do iix = 1, nb15_cell(ij)

        ix   = nb15_list(iix,ij)
        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp = charge(ix,i)

        ini_nb15 = num_nb15 + 1
        fin_nb15 = num_nb15 + num_nb15_calc(ix,ij)
        num_nb15 = fin_nb15

        num_count = 0
        do k = ini_nb15, fin_nb15
          iy   = nb15_calc_list(k,ij)
          dij(1) = rtmp(1) - coord_pbc(iy,1,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(iy,2,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(iy,3,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1:3,num_count)  = dij(1:3)
            rij2_list(num_count) = max(rij2, minimum_contact)
            j_list(num_count)   = iy
          end if
        end do

        force_local(1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp

        iatmcls = atmcls(ix,i)

        ! FEP: flag for atom (ix,i)
        fg1 = fepgrp(ix,i)

!ocl norecurrence
!ocl swp
        do k = 1, num_count
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          iy   = j_list(k)
          jatmcls = atmcls(iy,j)
          lj6  = nonb_lj6 (iatmcls,jatmcls)
          lj12 = nonb_lj12(iatmcls,jatmcls)
          qxy = qtmp*charge(iy,j)

          ! FEP: flag for atom (iy,j)
          fg2 = fepgrp(iy,j)

          ! FEP: LJ with soft core
          rij2_inv = 1.0_wp / (rij2_list(k) + table_sclj(fg1,fg2))
          term_lj6  = rij2_inv * rij2_inv * rij2_inv
          term_lj12 = term_lj6 * term_lj6
          evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
          term_lj12 = -12.0_wp*term_lj12*rij2_inv
          term_lj6  =  -6.0_wp*term_lj6 *rij2_inv

          ! FEP: elec with soft core
          rij2_scel  = density_tmp / (rij2_list(k) + table_scel(fg1,fg2))
          L    = int(rij2_scel)
          R    = rij2_scel - L
          term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
          elec_temp = elec_temp + qxy*term_elec
          term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))

          grad_coef = term_lj12*lj12 - term_lj6*lj6 + qxy*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)

          force_localj(1,k) = work(1)
          force_localj(2,k) = work(2)
          force_localj(3,k) = work(3)

        end do

        do k = 1, num_count
          iy = j_list(k)
          force(iy,1,j,id+1) = force(iy,1,j,id+1) + force_localj(1,k)
          force(iy,2,j,id+1) = force(iy,2,j,id+1) + force_localj(2,k)
          force(iy,3,j,id+1) = force(iy,3,j,id+1) + force_localj(3,k)
        end do
        if (check_virial == 1) &
          virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
        force(ix,1:3,i,id+1) = force(ix,1:3,i,id+1) + force_local(1:3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond_notbl_chk_generic_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_notbl_generic_fep
  !> @brief        calculate nonbonded force without lookup table of vdw
  !                for FEP
  !! @authors      HO
  !! @param[in]    domain    : domain information
  !! @param[in]    enefunc   : potential energy functions
  !! @param[in]    pairlist  : interaction list in each domain
  !! @param[inout] coord_pbc : pbc oriented coordinates
  !! @param[inout] force     : forces for each cell
  !! @param[inout] virial    : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_notbl_generic_fep( &
                                                 domain, enefunc, pairlist, &
                                                 virial)

    ! formal arguments
    type(s_domain),   target, intent(inout) :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(dp),                 intent(inout) :: virial(:,:)

    ! local variables
    real(wp)                  :: dij(1:3), rij2, inv_rij2
    real(wp)                  :: R, lj12, lj6, factor
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: trans_x, trans_y, trans_z
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: force_local(3)
    real(wp)                  :: density_tmp
    integer                   :: i, ix, iy, j, k, ij, iix, iwater, list, L, L1
    integer                   :: num_nb15, ini_nb15, fin_nb15, num_count
    integer                   :: id, omp_get_thread_num
    integer                   :: iatmcls,jatmcls
    integer                   :: ncell, ncell_local
    integer                   :: check_virial

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_grad(:)
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: atmcls(:,:)
    integer(int2),    pointer :: cell_pairlist(:,:)
    integer,          pointer :: natom(:)
    integer,          pointer :: nb15_list(:,:)
    integer,          pointer :: nb15_cell(:)
    integer,          pointer :: num_nb15_calc1(:,:), num_nb15_calc(:,:)
    integer,          pointer :: nb15_calc_list1(:,:), nb15_calc_list(:,:)
    integer(int1),    pointer :: virial_check(:,:)

    ! FEP
    real(wp)                  :: rij2_sclj, rij2_scel
    integer                   :: fg1, fg2
    integer,          pointer :: fepgrp(:,:)
    real(wp),         pointer :: table_sclj(:,:)
    real(wp),         pointer :: table_scel(:,:)
    real(wp),         pointer :: coord_pbc(:,:,:)
    real(wp),         pointer :: force(:,:,:,:)

    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff

    ! FEP
    nb15_list       => pairlist%nb15_list_fep
    nb15_cell       => pairlist%nb15_cell_fep
    num_nb15_calc1  => pairlist%num_nb15_calc1_fep
    num_nb15_calc   => pairlist%num_nb15_calc_fep
    nb15_calc_list1 => pairlist%nb15_calc_list1_fep
    nb15_calc_list  => pairlist%nb15_calc_list_fep
    fepgrp          => domain%fepgrp
    table_sclj      => enefunc%table_sclj
    table_scel      => enefunc%table_scel
    coord_pbc       => domain%translated_fep
    force           => domain%force_pbc_fep

    density_tmp     = cutoff2 * density

    ! calculate energy and gradient
    !
    !$omp parallel default(shared)                                      &
    !$omp private(id, i, ix, num_nb15, rtmp, qtmp, ini_nb15, fin_nb15,  & 
    !$omp         k, iy, iwater, ij, j, iix, list, num_count, L,        &
    !$omp         iatmcls, jatmcls, rij2, R, term_lj12, term_lj6,       &
    !$omp         dij, lj12, lj6, grad_coef, work, term_elec, inv_rij2, &
    !$omp         trans_x, trans_y, trans_z, force_local, check_virial, &
    !$omp         factor,                                               &
    !$omp         fg1, fg2, rij2_sclj, rij2_scel)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell, nthread
      do ix = 1, natom(i)
        coord_pbc(ix,1,i) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(ix,2,i) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(ix,3,i) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    !$omp barrier

    ! energy within a cell(solute-solute)
    !
    do i = id+1, ncell_local, nthread

      num_nb15 = 0

      do ix = 1, natom(i) - 1

        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp      = charge(ix,i)
        iatmcls   = atmcls(ix,i)

        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc1(ix,i)
        num_nb15  = fin_nb15

        force_local(1:3) = 0.0_wp

        ! FEP: flag for atom (ix,i)
        fg1 = fepgrp(ix,i)

!dir$ simd
        do k = ini_nb15, fin_nb15

          iy = nb15_calc_list1(k,i)
          jatmcls = atmcls(iy,i)
          lj12 = nonb_lj12(jatmcls,iatmcls)
          lj6  = nonb_lj6(jatmcls,iatmcls)

          ! Compute distance
          dij(1) = rtmp(1) - coord_pbc(iy,1,i)
          dij(2) = rtmp(2) - coord_pbc(iy,2,i)
          dij(3) = rtmp(3) - coord_pbc(iy,3,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          factor = 0.0_wp
          factor = merge(0.0_wp,1.0_wp,rij2>cutoff2)

          ! FEP: flag for atom (iy,i)
          fg2 = fepgrp(iy,i)

          ! FEP: LJ with soft core
          inv_rij2 = 1.0_wp / (rij2 + table_sclj(fg1,fg2))
          term_lj6  = inv_rij2 * inv_rij2 * inv_rij2
          term_lj12 = term_lj6 * term_lj6
          term_lj12 = -12.0_wp * term_lj12 * inv_rij2
          term_lj6  = - 6.0_wp * term_lj6  * inv_rij2

          ! FEP: elec with soft core
          rij2_scel  = density_tmp / (rij2 + table_scel(fg1,fg2))
          L  = int(rij2_scel)
          R  = rij2_scel - L
          term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))

          grad_coef = term_lj12*lj12 - term_lj6*lj6   &
                    + qtmp*charge(iy,i)*term_elec
          grad_coef = grad_coef * factor

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force(iy,1,i,id+1) = force(iy,1,i,id+1) + work(1)
          force(iy,2,i,id+1) = force(iy,2,i,id+1) + work(2)
          force(iy,3,i,id+1) = force(iy,3,i,id+1) + work(3)

        end do

        force(ix,1:3,i,id+1) = force(ix,1:3,i,id+1) + force_local(1:3)

      end do

    end do

    ! interaction between different cells
    !
    do ij = id+1, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      num_nb15 = 0

      ! FEP: shift for generic & gpu kernels
      if (domain%nonbond_kernel /= NBK_Fugaku .and. &
          domain%nonbond_kernel /= NBK_Intel) then
        trans_x = real(cell_move(1,j,i),wp) * system_size(1)
        trans_y = real(cell_move(2,j,i),wp) * system_size(2)
        trans_z = real(cell_move(3,j,i),wp) * system_size(3)
      else
        trans_x = 0.0_wp
        trans_y = 0.0_wp
        trans_z = 0.0_wp
      end if

      check_virial = virial_check(j,i)

      do iix = 1, nb15_cell(ij)

        ix = nb15_list(iix,ij)
        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp = charge(ix,i)
        iatmcls = atmcls(ix,i)

        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc(ix,ij)
        num_nb15  = fin_nb15

        force_local(1:3) = 0.0_wp

        ! FEP: flag for atom (ix,i)
        fg1 = fepgrp(ix,i)

!dir$ simd
        do k = ini_nb15, fin_nb15

          iy = nb15_calc_list(k,ij)
          jatmcls = atmcls(iy,j)
          lj12 = nonb_lj12(jatmcls,iatmcls)
          lj6  = nonb_lj6(jatmcls,iatmcls)

          ! Compute distance
          dij(1) = rtmp(1) - coord_pbc(iy,1,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(iy,2,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(iy,3,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          factor = 0.0_wp
          factor = merge(0.0_wp,1.0_wp,rij2>cutoff2)

          ! FEP: flag for atom (iy,j)
          fg2 = fepgrp(iy,j)

          ! FEP: LJ with soft core
          inv_rij2 = 1.0_wp / (rij2 + table_sclj(fg1,fg2))
          term_lj6  = inv_rij2 * inv_rij2 * inv_rij2
          term_lj12 = term_lj6 * term_lj6
          term_lj12 = -12.0_wp * term_lj12 * inv_rij2
          term_lj6  = - 6.0_wp * term_lj6 * inv_rij2

          ! FEP: elec with soft core
          rij2_scel = density_tmp / (rij2 + table_scel(fg1,fg2))
          L  = int(rij2_scel)
          R  = rij2_scel - L
          term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))

          grad_coef = term_lj12*lj12 - term_lj6*lj6   &
                    + qtmp*charge(iy,j)*term_elec
          grad_coef = grad_coef * factor

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force(iy,1,j,id+1) = force(iy,1,j,id+1) + work(1)
          force(iy,2,j,id+1) = force(iy,2,j,id+1) + work(2)
          force(iy,3,j,id+1) = force(iy,3,j,id+1) + work(3)

        end do

        if (check_virial == 1) &
          virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
        force(ix,1:3,i,id+1) = force(ix,1:3,i,id+1) + force_local(1:3)
      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_force_nonbond_notbl_generic_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_tbl_lnr_gpu_fep
  !> @brief        calculate nonbonded energy with gpu for FEP
  !! @authors      HO
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[in]    npt        : logical to check virial calculation is required
  !! @param[inout] coord_pbc  : pbc oriented coordinates
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] eelec      : electrostatic energy of target systems
  !! @param[inout] evdw       : van der Waals energy of target systems
  !! @param[inout] ene_virial : energy&virial outputs
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_tbl_lnr_gpu_fep( &
                                                 domain, enefunc, pairlist, &
                                                 npt, coord_pbc, force,     &
                                                 virial, eelec, evdw,       &
                                                 ene_virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    logical,                  intent(in)    :: npt
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)
    real(dp),                 intent(inout) :: ene_virial(:)

#ifdef USE_GPU

    ! local variables
    real(wp)                  :: cutoff2, pairlistdist2
    integer                   :: ncell, ncell_local, ncell_bound

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: nonb_lj6_factor(:)
    real(wp),         pointer :: density, cutoff, pairlistdist
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: atmcls(:,:)
    integer,          pointer :: natom(:), start_atom(:)
    integer(int1),    pointer :: virial_check(:,:)

    ! for GPU
    integer,          pointer :: univ_cell_pairlist1(:,:)
    integer,          pointer :: univ_cell_pairlist2(:,:)
    integer,          pointer :: univ_ix_natom(:)
    integer(1),       pointer :: univ_ix_list(:,:)
    integer,          pointer :: univ_iy_natom(:)
    integer(1),       pointer :: univ_iy_list(:,:)
    integer,          pointer :: univ_ij_sort_list(:)
    integer(1),       pointer :: univ_mask2(:,:)
    integer,          pointer :: univ_mask2_size
    integer                   :: univ_ncell_nonzero
    integer                   :: univ_update
      
    integer                   :: ncell_max, ij, i, j, ix, ixx, start_i
    integer                   :: cutoff_int
    integer                   :: num_atom_cls
    integer                   :: ret_shape(1)
    integer                   :: univ_ij, univ_ij0, iiy, idx, base_idx

    ! for FEP
    integer,          pointer :: fepgrp(:,:)
    integer(1),       pointer :: fep_mask(:,:)
    real(wp),         pointer :: table_sclj(:,:)
    real(wp),         pointer :: table_scel(:,:)
    integer(1),       pointer :: fepgrp_pbc(:)

    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6
    nonb_lj6_factor => enefunc%nonb_lj6_factor
    pairlistdist    => enefunc%pairlistdist

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    ncell_bound     =  domain%num_cell_boundary
    cutoff2         =  cutoff * cutoff
    pairlistdist2   =  pairlistdist * pairlistdist

    ret_shape = shape(natom)
    ncell_max = ret_shape(1)
    cutoff_int = enefunc%table%cutoff_int
    num_atom_cls = enefunc%num_atom_cls

    !
    ! univ
    univ_cell_pairlist1 => domain%univ_cell_pairlist1
    univ_cell_pairlist2 => domain%univ_cell_pairlist2
    univ_ix_natom       => pairlist%univ_ix_natom
    univ_ix_list        => pairlist%univ_ix_list
    univ_iy_natom       => pairlist%univ_iy_natom
    univ_iy_list        => pairlist%univ_iy_list
    univ_ij_sort_list   => pairlist%univ_ij_sort_list
    univ_mask2          => pairlist%univ_mask2
    univ_mask2_size     => pairlist%univ_mask2_size
    univ_ncell_nonzero  =  pairlist%univ_ncell_nonzero
    univ_update         =  pairlist%univ_update

    ! FEP
    fepgrp              => domain%fepgrp
    fep_mask            => enefunc%fep_mask
    table_sclj          => enefunc%table_sclj
    table_scel          => enefunc%table_scel
    fepgrp_pbc          => domain%fepgrp_pbc

    ij = domain%num_atom_domain
    !$omp parallel do private(i,ix,ixx,start_i)
    do i = 1, ncell
      start_i = start_atom(i)
      do ix = 1, natom(i)
        coord_pbc(     start_i+ix,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(  ij+start_i+ix,1,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(2*ij+start_i+ix,1,1) = coord(3,ix,i) + trans1(3,ix,i)
        fepgrp_pbc(start_i+ix) = fepgrp(ix,i)
      end do
    end do
    !$omp end parallel do

    !
    ! launch GPU kernels
    call gpu_launch_compute_energy_nonbond_table_linear_univ_fep( &
         coord_pbc, force(:,:,:,1), ene_virial,               &
         cell_move, natom, start_atom,                        &
         nonb_lj12, nonb_lj6, nonb_lj6_factor, table_ene,     &
         table_grad, univ_cell_pairlist1, univ_mask2,         &
         univ_ix_natom, univ_ix_list, univ_iy_natom,          &
         univ_iy_list, univ_ij_sort_list, virial_check,       &
         fepgrp_pbc, fep_mask, table_sclj, table_scel,        &
         domain%num_atom_domain, MaxAtom, MaxAtomCls,         &
         num_atom_cls, ncell_local, ncell_bound, ncell_max,   &
         cutoff_int, univ_maxcell, univ_maxcell1,             &
         univ_ncell_nonzero, univ_ncell_near, univ_update,    &
         univ_mask2_size, univ_natom_max, maxcell, density,   &
         cutoff2, system_size(1), system_size(2), system_size(3) )

#endif

    return

  end subroutine compute_energy_nonbond_tbl_lnr_gpu_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_notbl_gpu_fep
  !> @brief        calculate nonbonded energy with gpu for FEP
  !! @authors      HO
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[in]    npt        : logical to check virial calculation is required
  !! @param[inout] coord_pbc  : pbc oriented coordinates
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] eelec      : electrostatic energy of target systems
  !! @param[inout] evdw       : van der Waals energy of target systems
  !! @param[inout] ene_virial : energy&virial outputs
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_notbl_gpu_fep( &
                                                 domain, enefunc, pairlist, &
                                                 npt, coord_pbc, force,     &
                                                 virial, eelec, evdw,       &
                                                 ene_virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    logical,                  intent(in)    :: npt
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)
    real(dp),                 intent(inout) :: ene_virial(:)

#ifdef USE_GPU

    ! local variables
    real(wp)                  :: cutoff2, pairlistdist2
    integer                   :: ncell, ncell_local, ncell_bound

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: nonb_lj6_factor(:)
    real(wp),         pointer :: density, cutoff, pairlistdist
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: atmcls(:,:)
    integer,          pointer :: natom(:), start_atom(:)
    integer(int1),    pointer :: virial_check(:,:)

    ! for GPU
    integer,          pointer :: univ_cell_pairlist1(:,:)
    integer,          pointer :: univ_cell_pairlist2(:,:)
    integer,          pointer :: univ_ix_natom(:)
    integer(1),       pointer :: univ_ix_list(:,:)
    integer,          pointer :: univ_iy_natom(:)
    integer(1),       pointer :: univ_iy_list(:,:)
    integer,          pointer :: univ_ij_sort_list(:)
    integer(1),       pointer :: univ_mask2(:,:)
    integer,          pointer :: univ_mask2_size
    integer                   :: univ_ncell_nonzero
    integer                   :: univ_update
      
    integer                   :: ncell_max, ij, i, j, ix, ixx, start_i
    integer                   :: cutoff_int
    integer                   :: num_atom_cls
    integer                   :: ret_shape(1)
    integer                   :: univ_ij, univ_ij0, iiy, idx, base_idx

    ! for FEP
    integer,          pointer :: fepgrp(:,:)
    integer(1),       pointer :: fep_mask(:,:)
    real(wp),         pointer :: table_sclj(:,:)
    real(wp),         pointer :: table_scel(:,:)
    integer(1),       pointer :: fepgrp_pbc(:)

    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6
    nonb_lj6_factor => enefunc%nonb_lj6_factor
    pairlistdist    => enefunc%pairlistdist

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    ncell_bound     =  domain%num_cell_boundary
    cutoff2         =  cutoff * cutoff
    pairlistdist2   =  pairlistdist * pairlistdist

    ret_shape = shape(natom)
    ncell_max = ret_shape(1)
    cutoff_int = enefunc%table%cutoff_int
    num_atom_cls = enefunc%num_atom_cls

    !
    ! univ
    univ_cell_pairlist1 => domain%univ_cell_pairlist1
    univ_cell_pairlist2 => domain%univ_cell_pairlist2
    univ_ix_natom       => pairlist%univ_ix_natom
    univ_ix_list        => pairlist%univ_ix_list
    univ_iy_natom       => pairlist%univ_iy_natom
    univ_iy_list        => pairlist%univ_iy_list
    univ_ij_sort_list   => pairlist%univ_ij_sort_list
    univ_mask2          => pairlist%univ_mask2
    univ_mask2_size     => pairlist%univ_mask2_size
    univ_ncell_nonzero  =  pairlist%univ_ncell_nonzero
    univ_update         =  pairlist%univ_update


    ! FEP
    fepgrp              => domain%fepgrp
    fep_mask            => enefunc%fep_mask
    table_sclj          => enefunc%table_sclj
    table_scel          => enefunc%table_scel
    fepgrp_pbc          => domain%fepgrp_pbc

    ij = domain%num_atom_domain
    !$omp parallel do private(i,ix,ixx,start_i)
    do i = 1, ncell
      start_i = start_atom(i)
      do ix = 1, natom(i)
        coord_pbc(     start_i+ix,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(  ij+start_i+ix,1,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(2*ij+start_i+ix,1,1) = coord(3,ix,i) + trans1(3,ix,i)
        fepgrp_pbc(start_i+ix) = fepgrp(ix,i)
      end do
    end do
    !$omp end parallel do

    ! launch GPU kernels
    call gpu_launch_compute_energy_nonbond_notable_univ_fep(  &
         coord_pbc, force(:,:,:,1), ene_virial,               &
         cell_move, natom, start_atom,                        &
         nonb_lj12, nonb_lj6, nonb_lj6_factor, table_ene,     &
         table_grad, univ_cell_pairlist1, univ_mask2,         &
         univ_ix_natom, univ_ix_list, univ_iy_natom,          &
         univ_iy_list, univ_ij_sort_list, virial_check,       &
         fepgrp_pbc, fep_mask, table_sclj, table_scel,        &
         domain%num_atom_domain, MaxAtom, MaxAtomCls,         &
         num_atom_cls, ncell_local, ncell_bound, ncell_max,   &
         cutoff_int, univ_maxcell, univ_maxcell1,             &
         univ_ncell_nonzero, univ_ncell_near, univ_update,    &
         univ_mask2_size, univ_natom_max, maxcell, density,   &
         cutoff2, system_size(1), system_size(2), system_size(3) )

#endif

    return

  end subroutine compute_energy_nonbond_notbl_gpu_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_tbl_lnr_gpu_fep
  !> @brief        calculate nonbonded force with gpu for FEP
  !! @authors      HO
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[in]    npt        : logical to check virial calculation is required
  !! @param[in]    cpu_calc   : flag for cpu calculation or not
  !! @param[inout] coord_pbc  : pbc oriented coordinates
  !! @param[inout] force_omp  : temporary forces of target system
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] ene_virial : energy&virial outputs
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_tbl_lnr_gpu_fep( &
                                                domain, enefunc, pairlist, &
                                                npt, cpu_calc, coord_pbc,  &
                                                force_omp, force, virial,  &
                                                ene_virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    logical,                  intent(in)    :: npt
    logical,                  intent(in)    :: cpu_calc
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force_omp(:,:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: ene_virial(:)

#ifdef USE_GPU

    ! local variables
    real(wp)                  :: cutoff2, pairlistdist2
    integer                   :: ncell, ncell_local, ncell_bound
    ! integer                   :: check_virial

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: nonb_lj6_factor(:)
    real(wp),         pointer :: density, cutoff, pairlistdist
    real(wp),         pointer :: table_grad(:)
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: atmcls(:,:)
    integer(int2),    pointer :: cell_pairlist(:,:)
    integer,          pointer :: natom(:), start_atom(:)
    integer(int1),    pointer :: virial_check(:,:)

    ! for GPU
    integer,          pointer :: univ_cell_pairlist1(:,:)
    integer,          pointer :: univ_cell_pairlist2(:,:)
    integer,          pointer :: univ_ix_natom(:)
    integer(1),       pointer :: univ_ix_list(:,:)
    integer,          pointer :: univ_iy_natom(:)
    integer(1),       pointer :: univ_iy_list(:,:)
    integer,          pointer :: univ_ij_sort_list(:)
    integer(1),       pointer :: univ_mask2(:,:)
    integer,          pointer :: univ_mask2_size
    integer                   :: univ_ncell_nonzero, univ_gpu_start_index
    integer                   :: univ_update
    integer                   :: ncell_max
    integer                   :: cutoff_int
    integer                   :: num_atom_cls
    integer                   :: ret_shape(1)
    integer                   :: univ_ij, univ_ij0, iiy, idx, base_idx
    integer                   :: ij, start_i, ix, i

    ! for FEP
    integer,          pointer :: fepgrp(:,:)
    integer(1),       pointer :: fep_mask(:,:)
    real(wp),         pointer :: table_sclj(:,:)
    real(wp),         pointer :: table_scel(:,:)
    integer(1),       pointer :: fepgrp_pbc(:)

    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6
    nonb_lj6_factor => enefunc%nonb_lj6_factor
    pairlistdist    => enefunc%pairlistdist

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    ncell_bound     =  domain%num_cell_boundary
    cutoff2         =  cutoff * cutoff
    pairlistdist2   =  pairlistdist * pairlistdist

    ret_shape = shape(natom)
    ncell_max = ret_shape(1)
    cutoff_int = enefunc%table%cutoff_int
    num_atom_cls = enefunc%num_atom_cls

    !
    ! univ
    univ_cell_pairlist1 => domain%univ_cell_pairlist1
    univ_cell_pairlist2 => domain%univ_cell_pairlist2
    univ_ix_natom       => pairlist%univ_ix_natom
    univ_ix_list        => pairlist%univ_ix_list
    univ_iy_natom       => pairlist%univ_iy_natom
    univ_iy_list        => pairlist%univ_iy_list
    univ_ij_sort_list   => pairlist%univ_ij_sort_list
    univ_mask2          => pairlist%univ_mask2
    univ_mask2_size     => pairlist%univ_mask2_size
    univ_ncell_nonzero  =  pairlist%univ_ncell_nonzero
    univ_update         =  pairlist%univ_update

    ! FEP
    fepgrp              => domain%fepgrp
    fep_mask            => enefunc%fep_mask
    table_sclj          => enefunc%table_sclj
    table_scel          => enefunc%table_scel
    fepgrp_pbc          => domain%fepgrp_pbc

    ij = domain%num_atom_domain
    !$omp parallel do private(i,ix,start_i)
    do i = 1, ncell
      start_i = start_atom(i)
      do ix = 1, natom(i)
        coord_pbc(     start_i+ix,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(  ij+start_i+ix,1,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(2*ij+start_i+ix,1,1) = coord(3,ix,i) + trans1(3,ix,i)
        fepgrp_pbc(start_i+ix) = fepgrp(ix,i)
      end do
    end do
    !$omp end parallel do

    !  launch GPU kernels (data transfer from CPU to GPU as well)
    !
    univ_gpu_start = ncell_local
    call gpu_launch_compute_force_nonbond_table_linear_univ_fep( &
         coord_pbc, force(:,:,:,1), ene_virial,              &            
         cell_move, natom, start_atom, nonb_lj12, nonb_lj6,  &
         nonb_lj6_factor, table_grad, univ_cell_pairlist1,   &
         univ_mask2, univ_ix_natom, univ_ix_list,            &
         univ_iy_natom, univ_iy_list, univ_ij_sort_list,     &
         virial_check,                                       &
         fepgrp_pbc, fep_mask, table_sclj, table_scel,       &
         domain%num_atom_domain, MaxAtom,                    &
         MaxAtomCls, num_atom_cls,                           &            
         ncell_local, ncell_bound, ncell_max, cutoff_int,    &           
         univ_maxcell, univ_maxcell1, univ_ncell_nonzero,    &
         univ_ncell_near, univ_update, univ_mask2_size,      &
         univ_natom_max, npt, cpu_calc, density, cutoff2,    &
         pairlistdist2, univ_gpu_start,                      &
         system_size(1), system_size(2), system_size(3))

#endif

    return

  end subroutine compute_force_nonbond_tbl_lnr_gpu_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_notbl_gpu_fep
  !> @brief        calculate nonbonded force with gpu for FEP
  !! @authors      HO
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[in]    npt        : logical to check virial calculation is required
  !! @param[in]    cpu_calc   : flag for cpu calculation or not
  !! @param[inout] coord_pbc  : pbc oriented coordinates
  !! @param[inout] force_omp  : temporary forces of target system
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] ene_virial : energy&virial outputs
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_notbl_gpu_fep( &
                                                domain, enefunc, pairlist, &
                                                npt, cpu_calc, coord_pbc,  &
                                                force_omp, force, virial,  &
                                                ene_virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    logical,                  intent(in)    :: npt
    logical,                  intent(in)    :: cpu_calc
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force_omp(:,:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: ene_virial(:)

#ifdef USE_GPU

    ! local variables
    real(wp)                  :: cutoff2, pairlistdist2
    integer                   :: ncell, ncell_local, ncell_bound
    ! integer                   :: check_virial

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: nonb_lj6_factor(:)
    real(wp),         pointer :: density, cutoff, pairlistdist
    real(wp),         pointer :: table_grad(:)
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: atmcls(:,:)
    integer(int2),    pointer :: cell_pairlist(:,:)
    integer,          pointer :: natom(:), start_atom(:)
    integer(int1),    pointer :: virial_check(:,:)

    ! for GPU
    integer,          pointer :: univ_cell_pairlist1(:,:)
    integer,          pointer :: univ_cell_pairlist2(:,:)
    integer,          pointer :: univ_ix_natom(:)
    integer(1),       pointer :: univ_ix_list(:,:)
    integer,          pointer :: univ_iy_natom(:)
    integer(1),       pointer :: univ_iy_list(:,:)
    integer,          pointer :: univ_ij_sort_list(:)
    integer(1),       pointer :: univ_mask2(:,:)
    integer,          pointer :: univ_mask2_size
    integer                   :: univ_ncell_nonzero, univ_gpu_start_index
    integer                   :: univ_update
    integer                   :: ncell_max
    integer                   :: cutoff_int
    integer                   :: num_atom_cls
    integer                   :: ret_shape(1)
    integer                   :: univ_ij, univ_ij0, iiy, idx, base_idx
    integer                   :: ij, start_i, ix, i

    ! for FEP
    integer,          pointer :: fepgrp(:,:)
    integer(1),       pointer :: fep_mask(:,:)
    real(wp),         pointer :: table_sclj(:,:)
    real(wp),         pointer :: table_scel(:,:)
    integer(1),       pointer :: fepgrp_pbc(:)

    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6
    nonb_lj6_factor => enefunc%nonb_lj6_factor
    pairlistdist    => enefunc%pairlistdist

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    ncell_bound     =  domain%num_cell_boundary
    cutoff2         =  cutoff * cutoff
    pairlistdist2   =  pairlistdist * pairlistdist

    ret_shape = shape(natom)
    ncell_max = ret_shape(1)
    cutoff_int = enefunc%table%cutoff_int
    num_atom_cls = enefunc%num_atom_cls

    !
    ! univ
    univ_cell_pairlist1 => domain%univ_cell_pairlist1
    univ_cell_pairlist2 => domain%univ_cell_pairlist2
    univ_ix_natom       => pairlist%univ_ix_natom
    univ_ix_list        => pairlist%univ_ix_list
    univ_iy_natom       => pairlist%univ_iy_natom
    univ_iy_list        => pairlist%univ_iy_list
    univ_ij_sort_list   => pairlist%univ_ij_sort_list
    univ_mask2          => pairlist%univ_mask2
    univ_mask2_size     => pairlist%univ_mask2_size
    univ_ncell_nonzero  =  pairlist%univ_ncell_nonzero
    univ_update         =  pairlist%univ_update

    ! FEP
    fepgrp              => domain%fepgrp
    fep_mask            => enefunc%fep_mask
    table_sclj          => enefunc%table_sclj
    table_scel          => enefunc%table_scel
    fepgrp_pbc          => domain%fepgrp_pbc

    ij = domain%num_atom_domain
    !$omp parallel do private(i,ix,start_i)
    do i = 1, ncell
      start_i = start_atom(i)
      do ix = 1, natom(i)
        coord_pbc(     start_i+ix,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(  ij+start_i+ix,1,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(2*ij+start_i+ix,1,1) = coord(3,ix,i) + trans1(3,ix,i)
        fepgrp_pbc(start_i+ix) = fepgrp(ix,i)
      end do
    end do
    !$omp end parallel do

    !  launch GPU kernels (data transfer from CPU to GPU as well)
    !
    univ_gpu_start = ncell_local
    call gpu_launch_compute_force_nonbond_notable_univ_fep(  &
         coord_pbc, force(:,:,:,1), ene_virial,              &  
         cell_move, natom, start_atom, nonb_lj12, nonb_lj6,  &
         nonb_lj6_factor, table_grad, univ_cell_pairlist1,   &
         univ_mask2, univ_ix_natom, univ_ix_list,            &
         univ_iy_natom, univ_iy_list, univ_ij_sort_list,     &
         virial_check,                                       &
         fepgrp_pbc, fep_mask, table_sclj, table_scel,       &
         domain%num_atom_domain, MaxAtom,                    &
         MaxAtomCls, num_atom_cls,                           &  
         ncell_local, ncell_bound, ncell_max, cutoff_int,    & 
         univ_maxcell, univ_maxcell1, univ_ncell_nonzero,    &
         univ_ncell_near, univ_update, univ_mask2_size,      &
         univ_natom_max, npt, cpu_calc, density, cutoff2,    &
         pairlistdist2, univ_gpu_start,                      &
         system_size(1), system_size(2), system_size(3))

#endif

    return

  end subroutine compute_force_nonbond_notbl_gpu_fep

end module sp_energy_nonbond_fep_mod
