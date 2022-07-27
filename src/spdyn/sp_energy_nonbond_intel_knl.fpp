!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_nonbond_intel_knl_mod
!> @brief   calculate nonbonded energy with table and with linear interpolation
!! @authors Jaewoon Jung(JJ)
!  
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_nonbond_intel_knl_mod

  use sp_energy_nonbond_generic_mod !TODO
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
  public  :: compute_energy_nonbond_tbl_lnr_intel_knl
  public  :: compute_energy_nonbond_tbl_ljpme_intel_knl
  public  :: compute_force_nonbond_tbl_lnr_intel_knl
  public  :: compute_energy_nonbond_notbl_intel_knl
  public  :: compute_force_nonbond_notbl_intel_knl
  public  :: compute_force_nonbond_tbl_ljpme_intel_knl

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_tbl_lnr
  !> @brief        calculate nonbonded energy with lookup table
  !! @authors      JJ, NT
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] virial   : virial term of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : van der Waals energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_tbl_lnr_intel_knl( &
                                                 domain, enefunc, pairlist, &
                                                 coord_pbc, force, virial,  &
                                                 eelec, evdw)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)


    ! local variables
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: lj6, lj12
    real(wp)                  :: R
    real(wp)                  :: term_lj12, term_lj6, term_elec, factor
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: trans_x, trans_y, trans_z
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: elec_temp, evdw_temp
    real(wp)                  :: force_local(3)
    integer                   :: id, omp_get_thread_num
    integer                   :: ij, i, j, ix, iy, k, L, L1
    integer                   :: ncell, ncell_local
    integer                   :: check_virial

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: charge(:,:)
    real(wp),         pointer :: trans1(:,:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: atmcls(:,:)
    integer(int2),    pointer :: cell_pairlist(:,:)
    integer,          pointer :: natom(:)
    integer,          pointer :: nb15_list(:,:)
    integer(int1),    pointer :: virial_check(:,:)
    integer(1),       pointer :: exclusion_mask1(:,:,:), exclusion_mask(:,:,:)
    integer,          pointer :: nb15_cell(:)

#ifdef PKTIMER
    call timer_sta(222)
#ifdef FJ_PROF_FAPP
    call fapp_start("compute_energy_nonbond_table_linear",222,0)
#endif
#endif

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
    exclusion_mask1 => enefunc%exclusion_mask1
    exclusion_mask  => enefunc%exclusion_mask

    nb15_cell       => pairlist%nb15_cell
    nb15_list       => pairlist%nb15_list

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff


    !$omp parallel default(shared)                                      &
    !$omp private(id, i, ix, k, iy, L, ij, j, L1, rtmp, qtmp,  rij2, R, &
    !$omp         term_lj12, term_lj6, term_elec, grad_coef, work,      &
    !$omp         trans_x, trans_y, trans_z, force_local, lj6, lj12,    &
    !$omp         elec_temp, evdw_temp, dij, check_virial, factor)
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

      do ix = 1, natom(i) - 1

        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp = charge(ix,i)
        force_local(1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp

        do iy = ix+1, natom(i)

          ! compute distance
          !
          factor = real(exclusion_mask1(iy,ix,i),wp)
          dij(1) = rtmp(1) - coord_pbc(iy,1,i)
          dij(2) = rtmp(2) - coord_pbc(iy,2,i)
          dij(3) = rtmp(3) - coord_pbc(iy,3,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2  = cutoff2*density/rij2
          lj6  = nonb_lj6(atmcls(ix,i),atmcls(iy,i))
          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,i))
          L    = int(rij2)
          R    = rij2 - L
          L1   = 3*L - 2
          term_lj12 = table_ene(L1) + R*(table_ene(L1+3)-table_ene(L1))
          term_lj6  = table_ene(L1+1) + R*(table_ene(L1+4)-table_ene(L1+1))
          term_elec = table_ene(L1+2) + R*(table_ene(L1+5)-table_ene(L1+2))
          evdw_temp = evdw_temp + factor*(term_lj12*lj12-term_lj6*lj6)

          elec_temp = elec_temp + qtmp*charge(iy,i)*term_elec*factor
          term_lj12 = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6 &
                    + qtmp*charge(iy,i)*term_elec
          grad_coef = grad_coef*factor
          grad_coef = grad_coef
          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force(iy,1,i,id+1)  = force(iy,1,i,id+1) + work(1)
          force(iy,2,i,id+1)  = force(iy,2,i,id+1) + work(2)
          force(iy,3,i,id+1)  = force(iy,3,i,id+1) + work(3)

        end do

        force(ix,1,i,id+1) = force(ix,1,i,id+1) + force_local(1)
        force(ix,2,i,id+1) = force(ix,2,i,id+1) + force_local(2)
        force(ix,3,i,id+1) = force(ix,3,i,id+1) + force_local(3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp

      end do
    end do

    ! interaction between different cells (close)
    !
    do ij = id+1, maxcell_near, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)
      trans_x = real(cell_move(1,j,i),wp) * system_size(1)
      trans_y = real(cell_move(2,j,i),wp) * system_size(2)
      trans_z = real(cell_move(3,j,i),wp) * system_size(3)
      check_virial = virial_check(j,i)

      do ix = 1, natom(i)

        if (nb15_list(2*ix-1,ij) == 0) cycle

        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp = charge(ix,i)
        force_local(1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp

        do iy = nb15_list(2*ix-1,ij), nb15_list(2*ix,ij)

          factor = real(exclusion_mask(iy,ix,ij),wp)
          dij(1) = rtmp(1) - coord_pbc(iy,1,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(iy,2,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(iy,3,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2  = cutoff2*density/rij2
          lj6  = nonb_lj6 (atmcls(ix,i),atmcls(iy,j))
          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,j))
          L    = int(rij2)
          R    = rij2 - L
          L1   = 3*L - 2
          term_lj12 = table_ene(L1) + R*(table_ene(L1+3)-table_ene(L1))
          term_lj6  = table_ene(L1+1) + R*(table_ene(L1+4)-table_ene(L1+1))
          term_elec = table_ene(L1+2) + R*(table_ene(L1+5)-table_ene(L1+2))
          evdw_temp = evdw_temp + factor*(term_lj12*lj12 - term_lj6*lj6)
          elec_temp = elec_temp + qtmp*charge(iy,j)*term_elec*factor
          term_lj12 = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6 &
                    + qtmp*charge(iy,j)*term_elec
          grad_coef = grad_coef*factor
          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force(iy,1,j,id+1)  = force(iy,1,j,id+1) + work(1)
          force(iy,2,j,id+1)  = force(iy,2,j,id+1) + work(2)
          force(iy,3,j,id+1)  = force(iy,3,j,id+1) + work(3)

        end do

        if (check_virial == 1) &
          virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
        force(ix,1,i,id+1) = force(ix,1,i,id+1) + force_local(1)
        force(ix,2,i,id+1) = force(ix,2,i,id+1) + force_local(2)
        force(ix,3,i,id+1) = force(ix,3,i,id+1) + force_local(3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp

      end do
    end do

    do ij = id+1+maxcell_near, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)
      trans_x = real(cell_move(1,j,i),wp) * system_size(1)
      trans_y = real(cell_move(2,j,i),wp) * system_size(2)
      trans_z = real(cell_move(3,j,i),wp) * system_size(3)
      check_virial = virial_check(j,i)

      if (nb15_cell(ij) == 0) cycle

      do ix = 1, natom(i)

        if (nb15_list(2*ix-1,ij) == 0) cycle

        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp = charge(ix,i)
        force_local(1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp

        do iy = nb15_list(2*ix-1,ij), nb15_list(2*ix,ij)

          dij(1) = rtmp(1) - coord_pbc(iy,1,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(iy,2,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(iy,3,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 <= cutoff2) then

            rij2  = cutoff2*density/rij2
            lj6  = nonb_lj6 (atmcls(ix,i),atmcls(iy,j))
            lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,j))
            L    = int(rij2)
            R    = rij2 - L
            L1   = 3*L - 2
            term_lj12 = table_ene(L1) + R*(table_ene(L1+3)-table_ene(L1))
            term_lj6  = table_ene(L1+1) + R*(table_ene(L1+4)-table_ene(L1+1))
            term_elec = table_ene(L1+2) + R*(table_ene(L1+5)-table_ene(L1+2))
            evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
            elec_temp = elec_temp + qtmp*charge(iy,j)*term_elec
            term_lj12 = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))
            term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
            term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
            grad_coef = term_lj12*lj12 - term_lj6*lj6 &
                      + qtmp*charge(iy,j)*term_elec
            work(1) = grad_coef*dij(1)
            work(2) = grad_coef*dij(2)
            work(3) = grad_coef*dij(3)
  
            ! store force
            !
            force_local(1) = force_local(1) - work(1)
            force_local(2) = force_local(2) - work(2)
            force_local(3) = force_local(3) - work(3)
            force(iy,1,j,id+1)  = force(iy,1,j,id+1) + work(1)
            force(iy,2,j,id+1)  = force(iy,2,j,id+1) + work(2)
            force(iy,3,j,id+1)  = force(iy,3,j,id+1) + work(3)

          end if

        end do

        if (check_virial == 1) &
          virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
        force(ix,1,i,id+1) = force(ix,1,i,id+1) + force_local(1)
        force(ix,2,i,id+1) = force(ix,2,i,id+1) + force_local(2)
        force(ix,3,i,id+1) = force(ix,3,i,id+1) + force_local(3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp

      end do
    end do
    !$omp end parallel

#ifdef PKTIMER
    call timer_end(222)
#ifdef FJ_PROF_FAPP
    call fapp_stop("compute_energy_nonbond_table_linear",222,0)
#endif
#endif

    return

  end subroutine compute_energy_nonbond_tbl_lnr_intel_knl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_tbl_ljpme_intel_knl
  !> @brief        calculate nonbonded energy with ljpme
  !! @authors      JJ, NT
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] virial   : virial term of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : van der Waals energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_tbl_ljpme_intel_knl( &
                                                 domain, enefunc, pairlist, &
                                                 coord_pbc, force, virial,  &
                                                 eelec, evdw)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)


    ! local variables
    real(wp)                  :: dij(1:3), rij2, inv_r2, inv_r6, inv_r12
    real(wp)                  :: inv_cutoff2, inv_cutoff6
    real(wp)                  :: lj6, lj12, lj6_i, lj6_j, lj6_ij
    real(wp)                  :: R
    real(wp)                  :: term_lj12, term_lj6, term_elec, term_temp
    real(wp)                  :: factor
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: trans_x, trans_y, trans_z
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: elec_temp, evdw_temp
    real(wp)                  :: force_local(3)
    integer                   :: id, omp_get_thread_num
    integer                   :: ij, i, j, ix, iy, k, L, L1
    integer                   :: ncell, ncell_local
    integer                   :: iatmcls, jatmcls
    integer                   :: check_virial

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: charge(:,:)
    real(wp),         pointer :: trans1(:,:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: nonb_lj6_factor(:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: atmcls(:,:)
    integer(int2),    pointer :: cell_pairlist(:,:)
    integer,          pointer :: natom(:)
    integer,          pointer :: nb15_list(:,:)
    integer(int1),    pointer :: virial_check(:,:)
    integer(1),       pointer :: exclusion_mask1(:,:,:), exclusion_mask(:,:,:)
    integer,          pointer :: nb15_cell(:)


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
    nonb_lj6_factor => enefunc%nonb_lj6_factor
    exclusion_mask1 => enefunc%exclusion_mask1
    exclusion_mask  => enefunc%exclusion_mask

    nb15_cell       => pairlist%nb15_cell
    nb15_list       => pairlist%nb15_list

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff
    inv_cutoff2     =  1.0_wp /cutoff2
    inv_cutoff6     =  inv_cutoff2 * inv_cutoff2 * inv_cutoff2


    !$omp parallel default(shared)                                      &
    !$omp private(id, i, ix, k, iy, L, ij, j, L1, rtmp, qtmp,  rij2, R, &
    !$omp         term_lj12, term_lj6, term_elec, term_temp, grad_coef, &
    !$omp         work, trans_x, trans_y, trans_z, force_local, lj6,    &
    !$omp         lj6_i, lj6_j, lj6_ij, lj12, inv_r2, inv_r6, inv_r12,  &
    !$omp         elec_temp, evdw_temp, dij, check_virial, factor,      &
    !$omp         iatmcls, jatmcls)
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

      do ix = 1, natom(i) - 1

        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp    = charge(ix,i)
        iatmcls = atmcls(ix,i)
        lj6_i   = nonb_lj6_factor(iatmcls)
        force_local(1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp

        do iy = ix+1, natom(i)

          ! compute distance
          !
          factor   = real(exclusion_mask1(iy,ix,i),wp)
          dij(1)   = rtmp(1) - coord_pbc(iy,1,i)
          dij(2)   = rtmp(2) - coord_pbc(iy,2,i)
          dij(3)   = rtmp(3) - coord_pbc(iy,3,i)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          inv_r2   = 1.0_wp / rij2 
          inv_r6   = inv_r2 * inv_r2 * inv_r2
          inv_r12  = inv_r6 * inv_r6
          rij2     = cutoff2 * density * inv_r2
          jatmcls  = atmcls(iy,i)
          lj6_j    = nonb_lj6_factor(jatmcls)
          lj6_ij   = lj6_i * lj6_j
          lj6      = nonb_lj6 (jatmcls, iatmcls)
          lj12     = nonb_lj12(jatmcls, iatmcls)
          lj6      = lj6 -lj6_ij
          L        = int(rij2)
          R        = rij2 - L
          L1       = 2*L - 1
          term_temp = inv_r6 - inv_cutoff6
!         term_temp = inv_r6 
          term_lj6  = table_ene(L1  ) + R*(table_ene(L1+2)-table_ene(L1  ))
          term_elec = table_ene(L1+1) + R*(table_ene(L1+3)-table_ene(L1+1))
          evdw_temp = evdw_temp  &
                    + factor*(inv_r12*lj12-term_temp*lj6-term_lj6*lj6_ij)
          elec_temp = elec_temp + qtmp*charge(iy,i)*term_elec*factor
          term_lj12 = -12.0_wp*inv_r12*inv_r2
          term_temp = -6.0_wp*inv_r6*inv_r2
          term_lj6  = table_grad(L1  ) + R*(table_grad(L1+2)-table_grad(L1  ))
          term_elec = table_grad(L1+1) + R*(table_grad(L1+3)-table_grad(L1+1))
          grad_coef = term_lj12*lj12 - term_temp*lj6 - term_lj6*lj6_ij &
                    + qtmp*charge(iy,i)*term_elec
          grad_coef = grad_coef*factor
          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force(iy,1,i,id+1)  = force(iy,1,i,id+1) + work(1)
          force(iy,2,i,id+1)  = force(iy,2,i,id+1) + work(2)
          force(iy,3,i,id+1)  = force(iy,3,i,id+1) + work(3)

        end do

        force(ix,1,i,id+1) = force(ix,1,i,id+1) + force_local(1)
        force(ix,2,i,id+1) = force(ix,2,i,id+1) + force_local(2)
        force(ix,3,i,id+1) = force(ix,3,i,id+1) + force_local(3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp

      end do
    end do

    ! interaction between different cells (close)
    !
    do ij = id+1, maxcell_near, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)
      trans_x = real(cell_move(1,j,i),wp) * system_size(1)
      trans_y = real(cell_move(2,j,i),wp) * system_size(2)
      trans_z = real(cell_move(3,j,i),wp) * system_size(3)
      check_virial = virial_check(j,i)

      do ix = 1, natom(i)

        if (nb15_list(2*ix-1,ij) == 0) cycle

        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp    = charge(ix,i)
        iatmcls = atmcls(ix,i)
        lj6_i   = nonb_lj6_factor(iatmcls)
        force_local(1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp

        do iy = nb15_list(2*ix-1,ij), nb15_list(2*ix,ij)

          factor   = real(exclusion_mask(iy,ix,ij),wp)
          dij(1)   = rtmp(1) - coord_pbc(iy,1,j) + trans_x
          dij(2)   = rtmp(2) - coord_pbc(iy,2,j) + trans_y
          dij(3)   = rtmp(3) - coord_pbc(iy,3,j) + trans_z
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          inv_r2   = 1.0_wp / rij2
          inv_r6   = inv_r2 * inv_r2 * inv_r2
          inv_r12  = inv_r6 * inv_r6
          rij2     = cutoff2 * density * inv_r2
          jatmcls  = atmcls(iy,j)
          lj6_j    = nonb_lj6_factor(jatmcls)
          lj6_ij   = lj6_i * lj6_j
          lj6      = nonb_lj6 (jatmcls, iatmcls)
          lj12     = nonb_lj12(jatmcls, iatmcls)
          lj6      = lj6 - lj6_ij
          L        = int(rij2)
          R        = rij2 - L
          L1       = 2*L - 1
          term_temp = inv_r6 - inv_cutoff6
!         term_temp = inv_r6 
          term_lj6  = table_ene(L1  ) + R*(table_ene(L1+2)-table_ene(L1  ))
          term_elec = table_ene(L1+1) + R*(table_ene(L1+3)-table_ene(L1+1))
          evdw_temp = evdw_temp &
                    + factor*(inv_r12*lj12-term_temp*lj6-term_lj6*lj6_ij)
          elec_temp = elec_temp + qtmp*charge(iy,j)*term_elec*factor
          term_lj12 = -12.0_wp*inv_r12*inv_r2
          term_temp = -6.0_wp*inv_r6*inv_r2
          term_lj6  = table_grad(L1  ) + R*(table_grad(L1+2)-table_grad(L1  ))
          term_elec = table_grad(L1+1) + R*(table_grad(L1+3)-table_grad(L1+1))
          grad_coef = term_lj12*lj12 - term_temp*lj6 - term_lj6*lj6_ij &
                    + qtmp*charge(iy,j)*term_elec

          grad_coef = grad_coef*factor
          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force(iy,1,j,id+1)  = force(iy,1,j,id+1) + work(1)
          force(iy,2,j,id+1)  = force(iy,2,j,id+1) + work(2)
          force(iy,3,j,id+1)  = force(iy,3,j,id+1) + work(3)

        end do

        if (check_virial == 1) &
          virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
        force(ix,1,i,id+1) = force(ix,1,i,id+1) + force_local(1)
        force(ix,2,i,id+1) = force(ix,2,i,id+1) + force_local(2)
        force(ix,3,i,id+1) = force(ix,3,i,id+1) + force_local(3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp

      end do
    end do

    do ij = id+1+maxcell_near, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)
      trans_x = real(cell_move(1,j,i),wp) * system_size(1)
      trans_y = real(cell_move(2,j,i),wp) * system_size(2)
      trans_z = real(cell_move(3,j,i),wp) * system_size(3)
      check_virial = virial_check(j,i)

      if (nb15_cell(ij) == 0) cycle

      do ix = 1, natom(i)

        if (nb15_list(2*ix-1,ij) == 0) cycle

        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        iatmcls = atmcls(ix,i)
        lj6_i   = nonb_lj6_factor(iatmcls)
        qtmp    = charge(ix,i)
        force_local(1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp

        do iy = nb15_list(2*ix-1,ij), nb15_list(2*ix,ij)

          dij(1) = rtmp(1) - coord_pbc(iy,1,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(iy,2,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(iy,3,j) + trans_z
          rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 <= cutoff2) then

            inv_r2   = 1.0_wp / rij2
            inv_r6   = inv_r2 * inv_r2 * inv_r2
            inv_r12  = inv_r6 * inv_r6
            rij2     = cutoff2 * density * inv_r2
            jatmcls  = atmcls(iy,j)
            lj6_j    = nonb_lj6_factor(jatmcls)
            lj6_ij   = lj6_i * lj6_j
            lj6      = nonb_lj6 (jatmcls,iatmcls)
            lj12     = nonb_lj12(jatmcls,iatmcls)
            lj6      = lj6 - lj6_ij
            L        = int(rij2)
            R        = rij2 - L
            L1       = 2*L - 1
            term_temp = inv_r6 - inv_cutoff6
!           term_temp = inv_r6
            term_lj6  = table_ene(L1  ) + R*(table_ene(L1+2)-table_ene(L1  ))
            term_elec = table_ene(L1+1) + R*(table_ene(L1+3)-table_ene(L1+1))
            evdw_temp = evdw_temp &
                      + inv_r12*lj12 - term_temp*lj6 - term_lj6*lj6_ij
            elec_temp = elec_temp + qtmp*charge(iy,j)*term_elec
            term_lj12 = -12.0_wp*inv_r12*inv_r2
            term_temp = -6.0_wp*inv_r6*inv_r2
            term_lj6  = table_grad(L1  ) + R*(table_grad(L1+2)-table_grad(L1  ))
            term_elec = table_grad(L1+1) + R*(table_grad(L1+3)-table_grad(L1+1))
            grad_coef = term_lj12*lj12 - term_temp*lj6 - term_lj6*lj6_ij &
                      + qtmp*charge(iy,j)*term_elec
            work(1) = grad_coef*dij(1)
            work(2) = grad_coef*dij(2)
            work(3) = grad_coef*dij(3)
  
            ! store force
            !
            force_local(1) = force_local(1) - work(1)
            force_local(2) = force_local(2) - work(2)
            force_local(3) = force_local(3) - work(3)
            force(iy,1,j,id+1)  = force(iy,1,j,id+1) + work(1)
            force(iy,2,j,id+1)  = force(iy,2,j,id+1) + work(2)
            force(iy,3,j,id+1)  = force(iy,3,j,id+1) + work(3)

          end if

        end do

        if (check_virial == 1) &
          virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
        force(ix,1,i,id+1) = force(ix,1,i,id+1) + force_local(1)
        force(ix,2,i,id+1) = force(ix,2,i,id+1) + force_local(2)
        force(ix,3,i,id+1) = force(ix,3,i,id+1) + force_local(3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp

      end do
    end do
    !$omp end parallel

    return

  end subroutine compute_energy_nonbond_tbl_ljpme_intel_knl


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_tbl_lnr_intel_knl
  !> @brief        calculate nonbonded force without solvents with lookup table
  !! @authors      JJ 
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] virial   : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_tbl_lnr_intel_knl( &
                                                 domain, enefunc, pairlist, &
                                                 coord_pbc, force, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:)

    ! local variables
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: R, lj12, lj6
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, grad_coef, factor
    real(wp)                  :: work(1:3)
    real(wp)                  :: trans_x, trans_y, trans_z
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: force_local(3)
    integer                   :: i, ix, iy, j, k, ij, L, L1
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
    integer(int1),    pointer :: virial_check(:,:)
    integer(1),       pointer :: exclusion_mask1(:,:,:), exclusion_mask(:,:,:)


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
    exclusion_mask1 => enefunc%exclusion_mask1
    exclusion_mask  => enefunc%exclusion_mask

    nb15_list       => pairlist%nb15_list
    nb15_cell       => pairlist%nb15_cell

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff

    ! calculate energy and gradient
    !
#ifdef FJ_TIMER_2
      call timer_sta(412)
#endif

    !$omp parallel default(shared)                                     &
    !$omp private(id, i, ix, k, iy, ij, j, L, L1, iatmcls, jatmcls,    &
    !$omp         rtmp, qtmp, rij2, R, term_lj12, term_lj6, grad_coef, &
    !$omp         work, term_elec, trans_x, trans_y, trans_z,          &
    !$omp         force_local, lj12, lj6, dij, check_virial, factor)
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

    do i = id+1, ncell_local, nthread

      do ix = 1, natom(i) - 1

        iatmcls = atmcls(ix,i)
        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp = charge(ix,i)
        force_local(1:3) = 0.0_wp

        do iy = ix+1, natom(i)

          ! compute distance
          !
          jatmcls = atmcls(iy,i)
          factor = real(exclusion_mask1(iy,ix,i),wp)
          dij(1) = rtmp(1) - coord_pbc(iy,1,i)
          dij(2) = rtmp(2) - coord_pbc(iy,2,i)
          dij(3) = rtmp(3) - coord_pbc(iy,3,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2  = cutoff2*density/rij2
          lj6  = nonb_lj6(jatmcls,iatmcls)
          lj12 = nonb_lj12(jatmcls,iatmcls)
          L    = int(rij2)
          R    = rij2 - L
          L1   = 3*L - 2
          term_lj12 = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6 &
                    + qtmp*charge(iy,i)*term_elec
          grad_coef = grad_coef * factor
          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force(iy,1,i,id+1) = force(iy,1,i,id+1) + work(1)
          force(iy,2,i,id+1) = force(iy,2,i,id+1) + work(2)
          force(iy,3,i,id+1) = force(iy,3,i,id+1) + work(3)

        end do

        force(ix,1,i,id+1) = force(ix,1,i,id+1) + force_local(1)
        force(ix,2,i,id+1) = force(ix,2,i,id+1) + force_local(2)
        force(ix,3,i,id+1) = force(ix,3,i,id+1) + force_local(3)

      end do
    end do

    ! interaction between different cells 
    !
    do ij = id+1, maxcell_near, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)
      trans_x = real(cell_move(1,j,i),wp) * system_size(1)
      trans_y = real(cell_move(2,j,i),wp) * system_size(2)
      trans_z = real(cell_move(3,j,i),wp) * system_size(3)
      check_virial = virial_check(j,i)

      do ix = 1, natom(i)

        if (nb15_list(2*ix-1,ij) == 0) cycle
        iatmcls = atmcls(ix,i)
        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp = charge(ix,i)
        force_local(1:3) = 0.0_wp

        do iy = nb15_list(2*ix-1,ij), nb15_list(2*ix,ij)

          jatmcls = atmcls(iy,j)
          factor = real(exclusion_mask(iy,ix,ij),wp)
          dij(1) = rtmp(1) - coord_pbc(iy,1,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(iy,2,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(iy,3,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2  = cutoff2*density/rij2
          lj6  = nonb_lj6(jatmcls,iatmcls)
          lj12 = nonb_lj12(jatmcls,iatmcls)
          L    = int(rij2)
          R    = rij2 - L
          L1   = 3*L - 2
          term_lj12 = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6 &
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
        force(ix,1,i,id+1) = force(ix,1,i,id+1) + force_local(1)
        force(ix,2,i,id+1) = force(ix,2,i,id+1) + force_local(2)
        force(ix,3,i,id+1) = force(ix,3,i,id+1) + force_local(3)

      end do
    end do

    do ij = id+1+maxcell_near, maxcell, nthread

      if (nb15_cell(ij) == 0) cycle

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)
      trans_x = real(cell_move(1,j,i),wp) * system_size(1)
      trans_y = real(cell_move(2,j,i),wp) * system_size(2)
      trans_z = real(cell_move(3,j,i),wp) * system_size(3)
      check_virial = virial_check(j,i)

      do ix = 1, natom(i)

        if (nb15_list(2*ix-1,ij) == 0) cycle

        iatmcls = atmcls(ix,i)
        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp = charge(ix,i)
        force_local(1:3) = 0.0_wp

        do iy = nb15_list(2*ix-1,ij), nb15_list(2*ix,ij)

          dij(1) = rtmp(1) - coord_pbc(iy,1,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(iy,2,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(iy,3,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 <= cutoff2) then

            jatmcls = atmcls(iy,j)
            rij2  = cutoff2*density/rij2
            lj6  = nonb_lj6(jatmcls,iatmcls)
            lj12 = nonb_lj12(jatmcls,iatmcls)
            L    = int(rij2)
            R    = rij2 - L
            L1   = 3*L - 2
            term_lj12 = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))
            term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
            term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
            grad_coef = term_lj12*lj12 - term_lj6*lj6 &
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

          end if

        end do

        if (check_virial == 1) &
          virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
        force(ix,1,i,id+1) = force(ix,1,i,id+1) + force_local(1)
        force(ix,2,i,id+1) = force(ix,2,i,id+1) + force_local(2)
        force(ix,3,i,id+1) = force(ix,3,i,id+1) + force_local(3)

      end do
    end do
    !$omp end parallel

#ifdef FJ_TIMER_2
      call timer_end(412)
#endif

    return

  end subroutine compute_force_nonbond_tbl_lnr_intel_knl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_notbl_lnr_intel_knl
  !> @brief        calculate nonbonded energy without lookup table of vdw
  !! @authors      JJ, NT
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] virial   : virial term of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : van der Waals energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_notbl_intel_knl( &
                                                 domain, enefunc, pairlist, &
                                                 coord_pbc, force, virial,  &
                                                 eelec, evdw)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)


    ! local variables
    real(wp)                  :: dij(1:3), rij2, rij2_inv
    real(wp)                  :: lj6, lj12
    real(wp)                  :: R
    real(wp)                  :: term_lj12, term_lj6, term_elec, factor
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: trans_x, trans_y, trans_z
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: elec_temp, evdw_temp
    real(wp)                  :: force_local(3)
    integer                   :: id, omp_get_thread_num
    integer                   :: ij, i, j, ix, iy, k, L, L1
    integer                   :: ncell, ncell_local
    integer                   :: check_virial

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: charge(:,:)
    real(wp),         pointer :: trans1(:,:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: atmcls(:,:)
    integer(int2),    pointer :: cell_pairlist(:,:)
    integer,          pointer :: natom(:)
    integer,          pointer :: nb15_list(:,:)
    integer(int1),    pointer :: virial_check(:,:)
    integer(1),       pointer :: exclusion_mask1(:,:,:), exclusion_mask(:,:,:)


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
    exclusion_mask1 => enefunc%exclusion_mask1
    exclusion_mask  => enefunc%exclusion_mask

    nb15_list       => pairlist%nb15_list

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff


    !$omp parallel default(shared)                                      &
    !$omp private(id, i, ix, k, iy, L, ij, j, L1, rtmp, qtmp,  rij2, R, &
    !$omp         term_lj12, term_lj6, term_elec, grad_coef, work,      &
    !$omp         trans_x, trans_y, trans_z, force_local, lj6, lj12,    &
    !$omp         elec_temp, evdw_temp, dij, check_virial, factor,      &
    !$omp         rij2_inv)
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

      do ix = 1, natom(i) - 1

        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp = charge(ix,i)
        force_local(1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp

        do iy = ix+1, natom(i)

          ! compute distance
          !
          factor = real(exclusion_mask1(iy,ix,i),wp)
          dij(1) = rtmp(1) - coord_pbc(iy,1,i)
          dij(2) = rtmp(2) - coord_pbc(iy,2,i)
          dij(3) = rtmp(3) - coord_pbc(iy,3,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          if (rij2 > cutoff2) factor = 0.0_wp
          rij2_inv = 1.0_wp / rij2
          rij2  = cutoff2*density*rij2_inv
          lj6  = nonb_lj6(atmcls(ix,i),atmcls(iy,i))
          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,i))
          L    = int(rij2)
          R    = rij2 - L
          term_lj6  = rij2_inv * rij2_inv * rij2_inv
          term_lj12 = term_lj6 * term_lj6
          term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
          evdw_temp = evdw_temp + factor*(term_lj12*lj12-term_lj6*lj6)
          elec_temp = elec_temp + qtmp*charge(iy,i)*term_elec*factor
          term_lj12 = -12.0_wp * term_lj12 * rij2_inv
          term_lj6  = - 6.0_wp * term_lj6  * rij2_inv
          term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
          grad_coef = term_lj12*lj12 - term_lj6*lj6 &
                    + qtmp*charge(iy,i)*term_elec
          grad_coef = grad_coef*factor
          grad_coef = grad_coef
          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force(iy,1,i,id+1)  = force(iy,1,i,id+1) + work(1)
          force(iy,2,i,id+1)  = force(iy,2,i,id+1) + work(2)
          force(iy,3,i,id+1)  = force(iy,3,i,id+1) + work(3)

        end do

        force(ix,1,i,id+1) = force(ix,1,i,id+1) + force_local(1)
        force(ix,2,i,id+1) = force(ix,2,i,id+1) + force_local(2)
        force(ix,3,i,id+1) = force(ix,3,i,id+1) + force_local(3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp

      end do
    end do

    ! interaction between different cells (close)
    !
    do ij = id+1, maxcell_near, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)
      trans_x = real(cell_move(1,j,i),wp) * system_size(1)
      trans_y = real(cell_move(2,j,i),wp) * system_size(2)
      trans_z = real(cell_move(3,j,i),wp) * system_size(3)
      check_virial = virial_check(j,i)

      do ix = 1, natom(i)

        if (nb15_list(2*ix-1,ij) == 0) cycle

        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp = charge(ix,i)
        force_local(1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp

        do iy = nb15_list(2*ix-1,ij), nb15_list(2*ix,ij)

          factor = real(exclusion_mask(iy,ix,ij),wp)
          dij(1) = rtmp(1) - coord_pbc(iy,1,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(iy,2,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(iy,3,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          if (rij2 > cutoff2) factor = 0.0_wp
          rij2_inv = 1.0_wp / rij2
          rij2  = cutoff2*density*rij2_inv
          lj6  = nonb_lj6 (atmcls(ix,i),atmcls(iy,j))
          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,j))
          L    = int(rij2)
          R    = rij2 - L
          term_lj6  = rij2_inv * rij2_inv * rij2_inv
          term_lj12 = term_lj6 * term_lj6
          term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
          evdw_temp = evdw_temp + factor*(term_lj12*lj12 - term_lj6*lj6)
          elec_temp = elec_temp + qtmp*charge(iy,j)*term_elec*factor
          term_lj12 = -12.0_wp * term_lj12 * rij2_inv
          term_lj6  = - 6.0_wp * term_lj6  * rij2_inv
          term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
          grad_coef = term_lj12*lj12 - term_lj6*lj6 &
                    + qtmp*charge(iy,j)*term_elec
          grad_coef = grad_coef*factor
          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force(iy,1,j,id+1)  = force(iy,1,j,id+1) + work(1)
          force(iy,2,j,id+1)  = force(iy,2,j,id+1) + work(2)
          force(iy,3,j,id+1)  = force(iy,3,j,id+1) + work(3)

        end do

        if (check_virial == 1) &
          virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
        force(ix,1,i,id+1) = force(ix,1,i,id+1) + force_local(1)
        force(ix,2,i,id+1) = force(ix,2,i,id+1) + force_local(2)
        force(ix,3,i,id+1) = force(ix,3,i,id+1) + force_local(3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp

      end do
    end do

    do ij = id+1+maxcell_near, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)
      trans_x = real(cell_move(1,j,i),wp) * system_size(1)
      trans_y = real(cell_move(2,j,i),wp) * system_size(2)
      trans_z = real(cell_move(3,j,i),wp) * system_size(3)
      check_virial = virial_check(j,i)

      do ix = 1, natom(i)

        if (nb15_list(2*ix-1,ij) == 0) cycle

        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp = charge(ix,i)
        force_local(1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp

        do iy = nb15_list(2*ix-1,ij), nb15_list(2*ix,ij)

          dij(1) = rtmp(1) - coord_pbc(iy,1,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(iy,2,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(iy,3,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 <= cutoff2) then

            rij2_inv = 1.0_wp / rij2
            rij2  = cutoff2*density*rij2_inv
            lj6  = nonb_lj6 (atmcls(ix,i),atmcls(iy,j))
            lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,j))
            L    = int(rij2)
            R    = rij2 - L
            term_lj6  = rij2_inv * rij2_inv * rij2_inv
            term_lj12 = term_lj6 * term_lj6
            term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
            evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
            elec_temp = elec_temp + qtmp*charge(iy,j)*term_elec
            term_lj12 = -12.0_wp * term_lj12 * rij2_inv
            term_lj6  = - 6.0_wp * term_lj6  * rij2_inv
            term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
            grad_coef = term_lj12*lj12 - term_lj6*lj6 &
                      + qtmp*charge(iy,j)*term_elec
            work(1) = grad_coef*dij(1)
            work(2) = grad_coef*dij(2)
            work(3) = grad_coef*dij(3)

            ! store force
            !
            force_local(1) = force_local(1) - work(1)
            force_local(2) = force_local(2) - work(2)
            force_local(3) = force_local(3) - work(3)
            force(iy,1,j,id+1)  = force(iy,1,j,id+1) + work(1)
            force(iy,2,j,id+1)  = force(iy,2,j,id+1) + work(2)
            force(iy,3,j,id+1)  = force(iy,3,j,id+1) + work(3)

          end if
  
        end do

        if (check_virial == 1) &
          virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
        force(ix,1,i,id+1) = force(ix,1,i,id+1) + force_local(1)
        force(ix,2,i,id+1) = force(ix,2,i,id+1) + force_local(2)
        force(ix,3,i,id+1) = force(ix,3,i,id+1) + force_local(3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp

      end do
    end do
    !$omp end parallel

    return

  end subroutine compute_energy_nonbond_notbl_intel_knl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_notbl_intel_knl
  !> @brief        calculate nonbonded force without solvents with lookup table
  !! @authors      JJ 
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] virial   : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_notbl_intel_knl( &
                                                 domain, enefunc, pairlist, &
                                                 coord_pbc, force, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:)

    ! local variables
    real(wp)                  :: dij(1:3), rij2, rij2_inv
    real(wp)                  :: R, lj12, lj6
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, grad_coef, factor
    real(wp)                  :: work(1:3)
    real(wp)                  :: trans_x, trans_y, trans_z
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: force_local(3)
    integer                   :: i, ix, iy, j, k, ij, L, L1
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
    integer(int1),    pointer :: virial_check(:,:)
    integer(1),       pointer :: exclusion_mask1(:,:,:), exclusion_mask(:,:,:)


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
    exclusion_mask1 => enefunc%exclusion_mask1
    exclusion_mask  => enefunc%exclusion_mask

    nb15_list       => pairlist%nb15_list
    nb15_cell       => pairlist%nb15_cell

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff

    ! calculate energy and gradient
    !

    !$omp parallel default(shared)                                     &
    !$omp private(id, i, ix, k, iy, ij, j, L, L1, iatmcls, jatmcls,    &
    !$omp         rtmp, qtmp, rij2, R, term_lj12, term_lj6, grad_coef, &
    !$omp         work, term_elec, trans_x, trans_y, trans_z,          &
    !$omp         force_local, lj12, lj6, dij, check_virial, factor,   &
    !$omp         rij2_inv)
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

    do i = id+1, ncell_local, nthread

      do ix = 1, natom(i) - 1

        iatmcls = atmcls(ix,i)
        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp = charge(ix,i)
        force_local(1:3) = 0.0_wp

        do iy = ix+1, natom(i)

          ! compute distance
          !
          jatmcls = atmcls(iy,i)
          factor = real(exclusion_mask1(iy,ix,i),wp)
          dij(1) = rtmp(1) - coord_pbc(iy,1,i)
          dij(2) = rtmp(2) - coord_pbc(iy,2,i)
          dij(3) = rtmp(3) - coord_pbc(iy,3,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          if (rij2 > cutoff2) factor = 0.0_wp
          rij2_inv = 1.0_wp / rij2
          rij2  = cutoff2*density*rij2_inv
          lj6  = nonb_lj6(jatmcls,iatmcls)
          lj12 = nonb_lj12(jatmcls,iatmcls)
          L    = int(rij2)
          R    = rij2 - L
          term_lj6  = rij2_inv * rij2_inv * rij2_inv
          term_lj12 = term_lj6 * term_lj6
          term_lj12 = -12.0_wp * term_lj12 * rij2_inv
          term_lj6  = - 6.0_wp * term_lj6  * rij2_inv
          term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
          grad_coef = term_lj12*lj12 - term_lj6*lj6 &
                    + qtmp*charge(iy,i)*term_elec
          grad_coef = grad_coef * factor
          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force(iy,1,i,id+1) = force(iy,1,i,id+1) + work(1)
          force(iy,2,i,id+1) = force(iy,2,i,id+1) + work(2)
          force(iy,3,i,id+1) = force(iy,3,i,id+1) + work(3)

        end do

        force(ix,1,i,id+1) = force(ix,1,i,id+1) + force_local(1)
        force(ix,2,i,id+1) = force(ix,2,i,id+1) + force_local(2)
        force(ix,3,i,id+1) = force(ix,3,i,id+1) + force_local(3)

      end do
    end do

    ! interaction between different cells 
    !
    do ij = id+1, maxcell_near, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)
      trans_x = real(cell_move(1,j,i),wp) * system_size(1)
      trans_y = real(cell_move(2,j,i),wp) * system_size(2)
      trans_z = real(cell_move(3,j,i),wp) * system_size(3)
      check_virial = virial_check(j,i)

      do ix = 1, natom(i)

        if (nb15_list(2*ix-1,ij) == 0) cycle
        iatmcls = atmcls(ix,i)
        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp = charge(ix,i)
        force_local(1:3) = 0.0_wp

        do iy = nb15_list(2*ix-1,ij), nb15_list(2*ix,ij)

          jatmcls = atmcls(iy,j)
          factor = real(exclusion_mask(iy,ix,ij),wp)
          dij(1) = rtmp(1) - coord_pbc(iy,1,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(iy,2,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(iy,3,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          if (rij2 > cutoff2) factor = 0.0_wp
          rij2_inv = 1.0_wp / rij2
          rij2  = cutoff2*density*rij2_inv
          lj6  = nonb_lj6(jatmcls,iatmcls)
          lj12 = nonb_lj12(jatmcls,iatmcls)
          L    = int(rij2)
          R    = rij2 - L
          term_lj6  = rij2_inv * rij2_inv * rij2_inv
          term_lj12 = term_lj6 * term_lj6
          term_lj12 = -12.0_wp * term_lj12 * rij2_inv
          term_lj6  = - 6.0_wp * term_lj6  * rij2_inv
          term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
          grad_coef = term_lj12*lj12 - term_lj6*lj6 &
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
        force(ix,1,i,id+1) = force(ix,1,i,id+1) + force_local(1)
        force(ix,2,i,id+1) = force(ix,2,i,id+1) + force_local(2)
        force(ix,3,i,id+1) = force(ix,3,i,id+1) + force_local(3)

      end do
    end do

    do ij = id+1+maxcell_near, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)
      trans_x = real(cell_move(1,j,i),wp) * system_size(1)
      trans_y = real(cell_move(2,j,i),wp) * system_size(2)
      trans_z = real(cell_move(3,j,i),wp) * system_size(3)
      check_virial = virial_check(j,i)

      do ix = 1, natom(i)

        if (nb15_list(2*ix-1,ij) == 0) cycle

        iatmcls = atmcls(ix,i)
        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp = charge(ix,i)
        force_local(1:3) = 0.0_wp

        do iy = nb15_list(2*ix-1,ij), nb15_list(2*ix,ij)

          jatmcls = atmcls(iy,j)
          dij(1) = rtmp(1) - coord_pbc(iy,1,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(iy,2,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(iy,3,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 <= cutoff2) then

            rij2_inv = 1.0_wp / rij2
            rij2  = cutoff2*density*rij2_inv
            lj6  = nonb_lj6(jatmcls,iatmcls)
            lj12 = nonb_lj12(jatmcls,iatmcls)
            L    = int(rij2)
            R    = rij2 - L
            L1   = 3*L - 2
            term_lj6  = rij2_inv * rij2_inv * rij2_inv
            term_lj12 = term_lj6 * term_lj6
            term_lj12 = -12.0_wp * term_lj12 * rij2_inv
            term_lj6  = - 6.0_wp * term_lj6  * rij2_inv
            term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
            grad_coef = term_lj12*lj12 - term_lj6*lj6 &
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

          end if

        end do

        if (check_virial == 1) &
          virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
        force(ix,1,i,id+1) = force(ix,1,i,id+1) + force_local(1)
        force(ix,2,i,id+1) = force(ix,2,i,id+1) + force_local(2)
        force(ix,3,i,id+1) = force(ix,3,i,id+1) + force_local(3)

      end do
    end do
    !$omp end parallel

    return

  end subroutine compute_force_nonbond_notbl_intel_knl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_tbl_ljpme_intel_knl
  !> @brief        calculate nonbonded force with ljpme
  !! @authors      JJ, NT
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] virial   : virial term of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : van der Waals energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_tbl_ljpme_intel_knl( &
                                                 domain, enefunc, pairlist, &
                                                 coord_pbc, force, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:)


    ! local variables
    real(wp)                  :: dij(1:3), rij2, inv_r2, inv_r6, inv_r12
    real(wp)                  :: lj6, lj12, lj6_i, lj6_j, lj6_ij
    real(wp)                  :: R
    real(wp)                  :: term_lj12, term_lj6, term_elec, term_temp
    real(wp)                  :: factor
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: trans_x, trans_y, trans_z
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: elec_temp, evdw_temp
    real(wp)                  :: force_local(3)
    integer                   :: id, omp_get_thread_num
    integer                   :: ij, i, j, ix, iy, k, L, L1
    integer                   :: ncell, ncell_local
    integer                   :: iatmcls, jatmcls
    integer                   :: check_virial

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: charge(:,:)
    real(wp),         pointer :: trans1(:,:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: nonb_lj6_factor(:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: atmcls(:,:)
    integer(int2),    pointer :: cell_pairlist(:,:)
    integer,          pointer :: natom(:)
    integer,          pointer :: nb15_list(:,:)
    integer(int1),    pointer :: virial_check(:,:)
    integer(1),       pointer :: exclusion_mask1(:,:,:), exclusion_mask(:,:,:)
    integer,          pointer :: nb15_cell(:)


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
    nonb_lj6_factor => enefunc%nonb_lj6_factor
    exclusion_mask1 => enefunc%exclusion_mask1
    exclusion_mask  => enefunc%exclusion_mask

    nb15_cell       => pairlist%nb15_cell
    nb15_list       => pairlist%nb15_list

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff


    !$omp parallel default(shared)                                      &
    !$omp private(id, i, ix, k, iy, L, ij, j, L1, rtmp, qtmp,  rij2, R, &
    !$omp         term_lj12, term_lj6, term_elec, term_temp, grad_coef, &
    !$omp         work, trans_x, trans_y, trans_z, force_local, lj6,    &
    !$omp         lj6_i, lj6_j, lj6_ij, lj12, inv_r2, inv_r6, inv_r12,  &
    !$omp         elec_temp, evdw_temp, dij, check_virial, factor,      &
    !$omp         iatmcls, jatmcls)
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

      do ix = 1, natom(i) - 1

        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp    = charge(ix,i)
        iatmcls = atmcls(ix,i)
        lj6_i   = nonb_lj6_factor(iatmcls)
        force_local(1:3) = 0.0_wp

        do iy = ix+1, natom(i)

          ! compute distance
          !
          factor   = real(exclusion_mask1(iy,ix,i),wp)
          dij(1)   = rtmp(1) - coord_pbc(iy,1,i)
          dij(2)   = rtmp(2) - coord_pbc(iy,2,i)
          dij(3)   = rtmp(3) - coord_pbc(iy,3,i)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          inv_r2   = 1.0_wp / rij2 
          inv_r6   = inv_r2 * inv_r2 * inv_r2
          inv_r12  = inv_r6 * inv_r6
          rij2     = cutoff2 * density * inv_r2
          jatmcls  = atmcls(iy,i)
          lj6_j    = nonb_lj6_factor(jatmcls)
          lj6_ij   = lj6_i * lj6_j
          lj6      = nonb_lj6 (jatmcls, iatmcls)
          lj12     = nonb_lj12(jatmcls, iatmcls)
          lj6      = lj6 -lj6_ij
          L        = int(rij2)
          R        = rij2 - L
          L1       = 2*L - 1
          term_lj12 = -12.0_wp*inv_r12*inv_r2
          term_temp = -6.0_wp*inv_r6*inv_r2
          term_lj6  = table_grad(L1  ) + R*(table_grad(L1+2)-table_grad(L1  ))
          term_elec = table_grad(L1+1) + R*(table_grad(L1+3)-table_grad(L1+1))
          grad_coef = term_lj12*lj12 - term_temp*lj6 - term_lj6*lj6_ij &
                    + qtmp*charge(iy,i)*term_elec
          grad_coef = grad_coef*factor
          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force(iy,1,i,id+1)  = force(iy,1,i,id+1) + work(1)
          force(iy,2,i,id+1)  = force(iy,2,i,id+1) + work(2)
          force(iy,3,i,id+1)  = force(iy,3,i,id+1) + work(3)

        end do

        force(ix,1,i,id+1) = force(ix,1,i,id+1) + force_local(1)
        force(ix,2,i,id+1) = force(ix,2,i,id+1) + force_local(2)
        force(ix,3,i,id+1) = force(ix,3,i,id+1) + force_local(3)

      end do
    end do

    ! interaction between different cells (close)
    !
    do ij = id+1, maxcell_near, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)
      trans_x = real(cell_move(1,j,i),wp) * system_size(1)
      trans_y = real(cell_move(2,j,i),wp) * system_size(2)
      trans_z = real(cell_move(3,j,i),wp) * system_size(3)
      check_virial = virial_check(j,i)

      do ix = 1, natom(i)

        if (nb15_list(2*ix-1,ij) == 0) cycle

        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp    = charge(ix,i)
        iatmcls = atmcls(ix,i)
        lj6_i   = nonb_lj6_factor(iatmcls)
        force_local(1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp

        do iy = nb15_list(2*ix-1,ij), nb15_list(2*ix,ij)

          factor   = real(exclusion_mask(iy,ix,ij),wp)
          dij(1)   = rtmp(1) - coord_pbc(iy,1,j) + trans_x
          dij(2)   = rtmp(2) - coord_pbc(iy,2,j) + trans_y
          dij(3)   = rtmp(3) - coord_pbc(iy,3,j) + trans_z
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          inv_r2   = 1.0_wp / rij2
          inv_r6   = inv_r2 * inv_r2 * inv_r2
          inv_r12  = inv_r6 * inv_r6
          rij2     = cutoff2 * density * inv_r2
          jatmcls  = atmcls(iy,j)
          lj6_j    = nonb_lj6_factor(jatmcls)
          lj6_ij   = lj6_i * lj6_j
          lj6      = nonb_lj6 (jatmcls, iatmcls)
          lj12     = nonb_lj12(jatmcls, iatmcls)
          lj6      = lj6 - lj6_ij
          L        = int(rij2)
          R        = rij2 - L
          L1       = 2*L - 1
          term_lj12 = -12.0_wp*inv_r12*inv_r2
          term_temp = -6.0_wp*inv_r6*inv_r2
          term_lj6  = table_grad(L1  ) + R*(table_grad(L1+2)-table_grad(L1  ))
          term_elec = table_grad(L1+1) + R*(table_grad(L1+3)-table_grad(L1+1))
          grad_coef = term_lj12*lj12 - term_temp*lj6 - term_lj6*lj6_ij &
                    + qtmp*charge(iy,j)*term_elec
          grad_coef = grad_coef*factor
          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force(iy,1,j,id+1)  = force(iy,1,j,id+1) + work(1)
          force(iy,2,j,id+1)  = force(iy,2,j,id+1) + work(2)
          force(iy,3,j,id+1)  = force(iy,3,j,id+1) + work(3)

        end do

        if (check_virial == 1) &
          virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
        force(ix,1,i,id+1) = force(ix,1,i,id+1) + force_local(1)
        force(ix,2,i,id+1) = force(ix,2,i,id+1) + force_local(2)
        force(ix,3,i,id+1) = force(ix,3,i,id+1) + force_local(3)

      end do
    end do

    do ij = id+1+maxcell_near, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)
      trans_x = real(cell_move(1,j,i),wp) * system_size(1)
      trans_y = real(cell_move(2,j,i),wp) * system_size(2)
      trans_z = real(cell_move(3,j,i),wp) * system_size(3)
      check_virial = virial_check(j,i)

      if (nb15_cell(ij) == 0) cycle

      do ix = 1, natom(i)

        if (nb15_list(2*ix-1,ij) == 0) cycle

        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        iatmcls = atmcls(ix,i)
        lj6_i   = nonb_lj6_factor(iatmcls)
        qtmp    = charge(ix,i)
        force_local(1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp

        do iy = nb15_list(2*ix-1,ij), nb15_list(2*ix,ij)

          dij(1) = rtmp(1) - coord_pbc(iy,1,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(iy,2,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(iy,3,j) + trans_z
          rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 <= cutoff2) then

            inv_r2   = 1.0_wp / rij2
            inv_r6   = inv_r2 * inv_r2 * inv_r2
            inv_r12  = inv_r6 * inv_r6
            rij2     = cutoff2 * density * inv_r2
            jatmcls  = atmcls(iy,j)
            lj6_j    = nonb_lj6_factor(jatmcls)
            lj6_ij   = lj6_i * lj6_j
            lj6      = nonb_lj6 (jatmcls,iatmcls)
            lj12     = nonb_lj12(jatmcls,iatmcls)
            lj6      = lj6 - lj6_ij
            L        = int(rij2)
            R        = rij2 - L
            L1       = 2*L - 1
            term_lj12 = -12.0_wp * inv_r12 * inv_r2
            term_temp = -6.0_wp * inv_r6 * inv_r2
            term_lj6  = table_grad(L1  ) + R*(table_grad(L1+2)-table_grad(L1  ))
            term_elec = table_grad(L1+1) + R*(table_grad(L1+3)-table_grad(L1+1))
            grad_coef = term_lj12*lj12 - term_temp*lj6 - term_lj6*lj6_ij &
                      + qtmp*charge(iy,j)*term_elec
            work(1) = grad_coef*dij(1)
            work(2) = grad_coef*dij(2)
            work(3) = grad_coef*dij(3)
  
            ! store force
            !
            force_local(1) = force_local(1) - work(1)
            force_local(2) = force_local(2) - work(2)
            force_local(3) = force_local(3) - work(3)
            force(iy,1,j,id+1)  = force(iy,1,j,id+1) + work(1)
            force(iy,2,j,id+1)  = force(iy,2,j,id+1) + work(2)
            force(iy,3,j,id+1)  = force(iy,3,j,id+1) + work(3)

          end if

        end do

        if (check_virial == 1) &
          virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
        force(ix,1,i,id+1) = force(ix,1,i,id+1) + force_local(1)
        force(ix,2,i,id+1) = force(ix,2,i,id+1) + force_local(2)
        force(ix,3,i,id+1) = force(ix,3,i,id+1) + force_local(3)

      end do
    end do
    !$omp end parallel

    return

  end subroutine compute_force_nonbond_tbl_ljpme_intel_knl


end module sp_energy_nonbond_intel_knl_mod
