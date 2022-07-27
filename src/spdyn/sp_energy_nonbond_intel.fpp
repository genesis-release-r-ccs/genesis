!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_nonbond_intel_mod
!> @brief   calculate nonbonded energy with table and with linear interpolation
!! @authors Jaewoon Jung(JJ)
!  
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_nonbond_intel_mod

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
  public  :: compute_energy_nonbond_tbl_lnr_intel
  public  :: compute_force_nonbond_tbl_lnr_intel
  public  :: compute_energy_nonbond_notbl_intel
  public  :: compute_force_nonbond_notbl_intel
  public  :: compute_energy_nonbond_tbl_ljpme_intel 
  public  :: compute_force_nonbond_tbl_ljpme_intel 

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_tbl_lnr_intel
  !> @brief        calculate nonbonded energy with lookup table
  !! @authors      JJ, NT
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[inout] atmcls_pbc : atom class number
  !! @param[inout] coord_pbc  : pbc oriented coordinates
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] eelec      : electrostatic energy of target systems
  !! @param[inout] evdw       : van der Waals energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_tbl_lnr_intel( &
                                                 domain, enefunc, pairlist, &
                                                 atmcls_pbc,                &
                                                 coord_pbc, force, virial,  &
                                                 eelec, evdw)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    integer,                  intent(inout) :: atmcls_pbc(:)
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: lj6, lj12
    real(wp)                  :: R
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: rtmp(1:3), qtmp, jqtmp
    real(wp)                  :: elec_temp, evdw_temp
    real(wp)                  :: force_local(3)
    real(wp)                  :: viri(3)
    integer                   :: i, ix, iy, j, k, ki, kki, ij, L, L1
    integer                   :: id, omp_get_thread_num
    integer                   :: iatmcls, jatmcls
    integer                   :: ncell, ncell_local

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: charge(:,:)
    real(wp),         pointer :: trans1(:,:,:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer,          pointer :: atmcls(:,:)
    integer,          pointer :: natom(:), start_atom(:)
    integer,          pointer :: num_nb15_calc(:,:)
    integer,          pointer :: nb15_calc_list(:,:,:)


    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6

    num_nb15_calc   => pairlist%num_nb15_calc
    nb15_calc_list  => pairlist%nb15_calc_list_fugaku

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff

    !$omp parallel default(shared)                                         &
    !$omp private(id, i, ix, k, ki, kki, iy, ij, j, rij2, L, L1, R, rtmp,  &
    !$omp         qtmp, jqtmp, iatmcls, jatmcls, term_lj12, term_lj6,      &
    !$omp         term_elec, grad_coef, work, force_local, lj6, lj12,      &
    !$omp         elec_temp, evdw_temp, dij, viri)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread
      ki = start_atom(i)
      do ix = 1, natom(i)
        kki = ki + ix
        coord_pbc(kki,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(kki,2,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(kki,3,1) = coord(3,ix,i) + trans1(3,ix,i)
        coord_pbc(kki,4,1) = charge(ix,i)
        atmcls_pbc(kki)    = atmcls(ix,i)
      end do
    end do

    !$omp barrier

    do kki = id+1, domain%num_atom_domain, nthread

      rtmp(1) = coord_pbc(kki,1,1)
      rtmp(2) = coord_pbc(kki,2,1)
      rtmp(3) = coord_pbc(kki,3,1)
      qtmp    = coord_pbc(kki,4,1)
      iatmcls = atmcls_pbc(kki)

      force_local(1) = 0.0_wp
      force_local(2) = 0.0_wp
      force_local(3) = 0.0_wp
      elec_temp      = 0.0_wp
      evdw_temp      = 0.0_wp
      viri(1)        = 0.0_wp
      viri(2)        = 0.0_wp
      viri(3)        = 0.0_wp

!ocl norecurrence
      do k = 1, num_nb15_calc(kki,1)

        ij = nb15_calc_list(k,kki,1)

        dij(1) = rtmp(1) - coord_pbc(ij,1,1)
        dij(2) = rtmp(2) - coord_pbc(ij,2,1)
        dij(3) = rtmp(3) - coord_pbc(ij,3,1)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        rij2    = cutoff2*density / rij2

        jatmcls = atmcls_pbc(ij)
        lj12    = nonb_lj12(jatmcls,iatmcls)
        lj6     = nonb_lj6 (jatmcls,iatmcls)
        jqtmp   = coord_pbc(ij,4,1)

        L  = int(rij2)
        R  = rij2 - L
        L1 = 3*L - 2

        term_lj12 = table_ene(L1)   + R*(table_ene(L1+3)-table_ene(L1))
        term_lj6  = table_ene(L1+1) + R*(table_ene(L1+4)-table_ene(L1+1))
        term_elec = table_ene(L1+2) + R*(table_ene(L1+5)-table_ene(L1+2))
        evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
        elec_temp = elec_temp + qtmp*jqtmp*term_elec
        term_lj12 = table_grad(L1)   + R*(table_grad(L1+3)-table_grad(L1))
        term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
        term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
        grad_coef = term_lj12*lj12 - term_lj6*lj6 + qtmp*jqtmp*term_elec
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(ij,1,1,id+1) = force(ij,1,1,id+1) + work(1)
        force(ij,2,1,id+1) = force(ij,2,1,id+1) + work(2)
        force(ij,3,1,id+1) = force(ij,3,1,id+1) + work(3)

      end do

      force(kki,1,1,id+1) = force(kki,1,1,id+1) + force_local(1)
      force(kki,2,1,id+1) = force(kki,2,1,id+1) + force_local(2)
      force(kki,3,1,id+1) = force(kki,3,1,id+1) + force_local(3)
      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
      eelec(id+1) = eelec(id+1) + elec_temp
      evdw(id+1) = evdw(id+1) + evdw_temp

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond_tbl_lnr_intel

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_tbl_lnr_intel
  !> @brief        calculate nonbonded force without solvents with lookup table
  !! @authors      JJ 
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[in]    npt        : flag for NPT or not
  !! @param[inout] atmcls_pbc : atom class number
  !! @param[inout] coord_pbc  : pbc oriented coordinates
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_tbl_lnr_intel( &
                                                 domain, enefunc, pairlist,   &
                                                 npt, atmcls_pbc,             &
                                                 coord_pbc, force, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    logical,                  intent(in)    :: npt
    integer,                  intent(inout) :: atmcls_pbc(:)
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)


    if (npt) then
      call compute_force_nonbond_tbl_lnr_intel_npt( &  
           domain, enefunc, pairlist, atmcls_pbc,    &
           coord_pbc, force, virial)
    else
      call compute_force_nonbond_tbl_lnr_intel_nonpt( & 
           domain, enefunc, pairlist, atmcls_pbc,      &
           coord_pbc, force)
    end if

    return

  end subroutine compute_force_nonbond_tbl_lnr_intel

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_tbl_lnr_intel_npt
  !> @brief        calculate nonbonded force with virial
  !! @authors      JJ 
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[inout] atmcls_pbc : atom class number
  !! @param[inout] coord_pbc  : pbc oriented coordinates
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_tbl_lnr_intel_npt( &
                                                 domain, enefunc, pairlist, &
                                                 atmcls_pbc,                &
                                                 coord_pbc, force, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    integer,                  intent(inout) :: atmcls_pbc(:)
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)

    ! local variables
    real(wp)                  :: dij1, dij2, dij3, rij2
    real(wp)                  :: R, lj12, lj6
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: rtmp(1:3), qtmp, jqtmp
    real(wp)                  :: force_local(3), viri(3)
    integer                   :: i, ix, ixx, iy, j, k, ki, ij, L, L1
    integer                   :: id, omp_get_thread_num
    integer                   :: iatmcls,jatmcls
    integer                   :: ncell, ncell_local

    real(wip),        pointer,contiguous :: coord(:,:,:)
    real(wp),         pointer,contiguous :: charge(:,:)
    real(wp),         pointer,contiguous :: trans1(:,:,:)
    real(wp),         pointer,contiguous :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer            :: density, cutoff
    real(wp),         pointer,contiguous :: table_grad(:)
    integer,          pointer,contiguous :: atmcls(:,:)
    integer,          pointer,contiguous :: natom(:), start_atom(:)
    integer,          pointer,contiguous :: num_nb15_calc(:,:)
    integer,          pointer,contiguous :: nb15_calc_list(:,:,:)

    real(dp)                  :: sas,eae


    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge

    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6

    num_nb15_calc   => pairlist%num_nb15_calc
    nb15_calc_list  => pairlist%nb15_calc_list_fugaku

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff

    !$omp parallel default(shared)                                          &
    !$omp private(id, i, ix, rtmp, qtmp, k, ki, ij, rij2, L, L1, R,         &
    !$omp         term_lj12, term_lj6, grad_coef, work, term_elec,          &
    !$omp         force_local, lj12, lj6, iatmcls, jatmcls, dij1, dij2,     &
    !$omp         dij3, jqtmp, viri, ixx)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell, nthread
      ki = start_atom(i)
      do ix = 1, natom(i)
        ixx = ki + ix
        coord_pbc(ixx,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(ixx,2,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(ixx,3,1) = coord(3,ix,i) + trans1(3,ix,i)
        coord_pbc(ixx,4,1) = charge(ix,i)
        atmcls_pbc(ixx)    = atmcls(ix,i)
      end do
    end do

    !$omp barrier

    do ixx = id+1, domain%num_atom_domain, nthread

      rtmp(1) = coord_pbc(ixx,1,1)
      rtmp(2) = coord_pbc(ixx,2,1)
      rtmp(3) = coord_pbc(ixx,3,1)
      qtmp    = coord_pbc(ixx,4,1)
      iatmcls = atmcls_pbc(ixx)
      force_local(1) = 0.0_wp
      force_local(2) = 0.0_wp
      force_local(3) = 0.0_wp
      viri(1) = 0.0_wp
      viri(2) = 0.0_wp
      viri(3) = 0.0_wp

      !dir$ simd
      do k = 1, num_nb15_calc(ixx,1)

        ij = nb15_calc_list(k,ixx,1)

        jatmcls = atmcls_pbc(ij)
        lj12    = nonb_lj12(jatmcls,iatmcls)
        lj6     = nonb_lj6 (jatmcls,iatmcls)

        dij1    = rtmp(1) - coord_pbc(ij,1,1)
        dij2    = rtmp(2) - coord_pbc(ij,2,1)
        dij3    = rtmp(3) - coord_pbc(ij,3,1)

        jqtmp   = coord_pbc(ij,4,1)

        rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3
        rij2   = cutoff2 * density / rij2
        L      = int(rij2)
        R      = rij2 - L
        L1     = 3*L - 2

        work(1)  = table_grad(L1)
        work(2)  = table_grad(L1+1)
        work(3)  = table_grad(L1+2)
        work(1)  = work(1) + R*(table_grad(L1+3)-work(1))
        work(2)  = work(2) + R*(table_grad(L1+4)-work(2))
        work(3)  = work(3) + R*(table_grad(L1+5)-work(3))

        grad_coef = work(1)*lj12 - work(2)*lj6 + qtmp*jqtmp*work(3)
        work(1) = grad_coef*dij1
        work(2) = grad_coef*dij2
        work(3) = grad_coef*dij3
        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(ij,1,1,id+1) = force(ij,1,1,id+1) + work(1)
        force(ij,2,1,id+1) = force(ij,2,1,id+1) + work(2)
        force(ij,3,1,id+1) = force(ij,3,1,id+1) + work(3)
        viri(1) = viri(1) + dij1*work(1)
        viri(2) = viri(2) + dij2*work(2)
        viri(3) = viri(3) + dij3*work(3)

      end do

      force(ixx,1,1,id+1) = force(ixx,1,1,id+1) + force_local(1)
      force(ixx,2,1,id+1) = force(ixx,2,1,id+1) + force_local(2)
      force(ixx,3,1,id+1) = force(ixx,3,1,id+1) + force_local(3)
      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)

    end do

    !$omp end parallel

    return

  end subroutine compute_force_nonbond_tbl_lnr_intel_npt

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_tbl_lnr_intel_nonpt
  !> @brief        calculate nonbonded force without solvents with charmm
  !! @authors      JJ 
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[inout] atmcls_pbc : atom class number
  !! @param[inout] coord_pbc  : pbc oriented coordinates
  !! @param[inout] force      : forces for each cell
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_tbl_lnr_intel_nonpt( &
                                                 domain, enefunc, pairlist, &
                                                 atmcls_pbc,                &
                                                 coord_pbc, force)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    integer,                  intent(inout) :: atmcls_pbc(:)
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)

    ! local variables
    real(wp)                  :: dij(3), rij2
    real(wp)                  :: R, lj12, lj6
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(3)
    real(wp)                  :: rtmp(3), qtmp, jqtmp
    real(wp)                  :: force_local(3)
    integer                   :: i, ix, iy, j, k, ki, kki, ij, L, L1
    integer                   :: id, omp_get_thread_num
    integer                   :: iatmcls,jatmcls
    integer                   :: ncell, ncell_local

    real(wip),        pointer,contiguous :: coord(:,:,:)
    real(wp),         pointer,contiguous :: charge(:,:)
    real(wp),         pointer,contiguous :: trans1(:,:,:)
    real(wp),         pointer,contiguous :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer            :: density, cutoff
    real(wp),         pointer,contiguous :: table_grad(:)
    integer,          pointer,contiguous :: atmcls(:,:)
    integer,          pointer,contiguous :: natom(:), start_atom(:)
    integer,          pointer,contiguous :: num_nb15_calc(:,:)
    integer,          pointer,contiguous :: nb15_calc_list(:,:,:)

    real(dp)                  :: sas,eae


    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge

    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6

    num_nb15_calc   => pairlist%num_nb15_calc
    nb15_calc_list  => pairlist%nb15_calc_list_fugaku

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff

    !$omp parallel default(shared)                                      &
    !$omp private(id, i, ix, rtmp, qtmp, k, ki, kki, ij, rij2, L, L1,   &
    !$omp         R, term_lj12, term_lj6, grad_coef, work, term_elec,   &
    !$omp         force_local, lj12, lj6, iatmcls, jatmcls, dij,        &
    !$omp         jqtmp)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell, nthread
      ki = start_atom(i)
      do ix = 1, natom(i)
        kki = ki + ix
        coord_pbc(kki,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(kki,2,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(kki,3,1) = coord(3,ix,i) + trans1(3,ix,i)
        coord_pbc(kki,4,1) = charge(ix,i)
        atmcls_pbc(kki)    = atmcls(ix,i)
      end do
    end do

    !$omp barrier

    do kki = id+1, domain%num_atom_domain, nthread

      rtmp(1) = coord_pbc(kki,1,1)
      rtmp(2) = coord_pbc(kki,2,1)
      rtmp(3) = coord_pbc(kki,3,1)
      qtmp    = coord_pbc(kki,4,1)
      iatmcls = atmcls_pbc(kki)
      force_local(1) = 0.0_wp
      force_local(2) = 0.0_wp
      force_local(3) = 0.0_wp

      !dir$ simd
      do k = 1, num_nb15_calc(kki,1)

        ij = nb15_calc_list(k,kki,1)

        jatmcls = atmcls_pbc(ij)
        lj12    = nonb_lj12(jatmcls,iatmcls)
        lj6     = nonb_lj6 (jatmcls,iatmcls)

        dij(1)  = rtmp(1) - coord_pbc(ij,1,1)
        dij(2)  = rtmp(2) - coord_pbc(ij,2,1)
        dij(3)  = rtmp(3) - coord_pbc(ij,3,1)
        jqtmp   = coord_pbc(ij,4,1)

        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        rij2   = cutoff2 * density / rij2
        L      = int(rij2)
        R      = rij2 - L
        L1     = 3*L - 2

        work(1)  = table_grad(L1)
        work(2)  = table_grad(L1+1)
        work(3)  = table_grad(L1+2)
        work(1)  = work(1) + R*(table_grad(L1+3)-work(1))
        work(2)  = work(2) + R*(table_grad(L1+4)-work(2))
        work(3)  = work(3) + R*(table_grad(L1+5)-work(3))
        grad_coef = work(1)*lj12 - work(2)*lj6 + qtmp*jqtmp*work(3)
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(ij,1,1,id+1) = force(ij,1,1,id+1) + work(1)
        force(ij,2,1,id+1) = force(ij,2,1,id+1) + work(2)
        force(ij,3,1,id+1) = force(ij,3,1,id+1) + work(3)

      end do

      force(kki,1,1,id+1) = force(kki,1,1,id+1) + force_local(1)
      force(kki,2,1,id+1) = force(kki,2,1,id+1) + force_local(2)
      force(kki,3,1,id+1) = force(kki,3,1,id+1) + force_local(3)

    end do

    !$omp end parallel

    return

  end subroutine compute_force_nonbond_tbl_lnr_intel_nonpt

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_notbl_intel
  !> @brief        calculate nonbonded energy with lookup table for elec
  !! @authors      JJ, NT
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[inout] atmcls_pbc : atom class number
  !! @param[inout] coord_pbc  : pbc oriented coordinates
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] eelec      : electrostatic energy of target systems
  !! @param[inout] evdw       : van der Waals energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_notbl_intel( &
                                                 domain, enefunc, pairlist, &
                                                 atmcls_pbc,                &
                                                 coord_pbc, force, virial,  &
                                                 eelec, evdw)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    integer,                  intent(inout) :: atmcls_pbc(:)
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                  :: dij(1:3), rij2, rij2_inv
    real(wp)                  :: lj6, lj12
    real(wp)                  :: R
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: rtmp(1:3), qtmp, jqtmp
    real(wp)                  :: elec_temp, evdw_temp
    real(wp)                  :: force_local(3)
    real(wp)                  :: viri(3), within_cutoff
    integer                   :: i, ix, iy, j, k, ki, kki, ij, L, L1
    integer                   :: id, omp_get_thread_num
    integer                   :: iatmcls, jatmcls
    integer                   :: ncell, ncell_local

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: charge(:,:)
    real(wp),         pointer :: trans1(:,:,:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer,          pointer :: atmcls(:,:)
    integer,          pointer :: natom(:), start_atom(:)
    integer,          pointer :: num_nb15_calc(:,:)
    integer,          pointer :: nb15_calc_list(:,:,:)


    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6

    num_nb15_calc   => pairlist%num_nb15_calc
    nb15_calc_list  => pairlist%nb15_calc_list_fugaku

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff

    !$omp parallel default(shared)                                         &
    !$omp private(id, i, ix, k, ki, kki, iy, ij, j, rij2, L, L1, R, rtmp,  &
    !$omp         qtmp, jqtmp, iatmcls, jatmcls, term_lj12, term_lj6,      &
    !$omp         term_elec, grad_coef, work, force_local, lj6, lj12,      &
    !$omp         elec_temp, evdw_temp, dij, viri, within_cutoff, rij2_inv)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread
      ki = start_atom(i)
      do ix = 1, natom(i)
        kki = ki + ix
        coord_pbc(kki,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(kki,2,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(kki,3,1) = coord(3,ix,i) + trans1(3,ix,i)
        coord_pbc(kki,4,1) = charge(ix,i)
        atmcls_pbc(kki)    = atmcls(ix,i)
      end do
    end do

    !$omp barrier

    do kki = id+1, domain%num_atom_domain, nthread

      rtmp(1) = coord_pbc(kki,1,1)
      rtmp(2) = coord_pbc(kki,2,1)
      rtmp(3) = coord_pbc(kki,3,1)
      qtmp    = coord_pbc(kki,4,1)
      iatmcls = atmcls_pbc(kki)

      force_local(1) = 0.0_wp
      force_local(2) = 0.0_wp
      force_local(3) = 0.0_wp
      elec_temp      = 0.0_wp
      evdw_temp      = 0.0_wp
      viri(1)        = 0.0_wp
      viri(2)        = 0.0_wp
      viri(3)        = 0.0_wp

!ocl norecurrence
      do k = 1, num_nb15_calc(kki,1)

        ij = nb15_calc_list(k,kki,1)

        dij(1) = rtmp(1) - coord_pbc(ij,1,1)
        dij(2) = rtmp(2) - coord_pbc(ij,2,1)
        dij(3) = rtmp(3) - coord_pbc(ij,3,1)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij2    = cutoff2*density*rij2_inv

        jatmcls = atmcls_pbc(ij)
        lj12    = nonb_lj12(jatmcls,iatmcls)
        lj6     = nonb_lj6 (jatmcls,iatmcls)
        jqtmp   = coord_pbc(ij,4,1)

        L  = int(rij2)
        R  = rij2 - L

        term_lj6  = rij2_inv * rij2_inv * rij2_inv * within_cutoff
        term_lj12 = term_lj6 * term_lj6 
        term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
        term_elec = term_elec * within_cutoff
        evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
        elec_temp = elec_temp + qtmp*jqtmp*term_elec
        term_lj12 = -12.0_wp * term_lj12 * rij2_inv
        term_lj6  = - 6.0_wp * term_lj6 * rij2_inv
        term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
        term_elec = term_elec * within_cutoff
        grad_coef = term_lj12*lj12 - term_lj6*lj6 + qtmp*jqtmp*term_elec
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(ij,1,1,id+1) = force(ij,1,1,id+1) + work(1)
        force(ij,2,1,id+1) = force(ij,2,1,id+1) + work(2)
        force(ij,3,1,id+1) = force(ij,3,1,id+1) + work(3)

      end do

      force(kki,1,1,id+1) = force(kki,1,1,id+1) + force_local(1)
      force(kki,2,1,id+1) = force(kki,2,1,id+1) + force_local(2)
      force(kki,3,1,id+1) = force(kki,3,1,id+1) + force_local(3)
      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
      eelec(id+1) = eelec(id+1) + elec_temp
      evdw(id+1) = evdw(id+1) + evdw_temp

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond_notbl_intel

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_notbl_intel
  !> @brief        calculate nonbonded force without solvents with lookup table
  !! @authors      JJ 
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[in]    npt        : flag for NPT or not
  !! @param[inout] atmcls_pbc : atom class number
  !! @param[inout] coord_pbc  : pbc oriented coordinates
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_notbl_intel( &
                                                 domain, enefunc, pairlist,   &
                                                 npt, atmcls_pbc,             &
                                                 coord_pbc, force, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    logical,                  intent(in)    :: npt
    integer,                  intent(inout) :: atmcls_pbc(:)
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)


    if (npt) then
      call compute_force_nonbond_notbl_intel_npt( &  
           domain, enefunc, pairlist, atmcls_pbc,    &
           coord_pbc, force, virial)

    else
      call compute_force_nonbond_notbl_intel_nonpt( & 
           domain, enefunc, pairlist, atmcls_pbc,      &
           coord_pbc, force)
    end if

    return

  end subroutine compute_force_nonbond_notbl_intel

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_notbl_intel_npt
  !> @brief        calculate nonbonded force with virial
  !! @authors      JJ 
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[inout] atmcls_pbc : atom class number
  !! @param[inout] coord_pbc  : pbc oriented coordinates
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_notbl_intel_npt( &
                                                 domain, enefunc, pairlist, &
                                                 atmcls_pbc,                &
                                                 coord_pbc, force, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    integer,                  intent(inout) :: atmcls_pbc(:)
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)

    ! local variables
    real(wp)                  :: dij1, dij2, dij3, rij2, rij2_inv
    real(wp)                  :: R, lj12, lj6
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: rtmp(1:3), qtmp, jqtmp
    real(wp)                  :: force_local(3), viri(3), within_cutoff
    integer                   :: i, ix, iy, j, k, ki, ij, L, L1, ixx
    integer                   :: id, omp_get_thread_num
    integer                   :: iatmcls,jatmcls
    integer                   :: ncell, ncell_local

    real(wip),        pointer,contiguous :: coord(:,:,:)
    real(wp),         pointer,contiguous :: charge(:,:)
    real(wp),         pointer,contiguous :: trans1(:,:,:)
    real(wp),         pointer,contiguous :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer            :: density, cutoff
    real(wp),         pointer,contiguous :: table_grad(:)
    integer,          pointer,contiguous :: atmcls(:,:)
    integer,          pointer,contiguous :: natom(:), start_atom(:)
    integer,          pointer,contiguous :: num_nb15_calc(:,:)
    integer,          pointer,contiguous :: nb15_calc_list(:,:,:)

    real(dp)                  :: sas,eae


    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge

    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6

    num_nb15_calc   => pairlist%num_nb15_calc
    nb15_calc_list  => pairlist%nb15_calc_list_fugaku

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff

    !$omp parallel default(shared)                                          &
    !$omp private(id, i, ix, rtmp, qtmp, k, ki, ij, rij2, L, L1, R,         &
    !$omp         term_lj12, term_lj6, grad_coef, work, term_elec,          &
    !$omp         force_local, lj12, lj6, iatmcls, jatmcls, dij1, dij2,     &
    !$omp         dij3, jqtmp, viri, within_cutoff, rij2_inv, ixx)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell, nthread
      ki = start_atom(i)
      do ix = 1, natom(i)
        ixx = ki + ix
        coord_pbc(ixx,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(ixx,2,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(ixx,3,1) = coord(3,ix,i) + trans1(3,ix,i)
        coord_pbc(ixx,4,1) = charge(ix,i)
        atmcls_pbc(ixx)    = atmcls(ix,i)
      end do
    end do

    !$omp barrier

    do ixx = id+1, domain%num_atom_domain, nthread

      rtmp(1) = coord_pbc(ixx,1,1)
      rtmp(2) = coord_pbc(ixx,2,1)
      rtmp(3) = coord_pbc(ixx,3,1)
      qtmp    = coord_pbc(ixx,4,1)
      iatmcls = atmcls_pbc(ixx)
      force_local(1) = 0.0_wp
      force_local(2) = 0.0_wp
      force_local(3) = 0.0_wp
      viri(1) = 0.0_wp
      viri(2) = 0.0_wp
      viri(3) = 0.0_wp

      do k = 1, num_nb15_calc(ixx,1)

        ij = nb15_calc_list(k,ixx,1)

        jatmcls = atmcls_pbc(ij)
        lj12    = nonb_lj12(jatmcls,iatmcls)
        lj6     = nonb_lj6 (jatmcls,iatmcls)

        dij1    = rtmp(1) - coord_pbc(ij,1,1)
        dij2    = rtmp(2) - coord_pbc(ij,2,1)
        dij3    = rtmp(3) - coord_pbc(ij,3,1)

        jqtmp   = coord_pbc(ij,4,1)

        rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij2   = cutoff2 * density * rij2_inv
        L      = int(rij2)
        R      = rij2 - L

        term_lj6  = rij2_inv * rij2_inv * rij2_inv
        term_lj12 = term_lj6 * term_lj6
        work(1)  = -12.0_wp * term_lj12 * rij2_inv * within_cutoff
        work(2)  = - 6.0_wp * term_lj6  * rij2_inv * within_cutoff
        work(3)  = table_grad(L)
        work(3)  = work(3) + R*(table_grad(L+1)-work(3))
        work(3)  = work(3) * within_cutoff

        grad_coef = work(1)*lj12 - work(2)*lj6 + qtmp*jqtmp*work(3)
        work(1) = grad_coef*dij1
        work(2) = grad_coef*dij2
        work(3) = grad_coef*dij3
        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(ij,1,1,id+1) = force(ij,1,1,id+1) + work(1)
        force(ij,2,1,id+1) = force(ij,2,1,id+1) + work(2)
        force(ij,3,1,id+1) = force(ij,3,1,id+1) + work(3)
        viri(1) = viri(1) + dij1*work(1)
        viri(2) = viri(2) + dij2*work(2)
        viri(3) = viri(3) + dij3*work(3)

      end do

      force(ixx,1,1,id+1) = force(ixx,1,1,id+1) + force_local(1)
      force(ixx,2,1,id+1) = force(ixx,2,1,id+1) + force_local(2)
      force(ixx,3,1,id+1) = force(ixx,3,1,id+1) + force_local(3)
      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)

    end do

    !$omp end parallel

    return

  end subroutine compute_force_nonbond_notbl_intel_npt

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_notbl_intel_nonpt
  !> @brief        calculate nonbonded force without solvents with charmm
  !! @authors      JJ 
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[inout] atmcls_pbc : atom class number
  !! @param[inout] coord_pbc  : pbc oriented coordinates
  !! @param[inout] force      : forces for each cell
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_notbl_intel_nonpt( &
                                                 domain, enefunc, pairlist, &
                                                 atmcls_pbc,                &
                                                 coord_pbc, force)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    integer,                  intent(inout) :: atmcls_pbc(:)
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)

    ! local variables
    real(wp)                  :: dij(3), rij2, rij2_inv
    real(wp)                  :: R, lj12, lj6
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(3)
    real(wp)                  :: rtmp(3), qtmp, jqtmp
    real(wp)                  :: force_local(3), within_cutoff
    integer                   :: i, ix, iy, j, k, ki, kki, ij, L, L1
    integer                   :: id, omp_get_thread_num
    integer                   :: iatmcls,jatmcls
    integer                   :: ncell, ncell_local

    real(wip),        pointer,contiguous :: coord(:,:,:)
    real(wp),         pointer,contiguous :: charge(:,:)
    real(wp),         pointer,contiguous :: trans1(:,:,:)
    real(wp),         pointer,contiguous :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer            :: density, cutoff
    real(wp),         pointer,contiguous :: table_grad(:)
    integer,          pointer,contiguous :: atmcls(:,:)
    integer,          pointer,contiguous :: natom(:), start_atom(:)
    integer,          pointer,contiguous :: num_nb15_calc(:,:)
    integer,          pointer,contiguous :: nb15_calc_list(:,:,:)

    real(dp)                  :: sas,eae


    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge

    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6

    num_nb15_calc   => pairlist%num_nb15_calc
    nb15_calc_list  => pairlist%nb15_calc_list_fugaku

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff

    !$omp parallel default(shared)                                      &
    !$omp private(id, i, ix, rtmp, qtmp, k, ki, kki, ij, rij2, L, L1,   &
    !$omp         R, term_lj12, term_lj6, grad_coef, work, term_elec,   &
    !$omp         force_local, lj12, lj6, iatmcls, jatmcls, dij,        &
    !$omp         jqtmp, rij2_inv, within_cutoff)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell, nthread
      ki = start_atom(i)
      do ix = 1, natom(i)
        kki = ki + ix
        coord_pbc(kki,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(kki,2,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(kki,3,1) = coord(3,ix,i) + trans1(3,ix,i)
        coord_pbc(kki,4,1) = charge(ix,i)
        atmcls_pbc(kki)    = atmcls(ix,i)
      end do
    end do

    !$omp barrier

    do kki = id+1, domain%num_atom_domain, nthread

      rtmp(1) = coord_pbc(kki,1,1)
      rtmp(2) = coord_pbc(kki,2,1)
      rtmp(3) = coord_pbc(kki,3,1)
      qtmp    = coord_pbc(kki,4,1)
      iatmcls = atmcls_pbc(kki)

      force_local(1) = 0.0_wp
      force_local(2) = 0.0_wp
      force_local(3) = 0.0_wp

      do k = 1, num_nb15_calc(kki,1)

        ij = nb15_calc_list(k,kki,1)

        jatmcls = atmcls_pbc(ij)
        lj12    = nonb_lj12(jatmcls,iatmcls)
        lj6     = nonb_lj6 (jatmcls,iatmcls)

        dij(1)  = rtmp(1) - coord_pbc(ij,1,1)
        dij(2)  = rtmp(2) - coord_pbc(ij,2,1)
        dij(3)  = rtmp(3) - coord_pbc(ij,3,1)
        jqtmp   = coord_pbc(ij,4,1)

        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij2   = cutoff2 * density * rij2_inv
        L      = int(rij2)
        R      = rij2 - L

        term_lj6  = rij2_inv * rij2_inv * rij2_inv * within_cutoff
        term_lj12 = term_lj6 * term_lj6
        work(1)  = -12.0_wp * term_lj12 * rij2_inv 
        work(2)  = - 6.0_wp * term_lj6  * rij2_inv
        work(3)  = table_grad(L)
        work(3)  = work(3) + R*(table_grad(L+1)-work(3))
        work(3)  = work(3) * within_cutoff

        grad_coef = work(1)*lj12 - work(2)*lj6 + qtmp*jqtmp*work(3)
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(ij,1,1,id+1) = force(ij,1,1,id+1) + work(1)
        force(ij,2,1,id+1) = force(ij,2,1,id+1) + work(2)
        force(ij,3,1,id+1) = force(ij,3,1,id+1) + work(3)

      end do

      force(kki,1,1,id+1) = force(kki,1,1,id+1) + force_local(1)
      force(kki,2,1,id+1) = force(kki,2,1,id+1) + force_local(2)
      force(kki,3,1,id+1) = force(kki,3,1,id+1) + force_local(3)

    end do

    !$omp end parallel

    return

  end subroutine compute_force_nonbond_notbl_intel_nonpt

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_tbl_ljpme_intel 
  !> @brief        calculate nonbonded energy with lookup table
  !! @authors      JJ, NT
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[inout] atmcls_pbc : atom class number
  !! @param[inout] coord_pbc  : pbc oriented coordinates
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] eelec      : electrostatic energy of target systems
  !! @param[inout] evdw       : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_tbl_ljpme_intel ( &
                                      domain, enefunc, pairlist, atmcls_pbc, &
                                      coord_pbc, force, virial, eelec, evdw)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    integer,                  intent(inout) :: atmcls_pbc(:)
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                  :: dij(1:3), rij2, inv_r2, inv_r6, inv_r12
    real(wp)                  :: inv_cutoff2, inv_cutoff6
    real(wp)                  :: lj6, lj12, lj6_i, lj6_j, lj6_ij
    real(wp)                  :: R
    real(wp)                  :: term_lj12, term_lj6, term_elec, term_temp
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: rtmp(1:3), qtmp, jqtmp
    real(wp)                  :: elec_temp, evdw_temp
    real(wp)                  :: force_local(3)
    real(wp)                  :: viri(3)
    integer                   :: i, ix, iy, j, k, ki, kki, ij, L, L1
    integer                   :: id, omp_get_thread_num
    integer                   :: iatmcls, jatmcls
    integer                   :: ncell, ncell_local

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: charge(:,:)
    real(wp),         pointer :: trans1(:,:,:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: nonb_lj6_factor(:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer,          pointer :: atmcls(:,:)
    integer,          pointer :: natom(:), start_atom(:)
    integer,          pointer :: num_nb15_calc(:,:)
    integer,          pointer :: nb15_calc_list(:,:,:)


    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6
    nonb_lj6_factor => enefunc%nonb_lj6_factor

    num_nb15_calc   => pairlist%num_nb15_calc
    nb15_calc_list  => pairlist%nb15_calc_list_fugaku

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff
    inv_cutoff2     =  1.0_wp /cutoff2
    inv_cutoff6     =  inv_cutoff2 * inv_cutoff2 * inv_cutoff2

    !$omp parallel default(shared)                                         &
    !$omp private(id, i, ix, k, ki, kki, iy, ij, j, rij2, L, L1, R, rtmp,  &
    !$omp         qtmp, jqtmp, iatmcls, jatmcls, term_lj12, term_lj6,      &
    !$omp         term_elec, grad_coef, work, force_local, lj6, lj12,      &
    !$omp         lj6_i, lj6_j, lj6_ij, inv_r2, inv_r6, inv_r12,           &
    !$omp         elec_temp, evdw_temp, dij, viri)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread
      ki = start_atom(i)
      do ix = 1, natom(i)
        kki = ki + ix
        coord_pbc(kki,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(kki,2,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(kki,3,1) = coord(3,ix,i) + trans1(3,ix,i)
        coord_pbc(kki,4,1) = charge(ix,i)
        atmcls_pbc(kki)    = atmcls(ix,i)
      end do
    end do

    !$omp barrier

    do kki = id+1, domain%num_atom_domain, nthread

      rtmp(1) = coord_pbc(kki,1,1)
      rtmp(2) = coord_pbc(kki,2,1)
      rtmp(3) = coord_pbc(kki,3,1)
      qtmp    = coord_pbc(kki,4,1)
      iatmcls = atmcls_pbc(kki)
      lj6_i   = nonb_lj6_factor(iatmcls)

      force_local(1) = 0.0_wp
      force_local(2) = 0.0_wp
      force_local(3) = 0.0_wp
      elec_temp      = 0.0_wp
      evdw_temp      = 0.0_wp
      viri(1)        = 0.0_wp
      viri(2)        = 0.0_wp
      viri(3)        = 0.0_wp

!ocl norecurrence
      do k = 1, num_nb15_calc(kki,1)

        ij = nb15_calc_list(k,kki,1)

        dij(1)  = rtmp(1) - coord_pbc(ij,1,1)
        dij(2)  = rtmp(2) - coord_pbc(ij,2,1)
        dij(3)  = rtmp(3) - coord_pbc(ij,3,1)
        rij2    = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        inv_r2  = 1.0_wp / rij2
        inv_r6  = inv_r2 * inv_r2 * inv_r2
        inv_r12 = inv_r6 * inv_r6
        rij2    = cutoff2 * density * inv_r2

        jatmcls = atmcls_pbc(ij)
        lj6_j   = nonb_lj6_factor(jatmcls)
        lj6_ij  = lj6_i * lj6_j
        lj12    = nonb_lj12(jatmcls,iatmcls)
        lj6     = nonb_lj6 (jatmcls,iatmcls)
        lj6     = lj6 - lj6_ij
        jqtmp   = coord_pbc(ij,4,1)

        L  = int(rij2)
        R  = rij2 - L
        L1 = 2*L - 1

        term_temp = inv_r6 - inv_cutoff6
        term_lj6  = table_ene(L1  ) + R*(table_ene(L1+2)-table_ene(L1  ))
        term_elec = table_ene(L1+1) + R*(table_ene(L1+3)-table_ene(L1+1))
        evdw_temp = evdw_temp + inv_r12*lj12 - term_temp*lj6 - term_lj6*lj6_ij
        elec_temp = elec_temp + qtmp*jqtmp*term_elec

        term_lj12 = -12.0_wp*inv_r12*inv_r2
        term_temp = -6.0_wp*inv_r6*inv_r2
        term_lj6  = table_grad(L1  ) + R*(table_grad(L1+2)-table_grad(L1  ))
        term_elec = table_grad(L1+1) + R*(table_grad(L1+3)-table_grad(L1+1))
        grad_coef = term_lj12*lj12 - term_temp*lj6 - term_lj6*lj6_ij &
                  + qtmp*jqtmp*term_elec

        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(ij,1,1,id+1) = force(ij,1,1,id+1) + work(1)
        force(ij,2,1,id+1) = force(ij,2,1,id+1) + work(2)
        force(ij,3,1,id+1) = force(ij,3,1,id+1) + work(3)

      end do

      force(kki,1,1,id+1) = force(kki,1,1,id+1) + force_local(1)
      force(kki,2,1,id+1) = force(kki,2,1,id+1) + force_local(2)
      force(kki,3,1,id+1) = force(kki,3,1,id+1) + force_local(3)
      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
      eelec(id+1) = eelec(id+1) + elec_temp
      evdw(id+1) = evdw(id+1) + evdw_temp

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond_tbl_ljpme_intel

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_tbl_ljpme_intel
  !> @brief        calculate nonbonded force without solvents with lookup table
  !! @authors      JJ
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[in]    npt        : flag for NPT or not
  !! @param[inout] atmcls_pbc : atom class number
  !! @param[inout] coord_pbc  : pbc oriented coordinates
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_tbl_ljpme_intel( &
                                                 domain, enefunc, pairlist,   &
                                                 npt, atmcls_pbc,             &
                                                 coord_pbc, force, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    logical,                  intent(in)    :: npt
    integer,                  intent(inout) :: atmcls_pbc(:)
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)

    if (npt) then

      call compute_force_nonbond_tbl_ljpme_intel_npt( &
                              domain, enefunc, pairlist, atmcls_pbc, &
                              coord_pbc, force, virial)

    else

      call compute_force_nonbond_tbl_ljpme_intel_nonpt( &
                              domain, enefunc, pairlist, atmcls_pbc, &
                              coord_pbc, force)

    end if

    return

  end subroutine compute_force_nonbond_tbl_ljpme_intel

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_tbl_ljpme_intel_npt
  !> @brief        calculate nonbonded force with virial
  !! @authors      JJ
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[inout] atmcls_pbc : atom class number
  !! @param[inout] coord_pbc  : pbc oriented coordinates
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_tbl_ljpme_intel_npt( &
                                                 domain, enefunc, pairlist, &
                                                 atmcls_pbc,                &
                                                 coord_pbc, force, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    integer,                  intent(inout) :: atmcls_pbc(:)
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)

    ! local variables
    real(wp)                  :: dij1, dij2, dij3, rij2, inv_r2, inv_r6, inv_r12
    real(wp)                  :: R, lj12, lj6, lj6_i, lj6_j, lj6_ij
    real(wp)                  :: term_lj12, term_lj6, term_elec, term_temp
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: rtmp(1:3), qtmp, jqtmp
    real(wp)                  :: force_local(3), viri(3)
    integer                   :: i, ix, ixx, iy, j, k, ki, ij, L, L1
    integer                   :: id, omp_get_thread_num
    integer                   :: iatmcls,jatmcls
    integer                   :: ncell, ncell_local

    real(wip),        pointer,contiguous :: coord(:,:,:)
    real(wp),         pointer,contiguous :: charge(:,:)
    real(wp),         pointer,contiguous :: trans1(:,:,:)
    real(wp),         pointer,contiguous :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer,contiguous :: nonb_lj6_factor(:)
    real(wp),         pointer            :: density, cutoff
    real(wp),         pointer,contiguous :: table_grad(:)
    integer,          pointer,contiguous :: atmcls(:,:)
    integer,          pointer,contiguous :: natom(:), start_atom(:)
    integer,          pointer,contiguous :: num_nb15_calc(:,:)
    integer,          pointer,contiguous :: nb15_calc_list(:,:,:)


    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge

    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6
    nonb_lj6_factor => enefunc%nonb_lj6_factor

    num_nb15_calc   => pairlist%num_nb15_calc
    nb15_calc_list  => pairlist%nb15_calc_list_fugaku

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff


    !$omp parallel default(shared)                                          &
    !$omp private(id, i, ix, rtmp, qtmp, k, ki, ij, rij2, L, L1, R,         &
    !$omp         term_lj12, term_lj6, grad_coef, work, term_elec,          &
    !$omp         force_local, lj12, lj6, lj6_i, lj6_j, lj6_ij, inv_r2,     &
    !$omp         inv_r6, inv_r12, iatmcls, jatmcls, dij1, dij2, dij3,      &
    !$omp         jqtmp, viri, ixx)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell, nthread
      ki = start_atom(i)
      do ix = 1, natom(i)
        ixx = ki + ix
        coord_pbc(ixx,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(ixx,2,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(ixx,3,1) = coord(3,ix,i) + trans1(3,ix,i)
        coord_pbc(ixx,4,1) = charge(ix,i)
        atmcls_pbc(ixx)    = atmcls(ix,i)
      end do
    end do

    !$omp barrier

    do ixx = id+1, domain%num_atom_domain, nthread

      rtmp(1) = coord_pbc(ixx,1,1)
      rtmp(2) = coord_pbc(ixx,2,1)
      rtmp(3) = coord_pbc(ixx,3,1)
      qtmp    = coord_pbc(ixx,4,1)
      iatmcls = atmcls_pbc(ixx)
      lj6_i   = nonb_lj6_factor(iatmcls)
      force_local(1) = 0.0_wp
      force_local(2) = 0.0_wp
      force_local(3) = 0.0_wp
      viri(1) = 0.0_wp
      viri(2) = 0.0_wp
      viri(3) = 0.0_wp

!ocl norecurrence

      do k = 1, num_nb15_calc(ixx,1)

        ij = nb15_calc_list(k,ixx,1)

        jatmcls = atmcls_pbc(ij)
        lj6_j   = nonb_lj6_factor(jatmcls)
        lj6_ij  = lj6_i * lj6_j
        lj12    = nonb_lj12(jatmcls,iatmcls)
        lj6     = nonb_lj6 (jatmcls,iatmcls)
        lj6     = lj6 - lj6_ij

        dij1    = rtmp(1) - coord_pbc(ij,1,1)
        dij2    = rtmp(2) - coord_pbc(ij,2,1)
        dij3    = rtmp(3) - coord_pbc(ij,3,1)

        jqtmp   = coord_pbc(ij,4,1)

        rij2    = dij1*dij1 + dij2*dij2 + dij3*dij3
        inv_r2  = 1.0_wp / rij2
        inv_r6  = inv_r2 * inv_r2 * inv_r2
        inv_r12 = inv_r6 * inv_r6
        rij2    = cutoff2 * density * inv_r2
        L       = int(rij2)
        R       = rij2 - L
        L1      = 2*L - 1

        term_lj12 = -12.0_wp*inv_r12*inv_r2
        term_temp = -6.0_wp*inv_r6*inv_r2
        work(1)  = table_grad(L1)
        work(2)  = table_grad(L1+1)
        work(1)  = work(1) + R*(table_grad(L1+2)-work(1))
        work(2)  = work(2) + R*(table_grad(L1+3)-work(2))

        grad_coef = term_lj12*lj12 - term_temp*lj6 - work(1)*lj6_ij &
                  + qtmp*jqtmp*work(2)
        work(1) = grad_coef*dij1
        work(2) = grad_coef*dij2
        work(3) = grad_coef*dij3
        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(ij,1,1,id+1) = force(ij,1,1,id+1) + work(1)
        force(ij,2,1,id+1) = force(ij,2,1,id+1) + work(2)
        force(ij,3,1,id+1) = force(ij,3,1,id+1) + work(3)
        viri(1) = viri(1) + dij1*work(1)
        viri(2) = viri(2) + dij2*work(2)
        viri(3) = viri(3) + dij3*work(3)

      end do

      force(ixx,1,1,id+1) = force(ixx,1,1,id+1) + force_local(1)
      force(ixx,2,1,id+1) = force(ixx,2,1,id+1) + force_local(2)
      force(ixx,3,1,id+1) = force(ixx,3,1,id+1) + force_local(3)
      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)

    end do

    !$omp end parallel

    return

  end subroutine compute_force_nonbond_tbl_ljpme_intel_npt

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_tbl_ljpme_intel_nonpt
  !> @brief        calculate nonbonded force with virial
  !! @authors      JJ
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[inout] atmcls_pbc : atom class number
  !! @param[inout] coord_pbc  : pbc oriented coordinates
  !! @param[inout] force      : forces for each cell
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_tbl_ljpme_intel_nonpt( &
                                                 domain, enefunc, pairlist, &
                                                 atmcls_pbc,                &
                                                 coord_pbc, force)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    integer,                  intent(inout) :: atmcls_pbc(:)
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)

    ! local variables
    real(wp)                  :: dij1, dij2, dij3, rij2, inv_r2, inv_r6, inv_r12
    real(wp)                  :: R, lj12, lj6, lj6_i, lj6_j, lj6_ij
    real(wp)                  :: term_lj12, term_lj6, term_elec, term_temp
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: rtmp(1:3), qtmp, jqtmp
    real(wp)                  :: force_local(3)
    integer                   :: i, ix, ixx, iy, j, k, ki, ij, L, L1
    integer                   :: id, omp_get_thread_num
    integer                   :: iatmcls,jatmcls
    integer                   :: ncell, ncell_local

    real(wip),        pointer,contiguous :: coord(:,:,:)
    real(wp),         pointer,contiguous :: charge(:,:)
    real(wp),         pointer,contiguous :: trans1(:,:,:)
    real(wp),         pointer,contiguous :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer,contiguous :: nonb_lj6_factor(:)
    real(wp),         pointer            :: density, cutoff
    real(wp),         pointer,contiguous :: table_grad(:)
    integer,          pointer,contiguous :: atmcls(:,:)
    integer,          pointer,contiguous :: natom(:), start_atom(:)
    integer,          pointer,contiguous :: num_nb15_calc(:,:)
    integer,          pointer,contiguous :: nb15_calc_list(:,:,:)


    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge

    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6
    nonb_lj6_factor => enefunc%nonb_lj6_factor

    num_nb15_calc   => pairlist%num_nb15_calc
    nb15_calc_list  => pairlist%nb15_calc_list_fugaku

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff


    !$omp parallel default(shared)                                          &
    !$omp private(id, i, ix, rtmp, qtmp, k, ki, ij, rij2, L, L1, R,         &
    !$omp         term_lj12, term_lj6, grad_coef, work, term_elec,          &
    !$omp         force_local, lj12, lj6, lj6_i, lj6_j, lj6_ij, inv_r2,     &
    !$omp         inv_r6, inv_r12, iatmcls, jatmcls, dij1, dij2, dij3,      &
    !$omp         jqtmp, ixx)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell, nthread
      ki = start_atom(i)
      do ix = 1, natom(i)
        ixx = ki + ix
        coord_pbc(ixx,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(ixx,2,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(ixx,3,1) = coord(3,ix,i) + trans1(3,ix,i)
        coord_pbc(ixx,4,1) = charge(ix,i)
        atmcls_pbc(ixx)    = atmcls(ix,i)
      end do
    end do

    !$omp barrier

    do ixx = id+1, domain%num_atom_domain, nthread

      rtmp(1) = coord_pbc(ixx,1,1)
      rtmp(2) = coord_pbc(ixx,2,1)
      rtmp(3) = coord_pbc(ixx,3,1)
      qtmp    = coord_pbc(ixx,4,1)
      iatmcls = atmcls_pbc(ixx)
      lj6_i   = nonb_lj6_factor(iatmcls)
      force_local(1) = 0.0_wp
      force_local(2) = 0.0_wp
      force_local(3) = 0.0_wp

!ocl norecurrence

      do k = 1, num_nb15_calc(ixx,1)

        ij = nb15_calc_list(k,ixx,1)

        jatmcls = atmcls_pbc(ij)
        lj6_j   = nonb_lj6_factor(jatmcls)
        lj6_ij  = lj6_i * lj6_j
        lj12    = nonb_lj12(jatmcls,iatmcls)
        lj6     = nonb_lj6 (jatmcls,iatmcls)
        lj6     = lj6 - lj6_ij

        dij1    = rtmp(1) - coord_pbc(ij,1,1)
        dij2    = rtmp(2) - coord_pbc(ij,2,1)
        dij3    = rtmp(3) - coord_pbc(ij,3,1)

        jqtmp   = coord_pbc(ij,4,1)

        rij2    = dij1*dij1 + dij2*dij2 + dij3*dij3
        inv_r2  = 1.0_wp / rij2
        inv_r6  = inv_r2 * inv_r2 * inv_r2
        inv_r12 = inv_r6 * inv_r6
        rij2    = cutoff2 * density * inv_r2
        L       = int(rij2)
        R       = rij2 - L
        L1      = 2*L - 1

        term_lj12 = -12.0_wp*inv_r12*inv_r2
        term_temp = -6.0_wp*inv_r6*inv_r2
        work(1)  = table_grad(L1)
        work(2)  = table_grad(L1+1)
        work(1)  = work(1) + R*(table_grad(L1+2)-work(1))
        work(2)  = work(2) + R*(table_grad(L1+3)-work(2))

        grad_coef = term_lj12*lj12 - term_temp*lj6 - work(1)*lj6_ij &
                  + qtmp*jqtmp*work(2)
        work(1) = grad_coef*dij1
        work(2) = grad_coef*dij2
        work(3) = grad_coef*dij3
        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(ij,1,1,id+1) = force(ij,1,1,id+1) + work(1)
        force(ij,2,1,id+1) = force(ij,2,1,id+1) + work(2)
        force(ij,3,1,id+1) = force(ij,3,1,id+1) + work(3)

      end do

      force(ixx,1,1,id+1) = force(ixx,1,1,id+1) + force_local(1)
      force(ixx,2,1,id+1) = force(ixx,2,1,id+1) + force_local(2)
      force(ixx,3,1,id+1) = force(ixx,3,1,id+1) + force_local(3)

    end do

    !$omp end parallel

    return

  end subroutine compute_force_nonbond_tbl_ljpme_intel_nonpt

end module sp_energy_nonbond_intel_mod
