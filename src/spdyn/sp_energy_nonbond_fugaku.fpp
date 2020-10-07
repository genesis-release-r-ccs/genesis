!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_nonbond_fugaku_mod
!> @brief   calculate nonbonded energy with table and with linear interpolation
!! @authors Jaewoon Jung(JJ)
!  
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_nonbond_fugaku_mod

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
  public  :: compute_energy_nonbond_tbl_lnr_fugaku
  public  :: compute_force_nonbond_tbl_lnr_fugaku
  public  :: compute_energy_nonbond_notbl_fugaku
  public  :: compute_force_nonbond_notbl_fugaku

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_tbl_lnr_fugaku
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

  subroutine compute_energy_nonbond_tbl_lnr_fugaku( &
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
    integer,          pointer :: natom(:)
    integer,          pointer :: num_nb15_calc(:,:)
    integer,          pointer :: nb15_calc_list(:,:,:)

    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
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

#ifdef PKTIMER
    call timer_sta(222)
#ifdef FJ_PROF_FAPP
    call fapp_start("compute_energy_nonbond_table_linear",222,0)
#endif
#endif

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
      ki = (i-1)*MaxAtom
      do ix = 1, natom(i)
        kki = ki + ix
        coord_pbc(1,kki,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(2,kki,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(3,kki,1) = coord(3,ix,i) + trans1(3,ix,i)
        coord_pbc(4,kki,1) = charge(ix,i)
        atmcls_pbc(kki)    = atmcls(ix,i)
      end do
    end do

    !$omp barrier

    do i = id+1, ncell, nthread

      ki = (i-1)*MaxAtom

      do ix = 1, natom(i)

        kki = ki + ix
        rtmp(1) = coord_pbc(1,kki,1)
        rtmp(2) = coord_pbc(2,kki,1)
        rtmp(3) = coord_pbc(3,kki,1)
        qtmp    = coord_pbc(4,kki,1)
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
        do k = 1, num_nb15_calc(ix,i)

          ij = nb15_calc_list(k,ix,i)

          dij(1) = rtmp(1) - coord_pbc(1,ij,1)
          dij(2) = rtmp(2) - coord_pbc(2,ij,1)
          dij(3) = rtmp(3) - coord_pbc(3,ij,1)
          rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2    = cutoff2*density / rij2

          jatmcls = atmcls_pbc(ij)
          lj12    = nonb_lj12(jatmcls,iatmcls)
          lj6     = nonb_lj6 (jatmcls,iatmcls)
          jqtmp   = coord_pbc(4,ij,1)

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
          force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
          force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
          force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

        end do

        force(1,ki+ix,1,id+1) = force(1,ki+ix,1,id+1) + force_local(1)
        force(2,ki+ix,1,id+1) = force(2,ki+ix,1,id+1) + force_local(2)
        force(3,ki+ix,1,id+1) = force(3,ki+ix,1,id+1) + force_local(3)
        virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
        virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
        virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
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

  end subroutine compute_energy_nonbond_tbl_lnr_fugaku

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_tbl_lnr_fugaku
  !> @brief        calculate nonbonded force without solvents with lookup table
  !! @authors      JJ 
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] virial   : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_tbl_lnr_fugaku( &
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
#ifdef FJ_TIMER_2
      call timer_sta(411)
#endif
      call compute_force_nonbond_tbl_lnr_fugaku_npt( &  
           domain, enefunc, pairlist, atmcls_pbc,    &
           coord_pbc, force, virial)
#ifdef FJ_TIMER_2
      call timer_end(411)
#endif

    else
#ifdef FJ_TIMER_2
      call timer_sta(412)
#endif
      call compute_force_nonbond_tbl_lnr_fugaku_nonpt( & 
           domain, enefunc, pairlist, atmcls_pbc,      &
           coord_pbc, force)
#ifdef FJ_TIMER_2
      call timer_end(412)
#endif
    end if

    return

  end subroutine compute_force_nonbond_tbl_lnr_fugaku

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_tbl_lnr_fugaku_npt
  !> @brief        calculate nonbonded force with virial
  !! @authors      JJ 
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] virial   : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_tbl_lnr_fugaku_npt( &
                                                 domain, enefunc, pairlist, &
                                                 atmcls_pbc,                &
                                                 coord_pbc, force, virial)

#ifdef PKTIMER
    use Ctim
#endif
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
    integer                   :: i, ix, iy, j, k, ki, ij, L, L1
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
    integer,          pointer,contiguous :: natom(:)
    integer,          pointer,contiguous :: num_nb15_calc(:,:)
    integer,          pointer,contiguous :: nb15_calc_list(:,:,:)

    real(dp)                  :: sas,eae

    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
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

#ifdef PKTIMER
    call gettod(sas)
#ifdef FJ_PROF_FAPP
    call fapp_start("Nonb15F",11,0)
#endif
    call timer_sta(11)
#endif

    !$omp parallel default(shared)                                          &
    !$omp private(id, i, ix, rtmp, qtmp, k, ki, ij, rij2, L, L1, R,         &
    !$omp         term_lj12, term_lj6, grad_coef, work, term_elec,          &
    !$omp         force_local, lj12, lj6, iatmcls, jatmcls, dij1, dij2,     &
    !$omp         dij3, jqtmp, viri)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell, nthread
      ki = (i-1)*MaxAtom
      do ix = 1, natom(i)
        coord_pbc(1,ki+ix,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(2,ki+ix,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(3,ki+ix,1) = coord(3,ix,i) + trans1(3,ix,i)
        coord_pbc(4,ki+ix,1) = charge(ix,i)
        atmcls_pbc(ki+ix)    = atmcls(ix,i)
      end do
    end do

    !$omp barrier

    do i = id+1, ncell, nthread

      ki = (i-1)*MaxAtom

      do ix = 1, natom(i)

        rtmp(1) = coord_pbc(1,ki+ix,1)
        rtmp(2) = coord_pbc(2,ki+ix,1)
        rtmp(3) = coord_pbc(3,ki+ix,1)
        qtmp    = coord_pbc(4,ki+ix,1)
        iatmcls = atmcls_pbc(ki+ix)
        force_local(1) = 0.0_wp
        force_local(2) = 0.0_wp
        force_local(3) = 0.0_wp
        viri(1) = 0.0_wp
        viri(2) = 0.0_wp
        viri(3) = 0.0_wp

!ocl norecurrence

        do k = 1, num_nb15_calc(ix,i)

          ! prefetch for this ix loop prefetch
          !
          if( mod(k,64)==1 .and. (k+128)<=num_nb15_calc(ix,i) ) then
!ocl prefetch_read(nb15_calc_list(k+128,  ix,i),level=1,strong=1)
          endif

          ! prefetch for next ix loop prefetch
          if( ix < natom(i) ) then
            if((num_nb15_calc(ix,i)<=64.and.k==1) .or. &
               k==(num_nb15_calc(ix,i)-64) ) then
!ocl prefetch_read(nb15_calc_list(1,  ix+1,i),level=1,strong=1)
!ocl prefetch_read(nb15_calc_list(65, ix+1,i),level=1,strong=1)
            endif
          endif

          ij = nb15_calc_list(k+64,ix,i)

!ocl prefetch_read(coord_pbc(1,ij,1)  ,level=1,strong=1)
!ocl prefetch_read(atmcls_pbc(ij)     ,level=1,strong=1)
!ocl prefetch_write(force(1,ij,1,id+1),level=1,strong=1)

          ij = nb15_calc_list(k,ix,i)

          jatmcls = atmcls_pbc(ij)
          lj12    = nonb_lj12(jatmcls,iatmcls)
          lj6     = nonb_lj6 (jatmcls,iatmcls)

          dij1    = rtmp(1) - coord_pbc(1,ij,1)
          dij2    = rtmp(2) - coord_pbc(2,ij,1)
          dij3    = rtmp(3) - coord_pbc(3,ij,1)

          jqtmp   = coord_pbc(4,ij,1)

          rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3
          rij2   = cutoff2 * density / rij2
          L      = int(rij2)
          R      = rij2 - L
          L1     = 3*L - 2

!         work(1)  = table_grad(L1)   + R*(table_grad(L1+3)-table_grad(L1))
!         work(2)  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
!         work(3)  = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
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
          force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
          force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
          force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)
          viri(1) = viri(1) + dij1*work(1)
          viri(2) = viri(2) + dij2*work(2)
          viri(3) = viri(3) + dij3*work(3)

        end do

        force(1,ki+ix,1,id+1) = force(1,ki+ix,1,id+1) + force_local(1)
        force(2,ki+ix,1,id+1) = force(2,ki+ix,1,id+1) + force_local(2)
        force(3,ki+ix,1,id+1) = force(3,ki+ix,1,id+1) + force_local(3)
        virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
        virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
        virial(3,3,id+1) = virial(3,3,id+1) - viri(3)

      end do

    end do

    !$omp end parallel

#ifdef PKTIMER
    call timer_end(11)
#ifdef FJ_PROF_FAPP
    call fapp_stop("Nonb15F",11,0)
#endif
    call gettod(eae)
    Timc(1)=Timc(1)+(eae-sas)
#endif

    return

  end subroutine compute_force_nonbond_tbl_lnr_fugaku_npt

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_tbl_lnr_fugaku_nonpt
  !> @brief        calculate nonbonded force without solvents with charmm
  !! @authors      JJ 
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] virial   : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_tbl_lnr_fugaku_nonpt( &
                                                 domain, enefunc, pairlist, &
                                                 atmcls_pbc,                &
                                                 coord_pbc, force)

#ifdef PKTIMER
    use Ctim
#endif
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
    integer,          pointer,contiguous :: natom(:)
    integer,          pointer,contiguous :: num_nb15_calc(:,:)
    integer,          pointer,contiguous :: nb15_calc_list(:,:,:)

    real(dp)                  :: sas,eae

    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
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

#ifdef PKTIMER
    call gettod(sas)
#ifdef FJ_PROF_FAPP
    call fapp_start("Nonb15F_novirial",12,0)
#endif
    call timer_sta(12)
#endif

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
      ki = (i-1)*MaxAtom
      do ix = 1, natom(i)
        kki = ki + ix
        coord_pbc(1,kki,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(2,kki,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(3,kki,1) = coord(3,ix,i) + trans1(3,ix,i)
        coord_pbc(4,kki,1) = charge(ix,i)
        atmcls_pbc(kki)    = atmcls(ix,i)
      end do
    end do

    !$omp barrier

    do i = id+1, ncell, nthread

      ki = (i-1)*MaxAtom

      do ix = 1, natom(i)

        kki = ki + ix
        rtmp(1) = coord_pbc(1,kki,1)
        rtmp(2) = coord_pbc(2,kki,1)
        rtmp(3) = coord_pbc(3,kki,1)
        qtmp    = coord_pbc(4,kki,1)
        iatmcls = atmcls_pbc(kki)
        force_local(1) = 0.0_wp
        force_local(2) = 0.0_wp
        force_local(3) = 0.0_wp

!ocl norecurrence       

        do k = 1, num_nb15_calc(ix,i)

          ! prefetch for this ix loop prefetch
          !
          if( mod(k,64)==1 .and. (k+128)<=num_nb15_calc(ix,i) ) then
!ocl prefetch_read(nb15_calc_list(k+128,  ix,i),level=1,strong=1)
          endif

          ! prefetch for next ix loop prefetch
          if( ix < natom(i) ) then
            if((num_nb15_calc(ix,i)<=64.and.k==1) .or. &
               k==(num_nb15_calc(ix,i)-64) ) then
!ocl prefetch_read(nb15_calc_list(1,  ix+1,i),level=1,strong=1)
!ocl prefetch_read(nb15_calc_list(65, ix+1,i),level=1,strong=1)
            endif
          endif

          ij = nb15_calc_list(k+64,ix,i)

!ocl prefetch_read(coord_pbc(1,ij,1)  ,level=1,strong=1)
!ocl prefetch_read(atmcls_pbc(ij)     ,level=1,strong=1)
!ocl prefetch_write(force(1,ij,1,id+1),level=1,strong=1)

          ij = nb15_calc_list(k,ix,i)

          jatmcls = atmcls_pbc(ij)
          lj12    = nonb_lj12(jatmcls,iatmcls)
          lj6     = nonb_lj6 (jatmcls,iatmcls)

          dij(1)  = rtmp(1) - coord_pbc(1,ij,1)
          dij(2)  = rtmp(2) - coord_pbc(2,ij,1)
          dij(3)  = rtmp(3) - coord_pbc(3,ij,1)
          jqtmp   = coord_pbc(4,ij,1)

          rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2   = cutoff2 * density / rij2
          L      = int(rij2)
          R      = rij2 - L
          L1     = 3*L - 2

!         work(1)  = table_grad(L1)   + R*(table_grad(L1+3)-table_grad(L1))
!         work(2)  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
!         work(3)  = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
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
          force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
          force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
          force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

        end do

        force(1,kki,1,id+1) = force(1,kki,1,id+1) + force_local(1)
        force(2,kki,1,id+1) = force(2,kki,1,id+1) + force_local(2)
        force(3,kki,1,id+1) = force(3,kki,1,id+1) + force_local(3)

      end do

    end do

    !$omp end parallel

#ifdef PKTIMER
    call timer_end(12)
#ifdef FJ_PROF_FAPP
    call fapp_stop("Nonb15F_novirial",12,0)
#endif
    call gettod(eae)
    Timc(1)=Timc(1)+(eae-sas)
#endif

    return

  end subroutine compute_force_nonbond_tbl_lnr_fugaku_nonpt

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_notbl_fugaku
  !> @brief        calculate nonbonded energy with lookup table for elec
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

  subroutine compute_energy_nonbond_notbl_fugaku( &
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
    integer,          pointer :: natom(:)
    integer,          pointer :: num_nb15_calc(:,:)
    integer,          pointer :: nb15_calc_list(:,:,:)

    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
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
      ki = (i-1)*MaxAtom
      do ix = 1, natom(i)
        kki = ki + ix
        coord_pbc(1,kki,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(2,kki,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(3,kki,1) = coord(3,ix,i) + trans1(3,ix,i)
        coord_pbc(4,kki,1) = charge(ix,i)
        atmcls_pbc(kki)    = atmcls(ix,i)
      end do
    end do

    !$omp barrier

    do i = id+1, ncell, nthread

      ki = (i-1)*MaxAtom

      do ix = 1, natom(i)

        kki = ki + ix
        rtmp(1) = coord_pbc(1,kki,1)
        rtmp(2) = coord_pbc(2,kki,1)
        rtmp(3) = coord_pbc(3,kki,1)
        qtmp    = coord_pbc(4,kki,1)
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
        do k = 1, num_nb15_calc(ix,i)

          ij = nb15_calc_list(k,ix,i)

          dij(1) = rtmp(1) - coord_pbc(1,ij,1)
          dij(2) = rtmp(2) - coord_pbc(2,ij,1)
          dij(3) = rtmp(3) - coord_pbc(3,ij,1)
          rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          if (rij2 < cutoff2) within_cutoff = 1.0_wp
          if (rij2 >= cutoff2) within_cutoff = 0.0_wp
          rij2_inv = 1.0_wp / rij2
          rij2    = cutoff2*density*rij2_inv

          jatmcls = atmcls_pbc(ij)
          lj12    = nonb_lj12(jatmcls,iatmcls)
          lj6     = nonb_lj6 (jatmcls,iatmcls)
          jqtmp   = coord_pbc(4,ij,1)

          L  = int(rij2)
          R  = rij2 - L

          term_lj6  = rij2_inv * rij2_inv * rij2_inv * within_cutoff
          term_lj12 = term_lj6 * term_lj6 * within_cutoff
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
          force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
          force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
          force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

        end do

        force(1,ki+ix,1,id+1) = force(1,ki+ix,1,id+1) + force_local(1)
        force(2,ki+ix,1,id+1) = force(2,ki+ix,1,id+1) + force_local(2)
        force(3,ki+ix,1,id+1) = force(3,ki+ix,1,id+1) + force_local(3)
        virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
        virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
        virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond_notbl_fugaku

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_notbl_fugaku
  !> @brief        calculate nonbonded force without solvents with lookup table
  !! @authors      JJ 
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] virial   : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_notbl_fugaku( &
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
      call compute_force_nonbond_notbl_fugaku_npt( &  
           domain, enefunc, pairlist, atmcls_pbc,    &
           coord_pbc, force, virial)

    else
      call compute_force_nonbond_notbl_fugaku_nonpt( & 
           domain, enefunc, pairlist, atmcls_pbc,      &
           coord_pbc, force)
    end if

    return

  end subroutine compute_force_nonbond_notbl_fugaku

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_tbl_lnr_fugaku_npt
  !> @brief        calculate nonbonded force with virial
  !! @authors      JJ 
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] virial   : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_notbl_fugaku_npt( &
                                                 domain, enefunc, pairlist, &
                                                 atmcls_pbc,                &
                                                 coord_pbc, force, virial)

#ifdef PKTIMER
    use Ctim
#endif
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
    integer                   :: i, ix, iy, j, k, ki, ij, L, L1
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
    integer,          pointer,contiguous :: natom(:)
    integer,          pointer,contiguous :: num_nb15_calc(:,:)
    integer,          pointer,contiguous :: nb15_calc_list(:,:,:)

    real(dp)                  :: sas,eae

    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
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
    !$omp         dij3, jqtmp, viri, within_cutoff, rij2_inv)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell, nthread
      ki = (i-1)*MaxAtom
      do ix = 1, natom(i)
        coord_pbc(1,ki+ix,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(2,ki+ix,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(3,ki+ix,1) = coord(3,ix,i) + trans1(3,ix,i)
        coord_pbc(4,ki+ix,1) = charge(ix,i)
        atmcls_pbc(ki+ix)    = atmcls(ix,i)
      end do
    end do

    !$omp barrier

    do i = id+1, ncell, nthread

      ki = (i-1)*MaxAtom

      do ix = 1, natom(i)

        rtmp(1) = coord_pbc(1,ki+ix,1)
        rtmp(2) = coord_pbc(2,ki+ix,1)
        rtmp(3) = coord_pbc(3,ki+ix,1)
        qtmp    = coord_pbc(4,ki+ix,1)
        iatmcls = atmcls_pbc(ki+ix)
        force_local(1) = 0.0_wp
        force_local(2) = 0.0_wp
        force_local(3) = 0.0_wp
        viri(1) = 0.0_wp
        viri(2) = 0.0_wp
        viri(3) = 0.0_wp

!ocl norecurrence

        do k = 1, num_nb15_calc(ix,i)

          ! prefetch for this ix loop prefetch
          !
          if( mod(k,64)==1 .and. (k+128)<=num_nb15_calc(ix,i) ) then
!ocl prefetch_read(nb15_calc_list(k+128,  ix,i),level=1,strong=1)
          endif

          ! prefetch for next ix loop prefetch
          if( ix < natom(i) ) then
            if((num_nb15_calc(ix,i)<=64.and.k==1) .or. &
               k==(num_nb15_calc(ix,i)-64) ) then
!ocl prefetch_read(nb15_calc_list(1,  ix+1,i),level=1,strong=1)
!ocl prefetch_read(nb15_calc_list(65, ix+1,i),level=1,strong=1)
            endif
          endif

          ij = nb15_calc_list(k+64,ix,i)

!ocl prefetch_read(coord_pbc(1,ij,1)  ,level=1,strong=1)
!ocl prefetch_read(atmcls_pbc(ij)     ,level=1,strong=1)
!ocl prefetch_write(force(1,ij,1,id+1),level=1,strong=1)

          ij = nb15_calc_list(k,ix,i)

          jatmcls = atmcls_pbc(ij)
          lj12    = nonb_lj12(jatmcls,iatmcls)
          lj6     = nonb_lj6 (jatmcls,iatmcls)

          dij1    = rtmp(1) - coord_pbc(1,ij,1)
          dij2    = rtmp(2) - coord_pbc(2,ij,1)
          dij3    = rtmp(3) - coord_pbc(3,ij,1)

          jqtmp   = coord_pbc(4,ij,1)

          rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3
          if (rij2 < cutoff2) within_cutoff = 1.0_wp
          if (rij2 >= cutoff2) within_cutoff = 0.0_wp
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
          force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
          force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
          force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)
          viri(1) = viri(1) + dij1*work(1)
          viri(2) = viri(2) + dij2*work(2)
          viri(3) = viri(3) + dij3*work(3)

        end do

        force(1,ki+ix,1,id+1) = force(1,ki+ix,1,id+1) + force_local(1)
        force(2,ki+ix,1,id+1) = force(2,ki+ix,1,id+1) + force_local(2)
        force(3,ki+ix,1,id+1) = force(3,ki+ix,1,id+1) + force_local(3)
        virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
        virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
        virial(3,3,id+1) = virial(3,3,id+1) - viri(3)

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_force_nonbond_notbl_fugaku_npt

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_tbl_lnr_fugaku_nonpt
  !> @brief        calculate nonbonded force without solvents with charmm
  !! @authors      JJ 
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] virial   : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_notbl_fugaku_nonpt( &
                                                 domain, enefunc, pairlist, &
                                                 atmcls_pbc,                &
                                                 coord_pbc, force)

#ifdef PKTIMER
    use Ctim
#endif
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
    integer,          pointer,contiguous :: natom(:)
    integer,          pointer,contiguous :: num_nb15_calc(:,:)
    integer,          pointer,contiguous :: nb15_calc_list(:,:,:)

    real(dp)                  :: sas,eae

    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
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
      ki = (i-1)*MaxAtom
      do ix = 1, natom(i)
        kki = ki + ix
        coord_pbc(1,kki,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(2,kki,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(3,kki,1) = coord(3,ix,i) + trans1(3,ix,i)
        coord_pbc(4,kki,1) = charge(ix,i)
        atmcls_pbc(kki)    = atmcls(ix,i)
      end do
    end do

    !$omp barrier

    do i = id+1, ncell, nthread

      ki = (i-1)*MaxAtom

      do ix = 1, natom(i)

        kki = ki + ix
        rtmp(1) = coord_pbc(1,kki,1)
        rtmp(2) = coord_pbc(2,kki,1)
        rtmp(3) = coord_pbc(3,kki,1)
        qtmp    = coord_pbc(4,kki,1)
        iatmcls = atmcls_pbc(kki)

        force_local(1) = 0.0_wp
        force_local(2) = 0.0_wp
        force_local(3) = 0.0_wp

!ocl norecurrence       

        do k = 1, num_nb15_calc(ix,i)

          ! prefetch for this ix loop prefetch
          !
          if( mod(k,64)==1 .and. (k+128)<=num_nb15_calc(ix,i) ) then
!ocl prefetch_read(nb15_calc_list(k+128,  ix,i),level=1,strong=1)
          endif

          ! prefetch for next ix loop prefetch
          if( ix < natom(i) ) then
            if((num_nb15_calc(ix,i)<=64.and.k==1) .or. &
               k==(num_nb15_calc(ix,i)-64) ) then
!ocl prefetch_read(nb15_calc_list(1,  ix+1,i),level=1,strong=1)
!ocl prefetch_read(nb15_calc_list(65, ix+1,i),level=1,strong=1)
            endif
          endif

          ij = nb15_calc_list(k+64,ix,i)

!ocl prefetch_read(coord_pbc(1,ij,1)  ,level=1,strong=1)
!ocl prefetch_read(atmcls_pbc(ij)     ,level=1,strong=1)
!ocl prefetch_write(force(1,ij,1,id+1),level=1,strong=1)

          ij = nb15_calc_list(k,ix,i)

          jatmcls = atmcls_pbc(ij)
          lj12    = nonb_lj12(jatmcls,iatmcls)
          lj6     = nonb_lj6 (jatmcls,iatmcls)

          dij(1)  = rtmp(1) - coord_pbc(1,ij,1)
          dij(2)  = rtmp(2) - coord_pbc(2,ij,1)
          dij(3)  = rtmp(3) - coord_pbc(3,ij,1)
          jqtmp   = coord_pbc(4,ij,1)

          rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          if (rij2 < cutoff2) within_cutoff = 1.0_wp
          if (rij2 >= cutoff2) within_cutoff = 0.0_wp
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
          force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
          force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
          force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

        end do

        force(1,kki,1,id+1) = force(1,kki,1,id+1) + force_local(1)
        force(2,kki,1,id+1) = force(2,kki,1,id+1) + force_local(2)
        force(3,kki,1,id+1) = force(3,kki,1,id+1) + force_local(3)

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_force_nonbond_notbl_fugaku_nonpt

end module sp_energy_nonbond_fugaku_mod
