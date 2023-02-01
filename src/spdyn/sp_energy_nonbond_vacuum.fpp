!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_nonbond_vacuum_mod
!> @brief   calculate nonbonded energy without table for vdw and elec
!! @authors Hiraku Oshima
!  
!  (c) Copyright 2022 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_nonbond_vacuum_mod

  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use mpi_parallel_mod
  use constants_mod
  use timers_mod

  implicit none
  private

  ! subroutines
  !
  public  :: compute_energy_nonbond14_vacuum
  public  :: compute_energy_nonbond_vacuum
  public  :: compute_force_nonbond_vacuum
  private :: compute_energy_nonbond14_vacuum_charmm
  private :: compute_energy_nonbond14_vacuum_gro_amber
  ! FEP
  public  :: compute_energy_nonbond14_vacuum_fep
  public  :: compute_energy_nonbond_vacuum_fep
  public  :: compute_force_nonbond_vacuum_fep
  private :: compute_energy_nonbond14_vacuum_charmm_fep
  private :: compute_energy_nonbond14_vacuum_gro_amber_fep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_vacuum
  !> @brief        calculate nonbonded14 energy using CUTOFF
  !                without lookup table
  !  @authors      HO
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_vacuum(domain, enefunc, &
                                            force, virial, eelec, evdw)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    ! ==> Type 4
    if (enefunc%forcefield == ForcefieldCHARMM) then

      call compute_energy_nonbond14_vacuum_charmm( &
                                    domain, enefunc, &
                                    force, virial, eelec, evdw)

    ! ==> Type 10
    else  if (enefunc%forcefield == ForcefieldAMBER    .or. &
              enefunc%forcefield == ForcefieldGROAMBER .or. &
              enefunc%forcefield == ForcefieldGROMARTINI) then

      call compute_energy_nonbond14_vacuum_gro_amber( &
                                    domain, enefunc, &
                                    force, virial, eelec, evdw)

    endif

    return

  end subroutine compute_energy_nonbond14_vacuum

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_vacuum
  !> @brief        calculate nonbonded energy using CUTOFF without lookup table
  !! @authors      HO
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] virial   : virial term of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : van der Waals energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_vacuum(domain, enefunc, pairlist, &
                                                 force, virial, eelec, evdw)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(3,3,nthread)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: lj6, lj12
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, term_lj, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: viri_local(1:3,1:3)
    real(wp)                  :: trans_x, trans_y, trans_z
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: elec_temp, evdw_temp
    real(wp)                  :: dij_list(3,MaxAtom)
    real(wp)                  :: force_local(3), rij2_list(MaxAtom)
    real(wp)                  :: force_localj(3,MaxAtom)
    integer                   :: i, ix, iy, j, k, ij, iix, iwater, list, L, L1
    integer                   :: num_nb15, ini_nb15, fin_nb15, num_count
    integer                   :: id, omp_get_thread_num, j_list(MaxAtom)
    integer                   :: ncell, ncell_local
    real(wp)                  :: el_fact, rij2_inv, rij_inv

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), trans2(:,:,:), charge(:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: cutoff
    integer(1),       pointer :: cell_move(:,:,:)
    integer(int2),    pointer :: cell_pairlist(:,:)
    integer,          pointer :: atmcls(:,:)
    integer,          pointer :: natom(:), nsolute(:), nwater(:)
    integer,          pointer :: nb15_list(:,:)
    integer,          pointer :: nb15_cell(:)
    integer,          pointer :: num_nb15_calc1(:,:), num_nb15_calc(:,:)
    integer,          pointer :: nb15_calc_list1(:,:), nb15_calc_list(:,:)


    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    nsolute         => domain%num_solute
    nwater          => domain%num_water
    coord           => domain%coord
    trans1          => domain%trans_vec
    trans2          => domain%translated
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size

    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6

    nb15_list       => pairlist%nb15_list
    nb15_cell       => pairlist%nb15_cell
    num_nb15_calc1  => pairlist%num_nb15_calc1
    num_nb15_calc   => pairlist%num_nb15_calc
    nb15_calc_list1 => pairlist%nb15_calc_list1
    nb15_calc_list  => pairlist%nb15_calc_list

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff

    el_fact = ELECOEF / enefunc%dielec_const

    ! calculate energy and gradient
    !

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, ix, num_nb15, rtmp, qtmp, ini_nb15, fin_nb15, k, iy,  &
    !$omp         rij2, term_lj12, term_lj6, term_lj, grad_coef, work,   &
    !$omp         term_elec, viri_local, iwater, ij, j, trans_x, trans_y,      &
    !$omp         trans_z, iix, list, num_count, j_list, dij, dij_list,        &
    !$omp         rij2_list, force_local, force_localj, lj6, lj12, L1,         &
    !$omp         elec_temp,evdw_temp, rij2_inv, rij_inv)             
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread
      do ix = 1, natom(i)
        trans2(ix,1,i) = coord(1,ix,i) + trans1(1,ix,i)
        trans2(ix,2,i) = coord(2,ix,i) + trans1(2,ix,i)
        trans2(ix,3,i) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    !$omp barrier

    ! energy within a cell(solute-solute)
    !
    do i = id+1, ncell_local, nthread

      num_nb15 = 0
      do ix = 1, natom(i) - 1

        rtmp(1:3) = trans2(ix,1:3,i)
        qtmp      = charge(ix,i)
        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc1(ix,i)
        num_nb15  = fin_nb15
        num_count = 0

        viri_local(1:3,1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp

        do k = ini_nb15, fin_nb15

          iy = nb15_calc_list1(k,i)

          ! compute distance
          !
          dij(1) = rtmp(1) - trans2(iy,1,i)
          dij(2) = rtmp(2) - trans2(iy,2,i)
          dij(3) = rtmp(3) - trans2(iy,3,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1:3,num_count)  = dij(1:3)
            rij2_list(num_count) = rij2
            j_list(num_count)   = iy
          end if
        end do

        force_local(1:3) = 0.0_wp

!ocl norecurrence(force)
        do k = 1, num_count
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          iy     = j_list(k)
          rij2_inv = 1.0_wp / rij2_list(k)
          rij_inv  = sqrt(rij2_inv)
          lj6    = nonb_lj6(atmcls(ix,i),atmcls(iy,i))
          lj12   = nonb_lj12(atmcls(ix,i),atmcls(iy,i))

          term_lj6 = rij2_inv * rij2_inv * rij2_inv
          term_lj12 = term_lj6 * term_lj6
          term_elec = el_fact*rij_inv
          evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
          elec_temp = elec_temp + qtmp*charge(iy,i)*term_elec

          term_lj12 = -12.0_wp*term_lj12*rij2_inv
          term_lj6  = - 6.0_wp*term_lj6 *rij2_inv
          term_elec = -term_elec*rij2_inv
          grad_coef = term_lj12*lj12 - term_lj6*lj6                           &
                     + qtmp*charge(iy,i)*term_elec

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

          ! virial
          !
          viri_local(1,1) = viri_local(1,1) + dij(1)*work(1)
          viri_local(2,2) = viri_local(2,2) + dij(2)*work(2)
          viri_local(3,3) = viri_local(3,3) + dij(3)*work(3)
        end do

        do k = 1, num_count
          iy = j_list(k)
          force(iy,1:3,i,id+1) = force(iy,1:3,i,id+1) + force_localj(1:3,k)
        end do
        force(ix,1:3,i,id+1) = force(ix,1:3,i,id+1) + force_local(1:3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp
        do k = 1, 3
          virial(k,k,id+1) = virial(k,k,id+1) - viri_local(k,k)
        end do
       
      end do

    end do


    ! interaction between different cells
    !
    do ij = id+1, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      num_nb15 = 0

      trans_x = cell_move(1,j,i) * system_size(1)
      trans_y = cell_move(2,j,i) * system_size(2)
      trans_z = cell_move(3,j,i) * system_size(3)

      do iix = 1, nb15_cell(ij)

        ix   = nb15_list(iix,ij)
        rtmp(1:3) = trans2(ix,1:3,i)
        qtmp = charge(ix,i)

        ini_nb15 = num_nb15 + 1
        fin_nb15 = num_nb15 + num_nb15_calc(ix,ij)
        num_nb15 = fin_nb15

        num_count = 0

        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp
        viri_local(1:3,1:3) = 0.0_wp

        do k = ini_nb15, fin_nb15
          iy   = nb15_calc_list(k,ij)
          dij(1) = rtmp(1) - trans2(iy,1,j) + trans_x
          dij(2) = rtmp(2) - trans2(iy,2,j) + trans_y
          dij(3) = rtmp(3) - trans2(iy,3,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1:3,num_count)  = dij(1:3)
            rij2_list(num_count) = rij2
            j_list(num_count)   = iy
          end if
        end do

        force_local(1:3) = 0.0_wp

!ocl norecurrence(force)
        do k = 1, num_count
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          iy   = j_list(k)
          rij2_inv = 1.0_wp / rij2_list(k)
          rij_inv  = sqrt(rij2_inv)
          lj6    = nonb_lj6(atmcls(ix,i),atmcls(iy,j))
          lj12   = nonb_lj12(atmcls(ix,i),atmcls(iy,j))

          term_lj6 = rij2_inv * rij2_inv * rij2_inv
          term_lj12 = term_lj6 * term_lj6
          term_elec = el_fact*rij_inv
          evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
          elec_temp = elec_temp + qtmp*charge(iy,j)*term_elec

          term_lj12 = -12.0_wp*term_lj12*rij2_inv
          term_lj6  = - 6.0_wp*term_lj6 *rij2_inv
          term_elec = -term_elec*rij2_inv
          grad_coef = term_lj12*lj12 - term_lj6*lj6                           &
                     + qtmp*charge(iy,j)*term_elec

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

          ! virial
          !
          viri_local(1,1) = viri_local(1,1) + dij(1)*work(1)
          viri_local(2,1) = viri_local(2,1) + dij(2)*work(1)
          viri_local(2,2) = viri_local(2,2) + dij(2)*work(2)
          viri_local(3,1) = viri_local(3,1) + dij(3)*work(1)
          viri_local(3,2) = viri_local(3,2) + dij(3)*work(2)
          viri_local(3,3) = viri_local(3,3) + dij(3)*work(3)
        end do

        do k = 1, num_count
          iy = j_list(k)
          force(iy,1:3,j,id+1) = force(iy,1:3,j,id+1) + force_localj(1:3,k)
        end do
        force(ix,1:3,i,id+1) = force(ix,1:3,i,id+1) + force_local(1:3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp
        do k = 1, 3
          virial(k,k,id+1) = virial(k,k,id+1) - viri_local(k,k)
        end do

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond_vacuum

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_vacuum
  !> @brief        calculate nonbonded force using CUTOFF without lookup table
  !! @authors      HO
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] virial   : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_vacuum(domain, enefunc, pairlist, &
                                                force, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(3,3,nthread)

    ! local variables
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: lj12, lj6
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: viri_local(1:3,1:3)
    real(wp)                  :: trans_x, trans_y, trans_z
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: force_local(3), force_local_iy(3,MaxAtom)
    real(wp)                  :: dij_list(4,MaxAtom)
    integer                   :: i, ix, iy, j, k, ij, iix, iwater, list, L, L1
    integer                   :: num_nb15, ini_nb15, fin_nb15, num_count
    integer                   :: id, omp_get_thread_num, j_list(MaxAtom)
    integer                   :: iatmcls,jatmcls
    integer                   :: ncell, ncell_local
    real(wp)                  :: el_fact, rij2_inv, rij_inv

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), trans2(:,:,:), charge(:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: cutoff
    integer(1),       pointer :: cell_move(:,:,:)
    integer(int2),    pointer :: cell_pairlist(:,:)
    integer,          pointer :: atmcls(:,:)
    integer,          pointer :: natom(:), nsolute(:), nwater(:)
    integer,          pointer :: nb15_list(:,:)
    integer,          pointer :: nb15_cell(:)
    integer,          pointer :: num_nb15_calc1(:,:), num_nb15_calc(:,:)
    integer,          pointer :: nb15_calc_list1(:,:), nb15_calc_list(:,:)


    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    nsolute         => domain%num_solute
    nwater          => domain%num_water
    coord           => domain%coord
    trans1          => domain%trans_vec
    trans2          => domain%translated
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size

    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6 

    nb15_list       => pairlist%nb15_list
    nb15_cell       => pairlist%nb15_cell
    num_nb15_calc1  => pairlist%num_nb15_calc1
    num_nb15_calc   => pairlist%num_nb15_calc 
    nb15_calc_list1 => pairlist%nb15_calc_list1
    nb15_calc_list  => pairlist%nb15_calc_list

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff

    el_fact = ELECOEF / enefunc%dielec_const

    ! calculate energy and gradient
    !

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, ix, num_nb15, rtmp, qtmp, ini_nb15, fin_nb15, k, iy,  &
    !$omp         rij2, term_lj12, term_lj6, grad_coef, work, term_elec, &
    !$omp         viri_local, iwater, ij, j, trans_x, trans_y, trans_z,        &
    !$omp         iix, list, num_count, j_list, force_local, force_local_iy,   &
    !$omp         lj12, lj6, L1, dij, dij_list, iatmcls, jatmcls, rij2_inv,    &
    !$omp         rij_inv) 
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell, nthread
      do ix = 1, natom(i)
        trans2(ix,1,i) = coord(1,ix,i) + trans1(1,ix,i)
        trans2(ix,2,i) = coord(2,ix,i) + trans1(2,ix,i)
        trans2(ix,3,i) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    !$omp barrier

    ! energy within a cell(solute-solute)
    !
    do i = id+1, ncell_local, nthread

      num_nb15 = 0

      do ix = 1, natom(i) - 1

        rtmp(1:3) = trans2(ix,1:3,i)
        qtmp = charge(ix,i)
        iatmcls = atmcls(ix,i)

        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc1(ix,i)
        num_nb15  = fin_nb15
        force_local(1:3) = 0.0_wp
        num_count = 0
        viri_local(1:3,1:3) = 0.0_wp

!ocl norecurrence(force)
        do k = ini_nb15, fin_nb15

          iy = nb15_calc_list1(k,i)
          jatmcls = atmcls(iy,i)
          dij(1) = rtmp(1) - trans2(iy,1,i)
          dij(2) = rtmp(2) - trans2(iy,2,i)
          dij(3) = rtmp(3) - trans2(iy,3,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1,num_count) = dij(1)
            dij_list(2,num_count) = dij(2)
            dij_list(3,num_count) = dij(3)
            dij_list(4,num_count) = rij2
            j_list(num_count) = iy
          end if
        end do
        force_local(1:3) = 0.0_wp

        do k = 1, num_count
          iy   = j_list(k)
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          rij2_inv = 1.0_wp / dij_list(4,k)
          rij_inv = sqrt(rij2_inv)

          jatmcls = atmcls(iy,i)
          lj12 = nonb_lj12(iatmcls,jatmcls)
          lj6  = nonb_lj6 (iatmcls,jatmcls)

          term_lj6 = rij2_inv * rij2_inv * rij2_inv
          term_lj12 = term_lj6 * term_lj6
          term_lj12 = -12.0_wp*term_lj12*rij2_inv
          term_lj6  = - 6.0_wp*term_lj6 *rij2_inv
          term_elec = -el_fact*rij_inv*rij2_inv
          grad_coef = term_lj12*lj12 - term_lj6*lj6                           &
                     + qtmp*charge(iy,i)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force_local_iy(1,k) = work(1)
          force_local_iy(2,k) = work(2)
          force_local_iy(3,k) = work(3)

          ! virial
          !
          viri_local(1,1) = viri_local(1,1) + dij(1)*work(1)
          viri_local(2,1) = viri_local(2,1) + dij(2)*work(1)
          viri_local(2,2) = viri_local(2,2) + dij(2)*work(2)
          viri_local(3,1) = viri_local(3,1) + dij(3)*work(1)
          viri_local(3,2) = viri_local(3,2) + dij(3)*work(2)
          viri_local(3,3) = viri_local(3,3) + dij(3)*work(3)
        end do

        force(ix,1:3,i,id+1) = force(ix,1:3,i,id+1) + force_local(1:3)
        do k = 1, 3
          virial(k,k,id+1) = virial(k,k,id+1) - viri_local(k,k)
        end do

        do k = 1, num_count
          iy = j_list(k) 
          num_count = k - ini_nb15 + 1
          force(iy,1:3,i,id+1) = force(iy,1:3,i,id+1) + force_local_iy(1:3,k)
        end do
      end do

    end do

    ! interaction between different cells
    !
    do ij = id+1, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      num_nb15 = 0

      trans_x = cell_move(1,j,i) * system_size(1)
      trans_y = cell_move(2,j,i) * system_size(2)
      trans_z = cell_move(3,j,i) * system_size(3)

      do iix = 1, nb15_cell(ij)

        ix = nb15_list(iix,ij)
        rtmp(1:3) = trans2(ix,1:3,i)
        qtmp = charge(ix,i)
        iatmcls = atmcls(ix,i)

        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc(ix,ij)
        num_nb15  = fin_nb15

        force_local(1:3) = 0.0_wp
        num_count = 0

        viri_local(1:3,1:3) = 0.0_wp

!ocl norecurrence(force)
        do k = ini_nb15, fin_nb15
          iy = nb15_calc_list(k,ij)

          dij(1) = rtmp(1) - trans2(iy,1,j) + trans_x
          dij(2) = rtmp(2) - trans2(iy,2,j) + trans_y
          dij(3) = rtmp(3) - trans2(iy,3,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1,num_count) = dij(1)
            dij_list(2,num_count) = dij(2)
            dij_list(3,num_count) = dij(3)
            dij_list(4,num_count) = rij2
            j_list(num_count) = iy
          end if
        end do

        force_local(1:3) = 0.0_wp
        do k = 1, num_count
          iy = j_list(k)
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          rij2_inv = 1.0_wp / dij_list(4,k)
          rij_inv = sqrt(rij2_inv)
          jatmcls = atmcls(iy,j)
          lj12 = nonb_lj12(iatmcls,jatmcls)
          lj6  = nonb_lj6 (iatmcls,jatmcls)

          term_lj6 = rij2_inv * rij2_inv * rij2_inv
          term_lj12 = term_lj6 * term_lj6
          term_lj12 = -12.0_wp*term_lj12*rij2_inv
          term_lj6  = - 6.0_wp*term_lj6 *rij2_inv
          term_elec = -el_fact*rij_inv*rij2_inv
          grad_coef = term_lj12*lj12 - term_lj6*lj6                           &
                     + qtmp*charge(iy,j)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force_local_iy(1,k) = work(1)
          force_local_iy(2,k) = work(2)
          force_local_iy(3,k) = work(3)
          ! virial
          !
          viri_local(1,1) = viri_local(1,1) + dij(1)*work(1)
          viri_local(2,1) = viri_local(2,1) + dij(2)*work(1)
          viri_local(2,2) = viri_local(2,2) + dij(2)*work(2)
          viri_local(3,1) = viri_local(3,1) + dij(3)*work(1)
          viri_local(3,2) = viri_local(3,2) + dij(3)*work(2)
          viri_local(3,3) = viri_local(3,3) + dij(3)*work(3)
        end do

        do k = 1, num_count
          iy = j_list(k)
          force(iy,1:3,j,id+1) = force(iy,1:3,j,id+1) + force_local_iy(1:3,k)
        end do
        force(ix,1:3,i,id+1) = force(ix,1:3,i,id+1) + force_local(1:3)
        do k = 1, 3
          virial(k,k,id+1) = virial(k,k,id+1) - viri_local(k,k)
        end do

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_force_nonbond_vacuum

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_vacuum_charmm
  !> @brief        calculate nonbonded14 energy using CUTOFF
  !                without lookup table
  !  @authors      HO
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_vacuum_charmm(domain, enefunc, &
                                            force, virial, eelec, evdw)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                 :: dij(1:3), rij2
    real(wp)                 :: term_lj12, term_lj6, term_elec
    real(wp)                 :: cutoff2, grad_coef, lj12, lj6
    real(wp)                 :: work(1:3), elec_temp, evdw_temp
    real(wp)                 :: viri(1:3)
    integer                  :: i, ix, iy, j, k, ij, m, L
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: ncell_local, id, omp_get_thread_num
    integer                  :: iatmcls, jatmcls
    real(wp)                 :: el_fact, rij2_inv, rij_inv

    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:), system_size(:)
    real(wp),        pointer :: nb14_lj12(:,:), nb14_lj6(:,:)
    real(wp),        pointer :: cutoff
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: atmcls(:,:), natom(:)
    integer,         pointer :: num_nb14_calc(:), nb14_calc_list(:,:,:)
    integer(1),      pointer :: cell_move(:,:,:)


    coord           => domain%coord
    charge          => domain%charge
    cell_pair       => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    cell_move       => domain%cell_move
    system_size     => domain%system_size

    cutoff          => enefunc%cutoffdist
    nb14_lj12       => enefunc%nb14_lj12
    nb14_lj6        => enefunc%nb14_lj6
    num_nb14_calc   => enefunc%num_nb14_calc
    nb14_calc_list  => enefunc%nb14_calc_list

    cutoff2         =  cutoff * cutoff
    ncell_local     =  domain%num_cell_local

    el_fact = ELECOEF / enefunc%dielec_const

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, j, k, ij, ix, iy, m, iatmcls, jatmcls,    &
    !$omp         dij, rij2, term_lj12, term_lj6, term_elec, grad_coef,  &
    !$omp         work, viri, lj12, lj6, rij2_inv, rij_inv,      &
    !$omp         elec_temp, evdw_temp)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    ! calculate energy and gradient
    !

    do ij = id+1, ncell_local, nthread

      viri(1:3) = 0.0_wp
      elec_temp = 0.0_wp
      evdw_temp = 0.0_wp

      do k = 1, num_nb14_calc(ij)

        i  = nb14_calc_list(1,k,ij)
        j  = nb14_calc_list(2,k,ij)
        ix = nb14_calc_list(3,k,ij)
        iy = nb14_calc_list(4,k,ij)

        iatmcls = atmcls(ix,i)
        jatmcls = atmcls(iy,j)
        lj6   = nb14_lj6 (iatmcls,jatmcls)
        lj12  = nb14_lj12(iatmcls,jatmcls)

        ! compute distance
        !
        dij(1) = coord(1,ix,i) - coord(1,iy,j)
        dij(2) = coord(2,ix,i) - coord(2,iy,j)
        dij(3) = coord(3,ix,i) - coord(3,iy,j)
        dij(1) = dij(1) + real(cell_move(1,j,i),wp) * system_size(1)
        dij(2) = dij(2) + real(cell_move(2,j,i),wp) * system_size(2)
        dij(3) = dij(3) + real(cell_move(3,j,i),wp) * system_size(3)
        rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        ! energy and gradient
        !
        rij2_inv = 1.0_wp / rij2
        rij_inv  = sqrt(rij2_inv)
        term_lj6 = rij2_inv * rij2_inv * rij2_inv
        term_lj12 = term_lj6 * term_lj6
        term_elec = el_fact*rij_inv
        evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
        elec_temp = elec_temp + charge(ix,i)*charge(iy,j)*term_elec

        term_lj12 = -12.0_wp*term_lj12*rij2_inv
        term_lj6  = - 6.0_wp*term_lj6 *rij2_inv
        term_elec = -term_elec*rij2_inv
        grad_coef = term_lj12*lj12 - term_lj6*lj6 + &
                    charge(ix,i)*charge(iy,j)*term_elec

        work(1:3) = grad_coef*dij(1:3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force(ix,1:3,i,id+1) = force(ix,1:3,i,id+1) - work(1:3)
        force(iy,1:3,j,id+1) = force(iy,1:3,j,id+1) + work(1:3)
      end do

      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
      eelec(id+1) = eelec(id+1) + elec_temp
      evdw(id+1) = evdw(id+1) + evdw_temp
      
    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond14_vacuum_charmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_vacuum_gro_amber
  !> @brief        calculate nonbonded14 energy using CUTOFF
  !                without lookup table
  !  @authors      HO
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_vacuum_gro_amber(domain, enefunc, &
                                            force, virial, eelec, evdw)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                 :: dij(1:3), rij2
    real(wp)                 :: term_lj12, term_lj6, term_elec
    real(wp)                 :: cutoff2, grad_coef, lj12, lj6
    real(wp)                 :: work(1:3), elec_temp, evdw_temp
    real(wp)                 :: viri(1:3)
    real(wp)                 :: lj_scale, qq_scale, cc
    integer                  :: i, ix, iy, j, k, ij, m, L
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: ncell_local, id, omp_get_thread_num
    real(wp)                 :: el_fact, rij2_inv, rij_inv
    integer                  :: iatmcls, jatmcls

    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:), system_size(:)
    real(wp),        pointer :: nb14_lj12(:,:), nb14_lj6(:,:)
    real(wp),        pointer :: cutoff
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: atmcls(:,:), natom(:)
    integer,         pointer :: num_nb14_calc(:), nb14_calc_list(:,:,:)
    integer(1),      pointer :: cell_move(:,:,:)


    coord           => domain%coord
    charge          => domain%charge
    cell_pair       => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    cell_move       => domain%cell_move
    system_size     => domain%system_size

    cutoff          => enefunc%cutoffdist
    nb14_lj12       => enefunc%nb14_lj12
    nb14_lj6        => enefunc%nb14_lj6
    num_nb14_calc   => enefunc%num_nb14_calc
    nb14_calc_list  => enefunc%nb14_calc_list

    cutoff2         =  cutoff * cutoff
    ncell_local     =  domain%num_cell_local

    el_fact = ELECOEF / enefunc%dielec_const

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, j, k, ij, ix, iy, m, iatmcls, jatmcls,    &
    !$omp         dij, rij2, term_lj12, term_lj6, term_elec, grad_coef,  &
    !$omp         work, viri, lj12, lj6, rij2_inv, rij_inv,      &
    !$omp         qq_scale, lj_scale, cc, elec_temp, evdw_temp)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    ! calculate energy and gradient
    !
    do ij = id+1, ncell_local, nthread

      viri(1:3) = 0.0_wp
      elec_temp = 0.0_wp
      evdw_temp = 0.0_wp

      do k = 1, num_nb14_calc(ij)

        i  = nb14_calc_list(1,k,ij)
        j  = nb14_calc_list(2,k,ij)
        ix = nb14_calc_list(3,k,ij)
        iy = nb14_calc_list(4,k,ij)

        iatmcls = atmcls(ix,i)
        jatmcls = atmcls(iy,j)
        qq_scale = enefunc%nb14_qq_scale(k,ij)
        lj_scale = enefunc%nb14_lj_scale(k,ij)
        lj6   = nb14_lj6 (iatmcls,jatmcls)
        lj12  = nb14_lj12(iatmcls,jatmcls)
        cc    = charge(ix,i)*charge(iy,j)*qq_scale

        ! compute distance
        !
        dij(1) = coord(1,ix,i) - coord(1,iy,j)
        dij(2) = coord(2,ix,i) - coord(2,iy,j)
        dij(3) = coord(3,ix,i) - coord(3,iy,j)
        rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        ! energy and gradient
        !
        rij2_inv = 1.0_wp / rij2
        rij_inv  = sqrt(rij2_inv)
        term_lj6 = rij2_inv * rij2_inv * rij2_inv
        term_lj12 = term_lj6 * term_lj6
        term_elec = el_fact*rij_inv
        evdw_temp = evdw_temp + (term_lj12*lj12 - term_lj6*lj6)*lj_scale
        elec_temp = elec_temp + cc*term_elec

        term_lj12 = -12.0_wp*term_lj12*rij2_inv
        term_lj6  = - 6.0_wp*term_lj6 *rij2_inv
        term_elec = -term_elec*rij2_inv
        grad_coef = (term_lj12*lj12 - term_lj6*lj6)*lj_scale + cc*term_elec

        work(1:3) = grad_coef*dij(1:3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force(ix,1:3,i,id+1) = force(ix,1:3,i,id+1) - work(1:3)
        force(iy,1:3,j,id+1) = force(iy,1:3,j,id+1) + work(1:3)

      end do

      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
      eelec(id+1) = eelec(id+1) + elec_temp
      evdw(id+1) = evdw(id+1) + evdw_temp

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond14_vacuum_gro_amber

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_vacuum_fep
  !> @brief        calculate nonbonded14 energy using CUTOFF
  !                without lookup table for FEP
  !  @authors      HO
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_vacuum_fep(domain, enefunc, &
                                            force, virial, eelec, evdw)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    ! ==> Type 4
    if (enefunc%forcefield == ForcefieldCHARMM) then

      call compute_energy_nonbond14_vacuum_charmm( &
                                    domain, enefunc, &
                                    force, virial, eelec, evdw)

      call compute_energy_nonbond14_vacuum_charmm_fep( &
                                    domain, enefunc, &
                                    force, virial, eelec, evdw)

    ! ==> Type 10
    else  if (enefunc%forcefield == ForcefieldAMBER .or. &
              enefunc%forcefield == ForcefieldGROAMBER) then

      call compute_energy_nonbond14_vacuum_gro_amber( &
                                    domain, enefunc, &
                                    force, virial, eelec, evdw)

      call compute_energy_nonbond14_vacuum_gro_amber_fep( &
                                    domain, enefunc, &
                                    force, virial, eelec, evdw)

    endif

    return

  end subroutine compute_energy_nonbond14_vacuum_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_vacuum_fep
  !> @brief        calculate nonbonded energy using CUTOFF
  !                without lookup table for FEP
  !! @authors      HO
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] virial   : virial term of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : van der Waals energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_vacuum_fep(domain, enefunc, pairlist, &
                                                 force, virial, eelec, evdw)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(3,3,nthread)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: lj6, lj12
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, term_lj, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: viri_local(1:3,1:3)
    real(wp)                  :: trans_x, trans_y, trans_z
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: elec_temp, evdw_temp
    real(wp)                  :: dij_list(3,MaxAtom)
    real(wp)                  :: force_local(3), rij2_list(MaxAtom)
    real(wp)                  :: force_localj(3,MaxAtom)
    integer                   :: i, ix, iy, j, k, ij, iix, iwater, list, L, L1
    integer                   :: num_nb15, ini_nb15, fin_nb15, num_count
    integer                   :: id, omp_get_thread_num, j_list(MaxAtom)
    integer                   :: ncell, ncell_local
    real(wp)                  :: el_fact, rij2_inv, rij_inv

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), trans2(:,:,:), charge(:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: cutoff
    integer(1),       pointer :: cell_move(:,:,:)
    integer(int2),    pointer :: cell_pairlist(:,:)
    integer,          pointer :: atmcls(:,:)
    integer,          pointer :: natom(:), nsolute(:), nwater(:)
    integer,          pointer :: nb15_list(:,:)
    integer,          pointer :: nb15_cell(:)
    integer,          pointer :: num_nb15_calc1(:,:), num_nb15_calc(:,:)
    integer,          pointer :: nb15_calc_list1(:,:), nb15_calc_list(:,:)

    ! FEP
    integer                   :: fg1, fg2
    integer,          pointer :: fepgrp(:,:)
    real(wp),         pointer :: table_sclj(:,:), table_scel(:,:)

    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    nsolute         => domain%num_solute
    nwater          => domain%num_water
    coord           => domain%coord
    trans1          => domain%trans_vec
    trans2          => domain%translated
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size

    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff

    el_fact = ELECOEF / enefunc%dielec_const

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

    ! calculate energy and gradient
    !

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, ix, num_nb15, rtmp, qtmp, ini_nb15, fin_nb15, k, iy,  &
    !$omp         rij2, term_lj12, term_lj6, term_lj, grad_coef, work,   &
    !$omp         term_elec, viri_local, iwater, ij, j, trans_x, trans_y,      &
    !$omp         trans_z, iix, list, num_count, j_list, dij, dij_list,        &
    !$omp         rij2_list, force_local, force_localj, lj6, lj12, L1,         &
    !$omp         elec_temp,evdw_temp, rij2_inv, rij_inv, fg1, fg2)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread
      do ix = 1, natom(i)
        trans2(ix,1,i) = coord(1,ix,i) + trans1(1,ix,i)
        trans2(ix,2,i) = coord(2,ix,i) + trans1(2,ix,i)
        trans2(ix,3,i) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    !$omp barrier

    ! energy within a cell(solute-solute)
    !
    do i = id+1, ncell_local, nthread

      num_nb15 = 0
      do ix = 1, natom(i) - 1

        rtmp(1:3) = trans2(ix,1:3,i)
        qtmp      = charge(ix,i)
        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc1(ix,i)
        num_nb15  = fin_nb15
        num_count = 0

        viri_local(1:3,1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp

        do k = ini_nb15, fin_nb15

          iy = nb15_calc_list1(k,i)

          ! compute distance
          !
          dij(1) = rtmp(1) - trans2(iy,1,i)
          dij(2) = rtmp(2) - trans2(iy,2,i)
          dij(3) = rtmp(3) - trans2(iy,3,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1:3,num_count)  = dij(1:3)
            rij2_list(num_count) = rij2
            j_list(num_count)   = iy
          end if
        end do

        force_local(1:3) = 0.0_wp

!ocl norecurrence(force)
        do k = 1, num_count
          iy     = j_list(k)
          lj6    = nonb_lj6(atmcls(ix,i),atmcls(iy,i))
          lj12   = nonb_lj12(atmcls(ix,i),atmcls(iy,i))

          ! FEP: Determine lamblj and lambel
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,i)

          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)

          rij2_inv = 1.0_wp / (rij2_list(k) + table_sclj(fg1,fg2))
          term_lj6 = rij2_inv * rij2_inv * rij2_inv
          term_lj12 = term_lj6 * term_lj6
          evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
          term_lj12 = -12.0_wp*term_lj12*rij2_inv
          term_lj6  = - 6.0_wp*term_lj6 *rij2_inv

          rij2_inv = 1.0_wp / (rij2_list(k) + table_scel(fg1,fg2))
          rij_inv  = sqrt(rij2_inv)
          term_elec = el_fact*rij_inv
          elec_temp = elec_temp + qtmp*charge(iy,i)*term_elec
          term_elec = -term_elec*rij2_inv

          grad_coef = term_lj12*lj12 - term_lj6*lj6                           &
                     + qtmp*charge(iy,i)*term_elec

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

          ! virial
          !
          viri_local(1,1) = viri_local(1,1) + dij(1)*work(1)
          viri_local(2,2) = viri_local(2,2) + dij(2)*work(2)
          viri_local(3,3) = viri_local(3,3) + dij(3)*work(3)
        end do

        do k = 1, num_count
          iy = j_list(k)
          force(iy,1:3,i,id+1) = force(iy,1:3,i,id+1) + force_localj(1:3,k)
        end do
        force(ix,1:3,i,id+1) = force(ix,1:3,i,id+1) + force_local(1:3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp
        do k = 1, 3
          virial(k,k,id+1) = virial(k,k,id+1) - viri_local(k,k)
        end do
       
      end do

    end do


    ! interaction between different cells
    !
    do ij = id+1, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      num_nb15 = 0

      trans_x = cell_move(1,j,i) * system_size(1)
      trans_y = cell_move(2,j,i) * system_size(2)
      trans_z = cell_move(3,j,i) * system_size(3)

      do iix = 1, nb15_cell(ij)

        ix   = nb15_list(iix,ij)
        rtmp(1:3) = trans2(ix,1:3,i)
        qtmp = charge(ix,i)

        ini_nb15 = num_nb15 + 1
        fin_nb15 = num_nb15 + num_nb15_calc(ix,ij)
        num_nb15 = fin_nb15

        num_count = 0

        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp
        viri_local(1:3,1:3) = 0.0_wp

        do k = ini_nb15, fin_nb15
          iy   = nb15_calc_list(k,ij)
          dij(1) = rtmp(1) - trans2(iy,1,j) + trans_x
          dij(2) = rtmp(2) - trans2(iy,2,j) + trans_y
          dij(3) = rtmp(3) - trans2(iy,3,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1:3,num_count)  = dij(1:3)
            rij2_list(num_count) = rij2
            j_list(num_count)   = iy
          end if
        end do

        force_local(1:3) = 0.0_wp

!ocl norecurrence(force)
        do k = 1, num_count
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          iy   = j_list(k)

          lj6    = nonb_lj6(atmcls(ix,i),atmcls(iy,j))
          lj12   = nonb_lj12(atmcls(ix,i),atmcls(iy,j))

          ! FEP: Determine lamblj and lambel
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,j)

          rij2_inv = 1.0_wp / (rij2_list(k) + table_sclj(fg1,fg2))
          term_lj6 = rij2_inv * rij2_inv * rij2_inv
          term_lj12 = term_lj6 * term_lj6
          evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
          term_lj12 = -12.0_wp*term_lj12*rij2_inv
          term_lj6  = - 6.0_wp*term_lj6 *rij2_inv

          rij2_inv = 1.0_wp / (rij2_list(k) + table_scel(fg1,fg2))
          rij_inv  = sqrt(rij2_inv)
          term_elec = el_fact*rij_inv
          elec_temp = elec_temp + qtmp*charge(iy,j)*term_elec
          term_elec = -term_elec*rij2_inv

          grad_coef = term_lj12*lj12 - term_lj6*lj6                           &
                     + qtmp*charge(iy,j)*term_elec

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

          ! virial
          !
          viri_local(1,1) = viri_local(1,1) + dij(1)*work(1)
          viri_local(2,1) = viri_local(2,1) + dij(2)*work(1)
          viri_local(2,2) = viri_local(2,2) + dij(2)*work(2)
          viri_local(3,1) = viri_local(3,1) + dij(3)*work(1)
          viri_local(3,2) = viri_local(3,2) + dij(3)*work(2)
          viri_local(3,3) = viri_local(3,3) + dij(3)*work(3)
        end do

        do k = 1, num_count
          iy = j_list(k)
          force(iy,1:3,j,id+1) = force(iy,1:3,j,id+1) + force_localj(1:3,k)
        end do
        force(ix,1:3,i,id+1) = force(ix,1:3,i,id+1) + force_local(1:3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp
        do k = 1, 3
          virial(k,k,id+1) = virial(k,k,id+1) - viri_local(k,k)
        end do

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond_vacuum_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_vacuum_fep
  !> @brief        calculate nonbonded force using CUTOFF
  !                without lookup table for FEP
  !! @authors      HO
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] virial   : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_vacuum_fep(domain, enefunc, pairlist, &
                                                force, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(3,3,nthread)

    ! local variables
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: lj12, lj6
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: viri_local(1:3,1:3)
    real(wp)                  :: trans_x, trans_y, trans_z
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: force_local(3), force_local_iy(3,MaxAtom)
    real(wp)                  :: dij_list(4,MaxAtom)
    integer                   :: i, ix, iy, j, k, ij, iix, iwater, list, L, L1
    integer                   :: num_nb15, ini_nb15, fin_nb15, num_count
    integer                   :: id, omp_get_thread_num, j_list(MaxAtom)
    integer                   :: iatmcls,jatmcls
    integer                   :: ncell, ncell_local
    real(wp)                  :: el_fact, rij2_inv, rij_inv

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), trans2(:,:,:), charge(:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: cutoff
    integer(1),       pointer :: cell_move(:,:,:)
    integer(int2),    pointer :: cell_pairlist(:,:)
    integer,          pointer :: atmcls(:,:)
    integer,          pointer :: natom(:), nsolute(:), nwater(:)
    integer,          pointer :: nb15_list(:,:)
    integer,          pointer :: nb15_cell(:)
    integer,          pointer :: num_nb15_calc1(:,:), num_nb15_calc(:,:)
    integer,          pointer :: nb15_calc_list1(:,:), nb15_calc_list(:,:)

    ! FEP
    integer                   :: fg1, fg2
    integer,          pointer :: fepgrp(:,:)
    real(wp),         pointer :: table_sclj(:,:), table_scel(:,:)

    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    nsolute         => domain%num_solute
    nwater          => domain%num_water
    coord           => domain%coord
    trans1          => domain%trans_vec
    trans2          => domain%translated
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size

    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6 

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff

    el_fact = ELECOEF / enefunc%dielec_const

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

    ! calculate energy and gradient
    !

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, ix, num_nb15, rtmp, qtmp, ini_nb15, fin_nb15, k, iy,  &
    !$omp         rij2, term_lj12, term_lj6, grad_coef, work, term_elec, &
    !$omp         viri_local, iwater, ij, j, trans_x, trans_y, trans_z,        &
    !$omp         iix, list, num_count, j_list, force_local, force_local_iy,   &
    !$omp         lj12, lj6, L1, dij, dij_list, iatmcls, jatmcls, rij2_inv,    &
    !$omp         rij_inv, fg1, fg2) 
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell, nthread
      do ix = 1, natom(i)
        trans2(ix,1,i) = coord(1,ix,i) + trans1(1,ix,i)
        trans2(ix,2,i) = coord(2,ix,i) + trans1(2,ix,i)
        trans2(ix,3,i) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    !$omp barrier

    ! energy within a cell(solute-solute)
    !
    do i = id+1, ncell_local, nthread

      num_nb15 = 0

      do ix = 1, natom(i) - 1

        rtmp(1:3) = trans2(ix,1:3,i)
        qtmp = charge(ix,i)
        iatmcls = atmcls(ix,i)

        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc1(ix,i)
        num_nb15  = fin_nb15
        force_local(1:3) = 0.0_wp
        num_count = 0
        viri_local(1:3,1:3) = 0.0_wp

!ocl norecurrence(force)
        do k = ini_nb15, fin_nb15

          iy = nb15_calc_list1(k,i)
          jatmcls = atmcls(iy,i)
          dij(1) = rtmp(1) - trans2(iy,1,i)
          dij(2) = rtmp(2) - trans2(iy,2,i)
          dij(3) = rtmp(3) - trans2(iy,3,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1,num_count) = dij(1)
            dij_list(2,num_count) = dij(2)
            dij_list(3,num_count) = dij(3)
            dij_list(4,num_count) = rij2
            j_list(num_count) = iy
          end if
        end do
        force_local(1:3) = 0.0_wp

        do k = 1, num_count
          iy   = j_list(k)
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)

          jatmcls = atmcls(iy,i)
          lj12 = nonb_lj12(iatmcls,jatmcls)
          lj6  = nonb_lj6 (iatmcls,jatmcls)

          ! FEP: Determine lamblj, lambel, and softcore
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,i)

          ! FEP: LJ with soft-core for reference
          rij2_inv = 1.0_wp / (dij_list(4,k) + table_sclj(fg1,fg2))
          term_lj6 = rij2_inv * rij2_inv * rij2_inv
          term_lj12 = term_lj6 * term_lj6
          term_lj12 = -12.0_wp*term_lj12*rij2_inv
          term_lj6  = - 6.0_wp*term_lj6 *rij2_inv

          ! FEP: EL with soft-core for reference
          rij2_inv  = 1.0_wp / (dij_list(4,k) + table_scel(fg1,fg2))
          rij_inv   = sqrt(rij2_inv)
          term_elec = -el_fact*rij_inv*rij2_inv

          grad_coef = term_lj12*lj12 - term_lj6*lj6                           &
                     + qtmp*charge(iy,i)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force_local_iy(1,k) = work(1)
          force_local_iy(2,k) = work(2)
          force_local_iy(3,k) = work(3)

          ! virial
          !
          viri_local(1,1) = viri_local(1,1) + dij(1)*work(1)
          viri_local(2,1) = viri_local(2,1) + dij(2)*work(1)
          viri_local(2,2) = viri_local(2,2) + dij(2)*work(2)
          viri_local(3,1) = viri_local(3,1) + dij(3)*work(1)
          viri_local(3,2) = viri_local(3,2) + dij(3)*work(2)
          viri_local(3,3) = viri_local(3,3) + dij(3)*work(3)
        end do

        force(ix,1:3,i,id+1) = force(ix,1:3,i,id+1) + force_local(1:3)
        do k = 1, 3
          virial(k,k,id+1) = virial(k,k,id+1) - viri_local(k,k)
        end do

        do k = 1, num_count
          iy = j_list(k) 
          num_count = k - ini_nb15 + 1
          force(iy,1:3,i,id+1) = force(iy,1:3,i,id+1) + force_local_iy(1:3,k)
        end do
      end do

    end do

    ! interaction between different cells
    !
    do ij = id+1, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      num_nb15 = 0

      trans_x = cell_move(1,j,i) * system_size(1)
      trans_y = cell_move(2,j,i) * system_size(2)
      trans_z = cell_move(3,j,i) * system_size(3)

      do iix = 1, nb15_cell(ij)

        ix = nb15_list(iix,ij)
        rtmp(1:3) = trans2(ix,1:3,i)
        qtmp = charge(ix,i)
        iatmcls = atmcls(ix,i)

        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc(ix,ij)
        num_nb15  = fin_nb15

        force_local(1:3) = 0.0_wp
        num_count = 0

        viri_local(1:3,1:3) = 0.0_wp

!ocl norecurrence(force)
        do k = ini_nb15, fin_nb15
          iy = nb15_calc_list(k,ij)

          dij(1) = rtmp(1) - trans2(iy,1,j) + trans_x
          dij(2) = rtmp(2) - trans2(iy,2,j) + trans_y
          dij(3) = rtmp(3) - trans2(iy,3,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1,num_count) = dij(1)
            dij_list(2,num_count) = dij(2)
            dij_list(3,num_count) = dij(3)
            dij_list(4,num_count) = rij2
            j_list(num_count) = iy
          end if
        end do

        force_local(1:3) = 0.0_wp
        do k = 1, num_count
          iy = j_list(k)
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)

          jatmcls = atmcls(iy,j)
          lj12 = nonb_lj12(iatmcls,jatmcls)
          lj6  = nonb_lj6 (iatmcls,jatmcls)

          ! FEP: Determine lamblj, lambel, and softcore
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,j)

          ! FEP: LJ with soft-core for reference
          rij2_inv = 1.0_wp / (dij_list(4,k) + table_sclj(fg1,fg2))
          term_lj6 = rij2_inv * rij2_inv * rij2_inv
          term_lj12 = term_lj6 * term_lj6
          term_lj12 = -12.0_wp*term_lj12*rij2_inv
          term_lj6  = - 6.0_wp*term_lj6 *rij2_inv

          ! FEP: EL with soft-core for reference
          rij2_inv = 1.0_wp / (dij_list(4,k) + table_scel(fg1,fg2))
          rij_inv = sqrt(rij2_inv)
          term_elec = -el_fact*rij_inv*rij2_inv

          ! gradient
          grad_coef = term_lj12*lj12 - term_lj6*lj6                           &
                     + qtmp*charge(iy,j)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force_local_iy(1,k) = work(1)
          force_local_iy(2,k) = work(2)
          force_local_iy(3,k) = work(3)
          ! virial
          !
          viri_local(1,1) = viri_local(1,1) + dij(1)*work(1)
          viri_local(2,1) = viri_local(2,1) + dij(2)*work(1)
          viri_local(2,2) = viri_local(2,2) + dij(2)*work(2)
          viri_local(3,1) = viri_local(3,1) + dij(3)*work(1)
          viri_local(3,2) = viri_local(3,2) + dij(3)*work(2)
          viri_local(3,3) = viri_local(3,3) + dij(3)*work(3)
        end do

        do k = 1, num_count
          iy = j_list(k)
          force(iy,1:3,j,id+1) = force(iy,1:3,j,id+1) + force_local_iy(1:3,k)
        end do
        force(ix,1:3,i,id+1) = force(ix,1:3,i,id+1) + force_local(1:3)
        do k = 1, 3
          virial(k,k,id+1) = virial(k,k,id+1) - viri_local(k,k)
        end do

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_force_nonbond_vacuum_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_vacuum_charmm_fep
  !> @brief        calculate nonbonded14 energy using CUTOFF
  !                without lookup table for FEP
  !  @authors      HO
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_vacuum_charmm_fep(domain, enefunc, &
                                            force, virial, eelec, evdw)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                 :: dij(1:3), rij2
    real(wp)                 :: term_lj12, term_lj6, term_elec
    real(wp)                 :: cutoff2, grad_coef, lj12, lj6, cc
    real(wp)                 :: work(1:3)
    real(wp)                 :: viri(1:3)
    integer                  :: i, ix, iy, j, k, ij, m, L
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: ncell_local, id, omp_get_thread_num
    integer                  :: iatmcls, jatmcls
    real(wp)                 :: el_fact, rij2_inv, rij_inv

    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:), system_size(:)
    real(wp),        pointer :: nb14_lj12(:,:), nb14_lj6(:,:)
    real(wp),        pointer :: cutoff
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: atmcls(:,:), natom(:)
    integer,         pointer :: num_nb14_calc(:), nb14_calc_list(:,:,:)
    integer(1),      pointer :: cell_move(:,:,:)

    ! FEP
    integer                  :: fg1, fg2
    integer,         pointer :: fepgrp(:,:)
    real(wp),        pointer :: table_sclj(:,:), table_scel(:,:)
    real(wp)                 :: elec_temp, evdw_temp

    coord           => domain%coord
    charge          => domain%charge
    cell_pair       => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    system_size     => domain%system_size
    cell_move       => domain%cell_move

    cutoff          => enefunc%cutoffdist
    nb14_lj12       => enefunc%nb14_lj12
    nb14_lj6        => enefunc%nb14_lj6
    cutoff2         =  cutoff * cutoff
    ncell_local     =  domain%num_cell_local

    ! FEP
    num_nb14_calc   => enefunc%num_nb14_calc_fep
    nb14_calc_list  => enefunc%nb14_calc_list_fep
    fepgrp          => domain%fepgrp
    table_sclj      => enefunc%table_sclj
    table_scel      => enefunc%table_scel

    el_fact = ELECOEF / enefunc%dielec_const

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, j, k, ij, ix, iy, m,     &
    !$omp         dij, rij2, term_lj12, term_lj6, term_elec, grad_coef,  &
    !$omp         work, viri, lj12, lj6, cc, iatmcls, jatmcls,    &
    !$omp         elec_temp, evdw_temp, rij2_inv, rij_inv,                &
    !$omp         fg1, fg2)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    ! calculate energy and gradient
    !

    do ij = id+1, ncell_local, nthread

      viri(1:3) = 0.0_wp
      elec_temp = 0.0_wp
      evdw_temp = 0.0_wp

      do k = 1, num_nb14_calc(ij)

        i  = nb14_calc_list(1,k,ij)
        j  = nb14_calc_list(2,k,ij)
        ix = nb14_calc_list(3,k,ij)
        iy = nb14_calc_list(4,k,ij)

        iatmcls = atmcls(ix,i)
        jatmcls = atmcls(iy,j)
        lj6   = nb14_lj6(iatmcls,jatmcls)
        lj12  = nb14_lj12(iatmcls,jatmcls)
        cc = charge(ix,i) * charge(iy,j)

        ! compute distance
        !
        dij(1) = coord(1,ix,i) - coord(1,iy,j)
        dij(2) = coord(2,ix,i) - coord(2,iy,j)
        dij(3) = coord(3,ix,i) - coord(3,iy,j)
        dij(1) = dij(1) + real(cell_move(1,j,i),wp) * system_size(1)
        dij(2) = dij(2) + real(cell_move(2,j,i),wp) * system_size(2)
        dij(3) = dij(3) + real(cell_move(3,j,i),wp) * system_size(3)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        ! FEP: Determine lamblj, lambel, and softcore
        fg1 = fepgrp(ix,i)
        fg2 = fepgrp(iy,j)

        ! FEP: LJ with soft core
        rij2_inv = 1.0_wp / (rij2 + table_sclj(fg1,fg2))
        term_lj6 = rij2_inv * rij2_inv * rij2_inv
        term_lj12 = term_lj6 * term_lj6
        evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
        term_lj12 = -12.0_wp*term_lj12*rij2_inv
        term_lj6  = - 6.0_wp*term_lj6 *rij2_inv

        ! FEP: EL with soft core
        rij2_inv = 1.0_wp / (rij2 + table_scel(fg1,fg2))
        rij_inv = sqrt(rij2_inv)
        term_elec = el_fact*rij_inv
        elec_temp = elec_temp + cc*term_elec
        term_elec = -term_elec*rij2_inv

        ! gradient
        grad_coef = term_lj12*lj12 - term_lj6*lj6 + cc*term_elec

        work(1:3) = grad_coef*dij(1:3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force(ix,1:3,i,id+1) = force(ix,1:3,i,id+1) - work(1:3)
        force(iy,1:3,j,id+1) = force(iy,1:3,j,id+1) + work(1:3)

      end do

      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
      eelec(id+1) = eelec(id+1) + elec_temp
      evdw(id+1)  = evdw(id+1)  + evdw_temp
  
    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond14_vacuum_charmm_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_vacuum_gro_amber_fep
  !> @brief        calculate nonbonded14 energy using CUTOFF
  !                without lookup table for FEP
  !  @authors      HO
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_vacuum_gro_amber_fep(domain, enefunc, &
                                            force, virial, eelec, evdw)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                 :: dij(1:3), rij2
    real(wp)                 :: term_lj12, term_lj6, term_elec
    real(wp)                 :: cutoff2, grad_coef, lj12, lj6
    real(wp)                 :: work(1:3)
    real(wp)                 :: viri(1:3)
    real(wp)                 :: lj_scale, qq_scale, cc
    integer                  :: i, ix, iy, j, k, ij, m, L
    integer                  :: ncell_local, id, omp_get_thread_num
    integer                  :: iatmcls, jatmcls
    real(wp)                 :: el_fact, rij2_inv, rij_inv

    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:), system_size(:)
    real(wp),        pointer :: nb14_lj12(:,:), nb14_lj6(:,:)
    real(wp),        pointer :: cutoff
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: atmcls(:,:), natom(:)
    integer,         pointer :: num_nb14_calc(:), nb14_calc_list(:,:,:)
    integer(1),      pointer :: cell_move(:,:,:)

    ! FEP
    real(wp)                 :: elec_temp, evdw_temp
    integer                  :: fg1, fg2
    integer,         pointer :: fepgrp(:,:)
    real(wp),        pointer :: table_sclj(:,:), table_scel(:,:)

    coord           => domain%coord
    charge          => domain%charge
    cell_pair       => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    system_size     => domain%system_size
    cell_move       => domain%cell_move

    cutoff          => enefunc%cutoffdist
    nb14_lj12       => enefunc%nb14_lj12
    nb14_lj6        => enefunc%nb14_lj6
    cutoff2         =  cutoff * cutoff
    ncell_local     =  domain%num_cell_local

    ! FEP
    num_nb14_calc   => enefunc%num_nb14_calc_fep
    nb14_calc_list  => enefunc%nb14_calc_list_fep
    fepgrp          => domain%fepgrp
    table_sclj      => enefunc%table_sclj
    table_scel      => enefunc%table_scel

    el_fact = ELECOEF / enefunc%dielec_const

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, j, k, ij, ix, iy, m,      &
    !$omp         dij, rij2, term_lj12, term_lj6, term_elec, grad_coef,  &
    !$omp         work, viri, lj12, lj6, iatmcls, jatmcls,     &
    !$omp         qq_scale, lj_scale, cc, elec_temp, evdw_temp,                &
    !$omp         rij2_inv, rij_inv, fg1, fg2)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    ! calculate energy and gradient
    !
    do ij = id+1, ncell_local, nthread

      elec_temp = 0.0_wp
      evdw_temp = 0.0_wp
      viri(1:3) = 0.0_wp

      do k = 1, num_nb14_calc(ij)

        i  = nb14_calc_list(1,k,ij)
        j  = nb14_calc_list(2,k,ij)
        ix = nb14_calc_list(3,k,ij)
        iy = nb14_calc_list(4,k,ij)

        iatmcls = atmcls(ix,i)
        jatmcls = atmcls(iy,j)
        qq_scale = enefunc%nb14_qq_scale_fep(k,ij)
        lj_scale = enefunc%nb14_lj_scale_fep(k,ij)
        lj6    = nb14_lj6(iatmcls,jatmcls)
        lj12   = nb14_lj12(iatmcls,jatmcls)
        cc    = charge(ix,i)*charge(iy,j)*qq_scale

        ! compute distance
        !
        dij(1) = coord(1,ix,i) - coord(1,iy,j)
        dij(2) = coord(2,ix,i) - coord(2,iy,j)
        dij(3) = coord(3,ix,i) - coord(3,iy,j)
        dij(1) = dij(1) + real(cell_move(1,j,i),wp) * system_size(1)
        dij(2) = dij(2) + real(cell_move(2,j,i),wp) * system_size(2)
        dij(3) = dij(3) + real(cell_move(3,j,i),wp) * system_size(3)
        rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        ! FEP: Determine lamblj, lambel, and softcore
        fg1 = fepgrp(ix,i)
        fg2 = fepgrp(iy,j)

        ! FEP: LJ with soft core
        rij2_inv = 1.0_wp / (rij2 + table_sclj(fg1,fg2))
        term_lj6 = rij2_inv * rij2_inv * rij2_inv
        term_lj12 = term_lj6 * term_lj6
        evdw_temp = evdw_temp + (term_lj12*lj12 - term_lj6*lj6)*lj_scale
        term_lj12 = -12.0_wp*term_lj12*rij2_inv
        term_lj6  = - 6.0_wp*term_lj6 *rij2_inv

        ! FEP: EL with soft core
        rij2_inv = 1.0_wp / (rij2 + table_scel(fg1,fg2))
        rij_inv = sqrt(rij2_inv)
        term_elec = el_fact*rij_inv
        elec_temp = elec_temp + cc*term_elec
        term_elec = -term_elec*rij2_inv

        ! gradient
        grad_coef = (term_lj12*lj12 - term_lj6*lj6)*lj_scale + &
                    cc*term_elec

        work(1:3) = grad_coef*dij(1:3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force(ix,1:3,i,id+1) = force(ix,1:3,i,id+1) - work(1:3)
        force(iy,1:3,j,id+1) = force(iy,1:3,j,id+1) + work(1:3)
      end do

      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
      eelec(id+1) = eelec(id+1) + elec_temp
      evdw(id+1)  = evdw(id+1)  + evdw_temp

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond14_vacuum_gro_amber_fep

end module sp_energy_nonbond_vacuum_mod
