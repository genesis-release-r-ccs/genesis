!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_table_cubic_mod
!> @brief   calculate nonbonded energy with table and with linear interpolation
!! @authors Jaewoon Jung(JJ)
!  
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_table_cubic_mod

  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  !
  public  :: compute_energy_nonbond14_table
  public  :: compute_energy_nonbond_table
  public  :: compute_force_nonbond_table

  private :: compute_energy_nonbond14_table_charmm
  private :: compute_energy_nonbond14_table_gro_amber

  ! FEP
  public  :: compute_energy_nonbond14_table_fep
  public  :: compute_energy_nonbond_table_fep
  public  :: compute_force_nonbond_table_fep
  private :: compute_energy_nonbond14_table_charmm_fep
  private :: compute_energy_nonbond14_table_gro_amber_fep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_table
  !> @brief        calculate nonbonded14 energy with PME and lookup table
  !!               (linear)
  !  @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_table(domain, enefunc, &
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

      call compute_energy_nonbond14_table_charmm(    &
                                    domain, enefunc, &
                                    force, virial, eelec, evdw)

    ! ==> Type 10
    else  if (enefunc%forcefield == ForcefieldAMBER    .or. &
              enefunc%forcefield == ForcefieldGROAMBER .or. &
              enefunc%forcefield == ForcefieldGROMARTINI) then

      call compute_energy_nonbond14_table_gro_amber( &
                                    domain, enefunc, &
                                    force, virial, eelec, evdw)

    end if

    return

  end subroutine compute_energy_nonbond14_table

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_table
  !> @brief        calculate nonbonded energy with lookup table
  !! @authors      JJ 
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] virial   : virial term of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : van der Waals energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_table(domain, enefunc, pairlist, &
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
    real(wp)                  :: R, h00, h01, h10, h11
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

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), trans2(:,:,:), charge(:,:)
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


    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    trans2          => domain%translated
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
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

    ! calculate energy and gradient
    !

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, ix, num_nb15, rtmp, qtmp, ini_nb15, fin_nb15, k, iy,  &
    !$omp         rij2, L, R, term_lj12, term_lj6, term_lj, grad_coef, work,   &
    !$omp         term_elec, viri_local, iwater, ij, j, trans_x, trans_y,      &
    !$omp         trans_z, iix, list, num_count, j_list, dij, dij_list,        &
    !$omp         rij2_list, force_local, force_localj, lj6, lj12, L1,         &
    !$omp         elec_temp,evdw_temp, h00, h01, h10, h11)             
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

        force_local(1:3) = 0.0_wp
!ocl norecurrence(force)
        do k = ini_nb15, fin_nb15

          iy = nb15_calc_list1(k,i)

          ! compute distance
          !
          dij(1) = rtmp(1) - trans2(iy,1,i)
          dij(2) = rtmp(2) - trans2(iy,2,i)
          dij(3) = rtmp(3) - trans2(iy,3,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2   = density * rij2
          lj6    = nonb_lj6(atmcls(ix,i),atmcls(iy,i))
          lj12   = nonb_lj12(atmcls(ix,i),atmcls(iy,i))

          L    = int(rij2)
          R    = rij2 - L
          h00  = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
          h10  = R*(1.0_wp-R)*(1.0_wp-R)
          h01  = R*R*(3.0_wp-2.0_wp*R)
          h11  = R*R*(R-1.0_wp)

          L1 = 6*L - 5
          term_lj12 = table_ene(L1)*h00 + table_ene(L1+1)*h10
          term_lj6  = table_ene(L1+2)*h00 + table_ene(L1+3)*h10
          term_elec = table_ene(L1+4)*h00 + table_ene(L1+5)*h10
          term_lj12 = term_lj12 + table_ene(L1+6)*h01 + table_ene(L1+7)*h11
          term_lj6  = term_lj6  + table_ene(L1+8)*h01 + table_ene(L1+9)*h11
          term_elec = term_elec + table_ene(L1+10)*h01 + table_ene(L1+11)*h11
          evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
          elec_temp = elec_temp + qtmp*charge(iy,i)*term_elec

          term_lj12 = table_grad(L1)*h00 + table_grad(L1+1)*h10
          term_lj6  = table_grad(L1+2)*h00 + table_grad(L1+3)*h10
          term_elec = table_grad(L1+4)*h00 + table_grad(L1+5)*h10
          term_lj12 = term_lj12 + table_grad(L1+6)*h01 + table_grad(L1+7)*h11
          term_lj6  = term_lj6 + table_grad(L1+8)*h01 + table_grad(L1+9)*h11
          term_elec = term_elec + table_grad(L1+10)*h01 + table_grad(L1+11)*h11
          grad_coef = term_lj12*lj12 - term_lj6*lj6                             &
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

          ! virial
          !
          viri_local(1,1) = viri_local(1,1) + dij(1)*work(1)
          viri_local(2,2) = viri_local(2,2) + dij(2)*work(2)
          viri_local(3,3) = viri_local(3,3) + dij(3)*work(3)
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

      trans_x = real(cell_move(1,j,i),wp) * system_size(1)
      trans_y = real(cell_move(2,j,i),wp) * system_size(2)
      trans_z = real(cell_move(3,j,i),wp) * system_size(3)

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

        force_local(1:3) = 0.0_wp
        do k = ini_nb15, fin_nb15
          iy   = nb15_calc_list(k,ij)
          dij(1) = rtmp(1) - trans2(iy,1,j) + trans_x
          dij(2) = rtmp(2) - trans2(iy,2,j) + trans_y
          dij(3) = rtmp(3) - trans2(iy,3,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2  = density * rij2
          lj6  = nonb_lj6 (atmcls(ix,i),atmcls(iy,j))
          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,j))

          L    = int(rij2)
          R    = rij2 - L
          h00  = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
          h10  = R*(1.0_wp-R)*(1.0_wp-R)
          h01  = R*R*(3.0_wp-2.0_wp*R)
          h11  = R*R*(R-1.0_wp)

          L1 = 6*L - 5
          term_lj12 = table_ene(L1)*h00 + table_ene(L1+1)*h10
          term_lj6  = table_ene(L1+2)*h00 + table_ene(L1+3)*h10
          term_elec = table_ene(L1+4)*h00 + table_ene(L1+5)*h10
          term_lj12 = term_lj12 + table_ene(L1+6)*h01 + table_ene(L1+7)*h11
          term_lj6  = term_lj6 + table_ene(L1+8)*h01 + table_ene(L1+9)*h11
          term_elec = term_elec + table_ene(L1+10)*h01 + table_ene(L1+11)*h11
          evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
          elec_temp = elec_temp + qtmp*charge(iy,j)*term_elec

          term_lj12 = table_grad(L1)*h00 + table_grad(L1+1)*h10
          term_lj6  = table_grad(L1+2)*h00 + table_grad(L1+3)*h10
          term_elec = table_grad(L1+4)*h00 + table_grad(L1+5)*h10
          term_lj12 = term_lj12 + table_grad(L1+6)*h01 + table_grad(L1+7)*h11
          term_lj6  = term_lj6 + table_grad(L1+8)*h01 + table_grad(L1+9)*h11
          term_elec = term_elec + table_grad(L1+10)*h01 + table_grad(L1+11)*h11
          grad_coef = term_lj12*lj12 - term_lj6*lj6                            &
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
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp
        do k = 1, 3
          virial(k,k,id+1) = virial(k,k,id+1) - viri_local(k,k)
        end do

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond_table

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_table
  !> @brief        calculate nonbonded force without solvents with lookup table
  !! @authors      JJ 
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] virial   : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_table(domain, enefunc, pairlist, &
                                         force, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(3,3,nthread)

    ! local variables
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: R, h00, h01, h10, h11, lj12, lj6
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

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), trans2(:,:,:), charge(:,:)
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


    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    trans2          => domain%translated
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size

    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
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

    ! calculate energy and gradient
    !

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, ix, num_nb15, rtmp, qtmp, ini_nb15, fin_nb15, k, iy,  &
    !$omp         rij2, L, R, term_lj12, term_lj6, grad_coef, work, term_elec, &
    !$omp         viri_local, iwater, ij, j, trans_x, trans_y, trans_z,        &
    !$omp         iix, list, num_count, j_list, force_local, force_local_iy,   &
    !$omp         lj12, lj6, L1, dij, dij_list, iatmcls, jatmcls, h00, h01,    &
    !$omp         h10, h11)
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
!ocl nosimd
        do k = ini_nb15, fin_nb15

          iy = nb15_calc_list1(k,i)
          jatmcls = atmcls(iy,i)
          dij(1) = rtmp(1) - trans2(iy,1,i)
          dij(2) = rtmp(2) - trans2(iy,2,i)
          dij(3) = rtmp(3) - trans2(iy,3,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2  = min(cutoff2, rij2)
          rij2  = density * rij2
          lj12 = nonb_lj12(iatmcls,jatmcls)
          lj6  = nonb_lj6 (iatmcls,jatmcls)


          L  = int(rij2)
          R  = rij2 - L
          h00   = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
          h10   = R*(1.0_wp-R)*(1.0_wp-R)
          h01   = R*R*(3.0_wp-2.0_wp*R)
          h11   = R*R*(R-1.0_wp)

          L1 = 6*L - 5
          term_lj12 = table_grad(L1)*h00 + table_grad(L1+1)*h10
          term_lj6  = table_grad(L1+2)*h00 + table_grad(L1+3)*h10
          term_elec = table_grad(L1+4)*h00 + table_grad(L1+5)*h10
          term_lj12 = term_lj12 + table_grad(L1+6)*h01 + table_grad(L1+7)*h11
          term_lj6  = term_lj6 + table_grad(L1+8)*h01 + table_grad(L1+9)*h11
          term_elec = term_elec + table_grad(L1+10)*h01 + table_grad(L1+11)*h11
          grad_coef = term_lj12*lj12 - term_lj6*lj6                            &
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

      end do

    end do

    ! interaction between different cells
    !
    do ij = id+1, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      num_nb15 = 0

      trans_x = real(cell_move(1,j,i),wp) * system_size(1)
      trans_y = real(cell_move(2,j,i),wp) * system_size(2)
      trans_z = real(cell_move(3,j,i),wp) * system_size(3)

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
!ocl nosimd
        do k = ini_nb15, fin_nb15
          iy = nb15_calc_list(k,ij)

          dij(1) = rtmp(1) - trans2(iy,1,j) + trans_x
          dij(2) = rtmp(2) - trans2(iy,2,j) + trans_y
          dij(3) = rtmp(3) - trans2(iy,3,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2 = density * rij2 

          jatmcls = atmcls(iy,j)
          lj12 = nonb_lj12(iatmcls,jatmcls)
          lj6  = nonb_lj6 (iatmcls,jatmcls)

          L  = int(rij2)
          R  = rij2 - L
          h00   = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
          h10   = R*(1.0_wp-R)*(1.0_wp-R)
          h01   = R*R*(3.0_wp-2.0_wp*R)
          h11   = R*R*(R-1.0_wp)

          L1 = 6*L - 5
          term_lj12 = table_grad(L1)*h00 + table_grad(L1+1)*h10
          term_lj6  = table_grad(L1+2)*h00 + table_grad(L1+3)*h10
          term_elec = table_grad(L1+4)*h00 + table_grad(L1+5)*h10
          term_lj12 = term_lj12 + table_grad(L1+6)*h01 + table_grad(L1+7)*h11
          term_lj6  = term_lj6 + table_grad(L1+8)*h01 + table_grad(L1+9)*h11
          term_elec = term_elec + table_grad(L1+10)*h01 + table_grad(L1+11)*h11
          grad_coef = term_lj12*lj12 - term_lj6*lj6                            &
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

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_force_nonbond_table

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_table_charmm
  !> @brief        calculate nonbonded14 energy with PME and lookup table
  !!               (linear)
  !  @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_table_charmm(domain, enefunc, &
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
    real(wp)                 :: R, h00, h10, h01, h11, table(12)
    real(wp)                 :: term_lj12, term_lj6, term_elec
    real(wp)                 :: cutoff2, grad_coef, lj12, lj6
    real(wp)                 :: work(1:3), elec_temp, evdw_temp
    real(wp)                 :: viri_local(1:3,1:3)
    integer                  :: i, ix, iy, j, k, ij, m, L
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: ncell_local, id, omp_get_thread_num

    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: nb14_lj12(:,:), nb14_lj6(:,:)
    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: atmcls(:,:), natom(:)
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: num_nb14_calc(:), nb14_calc_list(:,:,:)


    coord           => domain%coord
    charge          => domain%charge
    cell_pair       => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nb14_lj12       => enefunc%nb14_lj12
    nb14_lj6        => enefunc%nb14_lj6
    num_nb14_calc   => enefunc%num_nb14_calc
    nb14_calc_list  => enefunc%nb14_calc_list

    cutoff2         =  cutoff * cutoff
    ncell_local     =  domain%num_cell_local


    !$omp parallel default(shared)                                             &
    !$omp private(id, i, j, k, ij, ix, iy, m, ini_nb14, fin_nb14, num_nb14,    &
    !$omp         dij, rij2, L, R, term_lj12, term_lj6, term_elec, grad_coef,  &
    !$omp         work, viri_local, table, lj12, lj6, h00, h10, h01, h11,      &
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

      viri_local(1:3,1:3) = 0.0_wp
      elec_temp = 0.0_wp
      evdw_temp = 0.0_wp

      do k = 1, num_nb14_calc(ij)

        i  = nb14_calc_list(1,k,ij)
        j  = nb14_calc_list(2,k,ij)
        ix = nb14_calc_list(3,k,ij)
        iy = nb14_calc_list(4,k,ij)

        ! compute distance
        !
        dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,j)
        rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        ! energy and gradient
        !
        rij2  = density * rij2
        L     = int(rij2)
        R     = rij2 - L
        lj6   = nb14_lj6 (atmcls(ix,i),atmcls(iy,j))
        lj12  = nb14_lj12(atmcls(ix,i),atmcls(iy,j))
        h00   = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
        h10   = R*(1.0_wp-R)*(1.0_wp-R)
        h01   = R*R*(3.0_wp-2.0_wp*R)
        h11   = R*R*(R-1.0_wp)

        table(1:12) = table_ene(6*L-5:6*L+6)
        term_lj12 = table(1)*h00 + table(2)*h10 + table(7)*h01 + table(8)*h11
        term_lj6  = table(3)*h00 + table(4)*h10 + table(9)*h01 + table(10)*h11
        term_elec = table(5)*h00 + table(6)*h10 + table(11)*h01 +table(12)*h11
        evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
        elec_temp = elec_temp + charge(ix,i)*charge(iy,j)*term_elec

        table(1:12) = table_grad(6*L-5:6*L+6)
        term_lj12 = table(1)*h00 + table(2)*h10 + table(7)*h01 + table(8)*h11
        term_lj6  = table(3)*h00 + table(4)*h10 + table(9)*h01 + table(10)*h11
        term_elec = table(5)*h00 + table(6)*h10 + table(11)*h01 +table(12)*h11
        grad_coef = term_lj12*lj12 - term_lj6*lj6 + &
                    charge(ix,i)*charge(iy,j)*term_elec

        work(1:3) = grad_coef*dij(1:3)

        force(ix,1:3,i,id+1) = force(ix,1:3,i,id+1) - work(1:3)
        force(iy,1:3,j,id+1) = force(iy,1:3,j,id+1) + work(1:3)

        ! virial
        !
        do m = 1, 3
          viri_local(1:3,m) = viri_local(1:3,m) + dij(1:3)*work(m)
        end do
      end do

      eelec(id+1) = eelec(id+1) + elec_temp
      evdw(id+1) = evdw(id+1) + evdw_temp
      do m = 1, 3
        virial(m,m,id+1) = virial(m,m,id+1) - viri_local(m,m)
      end do
      
    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond14_table_charmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_table_gro_amber
  !> @brief        calculate nonbonded14 energy with PME and lookup table
  !!               (linear)
  !  @authors      JJ, CK
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_table_gro_amber(domain, enefunc, &
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
    real(wp)                 :: R, h00, h10, h01, h11, table(12)
    real(wp)                 :: term_lj12, term_lj6, term_elec
    real(wp)                 :: cutoff2, grad_coef, lj12, lj6
    real(wp)                 :: work(1:3), elec_temp, evdw_temp
    real(wp)                 :: viri_local(1:3,1:3)
    real(wp)                 :: lj_scale, qq_scale, cc
    integer                  :: i, ix, iy, j, k, ij, m, L
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: ncell_local, id, omp_get_thread_num

    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: nb14_lj12(:,:), nb14_lj6(:,:)
    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: atmcls(:,:)
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: natom(:)
    integer,         pointer :: num_nb14_calc(:), nb14_calc_list(:,:,:)


    coord           => domain%coord
    charge          => domain%charge
    cell_pair       => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nb14_lj12       => enefunc%nb14_lj12
    nb14_lj6        => enefunc%nb14_lj6
    num_nb14_calc   => enefunc%num_nb14_calc
    nb14_calc_list  => enefunc%nb14_calc_list

    cutoff2         =  cutoff * cutoff
    ncell_local     =  domain%num_cell_local


    !$omp parallel default(shared)                                             &
    !$omp private(id, i, j, k, ij, ix, iy, m, ini_nb14, fin_nb14, num_nb14,    &
    !$omp         dij, rij2, L, R, term_lj12, term_lj6, term_elec, grad_coef,  &
    !$omp         work, viri_local, table, lj12, lj6, h00, h10, h01, h11,      &
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

      elec_temp = 0.0_wp
      evdw_temp = 0.0_wp
      viri_local(1:3,1:3) = 0.0_wp

      do k = 1, num_nb14_calc(ij)

        i  = nb14_calc_list(1,k,ij)
        j  = nb14_calc_list(2,k,ij)
        ix = nb14_calc_list(3,k,ij)
        iy = nb14_calc_list(4,k,ij)

        ! compute distance
        !
        dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,j)
        rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        qq_scale = enefunc%nb14_qq_scale(k,ij)
        lj_scale = enefunc%nb14_lj_scale(k,ij)

        ! energy and gradient
        !
        rij2  = density * rij2
        L     = int(rij2)
        R     = rij2 - L
        lj6   = nb14_lj6 (atmcls(ix,i),atmcls(iy,j))
        lj12  = nb14_lj12(atmcls(ix,i),atmcls(iy,j))
        h00   = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
        h10   = R*(1.0_wp-R)*(1.0_wp-R)
        h01   = R*R*(3.0_wp-2.0_wp*R)
        h11   = R*R*(R-1.0_wp)
        cc    = charge(ix,i)*charge(iy,j)*qq_scale

        table(1:12) = table_ene(6*L-5:6*L+6)
        term_lj12 = table(1)*h00 + table(2)*h10 + table(7)*h01 + table(8)*h11
        term_lj6  = table(3)*h00 + table(4)*h10 + table(9)*h01 + table(10)*h11
        term_elec = table(5)*h00 + table(6)*h10 + table(11)*h01 +table(12)*h11
        evdw_temp = evdw_temp + (term_lj12*lj12 - term_lj6*lj6)*lj_scale
        elec_temp = elec_temp + cc*term_elec

        table(1:12) = table_grad(6*L-5:6*L+6)
        term_lj12 = table(1)*h00 + table(2)*h10 + table(7)*h01 + table(8)*h11
        term_lj6  = table(3)*h00 + table(4)*h10 + table(9)*h01 + table(10)*h11
        term_elec = table(5)*h00 + table(6)*h10 + table(11)*h01 +table(12)*h11
        grad_coef = (term_lj12*lj12 - term_lj6*lj6)*lj_scale + cc*term_elec

        work(1:3) = grad_coef*dij(1:3)

        force(ix,1:3,i,id+1) = force(ix,1:3,i,id+1) - work(1:3)
        force(iy,1:3,j,id+1) = force(iy,1:3,j,id+1) + work(1:3)

        ! virial
        !
        do m = 1, 3
          viri_local(1:3,m) = viri_local(1:3,m) + dij(1:3)*work(m)
        end do

      end do

      eelec(id+1) = eelec(id+1) + elec_temp
      evdw(id+1) = evdw(id+1) + evdw_temp

      do m = 1, 3
        virial(m,m,id+1) = virial(m,m,id+1) - viri_local(m,m)
      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond14_table_gro_amber

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_table_fep
  !> @brief        calculate nonbonded14 energy with PME and lookup table
  !!               (linear) for FEP
  !  @authors      HO
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_table_fep(domain, enefunc, &
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

      call compute_energy_nonbond14_table_charmm(    &
                                    domain, enefunc, &
                                    force, virial, eelec, evdw)

      if (domain%fep_use) then
        call compute_energy_nonbond14_table_charmm_fep( &
                                      domain, enefunc, &
                                      force, virial, eelec, evdw)
      end if

    ! ==> Type 10
    else  if (enefunc%forcefield == ForcefieldAMBER    .or. &
              enefunc%forcefield == ForcefieldGROAMBER .or. &
              enefunc%forcefield == ForcefieldGROMARTINI) then

      call compute_energy_nonbond14_table_gro_amber( &
                                    domain, enefunc, &
                                    force, virial, eelec, evdw)

      if (domain%fep_use) then
        call compute_energy_nonbond14_table_gro_amber_fep( &
                                      domain, enefunc, &
                                      force, virial, eelec, evdw)
      end if

    endif

    return

  end subroutine compute_energy_nonbond14_table_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_table_fep
  !> @brief        calculate nonbonded energy with lookup table for FEP
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

  subroutine compute_energy_nonbond_table_fep(domain, enefunc, pairlist, &
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
    real(wp)                  :: R, h00, h01, h10, h11
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

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), trans2(:,:,:), charge(:,:)
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

    ! FEP
    real(wp)                  :: rij2_sclj, rij2_scel
    integer                   :: fg1, fg2
    integer,          pointer :: fepgrp(:,:)
    real(wp),         pointer :: table_sclj(:,:)
    real(wp),         pointer :: table_scel(:,:)

    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    trans2          => domain%translated
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size

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

    ! calculate energy and gradient
    !

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, ix, num_nb15, rtmp, qtmp, ini_nb15, fin_nb15, k, iy,  &
    !$omp         rij2, L, R, term_lj12, term_lj6, term_lj, grad_coef, work,   &
    !$omp         term_elec, viri_local, iwater, ij, j, trans_x, trans_y,      &
    !$omp         trans_z, iix, list, num_count, j_list, dij, dij_list,        &
    !$omp         rij2_list, force_local, force_localj, lj6, lj12, L1,         &
    !$omp         elec_temp,evdw_temp, h00, h01, h10, h11,                     &
    !$omp         fg1, fg2, rij2_sclj, rij2_scel)
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

        force_local(1:3) = 0.0_wp
!ocl norecurrence(force)
        do k = ini_nb15, fin_nb15

          iy = nb15_calc_list1(k,i)

          ! compute distance
          !
          dij(1) = rtmp(1) - trans2(iy,1,i)
          dij(2) = rtmp(2) - trans2(iy,2,i)
          dij(3) = rtmp(3) - trans2(iy,3,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          ! FEP: Determine softcore
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,i)
          rij2_sclj = rij2 + table_sclj(fg1,fg2)
          rij2_scel = rij2 + table_scel(fg1,fg2)

          ! FEP: LJ with soft core
          rij2_sclj = density * rij2_sclj
          lj6    = nonb_lj6(atmcls(ix,i),atmcls(iy,i))
          lj12   = nonb_lj12(atmcls(ix,i),atmcls(iy,i))
          L    = int(rij2_sclj)
          R    = rij2_sclj - L
          h00  = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
          h10  = R*(1.0_wp-R)*(1.0_wp-R)
          h01  = R*R*(3.0_wp-2.0_wp*R)
          h11  = R*R*(R-1.0_wp)

          L1 = 6*L - 5
          term_lj12 = table_ene(L1)*h00 + table_ene(L1+1)*h10
          term_lj6  = table_ene(L1+2)*h00 + table_ene(L1+3)*h10
          term_lj12 = term_lj12 + table_ene(L1+6)*h01 + table_ene(L1+7)*h11
          term_lj6  = term_lj6  + table_ene(L1+8)*h01 + table_ene(L1+9)*h11
          evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
          term_lj12 = table_grad(L1)*h00 + table_grad(L1+1)*h10
          term_lj6  = table_grad(L1+2)*h00 + table_grad(L1+3)*h10
          term_lj12 = term_lj12 + table_grad(L1+6)*h01 + table_grad(L1+7)*h11
          term_lj6  = term_lj6 + table_grad(L1+8)*h01 + table_grad(L1+9)*h11

          ! FEP: EL with soft core
          rij2_scel = density * rij2_scel
          L    = int(rij2_scel)
          R    = rij2_scel - L
          h00  = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
          h10  = R*(1.0_wp-R)*(1.0_wp-R)
          h01  = R*R*(3.0_wp-2.0_wp*R)
          h11  = R*R*(R-1.0_wp)
          L1 = 6*L - 5
          term_elec = table_ene(L1+4)*h00 + table_ene(L1+5)*h10
          term_elec = term_elec + table_ene(L1+10)*h01 + table_ene(L1+11)*h11
          elec_temp = elec_temp + qtmp*charge(iy,i)*term_elec
          term_elec = table_grad(L1+4)*h00 + table_grad(L1+5)*h10
          term_elec = term_elec + table_grad(L1+10)*h01 + table_grad(L1+11)*h11

          grad_coef = term_lj12*lj12 - term_lj6*lj6                             &
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

          ! virial
          !
          viri_local(1,1) = viri_local(1,1) + dij(1)*work(1)
          viri_local(2,2) = viri_local(2,2) + dij(2)*work(2)
          viri_local(3,3) = viri_local(3,3) + dij(3)*work(3)
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

      trans_x = real(cell_move(1,j,i),wp) * system_size(1)
      trans_y = real(cell_move(2,j,i),wp) * system_size(2)
      trans_z = real(cell_move(3,j,i),wp) * system_size(3)

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

        force_local(1:3) = 0.0_wp
        do k = ini_nb15, fin_nb15
          iy   = nb15_calc_list(k,ij)
          dij(1) = rtmp(1) - trans2(iy,1,j) + trans_x
          dij(2) = rtmp(2) - trans2(iy,2,j) + trans_y
          dij(3) = rtmp(3) - trans2(iy,3,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          ! FEP: Determine softcore
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,j)
          rij2_sclj = rij2 + table_sclj(fg1,fg2)
          rij2_scel = rij2 + table_scel(fg1,fg2)

          ! FEP: LJ with soft core
          rij2_sclj  = density * rij2_sclj
          lj6  = nonb_lj6 (atmcls(ix,i),atmcls(iy,j))
          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,j))
          L    = int(rij2_sclj)
          R    = rij2_sclj - L
          h00  = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
          h10  = R*(1.0_wp-R)*(1.0_wp-R)
          h01  = R*R*(3.0_wp-2.0_wp*R)
          h11  = R*R*(R-1.0_wp)

          L1 = 6*L - 5
          term_lj12 = table_ene(L1)*h00 + table_ene(L1+1)*h10
          term_lj6  = table_ene(L1+2)*h00 + table_ene(L1+3)*h10
          term_lj12 = term_lj12 + table_ene(L1+6)*h01 + table_ene(L1+7)*h11
          term_lj6  = term_lj6 + table_ene(L1+8)*h01 + table_ene(L1+9)*h11
          evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
          term_lj12 = table_grad(L1)*h00 + table_grad(L1+1)*h10
          term_lj6  = table_grad(L1+2)*h00 + table_grad(L1+3)*h10
          term_lj12 = term_lj12 + table_grad(L1+6)*h01 + table_grad(L1+7)*h11
          term_lj6  = term_lj6 + table_grad(L1+8)*h01 + table_grad(L1+9)*h11

          ! FEP: elec with soft core
          rij2_scel = density * rij2_scel
          L    = int(rij2_scel)
          R    = rij2_scel - L
          h00  = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
          h10  = R*(1.0_wp-R)*(1.0_wp-R)
          h01  = R*R*(3.0_wp-2.0_wp*R)
          h11  = R*R*(R-1.0_wp)
          L1 = 6*L - 5
          term_elec = table_ene(L1+4)*h00 + table_ene(L1+5)*h10
          term_elec = term_elec + table_ene(L1+10)*h01 + table_ene(L1+11)*h11
          elec_temp = elec_temp + qtmp*charge(iy,j)*term_elec
          term_elec = table_grad(L1+4)*h00 + table_grad(L1+5)*h10
          term_elec = term_elec + table_grad(L1+10)*h01 + table_grad(L1+11)*h11

          grad_coef = term_lj12*lj12 - term_lj6*lj6                            &
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
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp
        do k = 1, 3
          virial(k,k,id+1) = virial(k,k,id+1) - viri_local(k,k)
        end do

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond_table_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_table
  !> @brief        calculate nonbonded force without solvents with lookup table
  !! @authors      JJ 
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] virial   : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_table_fep(domain, enefunc, pairlist, &
                                                force, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(3,3,nthread)

    ! local variables
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: R, h00, h01, h10, h11, lj12, lj6
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

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), trans2(:,:,:), charge(:,:)
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

    ! FEP
    real(wp)                  :: rij2_sclj, rij2_scel
    integer                   :: fg1, fg2
    integer,          pointer :: fepgrp(:,:)
    real(wp),         pointer :: table_sclj(:,:)
    real(wp),         pointer :: table_scel(:,:)

    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    trans2          => domain%translated
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size

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

    ! calculate energy and gradient
    !

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, ix, num_nb15, rtmp, qtmp, ini_nb15, fin_nb15, k, iy,  &
    !$omp         rij2, L, R, term_lj12, term_lj6, grad_coef, work, term_elec, &
    !$omp         viri_local, iwater, ij, j, trans_x, trans_y, trans_z,        &
    !$omp         iix, list, num_count, j_list, force_local, force_local_iy,   &
    !$omp         lj12, lj6, L1, dij, dij_list, iatmcls, jatmcls, h00, h01,    &
    !$omp         h10, h11, rij2_sclj, rij2_scel, fg1, fg2)
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

        ! FEP: flag for atom ix
        fg1 = fepgrp(ix,i)

!ocl norecurrence(force)
!ocl nosimd
        do k = ini_nb15, fin_nb15

          iy = nb15_calc_list1(k,i)
          jatmcls = atmcls(iy,i)
          dij(1) = rtmp(1) - trans2(iy,1,i)
          dij(2) = rtmp(2) - trans2(iy,2,i)
          dij(3) = rtmp(3) - trans2(iy,3,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2  = min(cutoff2, rij2)

          ! FEP: Determine softcore
          fg2 = fepgrp(iy,i)
          rij2_sclj = rij2 + table_sclj(fg1,fg2)
          rij2_scel = rij2 + table_scel(fg1,fg2)

          ! FEP: LJ with soft-core
          rij2_sclj = density * rij2_sclj
          lj12 = nonb_lj12(iatmcls,jatmcls)
          lj6  = nonb_lj6 (iatmcls,jatmcls)
          L  = int(rij2_sclj)
          R  = rij2_sclj - L
          h00   = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
          h10   = R*(1.0_wp-R)*(1.0_wp-R)
          h01   = R*R*(3.0_wp-2.0_wp*R)
          h11   = R*R*(R-1.0_wp)

          L1 = 6*L - 5
          term_lj12 = table_grad(L1)*h00 + table_grad(L1+1)*h10
          term_lj6  = table_grad(L1+2)*h00 + table_grad(L1+3)*h10
          term_lj12 = term_lj12 + table_grad(L1+6)*h01 + table_grad(L1+7)*h11
          term_lj6  = term_lj6 + table_grad(L1+8)*h01 + table_grad(L1+9)*h11

          ! FEP: LJ with soft-core
          rij2_scel = density * rij2_scel
          L  = int(rij2_scel)
          R  = rij2_scel - L
          h00   = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
          h10   = R*(1.0_wp-R)*(1.0_wp-R)
          h01   = R*R*(3.0_wp-2.0_wp*R)
          h11   = R*R*(R-1.0_wp)
          L1 = 6*L - 5
          term_elec = table_grad(L1+4)*h00 + table_grad(L1+5)*h10
          term_elec = term_elec + table_grad(L1+10)*h01 + table_grad(L1+11)*h11

          grad_coef = term_lj12*lj12 - term_lj6*lj6                            &
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

      end do

    end do

    ! interaction between different cells
    !
    do ij = id+1, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      num_nb15 = 0

      trans_x = real(cell_move(1,j,i),wp) * system_size(1)
      trans_y = real(cell_move(2,j,i),wp) * system_size(2)
      trans_z = real(cell_move(3,j,i),wp) * system_size(3)

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

        ! FEP: flag for atom ix
        fg1 = fepgrp(ix,i)

!ocl norecurrence(force)
!ocl nosimd
        do k = ini_nb15, fin_nb15
          iy = nb15_calc_list(k,ij)

          dij(1) = rtmp(1) - trans2(iy,1,j) + trans_x
          dij(2) = rtmp(2) - trans2(iy,2,j) + trans_y
          dij(3) = rtmp(3) - trans2(iy,3,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          ! FEP: Determine softcore
          fg2 = fepgrp(iy,j)
          rij2_sclj = rij2 + table_sclj(fg1,fg2)
          rij2_scel = rij2 + table_scel(fg1,fg2)

          ! FEP: LJ with soft-core
          rij2_sclj = density * rij2_sclj
          jatmcls = atmcls(iy,j)
          lj12 = nonb_lj12(iatmcls,jatmcls)
          lj6  = nonb_lj6 (iatmcls,jatmcls)

          L  = int(rij2_sclj)
          R  = rij2_sclj - L
          h00   = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
          h10   = R*(1.0_wp-R)*(1.0_wp-R)
          h01   = R*R*(3.0_wp-2.0_wp*R)
          h11   = R*R*(R-1.0_wp)

          L1 = 6*L - 5
          term_lj12 = table_grad(L1)*h00 + table_grad(L1+1)*h10
          term_lj6  = table_grad(L1+2)*h00 + table_grad(L1+3)*h10
          term_lj12 = term_lj12 + table_grad(L1+6)*h01 + table_grad(L1+7)*h11
          term_lj6  = term_lj6 + table_grad(L1+8)*h01 + table_grad(L1+9)*h11

          ! FEP: elec with soft-core
          rij2_scel = density * rij2_scel
          L  = int(rij2_scel)
          R  = rij2_scel - L
          h00   = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
          h10   = R*(1.0_wp-R)*(1.0_wp-R)
          h01   = R*R*(3.0_wp-2.0_wp*R)
          h11   = R*R*(R-1.0_wp)
          L1 = 6*L - 5
          term_elec = table_grad(L1+4)*h00 + table_grad(L1+5)*h10
          term_elec = term_elec + table_grad(L1+10)*h01 + table_grad(L1+11)*h11

          grad_coef = term_lj12*lj12 - term_lj6*lj6                            &
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

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_force_nonbond_table_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_table_charmm_fep
  !> @brief        calculate nonbonded14 energy with PME and lookup table
  !!               (linear) for FEP
  !  @authors      HO
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_table_charmm_fep(domain, enefunc, &
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
    real(wp)                 :: R, h00, h10, h01, h11, table(12)
    real(wp)                 :: term_lj12, term_lj6, term_elec
    real(wp)                 :: cutoff2, grad_coef, lj12, lj6
    real(wp)                 :: work(1:3), elec_temp, evdw_temp
    real(wp)                 :: viri_local(1:3,1:3)
    integer                  :: i, ix, iy, j, k, ij, m, L
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: ncell_local, id, omp_get_thread_num

    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: nb14_lj12(:,:), nb14_lj6(:,:)
    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: atmcls(:,:), natom(:)
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: num_nb14_calc(:), nb14_calc_list(:,:,:)

    ! FEP
    real(wp)                  :: rij2_sclj, rij2_scel
    integer                   :: fg1, fg2
    integer,          pointer :: fepgrp(:,:)
    real(wp),         pointer :: table_sclj(:,:)
    real(wp),         pointer :: table_scel(:,:)

    coord           => domain%coord
    charge          => domain%charge
    cell_pair       => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
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

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, j, k, ij, ix, iy, m, ini_nb14, fin_nb14, num_nb14,    &
    !$omp         dij, rij2, L, R, term_lj12, term_lj6, term_elec, grad_coef,  &
    !$omp         work, viri_local, table, lj12, lj6, h00, h10, h01, h11,      &
    !$omp         elec_temp, evdw_temp, rij2_sclj, rij2_scel,                  &
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

      viri_local(1:3,1:3) = 0.0_wp
      elec_temp = 0.0_wp
      evdw_temp = 0.0_wp

      do k = 1, num_nb14_calc(ij)

        i  = nb14_calc_list(1,k,ij)
        j  = nb14_calc_list(2,k,ij)
        ix = nb14_calc_list(3,k,ij)
        iy = nb14_calc_list(4,k,ij)

        ! compute distance
        !
        dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,j)
        rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        ! FEP: Determine softcore
        fg1 = fepgrp(ix,i)
        fg2 = fepgrp(iy,j)
        rij2_sclj = rij2 + table_sclj(fg1,fg2)
        rij2_scel = rij2 + table_scel(fg1,fg2)

        ! FEP: LJ with soft core for reference
        rij2_sclj = density * rij2_sclj
        L     = int(rij2_sclj)
        R     = rij2_sclj - L
        lj6   = nb14_lj6 (atmcls(ix,i),atmcls(iy,j))
        lj12  = nb14_lj12(atmcls(ix,i),atmcls(iy,j))
        h00   = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
        h10   = R*(1.0_wp-R)*(1.0_wp-R)
        h01   = R*R*(3.0_wp-2.0_wp*R)
        h11   = R*R*(R-1.0_wp)

        table(1:12) = table_ene(6*L-5:6*L+6)
        term_lj12 = table(1)*h00 + table(2)*h10 + table(7)*h01 + table(8)*h11
        term_lj6  = table(3)*h00 + table(4)*h10 + table(9)*h01 + table(10)*h11
        evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
        table(1:12) = table_grad(6*L-5:6*L+6)
        term_lj12 = table(1)*h00 + table(2)*h10 + table(7)*h01 + table(8)*h11
        term_lj6  = table(3)*h00 + table(4)*h10 + table(9)*h01 + table(10)*h11

        ! FEP: elec with soft core for reference
        rij2_scel = density * rij2_scel
        L     = int(rij2_scel)
        R     = rij2_scel - L
        h00   = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
        h10   = R*(1.0_wp-R)*(1.0_wp-R)
        h01   = R*R*(3.0_wp-2.0_wp*R)
        h11   = R*R*(R-1.0_wp)
        table(1:12) = table_ene(6*L-5:6*L+6)
        term_elec = table(5)*h00 + table(6)*h10 + table(11)*h01 +table(12)*h11
        elec_temp = elec_temp + charge(ix,i)*charge(iy,j)*term_elec
        table(1:12) = table_grad(6*L-5:6*L+6)
        term_elec = table(5)*h00 + table(6)*h10 + table(11)*h01 +table(12)*h11

        grad_coef = term_lj12*lj12 - term_lj6*lj6 + &
                    charge(ix,i)*charge(iy,j)*term_elec

        work(1:3) = grad_coef*dij(1:3)

        force(ix,1:3,i,id+1) = force(ix,1:3,i,id+1) - work(1:3)
        force(iy,1:3,j,id+1) = force(iy,1:3,j,id+1) + work(1:3)

        ! virial
        !
        do m = 1, 3
          viri_local(1:3,m) = viri_local(1:3,m) + dij(1:3)*work(m)
        end do
      end do

      eelec(id+1) = eelec(id+1) + elec_temp
      evdw(id+1) = evdw(id+1) + evdw_temp
      do m = 1, 3
        virial(m,m,id+1) = virial(m,m,id+1) - viri_local(m,m)
      end do
      
    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond14_table_charmm_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_table_gro_amber_fep
  !> @brief        calculate nonbonded14 energy with PME and lookup table
  !!               (linear) for FEP
  !  @authors      HO
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_table_gro_amber_fep(domain, enefunc, &
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
    real(wp)                 :: R, h00, h10, h01, h11, table(12)
    real(wp)                 :: term_lj12, term_lj6, term_elec
    real(wp)                 :: cutoff2, grad_coef, lj12, lj6
    real(wp)                 :: work(1:3), elec_temp, evdw_temp
    real(wp)                 :: viri_local(1:3,1:3)
    real(wp)                 :: lj_scale, qq_scale, cc
    integer                  :: i, ix, iy, j, k, ij, m, L
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: ncell_local, id, omp_get_thread_num

    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: nb14_lj12(:,:), nb14_lj6(:,:)
    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: atmcls(:,:)
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: natom(:)
    integer,         pointer :: num_nb14_calc(:), nb14_calc_list(:,:,:)

    ! FEP
    real(wp)                  :: rij2_sclj, rij2_scel
    integer                   :: fg1, fg2
    integer,          pointer :: fepgrp(:,:)
    real(wp),         pointer :: table_sclj(:,:)
    real(wp),         pointer :: table_scel(:,:)

    coord           => domain%coord
    charge          => domain%charge
    cell_pair       => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nb14_lj12       => enefunc%nb14_lj12
    nb14_lj6        => enefunc%nb14_lj6

    cutoff2         =  cutoff * cutoff
    ncell_local     =  domain%num_cell_local

    ! FEP
    num_nb14_calc  => enefunc%num_nb14_calc_fep
    nb14_calc_list => enefunc%nb14_calc_list_fep
    fepgrp         => domain%fepgrp
    table_sclj     => enefunc%table_sclj
    table_scel     => enefunc%table_scel

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, j, k, ij, ix, iy, m, ini_nb14, fin_nb14, num_nb14,    &
    !$omp         dij, rij2, L, R, term_lj12, term_lj6, term_elec, grad_coef,  &
    !$omp         work, viri_local, table, lj12, lj6, h00, h10, h01, h11,      &
    !$omp         qq_scale, lj_scale, cc, elec_temp, evdw_temp,                &
    !$omp         fg1, fg2, rij2_sclj, rij2_scel)
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
      viri_local(1:3,1:3) = 0.0_wp

      do k = 1, num_nb14_calc(ij)

        i  = nb14_calc_list(1,k,ij)
        j  = nb14_calc_list(2,k,ij)
        ix = nb14_calc_list(3,k,ij)
        iy = nb14_calc_list(4,k,ij)

        ! compute distance
        !
        dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,j)
        rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        qq_scale = enefunc%nb14_qq_scale_fep(k,ij)
        lj_scale = enefunc%nb14_lj_scale_fep(k,ij)

        ! FEP: Determine softcore
        fg1 = fepgrp(ix,i)
        fg2 = fepgrp(iy,j)
        rij2_sclj = rij2 + table_sclj(fg1,fg2)
        rij2_scel = rij2 + table_scel(fg1,fg2)

        ! FEP: LJ with soft core for reference
        rij2_sclj = density * rij2_sclj
        L     = int(rij2_sclj)
        R     = rij2_sclj - L
        h00   = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
        h10   = R*(1.0_wp-R)*(1.0_wp-R)
        h01   = R*R*(3.0_wp-2.0_wp*R)
        h11   = R*R*(R-1.0_wp)
        lj6   = nb14_lj6 (atmcls(ix,i),atmcls(iy,j))
        lj12  = nb14_lj12(atmcls(ix,i),atmcls(iy,j))
        table(1:12) = table_ene(6*L-5:6*L+6)
        term_lj12 = table(1)*h00 + table(2)*h10 + table(7)*h01 + table(8)*h11
        term_lj6  = table(3)*h00 + table(4)*h10 + table(9)*h01 + table(10)*h11
        evdw_temp = evdw_temp + (term_lj12*lj12 - term_lj6*lj6)*lj_scale
        table(1:12) = table_grad(6*L-5:6*L+6)
        term_lj12 = table(1)*h00 + table(2)*h10 + table(7)*h01 + table(8)*h11
        term_lj6  = table(3)*h00 + table(4)*h10 + table(9)*h01 + table(10)*h11

        ! FEP: elec with soft core for reference
        rij2_scel = density * rij2_scel
        L     = int(rij2_scel)
        R     = rij2_scel - L
        h00   = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
        h10   = R*(1.0_wp-R)*(1.0_wp-R)
        h01   = R*R*(3.0_wp-2.0_wp*R)
        h11   = R*R*(R-1.0_wp)
        cc    = charge(ix,i)*charge(iy,j)*qq_scale
        table(1:12) = table_ene(6*L-5:6*L+6)
        term_elec = table(5)*h00 + table(6)*h10 + table(11)*h01 +table(12)*h11
        elec_temp = elec_temp + cc*term_elec
        table(1:12) = table_grad(6*L-5:6*L+6)
        term_elec = table(5)*h00 + table(6)*h10 + table(11)*h01 +table(12)*h11

        grad_coef = (term_lj12*lj12 - term_lj6*lj6)*lj_scale + cc*term_elec

        work(1:3) = grad_coef*dij(1:3)

        force(ix,1:3,i,id+1) = force(ix,1:3,i,id+1) - work(1:3)
        force(iy,1:3,j,id+1) = force(iy,1:3,j,id+1) + work(1:3)

        ! virial
        !
        do m = 1, 3
          viri_local(1:3,m) = viri_local(1:3,m) + dij(1:3)*work(m)
        end do

      end do

      eelec(id+1) = eelec(id+1) + elec_temp
      evdw(id+1) = evdw(id+1) + evdw_temp

      do m = 1, 3
        virial(m,m,id+1) = virial(m,m,id+1) - viri_local(m,m)
      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond14_table_gro_amber_fep

end module sp_energy_table_cubic_mod
