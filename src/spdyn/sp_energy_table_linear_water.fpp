!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_table_linear_water_mod
!> @brief   calculate nonbonded energy with table and with linear interpolation
!! @authors Jaewoon Jung(JJ)
!  
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_table_linear_water_mod

  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  !
  public  :: compute_energy_nonbond_table_water_linear
  public  :: compute_force_nonbond_table_water_linear

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_table_water_linear
  !> @brief        calculate nonbonded energy among water molecules
  !! @authors      JJ 
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

  subroutine compute_energy_nonbond_table_water_linear(domain, enefunc, &
                                           pairlist, coord_pbc, force,  &
                                           virial, eelec, evdw)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(wip),                intent(inout) :: virial(:,:)
    real(wip),                intent(inout) :: eelec(nthread)
    real(wip),                intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: R
    real(wp)                  :: cutoff2, cutoff_small
    real(wp)                  :: term_lj, grad_coef, term_elec
    real(wp)                  :: work(1:3)
    real(wp)                  :: trans_x, trans_y, trans_z
    real(wp)                  :: rtmp(3,3)
    real(wp)                  :: force_local(3,3)
    real(wp)                  :: elec_temp, evdw_temp
    real(wp)                  :: force_localj(3,3,3,MaxWater)
    integer                   :: i, iy, j, k, ij, iix, iwater, jwater
    integer                   :: num_nb15, ini_nb15, fin_nb15
    integer                   :: id, omp_get_thread_num
    integer                   :: i1, j1, L, L1
    integer                   :: num_count(3,3), c
    integer                   :: small_count, j_small(MaxWater)
    integer                   :: ilist(3), jlist(3)
    integer                   :: ncell_local
    integer                   :: water_type(3,3), type
    integer                   :: check_virial

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), system_size(:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_de_WW(:,:)
    real(wp),         pointer :: table_ene_WW(:,:)
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: natom(:)
    integer(int2),    pointer :: cell_pairlist(:,:)
    integer,          pointer :: nwater(:), nsolute(:)
    integer,          pointer :: nb15_cellww(:), nb15_listww(:,:)
    integer,          pointer :: num_nb15_calcww(:,:), nb15_calc_water(:,:)
    integer,          pointer :: virial_check(:,:)


    cell_pairlist   => domain%cell_pairlist1
    natom           => domain%num_atom
    nsolute         => domain%num_solute
    nwater          => domain%num_water
    coord           => domain%coord
    trans1          => domain%trans_vec
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    density         => enefunc%table%density
    table_ene_WW    => enefunc%table%table_ene_WW
    table_de_WW     => enefunc%table%table_de_WW
    cutoff          => enefunc%cutoffdist

    nb15_cellww     => pairlist%nb15_cellww
    nb15_listww     => pairlist%nb15_listww
    num_nb15_calcww => pairlist%num_nb15_calcww
    nb15_calc_water => pairlist%nb15_calc_water

    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff
    cutoff_small    =  (cutoff+2.0_wp)*(cutoff+2.0_wp)

    water_type(1,1)     = 1
    water_type(2:3,1)   = 2
    water_type(1,2:3)   = 2
    water_type(2:3,2:3) = 3

    ! calculate energy and gradient
    !

    !$omp parallel default(shared)                                            &
    !$omp private(id, i, j, iwater, i1, ilist, rtmp, force_local,             &
    !$omp         jwater, j1, jlist, rij2, L, R, dij, term_lj, term_elec,     &
    !$omp         grad_coef, work, ij, ini_nb15, fin_nb15, num_nb15,          &
    !$omp         trans_x, trans_y, trans_z, iix, iy, num_count, k, c,        &
    !$omp         small_count, j_small, force_localj, L1, check_virial,       &
    !$omp         elec_temp, evdw_temp, type)                      
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif

    ! energy within a cell
    !
    do i = id+1, ncell_local, nthread

      do iwater = 1, nwater(i) - 1

        do i1 = 1,3
          ilist(i1) = nsolute(i) + 3*(iwater-1) + i1
          rtmp(1:3,i1) = coord_pbc(1:3,ilist(i1),i)
        end do

        force_local(1:3,1:3) = 0.0_wp
        num_count(1:3,1:3) = 0
        small_count = 0
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp

        do jwater = iwater+1, nwater(i)

          iy   = nsolute(i) + 3*(jwater-1) + 1
          dij(1) = rtmp(1,1) - coord_pbc(1,iy,i)
          dij(2) = rtmp(2,1) - coord_pbc(2,iy,i)
          dij(3) = rtmp(3,1) - coord_pbc(3,iy,i)
          rij2  = dij(1)**2 + dij(2)**2 + dij(3)**2

          if (rij2 < cutoff_small) then
            small_count = small_count + 1
            j_small(small_count) = jwater
          end if

        end do

        do k = 1, small_count
          jwater = j_small(k)

          do j1 = 1, 3
            iy   = nsolute(i) + 3*(jwater-1) + j1
            do i1 = 1, 3

              type = water_type(i1,j1)
              dij(1) = rtmp(1,i1) - coord_pbc(1,iy,i)
              dij(2) = rtmp(2,i1) - coord_pbc(2,iy,i)
              dij(3) = rtmp(3,i1) - coord_pbc(3,iy,i)
              rij2  = dij(1)**2 + dij(2)**2 + dij(3)**2
              rij2  = cutoff2*density/rij2
              L    = int(rij2)
              R    = rij2 - L

              L1   = 6*L - 5
              term_lj   = table_ene_WW(L1,type)
              term_elec = table_ene_WW(L1+1,type)
              grad_coef = table_ene_WW(L1+2,type)
              term_lj   = term_lj + R*table_ene_WW(L1+3,type)
              term_elec = term_elec + R*table_ene_WW(L1+4,type)
              grad_coef = grad_coef + R*table_ene_WW(L1+5,type)

              evdw_temp = evdw_temp + term_lj
              elec_temp = elec_temp + term_elec
              work(1) = grad_coef*dij(1)
              work(2) = grad_coef*dij(2)
              work(3) = grad_coef*dij(3)
              force_local(1,i1) = force_local(1,i1) - work(1)
              force_local(2,i1) = force_local(2,i1) - work(2)
              force_local(3,i1) = force_local(3,i1) - work(3)
              force_localj(1,i1,j1,k) = work(1)
              force_localj(2,i1,j1,k) = work(2)
              force_localj(3,i1,j1,k) = work(3)
            end do
          end do

        end do

        do k = 1, small_count
          jwater = j_small(k)
          do j1 = 1, 3
            iy   = nsolute(i) + 3*(jwater-1) + j1
            do i1 = 1, 3
              force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1)  &
                                   + force_localj(1:3,i1,j1,k)
            end do
          end do
        end do
        do i1 = 1, 3
          iy = ilist(i1)
          force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1) + force_local(1:3,i1)
        end do
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1)  = evdw(id+1) + evdw_temp

      end do
    end do

    ! interaction between different cells
    !
    do ij = id+1, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      num_nb15 = 0

      trans_x = real(cell_move(1,j,i),wp)*system_size(1)
      trans_y = real(cell_move(2,j,i),wp)*system_size(2)
      trans_z = real(cell_move(3,j,i),wp)*system_size(3)
      check_virial = virial_check(j,i)

      do iix = 1, nb15_cellww(ij)

        iwater = nb15_listww(iix,ij)
        force_local(1:3,1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp 

        do i1 = 1, 3
          ilist(i1) = nsolute(i) + 3*(iwater-1) + i1
          rtmp(1,i1)  = coord_pbc(1,ilist(i1),i) + trans_x
          rtmp(2,i1)  = coord_pbc(2,ilist(i1),i) + trans_y
          rtmp(3,i1)  = coord_pbc(3,ilist(i1),i) + trans_z
        end do

        ini_nb15 = num_nb15 + 1
        fin_nb15 = num_nb15 + num_nb15_calcww(iwater,ij)
        num_nb15 = fin_nb15

        num_count(1:3,1:3) = 0
        small_count = 0

        do k = ini_nb15, fin_nb15

          jwater = nb15_calc_water(k,ij)
          iy   = nsolute(j) + 3*(jwater-1) + 1
          dij(1) = rtmp(1,1) - coord_pbc(1,iy,j)
          dij(2) = rtmp(2,1) - coord_pbc(2,iy,j)
          dij(3) = rtmp(3,1) - coord_pbc(3,iy,j)
          rij2  = dij(1)**2 + dij(2)**2 + dij(3)**2

          if (rij2 < cutoff_small) then
            small_count = small_count + 1
            j_small(small_count) = jwater
          end if

        end do

        do k = 1, small_count
          jwater = j_small(k)

          do j1 = 1, 3
            iy   = nsolute(j) + 3*(jwater-1) + j1
            do i1 = 1, 3
              type = water_type(i1,j1)
              dij(1) = rtmp(1,i1) - coord_pbc(1,iy,j)
              dij(2) = rtmp(2,i1) - coord_pbc(2,iy,j)
              dij(3) = rtmp(3,i1) - coord_pbc(3,iy,j)
              rij2  = dij(1)**2 + dij(2)**2 + dij(3)**2
              rij2  = cutoff2*density/rij2
              L    = int(rij2)
              R    = rij2 - L

              L1   = 6*L - 5
              term_lj   = table_ene_WW(L1,type)
              term_elec = table_ene_WW(L1+1,type)
              grad_coef = table_ene_WW(L1+2,type)
              term_lj   = term_lj + R*table_ene_WW(L1+3,type)
              term_elec = term_elec + R*table_ene_WW(L1+4,type)
              grad_coef = grad_coef + R*table_ene_WW(L1+5,type)

              evdw_temp = evdw_temp + term_lj
              elec_temp = elec_temp + term_elec
              work(1) = grad_coef*dij(1)
              work(2) = grad_coef*dij(2)
              work(3) = grad_coef*dij(3)
              force_local(1,i1) = force_local(1,i1) - work(1)
              force_local(2,i1) = force_local(2,i1) - work(2)
              force_local(3,i1) = force_local(3,i1) - work(3)
              force_localj(1,i1,j1,k) = work(1)
              force_localj(2,i1,j1,k) = work(2)
              force_localj(3,i1,j1,k) = work(3)
            end do
          end do

        end do

        do k = 1, small_count
          jwater = j_small(k)
          do j1 = 1, 3
            iy   = nsolute(j) + 3*(jwater-1) + j1
            do i1 = 1, 3
              force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1)  &
                                   + force_localj(1:3,i1,j1,k)
            end do
          end do
        end do
        do i1 = 1, 3
          iy = ilist(i1)
          force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1) + force_local(1:3,i1)
        end do
        if (check_virial == 1) then
          do i1 = 1, 3
            virial(1:3,ij) = virial(1:3,ij) - force_local(1:3,i1)
          end do
        end if
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond_table_water_linear

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_table_water_linear
  !> @brief        calculate nonbonded force among water molecules using lookup 
  !!               table (linear interpolation)
  !! @authors      JJ 
  !! @param[in]    domain    : domain information
  !! @param[in]    enefunc   : potential energy functions
  !! @param[in]    pairlist  : interaction list in each domain
  !! @param[inout] coord_pbc : pbc oriented coordinates
  !! @param[inout] force     : forces for each cell
  !! @param[inout] virial    : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_table_water_linear(domain, enefunc, &
                                               pairlist, coord_pbc,    &
                                               force, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(wip),                intent(inout) :: virial(:,:)

    ! local variables
    real(wp)                  :: dij(1:3), rij2, R
    real(wp)                  :: cutoff2, cutoff_small, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: trans_x, trans_y, trans_z
    real(wp)                  :: rtmp(3,3)
    real(wp)                  :: force_local(1:3,1:3)
    real(wp)                  :: force_localj(3,3,3,MaxWater)
    integer                   :: i, iy, j, k, ij, iix, iwater, jwater
    integer                   :: num_nb15, ini_nb15, fin_nb15
    integer                   :: id, omp_get_thread_num
    integer                   :: i1, j1, L
    integer                   :: num_count(3,3), c
    integer                   :: small_count, j_small(MaxWater)
    integer                   :: ilist(3), jlist(3)
    integer                   :: ncell_local
    integer                   :: water_type(3,3), type
    integer                   :: check_virial

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), system_size(:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_de_WW(:,:)
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: natom(:)
    integer(int2),    pointer :: cell_pairlist(:,:)
    integer,          pointer :: nwater(:)
    integer,          pointer :: nsolute(:)
    integer,          pointer :: nb15_cellww(:), nb15_listww(:,:)
    integer,          pointer :: num_nb15_calcww(:,:), nb15_calc_water(:,:)
    integer,          pointer :: virial_check(:,:)


    cell_pairlist   => domain%cell_pairlist1
    natom           => domain%num_atom
    nwater          => domain%num_water
    nsolute         => domain%num_solute
    coord           => domain%coord
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    trans1          => domain%trans_vec
    virial_check    => domain%virial_check

    density         => enefunc%table%density
    table_de_WW     => enefunc%table%table_de_WW
    cutoff          => enefunc%cutoffdist

    nb15_cellww     => pairlist%nb15_cellww
    nb15_listww     => pairlist%nb15_listww
    num_nb15_calcww => pairlist%num_nb15_calcww
    nb15_calc_water => pairlist%nb15_calc_water

    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff
    cutoff_small    =  (cutoff+2.0_wp)*(cutoff+2.0_wp)

    water_type(1,1)     = 1
    water_type(1,2:3)   = 2
    water_type(2:3,1)   = 2
    water_type(2:3,2:3) = 3

    ! calculate energy and gradient
    !
    !$omp parallel default(shared)                                             &
    !$omp private(id, i, j, iwater, i1, ilist, rtmp, force_local,              &
    !$omp         jwater, j1, jlist, rij2, L, R, grad_coef, work, dij, ij,     &
    !$omp         ini_nb15, fin_nb15, num_nb15, trans_x, trans_y, trans_z,     &
    !$omp         iix, iy, num_count, k, c, small_count, j_small, type,        &
    !$omp         check_virial, force_localj)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    ! energy within a cell
    !
    do i = id+1, ncell_local, nthread

      do iwater = 1, nwater(i) - 1

        do i1 = 1,3
          ilist(i1) = nsolute(i) + 3*(iwater-1)+ i1
          rtmp(1:3,i1) = coord_pbc(1:3,ilist(i1),i)
        end do

        force_local(1:3,1:3) = 0.0_wp

        num_count(1:3,1:3) = 0
        small_count = 0
        do jwater = iwater+1, nwater(i)

          iy   = nsolute(i) + 3*(jwater-1) + 1
          dij(1) = rtmp(1,1) - coord_pbc(1,iy,i) 
          dij(2) = rtmp(2,1) - coord_pbc(2,iy,i)
          dij(3) = rtmp(3,1) - coord_pbc(3,iy,i)
          rij2  = dij(1)**2 + dij(2)**2 + dij(3)**2

          if (rij2 < cutoff_small) then
            small_count = small_count + 1
            j_small(small_count) = jwater
          end if

        end do

        do k = 1, small_count
          jwater = j_small(k)

          do j1 = 1, 3
            iy = nsolute(i) + 3*(jwater-1) + j1
            do i1 = 1, 3
              type = water_type(i1,j1)
              dij(1) = rtmp(1,i1) - coord_pbc(1,iy,i)
              dij(2) = rtmp(2,i1) - coord_pbc(2,iy,i)
              dij(3) = rtmp(3,i1) - coord_pbc(3,iy,i)
              rij2  = dij(1)**2 + dij(2)**2 + dij(3)**2
              rij2  = cutoff2 * density / rij2
              L    = int(rij2)
              R    = rij2 - L
              grad_coef = table_de_WW(L,type)   &
                        + R*(table_de_WW(L+1,type)-table_de_WW(L,type))
              work(1) = grad_coef*dij(1)
              work(2) = grad_coef*dij(2)
              work(3) = grad_coef*dij(3)
              force_local(1,i1) = force_local(1,i1) - work(1)
              force_local(2,i1) = force_local(2,i1) - work(2)
              force_local(3,i1) = force_local(3,i1) - work(3)
              force_localj(1,i1,j1,k) = work(1)
              force_localj(2,i1,j1,k) = work(2)
              force_localj(3,i1,j1,k) = work(3)
            end do
          end do
        end do
        
        do k = 1, small_count
          jwater = j_small(k)
          do j1 = 1, 3
            iy = nsolute(i) + 3*(jwater-1) + j1
            do i1 = 1, 3
              force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1)   &
                                   + force_localj(1:3,i1,j1,k)
            end do
          end do
        end do
        do i1 = 1, 3
          iy = ilist(i1) 
          force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1) + force_local(1:3,i1)
        end do

      end do
    end do

    ! interaction between different cells
    !
    do ij = id+1, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      num_nb15 = 0

      trans_x = real(cell_move(1,j,i),wp)*system_size(1)
      trans_y = real(cell_move(2,j,i),wp)*system_size(2)
      trans_z = real(cell_move(3,j,i),wp)*system_size(3)
      check_virial = virial_check(j,i)

      do iix = 1, nb15_cellww(ij)

        iwater = nb15_listww(iix,ij)
        force_local(1:3,1:3) = 0.0_wp

        do i1 = 1, 3
          ilist(i1) = nsolute(i) + 3*(iwater-1) + i1
          rtmp(1,i1)  = coord_pbc(1,ilist(i1),i) + trans_x
          rtmp(2,i1)  = coord_pbc(2,ilist(i1),i) + trans_y
          rtmp(3,i1)  = coord_pbc(3,ilist(i1),i) + trans_z
        end do

        ini_nb15 = num_nb15 + 1
        fin_nb15 = num_nb15 + num_nb15_calcww(iwater,ij)
        num_nb15 = fin_nb15

        num_count(1:3,1:3) = 0
        small_count = 0

        do k = ini_nb15, fin_nb15
          jwater = nb15_calc_water(k,ij)

          ! first check O-O distance
          !
          iy   = nsolute(j) + 3*(jwater-1) + 1
          dij(1) = rtmp(1,1) - coord_pbc(1,iy,j) 
          dij(2) = rtmp(2,1) - coord_pbc(2,iy,j)
          dij(3) = rtmp(3,1) - coord_pbc(3,iy,j) 
          rij2  = dij(1)**2 + dij(2)**2 + dij(3)**2

          if (rij2 < cutoff_small) then
            small_count = small_count + 1
            j_small(small_count) = jwater
          end if
            
        end do

        do k = 1, small_count
          jwater = j_small(k)

          do j1 = 1, 3
            iy   = nsolute(j) + 3*(jwater-1) + j1
            do i1 = 1, 3
              type = water_type(i1,j1)
              dij(1) = rtmp(1,i1) - coord_pbc(1,iy,j)
              dij(2) = rtmp(2,i1) - coord_pbc(2,iy,j)
              dij(3) = rtmp(3,i1) - coord_pbc(3,iy,j)
              rij2  = dij(1)**2 + dij(2)**2 + dij(3)**2
              rij2  = cutoff2 * density / rij2
              L    = int(rij2)
              R    = rij2 - L
              grad_coef = table_de_WW(L,type)   &
                        + R*(table_de_WW(L+1,type)-table_de_WW(L,type))
              work(1) = grad_coef*dij(1)
              work(2) = grad_coef*dij(2)
              work(3) = grad_coef*dij(3)
              force_local(1,i1) = force_local(1,i1) - work(1)
              force_local(2,i1) = force_local(2,i1) - work(2)
              force_local(3,i1) = force_local(3,i1) - work(3)
              force_localj(1,i1,j1,k) = work(1)
              force_localj(2,i1,j1,k) = work(2)
              force_localj(3,i1,j1,k) = work(3)
            end do
          end do
        end do

        do k = 1, small_count
          jwater = j_small(k)
          do j1 = 1, 3
            iy   = nsolute(j) + 3*(jwater-1) + j1
            do i1 = 1, 3
              force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1)   &
                                   + force_localj(1:3,i1,j1,k)
            end do
          end do
        end do
        do i1 = 1, 3
          iy = ilist(i1) 
          force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1) + force_local(1:3,i1)
        end do
        if (check_virial == 1) then
          do i1 = 1, 3
            virial(1:3,ij) = virial(1:3,ij) - force_local(1:3,i1)
          end do
        end if

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_force_nonbond_table_water_linear

end module sp_energy_table_linear_water_mod
