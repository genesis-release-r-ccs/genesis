!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_table_linear_bondcorr_mod
!> @brief   calculate bond correction with linear interpolation table
!! @authors Jaewoon Jung(JJ)
!  
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_table_linear_bondcorr_mod

  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_energy_str_mod
  use sp_domain_str_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  !
  public   :: pme_bond_corr_linear
  private  :: pme_bond_corr_general_intel
  private  :: pme_bond_corr_general_fugaku
  private  :: pme_bond_corr_general_generic
  private  :: pme_bond_corr_gro_amber_intel
  private  :: pme_bond_corr_gro_amber_fugaku
  private  :: pme_bond_corr_gro_amber_generic
  private  :: pme_bond_corr_general_check
  private  :: pme_bond_corr_gro_amber_check
  private  :: pme_bond_corr_ljpme_intel
  private  :: pme_bond_corr_ljpme_fugaku
  private  :: pme_bond_corr_ljpme_generic
  private  :: pme_bond_corr_ljpme_check
  ! FEP
  private  :: pme_bond_corr_gro_amber_generic_fep
  private  :: pme_bond_corr_gro_amber_check_fep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_linear
  !> @brief        calculate pme bond correction with linear lookup table
  !  @authors      JJ
  !! @param[in]    domain    : domain information
  !! @param[in]    enefunc   : potential energy functions
  !! @param[in]    atmcls    : atom class number
  !! @param[in]    coord     : coordinates for each cell
  !! @param[inout] force     : forces for each cell
  !! @param[inout] force_pbc : force for each cell
  !! @param[inout] virial    : virial term of target systems
  !! @param[inout] eelec     : electrostatic energy of target systems
  !! @param[inout] evdw      : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_linear(domain, enefunc, atmcls, coord, force, &
                                  force_pbc, virial, evdw, eelec)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    integer,                 intent(in)    :: atmcls(:)
    real(wp),                intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial(:,:,:)
    real(dp),                intent(inout) :: evdw (nthread)
    real(dp),                intent(inout) :: eelec(nthread)


    if (enefunc%vdw == VDWPME) then 

      if (enefunc%nonb_limiter) then

        call pme_bond_corr_ljpme_check( &
                              domain, enefunc, coord, force, virial, &
                              evdw, eelec)

      else

        if (domain%nonbond_kernel == NBK_Intel) then
          call pme_bond_corr_ljpme_intel( &
                              domain, enefunc, atmcls, coord,        &
                              force_pbc, virial, evdw, eelec)
        else if (domain%nonbond_kernel == NBK_Fugaku) then
          call pme_bond_corr_ljpme_fugaku( &
                              domain, enefunc, atmcls, coord,        &
                              force_pbc, virial, evdw, eelec)
        else if (domain%nonbond_kernel == NBK_GPU) then
          call pme_bond_corr_ljpme_gpu( &
                              domain, enefunc, atmcls, coord,        &
                              force, virial, evdw, eelec)
        else
          call pme_bond_corr_ljpme_generic( &
                              domain, enefunc, coord, force, virial, &
                              evdw, eelec)
        end if

      end if

    else

      if (enefunc%nonb_limiter) then

        call pme_bond_corr_general_check( &
                              domain, enefunc, coord, force, virial, &
                              eelec)

      else

        if (domain%nonbond_kernel == NBK_Intel) then
           call pme_bond_corr_general_intel( &
                              domain, enefunc, coord, force_pbc,     &
                              virial, eelec)
        else if (domain%nonbond_kernel == NBK_Fugaku) then
           call pme_bond_corr_general_fugaku( &
                              domain, enefunc, coord, force_pbc,     &
                              virial, eelec)
        else if (domain%nonbond_kernel == NBK_GPU) then
           call pme_bond_corr_general_gpu( &
                              domain, enefunc, coord, force,         &
                              virial, eelec)
        else
          call pme_bond_corr_general_generic( &
                              domain, enefunc, coord, force, virial, &
                              eelec)
        end if

      end if

    end if

    if (enefunc%forcefield /= ForcefieldCHARMM) then

      if (enefunc%nonb_limiter) then

        call pme_bond_corr_gro_amber_check( &
                              domain, enefunc, coord, force, virial, &
                              eelec)

        if (domain%fep_use) then
          call pme_bond_corr_gro_amber_check_fep( &
                                domain, enefunc, coord, force, virial, &
                                eelec)
        end if

      else

        if (domain%nonbond_kernel == NBK_Intel) then
          call pme_bond_corr_gro_amber_intel( &
                              domain, enefunc, coord, force_pbc,     &
                              virial, eelec)
        else if (domain%nonbond_kernel == NBK_Fugaku) then
          call pme_bond_corr_gro_amber_fugaku( &
                              domain, enefunc, coord, force_pbc,     &
                              virial, eelec)
        else if (domain%nonbond_kernel == NBK_GPU) then
          call pme_bond_corr_gro_amber_gpu( &
                              domain, enefunc, coord, force,         &
                              virial, eelec)
        else
          call pme_bond_corr_gro_amber_generic( &
                              domain, enefunc, coord, force, virial, &
                              eelec)
        end if

        if (domain%fep_use) then
          call pme_bond_corr_gro_amber_generic_fep( &
                              domain, enefunc, coord, force, virial, &
                              eelec)
        end if

      end if

    end if

    return

   end subroutine pme_bond_corr_linear

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_general_intel
  !> @brief        calculate bond correction term in PME (general, intel kernel)
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates for each cell
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_general_intel(domain, enefunc, coord, force, &
                                         virial, eelec)

    ! formal arguments 
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), viri(3), cutoff2
    real(wp)                 :: R, term_elec, ccr, coef, elec_temp
    real(wp)                 :: iqtmp, jqtmp
    integer                  :: i, ix, iy, j, k, ij, L
    integer                  :: id, omp_get_thread_num
    integer                  :: ncell_local

    real(wp),        pointer :: table_ecor(:), table_decor(:)
    real(wp),        pointer :: density, cutoff
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: natom(:), nwater(:), water_list(:,:,:)
    integer,         pointer :: num_nonb_excl(:), nonb_excl_list(:,:,:)


    cell_pair       => domain%cell_pairlist1
    natom           => domain%num_atom
    nwater          => domain%num_water
    water_list      => domain%water_list

    table_ecor      => enefunc%table%table_ecor
    table_decor     => enefunc%table%table_decor
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    num_nonb_excl   => enefunc%num_nonb_excl
    nonb_excl_list  => enefunc%nonb_excl_list

    ncell_local = domain%num_cell_local
    cutoff2     = cutoff*cutoff

    !$omp parallel default(shared)                                      &
    !$omp private(id, i, j, k, ij, ix, iy, dij, rij2, L, R, coef, work, &
    !$omp         viri, term_elec,  ccr, elec_temp, iqtmp, jqtmp)
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
      viri(1:3) = 0.0_wp

      do k = 1, num_nonb_excl(ij)

        i  = nonb_excl_list(1,k,ij)
        j  = nonb_excl_list(2,k,ij)

        ! compute distance
        !
        dij(1) = coord(i,1,1) - coord(j,1,1)
        dij(2) = coord(i,2,1) - coord(j,2,1)
        dij(3) = coord(i,3,1) - coord(j,3,1)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        iqtmp = coord(i,4,1)
        jqtmp = coord(j,4,1)

        ! energy and gradient
        !
        rij2 = cutoff2*density/rij2
        L    = int(rij2)
        R    = rij2 - L

        term_elec = table_ecor(L) + R*(table_ecor(L+1)-table_ecor(L))
        ccr       = iqtmp * jqtmp * term_elec
        elec_temp = elec_temp + ccr

        term_elec = table_decor(L) + R*(table_decor(L+1)-table_decor(L))
        coef      = iqtmp * jqtmp * term_elec

        work(1) = coef*dij(1)
        work(2) = coef*dij(2)
        work(3) = coef*dij(3)

        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force(i,1,1,id+1) = force(i,1,1,id+1) - work(1)
        force(i,2,1,id+1) = force(i,2,1,id+1) - work(2)
        force(i,3,1,id+1) = force(i,3,1,id+1) - work(3)
        force(j,1,1,id+1) = force(j,1,1,id+1) + work(1)
        force(j,2,1,id+1) = force(j,2,1,id+1) + work(2)
        force(j,3,1,id+1) = force(j,3,1,id+1) + work(3)

      end do

      virial(1,1,id+1) = virial(1,1,id+1) - viri(1) 
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2) 
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3) 
      eelec(id+1) = eelec(id+1) + elec_temp  

    end do

    !$omp end parallel

    return

  end subroutine pme_bond_corr_general_intel

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_general_fugaku
  !> @brief        calculate bond correction term in PME (general, intel kernel)
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates for each cell
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_general_fugaku(domain, enefunc, coord, force, &
                                          virial, eelec)

    ! formal arguments 
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), viri(3), cutoff2
    real(wp)                 :: R, term_elec, ccr, coef, elec_temp
    real(wp)                 :: iqtmp, jqtmp
    integer                  :: i, ix, iy, j, k, ij, L
    integer                  :: id, omp_get_thread_num
    integer                  :: ncell_local

    real(wp),        pointer :: table_ecor(:), table_decor(:)
    real(wp),        pointer :: density, cutoff
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: natom(:), nwater(:), water_list(:,:,:)
    integer,         pointer :: num_nonb_excl(:), nonb_excl_list(:,:,:)


    cell_pair       => domain%cell_pairlist1
    natom           => domain%num_atom
    nwater          => domain%num_water
    water_list      => domain%water_list

    table_ecor      => enefunc%table%table_ecor
    table_decor     => enefunc%table%table_decor
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    num_nonb_excl   => enefunc%num_nonb_excl
    nonb_excl_list  => enefunc%nonb_excl_list

    ncell_local = domain%num_cell_local
    cutoff2     = cutoff*cutoff

    !$omp parallel default(shared)                                      &
    !$omp private(id, i, j, k, ij, ix, iy, dij, rij2, L, R, coef, work, &
    !$omp         viri, term_elec,  ccr, elec_temp, iqtmp, jqtmp)
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
      viri(1:3) = 0.0_wp

      do k = 1, num_nonb_excl(ij)

        i  = nonb_excl_list(1,k,ij)
        j  = nonb_excl_list(2,k,ij)

        ! compute distance
        !
        dij(1) = coord(1,i,1) - coord(1,j,1)
        dij(2) = coord(2,i,1) - coord(2,j,1)
        dij(3) = coord(3,i,1) - coord(3,j,1)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        iqtmp = coord(4,i,1)
        jqtmp = coord(4,j,1)

        ! energy and gradient
        !
        rij2 = cutoff2*density/rij2
        L    = int(rij2)
        R    = rij2 - L

        term_elec = table_ecor(L) + R*(table_ecor(L+1)-table_ecor(L))
        ccr       = iqtmp * jqtmp * term_elec
        elec_temp = elec_temp + ccr

        term_elec = table_decor(L) + R*(table_decor(L+1)-table_decor(L))
        coef      = iqtmp * jqtmp * term_elec

        work(1) = coef*dij(1)
        work(2) = coef*dij(2)
        work(3) = coef*dij(3)

        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force(1,i,1,id+1) = force(1,i,1,id+1) - work(1)
        force(2,i,1,id+1) = force(2,i,1,id+1) - work(2)
        force(3,i,1,id+1) = force(3,i,1,id+1) - work(3)
        force(1,j,1,id+1) = force(1,j,1,id+1) + work(1)
        force(2,j,1,id+1) = force(2,j,1,id+1) + work(2)
        force(3,j,1,id+1) = force(3,j,1,id+1) + work(3)

      end do

      virial(1,1,id+1) = virial(1,1,id+1) - viri(1) 
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2) 
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3) 
      eelec(id+1) = eelec(id+1) + elec_temp  

    end do

    !$omp end parallel

    return

  end subroutine pme_bond_corr_general_fugaku

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_general_gpu
  !> @brief        calculate bond correction term in PME (general, intel kernel)
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates for each cell
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_general_gpu(domain, enefunc, coord, force, &
                                       virial, eelec)

    ! formal arguments 
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), viri(3), cutoff2
    real(wp)                 :: R, term_elec, ccr, coef, elec_temp
    real(wp)                 :: iqtmp, jqtmp
    integer                  :: i, ix, iy, j, k, ij, L, icel, jcel
    integer                  :: id, omp_get_thread_num
    integer                  :: ncell_local, num_atom

    real(wp),        pointer :: system_size(:)
    real(wp),        pointer :: table_ecor(:), table_decor(:)
    real(wp),        pointer :: density, cutoff
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: natom(:), nwater(:), water_list(:,:,:)
    integer,         pointer :: start_atom(:)
    integer,         pointer :: num_nonb_excl(:), nonb_excl_list(:,:,:)
    integer(1),      pointer :: cell_move(:,:,:)


    cell_pair       => domain%cell_pairlist1
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    nwater          => domain%num_water
    water_list      => domain%water_list
    system_size     => domain%system_size
    cell_move       => domain%cell_move

    table_ecor      => enefunc%table%table_ecor
    table_decor     => enefunc%table%table_decor
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    num_nonb_excl   => enefunc%num_nonb_excl
    nonb_excl_list  => enefunc%nonb_excl_list

    ncell_local     = domain%num_cell_local
    cutoff2         = cutoff*cutoff
    num_atom        = domain%num_atom_domain

    !$omp parallel default(shared)                                      &
    !$omp private(id, i, j, k, ij, ix, iy, dij, rij2, L, R, coef, work, &
    !$omp         viri, icel, jcel, term_elec, ccr, elec_temp, iqtmp, jqtmp)
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
      viri(1:3) = 0.0_wp

      do k = 1, num_nonb_excl(ij)

        icel  = nonb_excl_list(1,k,ij)
        jcel  = nonb_excl_list(2,k,ij)
        ix    = nonb_excl_list(3,k,ij)
        iy    = nonb_excl_list(4,k,ij)
        i     = start_atom(icel) + ix
        j     = start_atom(jcel) + iy

        ! compute distance
        !
        dij(1) = coord(           i,1,1) - coord(           j,1,1)
        dij(2) = coord(  num_atom+i,1,1) - coord(  num_atom+j,1,1)
        dij(3) = coord(2*num_atom+i,1,1) - coord(2*num_atom+j,1,1)
        dij(1) = dij(1) + cell_move(1,jcel,icel)*system_size(1)
        dij(2) = dij(2) + cell_move(2,jcel,icel)*system_size(2)
        dij(3) = dij(3) + cell_move(3,jcel,icel)*system_size(3)

        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        iqtmp = coord(3*num_atom+i,1,1)
        jqtmp = coord(3*num_atom+j,1,1)

        ! energy and gradient
        !
        rij2 = cutoff2*density/rij2
        L    = int(rij2)
        R    = rij2 - L

        term_elec = table_ecor(L) + R*(table_ecor(L+1)-table_ecor(L))
        ccr       = iqtmp * jqtmp * term_elec
        elec_temp = elec_temp + ccr

        term_elec = table_decor(L) + R*(table_decor(L+1)-table_decor(L))
        coef      = iqtmp * jqtmp * term_elec

        work(1) = coef*dij(1)
        work(2) = coef*dij(2)
        work(3) = coef*dij(3)

        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force(1,ix,icel,id+1) = force(1,ix,icel,id+1) - work(1)
        force(2,ix,icel,id+1) = force(2,ix,icel,id+1) - work(2)
        force(3,ix,icel,id+1) = force(3,ix,icel,id+1) - work(3)
        force(1,iy,jcel,id+1) = force(1,iy,jcel,id+1) + work(1)
        force(2,iy,jcel,id+1) = force(2,iy,jcel,id+1) + work(2)
        force(3,iy,jcel,id+1) = force(3,iy,jcel,id+1) + work(3)

      end do

      virial(1,1,id+1) = virial(1,1,id+1) - viri(1) 
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2) 
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3) 
      eelec(id+1) = eelec(id+1) + elec_temp  

    end do

    !$omp end parallel

    return

  end subroutine pme_bond_corr_general_gpu

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_general_generic
  !> @brief        calculate bond correction term in PME (general)
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates for each cell
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_general_generic(domain, enefunc, coord,  &
                                           force, virial, eelec)

    ! formal arguments 
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), viri(3), cutoff2
    real(wp)                 :: R, term_elec, ccr, coef, elec_temp
    integer                  :: i, ix, iy, j, k, ij, L
    integer                  :: id, omp_get_thread_num
    integer                  :: ncell_local

    real(wp),        pointer :: charge(:,:), system_size(:)
    real(wp),        pointer :: table_ecor(:), table_decor(:)
    real(wp),        pointer :: density, cutoff
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: natom(:), nwater(:), water_list(:,:,:)
    integer,         pointer :: num_nonb_excl(:), nonb_excl_list(:,:,:)
    integer(1),      pointer :: cell_move(:,:,:)


    cell_pair       => domain%cell_pairlist1
    natom           => domain%num_atom
    nwater          => domain%num_water
    water_list      => domain%water_list
    charge          => domain%charge
    system_size     => domain%system_size
    cell_move       => domain%cell_move

    table_ecor      => enefunc%table%table_ecor
    table_decor     => enefunc%table%table_decor
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    num_nonb_excl   => enefunc%num_nonb_excl
    nonb_excl_list  => enefunc%nonb_excl_list

    ncell_local = domain%num_cell_local
    cutoff2     = cutoff*cutoff

    !$omp parallel default(shared)                                      &
    !$omp private(id, i, j, k, ij, ix, iy, dij, rij2, L, R, coef, work, &
    !$omp         viri, term_elec,  ccr, elec_temp)
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
      viri(1:3) = 0.0_wp

      do k = 1, num_nonb_excl(ij)

        i  = nonb_excl_list(1,k,ij)
        j  = nonb_excl_list(2,k,ij)
        ix = nonb_excl_list(3,k,ij)
        iy = nonb_excl_list(4,k,ij)

        ! compute distance
        !
        dij(1) = coord(ix,1,i) - coord(iy,1,j)
        dij(2) = coord(ix,2,i) - coord(iy,2,j)
        dij(3) = coord(ix,3,i) - coord(iy,3,j)
        dij(1) = dij(1) + real(cell_move(1,j,i),wp)*system_size(1)
        dij(2) = dij(2) + real(cell_move(2,j,i),wp)*system_size(2)
        dij(3) = dij(3) + real(cell_move(3,j,i),wp)*system_size(3)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        ! energy and gradient
        !
        rij2 = cutoff2*density/rij2
        L    = int(rij2)
        R    = rij2 - L

        term_elec = table_ecor(L) + R*(table_ecor(L+1)-table_ecor(L))
        ccr       = charge(ix,i) * charge(iy,j) * term_elec
        elec_temp = elec_temp + ccr

        term_elec = table_decor(L) + R*(table_decor(L+1)-table_decor(L))
        coef      = charge(ix,i) * charge(iy,j) * term_elec

        work(1) = coef*dij(1)
        work(2) = coef*dij(2)
        work(3) = coef*dij(3)

        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force(1,ix,i,id+1) = force(1,ix,i,id+1) - work(1)
        force(2,ix,i,id+1) = force(2,ix,i,id+1) - work(2)
        force(3,ix,i,id+1) = force(3,ix,i,id+1) - work(3)
        force(1,iy,j,id+1) = force(1,iy,j,id+1) + work(1)
        force(2,iy,j,id+1) = force(2,iy,j,id+1) + work(2)
        force(3,iy,j,id+1) = force(3,iy,j,id+1) + work(3)

      end do

      virial(1,1,id+1) = virial(1,1,id+1) - viri(1) 
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2) 
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3) 
      eelec(id+1) = eelec(id+1) + elec_temp  

    end do

    !$omp end parallel

    return

  end subroutine pme_bond_corr_general_generic

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_gro_amber_intel
  !> @brief        calculate bodn correction relating with 14 scaling
  !  @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[in]    coord   : coordinates for each cell
  !! @param[inout] force   : forces for each cell
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_gro_amber_intel(domain, enefunc, coord, &
                                           force, virial, eelec)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), viri(3), cutoff2
    real(wp)                 :: R, term_elec, coef, cc, qq_scale
    real(wp)                 :: iqtmp, jqtmp
    real(wp)                 :: elec_temp
    integer                  :: i, ix, iy, j, k, ij, L
    integer                  :: ncell_local, id, omp_get_thread_num

    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: natom(:)
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: num_nb14_calc(:), nb14_calc_list(:,:,:)


    cell_pair       => domain%cell_pairlist1
    natom           => domain%num_atom

    table_ene       => enefunc%table%table_ecor
    table_grad      => enefunc%table%table_decor
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    num_nb14_calc   => enefunc%num_nb14_calc
    nb14_calc_list  => enefunc%nb14_calc_list

    cutoff2         =  cutoff * cutoff
    ncell_local     =  domain%num_cell_local

    !$omp parallel default(shared)                                      &
    !$omp private(id, i, j, k, ij, ix, iy, dij, rij2, L, R, coef, work, &
    !$omp         term_elec, cc, qq_scale, iqtmp, jqtmp, viri, elec_temp)
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
      viri(1:3) = 0.0_wp

      do k = 1, num_nb14_calc(ij)

        i  = nb14_calc_list(1,k,ij)
        j  = nb14_calc_list(2,k,ij)

        ! compute distance
        !
        dij(1) = coord(i,1,1) - coord(j,1,1)
        dij(2) = coord(i,2,1) - coord(j,2,1)
        dij(3) = coord(i,3,1) - coord(j,3,1)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        iqtmp  = coord(i,4,1)
        jqtmp  = coord(j,4,1)

        ! energy and gradient
        !
        rij2 = cutoff2*density/rij2
        L    = int(rij2)
        R    = rij2 - L

        qq_scale  = -enefunc%nb14_qq_scale(k,ij)+1.0_wp
        term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
        cc        = iqtmp*jqtmp*qq_scale
        elec_temp = elec_temp + term_elec*cc

        term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
        coef      = cc*term_elec

        work(1) = coef*dij(1)
        work(2) = coef*dij(2)
        work(3) = coef*dij(3)

        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force(i,1,1,id+1) = force(i,1,1,id+1) - work(1)
        force(i,2,1,id+1) = force(i,2,1,id+1) - work(2)
        force(i,3,1,id+1) = force(i,3,1,id+1) - work(3)
        force(j,1,1,id+1) = force(j,1,1,id+1) + work(1)
        force(j,2,1,id+1) = force(j,2,1,id+1) + work(2)
        force(j,3,1,id+1) = force(j,3,1,id+1) + work(3)

      end do

      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,1,id+1) = virial(3,3,id+1) - viri(3)
      eelec(id+1) = eelec(id+1) + elec_temp

    end do

    !$omp end parallel

    return

  end subroutine pme_bond_corr_gro_amber_intel

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_gro_amber_fugaku
  !> @brief        calculate bodn correction relating with 14 scaling
  !  @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[in]    coord   : coordinates for each cell
  !! @param[inout] force   : forces for each cell
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_gro_amber_fugaku(domain, enefunc, coord, &
                                            force, virial, eelec)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), viri(3), cutoff2
    real(wp)                 :: R, term_elec, coef, cc, qq_scale
    real(wp)                 :: iqtmp, jqtmp
    real(wp)                 :: elec_temp
    integer                  :: i, ix, iy, j, k, ij, L
    integer                  :: ncell_local, id, omp_get_thread_num

    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: natom(:)
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: num_nb14_calc(:), nb14_calc_list(:,:,:)


    cell_pair       => domain%cell_pairlist1
    natom           => domain%num_atom

    table_ene       => enefunc%table%table_ecor
    table_grad      => enefunc%table%table_decor
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    num_nb14_calc   => enefunc%num_nb14_calc
    nb14_calc_list  => enefunc%nb14_calc_list

    cutoff2         =  cutoff * cutoff
    ncell_local     =  domain%num_cell_local

    !$omp parallel default(shared)                                      &
    !$omp private(id, i, j, k, ij, ix, iy, dij, rij2, L, R, coef, work, &
    !$omp         term_elec, cc, qq_scale, iqtmp, jqtmp, viri, elec_temp)
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
      viri(1:3) = 0.0_wp

      do k = 1, num_nb14_calc(ij)

        i  = nb14_calc_list(1,k,ij)
        j  = nb14_calc_list(2,k,ij)

        ! compute distance
        !
        dij(1) = coord(1,i,1) - coord(1,j,1)
        dij(2) = coord(2,i,1) - coord(2,j,1)
        dij(3) = coord(3,i,1) - coord(3,j,1)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        iqtmp  = coord(4,i,1)
        jqtmp  = coord(4,j,1)

        ! energy and gradient
        !
        rij2 = cutoff2*density/rij2
        L    = int(rij2)
        R    = rij2 - L

        qq_scale  = -enefunc%nb14_qq_scale(k,ij)+1.0_wp
        term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
        cc        = iqtmp*jqtmp*qq_scale
        elec_temp = elec_temp + term_elec*cc

        term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
        coef      = cc*term_elec

        work(1) = coef*dij(1)
        work(2) = coef*dij(2)
        work(3) = coef*dij(3)

        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force(1,i,1,id+1) = force(1,i,1,id+1) - work(1)
        force(2,i,1,id+1) = force(2,i,1,id+1) - work(2)
        force(3,i,1,id+1) = force(3,i,1,id+1) - work(3)
        force(1,j,1,id+1) = force(1,j,1,id+1) + work(1)
        force(2,j,1,id+1) = force(2,j,1,id+1) + work(2)
        force(3,j,1,id+1) = force(3,j,1,id+1) + work(3)

      end do

      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,1,id+1) = virial(3,3,id+1) - viri(3)
      eelec(id+1) = eelec(id+1) + elec_temp

    end do

    !$omp end parallel

    return

  end subroutine pme_bond_corr_gro_amber_fugaku

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_gro_amber_gpu
  !> @brief        calculate bodn correction relating with 14 scaling
  !  @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[in]    coord   : coordinates for each cell
  !! @param[inout] force   : forces for each cell
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_gro_amber_gpu(domain, enefunc, coord, &
                                         force, virial, eelec)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), viri(3), cutoff2
    real(wp)                 :: R, term_elec, coef, cc, qq_scale
    real(wp)                 :: iqtmp, jqtmp
    real(wp)                 :: elec_temp
    integer                  :: i, ix, iy, j, k, ij, L, num_atom, icel, jcel
    integer                  :: ncell_local, id, omp_get_thread_num

    real(wp),        pointer :: system_size(:)
    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: natom(:), start_atom(:)
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: num_nb14_calc(:), nb14_calc_list(:,:,:)
    integer(1),      pointer :: cell_move(:,:,:)


    cell_pair       => domain%cell_pairlist1
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    system_size     => domain%system_size
    cell_move       => domain%cell_move

    table_ene       => enefunc%table%table_ecor
    table_grad      => enefunc%table%table_decor
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    num_nb14_calc   => enefunc%num_nb14_calc
    nb14_calc_list  => enefunc%nb14_calc_list

    cutoff2         =  cutoff * cutoff
    ncell_local     =  domain%num_cell_local
    num_atom        =  domain%num_atom_domain

    !$omp parallel default(shared)                                      &
    !$omp private(id, i, j, k, ij, ix, iy, dij, rij2, L, R, coef, work, &
    !$omp         term_elec, cc, qq_scale, icel, jcel, iqtmp, jqtmp,    &
    !$omp         viri, elec_temp)
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
      viri(1:3) = 0.0_wp

      do k = 1, num_nb14_calc(ij)

        icel  = nb14_calc_list(1,k,ij)
        jcel  = nb14_calc_list(2,k,ij)
        ix    = nb14_calc_list(3,k,ij)
        iy    = nb14_calc_list(4,k,ij)
        i     = start_atom(icel) + ix
        j     = start_atom(jcel) + iy

        ! compute distance
        !
        dij(1) = coord(           i,1,1) - coord(           j,1,1)
        dij(2) = coord(  num_atom+i,1,1) - coord(  num_atom+j,1,1)
        dij(3) = coord(2*num_atom+i,1,1) - coord(2*num_atom+j,1,1)
        dij(1) = dij(1) + cell_move(1,jcel,icel)*system_size(1)
        dij(2) = dij(2) + cell_move(2,jcel,icel)*system_size(2)
        dij(3) = dij(3) + cell_move(3,jcel,icel)*system_size(3)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        iqtmp  = coord(3*num_atom+i,1,1)
        jqtmp  = coord(3*num_atom+j,1,1)

        ! energy and gradient
        !
        rij2 = cutoff2*density/rij2
        L    = int(rij2)
        R    = rij2 - L

        qq_scale  = -enefunc%nb14_qq_scale(k,ij)+1.0_wp
        term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
        cc        = iqtmp*jqtmp*qq_scale
        elec_temp = elec_temp + term_elec*cc

        term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
        coef      = cc*term_elec

        work(1) = coef*dij(1)
        work(2) = coef*dij(2)
        work(3) = coef*dij(3)

        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force(1,ix,icel,id+1) = force(1,ix,icel,id+1) - work(1)
        force(2,ix,icel,id+1) = force(2,ix,icel,id+1) - work(2)
        force(3,ix,icel,id+1) = force(3,ix,icel,id+1) - work(3)
        force(1,iy,jcel,id+1) = force(1,iy,jcel,id+1) + work(1)
        force(2,iy,jcel,id+1) = force(2,iy,jcel,id+1) + work(2)
        force(3,iy,jcel,id+1) = force(3,iy,jcel,id+1) + work(3)

      end do

      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,1,id+1) = virial(3,3,id+1) - viri(3)
      eelec(id+1) = eelec(id+1) + elec_temp

    end do

    !$omp end parallel

    return

  end subroutine pme_bond_corr_gro_amber_gpu

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_gro_amber_generic
  !> @brief        calculate bodn correction relating with 14 scaling
  !  @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[in]    coord   : coordinates for each cell
  !! @param[inout] force   : forces for each cell
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_gro_amber_generic(domain, enefunc, coord, &
                                             force, virial, eelec)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), viri(3), cutoff2
    real(wp)                 :: R, term_elec, coef, cc, qq_scale
    real(wp)                 :: elec_temp
    integer                  :: i, ix, iy, j, k, ij, L
    integer                  :: ncell_local, id, omp_get_thread_num

    real(wp),        pointer :: charge(:,:), system_size(:)
    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: natom(:)
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: num_nb14_calc(:), nb14_calc_list(:,:,:)
    integer(1),      pointer :: cell_move(:,:,:)


    charge          => domain%charge
    cell_pair       => domain%cell_pairlist1
    natom           => domain%num_atom
    system_size     => domain%system_size
    cell_move       => domain%cell_move

    table_ene       => enefunc%table%table_ecor
    table_grad      => enefunc%table%table_decor
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    num_nb14_calc   => enefunc%num_nb14_calc
    nb14_calc_list  => enefunc%nb14_calc_list

    cutoff2         =  cutoff * cutoff
    ncell_local     =  domain%num_cell_local

    !$omp parallel default(shared)                                      &
    !$omp private(id, i, j, k, ij, ix, iy, dij, rij2, L, R, coef, work, &
    !$omp         term_elec, cc, qq_scale, viri, elec_temp)
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
      viri(1:3) = 0.0_wp

      do k = 1, num_nb14_calc(ij)

        i  = nb14_calc_list(1,k,ij)
        j  = nb14_calc_list(2,k,ij)
        ix = nb14_calc_list(3,k,ij)
        iy = nb14_calc_list(4,k,ij)

        ! compute distance
        !
        dij(1) = coord(ix,1,i) - coord(iy,1,j)
        dij(2) = coord(ix,2,i) - coord(iy,2,j)
        dij(3) = coord(ix,3,i) - coord(iy,3,j)
        dij(1) = dij(1) + real(cell_move(1,j,i),wp)*system_size(1)
        dij(2) = dij(2) + real(cell_move(2,j,i),wp)*system_size(2)
        dij(3) = dij(3) + real(cell_move(3,j,i),wp)*system_size(3)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        ! energy and gradient
        !
        rij2 = cutoff2*density/rij2
        L    = int(rij2)
        R    = rij2 - L

        qq_scale  = -enefunc%nb14_qq_scale(k,ij)+1.0_wp
        term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
        cc        = charge(ix,i)*charge(iy,j)*qq_scale
        elec_temp = elec_temp + term_elec*cc

        term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
        coef      = cc*term_elec

        work(1) = coef*dij(1)
        work(2) = coef*dij(2)
        work(3) = coef*dij(3)

        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force(1,ix,i,id+1) = force(1,ix,i,id+1) - work(1)
        force(2,ix,i,id+1) = force(2,ix,i,id+1) - work(2)
        force(3,ix,i,id+1) = force(3,ix,i,id+1) - work(3)
        force(1,iy,j,id+1) = force(1,iy,j,id+1) + work(1)
        force(2,iy,j,id+1) = force(2,iy,j,id+1) + work(2)
        force(3,iy,j,id+1) = force(3,iy,j,id+1) + work(3)

      end do

      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,1,id+1) = virial(3,3,id+1) - viri(3)
      eelec(id+1) = eelec(id+1) + elec_temp

    end do

    !$omp end parallel

    return

  end subroutine pme_bond_corr_gro_amber_generic

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_general_check
  !> @brief        calculate bond correction term in PME (general)
  !! @authors      JJ, CK
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[in]    coord   : coordinates for each cell
  !! @param[inout] force   : forces for each cell
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_general_check(domain, enefunc, coord, force, &
                                         virial, eelec)

    ! formal arguments 
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), viri(3), cutoff2
    real(wp)                 :: R, term_elec, ccr, coef, elec_temp
    real(wp)                 :: minimum_contact
    integer                  :: i, ix, iy, j, k, ij, L
    integer                  :: id, omp_get_thread_num
    integer                  :: ncell_local

    real(wp),        pointer :: charge(:,:), system_size(:)
    real(wp),        pointer :: table_ecor(:), table_decor(:)
    real(wp),        pointer :: density, cutoff
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: natom(:), nwater(:), water_list(:,:,:)
    integer,         pointer :: num_nonb_excl(:), nonb_excl_list(:,:,:)
    integer(1),      pointer :: cell_move(:,:,:)


    cell_pair       => domain%cell_pairlist1
    natom           => domain%num_atom
    nwater          => domain%num_water
    water_list      => domain%water_list
    charge          => domain%charge
    system_size     => domain%system_size
    cell_move       => domain%cell_move

    table_ecor      => enefunc%table%table_ecor
    table_decor     => enefunc%table%table_decor
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    num_nonb_excl   => enefunc%num_nonb_excl
    nonb_excl_list  => enefunc%nonb_excl_list

    ncell_local = domain%num_cell_local
    cutoff2     = cutoff*cutoff
    minimum_contact =  enefunc%minimum_contact

    !$omp parallel default(shared)                                       &
    !$omp private(id, i, j, k, ij, ix, iy, dij, rij2, L, R, coef, work,  &
    !$omp         term_elec,  ccr, viri, elec_temp)
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
      viri(1:3) = 0.0_wp

      do k = 1, num_nonb_excl(ij)

        i  = nonb_excl_list(1,k,ij)
        j  = nonb_excl_list(2,k,ij)
        ix = nonb_excl_list(3,k,ij)
        iy = nonb_excl_list(4,k,ij)

        ! compute distance
        !
        dij(1) = coord(ix,1,i) - coord(iy,1,j)
        dij(2) = coord(ix,2,i) - coord(iy,2,j)
        dij(3) = coord(ix,3,i) - coord(iy,3,j)
        dij(1) = dij(1) + real(cell_move(1,j,i),wp)*system_size(1)
        dij(2) = dij(2) + real(cell_move(2,j,i),wp)*system_size(2)
        dij(3) = dij(3) + real(cell_move(3,j,i),wp)*system_size(3)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        rij2   = max(rij2, minimum_contact)

        ! energy and gradient
        !
        rij2 = cutoff2*density/rij2
        L    = int(rij2)
        R    = rij2 - L

        term_elec = table_ecor(L) + R*(table_ecor(L+1)-table_ecor(L))
        ccr       = charge(ix,i) * charge(iy,j) * term_elec
        elec_temp = elec_temp + ccr

        term_elec = table_decor(L) + R*(table_decor(L+1)-table_decor(L))
        coef      = charge(ix,i) * charge(iy,j) * term_elec

        work(1:3) = coef*dij(1:3)
        viri(1:3) = viri(1:3) + dij(1:3)*work(1:3)

        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
        force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + work(1:3)

      end do
      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
      eelec(id+1) = eelec(id+1) + elec_temp

    end do

    !$omp end parallel

    return

  end subroutine pme_bond_corr_general_check

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_gro_amber_check
  !> @brief        calculate bodn correction relating with 14 scaling
  !  @authors      JJ, CK
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[in]    coord   : coordinates for each cell
  !! @param[inout] force   : forces for each cell
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_gro_amber_check(domain, enefunc, coord, &
                                           force, virial, eelec)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), viri(3), cutoff2
    real(wp)                 :: R, term_elec, coef, cc, qq_scale
    real(wp)                 :: elec_temp, minimum_contact
    integer                  :: i, ix, iy, j, k, ij, L
    integer                  :: ncell_local, id, omp_get_thread_num

    real(wp),        pointer :: charge(:,:), system_size(:)
    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: natom(:)
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: num_nb14_calc(:), nb14_calc_list(:,:,:)
    integer(1),      pointer :: cell_move(:,:,:)


    charge          => domain%charge
    cell_pair       => domain%cell_pairlist1
    natom           => domain%num_atom
    system_size     => domain%system_size
    cell_move       => domain%cell_move

    table_ene       => enefunc%table%table_ecor
    table_grad      => enefunc%table%table_decor
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    num_nb14_calc   => enefunc%num_nb14_calc
    nb14_calc_list  => enefunc%nb14_calc_list

    cutoff2         =  cutoff * cutoff
    ncell_local     =  domain%num_cell_local
    minimum_contact =  enefunc%minimum_contact

    !$omp parallel default(shared)                                      &
    !$omp private(id, i, j, k, ij, ix, iy, dij, rij2, L, R, coef, work, &
    !$omp         term_elec, cc, qq_scale, viri, elec_temp)
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
      viri(1:3) = 0.0_wp

      do k = 1, num_nb14_calc(ij)

        i  = nb14_calc_list(1,k,ij)
        j  = nb14_calc_list(2,k,ij)
        ix = nb14_calc_list(3,k,ij)
        iy = nb14_calc_list(4,k,ij)

        ! compute distance
        !
        dij(1:3) = coord(ix,1:3,i) - coord(iy,1:3,j)
        dij(1:3) = dij(1:3) + real(cell_move(1:3,j,i),wp)*system_size(1:3)
        rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        rij2     = max(rij2, minimum_contact)

        ! energy and gradient
        !
        rij2 = cutoff2*density/rij2
        L    = int(rij2)
        R    = rij2 - L

        qq_scale  = -enefunc%nb14_qq_scale(k,ij)+1.0_wp
        term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
        cc        = charge(ix,i)*charge(iy,j)*qq_scale
        elec_temp = elec_temp + term_elec*cc

        term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
        coef      = cc*term_elec

        work(1:3) = coef*dij(1:3)
        viri(1:3) = viri(1:3) + dij(1:3)*work(1:3)

        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
        force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + work(1:3)


      end do

      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
      eelec(id+1) = eelec(id+1) + elec_temp

    end do

    !$omp end parallel

    return

  end subroutine pme_bond_corr_gro_amber_check

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_ljpme_intel
  !> @brief        calculate bond correction term in PME (general)
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[in]    atmcls  : atom class number
  !! @param[in]    coord   : coordinates for each cell
  !! @param[inout] force   : forces for each cell
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_ljpme_intel(domain, enefunc, atmcls, coord, force, &
                                       virial, evdw, eelec)

    ! formal arguments 
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    integer,                 intent(in)    :: atmcls(:)
    real(wp),                intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(:,:,:)
    real(dp),                intent(inout) :: evdw (nthread)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), viri(3), cutoff2
    real(wp)                 :: R, term_elec, term_evdw, ccr, coef
    real(wp)                 :: lj6_i, lj6_j, iqtmp, jqtmp
    real(wp)                 :: elec_temp, evdw_temp
    integer                  :: i, ix, iy, j, k, ij, L
    integer                  :: id, omp_get_thread_num
    integer                  :: ncell_local
    integer                  :: iatmcls, jatmcls

    real(wp),        pointer :: nonb_lj6_factor(:)
    real(wp),        pointer :: table_ecor(:), table_decor(:)
    real(wp),        pointer :: table_vcor(:), table_dvcor(:)
    real(wp),        pointer :: density, cutoff
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: natom(:), nwater(:), water_list(:,:,:)
    integer,         pointer :: num_nonb_excl(:), nonb_excl_list(:,:,:)


    cell_pair       => domain%cell_pairlist1
    natom           => domain%num_atom
    nwater          => domain%num_water
    water_list      => domain%water_list

    table_ecor      => enefunc%table%table_ecor
    table_decor     => enefunc%table%table_decor
    table_vcor      => enefunc%table%table_vcor
    table_dvcor     => enefunc%table%table_dvcor
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    num_nonb_excl   => enefunc%num_nonb_excl
    nonb_excl_list  => enefunc%nonb_excl_list
    nonb_lj6_factor => enefunc%nonb_lj6_factor

    ncell_local = domain%num_cell_local
    cutoff2     = cutoff*cutoff

    !$omp parallel default(shared)                                      &
    !$omp private(id, i, j, k, ij, ix, iy, dij, rij2, L, R, coef, work, &
    !$omp         term_elec,  term_evdw, ccr, elec_temp, evdw_temp,     &
    !$omp         iatmcls, jatmcls, lj6_i, lj6_j, iqtmp, jqtmp, viri)
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

      do k = 1, num_nonb_excl(ij)

        i  = nonb_excl_list(1,k,ij)
        j  = nonb_excl_list(2,k,ij)
        iatmcls = atmcls(i)
        jatmcls = atmcls(j)
        lj6_i = nonb_lj6_factor(iatmcls)
        lj6_j = nonb_lj6_factor(jatmcls)

        ! compute distance
        !
        dij(1) = coord(i,1,1) - coord(j,1,1)
        dij(2) = coord(i,2,1) - coord(j,2,1)
        dij(3) = coord(i,3,1) - coord(j,3,1)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        iqtmp  = coord(i,4,1)
        jqtmp  = coord(j,4,1)

        ! energy and gradient
        !
        rij2 = cutoff2*density/rij2
        L    = int(rij2)
        R    = rij2 - L

        term_elec = table_ecor(L) + R*(table_ecor(L+1)-table_ecor(L))
        ccr       = iqtmp * jqtmp * term_elec
        elec_temp = elec_temp + ccr

        term_evdw = table_vcor(L) + R*(table_vcor(L+1)-table_vcor(L))
        ccr       = lj6_i * lj6_j * term_evdw
        evdw_temp = evdw_temp + ccr

        term_elec = table_decor(L) + R*(table_decor(L+1)-table_decor(L))
        coef      = iqtmp * jqtmp * term_elec

        term_evdw = table_dvcor(L) + R*(table_dvcor(L+1)-table_dvcor(L))
        coef      = coef + lj6_i*lj6_j*term_evdw
        work(1)   = coef*dij(1)
        work(2)   = coef*dij(2)
        work(3)   = coef*dij(3)

        viri(1)   = viri(1) + dij(1)*work(1)
        viri(2)   = viri(2) + dij(2)*work(2)
        viri(3)   = viri(3) + dij(3)*work(3)

        force(i,1,1,id+1) = force(i,1,1,id+1) - work(1)
        force(i,2,1,id+1) = force(i,2,1,id+1) - work(2)
        force(i,3,1,id+1) = force(i,3,1,id+1) - work(3)
        force(j,1,1,id+1) = force(j,1,1,id+1) + work(1)
        force(j,2,1,id+1) = force(j,2,1,id+1) + work(2)
        force(j,3,1,id+1) = force(j,3,1,id+1) + work(3)

      end do

      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
      eelec(id+1) = eelec(id+1) + elec_temp
      evdw (id+1) = evdw (id+1) + evdw_temp

    end do

    !$omp end parallel

    return

  end subroutine pme_bond_corr_ljpme_intel

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_ljpme_fugaku
  !> @brief        calculate bond correction term in PME (general)
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[in]    atmcls  : atom class number
  !! @param[in]    coord   : coordinates for each cell
  !! @param[inout] force   : forces for each cell
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_ljpme_fugaku(domain, enefunc, atmcls, coord, force, &
                                        virial, evdw, eelec)

    ! formal arguments 
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    integer,                 intent(in)    :: atmcls(:)
    real(wp),                intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(:,:,:)
    real(dp),                intent(inout) :: evdw (nthread)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), viri(3), cutoff2
    real(wp)                 :: R, term_elec, term_evdw, ccr, coef
    real(wp)                 :: lj6_i, lj6_j, iqtmp, jqtmp
    real(wp)                 :: elec_temp, evdw_temp
    integer                  :: i, ix, iy, j, k, ij, L
    integer                  :: id, omp_get_thread_num
    integer                  :: ncell_local
    integer                  :: iatmcls, jatmcls

    real(wp),        pointer :: nonb_lj6_factor(:)
    real(wp),        pointer :: table_ecor(:), table_decor(:)
    real(wp),        pointer :: table_vcor(:), table_dvcor(:)
    real(wp),        pointer :: density, cutoff
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: natom(:), nwater(:), water_list(:,:,:)
    integer,         pointer :: num_nonb_excl(:), nonb_excl_list(:,:,:)


    cell_pair       => domain%cell_pairlist1
    natom           => domain%num_atom
    nwater          => domain%num_water
    water_list      => domain%water_list

    table_ecor      => enefunc%table%table_ecor
    table_decor     => enefunc%table%table_decor
    table_vcor      => enefunc%table%table_vcor
    table_dvcor     => enefunc%table%table_dvcor
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    num_nonb_excl   => enefunc%num_nonb_excl
    nonb_excl_list  => enefunc%nonb_excl_list
    nonb_lj6_factor => enefunc%nonb_lj6_factor

    ncell_local = domain%num_cell_local
    cutoff2     = cutoff*cutoff

    !$omp parallel default(shared)                                      &
    !$omp private(id, i, j, k, ij, ix, iy, dij, rij2, L, R, coef, work, &
    !$omp         term_elec,  term_evdw, ccr, elec_temp, evdw_temp,     &
    !$omp         iatmcls, jatmcls, lj6_i, lj6_j, iqtmp, jqtmp, viri)
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

      do k = 1, num_nonb_excl(ij)

        i  = nonb_excl_list(1,k,ij)
        j  = nonb_excl_list(2,k,ij)
        iatmcls = atmcls(i)
        jatmcls = atmcls(j)
        lj6_i = nonb_lj6_factor(iatmcls)
        lj6_j = nonb_lj6_factor(jatmcls)

        ! compute distance
        !
        dij(1) = coord(1,i,1) - coord(1,j,1)
        dij(2) = coord(2,i,1) - coord(2,j,1)
        dij(3) = coord(3,i,1) - coord(3,j,1)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        iqtmp  = coord(4,i,1)
        jqtmp  = coord(4,j,1)

        ! energy and gradient
        !
        rij2 = cutoff2*density/rij2
        L    = int(rij2)
        R    = rij2 - L

        term_elec = table_ecor(L) + R*(table_ecor(L+1)-table_ecor(L))
        ccr       = iqtmp * jqtmp * term_elec
        elec_temp = elec_temp + ccr

        term_evdw = table_vcor(L) + R*(table_vcor(L+1)-table_vcor(L))
        ccr       = lj6_i * lj6_j * term_evdw
        evdw_temp = evdw_temp + ccr

        term_elec = table_decor(L) + R*(table_decor(L+1)-table_decor(L))
        coef      = iqtmp * jqtmp * term_elec

        term_evdw = table_dvcor(L) + R*(table_dvcor(L+1)-table_dvcor(L))
        coef      = coef + lj6_i*lj6_j*term_evdw
        work(1)   = coef*dij(1)
        work(2)   = coef*dij(2)
        work(3)   = coef*dij(3)

        viri(1)   = viri(1) + dij(1)*work(1)
        viri(2)   = viri(2) + dij(2)*work(2)
        viri(3)   = viri(3) + dij(3)*work(3)

        force(1,i,1,id+1) = force(1,i,1,id+1) - work(1)
        force(2,i,1,id+1) = force(2,i,1,id+1) - work(2)
        force(3,i,1,id+1) = force(3,i,1,id+1) - work(3)
        force(1,j,1,id+1) = force(1,j,1,id+1) + work(1)
        force(2,j,1,id+1) = force(2,j,1,id+1) + work(2)
        force(3,j,1,id+1) = force(3,j,1,id+1) + work(3)

      end do

      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
      eelec(id+1) = eelec(id+1) + elec_temp
      evdw (id+1) = evdw (id+1) + evdw_temp

    end do

    !$omp end parallel

    return

  end subroutine pme_bond_corr_ljpme_fugaku

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_ljpme_gpu
  !> @brief        calculate bond correction term in PME (general)
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[in]    atmcls  : atom class number
  !! @param[in]    coord   : coordinates for each cell
  !! @param[inout] force   : forces for each cell
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_ljpme_gpu(domain, enefunc, atmcls, coord, force, &
                                     virial, evdw, eelec)

    ! formal arguments 
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    integer,                 intent(in)    :: atmcls(:)
    real(wp),                intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(:,:,:)
    real(dp),                intent(inout) :: evdw (nthread)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), viri(3), cutoff2
    real(wp)                 :: R, term_elec, term_evdw, ccr, coef
    real(wp)                 :: lj6_i, lj6_j, iqtmp, jqtmp
    real(wp)                 :: elec_temp, evdw_temp
    integer                  :: i, ix, iy, j, k, ij, L, icel, jcel
    integer                  :: id, omp_get_thread_num
    integer                  :: ncell_local, num_atom
    integer                  :: iatmcls, jatmcls

    real(wp),        pointer :: system_size(:)
    real(wp),        pointer :: nonb_lj6_factor(:)
    real(wp),        pointer :: table_ecor(:), table_decor(:)
    real(wp),        pointer :: table_vcor(:), table_dvcor(:)
    real(wp),        pointer :: density, cutoff
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: natom(:), nwater(:), water_list(:,:,:)
    integer,         pointer :: start_atom(:)
    integer,         pointer :: num_nonb_excl(:), nonb_excl_list(:,:,:)
    integer(1),      pointer :: cell_move(:,:,:)


    cell_pair       => domain%cell_pairlist1
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    nwater          => domain%num_water
    water_list      => domain%water_list
    system_size     => domain%system_size
    cell_move       => domain%cell_move

    table_ecor      => enefunc%table%table_ecor
    table_decor     => enefunc%table%table_decor
    table_vcor      => enefunc%table%table_vcor
    table_dvcor     => enefunc%table%table_dvcor
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    num_nonb_excl   => enefunc%num_nonb_excl
    nonb_excl_list  => enefunc%nonb_excl_list
    nonb_lj6_factor => enefunc%nonb_lj6_factor

    ncell_local     = domain%num_cell_local
    cutoff2         = cutoff*cutoff
    num_atom        = domain%num_atom_domain

    !$omp parallel default(shared)                                      &
    !$omp private(id, i, j, k, ij, ix, iy, dij, rij2, L, R, coef, work, &
    !$omp         term_elec,  term_evdw, ccr, elec_temp, evdw_temp,     &
    !$omp         icel, jcel, iatmcls, jatmcls, lj6_i, lj6_j,           &
    !$omp         iqtmp, jqtmp, viri)
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

      do k = 1, num_nonb_excl(ij)

        icel  = nonb_excl_list(1,k,ij)
        jcel  = nonb_excl_list(2,k,ij)
        ix    = nonb_excl_list(3,k,ij)
        iy    = nonb_excl_list(4,k,ij)
        i     = start_atom(icel) + ix
        j     = start_atom(jcel) + iy
        iatmcls = atmcls(i)
        jatmcls = atmcls(j)
        lj6_i = nonb_lj6_factor(iatmcls)
        lj6_j = nonb_lj6_factor(jatmcls)

        ! compute distance
        !
        dij(1) = coord(           i,1,1) - coord(           j,1,1)
        dij(2) = coord(  num_atom+i,1,1) - coord(  num_atom+j,1,1)
        dij(3) = coord(2*num_atom+i,1,1) - coord(2*num_atom+j,1,1)
        dij(1) = dij(1) + cell_move(1,jcel,icel)*system_size(1)
        dij(2) = dij(2) + cell_move(2,jcel,icel)*system_size(2)
        dij(3) = dij(3) + cell_move(3,jcel,icel)*system_size(3)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        iqtmp  = coord(3*num_atom+i,1,1)
        jqtmp  = coord(3*num_atom+j,1,1)

        ! energy and gradient
        !
        rij2 = cutoff2*density/rij2
        L    = int(rij2)
        R    = rij2 - L

        term_elec = table_ecor(L) + R*(table_ecor(L+1)-table_ecor(L))
        ccr       = iqtmp * jqtmp * term_elec
        elec_temp = elec_temp + ccr

        term_evdw = table_vcor(L) + R*(table_vcor(L+1)-table_vcor(L))
        ccr       = lj6_i * lj6_j * term_evdw
        evdw_temp = evdw_temp + ccr

        term_elec = table_decor(L) + R*(table_decor(L+1)-table_decor(L))
        coef      = iqtmp * jqtmp * term_elec

        term_evdw = table_dvcor(L) + R*(table_dvcor(L+1)-table_dvcor(L))
        coef      = coef + lj6_i*lj6_j*term_evdw
        work(1)   = coef*dij(1)
        work(2)   = coef*dij(2)
        work(3)   = coef*dij(3)

        viri(1)   = viri(1) + dij(1)*work(1)
        viri(2)   = viri(2) + dij(2)*work(2)
        viri(3)   = viri(3) + dij(3)*work(3)

        force(1,ix,icel,id+1) = force(1,ix,icel,id+1) - work(1)
        force(2,ix,icel,id+1) = force(2,ix,icel,id+1) - work(2)
        force(3,ix,icel,id+1) = force(3,ix,icel,id+1) - work(3)
        force(1,iy,jcel,id+1) = force(1,iy,jcel,id+1) + work(1)
        force(2,iy,jcel,id+1) = force(2,iy,jcel,id+1) + work(2)
        force(3,iy,jcel,id+1) = force(3,iy,jcel,id+1) + work(3)

      end do

      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
      eelec(id+1) = eelec(id+1) + elec_temp
      evdw (id+1) = evdw (id+1) + evdw_temp

    end do

    !$omp end parallel

    return

  end subroutine pme_bond_corr_ljpme_gpu

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_ljpme_generic
  !> @brief        calculate bond correction term in PME (general)
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[in]    coord   : coordinates for each cell
  !! @param[inout] force   : forces for each cell
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_ljpme_generic(domain, enefunc, coord, force, &
                                         virial, evdw, eelec)

    ! formal arguments 
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(:,:,:)
    real(dp),                intent(inout) :: evdw (nthread)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), viri(3), cutoff2
    real(wp)                 :: R, term_elec, term_evdw, ccr, coef
    real(wp)                 :: lj6_i, lj6_j
    real(wp)                 :: elec_temp, evdw_temp
    integer                  :: i, ix, iy, j, k, ij, L
    integer                  :: id, omp_get_thread_num
    integer                  :: ncell_local
    integer                  :: iatmcls, jatmcls

    real(wp),        pointer :: charge(:,:), system_size(:)
    real(wp),        pointer :: nonb_lj6_factor(:)
    real(wp),        pointer :: table_ecor(:), table_decor(:)
    real(wp),        pointer :: table_vcor(:), table_dvcor(:)
    real(wp),        pointer :: density, cutoff
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: natom(:), nwater(:), water_list(:,:,:)
    integer,         pointer :: num_nonb_excl(:), nonb_excl_list(:,:,:)
    integer,         pointer :: atmcls(:,:)
    integer(1),      pointer :: cell_move(:,:,:)


    cell_pair       => domain%cell_pairlist1
    natom           => domain%num_atom
    nwater          => domain%num_water
    water_list      => domain%water_list
    charge          => domain%charge
    atmcls          => domain%atom_cls_no
    system_size     => domain%system_size
    cell_move       => domain%cell_move

    table_ecor      => enefunc%table%table_ecor
    table_decor     => enefunc%table%table_decor
    table_vcor      => enefunc%table%table_vcor
    table_dvcor     => enefunc%table%table_dvcor
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    num_nonb_excl   => enefunc%num_nonb_excl
    nonb_excl_list  => enefunc%nonb_excl_list
    nonb_lj6_factor => enefunc%nonb_lj6_factor

    ncell_local = domain%num_cell_local
    cutoff2     = cutoff*cutoff

    !$omp parallel default(shared)                                      &
    !$omp private(id, i, j, k, ij, ix, iy, dij, rij2, L, R, coef, work, &
    !$omp         term_elec,  term_evdw, ccr, elec_temp, evdw_temp,     &
    !$omp         iatmcls, jatmcls, lj6_i, lj6_j, viri)
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

      do k = 1, num_nonb_excl(ij)

        i  = nonb_excl_list(1,k,ij)
        j  = nonb_excl_list(2,k,ij)
        ix = nonb_excl_list(3,k,ij)
        iy = nonb_excl_list(4,k,ij)
        iatmcls = atmcls(ix,i)
        jatmcls = atmcls(iy,j)
        lj6_i = nonb_lj6_factor(iatmcls)
        lj6_j = nonb_lj6_factor(jatmcls)

        ! compute distance
        !
        dij(1) = coord(ix,1,i) - coord(iy,1,j)
        dij(2) = coord(ix,2,i) - coord(iy,2,j)
        dij(3) = coord(ix,3,i) - coord(iy,3,j)
        dij(1) = dij(1) + real(cell_move(1,j,i),wp)*system_size(1)
        dij(2) = dij(2) + real(cell_move(2,j,i),wp)*system_size(2)
        dij(3) = dij(3) + real(cell_move(3,j,i),wp)*system_size(3)

        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        ! energy and gradient
        !
        rij2 = cutoff2*density/rij2
        L    = int(rij2)
        R    = rij2 - L

        term_elec = table_ecor(L) + R*(table_ecor(L+1)-table_ecor(L))
        ccr       = charge(ix,i) * charge(iy,j) * term_elec
        elec_temp = elec_temp + ccr

        term_evdw = table_vcor(L) + R*(table_vcor(L+1)-table_vcor(L))
        ccr       = lj6_i * lj6_j * term_evdw
        evdw_temp = evdw_temp + ccr

        term_elec = table_decor(L) + R*(table_decor(L+1)-table_decor(L))
        coef      = charge(ix,i) * charge(iy,j) * term_elec

        term_evdw = table_dvcor(L) + R*(table_dvcor(L+1)-table_dvcor(L))
        coef      = coef + lj6_i*lj6_j*term_evdw
        work(1)   = coef*dij(1)
        work(2)   = coef*dij(2)
        work(3)   = coef*dij(3)

        viri(1)   = viri(1) + dij(1)*work(1)
        viri(2)   = viri(2) + dij(2)*work(2)
        viri(3)   = viri(3) + dij(3)*work(3)

        force(1,ix,i,id+1) = force(1,ix,i,id+1) - work(1)
        force(2,ix,i,id+1) = force(2,ix,i,id+1) - work(2)
        force(3,ix,i,id+1) = force(3,ix,i,id+1) - work(3)
        force(1,iy,j,id+1) = force(1,iy,j,id+1) + work(1)
        force(2,iy,j,id+1) = force(2,iy,j,id+1) + work(2)
        force(3,iy,j,id+1) = force(3,iy,j,id+1) + work(3)

      end do

      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
      eelec(id+1) = eelec(id+1) + elec_temp
      evdw (id+1) = evdw (id+1) + evdw_temp

    end do

    !$omp end parallel

    return

  end subroutine pme_bond_corr_ljpme_generic

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_ljpme_check
  !> @brief        calculate bond correction term in PME (general)
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[in]    coord   : coordinates for each cell
  !! @param[inout] force   : forces for each cell
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_ljpme_check(domain, enefunc, coord, force, &
                                       virial, evdw, eelec)

    ! formal arguments 
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(:,:,:)
    real(dp),                intent(inout) :: evdw (nthread)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), viri(3), cutoff2
    real(wp)                 :: minimum_contact
    real(wp)                 :: R, term_elec, term_evdw, ccr, coef
    real(wp)                 :: lj6_i, lj6_j
    real(wp)                 :: elec_temp, evdw_temp
    integer                  :: i, ix, iy, j, k, ij, L
    integer                  :: id, omp_get_thread_num
    integer                  :: ncell_local
    integer                  :: iatmcls, jatmcls

    real(wp),        pointer :: charge(:,:), system_size(:)
    real(wp),        pointer :: nonb_lj6_factor(:)
    real(wp),        pointer :: table_ecor(:), table_decor(:)
    real(wp),        pointer :: table_vcor(:), table_dvcor(:)
    real(wp),        pointer :: density, cutoff
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: natom(:), nwater(:), water_list(:,:,:)
    integer,         pointer :: num_nonb_excl(:), nonb_excl_list(:,:,:)
    integer,         pointer :: atmcls(:,:)
    integer(1),      pointer :: cell_move(:,:,:)


    cell_pair       => domain%cell_pairlist1
    natom           => domain%num_atom
    nwater          => domain%num_water
    water_list      => domain%water_list
    charge          => domain%charge
    atmcls          => domain%atom_cls_no
    system_size     => domain%system_size
    cell_move       => domain%cell_move

    table_ecor      => enefunc%table%table_ecor
    table_decor     => enefunc%table%table_decor
    table_vcor      => enefunc%table%table_vcor
    table_dvcor     => enefunc%table%table_dvcor
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    num_nonb_excl   => enefunc%num_nonb_excl
    nonb_excl_list  => enefunc%nonb_excl_list
    nonb_lj6_factor => enefunc%nonb_lj6_factor

    ncell_local     = domain%num_cell_local
    cutoff2         = cutoff*cutoff
    minimum_contact = enefunc%minimum_contact

    !$omp parallel default(shared)                                      &
    !$omp private(id, i, j, k, ij, ix, iy, dij, rij2, L, R, coef, work, &
    !$omp         term_elec,  term_evdw, ccr, elec_temp, evdw_temp,     &
    !$omp         iatmcls, jatmcls, lj6_i, lj6_j, viri)
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

      do k = 1, num_nonb_excl(ij)

        i  = nonb_excl_list(1,k,ij)
        j  = nonb_excl_list(2,k,ij)
        ix = nonb_excl_list(3,k,ij)
        iy = nonb_excl_list(4,k,ij)
        iatmcls = atmcls(ix,i)
        jatmcls = atmcls(iy,j)
        lj6_i = nonb_lj6_factor(iatmcls)
        lj6_j = nonb_lj6_factor(jatmcls)

        ! compute distance
        !
        dij(1) = coord(ix,1,i) - coord(iy,1,j)
        dij(2) = coord(ix,2,i) - coord(iy,2,j)
        dij(3) = coord(ix,3,i) - coord(iy,3,j)
        dij(1) = dij(1) + real(cell_move(1,j,i),wp)*system_size(1)
        dij(2) = dij(2) + real(cell_move(2,j,i),wp)*system_size(2)
        dij(3) = dij(3) + real(cell_move(3,j,i),wp)*system_size(3)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        rij2   = max(rij2, minimum_contact)

        ! energy and gradient
        !
        rij2 = cutoff2*density/rij2
        L    = int(rij2)
        R    = rij2 - L

        term_elec = table_ecor(L) + R*(table_ecor(L+1)-table_ecor(L))
        ccr       = charge(ix,i) * charge(iy,j) * term_elec
        elec_temp = elec_temp + ccr

        term_evdw = table_vcor(L) + R*(table_vcor(L+1)-table_vcor(L))
        ccr       = lj6_i * lj6_j * term_evdw
        evdw_temp = evdw_temp + ccr

        term_elec = table_decor(L) + R*(table_decor(L+1)-table_decor(L))
        coef      = charge(ix,i) * charge(iy,j) * term_elec

        term_evdw = table_dvcor(L) + R*(table_dvcor(L+1)-table_dvcor(L))
        coef      = coef + lj6_i*lj6_j*term_evdw
        work(1)   = coef*dij(1)
        work(2)   = coef*dij(2)
        work(3)   = coef*dij(3)

        viri(1)   = viri(1) + dij(1)*work(1)
        viri(2)   = viri(2) + dij(2)*work(2)
        viri(3)   = viri(3) + dij(3)*work(3)

        force(1,ix,i,id+1) = force(1,ix,i,id+1) - work(1)
        force(2,ix,i,id+1) = force(2,ix,i,id+1) - work(2)
        force(3,ix,i,id+1) = force(3,ix,i,id+1) - work(3)
        force(1,iy,j,id+1) = force(1,iy,j,id+1) + work(1)
        force(2,iy,j,id+1) = force(2,iy,j,id+1) + work(2)
        force(3,iy,j,id+1) = force(3,iy,j,id+1) + work(3)

      end do

      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
      eelec(id+1) = eelec(id+1) + elec_temp
      evdw (id+1) = evdw (id+1) + evdw_temp

    end do

    !$omp end parallel

    return

  end subroutine pme_bond_corr_ljpme_check

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_gro_amber_generic_fep
  !> @brief        calculate bodn correction relating with 14 scaling for FEP
  !  @authors      HO
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[in]    coord   : coordinates for each cell
  !! @param[inout] force   : forces for each cell
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_gro_amber_generic_fep(domain, enefunc, coord, &
                                             force, virial, eelec)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), viri(3), cutoff2
    real(wp)                 :: R, term_elec, coef, cc, qq_scale
    real(wp)                 :: elec_temp
    integer                  :: i, ix, iy, j, k, ij, L
    integer                  :: ncell_local, id, omp_get_thread_num

    real(wp),        pointer :: charge(:,:), system_size(:)
    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: natom(:)
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: num_nb14_calc(:), nb14_calc_list(:,:,:)
    integer(1),      pointer :: cell_move(:,:,:)

    ! FEP
    real(wp),        pointer :: coord_pbc(:,:,:)
    real(wp)                  :: trans_x, trans_y, trans_z

    charge          => domain%charge
    cell_pair       => domain%cell_pairlist1
    natom           => domain%num_atom
    system_size     => domain%system_size
    cell_move       => domain%cell_move

    table_ene       => enefunc%table%table_ecor
    table_grad      => enefunc%table%table_decor
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    cutoff2         =  cutoff * cutoff
    ncell_local     =  domain%num_cell_local

    ! FEP
    num_nb14_calc   => enefunc%num_nb14_calc_fep
    nb14_calc_list  => enefunc%nb14_calc_list_fep
    coord_pbc      => domain%translated_fep

    !$omp parallel default(shared)                                      &
    !$omp private(id, i, j, k, ij, ix, iy, dij, rij2, L, R, coef, work, &
    !$omp         term_elec, cc, qq_scale, viri, elec_temp,             &
    !$omp         trans_x, trans_y, trans_z)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    ! FEP: shift coord for pbc
    do i = id+1, domain%num_cell_local + domain%num_cell_boundary, nthread
      do ix = 1, natom(i)
        coord_pbc(ix,1,i) = domain%coord(1,ix,i) + domain%trans_vec(1,ix,i)
        coord_pbc(ix,2,i) = domain%coord(2,ix,i) + domain%trans_vec(2,ix,i)
        coord_pbc(ix,3,i) = domain%coord(3,ix,i) + domain%trans_vec(3,ix,i)
      end do
    end do

    !$omp barrier

    ! calculate energy and gradient
    !
    do ij = id+1, ncell_local, nthread

      elec_temp = 0.0_wp
      viri(1:3) = 0.0_wp

      do k = 1, num_nb14_calc(ij)

        i  = nb14_calc_list(1,k,ij)
        j  = nb14_calc_list(2,k,ij)
        ix = nb14_calc_list(3,k,ij)
        iy = nb14_calc_list(4,k,ij)

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

        ! compute distance
        !
        dij(1) = coord_pbc(ix,1,i) - coord_pbc(iy,1,j)
        dij(2) = coord_pbc(ix,2,i) - coord_pbc(iy,2,j)
        dij(3) = coord_pbc(ix,3,i) - coord_pbc(iy,3,j)
        dij(1) = dij(1) + trans_x
        dij(2) = dij(2) + trans_y
        dij(3) = dij(3) + trans_z
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        ! energy and gradient
        !
        rij2 = cutoff2*density/rij2
        L    = int(rij2)
        R    = rij2 - L

        qq_scale  = -enefunc%nb14_qq_scale_fep(k,ij)+1.0_wp
        term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
        cc        = charge(ix,i)*charge(iy,j)*qq_scale
        elec_temp = elec_temp + term_elec*cc

        term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
        coef      = cc*term_elec

        work(1) = coef*dij(1)
        work(2) = coef*dij(2)
        work(3) = coef*dij(3)

        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force(1,ix,i,id+1) = force(1,ix,i,id+1) - work(1)
        force(2,ix,i,id+1) = force(2,ix,i,id+1) - work(2)
        force(3,ix,i,id+1) = force(3,ix,i,id+1) - work(3)
        force(1,iy,j,id+1) = force(1,iy,j,id+1) + work(1)
        force(2,iy,j,id+1) = force(2,iy,j,id+1) + work(2)
        force(3,iy,j,id+1) = force(3,iy,j,id+1) + work(3)

      end do

      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,1,id+1) = virial(3,3,id+1) - viri(3)
      eelec(id+1) = eelec(id+1) + elec_temp

    end do

    !$omp end parallel

    return

  end subroutine pme_bond_corr_gro_amber_generic_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_gro_amber_check_fep
  !> @brief        calculate bodn correction relating with 14 scaling for FEP
  !  @authors      HO
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[in]    coord   : coordinates for each cell
  !! @param[inout] force   : forces for each cell
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_gro_amber_check_fep(domain, enefunc, coord, &
                                           force, virial, eelec)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), viri(3), cutoff2
    real(wp)                 :: R, term_elec, coef, cc, qq_scale
    real(wp)                 :: elec_temp, minimum_contact
    integer                  :: i, ix, iy, j, k, ij, L
    integer                  :: ncell_local, id, omp_get_thread_num

    real(wp),        pointer :: charge(:,:), system_size(:)
    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: natom(:)
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: num_nb14_calc(:), nb14_calc_list(:,:,:)
    integer(1),      pointer :: cell_move(:,:,:)

    ! FEP
    real(wp),        pointer :: coord_pbc(:,:,:)
    real(wp)                 :: trans_x, trans_y, trans_z

    charge          => domain%charge
    cell_pair       => domain%cell_pairlist1
    natom           => domain%num_atom
    system_size     => domain%system_size
    cell_move       => domain%cell_move

    table_ene       => enefunc%table%table_ecor
    table_grad      => enefunc%table%table_decor
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    cutoff2         =  cutoff * cutoff
    ncell_local     =  domain%num_cell_local
    minimum_contact =  enefunc%minimum_contact

    ! FEP
    num_nb14_calc   => enefunc%num_nb14_calc_fep
    nb14_calc_list  => enefunc%nb14_calc_list_fep
    coord_pbc      => domain%translated_fep

    !$omp parallel default(shared)                                      &
    !$omp private(id, i, j, k, ij, ix, iy, dij, rij2, L, R, coef, work, &
    !$omp         term_elec, cc, qq_scale, viri, elec_temp,             &
    !$omp         trans_x, trans_y, trans_z)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    ! FEP: shift coord for pbc
    do i = id+1, domain%num_cell_local + domain%num_cell_boundary, nthread
      do ix = 1, natom(i)
        coord_pbc(ix,1,i) = domain%coord(1,ix,i) + domain%trans_vec(1,ix,i)
        coord_pbc(ix,2,i) = domain%coord(2,ix,i) + domain%trans_vec(2,ix,i)
        coord_pbc(ix,3,i) = domain%coord(3,ix,i) + domain%trans_vec(3,ix,i)
      end do
    end do

    !$omp barrier

    ! calculate energy and gradient
    !
    do ij = id+1, ncell_local, nthread

      elec_temp = 0.0_wp
      viri(1:3) = 0.0_wp

      do k = 1, num_nb14_calc(ij)

        i  = nb14_calc_list(1,k,ij)
        j  = nb14_calc_list(2,k,ij)
        ix = nb14_calc_list(3,k,ij)
        iy = nb14_calc_list(4,k,ij)

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

        ! compute distance
        !
        dij(1:3) = coord_pbc(ix,1:3,i) - coord_pbc(iy,1:3,j)
        dij(1)   = dij(1) + trans_x
        dij(2)   = dij(2) + trans_y
        dij(3)   = dij(3) + trans_z
        rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        rij2     = max(rij2, minimum_contact)

        ! energy and gradient
        !
        rij2 = cutoff2*density/rij2
        L    = int(rij2)
        R    = rij2 - L

        qq_scale  = -enefunc%nb14_qq_scale_fep(k,ij)+1.0_wp
        term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
        cc        = charge(ix,i)*charge(iy,j)*qq_scale
        elec_temp = elec_temp + term_elec*cc

        term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
        coef      = cc*term_elec

        work(1:3) = coef*dij(1:3)
        viri(1:3) = viri(1:3) + dij(1:3)*work(1:3)

        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
        force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + work(1:3)


      end do

      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
      eelec(id+1) = eelec(id+1) + elec_temp

    end do

    !$omp end parallel

    return

  end subroutine pme_bond_corr_gro_amber_check_fep

end module sp_energy_table_linear_bondcorr_mod
