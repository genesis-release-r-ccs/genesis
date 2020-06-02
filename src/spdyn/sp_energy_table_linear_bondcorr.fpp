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
  public   :: pme_bond_corr_linear_force
  private  :: pme_bond_corr_linear_general
  private  :: pme_bond_corr_linear_gro_amber
  private  :: pme_bond_corr_linear_general_check
  private  :: pme_bond_corr_linear_gro_amber_check

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_linear
  !> @brief        calculate pme bond correction with linear lookup table
  !  @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_linear(domain, enefunc, force, evdw, eelec)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: evdw (nthread)
    real(dp),                intent(inout) :: eelec(nthread)

    if (enefunc%nonb_limiter) then
      call pme_bond_corr_linear_general_check(domain, enefunc, force, eelec)
     
      if (enefunc%forcefield /= ForcefieldCHARMM) then
     
        call pme_bond_corr_linear_gro_amber_check(domain, enefunc, force, eelec)
     
      end if
    else
   
      if (enefunc%vdw == VDWPME) then 
        call pme_bond_corr_linear_general_lj(domain, enefunc, force, evdw, &
                                             eelec)
      else
        call pme_bond_corr_linear_general(domain, enefunc, force, eelec)
      end if
     
      if (enefunc%forcefield /= ForcefieldCHARMM) then
     
        call pme_bond_corr_linear_gro_amber(domain, enefunc, force, eelec)
     
      end if
    endif

    return

   end subroutine pme_bond_corr_linear

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_linear_force
  !> @brief        calculate pme bond correction force with linear lookup table
  !  @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_linear_force(domain, enefunc, force, evdw, eelec)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: evdw (nthread)
    real(dp),                intent(inout) :: eelec(nthread)

    if (enefunc%nonb_limiter) then
      call pme_bond_corr_linear_general_check(domain, enefunc, force, eelec)

      if (enefunc%forcefield /= ForcefieldCHARMM) then

        call pme_bond_corr_linear_gro_amber_check(domain, enefunc, force, eelec)

      end if
    else

      if (enefunc%vdw == VDWPME) then
        call pme_bond_corr_linear_general_lj(domain, enefunc, force, evdw, &
                                             eelec)
      else
        call pme_bond_corr_linear_general_force(domain, enefunc, force)
      end if

      if (enefunc%forcefield /= ForcefieldCHARMM) then

        call pme_bond_corr_linear_gro_amber(domain, enefunc, force, eelec)

      end if
    endif

    return

   end subroutine pme_bond_corr_linear_force

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_linear_general
  !> @brief        calculate bond correction term in PME (general)
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[inout] force   : forces of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_linear_general(domain, enefunc, force, eelec)

    ! formal arguments 
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), cutoff2
    real(wp)                 :: R, term_elec, ccr, coef, elec_temp
    integer                  :: i, ix, iy, j, k, ij, L
    integer                  :: id, omp_get_thread_num
    integer                  :: ncell_local

    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: table_ecor(:), table_decor(:)
    real(wp),        pointer :: density, cutoff
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: natom(:), nwater(:), water_list(:,:,:)
    integer,         pointer :: num_nonb_excl(:), nonb_excl_list(:,:,:)


    cell_pair       => domain%cell_pairlist1
    natom           => domain%num_atom
    nwater          => domain%num_water
    water_list      => domain%water_list
    coord           => domain%coord
    charge          => domain%charge

    table_ecor      => enefunc%table%table_ecor
    table_decor     => enefunc%table%table_decor
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    num_nonb_excl   => enefunc%num_nonb_excl
    nonb_excl_list  => enefunc%nonb_excl_list

    ncell_local = domain%num_cell_local
    cutoff2     = cutoff*cutoff

#ifdef PKTIMER
    call timer_sta(220)
#ifdef FJ_PROF_FAPP
    call fapp_start("pme_bond_corr_linear_general",220,0)
#endif
#endif

    !$omp parallel default(shared)                                      &
    !$omp private(id, i, j, k, ij, ix, iy, dij, rij2, L, R, coef, work, &
    !$omp         term_elec,  ccr, elec_temp)
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

      do k = 1, num_nonb_excl(ij)

        i  = nonb_excl_list(1,k,ij)
        j  = nonb_excl_list(2,k,ij)
        ix = nonb_excl_list(3,k,ij)
        iy = nonb_excl_list(4,k,ij)

        ! compute distance
        !
        dij(1) = coord(1,ix,i) - coord(1,iy,j)
        dij(2) = coord(2,ix,i) - coord(2,iy,j)
        dij(3) = coord(3,ix,i) - coord(3,iy,j)
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

        work(1:3) = coef*dij(1:3)

        force(1,ix,i,id+1) = force(1,ix,i,id+1) - work(1)
        force(2,ix,i,id+1) = force(2,ix,i,id+1) - work(2)
        force(3,ix,i,id+1) = force(3,ix,i,id+1) - work(3)
        force(1,iy,j,id+1) = force(1,iy,j,id+1) + work(1)
        force(2,iy,j,id+1) = force(2,iy,j,id+1) + work(2)
        force(3,iy,j,id+1) = force(3,iy,j,id+1) + work(3)

      end do
 
      eelec(id+1) = eelec(id+1) + elec_temp  

    end do

    !$omp end parallel

#ifdef PKTIMER
    call timer_end(220)
#ifdef FJ_PROF_FAPP
    call fapp_stop("pme_bond_corr_linear_general",220,0)
#endif
#endif

    return

  end subroutine pme_bond_corr_linear_general

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_linear_general_force
  !> @brief        calculate bond correction term in PME (general)
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[inout] force   : forces of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_linear_general_force(domain, enefunc, force)

    ! formal arguments 
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), cutoff2
    real(wp)                 :: R, term_elec, ccr, coef
    integer                  :: i, ix, iy, j, k, ij, L
    integer                  :: id, omp_get_thread_num
    integer                  :: ncell_local

    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: table_decor(:)
    real(wp),        pointer :: density, cutoff
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: natom(:), nwater(:), water_list(:,:,:)
    integer,         pointer :: num_nonb_excl(:), nonb_excl_list(:,:,:)


    cell_pair       => domain%cell_pairlist1
    natom           => domain%num_atom
    nwater          => domain%num_water
    water_list      => domain%water_list
    coord           => domain%coord
    charge          => domain%charge

    table_decor     => enefunc%table%table_decor
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    num_nonb_excl   => enefunc%num_nonb_excl
    nonb_excl_list  => enefunc%nonb_excl_list

    ncell_local = domain%num_cell_local
    cutoff2     = cutoff*cutoff

    !$omp parallel default(shared)                                      &
    !$omp private(id, i, j, k, ij, ix, iy, dij, rij2, L, R, coef, work, &
    !$omp         term_elec,  ccr)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    ! calculate energy and gradient
    !
    do ij = id+1, ncell_local, nthread

      do k = 1, num_nonb_excl(ij)

        i  = nonb_excl_list(1,k,ij)
        j  = nonb_excl_list(2,k,ij)
        ix = nonb_excl_list(3,k,ij)
        iy = nonb_excl_list(4,k,ij)

        ! compute distance
        !
        dij(1) = coord(1,ix,i) - coord(1,iy,j)
        dij(2) = coord(2,ix,i) - coord(2,iy,j)
        dij(3) = coord(3,ix,i) - coord(3,iy,j)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        ! energy and gradient
        !
        rij2 = cutoff2*density/rij2
        L    = int(rij2)
        R    = rij2 - L

        term_elec = table_decor(L) + R*(table_decor(L+1)-table_decor(L))
        coef      = charge(ix,i) * charge(iy,j) * term_elec

        work(1) = coef*dij(1)
        work(2) = coef*dij(2)
        work(3) = coef*dij(3)

        force(1,ix,i,id+1) = force(1,ix,i,id+1) - work(1)
        force(2,ix,i,id+1) = force(2,ix,i,id+1) - work(2)
        force(3,ix,i,id+1) = force(3,ix,i,id+1) - work(3)
        force(1,iy,j,id+1) = force(1,iy,j,id+1) + work(1)
        force(2,iy,j,id+1) = force(2,iy,j,id+1) + work(2)
        force(3,iy,j,id+1) = force(3,iy,j,id+1) + work(3)

      end do
 
    end do

    !$omp end parallel

    return

  end subroutine pme_bond_corr_linear_general_force


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_linear_general_lj
  !> @brief        calculate bond correction term in PME (general)
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[inout] force   : forces of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_linear_general_lj(domain, enefunc, force, evdw, &
                                             eelec)

    ! formal arguments 
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: evdw (nthread)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), cutoff2
    real(wp)                 :: R, term_elec, term_evdw, ccr, coef
    real(wp)                 :: lj6_i, lj6_j
    real(wp)                 :: elec_temp, evdw_temp
    integer                  :: i, ix, iy, j, k, ij, L
    integer                  :: id, omp_get_thread_num
    integer                  :: ncell_local
    integer                  :: iatmcls, jatmcls

    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: nonb_lj6_factor(:)
    real(wp),        pointer :: table_ecor(:), table_decor(:)
    real(wp),        pointer :: table_vcor(:), table_dvcor(:)
    real(wp),        pointer :: density, cutoff
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: natom(:), nwater(:), water_list(:,:,:)
    integer,         pointer :: num_nonb_excl(:), nonb_excl_list(:,:,:)
    integer,         pointer :: atmcls(:,:)


    cell_pair       => domain%cell_pairlist1
    natom           => domain%num_atom
    nwater          => domain%num_water
    water_list      => domain%water_list
    coord           => domain%coord
    charge          => domain%charge
    atmcls          => domain%atom_cls_no

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
    !$omp         iatmcls, jatmcls, lj6_i, lj6_j)
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
        dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,j)
        rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

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
        work(1:3) = coef*dij(1:3)

        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
        force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + work(1:3)

      end do

      eelec(id+1) = eelec(id+1) + elec_temp  
      evdw (id+1) = evdw (id+1) + evdw_temp  

    end do

    !$omp end parallel

    return

  end subroutine pme_bond_corr_linear_general_lj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_linear_gro_amber
  !> @brief        calculate bodn correction relating with 14 scaling
  !  @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_linear_gro_amber(domain, enefunc, &
                                          force, eelec)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), cutoff2
    real(wp)                 :: R, term_elec, coef, cc, qq_scale
    real(wp)                 :: elec_temp
    integer                  :: i, ix, iy, j, k, ij, L
    integer                  :: ncell_local, id, omp_get_thread_num

    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: natom(:)
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: num_nb14_calc(:), nb14_calc_list(:,:,:)


    coord           => domain%coord
    charge          => domain%charge
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
    !$omp         term_elec, cc, qq_scale, elec_temp)
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

        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
        force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + work(1:3)

      end do

      eelec(id+1) = eelec(id+1) + elec_temp

    end do

    !$omp end parallel

    return

  end subroutine pme_bond_corr_linear_gro_amber

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_linear_general_check
  !> @brief        calculate bond correction term in PME (general)
  !! @authors      JJ, CK
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[inout] force   : forces of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_linear_general_check(domain, enefunc, force, eelec)

    ! formal arguments 
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), cutoff2
    real(wp)                 :: R, term_elec, ccr, coef, elec_temp
    real(wp)                 :: minimum_contact
    integer                  :: i, ix, iy, j, k, ij, L
    integer                  :: id, omp_get_thread_num
    integer                  :: ncell_local

    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: table_ecor(:), table_decor(:)
    real(wp),        pointer :: density, cutoff
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: natom(:), nwater(:), water_list(:,:,:)
    integer,         pointer :: num_nonb_excl(:), nonb_excl_list(:,:,:)


    cell_pair       => domain%cell_pairlist1
    natom           => domain%num_atom
    nwater          => domain%num_water
    water_list      => domain%water_list
    coord           => domain%coord
    charge          => domain%charge

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
    !$omp         term_elec,  ccr, elec_temp)
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

      do k = 1, num_nonb_excl(ij)

        i  = nonb_excl_list(1,k,ij)
        j  = nonb_excl_list(2,k,ij)
        ix = nonb_excl_list(3,k,ij)
        iy = nonb_excl_list(4,k,ij)

        ! compute distance
        !
        dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,j)
        rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        rij2     = max(rij2, minimum_contact)

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

        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
        force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + work(1:3)

      end do
      eelec(id+1) = eelec(id+1) + elec_temp

    end do

    !$omp end parallel

    return

  end subroutine pme_bond_corr_linear_general_check

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_linear_gro_amber_check
  !> @brief        calculate bodn correction relating with 14 scaling
  !  @authors      JJ, CK
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_linear_gro_amber_check(domain, enefunc, &
                                                  force, eelec)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), cutoff2
    real(wp)                 :: R, term_elec, coef, cc, qq_scale
    real(wp)                 :: elec_temp, minimum_contact
    integer                  :: i, ix, iy, j, k, ij, L
    integer                  :: ncell_local, id, omp_get_thread_num

    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: natom(:)
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: num_nb14_calc(:), nb14_calc_list(:,:,:)


    coord           => domain%coord
    charge          => domain%charge
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
    minimum_contact =  enefunc%minimum_contact

    !$omp parallel default(shared)                                      &
    !$omp private(id, i, j, k, ij, ix, iy, dij, rij2, L, R, coef, work, &
    !$omp         term_elec, cc, qq_scale, elec_temp)
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

      do k = 1, num_nb14_calc(ij)

        i  = nb14_calc_list(1,k,ij)
        j  = nb14_calc_list(2,k,ij)
        ix = nb14_calc_list(3,k,ij)
        iy = nb14_calc_list(4,k,ij)

        ! compute distance
        !
        dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,j)
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

        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
        force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + work(1:3)


      end do

      eelec(id+1) = eelec(id+1) + elec_temp

    end do

    !$omp end parallel

    return

  end subroutine pme_bond_corr_linear_gro_amber_check

end module sp_energy_table_linear_bondcorr_mod
