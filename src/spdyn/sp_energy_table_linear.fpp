!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_table_linear_mod
!> @brief   calculate nonbonded energy with table and with linear interpolation
!! @authors Jaewoon Jung(JJ)
!  
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_table_linear_mod

  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_energy_str_mod
  use sp_domain_str_mod
  use mpi_parallel_mod
  use constants_mod
  use timers_mod

  implicit none
  private

  ! subroutines
  !
  public  :: compute_energy_nonbond14_table_linear
  public  :: compute_force_nonbond14_table_linear
  private :: compute_energy_nonbond14_table_linear_charmm
  private :: compute_energy_nonbond14_table_linear_gro_amber
  private :: compute_energy_nonbond14_table_linear_charmm_check
  private :: compute_energy_nonbond14_table_linear_gro_amber_check

contains
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_table_linear
  !> @brief        calculate nonbonded14 energy with PME and lookup table
  !!               (linear)
  !  @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_table_linear(domain, enefunc, &
                                               force, eelec, evdw)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    ! ==> Type 6/7
    if (enefunc%forcefield == ForcefieldCHARMM) then

      if (enefunc%nonb_limiter) then
        call compute_energy_nonbond14_table_linear_charmm_check(domain, &
                                            enefunc, force, eelec, evdw)
      else
        if (enefunc%vdw == VDWPME) then

          call compute_energy_nonbond14_ljpme_charmm(domain, enefunc, &
                                            force, eelec, evdw)

        else
          call compute_energy_nonbond14_table_linear_charmm(domain, enefunc, &
                                            force, eelec, evdw)
        end if
      endif

    ! ==> Type 12/13
    else ! ForcefieldAMBER, ForcefieldGROAMBER, ForcefieldGROMARTINI
    
      if (enefunc%nonb_limiter) then
        call compute_energy_nonbond14_table_linear_gro_amber_check(domain,  &
                                            enefunc, force, eelec, evdw)
      else
        call compute_energy_nonbond14_table_linear_gro_amber(domain, enefunc, &
                                            force, eelec, evdw)
      endif
    
    end if

    return

  end subroutine compute_energy_nonbond14_table_linear

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond14_table_linear
  !> @brief        calculate nonbonded14 force with PME and lookup table
  !!               (linear)
  !  @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond14_table_linear(domain, enefunc, &
                                                  force, eelec, evdw)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    ! ==> Type 6/7
    if (enefunc%forcefield == ForcefieldCHARMM) then

      if (enefunc%nonb_limiter) then
        call compute_energy_nonbond14_table_linear_charmm_check(domain, &
                                            enefunc, force, eelec, evdw)
      else
        if (enefunc%vdw == VDWPME) then

          call compute_energy_nonbond14_ljpme_charmm(domain, enefunc, &
                                            force, eelec, evdw)

        else
          call compute_force_nonbond14_table_linear_charmm(domain, enefunc, &
                                            force)
        end if
      endif

    ! ==> Type 12/13
    else ! ForcefieldAMBER, ForcefieldGROAMBER, ForcefieldGROMARTINI

      if (enefunc%nonb_limiter) then
        call compute_energy_nonbond14_table_linear_gro_amber_check(domain,  &
                                            enefunc, force, eelec, evdw)
      else
        call compute_energy_nonbond14_table_linear_gro_amber(domain, enefunc, &
                                            force, eelec, evdw)
      endif

    end if

    return

  end subroutine compute_force_nonbond14_table_linear

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond14_table_linear_charmm
  !> @brief        calculate nonbonded14 force with PME and lookup table
  !!               (linear)
  !  @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond14_table_linear_charmm(domain, enefunc, &
                                                         force)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)

    ! local variables
    real(wp)                 :: dij(1:3), rij2
    real(wp)                 :: R, table(6)
    real(wp)                 :: term_lj12, term_lj6, term_elec
    real(wp)                 :: cutoff2, grad_coef, lj12, lj6
    real(wp)                 :: work(1:3)
    integer                  :: i, ix, iy, j, k, ij, m, L, L1
    integer                  :: iatmcls, jatmcls
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: ncell_local, id, omp_get_thread_num

    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: nb14_lj12(:,:), nb14_lj6(:,:)
    real(wp),        pointer :: table_grad(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: natom(:), atmcls(:,:)
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: num_nb14_calc(:), nb14_calc_list(:,:,:)

#ifdef PKTIMER
    call timer_sta(224)
#ifdef FJ_PROF_FAPP
    call fapp_start("compute_energy_nonbond14_table_linear_charmm",224,0)
#endif
#endif

    coord           => domain%coord
    charge          => domain%charge
    cell_pair       => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom

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
    !$omp         work, table, lj12, lj6, iatmcls, jatmcls, L1)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    ! calculate energy and gradient
    !
    do ij = id+1, ncell_local, nthread

      do k = 1, num_nb14_calc(ij)

        i  = nb14_calc_list(1,k,ij) 
        j  = nb14_calc_list(2,k,ij) 
        ix = nb14_calc_list(3,k,ij) 
        iy = nb14_calc_list(4,k,ij) 
       
        iatmcls = atmcls(ix,i)
        jatmcls = atmcls(iy,j)

        ! compute distance
        !
        dij(1) = coord(1,ix,i) - coord(1,iy,j)
        dij(2) = coord(2,ix,i) - coord(2,iy,j)
        dij(3) = coord(3,ix,i) - coord(3,iy,j)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        ! energy and gradient
        !
        rij2  = cutoff2*density/rij2
        L     = int(rij2)
        R     = rij2 - L
        L1    = 3*L - 2

        lj6   = nb14_lj6 (iatmcls,jatmcls)
        lj12  = nb14_lj12(iatmcls,jatmcls)

        term_lj12 = table_grad(L1)   + R*(table_grad(L1+3)-table_grad(L1)  )
        term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
        term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
        grad_coef = term_lj12*lj12 - term_lj6*lj6 +     &
                    charge(ix,i)*charge(iy,j)*term_elec

        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)

        force(1,ix,i,id+1) = force(1,ix,i,id+1) - work(1)
        force(2,ix,i,id+1) = force(2,ix,i,id+1) - work(2)
        force(3,ix,i,id+1) = force(3,ix,i,id+1) - work(3)
        force(1,iy,j,id+1) = force(1,iy,j,id+1) + work(1)
        force(2,iy,j,id+1) = force(2,iy,j,id+1) + work(2)
        force(3,iy,j,id+1) = force(3,iy,j,id+1) + work(3)

      end do

    end do
    !$omp end parallel

#ifdef PKTIMER
    call timer_end(224)
#ifdef FJ_PROF_FAPP
    call fapp_stop("compute_energy_nonbond14_table_linear_charmm",224,0)
#endif
#endif

    return

  end subroutine compute_force_nonbond14_table_linear_charmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_table_linear_charmm
  !> @brief        calculate nonbonded14 energy with PME and lookup table
  !!               (linear)
  !  @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_table_linear_charmm(domain, enefunc, &
                                               force, eelec, evdw)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                 :: dij(1:3), rij2
    real(wp)                 :: R, table(6)
    real(wp)                 :: term_lj12, term_lj6, term_elec
    real(wp)                 :: cutoff2, grad_coef, lj12, lj6
    real(wp)                 :: work(1:3)
    real(wp)                 :: evdw_temp, elec_temp
    integer                  :: i, ix, iy, j, k, ij, m, L, L1
    integer                  :: iatmcls, jatmcls
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: ncell_local, id, omp_get_thread_num

    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: nb14_lj12(:,:), nb14_lj6(:,:)
    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: natom(:), atmcls(:,:)
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
    !$omp         work, table, lj12, lj6, iatmcls, jatmcls, L1,                &
    !$omp         evdw_temp, elec_temp)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    ! calculate energy and gradient
    !
    do ij = id+1, ncell_local, nthread

      evdw_temp = 0.0_wp
      elec_temp = 0.0_wp
     
      do k = 1, num_nb14_calc(ij)

        i  = nb14_calc_list(1,k,ij) 
        j  = nb14_calc_list(2,k,ij) 
        ix = nb14_calc_list(3,k,ij) 
        iy = nb14_calc_list(4,k,ij) 
       
        iatmcls = atmcls(ix,i)
        jatmcls = atmcls(iy,j)

        ! compute distance
        !
        dij(1) = coord(1,ix,i) - coord(1,iy,j)
        dij(2) = coord(2,ix,i) - coord(2,iy,j)
        dij(3) = coord(3,ix,i) - coord(3,iy,j)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        ! energy and gradient
        !
        rij2  = cutoff2*density/rij2
        L     = int(rij2)
        R     = rij2 - L
        L1    = 3*L - 2

        lj6   = nb14_lj6 (iatmcls,jatmcls)
        lj12  = nb14_lj12(iatmcls,jatmcls)

        term_lj12 = table_ene(L1)   + R*(table_ene(L1+3)-table_ene(L1)  )
        term_lj6  = table_ene(L1+1) + R*(table_ene(L1+4)-table_ene(L1+1))
        term_elec = table_ene(L1+2) + R*(table_ene(L1+5)-table_ene(L1+2))
        evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
        elec_temp = elec_temp + charge(ix,i)*charge(iy,j)*term_elec

        term_lj12 = table_grad(L1)   + R*(table_grad(L1+3)-table_grad(L1)  )
        term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
        term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
        grad_coef = term_lj12*lj12 - term_lj6*lj6 +     &
                    charge(ix,i)*charge(iy,j)*term_elec

        work(1:3) = grad_coef*dij(1:3)

        force(1,ix,i,id+1) = force(1,ix,i,id+1) - work(1)
        force(2,ix,i,id+1) = force(2,ix,i,id+1) - work(2)
        force(3,ix,i,id+1) = force(3,ix,i,id+1) - work(3)
        force(1,iy,j,id+1) = force(1,iy,j,id+1) + work(1)
        force(2,iy,j,id+1) = force(2,iy,j,id+1) + work(2)
        force(3,iy,j,id+1) = force(3,iy,j,id+1) + work(3)

      end do

      evdw(id+1)  = evdw(id+1) + evdw_temp
      eelec(id+1) = eelec(id+1) + elec_temp

    end do
    !$omp end parallel

    return

  end subroutine compute_energy_nonbond14_table_linear_charmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_ljpme_charmm
  !> @brief        calculate nonbonded14 energy with LJ-PME
  !!               (linear)
  !  @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_ljpme_charmm(domain, enefunc, &
                                               force, eelec, evdw)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                 :: dij(1:3), rij2, inv_r2, inv_r6, inv_r12
    real(wp)                 :: R, table(6)
    real(wp)                 :: lj6, lj12, lj6_i, lj6_j, lj6_ij
    real(wp)                 :: term_lj12, term_lj6, term_elec, term_temp
    real(wp)                 :: cutoff2, grad_coef
    real(wp)                 :: work(1:3)
    real(wp)                 :: evdw_temp, elec_temp
    integer                  :: i, ix, iy, j, k, ij, m, L, L1
    integer                  :: iatmcls, jatmcls
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: ncell_local, id, omp_get_thread_num

    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: nb14_lj12(:,:), nb14_lj6(:,:)
    real(wp),        pointer :: nonb_lj6_factor(:)
    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: natom(:), atmcls(:,:)
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
    nonb_lj6_factor => enefunc%nonb_lj6_factor
    num_nb14_calc   => enefunc%num_nb14_calc
    nb14_calc_list  => enefunc%nb14_calc_list

    cutoff2         =  cutoff * cutoff
    ncell_local     =  domain%num_cell_local

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, j, k, ij, ix, iy, m, ini_nb14, fin_nb14, num_nb14,    &
    !$omp         dij, rij2, L, R, term_lj12, term_lj6, term_elec, grad_coef,  &
    !$omp         work, table, lj12, lj6, lj6_i, lj6_j, lj6_ij, iatmcls,       &
    !$omp         jatmcls, L1, evdw_temp, elec_temp, inv_r2, inv_r6, inv_r12)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    ! calculate energy and gradient
    !
    do ij = id+1, ncell_local, nthread

      evdw_temp = 0.0_wp
      elec_temp = 0.0_wp
     
      do k = 1, num_nb14_calc(ij)

        i  = nb14_calc_list(1,k,ij) 
        j  = nb14_calc_list(2,k,ij) 
        ix = nb14_calc_list(3,k,ij) 
        iy = nb14_calc_list(4,k,ij) 
       
        iatmcls = atmcls(ix,i)
        jatmcls = atmcls(iy,j)
        lj6_i   = nonb_lj6_factor(iatmcls)
        lj6_j   = nonb_lj6_factor(jatmcls)
        lj6_ij  = lj6_i * lj6_j
        lj6     = nb14_lj6 (iatmcls,jatmcls)
        lj12    = nb14_lj12(iatmcls,jatmcls)
        lj6     = lj6 - lj6_ij

        ! compute distance
        !
        dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,j)
        rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        ! energy and gradient
        !
        inv_r2  = 1.0_wp / rij2
        inv_r6  = inv_r2 * inv_r2 * inv_r2
        inv_r12 = inv_r6 * inv_r6
        rij2    = cutoff2 * density * inv_r2
        L       = int(rij2)
        R       = rij2 - L
        L1      = 2*L - 1

        term_temp = inv_r6
        term_lj6  = table_ene(L1  ) + R*(table_ene(L1+2)-table_ene(L1  ))
        term_elec = table_ene(L1+1) + R*(table_ene(L1+3)-table_ene(L1+1))
        evdw_temp = evdw_temp &
                  + inv_r12*lj12 - term_temp*lj6 - term_lj6*lj6_ij
        elec_temp = elec_temp + charge(ix,i)*charge(iy,j)*term_elec

        term_lj12 = -12.0_wp * inv_r12 * inv_r2
        term_temp = -6.0_wp * inv_r6 * inv_r2
        term_lj6  = table_grad(L1  ) + R*(table_grad(L1+2)-table_grad(L1  ))
        term_elec = table_grad(L1+1) + R*(table_grad(L1+3)-table_grad(L1+1))
        grad_coef = term_lj12*lj12 - term_temp*lj6 - term_lj6*lj6_ij      &
                  + charge(ix,i)*charge(iy,j)*term_elec

        work(1:3) = grad_coef*dij(1:3)

        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
        force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + work(1:3)

      end do

      evdw(id+1)  = evdw(id+1) + evdw_temp
      eelec(id+1) = eelec(id+1) + elec_temp

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond14_ljpme_charmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_table_linear_gro_amber
  !> @brief        calculate nonbonded14 energy with PME and lookup table
  !!               (linear)
  !  @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_table_linear_gro_amber(domain, enefunc, &
                                                             force, eelec, evdw)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                 :: dij(1:3), rij2
    real(wp)                 :: R, table(6)
    real(wp)                 :: term_lj12, term_lj6, term_elec
    real(wp)                 :: cutoff2, grad_coef, lj12, lj6
    real(wp)                 :: work(1:3)
    real(wp)                 :: lj_scale, qq_scale, cc
    real(wp)                 :: inv_r12, inv_r121, inv_r123, inv_r126, inv_r1212
    real(wp)                 :: elec_temp, evdw_temp
    integer                  :: i, ix, iy, j, k, ij, m, L
    integer                  :: iatmcls, jatmcls
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: ncell_local, id, omp_get_thread_num

    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: nb14_lj12(:,:), nb14_lj6(:,:)
    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: natom(:), atmcls(:,:)
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
    !$omp         work, table, lj12, lj6, qq_scale, lj_scale, cc, inv_r12,     &
    !$omp         inv_r121, inv_r123, inv_r126, inv_r1212, iatmcls, jatmcls,   &
    !$omp         evdw_temp, elec_temp)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    ! calculate energy and gradient
    !
    do ij = id+1, ncell_local, nthread

      evdw_temp = 0.0_wp
      elec_temp = 0.0_wp

      do k = 1, num_nb14_calc(ij)

        i  = nb14_calc_list(1,k,ij)
        j  = nb14_calc_list(2,k,ij)
        ix = nb14_calc_list(3,k,ij)
        iy = nb14_calc_list(4,k,ij)

        iatmcls = atmcls(ix,i)
        jatmcls = atmcls(iy,j)

        ! compute distance
        !
        dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,j)
        rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        qq_scale = enefunc%nb14_qq_scale(k,ij)
        lj_scale = enefunc%nb14_lj_scale(k,ij)
        inv_r12  = 1.0_wp / rij2
        inv_r121 = sqrt(inv_r12)
        inv_r123 = inv_r12 * inv_r121
        inv_r126 = inv_r123 * inv_r123
        inv_r1212 = inv_r126 * inv_r126

        ! energy and gradient
        !
        rij2  = cutoff2*density/rij2
        L     = int(rij2)
        R     = rij2 - L
        lj6   = nb14_lj6 (iatmcls,jatmcls)
        lj12  = nb14_lj12(iatmcls,jatmcls)
        cc    = charge(ix,i)*charge(iy,j)*qq_scale

        term_lj12   = inv_r1212 
        term_lj6    = inv_r126
        term_elec   = table_ene(3*L) + R*(table_ene(3*L+3)-table_ene(3*L))
        evdw_temp   = evdw_temp  + (term_lj12*lj12 - term_lj6*lj6)*lj_scale
        elec_temp   = elec_temp  + cc*term_elec

        term_lj12   = -12.0_wp * inv_r1212 * inv_r12
        term_lj6    = -6.0_wp * inv_r126 * inv_r12
        term_elec   = table_grad(3*L) + R*(table_grad(3*L+3)-table_grad(3*L))
        grad_coef   = (term_lj12*lj12 - term_lj6*lj6)*lj_scale + cc*term_elec

        work(1:3) = grad_coef*dij(1:3)

        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
        force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + work(1:3)

      end do

      evdw(id+1)  = evdw(id+1) + evdw_temp
      eelec(id+1) = eelec(id+1) + elec_temp
 
    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond14_table_linear_gro_amber

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_table_linear_charmm_check
  !> @brief        calculate nonbonded14 energy with PME and lookup table
  !!               (linear)
  !  @authors      JJ, CK
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_table_linear_charmm_check &
                           (domain, enefunc, force, eelec, evdw)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                 :: dij(1:3), rij2
    real(wp)                 :: R, table(6)
    real(wp)                 :: term_lj12, term_lj6, term_elec
    real(wp)                 :: cutoff2, grad_coef, lj12, lj6
    real(wp)                 :: work(1:3)
    real(wp)                 :: minimum_contact
    real(wp)                 :: evdw_temp, elec_temp
    integer                  :: i, ix, iy, j, k, ij, m, L, L1
    integer                  :: iatmcls, jatmcls
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: ncell_local, id, omp_get_thread_num

    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: nb14_lj12(:,:), nb14_lj6(:,:)
    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: natom(:), atmcls(:,:)
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
    minimum_contact = enefunc%minimum_contact

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, j, k, ij, ix, iy, m, ini_nb14, fin_nb14, num_nb14,    &
    !$omp         dij, rij2, L, R, term_lj12, term_lj6, term_elec, grad_coef,  &
    !$omp         work, table, lj12, lj6, iatmcls, jatmcls, evdw_temp, elec_temp)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    ! calculate energy and gradient
    !
    do ij = id+1, ncell_local, nthread

      evdw_temp = 0.0_wp
      elec_temp = 0.0_wp
     
      do k = 1, num_nb14_calc(ij)

        i  = nb14_calc_list(1,k,ij) 
        j  = nb14_calc_list(2,k,ij) 
        ix = nb14_calc_list(3,k,ij) 
        iy = nb14_calc_list(4,k,ij) 
       
        iatmcls = atmcls(ix,i)
        jatmcls = atmcls(iy,j)

        ! compute distance
        !
        dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,j)
        rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        rij2     = max(rij2, minimum_contact)

        ! energy and gradient
        !
        rij2  = cutoff2*density/rij2
        L     = int(rij2)
        R     = rij2 - L
        L1    = 3*L - 2

        lj6   = nb14_lj6 (iatmcls,jatmcls)
        lj12  = nb14_lj12(iatmcls,jatmcls)

        term_lj12 = table_ene(L1)   + R*(table_ene(L1+3)-table_ene(L1)  )
        term_lj6  = table_ene(L1+1) + R*(table_ene(L1+4)-table_ene(L1+1))
        term_elec = table_ene(L1+2) + R*(table_ene(L1+5)-table_ene(L1+2))
        evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
        elec_temp = elec_temp +  charge(ix,i)*charge(iy,j)*term_elec

        term_lj12 = table_grad(L1)   + R*(table_grad(L1+3)-table_grad(L1)  )
        term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
        term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
        grad_coef = term_lj12*lj12 - term_lj6*lj6 +     &
                    charge(ix,i)*charge(iy,j)*term_elec

        work(1:3) = grad_coef*dij(1:3)

        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
        force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + work(1:3)

      end do

      evdw(id+1)  = evdw(id+1) + evdw_temp
      eelec(id+1) = eelec(id+1) + elec_temp

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond14_table_linear_charmm_check

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_table_linear_gro_amber_check
  !> @brief        calculate nonbonded14 energy with PME and lookup table
  !!               (linear)
  !  @authors      JJ, CK
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_table_linear_gro_amber_check &
                                   (domain, enefunc, force, eelec, evdw)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                 :: dij(1:3), rij2
    real(wp)                 :: R, table(6)
    real(wp)                 :: term_lj12, term_lj6, term_elec
    real(wp)                 :: cutoff2, grad_coef, lj12, lj6
    real(wp)                 :: work(1:3)
    real(wp)                 :: lj_scale, qq_scale, cc
    real(wp)                 :: inv_r12, inv_r121, inv_r123, inv_r126, inv_r1212
    real(wp)                 :: minimum_contact
    real(wp)                 :: evdw_temp, elec_temp
    integer                  :: i, ix, iy, j, k, ij, m, L
    integer                  :: iatmcls, jatmcls
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: ncell_local, id, omp_get_thread_num

    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: nb14_lj12(:,:), nb14_lj6(:,:)
    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: natom(:), atmcls(:,:)
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
    minimum_contact = enefunc%minimum_contact


    !$omp parallel default(shared)                                             &
    !$omp private(id, i, j, k, ij, ix, iy, m, ini_nb14, fin_nb14, num_nb14,    &
    !$omp         dij, rij2, L, R, term_lj12, term_lj6, term_elec, grad_coef,  &
    !$omp         work, table, lj12, lj6, qq_scale, lj_scale, cc, inv_r12,     &
    !$omp         inv_r121, inv_r123, inv_r126, inv_r1212, iatmcls, jatmcls,   &
    !$omp         evdw_temp, elec_temp)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    ! calculate energy and gradient
    !
    do ij = id+1, ncell_local, nthread

      evdw_temp = 0.0_wp
      elec_temp = 0.0_wp

      do k = 1, num_nb14_calc(ij)

        i  = nb14_calc_list(1,k,ij)
        j  = nb14_calc_list(2,k,ij)
        ix = nb14_calc_list(3,k,ij)
        iy = nb14_calc_list(4,k,ij)

        iatmcls = atmcls(ix,i)
        jatmcls = atmcls(iy,j)

        ! compute distance
        !
        dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,j)
        rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        rij2     = max(rij2, minimum_contact)
        qq_scale = enefunc%nb14_qq_scale(k,ij)
        lj_scale = enefunc%nb14_lj_scale(k,ij)
        inv_r12  = 1.0_wp / rij2
        inv_r121 = sqrt(inv_r12)
        inv_r123 = inv_r12 * inv_r121
        inv_r126 = inv_r123 * inv_r123
        inv_r1212 = inv_r126 * inv_r126

        ! energy and gradient
        !
        rij2  = cutoff2*density/rij2
        L     = int(rij2)
        R     = rij2 - L
        lj6   = nb14_lj6 (iatmcls,jatmcls)
        lj12  = nb14_lj12(iatmcls,jatmcls)
        cc    = charge(ix,i)*charge(iy,j)*qq_scale

        term_lj12   = inv_r1212 
        term_lj6    = inv_r126
        term_elec   = table_ene(3*L) + R*(table_ene(3*L+3)-table_ene(3*L))
        evdw_temp   = evdw_temp  + (term_lj12*lj12 - term_lj6*lj6)*lj_scale
        elec_temp   = elec_temp  + cc*term_elec

        term_lj12   = -12.0_wp * inv_r1212 * inv_r12
        term_lj6    = -6.0_wp * inv_r126 * inv_r12
        term_elec   = table_grad(3*L) + R*(table_grad(3*L+3)-table_grad(3*L))
        grad_coef   = (term_lj12*lj12 - term_lj6*lj6)*lj_scale + cc*term_elec

        work(1:3) = grad_coef*dij(1:3)

        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
        force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + work(1:3)

      end do

      evdw(id+1)  = evdw(id+1) + evdw_temp
      eelec(id+1) = eelec(id+1) + elec_temp
 
    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond14_table_linear_gro_amber_check


end module sp_energy_table_linear_mod
