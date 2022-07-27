!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_energy_eef1_mod
!> @brief   calculate nonbonded energy with EEF1/IMM1
!! @authors Takaharu Mori (TM)
! 
!  (c) Copyright 2016 RIKEN. All rights reserved.
! 
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_energy_eef1_mod

  use at_boundary_str_mod
  use at_pairlist_str_mod
  use at_enefunc_str_mod
  use molecules_str_mod
  use timers_mod
  use mpi_parallel_mod
  use messages_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  integer,               save :: istart, iend
  real(wp), allocatable, save :: FF(:)
  real(wp), allocatable, save :: dmfac(:,:)
  real(wp), allocatable, save :: coord_min(:,:), coord_tmp(:,:)

  real(wp), parameter         :: EEF1_CUTOFF_DIST = 9.0_wp

  ! subroutines
  public   :: compute_energy_nonbond_eef1
  private  :: compute_energy_reference
  private  :: compute_energy_nonbond14_eef1
  private  :: compute_energy_nonbond15_eef1
  private  :: compute_energy_nonbond14_imm1
  private  :: compute_energy_nonbond15_imm1
  private  :: compute_ellipsoid_depth
  private  :: g

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_eef1
  !> @brief        Calculate nonbond energy in no boundary condition
  !! @authors      TM
  !! @note         Shift function is not used in the ELECT term here.
  !!               Therefore, the ELECT energy differs from
  !!               that in the case of implicit_solvent = NONE.
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    nonb_ene : flag for calculate nonbonded energy
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !! @param[inout] esolv    : solvation free energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_eef1(enefunc, molecule, pairlist, nonb_ene,&
                                         coord, force, virial, eelec, evdw,    &
                                         esolv)

    ! formal arguments
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_molecule),        intent(in)    :: molecule
    type(s_pairlist),        intent(in)    :: pairlist
    logical,                 intent(in)    :: nonb_ene
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: eelec
    real(wp),                intent(inout) :: evdw
    real(wp),                intent(inout) :: esolv


    call timer(TimerNonBond, TimerOn)

    call timer(TimerPmeReal, TimerOn)


    call compute_energy_reference(enefunc, molecule, coord, force, esolv)

    if (.not. enefunc%imm1_use) then

      call compute_energy_nonbond14_eef1(enefunc, molecule,    &
                                         coord, force, virial, &
                                         eelec, evdw, esolv)

      call compute_energy_nonbond15_eef1(enefunc, molecule, pairlist, &
                                         coord, force, virial,        &
                                         eelec, evdw, esolv)

    else if (enefunc%imm1_use) then

      call compute_energy_nonbond14_imm1(enefunc, molecule,    &
                                         coord, force, virial, &
                                         eelec, evdw, esolv)
  
      call compute_energy_nonbond15_imm1(enefunc, molecule, pairlist, &
                                         coord, force, virial,        &
                                         eelec, evdw, esolv)

    end if


    call timer(TimerPmeReal, TimerOff)

    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_nonbond_eef1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_reference
  !> @brief        calculate self energy term
  !  @authors      TM
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !! @note         T.Lazaridis & M.Karplus, Proteins, 35, 133−152 (1999)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_reference(enefunc, molecule, coord, force, esolv)

    ! formal argumens
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_molecule), target, intent(in)    :: molecule
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: esolv

    ! local variables
    integer                   :: i, natom, n
    integer                   :: istart, iend, isd, id, omp_get_thread_num
    real(wp)                  :: ihalf_thick, zt, fz, hr
    real(wp)                  :: force_i, ddg, tmp, x2y2, pore_rad
    real(wp)                  :: rt, fr, x2y2z2, x1, y1, z1, x, y, z
    real(wp)                  :: a, b, c, m, s, m1, m2, fact, depth
    real(wp)                  :: ia, ib, ic, x2, y2, z2
    real(wp)                  :: nx, ny, nz, nn, d, x0, y0, z0
    real(wp)                  :: gg, grad_theta, grad_phi, diff, prev_d
    real(wp)                  :: r, theta, phi, theta1, theta2, f1, f2, phi1, phi2
    real(wp)                  :: s1, s2, t, s1_1, s2_1, t_1
    real(wp)                  :: fact_x, fact_y, fact_z, fact_xy
    logical                   :: make_pore
    integer,          pointer :: atmcls(:)
    real(wp),         pointer :: gref(:,:)


    natom       =  molecule%num_atoms
    atmcls      => molecule%atom_cls_no
    gref        => enefunc%eef1%gref_t
    ihalf_thick =  1.0_wp/(enefunc%eef1%imm1_memb_thick*0.5_wp)
    n           =  enefunc%eef1%imm1_exponent_n
    make_pore   =  enefunc%eef1%imm1_make_pore
    pore_rad    =  enefunc%eef1%imm1_pore_radius
    a           =  enefunc%eef1%imic_axis_a
    b           =  enefunc%eef1%imic_axis_b
    c           =  enefunc%eef1%imic_axis_c
    m1          =  enefunc%eef1%imic_exponent_m1
    m2          =  enefunc%eef1%imic_exponent_m2
    s           =  enefunc%eef1%imic_steepness
    istart      =  enefunc%eef1%istart
    iend        =  enefunc%eef1%iend

    esolv  = 0.0_wp
    ia     = 1.0_wp/a
    ib     = 1.0_wp/b
    ic     = 1.0_wp/c
    s1     = 2.0_wp/m1
    s2     = 2.0_wp/m2
    t      = m2/m1
    s1_1   = s1 - 1.0_wp
    s2_1   = s2 - 1.0_wp
    t_1    = t  - 1.0_wp
    fact_x = 2.0_wp/(a*m1)
    fact_y = 2.0_wp/(b*m1)
    fact_z = 2.0_wp/(c*m1)

    if (.not. enefunc%imm1_use) then

      !$omp parallel do       &
      !$omp private(i)        &
      !$omp reduction(+:esolv) 
      !
      do i = istart, iend
        esolv = esolv + gref(2,atmcls(i))
      end do
      !$omp end parallel do

    else if (enefunc%imm1_use) then

      if (.not. allocated(FF)) allocate(FF(natom), dmfac(3,natom))

      if (enefunc%imic_use) then

        ! implicit micelle model
        !
        call compute_ellipsoid_depth(a, b, c, m1, m2, natom, coord)

        !$omp parallel do default (none)                                       &
        !$omp private(i, fact, x, y, z, x1, y1, z1, x2, y2, z2, x0, y0, z0,    &
        !$omp         fact_xy, nx, ny, nz, nn, depth, gg)                      &
        !$omp shared (coord, coord_min, FF, dmfac, m1, ia, ib, ic, natom,      &
        !$omp         s, s1, s2, t, s1_1, s2_1, t_1, fact_x, fact_y, fact_z)
        !
        do i = 1, natom
          depth = sqrt((coord(1,i) - coord_min(1,i))**2 &
                     + (coord(2,i) - coord_min(2,i))**2 &
                     + (coord(3,i) - coord_min(3,i))**2)
          depth = s*depth

          x  = abs(coord(1,i))*ia
          y  = abs(coord(2,i))*ib
          z  = abs(coord(3,i))*ic
          gg = (x**s2 + y**s2)**t + z**s1

          if (gg < 1) depth = - depth
          FF(i) = 0.5_wp*(tanh(depth) + 1.0_wp)

          x0 = abs(coord_min(1,i))*ia
          y0 = abs(coord_min(2,i))*ib
          z0 = abs(coord_min(3,i))*ic
          x1 = x0**s2_1
          x2 = x1*x0
          y1 = y0**s2_1
          y2 = y1*y0
          z1 = z0**s1_1

          fact_xy = (x2 + y2)**t_1

          nx = sign(1.0_wp,coord_min(1,i)) * fact_x * x1 * fact_xy
          ny = sign(1.0_wp,coord_min(2,i)) * fact_y * y1 * fact_xy
          nz = sign(1.0_wp,coord_min(3,i)) * fact_z * z1
          nn = 1.0_wp/sqrt(nx*nx + ny*ny + nz*nz)

          fact = 0.5_wp*s/(cosh(depth))**2
          dmfac(1,i) = fact*nx*nn
          dmfac(2,i) = fact*ny*nn
          dmfac(3,i) = fact*nz*nn
        end do
        !$omp end parallel do

      else

        ! implicit membrane model

        if (make_pore) then
          if (abs(pore_rad - 0.0_wp) < EPS) make_pore = .false.
        end if

        !$omp parallel do                            &
        !$omp private(i, zt, fz, hr, ddg, x2y2, tmp)
        !
        do i = 1, natom

          zt  = abs(coord(3,i))*ihalf_thick
          fz  = zt**n / (1.0_wp + zt**n)
          ddg = gref(2,atmcls(i)) - gref(1,atmcls(i))

          if (make_pore) then
            x2y2  = coord(1,i)**2 + coord(2,i)**2
            tmp   = sqrt(x2y2)/pore_rad
            hr    = 1.0_wp - tmp**n / (1.0_wp + tmp**n)
            FF(i) = fz + hr - fz*hr

            if (abs(x2y2 - 0.0_wp) > EPS) then
              dmfac(1:2,i) = - dble(n)*hr*(1.0_wp-FF(i))*coord(1:2,i)/x2y2
            else
              dmfac(1:2,i) = 0.0_wp
            end if
          else
            FF(i) = fz
            dmfac(1:2,i) = 0.0_wp
          end if

          if (abs(coord(3,i) - 0.0_wp) > EPS) then
            dmfac(3,i) = dble(n)*fz*(1.0_wp-FF(i))/coord(3,i)
          else
            dmfac(3,i) = 0.0_wp
          end if

        end do
        !$omp end parallel do

      end if

      !$omp parallel                              &
      !$omp private(id, i, ddg, force_i)          &
      !$omp reduction(+:esolv) 
      !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
      do i = istart+id, iend, nthread
        ddg   = gref(2,atmcls(i)) - gref(1,atmcls(i))
        esolv = esolv + (gref(1,atmcls(i)) + FF(i)*ddg)
        force(1:3,i,id+1) = force(1:3,i,id+1) - dmfac(1:3,i)*ddg
      end do
      !$omp end parallel 

    end if

    return

  end subroutine compute_energy_reference

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_eef1
  !> @brief        calculate nonbonded14 energy
  !  @authors      TM
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !! @note         T.Lazaridis & M.Karplus, Proteins, 35, 133−152 (1999)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_eef1(enefunc, molecule, coord, force, &
                                           virial, eelec, evdw, esolv)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_molecule), target, intent(in)    :: molecule
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eelec
    real(wp),                 intent(inout) :: evdw
    real(wp),                 intent(inout) :: esolv

    ! local variables
    real(wp)                  :: dij(1:3), rij, rij2, inv_rij2
    real(wp)                  :: inv_rij, inv_r3, inv_r6, inv_r10, inv_r12
    real(wp)                  :: lj6, lj12, E14FAC
    real(wp)                  :: term_elec, shift, dshift, switch, dswitch
    real(wp)                  :: cutoff2, switchdist2
    real(wp)                  :: c1_switch, c2_switch, c4_switch, coef
    real(wp)                  :: c12, c6, term_lj12, term_lj6
    real(wp)                  :: factor_i, factor_j, alpha_4pi_i
    real(wp)                  :: inv_lambda_i, inv_lambda_j, x_i, x_j
    real(wp)                  :: work(1:3), force_local(1:3)
    real(wp)                  :: elec_i, esolv_tmp, rtmp(1:3), rvdw_i, vol_i
    integer                   :: i, j, k, l, natom, id, my_id
    integer                   :: num_nb14, ini_nb14, fin_nb14
    integer                   :: atmcls_i, atmcls_j

    real(wp),         pointer :: charge(:)
    real(wp),         pointer :: dielec_const
    real(wp),         pointer :: cutoff, switchdist
    real(wp),         pointer :: alpha_4pi(:,:)
    real(wp),         pointer :: rvdw(:), vol(:,:), inv_lambda(:,:)
    integer,          pointer :: atmcls(:)
    integer,          pointer :: num_nb14_calc(:), nb14_calc_list(:,:)

#ifdef OMP
    integer                   :: omp_get_thread_num
#endif

    charge         => molecule%charge
    atmcls         => molecule%atom_cls_no
    dielec_const   => enefunc%dielec_const
    cutoff         => enefunc%cutoffdist
    switchdist     => enefunc%switchdist
    num_nb14_calc  => enefunc%num_nb14_calc
    nb14_calc_list => enefunc%nb14_calc_list
    natom          =  molecule%num_atoms
    cutoff2        =  cutoff * cutoff
    switchdist2    =  switchdist * switchdist
    alpha_4pi      => enefunc%eef1%alpha_4pi
    rvdw           => enefunc%eef1%rvdw
    vol            => enefunc%eef1%volume
    inv_lambda     => enefunc%eef1%inv_lambda

    c1_switch = 1.0_wp/(cutoff2-switchdist2)
    c1_switch = c1_switch*c1_switch*c1_switch

    if (enefunc%forcefield == ForcefieldCHARMM19) then
      E14FAC   = 0.4_wp
    else
      E14FAC   = 1.0_wp
    end if

    ! calculate energy and gradient
    !
    num_nb14  = 0
    esolv_tmp = 0.0_wp

    !$omp parallel                                                            &
    !$omp firstprivate(num_nb14)                                              &
    !$omp private(id, my_id, ini_nb14, fin_nb14, i, k, j, rij, rij2, inv_rij, &
    !$omp         work, coef, lj6, lj12, term_elec, switch, dswitch, inv_rij2,&
    !$omp         inv_r6, inv_r12, atmcls_i, atmcls_j, factor_i, factor_j,    &
    !$omp         inv_lambda_i, inv_lambda_j, dij, c2_switch, c4_switch,      &
    !$omp         rvdw_i, vol_i, rtmp, elec_i, force_local, x_i, x_j,         &
    !$omp         alpha_4pi_i)                                                &
    !$omp reduction(+:eelec,evdw,esolv_tmp) 
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    my_id = my_city_rank * nthread + id

    do i = 1, natom-1

      if (mod(i-1,nproc_city*nthread) == my_id) then

        rtmp(1:3) = coord(1:3,i)
        elec_i    = charge(i) * ELECOEF * E14FAC
        atmcls_i  = atmcls(i)
        rvdw_i    = rvdw(atmcls_i)
        vol_i     = vol(2,atmcls_i)
        alpha_4pi_i  = alpha_4pi(2,atmcls_i)
        inv_lambda_i = inv_lambda(2,atmcls_i)

        force_local(1:3) = 0.0_wp

        do k = 1, num_nb14_calc(i)

          j = enefunc%nb14_calc_list(k,i)

          ! compute distance
          !
          dij(1:3) = rtmp(1:3) - coord(1:3,j)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij      = sqrt(rij2)
          inv_rij  = 1.0_wp / rij 

          ! lennard-jones energy and gradient
          !
          atmcls_j = atmcls(j)

          lj6  = enefunc%nb14_lj6 (atmcls_i,atmcls_j)
          lj12 = enefunc%nb14_lj12(atmcls_i,atmcls_j)

          if (rij2 > switchdist2) then
            c2_switch = cutoff2 - rij2
            c4_switch = c2_switch*(cutoff2 + 2.0_wp*rij2 - 3.0_wp*switchdist2)
            switch  = c2_switch*c4_switch*c1_switch
            dswitch = 4.0_wp*c1_switch*rij*(c2_switch*c2_switch - c4_switch)
          else
            switch  = 1.0_wp
            dswitch = 0.0_wp
          end if

          inv_rij2 = inv_rij * inv_rij
          inv_r6   = inv_rij2 * inv_rij2 * inv_rij2
          inv_r12  = inv_r6 * inv_r6
          evdw     = evdw + switch*(lj12 * inv_r12 - lj6 * inv_r6)
          coef     = ( dswitch * (lj12*inv_r12 - lj6*inv_r6) * rij  &
                     - switch * (12.0_wp*lj12*inv_r12-6.0_wp*lj6*inv_r6) )*inv_rij2

          ! electrostatic gradient
          !
          term_elec = elec_i * charge(j) * inv_rij2
          eelec = eelec + term_elec
          coef  = coef  + term_elec*(-2.0_wp*inv_rij2)

          ! solvation free energy in water
          !
          if (rij < EEF1_CUTOFF_DIST) then

            inv_lambda_j = inv_lambda(2,atmcls_j)

            x_i = inv_lambda_i*(rij - rvdw_i)
            x_j = inv_lambda_j*(rij - rvdw(atmcls_j))

            factor_i = exp(-(x_i*x_i))*inv_rij2*vol(2,atmcls_j)
            factor_j = exp(-(x_j*x_j))*inv_rij2*vol_i

            factor_i = factor_i*alpha_4pi_i
            factor_j = factor_j*alpha_4pi(2,atmcls_j)

            esolv_tmp = esolv_tmp + factor_i + factor_j

            coef  = coef + 2.0_wp*factor_i*inv_rij*(x_i*inv_lambda_i+inv_rij) &
                         + 2.0_wp*factor_j*inv_rij*(x_j*inv_lambda_j+inv_rij)

          end if

          ! store force
          !
          work(1:3) = coef*dij(1:3)
          force_local(1:3) = force_local(1:3) - work(1:3)
          force(1:3,j,id+1)= force(1:3,j,id+1)+ work(1:3)

        end do

        force(1:3,i,id+1) = force(1:3,i,id+1) + force_local(1:3)

      end if

    end do
    !$omp end parallel

    esolv = esolv - esolv_tmp 

    return

  end subroutine compute_energy_nonbond14_eef1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond15_eef1
  !> @brief        calculate nonbonded energy with pairlist (NOBC) with EEF1
  !! @authors      TM
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !! @note         T.Lazaridis & M.Karplus, Proteins, 35, 133−152 (1999)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond15_eef1(enefunc, molecule, pairlist, coord, &
                                           force, virial, eelec, evdw, esolv)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eelec
    real(wp),                 intent(inout) :: evdw
    real(wp),                 intent(inout) :: esolv

    ! local variables
    real(wp)                  :: dij(1:3), rij, rij2, inv_rij2
    real(wp)                  :: inv_rij, inv_r6, inv_r12, inv_r10, dr126
    real(wp)                  :: lj6, lj12
    real(wp)                  :: term_elec, shift, dshift, switch, dswitch
    real(wp)                  :: cutoff2, switchdist2
    real(wp)                  :: c1_switch, c2_switch, c4_switch, coef
    real(wp)                  :: inv_lambda_i, inv_lambda_j, x_i, x_j
    real(wp)                  :: work(1:3), force_local(1:3), esolv_tmp
    real(wp)                  :: alpha_4pi_i, factor_i, factor_j
    real(wp)                  :: rtmp(1:3), elec_i, rvdw_i, vol_i
    integer                   :: i, j, k, l, natom, id, my_id
    integer                   :: num_nb15, ini_nb15, fin_nb15
    integer                   :: atmcls_i, atmcls_j

    real(wp),         pointer :: charge(:)
    real(wp),         pointer :: dielec_const
    real(wp),         pointer :: cutoff, switchdist
    real(wp),         pointer :: alpha_4pi(:,:)
    real(wp),         pointer :: rvdw(:), vol(:,:), inv_lambda(:,:)
    integer,          pointer :: atmcls(:)
    integer,          pointer :: num_nb15_calc(:,:), nb15_calc_list(:,:)

#ifdef OMP
    integer                   :: omp_get_thread_num
#endif

    charge         => molecule%charge
    atmcls         => molecule%atom_cls_no
    dielec_const   => enefunc%dielec_const
    cutoff         => enefunc%cutoffdist
    switchdist     => enefunc%switchdist
    num_nb15_calc  => pairlist%num_nb15_calc
    nb15_calc_list => pairlist%nb15_calc_list
    natom          =  molecule%num_atoms
    cutoff2        =  cutoff * cutoff
    switchdist2    =  switchdist * switchdist

    c1_switch      =  1.0_wp/(cutoff2-switchdist2)
    c1_switch      =  c1_switch*c1_switch*c1_switch
    alpha_4pi      => enefunc%eef1%alpha_4pi
    rvdw           => enefunc%eef1%rvdw
    vol            => enefunc%eef1%volume
    inv_lambda     => enefunc%eef1%inv_lambda

    ! calculate energy and gradient
    !
    num_nb15 = 0
    esolv_tmp = 0.0_wp

    !$omp parallel                                                             &
    !$omp firstprivate(num_nb15)                                               &
    !$omp private(id, my_id, ini_nb15, fin_nb15, i, k, j, rij, rij2, inv_rij,  &
    !$omp         coef, dr126, lj6, lj12, switch, dswitch, term_elec, inv_rij2,&
    !$omp         inv_r6, inv_r12, atmcls_i, atmcls_j, factor_i, factor_j,     &
    !$omp         inv_lambda_i, inv_lambda_j, dij, c2_switch, c4_switch,       &
    !$omp         rvdw_i, vol_i, x_i, x_j, rtmp, elec_i, force_local, work,    &
    !$omp         alpha_4pi_i)                                                 &
    !$omp reduction(+:eelec,evdw,esolv_tmp) 
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    my_id = my_city_rank * nthread + id

    do i = 1, natom-1

      ini_nb15 = num_nb15 + 1
      fin_nb15 = num_nb15 + num_nb15_calc(i,id+1)
      num_nb15 = fin_nb15

      if (mod(i-1,nproc_city*nthread) /= my_id) cycle

      rtmp(1:3) = coord(1:3,i)
      elec_i    = charge(i) * ELECOEF
      atmcls_i  = atmcls(i)
      rvdw_i    = rvdw(atmcls_i)
      vol_i     = vol(2,atmcls_i)
      alpha_4pi_i  = alpha_4pi(2,atmcls_i)
      inv_lambda_i = inv_lambda(2,atmcls_i)

      force_local(1:3) = 0.0_wp

      do k = ini_nb15, fin_nb15

        j = nb15_calc_list(k,id+1)

        ! compute distance
        !
        dij(1:3) = rtmp(1:3) - coord(1:3,j)
        rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        ! cutoff
        !
        if (rij2 < cutoff2) then

          rij = sqrt(rij2)        

          if (rij2 > switchdist2) then
            c2_switch = cutoff2 - rij2
            c4_switch = c2_switch*(cutoff2 + 2.0_wp*rij2 - 3.0_wp*switchdist2)
            switch  = c2_switch*c4_switch*c1_switch
            dswitch = 4.0_wp*c1_switch*rij*(c2_switch*c2_switch - c4_switch)
          else
            switch  = 1.0_wp
            dswitch = 0.0_wp
          end if

          atmcls_j = atmcls(j)

          lj6  = enefunc%nonb_lj6 (atmcls_i,atmcls_j)
          lj12 = enefunc%nonb_lj12(atmcls_i,atmcls_j)

          ! lj energy and gradient
          !
          inv_rij  = 1.0_wp / rij
          inv_rij2 = inv_rij * inv_rij
          inv_r6   = inv_rij2 * inv_rij2 * inv_rij2
          inv_r12  = inv_r6 * inv_r6
          inv_r6   = lj6 * inv_r6
          inv_r12  = lj12 * inv_r12
          dr126    = inv_r12 - inv_r6
          evdw     = evdw + switch*dr126 
            
          ! gradient
          !
          coef = (  dswitch * dr126 * rij &
                   - switch * (12.0_wp*inv_r12-6.0_wp*inv_r6) ) * inv_rij2

          ! electrostatic energy with distance-dependent dielectric constant
          ! and gradient 
          !
          term_elec = elec_i * charge(j) * inv_rij2
          eelec = eelec + term_elec
          coef  = coef  + term_elec*(-2.0_wp*inv_rij2)

          ! solvation free energy in water and gradient
          !
          if (rij < EEF1_CUTOFF_DIST) then

            inv_lambda_j = inv_lambda(2,atmcls_j)

            x_i = inv_lambda_i*(rij - rvdw_i)
            x_j = inv_lambda_j*(rij - rvdw(atmcls_j))

            factor_i = exp(-(x_i*x_i))*inv_rij2*vol(2,atmcls_j)
            factor_j = exp(-(x_j*x_j))*inv_rij2*vol_i

            factor_i = factor_i*alpha_4pi_i
            factor_j = factor_j*alpha_4pi(2,atmcls_j)

            esolv_tmp = esolv_tmp + factor_i + factor_j

            coef  = coef + 2.0_wp*factor_i*inv_rij*(x_i*inv_lambda_i+inv_rij) &
                         + 2.0_wp*factor_j*inv_rij*(x_j*inv_lambda_j+inv_rij)

          end if

          ! store force
          !
          work(1:3) = coef*dij(1:3)
          force_local(1:3) = force_local(1:3) - work(1:3)
          force(1:3,j,id+1)= force(1:3,j,id+1)+ work(1:3)

        end if

      end do

      force(1:3,i,id+1) = force(1:3,i,id+1) + force_local(1:3)

    end do
    !$omp end parallel

    esolv = esolv - esolv_tmp

    return

  end subroutine compute_energy_nonbond15_eef1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_imm1
  !> @brief        calculate nonbonded14 energy with IMM1
  !  @authors      TM
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !! @note         T.Lazaridis & M.Karplus, Proteins, 35, 133−152 (1999)
  !!               T.Lazaridis, Proteins, 52, 176-192 (2003)
  !!               A. Rahaman & T. Lazaridis, BBA, 1838, 1440-1447 (2014)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_imm1(enefunc, molecule, coord, force, &
                                           virial, eelec, evdw, esolv)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_molecule), target, intent(in)    :: molecule
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eelec
    real(wp),                 intent(inout) :: evdw
    real(wp),                 intent(inout) :: esolv

    ! local variables
    real(wp)                  :: dij(1:3), rij, rij2, inv_rij2
    real(wp)                  :: inv_rij, inv_r3, inv_r6, inv_r10, inv_r12
    real(wp)                  :: lj6, lj12, elec_i, E14FAC
    real(wp)                  :: term_elec, shift, dshift, switch, dswitch
    real(wp)                  :: cutoff2, switchdist2
    real(wp)                  :: c1_switch, c2_switch, c4_switch, coef
    real(wp)                  :: c12, c6, term_lj12, term_lj6
    real(wp)                  :: factor_i, factor_j, alpha_4pi_i
    real(wp)                  :: inv_lambda_i, inv_lambda_j, x_i, x_j
    real(wp)                  :: work(1:3), force_local(1:3), rtmp(1:3)
    real(wp)                  :: esolv_tmp, esolv_i, esolv_j
    real(wp)                  :: fij, z, ihalf_thick, a
    real(wp)                  :: zt_i, zt_j, FF_i, FF_j, ddg_i, ddg_j
    real(wp)                  :: gfree_i, gfree_j, rvdw_i, vol_i
    real(wp)                  :: dqfac, sqfifj
    integer                   :: i, j, k, l, n, natom
    integer                   :: my_id, id
    integer                   :: num_nb14, ini_nb14, fin_nb14
    integer                   :: atmcls_i, atmcls_j

    real(wp),         pointer :: charge(:)
    real(wp),         pointer :: dielec_const
    real(wp),         pointer :: cutoff, switchdist
    real(wp),         pointer :: alpha_4pi(:,:), gfree_t(:,:)
    real(wp),         pointer :: rvdw(:), vol(:,:), inv_lambda(:,:)
    integer,          pointer :: atmcls(:)
    integer,          pointer :: num_nb14_calc(:), nb14_calc_list(:,:)

#ifdef OMP
    integer                   :: omp_get_thread_num
#endif

    charge         => molecule%charge
    atmcls         => molecule%atom_cls_no
    dielec_const   => enefunc%dielec_const
    cutoff         => enefunc%cutoffdist
    switchdist     => enefunc%switchdist
    num_nb14_calc  => enefunc%num_nb14_calc
    nb14_calc_list => enefunc%nb14_calc_list
    natom          =  molecule%num_atoms
    cutoff2        =  cutoff * cutoff
    switchdist2    =  switchdist * switchdist
    alpha_4pi      => enefunc%eef1%alpha_4pi
    gfree_t        => enefunc%eef1%gfree_t
    rvdw           => enefunc%eef1%rvdw
    vol            => enefunc%eef1%volume
    inv_lambda     => enefunc%eef1%inv_lambda

    c1_switch      =  1.0_wp/(cutoff2-switchdist2)
    c1_switch      =  c1_switch*c1_switch*c1_switch
    a              =  enefunc%eef1%imm1_factor_a
    n              =  enefunc%eef1%imm1_exponent_n
    ihalf_thick    =  1.0_wp/(enefunc%eef1%imm1_memb_thick*0.5_wp)

    if (enefunc%forcefield == ForcefieldCHARMM19) then
      E14FAC   = 0.4_wp
    else
      E14FAC   = 1.0_wp
    end if

    ! calculate energy and gradient
    !
    num_nb14  = 0
    esolv_tmp = 0.0_wp

    !$omp parallel                                                            &
    !$omp firstprivate(num_nb14)                                              &
    !$omp private(id, my_id, ini_nb14, fin_nb14, i, k, j, rij, rij2, inv_rij, &
    !$omp         work, coef, lj6, lj12, term_elec, switch, dswitch, inv_rij2,&
    !$omp         inv_r6, inv_r12, atmcls_i, atmcls_j, factor_i, factor_j,    &
    !$omp         inv_lambda_i, inv_lambda_j, dij, c2_switch, c4_switch,      &
    !$omp         rtmp, elec_i, force_local, x_i, x_j, gfree_i, gfree_j, fij, &
    !$omp         alpha_4pi_i, rvdw_i, vol_i, esolv_i, esolv_j, FF_i, FF_j,   &
    !$omp         ddg_i, ddg_j, sqfifj, dqfac)                                &
    !$omp reduction(+:eelec,evdw,esolv_tmp)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    my_id = my_city_rank * nthread + id
  
    do i = 1, natom-1
  
      if (mod(i-1,nproc_city*nthread) == my_id) then

        rtmp(1:3) = coord(1:3,i)
        elec_i    = charge(i) * ELECOEF * E14FAC
        atmcls_i  = atmcls(i)
        rvdw_i    = rvdw(atmcls_i)
        vol_i     = vol(2,atmcls_i)
        FF_i      = FF(i)
        gfree_i   = FF_i*gfree_t(2,atmcls_i) + (1.0_wp-FF_i)*gfree_t(1,atmcls_i)
        ddg_i     = gfree_t(2,atmcls_i) - gfree_t(1,atmcls_i)
        alpha_4pi_i  = alpha_4pi(2,atmcls_i)
        inv_lambda_i = inv_lambda(2,atmcls_i)

        force_local(1:3) = 0.0_wp

        do k = 1, num_nb14_calc(i)

          j = enefunc%nb14_calc_list(k,i)

          ! compute distance
          !
          dij(1:3) = rtmp(1:3) - coord(1:3,j)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij      = sqrt(rij2)
          inv_rij  = 1.0_wp / rij 

          ! lennard-jones energy and gradient
          !
          atmcls_j = atmcls(j)

          lj6  = enefunc%nb14_lj6 (atmcls_i,atmcls_j)
          lj12 = enefunc%nb14_lj12(atmcls_i,atmcls_j)

          if (rij2 > switchdist2) then
            c2_switch = cutoff2 - rij2
            c4_switch = c2_switch*(cutoff2 + 2.0_wp*rij2 - 3.0_wp*switchdist2)
            switch  = c2_switch*c4_switch*c1_switch
            dswitch = 4.0_wp*c1_switch*rij*(c2_switch*c2_switch - c4_switch)
          else
            switch  = 1.0_wp
            dswitch = 0.0_wp
          end if

          inv_rij2 = inv_rij * inv_rij
          inv_r6   = inv_rij2 * inv_rij2 * inv_rij2
          inv_r12  = inv_r6 * inv_r6
          evdw     = evdw + switch*(lj12 * inv_r12 - lj6 * inv_r6)
          coef     = ( dswitch * (lj12*inv_r12 - lj6*inv_r6) * rij  &
                      - switch * (12.0_wp*lj12*inv_r12-6.0_wp*lj6*inv_r6) )*inv_rij2

          ! solvation free energy in membrane/water
          !   Note that alpha_4pi does not include dGfree
          !
          FF_j = FF(j)

          if (rij < EEF1_CUTOFF_DIST) then

            inv_lambda_j = inv_lambda(2,atmcls_j)

            x_i = inv_lambda_i*(rij - rvdw_i)
            x_j = inv_lambda_j*(rij - rvdw(atmcls_j))

            gfree_j  = FF_j*gfree_t(2,atmcls_j) + (1.0_wp-FF_j)*gfree_t(1,atmcls_j)
            ddg_j    = gfree_t(2,atmcls_j) - gfree_t(1,atmcls_j)

            factor_i = exp(-(x_i*x_i))*inv_rij2*vol(2,atmcls_j)
            factor_j = exp(-(x_j*x_j))*inv_rij2*vol_i

            factor_i = factor_i*alpha_4pi_i
            factor_j = factor_j*alpha_4pi(2,atmcls_j)

            esolv_i  = factor_i*gfree_i
            esolv_j  = factor_j*gfree_j

            esolv_tmp = esolv_tmp + esolv_i + esolv_j

            coef  = coef +  2.0_wp*esolv_i*inv_rij*(x_i*inv_lambda_i+inv_rij) &
                         +  2.0_wp*esolv_j*inv_rij*(x_j*inv_lambda_j+inv_rij)

            force_local(1:3) = force_local(1:3) + dmfac(1:3,i)*factor_i*ddg_i
            force(1:3,j,id+1)= force(1:3,j,id+1)+ dmfac(1:3,j)*factor_j*ddg_j

          end if

          ! electrostatic energy and gradient
          !   Ref) A. Rahaman & T. Lazaridis, BBA, 1838, 1440–1447, 2014
          !
          if (enefunc%forcefield == ForcefieldCHARMM) then

            term_elec = elec_i * charge(j) * inv_rij2
            eelec = eelec + term_elec
            coef  = coef  + term_elec*(-2.0_wp*inv_rij2)

            work(1:3)        = coef*dij(1:3)
            force_local(1:3) = force_local(1:3) - work(1:3)
            force(1:3,j,id+1)= force(1:3,j,id+1)+ work(1:3)

          else if (enefunc%forcefield == ForcefieldCHARMM19) then

            sqfifj    = sqrt(FF_i*FF_j)
            fij       = a + (1.0_wp - a)*sqfifj
            term_elec = elec_i * charge(j) * inv_rij * inv_rij**(fij)
            eelec     = eelec + term_elec

            dqfac     = - term_elec*(1.0_wp - a)*log(rij2)/(4.0_wp*sqfifj)

            force_local(1:3) = force_local(1:3) - dqfac*FF_j*dmfac(1:3,i)
            force(1:3,j,id+1)= force(1:3,j,id+1)- dqfac*FF_i*dmfac(1:3,j)

          end if

        end do

        force(1:3,i,id+1) = force(1:3,i,id+1) + force_local(1:3)

      end if

    end do
    !$omp end parallel

    esolv = esolv - esolv_tmp 

    return

  end subroutine compute_energy_nonbond14_imm1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond15_imm1
  !> @brief        calculate nonbonded energy with pairlist (NOBC) with IMM1
  !! @authors      TM
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !! @note         T.Lazaridis & M.Karplus, Proteins, 35, 133−152 (1999)
  !!               T.Lazaridis, Proteins, 52, 176-192 (2003)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond15_imm1(enefunc, molecule, pairlist, coord, &
                                           force, virial, eelec, evdw, esolv)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eelec
    real(wp),                 intent(inout) :: evdw
    real(wp),                 intent(inout) :: esolv

    ! local variables
    real(wp)                  :: dij(1:3), rij, rij2
    real(wp)                  :: inv_rij, inv_r6, inv_r12, inv_r10, dr126
    real(wp)                  :: lj6, lj12, elec_i
    real(wp)                  :: term_elec, shift, dshift, switch, dswitch
    real(wp)                  :: cutoff2, switchdist2
    real(wp)                  :: c1_switch, c2_switch, c4_switch, coef
    real(wp)                  :: inv_lambda_i, inv_lambda_j, x_i, x_j
    real(wp)                  :: rtmp(1:3), force_local(1:3)
    real(wp)                  :: work(1:3), esolv_tmp, dqfac, sqfifj
    real(wp)                  :: factor_i, factor_j, alpha_4pi_i
    real(wp)                  :: fij, z, ihalf_thick, a
    real(wp)                  :: rvdw_i, vol_i, inv_rij2
    real(wp)                  :: gfree_i, gfree_j, zt_i, zt_j, ci, cj
    real(wp)                  :: esolv_i, esolv_j, FF_i, FF_j, ddg_i, ddg_j
    integer                   :: i, j, k, l, n, natom
    integer                   :: num_nb15, ini_nb15, fin_nb15
    integer                   :: id, my_id
    integer                   :: atmcls_i, atmcls_j

    real(wp),         pointer :: charge(:)
    real(wp),         pointer :: dielec_const
    real(wp),         pointer :: cutoff, switchdist
    real(wp),         pointer :: alpha_4pi(:,:), gfree_t(:,:)
    real(wp),         pointer :: rvdw(:), vol(:,:), inv_lambda(:,:)
    integer,          pointer :: atmcls(:)
    integer,          pointer :: num_nb15_calc(:,:), nb15_calc_list(:,:)

#ifdef OMP
    integer                   :: omp_get_thread_num
#endif

    charge         => molecule%charge
    atmcls         => molecule%atom_cls_no
    dielec_const   => enefunc%dielec_const
    cutoff         => enefunc%cutoffdist
    switchdist     => enefunc%switchdist
    num_nb15_calc  => pairlist%num_nb15_calc
    nb15_calc_list => pairlist%nb15_calc_list
    natom          =  molecule%num_atoms
    cutoff2        =  cutoff * cutoff
    switchdist2    =  switchdist * switchdist

    c1_switch      =  1.0_wp/(cutoff2-switchdist2)
    c1_switch      =  c1_switch*c1_switch*c1_switch
    alpha_4pi      => enefunc%eef1%alpha_4pi
    gfree_t        => enefunc%eef1%gfree_t
    rvdw           => enefunc%eef1%rvdw
    vol            => enefunc%eef1%volume
    inv_lambda     => enefunc%eef1%inv_lambda
    a              =  enefunc%eef1%imm1_factor_a
    n              =  enefunc%eef1%imm1_exponent_n
    ihalf_thick    =  1.0_wp/(enefunc%eef1%imm1_memb_thick*0.5_wp)

    ! calculate energy and gradient
    !
    num_nb15 = 0
    esolv_tmp = 0.0_wp

    !$omp parallel                                                             &
    !$omp firstprivate(num_nb15)                                               &
    !$omp private(id, my_id, ini_nb15, fin_nb15, i, k, j, rij, rij2, inv_rij,  &
    !$omp         coef, dr126, lj6, lj12, switch, dswitch, term_elec, dqfac,   &
    !$omp         inv_r6, inv_r12, atmcls_i, atmcls_j, factor_i, factor_j, dij,&
    !$omp         inv_lambda_i, inv_lambda_j, c2_switch, c4_switch, sqfifj,    &
    !$omp         x_i, x_j, rtmp, elec_i, force_local, work, gfree_i, gfree_j, &
    !$omp         alpha_4pi_i, rvdw_i, vol_i, fij, z, inv_rij2, FF_i, FF_j,    &
    !$omp         esolv_i, esolv_j, ddg_i, ddg_j, ci, cj)                      &
    !$omp reduction(+:eelec,evdw,esolv_tmp)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    my_id = my_city_rank * nthread + id
  
    do i = 1, natom-1

      ini_nb15 = num_nb15 + 1
      fin_nb15 = num_nb15 + num_nb15_calc(i,id+1)
      num_nb15 = fin_nb15
  
      if (mod(i-1,nproc_city*nthread) /= my_id) cycle

      rtmp(1:3) = coord(1:3,i)
      elec_i    = charge(i) * ELECOEF
      atmcls_i  = atmcls(i)
      rvdw_i    = rvdw(atmcls_i)
      vol_i     = vol(2,atmcls_i)
      FF_i      = FF(i)
      gfree_i   = FF_i*gfree_t(2,atmcls_i) + (1.0_wp-FF_i)*gfree_t(1,atmcls_i)
      ddg_i     = gfree_t(2,atmcls_i) - gfree_t(1,atmcls_i)
      alpha_4pi_i  = alpha_4pi(2,atmcls_i)
      inv_lambda_i = inv_lambda(2,atmcls_i)

      force_local(1:3) = 0.0_wp

      do k = ini_nb15, fin_nb15
  
        j = nb15_calc_list(k,id+1)
  
        ! compute distance
        !
        dij(1:3) = rtmp(1:3) - coord(1:3,j)
        rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        ! cutoff
        !
        if (rij2 < cutoff2) then

          rij = sqrt(rij2)        

          ! lj energy and gradient
          !
          if (rij2 > switchdist2) then
            c2_switch = cutoff2 - rij2
            c4_switch = c2_switch*(cutoff2 + 2.0_wp*rij2 - 3.0_wp*switchdist2)
            switch  = c2_switch*c4_switch*c1_switch
            dswitch = 4.0_wp*c1_switch*rij*(c2_switch*c2_switch - c4_switch)
          else
            switch  = 1.0_wp
            dswitch = 0.0_wp
          end if

          atmcls_j = atmcls(j)

          lj6  = enefunc%nonb_lj6 (atmcls_i,atmcls_j)
          lj12 = enefunc%nonb_lj12(atmcls_i,atmcls_j)

          inv_rij  = 1.0_wp / rij
          inv_rij2 = inv_rij * inv_rij
          inv_r6   = inv_rij2 * inv_rij2 * inv_rij2
          inv_r12  = inv_r6 * inv_r6
          inv_r6   = lj6 * inv_r6
          inv_r12  = lj12 * inv_r12
          dr126    = inv_r12 - inv_r6
          evdw     = evdw + switch*dr126 
            
          coef = (  dswitch * dr126 * rij &
                   - switch * (12.0_wp*inv_r12-6.0_wp*inv_r6) ) * inv_rij2

          ! solvation free energy in membrane/water
          !   Note that alpha_4pi does not include dGfree
          !
          FF_j = FF(j)

          if (rij < EEF1_CUTOFF_DIST) then

            inv_lambda_j = inv_lambda(2,atmcls_j)

            x_i = inv_lambda_i*(rij - rvdw_i)
            x_j = inv_lambda_j*(rij - rvdw(atmcls_j))

            gfree_j  = FF_j*gfree_t(2,atmcls_j) + (1.0_wp-FF_j)*gfree_t(1,atmcls_j)
            ddg_j    = gfree_t(2,atmcls_j) - gfree_t(1,atmcls_j)

            factor_i = exp(-(x_i*x_i))*inv_rij2*vol(2,atmcls_j)
            factor_j = exp(-(x_j*x_j))*inv_rij2*vol_i

            factor_i = factor_i*alpha_4pi_i
            factor_j = factor_j*alpha_4pi(2,atmcls_j)

            esolv_i  = factor_i*gfree_i
            esolv_j  = factor_j*gfree_j

            esolv_tmp = esolv_tmp + esolv_i + esolv_j
 
            coef  = coef + 2.0_wp*esolv_i*inv_rij*(x_i*inv_lambda_i+inv_rij) &
                         + 2.0_wp*esolv_j*inv_rij*(x_j*inv_lambda_j+inv_rij)

            force_local(1:3) = force_local(1:3) + dmfac(1:3,i)*factor_i*ddg_i
            force(1:3,j,id+1)= force(1:3,j,id+1)+ dmfac(1:3,j)*factor_j*ddg_j

          end if

          ! electrostatic energy and force
          !
!          ci        = a + (1.0_wp - a)*FF_i
!          cj        = a + (1.0_wp - a)*FF_j
!          fij       = sqrt(ci*cj)
!          term_elec = elec_i * charge(j) * inv_rij * inv_rij**(fij)
!          eelec     = eelec + term_elec
!
!          dqfac     = - term_elec*(1.0_wp - a)*log(rij2)/(4.0_wp*fij)
!
!          force_local(1:3) = force_local(1:3) - dqfac*cj*dmfac(1:3,i)
!          force(1:3,j)     = force(1:3,j)     - dqfac*ci*dmfac(1:3,j)

          sqfifj    = sqrt(FF_i*FF_j)
          fij       = a + (1.0_wp - a)*sqfifj
          term_elec = elec_i * charge(j) * inv_rij * inv_rij**(fij)
          eelec     = eelec + term_elec

          dqfac     = - term_elec*(1.0_wp - a)*log(rij2)/(4.0_wp*sqfifj)

          force_local(1:3) = force_local(1:3) - dqfac*FF_j*dmfac(1:3,i)
          force(1:3,j,id+1)= force(1:3,j,id+1)- dqfac*FF_i*dmfac(1:3,j)

          ! store force
          !
          coef = coef + (-1.0_wp - fij)*term_elec*inv_rij2

          work(1:3) = coef*dij(1:3)
          force_local(1:3) = force_local(1:3) - work(1:3)
          force(1:3,j,id+1)= force(1:3,j,id+1)+ work(1:3)

        end if

      end do

      force(1:3,i,id+1) = force(1:3,i,id+1) + force_local(1:3)

    end do
    !$omp end parallel

    esolv = esolv - esolv_tmp

    return

  end subroutine compute_energy_nonbond15_imm1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_ellipsoid_depth
  !> @brief        calculate minimum distance between a point and ellipsoid
  !! @authors      TM
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_ellipsoid_depth(a, b, c, m1, m2, natom, coord)

    ! formal arguments
    real(wp), intent(in)  :: a, b, c, m1, m2
    integer,  intent(in)  :: natom
    real(wp), intent(in)  :: coord(:,:)

    ! local variables
    integer               :: isd, ih, ig, il, is, i
    integer               :: ncycle, icycle, nlen, ixx
    real(wp)              :: x0, y0, z0, x, y, z, r
    real(wp)              :: x1, y1, z1, x2, y2, z2, d1, d2, ini_f
    real(wp)              :: theta, phi, dtheta, dphi
    real(wp)              :: tmp_theta(0:2), tmp_phi(0:2), f(0:2)
    real(wp)              :: com_theta, com_phi, hlamd1, hlamd2
    real(wp)              :: exp_theta, exp_phi, exp_f
    real(wp)              :: red_theta, red_phi, red_f
    real(wp)              :: ref_theta, ref_phi, ref_f
    real(wp)              :: s1, s2, t
    logical               :: middle(0:2)
    integer,  parameter   :: nsteps    = 10000
    real(wp), parameter   :: tolerance = 0.00001_wp
    real(wp), parameter   :: lamd1     = 0.33_wp
    real(wp), parameter   :: lamd2     = 1.70_wp


    if (.not. allocated(coord_min)) then
      allocate(coord_min(3,natom), coord_tmp(3,natom))
      call get_loop_index(natom, istart, iend)
    end if

    coord_min(:,:) = 0.0_wp

    hlamd1 = 0.5_wp*(1.0_wp - lamd1)
    hlamd2 = 0.5_wp*(1.0_wp - lamd2)
    s1     = 2.0_wp/m1
    s2     = 2.0_wp/m2
    t      = m2/m1

    !$omp parallel do default (none)                                         &
    !$omp private(i, x, y, z, theta, dtheta, phi, dphi, tmp_theta, tmp_phi,  &
    !$omp         f, isd, ih, ig, il, middle, com_theta, com_phi, ref_theta, &
    !$omp         ref_phi, ref_f, red_theta, red_phi, red_f, exp_theta,      &
    !$omp         exp_phi, exp_f, is, r, x0, y0, z0, d1, d2, x1, y1, z1,     &
    !$omp         x2, y2, z2, ini_f)                                         &
    !$omp shared (natom, a, b, c, m1, m2, coord, coord_min, istart, iend,    &
    !$omp         hlamd1, hlamd2, s1, s2, t)
    !
    do i = istart, iend

      x = coord(1,i)
      y = coord(2,i)
      z = coord(3,i)

      ini_f = ((abs(x)/a)**s2 + (abs(y)/b)**s2)**t + (abs(z)/c)**s1

      if (ini_f > 0.7_wp) then
        r = dsqrt(x*x + y*y + z*z)
        theta = dacos(z/r)
        if (abs(dsin(theta)) < 0.00001_wp .or. abs(y) < 0.00001_wp) then
          phi = (90.0_wp - sign(1.0_wp,x)*90.0_wp)*RAD
        else
          phi = sign(1.0_wp,y)*dacos((x/r)/dsin(theta))
        end if
        dtheta = sign(1.0_wp,z)*5.0_wp*RAD
        dphi   = 5.0_wp*RAD
      else
        phi    = (sign(1.0_wp,y)*90.0_wp - sign(1.0_wp,y)*sign(1.0_wp,x)*45.0_wp)*RAD
        theta  = (90.0_wp - sign(1.0_wp,z)*60.0_wp)*RAD
        dphi   = 15.0_wp*RAD
        dtheta = sign(1.0_wp,z)*30.0_wp*RAD
      end if

      tmp_theta(0) = theta
      tmp_phi  (0) = phi
      tmp_theta(1) = tmp_theta(0) + dtheta
      tmp_phi  (1) = tmp_phi  (0) - dphi
      tmp_theta(2) = tmp_theta(0) + dtheta
      tmp_phi  (2) = tmp_phi  (0) + dphi

      f(0) = g(tmp_theta(0), tmp_phi(0), a, b, c, m1, m2, x, y, z, x0, y0, z0)
      f(1) = g(tmp_theta(1), tmp_phi(1), a, b, c, m1, m2, x, y, z, x0, y0, z0)
      f(2) = g(tmp_theta(2), tmp_phi(2), a, b, c, m1, m2, x, y, z, x0, y0, z0)

      do isd = 1, nsteps

        ih = maxloc(f,dim=1) - 1
        il = minloc(f,dim=1) - 1

        middle(0:2) = .true.
        middle(ih)  = .false.
        middle(il)  = .false.
        do is = 0, 2
          if (middle(is)) ig = is
        end do

        if (2.0_wp*abs(f(ih)-f(il))/(abs(f(ih)) + abs(f(il))) < tolerance) then
          exit
        else if (isd == nsteps) then
          if (f(0) < 1.0e-15 .and. f(1) < 1.0e-15 .and. f(2) < 1.0e-15) then
            exit
          else
            write(MsgOut,*) x, y, z
            call error_msg('Compute_Ellipsoid_Depth> Depth calculation failed to converge')
          end if
        end if

        com_theta = tmp_theta(ig) + tmp_theta(il)
        com_phi   = tmp_phi  (ig) + tmp_phi  (il)

        ref_theta = com_theta - tmp_theta(ih)
        ref_phi   = com_phi   - tmp_phi  (ih)
        ref_f     = g(ref_theta, ref_phi, a, b, c, m1, m2, x, y, z, x0, y0, z0)

        if (ref_f < f(il))  then
          exp_theta = lamd2*ref_theta + hlamd2*com_theta
          exp_phi   = lamd2*ref_phi   + hlamd2*com_phi
          exp_f     = g(exp_theta, exp_phi, a, b, c, m1, m2, x, y, z, x0, y0, z0)
          if (exp_f < ref_f) then
            tmp_theta(ih) = exp_theta
            tmp_phi  (ih) = exp_phi
            f        (ih) = exp_f
          else
            tmp_theta(ih) = ref_theta
            tmp_phi  (ih) = ref_phi
            f        (ih) = ref_f
          end if
        else if (f(ig) < ref_f) then
          if (ref_f < f(ih)) then
            tmp_theta(ih) = ref_theta
            tmp_phi  (ih) = ref_phi
            f        (ih) = ref_f
          end if
          red_theta = lamd1*tmp_theta(ih) + hlamd1*com_theta
          red_phi   = lamd1*tmp_phi  (ih) + hlamd1*com_phi
          red_f     = g(red_theta, red_phi, a, b, c, m1, m2, x, y, z, x0, y0, z0)
          if (red_f < f(ih)) then
            tmp_theta(ih) = red_theta
            tmp_phi  (ih) = red_phi
            f        (ih) = red_f
          else
            tmp_theta(ih) = 0.5_wp*(tmp_theta(ih) + tmp_theta(il))
            tmp_phi  (ih) = 0.5_wp*(tmp_phi  (ih) + tmp_phi  (il))
            tmp_theta(ig) = 0.5_wp*(tmp_theta(ig) + tmp_theta(il))
            tmp_phi  (ig) = 0.5_wp*(tmp_phi  (ig) + tmp_phi  (il))
            f(ih) = g(tmp_theta(ih), tmp_phi(ih), a, b, c, m1, m2, x, y, z, x0, y0, z0)
            f(ig) = g(tmp_theta(ig), tmp_phi(ig), a, b, c, m1, m2, x, y, z, x0, y0, z0)
          end if
        else
          tmp_theta(ih) = ref_theta
          tmp_phi  (ih) = ref_phi
          f        (ih) = ref_f
        end if

      end do

      theta = (tmp_theta(ih) + tmp_theta(ig) + tmp_theta(il))/3.0_wp
      phi   = (tmp_phi  (ih) + tmp_phi  (ig) + tmp_phi  (il))/3.0_wp
      d1    = g(theta,         phi,         a, b, c, m1, m2, x, y, z, x1, y1, z1)
      d2    = g(tmp_theta(il), tmp_phi(il), a, b, c, m1, m2, x, y, z, x2, y2, z2)

      if (d1 < d2) then
        coord_min(1,i) = x1
        coord_min(2,i) = y1
        coord_min(3,i) = z1
      else
        coord_min(1,i) = x2
        coord_min(2,i) = y2
        coord_min(3,i) = z2
      end if

    end do
    !$omp end parallel do

#ifdef HAVE_MPI_GENESIS

    coord_tmp(1:3,1:natom) = 0.0_wp
    coord_tmp(1:3,1:natom) = coord_min(1:3,1:natom)

    ncycle = (natom - 1) / mpi_drain + 1
    nlen   = mpi_drain
    ixx    = 1

    do icycle = 1, ncycle
      if (icycle == ncycle) nlen = natom - (ncycle-1) * mpi_drain
      call mpi_allreduce(coord_tmp(1,ixx), coord_min(1,ixx), 3*nlen,    &
                         mpi_wp_real, mpi_sum, mpi_comm_country, ierror)
      ixx = ixx + nlen
    end do

#endif

    return

  end subroutine compute_ellipsoid_depth

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      g
  !> @brief        calculate distance between a point and ellipsoid
  !! @authors      TM
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function g(theta, phi, a, b, c, m1, m2, x, y, z, x0, y0, z0)

    ! formal arguments
    real(wp), intent(in)  :: theta, phi, a, b, c, m1, m2, x, y, z
    real(wp), intent(out) :: x0, y0, z0

    ! local variables
    real(wp)              :: g, f, r, sc, ss
    real(wp)              :: sin_theta, cos_theta, sin_phi, cos_phi


    sin_theta = sign(1.0_wp,sin(theta))*abs(sin(theta))**m1
    cos_phi   = sign(1.0_wp,cos(phi  ))*abs(cos(phi  ))**m2
    cos_theta = sign(1.0_wp,cos(theta))*abs(cos(theta))**m1
    sin_phi   = sign(1.0_wp,sin(phi  ))*abs(sin(phi  ))**m2

    x0 = a*sin_theta*cos_phi
    y0 = b*sin_theta*sin_phi
    z0 = c*cos_theta

    g = (x - x0)**2 + (y - y0)**2 + (z - z0)**2

    return

  end function g

end module at_energy_eef1_mod
