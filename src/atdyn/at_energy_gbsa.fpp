!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_energy_gbsa_mod
!> @brief   calculate nonbonded energy with GBSA
!! @authors Takaharu Mori (TM)
! 
!  (c) Copyright 2016 RIKEN. All rights reserved.
! 
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_energy_gbsa_mod

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

  ! Is this a risk of being wrong when multiple calcs are done?
  real(wp), allocatable, save :: Hij(:)
  real(wp), allocatable, save :: sgrad(:)
  real(wp), allocatable, save :: dfact(:)
  real(wp), allocatable, save :: tmp_allreduce(:)
  real(wp),              save :: gbcutoff2
  integer,  allocatable, save :: cell_pointer(:)
  integer,  allocatable, save :: head_atom(:)
  integer,  allocatable, save :: cell_list(:,:)
  integer,  allocatable, save :: cell_linked_list(:)
  integer,               save :: ncell_xyz_old

  real(wp), parameter         :: pa  =  1.0_wp/ 3.0_wp
  real(wp), parameter         :: pb  =  2.0_wp/ 5.0_wp
  real(wp), parameter         :: pc  =  3.0_wp/ 7.0_wp
  real(wp), parameter         :: pd  =  4.0_wp/ 9.0_wp
  real(wp), parameter         :: pe  =  5.0_wp/11.0_wp
  real(wp), parameter         :: paa =  4.0_wp/ 3.0_wp
  real(wp), parameter         :: pbb = 12.0_wp/ 5.0_wp
  real(wp), parameter         :: pcc = 24.0_wp/ 7.0_wp
  real(wp), parameter         :: pdd = 40.0_wp/ 9.0_wp
  real(wp), parameter         :: pee = 60.0_wp/11.0_wp

  ! subroutines
  public   :: compute_energy_gbsa
  public   :: compute_born_radius
  private  :: compute_Hij
  private  :: compute_energy_gb
  private  :: compute_energy_sa

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_gbsa
  !> @brief        Calculate solvation free energy with the GBSA model
  !! @authors      TM
  !! @param[inout] molecule : molecule information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] enefunc  : potential energy functions information
  !! @param[inout] force    : forces of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !! @param[inout] esolv    : solvation free energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_gbsa(enefunc, molecule, pairlist, coord, force, esolv)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_pairlist),        intent(in)    :: pairlist
    real(wp),                intent(in)    :: coord(:,:)
    type(s_enefunc),         intent(inout) :: enefunc
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: esolv

    ! local variables
    integer                  :: i

    call timer(TimerSolvation, TimerOn)


    ! Putting in QM charges
    if (enefunc%qmmm%do_qmmm) then
      do i = 1, enefunc%qmmm%qm_natoms
        molecule%charge(enefunc%qmmm%qmatom_id(i)) = enefunc%qmmm%qm_charge(i)
      end do
    end if

    ! calc GB term
    !
    call timer(TimerGB, TimerOn)

    call compute_energy_gb(enefunc, molecule, pairlist, coord, force, esolv)

    call timer(TimerGB, TimerOff)

    ! calc SA term
    !
    call timer(TimerSA, TimerOn)

    call compute_energy_sa(enefunc, molecule, coord, force, esolv)

    call timer(TimerSA, TimerOff)

    ! Unsetting QM charges
    if (enefunc%qmmm%do_qmmm) then
      do i = 1, enefunc%qmmm%qm_natoms
        molecule%charge(enefunc%qmmm%qmatom_id(i)) = 0.0D+00
      end do
    end if

    call timer(TimerSolvation, TimerOff)

    return

  end subroutine compute_energy_gbsa

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_Hij
  !> @brief        calculate Hij for the GB model
  !! @authors      
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_Hij(rij, cutoff, rho_i0, rho_js, Hij)

    ! formal arguments
    real(wp),                 intent(in)    :: rij
    real(wp),                 intent(in)    :: cutoff
    real(wp),                 intent(in)    :: rho_i0
    real(wp),                 intent(in)    :: rho_js
    real(wp),                 intent(inout) :: Hij

    ! local variables
    real(wp)                  :: rij2, rho_i02, rho_js2, cutoff2
    real(wp)                  :: term1, term2, term3, tmp

    rij2    = rij * rij
    rho_i02 = rho_i0 * rho_i0
    rho_js2 = rho_js * rho_js
    cutoff2 = cutoff * cutoff

    ! for Hij
    if (rij >= cutoff + rho_js) then
      Hij  = Hij + 0.0_wp
    else if (rij > cutoff - rho_js) then
      term1  = 1.0_wp + 2.0_wp*rij/(rij - rho_js)
      term2  = (rij2 - 4.0_wp*cutoff*rij - rho_js2)/cutoff2
      term3  = 2.0_wp*log((rij - rho_js)/cutoff)
      Hij  = Hij + 0.125_wp*(term1 + term2 + term3)/rij
    else if (rij > 4.0_wp*rho_js) then
      tmp    = rho_js2/rij2
      term1  = pa + tmp*(pb + tmp*(pc + tmp*(pd + tmp*pe)))
      term2  = tmp*tmp/rho_js
      Hij  = Hij + term1*term2
    else if (rij > rho_i0 + rho_js) then
      term1  = rho_js/(rij2 - rho_js2)
      term2  = 0.5_wp*log((rij - rho_js)/(rij + rho_js))/rij
      Hij  = Hij + 0.5_wp*(term1 + term2)
    else if (rij > abs(rho_i0 - rho_js)) then
      term1  = 2.0_wp - 1.0_wp/(2.0_wp*rij*rho_i0)*(rij2+rho_i02-rho_js2)
      term1  = term1/rho_i0
      term2  = - 1.0_wp/(rij + rho_js)
      term3  = 1.0_wp/rij*log(rho_i0/(rij+rho_js))
      Hij  = Hij + 0.25_wp*(term1 + term2 + term3)
    else if (rho_i0 < rho_js) then
      term1  = rho_js/(rij2 - rho_js2)
      term2  = 2.0_wp/rho_i0
      term3  = 1.0_wp/(2.0_wp*rij)*log((rho_js-rij)/(rho_js+rij))
      Hij  = Hij + 0.5_wp*(term1 + term2 + term3)
    else
      Hij  = Hij + 0.0_wp
    end if
  end subroutine compute_Hij



  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_born_radius
  !> @brief        calculate the born radius for the GB model
  !! @authors      TM
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    pairlist    : pairlist information
  !! @param[in]    coord       : coordinates of target systems
  !! @note         A.Onufriev et al., Proteins, 55, 383-394 (2004)
  !!               J.Srinivasan et al., Theor.Chem.Acc., 101, 426-434 (1999)
  !!               M.Schaefer & C.Froemmel, J.Mol.Biol., 216, 1045-1066 (1990)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_born_radius(enefunc, molecule, pairlist, coord)

    ! formal arguments
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_molecule), target, intent(in)    :: molecule
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)

    ! local variables
    integer                   :: i, j, k
    integer                   :: natom
    integer                   :: num_all, ini_all, fin_all
    integer                   :: id, my_id
    real(wp)                  :: drij(1:3), rij, rij2, tmp
    real(wp)                  :: a, b, c, rho_0
    real(wp)                  :: rho_i0, rho_is, rho_j0, rho_js
    real(wp)                  :: rho_js2, rho_is2, rho_i02, rho_j02
    real(wp)                  :: inv_alpha 
    real(wp)                  :: PSI1, PSI2, PSI3, b2, c3
    real(wp)                  :: Hij_i
    real(wp)                  :: term1, term2, term3
    real(wp)                  :: term_tanh, term_exp1
    real(wp)                  :: cutoff2
    real(wp)                  :: max_rho_is
    real(wp)                  :: coord_i(3)
    integer,          pointer :: num_all_calc(:,:), all_calc_list(:,:)
    real(wp),         pointer :: rho_k(:), scale_k(:)
    real(wp),         pointer :: cutoff
    real(wp),         pointer :: born_radius(:), dadrfac(:)
#ifdef OMP
    integer                   :: omp_get_thread_num
#endif


    natom    =  molecule%num_atoms
    cutoff   => enefunc%gbsa%cutoffdist
    rho_k    => enefunc%gbsa%vdw_radius
    scale_k  => enefunc%gbsa%scale_factor
    a        =  enefunc%gbsa%sfactor_alpha
    b        =  enefunc%gbsa%sfactor_beta
    c        =  enefunc%gbsa%sfactor_gamma
    rho_0    =  enefunc%gbsa%vdw_offset
    b2       =  2.0_wp*b
    c3       =  3.0_wp*c
    cutoff2  =  cutoff*cutoff
    num_all_calc  => pairlist%num_all_calc
    all_calc_list => pairlist%all_calc_list
    born_radius   => enefunc%gbsa%born_radius
    dadrfac       => enefunc%gbsa%dadrfac
    enefunc%gbsa%tmp_calc = .true.


    ! initial setup
    !
    if (.not. allocated(Hij)) then
      allocate(Hij(natom), sgrad(natom), dfact(natom), tmp_allreduce(natom))

      max_rho_is = -9999_wp
      do i = 1, natom
        rho_i0  = rho_k(i) - rho_0
        rho_is  = rho_i0 * scale_k(i)
        if (rho_is > max_rho_is) max_rho_is = rho_is
      end do
      gbcutoff2 = (cutoff + max_rho_is)**2

    end if

    ! calculate pairwise descreening function H
    !
    num_all = 0
    Hij(:)  = 0.0_wp

    !$omp parallel                                                   &
    !$omp firstprivate(num_all)                                      &
    !$omp private(id, my_id, ini_all, fin_all, i, k, j,              &
    !$omp         rho_i0, rho_is, rho_i02, rho_is2, Hij_i, coord_i,  &
    !$omp         rho_j0, rho_js, rho_j02, rho_js2, drij, rij2, rij, &
    !$omp         term1, term2, term3, tmp)                          &
    !$omp reduction(+:Hij)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    my_id = my_city_rank * nthread + id

    do i = 1, natom-1

      ini_all = num_all + 1
      fin_all = num_all + num_all_calc(i,id+1)
      num_all = fin_all

      if (mod(i-1,nproc_city*nthread) /= my_id) cycle

      rho_i0  = rho_k(i) - rho_0
      rho_is  = rho_i0 * scale_k(i)

      coord_i(1:3) = coord(1:3,i)

      do k = ini_all, fin_all

        j = all_calc_list(k,id+1)

        drij(1:3) = coord_i(1:3) - coord(1:3,j)
        rij2      = drij(1)*drij(1) + drij(2)*drij(2) + drij(3)*drij(3)

        if (rij2 < gbcutoff2) then
          rij     = sqrt(rij2)
          rho_j0  = rho_k(j) - rho_0
          rho_js  = rho_j0 * scale_k(j)
          call compute_Hij(rij, cutoff, rho_i0, rho_js, Hij(i))
          call compute_Hij(rij, cutoff, rho_j0, rho_is, Hij(j))
        end if

      end do

    end do
    !$omp end parallel

#ifdef HAVE_MPI_GENESIS
    tmp_allreduce(1:natom) = 0.0_wp
    tmp_allreduce(1:natom) = Hij(1:natom)
    call mpi_allreduce(tmp_allreduce, Hij, natom, &
                       mpi_wp_real, mpi_sum, mpi_comm_city, ierror)
#endif

    ! compute effective Born radius of atom i and gradient factor
    !
    !$omp parallel do                                                  &
    !$omp private(i, rho_i0, PSI1, PSI2, PSI3,                         &
    !$omp         term1, term2, term3, term_tanh, inv_alpha)
    !
    do i = 1, natom
      rho_i0           = rho_k(i) - rho_0
      PSI1             = rho_i0 * Hij(i)
      PSI2             = PSI1 * PSI1
      PSI3             = PSI1 * PSI2
      term_tanh        = tanh(a*PSI1 - b*PSI2 + c*PSI3)
      inv_alpha        = 1.0_wp/rho_i0 - 1.0_wp/rho_k(i)*term_tanh
      born_radius(i)   = 1.0_wp/inv_alpha

      term1      = born_radius(i)*rho_i0/rho_k(i)
      term2      = 1.0_wp - term_tanh*term_tanh
      term3      = a - b2*PSI1 + c3*PSI2
      dadrfac(i) = term1*term2*term3
    end do
    !$omp end parallel do
    return

  end subroutine compute_born_radius


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_gb
  !> @brief        calculate solvation free energy with the GB model
  !! @authors      TM
  !! @param[in]    molecule : molecule information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] enefunc  : potential energy functions information
  !! @param[inout] force    : forces of target systems
  !! @param[inout] esolv    : solvation free energy
  !! @note         A.Onufriev et al., Proteins, 55, 383-394 (2004)
  !!               J.Srinivasan et al., Theor.Chem.Acc., 101, 426-434 (1999)
  !!               M.Schaefer & C.Froemmel, J.Mol.Biol., 216, 1045-1066 (1990)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_gb(enefunc, molecule, pairlist, coord, force, esolv)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    type(s_enefunc),  target, intent(inout) :: enefunc
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: esolv

    ! local variables
    integer                   :: i, j, k
    integer                   :: natom
    integer                   :: num_all, ini_all, fin_all
    integer                   :: id, my_id, istart, iend
    real(wp)                  :: esolv_omp(nthread)
    real(wp)                  :: drij(1:3), rij, rij2, esolv_tmp, tmp
    real(wp)                  :: dfact_i
    real(wp)                  :: kappa, rho_0, ieps_p, ieps_w
    real(wp)                  :: rho_i0, rho_is, rho_j0, rho_js
    real(wp)                  :: rho_js2, rho_is2, rho_i02, rho_j02
    real(wp)                  :: alpha_i, alpha_j, alpha_ij
    real(wp)                  :: fij, Dij 
    real(wp)                  :: dHdr_i, dHdr_j, grad_coef
    real(wp)                  :: term1, term2, term3, term4
    real(wp)                  :: term_exp1, term_exp2
    real(wp)                  :: elec_term, cutoff2, ELECOEF2
    real(wp)                  :: charge_i, coord_i(3), force_i(3)
    integer,          pointer :: num_all_calc(:,:), all_calc_list(:,:)
    real(wp),         pointer :: charge(:), rho_k(:), scale_k(:)
    real(wp),         pointer :: cutoff
    real(wp),         pointer :: born_radius(:), dadrfac(:)
#ifdef OMP
    integer                   :: omp_get_thread_num
#endif

    natom    =  molecule%num_atoms
    charge   => molecule%charge
    cutoff   => enefunc%gbsa%cutoffdist
    rho_k    => enefunc%gbsa%vdw_radius
    scale_k  => enefunc%gbsa%scale_factor
    rho_0    =  enefunc%gbsa%vdw_offset
    istart   =  enefunc%gbsa%istart_gb
    iend     =  enefunc%gbsa%iend_gb
    ieps_w   =  1.0_wp/enefunc%gbsa%eps_solvent
    ieps_p   =  1.0_wp/enefunc%gbsa%eps_solute
    kappa    =  50.29216_wp*sqrt(enefunc%gbsa%salt_cons*ieps_w/enefunc%gbsa%temperature)
    cutoff2  =  cutoff*cutoff
    ELECOEF2 =  ELECOEF/2.0_wp
    num_all_calc  => pairlist%num_all_calc
    all_calc_list => pairlist%all_calc_list

    if (.not. enefunc%gbsa%tmp_calc) then
        call compute_born_radius(enefunc, molecule, pairlist, coord)
    end if
    enefunc%gbsa%tmp_calc = .false.
    born_radius => enefunc%gbsa%born_radius
    dadrfac => enefunc%gbsa%dadrfac

    !$omp parallel do                                                  &
    !$omp private(i, term1, term2, term_exp1)
    !
    do i = 1, natom
      term_exp1  = exp(-kappa*born_radius(i))*ieps_w
      term1      = ieps_p - term_exp1
      term2      = ELECOEF2*charge(i)*charge(i)/born_radius(i)
      sgrad(i)   = term1*term2 - kappa*term2*born_radius(i)*term_exp1
    end do
    !$omp end parallel do

    ! compute self energy
    !
    esolv_omp(1:nthread) = 0.0_wp
    !$omp parallel do              &
    !$omp private(i, term1, term2) &
    !$omp reduction(+:esolv_omp)
    !
    do i = istart, iend
      term1 = ieps_p - exp(-kappa*born_radius(i))*ieps_w
      term2 = ELECOEF2*charge(i)*charge(i)/born_radius(i)
      esolv_omp(1) = esolv_omp(1) - term1*term2
    end do
    !$omp end parallel do

    ! compute solvation free energy
    !
    dfact(:) = 0.0_wp
    num_all  = 0
    !$omp parallel                                                        &
    !$omp firstprivate(num_all)                                           &
    !$omp private(id, my_id, ini_all, fin_all, i, k, j,                   &
    !$omp         coord_i, drij, rij2, rij, term_exp1, term_exp2, term1,  &
    !$omp         term2, elec_term, fij, Dij, alpha_i, alpha_j, alpha_ij, &
    !$omp         tmp, grad_coef, charge_i, dfact_i, force_i, esolv_tmp)  &
    !$omp reduction(+:dfact)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    my_id = my_city_rank * nthread + id

    do i = 1, natom-1

      ini_all = num_all + 1
      fin_all = num_all + num_all_calc(i,id+1)
      num_all = fin_all

      if (mod(i-1,nproc_city*nthread) /= my_id) cycle

      alpha_i      = born_radius(i)
      charge_i     = charge(i)*ELECOEF
      coord_i(1:3) = coord(1:3,i)
      esolv_tmp    = 0.0_wp
      dfact_i      = 0.0_wp
      force_i(1:3) = 0.0_wp

      do k = ini_all, fin_all

        j = all_calc_list(k,id+1)

        drij(1:3) = coord_i(1:3) - coord(1:3,j)
        rij2      = drij(1)*drij(1) + drij(2)*drij(2) + drij(3)*drij(3)

        if (rij2 < gbcutoff2) then
          rij       = sqrt(rij2)
          alpha_j   = born_radius(j)
          alpha_ij  = alpha_i * alpha_j

          ! effective distance
          term_exp1 = exp(-rij2/(4.0_wp*alpha_ij))
          fij       = sqrt(rij2 + alpha_ij*term_exp1)

          ! solvation free energy
          term_exp2 = exp(-kappa*fij)
          Dij       = ieps_p - term_exp2*ieps_w
          elec_term = charge_i*charge(j)
          esolv_tmp = esolv_tmp - elec_term*Dij/fij

          ! gradient
          tmp       = - elec_term/(2.0_wp*fij*fij)
          term1     = tmp*(kappa*ieps_w*term_exp2 - Dij/fij)
          term2     = term1*(alpha_ij + 0.25_wp*rij2)*term_exp1

          dfact_i   = dfact_i  + term2
          dfact(j)  = dfact(j) + term2

          grad_coef         = term1*(2.0_wp - 0.5_wp*term_exp1)
          force_i(1:3)      = force_i(1:3)      - grad_coef*drij(1:3)
          force(1:3,j,id+1) = force(1:3,j,id+1) + grad_coef*drij(1:3)
        end if
      end do

      esolv_omp(id+1)   = esolv_omp(id+1)   + esolv_tmp
      force(1:3,i,id+1) = force(1:3,i,id+1) + force_i(1:3)
      dfact(i)          = dfact(i)          + dfact_i

    end do
    !$omp end parallel

    esolv = 0.0_wp
    do i = 1, nthread
      esolv = esolv + esolv_omp(i)
    end do

#ifdef HAVE_MPI_GENESIS
    tmp_allreduce(1:natom) = 0.0_wp
    tmp_allreduce(1:natom) = dfact(1:natom)
    call mpi_allreduce(tmp_allreduce, dfact, natom, &
                       mpi_wp_real, mpi_sum, mpi_comm_city, ierror)
#endif

    dfact(:) = dfact(:) + sgrad(:)

    ! compute solvation free energy and force
    !
    num_all = 0
    !$omp parallel                                                    &
    !$omp firstprivate(num_all)                                       &
    !$omp private(id, my_id, ini_all, fin_all, i, k, j,               &
    !$omp         rho_i0, rho_is, rho_i02, rho_is2, coord_i, force_i, &
    !$omp         rho_j0, rho_js, rho_j02, rho_js2, drij, rij2, rij,  &
    !$omp         term1, term2, term3, term4, term_exp1, term_exp2,   &
    !$omp         dHdr_i, dHdr_j, grad_coef, tmp)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    my_id = my_city_rank * nthread + id

    do i = 1, natom-1

      ini_all = num_all + 1
      fin_all = num_all + num_all_calc(i,id+1)
      num_all = fin_all

      if (mod(i-1,nproc_city*nthread) /= my_id) cycle

      rho_i0  = rho_k(i) - rho_0
      rho_is  = rho_i0 * scale_k(i)
      rho_i02 = rho_i0 * rho_i0
      rho_is2 = rho_is * rho_is
      coord_i(1:3) = coord(1:3,i)
      force_i(1:3) = 0.0_wp

      do k = ini_all, fin_all
        j = all_calc_list(k,id+1)

        drij(1:3) = coord_i(1:3) - coord(1:3,j)
        rij2      = drij(1)*drij(1) + drij(2)*drij(2) + drij(3)*drij(3)

        if (rij2 < gbcutoff2) then
          rij     = sqrt(rij2)
          rho_j0  = rho_k(j) - rho_0
          rho_js  = rho_j0 * scale_k(j)
          rho_j02 = rho_j0 * rho_j0
          rho_js2 = rho_js * rho_js

          ! dHij/drij
          if (rij > cutoff + rho_js) then
            dHdr_i = 0.0_wp
          else if (rij > cutoff - rho_js) then
            tmp    = rij - rho_js
            term1  = - (cutoff - tmp)*(cutoff + tmp)*(rho_js2 + rij2)
            term2  = 8.0_wp*cutoff2*rij2*tmp*tmp
            term3  = - 0.25_wp*log(tmp/cutoff)/rij2
            dHdr_i = term1/term2 + term3
          else if (rij > 4.0_wp*rho_js) then
            tmp    = rho_js2/rij2
            term1  = paa + tmp*(pbb + tmp*(pcc + tmp*(pdd + tmp*pee)))
            term2  = - tmp*tmp/rho_js/rij
            dHdr_i = term1*term2
          else if (rij > rho_i0 + rho_js) then
            term1  = - rho_js*(rij2 + rho_js2)/(rij*(rij2 - rho_js2)**2)
            term2  = - 0.5_wp/rij2 * log((rij - rho_js)/(rij + rho_js))
            dHdr_i = 0.5_wp * (term1 + term2)
          else if (rij > abs(rho_i0 - rho_js)) then
            term1  = - 0.5_wp/rho_i02
            term2  = (rij2+rho_js2)*(rho_i02-rho_js2) - 2.0_wp*rij*rho_js**3
            term2  = term2/(2.0_wp*rij2*rho_i02*(rij + rho_js)**2)
            term3  = - 1.0_wp/rij2 * log(rho_i0/(rij + rho_js))
            dHdr_i = 0.25_wp * (term1 + term2 + term3)
          else if (rho_i0 < rho_js) then
            term1  = - rho_js*(rij2 + rho_js2)/(rij*(rij2 - rho_js2)**2)
            term2  = - 0.5_wp/rij2 * log((rho_js - rij)/(rij + rho_js))
            dHdr_i = 0.5_wp * (term1 + term2)
          else
            dHdr_i = 0.0_wp
          end if

          ! dHji/drij
          if (rij > cutoff + rho_is) then
            dHdr_j = 0.0_wp
          else if (rij > cutoff - rho_is) then
            tmp    = rij - rho_is
            term1  = - (cutoff - tmp)*(cutoff + tmp)*(rho_is2 + rij2)
            term2  = 8.0_wp*cutoff2*rij2*tmp*tmp
            term3  = - 0.25_wp*log(tmp/cutoff)/rij2
            dHdr_j = term1/term2 + term3
          else if (rij > 4.0_wp*rho_is) then
            tmp    = rho_is2/rij2
            term1  = paa + tmp*(pbb + tmp*(pcc + tmp*(pdd + tmp*pee)))
            term2  = - tmp*tmp/rho_is/rij
            dHdr_j = term1*term2
          else if (rij > rho_j0 + rho_is) then
            term1  = - rho_is*(rij2 + rho_is2)/(rij*(rij2 - rho_is2)**2)
            term2  = - 0.5_wp/rij2 * log((rij - rho_is)/(rij + rho_is))
            dHdr_j = 0.5_wp * (term1 + term2)
          else if (rij > abs(rho_j0 - rho_is)) then
            term1  = - 0.5_wp/rho_j02
            term2  = (rij2+rho_is2)*(rho_j02-rho_is2) - 2.0_wp*rij*rho_is**3
            term2  = term2/(2.0_wp*rij2*rho_j02*(rij + rho_is)**2)
            term3  = - 1.0_wp/rij2 * log(rho_j0/(rij + rho_is))
            dHdr_j = 0.25_wp * (term1 + term2 + term3)
          else if (rho_j0 < rho_is) then
            term1  = - rho_is*(rij2 + rho_is2)/(rij*(rij2 - rho_is2)**2)
            term2  = - 0.5_wp/rij2 * log((rho_is - rij)/(rij + rho_is))
            dHdr_j = 0.5_wp * (term1 + term2)
          else
            dHdr_j = 0.0_wp
          end if

          term1 = dadrfac(i)*dHdr_i*dfact(i)
          term2 = dadrfac(j)*dHdr_j*dfact(j)
          grad_coef = (term1 + term2)/rij

          force_i(1:3)      = force_i(1:3)      - grad_coef*drij(1:3)
          force(1:3,j,id+1) = force(1:3,j,id+1) + grad_coef*drij(1:3)
        end if

      end do
      force(1:3,i,id+1) = force(1:3,i,id+1) + force_i(1:3)

    end do
    !$omp end parallel

    return

  end subroutine compute_energy_gb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_SA
  !> @brief        calculate solvation free energy with the SA model
  !! @authors      TM
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] esolv    : solvation free energy
  !! @note         J.Weiser et al., J.Comput.Chem., 20, 217-230 (1999)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_sa(enefunc, molecule, coord, force, esolv)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_enefunc),  target, intent(in)    :: enefunc
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: esolv

    ! local variables
    integer                   :: i, j, k, is, natom, natom_sasa
    integer                   :: c, m, n, ni, nj, mi, mj
    integer                   :: mcx, mcy, mcz, gmcx, gmcy, gmcz
    integer                   :: ncell_x, ncell_y, ncell_z, ncell_xyz
    integer                   :: ncell_yz
    integer                   :: id, my_id
    real(wp)                  :: P1, P2, P3, P4
    real(wp)                  :: term1, term2, term3, term4
    real(wp)                  :: Ri, Rj, Rk, sasa_i, sasa, tension
    real(wp)                  :: rij2, rjk2, rik2
    real(wp)                  :: Ri2, Rj2, Rk2, RjRk2, RiRk2, RiRj2
    real(wp)                  :: rij, rjk, Aij, Ajk
    real(wp)                  :: dxij(3), dxjk(3), dxik(3)
    real(wp)                  :: dAijdx(3), dAjkdx(3)
    real(wp)                  :: dfacij(3), dfacjk(3)
    real(wp)                  :: sum_Ajk, sum_dfacij(3), sum_dfacjk(3)
    real(wp)                  :: box_size_x, box_size_y, box_size_z
    real(wp)                  :: minx, miny, minz, maxx, maxy, maxz
    real(wp)                  :: csize_x, csize_y, csize_z
    real(wp)                  :: overlap_dist, fact
    real(wp)                  :: coord_i(3), coord_j(3)
    logical                   :: zero_cell
    integer,  pointer         :: atom_type(:), atom_list(:)
    real(wp), pointer         :: radius(:), param(:,:)
#ifdef OMP
    integer                   :: omp_get_thread_num
#endif


    natom      =  molecule%num_atoms
    natom_sasa =  enefunc%gbsa%num_sasa_atoms
    atom_type  => enefunc%gbsa%sasa_atom_type
    param      => enefunc%gbsa%sasa_parameters
    radius     => enefunc%gbsa%sasa_vdw_radius
    atom_list  => enefunc%gbsa%sasa_atom_list
    tension    =  enefunc%gbsa%surface_tension


    if (.not. allocated(head_atom)) then
      allocate(cell_linked_list(natom), cell_pointer(natom))
      ncell_xyz_old = 0
    end if

    overlap_dist = enefunc%gbsa%pairlist_distance

    maxx = maxval(coord(1,1:natom))
    minx = minval(coord(1,1:natom))
    maxy = maxval(coord(2,1:natom))
    miny = minval(coord(2,1:natom))
    maxz = maxval(coord(3,1:natom))
    minz = minval(coord(3,1:natom))

    box_size_x = (maxx - minx) + 2.0_wp*overlap_dist
    box_size_y = (maxy - miny) + 2.0_wp*overlap_dist
    box_size_z = (maxz - minz) + 2.0_wp*overlap_dist

    ncell_x = int(box_size_x/overlap_dist)
    ncell_y = int(box_size_y/overlap_dist)
    ncell_z = int(box_size_z/overlap_dist)

    csize_x = box_size_x/dble(ncell_x)
    csize_y = box_size_y/dble(ncell_y)
    csize_z = box_size_z/dble(ncell_z)

    minx = minx - csize_x
    miny = miny - csize_y
    minz = minz - csize_z

    ncell_xyz = ncell_x*ncell_y*ncell_z
    ncell_yz  = ncell_y*ncell_z

    if (ncell_xyz > ncell_xyz_old) then
      if (allocated(head_atom)) deallocate(cell_list, head_atom)
      allocate(cell_list(0:ncell_xyz,27), head_atom(0:ncell_xyz))
    end if
    ncell_xyz_old = ncell_xyz


    cell_list(:,:) = 0

    do mcx = 0, ncell_x-1
      do mcy = 0, ncell_y-1
        do mcz = 0, ncell_z-1
          c = 1 + mcx*ncell_yz + mcy*ncell_z + mcz
          m = 0
          do i = -1, 1
            do j = -1, 1
              do k = -1, 1
                m = m + 1
                zero_cell = .false.
                gmcx = mcx + i
                if (gmcx == -1)      zero_cell = .true.
                if (gmcx == ncell_x) zero_cell = .true.
                gmcy = mcy + j
                if (gmcy == -1)      zero_cell = .true.
                if (gmcy == ncell_y) zero_cell = .true.
                gmcz = mcz + k
                if (gmcz == -1)      zero_cell = .true.
                if (gmcz == ncell_z) zero_cell = .true.
                if (.not. zero_cell) then
                  n = 1 + gmcx*ncell_yz + gmcy*ncell_z + gmcz
                else
                  n = 0
                end if
                cell_list(c,m) = n
              end do
            end do
          end do
        end do
      end do
    end do

    head_atom(0:ncell_xyz) = 0
    do is = 1, natom_sasa
      i = atom_list(is)
      mcx = aint((coord(1,i) - minx)/csize_x)
      mcy = aint((coord(2,i) - miny)/csize_y)
      mcz = aint((coord(3,i) - minz)/csize_z)
      c = 1 + mcx*ncell_yz + mcy*ncell_z + mcz
      cell_pointer(i) = c
      cell_linked_list(i) = head_atom(c)
      head_atom(c) = i
    end do


    ! calculate SA term
    !
    sasa = 0.0_wp

    !$omp parallel                                                      &
    !$omp private(id, my_id, i, k, j, is, mi, mj, P1, P2, P3, P4,       &
    !$omp         Ri, Ri2, term1, term2, term3, term4, mcx, mcy, mcz,   &
    !$omp         ni, Rj, dxij, rij2, RiRj2, Rj2, rij, Aij, sum_Ajk,    &
    !$omp         sum_dfacjk, dxik, rik2, RiRk2, Rk, dxjk, rjk2, RjRk2, &
    !$omp         Rk2, rjk, fact, Ajk, dAjkdx, nj, dfacjk, sasa_i,      &
    !$omp         dfacij, dAijdx, sum_dfacij, coord_i, coord_j)         &
    !$omp reduction(+:sasa)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    my_id = my_city_rank * nthread + id

    do is = 1, natom_sasa

      if (mod(is-1,nproc_city*nthread) /= my_id) cycle

      i = atom_list(is)

      P1 = param(atom_type(i),1)
      P2 = param(atom_type(i),2)
      P3 = param(atom_type(i),3)
      P4 = param(atom_type(i),4)

      Ri  = radius(i)
      Ri2 = Ri*Ri

      term1 = 4.0_wp*PI*Ri2
      term2 = 0.0_wp
      term3 = 0.0_wp
      term4 = 0.0_wp

      coord_i(1) = coord(1,i)
      coord_i(2) = coord(2,i)
      coord_i(3) = coord(3,i)

      ni = cell_pointer(i)

      sum_dfacij(1:3) = 0.0_wp

      do mi = 1, 27

        j = head_atom(cell_list(ni,mi))

        do while (j /= 0)

          if (j == i) then
            j = cell_linked_list(j)
            cycle
          end if

          coord_j(1) = coord(1,j)
          coord_j(2) = coord(2,j)
          coord_j(3) = coord(3,j)
          Rj = radius(j)

          dxij(1) = coord_i(1) - coord_j(1)
          dxij(2) = coord_i(2) - coord_j(2)
          dxij(3) = coord_i(3) - coord_j(3)
          rij2  = dxij(1)*dxij(1) + dxij(2)*dxij(2) + dxij(3)*dxij(3)
          RiRj2 = Ri + Rj
          RiRj2 = RiRj2*RiRj2
          if (rij2 > RiRj2) then
            j = cell_linked_list(j)
            cycle
          end if

          Rj2   = Rj*Rj
          rij   = sqrt(rij2)
          Aij   = PI*Ri*(2.0_wp*Ri - rij - (Ri2 - Rj2)/rij)
          term2 = term2 + Aij

          dAijdx(1:3) = PI*Ri*((Ri2 - Rj2)/rij2 - 1.0_wp)*dxij(1:3)/rij

          nj = cell_pointer(j)

          sum_Ajk = 0.0_wp
          sum_dfacjk(1:3) = 0.0_wp

          do mj = 1, 27

            k = head_atom(cell_list(nj,mj))

            do while (k /= 0)

              if (k == j .or. k == i) then
                k = cell_linked_list(k)
                cycle
              end if

              Rk = radius(k)

              dxik(1) = coord_i(1) - coord(1,k)
              dxik(2) = coord_i(2) - coord(2,k)
              dxik(3) = coord_i(3) - coord(3,k)
              rik2  = dxik(1)*dxik(1) + dxik(2)*dxik(2) + dxik(3)*dxik(3)
              RiRk2 = Ri + Rk
              RiRk2 = RiRk2*RiRk2
              if (rik2 > RiRk2) then
                k = cell_linked_list(k)
                cycle
              end if

              dxjk(1) = coord_j(1) - coord(1,k)
              dxjk(2) = coord_j(2) - coord(2,k)
              dxjk(3) = coord_j(3) - coord(3,k)
              rjk2  = dxjk(1)*dxjk(1) + dxjk(2)*dxjk(2) + dxjk(3)*dxjk(3)
              RjRk2 = Rj + Rk
              RjRk2 = RjRk2*RjRk2
              if (rjk2 > RjRk2) then
                k = cell_linked_list(k)
                cycle
              end if

              Rk2  = Rk*Rk
              rjk  = sqrt(rjk2)
              fact = (Rj2 - Rk2)/rjk2
              Ajk  = PI*Rj*(2.0_wp*Rj - rjk - fact*rjk)

              sum_Ajk = sum_Ajk + Ajk

              dAjkdx(1:3) = PI*Rj*(fact - 1.0_wp)*dxjk(1:3)/rjk
              dfacjk(1:3) = (P3 + P4*Aij)*dAjkdx(1:3)

              sum_dfacjk(1:3) = sum_dfacjk(1:3) + dfacjk(1:3)

              force(1:3,k,id+1) = force(1:3,k,id+1) + tension*dfacjk(1:3)

              k = cell_linked_list(k)

            end do

          end do

          term3 = term3 + sum_Ajk
          term4 = term4 + sum_Ajk*Aij

          dfacij    (1:3) = P2*dAijdx(1:3) + P4*dAijdx(1:3)*sum_Ajk
          sum_dfacij(1:3) = sum_dfacij(1:3) + dfacij(1:3)

          force(1:3,j,id+1) = force(1:3,j,id+1) + tension*(dfacij(1:3) - sum_dfacjk(1:3))

          j = cell_linked_list(j)

        end do

      end do

      sasa_i = P1*term1 + P2*term2 + P3*term3 + P4*term4
      sasa = sasa + sasa_i

      force(1:3,i,id+1) = force(1:3,i,id+1) - tension*sum_dfacij(1:3)

    end do
    !$omp end parallel

    esolv = esolv + tension*sasa

    return

  end subroutine compute_energy_sa

end module at_energy_gbsa_mod
