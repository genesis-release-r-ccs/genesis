!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_pairlist_mod
!> @brief   set pairlist for nonbonded interactions
!! @authors Jaewoon Jung (JJ), Takaharu Mori (TM), Takashi Imai (TI), 
!!          Motoshi Kamiya (MK), Yuji Sugita (YS), Kiyoshi Yagi (KY),
!!          Cheng Tan (CT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_pairlist_mod

  use at_boundary_str_mod
  use at_pairlist_str_mod
  use at_enefunc_str_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  real(wp),         parameter :: TableWaterMargin = 2.0_wp
  real(wp),         parameter :: FactNumNb15   = 1.5_wp
  real(wp),         parameter :: FactThreshold = 0.95_wp

  ! subroutines
  public  :: setup_pairlist
  public  :: update_pairlist
  private :: update_pairlist_nobc
  private :: update_pairlist_nobc_gbsa
  private :: update_pairlist_nobc_cg
  private :: update_pairlist_nobc_cg_pwmcos
  private :: update_pairlist_nobc_cg_pwmcosns
  private :: update_pairlist_nobc_cg_IDR_HPS
  private :: update_pairlist_nobc_cg_IDR_KH
  private :: update_pairlist_pbc_cg_general
  private :: update_pairlist_pbc_cg_exv
  private :: update_pairlist_pbc_cg_ele
  private :: update_pairlist_pbc_cg_DNA_base
  private :: update_pairlist_pbc_cg_DNA_exv
  private :: update_pairlist_pbc_cg_pwmcos
  private :: update_pairlist_pbc_cg_pwmcosns
  private :: update_pairlist_pbc_cg_IDR_HPS
  private :: update_pairlist_pbc_cg_IDR_KH
  private :: update_pairlist_pbc_cg_KH
  private :: update_ecqm15_nonb
!  private :: update_pairlist_pbc
!  private :: update_pairlist_pbc_cell
  private :: update_pairlist_solute_solute
  private :: update_pairlist_solute_water
  private :: update_pairlist_water_water
  private :: check_pairlist_memory_size

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_pairlist
  !> @brief        initialize/allocate/setup pairlist for nonbonded interactions
  !! @authors      YS, TI, JJ, TM
  !! @param[in]    enefunc  : information of potential energy function
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    coord    : atomic coordinate
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_pairlist(enefunc, boundary, coord, trans, coord_pbc, &
                            pairlist)

    ! formal arguments
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_boundary),        intent(in)    :: boundary
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: trans(:,:)
    real(wp),                intent(inout) :: coord_pbc(:,:)
    type(s_pairlist),        intent(inout) :: pairlist
    
    ! local variables
    integer                  :: natom, nthread
    integer                  :: n_cg_DNA_base, n_cg_DNA_phos, n_cg_DNA_all
    integer                  :: n_cg_charged, n_cg_pro_charged
    integer                  :: n_cg_IDR_HPS
    integer                  :: n_cg_IDR_KH
    integer                  :: n_cg_KH
    integer                  :: n_cg_PWMcos
    integer                  :: n_cg_PWMcosns
    integer                  :: n_cg_cells
#ifdef OMP
    integer                  :: omp_get_max_threads
#endif


    ! do not make pairlist, if my_node does not compute 
    !   non-bonded energy in the real part
    !
    if (.not. real_calc) &
      return

    ! initialize pairlist
    !
    call init_pairlist(pairlist)

    pairlist%pairlistdist           = enefunc%pairlistdist
    ! CG
    pairlist%cg_pairlistdist_ele    = enefunc%cg_pairlistdist_ele
    pairlist%cg_pairlistdist_126    = enefunc%cg_pairlistdist_126
    pairlist%cg_pairlistdist_PWMcos = enefunc%cg_pairlistdist_PWMcos
    pairlist%cg_pairlistdist_DNAbp  = enefunc%cg_pairlistdist_DNAbp
    pairlist%cg_pairlistdist_exv    = enefunc%cg_pairlistdist_exv

    ! allocate pairlist
    !
#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif
    call alloc_pairlist(pairlist, PairListNthreads, nthread)
    if (enefunc%gbsa_use) then
      call alloc_pairlist(pairlist, PairListNthreadsGbsa, nthread)
    end if
    if (enefunc%forcefield == ForcefieldRESIDCG) then
      call alloc_pairlist(pairlist, PairListNthreadsCG, nthread)
      if (enefunc%num_pwmcos_resid > 0) then
        call alloc_pairlist(pairlist, PairListNthCGPWMcos, enefunc%num_pwmcos_resid)
      end if
      if (enefunc%num_pwmcosns_resid > 0) then
        call alloc_pairlist(pairlist, PairListNthCGPWMcosns, enefunc%num_pwmcosns_resid)
      end if
    end if
    
    natom = size(coord(1,:))

    select case (boundary%type)

    case (BoundaryTypeNOBC)

      call alloc_pairlist(pairlist, PairListAtomNobc, natom)
      if (enefunc%gbsa_use) then
        call alloc_pairlist(pairlist, PairListAtomNobcGbsa, natom)
      end if
      if (enefunc%forcefield == ForcefieldRESIDCG) then
        call alloc_pairlist(pairlist, PairListAtomNobcCG, natom)
      end if

    case (BoundaryTypePBC)

      if (enefunc%forcefield == ForcefieldRESIDCG) then
        n_cg_DNA_all     = enefunc%num_cg_particle_DNA_all
        n_cg_DNA_phos    = enefunc%num_cg_particle_DNA_phos
        n_cg_DNA_base    = enefunc%num_cg_particle_DNA_base
        n_cg_IDR_HPS     = enefunc%num_cg_particle_IDR_HPS
        n_cg_IDR_KH      = enefunc%num_cg_particle_IDR_KH
        n_cg_KH          = enefunc%num_cg_particle_KH
        n_cg_charged     = enefunc%num_cg_particle_charged
        n_cg_pro_charged = enefunc%num_cg_particle_pro_charged
        n_cg_PWMcos      = enefunc%num_pwmcos_terms
        n_cg_PWMcosns    = enefunc%num_pwmcosns_terms
        n_cg_cells       = boundary%num_cells
        call alloc_pairlist(pairlist, PairListNthreadsPbcCG, nthread)
        call alloc_pairlist(pairlist, PairListAtomPbcCGexv, natom)
        call alloc_pairlist(pairlist, PairListAtomPbcCGele, n_cg_charged, n_cg_DNA_phos)
        call alloc_pairlist(pairlist, PairListAtomPbcCGDNAbp, n_cg_DNA_base)
        call alloc_pairlist(pairlist, PairListAtomPbcCGDNAexv, n_cg_DNA_all)
        call alloc_pairlist(pairlist, PairListAtomPbcCGPWMcos, n_cg_PWMcos)
        call alloc_pairlist(pairlist, PairListAtomPbcCGPWMcosns, n_cg_PWMcosns)
        call alloc_pairlist(pairlist, PairListAtomPbcCGIDRHPS, n_cg_IDR_HPS)
        call alloc_pairlist(pairlist, PairListAtomPbcCGIDRKH, n_cg_IDR_KH)
        call alloc_pairlist(pairlist, PairListAtomPbcCGKH, n_cg_KH)
        call alloc_pairlist(pairlist, PairListCellsPbcCG, n_cg_cells)
      else
        call alloc_pairlist(pairlist, PairListPbcSolute, enefunc%table%num_solute)
        call alloc_pairlist(pairlist, PairListPbcWater,  enefunc%table%num_water)
      end if

    end select

    ! make pairlist
    !
    pairlist%allocation = .true.
    call update_pairlist(enefunc, boundary, coord, trans, coord_pbc, pairlist)
    pairlist%allocation = .false.

    return

  end subroutine setup_pairlist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist
  !> @brief        update pairlist for nonbonded interactions
  !! @authors      YS, TM, TI, JJ, CT
  !! @param[in]    enefunc  : information of potential energy function
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    coord    : atomic coordinate
  !! @param[inout] pairlist : information of pairlist
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist(enefunc, boundary, coord, trans, coord_pbc, &
                             pairlist)

    ! formal arguments
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_boundary),        intent(in)    :: boundary
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: trans(:,:)
    real(wp),                intent(inout) :: coord_pbc(:,:)
    type(s_pairlist),        intent(inout) :: pairlist


    call timer(TimerPairList, TimerOn)

    if (pairlist%allocation) then
      pairlist%allocate_nobc            = .true.
      pairlist%allocate_pbc             = .true.
      pairlist%allocate_solsol          = .true.
      pairlist%allocate_solwat          = .true.
      pairlist%allocate_watwat          = .true.
      pairlist%allocate_nobc_cg         = .true.
      pairlist%allocate_nobc_cg_pwmcos  = .true.
      pairlist%allocate_nobc_cg_pwmcosns= .true.
      pairlist%allocate_nobc_cg_IDR_HPS = .true.
      pairlist%allocate_nobc_cg_IDR_KH  = .true.
      pairlist%allocate_pbc_cg_exv      = .true.
      pairlist%allocate_pbc_cg_ele      = .true.
      pairlist%allocate_pbc_cg_DNA_bp   = .true.
      pairlist%allocate_pbc_cg_DNA_exv  = .true.
      pairlist%allocate_pbc_cg_pwmcos   = .true.
      pairlist%allocate_pbc_cg_pwmcosns = .true.
      pairlist%allocate_pbc_cg_IDR_HPS  = .true.
      pairlist%allocate_pbc_cg_IDR_KH   = .true.
      pairlist%allocate_pbc_cg_KH       = .true.
    end if

    select case (boundary%type)

    case (BoundaryTypeNOBC)

      if (enefunc%gbsa_use) then
        call update_pairlist_nobc_gbsa(enefunc, coord, pairlist)
      else if (enefunc%forcefield == ForcefieldRESIDCG) then
        call update_pairlist_nobc_cg  (enefunc, coord, pairlist)
        if (enefunc%num_pwmcos_terms > 0) then
          call update_pairlist_nobc_cg_pwmcos(enefunc, coord, pairlist)
        end if
        if (enefunc%num_pwmcosns_terms > 0) then
          call update_pairlist_nobc_cg_pwmcosns(enefunc, coord, pairlist)
        end if
        if (enefunc%cg_IDR_HPS_calc) then
          call update_pairlist_nobc_cg_IDR_HPS(enefunc, coord, pairlist)
        end if
        if (enefunc%cg_IDR_KH_calc) then
          call update_pairlist_nobc_cg_IDR_KH(enefunc, coord, pairlist)
        end if
      else
        call update_pairlist_nobc     (enefunc, coord, pairlist)
      end if

      if (enefunc%qmmm%do_qmmm .and. enefunc%qmmm%num_qmmmbonds > 0)  &
           call update_ecqm15_nonb(enefunc, coord, pairlist)

    case (BoundaryTypePBC)

      ! if (enefunc%table%table) then
      !   call update_pairlist_solute_solute(enefunc, boundary, coord, pairlist)
      !   call update_pairlist_solute_water (enefunc, boundary, coord, pairlist)
      !   call update_pairlist_water_water  (enefunc, boundary, coord, pairlist)
      ! else
      !   if (boundary%use_cell_linked_list) then
      !     call update_pairlist_pbc_cell(enefunc, boundary, coord, pairlist)
      !   else
      !     call update_pairlist_pbc(enefunc, boundary, coord, pairlist)
      !   end if
      ! end if
      if (enefunc%forcefield == ForcefieldRESIDCG) then
        call update_coord_pbc (enefunc, boundary, coord, trans, coord_pbc)

        call update_pairlist_pbc_cg_general    (enefunc, boundary, coord_pbc, pairlist)
        call update_pairlist_pbc_cg_exv        (enefunc, boundary, coord_pbc, pairlist)
        if (enefunc%cg_ele_calc) then
          call update_pairlist_pbc_cg_ele      (enefunc, boundary, coord_pbc, pairlist)
        end if
        if (enefunc%cg_DNA_base_pair_calc) then
          call update_pairlist_pbc_cg_DNA_base (enefunc, boundary, coord_pbc, pairlist)
        end if
        if (enefunc%cg_DNA_exv_calc) then
          call update_pairlist_pbc_cg_DNA_exv  (enefunc, boundary, coord_pbc, pairlist)
        end if
        if (enefunc%cg_pwmcos_calc) then
          call update_pairlist_pbc_cg_pwmcos   (enefunc, boundary, coord_pbc, pairlist)
        end if
        if (enefunc%cg_pwmcosns_calc) then
          call update_pairlist_pbc_cg_pwmcosns   (enefunc, boundary, coord_pbc, pairlist)
        end if
        if (enefunc%cg_IDR_HPS_calc) then
          call update_pairlist_pbc_cg_IDR_HPS  (enefunc, boundary, coord_pbc, pairlist)
        end if
        if (enefunc%cg_IDR_KH_calc) then
          call update_pairlist_pbc_cg_IDR_KH   (enefunc, boundary, coord_pbc, pairlist)
        end if
        if (enefunc%cg_KH_calc) then
          call update_pairlist_pbc_cg_KH       (enefunc, boundary, coord_pbc, pairlist)
        end if
      else
        call update_pairlist_solute_solute(enefunc, boundary, coord, pairlist)
        call update_pairlist_solute_water (enefunc, boundary, coord, pairlist)
        call update_pairlist_water_water  (enefunc, boundary, coord, pairlist)
      end if
    end select

    call timer(TimerPairList, TimerOff)

    return

  end subroutine update_pairlist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_nobc
  !> @brief        update pairlist when no boundary condition is applied
  !! @authors      YS, JJ, TM
  !! @param[in]    enefunc  : information of potential energy function
  !! @param[in]    coord    : atomic coordinate
  !! @param[inout] pairlist : information of pairlist
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_nobc(enefunc, coord, pairlist)
  
    ! formal arguments
    type(s_enefunc),          intent(in)    :: enefunc
    real(wp),                 intent(in)    :: coord(:,:)
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                 :: pairdist2, dij(1:3), rij2
    integer                  :: i, j, k, natom, n, nloops
    integer                  :: num_excl, ini_excl, fin_excl
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: num_nb15_max
    integer                  :: id, my_id, nthread
#ifdef OMP
    integer                  :: omp_get_thread_num, omp_get_max_threads
#endif
    logical                  :: nb15_calc
    logical                  :: proceed, do_allocate

    integer,     pointer     :: num_nb15_pre(:), num_nb15(:)


    num_nb15_pre => pairlist%num_nb15_pre
    num_nb15     => pairlist%num_nb15
    pairdist2    =  pairlist%pairlistdist * pairlist%pairlistdist
    natom        =  size(coord(1,:))

    ! initialize
    !
#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif
    pairlist%num_nb15_calc(1:natom,1:nthread) = 0

    ! if pairlist%allocation = .true. , allocation + make pairlist
    ! if pairlist%allocation = .false., make pairlist only
    !
    if (pairlist%allocate_nobc) then
      nloops = 2
      do_allocate = .true.
    else
      nloops = 1
      do_allocate = .false.
    end if
    
    do n = 1, nloops

      num_excl = 0
      num_nb14 = 0
      num_nb15(1:nthread) = 0

      if (.not. do_allocate) &
        num_nb15_pre(1:nthread) = 0

      !$omp parallel                                                       &
      !$omp private(id, my_id, i, ini_excl, fin_excl, ini_nb14, fin_nb14,  &
      !$omp         proceed, j, nb15_calc, k, dij, rij2)                   &
      !$omp firstprivate(num_excl, num_nb14, do_allocate)
      !
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      my_id = my_city_rank * nthread + id
      id = id + 1

      do i = 1, natom-1

        ini_excl = num_excl + 1
        fin_excl = num_excl + enefunc%num_nonb_excl(i)
        num_excl = fin_excl

        proceed = .true.
        if (mod(i-1,nproc_city*nthread) /= my_id) proceed = .false.

        if (proceed) then
          do j = i+1, natom

            ! compute distance
            !
            dij(1:3) = coord(1:3,i) - coord(1:3,j)
            rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            ! store interaction table
            !
            if (rij2 < pairdist2) then

              ! 1-2 and 1-3 interactions are removed
              !
              nb15_calc = .true.
              do k = ini_excl, fin_excl
                if (j == enefunc%nonb_excl_list(k)) then
                  nb15_calc = .false.
                  exit
                end if
              end do

              ! 1-4 interactions are removed
              !
              if (enefunc%forcefield/=ForcefieldKBGO) then
                do k = 1, enefunc%num_nb14_calc(i)
                  if (j == enefunc%nb14_calc_list(k,i)) then
                    nb15_calc = .false.
                    exit
                  end if
                end do
              end if
              if (.not. nb15_calc) cycle

              num_nb15(id) = num_nb15(id) + 1

              if (.not. do_allocate) &
                pairlist%nb15_calc_list(num_nb15(id),id) = j
            end if

          end do
        end if

        if (.not. do_allocate) then
          pairlist%num_nb15_calc(i,id) = num_nb15(id) - num_nb15_pre(id)
          num_nb15_pre(id)             = num_nb15(id)
        end if

      end do
      !$omp end parallel

      ! allocate memory of pairlist
      !
      if (do_allocate) then
        num_nb15_max = 0
        do i = 1, nthread
          num_nb15_max = max(num_nb15_max, num_nb15(i))
        end do

        num_nb15_max = int(real(num_nb15_max,wp)*FactNumNb15)
        call alloc_pairlist(pairlist, PairListIntNobc, num_nb15_max)
        pairlist%num_nb15_max = num_nb15_max

        do_allocate = .false.
      end if

    end do

    ! check memory size
    !
    call check_pairlist_memory_size(num_nb15, pairlist%num_nb15_max, &
                                    pairlist%allocate_nobc)

    return

  end subroutine update_pairlist_nobc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_nobc_gbsa
  !> @brief        update pairlist under NOBC for GBSA
  !! @authors      TM
  !! @param[in]    enefunc  : information of potential energy function
  !! @param[in]    coord    : atomic coordinate
  !! @param[inout] pairlist : information of pairlist
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_nobc_gbsa(enefunc, coord, pairlist)

    ! formal arguments
    type(s_enefunc),          intent(in)    :: enefunc
    real(wp),                 intent(in)    :: coord(:,:)
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                 :: pairdist2, dij(1:3), rij2
    integer                  :: i, j, k, natom, n, nloops
    integer                  :: num_excl, ini_excl, fin_excl
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: num_nb15_max
    integer                  :: num_all_max
    integer                  :: id, my_id, nthread
#ifdef OMP
    integer                  :: omp_get_thread_num, omp_get_max_threads
#endif
    logical                  :: nb15_calc
    logical                  :: proceed, do_allocate

    integer,     pointer     :: num_nb15_pre(:), num_nb15(:)
    integer,     pointer     :: num_all_pre(:), num_all(:)


    num_nb15_pre => pairlist%num_nb15_pre
    num_nb15     => pairlist%num_nb15
    num_all_pre  => pairlist%num_all_pre
    num_all      => pairlist%num_all
    pairdist2    =  pairlist%pairlistdist * pairlist%pairlistdist
    natom        =  size(coord(1,:))

    ! initialize
    !
#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif
    pairlist%num_nb15_calc(1:natom,1:nthread) = 0
    pairlist%num_all_calc (1:natom,1:nthread) = 0

    ! if pairlist%allocation = .true. , allocation + make pairlist
    ! if pairlist%allocation = .false., make pairlist only
    !
    if (pairlist%allocate_nobc) then
      nloops = 2
      do_allocate = .true.
    else
      nloops = 1
      do_allocate = .false.
    end if

    do n = 1, nloops

      num_excl = 0
      num_nb14 = 0
      num_nb15(1:nthread) = 0
      num_all (1:nthread) = 0

      if (.not. do_allocate) then
        num_nb15_pre(1:nthread) = 0
        num_all_pre (1:nthread) = 0
      end if

      !$omp parallel                                                       &
      !$omp private(id, my_id, i, ini_excl, fin_excl, ini_nb14, fin_nb14,  &
      !$omp         proceed, j, nb15_calc, k, dij, rij2)                   &
      !$omp firstprivate(num_excl, num_nb14, do_allocate)
      !
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      my_id = my_city_rank * nthread + id
      id = id + 1

      do i = 1, natom-1

        ini_excl = num_excl + 1
        fin_excl = num_excl + enefunc%num_nonb_excl(i)
        num_excl = fin_excl

        proceed = .true.
        if (mod(i-1,nproc_city*nthread) /= my_id) proceed = .false.

        if (proceed) then
          do j = i+1, natom

            ! compute distance
            !
            dij(1:3) = coord(1:3,i) - coord(1:3,j)
            rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            ! store interaction table
            !
            if (rij2 < pairdist2) then

              ! for GB
              !
              num_all(id) = num_all(id) + 1

              if (.not. do_allocate) &
                pairlist%all_calc_list(num_all(id),id) = j

              ! 1-2 and 1-3 interactions are removed
              !
              nb15_calc = .true.
              do k = ini_excl, fin_excl
                if (j == enefunc%nonb_excl_list(k)) then
                  nb15_calc = .false.
                  exit
                end if
              end do

              ! 1-4 interactions are removed
              !
              do k = 1, enefunc%num_nb14_calc(i)
                if (j == enefunc%nb14_calc_list(k,i)) then
                  nb15_calc = .false.
                  exit
                end if
              end do
              if (.not. nb15_calc) cycle

              num_nb15(id) = num_nb15(id) + 1

              if (.not. do_allocate) &
                pairlist%nb15_calc_list(num_nb15(id),id) = j
            end if

          end do
        end if

        if (.not. do_allocate) then
          pairlist%num_nb15_calc(i,id) = num_nb15(id) - num_nb15_pre(id)
          num_nb15_pre(id)             = num_nb15(id)

          pairlist%num_all_calc(i,id) = num_all(id) - num_all_pre(id)
          num_all_pre(id)             = num_all(id)
        end if

      end do
      !$omp end parallel

      ! allocate memory of pairlist
      !
      if (do_allocate) then

        num_nb15_max = 0
        do i = 1, nthread
          num_nb15_max = max(num_nb15_max, num_nb15(i))
        end do
        num_nb15_max = int(real(num_nb15_max,wp)*FactNumNb15)
        call alloc_pairlist(pairlist, PairListIntNobc, num_nb15_max)
        pairlist%num_nb15_max = num_nb15_max

        num_all_max = 0
        do i = 1, nthread
          num_all_max = max(num_all_max, num_all(i))
        end do
        num_all_max = int(real(num_all_max,wp)*FactNumNb15)
        call alloc_pairlist(pairlist, PairListIntNobcGbsa, num_all_max)
        pairlist%num_all_max = num_all_max

        do_allocate = .false.
      end if

    end do

    ! check memory size
    !
    call check_pairlist_memory_size(num_nb15, pairlist%num_nb15_max, &
                                    pairlist%allocate_nobc)

    call check_pairlist_memory_size(num_all, pairlist%num_all_max, &
                                    pairlist%allocate_nobc)

    return

  end subroutine update_pairlist_nobc_gbsa

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_nobc_cg
  !> @brief        update pairlist when nobc and AICG2P+3SPN2C 
  !! @authors      CT
  !! @param[in]    enefunc  : information of potential energy function
  !! @param[in]    coord    : atomic coordinate
  !! @param[inout] pairlist : information of pairlist
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_nobc_cg(enefunc, coord, pairlist)
  
    ! formal arguments
    type(s_enefunc),          intent(in)    :: enefunc
    real(wp),                 intent(in)    :: coord(:,:)
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                 :: pairdist2_ele
    real(wp)                 :: pairdist2_126
    real(wp)                 :: pairdist2_PWMcos
    real(wp)                 :: pairdist2_DNAbp
    real(wp)                 :: pairdist2_exv
    real(wp)                 :: pairdist2_longest
    real(wp)                 :: pairdist_longest
    real(wp)                 :: dij(1:3), rij2
    integer                  :: i, j, k, natom, n, nloops
    integer                  :: num_excl, ini_excl, fin_excl
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: num_DNA_bp_max
    integer                  :: num_DNA_exv_max
    integer                  :: num_ele_max
    integer                  :: num_exv_max
    integer                  :: num_kh_max
    integer                  :: i_base_type, j_base_type
    integer                  :: i_chain_id, j_chain_id
    integer                  :: id, my_id, nthread
#ifdef OMP
    integer                  :: omp_get_thread_num, omp_get_max_threads
#endif
    logical                  :: nb15_calc
    logical                  :: proceed, do_allocate
    logical                  :: is_nonlocal, is_WC_bp
    logical                  :: i_is_idr, j_is_idr
    logical                  :: i_is_DNA, j_is_DNA
    logical                  :: ij_is_KH_pair
    logical                  :: need_reallocate_1
    logical                  :: need_reallocate_2
    logical                  :: need_reallocate_3
    logical                  :: need_reallocate_4
    logical                  :: need_reallocate_5

    integer,     pointer     :: num_DNA_bp_pre(:), num_DNA_bp(:)
    integer,     pointer     :: num_DNA_exv_pre(:), num_DNA_exv(:)
    integer,     pointer     :: num_ele_pre(:), num_ele(:)
    integer,     pointer     :: num_exv_pre(:), num_exv(:)
    integer,     pointer     :: num_kh_pre(:),  num_kh(:)


    num_DNA_bp_pre  => pairlist%num_cg_DNA_basepair_pre
    num_DNA_bp      => pairlist%num_cg_DNA_basepair
    num_DNA_exv_pre => pairlist%num_cg_DNA_exv_pre
    num_DNA_exv     => pairlist%num_cg_DNA_exv
    num_ele_pre     => pairlist%num_cg_ele_pre
    num_ele         => pairlist%num_cg_ele
    num_exv_pre     => pairlist%num_cg_exv_pre
    num_exv         => pairlist%num_cg_exv
    num_kh_pre      => pairlist%num_cg_kh_pre
    num_kh          => pairlist%num_cg_kh

    pairdist2_ele    = pairlist%cg_pairlistdist_ele    * pairlist%cg_pairlistdist_ele
    pairdist2_126    = pairlist%cg_pairlistdist_126    * pairlist%cg_pairlistdist_126
    pairdist2_PWMcos = pairlist%cg_pairlistdist_PWMcos * pairlist%cg_pairlistdist_PWMcos
    pairdist2_DNAbp  = pairlist%cg_pairlistdist_DNAbp  * pairlist%cg_pairlistdist_DNAbp
    pairdist2_exv    = pairlist%cg_pairlistdist_exv    * pairlist%cg_pairlistdist_exv

    pairdist_longest = max(pairlist%cg_pairlistdist_ele,    &
                           pairlist%cg_pairlistdist_126,    &
                           pairlist%cg_pairlistdist_PWMcos, &
                           pairlist%cg_pairlistdist_DNAbp,  &
                           pairlist%cg_pairlistdist_exv)
    pairdist2_longest = pairdist_longest * pairdist_longest

    natom                =  size(coord(1,:))

    ! initialize
    !
#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif
    pairlist%num_cg_DNA_basepair_calc (1:natom,1:nthread) = 0
    pairlist%num_cg_DNA_exv_calc      (1:natom,1:nthread) = 0
    pairlist%num_cg_ele_calc          (1:natom,1:nthread) = 0
    pairlist%num_cg_exv_calc          (1:natom,1:nthread) = 0
    pairlist%num_cg_kh_calc           (1:natom,1:nthread) = 0

    if (pairlist%allocate_nobc_cg) then
      nloops      = 2
      do_allocate = .true.
    else
      nloops      = 1
      do_allocate = .false.
    end if
    
    do n = 1, nloops

      num_excl                     = 0
      num_nb14                     = 0
      num_DNA_bp(1:nthread)        = 0
      num_DNA_exv(1:nthread)       = 0
      num_ele(1:nthread)           = 0
      num_exv(1:nthread)           = 0
      num_kh(1:nthread)            = 0

      if (.not. do_allocate) then
        num_DNA_bp_pre(1:nthread)  = 0
        num_DNA_exv_pre(1:nthread) = 0
        num_ele_pre(1:nthread)     = 0
        num_exv_pre(1:nthread)     = 0
        num_kh_pre(1:nthread)      = 0
      end if

      !$omp parallel                                     &
      !$omp private(id, my_id, i, ini_excl, fin_excl,    &
      !$omp         ini_nb14, fin_nb14,                  &
      !$omp         proceed, j, nb15_calc, k, dij, rij2, &
      !$omp         is_nonlocal, is_WC_bp,               &
      !$omp         i_chain_id, j_chain_id,              &
      !$omp         i_is_idr, j_is_idr,                  &
      !$omp         i_is_DNA, j_is_DNA,                  &
      !$omp         ij_is_KH_pair,                       &
      !$omp         i_base_type, j_base_type)            &
      !$omp firstprivate(num_excl, num_nb14, do_allocate)
      !
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      my_id = my_city_rank * nthread + id
      id = id + 1

      do i = 1, natom-1

        ini_excl = num_excl + 1
        fin_excl = num_excl + enefunc%num_nonb_excl(i)
        num_excl = fin_excl

        proceed = .true.
        if (mod(i / 3 ,nproc_city*nthread) /= my_id) proceed = .false.

        if (proceed) then

          i_base_type = enefunc%NA_base_type(i)
          i_chain_id  = enefunc%mol_chain_id(i)
          i_is_idr    = enefunc%cg_IDR_HPS_is_IDR(i) .or. enefunc%cg_IDR_KH_is_IDR(i)
          i_is_DNA    = i_base_type <= NABaseTypeDBMAX .or. i_base_type == NABaseTypeDP .or. i_base_type == NABaseTypeDS

          do j = i + 1, natom

            dij(1:3) = coord(1:3,i) - coord(1:3,j)
            rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            if (rij2 > pairdist2_longest) then
              cycle
            end if

            ! ----------------------------------
            ! remove particles in exclusion list
            ! ----------------------------------
            ! 
            nb15_calc = .true.
            !
            do k = ini_excl, fin_excl
              if (j == enefunc%nonb_excl_list(k)) then
                nb15_calc = .false.
                exit
              end if
            end do
            if (.not. nb15_calc) cycle

            j_base_type = enefunc%NA_base_type(j)
            j_chain_id  = enefunc%mol_chain_id(j)
            j_is_idr    = enefunc%cg_IDR_HPS_is_IDR(j) .or. enefunc%cg_IDR_KH_is_IDR(j)
            j_is_DNA    = j_base_type <= NABaseTypeDBMAX .or. j_base_type == NABaseTypeDP .or. j_base_type == NABaseTypeDS

            ! ----------------------------
            ! test if i and j are KH-pairs
            ! ----------------------------
            ! 
            ij_is_KH_pair = .false.
            if (enefunc%cg_kh_mol_pair(i_chain_id, j_chain_id) > 0 .and. &
                enefunc%cg_pro_use_KH(i) .and. enefunc%cg_pro_use_KH(j)) then
              if (i_is_idr .and. j_is_idr) then
                ij_is_KH_pair = .false.
              else if (j_chain_id /= i_chain_id) then
                ij_is_KH_pair = .true.
              else if (i_is_idr .or. j_is_idr) then
                ij_is_KH_pair = .true.
              end if
            end if
            
            ! =============
            ! exv (general)
            ! =============
            ! 
            if (rij2 < pairdist2_exv) then
              if (i_is_DNA .and. j_is_DNA) then
                ! go to DNA-DNA and do nothing here
              else if (i_is_idr .and. j_is_idr) then
                ! go to idr-idr and do nothing here
              else if (ij_is_KH_pair) then
                ! go to KH-KH and do nothing here
              else
                num_exv(id) = num_exv(id) + 1
                if (.not. do_allocate) then
                  pairlist%cg_exv_list(num_exv(id), id) = j
                end if
              end if
            end if

            ! ================
            ! DNA (exv and bp)
            ! ================
            ! 
            if (i_is_DNA .and. j_is_DNA) then
              if (rij2 < pairdist2_DNAbp) then
                is_nonlocal = (j_chain_id /= i_chain_id) .or. j - i > 4
                is_WC_bp    = .false.
                if (i_base_type <= NABaseTypeDBMAX .and. &
                    j_base_type <= NABaseTypeDBMAX) then
                  if (i_chain_id /= j_chain_id .or. j - i > 10) then 
                    is_WC_bp = enefunc%base_pair_is_WC(i_base_type, j_base_type)
                  end if
                end if
                ! 
                if (is_WC_bp) then
                  num_DNA_bp(id) = num_DNA_bp(id) + 1
                  if (.not. do_allocate) &
                      pairlist%cg_DNA_basepair_list(num_DNA_bp(id), id) = j
                else if (is_nonlocal) then
                  if (rij2 < pairdist2_exv) then
                    num_DNA_exv(id) = num_DNA_exv(id) + 1
                    if (.not. do_allocate) &
                        pairlist%cg_DNA_exv_list(num_DNA_exv(id), id) = j
                  end if
                end if
              end if
            end if

            ! ======================
            ! 12-6 type interactions
            ! ======================
            if (rij2 < pairdist2_126) then
              !
              ! ---------------------
              ! KH model LJ potential
              ! ---------------------
              ! 
              if (ij_is_KH_pair) then
                num_kh(id) = num_kh(id) + 1
                if (.not. do_allocate) then
                  pairlist%cg_KH_list(num_kh(id), id) = j
                  pairlist%cg_KH_model_list(num_kh(id), id) = enefunc%cg_kh_mol_pair(i_chain_id, j_chain_id)
                end if
              end if
              ! 
            end if

            ! ===
            ! ele
            ! ===
            if (rij2 < pairdist2_ele) then
              ! 
              if (enefunc%cg_ele_mol_pair(i_chain_id, j_chain_id) == 1 .and. &
                  abs (enefunc%cg_charge(i)) > EPS .and. &
                  abs (enefunc%cg_charge(j)) > EPS) then
                ! 
                num_ele(id) = num_ele(id) + 1
                if (.not. do_allocate) then
                  pairlist%cg_ele_list(num_ele(id), id) = j
                  if ((i_base_type == NABaseTypeDP .and. &
                      j_base_type > NABaseTypeNAMAX) .or. &
                      (j_base_type == NABaseTypeDP .and. &
                      i_base_type > NABaseTypeNAMAX)) then
                    pairlist%cg_ele_scaling_list(num_ele(id), id) = &
                        - enefunc%cg_pro_DNA_ele_scale_Q / 0.6
                  else
                    pairlist%cg_ele_scaling_list(num_ele(id), id) = 1.0_wp
                  end if
                end if
                ! 
              end if
              ! 
            end if

          end do                ! j = i + 1, n_atom
        end if

        if (.not. do_allocate) then
          pairlist%num_cg_DNA_basepair_calc(i,id) = num_DNA_bp(id)  - num_DNA_bp_pre(id)
          num_DNA_bp_pre(id)                      = num_DNA_bp(id)
          pairlist%num_cg_DNA_exv_calc(i,id)      = num_DNA_exv(id) - num_DNA_exv_pre(id)
          num_DNA_exv_pre(id)                     = num_DNA_exv(id)
          pairlist%num_cg_ele_calc(i,id)          = num_ele(id)     - num_ele_pre(id)
          num_ele_pre(id)                         = num_ele(id)
          pairlist%num_cg_exv_calc(i,id)          = num_exv(id)     - num_exv_pre(id)
          num_exv_pre(id)                         = num_exv(id)
          pairlist%num_cg_kh_calc(i,id)           = num_kh(id)      - num_kh_pre(id)
          num_kh_pre(id)                          = num_kh(id)
        end if

      end do                    ! i = 1, n_atom - 1
      !$omp end parallel

      ! allocate memory of pairlist
      !
      if (do_allocate) then
        num_DNA_bp_max  = 0
        num_DNA_exv_max = 0
        num_ele_max     = 0
        num_exv_max     = 0
        num_kh_max      = 0
        do i = 1, nthread
          num_DNA_bp_max  = max(num_DNA_bp_max, num_DNA_bp(i))
          num_DNA_exv_max = max(num_DNA_exv_max, num_DNA_exv(i))
          num_ele_max     = max(num_ele_max, num_ele(i))
          num_exv_max     = max(num_exv_max, num_exv(i))
          num_kh_max      = max(num_kh_max, num_kh(i))
        end do

        num_DNA_bp_max  = int(real(num_DNA_bp_max,wp )*FactNumNb15)
        num_DNA_exv_max = int(real(num_DNA_exv_max,wp)*FactNumNb15)
        num_ele_max     = int(real(num_ele_max,wp    )*FactNumNb15)
        num_exv_max     = int(real(num_exv_max,wp    )*FactNumNb15)
        num_kh_max      = int(real(num_kh_max,wp     )*FactNumNb15)

        call alloc_pairlist (pairlist, PairListCGDNABP,  num_DNA_bp_max )
        call alloc_pairlist (pairlist, PairListCGDNAexv, num_DNA_exv_max)
        call alloc_pairlist (pairlist, PairListCGele,    num_ele_max    )
        call alloc_pairlist (pairlist, PairListCGexv,    num_exv_max    )
        call alloc_pairlist (pairlist, PairListCGKH,     num_kh_max     )

        pairlist%num_cg_DNA_basepair_max = num_DNA_bp_max
        pairlist%num_cg_DNA_exv_max      = num_DNA_exv_max
        pairlist%num_cg_ele_max          = num_ele_max
        pairlist%num_cg_exv_max          = num_exv_max
        pairlist%num_cg_kh_max           = num_kh_max

        do_allocate = .false.
      end if

    end do

    ! check memory size
    !
    call check_pairlist_memory_size(num_DNA_bp, &
         pairlist%num_cg_DNA_basepair_max, need_reallocate_1)
    call check_pairlist_memory_size(num_DNA_exv, &
         pairlist%num_cg_DNA_exv_max, need_reallocate_2)
    call check_pairlist_memory_size(num_ele, &
         pairlist%num_cg_ele_max, need_reallocate_3)
    call check_pairlist_memory_size(num_exv, &
         pairlist%num_cg_exv_max, need_reallocate_4)
    call check_pairlist_memory_size(num_kh, &
         pairlist%num_cg_kh_max, need_reallocate_5)
    pairlist%allocate_nobc_cg = need_reallocate_1 .or. need_reallocate_2 .or. &
        need_reallocate_3 .or. need_reallocate_4 .or. need_reallocate_5

    return

  end subroutine update_pairlist_nobc_cg

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_nobc_cg_pwmcos
  !> @brief        update pairlist for CG protein-DNA seq-(non)specific
  !!               interaction
  !! @authors      CT
  !! @param[in]    enefunc  : information of potential energy function
  !! @param[in]    coord    : atomic coordinate
  !! @param[inout] pairlist : information of pairlist
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_nobc_cg_pwmcos(enefunc, coord, pairlist)
  
    ! formal arguments
    type(s_enefunc),          intent(in)    :: enefunc
    real(wp),                 intent(in)    :: coord(:,:)
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                 :: pairdist2_PWMcos
    real(wp)                 :: rij2
    real(wp)                 :: dx, dy, dz
    real(wp)                 :: x1, y1, z1
    real(wp)                 :: x2, y2, z2
    integer                  :: i, j, k, n, nloops
    integer                  :: natom, npwmcos
    integer                  :: num_pwmcos_max
    integer                  :: i_pro_pwmcos
    integer                  :: specificity
    integer                  :: i_chain_id, j_chain_id
    integer                  :: i_base_type, j_base_type
    integer                  :: id, my_id, nthread
    integer                  :: alloc_stat
#ifdef OMP
    integer                  :: omp_get_thread_num, omp_get_max_threads
#endif
    logical                  :: proceed, do_allocate

    integer,     pointer     :: num_pwmcos(:)


    pairdist2_PWMcos = pairlist%cg_pairlistdist_PWMcos  &
        * pairlist%cg_pairlistdist_PWMcos

    num_pwmcos => pairlist%num_cg_pwmcos

    natom      =  size(coord(1,:))
    npwmcos    =  enefunc%num_pwmcos_resid

    ! initialize
    !
#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif

    if (pairlist%allocate_nobc_cg_pwmcos) then
      nloops      = 2
      do_allocate = .true.
    else
      nloops      = 1
      do_allocate = .false.
    end if
   
    do n = 1, nloops

      num_pwmcos(1:npwmcos)       = 0

      !$omp parallel                          &
      !$omp private(id, my_id, i, j, k,       &
      !$omp         proceed,                  &
      !$omp         i_pro_pwmcos,             &
      !$omp         x1, y1, z1, x2, y2, z2,   &
      !$omp         dx, dy, dz, rij2,         &
      !$omp         specificity,              &
      !$omp         alloc_stat,               &
      !$omp         i_chain_id, j_chain_id,   &
      !$omp         i_base_type, j_base_type) &
      !$omp firstprivate(do_allocate)
      !
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      my_id = my_city_rank * nthread + id
      id = id + 1
      
      do i = 1, npwmcos

        proceed = .true.
        if (mod(i - 1, nproc_city*nthread) /= my_id) proceed = .false.

        if (proceed) then

          i_pro_pwmcos = enefunc%pwmcos_involved_resid(i)
          specificity  = enefunc%pwmcos_involved_spec (i)

          i_base_type  = enefunc%NA_base_type(i_pro_pwmcos)
          i_chain_id   = enefunc%mol_chain_id(i_pro_pwmcos)

          x1 = coord(1, i_pro_pwmcos)
          y1 = coord(2, i_pro_pwmcos)
          z1 = coord(3, i_pro_pwmcos)

          do j = 1, natom

            j_base_type = enefunc%NA_base_type(j)

            if (j_base_type > NABaseTypeDBMAX) &
                cycle

            j_chain_id  = enefunc%mol_chain_id(j)

            if (enefunc%pwmcos_mol_pair(i_chain_id, j_chain_id) == 0) &
                cycle

            x2 = coord(1, j)
            y2 = coord(2, j)
            z2 = coord(3, j)

            dx = x2 - x1
            dy = y2 - y1
            dz = z2 - z1
            
            rij2 = dx * dx + dy * dy + dz * dz

            ! store interaction table
            !
            if (rij2 < pairdist2_PWMcos) then

              if (specificity == 1) then
                if (((j - 3 >= 1) .and. &
                    (enefunc%mol_chain_id(j - 3) == j_chain_id)) .and. &
                    ((j + 3 <= natom) .and. &
                    (enefunc%mol_chain_id(j + 3) == j_chain_id)) &
                   ) then
                  num_pwmcos(i) = num_pwmcos(i) + 1
                  if (.not. do_allocate) &
                      pairlist%cg_pwmcos_list(num_pwmcos(i), i) = j
                end if
              end if

            end if

          end do

        end if

      end do
      !$omp end parallel

      ! allocate memory of pairlist
      !
      if (do_allocate) then
        num_pwmcos_max  = 0
        do i = 1, npwmcos
          num_pwmcos_max  = max(num_pwmcos_max, num_pwmcos(i))
        end do

        num_pwmcos_max  = int(real(num_pwmcos_max, wp)*FactNumNb15)

        ! if (allocated(pairlist%cg_pwmcos_list)) &
        !     deallocate(pairlist%cg_pwmcos_list)
        ! allocate(pairlist%cg_pwmcos_list(num_pwmcos_max, npwmcos), stat=alloc_stat)
        ! if (alloc_stat /= 0) call error_msg_alloc
        call alloc_pairlist (pairlist, PairListCGPWMcos, num_pwmcos_max, npwmcos)

        pairlist%num_cg_pwmcos_max = num_pwmcos_max

        do_allocate = .false.
      end if

    end do

    ! check memory size
    !
    call check_pairlist_memory_size(num_pwmcos, pairlist%num_cg_pwmcos_max, &
        pairlist%allocate_nobc_cg_pwmcos)

    return

  end subroutine update_pairlist_nobc_cg_pwmcos

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_nobc_cg_pwmcosns
  !> @brief        update pairlist for CG protein-DNA seq-(non)specific
  !!               interaction
  !! @authors      CT
  !! @param[in]    enefunc  : information of potential energy function
  !! @param[in]    coord    : atomic coordinate
  !! @param[inout] pairlist : information of pairlist
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_nobc_cg_pwmcosns(enefunc, coord, pairlist)
  
    ! formal arguments
    type(s_enefunc),          intent(in)    :: enefunc
    real(wp),                 intent(in)    :: coord(:,:)
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                 :: pairdist2_PWMcos
    real(wp)                 :: rij2
    real(wp)                 :: dx, dy, dz
    real(wp)                 :: x1, y1, z1
    real(wp)                 :: x2, y2, z2
    integer                  :: i, j, k, n, nloops
    integer                  :: natom, npwmcosns
    integer                  :: num_pwmcosns_max
    integer                  :: i_pro_pwmcosns
    integer                  :: specificity
    integer                  :: i_chain_id, j_chain_id
    integer                  :: i_base_type, j_base_type
    integer                  :: id, my_id, nthread
    integer                  :: alloc_stat
#ifdef OMP
    integer                  :: omp_get_thread_num, omp_get_max_threads
#endif
    logical                  :: proceed, do_allocate

    integer,     pointer     :: num_pwmcosns(:)


    pairdist2_PWMcos = pairlist%cg_pairlistdist_PWMcos  &
        * pairlist%cg_pairlistdist_PWMcos

    num_pwmcosns => pairlist%num_cg_pwmcosns

    natom     = size(coord(1,:))
    npwmcosns = enefunc%num_pwmcosns_resid

    ! initialize
    !
#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif

    if (pairlist%allocate_nobc_cg_pwmcosns) then
      nloops      = 2
      do_allocate = .true.
    else
      nloops      = 1
      do_allocate = .false.
    end if
   
    do n = 1, nloops

      num_pwmcosns(1:npwmcosns)       = 0

      !$omp parallel                          &
      !$omp private(id, my_id, i, j, k,       &
      !$omp         proceed,                  &
      !$omp         i_pro_pwmcosns,           &
      !$omp         x1, y1, z1, x2, y2, z2,   &
      !$omp         dx, dy, dz, rij2,         &
      !$omp         specificity,              &
      !$omp         alloc_stat,               &
      !$omp         i_chain_id, j_chain_id,   &
      !$omp         i_base_type, j_base_type) &
      !$omp firstprivate(do_allocate)
      !
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      my_id = my_city_rank * nthread + id
      id = id + 1
      
      do i = 1, npwmcosns

        proceed = .true.
        if (mod(i - 1, nproc_city*nthread) /= my_id) proceed = .false.

        if (proceed) then

          i_pro_pwmcosns = enefunc%pwmcosns_involved_resid(i)
          specificity  = enefunc%pwmcosns_involved_spec (i)

          i_base_type  = enefunc%NA_base_type(i_pro_pwmcosns)
          i_chain_id   = enefunc%mol_chain_id(i_pro_pwmcosns)

          x1 = coord(1, i_pro_pwmcosns)
          y1 = coord(2, i_pro_pwmcosns)
          z1 = coord(3, i_pro_pwmcosns)

          do j = 1, natom

            j_base_type = enefunc%NA_base_type(j)

            if (j_base_type /= NABaseTypeDP) &
                cycle

            j_chain_id  = enefunc%mol_chain_id(j)

            if (enefunc%pwmcosns_mol_pair(i_chain_id, j_chain_id) == 0) &
                cycle

            x2 = coord(1, j)
            y2 = coord(2, j)
            z2 = coord(3, j)

            dx = x2 - x1
            dy = y2 - y1
            dz = z2 - z1
            
            rij2 = dx * dx + dy * dy + dz * dz

            ! store interaction table
            !
            if (rij2 < pairdist2_PWMcos) then
              num_pwmcosns(i) = num_pwmcosns(i) + 1
              if (.not. do_allocate) &
                  pairlist%cg_pwmcosns_list(num_pwmcosns(i), i) = j
            end if

          end do

        end if

      end do
      !$omp end parallel

      ! allocate memory of pairlist
      !
      if (do_allocate) then
        num_pwmcosns_max  = 0
        do i = 1, npwmcosns
          num_pwmcosns_max  = max(num_pwmcosns_max, num_pwmcosns(i))
        end do

        num_pwmcosns_max  = int(real(num_pwmcosns_max, wp)*FactNumNb15)

        call alloc_pairlist (pairlist, PairListCGPWMcosns, num_pwmcosns_max, npwmcosns)

        pairlist%num_cg_pwmcosns_max = num_pwmcosns_max

        do_allocate = .false.
      end if

    end do

    ! check memory size
    !
    call check_pairlist_memory_size(num_pwmcosns, pairlist%num_cg_pwmcosns_max, &
        pairlist%allocate_nobc_cg_pwmcosns)

    return

  end subroutine update_pairlist_nobc_cg_pwmcosns

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_nobc_cg_IDR_HPS
  !> @brief        update pairlist for nobc in IDR HPS model
  !! @authors      CT
  !! @param[in]    enefunc  : information of potential energy function
  !! @param[in]    coord    : atomic coordinate
  !! @param[inout] pairlist : information of pairlist
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_nobc_cg_IDR_HPS(enefunc, coord, pairlist)
  
    ! formal arguments
    type(s_enefunc),          intent(in)    :: enefunc
    real(wp),                 intent(in)    :: coord(:,:)
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                 :: pairdist2_126
    real(wp)                 :: dij(1:3), rij2
    integer                  :: i, j, natom, n, nloops
    integer                  :: i_chain_id, j_chain_id
    integer                  :: num_IDR_HPS_max
    integer                  :: id, my_id, nthread
#ifdef OMP
    integer                  :: omp_get_thread_num, omp_get_max_threads
#endif
    logical                  :: proceed, do_allocate

    integer,     pointer     :: num_IDR_HPS_pre(:), num_IDR_HPS(:)


    num_IDR_HPS_pre => pairlist%num_cg_IDR_HPS_pre
    num_IDR_HPS     => pairlist%num_cg_IDR_HPS

    pairdist2_126 = pairlist%cg_pairlistdist_126  &
        * pairlist%cg_pairlistdist_126
    natom         =  size(coord(1,:))

    ! initialize
    !
#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif
    pairlist%num_cg_IDR_HPS_calc(1:natom,1:nthread) = 0

    if (pairlist%allocate_nobc_cg_IDR_HPS) then
      nloops      = 2
      do_allocate = .true.
    else
      nloops      = 1
      do_allocate = .false.
    end if
    
    do n = 1, nloops

      num_IDR_HPS(1:nthread)       = 0

      if (.not. do_allocate) then
        num_IDR_HPS_pre(1:nthread) = 0
      end if

      !$omp parallel                    &
      !$omp private(id, my_id, i, j,    &
      !$omp         i_chain_id,         &
      !$omp         j_chain_id,         &
      !$omp         proceed, dij, rij2) &
      !$omp firstprivate(do_allocate)
      !
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      my_id = my_city_rank * nthread + id
      id = id + 1

      do i = 1, natom-1

        proceed = .true.
        if (mod(i, nproc_city*nthread) /= my_id) proceed = .false.
        if (.not. enefunc%cg_IDR_HPS_is_IDR(i))  proceed = .false.

        if (proceed) then

          i_chain_id  = enefunc%mol_chain_id(i)

          do j = i + 1, natom

            if (.not. enefunc%cg_IDR_HPS_is_IDR(j)) cycle

            j_chain_id  = enefunc%mol_chain_id(j)
            ! 
            if (i_chain_id == j_chain_id .and. i == j - 1) then
              cycle
            end if

            dij(1:3) = coord(1:3,i) - coord(1:3,j)
            rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            ! store interaction table
            !
            if (rij2 < pairdist2_126) then
              num_IDR_HPS(id) = num_IDR_HPS(id) + 1
              if (.not. do_allocate) then
                pairlist%cg_IDR_HPS_list(num_IDR_HPS(id), id) = j
              end if
            end if

          end do

        end if

        if (.not. do_allocate) then
          pairlist%num_cg_IDR_HPS_calc(i,id) = num_IDR_HPS(id) - num_IDR_HPS_pre(id)
          num_IDR_HPS_pre(id)                = num_IDR_HPS(id)
        end if

      end do
      !$omp end parallel

      ! allocate memory of pairlist
      !
      if (do_allocate) then
        num_IDR_HPS_max     = 0
        do i = 1, nthread
          num_IDR_HPS_max = max(num_IDR_HPS_max, num_IDR_HPS(i))
        end do

        num_IDR_HPS_max = int(real(num_IDR_HPS_max,wp)*FactNumNb15)

        call alloc_pairlist (pairlist, PairListCGIDRHPS, num_IDR_HPS_max)

        pairlist%num_cg_IDR_HPS_max      = num_IDR_HPS_max

        do_allocate = .false.
      end if

    end do

    ! check memory size
    !
    call check_pairlist_memory_size(num_IDR_HPS, pairlist%num_cg_IDR_HPS_max, &
        pairlist%allocate_nobc_cg_IDR_HPS)

    return

  end subroutine update_pairlist_nobc_cg_IDR_HPS

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_nobc_cg_IDR_KH
  !> @brief        update pairlist for nobc in IDR KH model
  !! @authors      CT
  !! @param[in]    enefunc  : information of potential energy function
  !! @param[in]    coord    : atomic coordinate
  !! @param[inout] pairlist : information of pairlist
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_nobc_cg_IDR_KH(enefunc, coord, pairlist)
  
    ! formal arguments
    type(s_enefunc),          intent(in)    :: enefunc
    real(wp),                 intent(in)    :: coord(:,:)
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                 :: pairdist2_126
    real(wp)                 :: dij(1:3), rij2
    integer                  :: i, j, natom, n, nloops
    integer                  :: i_chain_id, j_chain_id
    integer                  :: num_IDR_KH_max
    integer                  :: id, my_id, nthread
#ifdef OMP
    integer                  :: omp_get_thread_num, omp_get_max_threads
#endif
    logical                  :: proceed, do_allocate

    integer,     pointer     :: num_IDR_KH_pre(:), num_IDR_KH(:)


    num_IDR_KH_pre => pairlist%num_cg_IDR_KH_pre
    num_IDR_KH     => pairlist%num_cg_IDR_KH

    pairdist2_126 = pairlist%cg_pairlistdist_126  &
        * pairlist%cg_pairlistdist_126
    natom         =  size(coord(1,:))

    ! initialize
    !
#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif
    pairlist%num_cg_IDR_KH_calc(1:natom,1:nthread) = 0

    if (pairlist%allocate_nobc_cg_IDR_KH) then
      nloops      = 2
      do_allocate = .true.
    else
      nloops      = 1
      do_allocate = .false.
    end if
    
    do n = 1, nloops

      num_IDR_KH(1:nthread)       = 0

      if (.not. do_allocate) then
        num_IDR_KH_pre(1:nthread) = 0
      end if

      !$omp parallel                    &
      !$omp private(id, my_id, i, j,    &
      !$omp         i_chain_id,         &
      !$omp         j_chain_id,         &
      !$omp         proceed, dij, rij2) &
      !$omp firstprivate(do_allocate)
      !
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      my_id = my_city_rank * nthread + id
      id = id + 1

      do i = 1, natom - 1

        proceed = .true.
        if (mod(i, nproc_city*nthread) /= my_id) proceed = .false.
        if (.not. enefunc%cg_IDR_KH_is_IDR(i))  proceed = .false.

        if (proceed) then

          i_chain_id  = enefunc%mol_chain_id(i)

          do j = i + 1, natom

            if (.not. enefunc%cg_IDR_KH_is_IDR(j)) cycle
            
            j_chain_id  = enefunc%mol_chain_id(j)
            ! 
            if (i_chain_id == j_chain_id .and. i == j - 1) then
              cycle
            end if

            dij(1:3) = coord(1:3,i) - coord(1:3,j)
            rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            ! store interaction table
            !
            if (rij2 < pairdist2_126) then
              num_IDR_KH(id) = num_IDR_KH(id) + 1
              if (.not. do_allocate) then
                pairlist%cg_IDR_KH_list(num_IDR_KH(id), id) = j
              end if
            end if

          end do

        end if

        if (.not. do_allocate) then
          pairlist%num_cg_IDR_KH_calc(i,id) = num_IDR_KH(id) - num_IDR_KH_pre(id)
          num_IDR_KH_pre(id)                = num_IDR_KH(id)
        end if

      end do
      !$omp end parallel

      ! allocate memory of pairlist
      !
      if (do_allocate) then
        num_IDR_KH_max     = 0
        do i = 1, nthread
          num_IDR_KH_max = max(num_IDR_KH_max, num_IDR_KH(i))
        end do

        num_IDR_KH_max = int(real(num_IDR_KH_max,wp)*FactNumNb15)

        call alloc_pairlist (pairlist, PairListCGIDRKH, num_IDR_KH_max)

        pairlist%num_cg_IDR_KH_max      = num_IDR_KH_max

        do_allocate = .false.
      end if

    end do

    ! check memory size
    !
    call check_pairlist_memory_size(num_IDR_KH, pairlist%num_cg_IDR_KH_max, &
        pairlist%allocate_nobc_cg_IDR_KH)

    return

  end subroutine update_pairlist_nobc_cg_IDR_KH

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_ecqm15_nonb
  !> @brief        create ecqm_nb15_list
  !! @authors      KY
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    coord    : atomic coordinate
  !! @param[in]    pairlist : information of pairlist
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_ecqm15_nonb(enefunc, coord, pairlist)

    ! formal arguments
    type(s_enefunc),  target, intent(in)   :: enefunc
    real(wp),                 intent(in)   :: coord(:,:)
    type(s_pairlist), target, intent(inout):: pairlist

    ! local variables
    real(wp)                 :: pairdist2, dij(1:3), rij2
    integer                  :: i, j, k, ii, jj
    logical                  :: nb15_calc

    integer,     pointer     :: nonb_excl_list(:,:), nb14_calc_list(:,:)
    integer                  :: qm_natoms, ec_natoms, iqm, iec
    integer,     pointer     :: qmatom_id(:), ecatom_id(:)
    integer, allocatable     :: num_nb15(:), nb15_list(:,:)
    integer                  :: kalloc_stat, kdealloc_stat
    integer                  :: n_init, n_final

    integer                  :: id, nthread
#ifdef OMP
    integer                  :: omp_get_thread_num, omp_get_max_threads
#endif


    !dbg write(MsgOut,'("QM-EC_15 pair")')

    qm_natoms =  enefunc%qmmm%qm_natoms
    qmatom_id => enefunc%qmmm%qmatom_id
    ec_natoms =  enefunc%qmmm%ec_natoms
    ecatom_id => enefunc%qmmm%ecatom_id

    nonb_excl_list => enefunc%table%nonb_excl_list
    nb14_calc_list => enefunc%nb14_calc_list
    pairdist2      =  pairlist%pairlistdist * pairlist%pairlistdist

    if (.not. allocated(pairlist%ecqm_nb15_list)) then
      call alloc_pairlist(pairlist, PairListEcqm, qm_natoms*ec_natoms)
    end if

#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif

    !dbg write(MsgOut,'("qm_natom, ec_natom",2i4)') qm_natoms, ec_natoms

    allocate(num_nb15(nthread), stat=kalloc_stat)
    if (kalloc_stat /=0) call error_msg_alloc
    allocate(nb15_list(2, qm_natoms*ec_natoms), stat=kalloc_stat)
    if (kalloc_stat /=0) call error_msg_alloc

    !$omp parallel                                             &
    !$omp private(id, iqm, iec, i, j, k, dij, rij2, nb15_calc, &
    !$omp         nb15_list, n_init, n_final)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    id = id + 1

    num_nb15(id) = 0
    do iqm = 1, qm_natoms
      if (mod(iqm-1,nthread) /= id-1) cycle
      !dbg write(MsgOut,'("id, nthread, mod(iqm-1,nthread)",3i4)') id, nthread, mod(iqm-1,nthread)

      do iec = 1, ec_natoms
        if (qmatom_id(iqm) < ecatom_id(iec)) then
          i = qmatom_id(iqm)
          j = ecatom_id(iec)
        else
          i = ecatom_id(iec)
          j = qmatom_id(iqm)
        end if

        ! compute distance
        !
        dij(1:3) = coord(1:3,i) - coord(1:3,j)
        rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        ! store interaction table
        !
        if (rij2 < pairdist2) then

          ! 1-2 and 1-3 interactions are removed
          !
          nb15_calc = .true.
          do k = 1, enefunc%num_nonb_excl(i)
            if (j .eq. nonb_excl_list(k,i)) then
              nb15_calc = .false.
              exit
            end if
          end do
          if (.not. nb15_calc) cycle

          ! 1-4 interactions are removed
          !
          if (enefunc%forcefield/=ForcefieldKBGO) then
            do k = 1, enefunc%num_nb14_calc(i)
              if (j == nb14_calc_list(k,i)) then
                nb15_calc = .false.
                exit
              end if
            end do
          end if
          if (.not. nb15_calc) cycle

          ! NOTE: i < j
          num_nb15(id) = num_nb15(id) + 1
          nb15_list(1,num_nb15(id)) = i
          nb15_list(2,num_nb15(id)) = j
        end if

      end do
    end do

    !$omp barrier

    !dbg write(MsgOut,'("id,num_nb15",5i4)') id, num_nb15

    n_init  = 1
    do i = 1, id-1
       n_init = n_init + num_nb15(i)
    end do
    n_final = n_init + num_nb15(id) - 1
    !dbg write(MsgOut,'("id,n_init,n_final",3i4)') id, n_init, n_final

    if (id == nthread)  pairlist%ecqm_num_nb15 = n_final 
    if (num_nb15(id) > 0) then
      pairlist%ecqm_nb15_list(:,n_init:n_final) = nb15_list(:,1:num_nb15(id))
    end if

    !$omp end parallel

    deallocate(nb15_list, num_nb15, stat=kdealloc_stat)
    if (kdealloc_stat /=0) call error_msg_dealloc

    !dbg write(MsgOut,'("ecqm_num_nb15",i4)') pairlist%ecqm_num_nb15
    !dbg do i = 1, pairlist%ecqm_num_nb15
    !dbg   write(MsgOut,'(2i4)') pairlist%ecqm_nb15_list(:,i)
    !dbg end do

    return

  end subroutine update_ecqm15_nonb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_solute_solute
  !> @brief        update pairlist between solutes
  !! @authors      JJ, TM, MK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : boundary conditions information
  !! @param[in]    coord    : coordinates
  !! @param[inout] pairlist : pairlist information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_solute_solute(enefunc, boundary, coord, pairlist)

    ! formal arguments
    type(s_enefunc),          intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord(:,:)
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: dij(1:3), rij2, pairdist2
    real(wp)                  :: dij_pbc(1:3), origin(1:3)
    real(wp)                  :: r_shift(1:3)
    real(wp)                  :: csize_inv(1:3)
    real(wp)                  :: box_size(1:3), half_box_size(1:3)
    integer                   :: ii, jj, i, j, k, ik, jk, n, nloops
    integer                   :: natom
    integer                   :: ic(1:3), icel, inbc
    integer                   :: num_nb15_max
    integer                   :: id, my_id, nthread
    integer                   :: ncell(1:3)
#ifdef OMP
    integer                   :: omp_get_max_threads, omp_get_thread_num
#endif
    logical                   :: nb15_calc
    logical                   :: do_allocate

    integer,          pointer :: num_nb15_pre(:), num_nb15(:), num_nb15_calc(:)
    integer,          pointer :: list(:), ic_atm(:), head(:)
    integer,          pointer :: ncel


    natom   = enefunc%table%num_solute

    num_nb15_pre => pairlist%num_nb15_pre
    num_nb15     => pairlist%num_nb15
    num_nb15_calc => pairlist%table%num_nb15_calc
    list         => pairlist%table%cell_linked_list
    ic_atm       => pairlist%table%atom_cell_index
    head         => boundary%cell_head_atom
    ncel         => boundary%num_cells

    box_size(1)  = boundary%box_size_x
    box_size(2)  = boundary%box_size_y
    box_size(3)  = boundary%box_size_z
    origin(1)    = boundary%origin_x
    origin(2)    = boundary%origin_y
    origin(3)    = boundary%origin_z
    ncell(1)     = boundary%num_cells_x
    ncell(2)     = boundary%num_cells_y
    ncell(3)     = boundary%num_cells_z
    csize_inv(1) = 1.0_wp / boundary%cell_size_x
    csize_inv(2) = 1.0_wp / boundary%cell_size_y
    csize_inv(3) = 1.0_wp / boundary%cell_size_z
    half_box_size(1:3) = box_size(1:3)*0.5_wp

    pairdist2 = pairlist%pairlistdist * pairlist%pairlistdist

    ! initialize
    !
#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif


    ! make cell linked list for solute
    !
    head(1:ncel) = 0
    do ii = 1, natom

      ! coordinate shifted to the first quadrant and set into the boundary box
      !
      i = enefunc%table%solute_list(ii)
      r_shift(1:3) = coord(1:3,i) - origin(1:3)
      r_shift(1:3) = r_shift(1:3) + half_box_size(1:3) -  &
                     box_size(1:3)*anint(r_shift(1:3)/box_size(1:3))

      ! assign which cell
      !
      ic(1:3) = int(r_shift(1:3)*csize_inv(1:3))

      ! upper limit of ic[x,y,z] should be ncel_[x,y,z]
      !
      do k = 1, 3
        if (ic(k) >= ncell(k)) ic(k) = ncell(k)-1
      end do
      icel = 1 + ic(1) + ic(2)*ncell(1) + ic(3)*ncell(1)*ncell(2)

      ! store which cell atom i is placed
      !
      ic_atm(ii) = icel
      list(ii)   = head(icel)
      head(icel) = ii

    end do

    num_nb15_calc(1:natom) = 0

    ! if pairlist%allocation = .true. , allocation + make pairlist
    ! if pairlist%allocation = .false., make pairlist only
    !
    if (pairlist%allocate_solsol) then
      nloops = 2
      do_allocate = .true.
    else
      nloops = 1
      do_allocate = .false.
    end if

    do n = 1, nloops

      !$omp parallel                                      &
      !$omp private(id, my_id, ii, jj, i, j, k, inbc,     &
      !$omp         nb15_calc, dij, dij_pbc, rij2, icel) 
      !
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      my_id = my_city_rank * nthread + id
      id    = id + 1

      do icel = 1, ncel

        ii = head(icel)
        if (ii <= 0) cycle

        do while (ii > 0)

          i = enefunc%table%solute_list(ii)

          if (mod(ii-1,nproc_city*nthread) == my_id) then

            num_nb15_calc(ii) = 0

            ! serch only in the self and neighboring cells 
            !
            do inbc = 1, boundary%num_neighbor_cells
              jj = head(boundary%neighbor_cells(inbc,icel))

              do while (jj > 0)

                if (jj <= ii) then
                  jj = list(jj)
                  cycle
                end if

                ! 1-2 and 1-3 interactions are removed
                !
                nb15_calc = .true.
                do k = 1, enefunc%table%num_nonb_excl(ii)
                  if (jj == enefunc%table%nonb_excl_list(k,ii)) then
                    nb15_calc = .false.
                    exit
                  end if
                end do

                ! 1-4 interactions are removed
                !
                do k = 1, enefunc%num_nb14_calc(ii)
                  if (jj == enefunc%nb14_calc_list(k,ii)) then
                    nb15_calc = .false.
                    exit
                  end if
                end do

                if (.not.nb15_calc) then
                  jj = list(jj)
                  cycle
                end if

                j = enefunc%table%solute_list(jj)

                ! compute distance
                !
                dij(1:3) = coord(1:3,i) - coord(1:3,j)

                ! consider the periodic boundary
                !
                dij_pbc(1:3) = box_size(1:3)*anint(dij(1:3)/box_size(1:3))
                dij(1:3)     = dij(1:3) - dij_pbc(1:3)
                rij2         = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

                ! cout the number of interactions
                !
                if (rij2 < pairdist2) then
                  num_nb15_calc(ii) = num_nb15_calc(ii) + 1

                  if (.not. do_allocate) &
                    pairlist%table%nb15_calc_list(num_nb15_calc(ii),ii) = j
                end if

                jj = list(jj)

              end do
            end do
          end if

          ii = list(ii)

        end do
      end do
      !$omp end parallel

      ! allocate memory of pairlist
      !
      if (do_allocate) then
        num_nb15_max = max(1,maxval(num_nb15_calc(1:natom)))
        num_nb15_max = int(real(num_nb15_max,wp)*FactNumNb15)
        call alloc_pairlist2(pairlist, PairListPbcSoluteSolute,           &
                             num_nb15_max, natom)
        pairlist%num_nb15_max = num_nb15_max

        do_allocate = .false.
      end if

    end do

    ! check memory size
    !
    call check_pairlist_memory_size(num_nb15_calc, pairlist%num_nb15_max, &
                                    pairlist%allocate_solsol)

    return

  end subroutine update_pairlist_solute_solute

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_solute_water
  !> @brief        update pairlist between solute and water
  !! @authors      JJ, TM, MK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : boundary conditions information
  !! @param[in]    coord    : coordinates
  !! @param[inout] pairlist : pairlist information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_solute_water(enefunc, boundary, coord, pairlist)

    ! formal arguments
    type(s_enefunc),          intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord(:,:)
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: dij(1:3), rij2, pairdist2
    real(wp)                  :: dij_pbc(1:3), origin(1:3)
    real(wp)                  :: r_shift(1:3)
    real(wp)                  :: csize_inv(1:3)
    real(wp)                  :: box_size(1:3), half_box_size(1:3)
    integer                   :: ii, jj, i, j, k, ik, jk, n, nloops
    integer                   :: natom, nwater
    integer                   :: num_excl, ini_excl, fin_excl
    integer                   :: num_nb14, ini_nb14, fin_nb14
    integer                   :: ic(1:3), icel, inbc
    integer                   :: num_nb15_max
    integer                   :: id, my_id, nthread
    integer                   :: ncell(1:3)
#ifdef OMP
    integer                   :: omp_get_max_threads, omp_get_thread_num
#endif
    logical                   :: nb15_calc, do_allocate

    integer,          pointer :: num_nb15_pre(:), num_nb15_calcw(:)
    integer,          pointer :: ic_atm(:), head(:), list(:)
    integer,          pointer :: headw(:), listw(:)
    integer,          pointer :: ncel


    natom   = enefunc%table%num_solute
    nwater  = enefunc%table%num_water

    num_nb15_pre => pairlist%num_nb15_pre
    num_nb15_calcw => pairlist%table%num_nb15_calcw
    ic_atm       => pairlist%table%atom_cell_index
    list         => pairlist%table%cell_linked_list
    listw        => pairlist%table%cell_linked_listw
    head         => boundary%cell_head_atom
    headw        => boundary%cell_head_atomw
    ncel         => boundary%num_cells

    box_size(1)  = boundary%box_size_x
    box_size(2)  = boundary%box_size_y
    box_size(3)  = boundary%box_size_z
    origin(1)    = boundary%origin_x
    origin(2)    = boundary%origin_y
    origin(3)    = boundary%origin_z
    ncell(1)     = boundary%num_cells_x
    ncell(2)     = boundary%num_cells_y
    ncell(3)     = boundary%num_cells_z
    csize_inv(1) = 1.0_wp / boundary%cell_size_x
    csize_inv(2) = 1.0_wp / boundary%cell_size_y
    csize_inv(3) = 1.0_wp / boundary%cell_size_z
    half_box_size(1:3) = box_size(1:3)*0.5_wp

    pairdist2 = pairlist%pairlistdist * pairlist%pairlistdist

    ! initialize
    !
#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif
    num_nb15_calcw(1:natom) = 0


    ! make cell linked list of water oxygen
    !
    headw(1:ncel) = 0
    ik = 0
    do ii = 1, nwater

      ! water oxygen index
      !
      do k = 1, 3

        ik = ik + 1

        ! coordinate shifted to the first quadrant and set into the boundary box
        !
        i = enefunc%table%water_list(k,ii)
        r_shift(1:3) = coord(1:3,i) - origin(1:3)
        r_shift(1:3) = r_shift(1:3) + half_box_size(1:3) - &
                       box_size(1:3)*anint(r_shift(1:3)/box_size(1:3))

        ! assign which cell
        !
        ic(1:3) = int(r_shift(1:3)*csize_inv(1:3))

        ! upper limit of ic[x,y,z] should be ncel_[x,y,z]
        !
        if (ic(1) >= ncell(1)) ic(1) = ncell(1)-1
        if (ic(2) >= ncell(2)) ic(2) = ncell(2)-1
        if (ic(3) >= ncell(3)) ic(3) = ncell(3)-1
        icel = 1 + ic(1) + ic(2)*ncell(1) + ic(3)*ncell(1)*ncell(2)

        ! store which cell atom i is placed
        !
        listw(ik)   = headw(icel)
        headw(icel)  = ik
      end do

    end do

    num_nb15_calcw(1:natom) = 0

    ! if pairlist%allocation = .true. , allocation + make pairlist
    ! if pairlist%allocation = .false., make pairlist only
    !
    if (pairlist%allocate_solwat) then
      nloops = 2
      do_allocate = .true.
    else
      nloops = 1
      do_allocate = .false.
    end if

    do n = 1, nloops

      if (.not. do_allocate) &
        num_nb15_calcw(1:natom) = 0

      !$omp parallel                                                           &
      !$omp private(id, my_id, ii, jj, i, j, k, jk, inbc, ini_excl, fin_excl,  &
      !$omp         ini_nb14, fin_nb14, nb15_calc, dij, dij_pbc, rij2,         &
      !$omp         icel)&
      !$omp firstprivate(do_allocate)
      !
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      my_id = my_city_rank * nthread + id
      id    = id + 1

      do icel = 1, ncel

        ii = head(icel)
        if (ii <= 0) cycle

        do while (ii > 0)

          i = enefunc%table%solute_list(ii)

          if (mod(ii-1,nproc_city*nthread) == my_id) then

            ! serch only in the self and neighboring cells 
            !
            do inbc = 1, boundary%num_neighbor_cells
              jk = headw(boundary%neighbor_cells(inbc,icel))

              do while (jk > 0)

                ! water list
                !
                jj = (jk+2) / 3
                k = mod(jk+2,3) + 1
                j = enefunc%table%water_list(k,jj)

                ! compute distance
                !
                dij(1:3) = coord(1:3,i) - coord(1:3,j)

                ! consider the periodic boundary
                !
                dij_pbc(1:3) = box_size(1:3)*anint(dij(1:3)/box_size(1:3))
                dij(1:3)     = dij(1:3) - dij_pbc(1:3)
                rij2         = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

                ! cout the number of interactions
                !
                if (rij2 < pairdist2) then
                  num_nb15_calcw(ii) = num_nb15_calcw(ii) + 1
                  if (.not. do_allocate) &
                    pairlist%table%nb15_calc_listw(num_nb15_calcw(ii),ii) = j
                end if

                jk = listw(jk)

              end do
            end do
          end if

          ii = list(ii)

        end do
      end do
      !$omp end parallel

      ! allocate memory of pairlist
      !
      if (do_allocate) then
        num_nb15_max = max(1,maxval(num_nb15_calcw(1:natom)))
        num_nb15_max = int(real(num_nb15_max,wp)*FactNumNb15)
        call alloc_pairlist2(pairlist, PairListPbcSoluteWater,   &
                             num_nb15_max, max(1,natom))
        pairlist%table%num_nb15_max = num_nb15_max

        do_allocate = .false.
      end if

    end do

    ! check memory size
    !
    call check_pairlist_memory_size(num_nb15_calcw,              &
                                    pairlist%table%num_nb15_max, &
                                    pairlist%allocate_solwat)

    return

  end subroutine update_pairlist_solute_water

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_water_water
  !> @brief        update pairlist between water and water
  !! @authors      JJ, TM, MK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : boundary conditions information
  !! @param[in]    coord    : coordinates
  !! @param[inout] pairlist : pairlist information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_water_water(enefunc, boundary, coord, pairlist)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord(:,:)
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: dij(1:3), rij2, pairdist2
    real(wp)                  :: dij_pbc(1:3), origin(1:3)
    real(wp)                  :: r_shift(1:3)
    real(wp)                  :: csize_inv(1:3)
    real(wp)                  :: box_size(1:3), half_box_size(1:3)
    integer                   :: ii, jj, i, j, k, n, nloops
    integer                   :: nwater
    integer                   :: ic(1:3), ncell(1:3), icel, inbc
    integer                   :: num_nb15_max
    integer                   :: id, my_id, nthread
#ifdef OMP
    integer                   :: omp_get_max_threads, omp_get_thread_num
#endif
    logical                   :: do_allocate

    integer,          pointer :: num_nb15_calc_water(:)
    integer,          pointer :: listw(:), ic_atmw(:), headw(:)
    integer,          pointer :: ncel


    nwater       =  enefunc%table%num_water
    num_nb15_calc_water => pairlist%table%num_nb15_calc_water
    ic_atmw      => pairlist%table%atom_cell_index_water
    listw        => pairlist%table%cell_linked_list_water
    headw        => boundary%cell_head_atomw
    ncel         => boundary%num_cells

    box_size(1)  =  boundary%box_size_x
    box_size(2)  =  boundary%box_size_y
    box_size(3)  =  boundary%box_size_z
    csize_inv(1) =  1.0_wp / boundary%cell_size_x
    csize_inv(2) =  1.0_wp / boundary%cell_size_y
    csize_inv(3) =  1.0_wp / boundary%cell_size_z
    ncell(1)     =  boundary%num_cells_x
    ncell(2)     =  boundary%num_cells_y
    ncell(3)     =  boundary%num_cells_z
    origin(1)    =  boundary%origin_x
    origin(2)    =  boundary%origin_y
    origin(3)    =  boundary%origin_z
    half_box_size(1:3) = box_size(1:3)*0.5_wp

    pairdist2 = (pairlist%pairlistdist + TableWaterMargin) &
              * (pairlist%pairlistdist + TableWaterMargin)

    ! initialize
    !
#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif


    ! make cell linked list for water
    !
    headw(1:ncel) = 0
    do ii = 1, nwater

      ! coordinate shifted to the first quadrant and set into the boundary box
      !
      i = enefunc%table%water_list(1,ii)
      r_shift(1:3) = coord(1:3,i) - origin(1:3)
      r_shift(1:3) = r_shift(1:3) + half_box_size(1:3) - &
                     box_size(1:3)*anint(r_shift(1:3)/box_size(1:3))

      ! assign which cell
      !
      ic(1:3) = int(r_shift(1:3)*csize_inv(1:3))

      ! upper limit of ic[x,y,z] should be ncel_[x,y,z]
      !
      do k = 1, 3
        if (ic(k) >= ncell(k)) ic(k) = ncell(k)-1
      end do
      icel = 1 + ic(1) + ic(2)*ncell(1) + ic(3)*ncell(1)*ncell(2)

      ! store which cell atom i is placed
      !
      ic_atmw(ii) = icel
      listw(ii)   = headw(icel)
      headw(icel)  = ii

    end do

    num_nb15_calc_water(1:nwater) = 0


    ! if pairlist%allocation = .true. , allocation + make pairlist
    ! if pairlist%allocation = .false., make pairlist only
    !
    if (pairlist%allocate_watwat) then
      nloops = 2
      do_allocate = .true.
    else
      nloops = 1
      do_allocate = .false.
    end if

    do n = 1, nloops

      if (.not. do_allocate) &
        num_nb15_calc_water(1:nwater) = 0

      !$omp parallel                               &
      !$omp private(id, my_id, ii, jj, i, j, inbc, &
      !$omp         dij, dij_pbc, rij2, icel)      &
      !$omp firstprivate(do_allocate)
      !
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      my_id = my_city_rank * nthread + id
      id    = id + 1

      do icel = 1, ncel

        ii = headw(icel)
        if (ii <= 0) cycle

        do while (ii > 0)

          if (mod(ii-1,nproc_city*nthread) == my_id) then

            i = enefunc%table%water_list(1,ii)

            ! serch only in the self and neighboring cells
            !
            do inbc = 1, boundary%num_neighbor_cells

              jj = headw(boundary%neighbor_cells(inbc,icel))

              do while (jj > 0)

                if (jj <= ii) then
                  jj = listw(jj)
                  cycle
                end if

                j = enefunc%table%water_list(1,jj)

                ! compute distance
                !
                dij(1:3) = coord(1:3,i) - coord(1:3,j)

                ! consider the periodic boundary
                !
                dij_pbc(1:3) = box_size(1:3)*anint(dij(1:3)/box_size(1:3))
                dij(1:3)     = dij(1:3) - dij_pbc(1:3)
                rij2         = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

                ! cout the number of interactions
                !
                if (rij2 < pairdist2) then
                  num_nb15_calc_water(ii) = num_nb15_calc_water(ii) + 1

                  if (.not. do_allocate) &
                    pairlist%table%nb15_calc_list_water(num_nb15_calc_water(ii),ii) = jj
                end if

                jj = listw(jj)
              end do
            end do
          end if

          ii = listw(ii)

        end do
      end do
      !$omp end parallel

      ! allocate memory of pairlist
      !
      if (do_allocate) then
        num_nb15_max = max(1,maxval(num_nb15_calc_water(1:nwater)))
        num_nb15_max = int(real(num_nb15_max,wp)*FactNumNb15)
        call alloc_pairlist2(pairlist, PairListPbcWaterWater,          &
                             num_nb15_max, max(1,nwater))
        pairlist%table%water_nb15_max = num_nb15_max

        do_allocate = .false.
      end if

    end do

    ! check memory size
    !
    call check_pairlist_memory_size(num_nb15_calc_water,           &
                                    pairlist%table%water_nb15_max, &
                                    pairlist%allocate_watwat)

    return

  end subroutine update_pairlist_water_water

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_coord_pbc
  !> @brief        pbc oriented coordinates
  !! @authors      JJ
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : boundary conditions information
  !! @param[in]    coord    : coordinates
  !! @param[inout] coord_pbc: pbc oriented coordinates
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_coord_pbc(enefunc, boundary, coord, trans, coord_pbc)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: trans(:,:)
    real(wp),                 intent(inout) :: coord_pbc(:,:)

    ! local variables
    integer          :: i, j, k, num_atom
    real(wp)         :: box_size(1:3), half_box_size(1:3)
    real(wp)         :: move(3), shift(3)


    box_size(1) = boundary%box_size_x
    box_size(2) = boundary%box_size_y
    box_size(3) = boundary%box_size_z
    half_box_size(1:3) = box_size(1:3) * 0.5_dp
    num_atom = size(coord(1,:))
   
    do i = 1, num_atom
      shift(1) = coord(1,i) - boundary%origin_x
      shift(2) = coord(2,i) - boundary%origin_y
      shift(3) = coord(3,i) - boundary%origin_z
      move(1) = half_box_size(1) - box_size(1)*anint(shift(1)/box_size(1))
      move(2) = half_box_size(2) - box_size(2)*anint(shift(2)/box_size(2))
      move(3) = half_box_size(3) - box_size(3)*anint(shift(3)/box_size(3))
      trans(1,i) = move(1) - boundary%origin_x
      trans(2,i) = move(2) - boundary%origin_y
      trans(3,i) = move(3) - boundary%origin_z
      coord_pbc(1,i) = shift(1) + move(1)
      coord_pbc(2,i) = shift(2) + move(2)
      coord_pbc(3,i) = shift(3) + move(3)
    end do

    return

  end subroutine update_coord_pbc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_pbc_cg_general
  !> @brief        make linked cell lists
  !! @authors      CT
  !! @param[in]    enefunc   : potential energy functions information
  !! @param[in]    boundary  : boundary conditions information
  !! @param[in]    coord_pbc : coordinates
  !! @param[inout] pairlist  : pairlist information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_pbc_cg_general(enefunc, boundary, coord_pbc, &
                                            pairlist)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord_pbc(:,:)
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    integer          :: num_atom
    integer          :: num_dna, num_base, num_phos
    integer          :: num_pro_idr_hps
    integer          :: num_pro_idr_kh
    integer          :: num_pro_kh
    integer          :: num_charged
    ! 
    integer          :: i_atom
    integer          :: i_dna, i_base, i_phos
    integer          :: i_pro_idr_hps
    integer          :: i_pro_idr_kh
    integer          :: i_pro_kh
    integer          :: i_charged
    integer          :: i, j, k
    ! 
    real(wp)         :: r_shift(1:3)
    integer          :: icell_tmp(1:3), icell_index
    ! 
    integer          :: num_cell_all
    integer          :: num_cell(1:3)
    real(wp)         :: box_origin(1:3)
    real(wp)         :: box_size(1:3), half_box_size(1:3)
    real(wp)         :: csize_inv(1:3)
    ! 
    integer, pointer :: icell_atom(:)
    integer, pointer :: cell_list_all(:), cell_head_all(:)
    integer, pointer :: cell_list_DNA(:), cell_head_DNA(:)
    integer, pointer :: cell_list_DNA_base(:), cell_head_DNA_base(:)
    integer, pointer :: cell_list_DNA_phos(:), cell_head_DNA_phos(:)
    integer, pointer :: cell_list_IDR_HPS(:), cell_head_IDR_HPS(:)
    integer, pointer :: cell_list_IDR_KH(:), cell_head_IDR_KH(:)
    integer, pointer :: cell_list_KH(:), cell_head_KH(:)
    integer, pointer :: cell_list_charged(:), cell_head_charged(:)
    ! 
    integer          :: nthread
#ifdef OMP
    integer          :: omp_get_max_threads, omp_get_thread_num
#endif


    num_atom           = size(coord_pbc(1,:))
    num_dna            = enefunc%num_cg_particle_DNA_all
    num_base           = enefunc%num_cg_particle_DNA_base
    num_phos           = enefunc%num_cg_particle_DNA_phos
    num_pro_idr_hps    = enefunc%num_cg_particle_IDR_HPS
    num_pro_idr_kh     = enefunc%num_cg_particle_IDR_KH 
    num_pro_kh         = enefunc%num_cg_particle_KH 
    num_charged        = enefunc%num_cg_particle_charged

    num_cell_all       = boundary%num_cells
    num_cell(1)        = boundary%num_cells_x
    num_cell(2)        = boundary%num_cells_y
    num_cell(3)        = boundary%num_cells_z
    box_origin(1)      = boundary%origin_x
    box_origin(2)      = boundary%origin_y
    box_origin(3)      = boundary%origin_z
    box_size(1)        = boundary%box_size_x
    box_size(2)        = boundary%box_size_y
    box_size(3)        = boundary%box_size_z
    half_box_size(1:3) = box_size(1:3) * 0.5_wp
    csize_inv(1)       = 1.0_wp / boundary%cell_size_x
    csize_inv(2)       = 1.0_wp / boundary%cell_size_y
    csize_inv(3)       = 1.0_wp / boundary%cell_size_z

    icell_atom         => pairlist%cell_index_cg_all
    cell_list_all      => pairlist%cell_linked_list_cg_all
    cell_head_all      => pairlist%cell_head_index_cg_all
    cell_list_DNA      => pairlist%cell_linked_list_cg_DNA
    cell_head_DNA      => pairlist%cell_head_index_cg_DNA
    cell_list_DNA_base => pairlist%cell_linked_list_cg_DNA_base
    cell_head_DNA_base => pairlist%cell_head_index_cg_DNA_base
    cell_list_DNA_phos => pairlist%cell_linked_list_cg_DNA_phos
    cell_head_DNA_phos => pairlist%cell_head_index_cg_DNA_phos
    cell_list_IDR_HPS  => pairlist%cell_linked_list_cg_IDR_HPS
    cell_head_IDR_HPS  => pairlist%cell_head_index_cg_IDR_HPS
    cell_list_IDR_KH   => pairlist%cell_linked_list_cg_IDR_KH
    cell_head_IDR_KH   => pairlist%cell_head_index_cg_IDR_KH
    cell_list_KH       => pairlist%cell_linked_list_cg_KH
    cell_head_KH       => pairlist%cell_head_index_cg_KH
    cell_list_charged  => pairlist%cell_linked_list_cg_charged
    cell_head_charged  => pairlist%cell_head_index_cg_charged

#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif


    ! =======================================
    ! Make cell linked list for all particles
    ! =======================================
    !
    cell_head_all(1:num_cell_all) = 0
    do i_atom = 1, num_atom

      ! wrapping into first quadrant box image
      ! 
!     r_shift(1:3) = coord(1:3,i_atom) - box_origin(1:3)
!     r_shift(1:3) = r_shift(1:3) + half_box_size(1:3) &
!         - box_size(1:3)*anint(r_shift(1:3)/box_size(1:3))
      r_shift(1:3) = coord_pbc(1:3,i_atom)

      ! which cell am i in?  get the three indices
      ! 
      icell_tmp(1:3) = int(r_shift(1:3)*csize_inv(1:3))

      ! check the upper limit of icell index
      ! 
      do k = 1, 3
        if (icell_tmp(k) >= num_cell(k)) icell_tmp(k) = num_cell(k)-1
      end do

      ! convert 3-coor to 1-coor
      ! 
      icell_index = 1 + icell_tmp(1) + (icell_tmp(2) + icell_tmp(3) * num_cell(2))*num_cell(1)

      ! --------------------------
      ! KEY: make the linked list!
      ! --------------------------
      !
      icell_atom(i_atom)         = icell_index
      cell_list_all(i_atom)      = cell_head_all(icell_index)
      cell_head_all(icell_index) = i_atom

    end do

    ! =======================================
    ! Make cell linked list for DNA particles 
    ! =======================================
    !
    cell_head_DNA(1:num_cell_all) = 0
    do i_dna = 1, num_dna
      i_atom = enefunc%cg_particle_DNA_all(i_dna)
      icell_index = icell_atom(i_atom)
      cell_list_DNA(i_dna) = cell_head_DNA(icell_index)
      cell_head_DNA(icell_index) = i_dna
    end do

    ! ===================================
    ! Make cell linked list for DNA bases
    ! ===================================
    ! 
    cell_head_DNA_base(1:num_cell_all) = 0
    do i_base = 1, num_base
      i_atom = enefunc%cg_particle_DNA_base(i_base)
      icell_index = icell_atom(i_atom)
      cell_list_DNA_base(i_base) = cell_head_DNA_base(icell_index)
      cell_head_DNA_base(icell_index) = i_base
    end do

    ! ===================================
    ! Make cell linked list for DNA phoss
    ! ===================================
    ! 
    cell_head_DNA_phos(1:num_cell_all) = 0
    do i_phos = 1, num_phos
      i_atom = enefunc%cg_particle_DNA_phos(i_phos)
      icell_index = icell_atom(i_atom)
      cell_list_DNA_phos(i_phos) = cell_head_DNA_phos(icell_index)
      cell_head_DNA_phos(icell_index) = i_phos
    end do

    ! =================================
    ! Make cell linked list for IDR HPS
    ! =================================
    ! 
    cell_head_IDR_HPS(1:num_cell_all) = 0
    do i_pro_idr_hps = 1, num_pro_idr_hps
      i_atom = enefunc%cg_particle_IDR_HPS(i_pro_idr_hps)
      icell_index = icell_atom(i_atom)
      cell_list_IDR_HPS(i_pro_idr_hps) = cell_head_IDR_HPS(icell_index)
      cell_head_IDR_HPS(icell_index) = i_pro_idr_hps
    end do

    ! =================================
    ! Make cell linked list for IDR KH
    ! =================================
    ! 
    cell_head_IDR_KH(1:num_cell_all) = 0
    do i_pro_idr_kh = 1, num_pro_idr_kh
      i_atom = enefunc%cg_particle_IDR_KH(i_pro_idr_kh)
      icell_index = icell_atom(i_atom)
      cell_list_IDR_KH(i_pro_idr_kh) = cell_head_IDR_KH(icell_index)
      cell_head_IDR_KH(icell_index) = i_pro_idr_kh
    end do

    ! =============================
    ! Make cell linked list for KH
    ! =============================
    ! 
    cell_head_KH(1:num_cell_all) = 0
    do i_pro_kh = 1, num_pro_kh
      i_atom = enefunc%cg_particle_KH(i_pro_kh)
      icell_index = icell_atom(i_atom)
      cell_list_KH(i_pro_kh) = cell_head_KH(icell_index)
      cell_head_KH(icell_index) = i_pro_kh
    end do

    ! =================================
    ! Make cell linked list for charged
    ! =================================
    ! 
    cell_head_charged(1:num_cell_all) = 0
    do i_charged = 1, num_charged
      i_atom = enefunc%cg_particle_charged(i_charged)
      icell_index = icell_atom(i_atom)
      cell_list_charged(i_charged) = cell_head_charged(icell_index)
      cell_head_charged(icell_index) = i_charged
    end do


    return

  end subroutine update_pairlist_pbc_cg_general

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_pbc_cg_exv
  !> @brief        update pairlist for CG general excluded volume interactions
  !! @authors      CT
  !! @param[in]    enefunc   : potential energy functions information
  !! @param[in]    boundary  : boundary conditions information
  !! @param[in]    coord_pbc : coordinates
  !! @param[inout] pairlist  : pairlist information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_pbc_cg_exv(enefunc, boundary, coord_pbc, pairlist)
    
    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord_pbc(:,:)
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    integer          :: num_atom
    ! 
    integer          :: i, j, k, n, nloops
    integer          :: i_atom, j_atom
    integer          :: i_cell, j_cell
    integer          :: i_nbcell
    integer          :: i_base_type, j_base_type
    integer          :: i_chain_id, j_chain_id
    integer          :: pbc_int
    ! 
    integer          :: ini_excl, fin_excl
    integer          :: num_exv_max
    ! 
    real(wp)         :: dij(1:3)
    real(wp)         :: dij_pbc(1:3)
    real(wp)         :: rij_sqr
    real(wp)         :: pairdist_exv
    real(wp)         :: pairdist_exv_sqr
    ! 
    real(wp)         :: box_size(1:3), half_box_size(1:3)
    ! 
    integer, pointer :: icell_atom(:)
    integer, pointer :: cell_list_all(:), cell_head_all(:)
    integer, pointer :: num_exv_pre(:), num_exv(:)
    ! 
    logical          :: do_allocate
    logical          :: do_calc
    logical          :: is_dna_i
    logical          :: is_dna_j
    logical          :: is_idr_i
    logical          :: is_idr_j
    ! 
    integer          :: id, my_id, nthread
#ifdef OMP
    integer          :: omp_get_max_threads, omp_get_thread_num
#endif


    num_atom           = size(coord_pbc(1,:))

    box_size(1)        = boundary%box_size_x
    box_size(2)        = boundary%box_size_y
    box_size(3)        = boundary%box_size_z
    half_box_size(1:3) = box_size(1:3) * 0.5_wp

    pairdist_exv       = pairlist%cg_pairlistdist_exv
    pairdist_exv_sqr   = pairdist_exv * pairdist_exv

    icell_atom         => pairlist%cell_index_cg_all
    cell_list_all      => pairlist%cell_linked_list_cg_all
    cell_head_all      => pairlist%cell_head_index_cg_all
    num_exv_pre        => pairlist%num_cg_exv_pre
    num_exv            => pairlist%num_cg_exv

#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif


    pairlist%num_cg_exv_calc(1:num_atom,1:nthread) = 0

    if (pairlist%allocate_pbc_cg_exv) then
      nloops      = 2
      do_allocate = .true.
    else
      nloops      = 1
      do_allocate = .false.
    end if

    do n = 1, nloops

      num_exv(1:nthread) = 0

      if (.not. do_allocate) num_exv_pre(1:nthread) = 0

      !$omp parallel                          &
      !$omp private(id, my_id, i, j, k,       &
      !$omp         i_atom, j_atom,           &
      !$omp         i_cell, j_cell,           &
      !$omp         is_dna_i, is_dna_j,       &
      !$omp         is_idr_i, is_idr_j,       &
      !$omp         i_nbcell, do_calc,        &
      !$omp         ini_excl, fin_excl,       &
      !$omp         dij, dij_pbc, rij_sqr,    &
      !$omp         i_base_type, j_base_type, &
      !$omp         i_chain_id, j_chain_id,   &
      !$omp         pbc_int)                  &
      !$omp firstprivate(do_allocate)
      !
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      my_id = my_city_rank * nthread + id
      id    = id + 1
      
      do i_atom = 1, num_atom-1

        if (mod(i_atom - 1, nproc_city * nthread) == my_id) then

          i_cell = icell_atom(i_atom)

          i_chain_id  = enefunc%mol_chain_id(i_atom)
          i_base_type = enefunc%NA_base_type(i_atom)
          if (i_base_type <= NABaseTypeDBMAX &
              .or. i_base_type == NABaseTypeDP &
              .or. i_base_type == NABaseTypeDS) then
            is_dna_i = .true.
          else
            is_dna_i = .false.
          end if
          if (enefunc%cg_IDR_HPS_is_IDR(i_atom) .or. enefunc%cg_IDR_KH_is_IDR(i_atom)) then
            is_idr_i = .true.
          else
            is_idr_i = .false.
          end if

          do i_nbcell = 1, boundary%num_neighbor_cells_CG_exv

            j_cell = boundary%neighbor_cells_CG_exv(i_nbcell, i_cell)
            j_atom = cell_head_all(j_cell)

            do while (j_atom > i_atom)

              ! 
              ! exclude DNA-DNA interactions
              ! 
              if (is_dna_i) then
                j_base_type = enefunc%NA_base_type(j_atom)
                if (j_base_type <= NABaseTypeDBMAX &
                    .or. j_base_type == NABaseTypeDP &
                    .or. j_base_type == NABaseTypeDS) then
                  is_dna_j = .true.
                else
                  is_dna_j = .false.
                end if
                ! 
                if (is_dna_j) then
                  j_atom = cell_list_all(j_atom)
                  cycle
                end if
              end if

              ! exclude IDR-IDR interactions
              ! 
              if (enefunc%cg_IDR_HPS_is_IDR(j_atom) .or. &
                  enefunc%cg_IDR_KH_is_IDR(j_atom)) then
                is_idr_j = .true.
              else
                is_idr_j = .false.
              end if
              ! 
              if (is_idr_i .and. is_idr_j) then
                j_atom = cell_list_all(j_atom)
                cycle
              end if

              j_chain_id  = enefunc%mol_chain_id(j_atom)
              !
              ! exclude KH interaction
              ! 
              if (enefunc%cg_kh_mol_pair(i_chain_id, j_chain_id) > 0 .and. &
                  enefunc%cg_pro_use_KH(i_atom) .and. &
                  enefunc%cg_pro_use_KH(j_atom)) then
                if (j_chain_id /= i_chain_id) then
                  j_atom = cell_list_all(j_atom)
                  cycle
                else if (is_idr_i .or. is_idr_j) then
                  j_atom = cell_list_all(j_atom)
                  cycle
                end if
              end if

              do_calc = .true.
              ! 
              ! exclude interactions in exclusion list
              ! 
              ini_excl = enefunc%cg_istart_nonb_excl(i_atom)
              fin_excl = ini_excl + enefunc%num_nonb_excl(i_atom) - 1
              do k = ini_excl, fin_excl 
                if (j_atom == enefunc%nonb_excl_list(k)) then
                  do_calc = .false.
                  exit
                end if
              end do
              if (.not. do_calc) then
                j_atom = cell_list_all(j_atom)
                cycle
              end if

              ! 
              ! distance within pairlistdist?
              ! 
              dij(1:3) = coord_pbc(1:3,j_atom) - coord_pbc(1:3,i_atom)
              call check_pbc(box_size, dij, pbc_int)

              rij_sqr  = dij(1) * dij(1) + dij(2) * dij(2) + dij(3) * dij(3)
              
              ! --------------------
              ! core: count and fill
              ! --------------------
              ! 
              if (rij_sqr < pairdist_exv_sqr) then
                num_exv(id) = num_exv(id) + 1
                if (.not. do_allocate) &
                    pairlist%cg_exv_list(num_exv(id), id) = pbc_int + j_atom*27
              end if

              j_atom = cell_list_all(j_atom)
            end do              ! j_atom

          end do                ! i_nbcell

        end if                  ! compute on the current thread

        if (.not. do_allocate) then
          pairlist%num_cg_exv_calc(i_atom,id) = num_exv(id) - num_exv_pre(id)
          num_exv_pre(id) = num_exv(id)
        end if

      end do                    ! i = 1 to num_atom - 1
      !$omp end parallel

      ! allocate memory of pairlist
      !
      if (do_allocate) then
        num_exv_max = max(1, maxval(num_exv(1:nthread)))
        num_exv_max = int(real(num_exv_max,wp) * FactNumNb15)

        call alloc_pairlist(pairlist, PairListPbcCGexv, num_exv_max)

        pairlist%num_cg_exv_max = num_exv_max

        do_allocate = .false.
      end if
      
    end do                      ! nloops

    ! check memory size
    !
    call check_pairlist_memory_size(num_exv, pairlist%num_cg_exv_max, &
        pairlist%allocate_pbc_cg_exv)

    return

  end subroutine update_pairlist_pbc_cg_exv

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_pbc_cg_ele
  !> @brief        update pairlist for CG electrostatics
  !! @authors      CT
  !! @param[in]    enefunc   : potential energy functions information
  !! @param[in]    boundary  : boundary conditions information
  !! @param[in]    coord_pbc : coordinates
  !! @param[inout] pairlist  : pairlist information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_pbc_cg_ele(enefunc, boundary, coord_pbc, pairlist)
    
    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord_pbc(:,:)
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    integer          :: num_charge
    ! 
    integer          :: i, j, k, n, nloops
    integer          :: i_charge, j_charge
    integer          :: i_atom, j_atom
    integer          :: i_cell, j_cell
    integer          :: i_nbcell
    integer          :: i_base_type, j_base_type
    integer          :: i_chain_id, j_chain_id
    integer          :: pbc_int
    ! 
    integer          :: ini_excl, fin_excl
    integer          :: num_ele_max
    ! 
    real(wp)         :: dij(1:3)
    real(wp)         :: dij_pbc(1:3)
    real(wp)         :: rij_sqr
    real(wp)         :: pairdist_ele
    real(wp)         :: pairdist_ele_sqr
    ! 
    real(wp)         :: box_size(1:3), half_box_size(1:3)
    ! 
    integer, pointer :: cg_list_charged(:)
    integer, pointer :: icell_atom(:)
    integer, pointer :: cell_list_charged(:), cell_head_charged(:)
    integer, pointer :: num_ele_pre(:), num_ele(:)
    ! 
    logical          :: do_allocate
    logical          :: do_calc
    logical          :: is_phos_i
    logical          :: is_phos_j
    logical          :: is_pro_i
    logical          :: is_pro_j
    ! 
    integer          :: id, my_id, nthread
#ifdef OMP
    integer          :: omp_get_max_threads, omp_get_thread_num
#endif


    num_charge         = enefunc%num_cg_particle_charged

    box_size(1)        = boundary%box_size_x
    box_size(2)        = boundary%box_size_y
    box_size(3)        = boundary%box_size_z
    half_box_size(1:3) = box_size(1:3) * 0.5_wp

    pairdist_ele       = pairlist%cg_pairlistdist_ele
    pairdist_ele_sqr   = pairdist_ele * pairdist_ele

    icell_atom         => pairlist%cell_index_cg_all
    cell_list_charged  => pairlist%cell_linked_list_cg_charged
    cell_head_charged  => pairlist%cell_head_index_cg_charged
    num_ele_pre        => pairlist%num_cg_ele_pre
    num_ele            => pairlist%num_cg_ele
    cg_list_charged    => enefunc%cg_particle_charged

#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif


    pairlist%num_cg_ele_calc(1:num_charge,1:nthread) = 0

    if (pairlist%allocate_pbc_cg_ele) then
      nloops      = 2
      do_allocate = .true.
    else
      nloops      = 1
      do_allocate = .false.
    end if

    do n = 1, nloops

      num_ele(1:nthread) = 0

      if (.not. do_allocate) num_ele_pre(1:nthread) = 0

      !$omp parallel                          &
      !$omp private(id, my_id, i, j, k,       &
      !$omp         i_charge, j_charge,       &
      !$omp         i_atom, j_atom,           &
      !$omp         i_cell, j_cell,           &
      !$omp         is_phos_i, is_phos_j,     &
      !$omp         is_pro_i, is_pro_j,       &
      !$omp         i_nbcell, do_calc,        &
      !$omp         ini_excl, fin_excl,       &
      !$omp         dij, dij_pbc, rij_sqr,    &
      !$omp         i_chain_id, j_chain_id,   &
      !$omp         i_base_type, j_base_type, &
      !$omp         pbc_int)                  &
      !$omp firstprivate(do_allocate)
      !
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      my_id = my_city_rank * nthread + id
      id    = id + 1
      
      do i_charge = 1, num_charge - 1

        if (mod(i_charge - 1, nproc_city * nthread) == my_id) then

          i_atom = cg_list_charged(i_charge)
          i_cell      = icell_atom(i_atom)
          i_base_type = enefunc%NA_base_type(i_atom)
          i_chain_id  = enefunc%mol_chain_id(i_atom)
          is_phos_i   = .false.
          is_pro_i    = .false.
          if (i_base_type == NABaseTypeDP) then
            is_phos_i = .true.
          else if (i_base_type == NABaseTypeProtein) then
            is_pro_i  = .true.
          end if

          do i_nbcell = 1, boundary%num_neighbor_cells_CG_ele

            j_cell = boundary%neighbor_cells_CG_ele(i_nbcell, i_cell)
            j_charge = cell_head_charged(j_cell)

            do while (j_charge > i_charge)

              j_atom = cg_list_charged(j_charge)

              !
              ! get j properties
              ! 
              j_base_type = enefunc%NA_base_type(j_atom)
              j_chain_id  = enefunc%mol_chain_id(j_atom)
              is_phos_j   = .false.
              is_pro_j    = .false.
              if (j_base_type == NABaseTypeDP) then
                is_phos_j = .true.
              else if (j_base_type == NABaseTypeProtein) then
                is_pro_j  = .true.
              end if

              ! 
              ! check inter-molecule ele interactions
              ! 
              if (enefunc%cg_ele_mol_pair(i_chain_id, j_chain_id) == 0) then
                j_charge = cell_list_charged(j_charge)
                cycle
              end if

              do_calc = .true.
              !
              ! exclusion list
              ini_excl = enefunc%cg_istart_nonb_excl(i_atom)
              fin_excl = ini_excl + enefunc%num_nonb_excl(i_atom) - 1
              do k = ini_excl, fin_excl
                if (j_atom == enefunc%nonb_excl_list(k)) then
                  do_calc = .false.
                  exit
                end if
              end do
              !
              if (.not. do_calc) then
                j_charge = cell_list_charged(j_charge)
                cycle
              end if

              ! 
              ! distance within pairlistdist?
              ! 
              dij(1:3)     = coord_pbc(1:3,j_atom) - coord_pbc(1:3,i_atom)
              call check_pbc(box_size, dij, pbc_int)
!             dij_pbc(1:3) = box_size(1:3) * anint(dij(1:3) / box_size(1:3))
!             dij(1:3)     = dij(1:3) - dij_pbc(1:3)
              rij_sqr      = dij(1) * dij(1) + dij(2) * dij(2) + dij(3) * dij(3)
              
              ! --------------------
              ! core: count and fill
              ! --------------------
              ! 
              if (rij_sqr < pairdist_ele_sqr) then
                num_ele(id) = num_ele(id) + 1
                if (.not. do_allocate) then
                  pairlist%cg_ele_list(num_ele(id), id) = pbc_int + 27*j_atom

                  if ((is_phos_i .and. is_pro_j) .or. &
                      (is_phos_j .and. is_pro_i)) then
                    pairlist%cg_ele_scaling_list(num_ele(id), id) = &
                        - enefunc%cg_pro_DNA_ele_scale_Q / 0.6
                  else
                    pairlist%cg_ele_scaling_list(num_ele(id), id) = 1.0_wp
                  end if
                end if
              end if

              j_charge = cell_list_charged(j_charge)
            end do              ! j_atom

          end do                ! i_nbcell

        end if                  ! compute on the current thread

        if (.not. do_allocate) then
          pairlist%num_cg_ele_calc(i_charge,id) = num_ele(id) - num_ele_pre(id)
          num_ele_pre(id) = num_ele(id)
        end if

      end do                    ! icharge = 1 to num_charge
      !$omp end parallel

      ! allocate memory of pairlist
      !
      if (do_allocate) then
        num_ele_max = max(1, maxval(num_ele(1:nthread)))
        num_ele_max = int(real(num_ele_max,wp) * FactNumNb15)

        call alloc_pairlist(pairlist, PairListPbcCGele, num_ele_max)

        pairlist%num_cg_ele_max = num_ele_max

        do_allocate = .false.
      end if
      
    end do                      ! nloops

    ! check memory size
    !
    call check_pairlist_memory_size(num_ele, pairlist%num_cg_ele_max, &
        pairlist%allocate_pbc_cg_ele)

    return

  end subroutine update_pairlist_pbc_cg_ele

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_pbc_cg_DNA_base
  !> @brief        update pairlist for CG DNA basepairing and base crossstacking
  !! @authors      CT
  !! @param[in]    enefunc   : potential energy functions information
  !! @param[in]    boundary  : boundary conditions information
  !! @param[in]    coord_pbc : PBC coordinates
  !! @param[inout] pairlist  : pairlist information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_pbc_cg_DNA_base(enefunc, boundary, coord_pbc, &
                                             pairlist)
    
    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord_pbc(:,:)
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    integer          :: num_base
    ! 
    integer          :: i, j, k, n, nloops
    integer          :: i_base, j_base
    integer          :: i_atom, j_atom
    integer          :: i_cell, j_cell
    integer          :: i_nbcell
    integer          :: i_base_type, j_base_type
    integer          :: i_chain_id, j_chain_id
    integer          :: pbc_int
    ! 
    integer          :: ini_excl, fin_excl
    integer          :: num_DNA_bp_max
    ! 
    real(wp)         :: dij(1:3)
    real(wp)         :: dij_pbc(1:3)
    real(wp)         :: rij_sqr
    real(wp)         :: pairdist_DNAbp
    real(wp)         :: pairdist_DNAbp_sqr
    ! 
    real(wp)         :: box_size(1:3), half_box_size(1:3)
    ! 
    integer, pointer :: cg_list_base(:)
    integer, pointer :: icell_atom(:)
    integer, pointer :: cell_list_DNA_base(:), cell_head_DNA_base(:)
    integer, pointer :: num_DNA_bp_pre(:), num_DNA_bp(:)
    ! 
    logical          :: do_allocate
    logical          :: do_calc
    logical          :: is_WC_bp
    ! 
    integer          :: id, my_id, nthread
#ifdef OMP
    integer          :: omp_get_max_threads, omp_get_thread_num
#endif

    num_base           = enefunc%num_cg_particle_DNA_base

    box_size(1)        = boundary%box_size_x
    box_size(2)        = boundary%box_size_y
    box_size(3)        = boundary%box_size_z
    half_box_size(1:3) = box_size(1:3) * 0.5_wp

    pairdist_DNAbp     = pairlist%cg_pairlistdist_DNAbp
    pairdist_DNAbp_sqr = pairdist_DNAbp * pairdist_DNAbp

    icell_atom         => pairlist%cell_index_cg_all
    cell_list_DNA_base => pairlist%cell_linked_list_cg_DNA_base
    cell_head_DNA_base => pairlist%cell_head_index_cg_DNA_base
    num_DNA_bp_pre     => pairlist%num_cg_DNA_basepair_pre
    num_DNA_bp         => pairlist%num_cg_DNA_basepair
    cg_list_base       => enefunc%cg_particle_DNA_base

#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif


    pairlist%num_cg_DNA_basepair_calc(1:num_base,1:nthread) = 0

    if (pairlist%allocate_pbc_cg_DNA_bp) then
      nloops      = 2
      do_allocate = .true.
    else
      nloops      = 1
      do_allocate = .false.
    end if

    do n = 1, nloops

      num_DNA_bp(1:nthread) = 0

      if (.not. do_allocate) num_DNA_bp_pre(1:nthread) = 0

      !$omp parallel                          &
      !$omp private(id, my_id, i, j, k,       &
      !$omp         i_base, j_base,           &
      !$omp         i_atom, j_atom,           &
      !$omp         i_cell, j_cell,           &
      !$omp         i_nbcell, do_calc,        &
      !$omp         ini_excl, fin_excl,       &
      !$omp         dij, dij_pbc, rij_sqr,    &
      !$omp         is_WC_bp,                 &
      !$omp         i_chain_id, j_chain_id,   &
      !$omp         i_base_type, j_base_type, &
      !$omp         pbc_int)                  &
      !$omp firstprivate(do_allocate)
      !
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      my_id = my_city_rank * nthread + id
      id    = id + 1
      
      do i_base = 1, num_base - 1

        if (mod(i_base - 1, nproc_city * nthread) == my_id) then

          i_atom = cg_list_base(i_base)
          i_cell      = icell_atom(i_atom)
          i_base_type = enefunc%NA_base_type(i_atom)
          i_chain_id  = enefunc%mol_chain_id(i_atom)

          do i_nbcell = 1, boundary%num_neighbor_cells_CG_DNAbp

            j_cell = boundary%neighbor_cells_CG_DNAbp(i_nbcell, i_cell)
            j_base = cell_head_DNA_base(j_cell)

            do while (j_base > i_base)

              j_atom = cg_list_base(j_base)

              j_base_type = enefunc%NA_base_type(j_atom)
              j_chain_id  = enefunc%mol_chain_id(j_atom)

              do_calc = .true.
              !
              ! exclude interactions in exclusion list
              !
              ini_excl = enefunc%cg_istart_nonb_excl(i_atom)
              fin_excl = ini_excl + enefunc%num_nonb_excl(i_atom) - 1
              do k = ini_excl, fin_excl
                if (j_atom == enefunc%nonb_excl_list(k)) then
                  do_calc = .false.
                  exit
                end if
              end do
              !
              if (.not. do_calc) then
                j_base = cell_list_DNA_base(j_base)
                cycle
              end if
              ! 
              is_WC_bp    = .false.
              if (i_chain_id /= j_chain_id) then
                is_WC_bp = enefunc%base_pair_is_WC(i_base_type, j_base_type)
              else
                if (j_atom - i_atom > 10) &
                    is_WC_bp = enefunc%base_pair_is_WC(i_base_type, j_base_type)
              end if
              ! 
              if (.not. is_WC_bp) then
                j_base = cell_list_DNA_base(j_base)
                cycle
              end if

              ! 
              ! distance within pairlistdist?
              ! 
              dij(1:3)     = coord_pbc(1:3,j_atom) - coord_pbc(1:3,i_atom)
              call check_pbc(box_size, dij, pbc_int)

              rij_sqr      = dij(1) * dij(1) + dij(2) * dij(2) + dij(3) * dij(3)
              
              ! --------------------
              ! core: count and fill
              ! --------------------
              ! 
              if (rij_sqr < pairdist_DNAbp_sqr) then
                num_DNA_bp(id) = num_DNA_bp(id) + 1
                if (.not. do_allocate) then
                  pairlist%cg_DNA_basepair_list(num_DNA_bp(id), id) = j_atom * 27 + pbc_int
                end if
              end if

              j_base = cell_list_DNA_base(j_base)
            end do              ! j_base

          end do                ! i_nbcell

        end if                  ! compute on the current thread

        if (.not. do_allocate) then
          pairlist%num_cg_DNA_basepair_calc(i_base,id) = num_DNA_bp(id) - num_DNA_bp_pre(id)
          num_DNA_bp_pre(id) = num_DNA_bp(id)
        end if

      end do                    !
      !$omp end parallel

      ! allocate memory of pairlist
      !
      if (do_allocate) then
        num_DNA_bp_max = max(1, maxval(num_DNA_bp(1:nthread)))
        num_DNA_bp_max = int(real(num_DNA_bp_max,wp) * FactNumNb15)

        call alloc_pairlist(pairlist, PairListPbcCGDNABP, num_DNA_bp_max)

        pairlist%num_cg_DNA_basepair_max = num_DNA_bp_max

        do_allocate = .false.
      end if
      
    end do                      ! nloops

    ! check memory size
    !
    call check_pairlist_memory_size(num_DNA_bp, &
         pairlist%num_cg_DNA_basepair_max, pairlist%allocate_pbc_cg_DNA_bp)

    return

  end subroutine update_pairlist_pbc_cg_DNA_base

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_pbc_cg_DNA_exv
  !> @brief        update pairlist for CG DNA special excluded volume
  !! @authors      CT
  !! @param[in]    enefunc   : potential energy functions information
  !! @param[in]    boundary  : boundary conditions information
  !! @param[in]    coord_pbc : PBC coordinates
  !! @param[inout] pairlist  : pairlist information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_pbc_cg_DNA_exv(enefunc, boundary, coord_pbc, &
                                            pairlist)
    
    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord_pbc(:,:)
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    integer          :: num_dna
    ! 
    integer          :: i, j, k, n, nloops
    integer          :: i_dna, j_dna
    integer          :: i_atom, j_atom
    integer          :: i_cell, j_cell
    integer          :: i_nbcell
    integer          :: i_base_type, j_base_type
    integer          :: i_chain_id, j_chain_id
    integer          :: pbc_int
    ! 
    integer          :: ini_excl, fin_excl
    integer          :: num_DNA_exv_max
    ! 
    real(wp)         :: dij(1:3)
    real(wp)         :: dij_pbc(1:3)
    real(wp)         :: rij_sqr
    real(wp)         :: pairdist_exv
    real(wp)         :: pairdist_exv_sqr
    ! 
    real(wp)         :: box_size(1:3), half_box_size(1:3)
    ! 
    integer, pointer :: cg_list_dna(:)
    integer, pointer :: icell_atom(:)
    integer, pointer :: cell_list_DNA(:), cell_head_DNA(:)
    integer, pointer :: num_DNA_exv_pre(:), num_DNA_exv(:)
    ! 
    logical          :: do_allocate
    logical          :: do_calc
    logical          :: is_nonlocal
    logical          :: is_WC_bp
    ! 
    integer          :: id, my_id, nthread
#ifdef OMP
    integer          :: omp_get_max_threads, omp_get_thread_num
#endif

    num_dna            = enefunc%num_cg_particle_DNA_all

    box_size(1)        = boundary%box_size_x
    box_size(2)        = boundary%box_size_y
    box_size(3)        = boundary%box_size_z
    half_box_size(1:3) = box_size(1:3) * 0.5_wp

    pairdist_exv       = pairlist%cg_pairlistdist_exv
    pairdist_exv_sqr   = pairdist_exv * pairdist_exv

    icell_atom         => pairlist%cell_index_cg_all
    cell_list_DNA      => pairlist%cell_linked_list_cg_DNA
    cell_head_DNA      => pairlist%cell_head_index_cg_DNA
    num_DNA_exv_pre    => pairlist%num_cg_DNA_exv_pre
    num_DNA_exv        => pairlist%num_cg_DNA_exv
    cg_list_dna        => enefunc%cg_particle_DNA_all

#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif


    pairlist%num_cg_DNA_exv_calc(1:num_dna,1:nthread) = 0

    if (pairlist%allocate_pbc_cg_DNA_exv) then
      nloops      = 2
      do_allocate = .true.
    else
      nloops      = 1
      do_allocate = .false.
    end if

    do n = 1, nloops

      num_DNA_exv(1:nthread) = 0

      if (.not. do_allocate) num_DNA_exv_pre(1:nthread) = 0

      !$omp parallel                          &
      !$omp private(id, my_id, i, j, k,       &
      !$omp         i_dna, j_dna,             &
      !$omp         i_atom, j_atom,           &
      !$omp         i_cell, j_cell,           &
      !$omp         i_nbcell, do_calc,        &
      !$omp         ini_excl, fin_excl,       &
      !$omp         dij, dij_pbc, rij_sqr,    &
      !$omp         is_nonlocal, is_WC_bp,    &
      !$omp         i_chain_id, j_chain_id,   &
      !$omp         i_base_type, j_base_type, &
      !$omp         pbc_int)                  &
      !$omp firstprivate(do_allocate)
      !
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      my_id = my_city_rank * nthread + id
      id    = id + 1
      
      do i_dna = 1, num_dna - 1

        if (mod(i_dna - 1, nproc_city * nthread) == my_id) then

          i_atom      = cg_list_dna(i_dna)
          i_cell      = icell_atom(i_atom)
          i_base_type = enefunc%NA_base_type(i_atom)
          i_chain_id  = enefunc%mol_chain_id(i_atom)

          do i_nbcell = 1, boundary%num_neighbor_cells_CG_exv

            j_cell = boundary%neighbor_cells_CG_exv(i_nbcell, i_cell)
            j_dna  = cell_head_DNA(j_cell)

            do while (j_dna > i_dna)

              j_atom = cg_list_dna(j_dna)

              j_base_type = enefunc%NA_base_type(j_atom)
              j_chain_id  = enefunc%mol_chain_id(j_atom)

              do_calc = .true.
              !
              ! exclude interactions in exclusion list
              !
              ini_excl = enefunc%cg_istart_nonb_excl(i_atom)
              fin_excl = ini_excl + enefunc%num_nonb_excl(i_atom) - 1
              do k = ini_excl, fin_excl
                if (j_atom == enefunc%nonb_excl_list(k)) then
                  do_calc = .false.
                  exit
                end if
              end do
              !
              if (.not. do_calc) then
                j_dna = cell_list_DNA(j_dna)
                cycle
              end if
              ! 
              is_nonlocal = .false.
              is_WC_bp    = .false.
              if (i_chain_id /= j_chain_id) then
                is_nonlocal = .true.
                if (i_base_type <= NABaseTypeDBMAX .and. &
                    j_base_type <= NABaseTypeDBMAX) then
                  is_WC_bp = enefunc%base_pair_is_WC(i_base_type, j_base_type)
                end if
              else
                if (j_atom - i_atom > 4) &
                    is_nonlocal = .true.
                if (j_atom - i_atom > 10) then
                  if (i_base_type <= NABaseTypeDBMAX .and. &
                      j_base_type <= NABaseTypeDBMAX) then
                    is_WC_bp = enefunc%base_pair_is_WC(i_base_type, j_base_type)
                  end if
                end if
              end if
              ! 
              if (is_WC_bp) then
                j_dna = cell_list_DNA(j_dna)
                cycle
              end if
              if (.not. is_nonlocal) then
                j_dna = cell_list_DNA(j_dna)
                cycle
              end if

              ! 
              ! distance within pairlistdist?
              ! 
              dij(1:3)     = coord_pbc(1:3,j_atom) - coord_pbc(1:3,i_atom)
              call check_pbc(box_size, dij, pbc_int)

              rij_sqr      = dij(1) * dij(1) + dij(2) * dij(2) + dij(3) * dij(3)
              
              ! --------------------
              ! core: count and fill
              ! --------------------
              ! 
              if (rij_sqr < pairdist_exv_sqr) then
                num_DNA_exv(id) = num_DNA_exv(id) + 1
                if (.not. do_allocate) &
                    pairlist%cg_DNA_exv_list(num_DNA_exv(id), id) = j_atom * 27 + pbc_int
              end if

              j_dna = cell_list_DNA(j_dna)
            end do              ! j_atom

          end do                ! i_nbcell

        end if                  ! compute on the current thread

        if (.not. do_allocate) then
          pairlist%num_cg_DNA_exv_calc(i_dna,id) = num_DNA_exv(id) - num_DNA_exv_pre(id)
          num_DNA_exv_pre(id) = num_DNA_exv(id)
        end if

      end do                    !
      !$omp end parallel

      ! allocate memory of pairlist
      !
      if (do_allocate) then
        num_DNA_exv_max = max(1, maxval(num_DNA_exv(1:nthread)))
        num_DNA_exv_max = int(real(num_DNA_exv_max,wp) * FactNumNb15)

        call alloc_pairlist(pairlist, PairListPbcCGDNAexv, num_DNA_exv_max)

        pairlist%num_cg_DNA_exv_max = num_DNA_exv_max

        do_allocate = .false.
      end if
      
    end do                      ! nloops

    ! check memory size
    !
    call check_pairlist_memory_size(num_DNA_exv, pairlist%num_cg_DNA_exv_max, &
        pairlist%allocate_pbc_cg_DNA_exv)

    return

  end subroutine update_pairlist_pbc_cg_DNA_exv

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_pbc_cg_IDR_HPS
  !> @brief        update pairlist for CG protein HPS IDR model
  !! @authors      CT
  !! @param[in]    enefunc   : potential energy functions information
  !! @param[in]    boundary  : boundary conditions information
  !! @param[in]    coord_pbc : PBC coordinates
  !! @param[inout] pairlist  : pairlist information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine update_pairlist_pbc_cg_IDR_HPS(enefunc, boundary, coord_pbc, &
                                            pairlist)
    
    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord_pbc(:,:)
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    integer          :: num_pro_idr
    ! 
    integer          :: i, j, k, n, nloops
    integer          :: i_idr, j_idr
    integer          :: i_atom, j_atom
    integer          :: i_cell, j_cell
    integer          :: i_nbcell
    integer          :: i_chain_id, j_chain_id
    integer          :: pbc_int
    ! 
    integer          :: ini_excl, fin_excl
    integer          :: num_idr_max
    ! 
    real(wp)         :: dij(1:3)
    real(wp)         :: dij_pbc(1:3)
    real(wp)         :: rij_sqr
    real(wp)         :: pairdist_126
    real(wp)         :: pairdist_126_sqr
    ! 
    real(wp)         :: box_size(1:3), half_box_size(1:3)
    ! 
    integer, pointer :: cg_list_idr(:)
    integer, pointer :: icell_atom(:)
    integer, pointer :: cell_list_idr(:), cell_head_idr(:)
    integer, pointer :: num_idr_pre(:), num_idr(:)
    ! 
    logical          :: do_allocate
    logical          :: do_calc
    ! 
    integer          :: id, my_id, nthread
#ifdef OMP
    integer          :: omp_get_max_threads, omp_get_thread_num
#endif

    num_pro_idr        = enefunc%num_cg_particle_IDR_HPS

    box_size(1)        = boundary%box_size_x
    box_size(2)        = boundary%box_size_y
    box_size(3)        = boundary%box_size_z
    half_box_size(1:3) = box_size(1:3) * 0.5_wp

    pairdist_126       = pairlist%cg_pairlistdist_126
    pairdist_126_sqr   = pairdist_126 * pairdist_126

    icell_atom         => pairlist%cell_index_cg_all
    cell_list_idr      => pairlist%cell_linked_list_cg_IDR_HPS
    cell_head_idr      => pairlist%cell_head_index_cg_IDR_HPS
    num_idr_pre        => pairlist%num_cg_IDR_HPS_pre
    num_idr            => pairlist%num_cg_IDR_HPS
    cg_list_idr        => enefunc%cg_particle_IDR_HPS

#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif


    pairlist%num_cg_IDR_HPS_calc(1:num_pro_idr,1:nthread) = 0

    if (pairlist%allocate_pbc_cg_IDR_HPS) then
      nloops      = 2
      do_allocate = .true.
    else
      nloops      = 1
      do_allocate = .false.
    end if

    do n = 1, nloops

      num_idr(1:nthread) = 0

      if (.not. do_allocate) num_idr_pre(1:nthread) = 0

      !$omp parallel                        &
      !$omp private(id, my_id, i, j, k,     &
      !$omp         pbc_int,                &
      !$omp         i_idr, j_idr,           &
      !$omp         i_atom, j_atom,         &
      !$omp         i_cell, j_cell,         &
      !$omp         i_nbcell, do_calc,      &
      !$omp         ini_excl, fin_excl,     &
      !$omp         dij, dij_pbc, rij_sqr,  &
      !$omp         i_chain_id, j_chain_id) &
      !$omp firstprivate(do_allocate)
      !
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      my_id = my_city_rank * nthread + id
      id    = id + 1
      
      do i_idr = 1, num_pro_idr - 1


        if (mod(i_idr - 1, nproc_city * nthread) == my_id) then

          i_atom     = cg_list_idr(i_idr)
          i_cell     = icell_atom(i_atom)
          i_chain_id = enefunc%mol_chain_id(i_atom)

          do i_nbcell = 1, boundary%num_neighbor_cells_CG_126

            j_cell = boundary%neighbor_cells_CG_126(i_nbcell, i_cell)
            j_idr  = cell_head_idr(j_cell)

            do while (j_idr > i_idr)

              j_atom = cg_list_idr(j_idr)

              !
              ! get j properties
              ! 
              j_chain_id  = enefunc%mol_chain_id(j_atom)
              
              do_calc = .true.
              ! 
              if (i_chain_id == j_chain_id .and. i_atom == j_atom - 1) then
                do_calc = .false.
              end if
              ! 
              if (.not. do_calc) then
                j_idr = cell_list_idr(j_idr)
                cycle
              end if

              ! 
              ! distance within pairlistdist?
              ! 
              dij(1:3)     = coord_pbc(1:3,j_atom) - coord_pbc(1:3,i_atom)
              call check_pbc(box_size, dij, pbc_int)

              rij_sqr      = dij(1) * dij(1) + dij(2) * dij(2) + dij(3) * dij(3)
              
              ! --------------------
              ! core: count and fill
              ! --------------------
              ! 
              if (rij_sqr < pairdist_126_sqr) then
                num_idr(id) = num_idr(id) + 1
                if (.not. do_allocate) then
                  pairlist%cg_IDR_HPS_list(num_idr(id), id) = j_atom * 27 + pbc_int
                end if
              end if

              j_idr = cell_list_idr(j_idr)
            end do              ! j_atom

          end do                ! i_nbcell

        end if                  ! compute on the current thread

        if (.not. do_allocate) then
          pairlist%num_cg_IDR_HPS_calc(i_idr,id) = num_idr(id) - num_idr_pre(id)
          num_idr_pre(id) = num_idr(id)
        end if

      end do                    ! icharge = 1 to num_charge
      !$omp end parallel

      ! allocate memory of pairlist
      !
      if (do_allocate) then
        num_idr_max = max(1, maxval(num_idr(1:nthread)))
        num_idr_max = int(real(num_idr_max,wp) * FactNumNb15)

        call alloc_pairlist(pairlist, PairListPbcCGIDRHPS, num_idr_max)

        pairlist%num_cg_IDR_HPS_max = num_idr_max

        do_allocate = .false.
      end if
      
    end do                      ! nloops

    ! check memory size
    !
    call check_pairlist_memory_size(num_idr, pairlist%num_cg_IDR_HPS_max, &
        pairlist%allocate_pbc_cg_IDR_HPS)

    return

  end subroutine update_pairlist_pbc_cg_IDR_HPS

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_pbc_cg_IDR_KH
  !> @brief        update pairlist for CG protein KH IDR model
  !! @authors      CT
  !! @param[in]    enefunc   : potential energy functions information
  !! @param[in]    boundary  : boundary conditions information
  !! @param[in]    coord_pbc : PBC coordinates
  !! @param[inout] pairlist  : pairlist information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_pbc_cg_IDR_KH(enefunc, boundary, coord_pbc, &
                                           pairlist)
    
    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord_pbc(:,:)
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    integer          :: num_pro_idr
    ! 
    integer          :: i, j, k, n, nloops
    integer          :: i_idr, j_idr
    integer          :: i_atom, j_atom
    integer          :: i_cell, j_cell
    integer          :: i_nbcell
    integer          :: i_chain_id, j_chain_id
    integer          :: pbc_int
    ! 
    integer          :: ini_excl, fin_excl
    integer          :: num_idr_max
    ! 
    real(wp)         :: dij(1:3)
    real(wp)         :: dij_pbc(1:3)
    real(wp)         :: rij_sqr
    real(wp)         :: pairdist_126
    real(wp)         :: pairdist_126_sqr
    ! 
    real(wp)         :: box_size(1:3), half_box_size(1:3)
    ! 
    integer, pointer :: cg_list_idr(:)
    integer, pointer :: icell_atom(:)
    integer, pointer :: cell_list_idr(:), cell_head_idr(:)
    integer, pointer :: num_idr_pre(:), num_idr(:)
    ! 
    logical          :: do_allocate
    logical          :: do_calc
    ! 
    integer          :: id, my_id, nthread
#ifdef OMP
    integer          :: omp_get_max_threads, omp_get_thread_num
#endif

    num_pro_idr        = enefunc%num_cg_particle_IDR_KH

    box_size(1)        = boundary%box_size_x
    box_size(2)        = boundary%box_size_y
    box_size(3)        = boundary%box_size_z
    half_box_size(1:3) = box_size(1:3) * 0.5_wp

    pairdist_126       = pairlist%cg_pairlistdist_126
    pairdist_126_sqr   = pairdist_126 * pairdist_126

    icell_atom         => pairlist%cell_index_cg_all
    cell_list_idr      => pairlist%cell_linked_list_cg_IDR_KH
    cell_head_idr      => pairlist%cell_head_index_cg_IDR_KH
    num_idr_pre        => pairlist%num_cg_IDR_KH_pre
    num_idr            => pairlist%num_cg_IDR_KH
    cg_list_idr        => enefunc%cg_particle_IDR_KH

#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif

    pairlist%num_cg_IDR_KH_calc(1:num_pro_idr,1:nthread) = 0

    if (pairlist%allocate_pbc_cg_IDR_KH) then
      nloops      = 2
      do_allocate = .true.
    else
      nloops      = 1
      do_allocate = .false.
    end if

    do n = 1, nloops

      num_idr(1:nthread) = 0

      if (.not. do_allocate) num_idr_pre(1:nthread) = 0

      !$omp parallel                        &
      !$omp private(id, my_id, i, j, k,     &
      !$omp         pbc_int,                &
      !$omp         i_idr, j_idr,           &
      !$omp         i_atom, j_atom,         &
      !$omp         i_cell, j_cell,         &
      !$omp         i_nbcell, do_calc,      &
      !$omp         ini_excl, fin_excl,     &
      !$omp         dij, dij_pbc, rij_sqr,  &
      !$omp         i_chain_id, j_chain_id) &
      !$omp firstprivate(do_allocate)
      !
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      my_id = my_city_rank * nthread + id
      id    = id + 1
      
      do i_idr = 1, num_pro_idr - 1

        if (mod(i_idr - 1, nproc_city * nthread) == my_id) then

          i_atom = cg_list_idr(i_idr)
          i_cell      = icell_atom(i_atom)
          i_chain_id  = enefunc%mol_chain_id(i_atom)

          do i_nbcell = 1, boundary%num_neighbor_cells_CG_126

            j_cell = boundary%neighbor_cells_CG_126(i_nbcell, i_cell)
            j_idr = cell_head_idr(j_cell)

            do while (j_idr > i_idr)

              j_atom = cg_list_idr(j_idr)

              !
              ! get j properties
              ! 
              j_chain_id  = enefunc%mol_chain_id(j_atom)
              
              do_calc = .true.
              ! 
              if (i_chain_id == j_chain_id .and. i_atom == j_atom - 1) then
                do_calc = .false.
              end if
              ! 
              if (.not. do_calc) then
                j_idr = cell_list_idr(j_idr)
                cycle
              end if

              ! 
              ! distance within pairlistdist?
              ! 
              dij(1:3)     = coord_pbc(1:3,j_atom) - coord_pbc(1:3,i_atom)
              call check_pbc(box_size, dij, pbc_int)

              rij_sqr      = dij(1) * dij(1) + dij(2) * dij(2) + dij(3) * dij(3)
              
              ! --------------------
              ! core: count and fill
              ! --------------------
              ! 
              if (rij_sqr < pairdist_126_sqr) then
                num_idr(id) = num_idr(id) + 1
                if (.not. do_allocate) then
                  pairlist%cg_IDR_KH_list(num_idr(id), id) = j_atom * 27 + pbc_int
                end if
              end if

              j_idr = cell_list_idr(j_idr)
            end do              ! j_atom

          end do                ! i_nbcell

        end if                  ! compute on the current thread

        if (.not. do_allocate) then
          pairlist%num_cg_IDR_KH_calc(i_idr,id) = num_idr(id) - num_idr_pre(id)
          num_idr_pre(id) = num_idr(id)
        end if

      end do                    ! icharge = 1 to num_charge
      !$omp end parallel

      ! allocate memory of pairlist
      !
      if (do_allocate) then
        num_idr_max = max(1, maxval(num_idr(1:nthread)))
        num_idr_max = int(real(num_idr_max,wp) * FactNumNb15)

        call alloc_pairlist(pairlist, PairListPbcCGIDRKH, num_idr_max)

        pairlist%num_cg_IDR_KH_max = num_idr_max

        do_allocate = .false.
      end if
      
    end do                      ! nloops

    ! check memory size
    !
    call check_pairlist_memory_size(num_idr, pairlist%num_cg_IDR_KH_max, &
        pairlist%allocate_pbc_cg_IDR_KH)

    return

  end subroutine update_pairlist_pbc_cg_IDR_KH

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_pbc_cg_kh
  !> @brief        update pairlist for CG KH protein-protein interactions
  !! @authors      CT
  !! @param[in]    enefunc   : potential energy functions information
  !! @param[in]    boundary  : boundary conditions information
  !! @param[in]    coord_pbc : PBC coordinates
  !! @param[inout] pairlist  : pairlist information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_pbc_cg_kh(enefunc, boundary, coord_pbc, pairlist)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord_pbc(:,:)
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    integer          :: num_pro_kh
    ! 
    integer          :: i, j, k, n, nloops
    integer          :: i_kh, j_kh
    integer          :: i_atom, j_atom
    integer          :: i_cell, j_cell
    integer          :: i_nbcell
    integer          :: i_chain_id, j_chain_id
    integer          :: pbc_int
    ! 
    integer          :: ini_excl, fin_excl
    integer          :: num_kh_max
    ! 
    real(wp)         :: dij(1:3)
    real(wp)         :: dij_pbc(1:3)
    real(wp)         :: rij_sqr
    real(wp)         :: pairdist_126
    real(wp)         :: pairdist_126_sqr
    ! 
    real(wp)         :: box_size(1:3), half_box_size(1:3)
    ! 
    integer, pointer :: cg_list_kh(:)
    integer, pointer :: icell_atom(:)
    integer, pointer :: cell_list_kh(:), cell_head_kh(:)
    integer, pointer :: num_kh_pre(:), num_kh(:)
    ! 
    logical          :: do_allocate
    logical          :: do_calc
    logical          :: is_idr_i
    logical          :: is_idr_j
    ! 
    integer          :: id, my_id, nthread
#ifdef OMP
    integer          :: omp_get_max_threads, omp_get_thread_num
#endif


    num_pro_kh         = enefunc%num_cg_particle_KH

    box_size(1)        = boundary%box_size_x
    box_size(2)        = boundary%box_size_y
    box_size(3)        = boundary%box_size_z
    half_box_size(1:3) = box_size(1:3) * 0.5_wp

    pairdist_126       = pairlist%cg_pairlistdist_126
    pairdist_126_sqr   = pairdist_126 * pairdist_126

    icell_atom         => pairlist%cell_index_cg_all
    cell_list_kh       => pairlist%cell_linked_list_cg_KH
    cell_head_kh       => pairlist%cell_head_index_cg_KH
    num_kh_pre         => pairlist%num_cg_KH_pre
    num_kh             => pairlist%num_cg_KH
    cg_list_kh         => enefunc%cg_particle_KH

#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif


    pairlist%num_cg_KH_calc(1:num_pro_kh,1:nthread) = 0

    if (pairlist%allocate_pbc_cg_KH) then
      nloops      = 2
      do_allocate = .true.
    else
      nloops      = 1
      do_allocate = .false.
    end if

    do n = 1, nloops

      num_kh(1:nthread) = 0

      if (.not. do_allocate) num_kh_pre(1:nthread) = 0

      !$omp parallel                        &
      !$omp private(id, my_id, i, j, k,     &
      !$omp         pbc_int,                &
      !$omp         i_kh, j_kh,             &
      !$omp         i_atom, j_atom,         &
      !$omp         i_cell, j_cell,         &
      !$omp         is_idr_i, is_idr_j,     &
      !$omp         i_nbcell, do_calc,      &
      !$omp         ini_excl, fin_excl,     &
      !$omp         dij, dij_pbc, rij_sqr,  &
      !$omp         i_chain_id, j_chain_id) &
      !$omp firstprivate(do_allocate)
      !
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      my_id = my_city_rank * nthread + id
      id    = id + 1
      
      do i_kh = 1, num_pro_kh - 1


        if (mod(i_kh - 1, nproc_city * nthread) == my_id) then

          i_atom = cg_list_kh(i_kh)
          i_cell = icell_atom(i_atom)

          !
          ! get i properties
          ! 
          i_chain_id  = enefunc%mol_chain_id(i_atom)
          if (enefunc%cg_IDR_HPS_is_IDR(i_atom) .or. enefunc%cg_IDR_KH_is_IDR(i_atom)) then
            is_idr_i = .true.
          else
            is_idr_i = .false.
          end if

          do i_nbcell = 1, boundary%num_neighbor_cells_CG_126

            j_cell = boundary%neighbor_cells_CG_126(i_nbcell, i_cell)
            j_kh   = cell_head_kh (j_cell)

            do while (j_kh  > i_kh)

              j_atom = cg_list_kh (j_kh)

              !
              ! get j properties
              ! 
              j_chain_id  = enefunc%mol_chain_id(j_atom)
              if (enefunc%cg_IDR_HPS_is_IDR(j_atom) .or. enefunc%cg_IDR_KH_is_IDR(j_atom)) then
                is_idr_j = .true.
              else
                is_idr_j = .false.
              end if

              ! 
              ! check inter-molecule kh interactions
              ! 
              if (enefunc%cg_kh_mol_pair(i_chain_id, j_chain_id) == 0) then
                j_kh = cell_list_kh(j_kh)
                cycle
              end if

              do_calc = .true.
              ! 
              if (is_idr_i .and. is_idr_j) then
                do_calc = .false.
              else
                ! 
                ! exclude interactions in exclusion list
                ! 
                ini_excl = enefunc%cg_istart_nonb_excl(i_atom)
                fin_excl = ini_excl + enefunc%num_nonb_excl(i_atom) - 1
                do k = ini_excl, fin_excl 
                  if (j_atom == enefunc%nonb_excl_list(k)) then
                    do_calc = .false.
                    exit
                  end if
                end do
                !
                if (j_chain_id == i_chain_id) then
                  if (.not. (is_idr_i .or. is_idr_j)) then
                    do_calc = .false.
                  end if
                end if
                ! 
              end if
              ! 
              if (.not. do_calc) then
                j_kh  = cell_list_kh (j_kh)
                cycle
              end if

              ! 
              ! distance within pairlistdist?
              ! 
              dij(1:3) = coord_pbc(1:3,j_atom) - coord_pbc(1:3,i_atom)
              call check_pbc(box_size, dij, pbc_int)

              rij_sqr  = dij(1) * dij(1) + dij(2) * dij(2) + dij(3) * dij(3)
              
              ! --------------------
              ! core: count and fill
              ! --------------------
              ! 
              if (rij_sqr < pairdist_126_sqr) then
                num_kh (id) = num_kh (id) + 1
                if (.not. do_allocate) then
                  pairlist%cg_KH_list(num_kh(id), id) = j_atom * 27 + pbc_int
                  pairlist%cg_KH_model_list(num_kh(id), id) = enefunc%cg_kh_mol_pair(i_chain_id, j_chain_id)
                end if
              end if

              j_kh = cell_list_kh(j_kh)
            end do              ! j_atom

          end do                ! i_nbcell

        end if                  ! compute on the current thread

        if (.not. do_allocate) then
          pairlist%num_cg_KH_calc(i_kh,id) = num_kh(id) - num_kh_pre(id)
          num_kh_pre(id) = num_kh(id)
        end if

      end do                    ! icharge = 1 to num_charge
      !$omp end parallel

      ! allocate memory of pairlist
      !
      if (do_allocate) then
        num_kh_max = max(1, maxval(num_kh(1:nthread)))
        num_kh_max = int(real(num_kh_max,wp) * FactNumNb15)

        call alloc_pairlist(pairlist, PairListPbcCGKH, num_kh_max)

        pairlist%num_cg_KH_max = num_kh_max

        do_allocate = .false.
      end if
      
    end do                      ! nloops

    ! check memory size
    !
    call check_pairlist_memory_size(num_kh, pairlist%num_cg_KH_max, &
        pairlist%allocate_pbc_cg_KH)

    return

  end subroutine update_pairlist_pbc_cg_kh

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_pbc_cg_pwmcos
  !> @brief        update pairlist for CG protein-DNA PWMcos interactions
  !! @authors      CT
  !! @param[in]    enefunc   : potential energy functions information
  !! @param[in]    boundary  : boundary conditions information
  !! @param[in]    coord_pbc : PBC coordinates
  !! @param[inout] pairlist  : pairlist information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_pbc_cg_pwmcos(enefunc, boundary, coord_pbc, &
                                           pairlist)
    
    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord_pbc(:,:)
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    integer          :: num_pwmcos
    integer          :: num_atom
    ! 
    integer          :: i, k, n, nloops
    integer          :: i_pwmcos, j_base
    integer          :: i_atom, j_atom
    integer          :: i_cell, j_cell
    integer          :: i_nbcell
    integer          :: i_chain_id, j_chain_id
    integer          :: pbc_int
    ! 
    integer          :: num_pwmcos_max
    ! 
    real(wp)         :: dij(1:3)
    real(wp)         :: dij_pbc(1:3)
    real(wp)         :: rij_sqr
    real(wp)         :: pairdist_PWMcos
    real(wp)         :: pairdist_PWMcos_sqr
    ! 
    real(wp)         :: box_size(1:3), half_box_size(1:3)
    ! 
    integer, pointer :: cg_list_pwmcos(:)
    integer, pointer :: cg_list_base(:)
    integer, pointer :: icell_atom(:)
    integer, pointer :: cell_list_DNA_base(:), cell_head_DNA_base(:)
    integer, pointer :: num_pairlist_pwmcos_pre(:), num_pairlist_pwmcos(:)
    ! 
    logical          :: do_allocate
    logical          :: do_calc
    ! 
    integer          :: id, my_id, nthread
#ifdef OMP
    integer          :: omp_get_max_threads, omp_get_thread_num
#endif


    num_pwmcos              = enefunc%num_pwmcos_terms
    num_atom                = size(coord_pbc(1,:))

    box_size(1)             = boundary%box_size_x
    box_size(2)             = boundary%box_size_y
    box_size(3)             = boundary%box_size_z
    half_box_size(1:3)      = box_size(1:3) * 0.5_wp

    pairdist_PWMcos         = pairlist%cg_pairlistdist_PWMcos
    pairdist_PWMcos_sqr     = pairdist_PWMcos * pairdist_PWMcos

    icell_atom              => pairlist%cell_index_cg_all
    cell_list_DNA_base      => pairlist%cell_linked_list_cg_DNA_base
    cell_head_DNA_base      => pairlist%cell_head_index_cg_DNA_base
    num_pairlist_pwmcos_pre => pairlist%num_cg_pwmcos_pre
    num_pairlist_pwmcos     => pairlist%num_cg_pwmcos
    cg_list_pwmcos          => enefunc%pwmcos_protein_id
    cg_list_base            => enefunc%cg_particle_DNA_base

#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif


    pairlist%num_cg_pwmcos_calc(1:num_pwmcos, 1:nthread) = 0

    if (pairlist%allocate_pbc_cg_pwmcos) then
      nloops      = 2
      do_allocate = .true.
    else
      nloops      = 1
      do_allocate = .false.
    end if

    do n = 1, nloops

      num_pairlist_pwmcos(1:nthread) = 0

      if (.not. do_allocate) num_pairlist_pwmcos_pre(1:nthread) = 0

      !$omp parallel                        &
      !$omp private(id, my_id, i, k,        &
      !$omp         pbc_int,                &
      !$omp         i_pwmcos, j_base,       &
      !$omp         i_atom, j_atom,         &
      !$omp         i_cell, j_cell,         &
      !$omp         i_nbcell, do_calc,      &
      !$omp         dij, dij_pbc, rij_sqr,  &
      !$omp         i_chain_id, j_chain_id) &
      !$omp firstprivate(do_allocate)
      !
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      my_id = my_city_rank * nthread + id
      id    = id + 1
      
      do i_pwmcos = 1, num_pwmcos

        if (mod(i_pwmcos - 1, nproc_city * nthread) == my_id) then

          i_atom     = cg_list_pwmcos(i_pwmcos)
          i_cell     = icell_atom(i_atom)
          i_chain_id = enefunc%mol_chain_id(i_atom)

          do i_nbcell = 1, boundary%num_neighbor_cells_CG_PWMcos

            j_cell = boundary%neighbor_cells_CG_PWMcos(i_nbcell, i_cell)
            j_base = cell_head_DNA_base(j_cell)

            do while (j_base > 0)

              j_atom = cg_list_base(j_base)

              j_chain_id  = enefunc%mol_chain_id(j_atom)

              if (enefunc%pwmcos_mol_pair(i_chain_id, j_chain_id) == 0) then
                j_base = cell_list_DNA_base(j_base)
                cycle
              end if

              ! 
              ! distance within pairlistdist?
              ! 
              dij(1:3)     = coord_pbc(1:3,j_atom) - coord_pbc(1:3,i_atom)
              call check_pbc(box_size, dij, pbc_int)

              rij_sqr      = dij(1) * dij(1) + dij(2) * dij(2) + dij(3) * dij(3)
              
              ! --------------------
              ! core: count and fill
              ! --------------------
              ! 
              if (rij_sqr < pairdist_PWMcos_sqr) then
                if (((j_atom - 3 >= 1) .and. &
                    (enefunc%mol_chain_id(j_atom - 3) == j_chain_id)) .and. &
                    ((j_atom + 3 <= num_atom) .and. &
                    (enefunc%mol_chain_id(j_atom + 3) == j_chain_id))) then
                  num_pairlist_pwmcos(id) = num_pairlist_pwmcos(id) + 1
                  if (.not. do_allocate) then
                    pairlist%cg_pwmcos_list(num_pairlist_pwmcos(id), id) = j_atom * 27 + pbc_int
                  end if
                end if
              end if

              j_base = cell_list_DNA_base(j_base)
            end do              ! j_base

          end do                ! i_nbcell

        end if                  ! compute on the current thread

        if (.not. do_allocate) then
          pairlist%num_cg_pwmcos_calc(i_pwmcos,id) = num_pairlist_pwmcos(id) - num_pairlist_pwmcos_pre(id)
          num_pairlist_pwmcos_pre(id) = num_pairlist_pwmcos(id)
        end if

      end do
      !$omp end parallel

      ! allocate memory of pairlist
      !
      if (do_allocate) then
        num_pwmcos_max = max(1, maxval(num_pairlist_pwmcos(1:nthread)))
        num_pwmcos_max = int(real(num_pwmcos_max,wp) * FactNumNb15)

        call alloc_pairlist(pairlist, PairListPbcCGPWMcos, num_pwmcos_max)

        pairlist%num_cg_pwmcos_max = num_pwmcos_max

        do_allocate = .false.
      end if
      
    end do                      ! nloops

    ! check memory size
    !
    call check_pairlist_memory_size(num_pairlist_pwmcos, &
         pairlist%num_cg_pwmcos_max, pairlist%allocate_pbc_cg_pwmcos)

    return

  end subroutine update_pairlist_pbc_cg_pwmcos

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_pbc_cg_pwmcosns
  !> @brief        update pairlist for CG protein-DNA PWMcosns interactions
  !! @authors      CT
  !! @param[in]    enefunc   : potential energy functions information
  !! @param[in]    boundary  : boundary conditions information
  !! @param[in]    coord_pbc : PBC coordinates
  !! @param[inout] pairlist  : pairlist information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_pbc_cg_pwmcosns(enefunc, boundary, coord_pbc, &
                                             pairlist)
    
    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord_pbc(:,:)
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    integer          :: num_pwmcosns
    integer          :: num_atom
    ! 
    integer          :: i, k, n, nloops
    integer          :: i_pwmcosns, j_phos
    integer          :: i_atom, j_atom
    integer          :: i_cell, j_cell
    integer          :: i_nbcell
    integer          :: i_chain_id, j_chain_id
    integer          :: pbc_int
    ! 
    integer          :: num_pwmcosns_max
    ! 
    real(wp)         :: dij(1:3)
    real(wp)         :: dij_pbc(1:3)
    real(wp)         :: rij_sqr
    real(wp)         :: pairdist_PWMcos
    real(wp)         :: pairdist_PWMcos_sqr
    ! 
    real(wp)         :: box_size(1:3), half_box_size(1:3)
    ! 
    integer, pointer :: cg_list_pwmcosns(:)
    integer, pointer :: cg_list_phos(:)
    integer, pointer :: icell_atom(:)
    integer, pointer :: cell_list_DNA_phos(:), cell_head_DNA_phos(:)
    integer, pointer :: num_pairlist_pwmcosns_pre(:), num_pairlist_pwmcosns(:)
    ! 
    logical          :: do_allocate
    logical          :: do_calc
    ! 
    integer          :: id, my_id, nthread
#ifdef OMP
    integer          :: omp_get_max_threads, omp_get_thread_num
#endif


    num_pwmcosns              = enefunc%num_pwmcosns_terms
    num_atom                = size(coord_pbc(1,:))

    box_size(1)             = boundary%box_size_x
    box_size(2)             = boundary%box_size_y
    box_size(3)             = boundary%box_size_z
    half_box_size(1:3)      = box_size(1:3) * 0.5_wp

    pairdist_PWMcos         = pairlist%cg_pairlistdist_PWMcos
    pairdist_PWMcos_sqr     = pairdist_PWMcos * pairdist_PWMcos

    icell_atom              => pairlist%cell_index_cg_all
    cell_list_DNA_phos      => pairlist%cell_linked_list_cg_DNA_phos
    cell_head_DNA_phos      => pairlist%cell_head_index_cg_DNA_phos
    num_pairlist_pwmcosns_pre => pairlist%num_cg_pwmcosns_pre
    num_pairlist_pwmcosns     => pairlist%num_cg_pwmcosns
    cg_list_pwmcosns          => enefunc%pwmcosns_protein_id
    cg_list_phos            => enefunc%cg_particle_DNA_phos

#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif


    pairlist%num_cg_pwmcosns_calc(1:num_pwmcosns, 1:nthread) = 0

    if (pairlist%allocate_pbc_cg_pwmcosns) then
      nloops      = 2
      do_allocate = .true.
    else
      nloops      = 1
      do_allocate = .false.
    end if

    do n = 1, nloops

      num_pairlist_pwmcosns(1:nthread) = 0

      if (.not. do_allocate) num_pairlist_pwmcosns_pre(1:nthread) = 0

      !$omp parallel                        &
      !$omp private(id, my_id, i, k,        &
      !$omp         pbc_int,                &
      !$omp         i_pwmcosns, j_phos,     &
      !$omp         i_atom, j_atom,         &
      !$omp         i_cell, j_cell,         &
      !$omp         i_nbcell, do_calc,      &
      !$omp         dij, dij_pbc, rij_sqr,  &
      !$omp         i_chain_id, j_chain_id) &
      !$omp firstprivate(do_allocate)
      !
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      my_id = my_city_rank * nthread + id
      id    = id + 1
      
      do i_pwmcosns = 1, num_pwmcosns

        if (mod(i_pwmcosns - 1, nproc_city * nthread) == my_id) then

          i_atom     = cg_list_pwmcosns(i_pwmcosns)
          i_cell     = icell_atom(i_atom)
          i_chain_id = enefunc%mol_chain_id(i_atom)

          do i_nbcell = 1, boundary%num_neighbor_cells_CG_PWMcos

            j_cell = boundary%neighbor_cells_CG_PWMcos(i_nbcell, i_cell)
            j_phos = cell_head_DNA_phos(j_cell)

            do while (j_phos > 0)

              j_atom = cg_list_phos(j_phos)

              j_chain_id  = enefunc%mol_chain_id(j_atom)

              if (enefunc%pwmcosns_mol_pair(i_chain_id, j_chain_id) == 0) then
                j_phos = cell_list_DNA_phos(j_phos)
                cycle
              end if

              ! 
              ! distance within pairlistdist?
              ! 
              dij(1:3) = coord_pbc(1:3,j_atom) - coord_pbc(1:3,i_atom)
              call check_pbc(box_size, dij, pbc_int)

              rij_sqr  = dij(1) * dij(1) + dij(2) * dij(2) + dij(3) * dij(3)

              ! --------------------
              ! core: count and fill
              ! --------------------
              ! 
              if (rij_sqr < pairdist_PWMcos_sqr) then
                num_pairlist_pwmcosns(id) = num_pairlist_pwmcosns(id) + 1
                if (.not. do_allocate) then
                  pairlist%cg_pwmcosns_list(num_pairlist_pwmcosns(id), id) = j_atom * 27 + pbc_int
                end if
              end if

              j_phos = cell_list_DNA_phos(j_phos)
            end do              ! j_phos

          end do                ! i_nbcell

        end if                  ! compute on the current thread

        if (.not. do_allocate) then
          pairlist%num_cg_pwmcosns_calc(i_pwmcosns,id) = num_pairlist_pwmcosns(id) - num_pairlist_pwmcosns_pre(id)
          num_pairlist_pwmcosns_pre(id) = num_pairlist_pwmcosns(id)
        end if

      end do
      !$omp end parallel

      ! allocate memory of pairlist
      !
      if (do_allocate) then
        num_pwmcosns_max = max(1, maxval(num_pairlist_pwmcosns(1:nthread)))
        num_pwmcosns_max = int(real(num_pwmcosns_max,wp) * FactNumNb15)

        call alloc_pairlist(pairlist, PairListPbcCGPWMcosns, num_pwmcosns_max)

        pairlist%num_cg_pwmcosns_max = num_pwmcosns_max

        do_allocate = .false.
      end if
      
    end do                      ! nloops

    ! check memory size
    !
    call check_pairlist_memory_size(num_pairlist_pwmcosns, &
         pairlist%num_cg_pwmcosns_max, pairlist%allocate_pbc_cg_pwmcosns)

    return

  end subroutine update_pairlist_pbc_cg_pwmcosns

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_pairlist_memory_size
  !> @brief        check memory size of pairlist
  !! @authors      TM
  !! @param[in]    num_nb15       : number of 1-5 nonbonded pairs
  !! @param[in]    num_nb15_limit : already allocated size
  !! @param[out]   alloc_list     : allocate list or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_pairlist_memory_size(num_nb15_calc, num_nb15_limit, &
                                        alloc_list)

    ! formal arguments
    integer,          intent(in)    :: num_nb15_calc(:)
    integer,          intent(in)    :: num_nb15_limit
    logical,          intent(out)   :: alloc_list

    ! local variables
    integer                  :: threshold
    integer                  :: num_nb15_max

    num_nb15_max = 0
    num_nb15_max = max(num_nb15_max, maxval(num_nb15_calc(:)))

    ! check memory size
    !   if the num_nb15_max is larger than FactThreshold*(allocated size),
    !   allocate memory of pairlist at the next update
    !
    threshold = int(real(num_nb15_limit,wp)*FactThreshold)

    if (num_nb15_max >= threshold) then
      alloc_list = .true.
    else
      alloc_list = .false.
    end if

    return

  end subroutine check_pairlist_memory_size

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_pbc
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_pbc(box_size, dij, pbc_int)

    ! formal arguments
    real(wp),         intent(in)    :: box_size(:)
    real(wp),         intent(inout) :: dij(:)
    integer,          intent(inout) :: pbc_int

    ! local variables
    integer                  :: i, j, k


    if (dij(1) > box_size(1)/2.0_dp) then
      i = 0
      dij(1) = dij(1) - box_size(1)
    else if (dij(1) < -box_size(1)/2.0_dp) then
      i = 2
      dij(1) = dij(1) + box_size(1)
    else
      i = 1
    end if

    if (dij(2) > box_size(2)/2.0_dp) then
      j = 0
      dij(2) = dij(2) - box_size(2)
    else if (dij(2) < -box_size(2)/2.0_dp) then
      j = 2
      dij(2) = dij(2) + box_size(2)
    else
      j = 1
    end if

    if (dij(3) > box_size(3)/2.0_dp) then
      k = 0
      dij(3) = dij(3) - box_size(3)
    else if (dij(3) < -box_size(3)/2.0_dp) then
      k = 2
      dij(3) = dij(3) + box_size(3)
    else
      k = 1
    end if

    pbc_int = i + j*3 + k*9

    return

  end subroutine check_pbc

end module at_pairlist_mod
