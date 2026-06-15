!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_pairlist_mod
!> @brief   Sets pairlists for nonbonded interactions.
!! @authors Jaewoon Jung (JJ), Kiyotaka Sakamoto (KS)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_pairlist_mod

  use cg_pairlist_str_mod
  use cg_boundary_str_mod
  use cg_enefunc_str_mod
  use cg_domain_str_mod
  use cg_energy_str_mod
  use molecules_str_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif
#ifdef OMP
    use omp_lib, only: omp_get_thread_num
#endif

  implicit none
  private

  ! subroutines
  public  :: setup_pairlist
  public  :: update_pairlist_cg_alloc
  public  :: update_pairlist_cg
  public  :: update_pairlist_martini
  public  :: update_pairlist_cg_pwmcos_alloc
  public  :: update_pairlist_cg_pwmcos
  public  :: update_pairlist_cg_pwmcosns_alloc
  public  :: update_pairlist_cg_pwmcosns

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_pairlist
  !> @brief        Initializes, allocates, and sets up pairlists for nonbonded interactions.
  !! @authors      JJ
  !! @param[in]    boundary : boundary condition information
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_pairlist(boundary, enefunc, domain, pairlist)

    ! formal arguments
    type(s_boundary),         intent(in)    :: boundary
    type(s_enefunc),          intent(inout) :: enefunc 
    type(s_domain),   target, intent(inout) :: domain 
    type(s_pairlist),         intent(inout) :: pairlist

    integer                   :: ncel, i, alloc_num
    real(wip),        pointer :: coord(:,:)
    real(wp),         pointer :: trans(:,:), coord_pbc(:,:)

    coord                => domain%coord
    trans                => domain%trans_vec
    coord_pbc            => domain%translated



    pairlist%cg_pairlistdist_vdw    = enefunc%cg_pairlistdist_vdw
    pairlist%cg_pairlistdist_ele    = enefunc%cg_pairlistdist_ele
    pairlist%cg_pairlistdist_126    = enefunc%cg_pairlistdist_126
    pairlist%cg_pairlistdist_WF     = enefunc%cg_pairlistdist_WF
    pairlist%cg_pairlistdist_PWMcos = enefunc%cg_pairlistdist_PWMcos
    pairlist%cg_pairlistdist_DNAbp  = enefunc%cg_pairlistdist_DNAbp
    pairlist%cg_pairlistdist_exv    = enefunc%cg_pairlistdist_exv

    ncel      = domain%num_cell_local + domain%num_cell_boundary

    if (enefunc%forcefield == ForcefieldRESIDCG) then
      alloc_num = MaxAtom_domain
      call alloc_pairlist(pairlist, PairListCGExvNum,       alloc_num)
      if (enefunc%cg_ele_calc) then
        alloc_num = Max_cg_Elec
        call alloc_pairlist(pairlist, PairListCGEleNum,     alloc_num)
      end if
      if (enefunc%cg_DNA_exv_calc) then
        alloc_num = Max_cg_DNA
        call alloc_pairlist(pairlist, PairListCGDNAExvNum,  alloc_num)
      end if
      if (enefunc%cg_DNA_base_pair_calc) then
        alloc_num = Max_cg_Base
        call alloc_pairlist(pairlist, PairListCGDNABaseNum, alloc_num)
      end if
      if (enefunc%cg_KH_calc) then
        alloc_num = Max_cg_KH
        call alloc_pairlist(pairlist, PairListCGKHNum,      alloc_num)
      end if
      if (enefunc%cg_IDR_KH_calc) then
        alloc_num = Max_cg_IDR_KH
        call alloc_pairlist(pairlist, PairListCGIDRKHNum,   alloc_num)
      end if
      if (enefunc%cg_IDR_HPS_calc) then
        alloc_num = Max_cg_IDR_HPS
        call alloc_pairlist(pairlist, PairListCGIDRHPSNum,  alloc_num)
      end if
      if (enefunc%cg_IDR_MPIPI_calc) then
        alloc_num = Max_cg_IDR_MPIPI
        call alloc_pairlist(pairlist, PairListCGIDRMPIPINum,  alloc_num)
      end if
      if (enefunc%cg_pwmcos_calc) then
        alloc_num = MaxPwmCos
        call alloc_pairlist(pairlist, PairListPWMCOSNum,    alloc_num)
      end if
      if (enefunc%cg_pwmcosns_calc) then
        alloc_num = MaxPwmCosns
        call alloc_pairlist(pairlist, PairListPWMCOSnsNum,  alloc_num)
      end if
    else if (enefunc%forcefield == ForcefieldGROMARTINI) then
      alloc_num = MaxAtom_domain
      call alloc_pairlist(pairlist, PairListCGExvNum,       alloc_num)
      alloc_num = Max_cg_Elec
      call alloc_pairlist(pairlist, PairListCGEleNum,       alloc_num)
    end if

    call timer(TimerPairList, TimerOn)

    select case (boundary%type)

    case (BoundaryTypeNOBC)

      !$omp parallel do private(i)
      do i = 1, domain%num_atom_domain + domain%num_atom_boundary
        coord_pbc(i,1) = real(coord(i,1),wp)
        coord_pbc(i,2) = real(coord(i,2),wp)
        coord_pbc(i,3) = real(coord(i,3),wp)
      end do
      !$omp end parallel do

    case (boundaryTypePBC)

      !$omp parallel do private(i)
      do i = 1, domain%num_atom_domain + domain%num_atom_boundary
        coord_pbc(i,1) = real(coord(i,1),wp) + trans(i,1)
        coord_pbc(i,2) = real(coord(i,2),wp) + trans(i,2)
        coord_pbc(i,3) = real(coord(i,3),wp) + trans(i,3)
      end do
      !$omp end parallel do

    end select

    if (enefunc%forcefield == ForcefieldRESIDCG) then
      call update_pairlist_cg_alloc(coord_pbc, enefunc, domain, pairlist)
      call update_pairlist_cg(coord_pbc, enefunc, domain, pairlist)
      if (enefunc%cg_pwmcos_calc) then
        call update_pairlist_cg_pwmcos_alloc(coord_pbc, enefunc, domain, &
                                       pairlist)
        call update_pairlist_cg_pwmcos(coord_pbc, enefunc, domain, &
                                       pairlist)
      end if
      if (enefunc%cg_pwmcosns_calc) then
        call update_pairlist_cg_pwmcosns_alloc(coord_pbc, enefunc, domain, &
                                         pairlist)
        call update_pairlist_cg_pwmcosns(coord_pbc, enefunc, domain, &
                                         pairlist)
      end if
    else
      call update_pairlist_Martini_alloc(coord_pbc, enefunc, domain, pairlist)
      call update_pairlist_Martini(coord_pbc, enefunc, domain, pairlist)
    end if

    call timer(TimerPairList, TimerOff)

    return

  end subroutine setup_pairlist

 !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_cg_alloc
  !> @brief        Updates pairlists for AICG.
  !! @authors      JJ, CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_cg_alloc(coord, enefunc, domain, pairlist)

    ! formal arguments
    real(wp),                 intent(in)    :: coord(:,:)
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(in)    :: domain
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: pairdist2_ele, pairdist2_126
    real(wp)                  :: pairdist2_WF
    real(wp)                  :: pairdist2_DNAbp, pairdist2_exv
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: rtmp(1:3)

    integer                   :: i, j, ij, k, ix, ixx, iy, iyy, inbc
    integer                   :: k_dna, k_base, k_kh
    integer                   :: k_idr_kh, k_idr_hps, k_idr_mpipi
    integer                   :: chain_idx, chain_idy
    integer                   :: i_atom, j_atom, num_atom
    integer                   :: start_i, start_j
    integer                   :: i_base_type, j_base_type
    logical                   :: is_nonlocal, is_WC_bp
    integer                   :: num_nb15, num_nb15_dna, num_nb15_base
    integer                   :: num_nb15_kh, num_nb15_idr_kh
    integer                   :: num_nb15_idr_hps
    integer                   :: num_nb15_idr_mpipi
    integer                   :: i_dna, j_dna, ij_dna
    integer                   :: i_kh, j_kh, ij_kh
    integer                   :: i_idr_kh, j_idr_kh, ij_idr_kh
    integer                   :: i_idr_hps, j_idr_hps, ij_idr_hps
    integer                   :: i_idr_mpipi, j_idr_mpipi, ij_idr_mpipi
    integer                   :: i_idr, j_idr, ij_idr, ij_idr1
    integer                   :: i_126
    integer                   :: id
    integer                   :: ncell, nboundary
    logical                   :: nb15_calc

    integer,          pointer :: natom(:), ncharge(:)
    integer,          pointer :: start_atom(:)
    integer,          pointer :: near_cells_count(:), far_cells_count(:)
    integer,          pointer :: near_cells(:,:), far_cells(:,:)
    integer,          pointer :: id_l2g(:)
    integer(1),       pointer :: dna_check(:)
    integer,          pointer :: base_type(:)
    integer(1),       pointer :: cg_pro_use_KH(:), cg_IDR_KH(:), cg_IDR_HPS(:)
    integer(1),       pointer :: cg_IDR_MPIPI(:)
    integer,          pointer :: cg_elec_list_inv(:)
    integer,          pointer :: cg_DNA_list_inv(:)
    integer,          pointer :: cg_base_list_inv(:)
    integer,          pointer :: cg_kh_list_inv(:)
    integer,          pointer :: cg_idr_kh_list_inv(:)
    integer,          pointer :: cg_idr_hps_list_inv(:)
    integer,          pointer :: cg_idr_mpipi_list_inv(:)
    integer,          pointer :: num_cg_exv_calc(:)
    integer,          pointer :: num_cg_DNA_exv_calc(:)
    integer,          pointer :: num_cg_DNA_base_calc(:)
    integer,          pointer :: num_cg_ele_calc(:)
    integer,          pointer :: num_cg_kh_calc(:)
    integer,          pointer :: num_cg_idr_kh_calc(:)
    integer,          pointer :: num_cg_idr_hps_calc(:)
    integer,          pointer :: num_cg_idr_mpipi_calc(:)
    integer(1),       pointer :: exclusion_mask(:,:)
    integer,          pointer :: chain_id(:)
    integer(1),       pointer :: cg_ele_mol_pair(:,:), cg_kh_mol_pair(:,:)
    logical,          pointer :: base_pair_is_WC(:,:)

    ncell                =  domain%num_cell_local
    nboundary            =  domain%num_cell_boundary
    natom                => domain%num_atom
    start_atom           => domain%start_atom
    ncharge              => domain%num_charge
    near_cells_count     => domain%near_cells_count
    far_cells_count      => domain%far_cells_count
    near_cells           => domain%near_cells
    far_cells            => domain%far_cells
    dna_check            => domain%dna_check
    cg_pro_use_KH        => domain%cg_pro_use_KH
    cg_IDR_KH            => domain%cg_IDR_KH
    cg_IDR_HPS           => domain%cg_IDR_HPS
    cg_IDR_MPIPI         => domain%cg_IDR_MPIPI
    id_l2g               => domain%id_l2g
    chain_id             => domain%mol_chain_id
    base_type            => domain%NA_base_type

    exclusion_mask       => enefunc%exclusion_mask
    cg_ele_mol_pair      => enefunc%cg_ele_mol_pair
    cg_kh_mol_pair       => enefunc%cg_kh_mol_pair
    base_pair_is_WC      => enefunc%base_pair_is_WC
    cg_elec_list_inv     => enefunc%cg_elec_list_inv
    cg_dna_list_inv      => enefunc%cg_dna_list_inv
    cg_base_list_inv     => enefunc%cg_base_list_inv
    cg_kh_list_inv       => enefunc%cg_kh_list_inv
    cg_idr_kh_list_inv   => enefunc%cg_idr_kh_list_inv
    cg_idr_hps_list_inv  => enefunc%cg_idr_hps_list_inv
    cg_idr_mpipi_list_inv => enefunc%cg_idr_mpipi_list_inv

    num_cg_exv_calc      => pairlist%num_cg_exv_calc
    num_cg_DNA_exv_calc  => pairlist%num_cg_DNA_exv_calc
    num_cg_ele_calc      => pairlist%num_cg_ele_calc
    num_cg_DNA_base_calc => pairlist%num_cg_DNA_base_calc
    num_cg_kh_calc       => pairlist%num_cg_kh_calc
    num_cg_idr_kh_calc   => pairlist%num_cg_idr_kh_calc
    num_cg_idr_hps_calc  => pairlist%num_cg_idr_hps_calc
    num_cg_idr_mpipi_calc => pairlist%num_cg_idr_mpipi_calc

    pairdist2_ele        =  pairlist%cg_pairlistdist_ele**2
    pairdist2_126        =  pairlist%cg_pairlistdist_126**2
    pairdist2_WF         =  pairlist%cg_pairlistdist_WF**2
    pairdist2_DNAbp      =  pairlist%cg_pairlistdist_DNAbp**2
    pairdist2_exv        =  pairlist%cg_pairlistdist_exv**2

    num_atom = domain%num_atom_domain + domain%num_atom_boundary
    num_cg_exv_calc(1:num_atom) = 0
    if (enefunc%cg_DNA_exv_calc) then
      num_atom = enefunc%num_cg_dna
      num_cg_DNA_exv_calc(1:num_atom) = 0
    end if
    if (enefunc%cg_DNA_base_pair_calc) then
      num_atom = enefunc%num_cg_base
      num_cg_DNA_base_calc(1:num_atom) = 0
    end if
    if (enefunc%cg_ele_calc) then
      num_atom = enefunc%num_cg_elec
      num_cg_ele_calc(1:num_atom) = 0
    end if
    if (enefunc%cg_IDR_KH_calc) then
      num_atom = enefunc%num_cg_IDR_KH
      num_cg_idr_kh_calc(1:num_atom) = 0
    end if
    if (enefunc%cg_KH_calc) then
      num_atom = enefunc%num_cg_KH
      num_cg_kh_calc(1:num_atom) = 0
    end if
    if (enefunc%cg_IDR_HPS_calc) then
      num_atom = enefunc%num_cg_IDR_HPS
      num_cg_idr_hps_calc(1:num_atom) = 0
    end if
    if (enefunc%cg_IDR_MPIPI_calc) then
      num_atom = enefunc%num_cg_IDR_MPIPI
      num_cg_idr_mpipi_calc(1:num_atom) = 0
    end if

    !$omp parallel default(shared)                                        &
    !$omp private(id, i, ix, iy, k, ij, j, inbc, num_nb15, num_nb15_dna,  &
    !$omp         num_nb15_base, num_nb15_kh, num_nb15_idr_kh,            &
    !$omp         num_nb15_idr_hps, i_dna, j_dna, ij_dna, i_126,          &
    !$omp         num_nb15_idr_mpipi,                                     &
    !$omp         rtmp, dij, rij2, chain_idx, chain_idy, nb15_calc,       &
    !$omp         i_atom, j_atom, i_base_type, j_base_type, is_nonlocal,  &
    !$omp         is_WC_bp, start_i, start_j, ixx, iyy, i_kh, j_kh,       &
    !$omp         ij_kh, i_idr_kh, j_idr_kh, ij_idr_kh, i_idr_hps,        &
    !$omp         j_idr_hps, ij_idr_hps, i_idr, j_idr, ij_idr, ij_idr1,   &
    !$omp         j_idr_mpipi, ij_idr_mpipi, i_idr_mpipi, k_idr_mpipi,    &
    !$omp         k_base, k_dna, k_kh, k_idr_kh, k_idr_hps)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    ! make a pairlist in the same cell
    !
    do i = id+1, ncell+nboundary, nthread

      start_i = start_atom(i)

      do ix = 1, natom(i) 

        ixx = ix + start_i
        rtmp(1:3) = coord(ixx,1:3)

        num_nb15_dna = 0
        num_nb15_base = 0
        num_nb15 = 0
        num_nb15_kh = 0
        num_nb15_idr_kh = 0
        num_nb15_idr_hps = 0
        num_nb15_idr_mpipi = 0
        i_dna     = dna_check(ixx)
        i_kh      = cg_pro_use_KH(ixx)
        i_idr_kh  = cg_IDR_KH(ixx)
        i_idr_hps = cg_IDR_HPS(ixx)
        i_idr_mpipi = cg_IDR_MPIPI(ixx)
        i_idr     = i_idr_kh + i_idr_hps + i_idr_mpipi
        i_126     = i_idr + i_kh

        if (i_dna /= 0) then
          chain_idx = chain_id(ixx)
          i_base_type = base_type(ixx)
          i_atom = id_l2g(ixx)
        else
          chain_idx = chain_id(ixx)
          i_atom = id_l2g(ixx)
        end if

        if (i <= ncell) then

          do iy = ix + 1, natom(i)

            iyy = iy + start_i

            if (exclusion_mask(iyy,ixx) /= 1) then

              dij(1:3) = rtmp(1:3) - coord(iyy,1:3)
              rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

              ! case 1: rij < pairlistdist_exv : write all interactions
              !
              if (rij2 < pairdist2_exv) then

                call pairlist_exv_DNA_idr_check(dna_check, cg_pro_use_KH,     &
                                cg_IDR_KH, cg_IDR_HPS, cg_IDR_MPIPI,          &
                                chain_id, base_type,                          &
                                base_pair_is_WC, id_l2g, cg_kh_mol_pair,      &
                                iyy, i_atom, chain_idx, i_idr, i_dna,         &
                                i_base_type, i_kh, i_idr_kh, i_idr_hps,       &
                                i_idr_mpipi,                                  &
                                num_nb15_base, num_nb15_dna, num_nb15,        &
                                num_nb15_kh, num_nb15_idr_kh,                 &
                                num_nb15_idr_hps, num_nb15_idr_mpipi)

              ! case 1: pairlistdist_exv < rij < pairlistdist_DNAbp 
              !         interaction of Base pair and IDR
              !
              else if ((i_dna>0.or.i_126>0) .and. rij2 < pairdist2_DNAbp) then

                call pairlist_DNA_idr_check(dna_check, cg_pro_use_KH,         &
                                cg_IDR_KH, cg_IDR_HPS, cg_IDR_MPIPI,          &
                                chain_id, base_type,                          &
                                base_pair_is_WC, id_l2g, cg_kh_mol_pair,      &
                                iyy, i_atom, chain_idx, i_idr, i_dna,         &
                                i_base_type, i_kh, i_idr_kh, i_idr_hps,       &
                                i_idr_mpipi,                                  &
                                num_nb15_base, num_nb15_kh, num_nb15_idr_kh,  &
                                num_nb15_idr_hps, num_nb15_idr_mpipi)

              ! case 3: pairlist2_DNAbp < rij < pairlist2_126
              !
              else if (i_126 > 0 .and. rij2 < pairdist2_126) then

                call pairlist_idr_check(cg_pro_use_KH, cg_IDR_KH, cg_IDR_HPS, &
                                cg_IDR_MPIPI,                                 &
                                chain_id, id_l2g, cg_kh_mol_pair, iyy,        &
                                i_atom, chain_idx, i_idr, i_kh, i_idr_kh,     &
                                i_idr_hps, i_idr_mpipi,                       &
                                num_nb15_kh, num_nb15_idr_kh,                 &
                                num_nb15_idr_hps, num_nb15_idr_mpipi)

              end if

            end if
          end do

        end if

        do inbc = 1, near_cells_count(i)

          j = near_cells(inbc,i)
          start_j = start_atom(j)

          do iy = 1, natom(j)

            iyy = iy + start_j

            if (exclusion_mask(iyy,ixx) /= 1) then

              dij(1:3) = rtmp(1:3) - coord(iyy,1:3)
              rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

              ! case 1: rij < pairlistdist_exv : write all interactions
              !
              if (rij2 < pairdist2_exv) then

                call pairlist_exv_DNA_idr_check(dna_check, cg_pro_use_KH,     &
                                cg_IDR_KH, cg_IDR_HPS, cg_IDR_MPIPI,          &
                                chain_id, base_type,                          &
                                base_pair_is_WC, id_l2g, cg_kh_mol_pair,      &
                                iyy, i_atom, chain_idx, i_idr, i_dna,         &
                                i_base_type, i_kh, i_idr_kh, i_idr_hps,       &
                                i_idr_mpipi,                                  &
                                num_nb15_base, num_nb15_dna, num_nb15,        &
                                num_nb15_kh, num_nb15_idr_kh,                 &
                                num_nb15_idr_hps, num_nb15_idr_mpipi)

              ! case 1: pairlistdist_exv < rij < pairlistdist_DNAbp
              !         interaction of Base pair and IDR
              !
              else if ((i_dna>0.or.i_126>0) .and. rij2 < pairdist2_DNAbp) then

                call pairlist_DNA_idr_check(dna_check, cg_pro_use_KH,         &
                                cg_IDR_KH, cg_IDR_HPS, cg_IDR_MPIPI,          &
                                chain_id, base_type,                          &
                                base_pair_is_WC, id_l2g, cg_kh_mol_pair,      &
                                iyy, i_atom, chain_idx, i_idr, i_dna,         &
                                i_base_type, i_kh, i_idr_kh, i_idr_hps,       &
                                i_idr_mpipi,                                  &
                                num_nb15_base, num_nb15_kh, num_nb15_idr_kh,  &
                                num_nb15_idr_hps, num_nb15_idr_mpipi)

              ! case 3: pairlist2_DNAbp < rij < pairlist2_126
              !
              else if (i_126 > 0 .and. rij2 < pairdist2_126) then

                call pairlist_idr_check(cg_pro_use_KH, cg_IDR_KH, cg_IDR_HPS, &
                                cg_IDR_MPIPI,                                 &
                                chain_id, id_l2g, cg_kh_mol_pair, iyy,        &
                                i_atom, chain_idx, i_idr, i_kh, i_idr_kh,     &
                                i_idr_hps, i_idr_mpipi,                       &
                                num_nb15_kh, num_nb15_idr_kh,                 &
                                num_nb15_idr_hps, num_nb15_idr_mpipi)

              end if

            end if

          end do

        end do

        if (i_126 > 0) then

          do inbc = 1, far_cells_count(i)

            j = far_cells(inbc,i) 
            start_j = start_atom(j)

            do iy = 1, natom(j)

              iyy = iy + start_j
              dij(1:3) = rtmp(1:3) - coord(iyy,1:3)
              rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

              if (rij2 < pairdist2_126) then

                call pairlist_idr_check(cg_pro_use_KH, cg_IDR_KH, cg_IDR_HPS, &
                                cg_IDR_MPIPI,                                 &
                                chain_id, id_l2g, cg_kh_mol_pair, iyy,        &
                                i_atom, chain_idx, i_idr, i_kh, i_idr_kh,     &
                                i_idr_hps, i_idr_mpipi,                       &
                                num_nb15_kh, num_nb15_idr_kh,                 &
                                num_nb15_idr_hps, num_nb15_idr_mpipi)

              end if
            end do
          end do
        end if

        if (num_nb15 > 0) then
          num_cg_exv_calc(ixx) = num_nb15
        end if
        if (num_nb15_dna > 0) then
          k_dna = cg_DNA_list_inv(ixx)
          num_cg_DNA_exv_calc(k_dna) = num_nb15_dna
        end if
        if (num_nb15_base > 0) then
          k_base = cg_base_list_inv(ixx)
          num_cg_DNA_base_calc(k_base) = num_nb15_base
        end if
        if (num_nb15_kh > 0) then
          k_kh = cg_kh_list_inv(ixx)
          num_cg_kh_calc(k_kh) = num_nb15_kh
        end if
        if (num_nb15_idr_kh > 0) then
          k_idr_kh = cg_idr_kh_list_inv(ixx)
          num_cg_idr_kh_calc(k_idr_kh) = num_nb15_idr_kh
        end if
        if (num_nb15_idr_hps > 0) then
          k_idr_hps = cg_idr_hps_list_inv(ixx)
          pairlist%num_cg_idr_hps_calc(k_idr_hps) = num_nb15_idr_hps
        end if
        if (num_nb15_idr_mpipi > 0) then
          k_idr_mpipi = cg_idr_mpipi_list_inv(ixx)
          pairlist%num_cg_idr_mpipi_calc(k_idr_mpipi) = num_nb15_idr_mpipi
        end if

      end do

    end do

    ! electrostatic
    !
    if (enefunc%cg_ele_calc) then

      do i = id+1, ncell+nboundary, nthread

        start_i = start_atom(i)

        do ix = 1, ncharge(i) 

          num_nb15 = 0
          ixx = ix + start_i
          chain_idx = chain_id(ixx)
          rtmp(1:3) = coord(ixx,1:3)

          if (i <= ncell) then

            do iy = ix + 1, ncharge(i)

              iyy = iy + start_i

              if (exclusion_mask(iyy,ixx) /= 1) then

                chain_idy = chain_id(iyy)
                if (cg_ele_mol_pair(chain_idx,chain_idy) /= 0) then
                  dij(1:3) = rtmp(1:3) - coord(iyy,1:3)
                  rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
                  if (rij2 < pairdist2_ele) num_nb15 = num_nb15 + 1
                end if

              end if
            end do

          end if

          do inbc = 1, near_cells_count(i)

            j = near_cells(inbc,i)
            start_j = start_atom(j)

            do iy = 1, ncharge(j)

              iyy = iy + start_j

              if (exclusion_mask(iyy,ixx) /= 1) then

                chain_idy = chain_id(iyy)
                if (cg_ele_mol_pair(chain_idx,chain_idy) /= 0) then
                  dij(1:3) = rtmp(1:3) - coord(iyy,1:3)
                  rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
                  if (rij2 < pairdist2_ele) num_nb15 = num_nb15 + 1
                end if

              end if

            end do

          end do

          do inbc = 1, far_cells_count(i)

            j = far_cells(inbc,i)
            start_j = start_atom(j)

            do iy = 1, ncharge(j)

              iyy = iy + start_j
              chain_idy = chain_id(iyy)
              if (cg_ele_mol_pair(chain_idx,chain_idy) /= 0) then
                dij(1:3) = rtmp(1:3) - coord(iyy,1:3)
                rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
                if (rij2 < pairdist2_ele) num_nb15 = num_nb15 + 1
              end if

            end do

          end do

          k = cg_elec_list_inv(ixx)
          num_cg_ele_calc(k) = num_nb15

        end do
      end do

    end if

    !$omp end parallel

    num_atom = domain%num_atom_domain + domain%num_atom_boundary
    k = 0
    do i = 1, num_atom
      k = max(k,num_cg_exv_calc(i))
    end do
#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, k, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
#endif
    Max_exv_nb15 = max(k*2, min_atom)
    call alloc_pairlist(pairlist, PairListCGExvList, &
                        Max_exv_nb15, MaxAtom_domain)

    if (enefunc%cg_DNA_exv_calc) then
      num_atom = enefunc%num_cg_dna
      k = 0
      do i = 1, num_atom
        k = max(k,num_cg_DNA_exv_calc(i))
      end do
#ifdef HAVE_MPI_GENESIS
      call mpi_allreduce(mpi_in_place, k, 1, mpi_integer, &
                         mpi_max, mpi_comm_country, ierror)
#endif
      Max_dna_exv_nb15 = max(k*2, min_atom)
      call alloc_pairlist(pairlist, PairListCGDNAExvList, &
                          Max_dna_exv_nb15, Max_cg_DNA)
    else
      call alloc_pairlist(pairlist, PairListCGDNAExvList, &
                          1, 1)
    end if

    if (enefunc%cg_DNA_base_pair_calc) then
      num_atom = enefunc%num_cg_base
      k = 0
      do i = 1, num_atom
        k = max(k,num_cg_DNA_base_calc(i))
      end do
#ifdef HAVE_MPI_GENESIS
      call mpi_allreduce(mpi_in_place, k, 1, mpi_integer, &
                         mpi_max, mpi_comm_country, ierror)
#endif
      Max_dna_base_nb15 = max(k*2, min_atom)
      call alloc_pairlist(pairlist, PairListCGDNABaseList, &
                          Max_dna_base_nb15, Max_cg_Base)
    else
      call alloc_pairlist(pairlist, PairListCGDNABaseList, &
                          1, 1)
    end if

    if (enefunc%cg_ele_calc) then
      num_atom = enefunc%num_cg_elec
      k = 0
      do i = 1, num_atom
        k = max(k,num_cg_ele_calc(i))
      end do
#ifdef HAVE_MPI_GENESIS
      call mpi_allreduce(mpi_in_place, k, 1, mpi_integer, &
                         mpi_max, mpi_comm_country, ierror)
#endif
      Max_elec_nb15 = max(k*2,min_atom)
      call alloc_pairlist(pairlist, PairListCGEleList, &
                          Max_elec_nb15, Max_cg_Elec)
    else
      call alloc_pairlist(pairlist, PairListCGEleList, &
                          1, 1)
    end if

    if (enefunc%cg_IDR_KH_calc) then
      num_atom = enefunc%num_cg_IDR_KH
      k = 0
      do i = 1, num_atom
        k = max(k,num_cg_idr_kh_calc(i))
      end do
#ifdef HAVE_MPI_GENESIS
      call mpi_allreduce(mpi_in_place, k, 1, mpi_integer, &
                         mpi_max, mpi_comm_country, ierror)
#endif
      Max_IDR_KH_nb15 = max(k*2, min_atom)
      call alloc_pairlist(pairlist, PairListCGIDRKHList, &
                          Max_IDR_KH_nb15, Max_cg_IDR_KH)
    else
      call alloc_pairlist(pairlist, PairListCGIDRKHList, &
                          1, 1)
    end if

    if (enefunc%cg_KH_calc) then
      num_atom = enefunc%num_cg_KH
      k = 0
      do i = 1, num_atom
        k = max(k,num_cg_kh_calc(i))
      end do
#ifdef HAVE_MPI_GENESIS
      call mpi_allreduce(mpi_in_place, k, 1, mpi_integer, &
                         mpi_max, mpi_comm_country, ierror)
#endif
      Max_KH_nb15 = max(k*2, min_atom)
      call alloc_pairlist(pairlist, PairListCGKHList, &
                          Max_KH_nb15, Max_cg_KH)
    else
      call alloc_pairlist(pairlist, PairListCGKHList, &
                          1, 1)
    end if

    if (enefunc%cg_IDR_HPS_calc) then
      num_atom = enefunc%num_cg_IDR_HPS
      k = 0
      do i = 1, num_atom
        k = max(k,num_cg_idr_hps_calc(i))
      end do
#ifdef HAVE_MPI_GENESIS
      call mpi_allreduce(mpi_in_place, k, 1, mpi_integer, &
                         mpi_max, mpi_comm_country, ierror)
#endif
      Max_IDR_HPS_nb15 = max(k*2, min_atom)
      call alloc_pairlist(pairlist, PairListCGIDRHPSList, &
                          Max_IDR_HPS_nb15, Max_cg_IDR_HPS)
    else
      call alloc_pairlist(pairlist, PairListCGIDRHPSList, &
                          1, 1)
    end if

    if (enefunc%cg_IDR_MPIPI_calc) then
      num_atom = enefunc%num_cg_IDR_MPIPI
      k = 0
      do i = 1, num_atom
        k = max(k,num_cg_idr_mpipi_calc(i))
      end do
#ifdef HAVE_MPI_GENESIS
      call mpi_allreduce(mpi_in_place, k, 1, mpi_integer, &
          mpi_max, mpi_comm_country, ierror)
#endif
      Max_IDR_MPIPI_nb15 = k * 2
      call alloc_pairlist(pairlist, PairListCGIDRMPIPIList, &
          Max_IDR_MPIPI_nb15, Max_cg_IDR_MPIPI)
    else
      call alloc_pairlist(pairlist, PairListCGIDRMPIPIList, &
          1, 1)
    end if

    return

  end subroutine update_pairlist_cg_alloc

 !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_cg
  !> @brief        Updates pairlists for AICG.
  !! @authors      JJ
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_cg(coord, enefunc, domain, pairlist)

    ! formal arguments
    real(wp),                 intent(in)    :: coord(:,:)
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(in)    :: domain
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: pairdist2_ele, pairdist2_126
    real(wp)                  :: pairdist2_WF
    real(wp)                  :: pairdist2_DNAbp, pairdist2_exv
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: rtmp(1:3)

    integer                   :: i, j, ij, k, ix, ixx, iy, iyy, inbc
    integer                   :: k_exv, k_dna, k_base, k_kh
    integer                   :: k_idr_kh, k_idr_hps, k_idr_mpipi
    integer                   :: k_elec
    integer                   :: chain_idx, chain_idy
    integer                   :: i_atom, j_atom, num_atom
    integer                   :: start_i, start_j
    integer                   :: i_base_type, j_base_type
    logical                   :: is_nonlocal, is_WC_bp
    integer                   :: num_nb15, num_nb15_dna, num_nb15_base
    integer                   :: num_nb15_kh, num_nb15_idr_kh
    integer                   :: num_nb15_idr_hps
    integer                   :: num_nb15_idr_mpipi
    integer                   :: i_dna, j_dna, ij_dna
    integer                   :: i_kh, j_kh, ij_kh
    integer                   :: i_idr_kh, j_idr_kh, ij_idr_kh
    integer                   :: i_idr_hps, j_idr_hps, ij_idr_hps
    integer                   :: i_idr_mpipi, j_idr_mpipi, ij_idr_mpipi
    integer                   :: i_idr, j_idr, ij_idr, ij_idr1
    integer                   :: i_126
    integer                   :: id
    integer                   :: ncell, nboundary
    logical                   :: nb15_calc

    integer,          pointer :: natom(:), ncharge(:)
    integer,          pointer :: start_atom(:)
    integer,          pointer :: near_cells_count(:), far_cells_count(:)
    integer,          pointer :: near_cells(:,:), far_cells(:,:)
    integer,          pointer :: id_l2g(:)
    integer(1),       pointer :: dna_check(:)
    integer,          pointer :: base_type(:)
    integer(1),       pointer :: cg_pro_use_KH(:), cg_IDR_KH(:), cg_IDR_HPS(:)
    integer(1),       pointer :: cg_IDR_MPIPI(:)
    integer,          pointer :: cg_exv_list(:,:)
    integer,          pointer :: cg_DNA_exv_list(:,:)
    integer,          pointer :: cg_DNA_base_list(:,:)
    integer,          pointer :: cg_ele_list(:,:)
    integer,          pointer :: cg_elec_list_inv(:)
    integer,          pointer :: cg_DNA_list_inv(:)
    integer,          pointer :: cg_base_list_inv(:)
    integer,          pointer :: cg_kh_list_inv(:)
    integer,          pointer :: cg_idr_kh_list_inv(:)
    integer,          pointer :: cg_idr_hps_list_inv(:)
    integer,          pointer :: cg_idr_mpipi_list_inv(:)
    integer,          pointer :: num_cg_exv_calc(:)
    integer,          pointer :: num_cg_DNA_exv_calc(:)
    integer,          pointer :: num_cg_DNA_base_calc(:)
    integer,          pointer :: num_cg_ele_calc(:)
    integer,          pointer :: num_cg_kh_calc(:)
    integer,          pointer :: num_cg_idr_kh_calc(:)
    integer,          pointer :: num_cg_idr_hps_calc(:)
    integer,          pointer :: num_cg_idr_mpipi_calc(:)
    integer(1),       pointer :: exclusion_mask(:,:)
    integer,          pointer :: chain_id(:)
    integer(1),       pointer :: cg_ele_mol_pair(:,:), cg_kh_mol_pair(:,:)
    logical,          pointer :: base_pair_is_WC(:,:)
    integer,          pointer :: cg_kh_list(:,:), cg_idr_kh_list(:,:)
    integer,          pointer :: cg_idr_hps_list(:,:)
    integer,          pointer :: cg_idr_mpipi_list(:,:)
    integer,          pointer :: atom_2_cell(:)

    ncell                =  domain%num_cell_local
    nboundary            =  domain%num_cell_boundary
    natom                => domain%num_atom
    start_atom           => domain%start_atom
    ncharge              => domain%num_charge
    near_cells_count     => domain%near_cells_count
    far_cells_count      => domain%far_cells_count
    near_cells           => domain%near_cells      
    far_cells            => domain%far_cells      
    dna_check            => domain%dna_check
    cg_pro_use_KH        => domain%cg_pro_use_KH
    cg_IDR_KH            => domain%cg_IDR_KH
    cg_IDR_HPS           => domain%cg_IDR_HPS
    cg_IDR_MPIPI         => domain%cg_IDR_MPIPI
    id_l2g               => domain%id_l2g
    chain_id             => domain%mol_chain_id
    base_type            => domain%NA_base_type
    atom_2_cell          => domain%atom_2_cell

    exclusion_mask       => enefunc%exclusion_mask
    cg_ele_mol_pair      => enefunc%cg_ele_mol_pair
    cg_kh_mol_pair       => enefunc%cg_kh_mol_pair
    base_pair_is_WC      => enefunc%base_pair_is_WC
    cg_elec_list_inv     => enefunc%cg_elec_list_inv
    cg_dna_list_inv      => enefunc%cg_dna_list_inv
    cg_base_list_inv     => enefunc%cg_base_list_inv
    cg_kh_list_inv       => enefunc%cg_kh_list_inv
    cg_idr_kh_list_inv   => enefunc%cg_idr_kh_list_inv
    cg_idr_hps_list_inv  => enefunc%cg_idr_hps_list_inv
    cg_idr_mpipi_list_inv => enefunc%cg_idr_mpipi_list_inv


    num_cg_exv_calc      => pairlist%num_cg_exv_calc
    cg_exv_list          => pairlist%cg_exv_list
    num_cg_DNA_exv_calc  => pairlist%num_cg_DNA_exv_calc
    cg_DNA_exv_list      => pairlist%cg_DNA_exv_list
    num_cg_ele_calc      => pairlist%num_cg_ele_calc
    cg_ele_list          => pairlist%cg_ele_list
    num_cg_DNA_base_calc => pairlist%num_cg_DNA_base_calc
    cg_DNA_base_list     => pairlist%cg_DNA_base_list
    num_cg_kh_calc       => pairlist%num_cg_kh_calc
    cg_kh_list           => pairlist%cg_kh_list
    num_cg_idr_kh_calc   => pairlist%num_cg_idr_kh_calc
    cg_idr_kh_list       => pairlist%cg_idr_kh_list
    num_cg_idr_hps_calc  => pairlist%num_cg_idr_hps_calc
    cg_idr_hps_list      => pairlist%cg_idr_hps_list
    num_cg_idr_mpipi_calc => pairlist%num_cg_idr_mpipi_calc
    cg_idr_mpipi_list     => pairlist%cg_idr_mpipi_list

    pairdist2_ele        =  pairlist%cg_pairlistdist_ele**2
    pairdist2_126        =  pairlist%cg_pairlistdist_126**2
    pairdist2_WF         =  pairlist%cg_pairlistdist_WF**2
    pairdist2_DNAbp      =  pairlist%cg_pairlistdist_DNAbp**2
    pairdist2_exv        =  pairlist%cg_pairlistdist_exv**2

    ! initialize the key of reallocation
    pairlist%realloc_exv = 0
    pairlist%realloc_DNA_exv = 0
    pairlist%realloc_DNA_base = 0
    pairlist%realloc_elec = 0
    pairlist%realloc_IDR_KH = 0
    pairlist%realloc_KH = 0
    pairlist%realloc_IDR_HPS = 0
    pairlist%realloc_IDR_MPIPI = 0

    num_atom = domain%num_atom_domain + domain%num_atom_boundary
    num_cg_exv_calc(1:num_atom) = 0
    if (enefunc%cg_DNA_exv_calc) then
      num_atom = enefunc%num_cg_dna
      num_cg_DNA_exv_calc(1:num_atom) = 0
    end if
    if (enefunc%cg_DNA_base_pair_calc) then
      num_atom = enefunc%num_cg_base
      num_cg_DNA_base_calc(1:num_atom) = 0
    end if
    if (enefunc%cg_ele_calc) then
      num_atom = enefunc%num_cg_elec
      num_cg_ele_calc(1:num_atom) = 0
    end if
    if (enefunc%cg_IDR_KH_calc) then
      num_atom = enefunc%num_cg_IDR_KH
      pairlist%num_cg_idr_kh_calc(1:num_atom) = 0
    end if
    if (enefunc%cg_KH_calc) then
      num_atom = enefunc%num_cg_KH
      pairlist%num_cg_kh_calc(1:num_atom) = 0
    end if
    if (enefunc%cg_IDR_HPS_calc) then
      num_atom = enefunc%num_cg_IDR_HPS
      pairlist%num_cg_idr_hps_calc(1:num_atom) = 0
    end if
    if (enefunc%cg_IDR_MPIPI_calc) then
      num_atom = enefunc%num_cg_IDR_MPIPI
      pairlist%num_cg_idr_mpipi_calc(1:num_atom) = 0
    end if

    !$omp parallel default(shared)                                        &
    !$omp private(id, i, ix, iy, k, ij, j, inbc, num_nb15, num_nb15_dna,  &
    !$omp         num_nb15_base, num_nb15_kh, num_nb15_idr_kh,            &
    !$omp         num_nb15_idr_hps, i_dna, j_dna, ij_dna, i_126,          &
    !$omp         num_nb15_idr_mpipi,                                     &
    !$omp         rtmp, dij, rij2, chain_idx, chain_idy, nb15_calc,       &
    !$omp         i_atom, j_atom, i_base_type, j_base_type, is_nonlocal,  &
    !$omp         is_WC_bp, start_i, start_j, ixx, iyy, i_kh, j_kh,       &
    !$omp         ij_kh, i_idr_kh, j_idr_kh, ij_idr_kh, i_idr_hps,        &
    !$omp         j_idr_hps, ij_idr_hps, i_idr, j_idr, ij_idr, ij_idr1,   &
    !$omp         j_idr_mpipi, ij_idr_mpipi, i_idr_mpipi, k_idr_mpipi,    &
    !$omp         k_exv, k_base, k_dna, k_kh, k_idr_kh, k_idr_hps, k_elec)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    ! make a pairlist in the same cell
    !
    do ixx = id+1, domain%num_atom_domain+domain%num_atom_boundary, nthread

      i = atom_2_cell(ixx)
      rtmp(1:3) = coord(ixx,1:3)

      num_nb15_dna = 0
      num_nb15_base = 0
      num_nb15 = 0
      num_nb15_kh = 0
      num_nb15_idr_kh = 0
      num_nb15_idr_hps = 0
      num_nb15_idr_mpipi = 0
      i_dna     = dna_check(ixx)
      i_kh      = cg_pro_use_KH(ixx)
      i_idr_kh  = cg_IDR_KH(ixx)
      i_idr_hps = cg_IDR_HPS(ixx)
      i_idr_mpipi = cg_IDR_MPIPI(ixx)
      i_idr     = i_idr_kh + i_idr_hps + i_idr_mpipi
      i_126     = i_idr + i_kh

      if (i_dna /= 0) then
        chain_idx = chain_id(ixx)
        i_base_type = base_type(ixx)
        i_atom = id_l2g(ixx)
      else
        chain_idx = chain_id(ixx)
        i_atom = id_l2g(ixx)
      end if

      k_exv = ixx
      if (enefunc%cg_DNA_exv_calc) k_dna = cg_dna_list_inv(ixx)
      if (enefunc%cg_DNA_base_pair_calc) k_base = cg_base_list_inv(ixx)
      if (enefunc%cg_KH_calc) k_kh = cg_kh_list_inv(ixx)
      if (enefunc%cg_IDR_KH_calc) k_idr_kh = cg_idr_kh_list_inv(ixx)
      if (enefunc%cg_IDR_HPS_calc) k_idr_hps = cg_idr_hps_list_inv(ixx)
      if (enefunc%cg_IDR_MPIPI_calc) k_idr_mpipi = cg_idr_mpipi_list_inv(ixx)

      if (i <= ncell) then

        do iyy = ixx+1, start_atom(i+1) 

          if (exclusion_mask(iyy,ixx) /= 1) then

            dij(1:3) = rtmp(1:3) - coord(iyy,1:3)
            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            ! case 1: rij < pairlistdist_exv : write all interactions
            !
            if (rij2 < pairdist2_exv) then

              call pairlist_exv_DNA_idr(k_dna, k_base, k_kh, k_idr_kh,        &
                                k_idr_hps, k_idr_mpipi,                       &
                                dna_check, cg_pro_use_KH,                     &
                                cg_IDR_KH, cg_IDR_HPS, cg_IDR_MPIPI,          &
                                chain_id, base_type,                          &
                                base_pair_is_WC, id_l2g, cg_kh_mol_pair, ixx, &
                                iyy, i_atom, chain_idx, i_idr, i_dna,         &
                                i_base_type, i_kh, i_idr_kh, i_idr_hps,       &
                                i_idr_mpipi,                                  &
                                num_nb15_base, num_nb15_dna, num_nb15,        &
                                num_nb15_kh, num_nb15_idr_kh,                 &
                                num_nb15_idr_hps, num_nb15_idr_mpipi,         &
                                cg_DNA_base_list,                             &
                                cg_DNA_exv_list, cg_exv_list, cg_kh_list,     &
                                cg_idr_kh_list, cg_idr_hps_list,              &
                                cg_idr_mpipi_list)

            ! case 1: pairlistdist_exv < rij < pairlistdist_DNAbp 
            !         interaction of Base pair and IDR
            !
            else if ((i_dna>0.or.i_126>0) .and. rij2 < pairdist2_DNAbp) then

              call pairlist_DNA_idr(k_base, k_kh, k_idr_kh, k_idr_hps,       &
                                k_idr_mpipi,                                 &
                                dna_check, cg_pro_use_KH, cg_IDR_KH,         &
                                cg_IDR_HPS, cg_IDR_MPIPI,                    &
                                chain_id, base_type,                         &
                                base_pair_is_WC, id_l2g, cg_kh_mol_pair,     &
                                iyy, i_atom, chain_idx, i_idr, i_dna,        &
                                i_base_type, i_kh, i_idr_kh, i_idr_hps,      &
                                i_idr_mpipi,                                 &
                                num_nb15_base, num_nb15_kh, num_nb15_idr_kh, &
                                num_nb15_idr_hps, num_nb15_idr_mpipi,        &
                                cg_DNA_base_list,                            &
                                cg_kh_list, cg_idr_kh_list, cg_idr_hps_list, &
                                cg_idr_mpipi_list)

            ! case 3: pairlist2_DNAbp < rij < pairlist2_126
            !
            else if (i_126 > 0 .and. rij2 < pairdist2_126) then

              call pairlist_idr(k_kh, k_idr_kh, k_idr_hps, k_idr_mpipi,      &
                                cg_pro_use_KH,                               &
                                cg_IDR_KH, cg_IDR_HPS, cg_IDR_MPIPI,         &
                                chain_id, id_l2g,                            &
                                cg_kh_mol_pair, iyy, i_atom, chain_idx,      &
                                i_idr, i_kh, i_idr_kh, i_idr_hps,            &
                                i_idr_mpipi, num_nb15_kh,                    &
                                num_nb15_idr_kh, num_nb15_idr_hps,           &
                                num_nb15_idr_mpipi,                          &
                                cg_kh_list, cg_idr_kh_list, cg_idr_hps_list, &
                                cg_idr_mpipi_list)

            end if

          end if

        end do

      end if

      do inbc = 1, near_cells_count(i)

        j = near_cells(inbc,i)
        start_j = start_atom(j)

        do iy = 1, natom(j)
     
          iyy = iy + start_j

          if (exclusion_mask(iyy,ixx) /= 1) then

            dij(1:3) = rtmp(1:3) - coord(iyy,1:3)
            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            ! case 1: rij < pairlistdist_exv : write all interactions
            !
            if (rij2 < pairdist2_exv) then

              call pairlist_exv_DNA_idr(k_dna, k_base, k_kh, k_idr_kh,        &
                                k_idr_hps, k_idr_mpipi,                       &
                                dna_check, cg_pro_use_KH,                     &
                                cg_IDR_KH, cg_IDR_HPS, cg_IDR_MPIPI,          &
                                chain_id, base_type,                          &
                                base_pair_is_WC, id_l2g, cg_kh_mol_pair, ixx, &
                                iyy, i_atom, chain_idx, i_idr, i_dna,         &
                                i_base_type, i_kh, i_idr_kh, i_idr_hps,       &
                                i_idr_mpipi,                                  &
                                num_nb15_base, num_nb15_dna, num_nb15,        &
                                num_nb15_kh, num_nb15_idr_kh,                 &
                                num_nb15_idr_hps, num_nb15_idr_mpipi,         &
                                cg_DNA_base_list,                             &
                                cg_DNA_exv_list, cg_exv_list, cg_kh_list,     &
                                cg_idr_kh_list, cg_idr_hps_list,              &
                                cg_idr_mpipi_list)


            ! case 1: pairlistdist_exv < rij < pairlistdist_DNAbp 
            !         interaction of Base pair and IDR
            !
            else if ((i_dna>0.or.i_126>0) .and. rij2 < pairdist2_DNAbp) then

              call pairlist_DNA_idr(k_base, k_kh, k_idr_kh, k_idr_hps,       &
                                k_idr_mpipi,                                 &
                                dna_check, cg_pro_use_KH, cg_IDR_KH,         &
                                cg_IDR_HPS, cg_IDR_MPIPI,                    &
                                chain_id, base_type,                         &
                                base_pair_is_WC, id_l2g, cg_kh_mol_pair,     &
                                iyy, i_atom, chain_idx, i_idr, i_dna,        &
                                i_base_type, i_kh, i_idr_kh, i_idr_hps,      &
                                i_idr_mpipi,                                 &
                                num_nb15_base, num_nb15_kh, num_nb15_idr_kh, &
                                num_nb15_idr_hps, num_nb15_idr_mpipi,        &
                                cg_DNA_base_list,                            &
                                cg_kh_list, cg_idr_kh_list, cg_idr_hps_list, &
                                cg_idr_mpipi_list)


            ! case 3: pairlist2_DNAbp < rij < pairlist2_126
            !
            else if (i_126 > 0 .and. rij2 < pairdist2_126) then

              call pairlist_idr(k_kh, k_idr_kh, k_idr_hps, k_idr_mpipi,      &
                                cg_pro_use_KH,                               &
                                cg_IDR_KH, cg_IDR_HPS, cg_IDR_MPIPI,         &
                                chain_id, id_l2g,                            &
                                cg_kh_mol_pair, iyy, i_atom, chain_idx,      &
                                i_idr, i_kh, i_idr_kh, i_idr_hps,            &
                                i_idr_mpipi, num_nb15_kh,                    &
                                num_nb15_idr_kh, num_nb15_idr_hps,           &
                                num_nb15_idr_mpipi,                          &
                                cg_kh_list, cg_idr_kh_list, cg_idr_hps_list, &
                                cg_idr_mpipi_list)

            end if

          end if

        end do

      end do

      if (i_126 > 0) then

        do inbc = 1, far_cells_count(i)

          j = far_cells(inbc,i)
          start_j = start_atom(j)

          do iy = 1, natom(j)

            iyy = start_j + iy
            dij(1:3) = rtmp(1:3) - coord(iyy,1:3)
            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            if (rij2 < pairdist2_126) then

              call pairlist_idr(k_kh, k_idr_kh, k_idr_hps, k_idr_mpipi,      &
                                cg_pro_use_KH,                               &
                                cg_IDR_KH, cg_IDR_HPS, cg_IDR_MPIPI,         &
                                chain_id, id_l2g,                            &
                                cg_kh_mol_pair, iyy, i_atom, chain_idx,      &
                                i_idr, i_kh, i_idr_kh, i_idr_hps,            &
                                i_idr_mpipi, num_nb15_kh,                    &
                                num_nb15_idr_kh, num_nb15_idr_hps,           &
                                num_nb15_idr_mpipi,                          &
                                cg_kh_list, cg_idr_kh_list, cg_idr_hps_list, &
                                cg_idr_mpipi_list)

            end if

          end do

        end do

      end if

      if (num_nb15 > 0) num_cg_exv_calc(ixx) = num_nb15
      if (num_nb15 > Max_exv_nb15*0.7) pairlist%realloc_exv = 1
      if (num_nb15_dna > 0) num_cg_DNA_exv_calc(k_dna) = num_nb15_dna
      if (num_nb15_dna > Max_dna_exv_nb15*0.7) pairlist%realloc_DNA_exv = 1
      if (num_nb15_base > 0) num_cg_DNA_base_calc(k_base) = num_nb15_base
      if (num_nb15_base > Max_dna_base_nb15*0.7) pairlist%realloc_DNA_base = 1
      if (num_nb15_kh > 0) num_cg_kh_calc(k_kh) = num_nb15_kh
      if (num_nb15_kh > Max_KH_nb15*0.7) pairlist%realloc_KH = 1
      if (num_nb15_idr_kh > 0) num_cg_idr_kh_calc(k_idr_kh) = num_nb15_idr_kh
      if (num_nb15_idr_kh > Max_IDR_KH_nb15*0.7) pairlist%realloc_IDR_KH = 1
      if (num_nb15_idr_hps > 0) &
          num_cg_idr_hps_calc(k_idr_hps) = num_nb15_idr_hps
      if (num_nb15_idr_hps > Max_IDR_HPS_nb15*0.7) pairlist%realloc_IDR_HPS = 1
      if (num_nb15_idr_mpipi > 0) &
          num_cg_idr_mpipi_calc(k_idr_mpipi) = num_nb15_idr_mpipi
      if (num_nb15_idr_mpipi > Max_IDR_MPIPI_nb15*0.7) pairlist%realloc_IDR_MPIPI = 1

    end do

    ! electrostatic
    !
    if (enefunc%cg_ele_calc) then

      do i = id+1, ncell+nboundary, nthread

        start_i = start_atom(i)

        do ix = 1, ncharge(i) 

          num_nb15 = 0
          ixx = ix + start_i
          k_elec = cg_elec_list_inv(ixx)
          chain_idx = chain_id(ixx)
          rtmp(1:3) = coord(ixx,1:3)

          if (i <= ncell) then

            do iy = ix + 1, ncharge(i)
  
              iyy = iy + start_i
  
              if (exclusion_mask(iyy,ixx) /= 1) then
  
                chain_idy = chain_id(iyy)
                if (cg_ele_mol_pair(chain_idx,chain_idy) /= 0) then
  
                  dij(1:3) = rtmp(1:3) - coord(iyy,1:3)
                  rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
   
                  ! store interaction table
                  !
                  if (rij2 < pairdist2_ele) then
  
                    num_nb15 = num_nb15 + 1
                    pairlist%cg_ele_list(num_nb15,k_elec) = iyy
  
                  end if
                end if
              end if
            end do
          end if

          do inbc = 1, near_cells_count(i)

            j = near_cells(inbc,i)
            start_j = start_atom(j)

            do iy = 1, ncharge(j)

              iyy = iy + start_j

              if (exclusion_mask(iyy,ixx) /= 1) then

                chain_idy = chain_id(iyy)
                if (cg_ele_mol_pair(chain_idx,chain_idy) /= 0) then

                  dij(1:3) = rtmp(1:3) - coord(iyy,1:3)
                  rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

                  ! store interaction table
                  !
                  if (rij2 < pairdist2_ele) then

                    num_nb15 = num_nb15 + 1
                    pairlist%cg_ele_list(num_nb15,k_elec) = iyy

                  end if
                end if
              end if
            end do
          end do

          do inbc = 1, far_cells_count(i)

            j = far_cells(inbc,i)
            start_j = start_atom(j)

            do iy = 1, ncharge(j)

              iyy = iy + start_j
              chain_idy = chain_id(iyy)
              if (cg_ele_mol_pair(chain_idx,chain_idy) /= 0) then

                dij(1:3) = rtmp(1:3) - coord(iyy,1:3)
                rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

                ! store interaction table
                !
                if (rij2 < pairdist2_ele) then

                  num_nb15 = num_nb15 + 1
                  pairlist%cg_ele_list(num_nb15,k_elec) = iyy

                end if
              end if
            end do
          end do

          pairlist%num_cg_ele_calc(k_elec) = num_nb15
          if (num_nb15 > Max_elec_nb15*0.7) pairlist%realloc_elec = 1

        end do

      end do

    end if

    !$omp end parallel

    return

  end subroutine update_pairlist_cg

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_cg_pwmcos
  !> @brief        update pairlist for PWMCOS
  !! @authors      JJ
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_cg_pwmcos_alloc(coord, enefunc, domain, pairlist)

    ! formal arguments
    real(wp),                 intent(in)    :: coord(:,:)
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(in)    :: domain
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: pairdist2_PWMCos
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: rtmp(1:3)

    integer                   :: i, j, ij, k, i1, j1, ix, iy, ig
    integer                   :: icel, i_atom
    integer                   :: i_chain_id, j_chain_id
    integer                   :: start_i, start_j, start_icel
    integer                   :: iter, num_nb15
    integer                   :: id
    integer                   :: ncell, nboundary

    integer,          pointer :: nbase(:)
    integer,          pointer :: start_atom(:)
    integer,          pointer :: base_list(:), pwmcos_id(:)
    integer,          pointer :: chain_id(:)
    integer,          pointer :: id_l2g(:), id_g2l(:)
    integer,          pointer :: atom_2_cell(:)
    integer,          pointer :: num_cg_pwmcos_calc(:)
    integer(1),       pointer :: pwmcos_mol_pair(:,:)

    ncell                =  domain%num_cell_local
    nboundary            =  domain%num_cell_boundary
    nbase                => domain%num_base
    base_list            => domain%base_list
    start_atom           => domain%start_atom
    id_l2g               => domain%id_l2g
    id_g2l               => domain%id_g2l
    chain_id             => domain%mol_chain_id
    atom_2_cell          => domain%atom_2_cell
    pwmcos_id            => enefunc%pwmcos_protein_id
    pwmcos_mol_pair      => enefunc%pwmcos_mol_pair

    num_cg_pwmcos_calc   => pairlist%num_cg_pwmcos_calc

    pairdist2_PWMcos     =  pairlist%cg_pairlistdist_PWMcos ** 2 

    num_cg_pwmcos_calc(1:enefunc%num_pwmcos_domain) = 0

    !$omp parallel                                                       &
    !$omp private(id, i, ix, iy, k, ij, j, i1, j1, icel, i_atom, ig,     &
    !$omp         num_nb15, rtmp, dij, rij2, i_chain_id, j_chain_id,     &
    !$omp         iter, start_i, start_j, start_icel)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i1 = id+1, enefunc%num_pwmcos_domain, nthread

      num_nb15 = 0
      ix = id_g2l(pwmcos_id(i1))
      rtmp(1:3) = coord(ix,1:3)
      i_chain_id = chain_id(ix)
      i = atom_2_cell(ix)

      do ij = 1, 27

        j = domain%near_neighbor_cells(ij,i)
        start_j = start_atom(j)

        do j1 = 1, nbase(j)

          iy = base_list(j1+start_j)
          j_chain_id = chain_id(iy)

          if (pwmcos_mol_pair(i_chain_id,j_chain_id) /= 0) then

            dij(1:3) = rtmp(1:3) - coord(iy,1:3)
            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            if (rij2 < pairdist2_PWMcos) then

              ig = id_l2g(iy)
              if ((ig-3 >= 1) .and. (ig+3 <= domain%num_atom_all)) then
                i_atom = id_g2l(ig-3)
                if (chain_id(i_atom) == j_chain_id) then
                  i_atom = id_g2l(ig+3)
                  if (chain_id(i_atom) == j_chain_id) then
                    num_nb15 = num_nb15 + 1
                  end if
                end if
              end if
            end if
          end if

        end do

      end do

      num_cg_pwmcos_calc(i1) = num_nb15

    end do

    !$omp end parallel

    k = 0
    do i = 1, enefunc%num_pwmcos_domain
      k = max(k,num_cg_pwmcos_calc(i))
    end do
#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, k, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
#endif
    Max_pwmcos_nb15 = 2*k
    k = Max_pwmcos_nb15
    call alloc_pairlist(pairlist, PairListPwmCosList, &
                        k, MaxPwmCos)

    return

  end subroutine update_pairlist_cg_pwmcos_alloc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_cg_pwmcos
  !> @brief        update pairlist for PWMCOS
  !! @authors      JJ
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_cg_pwmcos(coord, enefunc, domain, pairlist)

    ! formal arguments
    real(wp),                 intent(in)    :: coord(:,:)
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(in)    :: domain
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: pairdist2_PWMCos
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: rtmp(1:3)

    integer                   :: i, j, ij, k, i1, i2, j1, ix, iy, ig
    integer                   :: icel, i_atom
    integer                   :: i_chain_id, j_chain_id
    integer                   :: start_i, start_j, start_icel
    integer                   :: iter, num_nb15
    integer                   :: id
    integer                   :: ncell, nboundary

    integer,          pointer :: nbase(:)
    integer,          pointer :: start_atom(:)
    integer,          pointer :: base_list(:), pwmcos_id(:)
    integer,          pointer :: chain_id(:)
    integer,          pointer :: id_l2g(:), id_g2l(:)
    integer,          pointer :: atom_2_cell(:)
    integer,          pointer :: cg_pwmcos_list(:,:)
    integer,          pointer :: num_cg_pwmcos_calc(:)
    integer(1),       pointer :: pwmcos_mol_pair(:,:)

    ncell                =  domain%num_cell_local
    nboundary            =  domain%num_cell_boundary
    nbase                => domain%num_base
    base_list            => domain%base_list
    start_atom           => domain%start_atom
    id_l2g               => domain%id_l2g
    id_g2l               => domain%id_g2l
    atom_2_cell          => domain%atom_2_cell
    chain_id             => domain%mol_chain_id
    pwmcos_id            => enefunc%pwmcos_protein_id
    pwmcos_mol_pair      => enefunc%pwmcos_mol_pair

    num_cg_pwmcos_calc   => pairlist%num_cg_pwmcos_calc
    cg_pwmcos_list       => pairlist%cg_pwmcos_list

    pairdist2_PWMcos     =  pairlist%cg_pairlistdist_PWMcos ** 2 
    pairlist%realloc_PWMcos = 0

    k = enefunc%num_pwmcos_domain
    num_cg_pwmcos_calc(1:k) = 0

    !$omp parallel                                                       &
    !$omp private(id, i, ix, iy, k, ij, j, i1, i2, j1, icel, i_atom, ig, &
    !$omp         num_nb15, rtmp, dij, rij2, i_chain_id, j_chain_id,     &
    !$omp         iter, start_i, start_j, start_icel)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i1 = id+1, enefunc%num_pwmcos_domain, nthread

      num_nb15 = 0
      ix = id_g2l(pwmcos_id(i1))
      rtmp(1:3) = coord(ix,1:3)
      i_chain_id = chain_id(ix)
      i = atom_2_cell(ix)

      do ij = 1, 27

        j = domain%near_neighbor_cells(ij,i)
        start_j = start_atom(j)

        do j1 = 1, nbase(j)

          iy = base_list(j1+start_j)
          j_chain_id = chain_id(iy)

          if (pwmcos_mol_pair(i_chain_id,j_chain_id) /= 0) then

            dij(1:3) = rtmp(1:3) - coord(iy,1:3)
            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

! store interaction table
!
            if (rij2 < pairdist2_PWMcos) then

              ig = id_l2g(iy)
              if ((ig-3 >= 1) .and. (ig+3 <= domain%num_atom_all)) then
                i_atom = id_g2l(ig-3)
                if (chain_id(i_atom) == j_chain_id) then
                  i_atom = id_g2l(ig+3)
                  if (chain_id(i_atom) == j_chain_id) then
                    num_nb15 = num_nb15 + 1
                    cg_pwmcos_list(num_nb15,i1) = iy
                  end if
                end if
              end if
            end if
          end if

        end do

      end do

      num_cg_pwmcos_calc(i1) = num_nb15
      if (num_nb15 > Max_pwmcos_nb15*3/4) pairlist%realloc_PWMcos = 1

    end do

    !$omp end parallel

    return

  end subroutine update_pairlist_cg_pwmcos

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_cg_pwmcosns
  !> @brief        update pairlist for PWMCOSns
  !! @authors      JJ
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_cg_pwmcosns_alloc(coord, enefunc, domain, pairlist)

    ! formal arguments
    real(wp),                 intent(in)    :: coord(:,:)
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(in)    :: domain
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: pairdist2_PWMcos
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: rtmp(1:3)

    integer                   :: i, j, ij, k, i1, j1, ix, iy, ig
    integer                   :: icel, i_atom
    integer                   :: i_chain_id, j_chain_id
    integer                   :: start_i, start_j, start_icel
    integer                   :: iter, num_nb15
    integer                   :: id
    integer                   :: ncell, nboundary

    integer,          pointer :: nphos(:)
    integer,          pointer :: start_atom(:)
    integer,          pointer :: phos_list(:), pwmcosns_id(:)
    integer,          pointer :: chain_id(:)
    integer,          pointer :: id_l2g(:), id_g2l(:)
    integer,          pointer :: atom_2_cell(:)
    integer,          pointer :: cg_pwmcosns_list(:,:)
    integer,          pointer :: num_cg_pwmcosns_calc(:)
    integer(1),       pointer :: pwmcosns_mol_pair(:,:)

    ncell                =  domain%num_cell_local
    nboundary            =  domain%num_cell_boundary
    nphos                => domain%num_phos
    phos_list            => domain%phos_list
    start_atom           => domain%start_atom
    id_l2g               => domain%id_l2g
    id_g2l               => domain%id_g2l
    chain_id             => domain%mol_chain_id
    atom_2_cell          => domain%atom_2_cell
    pwmcosns_id          => enefunc%pwmcosns_protein_id
    pwmcosns_mol_pair    => enefunc%pwmcosns_mol_pair

    num_cg_pwmcosns_calc => pairlist%num_cg_pwmcosns_calc
    cg_pwmcosns_list     => pairlist%cg_pwmcosns_list

    pairdist2_PWMcos     =  pairlist%cg_pairlistdist_PWMcos ** 2

    num_cg_pwmcosns_calc(1:enefunc%num_pwmcosns_domain) = 0

    !$omp parallel                                                       &
    !$omp private(id, i, ix, iy, k, ij, j, i1, j1, icel, i_atom, ig,     &
    !$omp         num_nb15, rtmp, dij, rij2, i_chain_id, j_chain_id,     &
    !$omp         iter, start_i, start_j, start_icel)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i1 = id+1, enefunc%num_pwmcosns_domain, nthread

      ix = id_g2l(pwmcosns_id(i1))
      rtmp(1:3) = coord(ix,1:3)
      num_nb15 = 0
      i_chain_id = chain_id(ix)
      i = atom_2_cell(ix)

      do ij = 1, 27

        j = domain%near_neighbor_cells(ij,i)
        start_j = start_atom(j)

        do j1 = 1, nphos(j)

          iy = phos_list(j1+start_j)
          j_chain_id = chain_id(iy)

          if (pwmcosns_mol_pair(i_chain_id,j_chain_id) /= 0) then

            dij(1:3) = rtmp(1:3) - coord(iy,1:3)
            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

! store interaction table
!
            if (rij2 < pairdist2_PWMcos) then
              num_nb15 = num_nb15 + 1
            end if
          end if

        end do

      end do

      num_cg_pwmcosns_calc(i1) = num_nb15

    end do

    !$omp end parallel

    k = 0
    do i = 1, enefunc%num_pwmcosns_domain
      k = max(k,num_cg_pwmcosns_calc(i))
    end do
#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, k, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
#endif
    Max_pwmcosns_nb15 = 2*k
    k = Max_pwmcosns_nb15
    call alloc_pairlist(pairlist, PairListPwmCosnsList, &
                        k, MaxPwmCosns)

    return

  end subroutine update_pairlist_cg_pwmcosns_alloc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_cg_pwmcosns
  !> @brief        update pairlist for PWMCOSns
  !! @authors      JJ
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_cg_pwmcosns(coord, enefunc, domain, pairlist)

    ! formal arguments
    real(wp),                 intent(in)    :: coord(:,:)
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(in)    :: domain
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: pairdist2_PWMcos
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: rtmp(1:3)

    integer                   :: i, j, ij, k, i1, i2, j1, ix, iy, ig
    integer                   :: icel, i_atom
    integer                   :: i_chain_id, j_chain_id
    integer                   :: start_i, start_j, start_icel
    integer                   :: iter, num_nb15
    integer                   :: id
    integer                   :: ncell, nboundary

    integer,          pointer :: nphos(:)
    integer,          pointer :: start_atom(:)
    integer,          pointer :: phos_list(:), pwmcosns_id(:)
    integer,          pointer :: chain_id(:)
    integer,          pointer :: id_l2g(:), id_g2l(:)
    integer,          pointer :: atom_2_cell(:)
    integer,          pointer :: cg_pwmcosns_list(:,:)
    integer,          pointer :: num_cg_pwmcosns_calc(:)
    integer(1),       pointer :: pwmcosns_mol_pair(:,:)

    ncell                =  domain%num_cell_local
    nboundary            =  domain%num_cell_boundary
    nphos                => domain%num_phos
    phos_list            => domain%phos_list
    start_atom           => domain%start_atom
    id_l2g               => domain%id_l2g
    id_g2l               => domain%id_g2l
    atom_2_cell          => domain%atom_2_cell
    chain_id             => domain%mol_chain_id
    pwmcosns_id          => enefunc%pwmcosns_protein_id
    pwmcosns_mol_pair    => enefunc%pwmcosns_mol_pair

    num_cg_pwmcosns_calc => pairlist%num_cg_pwmcosns_calc
    cg_pwmcosns_list     => pairlist%cg_pwmcosns_list

    pairdist2_PWMcos     =  pairlist%cg_pairlistdist_PWMcos ** 2
    pairlist%realloc_PWMcosns = 0


    num_cg_pwmcosns_calc(1:enefunc%num_pwmcosns_domain) = 0

    !$omp parallel                                                       &
    !$omp private(id, i, ix, iy, k, ij, j, i1, i2, j1, icel, i_atom, ig, &
    !$omp         num_nb15, rtmp, dij, rij2, i_chain_id, j_chain_id,     &
    !$omp         iter, start_i, start_j, start_icel)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i1 = id+1, enefunc%num_pwmcosns_domain, nthread

      ix = id_g2l(pwmcosns_id(i1))
      rtmp(1:3) = coord(ix,1:3)
      num_nb15 = 0
      i_chain_id = chain_id(ix)
      i = atom_2_cell(ix)

      do ij = 1, 27

        j = domain%near_neighbor_cells(ij,i)
        start_j = start_atom(j)

        do j1 = 1, nphos(j)

          iy = phos_list(j1+start_j)
          j_chain_id = chain_id(iy)

          if (pwmcosns_mol_pair(i_chain_id,j_chain_id) /= 0) then

            dij(1:3) = rtmp(1:3) - coord(iy,1:3)
            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

! store interaction table
!
            if (rij2 < pairdist2_PWMcos) then
              num_nb15 = num_nb15 + 1
              cg_pwmcosns_list(num_nb15,i1) = iy
            end if
          end if

        end do

      end do

      num_cg_pwmcosns_calc(i1) = num_nb15
      if (num_nb15 > Max_pwmcosns_nb15*3/4) pairlist%realloc_pwmcosns = 1

    end do

    !$omp end parallel

    return

  end subroutine update_pairlist_cg_pwmcosns

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_martini_alloc
  !> @brief        update pairlist for Martini (allocation)
  !! @authors      JJ
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_martini_alloc(coord, enefunc, domain, pairlist)

    ! formal arguments
    real(wp),                 intent(in)    :: coord(:,:)
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_domain),   target, intent(in)    :: domain
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: pairdist2_ele, pairdist2_vdw
    real(wp)                  :: dij(1:3), rij2, qtmp, term_elec, crf
    real(wp)                  :: rtmp(1:3), cutoff, el_fact
    real(dp)                  :: rf_self(nthread)

    integer                   :: i, j, ij, k, ix, ixx, iy, iyy, inbc
    integer                   :: i_atom, j_atom, num_atom
    integer                   :: start_i, start_j
    integer                   :: i_base_type, j_base_type
    integer                   :: num_nb15, num_nb15_ele
    integer                   :: id
    integer                   :: ncell, nboundary
    logical                   :: nb15_calc

    integer,          pointer :: natom(:), ncharge(:)
    integer,          pointer :: start_atom(:)
    integer,          pointer :: near_cells_count(:), far_cells_count(:)
    integer,          pointer :: near_cells(:,:), far_cells(:,:)
    integer,          pointer :: id_l2g(:)
    integer,          pointer :: cg_elec_list_inv(:)
    integer,          pointer :: num_cg_vdw_calc(:)
    integer,          pointer :: num_cg_ele_calc(:)
    integer(1),       pointer :: exclusion_mask(:,:)
    real(wp),         pointer :: charge(:)

    ncell                =  domain%num_cell_local
    nboundary            =  domain%num_cell_boundary
    natom                => domain%num_atom
    start_atom           => domain%start_atom
    ncharge              => domain%num_charge
    near_cells_count     => domain%near_cells_count
    far_cells_count      => domain%far_cells_count
    near_cells           => domain%near_cells
    far_cells            => domain%far_cells
    id_l2g               => domain%id_l2g
    charge               => domain%charge

    exclusion_mask       => enefunc%exclusion_mask
    cg_elec_list_inv     => enefunc%cg_elec_list_inv

    num_cg_vdw_calc      => pairlist%num_cg_exv_calc
    num_cg_ele_calc      => pairlist%num_cg_ele_calc

    pairdist2_ele        =  pairlist%cg_pairlistdist_ele**2
    pairdist2_vdw        =  pairlist%cg_pairlistdist_vdw**2

    num_atom = domain%num_atom_domain + domain%num_atom_boundary
    num_cg_vdw_calc(1:num_atom) = 0
    num_atom = enefunc%num_cg_elec
    num_cg_ele_calc(1:num_atom) = 0

    ! self correction term in reaction-field
    !
    cutoff = enefunc%cg_cutoffdist_ele
    if (enefunc%epsilon_rf < EPS) then
      crf = 1.0_wp/(2.0_wp*cutoff*cutoff*cutoff)
    else
      crf = (enefunc%epsilon_rf-enefunc%dielec_const) / &
            (2.0_wp*enefunc%epsilon_rf+enefunc%dielec_const)
      crf = crf / (cutoff*cutoff*cutoff)
    end if
    crf = 1.0_wp/cutoff + crf*cutoff*cutoff
    el_fact = ELECOEF / enefunc%dielec_const

    !$omp parallel default(shared)                                        &
    !$omp private(id, i, start_i, ix, ixx, rtmp, num_nb15, iy, iyy, dij,  &
    !$omp         rij2, inbc, j, start_j, qtmp, term_elec)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    ! make a pairlist in the same cell
    !
    do i = id+1, ncell+nboundary, nthread

      start_i = start_atom(i)

      do ix = 1, natom(i)  

        ixx = ix + start_i
        rtmp(1) = coord(ixx,1)
        rtmp(2) = coord(ixx,2)
        rtmp(3) = coord(ixx,3)
        num_nb15 = 0
      
        if (i <= ncell) then
          do iy = ix + 1, natom(i)
            iyy = iy + start_i
            if (exclusion_mask(iyy,ixx) /= 1) then
              dij(1) = rtmp(1) - coord(iyy,1)
              dij(2) = rtmp(2) - coord(iyy,2)
              dij(3) = rtmp(3) - coord(iyy,3)
              rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
              if (rij2 < pairdist2_vdw) num_nb15 = num_nb15 + 1
            end if
          end do
        end if

        do inbc = 1, near_cells_count(i)
          j = near_cells(inbc,i)
          start_j = start_atom(j)
          do iy = 1, natom(j)
            iyy = iy + start_j
            if (exclusion_mask(iyy,ixx) /= 1) then
              dij(1) = rtmp(1) - coord(iyy,1)
              dij(2) = rtmp(2) - coord(iyy,2)
              dij(3) = rtmp(3) - coord(iyy,3)
              rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
              if (rij2 < pairdist2_vdw) num_nb15 = num_nb15 + 1
            end if
          end do
        end do

        if (num_nb15 > 0) then
          num_cg_vdw_calc(ixx) = num_nb15
        end if

      end do

    end do

    ! electrostatic
    !
    rf_self(id+1) = 0.0_dp

    do i = id+1, ncell+nboundary, nthread

      start_i = start_atom(i)

      do ix = 1, ncharge(i) 

        num_nb15 = 0
        ixx = ix + start_i
        rtmp(1) = coord(ixx,1)
        rtmp(2) = coord(ixx,2)
        rtmp(3) = coord(ixx,3)

        if (i <= ncell) then

          if (enefunc%electrostatic == ElectrostaticRF) then
            qtmp = charge(ixx)
            term_elec = - qtmp*qtmp*crf
            term_elec = term_elec*el_fact
            rf_self(id+1) = rf_self(id+1) + 0.5_wp*term_elec
          end if

          do iy = ix + 1, ncharge(i)
            iyy = iy + start_i
            if (exclusion_mask(iyy,ixx) /= 1) then
              dij(1) = rtmp(1) - coord(iyy,1)
              dij(2) = rtmp(2) - coord(iyy,2)
              dij(3) = rtmp(3) - coord(iyy,3)
              rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
              if (rij2 < pairdist2_ele) num_nb15 = num_nb15 + 1
            end if
          end do
        end if

        do inbc = 1, near_cells_count(i)
          j = near_cells(inbc,i)
          start_j = start_atom(j)
          do iy = 1, ncharge(j)
            iyy = iy + start_j
            if (exclusion_mask(iyy,ixx) /= 1) then
              dij(1) = rtmp(1) - coord(iyy,1)
              dij(2) = rtmp(2) - coord(iyy,2)
              dij(3) = rtmp(3) - coord(iyy,3)
              rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
              if (rij2 < pairdist2_ele) num_nb15 = num_nb15 + 1
            end if
          end do
        end do

        k = cg_elec_list_inv(ixx)
        num_cg_ele_calc(k) = num_nb15

      end do
    end do

    !$omp end parallel

    do id = 1, nthread
      enefunc%rf_self = enefunc%rf_self + rf_self(id)
    end do

    num_atom = domain%num_atom_domain + domain%num_atom_boundary
    k = 0
    do i = 1, num_atom
      k = max(k,num_cg_vdw_calc(i))
    end do
#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, k, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
#endif
    Max_exv_nb15 = k * 2
    call alloc_pairlist(pairlist, PairListCGExvList, &
                        Max_exv_nb15, MaxAtom_domain)

    num_atom = enefunc%num_cg_elec
    k = 0
    do i = 1, num_atom
      k = max(k,num_cg_ele_calc(i))
    end do
#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, k, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
#endif
    if (Max_elec_nb15 > 0) Max_elec_nb15 = max(k*2,1000)
    call alloc_pairlist(pairlist, PairListCGEleList, &
                        Max_elec_nb15, Max_cg_Elec)

    return

  end subroutine update_pairlist_martini_alloc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_martini
  !> @brief        update pairlist for Martini
  !! @authors      JJ
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_martini(coord, enefunc, domain, pairlist)

    ! formal arguments
    real(wp),                 intent(in)    :: coord(:,:)
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(in)    :: domain
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: pairdist2_ele, pairdist2_vdw
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: rtmp(1:3)

    integer                   :: i, j, ij, k, ix, ixx, iy, iyy, inbc
    integer                   :: i_atom, j_atom, num_atom
    integer                   :: start_i, start_j
    integer                   :: i_base_type, j_base_type
    integer                   :: num_nb15, num_nb15_ele
    integer                   :: id
    integer                   :: ncell, nboundary
    logical                   :: nb15_calc

    integer,          pointer :: natom(:), ncharge(:)
    integer,          pointer :: start_atom(:), atom_2_cell(:)
    integer,          pointer :: near_cells_count(:), far_cells_count(:)
    integer,          pointer :: near_cells(:,:), far_cells(:,:)
    integer,          pointer :: id_l2g(:)
    integer,          pointer :: cg_elec_list_inv(:)
    integer,          pointer :: num_cg_vdw_calc(:)
    integer,          pointer :: num_cg_ele_calc(:)
    integer,          pointer :: cg_vdw_list(:,:), cg_ele_list(:,:)
    integer(1),       pointer :: exclusion_mask(:,:)

    ncell                =  domain%num_cell_local
    nboundary            =  domain%num_cell_boundary
    natom                => domain%num_atom
    start_atom           => domain%start_atom
    ncharge              => domain%num_charge
    atom_2_cell          => domain%atom_2_cell
    near_cells_count     => domain%near_cells_count
    far_cells_count      => domain%far_cells_count
    near_cells           => domain%near_cells
    far_cells            => domain%far_cells
    id_l2g               => domain%id_l2g

    exclusion_mask       => enefunc%exclusion_mask
    cg_elec_list_inv     => enefunc%cg_elec_list_inv

    num_cg_vdw_calc      => pairlist%num_cg_exv_calc
    num_cg_ele_calc      => pairlist%num_cg_ele_calc
    cg_vdw_list          => pairlist%cg_exv_list
    cg_ele_list          => pairlist%cg_ele_list

    pairdist2_ele        =  pairlist%cg_pairlistdist_ele**2
    pairdist2_vdw        =  pairlist%cg_pairlistdist_vdw**2

    num_atom = domain%num_atom_domain + domain%num_atom_boundary
    num_cg_vdw_calc(1:num_atom) = 0
    num_atom = enefunc%num_cg_elec
    num_cg_ele_calc(1:num_atom) = 0

    !$omp parallel default(shared)                                          &
    !$omp private(id, i, ixx, rtmp, num_nb15, iy, iyy, dij, rij2, inbc, j,  &
    !$omp         start_j, start_i, k)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do ixx = id+1, domain%num_atom_domain+domain%num_atom_boundary, nthread

      i = atom_2_cell(ixx)
      rtmp(1) = coord(ixx,1)
      rtmp(2) = coord(ixx,2)
      rtmp(3) = coord(ixx,3)

      num_nb15 = 0

      if (i <= ncell) then
 
        do iyy = ixx + 1, start_atom(i+1)
          if (exclusion_mask(iyy,ixx) /= 1) then
            dij(1) = rtmp(1) - coord(iyy,1)
            dij(2) = rtmp(2) - coord(iyy,2)
            dij(3) = rtmp(3) - coord(iyy,3)
            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
            if (rij2 < pairdist2_vdw) then
              num_nb15 = num_nb15 + 1
              cg_vdw_list(num_nb15,ixx) = iyy
            end if
          end if
        end do
      end if

      do inbc = 1, near_cells_count(i)

        j = near_cells(inbc,i)
        start_j = start_atom(j)
        do iy = 1, natom(j)
          iyy = iy + start_j
          if (exclusion_mask(iyy,ixx) /= 1) then
            dij(1) = rtmp(1) - coord(iyy,1)
            dij(2) = rtmp(2) - coord(iyy,2)
            dij(3) = rtmp(3) - coord(iyy,3)
            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
            if (rij2 < pairdist2_vdw) then
              num_nb15 = num_nb15 + 1
              cg_vdw_list(num_nb15,ixx) = iyy
            end if
          end if
        end do
      end do

      num_cg_vdw_calc(ixx) = num_nb15
      if (num_nb15 > Max_exv_nb15*0.7) pairlist%realloc_exv = 1

    end do

    ! electrostatic
    !
    do i = id+1, ncell+nboundary, nthread

      start_i = start_atom(i)

      do ix = 1, ncharge(i) 

        num_nb15 = 0
        ixx = ix + start_i
        k = cg_elec_list_inv(ixx)
        rtmp(1) = coord(ixx,1)
        rtmp(2) = coord(ixx,2)
        rtmp(3) = coord(ixx,3)

        if (i <= ncell) then
          do iy = ix + 1, ncharge(i)
            iyy = iy + start_i
            if (exclusion_mask(iyy,ixx) /= 1) then
              dij(1) = rtmp(1) - coord(iyy,1)
              dij(2) = rtmp(2) - coord(iyy,2)
              dij(3) = rtmp(3) - coord(iyy,3)
              rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
              if (rij2 < pairdist2_ele) then
                num_nb15 = num_nb15 + 1
                cg_ele_list(num_nb15,k) = iyy
              end if
            end if
          end do
        end if

        do inbc = 1, near_cells_count(i)
          j = near_cells(inbc,i)
          start_j = start_atom(j)
          do iy = 1, ncharge(j)
            iyy = iy + start_j
            if (exclusion_mask(iyy,ixx) /= 1) then
              dij(1) = rtmp(1) - coord(iyy,1)
              dij(2) = rtmp(2) - coord(iyy,2)
              dij(3) = rtmp(3) - coord(iyy,3)
              rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
              if (rij2 < pairdist2_ele) then
                num_nb15 = num_nb15 + 1
                cg_ele_list(num_nb15,k) = iyy
              end if
            end if
          end do
        end do

        num_cg_ele_calc(k) = num_nb15
        if (num_nb15 > Max_elec_nb15*0.7) pairlist%realloc_elec = 1

      end do
    end do

    !$omp end parallel

    return

  end subroutine update_pairlist_martini

  subroutine pairlist_exv_DNA_idr(k_dna, k_base, k_kh, k_idr_kh,             &
                         k_idr_hps, k_idr_mpipi,                             &
                         dna_check, cg_pro_use_KH, cg_IDR_KH,                &
                         cg_IDR_HPS, CG_IDR_MPIPI,                           &
                         chain_id, base_type, base_pair_is_WC,               &
                         id_l2g, cg_kh_mol_pair, ixx, iyy, i_atom,           &
                         chain_idx, i_idr, i_dna, i_base_type, i_kh,         &
                         i_idr_kh, i_idr_hps, i_idr_mpipi,                   &
                         num_nb15_base, num_nb15_dna,                        &
                         num_nb15, num_nb15_kh, num_nb15_idr_kh,             &
                         num_nb15_idr_hps, num_nb15_idr_mpipi,               &
                         cg_DNA_base_list,                                   &
                         cg_DNA_exv_list, cg_exv_list, cg_kh_list,           &
                         cg_idr_kh_list, cg_idr_hps_list, cg_idr_mpipi_list)

    integer,                  intent(in   ) :: k_dna, k_base, k_kh
    integer,                  intent(in   ) :: k_idr_kh, k_idr_hps, k_idr_mpipi
    integer(1),               intent(in   ) :: dna_check(:)
    integer(1),               intent(in   ) :: cg_pro_use_KH(:)
    integer(1),               intent(in   ) :: cg_IDR_KH(:), cg_IDR_HPS(:)
    integer(1),               intent(in   ) :: cg_IDR_MPIPI(:)
    integer,                  intent(in   ) :: chain_id(:), base_type(:)
    logical,                  intent(in   ) :: base_pair_is_WC(:,:)
    integer,                  intent(in   ) :: id_l2g(:)
    integer(1),               intent(in   ) :: cg_kh_mol_pair(:,:)
    integer,                  intent(in   ) :: chain_idx, i_atom, ixx, iyy
    integer,                  intent(in   ) :: i_idr, i_dna, i_kh, i_base_type
    integer,                  intent(in   ) :: i_idr_kh, i_idr_hps, i_idr_mpipi
    integer,                  intent(inout) :: num_nb15_base, num_nb15_dna
    integer,                  intent(inout) :: num_nb15, num_nb15_kh
    integer,                  intent(inout) :: num_nb15_idr_kh, num_nb15_idr_hps
    integer,                  intent(inout) :: num_nb15_idr_mpipi
    integer,                  intent(inout) :: cg_DNA_base_list(:,:)
    integer,                  intent(inout) :: cg_DNA_exv_list(:,:)
    integer,                  intent(inout) :: cg_exv_list(:,:)
    integer,                  intent(inout) :: cg_kh_list(:,:)
    integer,                  intent(inout) :: cg_idr_kh_list(:,:)
    integer,                  intent(inout) :: cg_idr_hps_list(:,:)
    integer,                  intent(inout) :: cg_idr_mpipi_list(:,:)
    
    integer                   :: j_atom, chain_idy
    integer                   :: j_dna, j_idr, j_kh
    integer                   :: j_idr_kh, j_idr_hps, j_idr_mpipi
    integer                   :: ij_dna, ij_kh
    integer                   :: ij_idr_kh, ij_idr_hps, ij_idr_mpipi
    integer                   :: ij_idr, ij_idr1
    integer                   :: j_base_type
    logical                   :: is_nonlocal, is_WC_bp

    j_dna      = dna_check(iyy)
    j_kh       = cg_pro_use_KH(iyy)
    j_idr_kh   = cg_IDR_KH(iyy)
    j_idr_hps  = cg_IDR_HPS(iyy)
    j_idr_mpipi= cg_IDR_MPIPI(iyy)
    j_idr      = j_idr_kh + j_idr_hps + j_idr_mpipi
    ij_dna     = i_dna * j_dna
    ij_kh      = i_kh * j_kh
    ij_idr_kh  = i_idr_kh * j_idr_kh
    ij_idr_hps = i_idr_hps * j_idr_hps
    ij_idr_mpipi = i_idr_mpipi * j_idr_mpipi
    ij_idr     = i_idr * j_idr
    ij_idr1    = i_idr + j_idr

    if (ij_dna /= 0) then

      is_nonlocal = .false.
      is_WC_bp    = .false.
      chain_idy = chain_id(iyy)
      j_base_type = base_type(iyy)
      if (chain_idx /= chain_idy) then
        is_nonlocal = .true.
        if (ij_dna == 1) is_WC_bp = base_pair_is_WC(i_base_type,j_base_type)
      else
        j_atom = id_l2g(iyy)
        if (abs(j_atom-i_atom) > 4) is_nonlocal = .true.
        if (abs(j_atom-i_atom) > 10) then
          if (ij_dna == 1) &
            is_WC_bp = base_pair_is_WC(i_base_type,j_base_type)
        end if
      end if
      if (is_WC_bp) then
        num_nb15_base = num_nb15_base + 1
        cg_DNA_base_list(num_nb15_base,k_base) = iyy
      else if (is_nonlocal .and. .not.is_WC_bp) then
        num_nb15_dna = num_nb15_dna + 1
        cg_DNA_exv_list(num_nb15_dna,k_dna) = iyy
      end if

    else

      chain_idy = chain_id(iyy)
      j_atom = id_l2g(iyy)
      if (ij_kh /= 0) then
        if (cg_kh_mol_pair(chain_idx,chain_idy)==0) ij_kh = 0
        if (ij_idr /= 0) ij_kh = 0
        if (ij_kh /= 0) then
          if (chain_idx == chain_idy) then
            if (ij_idr1 == 0) ij_kh = 0
          end if
        end if
      end if

      if (ij_kh == 0 .and. ij_idr_kh == 0 .and. &
          ij_idr_hps == 0 .and. ij_idr_mpipi == 0) then
        num_nb15 = num_nb15 + 1
        cg_exv_list(num_nb15,ixx) = iyy
      end if
      if (ij_kh /= 0) then
        if (.not. (chain_idx == chain_idy .and. abs(j_atom-i_atom)<=3)) then
          num_nb15_kh = num_nb15_kh + 1
          cg_kh_list(num_nb15_kh,k_kh) = iyy
        end if
      end if
      if (ij_idr_kh /= 0) then
        if (.not. (chain_idx == chain_idy .and. abs(i_atom-j_atom)==1)) then
          num_nb15_idr_kh = num_nb15_idr_kh + 1
          cg_idr_kh_list(num_nb15_idr_kh,k_idr_kh) = iyy
        end if
      end if
      if (ij_idr_hps /= 0) then
        if (.not. (chain_idx == chain_idy .and. abs(i_atom-j_atom)==1)) then
          num_nb15_idr_hps = num_nb15_idr_hps + 1
          cg_idr_hps_list(num_nb15_idr_hps,k_idr_hps) = iyy
        end if
      end if
      if (ij_idr_mpipi /= 0) then
        if (.not. (chain_idx == chain_idy .and. abs(i_atom-j_atom)==1)) then
          num_nb15_idr_mpipi = num_nb15_idr_mpipi + 1
          cg_idr_mpipi_list(num_nb15_idr_mpipi,k_idr_mpipi) = iyy
        end if
      end if

    end if

    return

  end subroutine pairlist_exv_DNA_idr

  subroutine pairlist_DNA_idr(k_base, k_kh, k_idr_kh, k_idr_hps, k_idr_mpipi, &
                          dna_check,                                          &
                          cg_pro_use_KH, cg_IDR_KH, cg_IDR_HPS, cg_IDR_MPIPI, &
                          chain_id, base_type, base_pair_is_WC, id_l2g,       &
                          cg_kh_mol_pair, iyy, i_atom, chain_idx, i_idr,      &
                          i_dna, i_base_type, i_kh, i_idr_kh, i_idr_hps,      &
                          i_idr_mpipi,                                        &
                          num_nb15_base, num_nb15_kh, num_nb15_idr_kh,        &
                          num_nb15_idr_hps, num_nb15_idr_mpipi,               &
                          cg_DNA_base_list, cg_kh_list,                       &
                          cg_idr_kh_list, cg_idr_hps_list, cg_idr_mpipi_list)

    integer,                  intent(in   ) :: k_base, k_kh
    integer,                  intent(in   ) :: k_idr_kh, k_idr_hps, k_idr_mpipi
    integer(1),               intent(in   ) :: dna_check(:)
    integer(1),               intent(in   ) :: cg_pro_use_KH(:)
    integer(1),               intent(in   ) :: cg_IDR_KH(:), cg_IDR_HPS(:)
    integer(1),               intent(in   ) :: cg_IDR_MPIPI(:)
    integer,                  intent(in   ) :: chain_id(:), base_type(:)
    logical,                  intent(in   ) :: base_pair_is_WC(:,:)
    integer,                  intent(in   ) :: id_l2g(:)
    integer(1),               intent(in   ) :: cg_kh_mol_pair(:,:)
    integer,                  intent(in   ) :: chain_idx, i_atom, iyy
    integer,                  intent(in   ) :: i_idr, i_dna, i_kh, i_base_type
    integer,                  intent(in   ) :: i_idr_kh, i_idr_hps, i_idr_mpipi
    integer,                  intent(inout) :: num_nb15_base
    integer,                  intent(inout) :: num_nb15_kh
    integer,                  intent(inout) :: num_nb15_idr_kh, num_nb15_idr_hps
    integer,                  intent(inout) :: num_nb15_idr_mpipi
    integer,                  intent(inout) :: cg_DNA_base_list(:,:)
    integer,                  intent(inout) :: cg_kh_list(:,:)
    integer,                  intent(inout) :: cg_idr_kh_list(:,:)
    integer,                  intent(inout) :: cg_idr_hps_list(:,:)
    integer,                  intent(inout) :: cg_idr_mpipi_list(:,:)

    integer                   :: j_atom, chain_idy
    integer                   :: j_dna, j_idr, j_kh
    integer                   :: j_idr_kh, j_idr_hps, j_idr_mpipi
    integer                   :: ij_dna, ij_kh
    integer                   :: ij_idr_kh, ij_idr_hps, ij_idr_mpipi
    integer                   :: ij_idr, ij_idr1
    integer                   :: j_base_type
    logical                   :: is_nonlocal, is_WC_bp

    j_dna      = dna_check(iyy)
    j_kh       = cg_pro_use_KH(iyy)
    j_idr_kh   = cg_IDR_KH(iyy)
    j_idr_hps  = cg_IDR_HPS(iyy)
    j_idr_mpipi = cg_IDR_MPIPI(iyy)
    j_idr      = j_idr_kh + j_idr_hps + j_idr_mpipi
    ij_dna     = i_dna * j_dna
    ij_kh      = i_kh * j_kh
    ij_idr_kh  = i_idr_kh * j_idr_kh
    ij_idr_hps = i_idr_hps * j_idr_hps
    ij_idr_mpipi = i_idr_mpipi * j_idr_mpipi
    ij_idr     = i_idr * j_idr
    ij_idr1    = i_idr + j_idr

    if (ij_dna /= 0) then

      is_nonlocal = .false.
      is_WC_bp    = .false.
      chain_idy = chain_id(iyy)
      j_base_type = base_type(iyy)
      if (chain_idx /= chain_idy) then
        is_nonlocal = .true.
        if (ij_dna == 1) is_WC_bp = base_pair_is_WC(i_base_type,j_base_type)
      else
        j_atom = id_l2g(iyy)
        if (abs(j_atom-i_atom) > 4) is_nonlocal = .true.
        if (abs(j_atom-i_atom) > 10) then
          if (ij_dna == 1) &
            is_WC_bp = base_pair_is_WC(i_base_type,j_base_type)
        end if
      end if
      if (is_WC_bp) then
        num_nb15_base = num_nb15_base + 1
        cg_DNA_base_list(num_nb15_base,k_base) = iyy
      end if

    else

      chain_idy = chain_id(iyy)
      j_atom = id_l2g(iyy)
      if (ij_kh /= 0) then
        if (cg_kh_mol_pair(chain_idx,chain_idy)==0) ij_kh = 0
        if (ij_idr /= 0) ij_kh = 0
        if (ij_kh /= 0) then
          if (chain_idx == chain_idy) then
            if (ij_idr1 == 0) ij_kh = 0
          end if
        end if
      end if

      if (ij_kh /= 0) then
        if (.not. (chain_idx == chain_idy .and. abs(j_atom-i_atom)<=3)) then
          num_nb15_kh = num_nb15_kh + 1
          cg_kh_list(num_nb15_kh,k_kh) = iyy
        end if
      end if
      if (ij_idr_kh /= 0) then
        if (.not. (chain_idx == chain_idy .and. abs(i_atom-j_atom)==1)) then
          num_nb15_idr_kh = num_nb15_idr_kh + 1
          cg_idr_kh_list(num_nb15_idr_kh,k_idr_kh) = iyy
        end if
      end if
      if (ij_idr_hps /= 0) then
        if (.not. (chain_idx == chain_idy .and. abs(i_atom-j_atom)==1)) then
          num_nb15_idr_hps = num_nb15_idr_hps + 1
          cg_idr_hps_list(num_nb15_idr_hps,k_idr_hps) = iyy
        end if
      end if
      if (ij_idr_mpipi /= 0) then
        if (.not. (chain_idx == chain_idy .and. abs(i_atom-j_atom)==1)) then
          num_nb15_idr_mpipi = num_nb15_idr_mpipi + 1
          cg_idr_mpipi_list(num_nb15_idr_mpipi,k_idr_mpipi) = iyy
        end if
      end if

    end if

    return

  end subroutine pairlist_DNA_idr

  subroutine pairlist_idr(k_kh, k_idr_kh, k_idr_hps, k_idr_mpipi,        &
                          cg_pro_use_KH,                                 &
                          cg_IDR_KH, cg_IDR_HPS, cg_IDR_MPIPI,           &
                          chain_id, id_l2g,                              &
                          cg_kh_mol_pair, iyy, i_atom, chain_idx, i_idr, &
                          i_kh, i_idr_kh, i_idr_hps, i_idr_mpipi,        &
                          num_nb15_kh,                                   &
                          num_nb15_idr_kh, num_nb15_idr_hps,             &
                          num_nb15_idr_mpipi, cg_kh_list,                &
                          cg_idr_kh_list, cg_idr_hps_list,               &
                          cg_idr_mpipi_list)

    integer,                  intent(in   ) :: k_kh, k_idr_kh, k_idr_hps
    integer,                  intent(in   ) :: k_idr_mpipi
    integer(1),               intent(in   ) :: cg_pro_use_KH(:)
    integer(1),               intent(in   ) :: cg_IDR_KH(:), cg_IDR_HPS(:)
    integer(1),               intent(in   ) :: cg_IDR_MPIPI(:)
    integer,                  intent(in   ) :: chain_id(:)
    integer,                  intent(in   ) :: id_l2g(:)
    integer(1),               intent(in   ) :: cg_kh_mol_pair(:,:)
    integer,                  intent(in   ) :: chain_idx, i_atom, iyy
    integer,                  intent(in   ) :: i_idr, i_kh
    integer,                  intent(in   ) :: i_idr_kh, i_idr_hps, i_idr_mpipi
    integer,                  intent(inout) :: num_nb15_kh
    integer,                  intent(inout) :: num_nb15_idr_kh, num_nb15_idr_hps
    integer,                  intent(inout) :: num_nb15_idr_mpipi
    integer,                  intent(inout) :: cg_kh_list(:,:)
    integer,                  intent(inout) :: cg_idr_kh_list(:,:)
    integer,                  intent(inout) :: cg_idr_hps_list(:,:)
    integer,                  intent(inout) :: cg_idr_mpipi_list(:,:)

    integer                   :: j_atom, chain_idy
    integer                   :: j_idr, j_kh
    integer                   :: j_idr_kh, j_idr_hps, j_idr_mpipi
    integer                   :: ij_kh, ij_idr_kh, ij_idr_hps, ij_idr_mpipi
    integer                   :: ij_idr, ij_idr1

    j_kh       = cg_pro_use_KH(iyy)
    j_idr_kh   = cg_IDR_KH(iyy)
    j_idr_hps  = cg_IDR_HPS(iyy)
    j_idr_mpipi  = cg_IDR_MPIPI(iyy)
    j_idr      = j_idr_kh + j_idr_hps + j_idr_mpipi
    ij_kh      = i_kh * j_kh
    ij_idr_kh  = i_idr_kh * j_idr_kh
    ij_idr_hps = i_idr_hps * j_idr_hps
    ij_idr_mpipi = i_idr_mpipi * j_idr_mpipi
    ij_idr     = i_idr * j_idr
    ij_idr1    = i_idr + j_idr

    chain_idy = chain_id(iyy)
    j_atom = id_l2g(iyy)
    if (ij_kh /= 0) then
      if (cg_kh_mol_pair(chain_idx,chain_idy)==0) ij_kh = 0
      if (ij_idr /= 0) ij_kh = 0
      if (ij_kh /= 0) then
        if (chain_idx == chain_idy) then
          if (ij_idr1 == 0) ij_kh = 0
        end if
      end if
    end if

    if (ij_kh /= 0) then
      if (.not. (chain_idx == chain_idy .and. abs(j_atom-i_atom)<=3)) then
        num_nb15_kh = num_nb15_kh + 1
        cg_kh_list(num_nb15_kh,k_kh) = iyy
      end if
    end if
    if (ij_idr_kh /= 0) then
      if (.not. (chain_idx == chain_idy .and. abs(i_atom-j_atom)==1)) then
        num_nb15_idr_kh = num_nb15_idr_kh + 1
        cg_idr_kh_list(num_nb15_idr_kh,k_idr_kh) = iyy
      end if
    end if
    if (ij_idr_hps /= 0) then
      if (.not. (chain_idx == chain_idy .and. abs(i_atom-j_atom)==1)) then
        num_nb15_idr_hps = num_nb15_idr_hps + 1
        cg_idr_hps_list(num_nb15_idr_hps,k_idr_hps) = iyy
      end if
    end if
    if (ij_idr_mpipi /= 0) then
      if (.not. (chain_idx == chain_idy .and. abs(i_atom-j_atom)==1)) then
        num_nb15_idr_mpipi = num_nb15_idr_mpipi + 1
        cg_idr_mpipi_list(num_nb15_idr_mpipi,k_idr_mpipi) = iyy
      end if
    end if

    return

  end subroutine pairlist_idr

  subroutine pairlist_exv_DNA_idr_check(dna_check, cg_pro_use_KH, cg_IDR_KH, &
                                cg_IDR_HPS, cg_IDR_MPIPI,                    &
                                chain_id, base_type,                         &
                                base_pair_is_WC, id_l2g, cg_kh_mol_pair,     &
                                iyy, i_atom, chain_idx, i_idr, i_dna,        &
                                i_base_type, i_kh, i_idr_kh,                 &
                                i_idr_hps, i_idr_mpipi,                      &
                                num_nb15_base, num_nb15_dna, num_nb15,       &
                                num_nb15_kh, num_nb15_idr_kh,                &
                                num_nb15_idr_hps, num_nb15_idr_mpipi)

    integer(1),               intent(in   ) :: dna_check(:)
    integer(1),               intent(in   ) :: cg_pro_use_KH(:)
    integer(1),               intent(in   ) :: cg_IDR_KH(:), cg_IDR_HPS(:)
    integer(1),               intent(in   ) :: cg_IDR_MPIPI(:)
    integer,                  intent(in   ) :: chain_id(:), base_type(:)
    logical,                  intent(in   ) :: base_pair_is_WC(:,:)
    integer,                  intent(in   ) :: id_l2g(:)
    integer(1),               intent(in   ) :: cg_kh_mol_pair(:,:)
    integer,                  intent(in   ) :: chain_idx, i_atom, iyy
    integer,                  intent(in   ) :: i_idr, i_dna, i_kh, i_base_type
    integer,                  intent(in   ) :: i_idr_kh, i_idr_hps, i_idr_mpipi
    integer,                  intent(inout) :: num_nb15_base, num_nb15_dna
    integer,                  intent(inout) :: num_nb15, num_nb15_kh
    integer,                  intent(inout) :: num_nb15_idr_kh, num_nb15_idr_hps
    integer,                  intent(inout) :: num_nb15_idr_mpipi

    integer                   :: j_atom, chain_idy
    integer                   :: j_dna, j_idr, j_kh
    integer                   :: j_idr_kh, j_idr_hps, j_idr_mpipi
    integer                   :: ij_dna, ij_kh
    integer                   :: ij_idr_kh, ij_idr_hps, ij_idr_mpipi
    integer                   :: ij_idr, ij_idr1
    integer                   :: j_base_type
    logical                   :: is_nonlocal, is_WC_bp

    j_dna      = dna_check(iyy)
    j_kh       = cg_pro_use_KH(iyy)
    j_idr_kh   = cg_IDR_KH(iyy)
    j_idr_hps  = cg_IDR_HPS(iyy)
    j_idr_mpipi  = cg_IDR_MPIPI(iyy)
    j_idr      = j_idr_kh + j_idr_hps + j_idr_mpipi
    ij_dna     = i_dna * j_dna
    ij_kh      = i_kh * j_kh
    ij_idr_kh  = i_idr_kh * j_idr_kh
    ij_idr_hps = i_idr_hps * j_idr_hps
    ij_idr_mpipi = i_idr_mpipi * j_idr_mpipi
    ij_idr     = i_idr * j_idr
    ij_idr1    = i_idr + j_idr

    if (ij_dna /= 0) then

      is_nonlocal = .false.
      is_WC_bp    = .false.
      chain_idy = chain_id(iyy)
      j_base_type = base_type(iyy)
      if (chain_idx /= chain_idy) then
        is_nonlocal = .true.
        if (ij_dna == 1) is_WC_bp = base_pair_is_WC(i_base_type,j_base_type)
      else
        j_atom = id_l2g(iyy)
        if (abs(j_atom-i_atom) > 4) is_nonlocal = .true.
        if (abs(j_atom-i_atom) > 10) then
          if (ij_dna == 1) &
            is_WC_bp = base_pair_is_WC(i_base_type,j_base_type)
        end if
      end if
      if (is_WC_bp) then
        num_nb15_base = num_nb15_base + 1
      else if (is_nonlocal .and. .not.is_WC_bp) then
        num_nb15_dna = num_nb15_dna + 1
      end if

    else

      chain_idy = chain_id(iyy)
      j_atom = id_l2g(iyy)
      if (ij_kh /= 0) then
        if (cg_kh_mol_pair(chain_idx,chain_idy)==0) ij_kh = 0
        if (ij_idr /= 0) ij_kh = 0
        if (ij_kh /= 0) then
          if (chain_idx == chain_idy) then
            if (ij_idr1 == 0) ij_kh = 0
          end if
        end if
      end if

      if (ij_kh == 0 .and. ij_idr_kh == 0 &
        .and. ij_idr_hps == 0 .and. ij_idr_mpipi == 0) then
        num_nb15 = num_nb15 + 1
      end if
      if (ij_kh /= 0) then
        if (.not. (chain_idx == chain_idy .and. abs(j_atom-i_atom)<=3)) then
          num_nb15_kh = num_nb15_kh + 1
        end if
      end if
      if (ij_idr_kh /= 0) then
        if (.not. (chain_idx == chain_idy .and. abs(i_atom-j_atom)==1)) then
          num_nb15_idr_kh = num_nb15_idr_kh + 1
        end if
      end if
      if (ij_idr_hps /= 0) then
        if (.not. (chain_idx == chain_idy .and. abs(i_atom-j_atom)==1)) then
          num_nb15_idr_hps = num_nb15_idr_hps + 1
        end if
      end if
      if (ij_idr_mpipi /= 0) then
        if (.not. (chain_idx == chain_idy .and. abs(i_atom-j_atom)==1)) then
          num_nb15_idr_mpipi = num_nb15_idr_mpipi + 1
        end if
      end if

    end if

    return

  end subroutine pairlist_exv_DNA_idr_check

  subroutine pairlist_DNA_idr_check(dna_check, cg_pro_use_KH, cg_IDR_KH,     &
                                cg_IDR_HPS, cg_IDR_MPIPI,                    &
                                chain_id, base_type,                         &
                                base_pair_is_WC, id_l2g, cg_kh_mol_pair,     &
                                iyy, i_atom, chain_idx, i_idr, i_dna,        &
                                i_base_type, i_kh,                           &
                                i_idr_kh, i_idr_hps, i_idr_mpipi,            &
                                num_nb15_base, num_nb15_kh, num_nb15_idr_kh, &
                                num_nb15_idr_hps, num_nb15_idr_mpipi)

    integer(1),               intent(in   ) :: dna_check(:)
    integer(1),               intent(in   ) :: cg_pro_use_KH(:)
    integer(1),               intent(in   ) :: cg_IDR_KH(:)
    integer(1),               intent(in   ) :: cg_IDR_HPS(:), cg_IDR_MPIPI(:)
    integer,                  intent(in   ) :: chain_id(:), base_type(:)
    logical,                  intent(in   ) :: base_pair_is_WC(:,:)
    integer,                  intent(in   ) :: id_l2g(:)
    integer(1),               intent(in   ) :: cg_kh_mol_pair(:,:)
    integer,                  intent(in   ) :: chain_idx, i_atom, iyy
    integer,                  intent(in   ) :: i_idr, i_dna, i_kh, i_base_type
    integer,                  intent(in   ) :: i_idr_kh, i_idr_hps, i_idr_mpipi
    integer,                  intent(inout) :: num_nb15_base
    integer,                  intent(inout) :: num_nb15_kh
    integer,                  intent(inout) :: num_nb15_idr_kh, num_nb15_idr_hps
    integer,                  intent(inout) :: num_nb15_idr_mpipi

    integer                   :: j_atom, chain_idy
    integer                   :: j_dna, j_idr, j_kh
    integer                   :: j_idr_kh, j_idr_hps, j_idr_mpipi
    integer                   :: ij_dna, ij_kh
    integer                   :: ij_idr_kh, ij_idr_hps, ij_idr_mpipi
    integer                   :: ij_idr, ij_idr1
    integer                   :: j_base_type
    logical                   :: is_nonlocal, is_WC_bp

    j_dna      = dna_check(iyy)
    j_kh       = cg_pro_use_KH(iyy)
    j_idr_kh   = cg_IDR_KH(iyy)
    j_idr_hps  = cg_IDR_HPS(iyy)
    j_idr_mpipi  = cg_IDR_MPIPI(iyy)
    j_idr      = j_idr_kh + j_idr_hps + j_idr_mpipi
    ij_dna     = i_dna * j_dna
    ij_kh      = i_kh * j_kh
    ij_idr_kh  = i_idr_kh * j_idr_kh
    ij_idr_hps = i_idr_hps * j_idr_hps
    ij_idr_mpipi = i_idr_mpipi * j_idr_mpipi
    ij_idr     = i_idr * j_idr
    ij_idr1    = i_idr + j_idr

    if (ij_dna /= 0) then

      is_nonlocal = .false.
      is_WC_bp    = .false.
      chain_idy = chain_id(iyy)
      j_base_type = base_type(iyy)
      if (chain_idx /= chain_idy) then
        is_nonlocal = .true.
        if (ij_dna == 1) is_WC_bp = base_pair_is_WC(i_base_type,j_base_type)
      else
        j_atom = id_l2g(iyy)
        if (abs(j_atom-i_atom) > 4) is_nonlocal = .true.
        if (abs(j_atom-i_atom) > 10) then
          if (ij_dna == 1) &
            is_WC_bp = base_pair_is_WC(i_base_type,j_base_type)
        end if
      end if
      if (is_WC_bp) then
        num_nb15_base = num_nb15_base + 1
      end if

    else

      chain_idy = chain_id(iyy)
      j_atom = id_l2g(iyy)
      if (ij_kh /= 0) then
        if (cg_kh_mol_pair(chain_idx,chain_idy)==0) ij_kh = 0
        if (ij_idr /= 0) ij_kh = 0
        if (ij_kh /= 0) then
          if (chain_idx == chain_idy) then
            if (ij_idr1 == 0) ij_kh = 0
          end if
        end if
      end if

      if (ij_kh /= 0) then
        if (.not. (chain_idx == chain_idy .and. abs(j_atom-i_atom)<=3)) then
          num_nb15_kh = num_nb15_kh + 1
        end if
      end if
      if (ij_idr_kh /= 0) then
        if (.not. (chain_idx == chain_idy .and. abs(i_atom-j_atom)==1)) then
          num_nb15_idr_kh = num_nb15_idr_kh + 1
        end if
      end if
      if (ij_idr_hps /= 0) then
        if (.not. (chain_idx == chain_idy .and. abs(i_atom-j_atom)==1)) then
          num_nb15_idr_hps = num_nb15_idr_hps + 1
        end if
      end if
      if (ij_idr_mpipi /= 0) then
        if (.not. (chain_idx == chain_idy .and. abs(i_atom-j_atom)==1)) then
          num_nb15_idr_mpipi = num_nb15_idr_mpipi + 1
        end if
      end if

    end if

    return

  end subroutine pairlist_DNA_idr_check

  subroutine pairlist_idr_check(cg_pro_use_KH, cg_IDR_KH,                 &
                                cg_IDR_HPS, cg_IDR_MPIPI,                 &
                                chain_id, id_l2g, cg_kh_mol_pair, iyy,    &
                                i_atom, chain_idx, i_idr, i_kh, i_idr_kh, &
                                i_idr_hps, i_idr_mpipi,                   &
                                num_nb15_kh, num_nb15_idr_kh,             &
                                num_nb15_idr_hps, num_nb15_idr_mpipi)

    integer(1),               intent(in   ) :: cg_pro_use_KH(:)
    integer(1),               intent(in   ) :: cg_IDR_KH(:), cg_IDR_HPS(:)
    integer(1),               intent(in   ) :: cg_IDR_MPIPI(:)
    integer,                  intent(in   ) :: chain_id(:)
    integer,                  intent(in   ) :: id_l2g(:)
    integer(1),               intent(in   ) :: cg_kh_mol_pair(:,:)
    integer,                  intent(in   ) :: chain_idx, i_atom, iyy
    integer,                  intent(in   ) :: i_idr, i_kh
    integer,                  intent(in   ) :: i_idr_kh, i_idr_hps, i_idr_mpipi
    integer,                  intent(inout) :: num_nb15_kh
    integer,                  intent(inout) :: num_nb15_idr_kh, num_nb15_idr_hps
    integer,                  intent(inout) :: num_nb15_idr_mpipi

    integer                   :: j_atom, chain_idy
    integer                   :: j_idr, j_kh
    integer                   :: j_idr_kh, j_idr_hps, j_idr_mpipi
    integer                   :: ij_kh, ij_idr_kh, ij_idr_hps, ij_idr_mpipi
    integer                   :: ij_idr, ij_idr1

    j_kh       = cg_pro_use_KH(iyy)
    j_idr_kh   = cg_IDR_KH(iyy)
    j_idr_hps  = cg_IDR_HPS(iyy)
    j_idr_mpipi  = cg_IDR_MPIPI(iyy)
    j_idr      = j_idr_kh + j_idr_hps + j_idr_mpipi
    ij_kh      = i_kh * j_kh
    ij_idr_kh  = i_idr_kh * j_idr_kh
    ij_idr_hps = i_idr_hps * j_idr_hps
    ij_idr_mpipi = i_idr_mpipi * j_idr_mpipi
    ij_idr     = i_idr * j_idr
    ij_idr1    = i_idr + j_idr

    chain_idy = chain_id(iyy)
    j_atom = id_l2g(iyy)
    if (ij_kh /= 0) then
      if (cg_kh_mol_pair(chain_idx,chain_idy)==0) ij_kh = 0
      if (ij_idr /= 0) ij_kh = 0
      if (ij_kh /= 0) then
        if (chain_idx == chain_idy) then
          if (ij_idr1 == 0) ij_kh = 0
        end if
      end if
    end if

    if (ij_kh /= 0) then
      if (.not. (chain_idx == chain_idy .and. abs(j_atom-i_atom)<=3)) then
        num_nb15_kh = num_nb15_kh + 1
      end if
    end if
    if (ij_idr_kh /= 0) then
      if (.not. (chain_idx == chain_idy .and. abs(i_atom-j_atom)==1)) then
        num_nb15_idr_kh = num_nb15_idr_kh + 1
      end if
    end if
    if (ij_idr_hps /= 0) then
      if (.not. (chain_idx == chain_idy .and. abs(i_atom-j_atom)==1)) then
        num_nb15_idr_hps = num_nb15_idr_hps + 1
      end if
    end if
    if (ij_idr_mpipi /= 0) then
      if (.not. (chain_idx == chain_idy .and. abs(i_atom-j_atom)==1)) then
        num_nb15_idr_mpipi = num_nb15_idr_mpipi + 1
      end if
    end if

    return

  end subroutine pairlist_idr_check

end module cg_pairlist_mod
