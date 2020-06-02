!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_prmtop_mod
!> @brief   read AMBER prmtop file
!! @authors Yusuke Sugihara FJ (YSFJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fileio_prmtop_mod

  use fileio_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  !
  type, public :: s_prmtop

    ! title
    character(80)    :: version
    character(80)    :: title

    ! atom
    character(4),     allocatable :: atom_name(:)
    character(4),     allocatable :: atom_cls_name(:)
    integer,          allocatable :: atom_cls_no(:)
    real(wp),         allocatable :: charge(:)
    real(wp),         allocatable :: mass(:)

    ! numlp

    ! amber prmtop pointers information
    integer          :: num_atoms          ! NATOM in AMBER
    integer          :: num_types          ! NTYPES 
    integer          :: num_bondh          ! NBONH 
    integer          :: num_mbonda         ! MBONA 
    integer          :: num_anglh          ! NTHETH 
    integer          :: num_mangla         ! MTHETA 
    integer          :: num_diheh          ! NPHIH 
    integer          :: num_mdihea         ! MPHIA 
    integer          :: num_hparm          ! NHPARM (currently not used)
    integer          :: num_parm           ! NPARM 
    integer          :: num_nb             ! NNB 
    integer          :: num_residues       ! NRES 
    integer          :: num_bona           ! NBONA 
    integer          :: num_theta          ! NTHETA 
    integer          :: num_phia           ! NPHIA 
    integer          :: num_uniqbond       ! NUMBND 
    integer          :: num_uniqangl       ! NUMANG
    integer          :: num_uniqdihe       ! NPTRA
    integer          :: num_types_prm      ! NATYP
    integer          :: num_bondtypes      ! NPHB
    integer          :: ifpert             ! IFPERT
    integer          :: num_bond_perturbed ! NBPER
    integer          :: num_angl_perturbed ! NGPER
    integer          :: num_dihe_perturbed ! NDPER
    integer          :: num_mbper          ! MBPER
    integer          :: num_mgper          ! MGPER
    integer          :: num_mdper          ! MDPER
    integer          :: ifbox              ! IFBOX
    integer          :: num_max_res_atoms  ! NMXRS
    integer          :: ifcap              ! IFCAP
    integer          :: num_extra_points   ! NUMEXTRA
    integer          :: num_copy           ! NCOPY
    integer          :: iptres             ! IPTRES
    integer          :: nspm               ! NSPM
    integer          :: nspsol             ! NSPSOL
    real(wp)         :: oldbeta            ! OLDBETA
    real(wp)         :: box(3)             ! BOX

    ! amber prmtop other information (working variables)
    integer,         allocatable :: excl_atom(:)
    integer,         allocatable :: nb_par_idx(:)
    character(4),    allocatable :: res_label(:)
    integer,         allocatable :: res_point(:)
    real(wp),        allocatable :: bond_fcons_uniq(:)
    real(wp),        allocatable :: bond_equil_uniq(:)
    real(wp),        allocatable :: angl_fcons_uniq(:)
    real(wp),        allocatable :: angl_equil_uniq(:)
    real(wp),        allocatable :: dihe_fcons_uniq(:)
    real(wp),        allocatable :: dihe_perio_uniq(:)
    real(wp),        allocatable :: dihe_phase_uniq(:)
    real(wp),        allocatable :: scee_scale_fact(:)
    real(wp),        allocatable :: scnb_scale_fact(:)
    real(wp),        allocatable :: solty(:)
    real(wp),        allocatable :: lennarda(:)
    real(wp),        allocatable :: lennardb(:)
    integer,         allocatable :: bond_inc_hy(:,:)
    integer,         allocatable :: bond_wo_hy(:,:)
    integer,         allocatable :: angl_inc_hy(:,:)
    integer,         allocatable :: angl_wo_hy(:,:)
    integer,         allocatable :: dihe_inc_hy(:,:)
    integer,         allocatable :: dihe_wo_hy(:,:)

    integer,         allocatable :: inb(:)

    real(wp),        allocatable :: hbond_acoef(:)
    real(wp),        allocatable :: hbond_bcoef(:)
    real(wp),        allocatable :: hb_cut(:)

    character(4),    allocatable :: classify(:)
    integer,         allocatable :: join_array(:)
    integer,         allocatable :: irotate(:)

    real(wp),        allocatable :: radi_born(:)
    real(wp),        allocatable :: fs_born(:)
    integer,        allocatable :: nsp(:)

    character(80)    :: radius_set

    ! flag for read data or not
    logical          :: lversion                    = .false.
    logical          :: ltitle                      = .false.
    logical          :: lpointers                   = .false.
    logical          :: latom_name                  = .false.
    logical          :: lcharge                     = .false.
    logical          :: lmass                       = .false.
    logical          :: latom_type_index            = .false.
    logical          :: lnumber_excluded_atoms      = .false.
    logical          :: lnonbonded_parm_index       = .false.
    logical          :: lresidue_label              = .false.
    logical          :: lresidue_pointer            = .false.
    logical          :: lbond_force_constant        = .false.
    logical          :: lbond_equil_value           = .false.
    logical          :: langle_force_constant       = .false.
    logical          :: langle_equil_value          = .false.
    logical          :: ldihedral_force_constant    = .false.
    logical          :: ldihedral_periodicity       = .false.
    logical          :: ldihedral_phase             = .false.
    logical          :: lscee_scale_factor          = .false.
    logical          :: lscnb_scale_factor          = .false.
    logical          :: lsolty                      = .false.
    logical          :: llennard_jones_acoef        = .false.
    logical          :: llennard_jones_bcoef        = .false.
    logical          :: lbonds_inc_hydrogen         = .false.
    logical          :: lbonds_without_hydrogen     = .false.
    logical          :: langles_inc_hydrogen        = .false.
    logical          :: langles_without_hydrogen    = .false.
    logical          :: ldihedrals_inc_hydrogen     = .false.
    logical          :: ldihedrals_without_hydrogen = .false.
    logical          :: lexcluded_atoms_list        = .false.
    logical          :: lhbond_acoef                = .false.
    logical          :: lhbond_bcoef                = .false.
    logical          :: lhbcut                      = .false.
    logical          :: lamber_atom_type            = .false.
    logical          :: ltree_chain_classification  = .false.
    logical          :: ljoin_array                 = .false.
    logical          :: lirotat                     = .false.
    logical          :: lradius_set                 = .false.
    logical          :: lradii                      = .false.
    logical          :: lscreen                     = .false.
    logical          :: lsolvent_pointers           = .false.
    logical          :: latoms_per_molecule         = .false.
    logical          :: lbox_dimensions             = .false.
    logical          :: lcap_info                   = .false.
    logical          :: lcap_info2                  = .false.
    logical          :: lpolarizability             = .false.

  end type s_prmtop

  ! parameters for allocatable variables
  integer, public, parameter :: PrmtopAtom   = 1
  integer, public, parameter :: PrmtopNbPrm  = 2
  integer, public, parameter :: PrmtopRes    = 3
  integer, public, parameter :: PrmtopBndUnq = 4
  integer, public, parameter :: PrmtopAngUnq = 5
  integer, public, parameter :: PrmtopDihUnq = 6
  integer, public, parameter :: PrmtopSolty  = 7
  integer, public, parameter :: PrmtopLennrd = 8
  integer, public, parameter :: PrmtopBonh   = 9
  integer, public, parameter :: PrmtopBona   = 10
  integer, public, parameter :: PrmtopAngh   = 11
  integer, public, parameter :: PrmtopAnga   = 12
  integer, public, parameter :: PrmtopDihh   = 13
  integer, public, parameter :: PrmtopDiha   = 14
  integer, public, parameter :: PrmtopHbond  = 15
  integer, public, parameter :: PrmtopExclud = 16
  integer, public, parameter :: PrmtopNSPM   = 17

  ! subroutines
  public  :: init_prmtop
  public  :: alloc_prmtop
  public  :: dealloc_prmtop
  public  :: dealloc_prmtop_all
  public  :: input_prmtop
  public  :: output_prmtop
  private :: read_prmtop
  private :: read_prmtop_main
  private :: read_prmtop_pointers
  private :: read_prmtop_char_format
  private :: read_prmtop_int_format
  private :: read_prmtop_real_format
  private :: write_prmtop

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_prmtop
  !> @brief        initialize PRMTOP information
  !! @authors      YSFJ
  !! @param[inout] prmtop : structure of PRMTOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_prmtop(prmtop)

    ! formal arguments
    type(s_prmtop),          intent(inout) :: prmtop


    prmtop%num_atoms          = 0 
    prmtop%num_residues       = 0
    prmtop%num_types          = 0
    prmtop%num_bondh          = 0
    prmtop%num_mbonda         = 0
    prmtop%num_anglh          = 0
    prmtop%num_mangla         = 0
    prmtop%num_diheh          = 0
    prmtop%num_mdihea         = 0
    prmtop%num_hparm          = 0
    prmtop%num_parm           = 0
    prmtop%num_nb             = 0
    prmtop%num_bona           = 0
    prmtop%num_theta          = 0
    prmtop%num_phia           = 0
    prmtop%num_uniqbond       = 0
    prmtop%num_uniqangl       = 0
    prmtop%num_uniqdihe       = 0
    prmtop%num_types_prm      = 0
    prmtop%num_bondtypes      = 0
    prmtop%ifpert             = 0
    prmtop%num_bond_perturbed = 0
    prmtop%num_angl_perturbed = 0
    prmtop%num_dihe_perturbed = 0
    prmtop%num_mbper          = 0
    prmtop%num_mgper          = 0
    prmtop%num_mdper          = 0
    prmtop%ifbox              = 0
    prmtop%num_max_res_atoms  = 0
    prmtop%ifcap              = 0
    prmtop%num_extra_points   = 0
    prmtop%num_copy           = 0
    prmtop%iptres             = 0
    prmtop%nspm               = 0
    prmtop%nspsol             = 0
    
    return

  end subroutine init_prmtop

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_prmtop
  !> @brief        allocate PRMTOP information
  !! @authors      YSFJ
  !! @param[inout] prmtop   : structure of PRMTOP information
  !! @param[in]    variable : selected variable
  !! @param[in]    var_size : size of the selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
 
  subroutine alloc_prmtop(prmtop, variable, var_size)

    ! formal arguments
    type(s_prmtop),          intent(inout) :: prmtop
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat


    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables

    select case (variable)

    case(PrmtopAtom)

      if (allocated(prmtop%atom_name)) then
        if (size(prmtop%atom_name) == var_size) return
        deallocate(prmtop%atom_name,                 &
                   prmtop%atom_cls_name,             &
                   prmtop%atom_cls_no,               &
                   prmtop%charge,                    &
                   prmtop%mass,                      &
                   prmtop%excl_atom,                 &
                   prmtop%classify,                  &
                   prmtop%join_array,                &
                   prmtop%irotate,                   &
                   prmtop%radi_born,                 &
                   prmtop%fs_born,                   &
                   stat = dealloc_stat)
      end if

      allocate(prmtop%atom_name(var_size),           &
               prmtop%atom_cls_name(var_size),       &
               prmtop%atom_cls_no(var_size),         &
               prmtop%charge(var_size),              &
               prmtop%mass(var_size),                &
               prmtop%excl_atom(var_size),           &
               prmtop%classify(var_size),            &
               prmtop%join_array(var_size),          &
               prmtop%irotate(var_size),             &
               prmtop%radi_born(var_size),           &
               prmtop%fs_born(var_size),             &
               stat = alloc_stat)

      prmtop%atom_name(:) = ''

    case(PrmtopNbPrm)

      if (allocated(prmtop%nb_par_idx)) then
        if (size(prmtop%nb_par_idx) == var_size) return
        deallocate(prmtop%nb_par_idx,                &
                   stat = dealloc_stat)
      end if

      allocate(prmtop%nb_par_idx(var_size),          &
               stat = alloc_stat)

    case(PrmtopRes)

      if (allocated(prmtop%res_label)) then
        if (size(prmtop%res_label) == var_size) return
        deallocate(prmtop%res_label,                 &
                   prmtop%res_point,                 &
                   stat = dealloc_stat)
      end if

      allocate(prmtop%res_label(var_size),           &
               prmtop%res_point(var_size),           &
               stat = alloc_stat)

      prmtop%res_label(:) = ''

    case(PrmtopBndUnq)

      if (allocated(prmtop%bond_fcons_uniq)) then
        if (size(prmtop%bond_fcons_uniq) == var_size) return
        deallocate(prmtop%bond_fcons_uniq,           &
                   prmtop%bond_equil_uniq,           &
                   stat = dealloc_stat)
      end if

      allocate(prmtop%bond_fcons_uniq(var_size),     &
               prmtop%bond_equil_uniq(var_size),     &
               stat = alloc_stat)

    case(PrmtopAngUnq)

      if (allocated(prmtop%angl_fcons_uniq)) then
        if (size(prmtop%angl_fcons_uniq) == var_size) return
        deallocate(prmtop%angl_fcons_uniq,           &
                   prmtop%angl_equil_uniq,           &
                   stat = dealloc_stat)
      end if

      allocate(prmtop%angl_fcons_uniq(var_size),     &
               prmtop%angl_equil_uniq(var_size),     &
               stat = alloc_stat)

    case(PrmtopDihUnq)

      if (allocated(prmtop%dihe_fcons_uniq)) then
        if (size(prmtop%dihe_fcons_uniq) == var_size) return
        deallocate(prmtop%dihe_fcons_uniq,           &
                   prmtop%dihe_perio_uniq,           &
                   prmtop%dihe_phase_uniq,           &
                   prmtop%scee_scale_fact,           &
                   prmtop%scnb_scale_fact,           &
                   stat = dealloc_stat)
      end if

      allocate(prmtop%dihe_fcons_uniq(var_size),     &
               prmtop%dihe_perio_uniq(var_size),     &
               prmtop%dihe_phase_uniq(var_size),     &
               prmtop%scee_scale_fact(var_size),     &
               prmtop%scnb_scale_fact(var_size),     &
               stat = alloc_stat)

    case(PrmtopSolty)

      if (allocated(prmtop%solty)) then
        if (size(prmtop%solty) == var_size) return
        deallocate(prmtop%solty,                     &
                   stat = dealloc_stat)
      end if

      allocate(prmtop%solty(var_size),               &
               stat = alloc_stat)

    case(PrmtopLennrd)

      if (allocated(prmtop%lennarda)) then
        if (size(prmtop%lennarda) == var_size) return
        deallocate(prmtop%lennarda,                  &
                   prmtop%lennardb,                  &
                   stat = dealloc_stat)
      end if

      allocate(prmtop%lennarda(var_size),            &
               prmtop%lennardb(var_size),            &
               stat = alloc_stat)

    case(PrmtopBonh)

      if (allocated(prmtop%bond_inc_hy)) then
        if (size(prmtop%bond_inc_hy(1,:)) == var_size) return
        deallocate(prmtop%bond_inc_hy,               &
                   stat = dealloc_stat)
      end if

      allocate(prmtop%bond_inc_hy(3,var_size),       &
               stat = alloc_stat)

    case(PrmtopBona)

      if (allocated(prmtop%bond_wo_hy)) then
        if (size(prmtop%bond_wo_hy(1,:)) == var_size) return
        deallocate(prmtop%bond_wo_hy,                &
                   stat = dealloc_stat)
      end if

      allocate(prmtop%bond_wo_hy(3,var_size),        &
               stat = alloc_stat)

    case(PrmtopAngh)

      if (allocated(prmtop%angl_inc_hy)) then
        if (size(prmtop%angl_inc_hy(1,:)) == var_size) return
        deallocate(prmtop%angl_inc_hy,               &
                   stat = dealloc_stat)
      end if

      allocate(prmtop%angl_inc_hy(4,var_size),       &
               stat = alloc_stat)

    case(PrmtopAnga)

      if (allocated(prmtop%angl_wo_hy)) then
        if (size(prmtop%angl_wo_hy(1,:)) == var_size) return
        deallocate(prmtop%angl_wo_hy,                &
                   stat = dealloc_stat)
      end if

      allocate(prmtop%angl_wo_hy(4,var_size),        &
               stat = alloc_stat)

    case(PrmtopDihh)

      if (allocated(prmtop%dihe_inc_hy)) then
        if (size(prmtop%dihe_inc_hy(1,:)) == var_size) return
        deallocate(prmtop%dihe_inc_hy,               &
                   stat = dealloc_stat)
      end if

      allocate(prmtop%dihe_inc_hy(5,var_size),       &
               stat = alloc_stat)

    case(PrmtopDiha)

      if (allocated(prmtop%dihe_wo_hy)) then
        if (size(prmtop%dihe_wo_hy(1,:)) == var_size) return
        deallocate(prmtop%dihe_wo_hy,                &
                   stat = dealloc_stat)
      end if

      allocate(prmtop%dihe_wo_hy(5,var_size),        &
               stat = alloc_stat)

    case (PrmtopHbond)

      if (allocated(prmtop%hbond_acoef)) then
        if (size(prmtop%hbond_acoef) == var_size) return
        deallocate(prmtop%hbond_acoef,               &
                   prmtop%hbond_bcoef,               &
                   prmtop%hb_cut,                    &
                   stat = dealloc_stat)
      end if

      allocate(prmtop%hbond_acoef(var_size),         &
               prmtop%hbond_bcoef(var_size),         &
               prmtop%hb_cut(var_size),              &
               stat = alloc_stat)

    case (PrmtopExclud)

      if (allocated(prmtop%inb)) then
        if (size(prmtop%inb) == var_size) return
        deallocate(prmtop%inb,                       &
                   stat = dealloc_stat)
      end if

      allocate(prmtop%inb(var_size),                 &
               stat = alloc_stat)

    case (PrmtopNSPM)

      if (allocated(prmtop%nsp)) then
        if (size(prmtop%nsp) == var_size) return
        deallocate(prmtop%nsp,                       &
                   stat = dealloc_stat)
      end if

      allocate(prmtop%nsp(var_size),                 &
               stat = alloc_stat)

    case default

      call error_msg('Alloc_Prmtop> bad variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_prmtop

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_prmtop
  !> @brief        deallocate PRMTOP information
  !! @authors      YSFJ
  !! @param[inout] prmtop   : structure of PRMTOP information
  !! @param[in]    variable : selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_prmtop(prmtop, variable)

    ! formal arguments
    type(s_prmtop),          intent(inout) :: prmtop
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat
   

    dealloc_stat = 0

    ! deallocate selected variable
    !
    select case (variable)

    case(PrmtopAtom)

      if (allocated(prmtop%atom_name)) then
        deallocate(prmtop%atom_name,             &
                   prmtop%atom_cls_name,         &
                   prmtop%atom_cls_no,           &
                   prmtop%charge,                &
                   prmtop%mass,                  &
                   prmtop%excl_atom,             &
                   prmtop%classify,              &
                   prmtop%join_array,            &
                   prmtop%irotate,               &
                   prmtop%radi_born,             &
                   prmtop%fs_born,               &
                   stat = dealloc_stat)
      end if

    case (PrmtopNbPrm)

      if (allocated(prmtop%nb_par_idx)) then
        deallocate(prmtop%nb_par_idx,            &
                   stat = dealloc_stat)
      end if

    case (PrmtopRes)

      if (allocated(prmtop%res_label)) then
        deallocate(prmtop%res_label,             &
                   prmtop%res_point,             &
                   stat = dealloc_stat)
       end if

    case (PrmtopBndUnq)

      if (allocated(prmtop%bond_fcons_uniq)) then
        deallocate(prmtop%bond_fcons_uniq,       &
                   prmtop%bond_equil_uniq,       &
                   stat = dealloc_stat)
      end if

    case (PrmtopAngUnq)

      if (allocated(prmtop%angl_fcons_uniq)) then
        deallocate(prmtop%angl_fcons_uniq,       &
                   prmtop%angl_equil_uniq,       &
                   stat = dealloc_stat)
      end if

    case (PrmtopDihUnq)

      if (allocated(prmtop%dihe_fcons_uniq)) then
        deallocate(prmtop%dihe_fcons_uniq,       &
                   prmtop%dihe_perio_uniq,       &
                   prmtop%dihe_phase_uniq,       &
                   prmtop%scee_scale_fact,       &
                   prmtop%scnb_scale_fact,       &
                   stat = dealloc_stat)
      end if

    case (PrmtopSolty)

      if (allocated(prmtop%solty)) then
        deallocate(prmtop%solty,                 &
                   stat = dealloc_stat)
      end if

    case(PrmtopLennrd)

      if (allocated(prmtop%lennarda)) then
        deallocate(prmtop%lennarda,              &
                   prmtop%lennardb,              &
                   stat = dealloc_stat)
      end if

    case(PrmtopBonh)

      if (allocated(prmtop%bond_inc_hy)) then
        deallocate(prmtop%bond_inc_hy,           &
                   stat = dealloc_stat)
      end if

    case(PrmtopBona)

      if (allocated(prmtop%bond_wo_hy)) then
        deallocate(prmtop%bond_wo_hy,            &
                   stat = dealloc_stat)
      end if

    case(PrmtopAngh)

      if (allocated(prmtop%angl_inc_hy)) then
        deallocate(prmtop%angl_inc_hy,           &
                   stat = dealloc_stat)
      end if

    case(PrmtopAnga)

      if (allocated(prmtop%angl_wo_hy)) then
        deallocate(prmtop%angl_wo_hy,            &
                   stat = dealloc_stat)
      end if

    case(PrmtopDihh)

      if (allocated(prmtop%dihe_inc_hy)) then
        deallocate(prmtop%dihe_inc_hy,           &
                   stat = dealloc_stat)
      end if

    case(PrmtopDiha)

      if (allocated(prmtop%dihe_wo_hy)) then
        deallocate(prmtop%dihe_wo_hy,            &
                   stat = dealloc_stat)
      end if

    case (PrmtopHbond)

      if (allocated(prmtop%hbond_acoef)) then
        deallocate(prmtop%hbond_acoef,           &
                   prmtop%hbond_bcoef,           &
                   prmtop%hb_cut,                &
                   stat = dealloc_stat)
      end if

    case (PrmtopExclud)

      if (allocated(prmtop%inb)) then
        deallocate(prmtop%inb,                   &
                   stat = dealloc_stat)
      end if

    case (PrmtopNSPM)

      if (allocated(prmtop%nsp)) then
        deallocate(prmtop%nsp,                   &
                   stat = dealloc_stat)
      end if

    case default

      call error_msg('Dealloc_Prmtop> bad variable')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_prmtop

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_prmtop_all
  !> @brief        deallocate all PRMTOP information
  !! @authors      YSFJ
  !! @param[inout] prmtop : structure of PRMTOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_prmtop_all(prmtop)

    ! formal arguments
    type(s_prmtop),          intent(inout) :: prmtop


    call dealloc_prmtop(prmtop, PrmtopAtom  )
    call dealloc_prmtop(prmtop, PrmtopNbPrm )
    call dealloc_prmtop(prmtop, PrmtopRes   )
    call dealloc_prmtop(prmtop, PrmtopBndUnq)
    call dealloc_prmtop(prmtop, PrmtopAngUnq)
    call dealloc_prmtop(prmtop, PrmtopDihUnq)
    call dealloc_prmtop(prmtop, PrmtopSolty )
    call dealloc_prmtop(prmtop, PrmtopLennrd)
    call dealloc_prmtop(prmtop, PrmtopBonh  )
    call dealloc_prmtop(prmtop, PrmtopBona  )
    call dealloc_prmtop(prmtop, PrmtopAngh  )
    call dealloc_prmtop(prmtop, PrmtopAnga  )
    call dealloc_prmtop(prmtop, PrmtopDihh  )
    call dealloc_prmtop(prmtop, PrmtopDiha  )
    call dealloc_prmtop(prmtop, PrmtopExclud)
    call dealloc_prmtop(prmtop, PrmtopHbond )
    call dealloc_prmtop(prmtop, PrmtopNSPM  )

    return

  end subroutine dealloc_prmtop_all

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_prmtop
  !> @brief        a driver subroutine for reading PRMTOP file
  !! @authors      YSFJ
  !! @param[in]    prmtop_filename : filename of PRMTOP file
  !! @param[inout] prmtop          : structure of PRMTOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_prmtop(prmtop_filename, prmtop)

    ! formal arguments
    character(*),            intent(in)    :: prmtop_filename
    type(s_prmtop),          intent(inout) :: prmtop
    
    ! local variables
    integer                  :: file


    ! open PRMTOP file
    !
    call open_file(file, prmtop_filename, IOFileInput)

    ! read coordinate data from PRMTOP file
    !
    call read_prmtop(file, prmtop)

    ! close PRMTOP file
    !
    call close_file(file)


    return

  end subroutine input_prmtop

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_prmtop
  !> @brief        a driver subroutine for writing PRMTOP file
  !! @authors      YSFJ
  !! @param[in]    psf_filename : filename of PRMTOP file
  !! @param[in]    psf          : structure of PRMTOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_prmtop(prmtop_filename, prmtop)

    ! formal arguments
    character(*),            intent(in)    :: prmtop_filename
    type(s_prmtop),          intent(in)    :: prmtop

    ! local variables
    integer                  :: file


    ! open PRMTOP file
    !
    call open_file(file, prmtop_filename, IOFileOutputNew)

    ! write coordinate data from PRMTOP file
    !
    call write_prmtop(file, prmtop)

    ! close PRMTOP file
    !
    call close_file(file)

    return

  end subroutine output_prmtop

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_prmtop
  !> @brief        read data from PRMTOP file
  !! @authors      YSFJ
  !! @param[in]    file    : unit number of PRMTOP file
  !! @param[inout] prmtop  : prmtop data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_prmtop(file, prmtop)

    ! formal arugments
    integer,                 intent(in)   :: file
    type(s_prmtop),          intent(inout):: prmtop
 
    ! local variables
    integer                  :: i, j


    ! deallocate old data
    !
    call init_prmtop(prmtop)

    ! read prmtop main section
    !
    call read_prmtop_main(file, prmtop)


    ! write summary of PRMTOP information
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Prmtop> Summary of PRMTOP file'

      write(MsgOut,'(A20,I10,A20,I10)')                       &
           '  NATOM           = ', prmtop%num_atoms,          &
           '  NTYPES          = ', prmtop%num_types
      write(MsgOut,'(A20,I10,A20,I10)')                       &
           '  NBONH           = ', prmtop%num_bondh,          &
           '  MBONA           = ', prmtop%num_mbonda
      write(MsgOut,'(A20,I10,A20,I10)')                       &
           '  NTHETH          = ', prmtop%num_anglh,          &
           '  MTHETA          = ', prmtop%num_mangla
      write(MsgOut,'(A20,I10,A20,I10)')                       &
           '  NPHIH           = ', prmtop%num_diheh,          &
           '  MPHIA           = ', prmtop%num_mdihea
      write(MsgOut,'(A20,I10,A20,I10)')                       &
           '  NPARM           = ', prmtop%num_parm,           &
           '  NNB             = ', prmtop%num_nb
      write(MsgOut,'(A20,I10,A20,I10)')                       &
           '  NRES            = ', prmtop%num_residues,       &
           '  NBONA           = ', prmtop%num_bona
      write(MsgOut,'(A20,I10,A20,I10)')                       &
           '  NTHETA          = ', prmtop%num_theta,          &
           '  NPHIA           = ', prmtop%num_phia
      write(MsgOut,'(A20,I10,A20,I10)')                       &
           '  NUMBND          = ', prmtop%num_uniqbond,       &
           '  NUMANG          = ', prmtop%num_uniqangl
      write(MsgOut,'(A20,I10,A20,I10)')                       &
           '  NPTRA           = ', prmtop%num_uniqdihe,       &
           '  NATYP           = ', prmtop%num_types_prm
      write(MsgOut,'(A20,I10)')                               &
           '  NPHB            = ', prmtop%num_bondtypes
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine read_prmtop

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_prmtop_main
  !> @brief        read data from PRMTOP file
  !! @authors      YSFJ
  !! @param[in]    file    : unit number of PRMTOP file
  !! @param[inout] prmtop  : prmtop data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_prmtop_main(file, prmtop)

    ! formal arugments
    integer,                 intent(in)    :: file
    type(s_prmtop),          intent(inout) :: prmtop

    ! local variables
    integer                  :: i
    integer                  :: size
    character(80)            :: line
    character(80)            :: key, temp
    character(80)            :: fmt01, fmt02, fmt03
 

    do while(.true.)

      read(file, '(A80)',end=100,err=100) line

      if      (line(1:8) == '%VERSION') then
        prmtop%lversion = .true.

        prmtop%version = line(1:80)

      else if (line(1:5) == '%FLAG') then

        backspace(file)
        read(file,'(A6,A)',end=100,err=100) temp, key

        read(file,'(A80)') line
        if (line(1:7) == '%FORMAT') then

          read(line(9:),'(A)') fmt03
          do i = 1, len_trim(fmt03)
            if (fmt03(i:i) == ')') then
              fmt03(i:i) = ' '
            end if
          end do

          fmt03 = trim(fmt03)

          do i = 1, len_trim(fmt03)
            if (fmt03(i:i) == 'A' .or. fmt03(i:i) == 'a' .or. &
                fmt03(i:i) == 'I' .or. fmt03(i:i) == 'i' .or. &
                fmt03(i:i) == 'E' .or. fmt03(i:i) == 'e') then

              read(fmt03(1:i-1),'(A)') fmt01
              read(fmt03(i+1:),'(A)')  fmt02

              exit
            end if
          end do
          fmt02 = trim(fmt02)

        else
          call error_msg( &
               'Read_Prmtop_Main> error! Please set %FORMAT line for all term.')
        end if

        if      (key(1: 5) == 'TITLE') then
          prmtop%ltitle = .true.

          read(file,'(A80)') line
          if (line(1:7) == '%FORMAT') then
            read(file,'(A80)') line
          end if
          prmtop%title = trim(line)

        else if (key(1: 8) == 'POINTERS') then
          prmtop%lpointers = .true.
          call read_prmtop_pointers(file, fmt01, fmt02, prmtop)

          ! allocate variables
          !
          call alloc_prmtop(prmtop, PrmtopAtom,   prmtop%num_atoms)
          size = prmtop%num_types * prmtop%num_types
          call alloc_prmtop(prmtop, PrmtopNbPrm, size)
          call alloc_prmtop(prmtop, PrmtopRes,    prmtop%num_residues)
          call alloc_prmtop(prmtop, PrmtopBndUnq, prmtop%num_uniqbond)
          call alloc_prmtop(prmtop, PrmtopAngUnq, prmtop%num_uniqangl)
          call alloc_prmtop(prmtop, PrmtopDihUnq, prmtop%num_uniqdihe)
          call alloc_prmtop(prmtop, PrmtopSolty,  prmtop%num_types_prm)
          size = prmtop%num_types * (prmtop%num_types + 1) * 0.5_wp
          call alloc_prmtop(prmtop, PrmtopLennrd, size)
          call alloc_prmtop(prmtop, PrmtopBonh,   prmtop%num_bondh)
          call alloc_prmtop(prmtop, PrmtopBona,   prmtop%num_mbonda)
          call alloc_prmtop(prmtop, PrmtopAngh,   prmtop%num_anglh)
          call alloc_prmtop(prmtop, PrmtopAnga,   prmtop%num_mangla)
          call alloc_prmtop(prmtop, PrmtopDihh,   prmtop%num_diheh)
          call alloc_prmtop(prmtop, PrmtopDiha,   prmtop%num_mdihea)
          call alloc_prmtop(prmtop, PrmtopExclud, prmtop%num_nb)
          call alloc_prmtop(prmtop, PrmtopHbond,  prmtop%num_bondtypes)


        else if (key(1: 9) == 'ATOM_NAME') then
          prmtop%latom_name = .true.
          call read_prmtop_char_format(file, fmt01, fmt02,          &
                                       'ATOM_NAME',                 &
                                       prmtop%num_atoms, prmtop) 

        else if (key(1: 6) == 'CHARGE') then
          prmtop%lcharge = .true.
          call read_prmtop_real_format(file, fmt01, fmt02,          &
                                       'CHARGE',                    &
                                       prmtop%num_atoms, prmtop) 

        else if (key(1: 4) == 'MASS') then
          prmtop%lmass = .true.
          call read_prmtop_real_format(file, fmt01, fmt02,          &
                                       'MASS',                      &
                                       prmtop%num_atoms, prmtop) 

        else if (key(1:15) == 'ATOM_TYPE_INDEX') then
          prmtop%latom_type_index = .true.
          call read_prmtop_int_format(file, fmt01, fmt02,           &
                                       'ATOM_TYPE_INDEX',           &
                                       prmtop%num_atoms, prmtop) 

        else if (key(1:21) == 'NUMBER_EXCLUDED_ATOMS') then
          prmtop%lnumber_excluded_atoms = .true.
          call read_prmtop_int_format(file, fmt01, fmt02,           &
                                       'NUMBER_EXCLUDED_ATOMS',     &
                                       prmtop%num_atoms, prmtop) 

        else if (key(1:20) == 'NONBONDED_PARM_INDEX') then
          prmtop%lnonbonded_parm_index = .true.
          size = prmtop%num_types * prmtop%num_types
          call read_prmtop_int_format(file, fmt01, fmt02,           &
                                       'NONBONDED_PARM_INDEX',      &
                                       size, prmtop) 

        else if (key(1:13) == 'RESIDUE_LABEL') then
          prmtop%lresidue_label = .true.
          call read_prmtop_char_format(file, fmt01, fmt02,          &
                                       'RESIDUE_LABEL',             &
                                       prmtop%num_residues, prmtop) 

        else if (key(1:15) == 'RESIDUE_POINTER') then
          prmtop%lresidue_pointer = .true.
          call read_prmtop_int_format(file, fmt01, fmt02,           &
                                       'RESIDUE_POINTER',           &
                                       prmtop%num_residues, prmtop) 

        else if (key(1:19) == 'BOND_FORCE_CONSTANT') then
          prmtop%lbond_force_constant = .true.
          call read_prmtop_real_format(file, fmt01, fmt02,          &
                                       'BOND_FORCE_CONSTANT',       &
                                       prmtop%num_uniqbond, prmtop) 

        else if (key(1:16) == 'BOND_EQUIL_VALUE') then
          prmtop%lbond_equil_value = .true.
          call read_prmtop_real_format(file, fmt01, fmt02,          &
                                       'BOND_EQUIL_VALUE',          &
                                       prmtop%num_uniqbond, prmtop) 

        else if (key(1:20) == 'ANGLE_FORCE_CONSTANT') then
          prmtop%langle_force_constant = .true.
          call read_prmtop_real_format(file, fmt01, fmt02,          &
                                       'ANGLE_FORCE_CONSTANT',      &
                                       prmtop%num_uniqangl, prmtop) 

        else if (key(1:17) == 'ANGLE_EQUIL_VALUE') then
          prmtop%langle_equil_value = .true.
          call read_prmtop_real_format(file, fmt01, fmt02,          &
                                       'ANGLE_EQUIL_VALUE',         &
                                       prmtop%num_uniqangl, prmtop) 

        else if (key(1:23) == 'DIHEDRAL_FORCE_CONSTANT') then
          prmtop%ldihedral_force_constant = .true.
          call read_prmtop_real_format(file, fmt01, fmt02,          &
                                       'DIHEDRAL_FORCE_CONSTANT',   &
                                       prmtop%num_uniqdihe, prmtop) 

        else if (key(1:20) == 'DIHEDRAL_PERIODICITY') then
          prmtop%ldihedral_periodicity = .true.
          call read_prmtop_real_format(file, fmt01, fmt02,          &
                                       'DIHEDRAL_PERIODICITY',      &
                                       prmtop%num_uniqdihe, prmtop) 

        else if (key(1:14) == 'DIHEDRAL_PHASE') then
          prmtop%ldihedral_phase = .true.
          call read_prmtop_real_format(file, fmt01, fmt02,          &
                                       'DIHEDRAL_PHASE',            &
                                       prmtop%num_uniqdihe, prmtop) 

        else if (key(1: 5) == 'SOLTY') then  ! currently unused
          prmtop%lsolty = .true.
          call read_prmtop_real_format(file, fmt01, fmt02,          &
                                       'SOLTY',                     &
                                       prmtop%num_types_prm, prmtop) 

        else if (key(1:17) == 'SCEE_SCALE_FACTOR') then
          prmtop%lscee_scale_factor = .true.
          call read_prmtop_real_format(file, fmt01, fmt02,          &
                                       'SCEE_SCALE_FACTOR',         &
                                       prmtop%num_uniqdihe, prmtop) 

        else if (key(1:17) == 'SCNB_SCALE_FACTOR') then
          prmtop%lscnb_scale_factor = .true.
          call read_prmtop_real_format(file, fmt01, fmt02,          &
                                       'SCNB_SCALE_FACTOR',         &
                                       prmtop%num_uniqdihe, prmtop) 

        else if (key(1:19) == 'LENNARD_JONES_ACOEF') then
          prmtop%llennard_jones_acoef = .true.
          size = prmtop%num_types * (prmtop%num_types + 1) * 0.5_wp
          call read_prmtop_real_format(file, fmt01, fmt02,          &
                                       'LENNARD_JONES_ACOEF',       &
                                       size, prmtop) 

        else if (key(1:19) == 'LENNARD_JONES_BCOEF') then
          prmtop%llennard_jones_bcoef = .true.
          size = prmtop%num_types * (prmtop%num_types + 1) * 0.5_wp
          call read_prmtop_real_format(file, fmt01, fmt02,          &
                                       'LENNARD_JONES_BCOEF',       &
                                       size, prmtop) 

        else if (key(1:18) == 'BONDS_INC_HYDROGEN') then
          prmtop%lbonds_inc_hydrogen = .true.
          call read_prmtop_int_format(file, fmt01, fmt02,           &
                                       'BONDS_INC_HYDROGEN',        &
                                       prmtop%num_bondh, prmtop) 

        else if (key(1:22) == 'BONDS_WITHOUT_HYDROGEN') then
          prmtop%lbonds_without_hydrogen = .true.
          call read_prmtop_int_format(file, fmt01, fmt02,           &
                                       'BONDS_WITHOUT_HYDROGEN',    &
                                       prmtop%num_mbonda, prmtop) 

        else if (key(1:19) == 'ANGLES_INC_HYDROGEN') then
          prmtop%langles_inc_hydrogen = .true.
          call read_prmtop_int_format(file, fmt01, fmt02,           &
                                       'ANGLES_INC_HYDROGEN',       &
                                       prmtop%num_anglh, prmtop) 

        else if (key(1:23) == 'ANGLES_WITHOUT_HYDROGEN') then
          prmtop%langles_without_hydrogen = .true.
          call read_prmtop_int_format(file, fmt01, fmt02,           &
                                       'ANGLES_WITHOUT_HYDROGEN',   &
                                       prmtop%num_mangla, prmtop) 

        else if (key(1:22) == 'DIHEDRALS_INC_HYDROGEN') then
          prmtop%ldihedrals_inc_hydrogen = .true.
          call read_prmtop_int_format(file, fmt01, fmt02,           &
                                       'DIHEDRALS_INC_HYDROGEN',    &
                                       prmtop%num_diheh, prmtop) 

        else if (key(1:26) == 'DIHEDRALS_WITHOUT_HYDROGEN') then
          prmtop%ldihedrals_without_hydrogen = .true.
          call read_prmtop_int_format(file, fmt01, fmt02,           &
                                       'DIHEDRALS_WITHOUT_HYDROGEN',&
                                       prmtop%num_mdihea, prmtop) 

        else if (key(1:19) == 'EXCLUDED_ATOMS_LIST') then
          prmtop%lexcluded_atoms_list = .true.
          call read_prmtop_int_format(file, fmt01, fmt02,           &
                                       'EXCLUDED_ATOMS_LIST',       &
                                       prmtop%num_nb, prmtop) 

        else if (key(1:11) == 'HBOND_ACOEF') then
          prmtop%lhbond_acoef = .true.
          call read_prmtop_real_format(file, fmt01, fmt02,          &
                                       'HBOND_ACOEF',               &
                                       prmtop%num_bondtypes, prmtop) 

        else if (key(1:11) == 'HBOND_BCOEF') then
          prmtop%lhbond_bcoef = .true.
          call read_prmtop_real_format(file, fmt01, fmt02,          &
                                       'HBOND_BCOEF',               &
                                       prmtop%num_bondtypes, prmtop) 

        else if (key(1: 5) == 'HBCUT') then ! no longer in use
          prmtop%lhbcut = .true.
          call read_prmtop_real_format(file, fmt01, fmt02,          &
                                       'HBCUT',                     &
                                       prmtop%num_bondtypes, prmtop) 

        else if (key(1:15) == 'AMBER_ATOM_TYPE') then
          prmtop%lamber_atom_type = .true.
          call read_prmtop_char_format(file, fmt01, fmt02,          &
                                       'AMBER_ATOM_TYPE',           &
                                       prmtop%num_atoms, prmtop) 

        else if (key(1:25) == 'TREE_CHAIN_CLASSIFICATION') then
          prmtop%ltree_chain_classification = .true.
          call read_prmtop_char_format(file, fmt01, fmt02,          &
                                       'TREE_CHAIN_CLASSIFICATION', &
                                       prmtop%num_atoms, prmtop) 

        else if (key(1:10) == 'JOIN_ARRAY') then ! Currently unused in sander
          prmtop%ljoin_array = .true.
          call read_prmtop_int_format(file, fmt01, fmt02,           &
                                       'JOIN_ARRAY',                &
                                       prmtop%num_atoms, prmtop) 

        else if (key(1: 6) == 'IROTAT') then
          prmtop%lirotat = .true.
          call read_prmtop_int_format(file, fmt01, fmt02,           &
                                       'IROTAT',                    &
                                       prmtop%num_atoms, prmtop) 

        else if (key(1:10) == 'RADIUS_SET') then
          prmtop%lradius_set = .true.
          call read_prmtop_char_format(file, fmt01, fmt02,          &
                                       'RADIUS_SET',                &
                                       1, prmtop) 

        else if (key(1: 5) == 'RADII') then
          prmtop%lradii = .true.
          call read_prmtop_real_format(file, fmt01, fmt02,          &
                                       'RADII',                     &
                                       prmtop%num_atoms, prmtop) 

        else if (key(1: 6) == 'SCREEN') then
          prmtop%lscreen = .true.
          call read_prmtop_real_format(file, fmt01, fmt02,          &
                                       'SCREEN',                    &
                                       prmtop%num_atoms, prmtop) 

          ! The following are only present if IFBOX .gt. 0
        else if (key(1:16) == 'SOLVENT_POINTERS') then
          prmtop%lsolvent_pointers = .true.
          call read_prmtop_int_format(file, fmt01, fmt02,           &
                                       'SOLVENT_POINTERS',          &
                                       3, prmtop) 
          call alloc_prmtop(prmtop, PrmtopNSPM,  prmtop%nspm)

        else if (key(1:18) == 'ATOMS_PER_MOLECULE') then
          prmtop%latoms_per_molecule = .true.
          if (prmtop%lsolvent_pointers) then
            call read_prmtop_int_format(file, fmt01, fmt02,         &
                                         'ATOMS_PER_MOLECULE',      &
                                         prmtop%nspm, prmtop) 
          end if

        else if (key(1:14) == 'BOX_DIMENSIONS') then
          prmtop%lbox_dimensions = .true.
          call read_prmtop_real_format(file, fmt01, fmt02,          &
                                       'BOX_DIMENSIONS',            &
               4, prmtop)

          ! The following are only present if IFCAP .gt. 0
        else if (key(1:8) == 'CAP_INFO') then
          prmtop%lcap_info = .true.

        else if (key(1:9) == 'CAP_INFO2') then
          prmtop%lcap_info2 = .true.

          ! The following are only present if IPOL .eq. 1
        else if (key(1:14) == 'POLARIZABILITY') then
          prmtop%lpolarizability = .true.

        end if
      end if
    end do

100 continue

    return
  end subroutine read_prmtop_main

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_prmtop_pointers
  !> @brief        read POINTERS information from PRMTOP file
  !! @authors      YSFJ
  !! @param[in]    file     : unit number of PRMTOP file
  !! @param[in]    fmt01    : number variables per line, 10 for10I8 
  !! @param[in]    fmt02    : column length per one variable, 8 for 10I8
  !! @param[inout] prmtop   : prmtop data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_prmtop_pointers(file, fmt01, fmt02, prmtop)

    ! parameters
    integer,                parameter     :: NUM_VARS = 40

    ! formal arguments
    integer,                intent(in)    :: file
    character(*),           intent(in)    :: fmt01
    character(*),           intent(in)    :: fmt02
    type(s_prmtop),         intent(inout) :: prmtop

    ! local variables
    integer                               :: i, j
    integer                               :: iwork(NUM_VARS) =0
    integer                               :: len
    integer                               :: count
    integer                               :: ifmt01, ifmt02
    integer                               :: start, end
    character(100)                        :: line
    character(10)                         :: fmt10


    read(fmt01,*) ifmt01
    read(fmt02,*) ifmt02

    len = ifmt01 * ifmt02

    if (ifmt02 < 10) then
      fmt10 = '(Ix)'
      write(fmt10(3:3),'(A1)') fmt02
    else if (ifmt02 < 100) then
      fmt10 = '(Ixx)'
      write(fmt10(3:4),'(A2)') fmt02
    end if

    count = 1
    do j = 1, 4
      read(file,'(A)') line
      len = len_trim(line)
      do i = 1, ifmt01
        start = (i - 1) * ifmt02 + 1
        end = start + ifmt02 - 1
        read(line(start:end),fmt10,end=100,err=100) iwork(count)
        count = count + 1
      end do
    end do
100 continue

    prmtop%num_atoms          = iwork(1)
    prmtop%num_types          = iwork(2)
    prmtop%num_bondh          = iwork(3)
    prmtop%num_mbonda         = iwork(4)
    prmtop%num_anglh          = iwork(5)
    prmtop%num_mangla         = iwork(6)
    prmtop%num_diheh          = iwork(7)
    prmtop%num_mdihea         = iwork(8)
    prmtop%num_hparm          = iwork(9)
    prmtop%num_parm           = iwork(10)
    prmtop%num_nb             = iwork(11)
    prmtop%num_residues       = iwork(12)
    prmtop%num_bona           = iwork(13)
    prmtop%num_theta          = iwork(14)
    prmtop%num_phia           = iwork(15)
    prmtop%num_uniqbond       = iwork(16)
    prmtop%num_uniqangl       = iwork(17)
    prmtop%num_uniqdihe       = iwork(18)
    prmtop%num_types_prm      = iwork(19)
    prmtop%num_bondtypes      = iwork(20)
    prmtop%ifpert             = iwork(21)
    prmtop%num_bond_perturbed = iwork(22)
    prmtop%num_angl_perturbed = iwork(23)
    prmtop%num_dihe_perturbed = iwork(24)
    prmtop%num_mbper          = iwork(25)
    prmtop%num_mgper          = iwork(26)
    prmtop%num_mdper          = iwork(27)
    prmtop%ifbox              = iwork(28)
    prmtop%num_max_res_atoms  = iwork(29)
    prmtop%ifcap              = iwork(30)
    prmtop%num_extra_points   = iwork(31)
    prmtop%num_copy           = iwork(32)
    
    return

  end subroutine read_prmtop_pointers

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_prmtop_char_format
  !> @brief        read character format style information from PRMTOP file
  !! @authors      YSFJ
  !! @param[in]    file    : unit number of PRMTOP file
  !! @param[in]    fmt01   : number variables per line, 10 for10I8 
  !! @param[in]    fmt02   : column length per one variable, 8 for 10I8
  !! @param[in]    keyword : input keyword information
  !! @param[in]    num     : number of variables
  !! @param[inout] prmtop  : prmtop data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_prmtop_char_format(file, fmt01, fmt02, keyword, num, prmtop) 

    ! formal arguments
    integer,                intent(in)    :: file
    character(*),           intent(in)    :: fmt01
    character(*),           intent(in)    :: fmt02
    character(*),           intent(in)    :: keyword
    integer,                intent(in)    :: num
    type(s_prmtop),         intent(inout) :: prmtop
    
    ! local variables
    integer                 :: i, j
    integer                 :: len
    integer                 :: count
    integer                 :: start, end
    integer                 :: ifmt01, ifmt02
    integer                 :: num_reads
    character(100)          :: line
    character(10)           :: fmt10


    read(fmt01,*) ifmt01
    read(fmt02,*) ifmt02
    if (ifmt02 < 10) then
      fmt10 = '(Ax)'
      write(fmt10(3:3),'(A1)') fmt02
    else if (ifmt02 < 100) then
      fmt10 = '(Axx)'
      write(fmt10(3:4),'(A2)') fmt02
    end if

    ! Read atom
    !
    len = ifmt01 * ifmt02

    num_reads = num / ifmt01
    if (mod(num,ifmt01) /= 0) then
      num_reads = num_reads + 1
    end if

    count = 1
    do j = 1, num_reads
      read(file,'(A)',err=100,end=100) line

      do i = 1, ifmt01
        if (count > num) then
          exit
        end if
        start = (i - 1) * ifmt02 + 1
        end = start + ifmt02 - 1

        if (keyword == 'ATOM_NAME') then
          read(line(start:end),fmt10,end=200,err=200) &
               prmtop%atom_name(count)
        else if (keyword == 'RESIDUE_LABEL') then
          read(line(start:end),fmt10,end=200,err=200) &
               prmtop%res_label(count)
        else if (keyword == 'AMBER_ATOM_TYPE') then
          read(line(start:end),fmt10,end=200,err=200) &
               prmtop%atom_cls_name(count)
        else if (keyword == 'TREE_CHAIN_CLASSIFICATION') then
          read(line(start:end),fmt10,end=200,err=200) &
               prmtop%classify(count)
        else if (keyword == 'RADIUS_SET') then
          read(line(start:end),fmt10,end=200,err=200) &
               prmtop%radius_set
        end if
        count = count + 1
      end do
200   continue

    end do
100 continue

    return 

  end subroutine read_prmtop_char_format

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_prmtop_int_format
  !> @brief        read int fomat style information from PRMTOP file
  !! @authors      YSFJ
  !! @param[in]    file    : unit number of PRMTOP file
  !! @param[in]    fmt01   : number variables per line, 10 for10I8 
  !! @param[in]    fmt02   : column length per one variable, 8 for 10I8
  !! @param[in]    keyword : input keyword information
  !! @param[in]    num     : number of variables
  !! @param[inout] prmtop  : prmtop data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_prmtop_int_format(file, fmt01, fmt02, keyword, num, prmtop) 

    ! formal arguments
    integer,                intent(in)    :: file
    character(*),           intent(in)    :: fmt01
    character(*),           intent(in)    :: fmt02
    character(*),           intent(in)    :: keyword
    integer,                intent(in)    :: num
    type(s_prmtop),         intent(inout) :: prmtop

    ! local variables
    integer                 :: i, j
    integer                 :: len
    integer                 :: count
    integer                 :: start, end
    integer                 :: ifmt01, ifmt02
    integer                 :: num_reads
    integer                 :: sub_count
    integer                 :: iwork
    integer                 :: iwork_num
    character(100)          :: line
    character(10)           :: fmt10


    read(fmt01,*) ifmt01
    read(fmt02,*) ifmt02
    len = len_trim(fmt02)

    if (ifmt02 < 10) then
      fmt10 = '(Ix)'
      write(fmt10(3:3),'(A1)') fmt02
    else if (ifmt02 < 100) then
      fmt10 = '(Ixx)'
      write(fmt10(3:4),'(A2)') fmt02
    end if

    iwork_num = num
    if      (keyword == 'BONDS_INC_HYDROGEN'       .or.       &
             keyword == 'BONDS_WITHOUT_HYDROGEN')     then
      iwork_num = num * 3
    else if (keyword == 'ANGLES_INC_HYDROGEN'      .or.       &
             keyword == 'ANGLES_WITHOUT_HYDROGEN')    then
      iwork_num = num * 4
    else if (keyword == 'DIHEDRALS_INC_HYDROGEN'   .or.       &
             keyword == 'DIHEDRALS_WITHOUT_HYDROGEN') then
      iwork_num = num * 5
    end if

    ! Read atom
    !
    len = ifmt01 * ifmt02

    num_reads = iwork_num / ifmt01
    if (mod(iwork_num,ifmt01) /= 0) then
      num_reads = num_reads + 1
    end if

    count = 1
    sub_count = 0

    do j = 1, num_reads
      read(file,'(A)',err=100,end=100) line

      do i = 1, ifmt01
        if (count > num) then
          exit
        end if
        start = (i - 1) * ifmt02 + 1
        end   = start + ifmt02 - 1

        if (keyword == 'ATOM_TYPE_INDEX') then
          read(line(start:end),fmt10) prmtop%atom_cls_no(count)

          count = count + 1
        else if (keyword == 'NUMBER_EXCLUDED_ATOMS') then
          read(line(start:end),fmt10) prmtop%excl_atom(count)

          count = count + 1
        else if (keyword == 'NONBONDED_PARM_INDEX') then
          read(line(start:end),fmt10) prmtop%nb_par_idx(count)

          count = count + 1

        else if (keyword == 'RESIDUE_POINTER') then
          read(line(start:end),fmt10) prmtop%res_point(count)

          count = count + 1

        else if (keyword == 'JOIN_ARRAY') then
          read(line(start:end),fmt10) prmtop%join_array(count)

          count = count + 1

        else if (keyword == 'IROTAT') then
          read(line(start:end),fmt10) prmtop%irotate(count)

          count = count + 1

        else if (keyword == 'EXCLUDED_ATOMS_LIST') then
          read(line(start:end),fmt10) prmtop%inb(count)

          count = count + 1

        else if (keyword == 'ATOMS_PER_MOLECULE') then
          read(line(start:end),fmt10) prmtop%nsp(count)

          count = count + 1

        else if (keyword == 'BONDS_INC_HYDROGEN') then
          read(line(start:end),fmt10) iwork
          sub_count = sub_count + 1
          prmtop%bond_inc_hy(sub_count,count) = iwork

          if (sub_count == 3) then
            count = count + 1
            sub_count = 0
          end if

        else if (keyword == 'BONDS_WITHOUT_HYDROGEN') then
          read(line(start:end),fmt10) iwork
          sub_count = sub_count + 1
          prmtop%bond_wo_hy(sub_count,count) = iwork

          if (sub_count == 3) then
            count = count + 1
            sub_count = 0
          end if

        else if (keyword == 'ANGLES_INC_HYDROGEN') then
          read(line(start:end),fmt10) iwork
          sub_count = sub_count + 1
          prmtop%angl_inc_hy(sub_count,count) = iwork

          if (sub_count == 4) then
            count = count + 1
            sub_count = 0
          end if

        else if (keyword == 'ANGLES_WITHOUT_HYDROGEN') then
          read(line(start:end),fmt10) iwork
          sub_count = sub_count + 1
          prmtop%angl_wo_hy(sub_count,count) = iwork

          if (sub_count == 4) then
            count = count + 1
            sub_count = 0
          end if

        else if (keyword == 'DIHEDRALS_INC_HYDROGEN') then
          read(line(start:end),fmt10) iwork
          sub_count = sub_count + 1
          prmtop%dihe_inc_hy(sub_count,count) = iwork

          if (sub_count == 5) then
            count = count + 1
            sub_count = 0
          end if

        else if (keyword == 'DIHEDRALS_WITHOUT_HYDROGEN') then
          read(line(start:end),fmt10) iwork
          sub_count = sub_count + 1
          prmtop%dihe_wo_hy(sub_count,count) = iwork

          if (sub_count == 5) then
            count = count + 1
            sub_count = 0
          end if

        else if (keyword == 'SOLVENT_POINTERS') then
          read(line(start:end),fmt10) iwork

          if (count == 1) then
            prmtop%iptres = iwork
          else if (count == 2) then
            prmtop%nspm = iwork
          else if (count == 3) then
            prmtop%nspsol = iwork
          end if

          count = count + 1

        end if

      end do
    end do
100 continue

    return

  end subroutine read_prmtop_int_format

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_prmtop_real_format
  !> @brief        read int real style information from PRMTOP file
  !! @authors      YSFJ
  !! @param[in]    file    : unit number of PRMTOP file
  !! @param[in]    fmt01   : number variables per line, 10 for10I8 
  !! @param[in]    fmt02   : column length per one variable, 8 for 10I8
  !! @param[in]    keyword : input keyword information
  !! @param[in]    num     : number of variables
  !! @param[inout] prmtop  : prmtop data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_prmtop_real_format(file, fmt01, fmt02, keyword, num, prmtop) 

    ! formal arguments
    integer,                 intent(in)    :: file
    character(*),            intent(in)    :: fmt01
    character(*),            intent(in)    :: fmt02
    character(*),            intent(in)    :: keyword
    integer,                 intent(in)    :: num
    type(s_prmtop),          intent(inout) :: prmtop

    ! local variables
    real(wp)                 :: dwork
    integer                  :: i, j
    integer                  :: len
    integer                  :: count
    integer                  :: start, end
    integer                  :: ifmt01
    integer                  :: ifmt02_int_part
    integer                  :: ifmt02_deci_part
    integer                  :: num_reads
    character(100)           :: line
    character(10)            :: fmt10
    character(10)            :: fmt02_int_part
    character(10)            :: fmt02_deci_part


    read(fmt01,*) ifmt01
    len = len_trim(fmt02)

    do i = 1, len
      if (fmt02(i:i) == '.') then
        read(fmt02(1:i-1),*)   ifmt02_int_part
        read(fmt02(i+1:len),*) ifmt02_deci_part
        read(fmt02(1:i-1),'(A)')   fmt02_int_part
        read(fmt02(i+1:len),'(A)') fmt02_deci_part
        exit
      end if
    end do

    if (ifmt02_int_part < 10) then
      if (ifmt02_deci_part < 10) then
        fmt10 = '(Ex.y)'
        read(fmt10(3:3),'(A1)') fmt02_int_part
        read(fmt10(5:5),'(A1)') fmt02_deci_part
      end if
    else if (ifmt02_int_part < 100) then
      if (ifmt02_deci_part < 10) then
        fmt10 = '(Exx.y)'
        write(fmt10(3:4),'(A2)') fmt02_int_part
        write(fmt10(6:6),'(A1)') fmt02_deci_part
      else if (ifmt02_deci_part < 100) then
        fmt10 = '(Exx.yy)'
        write(fmt10(3:4),'(A2)') fmt02_int_part
        write(fmt10(6:7),'(A2)') fmt02_deci_part
      end if
    end if

    ! Read atom
    !
    len = ifmt01 * ifmt02_int_part

    num_reads = num / ifmt01
    if (mod(num,ifmt01) /= 0) then
      num_reads = num_reads + 1
    end if

    count = 1
    do j = 1, num_reads
      read(file,'(A)',err=100,end=100) line

      do i = 1, ifmt01
        if (count > num) then
          exit
        end if
        start = (i - 1) * ifmt02_int_part + 1
        end = start + ifmt02_int_part - 1

        if (keyword == 'CHARGE') then
          read(line(start:end),fmt10) prmtop%charge(count)
          
        else if (keyword == 'MASS') then
          read(line(start:end),fmt10) prmtop%mass(count)

        else if (keyword == 'BOND_FORCE_CONSTANT') then
          read(line(start:end),fmt10) prmtop%bond_fcons_uniq(count)

        else if (keyword == 'BOND_EQUIL_VALUE') then
          read(line(start:end),fmt10) prmtop%bond_equil_uniq(count)

        else if (keyword == 'ANGLE_FORCE_CONSTANT') then
          read(line(start:end),fmt10) prmtop%angl_fcons_uniq(count)

        else if (keyword == 'ANGLE_EQUIL_VALUE') then
          read(line(start:end),fmt10) prmtop%angl_equil_uniq(count)

        else if (keyword == 'DIHEDRAL_FORCE_CONSTANT') then
          read(line(start:end),fmt10) prmtop%dihe_fcons_uniq(count)

        else if (keyword == 'DIHEDRAL_PERIODICITY') then
          read(line(start:end),fmt10) prmtop%dihe_perio_uniq(count)

        else if (keyword == 'DIHEDRAL_PHASE') then
          read(line(start:end),fmt10) prmtop%dihe_phase_uniq(count)

        else if (keyword == 'SOLTY') then
          read(line(start:end),fmt10) prmtop%solty(count)

        else if (keyword == 'SCEE_SCALE_FACTOR') then
          read(line(start:end),fmt10) prmtop%scee_scale_fact(count)

        else if (keyword == 'SCNB_SCALE_FACTOR') then
          read(line(start:end),fmt10) prmtop%scnb_scale_fact(count)

        else if (keyword == 'LENNARD_JONES_ACOEF') then
          read(line(start:end),fmt10) prmtop%lennarda(count)

        else if (keyword == 'LENNARD_JONES_BCOEF') then
          read(line(start:end),fmt10) prmtop%lennardb(count)

        else if (keyword == 'HBOND_ACOEF') then
          read(line(start:end),fmt10) prmtop%hbond_acoef(count)

        else if (keyword == 'HBOND_BCOEF') then
          read(line(start:end),fmt10) prmtop%hbond_bcoef(count)

        else if (keyword == 'HBCUT') then
          read(line(start:end),fmt10) prmtop%hb_cut(count)

        else if (keyword == 'RADII') then
          read(line(start:end),fmt10) prmtop%radi_born(count)

        else if (keyword == 'SCREEN') then
          read(line(start:end),fmt10) prmtop%fs_born(count)

        else if (keyword == 'BOX_DIMENSIONS') then
          read(line(start:end),fmt10) dwork

          if      (count == 1) then
            prmtop%oldbeta = dwork
          else if (count == 2) then
            prmtop%box(1) = dwork
          else if (count == 3) then
            prmtop%box(2) = dwork
          else if (count == 4) then
            prmtop%box(3) = dwork
          end if

        end if
        count = count + 1

      end do
200   continue

    end do
100 continue

    return

  end subroutine read_prmtop_real_format

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_prmtop
  !> @brief        write prmtop information
  !! @authors      YSFJ
  !! @param[in]    file   : unit number of file
  !! @param[in]    prmtop : prmtop data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_prmtop(file, prmtop)

    ! formal arguments
    integer,                 intent(in) :: file
    type(s_prmtop),          intent(in) :: prmtop

    ! local variable
    integer                  :: i, j
    integer                  :: size
    integer                  :: sub_count


    if (prmtop%lversion) then                   
      write(file,'(A80)') prmtop%version

    end if

    if (prmtop%ltitle) then
      write(file,'(A)') '%FLAG TITLE'
      write(file,'(A)') '%FORMAT(20a4)'
      write(file,'(A)') prmtop%title
    end if

    if (prmtop%lpointers) then                
      write(file,'(A)') '%FLAG POINTERS'
      write(file,'(A)') '%FORMAT(10I8)'
      write(file,'(10I8)') &
        prmtop%num_atoms, prmtop%num_types,  &
        prmtop%num_bondh, prmtop%num_mbonda, &
        prmtop%num_anglh, prmtop%num_mangla, &       
        prmtop%num_diheh, prmtop%num_mdihea, &     
        prmtop%num_hparm, prmtop%num_parm          
  
      write(file,'(10I8)') &
        prmtop%num_nb, prmtop%num_residues,         &
        prmtop%num_bona, prmtop%num_theta,          &
        prmtop%num_phia, prmtop%num_uniqbond,       &
        prmtop%num_uniqangl, prmtop%num_uniqdihe,   &
        prmtop%num_types_prm, prmtop%num_bondtypes

      write(file,'(10I8)') &
        prmtop%ifpert,             prmtop%num_bond_perturbed, &
        prmtop%num_angl_perturbed, prmtop%num_dihe_perturbed, &
        prmtop%num_mbper,          prmtop%num_mgper, &         
        prmtop%num_mdper,          prmtop%ifbox, &             
        prmtop%num_max_res_atoms,  prmtop%ifcap             

      write(file,'(2I8)') &
        prmtop%num_extra_points, prmtop%num_copy          

    end if

    if (prmtop%latom_name) then
      write(file,'(A)') '%FLAG ATOM_NAME'
      write(file,'(A)') '%FORMAT(20a4)'

      do i = 1, prmtop%num_atoms
        if (i == prmtop%num_atoms .or.  mod(i,20) == 0) then
          write(file,'(A4)') prmtop%atom_name(i)
        else
          write(file,'(A4)',advance='no') prmtop%atom_name(i)
        end if
      end do

    end if

    if (prmtop%lcharge) then
      write(file,'(A)') '%FLAG CHARGE'
      write(file,'(A)') '%FORMAT(5E16.8)'

      do i = 1, prmtop%num_atoms
        if (i == prmtop%num_atoms .or. mod(i,5) == 0) then
          write(file,'(ES16.8)') prmtop%charge(i)
        else
          write(file,'(ES16.8)',advance='no') prmtop%charge(i)
        end if
      end do

    end if

    if (prmtop%lmass) then
      write(file,'(A)') '%FLAG MASS'
      write(file,'(A)') '%FORMAT(5E16.8)'

      do i = 1, prmtop%num_atoms
        if (i == prmtop%num_atoms .or. mod(i,5) == 0) then
          write(file,'(ES16.8)') prmtop%mass(i)
        else
          write(file,'(ES16.8)',advance='no') prmtop%mass(i)
        end if
      end do

    end if

    if (prmtop%latom_type_index) then
      write(file,'(A)') '%FLAG ATOM_TYPE_INDEX'
      write(file,'(A)') '%FORMAT(10I8)'

      do i = 1, prmtop%num_atoms
        if (i == prmtop%num_atoms .or. mod(i,10) == 0) then
          write(file,'(I8)') prmtop%atom_cls_no(i)
        else
          write(file,'(I8)',advance='no') prmtop%atom_cls_no(i)
        end if
      end do

    end if

    if (prmtop%lnumber_excluded_atoms) then
      write(file,'(A)') '%FLAG NUMBER_EXCLUDED_ATOMS'
      write(file,'(A)') '%FORMAT(10I8)'

      do i = 1, prmtop%num_atoms
        if (i == prmtop%num_atoms .or. mod(i,10) == 0) then
          write(file,'(I8)') prmtop%excl_atom(i)
        else
          write(file,'(I8)',advance='no') prmtop%excl_atom(i)
        end if
      end do

    end if

    if (prmtop%lnonbonded_parm_index) then
      write(file,'(A)') '%FLAG NONBONDED_PARM_INDEX'
      write(file,'(A)') '%FORMAT(10I8)'

      size = prmtop%num_types * prmtop%num_types

      do i = 1, size
        if (i == size .or. mod(i,10) == 0) then
          write(file,'(I8)') prmtop%nb_par_idx(i)
        else
          write(file,'(I8)',advance='no') prmtop%nb_par_idx(i)
        end if
      end do

    end if

    if (prmtop%lresidue_label) then
      write(file,'(A)') '%FLAG RESIDUE_LABEL'
      write(file,'(A)') '%FORMAT(20a4)'

      do i = 1, prmtop%num_residues
        if (i == prmtop%num_residues .or. mod(i,20) == 0) then
          write(file,'(A4)') prmtop%res_label(i)
        else
          write(file,'(A4)',advance='no') prmtop%res_label(i)
        end if
      end do

    end if

    if (prmtop%lresidue_pointer) then
      write(file,'(A)') '%FLAG RESIDUE_POINTER'
      write(file,'(A)') '%FORMAT(10I8)'

      do i = 1, prmtop%num_residues
        if (i == prmtop%num_residues .or. mod(i,10) == 0) then
          write(file,'(I8)') prmtop%res_point(i)
        else
          write(file,'(I8)',advance='no') prmtop%res_point(i)
        end if
      end do

    end if

    if (prmtop%lbond_force_constant) then
      write(file,'(A)') '%FLAG BOND_FORCE_CONSTANT'
      write(file,'(A)') '%FORMAT(5E16.8)'

      do i = 1, prmtop%num_uniqbond
        if (i == prmtop%num_uniqbond .or. mod(i,5) == 0) then
          write(file,'(ES16.8)') prmtop%bond_fcons_uniq(i)
        else
          write(file,'(ES16.8)',advance='no') prmtop%bond_fcons_uniq(i)
        end if
      end do

    end if

    if (prmtop%lbond_equil_value) then
      write(file,'(A)') '%FLAG BOND_EQUIL_VALUE'
      write(file,'(A)') '%FORMAT(5E16.8)'

      do i = 1, prmtop%num_uniqbond
        if (i == prmtop%num_uniqbond .or. mod(i,5) == 0) then
          write(file,'(ES16.8)') prmtop%bond_equil_uniq(i)
        else
          write(file,'(ES16.8)',advance='no') prmtop%bond_equil_uniq(i)
        end if
      end do

    end if

    if (prmtop%langle_force_constant) then
      write(file,'(A)') '%FLAG ANGLE_FORCE_CONSTANT'
      write(file,'(A)') '%FORMAT(5E16.8)'

      do i = 1, prmtop%num_uniqangl
        if (i == prmtop%num_uniqangl .or. mod(i,5) == 0) then
          write(file,'(ES16.8)') prmtop%angl_fcons_uniq(i)
        else
          write(file,'(ES16.8)',advance='no') prmtop%angl_fcons_uniq(i)
        end if
      end do

    end if

    if (prmtop%langle_equil_value) then
      write(file,'(A)') '%FLAG ANGLE_EQUIL_VALUE'
      write(file,'(A)') '%FORMAT(5E16.8)'

      do i = 1, prmtop%num_uniqangl
        if (i == prmtop%num_uniqangl .or. mod(i,5) == 0) then
          write(file,'(ES16.8)') prmtop%angl_equil_uniq(i)
        else
          write(file,'(ES16.8)',advance='no') prmtop%angl_equil_uniq(i)
        end if
      end do

    end if

    if (prmtop%ldihedral_force_constant) then
      write(file,'(A)') '%FLAG DIHEDRAL_FORCE_CONSTANT'
      write(file,'(A)') '%FORMAT(5E16.8)'

      do i = 1, prmtop%num_uniqdihe
        if (i == prmtop%num_uniqdihe .or. mod(i,5) == 0) then
          write(file,'(ES16.8)') prmtop%dihe_fcons_uniq(i)
        else
          write(file,'(ES16.8)',advance='no') prmtop%dihe_fcons_uniq(i)
        end if
      end do

    end if

    if (prmtop%ldihedral_periodicity) then
      write(file,'(A)') '%FLAG DIHEDRAL_PERIODICITY'
      write(file,'(A)') '%FORMAT(5E16.8)'

      do i = 1, prmtop%num_uniqdihe
        if (i == prmtop%num_uniqdihe .or. mod(i,5) == 0) then
          write(file,'(ES16.8)') prmtop%dihe_perio_uniq(i)
        else
          write(file,'(ES16.8)',advance='no') prmtop%dihe_perio_uniq(i)
        end if
      end do

    end if

    if (prmtop%ldihedral_phase) then
      write(file,'(A)') '%FLAG DIHEDRAL_PHASE'
      write(file,'(A)') '%FORMAT(5E16.8)'

      do i = 1, prmtop%num_uniqdihe
        if (i == prmtop%num_uniqdihe .or. mod(i,5) == 0) then
          write(file,'(ES16.8)') prmtop%dihe_phase_uniq(i)
        else
          write(file,'(ES16.8)',advance='no') prmtop%dihe_phase_uniq(i)
        end if
      end do

    end if

    if (prmtop%lscee_scale_factor) then
      write(file,'(A)') '%FLAG SCEE_SCALE_FACTOR'
      write(file,'(A)') '%FORMAT(5E16.8)'

      do i = 1, prmtop%num_uniqdihe
        if (i == prmtop%num_uniqdihe .or. mod(i,5) == 0) then
          write(file,'(ES16.8)') prmtop%scee_scale_fact(i)
        else
          write(file,'(ES16.8)',advance='no') prmtop%scee_scale_fact(i)
        end if
      end do

    end if

    if (prmtop%lscnb_scale_factor) then
      write(file,'(A)') '%FLAG SCNB_SCALE_FACTOR'
      write(file,'(A)') '%FORMAT(5E16.8)'

      do i = 1, prmtop%num_uniqdihe
        if (i == prmtop%num_uniqdihe .or. mod(i,5) == 0) then
          write(file,'(ES16.8)') prmtop%scnb_scale_fact(i)
        else
          write(file,'(ES16.8)',advance='no') prmtop%scnb_scale_fact(i)
        end if
      end do

    end if

    if (prmtop%lsolty) then
      write(file,'(A)') '%FLAG SOLTY'
      write(file,'(A)') '%FORMAT(5E16.8)'

      do i = 1, prmtop%num_types_prm
        if (i == prmtop%num_types_prm .or. mod(i,5) == 0) then
          write(file,'(ES16.8)') prmtop%solty(i)
        else
          write(file,'(ES16.8)',advance='no') prmtop%solty(i)
        end if
      end do

    end if

    if (prmtop%llennard_jones_acoef) then
      write(file,'(A)') '%FLAG LENNARD_JONES_ACOEF'
      write(file,'(A)') '%FORMAT(5E16.8)'

      size = prmtop%num_types * (prmtop%num_types + 1) * 0.5

      do i = 1, size
        if (i == size .or. mod(i,5) == 0) then
          write(file,'(ES16.8)') prmtop%lennarda(i)
        else
          write(file,'(ES16.8)',advance='no') prmtop%lennarda(i)
        end if
      end do

    end if

    if (prmtop%llennard_jones_bcoef) then
      write(file,'(A)') '%FLAG LENNARD_JONES_BCOEF'
      write(file,'(A)') '%FORMAT(5E16.8)'

      size = prmtop%num_types * (prmtop%num_types + 1) * 0.5

      do i = 1, size
        if (i == size .or. mod(i,5) == 0) then
          write(file,'(ES16.8)') prmtop%lennardb(i)
        else
          write(file,'(ES16.8)',advance='no') prmtop%lennardb(i)
        end if
      end do

    end if

    if (prmtop%lbonds_inc_hydrogen) then
      write(file,'(A)') '%FLAG BONDS_INC_HYDROGEN'
      write(file,'(A)') '%FORMAT(10I8)'

      sub_count = 0

      do i = 1, prmtop%num_bondh

        do j = 1, 3
          sub_count = sub_count + 1

          if (sub_count == 10 .or. &
              (i == prmtop%num_bondh .and. j == 3)) then
            write(file,'(I8)') prmtop%bond_inc_hy(j,i)
            sub_count = 0
          else
            write(file,'(I8)',advance='no') prmtop%bond_inc_hy(j,i)

          end if
        end do
      end do

    end if

    if (prmtop%lbonds_without_hydrogen) then
      write(file,'(A)') '%FLAG BONDS_WITHOUT_HYDROGEN'
      write(file,'(A)') '%FORMAT(10I8)'

      sub_count = 0

      do i = 1, prmtop%num_mbonda

        do j = 1, 3
          sub_count = sub_count + 1

          if (sub_count == 10 .or. &
              (i == prmtop%num_mbonda .and. j == 3)) then
            write(file,'(I8)') prmtop%bond_wo_hy(j,i)
            sub_count = 0
          else
            write(file,'(I8)',advance='no') prmtop%bond_wo_hy(j,i)

          end if
        end do
      end do

    end if

    if (prmtop%langles_inc_hydrogen) then
      write(file,'(A)') '%FLAG ANGLES_INC_HYDROGEN'
      write(file,'(A)') '%FORMAT(10I8)'

      sub_count = 0

      do i = 1, prmtop%num_anglh

        do j = 1, 4
          sub_count = sub_count + 1

          if (sub_count == 10 .or. &
              (i == prmtop%num_anglh .and. j == 4)) then
            write(file,'(I8)') prmtop%angl_inc_hy(j,i)
            sub_count = 0
          else
            write(file,'(I8)',advance='no') prmtop%angl_inc_hy(j,i)

          end if
        end do
      end do

    end if

    if (prmtop%langles_without_hydrogen) then
      write(file,'(A)') '%FLAG ANGLES_WITHOUT_HYDROGEN'
      write(file,'(A)') '%FORMAT(10I8)'

      sub_count = 0

      do i = 1, prmtop%num_mangla

        do j = 1, 4
          sub_count = sub_count + 1

          if (sub_count == 10 .or. &
              (i == prmtop%num_mangla .and. j == 4)) then
            write(file,'(I8)') prmtop%angl_wo_hy(j,i)
            sub_count = 0
          else
            write(file,'(I8)',advance='no') prmtop%angl_wo_hy(j,i)

          end if
        end do
      end do

    end if

    if (prmtop%ldihedrals_inc_hydrogen) then
      write(file,'(A)') '%FLAG DIHEDRALS_INC_HYDROGEN'
      write(file,'(A)') '%FORMAT(10I8)'

      sub_count = 0

      do i = 1, prmtop%num_diheh

        do j = 1, 5
          sub_count = sub_count + 1

          if (sub_count == 10 .or. &
              (i == prmtop%num_diheh .and. j == 5)) then
            write(file,'(I8)') prmtop%dihe_inc_hy(j,i)
            sub_count = 0
          else
            write(file,'(I8)',advance='no') prmtop%dihe_inc_hy(j,i)

          end if
        end do
      end do

    end if

    if (prmtop%ldihedrals_without_hydrogen) then
      write(file,'(A)') '%FLAG DIHEDRALS_WITHOUT_HYDROGEN'
      write(file,'(A)') '%FORMAT(10I8)'

      sub_count = 0

      do i = 1, prmtop%num_mdihea

        do j = 1, 5
          sub_count = sub_count + 1

          if (sub_count == 10 .or. &
              (i == prmtop%num_mdihea .and. j == 5)) then
            write(file,'(I8)') prmtop%dihe_wo_hy(j,i)
            sub_count = 0
          else
            write(file,'(I8)',advance='no') prmtop%dihe_wo_hy(j,i)

          end if
        end do
      end do

    end if

    if (prmtop%lexcluded_atoms_list) then
      write(file,'(A)') '%FLAG EXCLUDED_ATOMS_LIST'
      write(file,'(A)') '%FORMAT(10I8)'

      do i = 1, prmtop%num_nb
        if (i == prmtop%num_nb .or. mod(i,10) == 0) then
          write(file,'(I8)') prmtop%inb(i)
        else
          write(file,'(I8)',advance='no') prmtop%inb(i)
        end if
      end do

    end if

    if (prmtop%lhbond_acoef) then
      write(file,'(A)') '%FLAG HBOND_ACOEF'
      write(file,'(A)') '%FORMAT(5E16.8)'

      do i = 1, prmtop%num_bondtypes
        if (i == prmtop%num_bondtypes .or. mod(i,10) == 0) then
          write(file,'(ES16.8)') prmtop%hbond_acoef(i)
        else
          write(file,'(ES16.8)',advance='no') prmtop%hbond_acoef(i)
        end if
      end do

    end if

    if (prmtop%lhbond_bcoef) then
      write(file,'(A)') '%FLAG HBOND_BCOEF'
      write(file,'(A)') '%FORMAT(5E16.8)'

      do i = 1, prmtop%num_bondtypes
        if (i == prmtop%num_bondtypes .or. mod(i,10) == 0) then
          write(file,'(ES16.8)') prmtop%hbond_bcoef(i)
        else
          write(file,'(ES16.8)',advance='no') prmtop%hbond_bcoef(i)
        end if
      end do

    end if

    if (prmtop%lhbcut) then
      write(file,'(A)') '%FLAG HBCUT'
      write(file,'(A)') '%FORMAT(5E16.8)'

      do i = 1, prmtop%num_bondtypes
        if (i == prmtop%num_bondtypes .or. mod(i,10) == 0) then
          write(file,'(ES16.8)') prmtop%hb_cut(i)
        else
          write(file,'(ES16.8)',advance='no') prmtop%hb_cut(i)
        end if
      end do

    end if

    if (prmtop%lamber_atom_type) then
      write(file,'(A)') '%FLAG AMBER_ATOM_TYPE'
      write(file,'(A)') '%FORMAT(20a4)'

      do i = 1, prmtop%num_atoms
        if (i == prmtop%num_atoms .or.  mod(i,20) == 0) then
          write(file,'(A4)') prmtop%atom_cls_name(i)
        else
          write(file,'(A4)',advance='no') prmtop%atom_cls_name(i)
        end if
      end do

    end if

    if (prmtop%ltree_chain_classification) then
      write(file,'(A)') '%FLAG TREE_CHAIN_CLASSIFICATION'
      write(file,'(A)') '%FORMAT(20a4)'

      do i = 1, prmtop%num_atoms
        if (i == prmtop%num_atoms .or.  mod(i,20) == 0) then
          write(file,'(A4)') prmtop%classify(i)
        else
          write(file,'(A4)',advance='no') prmtop%classify(i)
        end if
      end do

    end if

    if (prmtop%ljoin_array) then
      write(file,'(A)') '%FLAG JOIN_ARRAY'
      write(file,'(A)') '%FORMAT(10I8)'

      do i = 1, prmtop%num_atoms
        if (i == prmtop%num_atoms .or.  mod(i,10) == 0) then
          write(file,'(I8)') prmtop%join_array(i)
        else
          write(file,'(I8)',advance='no') prmtop%join_array(i)
        end if
      end do

    end if

    if (prmtop%lirotat) then
      write(file,'(A)') '%FLAG IROTAT'
      write(file,'(A)') '%FORMAT(10I8)'

      do i = 1, prmtop%num_atoms
        if (i == prmtop%num_atoms .or.  mod(i,10) == 0) then
          write(file,'(I8)') prmtop%irotate(i)
        else
          write(file,'(I8)',advance='no') prmtop%irotate(i)
        end if
      end do

    end if

    if (prmtop%lsolvent_pointers) then
      write(file,'(A)') '%FLAG SOLVENT_POINTERS'
      write(file,'(A)') '%FORMAT(3I8)'
      write(file,'(3I8)') prmtop%iptres, prmtop%nspm, prmtop%nspsol

    end if

    if (prmtop%latoms_per_molecule) then
      write(file,'(A)') '%FLAG ATOMS_PER_MOLECULE'
      write(file,'(A)') '%FORMAT(10I8)'

      do i = 1, prmtop%nspm
        if (i == prmtop%nspm .or. mod(i,10) == 0) then
          write(file,'(I8)') prmtop%nsp(i)
        else
          write(file,'(I8)',advance='no') prmtop%nsp(i)
        end if
      end do

    end if

    if (prmtop%lbox_dimensions) then
      write(file,'(A)') '%FLAG BOX_DIMENSIONS'
      write(file,'(A)') '%FORMAT(5E16.8)'

      write(file,'(4ES16.8)') prmtop%oldbeta,    &
                            prmtop%box(1),     &
                            prmtop%box(2),     &
                            prmtop%box(3)

    end if

    if (prmtop%lradius_set) then
      write(file,'(A)') '%FLAG RADIUS_SET'
      write(file,'(A)') '%FORMAT(1a80)'
      write(file,'(A80)') prmtop%radius_set

    end if

    if (prmtop%lradii) then
      write(file,'(A)') '%FLAG RADII'
      write(file,'(A)') '%FORMAT(5E16.8)'

      do i = 1, prmtop%num_atoms
        if (i == prmtop%num_atoms .or.  mod(i,5) == 0) then
          write(file,'(ES16.8)') prmtop%radi_born(i)
        else
          write(file,'(ES16.8)',advance='no') prmtop%radi_born(i)
        end if
      end do

    end if

    if (prmtop%lscreen) then
      write(file,'(A)') '%FLAG SCREEN'
      write(file,'(A)') '%FORMAT(5E16.8)'

      do i = 1, prmtop%num_atoms
        if (i == prmtop%num_atoms .or.  mod(i,5) == 0) then
          write(file,'(ES16.8)') prmtop%fs_born(i)
        else
          write(file,'(ES16.8)',advance='no') prmtop%fs_born(i)
        end if
      end do

    end if

    return

  end subroutine write_prmtop

end module fileio_prmtop_mod
