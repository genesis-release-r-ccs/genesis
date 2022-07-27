!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_par_mod
!> @brief   read and write a CHARMM-style parameter file
!! @authors Yuji Sugita (YS), Norio Takase (NT), Jaewoon Jung (JJ), 
!!          Naoyuki Miyashita (NM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fileio_par_mod

  use fileio_top_mod
  use fileio_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_par
    
    integer           :: type          = 0

    integer           :: num_bonds     = 0
    integer           :: num_angles    = 0
    integer           :: num_dihedrals = 0
    integer           :: num_impropers = 0
    integer           :: num_atom_cls  = 0
    integer           :: num_nbfix     = 0
    integer           :: num_cmaps     = 0

    ! bond
    character(6),     allocatable  :: bond_atom_cls(:,:)
    real(wp),         allocatable  :: bond_force_const(:)
    real(wp),         allocatable  :: bond_dist_min(:)

    !SPICA ENM
    integer           :: num_enmt      = 0
    integer           :: num_enmp      = 0
    integer,          allocatable  :: bond_enm_param_type(:)
    integer,          allocatable  :: bond_enm_atom_id(:,:)
    integer,          allocatable  :: bond_enm_ij_type(:)
    real(wp),         allocatable  :: bond_enm_force_const(:)
    real(wp),         allocatable  :: bond_enm_dist_min(:)

    integer           :: num_enmu     = 0
    integer           :: num_enmq     = 0
    integer,          allocatable  :: angle_enm_param_type(:)
    integer,          allocatable  :: angle_enm_atom_id(:,:)
    integer,          allocatable  :: angle_enm_ij_type(:)
    real(wp),         allocatable  :: angle_enm_force_const(:)
    real(wp),         allocatable  :: angle_enm_dist_min(:)
    real(wp),         allocatable  :: angle_enm_lj_eps(:)
    real(wp),         allocatable  :: angle_enm_lj_sigma(:)

    ! angle
    character(6),     allocatable  :: angl_atom_cls(:,:)
    real(wp),         allocatable  :: angl_force_const(:)
    real(wp),         allocatable  :: angl_theta_min(:)
    real(wp),         allocatable  :: angl_ub_force_const(:)
    real(wp),         allocatable  :: angl_ub_rmin(:)

    ! dihedral
    character(6),     allocatable  :: dihe_atom_cls(:,:)
    real(wp),         allocatable  :: dihe_force_const(:)
    integer,          allocatable  :: dihe_periodicity(:)
    real(wp),         allocatable  :: dihe_phase(:)

    ! improper
    character(6),     allocatable  :: impr_atom_cls(:,:)
    real(wp),         allocatable  :: impr_force_const(:)
    integer,          allocatable  :: impr_periodicity(:)
    real(wp),         allocatable  :: impr_phase(:)

    ! nonbond-defaults
    character(6),     allocatable  :: nonb_atom_cls(:)
    real(wp),         allocatable  :: nonb_polar(:)
    real(wp),         allocatable  :: nonb_eps(:)
    real(wp),         allocatable  :: nonb_rmin(:)
    real(wp),         allocatable  :: nonb_polar_14(:)
    real(wp),         allocatable  :: nonb_eps_14(:)
    real(wp),         allocatable  :: nonb_rmin_14(:)
    integer,          allocatable  :: nonb_lj_type_repuls(:)
    integer,          allocatable  :: nonb_lj_type_attract(:)

    ! nbfix
    character(6),     allocatable  :: nbfi_atom_cls(:,:)
    real(wp),         allocatable  :: nbfi_eps(:)
    real(wp),         allocatable  :: nbfi_rmin(:)
    real(wp),         allocatable  :: nbfi_eps_14(:)
    real(wp),         allocatable  :: nbfi_rmin_14(:)

    integer,          allocatable  :: nbfi_repuls(:)
    integer,          allocatable  :: nbfi_attract(:)

    ! cmap
    character(6),     allocatable  :: cmap_atom_cls(:,:)
    integer,          allocatable  :: cmap_resolution(:)
    real(wp),         allocatable  :: cmap_data(:,:,:)

  end type s_par

  integer,       parameter :: MAXCOLUMN_PAR   = 150

  ! parameters for allocatable variables
  integer,      public,  parameter :: ParBond = 1
  integer,      public,  parameter :: ParAngl = 2
  integer,      public,  parameter :: ParDihe = 3
  integer,      public,  parameter :: ParImpr = 4
  integer,      public,  parameter :: ParNbon = 5
  integer,      public,  parameter :: ParNbfi = 6
  integer,      public,  parameter :: ParCmap = 7

  integer,      public,  parameter :: ParENMT = 8
  integer,      public,  parameter :: ParENMP = 9
  integer,      public,  parameter :: ParENMU = 10
  integer,      public,  parameter :: ParENMQ = 11

  ! parameters
  integer,      public,  parameter :: ParTypeCHARMM = 1
  integer,      public,  parameter :: ParTypeXPLOR  = 2

  character(6), public,  parameter :: WildCard      = 'X     '

  integer,      private, parameter :: IPRE_READ     = 1
  integer,      private, parameter :: IREAD         = 2

  ! local variables
  logical,                 private :: vervose = .true. 

  ! subroutines
  public  :: input_par
  public  :: output_par
  public  :: init_par
  public  :: alloc_par
  public  :: dealloc_par
  public  :: dealloc_par_all
  public  :: read_par
  private :: read_par_bond
  private :: read_par_angl
  private :: read_par_dihe
  private :: read_par_impr
  private :: read_par_nonb
  private :: read_par_nbfi
  private :: read_par_cmap
  private :: skip_toppar_resi
  private :: write_par
  private :: write_par_bond
  private :: write_par_angl
  private :: write_par_dihe
  private :: write_par_impr
  private :: write_par_nonb
  private :: write_par_nbfi
  private :: write_par_cmap
  public  :: merge_par
  private :: copy_par
  public  :: resolve_wildcard

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_par
  !> @brief        a driver subroutine for reading CHARMM par file
  !! @authors      YS, NT
  !! @param[in]    par_filename : filename of CHARMM parameter file
  !! @param[out]   par          : CHARMM PAR information
  !! @param[in]    top          : CHARMM TOP information (optional)
  !! @note         TOP information needs resolving nonbonded wild-card atom name
  !!               in CHARMM19 force fields
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_par(par_filename, par, top)

    ! formal arguments
    character(*),            intent(in)    :: par_filename
    type(s_par),             intent(inout) :: par
    type(s_top),   optional, intent(in)    :: top

    ! local variables
    type(s_par)              :: par0
    integer                  :: file, i
    character(MaxFilename)   :: filename


    call init_par(par)

    i = 0
    do while(extract(par_filename, i, filename))

      ! open CHARMM parameter file
      !
      call open_file(file, filename, IOFileInput)

      ! read CHARMM parameter file
      !
      call read_par(file, par0)

      ! close CHARMM parameter file
      !
      call close_file(file)
 
      ! resolve CHARMM19 nonbonded wild-card atom name
      !
      call resolve_wildcard(par0, top)

      ! merge parameter information
      !
      call merge_par(par0, par)

    end do

    if (main_rank .and. vervose) then
      write(MsgOut,'(A)') 'Input_Par> Summary of Parfile'
      write(MsgOut,'(A20,I10,A20,I10)')               &
           '  num_bonds       = ', par%num_bonds,     & 
           '  num_angles      = ', par%num_angles
      write(MsgOut,'(A20,I10,A20,I10)')               &
           '  num_dihedrals   = ', par%num_dihedrals, &
           '  num_impropers   = ', par%num_impropers
      write(MsgOut,'(A20,I10,A20,I10)')               & 
           '  num_atom_cls    = ', par%num_atom_cls,  &
           '  num_nbfix       = ', par%num_nbfix
      write(MsgOut,'(A20,I10)')                       & 
           '  num_cmap_terms  = ', par%num_cmaps
      write(MSgOut,'(A)') ' '
      vervose = .false.
    end if

    return

  end subroutine input_par

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_par
  !> @brief        a driver subroutine for writing CHARMM par file
  !! @authors      YS
  !! @param[in]    par_filename : filename of CHARMM parameter file
  !! @param[in]    par          : CHARMM PAR information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_par(par_filename, par)

    ! formal arguments
    character(*),            intent(in)    :: par_filename
    type(s_par),             intent(in)    :: par

    ! local variables
    integer                  :: file


    ! open CHARMM parameter file
    !
    call open_file(file, par_filename, IOFileOutputNew)

    ! write CHARMM parameter file
    !
    call write_par(file, par)

    ! close CHARMM parameter file
    !
    call close_file(file)
 
    return

  end subroutine output_par

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_par
  !> @brief        initialize CHARMM PAR information
  !! @authors      YS, NT
  !! @param[out]   par : CHARMM PAR information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_par(par)

    ! formal arguments
    type(s_par),             intent(inout) :: par


    par%type          = 0
    par%num_bonds     = 0 
    par%num_angles    = 0
    par%num_dihedrals = 0
    par%num_impropers = 0
    par%num_atom_cls  = 0
    par%num_nbfix     = 0
    par%num_cmaps     = 0

    par%num_enmt = 0
    par%num_enmp = 0

    par%num_enmu = 0
    par%num_enmq = 0

    return

  end subroutine init_par

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_par
  !> @brief        allocate CHARMM PAR information
  !! @authors      YS, NT
  !! @param[inout] par      : CHARMM PAR information
  !! @param[in]    variable : a selected variable
  !! @param[in]    var_size : size of the selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_par(par, variable, var_size)

    ! formal arguments
    type(s_par),             intent(inout) :: par
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size
  
    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat

  
    alloc_stat = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case(ParBond)

      if (allocated(par%bond_atom_cls)) then
        if (size(par%bond_atom_cls(1,:)) == var_size) &
          return
        deallocate(par%bond_atom_cls,    &
                   par%bond_force_const, &
                   par%bond_dist_min,    &
                   stat = dealloc_stat)
      end if

      allocate(par%bond_atom_cls(2, var_size), &
               par%bond_force_const(var_size), &
               par%bond_dist_min(var_size),    &
               stat = alloc_stat)

      par%bond_atom_cls(1:2, 1:var_size) = ''
      par%bond_force_const  (1:var_size) = 0.0_wp
      par%bond_dist_min     (1:var_size) = 0.0_wp

    case(ParENMT)
      if (allocated(par%bond_enm_ij_type)) then
        if (size(par%bond_enm_ij_type(:)) == var_size) &
          return
        deallocate(par%bond_enm_atom_id,    &
                   par%bond_enm_ij_type,    &
                   stat = dealloc_stat)
      end if

      allocate(par%bond_enm_atom_id(2,var_size),  &
               par%bond_enm_ij_type(var_size),    &
               stat = alloc_stat)

      par%bond_enm_atom_id(1:2,1:var_size) = -1
      par%bond_enm_ij_type(1:var_size) = -1

    case(ParENMP)
      if (allocated(par%bond_enm_param_type)) then
        if (size(par%bond_enm_param_type(:)) == var_size) &
          return
        deallocate(par%bond_enm_param_type,    &
                   par%bond_enm_force_const,   &
                   par%bond_enm_dist_min,      &
                   stat = dealloc_stat)
      end if


      allocate(par%bond_enm_param_type(var_size),   &
               par%bond_enm_force_const(var_size),  &
               par%bond_enm_dist_min(var_size),     &
               stat = alloc_stat)


      par%bond_enm_param_type(1:var_size) = -1
      par%bond_enm_force_const(1:var_size) = 0.0_wp
      par%bond_enm_dist_min(1:var_size) = 0.0_wp

    case(ParENMU)
      if (allocated(par%angle_enm_ij_type)) then
        if (size(par%angle_enm_ij_type(:)) == var_size) &
          return
        deallocate(par%angle_enm_atom_id,    &
                   par%angle_enm_ij_type,    &
                   stat = dealloc_stat)
      end if

      allocate(par%angle_enm_atom_id(3,var_size),  &
               par%angle_enm_ij_type(var_size),    &
               stat = alloc_stat)

      par%angle_enm_atom_id(1:3,1:var_size) = -1
      par%angle_enm_ij_type(1:var_size) = -1

    case(ParENMQ)
      if (allocated(par%angle_enm_param_type)) then
        if (size(par%angle_enm_param_type(:)) == var_size) &
          return
        deallocate(par%angle_enm_param_type,    &
                   par%angle_enm_force_const,   &
                   par%angle_enm_dist_min,      &
                   par%angle_enm_lj_eps,        &
                   par%angle_enm_lj_sigma,      &
                   stat = dealloc_stat)
      end if

      allocate(par%angle_enm_param_type(var_size),  &
               par%angle_enm_force_const(var_size), &
               par%angle_enm_dist_min(var_size),    &
               par%angle_enm_lj_eps(var_size),      &
               par%angle_enm_lj_sigma(var_size),    &
               stat = alloc_stat)


      par%angle_enm_param_type(1:var_size) = -1
      par%angle_enm_force_const(1:var_size) = 0.0_wp
      par%angle_enm_dist_min(1:var_size) = 0.0_wp
      par%angle_enm_lj_eps(1:var_size) = 0.0_wp
      par%angle_enm_lj_sigma(1:var_size) = 0.0_wp

    case(ParAngl)

      if (allocated(par%angl_atom_cls)) then
        if (size(par%angl_atom_cls(1,:)) == var_size) &
          return
        deallocate(par%angl_atom_cls,       &
                   par%angl_force_const,    &
                   par%angl_theta_min,      &
                   par%angl_ub_force_const, &
                   par%angl_ub_rmin,        &
                   stat = dealloc_stat)
      end if

      allocate(par%angl_atom_cls(3, var_size),    &
               par%angl_force_const(var_size),    &
               par%angl_theta_min(var_size),      &
               par%angl_ub_force_const(var_size), &
               par%angl_ub_rmin(var_size),        &
               stat = alloc_stat)

      par%angl_atom_cls (1:3, 1:var_size) = ''
      par%angl_force_const   (1:var_size) = 0.0_wp
      par%angl_theta_min     (1:var_size) = 0.0_wp
      par%angl_ub_force_const(1:var_size) = 0.0_wp
      par%angl_ub_rmin       (1:var_size) = 0.0_wp

    case(ParDihe)

      if (allocated(par%dihe_atom_cls)) then
        if (size(par%dihe_atom_cls(1,:)) == var_size) &
          return
        deallocate(par%dihe_atom_cls,    &
                   par%dihe_force_const, &
                   par%dihe_periodicity, &
                   par%dihe_phase,       &
                   stat = dealloc_stat)
      end if

      allocate(par%dihe_atom_cls(4, var_size), &
               par%dihe_force_const(var_size), &
               par%dihe_periodicity(var_size), &
               par%dihe_phase(var_size),       &
               stat = alloc_stat)

      par%dihe_atom_cls(1:4, 1:var_size) = ''
      par%dihe_force_const  (1:var_size) = 0.0_wp
      par%dihe_periodicity  (1:var_size) = 0
      par%dihe_phase        (1:var_size) = 0.0_wp

    case(ParImpr)
      
      if (allocated(par%impr_atom_cls)) then
        if (size(par%impr_atom_cls(1,:)) == var_size) &
          return
        deallocate(par%impr_atom_cls,    &
                   par%impr_force_const, &
                   par%impr_periodicity, &
                   par%impr_phase,       &
                   stat = dealloc_stat)
      end if

      allocate(par%impr_atom_cls(4, var_size),  &
               par%impr_force_const(var_size),  &
               par%impr_periodicity(var_size),  &
               par%impr_phase(var_size),        &
               stat = alloc_stat)

      par%impr_atom_cls(1:4, 1:var_size) = ''
      par%impr_force_const  (1:var_size) = 0.0_wp
      par%impr_periodicity  (1:var_size) = 0
      par%impr_phase        (1:var_size) = 0.0_wp

    case(ParNbon)

      if (allocated(par%nonb_atom_cls)) then
        if (size(par%nonb_atom_cls) == var_size) &
             return
        deallocate(par%nonb_atom_cls,          &
                   par%nonb_polar,             &
                   par%nonb_eps,               &
                   par%nonb_rmin,              &
                   par%nonb_polar_14,          &
                   par%nonb_eps_14,            &
                   par%nonb_rmin_14,           &
                   par%nonb_lj_type_repuls,    &
                   par%nonb_lj_type_attract,   &
                   stat = dealloc_stat)
      end if

      allocate(par%nonb_atom_cls(var_size),          &
               par%nonb_polar(var_size),             &
               par%nonb_eps(var_size),               &
               par%nonb_rmin(var_size),              &
               par%nonb_polar_14(var_size),          &
               par%nonb_eps_14(var_size),            &
               par%nonb_rmin_14(var_size),           &
               par%nonb_lj_type_repuls(var_size),    &
               par%nonb_lj_type_attract(var_size),   &
               stat = alloc_stat)

      par%nonb_atom_cls(1:var_size) = ''
      par%nonb_polar   (1:var_size) = 0.0_wp
      par%nonb_eps     (1:var_size) = 0.0_wp
      par%nonb_rmin    (1:var_size) = 0.0_wp
      par%nonb_polar_14(1:var_size) = 0.0_wp
      par%nonb_eps_14  (1:var_size) = 0.0_wp
      par%nonb_rmin_14 (1:var_size) = 0.0_wp
      par%nonb_lj_type_repuls(1:var_size)  = -1
      par%nonb_lj_type_attract(1:var_size) = -1

    case(ParNbfi)

      if (allocated(par%nbfi_atom_cls)) then
        if (size(par%nbfi_atom_cls) == var_size) &
             return
        deallocate(par%nbfi_atom_cls,  &
                   par%nbfi_eps,       &
                   par%nbfi_rmin,      &
                   par%nbfi_eps_14,    &
                   par%nbfi_rmin_14,   &
                   par%nbfi_repuls,    &
                   par%nbfi_attract,   &
                   stat = dealloc_stat)
      end if

      allocate(par%nbfi_atom_cls(2,var_size), &
               par%nbfi_eps(var_size),        &
               par%nbfi_rmin(var_size),       &
               par%nbfi_eps_14(var_size),     &
               par%nbfi_rmin_14(var_size),    &
               par%nbfi_repuls(var_size),     &
               par%nbfi_attract(var_size),    &
               stat = alloc_stat)

      par%nbfi_atom_cls(1:2,1:var_size) = ''
      par%nbfi_eps         (1:var_size) = 0.0_wp
      par%nbfi_rmin        (1:var_size) = 0.0_wp
      par%nbfi_eps_14      (1:var_size) = 0.0_wp
      par%nbfi_rmin_14     (1:var_size) = 0.0_wp
      par%nbfi_repuls      (1:var_size) = -1
      par%nbfi_attract     (1:var_size) = -1

    case(ParCmap)

      if (allocated(par%cmap_atom_cls)) then
        if (size(par%cmap_atom_cls(1,:)) == var_size) &
          return
        deallocate(par%cmap_atom_cls,   &
                   par%cmap_resolution, &
                   par%cmap_data,       &
                   stat = dealloc_stat)
      end if

      allocate(par%cmap_atom_cls(8,var_size), &
               par%cmap_resolution(var_size), &
               par%cmap_data(24,24,var_size), &
               stat = alloc_stat)

      par%cmap_atom_cls  (1:8,1:var_size) = ''
      par%cmap_resolution    (1:var_size) = 0
      par%cmap_data(1:24,1:24,1:var_size) = 0.0_wp

    case default

      call error_msg('Alloc_Par> bad variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_par

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_par
  !> @brief        deallocate CHARMM PAR information
  !! @authors      YS, NT
  !! @param[inout] par      : CHARMM PAR information
  !! @param[in]    variable : a selected variables
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_par(par, variable)

    ! formal arguments
    type(s_par),             intent(inout) :: par
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat
   

    dealloc_stat = 0

    select case (variable)
      
    case(ParBond)

      if (allocated(par%bond_atom_cls)) then
        deallocate( par%bond_atom_cls,    &
                    par%bond_force_const, &
                    par%bond_dist_min,    &
                    stat = dealloc_stat)
      end if

    case(ParENMT)

      if (allocated(par%bond_enm_ij_type)) then
        deallocate(par%bond_enm_atom_id,    &
                   par%bond_enm_ij_type,    &
                   stat = dealloc_stat)
      end if

    case(ParENMP)

      if (allocated(par%bond_enm_param_type)) then
        deallocate(par%bond_enm_param_type,    &
                   par%bond_enm_force_const,    &
                   par%bond_enm_dist_min,    &
                   stat = dealloc_stat)
      end if

    case(ParENMU)

      if (allocated(par%angle_enm_ij_type)) then
        deallocate(par%angle_enm_atom_id,    &
                   par%angle_enm_ij_type,    &
                   stat = dealloc_stat)
      end if

    case(ParENMQ)

      if (allocated(par%angle_enm_param_type)) then
        deallocate(par%angle_enm_param_type,    &
                   par%angle_enm_force_const,    &
                   par%angle_enm_dist_min,    &
                   par%angle_enm_lj_eps,    &
                   par%angle_enm_lj_sigma,    &
                   stat = dealloc_stat)
      end if

    case(ParAngl)

      if (allocated(par%angl_atom_cls)) then
        deallocate( par%angl_atom_cls,       &
                    par%angl_force_const,    &
                    par%angl_theta_min,      &
                    par%angl_ub_force_const, &
                    par%angl_ub_rmin,        &
                    stat = dealloc_stat)
      end if

    case(ParDihe)

      if (allocated(par%dihe_atom_cls)) then
        deallocate( par%dihe_atom_cls,    &
                    par%dihe_force_const, &
                    par%dihe_periodicity, &
                    par%dihe_phase,       &
                    stat = dealloc_stat)
      end if

    case(ParImpr)

      if (allocated(par%impr_atom_cls)) then
        deallocate( par%impr_atom_cls,    &
                    par%impr_force_const, &
                    par%impr_periodicity, &
                    par%impr_phase,       &
                    stat = dealloc_stat)
      end if

    case(ParNbon)

      if (allocated(par%nonb_atom_cls)) then
        deallocate( par%nonb_atom_cls,        &
                    par%nonb_polar,           &
                    par%nonb_eps,             &
                    par%nonb_rmin,            &
                    par%nonb_polar_14,        &
                    par%nonb_eps_14,          &
                    par%nonb_rmin_14,         &
                    par%nonb_lj_type_repuls,  &
                    par%nonb_lj_type_attract, &
                    stat = dealloc_stat)
      end if

    case(ParNbfi)

      if (allocated(par%nbfi_atom_cls)) then
        deallocate( par%nbfi_atom_cls,  &
                    par%nbfi_eps,       &
                    par%nbfi_rmin,      &
                    par%nbfi_eps_14,    &
                    par%nbfi_rmin_14,   &
                    par%nbfi_repuls,    &
                    par%nbfi_attract,   &
                    stat = dealloc_stat)
      end if

    case(ParCmap)

      if (allocated(par%cmap_atom_cls)) then
        deallocate( par%cmap_atom_cls,   &
                    par%cmap_resolution, &
                    par%cmap_data,       &
                    stat = dealloc_stat)
      end if

    case default

      call error_msg('Dealloc_Par> bad variable')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_par

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_par_all
  !> @brief        deallocate PAR all information
  !! @authors      YS, NT
  !! @param[inout] par :  CHARMM PAR information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_par_all(par)

    ! formal arguments
    type(s_par),             intent(inout) :: par


    call dealloc_par(par, ParBond)
    call dealloc_par(par, ParENMP)
    call dealloc_par(par, ParENMT)
    call dealloc_par(par, ParAngl)
    call dealloc_par(par, ParDihe)
    call dealloc_par(par, ParImpr)
    call dealloc_par(par, ParNbon)
    call dealloc_par(par, ParNbfi)
    call dealloc_par(par, ParCmap)

    return

  end subroutine dealloc_par_all

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_par
  !> @brief        subroutine for reading CHARMM par file
  !! @authors      YS, NT, JJ, NM
  !! @param[in]    file : file unit number of CHARMM parameter file
  !! @param[out]   par  : CHARMM PAR information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_par(file, par)

    ! parameter
!    integer,       parameter :: MAXROW = 80
    integer,       parameter :: MAXKEY = 18

    integer,       parameter :: IBOND= 1,IANGL= 2,ITHET= 3,IDIHE= 4,IPHI = 5
    integer,       parameter :: IIMPR= 6,IIMPH= 7,INONB= 8,INBFI= 9,ICMAP=10
    integer,       parameter :: IHBON=11,IRESI=12,IPRES=13,IEND =14
    integer,       parameter :: ENMT =15,ENMP =16,ENMU =17,ENMQ =18  ! SPICA
    character(4),  parameter :: keywd(MAXKEY) = &
                                      (/'BOND','ANGL','THET','DIHE','PHI ', &
                                        'IMPR','IMPH','NONB','NBFI','CMAP', &
                                        'HBON','RESI','PRES','END ',        &
                                        'ENMT','ENMP','ENMU','ENMQ'/)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_par),             intent(out)   :: par

    ! local variables
    character(MAXCOLUMN_PAR) :: line
    character(4)             :: ckeywd
    character(6)             :: frmt
    integer                  :: i, j, ikeywd, nsta, nend, nchar, read_cnt
    integer                  :: nb, na, nd, ni, nn, nx, nc, mapi, mapi2
    integer                  :: nbenmt, nbenmp, nbenmu, nbenmq


    ! Check par file type
    !
    par%type = ParTypeCHARMM

    frmt = "      "
    if (MAXCOLUMN_PAR < 100) then
      write(frmt(1:5),'(A2,I2,A)') '(A',MAXCOLUMN_PAR,')'
    else
      write(frmt,'(A2,I3,A)') '(A',MAXCOLUMN_PAR,')'
    end if


    !  Read parameter file
    !
    do i = IPRE_READ, IREAD

      !  Initialize variables
      !
      ikeywd = 0
      nb = 0
      na = 0
      nd = 0
      ni = 0
      nn = 0
      nx = 0
      nc = 0
      read_cnt = 0

      nbenmt = 0
      nbenmp = 0
      nbenmu = 0
      nbenmq = 0

      do while(.true.)

        read(file, frmt) line
        read_cnt = read_cnt + 1

        call char_line(MAXCOLUMN_PAR, line, nchar)

        if (nchar <= 0) cycle

        !  Select section in the parameter file
        !
        nsta = 1
        nend = nchar

        call read_word(line, nsta, nend, 4, ckeywd)

        call toupper(ckeywd)

        do j = 1, MAXKEY
          if (ckeywd == keywd(j)) then
            ikeywd = j
            exit
          end if
        end do

        !  Read parameters
        !
        if (ikeywd == IBOND) then
          call read_par_bond(line, nchar, i, par, nb)

        else if (ikeywd == IANGL) then
          call read_par_angl(line, nchar, i, par, na)

        else if (ikeywd == ITHET) then
          call read_par_angl(line, nchar, i, par, na)

        else if (ikeywd == IDIHE) then
          call read_par_dihe(line, nchar, i, par, nd)

        else if (ikeywd == IPHI ) then
          call read_par_dihe(line, nchar, i, par, nd)

        else if (ikeywd == IIMPR) then
          call read_par_impr(line, nchar, i, par, ni)

        else if (ikeywd == IIMPH) then
          call read_par_impr(line, nchar, i, par, ni)

        else if (ikeywd == INONB) then
          call read_par_nonb(line, nchar, i, par, nn)

        else if (ikeywd == INBFI) then
          call read_par_nbfi(line, nchar, i, par, nx)

        else if (ikeywd == ICMAP) then
          call read_par_cmap(line, nchar, i, par, nc, mapi, mapi2)

        else if (ikeywd == IHBON) then
          ! Not implemented

        else if (ikeywd == IRESI .or. ikeywd == IPRES) then
          ikeywd = 0
          call skip_toppar_resi(file, i)

        else if (ikeywd == ENMT) then
           call read_par_enmt(line, nchar, i, par, nbenmt)

        else if (ikeywd == ENMP) then
           call read_par_enmp(line, nchar, i, par, nbenmp)

        else if (ikeywd == ENMU) then
           call read_par_enmu(line, nchar, i, par, nbenmu)

        else if (ikeywd == ENMQ) then
           call read_par_enmQ(line, nchar, i, par, nbenmq)

        else if (ikeywd == IEND) then

          if (i == IPRE_READ) then

            call init_par (par)
            call alloc_par(par, ParBond, nb)
            call alloc_par(par, ParAngl, na)
            call alloc_par(par, ParDihe, nd)
            call alloc_par(par, ParImpr, ni)
            call alloc_par(par, ParNbon, nn)
            call alloc_par(par, ParNbfi, nx)
            call alloc_par(par, ParCmap, nc)
            call alloc_par(par, ParENMT, nbenmt)
            call alloc_par(par, ParENMP, nbenmp)
            call alloc_par(par, ParENMU, nbenmu)
            call alloc_par(par, ParENMQ, nbenmq)

            do j = 1, read_cnt
              backspace(file)
            end do

          end if

          exit

        else
          !  Do Nothing (just Skipping header)
        end if

      end do

    end do

    par%num_bonds     = size(par%bond_atom_cls(1,:))
    par%num_enmt      = size(par%bond_enm_ij_type(:))
    par%num_enmp      = size(par%bond_enm_param_type(:))
    par%num_enmu      = size(par%angle_enm_ij_type(:))
    par%num_enmq      = size(par%angle_enm_param_type(:))
    par%num_angles    = size(par%angl_atom_cls(1,:))
    par%num_dihedrals = size(par%dihe_atom_cls(1,:))
    par%num_impropers = size(par%impr_atom_cls(1,:))
    par%num_atom_cls  = size(par%nonb_atom_cls(:))
    par%num_nbfix     = size(par%nbfi_atom_cls(1,:))
    par%num_cmaps     = size(par%cmap_atom_cls(1,:))


#ifdef DEBUG
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Par> Summary of Parfile'
      write(MsgOut,'(A20,I10,A20,I10)')               &
           '  num_bonds       = ', par%num_bonds,     & 
           '  num_angles      = ', par%num_angles
      write(MsgOut,'(A20,I10,A20,I10)')               &
           '  num_dihedrals   = ', par%num_dihedrals, &
           '  num_impropers   = ', par%num_impropers
      write(MsgOut,'(A20,I10,A20,I10)')               & 
           '  num_atom_cls    = ', par%num_atom_cls,  &
           '  num_nbfix       = ', par%num_nbfix
      write(MsgOut,'(A20,I10)')                       & 
           '  num_cmap_terms  = ', par%num_cmaps
      write(MSgOut,'(A)') ' '
    end if
#endif

    return

  end subroutine read_par

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_par_bond
  !> @brief        read bond information from PAR file
  !! @authors      YS
  !! @param[in]    line  : read file line
  !! @param[in]    nchar : number of characters to read
  !! @param[in]    mode  : mode (IPRE_READ or IREAD)
  !! @param[inout] par   : CHARMM PAR information
  !! @param[inout] nb    : number of read bonds
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_par_bond(line, nchar, mode, par, nb)

    ! formal arguments
    character(*),            intent(in)    :: line
    integer,                 intent(in)    :: nchar
    integer,                 intent(in)    :: mode
    type(s_par),             intent(inout) :: par
    integer,                 intent(inout) :: nb

    ! local variables
    real(wp)                 :: rbond1, rbond2
    integer                  :: nsta, nend, ndata
    character(6)             :: catm(2)


    nsta = 1
    nend = nchar

    call read_ndata(line, nsta, nend, ndata)

    if (ndata == 1) then
      
      !  Do nothing
      !
      
    else if (ndata == 4) then
      
      nb = nb + 1

      if (mode == IREAD) then

        !  Read information: lb(1), lb2(2), kb, b0
        !
        nsta = 1
        nend = nchar

        call read_word(line, nsta, nend, 6, catm(1))
        call read_word(line, nsta, nend, 6, catm(2))
        call read_real(line, nsta, nend, rbond1)
        call read_real(line, nsta, nend, rbond2)

        !  Store parameters
        !
        par%bond_atom_cls(1:2, nb) = catm(1:2)
        par%bond_force_const(nb)   = rbond1
        par%bond_dist_min(nb)      = rbond2

      end if

    else

      call error_msg('  Read_Par_Bond> ERROR: Format is not correct')

    end if

    return

  end subroutine read_par_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_par_angl
  !> @brief        read angle information from PAR file
  !! @authors      YS
  !! @param[in]    line  : read file line
  !! @param[in]    nchar : number of characters to read
  !! @param[in]    mode  : mode (IPRE_READ or IREAD)
  !! @param[inout] par   : CHARMM PAR information
  !! @param[inout] na    : number of read angles
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_par_angl(line, nchar, mode, par, na)

    ! formal arguments
    character(*),            intent(in)    :: line
    integer,                 intent(in)    :: nchar
    integer,                 intent(in)    :: mode
    type(s_par),             intent(inout) :: par
    integer,                 intent(inout) :: na

    ! local variables
    real(wp)                 :: rangl1, rangl2, rurey1, rurey2
    integer                  :: nsta, nend, ndata
    character(6)             :: catm(3)


    nsta = 1
    nend = nchar

    call read_ndata(line, nsta, nend, ndata)

    if (ndata == 1) then

      !  Do nothing
      !
    else if (ndata == 5) then

      na = na + 1

      if (mode == IREAD) then

        !  Read information: ia(1), ia(2), ia(3), ka, a0
        !
        nsta = 1
        nend = nchar

        call read_word(line, nsta, nend, 6, catm(1))
        call read_word(line, nsta, nend, 6, catm(2))
        call read_word(line, nsta, nend, 6, catm(3))
        call read_real(line, nsta, nend, rangl1)
        call read_real(line, nsta, nend, rangl2)

        !  Store the parameter
        !
        par%angl_atom_cls(1:3, na)  = catm(1:3)
        par%angl_force_const(na)    = rangl1
        par%angl_theta_min(na)      = rangl2
        par%angl_ub_force_const(na) = 0.0_wp
        par%angl_ub_rmin(na)        = 0.0_wp

      end if

    else if (ndata == 7) then
      
      na = na + 1

      if (mode == IREAD) then

        !  Read information: ia(1), ia(2), ia(3), ka, a0
        !
        nsta = 1
        nend = nchar

        call read_word(line, nsta, nend, 6, catm(1))
        call read_word(line, nsta, nend, 6, catm(2))
        call read_word(line, nsta, nend, 6, catm(3))
        call read_real(line, nsta, nend, rangl1)
        call read_real(line, nsta, nend, rangl2)
        call read_real(line, nsta, nend, rurey1)
        call read_real(line, nsta, nend, rurey2)


        !  Store the parameter
        !
        par%angl_atom_cls(1:3, na)  = catm(1:3)
        par%angl_force_const(na)    = rangl1
        par%angl_theta_min(na)      = rangl2
        par%angl_ub_force_const(na) = rurey1
        par%angl_ub_rmin(na)        = rurey2

      end if

    else

      call error_msg(' Read_Par_Angl> ERROR: Format is not correct')

    end if

    return

  end subroutine read_par_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_par_dihe
  !> @brief        read dihedral angle information from PAR file
  !! @authors      YS
  !! @param[in]    line  : read file line
  !! @param[in]    nchar : number of characters to read
  !! @param[in]    mode  : mode (IPRE_READ or IREAD)
  !! @param[inout] par   : CHARMM PAR information
  !! @param[inout] nd    : number of read dihedrals
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_par_dihe(line, nchar, mode, par, nd)

    ! formal arguments
    character(*),            intent(in)    :: line
    integer,                 intent(in)    :: nchar
    integer,                 intent(in)    :: mode
    type(s_par),             intent(inout) :: par
    integer,                 intent(inout) :: nd

    ! local variables
    real(wp)                 :: rdihe1, rdihe2
    integer                  :: nsta, nend, ndata, idihe1
    character(6)             :: catm(4)


    nsta = 1
    nend = nchar

    call read_ndata(line, nsta, nend, ndata)

    if (ndata == 1) then

      !  Do nothing
      !
    else if (ndata == 7) then

      nd = nd + 1

      if (mode == IREAD) then

        !  Read information: id(1), id(2), id(3), id(4), kd, md, dd
        !
        nsta = 1
        nend = nchar

        call read_word(line, nsta, nend, 6, catm(1))
        call read_word(line, nsta, nend, 6, catm(2))
        call read_word(line, nsta, nend, 6, catm(3))
        call read_word(line, nsta, nend, 6, catm(4))
        call read_real(line, nsta, nend, rdihe1)
        call read_int (line, nsta, nend, idihe1)
        call read_real(line, nsta, nend, rdihe2)

        !  Store the parameters
        !
        par%dihe_atom_cls(1:4, nd) = catm(1:4)
        par%dihe_force_const(nd)   = rdihe1
        par%dihe_periodicity(nd)   = idihe1
        par%dihe_phase(nd)         = rdihe2

      end if

    else

      call error_msg(' Read_Par_Dihe> ERROR: Format is not correct')

    end if

    return

  end subroutine read_par_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_par_impr
  !> @brief        read improper angle information from PAR file
  !! @authors      YS
  !! @param[in]    line  : read file line
  !! @param[in]    nchar : number of characters to read
  !! @param[in]    mode  : mode (IPRE_READ or IREAD)
  !! @param[inout] par   : CHARMM PAR information
  !! @param[inout] ni    : number of read impropers
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_par_impr(line, nchar, mode, par, ni)

    ! formal arguments
    character(*),            intent(in)    :: line
    integer,                 intent(in)    :: nchar
    integer,                 intent(in)    :: mode
    type(s_par),             intent(inout) :: par
    integer,                 intent(inout) :: ni

    ! local variables
    real(wp)                 :: rimpr1, rimpr2
    integer                  :: nsta, nend, ndata, iimpr1
    character(6)             :: catm(4)


    nsta = 1
    nend = nchar

    call read_ndata(line, nsta, nend, ndata)

    if (ndata == 1) then
      
      !  Do nothing
      !
    else if (ndata == 7) then

      ni = ni + 1

      if (mode == IREAD) then

        !  Read information: id(1), id(2), id(3), id(4), kd, md, dd
        !
        nsta = 1
        nend = nchar

        call read_word(line, nsta, nend, 6, catm(1))
        call read_word(line, nsta, nend, 6, catm(2))
        call read_word(line, nsta, nend, 6, catm(3))
        call read_word(line, nsta, nend, 6, catm(4))
        call read_real(line, nsta, nend, rimpr1)
        call read_int (line, nsta, nend, iimpr1)
        call read_real(line, nsta, nend, rimpr2)

        !  Store the parameters
        !
        par%impr_atom_cls(1:4, ni) = catm(1:4)
        par%impr_force_const(ni)   = rimpr1
        par%impr_periodicity(ni)   = iimpr1
        par%impr_phase(ni)         = rimpr2

      end if

    else

      call error_msg(' Read_Par_Impr> ERROR: Format is not correct')

    end if

    return

  end subroutine read_par_impr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_par_nonb
  !> @brief        read nb information from PAR file
  !! @authors      YS
  !! @param[in]    line  : read file line
  !! @param[in]    nchar : number of characters to read
  !! @param[in]    mode  : mode (IPRE_READ or IREAD)
  !! @param[inout] par   : CHARMM PAR information
  !! @param[inout] nn    : number of read non-bonded
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_par_nonb(line, nchar, mode, par, nn)

    ! formal arguments
    character(*),            intent(in)    :: line
    integer,                 intent(in)    :: nchar
    integer,                 intent(in)    :: mode
    type(s_par),             intent(inout) :: par
    integer,                 intent(inout) :: nn

    ! local variables
    real(wp)                 :: polar, eps, rmin, polar14, eps14, rmin14
    integer                  :: nsta, nend, ndata
    character(6)             :: catm1
    character(4)             :: key4


    nsta = 1
    nend = nchar

    call read_ndata(line, nsta, nend, ndata)
    key4 = line(1:4)
    call toupper(key4)

    if (ndata == 1 .or. ndata >= 8) then

      !  Do nothing
      !
    else if (key4 .eq. 'NONB') then
      !  Do nothing
      !

    else if (ndata == 4) then

      nn = nn + 1

      if (mode == IREAD) then

        !  Read atom, polar, epsilon, Rmin/2
        !
        nsta = 1
        nend = nchar

        call read_word(line, nsta, nend, 6, catm1)
        call read_real(line, nsta, nend, polar)
        call read_real(line, nsta, nend, eps)
        call read_real(line, nsta, nend, rmin)

        !  Store the parameter
        !
        par%nonb_atom_cls(nn) = catm1
        par%nonb_polar(nn)    = polar
        par%nonb_eps(nn)      = eps
        par%nonb_rmin(nn)     = rmin
        par%nonb_polar_14(nn) = polar
        par%nonb_eps_14(nn)   = eps  
        par%nonb_rmin_14(nn)  = rmin 

      end if

    else if (ndata == 7) then

      nn = nn + 1

      if (mode == IREAD) then

        !  Read atom, ignored, epsilon, Rmin/2, ignored, eps_1-4, Rmin/2_1-4
        !
        nsta = 1
        nend = nchar

        call read_word(line, nsta, nend, 6, catm1)
        call read_real(line, nsta, nend, polar)
        call read_real(line, nsta, nend, eps)
        call read_real(line, nsta, nend, rmin)

        call read_real(line, nsta, nend, polar14)
        call read_real(line, nsta, nend, eps14)
        call read_real(line, nsta, nend, rmin14)

        !  Store the parameter
        !
        par%nonb_atom_cls(nn) = catm1
        par%nonb_polar(nn)    = polar
        par%nonb_eps(nn)      = eps
        par%nonb_rmin(nn)     = rmin
        par%nonb_polar_14(nn) = polar14
        par%nonb_eps_14(nn)   = eps14
        par%nonb_rmin_14(nn)  = rmin14

      end if

    else

      call error_msg(' Read_Par_Nonb> ERROR: Data format is not correct')

    end if

    return

  end subroutine read_par_nonb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_par_nbfi
  !> @brief        read nbfix information from PAR file
  !! @authors      NT
  !! @param[in]    line  : read file line
  !! @param[in]    nchar : number of characters to read
  !! @param[in]    mode  : mode (IPRE_READ or IREAD)
  !! @param[inout] par   : CHARMM PAR information
  !! @param[inout] nx    : number of read nbfix
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_par_nbfi(line, nchar, mode, par, nx)

    ! formal arguments
    character(*),            intent(in)    :: line
    integer,                 intent(in)    :: nchar
    integer,                 intent(in)    :: mode
    type(s_par),             intent(inout) :: par
    integer,                 intent(inout) :: nx

    ! local variables
    real(wp)                 :: v1, v2
    integer                  :: nsta, nend, ndata
    character(6)             :: a1, a2
    integer                  :: pow1,pow2

    nsta = 1
    nend = nchar

    call read_ndata(line, nsta, nend, ndata)

    if (ndata == 1) then
      
      !  Do nothing
      !
      
    else if (ndata == 4) then
      
      nx = nx + 1

      if (mode == IREAD) then

        !  Read information: atom1, atom2, eps, rmin
        !
        nsta = 1
        nend = nchar

        call read_word(line, nsta, nend, 6, a1)
        call read_word(line, nsta, nend, 6, a2)
        call read_real(line, nsta, nend, v1)
        call read_real(line, nsta, nend, v2)

        !  Store parameters
        !
        par%nbfi_atom_cls(1, nx) = a1
        par%nbfi_atom_cls(2, nx) = a2
        par%nbfi_eps     (nx)    = v1
        par%nbfi_rmin    (nx)    = v2
        par%nbfi_eps_14  (nx)    = v1
        par%nbfi_rmin_14 (nx)    = v2

      end if

    else if (ndata == 6) then

      nx = nx + 1

      if (mode == IREAD) then

        !  Read information: atom1, atom2, eps, rmin
        !
        nsta = 1
        nend = nchar

        call read_word(line, nsta, nend, 6, a1)
        call read_word(line, nsta, nend, 6, a2)
        call read_real(line, nsta, nend, v1)
        call read_real(line, nsta, nend, v2)
        call read_int(line, nsta, nend, pow1)
        call read_int(line, nsta, nend, pow2)

        !  Store parameters
        !
        par%nbfi_atom_cls(1, nx) = a1
        par%nbfi_atom_cls(2, nx) = a2
        par%nbfi_eps     (nx)    = v1
        par%nbfi_rmin    (nx)    = v2
        par%nbfi_eps_14  (nx)    = v1
        par%nbfi_rmin_14 (nx)    = v2
        par%nbfi_repuls (nx)     = pow1
        par%nbfi_attract (nx)    = pow2

      end if

    else

      call error_msg('  Read_Par_Nbfi> ERROR: Format is not correct')

    end if

    return

  end subroutine read_par_nbfi

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_par_cmap
  !> @brief        read cmap information from PAR file
  !! @authors      YS
  !! @param[inout] line  : read file line
  !! @param[in]    nchar : number of characters to read
  !! @param[in]    mode  : mode (IPRE_READ or IREAD)
  !! @param[in]    par   : CHARMM PAR information
  !! @param[inout] nc    : number of read cmaps
  !! @param[inout] mapi  : index for cmap data 
  !! @param[inout] mapi2 : index for cmap data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_par_cmap(line, nchar, mode, par, nc, mapi, mapi2)

    ! formal arguments
    character(*),            intent(inout) :: line
    integer,                 intent(in)    :: nchar
    integer,                 intent(in)    :: mode
    type(s_par),             intent(inout) :: par
    integer,                 intent(inout) :: nc
    integer,                 intent(inout) :: mapi
    integer,                 intent(inout) :: mapi2

    ! local variables
    real(wp)                 :: rval
    integer                  :: nsta, nend, ndata, i
    character(MAXCOLUMN_PAR) :: temp_line
    character(6)             :: str


    nsta = 1
    nend = nchar + 3

    temp_line(1:nchar) = line(1:nchar)
    line = ''
    line(1:nchar) = temp_line(1:nchar)

    call read_ndata(line, nsta, nend, ndata)

    if (ndata == 9) then
      
      nc = nc + 1

      if (mode == IREAD) then

        !  Read cmap atom names
        !
        nsta = 1
        nend = nchar

        do i = 1, 8
          call read_word(line, nsta, nend, 6, str)
          par%cmap_atom_cls(i,nc) = str(1:6)
        end do

        !  Read cmap resolution
        !
        call read_int(line, nsta, nend, &
             par%cmap_resolution(nc))

        mapi  = 0
        mapi2 = 1

      end if

    else if (ndata == 5) then
      
      if (mode == IREAD) then

        !  Read cmap data
        !
        nsta = 1
        nend = nchar

        do i = 1, 5
          call read_real(line, nsta, nend, rval)
          mapi = mapi + 1
          par%cmap_data(mapi,mapi2,nc) = rval
        end do


      end if

    else if (ndata == 4) then
      
      if (mode == IREAD) then

        !  Read cmap data
        !
        nsta = 1
        nend = nchar

        do i = 1, 4
          call read_real(line, nsta, nend, rval)
          mapi = mapi + 1
          par%cmap_data(mapi,mapi2,nc) = rval
        end do

        mapi  = 0
        mapi2 = mapi2 + 1

      end if

    end if

    return

  end subroutine read_par_cmap

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    skip_toppar_resi
  !> @brief        skip the toppar file RESI or PRES block
  !! @authors      NT
  !! @param[in]    file : file unit number
  !! @param[in]    mode : mode (IPRE_READ or IREAD)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine skip_toppar_resi(file, mode)

    ! formal arguments
    integer,                 intent(in)    :: file
    integer,                 intent(in)    :: mode

    ! local variables
    character(MAXCOLUMN_PAR) :: line, str
    character(6)             :: frmt


    if (mode == IPRE_READ .and. main_rank) then
      write(MsgOut,'(A)') 'Read_Par> Skip RTF section in stream file.'
      write(MsgOut,'(A)') ' '
    end if

    frmt = "      "
    if (MAXCOLUMN_PAR < 100) then
      write(frmt(1:5),'(A2,I2,A)') '(A',MAXCOLUMN_PAR,')'
    else
      write(frmt,'(A2,I3,A)') '(A',MAXCOLUMN_PAR,')'
    endif

    do while (.true.)

      read(file, frmt) line
      str = adjustl(line)
      if (str(1:3) == 'end' .or. str(1:3) == 'END') &
        exit
    
    end do

    return

  end subroutine skip_toppar_resi

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_par
  !> @brief        subroutine for writing CHARMM par file
  !! @authors      YS, NT
  !! @param[in]    file : file unit number of CHARMM parameter file
  !! @param[in]    par  : CHARMM PAR information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_par(file, par)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_par),             intent(in)    :: par

    ! local variables
    integer                  :: i, j, k
    integer                  :: nb, na, nd, ni, nn, nx, nc, mapi, mapi2


    if (par%type == ParTypeXPLOR) then
      write(MsgOut,*) 'Write_Par> write a parameter file for X-PLOR'
    else ! ParTypeCHARMM
      write(MsgOut,*) 'Write_Par> write a parameter file for CHARMM'
    end if

    !  Initialize variables
    !
    nb = par%num_bonds
    na = par%num_angles
    nd = par%num_dihedrals
    ni = par%num_impropers
    nn = par%num_atom_cls
    nx = par%num_nbfix
    nc = par%num_cmaps

    write(MsgOut,*) 'Write_Par> Number of Bonds      = ', nb
    write(MsgOut,*) 'Write_Par> Number of Angles     = ', na
    write(MsgOut,*) 'Write_Par> Number of Dihedrals  = ', nd
    write(MsgOut,*) 'Write_Par> Number of Impropers  = ', ni
    write(MsgOut,*) 'Write_Par> Number of Non-bonded = ', nn
    write(MsgOut,*) 'Write_Par> Number of Nbfix      = ', nx

    if (par%type == ParTypeCHARMM) then
      write(MSgOut,*) 'Write_Par> number of Cross-term = ', nc
    else
      write(MsgOut,*) 'Write_Par> Cross-term data are not output.'
    end if

    write(MsgOut,*) ' '

    !  Write parameters
    !
    if (par%type == ParTypeXPLOR) &
      write (file, '(A/)') 'SET ECHO=FALSE END'

    if (nb > 0) then
      if (par%type == ParTypeCHARMM) &
        write (file, '(A)') 'BONDS'
      do i = 1, nb
        call write_par_bond(file, par, i)
      end do
    end if

    if (na > 0) then
      if (par%type == ParTypeCHARMM) &
        write (file, '(/A)') 'ANGLES'
      do i = 1, na
        call write_par_angl(file, par, i)
      end do
    end if

    if (nd > 0) then
      if (par%type == ParTypeCHARMM) &
        write (file, '(/A)') 'DIHEDRALS'
      do i = 1, nd
        call write_par_dihe(file, par, i)
      end do
    end if

    if (ni > 0) then
      if (par%type == ParTypeCHARMM) &
        write (file, '(/A)') 'IMPROPER'
      do i = 1, ni
        call write_par_impr(file, par, i)
      end do
    end if

    if (nn > 0) then
      if (par%type == ParTypeCHARMM) &
        write (file, '(/A)') 'NONBONDED'
      do i = 1, nn
        call write_par_nonb(file, par, i)
      end do
    end if

    if (nx > 0) then
      if (par%type == ParTypeCHARMM) &
        write (file, '(/A)') 'NBFIX'
      do i = 1, nx
        call write_par_nbfi(file, par, i)
      end do
    end if

    if (par%type == ParTypeCHARMM .and. nc > 0) then
      write (file, '(/A)') 'CMAP'
      do i = 1, nc
        call write_par_cmap(file, par, i, 0, 1)
        mapi2 = size(par%cmap_data(1,:,nc))
        do j = 1, mapi2
          mapi = size(par%cmap_data(:,mapi2,nc))
          do k = 1, mapi, 5
            call write_par_cmap(file, par, i, k, j)
          end do
          write(file, '(A)') ''
        end do
      end do
    end if

    if (par%type == ParTypeXPLOR) then
       write (file, '(/A)') 'SET ECHO=TRUE END'
    else ! ParTypeCHARMM
       write (file, '(A)') 'END'
    end if

    return

  end subroutine write_par

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_par_bond
  !> @brief        write bond information to PAR file
  !! @authors      YS
  !! @param[in]    file : file unit number of PAR file
  !! @param[in]    par  : CHARMM PAR information
  !! @param[in]    nb   : number of bonds
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_par_bond(file, par, nb)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_par),             intent(in)    :: par
    integer,                 intent(in)    :: nb

    ! local variables
    real(wp)                 :: rbond1, rbond2
    character(6)             :: catm(2)


    ! Restore parameters
    !
    catm(1:2) = par%bond_atom_cls(1:2, nb)
    rbond1    = par%bond_force_const(nb)
    rbond2    = par%bond_dist_min(nb)

    ! Write information to line: lb(1), lb(2), kb, b0
    !
    if (par%type == ParTypeXPLOR) then
      write(file, '(3(a6,1x),f8.3,f11.4)') &
                  'BOND',                  &
                  catm(1),                 &
                  catm(2),                 &
                  rbond1,                  &
                  rbond2
    else ! ParTypeCHARMM
      write(file, '(2(a6,1x),f8.3,f11.4)') &
                  catm(1),                 &
                  catm(2),                 &
                  rbond1,                  &
                  rbond2
    end if

    return

  end subroutine write_par_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_par_angl
  !> @brief        write angle information to PAR file
  !! @authors      YS
  !! @param[in]    file : file unit number of PAR file
  !! @param[in]    par  : CHARMM PAR information
  !! @param[in]    na   : number of angles
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_par_angl(file, par, na)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_par),             intent(in)    :: par
    integer,                 intent(in)    :: na

    ! local variables
    real(wp)                 :: rangl1, rangl2, rurey1, rurey2
    character(6)             :: catm(3)


    !  Restore the parameter
    !
    catm(1:3) = par%angl_atom_cls(1:3, na)
    rangl1    = par%angl_force_const(na)
    rangl2    = par%angl_theta_min(na)
    rurey1    = par%angl_ub_force_const(na)
    rurey2    = par%angl_ub_rmin(na)

    if (rurey1 == 0.0_wp .and. rurey2 == 0.0_wp) then

      !  Read information: ia(1), ia(2), ia(3), ka, a0
      !
      if (par%type == ParTypeXPLOR) then
        write (file, '(a5,1x,3(a6,1x),f8.3,f11.4)') &
                     'ANGLE',                       &
                     catm(1),                       &
                     catm(2),                       &
                     catm(3),                       &
                     rangl1,                        &
                     rangl2
      else ! ParTypeCHARMM
        write (file, '(3(a6,1x),f8.3,f11.4)') &
                     catm(1),                 &
                     catm(2),                 &
                     catm(3),                 &
                     rangl1,                  &
                     rangl2
      end if

    else

      !  Read information: ia(1), ia(2), ia(3), ka, a0
      !
      if (par%type == ParTypeXPLOR) then
        write (file, '(a5,1x,3(a6,1x),f8.3,f11.4,f8.2,f10.5)') &
                     'ANGLE',                                  &
                     catm(1),                                  &
                     catm(2),                                  &
                     catm(3),                                  &
                     rangl1,                                   &
                     rangl2,                                   &
                     rurey1,                                   &
                     rurey2
      else ! ParTypeCHARMM
        write (file, '(3(a6,1x),f8.3,f11.4,f8.2,f10.5)') &
                     catm(1),                            &
                     catm(2),                            &
                     catm(3),                            &
                     rangl1,                             &
                     rangl2,                             &
                     rurey1,                             &
                     rurey2
      end if

    end if

    return

  end subroutine write_par_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_par_dihe
  !> @brief        write dihedral angle information to PAR file
  !! @authors      YS
  !! @param[in]    file : file unit number of PAR file
  !! @param[in]    par  : CHARMM PAR information
  !! @param[in]    nd   : number of dihedrals
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_par_dihe(file, par, nd)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_par),             intent(in)    :: par
    integer,                 intent(in)    :: nd

    ! local variables
    real(wp)                 :: rdihe1, rdihe2
    integer                  :: idihe1
    character(6)             :: catm(4)


    !  Store the parameters
    !
    catm(1:4) = par%dihe_atom_cls(1:4, nd)
    rdihe1    = par%dihe_force_const(nd)
    idihe1    = par%dihe_periodicity(nd)
    rdihe2    = par%dihe_phase(nd)

    !  Read information: id(1), id(2), id(3), id(4), kd, md, dd
    !
    if (par%type == ParTypeXPLOR) then
      write (file, '(a8,1x,4(a6,1x),f10.4,i3,f9.2)') &
                   'DIHEDRAL',                       &
                   catm(1),                          &
                   catm(2),                          &
                   catm(3),                          &
                   catm(4),                          &
                   rdihe1,                           &
                   idihe1,                           &
                   rdihe2
    else ! ParTypeCHARMM
      write (file, '(4(a6,1x),f10.4,i3,f9.2)') &
                   catm(1),                    &
                   catm(2),                    &
                   catm(3),                    &
                   catm(4),                    &
                   rdihe1,                     &
                   idihe1,                     &
                   rdihe2
    end if

    return

  end subroutine write_par_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_par_impr
  !> @brief        write improper angle information to PAR file
  !! @authors      YS
  !! @param[in]    file : file unit number of PAR file
  !! @param[in]    par  : CHARMM PAR information
  !! @param[in]    ni   : number of impropers
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_par_impr(file, par, ni)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_par),             intent(in)    :: par
    integer,                 intent(in)    :: ni

    ! local variables
    real(wp)                 :: rimpr1, rimpr2
    integer                  :: iimpr1
    character(6)             :: catm(4)


    !  Store the parameters
    !
    catm(1:4) = par%impr_atom_cls(1:4, ni)
    rimpr1    = par%impr_force_const(ni)   
    iimpr1    = par%impr_periodicity(ni)   
    rimpr2    = par%impr_phase(ni)         

    !  Read information: id(1), id(2), id(3), id(4), kd, md, dd
    !
    if (par%type == ParTypeXPLOR) then
      write (file, '(a8,1x,4(a6,1x),f10.4,i10,f12.4)') &
                   'IMPROPER',                         &
                   catm(1),                            &
                   catm(2),                            &
                   catm(3),                            &
                   catm(4),                            &
                   rimpr1,                             &
                   iimpr1,                             &
                   rimpr2
    else ! ParTypeCHARMM
      write (file, '(4(a6,1x),f10.4,i10,f12.4)') &
                   catm(1),                      &
                   catm(2),                      &
                   catm(3),                      &
                   catm(4),                      &
                   rimpr1,                       &
                   iimpr1,                       &
                   rimpr2
    end if

    return

  end subroutine write_par_impr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_par_nonb
  !> @brief        write nb information to PAR file
  !! @authors      YS
  !! @param[in]    file : file unit number of PAR file
  !! @param[in]    par  : CHARMM PAR information
  !! @param[in]    nn   : number of non-bonded
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_par_nonb(file, par, nn)
    
    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_par),             intent(in)    :: par
    integer,                 intent(in)    :: nn

    ! local variables
    real(wp)                 :: polar, eps, rmin, polar14, eps14, rmin14
    character(6)             :: catm1


    !  Store the parameter
    !
    catm1   = par%nonb_atom_cls(nn)
    polar   = par%nonb_polar(nn)    
    eps     = par%nonb_eps(nn)      
    rmin    = par%nonb_rmin(nn)     
    polar14 = par%nonb_polar_14(nn)
    eps14   = par%nonb_eps_14(nn)  
    rmin14  = par%nonb_rmin_14(nn) 

    if (par%type == ParTypeXPLOR) then

      write (file, '(a9,1x,a6,4(f11.6))') &
                   'NONBONDED',           &
                   catm1,                 &
                   eps,                   &
                   rmin,                  &
                   eps14,                 &
                   rmin14

    else ! ParTypeCHARMM

      if (polar == polar14 .and. eps ==  eps14 .and. rmin == rmin14) then

        !  Write atom, polar, epsilon, Rmin/2
        !
        write (file, '(a6,f11.6,f11.6,f13.6)') &
                     catm1,                    &
                     polar,                    &
                     eps,                      &
                     rmin

      else

        !  Write atom, ignored, epsilon, Rmin/2, ignored, eps_1-4, Rmin/2_1-4
        !
        write (file, '(a6,2(f11.6,f11.6,f13.6))') &
                     catm1,                       &
                     polar,                       &
                     eps,                         &
                     rmin,                        &
                     polar14,                     &
                     eps14,                       &
                     rmin14
      end if

    end if

    return

  end subroutine write_par_nonb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_par_nbfi
  !> @brief        write nbfix information to PAR file
  !! @authors      NT
  !! @param[in]    file : file unit number of PAR file
  !! @param[in]    par  : CHARMM PAR information
  !! @param[in]    nx   : number of nbfix
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_par_nbfi(file, par, nx)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_par),             intent(in)    :: par
    integer,                 intent(in)    :: nx

    ! local variables
    real(wp)                 :: eps, rmin, eps14, rmin14
    character(6)             :: a1, a2


    ! Restore parameters
    !
    a1     = par%nbfi_atom_cls(1, nx)
    a2     = par%nbfi_atom_cls(2, nx)
    eps    = par%nbfi_eps     (nx)
    rmin   = par%nbfi_rmin    (nx)
    eps14  = par%nbfi_eps_14  (nx)
    rmin14 = par%nbfi_rmin_14 (nx)

    if (par%type == ParTypeXPLOR) then

      ! Write information to line: atom1, atom2, eps, rmin
      !
      if (eps == eps14 .and. rmin == rmin14) then

        write(file, '(3(a6,1x),f8.3,f11.4)') &
                              'NBFIX', a1, a2, eps, rmin

      else

      ! Write information to line: atom1, atom2, eps, rmin, eps-14, rmin-14
      !
        write(file, '(3(a6,1x),f8.3,f11.4)') &
                              'NBFIX', a1, a2, eps, rmin, eps14, rmin14

      end if

    else ! ParTypeCHARMM

      ! Write information to line: atom1, atom2, eps, rmin
      !
      if (eps == eps14 .and. rmin == rmin14) then

        write(file, '(2(a6,1x),f8.3,f11.4)') a1, a2, eps, rmin

      else

      ! Write information to line: atom1, atom2, eps, rmin, eps-14, rmin-14
      !
        write(file, '(2(a6,1x),f8.3,f11.4)') a1, a2, eps, rmin, eps14, rmin14

      end if

    end if

    return

  end subroutine write_par_nbfi

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_par_cmap
  !> @brief        write cmap information to PAR file
  !! @authors      YS
  !! @param[in]    file  : file unit number of PAR file
  !! @param[in]    par   : CHARMM PAR information
  !! @param[in]    nc    : number of cmaps
  !! @param[in]    mapi  : index for cmap data
  !! @param[in]    mapi2 : index for cmap data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_par_cmap(file, par, nc, mapi, mapi2)
    
    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_par),             intent(in)    :: par
    integer,                 intent(in)    :: nc
    integer,                 intent(in)    :: mapi, mapi2

    ! local variables
    integer                  :: i, nmapi


    nmapi = size(par%cmap_data(:,mapi2,nc))

    if (mapi == 0 .and. mapi2 == 1) then

      ! Write cmap atom names and resolution
      !
      write(file, '(8(a6,1x), i4)') &
            (par%cmap_atom_cls(i,nc), i = 1, 8), &
             par%cmap_resolution(nc)

    else if (mapi + 4 <= nmapi) then

      ! Write cmap data
      !
      write(file, '(4(f10.6, 2x), f10.6)') &
            (par%cmap_data(i,mapi2,nc), i = mapi, mapi+4)

    else

      ! Write cmap data
      !
      write(file, '(3(f10.6, 2x), f10.6)') &
            (par%cmap_data(i,mapi2,nc), i = mapi, nmapi)

    end if

    return

  end subroutine write_par_cmap

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    merge_par
  !> @brief        merge two parameter informations
  !! @authors      NT
  !! @param[in]    par0 : CHARMM PAR information of the source
  !! @param[inout] par  : CHARMM PAR information of the destination
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine merge_par(par0, par)

    ! formal arguments
    type(s_par),             intent(in)    :: par0
    type(s_par),             intent(inout) :: par

    ! local variables
    type(s_par)              :: parw
    integer                  :: i, j, ii, n


    if (par%num_bonds == 0) then

      call alloc_par(par, ParBond, par0%num_bonds)
      call alloc_par(par, ParENMT, par0%num_enmt)
      call alloc_par(par, ParENMP, par0%num_enmp)
      call alloc_par(par, ParENMU, par0%num_enmu)
      call alloc_par(par, ParENMQ, par0%num_enmq)
      call alloc_par(par, ParAngl, par0%num_angles)
      call alloc_par(par, ParDihe, par0%num_dihedrals)
      call alloc_par(par, ParImpr, par0%num_impropers)
      call alloc_par(par, ParNbon, par0%num_atom_cls)
      call alloc_par(par, ParNbfi, par0%num_nbfix)
      call alloc_par(par, ParCmap, par0%num_cmaps)

      par%num_bonds     = par0%num_bonds
      par%num_enmt      = par0%num_enmt
      par%num_enmp      = par0%num_enmp
      par%num_enmu      = par0%num_enmu
      par%num_enmq      = par0%num_enmq
      par%num_angles    = par0%num_angles
      par%num_dihedrals = par0%num_dihedrals
      par%num_impropers = par0%num_impropers
      par%num_atom_cls  = par0%num_atom_cls
      par%num_nbfix     = par0%num_nbfix
      par%num_cmaps     = par0%num_cmaps

      call copy_par(par0, par)

      return

    end if

    parw%num_bonds     = par%num_bonds     + par0%num_bonds
    parw%num_enmt      = par%num_enmt      + par0%num_enmt
    parw%num_enmp      = par%num_enmp      + par0%num_enmp
    parw%num_enmu      = par%num_enmu      + par0%num_enmu
    parw%num_enmq      = par%num_enmq      + par0%num_enmq
    parw%num_angles    = par%num_angles    + par0%num_angles
    parw%num_dihedrals = par%num_dihedrals + par0%num_dihedrals
    parw%num_impropers = par%num_impropers + par0%num_impropers
    parw%num_atom_cls  = par%num_atom_cls  + par0%num_atom_cls
    parw%num_nbfix     = par%num_nbfix     + par0%num_nbfix
    parw%num_cmaps     = par%num_cmaps     + par0%num_cmaps

    call alloc_par(parw, ParBond, parw%num_bonds)
    call alloc_par(parw, ParENMT, parw%num_enmt)
    call alloc_par(parw, ParENMP, parw%num_enmp)
    call alloc_par(parw, ParENMU, parw%num_enmu)
    call alloc_par(parw, ParENMQ, parw%num_enmq)
    call alloc_par(parw, ParAngl, parw%num_angles)
    call alloc_par(parw, ParDihe, parw%num_dihedrals)
    call alloc_par(parw, ParImpr, parw%num_impropers)
    call alloc_par(parw, ParNbon, parw%num_atom_cls)
    call alloc_par(parw, ParNbfi, parw%num_nbfix)
    call alloc_par(parw, ParCmap, parw%num_cmaps)

    call copy_par(par, parw)


    ! check and copy bonds
    !
    ii = 0
    do i = 1, par0%num_bonds
      do j = 1, par%num_bonds
        if (par0%bond_atom_cls(1,i) /= par%bond_atom_cls(1,j) .or. &
            par0%bond_atom_cls(2,i) /= par%bond_atom_cls(2,j)) &
          cycle
        if (main_rank .and. vervose) &
          write(MsgOut,'(5a)') ' Merge_Par> WARNING: BOND: Multiple definition: "', &
          trim(par%bond_atom_cls(1,j)), '"-"', &
          trim(par%bond_atom_cls(2,j)), '"'
        exit
      end do
      if (j <= par%num_bonds) &
        cycle
      
      ii = ii + 1
      parw%bond_atom_cls (:,ii+par%num_bonds) = par0%bond_atom_cls (:,i)
      parw%bond_force_const(ii+par%num_bonds) = par0%bond_force_const(i)
      parw%bond_dist_min   (ii+par%num_bonds) = par0%bond_dist_min   (i)
    end do

    n = ii + par%num_bonds
    call dealloc_par(par, ParBond)
    call alloc_par  (par, ParBond, n)

    ii = 0
    do i = 1, par0%num_enmt
      do j = 1, par%num_enmt
        if (par0%bond_enm_atom_id(1,i) /= par%bond_enm_atom_id(1,j) .or. &
            par0%bond_enm_atom_id(2,i) /= par%bond_enm_atom_id(2,j)) &
          cycle
        if (main_rank) &
          write(MsgOut,'(5a)') ' Merge_Par> WARNING: ENMBOND: Multiple definition: "', &
          par%bond_enm_atom_id(1,j), '"-"', &
          par%bond_enm_atom_id(2,j), '"'
        exit
     end do
      if (j <= par%num_enmt) &
        cycle

      ii = ii + 1
      parw%bond_enm_atom_id(:,ii+par%num_enmt)= par0%bond_enm_atom_id(:,i)
      parw%bond_enm_ij_type(ii+par%num_enmt)= par0%bond_enm_ij_type(i)

    end do

    n = ii + par%num_enmt
    call dealloc_par(par, ParENMT)
    call alloc_par  (par, ParENMT, n)

    ii = 0
    do i = 1, par0%num_enmp
      do j = 1, par%num_enmp
        if (par0%bond_enm_param_type(i) /= par%bond_enm_param_type(j) ) &
          cycle
        if (main_rank) &
          write(MsgOut,'(5a)') ' Merge_Par> WARNING: ENMBOND: Multiple definition: PARAM"', &
          par%bond_enm_param_type(j),'"'
        exit
     end do
      if (j <= par%num_enmp) &
        cycle

      ii = ii + 1
      parw%bond_enm_param_type(ii+par%num_enmp)=par0%bond_enm_param_type(i)
      parw%bond_enm_force_const(ii+par%num_enmp)=par0%bond_enm_force_const(i)
      parw%bond_enm_dist_min(ii+par%num_enmp)=par0%bond_enm_dist_min(i)

    end do

    n = ii + par%num_enmp
    call dealloc_par(par, ParENMP)
    call alloc_par  (par, ParENMP, n)

    ii = 0
    do i = 1, par0%num_enmu
      do j = 1, par%num_enmu
        if (par0%angle_enm_atom_id(1,i) /= par%angle_enm_atom_id(1,j) .or. &
            par0%angle_enm_atom_id(2,i) /= par%angle_enm_atom_id(2,j) .or. &
            par0%angle_enm_atom_id(3,i) /= par%angle_enm_atom_id(3,j) ) &
          cycle
        if (main_rank) &
          write(MsgOut,'(5a)') ' Merge_Par> WARNING: ENMANGLE: Multiple definition: "', &
          par%angle_enm_atom_id(1,j), '"-"', &
          par%angle_enm_atom_id(2,j), '"-"', &
          par%angle_enm_atom_id(3,j), '"'
        exit
     end do
      if (j <= par%num_enmu) &
        cycle

      ii = ii + 1
      parw%angle_enm_atom_id(:,ii+par%num_enmu)= par0%angle_enm_atom_id(:,i)
      parw%angle_enm_ij_type(ii+par%num_enmu)= par0%angle_enm_ij_type(i)

    end do

    n = ii + par%num_enmu
    call dealloc_par(par, ParENMU)
    call alloc_par  (par, ParENMU, n)

    ii = 0
    do i = 1, par0%num_enmq
      do j = 1, par%num_enmq
        if (par0%angle_enm_param_type(i) /= par%angle_enm_param_type(j) ) &
          cycle
        if (main_rank) &
          write(MsgOut,'(5a)') ' Merge_Par> WARNING: ENMANGLE: Multiple definition: PARAM"', &
          par%angle_enm_param_type(j),'"'
        exit
     end do
      if (j <= par%num_enmq) &
        cycle

      ii = ii + 1
      parw%angle_enm_param_type(ii+par%num_enmq)=par0%angle_enm_param_type(i)
      parw%angle_enm_force_const(ii+par%num_enmq)=par0%angle_enm_force_const(i)
      parw%angle_enm_dist_min(ii+par%num_enmq)=par0%angle_enm_dist_min(i)
      parw%angle_enm_lj_eps(ii+par%num_enmq)=par0%angle_enm_lj_eps(i)
      parw%angle_enm_lj_sigma(ii+par%num_enmq)=par0%angle_enm_lj_sigma(i)

    end do

    n = ii + par%num_enmq
    call dealloc_par(par, ParENMQ)
    call alloc_par  (par, ParENMQ, n)

    ! check and copy angles
    !
    ii = 0
    do i = 1, par0%num_angles
      do j = 1, par%num_angles
        if (par0%angl_atom_cls(1,i) /= par%angl_atom_cls(1,j) .or. &
            par0%angl_atom_cls(2,i) /= par%angl_atom_cls(2,j) .or. &
            par0%angl_atom_cls(3,i) /= par%angl_atom_cls(3,j)) &
          cycle
        if (main_rank .and. vervose) &
          write(MsgOut,'(7a)') ' Merge_Par> WARNING: ANGL: Multiple definition: "', &
          trim(par%angl_atom_cls(1,j)), '"-"', &
          trim(par%angl_atom_cls(2,j)), '"-"', &
          trim(par%angl_atom_cls(3,j)), '"'
        exit
      end do
      if (j <= par%num_angles) &
        cycle
      
      ii = ii + 1
      parw%angl_atom_cls      (:,ii+par%num_angles)= par0%angl_atom_cls    (:,i)
      parw%angl_force_const   (ii+par%num_angles)  = par0%angl_force_const   (i)
      parw%angl_theta_min     (ii+par%num_angles)  = par0%angl_theta_min     (i)
      parw%angl_ub_force_const(ii+par%num_angles)  = par0%angl_ub_force_const(i)
      parw%angl_ub_rmin       (ii+par%num_angles)  = par0%angl_ub_rmin       (i)
    end do

    n = ii + par%num_angles
    call dealloc_par(par, ParAngl)
    call alloc_par  (par, ParAngl, n)


    ! check and copy dihedrals
    !
    ii = 0
    do i = 1, par0%num_dihedrals
      do j = 1, par%num_dihedrals
        if (par0%dihe_atom_cls(1,i) /= par%dihe_atom_cls(1,j) .or. &
            par0%dihe_atom_cls(2,i) /= par%dihe_atom_cls(2,j) .or. &
            par0%dihe_atom_cls(3,i) /= par%dihe_atom_cls(3,j) .or. &
            par0%dihe_atom_cls(4,i) /= par%dihe_atom_cls(4,j)) &
          cycle
        if (main_rank .and. vervose) &
          write(MsgOut,'(9a)') ' Merge_Par> WARNING: DIHE: Multiple definition: "', &
          trim(par%dihe_atom_cls(1,j)), '"-"', &
          trim(par%dihe_atom_cls(2,j)), '"-"', &
          trim(par%dihe_atom_cls(3,j)), '"-"', &
          trim(par%dihe_atom_cls(4,j)), '"'
        exit
      end do
      if (j <= par%num_dihedrals) &
        cycle
      
      ii = ii + 1
      parw%dihe_atom_cls (:,ii+par%num_dihedrals) = par0%dihe_atom_cls (:,i)
      parw%dihe_force_const(ii+par%num_dihedrals) = par0%dihe_force_const(i)
      parw%dihe_periodicity(ii+par%num_dihedrals) = par0%dihe_periodicity(i)
      parw%dihe_phase      (ii+par%num_dihedrals) = par0%dihe_phase      (i)
    end do

    n = ii + par%num_dihedrals
    call dealloc_par(par, ParDihe)
    call alloc_par  (par, ParDihe, n)


    ! check and copy impropers
    !
    ii = 0
    do i = 1, par0%num_impropers
      do j = 1, par%num_impropers
        if (par0%impr_atom_cls(1,i) /= par%impr_atom_cls(1,j) .or. &
            par0%impr_atom_cls(2,i) /= par%impr_atom_cls(2,j) .or. &
            par0%impr_atom_cls(3,i) /= par%impr_atom_cls(3,j) .or. &
            par0%impr_atom_cls(4,i) /= par%impr_atom_cls(4,j)) &
          cycle
        if (main_rank .and. vervose) &
          write(MsgOut,'(9a)') ' Merge_Par> WARNING: IMPR: Multiple definition: "', &
          trim(par%impr_atom_cls(1,j)), '"-"', &
          trim(par%impr_atom_cls(2,j)), '"-"', &
          trim(par%impr_atom_cls(3,j)), '"-"', &
          trim(par%impr_atom_cls(4,j)), '"'
        exit
      end do
      if (j <= par%num_impropers) &
        cycle
      
      ii = ii + 1
      parw%impr_atom_cls (:,ii+par%num_impropers) = par0%impr_atom_cls (:,i)
      parw%impr_force_const(ii+par%num_impropers) = par0%impr_force_const(i)
      parw%impr_periodicity(ii+par%num_impropers) = par0%impr_periodicity(i)
      parw%impr_phase      (ii+par%num_impropers) = par0%impr_phase      (i)
    end do

    n = ii + par%num_impropers
    call dealloc_par(par, ParImpr)
    call alloc_par  (par, ParImpr, n)


    ! check and copy nonbonds
    !
    ii = 0
    do i = 1, par0%num_atom_cls
      do j = 1, par%num_atom_cls
        if (par0%nonb_atom_cls(i) /= par%nonb_atom_cls(j)) &
          cycle
        if (main_rank .and. vervose) &
          write(MsgOut,'(3a)') ' Merge_Par> WARNING: NONB: Multiple definition: "', &
          trim(par%nonb_atom_cls(j)), '"'
        exit
      end do
      if (j <= par%num_atom_cls) &
        cycle
      
      ii = ii + 1
      parw%nonb_atom_cls(ii+par%num_atom_cls) = par0%nonb_atom_cls(i)
      parw%nonb_polar   (ii+par%num_atom_cls) = par0%nonb_polar   (i)
      parw%nonb_eps     (ii+par%num_atom_cls) = par0%nonb_eps     (i)
      parw%nonb_rmin    (ii+par%num_atom_cls) = par0%nonb_rmin    (i)
      parw%nonb_polar_14(ii+par%num_atom_cls) = par0%nonb_polar_14(i)
      parw%nonb_eps_14  (ii+par%num_atom_cls) = par0%nonb_eps_14  (i)
      parw%nonb_rmin_14 (ii+par%num_atom_cls) = par0%nonb_rmin_14 (i)
      parw%nonb_lj_type_repuls(ii+par%num_atom_cls) = &
                                        par0%nonb_lj_type_repuls  (i)
      parw%nonb_lj_type_attract (ii+par%num_atom_cls) = &
                                        par0%nonb_lj_type_attract (i)
    end do

    n = ii + par%num_atom_cls
    call dealloc_par(par, ParNbon)
    call alloc_par  (par, ParNbon, n)


    ! check and copy nbfix
    !
    ii = 0
    do i = 1, par0%num_nbfix
      do j = 1, par%num_nbfix
        if (par0%nbfi_atom_cls(1,i) /= par%nbfi_atom_cls(1,j) .or. &
            par0%nbfi_atom_cls(2,i) /= par%nbfi_atom_cls(2,j)) &
          cycle
        if (main_rank .and. vervose) &
          write(MsgOut,'(5a)') ' Merge_Par> WARNING: NBFIX: Multiple definition: "',&
          trim(par%nbfi_atom_cls(1,j)), '"-"', &
          trim(par%nbfi_atom_cls(2,j)), '"'
        exit
      end do
      if (j <= par%num_nbfix) &
        cycle
      
      ii = ii + 1
      parw%nbfi_atom_cls(:,ii+par%num_nbfix) = par0%nbfi_atom_cls(:,i)
      parw%nbfi_eps       (ii+par%num_nbfix) = par0%nbfi_eps       (i)
      parw%nbfi_rmin      (ii+par%num_nbfix) = par0%nbfi_rmin      (i)
      parw%nbfi_eps_14    (ii+par%num_nbfix) = par0%nbfi_eps_14    (i)
      parw%nbfi_rmin_14   (ii+par%num_nbfix) = par0%nbfi_rmin_14   (i)
      parw%nbfi_repuls    (ii+par%num_nbfix) = par0%nbfi_repuls    (i)
      parw%nbfi_attract   (ii+par%num_nbfix) = par0%nbfi_attract   (i)
    end do

    n = ii + par%num_nbfix
    call dealloc_par(par, ParNbfi)
    call alloc_par  (par, ParNbfi, n)


    ! check and copy cmaps
    !
    ii = 0
    do i = 1, par0%num_cmaps
      do j = 1, par%num_cmaps
        if (par0%cmap_atom_cls(1,i) /= par%cmap_atom_cls(1,j) .or. &
            par0%cmap_atom_cls(2,i) /= par%cmap_atom_cls(2,j) .or. &
            par0%cmap_atom_cls(3,i) /= par%cmap_atom_cls(3,j) .or. &
            par0%cmap_atom_cls(4,i) /= par%cmap_atom_cls(4,j) .or. &
            par0%cmap_atom_cls(5,i) /= par%cmap_atom_cls(5,j) .or. &
            par0%cmap_atom_cls(6,i) /= par%cmap_atom_cls(6,j) .or. &
            par0%cmap_atom_cls(7,i) /= par%cmap_atom_cls(7,j) .or. &
            par0%cmap_atom_cls(8,i) /= par%cmap_atom_cls(8,j)) &
          cycle
        if (main_rank .and. vervose) &
          write(MsgOut,'(16a)')' Merge_Par> WARNING: CMAP: Multiple definition: "', &
          trim(par%cmap_atom_cls(1,j)), '","', &
          trim(par%cmap_atom_cls(2,j)), '","', &
          trim(par%cmap_atom_cls(3,j)), '","', &
          trim(par%cmap_atom_cls(4,j)), '","', &
          trim(par%cmap_atom_cls(5,j)), '","', &
          trim(par%cmap_atom_cls(6,j)), '","', &
          trim(par%cmap_atom_cls(7,j)), '","', &
          trim(par%cmap_atom_cls(8,j)), '"'
        exit
      end do
      if (j <= par%num_cmaps) &
        cycle
      
      ii = ii + 1
      parw%cmap_atom_cls(:,ii+par%num_cmaps) = par0%cmap_atom_cls(:,i)
      parw%cmap_resolution(ii+par%num_cmaps) = par0%cmap_resolution(i)
      parw%cmap_data  (:,:,ii+par%num_cmaps) = par0%cmap_data  (:,:,i)
    end do

    n = ii + par%num_cmaps
    call dealloc_par(par, ParCmap)
    call alloc_par  (par, ParCmap, n)

    par%num_bonds     = size(par%bond_atom_cls(1,:))
    par%num_enmt      = size(par%bond_enm_ij_type(:))
    par%num_enmp      = size(par%bond_enm_param_type(:))
    par%num_enmu      = size(par%angle_enm_ij_type(:))
    par%num_enmq      = size(par%angle_enm_param_type(:))
    par%num_angles    = size(par%angl_atom_cls(1,:))
    par%num_dihedrals = size(par%dihe_atom_cls(1,:))
    par%num_impropers = size(par%impr_atom_cls(1,:))
    par%num_atom_cls  = size(par%nonb_atom_cls(:))
    par%num_nbfix     = size(par%nbfi_atom_cls(1,:))
    par%num_cmaps     = size(par%cmap_atom_cls(1,:))

    call copy_par(parw, par)

    call dealloc_par_all(parw)

    return

  end subroutine merge_par

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    copy_par
  !> @brief        copy parameter informations
  !! @authors      NT
  !! @param[in]    par_src : CHARMM PAR information of the source
  !! @param[inout] par_dst : CHARMM PAR information of the destination
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine copy_par(par_src, par_dst)

    ! formal arguments
    type(s_par),             intent(in)    :: par_src
    type(s_par),             intent(inout) :: par_dst

    ! local variables
    integer                  :: n


    n = min(par_dst%num_bonds, par_src%num_bonds)
    par_dst%bond_atom_cls (:,1:n) = par_src%bond_atom_cls (:,1:n)
    par_dst%bond_force_const(1:n) = par_src%bond_force_const(1:n)
    par_dst%bond_dist_min   (1:n) = par_src%bond_dist_min   (1:n)

    n = min(par_dst%num_enmt, par_src%num_enmt)
    par_dst%bond_enm_atom_id(:,1:n) =    par_src%bond_enm_atom_id(:,1:n)
    par_dst%bond_enm_ij_type(1:n)   =    par_src%bond_enm_ij_type(1:n)

    n = min(par_dst%num_enmp, par_src%num_enmp)
    par_dst%bond_enm_param_type(1:n)=    par_src%bond_enm_param_type(1:n)
    par_dst%bond_enm_force_const(1:n)=    par_src%bond_enm_force_const(1:n)
    par_dst%bond_enm_dist_min(1:n)  =    par_src%bond_enm_dist_min(1:n)

    n = min(par_dst%num_enmu, par_src%num_enmu)
    par_dst%angle_enm_atom_id(:,1:n) =    par_src%angle_enm_atom_id(:,1:n)
    par_dst%angle_enm_ij_type(1:n)   =    par_src%angle_enm_ij_type(1:n)

    n = min(par_dst%num_enmq, par_src%num_enmq)
    par_dst%angle_enm_param_type(1:n)=    par_src%angle_enm_param_type(1:n)
    par_dst%angle_enm_force_const(1:n)=    par_src%angle_enm_force_const(1:n)
    par_dst%angle_enm_dist_min(1:n)  =    par_src%angle_enm_dist_min(1:n)
    par_dst%angle_enm_lj_eps(1:n)=    par_src%angle_enm_lj_eps(1:n)
    par_dst%angle_enm_lj_sigma(1:n)  =    par_src%angle_enm_lj_sigma(1:n)

    n = min(par_dst%num_angles, par_src%num_angles)
    par_dst%angl_atom_cls    (:,1:n) = par_src%angl_atom_cls    (:,1:n)
    par_dst%angl_force_const   (1:n) = par_src%angl_force_const   (1:n)
    par_dst%angl_theta_min     (1:n) = par_src%angl_theta_min     (1:n)
    par_dst%angl_ub_force_const(1:n) = par_src%angl_ub_force_const(1:n)
    par_dst%angl_ub_rmin       (1:n) = par_src%angl_ub_rmin       (1:n)

    n = min(par_dst%num_dihedrals, par_src%num_dihedrals)
    par_dst%dihe_atom_cls (:,1:n) = par_src%dihe_atom_cls (:,1:n)
    par_dst%dihe_force_const(1:n) = par_src%dihe_force_const(1:n)
    par_dst%dihe_periodicity(1:n) = par_src%dihe_periodicity(1:n)
    par_dst%dihe_phase      (1:n) = par_src%dihe_phase      (1:n)

    n = min(par_dst%num_impropers, par_src%num_impropers)
    par_dst%impr_atom_cls (:,1:n) = par_src%impr_atom_cls (:,1:n)
    par_dst%impr_force_const(1:n) = par_src%impr_force_const(1:n)
    par_dst%impr_periodicity(1:n) = par_src%impr_periodicity(1:n)
    par_dst%impr_phase      (1:n) = par_src%impr_phase      (1:n)

    n = min(par_dst%num_atom_cls, par_src%num_atom_cls)
    par_dst%nonb_atom_cls(1:n) = par_src%nonb_atom_cls(1:n)
    par_dst%nonb_polar   (1:n) = par_src%nonb_polar   (1:n)
    par_dst%nonb_eps     (1:n) = par_src%nonb_eps     (1:n)
    par_dst%nonb_rmin    (1:n) = par_src%nonb_rmin    (1:n)
    par_dst%nonb_polar_14(1:n) = par_src%nonb_polar_14(1:n)
    par_dst%nonb_eps_14  (1:n) = par_src%nonb_eps_14  (1:n)
    par_dst%nonb_rmin_14 (1:n) = par_src%nonb_rmin_14 (1:n)
    par_dst%nonb_lj_type_repuls (1:n) = par_src%nonb_lj_type_repuls (1:n)
    par_dst%nonb_lj_type_attract(1:n) = par_src%nonb_lj_type_attract(1:n)

    n = min(par_dst%num_nbfix, par_src%num_nbfix)
    par_dst%nbfi_atom_cls(:,1:n) = par_src%nbfi_atom_cls(:,1:n)
    par_dst%nbfi_eps       (1:n) = par_src%nbfi_eps       (1:n)
    par_dst%nbfi_rmin      (1:n) = par_src%nbfi_rmin      (1:n)
    par_dst%nbfi_eps_14    (1:n) = par_src%nbfi_eps_14    (1:n)
    par_dst%nbfi_rmin_14   (1:n) = par_src%nbfi_rmin_14   (1:n)
    par_dst%nbfi_repuls    (1:n) = par_src%nbfi_repuls    (1:n)
    par_dst%nbfi_attract   (1:n) = par_src%nbfi_attract   (1:n)

    n = min(par_dst%num_cmaps, par_src%num_cmaps)
    par_dst%cmap_atom_cls(:,1:n) = par_src%cmap_atom_cls(:,1:n)
    par_dst%cmap_resolution(1:n) = par_src%cmap_resolution(1:n)
    par_dst%cmap_data  (:,:,1:n) = par_src%cmap_data  (:,:,1:n)

    return

  end subroutine copy_par

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    resolve_wildcard
  !> @brief        resolve CHARMM19 wild-card nonbonded atom name
  !! @authors      NT
  !! @param[inout] par : CHARMM PAR information
  !! @param[in]    top : CHARMM TOP information (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine resolve_wildcard(par, top)

    ! parameters
    character(1),  parameter :: WildCardN = '*'
    character(1),  parameter :: WildCard1 = '%'

    ! formal arguments
    type(s_par),             intent(inout) :: par
    type(s_top),   optional, intent(in)    :: top

    ! local variables
    type(s_par)              :: par0
    integer                  :: i, j, k, m, i0, iw
    character(6)             :: catm, ctatm


    ! pre-check
    do i = 1, par%num_atom_cls
      if (scan(par%nonb_atom_cls(i), WildCardN//WildCard1) /= 0) &
        exit
    end do
    if (i > par%num_atom_cls) &
      return

    if (.not. present(top)) then
      if (main_rank .and. vervose) &
        write(MsgOut,'(a)') &
        'Resolve_Wildcard> WARNING: WildCard must be resolved.'
      return
    end if

    ! resolve wild-card
    if (main_rank .and. vervose) &
      write(MsgOut,'(a)') 'Resolve_Wildcard> '

    call alloc_par(par0, ParNbon, top%num_atom_cls)

    i0 = 0
    do i = 1, par%num_atom_cls

      catm = par%nonb_atom_cls(i)

      iw = scan(catm, WildCardN//WildCard1)

      if (iw == 0) then

        i0 = i0 + 1
        par0%nonb_atom_cls(i0) = catm
        par0%nonb_polar(i0)    = par%nonb_polar(i)   
        par0%nonb_eps(i0)      = par%nonb_eps(i)     
        par0%nonb_rmin(i0)     = par%nonb_rmin(i)    
        par0%nonb_polar_14(i0) = par%nonb_polar_14(i)
        par0%nonb_eps_14(i0)   = par%nonb_eps_14(i)  
        par0%nonb_rmin_14(i0)  = par%nonb_rmin_14(i) 
        par0%nonb_lj_type_repuls(i0)   = par%nonb_lj_type_repuls(i)
        par0%nonb_lj_type_attract(i0)  = par%nonb_lj_type_attract(i)

      else

        do j = 1, top%num_atom_cls

          ctatm = top%atom_cls_name(j)

          if (catm(iw:iw)  == WildCardN     .and. &
              catm(1:iw-1) == ctatm(1:iw-1) .or.  &
              catm(iw:iw)  == WildCard1     .and. &
              catm(1:iw-1) == ctatm(1:iw-1) .and. &
              len_trim(ctatm) == iw) then

            do k = 1, par%num_atom_cls
              if (ctatm == par%nonb_atom_cls(k)) &
                exit
            end do

            do m = 1, i0
              if (ctatm == par0%nonb_atom_cls(m)) &
                exit
            end do

            if (k > par%num_atom_cls .and. m > i0) then
              i0 = i0 + 1
              par0%nonb_atom_cls(i0) = ctatm
              par0%nonb_polar(i0)    = par%nonb_polar(i)   
              par0%nonb_eps(i0)      = par%nonb_eps(i)     
              par0%nonb_rmin(i0)     = par%nonb_rmin(i)    
              par0%nonb_polar_14(i0) = par%nonb_polar_14(i)
              par0%nonb_eps_14(i0)   = par%nonb_eps_14(i)  
              par0%nonb_rmin_14(i0)  = par%nonb_rmin_14(i) 
              par0%nonb_lj_type_repuls(i0)   = par%nonb_lj_type_repuls(i)
              par0%nonb_lj_type_attract(i0)  = par%nonb_lj_type_attract(i)

              if (main_rank) &
                write(MsgOut,'(4a)') '    ', catm, ' => ', ctatm

            end if

          end if

        end do
      end if

    end do
    
    call dealloc_par(par, ParNbon)
    call alloc_par(par, ParNbon, i0)

    par%nonb_atom_cls(1:i0) = par0%nonb_atom_cls(1:i0)
    par%nonb_polar(1:i0)    = par0%nonb_polar(1:i0)
    par%nonb_eps(1:i0)      = par0%nonb_eps(1:i0)
    par%nonb_rmin(1:i0)     = par0%nonb_rmin(1:i0)
    par%nonb_polar_14(1:i0) = par0%nonb_polar_14(1:i0)
    par%nonb_eps_14(1:i0)   = par0%nonb_eps_14(1:i0)
    par%nonb_rmin_14(1:i0)  = par0%nonb_rmin_14(1:i0)
    par%num_atom_cls = i0
    par%nonb_lj_type_repuls(1:i0)   = par0%nonb_lj_type_repuls(1:i0)
    par%nonb_lj_type_attract(1:i0)  = par0%nonb_lj_type_attract(1:i0)

    call dealloc_par(par0, ParNbon)

    if (main_rank) &
      write(MsgOut,'(a)') ' '

    return

  end subroutine resolve_wildcard

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_par_enmt
  !> @brief        read elastic network bond type information from PAR file
  !! @authors      RU
  !! @param[in]    line  : read file line
  !! @param[in]    nchar : number of characters to read
  !! @param[in]    mode  : mode (IPRE_READ or IREAD)
  !! @param[inout] par   : SPICA PAR information
  !! @param[inout] nb    :
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_par_enmt(line, nchar, mode, par, nb)

    ! formal arguments
    character(*),            intent(in)    :: line
    integer,                 intent(in)    :: nchar
    integer,                 intent(in)    :: mode
    type(s_par),             intent(inout) :: par
    integer,                 intent(inout) :: nb

    ! local variables
    real(wp)                 :: rbond1, rbond2
    integer                  :: nsta, nend, ndata
    character(6)             :: catm(2)
    integer                  :: atomi,atomj,enmtype


    nsta = 1
    nend = nchar

    call read_ndata(line, nsta, nend, ndata)

    if (ndata == 1) then

      !  Do nothing
      !

    else if (ndata == 3) then

      nb = nb + 1

      if (mode == IREAD) then

        !  Read information: lb(1), lb2(2), parmID
        !
        nsta = 1
        nend = nchar

        call read_int(line, nsta, nend,  atomi)
        call read_int(line, nsta, nend,  atomj)
        call read_int(line, nsta, nend, enmtype)
        ! call read_real(line, nsta, nend, rbond2)

        !  Store parameters
        !
        par%bond_enm_atom_id(1, nb) =atomi
        par%bond_enm_atom_id(2, nb) =atomj
        par%bond_enm_ij_type(nb)   = enmtype

      end if

    else

      call error_msg('  Read_Par_ENMT_Bond> ERROR: Format is not correct')

    end if

    return

  end subroutine read_par_enmt

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_par_enmp
  !> @brief        read elastic network bond param information from PAR file
  !! @authors      RU
  !! @param[in]    line  : read file line
  !! @param[in]    nchar : number of characters to read
  !! @param[in]    mode  : mode (IPRE_READ or IREAD)
  !! @param[inout] par   : CHARMM PAR information
  !! @param[inout] nb    :
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_par_enmp(line, nchar, mode, par, nb)

    ! formal arguments
    character(*),            intent(in)    :: line
    integer,                 intent(in)    :: nchar
    integer,                 intent(in)    :: mode
    type(s_par),             intent(inout) :: par
    integer,                 intent(inout) :: nb

    ! local variables
    real(wp)                 :: rbond1, rbond2
    integer                  :: nsta, nend, ndata
    character(6)             :: catm(2)
    integer                  :: enmtype


    nsta = 1
    nend = nchar

    call read_ndata(line, nsta, nend, ndata)

    if (ndata == 1) then

      !  Do nothing
      !

    else if (ndata == 3) then

      nb = nb + 1

      if (mode == IREAD) then

        !  Read information: lb(1), lb2(2), kb, b0
        !
        nsta = 1
        nend = nchar

        call read_int (line, nsta, nend,  enmtype)
        call read_real(line, nsta, nend, rbond1)
        call read_real(line, nsta, nend, rbond2)

        !  Store parameters
        !
        par%bond_enm_param_type(nb)  = enmtype
        par%bond_enm_force_const(nb) = rbond1
        par%bond_enm_dist_min(nb)    = rbond2

      end if

    else

      call error_msg('  Read_Par_ENMP_Bond> ERROR: Format is not correct')

    end if

    return

  end subroutine read_par_enmp

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_par_enmu
  !> @brief        read elastic network bond type information from PAR file
  !! @authors      RU
  !! @param[in]    line  : read file line
  !! @param[in]    nchar : number of characters to read
  !! @param[in]    mode  : mode (IPRE_READ or IREAD)
  !! @param[inout] par   : SPICA PAR information
  !! @param[inout] nb    :
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_par_enmu(line, nchar, mode, par, nb)

    ! formal arguments
    character(*),            intent(in)    :: line
    integer,                 intent(in)    :: nchar
    integer,                 intent(in)    :: mode
    type(s_par),             intent(inout) :: par
    integer,                 intent(inout) :: nb

    ! local variables
    real(wp)                 :: rbond1, rbond2
    integer                  :: nsta, nend, ndata
    character(6)             :: catm(2)
    integer                  :: atomi,atomj,atomk,enmtype


    nsta = 1
    nend = nchar

    call read_ndata(line, nsta, nend, ndata)

    if (ndata == 1) then

      !  Do nothing
      !

    else if (ndata == 4) then

      nb = nb + 1

      if (mode == IREAD) then

        !  Read information: lb(1), lb2(2),lb(3), paramID
        !
        nsta = 1
        nend = nchar

        call read_int(line, nsta, nend,  atomi)
        call read_int(line, nsta, nend,  atomj)
        call read_int(line, nsta, nend,  atomk)
        call read_int(line, nsta, nend, enmtype)

        !  Store parameters
        !
        par%angle_enm_atom_id(1, nb) =atomi
        par%angle_enm_atom_id(2, nb) =atomj
        par%angle_enm_atom_id(3, nb) =atomk
        par%angle_enm_ij_type(nb)   = enmtype

      end if

    else

      call error_msg('  Read_Par_ENMU_Angle> ERROR: Format is not correct')

    end if

    return

  end subroutine read_par_enmu

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_par_enmq
  !> @brief        read elastic network bond param information from PAR file
  !! @authors      RU
  !! @param[in]    line  : read file line
  !! @param[in]    nchar : number of characters to read
  !! @param[in]    mode  : mode (IPRE_READ or IREAD)
  !! @param[inout] par   : SPICA PAR information
  !! @param[inout] nb    :
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_par_enmq(line, nchar, mode, par, nb)

    ! formal arguments
    character(*),            intent(in)    :: line
    integer,                 intent(in)    :: nchar
    integer,                 intent(in)    :: mode
    type(s_par),             intent(inout) :: par
    integer,                 intent(inout) :: nb

    ! local variables
    real(wp)                 :: rbond1, rbond2,rlj1,rlj2
    integer                  :: nsta, nend, ndata
    character(6)             :: catm(2)
    integer                  :: enmtype


    nsta = 1
    nend = nchar

    call read_ndata(line, nsta, nend, ndata)

    if (ndata == 1) then

      !  Do nothing
      !

    else if (ndata == 5) then

      nb = nb + 1

      if (mode == IREAD) then

        !  Read information: lb(1), lb2(2), kb, b0
        !
        nsta = 1
        nend = nchar

        call read_int (line, nsta, nend,  enmtype)
        call read_real(line, nsta, nend, rbond1)
        call read_real(line, nsta, nend, rbond2)
        call read_real(line, nsta, nend, rlj1)
        call read_real(line, nsta, nend, rlj2)

        !  Store parameters
        !
        par%angle_enm_param_type(nb)  = enmtype
        par%angle_enm_force_const(nb) = rbond1
        par%angle_enm_dist_min(nb)    = rbond2
        par%angle_enm_lj_eps(nb)      = rlj1
        par%angle_enm_lj_sigma(nb)    = rlj2

      end if

    else

      call error_msg('  Read_Par_ENMQ_Angle> ERROR: Format is not correct')

    end if

    return

  end subroutine read_par_enmq

end module fileio_par_mod
