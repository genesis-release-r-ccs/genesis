!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pr_cached_molecule_mod
!> @brief   cached molecule information
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module pr_cached_molecule_mod

  use select_atoms_str_mod
  use fileio_top_mod
  use fileio_par_mod
  use fileio_psf_mod
  use fileio_mod
  use string_mod
  use messages_mod
  use constants_mod
  use omp_lib

  implicit none
  private

  ! public
  !

  ! parameters
  integer,        public, parameter :: CMMaxSelAtomsList = 10

  ! variables
  character(MaxFilename), public :: cm_cache_dir  = './'
  character(MaxFilename), public :: cm_cache_name = ''

  ! s_molecule
  integer,        public :: cm_num_atoms           = 0
  integer,        public :: cm_num_bonds           = 0
  integer,        public :: cm_num_angles          = 0
  integer,        public :: cm_num_dihedrals       = 0
  integer,        public :: cm_num_impropers       = 0
  integer,        public :: cm_num_cmaps           = 0
  integer,        public :: cm_num_molecule_no     = 0
  ! s_enefunc%table
  integer,        public :: cm_num_water_list      = 0
  integer,        public :: cm_num_solute_list     = 0
  integer,        public :: cm_num_solute_list_inv = 0
  integer,        public :: cm_num_solute          = 0
  ! s_constraints
  integer,        public :: cm_num_h_bond          = 0
  ! s_selatoms
  integer,        public :: cm_num_selatoms_list   = 0
  integer,        public :: cm_num_selatoms(1:CMMaxSelAtomsList) = 0

  ! subroutines
  public :: cm_initialize
  public :: cm_finalize

  ! s_molecule
  public :: cm_create_molecules
  public :: cm_delete_molecules
  public :: cm_regist_molecules
  public :: cm_atom_no
  public :: cm_atom_cls_no
  public :: cm_atom_cls_name
  public :: cm_atom_name
  public :: cm_atom_coord
  public :: cm_atom_refcoord
  public :: cm_atom_velocity
  public :: cm_charge
  public :: cm_mass
  public :: cm_residue_no
  public :: cm_residue_name
  public :: cm_segment_name
  public :: cm_molecule_no
  public :: cm_bond_list
  public :: cm_angl_list
  public :: cm_dihe_list
  public :: cm_impr_list
  public :: cm_cmap_list
  public :: cm_set_atom_coord
  public :: cm_set_atom_refcoord
  public :: cm_set_atom_velocity
  public :: cm_set_molecule_no
  public :: cm_sync_atom_coord
  public :: cm_sync_atom_refcoord
  public :: cm_sync_atom_velocity
  public :: cm_sync_molecule_no
  public :: cm_backup_atom_coord
  public :: cm_backup_atom_velocity
  public :: cm_recall_atom_coord
  public :: cm_recall_atom_velocity
  
  ! s_enefunc%table
  public :: cm_create_water_list
  public :: cm_create_solute_list
  public :: cm_create_solute_list_inv
  public :: cm_create_solute
  public :: cm_delete_water_list
  public :: cm_delete_solute_list
  public :: cm_delete_solute_list_inv
  public :: cm_delete_solute
  public :: cm_water_list
  public :: cm_solute_list
  public :: cm_solute_list_inv
  public :: cm_solute
  public :: cm_set_water_list
  public :: cm_set_solute_list
  public :: cm_set_solute_list_inv
  public :: cm_set_solute
  public :: cm_sync_water_list
  public :: cm_sync_solute_list
  public :: cm_sync_solute_list_inv
  public :: cm_sync_solute

  ! s_constraints
  public :: cm_create_duplicate
  public :: cm_create_H_index
  public :: cm_delete_duplicate
  public :: cm_delete_H_index
  public :: cm_duplicate
  public :: cm_H_index
  public :: cm_set_duplicate
  public :: cm_set_H_index
  public :: cm_sync_duplicate
  public :: cm_sync_H_index

  ! s_selatoms
  public :: cm_create_selatoms
  public :: cm_delete_selatoms
  public :: cm_selatoms


  ! private
  !

  ! parameters
  integer, public, parameter :: CMNumPages = 1000
  integer, public, parameter :: CMPageSize = 10000
  integer, public, parameter :: CMAnglSize = 10002

  ! structures
  type s_cm_page
    integer          :: mem_no
    integer          :: offset
    integer          :: size
    logical          :: dirty
    type(s_cm_page), pointer :: next
  end type s_cm_page

  type s_cm_selatoms
    integer,         allocatable :: v(:,:)
  end type s_cm_selatoms

  ! variables
  type(s_cm_page)    :: cmp_atom_no        (CMNumPages+1)
  type(s_cm_page)    :: cmp_atom_cls_no    (CMNumPages+1)
  type(s_cm_page)    :: cmp_atom_cls_name  (CMNumPages+1)
  type(s_cm_page)    :: cmp_atom_name      (CMNumPages+1)
  type(s_cm_page)    :: cmp_atom_coord     (CMNumPages+1)
  type(s_cm_page)    :: cmp_atom_velocity  (CMNumPages+1)
  type(s_cm_page)    :: cmp_atom_refcoord  (CMNumPages+1)
  type(s_cm_page)    :: cmp_charge         (CMNumPages+1)
  type(s_cm_page)    :: cmp_mass           (CMNumPages+1)
  type(s_cm_page)    :: cmp_residue_no     (CMNumPages+1)
  type(s_cm_page)    :: cmp_residue_name   (CMNumPages+1)
  type(s_cm_page)    :: cmp_segment_name   (CMNumPages+1)
  type(s_cm_page)    :: cmp_molecule_no    (CMNumPages+1)
  type(s_cm_page)    :: cmp_bond_list      (CMNumPages+1)
  type(s_cm_page)    :: cmp_angl_list      (CMNumPages+1)
  type(s_cm_page)    :: cmp_dihe_list      (CMNumPages+1)
  type(s_cm_page)    :: cmp_impr_list      (CMNumPages+1)
  type(s_cm_page)    :: cmp_cmap_list      (CMNumPages+1)
  type(s_cm_page)    :: cmp_water_list     (CMNumPages+1)
  type(s_cm_page)    :: cmp_solute_list    (CMNumPages+1)
  type(s_cm_page)    :: cmp_solute_list_inv(CMNumPages+1)
  type(s_cm_page)    :: cmp_solute         (CMNumPages+1)
  type(s_cm_page)    :: cmp_duplicate      (CMNumPages+1)
  type(s_cm_page)    :: cmp_H_index        (CMNumPages+1)
  type(s_cm_page)    :: cmp_selatoms       (CMNumPages+1, CMMaxSelAtomsList)

  !$omp threadprivate ( cmp_atom_no )
  !$omp threadprivate ( cmp_atom_cls_no )
  !$omp threadprivate ( cmp_atom_cls_name )
  !$omp threadprivate ( cmp_atom_name )
  !$omp threadprivate ( cmp_atom_coord )
  !$omp threadprivate ( cmp_atom_velocity )
  !$omp threadprivate ( cmp_atom_refcoord )
  !$omp threadprivate ( cmp_charge )
  !$omp threadprivate ( cmp_mass )
  !$omp threadprivate ( cmp_residue_no )
  !$omp threadprivate ( cmp_residue_name )
  !$omp threadprivate ( cmp_segment_name )
  !$omp threadprivate ( cmp_molecule_no )
  !$omp threadprivate ( cmp_bond_list )
  !$omp threadprivate ( cmp_angl_list )
  !$omp threadprivate ( cmp_dihe_list )
  !$omp threadprivate ( cmp_impr_list )
  !$omp threadprivate ( cmp_cmap_list )
  !$omp threadprivate ( cmp_water_list )
  !$omp threadprivate ( cmp_solute_list )
  !$omp threadprivate ( cmp_solute_list_inv )
  !$omp threadprivate ( cmp_solute )
  !$omp threadprivate ( cmp_duplicate )
  !$omp threadprivate ( cmp_H_index )
  !$omp threadprivate ( cmp_selatoms )

  integer,            target, allocatable :: cmm_atom_no        (:,:)
  integer,            target, allocatable :: cmm_atom_cls_no    (:,:)
  character(6),       target, allocatable :: cmm_atom_cls_name  (:,:)
  character(4),       target, allocatable :: cmm_atom_name      (:,:)
  real(wp),           target, allocatable :: cmm_atom_coord     (:,:,:)
  real(wp),           target, allocatable :: cmm_atom_velocity  (:,:,:)
  real(wp),           target, allocatable :: cmm_atom_refcoord  (:,:,:)
  real(wp),           target, allocatable :: cmm_charge         (:,:)
  real(wp),           target, allocatable :: cmm_mass           (:,:)
  integer,            target, allocatable :: cmm_residue_no     (:,:)
  character(6),       target, allocatable :: cmm_residue_name   (:,:)
  character(4),       target, allocatable :: cmm_segment_name   (:,:)
  integer,            target, allocatable :: cmm_molecule_no    (:,:)
  integer,            target, allocatable :: cmm_bond_list      (:,:,:)
  integer,            target, allocatable :: cmm_angl_list      (:,:,:)
  integer,            target, allocatable :: cmm_dihe_list      (:,:,:)
  integer,            target, allocatable :: cmm_impr_list      (:,:,:)
  integer,            target, allocatable :: cmm_cmap_list      (:,:,:)
  integer,            target, allocatable :: cmm_water_list     (:,:,:)
  integer,            target, allocatable :: cmm_solute_list    (:,:)
  integer,            target, allocatable :: cmm_solute_list_inv(:,:)
  logical,            target, allocatable :: cmm_solute         (:,:)
  integer,            target, allocatable :: cmm_duplicate      (:,:)
  integer,            target, allocatable :: cmm_H_index        (:,:,:)
  type(s_cm_selatoms),target              :: cmm_selatoms(CMMaxSelAtomsList)

  !$omp threadprivate ( cmm_atom_no )
  !$omp threadprivate ( cmm_atom_cls_no )
  !$omp threadprivate ( cmm_atom_cls_name )
  !$omp threadprivate ( cmm_atom_name )
  !$omp threadprivate ( cmm_atom_coord )
  !$omp threadprivate ( cmm_atom_velocity )
  !$omp threadprivate ( cmm_atom_refcoord )
  !$omp threadprivate ( cmm_charge )
  !$omp threadprivate ( cmm_mass )
  !$omp threadprivate ( cmm_residue_no )
  !$omp threadprivate ( cmm_residue_name )
  !$omp threadprivate ( cmm_segment_name )
  !$omp threadprivate ( cmm_molecule_no )
  !$omp threadprivate ( cmm_bond_list )
  !$omp threadprivate ( cmm_angl_list )
  !$omp threadprivate ( cmm_dihe_list )
  !$omp threadprivate ( cmm_impr_list )
  !$omp threadprivate ( cmm_cmap_list )
  !$omp threadprivate ( cmm_water_list )
  !$omp threadprivate ( cmm_solute_list )
  !$omp threadprivate ( cmm_solute_list_inv )
  !$omp threadprivate ( cmm_solute )
  !$omp threadprivate ( cmm_duplicate )
  !$omp threadprivate ( cmm_H_index )
  !$omp threadprivate ( cmm_selatoms )

  logical            :: cm_created   = .false.
  integer(8)         :: cm_fast      = 0
  integer(8)         :: cm_slow      = 0
  integer(8)         :: cm_page_in   = 0
  integer(8)         :: cm_page_out  = 0

  ! subroutines
  private :: cm_backup_cache
  private :: cm_recall_cache
  private :: cm_page_in_C1
  private :: cm_page_in_L1
  private :: cm_page_in_I1
  private :: cm_page_in_IN
  private :: cm_page_in_R1
  private :: cm_page_in_R3
  private :: cm_page_out_C1
  private :: cm_page_out_L1
  private :: cm_page_out_I1
  private :: cm_page_out_IN
  private :: cm_page_out_R1
  private :: cm_page_out_R3
  private :: cm_init_pages
  private :: cm_get_page
  private :: cm_dirty_page
  public  :: cm_get_cache_name

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_initialize
  !> @brief        initialize variables
  !! @authors      NT
  !! @param[in]    cache_path : path for cache files
  !! @param[in]    cache_name : cache file prefix
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_initialize(cache_path, cache_name)

    ! formal arguments
    character(*),            intent(in)    :: cache_path
    character(*),            intent(in)    :: cache_name


    cm_cache_dir  = cache_path
    cm_cache_name = cache_name

    write(MsgOut,'(a)') 'Cm_Initialize> cache file : '
    write(MsgOut,'(a)') '  cache directory   = ' // trim(cm_cache_dir)
    write(MsgOut,'(a)') '  cache file prefix = ' // trim(cm_cache_name)
    write(MsgOut,'(a)')   ' '

    return

  end subroutine cm_initialize

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_finalize
  !> @brief        delete cache files and print performance information
  !! @authors      NT
  !! @param[in]    cleanup : flag for clean-up cache file
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_finalize(cleanup)

    ! formal arguments
    logical,                 intent(in)    :: cleanup

    ! local variables
    integer                  :: iret
    integer                  :: system


    ! clean-up cache files
    !
    if (cleanup) then

      write(MsgOut,'(a)') 'Cm_Finalize> cleanup cache files.'
      write(MsgOut,'(a)') ' '

      iret =system('rm -f '//trim(cm_cache_dir)//'/'//trim(cm_cache_name)//'_*')

      if (iret /= 0) then
        write(MsgOut,'(a)') 'Cm_Delete_Cache> ERROR: Cannot delete cache file'
        write(MsgOut,'(a)') '   ** Below combination is often not work.     **'
        write(MsgOut,'(a)') '   ** [debug compile + OpenMP + system() call] **'
        write(MsgOut,'(a)') ' '
      end if

    end if
  

    ! print debug information
    !
    write(MsgOut,'(a,i14)') 'Cm_Finalize> Fast search : ', cm_fast
    write(MsgOut,'(a,i14)') '             Slow search : ', cm_slow
    write(MsgOut,'(a,i14)') '             Page-in     : ', cm_page_in
    write(MsgOut,'(a,i14)') '             Page-out    : ', cm_page_out
    write(MsgOut,'(a)')   ' '

    return

  end subroutine cm_finalize

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_create_molecules
  !> @brief        create cached molecule from physical storage
  !! @authors      NT
  !! @return       logical : flag for need to create cache files
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function cm_create_molecules()

    ! return
    logical                  :: cm_create_molecules

    ! local variables
    integer                  :: hdr_in
    character(MaxFilename)   :: file
    logical                  :: lexist


    if (cm_created) &
      call error_msg('Cm_Create_Molecules> Cached molecule is already created.')


    ! allocate page memory
    !

    !$omp parallel

    allocate(cmm_atom_no      (   CMPageSize, CMNumPages), &
             cmm_atom_cls_no  (   CMPageSize, CMNumPages), &
             cmm_atom_cls_name(   CMPageSize, CMNumPages), &
             cmm_atom_name    (   CMPageSize, CMNumPages), &
             cmm_atom_coord   (3, CMPageSize, CMNumPages), &
             cmm_atom_velocity(3, CMPageSize, CMNumPages), &
             cmm_atom_refcoord(3, CMPageSize, CMNumPages), &
             cmm_charge       (   CMPageSize, CMNumPages), &
             cmm_mass         (   CMPageSize, CMNumPages), &
             cmm_residue_no   (   CMPageSize, CMNumPages), &
             cmm_residue_name (   CMPageSize, CMNumPages), &
             cmm_segment_name (   CMPageSize, CMNumPages), &
             cmm_molecule_no  (   CMPageSize, CMNumPages), &
             cmm_bond_list    (2, CMPageSize, CMNumPages), &
             cmm_angl_list    (3, CMAnglSize, CMNumPages), &
             cmm_dihe_list    (4, CMPageSize, CMNumPages), &
             cmm_impr_list    (4, CMPageSize, CMNumPages), &
             cmm_cmap_list    (8, CMPageSize, CMNumPages))

    !$omp end parallel
    !$omp barrier


    ! initialize page information
    !
    
    !$omp parallel

    call cm_init_pages(cmp_atom_no,         CMPageSize)
    call cm_init_pages(cmp_atom_cls_no,     CMPageSize)
    call cm_init_pages(cmp_atom_cls_name,   CMPageSize)
    call cm_init_pages(cmp_atom_name,       CMPageSize)
    call cm_init_pages(cmp_atom_coord,      CMPageSize)
    call cm_init_pages(cmp_atom_velocity,   CMPageSize)
    call cm_init_pages(cmp_atom_refcoord,   CMPageSize)
    call cm_init_pages(cmp_charge,          CMPageSize)
    call cm_init_pages(cmp_mass,            CMPageSize)
    call cm_init_pages(cmp_residue_no,      CMPageSize)
    call cm_init_pages(cmp_residue_name,    CMPageSize)
    call cm_init_pages(cmp_segment_name,    CMPageSize)
    call cm_init_pages(cmp_molecule_no,     CMPageSize)
    call cm_init_pages(cmp_bond_list,       CMPageSize)
    call cm_init_pages(cmp_angl_list,       CMAnglSize)
    call cm_init_pages(cmp_dihe_list,       CMPageSize)
    call cm_init_pages(cmp_impr_list,       CMPageSize)
    call cm_init_pages(cmp_cmap_list,       CMPageSize)

    !$omp end parallel
    !$omp barrier
    

    ! create cache files
    !

    file = cm_get_cache_name(0, "HEADER")
    inquire(file=file, exist=lexist)

    if (lexist) then

      write(MsgOut,'(a)') 'Cm_Create_Molecules> cache files are already exist.'
      write(MsgOut,'(a)') '           Reuse cache file. '//&
           trim(cm_cache_name)//'...'
      write(msgOut,'(a)') ''

      call open_file(hdr_in, file, IOFileInput)
      read(hdr_in,*) cm_num_atoms
      read(hdr_in,*) cm_num_bonds
      read(hdr_in,*) cm_num_angles
      read(hdr_in,*) cm_num_dihedrals
      read(hdr_in,*) cm_num_impropers
      read(hdr_in,*) cm_num_cmaps
      call close_file(hdr_in)

    end if

    cm_created   = .true.
    cm_fast      = 0
    cm_slow      = 0
    cm_page_in   = 0
    cm_page_out  = 0

    cm_create_molecules = .not. lexist

    return

  end function cm_create_molecules

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_delete_molecules
  !> @brief        delete cached molecule
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_delete_molecules


    if (.not. cm_created) &
      call error_msg('Cm_delete_molecules> Cached molecule is not created yet.')


    ! synchronize page-mem and cache-file
    !
    call cm_sync_atom_coord
    call cm_sync_atom_refcoord
    call cm_sync_atom_velocity
    call cm_sync_molecule_no


    ! deallocate memory
    !
    !$omp parallel

    deallocate( &
         cmm_atom_no,       &
         cmm_atom_cls_no,   &
         cmm_atom_cls_name, &
         cmm_atom_name,     &
         cmm_atom_coord,    &
         cmm_atom_velocity, &
         cmm_atom_refcoord, &
         cmm_charge,        &
         cmm_mass,          &
         cmm_residue_no,    &
         cmm_residue_name,  &
         cmm_segment_name,  &
         cmm_molecule_no,   &
         cmm_bond_list,     &
         cmm_angl_list,     &
         cmm_dihe_list,     &
         cmm_impr_list,     &
         cmm_cmap_list)

    !$omp end parallel


    cm_created = .false.

    return

  end subroutine cm_delete_molecules

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_regist_molecules
  !> @brief        regist cached molecule to physical storage
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_regist_molecules

    ! local variables
    integer                  :: hdr_in


    call open_file(hdr_in, cm_get_cache_name(0, "HEADER"), IOFileOutputNew)
    write(hdr_in,*) cm_num_atoms
    write(hdr_in,*) cm_num_bonds
    write(hdr_in,*) cm_num_angles
    write(hdr_in,*) cm_num_dihedrals
    write(hdr_in,*) cm_num_impropers
    write(hdr_in,*) cm_num_cmaps
    call close_file(hdr_in)

    return

  end subroutine cm_regist_molecules

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_atom_no
  !> @brief        get molecule%atom_no
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @return       atom serial number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function cm_atom_no(i)

    ! return value
    integer                  :: cm_atom_no

    ! formal arguments
    integer,                 intent(in)    :: i

    ! local variables
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in


    call cm_get_page(cmp_atom_no, i, mem_no, mem_index, &
                     page_in, pin_offset)

    if (page_in) &
      call cm_page_in_I1(pin_offset, &
                     'atom_no', cmm_atom_no(:,mem_no))

    cm_atom_no = cmm_atom_no(mem_index, mem_no)

    return

  end function cm_atom_no

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_atom_cls_no
  !> @brief        get molecule%atom_cls_no
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @return       atom class number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function cm_atom_cls_no(i)

    ! return value
    integer                  :: cm_atom_cls_no

    ! formal arguments
    integer,                 intent(in)    :: i

    ! local variables
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in


    call cm_get_page(cmp_atom_cls_no, i, mem_no, mem_index, &
                     page_in, pin_offset)

    if (page_in) &
      call cm_page_in_I1(pin_offset, &
                     'atom_cls_no', cmm_atom_cls_no(:,mem_no))

    cm_atom_cls_no = cmm_atom_cls_no(mem_index, mem_no)

    return

  end function cm_atom_cls_no

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_atom_cls_name
  !> @brief        get molecule%atom_cls_name
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @return       atom class name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function cm_atom_cls_name(i)

    ! return value
    character(6)             :: cm_atom_cls_name

    ! formal arguments
    integer,                 intent(in)    :: i

    ! local variables
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in


    call cm_get_page(cmp_atom_cls_name, &
                     i, mem_no, mem_index, page_in, pin_offset)

    if (page_in) &
      call cm_page_in_C1(pin_offset, &
                     'atom_cls_name', cmm_atom_cls_name(:,mem_no))

    cm_atom_cls_name = cmm_atom_cls_name(mem_index, mem_no)

    return

  end function cm_atom_cls_name

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_atom_name
  !> @brief        get molecule%atom_name
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @return       atom name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function cm_atom_name(i)

    ! return value
    character(4)             :: cm_atom_name

    ! formal arguments
    integer,                 intent(in)    :: i

    ! local variables
    integer                  :: file
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in


    call cm_get_page(cmp_atom_name, &
                     i, mem_no, mem_index, page_in, pin_offset)

    if (page_in) &
      call cm_page_in_C1(pin_offset, &
                     'atom_name', cmm_atom_name(:,mem_no))

    cm_atom_name = cmm_atom_name(mem_index, mem_no)

    return

  end function cm_atom_name

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_atom_coord
  !> @brief        get molecule%atom_coord
  !! @authors      NT
  !! @param[in]    j : index of axis
  !! @param[in]    i : index of atom
  !! @return       atom coordinate
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function cm_atom_coord(j, i)

    ! return value
    real(wp)                 :: cm_atom_coord

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i

    ! local variables
    integer                  :: file
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in

    integer                  :: pout_mem_no, pout_offset
    logical                  :: page_out


    call cm_get_page(cmp_atom_coord, &
                     i, mem_no, mem_index, page_in, pin_offset, &
                     page_out, pout_mem_no, pout_offset)

    if (page_out) &
      call cm_page_out_R3(pout_offset, &
                          'atom_coord', cmm_atom_coord(:,:,pout_mem_no))

    if (page_in) &
      call cm_page_in_R3 (pin_offset, &
                          'atom_coord', cmm_atom_coord(:,:,mem_no))

    cm_atom_coord = cmm_atom_coord(j, mem_index, mem_no)

    return

  end function cm_atom_coord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_atom_refcoord
  !> @brief        get molecule%atom_refcoord
  !! @authors      NT
  !! @param[in]    j : index of axis
  !! @param[in]    i : index of atom
  !! @return       atom reference coordinate
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function cm_atom_refcoord(j, i)

    ! return value
    real(wp)                 :: cm_atom_refcoord

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i

    ! local variables
    integer                  :: file
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in


    call cm_get_page(cmp_atom_refcoord, &
                     i, mem_no, mem_index, page_in, pin_offset)

    if (page_in) &
      call cm_page_in_R3 (pin_offset, &
                          'atom_refcoord', cmm_atom_refcoord(:,:,mem_no))

    cm_atom_refcoord = cmm_atom_refcoord(j, mem_index, mem_no)

    return

  end function cm_atom_refcoord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_atom_velocity
  !> @brief        get molecule%atom_velocity
  !! @authors      NT
  !! @param[in]    j : index of axis
  !! @param[in]    i : index of atom
  !! @return       atom velocity
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function cm_atom_velocity(j, i)

    ! return value
    real(wp)                 :: cm_atom_velocity

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i

    ! local variables
    integer                  :: file
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in

    integer                  :: pout_mem_no, pout_offset
    logical                  :: page_out


    call cm_get_page(cmp_atom_velocity, &
                     i, mem_no, mem_index, page_in, pin_offset, &
                     page_out, pout_mem_no, pout_offset)

    if (page_out) &
      call cm_page_out_R3(pout_offset, &
                          'atom_velocity', cmm_atom_velocity(:,:,pout_mem_no))

    if (page_in) &
      call cm_page_in_R3(pin_offset, &
                         'atom_velocity', cmm_atom_velocity(:,:,mem_no))

    cm_atom_velocity = cmm_atom_velocity(j, mem_index, mem_no)

    return

  end function cm_atom_velocity

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_charge
  !> @brief        get molecule%charge
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @return       atom charge
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function cm_charge(i)

    ! return value
    real(wp)                 :: cm_charge

    ! formal arguments
    integer,                 intent(in)    :: i

    ! local variables
    integer                  :: file
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in


    call cm_get_page(cmp_charge, &
                     i, mem_no, mem_index, page_in, pin_offset)

    if (page_in) &
      call cm_page_in_R1(pin_offset, &
                         'charge', cmm_charge(:,mem_no))

    cm_charge = cmm_charge(mem_index, mem_no)

    return

  end function cm_charge

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_mass
  !> @brief        get molecule%mass
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @return       atom mass
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function cm_mass(i)

    ! return value
    real(wp)                 :: cm_mass

    ! formal arguments
    integer,                 intent(in)    :: i

    ! local variables
    integer                  :: file
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in


    call cm_get_page(cmp_mass, &
                     i, mem_no, mem_index, page_in, pin_offset)

    if (page_in) &
      call cm_page_in_R1(pin_offset, &
                         'mass', cmm_mass(:,mem_no))

    cm_mass = cmm_mass(mem_index, mem_no)

    return

  end function cm_mass

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_residue_no
  !> @brief        get molecule%residue_no
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @return       residue number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function cm_residue_no(i)

    ! return value
    integer                  :: cm_residue_no

    ! formal arguments
    integer,                 intent(in)    :: i

    ! local variables
    integer                  :: file
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in


    call cm_get_page(cmp_residue_no, &
                     i, mem_no, mem_index, page_in, pin_offset)

    if (page_in) &
      call cm_page_in_I1(pin_offset, &
                         'residue_no', cmm_residue_no(:,mem_no))

    cm_residue_no = cmm_residue_no(mem_index, mem_no)

    return

  end function cm_residue_no

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_residue_name
  !> @brief        get molecule%residue_name
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @return       residue name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function cm_residue_name(i)

    ! return value
    character(6)             :: cm_residue_name

    ! formal arguments
    integer,                 intent(in)    :: i

    ! local variables
    integer                  :: file
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in


    call cm_get_page(cmp_residue_name, &
                     i, mem_no, mem_index, page_in, pin_offset)

    if (page_in) &
      call cm_page_in_C1(pin_offset, &
                         'residue_name', cmm_residue_name(:,mem_no))

    cm_residue_name = cmm_residue_name(mem_index, mem_no)

    return

  end function cm_residue_name

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_segment_name
  !> @brief        get molecule%segment_name
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @return       segment name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function cm_segment_name(i)

    ! return value
    character(4)             :: cm_segment_name

    ! formal arguments
    integer,                 intent(in)    :: i

    ! local variables
    integer                  :: file
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in


    call cm_get_page(cmp_segment_name, &
                     i, mem_no, mem_index, page_in, pin_offset)

    if (page_in) &
      call cm_page_in_C1(pin_offset, &
                         'segment_name', cmm_segment_name(:,mem_no))

    cm_segment_name = cmm_segment_name(mem_index, mem_no)

    return

  end function cm_segment_name

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_molecule_no
  !> @brief        get molecule%molecule_no
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @return       molecule number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function cm_molecule_no(i)

    ! return value
    integer                  :: cm_molecule_no

    ! formal arguments
    integer,                 intent(in)    :: i

    ! local variables
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in

    integer                  :: pout_mem_no, pout_offset
    logical                  :: page_out


    call cm_get_page(cmp_molecule_no, &
                     i, mem_no, mem_index, page_in, pin_offset, &
                     page_out, pout_mem_no, pout_offset)

    if (page_out) &
      call cm_page_out_I1(pout_offset, &
                          'molecule_no', cmm_molecule_no(:,pout_mem_no))

    if (page_in) &
      call cm_page_in_I1 (pin_offset, &
                          'molecule_no', cmm_molecule_no(:,mem_no))

    cm_molecule_no = cmm_molecule_no(mem_index, mem_no)

    return

  end function cm_molecule_no

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_bond_list
  !> @brief        get molecule%bond_list
  !! @authors      NT
  !! @param[in]    j : index of bond atom
  !! @param[in]    i : index of bond
  !! @return       index of atom
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function cm_bond_list(j, i)

    ! return value
    integer                  :: cm_bond_list

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i

    ! local variables
    integer                  :: file
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in


    call cm_get_page(cmp_bond_list, i, mem_no, mem_index, page_in, pin_offset)

    if (page_in) &
      call cm_page_in_IN(2, pin_offset, &
                         'bond_list', cmm_bond_list(:,:,mem_no))

    cm_bond_list = cmm_bond_list(j, mem_index, mem_no)

    return

  end function cm_bond_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_angl_list
  !> @brief        get molecule%angl_list
  !! @authors      NT
  !! @param[in]    j : index of angle atom
  !! @param[in]    i : index of angle
  !! @return       index of atom
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function cm_angl_list(j, i)

    ! return value
    integer                  :: cm_angl_list

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i

    ! local variables
    integer                  :: file
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in


    call cm_get_page(cmp_angl_list, i, mem_no, mem_index, page_in, pin_offset)

    if (page_in) then

      !$omp critical
      call open_binary_file(file, &
           cm_get_cache_name(pin_offset, 'angl_list'), IOFileInput)
      !$omp end critical

      read(file) cmm_angl_list(1:3,1:CMAnglSize,mem_no)

      !$omp critical
      call close_file(file)
      !$omp end critical

    end if

    cm_angl_list = cmm_angl_list(j, mem_index, mem_no)

    return

  end function cm_angl_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_dihe_list
  !> @brief        get molecule%dihe_list
  !! @authors      NT
  !! @param[in]    j : index of dihedral atom
  !! @param[in]    i : index of dihedral
  !! @return       index of atom
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function cm_dihe_list(j, i)

    ! return value
    integer                  :: cm_dihe_list

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i

    ! local variables
    integer                  :: file
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in


    call cm_get_page(cmp_dihe_list, i, mem_no, mem_index, page_in, pin_offset)

    if (page_in) &
      call cm_page_in_IN(4, pin_offset, &
                         'dihe_list', cmm_dihe_list(:,:,mem_no))

    cm_dihe_list = cmm_dihe_list(j, mem_index, mem_no)

    return

  end function cm_dihe_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_impr_list
  !> @brief        get molecule%impr_list
  !! @authors      NT
  !! @param[in]    j : index of improper atom
  !! @param[in]    i : index of improper
  !! @return       index of atom
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function cm_impr_list(j, i)

    ! return value
    integer                  :: cm_impr_list

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i

    ! local variables
    integer                  :: file
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in


    call cm_get_page(cmp_impr_list, i, mem_no, mem_index, page_in, pin_offset)

    if (page_in) &
      call cm_page_in_IN(4, pin_offset, &
                         'impr_list', cmm_impr_list(:,:,mem_no))

    cm_impr_list = cmm_impr_list(j, mem_index, mem_no)

    return

  end function cm_impr_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_cmap_list
  !> @brief        get molecule%cmap_list
  !! @authors      NT
  !! @param[in]    j : index of cmap atom
  !! @param[in]    i : index of cmap
  !! @return       index of atom
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function cm_cmap_list(j, i)

    ! return value
    integer                  :: cm_cmap_list

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i

    ! local variables
    integer                  :: file
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in


    call cm_get_page(cmp_cmap_list, i, mem_no, mem_index, page_in, pin_offset)

    if (page_in) &
      call cm_page_in_IN(8, pin_offset, &
                         'cmap_list', cmm_cmap_list(:,:,mem_no))

    cm_cmap_list = cmm_cmap_list(j, mem_index, mem_no)

    return

  end function cm_cmap_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_set_atom_coord
  !> @brief        set molecule%atom_coord
  !! @authors      NT
  !! @param[in]    j   : index of axis
  !! @param[in]    i   : index of atom
  !! @param[in]    val : atom coordinate
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_set_atom_coord(j, i, val)

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i
    real(wp),                intent(in)    :: val

    ! local variables
    integer                  :: file
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in

    integer                  :: pout_mem_no, pout_offset
    logical                  :: page_out


#ifdef OMP
    if (omp_get_thread_num() /= 0) &
      call error_msg("Cm_Set_*** () ERROR : multi-thread is not supported.")
#endif

    call cm_get_page(cmp_atom_coord, &
                     i, mem_no, mem_index, page_in, pin_offset, &
                     page_out, pout_mem_no, pout_offset)

    if (page_out) &
      call cm_page_out_R3(pout_offset, &
                          'atom_coord', cmm_atom_coord(:,:,pout_mem_no))

    if (page_in) &
      call cm_page_in_R3 (pin_offset, &
                          'atom_coord', cmm_atom_coord(:,:,mem_no))

    cmm_atom_coord(j, mem_index, mem_no) = val

    call cm_dirty_page(cmp_atom_coord)

    return

  end subroutine cm_set_atom_coord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_set_atom_refcoord
  !> @brief        set molecule%atom_refcoord
  !! @authors      NT
  !! @param[in]    j   : index of axis
  !! @param[in]    i   : index of atom
  !! @param[in]    val : atom coordinate
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_set_atom_refcoord(j, i, val)

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i
    real(wp),                intent(in)    :: val

    ! local variables
    integer                  :: file
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in

    integer                  :: pout_mem_no, pout_offset
    logical                  :: page_out


#ifdef OMP
    if (omp_get_thread_num() /= 0) &
      call error_msg("Cm_Set_*** () ERROR : multi-thread is not supported.")
#endif

    call cm_get_page(cmp_atom_refcoord, &
                     i, mem_no, mem_index, page_in, pin_offset, &
                     page_out, pout_mem_no, pout_offset)

    if (page_out) &
      call cm_page_out_R3(pout_offset, &
                          'atom_refcoord', cmm_atom_refcoord(:,:,pout_mem_no))

    if (page_in) &
      call cm_page_in_R3 (pin_offset, &
                          'atom_refcoord', cmm_atom_refcoord(:,:,mem_no))

    cmm_atom_refcoord(j, mem_index, mem_no) = val

    call cm_dirty_page(cmp_atom_refcoord)

    return

  end subroutine cm_set_atom_refcoord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_set_atom_velocity
  !> @brief        set molecule%atom_velocity
  !! @authors      NT
  !! @param[in]    j   : index of axis
  !! @param[in]    i   : index of atom
  !! @param[in]    val : atom velocity
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_set_atom_velocity(j, i, val)

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i
    real(wp),                intent(in)    :: val

    ! local variables
    integer                  :: file
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in

    integer                  :: pout_mem_no, pout_offset
    logical                  :: page_out


#ifdef OMP
    if (omp_get_thread_num() /= 0) &
      call error_msg("Cm_Set_*** () ERROR : multi-thread is not supported.")
#endif

    call cm_get_page(cmp_atom_velocity, &
                     i, mem_no, mem_index, page_in, pin_offset, &
                     page_out, pout_mem_no, pout_offset)

    if (page_out) &
      call cm_page_out_R3(pout_offset, &
                          'atom_velocity', cmm_atom_velocity(:,:,pout_mem_no))

    if (page_in) &
      call cm_page_in_R3 (pin_offset, &
                          'atom_velocity', cmm_atom_velocity(:,:,mem_no))

    cmm_atom_velocity(j, mem_index, mem_no) = val

    call cm_dirty_page(cmp_atom_velocity)

    return

  end subroutine cm_set_atom_velocity

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_set_molecule_no
  !> @brief        set molecule%molecule_no
  !! @authors      NT
  !! @param[in]    i   : index of atom
  !! @param[in]    val : molecule number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_set_molecule_no(i, val)

    ! formal arguments
    integer,                 intent(in)    :: i
    integer,                 intent(in)    :: val

    ! local variables
    integer                  :: file
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in

    integer                  :: pout_mem_no, pout_offset
    logical                  :: page_out


#ifdef OMP
    if (omp_get_thread_num() /= 0) &
      call error_msg("Cm_Set_*** () ERROR : multi-thread is not supported.")
#endif

    call cm_get_page(cmp_molecule_no, &
                     i, mem_no, mem_index, page_in, pin_offset, &
                     page_out, pout_mem_no, pout_offset)

    if (page_out) &
      call cm_page_out_I1(pout_offset, &
                          'molecule_no', cmm_molecule_no(:,pout_mem_no))

    if (page_in) &
      call cm_page_in_I1 (pin_offset, &
                          'molecule_no', cmm_molecule_no(:,mem_no))

    cmm_molecule_no(mem_index, mem_no) = val

    call cm_dirty_page(cmp_molecule_no)

    return

  end subroutine cm_set_molecule_no

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_sync_atom_coord
  !> @brief        synchronize page-mem and cache-file about atom_coord
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_sync_atom_coord

    ! local variables
    integer                  :: i, mem_no, offset, file


    ! page-out last modified data
    !
    do i = 2, CMNumPages+1
      if (.not. cmp_atom_coord(i)%dirty) &
        cycle

      mem_no = cmp_atom_coord(i)%mem_no
      offset = cmp_atom_coord(i)%offset

      !$omp critical
      call open_binary_file(file, &
                           cm_get_cache_name(offset, 'atom_coord'), &
                           IOFileOutputReplace)
      !$omp end critical

      write(file) cmm_atom_coord(1:3, 1:CMPageSize, mem_no)

      !$omp critical
      call close_file(file)
      !$omp end critical

    end do

    return

  end subroutine cm_sync_atom_coord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_sync_atom_refcoord
  !> @brief        synchronize page-mem and cache-file about atom_refcoord
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_sync_atom_refcoord

    ! local variables
    integer                  :: i, mem_no, offset, file


    ! page-out last modified data
    !
    do i = 2, CMNumPages+1
      if (.not. cmp_atom_refcoord(i)%dirty) &
        cycle

      mem_no = cmp_atom_refcoord(i)%mem_no
      offset = cmp_atom_refcoord(i)%offset

      !$omp critical
      call open_binary_file(file, &
                           cm_get_cache_name(offset, 'atom_refcoord'), &
                           IOFileOutputReplace)
      !$omp end critical

      write(file) cmm_atom_refcoord(1:3, 1:CMPageSize, mem_no)

      !$omp critical
      call close_file(file)
      !$omp end critical

    end do

    return

  end subroutine cm_sync_atom_refcoord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_sync_atom_velocity
  !> @brief        synchronize page-mem and cache-file about atom_velocity
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_sync_atom_velocity

    ! local variables
    integer                  :: i, mem_no, offset, file


    ! page-out last modified data
    !
    do i = 2, CMNumPages+1
      if (.not. cmp_atom_velocity(i)%dirty) &
        cycle

      mem_no = cmp_atom_velocity(i)%mem_no
      offset = cmp_atom_velocity(i)%offset

      !$omp critical
      call open_binary_file(file, &
                           cm_get_cache_name(offset, 'atom_velocity'), &
                           IOFileOutputReplace)
      !$omp end critical

      write(file) cmm_atom_velocity(1:3, 1:CMPageSize, mem_no)

      !$omp critical
      call close_file(file)
      !$omp end critical

    end do

    return

  end subroutine cm_sync_atom_velocity

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_sync_molecule_no
  !> @brief        synchronize page-mem and cache-file about molecule_no
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_sync_molecule_no

    ! local variables
    integer                  :: i, mem_no, offset, file


    ! page-out last modified data
    !
    do i = 2, CMNumPages+1
      if (.not. cmp_molecule_no(i)%dirty) &
        cycle

      mem_no = cmp_molecule_no(i)%mem_no
      offset = cmp_molecule_no(i)%offset

      !$omp critical
      call open_binary_file(file, &
                           cm_get_cache_name(offset, 'molecule_no'), &
                           IOFileOutputReplace)
      !$omp end critical

      write(file) cmm_molecule_no(1:CMPageSize, mem_no)

      !$omp critical
      call close_file(file)
      !$omp end critical

    end do

    return

  end subroutine cm_sync_molecule_no

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_backup_atom_coord
  !> @brief        backup molecule%atom_coord
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_backup_atom_coord


    call cm_backup_cache('atom_coord')
    return

  end subroutine cm_backup_atom_coord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_backup_atom_velocity
  !> @brief        backup molecule%atom_velocity
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_backup_atom_velocity


    call cm_backup_cache('atom_velocity')
    return

  end subroutine cm_backup_atom_velocity

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_recall_atom_coord
  !> @brief        recall molecule%atom_coord
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_recall_atom_coord


    call cm_recall_cache('atom_coord')
    return

  end subroutine cm_recall_atom_coord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_recall_atom_velocity
  !> @brief        recall molecule%atom_velocity
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_recall_atom_velocity


    call cm_recall_cache('atom_velocity')
    return

  end subroutine cm_recall_atom_velocity

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_create_water_list
  !> @brief        allocate enefunc%table%water_list
  !! @authors      NT
  !! @param[in]    var_size : array size
  !! @param[in]    lexist   : flag for cache file is exist or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_create_water_list(var_size, lexist)

    ! formal arguments
    integer,                 intent(in)    :: var_size
    logical,                 intent(out)   :: lexist

    ! local variables
    integer                  :: i, file_out, num
    integer,    allocatable  :: water_list(:,:)


    cm_num_water_list = var_size


    ! allocate page memory
    !
    !$omp parallel
    allocate(cmm_water_list(3, CMPageSize, CMNumPages))
    call cm_init_pages(cmp_water_list,     CMPageSize)
    !$omp end parallel
    !$omp barrier


    ! create cache file
    !
    inquire(file=cm_get_cache_name(1, 'water_list'), exist=lexist)
    if (lexist) then
      write(MsgOut,'(a)') 'Cm_Create_Water_List> Re-use cache file. '
      return
    end if

    allocate(water_list(3, CMPageSize))
    water_list(1:3,1:CMPageSize)   = 0

    do i = 0, (var_size / CMPageSize)

      num = i * CMPageSize + 1

      call open_binary_file(file_out, &
                            cm_get_cache_name(num, 'water_list'), &
                            IOFileOutputNew)
      write(file_out) water_list(1:3, 1:CMPageSize)
      call close_file(file_out)

    end do

    deallocate(water_list)

    return

  end subroutine cm_create_water_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_create_solute_list
  !> @brief        allocate enefunc%table%solute_list
  !! @authors      NT
  !! @param[in]    var_size : array size
  !! @param[in]    lexist   : flag for cache file is exist or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_create_solute_list(var_size, lexist)

    ! formal arguments
    integer,                 intent(in)    :: var_size
    logical,                 intent(out)   :: lexist

    ! local variables
    integer                  :: i, file_out, num
    integer,    allocatable  :: solute_list(:)


    cm_num_solute_list = var_size


    ! allocate page memory
    !
    !$omp parallel
    allocate(cmm_solute_list(CMPageSize, CMNumPages))
    call cm_init_pages(cmp_solute_list,  CMPageSize)
    !$omp end parallel
    !$omp barrier


    ! create cache file
    !
    inquire(file=cm_get_cache_name(1, 'solute_list'), exist=lexist)
    if (lexist) then
      write(MsgOut,'(a)') 'Cm_Create_Solute_List> Re-use cache file. '
      return
    end if

    allocate(solute_list(CMPageSize))
    solute_list(1:CMPageSize)   = 0

    do i = 0, (var_size / CMPageSize)

      num = i * CMPageSize + 1

      call open_binary_file(file_out, &
                            cm_get_cache_name(num, 'solute_list'), &
                            IOFileOutputNew)
      write(file_out) solute_list(1:CMPageSize)
      call close_file(file_out)

    end do

    deallocate(solute_list)

    return

  end subroutine cm_create_solute_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_create_solute_list_inv
  !> @brief        allocate enefunc%table%solute_list_inv
  !! @authors      NT
  !! @param[in]    var_size : array size
  !! @param[in]    lexist   : flag for cache file is exist or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_create_solute_list_inv(var_size, lexist)

    ! formal arguments
    integer,                 intent(in)    :: var_size
    logical,                 intent(out)   :: lexist

    ! local variables
    integer                  :: i, file_out, num
    integer,    allocatable  :: solute_list_inv(:)


    cm_num_solute_list_inv = var_size


    ! allocate page memory
    !
    !$omp parallel
    allocate(cmm_solute_list_inv(CMPageSize, CMNumPages))
    call cm_init_pages(cmp_solute_list_inv,  CMPageSize)
    !$omp end parallel
    !$omp barrier


    ! create cache file
    !
    inquire(file=cm_get_cache_name(1, 'solute_list_inv'), exist=lexist)
    if (lexist) then
      write(MsgOut,'(a)') 'Cm_Create_Solute_List_Inv> Re-use cache file. '
      return
    end if

    allocate(solute_list_inv(CMPageSize))
    solute_list_inv(1:CMPageSize)   = 0

    do i = 0, (var_size / CMPageSize)

      num = i * CMPageSize + 1

      call open_binary_file(file_out, &
                            cm_get_cache_name(num, 'solute_list_inv'), &
                            IOFileOutputNew)
      write(file_out) solute_list_inv(1:CMPageSize)
      call close_file(file_out)

    end do

    deallocate(solute_list_inv)

    return

  end subroutine cm_create_solute_list_inv

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_create_solute
  !> @brief        allocate solute
  !! @authors      NT
  !! @param[in]    var_size : array size
  !! @param[in]    lexist   : flag for cache file is exist or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_create_solute(var_size, lexist)

    ! formal arguments
    integer,                 intent(in)    :: var_size
    logical,                 intent(out)   :: lexist

    ! local variables
    integer                  :: i, file_out, num
    logical,    allocatable  :: solute(:)


    cm_num_solute = var_size


    ! allocate page memory
    !
    !$omp parallel
    allocate(cmm_solute(CMPageSize, CMNumPages))
    call cm_init_pages(cmp_solute,  CMPageSize)
    !$omp end parallel
    !$omp barrier


    ! create cache file
    !
    inquire(file=cm_get_cache_name(1, 'solute'), exist=lexist)
    if (lexist) then
      write(MsgOut,'(a)') 'Cm_Create_Solute> Re-use cache file. '
      return
    end if

    allocate(solute(CMPageSize))
    solute(1:CMPageSize)   = .false.

    do i = 0, (var_size / CMPageSize)

      num = i * CMPageSize + 1

      call open_binary_file(file_out, &
                            cm_get_cache_name(num, 'solute'), &
                            IOFileOutputNew)
      write(file_out) solute(1:CMPageSize)
      call close_file(file_out)

    end do

    deallocate(solute)

    return

  end subroutine cm_create_solute

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_delete_water_list
  !> @brief        deallocate enefunc%table%water_list
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_delete_water_list

    ! local variables
    integer                  :: i, mem_no, offset, file


    if (.not. allocated(cmm_water_list)) &
      return

    ! synchronize page-mem and cache-file
    !
    call cm_sync_water_list


    ! deallocate page memory
    !
    !$omp parallel
    deallocate(cmm_water_list)
    !$omp end parallel

    return

  end subroutine cm_delete_water_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_delete_solute_list
  !> @brief        deallocate enefunc%table%solute_list
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_delete_solute_list

    ! local variables
    integer                  :: i, mem_no, offset, file


    if (.not. allocated(cmm_solute_list)) &
      return

    ! synchronize page-mem and cache-file
    !
    call cm_sync_solute_list


    ! deallocate page memory
    !
    !$omp parallel
    deallocate(cmm_solute_list)
    !$omp end parallel

    return

  end subroutine cm_delete_solute_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_delete_solute_list_inv
  !> @brief        deallocate enefunc%table%solute_list_inv
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_delete_solute_list_inv

    ! local variables
    integer                  :: i, mem_no, offset, file


    if (.not. allocated(cmm_solute_list_inv)) &
      return

    ! synchronize page-mem and cache-file
    !
    call cm_sync_solute_list_inv


    ! deallocate page memory
    !
    !$omp parallel
    if (allocated(cmm_solute_list_inv)) &
      deallocate(cmm_solute_list_inv)
    !$omp end parallel

    return

  end subroutine cm_delete_solute_list_inv

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_delete_solute
  !> @brief        deallocate solute
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_delete_solute

    ! local variables
    integer                  :: i, mem_no, offset, file


    if (.not. allocated(cmm_solute)) &
      return

    ! synchronize page-mem and cache-file
    !
    call cm_sync_solute


    ! deallocate page memory
    !
    !$omp parallel
    deallocate(cmm_solute)
    !$omp end parallel

    return

  end subroutine cm_delete_solute

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_water_list
  !> @brief        get enefunc%table%water_list
  !! @authors      NT
  !! @param[in]    j : index of water list atom
  !! @param[in]    i : index of water list
  !! @return       index of atom
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function cm_water_list(j, i)

    ! return value
    integer                  :: cm_water_list

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i

    ! local variables
    integer                  :: file
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in

    integer                  :: pout_mem_no, pout_offset
    logical                  :: page_out


    call cm_get_page(cmp_water_list, &
                     i, mem_no, mem_index, page_in, pin_offset, &
                     page_out, pout_mem_no, pout_offset)

    if (page_out) &
      call cm_page_out_IN(3, pout_offset, &
                          'water_list', cmm_water_list(:,:,pout_mem_no))

    if (page_in) &
      call cm_page_in_IN (3, pin_offset, &
                          'water_list', cmm_water_list(:,:,mem_no))

    cm_water_list = cmm_water_list(j, mem_index, mem_no)

    return

  end function cm_water_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_solute_list
  !> @brief        get enefunc%table%solute_list
  !! @authors      NT
  !! @param[in]    i : index of solute list
  !! @return       index of atom
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function cm_solute_list(i)

    ! return value
    integer                  :: cm_solute_list

    ! formal arguments
    integer,                 intent(in)    :: i

    ! local variables
    integer                  :: file
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in

    integer                  :: pout_mem_no, pout_offset
    logical                  :: page_out


    call cm_get_page(cmp_solute_list, &
                     i, mem_no, mem_index, page_in, pin_offset, &
                     page_out, pout_mem_no, pout_offset)

    if (page_out) &
      call cm_page_out_I1(pout_offset, &
                          'solute_list', cmm_solute_list(:,pout_mem_no))

    if (page_in) &
      call cm_page_in_I1 (pin_offset, &
                          'solute_list', cmm_solute_list(:,mem_no))

    cm_solute_list = cmm_solute_list(mem_index, mem_no)

    return

  end function cm_solute_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_solute_list_inv
  !> @brief        get enefunc%table%solute_list_inv
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @return       index of solute list
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function cm_solute_list_inv(i)

    ! return value
    integer                  :: cm_solute_list_inv

    ! formal arguments
    integer,                 intent(in)    :: i

    ! local variables
    integer                  :: file
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in

    integer                  :: pout_mem_no, pout_offset
    logical                  :: page_out


    call cm_get_page(cmp_solute_list_inv, &
                     i, mem_no, mem_index, page_in, pin_offset, &
                     page_out, pout_mem_no, pout_offset)

    if (page_out) &
      call cm_page_out_I1(pout_offset, &
                    'solute_list_inv', cmm_solute_list_inv(:,pout_mem_no))

    if (page_in) &
      call cm_page_in_I1 (pin_offset, &
                    'solute_list_inv', cmm_solute_list_inv(:,mem_no))

    cm_solute_list_inv = cmm_solute_list_inv(mem_index, mem_no)

    return

  end function cm_solute_list_inv

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_solute
  !> @brief        get solute
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @return       flag for solute or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function cm_solute(i)

    ! return value
    logical                  :: cm_solute

    ! formal arguments
    integer,                 intent(in)    :: i

    ! local variables
    integer                  :: file
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in

    integer                  :: pout_mem_no, pout_offset
    logical                  :: page_out


    call cm_get_page(cmp_solute, &
                     i, mem_no, mem_index, page_in, pin_offset, &
                     page_out, pout_mem_no, pout_offset)

    if (page_out) &
      call cm_page_out_L1(pout_offset, &
                          'solute', cmm_solute(:,pout_mem_no))

    if (page_in) &
      call cm_page_in_L1 (pin_offset, &
                          'solute', cmm_solute(:,mem_no))

    cm_solute = cmm_solute(mem_index, mem_no)

    return

  end function cm_solute

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_set_water_list
  !> @brief        set enefunc%table%water_list
  !! @authors      NT
  !! @param[in]    j   : index of water list atom
  !! @param[in]    i   : index of water list
  !! @param[in]    val : index of atom
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_set_water_list(j, i, val)

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i
    integer,                 intent(in)    :: val

    ! local variables
    integer                  :: file
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in

    integer                  :: pout_mem_no, pout_offset
    logical                  :: page_out


#ifdef OMP
    if (omp_get_thread_num() /= 0) &
      call error_msg("Cm_Set_*** () ERROR : multi-thread is not supported.")
#endif

    call cm_get_page(cmp_water_list, &
                     i, mem_no, mem_index, page_in, pin_offset, &
                     page_out, pout_mem_no, pout_offset)

    if (page_out) &
      call cm_page_out_IN(3, pout_offset, &
                          'water_list', cmm_water_list(:,:,pout_mem_no))

    if (page_in) &
      call cm_page_in_IN (3, pin_offset, &
                          'water_list', cmm_water_list(:,:,mem_no))

    cmm_water_list(j, mem_index, mem_no) = val

    call cm_dirty_page(cmp_water_list)

    return

  end subroutine cm_set_water_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_set_solute_list
  !> @brief        set enefunc%table%solute_list
  !! @authors      NT
  !! @param[in]    i   : index of solute list
  !! @param[in]    val : index of atom
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_set_solute_list(i, val)

    ! formal arguments
    integer,                 intent(in)    :: i
    integer,                 intent(in)    :: val

    ! local variables
    integer                  :: file
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in

    integer                  :: pout_mem_no, pout_offset
    logical                  :: page_out


#ifdef OMP
    if (omp_get_thread_num() /= 0) &
      call error_msg("Cm_Set_*** () ERROR : multi-thread is not supported.")
#endif

    call cm_get_page(cmp_solute_list, &
                     i, mem_no, mem_index, page_in, pin_offset, &
                     page_out, pout_mem_no, pout_offset)

    if (page_out) &
      call cm_page_out_I1(pout_offset, &
                          'solute_list', cmm_solute_list(:,pout_mem_no))

    if (page_in) &
      call cm_page_in_I1 (pin_offset, &
                          'solute_list', cmm_solute_list(:,mem_no))

    cmm_solute_list(mem_index, mem_no) = val

    call cm_dirty_page(cmp_solute_list)

    return

  end subroutine cm_set_solute_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_set_solute_list_inv
  !> @brief        set enefunc%table%solute_list_inv
  !! @authors      NT
  !! @param[in]    i   : index of atom
  !! @param[in]    val : index of solute list
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_set_solute_list_inv(i, val)

    ! formal arguments
    integer,                 intent(in)    :: i
    integer,                 intent(in)    :: val

    ! local variables
    integer                  :: file
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in

    integer                  :: pout_mem_no, pout_offset
    logical                  :: page_out


#ifdef OMP
    if (omp_get_thread_num() /= 0) &
      call error_msg("Cm_Set_*** () ERROR : multi-thread is not supported.")
#endif

    call cm_get_page(cmp_solute_list_inv, &
                     i, mem_no, mem_index, page_in, pin_offset, &
                     page_out, pout_mem_no, pout_offset)

    if (page_out) &
      call cm_page_out_I1(pout_offset, &
                    'solute_list_inv', cmm_solute_list_inv(:,pout_mem_no))

    if (page_in) &
      call cm_page_in_I1 (pin_offset, &
                    'solute_list_inv', cmm_solute_list_inv(:,mem_no))

    cmm_solute_list_inv(mem_index, mem_no) = val

    call cm_dirty_page(cmp_solute_list_inv)

    return

  end subroutine cm_set_solute_list_inv

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_set_solute
  !> @brief        set enefunc%table%solute
  !! @authors      NT
  !! @param[in]    i   : index of atom
  !! @param[in]    val : flag for solute or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_set_solute(i, val)

    ! formal arguments
    integer,                 intent(in)    :: i
    logical,                 intent(in)    :: val

    ! local variables
    integer                  :: file
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in

    integer                  :: pout_mem_no, pout_offset
    logical                  :: page_out


#ifdef OMP
    if (omp_get_thread_num() /= 0) &
      call error_msg("Cm_Set_*** () ERROR : multi-thread is not supported.")
#endif

    call cm_get_page(cmp_solute, &
                     i, mem_no, mem_index, page_in, pin_offset, &
                     page_out, pout_mem_no, pout_offset)

    if (page_out) &
      call cm_page_out_L1(pout_offset, &
                          'solute', cmm_solute(:,pout_mem_no))

    if (page_in) &
      call cm_page_in_L1 (pin_offset, &
                          'solute', cmm_solute(:,mem_no))

    cmm_solute(mem_index, mem_no) = val

    call cm_dirty_page(cmp_solute)

    return

  end subroutine cm_set_solute

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_sync_water_list
  !> @brief        synchronize page-mem and cache-file about water_list
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_sync_water_list

    ! local variables
    integer                  :: i, mem_no, offset, file


    ! page-out last modified data
    !
    do i = 2, CMNumPages+1
      if (.not. cmp_water_list(i)%dirty) &
        cycle

      mem_no = cmp_water_list(i)%mem_no
      offset = cmp_water_list(i)%offset

      !$omp critical
      call open_binary_file(file, &
                           cm_get_cache_name(offset, 'water_list'), &
                           IOFileOutputReplace)
      !$omp end critical

      write(file) cmm_water_list(1:3, 1:CMPageSize, mem_no)

      !$omp critical
      call close_file(file)
      !$omp end critical

    end do

    return

  end subroutine cm_sync_water_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_sync_solute_list
  !> @brief        synchronize page-mem and cache-file about solute_list
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_sync_solute_list

    ! local variables
    integer                  :: i, mem_no, offset, file


    ! page-out last modified data
    !
    do i = 2, CMNumPages+1
      if (.not. cmp_solute_list(i)%dirty) &
        cycle

      mem_no = cmp_solute_list(i)%mem_no
      offset = cmp_solute_list(i)%offset

      !$omp critical
      call open_binary_file(file, &
                           cm_get_cache_name(offset, 'solute_list'), &
                           IOFileOutputReplace)
      !$omp end critical

      write(file) cmm_solute_list(1:CMPageSize, mem_no)

      !$omp critical
      call close_file(file)
      !$omp end critical

    end do

    return

  end subroutine cm_sync_solute_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_sync_solute_list_inv
  !> @brief        synchronize page-mem and cache-file about solute_list_inv
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_sync_solute_list_inv

    ! local variables
    integer                  :: i, mem_no, offset, file


    ! page-out last modified data
    !
    do i = 2, CMNumPages+1
      if (.not. cmp_solute_list_inv(i)%dirty) &
        cycle

      mem_no = cmp_solute_list_inv(i)%mem_no
      offset = cmp_solute_list_inv(i)%offset

      !$omp critical
      call open_binary_file(file, &
                           cm_get_cache_name(offset, 'solute_list_inv'), &
                           IOFileOutputReplace)
      !$omp end critical

      write(file) cmm_solute_list_inv(1:CMPageSize, mem_no)

      !$omp critical
      call close_file(file)
      !$omp end critical

    end do

    return

  end subroutine cm_sync_solute_list_inv

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_sync_solute
  !> @brief        synchronize page-mem and cache-file about solute
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_sync_solute

    ! local variables
    integer                  :: i, mem_no, offset, file


    ! page-out last modified data
    !
    do i = 2, CMNumPages+1
      if (.not. cmp_solute_list(i)%dirty) &
        cycle

      mem_no = cmp_solute_list(i)%mem_no
      offset = cmp_solute_list(i)%offset

      !$omp critical
      call open_binary_file(file, &
                           cm_get_cache_name(offset, 'solute'), &
                           IOFileOutputReplace)
      !$omp end critical

      write(file) cmm_solute(1:CMPageSize, mem_no)

      !$omp critical
      call close_file(file)
      !$omp end critical

    end do

    return

  end subroutine cm_sync_solute

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_create_duplicate
  !> @brief        allocate domain_constraint%duplicate
  !! @authors      NT
  !! @param[in]    var_size : array size
  !! @param[in]    lexist   : flag for cache file is exist or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_create_duplicate(var_size, lexist)

    ! formal arguments
    integer,                 intent(in)    :: var_size
    logical,                 intent(out)   :: lexist

    ! local variables
    integer                  :: i, file_out, num
    integer,    allocatable  :: duplicate(:)


    cm_num_h_bond = var_size


    ! allocate page memory
    !
    !$omp parallel
    allocate(cmm_duplicate(CMPageSize, CMNumPages))
    call cm_init_pages(cmp_duplicate,  CMPageSize)
    !$omp end parallel
    !$omp barrier


    ! create cache file
    !
    inquire(file=cm_get_cache_name(1, 'duplicate'), exist=lexist)
    if (lexist) then
      write(MsgOut,'(a)') 'Cm_Create_Duplicate> Re-use cache file. '
      return
    end if

    allocate(duplicate(CMPageSize))
    duplicate(1:CMPageSize)   = 0

    do i = 0, (var_size / CMPageSize)

      num = i * CMPageSize + 1

      call open_binary_file(file_out, &
                            cm_get_cache_name(num, 'duplicate'), &
                            IOFileOutputNew)
      write(file_out) duplicate(1:CMPageSize)
      call close_file(file_out)

    end do

    deallocate(duplicate)

    return

  end subroutine cm_create_duplicate

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_create_H_index
  !> @brief        allocate domain_constraint%H_index
  !! @authors      NT
  !! @param[in]    var_size : array size
  !! @param[in]    lexist   : flag for cache file is exist or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_create_H_index(var_size, lexist)

    ! formal arguments
    integer,                 intent(in)    :: var_size
    logical,                 intent(out)   :: lexist

    ! local variables
    integer                  :: i, file_out, num
    integer,    allocatable  :: H_index(:,:)


    cm_num_h_bond = var_size


    ! allocate page memory
    !
    !$omp parallel
    allocate(cmm_H_index(8, CMPageSize, CMNumPages))
    call cm_init_pages(cmp_H_index,  CMPageSize)
    !$omp end parallel
    !$omp barrier


    ! create cache file
    !
    inquire(file=cm_get_cache_name(1, 'H_index'), exist=lexist)
    if (lexist) then
      write(MsgOut,'(a)') 'Cm_Create_H_index> Re-use cache file. '
      return
    end if

    allocate(H_index(8,CMPageSize))
    H_index(1:8, 1:CMPageSize)   = 0

    do i = 0, (var_size / CMPageSize)

      num = i * CMPageSize + 1

      call open_binary_file(file_out, &
                            cm_get_cache_name(num, 'H_index'), &
                            IOFileOutputNew)
      write(file_out) H_index(1:8, 1:CMPageSize)
      call close_file(file_out)

    end do

    deallocate(H_index)

    return

  end subroutine cm_create_H_index

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_delete_duplicate
  !> @brief        deallocate domain_constraint%duplicate
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_delete_duplicate

    ! local variables
    integer                  :: i, mem_no, offset, file


    if (.not. allocated(cmm_duplicate)) &
      return

    ! synchronize page-mem and cache-file
    !
    call cm_sync_duplicate


    ! deallocate page memory
    !
    !$omp parallel
    deallocate(cmm_duplicate)
    !$omp end parallel

    return

  end subroutine cm_delete_duplicate

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_delete_H_index
  !> @brief        deallocate domain_constraint%H_index
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_delete_H_index

    ! local variables
    integer :: i, mem_no, offset, file


    if (.not. allocated(cmm_H_index)) &
      return

    ! synchronize page-mem and cache-file
    !
    call cm_sync_H_index


    ! deallocate page memory
    !
    !$omp parallel
    deallocate(cmm_H_index)
    !$omp end parallel

    return

  end subroutine cm_delete_H_index

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_duplicate
  !> @brief        get domain_constraint%duplicate
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @return       duplicate
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function cm_duplicate(i)

    ! return value
    integer                  :: cm_duplicate

    ! formal arguments
    integer,                 intent(in)    :: i

    ! local variables
    integer                  :: file
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in

    integer                  :: pout_mem_no, pout_offset
    logical                  :: page_out


    call cm_get_page(cmp_duplicate, &
                     i, mem_no, mem_index, page_in, pin_offset, &
                     page_out, pout_mem_no, pout_offset)

    if (page_out) &
      call cm_page_out_I1(pout_offset, &
                          'duplicate', cmm_duplicate(:,pout_mem_no))

    if (page_in) &
      call cm_page_in_I1 (pin_offset, &
                          'duplicate', cmm_duplicate(:,mem_no))

    cm_duplicate = cmm_duplicate(mem_index, mem_no)

    return

  end function cm_duplicate

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_H_index
  !> @brief        get domain_constraint%H_index
  !! @authors      NT
  !! @param[in]    j : index of duplicate
  !! @param[in]    i : index of atom
  !! @return       H index
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function cm_H_index(j, i)

    ! return value
    integer                  :: cm_H_index

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i

    ! local variables
    integer                  :: file
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in

    integer                  :: pout_mem_no, pout_offset
    logical                  :: page_out


    call cm_get_page(cmp_H_index, &
                     i, mem_no, mem_index, page_in, pin_offset, &
                     page_out, pout_mem_no, pout_offset)

    if (page_out) &
      call cm_page_out_IN(8, pout_offset, &
                         'H_index', cmm_H_index(:,:,pout_mem_no))

    if (page_in) &
      call cm_page_in_IN (8, pin_offset, &
                         'H_index', cmm_H_index(:,:,mem_no))

    cm_H_index = cmm_H_index(j, mem_index, mem_no)

    return

  end function cm_H_index

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_set_duplicate
  !> @brief        set domain_constraint%duplicate
  !! @authors      NT
  !! @param[in]    i   : index of atom
  !! @param[in]    val : duplicate
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_set_duplicate(i, val)

    ! formal arguments
    integer,                 intent(in)    :: i
    integer,                 intent(in)    :: val

    ! local variables
    integer                  :: file
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in

    integer                  :: pout_mem_no, pout_offset
    logical                  :: page_out


#ifdef OMP
    if (omp_get_thread_num() /= 0) &
      call error_msg("Cm_Set_*** () ERROR : multi-thread is not supported.")
#endif

    call cm_get_page(cmp_duplicate, &
                     i, mem_no, mem_index, page_in, pin_offset, &
                     page_out, pout_mem_no, pout_offset)

    if (page_out) &
      call cm_page_out_I1(pout_offset, &
                          'duplicate', cmm_duplicate(:,pout_mem_no))

    if (page_in) &
      call cm_page_in_I1 (pin_offset, &
                          'duplicate', cmm_duplicate(:,mem_no))

    cmm_duplicate(mem_index, mem_no) = val

    call cm_dirty_page(cmp_duplicate)

    return

  end subroutine cm_set_duplicate

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_set_H_index
  !> @brief        set domain_constraint%solute
  !! @authors      NT
  !! @param[in]    j   : index of duplicate
  !! @param[in]    i   : index of atom
  !! @param[in]    val : H index
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_set_H_index(j, i, val)

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i
    integer,                 intent(in)    :: val

    ! local variables
    integer                  :: file
    integer                  :: mem_no, mem_index, pin_offset
    logical                  :: page_in

    integer                  :: pout_mem_no, pout_offset
    logical                  :: page_out


#ifdef OMP
    if (omp_get_thread_num() /= 0) &
      call error_msg("Cm_Set_*** () ERROR : multi-thread is not supported.")
#endif

    call cm_get_page(cmp_H_index, &
                     i, mem_no, mem_index, page_in, pin_offset, &
                     page_out, pout_mem_no, pout_offset)

    if (page_out) &
      call cm_page_out_IN(8, pout_offset, &
                         'H_index', cmm_H_index(:,:,pout_mem_no))

    if (page_in) &
      call cm_page_in_IN (8, pin_offset, &
                         'H_index', cmm_H_index(:,:,mem_no))

    cmm_H_index(j, mem_index, mem_no) = val

    call cm_dirty_page(cmp_H_index)

    return

  end subroutine cm_set_H_index

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_sync_duplicate
  !> @brief        synchronize page-mem and cache-file about duplicate
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_sync_duplicate

    ! local variables
    integer                  :: i, mem_no, offset, file


    ! page-out last modified data
    !
    do i = 2, CMNumPages+1
      if (.not. cmp_duplicate(i)%dirty) &
        cycle

      mem_no = cmp_duplicate(i)%mem_no
      offset = cmp_duplicate(i)%offset

      !$omp critical
      call open_binary_file(file, &
                           cm_get_cache_name(offset, 'duplicate'), &
                           IOFileOutputReplace)
      !$omp end critical

      write(file) cmm_duplicate(1:CMPageSize, mem_no)

      !$omp critical
      call close_file(file)
      !$omp end critical

    end do

    return

  end subroutine cm_sync_duplicate

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_sync_H_index
  !> @brief        synchronize page-mem and cache-file about H_index
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_sync_H_index

    ! local variables
    integer                  :: i, mem_no, offset, file


    ! page-out last modified data
    !
    do i = 2, CMNumPages+1
      if (.not. cmp_H_index(i)%dirty) &
        cycle

      mem_no = cmp_H_index(i)%mem_no
      offset = cmp_H_index(i)%offset

      !$omp critical
      call open_binary_file(file, &
                           cm_get_cache_name(offset, 'H_index'), &
                           IOFileOutputReplace)
      !$omp end critical

      write(file) cmm_H_index(1:8, 1:CMPageSize, mem_no)

      !$omp critical
      call close_file(file)
      !$omp end critical

    end do

    return

  end subroutine cm_sync_H_index

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_create_selatoms
  !> @brief        allocate domain_constraint%selatoms
  !! @authors      NT
  !! @param[in]    selatoms : select atom information
  !! @param[in]    lexist   : flag for cache file is exist or not
  !! @param[in]    sel_size : size of selatoms (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_create_selatoms(selatoms, lexist, sel_size)

    ! formal arguments
    type(s_selatoms),        intent(in)  :: selatoms
    logical,                 intent(out) :: lexist
    integer,       optional, intent(in)  :: sel_size

    ! local variables
    integer                  :: iselatoms, file_out
    integer                  :: i, idata, var_size, num
    character(30)            :: name

    integer,    allocatable  :: page_selatoms(:)


    cm_num_selatoms_list = cm_num_selatoms_list + 1
    iselatoms = cm_num_selatoms_list

    if (present(sel_size)) then
      cm_num_selatoms(iselatoms) = sel_size
    else
      cm_num_selatoms(iselatoms) = size(selatoms%idx)
    end if
    var_size = cm_num_selatoms(iselatoms)


    ! allocate page memory
    !
    !$omp parallel
    allocate(cmm_selatoms(iselatoms)%v(CMPageSize, CMNumPages))
    call cm_init_pages(cmp_selatoms(:,iselatoms),  CMPageSize)
    !$omp end parallel
    !$omp barrier


    ! create cache file
    !
    write(name,'(a,i0)') 'selatoms', iselatoms
    inquire(file=cm_get_cache_name(1, name), exist=lexist)
    if (lexist) then
      write(MsgOut,'(a)') 'Cm_Create_Selatoms> Re-use cache file. '//trim(name)
      return
    end if

    allocate(page_selatoms(CMPageSize))

    idata = 1

    do i = 1, var_size

      page_selatoms(idata) = selatoms%idx(i)
      idata = idata + 1
      
      if (idata > CMPageSize .or. i + 1 > var_size) then

        if (idata - 1 == CMPageSize) then
          num = i - CMPageSize + 1
        else
          num = i / CMPageSize * CMPageSize + 1
        end if

        call open_binary_file(file_out, &
                              cm_get_cache_name(num, name), &
                              IOFileOutputNew)
        write(file_out) page_selatoms(1:CMPageSize)
        call close_file(file_out)

        idata = 1

      end if
    end do

    deallocate(page_selatoms)

    return

  end subroutine cm_create_selatoms

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_delete_selatoms
  !> @brief        deallocate domain_constraint%selatoms
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_delete_selatoms

    ! local variables
    integer        :: iselatoms


    do iselatoms = 1, cm_num_selatoms_list

      if (.not. allocated(cmm_selatoms(iselatoms)%v)) &
        exit

      ! deallocate page memory
      !
      !$omp parallel
      deallocate(cmm_selatoms(iselatoms)%v)
      !$omp end parallel

    end do

    cm_num_selatoms_list = 0

    return

  end subroutine cm_delete_selatoms

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_selatoms
  !> @brief        get selatoms
  !! @authors      NT
  !! @param[in]    i         : index of atom selection
  !! @param[in]    iselatoms :  index of selection
  !! @return       index of atom
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function cm_selatoms(i, iselatoms)

    ! return value
    integer                  :: cm_selatoms

    ! formal arguments
    integer,                 intent(in) :: i
    integer,                 intent(in) :: iselatoms

    ! local variables
    integer                  :: mem_no, mem_index, pin_offset
    character(30)            :: name
    logical                  :: page_in


    call cm_get_page(cmp_selatoms(:,iselatoms), i, mem_no, mem_index, &
                     page_in, pin_offset)

    if (page_in) then
      write(name,'(a,i0)') 'selatoms', iselatoms
      call cm_page_in_I1(pin_offset, name, cmm_selatoms(iselatoms)%v(:, mem_no))
    end if

    cm_selatoms = cmm_selatoms(iselatoms)%v(mem_index, mem_no)

    return

  end function cm_selatoms

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_backup_cache
  !> @brief        backup cache file
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_backup_cache(var_name)

    ! formal arguments
    character(*),            intent(in)    :: var_name

    ! local varables
    integer                  :: offset
    character(MaxFilename)   :: srcname, dstname
    logical                  :: lexist


    offset = 1
    dstname = cm_get_cache_name(offset, trim(var_name)//'.org')

    inquire(file=dstname, exist=lexist)
    if (lexist) &
      return

    do while(.true.)

      srcname = cm_get_cache_name(offset, var_name)
      dstname = cm_get_cache_name(offset, trim(var_name)//'.org')

      inquire(file=srcname, exist=lexist)
      if (lexist) then
        call system('cp '//trim(srcname)//' '//trim(dstname))
      else
        exit
      end if

      offset = offset + CMPageSize
    end do

    return

  end subroutine cm_backup_cache

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_recall_cache
  !> @brief        recall cache file
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_recall_cache(var_name)

    ! formal arguments
    character(*),            intent(in)    :: var_name

    ! local varables
    integer                  :: offset
    character(Maxfilename)   :: srcname, dstname
    logical                  :: lexist


    offset = 1
    dstname = cm_get_cache_name(offset, trim(var_name)//'.org')

    inquire(file=dstname, exist=lexist)
    if (.not. lexist) &
      return

    do while(.true.)

      srcname = cm_get_cache_name(offset, trim(var_name)//'.org')
      dstname = cm_get_cache_name(offset, var_name)

      inquire(file=srcname, exist=lexist)
      if (lexist) then
        call system('mv '//trim(srcname)//' '//trim(dstname))
      else
        exit
      end if

      offset = offset + CMPageSize
    end do

    return

  end subroutine cm_recall_cache

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_page_in_C1
  !> @brief        page-in helper
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_page_in_C1(offset, var_name, var_array)

    ! formal arguments
    integer,                 intent(in)    :: offset
    character(*),            intent(in)    :: var_name
    character(*),            intent(inout) :: var_array(:)

    ! local variables
    integer                  :: file

    
    !$omp critical
    call open_binary_file(file, &
         cm_get_cache_name(offset, var_name), IOFileInput)
    !$omp end critical

    read(file) var_array(1:CMPageSize)

    !$omp critical
    call close_file(file)
    !$omp end critical

    return

  end subroutine cm_page_in_C1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_page_in_L1
  !> @brief        page-in helper
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_page_in_L1(offset, var_name, var_array)

    ! formal arguments
    integer,                 intent(in)    :: offset
    character(*),            intent(in)    :: var_name
    logical,                 intent(inout) :: var_array(:)

    ! local variables
    integer                  :: file

    
    !$omp critical
    call open_binary_file(file, &
         cm_get_cache_name(offset, var_name), IOFileInput)
    !$omp end critical

    read(file) var_array(1:CMPageSize)

    !$omp critical
    call close_file(file)
    !$omp end critical

    return

  end subroutine cm_page_in_L1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_page_in_I1
  !> @brief        page-in helper
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_page_in_I1(offset, var_name, var_array)

    ! formal arguments
    integer,                 intent(in)    :: offset
    character(*),            intent(in)    :: var_name
    integer,                 intent(inout) :: var_array(:)

    ! local variables
    integer                  :: file

    
    !$omp critical
    call open_binary_file(file, &
         cm_get_cache_name(offset, var_name), IOFileInput)
    !$omp end critical

    read(file) var_array(1:CMPageSize)

    !$omp critical
    call close_file(file)
    !$omp end critical

    return

  end subroutine cm_page_in_I1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_page_in_IN
  !> @brief        page-in helper
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_page_in_IN(n, offset, var_name, var_array)

    ! formal arguments
    integer,                 intent(in)    :: n
    integer,                 intent(in)    :: offset
    character(*),            intent(in)    :: var_name
    integer,                 intent(inout) :: var_array(:,:)

    ! local variables
    integer                  :: file

    
    !$omp critical
    call open_binary_file(file, &
         cm_get_cache_name(offset, var_name), IOFileInput)
    !$omp end critical

    read(file) var_array(1:n,1:CMPageSize)

    !$omp critical
    call close_file(file)
    !$omp end critical

    return

  end subroutine cm_page_in_IN

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_page_in_R1
  !> @brief        page-in helper
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_page_in_R1(offset, var_name, var_array)

    ! formal arguments
    integer,                 intent(in)    :: offset
    character(*),            intent(in)    :: var_name
    real(wp),                intent(inout) :: var_array(:)

    ! local variables
    integer                  :: file

    
    !$omp critical
    call open_binary_file(file, &
         cm_get_cache_name(offset, var_name), IOFileInput)
    !$omp end critical

    read(file) var_array(1:CMPageSize)

    !$omp critical
    call close_file(file)
    !$omp end critical

    return

  end subroutine cm_page_in_R1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_page_in_R3
  !> @brief        page-in helper
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_page_in_R3(offset, var_name, var_array)

    ! formal arguments
    integer,                 intent(in)    :: offset
    character(*),            intent(in)    :: var_name
    real(wp),                intent(inout) :: var_array(:,:)

    ! local variables
    integer                  :: file

    
    !$omp critical
    call open_binary_file(file, &
         cm_get_cache_name(offset, var_name), IOFileInput)
    !$omp end critical

    read(file) var_array(1:3,1:CMPageSize)

    !$omp critical
    call close_file(file)
    !$omp end critical

    return

  end subroutine cm_page_in_R3

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_page_out_C1
  !> @brief        page-out helper (for integer array)
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_page_out_C1(offset, var_name, var_array)

    ! formal arguments
    integer,                 intent(in)    :: offset
    character(*),            intent(in)    :: var_name
    character(*),            intent(inout) :: var_array(:)

    ! local variables
    integer                  :: file
    

    !$omp critical
    call open_binary_file(file, &
         cm_get_cache_name(offset, var_name), IOFileOutputReplace)
    !$omp end critical

    write(file) var_array(1:CMPageSize)

    !$omp critical
    call close_file(file)
    !$omp end critical

    return

  end subroutine cm_page_out_C1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_page_out_L1
  !> @brief        page-out helper (for integer array)
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_page_out_L1(offset, var_name, var_array)

    ! formal arguments
    integer,                 intent(in)    :: offset
    character(*),            intent(in)    :: var_name
    logical,                 intent(inout) :: var_array(:)

    ! local variables
    integer                  :: file
    

    !$omp critical
    call open_binary_file(file, &
         cm_get_cache_name(offset, var_name), IOFileOutputReplace)
    !$omp end critical

    write(file) var_array(1:CMPageSize)

    !$omp critical
    call close_file(file)
    !$omp end critical

    return

  end subroutine cm_page_out_L1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_page_out_I1
  !> @brief        page-out helper (for integer array)
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_page_out_I1(offset, var_name, var_array)

    ! formal arguments
    integer,                 intent(in)    :: offset
    character(*),            intent(in)    :: var_name
    integer,                 intent(inout) :: var_array(:)

    ! local variables
    integer                  :: file

    
    !$omp critical
    call open_binary_file(file, &
         cm_get_cache_name(offset, var_name), IOFileOutputReplace)
    !$omp end critical

    write(file) var_array(1:CMPageSize)

    !$omp critical
    call close_file(file)
    !$omp end critical

    return

  end subroutine cm_page_out_I1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_page_out_IN
  !> @brief        page-out helper (for integer array)
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_page_out_IN(n, offset, var_name, var_array)

    ! formal arguments
    integer,                 intent(in)    :: n
    integer,                 intent(in)    :: offset
    character(*),            intent(in)    :: var_name
    integer,                 intent(inout) :: var_array(:,:)

    ! local variables
    integer                  :: file
    
    !$omp critical
    call open_binary_file(file, &
         cm_get_cache_name(offset, var_name), IOFileOutputReplace)
    !$omp end critical

    write(file) var_array(1:n,1:CMPageSize)

    !$omp critical
    call close_file(file)
    !$omp end critical

    return

  end subroutine cm_page_out_IN

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_page_out_R1
  !> @brief        page-out helper (for integer array)
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_page_out_R1(offset, var_name, var_array)

    ! formal arguments
    integer,                 intent(in)    :: offset
    character(*),            intent(in)    :: var_name
    real(wp),                intent(inout) :: var_array(:)

    ! local variables
    integer                  :: file
    

    !$omp critical
    call open_binary_file(file, &
         cm_get_cache_name(offset, var_name), IOFileOutputReplace)
    !$omp end critical

    write(file) var_array(1:CMPageSize)

    !$omp critical
    call close_file(file)
    !$omp end critical

    return

  end subroutine cm_page_out_R1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_page_out_R3
  !> @brief        page-out helper (for integer array)
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_page_out_R3(offset, var_name, var_array)

    ! formal arguments
    integer,                 intent(in)    :: offset
    character(*),            intent(in)    :: var_name
    real(wp),                intent(inout) :: var_array(:,:)

    ! local variables
    integer                  :: file
    

    !$omp critical
    call open_binary_file(file, &
         cm_get_cache_name(offset, var_name), IOFileOutputReplace)
    !$omp end critical

    write(file) var_array(1:3,1:CMPageSize)

    !$omp critical
    call close_file(file)
    !$omp end critical

    return

  end subroutine cm_page_out_R3

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_init_pages
  !> @brief        initialize page information
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_init_pages(pages, size)

    ! formal arguments
    type(s_cm_page), target, intent(inout) :: pages(:)
    integer,                 intent(in)    :: size

    ! local variables
    integer                  :: i


    do i = 1, CMNumPages + 1

      pages(i)%mem_no = i - 1
      pages(i)%offset = 0
      pages(i)%size   = size
      pages(i)%dirty  = .false.
      if (i == CMNumPages + 1) then
        pages(i)%next => null()
      else
        pages(i)%next => pages(i+1)
      end if

    end do

    return

  end subroutine cm_init_pages

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_get_page
  !> @brief        get page information
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_get_page(pages, index, &
                         mem_no,       &
                         mem_index,    &
                         page_in,      &
                         pin_offset,   &
                         page_out,     &
                         pout_mem_no,  &
                         pout_offset)

    ! formal arguments
    type(s_cm_page), target, intent(inout) :: pages(:)
    integer,                 intent(in)    :: index
    integer,                 intent(out)   :: mem_no
    integer,                 intent(out)   :: mem_index
    logical,                 intent(out)   :: page_in
    integer,                 intent(out)   :: pin_offset
    logical,       optional, intent(out)   :: page_out
    integer,       optional, intent(out)   :: pout_mem_no
    integer,       optional, intent(out)   :: pout_offset

    ! local variables
    type(s_cm_page), pointer :: page, page_prv


    page_in  = .false.

    if (present(page_out)) &
      page_out = .false.

    if (index == 0) &
      call error_msg('Cm_Get_Page> index is ZERO.')

    ! check head
    !
    page => pages(1)%next

    ! empty page
    if (page%offset == 0) then

      page%offset = ((index - 1) / page%size) * page%size + 1

      page_in     = .true.
      pin_offset  = page%offset

      mem_no      = page%mem_no
      mem_index   = index - page%offset + 1

      cm_fast     = cm_fast + 1
      cm_page_in  = cm_page_in + 1
      return

    ! matched page
    else if (page%offset <= index .and. index < page%offset + page%size) then

      mem_no      = page%mem_no
      mem_index   = index - page%offset + 1

      cm_fast     = cm_fast + 1
      return

    end if


    ! check tails
    !

    do while(.true.)

      page_prv => page
      page     => page%next

      ! empty page
      if (page%offset == 0) then

        page%offset = ((index - 1) / page%size) * page%size + 1

        page_in     = .true.
        pin_offset  = page%offset

        mem_no      = page%mem_no
        mem_index   = index - page%offset + 1

        page_prv%next => page%next

        page%next     => pages(1)%next
        pages(1)%next => page

        cm_slow     = cm_slow + 1
        cm_page_in  = cm_page_in + 1
        return

      ! matched page
      else if (page%offset <= index .and. index < page%offset + page%size) then

        mem_no    = page%mem_no
        mem_index = index - page%offset + 1

        page_prv%next => page%next

        page%next     => pages(1)%next
        pages(1)%next => page

        cm_slow = cm_slow + 1
        return

      end if

      if (.not. associated(page%next)) &
        exit

    end do


    ! check tip
    !

    if (present(page_out)) then

      if (page%dirty) then
        page_out    = .true.
        pout_mem_no = page%mem_no
        pout_offset = page%offset
        page%dirty  = .false.
        
        cm_page_out = cm_page_out + 1

      end if

    end if

    page%offset = ((index - 1) / page%size) * page%size + 1

    page_in    = .true.
    pin_offset = page%offset

    mem_no     = page%mem_no
    mem_index  = index - page%offset + 1

    page_prv%next => page%next

    page%next     => pages(1)%next
    pages(1)%next => page

    cm_slow    = cm_slow + 1
    cm_page_in = cm_page_in + 1

    return

  end subroutine cm_get_page

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_dirty_page
  !> @brief        mark dirty to head page 
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_dirty_page(pages)

    ! formal arguments
    type(s_cm_page), target, intent(inout) :: pages(:)


    pages(1)%next%dirty = .true.

    return

  end subroutine cm_dirty_page

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      cm_get_cache_name
  !> @brief        get cache filename
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function cm_get_cache_name(num, type_str)

    ! return values
    character(MaxFilename)        :: cm_get_cache_name

    ! formal arguments
    integer,      intent(in) :: num
    character(*), intent(in) :: type_str


    write(cm_get_cache_name, '(a,a,a,a,i10.10,a,a)') &
         trim(cm_cache_dir),  '/', &
         trim(cm_cache_name), '_', &
         num, '_', trim(type_str)

    return

  end function cm_get_cache_name

end module pr_cached_molecule_mod
