!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pr_huge_molecule_mod
!> @brief   huge molecule information (memory version)
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

module pr_huge_molecule_mod

  use select_atoms_str_mod
  use fileio_mod
  use string_mod
  use messages_mod
  use constants_mod
  use mpi_parallel_mod
  use mpi

  implicit none
  private

  ! public
  !

  ! parameters
  integer,        public, parameter :: HMMaxSelAtomsList = 10

  ! variables
  character(MaxFilename), public :: hm_work_dir  = './'
  character(MaxFilename), public :: hm_work_name = ''

  integer, public :: my_prst_rank = 0
  integer, public :: nproc_prst   = 0
  logical, public :: hm_atom_ref_set_flag = .false.

  ! s_molecule
  integer,        public :: hm_num_atoms           = 0
  integer,        public :: hm_num_bonds           = 0
  integer,        public :: hm_num_angles          = 0
  integer,        public :: hm_num_dihedrals       = 0
  integer,        public :: hm_num_impropers       = 0
  integer,        public :: hm_num_cmaps           = 0
  ! s_enefunc%table
  integer,        public :: hm_num_water_list      = 0
  integer,        public :: hm_num_solute_list     = 0
  ! s_selatoms
  integer,        public :: hm_num_selatoms_list   = 0
  integer,        public :: hm_num_selatoms(1:HMMaxSelAtomsList) = 0

  ! subroutines
  public :: hm_initialize
  public :: hm_finalize

  ! s_molecule
  public :: hm_alloc_atoms
  public :: hm_atom_no
  public :: hm_atom_cls_no
  public :: hm_atom_cls_name
  public :: hm_atom_name
  public :: hm_light_atom_mass
  public :: hm_light_atom_name
  public :: hm_atom_coord
  public :: hm_atom_refcoord
  public :: hm_atom_velocity
  public :: hm_charge
  public :: hm_mass
  public :: hm_residue_no
  public :: hm_residue_name
  public :: hm_segment_name
  public :: hm_molecule_no
  public :: hm_set_atom_no
  public :: hm_set_atom_cls_no
  public :: hm_set_atom_cls_name
  public :: hm_set_atom_name
  public :: hm_set_light_atom_mass
  public :: hm_set_light_atom_name
  public :: hm_set_atom_coord
  public :: hm_set_atom_refcoord
  public :: hm_set_atom_velocity
  public :: hm_set_charge
  public :: hm_set_mass
  public :: hm_set_residue_no
  public :: hm_set_residue_name
  public :: hm_set_segment_name
  public :: hm_set_molecule_no

  public :: hm_alloc_bonds
  public :: hm_bond_list
  public :: hm_set_bond_list

  public :: hm_alloc_angls
  public :: hm_angl_list
  public :: hm_set_angl_list

  public :: hm_alloc_dihes
  public :: hm_dihe_list
  public :: hm_set_dihe_list

  public :: hm_alloc_imprs
  public :: hm_impr_list
  public :: hm_set_impr_list

  public :: hm_alloc_cmaps
  public :: hm_cmap_list
  public :: hm_set_cmap_list

  ! s_enefunc%table
  public :: hm_alloc_water_list
  public :: hm_water_list
  public :: hm_set_water_list

  public :: hm_alloc_solute_list
  public :: hm_solute_list
  public :: hm_set_solute_list

  public :: hm_alloc_solute_list_inv
  public :: hm_solute_list_inv
  public :: hm_set_solute_list_inv

  public :: hm_alloc_solute
  public :: hm_solute
  public :: hm_set_solute

  ! s_constraints
  public :: hm_alloc_duplicate
  public :: hm_duplicate
  public :: hm_set_duplicate

  public :: hm_alloc_H_index
  public :: hm_H_index
  public :: hm_set_H_index

  public :: hm_alloc_ring
  public :: hm_set_ring
  public :: hm_ring
  public :: hm_set_num_connect
  public :: hm_num_connect
  public :: hm_set_connectivity
  public :: hm_connectivity

  ! s_selatoms
  public :: hm_alloc_selatoms
  public :: hm_selatoms
  public :: hm_set_selatoms


  ! private
  !

  type s_hm_selatoms
    integer,    allocatable :: v(:)
  end type s_hm_selatoms

  integer,      allocatable :: hmm_atom_no        (:)
  integer,      allocatable :: hmm_atom_cls_no    (:)
  character(6), allocatable :: hmm_atom_cls_name  (:)
  character(4), allocatable :: hmm_atom_name      (:)
  integer,      allocatable :: hmm_light_atom_mass(:)
  integer,      allocatable :: hmm_light_atom_name(:)
  real(wip),    allocatable :: hmm_atom_coord     (:,:)
  real(wip),    allocatable :: hmm_atom_refcoord  (:,:)
  real(wip),    allocatable :: hmm_atom_velocity  (:,:)
  real(wp),     allocatable :: hmm_charge         (:)
  real(wip),    allocatable :: hmm_mass           (:)
  integer,      allocatable :: hmm_residue_no     (:)
  character(6), allocatable :: hmm_residue_name   (:)
  character(4), allocatable :: hmm_segment_name   (:)
  integer,      allocatable :: hmm_molecule_no    (:)
  integer,      allocatable :: hmm_bond_list      (:,:)
  integer,      allocatable :: hmm_angl_list      (:,:)
  integer,      allocatable :: hmm_dihe_list      (:,:)
  integer,      allocatable :: hmm_impr_list      (:,:)
  integer,      allocatable :: hmm_cmap_list      (:,:)
  integer,      allocatable :: hmm_water_list     (:,:)
  integer,      allocatable :: hmm_solute_list    (:)
  integer,      allocatable :: hmm_solute_list_inv(:)
  logical,      allocatable :: hmm_solute         (:)
  integer,      allocatable :: hmm_duplicate      (:)
  integer,      allocatable :: hmm_H_index        (:,:)
  integer,      allocatable :: hmm_ring           (:)
  integer,      allocatable :: hmm_num_connect    (:,:)
  integer,      allocatable :: hmm_connectivity   (:,:,:)
  type(s_hm_selatoms)       :: hmm_selatoms(HMMaxSelAtomsList)

  ! subroutines
  private :: hm_alloc_error
  private :: hm_check_memory

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_initialize
  !> @brief        initialize variables
  !! @authors      NT
  !! @param[in]    work_path : path for working files
  !! @param[in]    work_name : working file prefix
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_initialize(work_path, work_name)

    ! formal arguments
    character(*),            intent(in)    :: work_path
    character(*),            intent(in)    :: work_name


    hm_work_dir  = work_path
    hm_work_name = work_name

    if (main_rank) then

      write(MsgOut,'(a)') 'Hm_Initialize> '
      write(MsgOut,'(a)')   ' '

    end if

    return

  end subroutine hm_initialize

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_finalize
  !> @brief        delete working files and print performance information
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_finalize

    ! local variables
    integer                  :: i

    if (main_rank) then

      write(MsgOut,'(a)') 'Hm_Finalize> cleanup memory.'
      write(MsgOut,'(a)') ' '

    end if

    if (allocated(hmm_atom_no))         deallocate(hmm_atom_no)
    if (allocated(hmm_atom_cls_no))     deallocate(hmm_atom_cls_no)
    if (allocated(hmm_atom_cls_name))   deallocate(hmm_atom_cls_name)
    if (allocated(hmm_atom_name))       deallocate(hmm_atom_name)
    if (allocated(hmm_light_atom_mass)) deallocate(hmm_light_atom_mass)
    if (allocated(hmm_light_atom_name)) deallocate(hmm_light_atom_name)
    if (allocated(hmm_atom_coord))      deallocate(hmm_atom_coord)
    if (allocated(hmm_atom_refcoord))   deallocate(hmm_atom_refcoord)
    if (allocated(hmm_atom_velocity))   deallocate(hmm_atom_velocity)
    if (allocated(hmm_charge))          deallocate(hmm_charge)
    if (allocated(hmm_mass))            deallocate(hmm_mass)
    if (allocated(hmm_residue_no))      deallocate(hmm_residue_no)
    if (allocated(hmm_residue_name))    deallocate(hmm_residue_name)
    if (allocated(hmm_segment_name))    deallocate(hmm_segment_name)
    if (allocated(hmm_molecule_no))     deallocate(hmm_molecule_no)
    if (allocated(hmm_bond_list))       deallocate(hmm_bond_list)
    if (allocated(hmm_angl_list))       deallocate(hmm_angl_list)
    if (allocated(hmm_dihe_list))       deallocate(hmm_dihe_list)
    if (allocated(hmm_impr_list))       deallocate(hmm_impr_list)
    if (allocated(hmm_cmap_list))       deallocate(hmm_cmap_list)
    if (allocated(hmm_water_list))      deallocate(hmm_water_list)
    if (allocated(hmm_solute_list))     deallocate(hmm_solute_list)
    if (allocated(hmm_solute_list_inv)) deallocate(hmm_solute_list_inv)
    if (allocated(hmm_solute))          deallocate(hmm_solute)
    if (allocated(hmm_duplicate))       deallocate(hmm_duplicate)
    if (allocated(hmm_H_index))         deallocate(hmm_H_index)
    do i = 1, HMMaxSelAtomsList
      if (allocated(hmm_selatoms(i)%v)) deallocate(hmm_selatoms(i)%v)
    end do

    return

  end subroutine hm_finalize

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_alloc_atoms
  !> @brief        allocate atom informations
  !! @authors      NT
  !! @param[in]    num_atoms : number of atoms
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_alloc_atoms(num_atoms)

    ! formal arguments
    integer,                 intent(in)    :: num_atoms

    ! local variables
    integer                  :: stat, i, ierr


    hm_num_atoms = num_atoms

    call hm_check_memory

    allocate(hmm_atom_no        (num_atoms), stat=stat)
    if (stat /= 0) call hm_alloc_error('atom_no')
    allocate(hmm_atom_cls_no    (num_atoms), stat=stat)
    if (stat /= 0) call hm_alloc_error('atom_cls_no')
    allocate(hmm_atom_cls_name  (num_atoms), stat=stat)
    if (stat /= 0) call hm_alloc_error('atom_cls_name')
    allocate(hmm_atom_name      (num_atoms), stat=stat)
    if (stat /= 0) call hm_alloc_error('atom_name')
    allocate(hmm_light_atom_mass(num_atoms), stat=stat)
    if (stat /= 0) call hm_alloc_error('light_atom_mass')
    allocate(hmm_light_atom_name(num_atoms), stat=stat)
    if (stat /= 0) call hm_alloc_error('light_atom_name')
    allocate(hmm_atom_coord     (3,num_atoms), stat=stat)
    if (stat /= 0) call hm_alloc_error('atom_coord')
    allocate(hmm_atom_refcoord  (3,num_atoms), stat=stat)
    if (stat /= 0) call hm_alloc_error('atom_refcoord')
    allocate(hmm_atom_velocity  (3,num_atoms), stat=stat)
    if (stat /= 0) call hm_alloc_error('atom_velocity')
    allocate(hmm_charge         (num_atoms), stat=stat)
    if (stat /= 0) call hm_alloc_error('charge')
    allocate(hmm_mass           (num_atoms), stat=stat)
    if (stat /= 0) call hm_alloc_error('mass')
    allocate(hmm_residue_no     (num_atoms), stat=stat)
    if (stat /= 0) call hm_alloc_error('residue_no')
    allocate(hmm_residue_name   (num_atoms), stat=stat)
    if (stat /= 0) call hm_alloc_error('residue_name')
    allocate(hmm_segment_name   (num_atoms), stat=stat)
    if (stat /= 0) call hm_alloc_error('segment_name')
    allocate(hmm_molecule_no    (num_atoms), stat=stat)
    if (stat /= 0) call hm_alloc_error('molecule_no')

    do i = 1, num_atoms
      hmm_atom_no        (i) = 0
      hmm_atom_cls_no    (i) = 0
      hmm_atom_cls_name  (i) = ''
      hmm_atom_name      (i) = ''
      hmm_light_atom_mass(i) = 0
      hmm_light_atom_name(i) = 0
      hmm_atom_coord     (1:3,i) = 0.0_wip
      hmm_atom_refcoord  (1:3,i) = 0.0_wip
      hmm_atom_velocity  (1:3,i) = 0.0_wip
      hmm_charge         (i) = 0.0_wp
      hmm_mass           (i) = 0.0_wip
      hmm_residue_no     (i) = 0
      hmm_residue_name   (i) = ''
      hmm_segment_name   (i) = ''
      hmm_molecule_no    (i) = 0
    end do

    return

  end subroutine hm_alloc_atoms

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_atom_no
  !> @brief        get molecule%atom_no
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @return       atom serial number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_atom_no(i)

    ! return value
    integer                  :: hm_atom_no

    ! formal arguments
    integer,                 intent(in)    :: i


    hm_atom_no = hmm_atom_no(i)

    return

  end function hm_atom_no

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_atom_cls_no
  !> @brief        get molecule%atom_cls_no
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @return       atom class number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_atom_cls_no(i)

    ! return value
    integer                  :: hm_atom_cls_no

    ! formal arguments
    integer,                 intent(in)    :: i


    hm_atom_cls_no = hmm_atom_cls_no(i)

    return

  end function hm_atom_cls_no

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_atom_cls_name
  !> @brief        get molecule%atom_cls_name
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @return       atom class name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_atom_cls_name(i)

    ! return value
    character(6)             :: hm_atom_cls_name

    ! formal arguments
    integer,                 intent(in)    :: i


    hm_atom_cls_name = hmm_atom_cls_name(i)

    return

  end function hm_atom_cls_name

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_atom_name
  !> @brief        get molecule%atom_name
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @return       atom name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_atom_name(i)

    ! return value
    character(4)             :: hm_atom_name

    ! formal arguments
    integer,                 intent(in)    :: i


    hm_atom_name = hmm_atom_name(i)

    return

  end function hm_atom_name

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_light_atom_mass
  !> @brief        get molecule%light_atom_mass
  !! @authors      JJ
  !! @param[in]    i : index of atom
  !! @return       atom name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_light_atom_mass(i)

    ! return value
    integer                  :: hm_light_atom_mass

    ! formal arguments
    integer,                 intent(in)    :: i


    hm_light_atom_mass = hmm_light_atom_mass(i)

    return

  end function hm_light_atom_mass

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_light_atom_name
  !> @brief        get molecule%light_atom_name
  !! @authors      JJ
  !! @param[in]    i : index of atom
  !! @return       atom name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_light_atom_name(i)

    ! return value
    integer                  :: hm_light_atom_name

    ! formal arguments
    integer,                 intent(in)    :: i


    hm_light_atom_name = hmm_light_atom_name(i)

    return

  end function hm_light_atom_name

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_atom_coord
  !> @brief        get molecule%atom_coord
  !! @authors      NT
  !! @param[in]    j : index of axis
  !! @param[in]    i : index of atom
  !! @return       atom coordinate
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_atom_coord(j, i)

    ! return value
    real(wip)                :: hm_atom_coord

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i


    hm_atom_coord = hmm_atom_coord(j, i)

    return

  end function hm_atom_coord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_atom_refcoord
  !> @brief        get molecule%atom_refcoord
  !! @authors      NT
  !! @param[in]    j : index of axis
  !! @param[in]    i : index of atom
  !! @return       atom reference coordinate
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_atom_refcoord(j, i)

    ! return value
    real(wip)                :: hm_atom_refcoord

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i


    hm_atom_refcoord = hmm_atom_refcoord(j, i)

    return

  end function hm_atom_refcoord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_atom_velocity
  !> @brief        get molecule%atom_velocity
  !! @authors      NT
  !! @param[in]    j : index of axis
  !! @param[in]    i : index of atom
  !! @return       atom velocity
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_atom_velocity(j, i)

    ! return value
    real(wip)                :: hm_atom_velocity

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i


    hm_atom_velocity = hmm_atom_velocity(j, i)

    return

  end function hm_atom_velocity

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_charge
  !> @brief        get molecule%charge
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @return       atom charge
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_charge(i)

    ! return value
    real(wp)                 :: hm_charge

    ! formal arguments
    integer,                 intent(in)    :: i


    hm_charge = hmm_charge(i)

    return

  end function hm_charge

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_mass
  !> @brief        get molecule%mass
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @return       atom mass
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_mass(i)

    ! return value
    real(wip)                :: hm_mass

    ! formal arguments
    integer,                 intent(in)    :: i


    hm_mass = hmm_mass(i)

    return

  end function hm_mass

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_residue_no
  !> @brief        get molecule%residue_no
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @return       residue number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_residue_no(i)

    ! return value
    integer                  :: hm_residue_no

    ! formal arguments
    integer,                 intent(in)    :: i


    hm_residue_no = hmm_residue_no(i)

    return

  end function hm_residue_no

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_residue_name
  !> @brief        get molecule%residue_name
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @return       residue name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_residue_name(i)

    ! return value
    character(6)             :: hm_residue_name

    ! formal arguments
    integer,                 intent(in)    :: i


    hm_residue_name = hmm_residue_name(i)

    return

  end function hm_residue_name

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_segment_name
  !> @brief        get molecule%segment_name
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @return       segment name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_segment_name(i)

    ! return value
    character(4)             :: hm_segment_name

    ! formal arguments
    integer,                 intent(in)    :: i


    hm_segment_name = hmm_segment_name(i)

    return

  end function hm_segment_name

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_molecule_no
  !> @brief        get molecule%molecule_no
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @return       molecule number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_molecule_no(i)

    ! return value
    integer                  :: hm_molecule_no

    ! formal arguments
    integer,                 intent(in)    :: i


    hm_molecule_no = hmm_molecule_no(i)

    return

  end function hm_molecule_no

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_atom_no
  !> @brief        set molecule%atom_no
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @param[in]    no : atom serial number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_atom_no(i, no)

    ! formal arguments
    integer,                 intent(in)    :: i
    integer,                 intent(in)    :: no


    hmm_atom_no(i) = no

    return

  end subroutine hm_set_atom_no

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_atom_cls_no
  !> @brief        set molecule%atom_cls_no
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @param[in]    no : atom class number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_atom_cls_no(i, no)

    ! formal arguments
    integer,                 intent(in)    :: i
    integer,                 intent(in)    :: no


    hmm_atom_cls_no(i) = no

    return

  end subroutine hm_set_atom_cls_no

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_atom_cls_name
  !> @brief        set molecule%atom_cls_name
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @param[in]    name : atom class name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_atom_cls_name(i, name)

    ! formal arguments
    integer,                 intent(in)    :: i
    character(*),            intent(in)    :: name


    hmm_atom_cls_name(i) = name

    return

  end subroutine hm_set_atom_cls_name

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_atom_name
  !> @brief        set molecule%atom_name
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @param[in]    name : atom name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_atom_name(i, name)

    ! formal arguments
    integer,                 intent(in)    :: i
    character(*),            intent(in)    :: name

    hmm_atom_name(i) = name

    return

  end subroutine hm_set_atom_name

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_light_atom_mass
  !> @brief        set molecule%light_atom_mass
  !! @authors      JJ
  !! @param[in]    i : index of atom
  !! @param[in]    hlm : key of light_atom_mass
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_light_atom_mass(i, hlm)

    ! formal arguments
    integer,                 intent(in)    :: i
    integer,                 intent(in)    :: hlm

    hmm_light_atom_mass(i) = hlm

    return

  end subroutine hm_set_light_atom_mass

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_light_atom_name
  !> @brief        set molecule%light_atom_name
  !! @authors      JJ
  !! @param[in]    i : index of atom
  !! @param[in]    hln : key of light_atom_name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_light_atom_name(i, hln)

    ! formal arguments
    integer,                 intent(in)    :: i
    integer,                 intent(in)    :: hln

    hmm_light_atom_name(i) = hln

    return

  end subroutine hm_set_light_atom_name

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_atom_coord
  !> @brief        set molecule%atom_coord
  !! @authors      NT
  !! @param[in]    j : index of axis
  !! @param[in]    i : index of atom
  !! @param[in]    crd : coordinate value
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_atom_coord(j, i, crd)

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i
    real(wip),               intent(in)    :: crd


    hmm_atom_coord(j, i) = crd

    return

  end subroutine hm_set_atom_coord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_atom_refcoord
  !> @brief        set molecule%atom_refcoord
  !! @authors      NT
  !! @param[in]    j : index of axis
  !! @param[in]    i : index of atom
  !! @param[in]    crd : coordinate value
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_atom_refcoord(j, i, crd)

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i
    real(wip),               intent(in)    :: crd


    hmm_atom_refcoord(j, i) = crd
    if (.not. hm_atom_ref_set_flag) hm_atom_ref_set_flag=.true.

    return

  end subroutine hm_set_atom_refcoord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_atom_velocity
  !> @brief        set molecule%atom_velocity
  !! @authors      NT
  !! @param[in]    j : index of axis
  !! @param[in]    i : index of atom
  !! @param[in]    crd : coordinate value
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_atom_velocity(j, i, crd)

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i
    real(wip),               intent(in)    :: crd


    hmm_atom_velocity(j, i) = crd

    return

  end subroutine hm_set_atom_velocity

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_charge
  !> @brief        set molecule%charge
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @param[in]    charge : charge value
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_charge(i, charge)

    ! formal arguments
    integer,                 intent(in)    :: i
    real(wp),                intent(in)    :: charge


    hmm_charge(i) = charge

    return

  end subroutine hm_set_charge

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_mass
  !> @brief        set molecule%mass
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @param[in]    mass : mass value
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_mass(i, mass)

    ! formal arguments
    integer,                 intent(in)    :: i
    real(wip),               intent(in)    :: mass


    hmm_mass(i) = mass

    return

  end subroutine hm_set_mass

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_residue_no
  !> @brief        set molecule%residue_no
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @param[in]    no : residue number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_residue_no(i, no)

    ! formal arguments
    integer,                 intent(in)    :: i
    integer,                 intent(in)    :: no


    hmm_residue_no(i) = no

    return

  end subroutine hm_set_residue_no

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_residue_name
  !> @brief        set molecule%residue_name
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @param[in]    name : residue name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_residue_name(i, name)

    ! formal arguments
    integer,                 intent(in)    :: i
    character(*),            intent(in)    :: name


    hmm_residue_name(i) = name

    return

  end subroutine hm_set_residue_name

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_segment_name
  !> @brief        set molecule%segment_name
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @param[in]    name : segment name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_segment_name(i, name)

    ! formal arguments
    integer,                 intent(in)    :: i
    character(*),            intent(in)    :: name


    hmm_segment_name(i) = name

    return

  end subroutine hm_set_segment_name

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_molecule_no
  !> @brief        set molecule%molecule_no
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @param[in]    no : molecule number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_molecule_no(i, no)

    ! formal arguments
    integer,                 intent(in)    :: i
    integer,                 intent(in)    :: no


    hmm_molecule_no(i) = no

    return

  end subroutine hm_set_molecule_no

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_alloc_bonds
  !> @brief        allocate bond informations
  !! @authors      NT
  !! @param[in]    num_bonds : number of bonds
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_alloc_bonds(num_bonds)

    ! formal arguments
    integer,                 intent(in)    :: num_bonds

    ! local variables
    integer                  :: stat, i


    hm_num_bonds = num_bonds

    allocate(hmm_bond_list(2,num_bonds), stat=stat)
    if (stat /= 0) call hm_alloc_error('bond_list')

    do i = 1, num_bonds
      hmm_bond_list(1:2,i) = 0
    end do

    return

  end subroutine hm_alloc_bonds

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_bond_list
  !> @brief        get molecule%bond_list
  !! @authors      NT
  !! @param[in]    j : index of bond atom
  !! @param[in]    i : index of bond
  !! @return       index of atom
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_bond_list(j, i)

    ! return value
    integer                  :: hm_bond_list

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i


    hm_bond_list = hmm_bond_list(j, i)

    return

  end function hm_bond_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_bond_list
  !> @brief        set molecule%bond_list
  !! @authors      NT
  !! @param[in]    j : index of bond atom
  !! @param[in]    i : index of bond
  !! @param[in]    idx : index of atom
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_bond_list(j, i, idx)

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i
    integer,                 intent(in)    :: idx


    hmm_bond_list(j, i) = idx

    return

  end subroutine hm_set_bond_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_alloc_angls
  !> @brief        allocate angle informations
  !! @authors      NT
  !! @param[in]    num_angls : number of angles
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_alloc_angls(num_angls)

    ! formal arguments
    integer,                 intent(in)    :: num_angls

    ! local variables
    integer                  :: stat, i


    hm_num_angles = num_angls

    allocate(hmm_angl_list(3,num_angls), stat=stat)
    if (stat /= 0) call hm_alloc_error('angl_list')

    do i = 1, num_angls
      hmm_angl_list(1:3,i) = 0
    end do

    return

  end subroutine hm_alloc_angls

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_angl_list
  !> @brief        get molecule%angl_list
  !! @authors      NT
  !! @param[in]    j : index of angle atom
  !! @param[in]    i : index of angle
  !! @return       index of atom
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_angl_list(j, i)

    ! return value
    integer                  :: hm_angl_list

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i


    hm_angl_list = hmm_angl_list(j, i)

    return

  end function hm_angl_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_angl_list
  !> @brief        set molecule%angl_list
  !! @authors      NT
  !! @param[in]    j : index of angle atom
  !! @param[in]    i : index of angle
  !! @param[in]    idx : index of atom
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_angl_list(j, i, idx)

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i
    integer,                 intent(in)    :: idx


    hmm_angl_list(j, i) = idx

    return

  end subroutine hm_set_angl_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_alloc_dihes
  !> @brief        allocate dihedral informations
  !! @authors      NT
  !! @param[in]    num_dihes : number of dihedrals
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_alloc_dihes(num_dihes)

    ! formal arguments
    integer,                 intent(in)    :: num_dihes

    ! local variables
    integer                  :: stat, i, ierr


    hm_num_dihedrals = num_dihes

    allocate(hmm_dihe_list(4,num_dihes), stat=stat)
    if (stat /= 0) call hm_alloc_error('dihe_list')

    do i = 1, num_dihes
      hmm_dihe_list(1:4,i) = 0
    end do

    return

  end subroutine hm_alloc_dihes

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_dihe_list
  !> @brief        get molecule%dihe_list
  !! @authors      NT
  !! @param[in]    j : index of dihedral atom
  !! @param[in]    i : index of dihedral
  !! @return       index of atom
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_dihe_list(j, i)

    ! return value
    integer                  :: hm_dihe_list

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i


    hm_dihe_list = hmm_dihe_list(j, i)

    return

  end function hm_dihe_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_dihe_list
  !> @brief        set molecule%dihe_list
  !! @authors      NT
  !! @param[in]    j : index of dihedral atom
  !! @param[in]    i : index of dihedral
  !! @param[in]    idx : index of atom
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_dihe_list(j, i, idx)

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i
    integer,                 intent(in)    :: idx


    hmm_dihe_list(j, i) = idx

    return

  end subroutine hm_set_dihe_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_alloc_imprs
  !> @brief        allocate improper dihedral informations
  !! @authors      NT
  !! @param[in]    num_imprs : number of improper dihedrals
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_alloc_imprs(num_imprs)

    ! formal arguments
    integer,                 intent(in)    :: num_imprs

    ! local variables
    integer                  :: stat, i


    hm_num_impropers = num_imprs

    allocate(hmm_impr_list(4,num_imprs), stat=stat)
    if (stat /= 0) call hm_alloc_error('impr_list')

    do i = 1, num_imprs
      hmm_impr_list(1:4,i) = 0
    end do

    return

  end subroutine hm_alloc_imprs

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_impr_list
  !> @brief        get molecule%impr_list
  !! @authors      NT
  !! @param[in]    j : index of improper atom
  !! @param[in]    i : index of improper
  !! @return       index of atom
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_impr_list(j, i)

    ! return value
    integer                  :: hm_impr_list

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i


    hm_impr_list = hmm_impr_list(j, i)

    return

  end function hm_impr_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_impr_list
  !> @brief        set molecule%impr_list
  !! @authors      NT
  !! @param[in]    j : index of improper atom
  !! @param[in]    i : index of improper
  !! @param[in]    idx : index of atom
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_impr_list(j, i, idx)

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i
    integer,                 intent(in)    :: idx


    hmm_impr_list(j, i) = idx

    return

  end subroutine hm_set_impr_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_alloc_cmaps
  !> @brief        allocate cmap informations
  !! @authors      NT
  !! @param[in]    num_cmaps : number of cmaps
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_alloc_cmaps(num_cmaps)

    ! formal arguments
    integer,                 intent(in)    :: num_cmaps

    ! local variables
    integer                  :: stat, i


    hm_num_cmaps = num_cmaps

    allocate(hmm_cmap_list(8,num_cmaps), stat=stat)
    if (stat /= 0) call hm_alloc_error('cmap_list')

    do i = 1, num_cmaps
      hmm_cmap_list(1:8,i) = 0
    end do

    return

  end subroutine hm_alloc_cmaps

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_cmap_list
  !> @brief        get molecule%cmap_list
  !! @authors      NT
  !! @param[in]    j : index of cmap atom
  !! @param[in]    i : index of cmap
  !! @return       index of atom
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_cmap_list(j, i)

    ! return value
    integer                  :: hm_cmap_list

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i


    hm_cmap_list = hmm_cmap_list(j, i)

    return

  end function hm_cmap_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_cmap_list
  !> @brief        set molecule%cmap_list
  !! @authors      NT
  !! @param[in]    j : index of cmap atom
  !! @param[in]    i : index of cmap
  !! @param[in]    idx : index of atom
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_cmap_list(j, i, idx)

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i
    integer,                 intent(in)    :: idx


    hmm_cmap_list(j, i) = idx

    return

  end subroutine hm_set_cmap_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_alloc_water_list
  !> @brief        allocate water list informations
  !! @authors      NT, JJ
  !! @param[in]    num_waters : number of waters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_alloc_water_list(num_waters, na)

    ! formal arguments
    integer,                 intent(in)    :: num_waters
    integer,                 intent(in)    :: na

    ! local variables
    integer                  :: stat, i


    hm_num_water_list = num_waters

    allocate(hmm_water_list(na,num_waters), stat=stat)
    if (stat /= 0) call hm_alloc_error('water_list')

    do i = 1, num_waters
      hmm_water_list(1:na,i) = 0
    end do

    return

  end subroutine hm_alloc_water_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_water_list
  !> @brief        get enefunc%table%water_list
  !! @authors      NT
  !! @param[in]    j : index of water list atom
  !! @param[in]    i : index of water list
  !! @return       index of atom
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_water_list(j, i)

    ! return value
    integer                  :: hm_water_list

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i


    hm_water_list = hmm_water_list(j, i)

    return

  end function hm_water_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_water_list
  !> @brief        set enefunc%table%water_list
  !! @authors      NT
  !! @param[in]    j : index of water list atom
  !! @param[in]    i : index of water list
  !! @param[in]    idx : index of atom
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_water_list(j, i, idx)

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i
    integer,                 intent(in)    :: idx


    hmm_water_list(j, i) = idx

    return

  end subroutine hm_set_water_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_alloc_solute_list
  !> @brief        allocate solute list informations
  !! @authors      NT
  !! @param[in]    num_solutes : number of solutes
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_alloc_solute_list(num_solutes)

    ! formal arguments
    integer,                 intent(in)    :: num_solutes

    ! local variables
    integer                  :: stat, i


    hm_num_solute_list = num_solutes

    allocate(hmm_solute_list(num_solutes), stat=stat)
    if (stat /= 0) call hm_alloc_error('solute_list')

    do i = 1, num_solutes
      hmm_solute_list(i) = 0
    end do

    return

  end subroutine hm_alloc_solute_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_solute_list
  !> @brief        get enefunc%table%solute_list
  !! @authors      NT
  !! @param[in]    i : index of solute list
  !! @return       index of atom
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_solute_list(i)

    ! return value
    integer                  :: hm_solute_list

    ! formal arguments
    integer,                 intent(in)    :: i


    hm_solute_list = hmm_solute_list(i)

    return

  end function hm_solute_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_solute_list
  !> @brief        set enefunc%table%solute_list
  !! @authors      NT
  !! @param[in]    i : index of solute list
  !! @param[in]    idx : index of atom
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_solute_list(i, idx)

    ! formal arguments
    integer,                 intent(in)    :: i
    integer,                 intent(in)    :: idx


    hmm_solute_list(i) = idx

    return

  end subroutine hm_set_solute_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_alloc_solute_list_inv
  !> @brief        allocate solute list invert informations
  !! @authors      NT
  !! @param[in]    num_atoms : number of atoms
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_alloc_solute_list_inv(num_atoms)

    ! formal arguments
    integer,                 intent(in)    :: num_atoms

    ! local variables
    integer                  :: stat, i


    allocate(hmm_solute_list_inv(num_atoms), stat=stat)
    if (stat /= 0) call hm_alloc_error('solute_list_inv')

    do i = 1, num_atoms
      hmm_solute_list_inv(i) = 0
    end do

    return

  end subroutine hm_alloc_solute_list_inv

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_solute_list_inv
  !> @brief        get enefunc%table%solute_list_inv
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @return       index of solute list
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_solute_list_inv(i)

    ! return value
    integer                  :: hm_solute_list_inv

    ! formal arguments
    integer,                 intent(in)    :: i


    hm_solute_list_inv = hmm_solute_list_inv(i)

    return

  end function hm_solute_list_inv

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_solute_list_inv
  !> @brief        set enefunc%table%solute_list_inv
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @param[in]    idx : index of solute list
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_solute_list_inv(i, idx)

    ! formal arguments
    integer,                 intent(in)    :: i
    integer,                 intent(in)    :: idx


    hmm_solute_list_inv(i) = idx

    return

  end subroutine hm_set_solute_list_inv

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_alloc_solute
  !> @brief        allocate solute informations
  !! @authors      NT
  !! @param[in]    num_atoms : number of atoms
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_alloc_solute(num_atoms)

    ! formal arguments
    integer,                 intent(in)    :: num_atoms

    ! local variables
    integer                  :: stat, i


    allocate(hmm_solute(num_atoms), stat=stat)
    if (stat /= 0) call hm_alloc_error('solute')

    do i = 1, num_atoms
      hmm_solute(i) = 0
    end do

    return

  end subroutine hm_alloc_solute

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_solute
  !> @brief        get enefunc%table%solute
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @return       flag for solute or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_solute(i)

    ! return value
    logical                  :: hm_solute

    ! formal arguments
    integer,                 intent(in)    :: i


    hm_solute = hmm_solute(i)

    return

  end function hm_solute

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_solute
  !> @brief        set enefunc%table%solute
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @param[in]    flag : flag for solute or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_solute(i, flag)

    ! formal arguments
    integer,                 intent(in)    :: i
    logical,                 intent(in)    :: flag


    hmm_solute(i) = flag

    return

  end subroutine hm_set_solute

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_alloc_duplicate
  !> @brief        allocate domain_constraint%duplicate
  !! @authors      NT
  !! @param[in]    num_atoms : number of atoms
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_alloc_duplicate(num_atoms)

    ! formal arguments
    integer,                 intent(in)    :: num_atoms

    ! local variables
    integer                  :: stat, i


    allocate(hmm_duplicate(num_atoms), stat=stat)
    if (stat /= 0) call hm_alloc_error('duplicate')

    do i = 1, num_atoms
      hmm_duplicate(i) = 0
    end do

    return

  end subroutine hm_alloc_duplicate

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_duplicate
  !> @brief        get domain_constraint%duplicate
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @return       duplicate
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_duplicate(i)

    ! return value
    integer                  :: hm_duplicate

    ! formal arguments
    integer,                 intent(in)    :: i


    hm_duplicate = hmm_duplicate(i)

    return

  end function hm_duplicate

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_duplicate
  !> @brief        set domain_constraint%duplicate
  !! @authors      NT
  !! @param[in]    i : index of atom
  !! @param[in]    dup : duplicate
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_duplicate(i, dup)

    ! formal arguments
    integer,                 intent(in)    :: i
    integer,                 intent(in)    :: dup


    hmm_duplicate(i) = dup

    return

  end subroutine hm_set_duplicate

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_alloc_H_index
  !> @brief        allocate domain_constraint%H_index
  !! @authors      NT
  !! @param[in]    num_atoms : number of atoms
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_alloc_H_index(num_atoms)

    ! formal arguments
    integer,                 intent(in)    :: num_atoms

    ! local variables
    integer                  :: stat, i


    allocate(hmm_H_index(8,num_atoms), stat=stat)
    if (stat /= 0) call hm_alloc_error('H_index')

    do i = 1, num_atoms
      hmm_H_index(1:8,i) = 0
    end do

    return

  end subroutine hm_alloc_H_index

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_H_index
  !> @brief        get domain_constrant%H_index
  !! @authors      NT
  !! @param[in]    j : index of duplicate
  !! @param[in]    i : index of atom
  !! @return       H index
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_H_index(j, i)

    ! return value
    integer                  :: hm_H_index

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i


    hm_H_index = hmm_H_index(j, i)

    return

  end function hm_H_index

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_H_index
  !> @brief        set domain_constrant%H_index
  !! @authors      NT
  !! @param[in]    j : index of duplicate
  !! @param[in]    i : index of atom
  !! @param[in]    idx : H index
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_H_index(j, i, idx)

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i
    integer,                 intent(in)    :: idx


    hmm_H_index(j, i) = idx

    return

  end subroutine hm_set_H_index

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_alloc_ring
  !> @brief        allocate ring
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_alloc_ring(num_atoms)

    ! formal arguments
    integer,                 intent(in)    :: num_atoms

    ! local variables
    integer                  :: stat, i, ierr


    hm_num_atoms = num_atoms

    allocate(hmm_ring           (num_atoms), stat=stat)
    if (stat /= 0) call hm_alloc_error('ring')
    allocate(hmm_num_connect    (3,num_atoms), stat=stat)
    if (stat /= 0) call hm_alloc_error('num_connect')
    allocate(hmm_connectivity   (8,3,num_atoms), stat=stat)
    if (stat /= 0) call hm_alloc_error('connectivity')

    return

  end subroutine hm_alloc_ring

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_ring
  !> @brief        set ring checker
  !! @authors      JJ
  !! @param[in]    j : index of duplicate
  !! @param[in]    i : index of atom
  !! @param[in]    idx : ring identifier
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_ring(i, idx)

    ! formal arguments
    integer,                 intent(in)    :: i
    integer,                 intent(in)    :: idx

    hmm_ring(i) = idx

    return

  end subroutine hm_set_ring

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_num_connect
  !> @brief        set num_connect
  !! @authors      JJ
  !! @param[in]    j : index of connectivity order
  !! @param[in]    i : index of atom
  !! @param[in]    idx : connectivity
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_num_connect(j, i, idx)

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i
    integer,                 intent(in)    :: idx


    hmm_num_connect(j, i) = idx

    return

  end subroutine hm_set_num_connect

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_connectivity
  !> @brief        set connectivity
  !! @authors      JJ
  !! @param[in]    k : index of atom index that connected
  !! @param[in]    j : index of connectivity order
  !! @param[in]    i : index of atom
  !! @param[in]    idx : connectivity
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_connectivity(k, j, i, idx)

    ! formal arguments
    integer,                 intent(in)    :: k
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i
    integer,                 intent(in)    :: idx


    hmm_connectivity(k, j, i) = idx

    return

  end subroutine hm_set_connectivity

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_ring
  !> @brief        get ring identifier
  !! @authors      JJ
  !! @param[in]    i : index of atom
  !! @return       connectivity number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_ring(i)

    ! return value
    integer                  :: hm_ring

    ! formal arguments
    integer,                 intent(in)    :: i


    hm_ring = hmm_ring(i)

    return

  end function hm_ring

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_num_connect
  !> @brief        get num_contact
  !! @authors      JJ
  !! @param[in]    j : index of connectivity order
  !! @param[in]    i : index of atom
  !! @return       connectivity number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_num_connect(j, i)

    ! return value
    integer                  :: hm_num_connect

    ! formal arguments
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i


    hm_num_connect = hmm_num_connect(j, i)

    return

  end function hm_num_connect

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_connectivity
  !> @brief        get connectivity
  !! @authors      JJ
  !! @param[in]    k : index of atom index that connected
  !! @param[in]    j : index of connectivity order
  !! @param[in]    i : index of atom
  !! @return       connectivity
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_connectivity(k, j, i)

    ! return value
    integer                  :: hm_connectivity

    ! formal arguments
    integer,                 intent(in)    :: k
    integer,                 intent(in)    :: j
    integer,                 intent(in)    :: i


    hm_connectivity = hmm_connectivity(k, j, i)

    return

  end function hm_connectivity

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_alloc_selatoms
  !> @brief        allocate restraints%atomlist
  !! @authors      NT
  !! @param[in]    num_sels : number of selects
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_alloc_selatoms(num_sels)

    ! formal arguments
    integer,                 intent(in)    :: num_sels

    ! local variables
    integer                  :: stat, isel, i
    character(24)            :: ssel


    hm_num_selatoms_list = hm_num_selatoms_list + 1

    isel = hm_num_selatoms_list
    write(ssel,'(i0)') isel

    hm_num_selatoms(isel) = num_sels

    allocate(hmm_selatoms(isel)%v(num_sels), stat=stat)
    if (stat /= 0) call hm_alloc_error('selatoms')

    do i = 1, num_sels
      hmm_selatoms(isel)%v(i) = 0
    end do

    return

  end subroutine hm_alloc_selatoms

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_selatoms
  !> @brief        get selatoms
  !! @authors      NT
  !! @param[in]    i    : index of atom selection
  !! @param[in]    isel : index of selection
  !! @return       index of atom
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_selatoms(i, isel)

    ! return value
    integer                  :: hm_selatoms

    ! formal arguments
    integer,                 intent(in)    :: i
    integer,                 intent(in)    :: isel


    hm_selatoms = hmm_selatoms(isel)%v(i)

    return

  end function hm_selatoms

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_selatoms
  !> @brief        set selatoms
  !! @authors      NT
  !! @param[in]    i    : index of atom selection
  !! @param[in]    isel : index of selection
  !! @param[in]    idx  : index of atom
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_selatoms(i, isel, idx)

    ! formal arguments
    integer,                 intent(in)    :: i
    integer,                 intent(in)    :: isel
    integer,                 intent(in)    :: idx


    hmm_selatoms(isel)%v(i) = idx

    return

  end subroutine hm_set_selatoms

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_alloc_error
  !> @brief        show allocation error message
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_alloc_error(var_name)

    ! formal arguments
    character(*),            intent(in)    :: var_name


    if (main_rank) &
      call error_msg('ERROR: Allocation is failed. :['//&
      var_name//']. Molecule is too big. Do not Use HM_USE_MEMORY.')

    return

  end subroutine hm_alloc_error

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_check_memory
  !> @brief        check the memory usage amount
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_check_memory

    ! constants
    integer(8),    parameter :: BytePerAtom = 208

    ! local variables
    integer(8)               :: total, used, free, shared, buffers, cached
    integer(8)               :: mem_per_proc, nproc_per_node, mem_per_node
    integer                  :: pname_len, ierr, unit_no, iproc
    character(256)           :: line
    character(128)           :: filename, tmp
    character(MPI_MAX_PROCESSOR_NAME) :: pname

    character(MPI_MAX_PROCESSOR_NAME),allocatable :: pnames(:)


    ! check the free memory size
    unit_no = get_unit_no()

    write(filename,'(a,i0)') '/tmp/_hm_check_memory', my_prst_rank
    call system('free -b > '//filename)
    open(unit_no, file=filename, status='old', form='formatted', iostat=ierr)
    if (ierr /= 0) then
      write(filename,'(a,i0)') '/var/tmp/_hm_check_memory', my_prst_rank
      call system('free -b > '//filename)
      open(unit_no, file=filename, status='old', form='formatted', iostat=ierr)
      if (ierr /= 0) goto 10
    end if

    read(unit_no,'(a)') line
    read(unit_no,'(a)') line

    read(line,*,err=10,end=10) tmp, total, used, free, shared, buffers, cached

    close(unit_no)
    call free_unit_no(unit_no)
    call system('rm -f '//filename)

    free = total - (used - cached)

    ! check the required memory size per process
    mem_per_proc = int(hm_num_atoms,8) * BytePerAtom

    ! check the number of process per node
    allocate(pnames(nproc_prst))

    call MPI_Get_processor_name(pname, pname_len, ierr)

    call MPI_Allgather(pname,  MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, &
                       pnames, MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, &
                       MPI_COMM_WORLD, ierr)

    nproc_per_node = 0
    do iproc = 1, nproc_prst
      if (pname == pnames(iproc)) &
        nproc_per_node = nproc_per_node + 1
    end do

    deallocate(pnames)

    ! last check
    mem_per_node = mem_per_proc * nproc_per_node

    if (free < mem_per_node) then
      write(MsgOut,*)''
      write(MsgOut,*)'ERROR: Out of Memory.'
      write(MsgOut,*)'  Rank:', my_prst_rank
      write(MsgOut,*)'  Free memory size: ',free/1024/1024, ' MB'
      write(MsgOut,*)'  Required memory : ',mem_per_node/1024/1024,' MB /node'
      stop
    end if

    return

10  write(MsgOut,*) ' ##ERROR: Memory check failed.'
    close(unit_no)
    call free_unit_no(unit_no)
    call system('rm -f '//filename)
    return

  end subroutine hm_check_memory

end module pr_huge_molecule_mod
