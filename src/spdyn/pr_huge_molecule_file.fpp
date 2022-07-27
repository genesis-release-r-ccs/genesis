!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pr_huge_molecule_mod
!> @brief   huge molecule information (file version)
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

  ! s_selatoms
  public :: hm_alloc_selatoms
  public :: hm_selatoms
  public :: hm_set_selatoms


  ! private
  !

  integer                   :: hmu_atom_no         = -1
  integer                   :: hmu_atom_cls_no     = -1
  integer                   :: hmu_atom_cls_name   = -1
  integer                   :: hmu_atom_name       = -1
  integer                   :: hmu_light_atom_mass = -1
  integer                   :: hmu_light_atom_name = -1
  integer                   :: hmu_atom_coord      = -1
  integer                   :: hmu_atom_refcoord   = -1
  integer                   :: hmu_atom_velocity   = -1
  integer                   :: hmu_charge          = -1
  integer                   :: hmu_mass            = -1
  integer                   :: hmu_residue_no      = -1
  integer                   :: hmu_residue_name    = -1
  integer                   :: hmu_segment_name    = -1
  integer                   :: hmu_molecule_no     = -1
  integer                   :: hmu_bond_list       = -1
  integer                   :: hmu_angl_list       = -1
  integer                   :: hmu_dihe_list       = -1
  integer                   :: hmu_impr_list       = -1
  integer                   :: hmu_cmap_list       = -1
  integer                   :: hmu_water_list      = -1
  integer                   :: hmu_solute_list     = -1
  integer                   :: hmu_solute_list_inv = -1
  integer                   :: hmu_solute          = -1
  integer                   :: hmu_duplicate       = -1
  integer                   :: hmu_H_index         = -1
  integer                   :: hmu_selatoms(HMMaxSelAtomsList) = -1

  ! subroutines
  private :: hm_get_work_name

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

    call system('mkdir -p '//trim(work_path))
    call system('rm -f '   //trim(work_path)//'/'//trim(hm_work_name)//'_*')

    if (main_rank) then

      write(MsgOut,'(a)') 'Hm_Initialize> working file : '
      write(MsgOut,'(a)') '  working directory   = ' // trim(hm_work_dir)
      write(MsgOut,'(a)') '  working file prefix = ' // trim(hm_work_name)
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

      write(MsgOut,'(a)') 'Hm_Finalize> cleanup working files.'
      write(MsgOut,'(a)') ' '

    end if

    call pr_close(hmu_atom_no)
    call pr_close(hmu_atom_cls_no)
    call pr_close(hmu_atom_cls_name)
    call pr_close(hmu_atom_name)
    call pr_close(hmu_light_atom_mass)
    call pr_close(hmu_light_atom_name)
    call pr_close(hmu_atom_coord)
    call pr_close(hmu_atom_refcoord)
    call pr_close(hmu_atom_velocity)
    call pr_close(hmu_charge)
    call pr_close(hmu_mass)
    call pr_close(hmu_residue_no)
    call pr_close(hmu_residue_name)
    call pr_close(hmu_segment_name)
    call pr_close(hmu_molecule_no)
    call pr_close(hmu_bond_list)
    call pr_close(hmu_angl_list)
    call pr_close(hmu_dihe_list)
    call pr_close(hmu_impr_list)
    call pr_close(hmu_cmap_list)
    call pr_close(hmu_water_list)
    call pr_close(hmu_solute_list)
    call pr_close(hmu_solute_list_inv)
    call pr_close(hmu_solute)
    call pr_close(hmu_duplicate)
    call pr_close(hmu_H_index)
    do i = 1, HMMaxSelAtomsList
      call pr_close(hmu_selatoms(i))
    end do

    call system('rm -f '//trim(hm_work_dir)//'/'//trim(hm_work_name)//'_*')

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


    hm_num_atoms = num_atoms

    call pr_open(hm_get_work_name('atom_no'),       &
         int(num_atoms * 4,8),   hmu_atom_no)
    call pr_open(hm_get_work_name('atom_cls_no'),   &
         int(num_atoms * 4,8),   hmu_atom_cls_no)
    call pr_open(hm_get_work_name('atom_cls_name'), &
         int(num_atoms * 6,8),   hmu_atom_cls_name)
    call pr_open(hm_get_work_name('atom_name'),     &
         int(num_atoms * 4,8),   hmu_atom_name)
    call pr_open(hm_get_work_name('light_atom_mass'), &
         int(num_atoms * 4,8),   hmu_light_atom_mass)
    call pr_open(hm_get_work_name('light_atom_name'), &
         int(num_atoms * 4,8),   hmu_light_atom_name)
    call pr_open(hm_get_work_name('atom_coord'),    &
         int(num_atoms * 8*3,8), hmu_atom_coord)
    call pr_open(hm_get_work_name('atom_refcoord'), &
         int(num_atoms * 8*3,8), hmu_atom_refcoord)
    call pr_open(hm_get_work_name('atom_velocity'), &
         int(num_atoms * 8*3,8), hmu_atom_velocity)
    call pr_open(hm_get_work_name('charge'),        &
         int(num_atoms * 8,8),   hmu_charge)
    call pr_open(hm_get_work_name('mass'),          &
         int(num_atoms * 8,8),   hmu_mass)
    call pr_open(hm_get_work_name('residue_no'),    &
         int(num_atoms * 4,8),   hmu_residue_no)
    call pr_open(hm_get_work_name('residue_name'),  &
         int(num_atoms * 6,8),   hmu_residue_name)
    call pr_open(hm_get_work_name('segment_name'),  &
         int(num_atoms * 4,8),   hmu_segment_name)
    call pr_open(hm_get_work_name('molecule_no'),   &
         int(num_atoms * 4,8),   hmu_molecule_no)

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


    call pr_read(hmu_atom_no, int((i-1)*4,8), hm_atom_no, 4)

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


    call pr_read(hmu_atom_cls_no, int((i-1)*4,8), hm_atom_cls_no, 4)

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


    call pr_read(hmu_atom_cls_name, int((i-1)*6,8), hm_atom_cls_name, 6)

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


    call pr_read(hmu_atom_name, int((i-1)*4,8), hm_atom_name, 4)

    return

  end function hm_atom_name

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_light_atom_mass
  !> @brief        get molecule%light_atom_mass
  !! @authors      JJ
  !! @param[in]    i : index of atom
  !! @return       light atom mass key
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_light_atom_mass(i)

    ! return value
    integer                  :: hm_light_atom_mass

    ! formal arguments
    integer,                 intent(in)    :: i


    call pr_read(hmu_light_atom_mass, int((i-1)*4,8), hm_light_atom_mass, 4)

    return

  end function hm_light_atom_mass

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_light_atom_name
  !> @brief        get molecule%light_atom_name
  !! @authors      JJ
  !! @param[in]    i : index of atom
  !! @return       light atom name key
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_light_atom_name(i)

    ! return value
    integer                  :: hm_light_atom_name

    ! formal arguments
    integer,                 intent(in)    :: i


    call pr_read(hmu_light_atom_name, int((i-1)*4,8), hm_light_atom_name, 4)

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

    ! local variables
    real(dp)                 :: v


    call pr_read(hmu_atom_coord, int(((i-1)*3+(j-1)),8) * 8, v, 8)
    hm_atom_coord = real(v,wip)

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

    ! local variables
    real(dp)                 :: v


    call pr_read(hmu_atom_refcoord, int(((i-1)*3+(j-1))*8,8), v, 8)
    hm_atom_refcoord = real(v,wip)

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

    ! local variables
    real(dp)                 :: v


    call pr_read(hmu_atom_velocity, int(((i-1)*3+(j-1))*8,8), v, 8)
    hm_atom_velocity = real(v,wip)

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


    call pr_read(hmu_charge, int((i-1) * 8,8), hm_charge, 8)

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

    ! local variables
    real(dp)                 :: v


    call pr_read(hmu_mass, int((i-1) * 8,8), v, 8)
    hm_mass = real(v,wip)

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


    call pr_read(hmu_residue_no, int((i-1) * 4,8), hm_residue_no, 4)

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


    call pr_read(hmu_residue_name, int((i-1) * 6,8), hm_residue_name, 6)

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


    call pr_read(hmu_segment_name, int((i-1) * 4,8), hm_segment_name, 4)

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


    call pr_read(hmu_molecule_no, int((i-1) * 4,8), hm_molecule_no, 4)

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


    call pr_write(hmu_atom_no, int((i-1)*4,8), no, 4)

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


    call pr_write(hmu_atom_cls_no, int((i-1)*4,8), no, 4)

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

    character(6)             :: tmp


    tmp = name
    call pr_write(hmu_atom_cls_name, int((i-1)*6,8), tmp, 6)

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

    character(4)             :: tmp


    tmp = name
    call pr_write(hmu_atom_name, int((i-1)*4,8), tmp, 4)

    return

  end subroutine hm_set_atom_name

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_light_atom_mass
  !> @brief        set molecule%light_atom_mass
  !! @authors      JJ
  !! @param[in]    i : index of atom
  !! @param[in]    lam : light atom mass
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_light_atom_mass(i, lam)

    ! formal arguments
    integer,                 intent(in)    :: i
    integer,                 intent(in)    :: lam


    call pr_write(hmu_light_atom_mass, int((i-1) * 4,8), lam, 4)

    return

  end subroutine hm_set_light_atom_mass

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_set_light_atom_name
  !> @brief        set molecule%light_atom_name
  !! @authors      JJ
  !! @param[in]    i : index of atom
  !! @param[in]    lan : light atom name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_set_light_atom_name(i, lan)

    ! formal arguments
    integer,                 intent(in)    :: i
    integer,                 intent(in)    :: lan


    call pr_write(hmu_light_atom_name, int((i-1) * 4,8), lan, 4)

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

    ! local variables
    real(dp)                 :: v


    v = real(crd,dp)
    call pr_write(hmu_atom_coord, int(((i-1)*3+(j-1)) * 8,8), v, 8)

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

    ! local variables
    real(dp)                 :: v


    v = real(crd,dp)
    call pr_write(hmu_atom_refcoord, int(((i-1)*3+(j-1)) * 8,8), v, 8)
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

    ! local variables
    real(dp)                 :: v


    v = real(crd,dp)
    call pr_write(hmu_atom_velocity, int(((i-1)*3+(j-1)) * 8,8), v, 8)

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


    call pr_write(hmu_charge, int((i-1) * 8,8), charge, 8)

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

    ! local variables
    real(dp)                 :: v


    v = real(mass,dp)
    call pr_write(hmu_mass, int((i-1) * 8,8), v, 8)

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


    call pr_write(hmu_residue_no, int((i-1) * 4,8), no, 4)

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

    character(6)             :: tmp


    tmp = name
    call pr_write(hmu_residue_name, int((i-1) * 6,8), tmp, 6)

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

    character(4)             :: tmp


    tmp = name
    call pr_write(hmu_segment_name, int((i-1) * 4,8), tmp, 4)

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


    call pr_write(hmu_molecule_no, int((i-1) * 4,8), no, 4)

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


    hm_num_bonds = num_bonds

    call pr_open(hm_get_work_name('bond_list'), &
         int(num_bonds * 4 * 2,8), hmu_bond_list)

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


    call pr_read(hmu_bond_list, int(((i-1)*2+(j-1)) * 4,8), hm_bond_list, 4)

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


    call pr_write(hmu_bond_list, int(((i-1)*2+(j-1)) * 4,8), idx, 4)

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


    hm_num_angles = num_angls

    call pr_open(hm_get_work_name('angl_list'), &
         int(num_angls * 4 * 3,8), hmu_angl_list)

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


    call pr_read(hmu_angl_list, int(((i-1)*3+(j-1)) * 4,8), hm_angl_list, 4)

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


    call pr_write(hmu_angl_list, int(((i-1)*3+(j-1)) * 4,8), idx, 4)

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


    hm_num_dihedrals = num_dihes

    call pr_open(hm_get_work_name('dihe_list'), &
         int(num_dihes * 4 * 4,8), hmu_dihe_list)

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


    call pr_read(hmu_dihe_list, int(((i-1)*4+(j-1)) * 4,8), hm_dihe_list, 4)

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


    call pr_write(hmu_dihe_list, int(((i-1)*4+(j-1)) * 4,8), idx, 4)

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


    hm_num_impropers = num_imprs

    call pr_open(hm_get_work_name('impr_list'), &
         int(num_imprs * 4 * 4,8), hmu_impr_list)

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


    call pr_read(hmu_impr_list, int(((i-1)*4+(j-1)) * 4,8), hm_impr_list, 4)

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


    call pr_write(hmu_impr_list, int(((i-1)*4+(j-1)) * 4,8), idx, 4)

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


    hm_num_cmaps = num_cmaps

    call pr_open(hm_get_work_name('cmap_list'), &
         int(num_cmaps * 4 * 8,8), hmu_cmap_list)

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


    call pr_read(hmu_cmap_list, int(((i-1)*8+(j-1)) * 4,8), hm_cmap_list, 4)

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


    call pr_write(hmu_cmap_list, int(((i-1)*8+(j-1)) * 4,8), idx, 4)

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


    hm_num_water_list = num_waters

    call pr_open(hm_get_work_name('water_list'), &
         int(num_waters * 4 * na,8), hmu_water_list)

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


    call pr_read(hmu_water_list, int(((i-1)*3+(j-1)) * 4,8), hm_water_list, 4)

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


    call pr_write(hmu_water_list, int(((i-1)*3+(j-1)) * 4,8), idx, 4)

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


    hm_num_solute_list = num_solutes

    call pr_open(hm_get_work_name('solute_list'), &
         int(num_solutes * 4,8), hmu_solute_list)

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


    call pr_read(hmu_solute_list, int((i-1) * 4,8), hm_solute_list, 4)

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


    call pr_write(hmu_solute_list, int((i-1) * 4,8), idx, 4)

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


    call pr_open(hm_get_work_name('solute_list_inv'), &
         int(num_atoms * 4,8), hmu_solute_list_inv)

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


    call pr_read(hmu_solute_list_inv, int((i-1) * 4,8), hm_solute_list_inv, 4)

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


    call pr_write(hmu_solute_list_inv, int((i-1) * 4,8), idx, 4)

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


    call pr_open(hm_get_work_name('solute'), &
         int(num_atoms,8), hmu_solute)

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

    ! local variables
    integer(1)               :: flag


    call pr_read(hmu_solute, int(i-1,8), flag, 1)
    hm_solute = flag > 0

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

    ! local variables
    integer(1)               :: iflag


    if (flag) then
      iflag = 1
    else
      iflag = 0
    end if
    call pr_write(hmu_solute, int(i-1,8), iflag, 1)

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


    call pr_open(hm_get_work_name('duplicate'), &
         int(num_atoms * 4,8), hmu_duplicate)

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


    call pr_read(hmu_duplicate, int((i-1) * 4,8), hm_duplicate, 4)

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


    call pr_write(hmu_duplicate, int((i-1) * 4,8), dup, 4)

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


    call pr_open(hm_get_work_name('H_index'), &
         int(num_atoms * 4 * 8,8), hmu_H_index)

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


    call pr_read(hmu_H_index, int(((i-1)*8+(j-1)) * 4,8), hm_H_index, 4)

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


    call pr_write(hmu_H_index, int(((i-1)*8+(j-1)) * 4,8), idx, 4)

    return

  end subroutine hm_set_H_index

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
    integer                  :: isel
    character(24)            :: ssel


    hm_num_selatoms_list = hm_num_selatoms_list + 1

    isel = hm_num_selatoms_list
    write(ssel,'(i0)') isel

    hm_num_selatoms(isel) = num_sels

    call pr_open(hm_get_work_name('selatoms'//ssel), &
         int(num_sels * 4,8), hmu_selatoms(isel))

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


    call pr_read(hmu_selatoms(isel), int((i-1)*4,8), hm_selatoms, 4)

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


    call pr_write(hmu_selatoms(isel), int((i-1)*4,8), idx, 4)

    return

  end subroutine hm_set_selatoms

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      hm_get_work_name
  !> @brief        get working filename
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function hm_get_work_name(type_str)

    ! return values
    character(MaxFilename)        :: hm_get_work_name

    ! formal arguments
    character(*), intent(in) :: type_str


    write(hm_get_work_name, '(a,a,a,a,i0,a,a)') &
         trim(hm_work_dir),  '/', &
         trim(hm_work_name), '_', &
         my_prst_rank, '_', &
         trim(type_str)

    return

  end function hm_get_work_name

end module pr_huge_molecule_mod
