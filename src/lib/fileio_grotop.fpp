!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_grotop_mod
!> @brief   read GROMACS topology file
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fileio_grotop_mod

  use fileio_gropp_mod
  use fileio_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private


  ! GROMACS TOP molecule structures

  ! atom data
  type, public :: s_atom
    integer                        :: atom_no       = 0
    character(6)                   :: atom_type     = ''
    integer                        :: residue_no    = 0
    character(6)                   :: residue_name  = ''
    character(4)                   :: atom_name     = ''
    integer                        :: charge_group  = 0
    real(wp)                       :: charge        = 0.0_wp
    real(wp)                       :: mass          = 0.0_wp
    real(wp)                       :: stokes_r      = 0.0_wp
  end type s_atom

  ! bond data
  type, public :: s_bond
    integer                        :: atom_idx1     = 0
    integer                        :: atom_idx2     = 0
    integer                        :: func          = 0
    real(wp)                       :: b0            = 0.0_wp !nm
    real(wp)                       :: kb            = 0.0_wp
  end type s_bond

  ! angle data
  type, public :: s_angl
    integer                        :: atom_idx1     = 0
    integer                        :: atom_idx2     = 0
    integer                        :: atom_idx3     = 0
    integer                        :: func          = 0
    real(wp)                       :: theta_0       = 0.0_wp !deg.
    real(wp)                       :: kt            = 0.0_wp
    real(wp)                       :: r13           = 0.0_wp !nm
    real(wp)                       :: kub           = 0.0_wp
  end type s_angl

  ! dihedral data
  type, public :: s_dihe
    integer                        :: atom_idx1     = 0
    integer                        :: atom_idx2     = 0
    integer                        :: atom_idx3     = 0
    integer                        :: atom_idx4     = 0
    integer                        :: func          = 0
    real(wp)                       :: ps            = 0.0_wp !deg.
    real(wp)                       :: kp            = 0.0_wp
    integer                        :: multiplicity  = 0
    real(wp)                       :: c(6)          = 0.0_wp
  end type s_dihe

  ! cmap data
  type, public :: s_cmap
    integer                        :: atom_idx1     = 0
    integer                        :: atom_idx2     = 0
    integer                        :: atom_idx3     = 0
    integer                        :: atom_idx4     = 0
    integer                        :: atom_idx5     = 0
    integer                        :: func          = 0
    real(wp),          allocatable :: map(:,:)
  end type s_cmap

  ! exclusions data
  type, public :: s_excl
    integer                        :: atom_idx1     = 0
    integer                        :: atom_idx2     = 0
  end type s_excl

  ! constraints data
  type, public :: s_constr
    integer                        :: atom_idx1     = 0
    integer                        :: atom_idx2     = 0
    integer                        :: func          = 0
    real(wp)                       :: b0            = 0.0_wp !nm
  end type s_constr

  ! pair data
  !    func_type = 1  : LJ
  type, public :: s_pair
    integer                        :: atom_idx1     = 0
    integer                        :: atom_idx2     = 0
    integer                        :: func          = 0
    real(wp)                       :: v             = 0.0_wp
    real(wp)                       :: w             = 0.0_wp
    real(wp)                       :: khh           = 0.0_wp
  end type s_pair

  ! settle data
  type, public :: s_settles
    integer                        :: ow            = 0
    integer                        :: func          = 0
    real(wp)                       :: doh           = 0.0_wp
    real(wp)                       :: dhh           = 0.0_wp
  end type s_settles

  ! virtual site 2 data
  type, public :: s_vsite2
    integer                        :: atom_idx1     = 0 !site
    integer                        :: atom_idx2     = 0
    integer                        :: atom_idx3     = 0
    integer                        :: func          = 0
    real(wp)                       :: a             = 0.0_wp
  end type s_vsite2

  ! virtual site 3 data
  type, public :: s_vsite3
    integer                        :: atom_idx1     = 0 !site
    integer                        :: atom_idx2     = 0
    integer                        :: atom_idx3     = 0
    integer                        :: atom_idx4     = 0
    integer                        :: func          = 0
    real(wp)                       :: a             = 0.0_wp
    real(wp)                       :: b             = 0.0_wp
    real(wp)                       :: c             = 0.0_wp
    real(wp)                       :: d             = 0.0_wp
    real(wp)                       :: theta         = 0.0_wp
  end type s_vsite3

  ! virtual site 4 data
  type, public :: s_vsite4
    integer                        :: atom_idx1     = 0 !site
    integer                        :: atom_idx2     = 0
    integer                        :: atom_idx3     = 0
    integer                        :: atom_idx4     = 0
    integer                        :: atom_idx5     = 0
    integer                        :: func          = 0
    real(wp)                       :: a             = 0.0_wp
    real(wp)                       :: b             = 0.0_wp
    real(wp)                       :: c             = 0.0_wp
  end type s_vsite4

  ! virtual site n data
  type, public :: s_vsiten
    integer                        :: atom_idx      = 0 !site
    integer                        :: func          = 0
    integer,           allocatable :: atom_idcs(:)
  end type s_vsiten

  ! position restraints data
  type, public :: s_posres
    integer                        :: atom_idx      = 0
    integer                        :: func          = 0
    real(wp)                       :: kx            = 0.0_wp
    real(wp)                       :: ky            = 0.0_wp
    real(wp)                       :: kz            = 0.0_wp
  end type s_posres

  ! molecule definition data
  type, public :: s_grotop_mol
    integer                        :: num_atoms     = 0
    integer                        :: num_bonds     = 0
    integer                        :: num_angls     = 0
    integer                        :: num_dihes     = 0
    integer                        :: num_cmaps     = 0
    integer                        :: num_excls     = 0
    integer                        :: num_constrs   = 0
    integer                        :: num_pairs     = 0
    integer                        :: num_vsites2   = 0
    integer                        :: num_vsites3   = 0
    integer                        :: num_vsites4   = 0
    integer                        :: num_vsitesn   = 0
    integer                        :: num_posress   = 0
    type(s_atom),      allocatable :: atoms(:)
    type(s_bond),      allocatable :: bonds(:)
    type(s_angl),      allocatable :: angls(:)
    type(s_dihe),      allocatable :: dihes(:)
    type(s_cmap),      allocatable :: cmaps(:)
    type(s_excl),      allocatable :: excls(:)
    type(s_constr),    allocatable :: constrs(:)
    type(s_pair),      allocatable :: pairs(:)
    type(s_settles)                :: settles
    type(s_vsite2),    allocatable :: vsites2(:)
    type(s_vsite3),    allocatable :: vsites3(:)
    type(s_vsite4),    allocatable :: vsites4(:)
    type(s_vsiten),    allocatable :: vsitesn(:)
    type(s_posres),    allocatable :: posress(:)
  end type s_grotop_mol


  ! GROMACS TOP structures

  ! defaults data
  type, public :: s_defaults
    integer                        :: nonb_func    = 0
    integer                        :: combi_rule   = 0
    logical                        :: gen_pairs    = .false.
    real(wp)                       :: fudge_lj     = 1.0_wp
    real(wp)                       :: fudge_qq     = 1.0_wp
  end type s_defaults

  ! atom type data
  type, public :: s_atomtype
    character(6)                   :: type_name    = ''
    integer                        :: type_no      = 0
    real(wp)                       :: mass         = 0.0_wp
    real(wp)                       :: charge       = 0.0_wp
    character(4)                   :: ptype        = ''
    real(wp)                       :: v            = 0.0_wp
    real(wp)                       :: w            = 0.0_wp
    real(wp)                       :: khh          = 0.0_wp
    real(wp)                       :: stokes_r     = 0.0_wp
    character(20)                  :: alias        = ''
  end type s_atomtype

  ! bond type data
  !    func_type = 1  : bond
  !    func_type = 2  : G96 bond
  type, public :: s_bondtype
    character(6)                   :: atom_type1   = ''
    character(6)                   :: atom_type2   = ''
    integer                        :: func         = 0
    real(wp)                       :: b0           = 0.0_wp !nm
    real(wp)                       :: kb           = 0.0_wp
  end type s_bondtype

  ! angle type data
  !    func_type = 1  : angle
  !    func_type = 2  : G96 angle
  type, public :: s_angltype
    character(6)                   :: atom_type1   = ''
    character(6)                   :: atom_type2   = ''
    character(6)                   :: atom_type3   = ''
    integer                        :: func         = 0
    real(wp)                       :: theta_0      = 0.0_wp !deg.
    real(wp)                       :: kt           = 0.0_wp
    real(wp)                       :: r13          = 0.0_wp !nm
    real(wp)                       :: kub          = 0.0_wp
  end type s_angltype

  ! dihedral type data
  !    func_type = 1  : proper dih.
  !    func_type = 2  : improper dih.
  !    func_type = 3  : Ryckaert-Bellemans dih.
  type, public :: s_dihetype
    character(6)                   :: atom_type1   = ''
    character(6)                   :: atom_type2   = ''
    character(6)                   :: atom_type3   = ''
    character(6)                   :: atom_type4   = ''
    integer                        :: func         = 0
    real(wp)                       :: ps           = 0.0_wp !deg.
    real(wp)                       :: kp           = 0.0_wp
    integer                        :: multiplicity = 0
    real(wp)                       :: c(6)         = 0.0_wp
    logical                        :: four_type    = .false.
  end type s_dihetype

  ! constraint type data
  !    func_type = 1  : bond
  !    func_type = 2  : G96 bond
  type, public :: s_constrtype
    character(6)                   :: atom_type1   = ''
    character(6)                   :: atom_type2   = ''
    integer                        :: func         = 0
    real(wp)                       :: b0           = 0.0_wp !nm
  end type s_constrtype

  ! cmap type data
  !    func_type = 1  : ?
  type, public :: s_cmaptype
    character(6)                   :: atom_type1   = ''
    character(6)                   :: atom_type2   = ''
    character(6)                   :: atom_type3   = ''
    character(6)                   :: atom_type4   = ''
    character(6)                   :: atom_type5   = ''
    integer                        :: func         = 0
    real(wp),          allocatable :: map(:,:)
  end type s_cmaptype

  ! non-bonded parameter data
  !    func_type = 1  : LJ
  type, public :: s_nbonparm
    character(6)                   :: atom_type1   = ''
    character(6)                   :: atom_type2   = ''
    integer                        :: func         = 0
    real(wp)                       :: v            = 0.0_wp
    real(wp)                       :: w            = 0.0_wp
    real(wp)                       :: khh          = 0.0_wp
  end type s_nbonparm

  ! molecule type data
  type, public :: s_moltype
    character(20)                  :: name         = ''
    integer                        :: exclude_nbon = 0
    type(s_grotop_mol),  pointer   :: mol          => null()
  end type s_moltype

  ! molecules
  type, public :: s_mols
    character(20)                  :: name
    integer                        :: count
    type(s_moltype),     pointer   :: moltype      => null()
  end type s_mols

  ! topology data
  type, public :: s_grotop
    integer                        :: num_atomtypes   = 0
    integer                        :: num_bondtypes   = 0
    integer                        :: num_angltypes   = 0
    integer                        :: num_dihetypes   = 0
    integer                        :: num_constrtypes = 0
    integer                        :: num_cmaptypes   = 0
    integer                        :: num_nbonparms   = 0
    integer                        :: num_moltypes    = 0
    integer                        :: num_molss       = 0
    logical                        :: alias_atomtype  = .false.
    type(s_defaults)               :: defaults
    type(s_atomtype),  allocatable :: atomtypes(:)
    type(s_bondtype),  allocatable :: bondtypes(:)
    type(s_angltype),  allocatable :: angltypes(:)
    type(s_dihetype),  allocatable :: dihetypes(:)
    type(s_constrtype),allocatable :: constrtypes(:)
    type(s_cmaptype),  allocatable :: cmaptypes(:)
    type(s_nbonparm),  allocatable :: nbonparms(:)
    type(s_moltype),   allocatable :: moltypes(:)
    type(s_mols),      allocatable :: molss(:)
  end type s_grotop

  ! parameters for TOP molecule structure allocatable variables
  integer,      public, parameter :: GMGroMolAtom     = 101
  integer,      public, parameter :: GMGroMolBond     = 102
  integer,      public, parameter :: GMGroMolAngl     = 103
  integer,      public, parameter :: GMGroMolDihe     = 104
  integer,      public, parameter :: GMGroMolCmap     = 105
  integer,      public, parameter :: GMGroMolExcl     = 106
  integer,      public, parameter :: GMGroMolConstr   = 107
  integer,      public, parameter :: GMGroMolPair     = 108
  integer,      public, parameter :: GMGroMolVsites2  = 109
  integer,      public, parameter :: GMGroMolVsites3  = 110
  integer,      public, parameter :: GMGroMolVsites4  = 111
  integer,      public, parameter :: GMGroMolVsitesn  = 112
  integer,      public, parameter :: GMGroMolPosres   = 113

  ! parameters for TOP structure allocatable variables
  integer,      public, parameter :: GroTopAtomType   = 1
  integer,      public, parameter :: GroTopBondType   = 2
  integer,      public, parameter :: GroTopAnglType   = 3
  integer,      public, parameter :: GroTopDiheType   = 4
  integer,      public, parameter :: GroTopConstrType = 5
  integer,      public, parameter :: GroTopCmapType   = 6
  integer,      public, parameter :: GroTopNbonParm   = 7
  integer,      public, parameter :: GroTopMolType    = 8
  integer,      public, parameter :: GroTopMols       = 9

  ! parameters for directive
  integer,     private, parameter :: DDefaults        = 1
  integer,     private, parameter :: DAtomTypes       = 2
  integer,     private, parameter :: DBondTypes       = 3
  integer,     private, parameter :: DAngleTypes      = 4
  integer,     private, parameter :: DDihedralTypes   = 5
  integer,     private, parameter :: DConstrTypes     = 6
  integer,     private, parameter :: DCmapTypes       = 7
  integer,     private, parameter :: DNonbondParams   = 8
  !    molecule directive
  integer,     private, parameter :: DMoleculeType    = 9
  integer,     private, parameter :: DAtoms           = 10
  integer,     private, parameter :: DBonds           = 11
  integer,     private, parameter :: DAngles          = 12
  integer,     private, parameter :: DDihedrals       = 13
  integer,     private, parameter :: DCmap            = 14
  integer,     private, parameter :: DExclusions      = 15
  integer,     private, parameter :: DConstraints     = 16
  integer,     private, parameter :: DPairs           = 17
  integer,     private, parameter :: DSettles         = 18
  integer,     private, parameter :: DVirtualSites2   = 19
  integer,     private, parameter :: DVirtualSites3   = 20
  integer,     private, parameter :: DVirtualSites4   = 21
  integer,     private, parameter :: DVirtualSitesN   = 22
  integer,     private, parameter :: DPosition_Restraints = 23
  !    system directive
  integer,     private, parameter :: DSystem          = 24
  integer,     private, parameter :: DMolecules       = 25
  !    unknown directive
  integer,     private, parameter :: DUnknown         = 26

  ! parameters
  logical,     private, parameter :: VerboseOn        = .false.

  ! subroutines
  public  :: input_grotop
  public  :: output_grotop

  public  :: init_grotop
  public  :: init_grotop_mol
  public  :: realloc_grotop
  public  :: realloc_grotop_mol
  public  :: dealloc_grotop_all
  public  :: dealloc_grotop_mol_all

  private :: read_grotop
  private :: read_defaults
  private :: read_atom_types
  private :: read_bond_types
  private :: read_angle_types
  private :: read_dihedral_types
  private :: read_constraint_types
  private :: read_cmap_types
  private :: read_nonbond_parms
  private :: read_molecule_type
  private :: read_atoms
  private :: read_bonds
  private :: read_angles
  private :: read_dihedrals
  private :: read_cmaps
  private :: read_exclusions
  private :: read_constraints
  private :: read_pairs
  private :: read_settles
  private :: read_virtual_sites2
  private :: read_virtual_sites3
  private :: read_virtual_sites4
  private :: read_virtual_sitesn
  private :: read_position_restraints
  private :: read_system
  private :: read_molecules

  private :: write_grotop
  private :: write_defaults
  private :: write_atom_types
  private :: write_bond_types
  private :: write_angle_types
  private :: write_dihedral_types
  private :: write_constraint_types
  private :: write_cmap_types
  private :: write_nonbond_parms
  private :: write_molecule_type
  private :: write_atoms
  private :: write_bonds
  private :: write_angles
  private :: write_dihedrals
  private :: write_cmaps
  private :: write_exclusions
  private :: write_constraints
  private :: write_pairs
  private :: write_settles
  private :: write_virtual_sites2
  private :: write_virtual_sites3
  private :: write_virtual_sites4
  private :: write_virtual_sitesn
  private :: write_position_restraints
  private :: write_system
  private :: write_molecules

  private :: check_section_count
  private :: check_directive
  private :: size_grotop
  private :: size_grotop_mol
  private :: error_msg_grotop
  private :: search_atom_mass
  private :: search_stokes_r
  private :: search_bond_param
  private :: search_angl_param
  private :: search_dihe_param
  private :: search_dihe_param_duplicate
  private :: search_constr_param
  private :: search_cmap_param
  private :: rename_atom_type_name
  private :: match_format
  private :: is_digit

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_grotop
  !> @brief        a driver subroutine for reading GRO TOP file
  !! @authors      NT
  !! @param[in]    grotop_filename : filename of GRO TOP file
  !! @param[inout] grotop          : structure of GRO TOP data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_grotop(grotop_filename, grotop)

    ! formal arguments
    character(*),            intent(in)    :: grotop_filename
    type(s_grotop),          intent(inout) :: grotop

    ! local variables
    integer                  :: file, nstr
    character(MaxLine)       :: error

    character(100), allocatable :: strs(:)


    ! parse filename string
    !

    nstr = split_num(grotop_filename)
    allocate(strs(nstr))

    call split(nstr, nstr, grotop_filename, strs)


    ! open GROMACS TOP file
    !
    file = gro_pp_open_file(strs(1), strs(1:nstr), error)
    if (file == InvalidUnitNo) &
      call error_msg('Input_Grotop> '//trim(error))


    ! read GROMACS TOP file
    !
    call read_grotop(file, grotop)


    ! close GROMACS TOP file
    !
    call gro_pp_close_file(file)


    deallocate(strs)

    return

  end subroutine input_grotop

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_grotop
  !> @brief        a driver subroutine for writing GRO TOP file
  !! @authors      NT
  !! @param[in]    grotop_filename : filename of GRO TOP file
  !! @param[in]    grotop          : structure of GRO TOP data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_grotop(grotop_filename, grotop)

    ! formal arguments
    character(*),            intent(in)    :: grotop_filename
    type(s_grotop),          intent(in)    :: grotop

    ! local variables
    integer                  :: file


    ! open GROMACS TOP file
    !
    call open_file(file, grotop_filename, IOFileOutputNew)


    ! write GROMACS TOP file
    !
    call write_grotop(file, grotop)


    ! close GROMACS TOP file
    !
    call close_file(file)

    return

  end subroutine output_grotop

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_grotop
  !> @brief        initialize GROMACS TOP information
  !! @authors      NT
  !! @param[inout] grotop : structure of GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_grotop(grotop)

    ! formal arguments
    type(s_grotop),          intent(inout) :: grotop


    grotop%num_atomtypes   = 0
    grotop%num_bondtypes   = 0
    grotop%num_angltypes   = 0
    grotop%num_dihetypes   = 0
    grotop%num_constrtypes = 0
    grotop%num_cmaptypes   = 0
    grotop%num_nbonparms   = 0
    grotop%num_moltypes    = 0
    grotop%num_molss       = 0
    grotop%alias_atomtype  = .false.

    return

  end subroutine init_grotop

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_grotop_mod
  !> @brief        initialize GROMACS TOP molecule information
  !! @authors      NT
  !! @param[inout] gromol : structure of GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_grotop_mol(gromol)

    ! formal arguments
    type(s_grotop_mol),  intent(inout) :: gromol


    gromol%num_atoms   = 0
    gromol%num_bonds   = 0
    gromol%num_angls   = 0
    gromol%num_dihes   = 0
    gromol%num_cmaps   = 0
    gromol%num_excls   = 0
    gromol%num_constrs = 0
    gromol%num_pairs   = 0
    gromol%num_vsites2 = 0
    gromol%num_vsites3 = 0
    gromol%num_vsites4 = 0
    gromol%num_vsitesn = 0
    gromol%num_posress = 0

    return

  end subroutine init_grotop_mol

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    realloc_grotop
  !> @brief        reallocate GROMACS TOP information
  !! @authors      NT
  !! @param[inout] grotop   : structure of GROMACS TOP information
  !! @param[in]    variable : selected variable
  !! @param[in]    var_size : size of the selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine realloc_grotop(grotop, variable, var_size)

    ! formal arguments
    type(s_grotop),          intent(inout) :: grotop
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size

    ! local variables
    integer                  :: alloc_stat
    integer                  :: old_size

    type(s_atomtype),        allocatable :: atomtypes(:)
    type(s_bondtype),        allocatable :: bondtypes(:)
    type(s_angltype),        allocatable :: angltypes(:)
    type(s_dihetype),        allocatable :: dihetypes(:)
    type(s_constrtype),      allocatable :: constrtypes(:)
    type(s_cmaptype),        allocatable :: cmaptypes(:)
    type(s_nbonparm),        allocatable :: nbonparms(:)
    type(s_moltype),         allocatable :: moltypes(:)
    type(s_mols),            allocatable :: molss(:)


    select case(variable)

    case(GroTopAtomType)
      if (allocated(grotop%atomtypes)) then

        old_size = size(grotop%atomtypes)
        if (old_size == var_size) &
             return
        allocate(atomtypes(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        atomtypes(:) = grotop%atomtypes(:)
        deallocate(grotop%atomtypes, stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_dealloc
        allocate(grotop%atomtypes(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        if (old_size > var_size) then
          grotop%atomtypes(1:var_size) = atomtypes(1:var_size)
        else
          grotop%atomtypes(1:old_size) = atomtypes(1:old_size)
        end if
        deallocate(atomtypes, stat = alloc_stat)
      else
        old_size = 0
        allocate(grotop%atomtypes(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
      end if

    case(GroTopBondType)
      if (allocated(grotop%bondtypes)) then

        old_size = size(grotop%bondtypes)
        if (old_size == var_size) &
          return
        allocate(bondtypes(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        bondtypes(:) = grotop%bondtypes(:)
        deallocate(grotop%bondtypes, stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_dealloc
        allocate(grotop%bondtypes(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        if (old_size > var_size) then
          grotop%bondtypes(1:var_size) = bondtypes(1:var_size)
        else
          grotop%bondtypes(1:old_size) = bondtypes(1:old_size)
        end if
        deallocate(bondtypes, stat = alloc_stat)
      else
        old_size = 0
        allocate(grotop%bondtypes(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
      end if

    case(GroTopAnglType)
      if (allocated(grotop%angltypes)) then

        old_size = size(grotop%angltypes)
        if (old_size == var_size) &
          return
        allocate(angltypes(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        angltypes(:) = grotop%angltypes(:)
        deallocate(grotop%angltypes, stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_dealloc
        allocate(grotop%angltypes(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        if (old_size > var_size) then
          grotop%angltypes(1:var_size) = angltypes(1:var_size)
        else
          grotop%angltypes(1:old_size) = angltypes(1:old_size)
        end if
        deallocate(angltypes, stat = alloc_stat)
      else
        old_size = 0
        allocate(grotop%angltypes(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
      end if

    case(GroTopDiheType)
      if (allocated(grotop%dihetypes)) then

        old_size = size(grotop%dihetypes)
        if (old_size == var_size) &
          return
        allocate(dihetypes(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        dihetypes(:) = grotop%dihetypes(:)
        deallocate(grotop%dihetypes, stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_dealloc
        allocate(grotop%dihetypes(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        if (old_size > var_size) then
          grotop%dihetypes(1:var_size) = dihetypes(1:var_size)
        else
          grotop%dihetypes(1:old_size) = dihetypes(1:old_size)
        end if
        deallocate(dihetypes, stat = alloc_stat)
      else
        old_size = 0
        allocate(grotop%dihetypes(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
      end if

    case(GroTopConstrType)
      if (allocated(grotop%constrtypes)) then

        old_size = size(grotop%constrtypes)
        if (old_size == var_size) &
          return
        allocate(constrtypes(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        constrtypes(:) = grotop%constrtypes(:)
        deallocate(grotop%constrtypes, stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_dealloc
        allocate(grotop%constrtypes(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        if (old_size > var_size) then
          grotop%constrtypes(1:var_size) = constrtypes(1:var_size)
        else
          grotop%constrtypes(1:old_size) = constrtypes(1:old_size)
        end if
        deallocate(constrtypes, stat = alloc_stat)
      else
        old_size = 0
        allocate(grotop%constrtypes(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
      end if

    case(GroTopCmapType)
      if (allocated(grotop%cmaptypes)) then

        old_size = size(grotop%cmaptypes)
        if (old_size == var_size) &
          return
        allocate(cmaptypes(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        cmaptypes(:) = grotop%cmaptypes(:)
        deallocate(grotop%cmaptypes, stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_dealloc
        allocate(grotop%cmaptypes(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        if (old_size > var_size) then
          grotop%cmaptypes(1:var_size) = cmaptypes(1:var_size)
        else
          grotop%cmaptypes(1:old_size) = cmaptypes(1:old_size)
        end if
        deallocate(cmaptypes, stat = alloc_stat)
      else
        old_size = 0
        allocate(grotop%cmaptypes(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
      end if

    case(GroTopNbonParm)
      if (allocated(grotop%nbonparms)) then

        old_size = size(grotop%nbonparms)
        if (old_size == var_size) &
          return
        allocate(nbonparms(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        nbonparms(:) = grotop%nbonparms(:)
        deallocate(grotop%nbonparms, stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_dealloc
        allocate(grotop%nbonparms(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        if (old_size > var_size) then
          grotop%nbonparms(1:var_size) = nbonparms(1:var_size)
        else
          grotop%nbonparms(1:old_size) = nbonparms(1:old_size)
        end if
        deallocate(nbonparms, stat = alloc_stat)
      else
        old_size = 0
        allocate(grotop%nbonparms(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
      end if

    case(GroTopMolType)
      if (allocated(grotop%moltypes)) then

        old_size = size(grotop%moltypes)
        if (old_size == var_size) &
          return
        allocate(moltypes(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        moltypes(:) = grotop%moltypes(:)
        deallocate(grotop%moltypes, stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_dealloc
        allocate(grotop%moltypes(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        if (old_size > var_size) then
          grotop%moltypes(1:var_size) = moltypes(1:var_size)
        else
          grotop%moltypes(1:old_size) = moltypes(1:old_size)
        end if
        deallocate(moltypes, stat = alloc_stat)
      else
        old_size = 0
        allocate(grotop%moltypes(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
      end if

    case(GroTopMols)
      if (allocated(grotop%molss)) then

        old_size = size(grotop%molss)
        if (old_size == var_size) &
          return
        allocate(molss(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        molss(:) = grotop%molss(:)
        deallocate(grotop%molss, stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_dealloc
        allocate(grotop%molss(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        if (old_size > var_size) then
          grotop%molss(1:var_size) = molss(1:var_size)
        else
          grotop%molss(1:old_size) = molss(1:old_size)
        end if
        deallocate(molss, stat = alloc_stat)
      else
        old_size = 0
        allocate(grotop%molss(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
      end if

    case default

      call error_msg('Realloc_Grotop> bad variable')

    end select

    return

  end subroutine realloc_grotop

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    realloc_grotop_mol
  !> @brief        reallocate GROMACS TOP molecule information
  !! @authors      NT
  !! @param[inout] gromol   : structure of GROMACS TOP molecule information
  !! @param[in]    variable : selected variable
  !! @param[in]    var_size : size of the selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine realloc_grotop_mol(gromol, variable, var_size)

    ! formal arguments
    type(s_grotop_mol),      intent(inout) :: gromol
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size

    ! local variables
    integer                  :: alloc_stat
    integer                  :: old_size

    type(s_atom),            allocatable   :: atoms(:)
    type(s_bond),            allocatable   :: bonds(:)
    type(s_angl),            allocatable   :: angls(:)
    type(s_dihe),            allocatable   :: dihes(:)
    type(s_cmap),            allocatable   :: cmaps(:)
    type(s_excl),            allocatable   :: excls(:)
    type(s_constr),          allocatable   :: constrs(:)
    type(s_pair),            allocatable   :: pairs(:)
    type(s_vsite2),          allocatable   :: vsites2(:)
    type(s_vsite3),          allocatable   :: vsites3(:)
    type(s_vsite4),          allocatable   :: vsites4(:)
    type(s_vsiten),          allocatable   :: vsitesn(:)
    type(s_posres),          allocatable   :: posress(:)


    select case(variable)

    case(GMGroMolAtom)
      if (allocated(gromol%atoms)) then

        old_size = size(gromol%atoms)
        if (old_size == var_size) &
          return
        allocate(atoms(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        atoms(:) = gromol%atoms(:)
        deallocate(gromol%atoms, stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_dealloc
        allocate(gromol%atoms(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        if (old_size > var_size) then
          gromol%atoms(1:var_size) = atoms(1:var_size)
        else
          gromol%atoms(1:old_size) = atoms(1:old_size)
        end if
        deallocate(atoms, stat = alloc_stat)
      else
        old_size = 0
        allocate(gromol%atoms(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
      end if

    case(GMGroMolBond)
      if (allocated(gromol%bonds)) then

        old_size = size(gromol%bonds)
        if (old_size == var_size) &
          return
        allocate(bonds(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        bonds(:) = gromol%bonds(:)
        deallocate(gromol%bonds, stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_dealloc
        allocate(gromol%bonds(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        if (old_size > var_size) then
          gromol%bonds(1:var_size) = bonds(1:var_size)
        else
          gromol%bonds(1:old_size) = bonds(1:old_size)
        end if
        deallocate(bonds, stat = alloc_stat)
      else
        old_size = 0
        allocate(gromol%bonds(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
      end if

    case(GMGroMolAngl)
      if (allocated(gromol%angls)) then

        old_size = size(gromol%angls)
        if (old_size == var_size) &
          return
        allocate(angls(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        angls(:) = gromol%angls(:)
        deallocate(gromol%angls, stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_dealloc
        allocate(gromol%angls(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        if (old_size > var_size) then
          gromol%angls(1:var_size) = angls(1:var_size)
        else
          gromol%angls(1:old_size) = angls(1:old_size)
        end if
        deallocate(angls, stat = alloc_stat)
      else
        old_size = 0
        allocate(gromol%angls(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
      end if

    case(GMGroMolDihe)
      if (allocated(gromol%dihes)) then

        old_size = size(gromol%dihes)
        if (old_size == var_size) &
          return
        allocate(dihes(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        dihes(:) = gromol%dihes(:)
        deallocate(gromol%dihes, stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_dealloc
        allocate(gromol%dihes(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        if (old_size > var_size) then
          gromol%dihes(1:var_size) = dihes(1:var_size)
        else
          gromol%dihes(1:old_size) = dihes(1:old_size)
        end if
        deallocate(dihes, stat = alloc_stat)
      else
        old_size = 0
        allocate(gromol%dihes(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
      end if

    case(GMGroMolCmap)
      if (allocated(gromol%cmaps)) then

        old_size = size(gromol%cmaps)
        if (old_size == var_size) &
          return
        allocate(cmaps(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        cmaps(:) = gromol%cmaps(:)
        deallocate(gromol%cmaps, stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_dealloc
        allocate(gromol%cmaps(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        if (old_size > var_size) then
          gromol%cmaps(1:var_size) = cmaps(1:var_size)
        else
          gromol%cmaps(1:old_size) = cmaps(1:old_size)
        end if
        deallocate(cmaps, stat = alloc_stat)
      else
        old_size = 0
        allocate(gromol%cmaps(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
      end if

    case(GMGroMolExcl)
      if (allocated(gromol%excls)) then

        old_size = size(gromol%excls)
        if (old_size == var_size) &
          return
        allocate(excls(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        excls(:) = gromol%excls(:)
        deallocate(gromol%excls, stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_dealloc
        allocate(gromol%excls(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        if (old_size > var_size) then
          gromol%excls(1:var_size) = excls(1:var_size)
        else
          gromol%excls(1:old_size) = excls(1:old_size)
        end if
        deallocate(excls, stat = alloc_stat)
      else
        old_size = 0
        allocate(gromol%excls(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
      end if

    case(GMGroMolConstr)
      if (allocated(gromol%constrs)) then

        old_size = size(gromol%constrs)
        if (old_size == var_size) &
          return
        allocate(constrs(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        constrs(:) = gromol%constrs(:)
        deallocate(gromol%constrs, stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_dealloc
        allocate(gromol%constrs(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        if (old_size > var_size) then
          gromol%constrs(1:var_size) = constrs(1:var_size)
        else
          gromol%constrs(1:old_size) = constrs(1:old_size)
        end if
        deallocate(constrs, stat = alloc_stat)
      else
        old_size = 0
        allocate(gromol%constrs(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
      end if

    case(GMGroMolPair)
      if (allocated(gromol%pairs)) then

        old_size = size(gromol%pairs)
        if (old_size == var_size) &
          return
        allocate(pairs(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        pairs(:) = gromol%pairs(:)
        deallocate(gromol%pairs, stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_dealloc
        allocate(gromol%pairs(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        if (old_size > var_size) then
          gromol%pairs(1:var_size) = pairs(1:var_size)
        else
          gromol%pairs(1:old_size) = pairs(1:old_size)
        end if
        deallocate(pairs, stat = alloc_stat)
      else
        old_size = 0
        allocate(gromol%pairs(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
      end if

    case(GMGroMolVsites2)
      if (allocated(gromol%vsites2)) then

        old_size = size(gromol%vsites2)
        if (old_size == var_size) &
          return
        allocate(vsites2(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        vsites2(:) = gromol%vsites2(:)
        deallocate(gromol%vsites2, stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_dealloc
        allocate(gromol%vsites2(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        if (old_size > var_size) then
          gromol%vsites2(1:var_size) = vsites2(1:var_size)
        else
          gromol%vsites2(1:old_size) = vsites2(1:old_size)
        end if
        deallocate(vsites2, stat = alloc_stat)
      else
        old_size = 0
        allocate(gromol%vsites2(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
      end if

    case(GMGroMolVsites3)
      if (allocated(gromol%vsites3)) then

        old_size = size(gromol%vsites3)
        if (old_size == var_size) &
          return
        allocate(vsites3(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        vsites3(:) = gromol%vsites3(:)
        deallocate(gromol%vsites3, stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_dealloc
        allocate(gromol%vsites3(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        if (old_size > var_size) then
          gromol%vsites3(1:var_size) = vsites3(1:var_size)
        else
          gromol%vsites3(1:old_size) = vsites3(1:old_size)
        end if
        deallocate(vsites3, stat = alloc_stat)
      else
        old_size = 0
        allocate(gromol%vsites3(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
      end if

    case(GMGroMolVsites4)
      if (allocated(gromol%vsites4)) then

        old_size = size(gromol%vsites4)
        if (old_size == var_size) &
          return
        allocate(vsites4(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        vsites4(:) = gromol%vsites4(:)
        deallocate(gromol%vsites4, stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_dealloc
        allocate(gromol%vsites4(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        if (old_size > var_size) then
          gromol%vsites4(1:var_size) = vsites4(1:var_size)
        else
          gromol%vsites4(1:old_size) = vsites4(1:old_size)
        end if
        deallocate(vsites4, stat = alloc_stat)
      else
        old_size = 0
        allocate(gromol%vsites4(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
      end if

    case(GMGroMolVsitesn)
      if (allocated(gromol%vsitesn)) then

        old_size = size(gromol%vsitesn)
        if (old_size == var_size) &
          return
        allocate(vsitesn(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        vsitesn(:) = gromol%vsitesn(:)
        deallocate(gromol%vsitesn, stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_dealloc
        allocate(gromol%vsitesn(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        if (old_size > var_size) then
          gromol%vsitesn(1:var_size) = vsitesn(1:var_size)
        else
          gromol%vsitesn(1:old_size) = vsitesn(1:old_size)
        end if
        deallocate(vsitesn, stat = alloc_stat)
      else
        old_size = 0
        allocate(gromol%vsitesn(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
      end if

    case(GMGroMolPosres)
      if (allocated(gromol%posress)) then

        old_size = size(gromol%posress)
        if (old_size == var_size) &
          return
        allocate(posress(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        posress(:) = gromol%posress(:)
        deallocate(gromol%posress, stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_dealloc
        allocate(gromol%posress(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        if (old_size > var_size) then
          gromol%posress(1:var_size) = posress(1:var_size)
        else
          gromol%posress(1:old_size) = posress(1:old_size)
        end if
        deallocate(posress, stat = alloc_stat)
      else
        old_size = 0
        allocate(gromol%posress(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
      end if

    case default

      call error_msg('Realloc_Grotop_Mol> bad variable')

    end select

    return

  end subroutine realloc_grotop_mol

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_grotop_all
  !> @brief        deallocate GROMACS TOP information
  !! @authors      NT
  !! @param[inout] grotop : structure of GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_grotop_all(grotop)

    ! formal arguments
    type(s_grotop),          intent(inout) :: grotop

    ! local variables
    integer                  :: i, dealloc_stat


    if (allocated(grotop%atomtypes)) then
      deallocate(grotop%atomtypes, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(grotop%bondtypes)) then
      deallocate(grotop%bondtypes, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(grotop%angltypes)) then
      deallocate(grotop%angltypes, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(grotop%dihetypes)) then
      deallocate(grotop%dihetypes, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(grotop%constrtypes)) then
      deallocate(grotop%constrtypes, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(grotop%cmaptypes)) then
      deallocate(grotop%cmaptypes, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(grotop%nbonparms)) then
      deallocate(grotop%nbonparms, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(grotop%moltypes)) then
      do i = 1, size(grotop%moltypes)
        call dealloc_grotop_mol_all(grotop%moltypes(i)%mol)
      end do
      deallocate(grotop%moltypes, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(grotop%molss)) then
      deallocate(grotop%molss, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    return

  end subroutine dealloc_grotop_all

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_grotop_mol_all
  !> @brief        deallocate GROMACS TOP molecule information
  !! @authors      NT
  !! @param[inout] gromol : structure of GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_grotop_mol_all(gromol)

    ! formal arguments
    type(s_grotop_mol),      intent(inout) :: gromol

    ! local variables
    integer                  :: dealloc_stat


    if (allocated(gromol%atoms)) then
      deallocate(gromol%atoms, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(gromol%bonds)) then
      deallocate(gromol%bonds, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(gromol%angls)) then
      deallocate(gromol%angls, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(gromol%dihes)) then
      deallocate(gromol%dihes, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(gromol%cmaps)) then
      deallocate(gromol%cmaps, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(gromol%excls)) then
      deallocate(gromol%excls, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(gromol%constrs)) then
      deallocate(gromol%constrs, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(gromol%pairs)) then
      deallocate(gromol%pairs, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(gromol%vsites2)) then
      deallocate(gromol%vsites2, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(gromol%vsites3)) then
      deallocate(gromol%vsites3, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(gromol%vsites4)) then
      deallocate(gromol%vsites4, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(gromol%vsitesn)) then
      deallocate(gromol%vsitesn, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(gromol%posress)) then
      deallocate(gromol%posress, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    return

  end subroutine dealloc_grotop_mol_all

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_grotop
  !> @brief        read a GROMACS TOP file
  !! @authors      NT
  !! @param[in]    file   : unit number of GROTOP file
  !! @param[out]   grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_grotop(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),  target, intent(out)   :: grotop

    ! local variables
    integer                     :: directive, i, j, nmol
    character(MaxLine)          :: line
    logical                     :: unk, omit

    type(s_grotop_mol), pointer :: gromol


    directive = DUnknown
    unk = .false.

    do while(.true.)
      read(file,'(A)',end=100) line

      ! check directive
      if (.not. check_directive(line, directive)) &
        cycle

      ! read per directive
      select case(directive)

      case (DDefaults)
        call read_defaults(file, grotop)

      case (DAtomTypes)
        call read_atom_types(file, grotop)

      case (DBondTypes)
        call read_bond_types(file, grotop)

      case (DAngleTypes)
        call read_angle_types(file, grotop)

      case (DDihedralTypes)
        call read_dihedral_types(file, grotop)

      case (DConstrTypes)
        call read_constraint_types(file, grotop)

      case (DCmapTypes)
        call read_cmap_types(file, grotop)

      case (DNonbondParams)
        call read_nonbond_parms(file, grotop)

      case (DMoleculeType)
        call read_molecule_type(file, grotop, gromol)

      case (DAtoms)
        call read_atoms(file, grotop, gromol)

      case (DBonds)
        call read_bonds(file, grotop, gromol)

      case (DAngles)
        call read_angles(file, grotop, gromol)

      case (DDihedrals)
        call read_dihedrals(file, grotop, gromol)

      case (DCmap)
        call read_cmaps(file, grotop, gromol)

      case (DExclusions)
        call read_exclusions(file, gromol)

      case (DConstraints)
        call read_constraints(file, grotop, gromol)

      case (DPairs)
        call read_pairs(file, gromol)

      case (DSettles)
        call read_settles(file, gromol)

      case (DVirtualSites2)
        call read_virtual_sites2(file, gromol)

      case (DVirtualSites3)
        call read_virtual_sites3(file, gromol)

      case (DVirtualSites4)
        call read_virtual_sites4(file, gromol)

      case (DVirtualSitesN)
        call read_virtual_sitesn(file, gromol)

      case (DPosition_Restraints)
        call read_position_restraints(file, gromol)

      case (DSystem)
        call read_system()

      case (DMolecules)
        call read_molecules(file, grotop)

      case (DUnknown)
        if (main_rank) &
          write(MsgOut,*) 'Read_Grotop> INFO. Unknown directive:'// &
               line(1:index(line,']',.true.))
        unk = .true.

      end select

    end do

100 continue

    grotop%num_atomtypes   = size_grotop(grotop, GroTopAtomType)
    grotop%num_bondtypes   = size_grotop(grotop, GroTopBondType)
    grotop%num_angltypes   = size_grotop(grotop, GroTopAnglType)
    grotop%num_dihetypes   = size_grotop(grotop, GroTopDiheType)
    grotop%num_constrtypes = size_grotop(grotop, GroTopConstrType)
    grotop%num_cmaptypes   = size_grotop(grotop, GroTopCmapType)
    grotop%num_nbonparms   = size_grotop(grotop, GroTopNbonParm)
    grotop%num_moltypes    = size_grotop(grotop, GroTopMolType)
    grotop%num_molss       = size_grotop(grotop, GroTopMols)

    ! bind molecule and type
    !
    do i = 1, grotop%num_molss
      grotop%molss(i)%moltype => null()
      do j = 1, grotop%num_moltypes
        if (grotop%molss(i)%name == grotop%moltypes(j)%name) then
          grotop%molss(i)%moltype => grotop%moltypes(j)
        end if
      end do
      if (.not. associated(grotop%molss(i)%moltype)) then
        if (main_rank) &
          call error_msg('Read_Grotop> molecule is not found. '// &
                         trim(grotop%molss(i)%name))
      end if
    end do

    ! write summary of GROTOP information
    !
    if (main_rank) then
      if (unk) &
        write(MsgOut,'(A)') ''

      write(MsgOut,'(A)') 'Read_Grotop> Summy of Grotopfile'
      write(MsgOut,'(A20,I10)') &
           '  num_moltypes    = ', grotop%num_moltypes

      omit = .false.
      do i = 1, grotop%num_moltypes
        do j = 1, grotop%num_molss
          if (grotop%molss(j)%name == grotop%moltypes(i)%name) &
            exit
        end do
        if (j <= grotop%num_molss) then
          write(MsgOut,'(4X,A20,A)') &
               grotop%moltypes(i)%name(1:20), '  :'
          write(MsgOut,'(6X,A14,I10)') &
               'num_atoms   = ', grotop%moltypes(i)%mol%num_atoms
          write(MsgOut,'(6X,A14,I10)') &
               'num_bonds   = ', grotop%moltypes(i)%mol%num_bonds
          write(MsgOut,'(6X,A14,I10)') &
               'num_angls   = ', grotop%moltypes(i)%mol%num_angls
          write(MsgOut,'(6X,A14,I10)') &
               'num_dihes   = ', grotop%moltypes(i)%mol%num_dihes
          write(MsgOut,'(6X,A14,I10)') &
               'num_cmaps   = ', grotop%moltypes(i)%mol%num_cmaps
          write(MsgOut,'(6X,A14,I10)') &
               'num_excls   = ', grotop%moltypes(i)%mol%num_excls
          write(MsgOut,'(6X,A14,I10)') &
               'num_constrs = ', grotop%moltypes(i)%mol%num_constrs
          write(MsgOut,'(6X,A14,I10)') &
               'num_pairs   = ', grotop%moltypes(i)%mol%num_pairs
          write(MsgOut,'(6X,A14,I10)') &
               'num_vsites2 = ', grotop%moltypes(i)%mol%num_vsites2
          write(MsgOut,'(6X,A14,I10)') &
               'num_vsites3 = ', grotop%moltypes(i)%mol%num_vsites3
          write(MsgOut,'(6X,A14,I10)') &
               'num_vsites4 = ', grotop%moltypes(i)%mol%num_vsites4
          write(MsgOut,'(6X,A14,I10)') &
               'num_vsitesn = ', grotop%moltypes(i)%mol%num_vsitesn
          write(MsgOut,'(6X,A14,I10)') &
               'num_posress = ', grotop%moltypes(i)%mol%num_posress
        else
          omit = .true.
        end if

      end do
      write(MsgOut,'(4X,A)') '.. not used molecule types were hidden.'
      write(MsgOut,'(A)') ''

      nmol = 0
      do i = 1, grotop%num_molss
        nmol = nmol + grotop%molss(i)%count
      end do
      write(MsgOut,'(A20,I10)') &
           '  num_molecules   = ', nmol
      do i = 1, grotop%num_molss
        write(MsgOut,'(4X,A20,A,I0)') &
             grotop%molss(i)%name(1:20), '  :  ', grotop%molss(i)%count
      end do
      write(MsgOut,'(A)') ''

      write(MsgOut,'(A20,I10,A20,I10)') &
           '  num_atomtypes   = ', grotop%num_atomtypes,   &
           '  num_bondtypes   = ', grotop%num_bondtypes
      write(MsgOut,'(A20,I10,A20,I10)') &
           '  num_angltypes   = ', grotop%num_angltypes,   &
           '  num_dihetypes   = ', grotop%num_dihetypes
      write(MsgOut,'(A20,I10,A20,I10)') &
           '  num_constrtypes = ', grotop%num_constrtypes, &
           '  num_cmaptypes   = ', grotop%num_cmaptypes
      write(MsgOut,'(A20,I10)') &
           '  num_nbonparms   = ', grotop%num_nbonparms
      write(MsgOut,'(A)') ''

    end if

    return

  end subroutine read_grotop

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_defaults
  !> @brief        read section [defaults]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_defaults(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),  target, intent(inout) :: grotop

    ! local variables
    character(100)           :: line
    character(10)            :: yesno


    yesno = 'no'
    grotop%defaults%fudge_lj = 1.0_wp
    grotop%defaults%fudge_qq = 1.0_wp

    read(file,'(a)') line
    read(line,*,err=1,end=1) grotop%defaults%nonb_func,  &
                             grotop%defaults%combi_rule, &
                             yesno,                      &
                             grotop%defaults%fudge_lj,   &
                             grotop%defaults%fudge_qq

1   call tolower(yesno)
    grotop%defaults%gen_pairs = (yesno == 'yes')

    return

  end subroutine read_defaults

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_atom_types
  !> @brief        read section [atomtypes]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_atom_types(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),  target, intent(inout) :: grotop

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line

    type(s_atomtype), pointer :: at


    ! check count
    !
    cnt = check_section_count(file)

    ! allocate memory
    !
    old_cnt = size_grotop(grotop, GroTopAtomType)
    call realloc_grotop(grotop, GroTopAtomType, old_cnt+cnt)

    if (VerboseOn .and. main_rank) &
      write(MsgOut,'(" Read_Grotop> [atomtypes] :"'// &
                   ',i5," (total:",i5,")")') cnt, old_cnt+cnt

    ! read data
    !
    do i = 1, cnt
      at => grotop%atomtypes(old_cnt+i)

      read(file,'(A)',err=900,end=900) line

      if (match_format(line, 'SNNSNNNN')) then

        at%type_no = 0
        read(line,*)       &
             at%type_name, &
             at%mass,      &
             at%charge,    &
             at%ptype,     &
             at%v,         &
             at%w,         &
             at%khh,       &
             at%stokes_r

      else if (match_format(line, 'SSNNNSNN')) then

        read(line,*)       &
             at%alias,     &
             at%type_name, &
             at%type_no,   &
             at%mass,      &
             at%charge,    &
             at%ptype,     &
             at%v,         &
             at%w

        if (.not. grotop%alias_atomtype) then
          if (VerboseOn .and. main_rank) &
            write(MsgOut,*) &
                 'Read_Grotop> [atomtypes] : Found alias name. (column:0)'
          grotop%alias_atomtype = .true.
        end if

      else if (match_format(line, 'SNNNSNN')) then

        read(line,*)       &
             at%type_name, &
             at%type_no,   &
             at%mass,      &
             at%charge,    &
             at%ptype,     &
             at%v,         &
             at%w

      else if (match_format(line, 'SSNNSNN')) then

        at%type_no = 0
        read(line,*)       &
             at%alias,     &
             at%type_name, &
             at%mass,      &
             at%charge,    &
             at%ptype,     &
             at%v,         &
             at%w

        if (.not. grotop%alias_atomtype) then
          if (VerboseOn .and. main_rank) &
            write(MsgOut,*) &
                 'Read_Grotop> [atomtypes] : Found alias name. (column:0)2'
          grotop%alias_atomtype = .true.
        end if

      else if (match_format(line, 'SNNSNN')) then

        at%type_no = 0
        read(line,*)       &
             at%type_name, &
             at%mass,      &
             at%charge,    &
             at%ptype,     &
             at%v,         &
             at%w

      else

        goto 900

      end if

    end do

    return

900 call error_msg_grotop(file, 'read_grotop> read error. [atomtypes]')

  end subroutine read_atom_types

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_bond_types
  !> @brief        read section [bondtypes]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_bond_types(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),  target, intent(inout) :: grotop

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line

    type(s_bondtype), pointer :: bt


    ! check count
    !
    cnt = check_section_count(file)

    ! allocate memory
    !
    old_cnt = size_grotop(grotop, GroTopBondType)
    call realloc_grotop(grotop, GroTopBondType, old_cnt+cnt)

    if (VerboseOn .and. main_rank) &
      write(MsgOut,'(" Read_Grotop> [bondtypes] :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt

    ! read data
    !
    do i = 1, cnt

      bt => grotop%bondtypes(old_cnt+i)

      read(file,'(A)',err=900,end=900) line

      if (match_format(line, 'SSNNN')) then

        read(line,*)      &
           bt%atom_type1, &
           bt%atom_type2, &
           bt%func,       &
           bt%b0,         &
           bt%kb

        if (bt%func /= 1 .and. bt%func /= 2 .and. main_rank)     &
          write(MsgOut,'(" Read_Grotop> WARNING: [bondtypes] '// &
             'not supported function type. [", i3, "]")') bt%func

      else

        goto 900

      end if

    end do

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [bondtypes]')

  end subroutine read_bond_types

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_angle_types
  !> @brief        read section [angletypes]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_angle_types(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),  target, intent(inout) :: grotop

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line

    type(s_angltype), pointer :: at


    ! check count
    !
    cnt = check_section_count(file)

    ! allocate memory
    !
    old_cnt = size_grotop(grotop, GroTopAnglType)
    call realloc_grotop(grotop, GroTopAnglType, old_cnt+cnt)

    if (VerboseOn .and. main_rank) &
      write(MsgOut,'(" Read_Grotop> [angletypes] :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt

    ! read data
    !
    do i = 1, cnt

      at => grotop%angltypes(old_cnt+i)

      read(file,'(A)',err=900,end=900) line

      if (match_format(line, 'SSSNNNNN')) then

        read(line,*)      &
           at%atom_type1, &
           at%atom_type2, &
           at%atom_type3, &
           at%func,       &
           at%theta_0,    &
           at%kt,         &
           at%r13,        &
           at%kub

        if (at%func /= 5 .and. main_rank)      &
          write(MsgOut,'(" Read_Grotop> WARNING: [angletypes] '// &
             'not supported function type. [", i3, "]")') at%func

      else if (match_format(line, 'SSSNNN')) then

        read(line,*)      &
           at%atom_type1, &
           at%atom_type2, &
           at%atom_type3, &
           at%func,       &
           at%theta_0,    &
           at%kt

        if (at%func /= 1 .and. at%func /= 2 .and. main_rank)      &
          write(MsgOut,'(" Read_Grotop> WARNING: [angletypes] '// &
             'not supported function type. [", i3, "]")') at%func

      else

        goto 900

      end if

    end do

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [angletypes]')

  end subroutine read_angle_types

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_dihedral_types
  !> @brief        read section [dihedraltypes]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_dihedral_types(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),  target, intent(inout) :: grotop

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line

    type(s_dihetype), pointer :: dt


    ! check count
    !
    cnt = check_section_count(file)

    ! allocate memory
    !
    old_cnt = size_grotop(grotop, GroTopDiheType)

    call realloc_grotop(grotop, GroTopDiheType, old_cnt+cnt)

    if (VerboseOn .and. main_rank) then
      write(MsgOut,'(" Read_Grotop> [dihedraltypes] proper   :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt
    end if

    ! read data
    !
    do i = 1, cnt

      dt => grotop%dihetypes(old_cnt+i)

      read(file,'(A)',err=900,end=900) line

      if (match_format(line, 'SS1NNN') .or. &
          match_format(line, 'SS4NNN') .or. &
          match_format(line, 'SS9NNN')) then

        ! proper dihedral (two atom type)
        dt%atom_type1 = ''
        dt%atom_type4 = ''

        read(line,*) &
             dt%atom_type2, &
             dt%atom_type3, &
             dt%func, &
             dt%ps,   &
             dt%kp,   &
             dt%multiplicity
        dt%four_type = .false.

      else if (match_format(line, 'SS2NN')) then

        ! improper dihedral (two atom type)
        dt%atom_type2 = ''
        dt%atom_type3 = ''

        read(line,*) &
             dt%atom_type1, &
             dt%atom_type4, &
             dt%func, &
             dt%ps,   &
             dt%kp
        dt%four_type = .false.

      else if (match_format(line, 'SS3NNNNNN')) then

        ! Ryckaert-Bellemans dihedral (two atom type)
        dt%atom_type1 = ''
        dt%atom_type4 = ''

        read(line,*) &
             dt%atom_type2, &
             dt%atom_type3, &
             dt%func, &
             dt%c(1:6)
        dt%four_type = .false.

      else if (match_format(line, 'SSN')) then

        ! unknown (two atom type)
        read(line,*) &
             dt%atom_type1, &
             dt%atom_type2, &
             dt%func

        write(MsgOut,'(" Read_Grotop> WARNING: [dihedraltypes] '// &
             'not supported function type. [", i3, "]")') dt%func

      else if (match_format(line, 'SSSS1NNN') .or. &
               match_format(line, 'SSSS4NNN') .or. &
               match_format(line, 'SSSS9NNN')) then

        ! proper dihedral (four atom type)
        read(line,*) &
             dt%atom_type1, &
             dt%atom_type2, &
             dt%atom_type3, &
             dt%atom_type4, &
             dt%func, &
             dt%ps,   &
             dt%kp,   &
             dt%multiplicity
        dt%four_type = .true.

      else if (match_format(line, 'SSSS2NN')) then

        ! improper dihedral (four atom type)
        read(line,*) &
             dt%atom_type1, &
             dt%atom_type2, &
             dt%atom_type3, &
             dt%atom_type4, &
             dt%func, &
             dt%ps, &
             dt%kp
        dt%four_type = .true.

      else if (match_format(line, 'SSSS3NNNNNN')) then

        ! Ryckaert-Bellemans dihedral (four atom type)
        read(line,*) &
             dt%atom_type1, &
             dt%atom_type2, &
             dt%atom_type3, &
             dt%atom_type4, &
             dt%func, &
             dt%c(1:6)
        dt%four_type = .true.

      else if (match_format(line, 'SSSSN')) then

        ! unknown dihedral (four atom type)
        read(line,*) &
             dt%atom_type1, &
             dt%atom_type2, &
             dt%atom_type3, &
             dt%atom_type4, &
             dt%func

        write(MsgOut,'(" Read_Grotop> WARNING: [dihedraltypes] '// &
             'not supported function type. [", i3, "]")') dt%func

      else

        goto 900

      end if

    end do

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [dihedraltypes]')

  end subroutine read_dihedral_types

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_constraint_types
  !> @brief        read section [constrainttypes]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_constraint_types(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),  target, intent(inout) :: grotop

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line

    type(s_constrtype), pointer :: ct


    ! allocate memory
    !
    cnt     = check_section_count(file)

    old_cnt = size_grotop(grotop, GroTopConstrType)
    call realloc_grotop(grotop, GroTopConstrType, old_cnt+cnt)

    if (VerboseOn .and. main_rank) &
      write(MsgOut,'(" Read_Grotop> [constrainttypes] :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt

    ! read data
    !
    do i = 1, cnt

      ct => grotop%constrtypes(old_cnt+i)

      read(file,'(A)',err=900,end=900) line

      if (match_format(line, 'SSNN')) then

        read(line,*) &
             ct%atom_type1, &
             ct%atom_type2, &
             ct%func, &
             ct%b0

        if (ct%func /= 1 .and. ct%func /= 2 .and. main_rank) &
          write(MsgOut,'(" Read_Grotop> WARNING: [bondtypes] '// &
               'not supported function type. [", i3, "]")') ct%func

      else

        goto 900

      end if

    end do

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [constrainttypes]')

  end subroutine read_constraint_types

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_cmap_types
  !> @brief        read section [cmaptypes]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_cmap_types(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),  target, intent(inout) :: grotop

    ! local variables
    integer                  :: i, j, k, cnt, old_cnt, nrow, ncol
    character(10000)         :: line

    type(s_cmaptype), pointer :: ct


    ! check count
    !
    cnt = check_section_count(file)

    ! allocate memory
    !
    old_cnt = size_grotop(grotop, GroTopCmapType)
    call realloc_grotop(grotop, GroTopCmapType, old_cnt+cnt)

    if (VerboseOn .and. main_rank) &
      write(MsgOut,'(" Read_Grotop> [cmaptypes] :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt

    ! read data
    !
    do i = 1, cnt

      ct => grotop%cmaptypes(old_cnt+i)

      read(file,'(A)',err=900,end=900) line

      if (match_format(line, 'SSSSSNNN')) then

        read(line,*)      &
           ct%atom_type1, &
           ct%atom_type2, &
           ct%atom_type3, &
           ct%atom_type4, &
           ct%atom_type5, &
           ct%func,       &
           nrow, ncol

        allocate(ct%map(ncol, nrow))

        read(line,*)      &
           ct%atom_type1, &
           ct%atom_type2, &
           ct%atom_type3, &
           ct%atom_type4, &
           ct%atom_type5, &
           ct%func,       &
           nrow, ncol,    &
           ((ct%map(k,j),k=1,ncol),j=1,nrow)

        if (ct%func /= 1 .and. main_rank)      &
          write(MsgOut,'(" Read_Grotop> WARNING: [cmaptypes] '// &
             'not supported function type. [", i3, "]")') ct%func

      else

        goto 900

      end if

    end do

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [cmaptypes]')

    return

  end subroutine read_cmap_types

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_nonbond_parms
  !> @brief        read section [nonbond_parms]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_nonbond_parms(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),  target, intent(inout) :: grotop

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line

    type(s_nbonparm), pointer :: np


    ! check count
    !
    cnt = check_section_count(file)

    ! allocate memory
    !
    old_cnt = size_grotop(grotop, GroTopNbonParm)
    call realloc_grotop(grotop, GroTopNbonParm, old_cnt+cnt)

    if (VerboseOn .and. main_rank) &
      write(MsgOut,'(" Read_Grotop> [nonbond_params] :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt

    ! read data
    !
    do i = 1, cnt

      np => grotop%nbonparms(old_cnt+i)

      read(file,'(A)',err=900,end=900) line

      if (match_format(line, 'SSNNN')) then

        read(line,*) &
             np%atom_type1, &
             np%atom_type2, &
             np%func, &
             np%v,    &
             np%w

        if (np%func /= 1 .and. main_rank) &
          write(MsgOut,'(" Read_Grotop> WARNING: [nonbond_parms] '// &
             'not supported function type. [", i3, "]")') np%func

      else

        goto 900

      end if

    end do

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [nonbond_parms]')

  end subroutine read_nonbond_parms

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_molecule_type
  !> @brief        read section [moleculetype]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !! @param[inout] gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_molecule_type(file, grotop, gromol)

    ! formal arguments
    integer,                    intent(in)    :: file
    type(s_grotop),     target, intent(inout) :: grotop
    type(s_grotop_mol), pointer               :: gromol

    ! local variables
    integer                     :: old_cnt, alloc_stat
    type(s_moltype),   pointer  :: mt


    ! allocate memory
    !
    old_cnt = size_grotop(grotop, GroTopMolType)
    call realloc_grotop(grotop, GroTopMolType, old_cnt+1)

    ! allocate gromacs top molecule data
    allocate(gromol, stat=alloc_stat)
    if (alloc_stat /= 0) &
      goto 900

    mt => grotop%moltypes(old_cnt+1)

    mt%mol => gromol

    ! read molecule name, exclude nbon
    read(file,*,err=900,end=900) mt%name, mt%exclude_nbon

    if (VerboseOn .and. main_rank) &
      write(MsgOut,'(" Read_Grotop> [moleculetype] "'// &
                   ',i5,": ",a20)') old_cnt+1, mt%name
    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [moleculetype]')

  end subroutine read_molecule_type

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_atoms
  !> @brief        read section [atoms]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !! @param[inout] gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_atoms(file, grotop, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(inout) :: grotop
    type(s_grotop_mol),      pointer       :: gromol

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line
    character(20)            :: read_name

    type(s_atom),    pointer :: at


    ! check count
    !
    cnt = check_section_count(file)

    ! allocate memory
    !
    old_cnt = size_grotop_mol(gromol, GMGroMolAtom)
    call realloc_grotop_mol(gromol, GMGroMolAtom, old_cnt+cnt)

    if (VerboseOn .and. main_rank) &
      write(MsgOut,'(" Read_Grotop>   [atoms] :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt

    ! read data
    !
    do i = 1, cnt

      at => gromol%atoms(old_cnt+i)

      read(file,'(A)',err=900,end=900) line

      if (match_format(line, 'NSNSSNNNN')) then

        read(line,*)          &
             at%atom_no,      &
             read_name,       &
             at%residue_no,   &
             at%residue_name, &
             at%atom_name,    &
             at%charge_group, &
             at%charge,       &
             at%mass,         &
             at%stokes_r

        call rename_atom_type_name(grotop, read_name, at%atom_type)

      else if (match_format(line, 'NSNSSNNN')) then

        read(line,*)          &
             at%atom_no,      &
             read_name,       &
             at%residue_no,   &
             at%residue_name, &
             at%atom_name,    &
             at%charge_group, &
             at%charge,       &
             at%mass

        call rename_atom_type_name(grotop, read_name, at%atom_type)

        if (.not.search_stokes_r(grotop%atomtypes, at%atom_type, at%stokes_r)) &
          call error_msg_grotop(file, 'Read_Grotop> [atoms] '// &
               'atom type not found: '//at%atom_type)

      else if (match_format(line, 'NSNSSNN')) then

        read(line,*)          &
             at%atom_no,      &
             read_name,       &
             at%residue_no,   &
             at%residue_name, &
             at%atom_name,    &
             at%charge_group, &
             at%charge

        call rename_atom_type_name(grotop, read_name, at%atom_type)

        if (.not. search_atom_mass(grotop%atomtypes, at%atom_type, at%mass)) &
          call error_msg_grotop(file, 'Read_Grotop> [atoms] '// &
               'atom type not found: '//at%atom_type)

        if (.not.search_stokes_r(grotop%atomtypes, at%atom_type, at%stokes_r)) &
          call error_msg_grotop(file, 'Read_Grotop> [atoms] '// &
               'atom type not found: '//at%atom_type)

      else

        goto 900

      end if

    end do

    gromol%num_atoms = size_grotop_mol(gromol, GMGroMolAtom)

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [atoms]')

  end subroutine read_atoms

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_bonds
  !> @brief        read section [bonds]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !! @param[inout] gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_bonds(file, grotop, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(inout) :: grotop
    type(s_grotop_mol),      pointer       :: gromol

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line
    character(6)             :: at1, at2

    type(s_bond),    pointer :: bd


    ! check count
    !
    cnt = check_section_count(file)

    ! allocate memory
    !
    old_cnt = size_grotop_mol(gromol, GMGroMolBond)
    call realloc_grotop_mol(gromol, GMGroMolBond, old_cnt+cnt)

    if (VerboseOn .and. main_rank) &
      write(MsgOut,'(" Read_Grotop>   [bonds] :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt

    ! read data
    !
    do i = 1, cnt

      bd => gromol%bonds(old_cnt+i)

      read(file,'(A)',err=900,end=900) line

      if (match_format(line, 'NN1NN') .or. &
          match_format(line, 'NN2NN')) then

        read(line,*)       &
             bd%atom_idx1, &
             bd%atom_idx2, &
             bd%func,      &
             bd%b0,        &
             bd%kb

      else if (match_format(line, 'NN1') .or. &
               match_format(line, 'NN2')) then

        read(line,*)       &
             bd%atom_idx1, &
             bd%atom_idx2, &
             bd%func

        at1 = gromol%atoms(bd%atom_idx1)%atom_type
        at2 = gromol%atoms(bd%atom_idx2)%atom_type

        if (.not.search_bond_param(grotop%bondtypes,    &
                      at1, at2, bd%func, bd%b0, bd%kb)) &
          call error_msg_grotop(file, 'Read_Grotop> [bonds] '// &
             'atom type not found: '//at1//at2)

      else if (match_format(line, 'NNN')) then

        read(line,*) &
             bd%atom_idx1, &
             bd%atom_idx1, &
             bd%func

        if (main_rank) &
          write(MsgOut,'(" Read_Grotop> WARNING: [bonds] '// &
              'not supported function type. [", i3, "]")') bd%func

      else

        goto 900

      end if

    end do

    gromol%num_bonds = size_grotop_mol(gromol, GMGroMolBond)

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [bonds]')

  end subroutine read_bonds

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_angles
  !> @brief        read section [angles]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !! @param[inout] gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_angles(file, grotop, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(inout) :: grotop
    type(s_grotop_mol),      pointer       :: gromol

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line
    character(6)             :: at1, at2, at3

    type(s_angl), pointer    :: an


    ! check count
    !
    cnt = check_section_count(file)

    ! allocate memory
    !
    old_cnt = size_grotop_mol(gromol, GMGroMolAngl)
    call realloc_grotop_mol(gromol, GMGroMolAngl, old_cnt+cnt)

    if (VerboseOn .and. main_rank) &
      write(MsgOut,'(" Read_Grotop>   [angles] :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt

    ! read data
    !
    do i = 1, cnt

      an => gromol%angls(old_cnt+i)

      read(file,'(A)',err=900,end=900) line

      if (match_format(line, 'NNN5NNNN')) then

        read(line,*)       &
             an%atom_idx1, &
             an%atom_idx2, &
             an%atom_idx3, &
             an%func,      &
             an%theta_0,   &
             an%kt,        &
             an%r13,       &
             an%kub

      else if (match_format(line, 'NNN5')) then

        read(line,*)       &
             an%atom_idx1, &
             an%atom_idx2, &
             an%atom_idx3, &
             an%func

        at1 = gromol%atoms(an%atom_idx1)%atom_type
        at2 = gromol%atoms(an%atom_idx2)%atom_type
        at3 = gromol%atoms(an%atom_idx3)%atom_type

        if (.not.search_angl_param(grotop%angltypes, &
             at1, at2, at3, an%func, an%theta_0, an%kt, an%r13, an%kub)) &
          call error_msg_grotop(file, 'Read_Grotop> [angles] '// &
               'atom type not found: '//at1//at2//at3)

      else if (match_format(line, 'NNN1NN') .or. &
               match_format(line, 'NNN2NN')) then

        read(line,*)       &
             an%atom_idx1, &
             an%atom_idx2, &
             an%atom_idx3, &
             an%func,      &
             an%theta_0,   &
             an%kt

      else if (match_format(line, 'NNN1') .or. &
               match_format(line, 'NNN2')) then

        read(line,*)       &
             an%atom_idx1, &
             an%atom_idx2, &
             an%atom_idx3, &
             an%func

        at1 = gromol%atoms(an%atom_idx1)%atom_type
        at2 = gromol%atoms(an%atom_idx2)%atom_type
        at3 = gromol%atoms(an%atom_idx3)%atom_type

        if (.not.search_angl_param(grotop%angltypes,     &
             at1, at2, at3, an%func, an%theta_0, an%kt)) &
          call error_msg_grotop(file, 'Read_Grotop> [angles] '// &
               'atom type not found: '//at1//at2//at3)

      else if (match_format(line, 'NNNN')) then

        read(line,*)       &
             an%atom_idx1, &
             an%atom_idx2, &
             an%atom_idx3, &
             an%func

        if (main_rank) &
          write(MsgOut,'(" Read_Grotop> WARNING: [angles] '// &
                'not supported function type. [", i3, "]")') an%func

      else

        goto 900

      end if

    end do

    gromol%num_angls = size_grotop_mol(gromol, GMGroMolAngl)

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [angles]')

  end subroutine read_angles

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_dihedrals
  !> @brief        read section [dihedrals]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !! @param[inout] gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_dihedrals(file, grotop, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(inout) :: grotop
    type(s_grotop_mol),      pointer       :: gromol

    ! local variables
    integer                  :: i, cnt, old_cnt
    integer                  :: cnt_tri, old_dihedcnt, cur_size, ii
    integer                  :: iold, ix
    character(MaxLine)       :: line
    character(6)             :: at1, at2, at3, at4

    type(s_dihe),    pointer :: dh


    ! check count
    !
    cnt = check_section_count(file)

    cnt_tri      = cnt*5
    old_dihedcnt = gromol%num_dihes

    ! allocate memory
    !
    old_cnt = size_grotop_mol(gromol, GMGroMolDihe)
    if (old_cnt < old_dihedcnt+cnt_tri) then
      call realloc_grotop_mol(gromol, GMGroMolDihe, old_cnt+cnt_tri)
    end if
    cur_size = size_grotop_mol(gromol, GMGroMolDihe)

    if (VerboseOn .and. main_rank) then
      write(MsgOut,'(" Read_Grotop>   [dihedrals] proper   :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt_tri
    end if

    ! read data
    !
    ii = 0
    do i = 1, cnt

      ii = ii+1
      if (ii+5 > cur_size) then
        call realloc_grotop_mol(gromol, GMGroMolDihe, cur_size+cnt_tri)
        cur_size = size_grotop_mol(gromol, GMGroMolDihe)
      end if

      dh => gromol%dihes(old_dihedcnt+ii)

      read(file,'(A)',err=900,end=900) line

      if (match_format(line, 'NNNN1NNN')) then

        ! proper dihedral (full)
        read(line,*) &
             dh%atom_idx1, &
             dh%atom_idx2, &
             dh%atom_idx3, &
             dh%atom_idx4, &
             dh%func, &
             dh%ps,   &
             dh%kp,   &
             dh%multiplicity

      else if (match_format(line, 'NNNN1')) then

        ! proper dihedral (no parameters)
        read(line,*) &
             dh%atom_idx1, &
             dh%atom_idx2, &
             dh%atom_idx3, &
             dh%atom_idx4, &
             dh%func

        at1 = gromol%atoms(dh%atom_idx1)%atom_type
        at2 = gromol%atoms(dh%atom_idx2)%atom_type
        at3 = gromol%atoms(dh%atom_idx3)%atom_type
        at4 = gromol%atoms(dh%atom_idx4)%atom_type

        if (.not. search_dihe_param(grotop%dihetypes,        &
                                    at1, at2, at3, at4, dh)) &
          call error_msg_grotop(file, 'Read_Grotop> [proper dihedrals] '// &
               'atom type not found: '//at1//at2//at3//at4)

      else if (match_format(line, 'NNNN2NN')) then

        ! improper dihedral (full)
        read(line,*) &
             dh%atom_idx1, &
             dh%atom_idx2, &
             dh%atom_idx3, &
             dh%atom_idx4, &
             dh%func, &
             dh%ps,   &
             dh%kp

      else if (match_format(line, 'NNNN2')) then

        ! improper dihedral (no parameters)
        read(line,*) &
             dh%atom_idx1, &
             dh%atom_idx2, &
             dh%atom_idx3, &
             dh%atom_idx4, &
             dh%func

        at1 = gromol%atoms(dh%atom_idx1)%atom_type
        at2 = gromol%atoms(dh%atom_idx2)%atom_type
        at3 = gromol%atoms(dh%atom_idx3)%atom_type
        at4 = gromol%atoms(dh%atom_idx4)%atom_type

        if (.not.search_dihe_param(grotop%dihetypes,        &
                                   at1, at2, at3, at4, dh)) &
          call error_msg_grotop(file, 'Read_Grotop> [improper dih.] '// &
               'atom type not found: '//at1//at2//at3//at4)

      else if (match_format(line, 'NNNN3NNNNNN')) then

        ! Ryckaert-Bellemans dihedral (full)
        read(line,*) &
             dh%atom_idx1, &
             dh%atom_idx2, &
             dh%atom_idx3, &
             dh%atom_idx4, &
             dh%func,      &
             dh%c(1:6)

      else if (match_format(line, 'NNNN3')) then

        ! Ryckaert-Bellemans dihedral (no parameters)
        read(line,*) &
             dh%atom_idx1, &
             dh%atom_idx2, &
             dh%atom_idx3, &
             dh%atom_idx4, &
             dh%func

        at1 = gromol%atoms(dh%atom_idx1)%atom_type
        at2 = gromol%atoms(dh%atom_idx2)%atom_type
        at3 = gromol%atoms(dh%atom_idx3)%atom_type
        at4 = gromol%atoms(dh%atom_idx4)%atom_type

        if (.not. search_dihe_param(grotop%dihetypes,        &
                                    at1, at2, at3, at4, dh)) &
          call error_msg_grotop(file, 'Read_Grotop> [proper dihedrals] '// &
               'atom type not found: '//at1//at2//at3//at4)

      else if (match_format(line, 'NNNN4NNN')) then

        ! periodic improper dihedral (full)
        read(line,*) &
             dh%atom_idx1, &
             dh%atom_idx2, &
             dh%atom_idx3, &
             dh%atom_idx4, &
             dh%func, &
             dh%ps,   &
             dh%kp,   &
             dh%multiplicity

      else if (match_format(line, 'NNNN4')) then

        ! periodic improper dihedrals
        read(line,*) &
             dh%atom_idx1, &
             dh%atom_idx2, &
             dh%atom_idx3, &
             dh%atom_idx4, &
             dh%func

        at1 = gromol%atoms(dh%atom_idx1)%atom_type
        at2 = gromol%atoms(dh%atom_idx2)%atom_type
        at3 = gromol%atoms(dh%atom_idx3)%atom_type
        at4 = gromol%atoms(dh%atom_idx4)%atom_type
        iold = ii

        if (.not.search_dihe_param_duplicate(old_dihedcnt, grotop%dihetypes, &
                                   at1, at2, at3, at4, gromol%dihes, ii)) &
          call error_msg_grotop(file, 'Read_Grotop> [improper dih.] '// &
               'atom type not found: '//at1//at2//at3//at4)
        do ix = iold+1, ii
          dh => gromol%dihes(old_dihedcnt+ix)
          dh%func      = gromol%dihes(old_dihedcnt+iold)%func
          dh%atom_idx1 = gromol%dihes(old_dihedcnt+iold)%atom_idx1
          dh%atom_idx2 = gromol%dihes(old_dihedcnt+iold)%atom_idx2
          dh%atom_idx3 = gromol%dihes(old_dihedcnt+iold)%atom_idx3
          dh%atom_idx4 = gromol%dihes(old_dihedcnt+iold)%atom_idx4
        end do

      else if (match_format(line, 'NNNN9NNN')) then

        ! proper dihedral (multiple) (full)
        read(line,*) &
             dh%atom_idx1, &
             dh%atom_idx2, &
             dh%atom_idx3, &
             dh%atom_idx4, &
             dh%func, &
             dh%ps,   &
             dh%kp,   &
             dh%multiplicity

      else if (match_format(line, 'NNNN9')) then

        ! proper dihedral (multiple)
        read(line,*) &
             dh%atom_idx1, &
             dh%atom_idx2, &
             dh%atom_idx3, &
             dh%atom_idx4, &
             dh%func

        at1 = gromol%atoms(dh%atom_idx1)%atom_type
        at2 = gromol%atoms(dh%atom_idx2)%atom_type
        at3 = gromol%atoms(dh%atom_idx3)%atom_type
        at4 = gromol%atoms(dh%atom_idx4)%atom_type
        iold = ii

        if (.not.search_dihe_param_duplicate(old_dihedcnt, grotop%dihetypes, &
                                   at1, at2, at3, at4, gromol%dihes, ii)) &
          call error_msg_grotop(file, 'Read_Grotop> [proper dih.] '// &
               'atom type not found: '//at1//at2//at3//at4)
        do ix = iold+1, ii
          dh => gromol%dihes(old_dihedcnt+ix)
          dh%func      = gromol%dihes(old_dihedcnt+iold)%func
          dh%atom_idx1 = gromol%dihes(old_dihedcnt+iold)%atom_idx1
          dh%atom_idx2 = gromol%dihes(old_dihedcnt+iold)%atom_idx2
          dh%atom_idx3 = gromol%dihes(old_dihedcnt+iold)%atom_idx3
          dh%atom_idx4 = gromol%dihes(old_dihedcnt+iold)%atom_idx4
        end do

      else if (match_format(line, 'NNNNN')) then

        ! unknown
        read(line,*) &
             dh%atom_idx1, &
             dh%atom_idx2, &
             dh%atom_idx3, &
             dh%atom_idx4, &
             dh%func

        write(MsgOut,'(" Read_Grotop> WARNING: [dihedrals] '// &
              'not supported function type. [", i3, "]")') dh%func

      else

        goto 900

      end if

    end do

    gromol%num_dihes = gromol%num_dihes+ii

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [dihedrals]')

  end subroutine read_dihedrals

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_cmaps
  !> @brief        read section [cmap]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !! @param[inout] gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_cmaps(file, grotop, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(inout) :: grotop
    type(s_grotop_mol),      pointer       :: gromol

    ! local variables
    integer                  :: i, j, k, cnt, old_cnt, nrow, ncol
    character(10000)         :: line
    character(6)             :: at1, at2, at3, at4, at5

    type(s_cmap), pointer    :: cm


    ! check count
    !
    cnt = check_section_count(file)

    ! allocate memory
    !
    old_cnt = size_grotop_mol(gromol, GMGroMolCmap)
    call realloc_grotop_mol(gromol, GMGroMolCmap, old_cnt+cnt)

    if (VerboseOn .and. main_rank) &
      write(MsgOut,'(" Read_Grotop>   [cmap] :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt

    ! read data
    !
    do i = 1, cnt

      cm => gromol%cmaps(old_cnt+i)

      read(file,'(A)',err=900,end=900) line

      if (match_format(line, 'NNNNN1NN')) then

        read(line,*)       &
             cm%atom_idx1, &
             cm%atom_idx2, &
             cm%atom_idx3, &
             cm%atom_idx4, &
             cm%atom_idx5, &
             cm%func,      &
             nrow, ncol

        allocate(cm%map(ncol,nrow))

        read(line,*)       &
             cm%atom_idx1, &
             cm%atom_idx2, &
             cm%atom_idx3, &
             cm%atom_idx4, &
             cm%atom_idx5, &
             cm%func,      &
             nrow, ncol,   &
             ((cm%map(k,j),k=1,ncol),j=1,nrow)

      else if (match_format(line, 'NNNNN1')) then

        read(line,*)       &
             cm%atom_idx1, &
             cm%atom_idx2, &
             cm%atom_idx3, &
             cm%atom_idx4, &
             cm%atom_idx5, &
             cm%func

        !at1 = gromol%atoms(cm%atom_idx1)%atom_type
        !at2 = gromol%atoms(cm%atom_idx2)%atom_type
        !at3 = gromol%atoms(cm%atom_idx3)%atom_type
        !at4 = gromol%atoms(cm%atom_idx4)%atom_type
        !at5 = gromol%atoms(cm%atom_idx5)%atom_type
        !
        !if (.not.search_cmap_param(grotop%cmaptypes,   &
        !    at1, at2, at3, at4, at5, cm%func, cm%map)) &
        !  call error_msg_grotop(file, 'Read_Grotop> [cmap] '// &
        !       'atom type not found: '//at1//at2//at3//at4//at5)

      else

        goto 900

      end if

    end do

    gromol%num_cmaps = size_grotop_mol(gromol, GMGroMolCmap)

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [cmap]')

  end subroutine read_cmaps

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_exclusions
  !> @brief        read section [exclusions]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[inout] gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_exclusions(file, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop_mol),      pointer       :: gromol

    ! local variables
    integer                  :: i, j, cnt, cnt2, cnt3, old_cnt
    character(1000)          :: line

    integer,     allocatable :: v(:)


    ! check count
    !

    cnt = check_section_count(file)

    cnt2 = 0
    do i = 1, cnt
      read(file,'(A)') line
      cnt2 = cnt2 + split_num(line) - 1
    end do

    do i = 1, cnt
      backspace(file)
    end do

    ! allocate memory
    !

    old_cnt = size_grotop_mol(gromol, GMGroMolExcl)
    call realloc_grotop_mol(gromol, GMGroMolExcl, old_cnt+cnt2)

    if (VerboseOn .and. main_rank) &
      write(MsgOut,'(" Read_Grotop>   [exclusions] :"'// &
           ',i5," (total:",i5,")")') cnt2, old_cnt+cnt2

    ! read data
    !

    cnt2 = 0

    do i = 1, cnt

      read(file,'(a)') line
      cnt3 = split_num(line)
      allocate(v(cnt3))

      call split(cnt3,cnt3,line,v)

      do j = 2, cnt3
        cnt2 = cnt2 + 1
        gromol%excls(cnt2)%atom_idx1 = v(1)
        gromol%excls(cnt2)%atom_idx2 = v(j)
      end do

      deallocate(v)

    end do

    gromol%num_excls = size_grotop_mol(gromol, GMGroMolExcl)

    return

  end subroutine read_exclusions

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_constraints
  !> @brief        read section [constrants]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !! @param[inout] gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_constraints(file, grotop, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(inout) :: grotop
    type(s_grotop_mol),      pointer       :: gromol

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line
    character(6)             :: at1, at2

    type(s_constr),  pointer :: cd


    ! check count
    !
    cnt = check_section_count(file)

    ! allocate memory
    !
    old_cnt = size_grotop_mol(gromol, GMGroMolConstr)
    call realloc_grotop_mol(gromol, GMGroMolConstr, old_cnt+cnt)

    if (VerboseOn .and. main_rank) &
      write(MsgOut,'(" Read_Grotop>   [constrants] :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt

    ! read data
    !
    do i = 1, cnt

      cd => gromol%constrs(old_cnt+i)

      read(file,'(A)',err=900,end=900) line

      if (match_format(line, 'NN1N') .or. &
          match_format(line, 'NN2N')) then

        read(line,*)       &
             cd%atom_idx1, &
             cd%atom_idx2, &
             cd%func,      &
             cd%b0

      else if (match_format(line, 'NN1') .or. &
               match_format(line, 'NN2')) then

        read(line,*)       &
             cd%atom_idx1, &
             cd%atom_idx2, &
             cd%func

        at1 = gromol%atoms(cd%atom_idx1)%atom_type
        at2 = gromol%atoms(cd%atom_idx2)%atom_type

        if (.not.search_constr_param(grotop%constrtypes, &
                                     at1, at2, cd%func,cd%b0))&
          call error_msg_grotop(file, 'Read_Grotop> [constraints] '// &
             'atom type not found: '//at1//at2)

      else if (match_format(line, 'NNN')) then

        read(line,*)       &
             cd%atom_idx1, &
             cd%atom_idx2, &
             cd%func

        if (main_rank) &
          write(MsgOut,'(" Read_Grotop> WARNING: [constrants] '// &
                'not supported function type. [", i3, "]")') cd%func

      else

        goto 900

      end if

    end do

    gromol%num_constrs = size_grotop_mol(gromol, GMGroMolConstr)

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [constraints]')

  end subroutine read_constraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_pairs
  !> @brief        read section [pairs]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[inout] gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_pairs(file, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop_mol),      pointer       :: gromol

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line

    type(s_pair),    pointer :: pd


    ! check count
    !
    cnt = check_section_count(file)

    ! allocate memory
    !
    old_cnt = size_grotop_mol(gromol, GMGroMolPair)
    call realloc_grotop_mol(gromol, GMGroMolPair, old_cnt+cnt)

    if (VerboseOn .and. main_rank) &
      write(MsgOut,'(" Read_Grotop>   [pairs] :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt

    ! read data
    !
    do i = 1, cnt

      pd => gromol%pairs(old_cnt+i)

      read(file,'(A)',err=900,end=900) line

      if (match_format(line, 'NNNNNN')) then

        read(line,*) &
             pd%atom_idx1, &
             pd%atom_idx2, &
             pd%func,      &
             pd%v,         &
             pd%w,         &
             pd%khh

      else if (match_format(line, 'NN1NN') .or. &
               match_format(line, 'NN2NN')) then

        read(line,*) &
             pd%atom_idx1, &
             pd%atom_idx2, &
             pd%func,      &
             pd%v,         &
             pd%w

      else if (match_format(line, 'NNNNN')) then

        read(line,*) &
             pd%atom_idx1, &
             pd%atom_idx2, &
             pd%func,      &
             pd%v,         &
             pd%w

      else if (match_format(line, 'NN1') .or. &
               match_format(line, 'NN2')) then

        read(line,*) &
             pd%atom_idx1, &
             pd%atom_idx2, &
             pd%func

      else if (match_format(line, 'NNN')) then

        read(line,*) &
             pd%atom_idx1, &
             pd%atom_idx2, &
             pd%func

        if (main_rank) &
          write(MsgOut,'(" Read_Grotop> WARNING: [pairs] '// &
                'not supported function type. [", i3, "]")') pd%func

      else

        goto 900

      end if

    end do

    gromol%num_pairs = size_grotop_mol(gromol, GMGroMolPair)

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [pairs]')

  end subroutine read_pairs

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_settles
  !> @brief        read section [settles]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[inout] gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_settles(file, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop_mol),      pointer       :: gromol


    read(file,*,err=900,end=900) gromol%settles%ow,   &
                                 gromol%settles%func, &
                                 gromol%settles%doh,  &
                                 gromol%settles%dhh

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [settles]')

  end subroutine read_settles

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_virtual_sites2
  !> @brief        read section [virtual_sites2]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[inout] gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_virtual_sites2(file, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop_mol),      pointer       :: gromol

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line

    type(s_vsite2),  pointer :: vs2


    ! check count
    !
    cnt = check_section_count(file)

    ! allocate memory
    !
    old_cnt = size_grotop_mol(gromol, GMGroMolVsites2)
    call realloc_grotop_mol(gromol, GMGroMolVsites2, old_cnt+cnt)

    if (VerboseOn .and. main_rank) &
      write(MsgOut,'(" Read_Grotop>   [virtial_sites2] :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt

    ! read data
    !
    do i = 1, cnt

      vs2 => gromol%vsites2(old_cnt+i)

      read(file,'(A)',err=900,end=900) line

      if (match_format(line, 'NNN1N')) then

        read(line,*) &
             vs2%atom_idx1, &
             vs2%atom_idx2, &
             vs2%atom_idx3, &
             vs2%func,      &
             vs2%a

      else if (match_format(line, 'NNNN')) then

        read(line,*) &
             vs2%atom_idx1, &
             vs2%atom_idx2, &
             vs2%atom_idx3, &
             vs2%func

        if (main_rank) &
          write(MsgOut,'(" Read_Grotop> WARNING: [virtual_sites2] '// &
                'not supported function type. [", i3, "]")') vs2%func

      else

        goto 900

      end if

    end do

    gromol%num_vsites2 = size_grotop_mol(gromol, GMGroMolVsites2)

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [virtual_sites2]')

  end subroutine read_virtual_sites2

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_virtual_sites3
  !> @brief        read section [virtual_sites3]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[inout] gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_virtual_sites3(file, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop_mol),      pointer       :: gromol

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line

    type(s_vsite3),  pointer :: vs3


    ! check count
    !
    cnt = check_section_count(file)

    ! allocate memory
    !
    old_cnt = size_grotop_mol(gromol, GMGroMolVsites3)
    call realloc_grotop_mol(gromol, GMGroMolVsites3, old_cnt+cnt)

    if (VerboseOn .and. main_rank) &
      write(MsgOut,'(" Read_Grotop>   [virtial_sites3] :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt

    ! read data
    !
    do i = 1, cnt

      vs3 => gromol%vsites3(old_cnt+i)

      read(file,'(A)',err=900,end=900) line

      if (match_format(line, 'NNNN1NN')) then

        read(line,*) &
             vs3%atom_idx1, &
             vs3%atom_idx2, &
             vs3%atom_idx3, &
             vs3%atom_idx4, &
             vs3%func,      &
             vs3%a,         &
             vs3%b

      else if (match_format(line, 'NNNN2NN')) then

        read(line,*) &
             vs3%atom_idx1, &
             vs3%atom_idx2, &
             vs3%atom_idx3, &
             vs3%atom_idx4, &
             vs3%func,      &
             vs3%a,         &
             vs3%d

      else if (match_format(line, 'NNNN3NN')) then

        read(line,*) &
             vs3%atom_idx1, &
             vs3%atom_idx2, &
             vs3%atom_idx3, &
             vs3%atom_idx4, &
             vs3%func,      &
             vs3%theta,     &
             vs3%d

      else if (match_format(line, 'NNNN4NNN')) then

        read(line,*) &
             vs3%atom_idx1, &
             vs3%atom_idx2, &
             vs3%atom_idx3, &
             vs3%atom_idx4, &
             vs3%func,      &
             vs3%a,         &
             vs3%b,         &
             vs3%c

      else if (match_format(line, 'NNNNN')) then

        read(line,*) &
             vs3%atom_idx1, &
             vs3%atom_idx2, &
             vs3%atom_idx3, &
             vs3%atom_idx4, &
             vs3%func

        if (main_rank) &
          write(MsgOut,'(" Read_Grotop> WARNING: [virtual_sites3] '// &
                'not supported function type. [", i3, "]")') vs3%func

      else

        goto 900

      end if

    end do

    gromol%num_vsites3 = size_grotop_mol(gromol, GMGroMolVsites3)

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [virtual_sites3]')

  end subroutine read_virtual_sites3

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_virtual_sites4
  !> @brief        read section [virtual_sites4]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[inout] gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_virtual_sites4(file, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop_mol),      pointer       :: gromol

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line

    type(s_vsite4),  pointer :: vs4


    ! check count
    !
    cnt = check_section_count(file)

    ! allocate memory
    !
    old_cnt = size_grotop_mol(gromol, GMGroMolVsites4)
    call realloc_grotop_mol(gromol, GMGroMolVsites4, old_cnt+cnt)

    if (VerboseOn .and. main_rank) &
      write(MsgOut,'(" Read_Grotop>   [virtial_sites4] :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt

    ! read data
    !
    do i = 1, cnt

      vs4 => gromol%vsites4(old_cnt+i)

      read(file,'(A)',err=900,end=900) line

      if (match_format(line, 'NNNNN2NNN')) then

        read(line,*) &
             vs4%atom_idx1, &
             vs4%atom_idx2, &
             vs4%atom_idx3, &
             vs4%atom_idx4, &
             vs4%atom_idx5, &
             vs4%func,      &
             vs4%a,         &
             vs4%b,         &
             vs4%c

      else if (match_format(line, 'NNNNNN')) then

        read(line,*) &
             vs4%atom_idx1, &
             vs4%atom_idx2, &
             vs4%atom_idx3, &
             vs4%atom_idx4, &
             vs4%atom_idx5, &
             vs4%func

        if (main_rank) &
          write(MsgOut,'(" Read_Grotop> WARNING: [virtual_sites4] '// &
                'not supported function type. [", i3, "]")') vs4%func

      else

        goto 900

      end if

    end do

    gromol%num_vsites4 = size_grotop_mol(gromol, GMGroMolVsites4)

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [virtual_sites4]')

  end subroutine read_virtual_sites4

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_virtual_sitesn
  !> @brief        read section [virtual_sitesn]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[inout] gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_virtual_sitesn(file, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop_mol),      pointer       :: gromol

    ! constants
    integer,       parameter :: MaxFroms = 100

    ! local variables
    integer                  :: i, cnt, old_cnt, fcnt
    integer                  :: froms(MaxFroms)
    character(MaxLine)       :: line

    type(s_vsiten),  pointer :: vsn


    ! check count
    !
    cnt = check_section_count(file)

    ! allocate memory
    !
    old_cnt = size_grotop_mol(gromol, GMGroMolVsitesn)
    call realloc_grotop_mol(gromol, GMGroMolVsitesn, old_cnt+cnt)

    if (VerboseOn .and. main_rank) &
      write(MsgOut,'(" Read_Grotop>   [virtial_sitesn] :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt

    ! read data
    !
    do i = 1, cnt

      vsn => gromol%vsitesn(old_cnt+i)

      read(file,'(A)',err=900,end=900) line

      if (match_format(line, 'N1N') .or. &
          match_format(line, 'N2N') .or. &
          match_format(line, 'N3N')) then

!       froms(1:MaxFroms) = ''
        froms(1:MaxFroms) = 0
        read(line,*,err=1,end=1) &
             vsn%atom_idx, &
             vsn%func,     &
             froms(1:MaxFroms)
1       do fcnt = 1, MaxFroms
!         if (froms(fcnt) == '') &
          if (froms(fcnt) == 0) &
            exit
        end do
        fcnt = fcnt - 1

        allocate(vsn%atom_idcs(fcnt))
        vsn%atom_idcs(1:fcnt) = froms(1:fcnt)

      else if (match_format(line, 'NNN')) then

        read(line,*) &
             vsn%atom_idx, &
             vsn%func

        if (main_rank) &
          write(MsgOut,'(" Read_Grotop> WARNING: [virtual_sitesn] '// &
                'not supported function type. [", i3, "]")') vsn%func

      else

        goto 900

      end if

    end do

    gromol%num_vsitesn = size_grotop_mol(gromol, GMGroMolVsitesn)

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [virtual_sitesn]')

  end subroutine read_virtual_sitesn

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_position_restraints
  !> @brief        read section [position_restraints]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[inout] gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_position_restraints(file, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop_mol),      pointer       :: gromol

    ! local variables
    integer                  :: aidx
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line

    type(s_posres),  pointer :: pd


    ! check count
    !
    cnt = check_section_count(file)

    ! allocate memory
    !
    old_cnt = size_grotop_mol(gromol, GMGroMolPosres)
    call realloc_grotop_mol(gromol, GMGroMolPosres, old_cnt+cnt)

    if (VerboseOn .and. main_rank) &
      write(MsgOut,'(" Read_Grotop>   [position_restraints] :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt

    ! read data
    !
    do i = 1, cnt

      pd => gromol%posress(old_cnt+i)

      read(file,'(A)',err=900,end=900) line

      if (match_format(line, 'N1NNN')) then

        read(line,*) &
             pd%atom_idx, &
             pd%func,     &
             pd%kx,       &
             pd%ky,       &
             pd%kz

      else if (match_format(line, 'N1')) then

        read(line,*) &
             pd%atom_idx, &
             pd%func

      else if (match_format(line, 'NN')) then

        read(line,*) &
             pd%atom_idx, &
             pd%func

        if (main_rank) &
          write(MsgOut,'(" Read_Grotop> WARNING: [position_restraints] '// &
              'not supported function type. [", i3, "]")') pd%func

      else

        goto 900

      end if

    end do

    gromol%num_posress = size_grotop_mol(gromol, GMGroMolPosres)

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error.'//&
                         ' [position_restraints]')

  end subroutine read_position_restraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_system
  !> @brief        read section [system]
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_system()

    return

  end subroutine read_system

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_molecules
  !> @brief        read section [molecules]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_molecules(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(inout) :: grotop

    ! local variables
    integer                  :: directive, cnt, sz
    character(MaxLine)       :: line
    character(20)            :: name


    if (VerboseOn .and. main_rank) &
      write(MsgOut,'(" Read_Grotop> [molecules]")')

    do while(.true.)
      read(file,'(A)',end=100) line
      if (check_directive(line, directive)) &
        exit
      read(line,*,err=900,end=900) name, cnt

      sz = size_grotop(grotop, GroTopMols) + 1
      call realloc_grotop(grotop, GroTopMols, sz)
      grotop%molss(sz)%name = name
      grotop%molss(sz)%count = cnt

      if (VerboseOn .and. main_rank) &
        write(MsgOut,'(" Read_Grotop>   ",a20,"  count:",i5)') name, cnt
    end do

100 backspace(file)

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [molecules]')

  end subroutine read_molecules

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_grotop
  !> @brief        write a GROMACS TOP file
  !! @authors      NT
  !! @param[in]    file   : unit number of GROTOP file
  !! @param[in]    grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_grotop(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(in)    :: grotop


    ! defaults
    call write_defaults(file, grotop)

    ! atom type
    call write_atom_types(file, grotop)

    ! bond type
    call write_bond_types(file, grotop)

    ! angle type
    call write_angle_types(file, grotop)

    ! dihedral type
    call write_dihedral_types(file, grotop)

    ! constraint type
    call write_constraint_types(file, grotop)

    ! cmap type
    call write_cmap_types(file, grotop)

    ! nonbond parameter
    call write_nonbond_parms(file, grotop)

    ! molecule type
    call write_molecule_type(file, grotop)

    ! system
    call write_system(file, grotop)

    ! molecules
    call write_molecules(file, grotop)

    return

  end subroutine write_grotop

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_defaults
  !> @brief        write section [defaults]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_defaults(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(in)    :: grotop

    ! local variables
    character(16)            :: yesno


    if (grotop%defaults%gen_pairs) then
      yesno = 'yes'
    else
      yesno = 'no'
    end if

    write(file, '(a)') &
         '[ defaults ]'
    write(file, '(a)') &
         '; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ'

    write(file, '(i10,1x,i10,1x,a4,1x,e16.9,1x,e16.9)') &
         grotop%defaults%nonb_func,  &
         grotop%defaults%combi_rule, &
         trim(yesno),                &
         grotop%defaults%fudge_lj,   &
         grotop%defaults%fudge_qq

    write(file, '(a)') ''

    return

  end subroutine write_defaults

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_atom_types
  !> @brief        write section [atomtypes]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_atom_types(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(in)    :: grotop

    ! local variables
    integer                  :: i, cnt
    logical                  :: alias, atnr, khhst
    character(60)            :: vw_str


    cnt = size_grotop(grotop, GroTopAtomType)
    if (cnt == 0) &
      return

    alias = .false.
    atnr  = .false.
    khhst = .false.

    if (grotop%atomtypes(1)%alias /= '') &
      alias = .true.
    if (grotop%atomtypes(1)%type_no /= 0) &
      atnr  = .true.

    if (grotop%atomtypes(1)%khh > 0.0_wp .or. &
        grotop%atomtypes(1)%stokes_r > 0.0_wp) &
      khhst = .true.

    if (grotop%defaults%combi_rule == 1) then
      vw_str = '  c6  c12'
    else
      vw_str = '  sigma  epsilon'
    end if

    write(file, '(a)') '[ atomtypes ]'

    if (alias .and. atnr) then

      write(file, '(a)') &
           ';alias  name  at.num  mass  charge  ptype'//trim(vw_str)

      do i = 1, cnt
        write(file, '(a8,1x,a8,1x,i8,1x,e16.9,1x,e16.9,1x,a8,1x,e16.9,1x,e16.9)')&
             grotop%atomtypes(i)%alias,     &
             grotop%atomtypes(i)%type_name, &
             grotop%atomtypes(i)%type_no,   &
             grotop%atomtypes(i)%mass,      &
             grotop%atomtypes(i)%charge,    &
             grotop%atomtypes(i)%ptype,     &
             grotop%atomtypes(i)%v,         &
             grotop%atomtypes(i)%w

      end do

    else if (alias .and. .not. atnr) then

      write(file, '(a)') &
           ';alias  name  mass  charge  ptype'//trim(vw_str)

      do i = 1, cnt
        write(file, '(a8,1x,a8,1x,e16.9,1x,e16.9,1x,a8,1x,e16.9,1x,e16.9)') &
             grotop%atomtypes(i)%alias,     &
             grotop%atomtypes(i)%type_name, &
             grotop%atomtypes(i)%mass,      &
             grotop%atomtypes(i)%charge,    &
             grotop%atomtypes(i)%ptype,     &
             grotop%atomtypes(i)%v,         &
             grotop%atomtypes(i)%w

      end do

    else if (.not. alias .and. atnr) then

      write(file, '(a)') &
           ';name  at.num  mass  charge  ptype'//trim(vw_str)

      do i = 1, cnt
        write(file, '(a8,1x,i8,1x,e16.9,1x,e16.9,1x,a8,1x,e16.9,1x,e16.9)') &
             grotop%atomtypes(i)%type_name, &
             grotop%atomtypes(i)%type_no,   &
             grotop%atomtypes(i)%mass,      &
             grotop%atomtypes(i)%charge,    &
             grotop%atomtypes(i)%ptype,     &
             grotop%atomtypes(i)%v,         &
             grotop%atomtypes(i)%w

      end do

    else if (khhst) then

      write(file, '(a)') &
           ';name  mass  charge  ptype'//trim(vw_str)//' khh  stokes_r'

      do i = 1, cnt
        write(file, &
             '(a8,1x,e16.9,1x,e16.9,1x,a8,1x,e16.9,1x,e16.9,1x,2e16.9)') &
             grotop%atomtypes(i)%type_name, &
             grotop%atomtypes(i)%mass,      &
             grotop%atomtypes(i)%charge,    &
             grotop%atomtypes(i)%ptype,     &
             grotop%atomtypes(i)%v,         &
             grotop%atomtypes(i)%w,         &
             grotop%atomtypes(i)%khh,       &
             grotop%atomtypes(i)%stokes_r

      end do

    else ! .not. alias .and. .not. atnr .and. .not. khhst

      write(file, '(a)') &
           ';name  mass  charge  ptype'//trim(vw_str)

      do i = 1, cnt
        write(file, &
             '(a8,1x,e16.9,1x,e16.9,1x,a8,1x,e16.9,1x,e16.9)') &
             grotop%atomtypes(i)%type_name, &
             grotop%atomtypes(i)%mass,      &
             grotop%atomtypes(i)%charge,    &
             grotop%atomtypes(i)%ptype,     &
             grotop%atomtypes(i)%v,         &
             grotop%atomtypes(i)%w

      end do

    end if

    write(file, '(a)') ''

    return

  end subroutine write_atom_types

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_bond_types
  !> @brief        write section [bondtypes]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_bond_types(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(in)    :: grotop

    ! local variables
    integer                  :: i, cnt, func
    character(20)            :: str


    cnt  = size_grotop(grotop, GroTopBondType)
    if (cnt == 0) &
      return

    func = -1

    do i = 1, cnt

      if (func /= grotop%bondtypes(i)%func) then

        func = grotop%bondtypes(i)%func

        write(file, '(a)') '[ bondtypes ]'

        select case (func)

        case (1, 2)
          write(file,'(a)') ';i  j  func  b0  kb'

        case default
          write(str,*) func
          call error_msg('Write_Bond_Types> Unknown function type.'//trim(str))

        end select

      end if

      write(file, '(a8,1x,a8,1x,i8,1x,e16.9,1x,e16.9)') &
           grotop%bondtypes(i)%atom_type1,  &
           grotop%bondtypes(i)%atom_type2,  &
           grotop%bondtypes(i)%func,        &
           grotop%bondtypes(i)%b0,          &
           grotop%bondtypes(i)%kb

    end do

    write(file, '(a)') ''

    return

  end subroutine write_bond_types

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_angle_types
  !> @brief        write section [angletypes]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_angle_types(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(in)    :: grotop

    ! local variables
    integer                  :: i, cnt, func
    character(20)            :: str


    cnt  = size_grotop(grotop, GroTopAnglType)
    if (cnt == 0) &
      return

    func = -1

    do i = 1, cnt

      if (func /= grotop%angltypes(i)%func) then

        func = grotop%angltypes(i)%func

        write(file, '(a)') '[ angletypes ]'

        select case (func)

        case (1, 2)
          write(file,'(a)') ';i  j  k  func  th0  cth'

        case (5)
          write(file,'(a)') ';i  j  k  func  th0  cth  ub0  cub'

        case default
          write(str,*) func
          call error_msg('Write_Angle_Types> Unknown function type.'//trim(str))

        end select

      end if

      select case (func)

      case (1, 2)

        write(file, '(a8,1x,a8,1x,a8,1x,i8,1x,e16.9,1x,e16.9)') &
             grotop%angltypes(i)%atom_type1,     &
             grotop%angltypes(i)%atom_type2,     &
             grotop%angltypes(i)%atom_type3,     &
             grotop%angltypes(i)%func,           &
             grotop%angltypes(i)%theta_0,        &
             grotop%angltypes(i)%kt

      case (5)

       write(file,'(a8,1x,a8,1x,a8,1x,i8,1x,e16.9,1x,e16.9,1x,e16.9,1x,e16.9)')&
             grotop%angltypes(i)%atom_type1,     &
             grotop%angltypes(i)%atom_type2,     &
             grotop%angltypes(i)%atom_type3,     &
             grotop%angltypes(i)%func,           &
             grotop%angltypes(i)%theta_0,        &
             grotop%angltypes(i)%kt,             &
             grotop%angltypes(i)%r13,            &
             grotop%angltypes(i)%kub

      end select

    end do

    write(file, '(a)') ''

    return

  end subroutine write_angle_types

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_dihedral_types
  !> @brief        write section [dihedraltypes]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_dihedral_types(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(in)    :: grotop

    ! local variables
    integer                  :: i, cnt, func
    character(20)            :: str


    cnt  = size_grotop(grotop, GroTopDiheType)
    if (cnt == 0) &
      return

    func = -1

    do i = 1, cnt

      if (func /= grotop%dihetypes(i)%func) then

        func = grotop%dihetypes(i)%func

        write(file, '(a)') '[ dihedraltypes ]'

        if (grotop%dihetypes(i)%four_type) then

          select case (func)

          case (1, 4, 9)
            ! dihedral
            write(file,'(a)') ';i  j  k  l  func  ps  kp  mult'
          case (2)
            ! improper dihedral
            write(file,'(a)') ';i  j  k  l  func  ps  kp'
          case (3)
            ! Ryckaert-Bellemans dihedral
            write(file,'(a)') ';i  j  k  l  func  c1  c2 .. c6'
          case default
            write(str,*) func
            call error_msg('Write_Dihedral_Types> Unknown function type.'//&
                 trim(str))

          end select

        else

          select case (func)

          case (1, 4, 9)
            ! dihedral
            write(file,'(a)') ';j  k  func  ps  kp  mult'
          case (2)
            ! improper dihedral
            write(file,'(a)') ';i  l  func  ps  kp'
          case (3)
            ! Ryckaert-Bellemans dihedral
            write(file,'(a)') ';j  k  func  c1  c2 .. c6'
          case default
            write(str,*) func
            call error_msg('Write_Dihedral_Types> Unknown function type.'//&
                 trim(str))

          end select

        end if

      end if

      if (grotop%dihetypes(i)%four_type) then

        select case (func)

        case (1, 4, 9)
          ! dihedral
          write(file,'(a8,1x,a8,1x,a8,1x,a8,1x,i8,1x,e16.9,1x,e16.9,1x,i8)') &
               grotop%dihetypes(i)%atom_type1, &
               grotop%dihetypes(i)%atom_type2, &
               grotop%dihetypes(i)%atom_type3, &
               grotop%dihetypes(i)%atom_type4, &
               grotop%dihetypes(i)%func,       &
               grotop%dihetypes(i)%ps,         &
               grotop%dihetypes(i)%kp,         &
               grotop%dihetypes(i)%multiplicity

        case (2)
          ! improper dihedral
          write(file,'(a8,1x,a8,1x,a8,1x,a8,1x,i8,1x,e16.9,1x,e16.9)') &
               grotop%dihetypes(i)%atom_type1, &
               grotop%dihetypes(i)%atom_type2, &
               grotop%dihetypes(i)%atom_type3, &
               grotop%dihetypes(i)%atom_type4, &
               grotop%dihetypes(i)%func,       &
               grotop%dihetypes(i)%ps,         &
               grotop%dihetypes(i)%kp

        case (3)
          ! Ryckaert-Bellemans dihedral
          write(file,'(a8,1x,a8,1x,a8,1x,a8,1x,i8,1x,6(e16.9,1x))') &
               grotop%dihetypes(i)%atom_type1, &
               grotop%dihetypes(i)%atom_type2, &
               grotop%dihetypes(i)%atom_type3, &
               grotop%dihetypes(i)%atom_type4, &
               grotop%dihetypes(i)%func,       &
               grotop%dihetypes(i)%c(1:6)

        end select

      else

        select case (func)

        case (1, 4, 9)
          ! dihedral
          write(file,'(a8,1x,a8,1x,i8,1x,e16.9,1x,e16.9,1x,i8)') &
               grotop%dihetypes(i)%atom_type2, &
               grotop%dihetypes(i)%atom_type3, &
               grotop%dihetypes(i)%func,       &
               grotop%dihetypes(i)%ps,         &
               grotop%dihetypes(i)%kp,         &
               grotop%dihetypes(i)%multiplicity

        case (2)
          ! improper dihedral
          write(file,'(a8,1x,a8,1x,i8,1x,e16.9,1x,e16.9)') &
               grotop%dihetypes(i)%atom_type1, &
               grotop%dihetypes(i)%atom_type4, &
               grotop%dihetypes(i)%func,       &
               grotop%dihetypes(i)%ps,         &
               grotop%dihetypes(i)%kp

        case (3)
          ! Ryckaert-Bellemans dihedral
          write(file,'(a8,1x,a8,1x,i8,1x,6(e16.9,1x))') &
               grotop%dihetypes(i)%atom_type2, &
               grotop%dihetypes(i)%atom_type3, &
               grotop%dihetypes(i)%func,       &
               grotop%dihetypes(i)%c(1:6)

        end select

      end if

    end do

    write(file, '(a)') ''

    return

  end subroutine write_dihedral_types

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_constraint_types
  !> @brief        write section [constrainttypes]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_constraint_types(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(in)    :: grotop

    ! local variables
    integer                  :: i, cnt, func
    character(20)            :: str


    cnt  = size_grotop(grotop, GroTopConstrType)
    if (cnt == 0) &
      return

    func = -1

    do i = 1, cnt

      if (func /= grotop%constrtypes(i)%func) then

        func = grotop%constrtypes(i)%func

        write(file, '(a)') '[ constrainttypes ]'

        select case (func)

        case (1, 2)
          write(file,'(a)') ';i  j  func  b0'

        case default
          write(str,*) func
          call error_msg('Write_Constraint_Types> Unknown function type.'//&
               trim(str))

        end select

      end if

      write(file, '(a8,1x,a8,1x,i8,1x,e16.9)') &
           grotop%constrtypes(i)%atom_type1, &
           grotop%constrtypes(i)%atom_type2, &
           grotop%constrtypes(i)%func,       &
           grotop%constrtypes(i)%b0

    end do

    write(file, '(a)') ''

    return

  end subroutine write_constraint_types

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_cmap_types
  !> @brief        write section [cmaptypes]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_cmap_types(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(in)    :: grotop

    ! local variables
    integer                  :: i, j, k, cnt, func, nrow, ncol
    character(20)            :: str, nstr
    character(2)             :: term


    cnt  = size_grotop(grotop, GroTopCmapType)
    if (cnt == 0) &
      return

    func = -1

    do i = 1, cnt

      if (func /= grotop%cmaptypes(i)%func) then

        func = grotop%cmaptypes(i)%func

        write(file, '(a)') '[ cmaptypes ]'

        select case (func)

        case (1)
          write(file,'(a)') ';i  j  k  l  m  func  nrow  ncol  map'

        case default
          write(str,*) func
          call error_msg('Write_Cmap_Types> Unknown function type.'//trim(str))

        end select

      end if

      nrow = size(grotop%cmaptypes(i)%map(1,:))
      ncol = size(grotop%cmaptypes(i)%map(:,1))

      write(file, '(a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,i8,1x,i8,1x,i8,1x,a)') &
           grotop%cmaptypes(i)%atom_type1,     &
           grotop%cmaptypes(i)%atom_type2,     &
           grotop%cmaptypes(i)%atom_type3,     &
           grotop%cmaptypes(i)%atom_type4,     &
           grotop%cmaptypes(i)%atom_type5,     &
           grotop%cmaptypes(i)%func,           &
           nrow, ncol, char(92)
      write(nstr,'(i0)') ncol
      do j = 1, nrow
        if (j < nrow) then
          term = char(92)
        else
          term = ''
        end if
        write(file,'('//trim(nstr)//'(e16.9,1x),a)') &
             (grotop%cmaptypes(i)%map(k,j),k=1,ncol),term
      end do

    end do

    write(file, '(a)') ''

    return

  end subroutine write_cmap_types

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_nonbond_parms
  !> @brief        write section [nonbond_params]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_nonbond_parms(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(in)    :: grotop

    ! local variables
    integer                  :: i, cnt, func
    character(20)            :: str


    cnt  = size_grotop(grotop, GroTopNbonParm)
    if (cnt == 0) &
      return

    func = -1

    do i = 1, cnt

      if (func /= grotop%nbonparms(i)%func) then

        func = grotop%nbonparms(i)%func

        write(file, '(a)') '[ nonbond_params ]'

        select case (func)

        case (1, 2)
          if (grotop%defaults%combi_rule == 1) then
            write(file,'(a)') ';i  j  c6  c12'
          else
            write(file,'(a)') ';i  j  sigma  epsilon'
          end if

        case default
          write(str,*) func
          call error_msg('Write_Nonbond_Parms> Unknown function type.'//&
               trim(str))

        end select

      end if

      write(file, '(a8,1x,a8,1x,i8,1x,e16.9,1x,e16.9)') &
           grotop%nbonparms(i)%atom_type1,  &
           grotop%nbonparms(i)%atom_type2,  &
           grotop%nbonparms(i)%func,        &
           grotop%nbonparms(i)%v,           &
           grotop%nbonparms(i)%w

    end do

    write(file, '(a)') ''

    return

  end subroutine write_nonbond_parms

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_molecule_type
  !> @brief        write section [moleculetype]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_molecule_type(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(in)    :: grotop

    ! local variables
    integer                  :: i, cnt


    cnt = size_grotop(grotop, GroTopMolType)
    if (cnt == 0) &
      return

    do i = 1, cnt

      write(file,'(a)') '[ moleculetype ]'
      write(file,'(a)') '; molname  nrexcl'

!      write(file,'(a20,1x,i)')      &
!           grotop%moltypes(i)%name, &
!           grotop%moltypes(i)%exclude_nbon
      write(file,'(a20,1x,$)')      &
           grotop%moltypes(i)%name
      write(file,*) grotop%moltypes(i)%exclude_nbon

      write(file,'(a)') ''

      ! atoms
      call write_atoms(file, grotop%moltypes(i)%mol)

      ! bonds
      call write_bonds(file, grotop%moltypes(i)%mol)

      ! angles
      call write_angles(file, grotop%moltypes(i)%mol)

      ! dihedrals
      call write_dihedrals(file, grotop%moltypes(i)%mol)

      ! cmaps
      call write_cmaps(file, grotop%moltypes(i)%mol)

      ! exclusions
      call write_exclusions(file, grotop%moltypes(i)%mol)

      ! constraints
      call write_constraints(file, grotop%moltypes(i)%mol)

      ! pairs
      call write_pairs(file, grotop%moltypes(i)%mol)

      ! settles
      call write_settles(file, grotop%moltypes(i)%mol)

      ! virtual sites 2
      call write_virtual_sites2(file, grotop%moltypes(i)%mol)

      ! virtual sites 3
      call write_virtual_sites3(file, grotop%moltypes(i)%mol)

      ! virtual sites 4
      call write_virtual_sites4(file, grotop%moltypes(i)%mol)

      ! virtual sites n
      call write_virtual_sitesn(file, grotop%moltypes(i)%mol)

      ! position restraints
      call write_position_restraints(file, grotop%moltypes(i)%mol)

    end do

    return

  end subroutine write_molecule_type

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_atoms
  !> @brief        write section [atoms]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_atoms(file, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop_mol),      intent(in)    :: gromol

    ! local variables
    integer                  :: i, cnt
    logical                  :: mass, stok


    cnt = size_grotop_mol(gromol, GMGroMolAtom)
    if (cnt == 0) &
      return

    mass = .false.
    stok = .false.

    if (gromol%atoms(1)%mass /= 0.0_wp) &
      mass = .true.

    if (gromol%atoms(1)%stokes_r > 0.0_wp) &
      stok = .true.

    write(file,'(a)') '[ atoms ]'

    if (mass) then
      if (stok) then
        write(file,'(a)') &
        '; id  at-type  res-nr  res-name at-name  cg-nr  charge  mass stokes_r'
      else
        write(file,'(a)') &
        '; id  at-type  res-nr  res-name at-name  cg-nr  charge  mass'
      end if
    else
      write(file,'(a)') &
        '; id  at-type  res-nr  res-name at-name  cg-nr  charge'
    end if

    do i = 1, cnt

      if (mass) then

        if (stok) then

          write(file, &
            '(i8,1x,a6,1x,i8,1x,a6,1x,a4,1x,i8,1x,e16.9,1x,e16.9,1x,e16.9)') &
             gromol%atoms(i)%atom_no,      &
             gromol%atoms(i)%atom_type,    &
             gromol%atoms(i)%residue_no,   &
             gromol%atoms(i)%residue_name, &
             gromol%atoms(i)%atom_name,    &
             gromol%atoms(i)%charge_group, &
             gromol%atoms(i)%charge,       &
             gromol%atoms(i)%mass,         &
             gromol%atoms(i)%stokes_r

        else

          write(file,'(i8,1x,a6,1x,i8,1x,a6,1x,a4,1x,i8,1x,e16.9,1x,e16.9)') &
             gromol%atoms(i)%atom_no,      &
             gromol%atoms(i)%atom_type,    &
             gromol%atoms(i)%residue_no,   &
             gromol%atoms(i)%residue_name, &
             gromol%atoms(i)%atom_name,    &
             gromol%atoms(i)%charge_group, &
             gromol%atoms(i)%charge,       &
             gromol%atoms(i)%mass

        end if

      else

        write(file,'(i8,1x,a6,1x,i8,1x,a6,1x,a4,1x,i8,1x,e16.9)') &
             gromol%atoms(i)%atom_no,      &
             gromol%atoms(i)%atom_type,    &
             gromol%atoms(i)%residue_no,   &
             gromol%atoms(i)%residue_name, &
             gromol%atoms(i)%atom_name,    &
             gromol%atoms(i)%charge_group, &
             gromol%atoms(i)%charge

      end if

    end do

    write(file,'(a)') ''

    return

  end subroutine write_atoms

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_bonds
  !> @brief        write section [bonds]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_bonds(file, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop_mol),      intent(in)    :: gromol

    ! local variables
    integer                  :: i, cnt


    cnt = size_grotop_mol(gromol, GMGroMolBond)
    if (cnt == 0) &
      return

    write(file,'(a)') '[ bonds ]'
    write(file,'(a)') '; i  j  funct  length  force.c.'

    do i = 1, cnt

      write(file,'(i8,1x,i8,1x,i8,1x,e16.9,1x,e16.9)') &
           gromol%bonds(i)%atom_idx1, &
           gromol%bonds(i)%atom_idx2, &
           gromol%bonds(i)%func,      &
           gromol%bonds(i)%b0,        &
           gromol%bonds(i)%kb

    end do

    write(file,'(a)') ''

    return

  end subroutine write_bonds

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_angles
  !> @brief        write section [angles]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_angles(file, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop_mol),      intent(in)    :: gromol

    ! local variables
    integer                  :: i, cnt


    cnt = size_grotop_mol(gromol, GMGroMolAngl)
    if (cnt == 0) &
      return

    write(file,'(a)') '[ angles ]'

    select case(gromol%angls(i)%func)

    case (1,2)
      write(file,'(a)') '; i  j  k  funct  angle  force.c.'

    case (5)
      write(file,'(a)') '; i  j  k  funct  angle  force.c.  r13  kub'

    end select

    do i = 1, cnt

      select case(gromol%angls(i)%func)

      case (1,2)

        write(file,'(i8,1x,i8,1x,i8,1x,i8,1x,e16.9,1x,e16.9)') &
             gromol%angls(i)%atom_idx1, &
             gromol%angls(i)%atom_idx2, &
             gromol%angls(i)%atom_idx3, &
             gromol%angls(i)%func,      &
             gromol%angls(i)%theta_0,   &
             gromol%angls(i)%kt

      case (5)

       write(file,'(i8,1x,i8,1x,i8,1x,i8,1x,e16.9,1x,e16.9,1x,e16.9,1x,e16.9)')&
             gromol%angls(i)%atom_idx1, &
             gromol%angls(i)%atom_idx2, &
             gromol%angls(i)%atom_idx3, &
             gromol%angls(i)%func,      &
             gromol%angls(i)%theta_0,   &
             gromol%angls(i)%kt,        &
             gromol%angls(i)%r13,       &
             gromol%angls(i)%kub

      end select

    end do

    write(file,'(a)') ''

    return

  end subroutine write_angles

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_dihedrals
  !> @brief        write section [dihedrals]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_dihedrals(file, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop_mol),      intent(in)    :: gromol

    ! local variables
    integer                  :: i, cnt
    character(20)            :: str


    cnt = size_grotop_mol(gromol, GMGroMolDihe)
    if (cnt == 0) &
      return

    write(file,'(a)') '[ dihedrals ]'

    select case (gromol%dihes(1)%func)

    case (1,4,9)
      ! dihedral
      write(file,'(a)') '; ai  aj  ak  al  funct  ph0  cp  mult'

    case (2)
      ! improper dihedral
      write(file,'(a)') '; ai  aj  ak  al  funct  ph0  cp'

    case (3)
      ! Ryckaert-Bellemans dihedral
      write(file,'(a)') '; ai  aj  ak  al  funct  c1 ... c6'

    case default
      write(str,*) gromol%dihes(1)%func
      call error_msg('Write_Dihedrals> Unknown function type.'//trim(str))

    end select

    do i = 1, cnt

      select case (gromol%dihes(i)%func)

      case (1,4,9)

        write(file,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,e16.9,1x,e16.9,1x,i8)') &
             gromol%dihes(i)%atom_idx1, &
             gromol%dihes(i)%atom_idx2, &
             gromol%dihes(i)%atom_idx3, &
             gromol%dihes(i)%atom_idx4, &
             gromol%dihes(i)%func,      &
             gromol%dihes(i)%ps,        &
             gromol%dihes(i)%kp,        &
             gromol%dihes(i)%multiplicity

      case (2)

        write(file,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,e16.9,1x,e16.9)') &
             gromol%dihes(i)%atom_idx1, &
             gromol%dihes(i)%atom_idx2, &
             gromol%dihes(i)%atom_idx3, &
             gromol%dihes(i)%atom_idx4, &
             gromol%dihes(i)%func,      &
             gromol%dihes(i)%ps,        &
             gromol%dihes(i)%kp

      case (3)

        write(file,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,6(e16.9,1x))') &
             gromol%dihes(i)%atom_idx1, &
             gromol%dihes(i)%atom_idx2, &
             gromol%dihes(i)%atom_idx3, &
             gromol%dihes(i)%atom_idx4, &
             gromol%dihes(i)%func,      &
             gromol%dihes(i)%c(1:6)

      end select

    end do

    write(file,'(a)') ''

    return

  end subroutine write_dihedrals

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_cmaps
  !> @brief        write section [cmap]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_cmaps(file, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop_mol),      intent(in)    :: gromol

    ! local variables
    integer                  :: i, j, k, cnt, nrow, ncol
    character(20)            :: nstr
    character(2)             :: term


    cnt = size_grotop_mol(gromol, GMGroMolCmap)
    if (cnt == 0) &
      return

    write(file,'(a)') '[ cmap ]'
    write(file,'(a)') '; i  j  k  l  m  funct  nrow  ncol  map'

    do i = 1, cnt

      nrow = size(gromol%cmaps(i)%map(1,:))
      ncol = size(gromol%cmaps(i)%map(:,1))

      write(file, '(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,a)') &
           gromol%cmaps(i)%atom_idx1,     &
           gromol%cmaps(i)%atom_idx2,     &
           gromol%cmaps(i)%atom_idx3,     &
           gromol%cmaps(i)%atom_idx4,     &
           gromol%cmaps(i)%atom_idx5,     &
           gromol%cmaps(i)%func,          &
           nrow, ncol, char(92)
      write(nstr,'(i0)') ncol
      do j = 1, nrow
        if (j < nrow) then
          term = char(92)
        else
          term = ''
        end if
        write(file,'('//trim(nstr)//'(e16.9,1x),a)') &
             (gromol%cmaps(i)%map(k,j),k=1,ncol),term
      end do

    end do

    write(file,'(a)') ''

    return

  end subroutine write_cmaps

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_exclusions
  !> @brief        write section [exclusions]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_exclusions(file, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop_mol),      intent(in)    :: gromol

    ! local variables
    integer                  :: i, cnt, old1
    character(1000)          :: line, linep


    cnt = size_grotop_mol(gromol, GMGroMolExcl)

    if (cnt == 0) &
      return

    write(file,'(a)') '[ exclusions ]'

    old1 = gromol%excls(1)%atom_idx1
    write(line,'(i10)') old1

    do i = 1, cnt

      if (old1 /= gromol%excls(i)%atom_idx1) then
        write(file,'(a)') trim(line)
        old1 = gromol%excls(i)%atom_idx1
        write(line,'(i10)') old1
      end if

      write(linep,'(1x,i10)') gromol%excls(i)%atom_idx2
      line = trim(line) // linep

    end do

    write(file,'(a)') trim(line)
    write(file,'(a)') ''

    return

  end subroutine write_exclusions

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_constraints
  !> @brief        write section [constraints]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_constraints(file, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop_mol),      intent(in)    :: gromol

    ! local variables
    integer                  :: i, cnt


    cnt = size_grotop_mol(gromol, GMGroMolConstr)
    if (cnt == 0) &
      return

    write(file,'(a)') '[ constraints ]'
    write(file,'(a)') '; i  j  funct  length'

    do i = 1, cnt

      write(file,'(i8,1x,i8,1x,i8,1x,e16.9)') &
           gromol%constrs(i)%atom_idx1, &
           gromol%constrs(i)%atom_idx2, &
           gromol%constrs(i)%func,      &
           gromol%constrs(i)%b0

    end do

    write(file,'(a)') ''

    return

  end subroutine write_constraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_pairs
  !> @brief        write section [pairs]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_pairs(file, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop_mol),      intent(in)    :: gromol

    ! local variables
    integer                  :: i, cnt
    logical                  :: khh


    cnt = size_grotop_mol(gromol, GMGroMolPair)
    if (cnt == 0) &
      return

    khh = .false.
    if (gromol%pairs(1)%khh /= 0.0_wp) &
      khh = .true.

    write(file,'(a)') '[ pairs ]'

    if (khh) then
      write(file,'(a)') '; i  j  func  v  w  khh'
    else
      write(file,'(a)') '; i  j  func  v  w'
    end if

    do i = 1, cnt

      if (khh) then

        write(file,'(i8,1x,i8,1x,i8,1x,e16.9,1x,e16.9,1x,e16.9)') &
             gromol%pairs(i)%atom_idx1, &
             gromol%pairs(i)%atom_idx2, &
             gromol%pairs(i)%func,      &
             gromol%pairs(i)%v,         &
             gromol%pairs(i)%w,         &
             gromol%pairs(i)%khh

      else

        write(file,'(i8,1x,i8,1x,i8,1x,e16.9,1x,e16.9)') &
             gromol%pairs(i)%atom_idx1, &
             gromol%pairs(i)%atom_idx2, &
             gromol%pairs(i)%func,      &
             gromol%pairs(i)%v,         &
             gromol%pairs(i)%w

      end if

    end do

    write(file,'(a)') ''

    return

  end subroutine write_pairs

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_settles
  !> @brief        write section [settles]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_settles(file, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop_mol),      intent(in)    :: gromol


    if (gromol%settles%func == 0) &
      return

    write(file,'(a)') '[ settles ]'
    write(file,'(a)') '; OW  func  dOH  dHH'

    write(file,'(i8,1x,i8,1x,e16.9,1x,e16.9)') &
         gromol%settles%ow,   &
         gromol%settles%func, &
         gromol%settles%doh,  &
         gromol%settles%dhh

    write(file,'(a)') ''

    return

  end subroutine write_settles

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_virtual_sites2
  !> @brief        write section [virtual_sites2]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_virtual_sites2(file, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop_mol),      intent(in)    :: gromol

    ! local variables
    integer                  :: i, cnt
    character(20)            :: str


    cnt = size_grotop_mol(gromol, GMGroMolVsites2)
    if (cnt == 0) &
      return

    write(file,'(a)') '[ virtual_sites2 ]'

    select case (gromol%vsites2(1)%func)

    case (1)
      ! 2-body virtual site
      write(file,'(a)') '; Site  from  funct  a'

    case default
      write(str,*) gromol%vsites2(1)%func
      call error_msg('Write_Virtual_Sites2> Unknown function type.'//trim(str))

    end select

    do i = 1, cnt

      select case (gromol%vsites2(i)%func)

      case (1)
        write(file,'(i8,1x,i8,1x,i8,1x,i8,1x,e16.9)') &
             gromol%vsites2(i)%atom_idx1, &
             gromol%vsites2(i)%atom_idx2, &
             gromol%vsites2(i)%atom_idx3, &
             gromol%vsites2(i)%func,      &
             gromol%vsites2(i)%a

      end select

    end do

    write(file,'(a)') ''

    return

  end subroutine write_virtual_sites2

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_virtual_sites3
  !> @brief        write section [virtual_sites3]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_virtual_sites3(file, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop_mol),      intent(in)    :: gromol

    ! local variables
    integer                  :: i, cnt
    character(20)            :: str


    cnt = size_grotop_mol(gromol, GMGroMolVsites3)
    if (cnt == 0) &
      return

    write(file,'(a)') '[ virtual_sites3 ]'

    select case (gromol%vsites3(1)%func)

    case (1)
      ! 3-body virtual site
      write(file,'(a)') '; Site  from  funct  a  b'

    case (2)
      ! 3-body virtual site (fd)
      write(file,'(a)') '; Site  from  funct  a  d'

    case (3)
      ! 3-body virtual site (fad)
      write(file,'(a)') '; Site  from  funct  theta  d'

    case (4)
      ! 3-body virtual site (fdn)
      write(file,'(a)') '; Site  from  funct  a  b  c'

    case default
      write(str,*) gromol%vsites3(1)%func
      call error_msg('Write_Virtual_Sites3> Unknown function type.'//trim(str))

    end select

    do i = 1, cnt

      select case (gromol%vsites3(i)%func)

      case (1)
        write(file,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,e16.9,1x,e16.9)') &
             gromol%vsites3(i)%atom_idx1, &
             gromol%vsites3(i)%atom_idx2, &
             gromol%vsites3(i)%atom_idx3, &
             gromol%vsites3(i)%atom_idx4, &
             gromol%vsites3(i)%func,      &
             gromol%vsites3(i)%a,         &
             gromol%vsites3(i)%b

      case (2)
        write(file,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,e16.9,1x,e16.9)') &
             gromol%vsites3(i)%atom_idx1, &
             gromol%vsites3(i)%atom_idx2, &
             gromol%vsites3(i)%atom_idx3, &
             gromol%vsites3(i)%atom_idx4, &
             gromol%vsites3(i)%func,      &
             gromol%vsites3(i)%a,         &
             gromol%vsites3(i)%d

      case (3)
        write(file,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,e16.9,1x,e16.9)') &
             gromol%vsites3(i)%atom_idx1, &
             gromol%vsites3(i)%atom_idx2, &
             gromol%vsites3(i)%atom_idx3, &
             gromol%vsites3(i)%atom_idx4, &
             gromol%vsites3(i)%func,      &
             gromol%vsites3(i)%theta,     &
             gromol%vsites3(i)%d

      case (4)
        write(file,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,e16.9,1x,e16.9,1x,e16.9)') &
             gromol%vsites3(i)%atom_idx1, &
             gromol%vsites3(i)%atom_idx2, &
             gromol%vsites3(i)%atom_idx3, &
             gromol%vsites3(i)%atom_idx4, &
             gromol%vsites3(i)%func,      &
             gromol%vsites3(i)%a,         &
             gromol%vsites3(i)%b,         &
             gromol%vsites3(i)%c

      end select

    end do

    write(file,'(a)') ''

    return

  end subroutine write_virtual_sites3

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_virtual_sites4
  !> @brief        write section [virtual_sites4]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_virtual_sites4(file, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop_mol),      intent(in)    :: gromol

    ! local variables
    integer                  :: i, cnt
    character(20)            :: str


    cnt = size_grotop_mol(gromol, GMGroMolVsites4)
    if (cnt == 0) &
      return

    write(file,'(a)') '[ virtual_sites4 ]'

    select case (gromol%vsites4(1)%func)

    case (2)
      ! 4-body virtual site (fdn)
      write(file,'(a)') '; Site  from  funct  a  b  c'

    case default
      write(str,*) gromol%vsites4(1)%func
      call error_msg('Write_Virtual_Sites4> Unknown function type.'//trim(str))

    end select

    do i = 1, cnt

      select case (gromol%vsites4(i)%func)

      case (2)
        write(file, &
        '(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,e16.9,1x,e16.9,1x,e16.9)') &
             gromol%vsites4(i)%atom_idx1, &
             gromol%vsites4(i)%atom_idx2, &
             gromol%vsites4(i)%atom_idx3, &
             gromol%vsites4(i)%atom_idx4, &
             gromol%vsites4(i)%atom_idx5, &
             gromol%vsites4(i)%func,      &
             gromol%vsites4(i)%a,         &
             gromol%vsites4(i)%b,         &
             gromol%vsites4(i)%c

      end select

    end do

    write(file,'(a)') ''

    return

  end subroutine write_virtual_sites4

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_virtual_sitesn
  !> @brief        write section [virtual_sitesn]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_virtual_sitesn(file, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop_mol),      intent(in)    :: gromol

    ! local variables
    integer                  :: i, cnt, nidx
    character(20)            :: str, nstr


    cnt = size_grotop_mol(gromol, GMGroMolVsitesn)
    if (cnt == 0) &
      return

    write(file,'(a)') '[ virtual_sitesn ]'

    select case (gromol%vsitesn(1)%func)

    case (1,2,3)
      ! N-body virtual site (COG)
      ! N-body virtual site (COM)
      ! N-body virtual site (COW)
      write(file,'(a)') '; Site  funct  from'

    case default
      write(str,*) gromol%vsitesn(1)%func
      call error_msg('Write_Virtual_SitesN> Unknown function type.'//trim(str))

    end select

    do i = 1, cnt

      select case (gromol%vsitesn(i)%func)

      case (1,2,3)

        nidx = size(gromol%vsitesn(i)%atom_idcs)
        write(nstr,'(i0)') nidx

        write(file, '(i8,1x,i8,1x,'//trim(nstr)//'(i8,1x))') &
             gromol%vsitesn(i)%atom_idx, &
             gromol%vsitesn(i)%func,     &
             gromol%vsitesn(i)%atom_idcs(1:nidx)

      end select

    end do

    write(file,'(a)') ''

    return

  end subroutine write_virtual_sitesn

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_position_restraints
  !> @brief        write section [position_restraints]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_position_restraints(file, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop_mol),      intent(in)    :: gromol

    ! local variables
    integer                  :: i, cnt


    cnt = size_grotop_mol(gromol, GMGroMolPosres)
    if (cnt == 0) &
      return

    write(file,'(a)') '[ position_restraints ]'
    write(file,'(a)') '; iatom type  fx  fy  fz'

    do i = 1, cnt

      write(file,'(i8,1x,i8,1x,e16.9,1x,e16.9,1x,e16.9)') &
           gromol%posress(i)%atom_idx, &
           gromol%posress(i)%func,     &
           gromol%posress(i)%kx,       &
           gromol%posress(i)%ky,       &
           gromol%posress(i)%kz

    end do

    write(file,'(a)') ''

    return

  end subroutine write_position_restraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_system
  !> @brief        write section [system]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_system(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(in)    :: grotop


    return

  end subroutine write_system

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_molecules
  !> @brief        write section [molecules]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_molecules(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(in)    :: grotop

    ! local variables
    integer                  :: i, cnt


    cnt = size_grotop(grotop, GroTopMols)
    if (cnt == 0) &
      return

    write(file,'(a)') '[ molecules ]'

    do i = 1, cnt

      write(file,'(a20,1x,i12)') &
           grotop%molss(i)%name, &
           grotop%molss(i)%count

    end do

    return

  end subroutine write_molecules

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      check_section
  !> @brief        check section line count from current file pointer
  !! @authors      NT
  !! @param[in]    file : file unit number
  !! @return       number of lines
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function check_section_count(file)

    ! return value
    integer :: check_section_count

    ! formal arguments
    integer,                 intent(in)    :: file

    ! local variables
    integer                  :: i, cnt, directive
    character(MaxLine)       :: line


    cnt = 0

    do while(.true.)
      read(file,'(A)',end=100) line
      if (check_directive(line, directive)) &
        exit
      cnt = cnt + 1
    end do

100 do i = 1, cnt + 1
      backspace(file)
    end do

    check_section_count = cnt
    return

  end function check_section_count

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      check_directive
  !> @brief        check the directive
  !! @authors      NT
  !! @param[in]    line      : read line
  !! @param[out]   directive : directive ID
  !! @return       flag for directive was found or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function check_directive(line, directive)

    ! return value
    logical  :: check_directive

    ! formal arguments
    character(*),            intent(in)    :: line
    integer,                 intent(out)   :: directive

    ! local variables
    integer                  :: idxb, idxe
    character(20)            :: dir_str


    check_directive = .false.

    idxb = index(line,'[')
    if (idxb == 0) &
       return

    idxe = index(line,']')
    if (idxe == 0) &
       idxe = len(line)

    read(line(idxb+1:idxe-1),*,err=100,end=100) dir_str
100 continue

    select case (dir_str)

    case ('defaults')
      directive = DDefaults
    case ('atomtypes')
      directive = DAtomTypes
    case ('bondtypes')
      directive = DBondTypes
    case ('angletypes')
      directive = DAngleTypes
    case ('dihedraltypes')
      directive = DDihedralTypes
    case ('constrainttypes')
      directive = DConstrTypes
    case ('cmaptypes')
      directive = DCmapTypes
    case ('nonbond_params')
      directive = DNonbondParams
    case ('moleculetype')
      directive = DMoleculeType
    case ('atoms')
      directive = DAtoms
    case ('bonds')
      directive = DBonds
    case ('angles')
      directive = DAngles
    case ('dihedrals')
      directive = DDihedrals
    case ('cmap')
      directive = DCmap
    case ('exclusions')
      directive = DExclusions
    case ('constraints')
      directive = DConstraints
    case ('pairs')
      directive = DPairs
    case ('settles')
      directive = DSettles
    case ('virtual_sites2')
      directive = DVirtualSites2
    case ('virtual_sites3')
      directive = DVirtualSites3
    case ('dummies3')
      directive = DVirtualSites3
    case ('virtual_sites4')
      directive = DVirtualSites4
    case ('virtual_sitesn')
      directive = DVirtualSitesn
    case ('position_restraints')
      directive = DPosition_Restraints
    case ('system')
      directive = DSystem
    case ('molecules')
      directive = DMolecules
    case default
      directive = DUnknown
    end select

    check_directive = .true.
    return

  end function check_directive

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      size_grotop
  !> @brief        get the GROTOP allocatable variable size
  !! @authors      NT
  !! @param[in]    grotop   : GROMACS TOP information
  !! @param[in]    variable : parameter for allocatable variable
  !! @return       size of allocatable variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function size_grotop(grotop, variable)

    ! return value
    integer :: size_grotop

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    integer,                 intent(in)    :: variable


    size_grotop = 0

    select case(variable)

    case(GroTopAtomType)
      if (allocated(grotop%atomtypes)) &
        size_grotop = size(grotop%atomtypes)

    case(GroTopBondType)
      if (allocated(grotop%bondtypes)) &
        size_grotop = size(grotop%bondtypes)

    case(GroTopAnglType)
      if (allocated(grotop%angltypes)) &
        size_grotop = size(grotop%angltypes)

    case(GroTopDiheType)
      if (allocated(grotop%dihetypes)) &
        size_grotop = size(grotop%dihetypes)

    case(GroTopConstrType)
      if (allocated(grotop%constrtypes)) &
        size_grotop = size(grotop%constrtypes)

    case(GroTopCmapType)
      if (allocated(grotop%cmaptypes)) &
        size_grotop = size(grotop%cmaptypes)

    case(GroTopNbonParm)
      if (allocated(grotop%nbonparms)) &
        size_grotop = size(grotop%nbonparms)

    case(GroTopMolType)
      if (allocated(grotop%moltypes)) &
        size_grotop = size(grotop%moltypes)

    case(GroTopMols)
      if (allocated(grotop%molss)) &
        size_grotop = size(grotop%molss)

    case default
      call error_msg('Size_Grotop> bad variable')

    end select

    return

  end function size_grotop

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      size_grotop_mol
  !> @brief        get the GROTOP molecule allocatable variable size
  !! @authors      NT
  !! @param[in]    gromol   : GROMACS TOP molecule information
  !! @param[in]    variable : parameter for allocatable variable
  !! @return       size of allocatable variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function size_grotop_mol(gromol, variable)

    ! return value
    integer :: size_grotop_mol

    ! formal arguments
    type(s_grotop_mol),      intent(in)    :: gromol
    integer,                 intent(in)    :: variable


    size_grotop_mol = 0

    select case (variable)

    case (GMGroMolAtom)
      if (allocated(gromol%atoms)) then
        size_grotop_mol = size(gromol%atoms)
      else
        size_grotop_mol =0
      end if

    case (GMGroMolBond)
      if (allocated(gromol%bonds)) then
        size_grotop_mol = size(gromol%bonds)
      else
        size_grotop_mol =0
      end if

    case (GMGroMolAngl)
      if (allocated(gromol%angls)) then
        size_grotop_mol = size(gromol%angls)
      else
        size_grotop_mol =0
      end if

    case (GMGroMolDihe)
      if (allocated(gromol%dihes)) then
        size_grotop_mol = size(gromol%dihes)
      else
        size_grotop_mol =0
      end if

    case (GMGroMolCmap)
      if (allocated(gromol%cmaps)) then
        size_grotop_mol = size(gromol%cmaps)
      else
        size_grotop_mol =0
      end if

    case (GMGroMolExcl)
      if (allocated(gromol%excls)) then
        size_grotop_mol = size(gromol%excls)
      else
        size_grotop_mol =0
      end if

    case (GMGroMolConstr)
      if (allocated(gromol%constrs)) then
        size_grotop_mol = size(gromol%constrs)
      else
        size_grotop_mol =0
      end if

    case (GMGroMolPair)
      if (allocated(gromol%pairs)) then
        size_grotop_mol = size(gromol%pairs)
      else
        size_grotop_mol =0
      end if

    case (GMGroMolVsites2)
      if (allocated(gromol%vsites2)) then
        size_grotop_mol = size(gromol%vsites2)
      else
        size_grotop_mol =0
      end if

    case (GMGroMolVsites3)
      if (allocated(gromol%vsites3)) then
        size_grotop_mol = size(gromol%vsites3)
      else
        size_grotop_mol =0
      end if

    case (GMGroMolVsites4)
      if (allocated(gromol%vsites4)) then
        size_grotop_mol = size(gromol%vsites4)
      else
        size_grotop_mol =0
      end if

    case (GMGroMolVsitesn)
      if (allocated(gromol%vsitesn)) then
        size_grotop_mol = size(gromol%vsitesn)
      else
        size_grotop_mol =0
      end if

    case (GMGroMolPosres)
      if (allocated(gromol%posress)) then
        size_grotop_mol = size(gromol%posress)
      else
        size_grotop_mol =0
      end if

    case default
      call error_msg('Size_Grotop_Mol> bad variable')

    end select

    return

  end function size_grotop_mol

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    error_msg_grotop
  !> @brief        output error message
  !! @authors      NT
  !! @param[in]    file    : file unit number
  !! @param[in]    message : message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine error_msg_grotop(file, message)

    ! formal arguments
    integer,                 intent(in) :: file
    character(*),            intent(in) :: message


    call gro_pp_close_file(file)
    call error_msg(message)

    return

  end subroutine error_msg_grotop

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      search_atom_mass
  !> @brief        search atom mass
  !! @authors      NT
  !! @param[in]    atomtypes : data of atomtypes
  !! @param[in]    at        : atom type name
  !! @param[out]   mass      : mass of atoms
  !! @return       flag for found or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function search_atom_mass(atomtypes, at, mass)

    ! return value
    logical :: search_atom_mass

    ! formal arguments
    type(s_atomtype),        intent(in)    :: atomtypes(:)
    character(*),            intent(in)    :: at
    real(wp),                intent(out)   :: mass

    ! local variables
    integer                  :: i


    do i = 1, size(atomtypes)
      if (atomtypes(i)%type_name == at) then
        mass = atomtypes(i)%mass
        search_atom_mass = .true.
        return
      end if
    end do

    search_atom_mass = .false.
    return

  end function search_atom_mass

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      search_stokes_r
  !> @brief        search stokes radius
  !! @authors      NT
  !! @param[in]    atomtypes : data of atomtypes
  !! @param[in]    at        : atom type name
  !! @param[out]   stokes_r  : stokes radius
  !! @return       flag for found or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function search_stokes_r(atomtypes, at, stokes_r)

    ! return value
    logical :: search_stokes_r

    ! formal arguments
    type(s_atomtype),        intent(in)    :: atomtypes(:)
    character(*),            intent(in)    :: at
    real(wp),                intent(out)   :: stokes_r

    ! local variables
    integer                  :: i


    do i = 1, size(atomtypes)
      if (atomtypes(i)%type_name == at) then
        stokes_r = atomtypes(i)%stokes_r
        search_stokes_r = .true.
        return
      end if
    end do

    search_stokes_r = .false.
    return

  end function search_stokes_r

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      search_bond_param
  !> @brief        search bond parameters
  !! @authors      NT
  !! @param[in]    bondtypes : data of bondtypes
  !! @param[in]    at1       : atom 1 type name
  !! @param[in]    at2       : atom 2 type name
  !! @param[in]    func      : function type
  !! @param[out]   b0        : param b0
  !! @param[out]   kb        : param kb
  !! @return       flag for found or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function search_bond_param(bondtypes, at1, at2, func, b0, kb)

    ! return value
    logical :: search_bond_param

    ! formal arguments
    type(s_bondtype), target, intent(in)    :: bondtypes(:)
    character(*),             intent(in)    :: at1
    character(*),             intent(in)    :: at2
    integer,                  intent(in)    :: func
    real(wp),                 intent(out)   :: b0
    real(wp),                 intent(out)   :: kb

    ! local variables
    integer                   :: i
    type(s_bondtype), pointer :: bt


    do i = 1, size(bondtypes)
      bt => bondtypes(i)
      if (bt%func /= func) &
        cycle
      if (bt%atom_type1 == at1 .and. bt%atom_type2 == at2 .or. &
          bt%atom_type1 == at2 .and. bt%atom_type2 == at1) then
        b0 = bt%b0
        kb = bt%kb
        search_bond_param = .true.
        return
      end if
    end do

    search_bond_param = .false.
    return

  end function search_bond_param

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      search_angl_param
  !> @brief        search angle parameters
  !! @authors      NT
  !! @param[in]    angltypes : data of angletypes
  !! @param[in]    at1       : atom 1 type name
  !! @param[in]    at2       : atom 2 type name
  !! @param[in]    at3       : atom 3 type name
  !! @param[in]    func      : function type
  !! @param[out]   theta_0   : param theta_0
  !! @param[out]   kt        : param kt
  !! @param[out]   r13       : param r13 (optional)
  !! @param[out]   kub       : param kUB (optional)
  !! @return       flag for found or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function search_angl_param(angltypes, at1, at2, at3, func, theta_0,kt,r13,kub)

    ! return value
    logical :: search_angl_param

    ! formal arguments
    type(s_angltype), target, intent(in)    :: angltypes(:)
    character(*),             intent(in)    :: at1
    character(*),             intent(in)    :: at2
    character(*),             intent(in)    :: at3
    integer,                  intent(in)    :: func
    real(wp),                 intent(out)   :: theta_0
    real(wp),                 intent(out)   :: kt
    real(wp),       optional, intent(out)   :: r13
    real(wp),       optional, intent(out)   :: kub

    ! local variables
    integer                   :: i
    type(s_angltype), pointer :: at


    do i = 1, size(angltypes)
      at => angltypes(i)
      if (at%func /= func) &
        cycle
      if (at%atom_type2 /= at2) &
        cycle
      if (at%atom_type1 == at1 .and. at%atom_type3 == at3 .or. &
          at%atom_type1 == at3 .and. at%atom_type3 == at1) then
        theta_0 = at%theta_0
        kt = at%kt
        if (present(r13)) r13 = at%r13
        if (present(kub)) kub = at%kub
        search_angl_param = .true.
        return
      end if
    end do

    search_angl_param = .false.
    return

  end function search_angl_param

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      search_dihe_param
  !> @brief        search dihedral parameters
  !! @authors      NT
  !! @param[in]    dihetypes : data of dihedraltypes
  !! @param[in]    at1       : atom 1 type name
  !! @param[in]    at2       : atom 2 type name
  !! @param[in]    at3       : atom 3 type name
  !! @param[in]    at4       : atom 4 type name
  !! @param[inout] dihe      : dihedral info.
  !! @return       flag for found or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function search_dihe_param(dihetypes, at1, at2, at3, at4, dihe)

    ! return value
    logical :: search_dihe_param

    ! formal arguments
    type(s_dihetype), target, intent(in)    :: dihetypes(:)
    character(*),             intent(in)    :: at1
    character(*),             intent(in)    :: at2
    character(*),             intent(in)    :: at3
    character(*),             intent(in)    :: at4
    type(s_dihe),             intent(inout) :: dihe

    ! local variables
    integer                   :: i
    type(s_dihetype), pointer :: dt


    do i = 1, size(dihetypes)
      dt => dihetypes(i)

      if (dt%func /= dihe%func) &
        cycle

      if (dt%four_type) then

        if (dt%atom_type1 == 'X' .and. dt%atom_type4 == 'X') then

          if ((dt%atom_type2 == at2 .and. &
               dt%atom_type3 == at3) .or. &
              (dt%atom_type2 == at3 .and. &
               dt%atom_type3 == at2)) then
            exit
          end if

        else if (dt%atom_type2 == 'X' .and. dt%atom_type3 == 'X') then

          if ((dt%atom_type1 == at1 .and. &
               dt%atom_type4 == at4) .or. &
              (dt%atom_type1 == at4 .and. &
               dt%atom_type4 == at1)) then
            exit
          end if

        else if (dt%atom_type1 == 'X' .and. dt%atom_type2 == 'X') then

          if ((dt%atom_type3 == at3 .and. &
               dt%atom_type4 == at4) .or. &
              (dt%atom_type3 == at2 .and. &
               dt%atom_type4 == at1)) then
            exit
          end if

        else if (dt%atom_type1 == 'X') then

          if ((dt%atom_type2 == at2 .and. &
               dt%atom_type3 == at3 .and. &
               dt%atom_type4 == at4) .or. &
              (dt%atom_type2 == at3 .and. &
               dt%atom_type3 == at2 .and. &
               dt%atom_type4 == at1)) then
            exit
          end if

        else if (dt%atom_type4 == 'X') then

          if ((dt%atom_type1 == at1 .and. &
               dt%atom_type2 == at2 .and. &
               dt%atom_type3 == at3) .or. &
              (dt%atom_type1 == at4 .and. &
               dt%atom_type2 == at3 .and. &
               dt%atom_type3 == at2)) then
            exit
          end if

        else

          if ((dt%atom_type1 == at1 .and. &
               dt%atom_type2 == at2 .and. &
               dt%atom_type3 == at3 .and. &
               dt%atom_type4 == at4) .or. &
              (dt%atom_type1 == at4 .and. &
               dt%atom_type2 == at3 .and. &
               dt%atom_type3 == at2 .and. &
               dt%atom_type4 == at1)) then
            exit
          end if

        end if

      else
        if (dt%func == 1 .or. dt%func == 3) then

          ! proper dihedral
          if ((dt%atom_type2 == at2 .and. &
               dt%atom_type3 == at3) .or. &
              (dt%atom_type2 == at3 .and. &
               dt%atom_type3 == at2)) then
            exit
          end if

        else ! if (dt%func == 2)

          ! improper dihedral
          if ((dt%atom_type1 == at1 .and. &
               dt%atom_type4 == at4) .or. &
              (dt%atom_type1 == at4 .and. &
               dt%atom_type4 == at1)) then
            exit
          end if

        end if
      end if
    end do

    if (i <= size(dihetypes)) then
      dihe%ps           = dt%ps
      dihe%kp           = dt%kp
      dihe%multiplicity = dt%multiplicity
      dihe%c(1:6)       = dt%c(1:6)
      search_dihe_param = .true.
    else
      search_dihe_param = .false.
    end if

    return

  end function search_dihe_param

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      search_dihe_param_duplicate
  !> @brief        search dihedral parameters with duplicate manner
  !! @authors      CK
  !! @param[in]    dihetypes : data of dihedraltypes
  !! @param[in]    at1       : atom 1 type name
  !! @param[in]    at2       : atom 2 type name
  !! @param[in]    at3       : atom 3 type name
  !! @param[in]    at4       : atom 4 type name
  !! @param[inout] dihes     : array of dihedral info.
  !! @param[inout] ind       : index of dihedral info. array
  !! @return       flag for found or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function search_dihe_param_duplicate(old, dihetypes,                       &
                                       at1, at2, at3, at4, dihes, ind)

    ! return value
    logical :: search_dihe_param_duplicate

    ! formal arguments
    integer,                  intent(in)    :: old
    type(s_dihetype), target, intent(in)    :: dihetypes(:)
    character(*),             intent(in)    :: at1
    character(*),             intent(in)    :: at2
    character(*),             intent(in)    :: at3
    character(*),             intent(in)    :: at4
    type(s_dihe),             intent(inout) :: dihes(:)
    integer,                  intent(inout) :: ind

    ! local variables
    integer                   :: i, icn, iold
    type(s_dihetype), pointer :: dt
    integer                   :: icnt(5)
    logical                   :: four_atoms

    icn = 0
    icnt(1:5) = 0
    iold = ind
    four_atoms = .false.
    do i = 1, size(dihetypes)
      dt => dihetypes(i)

      if (dt%func /= dihes(old+iold)%func) &
        cycle

      if (dt%four_type) then

        if (dt%atom_type1 == 'X' .and. dt%atom_type4 == 'X' .and. &
            .not. four_atoms) then

          if ((dt%atom_type2 == at2 .and. &
               dt%atom_type3 == at3) .or. &
              (dt%atom_type2 == at3 .and. &
               dt%atom_type3 == at2)) then
            icn = icn+1
            icnt(icn) = i
          end if

        else if (dt%atom_type1 == 'X' .and. dt%atom_type2 == 'X' .and. &
            .not. four_atoms) then

          if ((dt%atom_type3 == at3 .and. &
               dt%atom_type4 == at4) .or. &
              (dt%atom_type3 == at2 .and. &
               dt%atom_type4 == at1)) then
            icn = icn+1
            icnt(icn) = i
          end if

        else if (dt%atom_type1 == 'X' .and. &
            .not. four_atoms) then

          if ((dt%atom_type2 == at2 .and. &
               dt%atom_type3 == at3 .and. &
               dt%atom_type4 == at4) .or. &
              (dt%atom_type2 == at3 .and. &
               dt%atom_type3 == at2 .and. &
               dt%atom_type4 == at1)) then
            icn = icn+1
            icnt(icn) = i
          end if

        else if (dt%atom_type4 == 'X' .and. &
            .not. four_atoms) then

          if ((dt%atom_type1 == at1 .and. &
               dt%atom_type2 == at2 .and. &
               dt%atom_type3 == at3) .or. &
              (dt%atom_type1 == at4 .and. &
               dt%atom_type2 == at3 .and. &
               dt%atom_type3 == at2)) then
            icn = icn+1
            icnt(icn) = i
          end if

        else

          if ((dt%atom_type1 == at1 .and. &
               dt%atom_type2 == at2 .and. &
               dt%atom_type3 == at3 .and. &
               dt%atom_type4 == at4) .or. &
              (dt%atom_type1 == at4 .and. &
               dt%atom_type2 == at3 .and. &
               dt%atom_type3 == at2 .and. &
               dt%atom_type4 == at1)) then
            four_atoms = .true.
            icn = icn+1
            icnt(icn) = i
          end if

        end if

      else
        if (dt%func == 9) then

          ! proper dihedral
          if ((dt%atom_type2 == at2 .and. &
               dt%atom_type3 == at3) .or. &
              (dt%atom_type2 == at3 .and. &
               dt%atom_type3 == at2)) then
            icn = icn+1
            icnt(icn) = i
          end if

        else if (dt%func == 4) then ! if (dt%func == 2)

          ! improper dihedral
          if ((dt%atom_type1 == at1 .and. &
               dt%atom_type4 == at4) .or. &
              (dt%atom_type1 == at4 .and. &
               dt%atom_type4 == at1)) then
            icn = icn+1
            icnt(icn) = i
          end if

        end if
      end if
      if (icn == 5) &
        call error_msg('Search_dihe_param_duplicate> '// &
                       '[dihedraltypes] too much dihedral functions')
    end do

    if (icn >= 1) then
      search_dihe_param_duplicate = .true.
    else
      search_dihe_param_duplicate = .false.
    endif

    do i = 1, icn
      dt => dihetypes(icnt(i))
      dihes(old+ind)%ps           = dt%ps
      dihes(old+ind)%kp           = dt%kp
      dihes(old+ind)%multiplicity = dt%multiplicity
      dihes(old+ind)%c(1:6)       = dt%c(1:6)
      ind = ind + 1
    end do
    ind = ind-1

    return

  end function search_dihe_param_duplicate

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      search_constr_param
  !> @brief        search constraint parameters
  !! @authors      NT
  !! @param[in]    constrtypes : data of constrainttypes
  !! @param[in]    at1         : atom 1 type name
  !! @param[in]    at2         : atom 2 type name
  !! @param[in]    func        : function type
  !! @param[out]   b0          : param b0
  !! @return       flag for found or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function search_constr_param(constrtypes, at1, at2, func, b0)

    ! return value
    logical :: search_constr_param

    ! formal arguments
    type(s_constrtype), target, intent(in)    :: constrtypes(:)
    character(*),               intent(in)    :: at1
    character(*),               intent(in)    :: at2
    integer,                    intent(in)    :: func
    real(wp),                   intent(out)   :: b0

    ! local variables
    integer                     :: i
    type(s_constrtype), pointer :: ct


    do i = 1, size(constrtypes)
      ct => constrtypes(i)
      if (ct%func /= func) &
        cycle
      if (ct%atom_type1 == at1 .and. ct%atom_type2 == at2 .or. &
          ct%atom_type1 == at2 .and. ct%atom_type2 == at1) then
        b0 = ct%b0
        search_constr_param = .true.
        return
      end if
    end do

    search_constr_param = .false.
    return

  end function search_constr_param

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      search_cmap_param
  !> @brief        search cmap parameters
  !! @authors      NT
  !! @param[in]    cmaptypes : data of constrainttypes
  !! @param[in]    at1       : atom 1 type name
  !! @param[in]    at2       : atom 2 type name
  !! @param[in]    at3       : atom 3 type name
  !! @param[in]    at4       : atom 4 type name
  !! @param[in]    at5       : atom 5 type name
  !! @param[in]    func      : function type
  !! @param[out]   map       : param map
  !! @return       flag for found or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function search_cmap_param(cmaptypes, at1, at2, at3, at4, at5, func, map)

    ! return value
    logical :: search_cmap_param

    ! formal arguments
    type(s_cmaptype), target, intent(in)    :: cmaptypes(:)
    character(*),             intent(in)    :: at1
    character(*),             intent(in)    :: at2
    character(*),             intent(in)    :: at3
    character(*),             intent(in)    :: at4
    character(*),             intent(in)    :: at5
    integer,                  intent(in)    :: func
    real(wp),                 allocatable   :: map(:,:)

    ! local variables
    integer                   :: i, nrow, ncol
    type(s_cmaptype), pointer :: ct


    do i = 1, size(cmaptypes)
      ct => cmaptypes(i)
      if (ct%func /= func) &
        cycle
      if (ct%atom_type1 == at1 .and. &
          ct%atom_type2 == at2 .and. &
          ct%atom_type3 == at3 .and. &
          ct%atom_type4 == at4 .and. &
          ct%atom_type5 == at5) then

        nrow = size(ct%map(:,1))
        ncol = size(ct%map(1,:))
        allocate(map(nrow,ncol))
        map(:,:) = ct%map(:,:)
        search_cmap_param = .true.
        return
      end if
    end do

    search_cmap_param = .false.
    return

  end function search_cmap_param

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    rename_atom_type_name
  !> @brief        rename atom type name
  !! @authors      NT
  !! @param[in]    grotop    : GROMACS TOP information
  !! @param[in]    read_name : read original name
  !! @param[out]   type_name : type name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine rename_atom_type_name(grotop, read_name, type_name)

    ! formal arguments
    type(s_grotop),  target, intent(in)    :: grotop
    character(*),            intent(in)    :: read_name
    character(*),            intent(out)   :: type_name

    ! local variables
    integer                  :: i

    type(s_atomtype), pointer :: at(:)


    if (.not. grotop%alias_atomtype) then
      type_name = read_name
      return
    end if

    at => grotop%atomtypes

    do i = 1, size(at)
      if (at(i)%alias == read_name) then
        type_name = at(i)%type_name
        return
      end if
    end do

    type_name = read_name
    return

  end subroutine rename_atom_type_name

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      match_format
  !> @brief        check the line format
  !! @authors      NT
  !! @return       flag for matched or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function match_format(line, fmt)

    ! parameters
    integer,    parameter    :: MaxTokens = 20

    ! return value
    logical                  :: match_format

    ! formal arguments
    character(*),            intent(in)    :: line
    character(*),            intent(in)    :: fmt

    ! local variables
    integer                  :: i, fcnt, pcnt
    character(20)            :: tokens(MaxTokens)


    match_format = .false.

    fcnt = len_trim(fmt)

    tokens(1:MaxTokens) = ''
    read(line,*,err=1,end=1) tokens(1:MaxTokens)
1   do pcnt = 1, MaxTokens
      if (tokens(pcnt) == '') &
        exit
    end do
    pcnt = pcnt - 1

    if (pcnt < fcnt) &
      return

    do i = 1, fcnt

      select case(fmt(i:i))

      case ('S')
        ! accept anything

      case ('N')
        if (.not. is_digit(tokens(i))) &
          return

      case default
        if (.not. is_digit(tokens(i))) &
          return

        if (fmt(i:i) /= tokens(i)) &
          return

      end select
    end do

    match_format = .true.
    return

  end function match_format

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      is_digit
  !> @brief        check digit
  !! @authors      NT
  !! @return       flag for found or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function is_digit(str)

    ! return value
    logical                  :: is_digit

    ! formal arguments
    character(*),            intent(in)    :: str

    ! local variables
    character                :: c


    c = str(1:1)
    is_digit = (c >= '0' .and. c <= '9' .or. c == '.' .or. c == '-')
    return

  end function is_digit

end module fileio_grotop_mod
