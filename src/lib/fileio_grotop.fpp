! --------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_grotop_mod
!> @brief   read GROMACS topology file
!! @authors Norio Takase (NT), Cheng Tan (CT)
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
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

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
    real(wp)                       :: v             = 0.0_wp
    real(wp)                       :: w             = 0.0_wp
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
    !CK
    real(wp)                       :: theta_01      = 0.0_wp !deg.
    real(wp)                       :: kt1           = 0.0_wp
    real(wp)                       :: cgamma        = 0.0_wp !deg.
    real(wp)                       :: epsa          = 0.0_wp
    ! shinobu-edited
    real(wp)                       :: w             = 0.0_wp
    integer                        :: types         = 0
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
    !shinobu-edited
    real(wp)                       :: w             = 0.0_wp
    real(wp)                       :: theta_0       = 0.0_wp
    integer                        :: types         = 0
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
    !shinobu-edited
    real(wp)                       :: r0            = 0.0_wp
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

  ! pair data for multi
  !    func_type = 1  : LJ
  type, public :: s_mcont
    integer                        :: atom_idx1     = 0
    integer                        :: atom_idx2     = 0
    integer                        :: func          = 0
    integer                        :: model         = 0
    real(wp)                       :: v             = 0.0_wp
    real(wp)                       :: w             = 0.0_wp
  end type s_mcont

  ! pair data for multi
  !    func_type = 1  : LJ
  type, public :: s_morph_pair
    integer                        :: atom_idx1     = 0
    integer                        :: atom_idx2     = 0
    real(wp)                       :: rmin          = 0.0_wp
    real(wp)                       :: rmin_other    = 0.0_wp
  end type s_morph_pair

  ! CG: PWMcos pairs
  type, public :: s_pwmcos
    integer                        :: func          = 0
    integer                        :: protein_idx   = 0
    real(wp)                       :: r0            = 0.0_wp
    real(wp)                       :: theta1        = 0.0_wp
    real(wp)                       :: theta2        = 0.0_wp
    real(wp)                       :: theta3        = 0.0_wp
    real(wp)                       :: ene_A         = 0.0_wp
    real(wp)                       :: ene_C         = 0.0_wp
    real(wp)                       :: ene_G         = 0.0_wp
    real(wp)                       :: ene_T         = 0.0_wp
    real(wp)                       :: gamma         = 0.0_wp
    real(wp)                       :: eps_shift     = 0.0_wp
  end type s_pwmcos

  ! CG: PWMcos-ns pairs
  type, public :: s_pwmcosns
    integer                        :: func          = 0
    integer                        :: protein_idx   = 0
    real(wp)                       :: r0            = 0.0_wp
    real(wp)                       :: theta1        = 0.0_wp
    real(wp)                       :: theta2        = 0.0_wp
    real(wp)                       :: ene           = 0.0_wp
  end type s_pwmcosns

  ! ~CG~ protein IDR HPS region
  !
  type, public :: s_idr_hps
    integer                        :: grp_start  = 0
    integer                        :: grp_end    = 0
  end type s_idr_hps

  ! ~CG~ protein IDR KH region
  !
  type, public :: s_idr_kh
    integer                        :: grp_start  = 0
    integer                        :: grp_end    = 0
  end type s_idr_kh

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
    integer                        :: num_mcontacts = 0
    integer                        :: num_morph_bb  = 0
    integer                        :: num_morph_sc  = 0
    integer                        :: num_pwmcos    = 0
    integer                        :: num_pwmcosns  = 0
    integer                        :: num_idr_hps   = 0
    integer                        :: num_idr_kh   = 0
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
    type(s_mcont),     allocatable :: mcontact(:)
    type(s_morph_pair),allocatable :: morph_bb(:)
    type(s_morph_pair),allocatable :: morph_sc(:)
    type(s_pwmcos),    allocatable :: pwmcos(:)
    type(s_pwmcosns),  allocatable :: pwmcosns(:)
    type(s_idr_hps),   allocatable :: idr_hps(:)
    type(s_idr_kh),    allocatable :: idr_kh(:)
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
    integer                        :: wild_num
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

  !Ck made
  ! gomodel data
  type, public :: s_gomodel
    integer                        :: go_func         = 1
    integer                        :: num_basins      = 1
  end type s_gomodel

  !Ck made
  ! gomodel data
  type, public :: s_itpread
    character(MaxLine)             :: itp             = ''
  end type s_itpread

  ! Ezmembrane type data
  type, public :: s_ezmembrane
    character(6)                   :: atom_type    = ''
    integer                        :: func         = 0
    real(wp)                       :: de           = 0.0_wp
    real(wp)                       :: zmin         = 0.0_wp
    real(wp)                       :: polym        = 0.0_wp
  end type s_ezmembrane

  ! Ezmembrane type data
  type, public :: s_flangltype
    character(6)                   :: atom_type    = ''
    integer                        :: func         = 0
    real(wp),          allocatable :: theta(:)
    real(wp),          allocatable :: efunc(:)
    real(wp),          allocatable :: d2func(:)
  end type s_flangltype

  type, public :: s_fldihetype
    character(6)                   :: atom_type1   = ''
    character(6)                   :: atom_type2   = ''
    integer                        :: func         = 0
    real(wp),          allocatable :: coef(:)
  end type s_fldihetype

  ! ~CG~ 3SPN.2C DNA: intra-strand base stacking
  !
  type, public :: s_basestacktype
    character(6)                   :: base_type5 = '' ! 5' Base type
    character(6)                   :: base_type3 = '' ! 3' Base type
    integer                        :: func       = 0
    real(wp)                       :: epsilon    = 0.0_wp
    real(wp)                       :: sigma      = 0.0_wp
    real(wp)                       :: theta_bs   = 0.0_wp ! degree
  end type s_basestacktype

  ! ~CG~ 3SPN.2C DNA: intra/inter-strand base pairing
  !
  type, public :: s_basepairtype
    character(6)                   :: base_type_a = ''
    character(6)                   :: base_type_b = ''
    integer                        :: func        = 0
    real(wp)                       :: theta_1     = 0.0_wp ! degree
    real(wp)                       :: theta_2     = 0.0_wp ! degree
    real(wp)                       :: theta_3     = 0.0_wp ! degree
    real(wp)                       :: phi_1       = 0.0_wp ! degree
    real(wp)                       :: sigma       = 0.0_wp ! nm
    real(wp)                       :: epsilon     = 0.0_wp ! kJ/mol
  end type s_basepairtype

  ! ~CG~ 3SPN.2C DNA: base cross-stacking
  !
  type, public :: s_basecrosstype
    character(6)                   :: base_type_a = ''
    character(6)                   :: base_type_b = ''
    integer                        :: func        = 0
    real(wp)                       :: epsilon     = 0.0_wp
    real(wp)                       :: sigma       = 0.0_wp
    real(wp)                       :: theta_cs    = 0.0_wp ! degree
  end type s_basecrosstype

  ! ~CG~ 3SPN.2C DNA: excluded volume
  !
  type, public :: s_cgdnaexvtype
    character(6)                   :: base_type   = ''
    integer                        :: func        = 0
    real(wp)                       :: sigma       = 0.0_wp
  end type s_cgdnaexvtype

  ! ~CG~ debye-huckel ele interaction mol-mol pairs
  !
  type, public :: s_cg_ele_mol_pair
    integer                        :: func        = -1
    integer                        :: grp1_start  = 0
    integer                        :: grp1_end    = 0
    integer                        :: grp2_start  = 0
    integer                        :: grp2_end    = 0
    logical                        :: is_intermol = .true.
  end type s_cg_ele_mol_pair

  ! ~CG~ PWMcos mol-mol pairs
  !
  type, public :: s_pwmcos_mol_pair
    integer                        :: func        = -1
    integer                        :: grp1_start  = 0
    integer                        :: grp1_end    = 0
    integer                        :: grp2_start  = 0
    integer                        :: grp2_end    = 0
  end type s_pwmcos_mol_pair

  ! ~CG~ PWMcos-ns mol-mol pairs
  !
  type, public :: s_pwmcosns_mol_pair
    integer                        :: func        = -1
    integer                        :: grp1_start  = 0
    integer                        :: grp1_end    = 0
    integer                        :: grp2_start  = 0
    integer                        :: grp2_end    = 0
  end type s_pwmcosns_mol_pair

  ! ~CG~ KH mol-mol pairs
  !
  type, public :: s_cg_kh_mol_pair
    integer                        :: func        = -1
    integer                        :: grp1_start  = 0
    integer                        :: grp1_end    = 0
    integer                        :: grp2_start  = 0
    integer                        :: grp2_end    = 0
    logical                        :: is_intermol = .true.
  end type s_cg_kh_mol_pair

  ! ~CG~ protein IDR HPS model params
  !
  type, public :: s_cg_IDR_HPS_atomtype
    character(6)                   :: type_name    = ''
    real(wp)                       :: mass         = 0.0_wp
    real(wp)                       :: charge       = 0.0_wp
    real(wp)                       :: sigma        = 0.0_wp
    real(wp)                       :: lambda       = 0.0_wp
  end type s_cg_IDR_HPS_atomtype

  ! ~CG~ protein KH model params
  !
  type, public :: s_cg_KH_atomtype
    character(6)                   :: type_name    = ''
    real(wp)                       :: mass         = 0.0_wp
    real(wp)                       :: charge       = 0.0_wp
    real(wp)                       :: sigma        = 0.0_wp
  end type s_cg_KH_atomtype

  ! ~CG~ protein-protein Miyazawa-Jernigan model epsilons
  !
  type, public :: s_cg_pair_MJ_eps
    character(6)                   :: type_name_1    = ''
    character(6)                   :: type_name_2    = ''
    real(wp)                       :: epsilon        = 0.0_wp
  end type s_cg_pair_MJ_eps

  ! topology data
  type, public :: s_grotop
    integer                               :: num_atomtypes      = 0
    integer                               :: num_bondtypes      = 0
    integer                               :: num_angltypes      = 0
    integer                               :: num_flangltypes    = 0
    integer                               :: num_dihetypes      = 0
    integer                               :: num_fldihetypes    = 0
    integer                               :: num_constrtypes    = 0
    integer                               :: num_cmaptypes      = 0
    integer                               :: num_nbonparms      = 0
    integer                               :: num_moltypes       = 0
    integer                               :: num_molss          = 0
    integer                               :: num_membranetypes  = 0
    integer                               :: num_basestacktypes = 0
    integer                               :: num_basepairtypes  = 0
    integer                               :: num_basecrosstypes = 0
    integer                               :: num_cgdnaexvtypes  = 0
    integer                               :: num_cgelemolpairs  = 0
    integer                               :: num_cgkhmolpairs  = 0
    integer                               :: num_pwmcosmolpairs = 0
    integer                               :: num_pwmcosnsmolpairs = 0
    integer                               :: num_cg_IDR_HPS_atomtypes = 0
    integer                               :: num_cg_KH_atomtypes = 0
    integer                               :: num_cg_pair_MJ_eps = 0
    logical                               :: alias_atomtype     = .false.
    type(s_defaults)                      :: defaults
    type(s_atomtype),         allocatable :: atomtypes(:)
    type(s_bondtype),         allocatable :: bondtypes(:)
    type(s_angltype),         allocatable :: angltypes(:)
    type(s_dihetype),         allocatable :: dihetypes(:)
    type(s_constrtype),       allocatable :: constrtypes(:)
    type(s_cmaptype),         allocatable :: cmaptypes(:)
    type(s_nbonparm),         allocatable :: nbonparms(:)
    type(s_moltype),          allocatable :: moltypes(:)
    type(s_mols),             allocatable :: molss(:)
    type(s_ezmembrane),       allocatable :: ez_membrane(:)
    type(s_flangltype),       allocatable :: flangltypes(:)
    type(s_fldihetype),       allocatable :: fldihetypes(:)
    type(s_basestacktype),    allocatable :: basestacktypes(:)
    type(s_basepairtype),     allocatable :: basepairtypes(:)
    type(s_basecrosstype),    allocatable :: basecrosstypes(:)
    type(s_cgdnaexvtype),     allocatable :: cgdnaexvtypes(:)
    type(s_cg_ele_mol_pair),  allocatable :: cg_ele_mol_pairs(:)
    type(s_pwmcos_mol_pair),  allocatable :: pwmcos_mol_pairs(:)
    type(s_pwmcosns_mol_pair),  allocatable :: pwmcosns_mol_pairs(:)
    type(s_cg_kh_mol_pair),   allocatable :: cg_kh_mol_pairs(:)
    type(s_cg_IDR_HPS_atomtype),  allocatable :: cg_IDR_HPS_atomtypes(:)
    type(s_cg_KH_atomtype),   allocatable :: cg_KH_atomtypes(:)

    type(s_cg_pair_MJ_eps),   allocatable :: cg_pair_MJ_eps(:)
    type(s_gomodel)                       :: gomodel
    type(s_itpread)                       :: itpread
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
  integer,      public, parameter :: GMGroMolMcont    = 114
  integer,      public, parameter :: GMGroMolMorphBB  = 115
  integer,      public, parameter :: GMGroMolMorphSC  = 116
  integer,      public, parameter :: GMGroMolPWMcos   = 117
  integer,      public, parameter :: GMGroMolPWMcosns = 118
  integer,      public, parameter :: GMGroMolIDRHPS   = 119
  integer,      public, parameter :: GMGroMolIDRKH    = 120

  ! parameters for TOP structure allocatable variables
  integer,      public, parameter :: GroTopAtomType      = 1
  integer,      public, parameter :: GroTopBondType      = 2
  integer,      public, parameter :: GroTopAnglType      = 3
  integer,      public, parameter :: GroTopDiheType      = 4
  integer,      public, parameter :: GroTopConstrType    = 5
  integer,      public, parameter :: GroTopCmapType      = 6
  integer,      public, parameter :: GroTopNbonParm      = 7
  integer,      public, parameter :: GroTopMolType       = 8
  integer,      public, parameter :: GroTopMols          = 9
  integer,      public, parameter :: GroTopEzMembrane    = 10
  integer,      public, parameter :: GroTopFlAnglType    = 11
  integer,      public, parameter :: GroTopBaseStackType = 12
  integer,      public, parameter :: GroTopBasePairType  = 13
  integer,      public, parameter :: GroTopBaseCrossType = 14
  integer,      public, parameter :: GroTopCGDNAExvType  = 15
  integer,      public, parameter :: GroTopFlDiheType    = 16
  integer,      public, parameter :: GroTopCGEleMolPair  = 17
  integer,      public, parameter :: GroTopPWMcosMolPair = 18
  integer,      public, parameter :: GroTopPWMcosnsMolPair  = 19
  integer,      public, parameter :: GroTopCGIDRHPSAtomType = 20
  integer,      public, parameter :: GroTopCGKHAtomType  = 21
  integer,      public, parameter :: GroTopCGPairMJEpsilon  = 22
  integer,      public, parameter :: GroTopCGKHMolPair   = 23

  ! parameters for directive
  integer,     private, parameter :: DDefaults             = 1
  integer,     private, parameter :: DAtomTypes            = 2
  integer,     private, parameter :: DBondTypes            = 3
  integer,     private, parameter :: DAngleTypes           = 4
  integer,     private, parameter :: DDihedralTypes        = 5
  integer,     private, parameter :: DConstrTypes          = 6
  integer,     private, parameter :: DCmapTypes            = 7
  integer,     private, parameter :: DNonbondParams        = 8
  !    molecule directive
  integer,     private, parameter :: DMoleculeType         = 9
  integer,     private, parameter :: DAtoms                = 10
  integer,     private, parameter :: DBonds                = 11
  integer,     private, parameter :: DAngles               = 12
  integer,     private, parameter :: DDihedrals            = 13
  integer,     private, parameter :: DCmap                 = 14
  integer,     private, parameter :: DExclusions           = 15
  integer,     private, parameter :: DConstraints          = 16
  integer,     private, parameter :: DPairs                = 17
  integer,     private, parameter :: DSettles              = 18
  integer,     private, parameter :: DVirtualSites2        = 19
  integer,     private, parameter :: DVirtualSites3        = 20
  integer,     private, parameter :: DVirtualSites4        = 21
  integer,     private, parameter :: DVirtualSitesN        = 22
  integer,     private, parameter :: DPosition_Restraints  = 23
  !    system directive
  integer,     private, parameter :: DSystem               = 24
  integer,     private, parameter :: DMolecules            = 25
  !    CK made directive
  integer,     private, parameter :: DGomodel              = 26
  integer,     private, parameter :: DCg_Atoms             = 27
  integer,     private, parameter :: DMulticontact         = 28
  integer,     private, parameter :: DMorphBB              = 29
  integer,     private, parameter :: DMorphSC              = 30
  integer,     private, parameter :: DEzMembrane           = 31
  integer,     private, parameter :: DFlexible_Local_Angle = 32
  integer,     private, parameter :: DFlexible_Local_Dihedral_Angle = 33
  !    unknown directive
  integer,     private, parameter :: DUnknown              = 34
  !    3SPN.2C DNA interaction type directive
  integer,     private, parameter :: DBaseStackTypes       = 35
  integer,     private, parameter :: DBasePairTypes        = 36
  integer,     private, parameter :: DBaseCrossTypes       = 37
  integer,     private, parameter :: DCGDNAExvTypes        = 38
  integer,     private, parameter :: DCGEleMolPairs        = 39
  integer,     private, parameter :: DCGPWMcos             = 40
  integer,     private, parameter :: DCGPWMcosMolPairs     = 41
  integer,     private, parameter :: DCGPWMcosns           = 42
  integer,     private, parameter :: DCGPWMcosnsMolPairs   = 43
  integer,     private, parameter :: DCGIDRHPSAtomTypes    = 44
  integer,     private, parameter :: DCGIDRHPSRegion       = 45
  integer,     private, parameter :: DCGKHAtomTypes        = 46
  integer,     private, parameter :: DCGIDRKHRegion        = 47
  integer,     private, parameter :: DCGPairMJEpsilon      = 48
  integer,     private, parameter :: DCGKHMolPairs         = 49

  ! parameters
  logical,     private, parameter :: VerboseOn        = .false.

  ! for duplicated dihedrals
  integer,     private, parameter :: MAX_MULTIPLY_NUM = 15

  ! local variables
  logical,                private :: vervose = .true.

  ! subroutines
  public  :: input_grotop
  public  :: output_grotop
  public  :: output_grotop_aadome

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
  private :: read_flangle_types
  private :: read_dihedral_types
  private :: read_fldihe_types
  private :: read_constraint_types
  private :: read_cmap_types
  private :: read_base_stack_types
  private :: read_base_pair_types
  private :: read_base_cross_types
  private :: read_base_exv_types
  private :: read_cg_ele_mol_pairs
  private :: read_cg_kh_mol_pairs
  private :: read_pwmcos_mol_pairs
  private :: read_pwmcosns_mol_pairs
  private :: read_cg_IDR_HPS_atom_types
  private :: read_cg_KH_atom_types
  private :: read_cg_pair_MJ_eps
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
  private :: read_gomodel
  private :: read_cgatoms
  private :: read_multicontact
  private :: read_morph_bb
  private :: read_morph_sc
  private :: read_membrane_type
  private :: read_pwmcos
  private :: read_pwmcosns
  private :: read_idr_hps
  private :: read_idr_kh

  private :: write_grotop
  private :: write_defaults
  private :: write_atom_types
  private :: write_bond_types
  private :: write_bond_types_aa
  private :: write_angle_types
  private :: write_angle_types_aa
  private :: write_dihedral_types
  private :: write_dihedral_types_aa
  private :: write_constraint_types
  private :: write_cmap_types
  private :: write_base_stack_types
  private :: write_base_pair_types
  private :: write_base_cross_types
  private :: write_base_exv_types
  private :: write_cg_IDR_HPS_atom_types
  private :: write_cg_KH_atom_types
  private :: write_nonbond_parms
  private :: write_molecule_type
  private :: write_molecule_type_aa
  private :: write_atoms
  private :: write_bonds
  private :: write_bonds_aa
  private :: write_angles
  private :: write_angles_aa
  private :: write_dihedrals
  private :: write_dihedrals_aa
  private :: write_cmaps
  private :: write_exclusions
  private :: write_constraints
  private :: write_pairs
  private :: write_multicontact
  private :: write_morph_pairs_bb
  private :: write_morph_pairs_sc
  private :: write_pwmcos
  private :: write_pwmcosns
  private :: write_settles
  private :: write_virtual_sites2
  private :: write_virtual_sites3
  private :: write_virtual_sites4
  private :: write_virtual_sitesn
  private :: write_position_restraints
  private :: write_system
  private :: write_gomodel
  private :: write_molecules

  private :: check_section_count
  private :: check_directive
  private :: size_grotop
  private :: size_grotop_mol
  private :: error_msg_grotop
  private :: search_atom_mass
  private :: search_atom_vw
  private :: search_stokes_r
  private :: search_atom_charge
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
    integer                  :: ierr
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
    if (file == 0) &
      call error_msg('Input_Grotop> '//trim(error))

    ! read GROMACS TOP file
    !
    call read_grotop(file, grotop)


    ! close GROMACS TOP file
    !
#ifdef HAVE_MPI_GENESIS
    if (nproc_country > 1) then
      call mpi_barrier(mpi_comm_world, ierr) ! barrier is required to unlink.
    endif
#endif
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
  !  Subroutine    output_grotop_aadome
  !> @brief        a driver subroutine for writing GRO TOP file
  !! @authors      CK
  !! @param[in]    grotop_filename : filename of GRO TOP file
  !! @param[in]    grotop          : structure of GRO TOP data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_grotop_aadome(grotop_filename, grotop)

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
    call write_grotop_aadome(file, grotop)


    ! close GROMACS TOP file
    !
    call close_file(file)

    return

  end subroutine output_grotop_aadome


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


    grotop%num_atomtypes      = 0
    grotop%num_bondtypes      = 0
    grotop%num_angltypes      = 0
    grotop%num_flangltypes    = 0
    grotop%num_dihetypes      = 0
    grotop%num_fldihetypes    = 0
    grotop%num_constrtypes    = 0
    grotop%num_cmaptypes      = 0
    grotop%num_nbonparms      = 0
    grotop%num_moltypes       = 0
    grotop%num_molss          = 0
    grotop%num_membranetypes  = 0
    grotop%num_basestacktypes = 0
    grotop%num_basepairtypes  = 0
    grotop%num_basecrosstypes = 0
    grotop%num_cgdnaexvtypes  = 0
    grotop%num_cgelemolpairs  = 0
    grotop%num_cgkhmolpairs   = 0
    grotop%num_pwmcosmolpairs = 0
    grotop%num_pwmcosnsmolpairs = 0
    grotop%num_cg_IDR_HPS_atomtypes = 0
    grotop%num_cg_KH_atomtypes= 0
    grotop%num_cg_pair_MJ_eps = 0
    grotop%alias_atomtype     = .false.

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


    gromol%num_atoms     = 0
    gromol%num_bonds     = 0
    gromol%num_angls     = 0
    gromol%num_dihes     = 0
    gromol%num_cmaps     = 0
    gromol%num_excls     = 0
    gromol%num_constrs   = 0
    gromol%num_pairs     = 0
    gromol%num_vsites2   = 0
    gromol%num_vsites3   = 0
    gromol%num_vsites4   = 0
    gromol%num_vsitesn   = 0
    gromol%num_posress   = 0
    gromol%num_mcontacts = 0
    gromol%num_morph_bb  = 0
    gromol%num_morph_sc  = 0
    gromol%num_pwmcos    = 0
    gromol%num_pwmcosns  = 0
    gromol%num_idr_hps   = 0
    gromol%num_idr_kh    = 0

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

    type(s_atomtype),         allocatable :: atomtypes(:)
    type(s_bondtype),         allocatable :: bondtypes(:)
    type(s_angltype),         allocatable :: angltypes(:)
    type(s_flangltype),       allocatable :: flangltypes(:)
    type(s_dihetype),         allocatable :: dihetypes(:)
    type(s_fldihetype),       allocatable :: fldihetypes(:)
    type(s_constrtype),       allocatable :: constrtypes(:)
    type(s_cmaptype),         allocatable :: cmaptypes(:)
    type(s_nbonparm),         allocatable :: nbonparms(:)
    type(s_moltype),          allocatable :: moltypes(:)
    type(s_mols),             allocatable :: molss(:)
    type(s_ezmembrane),       allocatable :: ez_membrane(:)
    type(s_basestacktype),    allocatable :: basestacktypes(:)
    type(s_basepairtype),     allocatable :: basepairtypes(:)
    type(s_basecrosstype),    allocatable :: basecrosstypes(:)
    type(s_cgdnaexvtype),     allocatable :: cgdnaexvtypes(:)
    type(s_cg_ele_mol_pair),  allocatable :: cg_ele_mol_pairs(:)
    type(s_cg_kh_mol_pair),   allocatable :: cg_kh_mol_pairs(:)
    type(s_pwmcos_mol_pair),  allocatable :: pwmcos_mol_pairs(:)
    type(s_pwmcosns_mol_pair),allocatable :: pwmcosns_mol_pairs(:)
    type(s_cg_IDR_HPS_atomtype), allocatable :: cg_IDR_HPS_atomtypes(:)
    type(s_cg_KH_atomtype),   allocatable :: cg_KH_atomtypes(:)
    type(s_cg_pair_MJ_eps),   allocatable :: cg_pair_MJ_eps(:)


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

    case(GroTopFlAnglType)
      if (allocated(grotop%flangltypes)) then

        old_size = size(grotop%flangltypes)
        if (old_size == var_size) &
          return
        allocate(flangltypes(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        flangltypes(:) = grotop%flangltypes(:)
        deallocate(grotop%flangltypes, stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_dealloc
        allocate(grotop%flangltypes(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        if (old_size > var_size) then
          grotop%flangltypes(1:var_size) = flangltypes(1:var_size)
        else
          grotop%flangltypes(1:old_size) = flangltypes(1:old_size)
        end if
        deallocate(flangltypes, stat = alloc_stat)
      else
        old_size = 0
        allocate(grotop%flangltypes(var_size), stat = alloc_stat)
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

    case(GroTopFlDiheType)
      if (allocated(grotop%fldihetypes)) then

        old_size = size(grotop%fldihetypes)
        if (old_size == var_size) &
          return
        allocate(fldihetypes(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        fldihetypes(:) = grotop%fldihetypes(:)
        deallocate(grotop%fldihetypes, stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_dealloc
        allocate(grotop%fldihetypes(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        if (old_size > var_size) then
          grotop%fldihetypes(1:var_size) = fldihetypes(1:var_size)
        else
          grotop%fldihetypes(1:old_size) = fldihetypes(1:old_size)
        end if
        deallocate(fldihetypes, stat = alloc_stat)
      else
        old_size = 0
        allocate(grotop%fldihetypes(var_size), stat = alloc_stat)
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

    case(GroTopEzMembrane)
      if (allocated(grotop%ez_membrane)) then

        old_size = size(grotop%ez_membrane)
        if (old_size == var_size) &
          return
        allocate(ez_membrane(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        ez_membrane(:) = grotop%ez_membrane(:)
        deallocate(grotop%ez_membrane, stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_dealloc
        allocate(grotop%ez_membrane(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        if (old_size > var_size) then
          grotop%ez_membrane(1:var_size) = ez_membrane(1:var_size)
        else
          grotop%ez_membrane(1:old_size) = ez_membrane(1:old_size)
        end if
        deallocate(molss, stat = alloc_stat)
      else
        old_size = 0
        allocate(grotop%ez_membrane(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
      end if

    ! ~CG~ 3SPN.2C DNA: base-base interaction types...
    case(GroTopBaseStackType)
      if (allocated(grotop%basestacktypes)) then

        old_size = size(grotop%basestacktypes)
        if (old_size == var_size) &
              return
        allocate(basestacktypes(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
        basestacktypes(:) = grotop%basestacktypes(:)
        deallocate(grotop%basestacktypes, stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_dealloc
        allocate(grotop%basestacktypes(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
        if (old_size > var_size) then
          grotop%basestacktypes(1:var_size) = basestacktypes(1:var_size)
        else
          grotop%basestacktypes(1:old_size) = basestacktypes(1:old_size)
        end if
        deallocate(basestacktypes, stat = alloc_stat)
      else
        old_size = 0
        allocate(grotop%basestacktypes(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
      end if

    case(GroTopBasePairType)
      if (allocated(grotop%basepairtypes)) then

        old_size = size(grotop%basepairtypes)
        if (old_size == var_size) &
              return
        allocate(basepairtypes(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
        basepairtypes(:) = grotop%basepairtypes(:)
        deallocate(grotop%basepairtypes, stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_dealloc
        allocate(grotop%basepairtypes(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
        if (old_size > var_size) then
          grotop%basepairtypes(1:var_size) = basepairtypes(1:var_size)
        else
          grotop%basepairtypes(1:old_size) = basepairtypes(1:old_size)
        end if
        deallocate(basepairtypes, stat = alloc_stat)
      else
        old_size = 0
        allocate(grotop%basepairtypes(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
      end if

    case(GroTopBaseCrossType)
      if (allocated(grotop%basecrosstypes)) then

        old_size = size(grotop%basecrosstypes)
        if (old_size == var_size) &
              return
        allocate(basecrosstypes(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
        basecrosstypes(:) = grotop%basecrosstypes(:)
        deallocate(grotop%basecrosstypes, stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_dealloc
        allocate(grotop%basecrosstypes(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
        if (old_size > var_size) then
          grotop%basecrosstypes(1:var_size) = basecrosstypes(1:var_size)
        else
          grotop%basecrosstypes(1:old_size) = basecrosstypes(1:old_size)
        end if
        deallocate(basecrosstypes, stat = alloc_stat)
      else
        old_size = 0
        allocate(grotop%basecrosstypes(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
      end if

    case(GroTopCGDNAExvType)
      if (allocated(grotop%cgdnaexvtypes)) then

        old_size = size(grotop%cgdnaexvtypes)
        if (old_size == var_size) &
              return
        allocate(cgdnaexvtypes(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
        cgdnaexvtypes(:) = grotop%cgdnaexvtypes(:)
        deallocate(grotop%cgdnaexvtypes, stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_dealloc
        allocate(grotop%cgdnaexvtypes(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
        if (old_size > var_size) then
          grotop%cgdnaexvtypes(1:var_size) = cgdnaexvtypes(1:var_size)
        else
          grotop%cgdnaexvtypes(1:old_size) = cgdnaexvtypes(1:old_size)
        end if
        deallocate(cgdnaexvtypes, stat = alloc_stat)
      else
        old_size = 0
        allocate(grotop%cgdnaexvtypes(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
      end if

    case(GroTopCGEleMolPair)
      if (allocated(grotop%cg_ele_mol_pairs)) then

        old_size = size(grotop%cg_ele_mol_pairs)
        if (old_size == var_size) &
              return
        allocate(cg_ele_mol_pairs(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
        cg_ele_mol_pairs(:) = grotop%cg_ele_mol_pairs(:)
        deallocate(grotop%cg_ele_mol_pairs, stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_dealloc
        allocate(grotop%cg_ele_mol_pairs(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
        if (old_size > var_size) then
          grotop%cg_ele_mol_pairs(1:var_size) = cg_ele_mol_pairs(1:var_size)
        else
          grotop%cg_ele_mol_pairs(1:old_size) = cg_ele_mol_pairs(1:old_size)
        end if
        deallocate(cg_ele_mol_pairs, stat = alloc_stat)
      else
        old_size = 0
        allocate(grotop%cg_ele_mol_pairs(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
      end if

    case(GroTopCGKHMolPair)
      if (allocated(grotop%cg_kh_mol_pairs)) then

        old_size = size(grotop%cg_kh_mol_pairs)
        if (old_size == var_size) &
              return
        allocate(cg_kh_mol_pairs(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
        cg_kh_mol_pairs(:) = grotop%cg_kh_mol_pairs(:)
        deallocate(grotop%cg_kh_mol_pairs, stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_dealloc
        allocate(grotop%cg_kh_mol_pairs(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
        if (old_size > var_size) then
          grotop%cg_kh_mol_pairs(1:var_size) = cg_kh_mol_pairs(1:var_size)
        else
          grotop%cg_kh_mol_pairs(1:old_size) = cg_kh_mol_pairs(1:old_size)
        end if
        deallocate(cg_kh_mol_pairs, stat = alloc_stat)
      else
        old_size = 0
        allocate(grotop%cg_kh_mol_pairs(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
      end if

    case(GroTopPWMcosMolPair)
      if (allocated(grotop%pwmcos_mol_pairs)) then

        old_size = size(grotop%pwmcos_mol_pairs)
        if (old_size == var_size) &
              return
        allocate(pwmcos_mol_pairs(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
        pwmcos_mol_pairs(:) = grotop%pwmcos_mol_pairs(:)
        deallocate(grotop%pwmcos_mol_pairs, stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_dealloc
        allocate(grotop%pwmcos_mol_pairs(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
        if (old_size > var_size) then
          grotop%pwmcos_mol_pairs(1:var_size) = pwmcos_mol_pairs(1:var_size)
        else
          grotop%pwmcos_mol_pairs(1:old_size) = pwmcos_mol_pairs(1:old_size)
        end if
        deallocate(pwmcos_mol_pairs, stat = alloc_stat)
      else
        old_size = 0
        allocate(grotop%pwmcos_mol_pairs(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
      end if

    case(GroTopPWMcosnsMolPair)
      if (allocated(grotop%pwmcosns_mol_pairs)) then

        old_size = size(grotop%pwmcosns_mol_pairs)
        if (old_size == var_size) &
              return
        allocate(pwmcosns_mol_pairs(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
        pwmcosns_mol_pairs(:) = grotop%pwmcosns_mol_pairs(:)
        deallocate(grotop%pwmcosns_mol_pairs, stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_dealloc
        allocate(grotop%pwmcosns_mol_pairs(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
        if (old_size > var_size) then
          grotop%pwmcosns_mol_pairs(1:var_size) = pwmcosns_mol_pairs(1:var_size)
        else
          grotop%pwmcosns_mol_pairs(1:old_size) = pwmcosns_mol_pairs(1:old_size)
        end if
        deallocate(pwmcosns_mol_pairs, stat = alloc_stat)
      else
        old_size = 0
        allocate(grotop%pwmcosns_mol_pairs(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
      end if

    case(GroTopCGIDRHPSAtomType)
      if (allocated(grotop%cg_IDR_HPS_atomtypes)) then

        old_size = size(grotop%cg_IDR_HPS_atomtypes)
        if (old_size == var_size) &
              return
        allocate(cg_IDR_HPS_atomtypes(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
        cg_IDR_HPS_atomtypes(:) = grotop%cg_IDR_HPS_atomtypes(:)
        deallocate(grotop%cg_IDR_HPS_atomtypes, stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_dealloc
        allocate(grotop%cg_IDR_HPS_atomtypes(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
        if (old_size > var_size) then
          grotop%cg_IDR_HPS_atomtypes(1:var_size) = cg_IDR_HPS_atomtypes(1:var_size)
        else
          grotop%cg_IDR_HPS_atomtypes(1:old_size) = cg_IDR_HPS_atomtypes(1:old_size)
        end if
        deallocate(cg_IDR_HPS_atomtypes, stat = alloc_stat)
      else
        old_size = 0
        allocate(grotop%cg_IDR_HPS_atomtypes(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
      end if

    case(GroTopCGKHAtomType)
      if (allocated(grotop%cg_KH_atomtypes)) then

        old_size = size(grotop%cg_KH_atomtypes)
        if (old_size == var_size) &
              return
        allocate(cg_KH_atomtypes(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
        cg_KH_atomtypes(:) = grotop%cg_KH_atomtypes(:)
        deallocate(grotop%cg_KH_atomtypes, stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_dealloc
        allocate(grotop%cg_KH_atomtypes(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
        if (old_size > var_size) then
          grotop%cg_KH_atomtypes(1:var_size) = cg_KH_atomtypes(1:var_size)
        else
          grotop%cg_KH_atomtypes(1:old_size) = cg_KH_atomtypes(1:old_size)
        end if
        deallocate(cg_KH_atomtypes, stat = alloc_stat)
      else
        old_size = 0
        allocate(grotop%cg_KH_atomtypes(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
      end if

    case(GroTopCGPairMJEpsilon)
      if (allocated(grotop%cg_pair_MJ_eps)) then

        old_size = size(grotop%cg_pair_MJ_eps)
        if (old_size == var_size) return

        allocate(cg_pair_MJ_eps(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) call error_msg_alloc
        cg_pair_MJ_eps(:) = grotop%cg_pair_MJ_eps(:)

        deallocate(grotop%cg_pair_MJ_eps, stat = alloc_stat)
        if (alloc_stat /= 0) call error_msg_dealloc
        allocate(grotop%cg_pair_MJ_eps(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) call error_msg_alloc

        if (old_size > var_size) then
          grotop%cg_pair_MJ_eps(1:var_size) = cg_pair_MJ_eps(1:var_size)
        else
          grotop%cg_pair_MJ_eps(1:old_size) = cg_pair_MJ_eps(1:old_size)
        end if

        deallocate(cg_pair_MJ_eps, stat = alloc_stat)

      else

        old_size = 0
        allocate(grotop%cg_pair_MJ_eps(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) call error_msg_alloc

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
    type(s_mcont),           allocatable   :: mcontact(:)
    type(s_morph_pair),      allocatable   :: morph_bb(:)
    type(s_morph_pair),      allocatable   :: morph_sc(:)
    type(s_pwmcos),          allocatable   :: pwmcos(:)
    type(s_pwmcosns),        allocatable   :: pwmcosns(:)
    type(s_idr_hps),         allocatable   :: idr_hps(:)
    type(s_idr_kh),          allocatable   :: idr_kh(:)

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

    case(GMGroMolMcont)
      if (allocated(gromol%mcontact)) then

        old_size = size(gromol%mcontact)
        if (old_size == var_size) &
          return
        allocate(mcontact(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        mcontact(:) = gromol%mcontact(:)
        deallocate(gromol%mcontact, stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_dealloc
        allocate(gromol%mcontact(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        if (old_size > var_size) then
          gromol%mcontact(1:var_size) = mcontact(1:var_size)
        else
          gromol%mcontact(1:old_size) = mcontact(1:old_size)
        end if
        deallocate(mcontact, stat = alloc_stat)
      else
        old_size = 0
        allocate(gromol%mcontact(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
      end if

    case(GMGroMolMorphBB)
      if (allocated(gromol%morph_bb)) then

        old_size = size(gromol%morph_bb)
        if (old_size == var_size) &
          return
        allocate(morph_bb(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        morph_bb(:) = gromol%morph_bb(:)
        deallocate(gromol%morph_bb, stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_dealloc
        allocate(gromol%morph_bb(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        if (old_size > var_size) then
          gromol%morph_bb(1:var_size) = morph_bb(1:var_size)
        else
          gromol%morph_bb(1:old_size) = morph_bb(1:old_size)
        end if
        deallocate(morph_bb, stat = alloc_stat)
      else
        old_size = 0
        allocate(gromol%morph_bb(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
      end if

    case(GMGroMolMorphSC)
      if (allocated(gromol%morph_sc)) then

        old_size = size(gromol%morph_sc)
        if (old_size == var_size) &
          return
        allocate(morph_sc(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        morph_sc(:) = gromol%morph_sc(:)
        deallocate(gromol%morph_sc, stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_dealloc
        allocate(gromol%morph_sc(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
        if (old_size > var_size) then
          gromol%morph_sc(1:var_size) = morph_sc(1:var_size)
        else
          gromol%morph_sc(1:old_size) = morph_sc(1:old_size)
        end if
        deallocate(morph_sc, stat = alloc_stat)
      else
        old_size = 0
        allocate(gromol%morph_sc(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc
      end if

    case(GMGroMolPWMcos)
      if (allocated(gromol%pwmcos)) then

        old_size = size(gromol%pwmcos)
        if (old_size == var_size) &
            return

        allocate(pwmcos(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) call error_msg_alloc
        pwmcos(:) = gromol%pwmcos(:)

        deallocate(gromol%pwmcos, stat = alloc_stat)
        if (alloc_stat /= 0) call error_msg_dealloc

        allocate(gromol%pwmcos(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) call error_msg_alloc

        if (old_size > var_size) then
          gromol%pwmcos(1:var_size) = pwmcos(1:var_size)
        else
          gromol%pwmcos(1:old_size) = pwmcos(1:old_size)
        end if

        deallocate(pwmcos, stat = alloc_stat)

      else

        old_size = 0
        allocate(gromol%pwmcos(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) call error_msg_alloc

      end if

    case(GMGroMolPWMcosns)
      if (allocated(gromol%pwmcosns)) then

        old_size = size(gromol%pwmcosns)
        if (old_size == var_size) &
            return

        allocate(pwmcosns(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) call error_msg_alloc
        pwmcosns(:) = gromol%pwmcosns(:)

        deallocate(gromol%pwmcosns, stat = alloc_stat)
        if (alloc_stat /= 0) call error_msg_dealloc

        allocate(gromol%pwmcosns(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) call error_msg_alloc

        if (old_size > var_size) then
          gromol%pwmcosns(1:var_size) = pwmcosns(1:var_size)
        else
          gromol%pwmcosns(1:old_size) = pwmcosns(1:old_size)
        end if

        deallocate(pwmcosns, stat = alloc_stat)

      else

        old_size = 0
        allocate(gromol%pwmcosns(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) call error_msg_alloc

      end if

    case(GMGroMolIDRHPS)
      if (allocated(gromol%idr_hps)) then

        old_size = size(gromol%idr_hps)
        if (old_size == var_size) &
              return
        allocate(idr_hps(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
        idr_hps(:) = gromol%idr_hps(:)
        deallocate(gromol%idr_hps, stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_dealloc
        allocate(gromol%idr_hps(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
        if (old_size > var_size) then
          gromol%idr_hps(1:var_size) = idr_hps(1:var_size)
        else
          gromol%idr_hps(1:old_size) = idr_hps(1:old_size)
        end if
        deallocate(idr_hps, stat = alloc_stat)
      else
        old_size = 0
        allocate(gromol%idr_hps(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
              call error_msg_alloc
      end if

    case(GMGroMolIDRKH)
      if (allocated(gromol%idr_kh)) then

        old_size = size(gromol%idr_kh)
        if (old_size == var_size) &
            return
        allocate(idr_kh(old_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
            call error_msg_alloc
        idr_kh(:) = gromol%idr_kh(:)
        deallocate(gromol%idr_kh, stat = alloc_stat)
        if (alloc_stat /= 0) &
            call error_msg_dealloc
        allocate(gromol%idr_kh(var_size), stat = alloc_stat)
        if (alloc_stat /= 0) &
            call error_msg_alloc
        if (old_size > var_size) then
          gromol%idr_kh(1:var_size) = idr_kh(1:var_size)
        else
          gromol%idr_kh(1:old_size) = idr_kh(1:old_size)
        end if
        deallocate(idr_kh, stat = alloc_stat)
      else
        old_size = 0
        allocate(gromol%idr_kh(var_size), stat = alloc_stat)
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

    if (allocated(grotop%flangltypes)) then
      deallocate(grotop%flangltypes, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(grotop%dihetypes)) then
      deallocate(grotop%dihetypes, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(grotop%fldihetypes)) then
      deallocate(grotop%fldihetypes, stat = dealloc_stat)
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

    if (allocated(grotop%ez_membrane)) then
      deallocate(grotop%ez_membrane, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(grotop%basestacktypes)) then
      deallocate(grotop%basestacktypes, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
            call error_msg_dealloc
    end if

    if (allocated(grotop%basepairtypes)) then
      deallocate(grotop%basepairtypes, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
            call error_msg_dealloc
    end if

    if (allocated(grotop%basecrosstypes)) then
      deallocate(grotop%basecrosstypes, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
            call error_msg_dealloc
    end if

    if (allocated(grotop%cgdnaexvtypes)) then
      deallocate(grotop%cgdnaexvtypes, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
          call error_msg_dealloc
    end if

    if (allocated(grotop%cg_ele_mol_pairs)) then
      deallocate(grotop%cg_ele_mol_pairs, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
          call error_msg_dealloc
    end if

    if (allocated(grotop%cg_kh_mol_pairs)) then
      deallocate(grotop%cg_kh_mol_pairs, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
          call error_msg_dealloc
    end if

    if (allocated(grotop%pwmcos_mol_pairs)) then
      deallocate(grotop%pwmcos_mol_pairs, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
          call error_msg_dealloc
    end if

    if (allocated(grotop%pwmcosns_mol_pairs)) then
      deallocate(grotop%pwmcosns_mol_pairs, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
          call error_msg_dealloc
    end if

    if (allocated(grotop%cg_IDR_HPS_atomtypes)) then
      deallocate(grotop%cg_IDR_HPS_atomtypes, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
          call error_msg_dealloc
    end if

    if (allocated(grotop%cg_KH_atomtypes)) then
      deallocate(grotop%cg_KH_atomtypes, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
          call error_msg_dealloc
    end if

    if (allocated(grotop%cg_pair_MJ_eps)) then
      deallocate(grotop%cg_pair_MJ_eps, stat = dealloc_stat)
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

    if (allocated(gromol%posress)) then
      deallocate(gromol%posress, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(gromol%mcontact)) then
      deallocate(gromol%mcontact, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(gromol%morph_bb)) then
      deallocate(gromol%morph_bb, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(gromol%morph_sc)) then
      deallocate(gromol%morph_sc, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(gromol%pwmcos)) then
      deallocate(gromol%pwmcos, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
          call error_msg_dealloc
    end if

    if (allocated(gromol%pwmcosns)) then
      deallocate(gromol%pwmcosns, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
          call error_msg_dealloc
    end if

    if (allocated(gromol%idr_hps)) then
      deallocate(gromol%idr_hps, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
          call error_msg_dealloc
    end if

    if (allocated(gromol%idr_kh)) then
      deallocate(gromol%idr_kh, stat = dealloc_stat)
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
    character(100)              :: err
    logical                     :: unk, omit

    type(s_grotop_mol), pointer :: gromol


    directive = DUnknown
    unk = .false.

    do while(.true.)
      if (.not. gro_pp_next_line(file, line, err)) &
        goto 100

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

      case (DFlexible_Local_Angle)
        call read_flangle_types(file, grotop)

      case (DDihedralTypes)
        call read_dihedral_types(file, grotop)

      case (DFlexible_Local_Dihedral_Angle)
        call read_fldihe_types(file, grotop)

      case (DConstrTypes)
        call read_constraint_types(file, grotop)

      case (DCmapTypes)
        call read_cmap_types(file, grotop)

      case (DBaseStackTypes)
        call read_base_stack_types(file, grotop)

      case (DBasePairTypes)
        call read_base_pair_types(file, grotop)

      case (DBaseCrossTypes)
        call read_base_cross_types(file, grotop)

      case (DCGDNAExvTypes)
        call read_base_exv_types(file, grotop)

      case (DCGEleMolPairs)
        call read_cg_ele_mol_pairs(file, grotop)

      case (DCGKHMolPairs)
        call read_cg_kh_mol_pairs(file, grotop)

      case (DCGPWMcosMolPairs)
        call read_pwmcos_mol_pairs(file, grotop)

      case (DCGPWMcosnsMolPairs)
        call read_pwmcosns_mol_pairs(file, grotop)

      case (DCGPWMcos)
        call read_PWMcos(file, grotop, gromol)

      case (DCGPWMcosns)
        call read_PWMcosns(file, grotop, gromol)

      case (DCGIDRHPSAtomTypes)
        call read_cg_IDR_HPS_atom_types(file, grotop)

      case (DCGIDRHPSRegion)
        call read_idr_hps(file, grotop, gromol)

      case (DCGKHAtomTypes)
        call read_cg_KH_atom_types(file, grotop)

      case (DCGIDRKHRegion)
        call read_idr_kh(file, grotop, gromol)

      case (DCGPairMJEpsilon)
        call read_cg_pair_MJ_eps(file, grotop)

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
        call read_exclusions(file, grotop, gromol)

      case (DConstraints)
        call read_constraints(file, grotop, gromol)

      case (DPairs)
        call read_pairs(file, grotop, gromol)

      case (DSettles)
        call read_settles(file, grotop, gromol)

      case (DVirtualSites2)
        call read_virtual_sites2(file, gromol)

      case (DVirtualSites3)
        call read_virtual_sites3(file, gromol)

      case (DVirtualSites4)
        call read_virtual_sites4(file, gromol)

      case (DVirtualSitesN)
        call read_virtual_sitesn(file, gromol)

      case (DPosition_Restraints)
        call read_position_restraints(file, grotop, gromol)

      case (DSystem)
        call read_system()

      case (DMolecules)
        call read_molecules(file, grotop)

      case (DGomodel)
        call read_gomodel(file, grotop)

      case (DCg_Atoms)
        call read_cgatoms(file, grotop, gromol)

      case (DMulticontact)
        call read_multicontact(file, grotop, gromol)

      case (DMorphBB)
        call read_morph_bb(file, grotop, gromol)

      case (DMorphSC)
        call read_morph_sc(file, grotop, gromol)

      case (DEzMembrane)
        call read_membrane_type(file, grotop)

      case (DUnknown)
        if (main_rank) &
          write(MsgOut,*) 'Read_Grotop> INFO. Unknown directive:'// &
               line(1:index(line,']',.true.))
        unk = .true.

      end select
    end do

100 continue
    if (err /= 'end of file.') &
      call error_msg('Read_gro_top> '//err)

    grotop%num_atomtypes      = size_grotop(grotop, GroTopAtomType)
    grotop%num_bondtypes      = size_grotop(grotop, GroTopBondType)
    grotop%num_angltypes      = size_grotop(grotop, GroTopAnglType)
    grotop%num_flangltypes    = size_grotop(grotop, GroTopFlAnglType)
    grotop%num_dihetypes      = size_grotop(grotop, GroTopDiheType)
    grotop%num_fldihetypes    = size_grotop(grotop, GroTopFlDiheType)
    grotop%num_constrtypes    = size_grotop(grotop, GroTopConstrType)
    grotop%num_cmaptypes      = size_grotop(grotop, GroTopCmapType)
    grotop%num_nbonparms      = size_grotop(grotop, GroTopNbonParm)
    grotop%num_moltypes       = size_grotop(grotop, GroTopMolType)
    grotop%num_molss          = size_grotop(grotop, GroTopMols)
    grotop%num_membranetypes  = size_grotop(grotop, GroTopEzMembrane)
    grotop%num_basestacktypes = size_grotop(grotop, GroTopBaseStackType)
    grotop%num_basepairtypes  = size_grotop(grotop, GroTopBasePairType)
    grotop%num_basecrosstypes = size_grotop(grotop, GroTopBaseCrossType)
    grotop%num_cgdnaexvtypes  = size_grotop(grotop, GroTopCGDNAExvType)
    grotop%num_cgelemolpairs  = size_grotop(grotop, GroTopCGEleMolPair)
    grotop%num_cgkhmolpairs   = size_grotop(grotop, GroTopCGKHMolPair)
    grotop%num_pwmcosmolpairs = size_grotop(grotop, GroTopPWMcosMolPair)
    grotop%num_pwmcosnsmolpairs = size_grotop(grotop, GroTopPWMcosnsMolPair)
    grotop%num_cg_IDR_HPS_atomtypes = size_grotop(grotop, GroTopCGIDRHPSAtomType)
    grotop%num_cg_KH_atomtypes= size_grotop(grotop, GroTopCGKHAtomType)
    grotop%num_cg_pair_MJ_eps = size_grotop(grotop, GroTopCGPairMJEpsilon)

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

      !shinobu-edited
      write(MsgOut,'(A)') 'Read_Grotop> Summary of Grotopfile'
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
          write(MsgOut,'(6X,A14,I10)') &
               'num_PWMcos  = ', grotop%moltypes(i)%mol%num_pwmcos
          write(MsgOut,'(6X,A14,I10)') &
               'num_PWMcosns= ', grotop%moltypes(i)%mol%num_pwmcosns
          write(MsgOut,'(6X,A14,I10)') &
               'num_IDR_HPS = ', grotop%moltypes(i)%mol%num_idr_hps
          write(MsgOut,'(6X,A14,I10)') &
               'num_IDR_KH  = ', grotop%moltypes(i)%mol%num_idr_kh
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
      if (grotop%num_flangltypes > 0) then
        write(MsgOut,'(A20,I10)') &
           '  num_flangltypes   = ', grotop%num_flangltypes
      endif
      if (grotop%num_fldihetypes > 0) then
        write(MsgOut,'(A20,I10)') &
           '  num_fldihetypes   = ', grotop%num_fldihetypes
      endif
      write(MsgOut,'(A20,I10,A20,I10)') &
           '  num_constrtypes = ', grotop%num_constrtypes, &
           '  num_cmaptypes   = ', grotop%num_cmaptypes
      write(MsgOut,'(A20,I10)') &
           '  num_nbonparms   = ', grotop%num_nbonparms
      if (grotop%num_basestacktypes > 0) then
        write(MsgOut,'(A24,I6)') &
            '  num_basestacktypes  = ', grotop%num_basestacktypes
      endif
      if (grotop%num_basepairtypes > 0) then
        write(MsgOut,'(A24,I6)') &
            '  num_basepairtypes   = ', grotop%num_basepairtypes
      endif
      if (grotop%num_basecrosstypes > 0) then
        write(MsgOut,'(A24,I6)') &
            '  num_basecrosstypes  = ', grotop%num_basecrosstypes
      endif
      if (grotop%num_cgdnaexvtypes > 0) then
        write(MsgOut,'(A24,I6)') &
            '  num_cgdnaexvtypes   = ', grotop%num_cgdnaexvtypes
      endif
      if (grotop%num_cgelemolpairs > 0) then
        write(MsgOut,'(A24,I6)') &
            '  num_cgelemolpairs   = ', grotop%num_cgelemolpairs
      endif
      if (grotop%num_cgkhmolpairs > 0) then
        write(MsgOut,'(A24,I6)') &
            '  num_cgkhmolpairs   = ', grotop%num_cgkhmolpairs
      endif
      if (grotop%num_pwmcosmolpairs > 0) then
        write(MsgOut,'(A24,I6)') &
            '  num_pwmcosmolpairs  = ', grotop%num_pwmcosmolpairs
      endif
      if (grotop%num_pwmcosnsmolpairs > 0) then
        write(MsgOut,'(A24,I6)') &
            '  num_pwmcosnsmolpairs  = ', grotop%num_pwmcosnsmolpairs
      endif
      if (grotop%num_cg_IDR_HPS_atomtypes > 0) then
        write(MsgOut,'(A30,I6)') &
            '  num_cg_IDR_HPS_atomtypes  = ', grotop%num_cg_IDR_HPS_atomtypes
      endif
      if (grotop%num_cg_KH_atomtypes > 0) then
        write(MsgOut,'(A24,I6)') &
            '  num_cg_KH_atomtypes  = ', grotop%num_cg_KH_atomtypes
      endif
      write(MsgOut,'(A)') ''
      if (grotop%num_cg_pair_MJ_eps > 0) then
        write(MsgOut,'(A24,I6)') &
            '  num_cg_pair_MJ_eps   = ', grotop%num_cg_pair_MJ_eps
      endif
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

    call gro_pp_next_line_s(file, line)
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
    character(100)           :: error

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

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

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
    character(100)           :: error

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

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

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
    character(100)           :: error

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

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

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
  !  Subroutine    read_flangle_types
  !> @brief        read section [flexible_local_angle] #original in GENESIS-Cafemol
  !! @authors      CK
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_flangle_types(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),  target, intent(inout) :: grotop

    ! local variables
    integer                  :: i, j, cnt, old_cnt, loop
    character(MaxLine)       :: line
    character(100)           :: error

    type(s_flangltype), pointer :: at


    ! check count
    !
    cnt = check_section_count(file)

    ! allocate memory
    !
    old_cnt = size_grotop(grotop, GroTopFlAnglType)
    call realloc_grotop(grotop, GroTopFlAnglType, old_cnt+cnt)

    if (VerboseOn .and. main_rank) &
      write(MsgOut,'(" Read_Grotop> [flexible_local_angle] :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt

    ! read data
    !
    do i = 1, cnt

      at => grotop%flangltypes(old_cnt+i)

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

      if (match_format(line, 'SNNN')) then

        read(line,*)      &
           at%atom_type,  &
           at%func,       &
           loop

        if (at%func /= 1)  then
          if (main_rank) then
            write(MsgOut,'(" Read_Grotop> Error: [flexible_local_angle] '// &
               'not supported function type. [", i3, "]")') at%func
          endif
          goto 900
        endif

        allocate(at%theta(loop))
        allocate(at%efunc(loop))
        allocate(at%d2func(loop))

        read(line,*)      &
           at%atom_type,  &
           at%func,       &
           loop,          &
           (at%theta(j), &
            at%efunc(j), &
            at%d2func(j), j = 1, loop)

     else

        goto 900

      end if

    end do

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [flexible_local_angle]')

  end subroutine read_flangle_types

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
    character(100)           :: error

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

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

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

      dt%wild_num = 0
      if (dt%atom_type1 == 'X') dt%wild_num = dt%wild_num + 1
      if (dt%atom_type2 == 'X') dt%wild_num = dt%wild_num + 1
      if (dt%atom_type4 == 'X') dt%wild_num = dt%wild_num + 1

    end do

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [dihedraltypes]')

  end subroutine read_dihedral_types

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_fldihe_types
  !> @brief        read section [flexible_local_dihedral_angle]
  !                #original in GENESIS-Cafemol
  !! @authors      JJ
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_fldihe_types(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),  target, intent(inout) :: grotop

    ! local variables
    integer                  :: i, j, cnt, old_cnt, loop
    character(MaxLine)       :: line
    character(100)           :: error

    type(s_fldihetype), pointer :: at


    ! check count
    !
    cnt = check_section_count(file)

    ! allocate memory
    !
    old_cnt = size_grotop(grotop, GroTopFlDiheType)
    call realloc_grotop(grotop, GroTopFlDiheType, old_cnt+cnt)

    if (VerboseOn .and. main_rank) &
      write(MsgOut,'(" Read_Grotop> [flexible_local_dihedral_angle] :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt

    ! read data
    !
    do i = 1, cnt

      at => grotop%fldihetypes(old_cnt+i)

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

      if (match_format(line, 'SSNNN')) then

        read(line,*)      &
           at%atom_type1, &
           at%atom_type2, &
           at%func,       &
           loop

        if (at%func /= 1)  then
          if (main_rank) then
            write(MsgOut,'(" Read_Grotop> Error: [flexible_local_dihedral_angle] '// &
               'not supported function type. [", i3, "]")') at%func
          endif
          goto 900
        endif

        allocate(at%coef(loop))

        read(line,*)      &
           at%atom_type1, &
           at%atom_type2, &
           at%func,       &
           loop,          &
           (at%coef(j), j = 1, loop)

     else

        goto 900

      end if

    end do

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [flexible_local_dihedral_angle]')

  end subroutine read_fldihe_types

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
    character(100)           :: error

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

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

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
    character(100)           :: error

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

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

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
  !  Subroutine    read_base_stack_types
  !> @brief        read section [ basestacktypes ]
  !! @authors      CT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_base_stack_types(file, grotop)

    ! formal arguments
    integer,                       intent(in)             :: file
    type(s_grotop),                target, intent(inout)  :: grotop

    ! local variables
    integer                                               :: i, cnt, old_cnt
    character(MaxLine)                                    :: line
    character(100)                                        :: error

    type(s_basestacktype),         pointer                :: bs


    ! check count
    cnt = check_section_count(file)

    ! allocate memory
    old_cnt = size_grotop(grotop, GroTopBaseStackType)

    call realloc_grotop(grotop, GroTopBaseStackType, old_cnt+cnt)

    if (VerboseOn .and. main_rank) then
      write(MsgOut,'(" Read_Grotop> [ basestacktypes ] proper   :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt
    end if

    ! read data
    do i = 1, cnt
      bs => grotop%basestacktypes(old_cnt+i)
      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900
      if (match_format(line, 'SS1NNN')) then
        read(line,*)         &
          bs%base_type5, &
          bs%base_type3, &
          bs%func,       &
          bs%epsilon,    &
          bs%sigma,      &
          bs%theta_bs
      else
        call error_msg_grotop(file, 'Read_Grotop> read error. [ basestacktypes ]')
      end if
    end do

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [ basestacktypes ]')

  end subroutine read_base_stack_types

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_base_pair_types
  !> @brief        read section [ basepairtypes ]
  !! @authors      CT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_base_pair_types(file, grotop)

    ! formal arguments
    integer,                 intent(in)             :: file
    type(s_grotop),  target, intent(inout)  :: grotop

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line
    character(100)           :: error

    type(s_basepairtype), pointer :: bp


    ! check count
    cnt = check_section_count(file)

    ! allocate memory
    old_cnt = size_grotop(grotop, GroTopBasePairType)

    call realloc_grotop(grotop, GroTopBasePairType, old_cnt+cnt)

    if (VerboseOn .and. main_rank) then
      write(MsgOut,'(" Read_Grotop> [ basepairtypes ] proper   :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt
    end if

    ! read data
    !
    do i = 1, cnt

      bp => grotop%basepairtypes(old_cnt+i)

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

      if (match_format(line, 'SS1NNNNNN')) then
        read(line,*)        &
            bp%base_type_a, &
            bp%base_type_b, &
            bp%func,        &
            bp%theta_1,     &
            bp%theta_2,     &
            bp%phi_1,       &
            bp%theta_3,     &
            bp%sigma,       &
            bp%epsilon
      else
        call error_msg_grotop(file, 'Read_Grotop> read error. [ basepairtypes ]')
      end if
    end do

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [ basepairtypes ]')

  end subroutine read_base_pair_types

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_base_cross_types
  !> @brief        read section [ basecrosstypes ]
  !! @authors      CT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_base_cross_types(file, grotop)

    ! formal arguments
    integer,                 intent(in)             :: file
    type(s_grotop),  target, intent(inout)  :: grotop

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line
    character(100)           :: error

    type(s_basecrosstype), pointer :: cs


    ! check count
    cnt = check_section_count(file)

    ! allocate memory
    old_cnt = size_grotop(grotop, GroTopBaseCrossType)

    call realloc_grotop(grotop, GroTopBaseCrossType, old_cnt+cnt)

    if (VerboseOn .and. main_rank) then
      write(MsgOut,'(" Read_Grotop> [ basecrosstypes ] proper   :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt
    end if

    ! read data
    !
    do i = 1, cnt

      cs => grotop%basecrosstypes(old_cnt+i)

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

      if (match_format(line, 'SS1NNN') .or. &
          match_format(line, 'SS2NNN')) then
        read(line,*)        &
            cs%base_type_a, &
            cs%base_type_b, &
            cs%func,        &
            cs%epsilon,     &
            cs%sigma,       &
            cs%theta_cs
      else
        call error_msg_grotop(file, 'Read_Grotop> read error. [ basecrosstypes ]')
      end if
    end do

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [ basecrosstypes ]')

  end subroutine read_base_cross_types

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_base_exv_types
  !> @brief        read section [ cgdnaexvtypes ]
  !! @authors      CT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_base_exv_types(file, grotop)

    ! formal arguments
    integer,                 intent(in)     :: file
    type(s_grotop),  target, intent(inout)  :: grotop

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line
    character(100)           :: error

    type(s_cgdnaexvtype), pointer :: exv


    ! check count
    cnt = check_section_count(file)

    ! allocate memory
    old_cnt = size_grotop(grotop, GroTopCGDNAExvType)

    call realloc_grotop(grotop, GroTopCGDNAExvType, old_cnt+cnt)

    if (VerboseOn .and. main_rank) then
      write(MsgOut,'(" Read_Grotop> [ cgdnaexvtypes ] proper   :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt
    end if

    ! read data
    !
    do i = 1, cnt

      exv => grotop%cgdnaexvtypes(old_cnt+i)

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

      if (match_format(line, 'SNN')) then
        read(line,*)       &
            exv%base_type, &
            exv%func,      &
            exv%sigma
      else
        call error_msg_grotop(file, 'Read_Grotop> read error. [ cgdnaexvtypes ]')
      end if
    end do

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [ cgdnaexvtypes ]')

  end subroutine read_base_exv_types

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_cg_IDR_HPS_atomtypes
  !> @brief        read section [ cg_IDR_HPS_atomtypes ]
  !! @authors      CT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_cg_IDR_HPS_atom_types(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),  target, intent(inout) :: grotop

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line
    character(100)           :: error

    type(s_cg_IDR_HPS_atomtype), pointer :: hps


    ! check count
    cnt = check_section_count(file)

    ! allocate memory
    old_cnt = size_grotop(grotop, GroTopCGIDRHPSAtomType)

    call realloc_grotop(grotop, GroTopCGIDRHPSAtomType, old_cnt+cnt)

    if (VerboseOn .and. main_rank) then
      write(MsgOut,'(" Read_Grotop> [ cg_IDR_HPS_atomtypes ] proper   :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt
    end if

    ! read data
    !
    do i = 1, cnt

      hps => grotop%cg_IDR_HPS_atomtypes(old_cnt+i)

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

      if (match_format(line, 'SNNNN')) then
        read(line,*)       &
            hps%type_name, &
            hps%mass,      &
            hps%charge,    &
            hps%sigma,     &
            hps%lambda
      else
        call error_msg_grotop(file, 'Read_Grotop> read error. [ cg_IDR_HPS_atomtypes ]')
      end if
    end do

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [ cg_IDR_HPS_atomtypes ]')

  end subroutine read_cg_IDR_HPS_atom_types

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_cg_KH_atomtypes
  !> @brief        read section [ cg_KH_atomtypes ]
  !! @authors      CT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_cg_KH_atom_types(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),  target, intent(inout) :: grotop

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line
    character(100)           :: error

    type(s_cg_KH_atomtype), pointer :: kh


    ! check count
    cnt = check_section_count(file)

    ! allocate memory
    old_cnt = size_grotop(grotop, GroTopCGKHAtomType)

    call realloc_grotop(grotop, GroTopCGKHAtomType, old_cnt+cnt)

    if (VerboseOn .and. main_rank) then
      write(MsgOut,'(" Read_Grotop> [ cg_KH_atomtypes ] proper   :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt
    end if

    ! read data
    !
    do i = 1, cnt

      kh => grotop%cg_KH_atomtypes(old_cnt+i)

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

      if (match_format(line, 'SNNN')) then
        read(line,*)       &
            kh%type_name, &
            kh%mass,      &
            kh%charge,    &
            kh%sigma
      else
        call error_msg_grotop(file, 'Read_Grotop> read error. [ cg_KH_atomtypes ]')
      end if
    end do

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [ cg_KH_atomtypes ]')

  end subroutine read_cg_KH_atom_types


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_cg_pair_MJ_eps
  !> @brief        read section [ cg_pair_energy_MJ_96 ]
  !! @authors      CT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_cg_pair_MJ_eps(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),  target, intent(inout) :: grotop

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line
    character(100)           :: error

    type(s_cg_pair_MJ_eps), pointer :: mj


    ! check count
    cnt = check_section_count(file)

    ! allocate memory
    old_cnt = size_grotop(grotop, GroTopCGPairMJEpsilon)

    call realloc_grotop(grotop, GroTopCGPairMJEpsilon, old_cnt+cnt)

    if (VerboseOn .and. main_rank) then
      write(MsgOut,'(" Read_Grotop> [ cg_pair_energy_MJ_96 ] proper   :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt
    end if

    ! read data
    !
    do i = 1, cnt

      mj => grotop%cg_pair_MJ_eps(old_cnt+i)

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

      if (match_format(line, 'SSN')) then
        read(line,*)        &
            mj%type_name_1, &
            mj%type_name_2, &
            mj%epsilon
      else
        call error_msg_grotop(file, 'Read_Grotop> read error. [ cg_pair_energy_MJ_96 ]')
      end if
    end do

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [ cg_pair_energy_MJ_96 ]')

  end subroutine read_cg_pair_MJ_eps


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_cg_ele_mol_pairs
  !> @brief        read section [ cg_ele_mol_pairs ]
  !! @authors      CT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_cg_ele_mol_pairs(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),  target, intent(inout) :: grotop

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line
    character(100)           :: error
    character                :: chr_tmp_hyphen
    character                :: chr_tmp_intermol
    character(3)             :: chr_tmp_switch

    type(s_cg_ele_mol_pair), pointer :: elepair


    ! check count
    cnt = check_section_count(file)

    ! allocate memory
    old_cnt = size_grotop(grotop, GroTopCGEleMolPair)

    call realloc_grotop(grotop, GroTopCGEleMolPair, old_cnt+cnt)

    if (VerboseOn .and. main_rank) then
      write(MsgOut,'(" Read_Grotop> [ cg_ele_chain_pairs ] proper   :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt
    end if

    ! read data
    !
    do i = 1, cnt

      elepair => grotop%cg_ele_mol_pairs(old_cnt+i)

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

      if (match_format(line, 'SN-NSN-N')) then
        read(line,*)            &
            chr_tmp_switch,     &
            elepair%grp1_start, &
            chr_tmp_hyphen,     &
            elepair%grp1_end,   &
            chr_tmp_intermol,   &
            elepair%grp2_start, &
            chr_tmp_hyphen,     &
            elepair%grp2_end
        if ( chr_tmp_intermol /= ':' ) then
          call error_msg_grotop(file, 'Read_Grotop> interaction format error. [ cg_ele_chain_pairs ]')
        end if
        if ( chr_tmp_switch == 'ON' ) then
          elepair%func = 1
        else if ( chr_tmp_switch == 'OFF' ) then
          elepair%func = 0
        else
          call error_msg_grotop(file, 'Read_Grotop> switch format error. [ cg_ele_chain_pairs ]')
        end if
      else if (match_format(line, 'SN-N')) then
        elepair%is_intermol = .false.
        read(line,*)            &
            chr_tmp_switch,     &
            elepair%grp1_start, &
            chr_tmp_hyphen,     &
            elepair%grp1_end
        if ( chr_tmp_switch == 'ON' ) then
          elepair%func = 1
        else if ( chr_tmp_switch == 'OFF' ) then
          elepair%func = 0
        else
          call error_msg_grotop(file, 'Read_Grotop> switch format error. [ cg_ele_chain_pairs ]')
        end if
      else
        call error_msg_grotop(file, 'Read_Grotop> format error. [ cg_ele_chain_pairs ]')
      end if
    end do

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [ cg_ele_chain_pairs ]')

  end subroutine read_cg_ele_mol_pairs


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_cg_kh_mol_pairs
  !> @brief        read section [ cg_kh_mol_pairs ]
  !! @authors      CT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_cg_kh_mol_pairs(file, grotop)

    ! formal arguments
    integer,                intent(in)             :: file
    type(s_grotop),         target, intent(inout)  :: grotop

    ! local variables
    integer                 :: i, cnt, old_cnt
    character(MaxLine)      :: line
    character(100)          :: error
    character               :: chr_tmp_hyphen
    character               :: chr_tmp_intermol
    character(3)            :: chr_tmp_switch

    type(s_cg_kh_mol_pair), pointer :: khpair


    ! check count
    cnt = check_section_count(file)

    ! allocate memory
    old_cnt = size_grotop(grotop, GroTopCGKHMolPair)

    call realloc_grotop(grotop, GroTopCGKHMolPair, old_cnt+cnt)

    if (VerboseOn .and. main_rank) then
      write(MsgOut,'(" Read_Grotop> [ cg_kh_chain_pairs ] proper   :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt
    end if

    ! read data
    !
    do i = 1, cnt

      khpair => grotop%cg_kh_mol_pairs(old_cnt+i)

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

      if (match_format(line, 'SN-NSN-N')) then
        read(line,*)           &
            chr_tmp_switch,    &
            khpair%grp1_start, &
            chr_tmp_hyphen,    &
            khpair%grp1_end,   &
            chr_tmp_intermol,  &
            khpair%grp2_start, &
            chr_tmp_hyphen,    &
            khpair%grp2_end
        if ( chr_tmp_intermol /= ':' ) then
          call error_msg_grotop(file, 'Read_Grotop> interaction format error. [ cg_kh_chain_pairs ]')
        end if
        if ( chr_tmp_switch == 'A' ) then
          khpair%func = 1
        else if ( chr_tmp_switch == 'B' ) then
          khpair%func = 2
        else if ( chr_tmp_switch == 'C' ) then
          khpair%func = 3
        else if ( chr_tmp_switch == 'D' ) then
          khpair%func = 4
        else if ( chr_tmp_switch == 'E' ) then
          khpair%func = 5
        else if ( chr_tmp_switch == 'F' ) then
          khpair%func = 6
        else if ( chr_tmp_switch == 'OFF' ) then
          khpair%func = 0
        else
          call error_msg_grotop(file, 'Read_Grotop> switch format error. [ cg_kh_chain_pairs ]')
        end if
      else if (match_format(line, 'SN-N')) then
        khpair%is_intermol = .false.
        read(line,*)           &
            chr_tmp_switch,    &
            khpair%grp1_start, &
            chr_tmp_hyphen,    &
            khpair%grp1_end
        if ( chr_tmp_switch == 'A' ) then
          khpair%func = 1
        else if ( chr_tmp_switch == 'B' ) then
          khpair%func = 2
        else if ( chr_tmp_switch == 'C' ) then
          khpair%func = 3
        else if ( chr_tmp_switch == 'D' ) then
          khpair%func = 4
        else if ( chr_tmp_switch == 'E' ) then
          khpair%func = 5
        else if ( chr_tmp_switch == 'F' ) then
          khpair%func = 6
        else if ( chr_tmp_switch == 'OFF' ) then
          khpair%func = 0
        else
          call error_msg_grotop(file, 'Read_Grotop> switch format error. [ cg_kh_chain_pairs ]')
        end if
      else
        call error_msg_grotop(file, 'Read_Grotop> format error. [ cg_kh_chain_pairs ]')
      end if
    end do

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [ cg_kh_chain_pairs ]')

  end subroutine read_cg_kh_mol_pairs


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_pwmcos_mol_pairs
  !> @brief        read section [ pwmcos_mol_pairs ]
  !! @authors      CT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_pwmcos_mol_pairs(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),  target, intent(inout) :: grotop

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line
    character(100)           :: error
    character                :: chr_tmp_hyphen
    character                :: chr_tmp_intermol
    character(3)             :: chr_tmp_switch

    type(s_pwmcos_mol_pair), pointer :: pwmcospair


    ! check count
    cnt = check_section_count(file)

    ! allocate memory
    old_cnt = size_grotop(grotop, GroTopPWMcosMolPair)

    call realloc_grotop(grotop, GroTopPWMcosMolPair, old_cnt+cnt)

    if (VerboseOn .and. main_rank) then
      write(MsgOut,'(" Read_Grotop> [ pwmcos_chain_pairs ] proper   :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt
    end if

    ! read data
    !
    do i = 1, cnt

      pwmcospair => grotop%pwmcos_mol_pairs(old_cnt+i)

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

      if (match_format(line, 'SN-NSN-N')) then
        read(line,*)            &
            chr_tmp_switch,     &
            pwmcospair%grp1_start, &
            chr_tmp_hyphen,     &
            pwmcospair%grp1_end,   &
            chr_tmp_intermol,   &
            pwmcospair%grp2_start, &
            chr_tmp_hyphen,     &
            pwmcospair%grp2_end
        if ( chr_tmp_intermol /= ':' ) then
          call error_msg_grotop(file, 'Read_Grotop> interaction format error. [ pwmcos_chain_pairs ]')
        end if
        if ( chr_tmp_switch == 'ON' ) then
          pwmcospair%func = 1
        else if ( chr_tmp_switch == 'OFF' ) then
          pwmcospair%func = 0
        else
          call error_msg_grotop(file, 'Read_Grotop> switch format error. [ pwmcos_chain_pairs ]')
        end if
      else
        call error_msg_grotop(file, 'Read_Grotop> format error. [ pwmcos_chain_pairs ]')
      end if
    end do

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [ pwmcos_chain_pairs ]')

  end subroutine read_pwmcos_mol_pairs


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_pwmcosns_mol_pairs
  !> @brief        read section [ pwmcosns_mol_pairs ]
  !! @authors      CT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_pwmcosns_mol_pairs(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),  target, intent(inout) :: grotop

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line
    character(100)           :: error
    character                :: chr_tmp_hyphen
    character                :: chr_tmp_intermol
    character(3)             :: chr_tmp_switch

    type(s_pwmcosns_mol_pair), pointer :: pwmcosnspair


    ! check count
    cnt = check_section_count(file)

    ! allocate memory
    old_cnt = size_grotop(grotop, GroTopPWMcosnsMolPair)

    call realloc_grotop(grotop, GroTopPWMcosnsMolPair, old_cnt+cnt)

    if (VerboseOn .and. main_rank) then
      write(MsgOut,'(" Read_Grotop> [ pwmcosns_chain_pairs ] proper   :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt
    end if

    ! read data
    !
    do i = 1, cnt

      pwmcosnspair => grotop%pwmcosns_mol_pairs(old_cnt+i)

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

      if (match_format(line, 'SN-NSN-N')) then
        read(line,*)            &
            chr_tmp_switch,     &
            pwmcosnspair%grp1_start, &
            chr_tmp_hyphen,     &
            pwmcosnspair%grp1_end,   &
            chr_tmp_intermol,   &
            pwmcosnspair%grp2_start, &
            chr_tmp_hyphen,     &
            pwmcosnspair%grp2_end
        if ( chr_tmp_intermol /= ':' ) then
          call error_msg_grotop(file, 'Read_Grotop> interaction format error. [ pwmcosns_chain_pairs ]')
        end if
        if ( chr_tmp_switch == 'ON' ) then
          pwmcosnspair%func = 1
        else if ( chr_tmp_switch == 'OFF' ) then
          pwmcosnspair%func = 0
        else
          call error_msg_grotop(file, 'Read_Grotop> switch format error. [ pwmcosns_chain_pairs ]')
        end if
      else
        call error_msg_grotop(file, 'Read_Grotop> format error. [ pwmcosns_chain_pairs ]')
      end if
    end do

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [ pwmcosns_chain_pairs ]')

  end subroutine read_pwmcosns_mol_pairs


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
    character(100)           :: error

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

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

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
    character(MaxLine)          :: line
    character(100)              :: error
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
    if (.not. gro_pp_next_line(file, line, error)) &
      goto 900
    read(line,*) mt%name, mt%exclude_nbon

    if (VerboseOn .and. main_rank) &
      write(MsgOut,'(" Read_Grotop> [moleculetype] "'// &
                   ',i5,": ",a20)') old_cnt+1, mt%name
    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [moleculetype]')

  end subroutine read_molecule_type

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_membrane_type
  !> @brief        read section [ez_membrane]
  !! @authors      CK
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_membrane_type(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop), target,  intent(inout) :: grotop

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line
    character(100)           :: error
    character(20)            :: read_name

    type(s_ezmembrane),    pointer :: ez


    ! check count
    !
    cnt = check_section_count(file)

    ! allocate memory
    !
    old_cnt = size_grotop(grotop, GroTopEzMembrane)
    call realloc_grotop(grotop, GroTopEzMembrane, old_cnt+cnt)

    if (VerboseOn .and. main_rank) &
      write(MsgOut,'(" Read_Grotop>   [ez_membrane] :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt

    ! read data
    !
    do i = 1, cnt

      ez => grotop%ez_membrane(old_cnt+i)

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

      if (match_format(line, 'SNNNN')) then

        read(line,*)       &
             ez%atom_type, &
             ez%func,      &
             ez%de,        &
             ez%zmin,      &
             ez%polym

      else

        goto 900

      end if

    end do

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [ez_membrane]')

  end subroutine read_membrane_type

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
    character(100)           :: error
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

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

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

      else if (match_format(line, 'NSNSSN')) then

        read(line,*)          &
             at%atom_no,      &
             read_name,       &
             at%residue_no,   &
             at%residue_name, &
             at%atom_name,    &
             at%charge_group

        call rename_atom_type_name(grotop, read_name, at%atom_type)

        if (.not. search_atom_mass(grotop%atomtypes, at%atom_type, at%mass)) &
          call error_msg_grotop(file, 'Read_Grotop> [atoms] '// &
               'atom type not found: '//at%atom_type)

        if (.not.search_atom_charge(grotop%atomtypes, at%atom_type, at%charge)) &
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
    character(100)           :: error
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

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

      if (match_format(line, 'NN1NN') .or. &
          match_format(line, 'NN2NN') .or. &
          ! ~CG~ 3SPN.2C DNA: read in param also for quartic bonds
          match_format(line, 'NNNNN')) then

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
    character(100)           :: error
    character(6)             :: at1, at2, at3
    !shinobu-edited
    character(6)             :: fnerr

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

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

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

      else if (match_format(line, 'NNNNNNNNNN')) then

        read(line,*)       &
             an%atom_idx1, &
             an%atom_idx2, &
             an%atom_idx3, &
             an%func,      &
             an%theta_0,   &
             an%theta_01,  &
             an%kt,        &
             an%kt1,       &
             an%cgamma,    &
             an%epsa

        if (an%func /= 10) then
          write(fnerr,'(i3)') an%func
          call error_msg_grotop(file, 'Read_Grotop> [angles] '// &
               'not supported function type: '//fnerr)
        endif

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

      ! shinobu-edited
      else if (match_format(line, 'NNNNNNN')) then

        ! AICG2+ Gaussian (local) angle
        read(line,*)       &
             an%atom_idx1, &
             an%atom_idx2, &
             an%atom_idx3, &
             an%func

        if (an%func /= 21) then
          write(fnerr,'(i3)') an%func
          call error_msg_grotop(file, 'Read_Grotop> [angles] '// &
               'not supported function type: '//fnerr)
        endif

      read(line,*)       &
             an%atom_idx1, &
             an%atom_idx2, &
             an%atom_idx3, &
             an%func,      &
             an%theta_0,   &
             an%kt,        &
             an%w

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

      !shinobu-edited
      else if (match_format(line, 'NNNN')) then
        ! AICG2+ Flexible angle
        read(line,*)       &
             an%atom_idx1, &
             an%atom_idx2, &
             an%atom_idx3, &
             an%func

        if (an%func /= 22) then
          write(fnerr,'(i3)') an%func
          call error_msg_grotop(file, 'Read_Grotop> [angles] '// &
               'not supported function type: '//fnerr)
        endif

        at1 = gromol%atoms(an%atom_idx1)%atom_type
        at2 = gromol%atoms(an%atom_idx2)%atom_type
        at3 = gromol%atoms(an%atom_idx3)%atom_type

        if (.not.search_flangl_param(grotop%flangltypes,     &
             at1, at2, at3, an%func, an%types)) &
          call error_msg_grotop(file, 'Read_Grotop> [angles] '// &
               'flexible angle atom type not found: '//at2)


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
    character(100)           :: error
    character(6)             :: at1, at2, at3, at4
    !shinobu-edited
    character(6)             :: fnerr

    type(s_dihe),    pointer :: dh

    ! check count
    !
    cnt = check_section_count(file)

    cnt_tri      = cnt*MAX_MULTIPLY_NUM
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
      if (ii+MAX_MULTIPLY_NUM > cur_size) then
        call realloc_grotop_mol(gromol, GMGroMolDihe, cur_size+cnt_tri)
        cur_size = size_grotop_mol(gromol, GMGroMolDihe)
      end if

      dh => gromol%dihes(old_dihedcnt+ii)

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

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

      else if (match_format(line, 'NNNNNNNN')) then

        ! AICG2+ Gaussian (local) dihedral
        read(line,*) &
            dh%atom_idx1, &
            dh%atom_idx2, &
            dh%atom_idx3, &
            dh%atom_idx4, &
            dh%func

        if (dh%func /= 21 &
            .and. dh%func /= 31 &
            .and. dh%func /= 32 &
            .and. dh%func /= 33 &
            .and. dh%func /= 41 &
            .and. dh%func /= 43 &
            .and. dh%func /= 52 &
            ) then
          write(fnerr,'(i3)') dh%func
          call error_msg_grotop(file, 'Read_Grotop> [dihedrals] '// &
              'not supported function type: '//fnerr)
        end if

        if (dh%func == 21 .or. dh%func == 41 .or. dh%func == 43) then
          read(line,*) &
              dh%atom_idx1, &
              dh%atom_idx2, &
              dh%atom_idx3, &
              dh%atom_idx4, &
              dh%func,      &
              dh%theta_0,   &
              dh%kp,        &
              dh%w
        end if

        if (dh%func == 31 .or. dh%func == 32 .or. dh%func == 33) then
          read(line,*) &
              dh%atom_idx1, &
              dh%atom_idx2, &
              dh%atom_idx3, &
              dh%atom_idx4, &
              dh%func, &
              dh%ps,   &
              dh%kp,   &
              dh%multiplicity
          end if


      else if (match_format(line, 'NNNNN')) then
        ! AICG2+ Flexible dihedral angle
        read(line,*)       &
             dh%atom_idx1, &
             dh%atom_idx2, &
             dh%atom_idx3, &
             dh%atom_idx4, &
             dh%func

        if (dh%func /= 22 .and. dh%func /= 52) then
          write(fnerr,'(i3)') dh%func
          call error_msg_grotop(file, 'Read_Grotop> [dihedrals] '// &
               'not supported function type: '//fnerr)
        endif

        at1 = gromol%atoms(dh%atom_idx1)%atom_type
        at2 = gromol%atoms(dh%atom_idx2)%atom_type
        at3 = gromol%atoms(dh%atom_idx3)%atom_type
        at4 = gromol%atoms(dh%atom_idx4)%atom_type

        if (.not.search_fldihe_param(grotop%fldihetypes,     &
             at1, at2, at3, at4, dh%func, dh%types)) &
          call error_msg_grotop(file, 'Read_Grotop> [dihedral] '// &
               'flexible dihedral angle atom types not found: '//at2// &
               'and '//at3)

      else if (match_format(line, 'NNNNN')) then

        ! unknown
        read(line,*) &
             dh%atom_idx1, &
             dh%atom_idx2, &
             dh%atom_idx3, &
             dh%atom_idx4, &
             dh%func

       if (dh%func == 22) then
         write(MsgOut,'(" Read_Grotop> Underconstructed ")')

       else

         write(MsgOut,'(" Read_Grotop> WARNING: [dihedrals] '// &
              'not supported function type. [", i3, "]")') dh%func

       endif

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
    character(100)           :: error
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

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

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
  !! @param[inout] grotop : GROMACS TOP information
  !! @param[inout] gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_exclusions(file, grotop, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(inout) :: grotop
    type(s_grotop_mol),      pointer       :: gromol

    ! local variables
    integer                  :: i, j, cnt, cnt2, cnt3, old_cnt
!!!original
!!!!!!    character(1000)          :: line
!!!original
!!!develop
    character(MaxLine)       :: line
!!!develop

    integer,     allocatable :: v(:)


    ! check count
    !

    cnt = check_section_count(file)

    call gro_pp_record_pos(file)

    cnt2 = 0
    do i = 1, cnt
      call gro_pp_next_line_s(file, line)
      cnt2 = cnt2 + split_num(line) - 1
    end do

    call gro_pp_restore_pos(file)

    ! allocate memory
    !

    old_cnt = size_grotop_mol(gromol, GMGroMolExcl)
    call realloc_grotop_mol(gromol, GMGroMolExcl, old_cnt+cnt2)

    if (VerboseOn .and. main_rank) &
      write(MsgOut,'(" Read_Grotop>   [exclusions] :"'// &
           ',i5," (total:",i5,")")') cnt2, old_cnt+cnt2

    ! read data
    !

    cnt2 = old_cnt

    do i = 1, cnt

      call gro_pp_next_line_s(file, line)
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
    character(100)           :: error
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

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

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
  !! @param[inout] grotop : GROMACS TOP information
  !! @param[inout] gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_pairs(file, grotop, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(inout) :: grotop
    type(s_grotop_mol),      pointer       :: gromol

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line
    character(100)           :: error
    character(6)             :: at1, at2

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

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

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

!shinobu-edited
        read(line,*) &
             pd%atom_idx1, &
             pd%atom_idx2, &
             pd%func

        if(pd%func==21) then
             read(line,*) &
                pd%atom_idx1, &
                pd%atom_idx2, &
                pd%func,      &
                pd%r0,        &
                pd%khh
        else
             read(line,*) &
                pd%atom_idx1, &
                pd%atom_idx2, &
                pd%func,      &
                pd%v,         &
                pd%w
        endif

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
  !! @param[inout] grotop : GROMACS TOP information
  !! @param[inout] gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_settles(file, grotop, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(inout) :: grotop
    type(s_grotop_mol),      pointer       :: gromol

    ! local variables
    character(MaxLine)       :: line
    character(100)           :: error


    if (.not. gro_pp_next_line(file, line, error)) &
      goto 900

    read(line,*) gromol%settles%ow,   &
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
    character(100)           :: error

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

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

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
    character(MaxLine)       :: error

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

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

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
    character(100)           :: error

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

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

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
    character(100)           :: error

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

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

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
  !! @param[inout] grotop : GROMACS TOP information
  !! @param[inout] gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_position_restraints(file, grotop, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(inout) :: grotop
    type(s_grotop_mol),      pointer       :: gromol

    ! local variables
    integer                  :: aidx
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line
    character(100)           :: error
    character(6)             :: at1, at2

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

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

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
    integer                  :: directive, i, cnt, sz, mol_cnt
    character(MaxLine)       :: line
    character(100)           :: error
    character(20)            :: name


    ! checko count
    !
    cnt = check_section_count(file)

    if (VerboseOn .and. main_rank) &
      write(MsgOut,'(" Read_Grotop> [molecules]")')

    do i = 1, cnt

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

      read(line,*,err=900,end=900) name, mol_cnt

      sz = size_grotop(grotop, GroTopMols) + 1
      call realloc_grotop(grotop, GroTopMols, sz)
      grotop%molss(sz)%name = name
      grotop%molss(sz)%count = mol_cnt

      if (VerboseOn .and. main_rank) &
        write(MsgOut,'(" Read_Grotop>   ",a20,"  count:",i5)') name, mol_cnt
    end do

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [molecules]')

  end subroutine read_molecules

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_gomodel
  !> @brief        read section [gomodel]
  !! @authors      CK
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_gomodel(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),  target, intent(inout) :: grotop

    ! local variables
    character(100)           :: line, error


    if (.not. gro_pp_next_line(file, line, error)) &
      goto 900
    read(line,*,err=900,end=900) grotop%gomodel%go_func,  &
                                 grotop%gomodel%num_basins

    return

900 call error_msg_grotop(file, 'read_grotop> read error. [gomodel]')

  end subroutine read_gomodel

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_cgatoms
  !> @brief        read section [cg_atoms]
  !! @authors      CK
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !! @param[inout] gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_cgatoms(file, grotop, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(inout) :: grotop
    type(s_grotop_mol),      pointer       :: gromol

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line
    character(100)           :: error
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
      write(MsgOut,'(" Read_Grotop>   [cg_atoms] :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt

    ! read data
    !
    do i = 1, cnt

      at => gromol%atoms(old_cnt+i)

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

      if (match_format(line, 'NSNSSNNNNN')) then

        read(line,*)          &
             at%atom_no,      &
             read_name,       &
             at%residue_no,   &
             at%residue_name, &
             at%atom_name,    &
             at%charge_group, &
             at%charge,       &
             at%mass,         &
             at%v,            &
             at%w

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

        if (.not. search_atom_vw(grotop%atomtypes, at%atom_type,             &
                                 at%v, at%w)) &
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

        if (.not. search_atom_vw(grotop%atomtypes, at%atom_type,             &
                                 at%v, at%w)) &
          call error_msg_grotop(file, 'Read_Grotop> [atoms] '// &
               'atom type not found: '//at%atom_type)
      else

        goto 900

      end if

    end do

    gromol%num_atoms = size_grotop_mol(gromol, GMGroMolAtom)

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [cg_atoms]')

  end subroutine read_cgatoms


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_multicontact
  !> @brief        read section [multicontact]
  !! @authors      CK
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !! @param[inout] gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_multicontact(file, grotop, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(inout) :: grotop
    type(s_grotop_mol),      pointer       :: gromol

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line
    character(100)           :: error
    character(6)             :: at1, at2

    type(s_mcont),   pointer :: pd


    ! check count
    !
    cnt = check_section_count(file)

    ! allocate memory
    !
    old_cnt = size_grotop_mol(gromol, GMGroMolMcont)
    call realloc_grotop_mol(gromol, GMGroMolMcont, old_cnt+cnt)

    if (VerboseOn .and. main_rank) &
      write(MsgOut,'(" Read_Grotop>   [multicontact] :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt

    ! read data
    !
    do i = 1, cnt

      pd => gromol%mcontact(old_cnt+i)

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

      if (match_format(line, 'NN2NNN')) then

        read(line,*) &
             pd%atom_idx1, &
             pd%atom_idx2, &
             pd%func,      &
             pd%model,     &
             pd%v,         &
             pd%w

      else

        goto 900

      end if

    end do

    gromol%num_mcontacts = size_grotop_mol(gromol, GMGroMolMcont)

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [multicontact]')

  end subroutine read_multicontact

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_morph_bb
  !> @brief        read section [moprh_bb]
  !! @authors      CK
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !! @param[inout] gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_morph_bb(file, grotop, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(inout) :: grotop
    type(s_grotop_mol),      pointer       :: gromol

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line
    character(100)           :: error
    character(6)             :: at1, at2

    type(s_morph_pair), pointer :: pd


    ! check count
    !
    cnt = check_section_count(file)

    ! allocate memory
    !
    old_cnt = size_grotop_mol(gromol, GMGroMolMorphBB)
    call realloc_grotop_mol(gromol, GMGroMolMorphBB, old_cnt+cnt)

    if (VerboseOn .and. main_rank) &
      write(MsgOut,'(" Read_Grotop>   [morph_bb] :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt

    ! read data
    !
    do i = 1, cnt

      pd => gromol%morph_bb(old_cnt+i)

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

      if (match_format(line, 'NNN')) then

        read(line,*) &
             pd%atom_idx1, &
             pd%atom_idx2, &
             pd%rmin

      else if (match_format(line, 'NNNN')) then

        read(line,*) &
             pd%atom_idx1, &
             pd%atom_idx2, &
             pd%rmin,      &
             pd%rmin_other
      else

        goto 900

      end if

    end do

    gromol%num_morph_bb = size_grotop_mol(gromol, GMGroMolMorphBB)

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [morph_bb]')

  end subroutine read_morph_bb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_morph_sc
  !> @brief        read section [moprh_sc]
  !! @authors      CK
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !! @param[inout] gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_morph_sc(file, grotop, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(inout) :: grotop
    type(s_grotop_mol),      pointer       :: gromol

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line
    character(100)           :: error
    character(6)             :: at1, at2

    type(s_morph_pair), pointer :: pd


    ! check count
    !
    cnt = check_section_count(file)

    ! allocate memory
    !
    old_cnt = size_grotop_mol(gromol, GMGroMolMorphSC)
    call realloc_grotop_mol(gromol, GMGroMolMorphSC, old_cnt+cnt)

    if (VerboseOn .and. main_rank) &
      write(MsgOut,'(" Read_Grotop>   [morph_sc] :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt

    ! read data
    !
    do i = 1, cnt

      pd => gromol%morph_sc(old_cnt+i)

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

      if (match_format(line, 'NNN')) then

        read(line,*) &
             pd%atom_idx1, &
             pd%atom_idx2, &
             pd%rmin

      else if (match_format(line, 'NNNN')) then

        read(line,*) &
             pd%atom_idx1, &
             pd%atom_idx2, &
             pd%rmin,      &
             pd%rmin_other

      else

        goto 900

      end if

    end do

    gromol%num_morph_sc = size_grotop_mol(gromol, GMGroMolMorphSC)

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [morph_sc]')

  end subroutine read_morph_sc


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_pwmcos
  !> @brief        read section [ pwmcos ]
  !! @authors      CT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !! @param[inout] gromol : GROMACS itp information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_pwmcos(file, grotop, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),  target, intent(inout) :: grotop
    type(s_grotop_mol), pointer :: gromol

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line
    character(100)           :: error
    ! character                :: chr_tmp_hyphen
    ! character                :: chr_tmp_intermol
    ! character(3)             :: chr_tmp_switch

    type(s_pwmcos),         pointer              :: pc


    ! check count
    cnt = check_section_count(file)

    ! allocate memory

    old_cnt = size_grotop_mol(gromol, GMGroMolPWMcos)
    call realloc_grotop_mol(gromol, GMGroMolPWMcos, old_cnt+cnt)

    if (VerboseOn .and. main_rank) then
      write(MsgOut,'(" Read_Grotop> [ pwmcos ] :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt
    end if

    ! read data
    !
    do i = 1, cnt

      pc => gromol%pwmcos(old_cnt+i)

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

      if (match_format(line, 'N1NNNNNNNNNN')) then
        read(line,*)        &
            pc%protein_idx, &
            pc%func,        &
            pc%r0,          &
            pc%theta1,      &
            pc%theta2,      &
            pc%theta3,      &
            pc%ene_A,       &
            pc%ene_C,       &
            pc%ene_G,       &
            pc%ene_T,       &
            pc%gamma,       &
            pc%eps_shift
      ! else if (match_format(line, 'N2NNNNN')) then
      !   read(line,*)        &
      !       pc%protein_idx, &
      !       pc%func,        &
      !       pc%r0,          &
      !       pc%theta1,      &
      !       pc%theta2,      &
      !       pc%theta3,      &
      !       pc%ene_A,       &
      !       pc%ene_C,       &
      !       pc%ene_G,       &
      !       pc%ene_T,       &
      !       pc%gamma,       &
      !       pc%eps_shift
      else
        call error_msg_grotop(file, 'Read_Grotop> format error. [ pwmcos ]')
      end if
    end do

    gromol%num_pwmcos = size_grotop_mol(gromol, GMGroMolPWMcos)

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [ pwmcos ]')

  end subroutine read_pwmcos


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_pwmcosns
  !> @brief        read section [ pwmcosns ]
  !! @authors      CT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !! @param[inout] gromol : GROMACS itp information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_pwmcosns(file, grotop, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),  target, intent(inout) :: grotop
    type(s_grotop_mol),      pointer       :: gromol

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line
    character(100)           :: error
    ! character                :: chr_tmp_hyphen
    ! character                :: chr_tmp_intermol
    ! character(3)             :: chr_tmp_switch

    type(s_pwmcosns), pointer :: pc


    ! check count
    cnt = check_section_count(file)

    ! allocate memory

    old_cnt = size_grotop_mol(gromol, GMGroMolPWMcosns)
    call realloc_grotop_mol(gromol, GMGroMolPWMcosns, old_cnt+cnt)

    if (VerboseOn .and. main_rank) then
      write(MsgOut,'(" Read_Grotop> [ pwmcosns ] :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt
    end if

    ! read data
    !
    do i = 1, cnt

      pc => gromol%pwmcosns(old_cnt+i)

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

      if (match_format(line, 'N2NNNN')) then
        read(line,*)        &
            pc%protein_idx, &
            pc%func,        &
            pc%r0,          &
            pc%theta1,      &
            pc%theta2,      &
            pc%ene
      else
        call error_msg_grotop(file, 'Read_Grotop> format error. [ pwmcosns ]')
      end if
    end do

    gromol%num_pwmcosns = size_grotop_mol(gromol, GMGroMolPWMcosns)

    return

900 call error_msg_grotop(file, 'Read_Grotop> read error. [ pwmcosns ]')

  end subroutine read_pwmcosns


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_idr_hps
  !> @brief        read section [ idr_hps ]
  !! @authors      CT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !! @param[inout] gromol : GROMACS ITP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_idr_hps(file, grotop, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),  target, intent(inout) :: grotop
    type(s_grotop_mol),      pointer       :: gromol

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line
    character(100)           :: error
    character(3)             :: chr_tmp

    type(s_idr_hps), pointer :: hps


    ! check count
    cnt = check_section_count(file)

    ! allocate memory
    old_cnt = size_grotop_mol(gromol, GMGroMolIDRHPS)
    call realloc_grotop_mol(gromol, GMGroMolIDRHPS, old_cnt+cnt)

    if (VerboseOn .and. main_rank) then
      write(MsgOut,'(" Read_Gromol> [ cg_IDR_HPS_region ] proper   :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt
    end if

    ! read data
    !
    do i = 1, cnt

      hps => gromol%idr_hps(old_cnt+i)

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

      if (match_format(line, 'NN')) then
        read(line,*) hps%grp_start, hps%grp_end
      else if (match_format(line, 'N')) then
        read(line,*) hps%grp_start
        hps%grp_end = hps%grp_start
      else
        call error_msg_grotop(file, 'Read_Gromol> read error. [ cg_IDR_HPS_region ]')
      end if
    end do

    gromol%num_idr_hps = size_grotop_mol(gromol, GMGroMolIDRHPS)

    return

900 call error_msg_grotop(file, 'Read_Gromol> read error. [ cg_IDR_HPS_region ]')

  end subroutine read_idr_hps


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_idr_kh
  !> @brief        read section [ idr_kh ]
  !! @authors      CT
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !! @param[inout] gromol : GROMACS ITP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_idr_kh(file, grotop, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),  target, intent(inout) :: grotop
    type(s_grotop_mol),      pointer       :: gromol

    ! local variables
    integer                  :: i, cnt, old_cnt
    character(MaxLine)       :: line
    character(100)           :: error
    character(3)             :: chr_tmp

    type(s_idr_kh), pointer :: kh


    ! check count
    cnt = check_section_count(file)

    ! allocate memory
    old_cnt = size_grotop_mol(gromol, GMGroMolIDRKH)
    call realloc_grotop_mol(gromol, GMGroMolIDRKH, old_cnt+cnt)

    if (VerboseOn .and. main_rank) then
      write(MsgOut,'(" Read_Gromol> [ cg_IDR_KH_region ] proper   :"'// &
           ',i5," (total:",i5,")")') cnt, old_cnt+cnt
    end if

    ! read data
    !
    do i = 1, cnt

      kh => gromol%idr_kh(old_cnt+i)

      if (.not. gro_pp_next_line(file, line, error)) &
        goto 900

      if (match_format(line, 'NN')) then
        read(line,*) kh%grp_start, kh%grp_end
      else if (match_format(line, 'N')) then
        read(line,*) kh%grp_start
        kh%grp_end = kh%grp_start
      else
        call error_msg_grotop(file, 'Read_Gromol> read error. [ cg_IDR_KH_region ]')
      end if
    end do

    gromol%num_idr_kh = size_grotop_mol(gromol, GMGroMolIDRKH)

    return

900 call error_msg_grotop(file, 'Read_Gromol> read error. [ cg_IDR_KH_region ]')

  end subroutine read_idr_kh



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

    ! basestack type
    call write_base_stack_types(file, grotop)

    ! basepair type
    call write_base_pair_types(file, grotop)

    ! basecross type
    call write_base_cross_types(file, grotop)

    ! cgdnaexv type
    call write_base_exv_types(file, grotop)

    ! cg IDR HPS atom type
    call write_cg_IDR_HPS_atom_types(file, grotop)

    ! cg KH atom type
    call write_cg_KH_atom_types(file, grotop)

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
  !  Subroutine    write_grotop_aadome
  !> @brief        write a GROMACS TOP file
  !! @authors      CK
  !! @param[in]    file   : unit number of GROTOP file
  !! @param[in]    grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_grotop_aadome(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(in)    :: grotop

    ! itp files
    call write_itp(file, grotop)

    call write_gomodel(file, grotop)

    ! defaults
    call write_defaults(file, grotop)

!    ! atom type
!    call write_atom_types(file, grotop)

    ! bond type
    call write_bond_types_aa(file, grotop)

    ! angle type
    call write_angle_types_aa(file, grotop)

    ! dihedral type
    call write_dihedral_types_aa(file, grotop)

    ! constraint type
    call write_constraint_types(file, grotop)

    ! nonbond parameter
    call write_nonbond_parms(file, grotop)

    ! molecule type
    call write_molecule_type_aa(file, grotop)

    ! system
    call write_system(file, grotop)

    ! molecules
    call write_molecules(file, grotop)

    return

  end subroutine write_grotop_aadome

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
  !  Subroutine    write_gomodel
  !> @brief        write section [gomodel]
  !! @authors      CK
  !! @param[in]    file   : file unit number
  !! @param[inout] grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_gomodel(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),  target, intent(in)    :: grotop

    ! local variables
    character(100)           :: line


    write(file, '(a)') &
         '[ gomodel ]'
    write(file, '(a)') &
         '; nbfunc  num_basins'
    write(file,'(i10,1x,i10)') grotop%gomodel%go_func,  &
                               grotop%gomodel%num_basins

    write(file, '(a)') ''

    return

  end subroutine write_gomodel

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_itp
  !> @brief        write section for itp files
  !! @authors      CK
  !! @param[in]    file   : file unit number
  !! @param[in]    grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_itp(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(in)    :: grotop


    if (grotop%itpread%itp /= '') then
      write(file, '(a,a,a)') '#include "',trim(adjustl(grotop%itpread%itp)),'"'

      write(file, '(a)') ''
    end if

    return

  end subroutine write_itp

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
  !  Subroutine    write_bond_types_aa
  !> @brief        write section [bondtypes]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_bond_types_aa(file, grotop)

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

      write(file, '(a8,1x,a8,1x,i8)') &
           grotop%bondtypes(i)%atom_type1,  &
           grotop%bondtypes(i)%atom_type2,  &
           grotop%bondtypes(i)%func

    end do

    write(file, '(a)') ''

    return

  end subroutine write_bond_types_aa

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
  !  Subroutine    write_angle_types_aa
  !> @brief        write section [angletypes]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_angle_types_aa(file, grotop)

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

        case default
          write(str,*) func
          call error_msg('Write_Angle_Types> Unknown function type.'//trim(str))

        end select

      end if

      write(file, '(a8,1x,a8,1x,a8,1x,i8)') &
           grotop%angltypes(i)%atom_type1,     &
           grotop%angltypes(i)%atom_type2,     &
           grotop%angltypes(i)%atom_type3,     &
           grotop%angltypes(i)%func

    end do

    write(file, '(a)') ''

    return

  end subroutine write_angle_types_aa


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
  !  Subroutine    write_dihedral_types_aa
  !> @brief        write section [dihedraltypes]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_dihedral_types_aa(file, grotop)

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

        ! dihedral
        write(file,'(a8,1x,a8,1x,a8,1x,a8,1x,i8)') &
             grotop%dihetypes(i)%atom_type1, &
             grotop%dihetypes(i)%atom_type2, &
             grotop%dihetypes(i)%atom_type3, &
             grotop%dihetypes(i)%atom_type4, &
             grotop%dihetypes(i)%func

      else

        ! dihedral
        write(file,'(a8,1x,a8,1x,i8)') &
             grotop%dihetypes(i)%atom_type2, &
             grotop%dihetypes(i)%atom_type3, &
             grotop%dihetypes(i)%func

      end if

    end do

    write(file, '(a)') ''

    return

  end subroutine write_dihedral_types_aa

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
  !  Subroutine    write_base_stack_types
  !> @brief        write section [ basestacktypes ]
  !! @authors      CT
  !! @param[in]    file   : file unit number
  !! @param[in]    grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_base_stack_types(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(in)    :: grotop

    ! local variables
    integer                  :: i, cnt, func
    character(20)            :: str


    cnt  = size_grotop(grotop, GroTopBaseStackType)
    if (cnt == 0) &
      return

    func = -1

    do i = 1, cnt

      if (func /= grotop%basestacktypes(i)%func) then
        func = grotop%basestacktypes(i)%func
        write(file,'(A)') '[ basestacktypes ]'
        write(file,'(A)') '; DNA - Base Stacking'
        write(file,'(A)') ";   5\'    3\' func   epsilon     sigma thetha_bs"
      end if

      write(file,'(a6,a6,i5,f10.2,f10.3,f10.2)') &
        grotop%basestacktypes(i)%base_type5, &
        grotop%basestacktypes(i)%base_type3, &
        grotop%basestacktypes(i)%func,       &
        grotop%basestacktypes(i)%epsilon,    &
        grotop%basestacktypes(i)%sigma,      &
        grotop%basestacktypes(i)%theta_bs

    end do

    write(file, '(a)') ''

    return

  end subroutine write_base_stack_types

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_base_pair_types
  !> @brief        write section [ basepairtypes ]
  !! @authors      CT
  !! @param[in]    file   : file unit number
  !! @param[in]    grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_base_pair_types(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(in)    :: grotop

    ! local variables
    integer                  :: i, cnt, func
    character(20)            :: str


    cnt  = size_grotop(grotop, GroTopBasePairType)
    if (cnt == 0) &
      return

    func = -1

    do i = 1, cnt

      if (func /= grotop%basepairtypes(i)%func) then
        func = grotop%basepairtypes(i)%func
        write(file, '(a)') '[ basepairtypes ]'
        write(file,'(a)') '; DNA - Base Pairing'
        write(file,'(a)') ';    a     b func   thetha1   thetha2      phi1     sigma   epsilon'
      end if

      write(file,'(a6,a6,i5,f10.2,f10.2,f10.2,f10.2,f10.2)') &
          grotop%basepairtypes(i)%base_type_a,               &
          grotop%basepairtypes(i)%base_type_b,               &
          grotop%basepairtypes(i)%func,                      &
          grotop%basepairtypes(i)%theta_1,                   &
          grotop%basepairtypes(i)%theta_2,                   &
          grotop%basepairtypes(i)%phi_1,                     &
          grotop%basepairtypes(i)%sigma,                     &
          grotop%basepairtypes(i)%epsilon

    end do

    write(file, '(a)') ''

    return

  end subroutine write_base_pair_types

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_base_cross_types
  !> @brief        write section [ basecrosstypes ]
  !! @authors      CT
  !! @param[in]    file   : file unit number
  !! @param[in]    grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_base_cross_types(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(in)    :: grotop

    ! local variables
    integer                  :: i, cnt, func
    character(20)            :: str


    cnt  = size_grotop(grotop, GroTopBaseCrossType)
    if (cnt == 0) &
      return

    func = -1

    write(file, '(a)') '[ basecrosstypes ]'

    do i = 1, cnt

      if (func /= grotop%basecrosstypes(i)%func) then
        func = grotop%basecrosstypes(i)%func
        select case (func)
        case (1)
          write(file,'(a)') ";DNA - Base Cross-Stacking up5\'-down5\'"
        case (2)
          write(file,'(a)') ";DNA - Base Cross-Stacking down3\'-up3\'"
        case default
          call error_msg('Write_Base_Cross_Types> Unknown function type.')
        end select
        write(file,'(a)') ';    a     b func   epsilon     sigma   thetha3 thetha_cs'
      end if

      write(file,'(a6,a6,i5,f10.3,f10.3,f10.2,f10.2)') &
        grotop%basecrosstypes(i)%base_type_a,      &
        grotop%basecrosstypes(i)%base_type_b,      &
        grotop%basecrosstypes(i)%func,             &
        grotop%basecrosstypes(i)%epsilon,          &
        grotop%basecrosstypes(i)%sigma,            &
        grotop%basecrosstypes(i)%theta_cs

    end do

    write(file, '(a)') ''

    return

  end subroutine write_base_cross_types

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_base_exv_types
  !> @brief        write section [ cgdnaexvtypes ]
  !! @authors      CT
  !! @param[in]    file   : file unit number
  !! @param[in]    grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_base_exv_types(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(in)    :: grotop

    ! local variables
    integer                  :: i, cnt, func
    character(20)            :: str


    cnt  = size_grotop(grotop, GroTopCGDNAExvType)
    if (cnt == 0) &
      return

    func = -1

    do i = 1, cnt

      if (func /= grotop%cgdnaexvtypes(i)%func) then
        func = grotop%cgdnaexvtypes(i)%func
        write(file, '(a)') '[ cgdnaexvtypes ]'
        write(file,'(a)') '; DNA - Base Exving'
        write(file,'(a)') ';    a func     sigma'
      end if

      write(file,'(a6,i5,f10.2)')           &
          grotop%cgdnaexvtypes(i)%base_type, &
          grotop%cgdnaexvtypes(i)%func,      &
          grotop%cgdnaexvtypes(i)%sigma

    end do

    write(file, '(a)') ''

    return

  end subroutine write_base_exv_types

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_cg_IDR_HPS_atom_types
  !> @brief        write section [ cg_IDR_HPS_atomtypes ]
  !! @authors      CT
  !! @param[in]    file   : file unit number
  !! @param[in]    grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_cg_IDR_HPS_atom_types(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(in)    :: grotop

    ! local variables
    integer                  :: i, cnt, func
    character(20)            :: str


    cnt  = size_grotop(grotop, GroTopCGIDRHPSAtomType)
    if (cnt == 0) &
      return

    func = -1

    write(file, '(a)') '[ cg_IDR_HPS_atomtypes ]'

    do i = 1, cnt

      if (func == -1) then
        write(file,'(a)') ';name      mass  charge  sigma   lambda'
        func = 1
      end if

      write(file,'(a5,f10.3,f8.1,f7.2,f8.3)')     &
        grotop%cg_IDR_HPS_atomtypes(i)%type_name, &
        grotop%cg_IDR_HPS_atomtypes(i)%mass,      &
        grotop%cg_IDR_HPS_atomtypes(i)%charge,    &
        grotop%cg_IDR_HPS_atomtypes(i)%sigma,     &
        grotop%cg_IDR_HPS_atomtypes(i)%lambda

    end do

    write(file, '(a)') ''

    return

  end subroutine write_cg_IDR_HPS_atom_types


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_cg_KH_atom_types
  !> @brief        write section [ cg_KH_atomtypes ]
  !! @authors      CT
  !! @param[in]    file   : file unit number
  !! @param[in]    grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_cg_KH_atom_types(file, grotop)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop),          intent(in)    :: grotop

    ! local variables
    integer                  :: i, cnt, func
    character(20)            :: str


    cnt  = size_grotop(grotop, GroTopCGKHAtomType)
    if (cnt == 0) &
      return

    func = -1

    write(file, '(a)') '[ cg_KH_atomtypes ]'

    do i = 1, cnt

      if (func == -1) then
        write(file,'(a)') ';name      mass  charge  sigma   lambda'
        func = 1
      end if

      write(file,'(a5,f10.3,f8.1,f8.1)')     &
        grotop%cg_KH_atomtypes(i)%type_name, &
        grotop%cg_KH_atomtypes(i)%mass,      &
        grotop%cg_KH_atomtypes(i)%charge,    &
        grotop%cg_KH_atomtypes(i)%sigma

    end do

    write(file, '(a)') ''

    return

  end subroutine write_cg_KH_atom_types


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

    ! constants
    integer,       parameter :: MaxFroms = 100

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

      ! morph pairs_bb
      call write_morph_pairs_bb(file, grotop%moltypes(i)%mol)

      ! morph pairs_sc
      call write_morph_pairs_sc(file, grotop%moltypes(i)%mol)

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
  !  Subroutine    write_molecule_type_aa
  !> @brief        write section [moleculetype]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    grotop : GROMACS TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_molecule_type_aa(file, grotop)

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
      call write_bonds_aa(file, grotop%moltypes(i)%mol)

      ! angles
      call write_angles_aa(file, grotop%moltypes(i)%mol)

      ! dihedrals
      call write_dihedrals_aa(file, grotop%moltypes(i)%mol)

      ! exclusions
      call write_exclusions(file, grotop%moltypes(i)%mol)

      ! constraints
      call write_constraints(file, grotop%moltypes(i)%mol)

      ! pairs
      call write_pairs(file, grotop%moltypes(i)%mol)

      ! multi-pairs
      call write_multicontact(file, grotop%moltypes(i)%mol)

      ! morph pairs bb
      call write_morph_pairs_bb(file, grotop%moltypes(i)%mol)

      ! morph pairs sc
      call write_morph_pairs_sc(file, grotop%moltypes(i)%mol)

      ! settles
      call write_settles(file, grotop%moltypes(i)%mol)

      ! position restraints
      call write_position_restraints(file, grotop%moltypes(i)%mol)

    end do

    return

  end subroutine write_molecule_type_aa

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

!!!original
!!!!!!    if (mass) then
!!!original
!!!develop
    if (gromol%atoms(1)%mass /= 0.0_wp) then
!!!develop
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

!!!original
!!!!!!      if (mass) then
!!!original
!!!develop
      if (gromol%atoms(i)%mass /= 0.0_wp) then
!!!develop

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
  !  Subroutine    write_bonds_aa
  !> @brief        write section [bonds]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_bonds_aa(file, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop_mol),      intent(in)    :: gromol

    ! local variables
    integer                  :: i, cnt


    cnt = size_grotop_mol(gromol, GMGroMolBond)
    if (cnt == 0) &
      return

    write(file,'(a)') '[ bonds ]'
    write(file,'(a)') '; i  j  funct'

    do i = 1, cnt

      write(file,'(i8,1x,i8,1x,i8)') &
           gromol%bonds(i)%atom_idx1, &
           gromol%bonds(i)%atom_idx2, &
           gromol%bonds(i)%func

    end do

    write(file,'(a)') ''

    return

  end subroutine write_bonds_aa

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
!!! value i is not initilized in GET_CG_2018

    i = 1
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
  !  Subroutine    write_angles_aa
  !> @brief        write section [angles]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_angles_aa(file, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop_mol),      intent(in)    :: gromol

    ! local variables
    integer                  :: i, cnt


    cnt = size_grotop_mol(gromol, GMGroMolAngl)
    if (cnt == 0) &
      return

    write(file,'(a)') '[ angles ]'
    write(file,'(a)') '; i  j  k  funct'

    do i = 1, cnt

      write(file,'(i8,1x,i8,1x,i8,1x,i8)') &
           gromol%angls(i)%atom_idx1, &
           gromol%angls(i)%atom_idx2, &
           gromol%angls(i)%atom_idx3, &
           gromol%angls(i)%func

    end do

    write(file,'(a)') ''

    return

  end subroutine write_angles_aa

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
  !  Subroutine    write_dihedrals_aa
  !> @brief        write section [dihedrals]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_dihedrals_aa(file, gromol)

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

      ! dihedral
    write(file,'(a)') '; ai  aj  ak  al  funct'

    do i = 1, cnt

      write(file,'(i8,1x,i8,1x,i8,1x,i8,1x,i8)') &
           gromol%dihes(i)%atom_idx1, &
           gromol%dihes(i)%atom_idx2, &
           gromol%dihes(i)%atom_idx3, &
           gromol%dihes(i)%atom_idx4, &
           gromol%dihes(i)%func

    end do

    write(file,'(a)') ''

    return

  end subroutine write_dihedrals_aa

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
    logical                  :: mass, stok


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
  !  Subroutine    write_multicontact
  !> @brief        write section [pairs]
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_multicontact(file, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop_mol),      intent(in)    :: gromol

    ! local variables
    integer                  :: i, cnt


    cnt = size_grotop_mol(gromol, GMGroMolMcont)
    if (cnt == 0) &
      return

    write(file,'(a)') '[ multicontact ]'
    write(file,'(a)') '; i  j  func  model v  w'

    do i = 1, cnt

      write(file,'(i8,1x,i8,1x,i8,1x,i8,1x,e16.9,1x,e16.9)') &
           gromol%mcontact(i)%atom_idx1, &
           gromol%mcontact(i)%atom_idx2, &
           gromol%mcontact(i)%func,      &
           gromol%mcontact(i)%model,     &
           gromol%mcontact(i)%v,         &
           gromol%mcontact(i)%w

    end do

    write(file,'(a)') ''

    return

  end subroutine write_multicontact

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_morph_pairs_bb
  !> @brief        write section [morph_bb]
  !! @authors      CK
  !! @param[in]    file   : file unit number
  !! @param[in]    gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_morph_pairs_bb(file, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop_mol),      intent(in)    :: gromol

    ! local variables
    integer                  :: i, cnt


    cnt = size_grotop_mol(gromol, GMGroMolMorphBB)
    if (cnt == 0) &
      return

    write(file,'(a)') '[ morph_bb ]'
    write(file,'(a)') '; i  j  length '

    do i = 1, cnt

      write(file,'(i8,1x,i8,1x,e16.9)') &
           gromol%morph_bb(i)%atom_idx1, &
           gromol%morph_bb(i)%atom_idx2, &
           gromol%morph_bb(i)%rmin,    &
           gromol%morph_bb(i)%rmin_other

    end do

    write(file,'(a)') ''

    return

  end subroutine write_morph_pairs_bb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_morph_pairs_sc
  !> @brief        write section [morph_sc]
  !! @authors      CK
  !! @param[in]    file   : file unit number
  !! @param[in]    gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_morph_pairs_sc(file, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop_mol),      intent(in)    :: gromol

    ! local variables
    integer                  :: i, cnt
    logical                  :: mass, stok


    cnt = size_grotop_mol(gromol, GMGroMolMorphSC)
    if (cnt == 0) &
      return

    write(file,'(a)') '[ morph_sc ]'
    write(file,'(a)') '; i  j  length '

    do i = 1, cnt

      write(file,'(i8,1x,i8,1x,e16.9)') &
           gromol%morph_sc(i)%atom_idx1, &
           gromol%morph_sc(i)%atom_idx2, &
           gromol%morph_sc(i)%rmin,      &
           gromol%morph_bb(i)%rmin_other

    end do

    write(file,'(a)') ''

    return

  end subroutine write_morph_pairs_sc


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_pwmcos
  !> @brief        write section [ pwmcos ]
  !! @authors      CT
  !! @param[in]    file   : file unit number
  !! @param[in]    gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_pwmcos(file, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop_mol),      intent(in)    :: gromol

    ! local variables
    integer                  :: i, cnt


    cnt = size_grotop_mol(gromol, GMGroMolPWMcos)
    if (cnt == 0) &
      return

    write(file,'(a)') '[ pwmcos ]'
    write(file,'(a)') '; i  j  length '

    do i = 1, cnt

      write(file,'(i8,1x,i8,1x,10e15.4)') &
          gromol%pwmcos(i)%protein_idx,   &
          gromol%pwmcos(i)%func,          &
          gromol%pwmcos(i)%r0,            &
          gromol%pwmcos(i)%theta1,        &
          gromol%pwmcos(i)%theta2,        &
          gromol%pwmcos(i)%theta3,        &
          gromol%pwmcos(i)%ene_A,         &
          gromol%pwmcos(i)%ene_C,         &
          gromol%pwmcos(i)%ene_G,         &
          gromol%pwmcos(i)%ene_T,         &
          gromol%pwmcos(i)%gamma,         &
          gromol%pwmcos(i)%eps_shift

    end do

    write(file,'(a)') ''

    return

  end subroutine write_pwmcos


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_pwmcosns
  !> @brief        write section [ pwmcosns ]
  !! @authors      CT
  !! @param[in]    file   : file unit number
  !! @param[in]    gromol : GROMACS TOP molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_pwmcosns(file, gromol)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grotop_mol),      intent(in)    :: gromol

    ! local variables
    integer                  :: i, cnt


    cnt = size_grotop_mol(gromol, GMGroMolPWMcosns)
    if (cnt == 0) &
      return

    write(file,'(a)') '[ pwmcosns ]'
    write(file,'(a)') '; i  j  length '

    do i = 1, cnt

      write(file,'(i8,1x,i8,1x,4e15.4)') &
          gromol%pwmcosns(i)%protein_idx,   &
          gromol%pwmcosns(i)%func,          &
          gromol%pwmcosns(i)%r0,            &
          gromol%pwmcosns(i)%theta1,        &
          gromol%pwmcosns(i)%theta2,        &
          gromol%pwmcosns(i)%ene

    end do

    write(file,'(a)') ''

    return

  end subroutine write_pwmcosns


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
    character(100)           :: error


    cnt = 0
    call gro_pp_record_pos(file)

    do while(.true.)
      if (.not. gro_pp_next_line(file, line, error)) &
        exit
      if (check_directive(line, directive)) &
        exit
      cnt = cnt + 1
    end do

    call gro_pp_restore_pos(file)

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
    character(30)            :: dir_str


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
    case ('flexible_local_angle')
      directive = DFlexible_Local_Angle
    case ('dihedraltypes')
      directive = DDihedralTypes
    case ('flexible_local_dihedral_angle')
      directive = DFlexible_Local_Dihedral_Angle
    case ('constrainttypes')
      directive = DConstrTypes
    case ('cmaptypes')
      directive = DCmapTypes
    case ('basestacktypes')
      directive = DBaseStackTypes
    case ('basepairtypes')
      directive = DBasePairTypes
    case ('basecrosstypes')
      directive = DBaseCrossTypes
    case ('cgdnaexvtypes')
      directive = DCGDNAExvTypes
    case ('cg_ele_chain_pairs')
      directive = DCGEleMolPairs
    case ('cg_KH_chain_pairs')
      directive = DCGKHMolPairs
    case ('pwmcos')
      directive = DCGPWMcos
    case ('pwmcosns')
      directive = DCGPWMcosns
    case ('cg_IDR_HPS_atomtypes')
      directive = DCGIDRHPSAtomTypes
    case ('cg_IDR_HPS_region')
      directive = DCGIDRHPSRegion
    case ('cg_KH_atomtypes')
      directive = DCGKHAtomTypes
    case ('cg_IDR_KH_region')
      directive = DCGIDRKHRegion
    case ('cg_pair_energy_MJ_96')
      directive = DCGPairMJEpsilon
    case ('pwmcos_chain_pairs')
      directive = DCGPWMcosMolPairs
    case ('pwmcosns_chain_pairs')
      directive = DCGPWMcosnsMolPairs
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
    case ('gomodel')
      directive = DGomodel
    case ('cg_atoms')
      directive = DCg_Atoms
    case ('multicontact')
      directive = DMulticontact
    case ('morph_bb')
      directive = DMorphBB
    case ('morph_sc')
      directive = DMorphSC
    case ('ez_membrane')
      directive = DEzMembrane
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

    case(GroTopFlAnglType)
      if (allocated(grotop%flangltypes)) &
        size_grotop = size(grotop%flangltypes)

    case(GroTopDiheType)
      if (allocated(grotop%dihetypes)) &
        size_grotop = size(grotop%dihetypes)

    case(GroTopFlDiheType)
      if (allocated(grotop%fldihetypes)) &
        size_grotop = size(grotop%fldihetypes)

    case(GroTopConstrType)
      if (allocated(grotop%constrtypes)) &
        size_grotop = size(grotop%constrtypes)

    case(GroTopCmapType)
      if (allocated(grotop%cmaptypes)) &
        size_grotop = size(grotop%cmaptypes)

    case(GroTopBaseStackType)
      if (allocated(grotop%basestacktypes)) &
        size_grotop = size(grotop%basestacktypes)

    case(GroTopBasePairType)
      if (allocated(grotop%basepairtypes)) &
        size_grotop = size(grotop%basepairtypes)

    case(GroTopBaseCrossType)
      if (allocated(grotop%basecrosstypes)) &
        size_grotop = size(grotop%basecrosstypes)

    case(GroTopCGDNAExvType)
      if (allocated(grotop%cgdnaexvtypes)) &
          size_grotop = size(grotop%cgdnaexvtypes)

    case(GroTopCGEleMolPair)
      if (allocated(grotop%cg_ele_mol_pairs)) &
          size_grotop = size(grotop%cg_ele_mol_pairs)

    case(GroTopCGKHMolPair)
      if (allocated(grotop%cg_kh_mol_pairs)) &
          size_grotop = size(grotop%cg_kh_mol_pairs)

    case(GroTopPWMcosMolPair)
      if (allocated(grotop%pwmcos_mol_pairs)) &
          size_grotop = size(grotop%pwmcos_mol_pairs)

    case(GroTopPWMcosnsMolPair)
      if (allocated(grotop%pwmcosns_mol_pairs)) &
          size_grotop = size(grotop%pwmcosns_mol_pairs)

    case(GroTopCGIDRHPSAtomType)
      if (allocated(grotop%cg_IDR_HPS_atomtypes)) &
          size_grotop = size(grotop%cg_IDR_HPS_atomtypes)

    case(GroTopCGKHAtomType)
      if (allocated(grotop%cg_KH_atomtypes)) &
          size_grotop = size(grotop%cg_KH_atomtypes)

    case(GroTopCGPairMJEpsilon)
      if (allocated(grotop%cg_pair_MJ_eps)) &
          size_grotop = size(grotop%cg_pair_MJ_eps)

    case(GroTopNbonParm)
      if (allocated(grotop%nbonparms)) &
        size_grotop = size(grotop%nbonparms)

    case(GroTopMolType)
      if (allocated(grotop%moltypes)) &
        size_grotop = size(grotop%moltypes)

    case(GroTopMols)
      if (allocated(grotop%molss)) &
        size_grotop = size(grotop%molss)

    case(GroTopEzMembrane)
      if (allocated(grotop%ez_membrane)) &
        size_grotop = size(grotop%ez_membrane)

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

    case (GMGroMolMcont)
      if (allocated(gromol%mcontact)) then
        size_grotop_mol = size(gromol%mcontact)
      else
        size_grotop_mol = 0
      end if

    case (GMGroMolPosres)
      if (allocated(gromol%posress)) then
        size_grotop_mol = size(gromol%posress)
      else
        size_grotop_mol =0
      end if

    case (GMGroMolMorphBB)
      if (allocated(gromol%morph_bb)) then
        size_grotop_mol = size(gromol%morph_bb)
      else
        size_grotop_mol = 0
      end if

    case (GMGroMolMorphSC)
      if (allocated(gromol%morph_sc)) then
        size_grotop_mol = size(gromol%morph_sc)
      else
        size_grotop_mol = 0
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
        size_grotop_mol = 0
      end if

    case (GMGroMolPWMcos)
      if (allocated(gromol%pwmcos)) then
        size_grotop_mol = size(gromol%pwmcos)
      else
        size_grotop_mol = 0
      end if

    case (GMGroMolPWMcosns)
      if (allocated(gromol%pwmcosns)) then
        size_grotop_mol = size(gromol%pwmcosns)
      else
        size_grotop_mol = 0
      end if

    case(GMGroMolIDRHPS)
      if (allocated(gromol%idr_hps)) then
        size_grotop_mol = size(gromol%idr_hps)
      else
        size_grotop_mol = 0
      end if

    case(GMGroMolIDRKH)
      if (allocated(gromol%idr_kh)) then
        size_grotop_mol = size(gromol%idr_kh)
      else
        size_grotop_mol = 0
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
  !  Function      search_atom_vw
  !> @brief        search atom vdw
  !! @authors      CK
  !! @param[in]    atomtypes : data of atomtypes
  !! @param[in]    at        : atom type name
  !! @param[out]   mass      : mass of atoms
  !! @return       flag for found or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function search_atom_vw(atomtypes, at, v, w)

    ! return value
    logical                  :: search_atom_vw

    ! formal arguments
    type(s_atomtype),        intent(in)    :: atomtypes(:)
    character(*),            intent(in)    :: at
    real(wp),                intent(out)   :: v
    real(wp),                intent(out)   :: w

    ! local variables
    integer                  :: i


    do i = 1, size(atomtypes)
      if (atomtypes(i)%type_name == at) then
        v = atomtypes(i)%v
        w = atomtypes(i)%w
        search_atom_vw = .true.
        return
      end if
    end do

    search_atom_vw = .false.
    return

  end function search_atom_vw

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      search_atom_charge
  !> @brief        search atom charge
  !! @authors      CK
  !! @param[in]    atomtypes : data of atomtypes
  !! @param[in]    at        : atom type name
  !! @param[out]   charge    : mass of atoms
  !! @return       flag for found or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function search_atom_charge(atomtypes, at, charge)

    ! return value
    logical                  :: search_atom_charge

    ! formal arguments
    type(s_atomtype),        intent(in)    :: atomtypes(:)
    character(*),            intent(in)    :: at
    real(wp),                intent(out)   :: charge

    ! local variables
    integer                  :: i


    do i = 1, size(atomtypes)
      if (atomtypes(i)%type_name == at) then
        charge = atomtypes(i)%charge
        search_atom_charge = .true.
        return
      end if
    end do

    search_atom_charge = .false.
    return

  end function search_atom_charge

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
  !  Function      search_flangl_param
  !> @brief        search flexible angle parameters
  !! @authors      CK
  !! @param[in]    flangltypes : data of flangletypes
  !! @param[in]    at1       : atom 1 type name
  !! @param[in]    at2       : atom 2 type name
  !! @param[in]    at3       : atom 3 type name
  !! @param[in]    func      : function type
  !! @param[out]   types     : param types
  !! @return       flag for found or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function search_flangl_param(flangltypes, at1, at2, at3, func, types)

    ! return value
    logical :: search_flangl_param

    ! formal arguments
    type(s_flangltype), target, intent(in)    :: flangltypes(:)
    character(*),               intent(in)    :: at1
    character(*),               intent(in)    :: at2
    character(*),               intent(in)    :: at3
    integer,                    intent(in)    :: func
    integer,                    intent(out)   :: types

    ! local variables
    integer                   :: i
    type(s_flangltype), pointer :: at


    do i = 1, size(flangltypes)
      at => flangltypes(i)
!      if (at%func /= func) &
!        cycle
      if (at%atom_type .ne. at2) &
        cycle

      types = i
      search_flangl_param = .true.
      return
    end do

    search_flangl_param = .false.
    return

  end function search_flangl_param

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
  !  Function      search_fldihe_param
  !> @brief        search flexible dihedral angle parameters
  !! @authors      JJ
  !! @param[in]    flangltypes : data of flangletypes
  !! @param[in]    at1       : atom 1 type name
  !! @param[in]    at2       : atom 2 type name
  !! @param[in]    at3       : atom 3 type name
  !! @param[in]    at4       : atom 4 type name
  !! @param[in]    func      : function type
  !! @param[out]   types     : param types
  !! @return       flag for found or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function search_fldihe_param(fldihetypes, at1, at2, at3, at4, func, types)

    ! return value
    logical :: search_fldihe_param

    ! formal arguments
    type(s_fldihetype), target, intent(in)    :: fldihetypes(:)
    character(*),               intent(in)    :: at1
    character(*),               intent(in)    :: at2
    character(*),               intent(in)    :: at3
    character(*),               intent(in)    :: at4
    integer,                    intent(in)    :: func
    integer,                    intent(out)   :: types

    ! local variables
    integer                   :: i
    type(s_fldihetype), pointer :: dh

    do i = 1, size(fldihetypes)
      dh => fldihetypes(i)
!      if (at%func /= func) &
!        cycle
      if (dh%atom_type1 .ne. at2 .or. dh%atom_type2 .ne. at3) &
        cycle

      types = i
      search_fldihe_param = .true.
      return
    end do

    search_fldihe_param = .false.
    return

  end function search_fldihe_param

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
    integer                   :: icnt(MAX_MULTIPLY_NUM)
    logical                   :: four_atoms

    icn = 0
    icnt(1:MAX_MULTIPLY_NUM) = 0
    iold = ind
    four_atoms = .false.

    do i = 1, size(dihetypes)
      dt => dihetypes(i)

      if (dt%func /= dihes(old+iold)%func) &
        cycle

      if (dt%wild_num == 0) then

        if (dt%four_type) then

          if ((dt%atom_type1 == at1 .and. &
               dt%atom_type2 == at2 .and. &
               dt%atom_type3 == at3 .and. &
               dt%atom_type4 == at4) .or. &
              (dt%atom_type1 == at4 .and. &
               dt%atom_type2 == at3 .and. &
               dt%atom_type3 == at2 .and. &
               dt%atom_type4 == at1)) then
            icn = icn+1
            icnt(icn) = i
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
      end if

      if (icn > MAX_MULTIPLY_NUM) &
        call error_msg('Search_dihe_param_duplicate> '// &
                       '[dihedraltypes] too much dihedral functions')

    end do

    if (icn == 0) then

      do i = 1, size(dihetypes)
        dt => dihetypes(i)

        if (dt%func /= dihes(old+iold)%func) &
          cycle

        if (dt%four_type .and. dt%wild_num == 1) then

          if (dt%atom_type1 == 'X' .and. &
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
      
          end if
        end if
        if (icn > MAX_MULTIPLY_NUM) &
          call error_msg('Search_dihe_param_duplicate> '// &
                         '[dihedraltypes] too much dihedral functions')
      end do

    end if

    if (icn == 0) then

      do i = 1, size(dihetypes)
        dt => dihetypes(i)

        if (dt%func /= dihes(old+iold)%func) &
          cycle

        if (dt%four_type .and. dt%wild_num == 2) then

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

          end if

        end if

        if (icn > MAX_MULTIPLY_NUM) &
          call error_msg('Search_dihe_param_duplicate> '// &
                         '[dihedraltypes] too much dihedral functions')

      end do

    end if

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
