!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_enefunc_str_mod
!> @brief   structure of energy functions
!! @authors Yuji Sugita (YS), Takashi Imai(TI), Chigusa Kobayashi (CK), 
!!          Takao Yoda (TY), Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_enefunc_str_mod

  use constants_mod
  use messages_mod
  use string_mod
#ifdef QSIMULATE
  use json_module
#endif

  implicit none
  private

  ! structures

  type, public :: s_pme_enefunc
    real(wp)                      :: alpha      ! alpha
    real(wp)                      :: alpha2m    ! -alpha^2
    real(wp)                      :: alpha2sp   ! 2*alpha/sqrt(pi)
    real(wp)                      :: u_self     ! Ewald self energy
    real(wp)                      :: bs_fact    ! B-spline factor (1/(n-1)...2)
    real(wp)                      :: bs_fact3   ! bs_fact^3
    real(wp)                      :: bs_fact3d  ! bs_fact^3*(n-1)
    integer                       :: n_bs       ! Order of B-spline
    integer                       :: ngrid(3)   ! Number of grid
    integer                       :: istart_atom
    integer                       :: iend_atom
    real(wp),         allocatable :: b2(:,:)         ! b^2(hx,hy,hz)
    real(wp),         allocatable :: gx(:)
    real(wp),         allocatable :: gy(:)
    real(wp),         allocatable :: gz(:)
    real(wp),         allocatable :: vir_fact(:,:,:) ! -2*(1+G^2/4a)/G^2
    real(wp),         allocatable :: theta(:,:,:)    ! F^-1[Theta](hx,hy,hz)
    real(wp),         allocatable :: f(:,:)
    real(wp),         allocatable :: v(:,:)
    real(wp),         allocatable :: qdf(:,:,:)
    real(wp),         allocatable :: force_wk(:,:)
    real(wp),         allocatable :: qdf_recv(:,:,:,:)
    real(wp),         allocatable :: qdf_send(:,:,:,:)
    integer,          allocatable :: i_list(:)
    integer,          allocatable :: x_start1(:)
    integer,          allocatable :: x_end1(:)
    integer,          allocatable :: y_start(:)
    integer,          allocatable :: y_end(:)
    integer,          allocatable :: y_local(:)
    integer,          allocatable :: y_start1(:)
    integer,          allocatable :: y_end1(:)
    integer,          allocatable :: z_start(:)
    integer,          allocatable :: z_end(:)
    integer,          allocatable :: z_local(:)
    complex(wp),      allocatable :: ftqdf_localA(:,:,:)
    complex(wp),      allocatable :: ftqdf_localB(:,:,:)
    real(wp)                      :: pme_max_spacing
  end type s_pme_enefunc

  type, public :: s_gbsa_enefunc
    integer                       :: num_sasa_atoms
    integer                       :: istart_gb
    integer                       :: iend_gb
    real(wp)                      :: cutoffdist
    real(wp)                      :: eps_solvent
    real(wp)                      :: eps_solute
    real(wp)                      :: sfactor_alpha
    real(wp)                      :: sfactor_beta
    real(wp)                      :: sfactor_gamma
    real(wp)                      :: salt_cons
    real(wp)                      :: vdw_offset
    real(wp)                      :: surface_tension
    real(wp)                      :: temperature
    real(wp)                      :: pairlist_distance
    logical                       :: tmp_calc
    real(wp),         allocatable :: born_radius(:)
    real(wp),         allocatable :: dadrfac(:)
    real(wp),         allocatable :: vdw_radius(:)
    real(wp),         allocatable :: scale_factor(:)
    real(wp),         allocatable :: sasa_parameters(:,:)
    real(wp),         allocatable :: sasa_vdw_radius(:)
    integer,          allocatable :: sasa_atom_type(:)
    integer,          allocatable :: sasa_atom_list(:)
  end type s_gbsa_enefunc

  type, public :: s_eef1_enefunc
    integer                       :: num_atoms
    integer                       :: istart
    integer                       :: iend
    character(6),     allocatable :: atom_name(:)
    real(wp)                      :: temperature
    real(wp)                      :: imm1_memb_thick
    real(wp)                      :: imm1_factor_a
    integer                       :: imm1_exponent_n
    logical                       :: imm1_make_pore
    real(wp)                      :: imm1_pore_radius
    real(wp)                      :: imic_axis_a
    real(wp)                      :: imic_axis_b
    real(wp)                      :: imic_axis_c
    real(wp)                      :: imic_exponent_m1
    real(wp)                      :: imic_exponent_m2
    real(wp)                      :: imic_steepness
    real(wp),         allocatable :: rvdw(:)
    real(wp),         allocatable :: volume(:,:)
    real(wp),         allocatable :: gref_0(:,:)
    real(wp),         allocatable :: gfree_0(:,:)
    real(wp),         allocatable :: href_0(:,:)
    real(wp),         allocatable :: cpref_0(:,:)
    real(wp),         allocatable :: gref_t(:,:)
    real(wp),         allocatable :: gfree_t(:,:)
    real(wp),         allocatable :: inv_lambda(:,:)
    real(wp),         allocatable :: alpha_4pi(:,:)
  end type s_eef1_enefunc

  type, public :: s_table_enefunc

    ! General table
    logical                       :: table
    integer                       :: table_order
    real(wp)                      :: density
    
    integer                       :: atom_cls_no_O
    integer                       :: atom_cls_no_H
    real(wp)                      :: charge_O
    real(wp)                      :: charge_H
    real(wp)                      :: mass_O
    real(wp)                      :: mass_H

    real(wp),         allocatable :: table_ene(:)
    real(wp),         allocatable :: table_grad(:)
    real(wp),         allocatable :: table_ecor(:)
    real(wp),         allocatable :: table_decor(:)
    real(wp),         allocatable :: dtable_ecor(:)
    real(wp),         allocatable :: dtable_decor(:)

    ! Water
    character(5)                  :: water_model
    integer                       :: num_water
    integer                       :: num_solute
    integer,          allocatable :: water_list(:,:)
    integer,          allocatable :: solute_list(:)
    integer,          allocatable :: solute_list_inv(:)

    real(wp),         allocatable :: table_ene_WW(:,:)
    real(wp),         allocatable :: table_de_WW(:,:)

    ! exclustion list
    !
    integer,          allocatable :: num_nonb_excl(:)
    integer,          allocatable :: nonb_excl_list(:,:)

    ! QM-EC list
    integer                   :: ecqm_num_nb14
    integer, allocatable      :: ecqm_nb14_list(:,:)

  end type s_table_enefunc

  type, public :: s_qmmm
    logical                   :: do_qmmm       = .false.
    character(MaxLine)        :: qmexe         = ''
    character(MaxLine)        :: qminp         = ''
    character(MaxLine)        :: qmout         = ''
    character(MaxLine)        :: qminfo        = ''
    character(MaxLine)        :: workdir       = ''
    character(MaxLine)        :: savedir       = ''
    character(MaxLine)        :: qmbase0       = ''
    character(MaxLine)        :: qmbasename    = ''
    character(MaxLine)        :: qmindex       = ''
    character(MaxLine)        :: qm_exemode    = ''
    integer                   :: qmsave_period = 0
    logical                   :: savefile      = .false.
    logical                   :: dryrun        = .false.
    logical                   :: save_qminfo   = .false.
    integer                   :: qmmaxtrial    = 0
    integer                   :: qmtrial       = 0
    integer                   :: qm_natoms     = 0
    integer                   :: ec_natoms     = 0
    integer, allocatable      :: qmatom_id(:)
    integer, allocatable      :: ecatom_id(:)
    real(wp)                  :: mm_cutoffdist = -1.0_wp
    logical                   :: mm_cutoff_byresidue = .true.
    integer                   :: num_qmmmbonds = 0
    integer, allocatable      :: qmmmbond_list(:,:)
    integer                   :: qm_nbonds     = 0
    integer                   :: qm_nangles    = 0
    integer                   :: qm_ndihedrals = 0
    integer                   :: qm_nimpropers = 0
    integer                   :: qm_ncmaps     = 0
    integer                   :: qm_count      = 0
    real(wp)                  :: qm_dipole(3)  = 0.0_wp
    logical                   :: ene_only      = .false.
    logical                   :: qm_classical  = .false.  ! Calculate ESP/MM forces
    logical                   :: qm_get_esp    = .false.  ! Read and store ESP charges
    real(wp), allocatable     :: qm_charge(:)
    real(wp), allocatable     :: qm_charge_save(:)
    logical                   :: qm_debug      = .false.
    real(wp)                  :: qm_total_charge = 0.0_wp ! For AMBER ff

    integer                   :: qm_nprocs     = 0
    integer                   :: inter_comm    = 0

    ! QMMM internal parameters -> at_qmmm.fpp
    integer                   :: qmtyp          = 1
    character(MaxLine)        :: qmcnt          = ''
    integer                   :: exclude_charge = 0
    integer,  allocatable     :: qm_atomic_no(:)
    real(wp), allocatable     :: qm_force(:,:)
    real(wp), allocatable     :: mm_force(:,:)
    integer                   :: mm_natoms      = 0
    integer,  allocatable     :: mmatom_id(:)
    integer                   :: mm_all_natoms  = 0
    integer                   :: mm_nres        = 0
    integer,  allocatable     :: mm_natoms_res(:)
    integer,  allocatable     :: mmatom_id_res(:,:)
    real(wp)                  :: linkbond
    real(wp), allocatable     :: linkatom_coord(:,:)
    real(wp), allocatable     :: linkatom_force(:,:)
    real(wp), allocatable     :: linkatom_charge(:)
    integer, allocatable      :: linkatom_global_address(:)
    logical                   :: is_qm_charge    = .false.
    logical                   :: ignore_qm_error = .false.
    logical                   :: qmmm_error      = .false.
    logical                   :: charges_bin     = .true.
    character(MaxLine)        :: tcscr   = ''
    character(MaxLine)        :: tcscr0  = ''
    logical                   :: defined_gen = .false.
#ifdef QSIMULATE
    type(json_file)           :: qs_input
    type(json_file)           :: qs_output
!    character(kind=c_char,len=:), allocatable :: bagel_output
#endif

  end type s_qmmm

  ! GaMD
  type, public :: s_enefunc_gamd
    ! gamd
    logical               :: gamd_stat
    logical               :: gamd_boost
    logical               :: boost_pot
    logical               :: boost_dih
    logical               :: boost_dual
    logical               :: thresh_lower
    logical               :: thresh_upper
    real(wp)              :: ene_pot_max
    real(wp)              :: ene_pot_min
    real(wp)              :: ene_pot_ave
    real(wp)              :: ene_pot_ave2
    real(wp)              :: ene_pot_dev
    real(wp)              :: ene_dih_max
    real(wp)              :: ene_dih_min
    real(wp)              :: ene_dih_ave
    real(wp)              :: ene_dih_ave2
    real(wp)              :: ene_dih_dev
    integer               :: count_pot
    integer               :: count_dih
    real(wp)              :: ene_pot_th
    real(wp)              :: ene_dih_th
    real(wp)              :: k0_pot
    real(wp)              :: k0_dih
    real(wp)              :: k_pot
    real(wp)              :: k_dih
    real(wp)              :: sigma0_pot
    real(wp)              :: sigma0_dih
    integer               :: update_period

    real(wp), allocatable :: f_dihe(:,:)
    real(wp), allocatable :: f_rest(:,:)
    real(wp)              :: v_dihe(3,3)
    real(wp)              :: v_rest(3,3)
  end type s_enefunc_gamd

  type, public :: s_enefunc

    integer                       :: forcefield
    integer                       :: output_style

    integer                       :: settle_func
    integer                       :: num_bonds
    integer                       :: num_angles
    integer                       :: num_ureys
    !shinobu-edited
    integer                       :: num_angflex
    !shinobu-edited
    integer                       :: num_dihedrals
    integer                       :: num_dihedflex
    integer                       :: num_rb_dihedrals
    integer                       :: num_impropers
    integer                       :: num_cmaps
    integer                       :: num_atom_cls
    integer                       :: num_contacts
    integer                       :: num_multi_contacts
    integer                       :: num_vsite2
    integer                       :: num_vsite3
    integer                       :: num_vsite3fd
    integer                       :: num_vsite3fad
    integer                       :: num_vsite3out
    integer                       :: num_vsite4fdn
    integer                       :: num_vsiten
    integer                       :: num_atoms_ref 
    integer                       :: num_restraintgroups
    integer                       :: num_restraintfuncs
    integer                       :: max_restraint_numatoms
    integer                       :: max_restraint_numgrps
    integer                       :: num_basins
    ! ~CG~ 3SPN.2C DNA: number of quartic bonds
    integer                       :: num_bonds_quartic
    ! ~CG~ 3SPN.2C DNA: number of base stacking
    integer                       :: num_base_stack
    logical                       :: cg_safe_dihedral_calc
    logical                       :: cg_DNA_base_pair_calc
    logical                       :: cg_DNA_exv_calc
    logical                       :: cg_ele_calc
    logical                       :: cg_pwmcos_calc
    logical                       :: cg_pwmcosns_calc
    logical                       :: cg_KH_calc
    logical                       :: cg_IDR_HPS_calc
    logical                       :: cg_IDR_KH_calc
    integer                       :: num_cg_particle_DNA_base
    integer                       :: num_cg_particle_DNA_phos
    integer                       :: num_cg_particle_DNA_all
    integer                       :: num_cg_particle_IDR_HPS
    integer                       :: num_cg_particle_IDR_KH
    integer                       :: num_cg_particle_KH
    integer                       :: num_cg_particle_charged
    integer                       :: num_cg_particle_pro_charged
    integer                       :: num_pwmcos_terms
    integer                       :: num_pwmcos_resid
    integer                       :: num_pwmcosns_terms
    integer                       :: num_pwmcosns_resid

    integer                       :: istart_bond,          iend_bond
    integer                       :: istart_angle,         iend_angle
    integer                       :: istart_urey,          iend_urey
    !shinobu-edited
    integer                       :: istart_angflex,       iend_angflex
    integer                       :: istart_dihedral,      iend_dihedral
    !shinobu-edited
    integer                       :: istart_dihedflex,     iend_dihedflex
    integer                       :: istart_rb_dihed,      iend_rb_dihed
    integer                       :: istart_improper,      iend_improper
    integer                       :: istart_cmap,          iend_cmap
    integer                       :: istart_contact,       iend_contact
    integer                       :: istart_vsite2,        iend_vsite2
    integer                       :: istart_vsite3,        iend_vsite3
    integer                       :: istart_vsite3fd,      iend_vsite3fd
    integer                       :: istart_vsite3fad,     iend_vsite3fad
    integer                       :: istart_vsite3out,     iend_vsite3out
    integer                       :: istart_vsite4fdn,     iend_vsite4fdn
    integer                       :: istart_vsiten,        iend_vsiten
    integer                       :: istart_restraint,     iend_restraint
    integer                       :: istart_multi_contact, iend_multi_contact
    integer                       :: istart_morph_bb,      iend_morph_bb
    integer                       :: istart_morph_sc,      iend_morph_sc
    ! ~CG~ 3SPN.2C DNA: quartic bonds
    integer                       :: istart_bond_quartic,  iend_bond_quartic
    ! ~CG~ 3SPN.2C DNA: base stacking
    integer                       :: istart_base_stack,    iend_base_stack

    ! bond (size = num_bonds)
    integer,          allocatable :: bond_list(:,:)
    real(wp),         allocatable :: bond_force_const(:)
    real(wp),         allocatable :: bond_dist_min(:)

    ! ~CG~ 3SPN.2C DNA : uncommon bond list for k(r-r0)^2 + 100 k(r-r0)^4
    integer,          allocatable :: bond_quartic_list(:,:)
    real(wp),         allocatable :: bond_quartic_force_const(:)
    real(wp),         allocatable :: bond_quartic_dist_min(:)

    ! angle (size = num_angles)
    integer,          allocatable :: angl_list(:,:)
    real(wp),         allocatable :: angl_force_const(:)
    real(wp),         allocatable :: angl_theta_min(:)
    ! ~CG~ AICG2+ angle
    integer,          allocatable :: angl_func(:)
    !CK
    logical                       :: multi_angle
    real(wp),         allocatable :: angl_force_const1(:)
    real(wp),         allocatable :: angl_theta_min1(:)
    real(wp),         allocatable :: angl_gamma(:)
    real(wp),         allocatable :: angl_epsa(:)
    !shinobu-edited (aicg2p local angl)
    real(wp),         allocatable :: angl_w(:)

    ! urey-bradley (size = num_ureys)
    integer,          allocatable :: urey_list(:,:)
    real(wp),         allocatable :: urey_force_const(:)
    real(wp),         allocatable :: urey_rmin(:)

    !shinobu-edited
    integer,          allocatable :: anglflex_list(:,:)
    integer,          allocatable :: anglflex_type(:)
    real(wp),         allocatable :: anglflex_theta(:,:)
    real(wp),         allocatable :: anglflex_efunc(:,:)
    real(wp),         allocatable :: anglflex_d2func(:,:)
    real(wp),         allocatable :: anglflex_min_th(:,:)
    real(wp),         allocatable :: anglflex_max_th(:,:)
    real(wp),         allocatable :: anglflex_ener_corr(:)

    ! dihedral (size = num_dihedrals)
    integer,          allocatable :: dihe_list(:,:)
    real(wp),         allocatable :: dihe_force_const(:)
    integer,          allocatable :: dihe_periodicity(:)
    real(wp),         allocatable :: dihe_phase(:)
    real(wp),         allocatable :: dihe_scee(:)
    real(wp),         allocatable :: dihe_scnb(:)
    !shinobu-edited
    real(wp),         allocatable :: dihe_theta_min(:)
    real(wp),         allocatable :: dihe_w(:)
    ! ~CG~ 3SPN.2C DNA: dihedral function types
    integer,          allocatable :: dihe_func(:)

    !shinobu-edited
    integer,          allocatable :: diheflex_list(:,:)
    integer,          allocatable :: diheflex_type(:)
    integer,          allocatable :: diheflex_func(:)
    real(wp),         allocatable :: diheflex_coef(:,:)
    real(wp),         allocatable :: diheflex_ener_corr(:)
    real(wp)                      :: cg_safe_dih_ene_shift

    ! dihedral (size = num_rb_dihedrals)
    integer,          allocatable :: rb_dihe_list(:,:)
    real(wp),         allocatable :: rb_dihe_c(:,:)

    ! improper (size = num_impropers)
    integer,          allocatable :: impr_list(:,:)
    real(wp),         allocatable :: impr_force_const(:)
    integer,          allocatable :: impr_periodicity(:)
    real(wp),         allocatable :: impr_phase(:)

    ! cmap (size = num_cmaps)
    !
    !   cmap_list(i,n) : atom list for cmap where (1 <= i <= 8 and
    !                        1 <= n <= psf%num_cross_terms)
    !
    !   cmap_resolution(n) : the number of cells in one direction
    !                        1 <= n <= psf%num_cross_terms)
    !
    !   cmap_coef(i,j,idih1,idih2,n) : the c_ij value for n-th cmap.
    !                    where (1 <= n <= psf%num_cross_terms)
    !                      ,(1 <= i, j <= 4) and 
    !                       (1 <= idih1, idih2 <= 24).
    integer,          allocatable :: cmap_list(:,:)
    integer,          allocatable :: cmap_resolution(:)
    integer,          allocatable :: cmap_type(:)
    real(wp),         allocatable :: cmap_coef(:,:,:,:,:)
    real(wp),         allocatable :: cmap_force(:,:,:)

    ! base stacking (size = num_base_stack)
    ! ~CG~ 3SPN.2C DNA: Base Stacking Data structures
    integer,          allocatable :: base_stack_list(:,:)
    integer,          allocatable :: base_stack_func(:)
    real(wp),         allocatable :: base_stack_epsilon(:)
    real(wp),         allocatable :: base_stack_sigma(:)
    real(wp),         allocatable :: base_stack_theta_bs(:)
    real(wp)                      :: base_stack_alpha
    real(wp)                      :: base_stack_K
    
    ! base pairing (size = num_base_type)
    ! ~CG~ 3SPN.2C DNA: Base Pairing Data structures
    real(wp),         allocatable :: base_pair_theta_1(:)
    real(wp),         allocatable :: base_pair_theta_2(:)
    real(wp),         allocatable :: base_pair_theta_3(:)
    real(wp),         allocatable :: base_pair_phi_1(:)
    real(wp),         allocatable :: base_pair_sigma(:)
    real(wp),         allocatable :: base_pair_epsilon(:)
    real(wp)                      :: base_pair_alpha
    real(wp)                      :: base_pair_K
    ! real(wp)                      :: base_pair_sigma_AT
    ! real(wp)                      :: base_pair_sigma_GC
    ! real(wp)                      :: base_pair_epsilon_AT
    ! real(wp)                      :: base_pair_epsilon_GC
    ! ~CG~ 3SPN.2C DNA: Base Crossing Data structures
    real(wp),         allocatable :: base_cross_1_epsilon(:,:)
    real(wp),         allocatable :: base_cross_1_sigma(:,:)
    real(wp),         allocatable :: base_cross_1_theta_cs(:,:)
    real(wp),         allocatable :: base_cross_2_epsilon(:,:)
    real(wp),         allocatable :: base_cross_2_sigma(:,:)
    real(wp),         allocatable :: base_cross_2_theta_cs(:,:)
    logical,          allocatable :: base_pair_is_WC(:,:)
    real(wp)                      :: base_cross_alpha
    real(wp)                      :: base_cross_K
    integer,          allocatable :: NA_base_type(:)
    integer,          allocatable :: mol_chain_id(:)
    ! integer,          allocatable :: mol_atom_id(:)
    integer,          allocatable :: atom_cls(:)

    ! ~CG~ 3spn.2c DNA exv (size = num_DNA_particle_types)
    real(wp),         allocatable :: cgDNA_exv_sigma(:,:)
    real(wp)                      :: cgDNA_exv_epsilon
    logical                       :: cg_infinite_DNA

    ! ~CG~ : debye-huckel ele
    real(wp),         allocatable :: cg_charge(:)
    real(wp)                      :: cg_pro_DNA_ele_scale_Q
    integer,          allocatable :: cg_ele_mol_pair(:,:)

    ! ~CG~ : KH model
    integer,          allocatable :: cg_KH_mol_pair(:,:)
    logical,          allocatable :: cg_pro_use_KH(:)
    real(wp),         allocatable :: cg_KH_sigma_half(:)
    real(wp),         allocatable :: cg_KH_epsilon(:, :)

    ! ~CG~ : IDR HPS model
    logical,          allocatable :: cg_IDR_HPS_is_IDR(:)
    real(wp),         allocatable :: cg_IDR_HPS_lambda_half(:)
    real(wp),         allocatable :: cg_IDR_HPS_sigma_half(:)
    real(wp)                      :: cg_IDR_HPS_epsilon

    ! ~CG~ : IDR KH model
    logical,          allocatable :: cg_IDR_KH_is_IDR(:)
    real(wp),         allocatable :: cg_IDR_KH_sigma_half(:)
    real(wp),         allocatable :: cg_IDR_KH_epsilon_D(:, :)

    ! ~CG~ : KH-MJ model
    real(wp)                      :: cg_KH_mod_A_lambda = 0.159
    real(wp)                      :: cg_KH_mod_B_lambda = 0.186
    real(wp)                      :: cg_KH_mod_C_lambda = 0.192
    real(wp)                      :: cg_KH_mod_D_lambda = 0.228
    real(wp)                      :: cg_KH_mod_E_lambda = 0.194
    real(wp)                      :: cg_KH_mod_F_lambda = 0.223
    real(wp)                      :: cg_KH_mod_A_eps_0  = -2.27
    real(wp)                      :: cg_KH_mod_B_eps_0  = -1.95
    real(wp)                      :: cg_KH_mod_C_eps_0  = -1.85
    real(wp)                      :: cg_KH_mod_D_eps_0  = -1.67
    real(wp)                      :: cg_KH_mod_E_eps_0  = -2.00
    real(wp)                      :: cg_KH_mod_F_eps_0  = -1.96
    
    ! ~CG~ PWMcos: CG protein-DNA seq-specific
    real(wp)                      :: pwmcos_sigma
    real(wp)                      :: pwmcos_phi
    integer,          allocatable :: pwmcos_protein_id(:)
    integer,          allocatable :: pwmcos_protein_id_N(:)
    integer,          allocatable :: pwmcos_protein_id_C(:)
    real(wp),         allocatable :: pwmcos_r0(:)
    real(wp),         allocatable :: pwmcos_theta1(:)
    real(wp),         allocatable :: pwmcos_theta2(:)
    real(wp),         allocatable :: pwmcos_theta3(:)
    real(wp),         allocatable :: pwmcos_ene_A(:)
    real(wp),         allocatable :: pwmcos_ene_C(:)
    real(wp),         allocatable :: pwmcos_ene_G(:)
    real(wp),         allocatable :: pwmcos_ene_T(:)
    real(wp),         allocatable :: pwmcos_gamma(:)
    real(wp),         allocatable :: pwmcos_eps(:)
    integer,          allocatable :: pwmcos_specificity(:)
    integer,          allocatable :: pwmcos_mol_pair(:,:)
    integer,          allocatable :: pwmcos_to_pairlist_id(:)
    integer,          allocatable :: pwmcos_involved_resid(:)
    integer,          allocatable :: pwmcos_involved_spec(:)

    ! ~CG~ PWMcosns: CG protein-DNA seq-nonspecific
    real(wp)                      :: pwmcosns_sigma
    real(wp)                      :: pwmcosns_phi
    integer,          allocatable :: pwmcosns_protein_id(:)
    integer,          allocatable :: pwmcosns_protein_id_N(:)
    integer,          allocatable :: pwmcosns_protein_id_C(:)
    real(wp),         allocatable :: pwmcosns_r0(:)
    real(wp),         allocatable :: pwmcosns_theta1(:)
    real(wp),         allocatable :: pwmcosns_theta2(:)
    real(wp),         allocatable :: pwmcosns_ene(:)
    integer,          allocatable :: pwmcosns_specificity(:)
    integer,          allocatable :: pwmcosns_mol_pair(:,:)
    integer,          allocatable :: pwmcosns_to_pairlist_id(:)
    integer,          allocatable :: pwmcosns_involved_resid(:)
    integer,          allocatable :: pwmcosns_involved_spec(:)

    ! ~CG~ : list of CG particles with special properties
    integer,          allocatable :: cg_particle_DNA_base(:)
    integer,          allocatable :: cg_particle_DNA_phos(:)
    integer,          allocatable :: cg_particle_DNA_all(:)
    integer,          allocatable :: cg_particle_IDR_HPS(:)
    integer,          allocatable :: cg_particle_IDR_KH(:)
    integer,          allocatable :: cg_particle_KH(:)
    integer,          allocatable :: cg_particle_charged(:)
    integer,          allocatable :: cg_particle_pro_charged(:)
    integer,          allocatable :: cg_istart_nonb_excl(:)

    ! nonbonded 1-4 interactions (size = num_nb14)

    ! non-bonded (size = num_atom_cls)
    integer,          allocatable :: nonb_atom_cls(:)
    real(wp),         allocatable :: nb14_lj6(:,:)
    real(wp),         allocatable :: nb14_lj10(:,:)
    real(wp),         allocatable :: nb14_lj12(:,:)
    real(wp),         allocatable :: nonb_lj6(:,:)
    real(wp),         allocatable :: nonb_lj10(:,:)
    real(wp),         allocatable :: nonb_lj12(:,:)
    real(wp),         allocatable :: nonb_aicg_eps(:,:)
    real(wp),         allocatable :: nonb_aicg_sig(:,:)
    real(wp),         allocatable :: cg_exv_eps_sqrt(:)
    real(wp),         allocatable :: cg_exv_sig_half(:)
    real(wp),         allocatable :: param_epsilon(:)
    real(wp),         allocatable :: param_sigma(:)
    ! for GO (size = num_atom)
    real(wp),         allocatable :: nonb_eps(:)
    real(wp),         allocatable :: nonb_rmin(:)

    ! non-bonded (size = num_atoms)
    integer,          allocatable :: num_nonb_excl(:)
    integer,          allocatable :: num_nb14_calc(:)
    integer,          allocatable :: nonb_excl_list(:)
    integer,          allocatable :: nb14_calc_list(:,:)
    real(wp),         allocatable :: nb14_qq_scale(:,:)
    real(wp),         allocatable :: nb14_lj_scale(:,:)
    real(wp),         allocatable :: work(:,:)

    real(wp)                      :: nb14_qq_scale_c19

    ! native contact (size = num_contacts)
    integer,          allocatable :: contact_list(:,:)
    integer,          allocatable :: contact_func(:)
    real(wp),         allocatable :: contact_lj6(:)
    real(wp),         allocatable :: contact_lj10(:)
    real(wp),         allocatable :: contact_lj12(:)
    !CK
    integer,          allocatable :: multi_contact_list(:,:)
    integer,          allocatable :: multi_contact_model(:)
    real(wp),         allocatable :: multi_contact_lj6(:)
    real(wp),         allocatable :: multi_contact_lj10(:)
    real(wp),         allocatable :: multi_contact_lj12(:)

    ! virtual site 2
    integer,          allocatable :: vsite2_list(:,:)
    real(wp),         allocatable :: vsite2_a(:)
    ! virtual site 3
    integer,          allocatable :: vsite3_list(:,:)

    real(wp),         allocatable :: vsite3_a(:)
    real(wp),         allocatable :: vsite3_b(:)
    ! virtual site 3fd
    integer,          allocatable :: vsite3fd_list(:,:)
    real(wp),         allocatable :: vsite3fd_a(:)
    real(wp),         allocatable :: vsite3fd_d(:)
    ! virtual site 3fad
    integer,          allocatable :: vsite3fad_list(:,:)
    real(wp),         allocatable :: vsite3fad_theta(:)
    real(wp),         allocatable :: vsite3fad_d(:)
    ! virtual site 3out
    integer,          allocatable :: vsite3out_list(:,:)
    real(wp),         allocatable :: vsite3out_a(:)
    real(wp),         allocatable :: vsite3out_b(:)
    real(wp),         allocatable :: vsite3out_c(:)
    ! virtual site 4fdn
    integer,          allocatable :: vsite4fdn_list(:,:)
    real(wp),         allocatable :: vsite4fdn_a(:)
    real(wp),         allocatable :: vsite4fdn_b(:)
    real(wp),         allocatable :: vsite4fdn_c(:)
    ! virtual site n
    integer,          allocatable :: vsiten_list(:,:)
    integer,          allocatable :: vsiten_n(:)

    ! non contact pairs
    real(wp)                      :: noncontact_lj12
    real(wp)                      :: mix_temperature
    real(wp)                      :: mix_beta

    real(wp)                      :: switchdist
    real(wp)                      :: cutoffdist
    real(wp)                      :: pairlistdist
    real(wp)                      :: dielec_const
    real(wp)                      :: debye
    real(wp)                      :: cg_cutoffdist_ele
    real(wp)                      :: cg_cutoffdist_126
    real(wp)                      :: cg_cutoffdist_DNAbp
    real(wp)                      :: cg_pairlistdist_ele
    real(wp)                      :: cg_pairlistdist_126
    real(wp)                      :: cg_pairlistdist_PWMcos
    real(wp)                      :: cg_pairlistdist_DNAbp
    real(wp)                      :: cg_pairlistdist_exv
    real(wp)                      :: cg_ele_coef
    real(wp)                      :: cg_ele_sol_T
    real(wp)                      :: cg_ele_sol_IC
    real(wp)                      :: cg_dielec_const
    real(wp)                      :: cg_debye_length

    ! flag for position restraint 
    logical                       :: restraint_posi
    ! flag for restraint calculation
    logical                       :: restraint_flag
    integer                       :: target_function
    integer                       :: steered_function
    real(wp)                      :: target_value
    real(wp)                      :: rmsd_force
    ! restraint group (size = num_restraintgroups)
    integer,          allocatable :: restraint_numatoms(:)
    integer,          allocatable :: restraint_atomlist(:,:)
    real(wp),         allocatable :: restraint_totmass_group(:)
    real(wp),         allocatable :: restraint_masscoef(:,:)
    real(wp),         allocatable :: restraint_masstmp(:)
    integer,          allocatable :: restraint_istartend(:,:)

    ! restraint func (size = num_restraintfuncs)
    integer,          allocatable :: restraint_kind(:)
    integer,          allocatable :: restraint_grouplist(:,:)
    integer,          allocatable :: restraint_funcgrp(:)
    integer,          allocatable :: restraint_exponent_func(:)
    integer,          allocatable :: restraint_exponent_dist(:,:)
    integer,          allocatable :: restraint_mode(:)
    integer,          allocatable :: restraint_rpath_func(:)
    integer,          allocatable :: restraint_replica_index(:)
    logical,          allocatable :: restraint_caging(:)
    real(wp),         allocatable :: restraint_flat_radius(:)
    real(wp),         allocatable :: restraint_weight_dist(:,:)
    real(wp),         allocatable :: restraint_const(:,:)
    real(wp),         allocatable :: restraint_ref(:,:)
    real(wp),         allocatable :: restraint_wcom1(:,:)
    real(wp),         allocatable :: restraint_wcom2(:,:)
    real(wp),         allocatable :: restraint_wcom3(:,:)
    real(wp),         allocatable :: restraint_wcom4(:,:)
    real(wp),         allocatable :: restraint_wcom5(:,:)
    real(wp),         allocatable :: restraint_wtmp(:,:)
    real(wp),         allocatable :: restraint_wdrt(:)
    real(wp),         allocatable :: restraint_wfrc(:,:)
    real(wp),         allocatable :: restraint_rgbuffer(:)
    real(wp),         allocatable :: restraint_wall_z(:,:)
    ! for repul & fb
    real(wp),         allocatable :: restraint_rcom1(:,:)
    real(wp),         allocatable :: restraint_rcom2(:,:)
    real(wp),         allocatable :: restraint_rdrt(:)
 
    ! restraint func (size = num_atoms_ref)
    real(wp),         allocatable :: restraint_refcoord(:,:)

    ! restraint func (size = num_restraintfunct x ndata)
    real(wp),         allocatable :: restraint_const_replica(:,:)
    real(wp),         allocatable :: restraint_ref_replica(:,:)

    ! multibasin work
    real(wp),         allocatable :: force_mb_work(:,:,:)
    real(wp),         allocatable :: energy_mb_work(:,:)
    real(wp),         allocatable :: virial_mb_work(:,:,:)
    real(wp),         allocatable :: basinenergy(:)

    ! membrane potential
    logical                       :: ez_membrane_flag
    real(wp),         allocatable :: ez_membrane_const(:)
    real(wp),         allocatable :: ez_membrane_zmin(:)
    real(wp),         allocatable :: ez_membrane_polym(:)
    real(wp),         allocatable :: ez_membrane_sigma(:)
    integer,          allocatable :: ez_membrane_func(:)
    ! principal component mode
    integer                       :: num_pc_modes
    real(wp),         allocatable :: pc_mode(:)
    real(wp),         allocatable :: pc_crd(:)
    real(wp),         allocatable :: pc_crd_ref(:)
    real(wp),         allocatable :: pc_crd_ref_fit(:)
    real(wp),         allocatable :: pc_mode_temp(:)
    real(wp),         allocatable :: pc_mode_fit(:)
    real(wp),         allocatable :: grad_pc(:,:)
     
    type(s_table_enefunc)         :: table
    logical                       :: force_switch
    logical                       :: vdw_shift

    integer                       :: nonb_func
    real(wp)                      :: fudge_lj
    real(wp)                      :: fudge_qq
    integer                       :: excl_level

    integer                       :: dispersion_corr
    real(wp)                      :: eswitch
    real(wp)                      :: vswitch
    real(wp)                      :: dispersion_energy
    real(wp)                      :: dispersion_virial

    ! pme
    logical                       :: pme_use
    type(s_pme_enefunc)           :: pme

    ! eef1/imm1/gbsa
    logical                       :: eef1_use
    logical                       :: imm1_use
    logical                       :: imic_use
    logical                       :: gbsa_use
    type(s_eef1_enefunc)          :: eef1
    type(s_gbsa_enefunc)          :: gbsa

    ! QM/MM
    type(s_qmmm)                  :: qmmm

    ! spherical boundary condition
    logical                       :: spot_use

    ! statistical variables
    logical                       :: rpath_flag
    logical                       :: rpath_sum_mf_flag
    integer                       :: rpath_pos_func
    integer                       :: stats_count
    integer                       :: stats_natom
    integer                       :: stats_dimension
    real(wp),         allocatable :: stats_delta(:)
    real(wp),         allocatable :: stats_grad(:, :, :)
    real(wp),         allocatable :: stats_force(:)
    real(wp),         allocatable :: stats_metric(:,:)
    integer,          allocatable :: stats_atom(:,:)
    real(wp),         allocatable :: stats_mass(:,:)

    real(wp),         allocatable :: rotated_coord(:,:)

    ! fit
    integer                       :: num_fitting
    integer                       :: fitting_method
    logical                       :: do_fitting
    integer                       :: fitting_move
    integer                       :: fitting_file
    logical                       :: mass_weight
    integer,          allocatable :: fitting_atom(:)
    real(wp),         allocatable :: fit_mass(:)
    real(wp),         allocatable :: fit_refcoord(:,:)
    real(wp),         allocatable :: fit_work(:,:)

    ! go
    integer                       :: go_electrostatic 

    ! morphing
    logical                       :: morph_flag = .false.
    logical                       :: morph_restraint_flag = .false.
    integer                       :: morph_ene_flag 
    integer                       :: num_morph_bb
    integer                       :: num_morph_sc
    integer,allocatable           :: morph_list_bb(:,:)
    integer,allocatable           :: morph_list_sc(:,:)
    real(wp),allocatable          :: morph_dist_bb(:)
    real(wp),allocatable          :: morph_dist_sc(:)
    real(wp),allocatable          :: morph_dist_bb_other(:)
    real(wp),allocatable          :: morph_dist_sc_other(:)
    real(wp)                      :: morph_spring = 5.0d0
    real(wp)                      :: morph_spring_max
    real(wp)                      :: morph_linear = 1.0d0

    logical                       :: contact_check
    logical                       :: nonb_limiter
    real(wp)                      :: minimum_contact
    logical                       :: pressure_position
    logical                       :: pressure_rmsd

    ! gamd
    logical                       :: gamd_use
    type(s_enefunc_gamd)          :: gamd

    ! mepmd
    logical                       :: rpath_sum_avforce
    integer                       :: mep_natoms
    integer,          pointer     :: mepatom_id(:)
    real(wp),         allocatable :: mepmd_avforce(:,:)


  end type s_enefunc

  ! parameter for allocatable variables
  integer,      public, parameter :: EneFuncBond          = 1
  integer,      public, parameter :: EneFuncAngl          = 2
  integer,      public, parameter :: EneFuncUrey          = 3
  integer,      public, parameter :: EneFuncDihe          = 4
  integer,      public, parameter :: EneFuncRBDihe        = 5
  integer,      public, parameter :: EneFuncImpr          = 6
  integer,      public, parameter :: EneFuncCmap          = 7 
  integer,      public, parameter :: EneFuncNbon          = 8 
  integer,      public, parameter :: EneFuncCntc          = 9
  integer,      public, parameter :: EneFuncRefg          = 10
  integer,      public, parameter :: EneFuncReff          = 11
  integer,      public, parameter :: EneFuncRefw          = 12
  integer,      public, parameter :: EneFuncRefc          = 13
  integer,      public, parameter :: EneFuncRefr          = 14
  integer,      public, parameter :: EneFuncTable         = 15
  integer,      public, parameter :: EneFuncTableWater    = 16
  integer,      public, parameter :: EneFuncTableSolute   = 17
  integer,      public, parameter :: EneFuncCmapType      = 18
  integer,      public, parameter :: EneFuncPme           = 19
  integer,      public, parameter :: EneFuncMode          = 20
  integer,      public, parameter :: EneFuncFitc          = 21
  integer,      public, parameter :: EneFuncQmmm          = 22
  integer,      public, parameter :: EneFuncModeRef       = 23
  integer,      public, parameter :: EneFuncEef1          = 24
  integer,      public, parameter :: EneFuncGbsa          = 25
  integer,      public, parameter :: EneFuncVsite2        = 26
  integer,      public, parameter :: EneFuncVsite3        = 27
  integer,      public, parameter :: EneFuncVsite3fd      = 28
  integer,      public, parameter :: EneFuncVsite3fad     = 29
  integer,      public, parameter :: EneFuncVsite3out     = 30
  integer,      public, parameter :: EneFuncVsite4fdn     = 31
  integer,      public, parameter :: EneFuncVsiten        = 32
  integer,      public, parameter :: EneFuncEzMembrane    = 33
  integer,      public, parameter :: EneFuncNonbGO        = 34
  integer,      public, parameter :: EneFuncMultiCntc     = 35
  integer,      public, parameter :: EneFuncMultiWork     = 36
  integer,      public, parameter :: EneFuncMorph         = 37
  !shinobu-edited
  integer,      public, parameter :: EneFuncAngFlex       = 38
  integer,      public, parameter :: EneFuncDiheFlex      = 39
  integer,      public, parameter :: EneFuncAngFlexTbl    = 40
  integer,      public, parameter :: EneFuncDiheFlexTbl   = 41
  ! ~CG~ 
  integer,      public, parameter :: EneFuncCGGeneral     = 42
  integer,      public, parameter :: EneFuncBondQuartic   = 43
  integer,      public, parameter :: EneFuncBaseStack     = 44
  integer,      public, parameter :: EneFuncBasePair      = 45
  integer,      public, parameter :: EneFuncCGDNAExv      = 46
  integer,      public, parameter :: EneFuncCGele         = 47
  integer,      public, parameter :: EneFuncPWMcos        = 48

  integer,      public, parameter :: EneFuncGamdDih       = 49
  integer,      public, parameter :: EneFuncGamdRest      = 50

  integer,      public, parameter :: EneFuncPWMcosns      = 51
  integer,      public, parameter :: EneFuncCGIDRHPS      = 52
  integer,      public, parameter :: EneFuncCGIDRKH       = 53
  integer,      public, parameter :: EneFuncCGKH          = 54
  integer,      public, parameter :: EneFuncCGKHmol       = 55

  ! parameters
  integer,      public, parameter :: ForcefieldCHARMM     = 1
  integer,      public, parameter :: ForcefieldCHARMM19   = 2
  integer,      public, parameter :: ForcefieldAAGO       = 3
  integer,      public, parameter :: ForcefieldCAGO       = 4
  integer,      public, parameter :: ForcefieldKBGO       = 5
  integer,      public, parameter :: ForcefieldAMBER      = 6
  integer,      public, parameter :: ForcefieldGROAMBER   = 7
  integer,      public, parameter :: ForcefieldGROMARTINI = 8
  integer,      public, parameter :: ForcefieldSOFT       = 9
  !shinobu-edited
  integer,      public, parameter :: ForcefieldRESIDCG     = 10

  integer,      public, parameter :: OutputStyleGENESIS   = 1
  integer,      public, parameter :: OutputStyleCHARMM    = 2
  integer,      public, parameter :: OutputStyleNAMD      = 3
  integer,      public, parameter :: OutputStyleGROMACS   = 4
  integer,      public, parameter :: MBeneNonb            = 1
  integer,      public, parameter :: MBeneCntc            = 2
  integer,      public, parameter :: MBeneTotal           = 3
  integer,      public, parameter :: MBeneExp             = 4 ! should be last

  integer,      public, parameter :: MorphEneFlagTot      = 1
  integer,      public, parameter :: MorphEneFlagBB       = 2
  integer,      public, parameter :: MorphEneFlagSC       = 3

  integer,      public, parameter :: EzMembraneNo         = 0
  integer,      public, parameter :: EzMembraneSigmoid    = 1
  integer,      public, parameter :: EzMembraneGaussian   = 2

  ! BaseType for bases are dynamicvally determined from .top files.
  ! 
  integer,      public, parameter :: NABaseTypeDBA        = 1
  integer,      public, parameter :: NABaseTypeDBC        = 2
  integer,      public, parameter :: NABaseTypeDBG        = 3
  integer,      public, parameter :: NABaseTypeDBT        = 4
  integer,      public, parameter :: NABaseTypeDBMAX      = 4
  integer,      public, parameter :: NABaseTypeRBA        = 5
  integer,      public, parameter :: NABaseTypeRBC        = 6
  integer,      public, parameter :: NABaseTypeRBG        = 7
  integer,      public, parameter :: NABaseTypeRBU        = 8
  integer,      public, parameter :: NABaseTypeRBMAX      = 8
  ! 
  ! max value for DNA/RNA Base type
  ! 
  integer,      public, parameter :: NABaseTypeBMAX       = 10 
  integer,      public, parameter :: NABaseTypeDP         = 11
  integer,      public, parameter :: NABaseTypeDS         = 12
  integer,      public, parameter :: NABaseTypeRP         = 13
  integer,      public, parameter :: NABaseTypeRS         = 14
  integer,      public, parameter :: NABaseTypeNAMAX      = 15
  integer,      public, parameter :: NABaseTypeProtein    = 21
  
  character(*), public, parameter :: ForceFieldTypes(10) = (/'CHARMM    ', &
                                                             'CHARMM19  ', &
                                                             'AAGO      ', &
                                                             'CAGO      ', &
                                                             'KBGO      ', &
                                                             'AMBER     ', &
                                                             'GROAMBER  ', &
                                                             'GROMARTINI', &
                                                             'SOFT      ', &
                                                             'RESIDCG   '/)
!shinobu-edited

  character(*), public, parameter :: OutputStyleTypes(4) = (/'GENESIS', &
                                                             'CHARMM ', &
                                                             'NAMD   ', &
                                                             'GROMACS'/)

  ! parameters (Dispersion Correction)
  integer,      public, parameter :: Disp_corr_NONE       = 1
  integer,      public, parameter :: Disp_corr_Energy     = 2
  integer,      public, parameter :: Disp_corr_EPress     = 3
  character(*), public, parameter :: Disp_corr_Types(3)   = (/'NONE  ', &
                                                              'ENERGY', &
                                                              'EPRESS'/)

  ! parameters (GO model)
  integer,      public, parameter :: GoElectrostaticNONE  = 1
  integer,      public, parameter :: GoElectrostaticALL   = 2 
  integer,      public, parameter :: GoElectrostaticINTER = 3 
  character(*), public, parameter :: GoElectrostaticTypes(3) = (/'NONE  ', &
                                                                 'ALL   ', &
                                                                 'INTER '/)

  ! parameters (qmtyp for QMMM)
  integer,      public, parameter :: QMtypQCHEM     = 1
  integer,      public, parameter :: QMtypG09       = 2
  integer,      public, parameter :: QMtypTERACHEM  = 3
  integer,      public, parameter :: QMtypDFTBPLUS  = 4
  integer,      public, parameter :: QMtypMOLPRO    = 5
  integer,      public, parameter :: QMtypG09_FRWF  = 6
  integer,      public, parameter :: QMtypQSimulate = 7
  integer,      public, parameter :: QMtypORCA      = 8
  character(*), public, parameter :: QMtypTypes(8)  = (/'qchem        ',  &
                                                        'gaussian     ',  &
                                                        'terachem     ',  &
                                                        'dftb+        ',  &
                                                        'molpro       ',  &
                                                        'gaussian_frwf',  &
                                                        'qsimulate    ',  &
                                                        'orca         '/)

  ! parameters (Exclude charge for QMMM)
  integer,      public, parameter :: ExcludeChargeAtom    = 1
  integer,      public, parameter :: ExcludeChargeGroup   = 2
  integer,      public, parameter :: ExcludeChargeAtomAndDistribute   = 3
  character(*), public, parameter :: ExcludeChargeTypes(3)  = (/'ATOM ', &
                                                                'GROUP', &
                                                                'AMBER'/)

  ! parameters (linlbond for QMMM)
  real(wp),     public, parameter :: LinkBondCHARMM = 1.0_wp
  real(wp),     public, parameter :: LinkBondAMBER = 1.09_wp

  ! The length of command for system call
  integer,      public, parameter :: len_exec = 1500

  ! subroutines
  public  :: init_enefunc
  public  :: alloc_enefunc
  public  :: dealloc_enefunc
  public  :: dealloc_enefunc_all

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_enefunc
  !> @brief        initialize energy functions information
  !! @authors      YS, CK
  !! @param[out]   enefunc  : structure of potential energy function
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_enefunc(enefunc)

    ! formal arguments
    type(s_enefunc),         intent(inout) :: enefunc

    
    enefunc%num_bonds               = 0 
    enefunc%num_angles              = 0
    enefunc%num_ureys               = 0
    enefunc%num_dihedrals           = 0 
    enefunc%num_rb_dihedrals        = 0
    enefunc%num_impropers           = 0 
    enefunc%num_cmaps               = 0
    enefunc%num_atom_cls            = 0
    enefunc%num_contacts            = 0
    enefunc%num_multi_contacts      = 0
    enefunc%num_vsite2              = 0
    enefunc%num_vsite3              = 0
    enefunc%num_vsite3fd            = 0
    enefunc%num_vsite3fad           = 0
    enefunc%num_vsite3out           = 0
    enefunc%num_vsite4fdn           = 0
    enefunc%num_vsiten              = 0
    enefunc%num_atoms_ref           = 0
    enefunc%num_restraintgroups     = 0 
    enefunc%num_restraintfuncs      = 0 
    enefunc%max_restraint_numatoms  = 0 
    enefunc%max_restraint_numgrps   = 0 
    enefunc%num_bonds_quartic       = 0
    enefunc%num_base_stack          = 0
    enefunc%num_pwmcos_terms        = 0
    enefunc%num_pwmcos_resid        = 0
    enefunc%num_pwmcosns_terms      = 0
    enefunc%num_pwmcosns_resid      = 0
    
    enefunc%istart_bond             = 0
    enefunc%iend_bond               = 0
    enefunc%istart_bond_quartic     = 0
    enefunc%iend_bond_quartic       = 0
    enefunc%istart_angle            = 0
    enefunc%iend_angle              = 0
    enefunc%istart_urey             = 0
    enefunc%iend_urey               = 0
    enefunc%istart_dihedral         = 0
    enefunc%iend_dihedral           = 0
    enefunc%istart_rb_dihed         = 0
    enefunc%iend_rb_dihed           = 0
    enefunc%istart_improper         = 0
    enefunc%iend_improper           = 0
    enefunc%istart_cmap             = 0
    enefunc%iend_cmap               = 0
    enefunc%istart_base_stack       = 0
    enefunc%iend_base_stack         = 0
    enefunc%istart_contact          = 0
    enefunc%iend_contact            = 0
    enefunc%istart_vsite2           = 0
    enefunc%iend_vsite2             = 0
    enefunc%istart_vsite3           = 0
    enefunc%iend_vsite3             = 0
    enefunc%istart_vsite3fd         = 0
    enefunc%iend_vsite3fd           = 0
    enefunc%istart_vsite3fad        = 0
    enefunc%iend_vsite3fad          = 0
    enefunc%istart_vsite3out        = 0
    enefunc%iend_vsite3out          = 0
    enefunc%istart_vsite4fdn        = 0
    enefunc%iend_vsite4fdn          = 0
    enefunc%istart_vsiten           = 0
    enefunc%iend_vsiten             = 0
    enefunc%istart_restraint        = 0
    enefunc%iend_restraint          = 0

    enefunc%forcefield              = ForcefieldCHARMM
    enefunc%output_style            = OutputStyleCHARMM

    enefunc%noncontact_lj12         = 0.0_wp

    enefunc%switchdist              = 0.0_wp
    enefunc%cutoffdist              = 0.0_wp
    enefunc%pairlistdist            = 0.0_wp
    enefunc%dielec_const            = 0.0_wp
    enefunc%debye                   = 0.0_wp

    enefunc%pme_use                 = .false.

    enefunc%target_function         = 0
    enefunc%steered_function        = 0
    enefunc%rmsd_force              = 0.0_wp
    enefunc%target_value            = 0.0_wp
    enefunc%restraint_posi          = .false.
    enefunc%restraint_flag          = .false.
    
    enefunc%table%table             = .false.
    enefunc%table%density           = 0.0_wp
    enefunc%table%water_model       = ''
    enefunc%table%num_water         = 0
    enefunc%table%num_solute        = 0
    enefunc%table%atom_cls_no_O     = 0
    enefunc%table%atom_cls_no_H     = 0
    enefunc%table%charge_O          = 0.0_wp
    enefunc%table%charge_H          = 0.0_wp
    enefunc%table%mass_O            = 0.0_wp
    enefunc%table%mass_H            = 0.0_wp

    enefunc%force_switch            = .false.
    enefunc%vdw_shift               = .false.

    enefunc%nonb_func               = 0
    enefunc%fudge_lj                = 1.0_wp
    enefunc%fudge_qq                = 1.0_wp
    enefunc%excl_level              = 3

    enefunc%multi_angle             = .false.
    enefunc%pme_use                 = .false.
    enefunc%num_basins              = 1

    enefunc%rpath_flag              = .false.
    enefunc%rpath_sum_mf_flag       = .false.
    enefunc%stats_count             = 0
    enefunc%stats_natom             = 0
    enefunc%stats_dimension         = 0
    enefunc%morph_ene_flag          = MorphEneFlagTot

    enefunc%num_fitting             = 0
    enefunc%fitting_method          = 0

    enefunc%fitting_file            = 0
    enefunc%fitting_move            = 0

    enefunc%pressure_rmsd           = .false.
    enefunc%pressure_position       = .false.
    enefunc%ez_membrane_flag        = .false.

    enefunc%eef1_use                = .false.
    enefunc%imm1_use                = .false.
    enefunc%imic_use                = .false.
    enefunc%gbsa_use                = .false.

    enefunc%num_morph_bb            = 0
    enefunc%num_morph_sc            = 0

    enefunc%cg_pwmcos_calc          = .false.
    enefunc%cg_pwmcosns_calc        = .false.
    enefunc%cg_safe_dihedral_calc   = .false.
    enefunc%cg_ele_calc             = .false.
    enefunc%cg_DNA_base_pair_calc   = .false.
    enefunc%cg_DNA_exv_calc         = .false.
    enefunc%cg_IDR_HPS_calc         = .false.
    enefunc%cg_IDR_KH_calc          = .false.
    enefunc%cg_KH_calc              = .false.

    enefunc%gamd_use                = .false.

    enefunc%eef1_use                = .false.
    enefunc%imm1_use                = .false.
    enefunc%imic_use                = .false.
    enefunc%gbsa_use                = .false.

    enefunc%gamd_use                = .false.

    enefunc%eef1_use                = .false.
    enefunc%imm1_use                = .false.
    enefunc%imic_use                = .false.
    enefunc%gbsa_use                = .false.
    
    enefunc%cg_safe_dih_ene_shift   = 0.0_wp

    enefunc%num_cg_particle_DNA_all     = 0
    enefunc%num_cg_particle_DNA_phos    = 0
    enefunc%num_cg_particle_DNA_base    = 0
    enefunc%num_cg_particle_IDR_HPS     = 0
    enefunc%num_cg_particle_IDR_KH      = 0
    enefunc%num_cg_particle_KH          = 0
    enefunc%num_cg_particle_charged     = 0
    enefunc%num_cg_particle_pro_charged = 0

    enefunc%rpath_sum_avforce       = .false.
    enefunc%mep_natoms              = 0

    return

  end subroutine init_enefunc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_enefunc
  !> @brief        allocate energy functions information
  !! @authors      YS, TY, CK
  !! @param[inout] enefunc   : structure of EneFunc information
  !! @param[in]    variable  : selected variable
  !! @param[in]    var_size  : size of the selected variable
  !! @param[in]    var_size2 : 2nd size of the selected variable (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_enefunc(enefunc, variable, var_size, var_size2)

    ! formal arguments
    type(s_enefunc),         intent(inout) :: enefunc
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size
    integer,       optional, intent(in)    :: var_size2

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat
    integer                  :: var_size3, var_size4


    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case(EneFuncBond)

      if (allocated(enefunc%bond_list)) then
        if (size(enefunc%bond_list(1,:)) == var_size) return
        deallocate(enefunc%bond_list,        &
                   enefunc%bond_force_const, &
                   enefunc%bond_dist_min,    &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%bond_list(2,var_size),      &
               enefunc%bond_force_const(var_size), &
               enefunc%bond_dist_min(var_size),    &
               stat = alloc_stat)

      enefunc%bond_list   (1:2,1:var_size) = 0
      enefunc%bond_force_const(1:var_size) = 0.0_wp
      enefunc%bond_dist_min   (1:var_size) = 0.0_wp

    ! ~CG~ 3SPN.2C DNA: quartic-form bond potential
    case(EneFuncBondQuartic)

      if (allocated(enefunc%bond_quartic_list)) then
         if (size(enefunc%bond_quartic_list(1,:)) == var_size) return
         deallocate(enefunc%bond_quartic_list,        &
              enefunc%bond_quartic_force_const, &
              enefunc%bond_quartic_dist_min,    &
              stat = dealloc_stat)
      end if

      allocate(enefunc%bond_quartic_list(2,var_size),      &
           enefunc%bond_quartic_force_const(var_size), &
           enefunc%bond_quartic_dist_min(var_size),    &
           stat = alloc_stat)

      enefunc%bond_quartic_list   (1:2,1:var_size) = 0
      enefunc%bond_quartic_force_const(1:var_size) = 0.0_wp
      enefunc%bond_quartic_dist_min   (1:var_size) = 0.0_wp

    case(EneFuncAngl)

      if (allocated(enefunc%angl_list)) then
        if (size(enefunc%angl_list(1,:)) == var_size) return
        deallocate(enefunc%angl_list,         &
                   enefunc%angl_force_const,  &
                   enefunc%angl_theta_min,    &
                   enefunc%angl_force_const1, &
                   enefunc%angl_theta_min1,   &
                   enefunc%angl_gamma,        &
                   enefunc%angl_epsa,         &
                   ! ~CG~ AICG2+ angle
                   enefunc%angl_func,         &
                   !shinobu-edited
                   enefunc%angl_w,            &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%angl_list(3, var_size),      &
               enefunc%angl_force_const(var_size),  &
               enefunc%angl_theta_min(var_size),    &
               ! ~CG~ AICG2+ angle
               enefunc%angl_func(var_size),         &
               enefunc%angl_force_const1(var_size), &
               enefunc%angl_theta_min1(var_size),   &
               enefunc%angl_gamma(var_size),        &
               enefunc%angl_epsa(var_size),         &
               !shinobu-edited
               enefunc%angl_w(var_size),            &
               stat = alloc_stat)

      enefunc%angl_list  (1:3, 1:var_size)  = 0
      enefunc%angl_force_const(1:var_size)  = 0.0_wp
      enefunc%angl_theta_min  (1:var_size)  = 0.0_wp
      enefunc%angl_force_const1(1:var_size) = 0.0_wp
      enefunc%angl_theta_min1  (1:var_size) = 0.0_wp
      enefunc%angl_gamma(1:var_size)        = 0.0_wp
      enefunc%angl_epsa(1:var_size)         = 0.0_wp
      !shinobu-edited
      enefunc%angl_w(1:var_size)            = 0.0_wp
      ! ~CG~ AICG2+ angle
      ! Default: 1???
      enefunc%angl_func(1:var_size)         = 1

    case(EneFuncUrey)

      if (allocated(enefunc%urey_list)) then
        if (size(enefunc%urey_list(1,:)) == var_size) return
        deallocate(enefunc%urey_list,        &
                   enefunc%urey_force_const, &
                   enefunc%urey_rmin,        &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%urey_list(2, var_size),     &
               enefunc%urey_force_const(var_size), &
               enefunc%urey_rmin(var_size),        &
               stat = alloc_stat)

      enefunc%urey_list  (1:2, 1:var_size) = 0
      enefunc%urey_force_const(1:var_size) = 0.0_wp
      enefunc%urey_rmin       (1:var_size) = 0.0_wp

    !shinobu-edited
    case(EneFuncAngFlex)
      if (allocated(enefunc%anglflex_list)) then
        if (size(enefunc%anglflex_list(1,:)) == var_size) return
        deallocate(enefunc%anglflex_list,        &
                   enefunc%anglflex_type,        &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%anglflex_list(3, var_size),     &
               enefunc%anglflex_type(var_size),        &
               stat = alloc_stat)

      enefunc%anglflex_list  (1:3, 1:var_size) = 0
      enefunc%anglflex_type  (1:var_size) = 0

    case(EneFuncAngFlexTbl)
      if (allocated(enefunc%anglflex_theta)) then
        if (size(enefunc%anglflex_theta(1,:)) == var_size) return
        deallocate(enefunc%anglflex_theta,            &
                   enefunc%anglflex_efunc,            &
                   enefunc%anglflex_d2func,           &
                   enefunc%anglflex_min_th,           &
                   enefunc%anglflex_max_th,           &
                   enefunc%anglflex_ener_corr,        &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%anglflex_theta(var_size2, var_size),      &
               enefunc%anglflex_efunc(var_size2, var_size),      &
               enefunc%anglflex_d2func(var_size2, var_size),     &
               enefunc%anglflex_min_th(1:2, var_size),           &
               enefunc%anglflex_max_th(1:2, var_size),           &
               enefunc%anglflex_ener_corr(var_size),             &
               stat = alloc_stat)

      enefunc%anglflex_theta (1:var_size2,1:var_size) = 0.0_wp
      enefunc%anglflex_efunc (1:var_size2,1:var_size) = 0.0_wp
      enefunc%anglflex_d2func(1:var_size2,1:var_size) = 0.0_wp
      enefunc%anglflex_min_th(1:2,1:var_size)         = 0.0_wp
      enefunc%anglflex_max_th(1:2,1:var_size)         = 0.0_wp
      enefunc%anglflex_ener_corr(1:var_size)          = 0.0_wp

    case(EneFuncDihe)

      if (allocated(enefunc%dihe_list)) then
        if (size(enefunc%dihe_list(1,:)) == var_size) return
        deallocate(enefunc%dihe_list,        &
                   enefunc%dihe_force_const, &
                   enefunc%dihe_periodicity, &
                   enefunc%dihe_phase,       &
                   enefunc%dihe_scee,        &
                   enefunc%dihe_scnb,        &
                   !shinobu-edited
                   enefunc%dihe_theta_min,   &
                   enefunc%dihe_w,           &
                   ! ~CG~ 3SPN.2C DNA: dihedral functype
                   enefunc%dihe_func,        &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%dihe_list(4,var_size),      &
               enefunc%dihe_force_const(var_size), &
               enefunc%dihe_periodicity(var_size), &
               enefunc%dihe_phase(var_size),       &
               enefunc%dihe_scee(var_size),        &
               enefunc%dihe_scnb(var_size),        &
               !shinobu-edited
               enefunc%dihe_theta_min(var_size),   &
               enefunc%dihe_w(var_size),           &
               ! ~CG~ 3SPN.2C DNA: dihedral functype
               enefunc%dihe_func(var_size),        &
               stat = alloc_stat)

      enefunc%dihe_list   (1:4,1:var_size) = 0
      enefunc%dihe_force_const(1:var_size) = 0.0_wp
      enefunc%dihe_periodicity(1:var_size) = 0
      enefunc%dihe_phase      (1:var_size) = 0.0_wp
      enefunc%dihe_scee       (1:var_size) = 0.0_wp
      enefunc%dihe_scnb       (1:var_size) = 0.0_wp
      !shinobu-edited
      enefunc%dihe_theta_min  (1:var_size) = 0.0_wp
      enefunc%dihe_w          (1:var_size) = 0.0_wp
      ! ~CG~ 3SPN.2C DNA: dihedral func type...
      enefunc%dihe_func       (1:var_size) = 1

    !shinobu-edited
    case(EneFuncDiheFlex)
      if (allocated(enefunc%diheflex_list)) then
        if (size(enefunc%diheflex_list(1,:)) == var_size) return
        deallocate(enefunc%diheflex_list,      & 
                   enefunc%diheflex_type,      &
                   enefunc%diheflex_func,      &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%diheflex_list(4,var_size),    &
               enefunc%diheflex_type(var_size),      &
               enefunc%diheflex_func(var_size),      &
               stat = alloc_stat)

      enefunc%diheflex_list   (1:4,1:var_size) = 0
      enefunc%diheflex_type   (    1:var_size) = 0
      enefunc%diheflex_func   (    1:var_size) = 0

    case(EneFuncDiheFlexTbl)
      if (allocated(enefunc%diheflex_coef)) then
        if (size(enefunc%diheflex_coef(1,:)) == var_size) return
        deallocate(enefunc%diheflex_coef,            &
                   enefunc%diheflex_ener_corr,       &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%diheflex_coef(var_size2, var_size),      &
               enefunc%diheflex_ener_corr(var_size),            &
               stat = alloc_stat)

      enefunc%diheflex_coef  (1:var_size2,1:var_size) = 0.0_wp
      enefunc%diheflex_ener_corr(1:var_size)          = 0.0_wp

    case(EneFuncRBDihe)

      if (allocated(enefunc%rb_dihe_list)) then
        if (size(enefunc%rb_dihe_list(1,:)) == var_size) return
        deallocate(enefunc%rb_dihe_list,     &
                   enefunc%rb_dihe_c,        &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%rb_dihe_list(4,var_size),   &
               enefunc%rb_dihe_c   (6,var_size),   &
               stat = alloc_stat)

      enefunc%rb_dihe_list(1:4,1:var_size) = 0
      enefunc%rb_dihe_c   (1:6,1:var_size) = 0.0_wp

    case(EneFuncImpr)

      if (allocated(enefunc%impr_list)) then
        if (size(enefunc%impr_list(1,:)) == var_size) return
        deallocate(enefunc%impr_list,        &
                   enefunc%impr_force_const, &
                   enefunc%impr_periodicity, &
                   enefunc%impr_phase,       &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%impr_list(4,var_size),      &
               enefunc%impr_force_const(var_size), &
               enefunc%impr_periodicity(var_size), &
               enefunc%impr_phase(var_size),       &
               stat = alloc_stat)

      enefunc%impr_list   (1:4,1:var_size) = 0
      enefunc%impr_force_const(1:var_size) = 0.0_wp
      enefunc%impr_periodicity(1:var_size) = 0
      enefunc%impr_phase      (1:var_size) = 0.0_wp

    case(EneFuncCmap)

      if (allocated(enefunc%cmap_list)) then
        if (size(enefunc%cmap_list(1,:)) == var_size) return
        deallocate(enefunc%cmap_list,       &
                   enefunc%cmap_type,       &
                   enefunc%cmap_force,      &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%cmap_list(8,var_size),    &
               enefunc%cmap_type(var_size),      &
               enefunc%cmap_force(3,8,var_size), &
               stat = alloc_stat)

      enefunc%cmap_list(1:8,1:var_size)   = 0
      enefunc%cmap_type(1:var_size) = 0
      enefunc%cmap_force(1:3,1:8,1:var_size) = 0.0_wp
      
    case(EneFuncCGGeneral)
      if (allocated(enefunc%NA_base_type)) then
        if (size(enefunc%NA_base_type)  == var_size .and. &
            size(enefunc%mol_chain_id)  == var_size .and. &
            size(enefunc%atom_cls)      == var_size .and. &
            size(enefunc%param_epsilon) == var_size .and. &
            size(enefunc%param_sigma)   == var_size) then
          return
        end if

        deallocate(enefunc%NA_base_type,       &
            enefunc%mol_chain_id,              &
            enefunc%atom_cls,                  &
            enefunc%param_epsilon,             &
            enefunc%param_sigma,               &
            stat = dealloc_stat)
      end if

      allocate(enefunc%NA_base_type(var_size), &
          enefunc%mol_chain_id(var_size),      &
          enefunc%atom_cls(var_size),          &
          enefunc%param_epsilon(var_size),     &
          enefunc%param_sigma(var_size),       &
          stat = alloc_stat)

      enefunc%NA_base_type (1:var_size) = 0
      enefunc%mol_chain_id (1:var_size) = 0
      enefunc%atom_cls     (1:var_size) = 0
      enefunc%param_epsilon(1:var_size) = 0
      enefunc%param_sigma  (1:var_size) = 0

    case(EneFuncBaseStack)
      ! ~CG~ 3SPN.2C DNA: base stacking
      if (allocated(enefunc%base_stack_list)) then
        if (size(enefunc%base_stack_list(1,:)) == var_size) return
        deallocate(enefunc%base_stack_list,           &
                   enefunc%base_stack_func,           &
                   enefunc%base_stack_epsilon,        &
                   enefunc%base_stack_sigma,          &
                   enefunc%base_stack_theta_bs,       &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%base_stack_list(3,var_size),   &
               enefunc%base_stack_func(var_size),     &
               enefunc%base_stack_epsilon(var_size),  &
               enefunc%base_stack_sigma(var_size),    &
               enefunc%base_stack_theta_bs(var_size), &
               stat = alloc_stat)

      enefunc%base_stack_list     (1:3, 1:var_size) = 0
      enefunc%base_stack_func     (1:var_size     ) = 0
      enefunc%base_stack_epsilon  (1:var_size     ) = 0.0_wp
      enefunc%base_stack_sigma    (1:var_size     ) = 0.0_wp
      enefunc%base_stack_theta_bs (1:var_size     ) = 0.0_wp

    case(EneFuncBasePair)
      ! ~CG~ 3SPN.2C DNA: base pairing
      ! var_size = 4 for DNA or RNA
      ! var_size = 8 for DNA + RNA
      if (allocated(enefunc%base_pair_theta_1)) then
        if (size(enefunc%base_pair_theta_1) == var_size) return
        deallocate(enefunc%base_pair_theta_1,     &
                   enefunc%base_pair_theta_2,     &
                   enefunc%base_pair_theta_3,     &
                   enefunc%base_pair_phi_1,       &
                   enefunc%base_pair_sigma,       &
                   enefunc%base_pair_epsilon,     &
                   enefunc%base_cross_1_epsilon,  &
                   enefunc%base_cross_1_sigma,    &
                   enefunc%base_cross_1_theta_cs, &
                   enefunc%base_cross_2_epsilon,  &
                   enefunc%base_cross_2_sigma,    &
                   enefunc%base_cross_2_theta_cs, &
                   enefunc%base_pair_is_WC,       &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%base_pair_theta_1     (var_size),          &
               enefunc%base_pair_theta_2     (var_size),          &
               enefunc%base_pair_theta_3     (var_size),          &
               enefunc%base_pair_phi_1       (var_size),          &
               enefunc%base_pair_sigma       (var_size),          &
               enefunc%base_pair_epsilon     (var_size),          &
               enefunc%base_cross_1_epsilon  (var_size,var_size), &
               enefunc%base_cross_1_sigma    (var_size,var_size), &
               enefunc%base_cross_1_theta_cs (var_size,var_size), &
               enefunc%base_cross_2_epsilon  (var_size,var_size), &
               enefunc%base_cross_2_sigma    (var_size,var_size), &
               enefunc%base_cross_2_theta_cs (var_size,var_size), &
               enefunc%base_pair_is_WC       (var_size,var_size), &
               stat = alloc_stat)

      enefunc%base_pair_theta_1     (1:var_size) = 0.0_wp
      enefunc%base_pair_theta_2     (1:var_size) = 0.0_wp
      enefunc%base_pair_theta_3     (1:var_size) = 0.0_wp
      enefunc%base_pair_phi_1       (1:var_size) = 0.0_wp
      enefunc%base_pair_sigma       (1:var_size) = 0.0_wp
      enefunc%base_pair_epsilon     (1:var_size) = 0.0_wp
      enefunc%base_cross_1_epsilon  (1:var_size,1:var_size) = 0.0_wp
      enefunc%base_cross_1_sigma    (1:var_size,1:var_size) = 0.0_wp
      enefunc%base_cross_1_theta_cs (1:var_size,1:var_size) = 0.0_wp
      enefunc%base_cross_2_epsilon  (1:var_size,1:var_size) = 0.0_wp
      enefunc%base_cross_2_sigma    (1:var_size,1:var_size) = 0.0_wp
      enefunc%base_cross_2_theta_cs (1:var_size,1:var_size) = 0.0_wp
      enefunc%base_pair_is_WC       (1:var_size,1:var_size) = .true.

    case(EneFuncCGDNAExv)

      if (allocated(enefunc%cgDNA_exv_sigma)) then
        if (size(enefunc%cgDNA_exv_sigma(1,:)) == var_size) return
        deallocate(enefunc%cgDNA_exv_sigma, &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%cgDNA_exv_sigma(var_size,var_size), &
               stat = alloc_stat)

      enefunc%cgDNA_exv_sigma(1:var_size,1:var_size) = 0.0_wp
      
    case(EneFuncCGele)

      if (allocated(enefunc%cg_charge) .and. &
          allocated(enefunc%cg_ele_mol_pair)) then
        if (size(enefunc%cg_charge(:)) == var_size .and.           &
            size(enefunc%cg_ele_mol_pair(1, :)) == var_size2 .and. &
            size(enefunc%cg_ele_mol_pair(:, 1)) == var_size2) return
        deallocate(enefunc%cg_charge, &
            enefunc%cg_ele_mol_pair,  &
            stat = dealloc_stat)
      end if

      allocate(enefunc%cg_charge(var_size),              &
          enefunc%cg_ele_mol_pair(var_size2, var_size2), &
          stat = alloc_stat)

      enefunc%cg_charge(1:var_size) = 0.0_wp
      enefunc%cg_ele_mol_pair(:, :) = 0

    case(EneFuncCGKH)
      ! ~CG~ protein KH model
      if (allocated(enefunc%cg_KH_sigma_half)) then
        if ((size(enefunc%cg_KH_sigma_half) == var_size) .and. &
            size(enefunc%cg_KH_epsilon(1, :)) == var_size2 .and. &
            size(enefunc%cg_KH_epsilon(:, 1)) == var_size2) return
        deallocate(enefunc%cg_KH_sigma_half, &
            enefunc%cg_pro_use_KH,           &
            enefunc%cg_KH_epsilon,           &
            stat = dealloc_stat)
      end if

      allocate(enefunc%cg_KH_sigma_half(var_size),     &
          enefunc%cg_pro_use_KH(var_size),             &
          enefunc%cg_KH_epsilon(var_size2, var_size2), &
          stat = alloc_stat)

      enefunc%cg_pro_use_KH(1:var_size)                  = .false.
      enefunc%cg_KH_sigma_half(1:var_size)               = 0.0
      enefunc%cg_KH_epsilon   (1:var_size2, 1:var_size2) = 0.0

    case(EneFuncCGKHmol)

      if (allocated(enefunc%cg_KH_mol_pair)) then
        if (size(enefunc%cg_KH_mol_pair(1, :)) == var_size .and. &
            size(enefunc%cg_KH_mol_pair(:, 1)) == var_size) return
        deallocate(enefunc%cg_KH_mol_pair,  &
            stat = dealloc_stat)
      end if

      allocate(enefunc%cg_KH_mol_pair(var_size, var_size), &
          stat = alloc_stat)

      enefunc%cg_KH_mol_pair(:, :) = 0

    case(EneFuncPWMcos)
      ! ~CG~ protein-DNA interaction
      if (allocated(enefunc%pwmcos_protein_id)) then
        if (size(enefunc%pwmcos_protein_id) == var_size) return
        deallocate(enefunc%pwmcos_protein_id, &
            enefunc%pwmcos_protein_id_N,      &
            enefunc%pwmcos_protein_id_C,      &
            enefunc%pwmcos_r0,                &
            enefunc%pwmcos_theta1,            &
            enefunc%pwmcos_theta2,            &
            enefunc%pwmcos_theta3,            &
            enefunc%pwmcos_ene_A,             &
            enefunc%pwmcos_ene_C,             &
            enefunc%pwmcos_ene_G,             &
            enefunc%pwmcos_ene_T,             &
            enefunc%pwmcos_gamma,             &
            enefunc%pwmcos_eps,               &
            enefunc%pwmcos_specificity,       &
            enefunc%pwmcos_to_pairlist_id,    &
            stat = dealloc_stat)
      end if

      allocate(enefunc%pwmcos_protein_id (var_size), &
          enefunc%pwmcos_protein_id_N   (var_size),  &
          enefunc%pwmcos_protein_id_C   (var_size),  &
          enefunc%pwmcos_r0             (var_size),  &
          enefunc%pwmcos_theta1         (var_size),  &
          enefunc%pwmcos_theta2         (var_size),  &
          enefunc%pwmcos_theta3         (var_size),  &
          enefunc%pwmcos_ene_A          (var_size),  &
          enefunc%pwmcos_ene_C          (var_size),  &
          enefunc%pwmcos_ene_G          (var_size),  &
          enefunc%pwmcos_ene_T          (var_size),  &
          enefunc%pwmcos_gamma          (var_size),  &
          enefunc%pwmcos_eps            (var_size),  &
          enefunc%pwmcos_specificity    (var_size),  &
          enefunc%pwmcos_to_pairlist_id (var_size),  &
          stat = alloc_stat)

      enefunc%pwmcos_protein_id     (1:var_size) = 0
      enefunc%pwmcos_protein_id_N   (1:var_size) = 0
      enefunc%pwmcos_protein_id_C   (1:var_size) = 0
      enefunc%pwmcos_r0             (1:var_size) = 0.0_wp
      enefunc%pwmcos_theta1         (1:var_size) = 0.0_wp
      enefunc%pwmcos_theta2         (1:var_size) = 0.0_wp
      enefunc%pwmcos_theta3         (1:var_size) = 0.0_wp
      enefunc%pwmcos_ene_A          (1:var_size) = 0.0_wp
      enefunc%pwmcos_ene_C          (1:var_size) = 0.0_wp
      enefunc%pwmcos_ene_G          (1:var_size) = 0.0_wp
      enefunc%pwmcos_ene_T          (1:var_size) = 0.0_wp
      enefunc%pwmcos_gamma          (1:var_size) = 0.0_wp
      enefunc%pwmcos_eps            (1:var_size) = 0.0_wp
      enefunc%pwmcos_specificity    (1:var_size) = 1
      enefunc%pwmcos_to_pairlist_id (1:var_size) = 1

    case(EneFuncPWMcosns)
      ! ~CG~ protein-DNA interaction
      if (allocated(enefunc%pwmcosns_protein_id)) then
        if (size(enefunc%pwmcosns_protein_id) == var_size) return
        deallocate(enefunc%pwmcosns_protein_id, &
            enefunc%pwmcosns_protein_id_N,      &
            enefunc%pwmcosns_protein_id_C,      &
            enefunc%pwmcosns_r0,                &
            enefunc%pwmcosns_theta1,            &
            enefunc%pwmcosns_theta2,            &
            enefunc%pwmcosns_ene,               &
            enefunc%pwmcosns_specificity,       &
            enefunc%pwmcosns_to_pairlist_id,    &
            stat = dealloc_stat)
      end if

      allocate(enefunc%pwmcosns_protein_id (var_size), &
          enefunc%pwmcosns_protein_id_N   (var_size),  &
          enefunc%pwmcosns_protein_id_C   (var_size),  &
          enefunc%pwmcosns_r0             (var_size),  &
          enefunc%pwmcosns_theta1         (var_size),  &
          enefunc%pwmcosns_theta2         (var_size),  &
          enefunc%pwmcosns_ene            (var_size),  &
          enefunc%pwmcosns_specificity    (var_size),  &
          enefunc%pwmcosns_to_pairlist_id (var_size),  &
          stat = alloc_stat)

      enefunc%pwmcosns_protein_id     (1:var_size) = 0
      enefunc%pwmcosns_protein_id_N   (1:var_size) = 0
      enefunc%pwmcosns_protein_id_C   (1:var_size) = 0
      enefunc%pwmcosns_r0             (1:var_size) = 0.0_wp
      enefunc%pwmcosns_theta1         (1:var_size) = 0.0_wp
      enefunc%pwmcosns_theta2         (1:var_size) = 0.0_wp
      enefunc%pwmcosns_ene            (1:var_size) = 0.0_wp
      enefunc%pwmcosns_specificity    (1:var_size) = 1
      enefunc%pwmcosns_to_pairlist_id (1:var_size) = 1

    case(EneFuncCGIDRHPS)
      ! ~CG~ protein IDR HPS model
      if (allocated(enefunc%cg_IDR_HPS_is_IDR)) then
        if (size(enefunc%cg_IDR_HPS_is_IDR) == var_size) return
        deallocate(enefunc%cg_IDR_HPS_is_IDR, &
            enefunc%cg_IDR_HPS_lambda_half,   &
            enefunc%cg_IDR_HPS_sigma_half,    &
            stat = dealloc_stat)
      end if

      allocate(enefunc%cg_IDR_HPS_is_IDR(var_size), &
          enefunc%cg_IDR_HPS_lambda_half(var_size), &
          enefunc%cg_IDR_HPS_sigma_half(var_size),  &
          stat = alloc_stat)

      enefunc%cg_IDR_HPS_is_IDR     (1:var_size) = .false.
      enefunc%cg_IDR_HPS_lambda_half(1:var_size) = 0.0
      enefunc%cg_IDR_HPS_sigma_half(1:var_size)  = 0.0

    case(EneFuncCGIDRKH)
      ! ~CG~ protein IDR KH model
      if (allocated(enefunc%cg_IDR_KH_is_IDR)) then
        if ((size(enefunc%cg_IDR_KH_is_IDR) == var_size) .and. &
            size(enefunc%cg_IDR_KH_epsilon_D(1, :)) == var_size2 .and. &
            size(enefunc%cg_IDR_KH_epsilon_D(:, 1)) == var_size2) return
        deallocate(enefunc%cg_IDR_KH_is_IDR, &
            enefunc%cg_IDR_KH_sigma_half,    &
            enefunc%cg_IDR_KH_epsilon_D,     &
            stat = dealloc_stat)
      end if

      allocate(enefunc%cg_IDR_KH_is_IDR(var_size), &
          enefunc%cg_IDR_KH_sigma_half(var_size),  &
          enefunc%cg_IDR_KH_epsilon_D(var_size2, var_size2),  &
          stat = alloc_stat)

      enefunc%cg_IDR_KH_is_IDR    (1:var_size)               = .false.
      enefunc%cg_IDR_KH_sigma_half(1:var_size)               = 0.0
      enefunc%cg_IDR_KH_epsilon_D (1:var_size2, 1:var_size2) = 0.0

    case(EneFuncNbon)

      if (allocated(enefunc%nonb_atom_cls)) then
        if (size(enefunc%nonb_atom_cls(:)) == var_size) return
        deallocate(enefunc%nonb_atom_cls, &
                   enefunc%nb14_lj6,      &
                   enefunc%nb14_lj10,     &
                   enefunc%nb14_lj12,     &
                   enefunc%nonb_lj6,      &
                   enefunc%nonb_lj10,     &
                   enefunc%nonb_lj12,     &
                   enefunc%nonb_aicg_sig, &
                   enefunc%nonb_aicg_eps, &
                   enefunc%cg_exv_sig_half, &
                   enefunc%cg_exv_eps_sqrt, &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%nonb_atom_cls(var_size),          &
               enefunc%nb14_lj6(var_size,var_size),      &
               enefunc%nb14_lj10(var_size,var_size),     &
               enefunc%nb14_lj12(var_size,var_size),     &
               enefunc%nonb_lj6(var_size,var_size),      &
               enefunc%nonb_lj10(var_size,var_size),     &
               enefunc%nonb_lj12(var_size,var_size),     &
               enefunc%nonb_aicg_sig(var_size,var_size), &
               enefunc%nonb_aicg_eps(var_size,var_size), &
               enefunc%cg_exv_sig_half(var_size),        &
               enefunc%cg_exv_eps_sqrt(var_size),        &
               stat = alloc_stat)

      enefunc%nonb_atom_cls(1:var_size) = 0
      enefunc%nb14_lj6     (1:var_size,1:var_size) = 0.0_wp
      enefunc%nb14_lj10    (1:var_size,1:var_size) = 0.0_wp
      enefunc%nb14_lj12    (1:var_size,1:var_size) = 0.0_wp
      enefunc%nonb_lj6     (1:var_size,1:var_size) = 0.0_wp
      enefunc%nonb_lj10    (1:var_size,1:var_size) = 0.0_wp
      enefunc%nonb_lj12    (1:var_size,1:var_size) = 0.0_wp
      enefunc%nonb_aicg_sig(1:var_size,1:var_size) = 0.0_wp
      enefunc%nonb_aicg_eps(1:var_size,1:var_size) = 0.0_wp
      enefunc%cg_exv_eps_sqrt(1:var_size) = 0.0_wp
      enefunc%cg_exv_sig_half(1:var_size) = 2.0_wp

    case(EneFuncCntc)

      if (allocated(enefunc%contact_list)) then
        if (size(enefunc%contact_list(1,:)) == var_size) return
        deallocate(enefunc%contact_list,  &
                   enefunc%contact_func,  &
                   enefunc%contact_lj6,   &
                   enefunc%contact_lj10,  &
                   enefunc%contact_lj12,  &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%contact_list(2,var_size),   &
               enefunc%contact_func(var_size),     &
               enefunc%contact_lj6(var_size),      &
               enefunc%contact_lj10(var_size),     &
               enefunc%contact_lj12(var_size),     &
               stat = alloc_stat)

      enefunc%contact_list(1:2,1:var_size) = 0
      enefunc%contact_func(1:var_size)  = 0
      enefunc%contact_lj6 (1:var_size)  = 0.0_wp
      enefunc%contact_lj10(1:var_size)  = 0.0_wp
      enefunc%contact_lj12(1:var_size)  = 0.0_wp

    case(EneFuncVsite2)

      if (allocated(enefunc%vsite2_list)) then
        if (size(enefunc%vsite2_list(1,:)) == var_size) return
        deallocate(enefunc%vsite2_list, &
                   enefunc%vsite2_a,    &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%vsite2_list(3,var_size), &
               enefunc%vsite2_a   (  var_size), &
               stat = alloc_stat)

      enefunc%vsite2_list(1:3,1:var_size) = 0
      enefunc%vsite2_a   (    1:var_size) = 0.0_wp

    case(EneFuncVsite3)

      if (allocated(enefunc%vsite3_list)) then
        if (size(enefunc%vsite3_list(1,:)) == var_size) return
        deallocate(enefunc%vsite3_list, &
                   enefunc%vsite3_a,    &
                   enefunc%vsite3_b,    &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%vsite3_list(4,var_size), &
               enefunc%vsite3_a   (  var_size), &
               enefunc%vsite3_b   (  var_size), &
               stat = alloc_stat)

      enefunc%vsite3_list(1:4,1:var_size) = 0
      enefunc%vsite3_a   (    1:var_size) = 0.0_wp
      enefunc%vsite3_b   (    1:var_size) = 0.0_wp

    case(EneFuncVsite3fd)

      if (allocated(enefunc%vsite3fd_list)) then
        if (size(enefunc%vsite3fd_list(1,:)) == var_size) return
        deallocate(enefunc%vsite3fd_list, &
                   enefunc%vsite3fd_a,    &
                   enefunc%vsite3fd_d,    &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%vsite3fd_list(4,var_size), &
               enefunc%vsite3fd_a   (  var_size), &
               enefunc%vsite3fd_d   (  var_size), &
               stat = alloc_stat)

      enefunc%vsite3fd_list(1:4,1:var_size) = 0
      enefunc%vsite3fd_a   (    1:var_size) = 0.0_wp
      enefunc%vsite3fd_d   (    1:var_size) = 0.0_wp

    case(EneFuncVsite3fad)

      if (allocated(enefunc%vsite3fad_list)) then
        if (size(enefunc%vsite3fad_list(1,:)) == var_size) return
        deallocate(enefunc%vsite3fad_list,  &
                   enefunc%vsite3fad_theta, &
                   enefunc%vsite3fad_d,     &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%vsite3fad_list (4,var_size), &
               enefunc%vsite3fad_theta(  var_size), &
               enefunc%vsite3fad_d    (  var_size), &
               stat = alloc_stat)

      enefunc%vsite3fad_list (1:4,1:var_size) = 0
      enefunc%vsite3fad_theta(    1:var_size) = 0.0_wp
      enefunc%vsite3fad_d    (    1:var_size) = 0.0_wp

    case(EneFuncVsite3out)

      if (allocated(enefunc%vsite3out_list)) then
        if (size(enefunc%vsite3out_list(1,:)) == var_size) return
        deallocate(enefunc%vsite3out_list, &
                   enefunc%vsite3out_a,    &
                   enefunc%vsite3out_b,    &
                   enefunc%vsite3out_c,    &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%vsite3out_list(4,var_size), &
               enefunc%vsite3out_a   (  var_size), &
               enefunc%vsite3out_b   (  var_size), &
               enefunc%vsite3out_c   (  var_size), &
               stat = alloc_stat)

      enefunc%vsite3out_list(1:4,1:var_size) = 0
      enefunc%vsite3out_a   (    1:var_size) = 0.0_wp
      enefunc%vsite3out_b   (    1:var_size) = 0.0_wp
      enefunc%vsite3out_c   (    1:var_size) = 0.0_wp

    case(EneFuncVsite4fdn)

      if (allocated(enefunc%vsite4fdn_list)) then
        if (size(enefunc%vsite4fdn_list(1,:)) == var_size) return
        deallocate(enefunc%vsite4fdn_list, &
                   enefunc%vsite4fdn_a,    &
                   enefunc%vsite4fdn_b,    &
                   enefunc%vsite4fdn_c,    &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%vsite4fdn_list(5,var_size), &
               enefunc%vsite4fdn_a   (  var_size), &
               enefunc%vsite4fdn_b   (  var_size), &
               enefunc%vsite4fdn_c   (  var_size), &
               stat = alloc_stat)

      enefunc%vsite4fdn_list(1:5,1:var_size) = 0
      enefunc%vsite4fdn_a   (    1:var_size) = 0.0_wp
      enefunc%vsite4fdn_b   (    1:var_size) = 0.0_wp
      enefunc%vsite4fdn_c   (    1:var_size) = 0.0_wp

    case(EneFuncVsiten)

      if (allocated(enefunc%vsiten_list)) then
        if (size(enefunc%vsiten_list(1,:)) == var_size .and. &
            size(enefunc%vsiten_list(:,1)) == var_size2) return
        deallocate(enefunc%vsiten_list, &
                   enefunc%vsiten_n,    &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%vsiten_list(var_size2,var_size), &
               enefunc%vsiten_n   (          var_size), &
               stat = alloc_stat)

      enefunc%vsiten_list(1:var_size2,1:var_size) = 0
      enefunc%vsiten_n   (            1:var_size) = 0

    case(EneFuncRefg)

      if (allocated(enefunc%restraint_numatoms)) then
        if (size(enefunc%restraint_numatoms) == var_size) return
        deallocate(enefunc%restraint_numatoms,      &
                   enefunc%restraint_atomlist,      &
                   enefunc%restraint_masscoef,      &
                   enefunc%restraint_totmass_group, &
                   enefunc%restraint_wcom3,         &
                   enefunc%restraint_wcom4,         &
                   enefunc%restraint_wcom5,         &
                   enefunc%restraint_wtmp,          &
                   enefunc%restraint_istartend, &
                   stat = dealloc_stat)
      end if
      allocate(enefunc%restraint_numatoms(1:var_size),             &
               enefunc%restraint_atomlist(1:var_size2,1:var_size), &
               enefunc%restraint_masscoef(1:var_size2,1:var_size), &
               enefunc%restraint_totmass_group(1:var_size),        &
               enefunc%restraint_wcom3(1:3,1:var_size2),           &
               enefunc%restraint_wcom4(1:3,1:var_size2),           &
               enefunc%restraint_wcom5(1:3,1:var_size2),           &
               enefunc%restraint_wtmp(1:var_size2,1:var_size),     &
               enefunc%restraint_istartend(1:2,1:var_size),        &
               stat = alloc_stat)

      enefunc%restraint_numatoms(1:var_size) = 0
      enefunc%restraint_atomlist(1:var_size2,1:var_size) = 0
      enefunc%restraint_masscoef(1:var_size2,1:var_size) = 0.0_wp
      enefunc%restraint_totmass_group(1:var_size)    = 0.0_wp
      enefunc%restraint_wcom3(1:3,1:var_size2)       = 0.0_wp
      enefunc%restraint_wcom4(1:3,1:var_size2)       = 0.0_wp
      enefunc%restraint_wcom5(1:3,1:var_size2)       = 0.0_wp
      enefunc%restraint_wtmp(1:var_size2,1:var_size) = 0.0_wp
      enefunc%restraint_istartend(1:2,1:var_size)    = 0

    case(EneFuncReff)

      if (allocated(enefunc%restraint_kind)) then
        if (size(enefunc%restraint_kind) == var_size) return
        deallocate(enefunc%restraint_kind,          &
                   enefunc%restraint_grouplist,     &
                   enefunc%restraint_const,         &
                   enefunc%restraint_ref,           &
                   enefunc%restraint_funcgrp,       &
                   enefunc%restraint_exponent_func, &
                   enefunc%restraint_exponent_dist, &
                   enefunc%restraint_mode,          &
                   enefunc%restraint_caging,        &
                   enefunc%restraint_flat_radius,   &
                   enefunc%restraint_replica_index, &
                   enefunc%restraint_weight_dist,   &
                   enefunc%restraint_wcom1,         &
                   enefunc%restraint_wcom2,         &
                   enefunc%restraint_wdrt,          &
                   enefunc%restraint_rgbuffer,      &
                   enefunc%restraint_wall_z,        &
                   enefunc%restraint_rpath_func,    &
                   enefunc%restraint_rcom1,         &
                   enefunc%restraint_rcom2,         &
                   enefunc%restraint_rdrt,          &
                   stat = dealloc_stat)
      end if

      var_size3 = max(int(var_size2/2),1)
      var_size4 = int(var_size2*(var_size2-1)/2)
      allocate(enefunc%restraint_kind(var_size),                        &
               enefunc%restraint_grouplist(1:var_size2,1:var_size),     &
               enefunc%restraint_const(1:4,1:var_size),                 &
               enefunc%restraint_ref(1:2,1:var_size),                   &
               enefunc%restraint_funcgrp(1:var_size),                   &
               enefunc%restraint_exponent_func(1:var_size),             &
               enefunc%restraint_exponent_dist(1:var_size3,1:var_size), &
               enefunc%restraint_mode(1:var_size),                      &
               enefunc%restraint_replica_index(1:var_size),             &
               enefunc%restraint_caging(1:var_size),                    &
               enefunc%restraint_flat_radius(1:var_size),               &
               enefunc%restraint_weight_dist(1:var_size3,1:var_size),   &
               enefunc%restraint_wcom1(1:3,1:var_size2),                &
               enefunc%restraint_wcom2(1:3,1:var_size2),                &
               enefunc%restraint_wdrt(1:var_size3),                     &
               enefunc%restraint_rgbuffer(1:var_size),                  &
               enefunc%restraint_wall_z(1:2,1:var_size),                &
               enefunc%restraint_rpath_func(1:var_size),                &
               enefunc%restraint_rcom1(1:3,1:var_size2),                &
               enefunc%restraint_rcom2(1:3,1:var_size4),                &
               enefunc%restraint_rdrt(1:var_size4),                     &
               stat = alloc_stat)

      enefunc%restraint_kind(1:var_size) = 0
      enefunc%restraint_grouplist(1:var_size2,1:var_size) = 0
      enefunc%restraint_const(1:4,1:var_size) = 0.0_wp
      enefunc%restraint_ref(1:2,1:var_size) = 0.0_wp
      enefunc%restraint_funcgrp(1:var_size) = 0
      enefunc%restraint_exponent_func(1:var_size) = 0
      enefunc%restraint_exponent_dist(1:var_size3,1:var_size) = 0
      enefunc%restraint_mode(1:var_size) = 0
      enefunc%restraint_replica_index(1:var_size) = 0
      enefunc%restraint_caging(1:var_size) = .false.
      enefunc%restraint_flat_radius(1:var_size) = 0.0_wp
      enefunc%restraint_weight_dist(1:var_size3,1:var_size) = 0.0_wp
      enefunc%restraint_wcom1(1:3,1:var_size2) = 0.0_wp
      enefunc%restraint_wcom2(1:3,1:var_size2) = 0.0_wp
      enefunc%restraint_wdrt(1:var_size3) = 0.0_wp
      enefunc%restraint_rgbuffer(1:var_size)   = 0.0_wp
      enefunc%restraint_wall_z(1:2,1:var_size) = 0.0_wp
      enefunc%restraint_rpath_func(1:var_size) = 0
      enefunc%restraint_rcom1(1:3,1:var_size2) = 0.0_wp
      enefunc%restraint_rcom2(1:3,1:var_size4) = 0.0_wp
      enefunc%restraint_rdrt(1:var_size4) = 0.0_wp

    case(EneFuncRefw)

      if (allocated(enefunc%restraint_wfrc)) then
        if (size(enefunc%restraint_wfrc(1,:)) == var_size) return
        deallocate(enefunc%restraint_wfrc, &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%restraint_wfrc(1:3,var_size), &
               stat = alloc_stat)

      enefunc%restraint_wfrc(1:3,1:var_size) = 0.0_wp

    case(EneFuncRefc)

      if (allocated(enefunc%restraint_refcoord)) then
        if (size(enefunc%restraint_refcoord) == var_size) return
        deallocate(enefunc%restraint_refcoord, &
                   enefunc%restraint_masstmp,  &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%restraint_refcoord(1:3,var_size), &
               enefunc%restraint_masstmp(var_size),      &
               stat = alloc_stat)

      enefunc%restraint_refcoord(1:3,1:var_size) = 0.0_wp
      enefunc%restraint_masstmp(1:var_size)      = 0.0_wp

    case(EneFuncRefr)

      if (allocated(enefunc%restraint_const_replica)) then
        if (size(enefunc%restraint_const_replica(1,:)) == var_size) return
        deallocate(enefunc%restraint_const_replica, &
                   enefunc%restraint_ref_replica,  &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%restraint_const_replica(var_size2,var_size), &
               enefunc%restraint_ref_replica(var_size2,var_size), &
               stat = alloc_stat)

      enefunc%restraint_const_replica(1:var_size2,1:var_size) = 0.0_wp
      enefunc%restraint_ref_replica(1:var_size2,1:var_size)   = 0.0_wp

    case(EneFuncMode)

      if (allocated(enefunc%pc_mode)) then
        if (size(enefunc%pc_mode) == var_size) return
        deallocate(enefunc%pc_mode, &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%pc_mode(1:var_size), &
               stat = alloc_stat)
      enefunc%num_pc_modes        = 0
      enefunc%pc_mode(1:var_size) = 0.0_wp
   
    case(EneFuncTable)

      if (allocated(enefunc%table%table_ene_WW)) then
        if (size(enefunc%table%table_ene(:)) /= 6*var_size) &
          deallocate(enefunc%table%table_ene,    &
                     enefunc%table%table_grad,   &
                     enefunc%table%table_ecor,   &
                     enefunc%table%table_decor,  &
                     enefunc%table%dtable_ecor,  &
                     enefunc%table%dtable_decor, &
                     enefunc%table%table_ene_WW, &
                     enefunc%table%table_de_WW,  &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%table%table_ene_WW)) &
        allocate(enefunc%table%table_ene   (6*var_size),    &
                 enefunc%table%table_grad  (6*var_size),    &
                 enefunc%table%table_ecor  (var_size),      &
                 enefunc%table%table_decor (var_size),      &
                 enefunc%table%dtable_ecor (var_size),      &
                 enefunc%table%dtable_decor(var_size),      &
                 enefunc%table%table_ene_WW(6*var_size, 3), &
                 enefunc%table%table_de_WW (2*var_size, 3), &
                 stat = alloc_stat)

      enefunc%table%table_ene   (1:6*var_size)      = 0.0_wp
      enefunc%table%table_grad  (1:6*var_size)      = 0.0_wp
      enefunc%table%table_ecor  (1:var_size)        = 0.0_wp
      enefunc%table%table_decor (1:var_size)        = 0.0_wp
      enefunc%table%dtable_ecor (1:var_size)        = 0.0_wp
      enefunc%table%dtable_decor(1:var_size)        = 0.0_wp
      enefunc%table%table_ene_WW(1:6*var_size, 1:3) = 0.0_wp
      enefunc%table%table_de_WW (1:2*var_size, 1:3) = 0.0_wp

    case(EneFuncTableWater)

      if (allocated(enefunc%table%water_list)) then
        if (size(enefunc%table%water_list(1,:)) == var_size) return
        deallocate(enefunc%table%water_list, &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%table%water_list(3,var_size), &
               stat = alloc_stat)

      enefunc%table%water_list(1:3,1:var_size) = 0

    case(EneFuncTableSolute)

      if (allocated(enefunc%table%solute_list)) then
        if (size(enefunc%table%solute_list(:)) == var_size) return
        deallocate(enefunc%table%solute_list,     &
                   enefunc%table%solute_list_inv, &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%table%solute_list(var_size),      &
               enefunc%table%solute_list_inv(var_size2), &
               stat = alloc_stat)

      enefunc%table%solute_list    (1:var_size)  = 0
      enefunc%table%solute_list_inv(1:var_size2) = 0

    case(EneFuncCmapType)

      if (allocated(enefunc%cmap_resolution)) then
        if (size(enefunc%cmap_resolution(:)) == var_size) return
        deallocate(enefunc%cmap_resolution, &
                   enefunc%cmap_coef,       &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%cmap_resolution(var_size),                   &
               enefunc%cmap_coef(4,4,var_size2,var_size2,var_size), &
               stat = alloc_stat)

      enefunc%cmap_resolution(1:var_size) = 0
      enefunc%cmap_coef(1:4,1:4,1:var_size2,1:var_size2,1:var_size) = 0.0_wp

    case(EneFuncFitc)

      if (allocated(enefunc%fit_refcoord)) then
        if (size(enefunc%fit_refcoord(1,:)) == var_size) return
        deallocate(enefunc%fit_refcoord,    &
                   enefunc%fit_work,        &
                   enefunc%fit_mass,        &
                   enefunc%fitting_atom,    &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%fit_refcoord(1:3,1:var_size),                   &
               enefunc%fit_work(1:3,1:var_size),                       &
               enefunc%fit_mass(1:var_size),                       &
               enefunc%fitting_atom(1:var_size2),                      &
               stat = alloc_stat)

      enefunc%fit_refcoord(1:3,1:var_size) = 0.0_wp
      enefunc%fit_work(1:3,1:var_size)     = 0.0_wp
      enefunc%fit_mass(1:var_size)         = 0.0_wp
      enefunc%fitting_atom(1:var_size2)    = 0

    case(EneFuncModeRef)

      if (allocated(enefunc%pc_crd)) then
        if (size(enefunc%pc_crd) == var_size*3) return
        deallocate(enefunc%pc_crd,         &
                   enefunc%pc_crd_ref,     &
                   enefunc%pc_crd_ref_fit, &
                   enefunc%pc_mode_temp,   &
                   enefunc%pc_mode_fit,    &
                   enefunc%grad_pc,        &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%pc_crd(1:var_size*3),          &
               enefunc%pc_crd_ref(1:var_size*3),      &
               enefunc%pc_crd_ref_fit(1:var_size*3),  &
               enefunc%pc_mode_temp(1:var_size*3),    &
               enefunc%pc_mode_fit(1:var_size*3),     &
               enefunc%grad_pc(1:3,1:var_size),       &
               stat = alloc_stat)

      enefunc%pc_crd(1:var_size*3)         = 0.0_wp
      enefunc%pc_crd_ref(1:var_size*3)     = 0.0_wp
      enefunc%pc_crd_ref_fit(1:var_size*3) = 0.0_wp
      enefunc%pc_mode_temp(1:var_size*3)   = 0.0_wp
      enefunc%pc_mode_fit(1:var_size*3)    = 0.0_wp
      enefunc%grad_pc(1:3,1:var_size)      = 0.0_wp

    case(EneFuncEef1)

      if (allocated(enefunc%eef1%atom_name)) then
        if (size(enefunc%eef1%atom_name) == var_size) return
        deallocate(enefunc%eef1%atom_name,      &
                   enefunc%eef1%volume,         &
                   enefunc%eef1%gref_0,         &
                   enefunc%eef1%gfree_0,        &
                   enefunc%eef1%href_0,         &
                   enefunc%eef1%cpref_0,        &
                   enefunc%eef1%gref_t,         &
                   enefunc%eef1%gfree_t,        &
                   enefunc%eef1%alpha_4pi,      &
                   enefunc%eef1%inv_lambda,     &
                   enefunc%eef1%rvdw,           &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%eef1%atom_name (var_size),   &
               enefunc%eef1%volume    (2,var_size), &
               enefunc%eef1%gref_0    (2,var_size), &
               enefunc%eef1%gfree_0   (2,var_size), &
               enefunc%eef1%href_0    (2,var_size), &
               enefunc%eef1%cpref_0   (2,var_size), &
               enefunc%eef1%gref_t    (2,var_size), &
               enefunc%eef1%gfree_t   (2,var_size), &
               enefunc%eef1%inv_lambda(2,var_size), &
               enefunc%eef1%alpha_4pi (2,var_size), &
               enefunc%eef1%rvdw      (var_size),   &
               stat = alloc_stat)

      enefunc%eef1%atom_name (1:var_size)     = ''
      enefunc%eef1%volume    (1:2,1:var_size) = 0.0_wp
      enefunc%eef1%gref_0    (1:2,1:var_size) = 0.0_wp
      enefunc%eef1%gfree_0   (1:2,1:var_size) = 0.0_wp
      enefunc%eef1%href_0    (1:2,1:var_size) = 0.0_wp
      enefunc%eef1%cpref_0   (1:2,1:var_size) = 0.0_wp
      enefunc%eef1%gref_t    (1:2,1:var_size) = 0.0_wp
      enefunc%eef1%gfree_t   (1:2,1:var_size) = 0.0_wp
      enefunc%eef1%inv_lambda(1:2,1:var_size) = 0.0_wp
      enefunc%eef1%alpha_4pi (1:2,1:var_size) = 0.0_wp
      enefunc%eef1%rvdw      (1:var_size)     = 0.0_wp

    case(EneFuncGbsa)

      if (allocated(enefunc%gbsa%vdw_radius)) then
        if (size(enefunc%gbsa%vdw_radius) == var_size) return
        deallocate(enefunc%gbsa%vdw_radius,      &
                   enefunc%gbsa%scale_factor,    &
                   enefunc%gbsa%sasa_atom_type,  &
                   enefunc%gbsa%sasa_parameters, &
                   enefunc%gbsa%sasa_vdw_radius, &
                   enefunc%gbsa%sasa_atom_list,  &
                   enefunc%gbsa%born_radius,     &
                   enefunc%gbsa%dadrfac,         &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%gbsa%vdw_radius     (var_size),  &
               enefunc%gbsa%scale_factor   (var_size),  &
               enefunc%gbsa%sasa_atom_type (var_size),  &
               enefunc%gbsa%sasa_parameters(21,4),      &
               enefunc%gbsa%sasa_vdw_radius(var_size),  &
               enefunc%gbsa%sasa_atom_list (var_size2), &
               enefunc%gbsa%born_radius    (var_size),  &
               enefunc%gbsa%dadrfac        (var_size),  &
               stat = alloc_stat)

      enefunc%gbsa%vdw_radius     (1:var_size)  = 0.0_wp
      enefunc%gbsa%scale_factor   (1:var_size)  = 0.0_wp
      enefunc%gbsa%sasa_atom_type (1:var_size)  = 0
      enefunc%gbsa%sasa_parameters(1:21,1:4)    = 0.0_wp
      enefunc%gbsa%sasa_vdw_radius(1:var_size)  = 0.0_wp
      enefunc%gbsa%sasa_atom_list (1:var_size2) = 0
      enefunc%gbsa%born_radius    (1:var_size)  = 0.0_wp
      enefunc%gbsa%dadrfac        (1:var_size)  = 0.0_wp

    case(EneFuncEzMembrane)

      if (allocated(enefunc%ez_membrane_const)) then
        if (size(enefunc%ez_membrane_const) == var_size) return
        deallocate(enefunc%ez_membrane_const, &
                   enefunc%ez_membrane_zmin, &
                   enefunc%ez_membrane_sigma, &
                   enefunc%ez_membrane_polym, &
                   enefunc%ez_membrane_func, &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%ez_membrane_const(1:var_size), &
               enefunc%ez_membrane_zmin(1:var_size),  &
               enefunc%ez_membrane_sigma(1:var_size), &
               enefunc%ez_membrane_polym(1:var_size), &
               enefunc%ez_membrane_func(1:var_size), &
               stat = alloc_stat)
      enefunc%ez_membrane_const(1:var_size) = 0.0_wp
      enefunc%ez_membrane_zmin(1:var_size)  = 0.0_wp
      enefunc%ez_membrane_sigma(1:var_size) = 0.0_wp
      enefunc%ez_membrane_polym(1:var_size) = 0.0_wp
      enefunc%ez_membrane_func(1:var_size)  = EzMembraneNo


    case(EneFuncNonbGO)

      if (allocated(enefunc%nonb_eps)) then
        if (size(enefunc%nonb_eps(:)) == var_size) return
        deallocate(enefunc%nonb_eps,      &
                   enefunc%nonb_rmin,     &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%nonb_eps(var_size),      &
               enefunc%nonb_rmin(var_size),     &
               stat = alloc_stat)

      enefunc%nonb_eps     (1:var_size) = 0.0_wp
      enefunc%nonb_rmin    (1:var_size) = 0.0_wp

    case(EneFuncMultiCntc)

      if (allocated(enefunc%multi_contact_list)) then
        if (size(enefunc%multi_contact_list(1,:)) == var_size) return
        deallocate(enefunc%multi_contact_list,  &
                   enefunc%multi_contact_model, &
                   enefunc%multi_contact_lj6,   &
                   enefunc%multi_contact_lj10,  &
                   enefunc%multi_contact_lj12,  &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%multi_contact_list(2,var_size),   &
               enefunc%multi_contact_model(var_size),    &
               enefunc%multi_contact_lj6(var_size),      &
               enefunc%multi_contact_lj10(var_size),     &
               enefunc%multi_contact_lj12(var_size),     &
               stat = alloc_stat)

      enefunc%multi_contact_list(1:2,1:var_size) = 0
      enefunc%multi_contact_model(1:var_size)    = 1
      enefunc%multi_contact_lj6 (1:var_size)     = 0.0_wp
      enefunc%multi_contact_lj10(1:var_size)     = 0.0_wp
      enefunc%multi_contact_lj12(1:var_size)     = 0.0_wp

    case(EneFuncMultiWork)

      if (allocated(enefunc%force_mb_work)) then
        if (size(enefunc%energy_mb_work(1,:)) == var_size) return
        deallocate(enefunc%force_mb_work,  &
                   enefunc%energy_mb_work, &
                   enefunc%virial_mb_work, &
                   enefunc%basinenergy,    &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%force_mb_work(3, var_size2, var_size),   &
               enefunc%energy_mb_work(10, var_size),    &
               enefunc%virial_mb_work(3, 3, var_size),  &
               enefunc%basinenergy(var_size),           &
               stat = alloc_stat)

      enefunc%force_mb_work(1:3,1:var_size2,1:var_size) = 0.0_wp
      enefunc%energy_mb_work(1:10,1:var_size) = 0.0_wp
      enefunc%virial_mb_work(1:3,1:3,1:var_size) = 0.0_wp
      enefunc%basinenergy(1:var_size) = 0.0_wp

    case(EneFuncMorph)

      if (allocated(enefunc%morph_dist_bb)) then
        if (size(enefunc%morph_dist_bb) == var_size) return
        deallocate(enefunc%morph_dist_bb, &
                   enefunc%morph_list_bb, &
                   enefunc%morph_dist_bb_other, &
                   stat = dealloc_stat)
      end if
      if (var_size2 > 0) then
        if (allocated(enefunc%morph_dist_sc)) then
          if (size(enefunc%morph_dist_sc) /= var_size2) then
            deallocate(enefunc%morph_dist_sc, &
                       enefunc%morph_list_sc, &
                       enefunc%morph_dist_sc_other, &
                       stat = dealloc_stat)
          end if
        end if
      end if

      allocate(enefunc%morph_dist_bb(1:var_size), &
               enefunc%morph_list_bb(1:2, 1:var_size), &
               enefunc%morph_dist_bb_other(1:var_size), &
               stat = alloc_stat)

      if (var_size2 > 0) then
        allocate(enefunc%morph_dist_sc(1:var_size2), &
                 enefunc%morph_list_sc(1:2, 1:var_size2), &
                 enefunc%morph_dist_sc_other(1:var_size2), &
                 stat = alloc_stat)
         enefunc%morph_dist_sc(1:var_size2) = 0.0_wp
         enefunc%morph_dist_sc_other(1:var_size2) = 0.0_wp
         enefunc%morph_list_sc(1:2,1:var_size2) = 0
      end if
      enefunc%morph_dist_bb(1:var_size) = 0.0_wp
      enefunc%morph_dist_bb_other(1:var_size) = 0.0_wp
      enefunc%morph_list_bb(1:2,1:var_size) = 0

    case(EneFuncGamdDih)

      if (allocated(enefunc%gamd%f_dihe)) then
        if (size(enefunc%gamd%f_dihe(1,:)) == var_size) return
        deallocate(enefunc%gamd%f_dihe, stat = dealloc_stat)
      end if

      allocate(enefunc%gamd%f_dihe(3,var_size), stat = alloc_stat)

      enefunc%gamd%f_dihe(1:3,1:var_size)  = 0.0_wp

    case(EneFuncGamdRest)

      if (allocated(enefunc%gamd%f_rest)) then
        if (size(enefunc%gamd%f_rest(1,:)) == var_size) return
        deallocate(enefunc%gamd%f_rest, stat = dealloc_stat)
      end if

      allocate(enefunc%gamd%f_rest(3,var_size), stat = alloc_stat)

      enefunc%gamd%f_rest(1:3,1:var_size)  = 0.0_wp

    case default

      call error_msg('Alloc_Enefunc> bad variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc    

    return

  end subroutine alloc_enefunc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_enefunc
  !> @brief        deallocate energy functions information
  !! @authors      YS, TY, CK
  !! @param[inout] enefunc  : structure of EneFunc information
  !! @param[in]    variable : selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_enefunc(enefunc, variable)

    ! formal arguments
    type(s_enefunc),         intent(inout) :: enefunc
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat

  
    dealloc_stat = 0

    select case (variable)

    case(EneFuncBond)
      
      if (allocated(enefunc%bond_list)) then
        deallocate(enefunc%bond_list,        &
                   enefunc%bond_force_const, &
                   enefunc%bond_dist_min,    &
                   stat = dealloc_stat)
      end if


    ! ~CG~ 3SPN.2C DNA: quartic-form bond potential
    case(EneFuncBondQuartic)

      if (allocated(enefunc%bond_quartic_list)) then
         deallocate(enefunc%bond_quartic_list,        &
              enefunc%bond_quartic_force_const, &
              enefunc%bond_quartic_dist_min,    &
              stat = dealloc_stat)
      end if

    case(EneFuncAngl)

      if (allocated(enefunc%angl_list)) then
        deallocate(enefunc%angl_list,           &
                   enefunc%angl_force_const,    &
                   enefunc%angl_theta_min,      &
                   enefunc%angl_func,           &
                   enefunc%angl_force_const1,   &
                   enefunc%angl_theta_min1,     &
                   enefunc%angl_gamma,          &
                   enefunc%angl_epsa,           &
                   !shinobu-edited
                   enefunc%angl_w,              &
                   stat = dealloc_stat)
      end if

    case(EneFuncUrey)

      if (allocated(enefunc%urey_list)) then
        deallocate(enefunc%urey_list,        &
                   enefunc%urey_force_const, &
                   enefunc%urey_rmin,        &
                   stat = dealloc_stat)
      end if

    !shinobu-edited
    case(EneFuncAngFlex)
      if (allocated(enefunc%anglflex_list)) then
        deallocate(enefunc%anglflex_list,     &
                   enefunc%anglflex_type,     &
                   stat = dealloc_stat)
      end if

    case(EneFuncAngFlexTbl)
      if (allocated(enefunc%anglflex_theta)) then
        deallocate(enefunc%anglflex_theta,            &
                   enefunc%anglflex_efunc,            &
                   enefunc%anglflex_d2func,           &
                   enefunc%anglflex_min_th,           &
                   enefunc%anglflex_max_th,           &
                   enefunc%anglflex_ener_corr,        &
                   stat = dealloc_stat)
      end if

    case(EneFuncDihe)

      if (allocated(enefunc%dihe_list)) then
        deallocate(enefunc%dihe_list,        &
                   enefunc%dihe_force_const, &
                   enefunc%dihe_periodicity, &
                   enefunc%dihe_phase,       &
                   enefunc%dihe_scee,        &
                   enefunc%dihe_scnb,        &
                   !shinobu-edited
                   enefunc%dihe_theta_min,   &
                   enefunc%dihe_w,           &
                   ! ~CG~ 3SPN.2C DNA: dihedral functype
                   enefunc%dihe_func,        &
                   stat = dealloc_stat)
      end if

    case(EneFuncDiheFlex)
      if (allocated(enefunc%diheflex_list)) then
        deallocate(enefunc%diheflex_list,      & 
                   enefunc%diheflex_type,      &
                   enefunc%diheflex_func,      &
                   stat = dealloc_stat)
      end if

    case(EneFuncDiheFlexTbl)
      if (allocated(enefunc%diheflex_coef)) then
        deallocate(enefunc%diheflex_coef,            &
                   enefunc%diheflex_ener_corr,       &
                   stat = dealloc_stat)
      end if

    case(EneFuncRBDihe)

      if (allocated(enefunc%rb_dihe_list)) then
        deallocate(enefunc%rb_dihe_list,     &
                   enefunc%rb_dihe_c,        &
                   stat = dealloc_stat)
      end if

    case(EneFuncImpr)

      if (allocated(enefunc%impr_list)) then
        deallocate(enefunc%impr_list,        &
                   enefunc%impr_force_const, &
                   enefunc%impr_periodicity, &
                   enefunc%impr_phase,       &
                   stat = dealloc_stat)
      end if

    case(EneFuncCmap)

      if (allocated(enefunc%cmap_list)) then
        deallocate(enefunc%cmap_list,        &
                   enefunc%cmap_resolution,  &
                   enefunc%cmap_coef,        &
                   enefunc%cmap_force,       &
                   stat = dealloc_stat )
      end if
      
    case(EneFuncCGGeneral)
      if (allocated(enefunc%NA_base_type)) then
        deallocate(enefunc%NA_base_type,  &
            stat = dealloc_stat)
      end if
      if (allocated(enefunc%mol_chain_id)) then
        deallocate(enefunc%mol_chain_id,  &
            stat = dealloc_stat)
      end if
      if (allocated(enefunc%atom_cls)) then
        deallocate(enefunc%atom_cls,  &
            stat = dealloc_stat)
      end if
      if (allocated(enefunc%param_epsilon)) then
        deallocate(enefunc%param_epsilon,  &
            stat = dealloc_stat)
      end if
      if (allocated(enefunc%param_sigma)) then
        deallocate(enefunc%param_sigma,  &
            stat = dealloc_stat)
      end if

    case(EneFuncBaseStack)
      ! ~CG~ 3SPN.2C DNA: bond stacking
      if (allocated(enefunc%base_stack_list)) then
        deallocate(enefunc%base_stack_list,     &
                   enefunc%base_stack_func,     &
                   enefunc%base_stack_epsilon,  &
                   enefunc%base_stack_sigma,    &
                   enefunc%base_stack_theta_bs, &
                   stat = dealloc_stat)
      end if

    case(EneFuncBasePair)
      ! ~CG~ 3SPN.2C DNA: bond pairing
      if (allocated(enefunc%base_pair_theta_1)) then
        deallocate(enefunc%base_pair_theta_1,     &
                   enefunc%base_pair_theta_2,     &
                   enefunc%base_pair_theta_3,     &
                   enefunc%base_pair_phi_1,       &
                   enefunc%base_pair_sigma,       &
                   enefunc%base_pair_epsilon,     &
                   enefunc%base_cross_1_epsilon,  &
                   enefunc%base_cross_1_sigma,    &
                   enefunc%base_cross_1_theta_cs, &
                   enefunc%base_cross_2_epsilon,  &
                   enefunc%base_cross_2_sigma,    &
                   enefunc%base_cross_2_theta_cs, &
                   enefunc%base_pair_is_WC,       &
                   stat = dealloc_stat)
      end if
      
      if (allocated(enefunc%cg_particle_DNA_all)) then
        deallocate(enefunc%cg_particle_DNA_all,  &
            stat = dealloc_stat)
      end if
      if (allocated(enefunc%cg_particle_DNA_base)) then
        deallocate(enefunc%cg_particle_DNA_base, &
            stat = dealloc_stat)
      end if
      if (allocated(enefunc%cg_particle_DNA_phos)) then
        deallocate(enefunc%cg_particle_DNA_phos, &
            stat = dealloc_stat)
      end if

    case(EneFuncCGDNAExv)

      if (allocated(enefunc%cgDNA_exv_sigma)) then
        deallocate(enefunc%cgDNA_exv_sigma,     &
                   stat = dealloc_stat)
      end if

    case(EneFuncCGele)

      if (allocated(enefunc%cg_charge)) then
        deallocate(enefunc%cg_charge,           &
            stat = dealloc_stat)
      end if
      if (allocated(enefunc%cg_ele_mol_pair)) then
        deallocate(enefunc%cg_ele_mol_pair,     &
            stat = dealloc_stat)
      end if
      if (allocated(enefunc%cg_particle_charged)) then
        deallocate(enefunc%cg_particle_charged, &
            stat = dealloc_stat)
      end if
      if (allocated(enefunc%cg_particle_pro_charged)) then
        deallocate(enefunc%cg_particle_pro_charged, &
            stat = dealloc_stat)
      end if

    case(EneFuncCGKH)

      if (allocated(enefunc%cg_KH_sigma_half)) then
        deallocate(enefunc%cg_KH_sigma_half, &
            stat = dealloc_stat)
      end if

      if (allocated(enefunc%cg_pro_use_KH)) then
        deallocate(enefunc%cg_pro_use_KH, &
            stat = dealloc_stat)
      end if

      if (allocated(enefunc%cg_particle_KH)) then
        deallocate(enefunc%cg_particle_KH, &
            stat = dealloc_stat)
      end if

      if (allocated(enefunc%cg_KH_epsilon)) then
        deallocate(enefunc%cg_KH_epsilon,    &
            stat = dealloc_stat)
      end if

    case(EneFuncCGKHmol)

      if (allocated(enefunc%cg_KH_mol_pair)) then
        deallocate(enefunc%cg_KH_mol_pair,   &
            stat = dealloc_stat)
      end if

    case(EneFuncPWMcos)
      ! ~CG~ protein-DNA interactions with PWMcos model
      if (allocated(enefunc%pwmcos_protein_id)) then
        deallocate(enefunc%pwmcos_protein_id, &
            enefunc%pwmcos_protein_id_N,      &
            enefunc%pwmcos_protein_id_C,      &
            enefunc%pwmcos_r0,                &
            enefunc%pwmcos_theta1,            &
            enefunc%pwmcos_theta2,            &
            enefunc%pwmcos_theta3,            &
            enefunc%pwmcos_ene_A,             &
            enefunc%pwmcos_ene_C,             &
            enefunc%pwmcos_ene_G,             &
            enefunc%pwmcos_ene_T,             &
            enefunc%pwmcos_gamma,             &
            enefunc%pwmcos_eps,               &
            enefunc%pwmcos_specificity,       &
            enefunc%pwmcos_to_pairlist_id,    &
            stat = dealloc_stat)
      end if
      if (allocated(enefunc%pwmcos_mol_pair)) then
        deallocate(enefunc%pwmcos_mol_pair,     &
            stat = dealloc_stat)
      end if
      if (allocated(enefunc%pwmcos_involved_resid)) then
        deallocate(enefunc%pwmcos_involved_resid,     &
            enefunc%pwmcos_involved_spec,             &
            stat = dealloc_stat)
      end if
      
    case(EneFuncPWMcosns)
      ! ~CG~ protein-DNA interactions with PWMcosns model
      if (allocated(enefunc%pwmcosns_protein_id)) then
        deallocate(enefunc%pwmcosns_protein_id, &
            enefunc%pwmcosns_protein_id_N,      &
            enefunc%pwmcosns_protein_id_C,      &
            enefunc%pwmcosns_r0,                &
            enefunc%pwmcosns_theta1,            &
            enefunc%pwmcosns_theta2,            &
            enefunc%pwmcosns_ene,               &
            enefunc%pwmcosns_specificity,       &
            enefunc%pwmcosns_to_pairlist_id,    &
            stat = dealloc_stat)
      end if
      if (allocated(enefunc%pwmcosns_mol_pair)) then
        deallocate(enefunc%pwmcosns_mol_pair,     &
            stat = dealloc_stat)
      end if
      if (allocated(enefunc%pwmcosns_involved_resid)) then
        deallocate(enefunc%pwmcosns_involved_resid,     &
            enefunc%pwmcosns_involved_spec,             &
            stat = dealloc_stat)
      end if
 
    case(EneFuncCGIDRHPS)

      if (allocated(enefunc%cg_IDR_HPS_is_IDR)) then
        deallocate(enefunc%cg_IDR_HPS_is_IDR, &
            enefunc%cg_IDR_HPS_lambda_half,   &
            enefunc%cg_IDR_HPS_sigma_half,    &
            stat = dealloc_stat)
      end if
      if (allocated(enefunc%cg_particle_IDR_HPS)) then
        deallocate(enefunc%cg_particle_IDR_HPS, &
            stat = dealloc_stat)
      end if

    case(EneFuncCGIDRKH)

      if (allocated(enefunc%cg_IDR_KH_is_IDR)) then
        deallocate(enefunc%cg_IDR_KH_is_IDR, &
            enefunc%cg_IDR_KH_sigma_half,    &
            stat = dealloc_stat)
      end if
      if (allocated(enefunc%cg_IDR_KH_epsilon_D)) then
        deallocate(enefunc%cg_IDR_KH_epsilon_D, &
            stat = dealloc_stat)
      end if
      if (allocated(enefunc%cg_particle_IDR_KH)) then
        deallocate(enefunc%cg_particle_IDR_KH, &
            stat = dealloc_stat)
      end if

    case(EneFuncNbon)

      if (allocated(enefunc%nonb_atom_cls)) then
        deallocate(enefunc%nonb_atom_cls, &
                   enefunc%nb14_lj6,      &
                   enefunc%nb14_lj10,     &
                   enefunc%nb14_lj12,     &
                   enefunc%nonb_lj6,      &
                   enefunc%nonb_lj10,     &
                   enefunc%nonb_lj12,     &
                   enefunc%nonb_aicg_sig, &
                   enefunc%nonb_aicg_eps, &
                   enefunc%cg_exv_eps_sqrt, &
                   enefunc%cg_exv_sig_half, &
                   stat = dealloc_stat)
      end if

    case(EneFuncCntc)

      if (allocated(enefunc%contact_list)) then
        deallocate(enefunc%contact_list,       &
                   enefunc%contact_func,       &
                   enefunc%contact_lj6,        &
                   enefunc%contact_lj10,       &
                   enefunc%contact_lj12,       &
                   stat = dealloc_stat)
      end if

    case(EneFuncVsite2)

      if (allocated(enefunc%vsite2_list)) then
        deallocate(enefunc%vsite2_list, &
                   enefunc%vsite2_a,    &
                   stat = dealloc_stat)
      end if

    case(EneFuncVsite3)

      if (allocated(enefunc%vsite3_list)) then
        deallocate(enefunc%vsite3_list, &
                   enefunc%vsite3_a,    &
                   enefunc%vsite3_b,    &
                   stat = dealloc_stat)
      end if

    case(EneFuncVsite3fd)

      if (allocated(enefunc%vsite3fd_list)) then
        deallocate(enefunc%vsite3fd_list, &
                   enefunc%vsite3fd_a,    &
                   enefunc%vsite3fd_d,    &
                   stat = dealloc_stat)
      end if

    case(EneFuncVsite3fad)

      if (allocated(enefunc%vsite3fad_list)) then
        deallocate(enefunc%vsite3fad_list,  &
                   enefunc%vsite3fad_theta, &
                   enefunc%vsite3fad_d,     &
                   stat = dealloc_stat)
      end if

    case(EneFuncVsite3out)

      if (allocated(enefunc%vsite3out_list)) then
        deallocate(enefunc%vsite3out_list, &
                   enefunc%vsite3out_a,    &
                   enefunc%vsite3out_b,    &
                   enefunc%vsite3out_c,    &
                   stat = dealloc_stat)
      end if

    case(EneFuncVsite4fdn)

      if (allocated(enefunc%vsite4fdn_list)) then
        deallocate(enefunc%vsite4fdn_list, &
                   enefunc%vsite4fdn_a,    &
                   enefunc%vsite4fdn_b,    &
                   enefunc%vsite4fdn_c,    &
                   stat = dealloc_stat)
      end if

    case(EneFuncVsiten)

      if (allocated(enefunc%vsiten_list)) then
        deallocate(enefunc%vsiten_list, &
                   enefunc%vsiten_n,    &
                   stat = dealloc_stat)
      end if

    case(EneFuncRefg)

      if (allocated(enefunc%restraint_numatoms)) then
        deallocate(enefunc%restraint_numatoms,      &
                   enefunc%restraint_atomlist,      &
                   enefunc%restraint_masscoef,      &
                   enefunc%restraint_totmass_group, &
                   enefunc%restraint_wcom3,         &
                   enefunc%restraint_wcom4,         &
                   enefunc%restraint_wcom5,         &
                   enefunc%restraint_wtmp,          &
                   enefunc%restraint_istartend, &
                   stat = dealloc_stat)
      end if

    case(EneFuncReff)

      if (allocated(enefunc%restraint_kind)) then
        deallocate(enefunc%restraint_kind,           &
                   enefunc%restraint_grouplist,      &
                   enefunc%restraint_const,          &
                   enefunc%restraint_ref,            &
                   enefunc%restraint_funcgrp,        &
                   enefunc%restraint_exponent_func,  &
                   enefunc%restraint_exponent_dist,  &
                   enefunc%restraint_mode,           &
                   enefunc%restraint_replica_index,  &
                   enefunc%restraint_caging,         &
                   enefunc%restraint_flat_radius,    &
                   enefunc%restraint_weight_dist,    &
                   enefunc%restraint_wcom1,          &
                   enefunc%restraint_wcom2,          &
                   enefunc%restraint_wdrt,           &
                   enefunc%restraint_rgbuffer,       &
                   enefunc%restraint_wall_z,         &
                   enefunc%restraint_rpath_func,     &
                   enefunc%restraint_rcom1,          &
                   enefunc%restraint_rcom2,          &
                   enefunc%restraint_rdrt,           &
                   stat = dealloc_stat)
      end if

    case(EneFuncRefw)

      if (allocated(enefunc%restraint_wfrc)) then
        deallocate(enefunc%restraint_wfrc, &
                   stat = dealloc_stat)
      end if

    case(EneFuncRefc)

      if (allocated(enefunc%restraint_refcoord)) then
        deallocate(enefunc%restraint_refcoord, &
                   enefunc%restraint_masstmp,  &
                   stat = dealloc_stat)
      end if

    case(EneFuncRefr)

      if (allocated(enefunc%restraint_const_replica)) then
        deallocate(enefunc%restraint_const_replica, &
                   enefunc%restraint_ref_replica, &
                   stat = dealloc_stat)
      end if

    case(EneFuncMode)

      if (allocated(enefunc%pc_mode)) then
        deallocate(enefunc%pc_mode, &
                   stat = dealloc_stat)
      end if

    case(EneFuncTable)

      if (allocated(enefunc%table%table_ene_WW)) then
        deallocate(enefunc%table%table_ene,    &
                   enefunc%table%table_grad,   &
                   enefunc%table%table_ecor,   &
                   enefunc%table%table_decor,  &
                   enefunc%table%dtable_ecor,  &
                   enefunc%table%dtable_decor, &
                   enefunc%table%table_ene_WW, &
                   enefunc%table%table_de_WW,  &
                   stat = dealloc_stat)
      end if

    case(EneFuncTableWater)

      if (allocated(enefunc%table%water_list)) then
        deallocate(enefunc%table%water_list, &
                   stat = dealloc_stat)
      end if

    case(EneFuncTableSolute)

      if (allocated(enefunc%table%solute_list)) then
        deallocate(enefunc%table%solute_list,     &
                   enefunc%table%solute_list_inv, &
                   stat = dealloc_stat)
      end if

    case(EneFuncCmapType)

      if (allocated(enefunc%cmap_resolution)) then
        deallocate(enefunc%cmap_resolution, &
                   enefunc%cmap_coef,       &
                   stat = dealloc_stat)
      end if

    case(EneFuncPme) 

      if (allocated(enefunc%pme%gx)) then
        deallocate(enefunc%pme%gx,           &
                   enefunc%pme%gy,           &
                   enefunc%pme%gz,           &
                   enefunc%pme%b2,           &
                   enefunc%pme%vir_fact,     &
                   enefunc%pme%theta,        &
                   enefunc%pme%ftqdf_localA, &
                   enefunc%pme%ftqdf_localB, &
                   enefunc%pme%qdf,          &
                   enefunc%pme%qdf_recv,     &
                   enefunc%pme%qdf_send,     &
                   enefunc%pme%f,            &
                   enefunc%pme%v,            &
                   enefunc%pme%i_list,       &
                   enefunc%pme%force_wk,     &
                   enefunc%pme%x_start1,     &
                   enefunc%pme%x_end1,       &
                   enefunc%pme%y_start,      &
                   enefunc%pme%y_end,        &
                   enefunc%pme%y_local,      &
                   enefunc%pme%y_start1,     &
                   enefunc%pme%y_end1,       &
                   enefunc%pme%z_start,      &
                   enefunc%pme%z_end,        &
                   enefunc%pme%z_local,      &
                   stat = dealloc_stat)
      end if

      if (dealloc_stat /= 0) call error_msg_dealloc

      if (allocated(enefunc%pme%i_list)) then
        deallocate( enefunc%pme%i_list, &
                    enefunc%pme%force_wk, &
                    stat = dealloc_stat)
      end if


    case(EneFuncFitc) 

      if (allocated(enefunc%fit_refcoord)) then

        deallocate(enefunc%fit_refcoord,    &
                   enefunc%fit_work,        &
                   enefunc%fit_mass,        &
                   enefunc%fitting_atom,    &
                   stat = dealloc_stat)
      end if

    case(EneFuncQmmm)

      if (allocated(enefunc%qmmm%qm_force)) then
        deallocate(enefunc%qmmm%qm_force,           &
                   enefunc%qmmm%mm_force,           &
                   enefunc%qmmm%qmatom_id,          &
                   enefunc%qmmm%mmatom_id,          &
                   enefunc%qmmm%qm_atomic_no,       &
                   enefunc%qmmm%qm_charge,          &
                   stat = dealloc_stat)

        if(enefunc%qmmm%num_qmmmbonds /= 0) then
          deallocate(enefunc%qmmm%linkatom_coord,     &
                     enefunc%qmmm%linkatom_force,     &
                     enefunc%qmmm%linkatom_charge,    &
                     enefunc%qmmm%linkatom_global_address, &
                     enefunc%qmmm%qmmmbond_list,      &
                     enefunc%qmmm%ecatom_id,          &
                     stat = dealloc_stat)
        end if

#ifdef QSIMULATE
        call enefunc%qmmm%qs_input%destroy()
        call enefunc%qmmm%qs_output%destroy()
#endif

      endif

    case(EneFuncModeRef)

      if (allocated(enefunc%pc_crd)) then
        deallocate(enefunc%pc_crd,         &
                   enefunc%pc_crd_ref,     &
                   enefunc%pc_crd_ref_fit, &
                   enefunc%pc_mode_temp,   &
                   enefunc%pc_mode_fit,    &
                   enefunc%grad_pc,        &
                   stat = dealloc_stat)
      end if

    case(EneFuncEef1)

      if (allocated(enefunc%eef1%atom_name)) then
        deallocate(enefunc%eef1%atom_name,  &
                   enefunc%eef1%volume,     &
                   enefunc%eef1%gref_0,     &
                   enefunc%eef1%gfree_0,    &
                   enefunc%eef1%href_0,     &
                   enefunc%eef1%cpref_0,    &
                   enefunc%eef1%gref_t,     &
                   enefunc%eef1%gfree_t,    &
                   enefunc%eef1%inv_lambda, &
                   enefunc%eef1%rvdw,       &
                   enefunc%eef1%alpha_4pi,  &
                   stat = dealloc_stat)
      end if

    case(EneFuncGbsa)

      if (allocated(enefunc%gbsa%vdw_radius)) then
        deallocate(enefunc%gbsa%vdw_radius,      &
                   enefunc%gbsa%scale_factor,    &
                   enefunc%gbsa%sasa_atom_type,  &
                   enefunc%gbsa%sasa_parameters, &
                   enefunc%gbsa%sasa_vdw_radius, &
                   enefunc%gbsa%sasa_atom_list,  &
                   enefunc%gbsa%born_radius,     &
                   enefunc%gbsa%dadrfac,         &
                   stat = dealloc_stat)
      end if

    case(EneFuncNonbGO)

      if (allocated(enefunc%nonb_eps)) then
        deallocate(enefunc%nonb_eps,      &
                   enefunc%nonb_rmin,     &
                   stat = dealloc_stat)
      end if

    case(EneFuncMultiCntc)

      if (allocated(enefunc%multi_contact_list)) then
        deallocate(enefunc%multi_contact_list,  &
                   enefunc%multi_contact_model, &
                   enefunc%multi_contact_lj6,   &
                   enefunc%multi_contact_lj10,  &
                   enefunc%multi_contact_lj12,  &
                   stat = dealloc_stat)
      end if

    case(EneFuncMultiWork)

      if (allocated(enefunc%force_mb_work)) then
        deallocate(enefunc%force_mb_work,  &
                   enefunc%energy_mb_work, &
                   enefunc%virial_mb_work,   &
                   stat = dealloc_stat)
      end if

    case(EneFuncMorph)

      if (allocated(enefunc%morph_dist_bb)) then
        deallocate(enefunc%morph_dist_bb, &
                   enefunc%morph_list_bb, &
                   enefunc%morph_dist_bb_other, &
                   stat = dealloc_stat)
      end if
      if (allocated(enefunc%morph_dist_sc)) then
        deallocate(enefunc%morph_dist_sc, &
                   enefunc%morph_list_sc, &
                   enefunc%morph_dist_sc_other, &
                   stat = dealloc_stat)
      end if

    case(EneFuncEzMembrane)

      if (allocated(enefunc%ez_membrane_const)) then
        deallocate(enefunc%ez_membrane_const, &
                   enefunc%ez_membrane_zmin, &
                   enefunc%ez_membrane_sigma, &
                   enefunc%ez_membrane_polym, &
                   enefunc%ez_membrane_func,  &
                   stat = dealloc_stat)
      end if

    case(EneFuncGamdDih)

      if (allocated(enefunc%gamd%f_dihe)) then
        deallocate(enefunc%gamd%f_dihe, stat = dealloc_stat)
      end if

    case(EneFuncGamdRest)

      if (allocated(enefunc%gamd%f_rest)) then
        deallocate(enefunc%gamd%f_rest, stat = dealloc_stat)
      end if

    case default

      call error_msg('Dealloc_Enefunc> bad variable')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_enefunc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_enefunc_all
  !> @brief        deallocate all energy functions information
  !! @authors      YS, CK
  !! @param[inout] enefunc : structure of potential energy function
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_enefunc_all(enefunc)

    ! formal arguments
    type(s_enefunc),         intent(inout) :: enefunc


    call dealloc_enefunc(enefunc, EneFuncBond)
    call dealloc_enefunc(enefunc, EneFuncAngl)
    call dealloc_enefunc(enefunc, EneFuncUrey)
    !shinobu-edited
    call dealloc_enefunc(enefunc, EneFuncAngFlex)
    call dealloc_enefunc(enefunc, EneFuncAngFlexTbl)
    call dealloc_enefunc(enefunc, EneFuncDihe)
    call dealloc_enefunc(enefunc, EneFuncDiheFlex)
    call dealloc_enefunc(enefunc, EneFuncRBDihe)
    call dealloc_enefunc(enefunc, EneFuncImpr)
    call dealloc_enefunc(enefunc, EneFuncCmap)
    call dealloc_enefunc(enefunc, EneFuncNbon)
    call dealloc_enefunc(enefunc, EneFuncCntc)
    call dealloc_enefunc(enefunc, EneFuncVsite2)
    call dealloc_enefunc(enefunc, EneFuncVsite3)
    call dealloc_enefunc(enefunc, EneFuncVsite3fd)
    call dealloc_enefunc(enefunc, EneFuncVsite3fad)
    call dealloc_enefunc(enefunc, EneFuncVsite3out)
    call dealloc_enefunc(enefunc, EneFuncVsite4fdn)
    call dealloc_enefunc(enefunc, EneFuncVsiten)
    call dealloc_enefunc(enefunc, EneFuncRefg)
    call dealloc_enefunc(enefunc, EneFuncReff)
    call dealloc_enefunc(enefunc, EneFuncRefw)
    call dealloc_enefunc(enefunc, EneFuncRefc)
    call dealloc_enefunc(enefunc, EneFuncRefr)
    call dealloc_enefunc(enefunc, EneFuncTable)
    call dealloc_enefunc(enefunc, EneFuncTableWater)
    call dealloc_enefunc(enefunc, EneFuncTableSolute)
    call dealloc_enefunc(enefunc, EneFuncCmapType)
    call dealloc_enefunc(enefunc, EneFuncPme)
    call dealloc_enefunc(enefunc, EneFuncMode)
    call dealloc_enefunc(enefunc, EneFuncFitc)
    call dealloc_enefunc(enefunc, EneFuncQmmm)
    call dealloc_enefunc(enefunc, EneFuncModeRef)
    call dealloc_enefunc(enefunc, EneFuncEef1)
    call dealloc_enefunc(enefunc, EneFuncGbsa)
    call dealloc_enefunc(enefunc, EneFuncEzMembrane)
    call dealloc_enefunc(enefunc, EneFuncNonbGO)
    call dealloc_enefunc(enefunc, EneFuncMultiCntc)
    call dealloc_enefunc(enefunc, EneFuncMultiWork)
    call dealloc_enefunc(enefunc, EneFuncMorph)
    ! ~CG~ 3SPN.2C DNA: quartic bond potential
    call dealloc_enefunc(enefunc, EneFuncBondQuartic)
    call dealloc_enefunc(enefunc, EneFuncBaseStack)
    call dealloc_enefunc(enefunc, EneFuncBasePair)
    call dealloc_enefunc(enefunc, EneFuncCGDNAExv)
    call dealloc_enefunc(enefunc, EneFuncCGele)
    call dealloc_enefunc(enefunc, EneFuncCGKH)
    call dealloc_enefunc(enefunc, EneFuncPWMcos)
    call dealloc_enefunc(enefunc, EneFuncPWMcosns)
    call dealloc_enefunc(enefunc, EneFuncGamdDih)
    call dealloc_enefunc(enefunc, EneFuncGamdRest)
    call dealloc_enefunc(enefunc, EneFuncCGIDRHPS)
    call dealloc_enefunc(enefunc, EneFuncCGIDRKH)

    return

  end subroutine dealloc_enefunc_all

end module at_enefunc_str_mod
