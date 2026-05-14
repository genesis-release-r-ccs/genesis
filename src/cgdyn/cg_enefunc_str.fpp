!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_enefunc_str_mod
!> @brief   Structure of energy functions.
!! @authors Yuji Sugita (YS), Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_enefunc_str_mod

  use cg_domain_str_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_enefunc

    integer                       :: forcefield
    integer                       :: electrostatic
    integer                       :: output_style

    integer                       :: num_bonds
    integer                       :: num_bonds_quartic
    integer                       :: num_angles
    integer                       :: num_angflex
    integer                       :: num_anglocal
    integer                       :: num_dihe_all
    integer                       :: num_base_stack_all
    integer                       :: num_contact_all
    integer                       :: num_excl
    integer                       :: num_nb14_all

    integer                       :: num_atom_cls
    integer                       :: num_atoms_ref
    integer                       :: num_restraintgroups
    integer                       :: num_restraintfuncs
    integer                       :: max_restraint_numatoms
    integer                       :: max_restraint_numgrps

    logical                       :: assign_force_max
    real(wp)                      :: upper_force_value

    ! AICG
    logical                       :: cg_safe_dihedral_calc
    logical                       :: cg_ele_calc
    logical                       :: cg_DNA_base_pair_calc
    logical                       :: cg_DNA_exv_calc
    logical                       :: cg_pwmcos_calc
    logical                       :: cg_pwmcosns_calc
    logical                       :: cg_KH_calc
    logical                       :: cg_IDR_KH_calc
    logical                       :: cg_IDR_HPS_calc
    logical                       :: cg_IDR_MPIPI_calc
    integer                       :: notation_14types
    integer,          allocatable :: atom_cls_2_base_type(:) 
    real(wp),         allocatable :: dihe_scnb(:), dihe_scee(:)

    ! bond (size = num_local_bond for each cell)
    integer                       :: num_bondsq_domain
    integer                       :: num_bond_domain
    integer,          allocatable :: bond_list(:,:)
    real(wp),         allocatable :: bond_force_const(:)
    real(wp),         allocatable :: bond_dist_min(:)
    integer,          allocatable :: bond_kind(:)

    ! angle (size = num_local_angle for each cell)
    integer                       :: num_angle_flexible_domain
    integer                       :: num_angle_local_domain
    integer                       :: num_angle_domain
    integer,          allocatable :: angl_list(:,:)
    real(wp),         allocatable :: angl_force_const(:)
    real(wp),         allocatable :: angl_theta_min(:)
    real(wp),         allocatable :: urey_force_const(:)
    real(wp),         allocatable :: urey_rmin(:)
    real(wp),         allocatable :: angl_width(:)
    integer,          allocatable :: angl_kind(:)

    ! for flexible angle
    real(wp),         allocatable :: anglflex_theta(:,:)
    real(wp),         allocatable :: anglflex_efunc(:,:)
    real(wp),         allocatable :: anglflex_d2func(:,:)
    real(wp),         allocatable :: anglflex_min_th(:,:)
    real(wp),         allocatable :: anglflex_max_th(:,:)
    real(wp),         allocatable :: anglflex_ener_corr(:)

    ! dihedral (size = num_local_dihedral for each cell)
    integer                       :: num_dihe_flexible_domain
    integer                       :: num_dihe_local_domain
    integer                       :: num_dihe_domain
    integer,          allocatable :: dihe_list(:,:)
    real(wp),         allocatable :: dihe_force_const(:)
    integer,          allocatable :: dihe_periodicity(:)
    real(wp),         allocatable :: dihe_phase(:)
    real(wp),         allocatable :: dihe_width(:)
    integer,          allocatable :: dihe_kind(:)

    ! for flexible dihedral angle
    !
    real(wp),         allocatable :: diheflex_coef(:,:)
    real(wp),         allocatable :: diheflex_ener_corr(:)
    real(wp)                      :: cg_safe_dih_ene_shift

    ! base stacking
    !
    integer                       :: num_stack_domain = 0
    integer,          allocatable :: base_stack_list(:,:)
    integer,          allocatable :: base_stack_func(:)
    real(wp),         allocatable :: base_stack_epsilon(:)
    real(wp),         allocatable :: base_stack_sigma(:)
    real(wp),         allocatable :: base_stack_theta_bs(:)
    real(wp)                      :: base_stack_alpha
    real(wp)                      :: base_stack_K

    ! contact (size = num_local_contact for each cell)
    integer                       :: num_contact_domain = 0
    integer                       :: num_contact_boundary = 0
    integer,          allocatable :: contact_list(:,:)
    integer,          allocatable :: contact_func(:)
    real(wp),         allocatable :: contact_lj12(:)
    real(wp),         allocatable :: contact_lj10(:)
    real(wp),         allocatable :: contact_lj6 (:)

    ! elec
    real(wp)                      :: cg_pro_DNA_ele_scale_Q
    integer(1),       allocatable :: cg_ele_mol_pair(:,:)

    ! ~CG~ PWMcos: CG protein-DNA seq-specific
    integer                       :: num_pwmcos_domain = 0
    real(wp)                      :: pwmcos_sigma
    real(wp)                      :: pwmcos_phi
    integer,          allocatable :: pwmcos_count(:)
    integer,          allocatable :: pwmcos_protein_id(:)
    integer,          allocatable :: pwmcos_protein_id_N(:)
    integer,          allocatable :: pwmcos_protein_id_C(:)
    real(wp),         allocatable :: pwmcos_r0(:,:)
    real(wp),         allocatable :: pwmcos_theta1(:,:)
    real(wp),         allocatable :: pwmcos_theta2(:,:)
    real(wp),         allocatable :: pwmcos_theta3(:,:)
    real(wp),         allocatable :: pwmcos_ene_A(:,:)
    real(wp),         allocatable :: pwmcos_ene_C(:,:)
    real(wp),         allocatable :: pwmcos_ene_G(:,:)
    real(wp),         allocatable :: pwmcos_ene_T(:,:)
    real(wp),         allocatable :: pwmcos_gamma(:,:)
    real(wp),         allocatable :: pwmcos_eps(:,:)
    integer,          allocatable :: pwmcos_specificity(:,:)
    integer(1),       allocatable :: pwmcos_mol_pair(:,:)
    integer,          allocatable :: pwmcos_involved_resid(:,:)
    integer,          allocatable :: pwmcos_involved_spec(:,:)

    ! ~CG~ PWMcos: CG protein-DNA seq-specific
    integer                       :: num_pwmcosns_domain = 0
    real(wp)                      :: pwmcosns_sigma
    real(wp)                      :: pwmcosns_phi
    integer,          allocatable :: pwmcosns_count(:)
    integer,          allocatable :: pwmcosns_protein_id(:)
    integer,          allocatable :: pwmcosns_protein_id_N(:)
    integer,          allocatable :: pwmcosns_protein_id_C(:)
    real(wp),         allocatable :: pwmcosns_r0(:,:)
    real(wp),         allocatable :: pwmcosns_theta1(:,:)
    real(wp),         allocatable :: pwmcosns_theta2(:,:)
    real(wp),         allocatable :: pwmcosns_ene(:,:)
    integer,          allocatable :: pwmcosns_specificity(:,:) 
    integer(1),       allocatable :: pwmcosns_mol_pair(:,:)
    integer,          allocatable :: pwmcosns_involved_resid(:,:)
    integer,          allocatable :: pwmcosns_involved_spec(:,:)

    ! excl
    integer,          allocatable :: excl_list(:,:)

    ! elec
    integer                       :: num_cg_elec
    real(wp)                      :: rf_self
    integer,          allocatable :: cg_elec_list(:)
    integer,          allocatable :: cg_elec_list_inv(:)

    ! base
    integer                       :: num_cg_base
    integer,          allocatable :: cg_base_list(:)
    integer,          allocatable :: cg_base_list_inv(:)

    !dna
    integer                       :: num_cg_dna
    integer,          allocatable :: DNA_end(:,:)
    integer,          allocatable :: cg_dna_list(:)
    integer,          allocatable :: cg_dna_list_inv(:)

    !KH
    integer                       :: num_cg_KH
    integer,          allocatable :: cg_KH_list(:)
    integer,          allocatable :: cg_KH_list_inv(:)
    integer(1),       allocatable :: cg_KH_mol_pair(:,:)
    real(wp),         allocatable :: cg_KH_epsilon(:,:)
    real(wp),         allocatable :: cg_KH_sigma(:,:)

    !IDR KH
    integer                       :: num_cg_IDR_KH
    integer,          allocatable :: cg_IDR_KH_list(:)
    integer,          allocatable :: cg_IDR_KH_list_inv(:)
    real(wp),         allocatable :: cg_IDR_KH_epsilon(:,:)
    real(wp),         allocatable :: cg_IDR_KH_sigma(:,:)

    !IDR HPS
    integer                       :: num_cg_IDR_HPS
    integer,          allocatable :: cg_IDR_HPS_list(:)
    integer,          allocatable :: cg_IDR_HPS_list_inv(:)
    real(wp),         allocatable :: cg_IDR_HPS_lambda(:,:)
    real(wp),         allocatable :: cg_IDR_HPS_sigma(:,:)
    real(wp)                      :: cg_IDR_HPS_epsilon

    !IDR MPIPI
    integer                       :: num_cg_IDR_MPIPI
    integer,          allocatable :: cg_IDR_MPIPI_list(:)
    integer,          allocatable :: cg_IDR_MPIPI_list_inv(:)
    ! ~CG~ : Wang-Frenkel potential
    real(wp),         allocatable :: cg_WangFrenkel_epsilon(:, :)
    real(wp),         allocatable :: cg_WangFrenkel_sigma(:, :)
    integer,          allocatable :: cg_WangFrenkel_mu(:, :)
    integer,          allocatable :: cg_WangFrenkel_nu(:, :)
    real(wp)                      :: cg_WangFrenkel_cutoff_scale

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

    ! ~CG~ 3spn.2c DNA exv (size = num_DNA_particle_types)
    real(wp),         allocatable :: cgDNA_exv_sigma(:,:)
    real(wp)                      :: cgDNA_exv_epsilon
    logical                       :: cg_infinite_DNA

    ! non-bonded (size = num_atoms)
    integer,          allocatable :: num_nonb_excl(:,:)
    integer,          allocatable :: num_nonb_excl1(:,:)
    integer,          allocatable :: num_nb14_calc(:,:)
    integer,          allocatable :: num_nb14_calc1(:,:)
    integer,          allocatable :: nonb_list(:,:,:)
    integer,          allocatable :: nonb_list1(:,:,:)
    integer,          allocatable :: nb14_list(:,:,:)
    integer,          allocatable :: nb14_list1(:,:,:)
    integer,          allocatable :: sc_list(:,:,:)
    integer,          allocatable :: sc_list1(:,:,:)
    integer,          allocatable :: num_nb14_total(:)
    integer,          allocatable :: num_nb14_total1(:)
    integer(1),       allocatable :: exclusion_mask(:,:)

    ! non-bonded list (size = num_atoms)
    real(wp)                      :: switchdist
    real(wp)                      :: cutoffdist
    real(wp)                      :: pairlistdist
    real(wp)                      :: dielec_const
    real(wp)                      :: epsilon_rf
    real(wp)                      :: debye
    real(wp)                      :: cg_cutoffdist_vdw
    real(wp)                      :: cg_switchdist_vdw
    real(wp)                      :: cg_cutoffdist_ele
    real(wp)                      :: cg_cutoffdist_126
    real(wp)                      :: cg_cutoffdist_DNAbp
    real(wp)                      :: cg_pairlistdist_vdw
    real(wp)                      :: cg_pairlistdist_ele
    real(wp)                      :: cg_pairlistdist_126
    real(wp)                      :: cg_pairlistdist_WF
    real(wp)                      :: cg_pairlistdist_PWMcos
    real(wp)                      :: cg_pairlistdist_DNAbp
    real(wp)                      :: cg_pairlistdist_exv
    real(wp)                      :: cg_ele_coef
    real(wp)                      :: cg_ele_sol_T
    real(wp)                      :: cg_ele_sol_IC
    real(wp)                      :: cg_dielec_const
    real(wp)                      :: cg_debye_length
    real(wp)                      :: buffer_min

    logical                       :: pme_use
    real(wp)                      :: pme_alpha
    integer                       :: pme_ngrid_x
    integer                       :: pme_ngrid_y
    integer                       :: pme_ngrid_z
    integer                       :: pme_nspline
    integer                       :: fft_scheme
    real(wp)                      :: pme_max_spacing

    ! flag for position restraint 
    logical                       :: restraint_posi
    logical                       :: restraint_rmsd
    logical                       :: restraint_rmsd_target
    logical                       :: restraint_emfit
    logical                       :: restraint_pc
    ! restraint group (size = num_restraintgroups)
    integer                       :: num_rest_domain = 0
    integer                       :: rest_n_stay
    integer                       :: rest_comm_domain
    integer,          allocatable :: restraint_numatoms(:)
    integer,          allocatable :: restraint_atomlist(:,:)
    integer,          allocatable :: restraint_bondslist(:,:)
    real(wp),         allocatable :: restraint_masscoef(:,:)

    ! restraint func (size = num_restraintfuncs)
    integer,          allocatable :: restraint_kind(:)
    integer,          allocatable :: restraint_grouplist(:,:)
    integer,          allocatable :: restraint_funcgrp(:)
    integer,          allocatable :: restraint_exponent_func(:)
    integer,          allocatable :: restraint_exponent_dist(:,:)
    integer,          allocatable :: restraint_mode(:)
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

    ! for repul & fb
    real(wp),         allocatable :: restraint_rcom1(:,:)
    real(wp),         allocatable :: restraint_rcom2(:,:)
    real(wp),         allocatable :: restraint_rdrt(:)

    ! restraint func (size = num_atoms_ref)
    real(wp),         allocatable :: restraint_refcoord(:,:)

    ! restraint func (size = num_restraintfunct x ndata)
    real(wp),         allocatable :: restraint_const_replica(:,:)
    real(wp),         allocatable :: restraint_ref_replica(:,:)

    ! for fitting
    integer                       :: num_fit_domain
    integer                       :: fit_n_stay
    integer                       :: fit_comm_domain
    real(wp),         allocatable :: fit_refcoord(:,:)
    real(wp),         allocatable :: fit_coord(:,:)

    ! principal component mode
    integer                       :: num_pc_modes
    real(wp),         allocatable :: pc_mode(:)
    real(wp),         allocatable :: pc_mode_fit(:)
    integer,          allocatable :: restraint_g2pc(:)
     
    ! restraint
    logical                       :: restraint
    logical                       :: local_restraint
    logical                       :: rmsd_withmass
    logical                       :: do_emfit
    integer                       :: num_atoms_bonds_restraint
    integer                       :: num_atoms_pc_restraint
    integer                       :: nrmsd
    integer                       :: num_restraint_domain
    integer,          allocatable :: restraint_atom(:)
    integer,          allocatable :: restraint_bondslist_to_atomlist(:)
    integer,          allocatable :: restraint_pclist_to_atomlist(:)
    real(wp),         allocatable :: restraint_force(:,:)
    real(wp),         allocatable :: restraint_coord(:,:)
    real(wp),         allocatable :: restraint_bonds_coord(:,:)
    real(wp),         allocatable :: restraint_bonds_force(:,:)
    real(wp)                      :: rmsd_force
    real(wp)                      :: rmsd_target
    real(wp)                      :: pc_force(100)
    real(wp)                      :: pc_target(100)
    real(wp),         allocatable :: rotated_coord(:,:)

    ! fit
    integer                       :: fitting_method
    integer                       :: fitting_move
    integer                       :: fitting_file
    logical                       :: mass_weight
    logical                       :: do_fitting
    integer,          allocatable :: fitting_atom(:)

    ! update (bond)
    integer                       :: bond_n_stay1
    integer                       :: bond_n_stay2
    integer                       :: bonds_comm_domain
    integer                       :: bondq_comm_domain

    integer                       :: angl_n_stay1
    integer                       :: angl_n_stay2
    integer                       :: angl_n_stay3
    integer                       :: anglf_comm_domain
    integer                       :: angll_comm_domain
    integer                       :: angl_comm_domain

    integer                       :: dihe_n_stay1
    integer                       :: dihe_n_stay2
    integer                       :: dihe_n_stay3
    integer                       :: dihef_comm_domain
    integer                       :: dihel_comm_domain
    integer                       :: dihe_comm_domain

    integer                       :: stack_n_stay
    integer                       :: stack_comm_domain

    integer                       :: pwmcos_n_stay
    integer                       :: pwmcos_comm_domain

    integer                       :: pwmcosns_n_stay
    integer                       :: pwmcosns_comm_domain

    integer,          allocatable :: contact_add(:)
    integer,          allocatable :: buf_contact_integer(:,:)
    real(wp),         allocatable :: buf_contact_real(:,:)

    logical                       :: force_switch
    logical                       :: vdw_shift

    real(wp)                      :: fudge_lj
    real(wp)                      :: fudge_qq
    integer                       :: excl_level

    integer                       :: dispersion_corr
    real(wp)                      :: eswitch
    real(wp)                      :: vswitch
    real(wp)                      :: dispersion_energy
    real(wp)                      :: dispersion_virial

    ! statistical variables
    integer                       :: stats_count
    integer                       :: stats_natom
    integer                       :: stats_dimension
    real(dp),         allocatable :: stats_delta(:)
    real(dp),         allocatable :: stats_grad(:, :, :)
    real(dp),         allocatable :: stats_force(:)
    real(dp),         allocatable :: stats_metric(:,:)
    integer,          allocatable :: stats_atom(:,:)
    real(dp),         allocatable :: stats_mass(:,:)
    integer,          allocatable :: stats_id_atm2cv(:)
    integer                       :: stats_icnt
    real(dp),         allocatable :: stats_force_save(:)

    logical                       :: contact_check
    logical                       :: nonb_limiter
    logical                       :: pairlist_check
    logical                       :: bonding_check
    real(wp)                      :: minimum_contact
    real(wp)                      :: err_minimum_contact
    logical                       :: pressure_position
    logical                       :: pressure_rmsd

  end type s_enefunc

  ! parameter for allocatable variables
  integer,      public, parameter :: EneFuncBase            = 1
  integer,      public, parameter :: EneFuncBond            = 2
  integer,      public, parameter :: EneFuncAngl            = 3
  integer,      public, parameter :: EneFuncAngFlexTbl      = 4
  integer,      public, parameter :: EneFuncDihe            = 5
  integer,      public, parameter :: EneFuncDiheFlexTbl     = 6
  integer,      public, parameter :: EneFuncImpr            = 7
  integer,      public, parameter :: EneFuncBaseStack       = 8
  integer,      public, parameter :: EneFuncNbon            = 9
  integer,      public, parameter :: EneFuncNonb            = 10
  integer,      public, parameter :: EneFuncNonbList        = 11
  integer,      public, parameter :: EneFuncRefg            = 12
  integer,      public, parameter :: EneFuncReff            = 13
  integer,      public, parameter :: EneFuncRefc            = 14
  integer,      public, parameter :: EneFuncRefr            = 15
  integer,      public, parameter :: EneFuncRest            = 16
  integer,      public, parameter :: EneFuncMode            = 17
  integer,      public, parameter :: EneFuncFitc            = 18
  integer,      public, parameter :: EneFuncFitd            = 19
  integer,      public, parameter :: EneFuncRestDomain      = 20
  integer,      public, parameter :: EneFuncContact         = 21
  integer,      public, parameter :: EneFuncPWMcos          = 22
  integer,      public, parameter :: EneFuncExcl            = 23
  integer,      public, parameter :: EneFuncCGele           = 24
  integer,      public, parameter :: EneFuncDNA             = 25
  integer,      public, parameter :: EneFuncBasePair        = 26
  integer,      public, parameter :: EneFuncCGDNAExv        = 27
  integer,      public, parameter :: EneFuncCGKHmol         = 28
  integer,      public, parameter :: EneFuncCGKHList        = 29
  integer,      public, parameter :: EneFuncCGKHInvList     = 30
  integer,      public, parameter :: EneFuncCGIDRKHList     = 31
  integer,      public, parameter :: EneFuncCGIDRKHInvList  = 32
  integer,      public, parameter :: EneFuncCGIDRHPSList    = 33
  integer,      public, parameter :: EneFuncCGIDRHPSInvList = 34
  integer,      public, parameter :: EneFuncPWMcosns        = 35
  integer,      public, parameter :: EneFuncCGBaseList      = 36
  integer,      public, parameter :: EneFuncCGBaseInvList   = 37
  integer,      public, parameter :: EneFuncCGDNAList       = 38
  integer,      public, parameter :: EneFuncCGDNAInvList    = 39
  integer,      public, parameter :: EneFuncCGElecList      = 40
  integer,      public, parameter :: EneFuncCGElecInvList   = 41
  integer,      public, parameter :: EneFuncCGDNAmol        = 42
  integer,      public, parameter :: EneFuncCGIDRMPIPIList    = 43
  integer,      public, parameter :: EneFuncCGIDRMPIPIInvList = 44

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
  integer,      public, parameter :: NABaseTypeKH         = 22
  integer,      public, parameter :: NABaseTypeIDRKH      = 23
  integer,      public, parameter :: NABaseTypeIDRHPS     = 24
  integer,      public, parameter :: NABaseTypeBothKH     = 25
  integer,      public, parameter :: NABaseTypeIDRMPIPI   = 26

  ! parameters (forcefield)
  integer,      public, parameter :: ForcefieldCHARMM     = 1
  integer,      public, parameter :: ForcefieldAMBER      = 2
  integer,      public, parameter :: ForcefieldGROAMBER   = 3
  integer,      public, parameter :: ForcefieldGROMARTINI = 4
  integer,      public, parameter :: ForcefieldAAGO       = 5
  integer,      public, parameter :: ForcefieldCAGO       = 6
  integer,      public, parameter :: ForcefieldKBGO       = 7
  integer,      public, parameter :: ForcefieldSOFT       = 8
  integer,      public, parameter :: ForcefieldRESIDCG    = 9
  character(*), public, parameter :: ForceFieldTypes(9)   = (/'CHARMM    ', &
                                                              'AMBER     ', &
                                                              'GROAMBER  ', &
                                                              'GROMARTINI', &
                                                              'AAGO      ', &
                                                              'CAGO      ', &
                                                              'KBGO      ', &
                                                              'SOFT      ', &
                                                              'RESIDCG   '/)
  ! FFT scheme
  integer,      public, parameter :: FFT_1dallgather      = 1
  integer,      public, parameter :: FFT_1dalltoall       = 2
  integer,      public, parameter :: FFT_2dalltoall       = 3
  character(*), public, parameter :: FFT_Types(3)         = (/'1DALLGATHER',&
                                                              '1DALLTOALL ',&
                                                              '2DALLTOALL '/)

  ! parameters (output style)
  integer,      public, parameter :: OutputStyleGENESIS   = 1
  integer,      public, parameter :: OutputStyleCHARMM    = 2
  integer,      public, parameter :: OutputStyleNAMD      = 3
  integer,      public, parameter :: OutputStyleGROMACS   = 4
  character(*), public, parameter :: OutputStyleTypes(4)  = (/'GENESIS ', &
                                                              'CHARMM  ', &
                                                              'NAMD    ', &
                                                              'GROMACS '/)

  ! parameters (Dispersion Correction)
  integer,      public, parameter :: Disp_corr_NONE       = 1
  integer,      public, parameter :: Disp_corr_Energy     = 2
  integer,      public, parameter :: Disp_corr_EPress     = 3
  character(*), public, parameter :: Disp_corr_Types(3)   = (/'NONE  ', &
                                                              'ENERGY', &
                                                              'EPRESS'/)

  ! variables for maximum numbers in one cell (these number will be updated)

  integer,      public            :: MaxBond       = 0
  integer,      public            :: MaxAngl       = 0
  integer,      public            :: MaxDihe       = 0
  integer,      public            :: MaxStack      = 0
  integer,      public            :: MaxRest       = 0
  integer,      public            :: MaxPosi       = 0
  integer,      public            :: MaxPwmCos     = 0
  integer,      public            :: MaxPwmCosns   = 0
  integer,      public            :: MaxExcl       = 0 
  integer,      public            :: MaxExcl_all   = 18
  integer,      public            :: MaxContact    = 0
  integer,      public            :: MaxFit        = 0

  integer,      public            :: BondMoveS     = 0
  integer,      public            :: BondMoveQ     = 0
  integer,      public            :: AnglMoveF     = 0
  integer,      public            :: AnglMoveL     = 0
  integer,      public            :: AnglMove      = 0
  integer,      public            :: DiheMoveF     = 0
  integer,      public            :: DiheMoveL     = 0
  integer,      public            :: DiheMove      = 0
  integer,      public            :: StackMove     = 0
  integer,      public            :: ContactMove   = 0
  integer,      public            :: PwmCosMove    = 0
  integer,      public            :: PwmCosnsMove  = 0
  integer,      public            :: ExclMove      = 0
  integer,      public            :: RestMove      = 50

  integer,      public, parameter :: MaxAtomCls = 1000
  integer,      public            :: max_class
  real(wp),     public            :: lj_coef(2,MaxAtomCls)

  ! subroutines
  public  :: init_enefunc
  public  :: alloc_enefunc
  public  :: dealloc_enefunc
  public  :: dealloc_enefunc_all

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_enefunc
  !> @brief        Initializes energy functions information.
  !! @authors      YS, CK
  !! @param[out]   enefunc  : structure of potential energy function
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_enefunc(enefunc)

    ! formal arguments
    type(s_enefunc),         intent(inout) :: enefunc

    
    enefunc%forcefield             = ForcefieldCHARMM
    enefunc%output_style           = OutputStyleCHARMM

    enefunc%num_atom_cls           = 0
    enefunc%num_atoms_ref          = 0
    enefunc%num_restraintgroups    = 0
    enefunc%num_restraintfuncs     = 0
    enefunc%max_restraint_numatoms = 0

    enefunc%switchdist             = 0.0_wp
    enefunc%cutoffdist             = 0.0_wp
    enefunc%pairlistdist           = 0.0_wp
    enefunc%dielec_const           = 1.0_wp
    enefunc%epsilon_rf             = 0.0_wp

    enefunc%pme_use                = .false.
    enefunc%pme_alpha              = 0.0_wp
    enefunc%pme_ngrid_x            = 0
    enefunc%pme_ngrid_y            = 0
    enefunc%pme_ngrid_z            = 0
    enefunc%pme_nspline            = 0
    enefunc%pme_max_spacing        = 1.2_wp

    enefunc%restraint_posi         = .false.
    enefunc%restraint_rmsd         = .false.
    enefunc%restraint_rmsd_target  = .false.
    enefunc%restraint_emfit        = .false.
    enefunc%restraint_pc           = .false.
    enefunc%local_restraint        = .false.
    
    enefunc%force_switch            = .false.
    enefunc%vdw_shift               = .false.

    enefunc%fudge_lj                = 1.0_wp
    enefunc%fudge_qq                = 1.0_wp
    enefunc%excl_level              = 3

    enefunc%notation_14types        = 0

    enefunc%stats_count             = 0
    enefunc%stats_natom             = 0
    enefunc%stats_dimension         = 0
    enefunc%stats_icnt              = 0

    enefunc%contact_check           = .false.
    enefunc%bonding_check           = .false.
    enefunc%pairlist_check          = .false.

    enefunc%fitting_file            = 0
    enefunc%fitting_move            = 0

    enefunc%pressure_rmsd           = .false.
    enefunc%pressure_position       = .false.

    enefunc%cg_DNA_base_pair_calc   = .false.
    enefunc%cg_DNA_exv_calc         = .false.
    enefunc%cg_pwmcos_calc          = .false.
    enefunc%cg_safe_dihedral_calc   = .false.
    enefunc%cg_pwmcosns_calc        = .false.
    enefunc%cg_ele_calc             = .false.
    enefunc%cg_KH_calc              = .false.
    enefunc%cg_IDR_KH_calc          = .false.
    enefunc%cg_IDR_HPS_calc         = .false.
    enefunc%cg_IDR_MPIPI_calc       = .false.

    enefunc%cg_safe_dih_ene_shift   = 0.0_wp

    enefunc%assign_force_max        = .false.
    enefunc%upper_force_value       = 100.0_wp

    return

  end subroutine init_enefunc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_enefunc
  !> @brief        Allocates energy functions information.
  !! @authors      YS, CK
  !! @param[out]   enefunc   : potential energy functions information
  !! @param[in]    variable  : selected variable
  !! @param[in]    var_size  : size of the selected variable
  !! @param[in]    var_size1 : 2nd size of the selected variable (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_enefunc(enefunc, variable, var_size, var_size1)

    ! formal arguments
    type(s_enefunc),         intent(inout) :: enefunc
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size
    integer,       optional, intent(in)    :: var_size1

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat
    integer                  :: var_size3, var_size4


    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case (EneFuncBond)

      if (allocated(enefunc%bond_list)) then
        if (size(enefunc%bond_list(1,:)) /= var_size) &
          deallocate(enefunc%bond_list,        &
                     enefunc%bond_force_const, &
                     enefunc%bond_dist_min,    &
                     enefunc%bond_kind,        &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%bond_list)) &
        allocate(enefunc%bond_list       (3, var_size), &
                 enefunc%bond_force_const(   var_size), &
                 enefunc%bond_dist_min   (   var_size), &
                 enefunc%bond_kind       (   var_size), &
                 stat = alloc_stat)

      enefunc%bond_list       (1:3, 1:var_size) = 0
      enefunc%bond_force_const(     1:var_size) = 0.0_wp
      enefunc%bond_dist_min   (     1:var_size) = 0.0_wp
      enefunc%bond_kind       (     1:var_size) = 0

    case (EneFuncContact)

      if (allocated(enefunc%contact_list)) then
        if (size(enefunc%contact_list(1,:)) /= var_size) &
          deallocate(enefunc%contact_list,  &
                     enefunc%contact_func,  &
                     enefunc%contact_lj12,  &
                     enefunc%contact_lj10,  &
                     enefunc%contact_lj6,   &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%contact_list)) &
        allocate(enefunc%contact_list(2, var_size), &
                 enefunc%contact_func(   var_size), &
                 enefunc%contact_lj12(   var_size), &
                 enefunc%contact_lj10(   var_size), &
                 enefunc%contact_lj6 (   var_size), &
                 stat = alloc_stat)

      enefunc%contact_list(1:2, 1:var_size) = 0
      enefunc%contact_func(     1:var_size) = 0
      enefunc%contact_lj12(     1:var_size) = 0.0_wp
      enefunc%contact_lj10(     1:var_size) = 0.0_wp
      enefunc%contact_lj6 (     1:var_size) = 0.0_wp

    case (EneFuncPWMcos)

      if (allocated(enefunc%pwmcos_count)) then
        if (size(enefunc%pwmcos_count(:)) /= var_size) &
          deallocate(enefunc%pwmcos_count,             &
                     enefunc%pwmcos_protein_id,        &
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
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%pwmcos_count)) &
        allocate(enefunc%pwmcos_count         (var_size), &
                 enefunc%pwmcos_protein_id    (var_size), &
                 enefunc%pwmcos_protein_id_N  (var_size), &
                 enefunc%pwmcos_protein_id_C  (var_size), &
                 enefunc%pwmcos_r0         (6, var_size), &
                 enefunc%pwmcos_theta1     (6, var_size), &
                 enefunc%pwmcos_theta2     (6, var_size), &
                 enefunc%pwmcos_theta3     (6, var_size), &
                 enefunc%pwmcos_ene_A      (6, var_size), &
                 enefunc%pwmcos_ene_C      (6, var_size), &
                 enefunc%pwmcos_ene_G      (6, var_size), &
                 enefunc%pwmcos_ene_T      (6, var_size), &
                 enefunc%pwmcos_gamma      (6, var_size), &
                 enefunc%pwmcos_eps        (6, var_size), &
                 enefunc%pwmcos_specificity(6, var_size), &
                 stat = alloc_stat)

      enefunc%pwmcos_count           (1:var_size) = 0
      enefunc%pwmcos_protein_id      (1:var_size) = 0
      enefunc%pwmcos_protein_id_N    (1:var_size) = 0
      enefunc%pwmcos_protein_id_C    (1:var_size) = 0
      enefunc%pwmcos_r0         (1:6, 1:var_size) = 0.0_wp
      enefunc%pwmcos_theta1     (1:6, 1:var_size) = 0.0_wp
      enefunc%pwmcos_theta2     (1:6, 1:var_size) = 0.0_wp
      enefunc%pwmcos_theta3     (1:6, 1:var_size) = 0.0_wp
      enefunc%pwmcos_ene_A      (1:6, 1:var_size) = 0.0_wp
      enefunc%pwmcos_ene_C      (1:6, 1:var_size) = 0.0_wp
      enefunc%pwmcos_ene_G      (1:6, 1:var_size) = 0.0_wp
      enefunc%pwmcos_ene_T      (1:6, 1:var_size) = 0.0_wp
      enefunc%pwmcos_gamma      (1:6, 1:var_size) = 0.0_wp
      enefunc%pwmcos_eps        (1:6, 1:var_size) = 0.0_wp
      enefunc%pwmcos_specificity(1:6, 1:var_size) = 1

    case (EneFuncPWMcosns)

      if (allocated(enefunc%pwmcosns_count)) then
        if (size(enefunc%pwmcosns_count(:)) /= var_size) &
          deallocate(enefunc%pwmcosns_count,           &
                     enefunc%pwmcosns_protein_id,      &
                     enefunc%pwmcosns_protein_id_N,    &
                     enefunc%pwmcosns_protein_id_C,    &
                     enefunc%pwmcosns_r0,              &
                     enefunc%pwmcosns_theta1,          &
                     enefunc%pwmcosns_theta2,          &
                     enefunc%pwmcosns_ene,             &
                     enefunc%pwmcosns_specificity,     &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%pwmcosns_count)) &
        allocate(enefunc%pwmcosns_count         (var_size), &
                 enefunc%pwmcosns_protein_id    (var_size), &
                 enefunc%pwmcosns_protein_id_N  (var_size), &
                 enefunc%pwmcosns_protein_id_C  (var_size), &
                 enefunc%pwmcosns_r0         (6, var_size), &
                 enefunc%pwmcosns_theta1     (6, var_size), &
                 enefunc%pwmcosns_theta2     (6, var_size), &
                 enefunc%pwmcosns_ene        (6, var_size), &
                 enefunc%pwmcosns_specificity(6, var_size), &
                 stat = alloc_stat)

      enefunc%pwmcosns_count           (1:var_size) = 0
      enefunc%pwmcosns_protein_id      (1:var_size) = 0
      enefunc%pwmcosns_protein_id_N    (1:var_size) = 0
      enefunc%pwmcosns_protein_id_C    (1:var_size) = 0
      enefunc%pwmcosns_r0         (1:6, 1:var_size) = 0.0_wp
      enefunc%pwmcosns_theta1     (1:6, 1:var_size) = 0.0_wp
      enefunc%pwmcosns_theta2     (1:6, 1:var_size) = 0.0_wp
      enefunc%pwmcosns_ene        (1:6, 1:var_size) = 0.0_wp
      enefunc%pwmcosns_specificity(1:6, 1:var_size) = 1

    case (EneFuncExcl)

      if (allocated(enefunc%excl_list)) then
        if (size(enefunc%excl_list(1,:)) /= var_size) &
          deallocate(enefunc%excl_list,  &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%excl_list))    &
        allocate(enefunc%excl_list(2, var_size), &
                 stat = alloc_stat)

      enefunc%excl_list(1:2, 1:var_size) = 0

    case (EneFuncAngl)

      if (allocated(enefunc%angl_list)) then
        if (size(enefunc%angl_list(1,:)) /= var_size) &
          deallocate(enefunc%angl_list,         &
                     enefunc%angl_force_const,  &
                     enefunc%angl_theta_min,    &
                     enefunc%angl_width,        &
                     enefunc%angl_kind,         & 
                     enefunc%urey_force_const,  &
                     enefunc%urey_rmin,         &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%angl_list)) &
        allocate(enefunc%angl_list        (4, var_size),  &
                 enefunc%angl_force_const (   var_size),  &
                 enefunc%angl_theta_min   (   var_size),  &
                 enefunc%angl_width       (   var_size),  &
                 enefunc%angl_kind        (   var_size),  &
                 enefunc%urey_force_const (   var_size),  &
                 enefunc%urey_rmin        (   var_size),  &
                 stat = alloc_stat)

      enefunc%angl_list        (1:4, 1:var_size) = 0
      enefunc%angl_force_const (     1:var_size) = 0.0_wp
      enefunc%angl_theta_min   (     1:var_size) = 0.0_wp
      enefunc%angl_width       (     1:var_size) = 0.0_wp
      enefunc%angl_kind        (     1:var_size) = 0
      enefunc%urey_force_const (     1:var_size) = 0.0_wp
      enefunc%urey_rmin        (     1:var_size) = 0.0_wp

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

      allocate(enefunc%anglflex_theta(var_size1, var_size),      &
               enefunc%anglflex_efunc(var_size1, var_size),      &
               enefunc%anglflex_d2func(var_size1, var_size),     &
               enefunc%anglflex_min_th(1:2, var_size),           &
               enefunc%anglflex_max_th(1:2, var_size),           &
               enefunc%anglflex_ener_corr(var_size),             &
               stat = alloc_stat)

      enefunc%anglflex_theta (1:var_size1,1:var_size) = 0.0_wp
      enefunc%anglflex_efunc (1:var_size1,1:var_size) = 0.0_wp
      enefunc%anglflex_d2func(1:var_size1,1:var_size) = 0.0_wp
      enefunc%anglflex_min_th(1:2,1:var_size)         = 0.0_wp
      enefunc%anglflex_max_th(1:2,1:var_size)         = 0.0_wp
      enefunc%anglflex_ener_corr(1:var_size)          = 0.0_wp

    case (EneFuncDihe)

      if (allocated(enefunc%dihe_list)) then
        if (size(enefunc%dihe_list(1,:)) /= var_size) &
          deallocate(enefunc%dihe_list,        &
                     enefunc%dihe_force_const, &
                     enefunc%dihe_periodicity, &
                     enefunc%dihe_phase,       &
                     enefunc%dihe_kind,        &
                     enefunc%dihe_width,       &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%dihe_list)) &
        allocate(enefunc%dihe_list       (5, var_size), &
                 enefunc%dihe_force_const(   var_size), &
                 enefunc%dihe_periodicity(   var_size), &
                 enefunc%dihe_phase      (   var_size), &
                 enefunc%dihe_kind       (   var_size), &
                 enefunc%dihe_width      (   var_size), &
                 stat = alloc_stat)

      enefunc%dihe_list       (1:5, 1:var_size) = 0
      enefunc%dihe_force_const(     1:var_size) = 0.0_wp
      enefunc%dihe_periodicity(     1:var_size) = 0
      enefunc%dihe_phase      (     1:var_size) = 0.0_wp
      enefunc%dihe_kind       (     1:var_size) = 0
      enefunc%dihe_width      (     1:var_size) = 0.0_wp

    case(EneFuncDiheFlexTbl)
      if (allocated(enefunc%diheflex_coef)) then
        if (size(enefunc%diheflex_coef(1,:)) == var_size) return
        deallocate(enefunc%diheflex_coef,            &
                   enefunc%diheflex_ener_corr,       &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%diheflex_coef(var_size1, var_size),      &
               enefunc%diheflex_ener_corr(var_size),            &
               stat = alloc_stat)

      enefunc%diheflex_coef  (1:var_size1,1:var_size) = 0.0_wp
      enefunc%diheflex_ener_corr(1:var_size)          = 0.0_wp

    case (EneFuncBaseStack)

      if (allocated(enefunc%base_stack_list)) then
        if (size(enefunc%base_stack_list(1,:)) /= var_size) &
          deallocate(enefunc%base_stack_list,     &
                     enefunc%base_stack_func,     &
                     enefunc%base_stack_epsilon,  &
                     enefunc%base_stack_sigma,    &
                     enefunc%base_stack_theta_bs, &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%base_stack_list)) &
        allocate(enefunc%base_stack_list    (3, var_size), &
                 enefunc%base_stack_func    (   var_size), &
                 enefunc%base_stack_epsilon (   var_size), &
                 enefunc%base_stack_sigma   (   var_size), &
                 enefunc%base_stack_theta_bs(   var_size), &
                 stat = alloc_stat)

      enefunc%base_stack_list    (1:3, 1:var_size) = 0
      enefunc%base_stack_func    (     1:var_size) = 0     
      enefunc%base_stack_epsilon (     1:var_size) = 0.0_wp
      enefunc%base_stack_sigma   (     1:var_size) = 0.0_wp
      enefunc%base_stack_theta_bs(     1:var_size) = 0.0_wp

    case(EneFuncCGele)

      if (allocated(enefunc%cg_ele_mol_pair)) then
        if (size(enefunc%cg_ele_mol_pair(:,1)) /= var_size) &
          deallocate(enefunc%cg_ele_mol_pair,               &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%cg_ele_mol_pair))           &
        allocate(enefunc%cg_ele_mol_pair(var_size, var_size), &
                 stat = alloc_stat)

      enefunc%cg_ele_mol_pair(:,:) = 0

    case(EneFuncCGBaseList)

      if (allocated(enefunc%cg_base_list)) then
        if (size(enefunc%cg_base_list(:)) /= var_size)        &
          deallocate(enefunc%cg_base_list,                    &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%cg_base_list))              &
        allocate(enefunc%cg_base_list(var_size),              &
                 stat = alloc_stat)

    case(EneFuncCGBaseInvList)

      if (allocated(enefunc%cg_base_list_inv)) then
        if (size(enefunc%cg_base_list_inv(:)) /= var_size)    &
          deallocate(enefunc%cg_base_list_inv,                &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%cg_base_list_inv))          &
        allocate(enefunc%cg_base_list_inv(var_size),          &
                 stat = alloc_stat)
        enefunc%cg_base_list_inv(1:var_size) = 0

    case(EneFuncCGDNAmol)

      if (allocated(enefunc%DNA_end)) then
        if (size(enefunc%DNA_end(1,:)) /= var_size)           &
          deallocate(enefunc%DNA_end,                         &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%DNA_end))                   &
        allocate(enefunc%DNA_end(2,var_size),                 &
                 stat = alloc_stat)
        enefunc%DNA_end(1:2,1:var_size) = 0

    case(EneFuncCGDNAList)

      if (allocated(enefunc%cg_dna_list)) then
        if (size(enefunc%cg_dna_list(:)) /= var_size)         &
          deallocate(enefunc%cg_dna_list,                     &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%cg_dna_list))               &
        allocate(enefunc%cg_dna_list(var_size),               &
                 stat = alloc_stat)

    case(EneFuncCGDNAInvList)

      if (allocated(enefunc%cg_dna_list_inv)) then
        if (size(enefunc%cg_dna_list_inv(:)) /= var_size)     &
          deallocate(enefunc%cg_dna_list_inv,                 &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%cg_dna_list_inv))           &
        allocate(enefunc%cg_dna_list_inv(var_size),           &
                 stat = alloc_stat)
        enefunc%cg_dna_list_inv(1:var_size) = 0

    case(EneFuncCGKHList)

      if (allocated(enefunc%cg_KH_list)) then
        if (size(enefunc%cg_KH_list(:)) /= var_size)          &
          deallocate(enefunc%cg_KH_list,                      &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%cg_KH_list))                &
        allocate(enefunc%cg_KH_list(var_size),                &
                 stat = alloc_stat)

    case(EneFuncCGKHInvList)

      if (allocated(enefunc%cg_KH_list_inv)) then
        if (size(enefunc%cg_KH_list_inv(:)) /= var_size)      &
          deallocate(enefunc%cg_KH_list_inv,                  &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%cg_KH_list_inv))            &
        allocate(enefunc%cg_KH_list_inv(var_size),            &
                 stat = alloc_stat)
        enefunc%cg_KH_list_inv(1:var_size) = 0

    case(EneFuncCGIDRKHList)

      if (allocated(enefunc%cg_IDR_KH_list)) then
        if (size(enefunc%cg_IDR_KH_list(:)) /= var_size)      &
          deallocate(enefunc%cg_IDR_KH_list,                  &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%cg_IDR_KH_list))            &
        allocate(enefunc%cg_IDR_KH_list(var_size),            &
                 stat = alloc_stat)

    case(EneFuncCGIDRKHInvList)

      if (allocated(enefunc%cg_IDR_KH_list_inv)) then
        if (size(enefunc%cg_IDR_KH_list_inv(:)) /= var_size)  &
          deallocate(enefunc%cg_IDR_KH_list_inv,              &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%cg_IDR_KH_list_inv))        &
        allocate(enefunc%cg_IDR_KH_list_inv(var_size),        &
                 stat = alloc_stat)
        enefunc%cg_IDR_KH_list_inv(1:var_size) = 0

    case(EneFuncCGIDRHPSList)

      if (allocated(enefunc%cg_IDR_HPS_list)) then
        if (size(enefunc%cg_IDR_HPS_list(:)) /= var_size)     &
          deallocate(enefunc%cg_IDR_HPS_list,                 &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%cg_IDR_HPS_list))           &
        allocate(enefunc%cg_IDR_HPS_list(var_size),           &
                 stat = alloc_stat)

    case(EneFuncCGIDRHPSInvList)

      if (allocated(enefunc%cg_IDR_HPS_list_inv)) then
        if (size(enefunc%cg_IDR_HPS_list_inv(:)) /= var_size) &
          deallocate(enefunc%cg_IDR_HPS_list_inv,             &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%cg_IDR_HPS_list_inv))       &
        allocate(enefunc%cg_IDR_HPS_list_inv(var_size),       &
                 stat = alloc_stat)
        enefunc%cg_IDR_HPS_list_inv(1:var_size) = 0

    case(EneFuncCGIDRMPIPIList)

      if (allocated(enefunc%cg_IDR_MPIPI_list)) then
        if (size(enefunc%cg_IDR_MPIPI_list(:)) /= var_size)     &
            deallocate(enefunc%cg_IDR_MPIPI_list,                 &
            stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%cg_IDR_MPIPI_list))           &
          allocate(enefunc%cg_IDR_MPIPI_list(var_size),           &
          stat = alloc_stat)

    case(EneFuncCGIDRMPIPIInvList)

      if (allocated(enefunc%cg_IDR_MPIPI_list_inv)) then
        if (size(enefunc%cg_IDR_MPIPI_list_inv(:)) /= var_size) &
            deallocate(enefunc%cg_IDR_MPIPI_list_inv,             &
            stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%cg_IDR_MPIPI_list_inv))       &
          allocate(enefunc%cg_IDR_MPIPI_list_inv(var_size),       &
          stat = alloc_stat)
      enefunc%cg_IDR_MPIPI_list_inv(1:var_size) = 0

    case(EneFuncCGElecList)

      if (allocated(enefunc%cg_elec_list)) then
        if (size(enefunc%cg_elec_list(:)) /= var_size)        &
          deallocate(enefunc%cg_elec_list,                    &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%cg_elec_list))              &
        allocate(enefunc%cg_elec_list(var_size),              &
                 stat = alloc_stat)

    case(EneFuncCGElecInvList)

      if (allocated(enefunc%cg_elec_list_inv)) then
        if (size(enefunc%cg_elec_list_inv(:)) /= var_size)    &
          deallocate(enefunc%cg_elec_list_inv,                &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%cg_elec_list_inv))          &
        allocate(enefunc%cg_elec_list_inv(var_size),          &
                 stat = alloc_stat)
        enefunc%cg_elec_list_inv(1:var_size) = 0

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

    case (EneFuncNbon)

      if (allocated(enefunc%nonb_atom_cls)) then
        if (size(enefunc%nonb_atom_cls(:)) /= var_size) then
          deallocate(enefunc%nonb_atom_cls, &
                     enefunc%nb14_lj6,      &
                     enefunc%nb14_lj10,     &
                     enefunc%nb14_lj12,     &
                     enefunc%nonb_lj6,      &
                     enefunc%nonb_lj10,     &
                     enefunc%nonb_lj12,     &
                     stat = dealloc_stat)
          if (enefunc%forcefield == ForcefieldRESIDCG) then
            deallocate(enefunc%nonb_aicg_sig, &
                       enefunc%nonb_aicg_eps, &
                       enefunc%atom_cls_2_base_type, &
                       stat = dealloc_stat)
            if (enefunc%cg_KH_calc) then
              deallocate(enefunc%cg_KH_epsilon, &
                         enefunc%cg_KH_sigma,   &
                         stat = dealloc_stat)
            end if
            if (enefunc%cg_IDR_KH_calc) then
              deallocate(enefunc%cg_IDR_KH_epsilon, &
                         enefunc%cg_IDR_KH_sigma,   &
                         stat = dealloc_stat)
            end if
            if (enefunc%cg_IDR_HPS_calc) then
              deallocate(enefunc%cg_IDR_HPS_lambda, &
                         enefunc%cg_IDR_HPS_sigma,  &
                         stat = dealloc_stat)
            end if
            if (enefunc%cg_IDR_MPIPI_calc) then
              deallocate(                             &
                  enefunc%cg_WangFrenkel_epsilon,     &
                  enefunc%cg_WangFrenkel_sigma,       &
                  enefunc%cg_WangFrenkel_mu,          &
                  enefunc%cg_WangFrenkel_nu,          &
                  stat = dealloc_stat)
            end if
          end if
        end if 
      end if

      if (.not. allocated(enefunc%nonb_atom_cls)) then
        allocate(enefunc%nonb_atom_cls(var_size),           &
                 enefunc%nb14_lj6     (var_size, var_size), &
                 enefunc%nb14_lj10    (var_size, var_size), &
                 enefunc%nb14_lj12    (var_size, var_size), &
                 enefunc%nonb_lj6     (var_size, var_size), &
                 enefunc%nonb_lj10    (var_size, var_size), &
                 enefunc%nonb_lj12    (var_size, var_size), &
                 stat = alloc_stat)
        if (enefunc%forcefield == ForcefieldRESIDCG) then
          allocate(enefunc%nonb_aicg_sig(var_size, var_size), &
                   enefunc%nonb_aicg_eps(var_size, var_size), &
                   enefunc%atom_cls_2_base_type(var_size),    &
                   stat = alloc_stat)
          if (enefunc%cg_KH_calc)                               &
            allocate(enefunc%cg_KH_epsilon(var_size, var_size), &
                     enefunc%cg_KH_sigma  (var_size, var_size), &
                     stat = alloc_stat)
          if (enefunc%cg_IDR_KH_calc)                               &
            allocate(enefunc%cg_IDR_KH_epsilon(var_size, var_size), &
                     enefunc%cg_IDR_KH_sigma  (var_size, var_size), &
                     stat = alloc_stat)

          if (enefunc%cg_IDR_HPS_calc)                              &
            allocate(enefunc%cg_IDR_HPS_lambda(var_size, var_size), &
                     enefunc%cg_IDR_HPS_sigma (var_size, var_size), &
                     stat = alloc_stat)

          if (enefunc%cg_IDR_MPIPI_calc)                                 &
            allocate(enefunc%cg_WangFrenkel_epsilon(var_size, var_size), &
                     enefunc%cg_WangFrenkel_sigma(var_size, var_size),   &
                     enefunc%cg_WangFrenkel_mu(var_size, var_size),      &
                     enefunc%cg_WangFrenkel_nu(var_size, var_size),      &
                     stat = alloc_stat)

        end if

      end if

      enefunc%nonb_atom_cls(1:var_size)             = 0
      enefunc%nb14_lj6     (1:var_size, 1:var_size) = 0.0_wp
      enefunc%nb14_lj10    (1:var_size, 1:var_size) = 0.0_wp
      enefunc%nb14_lj12    (1:var_size, 1:var_size) = 0.0_wp
      enefunc%nonb_lj6     (1:var_size, 1:var_size) = 0.0_wp
      enefunc%nonb_lj10    (1:var_size, 1:var_size) = 0.0_wp
      enefunc%nonb_lj12    (1:var_size, 1:var_size) = 0.0_wp
      if (enefunc%forcefield == ForcefieldRESIDCG) then
        enefunc%nonb_aicg_sig(1:var_size, 1:var_size) = 0.0_wp
        enefunc%nonb_aicg_eps(1:var_size, 1:var_size) = 0.0_wp
        enefunc%atom_cls_2_base_type(1:var_size) = 0.0_wp
        if (enefunc%cg_KH_calc) then
          enefunc%cg_KH_epsilon(1:var_size, 1:var_size) = 0.0_wp
          enefunc%cg_KH_sigma  (1:var_size, 1:var_size) = 0.0_wp
        end if
        if (enefunc%cg_IDR_KH_calc) then
          enefunc%cg_IDR_KH_epsilon(1:var_size, 1:var_size) = 0.0_wp
          enefunc%cg_IDR_KH_sigma  (1:var_size, 1:var_size) = 0.0_wp
        end if
        if (enefunc%cg_IDR_HPS_calc) then
          enefunc%cg_IDR_HPS_lambda(1:var_size, 1:var_size) = 0.0_wp
          enefunc%cg_IDR_HPS_sigma (1:var_size, 1:var_size) = 0.0_wp
        end if
        if (enefunc%cg_IDR_MPIPI_calc) then
          enefunc%cg_WangFrenkel_sigma   (1:var_size, 1:var_size) = 0.0_wp
          enefunc%cg_WangFrenkel_epsilon (1:var_size, 1:var_size) = 0.0_wp
          enefunc%cg_WangFrenkel_mu      (1:var_size, 1:var_size) = 1
          enefunc%cg_WangFrenkel_nu      (1:var_size, 1:var_size) = 1
        end if
      end if

    case (EneFuncNonb)

      if (allocated(enefunc%exclusion_mask)) then
        if (size(enefunc%exclusion_mask) /= var_size*var_size) &
          deallocate(enefunc%exclusion_mask,  &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%exclusion_mask)) &
        allocate(enefunc%exclusion_mask(var_size, var_size),   &
                 stat = alloc_stat)

    case (EneFuncRefg)

      if (allocated(enefunc%restraint_numatoms)) then
        if (size(enefunc%restraint_numatoms) /= var_size) &
          deallocate(enefunc%restraint_numatoms,  &
                     enefunc%restraint_atomlist,  &
                     enefunc%restraint_bondslist, &
                     enefunc%restraint_masscoef,  &
                     enefunc%restraint_wcom3,     &
                     enefunc%restraint_wcom4,     &
                     enefunc%restraint_wcom5,     &
                     enefunc%restraint_wtmp,      &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%restraint_numatoms)) &
        allocate(enefunc%restraint_numatoms (var_size),                &
                 enefunc%restraint_atomlist (1:var_size1, 1:var_size), &
                 enefunc%restraint_bondslist(1:var_size1, 1:var_size), &
                 enefunc%restraint_masscoef (1:var_size1, 1:var_size), &
                 enefunc%restraint_wcom3    (1:3, 1:var_size1),        &
                 enefunc%restraint_wcom4    (1:3, 1:var_size1),        &
                 enefunc%restraint_wcom5    (1:3, 1:var_size1),        &
                 enefunc%restraint_wtmp     (1:var_size1, 1:var_size), &
                 stat = alloc_stat)

      enefunc%restraint_numatoms (1:var_size)              = 0
      enefunc%restraint_atomlist (1:var_size1, 1:var_size) = 0
      enefunc%restraint_bondslist(1:var_size1, 1:var_size) = 0
      enefunc%restraint_masscoef (1:var_size1, 1:var_size) = 0.0_wp
      enefunc%restraint_wcom3    (1:3, 1:var_size1)        = 0.0_wp
      enefunc%restraint_wcom4    (1:3, 1:var_size1)        = 0.0_wp
      enefunc%restraint_wcom5    (1:3, 1:var_size1)        = 0.0_wp
      enefunc%restraint_wtmp     (1:var_size1, 1:var_size) = 0.0_wp

    case (EneFuncReff)

      if (allocated(enefunc%restraint_kind)) then
        if (size(enefunc%restraint_kind) /= var_size) &
          deallocate(enefunc%restraint_kind,          &
                     enefunc%restraint_grouplist,     &
                     enefunc%restraint_const,         &
                     enefunc%restraint_ref,           &
                     enefunc%restraint_funcgrp,       &
                     enefunc%restraint_exponent_func, &
                     enefunc%restraint_exponent_dist, &
                     enefunc%restraint_mode,          &
                     enefunc%restraint_weight_dist,   &
                     enefunc%restraint_wcom1,         &
                     enefunc%restraint_wcom2,         &
                     enefunc%restraint_wdrt,          &
                     enefunc%restraint_rcom1,         &
                     enefunc%restraint_rcom2,         &
                     enefunc%restraint_rdrt,          &
                     stat = dealloc_stat)
      end if

      var_size3 = max(int(var_size1/2),1)
      var_size4 = int(var_size1*(var_size1-1)/2)

      if (.not. allocated(enefunc%restraint_kind)) &
        allocate(enefunc%restraint_kind         (var_size),                &
                 enefunc%restraint_grouplist    (1:var_size1, 1:var_size), &
                 enefunc%restraint_const        (1:4, 1:var_size),         &
                 enefunc%restraint_ref          (1:2, 1:var_size),         &
                 enefunc%restraint_funcgrp      (1:var_size),              &
                 enefunc%restraint_exponent_func(1:var_size),              &
                 enefunc%restraint_exponent_dist(1:var_size3, 1:var_size), &
                 enefunc%restraint_mode(1:var_size),                       &
                 enefunc%restraint_weight_dist  (1:var_size3, 1:var_size), &
                 enefunc%restraint_wcom1        (1:3, 1:var_size1),        &
                 enefunc%restraint_wcom2        (1:3, 1:var_size1),        &
                 enefunc%restraint_wdrt         (1:var_size3),             &
                 enefunc%restraint_rcom1        (1:3, 1:var_size1),        &
                 enefunc%restraint_rcom2        (1:3, 1:var_size4),        &
                 enefunc%restraint_rdrt         (1:var_size4),             &
                 stat = alloc_stat)

      enefunc%restraint_kind         (1:var_size)              = 0
      enefunc%restraint_grouplist    (1:var_size1, 1:var_size) = 0
      enefunc%restraint_const        (1:4, 1:var_size)         = 0.0_wp
      enefunc%restraint_ref          (1:2, 1:var_size)         = 0.0_wp
      enefunc%restraint_funcgrp      (1:var_size)              = 0
      enefunc%restraint_exponent_func(1:var_size)              = 0
      enefunc%restraint_exponent_dist(1:var_size3, 1:var_size) = 0
      enefunc%restraint_mode(1:var_size)                       = 0
      enefunc%restraint_weight_dist  (1:var_size3, 1:var_size) = 0.0_wp
      enefunc%restraint_wcom1        (1:3, 1:var_size1)        = 0.0_wp
      enefunc%restraint_wcom2        (1:3, 1:var_size1)        = 0.0_wp
      enefunc%restraint_wdrt         (1:var_size3)             = 0.0_wp
      enefunc%restraint_rcom1        (1:3, 1:var_size1)        = 0.0_wp
      enefunc%restraint_rcom2        (1:3, 1:var_size4)        = 0.0_wp
      enefunc%restraint_rdrt         (1:var_size4)             = 0.0_wp

    case (EneFuncRefc)

      if (allocated(enefunc%restraint_refcoord)) then
        if (size(enefunc%restraint_refcoord) /= var_size) &
          deallocate(enefunc%restraint_refcoord, stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%restraint_refcoord)) &
        allocate(enefunc%restraint_refcoord(1:3, var_size), stat = alloc_stat)

      enefunc%restraint_refcoord(1:3, 1:var_size) = 0.0_wp

    case(EneFuncRefr)

      if (allocated(enefunc%restraint_const_replica)) then
        if (size(enefunc%restraint_const_replica(1,:)) == var_size) return
        deallocate(enefunc%restraint_const_replica, &
                   enefunc%restraint_ref_replica,  &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%restraint_const_replica(var_size1,var_size), &
               enefunc%restraint_ref_replica(var_size1,var_size), &
               stat = alloc_stat)

      enefunc%restraint_const_replica(1:var_size1,1:var_size) = 0.0_wp
      enefunc%restraint_ref_replica(1:var_size1,1:var_size)   = 0.0_wp

    case (EneFuncRest)

      if (allocated(enefunc%restraint_atom)) then
        if (size(enefunc%restraint_bondslist_to_atomlist(:)) /= var_size) &
          deallocate(enefunc%restraint_bondslist_to_atomlist, &
                     enefunc%restraint_bonds_coord,           &
                     enefunc%restraint_bonds_force,           &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%restraint_bondslist_to_atomlist))  &
        allocate(enefunc%restraint_bondslist_to_atomlist(var_size),  &
                 enefunc%restraint_bonds_coord(3, var_size),         &
                 enefunc%restraint_bonds_force(3, var_size),         &
                 stat = alloc_stat)

      enefunc%restraint_bondslist_to_atomlist(1:var_size) = 0
      enefunc%restraint_bonds_coord(1:3, 1:var_size)      = 0.0_wp
      enefunc%restraint_bonds_force(1:3, 1:var_size)      = 0.0_wp

    case (EneFuncRestDomain)

      if (allocated(enefunc%restraint_atom)) then
        if (size(enefunc%restraint_atom(:)) /= var_size) &
          deallocate(enefunc%restraint_atom,                  &
                     enefunc%restraint_force,                 &
                     enefunc%restraint_coord,                 &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%restraint_atom))            &
        allocate(enefunc%restraint_atom (var_size),           &
                 enefunc%restraint_force(4, var_size),        &
                 enefunc%restraint_coord(3, var_size),        &
                 stat = alloc_stat)

      enefunc%restraint_atom (1:var_size)                  = 0
      enefunc%restraint_force(1:4, 1:var_size)             = 0.0_wp
      enefunc%restraint_coord(1:3, 1:var_size)             = 0.0_wp

    case(EneFuncMode)

      if (allocated(enefunc%pc_mode)) then
        if (size(enefunc%pc_mode) == var_size) return
        deallocate(enefunc%pc_mode,        &
                   enefunc%pc_mode_fit,    &
                   enefunc%restraint_g2pc, &
                   stat = dealloc_stat)
      end if

      allocate(enefunc%pc_mode(1:var_size),         &
               enefunc%pc_mode_fit(1:var_size),     &
               enefunc%restraint_g2pc(1:var_size1), &
               stat = alloc_stat)
      enefunc%pc_mode(1:var_size) = 0.0_dp
      enefunc%pc_mode_fit(1:var_size) = 0.0_dp
      enefunc%restraint_g2pc(1:var_size1) = 0

    case (EneFuncFitc)

      if (allocated(enefunc%fit_refcoord)) then
        if (size(enefunc%fit_refcoord(1,:)) /= var_size) &
          deallocate(enefunc%fit_refcoord, stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%fit_refcoord)) &
        allocate(enefunc%fit_refcoord(1:3, 1:var_size), stat = alloc_stat)

      enefunc%fit_refcoord(1:3, 1:var_size) = 0.0_wp

    case (EneFuncFitd)

      if (allocated(enefunc%fit_coord)) then
        if (size(enefunc%fit_coord(1,:)) /= var_size) then
          deallocate(enefunc%fitting_atom, &
                     enefunc%fit_coord,    &
                     stat = dealloc_stat)
        end if
      end if

      if (.not. allocated(enefunc%fit_coord))        &
        allocate(enefunc%fitting_atom(1:var_size),   &
                 enefunc%fit_coord(1:3, 1:var_size), &
                 stat = alloc_stat)

      enefunc%fitting_atom(1:var_size)    = 0
      enefunc%fit_coord(1:3, 1:var_size)  = 0.0_wp

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

    case default

      call error_msg('Alloc_Enefunc> Invalid variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc    

    return

  end subroutine alloc_enefunc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_enefunc
  !> @brief        Deallocates energy functions information.
  !! @authors      YS, CK
  !! @param[inout] enefunc  : potential energy functions information
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

    case (EneFuncBond)

      if (allocated(enefunc%bond_list)) then
        deallocate(enefunc%bond_list,        &
                   enefunc%bond_force_const, &
                   enefunc%bond_dist_min,    &
                   enefunc%bond_kind,        &
                   stat = dealloc_stat)
      end if

    case (EneFuncAngl)

      if (allocated(enefunc%angl_list)) then
        deallocate(enefunc%angl_list,         &
                   enefunc%angl_force_const,  &
                   enefunc%angl_theta_min,    &
                   enefunc%angl_width,        &
                   enefunc%angl_kind,         &
                   enefunc%urey_force_const,  &
                   enefunc%urey_rmin,         &
                   stat = dealloc_stat)
      end if

    case (EneFuncDihe)

      if (allocated(enefunc%dihe_list)) then
        deallocate(enefunc%dihe_list,        &
                   enefunc%dihe_force_const, &
                   enefunc%dihe_periodicity, &
                   enefunc%dihe_phase,       &
                   enefunc%dihe_kind,        &
                   enefunc%dihe_width,       &
                   stat = dealloc_stat)
      end if

    case (EneFuncContact)

      if (allocated(enefunc%contact_list)) then
        deallocate(enefunc%contact_list,  &
                   enefunc%contact_func,  &
                   enefunc%contact_lj12,  &
                   enefunc%contact_lj10,  &
                   enefunc%contact_lj6,   &
                   stat = dealloc_stat)
      end if

    case (EneFuncNbon)

      if (allocated(enefunc%nonb_atom_cls)) then
        deallocate(enefunc%nonb_atom_cls, &
                   enefunc%nb14_lj6,      &
                   enefunc%nb14_lj10,     &
                   enefunc%nb14_lj12,     &
                   enefunc%nonb_lj6,      &
                   enefunc%nonb_lj10,     &
                   enefunc%nonb_lj12,     &
                   stat = dealloc_stat)
        if (enefunc%forcefield == ForcefieldRESIDCG) then
          deallocate(enefunc%nonb_aicg_sig, &
                     enefunc%nonb_aicg_eps, &
                     enefunc%atom_cls_2_base_type, &
                     stat = dealloc_stat)
          if (enefunc%cg_KH_calc) then
            deallocate(enefunc%cg_KH_epsilon, &
                       enefunc%cg_KH_sigma,   &
                       stat = dealloc_stat)
          end if
          if (enefunc%cg_IDR_KH_calc) then
            deallocate(enefunc%cg_IDR_KH_epsilon, &
                       enefunc%cg_IDR_KH_sigma,   &
                       stat = dealloc_stat)
          end if
          if (enefunc%cg_IDR_HPS_calc) then
            deallocate(enefunc%cg_IDR_HPS_lambda, &
                       enefunc%cg_IDR_HPS_sigma,  &
                       stat = dealloc_stat)
          end if
        end if
      end if

      if (enefunc%cg_IDR_MPIPI_calc.and.allocated(enefunc%cg_WangFrenkel_epsilon)) then
        deallocate(                             &
            enefunc%cg_WangFrenkel_epsilon,     &
            enefunc%cg_WangFrenkel_sigma,       &
            enefunc%cg_WangFrenkel_mu,          &
            enefunc%cg_WangFrenkel_nu,          &
            stat = dealloc_stat)
      end if

    case (EneFuncNonb)

      if (allocated(enefunc%exclusion_mask)) then
        deallocate(enefunc%exclusion_mask  ,   &
                   stat = dealloc_stat)
      end if

    case (EneFuncRest)

      if (allocated(enefunc%restraint_bondslist_to_atomlist)) then
        deallocate(enefunc%restraint_bondslist_to_atomlist, &
                   enefunc%restraint_bonds_coord,           &
                   enefunc%restraint_bonds_force,           &
                   stat = dealloc_stat)
      end if

    case (EneFuncRestDomain)

      if (allocated(enefunc%restraint_atom)) then
        deallocate(enefunc%restraint_atom,                  &
                   enefunc%restraint_force,                 &
                   enefunc%restraint_coord,                 &
                   stat = dealloc_stat)
      end if

    case (EneFuncExcl)

      if (allocated(enefunc%excl_list)) then
        deallocate(enefunc%excl_list,  &
                   stat = dealloc_stat)
      end if

    case default

      call error_msg('Error: Dealloc_Enefunc - Invalid variable specified.')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_enefunc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_enefunc_all
  !> @brief        Deallocates all energy functions information.
  !! @authors      YS, CK
  !! @param[out]   enefunc : structure of potential energy function
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_enefunc_all(enefunc)

    ! formal arguments
    type(s_enefunc),         intent(inout) :: enefunc


    call dealloc_enefunc(enefunc, EneFuncBond)
    call dealloc_enefunc(enefunc, EneFuncAngl)
    call dealloc_enefunc(enefunc, EneFuncDihe)
    call dealloc_enefunc(enefunc, EneFuncExcl)
    call dealloc_enefunc(enefunc, EneFuncContact)
    call dealloc_enefunc(enefunc, EneFuncNbon)
    call dealloc_enefunc(enefunc, EneFuncNonb)
    call dealloc_enefunc(enefunc, EneFuncRest)

    return

  end subroutine dealloc_enefunc_all

end module cg_enefunc_str_mod
