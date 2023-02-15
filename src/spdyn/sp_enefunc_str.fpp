!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_enefunc_str_mod
!> @brief   structure of energy functions
!! @authors Yuji Sugita (YS), Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_enefunc_str_mod

  use sp_domain_str_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_table_enefunc

    ! General table
    logical                       :: table
    logical                       :: water_table
    logical                       :: water_bond_calc = .false.
    logical                       :: water_angle_calc = .false.
    logical                       :: water_bond_calc_OH = .false.
    logical                       :: water_bond_calc_HH = .false.
    logical                       :: TIP2_bond_calc     = .false.
    integer                       :: atom_cls_no_1
    integer                       :: atom_cls_no_2
    integer                       :: atom_cls_no_O
    integer                       :: atom_cls_no_H
    integer                       :: atom_cls_no_D
    integer                       :: table_order
    real(wp)                      :: density
    real(wp)                      :: charge_1
    real(wp)                      :: charge_2
    real(wp)                      :: charge_O
    real(wp)                      :: charge_H
    real(wp)                      :: charge_D
    real(wp)                      :: mass_1
    real(wp)                      :: mass_2
    real(wp)                      :: mass_O
    real(wp)                      :: mass_H
    real(wp)                      :: mass_D
    real(wp)                      :: OH_bond
    real(wp)                      :: OH_force = 0.0_wp
    real(wp)                      :: HH_bond
    real(wp)                      :: HH_force = 0.0_wp
    real(wp)                      :: HOH_angle
    real(wp)                      :: HOH_force = 0.0_wp
    real(wp)                      :: TIP2_bond
    real(wp)                      :: TIP2_force = 0.0_wp

    real(wp),         allocatable :: table_ene(:)
    real(wp),         allocatable :: table_grad(:)
    real(wp),         allocatable :: table_ecor(:)
    real(wp),         allocatable :: table_decor(:)
    real(wp),         allocatable :: table_vcor(:)
    real(wp),         allocatable :: table_dvcor(:)

    ! Water
    character(5)                  :: water_model
    integer                       :: num_water
    integer                       :: num_solute
    integer                       :: num_all
    integer,          allocatable :: water_list(:,:)
    integer,          allocatable :: solute_list(:)
    integer,          allocatable :: solute_list_inv(:)

    real(wp),         allocatable :: table_ene_WW(:,:)
    real(wp),         allocatable :: table_de_WW(:,:)

    ! for allocation on GPU
    integer                       :: cutoff_int

  end type s_table_enefunc

  ! GaMD
  type, public :: s_enefunc_gamd
    ! gamd
    logical                       :: gamd_stat
    logical                       :: gamd_boost
    logical                       :: boost_pot
    logical                       :: boost_dih
    logical                       :: boost_dual
    logical                       :: thresh_lower
    logical                       :: thresh_upper
    real(dp)                      :: ene_pot_max
    real(dp)                      :: ene_pot_min
    real(dp)                      :: ene_pot_ave
    real(dp)                      :: ene_pot_ave2
    real(dp)                      :: ene_pot_dev
    real(dp)                      :: ene_dih_max
    real(dp)                      :: ene_dih_min
    real(dp)                      :: ene_dih_ave
    real(dp)                      :: ene_dih_ave2
    real(dp)                      :: ene_dih_dev
    integer                       :: count_pot
    integer                       :: count_dih
    real(dp)                      :: ene_pot_th
    real(dp)                      :: ene_dih_th
    real(dp)                      :: k0_pot
    real(dp)                      :: k0_dih
    real(dp)                      :: k_pot
    real(dp)                      :: k_dih
    real(dp)                      :: sigma0_pot
    real(dp)                      :: sigma0_dih
    integer                       :: update_period
    real(wp), allocatable         :: f_dihe_omp(:,:,:,:)
    real(wp), allocatable         :: f_rest_omp(:,:,:,:)
    real(dp), allocatable         :: v_dihe_omp(:,:,:)
    real(dp), allocatable         :: v_rest_omp(:,:,:)
  end type s_enefunc_gamd

  type, public :: s_enefunc

    integer                       :: forcefield
    integer                       :: output_style
    integer                       :: vdw          
    logical                       :: enm_use      

    integer                       :: num_bond_all
    integer                       :: num_enm_all
    integer                       :: num_angl_all
    integer                       :: num_dihe_all
    integer                       :: num_rb_dihe_all
    integer                       :: num_impr_all
    integer                       :: num_cmap_all
    integer                       :: num_contact_all
    integer                       :: num_excl_all
    integer                       :: num_nb14_all

    integer                       :: num_atom_cls
    integer                       :: num_atoms_ref
    integer                       :: num_restraintgroups
    integer                       :: num_restraintfuncs
    integer                       :: max_restraint_numatoms
    integer                       :: max_restraint_numgrps

    integer                       :: notation_14types
    real(wp),         allocatable :: dihe_scnb(:), dihe_scee(:)

    ! base
    integer,          allocatable :: num_bond(:)
    integer,          allocatable :: num_enm(:)
    integer,          allocatable :: num_angle(:)
    integer,          allocatable :: num_dihedral(:)
    integer,          allocatable :: num_rb_dihedral(:)
    integer,          allocatable :: num_improper(:)
    integer,          allocatable :: num_cmap(:)
    integer,          allocatable :: num_restraint(:)
    integer,          allocatable :: num_contact(:)
    integer,          allocatable :: num_excl(:)
    integer,          allocatable :: num_excl14(:)

    ! base for parallel I/O input
    integer,          allocatable :: num_bond_pio(:,:)
    integer,          allocatable :: num_angle_pio(:,:)
    integer,          allocatable :: num_dihedral_pio(:,:)
    integer,          allocatable :: num_rb_dihedral_pio(:,:)
    integer,          allocatable :: num_improper_pio(:,:)
    integer,          allocatable :: num_cmap_pio(:,:)
    integer,          allocatable :: num_restraint_pio(:,:)

    ! bond (size = num_local_bond for each cell)
    integer,          allocatable :: bond_list(:,:,:)
    real(wp),         allocatable :: bond_force_const(:,:)
    real(wp),         allocatable :: bond_dist_min(:,:)
    integer(1),       allocatable :: bond_kind(:,:)
    integer,          allocatable :: bond_pbc (:,:)

    ! bond (size = num_local_bond for each cell)
    integer,          allocatable :: enm_list(:,:,:)
    real(wp),         allocatable :: enm_force_const(:,:)
    real(wp),         allocatable :: enm_dist_min(:,:)

    ! excl and excl14
    integer,          allocatable :: excl_list(:,:,:)
    integer,          allocatable :: excl14_list(:,:,:)

    ! bond for parallel I/O input
    integer,          allocatable :: bond_list_pio(:,:,:,:)
    real(wp),         allocatable :: bond_force_const_pio(:,:,:)
    real(wp),         allocatable :: bond_dist_min_pio(:,:,:)
    integer,          allocatable :: bond_pbc_pio(:,:,:)

    ! angle (size = num_local_angle for each cell)
    integer,          allocatable :: angle_list(:,:,:)
    real(wp),         allocatable :: angle_force_const(:,:)
    real(wp),         allocatable :: angle_theta_min(:,:)
    real(wp),         allocatable :: urey_force_const(:,:)
    real(wp),         allocatable :: urey_rmin(:,:)
    integer(1),       allocatable :: angle_kind(:,:)
    integer,          allocatable :: angle_pbc(:,:,:)

    ! angle for parallel I/O input
    integer,          allocatable :: angle_list_pio(:,:,:,:)
    real(wp),         allocatable :: angle_force_const_pio(:,:,:)
    real(wp),         allocatable :: angle_theta_min_pio(:,:,:)
    real(wp),         allocatable :: urey_force_const_pio(:,:,:)
    real(wp),         allocatable :: urey_rmin_pio(:,:,:)
    integer,          allocatable :: angle_pbc_pio(:,:,:,:)
  
    ! dihedral (size = num_local_dihedral for each cell)
    integer,          allocatable :: dihe_list(:,:,:)
    real(wp),         allocatable :: dihe_force_const(:,:)
    integer,          allocatable :: dihe_periodicity(:,:)
    real(wp),         allocatable :: dihe_phase(:,:)
    integer(1),       allocatable :: dihe_kind(:,:)
    integer,          allocatable :: dihe_pbc (:,:,:)

    ! dihedral for parallel I/O input
    integer,          allocatable :: dihe_list_pio(:,:,:,:)
    real(wp),         allocatable :: dihe_force_const_pio(:,:,:)
    integer,          allocatable :: dihe_periodicity_pio(:,:,:)
    real(wp),         allocatable :: dihe_phase_pio(:,:,:)
    integer,          allocatable :: dihe_pbc_pio(:,:,:,:)

    ! dihedral (size = num_local_rb_dihedrals for each cell)
    integer,          allocatable :: rb_dihe_list(:,:,:)
    real(wp),         allocatable :: rb_dihe_c(:,:,:)
    integer,          allocatable :: rb_dihe_pbc(:,:,:)

    ! dihedral for parallel I/O input
    integer,          allocatable :: rb_dihe_list_pio(:,:,:,:)
    real(wp),         allocatable :: rb_dihe_c_pio(:,:,:,:)
    integer,          allocatable :: rb_dihe_pbc_pio(:,:,:,:)

    ! improper (size = num_local_improper for each cell)
    integer,          allocatable :: impr_list(:,:,:)
    real(wp),         allocatable :: impr_force_const(:,:)
    integer,          allocatable :: impr_periodicity(:,:)
    real(wp),         allocatable :: impr_phase(:,:)
    integer,          allocatable :: impr_pbc(:,:,:)

    ! improper for parallel I/O input
    integer,          allocatable :: impr_list_pio(:,:,:,:)
    real(wp),         allocatable :: impr_force_const_pio(:,:,:)
    integer,          allocatable :: impr_periodicity_pio(:,:,:)
    real(wp),         allocatable :: impr_phase_pio(:,:,:)
    integer,          allocatable :: impr_pbc_pio(:,:,:,:)

    ! cmap
    integer                       :: cmap_ngrid0
    integer                       :: cmap_ncmap_p
    integer,          allocatable :: cmap_list(:,:,:)
    integer,          allocatable :: cmap_resolution(:)
    integer,          allocatable :: cmap_type(:,:)
    real(wp),         allocatable :: cmap_coef(:,:,:,:,:)
    integer,          allocatable :: cmap_pbc(:,:,:)

    ! cmap for parallel I/O input
    integer,          allocatable :: cmap_list_pio(:,:,:,:)
    integer,          allocatable :: cmap_resolution_pio(:)
    integer,          allocatable :: cmap_type_pio(:,:,:)
    real(wp),         allocatable :: cmap_coef_pio(:,:,:,:,:)
    integer,          allocatable :: cmap_pbc_pio(:,:,:,:)

    ! contact (size = num_local_contact for each cell)
    integer,          allocatable :: contact_list(:,:,:)
    real(wp),         allocatable :: contact_lj12(:,:)
    real(wp),         allocatable :: contact_lj6(:,:)

    ! non-bonded (size = num_atom_cls)
    integer,          allocatable :: nonb_atom_cls(:)
    integer,          allocatable :: nonb_atom_cls_no(:)
    integer,          allocatable :: nonb_ljtype(:,:)
    real(wp),         allocatable :: nb14_lj6(:,:)
    real(wp),         allocatable :: nb14_lj12(:,:)
    real(wp),         allocatable :: nonb_lj6(:,:)
    real(wp),         allocatable :: nonb_lj12(:,:)
    real(wp),         allocatable :: nonb_lj6_factor(:)

    ! 1-4 (size = num_atoms)
    integer,          allocatable :: num_nonb_excl(:)
    integer,          allocatable :: num_nb14_calc(:)
    integer(1),       allocatable :: exclusion_mask(:,:,:)
    integer(1),       allocatable :: exclusion_mask1(:,:,:)

    ! exclusion
    integer,          allocatable :: nonb_excl_list(:,:,:)
    integer,          allocatable :: nb14_calc_list(:,:,:)
    integer,          allocatable :: sc_calc_list(:,:)
    real(wp),         allocatable :: nb14_qq_scale(:,:)
    real(wp),         allocatable :: nb14_lj_scale(:,:)

    real(wp)                      :: switchdist
    real(wp)                      :: cutoffdist
    real(wp)                      :: pairlistdist
    real(wp)                      :: dielec_const

    logical                       :: pme_use
    real(wp)                      :: pme_alpha
    integer                       :: pme_ngrid_x
    integer                       :: pme_ngrid_y
    integer                       :: pme_ngrid_z
    integer                       :: pme_nspline
    integer                       :: pme_scheme
    real(wp)                      :: pme_max_spacing

    ! external electric field
    !
    logical                       :: use_efield
    real(wp)                      :: efield(3)
    logical                       :: efield_virial
    logical                       :: efield_normal

    ! flag for position restraint 
    logical                       :: restraint_posi
    logical                       :: restraint_rmsd
    logical                       :: restraint_rmsd_target
    logical                       :: restraint_emfit
    logical                       :: restraint_pc
    ! restraint group (size = num_restraintgroups)
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
    integer,          allocatable :: restraint_rpath_func(:)
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
    real(wp),         allocatable :: fit_refcoord(:,:)
    real(wp),         allocatable :: fit_coord(:,:,:)

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
    integer,          allocatable :: restraint_atom(:,:)
    integer,          allocatable :: restraint_bondslist_to_atomlist(:)
    integer,          allocatable :: restraint_pclist_to_atomlist(:)
    real(wp),         allocatable :: restraint_force(:,:,:)
    real(wp),         allocatable :: restraint_coord(:,:,:)
    real(wp),         allocatable :: restraint_bonds_coord(:,:)
    real(wp),         allocatable :: restraint_bonds_force(:,:)
    real(wp)                      :: rmsd_force
    real(wp)                      :: rmsd_target
    real(wp)                      :: pc_force(100)
    real(wp)                      :: pc_target(100)
    real(wp),         allocatable :: rotated_coord(:,:,:)

    ! restraint for parallel I/O input
    integer,          allocatable :: restraint_atom_pio(:,:,:)
    real(wp),         allocatable :: restraint_force_pio(:,:,:,:)
    real(wp),         allocatable :: restraint_coord_pio(:,:,:,:)

    ! fit
    integer                       :: fitting_method
    integer                       :: fitting_move
    integer                       :: fitting_file
    logical                       :: mass_weight
    logical                       :: do_fitting
    integer,          allocatable :: nfitting(:)
    integer,          allocatable :: fitting_atom(:,:)
    integer,          allocatable :: fitting_exit(:)
    integer,          allocatable :: fitting_add(:)
    integer,          allocatable :: fitting_exit_index(:,:)
    integer,          allocatable :: buf_fitting_integer(:,:)
    real(wp),         allocatable :: buf_fitting_real(:,:,:)

    ! update
    integer,          allocatable :: bond_exit(:)
    integer,          allocatable :: bond_add(:)
    integer,          allocatable :: bond_exit_index(:,:)
    integer,          allocatable :: buf_bond_integer(:,:,:)
    real(wp),         allocatable :: buf_bond_real(:,:,:)

    ! update
    integer,          allocatable :: enm_add(:)
    integer,          allocatable :: buf_enm_integer(:,:,:)
    real(wp),         allocatable :: buf_enm_real(:,:,:)

    integer,          allocatable :: angle_exit(:)
    integer,          allocatable :: angle_add(:)
    integer,          allocatable :: angle_exit_index(:,:)
    integer,          allocatable :: buf_angle_integer(:,:,:)
    real(wp),         allocatable :: buf_angle_real(:,:,:)

    integer,          allocatable :: dihed_exit(:)
    integer,          allocatable :: dihed_add(:)
    integer,          allocatable :: dihed_exit_index(:,:)
    integer,          allocatable :: buf_dihed_integer(:,:,:)
    real(wp),         allocatable :: buf_dihed_real(:,:,:)

    integer,          allocatable :: rb_dihed_exit(:)
    integer,          allocatable :: rb_dihed_add(:)
    integer,          allocatable :: rb_dihed_exit_index(:,:)
    integer,          allocatable :: buf_rb_dihed_integer(:,:,:)
    real(wp),         allocatable :: buf_rb_dihed_real(:,:,:)

    integer,          allocatable :: impr_exit(:)
    integer,          allocatable :: impr_add(:)
    integer,          allocatable :: impr_exit_index(:,:)
    integer,          allocatable :: buf_impr_integer(:,:,:)
    real(wp),         allocatable :: buf_impr_real(:,:,:)

    integer,          allocatable :: cmap_exit(:)
    integer,          allocatable :: cmap_add(:)
    integer,          allocatable :: cmap_exit_index(:,:)
    integer,          allocatable :: buf_cmap_integer(:,:,:)

    integer,          allocatable :: contact_exit(:)
    integer,          allocatable :: contact_add(:)
    integer,          allocatable :: contact_exit_index(:,:)
    integer,          allocatable :: buf_contact_integer(:,:,:)
    real(wp),         allocatable :: buf_contact_real(:,:,:)

    integer,          allocatable :: restraint_exit(:)
    integer,          allocatable :: restraint_add(:)
    integer,          allocatable :: restraint_exit_index(:,:)
    integer,          allocatable :: buf_restraint_integer(:,:)
    real(wp),         allocatable :: buf_restraint_real(:,:,:)

    type(s_table_enefunc)         :: table

    logical                       :: force_switch
    logical                       :: vdw_shift
    logical                       :: vdw_no_switch

    real(wp)                      :: fudge_lj
    real(wp)                      :: fudge_qq
    integer                       :: excl_level

    integer                       :: dispersion_corr
    real(wp)                      :: eswitch
    real(wp)                      :: vswitch
    real(wp)                      :: dispersion_energy
    real(wp)                      :: dispersion_virial
    real(wp)                      :: pme_dispersion_self1
    real(wp)                      :: pme_dispersion_self2

    ! statistical variables
    logical                       :: rpath_flag
    logical                       :: rpath_sum_mf_flag
    integer                       :: rpath_pos_func
    integer                       :: stats_count
    integer                       :: stats_natom
    integer                       :: stats_dimension
    real(wp),         allocatable :: stats_delta(:)
    real(wp),         allocatable :: stats_grad(:, :, :)
    real(dp),         allocatable :: stats_force(:)
    real(dp),         allocatable :: stats_metric(:,:)
    integer,          allocatable :: stats_atom(:,:)
    real(wp),         allocatable :: stats_mass(:,:)
    integer,          allocatable :: stats_id_atm2cv(:)
    integer                       :: stats_icnt
    real(dp),         allocatable :: stats_force_save(:)
    integer,          allocatable :: rpath_rest_function(:)

    logical                       :: contact_check
    logical                       :: nonb_limiter
    logical                       :: pairlist_check
    logical                       :: bonding_check
    real(wp)                      :: minimum_contact
    real(wp)                      :: err_minimum_contact
    logical                       :: pressure_position
    logical                       :: pressure_rmsd

    ! gamd
    logical                       :: gamd_use
    type(s_enefunc_gamd)          :: gamd

    ! For vacuum
    logical                       :: vacuum

    !FEP
    real(wp)                      :: sc_alpha
    real(wp)                      :: sc_beta
    integer                       :: num_fep_neighbor
    integer                       :: fep_direction
    real(wp)                      :: lambljA
    real(wp)                      :: lambljB
    real(wp)                      :: lambelA
    real(wp)                      :: lambelB
    real(wp)                      :: lambbondA
    real(wp)                      :: lambbondB
    real(wp)                      :: lambrest
    integer                       :: fep_topology
    ! table_nonb_lambda is for parameters of lambda, softcore, and exclusion.
    ! table_nonb_lambda(index, fepgrp1, fepgrp2)
    ! index is an index for lambda, softcore, and exclusion (1-5). 1, 2, 3, and 4
    ! correspond to lambdaA, lambdaB, softcore for A, and softcore for B.
    ! 5 represents the excludion of interactions between partA-partB.
    ! fepgrp1 and fepgrp2 are indicies for perturbed states of each atom
    ! pair (1-5). If fepgrp has 1, 2, 3, 4, and 5, the perturbed state of 
    ! an atom correpond to singleA, singleB, dualA, dualB, and preserved,
    ! respectively.
    integer(1)                    :: fep_mask(5,5)
    real(wp)                      :: table_sclj(5,5)
    real(wp)                      :: table_scel(5,5)
    ! table_bond_lambda(fepgrp1,fepgrp2)
    ! table_angl_lambda(fepgrp1,fepgrp2,fepgrp3)
    ! table_dihe_lambda(fepgrp1,fepgrp2,fepgrp3,fepgrp4)
    real(wp)                      :: table_bond_lambda(5,5)
    real(wp)                      :: table_angl_lambda(5,5,5)
    real(wp)                      :: table_dihe_lambda(5,5,5,5)
    integer                       :: fepgrp_nonb(5,5)
    !
    real(wp)                      :: dispersion_energy_CC
    real(wp)                      :: dispersion_energy_AC
    real(wp)                      :: dispersion_energy_BC
    real(wp)                      :: dispersion_energy_AA
    real(wp)                      :: dispersion_energy_BB
    real(wp)                      :: dispersion_energy_AB
    real(wp)                      :: dispersion_virial_CC
    real(wp)                      :: dispersion_virial_AC
    real(wp)                      :: dispersion_virial_BC
    real(wp)                      :: dispersion_virial_AA
    real(wp)                      :: dispersion_virial_BB
    real(wp)                      :: dispersion_virial_AB
    !
    integer                       :: num_nb14_all_fep
    integer,          allocatable :: num_nb14_calc_fep(:)
    integer,          allocatable :: nb14_calc_list_fep(:,:,:)
    integer,          allocatable :: sc_calc_list_fep(:,:)
    real(wp),         allocatable :: nb14_qq_scale_fep(:,:)
    real(wp),         allocatable :: nb14_lj_scale_fep(:,:)
    ! Flag for singleB bonded
    integer,          allocatable :: bond_singleB(:,:)
    integer,          allocatable :: angl_singleB(:,:)
    integer,          allocatable :: dihe_singleB(:,:)
    integer,          allocatable :: rb_dihe_singleB(:,:)
    integer,          allocatable :: impr_singleB(:,:)
    integer,          allocatable :: cmap_singleB(:,:)
  end type s_enefunc

  ! parameter for allocatable variables
  integer,      public, parameter :: EneFuncBase          = 1
  integer,      public, parameter :: EneFuncBond          = 2
  integer,      public, parameter :: EneFuncAngl          = 3
  integer,      public, parameter :: EneFuncDihe          = 4
  integer,      public, parameter :: EneFuncRBDihe        = 5
  integer,      public, parameter :: EneFuncImpr          = 6
  integer,      public, parameter :: EneFuncCmap          = 7
  integer,      public, parameter :: EneFuncNbon          = 8
  integer,      public, parameter :: EneFuncNonb          = 9
  integer,      public, parameter :: EneFuncNonbList      = 10
  integer,      public, parameter :: EneFuncRefg          = 11
  integer,      public, parameter :: EneFuncReff          = 12
  integer,      public, parameter :: EneFuncRefc          = 13
  integer,      public, parameter :: EneFuncRefr          = 14
  integer,      public, parameter :: EneFuncRest          = 15
  integer,      public, parameter :: EneFuncBondCell      = 16
  integer,      public, parameter :: EneFuncTableWat      = 17
  integer,      public, parameter :: EneFuncTableSol      = 18
  integer,      public, parameter :: EneFuncTblWatDomain  = 19
  integer,      public, parameter :: EneFuncTableDomain   = 20
  integer,      public, parameter :: EneFuncAMBERScale    = 21
  integer,      public, parameter :: EneFuncMode          = 22
  integer,      public, parameter :: EneFuncFitc          = 23
  integer,      public, parameter :: EneFuncFitd          = 24
  integer,      public, parameter :: EneFuncRestDomain    = 25
  integer,      public, parameter :: EneFuncContact       = 26
  integer,      public, parameter :: EneFuncBase_pio      = 27
  integer,      public, parameter :: EneFuncBond_pio      = 28
  integer,      public, parameter :: EneFuncAngl_pio      = 29
  integer,      public, parameter :: EneFuncDihe_pio      = 30
  integer,      public, parameter :: EneFuncRBDihe_pio    = 31
  integer,      public, parameter :: EneFuncImpr_pio      = 32
  integer,      public, parameter :: EneFuncCmap_pio      = 33
  integer,      public, parameter :: EneFuncCmapcoef_pio  = 34
  integer,      public, parameter :: EneFuncRest_pio      = 35
  integer,      public, parameter :: EneFuncGamdDih       = 36
  integer,      public, parameter :: EneFuncGamdRest      = 37
  integer,      public, parameter :: EneFuncEnm           = 38
  integer,      public, parameter :: EneFuncExcl          = 39
  integer,      public, parameter :: EneFuncNonbFEPList   = 46
  integer,      public, parameter :: EneFuncFEPBonded     = 47
  integer,      public, parameter :: EneFuncBondCell_FEP  = 48

  ! water model type
  integer,      public, parameter :: TIP1                 = 1
  integer,      public, parameter :: TIP2                 = 2
  integer,      public, parameter :: TIP3                 = 3
  integer,      public, parameter :: TIP4                 = 4

  ! parameters (forcefield)
  integer,      public, parameter :: ForcefieldCHARMM     = 1
  integer,      public, parameter :: ForcefieldAMBER      = 2
  integer,      public, parameter :: ForcefieldGROAMBER   = 3
  integer,      public, parameter :: ForcefieldGROMARTINI = 4
  character(*), public, parameter :: ForceFieldTypes(4)   = (/'CHARMM    ', &
                                                              'AMBER     ', &
                                                              'GROAMBER  ', &
                                                              'GROMARTINI'/)
  ! FFT scheme
  integer,      public, parameter :: FFT_AutoSelect       = 1
  integer,      public, parameter :: FFT_opt_1dalltoall   = 2
  integer,      public, parameter :: FFT_noopt_1dalltoall = 3
  integer,      public, parameter :: FFT_opt_2dalltoall   = 4
  integer,      public, parameter :: FFT_noopt_2dalltoall = 5
  integer,      public, parameter :: FFT_opt_slab         = 6
  integer,      public, parameter :: FFT_noopt_slab       = 7
  character(*), public, parameter :: FFT_Types(7)         = &
                                                        (/'AUTOSELECT      ',&
                                                          'OPT_1DALLTOALL  ',&
                                                          'NOOPT_1DALLTOALL',&
                                                          'OPT_2DALLTOALL  ',&
                                                          'NOOPT_2DALLTOALL',&
                                                          'OPT_SLAB        ',&
                                                          'NOOPT_SLAB      '/)

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

  ! parameters (FEP calculation)
  integer,      public, parameter :: FEP_PRESERVE         = 0
  integer,      public, parameter :: FEP_APPEAR           = 1
  integer,      public, parameter :: FEP_VANISH           = -1

  ! variables for maximum numbers in one cell (these number will be updated)
  integer,      public            :: MaxBond   = 0
  integer,      public            :: MaxAngle  = 0
  integer,      public            :: MaxDihe   = 0
  integer,      public            :: MaxImpr   = 0
  integer,      public            :: MaxCmap   = 0
  integer,      public            :: MaxRest   = 100
  integer,      public            :: MaxExcl   = 36
  integer,      public            :: MaxContact= 0

  integer,      public            :: BondMove  = 0
  integer,      public            :: AngleMove = 0
  integer,      public            :: DiheMove  = 0
  integer,      public            :: ImprMove  = 0
  integer,      public            :: CmapMove  = 10
  integer,      public            :: RestMove  = 50
  integer,      public            :: ContactMove  = 0

  integer,      public, parameter :: MaxAtomCls = 1000
  integer,      public            :: max_class
! real(wp),     public            :: lj_coef(2,MaxAtomCls)


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
    enefunc%dielec_const           = 0.0_wp

    enefunc%enm_use                = .false.

    enefunc%pme_use                = .false.
    enefunc%pme_alpha              = 0.0_wp
    enefunc%pme_ngrid_x            = 0
    enefunc%pme_ngrid_y            = 0
    enefunc%pme_ngrid_z            = 0
    enefunc%pme_nspline            = 0
    enefunc%pme_max_spacing        = 1.2_wp

    enefunc%use_efield             = .false.
    enefunc%efield(1:3)            = 0.0_wp
    enefunc%efield_normal          = .false.
    enefunc%efield_virial          = .false.

    enefunc%restraint_posi         = .false.
    enefunc%restraint_rmsd         = .false.
    enefunc%restraint_rmsd_target  = .false.
    enefunc%restraint_emfit        = .false.
    enefunc%restraint_pc           = .false.
    enefunc%local_restraint        = .false.
    
    enefunc%table%table            = .false.
    enefunc%table%water_table      = .false.
    enefunc%table%atom_cls_no_O    = 0
    enefunc%table%atom_cls_no_H    = 0
    enefunc%table%table_order      = 0
    enefunc%table%density          = 0.0_wp
    enefunc%table%charge_O         = 0.0_wp
    enefunc%table%charge_H         = 0.0_wp
    enefunc%table%mass_O           = 0.0_wp
    enefunc%table%mass_H           = 0.0_wp
    enefunc%table%water_model      = ''
    enefunc%table%num_water        = 0
    enefunc%table%num_solute       = 0

    enefunc%force_switch           = .false.
    enefunc%vdw_shift              = .false.
    enefunc%vdw_no_switch          = .false.

    enefunc%fudge_lj                = 1.0_wp
    enefunc%fudge_qq                = 1.0_wp
    enefunc%excl_level              = 3

    enefunc%notation_14types        = 0

    enefunc%rpath_flag              = .false.
    enefunc%rpath_sum_mf_flag       = .false.
    enefunc%rpath_pos_func          = 0
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

    enefunc%gamd_use                = .false.

    return

  end subroutine init_enefunc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_enefunc
  !> @brief        allocate energy functions information
  !! @authors      YS, CK
  !! @param[out]   enefunc   : potential energy functions information
  !! @param[in]    variable  : selected variable
  !! @param[in]    var_size  : size of the selected variable
  !! @param[in]    var_size1 : 2nd size of the selected variable (optional)
  !! @param[in]    var_size2 : 3rd size of the selected variable (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_enefunc(enefunc, variable, var_size, var_size1, var_size2)

    ! formal arguments
    type(s_enefunc),         intent(inout) :: enefunc
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size
    integer,       optional, intent(in)    :: var_size1
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

    case (EneFuncBase)

      if (allocated(enefunc%num_bond)) then
        if (size(enefunc%num_bond(:)) /= var_size) &
          deallocate(enefunc%num_bond,        &
                     enefunc%num_angle,       &
                     enefunc%num_dihedral,    &
                     enefunc%num_rb_dihedral, &
                     enefunc%num_improper,    &
                     enefunc%num_cmap,        &
                     enefunc%num_contact,     &
                     enefunc%num_restraint,   &
                     enefunc%num_excl,        &
                     enefunc%num_excl14,      &
                     enefunc%nfitting)
      end if

      if (.not. allocated(enefunc%num_bond)) &
        allocate(enefunc%num_bond       (var_size), &
                 enefunc%num_angle      (var_size), &
                 enefunc%num_dihedral   (var_size), &
                 enefunc%num_rb_dihedral(var_size), &
                 enefunc%num_improper   (var_size), &
                 enefunc%num_cmap       (var_size), &
                 enefunc%num_contact    (var_size), &
                 enefunc%num_restraint  (var_size), &
                 enefunc%nfitting       (var_size), &
                 enefunc%num_excl       (var_size), &
                 enefunc%num_excl14     (var_size), &
                 stat = alloc_stat)

      enefunc%num_bond       (1:var_size) = 0
      enefunc%num_angle      (1:var_size) = 0
      enefunc%num_dihedral   (1:var_size) = 0
      enefunc%num_rb_dihedral(1:var_size) = 0
      enefunc%num_improper   (1:var_size) = 0
      enefunc%num_cmap       (1:var_size) = 0
      enefunc%num_contact    (1:var_size) = 0
      enefunc%num_restraint  (1:var_size) = 0
      enefunc%nfitting       (1:var_size) = 0
      enefunc%num_excl       (1:var_size) = 0
      enefunc%num_excl14     (1:var_size) = 0

    case (EneFuncBase_pio)

      if (allocated(enefunc%num_bond_pio)) then
        if (size(enefunc%num_bond_pio) /= var_size*var_size1) &
          deallocate(enefunc%num_bond_pio,        &
                     enefunc%num_angle_pio,       &
                     enefunc%num_dihedral_pio,    &
                     enefunc%num_rb_dihedral_pio, &
                     enefunc%num_improper_pio,    &
                     enefunc%num_cmap_pio,        &
                     enefunc%num_restraint_pio)
      end if

      if (.not. allocated(enefunc%num_bond_pio)) &
        allocate(enefunc%num_bond_pio       (var_size, var_size1), &
                 enefunc%num_angle_pio      (var_size, var_size1), &
                 enefunc%num_dihedral_pio   (var_size, var_size1), &
                 enefunc%num_rb_dihedral_pio(var_size, var_size1), &
                 enefunc%num_improper_pio   (var_size, var_size1), &
                 enefunc%num_cmap_pio       (var_size, var_size1), &
                 enefunc%num_restraint_pio  (var_size, var_size1), &
                 stat = alloc_stat)

      enefunc%num_bond_pio       (1:var_size, 1:var_size1) = 0
      enefunc%num_angle_pio      (1:var_size, 1:var_size1) = 0
      enefunc%num_dihedral_pio   (1:var_size, 1:var_size1) = 0
      enefunc%num_rb_dihedral_pio(1:var_size, 1:var_size1) = 0
      enefunc%num_improper_pio   (1:var_size, 1:var_size1) = 0
      enefunc%num_cmap_pio       (1:var_size, 1:var_size1) = 0
      enefunc%num_restraint_pio  (1:var_size, 1:var_size1) = 0

    case (EneFuncBond)

      if (allocated(enefunc%bond_list)) then
        if (size(enefunc%bond_list(1,1,:)) /= var_size) &
          deallocate(enefunc%bond_list,        &
                     enefunc%bond_force_const, &
                     enefunc%bond_dist_min,    &
                     enefunc%bond_kind,        &
                     enefunc%bond_pbc,         &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%bond_list)) &
        allocate(enefunc%bond_list       (2, MaxBond, var_size), &
                 enefunc%bond_force_const(MaxBond, var_size),    &
                 enefunc%bond_dist_min   (MaxBond, var_size),    &
                 enefunc%bond_kind       (MaxBond, var_size),    &
                 enefunc%bond_pbc        (MaxBond, var_size),    &
                 stat = alloc_stat)

      enefunc%bond_list       (1:2, 1:MaxBond, 1:var_size) = 0
      enefunc%bond_force_const(1:MaxBond, 1:var_size)      = 0.0_wp
      enefunc%bond_dist_min   (1:MaxBond, 1:var_size)      = 0.0_wp
      enefunc%bond_kind       (1:MaxBond, 1:var_size)      = 0
      enefunc%bond_pbc        (1:MaxBond, 1:var_size)      = 13

    case (EneFuncBond_pio)

      if (allocated(enefunc%bond_list_pio)) then
        if (size(enefunc%bond_list_pio(1,1,:,:)) /= var_size*var_size1) &
          deallocate(enefunc%bond_list_pio,        &
                     enefunc%bond_force_const_pio, &
                     enefunc%bond_dist_min_pio,    &
                     enefunc%bond_pbc_pio,         &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%bond_list_pio)) &
        allocate(enefunc%bond_list_pio       (2, MaxBond, var_size, var_size1),&
                 enefunc%bond_force_const_pio(MaxBond, var_size, var_size1),   &
                 enefunc%bond_dist_min_pio   (MaxBond, var_size, var_size1),   &
                 enefunc%bond_pbc_pio        (MaxBond, var_size, var_size1),   &
                 stat = alloc_stat)

      enefunc%bond_list_pio  (1:2, 1:MaxBond, 1:var_size, 1:var_size1) = 0
      enefunc%bond_force_const_pio(1:MaxBond, 1:var_size, 1:var_size1) = 0.0_wp
      enefunc%bond_dist_min_pio   (1:MaxBond, 1:var_size, 1:var_size1) = 0.0_wp
      enefunc%bond_pbc_pio        (1:MaxBond, 1:var_size, 1:var_size1) = 0.0_wp

    case (EneFuncEnm)

      if (allocated(enefunc%num_enm)) then
        if (size(enefunc%num_enm(:)) /= var_size) &
          deallocate(enefunc%num_enm,             &
                     enefunc%enm_list,            &
                     enefunc%enm_force_const,     &
                     enefunc%enm_dist_min,        &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%num_enm)) &
        allocate(enefunc%num_enm        (var_size),             &
                 enefunc%enm_list       (2, MaxBond, var_size), &
                 enefunc%enm_force_const(MaxBond, var_size),    &
                 enefunc%enm_dist_min   (MaxBond, var_size),    &
                 stat = alloc_stat)

      enefunc%num_enm        (1:var_size)                 = 0
      enefunc%enm_list       (1:2, 1:MaxBond, 1:var_size) = 0
      enefunc%enm_force_const(1:MaxBond, 1:var_size)      = 0.0_wp
      enefunc%enm_dist_min   (1:MaxBond, 1:var_size)      = 0.0_wp

    case (EneFuncExcl)

      if (allocated(enefunc%excl_list)) then
        if (size(enefunc%excl_list(1,1,:)) /= var_size) &
          deallocate(enefunc%excl_list,        &
                     enefunc%excl14_list,      &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%excl_list)) &
        allocate(enefunc%excl_list       (2, MaxAngle, var_size), &
                 enefunc%excl14_list     (2, MaxDihe , var_size), &
                 stat = alloc_stat)

      enefunc%excl_list  (1:2, 1:MaxAngle, 1:var_size) = 0
      enefunc%excl14_list(1:2, 1:MaxDihe , 1:var_size) = 0

    case (EneFuncContact)

      if (allocated(enefunc%contact_list)) then
        if (size(enefunc%contact_list(1,1,:)) /= var_size) &
          deallocate(enefunc%contact_list,   &
                     enefunc%contact_lj12,   &
                     enefunc%contact_lj6,    &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%contact_list)) &
        allocate(enefunc%contact_list(2, MaxContact, var_size), &
                 enefunc%contact_lj12(MaxContact, var_size),    &
                 enefunc%contact_lj6(MaxContact, var_size),     &
                 stat = alloc_stat)

      enefunc%contact_list(1:2, 1:MaxContact, 1:var_size) = 0
      enefunc%contact_lj12(1:MaxContact, 1:var_size)      = 0.0_wp
      enefunc%contact_lj6 (1:MaxContact, 1:var_size)      = 0.0_wp

    case (EneFuncAngl)

      if (allocated(enefunc%angle_list)) then
        if (size(enefunc%angle_list(1,1,:)) /= var_size) &
          deallocate(enefunc%angle_list,        &
                     enefunc%angle_force_const, &
                     enefunc%angle_theta_min,   &
                     enefunc%urey_force_const,  &
                     enefunc%urey_rmin,         &
                     enefunc%angle_kind,        &
                     enefunc%angle_pbc,         &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%angle_list)) &
        allocate(enefunc%angle_list       (3, MaxAngle, var_size),    &
                 enefunc%angle_force_const(   MaxAngle, var_size),    &
                 enefunc%angle_theta_min  (   MaxAngle, var_size),    &
                 enefunc%urey_force_const (   MaxAngle, var_size),    &
                 enefunc%urey_rmin        (   MaxAngle, var_size),    &
                 enefunc%angle_kind       (   MaxAngle, var_size),    &
                 enefunc%angle_pbc        (3, MaxAngle, var_size),    &
                 stat = alloc_stat)

      enefunc%angle_list       (1:3, 1:MaxAngle, 1:var_size) = 0
      enefunc%angle_force_const(     1:MaxAngle, 1:var_size) = 0.0_wp
      enefunc%angle_theta_min  (     1:MaxAngle, 1:var_size) = 0.0_wp
      enefunc%urey_force_const (     1:MaxAngle, 1:var_size) = 0.0_wp
      enefunc%urey_rmin        (     1:MaxAngle, 1:var_size) = 0.0_wp
      enefunc%angle_kind       (     1:MaxAngle, 1:var_size) = 0
      enefunc%angle_pbc        (1:3, 1:MaxAngle, 1:var_size) = 13

    case (EneFuncAngl_pio)

      if (allocated(enefunc%angle_list_pio)) then
        if (size(enefunc%angle_list_pio(1,1,:,:)) /= var_size*var_size1) &
          deallocate(enefunc%angle_list_pio,        &
                     enefunc%angle_force_const_pio, &
                     enefunc%angle_theta_min_pio,   &
                     enefunc%urey_force_const_pio,  &
                     enefunc%urey_rmin_pio,         &
                     enefunc%angle_pbc_pio,         &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%angle_list_pio)) &
        allocate(enefunc%angle_list_pio    (3, MaxAngle, var_size, var_size1), &
                 enefunc%angle_force_const_pio(MaxAngle, var_size, var_size1), &
                 enefunc%angle_theta_min_pio  (MaxAngle, var_size, var_size1), &
                 enefunc%urey_force_const_pio (MaxAngle, var_size, var_size1), &
                 enefunc%urey_rmin_pio        (MaxAngle, var_size, var_size1), &
                 enefunc%angle_pbc_pio     (3, MaxAngle, var_size, var_size1), &
                 stat = alloc_stat)

      enefunc%angle_list_pio  (1:3, 1:MaxAngle, 1:var_size, 1:var_size1)= 0
      enefunc%angle_force_const_pio(1:MaxAngle, 1:var_size, 1:var_size1)= 0.0_wp
      enefunc%angle_theta_min_pio  (1:MaxAngle, 1:var_size, 1:var_size1)= 0.0_wp
      enefunc%urey_force_const_pio (1:MaxAngle, 1:var_size, 1:var_size1)= 0.0_wp
      enefunc%urey_rmin_pio        (1:MaxAngle, 1:var_size, 1:var_size1)= 0.0_wp
      enefunc%angle_pbc_pio   (1:3, 1:MaxAngle, 1:var_size, 1:var_size1)= 13

    case (EneFuncDihe)

      if (allocated(enefunc%dihe_list)) then
        if (size(enefunc%dihe_list(1,1,:)) /= var_size) &
          deallocate(enefunc%dihe_list,        &
                     enefunc%dihe_force_const, &
                     enefunc%dihe_periodicity, &
                     enefunc%dihe_phase,       &
                     enefunc%dihe_kind,        &
                     enefunc%dihe_pbc,         &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%dihe_list)) &
        allocate(enefunc%dihe_list       (4, MaxDihe, var_size),    &
                 enefunc%dihe_force_const(   MaxDihe, var_size),    &
                 enefunc%dihe_periodicity(   MaxDihe, var_size),    &
                 enefunc%dihe_phase      (   MaxDihe, var_size),    &
                 enefunc%dihe_kind       (   MaxDihe, var_size),    &
                 enefunc%dihe_pbc        (3, MaxDihe, var_size),    &
                 stat = alloc_stat)

      enefunc%dihe_list       (1:4, 1:MaxDihe, 1:var_size) = 0
      enefunc%dihe_force_const(1:MaxDihe, 1:var_size)      = 0.0_wp
      enefunc%dihe_periodicity(1:MaxDihe, 1:var_size)      = 0
      enefunc%dihe_phase      (1:MaxDihe, 1:var_size)      = 0.0_wp
      enefunc%dihe_kind       (1:MaxDihe, 1:var_size)      = 0

    case (EneFuncDihe_pio)

      if (allocated(enefunc%dihe_list_pio)) then
        if (size(enefunc%dihe_list_pio(1,1,:,:)) /= var_size*var_size1) &
          deallocate(enefunc%dihe_list_pio,        &
                     enefunc%dihe_force_const_pio, &
                     enefunc%dihe_periodicity_pio, &
                     enefunc%dihe_phase_pio,       &
                     enefunc%dihe_pbc_pio,         &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%dihe_list_pio)) &
        allocate(enefunc%dihe_list_pio    (4, MaxDihe, var_size, var_size1), &
                 enefunc%dihe_force_const_pio(MaxDihe, var_size, var_size1), &
                 enefunc%dihe_periodicity_pio(MaxDihe, var_size, var_size1), &
                 enefunc%dihe_phase_pio      (MaxDihe, var_size, var_size1), &
                 enefunc%dihe_pbc_pio     (3, MaxDihe, var_size, var_size1), &
                 stat = alloc_stat)

      enefunc%dihe_list_pio  (1:4, 1:MaxDihe, 1:var_size, 1:var_size1) = 0
      enefunc%dihe_force_const_pio(1:MaxDihe, 1:var_size, 1:var_size1) = 0.0_wp
      enefunc%dihe_periodicity_pio(1:MaxDihe, 1:var_size, 1:var_size1) = 0
      enefunc%dihe_phase_pio      (1:MaxDihe, 1:var_size, 1:var_size1) = 0.0_wp
      enefunc%dihe_pbc_pio   (1:3, 1:MaxDihe, 1:var_size, 1:var_size1) = 13

    case (EneFuncRBDihe)

      if (allocated(enefunc%rb_dihe_list)) then

        if (size(enefunc%rb_dihe_list(1,1,:)) /= var_size) &
          deallocate(enefunc%rb_dihe_list, &
                     enefunc%rb_dihe_c,    &
                     enefunc%rb_dihe_pbc,  &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%rb_dihe_list)) &
        allocate(enefunc%rb_dihe_list    (4, MaxDihe, var_size), &
                 enefunc%rb_dihe_c       (6, MaxDihe, var_size), &
                 enefunc%rb_dihe_pbc     (3, MaxDihe, var_size), &
                 stat = alloc_stat)

      enefunc%rb_dihe_list(1:4, 1:MaxDihe, 1:var_size) = 0
      enefunc%rb_dihe_c   (1:6, 1:MaxDihe, 1:var_size) = 0.0_wp

    case (EneFuncRBDihe_pio)

      if (allocated(enefunc%rb_dihe_list_pio)) then

        if (size(enefunc%rb_dihe_list_pio(1,1,:,:)) /= var_size*var_size1) &
          deallocate(enefunc%rb_dihe_list_pio, &
                     enefunc%rb_dihe_c_pio,    &
                     enefunc%rb_dihe_pbc_pio,  &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%rb_dihe_list_pio)) &
        allocate(enefunc%rb_dihe_list_pio(4, MaxDihe, var_size, var_size1), &
                 enefunc%rb_dihe_c_pio   (6, MaxDihe, var_size, var_size1), &
                 enefunc%rb_dihe_pbc_pio (3, MaxDihe, var_size, var_size1), &
                 stat = alloc_stat)

      enefunc%rb_dihe_list_pio(1:4, 1:MaxDihe, 1:var_size, 1:var_size1) = 0
      enefunc%rb_dihe_c_pio   (1:6, 1:MaxDihe, 1:var_size, 1:var_size1) = 0.0_wp
      enefunc%rb_dihe_pbc_pio (1:3, 1:MaxDihe, 1:var_size, 1:var_size1) = 13

    case (EneFuncImpr)

      if (allocated(enefunc%impr_list)) then
        if (size(enefunc%impr_list(1,1,:)) /= var_size) &
          deallocate(enefunc%impr_list,        &
                     enefunc%impr_force_const, &
                     enefunc%impr_periodicity, &
                     enefunc%impr_phase,       &
                     enefunc%impr_pbc,         &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%impr_list)) &
        allocate(enefunc%impr_list       (4, MaxImpr, var_size), &
                 enefunc%impr_force_const(   MaxImpr, var_size), &
                 enefunc%impr_periodicity(   MaxImpr, var_size), &
                 enefunc%impr_phase      (   MaxImpr, var_size), &
                 enefunc%impr_pbc        (3, MaxImpr, var_size), &
                 stat = alloc_stat)

      enefunc%impr_list       (1:4, 1:MaxImpr, 1:var_size) = 0
      enefunc%impr_force_const(1:MaxImpr, 1:var_size)      = 0.0_wp
      enefunc%impr_periodicity(1:MaxImpr, 1:var_size)      = 0
      enefunc%impr_phase      (1:MaxImpr, 1:var_size)      = 0.0_wp

    case (EneFuncImpr_pio)

      if (allocated(enefunc%impr_list_pio)) then
        if (size(enefunc%impr_list_pio(1,1,:,:)) /= var_size*var_size1) &
          deallocate(enefunc%impr_list_pio,        &
                     enefunc%impr_force_const_pio, &
                     enefunc%impr_periodicity_pio, &
                     enefunc%impr_phase_pio,       &
                     enefunc%impr_pbc_pio,         &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%impr_list_pio)) &
        allocate(enefunc%impr_list_pio    (4, MaxImpr, var_size, var_size1), &
                 enefunc%impr_force_const_pio(MaxImpr, var_size, var_size1), &
                 enefunc%impr_periodicity_pio(MaxImpr, var_size, var_size1), &
                 enefunc%impr_phase_pio      (MaxImpr, var_size, var_size1), &
                 enefunc%impr_pbc_pio     (3, MaxImpr, var_size, var_size1), &
                 stat = alloc_stat)

      enefunc%impr_list_pio  (1:4, 1:MaxImpr, 1:var_size, 1:var_size1) = 0
      enefunc%impr_force_const_pio(1:MaxImpr, 1:var_size, 1:var_size1) = 0.0_wp
      enefunc%impr_periodicity_pio(1:MaxImpr, 1:var_size, 1:var_size1) = 0
      enefunc%impr_phase_pio      (1:MaxImpr, 1:var_size, 1:var_size1) = 0.0_wp
      enefunc%impr_pbc_pio   (1:3, 1:MaxImpr, 1:var_size, 1:var_size1) = 13

    case (EneFuncCmap)

      if (allocated(enefunc%cmap_resolution)) then
        if (size(enefunc%cmap_list(1,1,:))   /= var_size .or. &
            size(enefunc%cmap_resolution(:)) /= var_size2)    &
          deallocate(enefunc%cmap_resolution, &
                     enefunc%cmap_list,       &
                     enefunc%cmap_type,       &
                     enefunc%cmap_coef,       &
                     enefunc%cmap_pbc,        &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%cmap_resolution))      &
        allocate(enefunc%cmap_resolution(var_size2),     &
                 enefunc%cmap_list(8,MaxCmap, var_size), &
                 enefunc%cmap_type(  MaxCmap, var_size), &
                 enefunc%cmap_coef(4, 4, var_size1, var_size1, var_size2), &
                 enefunc%cmap_pbc (6,MaxCmap, var_size), &
                 stat = alloc_stat)

      enefunc%cmap_resolution(1:var_size2)          = 0
      enefunc%cmap_list(1:8, 1:MaxCmap, 1:var_size) = 0
      enefunc%cmap_type(1:MaxCmap, 1:var_size)      = 0
      enefunc%cmap_coef(1:4, 1:4, 1:var_size1, 1:var_size1, 1:var_size2) = &
                                                  0.0_wp
      enefunc%cmap_pbc (1:6, 1:MaxCmap, 1:var_size) = 13

    case (EneFuncCmapcoef_pio)

      if (allocated(enefunc%cmap_resolution_pio)) then
        if (size(enefunc%cmap_resolution_pio) /= var_size1)  &
          deallocate(enefunc%cmap_resolution_pio, &
                     enefunc%cmap_coef_pio,       &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%cmap_resolution_pio)) &
        allocate(enefunc%cmap_resolution_pio(var_size1),     &
                 enefunc%cmap_coef_pio(4,4,var_size,var_size,var_size1), &
                 stat = alloc_stat)

      enefunc%cmap_resolution_pio(1:var_size1)                         = 0
      enefunc%cmap_coef_pio(1:4,1:4,1:var_size,1:var_size,1:var_size1) = 0.0_wp

    case (EneFuncCmap_pio)

      if (allocated(enefunc%cmap_list_pio)) then
        if (size(enefunc%cmap_list_pio(1,1,:,:))   /= var_size*var_size1)  &
          deallocate(enefunc%cmap_list_pio,       &
                     enefunc%cmap_type_pio,       &
                     enefunc%cmap_pbc_pio,        &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%cmap_list_pio)) &
        allocate(enefunc%cmap_list_pio(8,MaxCmap, var_size, var_size1), &
                 enefunc%cmap_type_pio(MaxCmap, var_size, var_size1),   &
                 enefunc%cmap_pbc_pio (6,MaxCmap, var_size, var_size1), &
                 stat = alloc_stat)

      enefunc%cmap_list_pio(1:8, 1:MaxCmap, 1:var_size, 1:var_size1) = 0
      enefunc%cmap_type_pio(1:MaxCmap, 1:var_size, 1:var_size1)      = 0
      enefunc%cmap_pbc_pio (1:6, 1:MaxCmap, 1:var_size, 1:var_size1) = 13
    case (EneFuncFEPBonded)

      if (allocated(enefunc%bond_singleB)) then
        if (size(enefunc%bond_singleB(1,:)) /= var_size) &
          deallocate(enefunc%bond_singleB,        &
                     enefunc%angl_singleB,        &
                     enefunc%dihe_singleB,        &
                     enefunc%rb_dihe_singleB,     &
                     enefunc%impr_singleB,        &
                     enefunc%cmap_singleB,        &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%bond_singleB)) &
        allocate(enefunc%bond_singleB   (MaxBond, var_size),  &
                 enefunc%angl_singleB   (MaxAngle, var_size), &
                 enefunc%dihe_singleB   (MaxDihe, var_size),  &
                 enefunc%rb_dihe_singleB(MaxDihe, var_size),  &
                 enefunc%impr_singleB   (MaxImpr, var_size),  &
                 enefunc%cmap_singleB   (MaxCmap, var_size),  &
                 stat = alloc_stat)

      enefunc%bond_singleB   (1:MaxBond, 1:var_size)  = 0
      enefunc%angl_singleB   (1:MaxAngle, 1:var_size) = 0
      enefunc%dihe_singleB   (1:MaxDihe, 1:var_size)  = 0
      enefunc%rb_dihe_singleB(1:MaxDihe, 1:var_size)  = 0
      enefunc%impr_singleB   (1:MaxImpr, 1:var_size)  = 0
      enefunc%cmap_singleB   (1:MaxCmap, 1:var_size)  = 0

    case (EneFuncNbon)

      if (allocated(enefunc%nonb_atom_cls)) then
        if (size(enefunc%nonb_atom_cls(:)) /= var_size) &
          deallocate(enefunc%nonb_atom_cls,     &
                     enefunc%nonb_atom_cls_no,  &
                     enefunc%nonb_ljtype,       &
                     enefunc%nb14_lj6,          &
                     enefunc%nb14_lj12,         &
                     enefunc%nonb_lj6,          &
                     enefunc%nonb_lj12,         &
                     enefunc%nonb_lj6_factor,   &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%nonb_atom_cls)) &
        allocate(enefunc%nonb_atom_cls   (var_size),           &
                 enefunc%nonb_atom_cls_no(var_size),           &
                 enefunc%nonb_ljtype     (var_size, var_size), &
                 enefunc%nb14_lj6        (var_size, var_size), &
                 enefunc%nb14_lj12       (var_size, var_size), &
                 enefunc%nonb_lj6        (var_size, var_size), &
                 enefunc%nonb_lj12       (var_size, var_size), &
                 enefunc%nonb_lj6_factor (var_size),           &
                 stat = alloc_stat)

      enefunc%nonb_atom_cls   (1:var_size)             = 0
      enefunc%nonb_atom_cls_no(1:var_size)             = 0
      enefunc%nonb_ljtype     (1:var_size, 1:var_size) = 0
      enefunc%nb14_lj6        (1:var_size, 1:var_size) = 0.0_wp
      enefunc%nb14_lj12       (1:var_size, 1:var_size) = 0.0_wp
      enefunc%nonb_lj6        (1:var_size, 1:var_size) = 0.0_wp
      enefunc%nonb_lj12       (1:var_size, 1:var_size) = 0.0_wp
      enefunc%nonb_lj6_factor (1:var_size)             = 0.0_wp

    case (EneFuncNonb)

      if (allocated(enefunc%num_nonb_excl)) then
        if (size(enefunc%num_nonb_excl(:)) /= var_size) &
          deallocate(enefunc%num_nonb_excl,   &
                     enefunc%num_nb14_calc,   &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%num_nonb_excl))    &
        allocate(enefunc%num_nonb_excl(var_size),  &
                 enefunc%num_nb14_calc(var_size), &
                 stat = alloc_stat)

      enefunc%num_nonb_excl(1:var_size) = 0
      enefunc%num_nb14_calc(1:var_size) = 0

    case (EneFuncNonbList)

      if (allocated(enefunc%nonb_excl_list)) then
        if (size(enefunc%nonb_excl_list(1,1,:)) /= var_size) &
          deallocate(enefunc%nonb_excl_list,  &
                     enefunc%nb14_calc_list,  &
                     enefunc%sc_calc_list,    &
                     enefunc%nb14_qq_scale,   &
                     enefunc%nb14_lj_scale,   &
                     enefunc%exclusion_mask,  &
                     enefunc%exclusion_mask1, &
                     stat = dealloc_stat)
      end if

      var_size3 = MaxBond+MaxAngle*2+MaxWater*3
      var_size4 = MaxDihe*2+MaxImpr
      if (.not. allocated(enefunc%nonb_excl_list)) &
        allocate(enefunc%nonb_excl_list (4, var_size3, var_size),   &
                 enefunc%nb14_calc_list (4, var_size4, var_size),   &
                 enefunc%sc_calc_list   (var_size4, var_size),      &
                 enefunc%nb14_qq_scale  (var_size4, var_size),      &
                 enefunc%nb14_lj_scale  (var_size4, var_size),      &
                 enefunc%exclusion_mask(MaxAtom,MaxAtom,var_size1), &
                 enefunc%exclusion_mask1(MaxAtom,MaxAtom,var_size), &
                 stat = alloc_stat)

      enefunc%nonb_excl_list (1:4, 1:var_size3, 1:var_size) = 0
      enefunc%nb14_calc_list (1:4, 1:var_size4, 1:var_size) = 0
      enefunc%sc_calc_list   (1:var_size4, 1:var_size)      = 0
      enefunc%nb14_qq_scale  (1:var_size4, 1:var_size)      = 0.0_wp
      enefunc%nb14_lj_scale  (1:var_size4, 1:var_size)      = 0.0_wp

    case (EneFuncNonbFEPList)

      if (allocated(enefunc%nb14_calc_list_fep)) then
        if (size(enefunc%nb14_calc_list_fep(1,1,:)) /= var_size) &
          deallocate(enefunc%nb14_calc_list_fep,  &
                     enefunc%sc_calc_list_fep,    &
                     enefunc%nb14_qq_scale_fep,   &
                     enefunc%nb14_lj_scale_fep,   &
                     enefunc%num_nb14_calc_fep,   &
                     stat = dealloc_stat)
      end if

      var_size3 = MaxDihe*2+MaxImpr
      if (.not. allocated(enefunc%nb14_calc_list_fep)) &
        allocate(enefunc%nb14_calc_list_fep (4, var_size3, var_size),   &
                 enefunc%sc_calc_list_fep   (var_size3, var_size),      &
                 enefunc%nb14_qq_scale_fep  (var_size3, var_size),      &
                 enefunc%nb14_lj_scale_fep  (var_size3, var_size),      &
                 enefunc%num_nb14_calc_fep  (var_size),             &
                 stat = alloc_stat)

      enefunc%nb14_calc_list_fep (1:4, 1:var_size3, 1:var_size) = 0
      enefunc%sc_calc_list_fep   (1:var_size3, 1:var_size)      = 0
      enefunc%nb14_qq_scale_fep  (1:var_size3, 1:var_size)      = 0.0_wp
      enefunc%nb14_lj_scale_fep  (1:var_size3, 1:var_size)      = 0.0_wp
      enefunc%num_nb14_calc_fep(1:var_size)                     = 0

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
                     enefunc%restraint_rpath_func,    &
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
                 enefunc%restraint_rpath_func   (1:var_size),              &
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
      enefunc%restraint_rpath_func   (1:var_size)              = 0

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
        if (size(enefunc%restraint_atom(1,:)) /= var_size) &
          deallocate(enefunc%restraint_atom,                  &
                     enefunc%restraint_force,                 &
                     enefunc%restraint_coord,                 &
                     enefunc%restraint_bondslist_to_atomlist, &
                     enefunc%restraint_bonds_coord,           &
                     enefunc%restraint_bonds_force,           &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%restraint_atom)) &
        allocate(enefunc%restraint_atom (MaxRest, var_size),         &
                 enefunc%restraint_force(4, MaxRest, var_size),      &
                 enefunc%restraint_coord(3, MaxRest, var_size),      &
                 enefunc%restraint_bondslist_to_atomlist(var_size1), &
                 enefunc%restraint_bonds_coord(3, var_size1),        &
                 enefunc%restraint_bonds_force(3, var_size1),        &
                 stat = alloc_stat)

      enefunc%restraint_atom (1:MaxRest, 1:var_size)       = 0
      enefunc%restraint_force(1:4, 1:MaxRest, 1:var_size)  = 0.0_wp
      enefunc%restraint_coord(1:3, 1:MaxRest, 1:var_size)  = 0.0_wp
      enefunc%restraint_bondslist_to_atomlist(1:var_size1) = 0
      enefunc%restraint_bonds_coord(1:3, 1:var_size1)      = 0.0_wp
      enefunc%restraint_bonds_force(1:3, 1:var_size1)      = 0.0_wp

    case (EneFuncRest_pio)

      if (allocated(enefunc%restraint_atom_pio)) then
        if (size(enefunc%restraint_atom_pio) /= var_size*var_size1) &
          deallocate(enefunc%restraint_atom_pio,                  &
                     enefunc%restraint_force_pio,                 &
                     enefunc%restraint_coord_pio,                 &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%restraint_atom_pio)) &
        allocate(enefunc%restraint_atom_pio    (MaxRest, var_size, var_size1), &
                 enefunc%restraint_force_pio(4, MaxRest, var_size, var_size1), &
                 enefunc%restraint_coord_pio(3, MaxRest, var_size, var_size1), &
                 stat = alloc_stat)

      enefunc%restraint_atom_pio (    1:MaxRest,1:var_size,1:var_size1) = 0
      enefunc%restraint_force_pio(1:4,1:MaxRest,1:var_size,1:var_size1) = 0.0_wp
      enefunc%restraint_coord_pio(1:3,1:MaxRest,1:var_size,1:var_size1) = 0.0_wp

    case (EneFuncBondCell)

      if (allocated(enefunc%bond_exit)) then
        if (size(enefunc%bond_exit(:)) /= var_size) &
          deallocate(enefunc%bond_exit,             &
                     enefunc%bond_add,              &
                     enefunc%bond_exit_index,       &
                     enefunc%buf_bond_integer,      &
                     enefunc%buf_bond_real,         &
                     enefunc%enm_add,               &
                     enefunc%buf_enm_integer,       &
                     enefunc%buf_enm_real,          &
                     enefunc%angle_exit,            &
                     enefunc%angle_add,             &
                     enefunc%angle_exit_index,      &
                     enefunc%buf_angle_integer,     &
                     enefunc%buf_angle_real,        &
                     enefunc%dihed_exit,            &
                     enefunc%dihed_add,             &
                     enefunc%dihed_exit_index,      &
                     enefunc%buf_dihed_integer,     &
                     enefunc%buf_dihed_real,        &
                     enefunc%rb_dihed_exit,         &
                     enefunc%rb_dihed_add,          &
                     enefunc%rb_dihed_exit_index,   &
                     enefunc%buf_rb_dihed_integer,  &
                     enefunc%buf_rb_dihed_real,     &
                     enefunc%impr_exit,             &
                     enefunc%impr_add,              &
                     enefunc%impr_exit_index,       &
                     enefunc%buf_impr_integer,      &
                     enefunc%buf_impr_real,         &
                     enefunc%cmap_exit,             &
                     enefunc%cmap_add,              &
                     enefunc%cmap_exit_index,       &
                     enefunc%buf_cmap_integer,      &
                     enefunc%contact_exit,          &
                     enefunc%contact_add,           &
                     enefunc%contact_exit_index,    &
                     enefunc%buf_contact_integer,   &
                     enefunc%buf_contact_real,      &
                     enefunc%restraint_exit,        &
                     enefunc%restraint_add,         &
                     enefunc%restraint_exit_index,  &
                     enefunc%buf_restraint_integer, &
                     enefunc%buf_restraint_real,    &
                     enefunc%fitting_exit,          &
                     enefunc%fitting_add,           &
                     enefunc%fitting_exit_index,    &
                     enefunc%buf_fitting_integer,   &
                     enefunc%buf_fitting_real,      &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%bond_exit)) &
        allocate(enefunc%bond_exit            (var_size),               &
                 enefunc%bond_add             (var_size1),              &
                 enefunc%bond_exit_index      (BondMove, var_size),     &
                 enefunc%buf_bond_integer     (4, BondMove, var_size1), &
                 enefunc%buf_bond_real        (2, BondMove, var_size1), &
                 enefunc%enm_add              (var_size1),              &
                 enefunc%buf_enm_integer      (2, MaxBond , var_size1), &
                 enefunc%buf_enm_real         (2, MaxBond , var_size1), &
                 enefunc%angle_exit           (var_size),               &
                 enefunc%angle_add            (var_size1),              &
                 enefunc%angle_exit_index     (AngleMove, var_size),    &
                 enefunc%buf_angle_integer    (7, AngleMove, var_size1),&
                 enefunc%buf_angle_real       (4, AngleMove, var_size1),&
                 enefunc%dihed_exit           (var_size),               &
                 enefunc%dihed_add            (var_size1),              &
                 enefunc%dihed_exit_index     (DiheMove, var_size),     &
                 enefunc%buf_dihed_integer    (9, DiheMove, var_size1), &
                 enefunc%buf_dihed_real       (2, DiheMove, var_size1), &
                 enefunc%rb_dihed_exit        (var_size),               &
                 enefunc%rb_dihed_add         (var_size1),              &
                 enefunc%rb_dihed_exit_index  (DiheMove, var_size),     &
                 enefunc%buf_rb_dihed_integer (7, DiheMove, var_size1), &
                 enefunc%buf_rb_dihed_real    (6, DiheMove, var_size1), &
                 enefunc%impr_exit            (var_size),               &
                 enefunc%impr_add             (var_size1),              &
                 enefunc%impr_exit_index      (ImprMove, var_size),     &
                 enefunc%buf_impr_integer     (8, ImprMove, var_size1), &
                 enefunc%buf_impr_real        (2, ImprMove, var_size1), &
                 enefunc%cmap_exit            (var_size),               &
                 enefunc%cmap_add             (var_size1),              &
                 enefunc%cmap_exit_index      (CmapMove, var_size),     &
                 enefunc%buf_cmap_integer     (15,CmapMove, var_size1), &
                 enefunc%contact_exit         (var_size),               &
                 enefunc%contact_add          (var_size1),              &
                 enefunc%contact_exit_index   (ContactMove, var_size),  &
                 enefunc%buf_contact_integer  (2, ContactMove, var_size1), &
                 enefunc%buf_contact_real     (2, ContactMove, var_size1), &
                 enefunc%restraint_exit       (var_size),               &
                 enefunc%restraint_add        (var_size1),              &
                 enefunc%restraint_exit_index (RestMove, var_size),     &
                 enefunc%buf_restraint_integer(RestMove, var_size1),    &
                 enefunc%buf_restraint_real   (7, RestMove, var_size1), &
                 enefunc%fitting_exit         (var_size),               &
                 enefunc%fitting_add          (var_size1),              &
                 enefunc%fitting_exit_index   (RestMove, var_size),     &
                 enefunc%buf_fitting_integer  (RestMove, var_size1),    &
                 enefunc%buf_fitting_real     (7, RestMove, var_size1), &
                 stat = alloc_stat)

      enefunc%bond_exit            (1:var_size)                    = 0
      enefunc%bond_add             (1:var_size1)                   = 0
      enefunc%bond_exit_index      (1:BondMove, 1:var_size)        = 0
      enefunc%buf_bond_integer     (1:4, 1:BondMove, 1:var_size1)  = 0
      enefunc%buf_bond_real        (1:2, 1:BondMove, 1:var_size1)  = 0.0_wp
      enefunc%enm_add              (1:var_size1)                   = 0
      enefunc%buf_enm_integer      (1:2, 1:MaxBond , 1:var_size1)  = 0
      enefunc%buf_enm_real         (1:2, 1:MaxBond , 1:var_size1)  = 0.0_wp
      enefunc%angle_exit           (1:var_size)                    = 0
      enefunc%angle_add            (1:var_size1)                   = 0
      enefunc%angle_exit_index     (1:AngleMove, 1:var_size)       = 0
      enefunc%buf_angle_integer    (1:7, 1:AngleMove, 1:var_size1) = 0
      enefunc%buf_angle_real       (1:4, 1:AngleMove, 1:var_size1) = 0.0_wp
      enefunc%dihed_exit           (1:var_size)                    = 0
      enefunc%dihed_add            (1:var_size1)                   = 0
      enefunc%dihed_exit_index     (1:DiheMove, 1:var_size)        = 0
      enefunc%buf_dihed_integer    (1:9, 1:DiheMove, 1:var_size1)  = 0
      enefunc%buf_dihed_real       (1:2, 1:DiheMove, 1:var_size1)  = 0.0_wp
      enefunc%rb_dihed_exit        (1:var_size)                    = 0
      enefunc%rb_dihed_add         (1:var_size1)                   = 0
      enefunc%rb_dihed_exit_index  (1:DiheMove, 1:var_size)        = 0
      enefunc%buf_rb_dihed_integer (1:7, 1:DiheMove, 1:var_size1)  = 0
      enefunc%buf_rb_dihed_real    (1:6, 1:DiheMove, 1:var_size1)  = 0.0_wp
      enefunc%impr_exit            (1:var_size)                    = 0
      enefunc%impr_add             (1:var_size1)                   = 0
      enefunc%impr_exit_index      (1:ImprMove, 1:var_size)        = 0
      enefunc%buf_impr_integer     (1:8, 1:ImprMove, 1:var_size1)  = 0
      enefunc%buf_impr_real        (1:2, 1:ImprMove, 1:var_size1)  = 0.0_wp
      enefunc%cmap_exit            (1:var_size)                    = 0
      enefunc%cmap_add             (1:var_size1)                   = 0
      enefunc%cmap_exit_index      (1:CmapMove, 1:var_size)        = 0
      enefunc%buf_cmap_integer     (1:9, 1:CmapMove, 1:var_size1)  = 0
      enefunc%buf_cmap_integer     (1:9, 1:CmapMove, 1:var_size1)  = 0
      enefunc%contact_exit         (1:var_size)                    = 0
      enefunc%contact_add          (1:var_size1)                   = 0
      enefunc%contact_exit_index   (1:ContactMove, 1:var_size)        = 0
      enefunc%buf_contact_integer  (1:2, 1:ContactMove, 1:var_size1)  = 0
      enefunc%buf_contact_real     (1:2, 1:ContactMove, 1:var_size1)  = 0.0_wp
      enefunc%restraint_exit       (1:var_size)                    = 0
      enefunc%restraint_add        (1:var_size1)                   = 0
      enefunc%restraint_exit_index (1:RestMove, 1:var_size)        = 0
      enefunc%buf_restraint_integer(1:RestMove, 1:var_size1)       = 0
      enefunc%buf_restraint_real   (1:7, RestMove, 1:var_size1)    = 0.0_wp
      enefunc%fitting_exit         (1:var_size)                    = 0
      enefunc%fitting_add          (1:var_size1)                   = 0
      enefunc%fitting_exit_index   (1:RestMove, 1:var_size)        = 0
      enefunc%buf_fitting_integer  (1:RestMove, 1:var_size1)       = 0
      enefunc%buf_fitting_real     (1:7, RestMove, 1:var_size1)    = 0.0_wp

    case (EneFuncBondCell_FEP)
      ! Deallocate bond, angle, dihed, rb_dihed, impr, cmap for FEP
      if (allocated(enefunc%buf_bond_integer)) then
        deallocate(enefunc%buf_bond_integer,      &
                   enefunc%buf_angle_integer,     &
                   enefunc%buf_dihed_integer,     &
                   enefunc%buf_rb_dihed_integer,  &
                   enefunc%buf_impr_integer,      &
                   enefunc%buf_cmap_integer,      &
                   stat = dealloc_stat)
      end if

      ! Reallocate bond, angle, dihed, rb_dihed, impr, cmap for FEP
      if (.not. allocated(enefunc%buf_bond_integer)) &
        allocate(enefunc%buf_bond_integer     (5, BondMove, var_size1), &
                 enefunc%buf_angle_integer    (8, AngleMove, var_size1),&
                 enefunc%buf_dihed_integer    (10, DiheMove, var_size1), &
                 enefunc%buf_rb_dihed_integer (8, DiheMove, var_size1), &
                 enefunc%buf_impr_integer     (9, ImprMove, var_size1), &
                 enefunc%buf_cmap_integer     (16,CmapMove, var_size1), &
                 stat = alloc_stat)

      enefunc%buf_bond_integer     (1:5, 1:BondMove, 1:var_size1)  = 0
      enefunc%buf_angle_integer    (1:8, 1:AngleMove, 1:var_size1) = 0
      enefunc%buf_dihed_integer    (1:10, 1:DiheMove, 1:var_size1)  = 0
      enefunc%buf_rb_dihed_integer (1:8, 1:DiheMove, 1:var_size1)  = 0
      enefunc%buf_impr_integer     (1:9, 1:ImprMove, 1:var_size1)  = 0
      enefunc%buf_cmap_integer     (1:16, 1:CmapMove, 1:var_size1)  = 0

    case (EneFuncTableWat)

      if (allocated(enefunc%table%water_list)) then
         if (size(enefunc%table%water_list(1,:)) /= var_size) &
           deallocate(enefunc%table%water_list, stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%table%water_list)) &
        allocate(enefunc%table%water_list(4, var_size), stat = alloc_stat)

      enefunc%table%water_list(1:4, 1:var_size) = 0

    case (EneFuncTableSol)

      if (allocated(enefunc%table%solute_list)) then
         if (size(enefunc%table%solute_list(:)) /= var_size) &
           deallocate(enefunc%table%solute_list,      &
                      enefunc%table%solute_list_inv,  &
                      stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%table%solute_list)) &
        allocate(enefunc%table%solute_list    (var_size),  &
                 enefunc%table%solute_list_inv(var_size1), &
                 stat = alloc_stat)

      enefunc%table%solute_list    (1:var_size)  = 0
      enefunc%table%solute_list_inv(1:var_size1) = 0

    case (EneFuncTblWatDomain)

      if (allocated(enefunc%table%table_ene_WW)) then
        if (size(enefunc%table%table_ene(:)) /= 6*var_size) &
          deallocate(enefunc%table%table_ene,    &
                     enefunc%table%table_grad,   &
                     enefunc%table%table_ecor,   &
                     enefunc%table%table_decor,  &
                     enefunc%table%table_vcor,   &
                     enefunc%table%table_dvcor,  &
                     enefunc%table%table_ene_WW, &
                     enefunc%table%table_de_WW,  &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%table%table_ene_WW)) &
        allocate(enefunc%table%table_ene   (6*var_size),    &
                 enefunc%table%table_grad  (6*var_size),    &
                 enefunc%table%table_ecor  (var_size),      &
                 enefunc%table%table_decor (var_size),      &
                 enefunc%table%table_vcor  (var_size),      &
                 enefunc%table%table_dvcor (var_size),      &
                 enefunc%table%table_ene_WW(6*var_size, 3), &
                 enefunc%table%table_de_WW (2*var_size, 3), &
                 stat = alloc_stat)

      enefunc%table%table_ene   (1:6*var_size)      = 0.0_wp
      enefunc%table%table_grad  (1:6*var_size)      = 0.0_wp
      enefunc%table%table_ecor  (1:var_size)        = 0.0_wp
      enefunc%table%table_decor (1:var_size)        = 0.0_wp
      enefunc%table%table_vcor  (1:var_size)        = 0.0_wp
      enefunc%table%table_dvcor (1:var_size)        = 0.0_wp
      enefunc%table%table_ene_WW(1:6*var_size, 1:3) = 0.0_wp
      enefunc%table%table_de_WW (1:2*var_size, 1:3) = 0.0_wp

    case (EneFuncTableDomain)

      if (allocated(enefunc%table%table_ene)) then
        if (size(enefunc%table%table_ene(:)) /= 6*var_size) &
          deallocate(enefunc%table%table_ene,    &
                     enefunc%table%table_grad,   &
                     enefunc%table%table_ecor,   &
                     enefunc%table%table_decor,  &
                     enefunc%table%table_vcor,  &
                     enefunc%table%table_dvcor, &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%table%table_ene)) &
        allocate(enefunc%table%table_ene   (6*var_size), &
                 enefunc%table%table_grad  (6*var_size), &
                 enefunc%table%table_ecor  (var_size),   &
                 enefunc%table%table_decor (var_size),   &
                 enefunc%table%table_vcor  (var_size),   &
                 enefunc%table%table_dvcor (var_size),   &
                 stat = alloc_stat)

      enefunc%table%table_ene   (1:6*var_size) = 0.0_wp
      enefunc%table%table_grad  (1:6*var_size) = 0.0_wp
      enefunc%table%table_ecor  (1:var_size)   = 0.0_wp
      enefunc%table%table_decor (1:var_size)   = 0.0_wp
      enefunc%table%table_vcor  (1:var_size)   = 0.0_wp
      enefunc%table%table_dvcor (1:var_size)   = 0.0_wp

    case (EneFuncAMBERScale)

      if (allocated(enefunc%dihe_scnb)) then
        if (size(enefunc%dihe_scnb(:)) /= var_size) &
          deallocate(enefunc%dihe_scnb,    &
                     enefunc%dihe_scee,    &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(enefunc%dihe_scnb)) &
        allocate(enefunc%dihe_scnb  (0:var_size), &
                 enefunc%dihe_scee  (0:var_size), &
                 stat = alloc_stat)

      enefunc%dihe_scnb(0:var_size) = 0.0_wp
      enefunc%dihe_scee(0:var_size) = 0.0_wp

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
      enefunc%pc_mode(1:var_size) = 0.0_wip
      enefunc%pc_mode_fit(1:var_size) = 0.0_wip
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
        if (size(enefunc%fit_coord(1,MaxRest,:)) /= var_size) then
          deallocate(enefunc%fitting_atom, &
                     enefunc%fit_coord,    &
                     stat = dealloc_stat)
        end if
      end if

      if (.not. allocated(enefunc%fit_coord)) &
        allocate(enefunc%fitting_atom(1:MaxRest, 1:var_size),       &
                 enefunc%fit_coord(1:3, 1:MaxRest, 1:var_size),       &
                 stat = alloc_stat)

      enefunc%fitting_atom(1:MaxRest, 1:var_size)    = 0
      enefunc%fit_coord(1:3, 1:MaxRest, 1:var_size)  = 0.0_wp

    case (EneFuncRestDomain)

      if (allocated(enefunc%rotated_coord)) then
        if (size(enefunc%rotated_coord(1,MaxRest,:)) /= var_size) then
          deallocate(enefunc%rotated_coord,    &
                     stat = dealloc_stat)
        end if
      end if

      if (.not. allocated(enefunc%rotated_coord)) &
        allocate(enefunc%rotated_coord(1:3, 1:MaxRest, var_size),       &
                 stat = alloc_stat)

      enefunc%rotated_coord(1:3, 1:MaxRest, 1:var_size)  = 0.0_wp

    case (EneFuncGamdDih)

      if (allocated(enefunc%gamd%f_dihe_omp)) then
        if (size(enefunc%gamd%f_dihe_omp(1,:,1,1)) /= var_size) then
          deallocate(enefunc%gamd%f_dihe_omp,    &
                     enefunc%gamd%v_dihe_omp,    &
                     stat = dealloc_stat)
        end if
      end if

      if (.not. allocated(enefunc%gamd%f_dihe_omp)) &
        allocate(enefunc%gamd%f_dihe_omp(1:3,1:var_size,1:var_size1,&
                 1:var_size2), &
                 enefunc%gamd%v_dihe_omp(1:3,1:3,1:var_size2), &
                 stat = alloc_stat)

      enefunc%gamd%f_dihe_omp(1:3,1:var_size,1:var_size1,1:var_size2) = 0.0_wp
      enefunc%gamd%v_dihe_omp(1:3,1:3,1:var_size2) = 0.0_dp

    case (EneFuncGamdRest)

      if (allocated(enefunc%gamd%f_rest_omp)) then
        if (size(enefunc%gamd%f_rest_omp(1,:,1,1)) /= var_size) then
          deallocate(enefunc%gamd%f_rest_omp,    &
                     enefunc%gamd%v_rest_omp,    &
                     stat = dealloc_stat)
        end if
      end if

      if (.not. allocated(enefunc%gamd%f_rest_omp)) &
        allocate(enefunc%gamd%f_rest_omp(1:3,1:var_size,1:var_size1,&
                 1:var_size2), &
                 enefunc%gamd%v_rest_omp(1:3,1:3,1:var_size2), &
                 stat = alloc_stat)

      enefunc%gamd%f_rest_omp(1:3,1:var_size,1:var_size1,1:var_size2) = 0.0_wp
      enefunc%gamd%v_rest_omp(1:3,1:3,1:var_size2) = 0.0_dp

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

    case (EneFuncBase)

      if (allocated(enefunc%num_bond)) then
        deallocate(enefunc%num_bond,        &
                   enefunc%num_angle,       &
                   enefunc%num_dihedral,    &
                   enefunc%num_rb_dihedral, &
                   enefunc%num_improper,    &
                   enefunc%num_cmap,        &
                   enefunc%num_restraint,   &
                   enefunc%nfitting,        &
                   enefunc%num_excl,        &
                   enefunc%num_excl14,      &
                   stat = dealloc_stat)
      end if

    case (EneFuncBase_pio)

      if (allocated(enefunc%num_bond_pio)) then
        deallocate(enefunc%num_bond_pio,        &
                   enefunc%num_angle_pio,       &
                   enefunc%num_dihedral_pio,    &
                   enefunc%num_rb_dihedral_pio, &
                   enefunc%num_improper_pio,    &
                   enefunc%num_cmap_pio,        &
                   enefunc%num_restraint_pio)
      end if

    case (EneFuncBond)

      if (allocated(enefunc%bond_list)) then
        deallocate(enefunc%bond_list,        &
                   enefunc%bond_force_const, &
                   enefunc%bond_dist_min,    &
                   enefunc%bond_kind,        &
                   stat = dealloc_stat)
      end if

    case (EneFuncBond_pio)

      if (allocated(enefunc%bond_list_pio)) then
        deallocate(enefunc%bond_list_pio,        &
                   enefunc%bond_force_const_pio, &
                   enefunc%bond_dist_min_pio,    &
                   enefunc%bond_pbc_pio,         &
                   stat = dealloc_stat)
      end if

    case (EneFuncEnm)

      if (allocated(enefunc%num_enm)) then
        deallocate(enefunc%num_enm,         &
                   enefunc%enm_list,        &
                   enefunc%enm_force_const, &
                   enefunc%enm_dist_min,    &
                   stat = dealloc_stat)
      end if

    case (EneFuncExcl)

      if (allocated(enefunc%excl_list)) then
        deallocate(enefunc%excl_list,        &
                   enefunc%excl14_list,      &
                   stat = dealloc_stat)
      end if

    case (EneFuncAngl)

      if (allocated(enefunc%angle_list)) then
        deallocate(enefunc%angle_list,        &
                   enefunc%angle_force_const, &
                   enefunc%angle_theta_min,   &
                   enefunc%urey_force_const,  &
                   enefunc%urey_rmin,         &
                   enefunc%angle_kind,        &
                   stat = dealloc_stat)
      end if

    case (EneFuncAngl_pio)

      if (allocated(enefunc%angle_list_pio)) then
        deallocate(enefunc%angle_list_pio,        &
                   enefunc%angle_force_const_pio, &
                   enefunc%angle_theta_min_pio,   &
                   enefunc%urey_force_const_pio,  &
                   enefunc%urey_rmin_pio,         &
                   enefunc%angle_pbc_pio,         &
                   stat = dealloc_stat)
      end if

    case (EneFuncDihe)

      if (allocated(enefunc%dihe_list)) then
        deallocate(enefunc%dihe_list,        &
                   enefunc%dihe_force_const, &
                   enefunc%dihe_periodicity, &
                   enefunc%dihe_phase,       &
                   enefunc%dihe_kind,        &
                   stat = dealloc_stat)
      end if

    case (EneFuncDihe_pio)

      if (allocated(enefunc%dihe_list_pio)) then
        deallocate(enefunc%dihe_list_pio,        &
                   enefunc%dihe_force_const_pio, &
                   enefunc%dihe_periodicity_pio, &
                   enefunc%dihe_phase_pio,       &
                   enefunc%dihe_pbc_pio,         &
                   stat = dealloc_stat)
      end if

    case (EneFuncRBDihe)

      if (allocated(enefunc%rb_dihe_list)) then
        deallocate(enefunc%rb_dihe_list, &
                   enefunc%rb_dihe_c,    &
                   stat = dealloc_stat)
      end if

    case (EneFuncRBDihe_pio)

      if (allocated(enefunc%rb_dihe_list_pio)) then
        deallocate(enefunc%rb_dihe_list_pio, &
                   enefunc%rb_dihe_c_pio,    &
                   enefunc%rb_dihe_pbc_pio,  &
                   stat = dealloc_stat)
      end if

    case (EneFuncImpr)

      if (allocated(enefunc%impr_list)) then
        deallocate(enefunc%impr_list,        &
                   enefunc%impr_force_const, &
                   enefunc%impr_periodicity, &
                   enefunc%impr_phase,       &
                   stat = dealloc_stat)
      end if

    case (EneFuncImpr_pio)

      if (allocated(enefunc%impr_list_pio)) then
        deallocate(enefunc%impr_list_pio,        &
                   enefunc%impr_force_const_pio, &
                   enefunc%impr_periodicity_pio, &
                   enefunc%impr_phase_pio,       &
                   enefunc%impr_pbc_pio,         &
                   stat = dealloc_stat)
      end if

    case (EneFuncCmap)

      if (allocated(enefunc%cmap_resolution)) then
        deallocate(enefunc%cmap_resolution, &
                   enefunc%cmap_list, &
                   enefunc%cmap_type, &
                   enefunc%cmap_coef, &
                   stat = dealloc_stat)
      end if

    case (EneFuncCmap_pio)

      if (allocated(enefunc%cmap_list_pio)) then
        deallocate(enefunc%cmap_list_pio,       &
                   enefunc%cmap_resolution_pio, &
                   enefunc%cmap_type_pio,       &
                   enefunc%cmap_coef_pio,       &
                   enefunc%cmap_pbc_pio,        &
                   stat = dealloc_stat)
      end if

    case (EneFuncContact)

      if (allocated(enefunc%contact_list)) then
        deallocate(enefunc%contact_list,   &
                   enefunc%contact_lj12,   &
                   enefunc%contact_lj6,    &
                   stat = dealloc_stat)
      end if

    case (EneFuncFEPBonded)

      if (allocated(enefunc%bond_singleB)) then
        deallocate(enefunc%bond_singleB, &
                   enefunc%angl_singleB, &
                   enefunc%dihe_singleB, &
                   enefunc%rb_dihe_singleB, &
                   enefunc%impr_singleB, &
                   enefunc%cmap_singleB, &
                   stat = dealloc_stat)
      end if

    case (EneFuncNbon)

      if (allocated(enefunc%nonb_atom_cls)) then
        deallocate(enefunc%nonb_atom_cls,    &
                   enefunc%nonb_atom_cls_no, &
                   enefunc%nonb_ljtype,      &
                   enefunc%nb14_lj6,         &
                   enefunc%nb14_lj12,        &
                   enefunc%nonb_lj6,         &
                   enefunc%nonb_lj12,        &
                   enefunc%nonb_lj6_factor,  &
                   stat = dealloc_stat)
      end if

    case (EneFuncNonb)

      if (allocated(enefunc%num_nonb_excl)) then
        deallocate(enefunc%num_nonb_excl,   &
                   enefunc%num_nb14_calc,   &
                   stat = dealloc_stat)
      end if

    case (EneFuncNonbList)

      if (allocated(enefunc%nonb_excl_list)) then
        deallocate(enefunc%nonb_excl_list,  &
                   enefunc%nb14_calc_list,  &
                   enefunc%sc_calc_list,    &
                   enefunc%nb14_qq_scale,   &
                   enefunc%nb14_lj_scale,   &
                   stat = dealloc_stat)
      end if

    case (EneFuncNonbFEPList)

      if (allocated(enefunc%nb14_calc_list_fep)) then
        deallocate(enefunc%nb14_calc_list_fep,  &
                   enefunc%sc_calc_list_fep,    &
                   enefunc%nb14_qq_scale_fep,   &
                   enefunc%nb14_lj_scale_fep,   &
                   stat = dealloc_stat)
      end if

    case (EneFuncRefg)

      if (allocated(enefunc%restraint_numatoms)) then
        deallocate(enefunc%restraint_numatoms, &
                   enefunc%restraint_atomlist, &
                   enefunc%restraint_masscoef, &
                   enefunc%restraint_wcom3,    &
                   enefunc%restraint_wcom4,    &
                   enefunc%restraint_wcom5,    &
                   enefunc%restraint_wtmp,     &
                   stat = dealloc_stat)
      end if

    case (EneFuncReff)

      if (allocated(enefunc%restraint_kind)) then
        deallocate(enefunc%restraint_kind,           &
                   enefunc%restraint_grouplist,      &
                   enefunc%restraint_const,          &
                   enefunc%restraint_ref,            &
                   enefunc%restraint_funcgrp,        &
                   enefunc%restraint_exponent_func,  &
                   enefunc%restraint_exponent_dist,  &
                   enefunc%restraint_mode,           &
                   enefunc%restraint_weight_dist,    &
                   enefunc%restraint_wcom1,          &
                   enefunc%restraint_wcom2,          &
                   enefunc%restraint_wdrt,           &
                   enefunc%restraint_rcom1,          &
                   enefunc%restraint_rcom2,          &
                   enefunc%restraint_rdrt,           &
                   enefunc%restraint_rpath_func,     &
                   stat = dealloc_stat)
       end if

    case (EneFuncRefc)

      if (allocated(enefunc%restraint_refcoord)) then
        deallocate(enefunc%restraint_refcoord, &
                   stat = dealloc_stat)
       end if

    case(EneFuncRefr)

      if (allocated(enefunc%restraint_const_replica)) then
        deallocate(enefunc%restraint_const_replica, &
                   enefunc%restraint_ref_replica, &
                   stat = dealloc_stat)
      end if

    case (EneFuncRest)

      if (allocated(enefunc%restraint_atom)) then
        deallocate(enefunc%restraint_atom,                  &
                   enefunc%restraint_force,                 &
                   enefunc%restraint_coord,                 &
                   enefunc%restraint_bondslist_to_atomlist, &
                   enefunc%restraint_bonds_coord,           &
                   enefunc%restraint_bonds_force,           &
                   stat = dealloc_stat)
      end if

    case (EneFuncBondCell)

      if (allocated(enefunc%bond_exit)) then
        deallocate(enefunc%bond_exit,             &
                   enefunc%bond_add,              &
                   enefunc%bond_exit_index,       &
                   enefunc%buf_bond_integer,      &
                   enefunc%buf_bond_real,         &
                   enefunc%enm_add,               &
                   enefunc%buf_enm_integer,       &
                   enefunc%buf_enm_real,          &
                   enefunc%angle_exit,            &
                   enefunc%angle_add,             &
                   enefunc%angle_exit_index,      &
                   enefunc%buf_angle_integer,     &
                   enefunc%buf_angle_real,        &
                   enefunc%dihed_exit,            &
                   enefunc%dihed_add,             &
                   enefunc%dihed_exit_index,      &
                   enefunc%buf_dihed_integer,     &
                   enefunc%buf_dihed_real,        &
                   enefunc%rb_dihed_exit,         &
                   enefunc%rb_dihed_add,          &
                   enefunc%rb_dihed_exit_index,   &
                   enefunc%buf_rb_dihed_integer,  &
                   enefunc%buf_rb_dihed_real,     &
                   enefunc%impr_exit,             &
                   enefunc%impr_add,              &
                   enefunc%impr_exit_index,       &
                   enefunc%buf_impr_integer,      &
                   enefunc%buf_impr_real,         &
                   enefunc%cmap_exit,             &
                   enefunc%cmap_add,              &
                   enefunc%cmap_exit_index,       &
                   enefunc%buf_cmap_integer,      &
                   enefunc%contact_exit,          & 
                   enefunc%contact_add,           &
                   enefunc%contact_exit_index,    &
                   enefunc%buf_contact_integer,   &
                   enefunc%buf_contact_real,      &
                   enefunc%restraint_exit,        &
                   enefunc%restraint_add,         &
                   enefunc%restraint_exit_index,  &
                   enefunc%buf_restraint_integer, &
                   enefunc%buf_restraint_real,    &
                   enefunc%fitting_exit,          &
                   enefunc%fitting_add,           &
                   enefunc%fitting_exit_index,    &
                   enefunc%buf_fitting_integer,   &
                   enefunc%buf_fitting_real,      &
                   stat = dealloc_stat)
      end if

    case (EneFuncTableWat)

      if (allocated(enefunc%table%water_list)) then
         deallocate(enefunc%table%water_list, &
                    stat = dealloc_stat)
      end if

    case (EneFuncTableSol)

      if (allocated(enefunc%table%solute_list)) then
         deallocate(enefunc%table%solute_list,     &
                    enefunc%table%solute_list_inv, &
                    stat = dealloc_stat)
      end if
    case (EneFuncTblWatDomain)

      if (allocated(enefunc%table%table_ene_WW)) then
        deallocate(enefunc%table%table_ene,    &
                   enefunc%table%table_grad,   &
                   enefunc%table%table_ecor,   &
                   enefunc%table%table_decor,  &
                   enefunc%table%table_ene_WW, &
                   enefunc%table%table_de_WW,  &
                   stat = dealloc_stat)
       end if

    case (EneFuncTableDomain)

      if (allocated(enefunc%table%table_ene)) then
        deallocate(enefunc%table%table_ene,   &
                   enefunc%table%table_grad,  &
                   enefunc%table%table_ecor,  &
                   enefunc%table%table_decor, &
                   stat = dealloc_stat)
       end if

    case (EneFuncAMBERScale)

      if (allocated(enefunc%dihe_scnb)) then
        deallocate(enefunc%dihe_scnb,    &
                   enefunc%dihe_scee,    &
                   stat = dealloc_stat)
      end if

    case(EneFuncMode)

      if (allocated(enefunc%pc_mode)) then
        deallocate(enefunc%pc_mode,        &
                   enefunc%pc_mode_fit,    &
                   enefunc%restraint_g2pc, &
                   stat = dealloc_stat)
      end if

    case (EneFuncFitc)

      if (allocated(enefunc%fit_refcoord)) then
        deallocate(enefunc%fit_refcoord, &
                   stat = dealloc_stat)
       end if

    case (EneFuncFitd)

      if (allocated(enefunc%fit_coord)) then
        deallocate(enefunc%fitting_atom, &
                   enefunc%fit_coord,    &
                   stat = dealloc_stat)
      end if

    case (EneFuncRestDomain)

      if (allocated(enefunc%rotated_coord)) then
          deallocate(enefunc%rotated_coord,    &
                     stat = dealloc_stat)
      endif

    case (EneFuncGamdDih)

      if (allocated(enefunc%gamd%f_dihe_omp)) then
          deallocate(enefunc%gamd%f_dihe_omp,    &
                     enefunc%gamd%v_dihe_omp,    &
                     stat = dealloc_stat)
      endif

    case (EneFuncGamdRest)

      if (allocated(enefunc%gamd%f_rest_omp)) then
          deallocate(enefunc%gamd%f_rest_omp,    &
                     enefunc%gamd%v_rest_omp,    &
                     stat = dealloc_stat)
      endif
 
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
  !! @param[out]   enefunc : structure of potential energy function
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_enefunc_all(enefunc)

    ! formal arguments
    type(s_enefunc),         intent(inout) :: enefunc

    call dealloc_enefunc(enefunc, EneFuncBase)
    call dealloc_enefunc(enefunc, EneFuncBond)
    call dealloc_enefunc(enefunc, EneFuncEnm)
    call dealloc_enefunc(enefunc, EneFuncExcl)
    call dealloc_enefunc(enefunc, EneFuncAngl)
    call dealloc_enefunc(enefunc, EneFuncDihe)
    call dealloc_enefunc(enefunc, EneFuncRBDihe)
    call dealloc_enefunc(enefunc, EneFuncImpr)
    call dealloc_enefunc(enefunc, EneFuncCmap)
    call dealloc_enefunc(enefunc, EneFuncContact)
    call dealloc_enefunc(enefunc, EneFuncNbon)
    call dealloc_enefunc(enefunc, EneFuncNonb)
    call dealloc_enefunc(enefunc, EneFuncNonbList)
    call dealloc_enefunc(enefunc, EneFuncRefg)
    call dealloc_enefunc(enefunc, EneFuncReff)
    call dealloc_enefunc(enefunc, EneFuncRefc)
    call dealloc_enefunc(enefunc, EneFuncRefr)
    call dealloc_enefunc(enefunc, EneFuncRest)
    call dealloc_enefunc(enefunc, EneFuncBondCell)
    call dealloc_enefunc(enefunc, EneFuncTableWat)
    call dealloc_enefunc(enefunc, EneFuncTableSol)
    call dealloc_enefunc(enefunc, EneFuncTblWatDomain)
    call dealloc_enefunc(enefunc, EneFuncTableDomain)
    call dealloc_enefunc(enefunc, EneFuncAMBERScale)
    call dealloc_enefunc(enefunc, EneFuncMode)
    call dealloc_enefunc(enefunc, EneFuncFitc)
    call dealloc_enefunc(enefunc, EneFuncFitd)
    call dealloc_enefunc(enefunc, EneFuncGamdDih)
    call dealloc_enefunc(enefunc, EneFuncGamdRest)
    call dealloc_enefunc(enefunc, EneFuncNonbFEPList)
    call dealloc_enefunc(enefunc, EneFuncFEPBonded)

    return

  end subroutine dealloc_enefunc_all

end module sp_enefunc_str_mod
