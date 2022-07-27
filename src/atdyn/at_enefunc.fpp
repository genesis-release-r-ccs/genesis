!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_enefunc_mod
!> @brief   define potential energy functions
!! @authors Yuji Sugita (YS), Takaharu Mori (TM), Chigusa Kobayashi (CK)
!!          Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_enefunc_mod

  use at_enefunc_gromacs_mod
  use at_enefunc_amber_mod
  use at_enefunc_charmm_mod
  use at_enefunc_go_mod
  use at_energy_mod
  use at_restraints_str_mod
  use at_enefunc_str_mod
  use at_energy_str_mod
  use at_boundary_str_mod
  use at_constraints_str_mod
  use at_boundary_mod
  use molecules_str_mod
  use select_mod
  use select_atoms_mod
  use select_atoms_str_mod
  use fitting_mod
  use fitting_str_mod
  use fileio_table_mod
  use fileio_grotop_mod
  use fileio_prmtop_mod
  use fileio_par_mod
  use fileio_gpr_mod
  use fileio_eef1_mod
  use messages_mod
  use constants_mod
  use mpi_parallel_mod

  implicit none
  private

  ! subroutines
  public  :: define_enefunc
  public  :: setup_fitting_atdyn
  public  :: setup_enefunc_dispcorr
  private :: setup_enefunc_fit_refcoord

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc
  !> @brief        a driver subroutine for defining potential energy functions
  !! @authors      TM, CK
  !! @param[in]    ene_info   : ENERGY section control parameters
  !! @param[in]    boundary   : boundary conditions information
  !! @param[in]    par        : CHARMM PAR information
  !! @param[in]    gpr        : GO model parameter information
  !! @param[in]    prmtop     : AMBER parameter topology information
  !! @param[in]    grotop     : GROMACS parameter topology information
  !! @param[in]    table      : loolup table information
  !! @param[in]    molecule   : molecule information
  !! @param[in]    restraints : restraints information
  !! @param[out]   enefunc    : structure of enefunc
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc(ene_info, boundary, par, gpr, prmtop, grotop, eef1,&
                            table, molecule, restraints, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info 
    type(s_boundary),        intent(in)    :: boundary
    type(s_par),             intent(in)    :: par
    type(s_gpr),             intent(in)    :: gpr
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_grotop),          intent(in)    :: grotop
    type(s_eef1),            intent(in)    :: eef1
    type(s_table),           intent(in)    :: table
    type(s_molecule),        intent(in)    :: molecule
    type(s_restraints),      intent(in)    :: restraints
    type(s_enefunc),         intent(inout) :: enefunc


    call init_enefunc(enefunc)

    enefunc%forcefield       = ene_info%forcefield
    enefunc%output_style     = ene_info%output_style
                             
    enefunc%switchdist       = ene_info%switchdist
    enefunc%cutoffdist       = ene_info%cutoffdist
    enefunc%pairlistdist     = ene_info%pairlistdist
    enefunc%go_electrostatic = ene_info%go_electrostatic
    enefunc%dielec_const     = ene_info%dielec_const
    enefunc%debye            = ene_info%debye
    enefunc%force_switch     = ene_info%vdw_force_switch
    enefunc%vdw_shift        = ene_info%vdw_shift
    enefunc%dispersion_corr  = ene_info%dispersion_corr
    enefunc%nonb_limiter     = ene_info%nonb_limiter
    enefunc%contact_check    = ene_info%contact_check
    enefunc%minimum_contact  = ene_info%minimum_contact

    if (par%num_bonds > 0) then
      if (ene_info%forcefield /= ForcefieldCHARMM .and.       &
          ene_info%forcefield /= ForcefieldCHARMM19 .and.     &
          ene_info%forcefield /= ForcefieldKBGO) &
        call error_msg('Define_Enefunc> ERROR: Bad combination between the input files and force field type (see Chapter "INPUT SECTION" in the user manual).')

      call define_enefunc_charmm &
                         (ene_info, boundary, par, eef1, table, molecule, &
                          restraints, enefunc)

    else if (gpr%num_bonds > 0) then
      if (ene_info%forcefield /= ForcefieldAAGO .and.         &
          ene_info%forcefield /= ForcefieldCAGO .and.         &
          ene_info%forcefield /= ForcefieldKBGO) &
        call error_msg('Define_Enefunc> ERROR: Bad combination between the input files and force field type (see Chapter "INPUT SECTION" in the user manual).')

      call define_enefunc_go &
                         (ene_info, gpr, molecule, restraints, enefunc)

    else if (prmtop%num_atoms > 0) then
      if (ene_info%forcefield /= ForcefieldAMBER) &
        call error_msg('Define_Enefunc> ERROR: Bad combination between the input files and force field type (see Chapter "INPUT SECTION" in the user manual).')

      call define_enefunc_amber &
                          (ene_info, boundary, prmtop, table, molecule, &
                           restraints, enefunc)

    else if (grotop%num_atomtypes > 0) then
      if (ene_info%forcefield /= ForcefieldGROAMBER .and.      &
          ene_info%forcefield /= ForcefieldGROMARTINI .and.    &
          ene_info%forcefield /= ForcefieldAAGO .and.          &
          ene_info%forcefield /= ForcefieldKBGO .and.          &
          ene_info%forcefield /= ForcefieldCAGO .and.          &
          ene_info%forcefield /= ForcefieldSOFT .and.          &
          ene_info%forcefield /= ForcefieldRESIDCG) &
        call error_msg('Define_Enefunc> ERROR: Bad combination between the input files and force field type (see Chapter "INPUT SECTION" in the user manual).')

      call define_enefunc_gromacs &
                          (ene_info, boundary, grotop, table, molecule, &
                           restraints, enefunc)

    end if

    ! spherical potential
    !
    if (boundary%sph_pot) then
      enefunc%spot_use = .true.

    else
      enefunc%spot_use = .false.

    end if

    return

  end subroutine define_enefunc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_fitting_atdyn
  !> @brief        setup fitting in atdyn
  !! @authors      CK
  !! @param[in]    is_rpath : flag if rpath or not
  !! @param[in]    fit_info : fitting information
  !! @param[in]    sel_info : selection information
  !! @param[inout] molecule : molecule information
  !! @param[inout] enefunc  : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_fitting_atdyn(is_rpath, fit_info, sel_info, molecule,enefunc)

    ! formal arguments
    logical,                 intent(in)      :: is_rpath
    type(s_fit_info),        intent(in)      :: fit_info
    type(s_sel_info),        intent(in)      :: sel_info
    type(s_molecule),        intent(in)      :: molecule
    type(s_enefunc), target, intent(inout)   :: enefunc

    ! local variables
    integer                                :: i
    integer                                :: fitting_atom_idx
    integer                                :: fitting_move
    logical                                :: fit_check, fit_check_rpath
    type(s_selatoms)                       :: selatoms


    enefunc%fitting_method  = fit_info%fitting_method

    if (enefunc%pressure_rmsd) then
      if (enefunc%fitting_method /= FittingMethodTR .and.      &
          enefunc%fitting_method /= FittingMethodTR_ROT .and.  &
          enefunc%fitting_method /= FittingMethodTR_ZROT) then
          call error_msg('Setup_Fitting_Atdyn> pressure_rmsd option is '//&
                      'allowed only in Translation fitting')
      end if
    end if

    if (enefunc%fitting_method == FittingMethodNO) then
      if (main_rank) then
        fit_check = .false.
        do i = 1,enefunc%num_restraintfuncs
          if (enefunc%restraint_kind(i) == RestraintsFuncRMSD .or. &
              enefunc%restraint_kind(i) == RestraintsFuncRMSDCOM) then
            fit_check = .true.
          end if
        end do
        if (fit_check .and. .not. fit_info%force_no_fitting) then
          call error_msg('Setup_Fitting_Atdyn> No fit is not allowed '//&
                      'in RMSD restraint')
        end if
        if (fit_check .and. fit_info%force_no_fitting) then
          write(MsgOut,'(A)') "Setup_Fitting_Atdyn> RMSD restraint without FITTING"
          write(MsgOut,*) 
        end if
      end if
      do i = 1,enefunc%num_restraintfuncs
        if (enefunc%restraint_kind(i) == RestraintsFuncPC .or. &
            enefunc%restraint_kind(i) == RestraintsFuncPCCOM) then
           call error_msg('Setup_Fitting_Atdyn> No fit is not allowed '//&
                       'in PC/PCCOM ')
        end if
      end do
      return
    end if

    fit_check       = .false.
    do i = 1,enefunc%num_restraintfuncs
      if (enefunc%restraint_kind(i) == RestraintsFuncRMSD .or. &
          enefunc%restraint_kind(i) == RestraintsFuncRMSDCOM .or. &
          enefunc%restraint_kind(i) == RestraintsFuncPC .or. &
          enefunc%restraint_kind(i) == RestraintsFuncPCCOM) then
          fit_check = .true.
          exit
      end if
    end do

    fit_check_rpath = .false.
    do i = 1,enefunc%num_restraintfuncs
      if (enefunc%restraint_kind(i) == RestraintsFuncPOSI .and. is_rpath) then
        fit_check_rpath = .true.
        exit
      end if
    end do

    if (.not. fit_check .and. .not. fit_check_rpath) then
       if (main_rank) then
         write(MsgOut,'(A)') "Setup_Fitting_Atdyn> NO fitting is applied, skip"
         write(MsgOut,'(A)') "  fitting method  =  NO"
         write(MsgOut,*) 
       end if
       enefunc%fitting_method=FittingMethodNO
       return
    end if

    if (enefunc%fitting_method /= FittingMethodTR_ROT .and. &
        enefunc%fitting_method /= FittingMethodXYTR_ZROT) &
      call error_msg('Setup_Fitting_Atdyn> NO/TR+ROT/XYTR+ZROT is allowed')

    fitting_atom_idx          = fit_info%fitting_atom

    call select_atom(molecule, sel_info%groups(fitting_atom_idx), selatoms)
    enefunc%num_fitting = size(selatoms%idx)
    enefunc%mass_weight = fit_info%mass_weight

    if (is_rpath) then
      enefunc%fitting_file = FittingFileFIT
      enefunc%fitting_move = FittingMoveSYS
    else
      enefunc%fitting_file = FittingFileREF
      enefunc%fitting_move = FittingMoveREF
    end if

    do i = 1,enefunc%num_restraintfuncs
      if (enefunc%restraint_kind(i) == RestraintsFuncRMSDCOM .and. &
        (.not. enefunc%mass_weight) ) &
        call error_msg('Setup_Fitting_Atdyn> RESTRAINTS and FITTING '//&
                       'is inconsistent. Please check [FITTING] section.')
      if (enefunc%restraint_kind(i) == RestraintsFuncPCCOM .and. &
        (.not. enefunc%mass_weight) ) &
        call error_msg('Setup_Fitting_Atdyn> RESTRAINTS and FITTING '//&
                       'is inconsistent.  Please check [FITTING] section.')
      if (enefunc%restraint_kind(i) == RestraintsFuncRMSD .and. &
        (enefunc%mass_weight) ) &
        call error_msg('Setup_Fitting_Atdyn> RESTRAINTS and FITTING'//&
                       'is inconsistent.  Please check [FITTING] section.')
      if (enefunc%restraint_kind(i) == RestraintsFuncPC .and. &
        (enefunc%mass_weight) ) &
        call error_msg('Setup_Fitting_Atdyn> RESTRAINTS and FITTING'//&
                       'is inconsistent.  Please check [FITTING] section.')
      if (enefunc%restraint_kind(i) == RestraintsFuncRMSD .or.    &
          enefunc%restraint_kind(i) == RestraintsFuncRMSDCOM .or. &
          enefunc%restraint_kind(i) == RestraintsFuncPC .or.      &
          enefunc%restraint_kind(i) == RestraintsFuncPCCOM) then
        if (fitting_atom_idx .ne. enefunc%restraint_grouplist(1,i)) then
          call error_msg('Setup_Fitting_Atdyn> Fitting group should be'//&
                       'idenical to the group in restraint')
        end if
      end if
    end do

    call setup_enefunc_fit_refcoord(enefunc%fitting_file, molecule, enefunc)
    enefunc%fitting_atom(1:enefunc%num_fitting) =     &
        selatoms%idx(1:enefunc%num_fitting)

    return

  end subroutine setup_fitting_atdyn

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_dispcorr
  !> @brief        define dispersion correction term 
  !! @authors      CK, JJ
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    constraints : constraints information
  !! @param[inout] enefunc     : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_dispcorr(ene_info, molecule, constraints, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)      :: ene_info
    type(s_molecule),        intent(in)      :: molecule
    type(s_constraints),     intent(in)      :: constraints
    type(s_enefunc), target, intent(inout)   :: enefunc

    ! local variables
    integer                  :: i, j, iatmcls, ntypes
    integer                  :: i1, i2, i3, i4
    integer                  :: natom2, nexpair
    real(wp)                 :: lj6_tot, lj6_diff, lj6_ex
    real(wp)                 :: factor, rpair
    real(wp)                 :: diff_cs, diff_cs2, diff_cs3, diff_cs4
    real(wp)                 :: cutoff , cutoff2, cutoff3, cutoff4
    real(wp)                 :: cutoff5, cutoff6, cutoff7, cutoff8
    real(wp)                 :: cutoff14
    real(wp)                 :: inv_cutoff3, inv_cutoff6, inv_cutoff12
    real(wp)                 :: switchdist , switchdist2, switchdist3
    real(wp)                 :: switchdist4, switchdist5
    real(wp)                 :: switchdist6, switchdist7, switchdist8
    real(wp)                 :: shift_a, shift_b, shift_c
    real(wp)                 :: vswitch, eswitch, vlong

    integer,         pointer :: bondlist(:,:),anglelist(:,:), ureylist(:,:)
    integer,         pointer :: dihelist(:,:),rb_dihelist(:,:)
    integer,         pointer :: imprlist(:,:)

    integer,     allocatable :: atype(:)


    if (ene_info%dispersion_corr == Disp_corr_NONE) &
      return

    bondlist    => enefunc%bond_list
    anglelist   => enefunc%angl_list
    ureylist    => enefunc%urey_list
    dihelist    => enefunc%dihe_list
    rb_dihelist => enefunc%rb_dihe_list
    imprlist    => enefunc%impr_list

    ntypes = enefunc%num_atom_cls
    allocate(atype(1:ntypes))

    atype(1:ntypes) = 0

    do i = 1, molecule%num_atoms
      iatmcls = molecule%atom_cls_no(i)
      atype(iatmcls) = atype(iatmcls)+1
    end do

    lj6_tot = 0.0_wp
    do i = 1, ntypes
      do j = 1, ntypes
        lj6_tot = lj6_tot + enefunc%nonb_lj6(i,j)*atype(i)*atype(j)
      end do
    end do

    deallocate(atype)

    cutoff      = enefunc%cutoffdist
    cutoff2     = cutoff*cutoff
    cutoff3     = cutoff2*cutoff
    inv_cutoff3 = 1.0_wp/cutoff3

    eswitch = 0.0_wp
    vswitch = 0.0_wp
    vlong   = inv_cutoff3/3.0_wp

    if (enefunc%forcefield == ForcefieldAMBER) then

      factor = 2.0_wp*PI*lj6_tot
      enefunc%dispersion_energy = -factor*vlong
      enefunc%dispersion_virial = -2.0_wp*factor*vlong

    else if (enefunc%forcefield == ForcefieldGROAMBER .or.             &
             enefunc%forcefield == ForcefieldGROMARTINI) then

      ! remove exclusion
      !
      lj6_ex = 0.0_wp
      nexpair = 0
      do i = 1, molecule%num_atoms
        iatmcls = molecule%atom_cls_no(i)
        lj6_ex  = lj6_ex + enefunc%nonb_lj6(iatmcls,iatmcls)
      end do

      ! constraints
      !
      do i = 1, constraints%num_bonds
        i1 = constraints%bond_list(1,i)
        i2 = constraints%bond_list(2,i)
        i1 = molecule%atom_cls_no(i1)
        i2 = molecule%atom_cls_no(i2)
        lj6_ex = lj6_ex + 2.0_wp*enefunc%nonb_lj6(i1,i2)
        nexpair = nexpair + 2
      end do

      if (enefunc%settle_func /= 0) then
        do i = 1, constraints%num_water
          i1 = constraints%water_list(1,i)
          i2 = constraints%water_list(2,i)
          i3 = constraints%water_list(3,i)
          i1 = molecule%atom_cls_no(i1)
          i2 = molecule%atom_cls_no(i2)
          i3 = molecule%atom_cls_no(i3)
          lj6_ex = lj6_ex + 2.0_wp*enefunc%nonb_lj6(i1,i2)
          lj6_ex = lj6_ex + 2.0_wp*enefunc%nonb_lj6(i1,i3)
          lj6_ex = lj6_ex + 2.0_wp*enefunc%nonb_lj6(i2,i3)
          nexpair = nexpair + 6
        end do
      end if

      ! bonds
      do i = 1, enefunc%num_bonds
        i1     = molecule%atom_cls_no(bondlist(1,i))
        i2     = molecule%atom_cls_no(bondlist(2,i))
        lj6_ex = lj6_ex + 2.0_wp*enefunc%nb14_lj6(i1,i2)
      end do

      ! angles
      do i = 1, enefunc%num_angles
        i1     = molecule%atom_cls_no(anglelist(1,i))
        i3     = molecule%atom_cls_no(anglelist(3,i))
        lj6_ex = lj6_ex + 2.0_wp*enefunc%nb14_lj6(i1,i3)
      end do

      ! dihedral
      do i = 1, enefunc%num_dihedrals
        i1     = molecule%atom_cls_no(dihelist(1,i))
        i4     = molecule%atom_cls_no(dihelist(4,i))
        lj6_ex = lj6_ex + 2.0_wp*enefunc%nb14_lj6(i1,i4)
      end do

      ! RB dihedral
      do i = 1, enefunc%num_rb_dihedrals
        i1     = molecule%atom_cls_no(rb_dihelist(1,i))
        i4     = molecule%atom_cls_no(rb_dihelist(4,i))
        lj6_ex = lj6_ex + 2.0_wp*enefunc%nb14_lj6(i1,i4)
      end do

      ! improper
      do i = 1, enefunc%num_impropers
        i1     = molecule%atom_cls_no(imprlist(1,i))
        i4     = molecule%atom_cls_no(imprlist(4,i))
        lj6_ex = lj6_ex + 2.0_wp*enefunc%nb14_lj6(i1,i4)
      end do

      nexpair = nexpair + molecule%num_atoms       &
                        + 2*enefunc%num_bonds        &
                        + 2*enefunc%num_angles       &
                        + 2*enefunc%num_dihedrals    &
                        + 2*enefunc%num_rb_dihedrals &
                        + 2*enefunc%num_impropers

      lj6_diff = (lj6_tot - lj6_ex)

      natom2 = molecule%num_atoms*molecule%num_atoms
      rpair  = real(natom2,wp)/real(natom2-nexpair,wp)
      factor = 2.0_wp*PI*rpair*lj6_diff

      switchdist   = enefunc%switchdist
      diff_cs      = (cutoff - switchdist)

      if (diff_cs > EPS) then

        if (enefunc%vdw_shift) then

          cutoff4      = cutoff3*cutoff
          cutoff5      = cutoff4*cutoff
          cutoff6      = cutoff5*cutoff
          cutoff7      = cutoff6*cutoff
          cutoff8      = cutoff7*cutoff
          cutoff14     = cutoff7*cutoff7
          inv_cutoff6  = inv_cutoff3*inv_cutoff3
          inv_cutoff12 = inv_cutoff6*inv_cutoff6
  
          diff_cs2     = diff_cs*diff_cs
          diff_cs3     = diff_cs2*diff_cs
          diff_cs4     = diff_cs3*diff_cs
  
          switchdist2  = switchdist*switchdist
          switchdist3  = switchdist2*switchdist
          switchdist4  = switchdist3*switchdist
          switchdist5  = switchdist4*switchdist
          switchdist6  = switchdist5*switchdist
          switchdist7  = switchdist6*switchdist
          switchdist8  = switchdist7*switchdist
  
          ! LJ6
          ! 
          shift_a = -(10.0_wp*cutoff - 7.0_wp*switchdist)/(cutoff8*diff_cs2)
          shift_b =  ( 9.0_wp*cutoff - 7.0_wp*switchdist)/(cutoff8*diff_cs3)
  
          shift_c = inv_cutoff6 - 2.0_wp * shift_a * diff_cs3  &
                    - 1.5_wp * shift_b * diff_cs4
  
          eswitch = -2.0_wp * shift_a * ((1.0_wp/6.0_wp)*cutoff6               &
                                        -(3.0_wp/5.0_wp)*cutoff5*switchdist    &
                                        +(3.0_wp/4.0_wp)*cutoff4*switchdist2   &
                                        -(1.0_wp/3.0_wp)*cutoff3*switchdist3   &
                                        +(1.0_wp/6.0e1_wp)*switchdist6)        &
                    -1.5_wp * shift_b * ((1.0_wp/7.0_wp)*cutoff7               &
                                        -(2.0_wp/3.0_wp)*cutoff6*switchdist    &
                                        +(6.0_wp/5.0_wp)*cutoff5*switchdist2   &
                                        -                cutoff4*switchdist3   &
                                        +(1.0_wp/3.0_wp)*cutoff3*switchdist4   &
                                        -(1.0_wp/1.05e2_wp)*switchdist7)       &
                    -(1.0_wp/3.0_wp) * shift_c * (cutoff3-switchdist3)
  
          ! LJ12
          !
          shift_a = -(16.0_wp*cutoff - 13.0_wp*switchdist)/(cutoff14*diff_cs2)
          shift_b =  (15.0_wp*cutoff - 13.0_wp*switchdist)/(cutoff14*diff_cs3)
          shift_c = inv_cutoff12 - 2.0_wp * shift_a * diff_cs3  &
                    - 1.5_wp * shift_b * diff_cs4
  
          shift_a = -(10.0_wp*cutoff - 7.0_wp*switchdist)/(cutoff8*diff_cs2)
          shift_b =  ( 9.0_wp*cutoff - 7.0_wp*switchdist)/(cutoff8*diff_cs3)

          vswitch = shift_a * ( (1.0_wp/6.0_wp)*cutoff6                        &
                               -(2.0_wp/5.0_wp)*cutoff5*switchdist             &
                               +(1.0_wp/4.0_wp)*cutoff4*switchdist2            &
                               -(1.0_wp/6.0e1_wp)*switchdist6)                 &
                   +shift_b * ( (1.0_wp/7.0_wp)*cutoff7                        &
                               -(1.0_wp/2.0_wp)*cutoff6*switchdist             &
                               +(3.0_wp/5.0_wp)*cutoff5*switchdist2            &
                               -(1.0_wp/4.0_wp)*cutoff4*switchdist3            &
                               +(1.0_wp/1.4e2_wp)*switchdist7)
          enefunc%dispersion_energy = factor*(eswitch-vlong)
          enefunc%dispersion_virial = -2.0_wp*factor*(-vswitch+vlong)

        else

          eswitch = enefunc%eswitch
          vswitch = enefunc%vswitch
          enefunc%dispersion_energy = factor*(eswitch-vlong)
          enefunc%dispersion_virial = -factor*(vswitch+vlong)

        end if

      else

        enefunc%dispersion_energy = factor*(eswitch-vlong)
        enefunc%dispersion_virial = -2.0_wp*factor*(-vswitch+vlong)

      end if
    end if

  end subroutine setup_enefunc_dispcorr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_fit_refcoord
  !> @brief        setup enefunc fit reference coords
  !! @authors      KT, CK
  !! @param[in]    fitting_file : fitfile or reffile
  !! @param[in]    molecule     : molecule information
  !! @param[inout] enefunc      : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_fit_refcoord(fitting_file, molecule, enefunc)

    ! formal arguments
    integer,                 intent(in)    :: fitting_file
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: j


    if (fitting_file == FittingFileFIT) then
      if (size(molecule%atom_fitcoord(1,:)) /= molecule%num_atoms) &
      call error_msg('Setup_Enefunc_Fit_Refcoord> bad fitfile in [INPUT]')

    else if (fitting_file == FittingFileREF) then
      if (size(enefunc%restraint_refcoord(1,:)) /= molecule%num_atoms) &
      call error_msg('Setup_Enefunc_Fit_Refcoord> bad reffile in [INPUT]')
     
    end if

    enefunc%do_fitting = .true.

    call alloc_enefunc(enefunc, EneFuncFitc, molecule%num_atoms,  &
                       enefunc%num_fitting)

    if (fitting_file == FittingFileFIT) then

      do j = 1, molecule%num_atoms
        enefunc%fit_refcoord(1,j) = molecule%atom_fitcoord(1,j)
        enefunc%fit_refcoord(2,j) = molecule%atom_fitcoord(2,j)
        enefunc%fit_refcoord(3,j) = molecule%atom_fitcoord(3,j)
      end do

    else

      do j = 1, molecule%num_atoms
        enefunc%fit_refcoord(1,j) = enefunc%restraint_refcoord(1,j)
        enefunc%fit_refcoord(2,j) = enefunc%restraint_refcoord(2,j)
        enefunc%fit_refcoord(3,j) = enefunc%restraint_refcoord(3,j)
      end do

    end if

    do j = 1, molecule%num_atoms
      if (enefunc%mass_weight) then
        enefunc%fit_mass(j) = molecule%mass(j)
      else
        enefunc%fit_mass(j) = 1.0_wp
      end if
    end do

    return

  end subroutine setup_enefunc_fit_refcoord

end module at_enefunc_mod
