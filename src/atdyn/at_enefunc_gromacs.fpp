!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_enefunc_gromacs_mod
!> @brief   define potential energy functions
!! @authors Norio Takase (NT), Cheng Tan (CT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_enefunc_gromacs_mod

  use at_enefunc_restraints_mod
  use at_enefunc_table_mod
  use at_enefunc_pme_mod
  use at_energy_mod
  use at_boundary_str_mod
  use at_restraints_str_mod
  use at_enefunc_str_mod
  use at_energy_str_mod
  use nbond_list_mod
  use math_libs_mod
  use molecules_str_mod
  use table_libs_mod
  use fileio_table_mod
  use fileio_grotop_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! constants
  integer,   parameter :: MAXWRK   = 12
  integer,   parameter :: MAX_EXCL = 16 ! = 4(1-2) + 4x3 (1-3)
  integer,   parameter :: MAX_NB14 = 36 ! = 3x4x3 (1-4)

  ! subroutines
  public  :: define_enefunc_gromacs
  private :: setup_enefunc_bond
  private :: setup_enefunc_angl
  private :: setup_enefunc_dihe
  private :: setup_enefunc_rb_dihe
  private :: setup_enefunc_impr
  private :: setup_enefunc_cg_general
  private :: setup_enefunc_cgDNA_base_stack
  private :: setup_enefunc_cgDNA_nonb
  private :: setup_enefunc_cg_ele
  private :: setup_enefunc_cg_KH
  private :: setup_enefunc_cg_PWMcos
  private :: setup_enefunc_cg_PWMcosns
  private :: setup_enefunc_cg_IDR_HPS
  private :: setup_enefunc_cg_IDR_KH
  private :: setup_enefunc_cg_contact
  private :: setup_enefunc_cg_exv
  private :: setup_enefunc_cg_nonlocal_exclusion
  private :: setup_enefunc_nonb
  private :: setup_enefunc_morph
  private :: setup_enefunc_gro_restraints
  private :: setup_enefunc_vsite2
  private :: setup_enefunc_vsite3
  private :: setup_enefunc_vsite3fd
  private :: setup_enefunc_vsite3fad
  private :: setup_enefunc_vsite3out
  private :: setup_enefunc_vsite4fdn
  private :: setup_enefunc_vsiten
  private :: setup_enefunc_ez_membrane
  private :: count_nonb_excl
  private :: count_nonb_excl_go

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc_gromacs
  !> @brief        a driver subroutine for defining potential energy functions
  !! @authors      NT
  !! @param[in]    ene_info   : ENERGY section control parameters information
  !! @param[in]    boundary   : boundary condition
  !! @param[in]    grotop     : GROMACS parameter topology information
  !! @param[in]    table      : lookup table information
  !! @param[in]    molecule   : molecule information
  !! @param[in]    restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc_gromacs(ene_info, boundary, grotop, table, &
                                    molecule, restraints, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_boundary),        intent(in)    :: boundary
    type(s_grotop),          intent(in)    :: grotop
    type(s_table),           intent(in)    :: table
    type(s_molecule),        intent(in)    :: molecule
    type(s_restraints),      intent(in)    :: restraints
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: alloc_stat
    integer                  :: i, num_nonb_excl, num_nb14_calc
    integer                  :: nwork


    ! bond
    !
    call setup_enefunc_bond(grotop, enefunc)

    ! angle
    !
    call setup_enefunc_angl(grotop, enefunc)

    ! dihedral
    !
    call setup_enefunc_dihe(grotop, enefunc)

    ! Ryckaert-Bellemans dihedral
    !
    call setup_enefunc_rb_dihe(grotop, enefunc)

    ! improper
    !
    call setup_enefunc_impr(grotop, enefunc)

    ! other local interactions
    ! 
    if (enefunc%forcefield == ForcefieldRESIDCG) then

      ! General parameters for CG models
      !
      call setup_enefunc_cg_general(ene_info, grotop, molecule, enefunc)

      ! DNA base stacking
      !
      call setup_enefunc_cgDNA_base_stack(ene_info, grotop, molecule, enefunc)

    end if


    ! nonbonded
    !
    if (enefunc%forcefield == ForcefieldRESIDCG) then

      ! DNA base pairing + exv
      !
      call setup_enefunc_cgDNA_nonb(ene_info, grotop, molecule, enefunc)

      ! CG ele: Debye-Huckle
      !
      call setup_enefunc_cg_ele(ene_info, grotop, molecule, enefunc)

      ! CG protein-protein KH model
      !
      call setup_enefunc_cg_KH(ene_info, grotop, molecule, enefunc)

      ! CG PWMcos: protein-DNA interactions
      !
      call setup_enefunc_cg_PWMcos(ene_info, grotop, molecule, enefunc)

      ! CG PWMcosns: protein-DNA interactions
      !
      call setup_enefunc_cg_PWMcosns(ene_info, grotop, molecule, enefunc)

      ! CG IDR: HPS model
      !
      call setup_enefunc_cg_IDR_HPS(ene_info, grotop, molecule, enefunc)

      ! CG IDR: KH model
      !
      call setup_enefunc_cg_IDR_KH(ene_info, grotop, molecule, enefunc)

      ! CG : native contact
      !
      call setup_enefunc_cg_contact(ene_info, grotop, molecule, enefunc)

      ! CG : general exv
      !
      call setup_enefunc_cg_exv(ene_info, grotop, molecule, enefunc)

      ! CG : construct nonlocal excusion list
      !
      call setup_enefunc_cg_nonlocal_exclusion(ene_info, grotop, molecule, enefunc)

    else

      call setup_enefunc_nonb(ene_info, grotop, molecule, enefunc)
      if (enefunc%num_basins > 1 .and. boundary%type /= BoundaryTypeNOBC) then
        call error_msg( &
             'Define_Enefunc_Gromacs> multibasin is allowed only in NOBC')
      endif

    end if


    ! lookup table
    !
    call setup_enefunc_table(ene_info, table, molecule, enefunc)

    ! PME
    !
    call define_enefunc_pme(ene_info, boundary, molecule, enefunc)

    ! restraints
    !
    call setup_enefunc_restraints(molecule, restraints, enefunc)

    ! restraints (gromacs)
    !
    call setup_enefunc_gro_restraints(molecule, grotop, enefunc)

    ! morph
    if (enefunc%restraint_flag .and. restraints%climber_flag) &
      enefunc%morph_flag = .true.
    call setup_enefunc_morph(grotop, enefunc)

    call setup_enefunc_ez_membrane(molecule, grotop, enefunc)

    ! virtual site 2
    !
    call setup_enefunc_vsite2(grotop, enefunc)

    ! virtual site 3
    !
    call setup_enefunc_vsite3(grotop, enefunc)

    ! virtual site 3fd
    !
    call setup_enefunc_vsite3fd(grotop, enefunc)

    ! virtual site 3fad
    !
    call setup_enefunc_vsite3fad(grotop, enefunc)

    ! virtual site 3out
    !
    call setup_enefunc_vsite3out(grotop, enefunc)

    ! virtual site 4fdn
    !
    call setup_enefunc_vsite4fdn(grotop, enefunc)

    ! virtual site n
    !
    call setup_enefunc_vsiten(grotop, enefunc)

    ! allocate working array
    !
    nwork = max(molecule%num_atoms,         &
                enefunc%num_bonds,          &
                enefunc%num_bonds_quartic,  &
                enefunc%num_angles,         &
                enefunc%num_angflex,        &
                enefunc%num_ureys,          &
                enefunc%num_dihedrals,      &
                enefunc%num_dihedflex,      &
                enefunc%num_rb_dihedrals,   &
                enefunc%num_impropers,      &
                enefunc%num_cmaps*2,        &
                enefunc%num_base_stack,     &
                enefunc%num_vsite2,         &
                enefunc%num_vsite3,         &
                enefunc%num_vsite3fd,       &
                enefunc%num_vsite3fad,      &
                enefunc%num_vsite3out,      &
                enefunc%num_vsite4fdn,      &
                enefunc%num_vsiten,         &
                enefunc%num_restraintfuncs, &
                enefunc%num_restraintgroups,&
                enefunc%num_morph_sc,       &
                enefunc%num_morph_bb)

    allocate(enefunc%work(MAXWRK,nwork), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    ! write summary of energy function
    !
    if (main_rank) then

      if (.not. ene_info%table) then

        num_nonb_excl = 0
        num_nb14_calc = 0

        do i = 1, molecule%num_atoms
          num_nonb_excl = num_nonb_excl + enefunc%num_nonb_excl(i)
          num_nb14_calc = num_nb14_calc + enefunc%num_nb14_calc(i)
        end do

      end if

      write(MsgOut,'(A)') &
           'Define_Enefunc_Gromacs> Number of Interactions in Each Term'
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           '  bond_ene        = ', enefunc%num_bonds,           &
           '  angle_ene       = ', enefunc%num_angles
      if (enefunc%num_bonds_quartic > 0)  then
         write(MsgOut,'(A20,I10)')                         &
              '  bond_ene_cgDNA  = ', enefunc%num_bonds_quartic
      end if
      if (enefunc%num_ureys > 0)  then
        write(MsgOut,'(A20,I10)')                         &
           '  urey_ene        = ', enefunc%num_ureys
      end if
      if (enefunc%num_angflex > 0)  then
        write(MsgOut,'(A20,I10)')                         &
           '  flex_angle_ene  = ', enefunc%num_angflex
      end if
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           '  torsion_ene     = ', enefunc%num_dihedrals,       &
           '  rb_torsion_ene  = ', enefunc%num_rb_dihedrals
      if (enefunc%num_dihedflex > 0)  then
        write(MsgOut,'(A20,I10)')                         &
           '  flex_dihed_ene  = ', enefunc%num_dihedflex
      end if
      if (enefunc%num_base_stack > 0)  then
        write(MsgOut,'(A20,I10)')                         &
           '  base_stack_ene  = ', enefunc%num_base_stack
      end if
      write(MsgOut,'(A20,I10)')                                 &
           '  improper_ene    = ', enefunc%num_impropers
      if (.not. ene_info%table) then
        write(MsgOut,'(A20,I10,A20,I10)')                       &
           '  nb_exclusions   = ', num_nonb_excl,               &
           '  nb14_calc       = ', num_nb14_calc
      end if
      if (enefunc%num_contacts > 0) then
        write(MsgOut,'(A20,I10)')                               &
           '  contact_ene     = ', enefunc%num_contacts
      end if
      if (enefunc%num_multi_contacts > 0) then
        write(MsgOut,'(A20,I10)')                                 &
           '  multi_contact   = ', enefunc%num_multi_contacts
      end if
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           '  vsite2_ene      = ', enefunc%num_vsite2,          &
           '  vsite3_ene      = ', enefunc%num_vsite3
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           '  vsite3fd_ene    = ', enefunc%num_vsite3fd,        &
           '  vsite3fad_ene   = ', enefunc%num_vsite3fad
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           '  vsite3out_ene   = ', enefunc%num_vsite3out,       &
           '  vsite4fdn_ene   = ', enefunc%num_vsite4fdn
      write(MsgOut,'(A20,I10)')                                 &
           '  vsiten_ene      = ', enefunc%num_vsiten
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           ' restraint_groups = ', enefunc%num_restraintgroups, &
           ' restraint_funcs  = ', enefunc%num_restraintfuncs
      if (enefunc%morph_flag) then
        write(MsgOut,'(A20,I10,A20,I10)')                         &
           '  morphing_bb     = ', enefunc%num_morph_bb,          &
           '  morphing_sc     = ', enefunc%num_morph_sc
      end if
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine define_enefunc_gromacs

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_bond
  !> @brief        define BOND term in potential energy function
  !! @authors      NT
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond(grotop, enefunc)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: natom, nbond, ioffset
    integer                  :: istart, iend
    ! ~CG~ 3SPN.2C DNA: add local num for quartic bonds
    integer                  :: nbond_quartic
    integer                  :: istart_quartic, iend_quartic

    type(s_grotop_mol), pointer :: gromol


    nbond = 0
    ! ~CG~ 3SPN.2C DNA: 
    nbond_quartic = 0

    enefunc%settle_func = 0
    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      if (gromol%settles%func /= 0) enefunc%settle_func = gromol%settles%func
      do j = 1, grotop%molss(i)%count

        if (gromol%settles%func == 0) then
           ! ~CG~ 3SPN.2C DNA: add to num_bond and num_bond_quartic
           do k = 1, gromol%num_bonds
              if (gromol%bonds(k)%func == 21) then
                 nbond_quartic = nbond_quartic + 1
              else
                 nbond = nbond + 1
              end if
           end do
          ! nbond = nbond + gromol%num_bonds
        else
          nbond = nbond + 3
        end if

      end do
    end do

    call alloc_enefunc(enefunc, EneFuncBond, nbond)
    ! ~CG~ 3SPN.2C DNA: allocate arrays for quartic bonds
    call alloc_enefunc(enefunc, EneFuncBondQuartic, nbond_quartic)

    natom = 0
    nbond = 0
    ! ~CG~ 3SPN.2C DNA: 
    nbond_quartic = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        ioffset = natom
        natom   = natom + gromol%num_atoms

        if (gromol%settles%func == 0) then

          do k = 1, gromol%num_bonds
            ! ~CG~ 3SPN.2C DNA: parse parameters for quartic bonds...
            if (gromol%bonds(k)%func == 21) then
               nbond_quartic = nbond_quartic + 1
               enefunc%bond_quartic_list(1, nbond_quartic) = &
                    gromol%bonds(k)%atom_idx1 + ioffset
               enefunc%bond_quartic_list(2, nbond_quartic) = &
                    gromol%bonds(k)%atom_idx2 + ioffset
               enefunc%bond_quartic_force_const(nbond_quartic) = &
                    gromol%bonds(k)%kb * 0.01_wp * JOU2CAL * 0.5_wp
               enefunc%bond_quartic_dist_min   (nbond_quartic) = &
                    gromol%bonds(k)%b0 * 10.0_wp
            else
               nbond = nbond + 1
               enefunc%bond_list(1, nbond) = gromol%bonds(k)%atom_idx1 + ioffset
               enefunc%bond_list(2, nbond) = gromol%bonds(k)%atom_idx2 + ioffset
               if (gromol%bonds(k)%func == 1) then
                  enefunc%bond_force_const(nbond) = &
                       gromol%bonds(k)%kb * 0.01_wp * JOU2CAL * 0.5_wp
                  enefunc%bond_dist_min   (nbond) = &
                       gromol%bonds(k)%b0 * 10.0_wp
               else ! func == 2
                  enefunc%bond_force_const(nbond) = &
                       gromol%bonds(k)%kb * 0.0001_wp*JOU2CAL * 0.25_wp
                  enefunc%bond_dist_min   (nbond) = &
                       gromol%bonds(k)%b0 * 10.0_wp
               end if
            end if
            ! ~CG~ 3SPN.2C DNA: parse parameters END ---------------------
          end do

        else

          ! settle O-H (1)
          nbond = nbond + 1
          enefunc%bond_list(1, nbond)     = 1 + ioffset
          enefunc%bond_list(2, nbond)     = 2 + ioffset
          enefunc%bond_force_const(nbond) = 0.0_wp
          enefunc%bond_dist_min(nbond)    = gromol%settles%doh * 10.0_wp

          ! settle O-H (2)
          nbond = nbond + 1
          enefunc%bond_list(1, nbond)     = 1 + ioffset
          enefunc%bond_list(2, nbond)     = 3 + ioffset
          enefunc%bond_force_const(nbond) = 0.0_wp
          enefunc%bond_dist_min(nbond)    = gromol%settles%doh * 10.0_wp

          ! settle H-H
          nbond = nbond + 1
          enefunc%bond_list(1, nbond)     = 2 + ioffset
          enefunc%bond_list(2, nbond)     = 3 + ioffset
          enefunc%bond_force_const(nbond) = 0.0_wp
          enefunc%bond_dist_min(nbond)    = gromol%settles%dhh * 10.0_wp

        end if

      end do
    end do


    enefunc%num_bonds = nbond

    call get_loop_index(enefunc%num_bonds, istart, iend)

    enefunc%istart_bond = istart
    enefunc%iend_bond   = iend

    ! ~CG~ 3SPN.2C DNA: add back to num_quartic
    enefunc%num_bonds_quartic = nbond_quartic
    call get_loop_index(enefunc%num_bonds_quartic, istart_quartic, iend_quartic)
    enefunc%istart_bond_quartic = istart_quartic
    enefunc%iend_bond_quartic   = iend_quartic

    return

  end subroutine setup_enefunc_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_angl
  !> @brief        define ANGLE term in potential energy function
  !! @authors      NT
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_angl(grotop, enefunc)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: natom, nangl, nurey, ioffset
    !shinobu-edited
    integer                  :: nanglflex, nanglflextypes, ntable
    integer                  :: istart, iend
    real(wp)                 :: min_th, max_th, min_th_ener, max_th_ener
    real(wp)                 :: min_energy, etmp, center, t123, gradient

    type(s_grotop_mol), pointer :: gromol


    nangl = 0
    nurey = 0
    !shinobu-edited
    nanglflex = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count

        if (gromol%settles%func == 0) then
          do k = 1, gromol%num_angls
            if (gromol%angls(k)%func == 5) then
              nurey = nurey + 1
            !shinobu-edited
            else if (gromol%angls(k)%func == 22) then
              nanglflex = nanglflex + 1
            else
              ! nangl = nangl + gromol%num_angls
              nangl = nangl + 1
            end if
          end do

        else
          nangl = nangl + 1
        end if

      end do
    end do

    call alloc_enefunc(enefunc, EneFuncAngl, nangl)
    call alloc_enefunc(enefunc, EneFuncUrey, nurey)
    !shinobu-edited
    nanglflextypes = grotop%num_flangltypes
    if (nanglflextypes > 0) then

      call alloc_enefunc(enefunc, EneFuncAngFlex, nanglflex)

      ntable = size(grotop%flangltypes(1)%theta(:))
      call alloc_enefunc(enefunc, EneFuncAngFlexTbl, nanglflextypes, ntable)

      center = (AICG2P_FBA_MAX_ANG - AICG2P_FBA_MIN_ANG) * 0.5_wp
      do i = 1, nanglflextypes
        enefunc%anglflex_theta (1:ntable, i) = &
             grotop%flangltypes(i)%theta(1:ntable)
        enefunc%anglflex_efunc (1:ntable, i) = &
             grotop%flangltypes(i)%efunc(1:ntable) * JOU2CAL
        enefunc%anglflex_d2func(1:ntable, i) = &
             grotop%flangltypes(i)%d2func(1:ntable) * JOU2CAL

        t123   = AICG2P_FBA_MIN_ANG
        min_th = AICG2P_FBA_MIN_ANG
        max_th = AICG2P_FBA_MIN_ANG
        call table_flexibleangle(i, t123, enefunc%anglflex_theta,   &
                                enefunc%anglflex_efunc,             &
                                enefunc%anglflex_d2func,            &
                                etmp, gradient)
        min_th_ener = etmp
        min_energy  = etmp

        do while(t123 <= AICG2P_FBA_MAX_ANG)  
          t123 = t123 + AICG2P_FBA_DTHEATA
          call table_flexibleangle(i, t123, enefunc%anglflex_theta, &
                                   enefunc%anglflex_efunc,          &
                                   enefunc%anglflex_d2func,         &
                                   etmp, gradient)
          min_energy = min(min_energy, etmp)
          
          if (gradient < AICG2P_FBA_MIN_ANG_FORCE) then
            min_th      = t123
            min_th_ener =  etmp
          end if
          if (t123 > center .and.                        &
              abs(max_th-AICG2P_FBA_MIN_ANG) < EPS .and. &
              gradient > AICG2P_FBA_MAX_ANG_FORCE) then
            max_th      = t123
            max_th_ener =  etmp
          end if
        end do

        enefunc%anglflex_min_th(1, i) = min_th
        enefunc%anglflex_max_th(1, i) = max_th
        enefunc%anglflex_min_th(2, i) = min_th_ener
        enefunc%anglflex_max_th(2, i) = max_th_ener
        enefunc%anglflex_ener_corr(i) = min_energy
      end do

    else if (nanglflex > 0) then
        call error_msg( &
             'Setup_Enefunc_Angl> Flexible angle type is not defined')
    end if

    natom = 0
    nangl = 0
    nurey = 0
    !shinobu-edited
    nanglflex = 0

    do i = 1, grotop%num_molss

      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        ioffset = natom
        natom   = natom + gromol%num_atoms

        if (gromol%settles%func == 0) then

          do k = 1, gromol%num_angls
             if (gromol%angls(k)%func /= 5 .and. &
                   gromol%angls(k)%func /= 22) then
                nangl = nangl + 1
                enefunc%angl_list(1, nangl) = &
                     gromol%angls(k)%atom_idx1 + ioffset
                enefunc%angl_list(2, nangl) = &
                     gromol%angls(k)%atom_idx2 + ioffset
                enefunc%angl_list(3, nangl) = &
                     gromol%angls(k)%atom_idx3 + ioffset
                enefunc%angl_force_const(nangl) = &
                     gromol%angls(k)%kt * JOU2CAL * 0.5_wp
                enefunc%angl_theta_min  (nangl) = &
                     gromol%angls(k)%theta_0 * RAD
                enefunc%angl_func  (nangl) = &
                     gromol%angls(k)%func
             end if

            if (gromol%angls(k)%func == 10) then
              enefunc%multi_angle = .true.
              enefunc%angl_force_const1(nangl) = &
                   gromol%angls(k)%kt1 * JOU2CAL * 0.5_wp
              enefunc%angl_theta_min1 (nangl) = &
                   gromol%angls(k)%theta_01 * RAD
              enefunc%angl_gamma(nangl) = &
                   gromol%angls(k)%cgamma
              enefunc%angl_epsa(nangl) = &
                   gromol%angls(k)%epsa * JOU2CAL

           else if (gromol%angls(k)%func == 5) then
              nurey = nurey + 1
              enefunc%urey_list(1, nurey) = &
                   gromol%angls(k)%atom_idx1 + ioffset
              enefunc%urey_list(2, nurey) = &
                   gromol%angls(k)%atom_idx3 + ioffset
              enefunc%urey_force_const(nurey) = &
                   gromol%angls(k)%kub * JOU2CAL*0.5_wp*0.01_wp
              enefunc%urey_rmin(nurey) = &
                   gromol%angls(k)%r13 * 10.0_wp

           !shinobu-edited
           else if (gromol%angls(k)%func == 21) then
              enefunc%angl_theta_min  (nangl) = &
                   gromol%angls(k)%theta_0 * 10.0_wp
              enefunc%angl_force_const(nangl) = &
                   -gromol%angls(k)%kt * JOU2CAL
              ! temporary
              ! enefunc%angl_force_const(nangl) = &
              ! gromol%angls(k)%kt
              enefunc%angl_w(nangl) = &
                   gromol%angls(k)%w * 10.0_wp

           !shinobu-edited
           else if (gromol%angls(k)%func == 22) then
             nanglflex = nanglflex + 1
             enefunc%anglflex_list(1, nanglflex) = &
                  gromol%angls(k)%atom_idx1 + ioffset
             enefunc%anglflex_list(2, nanglflex) = &
                  gromol%angls(k)%atom_idx2 + ioffset
             enefunc%anglflex_list(3, nanglflex) = &
                  gromol%angls(k)%atom_idx3 + ioffset
             enefunc%anglflex_type(nanglflex)    = &
                  gromol%angls(k)%types

            end if

          end do

        else

          nangl = nangl + 1
          enefunc%angl_list(1, nangl) = 2 + ioffset
          enefunc%angl_list(2, nangl) = 1 + ioffset
          enefunc%angl_list(3, nangl) = 3 + ioffset
          enefunc%angl_force_const(nangl) = 0.0_wp
          enefunc%angl_theta_min  (nangl) = 0.0_wp

        end if

      end do
    end do

    enefunc%num_angles = nangl
    enefunc%num_ureys  = nurey
    enefunc%num_angflex  = nanglflex

    call get_loop_index(enefunc%num_angles, istart, iend)
    enefunc%istart_angle = istart
    enefunc%iend_angle   = iend

    call get_loop_index(enefunc%num_ureys, istart, iend)
    enefunc%istart_urey  = istart
    enefunc%iend_urey    = iend

    call get_loop_index(enefunc%num_angflex, istart, iend)
    enefunc%istart_angflex  = istart
    enefunc%iend_angflex    = iend

    return

  end subroutine setup_enefunc_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_dihe
  !> @brief        define proper DIHEDRAL term in potential energy function
  !! @authors      NT
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_dihe(grotop, enefunc)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: natom, ndihe, ioffset
    integer                  :: ndiheflextypes, ntable
    !shinobu-edited
    integer                  :: ndiheflex
    integer                  :: istart, iend
    real(wp)                 :: th, cos_dih, sin_dih, min_energy, etmp

    type(s_grotop_mol), pointer :: gromol


    ndihe = 0
    !shinobu-edited
    ndiheflex = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        do k = 1, gromol%num_dihes

        !shinobu-edited
          if (gromol%dihes(k)%func /= 1 .and. gromol%dihes(k)%func /= 4 &
              .and. gromol%dihes(k)%func /= 9 .and. gromol%dihes(k)%func /=21 &
              .and. gromol%dihes(k)%func /= 22 &
              .and. gromol%dihes(k)%func /= 31 &
              .and. gromol%dihes(k)%func /= 32 &
              .and. gromol%dihes(k)%func /= 33 &
              .and. gromol%dihes(k)%func /= 41 &
              .and. gromol%dihes(k)%func /= 43 &
              .and. gromol%dihes(k)%func /= 52 &
          ) &
            cycle
          ! ndihe = ndihe + 1
          !shinobu-edited
          if (gromol%dihes(k)%func == 22 .or. &
              gromol%dihes(k)%func == 52) then
             ndiheflex = ndiheflex + 1
          else
             ndihe = ndihe + 1
          end if

          if (gromol%dihes(k)%func == 32 .or. &
              gromol%dihes(k)%func == 33 .or. &
              gromol%dihes(k)%func == 41 .or. &
              gromol%dihes(k)%func == 43 .or. &
              gromol%dihes(k)%func == 52) then
            enefunc%cg_safe_dihedral_calc   = .true.
          end if

        end do
      end do
    end do

    call alloc_enefunc(enefunc, EneFuncDihe, ndihe)

    ndiheflextypes = grotop%num_fldihetypes

    ! parameter for each flexible dihedral angle type
    !
    if (ndiheflextypes > 0) then
      call alloc_enefunc(enefunc, EneFuncDiheFlex, ndiheflex)
      ntable = size(grotop%fldihetypes(1)%coef(:))
      call alloc_enefunc(enefunc, EneFuncDiheFlexTbl, ndiheflextypes, ntable)
      do i = 1, ndiheflextypes
        enefunc%diheflex_coef(1:ntable,i) = &
             grotop%fldihetypes(i)%coef(1:ntable) *JOU2CAL
        th = -PI
        cos_dih = cos(th)
        sin_dih = sin(th)
        min_energy = enefunc%diheflex_coef(1,i)                      &
                   + enefunc%diheflex_coef(2,i)*cos_dih              &
                   + enefunc%diheflex_coef(3,i)*sin_dih              &
                   + enefunc%diheflex_coef(4,i)                      &
                    *(2.0_wp*cos_dih*cos_dih-1.0_wp)                 &
                   + enefunc%diheflex_coef(5,i)                      &
                    *2.0_wp*cos_dih*sin_dih                          &
                   + enefunc%diheflex_coef(6,i)                      &
                    *(4.0_wp*cos_dih*cos_dih*cos_dih-3.0_wp*cos_dih) &
                   + enefunc%diheflex_coef(7,i)                      &
                    *(-4.0_wp*sin_dih*sin_dih*sin_dih+3.0_wp*sin_dih)
        do while (th <= PI)
          th = th + AICG2P_FBA_DTHEATA
          cos_dih = cos(th)
          sin_dih = sin(th)
          etmp = enefunc%diheflex_coef(1,i) &
                     + enefunc%diheflex_coef(2,i)*cos_dih              &
                     + enefunc%diheflex_coef(3,i)*sin_dih              &
                     + enefunc%diheflex_coef(4,i)                      &
                      *(2.0_wp*cos_dih*cos_dih-1.0_wp)                 &
                     + enefunc%diheflex_coef(5,i)                      &
                      *2.0_wp*cos_dih*sin_dih                          &
                     + enefunc%diheflex_coef(6,i)                      &
                      *(4.0_wp*cos_dih*cos_dih*cos_dih-3.0_wp*cos_dih) &
                     + enefunc%diheflex_coef(7,i)                      &
                      *(-4.0_wp*sin_dih*sin_dih*sin_dih+3.0_wp*sin_dih)
          min_energy = min(min_energy, etmp)
        end do
        enefunc%diheflex_ener_corr(i) = min_energy
      end do
    end if

    natom = 0
    ndihe = 0
    ndiheflex = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        ioffset = natom
        natom   = natom + gromol%num_atoms
        do k = 1, gromol%num_dihes
          !shinobu-edited
          if (gromol%dihes(k)%func /= 1 .and. gromol%dihes(k)%func /= 4 &
              .and. gromol%dihes(k)%func /= 9 .and. gromol%dihes(k)%func /= 21 &
              .and. gromol%dihes(k)%func /= 22 &
              .and. gromol%dihes(k)%func /= 31 &
              .and. gromol%dihes(k)%func /= 32 &
              .and. gromol%dihes(k)%func /= 33 &
              .and. gromol%dihes(k)%func /= 41 &
              .and. gromol%dihes(k)%func /= 43 &
              .and. gromol%dihes(k)%func /= 52 &
              ) &
            cycle
          if (gromol%dihes(k)%func /= 22 .and. gromol%dihes(k)%func /= 52) then
            ndihe = ndihe + 1
            enefunc%dihe_list    (1, ndihe) = gromol%dihes(k)%atom_idx1 + ioffset
            enefunc%dihe_list    (2, ndihe) = gromol%dihes(k)%atom_idx2 + ioffset
            enefunc%dihe_list    (3, ndihe) = gromol%dihes(k)%atom_idx3 + ioffset
            enefunc%dihe_list    (4, ndihe) = gromol%dihes(k)%atom_idx4 + ioffset
            enefunc%dihe_force_const(ndihe) = gromol%dihes(k)%kp * JOU2CAL
            enefunc%dihe_phase      (ndihe) = gromol%dihes(k)%ps * RAD
            enefunc%dihe_periodicity(ndihe) = gromol%dihes(k)%multiplicity
            enefunc%dihe_func       (ndihe) = gromol%dihes(k)%func
          end if
          !shinobu-edited
          if (gromol%dihes(k)%func == 21 .or. &
              gromol%dihes(k)%func == 41 .or. &
              gromol%dihes(k)%func == 43 ) then
            ! temporary
            ! enefunc%dihe_force_const(ndihe) = gromol%dihes(k)%kp
            enefunc%dihe_force_const(ndihe) = -gromol%dihes(k)%kp * JOU2CAL
            enefunc%dihe_theta_min(ndihe) = gromol%dihes(k)%theta_0 * RAD
            enefunc%dihe_w        (ndihe) = gromol%dihes(k)%w 
          end if
          !shinobu-edited
          if (gromol%dihes(k)%func == 22 &
              .or. gromol%dihes(k)%func == 52 &
                  ) then
            ndiheflex = ndiheflex + 1
            enefunc%diheflex_list (1, ndiheflex) = &
                 gromol%dihes(k)%atom_idx1 + ioffset
            enefunc%diheflex_list (2, ndiheflex) = &
                 gromol%dihes(k)%atom_idx2 + ioffset
            enefunc%diheflex_list (3, ndiheflex) = &
                 gromol%dihes(k)%atom_idx3 + ioffset
            enefunc%diheflex_list (4, ndiheflex) = &
                 gromol%dihes(k)%atom_idx4 + ioffset
            enefunc%diheflex_type (   ndiheflex) = &
                 gromol%dihes(k)%types
            enefunc%diheflex_func (   ndiheflex) = &
                 gromol%dihes(k)%func

          end if

        end do
      end do
    end do

    enefunc%num_dihedrals = ndihe
    !shinobu-edited
    enefunc%num_dihedflex = ndiheflex

    call get_loop_index(enefunc%num_dihedrals, istart, iend)
    enefunc%istart_dihedral = istart
    enefunc%iend_dihedral   = iend

    !shinobu-edited
    call get_loop_index(enefunc%num_dihedflex, istart, iend)
    enefunc%istart_dihedflex = istart
    enefunc%iend_dihedflex   = iend

    return

  end subroutine setup_enefunc_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_rb_dihe
  !> @brief        define Ryckaert-Bellemans DIHEDRAL term in potential energy
  !!               function
  !! @authors      NT
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_rb_dihe(grotop, enefunc)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: natom, ndihe, ioffset
    integer                  :: istart, iend

    type(s_grotop_mol), pointer :: gromol


    ndihe = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        do k = 1, gromol%num_dihes
          if (gromol%dihes(k)%func /= 3) &
            cycle
          ndihe = ndihe + 1
        end do
      end do
    end do

    call alloc_enefunc(enefunc, EneFuncRBDihe, ndihe)

    natom = 0
    ndihe = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        ioffset = natom
        natom   = natom + gromol%num_atoms
        do k = 1, gromol%num_dihes
          if (gromol%dihes(k)%func /= 3) &
            cycle
          ndihe = ndihe + 1
          enefunc%rb_dihe_list(1, ndihe) = gromol%dihes(k)%atom_idx1 + ioffset
          enefunc%rb_dihe_list(2, ndihe) = gromol%dihes(k)%atom_idx2 + ioffset
          enefunc%rb_dihe_list(3, ndihe) = gromol%dihes(k)%atom_idx3 + ioffset
          enefunc%rb_dihe_list(4, ndihe) = gromol%dihes(k)%atom_idx4 + ioffset
          enefunc%rb_dihe_c (1:6, ndihe) = gromol%dihes(k)%c(1:6) * JOU2CAL
        end do
      end do
    end do

    enefunc%num_rb_dihedrals = ndihe


    call get_loop_index(enefunc%num_rb_dihedrals, istart, iend)

    enefunc%istart_rb_dihed = istart
    enefunc%iend_rb_dihed   = iend

    return

  end subroutine setup_enefunc_rb_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_impr
  !> @brief        define improper DIHEDRAL term in potential energy function
  !! @authors      NT
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_impr(grotop, enefunc)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: natom, nimpr, ioffset
    integer                  :: istart, iend

    type(s_grotop_mol), pointer :: gromol


    nimpr = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        do k = 1, gromol%num_dihes
          if (gromol%dihes(k)%func /= 2) &
            cycle
          nimpr = nimpr + 1
        end do
      end do
    end do

    call alloc_enefunc(enefunc, EneFuncImpr, nimpr)

    natom = 0
    nimpr = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        ioffset = natom
        natom   = natom + gromol%num_atoms
        do k = 1, gromol%num_dihes
          if (gromol%dihes(k)%func /= 2) &
            cycle
          nimpr = nimpr + 1
          enefunc%impr_list    (1, nimpr) = gromol%dihes(k)%atom_idx1 + ioffset
          enefunc%impr_list    (2, nimpr) = gromol%dihes(k)%atom_idx2 + ioffset
          enefunc%impr_list    (3, nimpr) = gromol%dihes(k)%atom_idx3 + ioffset
          enefunc%impr_list    (4, nimpr) = gromol%dihes(k)%atom_idx4 + ioffset
          enefunc%impr_phase      (nimpr) = gromol%dihes(k)%ps * RAD
          enefunc%impr_force_const(nimpr) = gromol%dihes(k)%kp &
                                            * JOU2CAL * 0.5_wp
          enefunc%impr_periodicity(nimpr) = 0
        end do
      end do
    end do

    enefunc%num_impropers = nimpr


    call get_loop_index(enefunc%num_impropers, istart, iend)

    enefunc%istart_improper = istart
    enefunc%iend_improper   = iend

    return

  end subroutine setup_enefunc_impr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cg_general
  !> @brief        define parameters in CG potential energy function
  !! @authors      CT
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[in]    molecule : molecule including molecular information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine setup_enefunc_cg_general(ene_info, grotop, molecule, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer              :: i, k, l
    integer              :: n_atoms
    integer              :: n_phos
    integer              :: n_base
    integer              :: n_dna
    integer              :: alloc_stat

    integer, allocatable :: atom_cls_2_base_type(:)


    ! ===========================
    ! set cutoff and pairlistdist
    ! ===========================
    !
    enefunc%cg_cutoffdist_ele      = ene_info%cg_cutoffdist_ele
    enefunc%cg_cutoffdist_126      = ene_info%cg_cutoffdist_126
    enefunc%cg_cutoffdist_DNAbp    = ene_info%cg_cutoffdist_DNAbp
    enefunc%cg_pairlistdist_ele    = ene_info%cg_pairlistdist_ele
    enefunc%cg_pairlistdist_126    = ene_info%cg_pairlistdist_126
    enefunc%cg_pairlistdist_PWMcos = ene_info%cg_pairlistdist_PWMcos
    enefunc%cg_pairlistdist_DNAbp  = ene_info%cg_pairlistdist_DNAbp
    enefunc%cg_pairlistdist_exv    = ene_info%cg_pairlistdist_exv
    enefunc%cg_ele_sol_T           = ene_info%cg_sol_temperature
    enefunc%cg_ele_sol_IC          = ene_info%cg_sol_ionic_strength

    ! ==========================
    ! allocate arrays in enefunc
    ! ==========================
    ! 
    n_atoms = molecule%num_atoms
    call alloc_enefunc(enefunc, EneFuncCGGeneral, n_atoms)

    ! ------------
    ! NA_base_type
    ! ------------
    ! 
    allocate(atom_cls_2_base_type(grotop%num_atomtypes), &
        stat = alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    do i = 1, grotop%num_atomtypes
      if (grotop%atomtypes(i)%type_name .eq. 'DA') then
        atom_cls_2_base_type(i) = NABaseTypeDBA
      else if (grotop%atomtypes(i)%type_name .eq. 'DC') then
        atom_cls_2_base_type(i) = NABaseTypeDBC
      else if (grotop%atomtypes(i)%type_name .eq. 'DG') then
        atom_cls_2_base_type(i) = NABaseTypeDBG
      else if (grotop%atomtypes(i)%type_name .eq. 'DT') then
        atom_cls_2_base_type(i) = NABaseTypeDBT
      else if (grotop%atomtypes(i)%type_name .eq. 'DP') then
        atom_cls_2_base_type(i) = NABaseTypeDP
      else if (grotop%atomtypes(i)%type_name .eq. 'DS') then
        atom_cls_2_base_type(i) = NABaseTypeDS
      else if (grotop%atomtypes(i)%type_name .eq. 'RA') then
        atom_cls_2_base_type(i) = NABaseTypeRBA
      else if (grotop%atomtypes(i)%type_name .eq. 'RC') then
        atom_cls_2_base_type(i) = NABaseTypeRBC
      else if (grotop%atomtypes(i)%type_name .eq. 'RG') then
        atom_cls_2_base_type(i) = NABaseTypeRBG
      else if (grotop%atomtypes(i)%type_name .eq. 'RU') then
        atom_cls_2_base_type(i) = NABaseTypeRBU
      else if (grotop%atomtypes(i)%type_name .eq. 'RP') then
        atom_cls_2_base_type(i) = NABaseTypeRP
      else if (grotop%atomtypes(i)%type_name .eq. 'RS') then
        atom_cls_2_base_type(i) = NABaseTypeRS
      else
        atom_cls_2_base_type(i) = NABaseTypeProtein
      end if
    end do

    n_phos = 0
    n_base = 0
    n_dna  = 0
    do i = 1, n_atoms
      l = molecule%atom_cls_no(i)
      k = atom_cls_2_base_type(l)
      if (k <= NABaseTypeDBMAX .or. &
          k == NABaseTypeDP .or. k == NABaseTypeDS) then
        enefunc%cg_DNA_base_pair_calc = .true.
        enefunc%cg_DNA_exv_calc       = .true.

        n_dna = n_dna + 1
        if (k <= NABaseTypeDBMAX) then
          n_base = n_base + 1
        else if (k == NABaseTypeDP) then
          n_phos = n_phos + 1
        end if
      end if

      enefunc%NA_base_type(i) = k
      enefunc%atom_cls(i)     = l

    end do

    ! --------------------------------------------------------
    ! Set up look up tables for pairlist and energy/force calc 
    ! --------------------------------------------------------
    ! 
    enefunc%num_cg_particle_DNA_all  = n_dna
    enefunc%num_cg_particle_DNA_base = n_base
    enefunc%num_cg_particle_DNA_phos = n_phos

    if (allocated(enefunc%cg_particle_DNA_all)) &
        deallocate(enefunc%cg_particle_DNA_all)
    allocate(enefunc%cg_particle_DNA_all(n_dna), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    if (allocated(enefunc%cg_particle_DNA_base)) &
        deallocate(enefunc%cg_particle_DNA_base)
    allocate(enefunc%cg_particle_DNA_base(n_base), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    if (allocated(enefunc%cg_particle_DNA_phos)) &
        deallocate(enefunc%cg_particle_DNA_phos)
    allocate(enefunc%cg_particle_DNA_phos(n_phos), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    n_phos = 0
    n_base = 0
    n_dna  = 0
    do i = 1, n_atoms
      k = enefunc%NA_base_type(i)
      if (k <= NABaseTypeDBMAX .or. &
          k == NABaseTypeDP .or. k == NABaseTypeDS) then

        n_dna = n_dna + 1
        enefunc%cg_particle_DNA_all(n_dna) = i

        if (k <= NABaseTypeDBMAX) then
          n_base = n_base + 1
          enefunc%cg_particle_DNA_base(n_base) = i
        else if (k == NABaseTypeDP) then
          n_phos = n_phos + 1
          enefunc%cg_particle_DNA_phos(n_phos) = i
        end if
      end if
    end do


    ! ---------------------------------
    ! Assign mol chain id to every atom
    ! ---------------------------------
    !
    enefunc%mol_chain_id(1:n_atoms) = molecule%molecule_no(1:n_atoms)

    ! --------------------
    ! ending: deallocation
    ! --------------------
    ! 
    deallocate(atom_cls_2_base_type, stat = alloc_stat)

    return

  end subroutine setup_enefunc_cg_general

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cgDNA_base_stack
  !> @brief        define BASE_STACK term in potential energy function
  !! @authors      CT
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[in]    molecule : molecule including molecular information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cgDNA_base_stack(ene_info, grotop, molecule, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k, l
    integer                  :: natom, nbase_stack, ioffset
    integer                  :: istart, iend
    integer                  :: i_S5, i_B5, i_B3, i_B_last
    character(6)             :: basetype_5, basetype_3 ! base type of 5' and 3' bases

    type(s_grotop_mol), pointer :: gromol


    enefunc%cg_infinite_DNA  = ene_info%cg_infinite_DNA
    enefunc%base_stack_alpha = 3.0_wp
    enefunc%base_stack_K     = 6.0_wp

    ! count number of stacking bases
    natom = 0
    nbase_stack = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        ioffset = natom
        natom   = natom + gromol%num_atoms
        do k = 1, gromol%num_atoms
          if (gromol%atoms(k)%atom_name .eq. "DB") then
            if (.not. enefunc%cg_infinite_DNA) then
              if (k < 4) then
                cycle
              end if
              if (molecule%molecule_no(ioffset + k) /= &
                  molecule%molecule_no(ioffset + k - 3)) then
                cycle
              end if
            end if
            nbase_stack = nbase_stack + 1
          end if                ! atom_name == "DB"
        end do                  ! k
      end do                    ! j
    end do                      ! i

    call alloc_enefunc(enefunc, EneFuncBaseStack, nbase_stack)

    ! set parameters for every stacking base couples...
    natom = 0
    nbase_stack = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        ioffset = natom
        natom   = natom + gromol%num_atoms
        do k = 1, gromol%num_atoms
          if (gromol%atoms(k)%atom_name .eq. "DB") then

            if (.not. enefunc%cg_infinite_DNA) then
              if (k < 4) then
                cycle
              end if
              if (molecule%molecule_no(ioffset + k) /= &
                  molecule%molecule_no(ioffset + k - 3)) then
                cycle
              end if
            end if

            nbase_stack = nbase_stack + 1
            if (k > 3 .and. &
                molecule%molecule_no(ioffset + k) == &
                molecule%molecule_no(ioffset + k - 3)) then
              i_B5 = k - 3
            else
              do l = k + 1, gromol%num_atoms
                if (molecule%molecule_no(ioffset + k) == &
                    molecule%molecule_no(ioffset + l) .and. &
                    gromol%atoms(l)%atom_name .eq. "DB") then
                  i_B_last = l
                end if
              end do
              i_B5 = i_B_last
            end if
            i_S5 = i_B5 - 1
            i_B3 = k
            enefunc%base_stack_list(1, nbase_stack) = i_S5 + ioffset ! S
            enefunc%base_stack_list(2, nbase_stack) = i_B5 + ioffset ! 5' B
            enefunc%base_stack_list(3, nbase_stack) = i_B3 + ioffset ! 3' B
            basetype_5 = gromol%atoms(i_B5)%atom_type
            basetype_3 = gromol%atoms(i_B3)%atom_type
            do l = 1, grotop%num_basestacktypes
              if (grotop%basestacktypes(l)%base_type5 == basetype_5 .and. &
                  grotop%basestacktypes(l)%base_type3 == basetype_3) then
                enefunc%base_stack_func     (nbase_stack) = grotop%basestacktypes(l)%func
                enefunc%base_stack_epsilon  (nbase_stack) = grotop%basestacktypes(l)%epsilon  * JOU2CAL
                enefunc%base_stack_sigma    (nbase_stack) = grotop%basestacktypes(l)%sigma    * 10.0_wp
                enefunc%base_stack_theta_bs (nbase_stack) = grotop%basestacktypes(l)%theta_bs * RAD
              end if
            end do

          end if
        end do
      end do
    end do

    enefunc%num_base_stack = nbase_stack

    call get_loop_index(enefunc%num_base_stack, istart, iend)

    enefunc%istart_base_stack = istart
    enefunc%iend_base_stack   = iend

    return

  end subroutine setup_enefunc_cgDNA_base_stack

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cgDNA_nonb
  !> @brief        define nonbonded term for 3SPN.2C DNA model
  !! @authors      CT
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cgDNA_nonb(ene_info, grotop, molecule, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer              :: i, j, k
    integer              :: n_atoms
    integer              :: itype, jtype
    integer              :: alloc_stat
    logical              :: be_matched
    real(wp)             :: sigma_i, sigma_j

    integer, allocatable :: base_type_2_atom_cls(:)


    ! ----------------
    ! Model Parameters
    ! ----------------
    ! 
    enefunc%base_pair_alpha      = 2.0_wp
    enefunc%base_pair_K          = 12.0_wp
    enefunc%base_cross_alpha     = 4.0_wp
    enefunc%base_cross_K         = 8.0_wp
    enefunc%cgDNA_exv_epsilon    = 1.0_wp * JOU2CAL
    ! enefunc%cg_short_range_cutoffdist    = ene_info%cg_short_range_cutoffdist
    ! enefunc%cg_short_range_pairlistdist  = ene_info%cg_short_range_pairlistdist

    n_atoms = molecule%num_atoms

    ! =================================
    ! MAPPING base-pair <==> atom-class
    ! =================================
    ! count number of base-pair types
    ! and make the local mapping arrays
    ! 
    allocate(base_type_2_atom_cls(NABaseTypeNAMAX), &
        stat = alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    base_type_2_atom_cls(1:NABaseTypeNAMAX) = 0
    ! 
    do i = 1, grotop%num_atomtypes
      if (grotop%atomtypes(i)%type_name .eq. 'DA' .or. &
          grotop%atomtypes(i)%type_name .eq. 'DC' .or. &
          grotop%atomtypes(i)%type_name .eq. 'DG' .or. &
          grotop%atomtypes(i)%type_name .eq. 'DT' &
          ! Add something here if ~CG~ RNA will be used
          ) then
        if (grotop%atomtypes(i)%type_name .eq. 'DA') then
          base_type_2_atom_cls(NABaseTypeDBA) = i
        end if
        if (grotop%atomtypes(i)%type_name .eq. 'DC') then
          base_type_2_atom_cls(NABaseTypeDBC) = i
        end if
        if (grotop%atomtypes(i)%type_name .eq. 'DG') then
          base_type_2_atom_cls(NABaseTypeDBG) = i
        end if
        if (grotop%atomtypes(i)%type_name .eq. 'DT') then
          base_type_2_atom_cls(NABaseTypeDBT) = i
        end if
      else if (grotop%atomtypes(i)%type_name .eq. 'DP') then
        base_type_2_atom_cls(NABaseTypeDP) = i
      else if (grotop%atomtypes(i)%type_name .eq. 'DS') then
        base_type_2_atom_cls(NABaseTypeDS) = i
      end if
    end do

    ! --------------------------------------
    ! transfer params from grotop to enefunc
    ! --------------------------------------
    ! 
    call alloc_enefunc(enefunc, EneFuncBasePair, NABaseTypeBMAX)
    call alloc_enefunc(enefunc, EneFuncCGDNAExv, NABaseTypeNAMAX)

    ! enefunc%base_pair_xxx(:)
    do k = 1, grotop%num_basepairtypes
      do i = 1, NABaseTypeBMAX
        itype = base_type_2_atom_cls(i)
        if (itype == 0 ) &
            cycle

        if (grotop%atomtypes(itype)%type_name == &
            grotop%basepairtypes(k)%base_type_a) then

          enefunc%base_pair_theta_1(i) = &
              grotop%basepairtypes(k)%theta_1 * RAD
          enefunc%base_pair_theta_2(i) = &
              grotop%basepairtypes(k)%theta_2 * RAD
          enefunc%base_pair_theta_3(i) = &
              grotop%basepairtypes(k)%theta_3 * RAD
          enefunc%base_pair_phi_1(i)   = &
              grotop%basepairtypes(k)%phi_1 * RAD
          enefunc%base_pair_sigma(i)   = &
              grotop%basepairtypes(k)%sigma * 10.0_wp
          enefunc%base_pair_epsilon(i) = &
              grotop%basepairtypes(k)%epsilon * JOU2CAL
          exit

        end if
      end do
    end do

    ! enefunc%base_cross_xxx(:, :)
    do k = 1, grotop%num_basecrosstypes
      be_matched = .false.
      do i = 1, NABaseTypeBMAX
        itype = base_type_2_atom_cls(i)
        if (itype == 0) &
            cycle
        do j = 1, NABaseTypeBMAX
          jtype = base_type_2_atom_cls(j)
          if (jtype == 0) &
              cycle
          if (grotop%atomtypes(itype)%type_name ==    &
              grotop%basecrosstypes(k)%base_type_a .and. &
              grotop%atomtypes(jtype)%type_name ==    &
              grotop%basecrosstypes(k)%base_type_b) then

            if (grotop%basecrosstypes(k)%func == 1) then

              enefunc%base_cross_1_epsilon(i, j)  = &
                  grotop%basecrosstypes(k)%epsilon * JOU2CAL
              enefunc%base_cross_1_sigma(i, j)    = &
                  grotop%basecrosstypes(k)%sigma * 10.0_wp
              enefunc%base_cross_1_theta_cs(i, j) = &
                  grotop%basecrosstypes(k)%theta_cs * RAD
            else
              enefunc%base_cross_2_epsilon(i, j)  = &
                  grotop%basecrosstypes(k)%epsilon * JOU2CAL
              enefunc%base_cross_2_sigma(i, j)    = &
                  grotop%basecrosstypes(k)%sigma * 10.0_wp
              enefunc%base_cross_2_theta_cs(i, j) = &
                  grotop%basecrosstypes(k)%theta_cs * RAD

            end if
            be_matched = .true.
            exit
          end if
        end do                ! loop j
        if (be_matched ) exit
      end do                  ! loop i
    end do

    ! --------------------------
    ! Find out the base pairs!!!
    ! --------------------------
    ! 
    enefunc%base_pair_is_WC(1:NABaseTypeBMAX, 1:NABaseTypeBMAX) = .false.
    do i = 1, NABaseTypeBMAX
      itype = base_type_2_atom_cls(i)
      if (itype == 0) &
          cycle
      do j = 1, NABaseTypeBMAX
        jtype = base_type_2_atom_cls(j)
        if (jtype == 0) &
            cycle
        if (grotop%atomtypes(itype)%type_name == "DA" .and. &
            grotop%atomtypes(jtype)%type_name == "DT") then
          enefunc%base_pair_is_WC(i, j) = .true.
        else if (grotop%atomtypes(itype)%type_name == "DC" .and. &
            grotop%atomtypes(jtype)%type_name == "DG") then
          enefunc%base_pair_is_WC(i, j) = .true.
        else if (grotop%atomtypes(itype)%type_name == "DG" .and. &
            grotop%atomtypes(jtype)%type_name == "DC") then
          enefunc%base_pair_is_WC(i, j) = .true.
        else if (grotop%atomtypes(itype)%type_name == "DT" .and. &
            grotop%atomtypes(jtype)%type_name == "DA") then
          enefunc%base_pair_is_WC(i, j) = .true.
        end if
      end do                ! loop j
    end do                  ! loop i
    
    ! enefunc%cgDNA_exv_sigma(:)
    ! 
    do i = 1, NABaseTypeNAMAX
      itype = base_type_2_atom_cls(i)
      if (itype == 0) &
        cycle

      do k = 1, grotop%num_cgdnaexvtypes
        if (grotop%atomtypes(itype)%type_name == &
            grotop%cgdnaexvtypes(k)%base_type) then
          sigma_i = grotop%cgdnaexvtypes(k)%sigma * 10.0_wp
          exit
        end if
      end do

      do j = 1, NABaseTypeNAMAX
        jtype = base_type_2_atom_cls(j)
        if (jtype == 0) &
          cycle

        do k = 1, grotop%num_cgdnaexvtypes
          if (grotop%atomtypes(jtype)%type_name == &
              grotop%cgdnaexvtypes(k)%base_type) then
            sigma_j = grotop%cgdnaexvtypes(k)%sigma * 10.0_wp
            exit
          end if
        end do

        enefunc%cgDNA_exv_sigma(i, j) = 0.5_wp * (sigma_i + sigma_j)

      end do                    ! loop j
    end do                      ! loop i

    deallocate(base_type_2_atom_cls, stat = alloc_stat)

    return

  end subroutine setup_enefunc_cgDNA_nonb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cg_ele
  !> @brief        define electrostatic params
  !! @authors      CT
  !! @param[in]    ene_info : ENERGY section control parameters
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cg_ele(ene_info, grotop, molecule, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer              :: n_atoms
    integer              :: n_charged
    integer              :: n_pro_charged
    integer              :: n_mols
    real(wp)             :: e_T, a_C
    real(wp)             :: sol_T, sol_C
    integer              :: i, j, k
    integer              :: alloc_stat


    ! --------------------
    ! electrostatic params
    ! --------------------
    !
    if (grotop%num_cgelemolpairs > 0 ) then
      enefunc%cg_ele_calc       = .true.
    end if
    ! enefunc%cg_long_range_cutoffdist   = ene_info%cg_long_range_cutoffdist
    ! enefunc%cg_long_range_pairlistdist = ene_info%cg_long_range_pairlistdist

    enefunc%cg_ele_coef = ELEMENT_CHARGE * ELEMENT_CHARGE /   &
        (4.0_wp * PI * ELECTRIC_CONST) * AVOGADRO * JOU2CAL / &
        1e3_wp * 1e10_wp

    sol_T = enefunc%cg_ele_sol_T
    sol_C = enefunc%cg_ele_sol_IC
    e_T   = 2.494e2_wp - 7.88e-1_wp * sol_T &
        + 7.2e-4_wp * sol_T * sol_T
    a_C   = 1.0e0_wp - 2.551e-1_wp * sol_C  &
        + 5.151e-2_wp * sol_C * sol_C       &
        - 6.889e-3_wp * sol_C * sol_C * sol_C
    ! real calculations moved to compute_energy_CG_ele
    enefunc%cg_dielec_const = e_T * a_C

    ! real calculations moved to compute_energy_CG_ele
    enefunc%cg_debye_length = 1.0e10_wp             &
        * sqrt(                                     &
        (CAL2JOU * ELECTRIC_CONST                   &
        * enefunc%cg_dielec_const * KBOLTZ * sol_T) &
        /                                           &
        (2.0_wp * AVOGADRO * AVOGADRO               &
        * ELEMENT_CHARGE * ELEMENT_CHARGE * sol_C)  &
        )

    ! --------------------------------------
    ! Set charges of CG particles in enefunc
    ! --------------------------------------
    ! 
    n_atoms = molecule%num_atoms
    n_mols = molecule%num_molecules
    call alloc_enefunc(enefunc, EneFuncCGele, n_atoms, n_mols)

    enefunc%cg_charge(1:n_atoms)   = molecule%charge(1:n_atoms)
    enefunc%cg_pro_DNA_ele_scale_Q = ene_info%cg_pro_DNA_ele_scale_Q

    ! ---------------------------------
    ! set CG debye-huckel mol-mol pairs
    ! ---------------------------------
    do i = 1, grotop%num_cgelemolpairs

      if (grotop%cg_ele_mol_pairs(i)%is_intermol) then

        do j = grotop%cg_ele_mol_pairs(i)%grp1_start, grotop%cg_ele_mol_pairs(i)%grp1_end
          do k = grotop%cg_ele_mol_pairs(i)%grp2_start, grotop%cg_ele_mol_pairs(i)%grp2_end
            enefunc%cg_ele_mol_pair(j, k) = grotop%cg_ele_mol_pairs(i)%func
            enefunc%cg_ele_mol_pair(k, j) = grotop%cg_ele_mol_pairs(i)%func
          end do
        end do

      else

        do j = grotop%cg_ele_mol_pairs(i)%grp1_start, grotop%cg_ele_mol_pairs(i)%grp1_end
          enefunc%cg_ele_mol_pair(j, j) = grotop%cg_ele_mol_pairs(i)%func
        end do

      end if

    end do

    ! --------------------------------------------------------
    ! Set up look up tables for pairlist and energy/force calc 
    ! --------------------------------------------------------
    n_charged = 0
    n_pro_charged = 0
    do i = 1, n_atoms
      if (abs( enefunc%cg_charge(i) ) > EPS) then
        n_charged = n_charged + 1
        if (enefunc%NA_base_type(i) == NABaseTypeProtein) then
          n_pro_charged = n_pro_charged + 1
        end if
      end if
    end do
    enefunc%num_cg_particle_charged = n_charged
    enefunc%num_cg_particle_pro_charged = n_pro_charged

    if (allocated(enefunc%cg_particle_charged)) &
        deallocate(enefunc%cg_particle_charged)
    allocate(enefunc%cg_particle_charged(n_charged), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    
    if (allocated(enefunc%cg_particle_pro_charged)) &
        deallocate(enefunc%cg_particle_pro_charged)
    allocate(enefunc%cg_particle_pro_charged(n_pro_charged), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    n_charged = 0
    n_pro_charged = 0
    do i = 1, n_atoms
      if (abs( enefunc%cg_charge(i) ) > EPS ) then
        n_charged = n_charged + 1
        enefunc%cg_particle_charged(n_charged) = i
        if (enefunc%NA_base_type(i) == NABaseTypeProtein ) then
          n_pro_charged = n_pro_charged + 1
          enefunc%cg_particle_pro_charged(n_pro_charged) = i
        end if
      end if
    end do

    return

  end subroutine setup_enefunc_cg_ele

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cg_KH
  !> @brief        define params for cg protein KH model
  !! @authors      CT
  !! @param[in]    ene_info : ENERGY section control parameters
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cg_KH(ene_info, grotop, molecule, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer              :: n_atoms
    integer              :: n_atomtypes
    integer              :: n_mols
    integer              :: n_atom_KH
    integer              :: i, j, k
    integer              :: alloc_stat


    if (grotop%num_cgkhmolpairs > 0) then
      enefunc%cg_KH_calc       = .true.
    end if

    ! -----------------------
    ! set CG KH mol-mol pairs
    ! -----------------------
    ! 
    n_mols = molecule%num_molecules
    call alloc_enefunc(enefunc, EneFuncCGKHmol, n_mols)

    do i = 1, grotop%num_cgkhmolpairs

      if (grotop%cg_KH_mol_pairs(i)%is_intermol) then

        do j = grotop%cg_KH_mol_pairs(i)%grp1_start, grotop%cg_KH_mol_pairs(i)%grp1_end
          do k = grotop%cg_KH_mol_pairs(i)%grp2_start, grotop%cg_KH_mol_pairs(i)%grp2_end
            enefunc%cg_KH_mol_pair(j, k) = grotop%cg_KH_mol_pairs(i)%func
            enefunc%cg_KH_mol_pair(k, j) = grotop%cg_KH_mol_pairs(i)%func
          end do
        end do

      else

        do j = grotop%cg_KH_mol_pairs(i)%grp1_start, grotop%cg_KH_mol_pairs(i)%grp1_end
          enefunc%cg_KH_mol_pair(j, j) = grotop%cg_KH_mol_pairs(i)%func
        end do

      end if

    end do

    ! ---------------------------------------------------
    ! allocate sigma for atoms and epsilon for atom-types
    ! ---------------------------------------------------
    ! 
    n_atoms = molecule%num_atoms
    n_atomtypes = grotop%num_atomtypes
    call alloc_enefunc(enefunc, EneFuncCGKH, n_atoms, n_atomtypes)

    ! ------------------------------
    ! set epsilon for each atom type
    ! ------------------------------
    ! 
    do i = 1, grotop%num_atomtypes
      do j = 1, grotop%num_atomtypes
        do k = 1, grotop%num_cg_pair_MJ_eps
          if (grotop%atomtypes(i)%type_name == &
              grotop%cg_pair_MJ_eps(k)%type_name_1 .and. &
              grotop%atomtypes(j)%type_name == &
              grotop%cg_pair_MJ_eps(k)%type_name_2 ) then
            enefunc%cg_KH_epsilon(i, j) = grotop%cg_pair_MJ_eps(k)%epsilon
          end if
        end do                  !k
      end do                    !j
    end do                      !i

    ! ----------------------------------
    ! set sigma values for each particle
    ! ----------------------------------
    ! 
    n_atom_KH = 0
    do i = 1, n_atoms
      do j = 1, grotop%num_cg_KH_atomtypes
        if (grotop%cg_KH_atomtypes(j)%type_name == molecule%atom_cls_name(i)) then
          enefunc%cg_KH_sigma_half(i) = grotop%cg_KH_atomtypes(j)%sigma * 5.0_wp
          enefunc%cg_pro_use_KH(i) = .true.
          n_atom_KH = n_atom_KH + 1
        end if
      end do
    end do
    enefunc%num_cg_particle_KH = n_atom_KH

    if (allocated(enefunc%cg_particle_KH)) &
        deallocate(enefunc%cg_particle_KH)
    allocate(enefunc%cg_particle_KH(n_atom_KH), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    n_atom_KH = 0
    do i = 1, n_atoms
      if (enefunc%cg_pro_use_KH(i) ) then
        n_atom_KH = n_atom_KH + 1
        enefunc%cg_particle_KH(n_atom_KH) = i
      end if
    end do

    return

  end subroutine setup_enefunc_cg_KH

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cg_PWMcos
  !> @brief        define PWMcos params for protein-DNA interactions
  !! @authors      CT
  !! @param[in]    ene_info : ENERGY section control parameters
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cg_PWMcos(ene_info, grotop, molecule, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer              :: n_atoms
    integer              :: n_mols
    integer              :: n_pwmcos
    integer              :: n_pro_resid
    integer              :: ioffset
    integer              :: i, j, k, l
    integer              :: alloc_stat

    type(s_grotop_mol), pointer :: gromol


    ! -----------------------
    ! setup PWMcos parameters
    ! -----------------------
    !
    enefunc%pwmcos_sigma   = ene_info%cg_PWMcos_sigma
    enefunc%pwmcos_phi     = ene_info%cg_PWMcos_phi * RAD

    ! -----------------------------
    ! setup CG PWMcos mol-mol pairs
    ! -----------------------------
    n_mols = molecule%num_molecules
    allocate(enefunc%pwmcos_mol_pair(n_mols, n_mols), stat=alloc_stat)
    enefunc%pwmcos_mol_pair(:, :) = 0

    do i = 1, grotop%num_pwmcosmolpairs

      do j = grotop%pwmcos_mol_pairs(i)%grp1_start, grotop%pwmcos_mol_pairs(i)%grp1_end
        do k = grotop%pwmcos_mol_pairs(i)%grp2_start, grotop%pwmcos_mol_pairs(i)%grp2_end
          enefunc%pwmcos_mol_pair(j, k) = grotop%pwmcos_mol_pairs(i)%func
          enefunc%pwmcos_mol_pair(k, j) = grotop%pwmcos_mol_pairs(i)%func
        end do
      end do

    end do

    ! --------------------------------------------------
    ! setup pwmcos interaction for every protein residue
    ! --------------------------------------------------
    ! count number of pwmcos
    n_pwmcos    = 0
    n_pro_resid = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        l = 0
        do k = 1, gromol%num_pwmcos
          n_pwmcos = n_pwmcos + 1
          if (gromol%pwmcos(k)%protein_idx /= l) then
            n_pro_resid = n_pro_resid + 1
            l = gromol%pwmcos(k)%protein_idx
          end if
        end do
      end do
    end do

    if (n_pwmcos > 0) then
      enefunc%cg_pwmcos_calc = .true.
    end if
    
    call alloc_enefunc(enefunc, EneFuncPWMcos, n_pwmcos)

    if (allocated(enefunc%pwmcos_involved_resid)) &
        deallocate(enefunc%pwmcos_involved_resid)
    allocate(enefunc%pwmcos_involved_resid(n_pro_resid), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(enefunc%pwmcos_involved_spec (n_pro_resid), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    ! set parameters for every pwmcos term...
    n_atoms     = 0
    n_pwmcos    = 0
    n_pro_resid = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count

        ioffset = n_atoms
        n_atoms = n_atoms + gromol%num_atoms

        l = 0

        do k = 1, gromol%num_pwmcos

          n_pwmcos                             = n_pwmcos + 1

          if (gromol%pwmcos(k)%protein_idx /= l ) then
            l = gromol%pwmcos(k)%protein_idx
            n_pro_resid = n_pro_resid + 1
            enefunc%pwmcos_involved_resid(n_pro_resid) = l + ioffset
            enefunc%pwmcos_involved_spec (n_pro_resid) = gromol%pwmcos(k)%func
          end if

          enefunc%pwmcos_protein_id    (n_pwmcos) = gromol%pwmcos(k)%protein_idx + ioffset
          enefunc%pwmcos_r0            (n_pwmcos) = gromol%pwmcos(k)%r0 * 10.0_wp
          enefunc%pwmcos_theta1        (n_pwmcos) = gromol%pwmcos(k)%theta1 * RAD
          enefunc%pwmcos_theta2        (n_pwmcos) = gromol%pwmcos(k)%theta2 * RAD
          enefunc%pwmcos_theta3        (n_pwmcos) = gromol%pwmcos(k)%theta3 * RAD
          enefunc%pwmcos_ene_A         (n_pwmcos) = gromol%pwmcos(k)%ene_A * KBOLTZ * 300.0_wp
          enefunc%pwmcos_ene_C         (n_pwmcos) = gromol%pwmcos(k)%ene_C * KBOLTZ * 300.0_wp
          enefunc%pwmcos_ene_G         (n_pwmcos) = gromol%pwmcos(k)%ene_G * KBOLTZ * 300.0_wp
          enefunc%pwmcos_ene_T         (n_pwmcos) = gromol%pwmcos(k)%ene_T * KBOLTZ * 300.0_wp
          enefunc%pwmcos_gamma         (n_pwmcos) = gromol%pwmcos(k)%gamma
          enefunc%pwmcos_eps           (n_pwmcos) = gromol%pwmcos(k)%eps_shift
          enefunc%pwmcos_specificity   (n_pwmcos) = gromol%pwmcos(k)%func
          enefunc%pwmcos_to_pairlist_id(n_pwmcos) = n_pro_resid

          if (gromol%pwmcos(k)%protein_idx <= 1 ) then
            enefunc%pwmcos_protein_id_N(n_pwmcos) = gromol%pwmcos(k)%protein_idx + ioffset
          else
            enefunc%pwmcos_protein_id_N(n_pwmcos) = gromol%pwmcos(k)%protein_idx + ioffset - 1
          end if

          if (gromol%pwmcos(k)%protein_idx >= gromol%num_atoms ) then
            enefunc%pwmcos_protein_id_C(n_pwmcos) = gromol%pwmcos(k)%protein_idx + ioffset
          else
            enefunc%pwmcos_protein_id_C(n_pwmcos) = gromol%pwmcos(k)%protein_idx + ioffset + 1
          end if

        end do
      end do
    end do

    enefunc%num_pwmcos_terms = n_pwmcos
    enefunc%num_pwmcos_resid = n_pro_resid

    ! ~CG~ 3SPN.2C DNA DEBUG!!!!!!!!!!
    ! write (*,*) "n_pwmcos terms: ", enefunc%num_pwmcos_terms, "  n_resid : ", enefunc%num_pwmcos_resid
    ! do i = 1, enefunc%num_pwmcos_resid
    !   write (*,*) i, "    proid: ", enefunc%pwmcos_involved_resid(i)
    ! end do

    return

  end subroutine setup_enefunc_cg_PWMcos

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cg_PWMcosns
  !> @brief        define PWMcosns params for protein-DNA interactions
  !! @authors      CT
  !! @param[in]    ene_info : ENERGY section control parameters
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cg_PWMcosns(ene_info, grotop, molecule, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer              :: n_atoms
    integer              :: n_mols
    integer              :: n_pwmcosns
    integer              :: n_pro_resid
    integer              :: ioffset
    integer              :: i, j, k, l
    integer              :: alloc_stat

    type(s_grotop_mol), pointer :: gromol


    ! -----------------------
    ! setup PWMcosns parameters
    ! -----------------------
    !
    enefunc%pwmcosns_sigma   = ene_info%cg_PWMcosns_sigma
    enefunc%pwmcosns_phi     = ene_info%cg_PWMcosns_phi * RAD

    ! -----------------------------
    ! setup CG PWMcosns mol-mol pairs
    ! -----------------------------
    n_mols = molecule%num_molecules
    allocate(enefunc%pwmcosns_mol_pair(n_mols, n_mols), stat=alloc_stat)
    enefunc%pwmcosns_mol_pair(:, :) = 0

    do i = 1, grotop%num_pwmcosnsmolpairs

      do j = grotop%pwmcosns_mol_pairs(i)%grp1_start, grotop%pwmcosns_mol_pairs(i)%grp1_end
        do k = grotop%pwmcosns_mol_pairs(i)%grp2_start, grotop%pwmcosns_mol_pairs(i)%grp2_end
          enefunc%pwmcosns_mol_pair(j, k) = grotop%pwmcosns_mol_pairs(i)%func
          enefunc%pwmcosns_mol_pair(k, j) = grotop%pwmcosns_mol_pairs(i)%func
        end do
      end do

    end do

    ! --------------------------------------------------
    ! setup pwmcosns interaction for every protein residue
    ! --------------------------------------------------
    ! count number of pwmcosns
    n_pwmcosns  = 0
    n_pro_resid = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        l = 0
        do k = 1, gromol%num_pwmcosns
          n_pwmcosns = n_pwmcosns + 1
          if (gromol%pwmcosns(k)%protein_idx /= l) then
            n_pro_resid = n_pro_resid + 1
            l = gromol%pwmcosns(k)%protein_idx
          end if
        end do
      end do
    end do

    if (n_pwmcosns > 0) then
      enefunc%cg_pwmcosns_calc = .true.
    end if
    
    call alloc_enefunc(enefunc, EneFuncPWMcosns, n_pwmcosns)

    if (allocated(enefunc%pwmcosns_involved_resid)) &
        deallocate(enefunc%pwmcosns_involved_resid)
    allocate(enefunc%pwmcosns_involved_resid(n_pro_resid), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(enefunc%pwmcosns_involved_spec (n_pro_resid), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    ! set parameters for every pwmcosns term...
    n_atoms     = 0
    n_pwmcosns  = 0
    n_pro_resid = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count

        ioffset = n_atoms
        n_atoms = n_atoms + gromol%num_atoms

        l = 0

        do k = 1, gromol%num_pwmcosns

          n_pwmcosns = n_pwmcosns + 1

          if (gromol%pwmcosns(k)%protein_idx /= l ) then
            l = gromol%pwmcosns(k)%protein_idx
            n_pro_resid = n_pro_resid + 1
            enefunc%pwmcosns_involved_resid(n_pro_resid) = l + ioffset
            enefunc%pwmcosns_involved_spec (n_pro_resid) = gromol%pwmcosns(k)%func
          end if

          enefunc%pwmcosns_protein_id    (n_pwmcosns) = gromol%pwmcosns(k)%protein_idx + ioffset
          enefunc%pwmcosns_r0            (n_pwmcosns) = gromol%pwmcosns(k)%r0 * 10.0_wp
          enefunc%pwmcosns_theta1        (n_pwmcosns) = gromol%pwmcosns(k)%theta1 * RAD
          enefunc%pwmcosns_theta2        (n_pwmcosns) = gromol%pwmcosns(k)%theta2 * RAD
          enefunc%pwmcosns_ene           (n_pwmcosns) = gromol%pwmcosns(k)%ene
          enefunc%pwmcosns_specificity   (n_pwmcosns) = gromol%pwmcosns(k)%func
          enefunc%pwmcosns_to_pairlist_id(n_pwmcosns) = n_pro_resid

          if (gromol%pwmcosns(k)%protein_idx <= 1 ) then
            enefunc%pwmcosns_protein_id_N(n_pwmcosns) = gromol%pwmcosns(k)%protein_idx + ioffset
          else
            enefunc%pwmcosns_protein_id_N(n_pwmcosns) = gromol%pwmcosns(k)%protein_idx + ioffset - 1
          end if

          if (gromol%pwmcosns(k)%protein_idx >= gromol%num_atoms ) then
            enefunc%pwmcosns_protein_id_C(n_pwmcosns) = gromol%pwmcosns(k)%protein_idx + ioffset
          else
            enefunc%pwmcosns_protein_id_C(n_pwmcosns) = gromol%pwmcosns(k)%protein_idx + ioffset + 1
          end if

        end do
      end do
    end do

    enefunc%num_pwmcosns_terms = n_pwmcosns
    enefunc%num_pwmcosns_resid = n_pro_resid

    ! ~CG~ 3SPN.2C DNA DEBUG!!!!!!!!!!
    ! write (*,*) "n_pwmcosns terms: ", enefunc%num_pwmcosns_terms, "  n_resid : ", enefunc%num_pwmcosns_resid
    ! do i = 1, enefunc%num_pwmcosns_resid
    !   write (*,*) i, "    proid: ", enefunc%pwmcosns_involved_resid(i)
    ! end do

    return

  end subroutine setup_enefunc_cg_PWMcosns
 
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cg_IDR_HPS
  !> @brief        define params for HPS-IDR model
  !! @authors      CT
  !! @param[in]    ene_info : ENERGY section control parameters
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cg_IDR_HPS(ene_info, grotop, molecule, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer              :: n_atoms
    integer              :: n_idr_region
    integer              :: ioffset
    integer              :: i, j, k, l
    integer              :: alloc_stat
    integer              :: n_idr_particles
    integer              :: n_charged
    integer              :: n_pro_charged

    type(s_grotop_mol), pointer :: gromol


    ! ------------------------------------
    ! count IDR HPS "regions" in itp files
    ! ------------------------------------
    !
    n_idr_region = 0
    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        do k = 1, gromol%num_idr_hps
          n_idr_region = n_idr_region + 1
        end do
      end do
    end do

    ! switch on
    if (n_idr_region > 0) then
      enefunc%cg_IDR_HPS_calc = .true.
    end if

    enefunc%cg_IDR_HPS_epsilon = ene_info%cg_IDR_HPS_epsilon

    ! allocate everything
    n_atoms = molecule%num_atoms
    call alloc_enefunc(enefunc, EneFuncCGIDRHPS, n_atoms)

    ! -------------------------------------
    ! Set IDR regions if specified in input
    ! -------------------------------------
    ! 
    n_atoms    = 0
    n_idr_particles = 0
    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count

        ioffset = n_atoms
        n_atoms = n_atoms + gromol%num_atoms

        do k = 1, gromol%num_idr_hps
          do l = gromol%idr_hps(k)%grp_start, gromol%idr_hps(k)%grp_end
            n_idr_particles = n_idr_particles + 1
            enefunc%cg_IDR_HPS_is_IDR(l + ioffset) = .true.
          end do                ! l
        end do                  ! k

      end do                    ! j
    end do                      ! i


    ! ---------------------------------------------
    ! set sigma and lambda values for each particle
    ! ---------------------------------------------
    ! 
    do i = 1, n_atoms
      if (.not. enefunc%cg_IDR_HPS_is_IDR(i) ) then
        cycle
      end if
      do j = 1, grotop%num_cg_IDR_HPS_atomtypes
        if (grotop%cg_IDR_HPS_atomtypes(j)%type_name == molecule%atom_cls_name(i) ) then
          enefunc%cg_IDR_HPS_sigma_half(i)  = grotop%cg_IDR_HPS_atomtypes(j)%sigma * 5.0_wp
          enefunc%cg_IDR_HPS_lambda_half(i) = grotop%cg_IDR_HPS_atomtypes(j)%lambda * 0.5_wp
          enefunc%cg_charge(i)              = grotop%cg_IDR_HPS_atomtypes(j)%charge
        end if
      end do
    end do

    ! --------------------------------------------------------
    ! Set up look up tables for pairlist and energy/force calc 
    ! --------------------------------------------------------
    enefunc%num_cg_particle_IDR_HPS = n_idr_particles

    if (allocated(enefunc%cg_particle_IDR_HPS)) &
        deallocate(enefunc%cg_particle_IDR_HPS)
    allocate(enefunc%cg_particle_IDR_HPS(n_idr_particles), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    n_idr_particles = 0
    do i = 1, n_atoms
      if (enefunc%cg_IDR_HPS_is_IDR(i) ) then
        n_idr_particles = n_idr_particles + 1
        enefunc%cg_particle_IDR_HPS(n_idr_particles) = i
      end if
    end do

    ! --------------------------------------------------------
    ! Set up look up tables for pairlist and energy/force calc 
    ! --------------------------------------------------------
    n_charged = 0
    n_pro_charged = 0
    do i = 1, n_atoms
      if (abs(enefunc%cg_charge(i)) > EPS) then
        n_charged = n_charged + 1
        if (enefunc%NA_base_type(i) == NABaseTypeProtein) then
          n_pro_charged = n_pro_charged + 1
        end if
      end if
    end do
    enefunc%num_cg_particle_charged = n_charged
    enefunc%num_cg_particle_pro_charged = n_pro_charged

    if (allocated(enefunc%cg_particle_charged)) &
        deallocate(enefunc%cg_particle_charged)
    allocate(enefunc%cg_particle_charged(n_charged), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    
    if (allocated(enefunc%cg_particle_pro_charged)) &
        deallocate(enefunc%cg_particle_pro_charged)
    allocate(enefunc%cg_particle_pro_charged(n_pro_charged), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    n_charged = 0
    n_pro_charged = 0
    do i = 1, n_atoms
      if (abs(enefunc%cg_charge(i)) > EPS) then
        n_charged = n_charged + 1
        enefunc%cg_particle_charged(n_charged) = i
        if (enefunc%NA_base_type(i) == NABaseTypeProtein) then
          n_pro_charged = n_pro_charged + 1
          enefunc%cg_particle_pro_charged(n_pro_charged) = i
        end if
      end if
    end do
    ! --------------------------------------------------------
    ! Set up look up tables for pairlist and energy/force calc 
    ! --------------------------------------------------------

    return

  end subroutine setup_enefunc_cg_IDR_HPS

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cg_IDR_KH
  !> @brief        define params for KH-IDR model
  !! @authors      CT
  !! @param[in]    ene_info : ENERGY section control parameters
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cg_IDR_KH(ene_info, grotop, molecule, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer              :: n_atoms
    integer              :: n_atomtypes
    integer              :: n_idr_region
    integer              :: ioffset
    integer              :: i, j, k, l
    integer              :: alloc_stat
    integer              :: n_idr_particles
    integer              :: n_charged
    integer              :: n_pro_charged

    type(s_grotop_mol), pointer :: gromol


    ! ------------------------------------
    ! count IDR KH "regions" in itp files
    ! ------------------------------------
    !
    n_idr_region = 0
    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        do k = 1, gromol%num_idr_kh
          n_idr_region = n_idr_region + 1
        end do
      end do
    end do

    ! switch on
    if (n_idr_region > 0) then
      enefunc%cg_IDR_KH_calc = .true.
    end if

    ! allocate everything
    n_atoms = molecule%num_atoms
    n_atomtypes = grotop%num_atomtypes
    call alloc_enefunc(enefunc, EneFuncCGIDRKH, n_atoms, n_atomtypes)

    ! ------------------------------
    ! set epsilon for each atom type
    ! ------------------------------
    ! 
    do i = 1, grotop%num_atomtypes
      do j = 1, grotop%num_atomtypes
        do k = 1, grotop%num_cg_pair_MJ_eps
          if (grotop%atomtypes(i)%type_name == &
              grotop%cg_pair_MJ_eps(k)%type_name_1 .and. &
              grotop%atomtypes(j)%type_name == &
              grotop%cg_pair_MJ_eps(k)%type_name_2 ) then
            enefunc%cg_IDR_KH_epsilon_D(i, j) = enefunc%cg_KH_mod_D_lambda * &
                ( grotop%cg_pair_MJ_eps(k)%epsilon - enefunc%cg_KH_mod_D_eps_0 ) * 0.593_wp
          end if
        end do                  !k
      end do                    !j
    end do                      !i

    ! -------------------------------------
    ! Set IDR regions if specified in input
    ! -------------------------------------
    ! 
    n_atoms = 0
    n_idr_particles = 0
    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count

        ioffset = n_atoms
        n_atoms = n_atoms + gromol%num_atoms

        do k = 1, gromol%num_idr_kh
          do l = gromol%idr_kh(k)%grp_start, gromol%idr_kh(k)%grp_end
            n_idr_particles = n_idr_particles + 1
            enefunc%cg_IDR_KH_is_IDR(l + ioffset) = .true.
          end do                ! l
        end do                  ! k

      end do                    ! j
    end do                      ! i

    ! ----------------------------------
    ! set sigma values for each particle
    ! ----------------------------------
    ! 
    do i = 1, n_atoms
      if (.not. enefunc%cg_IDR_KH_is_IDR(i)) then
        cycle
      end if
      do j = 1, grotop%num_cg_KH_atomtypes
        if (grotop%cg_KH_atomtypes(j)%type_name .eq. molecule%atom_cls_name(i)) then
          enefunc%cg_IDR_KH_sigma_half(i)  = grotop%cg_KH_atomtypes(j)%sigma * 5.0_wp
          enefunc%cg_charge(i)             = grotop%cg_KH_atomtypes(j)%charge
        end if
      end do
    end do

    ! --------------------------------------------------------
    ! Set up look up tables for pairlist and energy/force calc 
    ! --------------------------------------------------------
    enefunc%num_cg_particle_IDR_KH = n_idr_particles

    if (allocated(enefunc%cg_particle_IDR_KH)) &
        deallocate(enefunc%cg_particle_IDR_KH)
    allocate(enefunc%cg_particle_IDR_KH(n_idr_particles), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    n_idr_particles = 0
    do i = 1, n_atoms
      if (enefunc%cg_IDR_KH_is_IDR(i)) then
        n_idr_particles = n_idr_particles + 1
        enefunc%cg_particle_IDR_KH(n_idr_particles) = i
      end if
    end do

    ! ----------------------------------------------------------------
    ! Set up look up tables for pairlist and energy/force calc (begin)
    ! ----------------------------------------------------------------
    n_charged = 0
    n_pro_charged = 0
    do i = 1, n_atoms
      if (abs( enefunc%cg_charge(i) ) > EPS) then
        n_charged = n_charged + 1
        if (enefunc%NA_base_type(i) == NABaseTypeProtein) then
          n_pro_charged = n_pro_charged + 1
        end if
      end if
    end do
    enefunc%num_cg_particle_charged = n_charged
    enefunc%num_cg_particle_pro_charged = n_pro_charged

    if (allocated(enefunc%cg_particle_charged)) &
        deallocate(enefunc%cg_particle_charged)
    allocate(enefunc%cg_particle_charged(n_charged), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    
    if (allocated(enefunc%cg_particle_pro_charged)) &
        deallocate(enefunc%cg_particle_pro_charged)
    allocate(enefunc%cg_particle_pro_charged(n_pro_charged), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    n_charged = 0
    n_pro_charged = 0
    do i = 1, n_atoms
      if (abs( enefunc%cg_charge(i) ) > EPS) then
        n_charged = n_charged + 1
        enefunc%cg_particle_charged(n_charged) = i
        if (enefunc%NA_base_type(i) == NABaseTypeProtein) then
          n_pro_charged = n_pro_charged + 1
          enefunc%cg_particle_pro_charged(n_pro_charged) = i
        end if
      end if
    end do
    ! --------------------------------------------------------------
    ! Set up look up tables for pairlist and energy/force calc (end)
    ! --------------------------------------------------------------

    return

  end subroutine setup_enefunc_cg_IDR_KH

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_nonb
  !> @brief        define NON-BOND term in potential energy function
  !! @authors      NT
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_nonb(ene_info, grotop, molecule, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: eps, sig, kij, ei, ej, si, sj, ki, kj
    real(wp)                 :: c6i, c6j, c12i, c12j, c6, c10, c12
    real(wp)                 :: vi, vj, wi, wj, vij, wij
    integer                  :: nnonb, i, j, k
    integer                  :: nexcl, natom, npairs, ioffset, excl_level
    integer                  :: istart, iend
    integer                  :: nmcont
    integer                  :: alloc_stat
    logical                  :: ljzero

    integer,        allocatable :: excls(:,:), nrexc(:)
    type(s_grotop_mol), pointer :: gromol


    enefunc%num_atom_cls = grotop%num_atomtypes
    enefunc%nonb_func    = grotop%defaults%nonb_func
    enefunc%fudge_lj     = grotop%defaults%fudge_lj
    enefunc%fudge_qq     = grotop%defaults%fudge_qq

    ELECOEF = ELECOEF_GROMACS

    ! set lennard-jones parameters
    !
    nnonb = enefunc%num_atom_cls

    call alloc_enefunc(enefunc, EneFuncNbon, nnonb)

    do i = 1, nnonb
      do j = 1, nnonb

        ! combination rule
        !

        vij = 0.0_wp
        wij = 0.0_wp
        c12 = 0.0_wp
        c10 = 0.0_wp
        c6  = 0.0_wp

        ljzero = .false.

        if (grotop%num_nbonparms > 0) then

          do k = 1, grotop%num_nbonparms
            if ((grotop%atomtypes(i)%type_name == &
                  grotop%nbonparms(k)%atom_type1 .and. &
                grotop%atomtypes(j)%type_name == &
                  grotop%nbonparms(k)%atom_type2) .or.  &
               (grotop%atomtypes(j)%type_name == &
                  grotop%nbonparms(k)%atom_type1 .and. &
                grotop%atomtypes(i)%type_name == &
                  grotop%nbonparms(k)%atom_type2)) then

              vij = grotop%nbonparms(k)%v
              wij = grotop%nbonparms(k)%w

              exit
            end if
          end do

          if (grotop%defaults%combi_rule == 2) then

            sig = vij * 10.0_wp
            eps = wij * JOU2CAL

            c6  = 4.0_wp * eps * (sig ** 6)
            c10 = 0.0_wp
            c12 = 4.0_wp * eps * (sig ** 12)
            if (eps == 0.0_wp) ljzero = .true.

          else ! combi_rule = 1 or 3

            c6  = vij * 1000000.0_wp * JOU2CAL
            c10 = 0.0_wp
            c12 = wij * 1000000.0_wp * 1000000.0_wp * JOU2CAL
            if (vij == 0.0_wp .and. wij == 0.0_wp) ljzero = .true.

          end if

        end if

        if (c6 == 0.0_wp .and. c10 == 0.0_wp .and. c12 == 0.0_wp) then

          vi = grotop%atomtypes(i)%v
          vj = grotop%atomtypes(j)%v
          wi = grotop%atomtypes(i)%w
          wj = grotop%atomtypes(j)%w

          if (grotop%defaults%combi_rule == 2) then

            si = vi * 10.0_wp
            sj = vj * 10.0_wp

            ei = wi * JOU2CAL
            ej = wj * JOU2CAL

            sig = (si + sj) * 0.5_wp
            eps = sqrt(ei * ej)
            if (eps == 0.0_wp) ljzero = .true.

            if (enefunc%forcefield == ForcefieldRESIDCG) then

              enefunc%nonb_aicg_eps(i,j) = eps
              enefunc%nonb_aicg_sig(i,j) = sig
              c6  = 4.0_wp * eps * (sig ** 6)
              c10 = 0.0_wp
              c12 = 4.0_wp * eps * (sig ** 12)
              if (eps == 0.0_wp) ljzero = .true.

              enefunc%cg_exv_eps_sqrt(i) = sqrt( ei )
              enefunc%cg_exv_eps_sqrt(j) = sqrt( ej )
              enefunc%cg_exv_sig_half(i) = si * 0.5_wp
              enefunc%cg_exv_sig_half(j) = sj * 0.5_wp

            else if (enefunc%forcefield /= ForcefieldSOFT) then

              c6  = 4.0_wp * eps * (sig ** 6)
              c10 = 0.0_wp
              c12 = 4.0_wp * eps * (sig ** 12)
              if (eps == 0.0_wp) ljzero = .true.

            else

              ki = grotop%atomtypes(i)%khh * JOU2CAL / 100.0_wp
              kj = grotop%atomtypes(j)%khh * JOU2CAL / 100.0_wp

              kij = sqrt(ki*kj)

              c6  = kij
              c10 = sig
              c12 = eps

            end if

          else ! combi_rule == 1 or 3

            c6i  = vi * 1.0E6_wp * JOU2CAL
            c6j  = vj * 1.0E6_wp * JOU2CAL

            c12i = wi * 1.0E12_wp * JOU2CAL
            c12j = wj * 1.0E12_wp * JOU2CAL

            c6  = sqrt(c6i  * c6j)
            c10 = 0.0_wp
            c12 = sqrt(c12i * c12j)

            if (c6 == 0.0_wp .and. c12 == 0.0_wp ) ljzero=.true.

            if (enefunc%forcefield == ForcefieldSOFT) then
              call error_msg( &
                   'Setup_Enefunc_Nonb> combination rule shuld be 2 for SOFT.')
            end if

          end if

        end if

        if (main_rank) then
          if (.not. ljzero .and. c6 == 0.0_wp .and. c10 == 0.0_wp  &
              .and. c12 == 0.0_wp) &
            write(MsgOut,'(A,A,A)') &
              'Setup_Enefunc_Nonb> WARNING, combination is not found.',&
               grotop%atomtypes(i)%type_name, grotop%atomtypes(j)%type_name 
        end if
!          call error_msg('Setup_Enefunc_Nonb> combination is not found.')

        ! set parameters
        !
        enefunc%nb14_lj12(i,j) = c12
        enefunc%nb14_lj10(i,j) = c10
        enefunc%nb14_lj6 (i,j) = c6

        enefunc%nonb_lj12(i,j) = c12
        enefunc%nonb_lj10(i,j) = c10
        enefunc%nonb_lj6 (i,j) = c6

      end do
    end do

    if (enefunc%forcefield == ForcefieldRESIDCG) then
      do i = 1, molecule%num_atoms
        enefunc%param_epsilon(i) = enefunc%cg_exv_eps_sqrt(enefunc%atom_cls(i))
        enefunc%param_sigma  (i) = enefunc%cg_exv_sig_half(enefunc%atom_cls(i))
      end do
    end if

    ! create native contact list
    if (enefunc%forcefield == ForcefieldKBGO .or. &
        enefunc%forcefield == ForcefieldAAGO .or. &
        enefunc%forcefield == ForcefieldCAGO .or. &
        enefunc%forcefield == ForcefieldSOFT .or. &
        enefunc%forcefield == ForcefieldRESIDCG) then
      !shinobu-edited

      npairs = 0
      natom  = 0
      do i = 1, grotop%num_molss
        gromol => grotop%molss(i)%moltype%mol
        do j = 1, grotop%molss(i)%count
          npairs = npairs + gromol%num_pairs
          natom  = natom + gromol%num_atoms
        end do
      end do

      enefunc%num_contacts = npairs
      call alloc_enefunc(enefunc, EneFuncNonbGO, natom)
      call alloc_enefunc(enefunc, EneFuncCntc, npairs)

      if (enefunc%forcefield == ForcefieldKBGO .or. &
          enefunc%forcefield == ForcefieldRESIDCG) then
        natom  = 0
        do i = 1, grotop%num_molss
          gromol => grotop%molss(i)%moltype%mol
          do j = 1, grotop%molss(i)%count
            do k = 1, gromol%num_atoms
              natom   = natom + 1
              sig = 10.0_wp * gromol%atoms(k)%v
              eps  = 4.0_wp * JOU2CAL * gromol%atoms(k)%w
              enefunc%nonb_eps(natom)  = eps
              enefunc%nonb_rmin(natom) = sig
            end do
          end do
        end do
      end if
      natom  = 0
      npairs = 0
      do i = 1, grotop%num_molss
        gromol => grotop%molss(i)%moltype%mol
        do j = 1, grotop%molss(i)%count
          ioffset = natom
          natom   = natom + gromol%num_atoms
          do k = 1, gromol%num_pairs
            npairs = npairs + 1

            enefunc%contact_list(1,npairs) = gromol%pairs(k)%atom_idx1  &
                                             + ioffset
            enefunc%contact_list(2,npairs) = gromol%pairs(k)%atom_idx2  &
                                             + ioffset

            if (gromol%pairs(k)%func == 2) then
              sig = 10.0_wp * gromol%pairs(k)%v
              eps  = JOU2CAL * gromol%pairs(k)%w
              !              write(6,*) k, eps

              if (enefunc%forcefield == ForcefieldKBGO) then
                enefunc%contact_lj12(npairs) = 13.0_wp * eps * (sig ** 12)
                enefunc%contact_lj10(npairs) = 18.0_wp * eps * (sig ** 10)
                enefunc%contact_lj6 (npairs) =  4.0_wp * eps * (sig ** 6)
              else if (enefunc%forcefield == ForcefieldCAGO .or. &
                       enefunc%forcefield == ForcefieldRESIDCG) then
                enefunc%contact_lj12(npairs) = 5.0_wp * eps * (sig ** 12)
                enefunc%contact_lj10(npairs) = 6.0_wp * eps * (sig ** 10)
                enefunc%contact_lj6 (npairs) = 0.0_wp
              else
                enefunc%contact_lj12(npairs) = eps * (sig ** 12)
                enefunc%contact_lj10(npairs) = 0.0_wp
                enefunc%contact_lj6 (npairs) = 2.0_wp * eps * (sig ** 6)
              end if

            ! (enefunc%forcefield == ForcefieldSOFT)
            else if ((gromol%pairs(k)%func == 10) .or. &
                     (gromol%pairs(k)%func == 11) .or. &
                     (gromol%pairs(k)%func == 14) .or. &
                     (gromol%pairs(k)%func == 15)) then

              enefunc%contact_func(npairs) = gromol%pairs(k)%func
              enefunc%contact_lj12(npairs) = JOU2CAL * gromol%pairs(k)%w
              enefunc%contact_lj10(npairs) = 10.0_wp * gromol%pairs(k)%v
              enefunc%contact_lj6 (npairs) = JOU2CAL * gromol%pairs(k)%khh / 100.0_wp

            ! (enefunc%forcefield == ForcefieldSOFT)
            else if ((gromol%pairs(k)%func == 12) .or. &
                     (gromol%pairs(k)%func == 13)) then

              call error_msg( &
                   'Setup_Enefunc_Nonb> func = 12 and 13 are not abailable.')

            !shinobu-edited
            else if (gromol%pairs(k)%func == 21) then
              sig = 10.0_wp * gromol%pairs(k)%r0
              eps  = JOU2CAL * gromol%pairs(k)%khh
              enefunc%contact_lj12(npairs) =  eps * (sig ** 12)
              enefunc%contact_lj10(npairs) =  eps * (sig ** 10)

            else ! combi_rule == 1 or 3
              if (enefunc%forcefield == ForcefieldCAGO) then
                enefunc%contact_lj12(npairs) = gromol%pairs(k)%w *  &
                                               1.0E12_wp * JOU2CAL
                enefunc%contact_lj10(npairs) = gromol%pairs(k)%v *  &
                                               1.0E10_wp * JOU2CAL
                enefunc%contact_lj6 (npairs) = 0.0_wp
              else
                enefunc%contact_lj12(npairs) = gromol%pairs(k)%w *  &
                                               1.0E12_wp * JOU2CAL
                enefunc%contact_lj10(npairs) = 0.0_wp
                enefunc%contact_lj6 (npairs) = gromol%pairs(k)%v *  &
                                               1.0E6_wp * JOU2CAL
              end if
            end if
          end do
        end do
      end do
      call get_loop_index(enefunc%num_contacts, istart, iend)
      enefunc%istart_contact = istart
      enefunc%iend_contact   = iend

      !CK20180225
      !
      ! multi-basin
      !
      if (ene_info%num_basins == grotop%gomodel%num_basins) then
        enefunc%num_basins = ene_info%num_basins
      else
        call error_msg( &
            'Setup_Enefunc_Nonb> number of basins is not correct.')
      end if
     
      enefunc%mix_temperature = ene_info%mix_temperature
      enefunc%mix_beta        = 1.0_wp/(ene_info%mix_temperature*KBOLTZ)

      if (enefunc%num_basins > 1) then
        call alloc_enefunc(enefunc, EneFuncMultiWork, enefunc%num_basins ,  &
                           natom)
        enefunc%basinenergy(1:enefunc%num_basins) =  &
                      ene_info%basinenergy(1:enefunc%num_basins)

        nmcont = 0
        do i = 1, grotop%num_molss
          gromol => grotop%molss(i)%moltype%mol
          do j = 1, grotop%molss(i)%count
            nmcont = nmcont + gromol%num_mcontacts
          end do
        end do
     
        enefunc%num_multi_contacts = nmcont
        call alloc_enefunc(enefunc, EneFuncMultiCntc, nmcont)

        nmcont = 0
        natom  = 0
        do i = 1, grotop%num_molss
          gromol => grotop%molss(i)%moltype%mol
          ioffset = natom
          natom   = natom + gromol%num_atoms
          do j = 1, grotop%molss(i)%count
            do k = 1, gromol%num_mcontacts
              nmcont = nmcont + 1
              enefunc%multi_contact_list(1,nmcont) = &
                                  gromol%mcontact(k)%atom_idx1 + ioffset
              enefunc%multi_contact_list(2,nmcont) = &
                                  gromol%mcontact(k)%atom_idx2 + ioffset

              enefunc%multi_contact_model(nmcont) = gromol%mcontact(k)%model

              if (gromol%mcontact(k)%func == 2) then
                sig = 10.0_wp * gromol%mcontact(k)%v
                eps = JOU2CAL * gromol%mcontact(k)%w
           
                if (enefunc%forcefield == ForcefieldKBGO) then
                  enefunc%multi_contact_lj12(nmcont) = 13.0_wp * eps * (sig**12)
                  enefunc%multi_contact_lj10(nmcont) = 18.0_wp * eps * (sig**10)
                  enefunc%multi_contact_lj6 (nmcont) =  4.0_wp * eps * (sig**6)
                else if (enefunc%forcefield == ForcefieldCAGO) then
                  enefunc%multi_contact_lj12(nmcont) = 5.0_wp * eps * (sig**12)
                  enefunc%multi_contact_lj10(nmcont) = 6.0_wp * eps * (sig**10)
                  enefunc%multi_contact_lj6 (nmcont) = 0.0_wp
                else 
                  enefunc%multi_contact_lj12(nmcont) = eps * (sig**12)
                  enefunc%multi_contact_lj10(nmcont) = 0.0_wp
                  enefunc%multi_contact_lj6 (nmcont) = 2.0_wp * eps * (sig**6)
                end if
             
              else ! combi_rule == 1 or 3
                if (enefunc%forcefield == ForcefieldCAGO) then
                  enefunc%multi_contact_lj12(nmcont) = gromol%mcontact(k)%w *  &
                                                 1.0E12_wp * JOU2CAL
                  enefunc%multi_contact_lj10(nmcont) = gromol%mcontact(k)%v *  &
                                                 1.0E10_wp * JOU2CAL
                  enefunc%multi_contact_lj6 (nmcont) = 0.0_wp
                else
                  enefunc%multi_contact_lj12(nmcont) = gromol%mcontact(k)%w *  &
                                                 1.0E12_wp * JOU2CAL
                  enefunc%multi_contact_lj10(nmcont) = 0.0_wp
                  enefunc%multi_contact_lj6 (nmcont) = gromol%mcontact(k)%v *  &
                                                 1.0E6_wp * JOU2CAL
                end if
              end if
            end do
          end do
        end do
        call get_loop_index(enefunc%num_multi_contacts, istart, iend)
        enefunc%istart_multi_contact = istart
        enefunc%iend_multi_contact   = iend

      end if

    end if

    ! check # of exclusion level
    !
    if (enefunc%forcefield == ForcefieldGROMARTINI) then

      !TODO

      enefunc%excl_level = -1

      do i = 1, grotop%num_molss
        excl_level = grotop%molss(i)%moltype%exclude_nbon
        if (enefunc%excl_level == -1) then
          enefunc%excl_level = excl_level
        else if (enefunc%excl_level /= excl_level) then
          call error_msg( &
               'Setup_Enefunc_Nonb> multiple "exclude_nbon" is not supported.')
        end if

      end do

    end if

    ! create exclusion list & exclude neighbours Nx
    !
    nexcl = 0
    natom = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        nexcl = nexcl + gromol%num_excls
        natom = natom + gromol%num_atoms
      end do
    end do

    ! ~CG~ 3SPN.2C DNA: add num of base-stacking pairs to nexcl
    ! 
    if (enefunc%forcefield == ForcefieldRESIDCG) then
      nexcl = nexcl + enefunc%num_base_stack
    end if

    allocate(excls(2,nexcl), nrexc(natom))

    nexcl = 0
    natom = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        ioffset = natom
        natom   = natom + gromol%num_atoms

        ! exclusion list
        do k = 1, gromol%num_excls

          nexcl = nexcl + 1
          excls(1,nexcl) = gromol%excls(k)%atom_idx1 + ioffset
          excls(2,nexcl) = gromol%excls(k)%atom_idx2 + ioffset
        end do

        ! exclude neighbours Nx
        nrexc(1+ioffset:gromol%num_atoms+ioffset) = &
             grotop%molss(i)%moltype%exclude_nbon

      end do
    end do

    ! ~CG~ 3SPN.2C DNA: Add base-steps into excl list
    ! 
    if (enefunc%forcefield == ForcefieldRESIDCG) then
      do i = 1, enefunc%num_base_stack
        nexcl = nexcl + 1
        excls(1, nexcl) = enefunc%base_stack_list(2, i)
        excls(2, nexcl) = enefunc%base_stack_list(3, i)
      end do
    end if

    ! treatment for 1-2, 1-3, 1-4 interactions
    !
    if (enefunc%forcefield == ForcefieldKBGO .or. &
        enefunc%forcefield == ForcefieldAAGO .or. &
        enefunc%forcefield == ForcefieldCAGO .or. &
        enefunc%forcefield == ForcefieldSOFT .or. &
        enefunc%forcefield == ForcefieldRESIDCG) then
      !shinobu-edited
      call count_nonb_excl_go(molecule, enefunc, excls, nrexc)

      ! set the starting index of nonb_excl list
      ! 
      allocate(enefunc%cg_istart_nonb_excl(molecule%num_atoms), stat=alloc_stat)
      if (alloc_stat /=0) &
          call error_msg_alloc
      j = 1
      do i = 1, molecule%num_atoms
        enefunc%cg_istart_nonb_excl(i) = j
        j = j + enefunc%num_nonb_excl(i)
      end do

    else if (.not.ene_info%table) then
      call count_nonb_excl(molecule, enefunc, excls, nrexc)
    end if

    deallocate(excls, nrexc)

    return

  end subroutine setup_enefunc_nonb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cg_contact
  !> @brief        define contact in CG energy functions
  !! @authors      CT
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cg_contact(ene_info, grotop, molecule, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: eps, sig
    integer                  :: i, j, k
    integer                  :: natom, npairs, ioffset
    integer                  :: istart, iend

    type(s_grotop_mol), pointer :: gromol


    ! ================================
    ! parameters for CG native-contact
    ! ================================
    !
    npairs = 0
    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        npairs  = npairs + gromol%num_pairs
      end do
    end do

    call alloc_enefunc(enefunc, EneFuncCntc, npairs)
    enefunc%num_contacts = npairs

    natom  = 0
    npairs = 0
    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        ioffset = natom
        natom   = natom + gromol%num_atoms
        do k = 1, gromol%num_pairs
          npairs = npairs + 1

          enefunc%contact_list(1,npairs) = gromol%pairs(k)%atom_idx1 + ioffset
          enefunc%contact_list(2,npairs) = gromol%pairs(k)%atom_idx2 + ioffset

          if (gromol%pairs(k)%func == 2) then

            sig = 10.0_wp * gromol%pairs(k)%v
            eps = JOU2CAL * gromol%pairs(k)%w

            enefunc%contact_lj12(npairs) = 5.0_wp * eps * (sig ** 12)
            enefunc%contact_lj10(npairs) = 6.0_wp * eps * (sig ** 10)
            enefunc%contact_lj6 (npairs) = 0.0_wp

          else

            call error_msg( &
                'Setup_Enefunc_CG > only func = 2 is abailable for AICG2P contact.')

          end if

        end do                ! k
      end do                  ! j
    end do                    ! i

    call get_loop_index(enefunc%num_contacts, istart, iend)
    enefunc%istart_contact = istart
    enefunc%iend_contact   = iend

    return

  end subroutine setup_enefunc_cg_contact

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cg_exv
  !> @brief        define exv term in CG energy functions
  !! @authors      CT
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cg_exv(ene_info, grotop, molecule, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: vi, wi
    integer                  :: nnonb, i

    type(s_grotop_mol), pointer :: gromol


    enefunc%num_atom_cls = grotop%num_atomtypes
    enefunc%nonb_func    = grotop%defaults%nonb_func

    ! set exv lennard-jones parameters
    !
    nnonb = enefunc%num_atom_cls

    call alloc_enefunc(enefunc, EneFuncNbon, nnonb)

    ! =============================
    ! parameters for CG general exv
    ! =============================
    ! 
    do i = 1, nnonb
      vi = grotop%atomtypes(i)%v
      wi = grotop%atomtypes(i)%w
      enefunc%cg_exv_sig_half(i) = vi * 5.0_wp * ene_info%cg_exv_sigma_scaling
      enefunc%cg_exv_eps_sqrt(i) = sqrt(wi  * JOU2CAL)
    end do
    !
    do i = 1, molecule%num_atoms
      enefunc%param_epsilon(i) = enefunc%cg_exv_eps_sqrt(enefunc%atom_cls(i))
      enefunc%param_sigma  (i) = enefunc%cg_exv_sig_half(enefunc%atom_cls(i))
    end do

    return

  end subroutine setup_enefunc_cg_exv

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cg_nonlocal_exclusion
  !> @brief        define exclusion list for CG nonlocal terms
  !! @authors      CT
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cg_nonlocal_exclusion(ene_info, grotop, molecule, &
                                                 enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                     :: natom
    integer                     :: i, j, k, l
    integer                     :: i1, i2
    integer                     :: n_tmp, ioffset
    integer                     :: total_num_excl
    integer                     :: max_num_excl
    integer, parameter          :: max_num_local = 20
    integer                     :: alloc_stat, dealloc_stat
    logical                     :: duplicate

    integer, allocatable        :: num_nonb_excl(:)
    integer, allocatable        :: nonb_excl_list(:,:)

    type(s_grotop_mol), pointer :: gromol


    natom = molecule%num_atoms

    allocate(num_nonb_excl(natom), stat=alloc_stat)
    if (alloc_stat /=0) call error_msg_alloc

    ! ==================================================
    ! Step 1: estimate exclusion list size for each atom
    ! ==================================================
    !
    total_num_excl = 0 
    num_nonb_excl(1:natom) = 0
    !
    ! count native contacts
    ! 
    do i = 1, enefunc%num_contacts
      i1 = enefunc%contact_list(1, i)
      i2 = enefunc%contact_list(2, i)
      if (i1 < i2) then
        num_nonb_excl(i1) = num_nonb_excl(i1) + 1
      else
        num_nonb_excl(i2) = num_nonb_excl(i2) + 1
      end if
    end do
    !
    ! count user-defined exclusion list
    ! 
    n_tmp = 0
    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        ioffset = n_tmp
        n_tmp   = n_tmp + gromol%num_atoms

        do k = 1, gromol%num_excls
          i1 = gromol%excls(k)%atom_idx1 + ioffset
          i2 = gromol%excls(k)%atom_idx2 + ioffset
          if (i1 < i2) then
            num_nonb_excl(i1) = num_nonb_excl(i1) + 1
          else
            num_nonb_excl(i2) = num_nonb_excl(i2) + 1
          end if
        end do

      end do
    end do
    !
    ! find out the maximum of num_nonb_excl
    ! 
    max_num_excl = 0
    do i = 1, natom
      if (max_num_excl < num_nonb_excl(i)) then
        max_num_excl = num_nonb_excl(i)
      end if
    end do
    max_num_excl = max_num_excl + max_num_local
    !
    ! -----------------------
    ! allocate nonb_excl_list
    ! -----------------------
    !
    allocate(nonb_excl_list(max_num_excl,natom),stat=alloc_stat)
    if (alloc_stat /=0) call error_msg_alloc

    ! =================================================
    ! Step 2: feed local and contacts to exclusion list 
    ! =================================================
    !
    num_nonb_excl(1:natom) = 0
    !
    ! ---------------------------
    ! add bonds to exclusion list
    ! ---------------------------
    ! 
    do i = 1, enefunc%num_bonds
      i1 = enefunc%bond_list(1, i)
      i2 = enefunc%bond_list(2, i)

      if (i1 < i2) then
        n_tmp = num_nonb_excl(i1) + 1
        num_nonb_excl(i1) = n_tmp
        nonb_excl_list(n_tmp, i1) = i2
      else
        n_tmp = num_nonb_excl(i2) + 1
        num_nonb_excl(i2) = n_tmp
        nonb_excl_list(n_tmp, i2) = i1
      end if
      total_num_excl = total_num_excl + 1
    end do
    ! 
    do i = 1, enefunc%num_bonds_quartic
      i1 = enefunc%bond_quartic_list(1, i)
      i2 = enefunc%bond_quartic_list(2, i)

      if (i1 < i2) then
        n_tmp = num_nonb_excl(i1) + 1
        num_nonb_excl(i1) = n_tmp
        nonb_excl_list(n_tmp, i1) = i2
      else
        n_tmp = num_nonb_excl(i2) + 1
        num_nonb_excl(i2) = n_tmp
        nonb_excl_list(n_tmp, i2) = i1
      end if
      total_num_excl = total_num_excl + 1
    end do
    !
    ! ----------------------------
    ! add angles to exclusion list
    ! ----------------------------
    !
    do i = 1, enefunc%num_angles
      i1 = enefunc%angl_list(1, i)
      i2 = enefunc%angl_list(3, i)

      if (i1 > i2) then
        j  = i1
        i1 = i2
        i2 = j
      end if

      duplicate = .false.
      do j = 1, num_nonb_excl(i1)
        if (i2 == nonb_excl_list(j, i1)) then
          duplicate = .true.
          exit
        end if
      end do

      if (.not. duplicate) then
        n_tmp = num_nonb_excl(i1) + 1
        num_nonb_excl(i1) = n_tmp
        nonb_excl_list(n_tmp, i1) = i2
        total_num_excl = total_num_excl + 1
      end if
    end do
    ! 
    do i = 1, enefunc%num_angflex
      i1 = enefunc%anglflex_list(1, i)
      i2 = enefunc%anglflex_list(3, i)

      if (i1 > i2) then
        j  = i1
        i1 = i2
        i2 = j
      end if

      duplicate = .false.
      do j = 1, num_nonb_excl(i1)
        if (i2 == nonb_excl_list(j, i1)) then
          duplicate = .true.
          exit
        end if
      end do

      if (.not. duplicate) then
        n_tmp = num_nonb_excl(i1) + 1
        num_nonb_excl(i1) = n_tmp
        nonb_excl_list(n_tmp, i1) = i2
        total_num_excl = total_num_excl + 1
      end if
    end do
    !
    ! -------------------------------
    ! add dihedrals to exclusion list
    ! -------------------------------
    !
    do i = 1, enefunc%num_dihedrals
      i1 = enefunc%dihe_list(1, i)
      i2 = enefunc%dihe_list(4, i)

      if (i1 > i2) then
        j  = i1
        i1 = i2
        i2 = j
      end if

      duplicate = .false.
      do j = 1, num_nonb_excl(i1)
        if (i2 == nonb_excl_list(j, i1)) then
          duplicate = .true.
          exit
        end if
      end do

      if (.not. duplicate) then
        n_tmp = num_nonb_excl(i1) + 1
        num_nonb_excl(i1) = n_tmp
        nonb_excl_list(n_tmp, i1) = i2
        total_num_excl = total_num_excl + 1
      end if
    end do
    !
    do i = 1, enefunc%num_dihedflex
      i1 = enefunc%diheflex_list(1, i)
      i2 = enefunc%diheflex_list(4, i)

      if (i1 > i2) then
        j  = i1
        i1 = i2
        i2 = j
      end if

      duplicate = .false.
      do j = 1, num_nonb_excl(i1)
        if (i2 == nonb_excl_list(j, i1)) then
          duplicate = .true.
          exit
        end if
      end do

      if (.not. duplicate) then
        n_tmp = num_nonb_excl(i1) + 1
        num_nonb_excl(i1) = n_tmp
        nonb_excl_list(n_tmp, i1) = i2
        total_num_excl = total_num_excl + 1
      end if
    end do
    !
    ! -----------------------------------------
    ! add base-stacking terms to exclusion list
    ! -----------------------------------------
    ! 
    do i = 1, enefunc%num_base_stack
      i1 = enefunc%base_stack_list(2, i)
      i2 = enefunc%base_stack_list(3, i)

      if (i1 > i2) then
        j  = i1
        i1 = i2
        i2 = j
      end if

      duplicate = .false.
      do j = 1, num_nonb_excl(i1)
        if (i2 == nonb_excl_list(j, i1)) then
          duplicate = .true.
          exit
        end if
      end do

      if (.not. duplicate) then
        n_tmp = num_nonb_excl(i1) + 1
        num_nonb_excl(i1) = n_tmp
        nonb_excl_list(n_tmp, i1) = i2
        total_num_excl = total_num_excl + 1
      end if
    end do
    ! 
    ! -------------------------------------
    ! add native contacts to exclusion list
    ! -------------------------------------
    ! 
    do i = 1, enefunc%num_contacts
      i1 = enefunc%contact_list(1, i)
      i2 = enefunc%contact_list(2, i)

      if (i1 > i2) then
        j  = i1
        i1 = i2
        i2 = j
      end if

      ! there is no overlap between contacts and local terms
      ! (hopefully...)
      ! 
      n_tmp = num_nonb_excl(i1) + 1
      num_nonb_excl(i1) = n_tmp
      nonb_excl_list(n_tmp, i1) = i2
      total_num_excl = total_num_excl + 1
    end do
    !
    ! -------------------------------
    ! add user-defined exclusion list
    ! -------------------------------
    ! 
    n_tmp = 0
    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        ioffset = n_tmp
        n_tmp   = n_tmp + gromol%num_atoms

        do k = 1, gromol%num_excls
          i1 = gromol%excls(k)%atom_idx1 + ioffset
          i2 = gromol%excls(k)%atom_idx2 + ioffset
          
          if (i1 > i2) then
            l  = i1
            i1 = i2
            i2 = l
          end if

          duplicate = .false.
          do l = 1, num_nonb_excl(i1)
            if (i2 == nonb_excl_list(l, i1)) then
              duplicate = .true.
              exit
            end if
          end do

          if (.not. duplicate) then
            l = num_nonb_excl(i1) + 1
            num_nonb_excl(i1) = l
            nonb_excl_list(l, i1) = i2
            total_num_excl = total_num_excl + 1
          end if
        end do                  ! do k
      end do                    ! do j
    end do                      ! do i

    ! ========================================
    ! Step 3: allocate enefunc exclusion lists
    ! ========================================
    !
    allocate(enefunc%num_nonb_excl(natom), stat=alloc_stat)
    if (alloc_stat /=0) call error_msg_alloc
    ! 
    allocate(enefunc%nonb_excl_list(total_num_excl), stat = alloc_stat)
    if (alloc_stat /=0) call error_msg_alloc
    ! 
    allocate(enefunc%cg_istart_nonb_excl(natom), stat=alloc_stat)
    if (alloc_stat /=0) call error_msg_alloc
    !
    enefunc%num_nonb_excl(1:natom) = num_nonb_excl(1:natom)
    !
    n_tmp = 0
    do i = 1, natom
      enefunc%cg_istart_nonb_excl(i) = n_tmp + 1
      do j = 1, num_nonb_excl(i)
        n_tmp = n_tmp + 1
        enefunc%nonb_excl_list(n_tmp) = nonb_excl_list(j, i)
      end do
    end do

    deallocate(num_nonb_excl,  stat = dealloc_stat)
    deallocate(nonb_excl_list, stat = dealloc_stat)


    ! ---------
    ! Other ...
    ! ---------
    !
    allocate(enefunc%num_nb14_calc(natom), stat=alloc_stat)
    enefunc%num_nb14_calc(1:natom) = 0

    return

  end subroutine setup_enefunc_cg_nonlocal_exclusion

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_morph
  !> @brief        define morphing term in potential energy function
  !! @authors      CK
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_morph(grotop, enefunc)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: natom, num_morph_bb, num_morph_sc
    integer                  :: ioffset, nmorph_bb, nmorph_sc
    integer                  :: istart, iend

    type(s_grotop_mol), pointer :: gromol


    if (.not. enefunc%morph_flag) return

    num_morph_bb=0
    num_morph_sc=0
    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        num_morph_bb = num_morph_bb + gromol%num_morph_bb
        num_morph_sc = num_morph_sc + gromol%num_morph_sc
      end do
    end do
    if (num_morph_bb == 0 .and. num_morph_bb == 0) return

    enefunc%num_morph_bb = num_morph_bb
    enefunc%num_morph_sc = num_morph_sc
    call alloc_enefunc(enefunc, EneFuncMorph, num_morph_bb, num_morph_sc)

    natom = 0
    nmorph_bb = 0
    nmorph_sc = 0
    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        ioffset = natom
        natom   = natom + gromol%num_atoms

        ! morph list
        do k = 1, gromol%num_morph_bb
          
          nmorph_bb = nmorph_bb + 1
          enefunc%morph_list_bb(1,nmorph_bb) = gromol%morph_bb(k)%atom_idx1 & 
                                            + ioffset
          enefunc%morph_list_bb(2,nmorph_bb) = gromol%morph_bb(k)%atom_idx2 &
                                            + ioffset
          enefunc%morph_dist_bb(nmorph_bb)   = gromol%morph_bb(k)%rmin * 10.0_wp
          enefunc%morph_dist_bb_other(nmorph_bb)= gromol%morph_bb(k)%rmin_other &
                                                      * 10.0_wp
        end do

        do k = 1, gromol%num_morph_sc
          
          nmorph_sc = nmorph_sc + 1
          enefunc%morph_list_sc(1,nmorph_sc) = gromol%morph_sc(k)%atom_idx1 &
                                            + ioffset
          enefunc%morph_list_sc(2,nmorph_sc) = gromol%morph_sc(k)%atom_idx2 &
                                            + ioffset
          enefunc%morph_dist_sc(nmorph_sc)   = gromol%morph_sc(k)%rmin * 10.0_wp

          enefunc%morph_dist_sc_other(nmorph_sc)= gromol%morph_sc(k)%rmin_other &
                                                      * 10.0_wp
        end do

      end do
    end do

    call get_loop_index(enefunc%num_morph_bb, istart, iend)
    enefunc%istart_morph_bb = istart
    enefunc%iend_morph_bb   = iend

    call get_loop_index(enefunc%num_morph_sc, istart, iend)
    enefunc%istart_morph_sc = istart
    enefunc%iend_morph_sc   = iend

    return

  end subroutine setup_enefunc_morph

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_gro_restraints
  !> @brief        define NON-BOND term in potential energy function
  !! @authors      NT
  !! @param[in]    molecule : molecule information
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_gro_restraints(molecule, grotop, enefunc)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_grotop),          intent(in)    :: grotop
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    type(s_enefunc)          :: ef0
    real(wp)                 :: kx, ky, kz
    integer                  :: nposres, npr_atom, max_pr_atom, natom, ioffset
    integer                  :: i, j, k, n, n2, n3, group0, func0, istart, iend

    type(s_grotop_mol), pointer :: gromol


    ! count # of position restraints
    !
    nposres     = 0
    max_pr_atom = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol

      if (gromol%num_posress > 0) &
        nposres = nposres + 1

      npr_atom = grotop%molss(i)%count * gromol%num_posress
      max_pr_atom = max(max_pr_atom,npr_atom)
    end do

    if (nposres == 0) &
      return

    enefunc%restraint_flag = .true.


    ! setup EneFuncRefc
    !
    if (.not. allocated(enefunc%restraint_refcoord)) then
      n = molecule%num_atoms
      call alloc_enefunc(enefunc, EneFuncRefc, n)
      enefunc%restraint_refcoord(1:3,1:n) = molecule%atom_coord(1:3,1:n)
    end if


    ! setup EneFuncRefg
    !
    if (allocated(enefunc%restraint_numatoms)) then

#ifdef QSIMULATE
      ! This is a work around for gcc 7.2 
      if (main_rank) then
        write(MsgOut,'("GROMACS restraint is not permitted with QSimulate")')
      end if
      call error_msg()
#else
      n  = size(enefunc%restraint_numatoms)
      n2 = size(enefunc%restraint_atomlist(:,1))
      call alloc_enefunc(ef0, EneFuncRefg, n, n2)
      ef0%restraint_numatoms     (1:n) = enefunc%restraint_numatoms     (1:n)
      ef0%restraint_atomlist(1:n2,1:n) = enefunc%restraint_atomlist(1:n2,1:n)
      ef0%restraint_masscoef(1:n2,1:n) = enefunc%restraint_masscoef(1:n2,1:n)

      call alloc_enefunc(enefunc, EneFuncRefg, n+nposres, max(n2,max_pr_atom))
      enefunc%restraint_numatoms     (1:n) = ef0%restraint_numatoms     (1:n)
      enefunc%restraint_atomlist(1:n2,1:n) = ef0%restraint_atomlist(1:n2,1:n)
      enefunc%restraint_masscoef(1:n2,1:n) = ef0%restraint_masscoef(1:n2,1:n)

      call dealloc_enefunc(ef0, EneFuncRefg)
#endif

    else

      n = 0
      call alloc_enefunc(enefunc, EneFuncRefg, nposres, max_pr_atom)

    end if


    group0 = n

    natom   = 0
    nposres = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol

      if (gromol%num_posress > 0) &
        nposres = nposres + 1

      npr_atom = 0
      do j = 1, grotop%molss(i)%count
        ioffset = natom
        natom   = natom + gromol%num_atoms

        do k = 1, gromol%num_posress
          npr_atom = npr_atom + 1
          enefunc%restraint_atomlist(npr_atom, group0+nposres) = &
               gromol%posress(k)%atom_idx + ioffset

          if (k > 1 .and. main_rank) then
            if (gromol%posress(k-1)%kx /= gromol%posress(k)%kx .or. &
                gromol%posress(k-1)%ky /= gromol%posress(k)%ky .or. &
                gromol%posress(k-1)%kz /= gromol%posress(k)%kz) then
              write(MsgOut,'(a)') 'Setup_Enefunc_Gro_Restraints> WARNING:'
              write(MsgOut,'(a,a,a)') &
   '   different restraint constant between foreach atoms is not supported. [',&
   trim(grotop%molss(i)%moltype%name), ']'
              write(MsgOut,'(a)') ' '
            end if
          end if

        end do

      end do

      if (npr_atom > 0) &
        enefunc%restraint_numatoms(group0+nposres) = npr_atom

    end do


    ! setup EneFuncReff
    !
    if (allocated(enefunc%restraint_kind)) then

#ifdef QSIMULATE
      ! This is a work around for gcc 7.2 
      if (main_rank) then
        write(MsgOut,'("GROMACS restraint is not permitted with QSimulate")')
      end if
      call error_msg()
#else
      n  = size(enefunc%restraint_kind)
      n2 = size(enefunc%restraint_grouplist(:,1))
      n3 = max(int(n2/2),1)
      call alloc_enefunc(ef0, EneFuncReff, n, n2)
      ef0%restraint_kind          (1:n) = enefunc%restraint_kind          (1:n)
      ef0%restraint_grouplist(1:n2,1:n) = enefunc%restraint_grouplist(1:n2,1:n)
      ef0%restraint_const    (1:4, 1:n) = enefunc%restraint_const    (1:4, 1:n)
      ef0%restraint_ref      (1:2, 1:n) = enefunc%restraint_ref      (1:2, 1:n)
      ef0%restraint_funcgrp       (1:n) = enefunc%restraint_funcgrp       (1:n)
      ef0%restraint_exponent_func (1:n) = enefunc%restraint_exponent_func (1:n)
      ef0%restraint_exponent_dist (1:n3,1:n) &
                                   = enefunc%restraint_exponent_dist (1:n3,1:n)
      ef0%restraint_weight_dist   (1:n3,1:n) &
                                   = enefunc%restraint_weight_dist   (1:n3,1:n)

      call alloc_enefunc(enefunc, EneFuncReff, n+nposres, max(n2,1))
      enefunc%restraint_kind          (1:n) = ef0%restraint_kind          (1:n)
      enefunc%restraint_grouplist(1:n2,1:n) = ef0%restraint_grouplist(1:n2,1:n)
      enefunc%restraint_const    (1:4, 1:n) = ef0%restraint_const    (1:4, 1:n)
      enefunc%restraint_ref      (1:2, 1:n) = ef0%restraint_ref      (1:2, 1:n)
      enefunc%restraint_funcgrp       (1:n) = ef0%restraint_funcgrp       (1:n)
      enefunc%restraint_exponent_func (1:n) = ef0%restraint_exponent_func (1:n)
      enefunc%restraint_exponent_dist (1:n3,1:n) &
                                       = ef0%restraint_exponent_dist (1:n3,1:n)
      enefunc%restraint_weight_dist   (1:n3,1:n) &
                                       = ef0%restraint_weight_dist   (1:n3,1:n)

      call dealloc_enefunc(ef0, EneFuncReff)
#endif

    else

      n = 0
      call alloc_enefunc(enefunc, EneFuncReff, nposres, 1)

    end if


    func0 = n

    nposres = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol

      if (gromol%num_posress == 0) &
        cycle

      nposres = nposres + 1

      kx = gromol%posress(1)%kx * 0.01_wp * JOU2CAL * 0.5_wp
      ky = gromol%posress(1)%ky * 0.01_wp * JOU2CAL * 0.5_wp
      kz = gromol%posress(1)%kz * 0.01_wp * JOU2CAL * 0.5_wp

      enefunc%restraint_kind         (func0+nposres) = RestraintsFuncPOSI
      enefunc%restraint_funcgrp      (func0+nposres) = 1
      enefunc%restraint_grouplist  (1,func0+nposres) = group0 + nposres
      enefunc%restraint_const      (1,func0+nposres) = 1.0_wp
      enefunc%restraint_const      (2,func0+nposres) = kx
      enefunc%restraint_const      (3,func0+nposres) = ky
      enefunc%restraint_const      (4,func0+nposres) = kz
      enefunc%restraint_exponent_func(func0+nposres) = 2

    end do

    enefunc%num_restraintgroups = group0+nposres
    enefunc%num_restraintfuncs  = func0+nposres


    ! define restraint functions for each processor
    !
    call get_loop_index(enefunc%num_restraintfuncs, istart, iend)

    enefunc%istart_restraint = istart
    enefunc%iend_restraint   = iend


    ! summary of setup enefunc_gro_restraint
    !
    if (main_rank) then

      write(MsgOut,'(A)')&
           'Setup_Enefunc_Gro_Restraints> Setup restraint functions'

      do i = func0+1, enefunc%num_restraintfuncs

        write(MsgOut,'(A,I5,A,I5)') &
          ' func  = ', i, ' kind  = ', enefunc%restraint_kind(i)

        ! summary for positional restraint
        !
        write(MsgOut,'(A,4F8.3)') &
             ' const(total, x, y, z) = ', enefunc%restraint_const(1:4,i)
        write(MsgOut,'(A,I5)') ' exponend of function = ', &
                                          enefunc%restraint_exponent_func(i)

        write(MsgOut,'(A,I5)') ' # of groups  = ', enefunc%restraint_funcgrp(i)
        write(MsgOut,'(" grouplist: ",$)')
        do j = 1, enefunc%restraint_funcgrp(i)
          write(MsgOut,'(i3,$)') enefunc%restraint_grouplist(j,i)
          if (mod(j,20) == 0 .and. j /= enefunc%restraint_funcgrp(i) )  &
            write(MsgOut,'(A)') ''
        end do
        write(MsgOut,'(A)') ''

      end do

      write(MsgOut,*) ''

    end if

    return

  end subroutine setup_enefunc_gro_restraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_vsite2
  !> @brief        define virtual site 2 term in potential energy function
  !! @authors      NT
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_vsite2(grotop, enefunc)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: natom, nsite, ioffset
    integer                  :: istart, iend

    type(s_grotop_mol), pointer :: gromol


    nsite = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        do k = 1, gromol%num_vsites2
          if (gromol%vsites2(k)%func /= 1) cycle
          nsite = nsite + 1
        end do
      end do
    end do

    call alloc_enefunc(enefunc, EneFuncVsite2, nsite)

    natom = 0
    nsite = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        ioffset = natom
        natom   = natom + gromol%num_atoms

        do k = 1, gromol%num_vsites2
          if (gromol%vsites2(k)%func /= 1) cycle
          nsite = nsite + 1
          enefunc%vsite2_list(1, nsite) = gromol%vsites2(k)%atom_idx1 + ioffset
          enefunc%vsite2_list(2, nsite) = gromol%vsites2(k)%atom_idx2 + ioffset
          enefunc%vsite2_list(3, nsite) = gromol%vsites2(k)%atom_idx3 + ioffset
          enefunc%vsite2_a   (   nsite) = gromol%vsites2(k)%a
        end do
      end do
    end do

    enefunc%num_vsite2 = nsite


    call get_loop_index(enefunc%num_vsite2, istart, iend)

    enefunc%istart_vsite2 = istart
    enefunc%iend_vsite2   = iend

    return

  end subroutine setup_enefunc_vsite2

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_vsite3
  !> @brief        define virtual site 3 term in potential energy function
  !! @authors      NT
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_vsite3(grotop, enefunc)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: natom, nsite, ioffset
    integer                  :: istart, iend

    type(s_grotop_mol), pointer :: gromol


    nsite = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        do k = 1, gromol%num_vsites3
          if (gromol%vsites3(k)%func /= 1) cycle
          nsite = nsite + 1
        end do
      end do
    end do

    call alloc_enefunc(enefunc, EneFuncVsite3, nsite)

    natom = 0
    nsite = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        ioffset = natom
        natom   = natom + gromol%num_atoms

        do k = 1, gromol%num_vsites3
          if (gromol%vsites3(k)%func /= 1) cycle
          nsite = nsite + 1
          enefunc%vsite3_list(1, nsite) = gromol%vsites3(k)%atom_idx1 + ioffset
          enefunc%vsite3_list(2, nsite) = gromol%vsites3(k)%atom_idx2 + ioffset
          enefunc%vsite3_list(3, nsite) = gromol%vsites3(k)%atom_idx3 + ioffset
          enefunc%vsite3_list(4, nsite) = gromol%vsites3(k)%atom_idx4 + ioffset
          enefunc%vsite3_a   (   nsite) = gromol%vsites3(k)%a
          enefunc%vsite3_b   (   nsite) = gromol%vsites3(k)%b
        end do
      end do
    end do

    enefunc%num_vsite3 = nsite


    call get_loop_index(enefunc%num_vsite3, istart, iend)

    enefunc%istart_vsite3 = istart
    enefunc%iend_vsite3   = iend

    return

  end subroutine setup_enefunc_vsite3

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_vsite3fd
  !> @brief        define virtual site 3fd term in potential energy function
  !! @authors      NT
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_vsite3fd(grotop, enefunc)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: natom, nsite, ioffset
    integer                  :: istart, iend

    type(s_grotop_mol), pointer :: gromol


    nsite = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        do k = 1, gromol%num_vsites3
          if (gromol%vsites3(k)%func /= 2) cycle
          nsite = nsite + 1
        end do
      end do
    end do

    call alloc_enefunc(enefunc, EneFuncVsite3fd, nsite)

    natom = 0
    nsite = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        ioffset = natom
        natom   = natom + gromol%num_atoms

        do k = 1, gromol%num_vsites3
          if (gromol%vsites3(k)%func /= 2) cycle
          nsite = nsite + 1
          enefunc%vsite3fd_list(1, nsite) =gromol%vsites3(k)%atom_idx1 + ioffset
          enefunc%vsite3fd_list(2, nsite) =gromol%vsites3(k)%atom_idx2 + ioffset
          enefunc%vsite3fd_list(3, nsite) =gromol%vsites3(k)%atom_idx3 + ioffset
          enefunc%vsite3fd_list(4, nsite) =gromol%vsites3(k)%atom_idx4 + ioffset
          enefunc%vsite3fd_a   (   nsite) =gromol%vsites3(k)%a
          enefunc%vsite3fd_d   (   nsite) =gromol%vsites3(k)%d * 10.0_wp
        end do
      end do
    end do

    enefunc%num_vsite3fd = nsite


    call get_loop_index(enefunc%num_vsite3fd, istart, iend)

    enefunc%istart_vsite3fd = istart
    enefunc%iend_vsite3fd   = iend

    return

  end subroutine setup_enefunc_vsite3fd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_vsite3fad
  !> @brief        define virtual site 3fad term in potential energy function
  !! @authors      NT
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_vsite3fad(grotop, enefunc)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: natom, nsite, ioffset
    integer                  :: istart, iend

    type(s_grotop_mol), pointer :: gromol


    nsite = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        do k = 1, gromol%num_vsites3
          if (gromol%vsites3(k)%func /= 3) cycle
          nsite = nsite + 1
        end do
      end do
    end do

    call alloc_enefunc(enefunc, EneFuncVsite3fad, nsite)

    natom = 0
    nsite = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        ioffset = natom
        natom   = natom + gromol%num_atoms

        do k = 1, gromol%num_vsites3
          if (gromol%vsites3(k)%func /= 3) cycle
          nsite = nsite + 1
          enefunc%vsite3fad_list(1, nsite) =gromol%vsites3(k)%atom_idx1 +ioffset
          enefunc%vsite3fad_list(2, nsite) =gromol%vsites3(k)%atom_idx2 +ioffset
          enefunc%vsite3fad_list(3, nsite) =gromol%vsites3(k)%atom_idx3 +ioffset
          enefunc%vsite3fad_list(4, nsite) =gromol%vsites3(k)%atom_idx4 +ioffset
          enefunc%vsite3fad_theta(  nsite) =gromol%vsites3(k)%theta * RAD
          enefunc%vsite3fad_d   (   nsite) =gromol%vsites3(k)%d * 10.0_wp
        end do
      end do
    end do

    enefunc%num_vsite3fad = nsite


    call get_loop_index(enefunc%num_vsite3fad, istart, iend)

    enefunc%istart_vsite3fad = istart
    enefunc%iend_vsite3fad   = iend

    return

  end subroutine setup_enefunc_vsite3fad

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_vsite3out
  !> @brief        define virtual site 3out term in potential energy function
  !! @authors      NT
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_vsite3out(grotop, enefunc)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: natom, nsite, ioffset
    integer                  :: istart, iend

    type(s_grotop_mol), pointer :: gromol


    nsite = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        do k = 1, gromol%num_vsites3
          if (gromol%vsites3(k)%func /= 4) cycle
          nsite = nsite + 1
        end do
      end do
    end do

    call alloc_enefunc(enefunc, EneFuncVsite3out, nsite)

    natom = 0
    nsite = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        ioffset = natom
        natom   = natom + gromol%num_atoms

        do k = 1, gromol%num_vsites3
          if (gromol%vsites3(k)%func /= 4) cycle
          nsite = nsite + 1
          enefunc%vsite3out_list(1, nsite) =gromol%vsites3(k)%atom_idx1 +ioffset
          enefunc%vsite3out_list(2, nsite) =gromol%vsites3(k)%atom_idx2 +ioffset
          enefunc%vsite3out_list(3, nsite) =gromol%vsites3(k)%atom_idx3 +ioffset
          enefunc%vsite3out_list(4, nsite) =gromol%vsites3(k)%atom_idx4 +ioffset
          enefunc%vsite3out_a   (   nsite) =gromol%vsites3(k)%a
          enefunc%vsite3out_b   (   nsite) =gromol%vsites3(k)%b
          enefunc%vsite3out_c   (   nsite) =gromol%vsites3(k)%c * 0.1_wp
        end do
      end do
    end do

    enefunc%num_vsite3out = nsite


    call get_loop_index(enefunc%num_vsite3out, istart, iend)

    enefunc%istart_vsite3out = istart
    enefunc%iend_vsite3out   = iend

    return

  end subroutine setup_enefunc_vsite3out

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_vsite4fdn
  !> @brief        define virtual site 4fdn term in potential energy function
  !! @authors      NT
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_vsite4fdn(grotop, enefunc)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: natom, nsite, ioffset
    integer                  :: istart, iend

    type(s_grotop_mol), pointer :: gromol


    nsite = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        do k = 1, gromol%num_vsites4
          if (gromol%vsites4(k)%func /= 2) cycle
          nsite = nsite + 1
        end do
      end do
    end do

    call alloc_enefunc(enefunc, EneFuncVsite4fdn, nsite)

    natom = 0
    nsite = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        ioffset = natom
        natom   = natom + gromol%num_atoms

        do k = 1, gromol%num_vsites4
          if (gromol%vsites4(k)%func /= 2) cycle
          nsite = nsite + 1
          enefunc%vsite4fdn_list(1, nsite) =gromol%vsites4(k)%atom_idx1 +ioffset
          enefunc%vsite4fdn_list(2, nsite) =gromol%vsites4(k)%atom_idx2 +ioffset
          enefunc%vsite4fdn_list(3, nsite) =gromol%vsites4(k)%atom_idx3 +ioffset
          enefunc%vsite4fdn_list(4, nsite) =gromol%vsites4(k)%atom_idx4 +ioffset
          enefunc%vsite4fdn_list(5, nsite) =gromol%vsites4(k)%atom_idx5 +ioffset
          enefunc%vsite4fdn_a   (   nsite) =gromol%vsites4(k)%a
          enefunc%vsite4fdn_b   (   nsite) =gromol%vsites4(k)%b
          enefunc%vsite4fdn_c   (   nsite) =gromol%vsites4(k)%c * 10.0_wp
        end do
      end do
    end do

    enefunc%num_vsite4fdn = nsite


    call get_loop_index(enefunc%num_vsite4fdn, istart, iend)

    enefunc%istart_vsite4fdn = istart
    enefunc%iend_vsite4fdn   = iend

    return

  end subroutine setup_enefunc_vsite4fdn

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_vsiten
  !> @brief        define virtual site N term in potential energy function
  !! @authors      NT
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_vsiten(grotop, enefunc)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k, m
    integer                  :: natom, nsite, ioffset, maxn
    integer                  :: istart, iend

    type(s_grotop_mol), pointer :: gromol


    nsite = 0
    maxn  = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        do k = 1, gromol%num_vsitesn
          if (gromol%vsitesn(k)%func /= 1 .and. &
              gromol%vsitesn(k)%func /= 2 .and. &
              gromol%vsitesn(k)%func /= 3) cycle
          nsite = nsite + 1
          maxn  = max(maxn, size(gromol%vsitesn(k)%atom_idcs)+1)
        end do
      end do
    end do

    call alloc_enefunc(enefunc, EneFuncVsiten, nsite, maxn)

    natom = 0
    nsite = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        ioffset = natom
        natom   = natom + gromol%num_atoms

        do k = 1, gromol%num_vsites4
          if (gromol%vsitesn(k)%func /= 1 .and. &
              gromol%vsitesn(k)%func /= 2 .and. &
              gromol%vsitesn(k)%func /= 3) cycle
          nsite = nsite + 1
          enefunc%vsiten_n(nsite) = size(gromol%vsitesn(k)%atom_idcs)+1
          enefunc%vsiten_list(1, nsite) = &
                 gromol%vsitesn(k)%atom_idx + ioffset
          do m = 1, enefunc%vsiten_n(nsite)
            enefunc%vsiten_list(m+1, nsite) = &
                 gromol%vsitesn(k)%atom_idcs(m) + ioffset
          end do
        end do
      end do
    end do

    enefunc%num_vsiten = nsite


    call get_loop_index(enefunc%num_vsiten, istart, iend)

    enefunc%istart_vsiten = istart
    enefunc%iend_vsiten   = iend

    return

  end subroutine setup_enefunc_vsiten

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_ez_membrane
  !> @brief        define morphing term in potential energy function
  !! @authors      CK
  !! @param[in]    molecule : molecule information
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_ez_membrane(molecule, grotop, enefunc)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_grotop),          intent(in)    :: grotop
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                     :: i, j, k, iit, itype
    integer                     :: num_memb_types, num_atoms, natom
    type(s_grotop_mol), pointer :: gromol

    num_memb_types = grotop%num_membranetypes
    if (num_memb_types > 0) enefunc%ez_membrane_flag = .true. 
    num_atoms = molecule%num_atoms

    call alloc_enefunc(enefunc, EneFuncEzMembrane, num_atoms)

    natom  = 0
    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        do k = 1, gromol%num_atoms
          natom   = natom + 1
          itype = 0
          iit = 0
          do while (itype == 0 .and. iit <  num_memb_types) 
            iit = iit + 1
            if (gromol%atoms(k)%atom_type .eq. &
                 grotop%ez_membrane(iit)%atom_type) then
              itype = iit
            end if
          end do
          if (itype /= 0) then
            enefunc%ez_membrane_const(natom) = grotop%ez_membrane(iit)%de
            enefunc%ez_membrane_zmin(natom)  = grotop%ez_membrane(iit)%zmin
            enefunc%ez_membrane_func(natom)  = grotop%ez_membrane(iit)%func
            enefunc%ez_membrane_polym(natom) = grotop%ez_membrane(iit)%polym
            if (enefunc%ez_membrane_func(natom) == EzMembraneGaussian) then
              enefunc%ez_membrane_polym(natom) = -0.5_wp/ &
                (grotop%ez_membrane(iit)%polym*grotop%ez_membrane(iit)%polym)
            end if
          else
            enefunc%ez_membrane_func(natom)  = EzMembraneNo
          end if
        end do
      end do
    end do

    return

  end subroutine setup_enefunc_ez_membrane

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    count_nonb_excl
  !> @brief        exclude 1-2, 1-3 interactions
  !! @authors      NT
  !! @param[in]    molecule : molecule information
  !! @param[out]   enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine count_nonb_excl(molecule, enefunc, excls, nrexc)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    integer,                 intent(in)    :: excls(:,:)
    integer,                 intent(in)    :: nrexc(:)

    ! local variables
    integer                  :: i, k, k1, k2, k3, k4
    integer                  :: num_excl, num_excl_total
    integer                  :: num_nb14, num_nb14_total
    integer                  :: max_nb14_num
    integer                  :: alloc_stat, dealloc_stat
    logical                  :: duplicate

    integer,     allocatable :: nonb_excl_list(:,:)
    integer,     allocatable :: nb14_calc_list(:,:)
    real(wp),    allocatable :: nb14_qq_scale (:,:)
    real(wp),    allocatable :: nb14_lj_scale (:,:)


    ! allocate nonbonded exclusion list and 1-4 interaction list
    !
    allocate(enefunc%num_nonb_excl(molecule%num_atoms), stat=alloc_stat)
    if (alloc_stat /=0) &
      call error_msg_alloc

    allocate(nonb_excl_list(max_excl,molecule%num_atoms), stat=alloc_stat)
    if (alloc_stat /=0) &
      call error_msg_alloc

    allocate(enefunc%num_nb14_calc(molecule%num_atoms), stat=alloc_stat)
    if (alloc_stat /=0) &
      call error_msg_alloc

    allocate(nb14_calc_list(MAX_NB14,molecule%num_atoms), &
             nb14_qq_scale (MAX_NB14,molecule%num_atoms), &
             nb14_lj_scale (MAX_NB14,molecule%num_atoms),stat=alloc_stat)
    if (alloc_stat /=0) &
      call error_msg_alloc

    enefunc%num_nonb_excl(1:molecule%num_atoms) = 0
    enefunc%num_nb14_calc(1:molecule%num_atoms) = 0

    ! exclude 1-2 interaction
    !
    num_excl_total = 0
    do k = 1, molecule%num_bonds
      k1 = molecule%bond_list(1,k)
      k2 = molecule%bond_list(2,k)

      if (nrexc(k1) < 1) then
        cycle
      end if

      if (k1 < k2) then
        num_excl = enefunc%num_nonb_excl(k1) + 1
        enefunc%num_nonb_excl(k1) = num_excl
        nonb_excl_list(num_excl,k1) = k2
        num_excl_total = num_excl_total + 1
      else
        num_excl = enefunc%num_nonb_excl(k2) + 1
        enefunc%num_nonb_excl(k2) = num_excl
        nonb_excl_list(num_excl,k2) = k1
        num_excl_total = num_excl_total + 1
      end if
    end do

    ! exclude 1-3 interaction
    !
    do k = 1, molecule%num_angles

      k1 = molecule%angl_list(1,k)
      k3 = molecule%angl_list(3,k)

      if (nrexc(k1) < 2) then
        cycle
      end if

      if (k1 < k3) then
        num_excl = enefunc%num_nonb_excl(k1)
        duplicate = .false.
        do i = 1, enefunc%num_nonb_excl(k1)
          if (k3 == nonb_excl_list(i,k1)) duplicate = .true.
        end do
        if (.not. duplicate) then
          num_excl = num_excl + 1
          enefunc%num_nonb_excl(k1) = num_excl
          nonb_excl_list(num_excl,k1) = k3
          num_excl_total = num_excl_total + 1
        end if
      else
        num_excl = enefunc%num_nonb_excl(k3)
        duplicate = .false.
        do i = 1, enefunc%num_nonb_excl(k3)
          if (k1 == nonb_excl_list(i,k3)) duplicate = .true.
        end do
        if (.not. duplicate) then
          num_excl = num_excl + 1
          enefunc%num_nonb_excl(k3) = num_excl
          nonb_excl_list(num_excl,k3) = k1
          num_excl_total = num_excl_total + 1
        end if
      end if
    end do

    ! count 1-4 interaction
    !
    num_nb14_total = 0
    do k = 1, molecule%num_dihedrals
      k1 = molecule%dihe_list(1,k)
      k4 = molecule%dihe_list(4,k)

      if (nrexc(k1) < 3) then
        cycle
      end if

      if (k1 < k4) then
        num_nb14 = enefunc%num_nb14_calc(k1)
        duplicate = .false.
        do i = 1, enefunc%num_nonb_excl(k1)
          if (k4 == nonb_excl_list(i,k1)) duplicate = .true.
        end do
        do i = 1, num_nb14
          if (k4 == nb14_calc_list(i,k1)) duplicate = .true.
        end do
        if (.not. duplicate) then
          num_nb14 = num_nb14 + 1
          enefunc%num_nb14_calc(k1) = num_nb14
          nb14_calc_list(num_nb14,k1) = k4
          nb14_qq_scale (num_nb14,k1) = enefunc%fudge_qq
          nb14_lj_scale (num_nb14,k1) = enefunc%fudge_lj
          num_nb14_total = num_nb14_total + 1
        end if
      else
        num_nb14 = enefunc%num_nb14_calc(k4)
        duplicate = .false.
        do i = 1, enefunc%num_nonb_excl(k4)
          if (k1 == nonb_excl_list(i,k4)) duplicate = .true.
        end do
        do i = 1, num_nb14
          if (k1 == nb14_calc_list(i,k4)) duplicate = .true.
        end do
        if (.not. duplicate) then
          num_nb14 = num_nb14 + 1
          enefunc%num_nb14_calc(k4) = num_nb14
          nb14_calc_list(num_nb14,k4) = k1
          nb14_qq_scale (num_nb14,k4) = enefunc%fudge_qq
          nb14_lj_scale (num_nb14,k4) = enefunc%fudge_lj
          num_nb14_total = num_nb14_total + 1
        end if
      end if
    end do

    ! exclude exclusions
    !
    do k = 1, size(excls(1,:))
      k1 = excls(1,k)
      k2 = excls(2,k)
      if (k1 < k2) then
        duplicate = .false.
        do i = 1, enefunc%num_nonb_excl(k1)
          if (k2 == nonb_excl_list(i,k1)) duplicate = .true.
        end do

        if (.not. duplicate) then
          num_excl = enefunc%num_nonb_excl(k1) + 1
          enefunc%num_nonb_excl(k1) = num_excl
          nonb_excl_list(num_excl,k1) = k2
          num_excl_total = num_excl_total + 1
        end if
      else
        duplicate = .false.
        do i = 1, enefunc%num_nonb_excl(k2)
          if (k1 == nonb_excl_list(i,k2)) duplicate = .true.
        end do

        if (.not. duplicate) then
          num_excl = enefunc%num_nonb_excl(k2) + 1
          enefunc%num_nonb_excl(k2) = num_excl
          nonb_excl_list(num_excl,k2) = k1
          num_excl_total = num_excl_total + 1
        end if
      end if
    end do


    ! pack 2D-array into 1D-array
    !
    allocate(enefunc%nonb_excl_list(num_excl_total), stat = alloc_stat)
      if (alloc_stat /=0) call error_msg_alloc

    call pack_array_i(molecule%num_atoms, enefunc%num_nonb_excl,  &
                      nonb_excl_list, enefunc%nonb_excl_list)
      deallocate(nonb_excl_list, stat = dealloc_stat)

    !allocate(enefunc%nb14_calc_list(num_nb14_total), &
    !         enefunc%nb14_qq_scale (num_nb14_total), &
    !         enefunc%nb14_lj_scale (num_nb14_total), stat = alloc_stat)
    max_nb14_num = max(1,maxval(enefunc%num_nb14_calc(1:molecule%num_atoms)))
    allocate(enefunc%nb14_calc_list(max_nb14_num,molecule%num_atoms), &
             enefunc%nb14_qq_scale (max_nb14_num,molecule%num_atoms), &
             enefunc%nb14_lj_scale (max_nb14_num,molecule%num_atoms), &
             stat = alloc_stat)
      if (alloc_stat /=0) call error_msg_alloc

    enefunc%nb14_calc_list(1:max_nb14_num,1:molecule%num_atoms) = 0
    enefunc%nb14_qq_scale (1:max_nb14_num,1:molecule%num_atoms) = 0.0_wp
    enefunc%nb14_lj_scale (1:max_nb14_num,1:molecule%num_atoms) = 0.0_wp
    do i = 1, molecule%num_atoms
      enefunc%nb14_calc_list(1:enefunc%num_nb14_calc(i),i) = &
             nb14_calc_list(1:enefunc%num_nb14_calc(i),i)
      enefunc%nb14_qq_scale (1:enefunc%num_nb14_calc(i),i) = &
             nb14_qq_scale (1:enefunc%num_nb14_calc(i),i)
      enefunc%nb14_lj_scale (1:enefunc%num_nb14_calc(i),i) = &
             nb14_lj_scale (1:enefunc%num_nb14_calc(i),i)
    end do
    !call pack_array_i(molecule%num_atoms, enefunc%num_nb14_calc,  &
    !                  nb14_calc_list, enefunc%nb14_calc_list)

    !call pack_array_r(molecule%num_atoms, enefunc%num_nb14_calc,  &
    !                  nb14_qq_scale, enefunc%nb14_qq_scale)

    !call pack_array_r(molecule%num_atoms, enefunc%num_nb14_calc,  &
    !                  nb14_lj_scale, enefunc%nb14_lj_scale)

    deallocate(nb14_calc_list, &
               nb14_qq_scale,  &
               nb14_lj_scale, stat = dealloc_stat)

    return

  end subroutine count_nonb_excl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    count_nonb_excl_go
  !> @brief        exclude 1-2, 1-3, 1-4 interactions
  !! @authors      NT
  !! @param[in]    molecule : molecule information
  !! @param[out]   enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine count_nonb_excl_go(molecule, enefunc, excls, nrexc)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc), target, intent(inout) :: enefunc
    integer,                 intent(in)    :: excls(:,:)
    integer,                 intent(in)    :: nrexc(:)


    ! local variables
    integer                  :: i, k, k1, k2, k3, k4, natom, ndat
    integer                  :: max_contact
    integer                  :: num_excl, num_excl_total
    integer                  :: num_nb14, num_nb14_total
    integer                  :: max_nonb_excl_num, max_nb14_num
    integer                  :: alloc_stat, dealloc_stat
    logical                  :: duplicate

    integer,     allocatable :: nonb_excl_list(:,:)
    integer,     allocatable :: nb14_calc_list(:,:)
    real(wp),    allocatable :: nb14_qq_scale (:,:)
    real(wp),    allocatable :: nb14_lj_scale (:,:)
    integer,     allocatable :: temporary_list(:,:)


    natom = molecule%num_atoms

    ! allocate nonbonded exclusion list and 1-4 interaction list
    !
    allocate(enefunc%table%num_nonb_excl(natom), stat=alloc_stat)
    if (alloc_stat /=0) &
      call error_msg_alloc

    enefunc%table%num_nonb_excl(1:natom) = 0
    do k = 1, size(excls(1,:))
      k1 = excls(1,k)
      k2 = excls(2,k)
      if (k1 < k2) then
        num_excl = enefunc%table%num_nonb_excl(k1) + 1
        enefunc%table%num_nonb_excl(k1) = num_excl
      else
        num_excl = enefunc%table%num_nonb_excl(k2) + 1
        enefunc%table%num_nonb_excl(k2) = num_excl
      end if
    end do
    max_contact = 0
    do k = 1, natom
      if (max_contact < enefunc%table%num_nonb_excl(k)) then
        max_contact = enefunc%table%num_nonb_excl(k)
      end if
    end do
    enefunc%table%num_nonb_excl(1:natom) = 0

    max_contact = max_excl+max_contact
    allocate(nonb_excl_list(max_contact,natom),stat=alloc_stat)
    if (alloc_stat /=0) &
      call error_msg_alloc

    allocate(enefunc%num_nonb_excl(natom), stat=alloc_stat)
    if (alloc_stat /=0) &
      call error_msg_alloc

    allocate(enefunc%num_nb14_calc(natom), stat=alloc_stat)
    if (alloc_stat /=0) &
      call error_msg_alloc

    allocate(nb14_calc_list(MAX_NB14,natom), &
             nb14_qq_scale (MAX_NB14,natom), &
             nb14_lj_scale (MAX_NB14,natom),stat=alloc_stat)
    if (alloc_stat /=0) &
      call error_msg_alloc

    enefunc%num_nonb_excl(1:natom) = 0
    enefunc%num_nb14_calc(1:natom) = 0

    ! exclude 1-2 interaction
    !
    num_excl_total = 0
    do k = 1, molecule%num_bonds
      k1 = molecule%bond_list(1,k)
      k2 = molecule%bond_list(2,k)
      
      if (nrexc(k1) < 1) then
        cycle
      end if

      if (k1 < k2) then
        num_excl = enefunc%table%num_nonb_excl(k1) + 1
        enefunc%table%num_nonb_excl(k1) = num_excl
        nonb_excl_list(num_excl,k1) = k2
        num_excl_total = num_excl_total + 1
      else
        num_excl = enefunc%table%num_nonb_excl(k2) + 1
        enefunc%table%num_nonb_excl(k2) = num_excl
        nonb_excl_list(num_excl,k2) = k1
        num_excl_total = num_excl_total + 1
      end if
    end do

    ! exclude 1-3 interaction
    !

    call create_nbond_list(molecule%bond_list, 3, temporary_list)

    ndat = 0
    if (allocated(temporary_list)) ndat = size(temporary_list(1,:))
    do k = 1, ndat

      k1 = temporary_list(1,k)
      k3 = temporary_list(3,k)

      if (nrexc(k1) < 2) then
        cycle
      end if

      if (k1 < k3) then
        num_excl = enefunc%table%num_nonb_excl(k1)
        duplicate = .false.
        do i = 1, enefunc%table%num_nonb_excl(k1)
          if (k3 == nonb_excl_list(i,k1)) duplicate = .true.
        end do
        if (.not. duplicate) then
          num_excl = num_excl + 1
          enefunc%table%num_nonb_excl(k1) = num_excl
          nonb_excl_list(num_excl,k1) = k3
          num_excl_total = num_excl_total + 1
        end if
      else
        num_excl = enefunc%table%num_nonb_excl(k3)
        duplicate = .false.
        do i = 1, enefunc%table%num_nonb_excl(k3)
          if (k1 == nonb_excl_list(i,k3)) duplicate = .true.
        end do
        if (.not. duplicate) then
          num_excl = num_excl + 1
          enefunc%table%num_nonb_excl(k3) = num_excl
          nonb_excl_list(num_excl,k3) = k1
          num_excl_total = num_excl_total + 1
        end if
      end if
    end do
    if (allocated(temporary_list)) &
       deallocate(temporary_list)

    ! count 1-4 interaction
    !
    call create_nbond_list(molecule%bond_list, 4, temporary_list)
    num_nb14_total = 0

    ndat = 0
    if (allocated(temporary_list)) ndat = size(temporary_list(1,:))

    do k = 1, ndat
      k1 = temporary_list(1,k)
      k4 = temporary_list(4,k)

      if (nrexc(k1) < 3) then
        cycle
      end if

      if (k1 < k4) then
        num_nb14 = enefunc%num_nb14_calc(k1)
        duplicate = .false.
        do i = 1, enefunc%table%num_nonb_excl(k1)
          if (k4 == nonb_excl_list(i,k1)) duplicate = .true.
        end do
        do i = 1, num_nb14
          if (k4 == nb14_calc_list(i,k1)) duplicate = .true.
        end do
        if (.not. duplicate) then
          num_nb14 = num_nb14 + 1
          enefunc%num_nb14_calc(k1) = num_nb14
          nb14_calc_list(num_nb14,k1) = k4
          nb14_qq_scale (num_nb14,k1) = enefunc%fudge_qq
          nb14_lj_scale (num_nb14,k1) = enefunc%fudge_lj
          num_nb14_total = num_nb14_total + 1
        end if
      else
        num_nb14 = enefunc%num_nb14_calc(k4)
        duplicate = .false.
        do i = 1, enefunc%table%num_nonb_excl(k4)
          if (k1 == nonb_excl_list(i,k4)) duplicate = .true.
        end do
        do i = 1, num_nb14
          if (k1 == nb14_calc_list(i,k4)) duplicate = .true.
        end do
        if (.not. duplicate) then
          num_nb14 = num_nb14 + 1
          enefunc%num_nb14_calc(k4) = num_nb14
          nb14_calc_list(num_nb14,k4) = k1
          nb14_qq_scale (num_nb14,k4) = enefunc%fudge_qq
          nb14_lj_scale (num_nb14,k4) = enefunc%fudge_lj
          num_nb14_total = num_nb14_total + 1
        end if
      end if
    end do
    if (allocated(temporary_list)) &
       deallocate(temporary_list)


    ! exclude exclusions
    !
    do k = 1, size(excls(1,:))
      k1 = excls(1,k)
      k2 = excls(2,k)
      if (k1 < k2) then
        duplicate = .false.
        do i = 1, enefunc%table%num_nonb_excl(k1)
          if (k2 == nonb_excl_list(i,k1)) duplicate = .true.
        end do

        if (.not. duplicate) then
          num_excl = enefunc%table%num_nonb_excl(k1) + 1
          enefunc%table%num_nonb_excl(k1) = num_excl
          nonb_excl_list(num_excl,k1) = k2
          num_excl_total = num_excl_total + 1
        end if
      else
        duplicate = .false.
        do i = 1, enefunc%table%num_nonb_excl(k2)
          if (k1 == nonb_excl_list(i,k2)) duplicate = .true.
        end do

        if (.not. duplicate) then
          num_excl = enefunc%table%num_nonb_excl(k2) + 1
          enefunc%table%num_nonb_excl(k2) = num_excl
          nonb_excl_list(num_excl,k2) = k1
          num_excl_total = num_excl_total + 1
        end if
      end if
    end do


    ! pack 2D-array into 1D-array
    !
    allocate(enefunc%nonb_excl_list(num_excl_total), stat = alloc_stat)
    if (alloc_stat /=0) call error_msg_alloc

    call pack_array_i(natom, enefunc%table%num_nonb_excl,  &
                      nonb_excl_list, enefunc%nonb_excl_list)
    enefunc%num_nonb_excl(1:natom) = enefunc%table%num_nonb_excl(1:natom)

    max_nonb_excl_num = max(1,maxval(enefunc%table%num_nonb_excl(1:natom)))
    allocate(enefunc%table%nonb_excl_list(max_nonb_excl_num,natom), &
             stat = alloc_stat)
    if (alloc_stat /=0) call error_msg_alloc

    enefunc%table%nonb_excl_list(1:max_nonb_excl_num,1:natom) = 0
    do i = 1, natom
      enefunc%table%nonb_excl_list(1:enefunc%table%num_nonb_excl(i),i) = &
             nonb_excl_list(1:enefunc%table%num_nonb_excl(i),i)
    end do
    deallocate(nonb_excl_list, stat = dealloc_stat)


    max_nb14_num = max(1,maxval(enefunc%num_nb14_calc(1:natom)))
    allocate(enefunc%nb14_calc_list(max_nb14_num,natom), &
             enefunc%nb14_qq_scale (max_nb14_num,natom), &
             enefunc%nb14_lj_scale (max_nb14_num,natom), &
             stat = alloc_stat)
      if (alloc_stat /=0) call error_msg_alloc

    enefunc%nb14_calc_list(1:max_nb14_num,1:natom) = 0
    enefunc%nb14_qq_scale (1:max_nb14_num,1:natom) = 0.0_wp
    enefunc%nb14_lj_scale (1:max_nb14_num,1:natom) = 0.0_wp
    do i = 1, natom
      enefunc%nb14_calc_list(1:enefunc%num_nb14_calc(i),i) = &
              nb14_calc_list(1:enefunc%num_nb14_calc(i),i)
      enefunc%nb14_qq_scale (1:enefunc%num_nb14_calc(i),i) = &
              nb14_qq_scale (1:enefunc%num_nb14_calc(i),i)
      enefunc%nb14_lj_scale (1:enefunc%num_nb14_calc(i),i) = &
              nb14_lj_scale (1:enefunc%num_nb14_calc(i),i)
    end do

    deallocate(nb14_calc_list, &
               nb14_qq_scale,  &
               nb14_lj_scale, stat = dealloc_stat)

    return

  end subroutine count_nonb_excl_go

end module at_enefunc_gromacs_mod
