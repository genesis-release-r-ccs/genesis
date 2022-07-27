!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_enefunc_amber_mod
!> @brief   define potential energy functions
!! @authors Norio Takase (NT), Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_enefunc_amber_mod

  use at_enefunc_restraints_mod
  use at_enefunc_table_mod
  use at_enefunc_gbsa_mod
  use at_enefunc_pme_mod
  use at_energy_mod
  use at_boundary_str_mod
  use at_restraints_str_mod
  use at_enefunc_str_mod
  use at_energy_str_mod
  use math_libs_mod
  use dihedral_libs_mod
  use molecules_str_mod
  use fileio_table_mod
  use fileio_prmtop_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! constants
  integer, private, parameter :: MAXWRK   = 12
  integer, private, parameter :: MAX_EXCL = 16 ! = 4(1-2) + 4x3 (1-3)
  integer, private, parameter :: MAX_NB14 = 36 ! = 3x4x3 (1-4) 

  ! subroutines
  public  :: define_enefunc_amber
  private :: setup_enefunc_bond
  private :: setup_enefunc_angl
  private :: setup_enefunc_dihe
  private :: setup_enefunc_impr
  private :: setup_enefunc_cmap
  private :: setup_enefunc_nonb
  private :: count_nonb_excl

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc_amber
  !> @brief        a driver subroutine for defining potential energy
  !! @authors      NT
  !! @param[in]    ene_info   : ENERGY section control parameters information
  !! @param[in]    boundary   : boundary condition
  !! @param[in]    prmtop     : AMBER parameter topology information
  !! @param[in]    table      : lookup table information
  !! @param[in]    molecule   : molecule information
  !! @param[in]    restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc_amber(ene_info, boundary, prmtop, table, &
                                  molecule, restraints, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info 
    type(s_boundary),        intent(in)    :: boundary
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_table),           intent(in)    :: table
    type(s_molecule),        intent(in)    :: molecule
    type(s_restraints),      intent(in)    :: restraints
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: alloc_stat
    integer                  :: i, found, found2
    integer                  :: nwork


    ! bond
    !
    call setup_enefunc_bond(prmtop, molecule, enefunc)

    ! angle
    !
    call setup_enefunc_angl(prmtop, molecule, enefunc)

    ! dihedral
    !
    call setup_enefunc_dihe(prmtop, molecule, enefunc)

    ! improper
    !
    call setup_enefunc_impr(prmtop, molecule, enefunc)

    ! cmap
    !
    call setup_enefunc_cmap(ene_info, prmtop, molecule, enefunc)

    ! nonbonded
    !
    call setup_enefunc_nonb(ene_info, prmtop, molecule, enefunc)

    ! lookup table
    !
    call setup_enefunc_table(ene_info, table, molecule, enefunc)

    ! PME
    !
    call define_enefunc_pme(ene_info, boundary, molecule, enefunc)

    ! Implicit solvent
    !
    call setup_enefunc_implicit_solvent_amber(ene_info, boundary, &
                                       molecule, prmtop, enefunc)

    ! restraints
    !
    call setup_enefunc_restraints(molecule, restraints, enefunc)

    ! allocate working array
    !
    nwork = max(molecule%num_atoms, &
                enefunc%num_bonds,  &
                enefunc%num_angles, &
                enefunc%num_dihedrals, &
                enefunc%num_impropers, &
                enefunc%num_cmaps*2,   &
                enefunc%num_restraintfuncs, &
                enefunc%num_restraintgroups,&
                enefunc%num_morph_sc)

    allocate(enefunc%work(MAXWRK,nwork), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    ! write summary of energy function
    !
    if (main_rank) then

      if (.not. ene_info%table) then

        found  = 0
        found2 = 0

        do i = 1, molecule%num_atoms
          found  = found  + enefunc%num_nonb_excl(i)
          found2 = found2 + enefunc%num_nb14_calc(i)
        end do

      end if

      write(MsgOut,'(A)') &
           'Define_Enefunc_AMBER> Number of Interactions in Each Term'
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           '  bond_ene        = ', enefunc%num_bonds,           &
           '  angle_ene       = ', enefunc%num_angles
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           '  torsion_ene     = ', enefunc%num_dihedrals,       &
           '  improper_ene    = ', enefunc%num_impropers
      write(MsgOut,'(A20,I10)')                                 &
           '  cmap_ene        = ', enefunc%num_cmaps
      if (.not. ene_info%table) then
        write(MsgOut,'(A20,I10,A20,I10)')                       &
           '  nb_exclusions   = ', found,                       &
           '  nb14_calc       = ', found2
      end if
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           '  restraint_groups= ', enefunc%num_restraintgroups, &
           '  restraint_funcs = ', enefunc%num_restraintfuncs
      write(MsgOut,'(A)') ' '

    end if

    return

  end subroutine define_enefunc_amber

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_bond
  !> @brief        define BOND term in potential energy function
  !! @authors      NT
  !! @param[in]    prmtop   : AMBER parameter topology information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond(prmtop, molecule, enefunc)

    ! formal arguments
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: nbond, i, j, icnt
    integer                  :: istart, iend
    integer,          allocatable :: bond_list(:,:)
    real(wp),         allocatable :: bond_force_const(:)
    real(wp),         allocatable :: bond_dist_min(:)
    logical,          allocatable :: bond_sele(:)


    nbond = prmtop%num_bondh + prmtop%num_mbonda

    call alloc_enefunc(enefunc, EneFuncBond, nbond)

    if (enefunc%qmmm%do_qmmm) then

      allocate(bond_list(2,nbond), bond_force_const(nbond), &
               bond_dist_min(nbond), bond_sele(nbond))

      ! Retreive the original data
      icnt = 0
      do i = 1, prmtop%num_bondh

        icnt = icnt + 1

        bond_list(1,icnt) = prmtop%bond_inc_hy(1,i) / 3 + 1
        bond_list(2,icnt) = prmtop%bond_inc_hy(2,i) / 3 + 1

        bond_force_const(icnt) = &
          prmtop%bond_fcons_uniq(prmtop%bond_inc_hy(3,i))
        bond_dist_min(icnt)    = &
          prmtop%bond_equil_uniq(prmtop%bond_inc_hy(3,i))
      end do

      do i = 1, prmtop%num_mbonda

        icnt = icnt + 1

        bond_list(1,icnt) = prmtop%bond_wo_hy(1,i) / 3 + 1
        bond_list(2,icnt) = prmtop%bond_wo_hy(2,i) / 3 + 1

        bond_force_const(icnt) = &
          prmtop%bond_fcons_uniq(prmtop%bond_wo_hy(3,i))
        bond_dist_min(icnt)    = &
          prmtop%bond_equil_uniq(prmtop%bond_wo_hy(3,i))

      end do

      ! Now copy the data to enefunc
      bond_sele = .false.
      do i = 1, nbond
        enefunc%bond_list(:,i) = molecule%bond_list(:,i)

        do j = 1, nbond
          if (bond_sele(j)) cycle
          if (bond_list(1,j) == molecule%bond_list(1,i) .and. &
              bond_list(2,j) == molecule%bond_list(2,i)) then
            enefunc%bond_force_const(i) = bond_force_const(j)
            enefunc%bond_dist_min(i)    = bond_dist_min(j)
            bond_sele(j) = .true.
            exit
          end if
        end do

      end do

      deallocate(bond_list, bond_force_const, bond_dist_min, bond_sele)

      ! Remove bond terms of QM region
      enefunc%num_bonds = nbond - enefunc%qmmm%qm_nbonds

    else
      icnt = 0
      do i = 1, prmtop%num_bondh

        icnt = icnt + 1

        enefunc%bond_list(1,icnt) = prmtop%bond_inc_hy(1,i) / 3 + 1
        enefunc%bond_list(2,icnt) = prmtop%bond_inc_hy(2,i) / 3 + 1

        enefunc%bond_force_const(icnt) = &
          prmtop%bond_fcons_uniq(prmtop%bond_inc_hy(3,i))
        enefunc%bond_dist_min(icnt)    = &
          prmtop%bond_equil_uniq(prmtop%bond_inc_hy(3,i))
      end do

      do i = 1, prmtop%num_mbonda

        icnt = icnt + 1

        enefunc%bond_list(1,icnt) = prmtop%bond_wo_hy(1,i) / 3 + 1
        enefunc%bond_list(2,icnt) = prmtop%bond_wo_hy(2,i) / 3 + 1

        enefunc%bond_force_const(icnt) = &
          prmtop%bond_fcons_uniq(prmtop%bond_wo_hy(3,i))
        enefunc%bond_dist_min(icnt)    = &
          prmtop%bond_equil_uniq(prmtop%bond_wo_hy(3,i))

      end do

      enefunc%num_bonds = nbond
    end if


    call get_loop_index(enefunc%num_bonds, istart, iend)

    enefunc%istart_bond = istart
    enefunc%iend_bond   = iend

    return

  end subroutine setup_enefunc_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_angl
  !> @brief        define ANGLE term in potential energy function
  !! @authors      NT, KY
  !! @param[in]    prmtop   : AMBER parameter topology information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_angl(prmtop, molecule, enefunc)

    ! formal arguments
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: nangl, i, j, icnt
    integer                  :: istart, iend
    integer,          allocatable :: angl_list(:,:)
    real(wp),         allocatable :: angl_force_const(:)
    real(wp),         allocatable :: angl_theta_min(:)
    logical,          allocatable :: angl_sele(:)


    nangl = prmtop%num_anglh + prmtop%num_mangla

    call alloc_enefunc(enefunc, EneFuncAngl, nangl)

    if (enefunc%qmmm%do_qmmm) then

      allocate(angl_list(3,nangl), angl_force_const(nangl), &
               angl_theta_min(nangl), angl_sele(nangl))

      ! Retreive the original data
      icnt = 0
      do i = 1, prmtop%num_anglh

        icnt = icnt + 1

        angl_list(1,icnt) = prmtop%angl_inc_hy(1,i) / 3 + 1
        angl_list(2,icnt) = prmtop%angl_inc_hy(2,i) / 3 + 1
        angl_list(3,icnt) = prmtop%angl_inc_hy(3,i) / 3 + 1

        angl_force_const(icnt) = &
          prmtop%angl_fcons_uniq(prmtop%angl_inc_hy(4,i))
        angl_theta_min(icnt)   = &
          prmtop%angl_equil_uniq(prmtop%angl_inc_hy(4,i))

      end do

      do i = 1, prmtop%num_mangla

        icnt = icnt + 1

        angl_list(1,icnt) = prmtop%angl_wo_hy(1,i) / 3 + 1
        angl_list(2,icnt) = prmtop%angl_wo_hy(2,i) / 3 + 1
        angl_list(3,icnt) = prmtop%angl_wo_hy(3,i) / 3 + 1

        angl_force_const(icnt) = &
          prmtop%angl_fcons_uniq(prmtop%angl_wo_hy(4,i))
        angl_theta_min(icnt)   = &
          prmtop%angl_equil_uniq(prmtop%angl_wo_hy(4,i))

      end do

      ! Now copy the data to enefunc
      angl_sele = .false.
      do i = 1, nangl
        enefunc%angl_list(:,i) = molecule%angl_list(:,i)

        do j = 1, nangl
          if (angl_sele(j)) cycle
          if (angl_list(1,j) == molecule%angl_list(1,i) .and. &
              angl_list(2,j) == molecule%angl_list(2,i) .and. &
              angl_list(3,j) == molecule%angl_list(3,i)) then
            enefunc%angl_force_const(i) = angl_force_const(j)
            enefunc%angl_theta_min(i)   = angl_theta_min(j)
            angl_sele(j) = .true.
            exit
          end if
        end do

      end do

      deallocate(angl_list, angl_force_const, angl_theta_min, angl_sele)

      ! Remove angl terms of QM region
      enefunc%num_angles = nangl - enefunc%qmmm%qm_nangles

    else
      icnt = 0
      do i = 1, prmtop%num_anglh

        icnt = icnt + 1

        enefunc%angl_list(1,icnt) = prmtop%angl_inc_hy(1,i) / 3 + 1
        enefunc%angl_list(2,icnt) = prmtop%angl_inc_hy(2,i) / 3 + 1
        enefunc%angl_list(3,icnt) = prmtop%angl_inc_hy(3,i) / 3 + 1

        enefunc%angl_force_const(icnt) = &
          prmtop%angl_fcons_uniq(prmtop%angl_inc_hy(4,i))
        enefunc%angl_theta_min(icnt)   = &
          prmtop%angl_equil_uniq(prmtop%angl_inc_hy(4,i))

      end do

      do i = 1, prmtop%num_mangla

        icnt = icnt + 1

        enefunc%angl_list(1,icnt) = prmtop%angl_wo_hy(1,i) / 3 + 1
        enefunc%angl_list(2,icnt) = prmtop%angl_wo_hy(2,i) / 3 + 1
        enefunc%angl_list(3,icnt) = prmtop%angl_wo_hy(3,i) / 3 + 1

        enefunc%angl_force_const(icnt) = &
          prmtop%angl_fcons_uniq(prmtop%angl_wo_hy(4,i))
        enefunc%angl_theta_min(icnt)   = &
          prmtop%angl_equil_uniq(prmtop%angl_wo_hy(4,i))

      end do

      enefunc%num_angles = nangl
    end if


    call get_loop_index(enefunc%num_angles, istart, iend)

    enefunc%istart_angle = istart
    enefunc%iend_angle   = iend

    return

  end subroutine setup_enefunc_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_dihe
  !> @brief        define DIHEDRAL term in potential energy function
  !! @authors      NT, KY
  !! @param[in]    prmtop   : AMBER parameter topology information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_dihe(prmtop, molecule, enefunc)

    ! formal arguments
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, idih
    integer                  :: istart, iend
    integer,     allocatable :: dihe_list(:,:)
    real(wp),    allocatable :: dihe_force_const(:)
    integer,     allocatable :: dihe_periodicity(:)
    real(wp),    allocatable :: dihe_phase(:)
    real(wp),    allocatable :: dihe_scee(:)
    real(wp),    allocatable :: dihe_scnb(:)
    logical,     allocatable :: dihe_sele(:)


    enefunc%num_dihedrals = 0

    do i = 1, prmtop%num_diheh
      if (prmtop%dihe_inc_hy(4,i) >= 0) then
        enefunc%num_dihedrals = enefunc%num_dihedrals + 1
      end if
    end do

    do i = 1, prmtop%num_mdihea
      if (prmtop%dihe_wo_hy(4,i) >= 0) then
        enefunc%num_dihedrals = enefunc%num_dihedrals + 1
      end if
    end do

    call alloc_enefunc(enefunc, EneFuncDihe, enefunc%num_dihedrals)

    if (enefunc%qmmm%do_qmmm) then

      allocate(dihe_list(4,enefunc%num_dihedrals),      &
               dihe_force_const(enefunc%num_dihedrals), &
               dihe_periodicity(enefunc%num_dihedrals), &
               dihe_phase(enefunc%num_dihedrals),       &
               dihe_scee(enefunc%num_dihedrals),        &
               dihe_scnb(enefunc%num_dihedrals),        &
               dihe_sele(enefunc%num_dihedrals))

      ! Retreive the original data
      idih = 0

      do i = 1, prmtop%num_diheh
      
        if (prmtop%dihe_inc_hy(4,i) >= 0) then

          idih = idih + 1
          dihe_list(1,idih) =      prmtop%dihe_inc_hy(1,i)  / 3 + 1
          dihe_list(2,idih) =      prmtop%dihe_inc_hy(2,i)  / 3 + 1
          dihe_list(3,idih) = iabs(prmtop%dihe_inc_hy(3,i)) / 3 + 1
          dihe_list(4,idih) =      prmtop%dihe_inc_hy(4,i)  / 3 + 1

          dihe_force_const(idih) = &
            prmtop%dihe_fcons_uniq(prmtop%dihe_inc_hy(5,i))
          dihe_periodicity(idih) = &
            prmtop%dihe_perio_uniq(prmtop%dihe_inc_hy(5,i))
          dihe_phase(idih)       = &
            prmtop%dihe_phase_uniq(prmtop%dihe_inc_hy(5,i))

          if (prmtop%lscee_scale_factor) then
            dihe_scee(idih) = &
              1.0_wp / prmtop%scee_scale_fact(prmtop%dihe_inc_hy(5,i))
          else
            dihe_scee(idih) = 1.0_wp / 1.2_wp
          end if

          if (prmtop%lscnb_scale_factor) then
            dihe_scnb(idih) = &
              1.0_wp / prmtop%scnb_scale_fact(prmtop%dihe_inc_hy(5,i))
          else
            dihe_scnb(idih) = 1.0_wp / 2.0_wp
          end if

        end if

      end do

      do i = 1, prmtop%num_mdihea

        if (prmtop%dihe_wo_hy(4,i) >= 0) then

          idih = idih + 1
          dihe_list(1,idih) =      prmtop%dihe_wo_hy(1,i)  / 3 + 1
          dihe_list(2,idih) =      prmtop%dihe_wo_hy(2,i)  / 3 + 1
          dihe_list(3,idih) = iabs(prmtop%dihe_wo_hy(3,i)) / 3 + 1
          dihe_list(4,idih) =      prmtop%dihe_wo_hy(4,i)  / 3 + 1

          dihe_force_const(idih) = &
            prmtop%dihe_fcons_uniq(prmtop%dihe_wo_hy(5,i))
          dihe_periodicity(idih) = &
            prmtop%dihe_perio_uniq(prmtop%dihe_wo_hy(5,i))
          dihe_phase(idih)       = &
            prmtop%dihe_phase_uniq(prmtop%dihe_wo_hy(5,i))

          if (prmtop%lscee_scale_factor) then
            dihe_scee(idih) = &
              1.0_wp / prmtop%scee_scale_fact(prmtop%dihe_wo_hy(5,i))
          else
            dihe_scee(idih) = 1.0_wp / 1.2_wp
          end if

          if (prmtop%lscnb_scale_factor) then
            dihe_scnb(idih) = &
              1.0_wp / prmtop%scnb_scale_fact(prmtop%dihe_wo_hy(5,i))
          else
            dihe_scnb(idih) = 1.0_wp / 2.0_wp
          end if

        end if

      end do

      ! Now copy the data to enefunc
      dihe_sele = .false.
      do i = 1, enefunc%num_dihedrals
        enefunc%dihe_list(:,i) = molecule%dihe_list(:,i)

        do j = 1, enefunc%num_dihedrals
          if (dihe_sele(j)) cycle
          if (dihe_list(1,j) == molecule%dihe_list(1,i) .and. &
              dihe_list(2,j) == molecule%dihe_list(2,i) .and. &
              dihe_list(3,j) == molecule%dihe_list(3,i) .and. &
              dihe_list(4,j) == molecule%dihe_list(4,i)) then
            enefunc%dihe_force_const(i) = dihe_force_const(j)
            enefunc%dihe_periodicity(i) = dihe_periodicity(j)
            enefunc%dihe_phase(i)       = dihe_phase(j)
            enefunc%dihe_scee(i)        = dihe_scee(j)
            enefunc%dihe_scnb(i)        = dihe_scnb(j)
            dihe_sele(j) = .true.
            exit
          end if
        end do

      end do

      deallocate(dihe_list, dihe_force_const, dihe_periodicity, &
                 dihe_phase, dihe_scee, dihe_scnb, dihe_sele)

      ! Remove dihe terms of QM region
      enefunc%num_dihedrals = molecule%num_dihedrals - enefunc%qmmm%qm_ndihedrals

    else
      idih = 0

      do i = 1, prmtop%num_diheh
      
        if (prmtop%dihe_inc_hy(4,i) >= 0) then

          idih = idih + 1
          enefunc%dihe_list(1,idih) =      prmtop%dihe_inc_hy(1,i)  / 3 + 1
          enefunc%dihe_list(2,idih) =      prmtop%dihe_inc_hy(2,i)  / 3 + 1
          enefunc%dihe_list(3,idih) = iabs(prmtop%dihe_inc_hy(3,i)) / 3 + 1
          enefunc%dihe_list(4,idih) =      prmtop%dihe_inc_hy(4,i)  / 3 + 1

          enefunc%dihe_force_const(idih) = &
            prmtop%dihe_fcons_uniq(prmtop%dihe_inc_hy(5,i))
          enefunc%dihe_periodicity(idih) = &
            prmtop%dihe_perio_uniq(prmtop%dihe_inc_hy(5,i))
          enefunc%dihe_phase(idih)       = &
            prmtop%dihe_phase_uniq(prmtop%dihe_inc_hy(5,i))

          if (prmtop%lscee_scale_factor) then
            enefunc%dihe_scee(idih) = &
              1.0_wp / prmtop%scee_scale_fact(prmtop%dihe_inc_hy(5,i))
          else
            enefunc%dihe_scee(idih) = 1.0_wp / 1.2_wp
          end if

          if (prmtop%lscnb_scale_factor) then
            enefunc%dihe_scnb(idih) = &
              1.0_wp / prmtop%scnb_scale_fact(prmtop%dihe_inc_hy(5,i))
          else
            enefunc%dihe_scnb(idih) = 1.0_wp / 2.0_wp
          end if

        end if

      end do

      ! check dihe_wo_hy array
      !
      do i = 1, prmtop%num_mdihea

        if (prmtop%dihe_wo_hy(4,i) >= 0) then

          idih = idih + 1
          enefunc%dihe_list(1,idih) =      prmtop%dihe_wo_hy(1,i)  / 3 + 1
          enefunc%dihe_list(2,idih) =      prmtop%dihe_wo_hy(2,i)  / 3 + 1
          enefunc%dihe_list(3,idih) = iabs(prmtop%dihe_wo_hy(3,i)) / 3 + 1
          enefunc%dihe_list(4,idih) =      prmtop%dihe_wo_hy(4,i)  / 3 + 1

          enefunc%dihe_force_const(idih) = &
            prmtop%dihe_fcons_uniq(prmtop%dihe_wo_hy(5,i))
          enefunc%dihe_periodicity(idih) = &
            prmtop%dihe_perio_uniq(prmtop%dihe_wo_hy(5,i))
          enefunc%dihe_phase(idih)       = &
            prmtop%dihe_phase_uniq(prmtop%dihe_wo_hy(5,i))

          if (prmtop%lscee_scale_factor) then
            enefunc%dihe_scee(idih) = &
              1.0_wp / prmtop%scee_scale_fact(prmtop%dihe_wo_hy(5,i))
          else
            enefunc%dihe_scee(idih) = 1.0_wp / 1.2_wp
          end if

          if (prmtop%lscnb_scale_factor) then
            enefunc%dihe_scnb(idih) = &
              1.0_wp / prmtop%scnb_scale_fact(prmtop%dihe_wo_hy(5,i))
          else
            enefunc%dihe_scnb(idih) = 1.0_wp / 2.0_wp
          end if

        end if

      end do
    end if


    call get_loop_index(enefunc%num_dihedrals, istart, iend)

    enefunc%istart_dihedral = istart
    enefunc%iend_dihedral   = iend

    return

  end subroutine setup_enefunc_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_impr
  !> @brief        define IMPROPER dihedral term in potential energy function
  !! @authors      NT, KY
  !! @param[in]    prmtop   : AMBER parameter topology information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_impr(prmtop,  molecule, enefunc)

    ! formal arguments
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, iimp
    integer                  :: istart, iend
    integer,     allocatable :: impr_list(:,:)
    real(wp),    allocatable :: impr_force_const(:)
    integer,     allocatable :: impr_periodicity(:)
    real(wp),    allocatable :: impr_phase(:)
    logical,     allocatable :: impr_sele(:)


    enefunc%num_impropers = 0

    do i = 1, prmtop%num_diheh
      if (prmtop%dihe_inc_hy(4,i) < 0) then
        enefunc%num_impropers = enefunc%num_impropers + 1
      end if
    end do

    do i = 1, prmtop%num_mdihea
      if (prmtop%dihe_wo_hy(4,i) < 0) then
        enefunc%num_impropers = enefunc%num_impropers + 1
      end if
    end do

    call alloc_enefunc(enefunc, EneFuncImpr, enefunc%num_impropers)

    if (enefunc%qmmm%do_qmmm) then

      allocate(impr_list(4,enefunc%num_impropers), &
               impr_force_const(enefunc%num_impropers), &
               impr_periodicity(enefunc%num_impropers), &
               impr_phase(enefunc%num_impropers),       &
               impr_sele(enefunc%num_impropers))

      ! Retreive the original data
      iimp = 0
      do i = 1, prmtop%num_diheh
        if (prmtop%dihe_inc_hy(4,i) < 0) then
          iimp = iimp + 1

          impr_list(1,iimp) =      prmtop%dihe_inc_hy(1,i)  / 3 + 1
          impr_list(2,iimp) =      prmtop%dihe_inc_hy(2,i)  / 3 + 1
          impr_list(3,iimp) = iabs(prmtop%dihe_inc_hy(3,i)) / 3 + 1
          impr_list(4,iimp) = iabs(prmtop%dihe_inc_hy(4,i)) / 3 + 1

          impr_force_const(iimp) = &
            prmtop%dihe_fcons_uniq(prmtop%dihe_inc_hy(5,i))
          impr_periodicity(iimp) = &
            prmtop%dihe_perio_uniq(prmtop%dihe_inc_hy(5,i))
          impr_phase(iimp)       = &
            prmtop%dihe_phase_uniq(prmtop%dihe_inc_hy(5,i))
      
        end if
      end do

      do i = 1, prmtop%num_mdihea
        if (prmtop%dihe_wo_hy(4,i) < 0) then
          iimp = iimp + 1

          impr_list(1,iimp) =      prmtop%dihe_wo_hy(1,i)  / 3 + 1
          impr_list(2,iimp) =      prmtop%dihe_wo_hy(2,i)  / 3 + 1
          impr_list(3,iimp) = iabs(prmtop%dihe_wo_hy(3,i)) / 3 + 1
          impr_list(4,iimp) = iabs(prmtop%dihe_wo_hy(4,i)) / 3 + 1

          impr_force_const(iimp) = &
            prmtop%dihe_fcons_uniq(prmtop%dihe_wo_hy(5,i))
          impr_periodicity(iimp) = &
            prmtop%dihe_perio_uniq(prmtop%dihe_wo_hy(5,i))
          impr_phase(iimp)       = &
            prmtop%dihe_phase_uniq(prmtop%dihe_wo_hy(5,i))

        end if
      end do

      ! Now copy the data to enefunc
      impr_sele = .false.
      do i = 1, enefunc%num_impropers
        enefunc%impr_list(:,i) = molecule%impr_list(:,i)

        do j = 1, enefunc%num_impropers
          if (impr_sele(j)) cycle
          if (impr_list(1,j) == molecule%impr_list(1,i) .and. &
              impr_list(2,j) == molecule%impr_list(2,i) .and. &
              impr_list(3,j) == molecule%impr_list(3,i) .and. &
              impr_list(4,j) == molecule%impr_list(4,i)) then
            enefunc%impr_force_const(i) = impr_force_const(j)
            enefunc%impr_periodicity(i) = impr_periodicity(j)
            enefunc%impr_phase(i)       = impr_phase(j)
            impr_sele(j) = .true.
            exit
          end if
        end do

      end do

      deallocate(impr_list, impr_force_const, impr_periodicity, impr_phase, impr_sele)

      ! Remove impr terms of QM region
      enefunc%num_impropers = molecule%num_impropers - enefunc%qmmm%qm_nimpropers

    else
      iimp = 0

      do i = 1, prmtop%num_diheh
      
        if (prmtop%dihe_inc_hy(4,i) < 0) then

          iimp = iimp + 1
          enefunc%impr_list(1,iimp) =      prmtop%dihe_inc_hy(1,i)  / 3 + 1
          enefunc%impr_list(2,iimp) =      prmtop%dihe_inc_hy(2,i)  / 3 + 1
          enefunc%impr_list(3,iimp) = iabs(prmtop%dihe_inc_hy(3,i)) / 3 + 1
          enefunc%impr_list(4,iimp) = iabs(prmtop%dihe_inc_hy(4,i)) / 3 + 1

          enefunc%impr_force_const(iimp) = &
            prmtop%dihe_fcons_uniq(prmtop%dihe_inc_hy(5,i))
          enefunc%impr_periodicity(iimp) = &
            prmtop%dihe_perio_uniq(prmtop%dihe_inc_hy(5,i))
          enefunc%impr_phase(iimp)       = &
            prmtop%dihe_phase_uniq(prmtop%dihe_inc_hy(5,i))

        end if

      end do

      ! check dihe_wo_hy array
      !
      do i = 1, prmtop%num_mdihea

        if (prmtop%dihe_wo_hy(4,i) < 0) then

          iimp = iimp + 1
          enefunc%impr_list(1,iimp) =      prmtop%dihe_wo_hy(1,i)  / 3 + 1
          enefunc%impr_list(2,iimp) =      prmtop%dihe_wo_hy(2,i)  / 3 + 1
          enefunc%impr_list(3,iimp) = iabs(prmtop%dihe_wo_hy(3,i)) / 3 + 1
          enefunc%impr_list(4,iimp) = iabs(prmtop%dihe_wo_hy(4,i)) / 3 + 1

          enefunc%impr_force_const(iimp) = &
            prmtop%dihe_fcons_uniq(prmtop%dihe_wo_hy(5,i))
          enefunc%impr_periodicity(iimp) = &
            prmtop%dihe_perio_uniq(prmtop%dihe_wo_hy(5,i))
          enefunc%impr_phase(iimp)       = &
            prmtop%dihe_phase_uniq(prmtop%dihe_wo_hy(5,i))

        end if

      end do

    end if


    call get_loop_index(enefunc%num_impropers, istart, iend)

    enefunc%istart_improper = istart
    enefunc%iend_improper   = iend

    return

  end subroutine setup_enefunc_impr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cmap
  !> @brief        define cmap term in potential energy function
  !! @authors      JJ, KY
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[in]    prmtop   : AMBER PRMTOP information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !!
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cmap(ene_info, prmtop, molecule, enefunc)

    ! formal variables
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k, idih1, idih2, ityp
    integer                  :: ncmap, ncmap_p, found, ngrid0
    integer                  :: alloc_stat, dealloc_stat
    integer                  :: istart, iend
    character(6)             :: ci1, ci2, ci3, ci4, ci5, ci6, ci7, ci8
    logical                  :: periodic

    real(wp),    allocatable :: c_ij(:,:,:,:)      ! cmap coeffs
    logical,     allocatable :: cmap_sele(:)



    periodic = ene_info%cmap_pspline
    ncmap    = molecule%num_cmaps
    ncmap_p  = prmtop%num_cmaptype
    enefunc%num_cmaps = ncmap

    ngrid0 = 0
    do i = 1, ncmap_p
      ngrid0 = max(ngrid0, prmtop%cmap_resolution(i))
    end do

    call alloc_enefunc(enefunc, EneFuncCmap, ncmap, ngrid0)
    call alloc_enefunc(enefunc, EneFuncCmapType, ncmap_p, ngrid0)

    if(enefunc%qmmm%do_qmmm) then

      allocate(cmap_sele(ncmap))
      cmap_sele = .false.

      ! Copy the data to enefunc in a rearranged order
      do i = 1, ncmap
        enefunc%cmap_list(1:8,i) = molecule%cmap_list(1:8,i)
        do j = 1, prmtop%num_cmap
          if (cmap_sele(j)) cycle
          if (prmtop%cmap_list(1,j) == molecule%cmap_list(1,i) .and. &
              prmtop%cmap_list(2,j) == molecule%cmap_list(2,i) .and. &
              prmtop%cmap_list(3,j) == molecule%cmap_list(3,i) .and. &
              prmtop%cmap_list(4,j) == molecule%cmap_list(4,i) .and. &
              prmtop%cmap_list(2,j) == molecule%cmap_list(5,i) .and. &
              prmtop%cmap_list(3,j) == molecule%cmap_list(6,i) .and. &
              prmtop%cmap_list(4,j) == molecule%cmap_list(7,i) .and. &
              prmtop%cmap_list(5,j) == molecule%cmap_list(8,i)) then

            enefunc%cmap_type(i) = prmtop%cmap_list(6,j)
            cmap_sele(j) = .true.
            exit

          end if
        end do
      end do

      ! Remove CMAP terms of QM region
      enefunc%num_cmaps = molecule%num_cmaps - enefunc%qmmm%qm_ncmaps

    else

      do i = 1, ncmap
    
        enefunc%cmap_list(1:8,i) = molecule%cmap_list(1:8,i)
        enefunc%cmap_type(    i) = prmtop%cmap_list(6,i)

      end do

    endif

    ! allocate local arrays
    !
    alloc_stat = 0
    allocate(c_ij(4,4,ngrid0,ngrid0), stat = alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    !  derive cmap coefficients by bicubic interpolation
    !
    do ityp = 1, ncmap_p
      enefunc%cmap_resolution(ityp) = prmtop%cmap_resolution(ityp)
    end do

    do ityp = 1, ncmap_p
      if (periodic) then
        call derive_cmap_amber_coefficients_p (ityp, prmtop, c_ij)
      else
        call derive_cmap_amber_coefficients_np(ityp, prmtop, c_ij)
      end if

      ! copy values of C_ij to enefunc%cmap_coef for output
      !
      do idih2 = 1, ngrid0
        do idih1 = 1, ngrid0
          do j = 1, 4
            do i = 1, 4
              enefunc%cmap_coef(i,j,idih1,idih2,ityp)     &
                = c_ij(i,j,idih1,idih2)
            end do
          end do
        end do
      end do
    end do

    !  deallocate local arrays
    !
    deallocate(c_ij, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc

    ! define cmap interaction list for each processor
    !
    call get_loop_index(enefunc%num_cmaps, istart, iend)
    enefunc%istart_cmap = istart
    enefunc%iend_cmap   = iend

    ! write summary
    !
    if (main_rank) then
      if (periodic) then
        write(MsgOut,'(A)')                                       &
        'Setup_Enefunc_Cmap_Par> '//                              &
        'Periodic-boundary spline is used to derive cmap coefs.'
        write(MsgOut,'(A)') ''
      else
        write(MsgOut,'(A)')                                       &
        'Setup_Enefunc_Cmap_Par> Natural spline is used'//        &
        ' to derive cmap coefs.'
        write(MsgOut,'(A)') ''
      end if
    end if

    return

  end subroutine setup_enefunc_cmap

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_nonb
  !> @brief        define NON-BOND term in potential energy function
  !! @authors      NT
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[in]    prmtop   : AMBER parameter topology information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_nonb(ene_info, prmtop, molecule, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k, nnonb, numi, numl, nb_par_idx


    nnonb = prmtop%num_types

    ELECOEF          = ELECOEF_AMBER

    call alloc_enefunc(enefunc, EneFuncNbon, nnonb)

    do i = 1, nnonb
      do j = 1, nnonb

        nb_par_idx = prmtop%nb_par_idx(nnonb*(i-1)+j)

        if (nb_par_idx < 0) then
          enefunc%nb14_lj12(i,j) = 0.0_wp
          enefunc%nb14_lj6 (i,j) = 0.0_wp
          enefunc%nonb_lj12(i,j) = 0.0_wp
          enefunc%nonb_lj6 (i,j) = 0.0_wp
        else
          enefunc%nb14_lj12(i,j) = prmtop%lennarda(nb_par_idx)
          enefunc%nb14_lj6 (i,j) = prmtop%lennardb(nb_par_idx)
          enefunc%nonb_lj12(i,j) = prmtop%lennarda(nb_par_idx)
          enefunc%nonb_lj6 (i,j) = prmtop%lennardb(nb_par_idx)
        end if

      end do
    end do

    enefunc%num_atom_cls = nnonb


    ! treatment for 1-2, 1-3, 1-4 interactions
    !
    if (.not.ene_info%table) &
      call count_nonb_excl(molecule, enefunc)

    return

  end subroutine setup_enefunc_nonb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    count_nonb_excl
  !> @brief        exclude 1-2, 1-3 interactions
  !! @authors      NT
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine count_nonb_excl(molecule, enefunc)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

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

    allocate(nonb_excl_list(MAX_EXCL,molecule%num_atoms),stat=alloc_stat)
    if (alloc_stat /=0) &
      call error_msg_alloc

    allocate(enefunc%num_nb14_calc(molecule%num_atoms), stat=alloc_stat)
    if (alloc_stat /=0) &
      call error_msg_alloc

    allocate(nb14_calc_list(MAX_NB14,molecule%num_atoms), &
             nb14_qq_scale (MAX_NB14,molecule%num_atoms), &
             nb14_lj_scale (MAX_NB14,molecule%num_atoms), stat=alloc_stat)
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
          nb14_qq_scale (num_nb14,k1) = enefunc%dihe_scee(k)
          nb14_lj_scale (num_nb14,k1) = enefunc%dihe_scnb(k)
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
          nb14_qq_scale (num_nb14,k4) = enefunc%dihe_scee(k)
          nb14_lj_scale (num_nb14,k4) = enefunc%dihe_scnb(k)
          num_nb14_total = num_nb14_total + 1
        end if

      end if
    end do
    
    ! pack 2D-array into 1D-array
    !
    allocate(enefunc%nonb_excl_list(num_excl_total), stat = alloc_stat) 
    if (alloc_stat /=0) &
      call error_msg_alloc

    call pack_array_i(molecule%num_atoms, enefunc%num_nonb_excl,  &
                      nonb_excl_list, enefunc%nonb_excl_list)

    deallocate(nonb_excl_list, stat = dealloc_stat)


    max_nb14_num = max(1,maxval(enefunc%num_nb14_calc(1:molecule%num_atoms)))
    allocate(enefunc%nb14_calc_list(max_nb14_num,molecule%num_atoms), &
             enefunc%nb14_qq_scale (max_nb14_num,molecule%num_atoms), &
             enefunc%nb14_lj_scale (max_nb14_num,molecule%num_atoms), &
             stat = alloc_stat)
    if (alloc_stat /=0) &
      call error_msg_alloc

    enefunc%nb14_calc_list(1:max_nb14_num,1:molecule%num_atoms) = 0
    enefunc%nb14_qq_scale (1:max_nb14_num,1:molecule%num_atoms) = 0.0_wp
    enefunc%nb14_lj_scale (1:max_nb14_num,1:molecule%num_atoms) = 0.0_wp

    do i = 1, molecule%num_atoms
      enefunc%nb14_calc_list(1:enefunc%num_nb14_calc(i),i) = &
                       nb14_calc_list(1:enefunc%num_nb14_calc(i),i)
      enefunc%nb14_lj_scale (1:enefunc%num_nb14_calc(i),i) = &
                       nb14_qq_scale (1:enefunc%num_nb14_calc(i),i)
      enefunc%nb14_lj_scale (1:enefunc%num_nb14_calc(i),i) = &
                       nb14_lj_scale (1:enefunc%num_nb14_calc(i),i)
    end do


    deallocate(nb14_calc_list, &
               nb14_qq_scale,  &
               nb14_lj_scale, stat = dealloc_stat)

    return

  end subroutine count_nonb_excl

end module at_enefunc_amber_mod
