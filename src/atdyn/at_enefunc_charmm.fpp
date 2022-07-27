!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_enefunc_charmm_mod
!> @brief   define potential energy functions
!! @authors Takaharu Mori (TM), Yuji Sugita (YS), Takao Yoda (TY), 
!!          Takashi Imai (TI), Jaewoon Jung (JJ), Chigusa Kobayashi (CK)
!!          Norio Takase (NT), Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_enefunc_charmm_mod

  use at_enefunc_restraints_mod
  use at_enefunc_table_mod
  use at_enefunc_pme_mod
  use at_enefunc_gbsa_mod
  use at_energy_mod
  use at_boundary_str_mod
  use at_restraints_str_mod
  use at_enefunc_str_mod
  use math_libs_mod
  use nbond_list_mod
  use dihedral_libs_mod
  use molecules_str_mod
  use fileio_table_mod
  use fileio_par_mod
  use fileio_gpr_mod
  use fileio_eef1_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! constants
  integer, parameter :: MAXWRK   = 12
  integer, parameter :: MAX_EXCL = 52 ! = 4(1-2) + 4x3 (1-3) + 3*4*3 (1-4:CH19) 
  integer, parameter :: MAX_NB14 = 36 ! = 3x4x3 (1-4) 

  ! subroutines
  public  :: define_enefunc_charmm
  private :: setup_enefunc_bond_par
  private :: setup_enefunc_angl_urey_par
  private :: setup_enefunc_dihe_par
  private :: setup_enefunc_impr_par
  private :: setup_enefunc_cmap_par  ! assign cmap parameters to each dihedral
  private :: setup_enefunc_nonb_par
  private :: count_nonb_excl_par
  private :: count_nonb_excl_par_c19

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc_charmm
  !> @brief        a driver subroutine for defining potential energy function
  !! @authors      YS, TI, JJ, CK
  !! @param[in]    ene_info   : ENERGY section control parameters information
  !! @param[in]    boundary   : boundary conditions information
  !! @param[in]    par        : CHARMM PAR information
  !! @param[in]    table      : lookup table information
  !! @param[in]    molecule   : molecule information
  !! @param[in]    restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc_charmm(ene_info, boundary, par, eelf1, table, &
                                   molecule, restraints, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info 
    type(s_boundary),        intent(in)    :: boundary
    type(s_par),             intent(in)    :: par
    type(s_eef1),            intent(in)    :: eelf1
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
    call setup_enefunc_bond_par(par, molecule, enefunc)

    ! angle
    !
    call setup_enefunc_angl_urey_par(par, molecule, enefunc)

    ! dihedral
    !
    call setup_enefunc_dihe_par(par, molecule, enefunc)

    ! improper
    !
    call setup_enefunc_impr_par(par, molecule, enefunc)

    ! cmap
    !
    call setup_enefunc_cmap_par(ene_info, par, molecule, enefunc)

    ! nonbonded
    !
    call setup_enefunc_nonb_par(ene_info, par, molecule, enefunc)

    ! lookup table
    !
    call setup_enefunc_table(ene_info, table, molecule, enefunc)

    ! PME
    !
    call define_enefunc_pme(ene_info, boundary, molecule, enefunc)

    ! Implicit solvent
    !
    call setup_enefunc_implicit_solvent_charmm(ene_info, boundary, &
                                    molecule, par, eelf1, enefunc)

    ! restraints
    !
    call setup_enefunc_restraints(molecule, restraints, enefunc)

    ! allocate working array
    !
    nwork = max(molecule%num_atoms, &
                enefunc%num_bonds,  &
                enefunc%num_angles, &
                enefunc%num_ureys,  &
                enefunc%num_dihedrals, &
                enefunc%num_impropers, &
                enefunc%num_cmaps*2,   &
                enefunc%num_restraintfuncs, &
                enefunc%num_restraintgroups)

    allocate(enefunc%work(MAXWRK,nwork),stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    ! write summary of energy function
    !
    if (main_rank) then

      if (.not. ene_info%table .or. &
          enefunc%forcefield==ForcefieldCHARMM19) then

        found  = 0
        found2 = 0

        do i = 1, molecule%num_atoms
          found  = found  + enefunc%num_nonb_excl(i)
          found2 = found2 + enefunc%num_nb14_calc(i)
        end do

      end if

      write(MsgOut,'(A)') &
           'Define_Enefunc_CHARMM> Number of Interactions in Each Term'
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           '  bond_ene        = ', enefunc%num_bonds,           &
           '  angle_ene       = ', enefunc%num_angles
      write(MsgOut,'(A20,I10)')  &
           '  urey_ene        = ', enefunc%num_ureys
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           '  torsion_ene     = ', enefunc%num_dihedrals,       &
           '  improper_ene    = ', enefunc%num_impropers
      write(MsgOut,'(A20,I10)')                                 &
           '  cmap_ene        = ', enefunc%num_cmaps
      if (.not. ene_info%table .or. enefunc%forcefield==ForcefieldCHARMM19) then
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

  end subroutine define_enefunc_charmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_bond_par
  !> @brief        define BOND term in potential energy function
  !! @authors      YS, JJ, TM
  !! @param[in]    par      : CHARMM PAR information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond_par(par, molecule, enefunc)

    ! formal arguments
    type(s_par),             intent(in)    :: par
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: nbond, nbond_p
    integer                  :: i, j, found
    integer                  :: istart, iend
    character(6)             :: ci1, ci2


    if(enefunc%qmmm%do_qmmm) then
      nbond   = molecule%num_bonds - enefunc%qmmm%qm_nbonds
    else
      nbond   = molecule%num_bonds
    end if
    nbond_p = par%num_bonds

    call alloc_enefunc(enefunc, EneFuncBond, nbond)

    found = 0
    do i = 1, nbond

      ci1 = molecule%atom_cls_name(molecule%bond_list(1, i))
      ci2 = molecule%atom_cls_name(molecule%bond_list(2, i))
      enefunc%bond_list(1:2,i)   = molecule%bond_list(1:2,i)

      do j = 1, nbond_p
        if ((ci1 == par%bond_atom_cls(1, j) .and. &
             ci2 == par%bond_atom_cls(2, j)) .or. &
            (ci1 == par%bond_atom_cls(2, j) .and. &
             ci2 == par%bond_atom_cls(1, j))) then
          enefunc%bond_force_const(i) = par%bond_force_const(j) 
          enefunc%bond_dist_min(i)    = par%bond_dist_min(j)
          found = found + 1
          exit
        end if
      end do

      if (j == nbond_p + 1) &
        write(MsgOut,*) 'Setup_Enefunc_Bond_Par> not found BOND: [', &
             ci1, ']-[', ci2, '] in parameter file. (ERROR)'

    end do
    enefunc%num_bonds = found

    if (found /= nbond) &
      call error_msg('Setup_Enefunc_Bond_Par> '//                    &
                     'Some Bond Paremeters are missing.')

    ! define bond interaction list for each processor
    !
    call get_loop_index(enefunc%num_bonds, istart, iend)
    enefunc%istart_bond = istart
    enefunc%iend_bond   = iend

    return
    
  end subroutine setup_enefunc_bond_par

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_angl_urey_par
  !> @brief        define ANGLE term in potential energy function
  !! @authors      YS, JJ, TM
  !! @param[in]    par      : CHARMM PAR information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_angl_urey_par(par, molecule, enefunc)

    ! formal arguments
    type(s_par),             intent(in)    :: par
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: nangl, nangl_p
    integer                  :: i, j, angl_found, urey_found
    integer                  :: istart, iend
    integer                  :: i1, i2, i3, j1, j2, j3
    character(6)             :: ci1, ci2, ci3
    integer                  :: num_virtual_angles
    logical                  :: done

    integer                  :: k
    integer                  :: num_virtual_angle_types
    integer, parameter       :: max_virtual_angle_types = 10
    character(6)             :: ignored_types(3,max_virtual_angle_types)


    if(enefunc%qmmm%do_qmmm) then
      nangl   = molecule%num_angles - enefunc%qmmm%qm_nangles
    else
      nangl   = molecule%num_angles
    end if
    nangl_p = par%num_angles

    done = .true.
    angl_found = 0
    urey_found = 0
    num_virtual_angles = 0
    num_virtual_angle_types = 0

    do i = 1, nangl

      i1 = molecule%angl_list(1,i)
      i2 = molecule%angl_list(2,i)
      i3 = molecule%angl_list(3,i)

      ci1 = molecule%atom_cls_name(i1)
      ci2 = molecule%atom_cls_name(i2)
      ci3 = molecule%atom_cls_name(i3)

      do j = 1, nangl_p
        if ((ci1 == par%angl_atom_cls(1,j) .and. &
             ci2 == par%angl_atom_cls(2,j) .and. &
             ci3 == par%angl_atom_cls(3,j)) .or. &
            (ci1 == par%angl_atom_cls(3,j) .and. &
             ci2 == par%angl_atom_cls(2,j) .and. &
             ci3 == par%angl_atom_cls(1,j))) then
          angl_found = angl_found + 1

          if (abs(par%angl_ub_force_const(j)) > EPS) then
            urey_found = urey_found + 1
          end if

          exit

        end if
      end do

      if (j == nangl_p + 1) then
        done = .false.
        do j = max(1,i-2), min(nangl,i+2)
          j1 = molecule%angl_list(1,j)
          j2 = molecule%angl_list(2,j)
          j3 = molecule%angl_list(3,j)
          if ((i1 == j1 .and. i2 == j3 .and. i3 == j2) .or. &
              (i1 == j2 .and. i2 == j1 .and. i3 == j3) .or. &
              (i1 == j2 .and. i2 == j3 .and. i3 == j1) .or. &
              (i1 == j3 .and. i2 == j1 .and. i3 == j2)) then
            ! this must be water! ... I hope.
            done = .true.
            num_virtual_angles = num_virtual_angles + 1

            ci1 = molecule%atom_cls_name(i1)
            ci2 = molecule%atom_cls_name(i2)
            ci3 = molecule%atom_cls_name(i3)

            do k = 1, num_virtual_angle_types
              if ((ci1 == ignored_types(1,k) .and. &
                   ci2 == ignored_types(2,k) .and. &
                   ci3 == ignored_types(3,k)) .or. &
                  (ci1 == ignored_types(3,k) .and. &
                   ci2 == ignored_types(2,k) .and. &
                   ci3 == ignored_types(1,k))) then
                exit
              end if
            end do
            if (k == num_virtual_angle_types + 1) then
              num_virtual_angle_types = num_virtual_angle_types + 1
              if (num_virtual_angle_types >= max_virtual_angle_types) then
                call error_msg('Setup_Enefunc_Angl_Urey_Par> Too many water like parameter types.')
              end if
              ignored_types(1,num_virtual_angle_types) = ci1
              ignored_types(2,num_virtual_angle_types) = ci2
              ignored_types(3,num_virtual_angle_types) = ci3
            end if
            exit
          end if
        end do

        if (.not. done) &
          write(MsgOut,*) 'Setup_Enefunc_Angl_Urey_Par> not found ANGL: [', &
               ignored_types(1,num_virtual_angle_types), ']-[', &
               ignored_types(2,num_virtual_angle_types), ']-[', &
               ignored_types(3,num_virtual_angle_types), &
               '] in parameter file. (ERROR)'
      end if

    end do

    if (angl_found /= nangl - num_virtual_angles) &
      call error_msg('Setup_Enefunc_Angl_Urey_Par> Some Angle Paremeters are missing.')

    do i = 1, num_virtual_angle_types
      if (main_rank) then
        if (i == 1) write(MsgOut,*)
        write(MsgOut,'("Setup_Enefunc_Angl_Urey_Par> [",a6,"]-[",a6,"]-[",a6, &
                     & "] type angles ignored.")') &
            ignored_types(1,i), ignored_types(2,i), ignored_types(3,i)
      end if
    end do
    if (num_virtual_angle_types > 0 .and. main_rank) then
      write(MsgOut,*)
    end if

    enefunc%num_angles = angl_found
    enefunc%num_ureys  = urey_found

    ! allocate enefunc
    !
    call alloc_enefunc(enefunc, EneFuncAngl, enefunc%num_angles)
    call alloc_enefunc(enefunc, EneFuncUrey, enefunc%num_ureys)

    angl_found = 0
    urey_found = 0
    do i = 1, nangl

      ci1 = molecule%atom_cls_name(molecule%angl_list(1,i))
      ci2 = molecule%atom_cls_name(molecule%angl_list(2,i))
      ci3 = molecule%atom_cls_name(molecule%angl_list(3,i))
      !enefunc%angl_list(1:3,i) = molecule%angl_list(1:3,i)

      do j = 1, nangl_p
        if ((ci1 == par%angl_atom_cls(1,j) .and. &
             ci2 == par%angl_atom_cls(2,j) .and. &
             ci3 == par%angl_atom_cls(3,j)) .or. &
            (ci1 == par%angl_atom_cls(3,j) .and. &
             ci2 == par%angl_atom_cls(2,j) .and. &
             ci3 == par%angl_atom_cls(1,j))) then
          angl_found = angl_found + 1
          enefunc%angl_list(1:3,angl_found)    = molecule%angl_list(1:3,i)
          enefunc%angl_force_const(angl_found) = par%angl_force_const(j)
          enefunc%angl_theta_min(angl_found)   = par%angl_theta_min(j) * RAD

          if (abs(par%angl_ub_force_const(j)) > EPS) then
            urey_found = urey_found + 1
            enefunc%urey_list(1,urey_found) = molecule%angl_list(1,i)
            enefunc%urey_list(2,urey_found) = molecule%angl_list(3,i)
            enefunc%urey_force_const(urey_found) = par%angl_ub_force_const(j)
            enefunc%urey_rmin(urey_found)        = par%angl_ub_rmin(j)
          end if

          exit

        end if
      end do

    end do


    ! define angle interaction list for each processor
    !
    call get_loop_index(enefunc%num_angles, istart, iend)
    enefunc%istart_angle = istart
    enefunc%iend_angle   = iend

    call get_loop_index(enefunc%num_ureys, istart, iend)
    enefunc%istart_urey = istart
    enefunc%iend_urey   = iend

    return

  end subroutine setup_enefunc_angl_urey_par

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_dihe_par
  !> @brief        define DIHEDRAL term in potential energy function
  !! @authors      YS, JJ, TM
  !! @param[in]    par      : CHARMM PAR information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_dihe_par(par, molecule, enefunc)

    ! formal arguments
    type(s_par),             intent(in)    :: par
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: ndihe, ndihe_p
    integer                  :: i, j, found, nw_found, alloc_stat
    integer                  :: istart, iend
    character(6)             :: ci1, ci2, ci3, ci4

    logical,     allocatable :: no_wild(:)


    if(enefunc%qmmm%do_qmmm) then
      ndihe   = molecule%num_dihedrals - enefunc%qmmm%qm_ndihedrals
    else
      ndihe   = molecule%num_dihedrals
    end if
    ndihe_p = par%num_dihedrals

    ! check usage of wild card
    !
    allocate(no_wild(ndihe_p), stat=alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    do i = 1, ndihe_p

      if ((par%dihe_atom_cls(1,i) /= WildCard) .and. &
          (par%dihe_atom_cls(4,i) /= WildCard)) then
        ! A-B-C-D type
        no_wild(i) = .true.
      else
        ! X-B-C-D type
        no_wild(i) = .false.
      end if

    end do

    ! find number of interactions
    !
    found = 0
    do i = 1, ndihe

      ci1 = molecule%atom_cls_name(molecule%dihe_list(1,i))
      ci2 = molecule%atom_cls_name(molecule%dihe_list(2,i))
      ci3 = molecule%atom_cls_name(molecule%dihe_list(3,i))
      ci4 = molecule%atom_cls_name(molecule%dihe_list(4,i))

      nw_found = 0
      do j = 1, ndihe_p
        if (no_wild(j)) then
          if ((( ci1 == par%dihe_atom_cls(1, j)) .and. &
               ( ci2 == par%dihe_atom_cls(2, j)) .and. &
               ( ci3 == par%dihe_atom_cls(3, j)) .and. &
               ( ci4 == par%dihe_atom_cls(4, j))) .or. &
              (( ci1 == par%dihe_atom_cls(4, j)) .and. &
               ( ci2 == par%dihe_atom_cls(3, j)) .and. &
               ( ci3 == par%dihe_atom_cls(2, j)) .and. &
               ( ci4 == par%dihe_atom_cls(1, j)))) then
            nw_found = nw_found + 1
            found = found + 1
          end if
        end if
      end do

      if (nw_found == 0) then
        do j = 1, ndihe_p
          if (.not. no_wild(j)) then
            if ((( ci2 == par%dihe_atom_cls(2, j)) .and. &
                 ( ci3 == par%dihe_atom_cls(3, j))) .or. &
                (( ci2 == par%dihe_atom_cls(3, j)) .and. &
                 ( ci3 == par%dihe_atom_cls(2, j)))) then
              found = found + 1
            end if
          end if
        end do
     end if

    end do

    enefunc%num_dihedrals = found
    call alloc_enefunc(enefunc, EneFuncDihe, found)

    ! setup parameters
    !
    found = 0

    do i = 1, ndihe

      ci1 = molecule%atom_cls_name(molecule%dihe_list(1,i))
      ci2 = molecule%atom_cls_name(molecule%dihe_list(2,i))
      ci3 = molecule%atom_cls_name(molecule%dihe_list(3,i))
      ci4 = molecule%atom_cls_name(molecule%dihe_list(4,i))

      nw_found = 0
      do j = 1, ndihe_p
        if (no_wild(j)) then
          if ((( ci1 == par%dihe_atom_cls(1, j)) .and. &
               ( ci2 == par%dihe_atom_cls(2, j)) .and. &
               ( ci3 == par%dihe_atom_cls(3, j)) .and. &
               ( ci4 == par%dihe_atom_cls(4, j))) .or. &
              (( ci1 == par%dihe_atom_cls(4, j)) .and. &
               ( ci2 == par%dihe_atom_cls(3, j)) .and. &
               ( ci3 == par%dihe_atom_cls(2, j)) .and. &
               ( ci4 == par%dihe_atom_cls(1, j)))) then
            nw_found = nw_found + 1
            found = found + 1
            enefunc%dihe_list(1:4,found)    = molecule%dihe_list(1:4,i)
            enefunc%dihe_force_const(found) = par%dihe_force_const(j)
            enefunc%dihe_periodicity(found) = par%dihe_periodicity(j)
            enefunc%dihe_phase(found)       = par%dihe_phase(j) * RAD
          end if
        end if
      end do

      if (nw_found == 0) then
        do j = 1, ndihe_p
          if (.not. no_wild(j)) then
            if ((( ci2 == par%dihe_atom_cls(2, j)) .and. &
                 ( ci3 == par%dihe_atom_cls(3, j))) .or. &
                (( ci2 == par%dihe_atom_cls(3, j)) .and. &
                 ( ci3 == par%dihe_atom_cls(2, j)))) then
              found = found + 1
              enefunc%dihe_list(1:4,found)    = molecule%dihe_list(1:4,i)
              enefunc%dihe_force_const(found) = par%dihe_force_const(j)
              enefunc%dihe_periodicity(found) = par%dihe_periodicity(j)
              enefunc%dihe_phase(found)       = par%dihe_phase(j) * RAD
            end if
          end if
        end do
      end if

    end do

    enefunc%num_dihedrals = found

    if (found < ndihe) &
       call error_msg('Setup_Enefunc_Dihe_Par> Some Dihedral '//            &
                      'Paremeters are missing.')

    ! define dihedral angle interaction list for each processor
    !
    call get_loop_index(enefunc%num_dihedrals, istart, iend)
    enefunc%istart_dihedral = istart
    enefunc%iend_dihedral   = iend

    return

  end subroutine setup_enefunc_dihe_par

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_impr_par
  !> @brief        define IMPROPER term in potential energy function
  !! @authors      YS
  !! @param[in]    par      : CHARMM PAR information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_impr_par(par, molecule, enefunc)

    ! formal variables
    type(s_par),             intent(in)    :: par
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: nimpr, nimpr_p
    integer                  :: i, j, found, alloc_stat
    integer                  :: istart, iend
    character(6)             :: ci1, ci2, ci3, ci4

    integer,     allocatable :: wc_type(:)
    logical,     allocatable :: used(:)


    if(enefunc%qmmm%do_qmmm) then
      nimpr   = molecule%num_impropers - enefunc%qmmm%qm_nimpropers
    else
      nimpr   = molecule%num_impropers
    end if
    nimpr_p = par%num_impropers

    ! check usage of wild card
    !
    allocate(wc_type(nimpr_p), used(nimpr), stat=alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    do i = 1, nimpr_p

      if ((par%impr_atom_cls(1,i) /= WildCard) .and. &
          (par%impr_atom_cls(2,i) /= WildCard) .and. &
          (par%impr_atom_cls(3,i) /= WildCard) .and. &
          (par%impr_atom_cls(4,i) /= WildCard)) then
        ! A-B-C-D type
        wc_type(i) = 0

      else if ((par%impr_atom_cls(1,i) == WildCard) .and. &
               (par%impr_atom_cls(2,i) /= WildCard) .and. &
               (par%impr_atom_cls(3,i) /= WildCard) .and. &
               (par%impr_atom_cls(4,i) /= WildCard)) then
        ! X-B-C-D type
        wc_type(i) = 1

      else if ((par%impr_atom_cls(1,i) == WildCard) .and. &
               (par%impr_atom_cls(2,i) == WildCard) .and. &
               (par%impr_atom_cls(3,i) /= WildCard) .and. &
               (par%impr_atom_cls(4,i) /= WildCard)) then
        ! X-X-C-D type
        wc_type(i) = 2

      else if ((par%impr_atom_cls(1,i) /= WildCard) .and. &
               (par%impr_atom_cls(2,i) == WildCard) .and. &
               (par%impr_atom_cls(3,i) == WildCard) .and. &
               (par%impr_atom_cls(4,i) /= WildCard)) then
        ! A-X-X-D type
        wc_type(i) = 3

      else
        call error_msg('Setup_Enefunc_Impr_Par> Undefined Wild Card')

      end if

    end do

    ! find number of interactions
    !
    found = 0
    do i = 1, nimpr

      ci1 = molecule%atom_cls_name(molecule%impr_list(1,i))
      ci2 = molecule%atom_cls_name(molecule%impr_list(2,i))
      ci3 = molecule%atom_cls_name(molecule%impr_list(3,i))
      ci4 = molecule%atom_cls_name(molecule%impr_list(4,i))

      do j = 1, nimpr_p

        ! A-B-C-D type
        if (wc_type(j) == 0) then
          if (((ci1 == par%impr_atom_cls(1, j)) .and. &
               (ci2 == par%impr_atom_cls(2, j)) .and. &
               (ci3 == par%impr_atom_cls(3, j)) .and. &
               (ci4 == par%impr_atom_cls(4, j))) .or. & 
              ((ci1 == par%impr_atom_cls(4, j)) .and. &
               (ci2 == par%impr_atom_cls(3, j)) .and. &
               (ci3 == par%impr_atom_cls(2, j)) .and. &
               (ci4 == par%impr_atom_cls(1, j)))) then
            found = found + 1
            exit
          end if

        ! X-B-C-D type
        else if (wc_type(j) == 1) then
          if (((ci2 == par%impr_atom_cls(2, j)) .and. &
               (ci3 == par%impr_atom_cls(3, j)) .and. &
               (ci4 == par%impr_atom_cls(4, j))) .or. &
              ((ci2 == par%impr_atom_cls(4, j)) .and. &
               (ci3 == par%impr_atom_cls(3, j)) .and. &
               (ci4 == par%impr_atom_cls(2, j)))) then
            found = found + 1
            exit
          end if

        ! X-X-C-D type
        else if (wc_type(j) == 2) then
          if (((ci3 == par%impr_atom_cls(3, j)) .and. &
               (ci4 == par%impr_atom_cls(4, j))) .or. &
              ((ci3 == par%impr_atom_cls(4, j)) .and. &
               (ci4 == par%impr_atom_cls(3, j)))) then
            found = found + 1
            exit
          end if

        ! A-X-X-D type
        else if (wc_type(j) == 3) then
          if (((ci1 == par%impr_atom_cls(1, j)) .and. &
               (ci4 == par%impr_atom_cls(4, j))) .or. &
              ((ci1 == par%impr_atom_cls(4, j)) .and. &
               (ci4 == par%impr_atom_cls(1, j)))) then
            found = found + 1
            exit
          end if

        else
          write(MsgOut,*) i, ci1, ci2, ci3, ci4
          call error_msg('Setup_Enefunc_Impr_Par> Unknown IMPR type.')

        end if

      end do

    end do
  
    enefunc%num_impropers = found
    call alloc_enefunc(enefunc, EneFuncImpr, found)

    ! setup parameters
    !
    found = 0

    do i = 1, nimpr

      ci1 = molecule%atom_cls_name(molecule%impr_list(1,i))
      ci2 = molecule%atom_cls_name(molecule%impr_list(2,i))
      ci3 = molecule%atom_cls_name(molecule%impr_list(3,i))
      ci4 = molecule%atom_cls_name(molecule%impr_list(4,i))

      used(i) = .false.

      ! A-B-C-D type
      do j = 1, nimpr_p
        if (wc_type(j) == 0) then
          if (((ci1 == par%impr_atom_cls(1, j)) .and. &
               (ci2 == par%impr_atom_cls(2, j)) .and. &
               (ci3 == par%impr_atom_cls(3, j)) .and. &
               (ci4 == par%impr_atom_cls(4, j))) .or. &
              ((ci1 == par%impr_atom_cls(4, j)) .and. &
               (ci2 == par%impr_atom_cls(3, j)) .and. &
               (ci3 == par%impr_atom_cls(2, j)) .and. &
               (ci4 == par%impr_atom_cls(1, j)))) then
            found = found + 1
            enefunc%impr_list(1:4,found)    = molecule%impr_list(1:4,i)
            enefunc%impr_force_const(found) = par%impr_force_const(j)
            enefunc%impr_phase(found)       = par%impr_phase(j) * RAD
            used(i) = .true.
            exit

          end if
        end if
      end do

      ! X-B-C-D type
      if (.not.used(i)) then
        do j = 1, nimpr_p
          if (wc_type(j) == 1) then
            if (((ci2 == par%impr_atom_cls(2, j)) .and. &
                 (ci3 == par%impr_atom_cls(3, j)) .and. &
                 (ci4 == par%impr_atom_cls(4, j))) .or. &
                ((ci2 == par%impr_atom_cls(4, j)) .and. &
                 (ci3 == par%impr_atom_cls(3, j)) .and. &
                 (ci4 == par%impr_atom_cls(2, j)))) then
              found = found + 1
              enefunc%impr_list(1:4,found)    = molecule%impr_list(1:4,i)
              enefunc%impr_force_const(found) = par%impr_force_const(j)
              enefunc%impr_phase(found)       = par%impr_phase(j) * RAD
              used(i) = .true.
              exit

            end if
          end if
        end do
      end if

      ! X-X-C-D type
      if (.not.used(i)) then
        do j = 1, nimpr_p
          if (wc_type(j) == 2) then
            if (((ci3 == par%impr_atom_cls(3, j)) .and. &
                 (ci4 == par%impr_atom_cls(4, j))) .or. &
                ((ci3 == par%impr_atom_cls(4, j)) .and. &
                 (ci4 == par%impr_atom_cls(3, j)))) then
              found = found + 1
              enefunc%impr_list(1:4,found)    = molecule%impr_list(1:4,i)
              enefunc%impr_force_const(found) = par%impr_force_const(j)
              enefunc%impr_phase(found)       = par%impr_phase(j) * RAD
              used(i) = .true.
              exit

            end if
          end if
        end do
      end if

      ! A-X-X-D type
      if (.not.used(i)) then
        do j = 1, nimpr_p
          if (wc_type(j) == 3) then
            if (((ci1 == par%impr_atom_cls(1, j)) .and. &
                 (ci4 == par%impr_atom_cls(4, j))) .or. &
                ((ci1 == par%impr_atom_cls(4, j)) .and. &
                 (ci4 == par%impr_atom_cls(1, j)))) then
              found = found + 1
              enefunc%impr_list(1:4,found)    = molecule%impr_list(1:4,i)
              enefunc%impr_force_const(found) = par%impr_force_const(j)
              enefunc%impr_phase(found)       = par%impr_phase(j) * RAD
              used(i) = .true.
              exit

            end if
          end if
        end do
      end if

      ! Error
      if (.not.used(i)) then
        write(MsgOut,*) i, ci1, ci2, ci3, ci4
        call error_msg('Setup_Enefunc_Impr_Par> Unknown IMPR type.')
      end if

    end do

    enefunc%num_impropers = found

    if (found < nimpr) &
       call error_msg( &
            'Setup_Enefunc_Impr_Par> Some Improper Paremeters are missing.')

    ! define improper dihedral angle interaction list for each processor
    !
    call get_loop_index(enefunc%num_impropers, istart, iend)
    enefunc%istart_improper = istart
    enefunc%iend_improper   = iend

    deallocate(wc_type,used)

    return

  end subroutine setup_enefunc_impr_par

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cmap_par
  !> @brief        define cmap term in potential energy function
  !! @authors      TY
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[in]    par      : CHARMM PAR information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !!
  !! @note         In str "par", following variables have been defined:
  !!   cmap_atom_cls(8,imap) (char) : 8 atom classes (4 for phi and 4 for psi)
  !!   cmap_resolution(imap) (int ) : = 24 (for all imap) = 360degree/15degree
  !!   cmap_data(i,j,imap)   (real) : Ecmap for grid points. i and j are the
  !!                                  1st (psi?) and the 2nd (phi?) grid IDs.
  !!                                  1<=i,j<=24.
  !!   Where imap is ID of cmap type (1 <= imap <= 6).
  !!
  !! @note         TY determined to use natural spline (periodic = .false.)
  !!               because it is similar to the CHARMM way.
  !!               However, I notice that force will be more accurately 
  !!               continuous at phi (or psi) = +-180 degree if periodic spline
  !!               was used.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cmap_par(ene_info, par, molecule, enefunc)

    ! formal variables
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_par),             intent(in)    :: par
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


    ! If 'periodic' is .true.,  
    !   then cubic spline with periodic (in dihedral-angle space) boundary 
    !   will be applied.
    ! If 'periodic' is .false., 
    !   then natural cubic spline without periodic boudnary 
    !   will be applied to expanded cross-term maps. 
    !   This is similar with CHARMM's source code.  
    !
    periodic = ene_info%cmap_pspline
    if(enefunc%qmmm%do_qmmm) then
      ncmap    = molecule%num_cmaps - enefunc%qmmm%qm_ncmaps
    else
      ncmap    = molecule%num_cmaps
    end if
    ncmap_p  = par%num_cmaps


    !  begin
    ! 
    ngrid0 = 0
    do i = 1, ncmap_p
      ngrid0 = max(ngrid0, par%cmap_resolution(i))
    end do

    call alloc_enefunc(enefunc, EneFuncCmap, ncmap, ngrid0)
    call alloc_enefunc(enefunc, EneFuncCmapType, ncmap_p, ngrid0)


    ! allocate local arrays
    !
    alloc_stat = 0
    allocate(c_ij(4,4,ngrid0,ngrid0), &
             stat = alloc_stat )
    if (alloc_stat /= 0) &
      call error_msg_alloc


    ! find Ecmap tables that will be used
    !
    enefunc%cmap_type(1:ncmap) = -1

    found = 0
    do i = 1, ncmap

      ! ci* will be atom-type strings
      !
      ci1 = molecule%atom_cls_name(molecule%cmap_list(1,i))
      ci2 = molecule%atom_cls_name(molecule%cmap_list(2,i))
      ci3 = molecule%atom_cls_name(molecule%cmap_list(3,i))
      ci4 = molecule%atom_cls_name(molecule%cmap_list(4,i))
      ci5 = molecule%atom_cls_name(molecule%cmap_list(5,i))
      ci6 = molecule%atom_cls_name(molecule%cmap_list(6,i))
      ci7 = molecule%atom_cls_name(molecule%cmap_list(7,i))
      ci8 = molecule%atom_cls_name(molecule%cmap_list(8,i))
      enefunc%cmap_list(1:8,i) = molecule%cmap_list(1:8,i)

      ! assign cmap type ID to each (psi,phi) pair
      !
      do j = 1, ncmap_p
        if (ci1 == par%cmap_atom_cls(1, j) .and. &
            ci2 == par%cmap_atom_cls(2, j) .and. &
            ci3 == par%cmap_atom_cls(3, j) .and. &
            ci4 == par%cmap_atom_cls(4, j) .and. &
            ci5 == par%cmap_atom_cls(5, j) .and. &
            ci6 == par%cmap_atom_cls(6, j) .and. &
            ci7 == par%cmap_atom_cls(7, j) .and. &
            ci8 == par%cmap_atom_cls(8, j)) then

          enefunc%cmap_type(i) = j
          found = found + 1
          exit

        end if
      end do

      ! if not found, print detail about the error.
      !
      if (enefunc%cmap_type(i) <= 0) then
        write(MsgOut,*) 'Setup_Enefunc_Cmap_Par> not found CMAP:'
        write(MsgOut,*) '[',ci1,']-[',ci2,']-[',ci3,']-[',ci4,']'
        write(MsgOut,*) '-[',ci5,']-[',ci6,']-[',ci7,']-[',ci8,']'
        write(MsgOut,*)  'in parameter file. (ERROR)'
      end if

    end do

    ! stop if parameter is missing
    !
    if (found < ncmap) &
      call error_msg('Setup_Enefunc_Cmap_Par> '//                       &
                     'Some Cmap Parameters are missing.')

    !  derive cmap coefficients by bicubic interpolation
    !    
    enefunc%num_cmaps = found

    do ityp = 1, ncmap_p
      enefunc%cmap_resolution(ityp) = par%cmap_resolution(ityp)
    end do

    do ityp = 1, ncmap_p
      if (periodic) then
        call derive_cmap_coefficients_p (ityp, par, c_ij)
      else
        call derive_cmap_coefficients_np(ityp, par, c_ij)
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
    if (dealloc_stat /= 0) &
      call error_msg_dealloc

    ! define improper dihedral angle interaction list for each processor
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

  end subroutine setup_enefunc_cmap_par

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_nonb_par
  !> @brief        define NON-BOND term in potential energy function
  !! @authors      YS, TI, NT
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[in]    par      : CHARMM PAR information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_nonb_par(ene_info, par, molecule, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_par),             intent(in)    :: par
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: eps14, rmin14, eps, rmin
    integer                  :: i, j, k
    integer                  :: nonb_p, nbfx_p
    integer                  :: istart, iend
    character(6)             :: ci1, ci2


    nonb_p = par%num_atom_cls
    enefunc%num_atom_cls = nonb_p

    ELECOEF = ELECOEF_CHARMM

    ! set lennard-jones parameters
    !
    call alloc_enefunc(enefunc, EneFuncNbon, nonb_p)

    do i = 1, nonb_p
      do j = 1, nonb_p

        ! combination rule
        !
        eps14  = sqrt(par%nonb_eps_14(i) * par%nonb_eps_14(j))
        rmin14 = par%nonb_rmin_14(i) + par%nonb_rmin_14(j)
        eps    = sqrt(par%nonb_eps(i) * par%nonb_eps(j))
        rmin   = par%nonb_rmin(i) + par%nonb_rmin(j)

        ! set parameters
        !
        if (enefunc%forcefield == ForcefieldKBGO) then
          enefunc%nb14_lj12(i,j) = eps14 * (rmin14 ** 12)
          enefunc%nb14_lj10(i,j) = 0.0_wp
          enefunc%nb14_lj6 (i,j) = 0.0_wp
          enefunc%nonb_lj12(i,j) = eps * (rmin ** 12)
          enefunc%nonb_lj10(i,j) = 0.0_wp
          enefunc%nonb_lj6 (i,j) = 0.0_wp
        else
          enefunc%nb14_lj12(i,j) = eps14 * (rmin14 ** 12)
          enefunc%nb14_lj10(i,j) = 0.0_wp
          enefunc%nb14_lj6 (i,j) = 2.0_wp * eps14 * (rmin14 ** 6) 
          enefunc%nonb_lj12(i,j) = eps * (rmin ** 12)
          enefunc%nonb_lj10(i,j) = 0.0_wp
          enefunc%nonb_lj6 (i,j) = 2.0_wp * eps * (rmin ** 6) 
        end if

      end do
    end do
    enefunc%nb14_qq_scale_c19 = 1.0_wp
    if (enefunc%forcefield == ForcefieldCHARMM19) then
      enefunc%nb14_qq_scale_c19 = 0.4_wp
    end if

    if (enefunc%forcefield == ForcefieldKBGO) then
      call alloc_enefunc(enefunc, EneFuncNonbGO, nonb_p)
      do i = 1, nonb_p
        enefunc%nonb_eps(i) = par%nonb_eps(i)
        enefunc%nonb_rmin(i) = par%nonb_rmin(i)
      end do

    end if


    ! overwrite lennard-jones parameters by NBFIX parameters
    !
    nbfx_p = par%num_nbfix
    if (enefunc%forcefield == ForcefieldKBGO) then
      enefunc%num_contacts = nbfx_p
      call alloc_enefunc(enefunc, EneFuncCntc, nbfx_p)
    end if

    do k = 1, nbfx_p
      
      ci1 = par%nbfi_atom_cls(1,k)
      ci2 = par%nbfi_atom_cls(2,k)

      do i = 1, nonb_p
        do j = 1, nonb_p

          if (ci1 .eq. par%nonb_atom_cls(i) .and. &
              ci2 .eq. par%nonb_atom_cls(j) .or.  &
              ci2 .eq. par%nonb_atom_cls(i) .and. &
              ci1 .eq. par%nonb_atom_cls(j)) then

            ! combination rule
            !
            eps14  = abs(par%nbfi_eps_14 (k)) !TODO CHECK
            rmin14 = par%nbfi_rmin_14(k)
            eps    = abs(par%nbfi_eps    (k))
            rmin   = par%nbfi_rmin   (k)

            ! set parameters
            !
            if (enefunc%forcefield == ForcefieldKBGO) then
              enefunc%contact_lj12(k) = 13.0_wp * eps * (rmin ** 12)
              enefunc%contact_lj10(k) = 18.0_wp * eps * (rmin ** 10)
              enefunc%contact_lj6 (k) =  4.0_wp * eps * (rmin ** 6)
              enefunc%contact_list(1,k) = i
              enefunc%contact_list(2,k) = j
            else
              enefunc%nb14_lj12(i,j) = eps14 * (rmin14 ** 12)
              enefunc%nb14_lj10(i,j) = 0.0_wp
              enefunc%nb14_lj6 (i,j) = 2.0_wp * eps14 * (rmin14 ** 6)
              enefunc%nonb_lj12(i,j) = eps * (rmin ** 12)
              enefunc%nonb_lj10(i,j) = 0.0_wp
              enefunc%nonb_lj6 (i,j) = 2.0_wp * eps * (rmin ** 6)
            end if

          end if

        end do
      end do

    end do
    if (enefunc%forcefield == ForcefieldKBGO) then
       call get_loop_index(enefunc%num_contacts, istart, iend)
       enefunc%istart_contact = istart
       enefunc%iend_contact   = iend
    end if


    ! treatment for 1-2, 1-3, 1-4 interactions
    !
    if (.not.ene_info%table) then
      if (enefunc%forcefield .eq. ForcefieldCHARMM19) then
        call count_nonb_excl_par_c19(molecule, enefunc)
      else
        call count_nonb_excl_par(molecule, enefunc)
      end if
    end if

    if (ene_info%num_basins > 1) then
        call error_msg( &
          'Setup_Enefunc_Nonb_Par> multibasin is not allowed in CHARMM-format')
    end if

    return

  end subroutine setup_enefunc_nonb_par

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    count_nonb_excl_par
  !> @brief        exclude 1-2, 1-3 interactions
  !! @authors      YS, TM
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine count_nonb_excl_par(molecule, enefunc)

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

    allocate(nb14_calc_list(MAX_NB14,molecule%num_atoms),stat=alloc_stat)
    if (alloc_stat /=0) &
      call error_msg_alloc

    enefunc%num_nonb_excl(1:molecule%num_atoms) = 0
    enefunc%num_nb14_calc(1:molecule%num_atoms) = 0

    num_excl_total = 0

    ! exclude native contacts
    !
    if (enefunc%forcefield == ForcefieldKBGO) then
      do k = 1, enefunc%num_contacts
        k1 = enefunc%contact_list(1,k)
        k2 = enefunc%contact_list(2,k)
   
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
    end if

    ! exclude 1-2 interaction
    !

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
             stat = alloc_stat)
    if (alloc_stat /=0) &
      call error_msg_alloc

    do i = 1, molecule%num_atoms
      enefunc%nb14_calc_list(1:enefunc%num_nb14_calc(i),i) = &
             nb14_calc_list(1:enefunc%num_nb14_calc(i),i)
    end do

    deallocate(nb14_calc_list, stat = dealloc_stat)

    return

  end subroutine count_nonb_excl_par

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    count_nonb_excl_par_c19
  !> @brief        exclude 1-2, 1-3 interactions
  !! @authors      YS, TM, CK
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine count_nonb_excl_par_c19(molecule, enefunc)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, k, k1, k2, k3, k4
    integer                  :: ik1, ik4
    integer                  :: num_excl, num_excl_total
    integer                  :: num_nb14, num_nb14_total
    integer                  :: max_nb14_num
    integer                  :: alloc_stat, dealloc_stat
    logical                  :: duplicate, improp

    integer,     allocatable :: nonb_excl_list(:,:)
    integer,     allocatable :: nb14_calc_list(:,:)
    integer,     allocatable :: temporary_list(:,:)


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

    allocate(nb14_calc_list(MAX_NB14,molecule%num_atoms),stat=alloc_stat)
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
    call create_nbond_list(molecule%bond_list, 3, temporary_list)

    do k = 1, size(temporary_list(1,:))

      k1 = temporary_list(1,k)
      k3 = temporary_list(3,k)

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
    deallocate(temporary_list)

    ! count 1-4 interaction
    !
    call create_nbond_list(molecule%bond_list, 4, temporary_list)
    num_nb14_total = 0
    do k = 1, size(temporary_list(1,:))

      k1 = temporary_list(1,k)
      k4 = temporary_list(4,k)
      improp = .false.

      do i = 1, enefunc%num_impropers
        ik1 = enefunc%impr_list(1,i)
        ik4 = enefunc%impr_list(4,i)
        if (min(k1,k4) == min(ik1,ik4) .and. &
            max(k1,k4) == max(ik1,ik4)) then
          improp = .true.
          exit
        end if
      end do

      if (k1 < k4) then

        num_nb14 = enefunc%num_nb14_calc(k1)
        num_excl = enefunc%num_nonb_excl(k1)
        duplicate = .false.

        do i = 1, enefunc%num_nonb_excl(k1)
          if (k4 == nonb_excl_list(i,k1)) duplicate = .true.
        end do

        do i = 1, num_nb14
          if (k4 == nb14_calc_list(i,k1)) duplicate = .true.
        end do

        if (.not. duplicate) then
          if (.not. improp) then
            num_nb14 = num_nb14 + 1
            enefunc%num_nb14_calc(k1) = num_nb14
            nb14_calc_list(num_nb14,k1) = k4
            num_nb14_total = num_nb14_total + 1
          else
            num_excl = num_excl + 1
            enefunc%num_nonb_excl(k1) = num_excl
            nonb_excl_list(num_excl,k1) = k4
            num_excl_total = num_excl_total + 1
          end if
        end if

      else

        num_nb14 = enefunc%num_nb14_calc(k4)
        num_excl = enefunc%num_nonb_excl(k4)
        duplicate = .false.

        do i = 1, enefunc%num_nonb_excl(k4)
          if (k1 == nonb_excl_list(i,k4)) duplicate = .true.
        end do

        do i = 1, num_nb14
          if (k1 == nb14_calc_list(i,k4)) duplicate = .true.
        end do

        if (.not. duplicate) then
          if (.not. improp) then
            num_nb14 = num_nb14 + 1
            enefunc%num_nb14_calc(k4) = num_nb14
            nb14_calc_list(num_nb14,k4) = k1
            num_nb14_total = num_nb14_total + 1
          else
            num_excl = num_excl + 1
            enefunc%num_nonb_excl(k4) = num_excl
            nonb_excl_list(num_excl,k4) = k1
            num_excl_total = num_excl_total + 1
          end if
        end if
      end if

    end do
    deallocate(temporary_list)
    
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
             stat = alloc_stat)
    if (alloc_stat /=0) &
      call error_msg_alloc
    do i = 1, molecule%num_atoms
      enefunc%nb14_calc_list(1:enefunc%num_nb14_calc(i),i) = &
             nb14_calc_list(1:enefunc%num_nb14_calc(i),i)
    end do

    deallocate(nb14_calc_list, stat = dealloc_stat)

    return

  end subroutine count_nonb_excl_par_c19

end module at_enefunc_charmm_mod
