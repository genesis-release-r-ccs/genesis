!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_enefunc_go_mod
!> @brief   define potential energy functions
!! @authors Takaharu Mori (TM), Chigusa Kobayashi (CK), Yuji Sugita (YS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_enefunc_go_mod

  use at_enefunc_restraints_mod
  use at_energy_mod
  use at_restraints_str_mod
  use at_enefunc_str_mod
  use math_libs_mod
  use molecules_str_mod
  use fileio_gpr_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! constants
  integer, parameter :: MAXWRK   = 12
  integer, parameter :: MAX_EXCL = 16 ! = 4(1-2) + 4x3 (1-3)
  integer, parameter :: MAX_NB14 = 36 ! = 3x4x3 (1-4) 

  ! subroutines
  public  :: define_enefunc_go
  private :: setup_enefunc_bond_gpr
  private :: setup_enefunc_angl_gpr
  private :: setup_enefunc_dihe_gpr
  private :: setup_enefunc_impr_gpr
  private :: setup_enefunc_nonb_gpr
  private :: count_nonb_excl_gpr

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc_go
  !> @brief        a driver subroutine for defining potential energy
  !! @authors      TM, CK
  !! @param[in]    ene_info   : ENERGY section control parameters information
  !! @param[in]    gpr        : GO model parameter information
  !! @param[in]    molecule   : molecule information
  !! @param[in]    restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc_go(ene_info, gpr, molecule, restraints, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_gpr),             intent(in)    :: gpr
    type(s_molecule),        intent(in)    :: molecule
    type(s_restraints),      intent(in)    :: restraints
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: alloc_stat
    integer                  :: i, found, found2
    integer                  :: nwork


    ! bond
    !
    call setup_enefunc_bond_gpr(gpr, enefunc)

    ! angle
    !
    call setup_enefunc_angl_gpr(gpr, enefunc)

    ! dihedral
    !
    call setup_enefunc_dihe_gpr(gpr, enefunc)

    ! improper
    !
    call setup_enefunc_impr_gpr(gpr, enefunc)

    ! nonbonded (contact and non-contact)
    !   Note that contact list is included in nonb_excl_list!!
    !
    call setup_enefunc_nonb_gpr(ene_info, gpr, molecule, enefunc)

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
                enefunc%num_restraintfuncs, &
                enefunc%num_restraintgroups,&
                enefunc%num_morph_sc)

    allocate(enefunc%work(MAXWRK,nwork),stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    ! write summary of energy function
    !
    found  = 0
    found2 = 0

    do i = 1, molecule%num_atoms
      found  = found  + enefunc%num_nonb_excl(i)
      found2 = found2 + enefunc%num_nb14_calc(i)
    end do

    write(MsgOut,'(A)') 'Define_Enefunc_Go> Number of Interactions in Each Term'
    write(MsgOut,'(A20,I10,A20,I10)')                         &
         '  bond_ene        = ', enefunc%num_bonds,           &
         '  angle_ene       = ', enefunc%num_angles
    write(MsgOut,'(A20,I10,A20,I10)')                         &
         '  torsion_ene     = ', enefunc%num_dihedrals,       &
         '  improper_ene    = ', enefunc%num_impropers
    write(MsgOut,'(A20,I10,A20,I10)')                         &
         '  contact         = ', enefunc%num_contacts,        &
         '  nc exclusions   = ', found + found2
    write(MsgOut,'(A20,I10,A20,I10)')                         &
         '  restraint_group = ', enefunc%num_restraintgroups, &
         '  restraint_func  = ', enefunc%num_restraintfuncs
    
    write(MsgOut,'(A)') ' '

    return

  end subroutine define_enefunc_go

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_bond_gpr
  !> @brief        define BOND term in potential energy function
  !! @authors      TM
  !! @param[in]    gpr     : GO model parameter information
  !! @param[inout] enefunc : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond_gpr(gpr, enefunc)

    ! formal arguments
    type(s_gpr),             intent(in)    :: gpr
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, nbond
    integer                  :: istart, iend


    nbond = gpr%num_bonds
    call alloc_enefunc(enefunc, EneFuncBond, nbond)

    do i = 1, nbond
      enefunc%bond_list(1:2,i)    = gpr%bond_list(1:2,i)
      enefunc%bond_dist_min(i)    = gpr%bond_dist_min(i)
      enefunc%bond_force_const(i) = gpr%bond_force_const(i) 
    end do
    enefunc%num_bonds = nbond

    ! define bond interaction list for each processor
    !
    call get_loop_index(enefunc%num_bonds, istart, iend)
    enefunc%istart_bond = istart
    enefunc%iend_bond   = iend

    return

  end subroutine setup_enefunc_bond_gpr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_angl_gpr
  !> @brief        define ANGLE term in potential energy function
  !! @authors      TM
  !! @param[in]    gpr     : GO model parameter information
  !! @param[inout] enefunc : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_angl_gpr(gpr, enefunc)

    ! formal arguments
    type(s_gpr),             intent(in)    :: gpr
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, nangl
    integer                  :: istart, iend


    nangl = gpr%num_angles
    call alloc_enefunc(enefunc, EneFuncAngl, nangl)

    do i = 1, nangl
      enefunc%angl_list(1:3,i)       = gpr%angl_list(1:3,i)
      enefunc%angl_force_const(i)    = gpr%angl_force_const(i)
      enefunc%angl_theta_min(i)      = gpr%angl_theta_min(i) * RAD
    end do
    enefunc%num_angles = nangl

    ! define angle interaction list for each processor
    !
    call get_loop_index(enefunc%num_angles, istart, iend)
    enefunc%istart_angle = istart
    enefunc%iend_angle   = iend

    return

  end subroutine setup_enefunc_angl_gpr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_dihe_gpr
  !> @brief        define DIHEDRAL term in potential energy function
  !! @authors      TM
  !! @param[in]    gpr     : GO model parameter information
  !! @param[inout] enefunc : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_dihe_gpr(gpr, enefunc)

    ! formal arguments
    type(s_gpr),             intent(in)    :: gpr
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, ndihe
    integer                  :: istart, iend


    ndihe = gpr%num_dihedrals
    call alloc_enefunc(enefunc, EneFuncDihe, ndihe)

    do i = 1, ndihe
      enefunc%dihe_list(1:4,i)    = gpr%dihe_list(1:4,i)
      enefunc%dihe_force_const(i) = gpr%dihe_force_const(i)
      enefunc%dihe_periodicity(i) = gpr%dihe_periodicity(i)
      enefunc%dihe_phase(i)       = gpr%dihe_phase(i) * RAD
    end do
    enefunc%num_dihedrals = ndihe

    ! define dihedral angle interaction list for each processor
    !
    call get_loop_index(enefunc%num_dihedrals, istart, iend)
    enefunc%istart_dihedral = istart
    enefunc%iend_dihedral   = iend

    return

  end subroutine setup_enefunc_dihe_gpr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_impr_gpr
  !> @brief        define IMPROPER term in potential energy function
  !! @authors      TM
  !! @param[in]    gpr     : GO model parameter information
  !! @param[inout] enefunc : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_impr_gpr(gpr, enefunc)

    ! formal arguments
    type(s_gpr),             intent(in)    :: gpr
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, nimpr
    integer                  :: istart, iend


    nimpr = gpr%num_impropers
    call alloc_enefunc(enefunc, EneFuncImpr, nimpr)

    do i = 1, nimpr
      enefunc%impr_list(1:4,i)    = gpr%impr_list(1:4,i)
      enefunc%impr_force_const(i) = gpr%impr_force_const(i)
      enefunc%impr_phase(i)       = gpr%impr_phase(i) * RAD
    end do
    enefunc%num_impropers = nimpr

    ! define improper dihedral angle interaction list for each processor
    !
    call get_loop_index(enefunc%num_impropers, istart, iend)
    enefunc%istart_improper = istart
    enefunc%iend_improper   = iend

    return

  end subroutine setup_enefunc_impr_gpr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_nonb_gpr
  !> @brief        define NON-BOND term in potential energy function
  !! @authors      TM
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[in]    gpr      : GO model parameter information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_nonb_gpr(ene_info, gpr, molecule, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_gpr),             intent(in)    :: gpr
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: eps, rmin
    integer                  :: i, npair
    integer                  :: istart, iend


    ! set contact pairs
    !
    npair = gpr%num_pairs
    call alloc_enefunc(enefunc, EneFuncCntc, npair)

    if (enefunc%forcefield == ForcefieldAAGO) then
      do i = 1, npair
        eps  = gpr%pair_eps(i)
        rmin = gpr%pair_dist_min(i)
        enefunc%contact_list(1:2,i) = gpr%pair_list(1:2,i)
        enefunc%contact_lj12(i)     =         eps * (rmin ** 12)
        enefunc%contact_lj6(i)      = 2.0_wp * eps * (rmin ** 6)
      end do

    else if (enefunc%forcefield == ForcefieldCAGO) then
      do i = 1, npair
        eps  = gpr%pair_eps(i)
        rmin = gpr%pair_dist_min(i)
        enefunc%contact_list(1:2,i) = gpr%pair_list(1:2,i)
        enefunc%contact_lj12(i)     = 5.0_wp * eps * (rmin ** 12)
        enefunc%contact_lj10(i)     = 6.0_wp * eps * (rmin ** 10)
      end do

    end if
    enefunc%num_contacts = npair


    ! define contact interaction list for each processor
    !
    call get_loop_index(enefunc%num_contacts, istart, iend)
    enefunc%istart_contact = istart
    enefunc%iend_contact   = iend


    ! set non-contact pairs
    !
    eps  = gpr%nonpair_eps
    rmin = gpr%nonpair_dist_min
    enefunc%noncontact_lj12 = eps * (rmin ** 12)


    ! set exclusion list including 1-2, 1-3, and contact pairs
    !   We also make 1-4 interaction list including dihedral and improper lists,
    !   but this is just used for exclusion for non-contact interactions.
    !
    call count_nonb_excl_gpr(gpr, molecule, enefunc)

    return

  end subroutine setup_enefunc_nonb_gpr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    count_nonb_excl_gpr
  !> @brief        exclude 1-2, 1-3 interactions
  !! @authors      TM
  !! @param[in]    gpr      : GO model parameter information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !! 
  !! @note         do NOT use this subroutine except for "ALL-ATOM" go model
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine count_nonb_excl_gpr(gpr, molecule, enefunc)

    ! formal arguments
    type(s_gpr),             intent(in)    :: gpr
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, k, k1, k2, k3, k4
    integer                  :: num_excl, num_excl_total
    integer                  :: num_nb14, num_nb14_total
    integer                  :: max_nb14_num
    integer                  :: max_contact, max_excl_go
    integer                  :: alloc_stat, dealloc_stat
    logical                  :: duplicate

    integer,     allocatable :: nonb_excl_list(:,:)
    integer,     allocatable :: nb14_calc_list(:,:)
    integer,     allocatable :: num_contact_excl(:)


    ! redefine max_excl
    !
    allocate(num_contact_excl(molecule%num_atoms), stat=alloc_stat)
    if (alloc_stat /=0) &
      call error_msg_alloc

    num_contact_excl(1:molecule%num_atoms) = 0

    do k = 1, gpr%num_pairs

      k1 = gpr%pair_list(1,k)
      k2 = gpr%pair_list(2,k)

      if (k1 < k2) then
        num_excl = num_contact_excl(k1) + 1
        num_contact_excl(k1) = num_excl
      else
        num_excl = num_contact_excl(k2) + 1
        num_contact_excl(k2) = num_excl
      end if

    end do

    max_contact = maxval(num_contact_excl)
    max_excl_go = MAX_EXCL + max_contact

    deallocate(num_contact_excl, stat = dealloc_stat)


    ! allocate exclusion list and 1-4 interaction list
    !
    allocate(enefunc%num_nonb_excl(molecule%num_atoms), stat=alloc_stat)
    if (alloc_stat /=0) &
      call error_msg_alloc

    allocate(nonb_excl_list(max_excl_go,molecule%num_atoms),stat=alloc_stat)
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


    ! exclude contact pairs
    !
    num_excl_total = 0

    do k = 1, gpr%num_pairs

      k1 = gpr%pair_list(1,k)
      k2 = gpr%pair_list(2,k)

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


    ! exclude 1-2 interaction
    !
    do k = 1, molecule%num_bonds

      k1 = molecule%bond_list(1,k)
      k2 = molecule%bond_list(2,k)

      if (k1 < k2) then

        num_excl = enefunc%num_nonb_excl(k1)
        duplicate = .false.

        do i = 1, enefunc%num_nonb_excl(k1)
          if (k2 == nonb_excl_list(i,k1)) duplicate = .true.
        end do

        if (.not. duplicate) then
          num_excl = num_excl + 1
          enefunc%num_nonb_excl(k1) = num_excl
          nonb_excl_list(num_excl,k1) = k2
          num_excl_total = num_excl_total + 1
        end if

      else

        num_excl = enefunc%num_nonb_excl(k2)
        duplicate = .false.

        do i = 1, enefunc%num_nonb_excl(k2)
          if (k1 == nonb_excl_list(i,k2)) duplicate = .true.
        end do 

        if (.not. duplicate) then
          num_excl = num_excl + 1
          enefunc%num_nonb_excl(k2) = num_excl
          nonb_excl_list(num_excl,k2) = k1
          num_excl_total = num_excl_total + 1
        end if

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

    do k = 1, molecule%num_impropers

      k1 = molecule%impr_list(1,k)
      k4 = molecule%impr_list(4,k)

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

    enefunc%nb14_calc_list(1:max_nb14_num,1:molecule%num_atoms) = 0
    do i = 1, molecule%num_atoms
      enefunc%nb14_calc_list(1:enefunc%num_nb14_calc(i),i) = &
             nb14_calc_list(1:enefunc%num_nb14_calc(i),i)
    end do

    deallocate(nb14_calc_list, stat = dealloc_stat)

    return

  end subroutine count_nonb_excl_gpr

end module at_enefunc_go_mod
