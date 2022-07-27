!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_enefunc_gbsa_mod
!> @brief   define enefuncs for implicit solvent models
!! @authors Takaharu Mori (TM)
!
!  (c) Copyright 2016 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_enefunc_gbsa_mod

  use atom_libs_mod
  use at_energy_mod
  use at_energy_str_mod
  use at_boundary_str_mod
  use at_enefunc_str_mod
  use molecules_mod
  use molecules_str_mod
  use fileio_prmtop_mod
  use fileio_par_mod
  use fileio_eef1_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! P1, P2, P3, P4 in LCPO (SASA)
  real(wp), parameter :: param1 (1:4) = (/7.7887e-01_wp, -2.8063e-01_wp, -1.2968e-03_wp,  3.9328e-04_wp/)
  real(wp), parameter :: param2 (1:4) = (/5.6482e-01_wp, -1.9608e-01_wp, -1.0219e-03_wp,  2.6580e-04_wp/)
  real(wp), parameter :: param3 (1:4) = (/2.3348e-01_wp, -7.2627e-02_wp, -2.0079e-04_wp,  7.9670e-05_wp/)
  real(wp), parameter :: param4 (1:4) = (/0.0000e+00_wp,  0.0000e+00_wp,  0.0000e+00_wp,  0.0000e+00_wp/)
  real(wp), parameter :: param5 (1:4) = (/5.1245e-01_wp, -1.5966e-01_wp, -1.9781e-04_wp,  1.6392e-04_wp/)
  real(wp), parameter :: param6 (1:4) = (/7.0344e-02_wp, -1.9015e-02_wp, -2.2009e-05_wp,  1.6875e-05_wp/)
  real(wp), parameter :: param7 (1:4) = (/7.7914e-01_wp, -2.5262e-01_wp, -1.6056e-03_wp,  3.5071e-04_wp/)
  real(wp), parameter :: param8 (1:4) = (/4.9392e-01_wp, -1.6038e-01_wp, -1.5512e-04_wp,  1.6453e-04_wp/)
  real(wp), parameter :: param9 (1:4) = (/6.8563e-01_wp, -1.8680e-01_wp, -1.3557e-03_wp,  2.3743e-04_wp/) ! P3 modified by Mori
  real(wp), parameter :: param10(1:4) = (/8.8857e-01_wp, -3.3421e-01_wp, -1.8683e-03_wp,  4.9372e-04_wp/)
  real(wp), parameter :: param11(1:4) = (/7.8602e-01_wp, -2.9198e-01_wp, -6.5370e-04_wp,  3.6247e-04_wp/) ! P1 modified by Mori
  real(wp), parameter :: param12(1:4) = (/2.2599e-01_wp, -3.6648e-02_wp, -1.2297e-03_wp,  8.0038e-05_wp/)
  real(wp), parameter :: param13(1:4) = (/5.1481e-02_wp, -1.2603e-02_wp, -3.2006e-04_wp,  2.4774e-05_wp/)
  real(wp), parameter :: param14(1:4) = (/7.3511e-01_wp, -2.2116e-01_wp, -8.9148e-04_wp,  2.5230e-04_wp/)
  real(wp), parameter :: param15(1:4) = (/4.1102e-01_wp, -1.2254e-01_wp, -7.5448e-05_wp,  1.1804e-04_wp/)
  real(wp), parameter :: param16(1:4) = (/6.2577e-02_wp, -1.7874e-02_wp, -8.3120e-05_wp,  1.9849e-05_wp/)
  real(wp), parameter :: param17(1:4) = (/7.7220e-01_wp, -2.6393e-01_wp,  1.0629e-03_wp,  2.1790e-04_wp/)
  real(wp), parameter :: param18(1:4) = (/5.4581e-01_wp, -1.9477e-01_wp, -1.2873e-03_wp,  2.9247e-04_wp/)
  real(wp), parameter :: param19(1:4) = (/3.8650e-01_wp, -1.8249e-01_wp, -3.6598e-03_wp,  4.2640e-04_wp/)
  real(wp), parameter :: param20(1:4) = (/3.8730e-02_wp, -8.9339e-03_wp,  8.3582e-06_wp,  3.0381e-06_wp/)
  real(wp), parameter :: param21(1:4) = (/9.8318e-01_wp, -4.0437e-01_wp,  1.1249e-04_wp,  4.9901e-04_wp/)

  ! subroutines
  public  :: setup_enefunc_implicit_solvent_charmm
  public  :: setup_enefunc_implicit_solvent_amber
  private :: setup_enefunc_gbsa
  private :: setup_enefunc_eef1
  public  :: setup_eef1_temperature

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_implicit_solvent_charmm
  !> @brief        define implicit solvent model with charmm
  !! @authors      TM
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[in]    boundary : boundary information
  !! @param[in]    molecule : molecule information
  !! @param[in]    par      : PAR information
  !! @param[in]    eef1     : EEF1 information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_implicit_solvent_charmm(ene_info, boundary, &
                                                   molecule, par, eef1, enefunc)

    ! formal arguments
    type(s_ene_info), intent(in)    :: ene_info
    type(s_boundary), intent(in)    :: boundary
    type(s_molecule), intent(in)    :: molecule
    type(s_par),      intent(in)    :: par
    type(s_eef1),     intent(in)    :: eef1
    type(s_enefunc),  intent(inout) :: enefunc

    ! local variables
    integer                :: natoms, i, j, ano
    logical                :: found, setuperr


    if (ene_info%implicit_solvent == ImplicitSolventNONE) then
      enefunc%eef1_use = .false.
      enefunc%imm1_use = .false.
      enefunc%imic_use = .false.
      enefunc%gbsa_use = .false.
    else if (ene_info%implicit_solvent == ImplicitSolventGBSA) then
      enefunc%eef1_use = .false.
      enefunc%imm1_use = .false.
      enefunc%imic_use = .false.
      enefunc%gbsa_use = .true.
    else if (ene_info%implicit_solvent == ImplicitSolventEEF1) then
      enefunc%eef1_use = .true.
      enefunc%imm1_use = .false.
      enefunc%imic_use = .false.
      enefunc%gbsa_use = .false.
    else if (ene_info%implicit_solvent == ImplicitSolventIMM1) then
      enefunc%eef1_use = .true.
      enefunc%imm1_use = .true.
      enefunc%imic_use = .false.
      enefunc%gbsa_use = .false.
    else if (ene_info%implicit_solvent == ImplicitSolventIMIC) then
      enefunc%eef1_use = .true.
      enefunc%imm1_use = .true.
      enefunc%imic_use = .true.
      enefunc%gbsa_use = .false.
    end if

    if (enefunc%gbsa_use) then
      call setup_enefunc_gbsa(ene_info, boundary, molecule, enefunc)
    end if

    if (enefunc%eef1_use) then
      call setup_enefunc_eef1(ene_info, boundary, par, eef1, molecule, enefunc)
    end if

    return

  end subroutine setup_enefunc_implicit_solvent_charmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_implicit_solvent_amber
  !> @brief        define implicit solvent model with amber
  !! @authors      TM
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[in]    boundary : boundary information
  !! @param[in]    molecule : molecule information
  !! @param[in]    prmtop   : PRMTOP information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_implicit_solvent_amber(ene_info, boundary, &
                                                  molecule, prmtop, enefunc)

    ! formal arguments
    type(s_ene_info), intent(in)    :: ene_info
    type(s_boundary), intent(in)    :: boundary
    type(s_molecule), intent(in)    :: molecule
    type(s_prmtop),   intent(in)    :: prmtop
    type(s_enefunc),  intent(inout) :: enefunc

    ! local variables
    integer                :: natoms, i, j, ano
    logical                :: found, setuperr


    if (ene_info%implicit_solvent == ImplicitSolventNONE) then
      enefunc%eef1_use = .false.
      enefunc%imm1_use = .false.
      enefunc%imic_use = .false.
      enefunc%gbsa_use = .false.
    else if (ene_info%implicit_solvent == ImplicitSolventGBSA) then
      enefunc%eef1_use = .false.
      enefunc%imm1_use = .false.
      enefunc%imic_use = .false.
      enefunc%gbsa_use = .true.
    else
      call error_msg('Setup_Enefunc_Implicit_Solvent_Amber> The specified implicit solvent model is not available in the AMBER force field')
    end if

    if (enefunc%gbsa_use) then
      call setup_enefunc_gbsa(ene_info, boundary, molecule, enefunc, prmtop)
    end if

    return

  end subroutine setup_enefunc_implicit_solvent_amber

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_gbsa
  !> @brief        define GBSA implicit solvent model
  !! @authors      TM
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[in]    boundary : boundary information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !! @param[in]    prmtop   : PRMTOP information (optional)
  !!
  !! @note         A.Bondi, J.Phys.Chem., 68, 441-451 (1964)
  !!               J.Srinivasan et al., Theor.Chem.Acc., 101, 426-434 (1999)
  !!               J.Weiser et al., J.Comput.Chem., 20, 217-230 (1999)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_gbsa(ene_info, boundary, molecule, enefunc, prmtop)

    integer, parameter     :: NO_CONNECTED_ATOM = -999

    ! formal arguments
    type(s_ene_info),         intent(in)    :: ene_info
    type(s_boundary),         intent(in)    :: boundary
    type(s_molecule),         intent(in)    :: molecule
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_prmtop), optional, intent(in)    :: prmtop

    ! local variables
    real(wp)               :: mass
    integer                :: natoms, i, j, k, m, n
    integer                :: z, z1, z2, z3, z4
    integer                :: ano1, ano2, hb, ifound
    character(6)           :: atmcls

    integer, allocatable   :: atomno(:)


    if (boundary%type /= BoundaryTypeNOBC) then
      call error_msg('Setup_Enefunc_Gbsa> GBSA is available in NOBC currently')
    end if

    if (main_rank) then
      write(MsgOut,'(a)') 'Setup_Enefunc_Gbsa> Setup parameters for GBSA'
      write(MsgOut,'(a)') ''
    end if

    natoms                       = molecule%num_atoms
    enefunc%gbsa%cutoffdist      = ene_info%cutoffdist - 1.0_wp
    enefunc%gbsa%eps_solvent     = ene_info%gbsa_eps_solvent
    enefunc%gbsa%eps_solute      = ene_info%gbsa_eps_solute
    enefunc%gbsa%sfactor_alpha   = ene_info%gbsa_alpha
    enefunc%gbsa%sfactor_beta    = ene_info%gbsa_beta
    enefunc%gbsa%sfactor_gamma   = ene_info%gbsa_gamma
    enefunc%gbsa%salt_cons       = ene_info%gbsa_salt_cons
    enefunc%gbsa%vdw_offset      = ene_info%gbsa_vdw_offset
    enefunc%gbsa%surface_tension = ene_info%gbsa_surf_tens
    enefunc%gbsa%temperature     = ene_info%gbsa_temperature
    enefunc%gbsa%tmp_calc        = .false.


    ! count the number of non-hydrogen atoms for SA term
    !
    allocate(atomno(molecule%num_atoms))

    ifound = 0
    do i = 1, molecule%num_atoms
      mass   = molecule%mass(i)
      atmcls = molecule%atom_cls_name(i)
      atomno(i) = atomic_number(mass, atmcls)
      if (atomno(i) > 1) ifound = ifound + 1
    end do
    enefunc%gbsa%num_sasa_atoms = ifound

    call alloc_enefunc(enefunc, EneFuncGbsa, natoms, ifound)


    ! setup GB term
    !
    if (present(prmtop)) then ! AMBER
      do i = 1, molecule%num_atoms
        enefunc%gbsa%vdw_radius  (i) = prmtop%radi_born(i)
        enefunc%gbsa%scale_factor(i) = prmtop%fs_born(i)
      end do
    else                      ! CHARMM
      do i = 1, molecule%num_atoms
        z = atomno(i)
        if      (z == 1)  then                  ! H
          enefunc%gbsa%vdw_radius  (i) = 1.20_wp
          enefunc%gbsa%scale_factor(i) = 0.85_wp
        else if (z == 6)  then                  ! C
          enefunc%gbsa%vdw_radius  (i) = 1.70_wp
          enefunc%gbsa%scale_factor(i) = 0.72_wp
        else if (z == 7)  then                  ! N
          enefunc%gbsa%vdw_radius  (i) = 1.55_wp
          enefunc%gbsa%scale_factor(i) = 0.79_wp
        else if (z == 8)  then                  ! O
          enefunc%gbsa%vdw_radius  (i) = 1.50_wp
          enefunc%gbsa%scale_factor(i) = 0.85_wp
        else if (z == 11)  then                 ! Na
          enefunc%gbsa%vdw_radius  (i) = 2.27_wp
          enefunc%gbsa%scale_factor(i) = 0.80_wp
        else if (z == 12)  then                 ! Mg
          enefunc%gbsa%vdw_radius  (i) = 1.73_wp
          enefunc%gbsa%scale_factor(i) = 0.80_wp
        else if (z == 15)  then                 ! P
          enefunc%gbsa%vdw_radius  (i) = 1.85_wp
          enefunc%gbsa%scale_factor(i) = 0.86_wp
        else if (z == 16)  then                 ! S
          enefunc%gbsa%vdw_radius  (i) = 1.80_wp
          enefunc%gbsa%scale_factor(i) = 0.96_wp
        else if (z == 17)  then                 ! Cl
          enefunc%gbsa%vdw_radius  (i) = 1.75_wp
          enefunc%gbsa%scale_factor(i) = 0.80_wp
        else if (z == 19)  then                 ! K
          enefunc%gbsa%vdw_radius  (i) = 2.75_wp
          enefunc%gbsa%scale_factor(i) = 0.80_wp
        else if (z == 30)  then                 ! Zn
          enefunc%gbsa%vdw_radius  (i) = 1.39_wp
          enefunc%gbsa%scale_factor(i) = 0.80_wp
        else                                    ! others
          write(MsgOut,'(a,a4)') 'Setup_Enefunc_Gbsa> WARNING in GB parameter setting! unknown atom:', molecule%atom_name(i)
          write(MsgOut,'(a)')    '                    General parameters were applied'
          write(MsgOut,'(a)')    ' '
          enefunc%gbsa%vdw_radius  (i) = 1.50_wp
          enefunc%gbsa%scale_factor(i) = 0.80_wp
        end if
      end do
    end if


    ! setup SA term
    !
    ifound = 0

    do i = 1, molecule%num_atoms

      z = atomno(i)

      if (z <= 1) cycle

      ifound = ifound + 1
      enefunc%gbsa%sasa_atom_list(ifound) = i

      ! first bonded neighbor atom
      z1 = NO_CONNECTED_ATOM
      do j = 1, molecule%num_bonds
        ano1 = molecule%bond_list(1,j)
        ano2 = molecule%bond_list(2,j)
        if (ano1 == i) then
          z1 = atomno(ano2)
          exit
        else if (ano2 == i) then
          z1 = atomno(ano1)
          exit
        end if
      end do

      ! second bonded neighbor atom
      z2 = NO_CONNECTED_ATOM
      do k = j+1, molecule%num_bonds
        ano1 = molecule%bond_list(1,k)
        ano2 = molecule%bond_list(2,k)
        if (ano1 == i) then
          z2 = atomno(ano2)
          exit
        else if (ano2 == i) then
          z2 = atomno(ano1)
          exit
        end if
      end do

      ! third bonded neighbor atom
      z3 = NO_CONNECTED_ATOM
      do m = k+1, molecule%num_bonds
        ano1 = molecule%bond_list(1,m)
        ano2 = molecule%bond_list(2,m)
        if (ano1 == i) then
          z3 = atomno(ano2)
          exit
        else if (ano2 == i) then
          z3 = atomno(ano1)
          exit
        end if
      end do

      ! fourth bonded neighbor atom
      z4 = NO_CONNECTED_ATOM
      do n = m+1, molecule%num_bonds
        ano1 = molecule%bond_list(1,n)
        ano2 = molecule%bond_list(2,n)
        if (ano1 == i) then
          z4 = atomno(ano2)
          exit
        else if (ano2 == i) then
          z4 = atomno(ano1)
          exit
        end if
      end do

      ! count the total number of H-bonds
      !
      hb = 0
      if (z1 == 1) hb = hb + 1
      if (z2 == 1) hb = hb + 1
      if (z3 == 1) hb = hb + 1
      if (z4 == 1) hb = hb + 1

      ! check atom type
      !
      enefunc%gbsa%sasa_atom_type(i) = 0

      if (z == 6) then  ! C
        if (z4 == NO_CONNECTED_ATOM) then
          ! sp2
          if (hb == 1) enefunc%gbsa%sasa_atom_type(i) = 5
          if (hb == 0) enefunc%gbsa%sasa_atom_type(i) = 6
        else
          ! sp3
          if (hb == 3) enefunc%gbsa%sasa_atom_type(i) = 1
          if (hb == 2) enefunc%gbsa%sasa_atom_type(i) = 2
          if (hb == 1) enefunc%gbsa%sasa_atom_type(i) = 3
          if (hb == 0) enefunc%gbsa%sasa_atom_type(i) = 4
        end if
        enefunc%gbsa%sasa_vdw_radius(i) = 1.70_wp

      else if (z == 8) then ! O
        if (z2 == NO_CONNECTED_ATOM) then
          ! sp2
          if (present(prmtop)) then ! AMBER
            if (prmtop%atom_cls_name(i) == 'O') then
              enefunc%gbsa%sasa_atom_type(i) = 9
            else
              enefunc%gbsa%sasa_atom_type(i) = 10
            end if
          else                      ! CHARMM
            if (molecule%atom_name(i) == 'O') then
              enefunc%gbsa%sasa_atom_type(i) = 9
            else
              enefunc%gbsa%sasa_atom_type(i) = 10
            end if
          end if
        else
          ! sp3 (hb = 2 for TIP3)
          if (hb >= 1) enefunc%gbsa%sasa_atom_type(i) = 7
          if (hb == 0) enefunc%gbsa%sasa_atom_type(i) = 8
        end if
        enefunc%gbsa%sasa_vdw_radius(i) = 1.60_wp

      else if (z == 7) then ! N
        if (z4 == NO_CONNECTED_ATOM) then
          ! sp2
          if (hb == 2) enefunc%gbsa%sasa_atom_type(i) = 14
          if (hb == 1) enefunc%gbsa%sasa_atom_type(i) = 15
          if (hb == 0) enefunc%gbsa%sasa_atom_type(i) = 16
        else
          ! sp3
          if (hb == 3) enefunc%gbsa%sasa_atom_type(i) = 11
          if (hb == 2) enefunc%gbsa%sasa_atom_type(i) = 12
          if (hb == 1) enefunc%gbsa%sasa_atom_type(i) = 13
        end if
        enefunc%gbsa%sasa_vdw_radius(i) = 1.65_wp

      else if (z == 16) then ! S
        if (hb == 1) enefunc%gbsa%sasa_atom_type(i) = 17
        if (hb == 0) enefunc%gbsa%sasa_atom_type(i) = 18
        enefunc%gbsa%sasa_vdw_radius(i) = 1.90_wp

      else if (z == 15) then ! P
        if (hb == 1) enefunc%gbsa%sasa_atom_type(i) = 19
        if (hb == 0) enefunc%gbsa%sasa_atom_type(i) = 20
        enefunc%gbsa%sasa_vdw_radius(i) = 1.90_wp

      else if (z == 17) then ! Cl
        enefunc%gbsa%sasa_atom_type (i) = 21
        enefunc%gbsa%sasa_vdw_radius(i) = 1.80_wp

      else ! others
        if      (z == 12) then  ! MG
          enefunc%gbsa%sasa_atom_type (i) = 8
          enefunc%gbsa%sasa_vdw_radius(i) = 1.185_wp
        else if (z == 30) then  ! ZN
          enefunc%gbsa%sasa_atom_type (i) = 8
          enefunc%gbsa%sasa_vdw_radius(i) = 1.09_wp
        else                    ! others
          write(MsgOut,'(a,a4)') 'Setup_Enefunc_Gbsa> WARNING in SA parameter setting! unknown atom:', molecule%atom_name(i)
          write(MsgOut,'(a)')    '                    General parameters were applied'
          write(MsgOut,'(a)')    ' '
          enefunc%gbsa%sasa_atom_type (i) = 8
          enefunc%gbsa%sasa_vdw_radius(i) = 1.2_wp
        end if

      end if

      if (enefunc%gbsa%sasa_atom_type(i) == 0) then
        write(MsgOut,'(a,a4)') 'Setup_Enefunc_Gbsa> Could not define the LCPO parameters for atom:', molecule%atom_name(i)
        call error_msg('Setup_Enefunc_Gbsa> Could not define the LCPO parameters for some atoms.')
      end if

    end do

    ! add solvent probe radius
    enefunc%gbsa%sasa_vdw_radius(:) = enefunc%gbsa%sasa_vdw_radius(:) + 1.40_wp

    enefunc%gbsa%pairlist_distance  = maxval(enefunc%gbsa%sasa_vdw_radius(:))*2.0_wp

    call get_loop_index(natoms, enefunc%gbsa%istart_gb, enefunc%gbsa%iend_gb)

    enefunc%gbsa%sasa_parameters( 1,1:4) = param1 (1:4)
    enefunc%gbsa%sasa_parameters( 2,1:4) = param2 (1:4)
    enefunc%gbsa%sasa_parameters( 3,1:4) = param3 (1:4)
    enefunc%gbsa%sasa_parameters( 4,1:4) = param4 (1:4)
    enefunc%gbsa%sasa_parameters( 5,1:4) = param5 (1:4)
    enefunc%gbsa%sasa_parameters( 6,1:4) = param6 (1:4)
    enefunc%gbsa%sasa_parameters( 7,1:4) = param7 (1:4)
    enefunc%gbsa%sasa_parameters( 8,1:4) = param8 (1:4)
    enefunc%gbsa%sasa_parameters( 9,1:4) = param9 (1:4)
    enefunc%gbsa%sasa_parameters(10,1:4) = param10(1:4)
    enefunc%gbsa%sasa_parameters(11,1:4) = param11(1:4)
    enefunc%gbsa%sasa_parameters(12,1:4) = param12(1:4)
    enefunc%gbsa%sasa_parameters(13,1:4) = param13(1:4)
    enefunc%gbsa%sasa_parameters(14,1:4) = param14(1:4)
    enefunc%gbsa%sasa_parameters(15,1:4) = param15(1:4)
    enefunc%gbsa%sasa_parameters(16,1:4) = param16(1:4)
    enefunc%gbsa%sasa_parameters(17,1:4) = param17(1:4)
    enefunc%gbsa%sasa_parameters(18,1:4) = param18(1:4)
    enefunc%gbsa%sasa_parameters(19,1:4) = param19(1:4)
    enefunc%gbsa%sasa_parameters(20,1:4) = param20(1:4)
    enefunc%gbsa%sasa_parameters(21,1:4) = param21(1:4)

    ! deallocate local array
    !
    deallocate(atomno)

    return

  end subroutine setup_enefunc_gbsa

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_eef1
  !> @brief        define EEF1 implicit solvent model
  !! @authors      TM
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[in]    boundary : boundary information
  !! @param[in]    par      : PAR information
  !! @param[in]    eef1     : EEF1 information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_eef1(ene_info, boundary, par, eef1, molecule, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_boundary),        intent(in)    :: boundary
    type(s_par),             intent(in)    :: par
    type(s_eef1),            intent(in)    :: eef1
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: natoms, i, j, ano
    logical                  :: found, setuperr
    character(6)             :: aname


    if (boundary%type /= BoundaryTypeNOBC) then
      call error_msg('Setup_Enefunc_Implicit_Solvent> EEF1/IMM1/IMIC is available in NOBC currently')
    end if

    if (main_rank) then
      write(MsgOut,'(a)') 'Setup_Enefunc_Eef1> Setup parameters for EEF1/IMM1/IMIC'
      write(MsgOut,'(a)') '  Lookup table was disabled'
      write(MsgOut,'(a)') ''
    end if

    call dealloc_enefunc(enefunc, EneFuncTable)
    call dealloc_enefunc(enefunc, EneFuncTableWater)
    call dealloc_enefunc(enefunc, EneFuncTableSolute)

    natoms = maxval(molecule%atom_cls_no(:))
    enefunc%eef1%num_atoms = natoms
    enefunc%eef1%imm1_memb_thick  = ene_info%imm1_memb_thick
    enefunc%eef1%imm1_exponent_n  = ene_info%imm1_exponent_n
    enefunc%eef1%imm1_factor_a    = ene_info%imm1_factor_a
    enefunc%eef1%imm1_make_pore   = ene_info%imm1_make_pore
    enefunc%eef1%imm1_pore_radius = ene_info%imm1_pore_radius
    enefunc%eef1%imic_axis_a      = ene_info%imic_axis_a
    enefunc%eef1%imic_axis_b      = ene_info%imic_axis_b
    enefunc%eef1%imic_axis_c      = ene_info%imic_axis_c
    enefunc%eef1%imic_exponent_m1 = ene_info%imic_exponent_m1
    enefunc%eef1%imic_exponent_m2 = ene_info%imic_exponent_m2
    enefunc%eef1%imic_steepness   = ene_info%imic_steepness

    call alloc_enefunc(enefunc, EneFuncEef1, natoms)

    setuperr = .false.
    do i = 1, molecule%num_atoms
      found = .false.
      do j = 1, eef1%num_atoms
        if (molecule%atom_cls_name(i) == eef1%atom_name(j)) then
          ano = molecule%atom_cls_no(i)
          enefunc%eef1%atom_name (    ano) = eef1%atom_name(j)
          enefunc%eef1%volume    (1:2,ano) = eef1%volume(1:2,j)
          enefunc%eef1%gref_0    (1:2,ano) = eef1%gref  (1:2,j)
          enefunc%eef1%gfree_0   (1:2,ano) = eef1%gfree (1:2,j)
          enefunc%eef1%href_0    (1:2,ano) = eef1%href  (1:2,j)
          enefunc%eef1%cpref_0   (1:2,ano) = eef1%cpref (1:2,j)
          enefunc%eef1%inv_lambda(1:2,ano) = 1.0_wp/eef1%sigw(1:2,j)
          enefunc%eef1%rvdw      (    ano) = par%nonb_rmin(ano)
          found = .true.
          exit
        end if
      end do

      aname = molecule%atom_cls_name(i)
      if (.not. found .and. aname(1:1) == 'H') then
        ano = molecule%atom_cls_no(i)
        enefunc%eef1%atom_name (    ano) = molecule%atom_cls_name(i)
        enefunc%eef1%volume    (1:2,ano) = 0.0_wp
        enefunc%eef1%gref_0    (1:2,ano) = 0.0_wp
        enefunc%eef1%gfree_0   (1:2,ano) = 0.0_wp
        enefunc%eef1%href_0    (1:2,ano) = 0.0_wp
        enefunc%eef1%cpref_0   (1:2,ano) = 0.0_wp
        enefunc%eef1%inv_lambda(1:2,ano) = 1.0_wp/3.5_wp
        enefunc%eef1%rvdw      (    ano) = par%nonb_rmin(ano)
        found = .true.
      end if

      if (.not. found) then
        write(MsgOut,*) 'Setup_Enefunc_Implicit_Solvent> not found EEF1 parameter:'
        write(MsgOut,*) '[',molecule%atom_cls_name(i),'] in eef1file. (ERROR)'
        setuperr = .true.
      end if

    end do

    if (setuperr) call error_msg('Setup_Enefunc_Implicit_Solvent> Setup error in EEF1/IMM1')

    call get_loop_index(molecule%num_atoms, enefunc%eef1%istart, enefunc%eef1%iend)

    call setup_eef1_temperature(ene_info%eef1_temperature, enefunc)

    return

  end subroutine setup_enefunc_eef1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_eef1_temperature
  !> @brief        get EEF1 parameters at any temperature
  !! @authors      TM
  !! @param[in]    temperature : temperature for EEF1 parameters
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_eef1_temperature(temperature, enefunc)

    ! formal arguments
    real(wp),                intent(in)    :: temperature
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer           :: natoms, i, j
    real(wp)          :: dT, lnTT0, dSref, ratio, gref

    real(wp), pointer :: gfree_0(:,:), gref_0(:,:), href_0(:,:), cpref_0(:,:)
    real(wp), pointer :: gfree_t(:,:), gref_t(:,:), alpha_4pi(:,:)
    real(wp), pointer :: inv_lambda(:,:)


    ! parameters for 298.15 K
    gfree_0    => enefunc%eef1%gfree_0
    gref_0     => enefunc%eef1%gref_0
    href_0     => enefunc%eef1%href_0
    cpref_0    => enefunc%eef1%cpref_0
    inv_lambda => enefunc%eef1%inv_lambda

    ! parameters for temperature T
    gfree_t    => enefunc%eef1%gfree_t
    gref_t     => enefunc%eef1%gref_t
    alpha_4pi  => enefunc%eef1%alpha_4pi


    dT     = temperature - 298.15_wp
    lnTT0  = log(temperature/298.15_wp)

    do i = 1, enefunc%eef1%num_atoms
      do j = 1, 2

        dSref = (href_0(j,i) - gref_0(j,i))/298.15_wp
        gref  = gref_0 (j,i) - dSref*dT &
              - cpref_0(j,i)*(temperature*lnTT0-dT)/1000.0_wp

        if (abs(gref_0(j,i) - 0.0_wp) < EPS) then
          ratio = 0.0_wp
        else
          ratio = gfree_0(j,i)/gref_0(j,i)
        end if

        gref_t (j,i) = gref
        gfree_t(j,i) = gref*ratio

        if (.not. enefunc%imm1_use) then
          alpha_4pi(j,i) = 2.0_wp*gfree_t(j,i)*inv_lambda(j,i)/(4.0_wp*PI*SQRT(PI))
        else if (enefunc%imm1_use) then
          alpha_4pi(j,i) = 2.0_wp*inv_lambda(j,i)/(4.0_wp*PI*SQRT(PI))
        end if

      end do
    end do

    return

  end subroutine setup_eef1_temperature

end module at_enefunc_gbsa_mod
