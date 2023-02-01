!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_constraints_mod
!> @brief   constraints module
!! @authors Jaewoon Jung (JJ), Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_constraints_mod

  use sp_constraints_str_mod
  use sp_dynamics_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use molecules_mod
  use molecules_str_mod
  use fileio_grotop_mod
  use fileio_prmtop_mod
  use fileio_par_mod
  use fileio_control_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
#ifdef HAVE_MPI_GENESIS
#ifdef MSMPI
!GCC$ ATTRIBUTES DLLIMPORT :: MPI_BOTTOM, MPI_IN_PLACE
#endif
#endif
  private

  ! structures
  type, public :: s_cons_info
    logical                  :: rigid_bond      = .false.
    logical                  :: fast_bond       = .false.
    logical                  :: fast_water      = .false.
    integer                  :: hydrogen_type   = ConstraintAtomName
    integer                  :: shake_iteration = 500
    real(wp)                 :: shake_tolerance = 1.0e-10_wip
    integer                  :: lincs_iteration = 1
    integer                  :: lincs_order     = 4
    character(5)             :: water_model     = 'TIP3'
    real(wp)                 :: hydrogen_mass_upper_bound = 2.1_wip
  end type s_cons_info

  ! parameters
  integer, public, parameter :: ConstraintModeLEAP  = 1
  integer, public, parameter :: ConstraintModeVVER1 = 2
  integer, public, parameter :: ConstraintModeVVER2 = 3

  ! subroutines
  public  :: show_ctrl_constraints
  public  :: read_ctrl_constraints
  public  :: setup_constraints
  public  :: setup_constraints_pio
  public  :: compute_constraints
  public  :: update_vel_group
  public  :: update_vel_group_3d
  public  :: compute_kin_group
  public  :: compute_group_deg_freedom
  public  :: compute_virial_group
  public  :: compute_vv1_group
  public  :: compute_vv1_coord_group
  public  :: update_vel_vv1_group
  public  :: update_coord_vv1_group
  public  :: update_coord_vv2_group
  public  :: compute_vv2_group
  public  :: water_force_redistribution
  public  :: decide_dummy
  public  :: compute_settle_min
  public  :: setup_mass_repartitioning
  public  :: setup_mass_repartitioning_back
  private :: setup_fast_water
  private :: setup_fast_water_tip4
  private :: setup_fast_water_pio
  private :: setup_rigid_bond
  private :: setup_rigid_bond_pio
  private :: compute_settle
  private :: compute_shake
  private :: compute_rattle_fast_vv1
  private :: compute_rattle_vv1
  private :: compute_rattle_fast_vv2
  private :: compute_rattle_vv2

contains
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_constraints
  !> @brief        show CONSTRAINTS section usage
  !! @authors      JJ
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "md", "min"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_constraints(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode


    if (show_all) then

      select case (run_mode)

      case ('md')

        write(MsgOut,'(A)') '[CONSTRAINTS]'
        write(MsgOut,'(A)') 'rigid_bond    = NO         # constraints all bonds involving hydrogen'
        write(MsgOut,'(A)') '# shake_iteration = 500      # max number of SHAKE/RATTLE iterations'
        write(MsgOut,'(A)') '# shake_tolerance = 1.0e-8   # SHAKE/RATTLE tolerance (Ang)'
        write(MsgOut,'(A)') '# water_model     = TIP3    # water model'
        write(MsgOut,'(A)') '# hydrogen_mass_upper_bound  = 2.1    # upper mass limit to define the hydrogen atom'
        write(MsgOut,'(A)') ' '

      end select

    else

      select case (run_mode)

      case ('md')

        write(MsgOut,'(A)') '[CONSTRAINTS]'
        write(MsgOut,'(A)') 'rigid_bond    = NO         # constraints all bonds involving hydrogen'
        write(MsgOut,'(A)') ' '

      end select

    end if

    return

  end subroutine show_ctrl_constraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_constraint
  !> @brief        read CONSTRAINTS section in the control file
  !! @authors      JJ
  !! @param[in]    handle    :unit number
  !! @param[out]   cons_info :CONSTRAINTS section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_constraints(handle, cons_info)
  
    ! parameters
    character(*),            parameter     :: Section = 'Constraints'

    ! formal arguments
    integer,                 intent(in)    :: handle
    type(s_cons_info),       intent(inout) :: cons_info
  

    ! read parameters from control file
    !
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_logical(handle, Section, 'rigid_bond',      &
                               cons_info%rigid_bond)
    call read_ctrlfile_logical(handle, Section, 'fast_water',      &
                               cons_info%fast_water)
    call read_ctrlfile_integer(handle, Section, 'shake_iteration', &
                               cons_info%shake_iteration)
    call read_ctrlfile_real   (handle, Section, 'shake_tolerance', &
                               cons_info%shake_tolerance)
    call read_ctrlfile_string (handle, Section, 'water_model',     &
                               cons_info%water_model)
    call read_ctrlfile_type   (handle, Section, 'hydrogen_type',   &
                               cons_info%hydrogen_type, ConstraintAtomType)
    call read_ctrlfile_real   (handle, Section, 'hydrogen_mass_upper_bound', &
                               cons_info%hydrogen_mass_upper_bound)

    call end_ctrlfile_section(handle)

    LIGHT_ATOM_MASS_LIMIT = cons_info%hydrogen_mass_upper_bound

    ! fast_water is true if rigid_bond is true
    !
    if (cons_info%rigid_bond) cons_info%fast_water = .true.

    ! write parameters to MsgOut
    !
    if (main_rank) then

      write(MsgOut,'(A)') 'Read_Ctrl_Constraints> Parameters for Constraints'

      if (cons_info%rigid_bond) then
        write(MsgOut,'(A30)')                                    &
              '  rigid_bond      =        yes'
        write(MsgOut,'(A20,I10,A20,E10.3)')                      &
              '  shake_iteration = ', cons_info%shake_iteration, &
              '  shake_tolerance = ', cons_info%shake_tolerance

        if (cons_info%fast_water) then
          write(MsgOut,'(A30,A20,A10)')                          &
                '  fast_water      =        yes',                &
                '  water_model     = ', trim(cons_info%water_model)
        else
          write(MsgOut,'(A30)')                                  &
                '  fast_water      =         no'
        end if

        if (cons_info%hydrogen_type==ConstraintAtomName) then
          write(MsgOut,'(A30)')                                    &
                '  hydrogen_type   =       name'
        else if (cons_info%hydrogen_type==ConstraintAtomMass) then
          write(MsgOut,'(A30)')                                    &
                '  hydrogen_type   =       mass'
        else 
          write(MsgOut,'(A30)')                                    &
                '  hydrogen_type   =  name|mass'
        end if

      else
        write(MsgOut,'(A30)') '  rigid_bond      =         no'
      end if

      write(MsgOut,'(A)') ' '

    end if

    return    

  end subroutine read_ctrl_constraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_constraints
  !> @brief        setup constraints for domain decomposition
  !! @authors      JJ
  !! @param[in]    cons_info :CONSTRAINTS section control parameters information
  !! @param[in]    par         : PAR information
  !! @param[in]    prmtop      : AMBER parameter topology information
  !! @param[in]    grotop      : GROMACS parameter topology information
  !! @param[inout] molecule    : molecule information
  !! @param[in]    domain      : domain information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_constraints(cons_info, par, prmtop, grotop, molecule, &
                               domain, enefunc, constraints)

    ! formal arguments
    type(s_cons_info),       intent(in)    :: cons_info
    type(s_par),             intent(in)    :: par
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(inout) :: molecule
    type(s_domain),          intent(in)    :: domain  
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_constraints),     intent(inout) :: constraints


    constraints%rigid_bond      = cons_info%rigid_bond
    constraints%fast_bond       = cons_info%fast_bond
    constraints%fast_water      = cons_info%fast_water
    constraints%shake_iteration = cons_info%shake_iteration
    constraints%shake_tolerance = cons_info%shake_tolerance
    constraints%lincs_iteration = cons_info%lincs_iteration
    constraints%lincs_order     = cons_info%lincs_order
    constraints%water_model     = cons_info%water_model
    constraints%hydrogen_type   = cons_info%hydrogen_type

    if (constraints%rigid_bond) then

      if (constraints%hydrogen_type == ConstraintAtomName  &
         .and. molecule%special_hydrogen) &
      call error_msg('Setup_Constraints> Non ordinary hydrogen name is not'//&
                         ' allowed. If you want, use hydrogen_type option.')

      ! setup SETTLE
      !
      if (constraints%fast_water .and. enefunc%table%num_water > 0) then

        if (constraints%water_type == TIP4) then
          call setup_fast_water_tip4(par, prmtop, grotop, &
                                     molecule, enefunc, constraints)
        else
          call setup_fast_water(par, prmtop, grotop, &
                                molecule, enefunc, constraints)
        end if

        ! update number of degrees of freedom
        !
        if (constraints%water_type == TIP4) then
          call update_num_deg_freedom('After setup of SETTLE',    &
                                      -6*enefunc%table%num_water  &
                                        *domain%num_duplicate,    &
                                      molecule%num_deg_freedom)
        else if (constraints%water_type == TIP3) then
          call update_num_deg_freedom('After setup of SETTLE',    &
                                      -3*enefunc%table%num_water  &
                                        *domain%num_duplicate,    &
                                      molecule%num_deg_freedom)
        else if (constraints%water_type == TIP2) then
          call update_num_deg_freedom('After setup of Water constraint',    &
                                      -enefunc%table%num_water  &
                                        *domain%num_duplicate,    &
                                      molecule%num_deg_freedom)
        end if

      end if

      ! setup SHAKE and RATTLE
      !
      call setup_rigid_bond(constraints)

      ! update number of degrees of freedom
      !
      if (constraints%num_bonds > 0) then
        call update_num_deg_freedom('After setup of SHAKE/RATTLE', &
                                    -constraints%num_bonds,        &
                                    molecule%num_deg_freedom)

      end if

    else if (constraints%fast_water .and. &
             constraints%water_type == TIP3 .and. &
             enefunc%table%num_water > 0) then

      call setup_fast_water(par, prmtop, grotop, &
                            molecule, enefunc, constraints)

    else if (constraints%water_type == TIP4 .and. &
             enefunc%table%num_water > 0) then

      if (main_rank) &
      write(MsgOut,'(A)') 'fast_water is defined as yes for TIP4 water model case'
      call update_num_deg_freedom('After setup of SETTLE',    &
                                  -6*enefunc%table%num_water  &
                                    *domain%num_duplicate,    &
                                  molecule%num_deg_freedom)
      call setup_fast_water_tip4(par, prmtop, grotop, &
                                 molecule, enefunc, constraints)
 
    end if
   
    return

  end subroutine setup_constraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_constraints_pio
  !> @brief        setup constraints for domain decomposition
  !! @authors      JJ
  !! @param[in]    cons_info :CONSTRAINTS section control parameters information
  !! @param[in]    pio_restart : flag for parallel I/O restart
  !! @param[inout  enefunc     : potential energy functions information
  !! @param[inout] constraints : constraints information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_constraints_pio(cons_info, pio_restart, enefunc,  &
                                   constraints, domain)

    ! formal arguments
    type(s_cons_info),       intent(in)    :: cons_info
    logical,                 intent(in)    :: pio_restart
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_constraints),     intent(inout) :: constraints
    type(s_domain),          intent(inout) :: domain


    constraints%rigid_bond      = cons_info%rigid_bond
    constraints%fast_bond       = cons_info%fast_bond
    constraints%fast_water      = cons_info%fast_water
    constraints%shake_iteration = cons_info%shake_iteration
    constraints%shake_tolerance = cons_info%shake_tolerance
    constraints%lincs_iteration = cons_info%lincs_iteration
    constraints%lincs_order     = cons_info%lincs_order
    constraints%water_model     = cons_info%water_model

    if (constraints%rigid_bond) then

      ! setup SETTLE
      !
      if (constraints%fast_water) then

        call setup_fast_water_pio(constraints)

        ! update number of degrees of freedom
        !
        call update_num_deg_freedom('After setup of SETTLE',    &
                                    -3*enefunc%table%num_water, &
                                    domain%num_deg_freedom)

      end if

      ! setup SHAKE and RATTLE
      !
      call setup_rigid_bond_pio(constraints)

      ! update number of degrees of freedom
      !
      if (constraints%num_bonds > 0) then
        call update_num_deg_freedom('After setup of SHAKE/RATTLE', &
                                    -constraints%num_bonds,        &
                                    domain%num_deg_freedom)

      end if

    else if (constraints%water_type == TIP4) then

      if (main_rank) &
      write(MsgOut,'(A)') 'fast_water is defined as yes for TIP4 water model case'
      call setup_fast_water_pio(constraints)
      call update_num_deg_freedom('After setup of SETTLE',    &
                                  -6*enefunc%table%num_water, &
                                  domain%num_deg_freedom)
    
    end if
   
    return

  end subroutine setup_constraints_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_constraints
  !> @brief        update coordinates according to constraints
  !! @authors      JJ
  !! @param[in]    cons_mode   : constraint mode
  !! @param[in]    vel_update  : flag for update velocity or not
  !! @param[in]    dt          : time step
  !! @param[inout] coord_ref   : reference coordinates
  !! @param[in]    domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] coord       : coordinates
  !! @param[inout] vel         : velocities
  !! @param[inout] virial      : virial of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_constraints(cons_mode, vel_update, dt, coord_ref, &
                                 domain, constraints, coord, vel, virial)

    ! formal arguments
    integer,                 intent(in)    :: cons_mode
    logical,                 intent(in)    :: vel_update
    real(wip),               intent(in)    :: dt
    real(wip),               intent(inout) :: coord_ref(:,:,:)
    type(s_domain),          intent(in)    :: domain
    type(s_constraints),     intent(inout) :: constraints
    real(wip),               intent(inout) :: coord(:,:,:)
    real(wip),               intent(inout) :: vel(:,:,:)
    real(dp ),               intent(inout) :: virial(:,:)


    call timer(TimerConstraint, TimerOn)

    ! constraints
    !
    select case (cons_mode)

    case (ConstraintModeLEAP)

      virial(1:3,1:3) = 0.0_dp 

      if (constraints%water_type >= TIP3) then
        if (constraints%num_water > 0) &
        call compute_settle(vel_update, dt, coord_ref, domain, &
                            constraints, coord, vel, virial)
        if (constraints%water_type == TIP4) &
          call decide_dummy(domain, constraints, coord)

      else if (constraints%water_type == TIP2) then
        call compute_settle_tip2(vel_update, dt, coord_ref, domain, &
                            constraints, coord, vel, virial)
      end if
      if (constraints%water_type /= TIP2) &
      call compute_shake (vel_update, dt, coord_ref, domain, &
                          constraints, coord, vel, virial)

      virial(1:3,1:3) = virial(1:3,1:3)/(dt*dt)

    case (ConstraintModeVVER1)

      call compute_rattle_fast_vv1(dt, coord_ref, &
                                   domain, constraints, coord, vel)

      call compute_rattle_vv1     (dt, coord_ref, &
                                   domain, constraints, coord, vel)

      if (constraints%water_type == TIP4) &
        call decide_dummy(domain, constraints, coord)

    case (ConstraintModeVVER2)

      if (constraints%water_type >= TIP3 .and. constraints%num_water > 0) then
        call compute_rattle_fast_vv2(domain, constraints, coord, vel)
      else if (constraints%water_type == TIP2) then
        call compute_settle_TIP2_vv2(domain, constraints, coord, vel)
      end if
      if (constraints%water_type /= TIP2) &
      call compute_rattle_vv2     (domain, constraints, coord, vel)

    end select

    call timer(TimerConstraint, TimerOff)

    return

  end subroutine compute_constraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_fast_water
  !> @brief        setup parameter for bond in water
  !! @authors      JJ
  !! @param[in]    par         : PAR information
  !! @param[in]    prmtop      : AMBER parameter topology information
  !! @param[in]    grotop      : GROMACS parameter topology information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_fast_water(par, prmtop, grotop, &
                              molecule, enefunc, constraints)

    ! formal arguments
    type(s_par),             intent(in)    :: par
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_constraints),     intent(inout) :: constraints

    ! local variables
    integer                  :: i, j, k, i1, i2, list(3), ioffset
    character(6)             :: ci1, ci2, ci3

    type(s_grotop_mol), pointer :: gromol


    if (constraints%water_type == TIP2) then

      constraints%water_mass1 = enefunc%table%mass_1
      constraints%water_mass2 = enefunc%table%mass_2
      list(1) = enefunc%table%water_list(1,1)
      list(2) = enefunc%table%water_list(2,1)

      if (par%num_bonds > 0) then
        ci1 = molecule%atom_cls_name(list(1))
        ci2 = molecule%atom_cls_name(list(2))
        do i = 1, par%num_bonds
          if (ci1 .eq. par%bond_atom_cls(1,i) .and. &
              ci2 .eq. par%bond_atom_cls(2,i) .or.  &
              ci1 .eq. par%bond_atom_cls(2,i) .and. &
              ci2 .eq. par%bond_atom_cls(1,i) ) then
            constraints%water_r12 = par%bond_dist_min(i)
            exit
          end if
        end do
      end if

    else

    ! mass
    !
    constraints%water_massO = enefunc%table%mass_O
    constraints%water_massH = enefunc%table%mass_H

    ! min distance
    !
    list(1:3) = enefunc%table%water_list(1:3,1)

    ! charmm
    if (par%num_bonds > 0) then

      ci1 = molecule%atom_cls_name(list(1))
      ci2 = molecule%atom_cls_name(list(2))
      ci3 = molecule%atom_cls_name(list(3))

      do i = 1, par%num_bonds
        if (ci1 .eq. par%bond_atom_cls(1,i) .and. &
            ci2 .eq. par%bond_atom_cls(2,i) .or.  &
            ci1 .eq. par%bond_atom_cls(2,i) .and. &
            ci2 .eq. par%bond_atom_cls(1,i) ) then

          constraints%water_rOH = par%bond_dist_min(i)
          exit

        end if
      end do

      do i = 1, par%num_bonds
        if (ci2 .eq. par%bond_atom_cls(1,i) .and. &
            ci3 .eq. par%bond_atom_cls(2,i)) then

          constraints%water_rHH = par%bond_dist_min(i)
          exit

        end if
      end do

    ! amber
    else if (prmtop%num_atoms > 0) then

      do i = 1, prmtop%num_bondh

        i1 = prmtop%bond_inc_hy(1,i) / 3 + 1
        i2 = prmtop%bond_inc_hy(2,i) / 3 + 1

        if (list(1) == i1 .and. list(2) == i2 .or.  &
            list(2) == i1 .and. list(1) == i2) then

          constraints%water_rOH = &
               prmtop%bond_equil_uniq(prmtop%bond_inc_hy(3,i))
          exit

        end if

      end do

      do i = 1, prmtop%num_bondh

        i1 = prmtop%bond_inc_hy(1,i) / 3 + 1
        i2 = prmtop%bond_inc_hy(2,i) / 3 + 1

        if (list(2) == i1 .and. list(3) == i2 .or. &
            list(3) == i1 .and. list(2) == i2) then

          constraints%water_rHH = &
               prmtop%bond_equil_uniq(prmtop%bond_inc_hy(3,i))
          exit

        end if

      end do

    ! gromacs
    else if (grotop%num_atomtypes > 0) then

      ioffset = 0

      do i = 1, grotop%num_molss
        gromol => grotop%molss(i)%moltype%mol

        if (gromol%settles%func == 0) then

          do j = 1, grotop%molss(i)%count

            do k = 1, gromol%num_bonds

              i1 = gromol%bonds(k)%atom_idx1 + ioffset
              i2 = gromol%bonds(k)%atom_idx2 + ioffset

              if (list(1) == i1 .and. list(2) == i2 .or.  &
                  list(2) == i1 .and. list(1) == i2) then

                constraints%water_rOH = gromol%bonds(k)%b0 * 10.0_wip
                goto 1

              end if
            end do

            ioffset = ioffset + gromol%num_atoms

          end do

        else

          constraints%water_rOH = gromol%settles%doh * 10.0_wip
          goto 1

        end if

      end do

1     ioffset = 0

      do i = 1, grotop%num_molss
        gromol => grotop%molss(i)%moltype%mol

        if (gromol%settles%func == 0) then

          do j = 1, grotop%molss(i)%count

            do k = 1, gromol%num_bonds

              i1 = gromol%bonds(k)%atom_idx1 + ioffset
              i2 = gromol%bonds(k)%atom_idx2 + ioffset

              if (list(2) == i1 .and. list(3) == i2 .or.  &
                  list(3) == i1 .and. list(2) == i2) then

                constraints%water_rHH = gromol%bonds(k)%b0 * 10.0_wip
                goto 2

              end if
            end do

            ioffset = ioffset + gromol%num_atoms

          end do

        else

          constraints%water_rHH = gromol%settles%dhh * 10.0_wip
          goto 2

        end if

      end do
2     continue

    end if

    end if

    ! write parameters to MsgOut
    !
    if (main_rank) then
      if (constraints%water_type >= TIP3) then
        write(MsgOut,'(A)') 'Setup_Fast_Water> Setup constraints for SETTLE'
        write(MsgOut,'(A20,F10.4,A20,F10.4)')                &
             '  r0(O-H)         = ', constraints%water_rOH,  &
             '  mass(O)         = ', constraints%water_massO
        write(MsgOut,'(A20,F10.4,A20,F10.4)')                &
             '  r0(H-H)         = ', constraints%water_rHH,  &
             '  mass(H)         = ', constraints%water_massH
        write(MsgOut,'(A)') ' '
      else if (constraints%water_type == TIP2) then
        write(MsgOut,'(A)') 'Setup_Fast_Water> Setup constraints for WATER'
        write(MsgOut,'(A20,F10.4,A20,F10.4)')                &
             '  r0(1-2)         = ', constraints%water_r12    
        write(MsgOut,'(A20,F10.4,A20,F10.4)')                &
             '  mass(1)         = ', constraints%water_mass1,&
             '  mass(2)         = ', constraints%water_mass2
      end if
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine setup_fast_water

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_fast_water_tip4
  !> @brief        setup parameter for bond in water (TIP4 model)
  !! @authors      JJ
  !! @param[in]    par         : PAR information
  !! @param[in]    prmtop      : AMBER parameter topology information
  !! @param[in]    grotop      : GROMACS parameter topology information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_fast_water_tip4(par, prmtop, grotop, &
                                   molecule, enefunc, constraints)

    ! formal arguments
    type(s_par),             intent(in)    :: par
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_constraints),     intent(inout) :: constraints

    ! local variables
    integer                  :: i, i1, i2, list(4)
    integer                  :: ioffset, j, k
    real(wp)                 :: a, b
    character(6)             :: ci1, ci2, ci3, ci4

    type(s_grotop_mol), pointer :: gromol


    ! mass
    !
    constraints%water_massO = enefunc%table%mass_O
    constraints%water_massH = enefunc%table%mass_H

    ! min distance
    !
    list(1:4) = enefunc%table%water_list(1:4,1)

    ! charmm
    if (par%num_bonds > 0) then

      ci1 = molecule%atom_cls_name(list(1))
      ci2 = molecule%atom_cls_name(list(2))
      ci3 = molecule%atom_cls_name(list(3))
      ci4 = molecule%atom_cls_name(list(4))

      do i = 1, par%num_bonds
        if (ci1 .eq. par%bond_atom_cls(1,i) .and. &
            ci2 .eq. par%bond_atom_cls(2,i) .or.  &
            ci1 .eq. par%bond_atom_cls(2,i) .and. &
            ci2 .eq. par%bond_atom_cls(1,i) ) then
          constraints%water_rOH = par%bond_dist_min(i)
          exit
        end if
      end do

      do i = 1, par%num_bonds
        if (ci2 .eq. par%bond_atom_cls(1,i) .and. &
            ci3 .eq. par%bond_atom_cls(2,i)) then
          constraints%water_rHH = par%bond_dist_min(i)
          exit
        end if
      end do

      do i = 1, par%num_bonds
        if (ci1 .eq. par%bond_atom_cls(1,i) .and. &
            ci4 .eq. par%bond_atom_cls(2,i) .or.  &
            ci1 .eq. par%bond_atom_cls(2,i) .and. &
            ci4 .eq. par%bond_atom_cls(1,i) ) then
          constraints%water_rOD = par%bond_dist_min(i)
          exit
        end if
      end do

    ! amber
    else if (prmtop%num_atoms > 0) then

      do i = 1, prmtop%num_bondh

        i1 = prmtop%bond_inc_hy(1,i) / 3 + 1
        i2 = prmtop%bond_inc_hy(2,i) / 3 + 1

        if (list(1) == i1 .and. list(2) == i2 .or.  &
            list(2) == i1 .and. list(1) == i2) then

          constraints%water_rOH = &
               prmtop%bond_equil_uniq(prmtop%bond_inc_hy(3,i))
          exit

        end if

      end do

      do i = 1, prmtop%num_bondh

        i1 = prmtop%bond_inc_hy(1,i) / 3 + 1
        i2 = prmtop%bond_inc_hy(2,i) / 3 + 1

        if (list(2) == i1 .and. list(3) == i2 .or. &
            list(3) == i1 .and. list(2) == i2) then

          constraints%water_rHH = &
               prmtop%bond_equil_uniq(prmtop%bond_inc_hy(3,i))
          exit

        end if

      end do

      do i = 1, prmtop%num_mbonda

        i1 = prmtop%bond_wo_hy(1,i) / 3 + 1
        i2 = prmtop%bond_wo_hy(2,i) / 3 + 1

        if (list(1) == i1 .and. list(4) == i2 .or. &
            list(4) == i1 .and. list(1) == i2) then

          constraints%water_rOD = &
               prmtop%bond_equil_uniq(prmtop%bond_wo_hy(3,i))
          exit

        end if

      end do

    ! gromacs
    else if (grotop%num_atomtypes > 0) then

      ioffset = 0

      do i = 1, grotop%num_molss
        gromol => grotop%molss(i)%moltype%mol

        if (gromol%settles%func == 0) then

          do j = 1, grotop%molss(i)%count

            do k = 1, gromol%num_bonds

              i1 = gromol%bonds(k)%atom_idx1 + ioffset
              i2 = gromol%bonds(k)%atom_idx2 + ioffset

              if (list(1) == i1 .and. list(2) == i2 .or.  &
                  list(2) == i1 .and. list(1) == i2) then

                constraints%water_rOH = gromol%bonds(k)%b0 * 10.0_dp
                goto 1

              end if
            end do

            ioffset = ioffset + gromol%num_atoms

          end do

        else

          constraints%water_rOH = gromol%settles%doh * 10.0_dp
          goto 1

        end if

      end do

1     ioffset = 0

      do i = 1, grotop%num_molss
        gromol => grotop%molss(i)%moltype%mol

        if (gromol%settles%func == 0) then

          do j = 1, grotop%molss(i)%count

            do k = 1, gromol%num_bonds

              i1 = gromol%bonds(k)%atom_idx1 + ioffset
              i2 = gromol%bonds(k)%atom_idx2 + ioffset

              if (list(2) == i1 .and. list(3) == i2 .or.  &
                  list(3) == i1 .and. list(2) == i2) then

                constraints%water_rHH = gromol%bonds(k)%b0 * 10.0_dp
                goto 2

              end if
            end do

            ioffset = ioffset + gromol%num_atoms

          end do

        else

          constraints%water_rHH = gromol%settles%dhh * 10.0_dp
          goto 2

        end if

      end do

2     ioffset = 0

      do i = 1, grotop%num_molss
        gromol => grotop%molss(i)%moltype%mol

        do j = 1, grotop%molss(i)%count

          do k = 1, gromol%num_vsites3
            i1 = gromol%vsites3(k)%func

            if (i1 == 1) then

              ioffset = 1
              a = gromol%vsites3(k)%a
              b = gromol%vsites3(k)%b
              if (a /= b) call error_msg('TIP4P molecule is not symmetric')
              constraints%water_rOD = 2.0_wp * a  &
                                * sqrt(constraints%water_rOH**2 &
                                       -0.25_wp*constraints%water_rHH**2)
            end if
            if (ioffset == 1) exit
          end do

          if (ioffset == 1) exit

        end do

        if (ioffset == 1) exit
      end do

      if (ioffset /= 1) call error_msg('Virtual site should be defined &
                                     & when using TIP4P')

    end if


    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Setup_Fast_Water> Setup constraints for SETTLE'
      write(MsgOut,'(A20,F10.4,A20,F10.4)')                &
           '  r0(O-H)         = ', constraints%water_rOH,  &
           '  mass(O)         = ', constraints%water_massO
      write(MsgOut,'(A20,F10.4)')                          &
           '  r0(O-D)         = ', constraints%water_rOD
      write(MsgOut,'(A20,F10.4,A20,F10.4)')                &
           '  r0(H-H)         = ', constraints%water_rHH,  &
           '  mass(H)         = ', constraints%water_massH
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine setup_fast_water_tip4

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_tip4_dummy
  !> @brief        r(O-D) in TIP4P
  !! @authors      JJ
  !! @param[in]    par         : PAR information
  !! @param[in]    prmtop      : AMBER parameter topology information
  !! @param[in]    grotop      : GROMACS parameter topology information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_tip4_dummy(par, prmtop, grotop, &
                              molecule, enefunc, constraints)

    ! formal arguments
    type(s_par),             intent(in)    :: par
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_constraints),     intent(inout) :: constraints

    ! local variables
    integer                  :: i, i1, i2, list(4), ioffset, j, k
    real(wp)                 :: a, b
    character(6)             :: ci1, ci2, ci3, ci4

    type(s_grotop_mol), pointer :: gromol

    ! min distance
    !
    list(1:4) = enefunc%table%water_list(1:4,1)

    ! charmm
    if (par%num_bonds > 0) then

      ci1 = molecule%atom_cls_name(list(1))
      ci2 = molecule%atom_cls_name(list(2))
      ci3 = molecule%atom_cls_name(list(3))
      ci4 = molecule%atom_cls_name(list(4))

      do i = 1, par%num_bonds
        if (ci1 .eq. par%bond_atom_cls(1,i) .and. &
            ci4 .eq. par%bond_atom_cls(2,i) .or.  &
            ci1 .eq. par%bond_atom_cls(2,i) .and. &
            ci4 .eq. par%bond_atom_cls(1,i) ) then
          constraints%water_rOD = par%bond_dist_min(i)
          exit
        end if
      end do

    ! amber
    else if (prmtop%num_atoms > 0) then

      do i = 1, prmtop%num_mbonda

        i1 = prmtop%bond_wo_hy(1,i) / 3 + 1
        i2 = prmtop%bond_wo_hy(2,i) / 3 + 1

        if (list(1) == i1 .and. list(4) == i2 .or. &
            list(4) == i1 .and. list(1) == i2) then

          constraints%water_rOD = &
               prmtop%bond_equil_uniq(prmtop%bond_wo_hy(3,i))
          exit

        end if

      end do

    else if (grotop%num_atomtypes > 0) then

      ioffset = 0

      do i = 1, grotop%num_molss
        gromol => grotop%molss(i)%moltype%mol

        if (gromol%settles%func == 0) then

          do j = 1, grotop%molss(i)%count

            do k = 1, gromol%num_bonds

              i1 = gromol%bonds(k)%atom_idx1 + ioffset
              i2 = gromol%bonds(k)%atom_idx2 + ioffset

              if (list(1) == i1 .and. list(2) == i2 .or.  &
                  list(2) == i1 .and. list(1) == i2) then

                constraints%water_rOH = gromol%bonds(k)%b0 * 10.0_dp
                goto 1

              end if
            end do

            ioffset = ioffset + gromol%num_atoms

          end do

        else

          constraints%water_rOH = gromol%settles%doh * 10.0_dp
          goto 1

        end if

      end do

1     ioffset = 0

      do i = 1, grotop%num_molss
        gromol => grotop%molss(i)%moltype%mol

        if (gromol%settles%func == 0) then

          do j = 1, grotop%molss(i)%count

            do k = 1, gromol%num_bonds

              i1 = gromol%bonds(k)%atom_idx1 + ioffset
              i2 = gromol%bonds(k)%atom_idx2 + ioffset

              if (list(2) == i1 .and. list(3) == i2 .or.  &
                  list(3) == i1 .and. list(2) == i2) then

                constraints%water_rHH = gromol%bonds(k)%b0 * 10.0_dp
                goto 2

              end if
            end do

            ioffset = ioffset + gromol%num_atoms

          end do

        else

          constraints%water_rHH = gromol%settles%dhh * 10.0_dp
          goto 2

        end if

      end do

2     ioffset = 0

      do i = 1, grotop%num_molss

        gromol => grotop%molss(i)%moltype%mol

        do j = 1, grotop%molss(i)%count

          do k = 1, gromol%num_vsites3
            i1 = gromol%vsites3(k)%func

            if (i1 == 1) then

              ioffset = 1
              a = gromol%vsites3(k)%a
              b = gromol%vsites3(k)%b
              if (a /= b) call error_msg('TIP4P molecule is not symmetric')
              constraints%water_rOD = 2.0_wp * a  &
                                * sqrt(constraints%water_rOH**2 &
                                       -0.25_wp*constraints%water_rHH**2)
            end if
            if (ioffset == 1) exit
          end do

          if (ioffset == 1) exit

        end do

      end do

      if (ioffset /= 1) call error_msg('Virtual site should be defined &
                                     & when using TIP4P')

    end if


    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A20,F10.4)')                          &
           '  r0(O-D)         = ', constraints%water_rOD
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine setup_tip4_dummy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_fast_water_pio
  !> @brief        setup parameter for bond in water
  !! @authors      JJ
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_fast_water_pio(constraints)

    ! formal arguments
    type(s_constraints),     intent(inout) :: constraints


    if (main_rank) then
      write(MsgOut,'(A)') &
           'Setup_Fast_Water_Pio> Constraints for SETTLE'
      write(MsgOut,'(A20,F10.4,A20,F10.4)')                &
           '  r0(O-H)         = ', constraints%water_rOH,  &
           '  mass(O)         = ', constraints%water_massO
      write(MsgOut,'(A20,F10.4,A20,F10.4)')                &
           '  r0(H-H)         = ', constraints%water_rHH,  &
           '  mass(H)         = ', constraints%water_massH
      if (constraints%water_type == TIP4)                  &
      write(MsgOut,'(A20,F10.4,A20)')                      &
           '  r0(O-D)         = ', constraints%water_rOD
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine setup_fast_water_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_rigid_bond
  !> @brief        setup rigid bonds for SHAKE and RATTLE in domain
  !!               decomposition
  !! @authors      JJ
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_rigid_bond(constraints)

    ! formal arguments
    type(s_constraints),     intent(inout) :: constraints


    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') &
           'Setup_Rigid_Bond> Setup constrains for SHAKE and RATTLE'
      write(MsgOut,'(A20,I10)') &
           '  num_rigid_bonds = ', constraints%num_bonds
      write(MsgOut,'(A)') &
           ' '
    end if

    return

  end subroutine setup_rigid_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_rigid_bond_pio
  !> @brief        setup rigid bonds for SHAKE and RATTLE in domain
  !!               decomposition
  !! @authors      JJ
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_rigid_bond_pio(constraints)

    ! formal arguments
    type(s_constraints),     intent(inout) :: constraints


    if (main_rank) then
      write(MsgOut,'(A)') &
           'Setup_Rigid_Bond_Pio> Constrains for SHAKE and RATTLE'
      write(MsgOut,'(A20,I10)') &
           '  num_rigid_bonds = ', constraints%num_bonds
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine setup_rigid_bond_pio

 !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_mass_repartitioning
  !> @brief        Mass repartitioning of hydrogen group
  !! @authors      JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    constraints : constraints information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : energy function information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_mass_repartitioning(dynamics, constraints, domain, enefunc)

    ! formal arguments
    type(s_dynamics),    target, intent(in   ) :: dynamics   
    type(s_constraints), target, intent(in   ) :: constraints
    type(s_domain),      target, intent(inout) :: domain
    type(s_enefunc),     target, intent(inout) :: enefunc

    ! local variables
    integer                      :: icel, i, j, k, ih, connect
    integer                      :: atm1, atm2
    real(wip)                    :: total_mass, total_massH
    real(wp)                     :: factor

    real(wip),           pointer :: mass(:,:), inv_mass(:,:)
    real(wp),            pointer :: ratio, ratio_xh1, ratio_xh2, ratio_xh3
    real(wp),            pointer :: ratio_ring
    integer,             pointer :: hmr_target
    integer,             pointer :: natom(:)
    integer,             pointer :: ring(:,:)
    integer,             pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)
    integer,             pointer :: ncell, nwater(:), water_list(:,:,:)


    hmr_target      => dynamics%hmr_target
    ratio           => dynamics%hmr_ratio
    ratio_xh1       => dynamics%hmr_ratio_xh1
    ratio_xh2       => dynamics%hmr_ratio_xh2
    ratio_xh3       => dynamics%hmr_ratio_xh3
    ratio_ring      => dynamics%hmr_ratio_ring

    mass            => domain%mass
    inv_mass        => domain%inv_mass
    ncell           => domain%num_cell_local
    natom           => domain%num_atom
    nwater          => domain%num_water
    water_list      => domain%water_list
    ring            => domain%ring

    HGr_local       => constraints%HGr_local
    HGr_bond_list   => constraints%HGr_bond_list
    connect         =  constraints%connect

    ! repartitioning water molecule mass
    !
    if (hmr_target == HmrTargetAll .or. hmr_target == HmrTargetSolvent) then

      if (ratio_xh2 > EPS) then
        factor = ratio_xh2
      else
        factor = ratio
      end if
      total_mass = enefunc%table%mass_O + 2.0_wip*enefunc%table%mass_H
      enefunc%table%mass_H = enefunc%table%mass_H * factor
      enefunc%table%mass_O = total_mass - 2.0_wip*enefunc%table%mass_H
      domain%water%mass(1) = enefunc%table%mass_O
      domain%water%mass(2) = enefunc%table%mass_H
      domain%water%mass(3) = enefunc%table%mass_H

      do icel = 1, ncell
        do k = 1, nwater(icel)
          atm1 = water_list(1,k,icel)
          mass(atm1,icel) = enefunc%table%mass_O
          atm2 = water_list(2,k,icel)
          mass(atm2,icel) = enefunc%table%mass_H
          atm2 = water_list(3,k,icel)
          mass(atm2,icel) = enefunc%table%mass_H
        end do
      end do

    end if

    ! mass repartitioning of XHn group
    !
    if (hmr_target == HmrTargetAll .or. hmr_target == HmrTargetSolute) then

      do icel = 1, ncell

        do j = 1, connect
          do k = 1, HGr_local(j,icel)
            total_mass  = 0.0_wip
            total_massH = 0.0_wip
            atm1 = HGr_bond_list(1,k,j,icel)
            total_mass = total_mass + mass(atm1,icel)
            factor = ratio
            if (j == 1 .and. ratio_xh1 > EPS) then
              factor = ratio_xh1
            end if
            if (j == 2 .and. ratio_xh2 > EPS) then
              factor = ratio_xh2
            end if
            if (j == 3 .and. ratio_xh3 > EPS) then
              factor = ratio_xh3
            end if
            if (ratio_ring > EPS .and. ring(atm1,icel) == 1) then
              factor = ratio_ring
            end if

            do ih = 1, j
              atm2 = HGr_bond_list(ih+1,k,j,icel)
              total_mass  = total_mass  + mass(atm2,icel)
              mass(atm2,icel) = mass(atm2,icel)*factor
              total_massH = total_massH + mass(atm2,icel)
            end do
            mass(atm1,icel) = total_mass - total_massH
          end do
        end do
 
      end do
 
    end if
      
    ! inverse mass
    !
    do icel = 1, ncell
      do i = 1, natom(icel)
        if (mass(i,icel) > EPS) then
          inv_mass(i,icel) = 1.0_wip / mass(i,icel)
        else
          inv_mass(i,icel) = 0.0_wip
        end if
      end do

    end do

    return

  end subroutine setup_mass_repartitioning

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_mass_repartitioning_back
  !> @brief        Mass repartitioning of hydrogen group is back to original
  !! @authors      JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    constraints : constraints information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : energy function information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_mass_repartitioning_back(dynamics, constraints, domain, &
                                            enefunc)

    ! formal arguments
    type(s_dynamics),    target, intent(in   ) :: dynamics   
    type(s_constraints), target, intent(in   ) :: constraints
    type(s_domain),      target, intent(inout) :: domain
    type(s_enefunc),     target, intent(inout) :: enefunc

    ! local variables
    integer                      :: icel, i, j, k, ih, connect
    integer                      :: atm1, atm2
    real(wip)                    :: total_mass, total_massH
    real(wp)                     :: factor

    real(wip),           pointer :: mass(:,:), inv_mass(:,:)
    real(wp),            pointer :: ratio, ratio_xh1, ratio_xh2, ratio_xh3
    real(wp),            pointer :: ratio_ring
    integer,             pointer :: hmr_target
    integer,             pointer :: natom(:)
    integer,             pointer :: ring(:,:)
    integer,             pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)
    integer,             pointer :: ncell, nwater(:), water_list(:,:,:)


    hmr_target      => dynamics%hmr_target
    ratio           => dynamics%hmr_ratio
    ratio_xh1       => dynamics%hmr_ratio_xh1
    ratio_xh2       => dynamics%hmr_ratio_xh2
    ratio_xh3       => dynamics%hmr_ratio_xh3
    ratio_ring      => dynamics%hmr_ratio_ring

    mass            => domain%mass
    inv_mass        => domain%inv_mass
    ncell           => domain%num_cell_local
    natom           => domain%num_atom
    nwater          => domain%num_water
    water_list      => domain%water_list
    ring            => domain%ring

    HGr_local       => constraints%HGr_local
    HGr_bond_list   => constraints%HGr_bond_list
    connect         =  constraints%connect

    ! repartitioning water molecule mass
    !
    if (hmr_target == HmrTargetAll .or. hmr_target == HmrTargetSolvent) then

      if (ratio_xh2 > EPS) then
        factor = ratio_xh2
      else
        factor = ratio
      end if
      total_mass = enefunc%table%mass_O + 2.0_wip*enefunc%table%mass_H
      enefunc%table%mass_H = enefunc%table%mass_H / factor
      enefunc%table%mass_O = total_mass - 2.0_wip*enefunc%table%mass_H
      domain%water%mass(1) = enefunc%table%mass_O
      domain%water%mass(2) = enefunc%table%mass_H
      domain%water%mass(3) = enefunc%table%mass_H

      do icel = 1, ncell
        do k = 1, nwater(icel)
          atm1 = water_list(1,k,icel)
          mass(atm1,icel) = enefunc%table%mass_O
          atm2 = water_list(2,k,icel)
          mass(atm2,icel) = enefunc%table%mass_H
          atm2 = water_list(3,k,icel)
          mass(atm2,icel) = enefunc%table%mass_H
        end do
      end do

    end if

    ! mass repartitioning of XHn group
    !
    if (hmr_target == HmrTargetAll .or. hmr_target == HmrTargetSolute) then

      do icel = 1, ncell

        do j = 1, connect
          do k = 1, HGr_local(j,icel)
            total_mass  = 0.0_wip
            total_massH = 0.0_wip
            atm1 = HGr_bond_list(1,k,j,icel)
            total_mass = total_mass + mass(atm1,icel)
            factor = ratio
            if (j == 1 .and. ratio_xh1 > EPS) then
                factor = ratio_xh1
            end if
            if (j == 2 .and. ratio_xh2 > EPS) then
                factor = ratio_xh2
            end if
            if (j == 3 .and. ratio_xh3 > EPS) then
                factor = ratio_xh3
            end if
            if (ratio_ring > EPS .and. ring(atm1,icel) == 1) then
              factor = ratio_ring
            end if

            do ih = 1, j
              atm2 = HGr_bond_list(ih+1,k,j,icel)
              total_mass  = total_mass  + mass(atm2,icel)
              mass(atm2,icel) = mass(atm2,icel)/factor
              total_massH = total_massH + mass(atm2,icel)
            end do
            mass(atm1,icel) = total_mass - total_massH
          end do
        end do
 
      end do
 
    end if
      
    ! inverse mass
    !
    do icel = 1, ncell
      do i = 1, natom(icel)
        if (mass(i,icel) > EPS) then
          inv_mass(i,icel) = 1.0_wip / mass(i,icel)
        else
          inv_mass(i,icel) = 0.0_wip
        end if
      end do

    end do

    return

  end subroutine setup_mass_repartitioning_back

 !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_vel_group
  !> @brief        Update group velocity
  !! @authors      JJ
  !! @param[in]    constraints : constraints information
  !! @param[in]    ncell       : number of cells
  !! @param[in]    nwater      : number of water
  !! @param[in]    water_list  : water molecule list
  !! @param[in]    scale_vel   : scale velocity
  !! @param[in]    mass        : mass
  !! @param[inout] vel         : velocities
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_vel_group(constraints, ncell, nwater, water_list, &
                              scale_vel, mass, vel)

    ! formal arguments
    type(s_constraints), target, intent(in)    :: constraints
    integer,                     intent(in)    :: ncell, nwater(:)
    integer,                     intent(in)    :: water_list(:,:,:)
    real(wip),                   intent(in)    :: scale_vel
    real(wip),                   intent(in)    :: mass(:,:)
    real(wip),                   intent(inout) :: vel(:,:,:)

    ! local variables
    integer                      :: icel, i, ix, j, k, ih, connect
    integer                      :: id, omp_get_thread_num
    integer                      :: iatm(1:8), water_atom
    real(wip)                    :: total_mass
    real(wip)                    :: vel_cm(1:3)

    integer,             pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)
    integer,             pointer :: nsolute(:)


    HGr_local       => constraints%HGr_local
    HGr_bond_list   => constraints%HGr_bond_list
    nsolute         => constraints%No_HGr
    connect         =  constraints%connect

    water_atom = constraints%water_type

    !$omp parallel private(id, icel, j, k, ix, total_mass, vel_cm, &
    !$omp                  iatm, ih)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do icel = id+1, ncell, nthread

      do ix = 1, nsolute(icel)
        vel(1:3,ix,icel) = vel(1:3,ix,icel) * scale_vel
      end do

!ocl nosimd
      do j = 1, connect
        do k = 1, HGr_local(j,icel)
          total_mass  = 0.0_wip
          vel_cm(1:3) = 0.0_wip
          do ih = 1, j+1
            iatm(ih)    = HGr_bond_list(ih,k,j,icel)
            total_mass  = total_mass + mass(iatm(ih),icel)
            vel_cm(1:3) = vel_cm(1:3)  &
                        + mass(iatm(ih),icel)*vel(1:3,iatm(ih),icel)
          end do
          vel_cm(1:3) = vel_cm(1:3) / total_mass
          do ih = 1, j+1
            iatm(ih)    = HGr_bond_list(ih,k,j,icel)
            vel(1:3,iatm(ih),icel) = vel(1:3,iatm(ih),icel) &
                                   + (scale_vel-1.0_wip)*vel_cm(1:3)
          end do
        end do
      end do

!ocl nosimd
      do ix = 1, nwater(icel)
        total_mass  = 0.0_wip
        vel_cm(1:3) = 0.0_wip
        do ih = 1, water_atom
          iatm(ih) = water_list(ih,ix,icel)
          total_mass = total_mass + mass(iatm(ih),icel)
          vel_cm(1:3) = vel_cm(1:3) &
                      + mass(iatm(ih),icel)*vel(1:3,iatm(ih),icel)
        end do
        vel_cm(1:3) = vel_cm(1:3) / total_mass
        do ih = 1, water_atom
          iatm(ih) = water_list(ih,ix,icel)
          vel(1:3,iatm(ih),icel) = vel(1:3,iatm(ih),icel) &
                                 + (scale_vel-1.0_wip)*vel_cm(1:3)
        end do
      end do

    end do
    !$omp end parallel

    return

  end subroutine update_vel_group

 !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_vel_group_3d
  !> @brief        Update group velocity
  !! @authors      JJ
  !! @param[in]    constraints : constraints information
  !! @param[in]    ncell       : number of cells
  !! @param[in]    nwater      : number of water
  !! @param[in]    water_list  : water molecule list
  !! @param[in]    scale_vel   : scale velocity
  !! @param[in]    mass        : mass
  !! @param[inout] vel         : velocities
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_vel_group_3d(constraints, ncell, nwater, water_list, &
                                 scale_vel, mass, vel)

    ! formal arguments
    type(s_constraints), target, intent(in)    :: constraints
    integer,                     intent(in)    :: ncell, nwater(:)
    integer,                     intent(in)    :: water_list(:,:,:)
    real(wip),                   intent(in)    :: scale_vel(:)
    real(wip),                   intent(in)    :: mass(:,:)
    real(wip),                   intent(inout) :: vel(:,:,:)

    ! local variables
    integer                      :: icel, i, ix, j, k, ih, connect
    integer                      :: id, omp_get_thread_num
    integer                      :: iatm(1:8), water_atom
    real(wip)                    :: total_mass
    real(wip)                    :: vel_cm(1:3)

    integer,             pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)
    integer,             pointer :: nsolute(:)


    HGr_local       => constraints%HGr_local
    HGr_bond_list   => constraints%HGr_bond_list
    nsolute         => constraints%No_HGr
    connect         =  constraints%connect

    water_atom = constraints%water_type

    !$omp parallel private(id, icel, j, k, ix, total_mass, vel_cm, &
    !$omp                  iatm, ih)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do icel = id+1, ncell, nthread

      do ix = 1, nsolute(icel)
        vel(1:3,ix,icel) = vel(1:3,ix,icel) * scale_vel(1:3)
      end do

!ocl nosimd
      do j = 1, connect
        do k = 1, HGr_local(j,icel)
          total_mass  = 0.0_wip
          vel_cm(1:3) = 0.0_wip
          do ih = 1, j+1
            iatm(ih)    = HGr_bond_list(ih,k,j,icel)
            total_mass  = total_mass + mass(iatm(ih),icel)
            vel_cm(1:3) = vel_cm(1:3)  &
                        + mass(iatm(ih),icel)*vel(1:3,iatm(ih),icel)
          end do
          vel_cm(1:3) = vel_cm(1:3) / total_mass
          do ih = 1, j+1
            iatm(ih)    = HGr_bond_list(ih,k,j,icel)
            vel(1:3,iatm(ih),icel) = vel(1:3,iatm(ih),icel) &
                                   + (scale_vel(1:3)-1.0_wip)*vel_cm(1:3)
          end do
        end do
      end do

!ocl nosimd
      do ix = 1, nwater(icel)
        total_mass  = 0.0_wip
        vel_cm(1:3) = 0.0_wip
        do ih = 1, water_atom
          iatm(ih) = water_list(ih,ix,icel)
          total_mass = total_mass + mass(iatm(ih),icel)
          vel_cm(1:3) = vel_cm(1:3) &
                      + mass(iatm(ih),icel)*vel(1:3,iatm(ih),icel)
        end do
        vel_cm(1:3) = vel_cm(1:3) / total_mass
        do ih = 1, water_atom
          iatm(ih) = water_list(ih,ix,icel)
          vel(1:3,iatm(ih),icel) = vel(1:3,iatm(ih),icel) &
                                 + (scale_vel(1:3)-1.0_wip)*vel_cm(1:3)
        end do
      end do

    end do
    !$omp end parallel

    return

  end subroutine update_vel_group_3d

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_kin_group
  !> @brief        Compute Group kinetic energy
  !! @authors      JJ
  !! @param[in]    constraints : constraints information
  !! @param[in]    ncell       : number of cells
  !! @param[in]    nwater      : number of water
  !! @param[in]    water_list  : water molecule list
  !! @param[in]    mass        : mass
  !! @param[in]    vel         : velocities
  !! @param[inout] kin         : component of kinetic energy
  !! @param[inout] ekin        : kinetic energy
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_kin_group(constraints, ncell, nwater, water_list, mass, &
                               vel, kin, ekin)

    ! formal arguments
    type(s_constraints), target, intent(in)    :: constraints
    integer,                     intent(in)    :: ncell, nwater(:)
    integer,                     intent(in)    :: water_list(:,:,:)
    real(wip),                   intent(in)    :: mass(:,:)
    real(wip),                   intent(in)    :: vel(:,:,:)
    real(dp),                    intent(inout) :: kin(:), ekin

    ! local variables
    integer                      :: icel, i, ix, j, k, ih, connect
    integer                      :: id, omp_get_thread_num
    integer                      :: iatm(1:8), water_atom
    real(wip)                    :: total_mass
    real(wip)                    :: vel_cm(1:3)
    real(dp)                     :: kinetic_omp(1:3,1:nthread)

    integer,             pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)
    integer,             pointer :: nsolute(:)


    HGr_local       => constraints%HGr_local
    HGr_bond_list   => constraints%HGr_bond_list
    nsolute         => constraints%No_HGr
    connect         =  constraints%connect

    kinetic_omp(1:3,1:nthread) = 0.0_dp
    kin(1:3) = 0.0_dp

    water_atom = constraints%water_type

    !$omp parallel private(id, icel, j, k, ix, total_mass, vel_cm, &
    !$omp                  iatm, ih) 
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do icel = id+1, ncell, nthread

      do ix = 1, nsolute(icel)
        kinetic_omp(1:3,id+1) = kinetic_omp(1:3,id+1) &
                              + mass(ix,icel)*vel(1:3,ix,icel)*vel(1:3,ix,icel)
      end do

!ocl nosimd
      do j = 1, connect
        do k = 1, HGr_local(j,icel)
          total_mass  = 0.0_wip
          vel_cm(1:3) = 0.0_wip
          do ih = 1, j+1
            iatm(ih)    = HGr_bond_list(ih,k,j,icel)
            total_mass  = total_mass + mass(iatm(ih),icel)
            vel_cm(1:3) = vel_cm(1:3)  &
                        + mass(iatm(ih),icel)*vel(1:3,iatm(ih),icel)
          end do
          vel_cm(1:3) = vel_cm(1:3) / total_mass
          kinetic_omp(1:3,id+1) = kinetic_omp(1:3,id+1) &
                                + total_mass*vel_cm(1:3)*vel_cm(1:3)
        end do
      end do

!ocl nosimd
      do ix = 1, nwater(icel)
        iatm(1:water_atom) = water_list(1:water_atom,ix,icel)
        total_mass  = 0.0_wip
        vel_cm(1:3) = 0.0_wip
        do ih = 1, water_atom
          total_mass  = total_mass + mass(iatm(ih),icel)
          vel_cm(1:3) = vel_cm(1:3)  &
                      + mass(iatm(ih),icel)*vel(1:3,iatm(ih),icel)
        end do
        vel_cm(1:3) = vel_cm(1:3) / total_mass
        kinetic_omp(1:3,id+1) = kinetic_omp(1:3,id+1) &
                              + total_mass*vel_cm(1:3)*vel_cm(1:3)
      end do
    end do
    !$omp end parallel

    do id = 1, nthread
      kin(1)      = kin(1)      + kinetic_omp(1,id)
      kin(2)      = kin(2)      + kinetic_omp(2,id)
      kin(3)      = kin(3)      + kinetic_omp(3,id)
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, kin, 3, mpi_real8, mpi_sum, &
                       mpi_comm_country, ierror)
#endif

    ekin = 0.5_dp*(kin(1)+kin(2)+kin(3))

    return

  end subroutine compute_kin_group

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_group_deg_freedom
  !> @brief        Compute degree of freedom for pressure
  !! @authors      JJ
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_group_deg_freedom(domain, constraints)

    ! formal arguments
    type(s_domain),      target, intent(inout) :: domain
    type(s_constraints), target, intent(inout) :: constraints

    ! local variables
    integer                      :: icel, i, ix, j, k, ih, connect
    integer                      :: id, omp_get_thread_num
    integer                      :: iatm(8)

    integer,             pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)
    integer,             pointer :: nsolute(:), nwater(:)
    integer,             pointer :: degree


    degree          => domain%num_group_freedom
    nwater          => domain%num_water
    HGr_local       => constraints%HGr_local
    HGr_bond_list   => constraints%HGr_bond_list
    nsolute         => constraints%No_HGr
    connect         =  constraints%connect

    degree = 0

    do icel = 1, domain%num_cell_local

      degree = degree + nsolute(icel)

      do j = 1, connect
        degree = degree + HGr_local(j,icel)
      end do

      degree = degree + nwater(icel)

    end do

#ifdef HAVE_MPI_GENESIS      
    call mpi_allreduce(mpi_in_place, degree, 1, mpi_integer, mpi_sum, &
                       mpi_comm_country, ierror)
#endif

    degree = degree * 3 - 3

    if (main_rank) then
      write(MsgOut,'(A)') &
           'Update_Num_Deg_Freedom> Number of degrees of freedom was updated'
      write(MsgOut,'(A20,I12,A34)') &
           '  num_deg_freedom = ', degree, ' (Group temperature/pressure case)'
      write(MsgOut,'(A)') 
    end if

    return

  end subroutine compute_group_deg_freedom

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_virial_group
  !> @brief        Compute Group virial
  !! @authors      JJ
  !! @param[in]    constraints : constraints information
  !! @param[in]    ncell       : number of cells
  !! @param[in]    nwater      : number of water
  !! @param[in]    water_list  : water molecule list
  !! @param[in]    mass        : mass
  !! @param[in]    coord       : coordinates
  !! @param[in]    force       : force
  !! @param[inout] virial      : virial
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_virial_group(constraints, ncell, nwater, water_list, &
                                  mass, coord, force, virial)

    ! formal arguments
    type(s_constraints), target, intent(in)    :: constraints
    integer,                     intent(in)    :: ncell, nwater(:)
    integer,                     intent(in)    :: water_list(:,:,:)
    real(wip),                   intent(in)    :: mass(:,:)
    real(wip),                   intent(in)    :: coord(:,:,:)
    real(wip),                   intent(in)    :: force(:,:,:)
    real(dp),                    intent(inout) :: virial(:,:)

    ! local variables
    integer                      :: icel, i, ix, j, k, ih, connect
    integer                      :: id, omp_get_thread_num
    integer                      :: iatm(1:8), water_atom
    real(wip)                    :: total_mass
    real(wip)                    :: r_cm(1:3), f_tot(1:3)
    real(dp)                     :: virial_omp(1:3,1:nthread)

    integer,             pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)
    integer,             pointer :: nsolute(:)


    HGr_local       => constraints%HGr_local
    HGr_bond_list   => constraints%HGr_bond_list
    nsolute         => constraints%No_HGr
    connect         =  constraints%connect

    virial_omp(1:3,1:nthread)  = 0.0_dp

    water_atom = constraints%water_type

    !$omp parallel private(id, icel, j, k, ix, total_mass, r_cm, &
    !$omp                  f_tot, iatm, ih)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do icel = id+1, ncell, nthread

      ! group virial from XHn group
      !
!ocl nosimd
      do j = 1, connect
        do k = 1, HGr_local(j,icel)
          total_mass  = 0.0_wip
          r_cm(1:3)   = 0.0_wip
          f_tot(1:3)  = 0.0_wip
          do ih = 1, j+1
            iatm(ih)    = HGr_bond_list(ih,k,j,icel)
            total_mass  = total_mass + mass(iatm(ih),icel)
            r_cm(1:3)   = r_cm(1:3)  &
                        + mass(iatm(ih),icel)*coord(1:3,iatm(ih),icel)
            f_tot(1:3)  = f_tot(1:3) + force(1:3,iatm(ih),icel)
          end do
          r_cm(1:3)   = r_cm(1:3) / total_mass
          do ih = 1, j+1
            virial_omp(1:3,id+1) = virial_omp(1:3,id+1) &
                               + (coord(1:3,iatm(ih),icel)-r_cm(1:3)) &
                                 *force(1:3,iatm(ih),icel)
          end do
        end do
      end do

      ! group virial from water
      !
!ocl nosimd
      do ix = 1, nwater(icel)
        iatm(1:water_atom) = water_list(1:water_atom,ix,icel)
        total_mass  = 0.0_wip
        r_cm(1:3)   = 0.0_wip
        f_tot(1:3)  = 0.0_wip
        do ih = 1, water_atom
          total_mass  = total_mass + mass(iatm(ih),icel)
          r_cm(1:3)   = r_cm(1:3)  &
                      + mass(iatm(ih),icel)*coord(1:3,iatm(ih),icel)
          f_tot(1:3)  = f_tot(1:3) + force(1:3,iatm(ih),icel)
        end do
        r_cm(1:3)  = r_cm(1:3) / total_mass
        do ih = 1, water_atom
          virial_omp(1:3,id+1) = virial_omp(1:3,id+1) &
                             + (coord(1:3,iatm(ih),icel)-r_cm(1:3)) &
                               *force(1:3,iatm(ih),icel)
        end do
      end do
    end do
    !$omp end parallel

    do id = 1, nthread
      virial(1,1) = virial(1,1) - virial_omp(1,id)
      virial(2,2) = virial(2,2) - virial_omp(2,id)
      virial(3,3) = virial(3,3) - virial_omp(3,id)
    end do

    return

  end subroutine compute_virial_group

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_vv1_group
  !> @brief        VV1 with group pressures
  !! @authors      JJ
  !! @param[in]    constraints : constraints information
  !! @param[in]    ncell       : number of cells
  !! @param[in]    natom       : number of atoms
  !! @param[in]    nwater      : number of water
  !! @param[in]    water_list  : water molecule list
  !! @param[in]    mass        : mass
  !! @param[in]    inv_mass    : inverse mass
  !! @param[in]    force       : force
  !! @param[in]    coord_ref   : reference coordinates
  !! @param[in]    size_scale  : size scale
  !! @param[in]    vel_scale   : velocity scale
  !! @param[in]    dt          : time step
  !! @param[in]    half_dt     : half time step
  !! @param[inout] coord       : coordinates
  !! @param[inout] vel         : velocities
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_vv1_group(constraints, ncell, natom, nwater, water_list, &
                               mass, inv_mass, force, coord_ref, size_scale,  &
                               vel_scale, dt, half_dt, coord, vel)

    ! formal arguments
    type(s_constraints), target, intent(in)    :: constraints
    integer,                     intent(in)    :: ncell, natom(:)
    integer,                     intent(in)    :: nwater(:), water_list(:,:,:)
    real(wip),                   intent(in)    :: mass(:,:), inv_mass(:,:)
    real(wip),                   intent(in)    :: force(:,:,:)
    real(wip),                   intent(in)    :: coord_ref(:,:,:)
    real(wip),                   intent(in)    :: size_scale(:)
    real(wip),                   intent(in)    :: vel_scale(:)
    real(wip),                   intent(in)    :: dt, half_dt
    real(wip),                   intent(inout) :: coord(:,:,:)
    real(wip),                   intent(inout) :: vel(:,:,:)

    ! local variables
    integer                      :: icel, i, ix, j, k, ih, connect
    integer                      :: id, omp_get_thread_num
    integer                      :: iatm(1:8), water_atom
    real(wip)                    :: factor
    real(wip)                    :: total_mass
    real(wip)                    :: r_cm(1:3), vel_cm(1:3)
    real(wip)                    :: r(1:3,1:8), v(1:3,1:8)

    integer,             pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)
    integer,             pointer :: nsolute(:)


    HGr_local       => constraints%HGr_local
    HGr_bond_list   => constraints%HGr_bond_list
    nsolute         => constraints%No_HGr
    connect         =  constraints%connect

    water_atom = constraints%water_type

    !$omp parallel private(id, icel, j, k, total_mass, r_cm, vel_cm, &
    !$omp                  iatm, ih, r, v, ix, factor)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do icel = id+1, ncell, nthread

      ! No group
      !
      do ix = 1, nsolute(icel)
        coord(1,ix,icel) = size_scale(1)*coord_ref(1,ix,icel)
        coord(2,ix,icel) = size_scale(2)*coord_ref(2,ix,icel)
        coord(3,ix,icel) = size_scale(3)*coord_ref(3,ix,icel)
        vel(1,ix,icel)   = vel_scale(1)*vel(1,ix,icel)
        vel(2,ix,icel)   = vel_scale(2)*vel(2,ix,icel)
        vel(3,ix,icel)   = vel_scale(3)*vel(3,ix,icel)
      end do

      ! XHn group
      !
!ocl nosimd
      do j = 1, connect
        do k = 1, HGr_local(j,icel)
          total_mass  = 0.0_wip
          r_cm(1:3)   = 0.0_wip
          vel_cm(1)   = 0.0_wip
          vel_cm(2)   = 0.0_wip
          vel_cm(3)   = 0.0_wip
          do ih = 1, j+1
            iatm(ih)    = HGr_bond_list(ih,k,j,icel)
            total_mass  = total_mass + mass(iatm(ih),icel)
            r_cm(1)     = r_cm(1)  &
                        + mass(iatm(ih),icel)*coord_ref(1,iatm(ih),icel)
            r_cm(2)     = r_cm(2)  &
                        + mass(iatm(ih),icel)*coord_ref(2,iatm(ih),icel)
            r_cm(3)     = r_cm(3)  &
                        + mass(iatm(ih),icel)*coord_ref(3,iatm(ih),icel)
            vel_cm(1)   = vel_cm(1)  &
                        + mass(iatm(ih),icel)*vel(1,iatm(ih),icel)
            vel_cm(2)   = vel_cm(2)  &
                        + mass(iatm(ih),icel)*vel(2,iatm(ih),icel)
            vel_cm(3)   = vel_cm(3)  &
                        + mass(iatm(ih),icel)*vel(3,iatm(ih),icel)
          end do
          r_cm(1)   = r_cm(1) / total_mass
          r_cm(2)   = r_cm(2) / total_mass
          r_cm(3)   = r_cm(3) / total_mass
          vel_cm(1) = vel_cm(1) / total_mass
          vel_cm(2) = vel_cm(2) / total_mass
          vel_cm(3) = vel_cm(3) / total_mass
          do ih = 1, j+1
            coord(1,iatm(ih),icel) = coord_ref(1,iatm(ih),icel) &
                                   + (size_scale(1)-1.0_wip)*r_cm(1)
            coord(2,iatm(ih),icel) = coord_ref(2,iatm(ih),icel) &
                                   + (size_scale(2)-1.0_wip)*r_cm(2)
            coord(3,iatm(ih),icel) = coord_ref(3,iatm(ih),icel) &
                                   + (size_scale(3)-1.0_wip)*r_cm(3)
            vel(1,iatm(ih),icel)   = vel(1,iatm(ih),icel)   &
                                   + (vel_scale(1)-1.0_wip)*vel_cm(1)
            vel(2,iatm(ih),icel)   = vel(2,iatm(ih),icel)   &
                                   + (vel_scale(2)-1.0_wip)*vel_cm(2)
            vel(3,iatm(ih),icel)   = vel(3,iatm(ih),icel)   &
                                   + (vel_scale(3)-1.0_wip)*vel_cm(3)
          end do

        end do
      end do

      ! Water
      !
!ocl nosimd
      do ix = 1, nwater(icel)
        iatm(1:water_atom) = water_list(1:water_atom,ix,icel)
        total_mass  = 0.0_wip
        r_cm(1)   = 0.0_wip
        r_cm(2)   = 0.0_wip
        r_cm(3)   = 0.0_wip
        vel_cm(1) = 0.0_wip
        vel_cm(2) = 0.0_wip
        vel_cm(3) = 0.0_wip
        do ih = 1, water_atom
          total_mass  = total_mass + mass(iatm(ih),icel)
          r_cm(1)   = r_cm(1)  &
                    + mass(iatm(ih),icel)*coord_ref(1,iatm(ih),icel)
          r_cm(2)   = r_cm(2)  &
                    + mass(iatm(ih),icel)*coord_ref(2,iatm(ih),icel)
          r_cm(3)   = r_cm(3)  &
                    + mass(iatm(ih),icel)*coord_ref(3,iatm(ih),icel)
          vel_cm(1) = vel_cm(1)  &
                    + mass(iatm(ih),icel)*vel(1,iatm(ih),icel)
          vel_cm(2) = vel_cm(2)  &
                    + mass(iatm(ih),icel)*vel(2,iatm(ih),icel)
          vel_cm(3) = vel_cm(3)  &
                    + mass(iatm(ih),icel)*vel(3,iatm(ih),icel)
        end do
        r_cm(1)  = r_cm(1) / total_mass
        r_cm(2)  = r_cm(2) / total_mass
        r_cm(3)  = r_cm(3) / total_mass
        vel_cm(1) = vel_cm(1) / total_mass
        vel_cm(2) = vel_cm(2) / total_mass
        vel_cm(3) = vel_cm(3) / total_mass
        do ih = 1, water_atom
          coord(1,iatm(ih),icel) = coord_ref(1,iatm(ih),icel) &
                                 + (size_scale(1)-1.0_wip)*r_cm(1)
          coord(2,iatm(ih),icel) = coord_ref(2,iatm(ih),icel) &
                                 + (size_scale(2)-1.0_wip)*r_cm(2)
          coord(3,iatm(ih),icel) = coord_ref(3,iatm(ih),icel) &
                                 + (size_scale(3)-1.0_wip)*r_cm(3)
          vel(1,iatm(ih),icel)   = vel(1,iatm(ih),icel)   &
                                 + (vel_scale(1)-1.0_wip)*vel_cm(1)
          vel(2,iatm(ih),icel)   = vel(2,iatm(ih),icel)   &
                                 + (vel_scale(2)-1.0_wip)*vel_cm(2)
          vel(3,iatm(ih),icel)   = vel(3,iatm(ih),icel)   &
                                 + (vel_scale(3)-1.0_wip)*vel_cm(3)

        end do
      end do
    end do

    do icel = id+1, ncell, nthread
      do ix = 1, natom(icel)
        factor = half_dt * inv_mass(ix,icel)
        vel(1,ix,icel) = vel(1,ix,icel) + factor*force(1,ix,icel)
        vel(2,ix,icel) = vel(2,ix,icel) + factor*force(2,ix,icel)
        vel(3,ix,icel) = vel(3,ix,icel) + factor*force(3,ix,icel)
        coord(1,ix,icel) = coord(1,ix,icel)+vel(1,ix,icel)*dt
        coord(2,ix,icel) = coord(2,ix,icel)+vel(2,ix,icel)*dt
        coord(3,ix,icel) = coord(3,ix,icel)+vel(3,ix,icel)*dt
      end do
    end do

    !$omp end parallel

    return

  end subroutine compute_vv1_group

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_vv1_coord_group
  !> @brief        VV1 with group T/P (coordinate only)
  !! @authors      JJ
  !! @param[in]    constraints : constraints information
  !! @param[in]    ncell       : number of cells
  !! @param[in]    natom       : number of atoms
  !! @param[in]    nwater      : number of water
  !! @param[in]    water_list  : water molecule list
  !! @param[in]    mass        : mass
  !! @param[in]    inv_mass    : inverse mass
  !! @param[in]    force       : force
  !! @param[in]    coord_ref   : reference coordinates
  !! @param[in]    vel         : velocities
  !! @param[in]    size_scale  : size scale
  !! @param[in]    dt          : time step
  !! @param[inout] coord       : coordinates
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_vv1_coord_group(constraints, ncell, natom, nwater,       &
                                     water_list, mass, inv_mass, force,       &
                                     coord_ref, vel, size_scale,  &
                                     dt, coord)

    ! formal arguments
    type(s_constraints), target, intent(in)    :: constraints
    integer,                     intent(in)    :: ncell, natom(:)
    integer,                     intent(in)    :: nwater(:), water_list(:,:,:)
    real(wip),                   intent(in)    :: mass(:,:), inv_mass(:,:)
    real(wip),                   intent(in)    :: force(:,:,:)
    real(wip),                   intent(in)    :: coord_ref(:,:,:)
    real(wip),                   intent(in)    :: vel(:,:,:)
    real(wip),                   intent(in)    :: size_scale(:)
    real(wip),                   intent(in)    :: dt
    real(wip),                   intent(inout) :: coord(:,:,:)

    ! local variables
    integer                      :: icel, i, ix, j, k, ih, connect
    integer                      :: id, omp_get_thread_num
    integer                      :: iatm(1:8), water_atom
    real(wip)                    :: factor
    real(wip)                    :: total_mass
    real(wip)                    :: r_cm(1:3), vel_cm(1:3)
    real(wip)                    :: r(1:3,1:8), v(1:3,1:8)

    integer,             pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)
    integer,             pointer :: nsolute(:)


    HGr_local       => constraints%HGr_local
    HGr_bond_list   => constraints%HGr_bond_list
    nsolute         => constraints%No_HGr
    connect         =  constraints%connect

    water_atom = constraints%water_type

    !$omp parallel private(id, icel, j, k, total_mass, r_cm, vel_cm, &
    !$omp                  iatm, ih, r, v, ix, factor)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    ! scale of coordinates
    !
    do icel = id+1, ncell, nthread

      ! No group
      !
!ocl nosimd
      do ix = 1, nsolute(icel)
        coord(1,ix,icel) = size_scale(1)*coord_ref(1,ix,icel)
        coord(2,ix,icel) = size_scale(2)*coord_ref(2,ix,icel)
        coord(3,ix,icel) = size_scale(3)*coord_ref(3,ix,icel)
      end do

      ! XHn group
      !
!ocl nosimd
      do j = 1, connect
        do k = 1, HGr_local(j,icel)
          total_mass  = 0.0_wip
          r_cm(1)     = 0.0_wip
          r_cm(2)     = 0.0_wip
          r_cm(3)     = 0.0_wip
          vel_cm(1)   = 0.0_wip
          vel_cm(2)   = 0.0_wip
          vel_cm(3)   = 0.0_wip
          iatm(1:j+1) = HGr_bond_list(1:j+1,k,j,icel)
          do ih = 1, j+1
            total_mass  = total_mass + mass(iatm(ih),icel)
            r_cm(1)     = r_cm(1)  &
                        + mass(iatm(ih),icel)*coord_ref(1,iatm(ih),icel)
            r_cm(2)     = r_cm(2)  &
                        + mass(iatm(ih),icel)*coord_ref(2,iatm(ih),icel)
            r_cm(3)     = r_cm(3)  &
                        + mass(iatm(ih),icel)*coord_ref(3,iatm(ih),icel)
            vel_cm(1)   = vel_cm(1)  &
                        + mass(iatm(ih),icel)*vel(1,iatm(ih),icel)
            vel_cm(2)   = vel_cm(2)  &
                        + mass(iatm(ih),icel)*vel(2,iatm(ih),icel)
            vel_cm(3)   = vel_cm(3)  &
                        + mass(iatm(ih),icel)*vel(3,iatm(ih),icel)
          end do
          r_cm(1)   = r_cm(1) / total_mass
          r_cm(2)   = r_cm(2) / total_mass
          r_cm(3)   = r_cm(3) / total_mass
          vel_cm(1) = vel_cm(1) / total_mass
          vel_cm(2) = vel_cm(2) / total_mass
          vel_cm(3) = vel_cm(3) / total_mass
          do ih = 1, j+1
            coord(1,iatm(ih),icel) = coord_ref(1,iatm(ih),icel) &
                                   + (size_scale(1)-1.0_wip)*r_cm(1)
            coord(2,iatm(ih),icel) = coord_ref(2,iatm(ih),icel) &
                                   + (size_scale(2)-1.0_wip)*r_cm(2)
            coord(3,iatm(ih),icel) = coord_ref(3,iatm(ih),icel) &
                                   + (size_scale(3)-1.0_wip)*r_cm(3)
          end do

        end do
      end do

      ! Water
      !
!ocl nosimd
      do ix = 1, nwater(icel)
        iatm(1:water_atom) = water_list(1:water_atom,ix,icel)
        total_mass  = 0.0_wip
        r_cm(1)   = 0.0_wip
        r_cm(2)   = 0.0_wip
        r_cm(3)   = 0.0_wip
        vel_cm(1) = 0.0_wip
        vel_cm(2) = 0.0_wip
        vel_cm(3) = 0.0_wip
        do ih = 1, water_atom
          total_mass  = total_mass + mass(iatm(ih),icel)
          r_cm(1)   = r_cm(1)  &
                    + mass(iatm(ih),icel)*coord_ref(1,iatm(ih),icel)
          r_cm(2)   = r_cm(2)  &
                    + mass(iatm(ih),icel)*coord_ref(2,iatm(ih),icel)
          r_cm(3)   = r_cm(3)  &
                    + mass(iatm(ih),icel)*coord_ref(3,iatm(ih),icel)
          vel_cm(1) = vel_cm(1)  &
                    + mass(iatm(ih),icel)*vel(1,iatm(ih),icel)
          vel_cm(2) = vel_cm(2)  &
                    + mass(iatm(ih),icel)*vel(2,iatm(ih),icel)
          vel_cm(3) = vel_cm(3)  &
                    + mass(iatm(ih),icel)*vel(3,iatm(ih),icel)
        end do
        r_cm(1)  = r_cm(1) / total_mass
        r_cm(2)  = r_cm(2) / total_mass
        r_cm(3)  = r_cm(3) / total_mass
        vel_cm(1) = vel_cm(1) / total_mass
        vel_cm(2) = vel_cm(2) / total_mass
        vel_cm(3) = vel_cm(3) / total_mass
        do ih = 1, water_atom 
          coord(1,iatm(ih),icel) = coord_ref(1,iatm(ih),icel) &
                                 + (size_scale(1)-1.0_wip)*r_cm(1)
          coord(2,iatm(ih),icel) = coord_ref(2,iatm(ih),icel) &
                                 + (size_scale(2)-1.0_wip)*r_cm(2)
          coord(3,iatm(ih),icel) = coord_ref(3,iatm(ih),icel) &
                                 + (size_scale(3)-1.0_wip)*r_cm(3)
        end do
      end do
    end do

    ! update coordiantes
    !
    do icel = id+1, ncell, nthread
      do ix = 1, natom(icel)
        coord(1,ix,icel) = coord(1,ix,icel)+vel(1,ix,icel)*dt
        coord(2,ix,icel) = coord(2,ix,icel)+vel(2,ix,icel)*dt
        coord(3,ix,icel) = coord(3,ix,icel)+vel(3,ix,icel)*dt
      end do
    end do

    ! scale of coordinates
    !
    do icel = id+1, ncell, nthread

      ! No group
      !
      do ix = 1, nsolute(icel)
        coord(1,ix,icel) = size_scale(1)*coord(1,ix,icel)
        coord(2,ix,icel) = size_scale(2)*coord(2,ix,icel)
        coord(3,ix,icel) = size_scale(3)*coord(3,ix,icel)
      end do

      ! XHn group
      !
!ocl nosimd
      do j = 1, connect
        do k = 1, HGr_local(j,icel)
          total_mass  = 0.0_wip
          r_cm(1)   = 0.0_wip
          r_cm(2)   = 0.0_wip
          r_cm(3)   = 0.0_wip
          iatm(1:j+1) = HGr_bond_list(1:j+1,k,j,icel)
          do ih = 1, j+1
            total_mass  = total_mass + mass(iatm(ih),icel)
            r_cm(1)   = r_cm(1)  &
                      + mass(iatm(ih),icel)*coord(1,iatm(ih),icel)
            r_cm(2)   = r_cm(2)  &
                      + mass(iatm(ih),icel)*coord(2,iatm(ih),icel)
            r_cm(3)   = r_cm(3)  &
                      + mass(iatm(ih),icel)*coord(3,iatm(ih),icel)
          end do
          r_cm(1)   = r_cm(1) / total_mass
          r_cm(2)   = r_cm(2) / total_mass
          r_cm(3)   = r_cm(3) / total_mass
          do ih = 1, j+1
            coord(1,iatm(ih),icel) = coord(1,iatm(ih),icel) &
                                   + (size_scale(1)-1.0_wip)*r_cm(1)
            coord(2,iatm(ih),icel) = coord(2,iatm(ih),icel) &
                                   + (size_scale(2)-1.0_wip)*r_cm(2)
            coord(3,iatm(ih),icel) = coord(3,iatm(ih),icel) &
                                   + (size_scale(3)-1.0_wip)*r_cm(3)
          end do

        end do
      end do

      ! Water
      !
!ocl nosimd
      do ix = 1, nwater(icel)
        iatm(1:water_atom) = water_list(1:water_atom,ix,icel)
        total_mass  = 0.0_wip
        r_cm(1)   = 0.0_wip
        r_cm(2)   = 0.0_wip
        r_cm(3)   = 0.0_wip
        do ih = 1, water_atom
          total_mass  = total_mass + mass(iatm(ih),icel)
          r_cm(1)   = r_cm(1)  &
                    + mass(iatm(ih),icel)*coord(1,iatm(ih),icel)
          r_cm(2)   = r_cm(2)  &
                    + mass(iatm(ih),icel)*coord(2,iatm(ih),icel)
          r_cm(3)   = r_cm(3)  &
                    + mass(iatm(ih),icel)*coord(3,iatm(ih),icel)
        end do
        r_cm(1)  = r_cm(1) / total_mass
        r_cm(2)  = r_cm(2) / total_mass
        r_cm(3)  = r_cm(3) / total_mass
        do ih = 1, water_atom
          coord(1,iatm(ih),icel) = coord(1,iatm(ih),icel) &
                                 + (size_scale(1)-1.0_wip)*r_cm(1)
          coord(2,iatm(ih),icel) = coord(2,iatm(ih),icel) &
                                 + (size_scale(2)-1.0_wip)*r_cm(2)
          coord(3,iatm(ih),icel) = coord(3,iatm(ih),icel) &
                                 + (size_scale(3)-1.0_wip)*r_cm(3)
        end do
      end do
    end do

    !$omp end parallel

    return

  end subroutine compute_vv1_coord_group

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_vv2_group
  !> @brief        VV2 with group pressures
  !! @authors      JJ
  !! @param[in]    constraints : constraints information
  !! @param[in]    ncell       : number of cells
  !! @param[in]    natom       : number of atoms
  !! @param[in]    nwater      : number of water
  !! @param[in]    water_list  : water molecule list
  !! @param[in]    mass        : mass
  !! @param[in]    force       : force
  !! @param[in]    vel_scale   : velocity scale
  !! @param[in]    half_dt     : half time step
  !! @param[inout] vel         : velocities
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_vv2_group(constraints, ncell, natom, nwater,    &
                               water_list, mass, imass, force,       &
                               vel_scale, half_dt, vel)

    ! formal arguments
    type(s_constraints), target, intent(in)    :: constraints
    integer,                     intent(in)    :: ncell, natom(:)
    integer,                     intent(in)    :: nwater(:), water_list(:,:,:)
    real(wip),                   intent(in)    :: mass(:,:)
    real(wip),                   intent(in)    :: imass(:,:)
    real(wip),                   intent(in)    :: force(:,:,:)
    real(wip),                   intent(in)    :: vel_scale(:)
    real(wip),                   intent(in)    :: half_dt
    real(wip),                   intent(inout) :: vel(:,:,:)

    ! local variables
    integer                      :: icel, i, ix, j, k, ih, connect
    integer                      :: id, omp_get_thread_num
    integer                      :: iatm(1:8), water_atom
    real(wip)                    :: factor
    real(wip)                    :: total_mass
    real(wip)                    :: vel_cm(1:3)
    real(wip)                    :: v(1:3,1:8)

    integer,             pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)
    integer,             pointer :: nsolute(:)


    HGr_local       => constraints%HGr_local
    HGr_bond_list   => constraints%HGr_bond_list
    nsolute         => constraints%No_HGr
    connect         =  constraints%connect

    water_atom = constraints%water_type

    !$omp parallel private(id, icel, j, k, total_mass, vel_cm, &
    !$omp                  iatm, ih, v, ix, factor)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do icel = id+1, ncell, nthread
      do ix = 1, natom(icel)
        factor = half_dt * imass(ix,icel)
        vel(1,ix,icel) = vel(1,ix,icel) + factor*force(1,ix,icel)
        vel(2,ix,icel) = vel(2,ix,icel) + factor*force(2,ix,icel)
        vel(3,ix,icel) = vel(3,ix,icel) + factor*force(3,ix,icel)
      end do
    end do

    do icel = id+1, ncell, nthread

      ! No group
      !
      do ix = 1, nsolute(icel)
        vel(1,ix,icel)   = vel(1,ix,icel)*vel_scale(1)
        vel(2,ix,icel)   = vel(2,ix,icel)*vel_scale(2)
        vel(3,ix,icel)   = vel(3,ix,icel)*vel_scale(3)
      end do

      ! XHn group
      !
!ocl nosimd
      do j = 1, connect
        do k = 1, HGr_local(j,icel)
          total_mass  = 0.0_wip
          vel_cm(1) = 0.0_wip
          vel_cm(2) = 0.0_wip
          vel_cm(3) = 0.0_wip
          do ih = 1, j+1
            iatm(ih)    = HGr_bond_list(ih,k,j,icel)
            total_mass  = total_mass + mass(iatm(ih),icel)
            vel_cm(1)   = vel_cm(1)  &
                        + mass(iatm(ih),icel)*vel(1,iatm(ih),icel)
            vel_cm(2)   = vel_cm(2)  &
                        + mass(iatm(ih),icel)*vel(2,iatm(ih),icel)
            vel_cm(3)   = vel_cm(3)  &
                        + mass(iatm(ih),icel)*vel(3,iatm(ih),icel)
          end do
          vel_cm(1) = vel_cm(1) / total_mass
          vel_cm(2) = vel_cm(2) / total_mass
          vel_cm(3) = vel_cm(3) / total_mass
          do ih = 1, j+1
            vel(1,iatm(ih),icel) = vel(1,iatm(ih),icel)   &
                                 + (vel_scale(1)-1.0_wip)*vel_cm(1)
            vel(2,iatm(ih),icel) = vel(2,iatm(ih),icel)   &
                                 + (vel_scale(2)-1.0_wip)*vel_cm(2)
            vel(3,iatm(ih),icel) = vel(3,iatm(ih),icel)   &
                                 + (vel_scale(3)-1.0_wip)*vel_cm(3)
          end do

        end do
      end do

      ! Water
      !
!ocl nosimd
      do ix = 1, nwater(icel)
        iatm(1:water_atom) = water_list(1:water_atom,ix,icel)
        total_mass  = 0.0_wip
        vel_cm(1) = 0.0_wip
        vel_cm(2) = 0.0_wip
        vel_cm(3) = 0.0_wip
        do ih = 1, water_atom
          total_mass  = total_mass + mass(iatm(ih),icel)
          vel_cm(1)   = vel_cm(1)  &
                      + mass(iatm(ih),icel)*vel(1,iatm(ih),icel)
          vel_cm(2)   = vel_cm(2)  &
                      + mass(iatm(ih),icel)*vel(2,iatm(ih),icel)
          vel_cm(3)   = vel_cm(3)  &
                      + mass(iatm(ih),icel)*vel(3,iatm(ih),icel)
        end do
        vel_cm(1) = vel_cm(1) / total_mass
        vel_cm(2) = vel_cm(2) / total_mass
        vel_cm(3) = vel_cm(3) / total_mass
        do ih = 1, water_atom
          vel(1,iatm(ih),icel) = vel(1,iatm(ih),icel)   &
                               + (vel_scale(1)-1.0_wip)*vel_cm(1)
          vel(2,iatm(ih),icel) = vel(2,iatm(ih),icel)   &
                               + (vel_scale(2)-1.0_wip)*vel_cm(2)
          vel(3,iatm(ih),icel) = vel(3,iatm(ih),icel)   &
                               + (vel_scale(3)-1.0_wip)*vel_cm(3)

        end do
      end do
    end do

    !$omp end parallel

    return

  end subroutine compute_vv2_group

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_vel_vv1_group
  !> @brief        VV2 with group pressures
  !! @authors      JJ
  !! @param[in]    constraints : constraints information
  !! @param[in]    ncell       : number of cells
  !! @param[in]    natom       : number of atoms
  !! @param[in]    nwater      : number of water
  !! @param[in]    water_list  : water molecule list
  !! @param[in]    mass        : mass
  !! @param[in]    force       : force
  !! @param[in]    vel_scale   : velocity scale
  !! @param[in]    half_dt     : half time step
  !! @param[inout] vel         : velocities
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_vel_vv1_group(constraints, ncell, natom, nwater,    &
                                  water_list, mass, imass, force,       &
                                  vel_scale, half_dt, vel)

    ! formal arguments
    type(s_constraints), target, intent(in)    :: constraints
    integer,                     intent(in)    :: ncell, natom(:)
    integer,                     intent(in)    :: nwater(:), water_list(:,:,:)
    real(wip),                   intent(in)    :: mass(:,:)
    real(wip),                   intent(in)    :: imass(:,:)
    real(wip),                   intent(in)    :: force(:,:,:)
    real(wip),                   intent(in)    :: vel_scale(:)
    real(wip),                   intent(in)    :: half_dt
    real(wip),                   intent(inout) :: vel(:,:,:)

    ! local variables
    integer                      :: icel, i, ix, j, k, ih, connect
    integer                      :: id, omp_get_thread_num
    integer                      :: iatm(1:8), water_atom
    real(wip)                    :: factor
    real(wip)                    :: total_mass
    real(wip)                    :: vel_cm(1:3)
    real(wip)                    :: v(1:3,1:8)

    integer,             pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)
    integer,             pointer :: nsolute(:)


    HGr_local       => constraints%HGr_local
    HGr_bond_list   => constraints%HGr_bond_list
    nsolute         => constraints%No_HGr
    connect         =  constraints%connect

    water_atom = constraints%water_type

    !$omp parallel private(id, icel, j, k, total_mass, vel_cm, &
    !$omp                  iatm, ih, v, ix, factor)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do icel = id+1, ncell, nthread

      ! No group
      !
      do ix = 1, nsolute(icel)
        vel(1,ix,icel)   = vel(1,ix,icel)*vel_scale(1)
        vel(2,ix,icel)   = vel(2,ix,icel)*vel_scale(2)
        vel(3,ix,icel)   = vel(3,ix,icel)*vel_scale(3)
      end do

      ! XHn group
      !
!ocl nosimd
      do j = 1, connect
        do k = 1, HGr_local(j,icel)
          total_mass  = 0.0_wip
          vel_cm(1)   = 0.0_wip
          do ih = 1, j+1
            iatm(ih)    = HGr_bond_list(ih,k,j,icel)
            total_mass  = total_mass + mass(iatm(ih),icel)
            vel_cm(1)   = vel_cm(1)  &
                        + mass(iatm(ih),icel)*vel(1,iatm(ih),icel)
            vel_cm(2)   = vel_cm(2)  &
                        + mass(iatm(ih),icel)*vel(2,iatm(ih),icel)
            vel_cm(3)   = vel_cm(3)  &
                        + mass(iatm(ih),icel)*vel(3,iatm(ih),icel)
          end do
          vel_cm(1) = vel_cm(1) / total_mass
          vel_cm(2) = vel_cm(2) / total_mass
          vel_cm(3) = vel_cm(3) / total_mass
          do ih = 1, j+1
            vel(1,iatm(ih),icel)= vel(1,iatm(ih),icel)   &
                                + (vel_scale(1)-1.0_wp)*vel_cm(1)
            vel(2,iatm(ih),icel)= vel(2,iatm(ih),icel)   &
                                + (vel_scale(2)-1.0_wp)*vel_cm(2)
            vel(3,iatm(ih),icel)= vel(3,iatm(ih),icel)   &
                                + (vel_scale(3)-1.0_wp)*vel_cm(3)
          end do

        end do
      end do

      ! Water
      !
!ocl nosimd
      do ix = 1, nwater(icel)
        iatm(1:water_atom) = water_list(1:water_atom,ix,icel)
        total_mass  = 0.0_wip
        vel_cm(1)   = 0.0_wip
        vel_cm(2)   = 0.0_wip
        vel_cm(3)   = 0.0_wip
        do ih = 1, water_atom
          total_mass  = total_mass + mass(iatm(ih),icel)
          vel_cm(1) = vel_cm(1)  &
                    + mass(iatm(ih),icel)*vel(1,iatm(ih),icel)
          vel_cm(2) = vel_cm(2)  &
                    + mass(iatm(ih),icel)*vel(2,iatm(ih),icel)
          vel_cm(3) = vel_cm(3)  &
                    + mass(iatm(ih),icel)*vel(3,iatm(ih),icel)
        end do
        vel_cm(1) = vel_cm(1) / total_mass
        vel_cm(2) = vel_cm(2) / total_mass
        vel_cm(3) = vel_cm(3) / total_mass
        do ih = 1, water_atom
          vel(1,iatm(ih),icel) = vel(1,iatm(ih),icel)   &
                               + (vel_scale(1)-1.0_wp)*vel_cm(1)
          vel(2,iatm(ih),icel) = vel(2,iatm(ih),icel)   &
                               + (vel_scale(2)-1.0_wp)*vel_cm(2)
          vel(3,iatm(ih),icel) = vel(3,iatm(ih),icel)   &
                               + (vel_scale(3)-1.0_wp)*vel_cm(3)
        end do
      end do
    end do

    do icel = id+1, ncell, nthread
      do ix = 1, natom(icel)
        factor = half_dt * imass(ix,icel)
        vel(1,ix,icel) = vel(1,ix,icel) + factor*force(1,ix,icel)
        vel(2,ix,icel) = vel(2,ix,icel) + factor*force(2,ix,icel)
        vel(3,ix,icel) = vel(3,ix,icel) + factor*force(3,ix,icel)
      end do
    end do

    !$omp end parallel

    return

  end subroutine update_vel_vv1_group

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_coord_vv1_group
  !> @brief        VV1 with group pressures
  !! @authors      JJ
  !! @param[in]    constraints : constraints information
  !! @param[in]    ncell       : number of cells
  !! @param[in]    natom       : number of atoms
  !! @param[in]    nwater      : number of water
  !! @param[in]    water_list  : water molecule list
  !! @param[in]    mass        : mass
  !! @param[in]    force       : force
  !! @param[in]    coord_ref   : reference coordinates
  !! @param[in]    vel_ref     : reference velocities
  !! @param[in]    size_scale  : size scale
  !! @param[in]    vel_scale   : velocity scale
  !! @param[in]    dt          : time step
  !! @param[in]    half_dt     : half time step
  !! @param[inout] coord       : coordinates
  !! @param[inout] vel         : velocities
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_coord_vv1_group(constraints, ncell, natom, nwater,    &
                                    water_list, mass, force, coord_ref,   &
                                    vel_ref, size_scale, vel_scale,       &
                                    dt, half_dt, coord, vel)

    ! formal arguments
    type(s_constraints), target, intent(in)    :: constraints
    integer,                     intent(in)    :: ncell, natom(:)
    integer,                     intent(in)    :: nwater(:), water_list(:,:,:)
    real(wip),                   intent(in)    :: mass(:,:)
    real(wip),                   intent(in)    :: force(:,:,:)
    real(wip),                   intent(in)    :: coord_ref(:,:,:)
    real(wip),                   intent(in)    :: vel_ref(:,:,:)
    real(wip),                   intent(in)    :: size_scale(:)
    real(wip),                   intent(in)    :: vel_scale(:)
    real(wip),                   intent(in)    :: dt, half_dt
    real(wip),                   intent(inout) :: coord(:,:,:)
    real(wip),                   intent(inout) :: vel(:,:,:)

    ! local variables
    integer                      :: icel, i, ix, j, k, ih, connect
    integer                      :: id, omp_get_thread_num
    integer                      :: iatm(1:8), water_atom
    real(wip)                    :: factor
    real(wip)                    :: total_mass
    real(wip)                    :: r_cm(1:3), vel_cm(1:3)
    real(wip)                    :: r(1:3,1:8), v(1:3,1:8)

    integer,             pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)
    integer,             pointer :: nsolute(:)


    HGr_local       => constraints%HGr_local
    HGr_bond_list   => constraints%HGr_bond_list
    nsolute         => constraints%No_HGr
    connect         =  constraints%connect

    water_atom = constraints%water_type

    !$omp parallel private(id, icel, j, k, total_mass, r_cm, vel_cm, &
    !$omp                  iatm, ih, r, v, ix, factor)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do icel = id+1, ncell, nthread

      ! No group
      !
      do ix = 1, nsolute(icel)
        coord(1,ix,icel) = size_scale(1)*coord_ref(1,ix,icel)
        coord(2,ix,icel) = size_scale(2)*coord_ref(2,ix,icel)
        coord(3,ix,icel) = size_scale(3)*coord_ref(3,ix,icel)
      end do

      ! XHn group
      !
!ocl nosimd
      do j = 1, connect
        do k = 1, HGr_local(j,icel)
          total_mass  = 0.0_dp
          r_cm(1)     = 0.0_dp
          r_cm(2)     = 0.0_dp
          r_cm(3)     = 0.0_dp
          do ih = 1, j+1
            iatm(ih)    = HGr_bond_list(ih,k,j,icel)
            total_mass  = total_mass + mass(iatm(ih),icel)
            r_cm(1)     = r_cm(1)  &
                        + mass(iatm(ih),icel)*coord_ref(1,iatm(ih),icel)
            r_cm(2)     = r_cm(2)  &
                        + mass(iatm(ih),icel)*coord_ref(2,iatm(ih),icel)
            r_cm(3)     = r_cm(3)  &
                        + mass(iatm(ih),icel)*coord_ref(3,iatm(ih),icel)
          end do
          r_cm(1)   = r_cm(1) / total_mass
          r_cm(2)   = r_cm(2) / total_mass
          r_cm(3)   = r_cm(3) / total_mass
          do ih = 1, j+1
            coord(1,iatm(ih),icel) = coord_ref(1,iatm(ih),icel) &
                                   + (size_scale(1)-1.0_wp)*r_cm(1)
            coord(2,iatm(ih),icel) = coord_ref(2,iatm(ih),icel) &
                                   + (size_scale(2)-1.0_wp)*r_cm(2)
            coord(3,iatm(ih),icel) = coord_ref(3,iatm(ih),icel) &
                                   + (size_scale(3)-1.0_wp)*r_cm(3)
          end do

        end do
      end do

      ! Water
      !
!ocl nosimd
      do ix = 1, nwater(icel)
        iatm(1:water_atom) = water_list(1:water_atom,ix,icel)
        total_mass  = 0.0_dp
        r_cm(1)     = 0.0_dp
        r_cm(2)     = 0.0_dp
        r_cm(3)     = 0.0_dp
        do ih = 1, water_atom
          total_mass  = total_mass + mass(iatm(ih),icel)
          r_cm(1)     = r_cm(1)  &
                      + mass(iatm(ih),icel)*coord_ref(1,iatm(ih),icel)
          r_cm(2)     = r_cm(2)  &
                      + mass(iatm(ih),icel)*coord_ref(2,iatm(ih),icel)
          r_cm(3)     = r_cm(3)  &
                      + mass(iatm(ih),icel)*coord_ref(3,iatm(ih),icel)
        end do
        r_cm(1)  = r_cm(1) / total_mass
        r_cm(2)  = r_cm(2) / total_mass
        r_cm(3)  = r_cm(3) / total_mass
        do ih = 1, water_atom
          coord(1,iatm(ih),icel) = coord_ref(1,iatm(ih),icel) &
                                 + (size_scale(1)-1.0_wp)*r_cm(1)
          coord(2,iatm(ih),icel) = coord_ref(2,iatm(ih),icel) &
                                 + (size_scale(2)-1.0_wp)*r_cm(2)
          coord(3,iatm(ih),icel) = coord_ref(3,iatm(ih),icel) &
                                 + (size_scale(3)-1.0_wp)*r_cm(3)
        end do
      end do
     end do

    do icel = id+1, ncell, nthread
      do ix = 1, natom(icel)
        coord(1,ix,icel) = coord(1,ix,icel)+vel(1,ix,icel)*half_dt
        coord(2,ix,icel) = coord(2,ix,icel)+vel(2,ix,icel)*half_dt
        coord(3,ix,icel) = coord(3,ix,icel)+vel(3,ix,icel)*half_dt
      end do
    end do

    !$omp end parallel

    return

  end subroutine update_coord_vv1_group

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_coord_vv2_group
  !> @brief        VV1 with group pressures
  !! @authors      JJ
  !! @param[in]    constraints : constraints information
  !! @param[in]    ncell       : number of cells
  !! @param[in]    natom       : number of atoms
  !! @param[in]    nwater      : number of water
  !! @param[in]    water_list  : water molecule list
  !! @param[in]    mass        : mass
  !! @param[in]    force       : force
  !! @param[in]    coord_ref   : reference coordinates
  !! @param[in]    vel_ref     : reference velocities
  !! @param[in]    size_scale  : size scale
  !! @param[in]    vel_scale   : velocity scale
  !! @param[in]    dt          : time step
  !! @param[in]    half_dt     : half time step
  !! @param[inout] coord       : coordinates
  !! @param[inout] vel         : velocities
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_coord_vv2_group(constraints, ncell, natom, nwater,    &
                                    water_list, mass, force, coord_ref,   &
                                    vel_ref, size_scale, vel_scale,       &
                                    dt, half_dt, coord, vel)

    ! formal arguments
    type(s_constraints), target, intent(in)    :: constraints
    integer,                     intent(in)    :: ncell, natom(:)
    integer,                     intent(in)    :: nwater(:), water_list(:,:,:)
    real(wip),                   intent(in)    :: mass(:,:)
    real(wip),                   intent(in)    :: force(:,:,:)
    real(wip),                   intent(in)    :: coord_ref(:,:,:)
    real(wip),                   intent(in)    :: vel_ref(:,:,:)
    real(wip),                   intent(in)    :: size_scale(:)
    real(wip),                   intent(in)    :: vel_scale(:)
    real(wip),                   intent(in)    :: dt, half_dt
    real(wip),                   intent(inout) :: coord(:,:,:)
    real(wip),                   intent(inout) :: vel(:,:,:)

    ! local variables
    integer                      :: icel, i, ix, j, k, ih, connect
    integer                      :: id, omp_get_thread_num
    integer                      :: iatm(1:8), water_atom
    real(wp)                     :: factor
    real(wp)                     :: total_mass
    real(wp)                     :: r_cm(1:3), vel_cm(1:3)
    real(wp)                     :: r(1:3,1:8), v(1:3,1:8)

    integer,             pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)
    integer,             pointer :: nsolute(:)


    HGr_local       => constraints%HGr_local
    HGr_bond_list   => constraints%HGr_bond_list
    nsolute         => constraints%No_HGr
    connect         =  constraints%connect

    water_atom = constraints%water_type

    !$omp parallel private(id, icel, j, k, total_mass, r_cm, vel_cm, &
    !$omp                  iatm, ih, r, v, ix, factor)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do icel = id+1, ncell, nthread
      do ix = 1, natom(icel)
        coord(1,ix,icel) = coord_ref(1,ix,icel)+vel(1,ix,icel)*half_dt
        coord(2,ix,icel) = coord_ref(2,ix,icel)+vel(2,ix,icel)*half_dt
        coord(3,ix,icel) = coord_ref(3,ix,icel)+vel(3,ix,icel)*half_dt
      end do
    end do

    do icel = id+1, ncell, nthread

      ! No group
      !
      do ix = 1, nsolute(icel)
        coord(1,ix,icel) = size_scale(1)*coord(1,ix,icel)
        coord(2,ix,icel) = size_scale(2)*coord(2,ix,icel)
        coord(3,ix,icel) = size_scale(3)*coord(3,ix,icel)
      end do

      ! XHn group
      !
!ocl nosimd
      do j = 1, connect
        do k = 1, HGr_local(j,icel)
          total_mass  = 0.0_dp
          r_cm(1)     = 0.0_dp
          do ih = 1, j+1
            iatm(ih)    = HGr_bond_list(ih,k,j,icel)
            total_mass  = total_mass + mass(iatm(ih),icel)
            r_cm(1)     = r_cm(1)  &
                        + mass(iatm(ih),icel)*coord(1,iatm(ih),icel)
            r_cm(2)     = r_cm(2)  &
                        + mass(iatm(ih),icel)*coord(2,iatm(ih),icel)
            r_cm(3)     = r_cm(3)  &
                        + mass(iatm(ih),icel)*coord(3,iatm(ih),icel)
          end do
          r_cm(1)   = r_cm(1) / total_mass
          r_cm(2)   = r_cm(2) / total_mass
          r_cm(3)   = r_cm(3) / total_mass
          do ih = 1, j+1
            coord(1,iatm(ih),icel) = coord(1,iatm(ih),icel) &
                                   + (size_scale(1)-1.0_wp)*r_cm(1)
            coord(2,iatm(ih),icel) = coord(2,iatm(ih),icel) &
                                   + (size_scale(2)-1.0_wp)*r_cm(2)
            coord(3,iatm(ih),icel) = coord(3,iatm(ih),icel) &
                                   + (size_scale(3)-1.0_wp)*r_cm(3)
          end do

        end do
      end do

      ! Water
      !
!ocl nosimd
      do ix = 1, nwater(icel)
        iatm(1:water_atom) = water_list(1:water_atom,ix,icel)
        total_mass  = 0.0_dp
        r_cm(1)     = 0.0_dp
        r_cm(2)     = 0.0_dp
        r_cm(3)     = 0.0_dp
        do ih = 1, water_atom
          total_mass  = total_mass + mass(iatm(ih),icel)
          r_cm(1)     = r_cm(1)  &
                      + mass(iatm(ih),icel)*coord(1,iatm(ih),icel)
          r_cm(2)     = r_cm(2)  &
                      + mass(iatm(ih),icel)*coord(2,iatm(ih),icel)
          r_cm(3)     = r_cm(3)  &
                      + mass(iatm(ih),icel)*coord(3,iatm(ih),icel)
        end do
        r_cm(1)  = r_cm(1) / total_mass
        r_cm(2)  = r_cm(2) / total_mass
        r_cm(3)  = r_cm(3) / total_mass
        do ih = 1, water_atom
          coord(1,iatm(ih),icel) = coord(1,iatm(ih),icel) &
                                 + (size_scale(1)-1.0_wp)*r_cm(1)
          coord(2,iatm(ih),icel) = coord(2,iatm(ih),icel) &
                                 + (size_scale(2)-1.0_wp)*r_cm(2)
          coord(3,iatm(ih),icel) = coord(3,iatm(ih),icel) &
                                 + (size_scale(3)-1.0_wp)*r_cm(3)
        end do
      end do
    end do

    !$omp end parallel

    return

  end subroutine update_coord_vv2_group

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_settle
  !> @brief        SETTLE for three-site water
  !! @authors      JJ
  !! @param[in]    vel_update  : flag for update velocity or not
  !! @param[in]    dt          : time step
  !! @param[in]    coord_old   : reference coordinates
  !! @param[in]    domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] coord       : coordinates
  !! @param[inout] vel         : velocities
  !! @param[inout] virial      : virial of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_settle(vel_update, dt, coord_old, domain, &
                            constraints, coord, vel, virial)

    ! formal arguments
    logical,                 intent(in)    :: vel_update
    real(wip),               intent(in)    :: dt
    real(wip),               intent(in)    :: coord_old(:,:,:)
    type(s_domain),  target, intent(in)    :: domain
    type(s_constraints),     intent(inout) :: constraints
    real(wip),               intent(inout) :: coord(:,:,:)
    real(wip),               intent(inout) :: vel(:,:,:)
    real(dp ),               intent(inout) :: virial(3,3)

    ! local variables
    real(dp)                 :: rHH, rOH, massO, massH
    real(dp)                 :: mass(1:3), mass_H2O, ra, rb, rc, inv_ra, inv_dt
    real(dp)                 :: x0(1:3,1:3), x1(1:3,1:3), x3(1:3,1:3)
    real(dp)                 :: delt(1:3,1:3), rf(1:3,1:3,1:3)
    real(dp)                 :: com0(1:3), com1(1:3)
    real(dp)                 :: oh1(1:3), oh2(1:3)
    real(dp)                 :: Xaxis(1:3), Yaxis(1:3), Zaxis(1:3)
    real(dp)                 :: Xaxis2(1:3), Yaxis2(1:3), Zaxis2(1:3)
    real(dp)                 :: mtrx(1:3,1:3)
    real(dp)                 :: xp0(1:3,1:3), xp1(1:3,1:3), xp2(1:3,1:3)
    real(dp)                 :: xp3(1:3,1:3)
    real(dp)                 :: dxp21(1:2), dxp23(1:2), dxp31(1:2)
    real(dp)                 :: sin_phi, cos_phi, sin_psi, cos_psi, tmp
    real(dp)                 :: rb_cosphi, rb_sinphi
    real(dp)                 :: rc_sinpsi_sinphi, rc_sinpsi_cosphi
    real(dp)                 :: alpha, beta, gamma, al2bt2, sin_theta, cos_theta
    real(dp)                 :: viri_local(1:3,1:3)
    real(dp)                 :: virial_omp(3,3,nthread)
    integer                  :: i, ix, iatom(1:3), i1, j1, k1
    integer                  :: id, omp_get_thread_num

    integer,         pointer :: nwater(:), water_list(:,:,:), ncell


    nwater     => domain%num_water
    water_list => domain%water_list
    ncell      => domain%num_cell_local

    rOH        =  constraints%water_rOH
    rHH        =  constraints%water_rHH
    massO      =  constraints%water_massO
    massH      =  constraints%water_massH

    mass(1)    =  massO
    mass(2)    =  massH
    mass(3)    =  massH
    mass_H2O   =  massO + 2.0_dp*massH
    rc         =  0.5_dp * rHH 
    ra         =  2.0_dp * sqrt(rOH*rOH-rc*rc) * massH / mass_H2O
    rb         =  sqrt(rOH*rOH-rc*rc) - ra
    inv_ra     =  1.0_dp / ra
    inv_dt     =  1.0_dp / dt

    virial_omp(1:3,1:3,1:nthread) = 0.0_dp

    !$omp parallel default(shared)                                             &
    !$omp private(i, ix, iatom, i1, j1, k1, x0, x1, x3, com0, com1, oh1, oh2,  &
    !$omp         delt, rf, Xaxis, Yaxis, Zaxis, Xaxis2, Yaxis2, Zaxis2,       &
    !$omp         mtrx, xp0, xp1, xp2, xp3, tmp, sin_phi, cos_phi, sin_psi,    &
    !$omp         cos_psi, rb_cosphi, rb_sinphi, rc_sinpsi_sinphi,             &
    !$omp         rc_sinpsi_cosphi, alpha, beta, gamma, dxp21, dxp23, dxp31,   &
    !$omp         al2bt2, sin_theta, cos_theta, viri_local, id)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    viri_local(1:3,1:3) = 0.0_dp

    do i = id+1, ncell, nthread

!ocl nosimd
      do ix = 1, nwater(i)

        iatom(1:3) = water_list(1:3,ix,i)
        com0(1:3)   = 0.0_dp
        com1(1:3)  = 0.0_dp
        do k1 = 1, 3
          x0(1:3,k1) = coord_old(1:3,iatom(k1),i)
          x1(1:3,k1) = coord    (1:3,iatom(k1),i)
          com0(1:3)  = com0(1:3) + x0(1:3,k1)*mass(k1)
          com1(1:3)  = com1(1:3) + x1(1:3,k1)*mass(k1)
        end do
        com0(1:3) = com0(1:3) / mass_H2O
        com1(1:3) = com1(1:3) / mass_H2O
  
        do k1 = 1, 3
          x0(1:3,k1) = x0(1:3,k1) - com0(1:3)
          x1(1:3,k1) = x1(1:3,k1) - com1(1:3)
        end do
  
        oh1(1:3) = x0(1:3,2) - x0(1:3,1)
        oh2(1:3) = x0(1:3,3) - x0(1:3,1)
  
        Zaxis(1) = oh1(2)*oh2(3) - oh1(3)*oh2(2)
        Zaxis(2) = oh1(3)*oh2(1) - oh1(1)*oh2(3)
        Zaxis(3) = oh1(1)*oh2(2) - oh1(2)*oh2(1)
        Xaxis(1) = x1(2,1)*Zaxis(3) - x1(3,1)*Zaxis(2)
        Xaxis(2) = x1(3,1)*Zaxis(1) - x1(1,1)*Zaxis(3)
        Xaxis(3) = x1(1,1)*Zaxis(2) - x1(2,1)*Zaxis(1)
        Yaxis(1) = Zaxis(2)*Xaxis(3) - Zaxis(3)*Xaxis(2)
        Yaxis(2) = Zaxis(3)*Xaxis(1) - Zaxis(1)*Xaxis(3)
        Yaxis(3) = Zaxis(1)*Xaxis(2) - Zaxis(2)*Xaxis(1)
  
        Xaxis2(1:3) = Xaxis(1:3) * Xaxis(1:3)
        Yaxis2(1:3) = Yaxis(1:3) * Yaxis(1:3)
        Zaxis2(1:3) = Zaxis(1:3) * Zaxis(1:3)
        mtrx(1:3,1) = Xaxis(1:3) / sqrt(Xaxis2(1)+Xaxis2(2)+Xaxis2(3))
        mtrx(1:3,2) = Yaxis(1:3) / sqrt(Yaxis2(1)+Yaxis2(2)+Yaxis2(3))
        mtrx(1:3,3) = Zaxis(1:3) / sqrt(Zaxis2(1)+Zaxis2(2)+Zaxis2(3))

        xp0(1:3,1:3) = 0.0_dp
        xp1(1:3,1:3) = 0.0_dp
        do i1 = 1, 3
          do k1 = 1, 3
            do j1 = 1, 3
              xp0(k1,i1) = xp0(k1,i1) + mtrx(j1,i1)*x0(j1,k1)
              xp1(k1,i1) = xp1(k1,i1) + mtrx(j1,i1)*x1(j1,k1)
            end do
          end do
        end do

        sin_phi = xp1(1,3)*inv_ra
        tmp = 1.0_dp - sin_phi*sin_phi
        if (tmp > 0.0_dp) then
          cos_phi = sqrt(tmp)
        else
          cos_phi = 0.0_dp
        end if

        sin_psi = (xp1(2,3) - xp1(3,3))/(rHH*cos_phi)
        tmp = 1.0_dp - sin_psi*sin_psi
        if (tmp > 0.0_dp) then
          cos_psi = sqrt(tmp)
        else
          cos_psi = 0.0_dp
        end if

        rb_cosphi = rb * cos_phi
        rb_sinphi = rb * sin_phi
        rc_sinpsi_sinphi = rc * sin_psi * sin_phi
        rc_sinpsi_cosphi = rc * sin_psi * cos_phi
  
        xp2(1,2) =   ra * cos_phi
        xp2(1,3) =   ra * sin_phi
        xp2(2,1) = - rc * cos_psi
        xp2(2,2) = - rb_cosphi - rc_sinpsi_sinphi
        xp2(2,3) = - rb_sinphi + rc_sinpsi_cosphi
        xp2(3,1) = - xp2(2,1)
        xp2(3,2) = - rb_cosphi + rc_sinpsi_sinphi
        xp2(3,3) = - rb_sinphi - rc_sinpsi_cosphi
  
        dxp21(1:2) = xp0(2,1:2) - xp0(1,1:2)
        dxp23(1:2) = xp0(2,1:2) - xp0(3,1:2)
        dxp31(1:2) = xp0(3,1:2) - xp0(1,1:2)
        alpha =   xp2(2,1)*dxp23(1) + dxp21(2)*xp2(2,2) + dxp31(2)*xp2(3,2)
        beta  = - xp2(2,1)*dxp23(2) + dxp21(1)*xp2(2,2) + dxp31(1)*xp2(3,2)
        gamma =   dxp21(1) * xp1(2,2) - xp1(2,1) * dxp21(2) &
                + dxp31(1) * xp1(3,2) - xp1(3,1) * dxp31(2)
  
        al2bt2 = alpha*alpha + beta*beta
        sin_theta = (alpha*gamma - beta*sqrt(al2bt2 - gamma*gamma))/al2bt2
        tmp = 1.0_dp - sin_theta*sin_theta
        if (tmp > 0.0_dp) then
          cos_theta = sqrt(tmp)
        else
          cos_theta = 0.0_dp
        end if

        xp3(1,1) = - xp2(1,2)*sin_theta
        xp3(1,2) =   xp2(1,2)*cos_theta
        xp3(1,3) =   xp2(1,3)
        xp3(2,1) =   xp2(2,1)*cos_theta - xp2(2,2)*sin_theta
        xp3(2,2) =   xp2(2,1)*sin_theta + xp2(2,2)*cos_theta
        xp3(2,3) =   xp2(2,3)
        xp3(3,1) =   xp2(3,1)*cos_theta - xp2(3,2)*sin_theta
        xp3(3,2) =   xp2(3,1)*sin_theta + xp2(3,2)*cos_theta
        xp3(3,3) =   xp2(3,3)
  
        x3(1:3,1:3) = 0.0_dp
        do k1 = 1, 3
          do i1 = 1, 3
            do j1 = 1, 3
              x3(i1,k1) = x3(i1,k1) + mtrx(i1,j1)*xp3(k1,j1)
            end do
          end do
          coord(1:3,iatom(k1),i) = real(x3(1:3,k1) + com1(1:3),wip)
        end do
  
        do k1 = 1, 3
          delt(1:3,k1) = x3(1:3,k1) - x1(1:3,k1)
          do i1 = 1, 3
            do j1 = 1, 3
              rf(i1,j1,k1) = - coord_old(i1,iatom(k1),i)*delt(j1,k1)*mass(k1)
            end do
          end do
        end do
  
        do i1 = 1, 3
          do j1 = 1, 3
            viri_local(1:3,j1) = viri_local(1:3,j1) + rf(1:3,j1,i1)
          end do
        end do
  
        if (vel_update) then
          do k1 = 1, 3
            vel(1:3,iatom(k1),i) = vel(1:3,iatom(k1),i) &
                                 + real(delt(1:3,k1)*inv_dt,wip)
          end do
        end if

      end do
  
    end do

    do k1 = 1, 3
      virial_omp(1:3,k1,id+1) = virial_omp(1:3,k1,id+1) - viri_local(1:3,k1)
    end do

    !$omp end parallel

    do i = 1, nthread
      virial(1:3,1:3) = virial(1:3,1:3) + virial_omp(1:3,1:3,i)
    end do

    return

  end subroutine compute_settle

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_settle_TIP2
  !> @brief        bond constraints for two-site water
  !! @authors      JJ
  !! @param[in]    vel_update  : flag for update velocity or not
  !! @param[in]    dt          : time step
  !! @param[in]    coord_old   : reference coordinates
  !! @param[in]    domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] coord       : coordinates
  !! @param[inout] vel         : velocities
  !! @param[inout] virial      : virial of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_settle_TIP2(vel_update, dt, coord_old, domain, &
                                 constraints, coord, vel, virial)

    ! formal arguments
    logical,                 intent(in)    :: vel_update
    real(wip),               intent(in)    :: dt
    real(wip),               intent(in)    :: coord_old(:,:,:)
    type(s_domain),  target, intent(in)    :: domain
    type(s_constraints),     intent(inout) :: constraints
    real(wip),               intent(inout) :: coord(:,:,:)
    real(wip),               intent(inout) :: vel(:,:,:)
    real(dp ),               intent(inout) :: virial(3,3)

    ! local variables
    real(dp)                 :: x0(3), x1(3), A, B, C, lambda
    real(dp)                 :: r12, mass1, mass2, imass1, imass2, imass
    real(dp)                 :: dcoord(3,2), dvel(3,2)
    real(dp)                 :: virial_omp(3,nthread)
    integer                  :: i, ix, list1, list2, id, omp_get_thread_num

    integer,         pointer :: nwater(:), water_list(:,:,:), ncell


    nwater     => domain%num_water
    water_list => domain%water_list
    ncell      => domain%num_cell_local

    r12        =  constraints%water_r12
    mass1      =  constraints%water_mass1
    mass2      =  constraints%water_mass2
    imass1     =  1.0_wp / mass1
    imass2     =  1.0_wp / mass2
    imass      = imass1 + imass2

    virial_omp(1:3,1:nthread) = 0.0_dp

    !$omp parallel default(shared)                                             &
    !$omp private(i, ix, list1, list2, x0, x1, A, B, C, lambda, dcoord, dvel,  &
    !$omp         id)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell, nthread

!ocl nosimd
      do ix = 1, nwater(i)

        list1 = water_list(1,ix,i)
        list2 = water_list(2,ix,i)
        x0(1:3) = coord_old(1:3,list1,i) - coord_old(1:3,list2,i)
        x1(1:3) = coord(1:3,list1,i) - coord(1:3,list2,i)
        A = imass*imass*(x0(1)*x0(1)+x0(2)*x0(2)+x0(3)*x0(3))
        B = imass      *(x0(1)*x1(1)+x0(2)*x1(2)+x0(3)*x1(3))
        C = x1(1)*x1(1)+x1(2)*x1(2)+x1(3)*x1(3) - r12*r12
        lambda = -B + sqrt(B*B-A*C)
        lambda = lambda / A

        dcoord(1:3,1) = lambda*x0(1:3)*imass1
        dcoord(1:3,2) = lambda*x0(1:3)*imass2 
        dvel  (1:3,1) = dcoord(1:3,1) / real(dt,dp)
        dvel  (1:3,2) = dcoord(1:3,2) / real(dt,dp)
        coord(1:3,list1,i) = coord(1:3,list1,i) + dcoord(1:3,1)
        coord(1:3,list2,i) = coord(1:3,list2,i) - dcoord(1:3,2)

        if (vel_update) then
          vel(1:3,list1,i) = vel(1:3,list1,i) + dvel(1:3,1)
          vel(1:3,list2,i) = vel(1:3,list2,i) - dvel(1:3,2)
        end if

        virial_omp(1,id+1) = virial_omp(1,id+1) + lambda*x0(1)*x0(1)
        virial_omp(2,id+1) = virial_omp(2,id+1) + lambda*x0(2)*x0(2)
        virial_omp(3,id+1) = virial_omp(3,id+1) + lambda*x0(3)*x0(3)
      end do
    end do

    !$omp end parallel

    do i = 1, nthread
      virial(1,1) = virial(1,1) + virial_omp(1,i)
      virial(2,2) = virial(2,2) + virial_omp(2,i)
      virial(3,3) = virial(3,3) + virial_omp(3,i)
    end do

    return

  end subroutine compute_settle_TIP2

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_settle_TIP2_vv2
  !> @brief        bond constraints for two-site water
  !! @authors      JJ
  !! @param[in]    domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] coord       : coordinates
  !! @param[inout] vel         : velocities
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_settle_TIP2_vv2(domain, constraints, coord, vel)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_constraints),     intent(inout) :: constraints
    real(wip),               intent(inout) :: coord(:,:,:)
    real(wip),               intent(inout) :: vel(:,:,:)

    ! local variables
    real(dp)                 :: x0(3), x1(3), A, B, lambda
    real(dp)                 :: r12, mass1, mass2, imass1, imass2, imass
    real(dp)                 :: dcoord(3,2), dvel(3,2)
    real(dp)                 :: virial_omp(3,nthread)
    integer                  :: i, ix, list1, list2, id, omp_get_thread_num

    integer,         pointer :: nwater(:), water_list(:,:,:), ncell


    nwater     => domain%num_water
    water_list => domain%water_list
    ncell      => domain%num_cell_local

    mass1      =  constraints%water_mass1
    mass2      =  constraints%water_mass2
    imass1     =  1.0_wp / mass1
    imass2     =  1.0_wp / mass2
    imass      = imass1 + imass2

    virial_omp(1:3,1:nthread) = 0.0_dp

    !$omp parallel default(shared)                                             &
    !$omp private(i, ix, list1, list2, x0, x1, A, B, lambda, dcoord, dvel,  &
    !$omp         id)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell, nthread

!ocl nosimd
      do ix = 1, nwater(i)

        list1 = water_list(1,ix,i)
        list2 = water_list(2,ix,i)
        x0(1:3) = coord(1:3,list1,i) - coord(1:3,list2,i)
        x1(1:3) = vel  (1:3,list1,i) - vel  (1:3,list2,i)
        A = x0(1)*x1(1) + x0(2)*x1(2) + x0(3)*x1(3)
        B = imass*(x0(1)*x0(1)+x0(2)*x0(2)+x0(3)*x0(3))
        lambda = -A/B

        dvel(1:3,1) =  lambda*x0(1:3)*imass1
        dvel(1:3,2) =  lambda*x0(1:3)*imass2

        vel(1:3,list1,i) = vel(1:3,list1,i) + dvel(1:3,1)
        vel(1:3,list2,i) = vel(1:3,list2,i) - dvel(1:3,2)

      end do
    end do

    !$omp end parallel

    return

  end subroutine compute_settle_TIP2_vv2

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_settle_min
  !> @brief        SETTLE for three-site water (only coordinaets)
  !! @authors      JJ
  !! @param[in]    coord_old   : reference coordinates
  !! @param[in]    domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] coord       : coordinates
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_settle_min(coord_old, domain, constraints, coord)

    ! formal arguments
    real(wip),               intent(in)    :: coord_old(:,:,:)
    type(s_domain),  target, intent(in)    :: domain
    type(s_constraints),     intent(inout) :: constraints
    real(wip),               intent(inout) :: coord(:,:,:)

    ! local variables
    real(dp)                 :: rHH, rOH, massO, massH
    real(dp)                 :: mass(1:3), mass_H2O, ra, rb, rc, inv_ra, inv_dt
    real(dp)                 :: x0(1:3,1:3), x1(1:3,1:3), x3(1:3,1:3)
    real(dp)                 :: delt(1:3,1:3), rf(1:3,1:3,1:3)
    real(dp)                 :: com0(1:3), com1(1:3)
    real(dp)                 :: oh1(1:3), oh2(1:3)
    real(dp)                 :: Xaxis(1:3), Yaxis(1:3), Zaxis(1:3)
    real(dp)                 :: Xaxis2(1:3), Yaxis2(1:3), Zaxis2(1:3)
    real(dp)                 :: mtrx(1:3,1:3)
    real(dp)                 :: xp0(1:3,1:3), xp1(1:3,1:3), xp2(1:3,1:3)
    real(dp)                 :: xp3(1:3,1:3)
    real(dp)                 :: dxp21(1:2), dxp23(1:2), dxp31(1:2)
    real(dp)                 :: sin_phi, cos_phi, sin_psi, cos_psi, tmp
    real(dp)                 :: rb_cosphi, rb_sinphi
    real(dp)                 :: rc_sinpsi_sinphi, rc_sinpsi_cosphi
    real(dp)                 :: alpha, beta, gamma, al2bt2, sin_theta, cos_theta
    real(dp)                 :: viri_local(1:3,1:3)
    real(dp)                 :: virial_omp(3,3,nthread)
    integer                  :: i, ix, iatom(1:3), i1, j1, k1
    integer                  :: id, omp_get_thread_num

    integer,         pointer :: nwater(:), water_list(:,:,:), ncell


    nwater     => domain%num_water
    water_list => domain%water_list
    ncell      => domain%num_cell_local

    rOH        =  constraints%water_rOH
    rHH        =  constraints%water_rHH
    massO      =  constraints%water_massO
    massH      =  constraints%water_massH

    mass(1)    =  massO
    mass(2)    =  massH
    mass(3)    =  massH
    mass_H2O   =  massO + 2.0_dp*massH
    rc         =  0.5_dp * rHH 
    ra         =  2.0_dp * dsqrt(rOH*rOH-rc*rc) * massH / mass_H2O
    rb         =  dsqrt(rOH*rOH-rc*rc) - ra
    inv_ra     =  1.0_dp / ra

    virial_omp(1:3,1:3,1:nthread) = 0.0_dp

    !$omp parallel default(shared)                                             &
    !$omp private(i, ix, iatom, i1, j1, k1, x0, x1, x3, com0, com1, oh1, oh2,  &
    !$omp         delt, rf, Xaxis, Yaxis, Zaxis, Xaxis2, Yaxis2, Zaxis2,       &
    !$omp         mtrx, xp0, xp1, xp2, xp3, tmp, sin_phi, cos_phi, sin_psi,    &
    !$omp         cos_psi, rb_cosphi, rb_sinphi, rc_sinpsi_sinphi,             &
    !$omp         rc_sinpsi_cosphi, alpha, beta, gamma, dxp21, dxp23, dxp31,   &
    !$omp         al2bt2, sin_theta, cos_theta, viri_local, id)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell, nthread
      do ix = 1, nwater(i)

      iatom(1:3) = water_list(1:3,ix,i)
      com0(1:3)   = 0.0_dp
      com1(1:3)  = 0.0_dp

      do k1 = 1, 3
        x0(1:3,k1) = coord_old(1:3,iatom(k1),i)
        x1(1:3,k1) = coord    (1:3,iatom(k1),i)
        com0(1:3)  = com0(1:3) + x0(1:3,k1)*mass(k1)
        com1(1:3)  = com1(1:3) + x1(1:3,k1)*mass(k1)
      end do
      com0(1:3) = com0(1:3) / mass_H2O
      com1(1:3) = com1(1:3) / mass_H2O

      do k1 = 1, 3
        x0(1:3,k1) = x0(1:3,k1) - com0(1:3)
        x1(1:3,k1) = x1(1:3,k1) - com1(1:3)
      end do

      oh1(1:3) = x0(1:3,2) - x0(1:3,1)
      oh2(1:3) = x0(1:3,3) - x0(1:3,1)

      Zaxis(1) = oh1(2)*oh2(3) - oh1(3)*oh2(2)
      Zaxis(2) = oh1(3)*oh2(1) - oh1(1)*oh2(3)
      Zaxis(3) = oh1(1)*oh2(2) - oh1(2)*oh2(1)
      Xaxis(1) = x1(2,1)*Zaxis(3) - x1(3,1)*Zaxis(2)
      Xaxis(2) = x1(3,1)*Zaxis(1) - x1(1,1)*Zaxis(3)
      Xaxis(3) = x1(1,1)*Zaxis(2) - x1(2,1)*Zaxis(1)
      Yaxis(1) = Zaxis(2)*Xaxis(3) - Zaxis(3)*Xaxis(2)
      Yaxis(2) = Zaxis(3)*Xaxis(1) - Zaxis(1)*Xaxis(3)
      Yaxis(3) = Zaxis(1)*Xaxis(2) - Zaxis(2)*Xaxis(1)

      Xaxis2(1:3) = Xaxis(1:3) * Xaxis(1:3)
      Yaxis2(1:3) = Yaxis(1:3) * Yaxis(1:3)
      Zaxis2(1:3) = Zaxis(1:3) * Zaxis(1:3)
      mtrx(1:3,1) = Xaxis(1:3) / dsqrt(Xaxis2(1)+Xaxis2(2)+Xaxis2(3))
      mtrx(1:3,2) = Yaxis(1:3) / dsqrt(Yaxis2(1)+Yaxis2(2)+Yaxis2(3))
      mtrx(1:3,3) = Zaxis(1:3) / dsqrt(Zaxis2(1)+Zaxis2(2)+Zaxis2(3))

      xp0(1:3,1:3) = 0.0_dp
      xp1(1:3,1:3) = 0.0_dp
      do i1 = 1, 3
        do k1 = 1, 3
          do j1 = 1, 3
            xp0(k1,i1) = xp0(k1,i1) + mtrx(j1,i1)*x0(j1,k1)
            xp1(k1,i1) = xp1(k1,i1) + mtrx(j1,i1)*x1(j1,k1)
          end do
        end do
      end do

      sin_phi = xp1(1,3)*inv_ra
      tmp = 1.0_dp - sin_phi*sin_phi
      if (tmp > 0.0_dp) then
        cos_phi = dsqrt(tmp)
      else
        cos_phi = 0.0_dp
      end if

      sin_psi = (xp1(2,3) - xp1(3,3))/(rHH*cos_phi)
      tmp = 1.0_dp - sin_psi*sin_psi
      if (tmp > 0.0_dp) then
        cos_psi = dsqrt(tmp)
      else
        cos_psi = 0.0_dp
      end if

      rb_cosphi = rb * cos_phi
      rb_sinphi = rb * sin_phi
      rc_sinpsi_sinphi = rc * sin_psi * sin_phi
      rc_sinpsi_cosphi = rc * sin_psi * cos_phi

      xp2(1,2) =   ra * cos_phi
      xp2(1,3) =   ra * sin_phi
      xp2(2,1) = - rc * cos_psi
      xp2(2,2) = - rb_cosphi - rc_sinpsi_sinphi
      xp2(2,3) = - rb_sinphi + rc_sinpsi_cosphi
      xp2(3,1) = - xp2(2,1)
      xp2(3,2) = - rb_cosphi + rc_sinpsi_sinphi
      xp2(3,3) = - rb_sinphi - rc_sinpsi_cosphi

      dxp21(1:2) = xp0(2,1:2) - xp0(1,1:2)
      dxp23(1:2) = xp0(2,1:2) - xp0(3,1:2)
      dxp31(1:2) = xp0(3,1:2) - xp0(1,1:2)
      alpha =   xp2(2,1)*dxp23(1) + dxp21(2)*xp2(2,2) + dxp31(2)*xp2(3,2)
      beta  = - xp2(2,1)*dxp23(2) + dxp21(1)*xp2(2,2) + dxp31(1)*xp2(3,2)
      gamma =   dxp21(1) * xp1(2,2) - xp1(2,1) * dxp21(2) &
              + dxp31(1) * xp1(3,2) - xp1(3,1) * dxp31(2)

      al2bt2 = alpha*alpha + beta*beta
      sin_theta = (alpha*gamma - beta*dsqrt(al2bt2 - gamma*gamma))/al2bt2
      tmp = 1.0_dp - sin_theta*sin_theta
      if (tmp > 0.0_dp) then
        cos_theta = dsqrt(tmp)
      else
        cos_theta = 0.0_dp
      end if

      xp3(1,1) = - xp2(1,2)*sin_theta
      xp3(1,2) =   xp2(1,2)*cos_theta
      xp3(1,3) =   xp2(1,3)
      xp3(2,1) =   xp2(2,1)*cos_theta - xp2(2,2)*sin_theta
      xp3(2,2) =   xp2(2,1)*sin_theta + xp2(2,2)*cos_theta
      xp3(2,3) =   xp2(2,3)
      xp3(3,1) =   xp2(3,1)*cos_theta - xp2(3,2)*sin_theta
      xp3(3,2) =   xp2(3,1)*sin_theta + xp2(3,2)*cos_theta
      xp3(3,3) =   xp2(3,3)

      x3(1:3,1:3) = 0.0_dp
      do k1 = 1, 3
        do i1 = 1, 3
          do j1 = 1, 3
            x3(i1,k1) = x3(i1,k1) + mtrx(i1,j1)*xp3(k1,j1)
          end do
        end do
        coord(1:3,iatom(k1),i) = x3(1:3,k1) + com1(1:3)
      end do

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_settle_min

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_shake
  !> @brief        SHAKE for rigid bonds
  !! @authors      JJ
  !! @param[in]    vel_update  : flag for update velocity or not
  !! @param[in]    dt          : time step
  !! @param[in]    coord_old   : reference coordinates
  !! @param[in]    domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] coord       : coordinates
  !! @param[inout] vel         : velocities
  !! @param[inout] virial      : virial of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_shake(vel_update, dt, coord_old, domain, &
                            constraints, coord, vel, virial)

    ! formal arguments
    logical,                     intent(in)    :: vel_update
    real(wip),                   intent(in)    :: dt
    real(wip),                   intent(in)    :: coord_old(:,:,:)
    type(s_domain),      target, intent(in)    :: domain
    type(s_constraints), target, intent(inout) :: constraints
    real(wip),                   intent(inout) :: coord(:,:,:)
    real(wip),                   intent(inout) :: vel(:,:,:)
    real(dp ),                   intent(inout) :: virial(3,3)

    ! local variables
    real(dp )                    :: tolerance, imass1, imass2
    real(dp )                    :: x12, y12, z12
    real(dp )                    :: x12_old, y12_old, z12_old
    real(dp )                    :: dist2, r
    real(dp )                    :: factor, g12, g12m1, g12m2, v12m1, v12m2
    real(dp )                    :: viri(1:3), fx, fy, fz
    real(dp )                    :: virial_omp(3,3,nthread)
    real(dp )                    :: coord_dtmp(1:3,1:8), vel_dtmp(1:3,1:8)
    logical                      :: shake_end
    integer                      :: icel, i, j, k, ih, connect, id
    integer                      :: atm1, atm2, iatm(1:8)
    integer                      :: iteration, omp_get_thread_num

    real(wip),           pointer :: HGr_bond_dist(:,:,:,:)
    real(dp ),           pointer :: HGr_bond_vector(:,:,:,:,:)
    real(wip),           pointer :: HGr_shake_force(:,:,:,:)
    real(wip),           pointer :: mass(:,:)
    integer,             pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)
    integer,             pointer :: ncell

    ! FEP
    real(dp)                     :: lambbond

    mass            => domain%mass
    ncell           => domain%num_cell_local

    HGr_local       => constraints%HGr_local
    HGr_bond_list   => constraints%HGr_bond_list 
    HGr_bond_dist   => constraints%HGr_bond_dist
    HGr_bond_vector => constraints%HGr_bond_vector
    HGr_shake_force => constraints%HGr_shake_force

    iteration       =  constraints%shake_iteration
    tolerance       =  real(constraints%shake_tolerance,dp)
    connect         =  constraints%connect

    virial_omp(1:3,1:3,1:nthread) = 0.0_dp

    !$omp parallel default(shared)                                             &
    !$omp private(icel, i, j, k, ih, atm1, atm2, shake_end, r, x12, y12, z12,  &
    !$omp         dist2, imass1, imass2, x12_old, y12_old, z12_old,            &
    !$omp         factor, g12, g12m1, g12m2, v12m1, v12m2, fx, fy, fz,         &
    !$omp         viri, id, coord_dtmp, vel_dtmp, iatm)

#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    ! initialize force and store old bond vector
    !
    do icel = id+1, ncell, nthread
      do j = 1, connect
        do k = 1, HGr_local(j,icel)
          atm1 = HGr_bond_list(1,k,j,icel)
!ocl nosimd
          do ih = 1, j
            atm2 = HGr_bond_list(ih+1,k,j,icel)
            HGr_shake_force(ih,k,j,icel) = 0.0_wip
            HGr_bond_vector(1:3,ih,k,j,icel) = &
              real(coord_old(1:3,atm1,icel)-coord_old(1:3,atm2,icel),dp)
          end do
        end do
      end do
    end do

    ! shake iteration
    !
    shake_end = .true.

    !$omp barrier

    do icel = id+1, ncell, nthread

!ocl nosimd
      do j = 1, connect
        do k = 1, HGr_local(j,icel)

          iatm(1:j+1) = HGr_bond_list(1:j+1,k,j,icel)
          imass1 = 1.0_dp / real(mass(iatm(1),icel),dp)

          coord_dtmp(1:3,1:j+1) = &
            real(coord(1:3,iatm(1:j+1),icel),dp)

          vel_dtmp  (1:3,1:j+1) = &
            real(vel  (1:3,iatm(1:j+1),icel),dp)

          do i = 1, iteration

            shake_end = .true.
            do ih = 1, j

              imass2 = 1.0_dp / real(mass(iatm(ih+1),icel),dp)

              r = real(HGr_bond_dist(ih+1,k,j,icel),dp)

              x12 = coord_dtmp(1,1) - coord_dtmp(1,ih+1)
              y12 = coord_dtmp(2,1) - coord_dtmp(2,ih+1)
              z12 = coord_dtmp(3,1) - coord_dtmp(3,ih+1)
              dist2 = x12**2 + y12**2 + z12**2

              if (abs(dist2 - r*r) >= 2.0_dp*tolerance*r) then 

                shake_end = .false.
                x12_old = HGr_bond_vector(1,ih,k,j,icel)
                y12_old = HGr_bond_vector(2,ih,k,j,icel)
                z12_old = HGr_bond_vector(3,ih,k,j,icel)

                factor = (x12*x12_old + y12*y12_old + z12*z12_old)* &
                         (imass1 + imass2)
                g12    = 0.5_dp*(dist2 - r*r)/factor
                g12m1  = g12 * imass1
                g12m2  = g12 * imass2
                v12m1  = g12m1 / real(dt,dp)
                v12m2  = g12m2 / real(dt,dp)

                HGr_shake_force(ih,k,j,icel) = &
                     HGr_shake_force(ih,k,j,icel) + real(g12,wip)

                coord_dtmp(1,1)    = coord_dtmp(1,1)    - g12m1 * x12_old
                coord_dtmp(2,1)    = coord_dtmp(2,1)    - g12m1 * y12_old
                coord_dtmp(3,1)    = coord_dtmp(3,1)    - g12m1 * z12_old
                coord_dtmp(1,ih+1) = coord_dtmp(1,ih+1) + g12m2 * x12_old
                coord_dtmp(2,ih+1) = coord_dtmp(2,ih+1) + g12m2 * y12_old
                coord_dtmp(3,ih+1) = coord_dtmp(3,ih+1) + g12m2 * z12_old

                if (vel_update) then
                  vel_dtmp(1,1)    = vel_dtmp(1,1)    - v12m1 * x12_old
                  vel_dtmp(2,1)    = vel_dtmp(2,1)    - v12m1 * y12_old
                  vel_dtmp(3,1)    = vel_dtmp(3,1)    - v12m1 * z12_old
                  vel_dtmp(1,ih+1) = vel_dtmp(1,ih+1) + v12m2 * x12_old
                  vel_dtmp(2,ih+1) = vel_dtmp(2,ih+1) + v12m2 * y12_old
                  vel_dtmp(3,ih+1) = vel_dtmp(3,ih+1) + v12m2 * z12_old
                end if 

              end if
            end do

            if (shake_end) exit

          end do

          if (.not.shake_end) then
            write(MsgOut, '(A,8i10)') &
              'Compute_Shake> SHAKE algorithm failed to converge: indexes', &
              domain%id_l2g(iatm(1:j+1), icel)
            call error_msg('')
          end if

          coord(1:3,iatm(1:j+1),icel) = &
               real(coord_dtmp(1:3,1:j+1),wip)

          vel  (1:3,iatm(1:j+1),icel) = &
               real(vel_dtmp  (1:3,1:j+1),wip)

        end do
      end do
    end do


    ! compute constraint virial
    !   Note: virial =>  virial/dt**2 in compute_constraint
    !
    viri(1:3) = 0.0_dp

    do icel = id+1, ncell, nthread
      do j = 1, connect
        do k = 1, HGr_local(j,icel)
!ocl nosimd
          do ih = 1, j

            x12_old = HGr_bond_vector(1,ih,k,j,icel)
            y12_old = HGr_bond_vector(2,ih,k,j,icel)
            z12_old = HGr_bond_vector(3,ih,k,j,icel)
            g12 = real(HGr_shake_force(ih,k,j,icel),dp)
      
            fx  = g12 * x12_old
            fy  = g12 * y12_old
            fz  = g12 * z12_old
            viri(1) = viri(1) + x12_old * fx
            viri(2) = viri(2) + y12_old * fy
            viri(3) = viri(3) + z12_old * fz

          end do
        end do
      end do
    end do

    ! FEP: virials for atoms in single topology are scaled by lambbond
    if (domain%fep_use) then
      do icel = id+1, ncell, nthread
        do j = 1, connect
          do k = 1, HGr_local(j,icel)
            iatm(1:j+1) = HGr_bond_list(1:j+1,k,j,icel)
            do ih = 1, j
              if ((domain%fepgrp(iatm(ih+1),icel) == 1) .and. &
                  (domain%fepgrp(iatm(1),icel) == 1)) then
                lambbond = real(domain%lambbondA+domain%lambbondB,dp)
              else
                lambbond = 1.0_dp
              end if
              x12_old = HGr_bond_vector(1,ih,k,j,icel)
              y12_old = HGr_bond_vector(2,ih,k,j,icel)
              z12_old = HGr_bond_vector(3,ih,k,j,icel)
              g12 = real(HGr_shake_force(ih,k,j,icel),dp) * (lambbond - 1.0_dp)
              fx  = g12 * x12_old
              fy  = g12 * y12_old
              fz  = g12 * z12_old
              viri(1) = viri(1) + x12_old * fx
              viri(2) = viri(2) + y12_old * fy
              viri(3) = viri(3) + z12_old * fz
            end do
          end do
        end do
      end do
    end if

    virial_omp(1,1,id+1) = virial_omp(1,1,id+1) - viri(1)
    virial_omp(2,2,id+1) = virial_omp(2,2,id+1) - viri(2)
    virial_omp(3,3,id+1) = virial_omp(3,3,id+1) - viri(3)

    !$omp end parallel

!ocl nosimd
    do i = 1, nthread
      virial(1:3,1:3) = virial(1:3,1:3) + real(virial_omp(1:3,1:3,i),wip)
    end do

    return

  end subroutine compute_shake

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_rattle_fast_vv1
  !> @brief        SETTLE for three-site water
  !! @authors      JJ
  !! @param[in]    dt          : time step
  !! @param[in]    coord_old   : reference coordinates
  !! @param[in]    domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] coord       : coordinates
  !! @param[inout] vel         : velocities
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_rattle_fast_vv1(dt, coord_old, domain, &
                                     constraints, coord, vel)

    ! formal arguments
    real(wip),                intent(in)    :: dt
    real(wip),                intent(in)    :: coord_old(:,:,:)
    type(s_domain),   target, intent(in)    :: domain
    type(s_constraints),      intent(inout) :: constraints
    real(wip),                intent(inout) :: coord(:,:,:)
    real(wip),                intent(inout) :: vel(:,:,:)

    ! local variables
    real(dp)                  :: rHH, rOH, massO, massH
    real(dp)                  :: mo, mh, mohh, ra, rb, rc, rc2, tmp, tmp2
    real(dp)                  :: alpa, beta, gama, sqr_ab
    real(dp)                  :: axlng, aylng, azlng
    real(dp)                  :: t1(1:3), t2(1:3), t3(1:3)
    real(dp)                  :: cosphi, costhe, sinphi
    real(dp)                  :: sinthe, cospsi, sinpsi
    real(dp)                  :: aksxd(1:3), aksyd(1:3), akszd(1:3)
    real(dp)                  :: rab(1:3), rac(1:3), rabc(1:3)
    real(dp)                  :: a1(1:3), b1(1:3), c1(1:3)
    real(dp)                  :: a3(1:3), b3(1:3), c3(1:3)
    real(dp)                  :: b0d(1:3), c0d(1:3)
    real(dp)                  :: a1d(1:3), b1d(1:3), c1d(1:3)
    real(dp)                  :: ya2d, yb2d, yc2d
    real(dp)                  :: a3d(1:3), b3d(1:3), c3d(1:3)
    real(dp)                  :: xb2d, xb2d2, hh2, deltx, hhhh
    integer                   :: i, ix, iatom(3)
    integer                   :: id, omp_get_thread_num

    real(wip),        pointer :: mass(:,:)
    integer,          pointer :: nwater(:), water_list(:,:,:)
    integer,          pointer :: ncell



    mass       => domain%mass
    nwater     => domain%num_water
    water_list => domain%water_list
    ncell      => domain%num_cell_local

    rOH        =  constraints%water_rOH
    rHH        =  constraints%water_rHH
    massO      =  constraints%water_massO
    massH      =  constraints%water_massH

    mo         =  massO
    mh         =  massH
    mohh       =  massO + 2.0_dp * massH
    rc         =  rHH / 2.0_dp
    ra         =  2.0_dp * mh * sqrt(rOH * rOH - rc * rc)/mohh
    rb         =  sqrt(rOH * rOH - rc * rc) - ra
    rc2        =  rHH
    hhhh       =  rHH * rHH
    mo         =  mo/mohh
    mh         =  mh/mohh

    !$omp parallel default(shared)                                             &
    !$omp private(i, ix, iatom, rab, rac, rabc, a1, b1, c1, a3, b3, c3,        &
    !$omp         aksXd, aksYd, aksZd, axlng, aylng, azlng, t1, t2, t3,        &
    !$omp         b0d, c0d, a1d, b1d, c1d, a3d, b3d, c3d, ya2d, yb2d, yc2d,    &
    !$omp         xb2d, xb2d2, sinphi, cosphi, sinpsi, cospsi, sinthe, costhe, &
    !$omp         tmp, tmp2, hh2, deltx, alpa, beta, gama, sqr_ab, id)
    !
    ! note: reduction cannot be used for "virial" because of the rounding error,
    ! which is very large especialy in the case of langevin NPT.
    !

#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

!ocl nosimd
    do i = id+1, ncell, nthread
      do ix = 1, nwater(i)

      iatom(1:3) = water_list(1:3,ix,i)

      ! A1'
      rab(1:3) = coord_old(1:3,iatom(2),i) - coord_old(1:3,iatom(1),i)
      rac(1:3) = coord_old(1:3,iatom(3),i) - coord_old(1:3,iatom(1),i)

      rabc(1:3) = coord(1:3,iatom(1),i)*mo + &
                 (coord(1:3,iatom(2),i) + coord(1:3,iatom(3),i))*mh

      a1(1:3) = coord(1:3,iatom(1),i) - rabc(1:3)
      b1(1:3) = coord(1:3,iatom(2),i) - rabc(1:3)
      c1(1:3) = coord(1:3,iatom(3),i) - rabc(1:3)

      aksZd(1) = rab(2)*rac(3) - rab(3)*rac(2)
      aksZd(2) = rab(3)*rac(1) - rab(1)*rac(3)
      aksZd(3) = rab(1)*rac(2) - rab(2)*rac(1)
      aksXd(1) = a1(2)*aksZd(3) - a1(3)*aksZd(2)
      aksXd(2) = a1(3)*aksZd(1) - a1(1)*aksZd(3)
      aksXd(3) = a1(1)*aksZd(2) - a1(2)*aksZd(1)
      aksYd(1) = aksZd(2)*aksXd(3) - aksZd(3)*aksXd(2)
      aksYd(2) = aksZd(3)*aksXd(1) - aksZd(1)*aksXd(3)
      aksYd(3) = aksZd(1)*aksXd(2) - aksZd(2)*aksXd(1)

      axlng=1.0_dp/sqrt(aksXd(1)*aksXd(1)+aksXd(2)*aksXd(2)+aksXd(3)*aksXd(3))
      aylng=1.0_dp/sqrt(aksYd(1)*aksYd(1)+aksYd(2)*aksYd(2)+aksYd(3)*aksYd(3))
      azlng=1.0_dp/sqrt(aksZd(1)*aksZd(1)+aksZd(2)*aksZd(2)+aksZd(3)*aksZd(3))

      t1(1:3) = aksXd(1:3) * axlng
      t2(1:3) = aksYd(1:3) * aylng
      t3(1:3) = aksZd(1:3) * azlng

      b0d(1) = t1(1)*rab(1) + t1(2)*rab(2) + t1(3)*rab(3)
      b0d(2) = t2(1)*rab(1) + t2(2)*rab(2) + t2(3)*rab(3)
      c0d(1) = t1(1)*rac(1) + t1(2)*rac(2) + t1(3)*rac(3)
      c0d(2) = t2(1)*rac(1) + t2(2)*rac(2) + t2(3)*rac(3)
      a1d(3) = t3(1)*a1(1) + t3(2)*a1(2) + t3(3)*a1(3)
      b1d(1) = t1(1)*b1(1) + t1(2)*b1(2) + t1(3)*b1(3)
      b1d(2) = t2(1)*b1(1) + t2(2)*b1(2) + t2(3)*b1(3)
      b1d(3) = t3(1)*b1(1) + t3(2)*b1(2) + t3(3)*b1(3)
      c1d(1) = t1(1)*c1(1) + t1(2)*c1(2) + t1(3)*c1(3)
      c1d(2) = t2(1)*c1(1) + t2(2)*c1(2) + t2(3)*c1(3)
      c1d(3) = t3(1)*c1(1) + t3(2)*c1(2) + t3(3)*c1(3)

      sinphi = a1d(3) / ra
      tmp    = 1.0_dp - sinphi*sinphi
      if (tmp <= 0.0_dp) then
        cosphi = 0.0_dp
      else
        cosphi = sqrt(tmp)
      end if

      sinpsi = ( b1d(3) - c1d(3) ) / (rc2 * cosphi)
      tmp2   = 1.0_dp - sinpsi*sinpsi
      if (tmp2 <= 0.0_dp ) then
        cospsi = 0.0_dp
      else
        cospsi = sqrt(tmp2)
      end if

      ya2d  =   ra * cosphi
      xb2d  = - rc * cospsi
      yb2d  = - rb * cosphi - rc *sinpsi * sinphi
      yc2d  = - rb * cosphi + rc *sinpsi * sinphi
      xb2d2 = xb2d * xb2d
      hh2   = 4.0_dp * xb2d2 + (yb2d-yc2d) * (yb2d-yc2d) &
                           + (b1d(3)-c1d(3)) * (b1d(3)-c1d(3))
      deltx = 2.0_dp * xb2d + sqrt(4.0_dp * xb2d2 - hh2 + hhhh)
      xb2d  = xb2d - deltx * 0.5_dp


      ! alpha, beta, gamma
      alpa = xb2d * (b0d(1)-c0d(1)) + b0d(2) * yb2d + c0d(2) * yc2d
      beta = xb2d * (c0d(2)-b0d(2)) + b0d(1) * yb2d + c0d(1) * yc2d
      gama = b0d(1)*b1d(2) - b1d(1)*b0d(2) + c0d(1)*c1d(2) - c1d(1)*c0d(2)

      sqr_ab = alpa * alpa + beta * beta
      sinthe = (alpa*gama - beta * sqrt(sqr_ab - gama * gama)) / sqr_ab


      ! A3'
      costhe = dsqrt(1.0_dp - sinthe * sinthe)
      a3d(1) = - ya2d * sinthe
      a3d(2) =   ya2d * costhe
      a3d(3) = a1d(3)
      b3d(1) =   xb2d * costhe - yb2d * sinthe
      b3d(2) =   xb2d * sinthe + yb2d * costhe
      b3d(3) = b1d(3)
      c3d(1) = - xb2d * costhe - yc2d * sinthe
      c3d(2) = - xb2d * sinthe + yc2d * costhe
      c3d(3) = c1d(3)


      ! A3
      a3(1:3) = t1(1:3)*a3d(1) + t2(1:3)*a3d(2) + t3(1:3)*a3d(3)
      b3(1:3) = t1(1:3)*b3d(1) + t2(1:3)*b3d(2) + t3(1:3)*b3d(3)
      c3(1:3) = t1(1:3)*c3d(1) + t2(1:3)*c3d(2) + t3(1:3)*c3d(3)

      coord(1:3,iatom(1),i) = real(rabc(1:3) + a3(1:3),wip)
      coord(1:3,iatom(2),i) = real(rabc(1:3) + b3(1:3),wip)
      coord(1:3,iatom(3),i) = real(rabc(1:3) + c3(1:3),wip)

      vel(1:3,iatom(1),i) = vel(1:3,iatom(1),i) + real((a3(1:3)-a1(1:3))/dt,wip)
      vel(1:3,iatom(2),i) = vel(1:3,iatom(2),i) + real((b3(1:3)-b1(1:3))/dt,wip)
      vel(1:3,iatom(3),i) = vel(1:3,iatom(3),i) + real((c3(1:3)-c1(1:3))/dt,wip)

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_rattle_fast_vv1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_rattle_vv1
  !> @brief        RATTLE_VV1 for rigid bonds
  !! @authors      JJ
  !! @param[in]    dt          : time step
  !! @param[in]    coord_old   : reference coordinates
  !! @param[in]    domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] coord       : coordinates
  !! @param[inout] vel         : velocities
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_rattle_vv1(dt, coord_old, domain, constraints, coord, vel)

    ! formal arguments
    real(wip),                   intent(in)    :: dt
    real(wip),                   intent(in)    :: coord_old(:,:,:)
    type(s_domain),      target, intent(in)    :: domain
    type(s_constraints), target, intent(inout) :: constraints
    real(wip),                   intent(inout) :: coord(:,:,:)
    real(wip),                   intent(inout) :: vel(:,:,:)

    ! local variables
    real(dp)                     :: tolerance, imass1, imass2
    real(dp)                     :: x12, y12, z12
    real(dp)                     :: x12_old, y12_old, z12_old
    real(dp)                     :: dist, dist2, r
    real(dp)                     :: factor, g12, g12m1, g12m2, v12m1, v12m2
    real(dp)                     :: coord_dtmp(1:3,1:8), vel_dtmp(1:3,1:8)
    integer                      :: icel, i, j, k, ih, connect, id
    integer                      :: atm1, atm2, iatm(1:8)
    integer                      :: iteration, omp_get_thread_num
    logical                      :: shake_end

    real(dp),            pointer :: HGr_bond_vector(:,:,:,:,:)
    real(wip),           pointer :: HGr_bond_dist(:,:,:,:)
    real(wip),           pointer :: mass(:,:)
    integer,             pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)
    integer,             pointer :: ncell


    mass            => domain%mass
    ncell           => domain%num_cell_local

    HGr_local       => constraints%HGr_local
    HGr_bond_list   => constraints%HGr_bond_list
    HGr_bond_dist   => constraints%HGr_bond_dist
    HGr_bond_vector => constraints%HGr_bond_vector

    iteration       =  constraints%shake_iteration
    tolerance       =  real(constraints%shake_tolerance,dp)
    connect         =  constraints%connect

    !$omp parallel default(shared)                                             &
    !$omp private(icel, i, j, k, ih, atm1, atm2, shake_end, r, x12, y12, z12,  &
    !$omp         dist2, dist, imass1, imass2, x12_old, y12_old, z12_old,      &
    !$omp         factor, g12, g12m1, g12m2, v12m1, v12m2,          &
    !$omp         id, coord_dtmp, vel_dtmp, iatm)

#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    ! shake iteration
    !
    ! initialize force and store old bond vector
    !
    do icel = id+1, ncell, nthread
      do j = 1, connect
        do k = 1, HGr_local(j,icel)
          atm1 = HGr_bond_list(1,k,j,icel)
          do ih = 1, j
            atm2 = HGr_bond_list(ih+1,k,j,icel)
            HGr_bond_vector(1:3,ih,k,j,icel) = &
              real(coord_old(1:3,atm1,icel) - coord_old(1:3,atm2,icel),dp)
          end do
        end do
      end do
    end do

    shake_end = .true.

    do icel = id+1, ncell, nthread
      do j = 1, connect
        do k = 1, HGr_local(j,icel)

          iatm(1:j+1) = HGr_bond_list(1:j+1,k,j,icel)
          imass1 = 1.0_dp / real(mass(iatm(1),icel),dp)

          coord_dtmp(1:3,1:j+1) = &
            real(coord(1:3,iatm(1:j+1),icel),dp)

          vel_dtmp  (1:3,1:j+1) = &
            real(vel  (1:3,iatm(1:j+1),icel),dp)

          do i = 1, iteration

            shake_end = .true.
            do ih = 1, j

              imass2 = 1.0_dp / real(mass(iatm(ih+1),icel),dp)

              r = real(HGr_bond_dist(ih+1,k,j,icel),dp)

              x12 = coord_dtmp(1,1) - coord_dtmp(1,ih+1)
              y12 = coord_dtmp(2,1) - coord_dtmp(2,ih+1)
              z12 = coord_dtmp(3,1) - coord_dtmp(3,ih+1)
              dist2 = x12**2 + y12**2 + z12**2

              if (abs(dist2 - r*r) >= 2.0_dp*tolerance*r) then

                shake_end = .false.
                x12_old = HGr_bond_vector(1,ih,k,j,icel)
                y12_old = HGr_bond_vector(2,ih,k,j,icel)
                z12_old = HGr_bond_vector(3,ih,k,j,icel)

                factor = (x12*x12_old + y12*y12_old + z12*z12_old)* &
                         (imass1 + imass2)
                g12    = 0.5_dp*(dist2 - r*r)/factor
                g12m1  = g12 * imass1
                g12m2  = g12 * imass2
                v12m1  = g12m1 / real(dt,dp)
                v12m2  = g12m2 / real(dt,dp)

                coord_dtmp(1,1)    = coord_dtmp(1,1)    - g12m1 * x12_old
                coord_dtmp(2,1)    = coord_dtmp(2,1)    - g12m1 * y12_old
                coord_dtmp(3,1)    = coord_dtmp(3,1)    - g12m1 * z12_old
                coord_dtmp(1,ih+1) = coord_dtmp(1,ih+1) + g12m2 * x12_old
                coord_dtmp(2,ih+1) = coord_dtmp(2,ih+1) + g12m2 * y12_old
                coord_dtmp(3,ih+1) = coord_dtmp(3,ih+1) + g12m2 * z12_old

                vel_dtmp(1,1)    = vel_dtmp(1,1)    - v12m1 * x12_old
                vel_dtmp(2,1)    = vel_dtmp(2,1)    - v12m1 * y12_old
                vel_dtmp(3,1)    = vel_dtmp(3,1)    - v12m1 * z12_old
                vel_dtmp(1,ih+1) = vel_dtmp(1,ih+1) + v12m2 * x12_old
                vel_dtmp(2,ih+1) = vel_dtmp(2,ih+1) + v12m2 * y12_old
                vel_dtmp(3,ih+1) = vel_dtmp(3,ih+1) + v12m2 * z12_old

              end if
            end do

            if (shake_end) exit

          end do

          coord(1:3,iatm(1:j+1),icel) = &
               real(coord_dtmp(1:3,1:j+1),wip)

          vel  (1:3,iatm(1:j+1),icel) = &
               real(vel_dtmp  (1:3,1:j+1),wip)

          if (.not. shake_end) then
            write(MsgOut, '(A,8i10)') &
              'Compute_Rattle_VV1> SHAKE algorithm failed to converge: '//&
              'indexes',domain%id_l2g(iatm(1:j+1), icel)
            call error_msg('')
          end if

        end do
      end do
    end do

    !$omp end parallel

    return

  end subroutine compute_rattle_vv1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_rattle_fast_vv2
  !> @brief        SETTLE for three-site water
  !! @authors      JJ
  !! @param[in]    domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] coord       : coordinates
  !! @param[inout] vel         : velocities
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_rattle_fast_vv2(domain, constraints, coord, vel)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_constraints),      intent(inout) :: constraints
    real(wip),                intent(inout) :: coord(:,:,:)
    real(wip),                intent(inout) :: vel(:,:,:)

    ! local variables
    real(dp)                  :: rHH, rOH
    real(dp)                  :: massO, massH, massOH, massHH, mass2OH
    real(dp)                  :: rab(1:3), rbc(1:3), rca(1:3)
    real(dp)                  :: vab(1:3), vbc(1:3), vca(1:3)
    real(dp)                  :: unit_ab(1:3), unit_bc(1:3), unit_ca(1:3)
    real(dp)                  :: length_ab, length_bc, length_ca
    real(dp)                  :: vab0, vbc0, vca0, cosA, cosB, cosC
    real(dp)                  :: MaBC, MbCA, McAB, T_AB, T_BC, T_CA, Det
    integer                   :: i, ix, iatom(1:3), k
    integer                   :: id, omp_get_thread_num

    integer,          pointer :: nwater(:), water_list(:,:,:)
    integer,          pointer :: ncell


    nwater     => domain%num_water
    water_list => domain%water_list
    ncell      => domain%num_cell_local

    rOH        = constraints%water_rOH
    rHH        = constraints%water_rHH
    massO      = constraints%water_massO
    massH      = constraints%water_massH

    massOH     = massO + massH
    massHH     = massH + massH
    mass2OH    = 2.0_dp * massOH

    !$omp parallel default(shared)                                             &
    !$omp private(i, ix, k, rab, rbc, rca, vab, vbc, vca, length_ab, length_bc,&
    !$omp         length_ca, unit_ab, unit_bc, unit_ca, vab0, vbc0, vca0,      &
    !$omp         cosA, cosB, cosC, MaBC, MbCA, McAB, T_AB, T_BC, T_CA, Det,   &
    !$omp         iatom, id)
    !

#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell, nthread
!ocl nosimd
      do ix = 1, nwater(i)

      iatom(1:3) = water_list(1:3,ix,i)

      rab(1:3) = coord(1:3,iatom(2),i) - coord(1:3,iatom(1),i)
      rbc(1:3) = coord(1:3,iatom(3),i) - coord(1:3,iatom(2),i)
      rca(1:3) = coord(1:3,iatom(1),i)  - coord(1:3,iatom(3),i)

      vab(1:3) = vel(1:3,iatom(2),i) - vel(1:3,iatom(1),i)
      vbc(1:3) = vel(1:3,iatom(3),i) - vel(1:3,iatom(2),i)
      vca(1:3) = vel(1:3,iatom(1),i)  - vel(1:3,iatom(3),i)

      length_ab = sqrt(rab(1)**2 + rab(2)**2 + rab(3)**2)
      length_bc = sqrt(rbc(1)**2 + rbc(2)**2 + rbc(3)**2)
      length_ca = sqrt(rca(1)**2 + rca(2)**2 + rca(3)**2)

      unit_ab(1:3) = rab(1:3)/length_ab
      unit_bc(1:3) = rbc(1:3)/length_bc
      unit_ca(1:3) = rca(1:3)/length_ca

      vab0 = 0.0_dp
      vbc0 = 0.0_dp
      vca0 = 0.0_dp

      do k = 1, 3
        vab0 = vab0 + vab(k)*unit_ab(k)
        vbc0 = vbc0 + vbc(k)*unit_bc(k)
        vca0 = vca0 + vca(k)*unit_ca(k)
      end do

      cosA = 0.0_dp
      cosB = 0.0_dp
      cosC = 0.0_dp

      do k = 1, 3
        cosA = cosA - unit_ab(k)*unit_ca(k)
        cosB = cosB - unit_bc(k)*unit_ab(k)
        cosC = cosC - unit_ca(k)*unit_bc(k)
      end do

      MbCA = massH*cosC*cosA - massOH*cosB
      MaBC = massO*cosB*cosC - massHH*cosA
      McAB = massH*cosA*cosB - massOH*cosC

      T_AB = vab0*(mass2OH-massO*cosC**2) + vbc0*MbCA + vca0*MaBC
      T_AB = T_AB * massO
      T_BC = vbc0*(massOH**2-(massH*cosA)**2)
      T_BC = T_BC + vca0*massO*McAB + vab0*massO*MbCA
      T_CA = vca0*(mass2OH-massO*cosB**2) + vab0*MaBC + vbc0*McAB
      T_CA = T_CA * massO
      Det  = 2.0_dp*massOH**2 + 2.0_dp*massO*massH*cosA*cosB*cosC
      Det  = Det - 2.0_dp*massH**2*cosA**2 - massO*massOH*(cosB**2+cosC**2)
      Det  = Det / massH

      T_AB = T_AB / Det
      T_BC = T_BC / Det
      T_CA = T_CA / Det

      vel(1:3,iatom(1),i)  = vel(1:3,iatom(1),i) +  &
                       real((T_AB*unit_ab(1:3)-T_CA*unit_ca(1:3))/massO,wip)
      vel(1:3,iatom(2),i) = vel(1:3,iatom(2),i) + &
                       real((T_BC*unit_bc(1:3)-T_AB*unit_ab(1:3))/massH,wip)
      vel(1:3,iatom(3),i) = vel(1:3,iatom(3),i) + &
                       real((T_CA*unit_ca(1:3)-T_BC*unit_bc(1:3))/massH,wip)

      end do
    end do

    !$omp end parallel

    return

  end subroutine compute_rattle_fast_vv2
 
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_rattle_vv2
  !> @brief        RATTLE_VV2 for rigid bonds
  !! @authors      JJ
  !! @param[in]    domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] coord       : coordinates
  !! @param[inout] vel         : velocities
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_rattle_vv2(domain, constraints, coord, vel)

    ! formal arguments
    type(s_domain),      target, intent(in)    :: domain
    type(s_constraints), target, intent(inout) :: constraints
    real(wip),                   intent(inout) :: coord(:,:,:)
    real(wip),                   intent(inout) :: vel(:,:,:)

    ! local variables
    real(dp)                     :: tolerance, imass1, imass2
    real(dp)                     :: x12, y12, z12
    real(dp)                     :: vel_x12, vel_y12, vel_z12
    real(dp)                     :: r, dot_d12_v12
    real(dp)                     :: factor, g12, g12m1, g12m2
    real(dp)                     :: coord_dtmp(1:3,1:8), vel_dtmp(1:3,1:8)
    integer                      :: icel, i, j, k, ih, connect, id
    integer                      :: iteration, omp_get_thread_num
    integer                      :: iatm(1:8)
    logical                      :: rattle_end

    real(wip),           pointer :: HGr_bond_dist(:,:,:,:)
    real(wip),           pointer :: mass(:,:)
    integer,             pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)
    integer,             pointer :: ncell


    mass          => domain%mass
    ncell         => domain%num_cell_local

    HGr_local     => constraints%HGr_local
    HGr_bond_list => constraints%HGr_bond_list
    HGr_bond_dist => constraints%HGr_bond_dist

    iteration     =  constraints%shake_iteration
    tolerance     =  real(constraints%shake_tolerance,dp)
    connect       =  constraints%connect

    !$omp parallel default(shared)                                             &
    !$omp private(icel, i, j, k, ih, rattle_end, r, x12, y12, z12,             &
    !$omp         dot_d12_v12, imass1, imass2, vel_x12, vel_y12, vel_z12,      &
    !$omp         factor, g12, g12m1, g12m2, id, coord_dtmp, vel_dtmp, iatm)

#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    ! rattle iteration
    !
    do icel = id+1, ncell, nthread
!ocl nosimd
      do j = 1, connect
        do k = 1, HGr_local(j,icel)

          iatm(1:j+1) = HGr_bond_list(1:j+1,k,j,icel)
          imass1 = 1.0_dp / mass(iatm(1),icel)

          coord_dtmp(1:3,1:j+1) = &
            real(coord(1:3,iatm(1:j+1),icel),dp)

          vel_dtmp  (1:3,1:j+1) = &
            real(vel  (1:3,iatm(1:j+1),icel),dp)

          do i = 1, iteration

            rattle_end = .true.

            do ih = 1, j

              imass2 = 1.0_dp / dble(mass(iatm(ih+1),icel))

              r = real(HGr_bond_dist(ih+1,k,j,icel),dp)

              x12     = coord_dtmp(1,1) - coord_dtmp(1,ih+1)
              y12     = coord_dtmp(2,1) - coord_dtmp(2,ih+1)
              z12     = coord_dtmp(3,1) - coord_dtmp(3,ih+1)
              vel_x12 = vel_dtmp(1,1)   - vel_dtmp(1,ih+1)
              vel_y12 = vel_dtmp(2,1)   - vel_dtmp(2,ih+1)
              vel_z12 = vel_dtmp(3,1)   - vel_dtmp(3,ih+1)

              dot_d12_v12 = x12*vel_x12 + y12*vel_y12 + z12*vel_z12

              if (abs(dot_d12_v12) >= tolerance) then

                rattle_end = .false.

                g12 = dot_d12_v12 / ( (imass1+imass2)*r**2)
                g12m1  = g12 * imass1
                g12m2  = g12 * imass2

                vel_dtmp(1,1)    = vel_dtmp(1,1)    - g12m1 * x12
                vel_dtmp(2,1)    = vel_dtmp(2,1)    - g12m1 * y12
                vel_dtmp(3,1)    = vel_dtmp(3,1)    - g12m1 * z12
                vel_dtmp(1,ih+1) = vel_dtmp(1,ih+1) + g12m2 * x12
                vel_dtmp(2,ih+1) = vel_dtmp(2,ih+1) + g12m2 * y12
                vel_dtmp(3,ih+1) = vel_dtmp(3,ih+1) + g12m2 * z12

              end if
            end do

            if (rattle_end) exit

          end do

          coord(1:3,iatm(1:j+1),icel) = &
               real(coord_dtmp(1:3,1:j+1),wip)

          vel  (1:3,iatm(1:j+1),icel) = &
               real(vel_dtmp  (1:3,1:j+1),wip)

          if (.not. rattle_end) then
            write(MsgOut, '(A,8i10)') &
              'Compute_Rattle_VV2> SHAKE algorithm failed to converge: '//&
              'indexes',domain%id_l2g(iatm(1:j+1), icel)
            call error_msg('')
          end if

        end do
      end do
    end do

    !$omp end parallel

    return

  end subroutine compute_rattle_vv2

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    decide_dummy
  !> @brief        Decide Dummy atom posiion
  !! @authors      JJ
  !! @param[in]    domain      : domain information
  !! @param[in]    constraints : constraints information
  !! @param[inout] coord       : coordinates
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine decide_dummy(domain, constraints, coord)

    ! formal arguments
    type(s_domain),      target, intent(in)    :: domain
    type(s_constraints),         intent(in)    :: constraints
    real(wip),                   intent(inout) :: coord(:,:,:)

    integer,       pointer :: ncell, nwater(:), water_list(:,:,:)

    ! local variables
    integer                :: i, j, index(4)
    integer                :: id, omp_get_thread_num
    real(wip)              :: rOD, dist_rijk, rij(3), rjk(3)
    real(wip)              :: rijk(3), uijk(3)


    ncell          => domain%num_cell_local
    nwater         => domain%num_water
    water_list     => domain%water_list

    rOD        =  constraints%water_rOD

    !$omp parallel default(shared)                                     &
    !$omp          private(id, i, j, index, rij, rjk, rijk, dist_rijk, &
    !$omp                  uijk)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread
      do j = 1, nwater(i)

        index(1:4) = water_list(1:4,j,i)
        rij(1:3)   = coord(1:3,index(2),i) - coord(1:3,index(1),i)
        rjk(1:3)   = coord(1:3,index(3),i) - coord(1:3,index(2),i)
        rijk(1:3)  = rij(1:3) + 0.5_wip*rjk(1:3)
        dist_rijk  = sqrt(rijk(1)*rijk(1) + rijk(2)*rijk(2) + rijk(3)*rijk(3))
        uijk(1:3)  = rijk(1:3) / dist_rijk
        coord(1:3,index(4),i) = coord(1:3,index(1),i) &
                              + rOD*uijk(1:3)
      end do
    end do
    !$omp end parallel

    return

  end subroutine decide_dummy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    water_force_redistribution
  !> @brief        Dummy atom force is redistributed to other OH2 atoms
  !! @authors      JJ
  !! @param[in]    constraints : constraints information
  !! @param[in]    domain      : domain information
  !! @param[inout] force       : forces
  !! @param[inout] virial      : virial
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine water_force_redistribution(constraints, domain, force, virial)

    ! formal arguments
    type(s_constraints),         intent(in)    :: constraints
    type(s_domain),      target, intent(in)    :: domain
    real(wip),                   intent(inout) :: force(:,:,:)
    real(dp ),                   intent(inout) :: virial(:,:)

    integer,       pointer :: ncell, nwater(:), water_list(:,:,:)
    real(wip),     pointer :: coord(:,:,:)

    ! local variables
    integer                :: i, j, index(4), k
    integer                :: id, omp_get_thread_num
    real(wip)              :: factor, rOD, dist_rijk, dist_rid
    real(wip)              :: rid(3), rij(3), rjk(3), rijk(3)
    real(wip)              :: Fd(3), F1(3), rid_Fd
    real(dp )              :: virial_sub(3), virial_add(3)


    ncell          => domain%num_cell_local
    nwater         => domain%num_water
    water_list     => domain%water_list
    coord          => domain%coord

    rOD        =  constraints%water_rOD

    virial_sub(1:3) = 0.0_dp 
    virial_add(1:3) = 0.0_dp 

    !$omp parallel default(shared)                                          &
    !$omp          private(id, i, j, index, rid, rij, rjk, rijk, dist_rijk, &
    !$omp                  dist_rid, factor, Fd, rid_Fd, F1)                &
    !$omp          reduction(+:virial_sub, virial_add)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread

      do j = 1, nwater(i)
        index(1:4) = water_list(1:4,j,i)
        do k = 1, 4
          virial_sub(1:3) = virial_sub(1:3) &
                          + force(1:3,index(k),i)*coord(1:3,index(k),i)
        end do
        rid(1:3)   = coord(1:3,index(4),i) - coord(1:3,index(1),i)
        rij(1:3)   = coord(1:3,index(2),i) - coord(1:3,index(1),i)
        rjk(1:3)   = coord(1:3,index(3),i) - coord(1:3,index(2),i)
        rijk(1:3)  = rij(1:3) + 0.5_wip*rjk(1:3)
        dist_rijk  = sqrt(rijk(1)*rijk(1) + rijk(2)*rijk(2) + rijk(3)*rijk(3))
        dist_rid   = rid(1)*rid(1) + rid(2)*rid(2) + rid(3)*rid(3)

        factor     = rOD / dist_rijk
        Fd(1:3)    = force(1:3,index(4),i)
        rid_Fd     = rid(1)*Fd(1) + rid(2)*Fd(2) + rid(3)*Fd(3)
        F1(1:3)    = rid(1:3) * (rid_Fd/dist_rid)
        F1(1:3)    = Fd(1:3) - F1(1:3)
        force(1:3,index(1),i) = force(1:3,index(1),i) &
                              + Fd(1:3) - factor*F1(1:3)
        force(1:3,index(2),i) = force(1:3,index(2),i) &
                              + 0.5_wip*factor*F1(1:3)
        force(1:3,index(3),i) = force(1:3,index(3),i) &
                              + 0.5_wip*factor*F1(1:3)
        force(1:3,index(4),i) = 0.0_wip
        do k = 1, 3
          virial_add(1:3) = virial_add(1:3) &
                          + force(1:3,index(k),i)*coord(1:3,index(k),i)
        end do
     end do

   end do
   !$omp end parallel

   virial(1,1) = virial(1,1) - virial_sub(1) + virial_add(1)
   virial(2,2) = virial(2,2) - virial_sub(2) + virial_add(2)
   virial(3,3) = virial(3,3) - virial_sub(3) + virial_add(3)

   return

 end subroutine water_force_redistribution

end module sp_constraints_mod
