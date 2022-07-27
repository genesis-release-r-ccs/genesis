!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_constraints_mod
!> @brief   constraints module
!! @authors Takaharu Mori (TM), Jaewoon Jung (JJ), Chigusa Kobayashi (CK), 
!!          Norio Takase (NT), Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_constraints_mod

  use at_dynamics_str_mod
  use at_dynvars_str_mod
  use at_constraints_str_mod
  use at_enefunc_str_mod
  use at_boundary_str_mod
  use molecules_mod
  use molecules_str_mod
  use fileio_control_mod
  use timers_mod
  use string_mod
  use select_mod
  use select_atoms_mod
  use select_atoms_str_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_cons_info
    logical                  :: rigid_bond      = .false.
    logical                  :: fast_bond       = .false.
    logical                  :: fast_water      = .true.
    integer                  :: hydrogen_type   = ConstraintAtomName
    integer                  :: shake_iteration = 500
    real(wp)                 :: shake_tolerance = 1.0e-10_wp
    integer                  :: lincs_iteration = 1
    integer                  :: lincs_order     = 4
    character(5)             :: water_model     = 'TIP3'
    real(wp)                 :: hydrogen_mass_upper_bound = 2.1_wp
    real(wp)                 :: water_rHH =  0.0_wp
    real(wp)                 :: water_rOH =  0.0_wp
    character(MaxLine)       :: noshake_index  = ''
    character(MaxLine)       :: fixatm_select_index = ''
  end type s_cons_info

  ! parameters
  integer, public, parameter :: ConstraintModeLEAP  = 1
  integer, public, parameter :: ConstraintModeVVER1 = 2
  integer, public, parameter :: ConstraintModeVVER2 = 3
  integer, public, parameter :: ConstraintModeBD    = 4


  ! subroutines
  public  :: show_ctrl_constraints
  public  :: read_ctrl_constraints
  public  :: setup_constraints
  public  :: compute_constraints
  private :: setup_shake
  private :: setup_shake_qmmm
  private :: setup_lincs
  private :: setup_settle
  private :: setup_fixatm
  private :: remove_fixatm
  private :: remove_fixatm_bond
  private :: remove_fixatm_angl
  private :: remove_fixatm_urey
  private :: remove_fixatm_dihe
  private :: remove_fixatm_rb_dihe
  private :: remove_fixatm_impr
  private :: remove_fixatm_cmap
  private :: exclude_light_atom
  private :: compute_shake
  private :: compute_settle
  private :: compute_lincs
  private :: compute_rattle
  private :: compute_settle_vv2
  public  :: clear_fixatm_component

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_constraints
  !> @brief        show CONSTRAINTS section usage
  !! @authors      NT, TM
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "md", "min", "remd", "rpath"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_constraints(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode


    if (show_all) then

      select case (run_mode)

      case ('md', 'remd', 'rpath', 'bd')

        write(MsgOut,'(A)') '[CONSTRAINTS]'
        write(MsgOut,'(A)') 'rigid_bond    = NO        # constraints all bonds involving hydrogen'
        write(MsgOut,'(A)') 'fast_water    = YES       # use SETTLE algorithm'
        write(MsgOut,'(A)') '# shake_iteration = 500     # max number of SHAKE/RATTLE iterations'
        write(MsgOut,'(A)') '# shake_tolerance = 1.0e-10 # SHAKE/RATTLE tolerance (Ang)'
        write(MsgOut,'(A)') '# water_model     = TIP3    # water model'
        write(MsgOut,'(A)') '# hydrogen_mass_upper_bound    = 2.1    # water model'
        write(MsgOut,'(A)') '# noshake_index = 1         # Hydrogen atoms without SHAKE'
        write(MsgOut,'(A)') '# fixatm_select_index = 1 2 # fixed atoms'
        write(MsgOut,'(A)') ' '


      end select

    else

      select case (run_mode)

      case ('md', 'remd', 'rpath', 'bd')

        write(MsgOut,'(A)') '[CONSTRAINTS]'
        write(MsgOut,'(A)') 'rigid_bond    = NO        # constraints all bonds involving hydrogen'
        write(MsgOut,'(A)') ' '


      end select

    end if

    return

  end subroutine show_ctrl_constraints
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_constraints
  !> @brief        read CONSTRAINTS section in the control file
  !! @authors      TM
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
    call read_ctrlfile_logical(handle, Section, 'fast_bond',       &
                               cons_info%fast_bond)
    call read_ctrlfile_integer(handle, Section, 'shake_iteration', &
                               cons_info%shake_iteration)
    call read_ctrlfile_real   (handle, Section, 'shake_tolerance', &
                               cons_info%shake_tolerance)
    call read_ctrlfile_string (handle, Section, 'water_model',     &
                               cons_info%water_model)
    call read_ctrlfile_integer(handle, Section, 'lincs_iteration', &
                               cons_info%lincs_iteration)
    call read_ctrlfile_integer(handle, Section, 'lincs_order',     &
                               cons_info%lincs_order)
    call read_ctrlfile_type(handle, Section, 'hydrogen_type',      &
                               cons_info%hydrogen_type, ConstraintAtomType)
    call read_ctrlfile_real   (handle, Section, 'hydrogen_mass_upper_bound', &
                               cons_info%hydrogen_mass_upper_bound)
    call read_ctrlfile_real   (handle, Section, 'water_rHH', &
                               cons_info%water_rHH)
    call read_ctrlfile_real   (handle, Section, 'water_rOH', &
                               cons_info%water_rOH)
    call read_ctrlfile_string (handle, Section, 'noshake_index', &
                               cons_info%noshake_index)
    call read_ctrlfile_string (handle, Section, 'fixatm_select_index',      &
                               cons_info%fixatm_select_index)

    call end_ctrlfile_section(handle)


    ! write parameters to MsgOut
    !
    if (main_rank) then

      write(MsgOut,'(A)') 'Read_Ctrl_Constraints> Parameters for Constraints'

      if (cons_info%rigid_bond) then
        write(MsgOut,'(A30)')                                      &
              '  rigid_bond      =        yes'
        write(MsgOut,'(A20,I10,A20,E10.3)')                        &
              '  shake_iteration = ', cons_info%shake_iteration,   &
              '  shake_tolerance = ', cons_info%shake_tolerance

        if (cons_info%fast_bond) then
          write(MsgOut,'(A30)')                                    &
                '  fast_bond       =        yes'
          write(MsgOut,'(A20,I10,A20,I10)')                        &
                '  lincs_iteration = ', cons_info%lincs_iteration, &
                '  lincs_order     = ', cons_info%lincs_order
        end if

        if (cons_info%fast_water) then
          write(MsgOut,'(A30,A20,A10)')                            &
                '  fast_water      =        yes',                  &
                '  water_model     = ', trim(cons_info%water_model)
          if (cons_info%water_rHH > 0.0_wp)                        &
            write(MsgOut,'(A20,F10.3)')                            &
                '  water_rHH       = ', cons_info%water_rHH
          if (cons_info%water_rOH > 0.0_wp)                        &
            write(MsgOut,'(A20,F10.3)')                            &
                '  water_rOH       = ', cons_info%water_rOH
        else
          write(MsgOut,'(A30)')                                    &
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

        if (.not. trim(cons_info%noshake_index) .eq. '') then
          write(MsgOut,'(a,a)') &
                '  noshake_index   = ', trim(cons_info%noshake_index)
        else
          write(MsgOut,'(a)') &
                '  noshake_index   =       none'
        end if

      else
        write(MsgOut,'(A30)') '  rigid_bond      =         no'
      end if

      if (len(trim(cons_info%fixatm_select_index)) /= 0) then
        write(MsgOut,'(A,A)')             &
            '  fixatm_select_index        = ', trim(cons_info%fixatm_select_index)
      else
        write(MsgOut,'(A)')               &
            '  fixatm_select_index        =       none'
      end if

      write(MsgOut,'(A)') ' '

    end if


    ! error check
    !
    if (main_rank) then

      if (cons_info%fast_bond .and. &
          .not. cons_info%rigid_bond) then
        call error_msg( &
          'Read_Ctrl_Constraints> rigid_bond must be YES, if fast_bond is YES')
      end if

    end if

    return    

  end subroutine read_ctrl_constraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_constraints
  !> @brief        setup constraints information
  !! @authors      TM, KY
  !! @param[in]    cons_info :CONSTRAINTS section control parameters information
  !! @param[in]    sel_info  : SELECTION section control parameters information
  !! @param[in]    dynamics    : dynamics information
  !! @param[inout] boundary    : information of boundary condition
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_constraints(cons_info, sel_info, dynamics, boundary, &
                               molecule, enefunc, constraints)

    ! formal arguments
    type(s_cons_info),       intent(in)    :: cons_info
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_boundary),        intent(inout) :: boundary
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc), target, intent(inout) :: enefunc
    type(s_constraints),     intent(inout) :: constraints

    ! local arguments
    type(s_qmmm), pointer :: qmmm


    qmmm => enefunc%qmmm

    ! initialize
    !
    call init_constraints(constraints)

    constraints%shake_iteration = cons_info%shake_iteration
    constraints%shake_tolerance = cons_info%shake_tolerance
    constraints%lincs_iteration = cons_info%lincs_iteration
    constraints%hydrogen_type   = cons_info%hydrogen_type
    constraints%lincs_order     = cons_info%lincs_order
    constraints%water_model     = cons_info%water_model
    constraints%water_rHH       = cons_info%water_rHH
    constraints%water_rOH       = cons_info%water_rOH
    
    ! setup fixed atoms
    !
    call setup_fixatm(molecule, boundary, dynamics, enefunc, &
                      sel_info, cons_info, constraints)


    ! decide hydrogen atom
    !
    call check_light_atom_name(cons_info%hydrogen_mass_upper_bound, &
                               molecule)
    if (.not. trim(cons_info%noshake_index) .eq. '') &
      call exclude_light_atom(cons_info%noshake_index, sel_info, molecule)

    if (cons_info%rigid_bond) then

      if (constraints%hydrogen_type .eq. ConstraintAtomName  &
         .and. molecule%special_hydrogen) &
      call error_msg('Setup_Constraints> Non ordinary hydrogen name is not'//&
                         ' allowed. If you want, use hydrogen_type option.')

      ! setup SETTLE
      !
      if (cons_info%fast_water) then

        call setup_settle(molecule, boundary, enefunc, constraints)

      end if

      ! setup SHAKE and RATTLE
      !
      call setup_shake(molecule, enefunc, constraints)
      if (qmmm%do_qmmm .and. .not. dynamics%esp_mm) then
        call setup_shake_qmmm(cons_info, molecule, enefunc, constraints)
      end if

      ! setup LINCS
      !
      if (cons_info%fast_bond) then

        if (dynamics%integrator == IntegratorVVER) then
          call error_msg('Setup_Constraints> fast_bond = YES with integrator'//&
                         ' = VVER is not supported')
        end if

        call setup_lincs(molecule, dynamics, constraints)

      end if

      ! setup constraints
      !
      if (constraints%num_water > 0 .or. constraints%num_bonds > 0) then
        constraints%rigid_bond = .true.
        if (cons_info%fast_water .and. constraints%num_water > 0) then
          constraints%fast_water = .true.
        else
          constraints%fast_water = .false.
        end if
        if (cons_info%fast_bond .and. constraints%num_bonds > 0) then
          constraints%fast_bond  = .true.
        else
          constraints%fast_bond  = .false.
        end if
      else
        constraints%rigid_bond = .false.
        constraints%fast_water = .false.
        constraints%fast_bond  = .false.
      end if

      ! update number of degrees of freedom
      !
      if (constraints%num_water > 0) then
        call update_num_deg_freedom('After setup of SETTLE',       &
                                    -3*constraints%num_water,      &
                                    molecule%num_deg_freedom)
      end if

      if (constraints%num_bonds > 0) then
        call update_num_deg_freedom('After setup of SHAKE/RATTLE', &
                                    -constraints%num_bonds,        &
                                    molecule%num_deg_freedom)
      end if

    end if
   
    return

  end subroutine setup_constraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_shake
  !> @brief        setup rigid bonds for SHAKE and RATTLE
  !! @authors      TM, KY
  !! @param[in]    molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_shake(molecule, enefunc, constraints)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_constraints),     intent(inout) :: constraints

    ! local variables
    integer                  :: i, nbond, nnsbond, nsbond
    integer                  :: ns_found, s_found
    integer                  :: istart, iend
    integer                  :: alloc_stat, dealloc_stat
    character(6)             :: ci1, ci2
    logical                  :: mi1, mi2
    logical                  :: cl1, cl2

    real(wp),  allocatable   :: temp_fc(:), temp_r0(:)
    integer,   allocatable   :: temp_list(:,:)


    nbond = enefunc%num_bonds

    ! count the number of bonds including hydrogen in enefunc%bond_list
    !
    s_found = 0
    do i = 1, nbond
      mi1 = molecule%light_atom_mass(enefunc%bond_list(1,i))
      mi2 = molecule%light_atom_mass(enefunc%bond_list(2,i))
      cl1 = molecule%light_atom_name(enefunc%bond_list(1,i))
      cl2 = molecule%light_atom_name(enefunc%bond_list(2,i))
      if (constraints%hydrogen_type == ConstraintAtomMass) then
        cl1 = mi1 
        cl2 = mi2 
      else if (constraints%hydrogen_type == ConstraintAtomBoth) then
        cl1 = (cl1 .or. mi1) 
        cl2 = (cl2 .or. mi2) 
      end if

      if (cl1 .or. cl2 .or. &
           (enefunc%forcefield == ForcefieldKBGO .or. &
            enefunc%forcefield == ForcefieldRESIDCG )) then
        s_found = s_found + 1
      end if
    end do

    nsbond  = s_found
    nnsbond = nbond - nsbond

    if (nsbond == 0) then
      if (main_rank) then
        write(MsgOut,'(A)') 'Setup_Shake> bond for SHAKE was not found'
        write(MsgOut,'(A)') ' '
      end if
      constraints%num_bonds  = 0
      return
    end if

  
    ! copy information from enefunc to temp array
    !
    allocate(temp_fc(nbond), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_r0(nbond), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_list(2,nbond), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    temp_fc(1:nbond)       = enefunc%bond_force_const(1:nbond)
    temp_r0(1:nbond)       = enefunc%bond_dist_min(1:nbond)
    temp_list(1:2,1:nbond) = enefunc%bond_list(1:2,1:nbond)


    ! re-allocate EneFuncBond and allocate ConstraintsShake
    !
    call dealloc_enefunc(enefunc, EneFuncBond)

    enefunc%num_bonds     = nnsbond
    constraints%num_bonds = nsbond

    call alloc_enefunc(enefunc, EneFuncBond, nnsbond)
    call alloc_constraints(constraints, ConstraintsShake, nsbond) 


    ! modify enefunc%bond and setup constraints%bond
    !
    s_found  = 0
    ns_found = 0

    do i = 1, nbond
      mi1 = molecule%light_atom_mass(temp_list(1,i))
      mi2 = molecule%light_atom_mass(temp_list(2,i))
      cl1 = molecule%light_atom_name(temp_list(1,i))
      cl2 = molecule%light_atom_name(temp_list(2,i))
      if (constraints%hydrogen_type == ConstraintAtomMass) then
        cl1 = mi1 
        cl2 = mi2 
      else if (constraints%hydrogen_type == ConstraintAtomBoth) then
        cl1 = (cl1 .or. mi1) 
        cl2 = (cl2 .or. mi2) 
      end if

      if (cl1 .or. cl2 .or. &
           enefunc%forcefield == ForcefieldKBGO .or.  &
           enefunc%forcefield == ForcefieldRESIDCG ) then
        s_found = s_found + 1
        constraints%bond_dist(s_found)     = temp_r0(i)
        constraints%bond_list(1:2,s_found) = temp_list(1:2,i)
      else
        ns_found = ns_found + 1
        enefunc%bond_force_const(ns_found) = temp_fc(i)
        enefunc%bond_dist_min(ns_found)    = temp_r0(i)
        enefunc%bond_list(1:2,ns_found)    = temp_list(1:2,i)
      end if

    end do
    
    ! re-setup loop index
    !
    call get_loop_index(ns_found, istart, iend)
    enefunc%istart_bond = istart
    enefunc%iend_bond   = iend


    ! deallocate temp array
    !
    deallocate(temp_fc, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_r0, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_list, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc


    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') &
           'Setup_Shake> Setup constrains for SHAKE and RATTLE'
      write(MsgOut,'(A20,I10,A20,I10)')               &  
           '  num_unconsbonds = ', enefunc%num_bonds, &
           '  num_rigid_bonds = ', constraints%num_bonds
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine setup_shake

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_shake_qmmm
  !> @brief        add rigid bonds for SHAKE and RATTLE in QM region
  !! @authors      KY
  !! @param[in]    cons_info: CONSTRAINTS section control parameters information
  !! @param[in]    molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_shake_qmmm(cons_info, molecule, enefunc, constraints)

    ! formal arguments
    type(s_cons_info),       intent(in)    :: cons_info
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_constraints),     intent(inout) :: constraints

    ! local variables
    integer                  :: i, j, nbond, nsbond_qm, nsbond_mm
    integer                  :: s_found
    integer                  :: istart, iend
    integer                  :: alloc_stat, dealloc_stat
    character(6)             :: ci1, ci2
    logical                  :: mi1, mi2
    logical                  :: cl1, cl2
    integer,   allocatable   :: qm_bond_list(:,:)
    real(wp)                 :: r1(3), r2(3), r12

    real(wp),  allocatable   :: temp_bond_dist(:)
    integer,   allocatable   :: temp_bond_list(:,:)

    character(5)             :: water_model


    nbond = enefunc%qmmm%qm_nbonds
    water_model = constraints%water_model

    allocate(qm_bond_list(2,nbond), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    qm_bond_list = molecule%bond_list(:,(molecule%num_bonds - nbond + 1):molecule%num_bonds)

    ! count the number of bonds including hydrogen in molecule%bond_list
    !
    s_found = 0
    do i = 1, nbond

      if (cons_info%fast_water .and. &
          trim(molecule%residue_name(qm_bond_list(1,i))) == trim(water_model)) cycle

      if (constraints%fixatm(qm_bond_list(1,i)) .and. &
          constraints%fixatm(qm_bond_list(2,i))) cycle

      mi1 = molecule%light_atom_mass(qm_bond_list(1,i))
      mi2 = molecule%light_atom_mass(qm_bond_list(2,i))
      cl1 = molecule%light_atom_name(qm_bond_list(1,i))
      cl2 = molecule%light_atom_name(qm_bond_list(2,i))
      if (constraints%hydrogen_type == ConstraintAtomMass) then
        cl1 = mi1 
        cl2 = mi2 
      else if (constraints%hydrogen_type == ConstraintAtomBoth) then
        cl1 = (cl1 .or. mi1) 
        cl2 = (cl2 .or. mi2) 
      end if

      if (cl1 .or. cl2 .or. &
           enefunc%forcefield == ForcefieldKBGO .or.  &
           enefunc%forcefield == ForcefieldRESIDCG ) then
        s_found = s_found + 1
      end if
    end do

    nsbond_qm  = s_found
    if (nsbond_qm == 0) then
      if (main_rank) then
        write(MsgOut,'(A)') '  QM bond for SHAKE was not found'
        write(MsgOut,'(A)') ' '
      end if
      return
    end if

  
    ! copy information of current constraints
    !
    nsbond_mm = constraints%num_bonds
    if (nsbond_mm /= 0) then
      allocate(temp_bond_dist(nsbond_mm), stat=alloc_stat)
      if (alloc_stat /= 0) call error_msg_alloc
      allocate(temp_bond_list(2,nsbond_mm), stat=alloc_stat)
      if (alloc_stat /= 0) call error_msg_alloc

      temp_bond_dist = constraints%bond_dist
      temp_bond_list = constraints%bond_list
    end if

    ! re-allocate
    call dealloc_constraints(constraints, ConstraintsShake)
    constraints%num_bonds = nsbond_mm + nsbond_qm
    call alloc_constraints(constraints, ConstraintsShake, constraints%num_bonds) 

    if (nsbond_mm /= 0) then
      constraints%bond_dist(1:nsbond_mm)   = temp_bond_dist(1:nsbond_mm)
      constraints%bond_list(:,1:nsbond_mm) = temp_bond_list(:,1:nsbond_mm)

      deallocate(temp_bond_dist, stat=dealloc_stat)
      if (dealloc_stat /= 0) call error_msg_dealloc
      deallocate(temp_bond_list, stat=dealloc_stat)
      if (dealloc_stat /= 0) call error_msg_dealloc
    end if

    ! setup constraints%bond
    !
    s_found  = 0

    do i = 1, nbond

      if (cons_info%fast_water .and. &
          trim(molecule%residue_name(qm_bond_list(1,i))) == trim(water_model)) cycle

      if (constraints%fixatm(qm_bond_list(1,i)) .and. &
          constraints%fixatm(qm_bond_list(2,i))) cycle

      mi1 = molecule%light_atom_mass(qm_bond_list(1,i))
      mi2 = molecule%light_atom_mass(qm_bond_list(2,i))
      cl1 = molecule%light_atom_name(qm_bond_list(1,i))
      cl2 = molecule%light_atom_name(qm_bond_list(2,i))
      if (constraints%hydrogen_type == ConstraintAtomMass) then
        cl1 = mi1 
        cl2 = mi2 
      else if (constraints%hydrogen_type == ConstraintAtomBoth) then
        cl1 = (cl1 .or. mi1) 
        cl2 = (cl2 .or. mi2) 
      endif

      if (cl1 .or. cl2 .or. &
           enefunc%forcefield == ForcefieldKBGO .or.  &
           enefunc%forcefield == ForcefieldRESIDCG ) then

        s_found = s_found + 1
        constraints%bond_list(1:2,nsbond_mm+s_found) = qm_bond_list(1:2,i)   !YA
        r1=molecule%atom_coord(:,qm_bond_list(1,i))
        r2=molecule%atom_coord(:,qm_bond_list(2,i))
        r12 = 0.0_wp
        do j = 1, 3
           r12 = r12 + (r1(j) - r2(j)) * (r1(j) - r2(j))
        end do
        constraints%bond_dist(nsbond_mm+s_found) = sqrt(r12)   !YA

      end if

    end do
    
    ! deallocate temp array
    !
    deallocate(qm_bond_list, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc

    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A23,I10,A20,I10)')               &  
           '  num_qm_rigid_bonds = ', nsbond_qm
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine setup_shake_qmmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_lincs
  !> @brief        setup fast shake (LINCS)
  !! @authors      TM
  !! @param[in]    molecule    : molecule information
  !! @param[in]    dynamics    : dynamics information
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_lincs(molecule, dynamics, constraints)

    ! formal arguments
    type(s_molecule),            intent(in)    :: molecule
    type(s_dynamics),            intent(in)    :: dynamics
    type(s_constraints), target, intent(inout) :: constraints

    ! local variables
    real(wp)                     :: imass1, imass2
    integer                      :: i, j, nsbond
    integer                      :: ncc, nncc(8)
    integer                      :: cmax, con, an1, an2, comm
    logical                      :: friction

    integer,             pointer :: list(:,:)


    ! use pointer
    !
    list   => constraints%bond_list
    nsbond =  constraints%num_bonds

    if (nsbond == 0) then
      if (main_rank) then
        write(MsgOut,'(A)') 'Setup_Lincs> bond for LINCS was not found'
        write(MsgOut,'(A)') ' '
      end if
      return
    end if

    friction = ((dynamics%integrator == IntegratorBDEM) .or. &
                (dynamics%integrator == IntegratorBD2N) .or. &
                (dynamics%integrator == IntegratorSDMP))

    ! define cmax
    !
    cmax = 0
    do i = 1, nsbond
      an1 = list(1,i)
      an2 = list(2,i)

      ncc = 0
      do j = 1, nsbond
        if (j /= i) then
          if (an1 == list(1,j) .or. an1 == list(2,j) .or.  &
              an2 == list(1,j) .or. an2 == list(2,j)) then
            ncc = ncc + 1
          end if
        end if
      end do
      if (ncc > cmax) cmax = ncc
    end do


    ! allocation for ncc[K], Sdiag[K], rhs[2,K], sol[K]
    !                con[cmax,K], coef[cmax,K], and A[cmax,K]
    !
    call alloc_constraints(constraints, ConstraintsLincs, nsbond, cmax)


    ! setup Sdiag[K], ncc[K], and con[cmax,K]
    !
    nncc(1:8) = 0
    do i = 1, nsbond
      an1 = list(1,i)
      an2 = list(2,i)

      ncc = 0
      do j = 1, nsbond
        if (j /= i) then
          if (an1 == list(1,j) .or. an1 == list(2,j) .or.  &
              an2 == list(1,j) .or. an2 == list(2,j)) then
            ncc = ncc + 1
            constraints%connected_cons_idx(ncc,i) = j
          end if
        end if
      end do

      constraints%num_connected_cons(i) = ncc
      nncc(ncc+1) = nncc(ncc+1) + 1

      if (friction) then
        imass1 = molecule%inv_stokes_radius(an1)
        imass2 = molecule%inv_stokes_radius(an2)
      else
        imass1 = molecule%inv_mass(an1)
        imass2 = molecule%inv_mass(an2)
      end if
      constraints%s_diagonal(i) = 1.0_wp / sqrt(imass1 + imass2)

    end do


    ! setup coef[cmax,K]
    !
    do i = 1, nsbond

      ! define comm (the atom coupling constraints i and con)
      !
      an1 = list(1,i)
      an2 = list(2,i)
      ncc = constraints%num_connected_cons(i)

      if (ncc == 0) then
        comm = an1
      else
        do j = 1, ncc
          con = constraints%connected_cons_idx(j,i)
          if (an1 == list(1,con) .or. an1 == list(2,con)) then
            comm = an1
          else if (an2 == list(1,con) .or. an2 == list(2,con)) then
            comm = an2
          end if
        end do
      end if

      do j = 1, ncc
        con = constraints%connected_cons_idx(j,i)
        if (friction) then
          constraints%lincs_coef(j,i) = constraints%s_diagonal(i)  &
                                     * constraints%s_diagonal(con) &
                                     * molecule%inv_stokes_radius(comm)
        else
          constraints%lincs_coef(j,i) = constraints%s_diagonal(i)  &
                                     * constraints%s_diagonal(con) &
                                     * molecule%inv_mass(comm)
        end if
        if (an1 == list(1,con) .or. an2 == list(2,con)) then
          constraints%lincs_coef(j,i) = -constraints%lincs_coef(j,i)
        end if
      end do

    end do


    ! LINCS
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Setup_Lincs> Re-setup constraints for LINCS'
      write(MsgOut,'(A20,I10,A20,I10)')       &
           '  num_2-body (CH) = ', nncc(1),   &
           '  num_3-body (CH2)= ', nncc(2)/2
      write(MsgOut,'(A20,I10,A20,I10)')       &
           '  num_4-body (CH3)= ', nncc(3)/3, &
           '  num_>4-body     = ', sum(nncc(4:8))
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine setup_lincs

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_settle
  !> @brief        setup water molecule list for SETTLE
  !! @authors      TM, CK, KY
  !! @param[in]    molecule    : molecule information
  !! @param[in]    boundary    : information of boundary condition
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] constraints : constraints information
  !! @note         H-H bond list is constructed based on the angle H-O-H,
  !!               because namd psffile does not contain H-H bonds.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_settle(molecule, boundary, enefunc, constraints)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_boundary),        intent(in)    :: boundary
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_constraints),     intent(inout) :: constraints

    ! local variables
    integer                  :: i, i1, i2, i3, natom, nwater, nbond, nangl
    integer                  :: istart, iend
    integer                  :: nthread, omp_get_num_threads, alloc_stat
    character(10)            :: water_model
    character(4)             :: an1, an2, an3

    logical,     allocatable :: lwater(:)
    integer,     allocatable :: temp_list(:,:)
    real(wp),    allocatable :: temp_fc(:), temp_min(:)


    water_model = constraints%water_model
    natom       = molecule%num_atoms


    ! count water molecule
    !
    nwater = 0

    do i = 1, natom
      if (molecule%residue_name(i) == water_model .and.  &
          .not. constraints%fixatm(i)) &
        nwater = nwater + 1
    end do

    if (mod(nwater, 3) /= 0) &
      call error_msg('Setup_Settle> # of water is incorrect.')

    nwater = nwater / 3
    constraints%num_water = nwater

    if (constraints%num_water == 0) then
      if (main_rank) then
        write(MsgOut,'(A)') 'Setup_Settle> water for SETTLE was not found'
        write(MsgOut,'(A)') ' '
      end if
      return
    end if


    ! setup constraints water list
    !
    call alloc_constraints(constraints, ConstraintsSettle, nwater)

    allocate(lwater(natom), stat=alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc
    lwater(1:natom) = .false.

    i = 1
    nwater = 0
    do while(.true.)

      if (i > natom) &
        exit

      if (molecule%residue_name(i) == water_model .and. &
          .not. constraints%fixatm(i)) then

        an1 = molecule%atom_name(i)
        an2 = molecule%atom_name(i+1)
        an3 = molecule%atom_name(i+2)

        if (an1(1:1) .eq. 'H' .and. &
            an2(1:1) .eq. 'O' .and. &
            an3(1:1) .eq. 'H' .or.  &
            an1(1:1) .eq. 'h' .and. &
            an2(1:1) .eq. 'o' .and. &
            an3(1:1) .eq. 'h') then

          nwater = nwater + 1
          constraints%water_list(1,nwater) = i+1
          constraints%water_list(2,nwater) = i
          constraints%water_list(3,nwater) = i+2
          lwater(i:i+2) = .true.

        else if ( &
            an1(1:1) .eq. 'O' .and. &
            an2(1:1) .eq. 'H' .and. &
            an3(1:1) .eq. 'H' .or.  &
            an1(1:1) .eq. 'o' .and. &
            an2(1:1) .eq. 'h' .and. &
            an3(1:1) .eq. 'h') then

          nwater = nwater + 1
          constraints%water_list(1,nwater) = i
          constraints%water_list(2,nwater) = i+1
          constraints%water_list(3,nwater) = i+2
          lwater(i:i+2) = .true.

        else if ( &
            an1(1:1) .eq. 'H' .and. &
            an2(1:1) .eq. 'H' .and. &
            an3(1:1) .eq. 'O' .or.  &
            an1(1:1) .eq. 'h' .and. &
            an2(1:1) .eq. 'h' .and. &
            an3(1:1) .eq. 'o') then

          nwater = nwater + 1
          constraints%water_list(1,nwater) = i+2
          constraints%water_list(2,nwater) = i
          constraints%water_list(3,nwater) = i+1
          lwater(i:i+2) = .true.

        end if

        i = i + 3

      else

        i = i + 1

      end if

    end do

    if (nwater /= constraints%num_water) &
      call error_msg('Setup_Settle> number of water is incorrect')

    constraints%water_massO = molecule%mass(constraints%water_list(1,1))
    constraints%water_massH = molecule%mass(constraints%water_list(2,1))


    ! reduce bond list
    !
    allocate(temp_list(2,enefunc%num_bonds), &
             temp_fc    (enefunc%num_bonds), &
             temp_min   (enefunc%num_bonds), stat=alloc_stat)

    nbond = 0

    do i = 1, enefunc%num_bonds

      i1 = enefunc%bond_list(1,i)
      i2 = enefunc%bond_list(2,i)

      if (.not. lwater(i1) .and. .not. lwater(i2)) then

        nbond = nbond + 1
        temp_list(1:2,nbond) = enefunc%bond_list(1:2,i)
        temp_fc  (    nbond) = enefunc%bond_force_const(i)
        temp_min (    nbond) = enefunc%bond_dist_min(i)

      else if (lwater(i1) .and. lwater(i2)) then

        an1 = molecule%atom_name(i1)
        an2 = molecule%atom_name(i2)

        if (an1(1:1) .eq. 'H' .and. an2(1:1) .eq. 'H' .or. &
            an1(1:1) .eq. 'h' .and. an2(1:1) .eq. 'h') then

          constraints%water_rHH = enefunc%bond_dist_min(i)

        else if &
           (an1(1:1) .eq. 'H' .and. an2(1:1) .eq. 'O' .or. &
            an1(1:1) .eq. 'O' .and. an2(1:1) .eq. 'H' .or. &
            an1(1:1) .eq. 'h' .and. an2(1:1) .eq. 'o' .or. &
            an1(1:1) .eq. 'o' .and. an2(1:1) .eq. 'h') then

          constraints%water_rOH = enefunc%bond_dist_min(i)

        end if

      else
        call error_msg('Setup_Settle> bond parameter is incorrect.')

      end if

    end do

    call alloc_enefunc(enefunc, EnefuncBond, nbond)

    enefunc%num_bonds = nbond
    enefunc%bond_list   (1:2,1:nbond) = temp_list(1:2,1:nbond)
    enefunc%bond_force_const(1:nbond) = temp_fc  (    1:nbond)
    enefunc%bond_dist_min   (1:nbond) = temp_min (    1:nbond)
    
    deallocate(temp_list, temp_fc, temp_min, stat=alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_dealloc


    ! reduce angle list
    !
    allocate(temp_list(3,enefunc%num_angles), &
             temp_fc    (enefunc%num_angles), &
             temp_min   (enefunc%num_angles), stat=alloc_stat)

    nangl = 0

    do i = 1, enefunc%num_angles

      i1 = enefunc%angl_list(1,i)
      i2 = enefunc%angl_list(2,i)
      i3 = enefunc%angl_list(3,i)

      if (.not. lwater(i1) .and. .not. lwater(i2) .and. .not. lwater(i3)) then

        nangl = nangl + 1
        temp_list(1:3,nangl) = enefunc%angl_list(1:3,i)
        temp_fc  (    nangl) = enefunc%angl_force_const(i)
        temp_min (    nangl) = enefunc%angl_theta_min(i)

      else if (lwater(i1) .and. lwater(i2) .and. lwater(i3)) then
        if (constraints%water_rHH < 1e-5_wp) then
          constraints%water_rHH = 2.0_wp*constraints%water_rOH   &
                                  * sin(0.5_wp*enefunc%angl_theta_min(i))
        endif

      else
        call error_msg('Setup_Settle> angle parameter is incorrect.')

      end if

    end do

    call alloc_enefunc(enefunc, EnefuncAngl, nangl)

    enefunc%num_angles = nangl
    enefunc%angl_list   (1:3,1:nangl) = temp_list(1:3,1:nangl)
    enefunc%angl_force_const(1:nangl) = temp_fc  (    1:nangl)
    enefunc%angl_theta_min  (1:nangl) = temp_min (    1:nangl)
    
    deallocate(temp_list, temp_fc, temp_min, lwater, stat=alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_dealloc


    ! re-setup loop index
    !
    call get_loop_index(nbond, istart, iend)
    enefunc%istart_bond = istart
    enefunc%iend_bond   = iend

    call get_loop_index(nangl, istart, iend)
    enefunc%istart_angle = istart
    enefunc%iend_angle   = iend


    ! allocate virial_tmp 
    !
#ifdef OMP
    !$omp parallel default(none) shared(nthread)
    nthread = omp_get_num_threads()
    !$omp end parallel
#else
    nthread = 1
#endif
    allocate(constraints%virial_tmp(3,3,nthread), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc


    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') &
           'Setup_Settle> Setup constraints for SETTLE'
      write(MsgOut,'(A20,I10,A20,I10)')                   &
           '  num_unconsbonds = ', enefunc%num_bonds,     &
           '  num_Settle  = ', constraints%num_water
      write(MsgOut,'(A20,F10.4,A20,F10.4)')               &
           '  r0(O-H)         = ', constraints%water_rOH, &
           '  mass(O)         = ', constraints%water_massO
      write(MsgOut,'(A20,F10.4,A20,F10.4)')               & 
           '  r0(H-H)         = ', constraints%water_rHH, &
           '  mass(H)         = ', constraints%water_massH
      write(MsgOut,'(A)') ' '
    end if

    ! error check
    !
    if (constraints%water_rOH   == 0.0_wp .or.  &
        constraints%water_rHH   == 0.0_wp .or.  &
        constraints%water_massO == 0.0_wp .or.  &
        constraints%water_massH == 0.0_wp) &
      call error_msg( &
           'Setup_Settle> Error in defining three-site water!')

    return

  end subroutine setup_settle

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_fixatm
  !> @brief        setup atoms fixed in space
  !! @authors      KY
  !! @param[in]    molecule    : molecule information
  !! @param[in]    boundary    : information of boundary condition
  !! @param[in]    dynamics    : dynamics information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[in]    sel_info    : selector input information
  !! @param[in]    cons_info   : CONSTRAINTS section control parameters
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_fixatm(molecule, boundary, dynamics, enefunc, &
                          sel_info, cons_info, constraints)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_boundary),        intent(inout) :: boundary
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_enefunc), target, intent(inout) :: enefunc
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_cons_info),       intent(in)    :: cons_info
    type(s_constraints),     intent(inout) :: constraints

    ! local variables
    integer               :: i, j, n
    type(s_qmmm), pointer :: qmmm

    integer                       :: ngroup, igroup
    integer,          allocatable :: group_list(:)
    type(s_selatoms), allocatable :: selatoms(:)


    qmmm => enefunc%qmmm

    call alloc_constraints(constraints, ConstraintsFixatm, molecule%num_atoms)

    ! fixed atoms at the boundary
    if (boundary%num_fixatm > 0) then
      constraints%num_fixatm = boundary%num_fixatm
      constraints%fixatm     = boundary%fixatm
      call update_num_deg_freedom('After setup of fixed atoms at boundary', &
                           -3*boundary%num_fixatm, molecule%num_deg_freedom)
    end if

    ! fix QM atoms for ESP/MM-MD
    if (dynamics%esp_mm) then
      n = 0
      do i = 1, qmmm%qm_natoms
        if (.not. constraints%fixatm(qmmm%qmatom_id(i))) then
          n = n + 1
          constraints%fixatm(qmmm%qmatom_id(i)) = .true.
          constraints%num_fixatm = constraints%num_fixatm + 1
        end if
      end do
      call update_num_deg_freedom('After setup of fixed QM atoms', &
                                   -3*n, molecule%num_deg_freedom)
    end if

    ! fixed atoms from the input
    !
    ngroup = split_num(trim(cons_info%fixatm_select_index))
    if (ngroup /= 0) then
      allocate(group_list(ngroup), selatoms(ngroup))
      call split(ngroup, ngroup, cons_info%fixatm_select_index, group_list)

      n = 0
      do i = 1, ngroup
        igroup = group_list(i)
        call select_atom(molecule, sel_info%groups(igroup), selatoms(i))
        do j = 1, size(selatoms(i)%idx)
          if (.not. constraints%fixatm(selatoms(i)%idx(j))) then
            !dbg if (main_rank) write(MsgOut,'("debug: ",i8)') selatoms(i)%idx(j)
            n = n + 1
            constraints%fixatm(selatoms(i)%idx(j)) = .true.
            constraints%num_fixatm = constraints%num_fixatm + 1
          end if
        end do
      end do
      deallocate(group_list, selatoms)

      call update_num_deg_freedom('After setup of fixed atoms', &
                                   -3*n, molecule%num_deg_freedom)

    end if


    if (constraints%num_fixatm > 0) &
      call remove_fixatm(molecule%num_atoms, constraints%fixatm, enefunc)

    return

  end subroutine setup_fixatm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    remove_fixatm
  !> @brief        remove potential functions for fixed atoms
  !! @authors      KY
  !! @param[in]    num_atoms : number of atoms
  !! @param[in]    fixatm    : information of fixed atom
  !! @param[inout] enefunc   : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine remove_fixatm(num_atoms, fixatm, enefunc)

    ! formal arguments
    integer,         intent(in)    :: num_atoms
    logical,         intent(in)    :: fixatm(num_atoms)
    type(s_enefunc), intent(inout) :: enefunc

    ! bond
    !
    if (enefunc%num_bonds > 0) &
      call remove_fixatm_bond(num_atoms, fixatm, enefunc)

    ! angle
    !
    if (enefunc%num_angles > 0) &
      call remove_fixatm_angl(num_atoms, fixatm, enefunc)

    ! Urey-Bradley
    !
    if (enefunc%num_ureys > 0) &
      call remove_fixatm_urey(num_atoms, fixatm, enefunc)

    ! dihedral
    !
    if (enefunc%num_dihedrals > 0) &
      call remove_fixatm_dihe(num_atoms, fixatm, enefunc)

    ! RB dihedral
    !
    if (enefunc%num_rb_dihedrals > 0) &
      call remove_fixatm_rb_dihe(num_atoms, fixatm, enefunc)

    ! improper
    !
    if (enefunc%num_impropers > 0) &
      call remove_fixatm_impr(num_atoms, fixatm, enefunc)

    ! cmap
    !
    if (enefunc%num_cmaps > 0) &
      call remove_fixatm_cmap(num_atoms, fixatm, enefunc)

    if (main_rank) then

      write(MsgOut,'(A)') &
           'Remove_Fixatm> Modified Interactions'
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           '  bond_ene        = ', enefunc%num_bonds,           &
           '  angle_ene       = ', enefunc%num_angles
      if (enefunc%num_ureys > 0)   &
        write(MsgOut,'(A20,I10)')  &
           '  urey_ene        = ', enefunc%num_ureys
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           '  torsion_ene     = ', enefunc%num_dihedrals,       &
           '  improper_ene    = ', enefunc%num_impropers
      if (enefunc%num_rb_dihedrals > 0) &
        write(MsgOut,'(A20,I10)')       &
           '  rb_dihed_ene    = ', enefunc%num_rb_dihedrals
      !write(MsgOut,'(A20,I10)')                                 &
      !     '  cmap_ene        = ', enefunc%num_cmaps
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine remove_fixatm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    remove_fixatm_bond
  !> @brief        remove bond functions for fixed atoms
  !! @authors      KY
  !! @param[in]    num_atoms : number of atoms
  !! @param[in]    fixatm    : information of fixed atom
  !! @param[inout] enefunc   : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine remove_fixatm_bond(num_atoms, fixatm, enefunc)

    ! formal arguments
    integer,         intent(in)    :: num_atoms
    logical,         intent(in)    :: fixatm(num_atoms)
    type(s_enefunc), intent(inout) :: enefunc

    ! local variables
    integer    :: i, j, k
    integer    :: nbonds, nbonds_new
    integer    :: found
    integer    :: istart, iend
    integer    :: alloc_stat, dealloc_stat

    real(wp),  allocatable   :: temp_fc(:), temp_r0(:)
    integer,   allocatable   :: temp_list(:,:)


    nbonds = enefunc%num_bonds

    ! count the number of bonds including fixed atoms
    !
    found = 0
    do i = 1, nbonds
      if (fixatm(enefunc%bond_list(1,i)) .and. &
          fixatm(enefunc%bond_list(2,i)))      &
        found = found + 1
    end do

    if (found == 0) return
    nbonds_new = nbonds - found

    ! copy information from enefunc to temp array
    !
    allocate(temp_fc(nbonds), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_r0(nbonds), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_list(2,nbonds), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    temp_fc   = enefunc%bond_force_const
    temp_r0   = enefunc%bond_dist_min
    temp_list = enefunc%bond_list

    ! re-allocate EneFuncBond
    !
    call dealloc_enefunc(enefunc, EneFuncBond)

    enefunc%num_bonds = nbonds_new
    call alloc_enefunc(enefunc, EneFuncBond, nbonds_new)

    ! modify enefunc%bond
    !
    found = 0
    do i = 1, nbonds
      if (.not. (fixatm(temp_list(1,i)) .and. &
                 fixatm(temp_list(2,i)))) then
        found = found + 1
        enefunc%bond_force_const(found) = temp_fc(i)
        enefunc%bond_dist_min(found)    = temp_r0(i)
        enefunc%bond_list(:,found)      = temp_list(:,i)
      end if

    end do

    ! re-setup loop index
    !
    call get_loop_index(enefunc%num_bonds, istart, iend)
    enefunc%istart_bond = istart
    enefunc%iend_bond   = iend

    ! deallocate temp array
    !
    deallocate(temp_fc, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_r0, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_list, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine remove_fixatm_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    remove_fixatm_angl
  !> @brief        remove angle and UB functions for fixed atoms
  !! @authors      KY
  !! @param[in]    num_atoms : number of atoms
  !! @param[in]    fixatm    : information of fixed atom
  !! @param[inout] enefunc  : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine remove_fixatm_angl(num_atoms, fixatm, enefunc)

    ! formal arguments
    integer,         intent(in)    :: num_atoms
    logical,         intent(in)    :: fixatm(num_atoms)
    type(s_enefunc), intent(inout) :: enefunc

    ! local variables
    integer    :: i, j, k
    integer    :: nangles, nangles_new
    integer    :: found
    integer    :: istart, iend
    integer    :: alloc_stat, dealloc_stat

    real(wp),  allocatable   :: temp_fc(:), temp_r0(:)
    integer,   allocatable   :: temp_list(:,:)


    nangles = enefunc%num_angles

    ! count the number of angles including fixed atoms
    !
    found = 0
    do i = 1, nangles
      if (fixatm(enefunc%angl_list(1,i)) .and. &
          fixatm(enefunc%angl_list(2,i)) .and. &
          fixatm(enefunc%angl_list(3,i)))      &
        found = found + 1
    end do

    if (found == 0) return

    nangles_new = nangles - found

    ! copy information from enefunc to temp array
    !
    allocate(temp_fc(nangles), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_r0(nangles), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_list(3,nangles), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    temp_fc   = enefunc%angl_force_const
    temp_r0   = enefunc%angl_theta_min
    temp_list = enefunc%angl_list

    ! re-allocate EneFuncAngl
    !
    call dealloc_enefunc(enefunc, EneFuncAngl)

    enefunc%num_angles = nangles_new
    call alloc_enefunc(enefunc, EneFuncAngl, nangles_new)

    ! modify enefunc%angle
    !
    found = 0
    do i = 1, nangles
      if (.not. (fixatm(temp_list(1,i)) .and. &
                 fixatm(temp_list(2,i)) .and. &
                 fixatm(temp_list(3,i)))) then
        found = found + 1
        enefunc%angl_force_const(found) = temp_fc(i)
        enefunc%angl_theta_min(found)   = temp_r0(i)
        enefunc%angl_list(:,found)      = temp_list(:,i)
      end if

    end do

    ! re-setup loop index
    !
    call get_loop_index(enefunc%num_angles, istart, iend)
    enefunc%istart_angle = istart
    enefunc%iend_angle   = iend

    ! deallocate temp array
    !
    deallocate(temp_fc, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_r0, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_list, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine remove_fixatm_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    remove_fixatm_urey
  !> @brief        remove Urey-Bradley functions for fixed atoms
  !! @authors      KY
  !! @param[in]    num_atoms : number of atoms
  !! @param[in]    fixatm    : information of fixed atom
  !! @param[inout] enefunc   : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine remove_fixatm_urey(num_atoms, fixatm, enefunc)

    ! formal arguments
    integer,         intent(in)    :: num_atoms
    logical,         intent(in)    :: fixatm(num_atoms)
    type(s_enefunc), intent(inout) :: enefunc

    ! local variables
    integer    :: i, j, k
    integer    :: nureys, nureys_new
    integer    :: found
    integer    :: istart, iend
    integer    :: alloc_stat, dealloc_stat

    real(wp),  allocatable   :: temp_fc(:), temp_r0(:)
    integer,   allocatable   :: temp_list(:,:)


    nureys = enefunc%num_ureys

    ! count the number of ureys including fixed atoms
    !
    found = 0
    do i = 1, nureys
      if (fixatm(enefunc%urey_list(1,i)) .and. &
          fixatm(enefunc%urey_list(2,i)))      &
        found = found + 1
    end do

    if (found == 0) return

    nureys_new = nureys - found

    ! copy information from enefunc to temp array
    !
    allocate(temp_fc(nureys), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_r0(nureys), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_list(3,nureys), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    temp_fc   = enefunc%urey_force_const
    temp_r0   = enefunc%urey_rmin
    temp_list = enefunc%urey_list

    ! re-allocate EneFuncUrey
    !
    call dealloc_enefunc(enefunc, EneFuncUrey)

    enefunc%num_ureys = nureys_new
    call alloc_enefunc(enefunc, EneFuncUrey, nureys_new)

    ! modify enefunc%urey
    !
    found = 0
    do i = 1, nureys
      if (.not. (fixatm(temp_list(1,i)) .and. &
                 fixatm(temp_list(2,i)))) then
        found = found + 1
        enefunc%urey_force_const(found) = temp_fc(i)
        enefunc%urey_rmin(found)        = temp_r0(i)
        enefunc%urey_list(:,found)      = temp_list(:,i)
      end if

    end do

    ! re-setup loop index
    !
    call get_loop_index(enefunc%num_ureys, istart, iend)
    enefunc%istart_urey = istart
    enefunc%iend_urey   = iend

    ! deallocate temp array
    !
    deallocate(temp_fc, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_r0, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_list, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine remove_fixatm_urey

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    remove_fixatm_dihe
  !> @brief        remove dihedral functions for fixed atoms
  !! @authors      KY
  !! @param[in]    num_atoms : number of atoms
  !! @param[in]    fixatm    : information of fixed atom
  !! @param[inout] enefunc   : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine remove_fixatm_dihe(num_atoms, fixatm, enefunc)

    ! formal arguments
    integer,         intent(in)    :: num_atoms
    logical,         intent(in)    :: fixatm(num_atoms)
    type(s_enefunc), intent(inout) :: enefunc

    ! local variables
    integer    :: i, j, k
    integer    :: ndihedrals, ndihedrals_new
    integer    :: found
    integer    :: istart, iend
    integer    :: alloc_stat, dealloc_stat

    real(wp),  allocatable   :: temp_fc(:), temp_phase(:)
    real(wp),  allocatable   :: temp_scee(:), temp_scnb(:)
    integer,   allocatable   :: temp_list(:,:), temp_period(:)


    ndihedrals = enefunc%num_dihedrals

    ! count the number of dihedrals including fixed atoms
    !
    found = 0
    do i = 1, ndihedrals
      if (fixatm(enefunc%dihe_list(1,i)) .and. &
          fixatm(enefunc%dihe_list(2,i)) .and. &
          fixatm(enefunc%dihe_list(3,i)) .and. &
          fixatm(enefunc%dihe_list(4,i)))      &
        found = found + 1
    end do

    if (found == 0) return

    ndihedrals_new = ndihedrals - found

    ! copy information from enefunc to temp array
    !
    allocate(temp_fc(ndihedrals), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_phase(ndihedrals), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_scee(ndihedrals), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_scnb(ndihedrals), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_period(ndihedrals), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_list(4,ndihedrals), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    temp_fc     = enefunc%dihe_force_const
    temp_phase  = enefunc%dihe_phase
    temp_scee   = enefunc%dihe_scee
    temp_scnb   = enefunc%dihe_scnb
    temp_period = enefunc%dihe_periodicity
    temp_list   = enefunc%dihe_list

    ! re-allocate EneFuncDihe
    !
    call dealloc_enefunc(enefunc, EneFuncDihe)

    enefunc%num_dihedrals = ndihedrals_new
    call alloc_enefunc(enefunc, EneFuncDihe, ndihedrals_new)

    ! modify enefunc%dihe
    !
    found = 0
    do i = 1, ndihedrals
      if (.not. (fixatm(temp_list(1,i)) .and. &
                 fixatm(temp_list(2,i)) .and. &
                 fixatm(temp_list(3,i)) .and. &
                 fixatm(temp_list(4,i)))) then
        found = found + 1
        enefunc%dihe_force_const(found) = temp_fc(i)
        enefunc%dihe_phase(found)       = temp_phase(i)
        enefunc%dihe_scee(found)        = temp_scee(i)
        enefunc%dihe_scnb(found)        = temp_scnb(i)
        enefunc%dihe_periodicity(found) = temp_period(i)
        enefunc%dihe_list(:,found)      = temp_list(:,i)
      end if

    end do

    ! re-setup loop index
    !
    call get_loop_index(enefunc%num_dihedrals, istart, iend)
    enefunc%istart_dihedral = istart
    enefunc%iend_dihedral   = iend

    ! deallocate temp array
    !
    deallocate(temp_fc, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_phase, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_scee, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_scnb, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_period, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_list, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine remove_fixatm_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    remove_fixatm_rb_dihe
  !> @brief        remove RB dihedral functions for fixed atoms
  !! @authors      KY
  !! @param[in]    num_atoms : number of atoms
  !! @param[in]    fixatm    : information of fixed atom
  !! @param[inout] enefunc   : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine remove_fixatm_rb_dihe(num_atoms, fixatm, enefunc)

    ! formal arguments
    integer,         intent(in)    :: num_atoms
    logical,         intent(in)    :: fixatm(num_atoms)
    type(s_enefunc), intent(inout) :: enefunc

    ! local variables
    integer    :: i, j, k
    integer    :: ndihedrals, ndihedrals_new
    integer    :: found
    integer    :: istart, iend
    integer    :: alloc_stat, dealloc_stat

    real(wp),  allocatable   :: temp_c(:,:)
    integer,   allocatable   :: temp_list(:,:)


    ndihedrals = enefunc%num_rb_dihedrals

    ! count the number of RB dihedrals including fixed atoms
    !
    found = 0
    do i = 1, ndihedrals
      if (fixatm(enefunc%rb_dihe_list(1,i)) .and. &
          fixatm(enefunc%rb_dihe_list(2,i)) .and. &
          fixatm(enefunc%rb_dihe_list(3,i)) .and. &
          fixatm(enefunc%rb_dihe_list(4,i)))      &
        found = found + 1
    end do

    if (found == 0) return

    ndihedrals_new = ndihedrals - found

    ! copy information from enefunc to temp array
    !
    allocate(temp_c(6,ndihedrals), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_list(4,ndihedrals), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    temp_c    = enefunc%rb_dihe_c
    temp_list = enefunc%rb_dihe_list

    ! re-allocate EneFuncRBDihe
    !
    call dealloc_enefunc(enefunc, EneFuncRBDihe)

    enefunc%num_rb_dihedrals = ndihedrals_new
    call alloc_enefunc(enefunc, EneFuncRBDihe, ndihedrals_new)

    ! modify RB enefunc%dihe
    !
    found = 0
    do i = 1, ndihedrals
      if (.not. (fixatm(temp_list(1,i)) .and. &
                 fixatm(temp_list(2,i)) .and. &
                 fixatm(temp_list(3,i)) .and. &
                 fixatm(temp_list(4,i)))) then
        found = found + 1
        enefunc%rb_dihe_c(:,found)    = temp_c(:,i)
        enefunc%rb_dihe_list(:,found) = temp_list(:,i)
      end if

    end do

    ! re-setup loop index
    !
    call get_loop_index(enefunc%num_rb_dihedrals, istart, iend)
    enefunc%istart_rb_dihed = istart
    enefunc%iend_rb_dihed   = iend

    ! deallocate temp array
    !
    deallocate(temp_c, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_list, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine remove_fixatm_rb_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    remove_fixatm_impr
  !> @brief        remove improper functions for fixed atoms
  !! @authors      KY
  !! @param[in]    num_atoms : number of atoms
  !! @param[in]    fixatm    : information of fixed atom
  !! @param[inout] enefunc   : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine remove_fixatm_impr(num_atoms, fixatm, enefunc)

    ! formal arguments
    integer,         intent(in)    :: num_atoms
    logical,         intent(in)    :: fixatm(num_atoms)
    type(s_enefunc), intent(inout) :: enefunc

    ! local variables
    integer    :: i, j, k
    integer    :: nimpropers, nimpropers_new
    integer    :: found
    integer    :: istart, iend
    integer    :: alloc_stat, dealloc_stat

    real(wp),  allocatable   :: temp_fc(:), temp_phase(:)
    integer,   allocatable   :: temp_list(:,:), temp_period(:)


    nimpropers = enefunc%num_impropers

    ! count the number of impropers including fixed atoms
    !
    found = 0
    do i = 1, nimpropers
      if (fixatm(enefunc%impr_list(1,i)) .and. &
          fixatm(enefunc%impr_list(2,i)) .and. &
          fixatm(enefunc%impr_list(3,i)) .and. &
          fixatm(enefunc%impr_list(4,i)))      &
        found = found + 1
    end do

    if (found == 0) return

    nimpropers_new = nimpropers - found

    ! copy information from enefunc to temp array
    !
    allocate(temp_fc(nimpropers), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_period(nimpropers), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_phase(nimpropers), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(temp_list(4,nimpropers), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    temp_fc     = enefunc%impr_force_const
    temp_period = enefunc%impr_periodicity
    temp_phase  = enefunc%impr_phase
    temp_list   = enefunc%impr_list

    ! re-allocate EneFuncDihe
    !
    call dealloc_enefunc(enefunc, EneFuncImpr)

    enefunc%num_impropers = nimpropers_new
    call alloc_enefunc(enefunc, EneFuncImpr, nimpropers_new)

    ! modify enefunc%dihe
    !
    found = 0
    do i = 1, nimpropers
      if (.not. (fixatm(temp_list(1,i)) .and. &
                 fixatm(temp_list(2,i)) .and. &
                 fixatm(temp_list(3,i)) .and. &
                 fixatm(temp_list(4,i)))) then
        found = found + 1
        enefunc%impr_force_const(found) = temp_fc(i)
        enefunc%impr_periodicity(found) = temp_period(i)
        enefunc%impr_phase(found)       = temp_phase(i)
        enefunc%impr_list(:,found)      = temp_list(:,i)
      end if

    end do

    ! re-setup loop index
    !
    call get_loop_index(enefunc%num_impropers, istart, iend)
    enefunc%istart_improper = istart
    enefunc%iend_improper   = iend

    ! deallocate temp array
    !
    deallocate(temp_fc, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_period, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_phase, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(temp_list, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine remove_fixatm_impr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    remove_fixatm_cmap
  !> @brief        remove CMAP functions for fixed atoms
  !! @authors      KY
  !! @param[in]    num_atoms : number of atoms
  !! @param[in]    fixatm    : information of fixed atom
  !! @param[inout] enefunc   : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine remove_fixatm_cmap(num_atoms, fixatm, enefunc)

    ! formal arguments
    integer,         intent(in)    :: num_atoms
    logical,         intent(in)    :: fixatm(num_atoms)
    type(s_enefunc), intent(inout) :: enefunc

    ! local variables
    integer    :: i, j, k
    integer    :: ncmaps, ncmaps_new
    integer    :: found
    integer    :: istart, iend
    integer    :: alloc_stat, dealloc_stat

    real(wp),  allocatable   :: temp_force(:,:,:), temp_coef(:,:,:,:,:)
    integer,   allocatable   :: temp_list(:,:), temp_res(:), temp_type(:)

    ! TODO

    return

  end subroutine remove_fixatm_cmap

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_constraints
  !> @brief        update coordinates and velocities
  !! @authors      TM, JJ
  !! @param[in]    cons_mode   : constraints mode
  !! @param[in]    dt          : time step
  !! @param[in]    molecule    : molecule information
  !! @param[inout] dynvars     : dynamic variables information
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_constraints(cons_mode, dt, molecule, dynvars, constraints)

    ! formal arguments
    integer,                 intent(in)    :: cons_mode
    real(wp),                intent(in)    :: dt
    type(s_molecule),        intent(in)    :: molecule
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_constraints),     intent(inout) :: constraints

    integer                  :: i, natom
    real(wp)                 :: vel_change(3), viri(3)

    call timer(TimerConstraint, TimerOn)

    ! constraints
    !
    select case (cons_mode)

    case (ConstraintModeLEAP)

      dynvars%virial_const(1:3,1:3) = 0.0_wp

      if (constraints%fast_water) then
        call compute_settle(dynvars%coord_ref, dynvars%coord,       &
                            dynvars%velocity, dynvars%virial_const, &
                            constraints, dt, .false.)
      end if

      if (constraints%fast_bond) then
        call compute_lincs(molecule%inv_mass,                       &
                           dynvars%coord_ref, dynvars%coord,        &
                           dynvars%virial_const, constraints)
      else
        call compute_shake(dt, molecule, molecule%inv_mass,         &
                           dynvars%coord_ref, dynvars%coord,        &
                           dynvars%velocity, dynvars%virial_const,  &
                           constraints, .false.)
      end if

      dynvars%virial_const(1:3,1:3) = dynvars%virial_const(1:3,1:3)/dt**2

    case (ConstraintModeVVER1)

      dynvars%virial_const(1:3,1:3) = 0.0_wp

      if (constraints%fast_water) then
        call compute_settle(dynvars%coord_ref, dynvars%coord,       &
                            dynvars%velocity, dynvars%virial_const, &
                            constraints, dt, .true.)
      end if

      call compute_shake(dt, molecule, molecule%inv_mass,           &
                         dynvars%coord_ref, dynvars%coord,          &
                         dynvars%velocity, dynvars%virial_const,    &
                         constraints, .true.)

      dynvars%virial_const(1:3,1:3) = dynvars%virial_const(1:3,1:3)/dt**2

    case (ConstraintModeVVER2)

      natom = molecule%num_atoms
      dynvars%velocity_ref(1:3,1:natom) = dynvars%velocity(1:3,1:natom)
      if (constraints%fast_water) then
        call compute_settle_vv2(dynvars%coord, dynvars%velocity,    &
                                constraints)
      end if
      call compute_rattle(molecule%mass, dynvars%coord,             &
                          dynvars%velocity, constraints)
      dynvars%virial_const(1:3,1:3) = 0.0_wp
      do i = 1, natom
        vel_change(1:3) = dynvars%velocity(1:3,i) - dynvars%velocity_ref(1:3,i)
        viri(1:3) = molecule%mass(i)*vel_change(1:3)*dynvars%coord(1:3,i)/dt
        viri(1:3) = viri(1:3) * 2.0_wp
        dynvars%virial_const(1,1) = dynvars%virial_const(1,1) + viri(1)
        dynvars%virial_const(2,2) = dynvars%virial_const(2,2) + viri(2)
        dynvars%virial_const(3,3) = dynvars%virial_const(3,3) + viri(3)
      end do

    case (ConstraintModeBD)

      dynvars%virial_const(1:3,1:3) = 0.0_wp

      if (constraints%fast_bond) then
        call compute_lincs(molecule%inv_stokes_radius,        &
                           dynvars%coord_ref, dynvars%coord,  &
                           dynvars%virial_const, constraints)
      end if

      if (constraints%rigid_bond) then
        call compute_shake(dt, molecule, molecule%inv_stokes_radius, &
                           dynvars%coord_ref, dynvars%coord,         &
                           dynvars%velocity, dynvars%virial_const,   &
                           constraints, .true.) !TODO check: velocity updt OK?
       end if

      dynvars%virial_const(1:3,1:3) = dynvars%virial_const(1:3,1:3)/dt**2

    end select


    call timer(TimerConstraint, TimerOff)

    return

  end subroutine compute_constraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    exclude_light_atom
  !> @brief        exclude atoms for SHAKE
  !! @authors      KY
  !! @param[in]    noshake_index : index of no-shake atoms
  !! @param[in]    sel_info  : SELECTION section control parameters information
  !! @param[inout] molecule  : molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine exclude_light_atom(noshake_index, sel_info, molecule)

    ! formal arguments
    character(MaxLine)       :: noshake_index
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_molecule),        intent(inout) :: molecule

    ! local variables
    integer                       :: i, j, k, idx, ncount
    integer                       :: ngroup, igroup, natom, iatom
    integer,          allocatable :: group_list(:)
    logical,          allocatable :: switch(:)
    type(s_selatoms), allocatable :: selatoms(:)
    integer                       :: kalloc_stat, kdealloc_stat


    ngroup = split_num(trim(noshake_index))
    allocate(group_list(ngroup), stat = kalloc_stat)
    if (kalloc_stat /= 0) call error_msg_alloc

    call split(ngroup, ngroup, noshake_index, group_list)
    allocate(selatoms(ngroup), stat = kalloc_stat)
    if (kalloc_stat /= 0) call error_msg_alloc

    natom = 0
    do i = 1, ngroup
      igroup = group_list(i)
      call select_atom(molecule, sel_info%groups(igroup), selatoms(i))
      natom  = natom + size(selatoms(i)%idx)
    end do

    if (natom == 0) then
      deallocate(selatoms, stat = kdealloc_stat)
      if (kdealloc_stat /= 0) call error_msg_dealloc
      deallocate(group_list, stat = kdealloc_stat)
      if (kdealloc_stat /= 0) call error_msg_dealloc
      return
    end if

    allocate(switch(natom), stat = kalloc_stat)
    switch = .false.

    ncount = 0
    k      = 1
    do i = 1, ngroup
      natom  = size(selatoms(i)%idx)
      do j = 1, natom
        idx = selatoms(i)%idx(j)

        if (molecule%light_atom_mass(idx) .or. &
            molecule%light_atom_name(idx) ) then
          molecule%light_atom_mass(idx) = .false.
          molecule%light_atom_name(idx) = .false.
          switch(k) = .true.
          ncount = ncount + 1
        end if
        k = k + 1

      end do
    end do

    if (main_rank .and. ncount > 0) then
      write(MsgOut,'(A)') &
        'Setup_Constraints> The following atoms are excluded from SHAKE:'
      k = 1
      do i = 1, ngroup
        natom  = size(selatoms(i)%idx)
        do j = 1, natom
          idx = selatoms(i)%idx(j)

          if (switch(k)) then
            write(MsgOut,'(i8,$)') idx
            if (mod(k,10) == 0) write(MsgOut,*)
          end if
          k = k + 1

        end do
      end do
      write(MsgOut,*)
    end if

    return

  end subroutine exclude_light_atom

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_shake
  !> @brief        SHAKE for rigid bonds
  !! @authors      TM, JJ
  !! @param[in]    dt          : delta time
  !! @param[in]    molecule    : molecule information
  !! @param[in]    inv_mass    : inverse of atomic mass
  !! @param[in]    coord_old   : reference coordinates at t
  !! @param[in]    coord       : unconstrained coordinates at t + dt
  !! @param[out]   coord       : constrained   coordinates at t + dt
  !! @param[inout] virial      : constraints virial 
  !! @param[inout] constraints : constraints information
  !! @param[in]    vel_update  : flag for update velocity
  !! @note         J.P.Ryckaert et al., J.Comput.Phys. 23, 327-341 (1977)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_shake(dt, molecule, inv_mass, coord_old, coord, vel, &
                           virial, constraints, vel_update)

    ! formal arguments
    real(wp),                    intent(in)    :: dt
    type(s_molecule), target,    intent(in)    :: molecule
    real(wp),                    intent(in)    :: inv_mass(:)
    real(wp),                    intent(in)    :: coord_old(:,:)
    real(wp),                    intent(inout) :: coord(:,:)
    real(wp),                    intent(inout) :: vel(:,:)
    real(wp),                    intent(inout) :: virial(3,3)
    type(s_constraints), target, intent(inout) :: constraints
    logical,                     intent(in)    :: vel_update

    ! local variables
    real(wp)                     :: inv_dt
    real(dp)                     :: tolerance, imass1, imass2
    real(dp)                     :: d12(3), d12_old(3)
    real(dp)                     :: dist, dist2, factor
    real(dp)                     :: g12, g12d12(3), g12d12m1(3), g12d12m2(3)
    real(dp)                     :: vxx, vyx, vzx, vyy, vzy, vzz, fx, fy, fz
    integer                      :: i, j
    integer                      :: atm1, atm2
    integer                      :: nbond, iteration
    logical                      :: shake_end
    logical, allocatable         :: shake_end_j(:)

    real(dp),            pointer :: vec(:,:)
    real(wp),            pointer :: r0(:), force(:)
    integer,             pointer :: bond_list(:,:)


    bond_list => constraints%bond_list
    r0        => constraints%bond_dist
    vec       => constraints%bond_vector
    force     => constraints%shake_force

    iteration =  constraints%shake_iteration
    tolerance =  real(constraints%shake_tolerance,dp)
    nbond     =  constraints%num_bonds

    inv_dt    = 1.0_wp/dt

    allocate(shake_end_j(nbond))

    ! initialize force and store old bond vector
    !
    do j = 1, nbond
      force(j) = 0.0_wp
      atm1 = bond_list(1,j)
      atm2 = bond_list(2,j)
      vec(1:3,j) = real(coord_old(1:3,atm1) - coord_old(1:3,atm2),dp)
    end do

    ! shake iteration
    !
    do i = 1, iteration

      shake_end   = .true.
      shake_end_j = .true.

      do j = 1, nbond

        atm1 = bond_list(1,j)
        atm2 = bond_list(2,j)

        d12(1:3) = real(coord(1:3,atm1) - coord(1:3,atm2),dp)
        dist2 = d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3)
        dist  = sqrt(dist2)

        if (abs(dist - real(r0(j),dp)) >= tolerance) then
          shake_end = .false.
          shake_end_j(j) = .false.

          imass1 = real(inv_mass(atm1),dp)
          imass2 = real(inv_mass(atm2),dp)
          d12_old(1:3) = vec(1:3,j)

          factor = ( d12(1)*d12_old(1)  &
                   + d12(2)*d12_old(2)  &
                   + d12(3)*d12_old(3)) &
                   * (imass1 + imass2)
          g12 = 0.5_dp*(dist2 - r0(j)*r0(j))/factor
          force(j) = force(j) + g12
          g12d12(1:3) = g12*d12_old(1:3)
          g12d12m1(1:3) = g12d12(1:3)*imass1
          g12d12m2(1:3) = g12d12(1:3)*imass2

          coord(1:3,atm1) = coord(1:3,atm1) - real(g12d12m1(1:3),wp)
          coord(1:3,atm2) = coord(1:3,atm2) + real(g12d12m2(1:3),wp) !TODO

          ! velocity update for VVER
          !
          if (vel_update) then
            vel(1:3,atm1) = vel(1:3,atm1) - real(g12d12m1(1:3)*inv_dt,wp)
            vel(1:3,atm2) = vel(1:3,atm2) + real(g12d12m2(1:3)*inv_dt,wp) !TODO
          end if

        end if

      end do

      if (shake_end) exit

    end do

    ! error check
    !
    if (.not. shake_end) then
      if (replica_main_rank) then
        write(MsgOut,'("Compute_Shake> SHAKE altorithm failed to converge")')
        do j = 1, nbond
          if (.not. shake_end_j(j)) then
            atm1 = bond_list(1,j)
            atm2 = bond_list(2,j)
            write(MsgOut,'(i8,x,a4,x,i6,x,a4,x,a6," - ",i8,x,a4,x,i6,x,a4,x,a6)')  &
              atm1,                         &
              molecule%segment_name(atm1),  &
              molecule%residue_no(atm1),    &
              molecule%residue_name(atm1),  &
              molecule%atom_name(atm1),     &
              atm2,                         &
              molecule%segment_name(atm2),  &
              molecule%residue_no(atm2),    &
              molecule%residue_name(atm2),  &
              molecule%atom_name(atm2)
          end if
        end do
      end if
      call error_msg('Compute_Shake> SHAKE algorithm failed to converge (see "Chapter: Trouble shooting" in the user manual)')

    end if

    do j = 1, nbond
      d12_old(1:3) = vec(1:3,j)
      g12 = force(j)

      fx  = g12 * d12_old(1)
      fy  = g12 * d12_old(2)
      fz  = g12 * d12_old(3)
      vxx = d12_old(1) * fx
      vyx = d12_old(2) * fx
      vzx = d12_old(3) * fx
      vyy = d12_old(2) * fy
      vzy = d12_old(3) * fy
      vzz = d12_old(3) * fz

      virial(1,1) = virial(1,1) - real(vxx,wp)
      virial(2,1) = virial(2,1) - real(vyx,wp)
      virial(3,1) = virial(3,1) - real(vzx,wp)
      virial(1,2) = virial(1,2) - real(vyx,wp)
      virial(2,2) = virial(2,2) - real(vyy,wp)
      virial(3,2) = virial(3,2) - real(vzy,wp)
      virial(1,3) = virial(1,3) - real(vzx,wp)
      virial(2,3) = virial(2,3) - real(vzy,wp)
      virial(3,3) = virial(3,3) - real(vzz,wp) !TODO
    end do

    return

  end subroutine compute_shake

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_settle
  !> @brief        SETTLE for rigid water
  !! @authors      TM, JJ
  !! @param[in]    coord_old   : reference coordinates at t
  !! @param[in]    coord       : unconstrained coordinates at t + dt
  !! @param[out]   coord       : constrained   coordinates at t + dt
  !! @param[inout] vel         : velocities
  !! @param[inout] virial      : constraints virial 
  !! @param[in]    constraints : constraints information
  !! @param[in]    dt          : delta t
  !! @param[in]    vel_update  : flag for velocity update
  !! @note         S.Miyamoto & P.A.Kollman, J.Comput.Chem. 13, 952-962 (1992)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_settle(coord_old, coord, vel, virial, constraints, &
                            dt, vel_update)

    ! formal arguments
    real(wp),                    intent(in)    :: coord_old(:,:)
    real(wp),                    intent(inout) :: coord(:,:)
    real(wp),                    intent(inout) :: vel(:,:)
    real(wp),                    intent(inout) :: virial(3,3)
    type(s_constraints), target, intent(in)    :: constraints
    real(wp),                    intent(in)    :: dt
    logical,                     intent(in)    :: vel_update

    ! local variables
    real(wp)          :: rOH, rHH, massO, massH
    real(wp)          :: ra, rb, rc, inv_ra, tmp, inv_massH2O, inv_dt
    real(wp)          :: Xaxis(3), Yaxis(3), Zaxis(3)
    real(wp)          :: Xaxis2(3), Yaxis2(3), Zaxis2(3)
    real(wp)          :: oh1(3), oh2(3), mtrx(3,3)
    real(wp)          :: x0(3,3), x1(3,3), x3(3,3)
    real(wp)          :: xp0(3,3), xp1(3,3), xp2(3,3), xp3(3,3)
    real(wp)          :: mass(3), com0(3), com1(3)
    real(wp)          :: dxp21(2), dxp23(2), dxp31(2)
    real(wp)          :: delt(3,3), rf(3,3,3)
    real(wp)          :: alpha, beta, gamma, al2bt2
    real(wp)          :: sin_theta, sin_phi, sin_psi
    real(wp)          :: cos_theta, cos_phi, cos_psi
    real(wp)          :: rb_cosphi, rb_sinphi
    real(wp)          :: rc_sinpsi_sinphi, rc_sinpsi_cosphi
    real(wp)          :: vxx, vyx, vzx, vxy, vyy, vzy, vxz, vyz, vzz
    integer           :: num_water
    integer           :: iw, i, j, k, iatom(3)
    integer           :: id, nthread
    integer           :: omp_get_num_threads, omp_get_thread_num
    integer           :: quotient, remainder
    integer           :: istart, iend

    real(wp), pointer :: virial_tmp(:,:,:)
    integer,  pointer :: water_list(:,:)


    water_list => constraints%water_list
    virial_tmp => constraints%virial_tmp

    num_water  =  constraints%num_water
    rOH        =  constraints%water_rOH
    rHH        =  constraints%water_rHH
    massO      =  constraints%water_massO
    massH      =  constraints%water_massH

    mass(1) = massO
    mass(2) = massH
    mass(3) = massH
    inv_massH2O = 1.0_wp/(massO + 2.0_wp*massH)
    rc = 0.5_wp * rHH
    ra = 2.0_wp * sqrt(rOH*rOH - rc*rc)*massH*inv_massH2O
    rb = sqrt(rOH*rOH - rc*rc) - ra
    inv_ra = 1.0_wp/ra
    inv_dt = 1.0_wp/dt

    virial_tmp(:,:,:) = 0.0_wp

    !$omp parallel default(none)                                               &
    !$omp private(iw, i, j, k, iatom, tmp, x0, x1, x3, xp0, xp1, xp2, xp3,     &
    !$omp         Xaxis, Yaxis, Zaxis, Xaxis2, Yaxis2, Zaxis2, oh1, oh2, mtrx, &
    !$omp         com0, com1, dxp21, dxp23, dxp31, alpha, beta, gamma,  al2bt2,&
    !$omp         sin_theta, sin_phi, sin_psi, cos_theta, cos_phi, cos_psi,    &
    !$omp         rb_cosphi, rb_sinphi, rc_sinpsi_sinphi, rc_sinpsi_cosphi,    &
    !$omp         delt, rf, vxx, vyx, vzx, vxy, vyy, vzy, vxz, vyz, vzz,       &
    !$omp         id, quotient, remainder, istart, iend)                       &
    !$omp shared (num_water, water_list, coord_old, coord, ra, rb, rc, mass,   &
    !$omp         inv_ra, inv_dt, rHH, vel_update, vel, virial_tmp, nthread,   &
    !$omp         inv_massH2O)
    ! 
    ! note:reduction cannot be used for "virial" because of the rounding error, 
    ! which is very large especialy in the case of langevin NPT.
    !
#ifdef OMP
    ! get loop index of each thread
    !
    nthread = omp_get_num_threads()
    id      = omp_get_thread_num()

    quotient  = num_water / nthread
    remainder = mod(num_water, nthread)

    if (id <= (remainder - 1)) then
      quotient = quotient + 1
      istart   = quotient * id + 1
      iend     = istart + quotient - 1
    else
      istart   = (quotient + 1)*remainder + quotient*(id - remainder) + 1
      iend     = istart + quotient - 1
    end if
#else
    istart  = 1
    iend    = num_water
    nthread = 1
    id      = 0
#endif

    do iw = istart, iend

      com0(1:3) = 0.0_wp
      com1(1:3) = 0.0_wp
      do k = 1, 3
        iatom(k)  = water_list(k,iw)
        x0(k,1:3) = coord_old(1:3,iatom(k))
        x1(k,1:3) = coord    (1:3,iatom(k))
        com0(1:3) = com0(1:3) + x0(k,1:3)*mass(k)
        com1(1:3) = com1(1:3) + x1(k,1:3)*mass(k)
      end do
      com0(1:3) = com0(1:3)*inv_massH2O
      com1(1:3) = com1(1:3)*inv_massH2O

      do k = 1, 3
        x0(k,1:3) = x0(k,1:3) - com0(1:3)
        x1(k,1:3) = x1(k,1:3) - com1(1:3)
      end do
      oh1(1:3) = x0(2,1:3) - x0(1,1:3)
      oh2(1:3) = x0(3,1:3) - x0(1,1:3)
      Zaxis(1) = oh1(2)*oh2(3) - oh1(3)*oh2(2)
      Zaxis(2) = oh1(3)*oh2(1) - oh1(1)*oh2(3)
      Zaxis(3) = oh1(1)*oh2(2) - oh1(2)*oh2(1)
      Xaxis(1) = x1(1,2)*Zaxis(3) - x1(1,3)*Zaxis(2)
      Xaxis(2) = x1(1,3)*Zaxis(1) - x1(1,1)*Zaxis(3)
      Xaxis(3) = x1(1,1)*Zaxis(2) - x1(1,2)*Zaxis(1)
      Yaxis(1) = Zaxis(2)*Xaxis(3) - Zaxis(3)*Xaxis(2)
      Yaxis(2) = Zaxis(3)*Xaxis(1) - Zaxis(1)*Xaxis(3)
      Yaxis(3) = Zaxis(1)*Xaxis(2) - Zaxis(2)*Xaxis(1)

      Xaxis2(1:3) = Xaxis(1:3)*Xaxis(1:3)
      Yaxis2(1:3) = Yaxis(1:3)*Yaxis(1:3)
      Zaxis2(1:3) = Zaxis(1:3)*Zaxis(1:3)
      mtrx(1:3,1) = Xaxis(1:3)/sqrt(Xaxis2(1) + Xaxis2(2) + Xaxis2(3))
      mtrx(1:3,2) = Yaxis(1:3)/sqrt(Yaxis2(1) + Yaxis2(2) + Yaxis2(3))
      mtrx(1:3,3) = Zaxis(1:3)/sqrt(Zaxis2(1) + Zaxis2(2) + Zaxis2(3))

      do k = 1, 3
        do i = 1, 3
          xp0(k,i) = 0.0_wp
          xp1(k,i) = 0.0_wp
          do j = 1, 3
            xp0(k,i) = xp0(k,i) + mtrx(j,i)*x0(k,j)
            xp1(k,i) = xp1(k,i) + mtrx(j,i)*x1(k,j)
          end do
        end do
      end do

      sin_phi = xp1(1,3)*inv_ra
      tmp = 1.0_wp - sin_phi*sin_phi
      if (tmp > 0.0_wp) then
        cos_phi = sqrt(tmp)
      else
        cos_phi = 0.0_wp
      end if

      sin_psi = (xp1(2,3) - xp1(3,3))/(rHH*cos_phi)
      tmp = 1.0_wp - sin_psi*sin_psi
      if (tmp > 0.0_wp) then
        cos_psi = sqrt(tmp)
      else
        cos_psi = 0.0_wp
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
      tmp = 1.0_wp - sin_theta*sin_theta
      if (tmp > 0.0_wp) then
        cos_theta = sqrt(tmp)
      else
        cos_theta = 0.0_wp
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

      do k = 1, 3
        do i = 1, 3
          x3(k,i) = 0.0_wp
          do j = 1, 3
            x3(k,i) = x3(k,i) + mtrx(i,j)*xp3(k,j)
          end do
        end do
        coord(1:3,iatom(k)) = x3(k,1:3) + com1(1:3)
      end do

      do k = 1, 3
        delt(k,1:3) = x3(k,1:3) - x1(k,1:3)
        do i = 1, 3
          do j = 1, 3
            rf(k,i,j) = - coord_old(i,iatom(k))*delt(k,j)*mass(k)
          end do
        end do
      end do

      vxx = rf(1,1,1) + rf(2,1,1) + rf(3,1,1)
      vyx = rf(1,2,1) + rf(2,2,1) + rf(3,2,1)
      vzx = rf(1,3,1) + rf(2,3,1) + rf(3,3,1)
      vxy = rf(1,1,2) + rf(2,1,2) + rf(3,1,2)
      vyy = rf(1,2,2) + rf(2,2,2) + rf(3,2,2)
      vzy = rf(1,3,2) + rf(2,3,2) + rf(3,3,2)
      vxz = rf(1,1,3) + rf(2,1,3) + rf(3,1,3)
      vyz = rf(1,2,3) + rf(2,2,3) + rf(3,2,3)
      vzz = rf(1,3,3) + rf(2,3,3) + rf(3,3,3)

      virial_tmp(1,1,id+1) = virial_tmp(1,1,id+1) - vxx
      virial_tmp(2,1,id+1) = virial_tmp(2,1,id+1) - vyx
      virial_tmp(3,1,id+1) = virial_tmp(3,1,id+1) - vzx
      virial_tmp(1,2,id+1) = virial_tmp(1,2,id+1) - vxy
      virial_tmp(2,2,id+1) = virial_tmp(2,2,id+1) - vyy
      virial_tmp(3,2,id+1) = virial_tmp(3,2,id+1) - vzy
      virial_tmp(1,3,id+1) = virial_tmp(1,3,id+1) - vxz
      virial_tmp(2,3,id+1) = virial_tmp(2,3,id+1) - vyz
      virial_tmp(3,3,id+1) = virial_tmp(3,3,id+1) - vzz

      ! velocity update for VVER
      !
      if (vel_update) then
        do k = 1, 3
          vel(1:3,iatom(k)) = vel(1:3,iatom(k)) + delt(k,1:3)*inv_dt
        end do
      end if

    end do
    !$omp end parallel

    do i = 1, nthread
      virial(1:3,1:3) = virial(1:3,1:3) + virial_tmp(1:3,1:3,i)
    end do

    return

  end subroutine compute_settle

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_lincs
  !> @brief        LINCS for rigid bonds
  !! @authors      TM
  !! @param[in]    inv_mass    : inversion of atomic mass
  !! @param[in]    coord_old   : reference coordinates at t
  !! @param[in]    coord       : unconstrained coordinates at t + dt
  !! @param[out]   coord       : constrained   coordinates at t + dt
  !! @param[inout] virial      : constraints virial 
  !! @param[inout] constraints : constraints information
  !! @note         B.Hess et al., J.Comput.Chem. 18, 1463-1472 (1997)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_lincs(inv_mass, coord_old, coord, virial, constraints)

    ! formal arguments
    real(wp),                    intent(in)    :: inv_mass(:)
    real(wp),                    intent(in)    :: coord_old(:,:)
    real(wp),                    intent(inout) :: coord(:,:)
    real(wp),                    intent(inout) :: virial(3,3)
    type(s_constraints), target, intent(inout) :: constraints

    ! local variables
    real(wp)                     :: vxx, vyx, vzx, vyy, vzy, vzz, fx, fy, fz
    real(wp)                     :: x12_old, y12_old, z12_old
    real(wp)                     :: p2, g12, frc, l2, d12(3)
    integer                      :: i, j, n, m, k, w
    integer                      :: a1, a2, nbond, nrec, niter

    real(dp),            pointer :: B(:,:) !TODO
    real(wp),            pointer :: blen(:), force(:), r0(:), coef(:,:)
    real(wp),            pointer :: sdiag(:), rhs(:,:), sol(:), A(:,:)
    integer,             pointer :: list(:,:), ncc(:), con(:,:)


    ! use pointers
    !
    nbond =  constraints%num_bonds
    niter =  constraints%lincs_iteration
    nrec  =  constraints%lincs_order
    list  => constraints%bond_list
    r0    => constraints%bond_dist
    ncc   => constraints%num_connected_cons
    con   => constraints%connected_cons_idx
    coef  => constraints%lincs_coef
    sdiag => constraints%s_diagonal
    rhs   => constraints%right_hand_side
    sol   => constraints%solution_array
    A     => constraints%cons_coupl_mat
    B     => constraints%bond_vector
    force => constraints%shake_force
    blen  => constraints%bond_length

    ! setup
    !
    do i = 1, nbond
      a1 = list(1,i)
      a2 = list(2,i)
      B(1:3,i) = coord_old(1:3,a1) - coord_old(1:3,a2)
      blen(i)  = sqrt(B(1,i)*B(1,i) + B(2,i)*B(2,i) + B(3,i)*B(3,i))
      B(1:3,i) = B(1:3,i)/blen(i)
    end do

    do i = 1, nbond
      do n = 1, ncc(i)
        k = con(n,i)
        A(n,i) = coef(n,i)*(B(1,i)*B(1,k) + B(2,i)*B(2,k) + B(3,i)*B(3,k))
      end do

      a1 = list(1,i)
      a2 = list(2,i)
      rhs(1,i) = sdiag(i)*( B(1,i)*(coord(1,a1) - coord(1,a2))          &
                          + B(2,i)*(coord(2,a1) - coord(2,a2))          &
                          + B(3,i)*(coord(3,a1) - coord(3,a2)) - r0(i))
      sol(i) = rhs(1,i)
    end do

    ! solve the equation
    !
    w = 2
    do m = 1, nrec
      do i = 1, nbond
        rhs(w,i) = 0.0_wp
        do n = 1, ncc(i)
          rhs(w,i) = rhs(w,i) + A(n,i)*rhs(3-w,con(n,i))
        end do
        sol(i) = sol(i) + rhs(w,i)
      end do
      w = 3 - w
    end do

    do i = 1, nbond
      a1 = list(1,i)
      a2 = list(2,i)
      force(i) = sdiag(i)*sol(i)
      coord(1:3,a1) = coord(1:3,a1) - inv_mass(a1)*B(1:3,i)*force(i)
      coord(1:3,a2) = coord(1:3,a2) + inv_mass(a2)*B(1:3,i)*force(i)
    end do

    ! Correction for rotational lengthening
    !
    do j = 1, niter
      do i = 1, nbond
        a1 = list(1,i)
        a2 = list(2,i)
        d12(1:3) = coord(1:3,a1) - coord(1:3,a2)

        l2 = d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3)
        p2 = 2.0_wp*r0(i)*r0(i) - l2
        if (p2 < 0.0_wp) p2 = 0.0_wp

        rhs(1,i) = sdiag(i)*(r0(i) - sqrt(p2))
        sol(i)   = rhs(1,i)
      end do

      ! solve the equation
      !
      w = 2
      do m = 1, nrec
        do i = 1, nbond
          rhs(w,i) = 0.0_wp
          do n = 1, ncc(i)
            rhs(w,i) = rhs(w,i) + A(n,i)*rhs(3-w,con(n,i))
          end do
          sol(i) = sol(i) + rhs(w,i)
        end do
        w = 3 - w
      end do

      do i = 1, nbond
        a1  = list(1,i)
        a2  = list(2,i)
        frc = sdiag(i)*sol(i)
        force(i) = force(i) + frc
        coord(1:3,a1) = coord(1:3,a1) - inv_mass(a1)*B(1:3,i)*frc
        coord(1:3,a2) = coord(1:3,a2) + inv_mass(a2)*B(1:3,i)*frc
      end do
    end do

    do i = 1, nbond
      g12 = force(i) * blen(i)
      fx  = g12 * B(1,i)
      fy  = g12 * B(2,i)
      fz  = g12 * B(3,i)

      vxx = B(1,i) * fx
      vyx = B(2,i) * fx
      vzx = B(3,i) * fx
      vyy = B(2,i) * fy
      vzy = B(3,i) * fy
      vzz = B(3,i) * fz

      virial(1,1) = virial(1,1) - vxx
      virial(2,1) = virial(2,1) - vyx
      virial(3,1) = virial(3,1) - vzx
      virial(1,2) = virial(1,2) - vyx
      virial(2,2) = virial(2,2) - vyy
      virial(3,2) = virial(3,2) - vzy
      virial(1,3) = virial(1,3) - vzx
      virial(2,3) = virial(2,3) - vzy
      virial(3,3) = virial(3,3) - vzz
    end do

    return

  end subroutine compute_lincs

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_rattle
  !> @brief        RATTLE for rigid bonds
  !! @authors      JJ, TM
  !! @param[in]    mass        : atomic mass
  !! @param[in]    coord       : coordinates at t + dt
  !! @param[inout] vel         : update velocities at t + dt
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_rattle(mass, coord, vel, constraints)

    ! formal arguments
    real(wp),                    intent(in)    :: mass(:)
    real(wp),                    intent(in)    :: coord(:,:)
    real(wp),                    intent(inout) :: vel(:,:)
    type(s_constraints), target, intent(inout) :: constraints

    ! local variables
    real(dp)                     :: imass1, imass2
    real(dp)                     :: dr12(3), dv12(3)
    real(dp)                     :: drdv, b12
    real(dp)                     :: tolerance
    integer                      :: i, j
    integer                      :: atm1, atm2
    integer                      :: nbond, iteration
    logical                      :: rattle_end

    real(wp),            pointer :: r0(:)
    integer,             pointer :: bond_list(:,:)


    bond_list => constraints%bond_list
    r0        => constraints%bond_dist

    iteration =  constraints%shake_iteration
    tolerance =  real(constraints%shake_tolerance,dp)
    nbond     =  constraints%num_bonds


    do i = 1, iteration

      rattle_end = .true.

      do j = 1, nbond
        atm1  = bond_list(1,j)
        atm2  = bond_list(2,j)

        dr12(1:3) = real(coord(1:3,atm1) - coord(1:3,atm2),dp)
        dv12(1:3) = real(vel  (1:3,atm1) - vel  (1:3,atm2),dp)

        drdv = dr12(1)*dv12(1) + dr12(2)*dv12(2) + dr12(3)*dv12(3)

        if (abs(drdv) >= tolerance) then

          imass1 = 1.0_dp/real(mass(atm1),dp)
          imass2 = 1.0_dp/real(mass(atm2),dp)

          b12 = drdv/((imass1 + imass2)*r0(j)*r0(j))

          vel(1:3,atm1) = vel(1:3,atm1) - real(b12*imass1*dr12(1:3),wp)
          vel(1:3,atm2) = vel(1:3,atm2) + real(b12*imass2*dr12(1:3),wp)

          rattle_end = .false.
        end if
      end do

      if (rattle_end) exit

    end do

    return

  end subroutine compute_rattle

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_settle_vv2
  !> @brief        SETTLE for rigid water (VV2)
  !! @authors      JJ, TM
  !! @param[inout] coord       : coordinates
  !! @param[inout] vel         : velocities
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_settle_vv2(coord, vel, constraints)

    ! formal arguments
    real(wp),                    intent(inout) :: coord(:,:)
    real(wp),                    intent(inout) :: vel(:,:)
    type(s_constraints), target, intent(inout) :: constraints

    ! local variables
    real(wp)                  :: massO, massH, massOH, massHH, mass2OH
    real(wp)                  :: massH2, massOH2
    real(wp)                  :: rab(1:3), rbc(1:3), rca(1:3)
    real(wp)                  :: vab(1:3), vbc(1:3), vca(1:3)
    real(wp)                  :: unit_ab(1:3), unit_bc(1:3), unit_ca(1:3)
    real(wp)                  :: length_ab, length_bc, length_ca
    real(wp)                  :: vab0, vbc0, vca0, cosA, cosB, cosC
    real(wp)                  :: cosA2, cosB2, cosC2
    real(wp)                  :: MaBC, MbCA, McAB, T_AB, T_BC, T_CA, Det
    integer                   :: num_water
    integer                   :: i, k, iatom_O, iatom_H1, iatom_H2
    integer                   :: id, nthread
    integer                   :: omp_get_num_threads, omp_get_thread_num
    integer                   :: quotient, remainder
    integer                   :: istart, iend

    integer,          pointer :: water_list(:,:)


    water_list => constraints%water_list

    num_water  = constraints%num_water
    massO      = constraints%water_massO
    massH      = constraints%water_massH

    massOH     = massO  + massH
    massHH     = massH  + massH
    mass2OH    = massOH + massOH
    massH2     = massH*massH
    massOH2    = massOH*massOH

    !$omp parallel default(none)                                               &
    !$omp private(i, k, rab, rbc, rca, vab, vbc, vca,                          &
    !$omp         length_ab, length_bc, length_ca, unit_ab, unit_bc, unit_ca,  &
    !$omp         vab0, vbc0, vca0, cosA, cosB, cosC, cosA2, cosB2, cosC2,     &
    !$omp         MaBC, MbCA, McAB, T_AB, T_BC, T_CA, Det, quotient, remainder,&
    !$omp         istart, iend, iatom_O, iatom_H1, iatom_H2, id)               &
    !$omp shared (water_list, coord, vel, num_water, massO, massH, massOH,     &
    !$omp         massHH, mass2OH, massH2, massOH2, nthread)
    !

#ifdef OMP
    nthread = omp_get_num_threads()
    id      = omp_get_thread_num()

    quotient = num_water / nthread
    remainder = mod(num_water, nthread)

    if (id <= (remainder-1)) then
      quotient = quotient + 1
      istart   = quotient*id + 1
      iend     = istart + quotient - 1
    else
      istart   = (quotient+1)*remainder + quotient*(id-remainder) + 1
      iend     = istart + quotient - 1
    end if
#else
    istart  = 1
    iend    = num_water
    nthread = 1
    id      = 0
#endif

    do i = istart, iend

      iatom_O  = water_list(1,i)
      iatom_H1 = water_list(2,i)
      iatom_H2 = water_list(3,i)

      rab(1:3) = coord(1:3,iatom_H1) - coord(1:3,iatom_O)
      rbc(1:3) = coord(1:3,iatom_H2) - coord(1:3,iatom_H1)
      rca(1:3) = coord(1:3,iatom_O)  - coord(1:3,iatom_H2)

      vab(1:3) = vel(1:3,iatom_H1) - vel(1:3,iatom_O)
      vbc(1:3) = vel(1:3,iatom_H2) - vel(1:3,iatom_H1)
      vca(1:3) = vel(1:3,iatom_O)  - vel(1:3,iatom_H2)

      length_ab = sqrt(rab(1)**2 + rab(2)**2 + rab(3)**2)
      length_bc = sqrt(rbc(1)**2 + rbc(2)**2 + rbc(3)**2)
      length_ca = sqrt(rca(1)**2 + rca(2)**2 + rca(3)**2)

      unit_ab(1:3) = rab(1:3)/length_ab
      unit_bc(1:3) = rbc(1:3)/length_bc
      unit_ca(1:3) = rca(1:3)/length_ca

      vab0 = 0.0_wp
      vbc0 = 0.0_wp
      vca0 = 0.0_wp

      do k = 1, 3
        vab0 = vab0 + vab(k)*unit_ab(k)
        vbc0 = vbc0 + vbc(k)*unit_bc(k)
        vca0 = vca0 + vca(k)*unit_ca(k)
      end do

      cosA = 0.0_wp
      cosB = 0.0_wp
      cosC = 0.0_wp

      do k = 1, 3
        cosA = cosA - unit_ab(k)*unit_ca(k)
        cosB = cosB - unit_bc(k)*unit_ab(k)
        cosC = cosC - unit_ca(k)*unit_bc(k)
      end do

      MbCA = massH*cosC*cosA - massOH*cosB
      MaBC = massO*cosB*cosC - massHH*cosA
      McAB = massH*cosA*cosB - massOH*cosC

      cosA2 = cosA*cosA
      cosB2 = cosB*cosB
      cosC2 = cosC*cosC

      T_AB = vab0*(mass2OH - massO *cosC2) + vbc0*MbCA + vca0*MaBC
      T_AB = T_AB*massO
      T_BC = vbc0*(massOH2 - massH2*cosA2)
      T_BC = T_BC + vca0*massO*McAB + vab0*massO*MbCA
      T_CA = vca0*(mass2OH - massO *cosB2) + vab0*MaBC + vbc0*McAB
      T_CA = T_CA*massO

      Det  = 2.0_wp*massOH2 + 2.0_wp*massO*massH*cosA*cosB*cosC
      Det  = Det - 2.0_wp*massH2*cosA2 - massO*massOH*(cosB2+cosC2)
      Det  = Det / massH

      T_AB = T_AB / Det
      T_BC = T_BC / Det
      T_CA = T_CA / Det

      vel(1:3,iatom_O)  =  vel(1:3,iatom_O)  &
                         + (T_AB*unit_ab(1:3) - T_CA*unit_ca(1:3))/massO
      vel(1:3,iatom_H1) = vel(1:3,iatom_H1)  &
                         + (T_BC*unit_bc(1:3) - T_AB*unit_ab(1:3))/massH
      vel(1:3,iatom_H2) = vel(1:3,iatom_H2)  &
                         + (T_CA*unit_ca(1:3) - T_BC*unit_bc(1:3))/massH

    end do
    !$omp end parallel

    return

  end subroutine compute_settle_vv2

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    clear_fixatm_component
  !> @brief        clear to zero the specified component
  !! @authors      KY
  !! @param[in]    constraints : boundary conditions information
  !! @param[in]    natom       : number of atoms
  !! @param[inout] comp        : velocity or force
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine clear_fixatm_component(constraints, natom, comp)

    ! formal arguments
    type(s_constraints), intent(in)    :: constraints
    integer,             intent(in)    :: natom
    real(wp),            intent(inout) :: comp(3,natom)

    ! local variable
    integer  :: i


    if (constraints%num_fixatm == 0) return

    !$omp parallel 
    !$omp do 
    do i = 1, natom
      if (constraints%fixatm(i)) comp(:,i) = 0.0_wp
    end do
    !$omp end do
    !$omp end parallel

    return

  end subroutine clear_fixatm_component

end module at_constraints_mod
