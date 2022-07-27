!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_boundary_mod
!> @brief   utilities for boundary conditions
!! @authors Takaharu Mori (TM), Takashi Imai (TI), Jaewoon Jung (JJ), 
!!          Norio Takase (NT), Motoshi Kamiya (MK), Kiyoshi Yagi (KY),
!!          Cheng Tan (CT), Yuji Sugita (YS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_boundary_mod

  use at_boundary_str_mod
  use molecules_str_mod
  use molecules_mod
  use fileio_rst_mod
  use fileio_spot_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use select_mod
  use select_atoms_mod
  use select_atoms_str_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  real(wp), parameter   :: TableWaterMargin = 2.0_wp
  real(wp), parameter   :: FakeDefault      = 999999_wp

  ! structures
  type, public :: s_pbc_info
    real(wp)            :: box_size_x       = FakeDefault
    real(wp)            :: box_size_y       = FakeDefault
    real(wp)            :: box_size_z       = FakeDefault
    real(wp)            :: pairlist_grid    = 3.0_wp
    logical             :: wrap_all         = .false.
    logical             :: calc_local_pbc   = .false.
  end type s_pbc_info

  type, public :: s_spot_info
    ! spherical potential
    logical                         :: spherical_pot  = .false.
    real(wp)                        :: const          = 10.0_wp
    integer                         :: exponent       = 2
    logical                         :: read_spot      = .false.
    real(wp)                        :: mod_ratio      = 1.0_wp
    integer                         :: nindex         = 0
    integer                         :: nfunctions     = 0
    character(MaxLine), allocatable :: center(:) 
    real(wp), allocatable           :: radius(:)
    logical                         :: fixatom        = .true.
    real(wp)                        :: fix_layer      = 1.0_wp
    character(MaxLine)              :: nospot_select_index = ''
    logical                         :: restart        = .true.
  end type s_spot_info

  type, public :: s_boundary_info
    integer             :: type             = BoundaryTypePBC
    type(s_pbc_info)    :: pbc_info
    type(s_spot_info)   :: spot_info
    real(wp)            :: origin_x         = 0.0_wp
    real(wp)            :: origin_y         = 0.0_wp
    real(wp)            :: origin_z         = 0.0_wp
    logical             :: shift_origin     = .true.
  end type s_boundary_info

  ! subroutines
  public  :: show_ctrl_boundary
  public  :: read_ctrl_boundary
  public  :: setup_boundary
  private :: setup_spot
  private :: setup_nospotlist
  private :: setup_fixatm_boundary
  public  :: update_boundary
  private :: update_boundary_pbc
  public  :: update_boundary_cg
  private :: update_boundary_pbc_cg
  public  :: update_spot_atomlist
  public  :: wrap_molecules
  private :: wrap_all

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_boundary
  !> @brief        show BOUNDARY section usage
  !! @authors      NT, TM, KY
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "md", "min", "remd", "rpath"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_boundary(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode


    if (show_all) then

      select case (run_mode)

      case ('md', 'min', 'vib', 'remd', 'rpath', 'bd')

        write(MsgOut,'(A)') '[BOUNDARY]'
        write(MsgOut,'(A)') 'type          = PBC       # [NOBC,PBC]'
        write(MsgOut,'(A)') 'box_size_x    = 999999    # box size (x)     in [PBC]'
        write(MsgOut,'(A)') 'box_size_y    = 999999    # box size (y)     in [PBC]'
        write(MsgOut,'(A)') 'box_size_z    = 999999    # box size (z)     in [PBC]'
        write(MsgOut,'(A)') 'spherical_pot = YES       # use sph. pot.    in [NOBC]'
        write(MsgOut,'(A)') 'constant      = 10.0      # force constants  in [NOBC]'
        write(MsgOut,'(A)') 'exponent      = 2         # exponent         in [NOBC]'
        write(MsgOut,'(A)') 'fixatom       = yes       # fixed atom at bounary  in [NOBC]'
        write(MsgOut,'(A)') 'fix_layer     = 1.0       # thickness of fix layer in [NOBC]'
        !write(MsgOut,'(A)') 'nospot_select_index = 1   # index for atoms without sph. pot.'
        write(MsgOut,'(A)') 'restart       = yes       # retrieve from rst in [NOBC]'
        write(MsgOut,'(A)') '# input center using spotfile'
        write(MsgOut,'(A)') '#modify_ratio  = 1.0       # change radius in spotfile    '
        write(MsgOut,'(A)') '# input center using selector'
        write(MsgOut,'(A)') '#nindex        = 1        # number of functions in [NOBC]'
        write(MsgOut,'(A)') '#center_select_index1 = 1 # selector index for center in [NOBC]'
        write(MsgOut,'(A)') '#radius1       = 20.0     # radius              in [NOBC]'
        write(MsgOut,'(A)') '# input center directly'
        write(MsgOut,'(A)') '#nfunctions    = 1             # number of functions in [NOBC]'
        write(MsgOut,'(A)') '#center1       = 0.0, 0.0, 0.0 # center x,y,z        in [NOBC]'
        write(MsgOut,'(A)') '#radius1       = 20.0          # radius              in [NOBC]'
        write(MsgOut,'(A)') ' '

      end select

    else

      select case (run_mode)

      case ('md', 'min', 'vib', 'remd', 'rpath', 'bd')

        write(MsgOut,'(A)') '[BOUNDARY]'
        write(MsgOut,'(A)') 'type          = PBC       # [NOBC,PBC]'
        write(MsgOut,'(A)') 'box_size_x    = 999999    # box size (x) in [PBC]'
        write(MsgOut,'(A)') 'box_size_y    = 999999    # box size (y) in [PBC]'
        write(MsgOut,'(A)') 'box_size_z    = 999999    # box size (z) in [PBC]'
        write(MsgOut,'(A)') 'spherical_pot = YES       # use sph. pot.    in [NOBC]'
        write(MsgOut,'(A)') 'constant      = 10.0      # force constants  in [NOBC]'
        write(MsgOut,'(A)') 'exponent      = 2         # exponent         in [NOBC]'
        write(MsgOut,'(A)') ' '

      end select

    end if


    return

  end subroutine show_ctrl_boundary
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_boundary
  !> @brief        read BOUNDARY section in the control file
  !! @authors      YS, TI, JJ, TM, KY
  !! @param[in]    handle     : unit number
  !! @param[out]   bound_info : BOUNDARY section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_boundary(handle, bound_info)

    ! parameters
    character(*),            parameter     :: Section = 'Boundary'

    ! formal arguments
    integer,                        intent(in)    :: handle
    type(s_boundary_info), target,  intent(inout) :: bound_info

    ! local variables
    integer                    :: i, nfunc, ierr
    character(20)              :: key1, key2
    type(s_spot_info), pointer :: spot_info


    ! use pointers
    !
    spot_info => bound_info%spot_info

    ! read parameters from control file
    ! 
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_type(handle, Section, 'type',          &
                            bound_info%type, BoundaryTypeTypes)

    call read_ctrlfile_logical(handle, Section, 'shift',      &
                               bound_info%shift_origin)

    select case (bound_info%type)

    case (BoundaryTypePBC)

      ! read parameters for PBC from control file
      !
      call read_ctrlfile_real(handle, Section, 'box_size_x',  &
                              bound_info%pbc_info%box_size_x)
      call read_ctrlfile_real(handle, Section, 'box_size_y',  &
                              bound_info%pbc_info%box_size_y)
      call read_ctrlfile_real(handle, Section, 'box_size_z',  &
                              bound_info%pbc_info%box_size_z)
      call read_ctrlfile_logical(handle, Section, 'wrapall',  &
                              bound_info%pbc_info%wrap_all)
      call read_ctrlfile_real(handle, Section, 'pairlist_grid',  &
                              bound_info%pbc_info%pairlist_grid)
      call read_ctrlfile_logical(handle, Section, 'local_pbc',  &
          bound_info%pbc_info%calc_local_pbc)
      
    case (BoundaryTypeNOBC)

      ! read parameters for spherical potential
      !
      spot_info%nindex     = 0
      spot_info%nfunctions = 0

      call read_ctrlfile_logical(handle, Section, 'spherical_pot',  &
                                 spot_info%spherical_pot)
      call read_ctrlfile_integer(handle, Section, 'nindex',       &
                                 spot_info%nindex)
      call read_ctrlfile_integer(handle, Section, 'nfunctions',   &
                                 spot_info%nfunctions)
      call read_ctrlfile_real   (handle, Section, 'constant',     &
                                 spot_info%const)
      call read_ctrlfile_integer(handle, Section, 'exponent',     &
                                 spot_info%exponent)
      call read_ctrlfile_real   (handle, Section, 'modify_ratio', &
                                 spot_info%mod_ratio)
      call read_ctrlfile_logical(handle, Section, 'fixatom',      &
                                 spot_info%fixatom)
      call read_ctrlfile_real   (handle, Section, 'fix_layer',    &
                                 spot_info%fix_layer)
      call read_ctrlfile_string (handle, Section, 'nospot_select_index', &
                                 spot_info%nospot_select_index)
      call read_ctrlfile_logical(handle, Section, 'restart',      &
                                 spot_info%restart)

      if (spot_info%spherical_pot) then
        if (spot_info%nindex == 0 .and. spot_info%nfunctions == 0)  &
          spot_info%read_spot = .true.

        if (spot_info%nindex > 0 .and. spot_info%nfunctions > 0)  &
          call error_msg('Read_Ctrl_Boundary> ERROR: Both nindex and &
                         &nfunctions are given in ctrl file.')

        nfunc = 0
        ierr  = 0
        if (spot_info%nindex > 0) then

          nfunc = spot_info%nindex

          allocate(spot_info%center(nfunc), &
                   spot_info%radius(nfunc))

          spot_info%center(1:nfunc) = ''
          spot_info%radius(1:nfunc) = -1.0_wp
          do i = 1, nfunc

            write(key1,'(a,i0)') 'center_select_index',i
            write(key2,'(a,i0)') 'radius',i
            call read_ctrlfile_string (handle, Section, key1, &
                                       spot_info%center(i))
            call read_ctrlfile_real   (handle, Section, key2, &
                                       spot_info%radius(i))

            if (spot_info%center(i) .eq. '') then
              ierr = -1
              if (main_rank) &
                write(MsgOut,'(a)') &
                  'Read_Ctrl_Boundary> ERROR: '//trim(key1)// &
                  ' is not specified in ctrl'
            end if

            if (spot_info%radius(i) < 0.0_wp) then
              ierr = -1
              if (main_rank) &
                write(MsgOut,'(a)') &
                  'Read_Ctrl_Boundary> ERROR: '//trim(key2)// &
                  ' is not specified in ctrl'
            end if

          end do

        else if (spot_info%nfunctions > 0) then

          nfunc = spot_info%nfunctions

          allocate(spot_info%center(nfunc), &
                   spot_info%radius(nfunc))

          spot_info%center(1:nfunc) = ''
          spot_info%radius(1:nfunc) = -1.0_wp
          do i = 1, nfunc

            write(key1,'(a,i0)') 'center',i
            write(key2,'(a,i0)') 'radius',i

            call read_ctrlfile_string(handle, Section, key1, &
                                      spot_info%center(i))
            call read_ctrlfile_real  (handle, Section, key2, &
                                      spot_info%radius(i))

            if (spot_info%center(i) .eq. '') then
              ierr = -1
              if(main_rank) &
                write(MsgOut,'(a)') &
                  'Read_Ctrl_Boundary> ERROR: '//trim(key1)// &
                  ' is not specified in ctrl'
            end if

            if (spot_info%radius(i) < 0.0_wp) then
              ierr = -1
              if(main_rank) &
                write(MsgOut,'(a)') &
                  'Read_Ctrl_Boundary> ERROR: '//trim(key2)// &
                  ' is not specified in ctrl'
            end if

          end do

        end if
        if (ierr == -1) &
          call error_msg('Read_Ctrl_Boundary> stop with error.')

      end if
      
    end select

    call end_ctrlfile_section(handle)


    ! write parameters to MsgOut
    !
    if (main_rank) then

      write(MsgOut,'(A)') 'Read_Ctrl_Boundary> Parameters of Boundary Condition'
      write(MsgOut,'(A20,A10)') &
            '  type            = ', BoundaryTypeTypes(bound_info%type)
      if (.not. bound_info%shift_origin) &
        write(MsgOut,'(A)') '  shift_origin    = no'

      select case (bound_info%type)

      case (BoundaryTypePBC)

        if (abs(bound_info%pbc_info%box_size_x - FakeDefault) > EPS .and. &
            abs(bound_info%pbc_info%box_size_y - FakeDefault) > EPS .and. &
            abs(bound_info%pbc_info%box_size_z - FakeDefault) > EPS) &
        write(MsgOut,'(A20,3F10.3)')                                  &
              '  box_size(x,y,z) = ', bound_info%pbc_info%box_size_x, &
                                      bound_info%pbc_info%box_size_y, &
                                      bound_info%pbc_info%box_size_z
        write(MsgOut,'(A20,F10.3)')                                   &
              '  pairlist_grid   = ', bound_info%pbc_info%pairlist_grid
        if (bound_info%pbc_info%calc_local_pbc) &
            write(MsgOut,'(A)') '  local_pbc       = yes'

      case (BoundaryTypeNOBC)

        if (spot_info%spherical_pot) then
          write(MsgOut,'(A30,A20,L10)') &
            '  spherical_pot   =        yes',  &
            '  restart         = ', spot_info%restart

          if (spot_info%read_spot) then
            write(MsgOut,'(A30,A20,F10.3)')     &
            '  center_info     =       read',   &
            '  modify_ratio    = ', spot_info%mod_ratio

          else if (spot_info%nindex /= 0) then
            write(MsgOut,'(A30,A20,I10)')     &
            '  center_info     =   selector',   &
            '  nindex          = ', spot_info%nindex

          else if (spot_info%nfunctions /= 0) then
            write(MsgOut,'(A30,A20,I10)')     &
            '  center_info     =    control',   &
            '  nfunctions      = ', spot_info%nfunctions

          end if

          write(MsgOut,'(A20,F10.3,A20,I10)')          &
            '  constant        = ', spot_info%const,   &
            '  exponent        = ', spot_info%exponent
          write(MsgOut,'(A20,L10,A20,F10.3)')          &
            '  fixatom         = ', spot_info%fixatom, &
            '  fix_layer       = ', spot_info%fix_layer
          !if (spot_info%nospot_select_index == '') then
          !  write(MsgOut,'(A)')        &
          !  '  nospot_select_index =   none'
          !else
          !  write(MsgOut,'(A)')        &
          !  '  nospot_select_index = '//trim(spot_info%nospot_select_index)
          !end if

        else
          write(MsgOut,'(A30,A20,L10)') &
            '  spherical_pot   =         no',  &
            '  restart         = ', spot_info%restart

        end if

      end select

      write(MsgOut,'(A)') ' '

    end if

    return

  end subroutine read_ctrl_boundary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_boundary
  !> @brief        set essential variables for boundary condition
  !! @authors      YS, TM, JJ, TI, KY, CT
  !! @param[in]    bound_info  : BOUNDARY section control parameters information
  !! @param[in]    use_table   : flag for use table or not
  !! @param[in]    pairlistdist: pair-list distance
  !! @param[in]    sel_info    : SELECTOR section in control parameters
  !! @param[in]    molecule    : molecule
  !! @param[in]    rst         : restart file information
  !! @param[in]    spot        : spherical potential information
  !! @param[out]   boundary    : boundary conditions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_boundary(bound_info, use_table, pairlistdist, sel_info, &
                            molecule, rst, spot, boundary)

    ! formal arguments
    type(s_boundary_info),   intent(in)    :: bound_info
    logical,                 intent(in)    :: use_table
    real(wp),                intent(in)    :: pairlistdist
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_molecule),        intent(inout) :: molecule
    type(s_rst),             intent(in)    :: rst
    type(s_spot),            intent(in)    :: spot
    type(s_boundary),        intent(inout) :: boundary

    ! local variables
    integer    :: nf, i, iatom


    ! initialize 
    !
    call init_boundary(boundary)

    ! allocate fixatm
    allocate(boundary%fixatm(molecule%num_atoms))
    boundary%fixatm = .false.

    ! setup variables
    !
    select case (bound_info%type)
    

    case (BoundaryTypePBC)

      boundary%type     = bound_info%type
      boundary%origin_x = bound_info%origin_x
      boundary%origin_y = bound_info%origin_y
      boundary%origin_z = bound_info%origin_z
      boundary%wrap_all = bound_info%pbc_info%wrap_all
      boundary%calc_local_pbc = bound_info%pbc_info%calc_local_pbc

      if (rst%rstfile_type == RstfileTypeUndef) then
        if (abs(bound_info%pbc_info%box_size_x - FakeDefault) < EPS) &
          call error_msg('Setup_Boundary> box_size_x is not specified in ctrl')
        if (abs(bound_info%pbc_info%box_size_y - FakeDefault) < EPS) &
          call error_msg('Setup_Boundary> box_size_y is not specified in ctrl')
        if (abs(bound_info%pbc_info%box_size_z - FakeDefault) < EPS) &
          call error_msg('Setup_Boundary> box_size_z is not specified in ctrl')

        boundary%box_size_x = bound_info%pbc_info%box_size_x
        boundary%box_size_y = bound_info%pbc_info%box_size_y
        boundary%box_size_z = bound_info%pbc_info%box_size_z
      else
        boundary%box_size_x = rst%box_size_x
        boundary%box_size_y = rst%box_size_y
        boundary%box_size_z = rst%box_size_z
      end if
      boundary%box_size_x_ref = boundary%box_size_x
      boundary%box_size_y_ref = boundary%box_size_y
      boundary%box_size_z_ref = boundary%box_size_z

      boundary%pairlist_grid = bound_info%pbc_info%pairlist_grid

    case (BoundaryTypeNOBC)

      boundary%type     = bound_info%type
      boundary%origin_x = bound_info%origin_x
      boundary%origin_y = bound_info%origin_y
      boundary%origin_z = bound_info%origin_z
      boundary%sph_pot  = bound_info%spot_info%spherical_pot

      if (boundary%sph_pot) then
        if (rst%rstfile_type == RstfileTypeUndef .or. &
            .not. bound_info%spot_info%restart) then

            call setup_spot(bound_info, sel_info, molecule, spot, boundary)

            boundary%fix_layer = bound_info%spot_info%fix_layer
            if (bound_info%spot_info%fixatom) &
              call setup_fixatm_boundary(molecule, boundary)

            !if (bound_info%spot_info%nospot_select_index /= '') &
            !  call setup_nospotlist(bound_info%spot_info,       &
            !                        sel_info, molecule, boundary)
            !
            !allocate(boundary%atomlist(molecule%num_atoms))
            !if (allocated(rst%coord)) then
            !  call update_spot_atomlist(molecule%num_atoms, & 
            !                            rst%coord, boundary)
            !else
            !  call update_spot_atomlist(molecule%num_atoms, & 
            !                            molecule%atom_coord, boundary)
            !end if

        else
          if (.not. rst%sph_pot) &
            call error_msg('Setup_Boundary> Information of spherical potential &
                           &is not found in rstfile.')
            
          call alloc_boundary(boundary, BoundarySphericalPot, rst%nfunctions)
          boundary%nfunctions = rst%nfunctions
          boundary%radius     = rst%radius    
          boundary%center     = rst%center    
          boundary%const      = rst%const     
          boundary%exponent   = rst%exponent  
          boundary%fix_layer  = rst%fix_layer 
          boundary%num_fixatm = rst%num_fixatm
          boundary%fixatm     = rst%fixatm    

        end if
      end if

    end select

    ! write setup info
    !
    if (main_rank) then
      write(MsgOut,'(A)') &
            'Setup_Boundary> Setup Variables for Boundary Condition'

      if (boundary%type == BoundaryTypePBC) then
        write(MsgOut,'(A20,3F10.3)')                       &
              '  box_size(x,y,z) = ', boundary%box_size_x, &
                                      boundary%box_size_y, &
                                      boundary%box_size_z
      end if

      write(MsgOut,'(A20,3F10.3)')                         &
              '  origin(x,y,z)   = ', boundary%origin_x,   &
                                      boundary%origin_y,   &
                                      boundary%origin_z

      if (boundary%type == BoundaryTypeNOBC) then

        if (boundary%sph_pot) then
          write(MsgOut,'(A20,I8)') &
              '  nfunctions      = ', boundary%nfunctions

          do nf = 1, boundary%nfunctions
            write(MsgOut,'(A20,3F10.3)')                           &
                    '    center(x,y,z) = ', boundary%center(1,nf),  &
                                            boundary%center(2,nf),  &
                                            boundary%center(3,nf)
            write(MsgOut,'(A20,F10.3)') &
              '    radius        = ',boundary%radius(nf)
            write(MsgOut,'(A20,F10.3)') &
              '    constant      = ',boundary%const(nf)
            write(MsgOut,'(A20,I10)') &
              '    exponent      = ',boundary%exponent(nf)
          end do
          write(MsgOut,*)

          write(MsgOut,'(A20,I8)') &
              '  num_fixatm      = ', boundary%num_fixatm
          write(MsgOut,*)

          if (boundary%num_nospot > 0) then
            write(MsgOut,'(A20,I8)') &
              '  num_nospot      = ', boundary%num_nospot
            if (boundary%num_nospot < 100) then
              do i = 1, boundary%num_nospot
                iatom = boundary%nospotlist(i)
                write(MsgOut,'(2x,i6,x,a4,x,i6,x,a4,x,a6,x,a6)') &
                  iatom,                         &
                  molecule%segment_name(iatom),  &
                  molecule%residue_no(iatom),    &
                  molecule%residue_name(iatom),  &
                  molecule%atom_name(iatom),     &
                  molecule%atom_cls_name(iatom)
              end do
            end if
          end if
        end if

      end if
      write(MsgOut,'(A)') ' '

    end if

    return

  end subroutine setup_boundary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_spot
  !> @brief        setup variables for spherical potential
  !! @authors      KY
  !! @param[in]    bound_info  : BOUNDARY section control parameters information
  !! @param[in]    sel_info    : SELECTOR section in control parameters
  !! @param[in]    molecule    : molecule
  !! @param[in]    spot        : spherical potential information
  !! @param[out]   boundary    : boundary conditions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_spot(bound_info, sel_info, molecule, spot, boundary)

    ! formal arguments
    type(s_boundary_info),   intent(in)    :: bound_info
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_molecule),        intent(in)    :: molecule
    type(s_spot),            intent(in)    :: spot
    type(s_boundary),        intent(inout) :: boundary

    ! local variables
    integer    :: ng, nf, nfunc
    integer    :: nindex
    integer    :: igroup, max_ngroup, i, j, offset
    integer, allocatable          :: ngroup(:)
    integer, allocatable          :: group_list(:)
    type(s_selatoms), allocatable :: selatoms(:,:)


    if (spot%ncenter > 0) then
      if (bound_info%spot_info%nindex > 0)     &
        call error_msg('Setup_Boundary> ERROR: Both spotfile and &
                       &nindex are given in ctrl file.')
      if (bound_info%spot_info%nfunctions > 0) &
        call error_msg('Setup_Boundary> ERROR: Both spotfile and &
                       &nfunctions are given in ctrl file.')
    end if

    if (bound_info%spot_info%read_spot) then
      ! read from spotfile

      if (spot%ncenter == 0) &
        call error_msg('Setup_Boundary> ERROR: Center of the spherical &
                       &potential must be specified by one of spotfile, nindex, nfunctions')
    
      nfunc = spot%ncenter
      boundary%nfunctions  = nfunc
      call alloc_boundary(boundary, BoundarySphericalPot, nfunc)

      boundary%const   = bound_info%spot_info%const
      boundary%exponent= bound_info%spot_info%exponent

      boundary%center  = spot%center
      boundary%radius  = spot%radius*bound_info%spot_info%mod_ratio

    else if (bound_info%spot_info%nindex > 0) then
      ! set from selector

      nindex = bound_info%spot_info%nindex

      allocate(ngroup(nindex))
      max_ngroup = 0
      do i = 1, nindex
        ngroup(i) = split_num(trim(bound_info%spot_info%center(i)))
        if(ngroup(i) > max_ngroup) max_ngroup = ngroup(i)
      end do
      allocate(selatoms(max_ngroup, nindex))

      nfunc  = 0
      do i = 1, nindex
        if (ngroup(i) /= 0) then
          allocate(group_list(ngroup(i)))
          call split(ngroup(i), ngroup(i), &
                     bound_info%spot_info%center(i), group_list)

          do ng = 1, ngroup(i)
            igroup = group_list(ng)
            call select_atom(molecule, sel_info%groups(igroup), selatoms(ng, i))
            nfunc = nfunc + size(selatoms(ng, i)%idx)
          end do

          deallocate(group_list)

        end if

      end do

      boundary%nfunctions = nfunc
      call alloc_boundary(boundary, BoundarySphericalPot, nfunc)

      boundary%const    = bound_info%spot_info%const
      boundary%exponent = bound_info%spot_info%exponent

      offset = 0
      do i = 1, nindex
        do ng = 1, ngroup(i)
          nfunc = size(selatoms(ng,i)%idx)
          boundary%radius(offset+1:offset+nfunc) = bound_info%spot_info%radius(i)

          do nf = 1, nfunc
            j = selatoms(ng,i)%idx(nf)
            boundary%center(:,offset+nf) = molecule%atom_coord(:,j)

            !write(MsgOut,'(i6,x,a4,x,i6,x,a4,x,a6,x,a6)') &
            !  j, &
            !  molecule%segment_name(j), &
            !  molecule%residue_no(j),   &
            !  molecule%residue_name(j), &
            !  molecule%atom_name(j),    &
            !  molecule%atom_cls_name(j)
            !write(MsgOut,'(3f12.4)') molecule%atom_coord(:,j)

          end do
          offset = offset + nfunc
        end do
      end do

      deallocate(selatoms, ngroup)

    else if(bound_info%spot_info%nfunctions > 0) then
      ! set from a control file

      nfunc = bound_info%spot_info%nfunctions
      boundary%nfunctions  = nfunc
      call alloc_boundary(boundary, BoundarySphericalPot, nfunc)

      boundary%const   = bound_info%spot_info%const
      boundary%exponent= bound_info%spot_info%exponent

      do nf = 1, nfunc
        call split(3,3,bound_info%spot_info%center(nf),boundary%center(:,nf))
      end do
      boundary%radius  = bound_info%spot_info%radius

    end if

    return

  end subroutine setup_spot

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_nospotlist
  !> @brief        setup list of atoms without the spherical potential
  !! @authors      KY
  !! @param[in]    spot_info : BOUNDARY section control parameters information
  !! @param[in]    sel_info  : SELECTOR section in control parameters
  !! @param[in]    molecule  : molecule
  !! @param[out]   boundary  : boundary conditions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_nospotlist(spot_info, sel_info, molecule, boundary)

    ! formal arguments
    type(s_spot_info),   intent(in)    :: spot_info
    type(s_sel_info),    intent(in)    :: sel_info
    type(s_molecule),    intent(in)    :: molecule
    type(s_boundary),    intent(inout) :: boundary

    ! local variables
    integer                       :: igroup, ngroup, natom, i, j, offset, temp
    integer,          allocatable :: group_list(:)
    type(s_selatoms), allocatable :: selatoms(:)
    integer                       :: kalloc_stat, kdealloc_stat


    ngroup = split_num(trim(spot_info%nospot_select_index))
    if (ngroup /= 0) then

      ! Number of atoms
      !
      allocate(group_list(ngroup))
      call split(ngroup, ngroup, spot_info%nospot_select_index, group_list)

      allocate(selatoms(ngroup))
      natom = 0
      do i = 1, ngroup
        igroup = group_list(i)
        call select_atom(molecule, sel_info%groups(igroup), selatoms(i))
        natom = natom + size(selatoms(i)%idx)
      end do

      boundary%num_nospot = natom

      ! List of atoms
      !
      allocate(boundary%nospotlist(natom), &
               stat = kalloc_stat)
      if(kalloc_stat /= 0) call error_msg_alloc

      offset = 0
      do i = 1, ngroup
        igroup = group_list(i)
        natom = size(selatoms(i)%idx)
        boundary%nospotlist(offset+1:offset+natom) = selatoms(i)%idx(1:natom)
        offset = offset + natom
      end do

      deallocate(selatoms)
      deallocate(group_list)

      ! Sort fixed atom indices in ascending order
      !
      do i = boundary%num_nospot, 2, -1
        do j = 1, i - 1
          if (boundary%nospotlist(j) > boundary%nospotlist(j+1)) then
            temp = boundary%nospotlist(j)
            boundary%nospotlist(j)   = boundary%nospotlist(j+1)
            boundary%nospotlist(j+1) = temp
          end if
        end do
      end do

    end if

    return

  end subroutine setup_nospotlist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_fixatm_boundary
  !> @brief        setup atoms to be fixed
  !! @authors      KY
  !! @param[in]    molecule : molecule information
  !! @param[inout] boundary : boundary conditions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_fixatm_boundary(molecule, boundary)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_boundary),         intent(inout) :: boundary

    ! local variables
    real(wp)     :: rr, dj(3), rj
    integer      :: natom, nres, nfunc, nfix
    integer      :: i, j, n, ierr
    logical      :: isFixed

    integer, allocatable :: jstart(:), jend(:)
    real(wp), pointer    :: coord(:,:)
 

    natom = molecule%num_atoms
    nres  = molecule%num_residues
    nfunc = boundary%nfunctions
    coord => molecule%atom_coord

    if (main_rank) then
      allocate(jstart(nres), jend(nres))
      jstart(1) = 1
      do i = 2, nres
        jstart(i) = jstart(i-1)
        do while (molecule%residue_c_no(jstart(i)) == i-1) 
          jstart(i) = jstart(i) + 1
        end do
      end do
      do i = 1, nres-1
        jend(i) = jstart(i+1) - 1
      end do
      jend(nres) = natom

      !$omp parallel &
      !$omp private(i, j, n, rr, dj, rj, isFixed)
      !$omp do
      do i = 1, nres
        do j = jstart(i), jend(i)

          isFixed = .true.
          do n = 1, nfunc
            rr  = boundary%radius(n) - boundary%fix_layer

            ! compute distance
            !
            dj(1:3) = coord(1:3,j) - boundary%center(1:3,n)
            rj      = sqrt(dj(1)*dj(1) + dj(2)*dj(2) + dj(3)*dj(3))

            if (rj < rr) then
              isFixed = .false.
              exit
            end if

          end do
          if (isFixed) exit

        end do

        if (isFixed) then
          do j = jstart(i), jend(i)
            boundary%fixatm(j) = .true.
          end do
        end if
      end do
      !$omp end do
      !$omp end parallel

      !dbg do i = 1, nres
      !dbg   if (boundary%fixatm(jstart(i))) &
      !dbg     write(MsgOut,'(i8,$)') i-1
      !dbg end do

      deallocate(jstart, jend)

      nfix = 0
      !$omp parallel 
      !$omp do reduction(+:nfix)
      do i = 1, natom
        if (boundary%fixatm(i)) &
          nfix = nfix + 1
      end do
      !$omp end do
      !$omp end parallel

      boundary%num_fixatm = nfix

    end if

#ifdef HAVE_MPI_GENESIS
    ! broadcast fixatom information
    call mpi_bcast(boundary%num_fixatm, 1, mpi_integer, 0, mpi_comm_world, ierr)
    call mpi_bcast(boundary%fixatm, natom, mpi_logical, 0, mpi_comm_world, ierr)
#endif

    return

  end subroutine setup_fixatm_boundary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_boundary
  !> @brief        update variables for boundary condition
  !! @authors      TM
  !! @param[in]    use_table    : flag for use table or not
  !! @param[in]    pairlistdist : pair list distance
  !! @param[inout] boundary     : information of boundary condition
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_boundary(use_table, pairlistdist, boundary)

    ! formal arguments
    logical,                 intent(in)    :: use_table
    real(wp),                intent(in)    :: pairlistdist
    type(s_boundary),        intent(inout) :: boundary


    select case (boundary%type)

    case (BoundaryTypeNOBC)

      ! do nothing

    case (BoundaryTypePBC)

      call update_boundary_pbc(use_table, pairlistdist, boundary)

    end select

    return

  end subroutine update_boundary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_boundary_pbc
  !> @brief        update variables for boundary condition
  !! @authors      TI, JJ, TM, MK
  !! @param[in]    use_table    : flag for use table or not
  !! @param[in]    pairlistdist : pair list distance
  !! @param[inout] boundary     : information of boundary condition
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_boundary_pbc(use_table, pairlistdist, boundary)

    ! formal arguments
    logical,                 intent(in)    :: use_table
    real(wp),                intent(in)    :: pairlistdist
    type(s_boundary),        intent(inout) :: boundary

    ! local variables
    real(wp)                 :: bsize_x, bsize_y, bsize_z
    real(wp)                 :: csize_x, csize_y, csize_z
    integer                  :: i, j, k
    integer                  :: inb, jnb, knb, inbs, jnbs, knbs, lc, lcnb, ln
    integer                  :: ncell_x, ncell_y, ncell_z, ncell
    integer                  :: ncell_x_check, ncell_y_check, ncell_z_check

    ! for small cell division
    real(wp)                 :: pairdist_w_margin, pairdist_w_margin2
    real(wp)                 :: pairdist_check
    real(wp)                 :: ci, cj, ck
    integer                  :: cellpair_max_x, cellpair_max_y, cellpair_max_z
    integer                  :: ik, jk, kk


    boundary%use_cell_linked_list = .true.
    bsize_x = boundary%box_size_x
    bsize_y = boundary%box_size_y
    bsize_z = boundary%box_size_z

    if (use_table) then

      ncell_x = int(bsize_x/boundary%pairlist_grid)
      ncell_y = int(bsize_y/boundary%pairlist_grid)
      ncell_z = int(bsize_z/boundary%pairlist_grid)

      pairdist_check = max(boundary%pairlist_grid,pairlistdist)

      ncell_x_check = int(bsize_x/pairdist_check)
      ncell_y_check = int(bsize_y/pairdist_check)
      ncell_z_check = int(bsize_z/pairdist_check)

      if (ncell_x_check < 3 .or. ncell_y_check < 3 .or. ncell_z_check < 3) &
        call error_msg('Update_Boundary_Pbc> too small boxsize/pairdist.'//&
                       ' larger boxsize or shorter pairdist should be used (see "Chapter: Trouble shooting" in the user manual).')

    else

      ncell_x = int(bsize_x/pairlistdist)
      ncell_y = int(bsize_y/pairlistdist)
      ncell_z = int(bsize_z/pairlistdist)

      if (ncell_x < 3 .or. ncell_y < 3 .or. ncell_z < 3) then
        if (pairlistdist < 0.5_wp*bsize_x .and. &
            pairlistdist < 0.5_wp*bsize_y .and. &
            pairlistdist < 0.5_wp*bsize_z) then
          boundary%use_cell_linked_list = .false.
          return
        else
          call error_msg('Update_Boundary_Pbc> too small boxsize/pairdist.'//&
                         ' larger boxsize or shorter pairdist should be used (see "Chapter: Trouble shooting" in the user manual).')
        end if
      end if

    end if

    csize_x = bsize_x/real(ncell_x,wp)
    csize_y = bsize_y/real(ncell_y,wp)
    csize_z = bsize_z/real(ncell_z,wp)
    ncell   = ncell_x*ncell_y*ncell_z

    ! MK (cellpairs; TableWaterMargin for water pairs)
    pairdist_w_margin = pairlistdist + TableWaterMargin
    pairdist_w_margin2 = pairdist_w_margin * pairdist_w_margin
    cellpair_max_x = pairdist_w_margin / csize_x + 1
    cellpair_max_y = pairdist_w_margin / csize_y + 1
    cellpair_max_z = pairdist_w_margin / csize_z + 1
    boundary%num_neighbor_cells = 0

    do k = -cellpair_max_z, cellpair_max_z

      ck = csize_z * real(max(0, abs(k) - 1), wp)
      ck = ck * ck

      do j = -cellpair_max_y, cellpair_max_y

        cj = csize_y * real(max(0, abs(j) - 1), wp)
        cj = cj * cj

        do i = -cellpair_max_x, cellpair_max_x

          ci = csize_x * real(max(0, abs(i) - 1), wp)
          ci = ci * ci

          if ((ci+cj+ck) < pairdist_w_margin2) then
            boundary%num_neighbor_cells = boundary%num_neighbor_cells + 1
          endif

        end do
      end do
    end do

    boundary%num_cells_x = ncell_x
    boundary%num_cells_y = ncell_y
    boundary%num_cells_z = ncell_z
    boundary%cell_size_x = csize_x
    boundary%cell_size_y = csize_y
    boundary%cell_size_z = csize_z
    boundary%num_cells   = ncell

    ! prepare cell neighbor list
    !
    call alloc_boundary(boundary, BoundaryCells, ncell)

    ln = 0
    do k = -cellpair_max_z, cellpair_max_z

      ck = csize_z * real(max(0, abs(k) - 1), wp)
      ck = ck * ck

      do j = -cellpair_max_y, cellpair_max_y

        cj = csize_y * real(max(0, abs(j) - 1), wp)
        cj = cj * cj

        do i = -cellpair_max_x, cellpair_max_x

          ci = csize_x * real(max(0, abs(i) - 1), wp)
          ci = ci * ci

          if ((ci+cj+ck) < pairdist_w_margin2) then
            ln = ln + 1
            boundary%neighbor_cell_common_x(ln) = i
            boundary%neighbor_cell_common_y(ln) = j
            boundary%neighbor_cell_common_z(ln) = k
          end if

        end do
      end do
    end do

    do k = 0, ncell_z-1
      do j = 0, ncell_y-1
        do i = 0, ncell_x-1
              
          lc = 1 + i + j*ncell_x + k*ncell_x*ncell_y

          do ln = 1, boundary%num_neighbor_cells

            inb = i + boundary%neighbor_cell_common_x(ln)
            if (inb < 0) then
              inbs = ncell_x + inb
            else if (inb >= ncell_x) then
              inbs = inb - ncell_x
            else
              inbs = inb
            end if

            jnb = j + boundary%neighbor_cell_common_y(ln)
            if (jnb < 0) then
              jnbs = ncell_y + jnb
            else if (jnb >= ncell_y) then
              jnbs = jnb - ncell_y
            else
              jnbs = jnb
            end if

            knb = k + boundary%neighbor_cell_common_z(ln)
            if (knb < 0) then
              knbs = ncell_z + knb
            else if (knb >= ncell_z) then
              knbs = knb - ncell_z
            else
              knbs = knb
            end if

            lcnb = 1 + inbs + jnbs*ncell_x + knbs*ncell_x*ncell_y
            boundary%neighbor_cells(ln,lc) = lcnb

          end do

        end do
      end do
    end do

    return

  end subroutine update_boundary_pbc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_boundary_cg
  !> @brief        update variables for boundary condition with CG models
  !! @authors      CT
  !! @param[in]    pairlistdist_ele    : cg electrostatics pairlist distance
  !! @param[in]    pairlistdist_126    : cg 12-6 type potential pairlist distance
  !! @param[in]    pairlistdist_PWMcos : cg PWMcos pairlist distance
  !! @param[in]    pairlistdist_DNAbp  : cg DNA basepairing pairlist distance
  !! @param[in]    pairlistdist_exv    : cg excluded volume pairlist distance
  !! @param[inout] boundary            : information of boundary condition
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_boundary_cg(pairlistdist_ele, pairlistdist_126, &
      pairlistdist_PWMcos, pairlistdist_DNAbp, pairlistdist_exv, boundary)

    ! formal arguments
    real(wp),         intent(in)    :: pairlistdist_ele
    real(wp),         intent(in)    :: pairlistdist_126
    real(wp),         intent(in)    :: pairlistdist_PWMcos
    real(wp),         intent(in)    :: pairlistdist_DNAbp
    real(wp),         intent(in)    :: pairlistdist_exv
    type(s_boundary), intent(inout) :: boundary


    select case (boundary%type)

    case (BoundaryTypeNOBC)

      ! do nothing

    case (BoundaryTypePBC)

      call update_boundary_pbc_cg(pairlistdist_ele, pairlistdist_126, &
          pairlistdist_PWMcos, pairlistdist_DNAbp, pairlistdist_exv, boundary)

    end select

    return

  end subroutine update_boundary_cg

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_boundary_pbc_cg
  !> @brief        update variables for boundary condition with CG models
  !! @authors      CT
  !! @param[in]    pairlistdist_ele    : cg electrostatics pairlist distance
  !! @param[in]    pairlistdist_126    : cg 12-6 type potential pairlist distance
  !! @param[in]    pairlistdist_PWMcos : cg PWMcos pairlist distance
  !! @param[in]    pairlistdist_DNAbp  : cg DNA basepairing pairlist distance
  !! @param[in]    pairlistdist_exv    : cg excluded volume pairlist distance
  !! @param[inout] boundary            : information of boundary condition
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_boundary_pbc_cg(pairlistdist_ele, pairlistdist_126, &
      pairlistdist_PWMcos, pairlistdist_DNAbp, pairlistdist_exv, boundary)

    ! formal arguments
    real(wp),         intent(in)    :: pairlistdist_ele
    real(wp),         intent(in)    :: pairlistdist_126
    real(wp),         intent(in)    :: pairlistdist_PWMcos
    real(wp),         intent(in)    :: pairlistdist_DNAbp
    real(wp),         intent(in)    :: pairlistdist_exv
    type(s_boundary), intent(inout) :: boundary

    ! local variables
    real(wp)                        :: bsize_x, bsize_y, bsize_z
    real(wp)                        :: csize_x, csize_y, csize_z
    integer                         :: i, j, k
    integer                         :: inb, jnb, knb, inbs, jnbs, knbs, lc, lcnb
    integer                         :: l_ele, l_126, l_PWMcos, l_DNAbp, l_exv
    integer                         :: ncell_x, ncell_y, ncell_z, ncell
    integer                         :: ncell_x_check, ncell_y_check, ncell_z_check

    ! for small cell division
    real(wp)                        :: pairlistdist_short, pairlistdist_long
    real(wp)                        :: pairdist_ele_sqr, pairdist_126_sqr
    real(wp)                        :: pairdist_PWMcos_sqr, pairdist_DNAbp_sqr
    real(wp)                        :: pairdist_exv_sqr
    real(wp)                        :: pairdist_check
    real(wp)                        :: ci, cj, ck
    integer                         :: cellpair_max_x, cellpair_max_y, cellpair_max_z
    integer                         :: ik, jk, kk


    boundary%use_cell_linked_list = .true.
    bsize_x = boundary%box_size_x
    bsize_y = boundary%box_size_y
    bsize_z = boundary%box_size_z

    ! -----------------------------------------------
    ! find out the shortest and longest pairlist dist
    ! -----------------------------------------------
    ! 
    pairlistdist_long  = max(pairlistdist_ele, pairlistdist_126, &
        pairlistdist_PWMcos, pairlistdist_DNAbp, pairlistdist_exv)
    ! 
    pairlistdist_short = min(pairlistdist_ele, pairlistdist_126, &
        pairlistdist_PWMcos, pairlistdist_DNAbp, pairlistdist_exv)

    ! -------------------------------
    ! check box size and cell numbers
    ! -------------------------------
    ! 
    ncell_x = int(bsize_x/pairlistdist_long)
    ncell_y = int(bsize_y/pairlistdist_long)
    ncell_z = int(bsize_z/pairlistdist_long)
    ! 
    if (ncell_x < 3 .or. ncell_y < 3 .or. ncell_z < 3) then
      ! if (pairlistdist_short < 0.5_wp*bsize_x .and. &
      !     pairlistdist_short < 0.5_wp*bsize_y .and. &
      !     pairlistdist_short < 0.5_wp*bsize_z) then
      !   boundary%use_cell_linked_list = .false.
      !   return
      ! else
      call error_msg('Update_Boundary_Pbc> too small boxsize/pairdist_long.'//&
          ' larger boxsize or shorter pairdist should be used.')
      ! end if
    end if

    ! ------------------------------------------
    ! reassign cells based on pairlistdist_short
    ! ------------------------------------------
    ! 
    if ( boundary%pairlist_grid < pairlistdist_short ) then
      ncell_x = int(bsize_x/pairlistdist_short)
      ncell_y = int(bsize_y/pairlistdist_short)
      ncell_z = int(bsize_z/pairlistdist_short)
    else
      ncell_x = int(bsize_x/boundary%pairlist_grid)
      ncell_y = int(bsize_y/boundary%pairlist_grid)
      ncell_z = int(bsize_z/boundary%pairlist_grid)
    end if

    ! -------------------------
    ! cell size and cell number
    ! -------------------------
    ! 
    csize_x = bsize_x/real(ncell_x,wp)
    csize_y = bsize_y/real(ncell_y,wp)
    csize_z = bsize_z/real(ncell_z,wp)
    ncell   = ncell_x*ncell_y*ncell_z

    boundary%num_cells_x = ncell_x
    boundary%num_cells_y = ncell_y
    boundary%num_cells_z = ncell_z
    boundary%cell_size_x = csize_x
    boundary%cell_size_y = csize_y
    boundary%cell_size_z = csize_z
    boundary%num_cells   = ncell


    ! ===================================
    ! assign cells into neighboring cells
    ! ===================================
    ! 
    pairdist_ele_sqr    = pairlistdist_ele    * pairlistdist_ele
    pairdist_126_sqr    = pairlistdist_126    * pairlistdist_126
    pairdist_PWMcos_sqr = pairlistdist_PWMcos * pairlistdist_PWMcos
    pairdist_DNAbp_sqr  = pairlistdist_DNAbp  * pairlistdist_DNAbp
    pairdist_exv_sqr    = pairlistdist_exv    * pairlistdist_exv

    cellpair_max_x = int(pairlistdist_long / csize_x) + 1
    cellpair_max_y = int(pairlistdist_long / csize_y) + 1
    cellpair_max_z = int(pairlistdist_long / csize_z) + 1
    boundary%num_neighbor_cells_CG_ele    = 0
    boundary%num_neighbor_cells_CG_126    = 0
    boundary%num_neighbor_cells_CG_PWMcos = 0
    boundary%num_neighbor_cells_CG_DNAbp  = 0
    boundary%num_neighbor_cells_CG_exv    = 0

    do k = -cellpair_max_z, cellpair_max_z

      ck = csize_z * real(max(0, abs(k) - 1), wp)
      ck = ck * ck

      do j = -cellpair_max_y, cellpair_max_y

        cj = csize_y * real(max(0, abs(j) - 1), wp)
        cj = cj * cj

        do i = -cellpair_max_x, cellpair_max_x

          ci = csize_x * real(max(0, abs(i) - 1), wp)
          ci = ci * ci

          if ((ci + cj + ck) < pairdist_ele_sqr) then
            boundary%num_neighbor_cells_CG_ele = boundary%num_neighbor_cells_CG_ele + 1
          end if
          if ((ci + cj + ck) < pairdist_126_sqr) then
            boundary%num_neighbor_cells_CG_126 = boundary%num_neighbor_cells_CG_126 + 1
          end if
          if ((ci + cj + ck) < pairdist_PWMcos_sqr) then
            boundary%num_neighbor_cells_CG_PWMcos = boundary%num_neighbor_cells_CG_PWMcos + 1
          end if
          if ((ci + cj + ck) < pairdist_DNAbp_sqr) then
            boundary%num_neighbor_cells_CG_DNAbp = boundary%num_neighbor_cells_CG_DNAbp + 1
          end if
          if ((ci + cj + ck) < pairdist_exv_sqr) then
            boundary%num_neighbor_cells_CG_exv = boundary%num_neighbor_cells_CG_exv + 1
          end if

        end do
      end do
    end do

    call alloc_boundary(boundary, BoundaryCellsCGele,    ncell)
    call alloc_boundary(boundary, BoundaryCellsCG126,    ncell)
    call alloc_boundary(boundary, BoundaryCellsCGPWMcos, ncell)
    call alloc_boundary(boundary, BoundaryCellsCGDNAbp,  ncell)
    call alloc_boundary(boundary, BoundaryCellsCGexv,    ncell)

    ! ------------------------------------------------
    ! count "common" neighbor cells for each potential
    ! ------------------------------------------------
    ! 
    l_ele    = 0
    l_126    = 0
    l_PWMcos = 0
    l_DNAbp  = 0
    l_exv    = 0
    ! 
    do k = -cellpair_max_z, cellpair_max_z

      ck = csize_z * real(max(0, abs(k) - 1), wp)
      ck = ck * ck

      do j = -cellpair_max_y, cellpair_max_y

        cj = csize_y * real(max(0, abs(j) - 1), wp)
        cj = cj * cj

        do i = -cellpair_max_x, cellpair_max_x

          ci = csize_x * real(max(0, abs(i) - 1), wp)
          ci = ci * ci

          if ((ci + cj + ck) < pairdist_ele_sqr) then
            l_ele = l_ele + 1
            boundary%neighbor_cell_common_x_ele(l_ele) = i
            boundary%neighbor_cell_common_y_ele(l_ele) = j
            boundary%neighbor_cell_common_z_ele(l_ele) = k
          end if

          if ((ci + cj + ck) < pairdist_126_sqr) then
            l_126 = l_126 + 1
            boundary%neighbor_cell_common_x_126(l_126) = i
            boundary%neighbor_cell_common_y_126(l_126) = j
            boundary%neighbor_cell_common_z_126(l_126) = k
          end if

          if ((ci + cj + ck) < pairdist_PWMcos_sqr) then
            l_PWMcos = l_PWMcos + 1
            boundary%neighbor_cell_common_x_PWMcos(l_PWMcos) = i
            boundary%neighbor_cell_common_y_PWMcos(l_PWMcos) = j
            boundary%neighbor_cell_common_z_PWMcos(l_PWMcos) = k
          end if

          if ((ci + cj + ck) < pairdist_DNAbp_sqr) then
            l_DNAbp = l_DNAbp + 1
            boundary%neighbor_cell_common_x_DNAbp(l_DNAbp) = i
            boundary%neighbor_cell_common_y_DNAbp(l_DNAbp) = j
            boundary%neighbor_cell_common_z_DNAbp(l_DNAbp) = k
          end if

          if ((ci + cj + ck) < pairdist_exv_sqr) then
            l_exv = l_exv + 1
            boundary%neighbor_cell_common_x_exv(l_exv) = i
            boundary%neighbor_cell_common_y_exv(l_exv) = j
            boundary%neighbor_cell_common_z_exv(l_exv) = k
          end if

        end do                  ! i
      end do                    ! j
    end do                      ! k

    ! ============================================
    ! Assign the real neighbor cells for each cell
    ! ============================================
    ! 
    do k = 0, ncell_z-1
      do j = 0, ncell_y-1
        do i = 0, ncell_x-1

          ! `lc` the index of each cell
          ! 
          lc = 1 + i + j*ncell_x + k*ncell_x*ncell_y

          ! ---
          ! ele
          ! ---
          ! 
          do l_ele = 1, boundary%num_neighbor_cells_CG_ele

            inb = i + boundary%neighbor_cell_common_x_ele(l_ele)
            if (inb < 0) then
              inbs = ncell_x + inb
            else if (inb >= ncell_x) then
              inbs = inb - ncell_x
            else
              inbs = inb
            end if

            jnb = j + boundary%neighbor_cell_common_y_ele(l_ele)
            if (jnb < 0) then
              jnbs = ncell_y + jnb
            else if (jnb >= ncell_y) then
              jnbs = jnb - ncell_y
            else
              jnbs = jnb
            end if

            knb = k + boundary%neighbor_cell_common_z_ele(l_ele)
            if (knb < 0) then
              knbs = ncell_z + knb
            else if (knb >= ncell_z) then
              knbs = knb - ncell_z
            else
              knbs = knb
            end if

            lcnb = 1 + inbs + jnbs*ncell_x + knbs*ncell_x*ncell_y
            boundary%neighbor_cells_CG_ele(l_ele,lc) = lcnb

          end do                ! l_ele

          ! ---
          ! 126
          ! ---
          ! 
          do l_126 = 1, boundary%num_neighbor_cells_CG_126

            inb = i + boundary%neighbor_cell_common_x_126(l_126)
            if (inb < 0) then
              inbs = ncell_x + inb
            else if (inb >= ncell_x) then
              inbs = inb - ncell_x
            else
              inbs = inb
            end if

            jnb = j + boundary%neighbor_cell_common_y_126(l_126)
            if (jnb < 0) then
              jnbs = ncell_y + jnb
            else if (jnb >= ncell_y) then
              jnbs = jnb - ncell_y
            else
              jnbs = jnb
            end if

            knb = k + boundary%neighbor_cell_common_z_126(l_126)
            if (knb < 0) then
              knbs = ncell_z + knb
            else if (knb >= ncell_z) then
              knbs = knb - ncell_z
            else
              knbs = knb
            end if

            lcnb = 1 + inbs + jnbs*ncell_x + knbs*ncell_x*ncell_y
            boundary%neighbor_cells_CG_126(l_126,lc) = lcnb

          end do                ! l_126

          ! ---
          ! PWMcos
          ! ---
          ! 
          do l_PWMcos = 1, boundary%num_neighbor_cells_CG_PWMcos

            inb = i + boundary%neighbor_cell_common_x_PWMcos(l_PWMcos)
            if (inb < 0) then
              inbs = ncell_x + inb
            else if (inb >= ncell_x) then
              inbs = inb - ncell_x
            else
              inbs = inb
            end if

            jnb = j + boundary%neighbor_cell_common_y_PWMcos(l_PWMcos)
            if (jnb < 0) then
              jnbs = ncell_y + jnb
            else if (jnb >= ncell_y) then
              jnbs = jnb - ncell_y
            else
              jnbs = jnb
            end if

            knb = k + boundary%neighbor_cell_common_z_PWMcos(l_PWMcos)
            if (knb < 0) then
              knbs = ncell_z + knb
            else if (knb >= ncell_z) then
              knbs = knb - ncell_z
            else
              knbs = knb
            end if

            lcnb = 1 + inbs + jnbs*ncell_x + knbs*ncell_x*ncell_y
            boundary%neighbor_cells_CG_PWMcos(l_PWMcos,lc) = lcnb

          end do                ! l_PWMcos

          ! ---
          ! DNAbp
          ! ---
          ! 
          do l_DNAbp = 1, boundary%num_neighbor_cells_CG_DNAbp

            inb = i + boundary%neighbor_cell_common_x_DNAbp(l_DNAbp)
            if (inb < 0) then
              inbs = ncell_x + inb
            else if (inb >= ncell_x) then
              inbs = inb - ncell_x
            else
              inbs = inb
            end if

            jnb = j + boundary%neighbor_cell_common_y_DNAbp(l_DNAbp)
            if (jnb < 0) then
              jnbs = ncell_y + jnb
            else if (jnb >= ncell_y) then
              jnbs = jnb - ncell_y
            else
              jnbs = jnb
            end if

            knb = k + boundary%neighbor_cell_common_z_DNAbp(l_DNAbp)
            if (knb < 0) then
              knbs = ncell_z + knb
            else if (knb >= ncell_z) then
              knbs = knb - ncell_z
            else
              knbs = knb
            end if

            lcnb = 1 + inbs + jnbs*ncell_x + knbs*ncell_x*ncell_y
            boundary%neighbor_cells_CG_DNAbp(l_DNAbp,lc) = lcnb

          end do                ! l_DNAbp

          ! ---
          ! exv
          ! ---
          ! 
          do l_exv = 1, boundary%num_neighbor_cells_CG_exv

            inb = i + boundary%neighbor_cell_common_x_exv(l_exv)
            if (inb < 0) then
              inbs = ncell_x + inb
            else if (inb >= ncell_x) then
              inbs = inb - ncell_x
            else
              inbs = inb
            end if

            jnb = j + boundary%neighbor_cell_common_y_exv(l_exv)
            if (jnb < 0) then
              jnbs = ncell_y + jnb
            else if (jnb >= ncell_y) then
              jnbs = jnb - ncell_y
            else
              jnbs = jnb
            end if

            knb = k + boundary%neighbor_cell_common_z_exv(l_exv)
            if (knb < 0) then
              knbs = ncell_z + knb
            else if (knb >= ncell_z) then
              knbs = knb - ncell_z
            else
              knbs = knb
            end if

            lcnb = 1 + inbs + jnbs*ncell_x + knbs*ncell_x*ncell_y
            boundary%neighbor_cells_CG_exv(l_exv,lc) = lcnb

          end do                ! l_exv

        end do                  ! i
      end do                    ! j
    end do                      ! k

    return

  end subroutine update_boundary_pbc_cg

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_spot_atomlist
  !> @brief        update list of atoms for spherical potential
  !!               This routine is currently not used.
  !! @authors      KY
  !! @param[in]    natom    : number of atoms
  !! @param[in]    coord    : coordinates
  !! @param[inout] boundary : boundary conditions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_spot_atomlist(natom, coord, boundary)

    ! formal arguments
    integer,                 intent(in)    :: natom
    real(wp),                intent(in)    :: coord(:,:)
    type(s_boundary),        intent(inout) :: boundary

    ! local variables
    integer  :: nfunc, nc
    integer  :: i, j, n
    integer  :: id, my_id
    real(wp) :: rn, rmin, din(3), rin
#ifdef OMP
    integer  :: omp_get_thread_num, omp_get_max_threads
#endif


#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif

    nfunc = boundary%nfunctions
    nc    = 0

    !$omp parallel                     &
    !$omp private(i, j, n,             &
    !$omp         rn, rmin, din, rin,  &
    !$omp         id, my_id)           &
    !$omp reduction(+:nc) 

#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    my_id = my_city_rank * nthread + id

    do i = 1, natom
      if (mod(i-1,nproc_city*nthread) /= my_id) cycle

      ! skip fix atoms
      if (boundary%fixatm(i)) then
        boundary%atomlist(i) = -1
        cycle
      end if

      ! find the nearest center
      rmin = -1.0_wp
      do n = 1, nfunc
        rn = boundary%radius(n)

        din = coord(:,i) - boundary%center(:,n)
        rin = sqrt(din(1)*din(1) + din(2)*din(2) + din(3)*din(3))

        if (rin < rn) then
          boundary%atomlist(i) = -1
          exit

        else
          if (rmin < 0.0_wp .or. (rin - rn) < rmin) then
            rmin  = rin - rn
            boundary%atomlist(i) = n
          end if

        end if

      end do

    end do
    !$omp end parallel

    return

  end subroutine update_spot_atomlist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    wrap molecules
  !> @brief        wrap molecules
  !! @authors      TM
  !! @param[in]    molecule : molecule information
  !! @param[in]    boundary : boundary conditions information
  !! @param[inout] coord    : coordinates
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine wrap_molecules(molecule, boundary, coord)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_boundary),        intent(in)    :: boundary
    real(wp),                intent(inout) :: coord(:,:)


    if (boundary%wrap_all) then
      call wrap_all(molecule, boundary, coord)
    else
      ! do nothing
    end if

    return

  end subroutine wrap_molecules

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    wrap_all
  !> @brief        wrap molecules
  !! @authors      TM
  !! @param[in]    molecule : molecule information
  !! @param[in]    boundary : boundary conditions information
  !! @param[inout] coord    : coordinates
  !! @note         both coord_ref and coord are wrapped for restart
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine wrap_all(molecule, boundary, coord)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_boundary),        intent(in)    :: boundary
    real(wp),                intent(inout) :: coord(:,:)

    ! local variables
    real(wp)                 :: com(3), box(3), ori(3)
    integer                  :: i, j, initial, final


    box(1) = boundary%box_size_x
    box(2) = boundary%box_size_y
    box(3) = boundary%box_size_z
    ori(1) = boundary%origin_x
    ori(2) = boundary%origin_y
    ori(3) = boundary%origin_z

    do i = 1, molecule%num_molecules

      initial  = molecule%molecule_atom_no(i)
      if (i /= molecule%num_molecules) then
        final = molecule%molecule_atom_no(i+1) - 1
      else
        final = molecule%num_atoms
      end if

      ! compute center of mass of a molecule
      !
      com(1:3) = 0.0_wp
      do j = initial, final
        coord(1:3,j) = coord(1:3,j) - ori(1:3)
        com(1:3) = com(1:3) + coord(1:3,j)*molecule%mass(j)
      end do
      com(1:3) = com(1:3)/molecule%molecule_mass(i)

      ! move molecule into the unit cell
      !
      do j = initial, final
        coord(1:3,j) = coord(1:3,j) - box(1:3)*nint(com(1:3)/box(1:3))
        coord(1:3,j) = coord(1:3,j) + ori(1:3)
      end do

    end do

    return

  end subroutine wrap_all

end module at_boundary_mod
