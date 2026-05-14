!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_boundary_mod
!> @brief   Provides utilities for boundary conditions.
!! @authors Jaewoon Jung (JJ), Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_boundary_mod

  use cg_ensemble_str_mod
  use cg_boundary_str_mod
  use cg_enefunc_str_mod
  use cg_energy_mod
  use molecules_str_mod
  use fileio_rst_mod
  use fileio_control_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_pbc_info
    real(wp)            :: box_size_x     = 0.0_wp
    real(wp)            :: box_size_y     = 0.0_wp
    real(wp)            :: box_size_z     = 0.0_wp
    integer             :: num_cells_x    = 0
    integer             :: num_cells_y    = 0
    integer             :: num_cells_z    = 0
    logical             :: calc_local_pbc = .false.
    logical             :: load_balance   = .true.
  end type s_pbc_info

  type, public :: s_nobc_info
    real(wp)            :: box_size_x    = 0.0_wp
    real(wp)            :: box_size_y    = 0.0_wp
    real(wp)            :: box_size_z    = 0.0_wp
    integer             :: num_cells_x   = 0
    integer             :: num_cells_y   = 0
    integer             :: num_cells_z   = 0
  end type s_nobc_info

  type, public :: s_boundary_info
    integer             :: type          = BoundaryTypePBC
    type(s_pbc_info)    :: pbc_info
    type(s_nobc_info)   :: nobc_info
    real(wp)            :: origin_x      = 0.0_wp
    real(wp)            :: origin_y      = 0.0_wp
    real(wp)            :: origin_z      = 0.0_wp
    integer             :: min_domain_cell = 2
    integer             :: domain_x      = 0
    integer             :: domain_y      = 0
    integer             :: domain_z      = 0
    real(wp)            :: box_size_x_max = 0.0_wp
    real(wp)            :: box_size_y_max = 0.0_wp
    real(wp)            :: box_size_z_max = 0.0_wp
    real(wp)            :: box_size_x_min = 0.0_wp
    real(wp)            :: box_size_y_min = 0.0_wp
    real(wp)            :: box_size_z_min = 0.0_wp
  end type s_boundary_info

  ! subroutines
  public  :: show_ctrl_boundary
  public  :: read_ctrl_boundary
  public  :: setup_boundary
  public  :: setup_processor_number
  public  :: setup_boundary_cell

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_boundary
  !> @brief        Displays usage information for the BOUNDARY section.
  !! @authors      NT
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "md", "min"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_boundary(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode


    if (show_all) then

      select case (run_mode)

      case ('md', 'min')

        write(MsgOut,'(A)') '[BOUNDARY]'
        write(MsgOut,'(A)') 'type          = PBC       # [PBC, NOBC]'
        write(MsgOut,'(A)') '# box_size_x    = 0.0       # Box size (x) in [PBC, NOBC].'
        write(MsgOut,'(A)') '# box_size_y    = 0.0       # Box size (y) in [PBC, NOBC].'
        write(MsgOut,'(A)') '# box_size_z    = 0.0       # Box size (z) in [PBC, NOBC].'
        write(MsgOut,'(A)') 'domain_x      = 0         # Domain size (x).'
        write(MsgOut,'(A)') 'domain_y      = 0         # Domain size (y).'
        write(MsgOut,'(A)') 'domain_z      = 0         # Domain size (z).'
        write(MsgOut,'(A)') ' '

      end select

    else

      select case (run_mode)

      case ('md', 'min')

        write(MsgOut,'(A)') '[BOUNDARY]'
        write(MsgOut,'(A)') 'type          = PBC       # [PBC, NOBC]'
        write(MsgOut,'(A)') ' '

      end select

    end if


    return

  end subroutine show_ctrl_boundary
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_boundary
  !> @brief        Reads the BOUNDARY section from the control file.
  !! @authors      JJ
  !! @param[in]    handle     : unit number
  !! @param[out]   bound_info : BOUNDARY section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_boundary(handle, bound_info)

    ! parameters
    character(*),            parameter     :: Section = 'Boundary'

    ! formal arguments
    integer,                 intent(in)    :: handle
    type(s_boundary_info),   intent(inout) :: bound_info


    ! read parameters from control file
    ! 
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_type(handle, Section, 'type', &
                            bound_info%type, BoundaryTypeTypes)


    select case (bound_info%type)

    case (BoundaryTypePBC)
      call read_ctrlfile_real (handle, Section, 'box_size_x',  &
                               bound_info%pbc_info%box_size_x)
      call read_ctrlfile_real (handle, Section, 'box_size_y',  &
                               bound_info%pbc_info%box_size_y)
      call read_ctrlfile_real (handle, Section, 'box_size_z',  &
                               bound_info%pbc_info%box_size_z)
      call read_ctrlfile_logical(handle, Section, 'local_pbc',  &
          bound_info%pbc_info%calc_local_pbc)
      call read_ctrlfile_logical(handle, Section, 'load_balance',  &
          bound_info%pbc_info%load_balance)

    case (BoundaryTypeNOBC)
      call read_ctrlfile_real (handle, Section, 'box_size_x',      &
                               bound_info%nobc_info%box_size_x)
      call read_ctrlfile_real (handle, Section, 'box_size_y',      &
                               bound_info%nobc_info%box_size_y)
      call read_ctrlfile_real (handle, Section, 'box_size_z',      &
                               bound_info%nobc_info%box_size_z)
      call read_ctrlfile_real (handle, Section, 'box_size_x_max',  &
                               bound_info%box_size_x_max)
      call read_ctrlfile_real (handle, Section, 'box_size_y_max',  &
                               bound_info%box_size_y_max)
      call read_ctrlfile_real (handle, Section, 'box_size_z_max',  &
                               bound_info%box_size_z_max)
      call read_ctrlfile_real (handle, Section, 'box_size_x_min',  &
                               bound_info%box_size_x_min)
      call read_ctrlfile_real (handle, Section, 'box_size_y_min',  &
                               bound_info%box_size_y_min)
      call read_ctrlfile_real (handle, Section, 'box_size_z_min',  &
                               bound_info%box_size_z_min)
      if (bound_info%box_size_x_max < EPS) &
      bound_info%box_size_x_max = max(bound_info%box_size_x_max, &
                                      bound_info%box_size_x_min, &
                                      bound_info%nobc_info%box_size_x)
      if (bound_info%box_size_y_max < EPS) &
      bound_info%box_size_y_max = max(bound_info%box_size_y_max, &
                                      bound_info%box_size_y_min, &
                                      bound_info%nobc_info%box_size_y)
      if (bound_info%box_size_z_max < EPS) &
      bound_info%box_size_z_max = max(bound_info%box_size_z_max, &
                                      bound_info%box_size_z_min, &
                                      bound_info%nobc_info%box_size_z)
      if (bound_info%box_size_x_min < EPS) &
      bound_info%box_size_x_min = min(bound_info%box_size_x_max, &
                                      bound_info%box_size_x_min, &
                                      bound_info%nobc_info%box_size_x)
      if (bound_info%box_size_y_min < EPS) &
      bound_info%box_size_y_min = min(bound_info%box_size_y_max, &
                                      bound_info%box_size_y_min, &
                                      bound_info%nobc_info%box_size_y)
      if (bound_info%box_size_z_min < EPS) &
      bound_info%box_size_z_min = min(bound_info%box_size_z_max, &
                                      bound_info%box_size_z_min, &
                                      bound_info%nobc_info%box_size_z)
    end select

    call read_ctrlfile_integer(handle, Section, 'min_domain_cell', &
                               bound_info%min_domain_cell)
    call read_ctrlfile_real (handle, Section, 'origin_x',  &
                             bound_info%origin_x)
    call read_ctrlfile_real (handle, Section, 'origin_y',  &
                             bound_info%origin_y)
    call read_ctrlfile_real (handle, Section, 'origin_z',  &
                             bound_info%origin_z)

    call read_ctrlfile_integer(handle, Section, 'domain_x', &
                               bound_info%domain_x)
    call read_ctrlfile_integer(handle, Section, 'domain_y', &
                               bound_info%domain_y)
    call read_ctrlfile_integer(handle, Section, 'domain_z', &
                               bound_info%domain_z)

    call end_ctrlfile_section(handle)

    ! write parameters to MsgOut
    !
    if (main_rank) then

      write(MsgOut,'(A)') &
        'Read_Ctrl_Boundary> Parameters for boundary conditions.'
      write(MsgOut,'(A20,A10)') &
            '  type            = ', BoundaryTypeTypes(bound_info%type)

      select case (bound_info%type)

      case (BoundaryTypePBC)

        write(MsgOut,'(A20,3F10.3)')                                &
            '  box_size(x,y,z) = ', bound_info%pbc_info%box_size_x, &
                                    bound_info%pbc_info%box_size_y, &
                                    bound_info%pbc_info%box_size_z

        if (bound_info%pbc_info%calc_local_pbc) &
            write(MsgOut,'(A)') '  local_pbc       = Yes'

        if (.not. bound_info%pbc_info%load_balance) &
            write(MsgOut,'(A)') '  load_balance    = No'

      case (BoundaryTypeNOBC)

        write(MsgOut,'(A24,3F10.3)')                             &
            'box_size_max(x,y,z) = ', bound_info%box_size_x_max, &
                                      bound_info%box_size_y_max, &
                                      bound_info%box_size_z_max
        write(MsgOut,'(A24,3F10.3)')                             &
            'box_size_min(x,y,z) = ', bound_info%box_size_x_min, &
                                      bound_info%box_size_y_min, &
                                      bound_info%box_size_z_min
      end select

      write(MsgOut,'(A20,I10)')  &
            '  min_domain_cell = ', bound_info%min_domain_cell

!     if (bound_info%min_domain_cell /= 1 .and. &
!         bound_info%min_domain_cell /= 2)      &
!       call error_msg('min_domain_cell should be either 1 or 2')
      
      if (bound_info%domain_x /= 0 .and. &
          bound_info%domain_y /= 0 .and. &
          bound_info%domain_z /= 0) &
        write(MsgOut,'(A20,3I10)')                         &
              '  domain (x,y,z)  = ', bound_info%domain_x, &
                                      bound_info%domain_y, &
                                      bound_info%domain_z

      write(MsgOut,'(A)') ' '

    end if

    return

  end subroutine read_ctrl_boundary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_boundary
  !> @brief        Sets essential variables for boundary conditions.
  !! @authors      JJ
  !! @param[in]    bound_info  : BOUNDARY section control parameters information
  !! @param[in]    table        : flag for use table or not
  !! @param[in]    pairlistdist : pair-list distance
  !! @param[in]    ensemble     : type of ensemble 
  !! @param[in]    rst          : restart file information
  !! @param[inout] boundary     : boundary information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_boundary(bound_info, ene_info, ensemble,          &
                            pairlistdist_ele, pairlistdist_126,      &
                            pairlistdist_WF,                         &
                            pairlistdist_PWMcos, pairlistdist_DNAbp, &
                            pairlistdist_exv, cellsize_min,          &
                            molecule, rst, boundary)

    ! formal arguments
    type(s_boundary_info),   intent(in)    :: bound_info
    type(s_ene_info),        intent(in)    :: ene_info
    integer,                 intent(in)    :: ensemble
    real(wp),                intent(in)    :: pairlistdist_ele
    real(wp),                intent(in)    :: pairlistdist_126
    real(wp),                intent(in)    :: pairlistdist_WF
    real(wp),                intent(in)    :: pairlistdist_PWMcos
    real(wp),                intent(in)    :: pairlistdist_DNAbp
    real(wp),                intent(in)    :: pairlistdist_exv
    real(wp),                intent(in)    :: cellsize_min
    type(s_rst),             intent(in)    :: rst
    type(s_molecule),        intent(in)    :: molecule
    type(s_boundary),        intent(inout) :: boundary

    ! local variables
    integer                  :: i
    real(dp)                 :: coord_min(3), coord_max(3), box_size(3)
    real(dp)                 :: center(3)

    call init_boundary(boundary)

    if (bound_info%domain_x == 1 .or. &
        bound_info%domain_y == 1 .or. &
        bound_info%domain_z == 1 ) then
      call error_msg('Error: Setup_Boundary - For domain[x,y,z], values '//&
           'other than (2,1,1), (2,2,1), or (1,1,1) must be greater than 1.')
    end if

    select case (bound_info%type)

    case (BoundaryTypeNOBC)

      boundary%type           = bound_info%type
      boundary%box_size_x     = bound_info%nobc_info%box_size_x    
      boundary%box_size_y     = bound_info%nobc_info%box_size_y    
      boundary%box_size_z     = bound_info%nobc_info%box_size_z    
      boundary%box_size_x_max = bound_info%box_size_x_max
      boundary%box_size_y_max = bound_info%box_size_y_max
      boundary%box_size_z_max = bound_info%box_size_z_max
      boundary%box_size_x_min = bound_info%box_size_x_min
      boundary%box_size_y_min = bound_info%box_size_y_min
      boundary%box_size_z_min = bound_info%box_size_z_min
      boundary%min_domain_cell= bound_info%min_domain_cell

      center(1:3) = 0.0_dp
      if (rst%rstfile_type /= RstfileTypeUndef) then
        do i = 1, rst%num_atoms
          center(1:3) = center(1:3) + rst%coord(1:3,i)
        end do
        center(1:3) = center(1:3) / real(rst%num_atoms,dp)
      else
        do i = 1, molecule%num_atoms
          center(1:3) = center(1:3) + molecule%atom_coord(1:3,i)
        end do
        center(1:3) = center(1:3) / real(molecule%num_atoms,dp)
      end if
      boundary%origin_x = center(1)
      boundary%origin_y = center(2)
      boundary%origin_z = center(3)

      ! Decidie the system size if it is not written
      !
      if (boundary%box_size_x < EPS) then
        coord_min(1) =  1000000000000.0_wp
        coord_max(1) = -1000000000000.0_wp
        if (rst%rstfile_type /= RstfileTypeUndef) then
          do i = 1, rst%num_atoms
            coord_min(1) = min(coord_min(1), rst%coord(1,i))
            coord_max(1) = max(coord_max(1), rst%coord(1,i))
          end do
        else
          do i = 1, molecule%num_atoms
            coord_min(1) = min(coord_min(1), molecule%atom_coord(1,i))
            coord_max(1) = max(coord_max(1), molecule%atom_coord(1,i))
          end do
        end if
        box_size(1) = max(-coord_min(1)+0.1_wp,coord_max(1)+0.1_wp)
        boundary%box_size_x     = box_size(1)*2.0_wp
        if (boundary%box_size_x > boundary%box_size_x_max) &
          call error_msg('Error: Calculated box_size_x is greater than '//&
                         'box_size_x_max.')
        if (boundary%box_size_x < boundary%box_size_x_min) &
          write(Msgout, '(A)') 'Warning: Calculated box_size_x is less '//&
                               'than box_size_x_min.'
      end if

      if (boundary%box_size_y < EPS) then
        coord_min(2) =  1000000000000.0_wp
        coord_max(2) = -1000000000000.0_wp
        if (rst%rstfile_type /= RstfileTypeUndef) then
          do i = 1, rst%num_atoms
            coord_min(2) = min(coord_min(2), rst%coord(2,i))
            coord_max(2) = max(coord_max(2), rst%coord(2,i))
          end do
        else
          do i = 1, molecule%num_atoms
            coord_min(2) = min(coord_min(2), molecule%atom_coord(2,i))
            coord_max(2) = max(coord_max(2), molecule%atom_coord(2,i))
          end do
        end if
        box_size(2) = max(-coord_min(2)+0.1_wp,coord_max(2)+0.1_wp)
        boundary%box_size_y     = box_size(2)*2.0_wp
        if (boundary%box_size_y > boundary%box_size_y_max) &
          call error_msg('Error: Calculated box_size_y is greater than '//&
                         'box_size_y_max.')
        if (boundary%box_size_y < boundary%box_size_y_min) &
          write(Msgout, '(A)') 'Warning: Calculated box_size_y is less '//&
                               'than box_size_y_min.'
      end if

      if (boundary%box_size_z < EPS) then
        coord_min(3) =  1000000000000.0_wp
        coord_max(3) = -1000000000000.0_wp
        if (rst%rstfile_type /= RstfileTypeUndef) then
          do i = 1, rst%num_atoms
            coord_min(3) = min(coord_min(3), rst%coord(3,i))
            coord_max(3) = max(coord_max(3), rst%coord(3,i))
          end do
        else
          do i = 1, molecule%num_atoms
            coord_min(3) = min(coord_min(3), molecule%atom_coord(3,i))
            coord_max(3) = max(coord_max(3), molecule%atom_coord(3,i))
          end do
        end if
        box_size(3) = max(-coord_min(3)+0.1_wp,coord_max(3)+0.1_wp)
        boundary%box_size_z     = box_size(3)*2.0_wp
        if (boundary%box_size_z > boundary%box_size_z_max) &
          call error_msg('Error: Calculated box_size_z is greater than '//&
                         'box_size_z_max.')
        if (boundary%box_size_z < boundary%box_size_z_min) &
          write(Msgout, '(A)') 'Warning: Calculated box_size_z is less '//&
                               'than box_size_z_min.'
      end if
 
      boundary%box_size_x_ref = boundary%box_size_x
      boundary%box_size_y_ref = boundary%box_size_y
      boundary%box_size_z_ref = boundary%box_size_z

    case (BoundaryTypePBC)

      boundary%type           = bound_info%type
      boundary%origin_x       = bound_info%origin_x
      boundary%origin_y       = bound_info%origin_y
      boundary%origin_z       = bound_info%origin_z
      boundary%box_size_x     = bound_info%pbc_info%box_size_x
      boundary%box_size_y     = bound_info%pbc_info%box_size_y
      boundary%box_size_z     = bound_info%pbc_info%box_size_z
      boundary%box_size_x_ref = boundary%box_size_x
      boundary%box_size_y_ref = boundary%box_size_y
      boundary%box_size_z_ref = boundary%box_size_z
      boundary%calc_local_pbc = bound_info%pbc_info%calc_local_pbc
      boundary%load_balance   = bound_info%pbc_info%load_balance
      boundary%min_domain_cell= bound_info%min_domain_cell
    
      if (rst%rstfile_type /= RstfileTypeUndef) then

        boundary%box_size_x     = rst%box_size_x
        boundary%box_size_y     = rst%box_size_y
        boundary%box_size_z     = rst%box_size_z
        boundary%box_size_x_ref = boundary%box_size_x
        boundary%box_size_y_ref = boundary%box_size_y
        boundary%box_size_z_ref = boundary%box_size_z

      end if

    end select

    call setup_processor_number(bound_info,          &
                                ene_info,            &
                                pairlistdist_ele,    &
                                pairlistdist_126,    &
                                pairlistdist_WF,     &
                                pairlistdist_PWMcos, &
                                pairlistdist_DNAbp,  &
                                pairlistdist_exv,    &
                                cellsize_min,        &
                                ensemble, boundary)

    call setup_boundary_cell   (ene_info,            &
                                ensemble,            &
                                pairlistdist_ele,    &
                                pairlistdist_126,    &
                                pairlistdist_WF,     &
                                pairlistdist_PWMcos, &
                                pairlistdist_DNAbp,  &
                                pairlistdist_exv,    &
                                cellsize_min,        &
                                boundary)

    return

  end subroutine setup_boundary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_processor_number
  !> @brief        Defines the number of processors in each dimension.
  !! @authors      JJ
  !! @param[in]    bound_info : BOUNDARY section control parameters information
  !! @param[in]    table        : flag for use table or not
  !! @param[in]    pairlistdist : pair-list distance
  !! @param[in]    ensemble     : type of ensemble 
  !! @param[inout] boundary     : boundary information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_processor_number(bound_info, ene_info, pairlistdist_ele,  &
                                    pairlistdist_126, pairlistdist_WF,       &
                                    pairlistdist_PWMcos,                     &
                                    pairlistdist_DNAbp, pairlistdist_exv,    &
                                    cellsize_min, ensemble, boundary)

    ! formal arguments
    type(s_boundary_info),   intent(in)    :: bound_info
    type(s_ene_info),        intent(in)    :: ene_info
    real(wp),                intent(in)    :: pairlistdist_ele
    real(wp),                intent(in)    :: pairlistdist_126
    real(wp),                intent(in)    :: pairlistdist_WF
    real(wp),                intent(in)    :: pairlistdist_PWMcos
    real(wp),                intent(in)    :: pairlistdist_DNAbp
    real(wp),                intent(in)    :: pairlistdist_exv
    real(wp),                intent(in)    :: cellsize_min
    integer,                 intent(in)    :: ensemble
    type(s_boundary),        intent(inout) :: boundary

    ! local variable
    real(wp)                 :: bsize_x, bsize_y, bsize_z
    real(wp)                 :: size_x, size_y, size_z
    real(wp)                 :: maxsize(0:500), cell_size(3,500)
    real(wp)                 :: cutoff
    logical                  :: accept_proc
    integer                  :: total_proc
    integer                  :: nx, ny, nz, i, j, k, itype
    integer                  :: nc(3,500)
    logical                  :: extend1

    ! check if nproc_city is the power of 2
    !
    total_proc = nproc_country

    accept_proc = .true.
    do while (total_proc > 1)
      if (mod(total_proc,2) /= 0) then
        accept_proc = .false.
        exit
      end if
      total_proc = total_proc / 2
    end do
    if (.not.accept_proc) &
      call error_msg('Error: Setup_Processor_Number - Number of processes '//&
                     'is not a power of 2.')
    if (nproc_country < 8) &
      call error_msg('Error: Setup_Processor_Number - Number of processes '//&
                     'should be greater than 8.')

    extend1 = .false.
    if (ensemble == EnsembleNPT  .or. &
        ensemble == EnsembleNPAT .or. &
        ensemble == EnsembleNPgT ) &
      extend1 = .true.

    ! check processor number based on the num of domains
    !
    if (bound_info%domain_x /= 0 .and. &
        bound_info%domain_y /= 0 .and. &
        bound_info%domain_z /= 0) then

      total_proc = bound_info%domain_x * &
                   bound_info%domain_y * &
                   bound_info%domain_z

      if (total_proc == nproc_country) then

        boundary%num_domain(1) = bound_info%domain_x
        boundary%num_domain(2) = bound_info%domain_y
        boundary%num_domain(3) = bound_info%domain_z
        return

      else
        call error_msg('Error: Setup_Processor_Number - Number of processes '//&
                     'does not match domain_x * domain_y * domain_z.')
      end if

    end if

    bsize_x = real(boundary%box_size_x,wp)
    bsize_y = real(boundary%box_size_y,wp)
    bsize_z = real(boundary%box_size_z,wp)

    if (ene_info%forcefield == ForcefieldGroMartini) then
      cutoff = max(pairlistdist_ele*2.0_wp, cellsize_min*2.0_wp)
    else
      cutoff = max(pairlistdist_ele,           pairlistdist_126,          &
                   pairlistdist_WF,                                       &
                   pairlistdist_PWMcos*2.0_wp, pairlistdist_DNAbp*2.0_wp, &
                   pairlistdist_exv*2.0_wp, cellsize_min*2.0_wp)
    end if

    if (extend1) then
      nx  = int(bsize_x / (cutoff+0.6_wp))
      ny  = int(bsize_y / (cutoff+0.6_wp))
      nz  = int(bsize_z / (cutoff+0.6_wp))
    else
      nx  = int(bsize_x / cutoff)
      ny  = int(bsize_y / cutoff)
      nz  = int(bsize_z / cutoff)
    end if

    if (nx*ny*nz >= nproc_city) then

      itype = 0
      do k = 2, nz
        do j = 2, ny
          do i = 2, nx
            if (i*j*k == nproc_city) then
              itype = itype + 1
              nc(1,itype) = i
              nc(2,itype) = j
              nc(3,itype) = k
              cell_size(1,itype) = nx/i
              cell_size(2,itype) = ny/j
              cell_size(3,itype) = nz/k
            end if
          end do
        end do
      end do

    else
      call error_msg('Error: Setup_Processor_Number - Cannot define domains '//&
                     'and cells. Smaller MPI processors, shorter '//&
                     'pairlistdist, or larger boxsize should be used.')
    end if

    if (itype == 0) then
      call error_msg('Error: Setup_Processor_Number - Cannot define domains '//&
                     'and cells. Smaller or adjusted MPI processors, '//&
                     'shorter pairlistdist, or larger boxsize should be used.')
    end if

    k = 0
    maxsize(0) = 10000000000000.0_wp
    do i = 1, itype
      size_x = bsize_x/cell_size(1,i)
      size_y = bsize_y/cell_size(2,i)
      size_z = bsize_z/cell_size(3,i)
      maxsize(i) = cell_size(1,i) + cell_size(2,i) + cell_size(3,i)
      if (maxsize(i) < maxsize(k)) &
        k = i
    end do

    boundary%num_domain(1) = nc(1,k)
    boundary%num_domain(2) = nc(2,k)
    boundary%num_domain(3) = nc(3,k)

    return

  end subroutine setup_processor_number

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_boundary_cell
  !> @brief        Sets up boundary cell information.
  !! @authors      NT
  !! @param[in]    table        : flag for use table or not
  !! @param[in]    pairlistdist : pair-list distance
  !! @param[in]    ensemble     : type of ensemble 
  !! @param[inout] boundary     : boundary information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_boundary_cell(ene_info, ensemble, pairlistdist_ele,         &
                                 pairlistdist_126, pairlistdist_WF,            &
                                 pairlistdist_PWMcos, pairlistdist_DNAbp,      &
                                 pairlistdist_exv, cellsize_min, boundary)

    ! formal arguments
    type(s_ene_info),         intent(in)    :: ene_info
    integer,                  intent(in)    :: ensemble
    real(wp),                 intent(in)    :: pairlistdist_ele
    real(wp),                 intent(in)    :: pairlistdist_126
    real(wp),                 intent(in)    :: pairlistdist_WF
    real(wp),                 intent(in)    :: pairlistdist_PWMcos
    real(wp),                 intent(in)    :: pairlistdist_DNAbp
    real(wp),                 intent(in)    :: pairlistdist_exv
    real(wp),                 intent(in)    :: cellsize_min
    type(s_boundary), target, intent(inout) :: boundary

    ! local variables
    real(wp)                 :: csize_x, csize_y, csize_z, cutoff
    real(wp)                 :: bsize_x, bsize_y, bsize_z
    integer                  :: i, j, k, ncell
    integer                  :: inb,jnb,knb,inbs,jnbs,knbs,lc,lcnb,ln
    integer                  :: ncell_x, ncell_y, ncell_z
    logical                  :: extend, extend1
    integer, pointer         :: num_domain(:)

    num_domain  => boundary%num_domain

    extend  = .false.

    extend1 = .false.
    if (ensemble == EnsembleNPT  .or. &
        ensemble == EnsembleNPAT .or. &
        ensemble == EnsembleNPgT) &
      extend1 = .true.

    if (ene_info%forcefield == ForcefieldGroMartini) then
      cutoff = max(pairlistdist_ele*2.0_wp, cellsize_min*2.0_wp)
    else
      cutoff = max(pairlistdist_ele,           pairlistdist_126,          &
                   pairlistdist_WF,                                       &
                   pairlistdist_PWMcos*2.0_wp, pairlistdist_DNAbp*2.0_wp, &
                   pairlistdist_exv*2.0_wp,    cellsize_min*2.0_wp)
    end if

    if (extend1) cutoff = cutoff + 1.2_wp

    bsize_x = real(boundary%box_size_x,wp)
    bsize_y = real(boundary%box_size_y,wp)
    bsize_z = real(boundary%box_size_z,wp)

    if ( (boundary%num_cells_x /= 0) .and. &
         (boundary%num_cells_y /= 0) .and. &
         (boundary%num_cells_z /= 0)) then
      ncell_x = boundary%num_cells_x
      ncell_y = boundary%num_cells_y
      ncell_z = boundary%num_cells_z
    else
      ncell_x = int(bsize_x/(cutoff/2.0_wp))
      ncell_y = int(bsize_y/(cutoff/2.0_wp))
      ncell_z = int(bsize_z/(cutoff/2.0_wp))
    end if

#ifdef DEBUG
    if (main_rank) then
      write(MsgOut,'(a,f15.8)')  'Debugging > cutoff', cutoff
      write(MsgOut,'(a,3f15.8)') 'Debugging > bsize_[x,y,z]', &
                                 bsize_x, bsize_y, bsize_z
      write(MsgOut,'(a,3i8)')    'Debugging > ncell_[x,y,z]', &
                                 ncell_x, ncell_y, ncell_z
    end if
#endif

    if (ncell_x <= 4 .or. ncell_y <= 4 .or. ncell_z <= 4) &
      call error_msg('Error: Setup_Boundary_Cell - Boxsize/pairlistdist is '//&
                     'too small. Shorter pairlistdist, larger boxsize, or '//&
                     'fewer MPI processors should be used.')
    if (.not. boundary%load_balance) then
      k = mod(ncell_x, num_domain(1))
      ncell_x = ncell_x - k
      k = mod(ncell_y, num_domain(2))
      ncell_y = ncell_y - k
      k = mod(ncell_z, num_domain(3))
      ncell_z = ncell_z - k
    end if

    csize_x = bsize_x/real(ncell_x, wp)
    csize_y = bsize_y/real(ncell_y, wp)
    csize_z = bsize_z/real(ncell_z, wp)
    ncell   = ncell_x*ncell_y*ncell_z

    boundary%num_cells_x      = ncell_x
    boundary%num_cells_y      = ncell_y
    boundary%num_cells_z      = ncell_z
    boundary%cell_size_x      = csize_x
    boundary%cell_size_y      = csize_y
    boundary%cell_size_z      = csize_z

    ! prepare cell neighbor list
    !
    if (ene_info%forcefield == ForcefieldGroMartini) then
      call alloc_boundary(boundary, BoundaryCells_Martini, ncell)
    else
      call alloc_boundary(boundary, BoundaryCells_AICG2P, ncell)
    end if

    do k = 0, ncell_z-1
    do j = 0, ncell_y-1
    do i = 0, ncell_x-1

      lc = 1 + i + j*ncell_x + k*ncell_x*ncell_y
      ln = 0

      do knb = k-1, k+1

        if (knb == -1) then
          knbs = ncell_z - 1
        else if (knb == ncell_z) then
          knbs = 0
        else
          knbs = knb
        end if

        do jnb = j-1, j+1

          if (jnb == -1) then
            jnbs = ncell_y - 1
          else if (jnb == ncell_y) then
           jnbs = 0
          else
            jnbs = jnb
          end if

          do inb = i-1, i+1
           if (inb == -1) then
              inbs = ncell_x - 1
            else if (inb == ncell_x) then
              inbs = 0
            else
              inbs = inb
            end if

            lcnb = 1 + inbs + jnbs*ncell_x + knbs*ncell_x*ncell_y
            ln   = ln + 1

            boundary%neighbor_cells(ln,lc) = lcnb
          end do

        end do
      end do

      if (ene_info%forcefield /= ForcefieldGroMartini) then

      do knb = k-2, k+2

        if (knb == -2) then
          knbs = ncell_z - 2
        else if (knb == -1) then
          knbs = ncell_z - 1
        else if (knb == ncell_z) then
          knbs = 0
        else if (knb == (ncell_z+1)) then
          knbs = 1
        else
          knbs = knb
        end if

        do jnb = j-2, j+2

          if (jnb == -2) then
            jnbs = ncell_y - 2
          else if (jnb == -1) then
            jnbs = ncell_y - 1
          else if (jnb == ncell_y) then
            jnbs = 0
          else if (jnb == (ncell_y+1)) then
            jnbs = 1
          else
            jnbs = jnb
          end if

          do inb = i-2, i+2
            if (inb == -2) then
              inbs = ncell_x - 2
            else if (inb == -1) then
              inbs = ncell_x - 1
            else if (inb == ncell_x) then
              inbs = 0
            else if (inb == (ncell_x+1)) then
              inbs = 1
            else
              inbs = inb
            end if

            if (abs(inb-i) >= 2 .or. abs(jnb-j) >= 2 .or. abs(knb-k) >= 2) then
              lcnb = 1 + inbs + jnbs*ncell_x + knbs*ncell_x*ncell_y
              ln   = ln + 1
              boundary%neighbor_cells(ln,lc) = lcnb
            end if

            boundary%neighbor_cells(ln,lc) = lcnb
          end do

        end do
      end do

      end if

    end do
    end do
    end do

    if (main_rank) then
      write(MsgOut,'(A)') &
           'Setup_Boundary_Cell> Set variables for boundary conditions.'

      write(MsgOut,'(A20,3I10)')                 &
           '  domains (x,y,z) = ', boundary%num_domain(1), &
                                   boundary%num_domain(2), &
                                   boundary%num_domain(3)

      write(MsgOut,'(A20,3I10)')            &
           '  ncells (x,y,z)  = ', ncell_x, &
                                   ncell_y, &
                                   ncell_z
      write(MsgOut,'(A)') ' '
    end if

    if (ncell_x <= boundary%num_domain(1) .or. &
        ncell_y <= boundary%num_domain(2) .or. &
        ncell_z <= boundary%num_domain(3))     &
      call error_msg('Error: Setup_Boundary_Cell - ncell_[x,y,z] should be '//&
                     'greater than or equal to 2*domain_[x,y,z]. Please '//&
                     'reduce MPI and increase OpenMP to use the same '//&
                     'number of processors.')
    return

  end subroutine setup_boundary_cell

end module cg_boundary_mod
