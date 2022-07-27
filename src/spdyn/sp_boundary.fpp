!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_boundary_mod
!> @brief   utilities for boundary conditions
!! @authors Jaewoon Jung (JJ), Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_boundary_mod

  use sp_ensemble_str_mod
  use sp_boundary_str_mod
  use fileio_rst_mod
  use fileio_control_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_pbc_info
    real(wp)            :: box_size_x    = 0.0_wp
    real(wp)            :: box_size_y    = 0.0_wp
    real(wp)            :: box_size_z    = 0.0_wp
    real(wp)            :: cell_size_buffer = -1.0_wp
    integer             :: num_cells_x   = 0
    integer             :: num_cells_y   = 0
    integer             :: num_cells_z   = 0
  end type s_pbc_info

  type, public :: s_boundary_info
    integer             :: type          = BoundaryTypePBC
    type(s_pbc_info)    :: pbc_info
    real(wp)            :: origin_x      = 0.0_wp
    real(wp)            :: origin_y      = 0.0_wp
    real(wp)            :: origin_z      = 0.0_wp
    integer             :: domain_x      = 0
    integer             :: domain_y      = 0
    integer             :: domain_z      = 0
    integer             :: domain_x_max  = 0
    integer             :: domain_y_max  = 0
    integer             :: domain_z_max  = 0
    integer             :: pio_domain_x  = 0
    integer             :: pio_domain_y  = 0
    integer             :: pio_domain_z  = 0
    integer             :: duplicate_x   = 1
    integer             :: duplicate_y   = 1
    integer             :: duplicate_z   = 1
    logical             :: shift_origin  = .true.
  end type s_boundary_info

  ! subroutines
  public  :: show_ctrl_boundary
  public  :: read_ctrl_boundary
  public  :: setup_boundary
  public  :: setup_boundary_pio
  public  :: setup_processor_number
  public  :: setup_boundary_cell

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_boundary
  !> @brief        show BOUNDARY section usage
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

      case ('md', 'min', 'remd')

        write(MsgOut,'(A)') '[BOUNDARY]'
        write(MsgOut,'(A)') 'type          = PBC       # [PBC]'
        write(MsgOut,'(A)') '# box_size_x    = 0.0       # box size (x) in [PBC]'
        write(MsgOut,'(A)') '# box_size_y    = 0.0       # box size (y) in [PBC]'
        write(MsgOut,'(A)') '# box_size_z    = 0.0       # box size (z) in [PBC]'
        write(MsgOut,'(A)') '# cell_size_buffer = 0.0    # buffer added in cell size'

        write(MsgOut,'(A)') 'domain_x      = 0         # domain size (x)'
        write(MsgOut,'(A)') 'domain_y      = 0         # domain size (y)'
        write(MsgOut,'(A)') 'domain_z      = 0         # domain size (z)'
        write(MsgOut,'(A)') 'duplicate_x   = 1         # domain size (x)'
        write(MsgOut,'(A)') 'duplicate_y   = 1         # domain size (x)'
        write(MsgOut,'(A)') 'duplicate_z   = 1         # domain size (x)'
        write(MsgOut,'(A)') ' '

      end select

    else

      select case (run_mode)

      case ('md', 'min', 'remd')

        write(MsgOut,'(A)') '[BOUNDARY]'
        write(MsgOut,'(A)') 'type          = PBC       # [PBC]'
        write(MsgOut,'(A)') ' '

      end select

    end if

    return

  end subroutine show_ctrl_boundary
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_boundary
  !> @brief        read BOUNDARY section in the control file
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

      call read_ctrlfile_real (handle, Section, 'box_size_x',   &
                               bound_info%pbc_info%box_size_x)
      call read_ctrlfile_real (handle, Section, 'box_size_y',   &
                               bound_info%pbc_info%box_size_y)
      call read_ctrlfile_real (handle, Section, 'box_size_z',   &
                               bound_info%pbc_info%box_size_z)
      call read_ctrlfile_real (handle, Section, 'cell_size_buffer',   &
                               bound_info%pbc_info%cell_size_buffer)

    end select

    call read_ctrlfile_integer(handle, Section, 'domain_x',     &
                               bound_info%domain_x)
    call read_ctrlfile_integer(handle, Section, 'domain_y',     &
                               bound_info%domain_y)
    call read_ctrlfile_integer(handle, Section, 'domain_z',     &
                               bound_info%domain_z)
    call read_ctrlfile_integer(handle, Section, 'domain_x_max', &
                               bound_info%domain_x_max)
    call read_ctrlfile_integer(handle, Section, 'domain_y_max', &
                               bound_info%domain_y_max)
    call read_ctrlfile_integer(handle, Section, 'domain_z_max', &
                               bound_info%domain_z_max)
    call read_ctrlfile_integer(handle, Section, 'pio_domain_x', &
                               bound_info%pio_domain_x)
    call read_ctrlfile_integer(handle, Section, 'pio_domain_y', &
                               bound_info%pio_domain_y)
    call read_ctrlfile_integer(handle, Section, 'pio_domain_z', &
                               bound_info%pio_domain_z)

    call read_ctrlfile_integer(handle, Section, 'duplicate_x',  &
                               bound_info%duplicate_x)
    call read_ctrlfile_integer(handle, Section, 'duplicate_y',  &
                               bound_info%duplicate_y)
    call read_ctrlfile_integer(handle, Section, 'duplicate_z',  &
                               bound_info%duplicate_z)

    call end_ctrlfile_section(handle)

    ! write parameters to MsgOut
    !
    if (main_rank) then

      write(MsgOut,'(A)') 'Read_Ctrl_Boundary> Parameters of Boundary Condition'
      write(MsgOut,'(A20,A10)') &
            '  type            = ', BoundaryTypeTypes(bound_info%type)

      select case (bound_info%type)

      case (BoundaryTypePBC)

        write(MsgOut,'(A20,3F10.3)')                                &
            '  box_size(x,y,z) = ', bound_info%pbc_info%box_size_x, &
                                    bound_info%pbc_info%box_size_y, &
                                    bound_info%pbc_info%box_size_z

      end select

      if (bound_info%domain_x /= 0 .and. &
          bound_info%domain_y /= 0 .and. &
          bound_info%domain_z /= 0) &
        write(MsgOut,'(A20,3I10)')                         &
              '  domain (x,y,z)  = ', bound_info%domain_x, &
                                      bound_info%domain_y, &
                                      bound_info%domain_z

      if (bound_info%domain_x_max /= 0 .and. &
          bound_info%domain_y_max /= 0 .and. &
          bound_info%domain_z_max /= 0) &
        write(MsgOut,'(A20,3I10)')                             &
              'domain_max(x,y,z) = ', bound_info%domain_x_max, &
                                      bound_info%domain_y_max, &
                                      bound_info%domain_z_max

      if (bound_info%pio_domain_x /= 0 .and. &
          bound_info%pio_domain_y /= 0 .and. &
          bound_info%pio_domain_z /= 0) &
        write(MsgOut,'(A20,3I10)')                             &
              'pio_domain(x,y,z) = ', bound_info%pio_domain_x, &
                                      bound_info%pio_domain_y, &
                                      bound_info%pio_domain_z

      if (bound_info%duplicate_x /= 1 .and. &
          bound_info%duplicate_y /= 1 .and. &
          bound_info%duplicate_z /= 1) &
        write(MsgOut,'(A20,3I10)')                            &
              'duplicate(x,y,z)  = ', bound_info%duplicate_x, &
                                      bound_info%duplicate_y, &
                                      bound_info%duplicate_z

      write(MsgOut,'(A)') ' '

    end if

    return

  end subroutine read_ctrl_boundary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_boundary
  !> @brief        set essential variables for boundary condition
  !! @authors      JJ
  !! @param[in]    bound_info  : BOUNDARY section control parameters information
  !! @param[in]    table        : flag for use table or not
  !! @param[in]    pairlistdist : pair-list distance
  !! @pmara[in]    water_model  : water model
  !! @param[in]    ensemble     : type of ensemble 
  !! @param[in]    rigid_bond   : flag for rigid-bond
  !! @param[in]    dsize_cg     : flag for reset domain size for CG-model
  !! @param[in]    dmin_size_cg : minimum domain size for CG-model
  !! @param[in]    rst          : restart file information
  !! @param[inout] boundary     : boundary information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_boundary(bound_info, table, pairlistdist, water_model, &
                            ensemble, rigid_bond, dsize_cg, dmin_size_cg, &
                            rst, boundary)

    ! formal arguments
    type(s_boundary_info),   intent(in)    :: bound_info
    logical,                 intent(in)    :: table
    real(wp),                intent(in)    :: pairlistdist
    character(*),            intent(in)    :: water_model
    integer,                 intent(in)    :: ensemble
    logical,                 intent(in)    :: rigid_bond
    type(s_rst),             intent(in)    :: rst
    logical,                 intent(in)    :: dsize_cg
    real(wp),                intent(in)    :: dmin_size_cg
    type(s_boundary),        intent(inout) :: boundary


    call init_boundary(boundary)

    if (bound_info%domain_x == 1 .or. bound_info%domain_y == 1 .or. &
        bound_info%domain_z == 1) then

      if (.not.((bound_info%domain_x == 2 .and. &
                 bound_info%domain_y == 1 .and. &
                 bound_info%domain_z == 1) .or. &
                (bound_info%domain_x == 2 .and. &
                 bound_info%domain_y == 2 .and. &
                 bound_info%domain_z == 1) .or. &
                (bound_info%domain_x == 1 .and. &
                 bound_info%domain_y == 1 .and. &
                 bound_info%domain_z == 1))) then

        call error_msg('Setup_Boundary> other than (2,1,1)/(2,2,1)/(1,1,1), '//&
                       'domain[x,y,z] should be larger than 1')

      end if

    end if

    boundary%type             = bound_info%type
    boundary%origin_x         = bound_info%origin_x
    boundary%origin_y         = bound_info%origin_y
    boundary%origin_z         = bound_info%origin_z

    boundary%box_size_x       = bound_info%pbc_info%box_size_x
    boundary%num_duplicate(1) = bound_info%duplicate_x
    boundary%box_size_orig(1) = boundary%box_size_x
    boundary%box_size_x       = boundary%box_size_x &
                               *real(bound_info%duplicate_x)

    boundary%box_size_y       = bound_info%pbc_info%box_size_y
    boundary%num_duplicate(2) = bound_info%duplicate_y
    boundary%box_size_orig(2) = boundary%box_size_y
    boundary%box_size_y       = boundary%box_size_y &
                               *real(bound_info%duplicate_y)

    boundary%box_size_z       = bound_info%pbc_info%box_size_z
    boundary%num_duplicate(3) = bound_info%duplicate_z
    boundary%box_size_orig(3) = boundary%box_size_z
    boundary%box_size_z       = boundary%box_size_z &
                               *real(bound_info%duplicate_z)

    boundary%box_size_x_ref = boundary%box_size_x
    boundary%box_size_y_ref = boundary%box_size_y
    boundary%box_size_z_ref = boundary%box_size_z
    
    if (rst%rstfile_type /= RstfileTypeUndef) then

      boundary%box_size_x       = rst%box_size_x
      boundary%box_size_y       = rst%box_size_y
      boundary%box_size_z       = rst%box_size_z

      boundary%box_size_orig(1) = boundary%box_size_x
      boundary%box_size_x       = boundary%box_size_x &
                                 *real(bound_info%duplicate_x)
      boundary%box_size_orig(2) = boundary%box_size_y
      boundary%box_size_y       = boundary%box_size_y &
                                 *real(bound_info%duplicate_y)
      boundary%box_size_orig(3) = boundary%box_size_z
      boundary%box_size_z       = boundary%box_size_z &
                                 *real(bound_info%duplicate_z)

      boundary%box_size_x_ref = boundary%box_size_x
      boundary%box_size_y_ref = boundary%box_size_y
      boundary%box_size_z_ref = boundary%box_size_z

    end if

    call setup_processor_number(bound_info, &
                                table, pairlistdist, water_model, &
                                ensemble, rigid_bond, boundary)

    call setup_boundary_cell   (bound_info%pbc_info%cell_size_buffer, &
                                table, pairlistdist, water_model,     &
                                ensemble, rigid_bond, dsize_cg,       &
                                dmin_size_cg, boundary)

    return

  end subroutine setup_boundary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_boundary_pio
  !> @brief        set essential variables for boundary condition
  !! @authors      NT
  !! @param[in]    bound_info  : BOUNDARY section control parameters information
  !! @param[inout] boundary     : boundary information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_boundary_pio(bound_info, boundary)

    ! formal arguments
    type(s_boundary_info),   intent(in)    :: bound_info
    type(s_boundary),        intent(inout) :: boundary


    call init_boundary(boundary)

    if (.not.((bound_info%domain_x == 2 .and. &
               bound_info%domain_y == 1 .and. &
               bound_info%domain_z == 1) .or. &
              (bound_info%domain_x == 2 .and. &
               bound_info%domain_y == 2 .and. &
               bound_info%domain_z == 1) .or. &
              (bound_info%domain_x == 1 .and. &
               bound_info%domain_y == 1 .and. &
               bound_info%domain_z == 1))) then

      if (bound_info%domain_x == 1 .or. &
          bound_info%domain_y == 1 .or. &
          bound_info%domain_z == 1 ) then

        call error_msg('Setup_Boundary_Pio> other than (2,1,1)/(2,2,1)/' &
                       //'(1,1,1), domain[x,y,z] should be larger than 1')

      end if

    end if

    call setup_processor_number_pio(bound_info, boundary)

    return

  end subroutine setup_boundary_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_processor_number
  !> @brief        define the processor number in each dimension
  !! @authors      JJ
  !! @param[in]    bound_info : BOUNDARY section control parameters information
  !! @param[in]    table        : flag for use table or not
  !! @param[in]    pairlistdist : pair-list distance
  !! @pmara[in]    water_model  : water model
  !! @param[in]    ensemble     : type of ensemble 
  !! @param[in]    rigid_bond   : flag for rigid-bond
  !! @param[inout] boundary     : boundary information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_processor_number(bound_info, table, pairlistdist, &
                                    water_model, ensemble, rigid_bond, boundary)

    ! formal arguments
    type(s_boundary_info),   intent(in)    :: bound_info
    logical,                 intent(in)    :: table
    real(wp),                intent(in)    :: pairlistdist
    character(*),            intent(in)    :: water_model
    integer,                 intent(in)    :: ensemble
    logical,                 intent(in)    :: rigid_bond
    type(s_boundary),        intent(inout) :: boundary

    ! local variable
    real(wp)                 :: bsize_x, bsize_y, bsize_z, buffer
    real(wp)                 :: size_x, size_y, size_z, maxsize(0:100)
    integer                  :: total_proc
    integer                  :: nx, ny, nz, nx1, ny1,nz1, i, j, k, itype
    integer                  :: nc(3,100), cell_size(3,100)
    logical                  :: defined_proc


    buffer  = bound_info%pbc_info%cell_size_buffer
    if (buffer < 0.0_wp) then
      if (ensemble == EnsembleNPT  .or. &
          ensemble == EnsembleNPAT .or. &
          ensemble == EnsembleNPgT) then
        buffer = 0.6_wp
      else
        buffer = 0.0_wp
      end if
    end if

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
        boundary%num_domain_max(1) = bound_info%domain_x_max
        boundary%num_domain_max(2) = bound_info%domain_y_max
        boundary%num_domain_max(3) = bound_info%domain_z_max

        if (boundary%num_domain_max(1) > 0 .and. &
            boundary%num_domain_max(2) > 0 .and. &
            boundary%num_domain_max(3) > 0) then

          if (boundary%num_domain_max(1) < boundary%num_domain(1))             &
            call error_msg('Setup_Processor_Number> domain_x_max should be '// &
                           'greater than domain_x')
          if (boundary%num_domain_max(2) < boundary%num_domain(2))             &
            call error_msg('Setup_Processor_Number> domain_y_max should be '// &
                           'greater than domain_y')
          if (boundary%num_domain_max(3) < boundary%num_domain(3))             &
            call error_msg('Setup_Processor_Number> domain_z_max should be '// &
                           'greater than domain_z')

          if (mod(boundary%num_domain_max(1),boundary%num_domain(1)) /= 0)     &
            call error_msg('Setup_Processor_Number> domain_x_max should be '// &
                           'multiple of domain_x')
          if (mod(boundary%num_domain_max(2),boundary%num_domain(2)) /= 0)     &
            call error_msg('Setup_Processor_Number> domain_y_max should be '// &
                           'multiple of domain_y')
          if (mod(boundary%num_domain_max(3),boundary%num_domain(3)) /= 0)     &
            call error_msg('Setup_Processor_Number> domain_z_max should be '// &
                           'multiple of domain_z')
        
        else

          boundary%num_domain_max(1:3) = boundary%num_domain(1:3)

        end if

        return

      else

        call error_msg('Setup_Processor_Number> # of process is not '// &
                       'domain_x * domain_y * domain_z ')

      end if

    end if

    bsize_x = boundary%box_size_x
    bsize_y = boundary%box_size_y
    bsize_z = boundary%box_size_z

    nx  = int (bsize_x / (pairlistdist+2.0_wp+buffer))
    ny  = int (bsize_y / (pairlistdist+2.0_wp+buffer))
    nz  = int (bsize_z / (pairlistdist+2.0_wp+buffer))
    nx1 = int (bsize_x / ((pairlistdist+2.0_wp+buffer)/2.0_wp))
    ny1 = int (bsize_y / ((pairlistdist+2.0_wp+buffer)/2.0_wp))
    nz1 = int (bsize_z / ((pairlistdist+2.0_wp+buffer)/2.0_wp))

    if (nx*ny*nz >= nproc_city) then

      itype = 0
      if (mod(nproc_city,8) == 0) then
        do k = 2, nz
          do j = 2, ny
            do i = 2, nx
              if (i*j*k == nproc_city) then
                itype = itype + 1
                nc(1,itype) = i
                nc(2,itype) = j
                nc(3,itype) = k
                cell_size(1,itype) = nx/i * i
                cell_size(2,itype) = ny/j * j
                cell_size(3,itype) = nz/k * k
              end if
            end do
          end do
        end do
      else if (nproc_city == 1) then
        itype = itype + 1
        nc(1,itype) = 1
        nc(2,itype) = 1
        nc(3,itype) = 1
        cell_size(1,itype) = nx
        cell_size(2,itype) = ny
        cell_size(3,itype) = nz
      else if (nproc_city == 2) then
        itype = itype + 1
        nc(1,itype) = 2
        nc(2,itype) = 1
        nc(3,itype) = 1
        cell_size(1,itype) = nx/2 * 2
        cell_size(2,itype) = ny
        cell_size(3,itype) = nz
      else if (nproc_city == 4) then
        itype = itype + 1
        nc(1,itype) = 2
        nc(2,itype) = 2
        nc(3,itype) = 1
        cell_size(1,itype) = nx/2 * 2
        cell_size(2,itype) = ny/2 * 2
        cell_size(3,itype) = nz
      else
        defined_proc = .false.
        do k = 2, nz
          do j = 2, ny
            do i = 2, nx
              if (i*j*k == nproc_city) then
                itype = itype + 1
                nc(1,itype) = i
                nc(2,itype) = j
                nc(3,itype) = k
                cell_size(1,itype) = nx/i * i
                cell_size(2,itype) = ny/j * j
                cell_size(3,itype) = nz/k * k
              end if
            end do
          end do
        end do
        if (.not. defined_proc) &
          call error_msg('Setup_Processor_Number> MPI Process number can not '//&
                         'be defined, please set them manualy')
      end if
    else
      call error_msg('Setup_Processor_Number> Cannot define domains'//         &
                     ' and cells. '//                                          &
                     'Smaller MPI processors, or shorter pairlistdist, or '//  &
                     'larger boxsize should be used.')
    end if

    if (itype == 0) then
      call error_msg('Setup_Processor_Number> Cannot define domains '//         &
                     'and cells. '//                                           &
                     'Smaller or adjusted MPI processors, or shorter'//        &
                     ' pairlistdist, or larger boxsize should be used.')
    end if

    k = 0
    maxsize(0) = 100000000000.0_wp
    do i = 1, itype
      size_x = bsize_x/cell_size(1,i)
      size_y = bsize_y/cell_size(2,i)
      size_z = bsize_z/cell_size(3,i)
      maxsize(i) = size_x*size_y*size_z
      if (maxsize(i) < maxsize(k)) &
        k = i
    end do

    boundary%num_domain(1) = nc(1,k)
    boundary%num_domain(2) = nc(2,k)
    boundary%num_domain(3) = nc(3,k)
    boundary%num_domain_max(1:3) = boundary%num_domain(1:3)

    return

  end subroutine setup_processor_number

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_processor_number_pio
  !> @brief        define the processor number in each dimension (parallel I/O)
  !! @authors      JJ
  !! @param[in]    bound_info : BOUNDARY section control parameters information
  !! @param[inout] boundary   : boundary information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_processor_number_pio(bound_info, boundary)

    ! formal arguments
    type(s_boundary_info), target,  intent(in)    :: bound_info
    type(s_boundary),      target,  intent(inout) :: boundary

    ! local variable
    real(wp)                 :: bsize_x, bsize_y, bsize_z
    real(wp)                 :: size_x, size_y, size_z
    integer                  :: total_proc
    integer                  :: multiple(3), iproc(3)
    logical                  :: extend, extend1
    logical                  :: defined_proc

    integer,         pointer :: num_domain(:), num_pio_domain(:)
    integer,         pointer :: domain_x, domain_y, domain_z
    integer,         pointer :: pio_domain_x, pio_domain_y, pio_domain_z


    num_domain     => boundary%num_domain
    num_pio_domain => boundary%num_pio_domain
    domain_x       => bound_info%domain_x
    domain_y       => bound_info%domain_y
    domain_z       => bound_info%domain_z
    pio_domain_x   => bound_info%pio_domain_x
    pio_domain_y   => bound_info%pio_domain_y
    pio_domain_z   => bound_info%pio_domain_z
    boundary%type  =  bound_info%type

    boundary%multiple_file(1:3) = 1

    if (domain_x /= 0 .and. domain_y /= 0 .and. domain_z /= 0) then

      total_proc = domain_x * domain_y * domain_z

      if (total_proc == nproc_country) then

        num_domain(1)     = domain_x
        num_domain(2)     = domain_y
        num_domain(3)     = domain_z
        num_pio_domain(1) = pio_domain_x
        num_pio_domain(2) = pio_domain_y
        num_pio_domain(3) = pio_domain_z

        if (num_pio_domain(1) /= 0 .or. num_pio_domain(2) /= 0 .or. &
            num_pio_domain(3) /= 0) then

          if (num_domain(1) >= num_pio_domain(1) .and. &
              num_domain(2) >= num_pio_domain(2) .and. &
              num_domain(3) >= num_pio_domain(3)) then

            ! num_domain should be the multiple of num_pio_domain
            !
            if (mod(num_domain(1), num_pio_domain(1)) /= 0)                   &
              call error_msg('Setup_Processor_Number_pio> domain_x should '// &
                             'be multiple of pio_domain_x')
            if (mod(num_domain(2), num_pio_domain(2)) /= 0)                   &
              call error_msg('Setup_Processor_Number_pio> domain_y should '// &
                             'be multiple of pio_domain_y')
            if (mod(num_domain(3), num_pio_domain(3)) /= 0)                   &
              call error_msg('Setup_Processor_Number_pio> domain_z should '// &
                             'be multiple of pio_domain_z')
  
            ! decide the rank of file to be open
            !
            multiple(1:3) = num_domain(1:3) / num_pio_domain(1:3)
            iproc(1)      = mod(my_country_rank, num_domain(1))
            iproc(2)      = mod(my_country_rank/num_domain(1), num_domain(2))
            iproc(3)      = my_country_rank/(num_domain(1)*num_domain(2))
            iproc(1:3) = iproc(1:3) / multiple(1:3)
            my_rank_pio = iproc(1) + iproc(2)*num_pio_domain(1) &
                        + iproc(3)*num_pio_domain(1)*num_pio_domain(2)

          else if (num_domain(1) <= num_pio_domain(1) .and. &
                   num_domain(2) <= num_pio_domain(2) .and. &
                   num_domain(3) <= num_pio_domain(3)) then

            ! num_pio_domain should be multiple of num_domain
            !
            if (mod(num_pio_domain(1), num_domain(1)) /= 0)                &
              call error_msg('Setup_Processor_Number_pio> pio_domain_x '// &
                             'should be multiple of domain_x')
            if (mod(num_pio_domain(2), num_domain(2)) /= 0)                &
              call error_msg('Setup_Processor_Number_pio> pio_domain_y '// &
                             'should be multiple of domain_y')
            if (mod(num_pio_domain(3), num_domain(3)) /= 0)                &
              call error_msg('Setup_Processor_Number_pio> pio_domain_z '// &
                             'should be multiple of domain_z')

            boundary%multiple_file_read = .true.
            boundary%multiple_file(1:3) = num_pio_domain(1:3) / num_domain(1:3)

          else
            call error_msg('Setup_Processor_Number_pio> pio_domain '// &
                           'should be greater/less than domain in '//  &
                           'all directoins')

          end if

        else

          ! pio_domain_[x,y,z] is not written => same number of MPIs
          !
          my_rank_pio = my_country_rank

        end if
 
      else

        call error_msg('Setup_Processor_Number> # of process is not '// &
                       'domain_x * domain_y * domain_z ')

      end if

    else

      ! domain_[x,y,z] is not written => same number of MPIs
      !
      my_rank_pio = my_country_rank
 
    end if

    return

  end subroutine setup_processor_number_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_boundary_cell
  !> @brief        setup boundary cell information
  !! @authors      NT
  !! @param[in]    buffer       : cell size buffer
  !! @param[in]    table        : flag for use table or not
  !! @param[in]    pairlistdist : pair-list distance
  !! @pmara[in]    water_model  : water model
  !! @param[in]    ensemble     : type of ensemble 
  !! @param[in]    rigid_bond   : flag for rigid-bond
  !! @param[in]    dsize_cg     : flag for reset domain size for CG-model
  !! @param[in]    dmin_size_cg : minimum domain size for CG-model
  !! @param[inout] boundary     : boundary information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_boundary_cell(buffer, table, pairlistdist,        &
                                 water_model, ensemble, rigid_bond,  &
                                 dsize_cg, dmin_size_cg,  boundary)

    ! formal arguments
    real(wp),                 intent(in)    :: buffer
    logical,                  intent(in)    :: table
    real(wp),                 intent(in)    :: pairlistdist
    character(*),             intent(in)    :: water_model
    integer,                  intent(in)    :: ensemble
    logical,                  intent(in)    :: rigid_bond
    logical,                  intent(in)    :: dsize_cg
    real(wp),                 intent(in)    :: dmin_size_cg
    type(s_boundary), target, intent(inout) :: boundary

    ! local variables
    real(wp)                 :: csize_x, csize_y, csize_z, cutoff
    real(wp)                 :: bsize_x, bsize_y, bsize_z
    real(wp)                 :: cell_size_buffer
    integer                  :: k, ncell
    integer                  :: ncell_x, ncell_y, ncell_z
    logical                  :: extend, extend1

    integer, pointer         :: num_domain(:), num_domain_max(:)


    num_domain     => boundary%num_domain
    num_domain_max => boundary%num_domain_max

    cell_size_buffer = buffer
    if (cell_size_buffer < 0.0_wp) then
      if (ensemble == EnsembleNPT  .or. &
          ensemble == EnsembleNPAT .or. &
          ensemble == EnsembleNPgT) then
        cell_size_buffer = 0.6_wp
      else
        cell_size_buffer = 0.0_wp
      end if
    end if

    cutoff = pairlistdist + 2.0_wp + cell_size_buffer

    if (dsize_cg) then
      cutoff=max(dmin_size_cg,pairlistdist)
    end if

    bsize_x = boundary%box_size_x
    bsize_y = boundary%box_size_y
    bsize_z = boundary%box_size_z

    if ((boundary%num_cells_x /= 0) .and. &
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
      write(MsgOut,'(a,f15.8)')  'Debuging > cutoff', cutoff
      write(MsgOut,'(a,3f15.8)') 'Debuging > bsize_[x,y,z]', &
                                 bsize_x, bsize_y, bsize_z
      write(MsgOut,'(a,3i8)')    'Debuging > ncell_[x,y,z]', &
                                 ncell_x, ncell_y, ncell_z
    end if
#endif

    k = mod(ncell_x, num_domain_max(1))
    if (k /= 0) ncell_x = ncell_x - k

    k = mod(ncell_y, num_domain_max(2))
    if (k /= 0) ncell_y = ncell_y - k

    k = mod(ncell_z, num_domain_max(3))
    if (k /= 0) ncell_z = ncell_z - k

    if (ncell_x < 5 .or. ncell_y < 5 .or. ncell_z < 5) &
      call error_msg('Setup_Boundary_Cell> too small boxsize/pairlistdist. '//&
                     'shorter pairlistdist or larger boxsize or less MPI processors'//&
                     ' should be used.')

    csize_x = bsize_x/real(ncell_x, wp)
    csize_y = bsize_y/real(ncell_y, wp)
    csize_z = bsize_z/real(ncell_z, wp)
    ncell   = ncell_x*ncell_y*ncell_z

    boundary%num_cells_x = ncell_x
    boundary%num_cells_y = ncell_y
    boundary%num_cells_z = ncell_z
    boundary%cell_size_x = csize_x
    boundary%cell_size_y = csize_y
    boundary%cell_size_z = csize_z

    if (main_rank) then
      write(MsgOut,'(A)') &
           'Setup_Boundary_Cell> Set Variables for Boundary Condition'

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
      call error_msg( &
          'Setup_Boundary_Cell> ncell_[x,y,z] should be greater than or equal to '//&
          '2*domain_[x,y,z]. Please reduce MPI and increase OpenMP to use the '//   &
          'same number of processors.')

    return

  end subroutine setup_boundary_cell

end module sp_boundary_mod
