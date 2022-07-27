!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sa_boundary_mod
!> @brief   utilities for boundary conditions
!! @authors Jaewoon Jung (JJ), Norio Takase (NT), Isseki Yu (IY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module sa_boundary_mod

  use sa_ensemble_str_mod
  use sa_boundary_str_mod
  use sa_option_str_mod
  use trajectory_str_mod
  use molecules_str_mod
  use fileio_control_mod
  use fileio_trj_mod
  use fileio_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_pbc_info
    real(wp)            :: box_size_x  = 0.0_wp
    real(wp)            :: box_size_y  = 0.0_wp
    real(wp)            :: box_size_z  = 0.0_wp
    integer             :: num_cells_x = 0
    integer             :: num_cells_y = 0
    integer             :: num_cells_z = 0
  end type s_pbc_info

  type, public :: s_boundary_info
    integer             :: type        = BoundaryTypePBC
    type(s_pbc_info)    :: pbc_info
    real(wp)            :: origin_x    = 0.0_wp
    real(wp)            :: origin_y    = 0.0_wp
    real(wp)            :: origin_z    = 0.0_wp
    integer             :: domain_x    = 0
    integer             :: domain_y    = 0
    integer             :: domain_z    = 0
  end type s_boundary_info

  ! subroutines
  public  :: show_ctrl_boundary
  public  :: read_ctrl_boundary
  public  :: setup_boundary
  public  :: setup_processor_number
  public  :: setup_boundary_cell
  public  :: refresh_boundary

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_boundary
  !> @brief        show BOUNDARY section usage
  !! @authors      IY, NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine show_ctrl_boundary

    write(MsgOut,'(A)') '[BOUNDARY]'
    write(MsgOut,'(A)') 'type          = PBC       # [PBC]'
    write(MsgOut,'(A)') 'box_size_x    = 0.0       # box size (x) in [PBC]'
    write(MsgOut,'(A)') 'box_size_y    = 0.0       # box size (y) in [PBC]'
    write(MsgOut,'(A)') 'box_size_z    = 0.0       # box size (z) in [PBC]'
    write(MsgOut,'(A)') 'domain_x      = 0         # domain size (x)'
    write(MsgOut,'(A)') 'domain_y      = 0         # domain size (y)'
    write(MsgOut,'(A)') 'domain_z      = 0         # domain size (z)'
    write(MsgOut,'(A)') 'num_cells_x   = 0         # number of cells (x)'
    write(MsgOut,'(A)') 'num_cells_y   = 0         # number of cells (x)'
    write(MsgOut,'(A)') 'num_cells_z   = 0         # number of cells (x)'
    write(MsgOut,'(A)') ' '

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

    ! should be remove case PBC
    ! currently box_size is read from .dcd (IY)
    select case (bound_info%type)

    case (BoundaryTypePBC)
      call read_ctrlfile_real (handle, Section, 'box_size_x',  &
                               bound_info%pbc_info%box_size_x)
      call read_ctrlfile_real (handle, Section, 'box_size_y',  &
                               bound_info%pbc_info%box_size_y)
      call read_ctrlfile_real (handle, Section, 'box_size_z',  &
                               bound_info%pbc_info%box_size_z)

      ! determin cell numbers
      call read_ctrlfile_integer (handle, Section, 'num_cells_x',  &
                               bound_info%pbc_info%num_cells_x)
      call read_ctrlfile_integer (handle, Section, 'num_cells_y',  &
                               bound_info%pbc_info%num_cells_y)
      call read_ctrlfile_integer (handle, Section, 'num_cells_z',  &
                               bound_info%pbc_info%num_cells_z)

    end select

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

      if (bound_info%domain_x .ne. 0 .and. &
          bound_info%domain_y .ne. 0 .and. &
          bound_info%domain_z .ne. 0) &
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
  !> @brief        set essential variables for boundary condition
  !! @authors      JJ, IY
  !! @param[in]    bound_info   : BOUNDARY section control parameters information
  !! @param[in]    bound_extend : extend the boundary or not
  !! @param[in]    pairlistdist : pair-list distance
  !! @param[in]    ensemble     : type of ensemble
  !! @param[in]    rigid_bond   : flag for rigid-bond
  !! @param[inout] boundary     : boundary information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_boundary(bound_info, bound_extend, pairlistdist, ensemble, &
                            molecule, boundary, trj_list, trajectory, option)

    ! formal arguments
    type(s_boundary_info),   intent(in)    :: bound_info
    logical,                 intent(in)    :: bound_extend
    real(wp),                intent(in)    :: pairlistdist
    integer,                 intent(in)    :: ensemble
    type(s_molecule),        intent(in)    :: molecule
    type(s_boundary),        intent(inout) :: boundary
    type(s_trajectory),      intent(inout) :: trajectory
    type(s_trj_list),        intent(in)    :: trj_list
    type(s_option),          intent(in)    :: option

    ! local variable
    integer                  :: i
    real(wp)                 :: coord_min(1:3), coord_max(1:3), box_size(1:3)
    type(s_trj_file)         :: trj_in


    select case (bound_info%type)

    case (BoundaryTypeNOBC)

     if (option%determine_box == DetermineBoxMax) then

       boundary%type     = bound_info%type
       boundary%origin_x = bound_info%origin_x
       boundary%origin_y = bound_info%origin_y
       boundary%origin_z = bound_info%origin_z

       ! Decide system size
       !
       coord_min(1:3) =  1000000.0_wp
       coord_max(1:3) = -1000000.0_wp
       do i = 1, molecule%num_atoms
        coord_min(1:3) = min(coord_min(1:3), molecule%atom_coord(1:3,i))
        coord_max(1:3) = max(coord_max(1:3), molecule%atom_coord(1:3,i))
       end do

       box_size(1:3) = max(-coord_min(1:3) + 0.1_wp, coord_max(1:3 ) + 0.1_wp)
       boundary%box_size_x     = box_size(1) * 2.0_wp
       boundary%box_size_y     = box_size(2) * 2.0_wp
       boundary%box_size_z     = box_size(3) * 2.0_wp

     else if (option%determine_box == DetermineboxManual) then

       boundary%box_size_x     = bound_info%pbc_info%box_size_x
       boundary%box_size_y     = bound_info%pbc_info%box_size_y
       boundary%box_size_z     = bound_info%pbc_info%box_size_z

     end if

    case (BoundaryTypePBC)

      if (option%determine_box == DetermineBoxTrajectory) then
        !--------read box size from trajectory (IY)
        call open_trj(trj_in, trj_list%filenames(1), &
                              trj_list%trj_format,   &
                              trj_list%trj_type, IOFileInput)

        call read_trj(trj_in, trajectory)

        boundary%type           = bound_info%type
        boundary%origin_x       = bound_info%origin_x
        boundary%origin_y       = bound_info%origin_y
        boundary%origin_z       = bound_info%origin_z
        boundary%box_size_x     = trajectory%pbc_box(1,1)
        boundary%box_size_y     = trajectory%pbc_box(2,2)
        boundary%box_size_z     = trajectory%pbc_box(3,3)

        call close_trj(trj_in)

      else if (option%determine_box == DetermineBoxManual) then

        boundary%type           = bound_info%type
        boundary%box_size_x     = bound_info%pbc_info%box_size_x
        boundary%box_size_y     = bound_info%pbc_info%box_size_y
        boundary%box_size_z     = bound_info%pbc_info%box_size_z

       end if

    end select

    !-------read cell size from control (IY)
    boundary%num_cells_x    = bound_info%pbc_info%num_cells_x
    boundary%num_cells_y    = bound_info%pbc_info%num_cells_y
    boundary%num_cells_z    = bound_info%pbc_info%num_cells_z

    call setup_processor_number(bound_info, bound_extend, &
                                pairlistdist, ensemble, boundary)

    call setup_boundary_cell   (bound_extend, pairlistdist, ensemble, boundary)

    return

  end subroutine setup_boundary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !  Subroutine    setup_processor_number
  !> @brief        define the processor number in each dimension
  !! @authors      JJ
  !! @param[in]    bound_info : BOUNDARY section control parameters information
  !! @param[in]    bound_extend : extend the boundary or not
  !! @param[in]    pairlistdist : pair-list distance
  !! @param[in]    ensemble     : type of ensemble
  !! @param[inout] boundary     : boundary information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_processor_number(bound_info, bound_extend,   &
                                    pairlistdist, ensemble, boundary)

    ! formal arguments
    type(s_boundary_info),   intent(in)    :: bound_info
    logical,                 intent(in)    :: bound_extend
    real(wp),                intent(in)    :: pairlistdist
    integer,                 intent(in)    :: ensemble
    type(s_boundary),        intent(inout) :: boundary

    ! local variable
    real(wp)                 :: bsize_x, bsize_y, bsize_z
    real(wp)                 :: size_x, size_y, size_z, maxsize(0:100)
    integer                  :: total_proc
    integer                  :: nx, ny, nz, nx1, ny1,nz1, i, j, k, itype
    integer                  :: nc(3,100), cell_size(3,100)
    logical                  :: extend1


    if (ensemble == EnsembleNPT  .or. &
        ensemble == EnsembleNPAT .or. &
        ensemble == EnsembleNPgT)  then
      extend1 = .true.

    else
      extend1 = .false.
    end if

    ! check processor number based on the num of domains
    !
    if (bound_info%domain_x /= 0 .and. &
        bound_info%domain_y /= 0 .and. &
        bound_info%domain_z /= 0) then

      total_proc = bound_info%domain_x * &
                   bound_info%domain_y * &
                   bound_info%domain_z

    end if

    bsize_x = boundary%box_size_x
    bsize_y = boundary%box_size_y
    bsize_z = boundary%box_size_z

    if (bound_extend) then
      if (extend1) then
        nx  = bsize_x / (pairlistdist + 2.6_wp)
        ny  = bsize_y / (pairlistdist + 2.6_wp)
        nz  = bsize_z / (pairlistdist + 2.6_wp)
        nx1 = bsize_x / ((pairlistdist + 2.6_wp) * half)
        ny1 = bsize_y / ((pairlistdist + 2.6_wp) * half)
        nz1 = bsize_z / ((pairlistdist + 2.6_wp) * half)
      else
        nx  = bsize_x / (pairlistdist + 2.0_wp)
        ny  = bsize_y / (pairlistdist + 2.0_wp)
        nz  = bsize_z / (pairlistdist + 2.0_wp)
        nx1 = bsize_x / ((pairlistdist + 2.0_wp) * half)
        ny1 = bsize_y / ((pairlistdist + 2.0_wp) * half)
        nz1 = bsize_z / ((pairlistdist + 2.0_wp) * half)
      end if
    else
      if (extend1) then
        nx  = bsize_x / (pairlistdist + 0.6_wp)
        ny  = bsize_y / (pairlistdist + 0.6_wp)
        nz  = bsize_z / (pairlistdist + 0.6_wp)
        nx1 = bsize_x / ((pairlistdist + 0.6_wp) * half)
        ny1 = bsize_y / ((pairlistdist + 0.6_wp) * half)
        nz1 = bsize_z / ((pairlistdist + 2.6_wp) * half)
      else
        nx  = bsize_x / pairlistdist
        ny  = bsize_y / pairlistdist
        nz  = bsize_z / pairlistdist
        nx1 = 2.0_wp * bsize_x / pairlistdist
        ny1 = 2.0_wp * bsize_y / pairlistdist
        nz1 = 2.0_wp * bsize_z / pairlistdist
      end if
    end if

    if (nx*ny*nz >= total_proc) then

      itype = 0
      if (mod(nproc_city,8) == 0) then
        do k = 2, nz
          do j = 2, ny
            do i = 2, nx
              if (i*j*k == total_proc) then
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
      else if (total_proc == 2) then
        itype = itype + 1
        nc(1,itype) = 2
        nc(2,itype) = 1
        nc(3,itype) = 1
        cell_size(1,itype) = nx/2 * 2
        cell_size(2,itype) = ny
        cell_size(3,itype) = nz
      else if (total_proc == 4) then
        itype = itype + 1
        nc(1,itype) = 2
        nc(2,itype) = 2
        nc(3,itype) = 1
        cell_size(1,itype) = nx/2 * 2
        cell_size(2,itype) = ny/2 * 2
        cell_size(3,itype) = nz
      else
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              if (i*j*k == total_proc) then
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
      end if
    else
      call error_msg('Setup_Processor_Number> MPI Process number should be smaller')
    end if

    if (itype == 0) then
      call error_msg('Setup_Processor_Number> MPI Process number should be smaller or adjusted')
    end if

    k = 0
    maxsize(0) = 100000000000.0_wp
    do i = 1, itype
      size_x = bsize_x / cell_size(1,i)
      size_y = bsize_y / cell_size(2,i)
      size_z = bsize_z / cell_size(3,i)
      maxsize(i) = size_x*size_y*size_z
      if (maxsize(i) < maxsize(k)) &
        k = i
    end do

     !    boundary%num_domain(1) = nc(1,k)
     !    boundary%num_domain(2) = nc(2,k)
     !    boundary%num_domain(3) = nc(3,k)

     ! fix the values of num_domain by control file (IY)
     boundary%num_domain(1) = bound_info%domain_x
     boundary%num_domain(2) = bound_info%domain_y
     boundary%num_domain(3) = bound_info%domain_z

    return

  end subroutine setup_processor_number

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_boundary_cell
  !> @brief        setup boundary cell information
  !! @authors      NT
  !! @param[in]    bound_extend : extend the boundary or not
  !! @param[in]    pairlistdist : pair-list distance
  !! @param[in]    ensemble     : type of ensemble
  !! @param[inout] boundary     : boundary information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_boundary_cell(bound_extend, pairlistdist, ensemble, boundary)

    ! formal arguments
    logical,                 intent(in)    :: bound_extend
    real(wp),                intent(in)    :: pairlistdist
    integer,                 intent(in)    :: ensemble
    type(s_boundary),        intent(inout) :: boundary

    ! local variables
    real(wp)                 :: csize_x, csize_y, csize_z, cutoff
    real(wp)                 :: bsize_x, bsize_y, bsize_z
    integer                  :: k, ncell
    integer                  :: ncell_x, ncell_y, ncell_z
    logical                  :: extend1


    if (ensemble == EnsembleNPT  .or. &
        ensemble == EnsembleNPAT .or. &
        ensemble == EnsembleNPgT)  then
      extend1 = .true.
    else
      extend1 = .false.
    end if

    if (bound_extend) then
      if (extend1) then
        cutoff = pairlistdist + 2.6_wp
      else
        cutoff = pairlistdist + 2.0_wp
      end if
    else
      if (extend1) then
        cutoff = pairlistdist + 0.6_wp
      else
        cutoff = pairlistdist
      end if
    end if

    bsize_x = boundary%box_size_x
    bsize_y = boundary%box_size_y
    bsize_z = boundary%box_size_z

    if ( (boundary%num_cells_x /= 0) .and. &
         (boundary%num_cells_y /= 0) .and. &
         (boundary%num_cells_z /= 0)) then
      ncell_x = boundary%num_cells_x
      ncell_y = boundary%num_cells_y
      ncell_z = boundary%num_cells_z
    else
      ncell_x = int(bsize_x / (cutoff * half))
      ncell_y = int(bsize_y / (cutoff * half))
      ncell_z = int(bsize_z / (cutoff * half))
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

    k = mod(ncell_x, boundary%num_domain(1))
    if (k .ne. 0) &
      ncell_x = ncell_x - k

    k = mod(ncell_y, boundary%num_domain(2))
    if (k .ne. 0) &
      ncell_y = ncell_y - k

    k = mod(ncell_z, boundary%num_domain(3))
    if (k .ne. 0) &
      ncell_z = ncell_z - k

    if (ncell_x < 5 .or. ncell_y < 5 .or. ncell_z < 5) &
      call error_msg('Setup_Boundary_Cell> too small boxsize/cutoff. '//&
                   'shorter cutoff or larger boxsize or less MPI processors'//&
                   ' should be used.')

    csize_x = bsize_x / real(ncell_x, wp)
    csize_y = bsize_y / real(ncell_y, wp)
    csize_z = bsize_z / real(ncell_z, wp)
    ncell   = ncell_x * ncell_y * ncell_z

    boundary%num_cells_x = ncell_x
    boundary%num_cells_y = ncell_y
    boundary%num_cells_z = ncell_z
    boundary%cell_size_x = csize_x
    boundary%cell_size_y = csize_y
    boundary%cell_size_z = csize_z
    ! prepare cell neighbor list
    !

!    call alloc_boundary(boundary, BoundaryCells, ncell)


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

    return

  end subroutine setup_boundary_cell

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    refresh_boundary
  !> @brief        refresh box size from trajectory
  !! @authors      IY
  !! @param[in]    trajectory : trajectory information
  !! @param[inout] boundary   : boundary information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine refresh_boundary(trajectory, boundary)

    type(s_trajectory),      intent(in)    :: trajectory
    type(s_boundary),        intent(inout) :: boundary

    real(wp)                 :: bsize_x, bsize_y, bsize_z
    real(wp)                 :: csize_x, csize_y, csize_z
    integer                  :: ncell_x, ncell_y, ncell_z


    boundary%box_size_x = trajectory%pbc_box(1,1)
    boundary%box_size_y = trajectory%pbc_box(2,2)
    boundary%box_size_z = trajectory%pbc_box(3,3)

    bsize_x = boundary%box_size_x
    bsize_y = boundary%box_size_y
    bsize_z = boundary%box_size_z

    ncell_x = boundary%num_cells_x
    ncell_y = boundary%num_cells_y
    ncell_z = boundary%num_cells_z

    csize_x = bsize_x / real(ncell_x, wp)
    csize_y = bsize_y / real(ncell_y, wp)
    csize_z = bsize_z / real(ncell_z, wp)

    boundary%cell_size_x = csize_x
    boundary%cell_size_y = csize_y
    boundary%cell_size_z = csize_z

    return

  end subroutine refresh_boundary

end module sa_boundary_mod
