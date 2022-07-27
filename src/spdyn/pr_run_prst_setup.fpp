!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pr_run_prst_setup_mod
!> @brief   run output parallel I/O restart files
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module pr_run_prst_setup_mod

  use pr_setup_spdyn_peer_mod
  use pr_control_mod
  use pr_domain_index_mod
  use pr_gromacs2hm_mod
  use pr_amber2hm_mod
  use pr_charmm2hm_mod
  use pr_huge_molecule_mod
  use sp_domain_mod
  use sp_parallel_io_mod
  use sp_dynamics_str_mod
  use sp_dynvars_str_mod
  use sp_restraints_str_mod
  use sp_constraints_str_mod
  use sp_boundary_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use fileio_grotop_mod
  use fileio_prmtop_mod
  use fileio_str_mod
  use fileio_top_mod
  use fileio_par_mod
  use fileio_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use timers_mod
  use mpi
 
  implicit none
  private

  ! subroutines
  public  :: run_prst_setup
  private :: setup_mpi_dummy
  private :: setup_huge_molecule
  private :: setup_molecule_variables
  private :: setup_domain_variables
  private :: setup_domain_rst_files
  private :: setup_domain_index_list
  private :: output_domain_restart_files
  private :: output_pdb_file
  private :: setup_the

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_prst_setup
  !> @brief        run output parallel I/O rstart files
  !! @authors      NT, JJ
  !! @param[in]    ctrl_data : information of control parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_prst_setup(ctrl_data)

    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data

    ! local variables
    type(s_boundary)         :: boundary
    type(s_boundary)         :: boundary_rst
    type(s_constraints)      :: constraints
    type(s_restraints)       :: restraints
    type(s_enefunc)          :: enefunc
    type(s_dynvars)          :: dynvars
    type(s_dynamics)         :: dynamics
    type(s_par)              :: par
    type(s_prmtop)           :: prmtop
    type(s_grotop)           :: grotop
    real(wp)                 :: ratio
    logical                  :: block_domain


    ! setup domain information
    !
    if (main_rank) then
      write(MsgOut,'(A)') '[STEP2-1] Setup MPI dummy'
      write(MsgOut,'(A)') ' '
    end if

    block_domain = .false.
    call setup_mpi_dummy(block_domain, ctrl_data)


    ! setup huge molecule data
    !
    if (main_rank) then
      write(MsgOut,'(A)') '[STEP2-2] Setup huge molecule data'
      write(MsgOut,'(A)') ' '
    end if

    call setup_huge_molecule(ctrl_data, par, prmtop, grotop)


    ! setup molecule variables
    !
    if (main_rank) then
      write(MsgOut,'(A)') '[STEP2-3] Setup molecule variables'
      write(MsgOut,'(A)') ' '
    end if

    call timer(TimerTest8, TimerOn)
    call setup_molecule_variables(ctrl_data)
    call timer(TimerTest8, TimerOff)


    ! setup domain restart files
    !
    if (main_rank) then
      write(MsgOut,'(A)') '[STEP2-4] Setup domain restart files'
      write(MsgOut,'(A)') ' '
    end if

    call setup_domain_rst_files(ctrl_data, boundary_rst)


    ! setup domain variables (once)
    !
    if (main_rank) then
      write(MsgOut,'(A)') '[STEP2-5] Setup domain variables'
      write(MsgOut,'(A)') ' '
    end if

    call timer(TimerTest9, TimerOn)
    call setup_domain_variables(ctrl_data, par, prmtop, grotop, &
                                boundary_rst, boundary,         &
                                enefunc, dynvars, dynamics,     &
                                restraints, constraints)
    call timer(TimerTest9, TimerOff)


    ! output PDB file
    !
    if (main_rank) then
      write(MsgOut,'(A)') '[STEP2-6] Output PDB file'
      write(MsgOut,'(A)') ' '
    end if

    call output_pdb_file(ctrl_data)


    ! setup domain index list
    !
    if (main_rank) then
      write(MsgOut,'(A)') '[STEP2-7] Setup domain index list'
      write(MsgOut,'(A)') ' '
    end if

    call setup_domain_index_list(ctrl_data, boundary)


    ! output domain restart files
    !
    if (main_rank) then
      write(MsgOut,'(A)') '[STEP2-8] Output domain restart files'
      write(MsgOut,'(A)') ' '
    end if

    call timer(TimerTest10, TimerOn)
    call output_domain_restart_files(block_domain, ctrl_data, &
                                     par, prmtop, grotop, &
                                     boundary, dynvars, dynamics,    &
                                     restraints, enefunc, constraints)
    call timer(TimerTest10, TimerOff)


    ! deallocate memory
    !
    if (main_rank) then
      write(MsgOut,'(A)') '[STEP2-9] Deallocate memory'
      write(MsgOut,'(A)') ' '
    end if

    call di_delete
    call hm_finalize

    call dealloc_restraints(restraints, RestraintsList)
    call dealloc_constraints(constraints, ConstraintsBondGroup)
    call dealloc_enefunc(enefunc, EneFuncReff)
    call dealloc_grotop_all(grotop)
    call dealloc_prmtop_all(prmtop)
    call dealloc_par_all(par)

    return

  end subroutine run_prst_setup

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_mpi_dummy
  !> @brief        setup MPI in dummy
  !! @authors      NT
  !! @param[in]    ctrl_data : information of control parameters
  !! @param[inout] boundary  : boundary condition information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_mpi_dummy(block_domain, ctrl_data)

    ! formal arguments
    logical,                 intent(inout) :: block_domain
    type(s_ctrl_data),       intent(in)    :: ctrl_data


    my_world_rank = 0
    my_city_rank  = 0

    if (ctrl_data%opt_info%domain_xyz /= 0) then
      nproc_country = ctrl_data%opt_info%domain_xyz 
    else
      nproc_country = ctrl_data%bou_info%domain_x * &
                      ctrl_data%bou_info%domain_y * &
                      ctrl_data%bou_info%domain_z
    end if
    nproc_city = nproc_country
 
    nproc_world = nproc_country

    return

  end subroutine setup_mpi_dummy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_huge_molecule
  !> @brief        setup huge molecule data
  !! @authors      NT
  !! @param[in]    ctrl_data : information of control parameters
  !! @param[inout] par       : CHARMM par information
  !! @param[inout] prmtop    : AMBER parameter / topology information
  !! @param[inout] grotop    : GROMACS topology information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_huge_molecule(ctrl_data, par, prmtop, grotop)

    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_par),             intent(inout) :: par
    type(s_prmtop),          intent(inout) :: prmtop
    type(s_grotop),          intent(inout) :: grotop

    ! local variables
    type(s_top)              :: top
    integer                  :: ispc, ierr
    character(200)           :: file, tfile, prefix, reffile
    character(200)           :: psffile, pdbfile, crdfile
    character(200)           :: prmtopfile, ambcrdfile, ambreffile
    character(200)           :: grotopfile, grocrdfile, groreffile

    ! initialization

    call init_top(top)
    call init_par(par)
    call init_prmtop(prmtop)
    call init_grotop(grotop)

    ! CHARMM input file
    !
    if (ctrl_data%inp_info%parfile /= '' .and. &
        ctrl_data%inp_info%topfile /= '' .and. &
        ctrl_data%inp_info%psffile /= '') then

      call input_top(ctrl_data%inp_info%topfile, top)
      call input_par(ctrl_data%inp_info%parfile, par)

      if (ctrl_data%inp_info%strfile /= '') &
      call input_str(ctrl_data%inp_info%strfile, top, par)
      

      ! setup prefix
      psffile = ctrl_data%inp_info%psffile
      pdbfile = ctrl_data%inp_info%pdbfile
      crdfile = ctrl_data%inp_info%crdfile
      reffile = ctrl_data%inp_info%reffile

      if (pdbfile /= '') then
        file = pdbfile
      else if (crdfile /= '') then
        file = crdfile
      else
        call error_msg('Setup_Huge_Molecule> coordinate is empty.')
      end if

      prefix = trim(psffile(scan(psffile,'/',.true.)+1:)) // '_' // &
               trim(file   (scan(file,   '/',.true.)+1:))


      ! create huge molecule data
      call hm_initialize(ctrl_data%opt_info%cachepath, prefix)

      call hm_create_charmm(top, par, psffile, pdbfile, crdfile, reffile)

      call dealloc_top_all(top)


    ! AMBER input file
    !
    else if (ctrl_data%inp_info%prmtopfile /= '') then

      call input_prmtop(ctrl_data%inp_info%prmtopfile, prmtop)


      ! setup prefix
      prmtopfile = ctrl_data%inp_info%prmtopfile
      ambcrdfile = ctrl_data%inp_info%ambcrdfile
      ambreffile = ctrl_data%inp_info%ambreffile
      pdbfile    = ctrl_data%inp_info%pdbfile
      reffile    = ctrl_data%inp_info%reffile

      if (ambcrdfile /= '') then
        file = ambcrdfile
      else if (pdbfile /= '') then
        file = pdbfile
      else
        call error_msg('Setup_Huge_Molecule> coordinate is empty.')
      end if


      do ispc = 1, len(prmtopfile)
        if (prmtopfile(ispc:ispc) == ' ' .or. &
            prmtopfile(ispc:ispc) == char(9)) &
          exit
      end do

      tfile = prmtopfile(:ispc-1)

      prefix = trim(tfile(scan(tfile,'/',.true.)+1:)) // '_' // &
               trim(file (scan(file, '/',.true.)+1:))


      ! create huge molecule data
      call hm_initialize(ctrl_data%opt_info%cachepath, prefix)

      call hm_create_amber(prmtop, ambcrdfile, ambreffile, pdbfile, reffile)


    ! GROMACS input file
    !
    else if (ctrl_data%inp_info%grotopfile /= '') then

      call input_grotop(ctrl_data%inp_info%grotopfile, grotop)


      ! setup prefix
      grotopfile = ctrl_data%inp_info%grotopfile
      grocrdfile = ctrl_data%inp_info%grocrdfile
      groreffile = ctrl_data%inp_info%groreffile
      pdbfile    = ctrl_data%inp_info%pdbfile
      reffile    = ctrl_data%inp_info%reffile
      
      if (grocrdfile /= '') then
        file = grocrdfile
      else if (pdbfile /= '') then
        file = pdbfile
      else
        call error_msg('Setup_Huge_Molecule> coordinate is empty.')
      end if


      do ispc = 1, len(grotopfile)
        if (grotopfile(ispc:ispc) == ' ' .or. &
            grotopfile(ispc:ispc) == char(9)) &
          exit
      end do

      tfile = grotopfile(:ispc-1)

      prefix = trim(tfile(scan(tfile,'/',.true.)+1:)) // '_' // &
               trim(file (scan(file, '/',.true.)+1:))


      ! create huge molecule data
      call hm_initialize(ctrl_data%opt_info%cachepath, prefix)

      call hm_create_gromacs(grotop, grocrdfile, groreffile, pdbfile, reffile)

    end if

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    return

  end subroutine setup_huge_molecule

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_molecule_variables
  !> @brief        setup molecule variables
  !! @authors      NT
  !! @param[in]    ctrl_data : information of control parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_molecule_variables(ctrl_data)

    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data

    ! local variables
    integer                  :: i, i1, i2, j, col, pcol, nmol, ierr
    integer                  :: src_col, dst_col, cur_col
    integer, allocatable     :: color(:), color_range(:,:)


    ! count the number of molecules
    !

    allocate(color(hm_num_atoms), color_range(2,hm_num_atoms/2))
    color(:)         = 0
    color_range(:,:) = 0
    
    ! coloring 
    cur_col = 0
    do i = 1, hm_num_bonds

      i1 = hm_bond_list(1,i)
      i2 = hm_bond_list(2,i)

      if (color(i1) == 0 .and. color(i2) == 0) then

        cur_col = cur_col + 1
        color(i1) = cur_col
        color(i2) = cur_col

        color_range(1,cur_col) = min(i1, i2)
        color_range(2,cur_col) = max(i1, i2)

      else if (color(i1) == 0 .or. color(i2) == 0) then

        if (color(i1) /= 0) then
          
          col = color(i1)
          color(i2) = col
          color_range(1,col) = min(color_range(1,col),i2)
          color_range(2,col) = max(color_range(2,col),i2)

        else ! color(i2) /= 0

          col = color(i2)
          color(i1) = col
          color_range(1,col) = min(color_range(1,col),i1)
          color_range(2,col) = max(color_range(2,col),i1)

        end if

      else if (color(i1) /= color(i2)) then

        if (color(i1) < color(i2)) then
          src_col = color(i1)
          dst_col = color(i2)
        else
          src_col = color(i2)
          dst_col = color(i1)
        end if

        color_range(1,src_col) = min(color_range(1,src_col), &
                                     color_range(1,dst_col))
        color_range(2,src_col) = max(color_range(2,src_col), &
                                     color_range(2,dst_col))

        do j = color_range(1,dst_col), color_range(2,dst_col)
          if (color(j) == dst_col) &
            color(j) = src_col
        end do

        color_range(1,dst_col) = 0
        color_range(2,dst_col) = 0

      end if

    end do

    ! define molecule number
    nmol = 0
    pcol = -1

    do i = 1, hm_num_atoms

      col = color(i)

      if (pcol /= col .or. col == 0) &
           nmol = nmol + 1

      call hm_set_molecule_no(i, nmol)

      pcol = col
    end do

    deallocate(color, color_range)

    if (main_rank) then
      write(MsgOut,'(a,i9)') &
           'Setup_Molecule_Variables> # of molecules : ', nmol
      write(MsgOut,'(a)') ''
    end if

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    return

  end subroutine setup_molecule_variables

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_domain_rst_files
  !> @brief        setup domain restart files
  !! @authors      NT
  !! @param[in]    ctrl_data : information of control parameters
  !! @param[inout] boundary  : boundary information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_domain_rst_files(ctrl_data, boundary)

    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_boundary),        intent(inout) :: boundary

    ! local variables
    type(s_domain)           :: domain
    integer                  :: idom, ndom, ierr
    integer                  :: i, ix, ig
    integer                  :: nplace
    character(200)           :: filename
    logical                  :: bex


    call init_boundary(boundary)

    ! check restart files are input
    !
    if (ctrl_data%inp_info%rstfile == '') then

      if (main_rank) then
        write(MsgOut,'(a)') '   ... skip'
        write(MsgOut,'(a)') ' '
      end if

      pio_restart = .false.

      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      return

    end if

    ! compute # of domains 
    !
    my_world_rank = 0
    my_city_rank  = 0

    do nplace = 1, 7
      filename = pio_get_ranked_filename(ctrl_data%inp_info%rstfile, &
                                         my_world_rank, nplace)
      inquire(file=filename, exist=bex)
      if (bex) &
        exit
    end do

    if (main_rank) then
      write(MsgOut,'(a,i10)')  'Restart file : # of places : ', nplace+1
      write(MsgOut,'(a)')      ' '
    end if

    call pio_read_domain_str(filename, boundary, domain,  &
                             ctrl_data%opt_info%convert)

    if (ctrl_data%opt_info%set_restart) pio_restart = .true.

    if (main_rank) then
      write(MsgOut,'(a)')      'Input domain information : '
      write(MsgOut,'(a,3i10)') '  domains (x,y,z) = ', boundary%num_domain(1), &
                                                       boundary%num_domain(2), &
                                                       boundary%num_domain(3)
    end if

    ndom = boundary%num_domain(1) * &
           boundary%num_domain(2) * &
           boundary%num_domain(3)


    ! update coordinates and velocities
    !

    do idom = 1, ndom

      if (main_rank) &
        write(MsgOut,'(a,i6,a,i6)') '  read ', idom, '/', ndom

      my_world_rank = idom - 1
      my_city_rank  = idom - 1
      filename      = pio_get_ranked_filename(ctrl_data%inp_info%rstfile, &
                                              my_world_rank, nplace)

      call pio_read_domain_str(filename, boundary, domain,  &
                             ctrl_data%opt_info%convert)

      do i = 1, domain%num_cell_local
        do ix = 1, domain%num_atom(i)
          ig = domain%id_l2g(ix,i)
          call hm_set_atom_coord(1, ig, domain%coord(1,ix, i))
          call hm_set_atom_coord(2, ig, domain%coord(2,ix, i))
          call hm_set_atom_coord(3, ig, domain%coord(3,ix, i))
          call hm_set_atom_velocity(1, ig, domain%velocity(1,ix, i))
          call hm_set_atom_velocity(2, ig, domain%velocity(2,ix, i))
          call hm_set_atom_velocity(3, ig, domain%velocity(3,ix, i))
        end do
      end do

    end do


    ! deallocate
    !

    call dealloc_domain_all(domain)

    my_world_rank = 0
    my_city_rank  = 0

    if (main_rank) &
      write(MsgOut,'(a)') ' '

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    return

  end subroutine setup_domain_rst_files

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_domain_variables
  !> @brief        setup domain variables
  !! @authors      NT
  !! @param[in]    ctrl_data   : information of control parameters
  !! @param[in]    par         : CHARMM PAR information
  !! @pmara[in]    prmtop      : AMBER parameter/topology information
  !! @param[in]    grotop      : GROMACS topology information
  !! @param[in]    boundary_rst: boundary condition information on restart
  !! @param[inout] boundary    : boundary condition information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamic variable information
  !! @param[inout] dynamics    : dynamics information
  !! @param[inout] constraints : domain constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_domain_variables(ctrl_data, par, prmtop, grotop, &
                                    boundary_rst, boundary, enefunc, dynvars, &
                                    dynamics, restraints, constraints)

    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_par),             intent(in)    :: par
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_grotop),          intent(in)    :: grotop
    type(s_boundary),        intent(in)    :: boundary_rst
    type(s_boundary),        intent(inout) :: boundary
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_restraints),      intent(inout) :: restraints
    type(s_constraints),     intent(inout) :: constraints

    integer                  :: ierr


    ! setups only to be related to the whole system
    !
    
    call setup_md_peer_0(ctrl_data,    &
                         par,          &
                         prmtop,       &
                         grotop,       &
                         boundary_rst, &
                         boundary,     &
                         enefunc,      &
                         dynvars,      &
                         dynamics,     &
                         constraints,  &
                         restraints)

    if (main_rank) &
      write(MsgOut,'(a)') ' '

    call dealloc_restraints(restraints, RestraintsGroup)
    call dealloc_restraints(restraints, RestraintsFunc)

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    return

  end subroutine setup_domain_variables

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_domain_index_list
  !> @brief        setup domain index list
  !! @authors      NT
  !! @param[in]    ctrl_data : information of control parameters
  !! @param[in]    boundary  : boundary condition information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_domain_index_list(ctrl_data, boundary)

    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_boundary),        intent(in)    :: boundary
    
    ! local variables
    integer                  :: ierr


    if (main_rank) &
      write(MsgOut,'(a)') 'Setup_Domain_Index_List> '

    di_file_dir = ctrl_data%opt_info%cachepath


    call di_create(boundary)

    if (main_rank) &
      write(MsgOut,'(a)') ''

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    return

  end subroutine setup_domain_index_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_domain_restart_files
  !> @brief        output domain restart files
  !! @authors      NT, JJ
  !! @param[in]    ctrl_data    : information of control parameters
  !! @param[in]    par          : CHARMM PAR information
  !! @param[in]    prmtop       : AMBER parameter/topology information
  !! @param[in]    grotop       : GROMACS topology information
  !! @param[in]    bundary      : boundary condition information
  !! @param[in]    dynvars      : dynamic variables information
  !! @param[in]    dynamics     : dynamics information
  !! @param[inout] enefuncO     : potential energy functions information
  !! @param[inout] constraintsO : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_domain_restart_files(block_domain, ctrl_data, par, prmtop, &
                                         grotop, boundary, dynvars, dynamics,  &
                                         restraints, enefuncO, constraintsO)

    ! formal arguments
    logical,                 intent(in)    :: block_domain
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_par),             intent(in)    :: par
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_grotop),          intent(in)    :: grotop
    type(s_boundary),        intent(in)    :: boundary
    type(s_dynvars),         intent(in)    :: dynvars
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_restraints),      intent(in)    :: restraints
    type(s_enefunc),         intent(inout) :: enefuncO
    type(s_constraints),     intent(inout) :: constraintsO
    
    ! local variables
    type(s_enefunc)          :: enefunc
    type(s_constraints)      :: constraints
    type(s_domain_index)     :: domain_index
    type(s_domain)           :: domain
    type(s_pio_t0_info)      :: t0_info

    integer                  :: local_my_rank, j, ierr
    integer                  :: idom, ndom
    character(200)           :: filename
    character(200)           :: filename1


    ! disable varbose standard output
    !
    main_rank = .false.


    ! start domain output
    !

    ndom = boundary%num_domain(1) * &
           boundary%num_domain(2) * &
           boundary%num_domain(3)


    do idom = 1+my_prst_rank, ndom, nproc_prst

      write(MsgOut,*) 'rank:', my_prst_rank, '  idom:', idom

      ! copy setup-once constraints, enefunc
      !
      enefunc     = enefuncO
      constraints = constraintsO

      ! get domain index information
      !
      call di_get_index(idom, domain_index)

      ! setup domain
      !
      local_my_rank  = idom - 1

      call setup_md_peer(local_my_rank, &
                         ctrl_data,     &
                         par,           &
                         prmtop,        &
                         grotop,        &
                         domain_index,  &
                         boundary,      &
                         enefunc,       &
                         constraints,   &
                         domain)

      ! create ranked filename
      !
      my_world_rank = local_my_rank
      my_city_rank  = local_my_rank
      filename      = pio_get_ranked_filename(ctrl_data%out_info%rstfile, &
                                              my_world_rank)
      filename1     = ctrl_data%out_info%selfile  

      ! output domain restart file
      !
  
      call setup_the(t0_info, ctrl_data, domain)

      if (local_my_rank == 0 .and. filename1 /= '')  &
        call pio_write_selection(filename1, restraints)
 
      call pio_write_domain_rst(filename,           &
                                boundary,           &
                                domain,             &
                                enefunc,            &
                                constraints,        &
                                dynvars,            &
                                dynamics)
    end do

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    if (main_rank) then
      write(MsgOut,'(a)') '  done.'
      write(MsgOut,'(a)') ''
    end if

    main_rank     = (my_prst_rank == 0)
    nproc_world   = 1
    nproc_city    = 1
    nproc_country = 1

    return

  end subroutine output_domain_restart_files

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output pdb file
  !> @brief        output molecule coordinates as PDB file
  !! @authors      NT
  !! @param[in]    ctrl_data : information of control parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_pdb_file(ctrl_data)
    
    ! formal arguments
    type(s_ctrl_data),       intent(in) :: ctrl_data

    integer                  :: old_mno
    integer                  :: i, pdb


    if (ctrl_data%out_info%pdbfile == '') then

      if (main_rank) then
        write(MsgOut,'(a)') '   ... skip'
        write(MsgOut,'(a)') ' '
      end if

      return

    end if

    if (main_rank) then

      old_mno = hm_molecule_no(1)

      call open_file(pdb, ctrl_data%out_info%pdbfile, IOFileOutputReplace)

      do i = 1, hm_num_atoms

        if (old_mno /= hm_molecule_no(i)) then
          write(pdb, '("TER ",i7,6x,a4,i6,48x,a1)') &
               mod(hm_atom_no(i)-1,1000000),  &
               hm_residue_name(i), &
               hm_residue_no(i),   &
               ' '
          old_mno = hm_molecule_no(i)
        end if

        write(pdb, '("ATOM",i7,1x,a4,1x,a4,i6,3x,3f8.3,f6.2,f6.2,6x,a4)') &
             mod(hm_atom_no(i)-1,1000000), &
             hm_atom_name(i),    &
             hm_residue_name(i), &
             hm_residue_no(i),   &
             hm_atom_coord(1,i), &
             hm_atom_coord(2,i), &
             hm_atom_coord(3,i), &
             1.0_wp,             &
             1.0_wp,             &
             hm_segment_name(i)
      end do

      call close_file(pdb)

    end if

    if (main_rank) then
      write(MsgOut,'(a)') '  done.'
      write(MsgOut,'(a)') ''
    end if

    return

  end subroutine output_pdb_file

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_the
  !> @brief        setup the parallel I/O local information
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_the(t0_info, ctrl_data, domain)

    ! formal arguments
    type(s_pio_t0_info),     intent(inout) :: t0_info
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_domain),  target, intent(in)    :: domain

    ! local variables
    integer                  :: nvar


    t0_info%input_topfile         = ctrl_data%inp_info%topfile
    t0_info%input_parfile         = ctrl_data%inp_info%parfile
    t0_info%energy_pairlistdist   = ctrl_data%ene_info%pairlistdist
    t0_info%energy_table          = ctrl_data%ene_info%table
    t0_info%energy_watermodel     = ctrl_data%ene_info%water_model
    t0_info%constraint_rigidbond  = ctrl_data%cons_info%rigid_bond
    t0_info%constraint_fastwater  = ctrl_data%cons_info%fast_water
    t0_info%constraint_watermodel = ctrl_data%cons_info%water_model
    t0_info%ensemble_type         = ctrl_data%ens_info%ensemble
    t0_info%boundary_boxsizex     = ctrl_data%bou_info%pbc_info%box_size_x
    t0_info%boundary_boxsizey     = ctrl_data%bou_info%pbc_info%box_size_y
    t0_info%boundary_boxsizez     = ctrl_data%bou_info%pbc_info%box_size_z
    t0_info%boundary_originx      = ctrl_data%bou_info%origin_x
    t0_info%boundary_originy      = ctrl_data%bou_info%origin_y
    t0_info%boundary_originz      = ctrl_data%bou_info%origin_z

    if (ctrl_data%opt_info%domain_xyz /= 0) then
      t0_info%boundary_domain_xyz = ctrl_data%opt_info%domain_xyz
    else
      t0_info%boundary_domain_xyz = &
                      ctrl_data%bou_info%domain_x * &
                      ctrl_data%bou_info%domain_y * &
                      ctrl_data%bou_info%domain_z
    end if

    if (allocated(ctrl_data%sel_info%groups)) then
      nvar = size(ctrl_data%sel_info%groups)
    else
      nvar = 0
    end if
    if (allocated(t0_info%selection_group)) then
      if (size(t0_info%selection_group) /= nvar) then
        deallocate(t0_info%selection_group)
        allocate  (t0_info%selection_group(nvar))
      end if
    else
      allocate(t0_info%selection_group(nvar))
    end if
    t0_info%selection_group(1:nvar) = ctrl_data%sel_info%groups(1:nvar)

    if (allocated(ctrl_data%res_info%function)) then
      nvar = size(ctrl_data%res_info%function)
    else
      nvar = 0
    end if
    if (allocated(t0_info%restraint_func)) then
      if (size(t0_info%restraint_func) /= nvar) then
        deallocate(t0_info%restraint_func,  &
                   t0_info%restraint_const, &
                   t0_info%restraint_index)
        allocate  (t0_info%restraint_func(nvar),  &
                   t0_info%restraint_const(nvar), &
                   t0_info%restraint_index(nvar))
      end if
    else
      allocate(t0_info%restraint_func(nvar),  &
               t0_info%restraint_const(nvar), &
               t0_info%restraint_index(nvar))
    end if
    t0_info%restraint_func(1:nvar)    = ctrl_data%res_info%function(1:nvar)
    t0_info%restraint_const(1:nvar)   = ctrl_data%res_info%constant(1:nvar)
    t0_info%restraint_index(1:nvar)   = ctrl_data%res_info%select_index(1:nvar)

    return

  end subroutine setup_the

end module pr_run_prst_setup_mod

