!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   rg_generate_mod
!> @brief   generate RPATH input files
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module rg_generate_mod

  use rg_option_mod
  use rg_option_str_mod
  use pbc_correct_mod
  use fileio_trj_mod
  use fitting_mod
  use fitting_str_mod
  use output_str_mod
  use input_str_mod
  use select_atoms_str_mod
  use molecules_str_mod
  use trajectory_str_mod
  use fileio_grocrd_mod
  use fileio_prmtop_mod
  use fileio_pdb_mod
  use fileio_rst_mod
  use fileio_mod
  use string_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type :: s_cv
    real(wp), allocatable :: ts(:)
    real(wp), allocatable :: cv(:,:)
    real(wp), allocatable :: cv2(:,:)
  end type s_cv

  type :: s_idxrep
    integer               :: idx
    integer               :: rep
  end type s_idxrep

  ! subroutine
  public  :: generate
  private :: open_cv
  private :: open_dcd
  private :: compute_path
  private :: reparametrize_path
  private :: output_snapshot
  private :: output_snapshot_pdb
  private :: output_snapshot_rst
  private :: output_control
  private :: define_pdb
  private :: get_distance
  private :: get_numbered_filename
  private :: sort_idxrep
  private :: cumsum
  private :: check_cvfile
  private :: check_dcdfile

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    generate
  !> @brief        generate molecule files
  !! @authors      NT
  !! @param[in]    molecule : molecule information
  !! @param[in]    input    : input information
  !! @param[in]    output   : output information
  !! @param[in]    option   : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine generate(molecule, input, output, fitting, option)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_input),           intent(in)    :: input
    type(s_output),          intent(in)    :: output
    type(s_fitting),         intent(inout) :: fitting
    type(s_option),          intent(in)    :: option

    ! local variables
    type(s_cv)               :: cv
    type(s_pdb)              :: pdb
    real(wp),   allocatable  :: path(:,:), path_equi(:,:)
    integer                  :: i


    if (option%is_cartesian) then 

      call open_dcd(input%dcdfile, molecule, fitting, option, cv)

    else

      call open_cv(input%cvfile, cv)

    end if


    call compute_path(option%nreplica, cv, path)


    do i = 1, option%iter_reparam
      call reparametrize_path(path, path_equi)
      path(:, :) = path_equi(:, :)
    end do


    call define_pdb(input%pdbfile,    &
                    input%prmtopfile, &
                    input%grocrdfile, &
                    pdb)


    call output_snapshot(pdb, cv, path_equi, &
                         input%dcdfile,    &
                         input%dcdvelfile, &
                         output%pdbfile,   &
                         output%rstfile,   &
                         molecule, fitting, option)


    call output_control(path_equi, input, output, fitting, option)


    call dealloc_pdb_all(pdb)
    deallocate(path, path_equi)

    return

  end subroutine generate

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine open_cv(cvfile, cv)

    ! formal arguments
    character(*),            intent(in)    :: cvfile
    type(s_cv),              intent(inout) :: cv

    ! local variables
    integer                  :: file, ndim, nstep, i, j
    character(1000)          :: line


    call open_file(file, cvfile, IOFileInput)

    nstep = 1

    read(file,'(a)') line
    ndim = split_num(line)
    ndim = ndim - 1

    do while(.true.)
      read(file,'(a)',end=10,err=10) line
      nstep = nstep + 1
    end do

10  write(MsgOut,'(a,i8)') 'Open_Cv> # of dimension : ', ndim
    write(MsgOut,'(a,i8)') '         # of steps     : ', nstep
    write(MsgOut,'(a)')    ' '

    allocate(cv%ts(nstep), cv%cv(ndim,nstep), cv%cv2(ndim,nstep))

    rewind(file)

    do i = 1, nstep
      read(file,*) cv%ts(i),(cv%cv(j,i),j=1,ndim)
      cv%cv2(:,i) = cv%cv(:,i)
    end do

    call close_file(file)

    return

  end subroutine open_cv

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine open_dcd(dcdfile, molecule, fitting, option, cv)

    ! formal arguments
    character(*),            intent(in)    :: dcdfile
    type(s_molecule),        intent(in)    :: molecule
    type(s_fitting),         intent(inout) :: fitting
    type(s_option),          intent(in)    :: option
    type(s_cv),              intent(inout) :: cv

    ! local variables
    type(s_trajectory)       :: trajectory
    type(s_trj_file)         :: file
    integer                  :: i, j, k, istep
    integer                  :: nbrella, ibrella, ndim, nstep, natom, nfunc
    character(MaxFilename)   :: filename
    integer                  :: fitting_method_old


    write(MsgOut,'(a)') 'Open_Dcd> '

    ndim = size(option%selatoms(option%cv_atom)%idx) * 3

    ! check nstep and allocate trajectory
    !
    call check_dcdfile(dcdfile, nstep, natom)


    ! allocate data_k
    !
    allocate(cv%ts(nstep), cv%cv(ndim,nstep), cv%cv2(ndim,nstep))


    if (natom /= option%num_atoms) &
      call error_msg( &
      'Open_Dcd> Dcd atom count is different from PDB/PSF/PRMTOP.')

    call alloc_trajectory(trajectory, natom)

    call open_trj(file, dcdfile, TrjFormatDCD, TrjTypeCoorBox, IOFileInput)

    do istep = 1, nstep
      cv%ts(istep) = istep

      call read_trj(file, trajectory)

      if (fitting%fitting_method /= FittingMethodNO) then
        if (fitting%mass_weight) then
          call run_fitting(fitting, &
                           molecule%atom_coord, &
                           trajectory%coord, &
                           trajectory%coord, &
                           molecule%mass)
        else
          call run_fitting(fitting, &
                           molecule%atom_coord, &
                           trajectory%coord, &
                           trajectory%coord)
        end if
      end if

      do i = 1, size(option%selatoms(option%cv_atom)%idx)
        cv%cv((i-1)*3+1, istep) = trajectory%coord(1, option%selatoms(option%cv_atom)%idx(i))
        cv%cv((i-1)*3+2, istep) = trajectory%coord(2, option%selatoms(option%cv_atom)%idx(i))
        cv%cv((i-1)*3+3, istep) = trajectory%coord(3, option%selatoms(option%cv_atom)%idx(i))
      end do

      if (fitting%fitting_method /= FittingMethodNO) then
        fitting_method_old = fitting%fitting_method

        if ((fitting_method_old == FittingMethodTR_ROT) .or. &
           (fitting_method_old == FittingMethodTR_ZROT)) then

          fitting%fitting_method = FittingMethodTR

        else if (fitting_method_old == FittingMethodXYTR_ZROT) then

          fitting%fitting_method = FittingMethodXYTR

        end if

        if (fitting%mass_weight) then
          call run_fitting(fitting, &
                           molecule%atom_coord, &
                           trajectory%coord, &
                           trajectory%coord, &
                           molecule%mass)
        else
          call run_fitting(fitting, &
                           molecule%atom_coord, &
                           trajectory%coord, &
                           trajectory%coord)
        end if

        fitting%fitting_method = fitting_method_old
      end if

      do i = 1, size(option%selatoms(option%cv_atom)%idx)
        cv%cv2((i-1)*3+1, istep) = trajectory%coord(1, option%selatoms(option%cv_atom)%idx(i))
        cv%cv2((i-1)*3+2, istep) = trajectory%coord(2, option%selatoms(option%cv_atom)%idx(i))
        cv%cv2((i-1)*3+3, istep) = trajectory%coord(3, option%selatoms(option%cv_atom)%idx(i))
      end do

    end do

    call close_trj(file)

    write(MsgOut,'(a)') '  ..done'
    write(MsgOut,'(a)') ''

    return

  end subroutine open_dcd

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_path(nreplica, cv, path)

    ! formal arguments
    integer,                 intent(in)    :: nreplica
    type(s_cv),              intent(inout) :: cv
    real(wp),  allocatable,  intent(inout) :: path(:,:)

    ! local variables
    integer                  :: ndim, nstep, ireplica, i, count
    real(wp),  allocatable   :: accum(:)


    ndim  = size(cv%cv(:,1))
    nstep = size(cv%cv(1,:))

    allocate(path(ndim, nreplica))
    allocate(accum(ndim))

    path(:,1)        = cv%cv(:,1)
    path(:,nreplica) = cv%cv(:,nstep)

    cv%ts(:) = (cv%ts(:) / maxval(cv%ts(:))) * (nreplica - 1)

    do ireplica = 1, nreplica-2
      accum(:) = 0
      count = 0
      do i = 1, nstep
        if (ireplica-0.5_wp <= cv%ts(i) .and. &
            cv%ts(i) < (ireplica + 0.5_wp)) then
          accum(:) = accum(:) + cv%cv(:,i)
          count = count + 1
        end if
      end do
      path(:,ireplica+1) = accum / real(count,wp)
    end do

    deallocate(accum)

    return

  end subroutine compute_path

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine reparametrize_path(path, path_equi)

    ! formal arguments
    real(wp),                intent(in)    :: path(:,:)
    real(wp), allocatable,   intent(inout) :: path_equi(:,:)

    ! local variables
    real(wp)                 :: d
    integer                  :: ndim, nreplica, i, q

    real(wp), allocatable    :: dist(:), dist_cumsum(:), dist_equi(:)


    ndim     = size(path(:,1))
    nreplica = size(path(1,:))

    if (.not. allocated(path_equi)) then
      allocate(path_equi(ndim, nreplica))
    end if
    allocate(dist(nreplica), dist_cumsum(nreplica), dist_equi(nreplica))

    dist(1) = 0.0_wp
    do i = 2, nreplica
      d = get_distance(path(:,i), path(:,i-1))
      dist(i) = d
    end do

    call cumsum(dist, dist_cumsum)

    do i = 1, nreplica
      d = dist_cumsum(nreplica)*real(i-1,wp)/real(nreplica-1)
      dist_equi(i) = d
    end do

    path_equi(:,1) = path(:,1)
    do i = 2, nreplica - 1
      q = 2
      do while (dist_cumsum(q) < dist_equi(i))
        q = q + 1
      end do
      path_equi(:,i) = path(:,q-1) + (path(:,q)-path(:,q-1)) * &
                                     (dist_equi(i) - dist_cumsum(q-1))/ &
                                     (dist_cumsum(q) - dist_cumsum(q-1))
      
    end do
    path_equi(:,nreplica) = path(:,nreplica)

    deallocate(dist, dist_cumsum, dist_equi)

    return

  end subroutine reparametrize_path

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_snapshot(pdb, cv, path, dcdfile, dcdvelfile, pdbfile, rstfile, &
                             molecule, fitting, option)

    ! formal arguments
    type(s_pdb),             intent(inout) :: pdb
    type(s_cv),              intent(in)    :: cv
    real(wp),                intent(in)    :: path(:,:)
    character(*),            intent(in)    :: dcdfile
    character(*),            intent(in)    :: dcdvelfile
    character(*),            intent(in)    :: pdbfile
    character(*),            intent(in)    :: rstfile
    type(s_molecule),        intent(in)    :: molecule
    type(s_fitting),         intent(inout) :: fitting
    type(s_option),          intent(in)    :: option
    
    ! local variables
    real(wp)                 :: d, dmin
    integer                  :: nstep, nreplica
    integer                  :: i, j, istep, ireplica, idx, prv_idx
    type(s_trj_file)         :: dcd, dvl

    type(s_trajectory)       :: trjcrd, trjvel
    type(s_idxrep), allocatable :: idxrep(:)

    integer                  :: fitting_method_old


    ! open trajectory
    !
    dcd%unit_no = 0
    dvl%unit_no = 0

    if (dcdfile /= '') then
      call open_trj(dcd, dcdfile, TrjFormatDCD, TrjTypeCoorBox, IOFileInput)
      call alloc_trajectory(trjcrd, molecule%num_atoms)
    else
      call error_msg('Output_Snapshot> Dcd file must be given.')
    end if

    if (dcdvelfile /= '') then
      call open_trj(dvl, dcdvelfile, TrjFormatDCD, TrjTypeCoorBox, IOFileInput)
      call alloc_trajectory(trjvel, molecule%num_atoms)
    end if


    ! select trajectory step for replica
    !
    nreplica = size(path(1,:))
    nstep    = size(cv%cv(1,:))

    allocate(idxrep(nreplica))

    do ireplica = 1, nreplica

      dmin = 100000.0_wp
      idx  = -1
      do istep = 1, nstep
        d = get_distance(path(:,ireplica),cv%cv2(:,istep))
        if (d < dmin) then
          dmin = d
          idx  = istep
        end if
      end do

      write(MsgOut,'(a,i5,a,i5,a)') &
           'Output_Snapshot> ', idx,      ' th snapshot is used for ', &
                                ireplica, ' th replica.'
      write(MsgOut,'(a,f10.5)') &
           '                   Euclid distance between the image and snapshot is ', dmin
      write(MsgOut,'(a)') ''

      idxrep(ireplica)%idx = idx
      idxrep(ireplica)%rep = ireplica

    end do
    write(MsgOut,'(a)') ''

    !call sort_idxrep(idxrep, 1, nreplica)


    ! read trajectory and write restart/pdb file
    !
    write(MsgOut,'(a)') 'Output_SnapShot> write restart/pdb file'

    prv_idx = 1

    do istep = 1, nstep

      if (dcd%unit_no /= 0) then
        call read_trj(dcd, trjcrd)
      end if

      if (fitting%fitting_method /= FittingMethodNO) then
        fitting_method_old = fitting%fitting_method

        if ((fitting_method_old == FittingMethodTR_ROT) .or. &
            (fitting_method_old == FittingMethodTR_ZROT)) then

          fitting%fitting_method = FittingMethodTR

        else if (fitting_method_old == FittingMethodXYTR_ZROT) then

          fitting%fitting_method = FittingMethodXYTR

        end if

        call run_fitting(fitting, &
                         molecule%atom_coord, &
                         trjcrd%coord,    &
                         trjcrd%coord)

        fitting%fitting_method = fitting_method_old
      end if

      if (dvl%unit_no /= 0) then
        call read_trj(dvl, trjvel)
      end if

      do ireplica = 1, nreplica
        if (istep .eq. idxrep(ireplica)%idx) then
          write(MsgOut,'(a,i6,a,i6)') &
            '   snapshot : ', istep, '  replica : ', ireplica
          if (rstfile /= '') &
            call output_snapshot_rst &
              (molecule%num_atoms, trjcrd, path, trjvel, ireplica, rstfile, option)
          if (trjcrd%pbc_box(3,3) > 0.0_wp) then
            call run_pbc_correct(PBCCModeMolecule,       &
                                 molecule,               &
                                 trjcrd)
          end if
          if (pdbfile /= '') &
            call output_snapshot_pdb &
              (pdb, trjcrd, path, ireplica, pdbfile, option)
        end if
      end do
    end do

    deallocate(idxrep)

    write(MsgOut,'(a)') ''


    ! close trjectory
    !
    call dealloc_trajectory(trjcrd)
    call dealloc_trajectory(trjvel)

    if (dcd%unit_no /= 0) &
      call close_trj(dcd)

    if (dvl%unit_no /= 0) &
      call close_trj(dvl)

    return

  end subroutine output_snapshot

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_snapshot_pdb(pdb, trjcrd, path, &
                                 irep, pdbfile, option)

    ! formal arguments
    type(s_pdb),             intent(inout) :: pdb
    type(s_trajectory),      intent(in)    :: trjcrd
    real(wp),                intent(in)    :: path(:,:)
    integer,                 intent(in)    :: irep
    character(*),            intent(in)    :: pdbfile
    type(s_option),          intent(in)    :: option

    ! local variables
    integer                  :: i
    

    pdb%atom_coord(1:3,1:pdb%num_atoms) = trjcrd%coord(1:3,1:pdb%num_atoms)

    if (option%is_cartesian) then
      do i = 1, size(option%selatoms(option%cv_atom)%idx)
        pdb%atom_coord(1, option%selatoms(option%cv_atom)%idx(i)) = path((i-1)*3+1, irep)
        pdb%atom_coord(2, option%selatoms(option%cv_atom)%idx(i)) = path((i-1)*3+2, irep)
        pdb%atom_coord(3, option%selatoms(option%cv_atom)%idx(i)) = path((i-1)*3+3, irep)
      end do
    end if

    call output_pdb(get_numbered_filename(pdbfile, irep), pdb)

    return

  end subroutine output_snapshot_pdb

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_snapshot_rst(num_atoms, trjcrd, path, trjvel, &
                                 irep, rstfile, option)

    ! formal arguments
    integer,                 intent(in)    :: num_atoms
    type(s_trajectory),      intent(in)    :: trjcrd
    real(wp),                intent(in)    :: path(:,:)
    type(s_trajectory),      intent(in)    :: trjvel
    integer,                 intent(in)    :: irep
    character(*),            intent(in)    :: rstfile
    type(s_option),          intent(in)    :: option
    
    ! local variables
    type(s_rst)              :: rst
    integer                  :: ireplica, idim


    rst%iseed_rpath = option%iseed + 362893
    rst%nreplicas   = option%nreplica
    rst%dimension   = size(path(:,1))

    call alloc_rst(rst, RestartAtom, num_atoms)
    call alloc_rst(rst, RestartRpath, rst%dimension, rst%nreplicas)

    do ireplica = 1, rst%nreplicas
      do idim = 1, rst%dimension
        rst%rest_reference(1:2,idim,ireplica) = path(idim, ireplica)
      end do
    end do

    if (allocated(trjvel%coord)) then

      ! md restart file

      rst%num_atoms                     = num_atoms
      rst%rstfile_type                  = RstfileTypeMd
      rst%iseed                         = option%iseed + irep
      rst%box_size_x                    = trjcrd%pbc_box(1,1)
      rst%box_size_y                    = trjcrd%pbc_box(2,2)
      rst%box_size_z                    = trjcrd%pbc_box(3,3)
      rst%thermostat_momentum           = 0.0_wp
      rst%barostat_momentum(1:3)        = 0.0_wp
      rst%coord   (1:3,1:num_atoms) = trjcrd%coord(1:3,1:num_atoms)
      rst%velocity(1:3,1:num_atoms) = trjvel%coord(1:3,1:num_atoms)

      call output_rst(get_numbered_filename(rstfile, irep), rst)

    else

      ! min restart file

      rst%num_atoms                     = num_atoms
      rst%rstfile_type                  = RstfileTypeMin
      rst%energy                        = 0.0_wp
      rst%delta_r                       = 0.0_wp
      rst%box_size_x                    = trjcrd%pbc_box(1,1)
      rst%box_size_y                    = trjcrd%pbc_box(2,2)
      rst%box_size_z                    = trjcrd%pbc_box(3,3)
      rst%coord(1:3,1:num_atoms)        = trjcrd%coord(1:3,1:num_atoms)

      call output_rst(get_numbered_filename(rstfile, irep), rst)

    end if

    call dealloc_rst_all(rst)


    return

  end subroutine output_snapshot_rst

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_control(path, input, output, fitting, option)

    ! formal arguments
    real(wp),                intent(in)    :: path(:,:)
    type(s_input),           intent(in)    :: input
    type(s_output),          intent(in)    :: output
    type(s_fitting),         intent(in)    :: fitting
    type(s_option),          intent(in)    :: option

    ! local variables
    integer                  :: nreplica, ndim, i, idim
    character(100)           :: fmt
    character(100)           :: fmt2


    ndim     = size(path(:,1))
    nreplica = size(path(1,:))

    write(fmt,'(a,i0,a)') '(a,i0,a,',nreplica,'(f0.2,1x))'
    write(fmt2,'(a,i0,a)') '(a,',nreplica,'(i0,1x))'

    if (option%is_cartesian) then
      write(MsgOut,'(a)') ''
      write(MsgOut,'(a)') 'copy the following for RPATH simulation'
      write(MsgOut,'(a)') '--------------------------------------'
      write(MsgOut,'(a)')           '[INPUT]'
      write(MsgOut,'(a,a)')         'fitfile       = ', trim(input%fitfile)
      write(MsgOut,'(a,a)')         'reffile       = ', trim(output%pdbfile)
      write(MsgOut,'(a,a)')         'rstfile       = ', trim(output%rstfile)
      write(MsgOut,'(a)')           ''
      write(MsgOut,'(a)')           '# for equilibration'
      write(MsgOut,'(a)')           '[RPATH]'
      write(MsgOut,'(a,i0)')        'nreplica      = ', nreplica
      write(MsgOut,'(a)')           'rpath_period  = 0'
      write(MsgOut,'(a,i0)')        'rest_function = ', 1
      write(MsgOut,'(a)')           ''
      write(MsgOut,'(a)')           '# for equilibration'
      write(MsgOut,'(a)')           '[RPATH]'
      write(MsgOut,'(a,i0)')        'nreplica       = ', nreplica
      write(MsgOut,'(a)')           'rpath_period   = 1000'
      write(MsgOut,'(a)')           'delta          = 0.001'
      write(MsgOut,'(a)')           'smooth         = 0.0'
      write(MsgOut,'(a,i0)')        'rest_function  = ', 1
      write(MsgOut,'(a)')           'fix_terminal   = NO'
      if (fitting%fitting_method == FittingMethodTR_ROT) then
        write(MsgOut,'(a)')         'fitting_method = TR+ROT'
        write(MsgOut,'(a,i0)')        'fitting_atom   = ', option%fitting_atom
      else if (fitting%fitting_method == FittingMethodTR_ZROT) then
        write(MsgOut,'(a)')         'fitting_method = TR+ZROT'
        write(MsgOut,'(a,i0)')      'fitting_atom   = ', option%fitting_atom
      else if (fitting%fitting_method == FittingMethodXYTR_ZROT) then
        write(MsgOut,'(a)')         'fitting_method = XYTR+ZROT'
        write(MsgOut,'(a,i0)')      'fitting_atom   = ', option%fitting_atom
      else 
        write(MsgOut,'(a)')         'fitting_method = NO'
      end if
      write(MsgOut,'(a)')           ''
      write(MsgOut,'(a)')           '[SELECTION]'
      do i = 1, size(option%groups)
        write(MsgOut,'(a,i0,a,a)') 'group', i, ' = ', trim(option%groups(i))
      end do
      write(MsgOut,'(a)')           ''
      write(MsgOut,'(a)')           '[RESTRAINTS]'
      write(MsgOut,'(a,i0)')        'nfunctions     = ', 1
      write(MsgOut,'(a,i0,a)')      'function', 1, '     = POSI'
      write(MsgOut,fmt)             'constant', 1, '     = ', (1.0_wp,i=1,nreplica)
      write(MsgOut,'(a,i0,a,i0)')   'select_index', 1, ' = ', option%cv_atom
      write(MsgOut,'(a)')           ''

    else
      write(MsgOut,'(a)') ''
      write(MsgOut,'(a)') 'copy the follwing for RPATH simulation'
      write(MsgOut,'(a)') '--------------------------------------'
      write(MsgOut,'(a)')           '[INPUT]'
      write(MsgOut,'(a,a)')         'rstfile       = ', trim(output%rstfile)
      write(MsgOut,'(a)')           ''
      write(MsgOut,'(a)')           '# for equilibration'
      write(MsgOut,'(a)')           '[RPATH]'
      write(MsgOut,'(a,i0)')        'nreplica      = ', nreplica
      write(MsgOut,'(a)')           'rpath_period  = 0'
      write(MsgOut,'(a,1000(i0,1x))') 'rest_function = ', (i,i=1,ndim)
      write(MsgOut,'(a)')           ''
      write(MsgOut,'(a)')           '# for string method'
      write(MsgOut,'(a)')           '[RPATH]'
      write(MsgOut,'(a,i0)')        'nreplica      = ', nreplica
      write(MsgOut,'(a)')           'rpath_period  = 1000'
      write(MsgOut,'(a)')           'delta         = 0.001'
      write(MsgOut,'(a)')           'smooth        = 0.0'
      write(MsgOut,fmt2)            'rest_function = ', (i,i=1,ndim)
      write(MsgOut,'(a)')           'fix_terminal  = NO'
      write(MsgOut,'(a)')           ''
      write(MsgOut,'(a)')           '[RESTRAINTS]'
      write(MsgOut,'(a,i0)')        'nfunctions     = ', ndim
      do idim = 1, ndim
        write(MsgOut,'(a,i0,a)')    'function', idim, '     = XXX'
        write(MsgOut,fmt)           'constant', idim, '     = ', (1.0_wp,i=1,nreplica)
        write(MsgOut,fmt)           'reference', idim, '    = ', (path(idim,i),i=1,nreplica)
        write(MsgOut,'(a,i0,a)')    'select_index', idim, ' = YYY'
        write(MsgOut,'(a)')         ''
      end do
    end if

    write(MsgOut,'(a)') '--------------------------------------'

    return

  end subroutine output_control

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_pdb(pdbfile, prmtopfile, grocrdfile, pdb)

    ! formal arguments
    character(*),            intent(in)    :: pdbfile
    character(*),            intent(in)    :: prmtopfile
    character(*),            intent(in)    :: grocrdfile
    type(s_pdb),             intent(inout) :: pdb

    ! local variables
    type(s_prmtop)           :: prmtop
    type(s_grocrd)           :: grocrd
    integer                  :: i


    if (pdbfile == '' .and. prmtopfile == '' .and. grocrdfile == '') &
      return
      !call error_msg('Check_Reference> reference file is not given.')

    ! with pdb
    if (pdbfile /= '') then

      call input_pdb(pdbfile, pdb)

    end if

    ! with prmtop
    if (prmtopfile /= '') then

      call input_prmtop(prmtopfile, prmtop)

      call init_pdb(pdb)
      pdb%hetatm_rec = .true.
      pdb%cryst_rec  = .true.
      pdb%atom_col7  = .true.
      pdb%res_col5   = .true.
      pdb%num_atoms  = prmtop%num_atoms

      call alloc_pdb(pdb, PdbAtom, pdb%num_atoms)

      ! atom_name
      pdb%atom_name(:) = prmtop%atom_name(:)

      ! atom no
      do i = 1, pdb%num_atoms
        pdb%atom_no(i) = i
      end do

      ! residue name
      ! residue no
      do i = 1, size(prmtop%res_label)-1
        pdb%residue_name(prmtop%res_point(i):prmtop%res_point(i+1)-1) = &
             prmtop%res_label(i)
        pdb%residue_no  (prmtop%res_point(i):prmtop%res_point(i+1)-1) = &
             i
      end do
      i = size(prmtop%res_label)
      pdb%residue_name(prmtop%res_point(i):) = &
           prmtop%res_label(i)
      pdb%residue_no  (prmtop%res_point(i):) = &
           i

    end if

    ! with grocrd
    if (grocrdfile /= '') then

      call input_grocrd(grocrdfile, grocrd)

      call init_pdb(pdb)
      pdb%hetatm_rec = .true.
      pdb%cryst_rec  = .true.
      pdb%atom_col7  = .true.
      pdb%res_col5   = .true.
      pdb%num_atoms  = grocrd%num_atoms

      call alloc_pdb(pdb, PdbAtom, pdb%num_atoms)

      ! atom_no
      pdb%atom_no(:) = grocrd%atom_no(:)

      ! atom_name
      pdb%atom_name(:) = grocrd%atom_name(:)

      ! residue_no
      pdb%residue_no(:) = grocrd%residue_no(:)

      ! residue_name
      pdb%residue_name(:) = grocrd%residue_name(:)

    end if

    return

  end subroutine define_pdb

  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_distance(p1, p2)

    ! return values
    real(wp)                 :: get_distance

    ! formal arguments
    real(wp),                intent(in)    :: p1(:)
    real(wp),                intent(in)    :: p2(:)

    ! local variables
    real(wp)                 :: accum, d
    integer                  :: i
    

    accum = 0.0_wp

    do i = 1, size(p1)
      d = p2(i) - p1(i)
      accum = accum + d * d
    end do
    get_distance = sqrt(accum)

    return

  end function get_distance

  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_numbered_filename(basename, num)

    ! return value
    character(MaxFilename)   :: get_numbered_filename

    ! formal arguments
    character(*),            intent(in)    :: basename
    integer,                 intent(in)    :: num

    ! local variables
    integer                  :: lidx, ridx
    character(MaxFilename)   :: lstr, rstr


    lidx = index(basename, '{', back=.true.)
    ridx = index(basename, '}', back=.true.)

    if (lidx == 0 .or. ridx == 0) then
      get_numbered_filename = basename
      return
    end if

    lstr = basename(:lidx-1)
    rstr = basename(ridx+1:)

    write(get_numbered_filename, '(a,i0,a)') trim(lstr), num, trim(rstr)

    return

  end function get_numbered_filename

  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine sort_idxrep(a, start, end)

    ! formal arguments
    type(s_idxrep),          intent(inout) :: a(*)
    integer,                 intent(in)    :: start
    integer,                 intent(in)    :: end

    ! local variables
    type(s_idxrep)           :: t, x
    integer                  :: i, j


    x = a((start + end) / 2)
    i = start
    j = end

    do
      do while (a(i)%idx < x%idx)
        i = i + 1
      end do

      do while (x%idx < a(j)%idx)
        j = j - 1
      end do

      if (i >= j) exit

      t = a(i)
      a(i) = a(j)
      a(j) = t
      i = i + 1
      j = j - 1

    end do

    if (start < i - 1) &
      call sort_idxrep(a, start, i - 1)
    if (j + 1 < end) &
      call sort_idxrep(a, j + 1, end)

    return

  end subroutine sort_idxrep

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cumsum(array_in, array_out)

    ! formal arguments
    real(wp),                intent(in)    :: array_in(:)
    real(wp),                intent(inout) :: array_out(:)

    ! local variables
    integer                  :: i
    real(wp)                 :: accum


    accum = 0.0_wp

    do i = 1, size(array_in)
      accum = accum + array_in(i)
      array_out(i) = accum
    end do

    return

  end subroutine cumsum

  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_replicate_name1(filename, no)

    ! return
    character(Maxfilename)   :: get_replicate_name1

    ! formal arguments
    character(*),            intent(in)    :: filename
    integer,                 intent(in)    :: no

    ! local variables
    integer                  :: bl, br


    bl = index(filename, '{', back=.true.)
    br = index(filename, '}', back=.true.)

    if (bl == 0 .or. br == 0 .or. bl > br) &
      call error_msg('Get_Replicate_Name1> Syntax error.')

    write(get_replicate_name1, '(a,i0,a)') &
         filename(:bl-1),no,filename(br+1:len_trim(filename))

    return

  end function get_replicate_name1

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_cvfile(cvfile, nsteps)

    ! formal arguments
    character(*),            intent(in)    :: cvfile
    integer,                 intent(inout) :: nsteps

    ! local variables
    integer                  :: file
    character(MaxFilename)   :: filename
    character(MaxLine)       :: line


    !filename = get_replicate_name1(cvfile, 1)
    filename = cvfile

    call open_file(file, filename, IOFileInput)

    nsteps = 0
    do while(.true.)
      read(file,'(a)',end=10) line
      nsteps = nsteps + 1
    end do

10  call close_file(file)

    write(MsgOut,'(a,i12)') '  number of time steps       : ', nsteps
    write(MsgOut,'(a)')   ''

    return

  end subroutine check_cvfile

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_dcdfile(dcdfile, nsteps, natom)

    ! formal arguments
    character(*),            intent(in)    :: dcdfile
    integer,                 intent(inout) :: nsteps
    integer,                 intent(inout) :: natom

    ! local variables
    type(s_trj_file)         :: file
    integer                  :: i, hdr_size, step_size
    integer(8)               :: file_size
    integer(4)               :: icntrl(20), ntitle
    character(MaxFilename)   :: filename
    character(80)            :: title(10)
    character(4)             :: hdr
    logical                  :: exist

    integer(8)               :: ftell


    !filename = get_replicate_name1(dcdfile, 1)
    filename = dcdfile

    call open_trj(file, filename, TrjFormatDCD, TrjTypeCoorBox, IOFileInput)

    read(file%unit_no) hdr, icntrl(1:20)
    read(file%unit_no) ntitle,(title(i),i=1,ntitle)
    read(file%unit_no) natom

    ! check header size
    hdr_size = &
         4 + 4 + 20*4      + 4 + &  ! => read() hdr, icntrl
         4 + 4 + 80*ntitle + 4 + &  ! => read() ntitle, title(:)
         4 + 4             + 4      ! => read() natom

    ! check trajectory step size
    step_size = (4 + 4 * natom + 4) * 3

    ! check file size
#ifdef KCOMP
    inquire(file%unit_no,flen=file_size)
#else
    inquire(file%unit_no,size=file_size)
#endif

    if (mod(file_size - hdr_size, step_size) /= 0) &
      step_size = step_size + 4 + 8 * 6 + 4   ! box

    nsteps = (file_size - hdr_size) / step_size

    write(MsgOut,'(a,i12)') '  number of trajectory steps       : ', nsteps
    write(MsgOut,'(a)')   ''

    call close_trj(file)

    return

  end subroutine check_dcdfile

end module rg_generate_mod
