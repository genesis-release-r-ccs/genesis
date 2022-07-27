!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pa_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module pa_analyze_mod

  use pa_option_str_mod
  use fitting_mod
  use fileio_trj_mod
  use fitting_str_mod
  use trajectory_str_mod
  use output_str_mod
  use input_str_mod
  use molecules_mod
  use molecules_str_mod
  use fileio_pdb_mod
  use fileio_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public  :: analyze
  private :: assign_mass
  private :: sum_project
  private :: out_project

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        run analyzing trajectories
  !! @authors      NT
  !! @param[inout] molecule   : molecule information
  !! @param[inout] trj_list   : trajectory file list information
  !! @param[inout] trajectory : trajectory information
  !! @param[inout] fitting    : fitting information
  !! @param[inout] option     : option information
  !! @param[inout] input      : input information
  !! @param[inout] output     : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(molecule, trj_list, trajectory, fitting, option, &
                     input, output)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_trj_list),        intent(inout) :: trj_list
    type(s_trajectory),      intent(inout) :: trajectory
    type(s_fitting),         intent(inout) :: fitting
    type(s_option),          intent(inout) :: option
    type(s_input),           intent(inout) :: input
    type(s_output),          intent(inout) :: output

    ! local variables
    type(s_pdb)              :: pdb_in
    type(s_trj_file)         :: trj_in
    real(wp)                 :: nstru_r, ttl_mass, rmsf
    integer                  :: natom, ntraj, nsela, nself, nstru
    integer                  :: iatom, itraj, isel, istep, nvcv
    integer                  :: val_in, vec_in, prj_out
    integer                  :: i, j, alloc_stat

    real(wp),   allocatable  :: sqrt_mass (:)
    real(wp),   allocatable  :: ave_coord (:,:)
    real(wp),   allocatable  :: aft_coord (:,:)
    real(wp),   allocatable  :: trj_coord (:,:)
    real(wp),   allocatable  :: eig_value (:)
    real(wp),   allocatable  :: eig_vector(:,:)
    real(wp),   allocatable  :: vector(:)
    real(wp),   allocatable  :: project(:)


    if (option%check_only) &
      return

    natom = molecule%num_atoms
    ntraj = size(trj_list%md_steps)
    nsela = size(option%analysis_atom%idx)
    nself = size(fitting%fitting_atom%idx)
    if (option%vcv_matrix == VcvMatrixGlobal) then
      nvcv = nsela * 3
    else
      nvcv = nself * 3
    end if


    ! allocate memory
    !
    allocate(sqrt_mass (natom),      &
             ave_coord (3,natom),    &
             aft_coord (3,natom),    &
             trj_coord (3,natom),    &
             eig_value (nvcv),       &
             eig_vector(nvcv, nvcv), &
             vector    (nvcv),       &
             project   (option%num_pca), stat=alloc_stat)

    if (alloc_stat /= 0) &
      call error_msg_alloc


    ! prepare data
    !

    ! setup mass
    if (fitting%mass_weight) then
      call assign_mass(molecule)
    end if

    if (fitting%mass_weight) then
      do iatom = 1, natom
        sqrt_mass(iatom) = sqrt(molecule%mass(iatom))
      end do
    else
      do iatom = 1, natom
        sqrt_mass(iatom) = 1.0_wp
      end do
    end if

    ! setup average coordintes
    call input_pdb(input%pdb_avefile, pdb_in)

    if (nsela /= size(pdb_in%atom_name)) &
      call error_msg('Analyze> different count between pdb_avefile atom '// &
                     'and analysis_atom selection.')

    do isel = 1, nsela
      iatom = option%analysis_atom%idx(isel)
      ave_coord(1:3,iatom) = pdb_in%atom_coord(1:3,isel) * sqrt_mass(iatom)
    end do

    call dealloc_pdb_all(pdb_in)

    ! setup fitted average coordinates
    call input_pdb(input%pdb_aftfile, pdb_in)

    if (nself /= size(pdb_in%atom_name)) &
      call error_msg('Analyze> different count between pdb_aftfile atom '// &
                     'and fitting_atom selection.')

    do isel = 1, nself
      iatom = fitting%fitting_atom%idx(isel)
      aft_coord(1:3,iatom) = pdb_in%atom_coord(1:3,isel) * sqrt_mass(iatom)
    end do

    call dealloc_pdb_all(pdb_in)

    ! setup eigen-vector and eigen-val
    call open_file(val_in, input%valfile, IOFileInput)
    call open_file(vec_in, input%vecfile, IOFileInput)

    do i = 1, nvcv
      read(val_in,*,err=900,end=900) eig_value(i)
      do j = 1, nvcv
        read(vec_in,*,err=900,end=900) eig_vector(j,i)
      end do
    end do

    call close_file(val_in)
    call close_file(vec_in)

    ! setup total mass
    if (fitting%mass_weight) then
      ttl_mass = 0.0_wp
      do  isel = 1, nsela
        iatom = option%analysis_atom%idx(isel)
        ttl_mass = ttl_mass + molecule%mass(iatom)
      end do
    else
      ttl_mass = 1.0_wp
    end if

    ! setup project
    do i = 1, option%num_pca
      project(i) = 0.0_wp
    end do


    ! compute rmsf
    !
    write(MsgOut,'(A)') 'Analyze> rmsf (A)'

    do i = 1, option%num_pca
      rmsf = sqrt(eig_value(i) / ttl_mass)
      write(MsgOut,*) i, rmsf
    end do
    write(MsgOut,*)    ' '


    ! open output file
    !
    call open_file(prj_out, output%prjfile, IOFileOutputNew)


    ! analysis loop
    !
    nstru = 0

    do itraj = 1, ntraj

      call open_trj(trj_in, &
                    trj_list%filenames(itraj), &
                    trj_list%trj_format,       &
                    trj_list%trj_type, IOFileInput)

      do istep = 1, trj_list%md_steps(itraj)

        call read_trj(trj_in, trajectory)

        if (mod(istep, trj_list%ana_periods(itraj)) == 0) then

          nstru = nstru + 1
          write(MsgOut,*) '      number of structures = ', nstru


          ! geometrical fitting
          !
          call run_fitting(fitting,             &
                           molecule%atom_coord, &
                           trajectory%coord,    &
                           trajectory%coord)

          ! mass-weighted fitting
          !
          do iatom = 1, natom
            trj_coord(1:3,iatom) = trajectory%coord(1:3,iatom)*sqrt_mass(iatom)
          end do

          call run_fitting(fitting,   &
                           aft_coord, &
                           trj_coord, &
                           trj_coord)

          if (option%vcv_matrix == VcvMatrixGlobal) then

            call sum_project(option%analysis_atom%idx, &
                             option%num_pca, ttl_mass, &
                             trj_coord, ave_coord,     &
                             eig_vector, vector, project)
          else

            call sum_project(fitting%fitting_atom%idx, &
                             option%num_pca, ttl_mass, &
                             trj_coord, aft_coord,     &
                             eig_vector, vector, project)
          end if

          call out_project(prj_out, option%num_pca, nstru, &
                             project)

        end if

      end do

      call close_trj(trj_in)

    end do


    ! close output file
    !
    call close_file(prj_out)

    ! Output summary
    !
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze> Detailed information in the output files'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '  [prjfile] ' // trim(output%prjfile)
    write(MsgOut,'(A)') '    Column 1: Snapshot index'
    write(MsgOut,'(A)') '    Column 2: Projected coordinates on the PC1 axis'
    write(MsgOut,'(A)') '       :                    :'
    write(MsgOut,'(A)') '    Column N: Projected coordinates on the PC(N-1) axis'
    write(MsgOut,'(A)') ''

    ! deallocate memory
    !
    deallocate(project, vector, eig_vector, eig_value, &
               trj_coord, aft_coord, ave_coord, sqrt_mass)

    return

900 call error_msg('Analyze> mismatch VEC/ VAL file from atom selection count')

  end subroutine analyze

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_mass
  !> @brief        assign mass
  !! @authors      NT
  !! @param[inout] molecule   : molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_mass(molecule)

    ! parameters
    real(wp),                parameter     :: MassH   =  1.008000_wp
    real(wp),                parameter     :: MassC   = 12.011000_wp
    real(wp),                parameter     :: MassN   = 14.007000_wp
    real(wp),                parameter     :: MassO   = 15.999000_wp
    real(wp),                parameter     :: MassS   = 32.060000_wp
    real(wp),                parameter     :: MassP   = 30.974000_wp
    real(wp),                parameter     :: MassMG  = 24.305000_wp
    real(wp),                parameter     :: MassZN  = 65.370000_wp

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule

    ! local variables
    integer                  :: i


    write(MsgOut,'(a)'),      'WARNING: atom mass is not assigned.'
    write(MsgOut,'(a)'),      '   uses default mass.'
    write(MsgOut,'(a,f9.6)')  '      1) H   : ', MassH
    write(MsgOut,'(a,f9.6)')  '      2) C   : ', MassC
    write(MsgOut,'(a,f9.6)')  '      3) N   : ', MassN
    write(MsgOut,'(a,f9.6)')  '      4) O   : ', MassO
    write(MsgOut,'(a,f9.6)')  '      5) S   : ', MassS
    write(MsgOut,'(a,f9.6)')  '      6) P   : ', MassP
    write(MsgOut,'(a,f9.6)')  '      7) MG  : ', MassMG
    write(MsgOut,'(a,f9.6)')  '      8) ZN  : ', MassZN


    do i = 1, molecule%num_atoms

      if (molecule%atom_name(i)(1:1) == 'H') then
        molecule%mass(i) = MassH
      else if (molecule%atom_name(i)(1:1) == 'C') then
        molecule%mass(i) = MassC
      else if (molecule%atom_name(i)(1:1) == 'N') then
        molecule%mass(i) = MassN
      else if (molecule%atom_name(i)(1:1) == 'O') then
        molecule%mass(i) = MassO
      else if (molecule%atom_name(i)(1:1) == 'S') then
        molecule%mass(i) = MassS
      else if (molecule%atom_name(i)(1:1) == 'P') then
        molecule%mass(i) = MassP
      else if (molecule%atom_name(i)(1:2) == 'MG') then
        molecule%mass(i) = MassMG
      else if (molecule%atom_name(i)(1:2) == 'ZN') then
        molecule%mass(i) = MassZN
      else
        write(MsgOut,'(a,a)') 'Assign_Mass> Unknown atom :', &
             molecule%atom_name(i)
      end if

    end do

    write(MsgOut,'(a)') ''

    return

  end subroutine assign_mass

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    sum_project
  !> @brief        sum project
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine sum_project(idx, num_pca, ttl_mass, trj_coord, ave_coord, &
                         eig_vector, vector, project)

    ! formal arguments
    integer,                 intent(in)    :: idx(:)
    integer,                 intent(in)    :: num_pca
    real(wp),                intent(in)    :: ttl_mass
    real(wp),                intent(in)    :: trj_coord(:,:)
    real(wp),                intent(in)    :: ave_coord(:,:)
    real(wp),                intent(in)    :: eig_vector(:,:)
    real(wp),                intent(inout) :: vector(:)
    real(wp),                intent(inout) :: project(:)

    ! local variables
    integer                  :: i, j, k, iatom, nsel, nvcv


    nsel = size(idx)
    nvcv  = nsel * 3

    do i = 1, nsel
      iatom = idx(i)
      vector(3*(i-1)+1) = trj_coord(1,iatom) - ave_coord(1,iatom)
      vector(3*(i-1)+2) = trj_coord(2,iatom) - ave_coord(2,iatom)
      vector(3*(i-1)+3) = trj_coord(3,iatom) - ave_coord(3,iatom)
    end do

    do i = 1, num_pca
      project(i) = 0.0_wp
      ! j = nvcv - i + 1
      do k = 1, nvcv
        project(i) = project(i) + vector(k) * eig_vector(k,i)
      end do
      project(i) = project(i) / sqrt(ttl_mass)
    end do

    return

  end subroutine sum_project

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    out_project
  !> @brief        out project
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine out_project(prj_out, num_pca, nstru, project)

    ! formal arguments
    integer,                 intent(in)    :: prj_out
    integer,                 intent(in)    :: num_pca
    integer,                 intent(in)    :: nstru
    real(wp),                intent(in)    :: project(:)
    
    ! local variables
    integer                  :: i


    write(prj_out, '(i9,100(f10.5,1x))') nstru, (project(i),i=1,num_pca)
    write(MsgOut,  '(i9,100(f10.5,1x))') nstru, (project(i),i=1,num_pca)

    return

  end subroutine out_project

end module pa_analyze_mod
