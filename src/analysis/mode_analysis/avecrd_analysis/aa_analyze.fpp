!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   aa_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module aa_analyze_mod

  use aa_option_str_mod
  use fitting_mod
  use fileio_trj_mod
  use fitting_str_mod
  use trajectory_str_mod
  use output_str_mod
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
  !! @param[inout] output     : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(molecule, trj_list, trajectory, fitting, option, output)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_trj_list),        intent(inout) :: trj_list
    type(s_trajectory),      intent(inout) :: trajectory
    type(s_fitting),         intent(inout) :: fitting
    type(s_option),          intent(inout) :: option
    type(s_output),          intent(inout) :: output

    ! local variables
    type(s_trj_file)         :: trj_in
    type(s_pdb)              :: pdb_out
    real(wp)                 :: nstru_inv
    integer                  :: natom, niter, ntraj, nstru
    integer                  :: iatom, iiter, itraj, istep
    integer                  :: rms_out, alloc_stat

    real(wp), allocatable  :: av0_coord(:,:)
    real(wp), allocatable  :: ave_coord(:,:)
    real(wp), allocatable  :: trj_coord(:,:)
    real(wp), allocatable  :: sqrt_mass(:)


    if (option%check_only) &
      return

    natom = molecule%num_atoms
    niter = option%num_iterations
    ntraj = size(trj_list%md_steps)


    ! allocate memory
    !
    allocate(sqrt_mass(natom),   &
             av0_coord(3,natom), &
             ave_coord(3,natom), &
             trj_coord(3,natom), stat=alloc_stat)
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

    ! setup initial average coordinates
    do iatom = 1, natom
      av0_coord(1:3,iatom) = molecule%atom_coord(1:3,iatom) * sqrt_mass(iatom)
    end do


    ! open output file
    !
    if (output%rmsfile /= '') &
      call open_file(rms_out, output%rmsfile, IOFileOutputNew)


    ! analysis loop
    !
    do iiter = 1, niter

      write(MsgOut,*) 'Analyze> number of iterations = ',iiter
      write(MsgOut,*) ' '

      ! reset average coordinates
      do iatom = 1, natom
        ave_coord(1:3,iatom) = 0.0_wp
      end do

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
                             av0_coord, &
                             trj_coord, &
                             trj_coord)

            if (output%rmsfile /= '') &
              call out_rmsd(rms_out, nstru, fitting)

            ! sum trajectory coordinates
            !
            do iatom = 1, natom
              ave_coord(1:3,iatom) = ave_coord(1:3,iatom) + trj_coord(1:3,iatom)
            end do

          end if

        end do

        call close_trj(trj_in)

      end do

      ! compute average coordinates
      nstru_inv = 1.0_wp / real(nstru, wp)

      do iatom = 1, natom
        ave_coord(1:3,iatom) = ave_coord(1:3,iatom) * nstru_inv
      end do

      ! execute fitting
      call run_fitting(fitting,   &
                       av0_coord, &
                       ave_coord, &
                       ave_coord)

      ! check convergence
      write(MsgOut,*) 'INFO> check the convergence: RMSD = ', fitting%rmsd
      write(MsgOut,*) ' '

      ! copy average coordinates as initial
      do iatom = 1, natom
        av0_coord(1:3,iatom) = ave_coord(1:3,iatom)
      end do

    end do


    ! output average coordinate data
    !
    do iatom = 1, natom
      molecule%atom_coord(1:3, iatom) = ave_coord(1:3,iatom) / sqrt_mass(iatom)
    end do

    if (output%pdb_avefile /= '') then
      call export_molecules(molecule, option%analysis_atom, pdb_out)
      call output_pdb(output%pdb_avefile, pdb_out)
      call dealloc_pdb_all(pdb_out)
    end if

    if (output%pdb_aftfile /= '') then
      call export_molecules(molecule, fitting%fitting_atom, pdb_out)
      call output_pdb(output%pdb_aftfile, pdb_out)
      call dealloc_pdb_all(pdb_out)
    end if


    ! close output file
    !
    if (output%rmsfile /= '') &
      call close_file(rms_out)


    ! Output summary
    !
    write(MsgOut,'(A)') 'Analyze> Detailed information in the output files'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '  [rmsfile] ' // trim(output%rmsfile)
    write(MsgOut,'(A)') '    Column 1: Snapshot index'
    write(MsgOut,'(A)') '    Column 2: Root-mean-square deviation (RMSD) with respect'
    write(MsgOut,'(A)') '              to the averaged structure (angstrom)'
    write(MsgOut,'(A)') ''


    ! deallocate memory
    !
    deallocate(av0_coord, ave_coord, trj_coord, sqrt_mass)

    return

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

end module aa_analyze_mod
