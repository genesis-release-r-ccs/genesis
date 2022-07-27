!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ca_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Takaharu Mori (TM), Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ca_analyze_mod

  use ca_option_str_mod
  use fileio_trj_mod
  use fitting_mod
  use fitting_str_mod
  use pbc_correct_mod
  use measure_mod
  use trajectory_str_mod
  use output_str_mod
  use molecules_str_mod
  use select_atoms_mod
  use select_atoms_str_mod
  use fileio_mod
  use messages_mod
  use constants_mod
 
  implicit none
  private


  ! subroutines
  public  :: analyze
  private :: define_molecule_range_mass
  private :: analyze_com

  real(wp), allocatable, save :: molecule_mass(:)
  logical,  allocatable, save :: stop_accumulation(:)

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        run analyzing trajectories
  !! @authors      TM
  !! @param[in]    trj_list   : trajectory file list information
  !! @param[in]    output     : output information
  !! @param[inout] option     : option information
  !! @param[inout] trajectory : trajectory information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(molecule, trj_list, output, fitting, option, trajectory)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_trj_list),        intent(in)    :: trj_list
    type(s_output),          intent(in)    :: output
    type(s_fitting),         intent(inout) :: fitting
    type(s_option),          intent(inout) :: option
    type(s_trajectory),      intent(inout) :: trajectory

    ! local variables
    type(s_trj_file)         :: trj_in
    integer                  :: nstru, ifile, istep, num_trjfiles
    integer                  :: iatom, trj_out, idx


    if (option%check_only) &
      return


    ! define analysis molecule range
    !
    call define_molecule_range_mass(option, molecule)


    ! open output file
    !
    if (output%trjfile /= '') &
      call open_file(trj_out, output%trjfile, IOFileOutputNew)


    ! analysis loop
    !
    nstru = 0
    num_trjfiles = size(trj_list%md_steps)

    do ifile = 1, num_trjfiles

      ! open trajectory file
      !
      call open_trj(trj_in, trj_list%filenames(ifile), &
                            trj_list%trj_format,       &
                            trj_list%trj_type, IOFileInput)

      do istep = 1, trj_list%md_steps(ifile)

        ! read trajectory
        !   coordinates of one MD snapshot are saved in trajectory%coord)
        !
        call read_trj(trj_in, trajectory)

        if (mod(istep, trj_list%ana_periods(ifile)) == 0) then

          nstru = nstru + 1
          write(MsgOut,*) '      number of structures = ', nstru

          ! pbc-correct
          !
          call run_pbc_correct(option%pbcc_mode, molecule, trajectory)

          ! fitting
          !
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

          ! calculate COM
          !
          call analyze_com(nstru, trj_out, option, molecule, trajectory)

        end if

      end do

      ! close trajectory file
      !
      call close_trj(trj_in)

    end do


    ! close output file
    !
    if (output%trjfile /= '') call close_file(trj_out)


    ! Output summary
    !
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze> Detailed information in the output files'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '  [trjfile] ' // trim(output%trjfile)
    write(MsgOut,'(A)') '    Column  1 : Snapshot index'
    write(MsgOut,'(A)') '    Columns 2-: Coordinates of the center of mass (COM) of the selected group'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '    1) analysis_for = ALL, output_coord = XYZ:'
    write(MsgOut,'(A)') '      Columns 2-4: X, Y, Z of the COM of the selected group'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '    2) analysis_for = ALL, output_coord = XY:'
    write(MsgOut,'(A)') '      Columns 2-3: X, Y of the COM of the selected group'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '    3) analysis_for = MOLECULE, output_coord = XYZ:'
    write(MsgOut,'(A)') '      Columns 2-4: X, Y, Z of the COM of the 1st molecule in the selected group'
    write(MsgOut,'(A)') '      Columns 5-7: X, Y, Z of the COM of the 2nd molecule in the selected group'
    write(MsgOut,'(A)') '               :                              :'
    write(MsgOut,'(A)') ''

    return

  end subroutine analyze

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_molecule_range_mass
  !> @brief        define molecule range and mass
  !! @authors      TM
  !! @param[in]    molecule   : molecule information
  !! @param[inout] trajectory : trajectory information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_molecule_range_mass(option, molecule)

    ! formal arguments
    type(s_option),             intent(inout) :: option
    type(s_molecule),           intent(inout) :: molecule

    ! local variables
    real(wp)             :: mol_mass
    integer              :: i, atom_index, natoms, nmols, imol
    integer              :: current_mol_id, next_mol_id
    integer              :: mol_first, mol_last

    natoms = size(option%analysis_atom%idx)

    allocate(stop_accumulation(natoms))


    select case (option%analysis_for)

    case (AnalysisForMolecule)

      nmols = 0
      do i = 1, natoms
        if (i /= natoms) then
          current_mol_id = molecule%molecule_no(option%analysis_atom%idx(i))
          next_mol_id    = molecule%molecule_no(option%analysis_atom%idx(i+1))

          if (current_mol_id == next_mol_id) then
            stop_accumulation(i) = .false.
          else
            stop_accumulation(i) = .true.
            nmols = nmols + 1
          end if
        else
          stop_accumulation(i) = .true.
          nmols = nmols + 1
        end if
      end do

    case (AnalysisForAll)

      nmols = 1
      stop_accumulation(1:natoms) = .false.
      stop_accumulation(natoms)   = .true.

    end select


    allocate(molecule_mass(nmols))

    imol = 0
    mol_mass = 0.0_wp
    do i = 1, natoms
      atom_index = option%analysis_atom%idx(i)
      mol_mass   = mol_mass + molecule%mass(atom_index)

      if (stop_accumulation(i)) then
        imol = imol + 1
        molecule_mass(imol) = mol_mass
        mol_mass = 0.0_wp
      end if
    end do


    ! write setup results
    !
    write(MsgOut,'(a)') '     molecule   first atom    last atom     molecule mass'

    imol = 0
    mol_first = option%analysis_atom%idx(1)

    do i = 1, natoms
      atom_index = option%analysis_atom%idx(i)
      mol_last   = atom_index

      if (stop_accumulation(i)) then
        imol = imol + 1
        write(MsgOut,'(3i13,f18.4)') imol, mol_first, mol_last, molecule_mass(imol)
        if (i /= natoms) then
          mol_first = option%analysis_atom%idx(i+1)
        end if
      end if

    end do

    write(MsgOut,'(a)') ''


    return

  end subroutine define_molecule_range_mass

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze_com_by_molecule
  !> @brief        analyze com of each molecule
  !! @authors      TM, KY
  !! @param[in]    molecule   : molecule information
  !! @param[inout] trajectory : trajectory information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze_com(nstru, trj_out, option, molecule, trajectory)

    ! formal arguments
    integer,                    intent(inout) :: nstru
    integer,                    intent(inout) :: trj_out
    type(s_option),             intent(inout) :: option
    type(s_molecule),           intent(inout) :: molecule
    type(s_trajectory), target, intent(inout) :: trajectory

    ! local variables
    real(wp)             :: com(3), mol_mass
    integer              :: i, k, atom_index, natoms, nmols, imol
    integer              :: j, jst, close_index, kk
    real(wp)             :: dd, dj
    integer              :: current_mol_id, next_mol_id
    real(wp), pointer    :: coord(:,:)
    character            :: frmt1*10, frmt2*8


    coord  => trajectory%coord
    natoms =  size(option%analysis_atom%idx)

    ! calculate com of each molecule
    !
    k = option%output_coord
    write(frmt1,'(a,i0,a)' ) '(',k,'F10.3,$)'
    write(frmt2,'(a,i0,a)' ) '(',k,'F10.3)'

    write(trj_out,'(I10,$)') nstru

    imol = 0
    jst  = 1
    do i = 1, natoms

      atom_index = option%analysis_atom%idx(i)
      com(1:k)   = com(1:k) + coord(1:k,atom_index) * molecule%mass(atom_index)

      if (stop_accumulation(i)) then
        imol = imol + 1
        com(1:k) = com(1:k)/molecule_mass(imol)
        write(trj_out,frmt1) com(1:k)

        if (option%output_atomno) then
          dd = 100.0_wp
          close_index = 0
          do j = jst, i

            atom_index = option%analysis_atom%idx(j)
            dj = 0.0_wp
            do kk = 1, k
              dj = dj + (com(kk) - coord(kk,atom_index))**2
            end do
            dj = sqrt(dj)

            if (dj < dd) then
              dd = dj
              close_index = j
            end if

          end do
          write(trj_out,'(2x,i8,$)') close_index

        end if

        com(1:k) = 0.0_wp
        if (i == natoms) write(trj_out,*)

      end if

    end do

    return

  end subroutine analyze_com

end module ca_analyze_mod
