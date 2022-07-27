!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   rg_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Motoshi Kamiya (MK), Takaharu Mori (TM)
!
!  (c) Copyright 2016 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module rg_analyze_mod

  use rg_option_str_mod
  use fileio_trj_mod
  use measure_mod
  use trajectory_str_mod
  use output_str_mod
  use molecules_str_mod
  use select_atoms_mod
  use fileio_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public  :: analyze

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        run analyzing trajectories
  !! @authors      MK, TM
  !! @param[in]    molecule   : molecule information
  !! @param[in]    trj_list   : trajectory file list information
  !! @param[in]    output     : output information
  !! @param[inout] option     : option information
  !! @param[inout] trajectory : trajectory information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(molecule, trj_list, output, option, trajectory)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_trj_list),        intent(in)    :: trj_list
    type(s_output),          intent(in)    :: output
    type(s_option),          intent(inout) :: option
    type(s_trajectory),      intent(inout) :: trajectory

    ! local variables
    type(s_trj_file)          :: trj_in
    integer                   :: nstru, ifile, istep, num_trjfiles
    integer                   :: i, iatom, rg_out, idx
    real(wp)                  :: com(3), weight, tot_weight, rg


    if (option%check_only) &
      return


    ! open output file
    !
    if (output%rgfile /= '' ) &
      call open_file(rg_out, output%rgfile, IOFileOutputNew)


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

          ! compute rg
          !
          com(1:3) = 0.0_wp
          tot_weight = 0.0_wp
          do iatom = 1, size(option%analysis_atom%idx)

            idx = option%analysis_atom%idx(iatom)
            weight = 1.0_wp
            if (option%mass_weighted) weight = molecule%mass(idx)
            com(:) = com(:) + weight * trajectory%coord(:,idx)

            tot_weight = tot_weight + weight

          end do
          com(:) = com(:) / tot_weight

          rg = 0.0_wp
          do iatom = 1, size(option%analysis_atom%idx)

            idx = option%analysis_atom%idx(iatom)
            weight = 1.0_wp
            if (option%mass_weighted) weight = molecule%mass(idx)

            do i = 1, 3
              rg = rg + weight * ( ( trajectory%coord(i,idx) - com(i) ) ** 2 )
            end do
          end do
          rg = sqrt( rg / tot_weight )

          ! output results
          !
          write(MsgOut,'(a,f10.5)') '              RG of analysis atoms = ',rg

          write(MsgOut,*) ''

          if (output%rgfile /= '') &
            write(rg_out, '(i10,1x,f10.5)') nstru, rg

        end if

      end do

      ! close trajectory file
      !
      call close_trj(trj_in)

    end do

    ! close output file
    !
    if (output%rgfile /= '') call close_file(rg_out)


    ! Output summary
    !
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze> Detailed information in the output files'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '  [rgfile] ' // trim(output%rgfile)
    write(MsgOut,'(A)') '    Column 1: Snapshot index'
    write(MsgOut,'(A)') '    Column 2: Radius of gyration (angstrom)'
    write(MsgOut,'(A)') ''

    return

  end subroutine analyze

end module rg_analyze_mod
