!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   qa_cg_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Daisuke Matsuoka (DM), Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module qa_cg_analyze_mod

  use qa_cg_option_str_mod
  use fileio_trj_mod
  use trajectory_str_mod
  use output_str_mod
  use molecules_str_mod
  use fileio_mod
  use measure_mod
  use messages_mod
  use constants_mod
  use atom_libs_mod

  implicit none
  private

  ! subroutines
  public  :: analyze
  private :: analyze_qnt
  private :: out_result

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        run analyzing trajectories
  !! @authors      DM, NT
  !! @param[in]    trj_list   : trajectory file list information
  !! @param[in]    output     : output information
  !! @param[inout] option     : option information
  !! @param[inout] contact    : contact information
  !! @param[inout] trajectory : trajectory information
  !! @note         PNAS (2013) 110, 17874 - 17879
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(molecule, trj_list, output, option, contact, trajectory)

    ! formal arguments
    type(s_molecule),   intent(in)    :: molecule
    type(s_trj_list),   intent(in)    :: trj_list
    type(s_output),     intent(in)    :: output
    type(s_option),     intent(inout) :: option
    type(s_contact),    intent(inout) :: contact
    type(s_trajectory), intent(inout) :: trajectory

    ! local variables
    type(s_trj_file)          :: trj_in
    integer                   :: nstru, ifile, istep, num_trjfiles, qnt_unit
    integer                   :: n_pair
    real(wp)                  :: qval


    if (option%check_only) &
      return

    ! open output file
    !
    if (output%qntfile /= '') &
      call open_file(qnt_unit, output%qntfile, IOFileOutputNew)

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

          call analyze_qnt(trajectory, option, contact, qval)
          call out_result (nstru, qnt_unit, qval)

        end if

      end do

      ! close trajectory file
      !
      call close_trj(trj_in)

    end do


    ! close output file
    !
    call close_file(qnt_unit)

    ! Output summary
    !
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze> Detailed information in the output files'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '  [qntfile] ' // trim(output%qntfile)
    write(MsgOut,'(A)') '    Column 1: Snapshot index'
    write(MsgOut,'(A)') '    Column 2: Fraction of the native contact'
    write(MsgOut,'(A)') ''

    ! deallocate(contact%cnt_pair, contact%r0_ij)

    return

  end subroutine analyze


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze_qnt
  !> @brief        analysis fraction of natice contact
  !! @authors      DM
  !! @param[in]    trajectory : trajectory informatin
  !! @param[in]    option     : option informatin
  !! @param[in]    contact    : contact information
  !! @param[out]   qval       : Q-value
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze_qnt(trajectory, option, contact, qval)

    ! formal argments
    type(s_trajectory),   intent(in)    :: trajectory
    type(s_option),       intent(in)    :: option
    type(s_contact),      intent(in)    :: contact
    real(wp),             intent(out)   :: qval

    ! local variables
    integer  :: i_pair
    integer  :: idx, jdx
    real(wp) :: dist
    real(wp) :: q_total


    ! examine contact
    !
    q_total = 0.0_wp

    do i_pair = 1, contact%n_pair
      idx = contact%cnt_pair(1, i_pair)
      jdx = contact%cnt_pair(2, i_pair)

      dist = compute_dis(trajectory%coord(:, idx), &
                         trajectory%coord(:, jdx))

      if ( dist < option%lambda * contact%r0_ij(i_pair) ) then
        q_total = q_total + 1.0_wp
      end if

    end do

    qval = q_total / real(contact%n_pair, wp)

    return

  end subroutine analyze_qnt


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    out_result
  !> @brief        output analysis results
  !! @authors      NT
  !! @param[in]    sturct_no : sturcture number
  !! @param[in]    unit_no   : file unit number
  !! @param[in]    qval      : analysis results
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine out_result(struct_no, unit_no, qval)

    ! formal arguments
    integer,                 intent(in)    :: struct_no
    integer,                 intent(in)    :: unit_no
    real(wp),                intent(in)    :: qval


    write(unit_no, '(i10,1x,f8.5)') struct_no, qval

    return

  end subroutine out_result

end module qa_cg_analyze_mod
