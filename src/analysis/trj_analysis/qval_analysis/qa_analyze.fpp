!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   qa_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Daisuke Matsuoka (DM), Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module qa_analyze_mod

  use qa_option_str_mod
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

  ! structure
  type, private :: s_contact
    integer               :: n_pair
    real(wp), allocatable :: r0_ij(:)
    integer,  allocatable :: cnt_pair(:,:)
  end type s_contact

  ! subroutines
  public  :: analyze
  private :: count_contact_pair
  private :: native_contact
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
  !! @param[inout] trajectory : trajectory information
  !! @note         PNAS (2013) 110, 17874 - 17879
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(molecule, trj_list, output, option, trajectory)

    ! formal arguments
    type(s_molecule),   intent(in)    :: molecule
    type(s_trj_list),   intent(in)    :: trj_list
    type(s_output),     intent(in)    :: output
    type(s_option),     intent(inout) :: option
    type(s_trajectory), intent(inout) :: trajectory

    ! local variables
    type(s_trj_file)          :: trj_in
    integer                   :: nstru, ifile, istep, num_trjfiles, qnt_unit
    integer                   :: n_pair
    type(s_contact)           :: contact
    real(wp)                  :: qval


    if (option%check_only) &
      return

    ! identify native contact
    !
    call count_contact_pair(molecule, option, n_pair)

    contact%n_pair = n_pair
    allocate(contact%r0_ij(n_pair), contact%cnt_pair(2, n_pair))

    call native_contact(molecule, option, contact)

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

    deallocate(contact%cnt_pair, contact%r0_ij)

    return

  end subroutine analyze

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    count_contact_pair
  !> @brief        analysis natice contact
  !! @authors      DM
  !! @param[in]    molecule   : molecule informatin
  !! @param[inout] option     : option information
  !! @param[out]   i_pair     : number of the native contact pair
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine count_contact_pair(molecule, option, i_pair)

    ! formal argments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_option),           intent(in)    :: option
    integer,                  intent(out)   :: i_pair

    ! local variables
    integer,          pointer :: numa(:)
    integer,          pointer :: numr(:)
    character(len=4), pointer :: nama(:)
    character(len=6), pointer :: namr(:)
    character(len=4), pointer :: seg(:)
    character(len=1), pointer :: chn(:)

    real(wp)                  :: dist
    integer                   :: natm
    integer                   :: iatm, jatm
    integer                   :: idx, jdx

    ! format
100 format(i8,' | ',i8,a5,a7,i5,1x,a1,1x,a5' ... ',i8,a5,a7,i5,1x,a1,1x,a5' | ',f5.3)

    ! molecule information
    !
    numa => molecule%atom_no
    nama => molecule%atom_name
    numr => molecule%residue_no
    namr => molecule%residue_name
    seg  => molecule%segment_name
    chn  => molecule%chain_id

    ! examine native contact
    !
    if (option%verbose) then
      write(MsgOut,'(A)') '[native contact list]'
    end if

    i_pair = 0

    natm = size(option%analysis_atom%idx)
    do iatm = 1, natm
      idx = option%analysis_atom%idx(iatm)

      do jatm = iatm + 1, natm
        jdx = option%analysis_atom%idx(jatm)

        if ((seg(idx) .ne. seg(jdx)) .or.  &
            (chn(idx) .ne. chn(jdx))) then

          dist = compute_dis(molecule%atom_coord(:,idx), &
                             molecule%atom_coord(:,jdx))

          if (dist >= option%maximum_distance) cycle

          i_pair = i_pair + 1

          if (option%verbose) then
            write(MsgOut, 100) i_pair,  &
                               numa(idx), nama(idx), namr(idx), numr(idx), chn(idx), seg(idx), &
                               numa(jdx), nama(jdx), namr(jdx), numr(jdx), chn(jdx), seg(jdx), &
                               dist
          end if

        else if (abs(numr(idx) - numr(jdx)) > 3) then
          dist = compute_dis(molecule%atom_coord(:,idx), &
                              molecule%atom_coord(:,jdx))

          if (dist >= option%maximum_distance) cycle

          i_pair = i_pair + 1

          if (option%verbose) then
            write(MsgOut, 100) i_pair,  &
                               numa(idx), nama(idx), namr(idx), numr(idx), chn(idx), seg(idx),  &
                               numa(jdx), nama(jdx), namr(jdx), numr(jdx), chn(jdx), seg(jdx),  &
                               dist
          end if
        end if
      end do
    end do

    if (option%verbose) then
      write(MsgOut,'(A)') '[end of native contact list]'
      write(MsgOut,'(A)') ''

    else
      write(MsgOut,'(A, i0)') 'number of native contact is ', i_pair
    end if

    nullify(numa, nama, numr, namr, seg, chn)

    return

  end subroutine count_contact_pair

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    native_contact
  !> @brief        analysis natice contact
  !! @authors      DM
  !! @param[in]    molecule   : molecule informatin
  !! @param[inout] option     : option information
  !! @param[out]   contact    : contact information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine native_contact(molecule, option, contact)

    ! formal argments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_option),           intent(in)    :: option
    type(s_contact),  target, intent(inout) :: contact

    ! local variables
    integer,          pointer :: numr(:)
    character(len=4), pointer :: seg(:)
    character(len=1), pointer :: chn(:)

    integer,          pointer :: n_pair
    real(wp),         pointer :: r0_ij(:)
    integer,          pointer :: cnt_pair(:,:)

    real(wp)                  :: dist
    integer                   :: natm
    integer                   :: iatm, jatm
    integer                   :: idx, jdx
    integer                   :: i_pair

    ! molecule information
    !
    numr => molecule%residue_no
    seg  => molecule%segment_name
    chn  => molecule%chain_id

    ! contact information
    n_pair   => contact%n_pair
    r0_ij    => contact%r0_ij
    cnt_pair => contact%cnt_pair


    ! obratin native contact pair information
    !
    i_pair = 0

    natm = size(option%analysis_atom%idx)
    do iatm = 1, natm
      idx = option%analysis_atom%idx(iatm)

      do jatm = iatm + 1, natm
        jdx = option%analysis_atom%idx(jatm)

        if ((seg(idx) .ne. seg(jdx)) .or.  &
            (chn(idx) .ne. chn(jdx))) then

          dist = compute_dis(molecule%atom_coord(:,idx), &
                             molecule%atom_coord(:,jdx))

          if (dist >= option%maximum_distance) cycle

          i_pair = i_pair + 1
          r0_ij(i_pair)       = dist
          cnt_pair(1, i_pair) = idx
          cnt_pair(2, i_pair) = jdx

        else if (abs(numr(idx) - numr(jdx)) > 3) then
          dist = compute_dis(molecule%atom_coord(:,idx), &
                             molecule%atom_coord(:,jdx))

          if (dist >= option%maximum_distance) cycle

          i_pair = i_pair + 1
          r0_ij(i_pair)       = dist
          cnt_pair(1, i_pair) = idx
          cnt_pair(2, i_pair) = jdx

        end if
      end do
    end do

    nullify(numr, seg, chn)
    nullify(n_pair, r0_ij, cnt_pair)

    return

  end subroutine native_contact

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
    real(wp) :: diff
    real(wp) :: qij
    real(wp) :: q_total


    ! examine contact
    !
    q_total = 0.0_wp

    do i_pair = 1, contact%n_pair
      idx = contact%cnt_pair(1, i_pair)
      jdx = contact%cnt_pair(2, i_pair)

      dist = compute_dis(trajectory%coord(:, idx), &
                         trajectory%coord(:, jdx))

      diff = option%beta * (dist - option%lambda * contact%r0_ij(i_pair))
      qij  = 1.0_wp / (1.0_wp + exp(diff))

      q_total = q_total + qij
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

end module qa_analyze_mod
