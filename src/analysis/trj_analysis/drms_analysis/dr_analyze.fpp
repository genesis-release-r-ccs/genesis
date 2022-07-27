!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   dr_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Chigusa Kobayashi (CK), Daisuke Matsuoka (DM), Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module dr_analyze_mod

  use dr_option_str_mod
  use fileio_trj_mod
  use fitting_mod
  use fitting_str_mod
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
  private :: compute_drms

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        run analyzing trajectories
  !! @authors      CK, DM, NT
  !! @param[in]    trj_list   : trajectory file list information
  !! @param[in]    output     : output information
  !! @param[inout] option     : option information
  !! @param[inout] trajectory : trajectory information
  !! @note         JPCB (2017) 121, 3364 - 3375
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
    type(s_trj_file)         :: trj_in
    integer                  :: nstru, ifile, istep, num_trjfiles
    integer                  :: rms_out
    real(wp)                 :: drms, drms_cur


    if (option%check_only) &
      return

    ! open output file
    !
    if (output%rmsfile /= '') &
      call open_file(rms_out, output%rmsfile, IOFileOutputNew)

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

          if (option%two_states) then
            call compute_drms_two_states(option, trajectory%coord,  &
                                         trajectory%pbc_box, drms, drms_cur)

            if (output%rmsfile /= '') &
              write(rms_out, '(i10,1x,2f8.3)') nstru, drms, drms_cur
          else
            call compute_drms(option, trajectory%coord, trajectory%pbc_box, drms)

            if (output%rmsfile /= '') &
              write(rms_out, '(i10,1x,f8.3)') nstru, drms
          endif

        end if

      end do

      ! close trajectory file
      !
      call close_trj(trj_in)

    end do


    ! close output file
    !
    if (output%rmsfile /= '') call close_file(rms_out)


    ! Output summary
    !
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze> Detailed information in the output files'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '  [rmsfile] ' // trim(output%rmsfile)
    write(MsgOut,'(A)') '    Column 1: Snapshot index'
    write(MsgOut,'(A)') '    Column 2: Distance RMSD (angstrom)'
    write(MsgOut,'(A)') ''

    return

  end subroutine analyze

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine  compute_drms
  !> @brief      calculate drms
  !! @authors    CK
  !! @param[in]  enefunc potential energy functions [str]
  !! @param[in]  coord coordinates of target systems [dble]
  !! @date       2018/07/19 (CK)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_drms(option, coord, pbc_box, drms)

    ! formal arguments
    type(s_option), target, intent(in)    :: option
    real(wp),       intent(in)    :: coord(:,:)
    real(wp),       intent(in)    :: pbc_box(3,3)
    real(wp),       intent(out)   :: drms

    ! local variables
    integer                       :: i, i1, i2,  num_contact
    real(wp)                      :: d12(1:3), r12
    real(wp)                      :: t, dn, tmp
    real(wp)                      :: box(1:3)
    integer,              pointer :: contact_list(:,:)
    real(wp),             pointer :: contact_dist(:)
    logical                       :: pbc_flag

    contact_list   => option%contact_list
    contact_dist   => option%contact_dist

    pbc_flag = option%pbc_correct

    box(1) = pbc_box(1,1)
    box(2) = pbc_box(2,2)
    box(3) = pbc_box(3,3)
    if (pbc_flag) then
      if (box(1) <= 0.0_wp .or.  &
        box(2) <= 0.0_wp .or.    &
        box(3) <= 0.0_wp  ) then
        call error_msg('Compute_Drms> Box is required for pbc_correct')
      endif
    endif

    num_contact = option%num_contact
    dn          = real(num_contact,wp)
    drms = 0.0_wp

    !$omp parallel do default(none)                          &
    !$omp private(i,  d12, r12,  tmp, i1, i2)    &
    !$omp shared(coord, contact_list, contact_dist, num_contact, box, pbc_flag)   &
    !$omp reduction(+:drms)
    !

    do i = 1, num_contact

      i1 = contact_list(1,i)
      i2 = contact_list(2,i)
      d12(1:3) = coord(1:3,i1) - coord(1:3,i2)
      if (pbc_flag) &
        d12(1:3) = d12(1:3)-anint(d12(1:3)/box(1:3))*box(1:3)
      r12      = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )
        
      tmp      = r12 - contact_dist(i)
      drms     = drms + tmp*tmp

    end do
    !$omp end parallel do

    if (dn > EPS) drms = sqrt(drms/dn)

    return
  end subroutine compute_drms


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine  compute_drms_two_states
  !> @brief      calculate drms
  !! @authors    CK
  !! @param[in]  enefunc potential energy functions [str]
  !! @param[in]  coord coordinates of target systems [dble]
  !! @date       2018/08/16 (CK)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_drms_two_states(option, coord, pbc_box, drms, drms_cur)

    ! formal arguments
    type(s_option), target, intent(in)    :: option
    real(wp),       intent(in)    :: coord(:,:)
    real(wp),       intent(in)    :: pbc_box(3,3)
    real(wp),       intent(out)   :: drms
    real(wp),       intent(out)   :: drms_cur

    ! local variables
    integer                       :: i, i1, i2,  num_contact
    real(wp)                      :: d12(1:3), r12
    real(wp)                      :: t, dn, tmp, tmp2
    real(wp)                      :: box(1:3)
    integer,              pointer :: contact_list(:,:)
    real(wp),             pointer :: contact_dist(:)
    real(wp),             pointer :: contact_cur_dist(:)
    logical                       :: pbc_flag

    contact_list     => option%contact_list
    contact_dist     => option%contact_dist
    contact_cur_dist => option%contact_cur_dist

    pbc_flag = option%pbc_correct

    box(1) = pbc_box(1,1)
    box(2) = pbc_box(2,2)
    box(3) = pbc_box(3,3)
    if (pbc_flag) then
      if (box(1) <= 0.0_wp .or.  &
        box(2) <= 0.0_wp .or.    &
        box(3) <= 0.0_wp  ) then
        call error_msg('Compute_Drms> Box is required for pbc_correct')
      endif
    endif

    num_contact = option%num_contact
    dn          = real(num_contact,wp)
    drms        = 0.0_wp
    drms_cur    = 0.0_wp

    !$omp parallel do default(none)                          &
    !$omp private(i,  d12, r12,  tmp, tmp2, i1, i2)    &
    !$omp shared(coord, contact_list, contact_dist, contact_cur_dist, num_contact,  &
    !$omp        box, pbc_flag)   &
    !$omp reduction(+:drms) reduction(+:drms_cur)
    !

    do i = 1, num_contact

      i1 = contact_list(1,i)
      i2 = contact_list(2,i)
      d12(1:3) = coord(1:3,i1) - coord(1:3,i2)
      if (pbc_flag) &
        d12(1:3) = d12(1:3)-anint(d12(1:3)/box(1:3))*box(1:3)
      r12      = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )
        
      tmp      = r12 - contact_dist(i)
      tmp2     = r12 - contact_cur_dist(i)
      drms     = drms     + tmp*tmp
      drms_cur = drms_cur + tmp2*tmp2

    end do
    !$omp end parallel do

    if (dn > EPS) drms     = sqrt(drms/dn)
    if (dn > EPS) drms_cur = sqrt(drms_cur/dn)

    return
  end subroutine compute_drms_two_states

end module dr_analyze_mod
