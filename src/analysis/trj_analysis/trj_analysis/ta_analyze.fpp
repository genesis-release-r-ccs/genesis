!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ta_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Norio Takase (NT), Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ta_analyze_mod

  use ta_option_str_mod
  use fileio_trj_mod
  use measure_mod
  use trajectory_str_mod
  use output_str_mod
  use molecules_str_mod
  use fileio_mod
  use messages_mod
  use constants_mod
 
  implicit none
  private

  ! subroutines
  public  :: analyze
  private :: analyze_dis
  private :: analyze_ang
  private :: analyze_tor
  private :: analyze_comdis
  private :: analyze_comang
  private :: analyze_comtor
  private :: out_result
  private :: print_output_info

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        run analyzing trajectories
  !! @authors      NT, TM
  !! @param[in]    trj_list   : trajectory file list information
  !! @param[in]    output     : output information
  !! @param[inout] option     : option information
  !! @param[inout] trajectory : trajectory information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(molecule, trj_list, output, option, trajectory)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_trj_list),        intent(in)    :: trj_list
    type(s_output),          intent(in)    :: output
    type(s_option),          intent(inout) :: option
    type(s_trajectory),      intent(inout) :: trajectory

    ! local variables
    type(s_trj_file)         :: trj_in
    integer                  :: nstru, ifile, istep, num_trjfiles
    integer                  :: dis_unit, ang_unit, tor_unit
    integer                  :: cdis_unit, cang_unit, ctor_unit


    if (option%check_only) &
      return

    ! open output file
    !
    if (option%out_dis)  &
      call open_file(dis_unit, output%disfile, IOFileOutputNew)
    if (option%out_ang)  &
      call open_file(ang_unit, output%angfile, IOFileOutputNew)
    if (option%out_tor)  &
      call open_file(tor_unit, output%torfile, IOFileOutputNew)
    if (option%out_cdis) &
      call open_file(cdis_unit, output%comdisfile, IOFileOutputNew)
    if (option%out_cang) &
      call open_file(cang_unit, output%comangfile, IOFileOutputNew)
    if (option%out_ctor) &
      call open_file(ctor_unit, output%comtorfile, IOFileOutputNew)


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

          if (option%out_dis) then
            call analyze_dis(trajectory, option)
            call out_result (nstru, dis_unit, option%distance)
          end if

          if (option%out_ang) then
            call analyze_ang(trajectory, option)
            call out_result (nstru, ang_unit, option%angle)
          end if

          if (option%out_tor) then
            call analyze_tor(trajectory, option)
            call out_result (nstru, tor_unit, option%torsion)
          end if

          if (option%out_cdis) then
            call analyze_comdis(molecule, trajectory, option)
            call out_result (nstru, cdis_unit, option%cdistance)
          end if

          if (option%out_cang) then
            call analyze_comang(molecule, trajectory, option)
            call out_result (nstru, cang_unit, option%cangle)
          end if

          if (option%out_ctor) then
            call analyze_comtor(molecule, trajectory, option)
            call out_result (nstru, ctor_unit, option%ctorsion)
          end if

        end if

      end do

      ! close trajectory file
      !
      call close_trj(trj_in)

    end do


    ! close output file
    !
    if (option%out_ctor) call close_file(ctor_unit)
    if (option%out_cang) call close_file(cang_unit)
    if (option%out_cdis) call close_file(cdis_unit)
    if (option%out_tor)  call close_file(tor_unit)
    if (option%out_ang)  call close_file(ang_unit)
    if (option%out_dis)  call close_file(dis_unit)


    ! Output summary
    !
    call print_output_info(output, option)

    return

  end subroutine analyze

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze_dis
  !> @brief        analyze distances
  !! @authors      NT, TM, SI
  !! @param[in]    trajectory : trajectory information
  !! @param[inout] option     : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze_dis(trajectory, option)

    ! formal argments
    type(s_trajectory),      intent(in)    :: trajectory
    type(s_option),          intent(inout) :: option

    ! local variables
    integer                  :: i, j, idx1, idx2


    do i = 1, size(option%distance)
      option%distance(i) = 0.0_wp
      
      do j = 1, option%dist_num(i)/2
        idx1 = option%dist_list(2*j-1,i)
        idx2 = option%dist_list(2*j,i)

        option%distance(i) = option%distance(i) + option%dist_weight(i,j) &
                           * compute_dis(trajectory%coord(:,idx1), &
                                          trajectory%coord(:,idx2))
      end do

    end do

    return

  end subroutine analyze_dis

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze_ang
  !> @brief        analize angles
  !! @authors      NT, TM
  !! @param[in]    trajectory : trajectory information
  !! @param[inout] option     : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze_ang(trajectory, option)

    ! formal argments
    type(s_trajectory),      intent(in)    :: trajectory
    type(s_option),          intent(inout) :: option

    ! local variables
    integer                  :: i, idx1, idx2, idx3


    do i = 1, size(option%angle)
      
      idx1 = option%angl_list(1, i)
      idx2 = option%angl_list(2, i)
      idx3 = option%angl_list(3, i)

      option%angle(i) = compute_ang(trajectory%coord(:,idx1), &
                                    trajectory%coord(:,idx2), &
                                    trajectory%coord(:,idx3))

    end do

    return

  end subroutine analyze_ang

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze_tor
  !> @brief        analize torsions
  !! @authors      NT, TM
  !! @param[in]    trajectory : trajectory informatin
  !! @param[inout] option     : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze_tor(trajectory, option)

    ! formal argments
    type(s_trajectory),      intent(in)    :: trajectory
    type(s_option),          intent(inout) :: option

    ! local variables
    integer                  :: i, idx1, idx2, idx3, idx4


    do i = 1, size(option%torsion)
      
      idx1 = option%tors_list(1, i)
      idx2 = option%tors_list(2, i)
      idx3 = option%tors_list(3, i)
      idx4 = option%tors_list(4, i)

      option%torsion(i) = compute_dih(trajectory%coord(:,idx1), &
                                      trajectory%coord(:,idx2), &
                                      trajectory%coord(:,idx3), &
                                      trajectory%coord(:,idx4))

    end do

    return

  end subroutine analyze_tor

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze_comdis
  !> @brief        analyze COM distances
  !! @authors      NT
  !! @param[in]    molecule   : molecule information
  !! @param[in]    trajectory : trajectory information
  !! @param[inout] option     : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze_comdis(molecule, trajectory, option)

    ! formal argments
    type(s_molecule),        intent(in)    :: molecule
    type(s_trajectory),      intent(in)    :: trajectory
    type(s_option),          intent(inout) :: option

    ! local variables
    real(wp)                 :: c1(3), c2(3)
    integer                  :: i


    do i = 1, size(option%cdist_group(1,:))

      ! atom group 1
      c1 = compute_com(trajectory%coord, molecule%mass, &
                       option%selatoms(option%cdist_group(1,i))%idx)

      ! atom group 2
      c2 = compute_com(trajectory%coord, molecule%mass, &
                       option%selatoms(option%cdist_group(2,i))%idx)

      ! compute distance
      option%cdistance(i) = compute_dis(c1, c2)

    end do

    return

  end subroutine analyze_comdis

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze_comang
  !> @brief        analyze COM angle
  !! @authors      NT
  !! @param[in]    molecule   : molecule information
  !! @param[in]    trajectory : trajectory information
  !! @param[inout] option     : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze_comang(molecule, trajectory, option)

    ! formal argments
    type(s_molecule),        intent(in)    :: molecule
    type(s_trajectory),      intent(in)    :: trajectory
    type(s_option), target,  intent(inout) :: option

    ! local variables
    real(wp)                 :: c1(3), c2(3), c3(3)
    integer                  :: i


    do i = 1, size(option%cangl_group(1,:))

      ! atom group 1
      c1 = compute_com(trajectory%coord, molecule%mass, &
                       option%selatoms(option%cangl_group(1,i))%idx)

      ! atom group 2
      c2 = compute_com(trajectory%coord, molecule%mass, &
                       option%selatoms(option%cangl_group(2,i))%idx)

      ! atom group 3
      c3 = compute_com(trajectory%coord, molecule%mass, &
                       option%selatoms(option%cangl_group(3,i))%idx)

      ! compute angle
      option%cangle(i) = compute_ang(c1, c2, c3)

    end do

    return

  end subroutine analyze_comang

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze_comtor
  !> @brief        analyze COM torsion
  !! @authors      NT
  !! @param[in]    molecule   : molecule information
  !! @param[in]    trajectory : trajectory information
  !! @param[inout] option     : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze_comtor(molecule, trajectory, option)

    ! formal argments
    type(s_molecule),        intent(in)    :: molecule
    type(s_trajectory),      intent(in)    :: trajectory
    type(s_option),          intent(inout) :: option

    ! local variables
    real(wp)                 :: c1(3), c2(3), c3(3), c4(3)
    integer                  :: i


    do i = 1, size(option%ctor_group(1,:))

      ! atom group 1
      c1 = compute_com(trajectory%coord, molecule%mass, &
                       option%selatoms(option%ctor_group(1,i))%idx)

      ! atom group 2
      c2 = compute_com(trajectory%coord, molecule%mass, &
                       option%selatoms(option%ctor_group(2,i))%idx)

      ! atom group 3
      c3 = compute_com(trajectory%coord, molecule%mass, &
                       option%selatoms(option%ctor_group(3,i))%idx)

      ! atom group 4
      c4 = compute_com(trajectory%coord, molecule%mass, &
                       option%selatoms(option%ctor_group(4,i))%idx)

      ! compute angle
      option%ctorsion(i) = compute_dih(c1, c2, c3, c4)

    end do

    return

  end subroutine analyze_comtor

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    out_result
  !> @brief        output analysis results
  !! @authors      NT
  !! @param[in]    sturct_no : sturcture number
  !! @param[in]    unit_no   : file unit number
  !! @param[in]    results   : analysis results
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine out_result(struct_no, unit_no, results)

    ! formal arguments
    integer,                 intent(in)    :: struct_no
    integer,                 intent(in)    :: unit_no
    real(wp),                intent(in)    :: results(:)

    ! local variables
    integer                  :: i


    write(unit_no, '(i10,1x,100(f8.3,1x))') &
         struct_no, (results(i), i=1, size(results))

    return

  end subroutine out_result

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    print_output_info
  !> @brief        print detailed output information
  !! @authors      TM
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine print_output_info(output, option)

    ! formal arguments
    type(s_output),          intent(in) :: output
    type(s_option),          intent(in) :: option

    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze> Detailed information in the output files'
    write(MsgOut,'(A)') ''

    if (option%out_dis)  then
      write(MsgOut,'(A)') '  [disfile] ' // trim(output%disfile)
      write(MsgOut,'(A)') '    Column 1: Snapshot index'
      write(MsgOut,'(A)') '    Column 2: Distance (angstrom)'
      write(MsgOut,'(A)') '    If multiple groups were specified in [OPTION],'
      write(MsgOut,'(A)') '    Column N+1: Distance for the N-th selected group'
      write(MsgOut,'(A)') ''
    end if

    if (option%out_ang)  then
      write(MsgOut,'(A)') '  [angfile] ' // trim(output%angfile)
      write(MsgOut,'(A)') '    Column 1: Snapshot index'
      write(MsgOut,'(A)') '    Column 2: Angle (degree)'
      write(MsgOut,'(A)') '    If multiple groups were specified in [OPTION],'
      write(MsgOut,'(A)') '    Column N+1: Angle for the N-th selected group'
      write(MsgOut,'(A)') ''
    end if

    if (option%out_tor)  then
      write(MsgOut,'(A)') '  [torfile] ' // trim(output%torfile)
      write(MsgOut,'(A)') '    Column 1: Snapshot index'
      write(MsgOut,'(A)') '    Column 2: Torsion angle (degree)'
      write(MsgOut,'(A)') '    If multiple groups were specified in [OPTION],'
      write(MsgOut,'(A)') '    Column N+1: Torsion angle for the N-th selected group'
      write(MsgOut,'(A)') ''
    end if

    if (option%out_cdis)  then
      write(MsgOut,'(A)') '  [comdisfile] ' // trim(output%comdisfile)
      write(MsgOut,'(A)') '    Column 1: Snapshot index'
      write(MsgOut,'(A)') '    Column 2: Distance between the centers of mass (angstrom)'
      write(MsgOut,'(A)') '    If multiple groups were specified in [OPTION],'
      write(MsgOut,'(A)') '    Column N+1: Distance for the N-th selected group'
      write(MsgOut,'(A)') ''
    end if

    if (option%out_cang)  then
      write(MsgOut,'(A)') '  [comangfile] ' // trim(output%comangfile)
      write(MsgOut,'(A)') '    Column 1: Snapshot index'
      write(MsgOut,'(A)') '    Column 2: Angle between the centers of mass (degree)'
      write(MsgOut,'(A)') '    If multiple groups were specified in [OPTION],'
      write(MsgOut,'(A)') '    Column N+1: Angle for the N-th selected group'
      write(MsgOut,'(A)') ''
    end if

    if (option%out_ctor)  then
      write(MsgOut,'(A)') '  [comtorfile] ' // trim(output%comtorfile)
      write(MsgOut,'(A)') '    Column 1: Snapshot index'
      write(MsgOut,'(A)') '    Column 2: Torsion angle between the centers of mass (degree)'
      write(MsgOut,'(A)') '    If multiple groups were specified in [OPTION],'
      write(MsgOut,'(A)') '    Column N+1: Torsion angle for the N-th selected group'
      write(MsgOut,'(A)') ''
    end if

    return

  end subroutine print_output_info

end module ta_analyze_mod
