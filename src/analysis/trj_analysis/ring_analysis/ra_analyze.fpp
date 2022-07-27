!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ra_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ra_analyze_mod

  use ra_option_str_mod
  use fileio_trj_mod
  use fitting_mod
  use fitting_str_mod
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

  integer    ,parameter    :: max_ring=1000
  integer    ,parameter    :: max_selected_ring=10

  ! subroutines
  public  :: analyze
  private :: compute_ring

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        run analyzing trajectories
  !! @authors      CK
  !! @param[in]    trj_list   : trajectory file list information
  !! @param[in]    output     : output information
  !! @param[inout] option     : option information
  !! @param[inout] trajectory : trajectory information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(molecule, output, option)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_output),          intent(in)    :: output
    type(s_option),          intent(inout) :: option

    ! local variables
    integer                  :: i, j , rot_out


    if (option%check_only) &
      return

    ! open output file
    !
    if (output%torfile /= '') &
      call open_file(rot_out, output%torfile, IOFileOutputNew)

    call compute_ring(molecule%atom_coord, molecule, rot_out, option) 


    ! close output file
    !
    if (output%torfile /= '') call close_file(rot_out)

    return

  end subroutine analyze

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine  compute_ring
  !> @brief      calculate quaternion
  !! @authors    CK
  !! @param[in]  enefunc potential energy functions [str]
  !! @param[in]  refcoord coordinates of target systems [dble]
  !! @param[in]  coord coordinates of target systems [dble]
  !! @date       2015/06/10 (CK)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_ring(coord, molecule, rot_out, option)

    ! formal arguments
    real(wp),               intent(in)    :: coord(:,:)
    type(s_molecule),       intent(in)    :: molecule
    integer,                intent(in)    :: rot_out
    type(s_option), target, intent(inout) :: option

    ! parameters

    ! local variables
    integer                  :: i, j, k, iatm, jatm, katm, icoun, kk
    integer                  :: iring, iother, kflag
    integer                  :: ires, iseg, jres, jseg, kres, kseg
    integer                  :: num_ring_penetrations
    real(wp)                 :: ci(3), dc(3)
    real(wp)                 :: dr
    integer                  :: rings_selected(1:max_selected_ring, 1:max_ring)
    integer                  :: num_rings_selected(1:max_ring)


    icoun = 0

    rings_selected(1:max_selected_ring,1:max_ring)=0
    num_rings_selected(1:max_ring)=0

    do i = 1, option%num_ring
      iatm = option%ring_atoms(i)
      ci(1:3) = coord(1:3,iatm)
      iseg = molecule%segment_no(iatm)
      ires = molecule%residue_no(iatm)

      do j = 1, option%num_other
        jatm = option%other_atoms(j)
        jseg = molecule%segment_no(jatm)
        jres = molecule%residue_no(jatm)
        if (iseg /= jseg .or. ires /= jres) then

          dc(1:3) = ci(1:3)-coord(1:3,jatm)
          dr = sqrt(dc(1)*dc(1)+dc(2)*dc(2)+dc(3)*dc(3))
          if (dr < option%min_contact) then
            write(6,*) iatm, jatm, dr
            icoun = icoun + 1
            if (icoun > max_ring) then
              write(0,*) "error, too much ring"
              stop
            end if
            rings_selected(1,icoun) = jatm
            rings_selected(2,icoun) = iatm
            exit
         
          end if
        end if

      end do

    end do

    num_ring_penetrations = icoun

    do i = 1, num_ring_penetrations

      iatm = rings_selected(2,i)
      iseg = molecule%segment_no(iatm)
      ires = molecule%residue_no(iatm)

      jatm = rings_selected(1,i)
      jseg = molecule%segment_no(jatm)
      jres = molecule%residue_no(jatm)

      kflag=0
      kk = 2
      do k = 1, option%num_all_ring
        katm = option%ring_all_atoms(k)
        kseg = molecule%segment_no(katm)
        kres = molecule%residue_no(katm)
        if (kseg .eq. iseg .and. kres == ires) then
          kflag=1
          kk = kk + 1
          if (kk > max_selected_ring) &
              call error_msg('Compute_Ring> too much selected ring atoms')
          rings_selected(kk, i) = katm
        else if (kflag == 1) then
          exit
        end if
      end do
      num_rings_selected(i)=kk

    end do

    write(rot_out,'(A,I0)') '     # of candidates = ', num_ring_penetrations

    do i = 1, num_ring_penetrations
      write(rot_out,'(A,I0)') "  Penetration",i

      iatm = rings_selected(1,i)

      ci(1:3)=coord(1:3,iatm)
      write(rot_out,'(A,x,I10,x,A,x,I4,x,A,x,A)')  &
          '  Atom',                              &
          iatm,                                  & 
          molecule%segment_name(iatm),           &
          molecule%residue_no(iatm),             &
          molecule%residue_name(iatm),           &
          molecule%atom_name(iatm)

      do j = 2, num_rings_selected(i)
        iatm = rings_selected(j,i)
        dc(1:3) = ci(1:3)-coord(1:3,iatm)
        dr = sqrt(dc(1)*dc(1)+dc(2)*dc(2)+dc(3)*dc(3))

        write(rot_out,'(A,x,I10,x,A,x,I4,x,A,x,A,x,F7.4)')     &
            '  Ring',                            &
            iatm,                                &  
            molecule%segment_name(iatm),         &
            molecule%residue_no(iatm),           &
            molecule%residue_name(iatm),         &
            molecule%atom_name(iatm), dr
      end do
    end do

    return

  end subroutine compute_ring

end module ra_analyze_mod
