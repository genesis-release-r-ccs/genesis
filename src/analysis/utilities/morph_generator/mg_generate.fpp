!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   mg_generate_mod
!> @brief   generate morph files for GENESIS
!! @authors Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module mg_generate_mod

  use mg_option_mod
  use mg_option_str_mod
  use output_str_mod
  use input_str_mod
  use molecules_str_mod
  use fileio_morph_mod
  use fileio_grotop_mod
  use select_atoms_str_mod
  use fileio_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! subroutine
  public  :: generate
  private :: generate_morph_pairs
  private :: calc_avoid_bonding !ported from drms_analysis

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    generate
  !> @brief        generate molecule files
  !! @authors      NT
  !! @param[in]    input    : input information
  !! @param[in]    option   : option information
  !! @param[in]    output   : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine generate(molecule, output, option)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_output),          intent(in)    :: output
    type(s_option),          intent(inout) :: option

    ! local variables
    integer                  :: dealloc_stat


    if (option%check_only) &
      return

    call generate_morph_pairs(molecule, option)

    call output_morph_in(trim(output%morphfile),option%morph_in) 


    ! deallocate memory
    !
    call dealloc_morph_in(option%morph_in)

    return

  end subroutine generate

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    generate_morph_pairs
  !> @brief        generate AAGO-model GROTOP file
  !! @authors      CK
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine generate_morph_pairs(molecule, option)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_option), target,  intent(inout) :: option

    ! local variables
    integer                  :: i, j, k, natom
    integer                  :: iatm, jatm
    integer                  :: imax, imin
    integer                  :: jmax, jmin
    integer                  :: ires, jres
    integer                  :: nrexcl
    integer                  :: iexcl, jexcl
    integer                  :: ncount, ncount_d
    real(wp)                 :: di(1:3)
    real(wp)                 :: dj(1:3)
    real(wp)                 :: d(1:3), dr
    real(wp)                 :: ei(1:3)
    real(wp)                 :: ej(1:3)
    real(wp)                 :: e(1:3), er
    real(wp)                 :: box_cur(1:3)
    real(wp)                 :: box_ref(1:3)

    character(20)            :: atom_id
    character(4)             :: seg_i, seg_j
    integer                  :: imol, jmol
    integer                  :: ifst, icount, itmp
    integer                  :: igrp1, igrp2

    integer                  :: alloc_stat, dealloc_stat
    logical                  :: segflag, molflag, contflag
    logical                  :: exclflag

    integer,  allocatable    :: temp_contact_list(:,:)
    integer,  allocatable    :: temp_conv(:)
    real(wp), allocatable    :: temp_contact_dist(:)
    real(wp), allocatable    :: temp_contact_dist_other(:)


    molflag = (molecule%molecule_no(1) > 0)
    segflag = (len_trim(molecule%segment_name(1)) > 0)

    box_cur(1:3) = option%box_size_cur(1:3)
    box_ref(1:3) = option%box_size_ref(1:3)

    natom  = molecule%num_atoms

    write(MsgOut,'(A)') 'Generate_Morph_Pairs> '

    dealloc_stat = 0
    alloc_stat = 0

    if (option%identical_group) then
      igrp1 = 1
      igrp2 = 1
    else
      igrp1 = 1
      igrp2 = 2
    end if

    do icount = 1, 2
      ncount = 0
      do i = 1, option%num_atoms_group(igrp1)
        iatm = option%contact_atoms(i,igrp1)
        ires = molecule%residue_no(iatm)
        seg_i = molecule%segment_name(iatm)
        imol  = molecule%molecule_no(iatm)
        iexcl = option%exclude_group_list(iatm)
        di(1:3) = molecule%atom_refcoord(1:3,iatm)
        ei(1:3) = molecule%atom_coord(1:3,iatm)
        ifst  = 1
        if (option%identical_group) ifst = i+1
        do j = ifst,option%num_atoms_group(igrp2)
          jatm = option%contact_atoms(j,igrp2)
          jres = molecule%residue_no(jatm)
          seg_j = molecule%segment_name(jatm)
          jmol  = molecule%molecule_no(jatm)
          jexcl = option%exclude_group_list(jatm)
          exclflag= (iexcl > 0 .and. iexcl .eq. jexcl) 
          if (exclflag) cycle
          if (iatm == jatm) cycle
     
          contflag=.true.
          if (segflag) contflag=(seg_i .eq. seg_j)
          if (molflag) contflag=(imol .eq. jmol)
          if (contflag .and. abs(ires-jres) < option%exclude_residues) cycle
          dj(1:3) = molecule%atom_refcoord(1:3,jatm)
          ej(1:3) = molecule%atom_coord(1:3,jatm)
          d(1:3) = di(1:3) - dj(1:3)
          e(1:3) = ei(1:3) - ej(1:3)
          if (option%pbc_correct_setup) &
            d(1:3) = d(1:3)-anint(d(1:3)/box_ref(1:3))*box_ref(1:3)
          if (option%pbc_correct_setup) &
            e(1:3) = e(1:3)-anint(e(1:3)/box_cur(1:3))*box_cur(1:3)
          dr = sqrt(d(1)*d(1)+d(2)*d(2)+d(3)*d(3))
          er = sqrt(e(1)*e(1)+e(2)*e(2)+e(3)*e(3))
          if (((dr >= option%minimum_distance .and.  &
              dr < option%maximum_distance) .or.     &
               (er >= option%minimum_distance .and.  &
                er < option%maximum_distance)) .and. &
              abs(dr-er) >= option%minimum_difference) then
            ncount = ncount + 1
            if (icount == 2) then
              temp_conv(ncount)               = ncount
              temp_contact_dist(ncount)       = dr
              temp_contact_dist_other(ncount) = er
              temp_contact_list(1,ncount)     = min(iatm,jatm)
              temp_contact_list(2,ncount)     = max(iatm,jatm)
            end if
          end if
        end do
      end do

      if (icount == 1) then
        if (ncount > 0) then
          allocate(temp_contact_dist(1:ncount), &
                   temp_contact_dist_other(1:ncount), &
                   temp_contact_list(1:2, 1:ncount),  &
                   temp_conv(1:ncount),               &
                   stat = alloc_stat)
          if (alloc_stat /= 0)   call error_msg_alloc
        else
          call error_msg('Setup_Contact_List> ERROR : no contact is defined.')
        end if
      end if
    end do

    if (option%avoid_bonding ) then
      if (molecule%num_bonds > 0) then
        call calc_avoid_bonding(molecule, ncount, temp_contact_list, ncount_d, temp_conv)
        option%num_contact = ncount_d
      else
        call error_msg('Setup_Contact_List> ERROR : bond/angle/dihedral information is required in avoid_bonding option')
      end if
    else
      option%num_contact = ncount
      do i = 1, option%num_contact
        temp_conv(i) = i
      end do
    end if
    call alloc_morph_in(option%morph_in, option%num_contact, 0)

    do i = 1, option%num_contact
      itmp = temp_conv(i)
      option%morph_in%morph_bb(i)%atom_idx1  = temp_contact_list(1,itmp)
      option%morph_in%morph_bb(i)%atom_idx2  = temp_contact_list(2,itmp)
      option%morph_in%morph_bb(i)%rmin       = temp_contact_dist(itmp)
      option%morph_in%morph_bb(i)%rmin_other = temp_contact_dist_other(itmp)
    end do

    write(MsgOut,'(A)') 'Setup_Contact_List> Contact List'
    write(MsgOut,'(A20,I10)') '  # of contacts   = ', option%num_contact

    do i = 1, option%num_contact
      write(MsgOut,'(A8,I0,A2,$)') ' Contact',i,': '

      iatm = option%morph_in%morph_bb(i)%atom_idx1
      if (molecule%segment_name(iatm) .ne. "") then
        write(atom_id,'(A,A,A,A,I0,A,A,A,A)') trim(molecule%segment_name(iatm)),':', &
                                          trim(molecule%residue_name(iatm)),':', &
                                          molecule%residue_no(iatm),':', &
                                          trim(molecule%atom_name(iatm))
      else
        write(atom_id,'(A,A,I0,A,A)') trim(molecule%residue_name(iatm)),':', &
                                          molecule%residue_no(iatm),':', &
                                          trim(molecule%atom_name(iatm))
      end if
      write(MsgOut,'(I0,A3,A,A3,$)') iatm, ' ( ',trim(atom_id),' ) '
      write(MsgOut,'(A3,$)') ' - '

      iatm = option%morph_in%morph_bb(i)%atom_idx2
      if (molecule%segment_name(iatm) .ne. "") then
        write(atom_id,'(A,A,A,A,I0,A,A,A,A)') trim(molecule%segment_name(iatm)),':', &
                                          trim(molecule%residue_name(iatm)),':', &
                                          molecule%residue_no(iatm),':', &
                                          trim(molecule%atom_name(iatm))
      else
        write(atom_id,'(A,A,I0,A,A)') trim(molecule%residue_name(iatm)),':', &
                                          molecule%residue_no(iatm),':', &
                                          trim(molecule%atom_name(iatm))
      end if
      write(MsgOut,'(I0,A3,A,A3,$)') iatm, ' ( ',trim(atom_id),' ) '


      write(MsgOut,'(2F10.2)') option%morph_in%morph_bb(i)%rmin, &
                               option%morph_in%morph_bb(i)%rmin_other

    end do

     deallocate(temp_contact_dist,       &
                temp_contact_list,       &
                temp_contact_dist_other, &
                temp_conv,               &
              stat = dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc    

    return

  end subroutine generate_morph_pairs

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calc_avoid_bonding
  !> @brief        calc_avoid_bonding
  !! @authors      CK
  !! @param[in]    molecule   : molecule information
  !! @param[in]    ncount     : number of members in temporary contact lists
  !! @param[in]    temp_contact_list : temporary contact lists
  !! @param[out]   ncount_d   : number of member in real contact lists
  !! @param[out]   temp_conv  : temporary convert lists
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine calc_avoid_bonding(molecule, ncount, temp_contact_list, &
                                ncount_d, temp_conv)

    ! formal arguments
    type(s_molecule), intent(in)    :: molecule
    integer,          intent(in)    :: ncount
    integer,          intent(in)    :: temp_contact_list(:,:)
    integer,          intent(out)   :: ncount_d
    integer,          intent(out)   :: temp_conv(:)

    logical                         :: duplicate
    integer                         :: i, j, k
    integer                         :: imax, imin
    integer                         :: jmax, jmin
    integer                         :: iatm, jatm

!
! skip contacts with bonding terms
!
    ncount_d = 0
    do i = 1, ncount
      imin = temp_contact_list(1,i)
      imax = temp_contact_list(2,i)
      duplicate = .false.

      do j = 1, molecule%num_bonds
        jmax = max(molecule%bond_list(1,j), molecule%bond_list(2,j))
        jmin = min(molecule%bond_list(1,j), molecule%bond_list(2,j))
        if (imax == jmax .and. imin == jmin) then
          duplicate = .true.
          exit
        end if
      end do
      if (.not. duplicate) then
        do j = 1, molecule%num_angles
          jmax = max(molecule%angl_list(1,j), molecule%angl_list(3,j))
          jmin = min(molecule%angl_list(1,j), molecule%angl_list(3,j))
          if (imax == jmax .and. imin == jmin) then
            duplicate = .true.
            exit
          end if
        end do
      end if
      if (.not. duplicate) then
        do j = 1, molecule%num_dihedrals
          jmax = max(molecule%dihe_list(1,j), molecule%dihe_list(4,j))
          jmin = min(molecule%dihe_list(1,j), molecule%dihe_list(4,j))
          if (imax == jmax .and. imin == jmin) then
            duplicate = .true.
            exit
          end if
        end do
      end if
      if (.not. duplicate) then
        do j = 1, molecule%num_impropers
          jmax = max(molecule%impr_list(1,j), molecule%impr_list(4,j))
          jmin = min(molecule%impr_list(1,j), molecule%impr_list(4,j))
          if (imax == jmax .and. imin == jmin) then
            duplicate = .true.
            exit
          end if
        end do
      end if
      if (.not. duplicate) then
        ncount_d = ncount_d + 1
        temp_conv(ncount_d) = i
      end if
    end do

    return

  end subroutine calc_avoid_bonding

end module mg_generate_mod
