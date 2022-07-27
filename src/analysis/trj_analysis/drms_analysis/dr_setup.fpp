!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   dr_setup_mod
!> @brief   setup variables and structures in DRMS_ANALYSIS
!! @authors Chigusa Kobayashi (CK)
!
!  (c) Copyright 2018 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module dr_setup_mod

  use dr_control_mod
  use dr_option_mod
  use dr_option_str_mod
  use trajectory_mod
  use output_mod
  use input_mod
  use trajectory_str_mod
  use output_str_mod
  use input_str_mod
  use select_mod
  use molecules_mod
  use molecules_str_mod
  use fileio_prmtop_mod
  use fileio_ambcrd_mod
  use fileio_grotop_mod
  use fileio_grocrd_mod
  use fileio_psf_mod
  use fileio_pdb_mod
  use constants_mod
 
  implicit none
  private

  ! subroutines
  public  :: setup
  private :: setup_contact_list
  private :: setup_contact_list_two_states
  private :: calc_avoid_bonding

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup
  !> @brief        setup variables and structures in DRMS_ANALYSIS
  !! @authors      NT
  !! @param[in]    ctrl_data  : information of control parameters
  !! @param[inout] molecule   : molecule information
  !! @param[inout] trj_list   : trajectory file list information
  !! @param[inout] option     : option information
  !! @param[inout] output     : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup(ctrl_data, molecule, trj_list, trajectory, output, option)

    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_molecule),        intent(inout) :: molecule
    type(s_trj_list),        intent(inout) :: trj_list
    type(s_trajectory),      intent(inout) :: trajectory
    type(s_output),          intent(inout) :: output
    type(s_option),          intent(inout) :: option

    ! local variables
    type(s_psf)              :: psf
    type(s_pdb)              :: pdb
    type(s_pdb)              :: ref
    type(s_prmtop)           :: prmtop
    type(s_ambcrd)           :: ambcrd
    type(s_ambcrd)           :: ambref
    type(s_grotop)           :: grotop
    type(s_grocrd)           :: grocrd
    type(s_grocrd)           :: groref


    ! input files
    !
    call input_files(ctrl_data%inp_info, &
                     psf=psf,            &
                     ref=ref,            &
                     pdb=pdb,            &
                     prmtop=prmtop,      &
                     ambcrd=ambcrd,      &
                     ambref=ambref,      &
                     grotop=grotop,      &
                     grocrd=grocrd,      &
                     groref=groref)


    ! define molecules
    !
    call define_molecules(molecule, pdb=pdb, ref=ref,       &
                                    psf=psf,       &
                                    prmtop=prmtop, &
                                    ambcrd=ambcrd, &
                                    ambref=ambref, &
                                    grotop=grotop, &
                                    grocrd=grocrd, &
                                    groref=groref)

    call dealloc_pdb_all(ref)
    call dealloc_pdb_all(pdb)
    call dealloc_psf_all(psf)
    call dealloc_prmtop_all(prmtop)
    call dealloc_ambcrd_all(ambref)
    call dealloc_ambcrd_all(ambcrd)
    call dealloc_grotop_all(grotop)
    call dealloc_grocrd_all(grocrd)
    call dealloc_grocrd_all(groref)


    ! setup trajectory
    !
    call setup_trajectory(ctrl_data%trj_info, molecule, trj_list, trajectory)

    ! setup output 
    !
    call setup_output(ctrl_data%out_info, output)

    ! setup selection
    !
    call setup_selection(ctrl_data%sel_info, molecule)

    ! setup option
    !
    call setup_option(ctrl_data%opt_info, ctrl_data%sel_info, &
                      molecule, option)

    ! setup contact list
    !
    if (option%two_states) then
      call setup_contact_list_two_states(molecule, option)
    else
      call setup_contact_list(molecule, option)
    endif

    return

  end subroutine setup

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_contact_list
  !> @brief        setup contact list
  !! @authors      CK
  !! @param[in]    molecule : molecule information
  !! @param[inout] option   : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine setup_contact_list(molecule, option)

    ! formal arguments
    type(s_molecule), intent(in)    :: molecule
    type(s_option),   intent(inout) :: option

    logical                         :: duplicate
    integer                         :: i, j, k
    integer                         :: imax, imin
    integer                         :: jmax, jmin
    integer                         :: iatm, jatm
    integer                         :: ires, jres
    integer                         :: iexcl, jexcl
    integer                         :: ncount, ncount_d
    real(wp)                        :: di(1:3)
    real(wp)                        :: dj(1:3)
    real(wp)                        :: d(1:3), dr
    real(wp)                        :: box_ref(1:3)

    character(20)                   :: atom_id
    character(4)                    :: seg_i, seg_j
    integer                         :: imol, jmol
    integer                         :: ifst, icount, itmp
    integer                         :: igrp1, igrp2

    integer                         :: alloc_stat, dealloc_stat
    logical                         :: segflag, molflag, contflag
    logical                         :: exclflag

    integer, allocatable            :: temp_contact_list(:,:)
    integer, allocatable            :: temp_conv(:)
    real(wp), allocatable           :: temp_contact_dist(:)

    molflag=(molecule%molecule_no(1) > 0)
    segflag=(len_trim(molecule%segment_name(1)) > 0)

    box_ref(1:3) = option%box_size_ref(1:3)

    if (option%identical_group) then
      igrp1 = 1
      igrp2 = 1
    else
      igrp1 = 1
      igrp2 = 2
    endif
    dealloc_stat = 0
    alloc_stat = 0

    do icount = 1, 2
      ncount = 0
      do i = 1, option%num_atoms_group(igrp1)
        iatm = option%contact_atoms(i,igrp1)
        ires = molecule%residue_no(iatm)
        seg_i = molecule%segment_name(iatm)
        imol  = molecule%molecule_no(iatm)
        iexcl = option%exclude_group_list(iatm)
        di(1:3) = molecule%atom_refcoord(1:3,iatm)
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
          d(1:3) = di(1:3) - dj(1:3)
          if (option%pbc_correct_setup) &
            d(1:3) = d(1:3)-anint(d(1:3)/box_ref(1:3))*box_ref(1:3)
          dr = sqrt(d(1)*d(1)+d(2)*d(2)+d(3)*d(3))
          if (dr >= option%minimum_distance .and.  &
              dr < option%maximum_distance) then
            ncount = ncount + 1
            if (icount == 2) then
              temp_contact_dist(ncount)   = dr
              temp_contact_list(1,ncount) = min(iatm,jatm)
              temp_contact_list(2,ncount) = max(iatm,jatm)
            endif
          endif
        end do
      end do

      if (icount == 1) then
        if (ncount > 0) then
          allocate(temp_contact_dist(1:ncount),      &
                   temp_contact_list(1:2, 1:ncount), &
                   temp_conv(1:ncount),              &
                   stat = alloc_stat)
          if (alloc_stat /= 0)   call error_msg_alloc
          call alloc_option(option, DA_Contact, ncount)
        else
          call error_msg('Setup_Contact_List> ERROR : no contact is defined.')
        endif
      endif
    end do

    if (option%avoid_bonding ) then
      if (molecule%num_bonds > 0) then
        call calc_avoid_bonding(molecule, ncount, temp_contact_list, ncount_d, temp_conv)
        option%num_contact = ncount_d
        do i = 1, option%num_contact
          itmp = temp_conv(i)
          option%contact_list(1:2,i) = temp_contact_list(1:2,itmp)
          option%contact_dist(i) = temp_contact_dist(itmp)
        end do
      else
        call error_msg('Setup_Contact_List> ERROR : bond/angle/dihedral information is required in avoid_bonding option')
      endif
    else
      option%num_contact = ncount
      option%contact_list(1:2,1:ncount) = temp_contact_list(1:2,1:ncount)
      option%contact_dist(1:ncount)     = temp_contact_dist(1:ncount)
    endif

    write(MsgOut,'(A)') 'Setup_Contact_List> Contact List'
    write(MsgOut,'(A20,I10)') '  # of contacts   = ', option%num_contact
    do i = 1, option%num_contact
      write(MsgOut,'(A8,I0,A2,$)') ' Contact',i,': '

      do k  = 1, 2
        iatm = option%contact_list(k,i)
        if (molecule%segment_name(iatm) .ne. "") then
          write(atom_id,'(A,A,A,A,I0,A,A,A,A)') trim(molecule%segment_name(iatm)),':', &
                                            trim(molecule%residue_name(iatm)),':', &
                                            molecule%residue_no(iatm),':', &
                                            trim(molecule%atom_name(iatm))
        else
          write(atom_id,'(A,A,I0,A,A)') trim(molecule%residue_name(iatm)),':', &
                                            molecule%residue_no(iatm),':', &
                                            trim(molecule%atom_name(iatm))
        endif
        write(MsgOut,'(I0,A3,A,A3,$)') option%contact_list(k,i), ' ( ',trim(atom_id),' ) '
        if (k == 1) then
          write(MsgOut,'(A3,$)') ' - '
        endif

      end do
      write(MsgOut,'(F10.2)') option%contact_dist(i)

    end do

     deallocate(temp_contact_dist, &
                temp_contact_list, &
                temp_conv,         &
              stat = dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc    

    return

  end subroutine setup_contact_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_contact_list_two_states
  !> @brief        setup contact list considering two states
  !! @authors      CK
  !! @param[in]    molecule : molecule information
  !! @param[inout] option   : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine setup_contact_list_two_states(molecule, option)

    ! formal arguments
    type(s_molecule), intent(in)    :: molecule
    type(s_option),   intent(inout) :: option

    logical                         :: duplicate
    integer                         :: i, j, k
    integer                         :: imax, imin
    integer                         :: jmax, jmin
    integer                         :: iatm, jatm
    integer                         :: ires, jres
    integer                         :: iexcl, jexcl
    integer                         :: ncount, ncount_d
    real(wp)                        :: di(1:3)
    real(wp)                        :: dj(1:3)
    real(wp)                        :: d(1:3), dr
    real(wp)                        :: ei(1:3)
    real(wp)                        :: ej(1:3)
    real(wp)                        :: e(1:3), er
    real(wp)                        :: box_cur(1:3)
    real(wp)                        :: box_ref(1:3)

    character(20)                   :: atom_id
    character(4)                    :: seg_i, seg_j
    integer                         :: imol, jmol
    integer                         :: ifst, icount, itmp
    integer                         :: igrp1, igrp2

    integer                         :: alloc_stat, dealloc_stat
    logical                         :: segflag, molflag, contflag
    logical                         :: exclflag

    integer, allocatable            :: temp_contact_list(:,:)
    integer, allocatable            :: temp_conv(:)
    real(wp), allocatable           :: temp_contact_dist(:)
    real(wp), allocatable           :: temp_contact_dist_other(:)

    molflag=(molecule%molecule_no(1) > 0)
    segflag=(len_trim(molecule%segment_name(1)) > 0)

    box_cur(1:3) = option%box_size_cur(1:3)
    box_ref(1:3) = option%box_size_ref(1:3)

    dealloc_stat = 0
    alloc_stat = 0

    if (option%identical_group) then
      igrp1 = 1
      igrp2 = 1
    else
      igrp1 = 1
      igrp2 = 2
    endif

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
            endif
          endif
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
          call alloc_option(option, DA_Contact, ncount)
        else
          call error_msg('Setup_Contact_List> ERROR : no contact is defined.')
        endif
      endif
    end do

    if (option%avoid_bonding ) then
      if (molecule%num_bonds > 0) then
        call calc_avoid_bonding(molecule, ncount, temp_contact_list, ncount_d, temp_conv)
        option%num_contact = ncount_d
        do i = 1, option%num_contact
          itmp = temp_conv(i)
          option%contact_list(1:2,i) = temp_contact_list(1:2,itmp)
          option%contact_dist(i)     = temp_contact_dist(itmp)
          option%contact_cur_dist(i) = temp_contact_dist_other(itmp)
        end do
      else
        call error_msg('Setup_Contact_List> ERROR : bond/angle/dihedral information is required in avoid_bonding option')
      endif
    else
      option%num_contact = ncount
      option%contact_list(1:2,1:ncount) = temp_contact_list(1:2,1:ncount)
      option%contact_dist(1:ncount)     = temp_contact_dist(1:ncount)
      option%contact_cur_dist(1:ncount) = temp_contact_dist_other(1:ncount)
    endif

    write(MsgOut,'(A)') 'Setup_Contact_List> Contact List'
    write(MsgOut,'(A20,I10)') '  # of contacts   = ', option%num_contact
    do i = 1, option%num_contact
      write(MsgOut,'(A8,I0,A2,$)') ' Contact',i,': '

      do k  = 1, 2
        iatm = option%contact_list(k,i)
        if (molecule%segment_name(iatm) .ne. "") then
          write(atom_id,'(A,A,A,A,I0,A,A,A,A)') trim(molecule%segment_name(iatm)),':', &
                                            trim(molecule%residue_name(iatm)),':', &
                                            molecule%residue_no(iatm),':', &
                                            trim(molecule%atom_name(iatm))
        else
          write(atom_id,'(A,A,I0,A,A)') trim(molecule%residue_name(iatm)),':', &
                                            molecule%residue_no(iatm),':', &
                                            trim(molecule%atom_name(iatm))
        endif
        write(MsgOut,'(I0,A3,A,A3,$)') option%contact_list(k,i), ' ( ',trim(atom_id),' ) '
        if (k == 1) then
          write(MsgOut,'(A3,$)') ' - '
        endif

      end do
      write(MsgOut,'(2F10.2)') option%contact_dist(i), option%contact_cur_dist(i)

    end do

     deallocate(temp_contact_dist,       &
                temp_contact_list,       &
                temp_contact_dist_other, &
                temp_conv,               &
              stat = dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc    

    return

  end subroutine setup_contact_list_two_states

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calc_avoid_bonding
  !> @brief        calc_avoid_bonding
  !! @authors      CK
  !! @param[in]    molecule          : molecule information
  !! @param[in]    ncount            : number of members in temporary contact lists
  !! @param[in]    temp_contact_list : temporary contact lists
  !! @param[out]   ncount_d          : number of member in real contact lists
  !! @param[out]   temp_conv         : temporary convert lists
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine calc_avoid_bonding(molecule, ncount, temp_contact_list, ncount_d, temp_conv)

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
        endif
      end do
      if (.not. duplicate) then
        do j = 1, molecule%num_angles
          jmax = max(molecule%angl_list(1,j), molecule%angl_list(3,j))
          jmin = min(molecule%angl_list(1,j), molecule%angl_list(3,j))
          if (imax == jmax .and. imin == jmin) then
            duplicate = .true.
            exit
          endif
        end do
      endif
      if (.not. duplicate) then
        do j = 1, molecule%num_dihedrals
          jmax = max(molecule%dihe_list(1,j), molecule%dihe_list(4,j))
          jmin = min(molecule%dihe_list(1,j), molecule%dihe_list(4,j))
          if (imax == jmax .and. imin == jmin) then
            duplicate = .true.
            exit
          endif
        end do
      endif
      if (.not. duplicate) then
        do j = 1, molecule%num_impropers
          jmax = max(molecule%impr_list(1,j), molecule%impr_list(4,j))
          jmin = min(molecule%impr_list(1,j), molecule%impr_list(4,j))
          if (imax == jmax .and. imin == jmin) then
            duplicate = .true.
            exit
          endif
        end do
      endif
      if (.not. duplicate) then
        ncount_d = ncount_d + 1
        temp_conv(ncount_d) = i
      endif
    end do

    return

  end subroutine calc_avoid_bonding

end module dr_setup_mod
