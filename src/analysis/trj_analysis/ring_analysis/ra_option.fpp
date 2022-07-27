!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ra_option_mod
!> @brief   module for analysis options
!! @authors Chigusa Kobayashi (CK)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ra_option_mod

  use ra_option_str_mod
  use select_mod
  use select_atoms_mod
  use molecules_str_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_opt_info
    logical                         :: check_only     = .true.
    integer                         :: num_ring_type
    real(wp)                        :: min_contact = 2.6_wp
    character(6), allocatable       :: ring_atoms_in(:)
  end type s_opt_info

  integer                 :: num_ring_type = 13

  character(6), parameter :: RingatomType(13) = (/'CN1   ', &
                                                  'CN1T  ', &
                                                  'CN2   ', &
                                                  'CN3   ', &
                                                  'CN3T  ', &
                                                  'CN4   ', &
                                                  'CN5   ', &
                                                  'CN5G  ', &
                                                  'CN6   ', &
                                                  'CA    ', &
                                                  'CAI   ', &
                                                  'CPT   ', &
                                                  'CP2   '/)

  ! subroutines
  public  :: show_ctrl_option
  public  :: read_ctrl_option
  public  :: setup_option
  private :: setup_atoms

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_option
  !> @brief        show control parameters in OPTION section
  !! @authors      CK
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_option

    write(MsgOut,'(A)') '[OPTION]'
    write(MsgOut,'(A)') 'check_only           = YES       # (YES/NO)'
    write(MsgOut,'(A)') 'min_contact          = 2.6       # minimum of contact'
    write(MsgOut,'(A)') ' '

    return

  end subroutine show_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_option
  !> @brief        read control parameters in OPTION section
  !! @authors      CK
  !! @param[in]    handle   : unit number of control file
  !! @param[out]   opt_info : OPTION section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine read_ctrl_option(handle, opt_info)
  
    ! parameters
    character(*),            parameter     :: Section = 'Option'

    ! formal argments
    integer,                 intent(in)    :: handle
    integer                  :: i
    character(MaxLine)       :: value
    character(20)            :: atom_type
    integer                  :: natomring
    type(s_opt_info),        intent(inout) :: opt_info


    ! read parameters
    !
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_logical(handle, Section, &
                              'check_only', opt_info%check_only)

    call read_ctrlfile_real(handle, Section, 'min_contact',       &
                                opt_info%min_contact)

    natomring = 0
    do while (.true.)
      value = ''
      write(atom_type,'(A9,I0)') 'atomtypes', natomring + 1
      call read_ctrlfile_string(handle, Section, atom_type, value)
      if (value == '') &
        exit
      natomring = natomring + 1
    end do

    if (natomring >  0) then
      allocate(opt_info%ring_atoms_in(natomring))
      do i = 1, natomring
        write(atom_type,'(A9,I0)') 'atomtypes', i
        call read_ctrlfile_string(handle, Section, atom_type, &
                                  value)
        if (len_trim(value) > 6) then
          opt_info%ring_atoms_in(i) = value(1:6)
        else
          opt_info%ring_atoms_in(i) = value
        endif
      end do
    else
      natomring=num_ring_type 
      allocate(opt_info%ring_atoms_in(natomring))
      do i = 1, natomring
         opt_info%ring_atoms_in(i) = RingatomType(i)
      end do
    end if
    opt_info%num_ring_type = natomring

    call end_ctrlfile_section(handle)


    ! write parameters to MsgOut
    !
    write(MsgOut,'(A)') 'Read_Ctrl_Option> Parameters of Options'
    if (opt_info%check_only) then
      write(MsgOut,'(A20,A3)') '  check only      = ', 'yes'
    else
      write(MsgOut,'(A20,A2)') '  check only      = ', 'no'
    end if
    write(MsgOut,'(A14,F10.4)') &
         ' min_contact= ',opt_info%min_contact

    write(MsgOut,'(A20)') '  ring_atoms_in     '
    do i = 1, opt_info%num_ring_type
      write(MsgOut,'(I3,2X,A)') i, opt_info%ring_atoms_in(i)
    end do

    write(MsgOut,'(A)') ' '

    return

  end subroutine read_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_option
  !> @brief        setup option information
  !! @authors      CK
  !! @param[in]    opt_info : OPTION section control parameters information
  !! @param[in]    molecule : molecule information
  !! @param[inout] option   : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine setup_option(opt_info, molecule, option)

    ! formal arguments
    type(s_opt_info),           intent(in)    :: opt_info
    type(s_molecule),           intent(in)    :: molecule
    type(s_option),             intent(inout) :: option

    integer                     :: i


    ! check only
    !
    option%check_only = opt_info%check_only

    option%min_contact  = opt_info%min_contact
    option%num_ring_type  = opt_info%num_ring_type
    allocate(option%ring_atoms_in(option%num_ring_type))
    do i = 1, option%num_ring_type
      option%ring_atoms_in(i)=opt_info%ring_atoms_in(i)
    end do

    allocate(option%atomflag(molecule%num_atoms))
    option%atomflag(1:molecule%num_atoms)=0

    call setup_atoms(molecule, option)

    return

  end subroutine setup_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_atoms
  !> @brief        setup option information
  !! @authors      CK
  !! @param[in]    molecule : molecule information
  !! @param[inout] option   : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine setup_atoms(molecule, option)

    ! formal arguments
    type(s_molecule),  intent(in)    :: molecule
    type(s_option),    intent(inout) :: option

    integer                          :: i, j, i1, i2
    integer                          :: num_atoms, num_ring, num_other
    integer                          :: num_all_ring
    integer                          :: oldsegno, oldresno, flag,anum
    character(6)                     :: atomcls, resname
    character(6)                     :: atomcls1, atomcls2
    integer, allocatable             :: bond_flag(:)

    num_atoms=molecule%num_atoms

    allocate(bond_flag(1:num_atoms))
    bond_flag(1:num_atoms) = 0
    do i = 1, molecule%num_bonds
      i1 = molecule%bond_list(1,i)
      atomcls1 = molecule%atom_cls_name(i1)
      atomcls1 = trim(atomcls1)
      i2 = molecule%bond_list(2,i)
      atomcls2 = molecule%atom_cls_name(i2)
      atomcls2 = trim(atomcls2)
      
      if (atomcls1(1:1) .ne. "H" .and. atomcls2 .ne. "H") then
        bond_flag(i1) = 1
        bond_flag(i2) = 1
      end if
    end do

    oldsegno = -1
    oldresno = -1
    num_ring = 0
    num_all_ring = 0
    num_other = 0

    do i = 1, num_atoms

      atomcls = molecule%atom_cls_name(i)
      atomcls = trim(atomcls)
      resname = molecule%residue_name(i)
      resname = trim(resname)
      flag = 0
      
      if (resname .ne. "WAT" .and. resname .ne. "TIP3" .and. &
          resname .ne. "SPC" .and. bond_flag(i) == 1) then

      if (atomcls(1:1) .eq. "C" .or.  &
          atomcls(1:1) .eq. "O" .or.  &
          atomcls(1:1) .eq. "N" .or.  &
          atomcls(1:1) .eq. "P" ) then
        anum=0
        if (atomcls(1:1) .eq. "C") then 
          do j = 1, option%num_ring_type
            if (atomcls .eq. option%ring_atoms_in(j)) then
              anum=j
              exit
            end if
          end do
        end if
        if (anum > 0) then
          if (oldsegno /= molecule%segment_no(i) .or. &
              oldresno /= molecule%residue_no(i)) then
              flag=2
              num_ring = num_ring+1
          else
              flag=3
              num_all_ring = num_all_ring+1
          end if
          oldsegno=molecule%segment_no(i)
          oldresno=molecule%residue_no(i)
        else
          flag=1
          num_other = num_other+1
        end if
      end if
      end if
      option%atomflag(i) = flag
     
    end do

    option%num_ring=num_ring
    option%num_all_ring=num_all_ring
    option%num_other=num_other
    write(MsgOut,'(A)') 'Setup_Atoms> Number of ring and others'
    write(MsgOut,'(A,I0, A, I0)') 'rings = ',num_ring,  &
                                  ' all_other_ring = ',num_all_ring, &
                                  ' other = ',num_other

    allocate(option%ring_atoms(1:num_ring), &
             option%ring_all_atoms(1:num_all_ring), &
             option%other_atoms(1:num_other))

    num_ring = 0
    num_all_ring = 0
    num_other = 0

    do i = 1, num_atoms

      if (option%atomflag(i) == 2) then
        num_ring=num_ring+1
        if (num_ring > option%num_ring) then
          write(0,*) 'error num_ring > num_ring0'
          stop
        end if
        option%ring_atoms(num_ring)=i

      else if (option%atomflag(i) == 1) then
        num_other=num_other+1
        if (num_other > option%num_other) then
          write(0,*) 'error num_other > num_other0'
          stop
        end if
        option%other_atoms(num_other)=i

      else if (option%atomflag(i) == 3) then
        num_all_ring=num_all_ring+1
        if (num_all_ring > option%num_all_ring) then
          write(0,*) 'error num_all_ring > num_all_ring0'
          stop
        end if
        option%ring_all_atoms(num_all_ring)=i

      end if
    end do

    deallocate(bond_flag)

    return

  end subroutine setup_atoms

end module ra_option_mod
