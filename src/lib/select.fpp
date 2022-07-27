!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   select_mod
!> @brief   utilities for atom selection
!! @authors Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module select_mod

  use molecules_str_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod

  implicit none
  private

  ! parameters
  integer,    private, parameter :: MaxRangeRes     = 50
  integer,    private, parameter :: MaxGroupRes     = 50
  integer,    private, parameter :: MaxGroupResName = 10

  integer,    private, parameter :: TblType         = 1
  integer,    private, parameter :: TblMolNo        = 2
  integer,    private, parameter :: TblMolID        = 3

  ! structures
  type, public :: s_sel_info
    character(MaxLineLong), allocatable :: groups(:)
    character(MaxLine),     allocatable :: mole_names(:)
  end type s_sel_info

  type s_range_res
    character(20)                :: name
    character(4)                 :: seg_name_b
    character(4)                 :: seg_name_e
    integer                      :: res_no_b
    integer                      :: res_no_e
    character(6)                 :: res_name_b
    character(6)                 :: res_name_e
  end type s_range_res

  type s_group_res
    character(20)                :: name
    integer                      :: num_res
    character(6)                 :: res_name(MaxGroupResName)
  end type s_group_res

  type s_mole_info
    integer                      :: num_range
    integer                      :: num_group
    type(s_range_res)            :: range_res(MaxRangeRes)
    type(s_group_res)            :: group_res(MaxGroupRes)
  end type s_mole_info
  
  ! subroutines
  public  :: show_ctrl_selection
  public  :: read_ctrl_selection
  public  :: setup_selection
  private :: parse_molecule_names
  private :: setup_molecule_no

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_selection
  !> @brief        show SELECTION section usage
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_selection(tool_name)

    ! formal arguments
    character(len=*), optional, intent(in) :: tool_name

    if (.not.present(tool_name)) then
      write(MsgOut,'(A)') '[SELECTION]'
      write(MsgOut,'(A)') '# group1         = all             # selection group 1'
      write(MsgOut,'(A)') '# group2         = molname:protein # selection group 2'
      write(MsgOut,'(A)') '# mole_name1     = protein  P1:1:TYR P1:5:MET'
      write(MsgOut,'(A)') '# mole_name2     = lipid    OLEO:PCGL:OLEO'
      write(MsgOut,'(A)') ' '

    else
      select case (tool_name)

      case ('qg') 
        write(MsgOut,'(A)') '[SELECTION]'
        write(MsgOut,'(A)') '# group1         = atomno:300-330              # selection group 1'
        write(MsgOut,'(A)') '# group2         = resno:501-502               # selection group 2'
        write(MsgOut,'(A)') '# group3         = segid:A                     # selection group 3'
        write(MsgOut,'(A)') '# group4         = (resname:TIP3|segid:ION) &\ # selection group 4'
        write(MsgOut,'(A)') '#                  (resno:502 around_mol:50.0) # selection group 4'
        write(MsgOut,'(A)') '# group5         = segid:A and an:CA           # selection group 5'
        !write(MsgOut,'(A)') '# group6         = molname:protein             # selection group 6'
        !write(MsgOut,'(A)') '# mole_name1     = protein  P1:1:TYR P1:5:MET'
        !write(MsgOut,'(A)') '# mole_name2     = lipid    OLEO:PCGL:OLEO'
        write(MsgOut,'(A)') ' '

      case default
        write(MsgOut,'(A)') '[SELECTION]'
        write(MsgOut,'(A)') '# group1         = all             # selection group 1'
        write(MsgOut,'(A)') '# group2         = molname:protein # selection group 2'
        write(MsgOut,'(A)') '# mole_name1     = protein  P1:1:TYR P1:5:MET'
        write(MsgOut,'(A)') '# mole_name2     = lipid    OLEO:PCGL:OLEO'
        write(MsgOut,'(A)') ' '

      end select
    end if  

    return

  end subroutine show_ctrl_selection
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_selection
  !> @brief        read control parameters in SELECTION section
  !! @authors      NT
  !! @param[in]    handle   : unit number of control file [int] 
  !! @param[out]   sel_info : information of selection section [str] 
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine read_ctrl_selection(handle, sel_info)
  
    ! parameters
    character(*),            parameter     :: Section = 'Selection'

    ! formal argments
    integer,                 intent(in)    :: handle
    type(s_sel_info),        intent(inout) :: sel_info

    ! local variables
    integer                 :: i, ngroup, nmole
    character(MaxLineLong)  :: value
    character(MaxLineLong)  :: mole_name, group_name


    ! read parameters 
    !

    call begin_ctrlfile_section(handle, Section)

    ! check readable selection group
    ngroup = 0
    do while (.true.)
      value = ''
      write(group_name,fmt='(A5,I0)') 'group', ngroup + 1
      call read_ctrlfile_string(handle, Section, group_name, value)
      if (value == '') then
        exit
      end if
      ngroup = ngroup + 1
    end do

    ! check readable molecule name
    nmole = 0
    do while (.true.)
      value = ''
      write(mole_name,fmt='(A9,I0)') 'mole_name', nmole + 1
      call read_ctrlfile_string(handle, Section, mole_name, value)
      if (value == '') then
        exit
      end if
      nmole = nmole + 1
    end do

    ! allocate for selection groups, molecule names
    allocate(sel_info%groups(ngroup), &
             sel_info%mole_names(nmole))
 
    ! read selection groups
    do i = 1, ngroup
      write(group_name,fmt='(A5,I0)') 'group', i
      call read_ctrlfile_string(handle, Section, group_name, &
                                sel_info%groups(i))
    end do

    ! read molecule names
    do i = 1, nmole
      write(mole_name,fmt='(A9,I0)') 'mole_name', i
      call read_ctrlfile_string(handle, Section, mole_name, &
                                sel_info%mole_names(i))
    end do

    call end_ctrlfile_section(handle)


    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Ctrl_Selection> Parameters of Selection'
      write(MsgOut,'(A20,I10)') '  # of groups     = ', ngroup
      do i = 1, ngroup
        write(MsgOut,'(A10,I0,A3,A)') &
             '    group ', i, ' = ', trim(sel_info%groups(i))
      end do
      write(MsgOut,'(A20,I10)') '  # of mole names = ', nmole
      do i = 1, nmole
        write(MsgOut,'(A10,I0,A3,A)') &
             '    name  ', i, ' = ', trim(sel_info%mole_names(i))
      end do
      write(MsgOut,'(A)') ' '
    endif

    return

  end subroutine read_ctrl_selection

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_selection
  !> @brief        setup selection
  !! @authors      NT
  !! @param[in]    sel_info  : information of selection section
  !! @param[inout] molecule  : information of molecule
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine setup_selection(sel_info, molecule)
  
    ! formal argments
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_molecule),        intent(inout) :: molecule

    ! local variables
    type(s_mole_info)        :: mole_info
    integer                  :: nmole


    nmole  = size(sel_info%mole_names)

    ! setup molecule number
    !
    if (nmole > 0) then

      if (allocated(molecule%molecule_atom_no) .and. main_rank) then
        write(MsgOut,'(A)') &
             'Setup_Selection> Original molecule information is overwrote.'
      end if

      call parse_molecule_names(sel_info%mole_names, mole_info)
      call setup_molecule_no   (mole_info, molecule)

      if (molecule%num_molecules > 0 .and. main_rank) then
        write(MsgOut,'(A)') &
             'Setup_Selection> Molecule information is defined newly.'
        write(MsgOut,'(A,I7)') &
             '                   number of molecules : ', &
             molecule%num_molecules
        write(MsgOut,'(A)') ' '
      end if

    end if

    return

  end subroutine setup_selection

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    parse_molecule_names
  !> @brief        parse molecule name definitions
  !! @authors      NT
  !! @param[in]    mole_names : molecule names
  !! @param[inout] mole_info  : molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine parse_molecule_names(mole_names, mole_info)

    ! formal arguments
    character(*),            intent(in)    :: mole_names(:)
    type(s_mole_info),       intent(inout) :: mole_info

    ! local variables
    integer                  :: res_no, i, j, nmole
    character(MaxLine)       :: mole_name
    character(MaxLine)       :: s1, s2, s3
    character(4)             :: seg_name
    character(6)             :: res_name
    

    nmole = size(mole_names)

    mole_info%num_range = 0
    mole_info%num_group = 0

    do j = 1, nmole

      mole_name = mole_names(j)
       
      s1 = ''
      s2 = ''
      s3 = ''
      read(mole_name,*,end=110,err=900) s1, s2, s3

110   if (s1 == '') &
        cycle

      if (s3 /= '') then

        ! residue by range 
        !
        mole_info%num_range = mole_info%num_range + 1
        if (mole_info%num_range > MaxRangeRes) &
          cycle

        do i = 1, len(s2)
          if (s2(i:i) == ':') s2(i:i) = ' '
        end do

        ! begin value
        seg_name = ''
        res_no   = -1
        res_name = ''

        read(s2,*,end=120,err=120) seg_name, res_no, res_name

120     if (seg_name == '') &
          call error_msg('Parse_Molecule_Range> ERROR : molecule name : '//&
          trim(s1))

        if (res_name == '') then

          seg_name = ''
          res_no   = -1

          read(s2,*,end=121,err=121) res_no, res_name

121       if (res_no == -1 .or. res_name == '') &
            call error_msg('Parse_Molecule_Range> ERROR : molecule name : '//&
            trim(s1))

        end if

        mole_info%range_res(mole_info%num_range)%name       = s1
        mole_info%range_res(mole_info%num_range)%seg_name_b = seg_name
        mole_info%range_res(mole_info%num_range)%res_no_b   = res_no
        mole_info%range_res(mole_info%num_range)%res_name_b = res_name

        do i = 1, len(s3)
          if (s3(i:i) == ':') s3(i:i) = ' '
        end do

        ! end value
        seg_name = ''
        res_no   = -1
        res_name = ''

        read(s3,*,end=130,err=130) seg_name, res_no, res_name

130     if (seg_name == '') &
          call error_msg('Parse_Molecule_Range> ERROR : molecule name : '//&
          trim(s1))

        if (res_name == '') then

          seg_name = ''
          res_no   = -1

          read(s3,*,end=131,err=131) res_no, res_name

131       if (res_no == -1 .or. res_name == '') &
            call error_msg('Parse_Molecule_Range> ERROR : molecule name : '//&
            trim(s1))

        end if

        mole_info%range_res(mole_info%num_range)%name       = s1
        mole_info%range_res(mole_info%num_range)%seg_name_e = seg_name
        mole_info%range_res(mole_info%num_range)%res_no_e   = res_no
        mole_info%range_res(mole_info%num_range)%res_name_e = res_name

      else
        ! residue by group
        !

        mole_info%num_group = mole_info%num_group + 1
        if (mole_info%num_group > MaxGroupRes) &
          cycle

        do i = 1, len(s2)
          if (s2(i:i) == ':') s2(i:i) = ' '
        end do

        mole_info%group_res(mole_info%num_group)%name = s1

        read(s2,*,end=140,err=900) &
          (mole_info%group_res(mole_info%num_group)%res_name(i), &
          i=1,MaxGroupResName)

140     continue

        if (i - 1 > 0) then
          mole_info%group_res(mole_info%num_group)%num_res = i - 1
        else
          mole_info%num_group = mole_info%num_group - 1
        end if

      end if

    end do

100 continue

    return

900 call error_msg_fileio
    return

  end subroutine parse_molecule_names

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_molecule_no
  !> @brief        setup molecule number
  !! @authors      NT
  !! @param[in]    mole_info : molecule information
  !! @param[inout] molecule  : molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine setup_molecule_no(mole_info, molecule)

    ! formal argments
    type(s_mole_info),        intent(in)    :: mole_info
    type(s_molecule), target, intent(inout) :: molecule

    ! local variables
    integer                   :: alloc_stat
    integer                   :: i, j, k, natom, nresi, nmole
    integer                   :: rnob, rnoe, si1, si2
    integer                   :: mid, ngr
    integer                   :: ty, rn
    character(4)              :: snb, sne
    character(6)              :: rnab, rnae
    character(6)              :: gr(MaxGroupResName)

    integer,      allocatable :: resno(:), atomidx(:,:), mol_tbl(:,:)
    character(4), allocatable :: segname(:)
    character(6), allocatable :: resname(:)

    integer,          pointer :: presno(:)
    character(4),     pointer :: psegname(:)
    character(6),     pointer :: presname(:)


    natom = molecule%num_atoms
    if (natom == 0) then
      call alloc_molecules(molecule, MoleculeMole, 0)
      return
    end if

    psegname => molecule%segment_name
    presname => molecule%residue_name
    presno   => molecule%residue_no

    ! allocate memory
    nresi = 1
    do i = 2, natom
      if (psegname(i-1) /= psegname(i) .or. &
          presno  (i-1) /= presno  (i)) &
          nresi = nresi + 1
    end do

    allocate(segname(nresi),   &
             resno  (nresi),   &
             resname(nresi),   &
             atomidx(2,nresi), &
             mol_tbl(3,natom), &
             stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    ! prepare residue lists
    nresi = 1
    segname(nresi)   = psegname(1)
    resno  (nresi)   = presno  (1)
    resname(nresi)   = presname(1)
    atomidx(1,nresi) = 1

    do i = 2, natom
      if (psegname(i-1) /= psegname(i) .or. &
          presno  (i-1) /= presno  (i)) then
        atomidx(2,nresi) = i - 1
        nresi = nresi + 1
        segname(nresi)   = psegname(i)
        resno  (nresi)   = presno  (i)
        resname(nresi)   = presname(i)
        atomidx(1,nresi) = i
      end if
    end do
    atomidx(2,nresi) = natom

    mid = 1
    mol_tbl(:,:) = 0

    ! prepare range residue
    do i = 1, mole_info%num_range

      snb  = mole_info%range_res(i)%seg_name_b
      sne  = mole_info%range_res(i)%seg_name_e
      rnob = mole_info%range_res(i)%res_no_b
      rnoe = mole_info%range_res(i)%res_no_e
      rnab = mole_info%range_res(i)%res_name_b
      rnae = mole_info%range_res(i)%res_name_e

      si1 = -1
      do j = 1, nresi

        if (('' == snb .or. segname(j) == snb) .and. &
             resno(j) == rnob .and. &
             resname(j) == rnab) then
          si1 = atomidx(1,j)
        end if

        if (si1 == -1) &
          cycle

        if (('' == sne .or. segname(j) == sne) .and. &
             resno(j) == rnoe .and. &
             resname(j) == rnae) then
          si2 = atomidx(2,j)

          mid = mid + 1
          do k = si1, si2
            mol_tbl(:,k) = (/1,i,mid/)
          end do

          si1 = -1
        end if

      end do
    end do

    ! prepare group residue
    do i = 1, mole_info%num_group
      ngr      = mole_info%group_res(i)%num_res
      gr(:ngr) = mole_info%group_res(i)%res_name(:ngr)

      j = 1
      do while(j <= nresi)

        do k = 1, ngr
          if (nresi < j+k-1) &
            exit
          if (resname(j+k-1) /= gr(k)) &
            exit
        end do

        if (k > ngr) then
          si1 = atomidx(1,j)
          si2 = atomidx(2,j+ngr-1)

          mid = mid + 1
          do k = si1, si2
            if (mol_tbl(TblType,k) == 0) &
              mol_tbl(:,k) = (/2,i,mid/)
          end do
          j = j + ngr
        else
          j = j + 1
        end if

      end do
    end do

    ! setup molecule number
    nmole = 1
    do i = 2, natom
      if (mol_tbl(TblMolID,i-1) /= mol_tbl(TblMolID,i)) &
        nmole = nmole + 1
    end do

    call alloc_molecules(molecule, MoleculeMole, nmole)
    
    nmole = 1
    ty = mol_tbl(TblType, 1)
    rn = mol_tbl(TblMolNo, 1)

    select case(ty)
    case (0)
      molecule%molecule_name(nmole) = '@UNDEFINED@'
    case (1)
      molecule%molecule_name(nmole) = mole_info%range_res(rn)%name
    case (2)
      molecule%molecule_name(nmole) = mole_info%group_res(rn)%name
    end select
    molecule%molecule_atom_no(nmole) = 1
    molecule%molecule_no(1) = nmole

    do i = 2, natom
      if (mol_tbl(TblMolID,i-1) /= mol_tbl(TblMolID,i)) then

        nmole = nmole + 1
        ty = mol_tbl(TblType, i)
        rn = mol_tbl(TblMolNo, i)
        select case(ty)
        case (0)
          molecule%molecule_name(nmole) = '@UNDEFINED@'
        case (1)
          molecule%molecule_name(nmole) = mole_info%range_res(rn)%name
        case (2)
          molecule%molecule_name(nmole) = mole_info%group_res(rn)%name
        end select
        molecule%molecule_atom_no(nmole) = i

      end if
      molecule%molecule_no(i) = nmole
    end do

    molecule%num_molecules = nmole

    ! deallocate memory
    deallocate(segname, resno, resname, atomidx, mol_tbl, &
               stat = alloc_stat)

    return

  end subroutine setup_molecule_no

end module select_mod
