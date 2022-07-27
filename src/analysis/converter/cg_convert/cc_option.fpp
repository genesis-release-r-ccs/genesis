!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cc_option_mod
!> @brief   module for trajectory conversion options
!! @authors Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module cc_option_mod

  use constants_mod
  use cc_option_str_mod
  use pbc_correct_mod
  use trajectory_str_mod
  use select_mod
  use select_atoms_mod
  use select_atoms_str_mod
  use molecules_str_mod
  use fileio_mod
  use fileio_control_mod
  use string_mod
  use messages_mod

  implicit none
  private

  ! structures
  type, public :: s_opt_info
    logical                         :: check_only     = .false.
    logical                         :: allow_backup   = .false.
    logical                         :: geometric_flag = .false.
    integer                         :: trjout_format  = TrjFormatDCD
    integer                         :: trjout_type    = TrjTypeCoorBox
    integer                         :: trjout_atom    = 1
    character(MaxLine)              :: trjout_atom_char  = ""
    integer                         :: pbc_correct    = PBCCModeNo
    logical                         :: split_trjpdb   = .false.
    logical                         :: centering      = .false.
    integer                         :: centering_atom = 1
    integer                         :: num_outgroups  = 1
    character(MaxLine)              :: center_coord   = ''
    character(MaxLine), allocatable :: rename_res(:)
  end type s_opt_info

  public  :: show_ctrl_option
  public  :: read_ctrl_option
  public  :: setup_option

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_option
  !> @brief        show control parameters in OPTION section
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_option

    write(MsgOut,'(A)') '[OPTION]'
    write(MsgOut,'(A)') 'check_only     = NO              # only checking input files (YES/NO)'
    write(MsgOut,'(A)') 'allow_backup   = NO              # backup existing output files (YES/NO)'
    write(MsgOut,'(A)') 'trjout_format  = DCD             # (PDB/DCD)'
    write(MsgOut,'(A)') 'trjout_type    = COOR+BOX        # (COOR/COOR+BOX)'
    write(MsgOut,'(A)') 'trjout_atom    = 1               # atom group'
    write(MsgOut,'(A)') 'geometric      = YES             # no mass weight'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '# Note: Coordinates are converted in the following order:'
    write(MsgOut,'(A)') '#       read traj -> write traj'
    write(MsgOut,'(A)') ''

    return

  end subroutine show_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_option
  !> @brief        read control parameters in OPTION section
  !! @authors      NT
  !! @param[in]    handle   : unit number of control file
  !! @param[out]   opt_info : OPTION section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine read_ctrl_option(handle, opt_info)
  
    ! parameters
    character(*),            parameter     :: Section = 'Option'

    ! formal argments
    integer,                 intent(in)    :: handle
    type(s_opt_info),        intent(inout) :: opt_info

    ! clocal variables
    integer                  :: i, nresi
    character(MaxLine)       :: value
    character(MaxLine)       :: rename_res


    ! read parameters
    !

    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_logical(handle, Section, 'check_only',     &
                               opt_info%check_only)

    call read_ctrlfile_logical(handle, Section, 'allow_backup',   &
                               opt_info%allow_backup)

    call read_ctrlfile_type   (handle, Section, 'trjout_format',  &
                               opt_info%trjout_format, TrjFormatTypes)

    call read_ctrlfile_type   (handle, Section, 'trjout_type',    &
                               opt_info%trjout_type, TrjTypeTypes)

    call read_ctrlfile_string(handle, Section, 'trjout_atom',    &
                               opt_info%trjout_atom_char)

    call read_ctrlfile_logical(handle, Section, 'centering',      &
                               opt_info%centering)

    call read_ctrlfile_logical(handle, Section, 'geometric',      &
                               opt_info%geometric_flag)

    call read_ctrlfile_integer(handle, Section, 'centering_atom', &
                               opt_info%centering_atom)

    call read_ctrlfile_string (handle, Section, 'center_coord',   &
                               opt_info%center_coord)

    call read_ctrlfile_type   (handle, Section, 'pbc_correct',    &
                               opt_info%pbc_correct, PBCCModeTypes)

    call read_ctrlfile_logical(handle, Section, 'split_trjpdb',   &
                               opt_info%split_trjpdb)


    ! read rename_res parameters
    nresi = 0

    do while (.true.)

      value = ''
      write(rename_res,'(A10,I0)') 'rename_res', nresi + 1
      call read_ctrlfile_string(handle, Section, rename_res, value)

      if (value == '') &
        exit

      nresi = nresi + 1

    end do

    allocate(opt_info%rename_res(nresi))
 
    do i = 1, nresi
      write(rename_res,'(A10,I0)') 'rename_res', i
      call read_ctrlfile_string(handle, Section, rename_res, &
                                opt_info%rename_res(i))
    end do

    call end_ctrlfile_section(handle)


    ! write parameters to MsgOut
    !
    write(MsgOut,'(A)') 'Read_Ctrl_Option> Parameters of Options'

    if (opt_info%check_only) then
      write(MsgOut,'(A20,A3)') '  check only      = ', 'yes'
    else
      write(MsgOut,'(A20,A2)') '  check only      = ', 'no'
    end if

    if (opt_info%allow_backup) then
      write(MsgOut,'(A20,A3)') '  allow backup    = ', 'yes'
    else
      write(MsgOut,'(A20,A2)') '  allow backup    = ', 'no'
    end if

    write(MsgOut,'(A20,A10)')     &
                 '  trjout format   = ', &
                 TrjFormatTypes(opt_info%trjout_format)
    write(MsgOut,'(A20,A10)')     &
                 '  trjout type     = ', &
                 TrjTypeTypes(opt_info%trjout_type)
    write(MsgOut,'(A20,A6,A)')   &
                 '  trjout atom     = ', &
                 'groups', trim(opt_info%trjout_atom_char)

    if (opt_info%centering) then
      write(MsgOut,'(A20,A3)')    &
                 '  centering       = ', 'yes'
      write(MsgOut,'(A20,A5,I0)') &
                 '  centering atom  = ', &
                 'group', opt_info%centering_atom
      write(MsgOut,'(A20,A)')     &
                 '  center_coord    = ', trim(opt_info%center_coord)
    else
      write(MsgOut,'(A20,A2)')    &
                 '  centering       = ', 'no'
    end if

    if (opt_info%split_trjpdb) then
      write(MsgOut,'(A20,A3)')    &
                 '  split_trjpdb    = ', 'yes'
    else
      write(MsgOut,'(A20,A2)')    &
                 '  split_trjpdb    = ', 'no'
    end if

    write(MsgOut,'(A20,A10)')     &
                 '  pbc correction  = ', &
                 PBCCModeTypes(opt_info%pbc_correct)

    write(MsgOut,'(A20,I10)') '  # of rename res = ', nresi
    do i = 1, nresi
      write(MsgOut,'(A10,I0,A3,A)')      &
                 '   rename ', i, ' = ', trim(opt_info%rename_res(i))
    end do

    write(MsgOut,'(A)') ' '


    return

  end subroutine read_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_option
  !> @brief        setup option information
  !! @authors      NT
  !! @param[in]    opt_info : OPTION section control parameters information
  !! @param[in]    sel_info : SELECTION section control parameters information
  !! @param[inout] molecule : molecule information
  !! @param[out]   option   : option information
  !! @param[inout] molecule_out : molecule information
  !! @param[inout] trajectory_out : trajectory information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine setup_option(opt_info, &
                          sel_info, &
                          molecule, &
                          option,   &
                          molecule_out,   &
                          trajectory_out)
  
    ! formal argments
    type(s_opt_info),        intent(in)    :: opt_info
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_molecule),        intent(inout) :: molecule
    type(s_option),          intent(inout) :: option
    type(s_molecule),        intent(inout) :: molecule_out
    type(s_trajectory),      intent(inout) :: trajectory_out

    ! local variables
    integer                  :: i, j, k
    integer                  :: nresi, nrenm, ndata, ngroup, max_natom, natom
    integer                  :: igroup
    integer                  :: alloc_stat
    character(MaxLine)       :: s1, s2
    character(6)             :: st
    integer                  :: iatm, ires, oldres
    character(4)             :: segname, oldsegnm
    character(6)             :: resname, oldresnm
    character(4)             :: atname

    type(s_selatoms), allocatable :: seleatoms(:)


    ! check only
    option%check_only    = opt_info%check_only

    ! backup output file
    call setup_backup(opt_info%allow_backup)

    ! trajectory format
    option%trjout_format = opt_info%trjout_format

    ! trajectory type
    option%trjout_type = opt_info%trjout_type

    ! trajectory atom
    option%num_groups = split_num(opt_info%trjout_atom_char)
    ngroup = size(sel_info%groups)
    if (ngroup > 0 .and. option%num_groups > 0) then
      allocate(seleatoms(1:ngroup),                            &
               option%trjout_groups(1:option%num_groups),      &
               trajectory_out%coord(1:3,1:option%num_groups),  &
               stat=alloc_stat)
      option%trjout_groups(1:option%num_groups) = 0

      call alloc_molecules(molecule_out, MoleculeAtom, option%num_groups)
      molecule_out%chain_id(1:option%num_groups)         = 'A'
      molecule_out%num_atoms                             = option%num_groups
      molecule_out%residue_name(1:option%num_groups)     = 'ALA   '
      molecule_out%atom_name(1:option%num_groups)        = ' CA '
      molecule_out%segment_name(1:option%num_groups)     = 'PROT'
      do i = 1, option%num_groups
        molecule_out%residue_no(i)     = i
      end do

      if (alloc_stat /= 0) &
        call error_msg('Setup_Option> Memory allocation error.')

       call split(option%num_groups, option%num_groups,  &
                  opt_info%trjout_atom_char, option%trjout_groups)

      max_natom = 0
      do i = 1, ngroup
        call select_atom(molecule, sel_info%groups(i), seleatoms(i))
        natom = size(seleatoms(i)%idx)
        max_natom = max (max_natom, natom)
      end do

      allocate(option%num_atoms(1:option%num_groups),              &
               option%atomlist(1:max_natom, 1:option%num_groups),  &
               stat=alloc_stat)
      option%atomlist(1:max_natom, 1:option%num_groups) = 0
      option%num_atoms(1:option%num_groups) = 0

      write(MsgOut,'(A,I0)') 'Setup_Option> number of groups: ', &
                                   option%num_groups
      do i = 1, option%num_groups
        igroup = option%trjout_groups(i)
        if (igroup > ngroup .or. igroup <= 0) &
          call error_msg('Setup_Option> trj-out atom selection is out of range.')
        option%num_atoms(i) = size(seleatoms(igroup)%idx)
        do j = 1, option%num_atoms(i)
          option%atomlist(j,i) = seleatoms(igroup)%idx(j)
        end do

        write(MsgOut,'(A,I5,A,I0,A)') ' group = ', i, ', "', option%trjout_groups(i),'"'
        write(MsgOut,'(A,I5)')   ' # of atoms = ', option%num_atoms(i)
        write(MsgOut,'(A)')      ' atomlist: '
        oldres   =-10000
        oldsegnm = "    "
        oldresnm = "      "
        do j = 1, option%num_atoms(i)
          iatm=option%atomlist(j,i)
          if (iatm <= 0 .or. iatm > molecule%num_atoms) &
            call error_msg('Setup_Option> atom idx is not correct.')
          segname = molecule%segment_name(iatm)
          resname = molecule%residue_name(iatm)
          atname  = molecule%atom_name(iatm)
          ires    = molecule%residue_no(iatm)
          if (ires .ne. oldres .or. resname .ne. oldresnm .or. &
             segname .ne. oldsegnm) then
            if (oldres .ne. -10000) then
              write(MsgOut,'(A)') ''
            endif
            write(MsgOut,'(5X,A4,A,A6,A,i6,A,$)') &
                trim(segname),'.',trim(resname),'.',ires,':'
            write(MsgOut,'(1X,A,$)') trim(atname)
            oldres = ires
            oldresnm = resname
            oldsegnm = segname
          else
            write(MsgOut,'(A,1X,A,$)') ',',trim(atname)
          endif
        end do
        write(MsgOut,'(A)') 
      end do
       
    else
      call error_msg('Setup_Option> selection is required.')
    endif

    if (opt_info%geometric_flag) then
      molecule%mass(1:molecule%num_atoms) = 1.0_wp
    endif

    write(MsgOut,'(A)') ' '

    return

900 call error_msg('Setup_Option> Rename residue : Bad format:'//&
                   trim(opt_info%rename_res(i)))

  end subroutine setup_option

end module cc_option_mod
