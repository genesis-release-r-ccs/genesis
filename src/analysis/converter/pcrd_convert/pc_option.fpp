!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pc_option_mod
!> @brief   module for conversion options
!! @authors Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module pc_option_mod

  use constants_mod
  use pc_option_str_mod
  use pbc_correct_mod
  use trajectory_str_mod
  use select_mod
  use select_atoms_mod
  use molecules_str_mod
  use fileio_mod
  use fileio_control_mod
  use string_mod
  use messages_mod

  implicit none
  private

  ! structures
  type, public :: s_opt_info
    logical                         :: check_only    = .false.
    logical                         :: allow_backup  = .false.
    integer                         :: trjout_format = TrjFormatDCD
    integer                         :: trjout_type   = TrjTypeCoorBox
    integer                         :: trjout_atom   = 1
    integer                         :: pbc_correct   = PBCCModeNo
    character(MaxLine), allocatable :: rename_res(:)
  end type s_opt_info

  ! subroutines
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
    write(MsgOut,'(A)') 'pbc_correct    = NO              # (NO/MOLECULE)'
    write(MsgOut,'(A)') '# rename_res1  = HSE HIS'
    write(MsgOut,'(A)') '# rename_res2  = HSD HIS'
    write(MsgOut,'(A)') ' '

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

    ! local variables
    integer                  :: i, nresi
    character(Maxline)       :: value
    character(Maxline)       :: rename_res


    ! read parameters
    !

    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_logical(handle, Section, 'check_only', &
                               opt_info%check_only)

    call read_ctrlfile_logical(handle, Section, 'allow_backup', &
                               opt_info%allow_backup)

    call read_ctrlfile_type   (handle, Section, 'trjout_format', &
                               opt_info%trjout_format, TrjFormatTypes)

    call read_ctrlfile_type   (handle, Section, 'trjout_type', &
                               opt_info%trjout_type, TrjTypeTypes)

    call read_ctrlfile_integer(handle, Section, 'trjout_atom', &
                               opt_info%trjout_atom)

    call read_ctrlfile_type   (handle, Section, 'pbc_correct', &
                               opt_info%pbc_correct, PBCCModeTypes)

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
    write(MsgOut,'(A20,A10)') &
                    '  trjout format   = ', &
                    TrjFormatTypes(opt_info%trjout_format)
    write(MsgOut,'(A20,A10)') &
                    '  trjout type     = ', &
                    TrjTypeTypes(opt_info%trjout_type)
    write(MsgOut,'(A20,A5,I0)') &
                    '  trjout atom     = ', 'group', opt_info%trjout_atom
    write(MsgOut,'(A20,A10)') &
                    '  pbc correction  = ', &
                    PBCCModeTypes(opt_info%pbc_correct)
    write(MsgOut,'(A20,I10)') '  # of rename res = ', nresi
    do i = 1, nresi
      write(MsgOut,'(A10,I0,A3,A)') &
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
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine setup_option(opt_info, &
                          sel_info, &
                          molecule, &
                          option)
  
    ! formal argments
    type(s_opt_info),        intent(in)    :: opt_info
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_molecule),        intent(inout) :: molecule
    type(s_option),          intent(inout) :: option

    ! local variables
    integer                 :: i, j
    integer                 :: nresi, nrenm
    character(MaxLine)      :: s1, s2
    character(6)            :: st


    ! check only
    option%check_only    = opt_info%check_only

    ! backup output file
    call setup_backup(opt_info%allow_backup)

    ! trajectory format
    option%trjout_format = opt_info%trjout_format

    ! trajectory type
    option%trjout_type = opt_info%trjout_type

    ! trajectory atom
    if (opt_info%trjout_atom > size(sel_info%groups)) &
        call error_msg('Setup_Option> trj-out atom selection is out of range.')

    option%trjout_atom_exp = sel_info%groups(opt_info%trjout_atom)

    call select_atom(molecule, &
                     option%trjout_atom_exp, &
                     option%trjout_atom)

    write(MsgOut,'(A,I8)') 'Setup_Option> trj-out atom count: ', &
         size(option%trjout_atom%idx)

    ! pbc correct
    option%pbcc_mode = opt_info%pbc_correct

    ! setup for PBC-correct mode "molecule"
    call setup_pbc_correct(option%pbcc_mode, molecule)

    ! rename residues of molecule informations
    nresi = size(opt_info%rename_res)
    nrenm = 0
    do i = 1, nresi
      
      s1 = ''
      s2 = ''
      read(opt_info%rename_res(i), *, err=900, end=900) s1, s2

      do j = 1, size(molecule%residue_name)
        st = molecule%residue_name(j)
        if (st == s1) then
          molecule%residue_name(j) = s2
          nrenm = nrenm + 1
        end if
      end do
    end do

    if (nresi > 0) &
      write(MsgOut,'(A,I8)') 'Setup_Option> re-named residues : ', nrenm

    write(MsgOut,'(A)') ' '

    return

900 call error_msg('Setup_Option> Rename residue : Bad format:'//&
                   trim(opt_info%rename_res(i)))

  end subroutine setup_option

end module pc_option_mod
