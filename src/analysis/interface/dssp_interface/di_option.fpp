!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   di_option_mod
!> @brief   module for analysis options
!! @authors Takaharu Mori (TM)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module di_option_mod

  use di_option_str_mod
  use select_mod
  use select_atoms_mod
  use molecules_str_mod
  use fileio_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_opt_info
    logical                         :: check_only     = .false.
    logical                         :: allow_backup   = .false.
    integer                         :: analysis_atom  = 1
    character(MaxLine)              :: dssp_exec      = '/home/user/bin/dssp.exe'
    character(MaxLine)              :: temporary_pdb  = 'temporary.pdb'
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
  !! @authors      TM
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_option

    write(MsgOut,'(A)') '[OPTION]'
    write(MsgOut,'(A)') 'check_only     = NO                      # only checking input files (YES/NO)'
    write(MsgOut,'(A)') 'allow_backup   = NO                      # backup existing output files (YES/NO)'
    write(MsgOut,'(A)') 'analysis_atom  = 1                       # protein atom group'
    write(MsgOut,'(A)') 'dssp_exec      = /home/user/bin/dssp.exe # dssp binary path'
    write(MsgOut,'(A)') 'temporary_pdb  = temporary.pdb           # temporary pdb for DSSP'
    write(MsgOut,'(A)') '# rename_res1    = HSE HIS'
    write(MsgOut,'(A)') '# rename_res2    = HSD HIS'
    write(MsgOut,'(A)') ''

    return

  end subroutine show_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_option
  !> @brief        read control parameters in OPTION section
  !! @authors      TM
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
    character(MaxLine)       :: value
    character(MaxLine)       :: rename_res


    ! read parameters
    !
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_logical(handle, Section, &
                              'check_only', opt_info%check_only)

    call read_ctrlfile_logical(handle, Section, &
                              'allow_backup', opt_info%allow_backup)

    call read_ctrlfile_integer(handle, Section, &
                              'analysis_atom', opt_info%analysis_atom)

    call read_ctrlfile_string (handle, Section, &
                              'dssp_exec', opt_info%dssp_exec)

    call read_ctrlfile_string (handle, Section, &
                              'temporary_pdb', opt_info%temporary_pdb)

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
    write(MsgOut,'(A20,A5,I0)') &
                    '  analysis atom   = ', 'group', opt_info%analysis_atom

    write(MsgOut,'(A20,A)')    '  dssp_exec       = ', trim(opt_info%dssp_exec)
    write(MsgOut,'(A20,A)')    '  temporary_pdb   = ', trim(opt_info%temporary_pdb)

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
  !! @authors      TM
  !! @param[in]    opt_info : OPTION section control parameters information
  !! @param[in]    sel_info : SELECTION section control parameters information
  !! @param[in]    molecule : molecule information
  !! @param[inout] option   : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine setup_option(opt_info, sel_info, molecule, option)

    ! formal arguments
    type(s_opt_info),        intent(in)    :: opt_info
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_molecule),        intent(inout) :: molecule
    type(s_option),          intent(inout) :: option

    ! local variables
    integer                  :: i, j
    integer                  :: nresi, nrenm, ndata
    character(MaxLine)       :: s1, s2
    character(6)             :: st

    ! check only
    !
    option%check_only = opt_info%check_only

    ! backup output file
    !
    call setup_backup(opt_info%allow_backup)

    ! analysis atom
    !
    if (opt_info%analysis_atom > size(sel_info%groups)) &
        call error_msg('Setup_Option> analysis atom selection is out of range.')

    option%analysis_atom_exp = sel_info%groups(opt_info%analysis_atom)

    call select_atom(molecule, &
                     option%analysis_atom_exp, &
                     option%analysis_atom)

    write(MsgOut,'(A,I8)') 'Setup_Option> analysis atom count: ', &
         size(option%analysis_atom%idx)

    ! dssp command
    !
    option%dssp_exec     = opt_info%dssp_exec
    option%temporary_pdb = opt_info%temporary_pdb

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

end module di_option_mod
