!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   qg_option_mod
!> @brief   module for trajectory extraction options
!! @authors Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module qg_option_mod

  use qg_option_str_mod
  use pbc_correct_mod
  use trajectory_str_mod
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
    logical                         :: check_only          = .false.
    logical                         :: allow_backup        = .false.
    integer                         :: coord_format        = CrdFormatCharmm
    character(MaxLine)              :: qmmm_atom_index     = ''
    character(MaxLine)              :: qm_atom_index       = ''
    logical                         :: reconcile_num_atoms = .false.
    character(MaxLine)              :: origin_atom_index   = ''
    character(MaxLine)              :: frame_number        = '1'
    character(MaxLine), allocatable :: rename_res(:)
  end type s_opt_info

  public  :: show_ctrl_option
  public  :: read_ctrl_option
  public  :: setup_option
  private :: assign_number
  private :: check_around
  private :: setup_wrap_qmmm

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_option
  !> @brief        show control parameters in OPTION section
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_option

    write(MsgOut,'(A)')  '[OPTION]'
    write(MsgOut,'(A)')  'check_only          = NO          # only checking input files (YES/NO)'
    write(MsgOut,'(A)')  'allow_backup        = NO          # backup existing output files (YES/NO)'
    write(MsgOut,'(A)')  'coord_format        = CHARMM      # (CHARMM CARD)'
    write(MsgOut,'(A)')  'qmmm_atom_index     = 1 2 3 4     # atom groups'
    write(MsgOut,'(A)')  '# qm_atom_index     = 1 2         # atom groups'
    write(MsgOut,'(A)')  'frame_number        = 1,3,5,7,9   # A. 1,3,5,7,9; B. 1:2:10; C. 1,3 / 5:2:9'
    write(MsgOut,'(A)')  'reconcile_num_atoms = NO          # (YES/NO)'
    write(MsgOut,'(A)')  'origin_atom_index   = 1 2         # atom groups'
    write(MsgOut,'(A)')  ' '

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

    ! formal arguments
    integer,                 intent(in)    :: handle
    type(s_opt_info),        intent(inout) :: opt_info

    ! clocal variables
    integer                  :: i, nresi
    character(MaxLine)       :: value
    character(MaxLine)       :: rename_res


    ! read parameters
    !

    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_logical(handle, Section, 'check_only',          &
                               opt_info%check_only)

    call read_ctrlfile_logical(handle, Section, 'allow_backup',        &
                               opt_info%allow_backup)

    call read_ctrlfile_type   (handle, Section, 'coord_format',        &
                               opt_info%coord_format, CrdFormatTypes)

    call read_ctrlfile_string (handle, Section, 'qmmm_atom_index',     &
                               opt_info%qmmm_atom_index)

    call read_ctrlfile_string (handle, Section, 'qm_atom_index',       &
                               opt_info%qm_atom_index)

    call read_ctrlfile_string (handle, Section, 'frame_number',        &
                               opt_info%frame_number)

    call read_ctrlfile_logical(handle, Section, 'reconcile_num_atoms', &
                               opt_info%reconcile_num_atoms)

    call read_ctrlfile_string (handle, Section, 'origin_atom_index',   &
                               opt_info%origin_atom_index)

    ! read rename_res parameters
    nresi = 0

    do while (.true.)

      value = ''
      write(rename_res,'(A10,I0)') 'rename_res', nresi + 1
      call read_ctrlfile_string(handle, Section, rename_res, value)

      if (value .eq. '') &
        exit

      nresi = nresi + 1

    end do

    if(allocated(opt_info%rename_res)) deallocate(opt_info%rename_res)
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
      write(MsgOut,'(A20,A3)') '  Check only      = ', 'YES'
    else
      write(MsgOut,'(A20,A2)') '  Check only      = ', 'NO'
    end if

    if (opt_info%allow_backup) then
      write(MsgOut,'(A20,A3)') '  allow backup    = ', 'yes'
    else
      write(MsgOut,'(A20,A2)') '  allow backup    = ', 'no'
    end if

    write(MsgOut,'(A20,A10)')            &
                 '  Output format   = ', &
                 CrdFormatTypes(opt_info%coord_format)
    write(MsgOut,'(A20,A)')              &
                 '  QM/MM atoms     = ', &
                 'group: '//trim(opt_info%qmmm_atom_index)

    if(len_trim(opt_info%qm_atom_index) /= 0) &
      write(MsgOut,'(A20,A)')                 &
                 '  QM atoms        = ',      &
                 'group: '//trim(opt_info%qm_atom_index)

    if(len_trim(opt_info%origin_atom_index) /= 0) &
      write(MsgOut,'(A20,A)')                     &
                 '  Origin defined  = ',          &
                 'group: '//trim(opt_info%origin_atom_index)

    if (opt_info%reconcile_num_atoms) then
      write(MsgOut,'(A20,A3)') '  Reconciliation  = ', 'YES'
    else
      write(MsgOut,'(A20,A2)') '  Reconciliation  = ', 'NO'
    end if

    write(MsgOut,'(A20,A)')     &
                 '  Frame number    = ', trim(opt_info%frame_number)

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
  !! @authors      NT, KYMD
  !! @param[in]    opt_info : OPTION section control parameters information
  !! @param[in]    trj_list : TRAJECTORY section control parameters information
  !! @param[in]    sel_info : SELECTION section control parameters information
  !! @param[inout] molecule : molecule information
  !! @param[out]   option   : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine setup_option(opt_info, &
                          sel_info, &
                          trj_list, &
                          molecule, &
                          option)
  
    ! formal arguments
    type(s_opt_info),        intent(in)    :: opt_info
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_trj_list),        intent(in)    :: trj_list
    type(s_molecule),        intent(inout) :: molecule
    type(s_option),          intent(inout) :: option

    ! local variables
    integer                  :: i, j
    integer                  :: nresi, nrenm, ndata
    character(MaxLine)       :: s1, s2
    character(6)             :: st

    integer                  :: nframe_all, nframe_save
    integer                  :: irun, nrun, itrj, ntrj
    integer                  :: ngrp_sel, ngrp_qmmm, ngrp_qm, ngrp_cent
    integer                  :: jcheck_around


    ! check only
    option%check_only   = opt_info%check_only

    ! backup output file
    call setup_backup(opt_info%allow_backup)

    ! trajectory format
    option%coord_format = opt_info%coord_format

    ! check if groups are defined
    if(size(sel_info%groups) == 0) &
      call error_msg('Setup_Option> index in [SELECTION] must be required')

    ! QM/MM atoms
    !
    !if (len_trim(opt_info%qmmm_atom_index) == 0) &
    !    call error_msg('Setup_Option> No QM/MM atoms defined')
    if (len_trim(opt_info%qmmm_atom_index) /= 0) then

      ! number of QM and MM atoms
      ngrp_sel  = size(sel_info%groups)
      ngrp_qmmm = split_num(trim(opt_info%qmmm_atom_index))
      if (ngrp_qmmm > ngrp_sel) &
          call error_msg('Setup_Option> QM/MM atom selection is out of range.')

      ! merge selection group
      call merge_selection(sel_info%groups, option%qmmm_atom_exp,  &
                           opt_info%qmmm_atom_index, ngrp_sel, ngrp_qmmm)

    else

      write(option%qmmm_atom_exp, '(a)') "all"

    end if
    
    ! assign selected atoms
    call select_atom(molecule, &
                     option%qmmm_atom_exp, &
                     option%qmmm_atom)

    write(MsgOut,'(A,I8)') 'Setup_Option> # of QMMM atoms = ', &
         size(option%qmmm_atom%idx)

    ! QM-subsystem atoms
    !
    !if (len_trim(opt_info%qm_atom_index) == 0) &
    !  call error_msg('Setup_Option> No QM atoms defined')
    if (len_trim(opt_info%qm_atom_index) /= 0) then

      ngrp_qm = split_num(trim(opt_info%qm_atom_index))
      if (ngrp_qm > ngrp_sel) &
        call error_msg('Setup_Option> QM atom selection is out of range.')

      ! merge selection group
      call merge_selection(sel_info%groups, option%qm_atom_exp, &
                           opt_info%qm_atom_index, ngrp_sel, ngrp_qm)

      ! assign selected atoms
      call select_atom(molecule, &
                       option%qm_atom_exp, &
                       option%qm_atom)

      write(MsgOut,'(A,I8)') 'Setup_Option> # of QM atoms   = ', &
        size(option%qm_atom%idx)

    end if

    if (len_trim(opt_info%origin_atom_index) /= 0) then

      ngrp_cent = split_num(trim(opt_info%origin_atom_index))
      if (ngrp_cent > ngrp_sel) &
        call error_msg(       &
            'Setup_Option> selection for origin definition is out of range.')

      ! merge selection group
      call merge_selection(sel_info%groups, option%origin_atom_exp, &
                           opt_info%origin_atom_index, ngrp_sel, ngrp_cent)

      ! assign selected atoms
      call select_atom(molecule, &
                       option%origin_atom_exp, &
                       option%origin_atom)

      write(MsgOut,'(A,I8)') 'Setup_Option> # of atoms to define center = ', &
        size(option%origin_atom%idx)

    else

      if (len_trim(opt_info%qm_atom_index) == 0) &
        call error_msg('Setup_Option> Neither QM atoms nor Origin atoms defined')

      ! assign selected atoms
      call select_atom(molecule, &
                       option%qm_atom_exp, &
                       option%origin_atom)

    end if


    ! Check whether around_xxx (atoms, residues, molecules) is used 
    !
    jcheck_around = 0
    call check_around(option%qmmm_atom_exp, jcheck_around)
    call check_around(option%qm_atom_exp, jcheck_around)
    if(mod(jcheck_around, 10) /= 0) &
      call error_msg('ERROR: Use of around and around_atoms is not permitted. &
                     &Use around_res or around_mol instead.')
    if(mod(jcheck_around, 100) / 10 /= 0) then
       write(MsgOut,'(A)')      'Setup_Option> **NOTICE** ' 
       write(MsgOut,'(14X, A)') 'Around_residues for protein is not recommended.'
       write(MsgOut,'(14X, A)') 'Please use qmmm_psffile at your own risk.'
       write(MsgOut,'(A)')      'Setup_Option> **NOTICE** ' 
    end if

    ! Reconsile
    option%reconcile_num_atoms = opt_info%reconcile_num_atoms
    if (jcheck_around / 10 /= 0 .and. option%reconcile_num_atoms) then
      write(MsgOut,'(A)') 'Setup_Option> Reconciliation of # of atoms : DO'
    else
      option%reconcile_num_atoms = .false.
      write(MsgOut,'(A)') 'Setup_Option> Reconciliation of # of atoms : SKIP'
    end if


    ! Frames to be saved
    !
    ! Currently CHARMM CRD is only permitted for output format 
    if (opt_info%frame_number /= '') then 
      if (option%coord_format /= CrdFormatCharmm) &
        call error_msg('Setup_Option> coord_format should be CHARMM .crd')
    end if

    ! Assign and count frames taken from trajectory in a form of CRD
    nframe_all = 0
    nrun = size(trj_list%md_steps)
    do irun = 1, nrun
      ntrj = trj_list%md_steps(irun)
      do itrj = 1, ntrj
        if (mod(itrj, trj_list%ana_periods(irun)) == 0) &
          nframe_all = nframe_all + 1
      end do
    end do
    call alloc_option(option, nframe_all)
    call assign_number(option%frame_no, opt_info%frame_number, &
                       nframe_save)
    write(MsgOut,'(A,I4)')                                     &
     'Setup_Option> # of frames to be saved = ', nframe_save


    ! setup for wrapping solvent molecules and so on
    !
    call setup_wrap_qmmm(molecule)


    ! rename residues of molecule informations
    !
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


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_number
  !> @brief        assign frame number to be saved
  !! @authors      KYMD
  !! @param[out]   frame_no     : frame number
  !! @param[in]    frame_number : frame number given in option
  !! @param[out]   nframe_save  : # of frames to be saved
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_number(frame_no, frame_number, nframe_save)

    ! formal
    integer, intent(inout)       :: frame_no(:)
    character(len=*), intent(in) :: frame_number
    integer, intent(out)         :: nframe_save

    ! local
    integer                      :: nframes
    character(MaxLine)           :: line_type1, line_type2
    integer                      :: numlen, iposf, iposf_t2, istart
    integer                      :: iframe
    integer                      :: ibegin, istride, iend 


    nframes    = size(frame_no(:))
    line_type1 = ''
    line_type2 = ''

    ! allowed formats are 
    !    1. 1, 2, 3, ..., 10                
    !    2. 1:1:10
    !    3. 1, 2, 3 / 4:1:10
    numlen = len_trim(frame_number)
    iposf  = index(frame_number(1: numlen), '/')
    if(iposf /= 0) then
      iposf_t2 = index(frame_number(1: iposf - 1), ':')
      if(iposf_t2 == 0) then
        write(line_type1, '(a)') trim(frame_number(1: iposf - 1))//","
        write(line_type2, '(a)') trim(frame_number(iposf + 1: numlen))//":"
      else
        write(line_type2, '(a)') trim(frame_number(1: iposf - 1))//":"
        write(line_type1, '(a)') trim(frame_number(iposf + 1: numlen))//","
      end if

    else
      iposf_t2 = index(frame_number(1: numlen), ':')
      if(iposf_t2 == 0) then
        write(line_type1, '(a)') frame_number(1: numlen)//","
      else
        write(line_type2, '(a)') frame_number(1: numlen)//":"
      end if

    end if

    ! frame_no with given frame number is set to 1 
    nframe_save = 0
    if (len_trim(line_type1) /= 0) then
      istart = 1
      numlen = len_trim(line_type1)
      do
        iposf = index(line_type1(istart: numlen), ",") 
        if (iposf == 0) iposf = index(line_type1(istart: numlen), " ") 
        if (iposf == 0) exit
        iposf = iposf + istart - 1
        read(line_type1(istart: iposf - 1), *) iframe
        if (iframe > nframes) then
          call error_msg('Assign_Number> '//&
           'a frame number should be smaller than md_steps')
        end if
        frame_no(iframe) = 1 
        nframe_save      = nframe_save + 1
        istart           = iposf + 1
        if (istart > numlen) exit
      end do
    end if
    if (len_trim(line_type2) /= 0) then
      istart = 1
      numlen = len_trim(line_type2)
      iposf = index(line_type2(istart: numlen), ":") + istart - 1
      read(line_type2(istart: iposf - 1), *) ibegin
      istart = iposf + 1
      iposf = index(line_type2(istart: numlen), ":") + istart - 1
      read(line_type2(istart: iposf - 1), *) istride
      istart = iposf + 1
      iposf = index(line_type2(istart: numlen), ":") + istart - 1
      read(line_type2(istart: iposf - 1), *) iend

      do iframe = ibegin, iend, istride
        if (iend > nframes) exit
        if (iframe > nframes) then
          call error_msg('Assign_Number> '//&
           'a frame number should be smaller than md_steps')
        end if
        frame_no(iframe) = 1
        nframe_save      = nframe_save + 1
      end do
    end if

    return

  end subroutine assign_number


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    merge_selection
  !> @brief        merge multiple expression defined in selection
  !! @authors      KYMD
  !! @param[in]    groups         : expressions defined in selection
  !! @param[out]   exp_qmmxx      : results
  !! @param[in]    str_index_qmxx : indexes of groups assigned in option
  !! @param[in]    ngrp_sel       : # of expressions
  !! @param[in]    ngrp_qmxx      : # of indexes
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine merge_selection(groups, exp_qmxx, str_index_qmxx, &
                             ngrp_sel, ngrp_qmxx)

    ! formal arguments
    character(len=*),        intent(in)    :: groups(:)
    character(len=*),        intent(out)   :: exp_qmxx
    character(len=*),        intent(in)    :: str_index_qmxx
    integer,                 intent(in)    :: ngrp_sel, ngrp_qmxx

    ! local arguments
    integer                  :: igrp, igposf, iglen, igposl
    !integer, allocatable     :: list_index_qmxx(:)
    !integer                  :: alloc_stat, dealloc_stat
    integer                  :: list_index_qmxx(1: ngrp_qmxx)


    !allocate(list_index_qmxx(1: ngrp_qmxx), stat = alloc_stat)
    !if (alloc_stat /= 0) call error_msg_alloc

    call split(ngrp_qmxx, ngrp_qmxx, str_index_qmxx, list_index_qmxx)

    if (any(list_index_qmxx(1: ngrp_qmxx) > ngrp_sel)) &
        call error_msg(                                &
        'Merge_Selection> QM/MM or QM atom selection includes unknown groups.')

    igposl = 0
    do igrp = 1, ngrp_qmxx
      igposf = igposl + 1
      iglen  = len_trim(groups(list_index_qmxx(igrp)))
      igposl = igposl + iglen 
      if (igposl > MaxLine) igposl = MaxLine
      write(exp_qmxx(igposf: igposl), '(a)') trim(groups(list_index_qmxx(igrp)))

      if (igrp /= ngrp_qmxx) then
        igposf = igposf + iglen
        igposl = igposl + 1
        write(exp_qmxx(igposf:igposl), '(a)') '|'
      end if
    end do

    !deallocate(list_index_qmxx, stat = dealloc_stat)
    !if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine merge_selection

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_around
  !> @brief        check whether expression exp_qmxx includes around or not
  !! @authors      KYMD
  !! @param[in   ] exp_qmmxx      : expression of selection
  !! @param[inout] jcehck_around  : results
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_around(exp_qmxx, jcheck_around)
    
    ! formal arguments
    character(len=*),        intent(in)    :: exp_qmxx
    integer,                 intent(inout) :: jcheck_around

    ! local arguments
    character(len=MaxLine)   :: lowexp_qmxx
    integer                  :: i, jexist_around_atom, &
                                   jexist_around_res,  &
                                   jexist_around_mol
    integer                  :: len_exp


    len_exp = len_trim(exp_qmxx)
    write(lowexp_qmxx, '(a)') exp_qmxx
    do i = 1, len_exp
      if (lowexp_qmxx(i:i) >= 'A' .and. lowexp_qmxx(i:i) <= 'Z') then
        lowexp_qmxx(i:i) = char(ichar(lowexp_qmxx(i:i)) + 32)
      end if
    end do

    jexist_around_atom = index(lowexp_qmxx(1: len_exp), 'around:')   &
                         + index(lowexp_qmxx(1: len_exp), 'around ') &
                         + index(lowexp_qmxx(1: len_exp), 'around_atoms') 
    jexist_around_res  = index(lowexp_qmxx(1: len_exp), 'around_res')   &
                         + index(lowexp_qmxx(1: len_exp), 'around_residues') 
    jexist_around_mol  = index(lowexp_qmxx(1: len_exp), 'around_mol')   &
                         + index(lowexp_qmxx(1: len_exp), 'around_molecules') 

    if (jexist_around_atom /= 0) jcheck_around = jcheck_around + 1
    if (jexist_around_res  /= 0) jcheck_around = jcheck_around + 10
    if (jexist_around_mol  /= 0) jcheck_around = jcheck_around + 100

    return

  end subroutine check_around

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_wrap_qmmm
  !> @brief        setup wrap procedure (this routine is just wrapper)
  !! @authors      KYMD
  !! @param[inout] molecule : molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_wrap_qmmm(molecule)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule

    call setup_pbc_correct(PBCCModeMolecule, molecule)

    return

  end subroutine setup_wrap_qmmm

end module qg_option_mod
