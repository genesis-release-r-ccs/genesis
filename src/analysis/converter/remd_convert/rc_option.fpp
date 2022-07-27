!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   rc_option_mod
!> @brief   module for REMD trajectory conversion options
!! @authors Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module rc_option_mod

  use constants_mod
  use rc_option_str_mod
  use pbc_correct_mod
  use trajectory_str_mod
  use output_mod
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
    logical             :: check_only      = .false.
    logical             :: allow_backup    = .false.
    integer             :: convert_type    = ConvertTypeParameter
    character(MaxLine)  :: convert_ids     = ''
    integer             :: num_replicas    = 0
    integer             :: nsteps          = 0
    integer             :: exchange_period = 0
    integer             :: crdout_period   = 0
    integer             :: eneout_period   = 0
    integer             :: logout_period   = 0
    integer             :: trjout_period   = 0
    integer             :: trjout_format   = TrjFormatDCD
    integer             :: trjout_type     = TrjTypeCoorBox
    integer             :: trjout_atom     = 1
    integer             :: pbc_correct     = PBCCModeNo
    logical             :: centering       = .false.
    integer             :: centering_atom  = 1
    character(MaxLine)  :: center_coord    = '0.0 0.0 0.0'
  end type s_opt_info

  public  :: show_ctrl_option
  public  :: read_ctrl_option
  public  :: setup_option
  private :: parse_convert_ids

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
    write(MsgOut,'(A)') 'check_only      = NO             # only checking input files (YES/NO)'
    write(MsgOut,'(A)') 'allow_backup    = NO             # backup existing output files (YES/NO)'
    write(MsgOut,'(A)') 'convert_type    = PARAMETER      # (REPLICA/PARAMETER)'
    write(MsgOut,'(A)') 'convert_ids     =                # selected index (empty = all)(example: 1 2 5-10)'
    write(MsgOut,'(A)') 'num_replicas    = 0              # total number of replicas used in the simulation'
    write(MsgOut,'(A)') 'nsteps          = 0              # nsteps in [DYNAMICS]'
    write(MsgOut,'(A)') 'exchange_period = 0              # exchange_period in [REMD]'
    write(MsgOut,'(A)') 'crdout_period   = 0              # crdout_period in [DYNAMICS]'
    write(MsgOut,'(A)') 'eneout_period   = 0              # eneout_period in [DYNAMICS]'
    write(MsgOut,'(A)') 'logout_period   = 0              # output frequency of logfile'
    write(MsgOut,'(A)') 'trjout_period   = 0              # output frequency of trjfile'
    write(MsgOut,'(A)') 'trjout_format   = DCD            # (PDB/DCD)'
    write(MsgOut,'(A)') 'trjout_type     = COOR+BOX       # (COOR/COOR+BOX)'
    write(MsgOut,'(A)') 'trjout_atom     = 1              # atom group'
    write(MsgOut,'(A)') 'centering       = NO             # shift center of mass (YES requres psf/prmtop/grotop)'
    write(MsgOut,'(A)') 'centering_atom  = 1              # atom group'
    write(MsgOut,'(A)') 'center_coord    = 0.0 0.0 0.0    # target center coordinates'
    write(MsgOut,'(A)') 'pbc_correct     = NO             # (NO/MOLECULE)'
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



    ! read parameters
    !

    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_logical(handle, Section, 'check_only',      &
                               opt_info%check_only)

    call read_ctrlfile_logical(handle, Section, 'allow_backup',    &
                               opt_info%allow_backup)

    call read_ctrlfile_type   (handle, Section, 'convert_type',    &
                               opt_info%convert_type, ConvertTypeTypes)

    call read_ctrlfile_string (handle, Section, 'convert_ids',     &
                               opt_info%convert_ids)

    call read_ctrlfile_integer(handle, Section, 'num_replicas',    &
                               opt_info%num_replicas)

    call read_ctrlfile_integer(handle, Section, 'nsteps',          &
                               opt_info%nsteps)

    call read_ctrlfile_integer(handle, Section, 'exchange_period', &
                               opt_info%exchange_period)

    call read_ctrlfile_integer(handle, Section, 'crdout_period',   &
                               opt_info%crdout_period)

    call read_ctrlfile_integer(handle, Section, 'eneout_period',   &
                               opt_info%eneout_period)

    call read_ctrlfile_integer(handle, Section, 'logout_period',   &
                               opt_info%logout_period)

    call read_ctrlfile_integer(handle, Section, 'trjout_period',   &
                               opt_info%trjout_period)

    call read_ctrlfile_type   (handle, Section, 'trjout_format',   &
                               opt_info%trjout_format, TrjFormatTypes)

    call read_ctrlfile_type   (handle, Section, 'trjout_type',     &
                               opt_info%trjout_type, TrjTypeTypes)

    call read_ctrlfile_integer(handle, Section, 'trjout_atom',     &
                               opt_info%trjout_atom)

    call read_ctrlfile_logical(handle, Section, 'centering',       &
                               opt_info%centering)

    call read_ctrlfile_integer(handle, Section, 'centering_atom',  &
                               opt_info%centering_atom)

    call read_ctrlfile_string (handle, Section, 'center_coord',    &
                               opt_info%center_coord)

    call read_ctrlfile_type   (handle, Section, 'pbc_correct',     &
                               opt_info%pbc_correct, PBCCModeTypes)

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

    write(MsgOut,'(A20,A10)')  '  convert type    = ', &
         ConvertTypeTypes(opt_info%convert_type)
    write(MsgOut,'(A20,A10)')  '  convert ids     = ', &
         opt_info%convert_ids
    write(MsgOut,'(A20,I10)')  '  num_replicas    = ', &
         opt_info%num_replicas
    write(MsgOut,'(A20,I10)')  '  nsteps          = ', &
         opt_info%nsteps
    write(MsgOut,'(A20,I10)')  '  exchange_period = ', &
         opt_info%exchange_period
    write(MsgOut,'(A20,I10)')  '  crdout_period   = ', &
         opt_info%crdout_period
    write(MsgOut,'(A20,I10)')  '  eneout_period   = ', &
         opt_info%eneout_period
    write(MsgOut,'(A20,I10)')  '  logout_period   = ', &
         opt_info%logout_period
    write(MsgOut,'(A20,I10)')  '  trjout_period   = ', &
         opt_info%trjout_period
    write(MsgOut,'(A20,A10)')  '  trjout format   = ', &
         TrjFormatTypes(opt_info%trjout_format)
    write(MsgOut,'(A20,A10)')  '  trjout type     = ', &
         TrjTypeTypes(opt_info%trjout_type)
    write(MsgOut,'(A20,A5,I0)')'  trjout atom     = ', &
         'group', opt_info%trjout_atom

    if (opt_info%centering) then
      write(MsgOut,'(A20,A3)')    '  centering       = ', &
                                  'yes'
      write(MsgOut,'(A20,A5,I0)') '  centering atom  = ', &
                                  'group', opt_info%centering_atom
      write(MsgOut,'(A20,A)')     '  center_coord    = ', &
                                  trim(opt_info%center_coord)
    else
      write(MsgOut,'(A20,A2)')    '  centering       = ', &
                                  'no'
    end if

    write(MsgOut,'(A20,A10)')  '  pbc correction  = ', &
         PBCCModeTypes(opt_info%pbc_correct)
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
  !! @param[in]    out_info : OUTPUT section control parameters information
  !! @param[inout] molecule : molecule information
  !! @param[out]   option   : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine setup_option(opt_info, &
                          sel_info, &
                          out_info, &
                          molecule, &
                          option)
  
    ! formal argments
    type(s_opt_info),        intent(in)    :: opt_info
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_out_info),        intent(in)    :: out_info
    type(s_molecule),        intent(inout) :: molecule
    type(s_option),          intent(inout) :: option

    ! local variables
    integer                  :: ndata

    ! check only
    option%check_only      = opt_info%check_only

    ! backup output file
    call setup_backup(opt_info%allow_backup)

    ! convert type
    option%convert_type    = opt_info%convert_type

    ! convert ids
    call parse_convert_ids(opt_info%convert_ids, option%convert_ids)

    ! REMD information
    option%num_replicas    = opt_info%num_replicas
    option%nsteps          = opt_info%nsteps
    option%exchange_period = opt_info%exchange_period
    option%crdout_period   = opt_info%crdout_period
    option%eneout_period   = opt_info%eneout_period

    ! logfile output frequency
    if (opt_info%logout_period == 0) then
      option%logout_period   = opt_info%eneout_period
      write(MsgOut,'(A)') 'Setup_Option> logout_period is set to eneout_period.'
    else
      option%logout_period   = opt_info%logout_period
    end if

    if (out_info%trjfile /= '') then

      ! trajectory output frequency
      if (opt_info%trjout_period == 0) then
        option%trjout_period   = option%crdout_period
        write(MsgOut,'(A)') 'Setup_Option> trjout_period is set to crdout_period.'
      else
        option%trjout_period   = opt_info%trjout_period
      end if

      ! trajectory format
      option%trjout_format   = opt_info%trjout_format

      ! trajectory type
      option%trjout_type     = opt_info%trjout_type

      ! trajectory atom
      if (opt_info%trjout_atom > size(sel_info%groups)) &
        call error_msg('Setup_Option> trj-out atom selection is out of range.')

      option%trjout_atom_exp = sel_info%groups(opt_info%trjout_atom)

      call select_atom(molecule, &
                       option%trjout_atom_exp, &
                       option%trjout_atom)

      write(MsgOut,'(A,I8)') 'Setup_Option> trj-out atom count: ', &
           size(option%trjout_atom%idx)

    end if

    ! centering
    option%centering = opt_info%centering

    if (option%centering) then
      if (opt_info%centering_atom > size(sel_info%groups)) &
        call error_msg('Setup_Option> centering atom selection is out of range.')

      option%centering_atom_exp = sel_info%groups(opt_info%centering_atom)

      call select_atom(molecule, &
                       option%centering_atom_exp, &
                       option%centering_atom)

      write(MsgOut,'(A,I8)') 'Setup_Option> centering atom count: ', &
        size(option%centering_atom%idx)


      ndata = split_num(opt_info%center_coord)
      if (ndata == 3) then
        call split(ndata, 3, opt_info%center_coord, option%center_coord)
        write(MsgOut,'(A,3F10.3)') 'Setup_Option> center coordinates: ', &
                                   option%center_coord(1:3)
      else
        call error_msg('Setup_Option> center_coord in [OPTION] lacks X or Y or Z')
      endif

    end if

    ! pbc correct
    option%pbcc_mode = opt_info%pbc_correct

    ! setup for PBC-correct mode "molecule"
    call setup_pbc_correct(option%pbcc_mode, molecule)

    write(MsgOut,'(A)') ' '

    return

  end subroutine setup_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    parse_convert_ids
  !> @brief        parse convert ids
  !! @authors      NT
  !! @param[in]    str     : ids information
  !! @param[out]   id_list : id list
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine parse_convert_ids(str, id_list)
  
    ! parameters
    integer,      parameter  :: MaxTok = 100

    ! formal arguments
    character(*),            intent(in)    :: str
    integer,    allocatable, intent(inout) :: id_list(:)
   
    ! local variables
    integer                  :: i, j, ntok, nids, rids, alloc_stat
    integer                  :: toki(2,MaxTok), ifr, ito
    character(20)            :: toks(  MaxTok)

    logical,    allocatable  :: ids(:)


    if (len_trim(str) == 0) return

    ! parse ids string
    toks(:) = ''
    read(str,*,end=10,err=900) toks

10  ntok = 0
    rids = 0

    do i = 1, MaxTok

      if (len_trim(toks(i)) == 0) &
        exit

      j = index(toks(i),'-')
      if (j /= 0) then

        ! range id
        toks(i)(j:j) = ' '
        read(toks(i),*,end=900,err=900) ifr, ito
        if (ifr > ito) ifr = ito

      else

        ! single id
        read(toks(i),*,end=900,err=900) ito
        ifr = ito

      end if

      toki(1,i) = ifr
      toki(2,i) = ito
      ntok = ntok + 1
      rids = MAX(rids, ito)

    end do

    ! check id 
    if (ntok == 0) &
      call error_msg('Parse_Convert_Ids> There is no IDs.')

    alloc_stat = 0
    allocate(ids(rids), stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg('Parse_Convert_Ids> Id is too large.')
    ids(:) = .false.

    do i = 1, ntok
      do j = toki(1,i), toki(2,i)
        ids(j) = .true.
      end do
    end do

    nids = 0
    do i = 1, rids
      if (ids(i)) nids = nids + 1
    end do

    ! make id list
    allocate(id_list(nids), stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    nids = 0
    do i = 1, rids
      if (ids(i)) then
        nids = nids + 1
        id_list(nids) = i
      end if
    end do

    ! deallocate
    deallocate(ids, stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_dealloc

    return

900 call error_msg('Parse_Convert_Ids> Read error. :'//trim(str))

  end subroutine parse_convert_ids

end module rc_option_mod
