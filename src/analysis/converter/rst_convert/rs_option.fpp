!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   rs_option_mod
!> @brief   module for conversion options
!! @authors Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module rs_option_mod

  use rs_option_str_mod
  use output_str_mod
  use fileio_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use constants_mod
  use pbc_correct_mod
  use molecules_str_mod

  implicit none
  private

  ! structures
  type, public :: s_opt_info

    logical                         :: check_only        = .false.
    logical                         :: allow_backup      = .false.
    integer                         :: rstin_format      = RstFormatGENESIS
    integer                         :: rstout_format     = RstFormatPDB
    integer                         :: rstout_type       = RstTypeAuto
    integer                         :: pbc_correct       = PBCCModeNO
    integer                         :: nreplicas         = -1
    character(MaxLine)              :: shift_coord       = ''
    character(MaxLine), allocatable :: rename_res(:)

    real(wp)                        :: min_energy        = 0.0_wp
    real(wp)                        :: min_delta_r       = 0.0_wp

    integer                         :: md_iseed          = 0
    integer                         :: step              = 0
    real(wp)                        :: md_thermo_moment  = 0.0_wp
    real(wp)                        :: md_baro_moment(3) = 0.0_wp
    real(wp)                        :: origin_x          = 0.0_wp
    real(wp)                        :: origin_y          = 0.0_wp
    real(wp)                        :: origin_z          = 0.0_wp

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
    write(MsgOut,'(A)') 'rstin_format   = GENESIS         # (NAMD/GENESIS)'
    write(MsgOut,'(A)') 'rstout_format  = PDB             # (NAMD/GENESIS/PDB/CHARMM/AMBER)'
    write(MsgOut,'(A)') 'rstout_type    = AUTO            # (AUTO/MIN/MD/REMD)'
    write(MsgOut,'(A)') 'pbc_correct    = NO              # (NO/MOLECULE)'
    write(MsgOut,'(A)') '# rename_res1    = HSE HIS       # change residue name'
    write(MsgOut,'(A)') '# rename_res2    = HSD HIS       # change residue name'
    write(MsgOut,'(A)') '# shift_coord    = 0.0 0.0 0.0   # These are added to the original coordinates'
    write(MsgOut,'(A)') '# nreplicas      = 16            # If multiple restart files (input{}.rst) generated from'
    write(MsgOut,'(A)') '                                 # REMD are converted, specify the number of replicas'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '# The following are for NAMD'
    write(MsgOut,'(A)') 'min_energy     = 0.0             # energy'
    write(MsgOut,'(A)') 'min_delta_r    = 0.0             # delta-r'
    write(MsgOut,'(A)') 'md_iseed       = 0               # random seed'
    write(MsgOut,'(A)') 'md_thermo_moment = 0.0           # thermostat momentum'
    write(MsgOut,'(A)') 'md_baro_moment   = 0.0 0.0 0.0   # barostat momentum'
    write(MsgOut,'(A)') 'namd_step        = 0             # step'
    write(MsgOut,'(A)') 'namd_originx     = 0.0           # box origin (x)'
    write(MsgOut,'(A)') 'namd_originy     = 0.0           # box origin (y)'
    write(MsgOut,'(A)') 'namd_originz     = 0.0           # box origin (z)'
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

    integer                                :: i, nresi
    character(Maxline)                     :: value, line
    character(MaxLine)                     :: rename_res


    ! read parameters
    !

    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_logical(handle, Section, 'check_only', &
                               opt_info%check_only)

    call read_ctrlfile_logical(handle, Section, 'allow_backup', &
                               opt_info%allow_backup)

    call read_ctrlfile_type   (handle, Section, 'rstin_format', &
                               opt_info%rstin_format, RstFormatTypes)

    call read_ctrlfile_type   (handle, Section, 'rstout_format', &
                               opt_info%rstout_format, RstFormatTypes)

    call read_ctrlfile_type   (handle, Section, 'rstout_type', &
                               opt_info%rstout_type, RstTypeTypes)

    call read_ctrlfile_type   (handle, Section, 'pbc_correct',    &
                               opt_info%pbc_correct, PBCCModeTypes)

    call read_ctrlfile_string (handle, Section, 'shift_coord', &
                               opt_info%shift_coord)

    call read_ctrlfile_integer(handle, Section, 'nreplicas', &
                               opt_info%nreplicas)


!     if (opt_info%rstout_format == opt_info%rstin_format)  then
!       call error_msg                                                         &
!            ('Read_Ctrl_Option> Formats of input/output should be different')
!     endif

    ! for rename residues
    !
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

    ! for NAMD
    !
    call read_ctrlfile_real   (handle, Section, 'min_energy', &
                               opt_info%min_energy)

    call read_ctrlfile_real   (handle, Section, 'min_delta_r', &
                               opt_info%min_delta_r)

    call read_ctrlfile_integer(handle, Section, 'md_iseed' , &
                               opt_info%md_iseed)

    call read_ctrlfile_real   (handle, Section, 'md_thermo_moment', &
                               opt_info%md_thermo_moment)

    if (opt_info%rstout_format == RstFormatGENESIS) then
      line = ''
      call read_ctrlfile_string(handle, Section, 'md_baro_moment', line)
  
      if (line /= '') &
        read(line,*,err=900,end=900) opt_info%md_baro_moment
    endif

    call read_ctrlfile_integer(handle, Section, 'namd_step', &
                               opt_info%step)

    call read_ctrlfile_real   (handle, Section, 'namd_originX', &
                               opt_info%origin_x)

    call read_ctrlfile_real   (handle, Section, 'namd_originY', &
                               opt_info%origin_y)

    call read_ctrlfile_real   (handle, Section, 'namd_originZ', &
                               opt_info%origin_z)

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

    write(MsgOut,'(A20,A10)')   '  rst in format   = ', &
         RstFormatTypes(opt_info%rstin_format)

    write(MsgOut,'(A20,A10)')   '  rst out format  = ', &
         RstFormatTypes(opt_info%rstout_format)

    write(MsgOut,'(A20,A10)')   '  rst out type    = ', &
         RstTypeTypes(opt_info%rstout_type)

    write(MsgOut,'(A20,I10)')     &
                 '  # of replicas   = ', opt_info%nreplicas

    write(MsgOut,'(A20,A10)')     &
                 '  pbc correction  = ', &
                 PBCCModeTypes(opt_info%pbc_correct)

    write(MsgOut,'(A20,A)') &
                 ' shift coordinate = ', trim(opt_info%shift_coord)

    write(MsgOut,'(A20,I10)') '  # of rename res = ', nresi
    do i = 1, nresi
      write(MsgOut,'(A10,I0,A3,A)')      &
                 '   rename ', i, ' = ', trim(opt_info%rename_res(i))
    end do

    ! for NAMD
    !
    write(MsgOut,'(A20,F10.3)') '  min: energy     = ', &
         opt_info%min_energy

    write(MsgOut,'(A20,F10.3)') '  min: delta-r    = ', &
         opt_info%min_delta_r

    write(MsgOut,'(A20,I16)')   '  md: iseed       = ', &
         opt_info%md_iseed

    write(MsgOut,'(A20,F16.9)') '  md: thermo moment = ', &
         opt_info%md_thermo_moment

    write(MsgOut,'(A20,F16.9)') '  md: baro moment   = ', &
         opt_info%md_baro_moment(1)

    write(MsgOut,'(A20,F16.9)') '                    ', &
         opt_info%md_baro_moment(2)

    write(MsgOut,'(A20,F16.9)') '                    ', &
         opt_info%md_baro_moment(3)

    write(MsgOut,'(A20,I16)')   'namd: step          = ', &
         opt_info%step

    write(MsgOut,'(A20,F16.9)') 'namd: origin_x      = ', &
         opt_info%origin_x

    write(MsgOut,'(A20,F16.9)') 'namd: origin_y      = ', &
         opt_info%origin_y

    write(MsgOut,'(A20,F16.9)') 'namd: origin_z      = ', &
         opt_info%origin_z

    write(MsgOut,'(A)') ' '

    return

900 call error_msg('Read_Ctrl_Option> format error : "md_baro_moment"')

  end subroutine read_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_option
  !> @brief        setup option information
  !! @authors      NT
  !! @param[in]    opt_info : OPTION section control parameters information
  !! @param[out]   option   : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine setup_option(opt_info, molecule, option)
  
    ! formal argments
    type(s_opt_info),        intent(inout) :: opt_info
    type(s_molecule),        intent(inout) :: molecule
    type(s_option),          intent(inout) :: option

    ! local variable
    integer                  :: i, j
    integer                  :: nresi, nrenm, ndata
    character(MaxLine)       :: s1, s2
    character(6)             :: st


    ! check only
    option%check_only    = opt_info%check_only

    ! backup output file
    call setup_backup(opt_info%allow_backup)

    ! restart input format
    option%rstin_format = opt_info%rstin_format

    ! restart output format
    option%rstout_format = opt_info%rstout_format

    ! restart output type
    option%rstout_type = opt_info%rstout_type

    ! restart output type
    option%nreplicas   = opt_info%nreplicas

    ! pbc correct
    option%pbcc_mode = opt_info%pbc_correct

    ! shift coordinates
    if (opt_info%shift_coord /= '') then
      ndata = split_num(opt_info%shift_coord)
      if (ndata == 3) then
        call split(ndata, 3, opt_info%shift_coord, option%shift_coord)
        write(MsgOut,'(A,3F10.3)') 'Setup_Option> shift coordinates: ', &
                                   option%shift_coord(1:3)
      else
        call error_msg('Setup_Option> shift_corod in [OPTION] lacks X or Y or Z')
      end if

    else
      option%shift_coord(1:3) = 0.0_wp
    end if

    ! setup for PBC-correct mode "molecule"
    call setup_pbc_correct(option%pbcc_mode, molecule)

    ! minimization parameters
    option%min_energy  = opt_info%min_energy
    option%min_delta_r = opt_info%min_delta_r

    ! MD parameters
    option%md_iseed         = opt_info%md_iseed
    option%md_thermo_moment = opt_info%md_thermo_moment
    option%md_baro_moment   = opt_info%md_baro_moment

    ! NAMD parameters
    option%step             = opt_info%step
    option%origin_x         = opt_info%origin_x
    option%origin_y         = opt_info%origin_y
    option%origin_z         = opt_info%origin_z

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

    deallocate(opt_info%rename_res)

    write(MsgOut,'(A)') ' '

    return

900 call error_msg('Setup_Option> Rename residue : Bad format:'//&
                   trim(opt_info%rename_res(i)))

  end subroutine setup_option

end module rs_option_mod
