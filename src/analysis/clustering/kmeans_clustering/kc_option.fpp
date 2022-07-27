!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   kc_option_mod
!> @brief   module for analysis options
!! @authors Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module kc_option_mod

  use kc_option_str_mod
  use select_mod
  use select_atoms_mod
  use molecules_str_mod
  use fileio_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use constants_mod
  use trajectory_str_mod
  use pbc_correct_mod

  implicit none
  private

  ! structures
  type, public :: s_opt_info
    logical                         :: check_only     = .false.
    logical                         :: allow_backup   = .false.
    integer                         :: analysis_atom  = 1
    integer                         :: num_clusters   = 10
    integer                         :: max_iteration  = 100
    real(wp)                        :: stop_threshold = 98.0
    integer                         :: num_iterations = 10
    integer                         :: iseed          = 3141592
    integer                         :: trjout_format  = TrjFormatPDB
    integer                         :: trjout_type    = TrjTypeCoor
    integer                         :: trjout_atom    = 1
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
    write(MsgOut,'(A)') 'check_only      = NO             # only checking input files (YES/NO)'
    write(MsgOut,'(A)') 'allow_backup    = NO             # backup existing output files (YES/NO)'
    write(MsgOut,'(A)') 'analysis_atom   = 1              # target atoms for the cluster analysis'
    write(MsgOut,'(A)') 'num_clusters    = 10             # number of clusters'
    write(MsgOut,'(A)') 'max_iteration   = 100            # max number of iterations for k-means algorithm'
    write(MsgOut,'(A)') 'stop_threshold  = 98.0           # stop threshold of convergence (%) for k-means algorithm'
    write(MsgOut,'(A)') 'num_iterations  = 10             # number of iterations for calculating averaged coordinates'
    write(MsgOut,'(A)') 'trjout_atom     = 1              # atom selection for pdbfile and trjfile'
    write(MsgOut,'(A)') 'trjout_format   = DCD            # (PDB/DCD)'
    write(MsgOut,'(A)') 'trjout_type     = COOR           # (COOR/COOR+BOX)'
    write(MsgOut,'(A)') 'iseed           = 3141592        # random number seed'
    write(MsgOut,'(A)') ' '

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


    ! read parameters
    !

    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_logical(handle, Section, 'check_only', &
                               opt_info%check_only)

    call read_ctrlfile_logical(handle, Section, 'allow_backup', &
                               opt_info%allow_backup)

    call read_ctrlfile_integer(handle, Section, 'analysis_atom', &
                               opt_info%analysis_atom)

    call read_ctrlfile_integer(handle, Section, 'num_clusters', &
                               opt_info%num_clusters)

    call read_ctrlfile_integer(handle, Section, 'max_iteration', &
                               opt_info%max_iteration)

    call read_ctrlfile_real(handle, Section, 'stop_threshold', &
                               opt_info%stop_threshold)

    call read_ctrlfile_integer(handle, Section, 'num_iterations', &
                               opt_info%num_iterations)

    call read_ctrlfile_integer(handle, Section, 'iseed', &
                               opt_info%iseed)

    call read_ctrlfile_type   (handle, Section, 'trjout_format', &
                               opt_info%trjout_format, TrjFormatTypes)

    call read_ctrlfile_type   (handle, Section, 'trjout_type',  &
                               opt_info%trjout_type, TrjTypeTypes)

    call read_ctrlfile_integer(handle, Section, 'trjout_atom', &
                               opt_info%trjout_atom)

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

    write(MsgOut,'(A20,I0)') &
                    '  num_clusters    = ', opt_info%num_clusters

    write(MsgOut,'(A20,I0)') &
                    '  max_iteration   = ', opt_info%max_iteration

    write(MsgOut,'(A20,F10.3)') &
                    '  stop_threshold  = ', opt_info%stop_threshold

    write(MsgOut,'(A20,I0)') &
                    '  num interations = ', opt_info%num_iterations

    write(MsgOut,'(A20,I0)') &
                    '  iseed           = ', opt_info%iseed

    write(MsgOut,'(A20,A10)') &
                    '  trjout format   = ', TrjFormatTypes(opt_info%trjout_format)

    write(MsgOut,'(A20,A10)') &
                    '  trjout type     = ', TrjTypeTypes(opt_info%trjout_type)

    write(MsgOut,'(A20,A5,I0)') &
                    '  trjout atom     = ', 'group', opt_info%trjout_atom

    write(MsgOut,'(A)') ''

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


    ! check only
    option%check_only = opt_info%check_only

    ! backup output file
    call setup_backup(opt_info%allow_backup)

    ! number of iterations
    option%num_iterations = opt_info%num_iterations
    option%num_clusters   = opt_info%num_clusters
    option%stop_threshold = opt_info%stop_threshold
    option%max_iteration  = opt_info%max_iteration
    option%iseed          = opt_info%iseed
    option%trjout_format  = opt_info%trjout_format
    option%trjout_type    = opt_info%trjout_type

    ! analysis target atom
    if (opt_info%analysis_atom > size(sel_info%groups)) &
        call error_msg('Setup_Option> analysis target atom selection is out of range.')

    call select_atom(molecule, &
                     sel_info%groups(opt_info%analysis_atom), &
                     option%analysis_atom)

    write(MsgOut,'(A,I8)') 'Setup_Option> analysis atom count: ', &
         size(option%analysis_atom%idx)


    ! trjout target atom
    if (opt_info%trjout_atom > size(sel_info%groups)) &
        call error_msg('Setup_Option> trjout target atom selection is out of range.')

    call select_atom(molecule, &
                     sel_info%groups(opt_info%trjout_atom), &
                     option%trjout_atom)

    write(MsgOut,'(A,I8)') 'Setup_Option> output atom count: ', &
         size(option%trjout_atom%idx)


    write(MsgOut,'(A)') ''

    return

  end subroutine setup_option

end module kc_option_mod
