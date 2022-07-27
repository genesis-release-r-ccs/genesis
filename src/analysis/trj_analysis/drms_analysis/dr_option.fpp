!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   dr_option_mod
!> @brief   module for analysis options
!! @authors Chigusa Kobayashi (CK), Daisuke Matsuoka (DM), Norio Takase (NT)
! 
!  (c) Copyright 2018 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module dr_option_mod

  use dr_option_str_mod
  use select_mod
  use select_atoms_mod
  use select_atoms_str_mod
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
    logical                     :: check_only           = .true.
    logical                     :: verbose              = .true.
    logical                     :: allow_backup         = .false.
    logical                     :: avoid_bonding        = .true.
    logical                     :: pbc_correct          = .false.
    logical                     :: pbc_correct_setup    = .false.
    logical                     :: ignore_hydrogen      = .false.
    logical                     :: two_states           = .false.
    integer                     :: exclude_residues     = 0
    real(wp)                    :: minimum_distance     = 1.0_wp
    real(wp)                    :: maximum_distance     = 10.0_wp
    real(wp)                    :: minimum_difference   = 5.0_wp
    real(wp)                    :: box_size_ref_x       = 0.0_wp
    real(wp)                    :: box_size_ref_y       = 0.0_wp
    real(wp)                    :: box_size_ref_z       = 0.0_wp
    real(wp)                    :: box_size_cur_x       = 0.0_wp
    real(wp)                    :: box_size_cur_y       = 0.0_wp
    real(wp)                    :: box_size_cur_z       = 0.0_wp
    character(Maxline)          :: contact_groups       = ''
    character(Maxline)          :: exclude_groups       = ''
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
  !! @authors      DM, NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_option


    write(MsgOut,'(A)') '[SELECTION]'
    write(MsgOut,'(A)') 'group1            = sid:PROA and not hydrogen   # atom group 1'
    write(MsgOut,'(A)') '# group2          = sid:PROB and not hydrogen   # atom group 2'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '[OPTION]'
    write(MsgOut,'(A)') 'check_only             = YES             # (YES/NO)'
    write(MsgOut,'(A)') 'allow_backup           = NO        # backup existing output files (YES/NO)'
    write(MsgOut,'(A)') ' '
    write(MsgOut,'(A)') 'pbc_correct            = NO      '
    write(MsgOut,'(A)') 'ignore_hydrogen        = NO      '
    write(MsgOut,'(A)') 'two_states             = NO      '
    write(MsgOut,'(A)') 'avoid_bonding          = YES     '
    write(MsgOut,'(A)') 'exclude_residues       = 0       '
    write(MsgOut,'(A)') 'minimum_distance       = 1.0     '
    write(MsgOut,'(A)') 'maximum_distance       = 10.0    '
    write(MsgOut,'(A)') 'minimum_difference     = 5.0     '
    write(MsgOut,'(A)') 'contact_groups         = 1 2     '
    write(MsgOut,'(A)') 'exclude_groups         = 3       '
    write(MsgOut,'(A)') ' '

    return

  end subroutine show_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_option
  !> @brief        read control parameters in OPTION section
  !! @authors      CK, DM, NT
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
    integer           :: i
    character(len=20) :: name


    ! read parameters
    !
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_logical(handle, Section, &
                              'check_only', opt_info%check_only)

    call read_ctrlfile_logical(handle, Section, &
                              'allow_backup', opt_info%allow_backup)

    call read_ctrlfile_logical(handle, Section, &
                              'avoid_bonding', opt_info%avoid_bonding)

    call read_ctrlfile_logical(handle, Section, &
                              'ignore_hydrogen', opt_info%ignore_hydrogen)

    call read_ctrlfile_logical(handle, Section, &
                              'pbc_correct', opt_info%pbc_correct)

    call read_ctrlfile_logical(handle, Section, &
                              'two_states', opt_info%two_states)

    call read_ctrlfile_integer(handle, Section, 'exclude_residues', &
                               opt_info%exclude_residues)

    call read_ctrlfile_real   (handle, Section, 'minimum_distance', &
                               opt_info%minimum_distance)

    call read_ctrlfile_real   (handle, Section, 'maximum_distance', &
                               opt_info%maximum_distance)

    call read_ctrlfile_real   (handle, Section, 'minimum_difference', &
                               opt_info%minimum_difference)

    call read_ctrlfile_real   (handle, Section, 'box_size_x', &
                               opt_info%box_size_ref_x)

    call read_ctrlfile_real   (handle, Section, 'box_size_y', &
                               opt_info%box_size_ref_y)

    call read_ctrlfile_real   (handle, Section, 'box_size_z', &
                               opt_info%box_size_ref_z)

    call read_ctrlfile_real   (handle, Section, 'box_size_cur_x', &
                               opt_info%box_size_cur_x)

    call read_ctrlfile_real   (handle, Section, 'box_size_cur_y', &
                               opt_info%box_size_cur_y)

    call read_ctrlfile_real   (handle, Section, 'box_size_cur_z', &
                               opt_info%box_size_cur_z)

    call read_ctrlfile_string (handle, Section, 'contact_groups', &
                               opt_info%contact_groups)

    call read_ctrlfile_string (handle, Section, 'exclude_groups', &
                               opt_info%exclude_groups)


    call read_ctrlfile_logical(handle, Section, 'verbose', opt_info%verbose)

    call end_ctrlfile_section(handle)

    ! Error check
    if (opt_info%minimum_distance <= 0.0_wp) then
      call error_msg('Read_Ctrl_Option> ERROR : minimum_distance was not set correctly.')
    endif

    if (opt_info%maximum_distance <= 0.0_wp .or.  &
        opt_info%maximum_distance < opt_info%minimum_distance ) then
      call error_msg('Read_Ctrl_Option> ERROR : maximum_distance was not set correctly.')
    endif

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

    write(MsgOut,'(A20,F10.3,A20,F10.3)') &
     ' minimum_distance = ', opt_info%minimum_distance, &
     ' maximum_distance = ', opt_info%maximum_distance

    write(MsgOut,'(A20,I5)') &
         ' exclude_residues = ', opt_info%exclude_residues

    if (opt_info%avoid_bonding) then
      write(MsgOut,'(A20,A3)') ' avoid_bonding    = ', 'yes'
    else
      write(MsgOut,'(A20,A2)') ' avoid_bonding    = ', 'no'
    end if

    if (opt_info%ignore_hydrogen) then
      write(MsgOut,'(A20,A3)') ' ignore_hydrogen  = ', 'yes'
    else
      write(MsgOut,'(A20,A2)') ' ignore_hydrogen  = ', 'no'
    end if

    if (opt_info%pbc_correct) then
      write(MsgOut,'(A20,A3)') ' pbc_correct      = ', 'yes'
    else                                          
      write(MsgOut,'(A20,A2)') ' pbc_correct      = ', 'no'
    end if

    if (opt_info%two_states) then
      write(MsgOut,'(A20,A3)') ' two_states       = ', 'yes'
      write(MsgOut,'(A20,F10.3)') &
       'minimum_difference= ', opt_info%minimum_difference
    else                                       
      write(MsgOut,'(A20,A2)') ' two_states       = ', 'no'
    end if

    if (opt_info%check_only) then
      write(MsgOut,'(A20,A3)') ' check only       = ', 'yes'
    else
      write(MsgOut,'(A20,A2)') ' check only       = ', 'no'
    end if

    if (len_trim(opt_info%contact_groups) <= 0) then
      call error_msg('Read_Ctrl_Option> ERROR : contact_groups should be defined.')
    else
      write(MsgOut,'(A20,A)') &
       ' contact_groups  = ', trim(opt_info%contact_groups)
    endif
    if (len_trim(opt_info%exclude_groups) > 0) then
      write(MsgOut,'(A20,A)') &
       ' exclude_groups  = ', trim(opt_info%exclude_groups)
    endif

    if (opt_info%pbc_correct .and. opt_info%box_size_ref_x > 0.0_wp .and. &
        opt_info%box_size_ref_y > 0.0_wp .and.  &
        opt_info%box_size_ref_z > 0.0_wp  ) then
      opt_info%pbc_correct_setup = .true.
      write(MsgOut,'(A20,3F10.3)') '  box_size(x,y,z) = ',     &
                                      opt_info%box_size_ref_x, &
                                      opt_info%box_size_ref_y, &
                                      opt_info%box_size_ref_z

      if (opt_info%two_states) then
         if (opt_info%box_size_cur_x > 0.0_wp .and. &
             opt_info%box_size_cur_y > 0.0_wp .and.  &
             opt_info%box_size_cur_z > 0.0_wp  ) then
           write(MsgOut,'(A24,3F10.3)') ' box_size_cur_(x,y,z) = ',  &
                                      opt_info%box_size_cur_x, &
                                      opt_info%box_size_cur_y, &
                                      opt_info%box_size_cur_z
         else
           call error_msg('Read_Ctrl_Option> ERROR : box_size_ref should be defined.')
         endif
      endif
    endif

    if (opt_info%pbc_correct .and. .not.  opt_info%pbc_correct_setup) then
      write(MsgOut,'(A)') 'WARNING: pbc_correct does not apply to initial structures' 
    endif

    write(MsgOut,'(A)') ' '

    if (opt_info%verbose) then
      write(MsgOut,'(A20,A3)') '  verbose         = ', 'yes'
    else
      write(MsgOut,'(A20,A2)') '  verbose         = ', 'no'
    end if


    write(MsgOut,'(A)') ' '

    return

  end subroutine read_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_option
  !> @brief        setup option information
  !! @authors      CK, DM, NT
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
    integer                         :: igroup(2)
    integer, allocatable            :: ex_group(:)
    integer                         :: ngroup, i, j, jatm
    integer                         :: icount
    integer                         :: alloc_stat
    type(s_selatoms)                :: selatoms_grp

    alloc_stat = 0

    ! check only
    !
    option%check_only = opt_info%check_only

    ! backup output file
    !
    call setup_backup(opt_info%allow_backup)

    ! verbose log output
    !
    option%verbose    = opt_info%verbose

    ! cutoff distance
    !
    option%minimum_distance       = opt_info%minimum_distance
    option%maximum_distance       = opt_info%maximum_distance
    option%minimum_difference     = opt_info%minimum_difference

    ! calc mode
    !
    option%avoid_bonding          = opt_info%avoid_bonding
    option%pbc_correct            = opt_info%pbc_correct
    option%pbc_correct_setup      = opt_info%pbc_correct_setup
    option%ignore_hydrogen        = opt_info%ignore_hydrogen
    option%exclude_residues       = opt_info%exclude_residues
    option%two_states             = opt_info%two_states   
    option%box_size_cur(1)        = opt_info%box_size_cur_x
    option%box_size_cur(2)        = opt_info%box_size_cur_y
    option%box_size_cur(3)        = opt_info%box_size_cur_z
    option%box_size_ref(1)        = opt_info%box_size_ref_x
    option%box_size_ref(2)        = opt_info%box_size_ref_y
    option%box_size_ref(3)        = opt_info%box_size_ref_z

    if (option%ignore_hydrogen) then
      if (len_trim(molecule%atom_cls_name(1)) == 0) &
      call error_msg('Setup_Option> atom class name is required in ignore_hydrogen')
      write(MsgOut, '(A)') "Setup_Option> ignore_hydrogen option is set"
    endif

    call alloc_option(option, DA_Atoms, molecule%num_atoms)

    ngroup = split_num(opt_info%contact_groups)
    option%identical_group = .false.
    if (ngroup == 1)  then
      option%identical_group = .true.
    else if (ngroup /= 2)  then
      call error_msg('Setup_Option> ERROR: # of groups must be 2.')
    endif

    call split(ngroup, ngroup, opt_info%contact_groups, igroup(1:ngroup))

    do i = 1, ngroup
      call select_atom(molecule, sel_info%groups(igroup(i)), selatoms_grp)
      option%num_atoms_group(i)=size(selatoms_grp%idx)
      icount=0
      do j = 1, option%num_atoms_group(i)
        jatm=selatoms_grp%idx(j)
        if (option%ignore_hydrogen) then
            if (molecule%atom_cls_name(i)(1:1) .eq. "H" .or. &
                molecule%atom_cls_name(i)(1:1) .eq. "D" .or. &
                molecule%atom_cls_name(i)(1:1) .eq. "h" .or. &
                molecule%atom_cls_name(i)(1:1) .eq. "d") &
              cycle
        endif
        icount = icount+1
        option%contact_atoms(icount,i)=jatm
      end do

      option%num_atoms_group(i) = icount
      write(MsgOut,'(a,i12,a,a)') 'Contact groups: ', i, ' :', &
                                trim(sel_info%groups(igroup(i)))
      write(MsgOut,'(a,i12)')     '    # of selected atom = ',  &
                                  option%num_atoms_group(i) 
      write(MsgOut,'(a)')       ''
    end do

    call alloc_option(option, DA_Exclude, molecule%num_atoms)

    if (len_trim(opt_info%exclude_groups) > 0) then
      ngroup = split_num(opt_info%exclude_groups)
     
      allocate(ex_group(1:ngroup), &
               stat = alloc_stat)
     
      if (alloc_stat /= 0)   call error_msg_alloc
     
      call split(ngroup, ngroup, opt_info%exclude_groups, ex_group(1:ngroup))
      do i = 1, ngroup
        call select_atom(molecule, sel_info%groups(ex_group(i)), selatoms_grp)
     
        do j = 1, size(selatoms_grp%idx(:))
          jatm=selatoms_grp%idx(j)
          if (option%exclude_group_list(jatm) /= 0)  then
            call error_msg('Setup_Option> ERROR : duplicate exclude_group')
          else
            option%exclude_group_list(jatm) = i
          endif
        end do
      
        write(MsgOut,'(a,i12,a,a)') 'Exclude groups: ', i, ' :', &
                                  trim(sel_info%groups(ex_group(i)))
        write(MsgOut,'(a,i12)')     '    # of selected atom = ', &
                                    size(selatoms_grp%idx(:))
        write(MsgOut,'(a)')       ''
      end do
    else
      write(MsgOut,'(a)') 'No exclude group is found'
    endif

    write(MsgOut,'(A)') ''

    return

  end subroutine setup_option

end module dr_option_mod
