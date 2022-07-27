!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ma_option_mod
!> @brief   module for analysis options
!! @authors Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ma_option_mod

  use ma_option_str_mod
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
    logical                         :: check_only      = .false.
    logical                         :: read_ref_pdb    = .false.
    logical                         :: read_ref_path   = .false.

    ! mbar variables
    integer                         :: num_replicas     = 0
    integer                         :: nreplica         = 0
    integer                         :: dimension        = 1
    integer                         :: input_type       = InputTypeUS
    integer                         :: nblocks          = 1
    integer                         :: self_iteration   = 5
    integer                         :: newton_iteration = 100
    character(MaxLine)              :: temperature
    real(wp)                        :: target_temperature = 0.0_wp
    real(wp)                        :: tolerance        = 10E-08_wp
    character(MaxLineLong_CV), allocatable :: rest_func_no(:)
    character(MaxLine), allocatable :: grids(:)
    character(MaxLine)              :: out_unit         = 'NONE'

    ! restraints variables
    integer,            allocatable :: rest_funcs(:)
    character(MaxLine), allocatable :: rest_sel_index(:)
    character(MaxLine), allocatable :: rest_constants(:)
    character(MaxLine), allocatable :: rest_references(:)
    logical,            allocatable :: is_periodic(:)
    real(wp),           allocatable :: box_size(:)

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

    write(MsgOut,'(A)') '[MBAR]'
    write(MsgOut,'(A)') 'nreplica         = 4'
    write(MsgOut,'(A)') 'dimension        = 1'
    write(MsgOut,'(A)') 'nblocks          = 1'
    write(MsgOut,'(A)') 'input_type       = US  # US REMD FEP REST MBGO'
    write(MsgOut,'(A)') 'self_iteration   = 5'
    write(MsgOut,'(A)') 'newton_iteration = 100'
    write(MsgOut,'(A)') 'temperature      = 300.0'
    write(MsgOut,'(A)') 'target_temperature = 300.0'
    write(MsgOut,'(A)') 'tolerance        = 10E-08'
    write(MsgOut,'(A)') 'rest_function1   = 1'
    write(MsgOut,'(A)') 'grids1           = 0.0 360.0 91 # min max num_grids, num_grids means the number of bins + 1'
    write(MsgOut,'(A)') 'output_unit      = NONE  # unit for fene and pmf [NONE, kcal/mol]. NONE corresponds to dimensionless'
    write(MsgOut,'(A)') ''
    ! write(MsgOut,'(A)') '[SELECTION]'
    ! write(MsgOut,'(A)') 'group1         = an:NL'
    ! write(MsgOut,'(A)') 'group2         = an:CA'
    ! write(MsgOut,'(A)') 'group3         = an:CRP'
    ! write(MsgOut,'(A)') 'group4         = an:NR'
    ! write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '[RESTRAINTS]'
    write(MsgOut,'(A)') 'function1      = DIHED'
    write(MsgOut,'(A)') 'select_index1  = 1 3 7 11'
    write(MsgOut,'(A)') 'constant1      = 1.0 1.0 1.0 1.0'
    write(MsgOut,'(A)') 'reference1     = 144.0 146.0 148.0 150.0'
    write(MsgOut,'(A)') 'is_periodic1   = YES'
    write(MsgOut,'(A)') 'box_size1      = 360.0'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '# Please, be careful about the unit of the force constant and reference.'
    write(MsgOut,'(A)') '# For example, unit of the angle used in the control file of ATDYN and SPDYN is "degree".'
    write(MsgOut,'(A)') '# However, it may have to be "radian" in this control file, which depends on the input data.'
    write(MsgOut,'(A)') ''

    return

  end subroutine show_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_option
  !> @brief        read control parameters in MBAR/RESTRAINTS section
  !! @authors      NT
  !! @param[in]    handle   : unit number of control file
  !! @param[out]   opt_info : OPTION section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine read_ctrl_option(handle, opt_info)
  
    ! parameters
    character(*),            parameter     :: SectionMbar = 'MBAR'
    character(*),            parameter     :: SectionRest = 'RESTRAINTS'

    ! formal argments
    integer,                 intent(in)    :: handle
    type(s_opt_info),        intent(inout) :: opt_info

    ! local variables
    integer                  :: i, nfunc, len_str
    character(MaxLineLong_CV)       :: value
    character(30)            :: cdim, cfun


    ! read control parameters
    !


    ! MBAR section
    call begin_ctrlfile_section(handle, SectionMbar)

    call read_ctrlfile_logical(handle, SectionMbar, &
                               'check_only', opt_info%check_only)

    call read_ctrlfile_logical(handle, SectionMbar, &
                               'read_ref_pdb', opt_info%read_ref_pdb)

    call read_ctrlfile_logical(handle, SectionMbar, &
                               'read_ref_path', opt_info%read_ref_path)

    call read_ctrlfile_integer(handle, SectionMbar, &
                               'num_replicas', opt_info%num_replicas)

    call read_ctrlfile_integer(handle, SectionMbar, &
                               'nreplica', opt_info%nreplica)

    if (opt_info%num_replicas == 0 .and. opt_info%nreplica /= 0) &
      opt_info%num_replicas = opt_info%nreplica
    
    call read_ctrlfile_integer(handle, SectionMbar, &
                               'dimension', opt_info%dimension)

    call read_ctrlfile_type   (handle, SectionMbar, &
                              'input_type', opt_info%input_type, InputType)

    call read_ctrlfile_integer(handle, SectionMbar, &
                               'self_iteration', opt_info%self_iteration)

    call read_ctrlfile_integer(handle, SectionMbar, &
                               'newton_iteration', opt_info%newton_iteration)

    call read_ctrlfile_integer(handle, SectionMbar, &
                               'nblocks', opt_info%nblocks)

    ! call read_ctrlfile_real   (handle, SectionMbar, &
    !                            'temperature', opt_info%temperature)
    call read_ctrlfile_string(handle, SectionMbar, &
                               'temperature', opt_info%temperature)

    call read_ctrlfile_real   (handle, SectionMbar, &
                               'target_temperature', opt_info%target_temperature)

    call read_ctrlfile_real   (handle, SectionMbar, &
                               'tolerance', opt_info%tolerance)

    if (opt_info%target_temperature == 0.0_wp) &
      call error_msg('Read_Ctrl_Option> target_temperature should be specified.')

    if (opt_info%dimension > 2) &
      call error_msg('Read_Ctrl_Option> Dimension must be 1 or 2.')

    if (opt_info%nblocks == 0) &
      call error_msg('Read_Ctrl_Option> Nblocks must not be ZERO')

    if (opt_info%read_ref_pdb .and. opt_info%read_ref_path) &
      call error_msg('Read_Ctrl_Option> read_ref_pdb and read_ref_path can not be set at same time')

    allocate(opt_info%rest_func_no(opt_info%dimension), &
             opt_info%grids(opt_info%dimension))

    opt_info%rest_func_no(1:opt_info%dimension) = ''
    opt_info%grids(1:opt_info%dimension) = ''

    do i = 1, opt_info%dimension
      write(cdim,'(i0)') i
      call read_ctrlfile_string (handle, SectionMbar, &
                                 'rest_function'//cdim, opt_info%rest_func_no(i))
    end do

    do i = 1, opt_info%dimension
      write(cdim,'(i0)') i
      call read_ctrlfile_string (handle, SectionMbar, &
                                 'grids'//cdim, opt_info%grids(i))
    end do

    call read_ctrlfile_string (handle, SectionMbar, &
                               'output_unit', opt_info%out_unit)

    call end_ctrlfile_section(handle)

    ! RESTRAINTS section

    call begin_ctrlfile_section(handle, SectionRest)

    nfunc = 0
    do while (.true.)
      value = ''
      write(cfun,'(i0)') nfunc + 1
      call read_ctrlfile_string(handle, SectionRest, 'constant'//cfun, value)
      if (value == '') &
        exit
      nfunc = nfunc + 1
    end do
    allocate(opt_info%rest_funcs     (nfunc), &
             opt_info%rest_sel_index (nfunc), &
             opt_info%rest_constants (nfunc), &
             opt_info%rest_references(nfunc), &
             opt_info%is_periodic(nfunc), &
             opt_info%box_size(nfunc))

    opt_info%rest_funcs     (1:nfunc) = 0
    opt_info%rest_sel_index (1:nfunc) = ''
    opt_info%rest_constants (1:nfunc) = ''
    opt_info%rest_references(1:nfunc) = ''
    opt_info%is_periodic(1:nfunc)     = .false.
    opt_info%box_size(1:nfunc)        = 0.0_wp

    do i = 1, nfunc
      write(cfun,'(i0)') i

      call read_ctrlfile_type  (handle, SectionRest, &
                               'function'//cfun, opt_info%rest_funcs(i), &
                               RestraintsFuncs)

      call read_ctrlfile_string(handle, SectionRest, &
                               'select_index'//cfun, opt_info%rest_sel_index(i))

      call read_ctrlfile_string(handle, SectionRest, &
                               'constant'//cfun, opt_info%rest_constants(i))

      call read_ctrlfile_string(handle, SectionRest, &
                               'reference'//cfun, opt_info%rest_references(i))

      call read_ctrlfile_logical(handle, SectionRest, &
                               'is_periodic'//cfun, opt_info%is_periodic(i))

      call read_ctrlfile_real(handle, SectionRest, &
                               'box_size'//cfun, opt_info%box_size(i))

    end do
    
    call end_ctrlfile_section(handle)


    ! write parameter to MsgOut
    !

    write(MsgOut,'(A)') 'Read_Ctrl_Option> Parameters of MBAR'
    write(MsgOut,'(A)') ''
    if (opt_info%check_only) then
      write(MsgOut,'(A20,A3)') '  check only      = ', 'yes'
    else
      write(MsgOut,'(A20,A2)') '  check only      = ', 'no'
    end if

    write(MsgOut,'(A,I0)') '  # of replicas : ', opt_info%num_replicas

    write(MsgOut,'(A,I0)') '  # of dimension : ', opt_info%dimension

    write(MsgOut,'(A20,I0)') '  self_iteration  = ', opt_info%self_iteration

    write(MsgOut,'(A20,I0)') '  newton_iteration= ', opt_info%newton_iteration

    write(MsgOut,'(A20,I0)') '  # of blocks     = ', opt_info%nblocks

    write(MsgOut,'(A20,A)') &
      '  temperature     = ', trim(opt_info%temperature)

    write(MsgOut,'(A,F10.3)') &
      '  target_temperature = ', opt_info%target_temperature

    write(MsgOut,'(A20,F10.3)') &
      '  tolerance       = ', opt_info%tolerance

    do i = 1, opt_info%dimension

      write(MsgOut,'(A20,I10)') &
        '  dimension :      ', i
      write(MsgOut,'(A20,A)') &
        '  rest func no    = ', trim(opt_info%rest_func_no(i))
      write(MsgOut,'(A20,A)') &
        '  grids           = ', trim(opt_info%grids(i))

    end do

    write(MsgOut,'(A20,A)') &
         ' output_unit     = ', trim(opt_info%out_unit)

    write(MsgOut,'(A)') ''


    write(MsgOut,'(A)') 'Read_Ctrl_Option> Parameters of Restraints'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A,I0)') '  # of functions : ', nfunc

    do i = 1, nfunc

      write(MsgOut,'(A20,I10)') &
        '   function :       ', i
      if (opt_info%rest_funcs(i) /= 0) then
        write(MsgOut,'(A20,A)')   &
        '     restraint    = ', RestraintsFuncs(opt_info%rest_funcs(i))
        write(MsgOut,'(A20,A)')   &
        '     select index = ', trim(opt_info%rest_sel_index(i))
      end if

      len_str = len_trim(opt_info%rest_constants(i))
      if (len_str > 50) then
        write(MsgOut,'(A20,A)') &
        '     constant     = ', opt_info%rest_constants(i)(1:25)//'...'// &
                                opt_info%rest_constants(i)(len_str-25:len_str)
      else
        write(MsgOut,'(A20,A)') &
        '     constant     = ', trim(opt_info%rest_constants(i))
      end if

      len_str = len_trim(opt_info%rest_references(i))
      if (len_str > 50) then
        write(MsgOut,'(A20,A)') &
        '     reference    = ', opt_info%rest_references(i)(1:25)//'...'// &
                                opt_info%rest_references(i)(len_str-25:len_str)
      else
        write(MsgOut,'(A20,A)') &
        '     reference    = ', trim(opt_info%rest_references(i))
      end if

      if (opt_info%is_periodic(i)) then
        write(MsgOut,'(A20,A3)') '     is_periodic  = ', 'yes'
        write(MsgOut,'(A20,F10.3)') &
          '     box_size     = ', opt_info%box_size(i)
      else
        write(MsgOut,'(A20,A2)') '     is_periodic  = ', 'no'
      end if
      write(MsgOut,'(A)') ''

    end do

    write(MsgOut,'(A)') ''

    return

  end subroutine read_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_option
  !> @brief        setup option information
  !! @authors      NT
  !! @param[in]    opt_info : OPTION section control parameters information
  !! @param[in]    sel_info : SELECTION section control parameters information
  !! @param[in]    molecule : molecule information
  !! @param[inout] option   : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine setup_option(opt_info, sel_info, molecule, option)

    ! params
    integer                  :: NumIndices(size(RestraintsFuncs)) = &
                                                          (/1,2,2,1,1,3,3,4,4/)

    ! formal arguments
    type(s_opt_info),        intent(in)    :: opt_info
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_molecule),        intent(in)    :: molecule
    type(s_option),          intent(inout) :: option

    ! local variables
    integer                    :: i, j, jj, max_nrep, nrep, nval, nfunc
    character(MaxLineLong_CV)  :: str
    character(MaxLineLong_CV)  :: stmp

    real(wp), allocatable    :: val(:)
    real(wp)                 :: tmp


    option%check_only    = opt_info%check_only
    option%read_ref_pdb  = opt_info%read_ref_pdb
    option%read_ref_path = opt_info%read_ref_path

    ! MBAR settings
    !

    ! number of replicas
    option%num_replicas = opt_info%num_replicas

    ! dimension
    option%dimension = opt_info%dimension

    ! input type
    option%input_type = opt_info%input_type

    ! nblocks
    option%nblocks = opt_info%nblocks

    ! self_iteration
    option%self_iteration = opt_info%self_iteration

    ! newton_iteration
    option%newton_iteration = opt_info%newton_iteration

    ! temperature
    !option%temperature = opt_info%temperature
    allocate(option%temperature(option%num_replicas))
    i = split_num(opt_info%temperature)
    if (i == option%num_replicas) then
      call split(option%num_replicas, option%num_replicas, &
                 opt_info%temperature, &
                 option%temperature)
    else if (i == 1) then
      read(opt_info%temperature, *) tmp
      option%temperature(:) = tmp
    end if
    !write(*, *) 'temperature = ', option%temperature(:)

    ! target temperature
    option%target_temperature = opt_info%target_temperature

    ! tolerance
    option%tolerance = opt_info%tolerance

    if (option%input_type == InputTypeCV .or. option%input_type == InputTypeUS) then

      ! rest function number
      allocate(option%rest_func_num(option%dimension))
      nval = 1
      do i = 1, option%dimension
        option%rest_func_num(i) = split_num(opt_info%rest_func_no(i))
        if (nval < option%rest_func_num(i)) then
          nval = option%rest_func_num(i)
        end if
      end do

      allocate(option%rest_func_no(option%dimension, nval))
      option%rest_func_no(:, :) = 0
      do i = 1, option%dimension
        nval = option%rest_func_num(i)
        if (option%rest_func_num(i) > 0 ) then
          stmp = opt_info%rest_func_no(i)
          call tolower(stmp)
          if (stmp .eq. "all") then
            do j = 1, nval
              option%rest_func_no(i, j) = j
            end do
          else
            call split(nval, nval,           &
                       opt_info%rest_func_no(i), option%rest_func_no(i, 1:nval))
          endif
        end if
      end do

      if ((option%read_ref_path .or. option%read_ref_pdb) .and. &
        option%dimension /= 1) &
        call error_msg('Setup_Option> dimension > 1 not supported in read_ref_path/pdb')

      ! !
      allocate(option%grid_min (option%dimension), &
               option%grid_max (option%dimension), &
               option%num_grids(option%dimension))

      do i = 1, option%dimension

        ! grids
        str = opt_info%grids(i)
        do j = 1, len(str)
          if (str(j:j) == '(' .or. &
            str(j:j) == ')' .or. &
            str(j:j) == ',') &
            str(j:j) = ' '
        end do
        if (str /= '') then
          if (split_num(str) /= 3) &
            call error_msg('Setup_Option> grids must be (min max num_grids)')
          read(str,*) option%grid_min(i), &
                      option%grid_max(i), &
                      option%num_grids(i)
        else
          option%grid_min(i)  = 0.0_wp
          option%grid_max(i)  = 0.0_wp
          option%num_grids(i) = 0
        end if

      end do

    end if

    ! unit for output
    option%out_unit = opt_info%out_unit
    if (option%out_unit == 'kcal/mol') then
      write(MsgOut,'(A)') &
        'Setup_Option> Unit for fene and pmf is set to kcal/mol'
    else if (option%out_unit == 'NONE') then
      write(MsgOut,'(A)') &
        'Setup_Option> Unit for fene and pmf is set to dimension less'
    else
      call error_msg('Setup_Option> Unknown unit. out_unit should be NONE or kcal/mol')
    end if
 
    ! selection settings
    !

    ! num_atoms
    option%num_atoms = molecule%num_atoms

    ! selatoms
    allocate(option%selatoms(1:size(sel_info%groups)))

    do i = 1, size(sel_info%groups)
      call select_atom(molecule, sel_info%groups(i), option%selatoms(i))

      write(MsgOut,'(a,i12,a,a)') &
           ' group: ', i, ' :', trim(sel_info%groups(i))
      write(MsgOut,'(a,i12)') &
           '    # of selected atom = ', size(option%selatoms(i)%idx)
    end do

    write(MsgOut,'(a)') ''


    if (option%input_type == InputTypeCV .or. option%input_type == InputTypeUS) then

      ! restraints settings
      !
      nfunc = size(opt_info%rest_funcs)

      if ((option%read_ref_path .or. option%read_ref_pdb) .and. nfunc /= 1) &
        call error_msg('Setup_Option> # of function should be 1 in read_ref_path/pdb')

      allocate(option%rest_nreplica(nfunc))
      do i = 1, nfunc
        option%rest_nreplica(i) = split_num(opt_info%rest_constants(i))
      end do
      max_nrep = maxval(option%rest_nreplica(:))

      allocate(option%rest_funcs     (          nfunc), &
               option%rest_sel_index (       4, nfunc), &
               option%rest_constants (max_nrep, nfunc), &
               option%rest_references(max_nrep, nfunc), &
               option%is_periodic(nfunc), &
               option%box_size(nfunc))

      do i = 1, nfunc

        ! rest_funcs
        option%rest_funcs(i) = opt_info%rest_funcs(i)

        ! rest_sel_index
        option%rest_sel_index(1:4,i) = 0

        if (option%rest_funcs(i) /= 0) then
          nval = split_num(opt_info%rest_sel_index(i))
          if (nval /= NumIndices(option%rest_funcs(i))) &
            call error_msg('Setup_Option> Number of restraint select index is incorrect.')

          call split(nval, nval, opt_info%rest_sel_index(i), &
                   option%rest_sel_index(:,i))
        end if

        ! rest_constants
        call split(option%rest_nreplica(i), option%rest_nreplica(i), &
                   opt_info%rest_constants(i), &
                   option%rest_constants(:,i))

        ! rest_references
        if (opt_info%rest_references(i) .ne. "") then
          call split(option%rest_nreplica(i), option%rest_nreplica(i), &
                     opt_info%rest_references(i), &
                     option%rest_references(:,i))
        else
          if (option%rest_funcs(i) .ne. RestraintsFuncPOSI) &
            call error_msg('Setup_Option> Please insert Reference')
        endif

        if (option%rest_funcs(i) .ne. RestraintsFuncPOSI) then
          if (option%read_ref_path .or. option%read_ref_pdb) &
            call error_msg('Setup_Option> read_ref_path/pdb is available only in POSI')
        else
          if (opt_info%is_periodic(i)) &
            call error_msg('Setup_Option> is_periodic is not allowed in POSI')
        endif
      end do

      option%is_periodic(1:nfunc) = opt_info%is_periodic(1:nfunc)
      option%box_size(1:nfunc) = opt_info%box_size(1:nfunc)

    end if

    return

  end subroutine setup_option
  
end module ma_option_mod
