!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_restraints_mod
!> @brief   restraints function
!! @authors Chigusa Kobayashi (CK), Takaharu Mori (TM), Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8
  
#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_restraints_mod

  use sp_restraints_str_mod
  use sp_enefunc_str_mod
  use select_mod
  use select_atoms_mod
  use select_atoms_str_mod
  use molecules_str_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod

  implicit none
  private

  ! structures
  type, public :: s_res_info
     integer                         :: nfunctions = 0
     integer,            allocatable :: function(:)
     character(MaxLine), allocatable :: constant(:) 
     character(MaxLine), allocatable :: reference(:) 
     character(MaxLine), allocatable :: select_index(:)
     integer,            allocatable :: direction(:)
     integer,            allocatable :: exponent(:)
     character(MaxLine), allocatable :: exponent_dist(:)
     character(MaxLine), allocatable :: weight_dist(:)
     character(MaxLine), allocatable :: vector(:) 
     logical                         :: pressure_position = .false.
     logical                         :: pressure_rmsd   = .false.
  end type s_res_info

  ! subroutines
  public  :: show_ctrl_restraints
  public  :: read_ctrl_restraints
  public  :: setup_restraints
  public  :: setup_restraints_pio

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_restraints
  !> @brief        show RESTRAINTS section usage
  !! @authors      NT
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "md", "min"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_restraints(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode


    if (show_all) then

      select case (run_mode)

      case ('md', 'min')

        write(MsgOut,'(A)') '[SELECTION]'
        write(MsgOut,'(A)') '# group1        = an:CA              # restraint group 1'
        write(MsgOut,'(A)') '# group2        = an:CA and resno:1  # restraint group 2'
        write(MsgOut,'(A)') '# group3        = an:CA and resno:5  # restraint group 3'
        write(MsgOut,'(A)') '# group4        = an:CA and resno:7  # restraint group 4'
        write(MsgOut,'(A)') '# group5        = an:CA and resno:12 # restraint group 5'
        write(MsgOut,'(A)') ' '

        write(MsgOut,'(A)') '[RESTRAINTS]'
        write(MsgOut,'(A)') '## Example of restraints'
        write(MsgOut,'(A)') '# nfunctions        = 2         # number of functions'
        write(MsgOut,'(A)') '#'                 
        write(MsgOut,'(A)') '# function1         = POSI      # restraint function type'
        write(MsgOut,'(A)') '# direction1        = ALL       # direction [ALL,X,Y,Z]'
        write(MsgOut,'(A)') '# constant1         = 10.0      # force constant'
        write(MsgOut,'(A)') '# select_index1     = 1         # restrained groups'
        write(MsgOut,'(A)') '#'                 
        write(MsgOut,'(A)') '# function2         = DIST      # restraint function type'
        write(MsgOut,'(A)') '# constant2         = 10.0      # force constants'
        write(MsgOut,'(A)') '# reference2        = 5.5       # references'
        write(MsgOut,'(A)') '# select_index2     = 2 3       # restrained groups'
        write(MsgOut,'(A)') '# pressure_position = NO        # including position restraint virial for pressure'
        write(MsgOut,'(A)') '# pressure_rmsd     = NO        # including rmsd restraint virial for pressure'
        write(MsgOut,'(A)') ' '


      end select

    end if

    return

  end subroutine show_ctrl_restraints
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      read_ctrl_restraints
  !> @brief        read RESTRAINTS section in the control file
  !! @authors      CK, TM
  !! @param[in]    handle   : unit number of control files
  !! @param[out]   res_info : RESTRAINTS section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_restraints(handle, res_info) 

    ! parameters
    character(*),            parameter     :: Section = 'Restraints'

    ! formal arguments
    integer,                 intent(in)    :: handle
    type(s_res_info),        intent(inout) :: res_info 

    ! local variables
    integer                  :: i, nfunc
    character(20)            :: name


    ! read parameters from control file
    ! 
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_integer(handle, Section,                           &
                       'nfunctions', res_info%nfunctions)


    call read_ctrlfile_logical(handle, Section,                           &
                       'pressure_position', res_info%pressure_position)

    call read_ctrlfile_logical(handle, Section,                           &
                       'pressure_rmsd', res_info%pressure_rmsd)

    ! check parameters
    !
    if (res_info%nfunctions == 0) then

      if (main_rank) then
        write(MsgOut,'(a)') 'Read_Ctrl_Restraints> No restraint in the system'
        write(MsgOut,'(a)') ''
      end if
      return

    end if


    ! allocate for res_info
    !
    nfunc = res_info%nfunctions

    allocate(res_info%function     (nfunc), &
             res_info%constant     (nfunc), &
             res_info%reference    (nfunc), &
             res_info%select_index (nfunc), &
             res_info%direction    (nfunc), &
             res_info%exponent     (nfunc), &
             res_info%exponent_dist(nfunc), &
             res_info%weight_dist  (nfunc), &
             res_info%vector       (nfunc))

    res_info%function     (1:nfunc) = RestraintsFuncPOSI
    res_info%constant     (1:nfunc) = ''
    res_info%reference    (1:nfunc) = ''
    res_info%select_index (1:nfunc) = ''
    res_info%direction    (1:nfunc) = RestraintsDirALL
    res_info%exponent     (1:nfunc) = 2
    res_info%exponent_dist(1:nfunc) = ''
    res_info%weight_dist  (1:nfunc) = ''
    res_info%vector       (1:nfunc) = ''


    ! read each restraint function
    !
    do i = 1,nfunc

      write(name,'(a,i0)') 'function', i
      call read_ctrlfile_type   (handle, Section, name, &
                                 res_info%function(i), RestraintsFuncTypes)

      write(name,'(a,i0)') 'constant', i
      call read_ctrlfile_string (handle, Section, name, &
                                 res_info%constant(i))

      write(name,'(a,i0)') 'reference', i
      call read_ctrlfile_string (handle, Section, name, &
                                 res_info%reference(i))

      write(name,'(a,i0)') 'select_index', i
      call read_ctrlfile_string (handle, Section, name, &
                                 res_info%select_index(i))

      write(name,'(a,i0)') 'direction', i
      call read_ctrlfile_type   (handle, Section, name, &
                                 res_info%direction(i), RestraintsDirTypes)

      write(name,'(a,i0)') 'exponent', i
      call read_ctrlfile_integer(handle, Section, name, &
                                 res_info%exponent(i))

      write(name,'(a,i0)') 'exponent_dist', i
      call read_ctrlfile_string (handle, Section, name, &
                                 res_info%exponent_dist(i))

      write(name,'(a,i0)') 'weight_dist', i
      call read_ctrlfile_string (handle, Section, name, &
                                 res_info%weight_dist(i))

      write(name,'(a,i0)') 'vector', i
      call read_ctrlfile_string (handle, Section, name, &
                                 res_info%vector(i))

      if (res_info%exponent(i) <= 0) &
        call error_msg('Read_Ctrl_Restraints> exponent is incorrect')

      if (res_info%select_index(i) == '') &
        call error_msg('Read_Ctrl_Restraints> No index is found')

      if (res_info%function(i) == RestraintsFuncPC .or. &
          res_info%function(i) == RestraintsFuncPCCOM ) &
        call error_msg('Read_Ctrl_Restraints> PC restraint is not allowed in spdyn')

    end do

    call end_ctrlfile_section(handle)


    ! write parameters to MsgOut
    !
    if (main_rank) then

      write(MsgOut,'(a)')  'Read_Ctrl_Restraints> Parameters of Restraint'

      write(MsgOut,'(a20,i3)')                            &
            '  numfunc         = ', nfunc

      do i = 1,nfunc
        write(MsgOut,'(a,i0,a,a)') &
              '  function', i, ' = ', &
                            RestraintsFuncTypes(res_info%function(i))

        if (res_info%function(i) == RestraintsFuncPOSI) then
          write(MsgOut,'(a,a)') &
              '    direction     = ',RestraintsDirTypes(res_info%direction(i))
        end if


        write(MsgOut,'(a,a)') &
              '    constant      = ', trim(res_info%constant(i))

        write(MsgOut,'(a,a)') &
              '    select_index  = ', trim(res_info%select_index(i))

        write(MsgOut,'(a,a)') &
              '    reference     = ', trim(res_info%reference(i))

        write(MsgOut,'(a,i3)') &
                '    exponent      = ', res_info%exponent(i)
  
        write(MsgOut,'(a,a)') &
                '    exponent_dist = ',trim(res_info%exponent_dist(i))
  
        write(MsgOut,'(a,a)') &
                '    weight_dist   = ',trim(res_info%weight_dist(i))

      end do

      if (res_info%pressure_position) then
        write(MsgOut,'(a)') '  pressure_position   = YES'
      else
        write(MsgOut,'(a)') '  pressure_position   = NO'
      end if

      if (res_info%pressure_rmsd) then
        write(MsgOut,'(a)') '  pressure_rmsd     = YES'
      else
        write(MsgOut,'(a)') '  pressure_rmsd     = NO'
      end if

      write(MsgOut,'(a)')

    end if

    return

  end subroutine read_ctrl_restraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      setup_restraints
  !> @brief        setup restraints information
  !! @authors      CK, TM
  !! @param[in]    res_info   : RESTRAINTS section control parameter information
  !! @param[in]    sel_info   : SELECTIONS section control parameter information
  !! @param[inout] molecule   : molecule information
  !! @param[inout] restraints : restraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine setup_restraints(res_info, sel_info, molecule, restraints)

    ! formal arguments
    type(s_res_info),        intent(in)    :: res_info
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_molecule),        intent(inout) :: molecule
    type(s_restraints),      intent(inout) :: restraints

    ! local variables
    integer                  :: nfunc, ngroup, natom
    integer                  :: max_natom
    integer                  :: alloc_stat
    integer                  :: dealloc_stat
    integer                  :: i, j

    type(s_selatoms),        allocatable :: selatoms(:)


    restraints%nfunctions  = res_info%nfunctions
    restraints%num_groups = size(sel_info%groups)

    nfunc  = restraints%nfunctions
    ngroup = restraints%num_groups


    ! allocate for restraints
    !
    call alloc_restraints(restraints, RestraintsFunc,  nfunc)
    call alloc_restraints(restraints, RestraintsGroup, ngroup)

    if (nfunc > 0) then
      restraints%function     (1:nfunc) = res_info%function     (1:nfunc)    
      restraints%constant     (1:nfunc) = res_info%constant     (1:nfunc)
      restraints%reference    (1:nfunc) = res_info%reference    (1:nfunc)
      restraints%select_index (1:nfunc) = res_info%select_index (1:nfunc)
      restraints%direction    (1:nfunc) = res_info%direction    (1:nfunc)
      restraints%exponent     (1:nfunc) = res_info%exponent     (1:nfunc)
      restraints%exponent_dist(1:nfunc) = res_info%exponent_dist(1:nfunc)
      restraints%weight_dist  (1:nfunc) = res_info%weight_dist  (1:nfunc)
      restraints%group       (1:ngroup) = sel_info%groups      (1:ngroup)
    end if

    restraints%restraint_flag     = (restraints%nfunctions > 0)
    restraints%pressure_position  = res_info%pressure_position
    restraints%pressure_rmsd      = res_info%pressure_rmsd

    ! set parameters for selection
    !
    call setup_selection(sel_info, molecule)


    ! selection
    !
    allocate(selatoms(ngroup), stat=alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg('Setup_Restraints> Memory allocation error.')

    max_natom = 0
    do i = 1, ngroup
      call select_atom(molecule, sel_info%groups(i), selatoms(i))
      natom     = size(selatoms(i)%idx)
      max_natom = max (max_natom, natom)
    end do

    call alloc_restraints(restraints, RestraintsList, ngroup, max_natom)

    restraints%max_atoms = max_natom

    do i = 1, ngroup

      natom = size(selatoms(i)%idx)
      restraints%num_atoms(i) = natom
      restraints%atomlist(1:max_natom,i) = 0

      do j = 1, natom
        restraints%atomlist(j,i) = selatoms(i)%idx(j)
      end do

    end do


    ! deallocate local array
    !
    do i = 1, ngroup
      call dealloc_selatoms(selatoms(i))
    end do
    deallocate(selatoms, stat=dealloc_stat)
    if (dealloc_stat /= 0) &
      call error_msg('Setup_Restraints> Memory deallocation error.')
    
    return

  end subroutine setup_restraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      setup_restraints_pio
  !> @brief        setup restraints information
  !! @authors      JJ
  !! @param[in]    res_info   : RESTRAINTS section control parameter information
  !! @param[in]    sel_info   : SELECTIONS section control parameter information
  !! @param[inout] restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_restraints_pio(res_info, sel_info, restraints, enefunc)

    ! formal arguments
    type(s_res_info),        intent(in)    :: res_info
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_restraints),      intent(inout) :: restraints
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: nfunc, ngroup, natom
    integer                  :: max_natom
    integer                  :: alloc_stat
    integer                  :: dealloc_stat
    integer                  :: i, j


    restraints%nfunctions  = res_info%nfunctions
    restraints%num_groups = size(sel_info%groups)

    nfunc  = restraints%nfunctions
    ngroup = restraints%num_groups

    ! allocate for restraints
    !
    call alloc_restraints(restraints, RestraintsFunc,  nfunc)
    call alloc_restraints(restraints, RestraintsGroup, ngroup)

    if (nfunc > 0) then
      restraints%function     (1:nfunc) = res_info%function     (1:nfunc)
      restraints%constant     (1:nfunc) = res_info%constant     (1:nfunc)
      restraints%reference    (1:nfunc) = res_info%reference    (1:nfunc)
      restraints%select_index (1:nfunc) = res_info%select_index (1:nfunc)
      restraints%direction    (1:nfunc) = res_info%direction    (1:nfunc)
      restraints%exponent     (1:nfunc) = res_info%exponent     (1:nfunc)
      restraints%exponent_dist(1:nfunc) = res_info%exponent_dist(1:nfunc)
      restraints%weight_dist  (1:nfunc) = res_info%weight_dist  (1:nfunc)
      restraints%group       (1:ngroup) = sel_info%groups      (1:ngroup)
    end if

    restraints%restraint_flag     = (restraints%nfunctions > 0)
    restraints%pressure_position  = res_info%pressure_position
    restraints%pressure_rmsd      = res_info%pressure_rmsd

    ! selection is already done by selfile

    return

  end subroutine setup_restraints_pio

end module sp_restraints_mod
