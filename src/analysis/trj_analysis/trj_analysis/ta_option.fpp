!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ta_option_mod
!> @brief   module for analysis options
!! @authors Norio Takase (NT), Takaharu Mori (TM)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ta_option_mod

  use ta_option_str_mod
  use output_str_mod
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
    integer                         :: n_distance     = 1     
    logical                         :: check_only     = .false.
    logical                         :: allow_backup   = .false.
    logical                         :: linear_sum     = .false.    
    character(Maxline), allocatable :: dis_atoms(:)
    character(Maxline), allocatable :: ang_atoms(:)
    character(Maxline), allocatable :: tor_atoms(:)
    character(Maxline), allocatable :: cdis_groups(:)
    character(Maxline), allocatable :: cang_groups(:)
    character(Maxline), allocatable :: ctor_groups(:)
    character(Maxline), allocatable :: ldis_groups(:)
    character(Maxline), allocatable :: dist_weight(:)
  end type s_opt_info

  ! subroutines
  public  :: show_ctrl_option
  public  :: read_ctrl_option
  public  :: setup_option
  private :: parse_atom_defs
  private :: parse_group_defs
  private :: get_atom_index

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_option
  !> @brief        show control parameters in OPTION section
  !! @authors      NT
  !! @param[in]    mode : show mode
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_option(mode)

    ! formal arguments
    character(*),            intent(in)    :: mode


    select case(mode)

    case ('com_distance','com_angle','com_torsion','')
      write(MsgOut,'(A)') '[SELECTION]'
      write(MsgOut,'(A)') '# group1         = an:CA              # COM group 1'
      write(MsgOut,'(A)') '# group2         = an:CA and resno:1  # COM group 2'
      write(MsgOut,'(A)') '# group3         = an:CA and resno:5  # COM group 3'
      write(MsgOut,'(A)') '# group4         = an:CA and resno:7  # COM group 4'
      write(MsgOut,'(A)') '# group5         = an:CA and resno:12 # COM group 5'
      write(MsgOut,'(A)') ''

    end select

    write(MsgOut,'(A)') '[OPTION]'
    write(MsgOut,'(A)') 'check_only     = NO              # only checking input files (YES/NO)'
    write(MsgOut,'(A)') 'allow_backup   = NO              # backup existing output files (YES/NO)'

    select case (mode)

    case ('com_distance')
      write(MsgOut,'(A)') '# com_distance1  = 1 2'
      write(MsgOut,'(A)') '# com_distance2  = 3 4'

    case ('com_angle')
      write(MsgOut,'(A)') '# com_angle1     = 1 2 3'
      write(MsgOut,'(A)') '# com_angle2     = 4 5 6'

    case ('com_torsion')
      write(MsgOut,'(A)') '# com_torsion1   = 1 2 3 4'
      write(MsgOut,'(A)') '# com_torsion2   = 5 6 7 8'

    case ('distance')
      write(MsgOut,'(A)') '# distance1      = P1:306:ALA:O   P1:768:ASN:HD21'
      write(MsgOut,'(A)') '# distance2      = P1:771:GLU:OE1 P1:768:ASN:ND2'

    case ('angle')
      write(MsgOut,'(A)') '# angle1         = P1:306:ALA:O   P1:768:ASN:HD21 P1:277:TIP3:H1'
      write(MsgOut,'(A)') '# angle2         = P1:771:GLU:OE1 P1:768:ASN:ND2  P1:272:TIP3:H1'

    case ('torsion')
      write(MsgOut,'(A)') '# torsion1       = P1:306:ALA:O   P1:768:ASN:HD21 P1:277:TIP3:H1 P1:271:TIP3:H2'
      write(MsgOut,'(A)') '# torsion2       = P1:771:GLU:OE1 P1:768:ASN:ND2  P1:272:TIP3:H1 P1:204:TIP3:OH2'

    case ('')
      write(MsgOut,'(A)') '# com_distance1  = 1 2'
      write(MsgOut,'(A)') '# com_distance2  = 3 4'
      write(MsgOut,'(A)') '# com_angle1     = 1 2 3'
      write(MsgOut,'(A)') '# com_angle2     = 4 5 6'
      write(MsgOut,'(A)') '# com_torsion1   = 1 2 3 4'
      write(MsgOut,'(A)') '# com_torsion2   = 5 6 7 8'
      write(MsgOut,'(A)') '# distance1      = P1:306:ALA:O   P1:768:ASN:HD21'
      write(MsgOut,'(A)') '# distance2      = P1:771:GLU:OE1 P1:768:ASN:ND2'
      write(MsgOut,'(A)') '# angle1         = P1:306:ALA:O   P1:768:ASN:HD21 P1:277:TIP3:H1'
      write(MsgOut,'(A)') '# angle2         = P1:771:GLU:OE1 P1:768:ASN:ND2  P1:272:TIP3:H1'
      write(MsgOut,'(A)') '# torsion1       = P1:306:ALA:O   P1:768:ASN:HD21 P1:277:TIP3:H1 P1:271:TIP3:H2'
      write(MsgOut,'(A)') '# torsion2       = P1:771:GLU:OE1 P1:768:ASN:ND2  P1:272:TIP3:H1 P1:204:TIP3:OH2'
    end select

    write(MsgOut,'(A)') ' '

    return

  end subroutine show_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_option
  !> @brief        read control parameters in OPTION section
  !! @authors      NT, TM, SI
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
    integer                  :: i, ndis, nang, ntor, ncdis, ncang, nctor, nwdis
    character(MaxLine)       :: value, dis_name, ang_name, tor_name, wgt_name


    ! read parameters
    !
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_logical(handle, Section, &
                              'check_only', opt_info%check_only)

    call read_ctrlfile_logical(handle, Section, &
                              'allow_backup', opt_info%allow_backup)

   

    ! read distance atom definitions
    !
    ndis = 0
    do while (.true.)
      value = ''
      write(dis_name,'(A8,I0)') 'distance', ndis + 1
      call read_ctrlfile_string(handle, Section, dis_name, value)
      if (value == '') &
        exit
      ndis = ndis + 1
    end do
    allocate(opt_info%dis_atoms(ndis))

    do i = 1, ndis
      write(dis_name,'(A8,I0)') 'distance', i
      call read_ctrlfile_string(handle, Section, dis_name, &
                                opt_info%dis_atoms(i))
    end do

    ! read weight for linear combination of distance 
    !
    nwdis = 0
    do while (.true.)
      value = ''
      write(wgt_name,'(A11,I0)') 'dist_weight', nwdis + 1
      call read_ctrlfile_string(handle, Section, wgt_name, value)
      if (value == '') exit
      nwdis = nwdis + 1
    end do
    if (nwdis == 0) then
      allocate(opt_info%dist_weight(ndis))
      opt_info%dist_weight = ''
    else  
      allocate(opt_info%dist_weight(nwdis))
      opt_info%dist_weight = ''
      do i = 1, nwdis
        write(wgt_name,'(A11,I0)') 'dist_weight', i
        call read_ctrlfile_string(handle, Section, wgt_name, &
                                  opt_info%dist_weight(i))
      end do 
    end if

    ! read angle atom definitions
    !
    nang = 0
    do while (.true.)
      value = ''
      write(ang_name,'(A5,I0)') 'angle', nang + 1
      call read_ctrlfile_string(handle, Section, ang_name, value)
      if (value == '') &
        exit
      nang = nang + 1
    end do
    allocate(opt_info%ang_atoms(nang))

    do i = 1, nang
      write(ang_name,'(A5,I0)') 'angle', i
      call read_ctrlfile_string(handle, Section, ang_name, &
                                opt_info%ang_atoms(i))
    end do

    ! read torsion atom definitions
    !
    ntor = 0
    do while (.true.)
      value = ''
      write(tor_name,'(A7,I0)') 'torsion', ntor + 1
      call read_ctrlfile_string(handle, Section, tor_name, value)
      if (value == '') &
        exit
      ntor = ntor + 1
    end do
    allocate(opt_info%tor_atoms(ntor))

    do i = 1, ntor
      write(tor_name,'(A7,I0)') 'torsion', i
      call read_ctrlfile_string(handle, Section, tor_name, &
                                opt_info%tor_atoms(i))
    end do

    ! read COM distance atom definitions
    !
    ncdis = 0
    do while (.true.)
      value = ''
      write(dis_name,'(A12,I0)') 'com_distance', ncdis + 1
      call read_ctrlfile_string(handle, Section, dis_name, value)
      if (value == '') &
        exit
      ncdis = ncdis + 1
    end do
    allocate(opt_info%cdis_groups(ncdis))

    do i = 1, ncdis
      write(dis_name,'(A12,I0)') 'com_distance', i
      call read_ctrlfile_string(handle, Section, dis_name, &
                                opt_info%cdis_groups(i))
    end do

    ! read COM angle atom definitions
    !
    ncang = 0
    do while (.true.)
      value = ''
      write(ang_name,'(A9,I0)') 'com_angle', ncang + 1
      call read_ctrlfile_string(handle, Section, ang_name, value)
      if (value == '') &
        exit
      ncang = ncang + 1
    end do
    allocate(opt_info%cang_groups(ncang))

    do i = 1, ncang
      write(ang_name,'(A9,I0)') 'com_angle', i
      call read_ctrlfile_string(handle, Section, ang_name, &
                                opt_info%cang_groups(i))
    end do

    ! read torsion atom definitions
    !
    nctor = 0
    do while (.true.)
      value = ''
      write(tor_name,'(A11,I0)') 'com_torsion', nctor + 1
      call read_ctrlfile_string(handle, Section, tor_name, value)
      if (value == '') &
        exit
      nctor = nctor + 1
    end do
    allocate(opt_info%ctor_groups(nctor))

    do i = 1, nctor
      write(tor_name,'(A11,I0)') 'com_torsion', i
      call read_ctrlfile_string(handle, Section, tor_name, &
                                opt_info%ctor_groups(i))
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

    if (nwdis > 0) then
      write(MsgOut,'(A31,I10)') '   # of weight for distance  = ', nwdis
      do i = 1, nwdis
        write(MsgOut,'(A25,I0,A3,A)') &
             '     weight for distance ', i, ' = ', trim(opt_info%dist_weight(i))
      end do
    end if

    if (ndis > 0) then
      write(MsgOut,'(A20,I10)') '   # of distance  = ', ndis
      do i = 1, ndis
        write(MsgOut,'(A14,I0,A3,A)') &
             '     distance ', i, ' = ', trim(opt_info%dis_atoms(i))
      end do
    end if

    if (nang > 0) then
      write(MsgOut,'(A20,I10)') '   # of angle     = ', nang
      do i = 1, nang
        write(MsgOut,'(A14,I0,A3,A)') &
             '        angle ', i, ' = ', trim(opt_info%ang_atoms(i))
      end do
    end if

    if (ntor > 0) then
      write(MsgOut,'(A20,I10)') '   # of torsion   = ', ntor
      do i = 1, ntor
        write(MsgOut,'(A14,I0,A3,A)') &
             '      torsion ', i, ' = ', trim(opt_info%tor_atoms(i))
      end do
    end if

    if (ncdis > 0) then
      write(MsgOut,'(A20,I10)') '   # of COM dist. = ', ncdis
      do i = 1, ncdis
        write(MsgOut,'(A14,I0,A3,A)') &
             '     distance ', i, ' = ', trim(opt_info%cdis_groups(i))
      end do
    end if

    if (ncang > 0) then
      write(MsgOut,'(A20,I10)') '   # of COM angle = ', ncang
      do i = 1, ncang
        write(MsgOut,'(A14,I0,A3,A)') &
             '        angle ', i, ' = ', trim(opt_info%cang_groups(i))
      end do
    end if

    if (nctor > 0) then
      write(MsgOut,'(A20,I10)') '   # of COM tor.  = ', nctor
      do i = 1, nctor
        write(MsgOut,'(A14,I0,A3,A)') &
             '      torsion ', i, ' = ', trim(opt_info%ctor_groups(i))
      end do
    end if

    write(MsgOut,'(A)') ' '

    return

  end subroutine read_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_option
  !> @brief        setup option information
  !! @authors      NT, SI
  !! @param[in]    opt_info : OPTION section control parameters information
  !! @param[in]    sel_info : SELECTION section control parameters information
  !! @param[in]    molecule : molecule information
  !! @param[inout] option   : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine setup_option(opt_info, sel_info, molecule, output, option)

    ! formal arguments
    type(s_opt_info),        intent(in)    :: opt_info
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_molecule),        intent(in)    :: molecule
    type(s_output),          intent(in)    :: output
    type(s_option),          intent(inout) :: option
    
    ! local variables
    integer                  :: i, ndis, nang, ntor, ncdis, ncang, nctor
    integer                  :: natm, nwdis, tnwdis

    ! check only
    !
    option%check_only = opt_info%check_only

    ! backup output file
    !
    call setup_backup(opt_info%allow_backup)

    ! linear combination of distance
    !
    ndis   = size(opt_info%dis_atoms)
    write(*,*) "ndis= ", ndis
    tnwdis = size(opt_info%dist_weight)

    !if ((ndis == 1) .and. (tnwdis > 1)) then
    !  call error_msg('Setup_Option> Total # of weight factor does match Toatal # of distance') 
    !else if ((ndis > 1) .and. (ndis /= tnwdis)) then
    !  call error_msg('Setup_Option> Total # of weight factor does match Toatal # of distance') 
    !end if  

    allocate(option%dist_weight(ndis, 50))
    option%dist_weight = 1.0_wp

    do i=1, ndis

      natm  = split_num(opt_info%dis_atoms(i))
      nwdis = split_num(opt_info%dist_weight(i))
      write(*,*) "natm= ", natm
      write(*,*) "nwdis= ", nwdis
      !if ((nwdis > 1) .and. (nwdis /= (natm/2))) &
      ! call error_msg('Setup_Option> # of weight factor does match # of distance')
  
      ! if nwdis == 0, weight should be defined.
      if (nwdis > 0) &
      call split(nwdis, nwdis, opt_info%dist_weight(i), option%dist_weight(i,:))

    end do       

    ! distance
    !
    !ndis = size(opt_info%dis_atoms)
    option%out_dis = ndis > 0

    if (option%out_dis) then

      if (output%disfile == '')  call error_msg('Setup_Option> ERROR : disfile name is not specified')
      if (natm > 50) call error_msg('Setup_Option> Maximum # of atom in one distance is 50')
      if (mod(natm,2) /= 0) call error_msg('Setup_Option> Maximum # of atom must be even number')

      allocate(option%dist_list(50,ndis), option%distance(ndis), option%dist_num(ndis))
      option%dist_list(50,:) = 0

      write(MsgOut,'(A)') 'Setup_Option> distance atom indices: '
      do i = 1, ndis
        call parse_atom_defs( &
               molecule, opt_info%dis_atoms(i), option%dist_list(:,i), option%dist_num(i))
        write(MsgOut,'(I5,A,50I7:)') i, ' ) ', option%dist_list(1:option%dist_num(i),i)
      end do

      write(MsgOut, '(A)') ' '

    end if

    ! angle
    !
    nang = size(opt_info%ang_atoms)
    option%out_ang = nang > 0

    if (option%out_ang) then

      if (output%angfile == '')  call error_msg('Setup_Option> ERROR : angfile name is not specified')

      allocate(option%angl_list(3,nang), option%angle(nang))

      write(MsgOut,'(A)') 'Setup_Option> angle atom indices: '
      do i = 1, nang
        call parse_atom_defs( &
             molecule, opt_info%ang_atoms(i), option%angl_list(:,i))
        write(MsgOut,'(I5,A,I7,I7,I7)') i, ' ) ', option%angl_list(:,i)
      end do
      write(MsgOut, '(A)') ' '

    end if

    ! torsion
    !
    ntor = size(opt_info%tor_atoms)
    option%out_tor = ntor > 0

    if (option%out_tor) then

      if (output%torfile == '')  call error_msg('Setup_Option> ERROR : torfile name is not specified')

      allocate(option%tors_list(4,ntor), option%torsion(ntor))

      write(MsgOut,'(A)') 'Setup_Option> torsion atom indices: '
      do i = 1, ntor
        call parse_atom_defs( &
             molecule, opt_info%tor_atoms(i), option%tors_list(:,i))
        write(MsgOut,'(I5,A,I7,I7,I7,I7)') i, ' ) ', option%tors_list(:,i)
      end do
      write(MsgOut, '(A)') ' '

    end if

    ! COM distance
    !
    ncdis = size(opt_info%cdis_groups)
    option%out_cdis = ncdis > 0

    if (option%out_cdis) then

      if (output%comdisfile == '')  call error_msg('Setup_Option> ERROR : comdisfile name is not specified')

      allocate(option%cdist_group(2,ncdis), option%cdistance(ncdis))

      write(MsgOut,'(A)') &
           'Setup_Option> COM distance atom selection expression: '
      do i = 1, ncdis
        call parse_group_defs( &
             sel_info, opt_info%cdis_groups(i), option%cdist_group(:,i))
        write(MsgOut,'(I5,A)') i, ' : '
        write(MsgOut,'(5X,A)') trim(sel_info%groups(option%cdist_group(1,i)))
        write(MsgOut,'(5X,A)') trim(sel_info%groups(option%cdist_group(2,i)))
      end do

      write(MsgOut, '(A)') ' '

    end if

    ! COM angle
    !
    ncang = size(opt_info%cang_groups)
    option%out_cang = ncang > 0

    if (option%out_cang) then

      if (output%comangfile == '')  call error_msg('Setup_Option> ERROR : comangfile name is not specified')

      allocate(option%cangl_group(3,ncang), option%cangle(ncang))

      write(MsgOut,'(A)') &
           'Setup_Option> COM angle atom selection expression: '
      do i = 1, ncang
        call parse_group_defs( &
             sel_info, opt_info%cang_groups(i), option%cangl_group(:,i))
        write(MsgOut,'(I5,A)') i, ' : '
        write(MsgOut,'(5X,A)') trim(sel_info%groups(option%cangl_group(1,i)))
        write(MsgOut,'(5X,A)') trim(sel_info%groups(option%cangl_group(2,i)))
        write(MsgOut,'(5X,A)') trim(sel_info%groups(option%cangl_group(3,i)))
      end do

      write(MsgOut, '(A)') ' '

    end if

    ! COM torsion
    !
    nctor = size(opt_info%ctor_groups)
    option%out_ctor = nctor > 0

    if (option%out_ctor) then

      if (output%comtorfile == '')  call error_msg('Setup_Option> ERROR : comtorfile name is not specified')

      allocate(option%ctor_group(4,nctor), option%ctorsion(nctor))

      write(MsgOut,'(A)') &
           'Setup_Option> COM torsion atom selection expression: '
      do i = 1, nctor
        call parse_group_defs( &
             sel_info, opt_info%ctor_groups(i), option%ctor_group(:,i))
        write(MsgOut,'(I5,A)') i, ' : '
        write(MsgOut,'(5X,A)') trim(sel_info%groups(option%ctor_group(1,i)))
        write(MsgOut,'(5X,A)') trim(sel_info%groups(option%ctor_group(2,i)))
        write(MsgOut,'(5X,A)') trim(sel_info%groups(option%ctor_group(3,i)))
        write(MsgOut,'(5X,A)') trim(sel_info%groups(option%ctor_group(4,i)))
      end do

      write(MsgOut, '(A)') ' '

    end if


    if (option%out_cdis .or. option%out_cang .or. option%out_ctor) then

      ! select atoms
      !
      allocate(option%selatoms(1:size(sel_info%groups)))

      do i = 1, size(sel_info%groups)
        call select_atom(molecule, sel_info%groups(i), option%selatoms(i))
      end do

      ! check molecule mass
      !
      do i = 1, molecule%num_atoms
        if (molecule%mass(i) == 0.0_wp) &
          exit
      end do

      if (i <= molecule%num_atoms) then
        call error_msg('Setup_Option> ERROR : molecule mass were not assigned.')
      end if

    end if

    return

  end subroutine setup_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    parse_atom_defs
  !> @brief        parse distance/angle/torsion atom definition line
  !! @authors      NT, SI
  !! @param[in]    molecule : molecule information
  !! @param[in]    defs     : atom definitions
  !! @param[inout] atoms    : atom indices
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine parse_atom_defs(molecule, defs, atoms, n_dist)

    ! formal argments
    type(s_molecule),        intent(in)    :: molecule
    character(*),            intent(in)    :: defs
    integer,                 intent(inout) :: atoms(:)
    integer, optional,       intent(inout) :: n_dist    

    ! local variables
    integer                  :: idx, i
    character(Maxline)       :: s(50)


    atoms(:) = -1

    s(:) = ''
    read(defs, *, err=900, end=100) (s(i), i=1, 50)

100 do i = 1, 50
      if (s(i) /= '') then
        idx = get_atom_index(molecule, s(i))
        if (idx == -1) &
          call error_msg('Parse_atom_defs> read error. atom not found : '//&
                         trim(s(i)))
        if (size(atoms) < i) &
          call error_msg('Parse_atom_defs> read error. bad format : '//&
                         trim(defs))
        atoms(i) = idx
      else
        exit
      end if
    end do

    do i = 1, size(atoms)
      if (atoms(i) == -1) then
        n_dist = i - 1
        exit
      end if
    end do

    return

900 call error_msg('Parse_atom_defs> read error.')

  end subroutine parse_atom_defs

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    parse_group_defs
  !> @brief        parse distance/angle/torsion group definition line
  !! @authors      NT
  !! @param[in]    sel_info : SELECTION section control parameters information
  !! @param[in]    defs     : group definitions
  !! @param[inout] groups   : atom selection groups
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine parse_group_defs(sel_info, defs, groups)

    ! formal argments
    type(s_sel_info),        intent(in)    :: sel_info
    character(*),            intent(in)    :: defs
    integer,                 intent(inout) :: groups(:)

    ! local variables
    integer                  :: i, ngroup
    character(100)           :: str


    ngroup = split_num(defs)
    if (ngroup /= size(groups)) then
      write(str,'(i0)') size(groups)
      call error_msg('Parse_Group_Defs> ERROR: group count must BE '//trim(str))
    end if

    call split(ngroup, ngroup, defs, groups)

    do i = 1, ngroup
      if (groups(i) > size(sel_info%groups)) &
        call error_msg('Parse_Group_Defs> ERROR: group index is out of bounds.')
    end do

    return

  end subroutine parse_group_defs

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      get_atom_index
  !> @brief        get atom index
  !! @authors      NT
  !! @return       atom index
  !! @param[in]    molecule : molecule information
  !! @param[in]    def      : atom definition
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  function get_atom_index(molecule, def)

    ! return value
    integer                  :: get_atom_index

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    character(*),            intent(in)    :: def

    ! local variables
    integer                  :: i, res_no, ntok
    character(Maxline)       :: cstr
    character(80)            :: an, sn, seg_name, res_name, atom_name


    get_atom_index = -1

    cstr = def
    ntok = 1
    do i = 1, len(cstr)
      if (cstr(i:i) == ':') then
        cstr(i:i) = ' '
        ntok = ntok + 1
      end if
    end do

    if (ntok < 3) then
      goto 900

    else if (ntok == 3) then
      
      read(cstr,*,err=900,end=900) res_no, res_name, atom_name
      do i = 1, molecule%num_atoms
        read(molecule%atom_name(i),*,err=900) an
        if (res_no    == molecule%residue_no(i)   .and. &
            res_name  == molecule%residue_name(i) .and. &
            atom_name == an) then
          get_atom_index = i
          return 
        end if
      end do

    else

      read(cstr,*,err=900,end=900) seg_name, res_no, res_name, atom_name
      do i = 1, molecule%num_atoms
        read(molecule%segment_name(i),*,err=900) sn
        read(molecule%atom_name(i),*,err=900) an
        if (seg_name  == sn                       .and. &
            res_no    == molecule%residue_no(i)   .and. &
            res_name  == molecule%residue_name(i) .and. &
            atom_name == an) then
          get_atom_index = i
          return 
        end if
      end do

    end if

    return 

900 call error_msg('Get_Atom_Index> read error '//def)

  end function get_atom_index

end module ta_option_mod
