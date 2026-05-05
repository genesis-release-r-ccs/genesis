!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_reac_coord_mod
!> @brief   reaction coordinate
!! @authors Junichiro Koga (JK), Shingo Ito (SI)
!
!  (c) Copyright 2020 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8
  
#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_reac_coord_mod
  use constants_mod
  use at_restraints_str_mod
  use at_reac_coord_str_mod
  !use at_abmd_mod
  !use at_bias_mod 
  !use at_bias_str_mod   
  use at_enefunc_str_mod
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
  type, public :: s_reac_info
     integer                         :: nreac_coord = 0
     integer,            allocatable :: function(:)
     character(MaxLine), allocatable :: select_index(:)
     integer,            allocatable :: direction(:)
     character(MaxLine), allocatable :: unitv(:)
     character(MaxLine), allocatable :: exponent_dist(:)
     character(MaxLine), allocatable :: weight_dist(:)

     integer,            allocatable :: nxi(:)
     real(wp),           allocatable :: ximin(:), ximax(:), dxi(:)
     real(wp),           allocatable :: uboundary(:)
     integer                         :: monitor_period
     logical,            allocatable :: periodicity(:)
  end type s_reac_info

  ! subroutines
  public  :: show_ctrl_reac_coord
  public  :: read_ctrl_reac_coord
  public  :: setup_reac_coord

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_reac_coord
  !> @brief        show ReacCoord section usage
  !! @authors      JK
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "md", "remd"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_reac_coord(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode


    if (show_all) then

      select case (run_mode)

      case ('md', 'remd')

        write(MsgOut,'(A)') '[SELECTION]'
        write(MsgOut,'(A)') '# group1        = an:CA              # ReacCoord group 1'
        write(MsgOut,'(A)') '# group2        = an:CA and resno:1  # ReacCoord group 2'
        write(MsgOut,'(A)') '# group3        = an:CA and resno:5  # ReacCoord group 3'
        write(MsgOut,'(A)') '# group4        = an:CA and resno:7  # ReacCoord group 4'
        write(MsgOut,'(A)') '# group5        = an:CA and resno:12 # ReacCoord group 5'
        write(MsgOut,'(A)') ' '

        write(MsgOut,'(A)') '[REAC_COORD]'
        write(MsgOut,'(A)') '## Example of reaction coordinates'
        write(MsgOut,'(A)') '# nfunctions        = 2         # number of reaction coordinates'
        write(MsgOut,'(A)') '#'                 
        write(MsgOut,'(A)') '# function1         = POSI      # reaction coordinate type'
        write(MsgOut,'(A)') '# direction1        = ALL       # direction [ALL,X,Y,Z]'
        write(MsgOut,'(A)') '# select_index1     = 1         # associated groups'
        write(MsgOut,'(A)') '#'                 
        write(MsgOut,'(A)') '# function2         = DIST      # reaction coordinate type'
        write(MsgOut,'(A)') '# select_index2     = 2 3       # assocaited groups'
        write(MsgOut,'(A)') '#'                 
        write(MsgOut,'(A)') '# function3         = DIHED     # reaction coordinate type'
        write(MsgOut,'(A)') '# select_index3     = 4 5 6 7   # assocaited groups'
        write(MsgOut,'(A)') '# periodicity3      = YES       # periodicity of reaction coordinates'
        write(MsgOut,'(A)') ' '

      end select

    end if

    return
  end subroutine show_ctrl_reac_coord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      read_ctrl_reac_coord
  !> @brief        read ReacCoord section in the control file
  !! @authors      JK
  !! @param[in]    handle   : unit number of control files
  !! @param[out]   reac_info : ReacCoord section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_reac_coord(handle, reac_info) 

    ! parameters
    character(*),            parameter     :: Section = 'ReacCoord'

    ! formal arguments
    integer,                 intent(in)    :: handle
    type(s_reac_info),       intent(inout) :: reac_info 

    ! local variables
    integer                  :: i, nfunc
    character(20)            :: name


    ! read parameters from control file
    ! 
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_integer(handle, Section,                           &
                       'nreac_coord', reac_info%nreac_coord)

    reac_info%monitor_period = 1000
    call read_ctrlfile_integer(handle, Section,                           &
                       'monitor_period', reac_info%monitor_period)
    if (reac_info%nreac_coord == 0) then

      if (main_rank) then
        write(MsgOut,'(a)') 'Read_Ctrl_ReacCoord> No reaction coordinate in the system'
        write(MsgOut,'(a)') ''
      end if
      return

    end if

    ! allocate for reac_info
    !
    nfunc = reac_info%nreac_coord

    allocate(reac_info%function      (nfunc), &
             reac_info%select_index  (nfunc), &
             reac_info%unitv         (nfunc), &
             reac_info%direction     (nfunc), &
             reac_info%exponent_dist (nfunc), &
             reac_info%weight_dist   (nfunc), &
             reac_info%nxi           (nfunc), &
             reac_info%ximin         (nfunc), &
             reac_info%ximax         (nfunc), &
             reac_info%uboundary     (nfunc), &
             reac_info%dxi           (nfunc), &
             reac_info%periodicity   (nfunc))

    reac_info%function     (1:nfunc) = ReacCoordsFuncPOSI
    reac_info%select_index (1:nfunc) = ''
    reac_info%exponent_dist(1:nfunc) = ''
    reac_info%weight_dist  (1:nfunc) = ''
    reac_info%nxi          (1:nfunc) = 0
    reac_info%ximin        (1:nfunc) = 0_wp
    reac_info%ximax        (1:nfunc) = 0_wp
    reac_info%uboundary    (1:nfunc) = 0_wp !200_wp
    reac_info%dxi          (1:nfunc) = 0_wp
    reac_info%direction    (1:nfunc) = 1
    reac_info%unitv        (1:nfunc) = ''
    reac_info%periodicity(1:nfunc) = .false.


    ! read each ReacCoord function
    !
    do i = 1,nfunc

      write(name,'(a,i0)') 'function', i
      call read_ctrlfile_type   (handle, Section, name, &
                                 reac_info%function(i), ReacCoordsFuncTypes)

      write(name,'(a,i0)') 'select_index', i
      call read_ctrlfile_string (handle, Section, name, &
                                 reac_info%select_index(i))

      write(name,'(a,i0)') 'direction', i
      call read_ctrlfile_type   (handle, Section, name, &
                                 reac_info%direction(i), RestraintsDirTypes)

      write(name,'(a,i0)') 'unit_vector', i
      call read_ctrlfile_string (handle, Section, name, &
                                 reac_info%unitv(i))

      write(name,'(a,i0)') 'exponent_dist', i
      call read_ctrlfile_string (handle, Section, name, &
                                 reac_info%exponent_dist(i))

      write(name,'(a,i0)') 'weight_dist', i
      call read_ctrlfile_string (handle, Section, name, &
                                 reac_info%weight_dist(i))

      write(name,'(a,i0)') 'nxi', i
      call read_ctrlfile_integer(handle, Section, name, &
                                 reac_info%nxi(i))

      write(name,'(a,i0)') 'ximin', i
      call read_ctrlfile_real   (handle, Section, name, &
                                 reac_info%ximin(i))

      write(name,'(a,i0)') 'ximax', i
      call read_ctrlfile_real   (handle, Section, name, &
                                 reac_info%ximax(i))

      write(name,'(a,i0)') 'dxi', i
      call read_ctrlfile_real   (handle, Section, name, &
                                 reac_info%dxi(i))

      if (reac_info%function(i) == ReacCoordsFuncANGLE    .or. &
          reac_info%function(i) == ReacCoordsFuncANGLECOM .or. &
          reac_info%function(i) == ReacCoordsFuncDIHED    .or. &
          reac_info%function(i) == ReacCoordsFuncDIHEDCOM ) then
        reac_info%uboundary(i) = 0_wp
      endif
      write(name,'(a,i0)') 'uboundary', i
      call read_ctrlfile_real   (handle, Section, name, &
                                 reac_info%uboundary(i))

      write(name,'(a,i0)') 'periodicity', i
      call read_ctrlfile_logical(handle, Section, name,  &
                                 reac_info%periodicity(i))

      if (reac_info%select_index(i) == '') &
        call error_msg('Read_Ctrl_ReacCoord> No index is found')


    end do

    call end_ctrlfile_section(handle)


    ! write parameters to MsgOut
    !
    if (main_rank) then

      write(MsgOut,'(a)')  'Read_Ctrl_ReacCoord> Parameters of ReactionCoordinate'

      write(MsgOut,'(a20,i3)')                            &
            '  nreac_coord       = ', nfunc
      write(MsgOut,'(a,i10)')   &
            '  monitor_period    = ', reac_info%monitor_period

      do i = 1,nfunc
        write(MsgOut,'(a,i0,a,a)') &
              '  function', i, ' = ', &
                            ReacCoordsFuncTypes(reac_info%function(i))

        write(MsgOut,'(a,a)') &
              '    select_index  = ', trim(reac_info%select_index(i))

        write(MsgOut,'(a,a)') &
              '    exponent_dist = ', trim(reac_info%exponent_dist(i))

        write(MsgOut,'(a,a)') &
              '    weight_dist   = ', trim(reac_info%weight_dist(i))

        write(MsgOut,'(a,a)') &
              '    unit_vector   = ', trim(reac_info%unitv(i))

        write(MsgOut,'(a,i0)') &
              '  nxi             = ', reac_info%nxi(i)

        write(MsgOut,'(a,f10.5)') &
              '  ximin           = ', reac_info%ximin(i)

        write(MsgOut,'(a,f10.5)') &
              '  ximax           = ', reac_info%ximax(i)

        write(MsgOut,'(a,f10.5)') &
              '  dxi             = ', reac_info%dxi(i)
              
        write(MsgOut,'(a,f10.5)') &
              '  uboundary       = ', reac_info%uboundary(i)

        write(MsgOut,'(a,l)') &
              '  periodicity     = ', reac_info%periodicity(i)      

      end do

      write(MsgOut,'(a)')

    end if

    return
  end subroutine read_ctrl_reac_coord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      setup_reac_coord
  !> @brief        setup reaction coordinate information
  !! @authors      JK
  !! @param[in]    reac_info  : ReacCoord section control parameter information
  !! @param[in]    sel_info   : SELECTIONS section control parameter information
  !! @param[inout] molecule   : molecule information
  !! @param[inout] reac_coord : reac_coord information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine setup_reac_coord(reac_info, sel_info, molecule, &
                              reac_coord)

    ! formal arguments
    type(s_reac_info),       intent(in)    :: reac_info
    !type(s_abmd_info),       intent(in)    :: abmd_info
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_molecule),        intent(inout) :: molecule
    type(s_reac_coord),      intent(inout) :: reac_coord

    ! local variables
    integer                  :: nfunc, ngroup, natom
    integer                  :: max_natom
    integer                  :: alloc_stat
    integer                  :: dealloc_stat
    integer                  :: i, j

    type(s_selatoms),        allocatable :: selatoms(:)

    reac_coord%nreac_coord    = reac_info%nreac_coord
    reac_coord%num_groups     = size(sel_info%groups)
    reac_coord%monitor_period = reac_info%monitor_period

    nfunc  = reac_coord%nreac_coord
    ngroup = reac_coord%num_groups

    !if (abmd_info%abmd .eqv. .false.) then
    !  call error_msg('Setup_ReacCoord> ABMD is not set')
    !end if

    ! allocate for reac coord
    !
    call alloc_reac_coord(reac_coord, RestraintsFunc,  nfunc)
    call alloc_reac_coord(reac_coord, RestraintsGroup, ngroup)

    if (nfunc > 0) then
      reac_coord%function     (1:nfunc) = reac_info%function     (1:nfunc)    
      reac_coord%select_index (1:nfunc) = reac_info%select_index (1:nfunc)
      reac_coord%direction    (1:nfunc) = reac_info%direction    (1:nfunc)
      reac_coord%unitv        (1:nfunc) = reac_info%unitv        (1:nfunc)
      reac_coord%exponent_dist(1:nfunc) = reac_info%exponent_dist(1:nfunc)
      reac_coord%weight_dist  (1:nfunc) = reac_info%weight_dist  (1:nfunc)
      reac_coord%group       (1:ngroup) = sel_info%groups       (1:ngroup)
    endif

    reac_coord%reac_coord_flag  = (reac_coord%nreac_coord > 0)

    ! set parameters for selection
    !
    call setup_selection(sel_info, molecule)

    ! selection
    !
    allocate(selatoms(ngroup), stat=alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg('Setup_ReacCoord> Memory allocation error.')

    max_natom = 0
    do i = 1, ngroup
      call select_atom(molecule, sel_info%groups(i), selatoms(i))
      natom     = size(selatoms(i)%idx)
      max_natom = max (max_natom, natom)
    end do

    call alloc_reac_coord(reac_coord, RestraintsList, ngroup, max_natom)

    reac_coord%max_atoms = max_natom

    do i = 1, ngroup

      natom = size(selatoms(i)%idx)
      reac_coord%num_atoms(i) = natom
      reac_coord%atomlist(1:max_natom,i) = 0

      do j = 1, natom
        reac_coord%atomlist(j,i) = selatoms(i)%idx(j)
      end do

    end do

    call alloc_reac_coord(reac_coord, ReacCoord, nfunc)

    reac_coord%nxi      (1:nfunc) = reac_info%nxi      (1:nfunc)
    reac_coord%ximin    (1:nfunc) = reac_info%ximin    (1:nfunc)
    reac_coord%ximax    (1:nfunc) = reac_info%ximax    (1:nfunc)
    reac_coord%dxi      (1:nfunc) = reac_info%dxi      (1:nfunc)
    reac_coord%uboundary(1:nfunc) = reac_info%uboundary(1:nfunc)
    reac_coord%periodicity(1:nfunc) = reac_info%periodicity(1:nfunc)

    ! deallocate local array
    !
    do i = 1, ngroup
      call dealloc_selatoms(selatoms(i))
    end do
    deallocate(selatoms, stat=dealloc_stat)
    if (dealloc_stat /= 0) &
      call error_msg('Setup_ReacCoord> Memory deallocation error.')
    
    return

  end subroutine setup_reac_coord

end module at_reac_coord_mod

