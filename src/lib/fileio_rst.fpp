!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_rst
!> @brief   Restart file I/O module
!! @authors Takaharu Mori (TM), Chigusa Kobayashi (CK), Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fileio_rst_mod

  use fileio_data_mod
  use fileio_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_rst
    integer                       :: rstfile_type    = 0
    integer                       :: num_atoms       = 0
    integer                       :: iseed           = 0
    real(wp),         allocatable :: coord(:,:)
    real(wp),         allocatable :: velocity(:,:)
    real(wp)                      :: thermostat_momentum  = 0.0_wp
    real(wp)                      :: barostat_momentum(3) = 0.0_wp
    real(wp)                      :: box_size_x      = 0.0_wp
    real(wp)                      :: box_size_y      = 0.0_wp
    real(wp)                      :: box_size_z      = 0.0_wp
    ! min
    real(wp)                      :: energy          = 0.0_wp
    real(wp)                      :: delta_r         = 0.0_wp
    ! remd
    integer                       :: iseed_remd      = 0
    integer                       :: iseed_rpath     = 0
    integer                       :: dimension       = 0
    integer                       :: nreplicas       = 0
    integer,          allocatable :: repid2parmsetid(:)
    integer,          allocatable :: num_criteria(:,:,:)
    integer,          allocatable :: num_exchanges(:,:,:)
    real(wp),         allocatable :: rest_reference(:,:,:)
    ! random
    character,        allocatable :: random(:)
    ! spherical potential
    logical                       :: sph_pot     = .false.
    integer                       :: nfunctions  = 0
    real(wp), allocatable         :: radius(:)
    real(wp), allocatable         :: center(:,:)
    real(wp), allocatable         :: const(:)
    integer, allocatable          :: exponent(:)
    real(wp)                      :: fix_layer   = 0.0_wp
    integer                       :: num_fixatm  = 0
    logical, allocatable          :: fixatm(:)
    ! QMMM
    real(wp), allocatable         :: qm_charge(:)
  end type s_rst

  ! parameters for allocatable variables
  integer,      public, parameter :: RestartAtom      = 1
  integer,      public, parameter :: RestartReplica   = 2
  integer,      public, parameter :: RestartRandom    = 3
  integer,      public, parameter :: RestartRpath     = 4
  integer,      public, parameter :: RestartSpot      = 5
  integer,      public, parameter :: RestartQMMM      = 6

  ! parameters
  integer,      public, parameter :: RstfileTypeUndef = 0
  integer,      public, parameter :: RstfileTypeMin   = 1
  integer,      public, parameter :: RstfileTypeMd    = 2
  integer,      public, parameter :: RstfileTypeRemd  = 3
  integer,      public, parameter :: RstfileTypeRpath = 4
  integer,      public, parameter :: RstfileTypeBd    = 5

  integer,              parameter :: HeaderLength     = 80

  ! local variables
  logical,                private :: vervose = .true.

  ! subroutines
  public  :: input_rst
  public  :: output_rst
  public  :: init_rst
  public  :: alloc_rst
  public  :: dealloc_rst
  public  :: dealloc_rst_all
  private :: read_rst_binary
  private :: read_rst_ascii
  private :: write_rst_binary
  private :: write_rst_ascii

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_rst
  !> @brief        a driver subroutine for reading Restart file
  !! @authors      TM
  !! @param[in]    rst_filename : filename of restart file
  !! @param[out]   rst          : structure of restart information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_rst(rst_filename, rst)

    ! formal arguments
    character(*),            intent(in)    :: rst_filename
    type(s_rst),             intent(inout) :: rst

    ! local variables
    integer                  :: file, endian, ival


    if (get_extension(rst_filename) .ne. 'rsa') then

      ! Open binary restart file
      !

      call open_data(rst_filename, IOFileDataRead, file)

      call read_rst_binary(file, rst)

      call close_data(file)

    else

      ! Open ascii restart file
      !

      call open_file(file, rst_filename, IOFileInput)

      call read_rst_ascii(file, rst)

      call close_file(file)

    end if

    return

  end subroutine input_rst

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_rst
  !> @brief        a driver subroutine for writing Restart file
  !! @authors      TM
  !! @param[in]    rst_filename : filename of restart file
  !! @param[in]    rst          : structure of restart information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_rst(rst_filename, rst)

    ! formal arguments
    character(*),            intent(in) :: rst_filename
    type(s_rst),             intent(in) :: rst

    ! local variables
    integer                  :: file


    if (get_extension(rst_filename) .ne. 'rsa') then

      ! open binary restart file
      !

      call open_data(rst_filename, IOFileDataWrite, file)

      call write_rst_binary(file, rst)

      call close_data(file)

    else

      ! Open ascii restart file
      !

      call open_file(file, rst_filename, IOFileOutputReplace)

      call write_rst_ascii(file, rst)

      call close_file(file)

    end if

    return

  end subroutine output_rst

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_rst
  !> @brief        initialize RST information
  !! @authors      NT
  !! @param[inout] rst : structure of restart information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_rst(rst)

    ! formal arguments
    type(s_rst),             intent(inout) :: rst


    rst%rstfile_type         = 0
    rst%num_atoms            = 0
    rst%iseed                = 0
    rst%thermostat_momentum  = 0.0_wp
    rst%barostat_momentum(3) = 0.0_wp
    rst%box_size_x           = 0.0_wp
    rst%box_size_y           = 0.0_wp
    rst%box_size_z           = 0.0_wp
    rst%energy               = 0.0_wp
    rst%delta_r              = 0.0_wp
    rst%iseed_remd           = 0
    rst%iseed_rpath          = 0
    rst%dimension            = 0
    rst%nreplicas            = 0

    return

  end subroutine init_rst

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_rst
  !> @brief        allocate restart information
  !! @authors      TM, CK, KY
  !! @param[inout] rst       : structure of restart information
  !! @param[in]    variable  : selected variable
  !! @param[in]    var_size  : size of the selected variable
  !! @param[in]    var_size2 : 2nd size of the selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_rst(rst, variable, var_size1, var_size2)

    ! formal arguments
    type(s_rst),             intent(inout) :: rst
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size1
    integer,       optional, intent(in)    :: var_size2

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat


    ! initialize
    !
    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case(RestartAtom)

      if (allocated(rst%coord)) then
        if (size(rst%coord(1,:)) == var_size1) return
        deallocate(rst%coord,    &
                   rst%velocity, &
                   stat = dealloc_stat)
      end if

      allocate(rst%coord(3,var_size1),    &
               rst%velocity(3,var_size1), &
               stat = alloc_stat)

      rst%coord   (1:3,1:var_size1) = 0.0_dp
      rst%velocity(1:3,1:var_size1) = 0.0_dp

    case(RestartReplica)

      if (allocated(rst%num_criteria)) then
        if (size(rst%num_criteria(:,1,1)) == var_size1 .and. &
            size(rst%num_criteria(1,:,1)) == var_size2 ) return
        deallocate(rst%repid2parmsetid, &
                   rst%num_criteria,    &
                   rst%num_exchanges,   &
                   stat = dealloc_stat)
      end if

      allocate(rst%repid2parmsetid(var_size1),              &
               rst%num_criteria (var_size1,var_size2,2),    &
               rst%num_exchanges(var_size1,var_size2,2),    &
               stat = alloc_stat)

      rst%repid2parmsetid(1:var_size1)                = 0
      rst%num_criteria (1:var_size1,1:var_size2,1:2)  = 0
      rst%num_exchanges(1:var_size1,1:var_size2,1:2)  = 0

    case(RestartRpath)

      if (allocated(rst%rest_reference)) then
        if (size(rst%rest_reference(1,:,1)) == var_size1 .and. &
            size(rst%rest_reference(1,1,:)) == var_size2 ) return
        deallocate(rst%rest_reference,  &
                   stat = dealloc_stat)
      end if

      allocate(rst%rest_reference(2,var_size1,var_size2),   &
               stat = alloc_stat)

      rst%rest_reference(1:2,1:var_size1,1:var_size2) = 0_wp

    case (RestartRandom)

      if (allocated(rst%random)) then
        if (size(rst%random) == var_size1) return
        deallocate(rst%random, &
                   stat = dealloc_stat)
      end if

      allocate(rst%random(var_size1), &
               stat = alloc_stat)

      rst%random(1:var_size1) = char(0)

    case (RestartSpot)

      if (allocated(rst%radius)) then
        if (size(rst%radius) == var_size1) return
        deallocate(rst%radius,   &
                   rst%center,   &
                   rst%const,    &
                   rst%exponent, &
                   rst%fixatm,   &
                   stat = dealloc_stat)
      end if

      allocate(rst%radius  (var_size1), &
               rst%center(3,var_size1), &
               rst%const   (var_size1), &
               rst%exponent(var_size1), &
               rst%fixatm  (var_size2), &
               stat = alloc_stat)

      rst%radius  (1:var_size1) = 0.0_wp
      rst%center(:,1:var_size1) = 0.0_wp
      rst%const   (1:var_size1) = 0.0_wp
      rst%exponent(1:var_size1) = 0
      rst%fixatm  (1:var_size2) = .false.

    case (RestartQMMM)

      if (allocated(rst%qm_charge)) then
        if (size(rst%qm_charge) == var_size1) return
        deallocate(rst%qm_charge, stat = dealloc_stat)
      end if

      allocate(rst%qm_charge(var_size1), &
               stat = alloc_stat)

      rst%qm_charge (1:var_size1) = 0.0_wp

    case default

      call error_msg('Alloc_Rst> bad variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_rst

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_rst
  !> @brief        deallocate restart information
  !! @authors      TM, CK, KY
  !! @param[inout] rst      : structure of restart information
  !! @param[in]    variable : selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_rst(rst, variable)

    ! formal arguments
    type(s_rst),             intent(inout) :: rst
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat


    ! initialize
    !
    dealloc_stat = 0

    ! deallocate selected variable
    !
    select case (variable)

    case (RestartAtom)

      if (allocated(rst%coord)) then
        deallocate(rst%coord,    &
                   rst%velocity, &
                   stat = dealloc_stat)
      end if

    case(RestartReplica)

      if (allocated(rst%repid2parmsetid)) then
        deallocate(rst%repid2parmsetid, &
                   rst%num_criteria,    &
                   rst%num_exchanges,   &
                   stat = dealloc_stat)
      end if

    case(RestartRpath)
      if (allocated(rst%rest_reference)) then
        deallocate(rst%rest_reference,  &
                   stat = dealloc_stat)
      end if

    case (RestartRandom)

      if (allocated(rst%random)) then
        deallocate(rst%random, &
                   stat = dealloc_stat)
      end if

    case (RestartSpot)

      if (allocated(rst%radius)) then
        deallocate(rst%radius,   &
                   rst%center,   &
                   rst%const,    &
                   rst%exponent, &
                   rst%fixatm,   &
                   stat = dealloc_stat)
      end if

    case (RestartQMMM)

      if (allocated(rst%qm_charge)) then
        deallocate(rst%qm_charge, &
                   stat = dealloc_stat)
      end if

    case default
      call error_msg('Dealloc_Restart> bad variable')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_rst

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_rst_all
  !> @brief        deallocate all restart information
  !! @authors      TM, CK, KY
  !! @param[inout] rst : structure of restart information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_rst_all(rst)

    ! format arguments
    type(s_rst),             intent(inout) :: rst


    call dealloc_rst(rst, RestartAtom)
    call dealloc_rst(rst, RestartReplica)
    call dealloc_rst(rst, RestartRandom)
    call dealloc_rst(rst, RestartRpath)
    call dealloc_rst(rst, RestartSpot)
    call dealloc_rst(rst, RestartQMMM)

    return

  end subroutine dealloc_rst_all

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_rst_binary
  !> @brief        read binary restart data
  !! @authors      TM, CK, KY
  !! @param[in]    file : unit number of restart file
  !! @param[out]   rst  : structure of restart information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_rst_binary(file, rst)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_rst),             intent(inout) :: rst

    ! local variables
    integer                  :: i, random_size, qm_size
    logical                  :: is_tag

    real(wp),    allocatable :: v(:)


    ! deallocate old data
    !
    call dealloc_rst_all(rst)

    ! read header
    !
    call read_data_integer &
         (file, 'num_atoms', rst%num_atoms)

    ! allocate buffer
    !
    !call alloc_rst(rst, RestartAtom, rst%num_atoms)

    ! restart file type
    !
    call read_data_integer &
         (file, 'rstfile_type', rst%rstfile_type)

    ! MIN
    if (find_tag(file, 'energy')) &
      call read_data_real_wp(file, 'energy', rst%energy)

    if (find_tag(file, 'delta_r')) &
      call read_data_real_wp(file, 'delta_r', rst%delta_r)

    ! MIN and MD
    if (find_tag(file, 'box_size_x')) &
      call read_data_real_wp(file, 'box_size_x', rst%box_size_x)

    if (find_tag(file, 'box_size_y')) &
      call read_data_real_wp(file, 'box_size_y', rst%box_size_y)

    if (find_tag(file, 'box_size_z')) &
      call read_data_real_wp(file, 'box_size_z', rst%box_size_z)

    if (find_tag(file, 'coord_x')) then
      call alloc_rst(rst, RestartAtom, rst%num_atoms)
      allocate(v(rst%num_atoms))

      call read_data_real_wp_array(file, 'coord_x', (/rst%num_atoms/), v(1:rst%num_atoms))
      do i = 1, rst%num_atoms
        rst%coord(1,i) = v(i)
      end do

      call read_data_real_wp_array(file, 'coord_y', (/rst%num_atoms/), v(1:rst%num_atoms))
      do i = 1, rst%num_atoms
        rst%coord(2,i) = v(i)
      end do

      call read_data_real_wp_array(file, 'coord_z', (/rst%num_atoms/), v(1:rst%num_atoms))
      do i = 1, rst%num_atoms
        rst%coord(3,i) = v(i)
      end do

      deallocate(v)
    end if

    ! MD
    if (find_tag(file, 'iseed')) &
      call read_data_integer(file, 'iseed', rst%iseed)

    if (find_tag(file, 'thermostat_momentum')) &
      call read_data_real_wp(file, 'thermostat_momentum', rst%thermostat_momentum)

    if (find_tag(file, 'barostat_momentum')) &
      call read_data_real_wp_array(file, 'barostat_momentum', (/3/), rst%barostat_momentum(1:3))

    if (find_tag(file, 'velocity_x')) then
      call alloc_rst(rst, RestartAtom, rst%num_atoms)
      allocate(v(rst%num_atoms))

      call read_data_real_wp_array(file, 'velocity_x', (/rst%num_atoms/), v(1:rst%num_atoms))
      do i = 1, rst%num_atoms
        rst%velocity(1,i) = v(i)
      end do

      call read_data_real_wp_array(file, 'velocity_y', (/rst%num_atoms/), v(1:rst%num_atoms))
      do i = 1, rst%num_atoms
        rst%velocity(2,i) = v(i)
      end do

      call read_data_real_wp_array(file, 'velocity_z', (/rst%num_atoms/), v(1:rst%num_atoms))
      do i = 1, rst%num_atoms
        rst%velocity(3,i) = v(i)
      end do

      deallocate(v)
    end if

    !call get_data_size(file, 'random', random_size)
    !if (random_size > 0) then
    !  call alloc_rst(rst, RestartRandom, random_size)
    !  call read_data_byte_array &
    !    (file, 'random', (/random_size/), rst%random(1:random_size))
    !end if

    if (find_tag(file, 'random')) then
      call get_data_size(file, 'random', random_size)
      if (random_size > 0) then
        call alloc_rst(rst, RestartRandom, random_size)
        call read_data_byte_array &
          (file, 'random', (/random_size/), rst%random(1:random_size))
      end if
    end if

    ! REMD and RPATH
    if (find_tag(file, 'nreplicas')) &
      call read_data_integer(file, 'nreplicas', rst%nreplicas)

    if (find_tag(file, 'dimension')) &
      call read_data_integer(file, 'dimension', rst%dimension)

    ! REMD
    if (find_tag(file, 'iseed_remd')) &
      call read_data_integer(file, 'iseed_remd', rst%iseed_remd)

    if (find_tag(file, 'repid2parmsetid')) then
      call alloc_rst(rst, RestartReplica, rst%nreplicas, rst%dimension)

      call read_data_integer_array &
        (file, 'repid2parmsetid', (/rst%nreplicas/), rst%repid2parmsetid(:))

      call read_data_integer_array &
        (file, 'num_criteria', (/rst%nreplicas,rst%dimension,2/), &
                                                       rst%num_criteria(:,:,:))
      call read_data_integer_array &
        (file, 'num_exchanges', (/rst%nreplicas,rst%dimension,2/), &
                                                      rst%num_exchanges(:,:,:))
    end if

    ! RPATH
    if (find_tag(file, 'iseed_rpath')) &
      call read_data_integer(file, 'iseed_rpath', rst%iseed_rpath)


    if (find_tag(file, 'rest_reference')) then
      call alloc_rst(rst, RestartRpath, rst%dimension, rst%nreplicas)

      call read_data_real_wp_array &
        (file, 'rest_reference', (/2,rst%dimension,rst%nreplicas/), &
                                                      rst%rest_reference(:,:,:))
    end if

    ! SPOT
    if (find_tag(file, 'sph_pot')) then
      call read_data_logical(file, 'sph_pot', rst%sph_pot)
    else
      rst%sph_pot = .false.
    end if

    if (rst%sph_pot) then
      call read_data_integer(file, 'nfunctions', rst%nfunctions)
      call alloc_rst(rst, RestartSpot, rst%nfunctions, rst%num_atoms)

      call read_data_real_wp_array(file, 'radius', (/rst%nfunctions/), rst%radius)
      call read_data_real_wp_array(file, 'center', (/3,rst%nfunctions/), rst%center)
      call read_data_real_wp_array(file, 'const', (/rst%nfunctions/), rst%const)
      call read_data_integer_array(file, 'exponent', (/rst%nfunctions/), rst%exponent)
      call read_data_real_wp(file, 'fix_layer', rst%fix_layer)
      call read_data_integer(file, 'num_fixatm', rst%num_fixatm)
      call read_data_logical_array(file, 'fixatm', (/rst%num_atoms/), rst%fixatm)
    end if

    ! QMMM
    if (find_tag(file, 'qm_charge')) then
      call get_data_size(file, 'qm_charge', qm_size)
      if (qm_size > 0) then
        call alloc_rst(rst, RestartQMMM, qm_size)
        call read_data_real_wp_array &
          (file, 'qm_charge', (/qm_size/), rst%qm_charge(1:qm_size))
      end if
    end if

    ! write the summary
    !
    if (nrep_per_proc > 1) vervose = .false.
    if (main_rank .and. vervose) then
      write(MsgOut,'(A)') 'Read_Rst_Binary> Summary of RST file'
      !write(MsgOut,'(A)') '  RstfileType     =         MD'
      write(MsgOut,'(A20,I10,A20,I10)')             &
        '  num_atoms       = ', rst%num_atoms,      &
        '  iseed           = ', rst%iseed
      write(MsgOut,'(A20,3F10.3)')                  &
        '  boxsize (x,y,z) = ', rst%box_size_x,     &
                                rst%box_size_y,     &
                                rst%box_size_z

      if (allocated(rst%repid2parmsetid)) then
        !write(MsgOut,'(A)') '  RstfileType     =       REMD'
        write(MsgOut,'(A20,I10,A20,I10)')              &
             '  nreplicas         = ', rst%nreplicas,  &
             '  dimension         = ', rst%dimension
        write(MsgOut,'(A20,I10)') &
             '  remd_iseed        = ', rst%iseed_remd
      end if

      if (allocated(rst%rest_reference)) then
        !write(MsgOut,'(A)') '  RstfileType     =       RPATH'
        write(MsgOut,'(A20,I10,A20,I10)')              &
             '  nreplicas         = ', rst%nreplicas,  &
             '  dimension         = ', rst%dimension
        write(MsgOut,'(A20,I10)') &
             '  rpath iseed       = ', rst%iseed_rpath
      end if

      if (rst%sph_pot) then
        write(MsgOut,'(A30)') &
             '  spherical_pot     =      yes'
        write(MsgOut,'(A20,I10,A20,I10)')              &
             '  nfunctions        = ', rst%nfunctions, &
             '  num_fixatm        = ', rst%num_fixatm
      end if

      if (allocated(rst%qm_charge)) then
        write(MsgOut,'(A20,I10,A20,I10)')              &
             '  num_of_qm_charge  = ', size(rst%qm_charge)
      end if
      write(MsgOut,'(A)') ''
      vervose = .false.
    end if

    return

  end subroutine read_rst_binary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_rst_ascii
  !> @brief        read ascii restart data
  !! @authors      TM, CK, KY
  !! @param[in]    file : unit number of restart file
  !! @param[out]   rst  : structure of restart information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_rst_ascii(file, rst)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_rst),             intent(inout) :: rst

    ! local variables
    integer                  :: i, random_size, qm_size
    character(HeaderLength)  :: head1, head2
    character(100)           :: tag

    integer, allocatable     :: irnd(:)


    ! deallocate old data
    !
    call dealloc_rst_all(rst)

    ! read header
    !
    rewind(file)

    read(file,*) head1
    read(file,*) head2
    read(file,'(a)') tag
    read(file,*) rst%num_atoms

    ! allocate buffer
    !
    call alloc_rst(rst, RestartAtom, rst%num_atoms)

    ! select restart file type
    !
    read(file,'(a)') tag
    read(file,*) rst%rstfile_type

    select case(rst%rstfile_type)

    case(RstfileTypeMin)

      read(file,'(a)') tag
      read(file,*) rst%box_size_x, rst%box_size_y, rst%box_size_z
      read(file,'(a)') tag
      read(file,*) rst%energy, rst%delta_r
      read(file,'(a)') tag
      read(file,*) rst%coord(1,:)
      read(file,'(a)') tag
      read(file,*) rst%coord(2,:)
      read(file,'(a)') tag
      read(file,*) rst%coord(3,:)

      ! write the summary
      !
      if (nrep_per_proc > 1) vervose = .false.
      if (main_rank .and. vervose) then

        write(MsgOut,'(A)') 'Read_Rst_Ascii> Summary of RST file'

        write(MsgOut,'(A)') '  RstfileType     =        MIN'
        write(MsgOut,'(A20,I10,A20,I10)')                &
             '  num_atoms       = ', rst%num_atoms
        write(MsgOut,'(A20,3F10.3)')                     &
             '  boxsize (x,y,z) = ', rst%box_size_x,     &
                                     rst%box_size_y,     &
                                     rst%box_size_z

        write(MsgOut,'(A)') ' '
        vervose = .false.
      end if

    case(RstfileTypeMd)

      read(file,'(a)') tag
      read(file,*) rst%iseed
      read(file,'(a)') tag
      read(file,*) rst%box_size_x, rst%box_size_y, rst%box_size_z
      read(file,'(a)') tag
      read(file,*) rst%thermostat_momentum
      read(file,'(a)') tag
      read(file,*) rst%barostat_momentum(1:3)
      read(file,'(a)') tag
      read(file,*) rst%coord(1,:)
      read(file,'(a)') tag
      read(file,*) rst%coord(2,:)
      read(file,'(a)') tag
      read(file,*) rst%coord(3,:)
      read(file,'(a)') tag
      read(file,*) rst%velocity(1,:)
      read(file,'(a)') tag
      read(file,*) rst%velocity(2,:)
      read(file,'(a)') tag
      read(file,*) rst%velocity(3,:)

      random_size = 0
      read(file,'(a)') tag
      read(file,*,end=10) random_size
10    if (random_size > 0) then
        call alloc_rst(rst, RestartRandom, random_size)
        allocate(irnd(random_size))
        read(file,*) irnd(1:random_size)
        rst%random(1:random_size) = char(irnd(1:random_size))
        deallocate(irnd)
      end if

      ! write the summary
      !
      if (main_rank) then
        write(MsgOut,'(A)') 'Read_Rst_Ascii> Summary of RST file'

        write(MsgOut,'(A)') '  RstfileType     =         MD'
        write(MsgOut,'(A20,I10,A20,I10)')                &
             '  num_atoms       = ', rst%num_atoms,      &
             '  iseed           = ', rst%iseed
        write(MsgOut,'(A20,3F10.3)')                     &
             '  boxsize (x,y,z) = ', rst%box_size_x,     &
                                     rst%box_size_y,     &
                                     rst%box_size_z
        write(MsgOut,'(A)') ''

      end if

    case(RstfileTypeRemd)

      read(file,'(a)') tag
      read(file,*) rst%iseed
      read(file,'(a)') tag
      read(file,*) rst%box_size_x, rst%box_size_y, rst%box_size_z
      read(file,'(a)') tag
      read(file,*) rst%thermostat_momentum
      read(file,'(a)') tag
      read(file,*) rst%barostat_momentum(1:3)
      read(file,'(a)') tag
      read(file,*) rst%coord(1,:)
      read(file,'(a)') tag
      read(file,*) rst%coord(2,:)
      read(file,'(a)') tag
      read(file,*) rst%coord(3,:)
      read(file,'(a)') tag
      read(file,*) rst%velocity(1,:)
      read(file,'(a)') tag
      read(file,*) rst%velocity(2,:)
      read(file,'(a)') tag
      read(file,*) rst%velocity(3,:)

      read(file,'(a)') tag
      read(file,*) rst%iseed_remd
      read(file,'(a)') tag
      read(file,*) rst%nreplicas
      read(file,'(a)') tag
      read(file,*) rst%dimension

      call alloc_rst(rst, RestartReplica, rst%nreplicas, rst%dimension)

      read(file,'(a)') tag
      read(file,*) rst%repid2parmsetid(:)
      do i = 1, rst%nreplicas
        read(file,'(a)') tag
        read(file,*) rst%num_criteria (i,:,1)
        read(file,*) rst%num_criteria (i,:,2)
        read(file,'(a)') tag
        read(file,*) rst%num_exchanges(i,:,1)
        read(file,*) rst%num_exchanges(i,:,2)
      end do

      random_size = 0
      read(file,'(a)') tag
      read(file,*,end=20) random_size
20    if (random_size > 0) then
        call alloc_rst(rst, RestartRandom, random_size)
        allocate(irnd(random_size))
        read(file,*) irnd(1:random_size)
        rst%random(1:random_size) = char(irnd(1:random_size))
        deallocate(irnd)
      end if

      if (main_rank) then
        write(MsgOut,'(A)') 'Read_Rst_Ascii> Summary of RST file'

        write(MsgOut,'(A)') '  RstfileType     =       REMD'
        write(MsgOut,'(A20,I10,A20,I10)')              &
             '  num_atoms       = ', rst%num_atoms,    &
             '  iseed           = ', rst%iseed
        write(MsgOut,'(A20,3F10.3)')                   &
             '  boxsize (x,y,z) = ', rst%box_size_x,   &
                                     rst%box_size_y,   &
                                     rst%box_size_z
        write(MsgOut,'(A20,I10,A20,I10)')              &
             '  nreplicas         = ', rst%nreplicas,  &
             '  dimension         = ', rst%dimension
        write(MsgOut,'(A20,I10)') &
             '  remd iseed        = ', rst%iseed_remd
      end if

    case(RstfileTypeRpath)

      read(file,'(a)') tag
      read(file,*) rst%iseed
      read(file,'(a)') tag
      read(file,*) rst%box_size_x, rst%box_size_y, rst%box_size_z
      read(file,'(a)') tag
      read(file,*) rst%thermostat_momentum
      read(file,'(a)') tag
      read(file,*) rst%barostat_momentum(1:3)
      read(file,'(a)') tag
      read(file,*) rst%coord(1,:)
      read(file,'(a)') tag
      read(file,*) rst%coord(2,:)
      read(file,'(a)') tag
      read(file,*) rst%coord(3,:)
      read(file,'(a)') tag
      read(file,*) rst%velocity(1,:)
      read(file,'(a)') tag
      read(file,*) rst%velocity(2,:)
      read(file,'(a)') tag
      read(file,*) rst%velocity(3,:)

      read(file,'(a)') tag
      read(file,*) rst%iseed_rpath
      read(file,'(a)') tag
      read(file,*) rst%nreplicas
      read(file,'(a)') tag
      read(file,*) rst%dimension

      call alloc_rst(rst, RestartRpath, rst%dimension, rst%nreplicas)

      do i = 1, rst%dimension
        read(file,'(a)') tag
        read(file,*) rst%rest_reference (1,i,:)
        read(file,*) rst%rest_reference (2,i,:)
      end do

      random_size = 0
      read(file,'(a)') tag
      read(file,*,end=30) random_size
30    if (random_size > 0) then
        call alloc_rst(rst, RestartRandom, random_size)
        allocate(irnd(random_size))
        read(file,*) irnd(1:random_size)
        rst%random(1:random_size) = char(irnd(1:random_size))
        deallocate(irnd)
      end if

      if (main_rank) then
        write(MsgOut,'(A)') 'Read_Rst_Ascii> Summary of RST file'

        write(MsgOut,'(A)') '  RstfileType     =       RPATH'
        write(MsgOut,'(A20,I10,A20,I10)')              &
             '  num_atoms       = ', rst%num_atoms,    &
             '  iseed           = ', rst%iseed
        write(MsgOut,'(A20,3F10.3)')                   &
             '  boxsize (x,y,z) = ', rst%box_size_x,   &
                                     rst%box_size_y,   &
                                     rst%box_size_z
        write(MsgOut,'(A20,I10,A20,I10)')              &
             '  nreplicas         = ', rst%nreplicas,  &
             '  dimension         = ', rst%dimension
        write(MsgOut,'(A20,I10)') &
             '  rpath iseed       = ', rst%iseed_rpath
      end if

    end select

    read(file,'(a)', end=999) tag
    read(file,*) rst%sph_pot
    if (rst%sph_pot) then
      read(file,'(a)') tag
      read(file,*) rst%nfunctions
      call alloc_rst(rst, RestartSpot, rst%nfunctions, rst%num_atoms)

      read(file,'(a)') tag
      read(file,*) rst%radius
      read(file,'(a)') tag
      read(file,*) rst%center
      read(file,'(a)') tag
      read(file,*) rst%const
      read(file,'(a)') tag
      read(file,*) rst%exponent
      read(file,'(a)') tag
      read(file,*) rst%fix_layer
      read(file,'(a)') tag
      read(file,*) rst%num_fixatm
      read(file,'(a)') tag
      read(file,*) rst%fixatm
    end if

    if (main_rank .and. rst%sph_pot) then
      write(MsgOut,'(A30)') &
           '  spherical_pot     =      yes'
      write(MsgOut,'(A20,I10,A20,I10)')              &
           '  nfunctions        = ', rst%nfunctions, &
           '  num_fixatm        = ', rst%num_fixatm
    end if

    read(file,'(a)', end=999) tag
    read(file,*,end=40) qm_size
40  if (qm_size > 0) then
      call alloc_rst(rst, RestartQMMM, qm_size)
      read(file,*)  rst%qm_charge
    end if

    if (main_rank .and. qm_size > 0) then
      write(MsgOut,'(A20,I10,A20,I10)')              &
           '  num_of_qm_charge  = ', size(rst%qm_charge)
    end if

999 continue

    return

  end subroutine read_rst_ascii

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_rst_binary
  !> @brief        write data to binary restart file
  !! @authors      TM, CK
  !! @param[in]    file : unit number of restart file
  !! @param[in]    rst  : structure of restart information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_rst_binary(file, rst)

    ! formal arguments
    integer,                 intent(in) :: file
    type(s_rst),             intent(in) :: rst

    ! local variables
    integer                  :: i
    character(HeaderLength)  :: head1, head2
    character(24)            :: name, date


    ! check memory allocation
    !
    if (.not. allocated(rst%coord))    &
      call error_msg('Write_Rst_Binary> not allocated: rst%coord')
    if (.not. allocated(rst%velocity)) &
      call error_msg('Write_Rst_Binary> not allocated: rst%velocity')

    ! make header
    !
    call fdate(date)
    call getlog(name)

    head1 = 'REMARKS CREATED BY GENESIS                                                      '
    head2 = 'REMARKS DATE: ' // date // ' CREATED BY USER: ' // name

    call write_data_byte_array &
         (file, 'header1', (/HeaderLength/), head1)
    call write_data_byte_array &
         (file, 'header2', (/HeaderLength/), head2)
    call write_data_integer &
         (file, 'num_atoms', rst%num_atoms)

    ! restart file type
    !
    call write_data_integer(file, 'rstfile_type', rst%rstfile_type)

    ! MIN
    call write_data_real_wp(file, 'energy', rst%energy)
    call write_data_real_wp(file, 'delta_r', rst%delta_r)

    ! MIN and MD
    call write_data_real_wp(file, 'box_size_x', rst%box_size_x)
    call write_data_real_wp(file, 'box_size_y', rst%box_size_y)
    call write_data_real_wp(file, 'box_size_z', rst%box_size_z)

    if (allocated(rst%coord)) then
      call write_data_real_wp_array(file, 'coord_x', (/rst%num_atoms/), rst%coord(1,:))
      call write_data_real_wp_array(file, 'coord_y', (/rst%num_atoms/), rst%coord(2,:))
      call write_data_real_wp_array(file, 'coord_z', (/rst%num_atoms/), rst%coord(3,:))
    end if

    ! MD
    call write_data_integer(file, 'iseed', rst%iseed)
    call write_data_real_wp(file, 'thermostat_momentum', rst%thermostat_momentum)
    call write_data_real_wp_array(file, 'barostat_momentum', (/3/), rst%barostat_momentum(1:3))

    if (allocated(rst%velocity)) then
      call write_data_real_wp_array(file, 'velocity_x', (/rst%num_atoms/), rst%velocity(1,:))
      call write_data_real_wp_array(file, 'velocity_y', (/rst%num_atoms/), rst%velocity(2,:))
      call write_data_real_wp_array(file, 'velocity_z', (/rst%num_atoms/), rst%velocity(3,:))
    end if

    if (allocated(rst%random)) then
      call write_data_byte_array(file, 'random',(/size(rst%random)/), rst%random(1:size(rst%random)))
    !else
    !  call write_data_byte_array(file, 'random', (/0/), rst%random)
    end if

    ! REMD and RPATH
    call write_data_integer(file, 'nreplicas', rst%nreplicas)
    call write_data_integer(file, 'dimension', rst%dimension)

    ! REMD
    call write_data_integer(file, 'iseed_remd', rst%iseed_remd)

    if (allocated(rst%repid2parmsetid)) then
      call write_data_integer_array(file, 'repid2parmsetid', (/rst%nreplicas/), rst%repid2parmsetid(:))
      call write_data_integer_array(file, 'num_criteria', (/rst%nreplicas,rst%dimension,2/), rst%num_criteria(:,:,:))
      call write_data_integer_array(file, 'num_exchanges', (/rst%nreplicas,rst%dimension,2/), rst%num_exchanges(:,:,:))
    end if

    ! RPATH
    call write_data_integer(file, 'iseed_rpath', rst%iseed_rpath)

    if (allocated(rst%rest_reference)) then
      call write_data_real_wp_array(file, 'rest_reference', (/2,rst%dimension,rst%nreplicas/), rst%rest_reference)
    end if

    ! SPOT
    call write_data_logical(file, 'sph_pot', rst%sph_pot)
    if (rst%sph_pot) then
      call write_data_integer(file, 'nfunctions', rst%nfunctions)
      call write_data_real_wp_array(file, 'radius', (/rst%nfunctions/), rst%radius)
      call write_data_real_wp_array(file, 'center', (/3,rst%nfunctions/), rst%center)
      call write_data_real_wp_array(file, 'const', (/rst%nfunctions/), rst%const)
      call write_data_integer_array(file, 'exponent', (/rst%nfunctions/), rst%exponent)
      call write_data_real_wp(file, 'fix_layer', rst%fix_layer)
      call write_data_integer(file, 'num_fixatm', rst%num_fixatm)
      call write_data_logical_array(file, 'fixatm', (/rst%num_atoms/), rst%fixatm)
    end if

    ! QMMM
    if (allocated(rst%qm_charge)) then
      call write_data_real_wp_array(file, 'qm_charge', (/size(rst%qm_charge)/), &
                                 rst%qm_charge)
    !else
    !  call write_data_real_wp_array(file, 'qm_charge', (/0/), rst%qm_charge)
    end if

    return

  end subroutine write_rst_binary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_rst_ascii
  !> @brief        write data to ascii restart file
  !! @authors      TM, CK
  !! @param[in]    file : unit number of restart file
  !! @param[in]    rst  : structure of restart information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_rst_ascii(file, rst)

    ! formal arguments
    integer,                 intent(in) :: file
    type(s_rst),             intent(in) :: rst

    ! local variables
    integer                  :: i, dummy
    character(HeaderLength)  :: head1, head2
    character(24)            :: name, date


    ! check memory allocation
    !
    if (.not. allocated(rst%coord))    &
      call error_msg('Write_Rst_Ascii> not allocated: rst%coord')
    if (.not. allocated(rst%velocity)) &
      call error_msg('Write_Rst_Ascii> not allocated: rst%velocity')

    ! make header
    !
    call fdate(date)
    call getlog(name)

    head1 = 'REMARKS CREATED BY GENESIS                                                      '
    head2 = 'REMARKS DATE: ' // date // ' CREATED BY USER: ' // name

    write(file,*) trim(head1)
    write(file,*) trim(head2)
    write(file,*) '# num_atoms'
    write(file,*) rst%num_atoms

    ! select restart file type
    !
    write(file,*) '# rstfile_type'
    write(file,*) rst%rstfile_type

    select case(rst%rstfile_type)

    case(RstfileTypeMin)

      write(file,*) '# box_size'
      write(file,*) rst%box_size_x, rst%box_size_y, rst%box_size_z
      write(file,*) '# energy, delta_r'
      write(file,*) rst%energy, rst%delta_r
      write(file,*) '# coord_x'
      write(file,*) rst%coord(1,:)
      write(file,*) '# coord_y'
      write(file,*) rst%coord(2,:)
      write(file,*) '# coord_z'
      write(file,*) rst%coord(3,:)

    case(RstfileTypeMd)

      write(file,*) '# iseed'
      write(file,*) rst%iseed
      write(file,*) '# box_size'
      write(file,*) rst%box_size_x, rst%box_size_y, rst%box_size_z
      write(file,*) '# thermostat_momentum'
      write(file,*) rst%thermostat_momentum
      write(file,*) '# barostat_momentum'
      write(file,*) rst%barostat_momentum(1:3)
      write(file,*) '# coord_x'
      write(file,*) rst%coord(1,:)
      write(file,*) '# coord_y'
      write(file,*) rst%coord(2,:)
      write(file,*) '# coord_z'
      write(file,*) rst%coord(3,:)
      write(file,*) '# velocity_x'
      write(file,*) rst%velocity(1,:)
      write(file,*) '# velocity_y'
      write(file,*) rst%velocity(2,:)
      write(file,*) '# velocity_z'
      write(file,*) rst%velocity(3,:)

      write(file,*) '# random'
      if (allocated(rst%random)) then
        write(file,*) size(rst%random)
        write(file,*) (transfer(rst%random(i),dummy),i=1,size(rst%random))
      else
        write(file,*) 0
      end if

    case(RstfileTypeRemd)

      write(file,*) '# iseed'
      write(file,*) rst%iseed
      write(file,*) '# box_size'
      write(file,*) rst%box_size_x, rst%box_size_y, rst%box_size_z
      write(file,*) '# thermostat_momentum'
      write(file,*) rst%thermostat_momentum
      write(file,*) '# barostat_momentum'
      write(file,*) rst%barostat_momentum(1:3)
      write(file,*) '# coord_x'
      write(file,*) rst%coord(1,:)
      write(file,*) '# coord_y'
      write(file,*) rst%coord(2,:)
      write(file,*) '# coord_z'
      write(file,*) rst%coord(3,:)
      write(file,*) '# velocity_x'
      write(file,*) rst%velocity(1,:)
      write(file,*) '# velocity_y'
      write(file,*) rst%velocity(2,:)
      write(file,*) '# velocity_z'
      write(file,*) rst%velocity(3,:)
      write(file,*) '# iseed_remd'
      write(file,*) rst%iseed_remd
      write(file,*) '# nreplicas'
      write(file,*) rst%nreplicas
      write(file,*) '# dimension'
      write(file,*) rst%dimension
      write(file,*) '# repid2parmsetid'
      write(file,*) rst%repid2parmsetid(:)
      do i = 1, rst%nreplicas
        write(file,*) '# num_criteria', i
        write(file,*) rst%num_criteria (i,:,1)
        write(file,*) rst%num_criteria (i,:,2)
        write(file,*) '# num_exchanges', i
        write(file,*) rst%num_exchanges(i,:,1)
        write(file,*) rst%num_exchanges(i,:,2)
      end do

      write(file,*) '# random'
      if (allocated(rst%random)) then
        write(file,*) size(rst%random)
        write(file,*) (transfer(rst%random(i),dummy),i=1,size(rst%random))
      else
        write(file) 0
      end if

    case(RstfileTypeRpath)

      write(file,*) '# iseed'
      write(file,*) rst%iseed
      write(file,*) '# box_size'
      write(file,*) rst%box_size_x, rst%box_size_y, rst%box_size_z
      write(file,*) '# thermostat_momentum'
      write(file,*) rst%thermostat_momentum
      write(file,*) '# barostat_momentum'
      write(file,*) rst%barostat_momentum(1:3)
      write(file,*) '# coord_x'
      write(file,*) rst%coord(1,:)
      write(file,*) '# coord_y'
      write(file,*) rst%coord(2,:)
      write(file,*) '# coord_z'
      write(file,*) rst%coord(3,:)
      write(file,*) '# velocity_x'
      write(file,*) rst%velocity(1,:)
      write(file,*) '# velocity_y'
      write(file,*) rst%velocity(2,:)
      write(file,*) '# velocity_z'
      write(file,*) rst%velocity(3,:)
      write(file,*) '# iseed_rpath'
      write(file,*) rst%iseed_rpath
      write(file,*) '# nreplicas'
      write(file,*) rst%nreplicas
      write(file,*) '# dimension'
      write(file,*) rst%dimension
      do i = 1, rst%dimension
        write(file,*) '# num_replica', i
        write(file,*) rst%rest_reference (1,i,:)
        write(file,*) rst%rest_reference (2,i,:)
      end do

      write(file,*) '# random'
      if (allocated(rst%random)) then
        write(file,*) size(rst%random)
        write(file,*) (transfer(rst%random(i),dummy),i=1,size(rst%random))
      else
        write(file) 0
      end if

    end select

    write(file,*) '# sph_pot'
    write(file,*) rst%sph_pot
    if (rst%sph_pot) then
      write(file,*) '# nfunctions'
      write(file,*) rst%nfunctions
      write(file,*) '# radius'
      write(file,*) rst%radius
      write(file,*) '# center'
      write(file,*) rst%center
      write(file,*) '# const'
      write(file,*) rst%const
      write(file,*) '# exponent'
      write(file,*) rst%exponent
      write(file,*) '# fix_layer'
      write(file,*) rst%fix_layer
      write(file,*) '# num_fixatm'
      write(file,*) rst%num_fixatm
      write(file,*) '# fixatm'
      write(file,*) rst%fixatm
    end if

    write(file,*) '# qm_charge'
    if (allocated(rst%qm_charge)) then
      write(file,*) size(rst%qm_charge)
      write(file,*) rst%qm_charge
    else
      write(file,*) 0
    end if

    return

  end subroutine write_rst_ascii

end module fileio_rst_mod
