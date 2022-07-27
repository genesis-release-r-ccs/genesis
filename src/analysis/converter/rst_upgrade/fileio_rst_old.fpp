!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_rst_old
!> @brief   Restart file I/O module
!! @authors Takaharu Mori (TM), Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module fileio_rst_old_mod

  use fileio_rst_mod
  use fileio_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  integer,              parameter :: HeaderLength     = 80

  ! subroutines
  public  :: input_rst_old
  private :: read_rst_binary
  private :: read_rst_ascii
  private :: check_endian

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_rst_old
  !> @brief        a driver subroutine for reading Restart file
  !! @authors      TM
  !! @param[in]    rst_filename : filename of restart file
  !! @param[out]   rst          : structure of restart information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_rst_old(rst_filename, rst)

    ! formal arguments
    character(*),            intent(in)    :: rst_filename
    type(s_rst),             intent(inout) :: rst

    ! local variables
    integer                  :: file, endian, ival


    if (get_extension(rst_filename) /= 'rsa') then

      ! Open binary restart file
      !

      call open_binary_file(file, rst_filename, IOFileInput, &
                            check_endian(rst_filename))

      call read_rst_binary(file, rst)

      call close_file(file)

    else

      ! Open ascii restart file
      !

      call open_file(file, rst_filename, IOFileInput)

      call read_rst_ascii(file, rst)

      call close_file(file)

    end if

    return

  end subroutine input_rst_old

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_rst_binary
  !> @brief        read binary restart data
  !! @authors      TM, CK
  !! @param[in]    file : unit number of restart file
  !! @param[out]   rst  : structure of restart information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_rst_binary(file, rst)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_rst),             intent(inout) :: rst

    ! local variables
    integer                  :: i, random_size
    integer                  :: num_deg_freedom
    character(HeaderLength)  :: head1, head2


    ! deallocate old data
    !
    call dealloc_rst_all(rst)

    ! read header
    !
    rewind(file)

    read(file) head1
    read(file) head2
    read(file) rst%num_atoms

    ! allocate buffer
    !
    call alloc_rst(rst, RestartAtom, rst%num_atoms)

    ! select restart file type
    !
    read(file) rst%rstfile_type

    select case(rst%rstfile_type)

    case(RstfileTypeMin)

      read(file) rst%box_size_x, rst%box_size_y, rst%box_size_z
      read(file) rst%energy, rst%delta_r
      read(file) rst%coord(1,:)
      read(file) rst%coord(2,:)
      read(file) rst%coord(3,:)

      ! write the summary
      !
      if (main_rank) then

        write(MsgOut,'(A)') 'Read_Rst_Binary> Summary of RST file'

        write(MsgOut,'(A)') '  RstfileType     =        MIN'
        write(MsgOut,'(A20,I10,A20,I10)')                &
             '  num_atoms       = ', rst%num_atoms
        write(MsgOut,'(A20,3F10.3)')                     &
             '  boxsize (x,y,z) = ', rst%box_size_x,     &
                                     rst%box_size_y,     &
                                     rst%box_size_z

        write(MsgOut,'(A)') ' '
      end if

    case(RstfileTypeMd)

      read(file) rst%iseed, num_deg_freedom
      read(file) rst%box_size_x, rst%box_size_y, rst%box_size_z
      read(file) rst%thermostat_momentum
      read(file) rst%barostat_momentum(1:3)
      read(file) rst%coord(1,:)
      read(file) rst%coord(2,:)
      read(file) rst%coord(3,:)
      read(file) rst%velocity(1,:)
      read(file) rst%velocity(2,:)
      read(file) rst%velocity(3,:)

      random_size = 0
      read(file,end=10) random_size
10    if (random_size > 0) then
        call alloc_rst(rst, RestartRandom, random_size)
        read(file) rst%random(1:random_size)
      end if

      ! write the summary
      !
      if (main_rank) then
        write(MsgOut,'(A)') 'Read_Rst_Binary> Summary of RST file'

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

      read(file) rst%iseed, num_deg_freedom
      read(file) rst%box_size_x, rst%box_size_y, rst%box_size_z
      read(file) rst%thermostat_momentum
      read(file) rst%barostat_momentum(1:3)
      read(file) rst%coord(1,:)
      read(file) rst%coord(2,:)
      read(file) rst%coord(3,:)
      read(file) rst%velocity(1,:)
      read(file) rst%velocity(2,:)
      read(file) rst%velocity(3,:)

      read(file) rst%iseed_remd
      read(file) rst%nreplicas
      read(file) rst%dimension

      call alloc_rst(rst, RestartReplica, rst%nreplicas, rst%dimension)

      read(file) rst%repid2parmsetid(:)
      do i = 1, rst%nreplicas
        read(file) rst%num_criteria (i,:,1)
        read(file) rst%num_criteria (i,:,2)
        read(file) rst%num_exchanges(i,:,1)
        read(file) rst%num_exchanges(i,:,2)
      end do

      random_size = 0
      read(file,end=20) random_size
20    if (random_size > 0) then
        call alloc_rst(rst, RestartRandom, random_size)
        read(file) rst%random(1:random_size)
      end if

      if (main_rank) then
        write(MsgOut,'(A)') 'Read_Rst_Binary> Summary of RST file'

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
    end select

    return

  end subroutine read_rst_binary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_rst_ascii
  !> @brief        read ascii restart data
  !! @authors      TM, CK
  !! @param[in]    file : unit number of restart file
  !! @param[out]   rst  : structure of restart information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_rst_ascii(file, rst)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_rst),             intent(inout) :: rst

    ! local variables
    integer                  :: i, random_size
    integer                  :: num_deg_freedom
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
      if (main_rank) then

        write(MsgOut,'(A)') 'Read_Rst_Ascii> Summary of RST file'

        write(MsgOut,'(A)') '  RstfileType     =        MIN'
        write(MsgOut,'(A20,I10,A20,I10)')                &
             '  num_atoms       = ', rst%num_atoms
        write(MsgOut,'(A20,3F10.3)')                     &
             '  boxsize (x,y,z) = ', rst%box_size_x,     &
                                     rst%box_size_y,     &
                                     rst%box_size_z

        write(MsgOut,'(A)') ' '
      end if

    case(RstfileTypeMd)

      read(file,'(a)') tag
      read(file,*) rst%iseed, num_deg_freedom
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
      read(file,*) rst%iseed, num_deg_freedom
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

      read(file) rst%iseed, num_deg_freedom
      read(file) rst%box_size_x, rst%box_size_y, rst%box_size_z
      read(file) rst%thermostat_momentum
      read(file) rst%barostat_momentum(1:3)
      read(file) rst%coord(1,:)
      read(file) rst%coord(2,:)
      read(file) rst%coord(3,:)
      read(file) rst%velocity(1,:)
      read(file) rst%velocity(2,:)
      read(file) rst%velocity(3,:)

      read(file) rst%iseed_rpath
      read(file) rst%nreplicas
      read(file) rst%dimension

      call alloc_rst(rst, RestartRpath, rst%nreplicas, rst%dimension)

      read(file) rst%repid2parmsetid(:)
      do i = 1, rst%nreplicas
        read(file) rst%num_criteria (i,:,1)
        read(file) rst%num_criteria (i,:,2)
        read(file) rst%num_exchanges(i,:,1)
        read(file) rst%num_exchanges(i,:,2)
      end do

      if (main_rank) then
        write(MsgOut,'(A)') 'Read_Rst> Summary of RST file'

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
             '  rpath iseed        = ', rst%iseed_rpath
      end if

    end select

    return

  end subroutine read_rst_ascii

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      check_endian
  !> @brief        check restart file endian
  !! @authors      NT
  !! @param[in]    rst_filename : filename of restart file
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function check_endian(rst_filename)

    ! return value
    integer                  :: check_endian

    ! formal arguments
    character(*),            intent(in)    :: rst_filename

    ! local variables
    integer                  :: file, ival


    file = get_unit_no()

    open(file, &
         file    = rst_filename,    &
         status  = 'old',           &
         form    = 'unformatted',   &
         convert = 'little_endian', &
         access  = 'direct',        &
         recl    = 4)
    read(file,rec=1) ival
    close(file)

    if (ival /= HeaderLength) then

      open(file, &
           file    = rst_filename,  &
           status  = 'old',         &
           form    = 'unformatted', &
           convert = 'big_endian',  &
           access  = 'direct',      &
           recl    = 4)
      read(file,rec=1) ival
      close(file)

      if (ival /= HeaderLength) &
        call error_msg('Input_Rst> unknown file format.')

      check_endian = IOFileBigEndian

    else

      check_endian = IOFileLittleEndian

    end if

    call free_unit_no(file)

    return

  end function check_endian

end module fileio_rst_old_mod
