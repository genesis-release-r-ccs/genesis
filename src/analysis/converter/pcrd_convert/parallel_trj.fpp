!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   parallel_trj_mod
!> @brief   module for read parallel trajectory files
!! @authors Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module parallel_trj_mod

  use trajectory_str_mod
  use fileio_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef OMP
  use omp_lib
#endif

  implicit none
  private

  ! parameters
  integer, parameter :: MaxPTInfo = 10
  integer, parameter :: MaxPlace  = 8

  ! structures
  type, public :: s_pt_file
    character(MaxLine)                    :: name
    integer                               :: pos
    integer                               :: endian
  end type s_pt_file

  type, public :: s_pt_buf
    integer                               :: natom
    integer,                  allocatable :: idx(:)
    real(sp),                 allocatable :: crd(:,:)
    real(sp)                              :: box(6)
  end type s_pt_buf

  type, public :: s_pt_info
    type(s_pt_file),          allocatable :: files(:)
    type(s_pt_buf),           allocatable :: bufs(:)
  end type s_pt_info

  type(s_pt_info),                 target :: g_pt_info(MaxPtInfo)

  ! subroutines
  public  :: open_parallel_trj
  public  :: close_parallel_trj
  public  :: read_parallel_trj
  private :: read_parallel_trj_
  private :: get_ranked_filename
  private :: get_file_endian

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    open_parallel_trj
  !> @brief        open parallel trajectory file
  !! @authors      NT
  !! @param[inout] handle   : handle for parallel trajectory
  !! @param[in]    filename : file name for parallel trajectory
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine open_parallel_trj(handle, filename)
  
    ! formal argments
    type(s_trj_file),        intent(inout) :: handle
    character(*),            intent(in)    :: filename

    ! local variables
    integer                  :: i, nfiles, nthreads, nplace
    logical                  :: bex

    type(s_pt_file), pointer :: file
    type(s_pt_buf),  pointer :: buf


    ! check opened
    do i = 1, MaxPTInfo
      if (.not. allocated(g_pt_info(i)%files)) &
        exit
    end do

    if (i > MaxPtInfo) &
      call error_msg('Open_Parallel_Trj> Open error.')

    ! check # of place
    nfiles = 0
    nplace = 1

    do i = 1, MaxPlace

      inquire(file  = get_ranked_filename(filename, nfiles, nplace), &
              exist = bex)
      if (bex) &
        exit
      nplace = nplace + 1

    end do

    if (nplace > MaxPlace) &
      call error_msg('Open_Parallel_Trj> ERROR: Rank0 File not found:'// &
      trim(filename))

    ! check # of parallel trajectries
    do while(.true.)
      inquire(file  = get_ranked_filename(filename, nfiles, nplace), &
              exist = bex)
      if (.not. bex) &
        exit
      nfiles = nfiles + 1
    end do

    ! check # of threads
#ifdef OMP
    nthreads = omp_get_max_threads()
#else
    nthreads = 1
#endif

    ! setup handle
    handle%name       = filename
    handle%trj_format = 0
    handle%trj_type   = 0
    handle%in_out     = IOFileInput
    handle%unit_no    = i
    handle%byte_swap  = .false.
    handle%rw_header  = .false.
    handle%step_no    = 0

    ! setup internal informations
    allocate(g_pt_info(handle%unit_no)%files(nfiles))
    
    do i = 1, nfiles
      file => g_pt_info(handle%unit_no)%files(i)
      file%name   = get_ranked_filename(filename, i-1, nplace)
      file%pos    = 0
      file%endian = get_file_endian(file%name)
    end do

    allocate(g_pt_info(handle%unit_no)%bufs(nthreads))

    do i = 1, nthreads
      buf  => g_pt_info(handle%unit_no)%bufs(i)
      buf%natom = 0
    end do

    write(MsgOut,'(a,a,a,I0,a,I0,a)') &
                        'Open_Parallel_Trj> File:', trim(filename), &
                        '  # of parallel:', nfiles, &
                        '  place(', nplace, ')'
    write(MsgOut,'(a)') ' '

    return

  end subroutine open_parallel_trj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    close_parallel_trj
  !> @brief        close parallel trajectory file
  !! @authors      NT
  !! @param[in]    handle : handle for parallel trajectory
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine close_parallel_trj(handle)
  
    ! formal argments
    type(s_trj_file),        intent(in) :: handle

    ! local variables
    integer                  :: i

    type(s_pt_buf),  pointer :: buf


    if (handle%unit_no == 0 .or. &
        handle%unit_no > MaxPTInfo) &
      call error_msg('Close_Parallel_Trj> Close error.[1]')

    if (.not. allocated(g_pt_info(handle%unit_no)%files)) &
      call error_msg('Close_Parallel_Trj> Close error.[2]')

    
    do i = 1, size(g_pt_info(handle%unit_no)%bufs)
      buf => g_pt_info(handle%unit_no)%bufs(i)
      if (allocated(buf%crd)) &
        deallocate(buf%crd, buf%idx)
    end do

    deallocate(g_pt_info(handle%unit_no)%bufs)
    deallocate(g_pt_info(handle%unit_no)%files)

    return

  end subroutine close_parallel_trj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_parallel_trj
  !> @brief        read parallel trajectory file
  !! @authors      NT
  !! @param[in]    handle     : handle for parallel trajectory
  !! @param[inout] trajectory : trajectory information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine read_parallel_trj(handle, trajectory)
  
    ! formal argments
    type(s_trj_file),        intent(in)    :: handle
    type(s_trajectory),      intent(inout) :: trajectory

    ! local variables
    integer                  :: i, ifile, ibuf, nt, it


    if (handle%unit_no <= 0 .or. handle%unit_no > MaxPTInfo) &
      call error_msg('Read_Parallel_Trj> read error.')


    !$omp parallel private(it, ifile, ibuf)
    !

#ifdef OMP
    nt = omp_get_num_threads()
#else
    nt = 1
#endif

    do i = 1, size(g_pt_info(handle%unit_no)%files), nt

#ifdef OMP
      it = omp_get_thread_num()
#else
      it = 0
#endif

      ifile = i + it
      if (ifile > size(g_pt_info(handle%unit_no)%files)) &
        cycle

      ibuf  = it + 1
      write(MsgOut,'(A,I5,A,I3)') '   File:', ifile, '  Buf:', ibuf

      call read_parallel_trj_(g_pt_info(handle%unit_no)%files(ifile), &
                              g_pt_info(handle%unit_no)%bufs (ibuf),  &
                              trajectory%coord,                       &
                              trajectory%pbc_box)

    end do

    !
    !$omp end parallel

    return

  end subroutine read_parallel_trj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_parallel_trj_
  !> @brief        read parallel trajectory file
  !! @authors      NT
  !! @param[inout] file  : file handle
  !! @param[inout] buf   : buffer data
  !! @param[inout] coord : coordinates
  !! @param[inout] box   : box of boundary condition
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine read_parallel_trj_(file, buf, coord, box)
  
    ! formal argments
    type(s_pt_file),         intent(inout) :: file
    type(s_pt_buf),          intent(inout) :: buf
    real(wp),                intent(inout) :: coord(:,:)
    real(wp),                intent(inout) :: box(:,:)

    ! local variables
    integer                  :: i, un, natom, natom_all, rank, nsnap
    character(4)             :: mark

    ! function
    integer                  :: ftell


    !$omp critical
    !
    call open_binary_file(un, file%name, IOFileInput,file%endian)
    !
    !$omp end critical

    call fseek(un, file%pos, 0)


    ! read header
    !
    if (file%pos == 0) then

      read(un) mark
      read(un) rank
      read(un) natom_all
      read(un) nsnap

    end if


    ! read body
    !
    read(un) buf%box

    read(un) natom
    
    if (buf%natom == 0) then

      allocate(buf%idx(1:natom), buf%crd(1:3,1:natom))
      buf%natom = natom

    else if (buf%natom < natom) then

      deallocate(buf%idx, buf%crd)
      allocate(buf%idx(1:natom), buf%crd(1:3,1:natom))
      buf%natom = natom

    end if

    read(un) buf%idx(1:natom)
    read(un) buf%crd(1:3,1:natom)

    file%pos = ftell(un)

    !$omp critical
    !
    call close_file(un)

    do i = 1, natom
      coord(1:3,buf%idx(i)) = buf%crd(1:3,i)
    end do

    box      = 0.0_wp
    box(1,1) = buf%box(1)
    box(2,2) = buf%box(3)
    box(3,3) = buf%box(6)

    !
    !$omp end critical

    return

  end subroutine read_parallel_trj_

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      get_ranked_filename
  !> @brief        get ranked filename
  !! @authors      NT
  !! @param[in]    filename : base file name
  !! @param[in]    rank_no  : rank number
  !! @param[in]    nplace   : number of places
  !! @return       ranked filename
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  function get_ranked_filename(filename, rank_no, nplace)
  
    ! return
    character(MaxFilename)   :: get_ranked_filename

    ! formal argments
    character(*),            intent(in)    :: filename
    integer,                 intent(in)    :: rank_no
    integer,                 intent(in)    :: nplace

    ! local variables
    integer                  :: ci1, ci2
    character(50)            :: fmt_str


    ! check filename
    !
    get_ranked_filename = filename

    ci1 = scan(filename, '(', .true.)
    ci2 = scan(filename, ')', .true.)

    if (ci1 == 0 .or. ci2 ==0 .or. ci1 > ci2) then
      call error_msg( &
        'Get_Ranked_Filename> Filename is not correctly ranked')
    end if

    write(fmt_str,'(A,I0,A,I0,A)') '(A,I',nplace,'.',nplace,',A)'

    write(get_ranked_filename,fmt=fmt_str) &
         trim(filename(1:ci1-1)), rank_no, &
         trim(filename(ci2+1:))

    return

  end function get_ranked_filename

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      get_file_endian
  !> @brief        get file endian
  !! @authors      NT
  !! @param[in]    filename : filename
  !! @return       endian
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_file_endian(filename)

    ! function
    integer                  :: get_file_endian

    ! formal arguments
    character(*),            intent(in) :: filename

    ! local variables
    integer                  :: file, ival


    ! check endianess
    !
    get_file_endian = IOFileNativeEndian

    ! GCC gfortan do not support raw binary format extension (form = "binary")
#ifndef __GFORTRAN__

    file = get_unit_no()

    open(file, &
         file    = filename,        &
         status  = 'old',           &
         form    = 'binary',        &
         access  = 'direct',        &
         convert = 'little_endian', &
         recl    = 4)
    read(file,rec=1) ival

    if (ival /= 4) then

      close(file)

      open(file, &
           file    = filename,      &
           status  = 'old',         &
           form    = 'binary',      &
           access  = 'direct',      &
           convert = 'big_endian',  &
           recl    = 4)
      read(file,rec=1) ival

      if (ival /= 4) &
        call error_msg('Get_File_Endian> unknown file format.')

      get_file_endian = IOFileBigEndian

    else

      get_file_endian = IOFileLittleEndian

    end if

    close(file)
    call free_unit_no(file)

#endif

    return

  end function get_file_endian

end module parallel_trj_mod
