!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   hardwareinfo_mod
!> @brief   utilities for get hardware information
!! @authors Chigusa Kobayashi (CK), Motoshi Kamiya (MK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module hardwareinfo_mod

  use fileio_str_mod
  use fileio_mod
  use string_mod
  use constants_mod
  use messages_mod

#ifdef USE_GPU
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif /* HAVE_MPI_GENESIS */
  use mpi_parallel_mod
  use, intrinsic :: iso_c_binding
#endif /* USE_GPU */

  implicit none
  private

  character(13), parameter  :: cpuinfo = "/proc/cpuinfo"

#ifdef USE_GPU
  interface get_device_count
    subroutine gpu_get_device_count(device_count) BIND(C)
      use, intrinsic :: iso_c_binding, only: c_int
      integer(c_int), intent(out) :: device_count
    end subroutine
  end interface get_device_count

  interface get_major_version
    subroutine gpu_get_major_version(device_id, major_version) BIND(C)
      use, intrinsic :: iso_c_binding, only: c_int
      integer(c_int), value       :: device_id
      integer(c_int), intent(out) :: major_version
    end subroutine
  end interface get_major_version

  interface get_minor_version
    subroutine gpu_get_minor_version(device_id, minor_version) BIND(C)
      use, intrinsic :: iso_c_binding, only: c_int
      integer(c_int), value       :: device_id
      integer(c_int), intent(out) :: minor_version
    end subroutine
  end interface get_minor_version

  interface set_device
    subroutine gpu_set_device(device_id) BIND(C)
      use, intrinsic :: iso_c_binding, only: c_int
      integer(c_int), value :: device_id
    end subroutine
  end interface set_device

  interface get_ecc_support
    subroutine gpu_get_ecc_support(device_id, ecc_support) BIND(C)
      use, intrinsic :: iso_c_binding, only: c_int
      integer(c_int), value       :: device_id
      integer(c_int), intent(out) :: ecc_support
    end subroutine
  end interface get_ecc_support

  interface get_device_name
    subroutine gpu_get_device_name(device_id, device_name) BIND(C)
      use, intrinsic :: iso_c_binding, only: c_int, c_char
      integer(c_int),    value       :: device_id
      character(c_char), intent(out) :: device_name(256)
    end subroutine
  end interface get_device_name

#endif /* USE_GPU */

  public  :: hw_information
  private :: runtime_information
  private :: build_information
  private :: genesis_information
  public  :: get_cpu_information
  public  :: get_cpu_flags
#ifdef USE_GPU
  public  :: assign_gpu
  private :: get_local_rank
#endif /* USE_GPU */
  private :: get_env_information

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hw_information
  !> @brief        write architecture and compile information
  !! @param[in]    get_arch_infomation
  !! @authors      CK
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hw_information

    call genesis_information

    call build_information

    call runtime_information

    return

  end subroutine hw_information

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    runtime_information
  !> @brief        write runtime information
  !! @authors      CK
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine runtime_information
    use, intrinsic :: iso_c_binding, only: c_int, c_char

    ! local variables
    character(MaxLine)       :: ldlibrary,  rt_user, rt_host
    character(MaxLine)       :: mpiversion
    character(80)            :: cpuname
    integer                  :: tval(8)
    integer                  :: leng, ierr
    integer(8)               :: year, month, day, hour, minute, sec
#ifdef USE_GPU
    integer                  :: itmp, my_device_id
    character(256)           :: gpu_modelname
    character(256)           :: gpu_major, gpu_minor
    logical                  :: gpu_ecc
    integer(c_int)           :: major, minor, ecc_support
    character(c_char)        :: device_name(256)
#endif /* USE_GPU */


#ifndef KCOMP
    write(MsgOut,'(A)') 'Runtime_Information> Machine and Library Information'

    call date_and_time(values=tval)
    year   = int(tval(1),kind=8)
    month  = int(tval(2),kind=8)
    day    = int(tval(3),kind=8)
    hour   = int(tval(5),kind=8)
    minute = int(tval(6),kind=8)
    sec    = int(tval(7),kind=8)
    write(MsgOut,'(A,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)') &
    '  date       = ',&
      year,"/",month,"/",day," ",hour,":",minute,":", sec

    call get_cpu_information(cpuname)
    write(MsgOut,*) ' cpu model    = ',trim(cpuname)
#endif
#ifdef USE_GPU
    major = -1
    minor = -1
    ecc_support = -1
    device_name = ""
    call assign_gpu(my_device_id)
    call get_major_version(my_device_id, major)
    call get_minor_version(my_device_id, minor)
    call get_ecc_support(my_device_id, ecc_support)
    call get_device_name(my_device_id, device_name)
    gpu_modelname(1:256) = ' '
    gpu_ecc = (ecc_support > 0)
    do itmp = 1, 256
      if (device_name(itmp) == CHAR(0)) then
        exit
      end if
      gpu_modelname(itmp:itmp) = device_name(itmp)
    end do
    write(gpu_major,'(i8)') major
    write(gpu_minor,'(i8)') minor
    write(MsgOut,*) ' gpu model    = ',trim(gpu_modelname), &
                    " (CC ", trim(adjustl(gpu_major)), ".", &
                    trim(adjustl(gpu_minor)), ")"
    write(MsgOut,*) ' gpu ECC      = ',gpu_ecc
#endif /* USE_GPU */
    call get_env_information(ldlibrary,rt_host,rt_user)
    write(MsgOut,*) ' exec. host   = ',trim(rt_user),'@',trim(rt_host)
    write(MsgOut,*) ' LD library   = ',trim(ldlibrary)

#ifdef HAVE_MPI_GENESIS
    call MPI_GET_LIBRARY_VERSION(mpiversion, leng, ierr)
    if (trim(mpiversion) .ne. "") then
      write(MsgOut,*) 'MPI Runtime = ',trim(mpiversion)
    endif
#endif
    write(MsgOut,'(A)')

    return

  end subroutine runtime_information

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    genesis_information
  !> @brief        write GENESIS information
  !! @authors      CK
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine genesis_information

    write(MsgOut,'(A)') 'GENESIS_Information> GENESIS Information'
    write(MsgOut,*) ' version      = ',VERSION
    write(MsgOut,*) ' commit ID    = ',COMPILE_GENESIS_VERSION
    write(MsgOut,*) ' precision    = ',precision_char
#ifdef USE_GPU
    write(MsgOut,*) ' nonbonding   = GPU'
#else
    write(MsgOut,*) ' nonbonding   = CPU'
#endif
    write(MsgOut,'(A)')

    return

  end subroutine genesis_information

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    build_information
  !> @brief        write build information
  !! @authors      CK
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine build_information

    write(MsgOut,'(A)') 'Build_Information> Compiler Information'
    write(MsgOut,*) ' build host   = ',COMPILE_USER,'@',COMPILE_HOST
    write(MsgOut,*) ' fortran      = ',COMPILE_FC_VER
    write(MsgOut,*) ' option       = ',COMPILE_FCFLAGS
    write(MsgOut,*) ' C            = ',COMPILE_CC_VER
    write(MsgOut,*) ' option       = ',COMPILE_CFLAGS
    write(MsgOut,*) ' defined var. = ',COMPILE_DEFINED_VARIABLES
    write(MsgOut,*) ' link option  = ',COMPILE_LDFLAGS
#ifdef HAVE_MPI_GENESIS
    write(MsgOut,*) ' MPI Compiler = ',MPITYPE,' MPI'
#endif
#ifdef CUDAGPU
    write(MsgOut,*) ' CUDA         = ',COMPILE_NVCC_VER
#endif
    write(MsgOut,'(A)')

    return

  end subroutine build_information

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_env_information
  !> @brief        write Enviroment information
  !! @authors      CK
  !! @param[out]   ldlibrary    : LD_LIBRARY_PATH
  !! @param[out]   rt_host      : host
  !! @param[out]   rt_user      : user
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine get_env_information(ldlibrary, rt_host, rt_user)

#if defined(INTEL)
    use ifport, only:getenv

#elif defined(KCOMP) || defined(RICC)
    use service_routines, only:getenv

#endif

    ! formal arguments
    character(*),            intent(out)    :: ldlibrary
    character(*),            intent(out)    :: rt_host
    character(*),            intent(out)    :: rt_user

    call getenv("LD_LIBRARY_PATH",ldlibrary)
    call getenv("USER",rt_user)
    call getenv("HOSTNAME",rt_host)

    return

  end subroutine get_env_information

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_cpu_information
  !> @brief        write CPU information
  !! @param[out]   cpuname    :  cpu model name
  !! @authors      CK
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine get_cpu_information(cpuname)

    ! formal arguments
    character(*),            intent(out)    :: cpuname

    ! local variables
    integer       :: i, inchar, file
    character(80) :: line
    logical       :: exists


    exists=.false.

#ifdef KCOMP
    cpuname="Fugaku"
    return
#endif
    inquire(file=cpuinfo, EXIST=exists)
    if (.not. exists) then
      cpuname="N/A"
      return
    endif

    call open_file(file, cpuinfo, IOFileInput)

    do while(.true.)
      read(file, '(A80)',end=100,err=100) line
      inchar = index(line,":")
      if (inchar > 0) then
        if (line(1:10) == 'model name') then
          read(line(inchar+2:80),'(A)') cpuname
          exit
        else if (inchar == 6 .and. line(1:3) == 'cpu') then
          read(line(inchar+2:80),'(A)') cpuname
          exit
        end if
      end if
    end do
100 continue

    call close_file(file)

    return

  end subroutine get_cpu_information

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_cpu_flags
  !> @brief        get CPU flags
  !! @param[out]   cpuflags    :  cpu flags
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine get_cpu_flags(cpuflags)

    ! formal arguments
    character(*),            allocatable   :: cpuflags(:)

    ! local variables
    integer                  :: i, inchar, file, cnt
    character(1000)          :: line, flags
    logical                  :: exists


    exists=.false.

    inquire(file=cpuinfo, EXIST=exists)
    if (.not. exists) then
      allocate(cpuflags(1))
      cpuflags(1) = "N/A"
      return
    end if

    call open_file(file, cpuinfo, IOFileInput)

    flags = ''
    do while(.true.)
      read(file, '(a)', end=100) line
      inchar = index(line,":")
      if (inchar > 0) then
        if (line(1:5) == 'flags') then
          read(line(inchar+2:),'(A)') flags
          exit
        end if
      end if
    end do

100 call close_file(file)

    if (flags == '') then
      allocate(cpuflags(1))
      cpuflags(1) = "N/A"
      return
    end if

    cnt = split_num(flags)
    allocate(cpuflags(cnt))
    call split(cnt, cnt, flags, cpuflags)

    return

  end subroutine get_cpu_flags

#ifdef USE_GPU
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_gpu
  !> @brief        assign GPU for each process
  !! @authors      MK
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_gpu(my_device_id)
    use, intrinsic :: iso_c_binding, only: c_int

    ! formal arguments
    integer, intent(out) :: my_device_id

    ! local variables
    integer :: ret
    integer :: my_count, major, minor
    integer :: counter, i, not_avail
    logical, allocatable :: avail_gpus(:)


    my_device_id = 0
    my_count = -1
    major = -1
    minor = -1
#ifndef HAVE_MPI_GENESIS
    ! do nothing if MPI is not available
    return
#else

    call get_device_count(my_count)

    not_avail = 0
    if (my_count > 0) then
      allocate(avail_gpus(my_count))
      avail_gpus(:) = .true.
      do i = 1, my_count
        ! check compute capability (is >= 3.5 or not)
        call get_major_version(i-1, major)
        call get_minor_version(i-1, minor)
        if (major < 3 .or. &
             (major == 3 .and. minor < 5)) then
          not_avail = not_avail + 1
          avail_gpus(i) = .false.
        end if
      end do
      my_count = my_count - not_avail
      if (not_avail > 0) then
        if (main_rank) then
          write(MsgOut,'(a,i5)') "  ignored GPUs =", not_avail
        end if
      end if
    end if

    if (main_rank) then
      write(MsgOut,'(a,i5)') "  # of GPUs    =", my_count
    end if
    ! if only 1 GPU is found, do nothing
    if (my_count == 1) return

    ! put error if no gpu found
    if (my_count == 0) then
      call error_msg('get_gpu_info> CUDA enabled GPU is not available')
    endif

    call get_local_rank()
    my_device_id = mod(my_node_local_rank,my_count)

    counter = -1
    do i = 1, my_count + not_avail
      if (avail_gpus(i)) counter = counter + 1
      if (counter == my_device_id) then
        my_device_id = i - 1
        exit
      end if
      if (i == (my_count + not_avail)) then
        call error_msg('get_gpu_info> Error: unexpected error')
      end if
    end do

    ! assign device
    call set_device(my_device_id)
#endif /* HAVE_MPI_GENESIS */

    return

  end subroutine assign_gpu

  !======1=========2=========3=========4=========5=========6=========7=========8  !
  !  Subroutine    get_local_rank
  !> @brief        get node local rank of this process
  !! @authors      MK
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine get_local_rank()

#if defined(INTEL)
    use ifport, only:getenv

#elif defined(KCOMP) || defined(RICC)
    use service_routines, only:getenv

#endif

    ! local variables
    character(MaxLine) :: myname
    integer            :: len_myname, i, ierr

    character(MaxLine), allocatable :: nodenames(:)
    integer,            allocatable :: ranks(:)


#ifdef HAVE_MPI_GENESIS
    ! check OpenMPI
    call getenv("OMPI_COMM_WORLD_LOCAL_RANK", myname)
    if (len_trim(myname) /= 0) then
      if (main_rank) then
        write(MsgOut,'(a)') &
          "get_local_rank> OpenMPI environment variable found."
        write(MsgOut,*)
      end if
      read(myname,*) my_node_local_rank
      return
    end if

    ! general case
    call MPI_Get_processor_name(myname, len_myname, ierr)

    allocate(nodenames(nproc_world),ranks(nproc_world))

    nodenames(1:nproc_world) = ""
    ranks(1:nproc_world)     = 0

    nodenames(my_world_rank+1) = myname
    ranks(my_world_rank+1)     = my_world_rank

    do i = 0, nproc_world - 1
      call MPI_Bcast(nodenames(i+1), MaxLine, MPI_CHARACTER, &
                     i, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(ranks(i+1), 1, MPI_INTEGER, i, MPI_COMM_WORLD, ierr)
    end do

    my_node_local_rank = 0
    do i = 0, nproc_world - 1
      if (trim(nodenames(my_world_rank+1)) == trim(nodenames(i+1))) then
        if (my_world_rank > ranks(i+1)) then
          my_node_local_rank = my_node_local_rank + 1
        end if
      end if
    end do

    deallocate(nodenames,ranks)
#endif /* HAVE_MPI_GENESIS */

    return

  end subroutine get_local_rank

#endif /* USE_GPU */

end module hardwareinfo_mod
