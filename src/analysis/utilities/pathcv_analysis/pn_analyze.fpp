!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pn_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module pn_analyze_mod

  use pn_option_str_mod
  use trajectory_str_mod
  use output_str_mod
  use molecules_mod
  use molecules_str_mod
  use input_str_mod
  use fileio_trj_mod
  use fileio_mod
  use string_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type s_data
    real(wp), allocatable :: v(:,:)
  end type s_data

  ! subroutines
  public  :: analyze
  private :: read_data
  private :: calc_pathcv
  private :: write_pathcv
  private :: get_replicate_name1

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        run analyzing trajectories
  !! @authors      NT, CK
  !! @param[in]    input    : input information
  !! @param[in]    molecule   : molecule information
  !! @param[in]    trj_list   : trajectory file list information
  !! @param[in]    output   : output information
  !! @param[in]    option   : option information
  !! @param[inout] trajectory : trajectory information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(input, output, option, molecule, trj_list, trajectory)

    ! formal arguments
    type(s_input),                     intent(in)    :: input
    type(s_output),                    intent(in)    :: output
    type(s_option),                    intent(in)    :: option
    type(s_molecule),    optional,     intent(in)    :: molecule
    type(s_trj_list),    optional,     intent(in)    :: trj_list
    type(s_trajectory),  optional,     intent(inout) :: trajectory

    ! local variables
    type(s_data)             :: path, cv
    integer                  :: ireplica

    real(wp),     allocatable :: progress(:)
    real(wp),     allocatable :: distance(:)


    if (option%check_only) &
      return


    ! read path file
    !
    call read_data(input%pathfile, path)


!    do ireplica = 1, option%nreplicas
    do ireplica = 1, option%nfiles

      ! read cv files
      !
      if (option%trajectory) then
        call read_trajectory_path(ireplica, molecule, trj_list, trajectory, cv)
      else
        if (option%nfiles .eq. 1) then
          call read_data(input%cvfile, cv)
        else
          call read_data(get_replicate_name1(input%cvfile, ireplica), cv)
        endif
      endif


      ! calculate path cv
      !
      call calc_pathcv(path, cv, progress, distance)


      ! write pathcv files
      !
      if (option%nfiles .eq. 1) then
        call write_pathcv(output%pathcvfile, progress, distance)
      else
        call write_pathcv(get_replicate_name1(output%pathcvfile, ireplica), &
                        progress, distance)
      endif

    end do

    return

  end subroutine analyze

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_trajectory_path(file_no, molecule, trj_list, trajectory, data)
    ! formal arguments
    integer,                 intent(in)    :: file_no
    type(s_molecule),        intent(in)    :: molecule
    type(s_trj_list),        intent(in)    :: trj_list
    type(s_trajectory),      intent(inout) :: trajectory
    type(s_data),            intent(inout) :: data

    ! local variables
    type(s_trj_file)         :: trj_in
    integer                  :: nstru, istep, natom3, jmol, jj

    call open_trj(trj_in, &
                trj_list%filenames(file_no), &
                trj_list%trj_format,       &
                trj_list%trj_type, IOFileInput)
    
    natom3 = molecule%num_atoms*3

    if (allocated(data%v)) &
      deallocate(data%v)
    allocate(data%v(trj_list%md_steps(file_no),natom3))

    nstru=0
    do istep = 1, trj_list%md_steps(file_no)
      call read_trj(trj_in, trajectory)
      if (mod(istep, trj_list%ana_periods(file_no)) == 0) then
        nstru = nstru + 1
        jj = 0
        do jmol = 1, molecule%num_atoms
          jj = jj+1
          data%v(nstru, jj)=trajectory%coord(1,jmol)
          jj = jj+1
          data%v(nstru, jj)=trajectory%coord(2,jmol)
          jj = jj+1
          data%v(nstru, jj)=trajectory%coord(3,jmol)
        end do
      endif
    end do
    call close_trj(trj_in)

  end subroutine read_trajectory_path

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_data(filename, data)

    ! formal arguments
    character(*),            intent(in)    :: filename
    type(s_data),            intent(inout) :: data

    ! local variables
    real(wp)                               :: dummy
    integer                                :: file, ndim, nstep, i, j
    character(MaxLineLong_CV)              :: line


    call open_file(file, filename, IOFileInput)

    read(file,'(a)') line
    ndim = split_num(line) - 1
    if (ndim == 0) &
      call error_msg('Read_Data> bad format. '//trim(filename))

    nstep = 1
    do while(.true.)
      read(file,*,end=10,err=10) line
      nstep = nstep + 1
    end do

10  rewind(file)
    write(MsgOut,'(a)')    'Read_Data> '
    write(MsgOut,'(a,a)')  '  file name = ', trim(filename)
    write(MsgOut,'(a,i0)') '    # of dimension : ', ndim
    write(MsgOut,'(a,i0)') '    # of steps     : ', nstep
    write(MsgOut,'(a)')    ''

    if (allocated(data%v)) &
      deallocate(data%v)
    allocate(data%v(nstep,ndim))
    do i = 1, nstep
      read(file,*) dummy, (data%v(i,j), j=1,ndim)
    end do

    call close_file(file)

    return

  end subroutine read_data

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calc_pathcv(path, cv, progress, distance)

    ! formal arguments
    type(s_data),            intent(in)    :: path
    type(s_data),            intent(in)    :: cv
    real(wp),                allocatable   :: progress(:)
    real(wp),                allocatable   :: distance(:)

    ! local variables
    real(wp)                 :: d_pathpoints, lambda
    integer                  :: ndim, nstep_p, nstep_c, istep_c, istep_p, idim

    real(wp),    allocatable :: v(:), log_m(:)
    real(wp),    allocatable :: dev(:,:), d(:,:), tmp(:,:)
    real(wp),    allocatable :: log_numerator(:), log_denominator(:)


    write(MsgOut,'(a)')    'Calc_Pathcv> '

    if (allocated(progress)) &
      deallocate(progress)

    if (allocated(distance)) &
      deallocate(distance)

    ndim    = size(path%v(1,:))
    nstep_p = size(path%v(:,1))
    nstep_c = size(cv%v(:,1))

    if (ndim /= size(cv%v(1,:))) &
      call error_msg('Calc_Pathcv> dimension of ref and data do not match.')

    allocate(progress(nstep_c))
    allocate(distance(nstep_c))

    allocate(v(ndim))

    ! determine lambda parameter
    d_pathpoints = 0
    do istep_p = 1, nstep_p-1
      v(1:ndim) = path%v(istep_p+1,1:ndim) - path%v(istep_p,1:ndim)
      v(1:ndim) = v(1:ndim) ** 2.0_wp
      d_pathpoints = d_pathpoints + sqrt(sum(v))
    end do
    d_pathpoints = d_pathpoints / real(nstep_p - 1,wp)
    lambda = 2.3_wp / (d_pathpoints ** 2.0_wp)

    ! calculate pathway CVs
    allocate(dev(nstep_c, ndim))
    allocate(d  (nstep_c, nstep_p))
    allocate(tmp(nstep_c, nstep_p))
    allocate(log_numerator  (nstep_c))
    allocate(log_denominator(nstep_c))

    do istep_p = 1, nstep_p
      do istep_c = 1, nstep_c
        dev(istep_c,1:ndim) = path%v(istep_p,1:ndim) - cv%v(istep_c,1:ndim)
        d(istep_c, istep_p) = sum(dev(istep_c,1:ndim)**2.0_wp)
      end do
    end do

    allocate(log_m(nstep_p))
    do istep_p = 1, nstep_p
      log_m(istep_p) = log(real(istep_p, wp))
    end do

    do istep_c = 1, nstep_c
      tmp(istep_c,1:nstep_p) = log_m(1:nstep_p) - lambda * d(istep_c,1:nstep_p)
    end do
    call logsumexp_array(tmp, log_numerator)

    do istep_c = 1, nstep_c
      tmp(istep_c,1:nstep_p) = - lambda * d(istep_c,1:nstep_p)
    end do
    call logsumexp_array(tmp, log_denominator)

    progress(1:nstep_c) = &
         exp(log_numerator(1:nstep_c) - log_denominator(1:nstep_c))
    
    distance(1:nstep_c) = &
         - (1.0_wp / lambda) * log_denominator(1:nstep_c)

    deallocate(log_m, &
               log_denominator, &
               log_numerator,   &
               tmp, d, dev, v)

    write(MsgOut,'(a)')    '  ..done'
    write(MsgOut,'(a)')    ''

    return

  end subroutine calc_pathcv

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_pathcv(filename, progress, distance)

    ! formal arguments
    character(*),            intent(in)    :: filename
    real(wp),                intent(in)    :: progress(:)
    real(wp),                intent(in)    :: distance(:)

    ! local variables
    integer                  :: file, istep, nstep
    real(wp)                 :: ave_prog, std_prog
    real(wp)                 :: ave_dist, std_dist, rstep


    call open_file(file, filename, IOFileOutputNew)

    nstep = size(progress)

    ave_prog = 0.0_wp
    ave_dist = 0.0_wp
    std_prog = 0.0_wp
    std_dist = 0.0_wp
    rstep    = 1.0_wp/real(nstep,wp)
    do istep = 1, nstep
      write(file,'(i8,2f25.16)') istep, progress(istep), distance(istep)
      ave_prog = ave_prog + progress(istep)
      std_prog = std_prog + progress(istep)*progress(istep)
      ave_dist = ave_dist + distance(istep)
      std_dist = std_dist + distance(istep)*distance(istep)
    end do
    ave_prog = ave_prog * rstep
    std_prog = std_prog * rstep
    ave_dist = ave_dist * rstep
    std_dist = std_dist * rstep
    std_prog = sqrt(std_prog-ave_prog*ave_prog)
    std_dist = sqrt(std_dist-ave_dist*ave_dist)

    write(MsgOut, '(a,F12.5,a,F12.5)') 'Progress : Ave',ave_prog, \
                                       ' Std',std_prog 
    write(MsgOut, '(a,F12.5,a,F12.5)') 'Distance : Ave',ave_dist, \
                                       ' Std',std_dist 

    call close_file(file)

    return

  end subroutine write_pathcv

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine logsumexp_array(x, s)

    ! formal arguments
    real(wp),                intent(in)    :: x(:,:)
    real(wp),                intent(inout) :: s(:)

    ! local variables
    real(wp)                 :: max_x
    integer                  :: nstep_c, nstep_p, istep_c

    real(wp),    allocatable :: v(:)


    nstep_p = size(x(1,:))
    nstep_c = size(x(:,1))

    allocate(v(nstep_p))

    do istep_c = 1, nstep_c

      max_x = maxval(x(istep_c,1:nstep_p))
      v(1:nstep_p) = exp(x(istep_c,1:nstep_p) - max_x)

      s(istep_c) = log(sum(v)) + max_x

    end do

    deallocate(v)

    return

  end subroutine logsumexp_array

  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_replicate_name1(filename, no)

    ! return
    character(Maxfilename)   :: get_replicate_name1

    ! formal arguments
    character(*),            intent(in)    :: filename
    integer,                 intent(in)    :: no

    ! local variables
    integer                  :: bl, br


    bl = index(filename, '{', back=.true.)
    br = index(filename, '}', back=.true.)

    if (bl == 0 .or. br == 0 .or. bl > br) &
      call error_msg('Get_Replicate_Name1> Syntax error.')

    write(get_replicate_name1, '(a,i0,a)') &
         filename(:bl-1),no,filename(br+1:)

    return

  end function get_replicate_name1

end module pn_analyze_mod
