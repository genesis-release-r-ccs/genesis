!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pm_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module pm_analyze_mod

  use pm_option_str_mod
  use output_str_mod
  use input_str_mod
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

  ! constants
  real(wp),        parameter   :: KB   = 0.00198719168260038_wp
  !real(wp),        parameter   :: PI   = 3.14159265358979323_wp

  ! subroutines
  public  :: analyze
  private :: calc_pmf1d
  private :: calc_pmf2d
  private :: check_file_lines
  private :: check_file_column
  private :: get_replicate_name1
  private :: logsumexp2
  private :: periodic

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        run analyzing trajectories
  !! @authors      NT
  !! @param[in]    input    : input information
  !! @param[in]    output   : output information
  !! @param[in]    option   : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(input, output, option)

    ! formal arguments
    type(s_input),           intent(in)    :: input
    type(s_output),          intent(in)    :: output
    type(s_option),          intent(in)    :: option

    ! local variables
    real(wp)                 :: w, d1, d2, KBT, distance
    integer                  :: ireplica, i, j, file_w, file_d, file_dist, tmp
    integer                  :: nline, nline_w, nline_d, nline_dist, ncol_d
    real(wp)                 :: tim

    real(wp),    allocatable :: weight(:), data(:, :), pmf1d(:,:), pmf2d(:, :)


    if (option%check_only) &
      return

    KBT = KB * option%temperature

    ! setup data
    !
    call check_file_lines(get_replicate_name1(input%weightfile, 1), nline_w)
    call check_file_lines(get_replicate_name1(input%cvfile, 1), nline_d)
    call check_file_column(get_replicate_name1(input%cvfile, 1), ncol_d)

    if (nline_w /= 0 .and. nline_w /= nline_d) &
      call error_msg( &
      'Analyze> # of weight file lines is different from cv file lines')

    if (ncol_d < (option%dimension+1)) &
      call error_msg( &
      'Analyze> # of column of cv file must be >= dimension+1')

    if (input%distfile /= '') then
      call check_file_lines(get_replicate_name1(input%distfile, 1), nline_dist)

      if (nline_w /= 0 .and. nline_w /= nline_dist) &
        call error_msg( &
        'Analyze> # of weight file lines is different from distance file lines')
    end if

    allocate(weight  (nline_d * option%nreplicas), &
             data(option%dimension, nline_d * option%nreplicas))

    nline = 0
    do ireplica = 1, option%nreplicas

      file_w = 0
      file_d = 0
      file_dist = 0

      ! open weight file
      if (input%weightfile /= '') &
        call open_file(file_w, get_replicate_name1(input%weightfile, ireplica), &
                       IOFileInput)

      ! open cv file
      call open_file(file_d, get_replicate_name1(input%cvfile, ireplica), &
                     IOFileInput)

      ! open distance file
      if (input%distfile /= '') &
        call open_file(file_dist, get_replicate_name1(input%distfile, ireplica), &
                       IOFileInput)

      do i = 1, nline_d

        ! read weight file
        if (file_w /= 0) &
          read(file_w,*) tim, w

        ! read cv file
        if (option%dimension == 1) then
          read(file_d,*) tim, d1
          if (file_dist /= 0) then
            read(file_dist,*) tim, distance
            if (option%cutoff == 0.0_wp .or. distance < option%cutoff) then
              nline = nline + 1
              weight(nline) = w
              data(1, nline) = d1
            end if
          else
            nline = nline + 1
            weight(nline) = w
            data(1, nline) = d1
          end if
        else if (option%dimension == 2) then
          read(file_d,*) tim, d1, d2
          if (file_dist /= 0) then
            read(file_dist,*) tim, distance
            if (option%cutoff == 0.0_wp .or. distance < option%cutoff) then
              nline = nline + 1
              weight(nline) = w
              data(1, nline) = d1
              data(2, nline) = d2
            end if
          else
            nline = nline + 1
            weight(nline) = w
            data(1, nline) = d1
            data(2, nline) = d2
          end if
        end if

      end do

      ! close cv file
      call close_file(file_d)

      ! close weight file
      if (file_w /= 0) &
        call close_file(file_w)

      ! close distance file
      if (file_dist /= 0) &
        call close_file(file_dist)

    end do

    if (input%weightfile == '') &
      weight(1:nline) = 1.0_wp / real(nline,wp)


    ! calc pmf
    !
    if (option%dimension == 1) then
      call calc_pmf1d(option, data(:, 1:nline), weight(1:nline), pmf1d)
    else if (option%dimension == 2) then
      call calc_pmf2d(option, data(:, 1:nline), weight(1:nline), pmf2d)
    end if

    ! output pmf
    !
    call open_file(file_d, output%pmffile, IOFileOutputNew)

    if (option%dimension == 1) then
      ! Periodic CV cannot print out standard type PMF
      if (option%is_periodic(1)) then
        do i = 1, size(pmf1d(1,:))
          write(file_d,*) option%center(1, i), pmf1d(2, i)
        end do   
      else
        do i = 1, size(pmf1d(1,:))
          write(file_d,*) option%center(1, i), pmf1d(1, i), pmf1d(2, i)
        end do
      end if
    else if (option%dimension == 2) then

      if(option%output_type == OutputTypeMATLAB) then
        do j = 1, size(pmf2d(1, :))
          do i = 1, size(pmf2d(:, 1))
            write(file_d,'(es25.16e3,$)') pmf2d(i, j)
          end do
          write(file_d, *)
        end do
      else if(option%output_type == OutputTypeGNUPLOT) then
        do i = 1, size(pmf2d(:, 1))
          do j = 1, size(pmf2d(1, :))
            write(file_d, *) option%center(1, i), option%center(2, j), pmf2d(i, j)
          end do
          write(file_d, *)
        end do
      end if

    end if

    call close_file(file_d)

    if (option%dimension == 1) then
      deallocate(pmf1d)
    else if (option%dimension == 2) then
      deallocate(pmf2d)
    end if

    deallocate(weight, data)


    ! Output summary
    !
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze> Detailed information in the output files'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '  [pmffile] ' // trim(output%pmffile)
    if (option%dimension == 1) then
      ! Periodic CV cannot print out standard type PMF now 
      if (option%is_periodic(1)) then
        write(MsgOut,'(A)') '    Column 1: coordinates of grid centers'
        write(MsgOut,'(A)') '    Column 2: Free energy profile at the corresponding bin'
        write(MsgOut,'(A)') ''
      else
        write(MsgOut,'(A)') '    Column 1: coordinates of grid centers'
        write(MsgOut,'(A)') '    Column 2: Free energy profile at the corresponding bin &
                                 by standard method'
        write(MsgOut,'(A)') '    Column 3: Free energy profile at the corresponding bin &
                                 by Gaussian Distribution'
        write(MsgOut,'(A)') ''
      end if
    else if (option%dimension == 2) then
      write(MsgOut,'(A)') '    Row X and Column Y: free energy profile at bin center (X,Y)'
      write(MsgOut,'(A)') ''
    end if

    return

  end subroutine analyze

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calc_pmf1d(option, data, weight, pmf)

    ! parameters
    real(wp),      parameter :: RMIN = 1.1755e-38

    ! formal arguments
    type(s_option),          intent(in)    :: option
    real(wp),                intent(in)    :: data(:, :)
    real(wp),                intent(in)    :: weight(:)
    real(wp),                allocatable   :: pmf(:,:)
    real(wp)                               :: KBT

    ! local variables
    integer                  :: nbin, ndata
    integer                  :: ibin, idata
    integer                  :: ipmf 

    real(wp)                 :: ZERO = 0.0_wp
    real(wp)                 :: NaN
    real(wp)                 :: ddata = 0.0_wp

    real(wp),    allocatable :: dx1(:,:), x(:,:), s(:)


    NaN = 0.0_wp/ZERO
    KBT = KB * option%temperature

    nbin  = size(option%center(1, :))
    ndata = size(data(1, :))

    allocate(pmf(2, nbin))
    allocate(dx1(nbin, ndata), x(ndata, nbin), s(nbin))

    ! impf = 1 > Standard PMF, 2 > Gaussian type PMF 
    do ipmf = 1, 2
      dx1(1:nbin, 1:ndata) = 0.0_wp

      ! Periodic CV cannot print out standard type PMF now
      if((ipmf == 1) .and. (option%is_periodic(1))) exit
      do ibin = 1, nbin
        if (option%is_periodic(1)) then
          do idata = 1, ndata
            dx1(ibin, idata) = (periodic(data(1, idata), option%center(1, ibin), &
              option%box_size(1)) / option%band_width(1))**2.0_wp
          end do
        else
          if(ipmf == 2) then
            ! For Gauusian type PMF
            dx1(ibin, 1:ndata) = ((data(1, 1:ndata) - option%center(1, ibin)) / &
              option%band_width(1)) ** 2.0_wp
          else
            ! For standard type PMF
            do idata = 1, ndata
              ddata = data(1, idata) - option%center(1, ibin)
              if((-option%delta_grid(1) < ddata) .and. (option%delta_grid(1) > ddata)) &
                dx1(ibin, idata) = dx1(ibin, idata) + 1.0_wp
            end do
          end if
        end if
        !x  (1:nbin,idata) = -0.5_wp * dx1(idata, 1:nbin) * log(weight(idata))

        if(ipmf == 2) dx1(ibin, 1:ndata) = exp(-0.5_wp * dx1(ibin, 1:ndata)) / &
                                           (sqrt(2.0_wp*PI)*option%band_width(1))
        do idata = 1, ndata
          dx1(ibin, idata) = dx1(ibin, idata)*weight(idata)
        end do
      end do

      if(ipmf == 1) then
        ! Check Histogram
        write(MsgOut,'(A,f5.3)') '    Summation P(r) in each = ', sum(dx1)
      end if   

      do ibin = 1, nbin
        pmf(ipmf, ibin) = 0.0_wp
        do idata = 1, ndata
          pmf(ipmf, ibin) = pmf(ipmf, ibin) + dx1(ibin, idata)
        end do
      end do

      do ibin = 1, nbin
        pmf(ipmf, ibin) = -KBT*log(pmf(ipmf, ibin))
      end do

      ! call logsumexp2(x, s)
  
      ! do ibin = 1, nbin
      !   if (s(ibin) < RMIN) &
      !     s(ibin) = NaN
      ! end do
  
      ! pmf(1:nbin) = -log(s(1:nbin))

      ! Lowest PMF is set to be 0.0
      pmf(ipmf, 1:nbin) = pmf(ipmf, 1:nbin) - minval(pmf(ipmf, 1:nbin))

    end do  

    deallocate(s, x, dx1)

    return

  end subroutine calc_pmf1d

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calc_pmf2d(option, data, weight, pmf)

    ! parameters
    real(wp),      parameter :: RMIN = 1.1755e-38

    ! formal arguments
    type(s_option),          intent(in)    :: option
    real(wp),                intent(in)    :: data(:, :)
    real(wp),                intent(in)    :: weight(:)
    real(wp),                allocatable   :: pmf(:, :)
    real(wp)                               :: KBT

    ! local variables
    integer                  :: ibin, nbin1, nbin2
    integer                  :: idata, ndata
    integer                  :: i, j

    real(wp)                 :: ZERO = 0.0_wp
    real(wp)                 :: NaN

    real(wp),    allocatable :: dx1(:,:), dx2(:,:), x(:,:), s(:)


    NaN = 0.0_wp/ZERO
    KBT = KB * option%temperature

    nbin1  = option%num_grids(1)-1
    nbin2  = option%num_grids(2)-1
    ndata  = size(data(1, :))

    ! write(*,*)'nbin1 = ', nbin1
    ! write(*,*)'nbin2 = ', nbin2
    ! write(*,*)'ndata = ', ndata
    ! write(*,*)'option%center(1, :) = ', option%center(1, :)
    ! write(*,*)'option%center(2, :) = ', option%center(2, :)

    allocate(pmf(nbin1, nbin2))
    allocate(dx1(nbin1, ndata), dx2(nbin2, ndata))

    do i = 1, nbin1
      if (option%is_periodic(1)) then
        do idata = 1, ndata
          dx1(i, idata) = (periodic(data(1, idata), option%center(1, i), option%box_size(1)) &
            /option%band_width(1))**2.0_wp
        end do
      else
        dx1(i, 1:ndata) = ((data(1, 1:ndata) - option%center(1, i)) / &
          option%band_width(1)) ** 2.0_wp
      end if
      dx1(i, 1:ndata) = exp(-0.5_wp * dx1(i, 1:ndata)) / &
                           (sqrt(2.0_wp*PI)*option%band_width(1))
    end do

    do j = 1, nbin2
      if (option%is_periodic(2)) then
        do idata = 1, ndata
          dx2(j, idata) = (periodic(data(2, idata), option%center(2, j), option%box_size(2)) &
            /option%band_width(2))**2.0_wp
        end do
      else
        dx2(j, 1:ndata) = ((data(2, 1:ndata) - option%center(2, j)) / &
          option%band_width(2)) ** 2.0_wp
      end if
      dx2(j, 1:ndata) = exp(-0.5_wp * dx2(j, 1:ndata)) / &
                           (sqrt(2.0_wp*PI)*option%band_width(2))
    end do

    pmf(:, :) = 0.0_wp
    do i = 1, nbin1
      do j = 1, nbin2
        do idata = 1, ndata
          pmf(i, j) = pmf(i, j) + dx1(i, idata)*dx2(j, idata)*weight(idata)
        end do
      end do
    end do

    do i = 1, nbin1
      do j = 1, nbin2
        pmf(i, j) = -KBT*log(pmf(i, j))
      end do
    end do

    pmf(1:nbin1, 1:nbin2) = pmf(1:nbin1, 1:nbin2) - minval(pmf(1:nbin1, 1:nbin2))

    deallocate(dx1, dx2)

    return

  end subroutine calc_pmf2d

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_file_lines(filename, nline)

    ! formal arguments
    character(*),            intent(in)    :: filename
    integer,                 intent(inout) :: nline

    ! local variables
    integer                  :: file
    character(100)           :: line


    if (filename == '') then
      nline = 0
      return
    end if

    call open_file(file, filename, IOFileInput)

    nline = 0
    do while(.true.)
      read(file,'(a)',end=10) line
      nline = nline + 1
    end do

10  call close_file(file)

    write(MsgOut,'(a)')    'Check_File_Lines>'
    write(MsgOut,'(a,a)')  '    file name   : ', trim(filename)
    write(MsgOut,'(a,i0)') '    # of line   : ', nline
    write(MsgOut,'(a)')    ''

    return

  end subroutine check_file_lines

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_file_column(filename, ncol)

    ! formal arguments
    character(*),            intent(in)    :: filename
    integer,                 intent(inout) :: ncol

    ! local variables
    integer                  :: file
    character(1000)          :: line


    call open_file(file, filename, IOFileInput)

    read(file,'(a)') line
    ncol = split_num(line)

    call close_file(file)

    write(MsgOut,'(a)')    'Check_File_Column>'
    write(MsgOut,'(a,a)')  '    file name   : ', trim(filename)
    write(MsgOut,'(a,i0)') '    # of column : ', ncol
    write(MsgOut,'(a)')    ''

    return

  end subroutine check_file_column

  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_replicate_name1(filename, no)

    ! return
    character(Maxfilename)   :: get_replicate_name1

    ! formal arguments
    character(*),            intent(in)    :: filename
    integer,                 intent(in)    :: no

    ! local variables
    integer                  :: bl, br


    if (filename == '') then
      get_replicate_name1 = ''
      return
    end if

    bl = index(filename, '{', back=.true.)
    br = index(filename, '}', back=.true.)

    if (bl == 0 .or. br == 0 .or. bl > br) &
      call error_msg('Get_Replicate_Name1> Syntax error.')

    write(get_replicate_name1, '(a,i0,a)') &
         filename(:bl-1),no,filename(br+1:len_trim(filename))

    return

  end function get_replicate_name1

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine logsumexp2(x, s)

    ! formal arguments
    real(wp),                intent(in)    :: x(:,:)
    real(wp),                intent(inout) :: s(:)

    ! local variables
    integer                  :: nbin, ibin, ndata

    real(wp),    allocatable :: max_x(:), exp_x(:,:)


    ndata = size(x(1,:))
    nbin  = size(s)

    allocate(max_x(nbin), exp_x(nbin,ndata))

    do ibin = 1, nbin
      max_x(ibin) = maxval(x(ibin,1:ndata))
    end do

    do ibin = 1, nbin
      exp_x(ibin,1:ndata) = exp(x(ibin,1:ndata) - max_x(ibin))
    end do

    do ibin = 1, nbin
      s(ibin) = log(sum(exp_x(ibin,1:ndata))) + max_x(ibin)
    end do

    deallocate(exp_x, max_x)

    return

  end subroutine logsumexp2

  !======1=========2=========3=========4=========5=========6=========7=========8

  function periodic(x, center, box_size)

    ! return value
    real(wp)                 :: periodic

    ! formal arguments
    real(wp)                 :: x
    real(wp)                 :: center
    real(wp)                 :: box_size


    periodic = x - center
    periodic = periodic - nint(periodic/box_size)*box_size

    return

  end function periodic

end module pm_analyze_mod
