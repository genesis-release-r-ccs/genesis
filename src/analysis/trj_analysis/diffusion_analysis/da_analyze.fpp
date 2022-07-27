!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   da_analyze_mod
!> @brief   analyze mean square displacement
!! @authors Donatas Surblys (DS), Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module da_analyze_mod

  use da_option_str_mod
  use input_str_mod
  use output_str_mod
  use fileio_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public  :: analyze
  private :: fit_least_squares
  private :: get_column_count
  private :: get_line_count

  ! fitting result for a + bx
  type :: s_fitting_result
    real(wp), dimension(2) :: coeff   ! coefficients of the fitting (a, b)
    real(wp), dimension(2) :: stderr  ! standard error of each coefficient
    real(wp)               :: corr    ! correlation coefficient, showing 
                                      ! the quality of the fit
  end type s_fitting_result

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        analyze translational diffusion coefficient
  !! @authors      DS, TM
  !! @param[in]    molecule   : molecule structure
  !! @param[in]    output     : output information
  !! @param[in]    option     : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(input, output, option)

    ! formal arguments
    type(s_input),           intent(in)    :: input
    type(s_output),          intent(in)    :: output
    type(s_option),          intent(in)    :: option

    ! local variables
    integer                                :: msd_data, fit_data, ndata, ncols
    integer                                :: iset, idata, nlines
    integer                                :: idx_start_fit, idx_stop_fit, iline
    integer,  allocatable, dimension(:)    :: ndofs
    real(wp)                               :: diffusion_coefficient
    real(wp), allocatable, dimension(:, :) :: xydata
    character(*), parameter                :: format_float = "es25.16e3"
    type(s_fitting_result), dimension(:), allocatable :: fittings


    if (option%check_only) &
      return


    ! check data size
    !
    nlines  = get_line_count (trim(input%msdfile))
    ndata   = nlines
    ncols   = get_column_count(trim(input%msdfile), 1)

    allocate(ndofs(ncols-1))
    if (size(option%ndofs) == 1) then
      ndofs(:) = option%ndofs(1)
    else if (size(option%ndofs) == size(ndofs)) then
      ndofs(:) = option%ndofs
    else
      call error_msg('Analyze> Number of degrees of freedom in "ndofs"'//&
        ' does not match the number of MSD sets found in data file')
    end if

    allocate(xydata(ncols, ndata))
    allocate(fittings(ncols-1))
    

    ! read MDS data
    !
    call open_file(msd_data, trim(input%msdfile), IOFileInput)
    read(msd_data, *) xydata
    call close_file(msd_data)

    ! convert to ps and angstroms
    !
    xydata(1, :)  = xydata(1, :)  * option%time_step
    xydata(2:, :) = xydata(2:, :) * option%distance_unit ** 2
    do iset = 2, ncols
      xydata(iset, :) = xydata(iset, :) / (2.0_wp * ndofs(iset-1))
    end do

    ! determine data index from where to fit
    !
    if (option%start_percent >= 0_wp) then
      idx_start_fit = max(nint(ndata*option%start_percent/100), 1)
    else if (option%start_time >= 0_wp) then
      idx_start_fit = max(minloc(abs(xydata(1, :)-option%start_time), 1), 1)
    else if (option%start_step > 0) then
      idx_start_fit = option%start_step
    else
      idx_start_fit = 1
    end if

    ! determine data index up to where to fit
    !
    if (option%stop_percent >= 0_wp) then
      idx_stop_fit = min(nint(ndata*option%stop_percent/100), ndata)
    else if (option%stop_time >= 0_wp) then
      idx_stop_fit = min(minloc(abs(xydata(1, :)-option%stop_time), 1), ndata)
    else if (option%stop_step > 0) then
      idx_stop_fit = option%stop_step
    else
      idx_stop_fit = ndata
    end if


    ! Analyze
    !
    write(MsgOut, '()')
    write(MsgOut, '("Analyze> Starting fit at",es9.2e2," ps and using ",i0,' &
      //'" out of ",i0," available sample points with ", i0, " data sets")') &
      xydata(1, idx_start_fit), idx_stop_fit - idx_start_fit + 1, ndata, ncols-1

    write(MsgOut, '("Analyze> Fitting function (A^2/ps): f(x) = b * x + a ")')
    write(MsgOut, '()')

    do iset = 1, ncols-1

      ! results are in A^2/ps
      !
      fittings(iset) = fit_least_squares(xydata(1,idx_start_fit:idx_stop_fit), &
        xydata(iset+1, idx_start_fit:idx_stop_fit), .true.)
      write(MsgOut, &
        '("Analyze> (Set ",i0,") Fitting coefficient:        a =",'&
        //format_float//')') iset, fittings(iset)%coeff(1)
      write(MsgOut, &
        '("Analyze> (Set ",i0,") Fitting coefficient:        b =",'&
        //format_float//')') iset, fittings(iset)%coeff(2)
      write(MsgOut, &
        '("Analyze> (Set ",i0,") Standard error:         SE(a) =",'&
        //format_float//')') iset, fittings(iset)%stderr(1)
      write(MsgOut, &
        '("Analyze> (Set ",i0,") Standard error:         SE(b) =",'&
        //format_float//')') iset, fittings(iset)%stderr(2)
      write(MsgOut, &
        '("Analyze> (Set ",i0,") Correleation coefficient:   r =",'&
        //format_float//')') iset, fittings(iset)%corr

      ! convert diffusion coefficient to cm/s
      !
      diffusion_coefficient = fittings(iset)%coeff(2) * 1e-4_wp
      write(MsgOut, &
        '("Analyze> (Set ",i0,") Diffusion coefficient (cm^2/s):",'&
        //format_float//')') iset, diffusion_coefficient
      write(MsgOut, '()')

    end do


    ! Output results
    !
    if (output%outfile /= '') then

      write(MsgOut, '("Analyze> Writing fitted data to ", A)') trim(output%outfile)

      call open_file(fit_data, trim(output%outfile), IOFileOutputNew)

      do idata = 1, ndata
        write(fit_data, '('//format_float//')', advance="no") xydata(1, idata)
        do iset = 1, ncols-1
          write(fit_data, '(2'//format_float//')', advance="no") &
            xydata(iset+1, idata), &
            fittings(iset)%coeff(1)+fittings(iset)%coeff(2)*xydata(1, idata)
        end do
        write(fit_data, '()')
      end do

      call close_file(fit_data)

    end if


    ! Output Summary
    !
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze> Detailed information in the output files'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '  [outfile] ' // trim(output%outfile)
    write(MsgOut,'(A)') '    Column 1: time (ps)'
    do iset = 1, ncols-1
      write(MsgOut,'(A,i0,A,i0,A,i0)') &
                        '    Column ', iset*2, ': MSD/', ndofs(iset)*2, &
                        ' (A) from data set ', iset
      write(MsgOut,'(A,i0,A,i0)')      &
                        '    Column ', iset*2+1, ': fitting of column ', iset*2
    end do
    write(MsgOut,'(A)') ''


    return

  end subroutine analyze

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      get_line_count
  !> @brief        count number of lines in file
  !! @authors      DS
  !! @param[in]    file_name  : file name
  !! @return       integer containing number of lines in file
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_line_count(file_name) result(n)

    character(*), intent(in) :: file_name

    integer :: n
    integer :: file_unit
    integer :: iostat

    call open_file(file_unit, file_name, IOFileInput)

    n = -1
    iostat = 0
    do while(iostat == 0)
      read(file_unit, *, iostat=iostat)
      n = n + 1
    end do

    call close_file(file_unit)

    return

  end function get_line_count

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      get_column_count
  !> @brief        count number of columns in file at the specified row
  !! @authors      DS
  !! @param[in]    file_name  : file name
  !! @param[in]    row        : row number
  !! @return       integer containing number of columns at the specified row
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_column_count(file_name, row) result(n)

    character(*),      intent(in) :: file_name
    integer, optional, intent(in) :: row
    integer                       :: n

    character(*), parameter :: characters_whitespace = "  "

    logical      :: in_whitespace, is_whitespace
    integer      :: irow, nrow
    integer      :: file_unit
    integer      :: iostat
    character(1) :: c


    if (present(row)) then
      nrow = row
    else
      nrow = 1
    end if

    n = 0

    call open_file(file_unit, file_name, IOFileInput)

    do irow = 1, nrow-1
      read(file_unit, *, iostat=iostat)
      if (iostat /= 0) then
        call close_file(file_unit)
        return
      end if
    end do

    in_whitespace = .true.

    do
      read(file_unit, '(A1)', advance="no", iostat=iostat) c

      if (iostat /= 0) exit

      is_whitespace = (scan(c, characters_whitespace) /= 0)

      if (in_whitespace .and. .not. is_whitespace) then
        n = n + 1
      end if

      in_whitespace = is_whitespace

    end do

    call close_file(file_unit)

    return

  end function get_column_count

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      fit_least_squares
  !> @brief        fit line to data using least squares method
  !! @authors      DS
  !! @param[in]    xdata      : one dimensinal array of floats
  !! @param[in]    ydata      : one dimensinal array of floats
  !! @param[in]    normalize  : whether to normalize date before fitting
  !! @return       a data type containing fitting results
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function fit_least_squares(xdata, ydata, normalize) result(res)

    ! Formulas were taken from
    ! http://mathworld.wolfram.com/LeastSquaresFitting.html

    real(wp), dimension(:),           intent(in) :: xdata
    real(wp), dimension(size(xdata)), intent(in) :: ydata
    logical, optional                            :: normalize

    type(s_fitting_result) :: res

    logical  :: norm
    integer  :: i, n
    real(wp) :: x_mean, y_mean, x_min, x_max, y_min, y_max, x_scale, y_scale
    real(wp) :: ss_xx, ss_xy, ss_yy, a, b, r2, s, se_a, se_b

    real(wp), dimension(:), allocatable :: xwork, ywork


    norm = .false.
    if (present(normalize)) then
      norm = normalize
    end if

    n = size(xdata)

    allocate(xwork(n), ywork(n))

    if (norm) then

      x_min = minval(xdata)
      x_max = maxval(xdata)
      x_scale = x_max - x_min
      y_min = minval(ydata)
      y_max = maxval(ydata)
      y_scale = y_max - y_min

      xwork = (xdata - x_min) / x_scale
      ywork = (ydata - y_min) / y_scale

    else

      xwork = xdata
      ywork = ydata

    end if

    x_mean = sum(xwork) / n
    y_mean = sum(ywork) / n

    ss_xx = sum([ ((xwork(i)-x_mean)**2,i=1,n) ])
    ss_yy = sum([ ((ywork(i)-y_mean)**2,i=1,n) ])
    ss_xy = sum([ ((xwork(i)-x_mean)*(ywork(i)-y_mean),i=1,n) ])

    b = ss_xy / ss_xx
    a = y_mean - b  * x_mean

    r2 = ss_xy**2/(ss_xx*ss_yy)

    s = sqrt((ss_yy - b * ss_xy)/(n - 2))

    se_b = s / sqrt(ss_xx)

    if (norm) then

      b = b * y_scale / x_scale
      a = y_scale * a - b * x_min + y_min

      se_b = se_b * y_scale / x_scale
      se_a = s * y_scale &
        * sqrt(1_wp / n + (x_scale * x_mean + x_min)**2 / (x_scale)**2 / ss_xx)
    else
      se_a = s * sqrt(1_wp / n + x_mean**2/ss_xx)
    end if

    res%coeff = [a, b]
    res%stderr = [se_a, se_b]
    res%corr = sqrt(r2)

    return

  end function fit_least_squares

end module da_analyze_mod
