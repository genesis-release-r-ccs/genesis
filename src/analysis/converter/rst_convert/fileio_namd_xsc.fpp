!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_namd_xsc_mod
!> @brief   NAMD xsc file I/O
!! @authors Norio Takase (NT), Chigusa Kobayashi (CK)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module fileio_namd_xsc_mod

  use fileio_mod
  use string_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_namd_xsc
    integer    :: step
    real(dp)   :: a(3)
    real(dp)   :: b(3)
    real(dp)   :: c(3)
    real(dp)   :: o(3)
    real(dp)   :: s(6)
  end type s_namd_xsc

  ! subroutines
  public  :: input_namd_xsc
  public  :: output_namd_xsc
  private :: output_format

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_namd_xsc
  !> @brief        input NAMD extended system
  !! @authors      NT
  !! @param[in]    filename : file name
  !! @param[inout] xsc      : extended system information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_namd_xsc(filename, xsc)

    ! formal arguments
    character(*),            intent(in)    :: filename
    type(s_namd_xsc),        intent(inout) :: xsc

    ! local varialbes
    integer                  :: xsc_file
    character(MaxLine)       :: line


    call open_file(xsc_file, filename, IOFileInput)

    do while (.true.)
      read(xsc_file,'(A)',end=900) line
      if (index(line, '#') == 0) &
        exit
    end do

    xsc%step = 0
    xsc%a(:) = 0.0_wp
    xsc%b(:) = 0.0_wp
    xsc%c(:) = 0.0_wp
    xsc%o(:) = 0.0_wp
    xsc%s(:) = 0.0_wp

    read(line,*,err=900,end=100) xsc%step, &
                                 xsc%a,    &
                                 xsc%b,    &
                                 xsc%c,    &
                                 xsc%o,    &
                                 xsc%s

100 write(MsgOut,'(A,I8)' )    'Input_Namd_Xsc> step : ', xsc%step
    write(MsgOut,'(A,3F8.3)')  '                   a : ', xsc%a
    write(MsgOut,'(A,3F8.3)')  '                   b : ', xsc%b
    write(MsgOut,'(A,3F8.3)')  '                   c : ', xsc%c
    write(MsgOut,'(A,3F8.3)')  '                   o : ', xsc%o
    write(MsgOut,'(A,3F16.9)') '                   s : ', xsc%s(1:3)
    write(MsgOut,'(A,3F16.9)') '                       ', xsc%s(4:6)
    write(MsgOut,'(A)')        ' '

    call close_file(xsc_file)

    return

900 call error_msg('Input_Namd_Xsc> File format error.:')

  end subroutine input_namd_xsc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_namd_xsc
  !> @brief        output NAMD extended system
  !! @authors      CK
  !! @param[in]    filename : file name
  !! @param[inout] xsc      : extended system information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_namd_xsc(filename, xsc)

    ! formal arguments
    character(*),            intent(in)    :: filename
    type(s_namd_xsc),        intent(in)    :: xsc

    ! local varialbes
    integer                  :: xsc_file

    call open_file(xsc_file, filename, IOFileOutputNew)

    write(xsc_file, '(A)')   &
      "# NAMD extended system configuration restart file"
    write(xsc_file, '(A)')   &
      "#$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z s_x s_y s_z s_u s_v s_w"
    write(xsc_file, '(i5,$)') xsc%step
    call output_format(xsc_file, xsc%a)
    call output_format(xsc_file, xsc%b)
    call output_format(xsc_file, xsc%c)
    call output_format(xsc_file, xsc%o)
    call output_format(xsc_file, xsc%s)
    call close_file(xsc_file)

    return

  end subroutine output_namd_xsc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_format
  !> @brief        output NAMD extended system with format
  !! @authors      CK
  !! @param[in]    file     : file number
  !! @param[inout] vector   : box information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_format(file, vector)
    integer,          intent(in) :: file
    real(wp),         intent(in) :: vector(:)
    integer                      :: i
    integer                      :: zero

    zero = 0

    do i = 1,size(vector(:))
      if (abs(vector(i)) > 1d-5) then 
        write(file,'(1x,f8.3,$)') vector(i)
      else
        write(file,'(1x,i1,$)')  zero
      endif
    end do

   return

  end subroutine output_format

end module fileio_namd_xsc_mod
