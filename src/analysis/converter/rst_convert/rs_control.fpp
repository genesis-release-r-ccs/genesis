!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   rs_control_mod
!> @brief   read parameters and data for restart file conversion
!! @authors Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module rs_control_mod

  use rs_option_mod
  use output_mod
  use input_mod
  use fileio_control_mod
  use string_mod
  use messages_mod

  implicit none
  private

  ! structures
  type, public :: s_ctrl_data

    ! data for section input
    type(s_inp_info)      :: inp_info

    ! data for section output
    type(s_out_info)      :: out_info

    ! data for section option
    type(s_opt_info)      :: opt_info

  end type s_ctrl_data

  ! subroutines
  public  :: usage
  public  :: control

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    usage
  !> @brief        show usage of RST_CONVERT
  !! @authors      NT
  !! @param[in]    arg1 : arguments in execution
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine usage(arg1)

    ! formal arguments
    character(*),            intent(out) :: arg1

    ! local variables
    character(Maxline)       :: arg2
    integer                  :: iargc


    call getarg(1, arg1)

    if (iargc() < 1      .or. &
         arg1 == '-help' .or. &
         arg1 == '-h'    .or. &
         arg1 == '-HELP' .or. &
         arg1 == '-H') then

      call getarg(2, arg2)
      call tolower(arg2)

      if (iargc() < 2 .or. arg2 /= 'ctrl') then
        write(MsgOut,'(A)') ' '
        write(MsgOut,'(A)') '# normal usage'
        write(MsgOut,'(A)') '  % ./rst_convert INPRST'
        write(MsgOut,'(A)') ' '
        write(MsgOut,'(A)') '# usage to see control parameters'
        write(MsgOut,'(A)') '  % ./rst_convert -h ctrl'
        write(MsgOut,'(A)') ' '

      else

        write(MsgOut,'(A)') '# control parameters in rst_convert'
        write(MsgOut,'(A)') ' '

        call show_ctrl_input ('coor,vel,xsc,rst,pdb,psf,ambcrd,prmtop,grotop,grocrd,RST')
        call show_ctrl_output('rst,coor,vel,xsc,pdb,PDB,crd,CRD,ambcrd,AMBCRD')
        call show_ctrl_option

      end if

      stop

    else
      write(MsgOut,'(A)')'****************************************************'
      write(MsgOut,'(A)')'*                                                  *'
      write(MsgOut,'(A)')'*        RST_CONVERT: convert restart files        *'
      write(MsgOut,'(A)')'*                                                  *'
      write(MsgOut,'(A)')'*              Developed by RIKEN                  *'
      write(MsgOut,'(A)')'*                                                  *'
      write(MsgOut,'(A)')'****************************************************'
      write(MsgOut,'(A)')' '
    end if

    return
  
  end subroutine usage

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    control
  !> @brief        open/read/close control files 
  !! @authors      NT
  !! @param[in]    filename  : file name of control file
  !! @param[inout] ctrl_data : information of control parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine control(filename, ctrl_data)
  
    ! formal arguments
    character(*),            intent(in)    :: filename
    type(s_ctrl_data),       intent(inout) :: ctrl_data
    
    ! local variables
    integer                  :: handle
    

    ! open control file
    ! 
    if ( len_trim(filename) > 0 ) then
      call open_ctrlfile(filename, handle)

      if (handle == 0) &
        call error_msg('Control_Rstcnv> File IO Error')


      ! read input section
      !
      call read_ctrl_input (handle, ctrl_data%inp_info)

      ! read output section
      !
      call read_ctrl_output(handle, ctrl_data%out_info)

      ! read option section
      !
      call read_ctrl_option(handle, ctrl_data%opt_info)

      ! close control file
      !
      call close_ctrlfile(handle)

    end if

    return

  end subroutine control

end module rs_control_mod
