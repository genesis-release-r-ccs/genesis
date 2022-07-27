!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_str_mod
!> @brief   read a CHARMM stream file
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fileio_str_mod

  use fileio_par_mod
  use fileio_top_mod
  use fileio_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! local variables
  logical,                private :: vervose = .true.  

  ! subroutines
  public  :: input_str
  private :: next_stream

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_str
  !> @brief        a driver subroutine for reading CHARMM stream file
  !! @authors      NT
  !! @param[in]    str_filename : filename of CHARMM stream file
  !! @param[in]    top          : CHARMM TOP information
  !! @param[out]   par          : CHARMM PAR information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_str(str_filename, top, par)

    ! formal arguments
    character(*),            intent(in)    :: str_filename
    type(s_top),             intent(inout) :: top
    type(s_par),             intent(inout) :: par
  
    ! local variables
    type(s_par)              :: par0
    type(s_top)              :: top0
    integer                  :: file, i
    character(MaxFilename)   :: filename, str_name


    i  = 0
    do while(extract(str_filename, i, filename))

      ! open CHARMM stream file
      !
      call open_file(file, filename, IOFileInput)


      do while(next_stream(file, str_name))

        if (str_name .eq. 'rtf') then

          ! read CHARMM topology file
          !
          call read_top(file, top0)

          ! merge topology information
          !
          call merge_top(top0, top)

        else if (str_name(1:4) .eq. 'para') then

          ! read CHARMM parameter file
          !
          call read_par(file, par0)

          ! resolve CHARMM19 nonbonded wild-card atom name
          !
          call resolve_wildcard(par0, top)

          ! merge parameter information
          !
          call merge_par(par0, par)

        else

          if (main_rank) &
            write(MsgOut,'(A)') 'Input_Str> Unknown stream :'//trim(str_name)

        end if
          
      end do


      ! close CHARMM stream file
      !
      call close_file(file)

    end do


    if (main_rank .and. vervose) then
      write(MsgOut,'(A)') 'Input_Str> Summary of Top information'
      write(MsgOut,'(A20,I10,A20,I10)') &
           '  num_atom_class  = ', top%num_atom_cls, &
           '  num_resi_type   = ', top%num_res_type
      write(MsgOut,'(A)') ' '
    end if

    if (main_rank .and. vervose) then
      write(MsgOut,'(A)') 'Input_Str> Summary of Par information'
      write(MsgOut,'(A20,I10,A20,I10)')               &
           '  num_bonds       = ', par%num_bonds,     & 
           '  num_angles      = ', par%num_angles
      write(MsgOut,'(A20,I10,A20,I10)')               &
           '  num_dihedrals   = ', par%num_dihedrals, &
           '  num_impropers   = ', par%num_impropers
      write(MsgOut,'(A20,I10,A20,I10)')               & 
           '  num_atom_cls    = ', par%num_atom_cls,  &
           '  num_nbfix       = ', par%num_nbfix
      write(MsgOut,'(A20,I10)')                       & 
           '  num_cmap_terms  = ', par%num_cmaps
      write(MSgOut,'(A)') ' '
    end if
    vervose = .false.


    return

  end subroutine input_str

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      next_stream
  !> @brief        seek to next stream
  !! @authors      NT
  !! @param[in]    unit_no  : file unit number
  !! @param[inout] str_name : stream name
  !! @return       flag for next stream was found
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function next_stream(unit_no, str_name)

    ! return
    logical                  :: next_stream

    ! formal arguments
    integer,                 intent(in)    :: unit_no
    character(*),            intent(inout) :: str_name

    ! local variables
    integer                  :: idx, stat
    character(MaxLine)       :: line
    character(20)            :: command


    do while(.true.)

      read(unit_no,'(a)',iostat=stat) line
      if (stat /= 0) &
        exit

      idx = scan(line, '!')
      if (idx /= 0) line(idx:) = ''

      command  = ''
      str_name = ''
      read(line,*,end=10,err=10) command, str_name
10    call tolower(command)
      call tolower(str_name)

      if (command .eq. 'read') then
        next_stream = .true.
        return
      end if

    end do

    next_stream = .false.
    return

  end function next_stream

end module fileio_str_mod
