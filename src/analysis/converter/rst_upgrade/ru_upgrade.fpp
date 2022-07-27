!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ru_upgrade_mod
!> @brief   upgrade restart files
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ru_upgrade_mod

  use sp_parallel_io_old_mod
  use sp_parallel_io_mod
  use sp_domain_str_mod
  use sp_enefunc_str_mod
  use sp_boundary_str_mod
  use sp_constraints_str_mod
  use sp_dynvars_str_mod
  use sp_dynamics_str_mod
  use output_str_mod
  use input_str_mod
  use fileio_rst_old_mod
  use fileio_rst_mod
  use fileio_mod
  use messages_mod
  use constants_mod
  use mpi_parallel_mod
 
  implicit none
  private

  ! subroutines
  public  :: upgrade
  private :: upgrade_single
  private :: upgrade_parallel
  private :: setup_num_cell_local

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    upgrade
  !> @brief        upgrade restart file
  !! @authors      NT
  !! @param[in]    input  : input information
  !! @param[in]    output : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine upgrade(input, output)

    ! formal arguments
    type(s_input),           intent(in)    :: input
    type(s_output),          intent(in)    :: output


    if (index(input%rstfile, '(') == 0) then

      call upgrade_single  (input%rstfile, output%rstfile)

    else

      call upgrade_parallel(input%rstfile, output%rstfile)

    end if

    return

  end subroutine upgrade

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine upgrade_single(infile, outfile)

    ! formal arguments
    character(*),            intent(in)    :: infile
    character(*),            intent(in)    :: outfile

    ! local variables
    type(s_rst)              :: rst


    write(MsgOut,'(a)')   'Upgrade_Single Restart file> '
    write(MsgOut,'(a,a)') '    input file  = ', trim(infile)
    write(MsgOut,'(a,a)') '    output file = ', trim(outfile)
    write(MsgOut,'(a)')   ' '


    call input_rst_old(infile, rst)

    
    if (outfile == '') then
      call output_rst(infile, rst)
    else
      call output_rst(outfile, rst)
    end if

    return

  end subroutine upgrade_single

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine upgrade_parallel(infile_rnk, outfile_rnk)

    ! formal arguments
    character(*),            intent(in)    :: infile_rnk
    character(*),            intent(in)    :: outfile_rnk

    ! local variables
    type(s_boundary)         :: boundary
    type(s_domain)           :: domain
    type(s_enefunc)          :: enefunc
    type(s_constraints)      :: constraints
    type(s_dynvars)          :: dynvars
    type(s_dynamics)         :: dynamics
    integer                  :: nplace
    character(1000)          :: infile, outfile
    logical                  :: bex


    write(MsgOut,'(a)')   'Upgrade_Parallel Restart file> '
    write(MsgOut,'(a)')   ' '


    my_world_rank = 0
    my_city_rank = 0

    do nplace = 1, 7
      infile = pio_get_ranked_filename(infile_rnk, nplace)
      inquire(file=infile, exist=bex)
      if (bex) &
        exit
    end do

    do while(.true.)

      call init_boundary(boundary)
      call init_domain(domain)
      call init_enefunc(enefunc)
      call init_constraints(constraints)
      call init_dynvars(dynvars)
      call init_dynamics(dynamics)

      ! read old restart file
      !
      infile = pio_get_ranked_filename(infile_rnk, nplace)
      inquire(file=infile, exist=bex)
      if (.not. bex) &
        exit

      call pio_read_domain_rst_old(infile,      &
                                   boundary,    &
                                   domain,      &
                                   enefunc,     &
                                   constraints, &
                                   dynvars,     &
                                   dynamics)


      ! setup num_cell_local
      !
      call setup_num_cell_local(boundary, domain)


      ! write new restart file

      outfile = pio_get_ranked_filename(outfile_rnk, nplace)

      call pio_write_domain_rst(   outfile,     &
                                   boundary,    &
                                   domain,      &
                                   enefunc,     &
                                   constraints, &
                                   dynvars,     &
                                   dynamics)


      write(MsgOut,'(a,i0)')'  rank : ', my_world_rank
      write(MsgOut,'(a,a)') '    input file  = ', trim(infile)
      write(MsgOut,'(a,a)') '    output file = ', trim(outfile)
      if (pio_restart) then
        write(MsgOut,'(a)') '        restart : Yes'
      else
        write(MsgOut,'(a)') '        restart : No'
      end if
      write(MsgOut,'(a)')   ' '

      my_world_rank = my_world_rank + 1
      my_city_rank = my_city_rank + 1

    end do

    return

  end subroutine upgrade_parallel

  !======1=========2=========3=========4=========5=========6=========7=========8

  !!!! note: part of sp_domain.fpp :: setup_processor_rank()

  subroutine setup_num_cell_local(boundary, domain)

    ! formal arguments
    type(s_boundary),target, intent(in)    :: boundary
    type(s_domain),  target, intent(inout) :: domain

    ! local variables
    integer                  :: iproc(3), cell(3)
    integer                  :: i, quotient, remainder

    integer, pointer         :: ncell_local, num_domain(:)


    num_domain     => boundary%num_domain
    ncell_local    => domain%num_cell_local

    ! Assign the rank of each dimension from my_rank
    !
    iproc(1) = mod(my_city_rank, num_domain(1))
    iproc(2) = mod(my_city_rank/num_domain(1), num_domain(2))
    iproc(3) = my_city_rank/(num_domain(1)*num_domain(2))

    ! Cell number of each dimension
    !
    cell(1) = boundary%num_cells_x
    cell(2) = boundary%num_cells_y
    cell(3) = boundary%num_cells_z

    ! Default value of the number of cell in each domain
    !
    ncell_local = 1

    ! Assign the cell index for each processor and the total number of cells
    !
    do i = 1, 3

      quotient = cell(i) / num_domain(i)
      remainder = mod(cell(i), num_domain(i))

      if (iproc(i) <= (remainder -1)) then
        quotient       = quotient + 1
        ncell_local    = ncell_local * quotient
      else
        ncell_local    = ncell_local * quotient
      end if

    end do

    return

  end subroutine setup_num_cell_local

end module ru_upgrade_mod
