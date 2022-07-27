!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_table_mod
!> @brief   read energy data from table file
!! @authors Yuji Sugita (YS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fileio_table_mod

  use fileio_data_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_table
    real(wp)                 :: cutoffdist
    real(wp)                 :: switchdist
    real(wp)                 :: dielec_const
    real(wp)                 :: table_density
    integer                  :: table_order
    real(wp)                 :: pme_alpha
    real(wp),    allocatable :: table_ene(:)
    real(wp),    allocatable :: table_grad(:)
    real(wp),    allocatable :: table_ene_lj4(:)
    real(wp),    allocatable :: table_grad_lj4(:)
    real(wp),    allocatable :: table_ecor(:)
    real(wp),    allocatable :: table_decor(:)
  end type s_table

  ! subroutines
  public  :: input_table
  public  :: output_table
  public  :: check_table
  public  :: get_table
  public  :: alloc_table
  public  :: dealloc_table

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_table
  !> @brief        a driver subroutine for reading table file
  !! @authors      NT
  !! @param[in]    table_filename : filename of table file
  !! @param[out]   table          : structure of table information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_table(table_filename, table)

    ! formal arguments
    character(*),            intent(in)    :: table_filename
    type(s_table),           intent(inout) :: table

    ! local variables
    integer                  :: f, var_size, alloc_stat


    ! open table file
    call open_data(table_filename, IOFileDataRead, f)

    ! allocate memory
    call get_data_size(f, 'table_decor', var_size)

    if (allocated(table%table_ene)) then
      deallocate(table%table_ene,      &
                 table%table_grad,     &
                 table%table_ene_lj4,  &
                 table%table_grad_lj4, &
                 table%table_ecor,     &
                 table%table_decor)
    end if

    allocate(table%table_ene     (6*var_size), &
             table%table_grad    (6*var_size), &
             table%table_ene_lj4 (2*var_size), &
             table%table_grad_lj4(2*var_size), &
             table%table_ecor    (  var_size), &
             table%table_decor   (  var_size), &
             stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg('Input_Table> memory allocation error.')

    ! read data
    call read_data_real_wp(f, 'cutoffdist',    table%cutoffdist)

    call read_data_real_wp(f, 'switchdist',    table%switchdist)

    call read_data_real_wp(f, 'dielec_const',  table%dielec_const)

    call read_data_real_wp(f, 'table_density', table%table_density)

    call read_data_integer(f, 'table_order',   table%table_order)

    call read_data_real_wp(f, 'pme_alpha',     table%pme_alpha)

    call read_data_real_wp_array(f, 'table_ene',  &
                               (/6*var_size/), table%table_ene)
    call read_data_real_wp_array(f, 'table_grad', &
                               (/6*var_size/), table%table_grad)
    call read_data_real_wp_array(f, 'table_ene_lj4',  &
                               (/2*var_size/), table%table_ene_lj4)
    call read_data_real_wp_array(f, 'table_grad_lj4', &
                               (/2*var_size/), table%table_grad_lj4)
    call read_data_real_wp_array(f, 'table_ecor',  &
                               (/  var_size/), table%table_ecor)
    call read_data_real_wp_array(f, 'table_decor', &
                               (/  var_size/), table%table_decor)

    ! close table file
    call close_data(f)

    return

  end subroutine input_table

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_table
  !> @brief        a driver subroutine for writing table file
  !! @authors      NT
  !! @param[in]    table_filename : filename of table file
  !! @param[in]    use_ascii      : flag for output in ascii format
  !! @param[in]    table          : structure of table information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_table(table_filename, use_ascii, table)

    ! formal arguments
    character(*),            intent(in)    :: table_filename
    logical,                 intent(in)    :: use_ascii
    type(s_table),           intent(in)    :: table

    ! local variables
    integer                  :: f, var_size, mode


    var_size = 0
    if (allocated(table%table_ecor)) var_size = size(table%table_ecor)

    ! open table file
    mode = IOFileDataWriteAscii
    if (.not. use_ascii) mode = IOFileDataWrite

    call open_data(table_filename, mode, f)

    ! write data
    call write_data_real_wp(f, 'cutoffdist',    table%cutoffdist)

    call write_data_real_wp(f, 'switchdist',    table%switchdist)

    call write_data_real_wp(f, 'dielec_const',  table%dielec_const)

    call write_data_real_wp(f, 'table_density', table%table_density)

    call write_data_integer(f, 'table_order',   table%table_order)

    call write_data_real_wp(f, 'pme_alpha',     table%pme_alpha)

    call write_data_real_wp_array(f, 'table_ene',  &
                                  (/6*var_size/), table%table_ene)
    call write_data_real_wp_array(f, 'table_grad', &
                                  (/6*var_size/), table%table_grad)
    call write_data_real_wp_array(f, 'table_ene_lj4',  &
                                  (/2*var_size/), table%table_ene_lj4)
    call write_data_real_wp_array(f, 'table_grad_lj4', &
                                  (/2*var_size/), table%table_grad_lj4)
    call write_data_real_wp_array(f, 'table_ecor',  &
                                  (/  var_size/), table%table_ecor)
    call write_data_real_wp_array(f, 'table_decor', &
                                  (/  var_size/), table%table_decor)

    ! close table file
    call close_data(f)

    return

  end subroutine output_table

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_table
  !> @brief        check table value
  !! @authors      NT
  !! @param[in]    table         : structure of table information
  !! @param[in]    cutoffdist    : cutoff distance
  !! @param[in]    switchdist    : switch distance
  !! @param[in]    dielec_const  : dielectric constant
  !! @param[in]    table_density : table density
  !! @param[in]    table_order   : table order
  !! @param[in]    pme alpha     : pme alpha
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function check_table(table,         &
                       cutoffdist,    &
                       switchdist,    &
                       dielec_const,  &
                       table_density, &
                       table_order,   &
                       pme_alpha)

    ! return value
    logical :: check_table

    ! formal arguments
    type(s_table),           intent(in)    :: table
    real(wp),                intent(in)    :: cutoffdist
    real(wp),                intent(in)    :: switchdist
    real(wp),                intent(in)    :: dielec_const
    real(wp),                intent(in)    :: table_density
    integer,                 intent(in)    :: table_order
    real(wp),                intent(in)    :: pme_alpha


    check_table = table%cutoffdist    == cutoffdist    .and. &
                  table%switchdist    == switchdist    .and. &
                  table%dielec_const  == dielec_const  .and. &
                  table%table_density == table_density .and. &
                  table%table_order   == table_order   .and. &
                  table%pme_alpha     == pme_alpha
    return

  end function check_table

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_table
  !> @brief        get table array variables
  !! @authors      NT
  !! @param[in]    table : structure of table information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine get_table(table,          &
                       table_ene,      &
                       table_grad,     &
                       table_ene_lj4,  &
                       table_grad_lj4, &
                       table_ecor,     &
                       table_decor)

    ! formal arguments
    type(s_table),           intent(in)    :: table
    real(wp),                intent(inout) :: table_ene(:)
    real(wp),                intent(inout) :: table_grad(:)
    real(wp),      optional, intent(inout) :: table_ene_lj4(:)
    real(wp),      optional, intent(inout) :: table_grad_lj4(:)
    real(wp),      optional, intent(inout) :: table_ecor(:)
    real(wp),      optional, intent(inout) :: table_decor(:)

    ! local variables
    integer                  :: var_size


    if (present(table_ecor) .and. present(table_decor)) then

      if (size(table%table_ene)      /= size(table_ene)      .or. &
          size(table%table_grad)     /= size(table_grad)     .or. &
          size(table%table_ene_lj4)  /= size(table_ene_lj4)  .or. &
          size(table%table_grad_lj4) /= size(table_grad_lj4) .or. &
          size(table%table_ecor)     /= size(table_ecor)     .or. &
          size(table%table_decor)    /= size(table_decor)) &
        call error_msg('Get_Table> table size is different.')


      var_size = size(table%table_ecor)

      table_ene     (1:6*var_size) = table%table_ene     (1:6*var_size)
      table_grad    (1:6*var_size) = table%table_grad    (1:6*var_size)
      table_ene_lj4 (1:2*var_size) = table%table_ene_lj4 (1:2*var_size)
      table_grad_lj4(1:2*var_size) = table%table_grad_lj4(1:2*var_size)
      table_ecor    (1:  var_size) = table%table_ecor    (1:  var_size)
      table_decor   (1:  var_size) = table%table_decor   (1:  var_size)

    else if (present(table_ene_lj4) .and. present(table_grad_lj4)) then

      if (size(table%table_ene)      /= size(table_ene)      .or. &
          size(table%table_grad)     /= size(table_grad)     .or. &
          size(table%table_ene_lj4)  /= size(table_ene_lj4)  .or. &
          size(table%table_grad_lj4) /= size(table_grad_lj4)) &
        call error_msg('Get_Table> table size is different.')


      var_size = size(table%table_ene)/6

      table_ene     (1:6*var_size) = table%table_ene     (1:6*var_size)
      table_grad    (1:6*var_size) = table%table_grad    (1:6*var_size)
      table_ene_lj4 (1:2*var_size) = table%table_ene_lj4 (1:2*var_size)
      table_grad_lj4(1:2*var_size) = table%table_grad_lj4(1:2*var_size)

    else 
      if (size(table%table_ene)      /= size(table_ene)      .or. &
          size(table%table_grad)     /= size(table_grad)) &
        call error_msg('Get_Table> table size is different.')


      var_size = size(table%table_ene)/6

      table_ene     (1:6*var_size) = table%table_ene     (1:6*var_size)
      table_grad    (1:6*var_size) = table%table_grad    (1:6*var_size)

    end if

    return

  end subroutine get_table

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_table
  !> @brief        allocate table information
  !! @authors      NT
  !! @param[inout] table    : structure of table information
  !! @param[in]    var_size : variable size
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_table(table, var_size)

    ! formal arguments
    type(s_table),           intent(inout) :: table
    integer,                 intent(in)    :: var_size

    ! local variables
    integer                  :: alloc_stat


    if (allocated(table%table_ene)) then
      deallocate(table%table_ene,      &
                 table%table_grad,     &
                 table%table_ene_lj4,  &
                 table%table_grad_lj4, &
                 table%table_ecor,     &
                 table%table_decor)
    endif

    allocate(table%table_ene     (6*var_size),  &
             table%table_grad    (6*var_size),  &
             table%table_ene_lj4 (2*var_size),  &
             table%table_grad_lj4(2*var_size),  &
             table%table_ecor    (  var_size),  &
             table%table_decor   (  var_size),  &
             stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg('Alloc_Table> memory allocation error.')

    table%table_ene     (1:6*var_size) = 0.0_wp
    table%table_grad    (1:6*var_size) = 0.0_wp
    table%table_ene_lj4 (1:2*var_size) = 0.0_wp
    table%table_grad_lj4(1:2*var_size) = 0.0_wp
    table%table_ecor    (1:  var_size) = 0.0_wp
    table%table_decor   (1:  var_size) = 0.0_wp

    return

  end subroutine alloc_table

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_table
  !> @brief        deallocate table information
  !! @authors      NT
  !! @param[inout] table : structure of table information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_table(table)

    ! formal arguments
    type(s_table),           intent(inout) :: table


    if (allocated(table%table_ene)) then
      deallocate(table%table_ene,      &
                 table%table_grad,     &
                 table%table_ene_lj4,  &
                 table%table_grad_lj4, &
                 table%table_ecor,     &
                 table%table_decor)
    endif

    return

  end subroutine dealloc_table

end module fileio_table_mod
