!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_top_mod
!> @brief   read a topology file (CHARMM-style)
!! @authors Yuji Sugita (YS), Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fileio_top_mod

  use fileio_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_top

    integer           :: num_atom_cls = 0
    integer           :: num_res_type = 0

    ! atom class
    integer,          allocatable :: atom_cls_no(:)
    character(6),     allocatable :: atom_cls_name(:)
    real(wp),         allocatable :: atom_cls_mass(:)
    character(2),     allocatable :: atom_ele_name(:)
    character(60),    allocatable :: atom_cls_comment(:)

    ! residue types
    character(6),     allocatable :: res_name(:)
    integer,          allocatable :: res_atom_cls_bgn(:)
    integer,          allocatable :: res_atom_cls_end(:)

  end type s_top

  ! parameters for allocatable variables
  integer,      public, parameter :: TopAtomCls = 1
  integer,      public, parameter :: TopResType = 2

  ! parameters
  integer,     private, parameter :: MAXROW = 80

  ! local variables
  logical,                private :: vervose = .true.

  ! subroutines
  public  :: input_top
  public  :: init_top
  public  :: alloc_top
  public  :: dealloc_top
  public  :: dealloc_top_all
  public  :: read_top
  public  :: merge_top
  private :: copy_top

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_top
  !> @brief        open, read, and close a topfile 
  !! @authors      YS, NT
  !! @param[in]    top_filename : filename of a topfile
  !! @param[out]   top          : CHARMM TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_top(top_filename, top)

    ! formal arguments
    character(*),            intent(in)    :: top_filename
    type(s_top),             intent(inout) :: top

    ! local variables
    type(s_top)              :: top0
    integer                  :: unit_no, i
    character(MaxFilename)   :: filename


    call init_top(top)

    i = 0
    do while(extract(top_filename, i, filename))

      ! open a topology file
      !
      call open_file(unit_no, filename, IOFileInput)

      ! read only the information of weight and class for each atomtype
      !
      call read_top(unit_no, top0)

      ! close a topology file
      !
      call close_file(unit_no)

      ! merge topology information
      !
      call merge_top(top0, top)

    end do

    if (main_rank .and. vervose) then 
      write(MsgOut,'(A)') 'Input_Top> Summary of Topfile'
      write(MsgOut,'(A20,I10,A20,I10)') &
           '  num_atom_class  = ', top%num_atom_cls, &
           '  num_resi_type   = ', top%num_res_type
      write(MsgOut,'(A)') ' '
      vervose = .false.
    end if

    return

  end subroutine input_top

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_top
  !> @brief        initialize CHARMM TOP information
  !! @authors      YS
  !! @param[out]   top : CHARMM TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_top(top)

    ! formal arguments
    type(s_top),             intent(inout) :: top


    top%num_atom_cls = 0
    top%num_res_type = 0

    return

  end subroutine init_top

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_top
  !> @brief        allocate CHARMM TOP information
  !! @authors      YS
  !! @param[inout] top      : CHARMM TOP information
  !! @param[in]    variable : an variable to be allocated (TopAtomCls)
  !! @param[in]    var_size : allocation size
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_top(top, variable, var_size)

    ! formal arguments
    type(s_top),             intent(inout) :: top
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size
    
    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat

    
    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case(TopAtomCls)

      if (allocated(top%atom_cls_no)) then
        if (size(top%atom_cls_no) == var_size) &
          return
        deallocate(top%atom_cls_no,      &
                   top%atom_cls_name,    &
                   top%atom_cls_mass,    &
                   top%atom_ele_name,    &
                   top%atom_cls_comment, &
                   stat = dealloc_stat)
      end if

      allocate(top%atom_cls_no(var_size),      &
               top%atom_cls_name(var_size),    &
               top%atom_cls_mass(var_size),    &
               top%atom_ele_name(var_size),    &
               top%atom_cls_comment(var_size), &
               stat = alloc_stat)

      top%atom_cls_no     (1:var_size) = 0
      top%atom_cls_name   (1:var_size) = ''
      top%atom_cls_mass   (1:var_size) = 0.0_wp
      top%atom_ele_name   (1:var_size) = ''
      top%atom_cls_comment(1:var_size) = ''

    case (TopResType)

      if (allocated(top%res_name)) then
        if (size(top%res_name) == var_size) &
          return
        deallocate(top%res_name,         &
                   top%res_atom_cls_bgn, &
                   top%res_atom_cls_end, &
                   stat = dealloc_stat)
      end if

      allocate(top%res_name(var_size),         &
               top%res_atom_cls_bgn(var_size), &
               top%res_atom_cls_end(var_size), &
               stat = alloc_stat)

      top%res_name        (1:var_size) = ''
      top%res_atom_cls_bgn(1:var_size) = 0
      top%res_atom_cls_end(1:var_size) = 0

    case default

      call error_msg('Alloc_Top> bad variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_top

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_top
  !> @brief        deallocate CHARMM TOP information
  !! @authors      YS
  !! @param[inout] top      : CHARMM TOP information
  !! @param[in]    variable : an variable to be allocated (TopAtomCls)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_top(top, variable)

    ! formal arguments
    type(s_top),             intent(inout) :: top
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat
    

    dealloc_stat = 0

    ! deallocate selected variable
    !
    select case (variable)

    case(TopAtomCls)

      if (allocated(top%atom_cls_no)) then
        deallocate (top%atom_cls_no,      &
                    top%atom_cls_name,    &
                    top%atom_cls_mass,    &
                    top%atom_ele_name,    &
                    top%atom_cls_comment, &
                    stat = dealloc_stat)
      end if

    case (TopResType)

      if (allocated(top%res_name)) then
        deallocate (top%res_name,         &
                    top%res_atom_cls_bgn, &
                    top%res_atom_cls_end, &
                    stat = dealloc_stat)
       end if

    case default

      call error_msg('Dealloc_Top> bad variable')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_top

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_top_all
  !> @brief        deallocate all CHARMM TOP iformation
  !! @authors      YS
  !! @param[inout] top : CHARMM TOP information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_top_all(top)

    ! formal arguments
    type(s_top),             intent(inout) :: top


    call dealloc_top(top, TopAtomCls)
    call dealloc_top(top, TopResType)

    return

  end subroutine dealloc_top_all

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_top
  !> @brief        read data from a topfile
  !! @authors      YS
  !! @param[in]    unit_no : unit number of a topfile
  !! @param[out]   top     : CHARMM TOP inforation
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_top(unit_no, top)

    ! formal arguments
    integer,                 intent(in)    :: unit_no
    type(s_top),             intent(out)   :: top

    ! local variables
    integer                  :: ia, ir, nchar, nsta, nend, ndata
    integer                  :: num_cls, atom_code, num_res_type
    integer                  :: i, read_cnt
    character(80)            :: line, linebk


    ! initialize variables
    !
    call init_top(top)

    read_cnt     = 0
    num_cls      = 0
    num_res_type = 0

    !  read a topfile and check the number of atom classes
    !
    do while(.true.)
      
      read(unit_no, '(A80)') line
      read_cnt = read_cnt + 1

      call char_line(MAXROW, line, nchar)

      if (line(1:4) .eq. 'MASS') then

        !  Read MASS section 
        !
        nsta = 5
        nend = nchar

        call read_int(line, nsta, nend, atom_code)

        num_cls = num_cls + 1

      else if (line(1:4) .eq. 'RESI' .or. line(1:4) .eq. 'PRES') then

        !  read RESI section
        !
        num_res_type = num_res_type + 1

      else if (line(1:3) .eq. 'end' .or. line(1:3) .eq. 'END') then

        exit

      end if

    end do


    !  allocate top   
    !
    call alloc_top(top, TopAtomCls, num_cls)
    call alloc_top(top, TopResType, num_res_type)


    !  read a topfile and build top
    !
    do i = 1, read_cnt
      backspace(unit_no)
    end do

    ia = 0
    ir = 0

    do while(.true.)
       
      read(unit_no, '(a80)') line

      linebk = line

      call char_line(MAXROW, line, nchar)

      if (line(1:4) .eq. 'MASS') then

        !  read atom_class section and store the information
        !
        call read_ndata(line, 1, nchar, ndata)

        ia   = ia + 1
        nsta = 5
        nend = nchar

        call read_int (line, nsta, nend, top%atom_cls_no(ia))
        call read_word(line, nsta, nend, 6, top%atom_cls_name(ia))
        call read_real(line, nsta, nend, top%atom_cls_mass(ia))
        if (ndata == 5) then
          call read_word(line, nsta, nend, 2, top%atom_ele_name(ia))
        else
          top%atom_ele_name(ia) = 'X '
        end if
        call read_comment(linebk, 1, MAXROW, top%atom_cls_comment(ia))

      else if (line(1:4) .eq. 'RESI' .or. line(1:4) .eq. 'PRES') then

        !  read RESI section
        !
        ir   = ir + 1
        nsta = 5
        nend = nchar

        call read_word(line, nsta, nend, 6, top%res_name(ir))

        top%res_atom_cls_bgn(ir) = 1
        top%res_atom_cls_end(ir) = num_cls

      else if (line(1:3) .eq. 'end' .or. line(1:3) .eq. 'END') then

#ifdef DEBUG
        if (main_rank) then
          write(MsgOut,'(A)') 'Read_Top> Summary of Topfile'
          write(MsgOut,'(A20,I10,A20,I10)')          &
                '  num_atom_class  = ', num_cls,     &
                '  num_resi_type   = ', num_res_type
          write(MsgOut,'(A)') ' '
        end if
#endif

        top%num_atom_cls = num_cls
        top%num_res_type = num_res_type
        exit

      end if

    end do

    return

  end subroutine read_top

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    merge_top
  !> @brief        merge two topology informations
  !! @authors      NT
  !! @param[in]    top0 : CHARMM TOP information of the source
  !! @param[inout] top  : CHARMM TOP information of the destination
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine merge_top(top0, top)

    ! formal arguments
    type(s_top),             intent(in)    :: top0
    type(s_top),             intent(inout) :: top

    ! local variables
    type(s_top)              :: topw
    integer                  :: i, j, ii, n

    
    if (top%num_atom_cls == 0) then

      call alloc_top(top, TopAtomCls, top0%num_atom_cls)
      call alloc_top(top, TopResType, top0%num_res_type)

      top%num_atom_cls = top0%num_atom_cls
      top%num_res_type = top0%num_res_type

      call copy_top(top0, top)

      return

    end if

    topw%num_atom_cls = top%num_atom_cls + top0%num_atom_cls
    topw%num_res_type = top%num_res_type + top0%num_res_type

    call alloc_top(topw, TopAtomCls, topw%num_atom_cls)
    call alloc_top(topw, TopResType, topw%num_res_type)

    call copy_top(top, topw)


    ! check and copy atom class
    !
    do i = 1, top0%num_atom_cls

      ! multiple definition check
      do j = 1, top%num_atom_cls
        if (top0%atom_cls_name(i) .ne. top%atom_cls_name(j)) &
          cycle

        if (top0%atom_cls_mass(i) /= top%atom_cls_mass(j) .or. &
             top0%atom_ele_name(i) .ne. top%atom_ele_name(j)) then

          if (main_rank .and. top0%atom_cls_mass(i) /= top%atom_cls_mass(j)) &
               write(MsgOut,'(3a,f8.4,a,f8.4)') &
               ' Merge_Top> WARNING: Multiple definition: "', &
               trim(top%atom_cls_name(j)), '" Mass: ', &
               top0%atom_cls_mass(i), ' <=> ', top%atom_cls_mass(j)

          if (main_rank .and. top0%atom_ele_name(i) .ne. top%atom_ele_name(j))&
               write(MsgOut,'(6a)') &
               ' Merge_Top> WARNING: Multiple definition: "', &
               trim(top%atom_cls_name(j)), '" Element: ', &
               top0%atom_ele_name(i), ' <=> ', top%atom_ele_name(j)

        else

          if (main_rank) &
            write(MsgOut,'(3a)') ' Merge_Top> WARNING: Multiple definition: "', &
            trim(top%atom_cls_name(j)), '"'

        end if

      end do

      topw%atom_cls_no  (i+top%num_atom_cls) = top0%atom_cls_no  (i)
      topw%atom_cls_name(i+top%num_atom_cls) = top0%atom_cls_name(i)
      topw%atom_cls_mass(i+top%num_atom_cls) = top0%atom_cls_mass(i)
      topw%atom_ele_name(i+top%num_atom_cls) = top0%atom_ele_name(i)

    end do
    
    n = top0%num_atom_cls + top%num_atom_cls
    call dealloc_top(top, TopAtomCls)
    call alloc_top  (top, TopAtomCls, n)


    ! check and copy residue types
    !

    ii = 0
    do i = 1, top0%num_res_type
      do j = 1, top%num_res_type
        if (top0%res_name(i) .ne. top%res_name(j)) &
          cycle
        if (main_rank) &
          write(MsgOut,'(3a)')' Merge_Top> WARNING: RESI: Multiple definition: "', &
          trim(top%res_name(j)), '" : will be ignored.'
        exit
      end do
      if (j <= top%num_res_type) &
        cycle

      ii = ii + 1
      topw%res_name        (ii+top%num_res_type) = top0%res_name(i)
      if (top0%num_atom_cls == 0) then
        !
        ! read all atom cls
        !
        topw%res_atom_cls_bgn(ii+top%num_res_type) = 1 
        topw%res_atom_cls_end(ii+top%num_res_type) = top%num_atom_cls
      else
        topw%res_atom_cls_bgn(ii+top%num_res_type) = &
                          top0%res_atom_cls_bgn(i) + top%num_atom_cls
        topw%res_atom_cls_end(ii+top%num_res_type) = &
                          top0%res_atom_cls_end(i) + top%num_atom_cls
      end if
    end do

    n = ii + top%num_res_type
    call dealloc_top(top, TopResType)
    call alloc_top  (top, TopResType, n)


    top%num_atom_cls = size(top%atom_cls_no)
    top%num_res_type = size(top%res_name)

    call copy_top(topw, top)

    call dealloc_top_all(topw)

    return

  end subroutine merge_top

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    copy_top
  !> @brief        copy topology informations
  !! @authors      NT
  !! @param[in]    top_src : CHARMM TOP information of the source
  !! @param[inout] top_dst : CHARMM TOP information of the destination
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine copy_top(top_src, top_dst)

    ! formal arguments
    type(s_top),             intent(in)    :: top_src
    type(s_top),             intent(inout) :: top_dst

    ! local variables
    integer                  :: n


    n = min(top_dst%num_atom_cls, top_src%num_atom_cls)
    top_dst%atom_cls_no  (1:n) = top_src%atom_cls_no  (1:n)
    top_dst%atom_cls_name(1:n) = top_src%atom_cls_name(1:n)
    top_dst%atom_cls_mass(1:n) = top_src%atom_cls_mass(1:n)
    top_dst%atom_ele_name(1:n) = top_src%atom_ele_name(1:n)

    n = min(top_dst%num_res_type, top_src%num_res_type)
    top_dst%res_name        (1:n) = top_src%res_name        (1:n)
    top_dst%res_atom_cls_bgn(1:n) = top_src%res_atom_cls_bgn(1:n)
    top_dst%res_atom_cls_end(1:n) = top_src%res_atom_cls_end(1:n)

    return

  end subroutine copy_top

end module fileio_top_mod
