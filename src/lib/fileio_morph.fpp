!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_morph_mod
!> @brief   read & write morphing information
!! @authors Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fileio_morph_mod

  use fileio_grotop_mod
  use fileio_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_morph_in

    integer               :: num_morph_bb   = 0
    integer               :: num_morph_sc   = 0

    type(s_morph_pair),allocatable :: morph_bb(:)
    type(s_morph_pair),allocatable :: morph_sc(:)

  end type s_morph_in

  integer,      public, parameter :: MorphBB = 1
  integer,      public, parameter :: MorphSC = 2

  ! subroutines
  public  :: input_morph_in
  public  :: output_morph_in
  public  :: init_morph_in
  public  :: alloc_morph_in
  public  :: dealloc_morph_in
  private :: read_morph
  private :: write_morph
  private :: check_morph

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_morph_in
  !> @brief        a driver subroutine for reading Morph file
  !! @authors      CK
  !! @param[in]    morph_filename : filename of Morph file
  !! @param[out]   morph_in       : structure of Morph information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_morph_in(morph_filename, morph_in)

    ! formal arguments
    character(*),            intent(in)    :: morph_filename
    type(s_morph_in),        intent(inout) :: morph_in
   
    ! local variables
    integer                  :: file
    

    ! open Morph file
    !
    call open_file(file, morph_filename, IOFileInput)

    ! read principal component vector from Morph file
    !
    call read_morph(file, morph_in)

    ! close Morph file
    !
    call close_file(file)

    return

  end subroutine input_morph_in

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_morph_in
  !> @brief        a driver subroutine for writing Morph file
  !! @authors      CK
  !! @param[in]    morph_filename : filename of Morph file
  !! @param[in]    morph_in       : structure of Morph information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_morph_in(morph_filename, morph_in)

    ! formal arguments
    character(*),            intent(in) :: morph_filename
    type(s_morph_in),        intent(in) :: morph_in
   
    ! local variables
    integer                  :: file
    

    ! open Morph file
    !
    call open_file(file, morph_filename, IOFileOutputNew)

    ! write coordinate data from Morph file
    !
    call write_morph(file, morph_in%morph_bb, MorphBB)
    if (allocated(morph_in%morph_sc)) then
      call write_morph(file, morph_in%morph_sc, MorphSC)
    endif

    ! close Morph file
    !
    call close_file(file)

    return

  end subroutine output_morph_in

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_morph_in
  !> @brief        initialize Morph information
  !! @authors      CK
  !! @param[out]   morph_in : structure of Morph information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_morph_in(morph_in)

    ! formal arguments
    type(s_morph_in),             intent(inout) :: morph_in


    morph_in%num_morph_bb         = 0
    morph_in%num_morph_sc         = 0

    return

  end subroutine init_morph_in

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_morph_in
  !> @brief        allocate Morph information
  !! @authors      CK
  !! @param[inout] morph_in  : structure of Morph information
  !! @param[in]    var_size1 : size of the selected variable
  !! @param[in]    var_size2 : size of the selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_morph_in(morph_in, var_size1,var_size2)

    ! formal arguments
    type(s_morph_in),        intent(inout) :: morph_in
    integer,                 intent(in)    :: var_size1
    integer,                 intent(in)    :: var_size2
  
    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat
  

    ! initialize
    !
    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate variable
    !
    if (var_size1 > 0) then

      if (allocated(morph_in%morph_bb)) then
        if (size(morph_in%morph_bb) == var_size1) return
        deallocate (morph_in%morph_bb,  &
                    stat = dealloc_stat)
      end if
      allocate(morph_in%morph_bb(1:var_size1),  &
             stat = alloc_stat)

    endif

    if (var_size2 > 0) then
      if (allocated(morph_in%morph_sc)) then
        if (size(morph_in%morph_sc) == var_size2) return
        deallocate (morph_in%morph_sc,  &
                    stat = dealloc_stat)
      end if
      allocate(morph_in%morph_sc(1:var_size2),  &
             stat = alloc_stat)
    endif

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_morph_in

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_morph_in
  !> @brief        deallocate Morph information
  !! @authors      CK
  !! @param[inout] morph     : structure of Morph information
  !! @param[in]    variable : selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_morph_in(morph_in)

    ! formal arguments
    type(s_morph_in),        intent(inout) :: morph_in

    ! local variables
    integer                  :: dealloc_stat


    ! initialize
    !
    dealloc_stat = 0

    ! deallocate selected variable
    !  
    if (allocated(morph_in%morph_bb)) then
      deallocate (morph_in%morph_bb,  &
                  stat = dealloc_stat)
    end if
    if (allocated(morph_in%morph_sc)) then
      deallocate (morph_in%morph_sc,  &
                  stat = dealloc_stat)
    end if

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_morph_in

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_morph
  !> @brief        read data from Morph file
  !! @authors      CK
  !! @param[in]    file : unit number of Morph file
  !! @param[out]   morph_in : morph data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_morph(file, morph_in)

    ! formal arugments
    integer,                 intent(in)    :: file
    type(s_morph_in),        intent(inout) :: morph_in

    ! local variables
    integer                  :: i
    integer                  :: num_morph_bb, num_morph_sc, num_char
    character(80)            :: line
    logical                  :: backbone
    integer                  :: idxb, idxe
    character(20)            :: dir_str


    ! deallocate old data
    !
    call dealloc_morph_in(morph_in)
    call init_morph_in(morph_in)

    ! check line count
    !
    call check_morph(file, num_morph_bb, num_morph_sc)

    ! allocate
    !
    call alloc_morph_in(morph_in, num_morph_bb, num_morph_sc)
    morph_in%num_morph_bb = num_morph_bb
    morph_in%num_morph_sc = num_morph_sc

    num_morph_bb  = 0
    num_morph_sc  = 0

    !
    ! read
    !
    backbone   = .true.

    ! check format of Morphinb data
    !
    do while(.true.)

      read(file, fmt='(a80)', end=100, err=910) line

      line = adjustl(line)

      idxb = index(line,'[')
      idxe = index(line,']')

      if (idxb > 0 .and. idxe > 0) then
        read(line(idxb+1:idxe-1),*) dir_str
        if (dir_str .eq. 'morph_bb') then
          backbone = .true.
        else
          backbone = .false.
        end if

      else if (line(1:1) .ne. ';') then

        num_char=split_num(line)
        if (num_char == 3) then
          if (backbone) then

            num_morph_bb = num_morph_bb + 1
            read(line,*) &
               morph_in%morph_bb(num_morph_bb)%atom_idx1, &
               morph_in%morph_bb(num_morph_bb)%atom_idx2, &
               morph_in%morph_bb(num_morph_bb)%rmin

          else

            num_morph_sc = num_morph_sc + 1
            read(line,*) &
               morph_in%morph_sc(num_morph_sc)%atom_idx1, &
               morph_in%morph_sc(num_morph_sc)%atom_idx2, &
               morph_in%morph_sc(num_morph_sc)%rmin
          end if
        else if (num_char == 4) then
          if (backbone) then

            num_morph_bb = num_morph_bb + 1
            read(line,*) &
               morph_in%morph_bb(num_morph_bb)%atom_idx1, &
               morph_in%morph_bb(num_morph_bb)%atom_idx2, &
               morph_in%morph_bb(num_morph_bb)%rmin,      &
               morph_in%morph_bb(num_morph_bb)%rmin_other

          else

            num_morph_sc = num_morph_sc + 1
            read(line,*) &
               morph_in%morph_sc(num_morph_sc)%atom_idx1, &
               morph_in%morph_sc(num_morph_sc)%atom_idx2, &
               morph_in%morph_sc(num_morph_sc)%rmin,      &
               morph_in%morph_bb(num_morph_bb)%rmin_other
          end if

        endif

      end if

    end do
100 return

910 call error_msg('Read_Morph> read error (read line)')

    return

  end subroutine read_morph

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_morph
  !> @brief        write section [morph_bb]
  !! @authors      CK
  !! @param[in]    file   : file unit number
  !! @param[in]    morph  : information morph pairs
  !! @param[in]    itype  : input type
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_morph(file, morph, itype)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_morph_pair),      intent(in)    :: morph(:)
    integer,                 intent(in)    :: itype

    ! local variables
    integer                  :: i


    if (itype == MorphBB) then
      write(file,'(a)') '[ morph_bb ]'
    else
      write(file,'(a)') '[ morph_sc ]'
    endif

    write(file,'(a)') '; i  j  length '

    do i = 1, size(morph)

      write(file,'(i8,1x,i8,1x,2e16.9)') &
           morph(i)%atom_idx1, &
           morph(i)%atom_idx2, &
           morph(i)%rmin,      &
           morph(i)%rmin_other

    end do

    write(file,'(a)') ''

    return

  end subroutine write_morph

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_morph
  !> @brief        check format of Morphing file
  !! @authors      CK
  !! @param[in]    file         : unit number of local restraint file
  !! @param[out]   num_morph_bb : number of backbone morph
  !! @param[out]   num_morph_sc : number of sidechain morph
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_morph(file, num_morph_bb, num_morph_sc)

    ! formal arugments
    integer,                 intent(in)    :: file
    integer,                 intent(out)   :: num_morph_bb
    integer,                 intent(out)   :: num_morph_sc
 
    ! local variables
    character(80)            :: line
    integer                  :: num_char
    logical                  :: backbone
    integer                  :: idxb, idxe
    character(20)            :: dir_str


    num_morph_bb  = 0
    num_morph_sc  = 0
    backbone   = .true.

    ! check format of Morphinb data
    !
    do while(.true.)

      read(file, fmt='(a80)', end=100, err=910) line

      line = adjustl(line)

      idxb = index(line,'[')
      idxe = index(line,']')

      if (idxb > 0 .and. idxe > 0) then
        read(line(idxb+1:idxe-1),*,err=100,end=100) dir_str
        if (dir_str .eq. 'morph_bb') then
          backbone = .true.
        else
          backbone = .false.
        end if

      else if (line(1:1) .ne. ';') then

        num_char=split_num(line)
        if (num_char == 3) then
          if (backbone) then
            num_morph_bb = num_morph_bb + 1
          else
            num_morph_sc = num_morph_sc + 1
          end if
        else if (num_char == 4) then
          if (backbone) then
            num_morph_bb = num_morph_bb + 1
          else
            num_morph_sc = num_morph_sc + 1
          end if
        end if

      end if

    end do

100 rewind(file)

    return

910 call error_msg('Check_Morph> read error (read line)')

  end subroutine check_morph

end module fileio_morph_mod
