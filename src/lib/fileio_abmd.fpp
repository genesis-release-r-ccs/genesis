!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_abmd_mod
!> @brief   read ABMD bias potential from file
!! @authors Junichiro Koga (JK)
!
!  (c) Copyright 2021 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fileio_abmd_mod

  use fileio_mod
  use bspline_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_abmd
    integer               :: nreac_coord
    integer               :: nreac_coord_comb
    integer,  allocatable :: nxi(:)
    real(wp), allocatable :: dxi(:)
    real(dp), allocatable :: ubuf(:)

  end type s_abmd

  public  :: input_abmd
  public  :: dealloc_abmd
  public  :: read_bias_from_file

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_abmd
  !> @brief        a driver subroutine for reading ABMD bias potential file
  !! @authors      JK
  !! @param[in]    mode_filename : filename of ABMD file
  !! @param[out]   mode          : structure of ABMD information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_abmd(abmd_filename, abmd)
    ! formal arguments
    character(*),            intent(in)    :: abmd_filename
    type(s_abmd),            intent(inout) :: abmd
   
    ! local variables
    integer                  :: file
    

    ! open ABMD bias potential file
    !
    call open_file(file, abmd_filename, IOFileInput)

    ! read bias potential from ABMD bias potential file
    !
    call read_abmd_bias_potential(file, abmd)

    ! close ABMD bias potential file
    !
    call close_file(file)

    return

  end subroutine input_abmd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_bias_from_file
  !> @brief        judge whether the bias potential has been read from file or
  !                not
  !! @authors      JK
  !! @param[inout] abmd     : structure of ABMD information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  logical function read_bias_from_file(abmd)
    type(s_abmd), intent(in) :: abmd
    read_bias_from_file = allocated(abmd%nxi) .and. &
                          allocated(abmd%dxi) .and. &
                          allocated(abmd%ubuf)
  end function read_bias_from_file

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_abmd
  !> @brief        deallocate ABMD information
  !! @authors      JK
  !! @param[inout] abmd     : structure of ABMD information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_abmd(abmd)
    type(s_abmd),            intent(inout) :: abmd

    ! local variables
    integer                  :: dealloc_stat

    ! initialize
    !
    dealloc_stat = 0

    ! deallocate selected variable
    !  
    if (allocated(abmd%nxi)) then
      deallocate (abmd%nxi,  &
                  abmd%dxi,  &
                  abmd%ubuf, &
                  stat = dealloc_stat)
    end if

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_abmd

  subroutine read_abmd_bias_potential(file, abmd)
    integer,      intent(in)    :: file
    type(s_abmd), intent(inout) :: abmd
    character(80)               :: strbuf
    integer                     :: intbuf
    real(kind=wp)               :: rbuf
    real(kind=wp), allocatable  :: rbufa(:)
    integer,       allocatable  :: indices(:)
    integer                     :: i, j, icount

    rewind(file)
    read(file,*,err=914) strbuf,strbuf,abmd%nreac_coord
    allocate(abmd%nxi(abmd%nreac_coord))
    allocate(abmd%dxi(abmd%nreac_coord))
    do i = 1, abmd%nreac_coord
      read(file,*,err=914) &
      strbuf, strbuf, strbuf, strbuf, strbuf, strbuf, strbuf, &
      intbuf, intbuf, abmd%nxi(i), rbuf, rbuf, abmd%dxi(i)
    enddo

    abmd%nreac_coord_comb = bspline_get_total_num_elements(&
                            abmd%nreac_coord, abmd%nxi)
    allocate(abmd%ubuf(abmd%nreac_coord_comb))
    allocate(rbufa(abmd%nreac_coord))
    allocate(indices(abmd%nreac_coord))

    abmd%ubuf = 0.0_wp

    do while(.true.) 
      read(file, *, end=1000, err=915)
      icount = 0
      combo:do i = 1, abmd%nreac_coord_comb
        call bspline_resolve_indices(i, abmd%nreac_coord, abmd%nxi, indices)
        do j = 1, abmd%nreac_coord
          if(indices(j)==-1 .or. indices(j) == abmd%nxi(j)) then
            cycle combo
          endif
        enddo
        read(file,*,end=1000, err=915) (rbufa(j),j=1, abmd%nreac_coord), &
                                       abmd%ubuf(i)
        icount = icount + 1
        if (icount == abmd%nxi(1)) then
          read(file, *, end=1000, err=915)
          icount = 0
        endif
      enddo combo
    enddo
1000 rewind(file)
    deallocate(rbufa)
    deallocate(indices)

    return

914 call error_msg('Read_Bias> read error (header)')
915 call error_msg('Read_Bias> read error (bias potential)')

  end subroutine read_abmd_bias_potential

end module fileio_abmd_mod

