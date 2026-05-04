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

module fileio_spline_mod

  use fileio_mod
  use bspline_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_spline
    integer               :: nreac_coord
    integer               :: nreac_coord_comb
    integer,  allocatable :: nxi(:)
    real(wp), allocatable :: dxi(:)
    real(dp), allocatable :: ubuf(:)

  end type s_spline

  public  :: input_spline
  public  :: dealloc_spline
  public  :: read_bias_from_file

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_spline
  !> @brief        a driver subroutine for reading ABMD bias potential file
  !! @authors      JK
  !! @param[in]    mode_filename : filename of ABMD file
  !! @param[out]   mode          : structure of ABMD information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_spline(spline_filename, spline)
    ! formal arguments
    character(*),            intent(in)    :: spline_filename
    type(s_spline),            intent(inout) :: spline
   
    ! local variables
    integer                  :: file
    

    ! open ABMD bias potential file
    !
    call open_file(file, spline_filename, IOFileInput)

    ! read bias potential from ABMD bias potential file
    !
    call read_spline_bias_potential(file, spline)

    ! close ABMD bias potential file
    !
    call close_file(file)

    return

  end subroutine input_spline

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_bias_from_file
  !> @brief        judge whether the bias potential has been read from file or
  !                not
  !! @authors      JK
  !! @param[inout] spline     : structure of ABMD information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  logical function read_bias_from_file(spline)
    type(s_spline), intent(in) :: spline
    read_bias_from_file = allocated(spline%nxi) .and. &
                          allocated(spline%dxi) .and. &
                          allocated(spline%ubuf)
  end function read_bias_from_file

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_spline
  !> @brief        deallocate ABMD information
  !! @authors      JK
  !! @param[inout] spline     : structure of ABMD information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_spline(spline)
    type(s_spline),            intent(inout) :: spline

    ! local variables
    integer                  :: dealloc_stat

    ! initialize
    !
    dealloc_stat = 0

    ! deallocate selected variable
    !  
    if (allocated(spline%nxi)) then
      deallocate (spline%nxi,  &
                  spline%dxi,  &
                  spline%ubuf, &
                  stat = dealloc_stat)
    end if

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_spline

  subroutine read_spline_bias_potential(file, spline)
    integer,      intent(in)    :: file
    type(s_spline), intent(inout) :: spline
    character(80)               :: strbuf
    integer                     :: intbuf
    real(kind=wp)               :: rbuf
    real(kind=wp), allocatable  :: rbufa(:)
    integer,       allocatable  :: indices(:)
    integer                     :: i, j, icount

    ! Check nreac_coord from bspline file
    read(file, *, err=914)
    read(file, *, err=914) strbuf, strbuf, strbuf, spline%nreac_coord

    allocate(spline%nxi(spline%nreac_coord))
    allocate(spline%dxi(spline%nreac_coord))

    rewind(file)
    read(file, *, err=914)
    read(file, *, err=914) strbuf, strbuf, strbuf, intbuf, (spline%nxi(i), i=1, spline%nreac_coord), (spline%dxi(i), i=1, spline%nreac_coord)

    
    !read(file,*,err=914) strbuf,strbuf,intbuf
    !!allocate(spline%nxi(spline%nreac_coord))
    !!allocate(spline%dxi(spline%nreac_coord))
    !do i = 1, intbuf !spline%nreac_coord
    !  !read(file,*,err=914) &
    !  !strbuf, strbuf, strbuf, strbuf, strbuf, strbuf, strbuf, &
    !  !intbuf, intbuf, spline%nxi(i), rbuf, rbuf, spline%dxi(i)
    !  read(file,*,err=914) strbuf
    !enddo

    spline%nreac_coord_comb = bspline_get_total_num_elements(&
                            spline%nreac_coord, spline%nxi)
    allocate(spline%ubuf(spline%nreac_coord_comb))
    allocate(rbufa(spline%nreac_coord))
    allocate(indices(spline%nreac_coord))

    spline%ubuf = 0.0_wp

    do while(.true.) 
      !read(file, *, end=1000, err=915)
      icount = 0
      combo:do i = 1, spline%nreac_coord_comb
        call bspline_resolve_indices(i, spline%nreac_coord, spline%nxi, indices)
        do j = 1, spline%nreac_coord
          if(indices(j)==-1 .or. indices(j) == spline%nxi(j)) then
            cycle combo
          endif
        enddo
        read(file,*,end=1000, err=915) (rbufa(j),j=1, spline%nreac_coord), &
                                       spline%ubuf(i)
        icount = icount + 1
        if (icount == spline%nxi(1)) then
          read(file, *, end=1000, err=915)
          icount = 0
        endif
      enddo combo
      !read(file, *, end=1000, err=915)
    enddo
1000 rewind(file)
    deallocate(rbufa)
    deallocate(indices)

    return

914 call error_msg('Read_Bias> read error (header)')
915 call error_msg('Read_Bias> read error (bias potential)')

  end subroutine read_spline_bias_potential

end module fileio_spline_mod

