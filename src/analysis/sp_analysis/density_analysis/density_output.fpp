!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   density_output_str_mod
!> @brief   subroutines for density output
!! @authors Daisuke Matsuoka (DM)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module density_output_mod

  use density_option_mod
  use density_option_str_mod
  use output_str_mod
  use molecules_mod
  use molecules_str_mod
  use fileio_pdb_mod
  use fileio_mod
  use string_mod
  use messages_mod
  use constants_mod
  
  implicit none
  private

  ! structure
  type, public :: s_density
    integer(kind=4)       :: ngrid_x
    integer(kind=4)       :: nxmin
    integer(kind=4)       :: nxmax
    integer(kind=4)       :: ngrid_y
    integer(kind=4)       :: nymin
    integer(kind=4)       :: nymax
    integer(kind=4)       :: ngrid_z
    integer(kind=4)       :: nzmin
    integer(kind=4)       :: nzmax
    real(dp), allocatable :: value(:,:,:)
    real(wp)              :: voxel_size_x
    real(wp)              :: voxel_size_y
    real(wp)              :: voxel_size_z
  end type s_density
  
  ! subroutines
  public  :: write_density
  private :: write_xplor
  private :: write_ccp4
  private :: write_dx
  private :: write_pdb
  public  :: alloc_density
  public  :: dealloc_density

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_density
  !> @brief        Write density to an output file
  !! @authors      DM
  !! @param[in]    output   : output information
  !! @param[in]    option   : option information
  !! @param[in]    density  : density information
  !! @param[in]    molecule : molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_density(output, option, density, molecule)

    ! formal arguments
    type(s_output),          intent(in)    :: output
    type(s_density_option),  intent(in)    :: option
    type(s_density),         intent(in)    :: density
    type(s_molecule),        intent(in)    :: molecule


    select case (option%output_format)

    case (DensityFormatXPLOR)
      call write_xplor(output, density, option)

    case (DensityFormatCCP4)
      call write_ccp4(output, density)

    case (DensityFormatDX)
      call write_dx(output, density, option)

    case default
      call error_msg('Write_Density> Unknown grid format type.')

    end select

    call write_pdb(output, molecule);

    return

  end subroutine write_density

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_xplor
  !> @brief        Write the density data in xplor format
  !! @authors      DM
  !! @param[in]    output     : output information
  !! @param[in]    density    : density information
  !! @param[in]    option     : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_xplor(output, density, option)

    ! formal arguments
    type(s_output),          intent(in)    :: output
    type(s_density),         intent(in)    :: density
    type(s_density_option),  intent(in)    :: option

    ! local varialbles
    integer                  :: ix, iy, iz, unit_no


    ! open density data file
    !
    call open_file(unit_no, output%mapfile, IOFileOutputNew)

    ! write density data
    !
    write(unit_no,'(/i8,a8)') 2,' !NTITLE'
    write(unit_no,'(A)')  ' REMARKS: This file was created with density_analysis.'
    if (option%density_type== DensityTypeElectron) then
      write(unit_no,'(A)')" REMARKS: density_option is 'ELECTRON'"
    else
      write(unit_no,'(A)')" REMARKS: density_option is 'NUMBER'"
    end if

    write(unit_no, '(9i8)') density%ngrid_x, density%nxmin, density%nxmax, &
                            density%ngrid_y, density%nymin, density%nymax, &
                            density%ngrid_z, density%nzmin, density%nzmax
                             
    write(unit_no, '(6e12.5)') real(density%ngrid_x, wp)*density%voxel_size_x, &
                               real(density%ngrid_y, wp)*density%voxel_size_y, &
                               real(density%ngrid_z, wp)*density%voxel_size_z, &
                               90.0_wp, 90.0_wp, 90.0_wp

    write(unit_no, '(A3)') 'ZYX'
    
    do iz = density%nzmin, density%nzmax
      write(unit_no, '(i8)') iz
      write(unit_no, '(6e12.5)') ((density%value(ix, iy, iz),  &
                                   ix = density%nxmin, density%nxmax), &
                                   iy = density%nymin, density%nymax)
    end do
    
    write(unit_no, '(i8)') -9999

    ! close the file
    !
    call close_file(unit_no)
    
    return

  end subroutine write_xplor

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_ccp4
  !> @brief        Write the density data in ccp4 binary format
  !! @authors      DM
  !! @param[in]    output     : output information
  !! @param[in]    density    : density information
  !! @note : http://www.ccp4.ac.uk/html/maplib.html
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_ccp4(output, density)

    ! formal arguments
    type(s_output),          intent(in) :: output
    type(s_density),         intent(in) :: density

    ! local varialbles
    real(sp)                 :: xlength, ylength, zlength
    real(sp)                 :: mindensity, maxdensity
    real(sp)                 :: avedensity, sigma
    real(sp)                 :: skew_mat(3,3), skew_trs(3)
    integer                  :: ix, iy, iz, idx, unit_no
    integer                  :: ndata
    character(MaxFilename)   :: filename

    ! constants
    integer (kind=4), parameter :: lmode = 2
    character(len=4), parameter :: mapstring = 'MAP '
    integer (kind=4), parameter :: ispg = 1
    integer (kind=4), parameter :: nlabel = 0
    real(sp),         parameter :: dummy = 0.0_sp
    integer (kind=4), parameter :: skew_flg = 0
    character(len=1), parameter :: dummychar = char(0)  ! NULL
    integer (kind=4), parameter :: mapc = 1
    integer (kind=4), parameter :: mapr = 2
    integer (kind=4), parameter :: maps = 3
    integer (kind=1)            :: machine_stamp(4)

    ! specify endianness used in an output binary file
    ! In this program, little endian is selected.
    data machine_stamp  /z'44', z'41', z'00', z'00'/

    ndata = density%ngrid_x * density%ngrid_y * density%ngrid_z

    skew_mat(1,:) = (/1.0_sp, 0.0_sp, 0.0_sp/)
    skew_mat(2,:) = (/0.0_sp, 1.0_sp, 0.0_sp/)
    skew_mat(3,:) = (/0.0_sp, 0.0_sp, 1.0_sp/)
    skew_trs(:)   = (/0.0_sp, 0.0_sp, 0.0_sp/)

    xlength = real(density%ngrid_x, sp) * real(density%voxel_size_x, sp)
    ylength = real(density%ngrid_y, sp) * real(density%voxel_size_y, sp)
    zlength = real(density%ngrid_z, sp) * real(density%voxel_size_z, sp)

    mindensity = real(minval(density%value), sp)
    maxdensity = real(maxval(density%value), sp)
    avedensity = real(sum(density%value), sp) / real(ndata, sp)

    sigma = 0.0_sp
    do iz = density%nzmin, density%nzmax
      do iy = density%nymin, density%nymax
        do ix = density%nxmin, density%nxmax
          sigma = sigma + (real(density%value(ix, iy, iz), sp) - avedensity) * &
                          (real(density%value(ix, iy, iz), sp) - avedensity)
        end do
      end do
    end do
    sigma = sqrt(sigma / real(ndata, sp))

    ! open density data file
    !
    unit_no  = get_unit_no()
    filename = output%mapfile
    open (unit    = unit_no,         &
          file    = filename,        &
          status  = 'new',           &
          access  = 'stream',        &
          form    = 'unformatted',   &
          convert = 'little_endian', &
          err     = 900)

    write(unit_no) density%ngrid_x   ! NC: number of Columns
    write(unit_no) density%ngrid_y   ! NR: number of Rows
    write(unit_no) density%ngrid_z   ! NS: number of Sections
    write(unit_no) lmode             ! MODE: Data type (stored as real)

    write(unit_no) density%nxmin     ! NCSTART: number of first COLUMN
    write(unit_no) density%nymin     ! NRSTART: number of first ROW
    write(unit_no) density%nzmin     ! NSSTART: number of first SECTION
    write(unit_no) density%ngrid_x   ! NX: numberf of intervals along X
    write(unit_no) density%ngrid_y   ! NY: numberf of intervals along y
    write(unit_no) density%ngrid_z   ! NZ: numberf of intervals along z
    write(unit_no) xlength           ! X_length (Cell dimensions)
    write(unit_no) ylength           ! Y_length (Cell dimensions)
    write(unit_no) zlength           ! Z_length (Cell dimensions)
    write(unit_no) 90.0_sp           ! Alpha    (Cell angles)
    write(unit_no) 90.0_sp           ! Beta     (Cell angles)
    write(unit_no) 90.0_sp           ! Gamma    (Cell angles)
    write(unit_no) mapc              ! MAPC: X axis corresponds to Columns
    write(unit_no) mapr              ! MAPR: Y axis corresponds to Rows
    write(unit_no) maps              ! MAPS: Z axis corresponds to Sections

    write(unit_no) mindensity        ! AMIN:  minimum density value
    write(unit_no) maxdensity        ! AMAX:  maximum density value
    write(unit_no) avedensity        ! AMEAN: average density value
    write(unit_no) ispg              ! ISPG:  space group number (1 -> P1)

    write(unit_no) dummy             ! NSYMBT: bytes used for storing symmetry operators
    write(unit_no) skew_flg          ! LSKGLH: skew matrix flag: 0:NONE
    write(unit_no) skew_mat(1,1)     ! SKWMAT: Skew matrix S11
    write(unit_no) skew_mat(1,2)     ! SKWMAT: Skew matrix S12
    write(unit_no) skew_mat(1,3)     ! SKWMAT: Skew matrix S13
    write(unit_no) skew_mat(2,1)     ! SKWMAT: Skew matrix S21
    write(unit_no) skew_mat(2,2)     ! SKWMAT: Skew matrix S22
    write(unit_no) skew_mat(2,3)     ! SKWMAT: Skew matrix S23
    write(unit_no) skew_mat(3,1)     ! SKWMAT: Skew matrix S31
    write(unit_no) skew_mat(3,2)     ! SKWMAT: Skew matrix S32
    write(unit_no) skew_mat(3,3)     ! SKWMAT: Skew matrix S33
    write(unit_no) skew_trs(1)       ! SKWTRN: Skew translation t1
    write(unit_no) skew_trs(2)       ! SKWTRN: Skew translation t2
    write(unit_no) skew_trs(3)       ! SKWTRN: Skew translation t3

    do ix = 38, 52
      write(unit_no) dummy
    end do

    write(unit_no) mapstring        ! MAP :   character string to identify file type
    write(unit_no) machine_stamp    ! MACHST: machine type which wrote file
    write(unit_no) sigma            ! ARMS:   RMS deviation of map from mean density
    write(unit_no) nlabel           ! NLABL:  number of labels being used

    do ix = 1, 800
      write(unit_no) dummychar
    end do

    ! density value
    do iz = density%nzmin, density%nzmax
      do iy = density%nymin, density%nymax
        do ix = density%nxmin, density%nxmax
          write(unit_no) real(density%value(ix, iy, iz), sp)
        end do
      end do
    end do

    ! close file
    !
    call close_file(unit_no)

    return

900 write(ErrOut,'(A)') 'Open_Binary_File> ', trim(filename)
    call error_msg('Open_Binary_file> Error in opening ' // trim(filename) // &
                   ' (Output file already exists, or other reasons)')

  end subroutine write_ccp4

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_dx
  !> @brief        Write the density in OpenDX grid data format
  !! @authors      DM
  !! @param[in]    output  : output information
  !! @param[in]    density : density information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_dx(output, density, option)

    ! formal arguments
    type(s_output),          intent(in)    :: output
    type(s_density),         intent(in)    :: density
    type(s_density_option),  intent(in)    :: option

    ! local varialbles
    integer                  :: ix, iy, iz, idx, unit_no
    integer                  :: idata, ndata


    ! open density data file
    !
    call open_file(unit_no, output%mapfile, IOFileOutputNew)

    ! write density data
    !
    ndata = density%ngrid_x * density%ngrid_y * density%ngrid_z

    write(unit_no,'(A)')   '# This file is created by density_analysis'
    if (option%density_type== DensityTypeElectron) then
      write(unit_no,'(A)') "# density_option is 'ELECTRON'"
    else
      write(unit_no,'(A)') "# density_option is 'NUMBER'"
    end if
    write(unit_no,'(A)')   '# '

    write(unit_no,'("object 1 class gridpositions counts ",i0,1x,i0,1x,i0)') &
         density%ngrid_x, density%ngrid_y, density%ngrid_z

    write(unit_no,'("origin",3(1x,f10.3))') &
         real(density%nxmin) * density%voxel_size_x, &
         real(density%nymin) * density%voxel_size_y, &
         real(density%nzmin) * density%voxel_size_z

    write(unit_no,'("delta ",3(1x,f10.3))') density%voxel_size_x, 0.0_wp, 0.0_wp
    write(unit_no,'("delta ",3(1x,f10.3))') 0.0_wp, density%voxel_size_y, 0.0_wp
    write(unit_no,'("delta ",3(1x,f10.3))') 0.0_wp, 0.0_wp, density%voxel_size_z

    write(unit_no,'("object 2 class gridconnections counts ",i0,1x,i0,1x,i0)') &
         density%ngrid_x, density%ngrid_y, density%ngrid_z

    write(unit_no,'("object 3 class array type double rank 0 items ",i0," data follows")') &
          ndata

    idata = 0
    do ix = density%nxmin, density%nxmax
      do iy = density%nymin, density%nymax
        do iz = density%nzmin, density%nzmax
          idata = idata + 1
          write(unit_no,'(e14.5,$)') density%value(ix, iy, iz)
          if (mod(idata, 3) == 0 .and. idata /= ndata) then
            write(unit_no,'(A)') ''
          end if
        end do
      end do
    end do
    write(unit_no,'(A)') ''

    write(unit_no,'(A)') 'attribute "dep" string "positions"'
    write(unit_no,'(A)') 'object "regular positions regular connections" class field'
    write(unit_no,'(A)') 'component "positions" value 1'
    write(unit_no,'(A)') 'component "connections" value 2'
    write(unit_no,'(A)') 'component "data" value 3'

    ! close the file
    !
    call close_file(unit_no)

    return

  end subroutine write_dx

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_pdb
  !> @brief        Write the reference molecule structure
  !! @authors      NT
  !! @param[in]    output   : output information
  !! @param[in]    molecule : molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_pdb(output, molecule)

    ! formal arguments
    type(s_output),          intent(in) :: output
    type(s_molecule),        intent(in) :: molecule

    ! local variables
    type(s_pdb)              :: pdb


    if (output%pdbfile .eq. '') &
      return

    call export_molecules(molecule, pdb = pdb)
    call output_pdb(output%pdbfile, pdb)

    call dealloc_pdb_all(pdb)

    return

  end subroutine write_pdb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_density
  !> @brief        allocate density
  !! @authors      DM
  !! @param[inout] density  : density informatoin
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_density(density)

    ! formal arguments
    type(s_density),         intent(inout) :: density

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat


    alloc_stat   = 0
    dealloc_stat = 0

    if (allocated(density%value)) &
      deallocate(density%value, stat = dealloc_stat)

    if (dealloc_stat /= 0) &
      call error_msg_dealloc

    if (.not. allocated(density%value)) then
      allocate(density%value(density%ngrid_x,  &
                             density%ngrid_y,  &
                             density%ngrid_z), &
               stat = alloc_stat)

      density%value(:,:,:) = 0.0_dp
    end if

    if (alloc_stat /= 0) &
      call error_msg_alloc

    return 

  end subroutine alloc_density

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_densityr
  !> @brief        deallocate density
  !! @authors      DM
  !! @param[inout] density  : density informatoin
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_density(density)

    ! formal arguments
    type(s_density),         intent(inout) :: density

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0

      if (allocated(density%value)) then
        deallocate(density%value, stat = dealloc_stat)
      end if

    if (dealloc_stat /= 0) call error_msg_dealloc

    return 

  end subroutine dealloc_density

end module density_output_mod
