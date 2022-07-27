!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   eg_generate_mod
!> @brief   generate trajectory files
!! @authors Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module eg_generate_mod

  use eg_option_mod
  use eg_option_str_mod
  use output_str_mod
  use molecules_str_mod
  use molecules_mod
  use fileio_mod
  use messages_mod
  use constants_mod
  use fileio_pdb_mod

  implicit none
  private

  ! subroutines
  public :: generate

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    generate
  !> @brief        generate trajectory files
  !! @authors      TM
  !! @param[inout] molecule   : molecule information
  !! @param[inout] trj_list   : trajectory list information
  !! @param[inout] trajectory : trajectory information
  !! @param[inout] fitting    : fitting information
  !! @param[inout] option     : option information
  !! @param[inout] output     : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

 subroutine generate(molecule, option, output)

    ! formal arguments
    type(s_molecule), target, intent(inout) :: molecule
    type(s_option),           intent(inout) :: option
    type(s_output),           intent(inout) :: output

    ! local variables
    integer                  :: i, j, k, n, itmp1, itmp2
    integer                  :: natom
    integer                  :: n_grid_cut_x, n_grid_cut_y, n_grid_cut_z
    integer                  :: nx, ny, nz, nb(3), max_nxyz
    integer                  :: map_out, pdb_out
    integer, allocatable     :: ig(:,:), ig_upper(:,:), ig_lower(:,:)
    real(wp)                 :: f, d, pi, r, y, vol, total_mass
    real(wp)                 :: mincrd(3), maxcrd(3)
    real(wp)                 :: coeff_rho, zfactor, yzfactor
    real(wp), allocatable    :: simulated_map(:,:,:)
    real(wp), allocatable    :: bound_x(:), bound_y(:), bound_z(:)
    real(wp), allocatable    :: derfa_x(:,:), derfa_y(:,:), derfa_z(:,:), erfa(:,:)
    real(wp), pointer        :: coord(:,:)
    type(s_pdb)              :: pdb,ref

    natom  =   molecule%num_atoms
    coord  =>  molecule%atom_coord


    ! check-only
    if (option%check_only) &
      return


    ! check molecule size
    !
    mincrd(1) = minval(coord(1,1:natom))
    mincrd(2) = minval(coord(2,1:natom))
    mincrd(3) = minval(coord(3,1:natom))
    maxcrd(1) = maxval(coord(1,1:natom))
    maxcrd(2) = maxval(coord(2,1:natom))
    maxcrd(3) = maxval(coord(3,1:natom))
    write(MsgOut,'(A)') 'Generate> Check molecule size'
    write(MsgOut,'(A,f10.3,A,f10.3)') '  min x           = ', mincrd(1), &
                                      '  max x           = ', maxcrd(1)
    write(MsgOut,'(A,f10.3,A,f10.3)') '  min y           = ', mincrd(2), &
                                      '  max y           = ', maxcrd(2)
    write(MsgOut,'(A,f10.3,A,f10.3)') '  min z           = ', mincrd(3), &
                                      '  max z           = ', maxcrd(3)
    write(MsgOut,'(A)') ''


    ! determine_cutoff
    !
    pi  = acos(-1.0_wp)
    r   = 0.0_wp
    y   = 0.0_wp

    do while (1.0_wp - y > option%tolerance)
      y = erf(sqrt(3.0_wp/2.0_wp)*r) - sqrt(6.0_wp/pi)*r*exp(-3.0_wp*r*r/2.0_wp)
      r = r + 0.01_wp
    end do

    n_grid_cut_x = ceiling(r*option%sigma/option%voxel_size)
    n_grid_cut_y = ceiling(r*option%sigma/option%voxel_size)
    n_grid_cut_z = ceiling(r*option%sigma/option%voxel_size)


    ! determine map size
    !
    if (option%auto_margin) then
      option%x0 = mincrd(1) - option%margin_size_x
      option%y0 = mincrd(2) - option%margin_size_y
      option%z0 = mincrd(3) - option%margin_size_z
      option%box_size_x = maxcrd(1) - mincrd(1) + 2.0_wp*option%margin_size_x
      option%box_size_y = maxcrd(2) - mincrd(2) + 2.0_wp*option%margin_size_y
      option%box_size_z = maxcrd(3) - mincrd(3) + 2.0_wp*option%margin_size_z
      write(MsgOut,'(A)') 'Generate> Auto_margin was done'
      write(MsgOut,'(A)') ''
    end if

    nx = int(option%box_size_x/option%voxel_size + 0.5_wp)
    ny = int(option%box_size_y/option%voxel_size + 0.5_wp)
    nz = int(option%box_size_z/option%voxel_size + 0.5_wp)

    allocate(simulated_map(0:nx,0:ny,0:nz))
    allocate(bound_x(0:nx))
    allocate(bound_y(0:ny))
    allocate(bound_z(0:nz))

    simulated_map(:,:,:) = 0.0_wp

    do i = 0, nx
      bound_x(i) = option%x0 + option%voxel_size * (dble(i) - 0.5_wp)
    end do
    do i = 0, ny
      bound_y(i) = option%y0 + option%voxel_size * (dble(i) - 0.5_wp)
    end do
    do i = 0, nz
      bound_z(i) = option%z0 + option%voxel_size * (dble(i) - 0.5_wp)
    end do


    ! calculate the index of the grid that each atom is located
    !
    allocate(ig(3,natom),ig_upper(3,natom),ig_lower(3,natom))

    do n = 1, natom
      ig(1,n) = int((coord(1,n) - option%x0)/option%voxel_size + 0.5_wp)
      ig(2,n) = int((coord(2,n) - option%y0)/option%voxel_size + 0.5_wp)
      ig(3,n) = int((coord(3,n) - option%z0)/option%voxel_size + 0.5_wp)

      if (ig(1,n) < n_grid_cut_x .or. ig(1,n) >= nx - n_grid_cut_x .or. &
          ig(2,n) < n_grid_cut_y .or. ig(2,n) >= ny - n_grid_cut_y .or. &
          ig(3,n) < n_grid_cut_z .or. ig(3,n) >= nz - n_grid_cut_z ) then
        write(MsgOut,'(A,I10)') 'Gaussian kernes is extending outside the map box for atom ', n
        write(MsgOut,'(A,3F12.3)') 'Coordinate: ', coord(1,n), coord(2,n), coord(3,n)
        write(MsgOut,*) 'grid index igx(n), igy(n), igz(n): ', ig(1,n), ig(2,n), ig(3,n)
        write(MsgOut,*) 'grid index range, xmin(n), xmax(n), ymin(n), ymax(n), zmin(n), zmax(n) ', &
          int(ig(1,n),8)-n_grid_cut_x, int(ig(1,n),8)+n_grid_cut_x, int(ig(2,n),8)-n_grid_cut_y,   &
          int(ig(2,n),8)+n_grid_cut_y, int(ig(3,n),8)-n_grid_cut_z, int(ig(3,n),8)+n_grid_cut_z
        call error_msg('Generate> Gaussian kernes is extending outside the map box')
      end if

      ig_lower(1,n) = ig(1,n) - n_grid_cut_x
      ig_lower(2,n) = ig(2,n) - n_grid_cut_y
      ig_lower(3,n) = ig(3,n) - n_grid_cut_z
      ig_upper(1,n) = ig(1,n) + n_grid_cut_x
      ig_upper(2,n) = ig(2,n) + n_grid_cut_y
      ig_upper(3,n) = ig(3,n) + n_grid_cut_z
    end do


    ! calculate error and exp functions and differences for integral
    !
    f = sqrt(3.0_wp / (2.0_wp * option%sigma**2))

    nb(1) = nx
    nb(2) = ny
    nb(3) = nz
    max_nxyz = maxval(nb)

    allocate(erfa   (0:max_nxyz,1:natom))
    allocate(derfa_x(0:nx-1,1:natom))
    allocate(derfa_y(0:ny-1,1:natom))
    allocate(derfa_z(0:nz-1,1:natom))
    erfa(:,:) = 0.0_wp
    derfa_x(:,:) = 0.0_wp
    derfa_y(:,:) = 0.0_wp
    derfa_z(:,:) = 0.0_wp

    do n = 1, natom
      itmp1 = ig(1,n) - n_grid_cut_x
      itmp2 = ig(1,n) + n_grid_cut_x + 1
      do i = itmp1, itmp2
        d = bound_x(i) - coord(1,n)
        erfa(i,n) = erf(f*d)
      end do
      do i = itmp1, itmp2 - 1
        derfa_x(i,n) = erfa(i+1,n) - erfa(i,n)
      end do

      itmp1 = ig(2,n) - n_grid_cut_y
      itmp2 = ig(2,n) + n_grid_cut_y + 1
      do j = itmp1, itmp2
        d = bound_y(j) - coord(2,n)
        erfa(j,n) = erf(f*d)
      end do
      do j = itmp1, itmp2 - 1
        derfa_y(j,n) = erfa(j+1,n) - erfa(j,n)
      end do

      itmp1 = ig(3,n) - n_grid_cut_z
      itmp2 = ig(3,n) + n_grid_cut_z + 1
      do k = itmp1, itmp2
        d = bound_z(k) - coord(3,n)
        erfa(k,n) = erf(f*d)
      end do
      do k = itmp1, itmp2 - 1
        derfa_z(k,n) = erfa(k+1,n) - erfa(k,n)
      end do
    end do


    ! calculate simulated density map
    !
    vol = option%voxel_size**3
    coeff_rho = (option%sigma**2 * pi / 6.0_wp)**(3.0_wp/2.0_wp) / vol

    do n = 1, natom
      do k = ig_lower(3,n), ig_upper(3,n)
        zfactor = coeff_rho*derfa_z(k,n)
        do j = ig_lower(2,n), ig_upper(2,n)
          yzfactor = derfa_y(j,n)*zfactor
          do i = ig_lower(1,n), ig_upper(1,n)
            simulated_map(i,j,k) = simulated_map(i,j,k) + derfa_x(i,n)*yzfactor
          end do
        end do
      end do
    end do


    ! open output files
    !
    if (output%mapfile /= '') &
      call open_file(map_out, output%mapfile, IOFileOutputNew)

      write(map_out,'(4f12.6,3i12)') option%voxel_size, &
                 option%x0, option%y0, option%z0, nx, ny, nz

      write(map_out,'(a)') ''

      write(map_out,'(10f12.6)') simulated_map(0:nx-1,0:ny-1,0:nz-1)

    if (output%mapfile /= '') call close_file(map_out)

    return

  end subroutine generate

end module eg_generate_mod
