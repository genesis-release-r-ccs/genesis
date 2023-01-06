!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_pme_opt_1dalltoall_mod
!> @brief   Smooth particle mesh ewald method
!! @authors Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_pme_opt_1dalltoall_mod

  use sp_boundary_str_mod
  use sp_enefunc_str_mod
  use sp_energy_str_mod
  use sp_domain_str_mod
  use math_libs_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use fft3d_1dalltoall_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  real(wp),         save :: el_fact      ! e^2/(4*pi*eps) A*kcal/mol
  real(wp),         save :: alpha        ! alpha
  real(wp),         save :: inv_alpha    ! 1/alpha
  real(wp),         save :: alpha2m      ! -alpha^2
  real(wp),         save :: alpha3       ! alpha^3
  real(wp),         save :: alpha6       ! alpha^6
  real(wp),         save :: vol_fact2    ! (2*pi/V)*el_fact
  real(wp),         save :: vol_fact4    ! (4*pi/V)*el_fact
  real(wp),         save :: vol_fact2_lj ! (sqrt(pi^3)/24V)
  real(wp),         save :: vol_fact4_lj ! (2*sqrt(pi^3)/24V)
  real(dp), public, save :: u_self_o1    ! Ewald self energy
  real(wp),         save :: inv_volume   ! inverse of volume
  real(wp),         save :: box(3)       ! box size
  real(wp),         save :: box_inv(3)   ! Inverse of box size
  real(wp),         save :: bs_fact      ! B-spline factor (1/(n-1)...2)
  real(wp),         save :: bs_fact3     ! bs_fact^3
  real(wp),         save :: bs_fact3d    ! bs_fact^3*(n-1)
  real(wp),         save :: r_scale(3)   ! coordinate-scaling factor (I/L)
  integer,          save :: n_bs         ! Order of B-spline
  integer,          save :: grid_bd      ! boundary grid number
  integer,          save :: ngrid(4)     ! Number of grid
  integer,          save :: nx           ! process number in x dimension
  integer,          save :: ny           ! process number in y dimension
  integer,          save :: nz           ! process number in z dimension
  integer,          save :: x_start
  integer,          save :: x_end
  integer,          save :: x_start1
  integer,          save :: x_end1
  integer,          save :: x_local1 
  integer,          save :: y_start
  integer,          save :: y_end
  integer,          save :: z_start
  integer,          save :: z_end
  integer,          save :: nlocalx
  integer,          save :: nlocalx1
  integer,          save :: nlocaly
  integer,          save :: nlocalz
  integer,          save :: maxproc
  integer,          save :: ngridmax

  integer,          save, allocatable :: grid_g2lx(:)
  integer,          save, allocatable :: grid_g2ly(:)
  integer,          save, allocatable :: grid_g2lz(:)
  real(wp),         save, allocatable :: b2(:,:)        ! b^2(hx,hy,hz)
  real(wp),         save, allocatable :: gx(:)
  real(wp),         save, allocatable :: gy(:)
  real(wp),         save, allocatable :: gz(:)
  real(wp),         save, allocatable :: vir_fact(:,:,:)! -2*(1+G^2/4a)/G^2
  real(wp),         save, allocatable :: vir_fact_lj(:,:,:)
  real(wp),         save, allocatable :: theta(:,:,:)   ! prefactor in recip energy (elec)
  real(wp),         save, allocatable :: theta_lj(:,:,:)   ! prefactor in recip energy (dispersion)
  real(wp),         save, allocatable :: f(:,:,:)
  integer,          save, allocatable :: vi(:,:,:)
  integer,          save, allocatable :: vii(:,:,:,:)
  real(wp),         save, allocatable :: bsc(:,:,:,:)
  real(wp),         save, allocatable :: bscd(:,:,:,:)
  real(wp),         save, allocatable :: qdf(:,:,:,:)
  real(wp),         save, allocatable :: ldf(:,:,:,:)
  real(wp),         save, allocatable :: buf_send(:,:), buf_recv(:,:)
  real(wp),         save, allocatable :: qdf_real(:)
  real(wp),         save, allocatable :: qdf_work(:,:)
  real(wp),         save, allocatable :: ldf_real(:)
  real(wp),         save, allocatable :: ldf_work(:,:)
  complex(wp),      save, allocatable :: ftqdf(:)
  complex(wp),      save, allocatable :: ftqdf_work(:)
  complex(wp),      save, allocatable :: ftldf(:)
  complex(wp),      save, allocatable :: ftldf_work(:)

  complex(wp),      save, allocatable :: ftqdf_work2(:,:,:)
  complex(wp),      save, allocatable :: ftldf_work2(:,:,:)

  ! parameters
  integer, parameter      :: NumIndex        = 45
  integer, parameter      :: Index(NumIndex) = (/ &
                              1,   2,   3,   4,   5,   6,   8,   9,  10,  12, &
                             15,  16,  18,  20,  24,  25,  27,  30,  32,  36, &
                             40,  45,  48,  50,  54,  60,  64,  72,  75,  80, &
                             81,  90,  96, 100, 120, 125, 128, 135, 144, 150, &
                            160, 162, 180, 192, 200/)

  ! subroutines
  public  :: setup_pme_opt_1dalltoall
  public  :: dealloc_pme_opt_1dalltoall
  public  :: pme_pre_opt_1dalltoall
  public  :: pme_pre_opt_1dalltoall_lj
  public  :: pme_recip_opt_1dalltoall_4
  public  :: pme_recip_opt_1dalltoall_lj_4
  public  :: pme_recip_opt_1dalltoall_6
  public  :: pme_recip_opt_1dalltoall_lj_6
  public  :: pme_recip_opt_1dalltoall_8
  public  :: pme_recip_opt_1dalltoall_lj_8

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_pme_opt_1dalltoall
  !> @brief        Setup for PME with domain decomposition
  !! @authors      JJ, NT
  !! @param[in]    domain   : domain information
  !! @param[in]    boundary : boundary information
  !! @param[in]    enefunc  : energy function
  !! @param[inout] calc     : flag of whether calculated
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_pme_opt_1dalltoall(domain, boundary, enefunc, calc)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_boundary),        intent(in)    :: boundary
    type(s_enefunc),         intent(in)    :: enefunc
    logical,                 intent(inout) :: calc

    ! local variables
    real(wp)                 :: bcos, bsin, fact, grid_space
    integer                  :: i, j, k, js
    integer                  :: localmax, ncell
    integer                  :: iproc(3), index_x, index_y, index_z
    integer                  :: remainder, quotient, expo
    integer                  :: alloc_stat

    real(wp),    allocatable :: bs(:)


    ncell  =  domain%num_cell_local + domain%num_cell_boundary

    ! Setting the number of grids (even number)
    !
    ngrid(1) = enefunc%pme_ngrid_x
    ngrid(2) = enefunc%pme_ngrid_y
    ngrid(3) = enefunc%pme_ngrid_z
    ngridmax = max(ngrid(1),ngrid(2),ngrid(3))

    if (ngrid(1) == 0 .and. ngrid(2) == 0 .and. ngrid(3) == 0) then

      grid_space = enefunc%pme_max_spacing
      ngrid(1) = int(boundary%box_size_x/grid_space)
      ngrid(2) = int(boundary%box_size_y/grid_space)
      ngrid(3) = int(boundary%box_size_z/grid_space)
      ngridmax = max(ngrid(1),ngrid(2),ngrid(3))

    end if

    ! Setting common parameters
    !
    el_fact    = ELECOEF/enefunc%dielec_const
    n_bs       = enefunc%pme_nspline
    alpha      = enefunc%pme_alpha
    alpha2m    = -alpha**2
    alpha3     = alpha*alpha*alpha
    alpha6     = alpha3*alpha3
    inv_alpha  = 1.0_wp/alpha

    ! grid partition
    !
    if (reciprocal_calc) then

      nx = boundary%num_domain(1)
      ny = boundary%num_domain(2)
      nz = boundary%num_domain(3)

      ! Check if grid numbers in each subdomain is sufficent for communication
      !
      ngridmax = ngrid(1) / nx
      ngridmax = min(ngridmax, ngrid(2)/ny)
      ngridmax = min(ngridmax, ngrid(3)/nz)
      if (ngridmax < 2*n_bs) &
        call error_msg('Cell size or PME grids should be larger')

      ! Check grid point in x direction
      !
      remainder = mod(ngrid(1),2*nx)
      if (remainder /= 0) ngrid(1) = ngrid(1) + 2*nx - remainder

      quotient = ngrid(1)/(2*nx)
      if (quotient <= Index(NumIndex)) then
        do i = 1, NumIndex
          if (quotient <= Index(i)) exit
        end do
        quotient = Index(i)
        ngrid(1) = (2*nx) * quotient
      else
        expo = int(log(real(quotient,wp))/log(real(2,wp)))
        if (2**expo >= quotient) then
          quotient = 2**expo
        else
          quotient = 2**(expo+1)
        end if
        ngrid(1) = (2*nx) * quotient
      end if

      ! check grid points in y direction
      !
      if (mod(nz,2) == 0) then

        remainder = mod(ngrid(2),ny*nz)
        if (remainder /= 0) ngrid(2) = ngrid(2) + nz*ny - remainder

        quotient = ngrid(2)/(nz*ny)
        if (quotient <= Index(NumIndex)) then
          do i = 1, NumIndex
            if (quotient <= Index(i)) exit
          end do
          quotient = Index(i)
          ngrid(2) = (nz*ny) * quotient
        else
          expo = int(log(real(quotient,wp))/log(real(2,wp)))
          if (2**expo >= quotient) then
            quotient = 2**expo
          else
            quotient = 2**(expo+1)
          end if
          ngrid(2) = (nz*ny) * quotient
        end if

      else

        remainder = mod(ngrid(2),ny*nz*2)
        if (remainder /= 0) ngrid(2) = ngrid(2) + nz*ny*2 - remainder

        quotient = ngrid(2)/(nz*ny*2)
        if (quotient <= Index(NumIndex)) then
          do i = 1, NumIndex
            if (quotient <= Index(i)) exit
          end do
          quotient = Index(i)
          ngrid(2) = (nz*ny*2) * quotient
        else
          expo = int(log(real(quotient,wp))/log(real(2,wp)))
          if (2**expo >= quotient) then
            quotient = 2**expo
          else
            quotient = 2**(expo+1)
          end if
          ngrid(2) = (nz*ny*2) * quotient
        end if

      end if

      ! Check grid point in z direction
      !
      k = lcm(nz*nx, nz*ny)
      remainder = mod(ngrid(3),k)
      if (remainder /= 0) ngrid(3) = ngrid(3) + k - remainder

      quotient = ngrid(3)/(k)
      if (quotient <= Index(NumIndex)) then
        do i = 1, NumIndex
          if (quotient <= Index(i)) exit
        end do
        quotient = Index(i)
        ngrid(3) = k * quotient
      else
        expo = int(log(real(quotient,wp))/log(real(2,wp)))
        if (2**expo >= quotient) then
          quotient = 2**expo
        else
          quotient = 2**(expo+1)
        end if
        ngrid(3) = k * quotient
      end if

      if ((enefunc%pme_ngrid_x /= ngrid(1)) .or. &
          (enefunc%pme_ngrid_y /= ngrid(2)) .or. &
          (enefunc%pme_ngrid_z /= ngrid(3))) then

        if (enefunc%pme_ngrid_x /= 0 .or. enefunc%pme_ngrid_y /= 0 .or. &
            enefunc%pme_ngrid_z /= 0) calc = .false.

        if (main_rank) then
          write(MsgOut,'(A)') &
            ''
          if (enefunc%pme_ngrid_x==0 .and. enefunc%pme_ngrid_y==0 .and. &
              enefunc%pme_ngrid_z==0) then
            write(MsgOut,'(A)') &
              'Setup_PME> Proper PME grid number was generated automatically'
          else
            write(MsgOut,'(A)') &
              '  WARNING: PME grid number is different from the input'
          end if
          write(MsgOut,'(A20,3I10)') &
            '  pme_ngrid(x,y,z)= ', ngrid(1), ngrid(2), ngrid(3)
          write(MsgOut,'(A)') &
            ''
        end if
      end if

      ngridmax = max(ngrid(1),ngrid(2),ngrid(3))

      ! index for each processor
      !

      iproc(1) = mod(my_city_rank, boundary%num_domain(1))
      iproc(2) = mod(my_city_rank/boundary%num_domain(1),boundary%num_domain(2))
      iproc(3) = my_city_rank/(boundary%num_domain(1)*boundary%num_domain(2))

      nlocalx  = int((ngrid(1)+nx-1)/nx)
      nlocaly  = int((ngrid(2)+ny-1)/ny)
      nlocalz  = int((ngrid(3)+nz-1)/nz)

      x_start  = (iproc(1))*nlocalx + 1
      y_start  = (iproc(2))*nlocaly + 1
      z_start  = (iproc(3))*nlocalz + 1

      x_end    = (iproc(1)+1)*nlocalx 
      y_end    = (iproc(2)+1)*nlocaly 
      z_end    = (iproc(3)+1)*nlocalz 

      localmax = max(nlocalx*nlocaly,nlocalx*nlocalz,nlocaly*nlocalz)
      index_x  = iproc(2)*nlocaly*nlocalz + iproc(3)
      index_y  = iproc(3)*nlocalx*nlocalz + iproc(1)
      index_z  = iproc(1)*nlocalx*nlocaly + iproc(2)

      maxproc = max(nprocx,nprocy,nprocz,1)

      nlocalx1 = nlocalx/2 + 1
      if (my_x_rank == 0) then
        x_local1 = nlocalx1
        x_start1 = 1
        x_end1   = x_local1
      else
        x_local1 = nlocalx1 - 1
        x_start1 = my_x_rank*nlocalx/2 + 1
        x_end1   = x_start1 + x_local1 - 1
      end if
    end if

    if (reciprocal_calc) then

      j = 1
      do i = 1, n_bs - 2
        j = j*(n_bs - i)
      end do

      bs_fact   = 1.0_wp/real(j,wp)
      bs_fact3  = bs_fact**3
      bs_fact3d = bs_fact3 * real(n_bs - 1,wp)

      ! Preparing b2=b(h)^2, h shifted
      !
      allocate(bs(n_bs), b2(ngridmax,3))

      call b_spline_coef(n_bs, 1.0_wp, bs)

      do i = 1, n_bs - 1
        bs(i) = bs(i)*bs_fact
      end do

      do k = 1, 3

        fact = 2.0_wp * PI/real(ngrid(k),wp)

        do j = 1, ngrid(k)/2 + 1

          js = j - 1
          bcos = 0.0_wp
          bsin = 0.0_wp

          do i = 0, n_bs - 2
            bcos = bcos + bs(i+1) * cos(real(js*i,wp)*fact)
            bsin = bsin + bs(i+1) * sin(real(js*i,wp)*fact)
          end do

          b2(j, k) = 1.0_wp/(bcos*bcos + bsin*bsin)

        end do

        do j = ngrid(k)/2 + 2, ngrid(k)

          js = j - 1 - ngrid(k)
          bcos = 0.0_wp
          bsin = 0.0_wp

          do i = 0, n_bs - 2
            bcos = bcos + bs(i+1) * cos(real(js*i,wp)*fact)
            bsin = bsin + bs(i+1) * sin(real(js*i,wp)*fact)
          end do

          b2(j, k) = 1.0_wp/(bcos*bcos + bsin*bsin)

        end do

      end do

      if (mod(n_bs,2) == 1) then
        do k = 1, 3
          b2(ngrid(k)/2 + 1,k) = 0.0_wp
        end do
      end if

      ! setup global/local grid indices
      !
      allocate(grid_g2lx(ngrid(1)), grid_g2ly(ngrid(2)),      &
               grid_g2lz(ngrid(3)), stat = alloc_stat)

      grid_bd = n_bs
      if (nprocx >= 2) then
        if (my_x_rank == 0) then
          do k = 1, grid_bd
            grid_g2lx(ngrid(1)-grid_bd+k)  = k
            grid_g2lx(x_end+k)             = grid_bd+nlocalx+k
          end do
          grid_g2lx(ngrid(1)-n_bs)         = 1
          grid_g2lx(x_end+n_bs+1)          = nlocalx+2*n_bs
        else if (my_x_rank == (nprocx-1)) then
          do k = 1, grid_bd
            grid_g2lx(x_start-1-grid_bd+k) = k
            grid_g2lx(k)                   = grid_bd+nlocalx+k
          end do
          grid_g2lx(x_start-n_bs-1)        = 1
          grid_g2lx(n_bs+1)                = nlocalx+2*n_bs
        else
          do k = 1, grid_bd
            grid_g2lx(x_start-1-grid_bd+k) = k
            grid_g2lx(x_end+k)             = grid_bd+nlocalx+k
          end do
          grid_g2lx(x_start-n_bs-1)        = 1
          grid_g2lx(x_end+n_bs+1)          = nlocalx+2*n_bs
        end if
      else
        do k = 1, grid_bd
          grid_g2lx(ngrid(1)-grid_bd+k)    = k
          grid_g2lx(k)                     = grid_bd+nlocalx+k
        end do
        grid_g2lx(ngrid(1)-n_bs)           = 1
        grid_g2lx(n_bs+1)                  = nlocalx+2*n_bs
      end if
      do k = x_start, x_end
        grid_g2lx(k)                       = grid_bd + k - x_start + 1
      end do
      if (nprocy >= 2) then
        if (my_y_rank == 0) then
          do k = 1, grid_bd
            grid_g2ly(ngrid(2)-grid_bd+k)  = k
            grid_g2ly(y_end+k)             = grid_bd+nlocaly+k
          end do
          grid_g2ly(ngrid(2)-n_bs)         = 1
          grid_g2ly(y_end+n_bs+1)          = nlocaly+2*n_bs
        else if (my_y_rank == (nprocy-1)) then
          do k = 1, grid_bd
            grid_g2ly(y_start-1-grid_bd+k) = k
            grid_g2ly(k)                   = grid_bd+nlocaly+k
          end do
          grid_g2ly(y_start-n_bs-1)        = 1
          grid_g2ly(n_bs+1)                = nlocaly+2*n_bs
        else
          do k = 1, grid_bd
            grid_g2ly(y_start-1-grid_bd+k) = k
            grid_g2ly(y_end+k)             = grid_bd+nlocaly+k
          end do
          grid_g2ly(y_start-n_bs-1)        = 1
          grid_g2ly(y_end+n_bs+1)          = nlocaly+2*n_bs
        end if
      else
        do k = 1, grid_bd
          grid_g2ly(ngrid(2)-grid_bd+k)    = k
          grid_g2ly(k)                     = grid_bd+nlocaly+k
        end do
        grid_g2ly(ngrid(2)-n_bs)           = 1
        grid_g2ly(n_bs+1)                  = nlocaly+2*n_bs
      end if
      do k = y_start, y_end
        grid_g2ly(k)                       = grid_bd + k - y_start + 1
      end do

      if (nprocz >= 2) then
        if (my_z_rank == 0) then
          do k = 1, grid_bd
            grid_g2lz(ngrid(3)-grid_bd+k)  = k
            grid_g2lz(z_end+k)             = grid_bd+nlocalz+k
          end do
          grid_g2lz(ngrid(3)-n_bs)         = 1
          grid_g2lz(z_end+n_bs+1)          = nlocalz+2*n_bs
        else if (my_z_rank == (nprocz-1)) then
          do k = 1, grid_bd
            grid_g2lz(z_start-1-grid_bd+k) = k
            grid_g2lz(k)                   = grid_bd+nlocalz+k
          end do
          grid_g2lz(z_start-n_bs-1)        = 1
          grid_g2lz(n_bs+1)                = nlocalz+2*n_bs
        else
          do k = 1, grid_bd
            grid_g2lz(z_start-1-grid_bd+k) = k
            grid_g2lz(z_end+k)             = grid_bd+nlocalz+k
          end do
          grid_g2lz(z_start-n_bs-1)        = 1
          grid_g2lz(z_end+n_bs+1)          = nlocalz+2*n_bs
        end if
      else
        do k = 1, grid_bd
          grid_g2lz(ngrid(3)-grid_bd+k)    = k
          grid_g2lz(k)                     = grid_bd+nlocalz+k
        end do
        grid_g2lz(ngrid(3)-n_bs)           = 1
        grid_g2lz(n_bs+1)                  = nlocalz+2*n_bs
      end if
      do k = z_start, z_end
        grid_g2lz(k) = grid_bd + k - z_start + 1
      end do

      ! Prepareing theta = F^-1[theta](h), h shifted
      !
      localmax = max((nlocalx+2*grid_bd)*(nlocaly+2*grid_bd), &
                     (nlocalx+2*grid_bd)*(nlocalz+2*grid_bd), &
                     (nlocaly+2*grid_bd)*(nlocalz+2*grid_bd))
      localmax = localmax*grid_bd

      ! Prepareing theta = F^-1[theta](h), h shifted
      !
      if (enefunc%vdw == VDWPME) then
        allocate(qdf(nlocalx+2*n_bs,nlocaly+2*n_bs,nlocalz+2*n_bs,nthread), &
                 ldf(nlocalx+2*n_bs,nlocaly+2*n_bs,nlocalz+2*n_bs,nthread), &
                 qdf_real(nlocalx*nlocaly*nlocalz),                         &
                 ldf_real(nlocalx*nlocaly*nlocalz),                         &
                 buf_send(2*localmax,2), buf_recv(2*localmax,2),            &
                 gx(x_local1), gy(nlocaly), gz(ngrid(3)),                   &
                 vir_fact(ngrid(3), nlocaly, x_local1),                     &
                 vir_fact_lj(ngrid(3), nlocaly, x_local1),                  &
                 theta(ngrid(3), nlocaly, x_local1),                        &
                 theta_lj(ngrid(3), nlocaly, x_local1),                     &
                 qdf_work(nlocalx*nlocaly*nlocalz, maxproc),                &
                 ldf_work(nlocalx*nlocaly*nlocalz, maxproc),                &
                 ftqdf(nlocalx1*nlocaly*nlocalz),                           &
                 ftqdf_work(nlocalx1*nlocaly*nlocalz*maxproc),              &
                 ftqdf_work2(ngrid(3),nlocalx1,nlocaly),                    &
                 ftldf(nlocalx1*nlocaly*nlocalz),                           &
                 ftldf_work(nlocalx1*nlocaly*nlocalz*maxproc),              &
                 ftldf_work2(ngrid(3),nlocalx1,nlocaly),                    &
                 stat = alloc_stat)
      else
        allocate(qdf(nlocalx+2*n_bs,nlocaly+2*n_bs,nlocalz+2*n_bs,nthread), &
                 qdf_real(nlocalx*nlocaly*nlocalz),                         &
                 buf_send(localmax,2), buf_recv(localmax,2),                &
                 gx(x_local1), gy(nlocaly), gz(ngrid(3)),                   &
                 vir_fact(ngrid(3), nlocaly, x_local1),                     &
                 theta(ngrid(3), nlocaly, x_local1),                        &
                 qdf_work(nlocalx*nlocaly*nlocalz, maxproc),                &
                 ftqdf(nlocalx1*nlocaly*nlocalz),                           &
                 ftqdf_work(nlocalx1*nlocaly*nlocalz*maxproc),              &
                 ftqdf_work2(ngrid(3),nlocalx1,nlocaly),                    &
                 stat = alloc_stat) 
      end if
      if (alloc_stat /= 0) call error_msg_alloc

      allocate(f(3,MaxAtom,ncell),                                   &
               bsc(n_bs,3,MaxAtom,ncell),bscd(n_bs,3,MaxAtom,ncell), &
               vi(4,MaxAtom,ncell), vii(n_bs,3,MaxAtom,ncell),       &
               stat = alloc_stat)
      if (alloc_stat /= 0)   call error_msg_alloc

      deallocate(bs)

    end if

    return

  end subroutine setup_pme_opt_1dalltoall

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_pme_opt_1dalltoall
  !> @brief        deallocate all arrays used in PME
  !! @authors      JJ, NT
  !! @param[in]    enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_pme_opt_1dalltoall(enefunc)

    type(s_enefunc),          intent(in) :: enefunc

    ! local variables
    integer :: dealloc_stat, ierr


    if (.not. allocated(f)) &
      return 

    deallocate(f, bsc, bscd, vi, vii, stat = dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc

    if (enefunc%vdw == VDWPME) then
      deallocate(qdf, qdf_real, ldf, ldf_real, buf_send, buf_recv, gx, gy, gz, &
                 vir_fact, vir_fact_lj, theta, theta_lj, qdf_work, ldf_work,   &
                 ftqdf, ftqdf_work, ftqdf_work2, ftldf, ftldf_work,            &
                 ftldf_work2,b2, stat = dealloc_stat)
    else
      deallocate(qdf, qdf_real, buf_send, buf_recv, gx, gy, gz, vir_fact, &
                 theta, qdf_work, ftqdf, ftqdf_work, ftqdf_work2, b2,     &
                 stat = dealloc_stat)
    end if
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_pme_opt_1dalltoall

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_pre_opt_1dalltoall
  !> @brief        Prepare functions for PME calculation
  !! @authors      JJ, NT
  !! @param[in]    domain   : domain information
  !! @param[in]    boundary : boundary information
  !! @note         Extracted from setup_pme for NPT calculation
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_pre_opt_1dalltoall(domain, boundary)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_boundary),        intent(in)    :: boundary

    ! local variables
    integer                  :: i, j, k, is, js, ks, ix, iy
    integer                  :: niy
    integer                  :: iy_start, iy_end
    real(wp)                 :: fact, gfact(3), g2


    box(1) = boundary%box_size_x
    box(2) = boundary%box_size_y
    box(3) = boundary%box_size_z

    do k = 1, 3
      box_inv(k) = 1.0_wp / box(k)
      r_scale(k) = real(ngrid(k),wp) * box_inv(k)
    end do

    vol_fact2 = 2.0_wp * PI * el_fact * box_inv(1)*box_inv(2)*box_inv(3)
    vol_fact4 = 2.0_wp * vol_fact2

    ! Prepareing theta = F^-1[theta](h), h shifted
    ! Gx, Gy, Gz, and vir_fact=-2(1+G^2/4a)/G^2 are also prepared 
    ! for virial calculation
    !
    do k = 1, 3
      gfact(k) = 2.0_wp * PI * box_inv(k)
    end do

    fact = 0.25_wp / alpha2m

    niy = nlocaly / nprocz
    iy_start = y_start + niy*my_z_rank
    iy_end   = iy_start + niy - 1

    do i = x_start1, x_end1

      ix = i - x_start1 + 1
      is = i - 1
      gx(ix) = gfact(1) * real(is,wp)

      do j = iy_start, iy_end

        iy = j - iy_start + 1
        if (j <= ngrid(2)/2+1) then
          js = j - 1
        else
          js = j - 1 - ngrid(2)
        end if
        gy(iy) = gfact(2) * real(js,wp)

        do k = 1, ngrid(3)
          if (k <= ngrid(3)/2+1) then
            ks = k - 1
          else
            ks = k - 1 - ngrid(3)
          end if
          gz(k) = gfact(3) * real(ks,wp)
          g2 = gx(ix)*gx(ix) + gy(iy)*gy(iy) + gz(k)*gz(k)
          if (g2 > EPS) then
            vir_fact(k,iy,ix) = -2.0_wp * (1.0_wp - g2 * fact)/g2
            if (g2*fact < -80.0_wp) then
              theta(k,iy,ix) = 0.0_wp
            else
              theta(k,iy,ix) = b2(i,1)*b2(j,2)*b2(k,3)*exp(g2 * fact)/g2
            end if
            if (abs(theta(k,iy,ix)) < 1.0e-15) theta(k,iy,ix) = 0.0_wp
          else
            vir_fact(k,iy,ix) = 0.0_wp
            theta(k,iy,ix) = 0.0_wp
          end if
        end do

      end do

    end do

    if (my_x_rank == 0) then

      i  = ngrid(1)/2 + 1
      ix = nlocalx/2 + 1
      is = i - 1
      gx(ix) = gfact(1) * real(is,wp)

      do j = iy_start, iy_end
        iy = j - iy_start + 1
        if (j <= ngrid(2)/2+1) then
          js = j - 1
        else
          js = j - 1 - ngrid(2)
        end if
        gy(iy) = gfact(2) * real(js,wp)
        do k = 1, ngrid(3)
          if (k <= ngrid(3)/2+1) then
            ks = k - 1
          else
            ks = k - 1 - ngrid(3)
          end if
          gz(k) = gfact(3) * real(ks,wp)
          g2 = gx(ix)**2 + gy(iy)**2 + gz(k)*gz(k)
          if (g2 > EPS) then
            vir_fact(k,iy,ix) = -2.0_wp * (1.0_wp - g2 * fact)/g2
            if (g2*fact < -80.0_wp) then
              theta(k,iy,ix) = 0.0_wp
            else
              theta(k,iy,ix) = b2(i,1) * b2(j,2) * b2(k,3) * exp(g2 * fact)/g2
            end if
            if (abs(theta(k,iy,ix)) < 1.0e-15) theta(k,iy,ix) = 0.0_wp
          else
            vir_fact(k,iy,ix) = 0.0_wp
            theta(k,iy,ix) = 0.0_wp
          end if
        end do
      end do

    end if

    if (my_x_rank == 0 .and. my_y_rank == 0 .and. my_z_rank == 0) &
      theta(1,1,1) = 0.0_wp

    if (my_city_rank == 0) &
      vir_fact(1,1,1) = 0.0_wp

    ! Calculating self energy
    !
    u_self_o1 = 0.0_dp 

    do i = 1, domain%num_cell_local
      do ix = 1, domain%num_atom(i)
        u_self_o1 = u_self_o1 + domain%charge(ix,i)*domain%charge(ix,i)
      end do 
    end do

    u_self_o1 = - u_self_o1 * el_fact * alpha/sqrt(PI)

    return

  end subroutine pme_pre_opt_1dalltoall

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_pre_opt_1dalltoall_lj
  !> @brief        Prepare functions for PME calculation (Elec and LJ)
  !! @authors      JJ, NT
  !! @param[in]    domain   : domain information
  !! @param[in]    boundary : boundary information
  !! @note         Extracted from setup_pme for NPT calculation
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_pre_opt_1dalltoall_lj(domain, boundary)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_boundary),        intent(in)    :: boundary

    ! local variables
    integer                  :: i, j, k, is, js, ks, ix, iy
    integer                  :: niy
    integer                  :: iy_start, iy_end
    real(wp)                 :: fact, gfact(3), g1, g2, g3, g_lj
    real(wp)                 :: fact1, fact2, fact3, fact4, fact5


    box(1) = boundary%box_size_x
    box(2) = boundary%box_size_y
    box(3) = boundary%box_size_z

    do k = 1, 3
      box_inv(k) = 1.0_wp / box(k)
      r_scale(k) = real(ngrid(k),wp) * box_inv(k)
    end do

    inv_volume   = box_inv(1)*box_inv(2)*box_inv(3)
    vol_fact2    = 2.0_wp * PI * el_fact * inv_volume
    vol_fact4    = 2.0_wp * vol_fact2
    vol_fact2_lj = PI * sqrt(PI) * inv_volume / (24.0_wp)
    vol_fact4_lj = 2.0_wp * vol_fact2_lj

    ! factor that converts from m to G=2*pi*m/L
    !
    gfact(1:3) = 2.0_wp * PI * box_inv(1:3)

    fact = 0.25_wp / alpha2m

    niy = nlocaly / nprocz
    iy_start = y_start + niy*my_z_rank
    iy_end   = iy_start + niy - 1

    do i = x_start1, x_end1

      ix = i - x_start1 + 1
      is = i - 1
      gx(ix) = gfact(1) * real(is,wp)

      do j = iy_start, iy_end

        iy = j - iy_start + 1
        if (j <= ngrid(2)/2+1) then
          js = j - 1
        else
          js = j - 1 - ngrid(2)
        end if
        gy(iy) = gfact(2) * real(js,wp)

        do k = 1, ngrid(3)
          if (k <= ngrid(3)/2+1) then
            ks = k - 1
          else
            ks = k - 1 - ngrid(3)
          end if
          gz(k) = gfact(3) * real(ks,wp)

          g2 = gx(ix)*gx(ix) + gy(iy)*gy(iy) + gz(k)*gz(k)
          g1 = sqrt(g2)
          g3 = g1 * g2

          ! -g^2/(4*alpha^2)
          !
          fact1 = g2*fact

          ! exp(-g^2/(4*alpha^2))
          !
          fact2 = exp(fact1)

          ! g / (2*alpha)
          !
          fact3 = 0.5*g1*inv_alpha

          ! erfc(g/(2*alpha))
          !
          fact4 = erfc(fact3)

          ! |bx(mx)|^2 * |by(my)|^2 * |bz(mz)|^2
          !
          fact5 = b2(i,1)*b2(j,2)*b2(k,3)

          ! calculate theta and vir_fact (prefactor for energy & virial)
          !
          g_lj = (4.0_wp*alpha3-2.0_wp*alpha*g2)*fact2 &
               + SQRT_PI*g3*fact4
          theta_lj(k,iy,ix) = fact5 * g_lj
          vir_fact_lj(k,iy,ix) = 3.0_wp*g1*SQRT_PI*fact4-6.0_wp*alpha*fact2
          vir_fact_lj(k,iy,ix) = vir_fact_lj(k,iy,ix) * fact5

          if (g2 > EPS) then
            vir_fact(k,iy,ix) = -2.0_wp * (1.0_wp-fact1)/g2
            theta(k,iy,ix) = fact5 * fact2 / g2
          else
            vir_fact(k,iy,ix) = 0.0_wp
            theta(k,iy,ix) = 0.0_wp
          end if
        end do

      end do

    end do

    if (my_x_rank == 0) then

      i  = ngrid(1)/2 + 1
      ix = nlocalx/2 + 1
      is = i - 1
      gx(ix) = gfact(1) * real(is,wp)

      do j = iy_start, iy_end
        iy = j - iy_start + 1
        if (j <= ngrid(2)/2+1) then
          js = j - 1
        else
          js = j - 1 - ngrid(2)
        end if
        gy(iy) = gfact(2) * real(js,wp)
        do k = 1, ngrid(3)
          if (k <= ngrid(3)/2+1) then
            ks = k - 1
          else
            ks = k - 1 - ngrid(3)
          end if
          gz(k) = gfact(3) * real(ks,wp)

          g2 = gx(ix)*gx(ix) + gy(iy)*gy(iy) + gz(k)*gz(k)
          g1 = sqrt(g2)
          g3 = g1 * g2

          ! -g^2/(4*alpha^2)
          !
          fact1 = g2*fact

          ! exp(-g^2/(4*alpha^2))
          !
          fact2 = exp(fact1)

          ! g / (2*alpha)
          !
          fact3 = 0.5*g1*inv_alpha

          ! erfc(g/(2*alpha))
          !
          fact4 = erfc(fact3)

          ! |bx(mx)|^2 * |by(my)|^2 * |bz(mz)|^2
          !
          fact5 = b2(i,1)*b2(j,2)*b2(k,3)

          ! prefactor for energy and virial
          !
          g_lj = (4.0_wp*alpha3-2.0_wp*alpha*g2)*fact2 &
               + SQRT_PI*g3*fact4
          theta_lj(k,iy,ix) = fact5 * g_lj
          vir_fact_lj(k,iy,ix) = 3.0_wp*g1*SQRT_PI*fact4-6.0_wp*alpha*fact2
          vir_fact_lj(k,iy,ix) = vir_fact_lj(k,iy,ix) * fact5

          if (g2 > EPS) then
            vir_fact(k,iy,ix) = -2.0_wp * (1.0_wp-fact1)/g2
            theta(k,iy,ix) = fact5 * fact2/g2
          else
            vir_fact(k,iy,ix) = 0.0_wp
            theta(k,iy,ix) = 0.0_wp
          end if
        end do
      end do

    end if

    if (my_city_rank == 0) then
      theta(1,1,1) = 0.0_wp
      theta_lj(1,1,1) = 0.0_wp
      vir_fact(1,1,1) = 0.0_wp
      vir_fact_lj(1,1,1) = 0.0_wp
    end if

    ! Calculating self energy
    !
    u_self_o1 = 0.0_dp 

    do i = 1, domain%num_cell_local
      do ix = 1, domain%num_atom(i)
        u_self_o1 = u_self_o1 + domain%charge(ix,i)*domain%charge(ix,i)
      end do 
    end do

    u_self_o1 = - u_self_o1 * el_fact * alpha/sqrt(PI)

    return

  end subroutine pme_pre_opt_1dalltoall_lj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_recip_opt_1dalltoall_4
  !> @brief        Calculate PME reciprocal part with domain decomposition
  !                (B spline order 4)                                     
  !! @authors      JJ, NT
  !! @param[in]    domain : domain information
  !! @param[inout] force  : forces of target systems
  !! @param[inout] virial : virial term of target systems
  !! @param[inout] eelec  : electrostatic energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_recip_opt_1dalltoall_4(domain, force, virial, eelec)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    real(wip),               intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: elec_temp
    real(wp)                 :: vr(4), dv(4), grid(4), half_grid(4)
    real(wp)                 :: bscx(4), bscy(4), bscz(4)
    real(wp)                 :: bscdx(4), bscdy(4), bscdz(4)
    real(wp)                 :: vxx, vyy, vzz, vyyd, vzzd
    real(wp)                 :: tqq, qtmp, tqq1, tqq2
    real(wp)                 :: force_local(3)
    real(wp)                 :: u_2(3), u_3(3)
    integer                  :: i, k, icel, iproc, id
    integer                  :: ix, iv, iy, iz, ixs, iys, izs
    integer                  :: ixyz, nix, nix1, niy, niz, nizx, nizy
    integer                  :: ixx, iyy, izz, iyyzz, ii(4)
    integer                  :: iix(4), iiy(4), iiz(4)
    integer                  :: start(3), end(3)
    integer                  :: ncell, ncell_local, omp_get_thread_num
    integer                  :: kk
    integer                  :: iprocx, iprocy
    integer                  :: ix_start, ix_end, iz_start, iz_end

    complex(wp), allocatable :: work1(:), work2(:)
    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    integer,         pointer :: natom(:)


    coord  => domain%coord
    charge => domain%charge
    natom  => domain%num_atom

    ncell  = domain%num_cell_local + domain%num_cell_boundary
    ncell_local  = domain%num_cell_local

    ! Initializing the energy and force
    !
    qdf_real(1:nlocalx*nlocaly*nlocalz) = 0.0_wp

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, iv, ix, iy, iz, icel, k, ii, dv, izz, izs, iyy, iys,  &
    !$omp         iz_start, iz_end, ix_start, ix_end, ixx, ixs, vxx, vyy, vzz, &
    !$omp         iproc, tqq, vr, work1, work2, elec_temp, force_local, ixyz,  &
    !$omp         tqq1, tqq2, nix, nix1, niy, niz, bscx, bscy, bscz, bscdx,    &
    !$omp         bscdy, bscdz, kk, iprocx, iprocy, qtmp, grid, half_grid,     &
    !$omp         start, end, nizx, nizy, iix, iiy, iiz, u_2, u_3, vyyd, vzzd, &
    !$omp         iyyzz)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    ! initialization
    !
    allocate(work1(ngridmax), work2(2*ngridmax))

    start(1) = x_start
    start(2) = y_start
    start(3) = z_start
    end(1) = x_end
    end(2) = y_end
    end(3) = z_end

    ! Initialization of Q-fanction
    !
    do iz = 1, nlocalz+8
      do iy = 1, nlocaly+8
        do ix = 1, nlocalx+8
          qdf(ix,iy,iz,id+1) = 0.0_wp
        end do
      end do
    end do

    ! Calculateing Q-fanction
    !
    grid(1:3) = real(ngrid(1:3),wp)
    half_grid(1:3) = real(ngrid(1:3)/2,wp)
    grid(4) = 2.0_wp
    half_grid(4) = 1.0_wp
    ngrid(4) = 100000

    do icel = id+1, ncell_local, nthread

      do i = 1, natom(icel)

        vr(1) = coord(1,i,icel) * r_scale(1)
        vr(2) = coord(2,i,icel) * r_scale(2)
        vr(3) = coord(3,i,icel) * r_scale(3)
        vr(1) = vr(1) + half_grid(1) - grid(1)*anint(vr(1)/grid(1))
        vr(2) = vr(2) + half_grid(2) - grid(2)*anint(vr(2)/grid(2))
        vr(3) = vr(3) + half_grid(3) - grid(3)*anint(vr(3)/grid(3))
        vr(1) = min(vr(1),grid(1)-0.0000001_wp)
        vr(2) = min(vr(2),grid(2)-0.0000001_wp)
        vr(3) = min(vr(3),grid(3)-0.0000001_wp)
        vi(1,i,icel) = min(int(vr(1)),ngrid(1)-1)
        vi(2,i,icel) = min(int(vr(2)),ngrid(2)-1)
        vi(3,i,icel) = min(int(vr(3)),ngrid(3)-1)
        dv(1) = vr(1) - real(vi(1,i,icel),wp)
        dv(2) = vr(2) - real(vi(2,i,icel),wp)
        dv(3) = vr(3) - real(vi(3,i,icel),wp)
 
        u_2(1) = dv(1)*dv(1)
        u_2(2) = dv(2)*dv(2)
        u_2(3) = dv(3)*dv(3)
        u_3(1) = u_2(1)*dv(1)
        u_3(2) = u_2(2)*dv(2)
        u_3(3) = u_2(3)*dv(3)
        bsc(1,1,i,icel) = u_3(1)
        bsc(1,2,i,icel) = u_3(2)
        bsc(1,3,i,icel) = u_3(3)
        bsc(2,1,i,icel) = -3.0_wp*u_3(1) + 3.0_wp*u_2(1) &
                          + 3.0_wp*dv(1) + 1.0_wp
        bsc(2,2,i,icel) = -3.0_wp*u_3(2) + 3.0_wp*u_2(2) &
                          + 3.0_wp*dv(2) + 1.0_wp
        bsc(2,3,i,icel) = -3.0_wp*u_3(3) + 3.0_wp*u_2(3) &
                          + 3.0_wp*dv(3) + 1.0_wp
        bsc(3,1,i,icel) = 3.0_wp*u_3(1) - 6.0_wp*u_2(1) + 4.0_wp
        bsc(3,2,i,icel) = 3.0_wp*u_3(2) - 6.0_wp*u_2(2) + 4.0_wp
        bsc(3,3,i,icel) = 3.0_wp*u_3(3) - 6.0_wp*u_2(3) + 4.0_wp
        bsc(4,1,i,icel) = (1.0_wp-dv(1))*(1.0_wp-dv(1))*(1.0_wp-dv(1))
        bsc(4,2,i,icel) = (1.0_wp-dv(2))*(1.0_wp-dv(2))*(1.0_wp-dv(2))
        bsc(4,3,i,icel) = (1.0_wp-dv(3))*(1.0_wp-dv(3))*(1.0_wp-dv(3))
        bscd(1,1,i,icel) = u_2(1)
        bscd(1,2,i,icel) = u_2(2)
        bscd(1,3,i,icel) = u_2(3)
        bscd(2,1,i,icel) = -3.0_wp*u_2(1) + 2.0_wp*dv(1) + 1.0_wp
        bscd(2,2,i,icel) = -3.0_wp*u_2(2) + 2.0_wp*dv(2) + 1.0_wp
        bscd(2,3,i,icel) = -3.0_wp*u_2(3) + 2.0_wp*dv(3) + 1.0_wp
        bscd(3,1,i,icel) = 3.0_wp*u_2(1) - 4.0_wp*dv(1)
        bscd(3,2,i,icel) = 3.0_wp*u_2(2) - 4.0_wp*dv(2)
        bscd(3,3,i,icel) = 3.0_wp*u_2(3) - 4.0_wp*dv(3)
        bscd(4,1,i,icel) = -(1.0_wp-dv(1))*(1.0_wp-dv(1))
        bscd(4,2,i,icel) = -(1.0_wp-dv(2))*(1.0_wp-dv(2))
        bscd(4,3,i,icel) = -(1.0_wp-dv(3))*(1.0_wp-dv(3))
  
      end do

    end do

    do icel = id+1, ncell_local, nthread
      do i = 1, natom(icel)

        qtmp = charge(i,icel)*bs_fact3
        ii(1) = vi(1,i,icel)
        ii(2) = vi(2,i,icel)
        ii(3) = vi(3,i,icel)

!ocl nosimd
        do ixyz = 1, 4
          iz = ii(3) - ixyz + 2
          if (iz <= 0) iz = iz + ngrid(3)
          vii(ixyz,3,i,icel) = grid_g2lz(iz)
          iy = ii(2) - ixyz + 2
          if (iy <= 0) iy = iy + ngrid(2)
          vii(ixyz,2,i,icel) = grid_g2ly(iy)
          ix = ii(1) - ixyz + 2
          if (ix <= 0) ix = ix + ngrid(1)
          vii(ixyz,1,i,icel) = grid_g2lx(ix)
        end do

        do izz = 1, 4
          izs = vii(izz,3,i,icel)
          vzz = bsc(izz,3,i,icel)
          do iyy = 1, 4
            iys = vii(iyy,2,i,icel)
            vyy = bsc(iyy,2,i,icel)
!ocl norecurrence(qdf)
!ocl nosimd
            do ixx = 1, 4
              ixs = vii(ixx,1,i,icel)
              qdf(ixs,iys,izs,id+1) = qdf(ixs,iys,izs,id+1)            &
                                    + bsc(ixx,1,i,icel)*vyy*vzz*qtmp
            end do
          end do
        end do
      end do
    end do

    !$omp barrier
    do iproc = 2, nthread
      do iz = id+1, nlocalz+8,nthread
        do iy = 1, nlocaly+8
          do ix = 1, nlocalx+8
            qdf(ix,iy,iz,1) = qdf(ix,iy,iz,1) + qdf(ix,iy,iz,iproc)
          end do
        end do
      end do
    end do

    !$omp barrier
    !$omp master
    call communicate_pme_pre
    !$omp end master
    !$omp barrier

    do iz = id+1, nlocalz, nthread
      izs = iz + 4
      do iy = 1, nlocaly
        k = (iy-1)*nlocalx + (iz-1)*nlocalx*nlocaly
        iys = iy + 4
        do ix = 1, nlocalx
          ixs = ix + 4
          qdf_real(k+ix) = qdf_real(k+ix) + qdf(ixs,iys,izs,1)
        end do
      end do
    end do

    !$omp barrier
    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdf_real, nlocalx*nlocaly*nlocalz/nprocx, mpi_wp_real, &
                      qdf_work, nlocalx*nlocaly*nlocalz/nprocx, mpi_wp_real, &
                      grid_commx, ierror)
#else
    qdf_work(1:nlocalx*nlocaly*nlocalz,1) = qdf_real(1:nlocalx*nlocaly*nlocalz)
#endif
    !$omp end master
    !$omp barrier

    nizx = nlocalz / nprocx
    nizy = nlocalz / nprocy
    call fft3d_1d_alltoall(qdf_work, ftqdf, ftqdf, ftqdf, ftqdf_work,        &
                           ftqdf_work, ftqdf_work, ftqdf_work, ftqdf_work2,  &
                           work1, work2, nlocalx, nlocaly, nlocalz,          &
                           nizx, nizy, nlocalx1, x_local1, ngrid(1),         &
                           ngrid(2), ngrid(3), nprocx, nprocy, nprocz, id,   &
                           nthread, grid_commx, grid_commy, grid_commz)

    ! Energy calculation
    !

    !$omp barrier

    iproc = my_z_rank + 1
    niy = nlocaly / nprocz

    do ix = 1, nlocalx/2
      do iy = id+1, niy, nthread
        elec_temp = 0.0_wp
        vxx = 0.0_wp
        vyy = 0.0_wp
        vzz = 0.0_wp
        do iz = 1, ngrid(3)
          tqq1 = real(ftqdf_work2(iz,ix,iy),wp)
          tqq2 = imag(ftqdf_work2(iz,ix,iy))
          tqq = tqq1*tqq1 + tqq2*tqq2
          tqq = tqq * theta(iz,iy,ix) * vol_fact4
          elec_temp = elec_temp + tqq

          ! virial
          !
          vxx = vxx + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gx(ix) * gx(ix))
          vyy = vyy + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gy(iy) * gy(iy))
          vzz = vzz + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gz(iz) * gz(iz))
        end do
        eelec(id+1) = eelec(id+1) + elec_temp
        virial(1,1,id+1) = virial(1,1,id+1) + vxx
        virial(2,2,id+1) = virial(2,2,id+1) + vyy
        virial(3,3,id+1) = virial(3,3,id+1) + vzz
      end do
    end do

    if (my_x_rank == 0) then
      do iy = id+1, niy, nthread
        elec_temp = 0.0_wp
        vxx = 0.0_wp
        vyy = 0.0_wp
        vzz = 0.0_wp
        do iz = 1, ngrid(3)
          tqq1 = real(ftqdf_work2(iz,1,iy),wp)
          tqq2 = imag(ftqdf_work2(iz,1,iy))
          tqq  = tqq1*tqq1 + tqq2*tqq2
          tqq = tqq * theta(iz,iy,1) * vol_fact2
          elec_temp = elec_temp - tqq
          vxx = vxx - tqq * (1.0_wp + vir_fact(iz,iy,1) * gx(1) * gx(1))
          vyy = vyy - tqq * (1.0_wp + vir_fact(iz,iy,1) * gy(iy) * gy(iy))
          vzz = vzz - tqq * (1.0_wp + vir_fact(iz,iy,1) * gz(iz) * gz(iz))
        end do
        eelec(id+1) = eelec(id+1) + elec_temp
        virial(1,1,id+1) = virial(1,1,id+1) + vxx
        virial(2,2,id+1) = virial(2,2,id+1) + vyy
        virial(3,3,id+1) = virial(3,3,id+1) + vzz
      end do

      ix = nlocalx/2 + 1
      do iy = id+1, niy, nthread
        elec_temp = 0.0_wp
        vxx = 0.0_wp
        vyy = 0.0_wp
        vzz = 0.0_wp
        do iz = 1, ngrid(3)
          tqq1 = real(ftqdf_work2(iz,ix,iy),wp)
          tqq2 = imag(ftqdf_work2(iz,ix,iy))
          tqq  = tqq1*tqq1 + tqq2*tqq2
          tqq = tqq * theta(iz,iy,ix) * vol_fact2
          elec_temp = elec_temp + tqq
          vxx = vxx + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gx(ix) * gx(ix))
          vyy = vyy + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gy(iy) * gy(iy))
          vzz = vzz + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gz(iz) * gz(iz))
        end do
        eelec(id+1) = eelec(id+1) + elec_temp
        virial(1,1,id+1) = virial(1,1,id+1) + vxx
        virial(2,2,id+1) = virial(2,2,id+1) + vyy
        virial(3,3,id+1) = virial(3,3,id+1) + vzz
      end do
    end if

    !$omp barrier

    ! F^-1[Th]*F^-1[Q] (=X)
    !
    do ix = 1, x_local1
      do iy = id+1, niy, nthread
        k = (iy-1)*ngrid(3) + (ix-1)*ngrid(3)*niy
        do iz = 1, ngrid(3)
          ftqdf_work(k+iz) = ftqdf_work2(iz,ix,iy) &
                            *cmplx(theta(iz,iy,ix),0.0_wp)
        end do
      end do
    end do

    !$omp barrier

    call bfft3d_1d_alltoall(qdf_work, qdf_real, ftqdf, ftqdf, ftqdf_work,    &
                    ftqdf_work, ftqdf_work2, ftqdf_work, work1, work2,       &
                    nlocalx, nlocaly, nlocalz, nlocalx1, x_local1, ngrid(1), &
                    ngrid(2), ngrid(3), nprocx, nprocy, nprocz, id, nthread, &
                    my_x_rank, grid_commx, grid_commy, grid_commz, niy,      &
                    nizy, nizx)

    !$omp barrier

    ! X is saved on qdf
    !
    do iz = id+1, nlocalz, nthread
      izs = iz + 4
      do iy = 1, nlocaly
        iys = iy + 4
        k = (iy-1)*nlocalx + (iz-1)*nlocalx*nlocaly
        do ix = 1, nlocalx
          ixs = ix + 4
          qdf(ixs,iys,izs,1) = qdf_real(k+ix)
        end do
      end do
    end do

    !$omp barrier

    !$omp barrier
    !$omp master
    call communicate_pme_post
    !$omp end master
    !$omp barrier

    do icel = id+1, ncell_local, nthread
      do i = 1, natom(icel)

        force_local(1) = 0.0_wp
        force_local(2) = 0.0_wp
        force_local(3) = 0.0_wp
        qtmp = charge(i,icel)*vol_fact4*bs_fact3d

!ocl nosimd
        do ixyz = 1, 4
          ii(ixyz)    = vi(ixyz,i,icel)
          bscx(ixyz)  = bsc(ixyz,1,i,icel)
          bscy(ixyz)  = bsc(ixyz,2,i,icel)
          bscz(ixyz)  = bsc(ixyz,3,i,icel)
          bscdx(ixyz) = bscd(ixyz,1,i,icel)
          bscdy(ixyz) = bscd(ixyz,2,i,icel)
          bscdz(ixyz) = bscd(ixyz,3,i,icel)
        end do

!ocl nosimd
        do ixyz = 1, 4
          ix = ii(1) - ixyz + 2
          if (ix <= 0) ix = ix + ngrid(1)
          iix(ixyz) = grid_g2lx(ix)
          iy = ii(2) - ixyz + 2
          if (iy <= 0) iy = iy + ngrid(2)
          iiy(ixyz) = grid_g2ly(iy)
          iz = ii(3) - ixyz + 2
          if (iz <= 0) iz = iz + ngrid(3)
          iiz(ixyz) = grid_g2lz(iz)
        end do

        do iyyzz = 1, 16
          izz = (iyyzz-1)/4 + 1
          iyy = iyyzz - 4*(izz-1)
          izs = iiz(izz)
          iys = iiy(iyy)
          vzz = bscz(izz)
          vzzd = bscdz(izz)
          vyy = bscy(iyy)
          vyyd = bscdy(iyy)
!ocl nosimd
          do ixx = 1, 4
            ixs = iix(ixx)
            force_local(1) = force_local(1) &
                           + bscdx(ixx)*vyy*vzz*qdf(ixs,iys,izs,1)
            force_local(2) = force_local(2)     &
                           + bscx(ixx)*vyyd*vzz*qdf(ixs,iys,izs,1)
            force_local(3) = force_local(3)     &
                           + bscx(ixx)*vyy*vzzd*qdf(ixs,iys,izs,1)
          end do
        end do

        force(1:3,i,icel) = - force_local(1:3)*qtmp*r_scale(1:3)

      end do
    end do

    deallocate(work1, work2)

    !$omp end parallel

    return

  end subroutine pme_recip_opt_1dalltoall_4

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_recip_opt_1dalltoall_lj_4
  !> @brief        Calculate PME reciprocal part with domain decomposition
  !                (B spline order 4)                                     
  !! @authors      JJ, NT
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_recip_opt_1dalltoall_lj_4(domain, enefunc, force, virial, &
                                           eelec, evdw)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wip),               intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw (nthread)

    ! local variables
    real(wp)                 :: elec_temp, evdw_temp, temp
    real(wp)                 :: vr(4), dv(4), grid(4), half_grid(4)
    real(wp)                 :: bscx(4), bscy(4), bscz(4)
    real(wp)                 :: bscdx(4), bscdy(4), bscdz(4)
    real(wp)                 :: vxx, vyy, vzz, vyyd, vzzd
    real(wp)                 :: efft_real, efft_imag
    real(wp)                 :: vfft_real, vfft_imag
    real(wp)                 :: tqq, tll, qtmp, ltmp
    real(wp)                 :: force_local(3), force_local_lj(3)
    real(wp)                 :: vol_factor
    real(wp)                 :: u_2(3), u_3(3)
    integer                  :: iatmcls
    integer                  :: i, k, icel, iproc, id
    integer                  :: ix, iv, iy, iz, ixs, iys, izs
    integer                  :: ixyz, nix, nix1, niy, niz, nizx, nizy
    integer                  :: ixx, iyy, izz, ii(4)
    integer                  :: iix(4), iiy(4), iiz(4)
    integer                  :: start(3), end(3)
    integer                  :: ncell, ncell_local, omp_get_thread_num
    integer                  :: kk
    integer                  :: iprocx, iprocy
    integer                  :: ix_start, ix_end, iz_start, iz_end

    complex(wp), allocatable :: work1(:), work2(:), work3(:), work4(:)
    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: nonb_lj6_factor(:)
    integer,         pointer :: natom(:)
    integer,         pointer :: atmcls(:,:)


    coord  => domain%coord
    charge => domain%charge
    natom  => domain%num_atom
    atmcls => domain%atom_cls_no
    nonb_lj6_factor => enefunc%nonb_lj6_factor

    ncell  = domain%num_cell_local + domain%num_cell_boundary
    ncell_local  = domain%num_cell_local

    ! Initializing the energy and force
    !
    qdf_real(1:nlocalx*nlocaly*nlocalz) = 0.0_wp
    ldf_real(1:nlocalx*nlocaly*nlocalz) = 0.0_wp

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, iv, ix, iy, iz, icel, k, ii, dv, izz, izs, iyy, iys,  &
    !$omp         iz_start, iz_end, ix_start, ix_end, ixx, ixs, vxx, vyy, vzz, &
    !$omp         iproc, tqq, vr, work1, work2, elec_temp, evdw_temp, temp,    &
    !$omp         force_local, force_local_lj, ixyz, nix, nix1, niy, niz,      &
    !$omp         bscx, bscy, bscz, bscdx, bscdy, bscdz, kk, iprocx, iprocy,   &
    !$omp         qtmp, grid, half_grid, start, end, nizx, nizy, iix, iiy,     &
    !$omp         iiz, u_2, u_3, vyyd, vzzd, iatmcls, ltmp, efft_real,         &
    !$omp         efft_imag, vfft_real, vfft_imag, tll, work3, work4)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    ! initialization
    !
    allocate(work1(ngridmax), work2(2*ngridmax),  &
             work3(ngridmax), work4(2*ngridmax))

    start(1) = x_start
    start(2) = y_start
    start(3) = z_start
    end(1) = x_end
    end(2) = y_end
    end(3) = z_end

    ! Initialization of Q-fanction
    !
    do iz = 1, nlocalz+8
      do iy = 1, nlocaly+8
        do ix = 1, nlocalx+8
          qdf(ix,iy,iz,id+1) = 0.0_wp
          ldf(ix,iy,iz,id+1) = 0.0_wp
        end do
      end do
    end do

    ! Calculateing Q-fanction
    !
    grid(1:3) = real(ngrid(1:3),wp)
    half_grid(1:3) = real(ngrid(1:3)/2,wp)
    grid(4) = 2.0_wp
    half_grid(4) = 1.0_wp
    ngrid(4) = 100000

    do icel = id+1, ncell_local, nthread

      do i = 1, natom(icel)

        vr(1) = coord(1,i,icel) * r_scale(1)
        vr(2) = coord(2,i,icel) * r_scale(2)
        vr(3) = coord(3,i,icel) * r_scale(3)
        vr(1) = vr(1) + half_grid(1) - grid(1)*anint(vr(1)/grid(1))
        vr(2) = vr(2) + half_grid(2) - grid(2)*anint(vr(2)/grid(2))
        vr(3) = vr(3) + half_grid(3) - grid(3)*anint(vr(3)/grid(3))
        vr(1) = min(vr(1),grid(1)-0.0000001_wp)
        vr(2) = min(vr(2),grid(2)-0.0000001_wp)
        vr(3) = min(vr(3),grid(3)-0.0000001_wp)
        vi(1,i,icel) = min(int(vr(1)),ngrid(1)-1)
        vi(2,i,icel) = min(int(vr(2)),ngrid(2)-1)
        vi(3,i,icel) = min(int(vr(3)),ngrid(3)-1)
        dv(1) = vr(1) - real(vi(1,i,icel),wp)
        dv(2) = vr(2) - real(vi(2,i,icel),wp)
        dv(3) = vr(3) - real(vi(3,i,icel),wp)
 
        u_2(1) = dv(1)*dv(1)
        u_2(2) = dv(2)*dv(2)
        u_2(3) = dv(3)*dv(3)
        u_3(1) = u_2(1)*dv(1)
        u_3(2) = u_2(2)*dv(2)
        u_3(3) = u_2(3)*dv(3)
        bsc(1,1,i,icel) = u_3(1)
        bsc(1,2,i,icel) = u_3(2)
        bsc(1,3,i,icel) = u_3(3)
        bsc(2,1,i,icel) = -3.0_wp*u_3(1) + 3.0_wp*u_2(1) &
                          + 3.0_wp*dv(1) + 1.0_wp
        bsc(2,2,i,icel) = -3.0_wp*u_3(2) + 3.0_wp*u_2(2) &
                          + 3.0_wp*dv(2) + 1.0_wp
        bsc(2,3,i,icel) = -3.0_wp*u_3(3) + 3.0_wp*u_2(3) &
                          + 3.0_wp*dv(3) + 1.0_wp
        bsc(3,1,i,icel) = 3.0_wp*u_3(1) - 6.0_wp*u_2(1) + 4.0_wp
        bsc(3,2,i,icel) = 3.0_wp*u_3(2) - 6.0_wp*u_2(2) + 4.0_wp
        bsc(3,3,i,icel) = 3.0_wp*u_3(3) - 6.0_wp*u_2(3) + 4.0_wp
        bsc(4,1,i,icel) = (1.0_wp-dv(1))*(1.0_wp-dv(1))*(1.0_wp-dv(1))
        bsc(4,2,i,icel) = (1.0_wp-dv(2))*(1.0_wp-dv(2))*(1.0_wp-dv(2))
        bsc(4,3,i,icel) = (1.0_wp-dv(3))*(1.0_wp-dv(3))*(1.0_wp-dv(3))
        bscd(1,1,i,icel) = u_2(1)
        bscd(1,2,i,icel) = u_2(2)
        bscd(1,3,i,icel) = u_2(3)
        bscd(2,1,i,icel) = -3.0_wp*u_2(1) + 2.0_wp*dv(1) + 1.0_wp
        bscd(2,2,i,icel) = -3.0_wp*u_2(2) + 2.0_wp*dv(2) + 1.0_wp
        bscd(2,3,i,icel) = -3.0_wp*u_2(3) + 2.0_wp*dv(3) + 1.0_wp
        bscd(3,1,i,icel) = 3.0_wp*u_2(1) - 4.0_wp*dv(1)
        bscd(3,2,i,icel) = 3.0_wp*u_2(2) - 4.0_wp*dv(2)
        bscd(3,3,i,icel) = 3.0_wp*u_2(3) - 4.0_wp*dv(3)
        bscd(4,1,i,icel) = -(1.0_wp-dv(1))*(1.0_wp-dv(1))
        bscd(4,2,i,icel) = -(1.0_wp-dv(2))*(1.0_wp-dv(2))
        bscd(4,3,i,icel) = -(1.0_wp-dv(3))*(1.0_wp-dv(3))
  
      end do

    end do

    do icel = id+1, ncell_local, nthread
      do i = 1, natom(icel)

        qtmp = charge(i,icel)*bs_fact3
        iatmcls = atmcls(i,icel)
        ltmp = nonb_lj6_factor(iatmcls)*bs_fact3
        ii(1:3) = vi(1:3,i,icel)

        do ixyz = 1, 4
          iz = ii(3) - ixyz + 2
          if (iz <= 0) iz = iz + ngrid(3)
          vii(ixyz,3,i,icel) = grid_g2lz(iz)
          iy = ii(2) - ixyz + 2
          if (iy <= 0) iy = iy + ngrid(2)
          vii(ixyz,2,i,icel) = grid_g2ly(iy)
          ix = ii(1) - ixyz + 2
          if (ix <= 0) ix = ix + ngrid(1)
          vii(ixyz,1,i,icel) = grid_g2lx(ix)
        end do

        do izz = 1, 4
          izs = vii(izz,3,i,icel)
          vzz = bsc(izz,3,i,icel)
          do iyy = 1, 4
            iys = vii(iyy,2,i,icel)
            vyy = bsc(iyy,2,i,icel)
!ocl norecurrence(qdf)
            do ixx = 1, 4
              ixs = vii(ixx,1,i,icel)
              qdf(ixs,iys,izs,id+1) = qdf(ixs,iys,izs,id+1)            &
                                    + bsc(ixx,1,i,icel)*vyy*vzz*qtmp
              ldf(ixs,iys,izs,id+1) = ldf(ixs,iys,izs,id+1)            &
                                    + bsc(ixx,1,i,icel)*vyy*vzz*ltmp
            end do
          end do
        end do
      end do
    end do

    !$omp barrier
    do iproc = 2, nthread
      do iz = id+1, nlocalz+8,nthread
        do iy = 1, nlocaly+8
          do ix = 1, nlocalx+8
            qdf(ix,iy,iz,1) = qdf(ix,iy,iz,1) + qdf(ix,iy,iz,iproc)
            ldf(ix,iy,iz,1) = ldf(ix,iy,iz,1) + ldf(ix,iy,iz,iproc)
          end do
        end do
      end do
    end do

    !$omp barrier
    !$omp master
    call communicate_pme_pre_lj
    !$omp end master
    !$omp barrier

    do iz = id+1, nlocalz, nthread
      izs = iz + 4
      do iy = 1, nlocaly
        k = (iy-1)*nlocalx + (iz-1)*nlocalx*nlocaly
        iys = iy + 4
        do ix = 1, nlocalx
          ixs = ix + 4
          qdf_real(k+ix) = qdf_real(k+ix) + qdf(ixs,iys,izs,1)
          ldf_real(k+ix) = ldf_real(k+ix) + ldf(ixs,iys,izs,1)
        end do
      end do
    end do

    !$omp barrier
    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdf_real, nlocalx*nlocaly*nlocalz/nprocx, mpi_wp_real, &
                      qdf_work, nlocalx*nlocaly*nlocalz/nprocx, mpi_wp_real, &
                      grid_commx, ierror)
    call mpi_alltoall(ldf_real, nlocalx*nlocaly*nlocalz/nprocx, mpi_wp_real, &
                      ldf_work, nlocalx*nlocaly*nlocalz/nprocx, mpi_wp_real, &
                      grid_commx, ierror)
#else
    qdf_work(1:nlocalx*nlocaly*nlocalz,1) = qdf_real(1:nlocalx*nlocaly*nlocalz)
    ldf_work(1:nlocalx*nlocaly*nlocalz,1) = ldf_real(1:nlocalx*nlocaly*nlocalz)
#endif
    !$omp end master
    !$omp barrier

    nizx = nlocalz / nprocx
    nizy = nlocalz / nprocy
    call fft3d_1d_alltoall_lj(qdf_work, ftqdf, ftqdf, ftqdf, ftqdf_work,       &
                              ftqdf_work, ftqdf_work, ftqdf_work, ftqdf_work2, &
                              ldf_work, ftldf, ftldf, ftldf, ftldf_work,       &
                              ftldf_work, ftldf_work, ftldf_work, ftldf_work2, &
                              work1, work2, work3, work4,                      &
                              nlocalx, nlocaly, nlocalz, nizx, nizy,           &
                              nlocalx1, x_local1, ngrid(1), ngrid(2),          &
                              ngrid(3), nprocx, nprocy, nprocz, id, nthread,   &
                              grid_commx, grid_commy, grid_commz)


    ! Energy calculation
    !
    iproc = my_z_rank + 1
    niy = nlocaly / nprocz

    do ix = 1, nlocalx/2
      do iy = id+1, niy, nthread
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp
        vxx = 0.0_wp
        vyy = 0.0_wp
        vzz = 0.0_wp

        do iz = 1, ngrid(3)
          efft_real = real(ftqdf_work2(iz,ix,iy),wp)
          efft_imag = imag(ftqdf_work2(iz,ix,iy))
          vfft_real = real(ftldf_work2(iz,ix,iy),wp)
          vfft_imag = imag(ftldf_work2(iz,ix,iy))
          tqq = efft_real*efft_real + efft_imag*efft_imag
          tqq = tqq * theta(iz,iy,ix) * vol_fact4
          elec_temp = elec_temp + tqq
          tll = vfft_real*vfft_real + vfft_imag*vfft_imag
          tll = tll * vol_fact4_lj
          evdw_temp = evdw_temp - tll*theta_lj(iz,iy,ix)

          ! virial
          !
          vxx = vxx + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gx(ix) * gx(ix))
          vyy = vyy + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gy(iy) * gy(iy))
          vzz = vzz + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gz(iz) * gz(iz))
          vxx = vxx  &
              - tll*(theta_lj(iz,iy,ix)+vir_fact_lj(iz,iy,ix)*gx(ix)*gx(ix))
          vyy = vyy  &
              - tll*(theta_lj(iz,iy,ix)+vir_fact_lj(iz,iy,ix)*gy(iy)*gy(iy))
          vzz = vzz  &
              - tll*(theta_lj(iz,iy,ix)+vir_fact_lj(iz,iy,ix)*gz(iz)*gz(iz))
        end do

        eelec(id+1) = eelec(id+1) + elec_temp
        evdw (id+1) = evdw (id+1) + evdw_temp
        virial(1,1,id+1) = virial(1,1,id+1) + vxx
        virial(2,2,id+1) = virial(2,2,id+1) + vyy
        virial(3,3,id+1) = virial(3,3,id+1) + vzz
      end do
    end do

    if (my_x_rank == 0) then

      do iy = id+1, niy, nthread

        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp
        vxx = 0.0_wp
        vyy = 0.0_wp
        vzz = 0.0_wp

        do iz = 1, ngrid(3)

          efft_real = real(ftqdf_work2(iz,1,iy),wp)
          efft_imag = imag(ftqdf_work2(iz,1,iy))
          vfft_real = real(ftldf_work2(iz,1,iy),wp)
          vfft_imag = imag(ftldf_work2(iz,1,iy))
          tqq = efft_real*efft_real + efft_imag*efft_imag
          tqq = tqq * theta(iz,iy,1) * vol_fact2
          elec_temp = elec_temp - tqq
          tll = vfft_real*vfft_real + vfft_imag*vfft_imag
          tll = tll * vol_fact2_lj
          evdw_temp = evdw_temp + tll*theta_lj(iz,iy,1)

          vxx = vxx - tqq * (1.0_wp + vir_fact(iz,iy,1) * gx(1) * gx(1))
          vyy = vyy - tqq * (1.0_wp + vir_fact(iz,iy,1) * gy(iy) * gy(iy))
          vzz = vzz - tqq * (1.0_wp + vir_fact(iz,iy,1) * gz(iz) * gz(iz))
          vxx = vxx  &
              + tll*(theta_lj(iz,iy,1)+vir_fact_lj(iz,iy,1)*gx(1 )*gx(1 ))
          vyy = vyy  &
              + tll*(theta_lj(iz,iy,1)+vir_fact_lj(iz,iy,1)*gy(iy)*gy(iy))
          vzz = vzz  &
              + tll*(theta_lj(iz,iy,1)+vir_fact_lj(iz,iy,1)*gz(iz)*gz(iz))

        end do
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw (id+1) = evdw (id+1) + evdw_temp
        virial(1,1,id+1) = virial(1,1,id+1) + vxx
        virial(2,2,id+1) = virial(2,2,id+1) + vyy
        virial(3,3,id+1) = virial(3,3,id+1) + vzz
      end do

      ix = nlocalx/2 + 1

      do iy = id+1, niy, nthread

        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp
        vxx = 0.0_wp
        vyy = 0.0_wp
        vzz = 0.0_wp

        do iz = 1, ngrid(3)
  
          efft_real = real(ftqdf_work2(iz,ix,iy),wp)
          efft_imag = imag(ftqdf_work2(iz,ix,iy))
          vfft_real = real(ftldf_work2(iz,ix,iy),wp)
          vfft_imag = imag(ftldf_work2(iz,ix,iy))
          tqq = efft_real*efft_real + efft_imag*efft_imag
          tqq = tqq * theta(iz,iy,ix) * vol_fact2
          elec_temp = elec_temp + tqq
          tll = vfft_real*vfft_real + vfft_imag*vfft_imag
          tll = tll * vol_fact2_lj
          evdw_temp = evdw_temp - tll*theta_lj(iz,iy,ix)

          vxx = vxx + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gx(ix) * gx(ix))
          vyy = vyy + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gy(iy) * gy(iy))
          vzz = vzz + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gz(iz) * gz(iz))
          vxx = vxx &
              - tll*(theta_lj(iz,iy,ix)+vir_fact_lj(iz,iy,ix)*gx(ix)*gx(ix))
          vyy = vyy &
              - tll*(theta_lj(iz,iy,ix)+vir_fact_lj(iz,iy,ix)*gy(iy)*gy(iy))
          vzz = vzz &
              - tll*(theta_lj(iz,iy,ix)+vir_fact_lj(iz,iy,ix)*gz(iz)*gz(iz))

        end do

        eelec(id+1) = eelec(id+1) + elec_temp
        evdw (id+1) = evdw (id+1) + evdw_temp
        virial(1,1,id+1) = virial(1,1,id+1) + vxx
        virial(2,2,id+1) = virial(2,2,id+1) + vyy
        virial(3,3,id+1) = virial(3,3,id+1) + vzz

      end do
    end if

    ! add self dispersion term
    !
    !$omp master
    if (my_city_rank == 0) then
      temp    = PI*sqrt(PI)*alpha3*inv_volume/(6.0_wp)*enefunc%pme_dispersion_self1
      evdw(1) = evdw(1) - temp
      evdw(1) = evdw(1) + alpha6/12.0_wp*enefunc%pme_dispersion_self2
      virial(1,1,1) = virial(1,1,1) - temp
      virial(2,2,1) = virial(2,2,1) - temp
      virial(3,3,1) = virial(3,3,1) - temp
    end if
    !$omp end master

    ! F^-1[Th]*F^-1[Q] (=X)
    !
    !$omp barrier
    do ix = 1, x_local1
      do iy = id+1, niy, nthread
        k = (iy-1)*ngrid(3) + (ix-1)*ngrid(3)*niy
        do iz = 1, ngrid(3)
          ftqdf_work(k+iz) = ftqdf_work2(iz,ix,iy) &
                            *cmplx(theta(iz,iy,ix),0.0_wp)
          ftldf_work(k+iz) = ftldf_work2(iz,ix,iy) &
                            *cmplx(theta_lj(iz,iy,ix),0.0_wp)
        end do
      end do
    end do

    !$omp barrier
    call bfft3d_1d_alltoall_lj(qdf_work, qdf_real, ftqdf, ftqdf, ftqdf_work, &
                               ftqdf_work, ftqdf_work2, ftqdf_work,          &
                               ldf_work, ldf_real, ftldf, ftldf, ftldf_work, &
                               ftldf_work, ftldf_work2, ftldf_work,          &
                               work1, work2, work3, work4,                   &
                               nlocalx, nlocaly, nlocalz, nlocalx1,          &
                               x_local1, ngrid(1), ngrid(2), ngrid(3),       &
                               nprocx, nprocy, nprocz, id, nthread,          &
                               my_x_rank, grid_commx, grid_commy,            &
                               grid_commz, niy, nizy, nizx)

    ! X is saved on qdf
    !
    do iz = id+1, nlocalz, nthread
      izs = iz + 4
      do iy = 1, nlocaly
        iys = iy + 4
        k = (iy-1)*nlocalx + (iz-1)*nlocalx*nlocaly
        do ix = 1, nlocalx
          ixs = ix + 4
          qdf(ixs,iys,izs,1) = qdf_real(k+ix)
          ldf(ixs,iys,izs,1) = ldf_real(k+ix)
        end do
      end do
    end do

    !$omp barrier
    !$omp master
    call communicate_pme_post_lj
    !$omp end master
    !$omp barrier

    do icel = id+1, ncell_local, nthread
      do i = 1, natom(icel)

        force_local(1:3) = 0.0_wp
        force_local_lj(1:3) = 0.0_wp
        qtmp = charge(i,icel)*vol_fact4*bs_fact3d
        iatmcls = atmcls(i,icel)
        ltmp = nonb_lj6_factor(iatmcls)*vol_fact4_lj*bs_fact3d

        do izz = 1, 4
          izs = vii(izz,3,i,icel)
          vzz = bsc(izz,3,i,icel)
          vzzd = bscd(izz,3,i,icel)
          do iyy = 1, 4
            iys = vii(iyy,2,i,icel)
            vyy = bsc(iyy,2,i,icel)
            vyyd = bscd(iyy,2,i,icel)
            do ixx = 1, 4
              ixs = vii(ixx,1,i,icel)
              force_local(1)    = force_local(1) &
                                + bscd(ixx,1,i,icel)*vyy*vzz*qdf(ixs,iys,izs,1)
              force_local(2)    = force_local(2)     &
                                + bsc(ixx,1,i,icel)*vyyd*vzz*qdf(ixs,iys,izs,1)
              force_local(3)    = force_local(3)     &
                                + bsc(ixx,1,i,icel)*vyy*vzzd*qdf(ixs,iys,izs,1)
              force_local_lj(1) = force_local_lj(1) &
                                + bscd(ixx,1,i,icel)*vyy*vzz*ldf(ixs,iys,izs,1)
              force_local_lj(2) = force_local_lj(2)     &
                                + bsc(ixx,1,i,icel)*vyyd*vzz*ldf(ixs,iys,izs,1)
              force_local_lj(3) = force_local_lj(3)     &
                                + bsc(ixx,1,i,icel)*vyy*vzzd*ldf(ixs,iys,izs,1)
            end do
          end do
        end do

        force(1:3,i,icel) = - force_local(1:3)*qtmp*r_scale(1:3)    &
                            + force_local_lj(1:3)*ltmp*r_scale(1:3)

      end do
    end do

    deallocate(work1, work2, work3, work4)

    !$omp end parallel
 
    return

  end subroutine pme_recip_opt_1dalltoall_lj_4

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_recip_opt_1dalltoall_6
  !> @brief        Calculate PME reciprocal part with domain decomposition
  !                (B spline order 6)                                     
  !! @authors      JJ, NT
  !! @param[in]    domain : domain information
  !! @param[inout] force  : forces of target systems
  !! @param[inout] virial : virial term of target systems
  !! @param[inout] eelec  : electrostatic energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_recip_opt_1dalltoall_6(domain, force, virial, eelec)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    real(wip),               intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: elec_temp
    real(wp)                 :: vr(4), dv(4), grid(4), half_grid(4)
    real(wp)                 :: bscx(4), bscy(4), bscz(4)
    real(wp)                 :: bscdx(4), bscdy(4), bscdz(4)
    real(wp)                 :: vxx, vyy, vzz, vyyd, vzzd, temp
    real(wp)                 :: tqq, qtmp, tqq1, tqq2
    real(wp)                 :: force_local(3)
    real(wp)                 :: u_2(3), u_3(3), u_4(3), u_5(3)
    real(wp)                 :: v_1(3), v_4(3), v_5(3)
    integer                  :: i, k, k1, icel, iproc, id
    integer                  :: ix, iv, iy, iz, ixs, iys, izs
    integer                  :: ixyz, nix, nix1, niy, niz, nizx, nizy
    integer                  :: ixx, iyy, izz, ii(4)
    integer                  :: iix(4), iiy(4), iiz(4)
    integer                  :: start(3), end(3)
    integer                  :: ncell, ncell_local, omp_get_thread_num
    integer                  :: kk, iorg, k_f, k_t
    integer                  :: iprocx, iprocy
    integer                  :: ix_start, ix_end, iz_start, iz_end

    complex(wp), allocatable :: work1(:), work2(:)
    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    integer,         pointer :: natom(:)


    coord  => domain%coord
    charge => domain%charge
    natom  => domain%num_atom

    ncell  = domain%num_cell_local + domain%num_cell_boundary
    ncell_local  = domain%num_cell_local

    ! Initializing the energy and force
    !
    qdf_real(1:nlocalx*nlocaly*nlocalz) = 0.0_wp

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, iv, ix, iy, iz, icel, k, ii, dv, izz, izs, iyy, iys,  &
    !$omp         iz_start, iz_end, ix_start, ix_end, iorg, k_f, k_t,          &
    !$omp         ixx, ixs, vxx, vyy, vzz, iproc, tqq, k1, temp, vr,           &
    !$omp         work1, work2, elec_temp, force_local, ixyz, tqq1, tqq2,      &
    !$omp         nix, nix1, niy, niz, bscx, bscy, bscz, bscdx, bscdy,         &
    !$omp         bscdz, kk, iprocx, iprocy, qtmp, grid, half_grid, start,     &
    !$omp         end, nizx, nizy, iix, iiy, iiz, u_2, u_3, vyyd, vzzd,        &
    !$omp         u_4, u_5, v_1, v_4, v_5)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    ! initialization
    !
    allocate(work1(ngridmax), work2(2*ngridmax))

    start(1) = x_start
    start(2) = y_start
    start(3) = z_start
    end(1) = x_end
    end(2) = y_end
    end(3) = z_end

    ! Initialization of Q-fanction
    !
    do iz = 1, nlocalz+12
      do iy = 1, nlocaly+12
        do ix = 1, nlocalx+12
          qdf(ix,iy,iz,id+1) = 0.0_wp
        end do
      end do
    end do

    ! Calculateing Q-fanction
    !
    grid(1:3) = real(ngrid(1:3),wp)
    half_grid(1:3) = real(ngrid(1:3)/2,wp)
    grid(4) = 2.0_wp
    half_grid(4) = 1.0_wp
    ngrid(4) = 100000

    do icel = id+1, ncell_local, nthread

      do i = 1, natom(icel)

        vr(1) = coord(1,i,icel) * r_scale(1)
        vr(2) = coord(2,i,icel) * r_scale(2)
        vr(3) = coord(3,i,icel) * r_scale(3)
        vr(1) = vr(1) + half_grid(1) - grid(1)*anint(vr(1)/grid(1))
        vr(2) = vr(2) + half_grid(2) - grid(2)*anint(vr(2)/grid(2))
        vr(3) = vr(3) + half_grid(3) - grid(3)*anint(vr(3)/grid(3))
        vr(1) = min(vr(1),grid(1)-0.0000001_wp)
        vr(2) = min(vr(2),grid(2)-0.0000001_wp)
        vr(3) = min(vr(3),grid(3)-0.0000001_wp)
        vi(1,i,icel) = min(int(vr(1)),ngrid(1)-1)
        vi(2,i,icel) = min(int(vr(2)),ngrid(2)-1)
        vi(3,i,icel) = min(int(vr(3)),ngrid(3)-1)
        dv(1) = vr(1) - real(vi(1,i,icel),wp)
        dv(2) = vr(2) - real(vi(2,i,icel),wp)
        dv(3) = vr(3) - real(vi(3,i,icel),wp)

        v_1(1) = 1.0_wp - dv(1) 
        v_1(2) = 1.0_wp - dv(2) 
        v_1(3) = 1.0_wp - dv(3) 
        v_4(1) = v_1(1)**4
        v_4(2) = v_1(2)**4
        v_4(3) = v_1(3)**4
        v_5(1) = v_4(1)*v_1(1)
        v_5(2) = v_4(2)*v_1(2)
        v_5(3) = v_4(3)*v_1(3)
        u_2(1) = dv(1)*dv(1)
        u_2(2) = dv(2)*dv(2)
        u_2(3) = dv(3)*dv(3)
        u_3(1) = u_2(1)*dv(1)
        u_3(2) = u_2(2)*dv(2)
        u_3(3) = u_2(3)*dv(3)
        u_4(1) = u_2(1)*u_2(1)
        u_4(2) = u_2(2)*u_2(2)
        u_4(3) = u_2(3)*u_2(3)
        u_5(1) = u_2(1)*u_3(1)
        u_5(2) = u_2(2)*u_3(2)
        u_5(3) = u_2(3)*u_3(3)
        bsc(1,1,i,icel) = u_5(1)
        bsc(1,2,i,icel) = u_5(2)
        bsc(1,3,i,icel) = u_5(3)
        bsc(2,1,i,icel) = -5.0_wp*u_5(1)+5.0_wp*u_4(1)+10.0_wp*u_3(1) &
                        + 10.0_wp*u_2(1)+5.0_wp*dv(1)+1.0_wp
        bsc(2,2,i,icel) = -5.0_wp*u_5(2)+5.0_wp*u_4(2)+10.0_wp*u_3(2) &
                        + 10.0_wp*u_2(2)+5.0_wp*dv(2)+1.0_wp
        bsc(2,3,i,icel) = -5.0_wp*u_5(3)+5.0_wp*u_4(3)+10.0_wp*u_3(3) &
                        + 10.0_wp*u_2(3)+5.0_wp*dv(3)+1.0_wp
        bsc(3,1,i,icel) = 10.0_wp*u_5(1)-20.0_wp*u_4(1)-20.0_wp*u_3(1) &
                        + 20.0_wp*u_2(1)+50.0_wp*dv(1)+26.0_wp
        bsc(3,2,i,icel) = 10.0_wp*u_5(2)-20.0_wp*u_4(2)-20.0_wp*u_3(2) &
                        + 20.0_wp*u_2(2)+50.0_wp*dv(2)+26.0_wp
        bsc(3,3,i,icel) = 10.0_wp*u_5(3)-20.0_wp*u_4(3)-20.0_wp*u_3(3) &
                        + 20.0_wp*u_2(3)+50.0_wp*dv(3)+26.0_wp
        bsc(4,1,i,icel) = -10.0_wp*u_5(1)+30.0_wp*u_4(1)-60.0_wp*u_2(1) &
                        + 66.0_wp
        bsc(4,2,i,icel) = -10.0_wp*u_5(2)+30.0_wp*u_4(2)-60.0_wp*u_2(2) &
                        + 66.0_wp
        bsc(4,3,i,icel) = -10.0_wp*u_5(3)+30.0_wp*u_4(3)-60.0_wp*u_2(3) &
                        + 66.0_wp
        bsc(5,1,i,icel) = 5.0_wp*u_5(1)-20.0_wp*u_4(1)+20.0_wp*u_3(1) &
                        + 20.0_wp*u_2(1)-50.0_wp*dv(1)+26.0_wp
        bsc(5,2,i,icel) = 5.0_wp*u_5(2)-20.0_wp*u_4(2)+20.0_wp*u_3(2) &
                        + 20.0_wp*u_2(2)-50.0_wp*dv(2)+26.0_wp
        bsc(5,3,i,icel) = 5.0_wp*u_5(3)-20.0_wp*u_4(3)+20.0_wp*u_3(3) &
                        + 20.0_wp*u_2(3)-50.0_wp*dv(3)+26.0_wp
        bsc(6,1,i,icel) = v_5(1)
        bsc(6,2,i,icel) = v_5(2)
        bsc(6,3,i,icel) = v_5(3)
        bscd(1,1,i,icel) = u_4(1)
        bscd(1,2,i,icel) = u_4(2)
        bscd(1,3,i,icel) = u_4(3)
        bscd(2,1,i,icel) = -5.0_wp*u_4(1)+4.0_wp*u_3(1)+6.0_wp*u_2(1) &
                         + 4.0_wp*dv(1)+1.0_wp
        bscd(2,2,i,icel) = -5.0_wp*u_4(2)+4.0_wp*u_3(2)+6.0_wp*u_2(2) &
                         + 4.0_wp*dv(2)+1.0_wp
        bscd(2,3,i,icel) = -5.0_wp*u_4(3)+4.0_wp*u_3(3)+6.0_wp*u_2(3) &
                         + 4.0_wp*dv(3)+1.0_wp
        bscd(3,1,i,icel) = 10.0_wp*u_4(1)-16.0_wp*u_3(1)-12.0_wp*u_2(1) &
                         + 8.0_wp*dv(1)+10.0_wp
        bscd(3,2,i,icel) = 10.0_wp*u_4(2)-16.0_wp*u_3(2)-12.0_wp*u_2(2) &
                         + 8.0_wp*dv(2)+10.0_wp
        bscd(3,3,i,icel) = 10.0_wp*u_4(3)-16.0_wp*u_3(3)-12.0_wp*u_2(3) &
                         + 8.0_wp*dv(3)+10.0_wp
        bscd(4,1,i,icel) = -10.0_wp*u_4(1)+24.0_wp*u_3(1)-24.0_wp*dv(1)
        bscd(4,2,i,icel) = -10.0_wp*u_4(2)+24.0_wp*u_3(2)-24.0_wp*dv(2)
        bscd(4,3,i,icel) = -10.0_wp*u_4(3)+24.0_wp*u_3(3)-24.0_wp*dv(3)
        bscd(5,1,i,icel) = 5.0_wp*u_4(1)-16.0_wp*u_3(1)+12.0_wp*u_2(1) &
                         + 8.0_wp*dv(1)-10.0_wp
        bscd(5,2,i,icel) = 5.0_wp*u_4(2)-16.0_wp*u_3(2)+12.0_wp*u_2(2) &
                         + 8.0_wp*dv(2)-10.0_wp
        bscd(5,3,i,icel) = 5.0_wp*u_4(3)-16.0_wp*u_3(3)+12.0_wp*u_2(3) &
                         + 8.0_wp*dv(3)-10.0_wp
        bscd(6,1,i,icel) = -v_4(1)
        bscd(6,2,i,icel) = -v_4(2)
        bscd(6,3,i,icel) = -v_4(3)
  
      end do

    end do

    do icel = id+1, ncell_local, nthread
      do i = 1, natom(icel)

        qtmp = charge(i,icel)*bs_fact3
        ii(1:3) = vi(1:3,i,icel)

        do ixyz = 1, 6
          iz = ii(3) - ixyz + 2
          if (iz <= 0) iz = iz + ngrid(3)
          vii(ixyz,3,i,icel) = grid_g2lz(iz)
          iy = ii(2) - ixyz + 2
          if (iy <= 0) iy = iy + ngrid(2)
          vii(ixyz,2,i,icel) = grid_g2ly(iy)
          ix = ii(1) - ixyz + 2
          if (ix <= 0) ix = ix + ngrid(1)
          vii(ixyz,1,i,icel) = grid_g2lx(ix)
        end do

        do izz = 1, 6
          izs = vii(izz,3,i,icel)
          vzz = bsc(izz,3,i,icel)
          do iyy = 1, 6
            iys = vii(iyy,2,i,icel)
            vyy = bsc(iyy,2,i,icel)
!ocl norecurrence(qdf)
            do ixx = 1, 6
              ixs = vii(ixx,1,i,icel)
              qdf(ixs,iys,izs,id+1) = qdf(ixs,iys,izs,id+1)            &
                                    + bsc(ixx,1,i,icel)*vyy*vzz*qtmp
            end do
          end do
        end do
      end do
    end do

    !$omp barrier
    do iproc = 2, nthread
      do iz = id+1, nlocalz+12,nthread
        do iy = 1, nlocaly+12
          do ix = 1, nlocalx+12
            qdf(ix,iy,iz,1) = qdf(ix,iy,iz,1) + qdf(ix,iy,iz,iproc)
          end do
        end do
      end do
    end do

    !$omp barrier
    !$omp master
    call communicate_pme_pre
    !$omp end master
    !$omp barrier

    do iz = id+1, nlocalz, nthread
      izs = iz + 6
      do iy = 1, nlocaly
        k = (iy-1)*nlocalx + (iz-1)*nlocalx*nlocaly
        iys = iy + 6
        do ix = 1, nlocalx
          ixs = ix + 6
          qdf_real(k+ix) = qdf_real(k+ix) + qdf(ixs,iys,izs,1)
        end do
      end do
    end do

    !$omp barrier
    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdf_real, nlocalx*nlocaly*nlocalz/nprocx, mpi_wp_real, &
                      qdf_work, nlocalx*nlocaly*nlocalz/nprocx, mpi_wp_real, &
                      grid_commx, ierror)
#else
    qdf_work(1:nlocalx*nlocaly*nlocalz,1) = qdf_real(1:nlocalx*nlocaly*nlocalz)
#endif
    !$omp end master
    !$omp barrier

    nizx = nlocalz / nprocx
    nizy = nlocalz / nprocy
    call fft3d_1d_alltoall(qdf_work, ftqdf, ftqdf, ftqdf, ftqdf_work,        &
                           ftqdf_work, ftqdf_work, ftqdf_work, ftqdf_work2,  &
                           work1, work2, nlocalx, nlocaly, nlocalz,          &
                           nizx, nizy, nlocalx1, x_local1, ngrid(1),         &
                           ngrid(2), ngrid(3), nprocx, nprocy, nprocz, id,   &
                           nthread, grid_commx, grid_commy, grid_commz)


    ! Energy calculation
    !
    iproc = my_z_rank + 1
    niy = nlocaly / nprocz

    do ix = 1, nlocalx/2
      do iy = id+1, niy, nthread
        elec_temp = 0.0_wp
        vxx = 0.0_wp
        vyy = 0.0_wp
        vzz = 0.0_wp
        do iz = 1, ngrid(3)
          tqq1 = real(ftqdf_work2(iz,ix,iy),wp)
          tqq2 = imag(ftqdf_work2(iz,ix,iy))
          tqq = tqq1*tqq1 + tqq2*tqq2
          tqq = tqq * theta(iz,iy,ix) * vol_fact4
          elec_temp = elec_temp + tqq

          ! virial
          !
          vxx = vxx + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gx(ix) * gx(ix))
          vyy = vyy + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gy(iy) * gy(iy))
          vzz = vzz + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gz(iz) * gz(iz))
        end do
        eelec(id+1) = eelec(id+1) + elec_temp
        virial(1,1,id+1) = virial(1,1,id+1) + vxx
        virial(2,2,id+1) = virial(2,2,id+1) + vyy
        virial(3,3,id+1) = virial(3,3,id+1) + vzz
      end do
    end do

    if (my_x_rank == 0) then
      do iy = id+1, niy, nthread
        elec_temp = 0.0_wp
        vxx = 0.0_wp
        vyy = 0.0_wp
        vzz = 0.0_wp
        do iz = 1, ngrid(3)
          tqq1 = real(ftqdf_work2(iz,1,iy),wp)
          tqq2 = imag(ftqdf_work2(iz,1,iy))
          tqq  = tqq1*tqq1 + tqq2*tqq2
          tqq = tqq * theta(iz,iy,1) * vol_fact2
          elec_temp = elec_temp - tqq
          vxx = vxx - tqq * (1.0_wp + vir_fact(iz,iy,1) * gx(1) * gx(1))
          vyy = vyy - tqq * (1.0_wp + vir_fact(iz,iy,1) * gy(iy) * gy(iy))
          vzz = vzz - tqq * (1.0_wp + vir_fact(iz,iy,1) * gz(iz) * gz(iz))
        end do
        eelec(id+1) = eelec(id+1) + elec_temp
        virial(1,1,id+1) = virial(1,1,id+1) + vxx
        virial(2,2,id+1) = virial(2,2,id+1) + vyy
        virial(3,3,id+1) = virial(3,3,id+1) + vzz
      end do

      ix = nlocalx/2 + 1
      do iy = id+1, niy, nthread
        elec_temp = 0.0_wp
        vxx = 0.0_wp
        vyy = 0.0_wp
        vzz = 0.0_wp
        do iz = 1, ngrid(3)
          tqq1 = real(ftqdf_work2(iz,ix,iy),wp)
          tqq2 = imag(ftqdf_work2(iz,ix,iy))
          tqq  = tqq1*tqq1 + tqq2*tqq2
          tqq = tqq * theta(iz,iy,ix) * vol_fact2
          elec_temp = elec_temp + tqq
          vxx = vxx + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gx(ix) * gx(ix))
          vyy = vyy + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gy(iy) * gy(iy))
          vzz = vzz + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gz(iz) * gz(iz))
        end do
        eelec(id+1) = eelec(id+1) + elec_temp
        virial(1,1,id+1) = virial(1,1,id+1) + vxx
        virial(2,2,id+1) = virial(2,2,id+1) + vyy
        virial(3,3,id+1) = virial(3,3,id+1) + vzz
      end do
    end if

    !$omp barrier

    ! F^-1[Th]*F^-1[Q] (=X)
    !
    do ix = 1, x_local1
      do iy = id+1, niy, nthread
        k = (iy-1)*ngrid(3) + (ix-1)*ngrid(3)*niy
        do iz = 1, ngrid(3)
          ftqdf_work(k+iz) = ftqdf_work2(iz,ix,iy)*cmplx(theta(iz,iy,ix),0.0_wp)
        end do
      end do
    end do

    !$omp barrier

    call bfft3d_1d_alltoall(qdf_work, qdf_real, ftqdf, ftqdf, ftqdf_work,    &
                    ftqdf_work, ftqdf_work2, ftqdf_work, work1, work2,       &
                    nlocalx, nlocaly, nlocalz, nlocalx1, x_local1, ngrid(1), &
                    ngrid(2), ngrid(3), nprocx, nprocy, nprocz, id, nthread, &
                    my_x_rank, grid_commx, grid_commy, grid_commz, niy,      &
                    nizy, nizx)

    ! X is saved on qdf
    !
    do iz = id+1, nlocalz, nthread
      izs = iz + 6
      do iy = 1, nlocaly
        iys = iy + 6
        k = (iy-1)*nlocalx + (iz-1)*nlocalx*nlocaly
        do ix = 1, nlocalx
          ixs = ix + 6
          qdf(ixs,iys,izs,1) = qdf_real(k+ix)
        end do
      end do
    end do

    !$omp barrier
    !$omp master
    call communicate_pme_post
    !$omp end master
    !$omp barrier

    do icel = id+1, ncell_local, nthread
      do i = 1, natom(icel)

        force_local(1:3) = 0.0_wp
        qtmp = charge(i,icel)*vol_fact4*bs_fact3d

        do izz = 1, 6
          izs = vii(izz,3,i,icel)
          vzz = bsc(izz,3,i,icel)
          vzzd = bscd(izz,3,i,icel)
          do iyy = 1, 6
            iys = vii(iyy,2,i,icel)
            vyy = bsc(iyy,2,i,icel)
            vyyd = bscd(iyy,2,i,icel)
            do ixx = 1, 6
              ixs = vii(ixx,1,i,icel)
              force_local(1) = force_local(1) &
                             + bscd(ixx,1,i,icel)*vyy*vzz*qdf(ixs,iys,izs,1)
              force_local(2) = force_local(2)     &
                             + bsc(ixx,1,i,icel)*vyyd*vzz*qdf(ixs,iys,izs,1)
              force_local(3) = force_local(3)     &
                             + bsc(ixx,1,i,icel)*vyy*vzzd*qdf(ixs,iys,izs,1)
            end do
          end do
        end do

        force(1:3,i,icel) = - force_local(1:3)*qtmp*r_scale(1:3)

      end do
    end do

    deallocate(work1, work2)

    !$omp end parallel
 
    return

  end subroutine pme_recip_opt_1dalltoall_6

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_recip_opt_1dalltoall_lj_6
  !> @brief        Calculate PME reciprocal part with domain decomposition
  !                (B spline order 6)                                     
  !! @authors      JJ, NT
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_recip_opt_1dalltoall_lj_6(domain, enefunc, force, virial, &
                                           eelec, evdw)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wip),               intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw (nthread)

    ! local variables
    real(wp)                 :: elec_temp, evdw_temp, temp
    real(wp)                 :: vr(4), dv(4), grid(4), half_grid(4)
    real(wp)                 :: bscx(4), bscy(4), bscz(4)
    real(wp)                 :: bscdx(4), bscdy(4), bscdz(4)
    real(wp)                 :: vxx, vyy, vzz, vyyd, vzzd
    real(wp)                 :: efft_real, efft_imag
    real(wp)                 :: vfft_real, vfft_imag
    real(wp)                 :: tqq, tll, qtmp, ltmp
    real(wp)                 :: force_local(3), force_local_lj(3)
    real(wp)                 :: vol_factor
    real(wp)                 :: u_2(3), u_3(3), u_4(3), u_5(3)
    real(wp)                 :: v_1(3), v_4(3), v_5(3)
    integer                  :: iatmcls
    integer                  :: i, k, icel, iproc, id
    integer                  :: ix, iv, iy, iz, ixs, iys, izs
    integer                  :: ixyz, nix, nix1, niy, niz, nizx, nizy
    integer                  :: ixx, iyy, izz, ii(4)
    integer                  :: iix(4), iiy(4), iiz(4)
    integer                  :: start(3), end(3)
    integer                  :: ncell, ncell_local, omp_get_thread_num
    integer                  :: kk
    integer                  :: iprocx, iprocy
    integer                  :: ix_start, ix_end, iz_start, iz_end

    complex(wp), allocatable :: work1(:), work2(:), work3(:), work4(:)
    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: nonb_lj6_factor(:)
    integer,         pointer :: natom(:)
    integer,         pointer :: atmcls(:,:)


    coord  => domain%coord
    charge => domain%charge
    natom  => domain%num_atom
    atmcls => domain%atom_cls_no
    nonb_lj6_factor => enefunc%nonb_lj6_factor

    ncell  = domain%num_cell_local + domain%num_cell_boundary
    ncell_local  = domain%num_cell_local

    ! Initializing the energy and force
    !
    qdf_real(1:nlocalx*nlocaly*nlocalz) = 0.0_wp
    ldf_real(1:nlocalx*nlocaly*nlocalz) = 0.0_wp

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, iv, ix, iy, iz, icel, k, ii, dv, izz, izs, iyy, iys,  &
    !$omp         iz_start, iz_end, ix_start, ix_end, ixx, ixs, vxx, vyy, vzz, &
    !$omp         iproc, tqq, vr, work1, work2, elec_temp, evdw_temp, temp,    &
    !$omp         force_local, force_local_lj, ixyz, nix, nix1, niy, niz,      &
    !$omp         bscx, bscy, bscz, bscdx, bscdy, bscdz, kk, iprocx, iprocy,   &
    !$omp         qtmp, grid, half_grid, start, end, nizx, nizy, iix, iiy,     &
    !$omp         iiz, u_2, u_3, vyyd, vzzd, iatmcls, ltmp, efft_real,         &
    !$omp         efft_imag, vfft_real, vfft_imag, u_4, u_5, v_1, v_4, v_5,    &
    !$omp         work3, work4, tll)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    ! initialization
    !
    allocate(work1(ngridmax), work2(2*ngridmax),  &
             work3(ngridmax), work4(2*ngridmax))

    start(1) = x_start
    start(2) = y_start
    start(3) = z_start
    end(1) = x_end
    end(2) = y_end
    end(3) = z_end

    ! Initialization of Q-fanction
    !
    do iz = 1, nlocalz+12
      do iy = 1, nlocaly+12
        do ix = 1, nlocalx+12
          qdf(ix,iy,iz,id+1) = 0.0_wp
          ldf(ix,iy,iz,id+1) = 0.0_wp
        end do
      end do
    end do

    ! Calculateing Q-fanction
    !
    grid(1:3) = real(ngrid(1:3),wp)
    half_grid(1:3) = real(ngrid(1:3)/2,wp)
    grid(4) = 2.0_wp
    half_grid(4) = 1.0_wp
    ngrid(4) = 100000

    do icel = id+1, ncell_local, nthread

      do i = 1, natom(icel)

        vr(1) = coord(1,i,icel) * r_scale(1)
        vr(2) = coord(2,i,icel) * r_scale(2)
        vr(3) = coord(3,i,icel) * r_scale(3)
        vr(1) = vr(1) + half_grid(1) - grid(1)*anint(vr(1)/grid(1))
        vr(2) = vr(2) + half_grid(2) - grid(2)*anint(vr(2)/grid(2))
        vr(3) = vr(3) + half_grid(3) - grid(3)*anint(vr(3)/grid(3))
        vr(1) = min(vr(1),grid(1)-0.0000001_wp)
        vr(2) = min(vr(2),grid(2)-0.0000001_wp)
        vr(3) = min(vr(3),grid(3)-0.0000001_wp)
        vi(1,i,icel) = min(int(vr(1)),ngrid(1)-1)
        vi(2,i,icel) = min(int(vr(2)),ngrid(2)-1)
        vi(3,i,icel) = min(int(vr(3)),ngrid(3)-1)
        dv(1) = vr(1) - real(vi(1,i,icel),wp)
        dv(2) = vr(2) - real(vi(2,i,icel),wp)
        dv(3) = vr(3) - real(vi(3,i,icel),wp)

        v_1(1) = 1.0_wp - dv(1)
        v_1(2) = 1.0_wp - dv(2)
        v_1(3) = 1.0_wp - dv(3)
        v_4(1) = v_1(1)**4
        v_4(2) = v_1(2)**4
        v_4(3) = v_1(3)**4
        v_5(1) = v_4(1)*v_1(1)
        v_5(2) = v_4(2)*v_1(2)
        v_5(3) = v_4(3)*v_1(3)
        u_2(1) = dv(1)*dv(1)
        u_2(2) = dv(2)*dv(2)
        u_2(3) = dv(3)*dv(3)
        u_3(1) = u_2(1)*dv(1)
        u_3(2) = u_2(2)*dv(2)
        u_3(3) = u_2(3)*dv(3)
        u_4(1) = u_2(1)*u_2(1)
        u_4(2) = u_2(2)*u_2(2)
        u_4(3) = u_2(3)*u_2(3)
        u_5(1) = u_2(1)*u_3(1)
        u_5(2) = u_2(2)*u_3(2)
        u_5(3) = u_2(3)*u_3(3)
        bsc(1,1,i,icel) = u_5(1)
        bsc(1,2,i,icel) = u_5(2)
        bsc(1,3,i,icel) = u_5(3)
        bsc(2,1,i,icel) = -5.0_wp*u_5(1)+5.0_wp*u_4(1)+10.0_wp*u_3(1) &
                        + 10.0_wp*u_2(1)+5.0_wp*dv(1)+1.0_wp
        bsc(2,2,i,icel) = -5.0_wp*u_5(2)+5.0_wp*u_4(2)+10.0_wp*u_3(2) &
                        + 10.0_wp*u_2(2)+5.0_wp*dv(2)+1.0_wp
        bsc(2,3,i,icel) = -5.0_wp*u_5(3)+5.0_wp*u_4(3)+10.0_wp*u_3(3) &
                        + 10.0_wp*u_2(3)+5.0_wp*dv(3)+1.0_wp
        bsc(3,1,i,icel) = 10.0_wp*u_5(1)-20.0_wp*u_4(1)-20.0_wp*u_3(1) &
                        + 20.0_wp*u_2(1)+50.0_wp*dv(1)+26.0_wp
        bsc(3,2,i,icel) = 10.0_wp*u_5(2)-20.0_wp*u_4(2)-20.0_wp*u_3(2) &
                        + 20.0_wp*u_2(2)+50.0_wp*dv(2)+26.0_wp
        bsc(3,3,i,icel) = 10.0_wp*u_5(3)-20.0_wp*u_4(3)-20.0_wp*u_3(3) &
                        + 20.0_wp*u_2(3)+50.0_wp*dv(3)+26.0_wp
        bsc(4,1,i,icel) = -10.0_wp*u_5(1)+30.0_wp*u_4(1)-60.0_wp*u_2(1) &
                        + 66.0_wp
        bsc(4,2,i,icel) = -10.0_wp*u_5(2)+30.0_wp*u_4(2)-60.0_wp*u_2(2) &
                        + 66.0_wp
        bsc(4,3,i,icel) = -10.0_wp*u_5(3)+30.0_wp*u_4(3)-60.0_wp*u_2(3) &
                        + 66.0_wp
        bsc(5,1,i,icel) = 5.0_wp*u_5(1)-20.0_wp*u_4(1)+20.0_wp*u_3(1) &
                        + 20.0_wp*u_2(1)-50.0_wp*dv(1)+26.0_wp
        bsc(5,2,i,icel) = 5.0_wp*u_5(2)-20.0_wp*u_4(2)+20.0_wp*u_3(2) &
                        + 20.0_wp*u_2(2)-50.0_wp*dv(2)+26.0_wp
        bsc(5,3,i,icel) = 5.0_wp*u_5(3)-20.0_wp*u_4(3)+20.0_wp*u_3(3) &
                        + 20.0_wp*u_2(3)-50.0_wp*dv(3)+26.0_wp
        bsc(6,1,i,icel) = v_5(1)
        bsc(6,2,i,icel) = v_5(2)
        bsc(6,3,i,icel) = v_5(3)
        bscd(1,1,i,icel) = u_4(1)
        bscd(1,2,i,icel) = u_4(2)
        bscd(1,3,i,icel) = u_4(3)
        bscd(2,1,i,icel) = -5.0_wp*u_4(1)+4.0_wp*u_3(1)+6.0_wp*u_2(1) &
                         + 4.0_wp*dv(1)+1.0_wp
        bscd(2,2,i,icel) = -5.0_wp*u_4(2)+4.0_wp*u_3(2)+6.0_wp*u_2(2) &
                         + 4.0_wp*dv(2)+1.0_wp
        bscd(2,3,i,icel) = -5.0_wp*u_4(3)+4.0_wp*u_3(3)+6.0_wp*u_2(3) &
                         + 4.0_wp*dv(3)+1.0_wp
        bscd(3,1,i,icel) = 10.0_wp*u_4(1)-16.0_wp*u_3(1)-12.0_wp*u_2(1) &
                         + 8.0_wp*dv(1)+10.0_wp
        bscd(3,2,i,icel) = 10.0_wp*u_4(2)-16.0_wp*u_3(2)-12.0_wp*u_2(2) &
                         + 8.0_wp*dv(2)+10.0_wp
        bscd(3,3,i,icel) = 10.0_wp*u_4(3)-16.0_wp*u_3(3)-12.0_wp*u_2(3) &
                         + 8.0_wp*dv(3)+10.0_wp
        bscd(4,1,i,icel) = -10.0_wp*u_4(1)+24.0_wp*u_3(1)-24.0_wp*dv(1)
        bscd(4,2,i,icel) = -10.0_wp*u_4(2)+24.0_wp*u_3(2)-24.0_wp*dv(2)
        bscd(4,3,i,icel) = -10.0_wp*u_4(3)+24.0_wp*u_3(3)-24.0_wp*dv(3)
        bscd(5,1,i,icel) = 5.0_wp*u_4(1)-16.0_wp*u_3(1)+12.0_wp*u_2(1) &
                         + 8.0_wp*dv(1)-10.0_wp
        bscd(5,2,i,icel) = 5.0_wp*u_4(2)-16.0_wp*u_3(2)+12.0_wp*u_2(2) &
                         + 8.0_wp*dv(2)-10.0_wp
        bscd(5,3,i,icel) = 5.0_wp*u_4(3)-16.0_wp*u_3(3)+12.0_wp*u_2(3) &
                         + 8.0_wp*dv(3)-10.0_wp
        bscd(6,1,i,icel) = -v_4(1)
        bscd(6,2,i,icel) = -v_4(2)
        bscd(6,3,i,icel) = -v_4(3)
 
      end do

    end do

    do icel = id+1, ncell_local, nthread
      do i = 1, natom(icel)

        qtmp = charge(i,icel)*bs_fact3
        iatmcls = atmcls(i,icel)
        ltmp = nonb_lj6_factor(iatmcls)*bs_fact3
        ii(1:3) = vi(1:3,i,icel)

        do ixyz = 1, 6
          iz = ii(3) - ixyz + 2
          if (iz <= 0) iz = iz + ngrid(3)
          vii(ixyz,3,i,icel) = grid_g2lz(iz)
          iy = ii(2) - ixyz + 2
          if (iy <= 0) iy = iy + ngrid(2)
          vii(ixyz,2,i,icel) = grid_g2ly(iy)
          ix = ii(1) - ixyz + 2
          if (ix <= 0) ix = ix + ngrid(1)
          vii(ixyz,1,i,icel) = grid_g2lx(ix)
        end do

        do izz = 1, 6
          izs = vii(izz,3,i,icel)
          vzz = bsc(izz,3,i,icel)
          do iyy = 1, 6
            iys = vii(iyy,2,i,icel)
            vyy = bsc(iyy,2,i,icel)
!ocl norecurrence(qdf)
            do ixx = 1, 6
              ixs = vii(ixx,1,i,icel)
              qdf(ixs,iys,izs,id+1) = qdf(ixs,iys,izs,id+1)            &
                                    + bsc(ixx,1,i,icel)*vyy*vzz*qtmp
              ldf(ixs,iys,izs,id+1) = ldf(ixs,iys,izs,id+1)            &
                                    + bsc(ixx,1,i,icel)*vyy*vzz*ltmp
            end do
          end do
        end do
      end do
    end do

    !$omp barrier
    do iproc = 2, nthread
      do iz = id+1, nlocalz+12,nthread
        do iy = 1, nlocaly+12
          do ix = 1, nlocalx+12
            qdf(ix,iy,iz,1) = qdf(ix,iy,iz,1) + qdf(ix,iy,iz,iproc)
            ldf(ix,iy,iz,1) = ldf(ix,iy,iz,1) + ldf(ix,iy,iz,iproc)
          end do
        end do
      end do
    end do

    !$omp barrier
    !$omp master
    call communicate_pme_pre_lj
    !$omp end master
    !$omp barrier

    do iz = id+1, nlocalz, nthread
      izs = iz + 6
      do iy = 1, nlocaly
        k = (iy-1)*nlocalx + (iz-1)*nlocalx*nlocaly
        iys = iy + 6
        do ix = 1, nlocalx
          ixs = ix + 6
          qdf_real(k+ix) = qdf_real(k+ix) + qdf(ixs,iys,izs,1)
          ldf_real(k+ix) = ldf_real(k+ix) + ldf(ixs,iys,izs,1)
        end do
      end do
    end do

    !$omp barrier
    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdf_real, nlocalx*nlocaly*nlocalz/nprocx, mpi_wp_real, &
                      qdf_work, nlocalx*nlocaly*nlocalz/nprocx, mpi_wp_real, &
                      grid_commx, ierror)
    call mpi_alltoall(ldf_real, nlocalx*nlocaly*nlocalz/nprocx, mpi_wp_real, &
                      ldf_work, nlocalx*nlocaly*nlocalz/nprocx, mpi_wp_real, &
                      grid_commx, ierror)
#else
    qdf_work(1:nlocalx*nlocaly*nlocalz,1) = qdf_real(1:nlocalx*nlocaly*nlocalz)
    ldf_work(1:nlocalx*nlocaly*nlocalz,1) = ldf_real(1:nlocalx*nlocaly*nlocalz)
#endif
    !$omp end master
    !$omp barrier

    nizx = nlocalz / nprocx
    nizy = nlocalz / nprocy
    call fft3d_1d_alltoall_lj(qdf_work, ftqdf, ftqdf, ftqdf, ftqdf_work,       &
                              ftqdf_work, ftqdf_work, ftqdf_work, ftqdf_work2, &
                              ldf_work, ftldf, ftldf, ftldf, ftldf_work,       &
                              ftldf_work, ftldf_work, ftldf_work, ftldf_work2, &
                              work1, work2, work3, work4,                      &
                              nlocalx, nlocaly, nlocalz, nizx, nizy,           &
                              nlocalx1, x_local1, ngrid(1), ngrid(2),          &
                              ngrid(3), nprocx, nprocy, nprocz, id, nthread,   &
                              grid_commx, grid_commy, grid_commz)


    ! Energy calculation
    !
    iproc = my_z_rank + 1
    niy = nlocaly / nprocz

    do ix = 1, nlocalx/2
      do iy = id+1, niy, nthread
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp
        vxx = 0.0_wp
        vyy = 0.0_wp
        vzz = 0.0_wp

        do iz = 1, ngrid(3)
          efft_real = real(ftqdf_work2(iz,ix,iy),wp)
          efft_imag = imag(ftqdf_work2(iz,ix,iy))
          vfft_real = real(ftldf_work2(iz,ix,iy),wp)
          vfft_imag = imag(ftldf_work2(iz,ix,iy))
          tqq = efft_real*efft_real + efft_imag*efft_imag
          tqq = tqq * theta(iz,iy,ix) * vol_fact4
          elec_temp = elec_temp + tqq
          tll = vfft_real*vfft_real + vfft_imag*vfft_imag
          tll = tll * vol_fact4_lj
          evdw_temp = evdw_temp - tll*theta_lj(iz,iy,ix)

          ! virial
          !
          vxx = vxx + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gx(ix) * gx(ix))
          vyy = vyy + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gy(iy) * gy(iy))
          vzz = vzz + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gz(iz) * gz(iz))
          vxx = vxx  &
              - tll*(theta_lj(iz,iy,ix)+vir_fact_lj(iz,iy,ix)*gx(ix)*gx(ix))
          vyy = vyy  &
              - tll*(theta_lj(iz,iy,ix)+vir_fact_lj(iz,iy,ix)*gy(ix)*gy(ix))
          vzz = vzz  &
              - tll*(theta_lj(iz,iy,ix)+vir_fact_lj(iz,iy,ix)*gz(ix)*gz(ix))
        end do

        eelec(id+1) = eelec(id+1) + elec_temp
        evdw (id+1) = evdw (id+1) + evdw_temp
        virial(1,1,id+1) = virial(1,1,id+1) + vxx
        virial(2,2,id+1) = virial(2,2,id+1) + vyy
        virial(3,3,id+1) = virial(3,3,id+1) + vzz
      end do
    end do

    if (my_x_rank == 0) then

      do iy = id+1, niy, nthread

        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp
        vxx = 0.0_wp
        vyy = 0.0_wp
        vzz = 0.0_wp

        do iz = 1, ngrid(3)

          efft_real = real(ftqdf_work2(iz,1,iy),wp)
          efft_imag = imag(ftqdf_work2(iz,1,iy))
          vfft_real = real(ftldf_work2(iz,1,iy),wp)
          vfft_imag = imag(ftldf_work2(iz,1,iy))
          tqq = efft_real*efft_real + efft_imag*efft_imag
          tqq = tqq * theta(iz,iy,1) * vol_fact2
          elec_temp = elec_temp - tqq
          tll = vfft_real*vfft_real + vfft_imag*vfft_imag
          tll = tll * vol_fact2_lj
          evdw_temp = evdw_temp + tll*theta_lj(iz,iy,1)

          vxx = vxx - tqq * (1.0_wp + vir_fact(iz,iy,1) * gx(1) * gx(1))
          vyy = vyy - tqq * (1.0_wp + vir_fact(iz,iy,1) * gy(iy) * gy(iy))
          vzz = vzz - tqq * (1.0_wp + vir_fact(iz,iy,1) * gz(iz) * gz(iz))
          vxx = vxx  &
              + tll*(theta_lj(iz,iy,1)+vir_fact_lj(iz,iy,1)*gx(1 )*gx(1 ))
          vyy = vyy  &
              + tll*(theta_lj(iz,iy,1)+vir_fact_lj(iz,iy,1)*gy(iy)*gy(iy))
          vzz = vzz  &
              + tll*(theta_lj(iz,iy,1)+vir_fact_lj(iz,iy,1)*gz(iz)*gz(iz))

        end do
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw (id+1) = evdw (id+1) + evdw_temp
        virial(1,1,id+1) = virial(1,1,id+1) + vxx
        virial(2,2,id+1) = virial(2,2,id+1) + vyy
        virial(3,3,id+1) = virial(3,3,id+1) + vzz
      end do

      ix = nlocalx/2 + 1

      do iy = id+1, niy, nthread

        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp
        vxx = 0.0_wp
        vyy = 0.0_wp
        vzz = 0.0_wp

        do iz = 1, ngrid(3)
  
          efft_real = real(ftqdf_work2(iz,ix,iy),wp)
          efft_imag = imag(ftqdf_work2(iz,ix,iy))
          vfft_real = real(ftldf_work2(iz,ix,iy),wp)
          vfft_imag = imag(ftldf_work2(iz,ix,iy))
          tqq = efft_real*efft_real + efft_imag*efft_imag
          tqq = tqq * theta(iz,iy,ix) * vol_fact2
          elec_temp = elec_temp + tqq
          tll = vfft_real*vfft_real + vfft_imag*vfft_imag
          tll = tll * vol_fact2_lj
          evdw_temp = evdw_temp - tll*theta_lj(iz,iy,ix)

          vxx = vxx + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gx(ix) * gx(ix))
          vyy = vyy + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gy(iy) * gy(iy))
          vzz = vzz + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gz(iz) * gz(iz))
          vxx = vxx &
              - tll*(theta_lj(iz,iy,ix)+vir_fact_lj(iz,iy,ix)*gx(ix)*gx(ix))
          vyy = vyy &
              - tll*(theta_lj(iz,iy,ix)+vir_fact_lj(iz,iy,ix)*gy(iy)*gy(iy))
          vzz = vzz &
              - tll*(theta_lj(iz,iy,ix)+vir_fact_lj(iz,iy,ix)*gz(iz)*gz(iz))

        end do

        eelec(id+1) = eelec(id+1) + elec_temp
        evdw (id+1) = evdw (id+1) + evdw_temp
        virial(1,1,id+1) = virial(1,1,id+1) + vxx
        virial(2,2,id+1) = virial(2,2,id+1) + vyy
        virial(3,3,id+1) = virial(3,3,id+1) + vzz

      end do
    end if

    ! add self dispersion term
    !
    !$omp master
    if (my_city_rank == 0) then
      temp    = PI*sqrt(PI)*alpha3*inv_volume/(6.0_wp)*enefunc%pme_dispersion_self1
      evdw(1) = evdw(1) - temp
      evdw(1) = evdw(1) + alpha6/12.0_wp*enefunc%pme_dispersion_self2
      virial(1,1,1) = virial(1,1,1) - temp
      virial(2,2,1) = virial(2,2,1) - temp
      virial(3,3,1) = virial(3,3,1) - temp
    end if
    !$omp end master

    ! F^-1[Th]*F^-1[Q] (=X)
    !
    !$omp barrier
    do ix = 1, x_local1
      do iy = id+1, niy, nthread
        k = (iy-1)*ngrid(3) + (ix-1)*ngrid(3)*niy
        do iz = 1, ngrid(3)
          ftqdf_work(k+iz) = ftqdf_work2(iz,ix,iy) &
                            *cmplx(theta(iz,iy,ix),0.0_wp)
          ftldf_work(k+iz) = ftldf_work2(iz,ix,iy) &
                            *cmplx(theta_lj(iz,iy,ix),0.0_wp)
        end do
      end do
    end do

    !$omp barrier
    call bfft3d_1d_alltoall_lj(qdf_work, qdf_real, ftqdf, ftqdf, ftqdf_work, &
                               ftqdf_work, ftqdf_work2, ftqdf_work,          &
                               ldf_work, ldf_real, ftldf, ftldf, ftldf_work, &
                               ftldf_work, ftldf_work2, ftldf_work,          &
                               work1, work2, work3, work4,                   &
                               nlocalx, nlocaly, nlocalz, nlocalx1,          &
                               x_local1, ngrid(1), ngrid(2), ngrid(3),       &
                               nprocx, nprocy, nprocz, id, nthread,          &
                               my_x_rank, grid_commx, grid_commy,            &
                               grid_commz, niy, nizy, nizx)

    ! X is saved on qdf
    !
    do iz = id+1, nlocalz, nthread
      izs = iz + 6
      do iy = 1, nlocaly
        iys = iy + 6
        k = (iy-1)*nlocalx + (iz-1)*nlocalx*nlocaly
        do ix = 1, nlocalx
          ixs = ix + 6
          qdf(ixs,iys,izs,1) = qdf_real(k+ix)
          ldf(ixs,iys,izs,1) = ldf_real(k+ix)
        end do
      end do
    end do

    !$omp barrier
    !$omp master
    call communicate_pme_post_lj
    !$omp end master
    !$omp barrier

    do icel = id+1, ncell_local, nthread
      do i = 1, natom(icel)

        force_local(1:3) = 0.0_wp
        force_local_lj(1:3) = 0.0_wp
        qtmp = charge(i,icel)*vol_fact4*bs_fact3d
        iatmcls = atmcls(i,icel)
        ltmp = nonb_lj6_factor(iatmcls)*vol_fact4_lj*bs_fact3d

        do izz = 1, 6
          izs = vii(izz,3,i,icel)
          vzz = bsc(izz,3,i,icel)
          vzzd = bscd(izz,3,i,icel)
          do iyy = 1, 6
            iys = vii(iyy,2,i,icel)
            vyy = bsc(iyy,2,i,icel)
            vyyd = bscd(iyy,2,i,icel)
            do ixx = 1, 6
              ixs = vii(ixx,1,i,icel)
              force_local(1)    = force_local(1) &
                                + bscd(ixx,1,i,icel)*vyy*vzz*qdf(ixs,iys,izs,1)
              force_local(2)    = force_local(2)     &
                                + bsc(ixx,1,i,icel)*vyyd*vzz*qdf(ixs,iys,izs,1)
              force_local(3)    = force_local(3)     &
                                + bsc(ixx,1,i,icel)*vyy*vzzd*qdf(ixs,iys,izs,1)
              force_local_lj(1) = force_local_lj(1) &
                                + bscd(ixx,1,i,icel)*vyy*vzz*ldf(ixs,iys,izs,1)
              force_local_lj(2) = force_local_lj(2)     &
                                + bsc(ixx,1,i,icel)*vyyd*vzz*ldf(ixs,iys,izs,1)
              force_local_lj(3) = force_local_lj(3)     &
                                + bsc(ixx,1,i,icel)*vyy*vzzd*ldf(ixs,iys,izs,1)
            end do
          end do
        end do

        force(1:3,i,icel) = - force_local(1:3)*qtmp*r_scale(1:3)    &
                            + force_local_lj(1:3)*ltmp*r_scale(1:3)

      end do
    end do

    deallocate(work1, work2, work3, work4)

    !$omp end parallel
 
    return

  end subroutine pme_recip_opt_1dalltoall_lj_6

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_recip_opt_1dalltoall_8
  !> @brief        Calculate PME reciprocal part with domain decomposition
  !                (B spline order 8)                                     
  !! @authors      JJ, NT
  !! @param[in]    domain : domain information
  !! @param[inout] force  : forces of target systems
  !! @param[inout] virial : virial term of target systems
  !! @param[inout] eelec  : electrostatic energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_recip_opt_1dalltoall_8(domain, force, virial, eelec)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    real(wip),               intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: elec_temp
    real(wp)                 :: vr(4), dv(4), grid(4), half_grid(4)
    real(wp)                 :: bscx(4), bscy(4), bscz(4)
    real(wp)                 :: bscdx(4), bscdy(4), bscdz(4)
    real(wp)                 :: vxx, vyy, vzz, vyyd, vzzd, temp
    real(wp)                 :: tqq, qtmp, tqq1, tqq2
    real(wp)                 :: force_local(3)
    real(wp)                 :: u_2(3), u_3(3), u_4(3), u_5(3)
    real(wp)                 :: u_6(3), u_7(3)
    real(wp)                 :: v_1(3), v_4(3), v_6(3), v_7(3)
    integer                  :: i, k, k1, icel, iproc, id
    integer                  :: ix, iv, iy, iz, ixs, iys, izs
    integer                  :: ixyz, nix, nix1, niy, niz, nizx, nizy
    integer                  :: ixx, iyy, izz, ii(4)
    integer                  :: iix(4), iiy(4), iiz(4)
    integer                  :: start(3), end(3)
    integer                  :: ncell, ncell_local, omp_get_thread_num
    integer                  :: kk, iorg, k_f, k_t
    integer                  :: iprocx, iprocy
    integer                  :: ix_start, ix_end, iz_start, iz_end

    complex(wp), allocatable :: work1(:), work2(:)
    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    integer,         pointer :: natom(:)


    coord  => domain%coord
    charge => domain%charge
    natom  => domain%num_atom

    ncell  = domain%num_cell_local + domain%num_cell_boundary
    ncell_local  = domain%num_cell_local

    ! Initializing the energy and force
    !
    qdf_real(1:nlocalx*nlocaly*nlocalz) = 0.0_wp

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, iv, ix, iy, iz, icel, k, ii, dv, izz, izs, iyy, iys,  &
    !$omp         iz_start, iz_end, ix_start, ix_end, iorg, k_f, k_t,          &
    !$omp         ixx, ixs, vxx, vyy, vzz, iproc, tqq, k1, temp, vr,           &
    !$omp         work1, work2, elec_temp, force_local, ixyz, tqq1, tqq2,      &
    !$omp         nix, nix1, niy, niz, bscx, bscy, bscz, bscdx, bscdy,         &
    !$omp         bscdz, kk, iprocx, iprocy, qtmp, grid, half_grid, start,     &
    !$omp         end, nizx, nizy, iix, iiy, iiz, u_2, u_3, vyyd, vzzd,        &
    !$omp         u_4, u_5, u_6, u_7, v_1, v_4, v_6, v_7)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    ! initialization
    !
    allocate(work1(ngridmax), work2(2*ngridmax))

    start(1) = x_start
    start(2) = y_start
    start(3) = z_start
    end(1) = x_end
    end(2) = y_end
    end(3) = z_end

    ! Initialization of Q-fanction
    !
    do iz = 1, nlocalz+16
      do iy = 1, nlocaly+16
        do ix = 1, nlocalx+16
          qdf(ix,iy,iz,id+1) = 0.0_wp
        end do
      end do
    end do

    ! Calculateing Q-fanction
    !
    grid(1:3) = real(ngrid(1:3),wp)
    half_grid(1:3) = real(ngrid(1:3)/2,wp)
    grid(4) = 2.0_wp
    half_grid(4) = 1.0_wp
    ngrid(4) = 100000

    do icel = id+1, ncell_local, nthread

      do i = 1, natom(icel)

        vr(1) = coord(1,i,icel) * r_scale(1)
        vr(2) = coord(2,i,icel) * r_scale(2)
        vr(3) = coord(3,i,icel) * r_scale(3)
        vr(1) = vr(1) + half_grid(1) - grid(1)*anint(vr(1)/grid(1))
        vr(2) = vr(2) + half_grid(2) - grid(2)*anint(vr(2)/grid(2))
        vr(3) = vr(3) + half_grid(3) - grid(3)*anint(vr(3)/grid(3))
        vr(1) = min(vr(1),grid(1)-0.0000001_wp)
        vr(2) = min(vr(2),grid(2)-0.0000001_wp)
        vr(3) = min(vr(3),grid(3)-0.0000001_wp)
        vi(1,i,icel) = min(int(vr(1)),ngrid(1)-1)
        vi(2,i,icel) = min(int(vr(2)),ngrid(2)-1)
        vi(3,i,icel) = min(int(vr(3)),ngrid(3)-1)
        dv(1) = vr(1) - real(vi(1,i,icel),wp)
        dv(2) = vr(2) - real(vi(2,i,icel),wp)
        dv(3) = vr(3) - real(vi(3,i,icel),wp)

        v_1(1) = 1.0_wp - dv(1) 
        v_1(2) = 1.0_wp - dv(2) 
        v_1(3) = 1.0_wp - dv(3) 
        v_4(1) = v_1(1)*v_1(1)*v_1(1)*v_1(1)
        v_4(2) = v_1(2)*v_1(2)*v_1(2)*v_1(2)
        v_4(3) = v_1(3)*v_1(3)*v_1(3)*v_1(3)
        v_6(1) = v_4(1)*v_1(1)*v_1(1)
        v_6(2) = v_4(1)*v_1(2)*v_1(2)
        v_6(3) = v_4(1)*v_1(3)*v_1(3)
        v_7(1) = v_6(1)*v_1(1)
        v_7(2) = v_6(2)*v_1(2)
        v_7(3) = v_6(3)*v_1(3)
        u_2(1) = dv(1)*dv(1)
        u_2(2) = dv(2)*dv(2)
        u_2(3) = dv(3)*dv(3)
        u_3(1) = u_2(1)*dv(1)
        u_3(2) = u_2(2)*dv(2)
        u_3(3) = u_2(3)*dv(3)
        u_4(1) = u_2(1)*u_2(1)
        u_4(2) = u_2(2)*u_2(2)
        u_4(3) = u_2(3)*u_2(3)
        u_5(1) = u_2(1)*u_3(1)
        u_5(2) = u_2(2)*u_3(2)
        u_5(3) = u_2(3)*u_3(3)
        u_6(1) = u_3(1)*u_3(1)
        u_6(2) = u_3(2)*u_3(2)
        u_6(3) = u_3(3)*u_3(3)
        u_7(1) = u_3(1)*u_4(1)
        u_7(2) = u_3(2)*u_4(2)
        u_7(3) = u_3(3)*u_4(3)
        bsc(1,1,i,icel) = u_7(1)
        bsc(1,2,i,icel) = u_7(2)
        bsc(1,3,i,icel) = u_7(3)
        bsc(2,1,i,icel) = -7.0_wp*u_7(1)+7.0_wp*u_6(1)+21.0_wp*u_5(1)  &
                        + 35.0_wp*u_4(1)+35.0_wp*u_3(1)+21.0_wp*u_2(1) &
                        + 7.0_wp*dv(1)+1.0_wp
        bsc(2,2,i,icel) = -7.0_wp*u_7(2)+7.0_wp*u_6(2)+21.0_wp*u_5(2)  &
                        + 35.0_wp*u_4(2)+35.0_wp*u_3(2)+21.0_wp*u_2(2) &
                        + 7.0_wp*dv(2)+1.0_wp
        bsc(2,3,i,icel) = -7.0_wp*u_7(3)+7.0_wp*u_6(3)+21.0_wp*u_5(3)  &
                        + 35.0_wp*u_4(3)+35.0_wp*u_3(3)+21.0_wp*u_2(3) &
                        + 7.0_wp*dv(3)+1.0_wp
        bsc(3,1,i,icel) = 21.0_wp*u_7(1)-42.0_wp*u_6(1)-84.0_wp*u_5(1) &
                        + 280.0_wp*u_3(1)+504.0_wp*u_2(1)              &
                        + 392.0_wp*dv(1)+120.0_wp
        bsc(3,2,i,icel) = 21.0_wp*u_7(2)-42.0_wp*u_6(2)-84.0_wp*u_5(2) &
                        + 280.0_wp*u_3(2)+504.0_wp*u_2(2)              &
                        + 392.0_wp*dv(2)+120.0_wp
        bsc(3,3,i,icel) = 21.0_wp*u_7(3)-42.0_wp*u_6(3)-84.0_wp*u_5(3) &
                        + 280.0_wp*u_3(3)+504.0_wp*u_2(3)              &
                        + 392.0_wp*dv(3)+120.0_wp
        bsc(4,1,i,icel) = -35.0_wp*u_7(1)+105.0_wp*u_6(1)+105.0_wp*u_5(1) &
                        - 315.0_wp*u_4(1)-665.0_wp*u_3(1)+315.0_wp*u_2(1) &
                        + 1715.0_wp*dv(1)+1191.0_wp
        bsc(4,2,i,icel) = -35.0_wp*u_7(2)+105.0_wp*u_6(2)+105.0_wp*u_5(2) &
                        - 315.0_wp*u_4(2)-665.0_wp*u_3(2)+315.0_wp*u_2(2) &
                        + 1715.0_wp*dv(2)+1191.0_wp
        bsc(4,3,i,icel) = -35.0_wp*u_7(3)+105.0_wp*u_6(3)+105.0_wp*u_5(3) &
                        - 315.0_wp*u_4(3)-665.0_wp*u_3(3)+315.0_wp*u_2(3) &
                        + 1715.0_wp*dv(3)+1191.0_wp
        bsc(5,1,i,icel) = 35.0_wp*u_7(1)-140.0_wp*u_6(1)+560.0_wp*u_4(1) &
                        - 1680.0_wp*u_2(1)+2416.0_wp
        bsc(5,2,i,icel) = 35.0_wp*u_7(2)-140.0_wp*u_6(2)+560.0_wp*u_4(2) &
                        - 1680.0_wp*u_2(2)+2416.0_wp
        bsc(5,3,i,icel) = 35.0_wp*u_7(3)-140.0_wp*u_6(3)+560.0_wp*u_4(3) &
                        - 1680.0_wp*u_2(3)+2416.0_wp
        bsc(6,1,i,icel) = -21.0_wp*u_7(1)+105.0_wp*u_6(1)-105.0_wp*u_5(1) &
                        - 315.0_wp*u_4(1)+665.0_wp*u_3(1)+315.0_wp*u_2(1) &
                        - 1715.0_wp*dv(1)+1191.0_wp
        bsc(6,2,i,icel) = -21.0_wp*u_7(2)+105.0_wp*u_6(2)-105.0_wp*u_5(2) &
                        - 315.0_wp*u_4(2)+665.0_wp*u_3(2)+315.0_wp*u_2(2) &
                        - 1715.0_wp*dv(2)+1191.0_wp
        bsc(6,3,i,icel) = -21.0_wp*u_7(3)+105.0_wp*u_6(3)-105.0_wp*u_5(3) &
                        - 315.0_wp*u_4(3)+665.0_wp*u_3(3)+315.0_wp*u_2(3) &
                        - 1715.0_wp*dv(3)+1191.0_wp
        bsc(7,1,i,icel) = 7.0_wp*u_7(1)-42.0_wp*u_6(1)+84.0_wp*u_5(1) &
                        - 280.0_wp*u_3(1)+504.0_wp*u_2(1)             &
                        - 392.0_wp*dv(1)+120.0_wp
        bsc(7,2,i,icel) = 7.0_wp*u_7(2)-42.0_wp*u_6(2)+84.0_wp*u_5(2) &
                        - 280.0_wp*u_3(2)+504.0_wp*u_2(2)             &
                        - 392.0_wp*dv(2)+120.0_wp
        bsc(7,3,i,icel) = 7.0_wp*u_7(3)-42.0_wp*u_6(3)+84.0_wp*u_5(3) &
                        - 280.0_wp*u_3(3)+504.0_wp*u_2(3)             &
                        - 392.0_wp*dv(3)+120.0_wp
        bsc(8,1,i,icel) = v_7(1)
        bsc(8,2,i,icel) = v_7(2)
        bsc(8,3,i,icel) = v_7(3)
        bscd(1,1,i,icel) = u_6(1)
        bscd(1,2,i,icel) = u_6(2)
        bscd(1,3,i,icel) = u_6(3)
        bscd(2,1,i,icel) = -7.0_wp*u_6(1)+6.0_wp*u_5(1)+15.0_wp*u_4(1) &
                         + 20.0_wp*u_3(1)+15.0_wp*u_2(1)+6.0_wp*dv(1)  &
                         + 1.0_wp
        bscd(2,2,i,icel) = -7.0_wp*u_6(2)+6.0_wp*u_5(2)+15.0_wp*u_4(2) &
                         + 20.0_wp*u_3(2)+15.0_wp*u_2(2)+6.0_wp*dv(2)  &
                         + 1.0_wp
        bscd(2,3,i,icel) = -7.0_wp*u_6(3)+6.0_wp*u_5(3)+15.0_wp*u_4(3) &
                         + 20.0_wp*u_3(3)+15.0_wp*u_2(3)+6.0_wp*dv(3)  &
                         + 1.0_wp
        bscd(3,1,i,icel) = 21.0_wp*u_6(1)-36.0_wp*u_5(1)-60.0_wp*u_4(1) &
                         + 120.0_wp*u_2(1)+144.0_wp*dv(1)+56.0_wp
        bscd(3,2,i,icel) = 21.0_wp*u_6(2)-36.0_wp*u_5(2)-60.0_wp*u_4(2) &
                         + 120.0_wp*u_2(2)+144.0_wp*dv(2)+56.0_wp
        bscd(3,3,i,icel) = 21.0_wp*u_6(3)-36.0_wp*u_5(3)-60.0_wp*u_4(3) &
                         + 120.0_wp*u_2(3)+144.0_wp*dv(3)+56.0_wp
        bscd(4,1,i,icel) = -35.0_wp*u_6(1)+90.0_wp*u_5(1)+75.0_wp*u_4(1) &
                         - 180.0_wp*u_3(1)-285.0_wp*u_2(1)+90.0_wp*dv(1) &
                         + 245.0_wp
        bscd(4,2,i,icel) = -35.0_wp*u_6(2)+90.0_wp*u_5(2)+75.0_wp*u_4(2) &
                         - 180.0_wp*u_3(2)-285.0_wp*u_2(2)+90.0_wp*dv(2) &
                         + 245.0_wp
        bscd(4,3,i,icel) = -35.0_wp*u_6(3)+90.0_wp*u_5(3)+75.0_wp*u_4(3) &
                         - 180.0_wp*u_3(3)-285.0_wp*u_2(3)+90.0_wp*dv(3) &
                         + 245.0_wp
        bscd(5,1,i,icel) = 35.0_wp*u_6(1)-120.0_wp*u_5(1)+320.0_wp*u_3(1) &
                         - 480.0_wp*dv(1)
        bscd(5,2,i,icel) = 35.0_wp*u_6(2)-120.0_wp*u_5(2)+320.0_wp*u_3(2) &
                         - 480.0_wp*dv(2)
        bscd(5,3,i,icel) = 35.0_wp*u_6(3)-120.0_wp*u_5(3)+320.0_wp*u_3(3) &
                         - 480.0_wp*dv(3)
        bscd(6,1,i,icel) = -21.0_wp*u_6(1)+90.0_wp*u_5(1)-75.0_wp*u_4(1) &
                         - 180.0_wp*u_3(1)+285.0_wp*u_2(1)+90.0_wp*dv(1) &
                         -245.0_wp
        bscd(6,2,i,icel) = -21.0_wp*u_6(2)+90.0_wp*u_5(2)-75.0_wp*u_4(2) &
                         - 180.0_wp*u_3(2)+285.0_wp*u_2(2)+90.0_wp*dv(2) &
                         -245.0_wp
        bscd(6,3,i,icel) = -21.0_wp*u_6(3)+90.0_wp*u_5(3)-75.0_wp*u_4(3) &
                         - 180.0_wp*u_3(3)+285.0_wp*u_2(3)+90.0_wp*dv(3) &
                         -245.0_wp
        bscd(7,1,i,icel) = 7.0_wp*u_6(1)-36.0_wp*u_5(1)+60.0_wp*u_4(1) &
                         -120.0_wp*u_2(1)+144.0_wp*dv(1)-56.0_wp
        bscd(7,2,i,icel) = 7.0_wp*u_6(2)-36.0_wp*u_5(2)+60.0_wp*u_4(2) &
                         -120.0_wp*u_2(2)+144.0_wp*dv(2)-56.0_wp
        bscd(7,3,i,icel) = 7.0_wp*u_6(3)-36.0_wp*u_5(3)+60.0_wp*u_4(3) &
                         -120.0_wp*u_2(3)+144.0_wp*dv(3)-56.0_wp
        bscd(8,1,i,icel) = -v_6(1)
        bscd(8,2,i,icel) = -v_6(2)
        bscd(8,3,i,icel) = -v_6(3)
  
      end do

    end do

    do icel = id+1, ncell_local, nthread
      do i = 1, natom(icel)

        qtmp = charge(i,icel)*bs_fact3
        ii(1:3) = vi(1:3,i,icel)

        do ixyz = 1, 8
          iz = ii(3) - ixyz + 2
          if (iz <= 0) iz = iz + ngrid(3)
          vii(ixyz,3,i,icel) = grid_g2lz(iz)
          iy = ii(2) - ixyz + 2
          if (iy <= 0) iy = iy + ngrid(2)
          vii(ixyz,2,i,icel) = grid_g2ly(iy)
          ix = ii(1) - ixyz + 2
          if (ix <= 0) ix = ix + ngrid(1)
          vii(ixyz,1,i,icel) = grid_g2lx(ix)
        end do

        do izz = 1, 8
          izs = vii(izz,3,i,icel)
          vzz = bsc(izz,3,i,icel)
          do iyy = 1, 8
            iys = vii(iyy,2,i,icel)
            vyy = bsc(iyy,2,i,icel)
!ocl norecurrence(qdf)
            do ixx = 1, 8
              ixs = vii(ixx,1,i,icel)
              qdf(ixs,iys,izs,id+1) = qdf(ixs,iys,izs,id+1)            &
                                    + bsc(ixx,1,i,icel)*vyy*vzz*qtmp
            end do
          end do
        end do
      end do
    end do

    !$omp barrier
    do iproc = 2, nthread
      do iz = id+1, nlocalz+16,nthread
        do iy = 1, nlocaly+16
          do ix = 1, nlocalx+16
            qdf(ix,iy,iz,1) = qdf(ix,iy,iz,1) + qdf(ix,iy,iz,iproc)
          end do
        end do
      end do
    end do

    !$omp barrier
    !$omp master
    call communicate_pme_pre
    !$omp end master
    !$omp barrier

    do iz = id+1, nlocalz, nthread
      izs = iz + 8
      do iy = 1, nlocaly
        k = (iy-1)*nlocalx + (iz-1)*nlocalx*nlocaly
        iys = iy + 8
        do ix = 1, nlocalx
          ixs = ix + 8
          qdf_real(k+ix) = qdf_real(k+ix) + qdf(ixs,iys,izs,1)
        end do
      end do
    end do

    !$omp barrier
    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdf_real, nlocalx*nlocaly*nlocalz/nprocx, mpi_wp_real, &
                      qdf_work, nlocalx*nlocaly*nlocalz/nprocx, mpi_wp_real, &
                      grid_commx, ierror)
#else
    qdf_work(1:nlocalx*nlocaly*nlocalz,1) = qdf_real(1:nlocalx*nlocaly*nlocalz)
#endif
    !$omp end master
    !$omp barrier

    nizx = nlocalz / nprocx
    nizy = nlocalz / nprocy
    call fft3d_1d_alltoall(qdf_work, ftqdf, ftqdf, ftqdf, ftqdf_work,        &
                           ftqdf_work, ftqdf_work, ftqdf_work, ftqdf_work2,  &
                           work1, work2, nlocalx, nlocaly, nlocalz,          &
                           nizx, nizy, nlocalx1, x_local1, ngrid(1),         &
                           ngrid(2), ngrid(3), nprocx, nprocy, nprocz, id,   &
                           nthread, grid_commx, grid_commy, grid_commz)


    ! Energy calculation
    !
    iproc = my_z_rank + 1
    niy = nlocaly / nprocz

    do ix = 1, nlocalx/2
      do iy = id+1, niy, nthread
        elec_temp = 0.0_wp
        vxx = 0.0_wp
        vyy = 0.0_wp
        vzz = 0.0_wp
        do iz = 1, ngrid(3)
          tqq1 = real(ftqdf_work2(iz,ix,iy),wp)
          tqq2 = imag(ftqdf_work2(iz,ix,iy))
          tqq = tqq1*tqq1 + tqq2*tqq2
          tqq = tqq * theta(iz,iy,ix) * vol_fact4
          elec_temp = elec_temp + tqq

          ! virial
          !
          vxx = vxx + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gx(ix) * gx(ix))
          vyy = vyy + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gy(iy) * gy(iy))
          vzz = vzz + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gz(iz) * gz(iz))
        end do
        eelec(id+1) = eelec(id+1) + elec_temp
        virial(1,1,id+1) = virial(1,1,id+1) + vxx
        virial(2,2,id+1) = virial(2,2,id+1) + vyy
        virial(3,3,id+1) = virial(3,3,id+1) + vzz
      end do
    end do

    if (my_x_rank == 0) then
      do iy = id+1, niy, nthread
        elec_temp = 0.0_wp
        vxx = 0.0_wp
        vyy = 0.0_wp
        vzz = 0.0_wp
        do iz = 1, ngrid(3)
          tqq1 = real(ftqdf_work2(iz,1,iy),wp)
          tqq2 = imag(ftqdf_work2(iz,1,iy))
          tqq  = tqq1*tqq1 + tqq2*tqq2
          tqq = tqq * theta(iz,iy,1) * vol_fact2
          elec_temp = elec_temp - tqq
          vxx = vxx - tqq * (1.0_wp + vir_fact(iz,iy,1) * gx(1) * gx(1))
          vyy = vyy - tqq * (1.0_wp + vir_fact(iz,iy,1) * gy(iy) * gy(iy))
          vzz = vzz - tqq * (1.0_wp + vir_fact(iz,iy,1) * gz(iz) * gz(iz))
        end do
        eelec(id+1) = eelec(id+1) + elec_temp
        virial(1,1,id+1) = virial(1,1,id+1) + vxx
        virial(2,2,id+1) = virial(2,2,id+1) + vyy
        virial(3,3,id+1) = virial(3,3,id+1) + vzz
      end do

      ix = nlocalx/2 + 1
      do iy = id+1, niy, nthread
        elec_temp = 0.0_wp
        vxx = 0.0_wp
        vyy = 0.0_wp
        vzz = 0.0_wp
        do iz = 1, ngrid(3)
          tqq1 = real(ftqdf_work2(iz,ix,iy),wp)
          tqq2 = imag(ftqdf_work2(iz,ix,iy))
          tqq  = tqq1*tqq1 + tqq2*tqq2
          tqq = tqq * theta(iz,iy,ix) * vol_fact2
          elec_temp = elec_temp + tqq
          vxx = vxx + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gx(ix) * gx(ix))
          vyy = vyy + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gy(iy) * gy(iy))
          vzz = vzz + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gz(iz) * gz(iz))
        end do
        eelec(id+1) = eelec(id+1) + elec_temp
        virial(1,1,id+1) = virial(1,1,id+1) + vxx
        virial(2,2,id+1) = virial(2,2,id+1) + vyy
        virial(3,3,id+1) = virial(3,3,id+1) + vzz
      end do
    end if

    !$omp barrier

    ! F^-1[Th]*F^-1[Q] (=X)
    !
    do ix = 1, x_local1
      do iy = id+1, niy, nthread
        k = (iy-1)*ngrid(3) + (ix-1)*ngrid(3)*niy
        do iz = 1, ngrid(3)
          ftqdf_work(k+iz) = ftqdf_work2(iz,ix,iy)*cmplx(theta(iz,iy,ix),0.0_wp)
        end do
      end do
    end do

    !$omp barrier

    call bfft3d_1d_alltoall(qdf_work, qdf_real, ftqdf, ftqdf, ftqdf_work,    &
                    ftqdf_work, ftqdf_work2, ftqdf_work, work1, work2,       &
                    nlocalx, nlocaly, nlocalz, nlocalx1, x_local1, ngrid(1), &
                    ngrid(2), ngrid(3), nprocx, nprocy, nprocz, id, nthread, &
                    my_x_rank, grid_commx, grid_commy, grid_commz, niy,      &
                    nizy, nizx)

    ! X is saved on qdf
    !
    do iz = id+1, nlocalz, nthread
      izs = iz + 8
      do iy = 1, nlocaly
        iys = iy + 8
        k = (iy-1)*nlocalx + (iz-1)*nlocalx*nlocaly
        do ix = 1, nlocalx
          ixs = ix + 8
          qdf(ixs,iys,izs,1) = qdf_real(k+ix)
        end do
      end do
    end do

    !$omp barrier
    !$omp master
    call communicate_pme_post
    !$omp end master
    !$omp barrier

    do icel = id+1, ncell_local, nthread
      do i = 1, natom(icel)

        force_local(1:3) = 0.0_wp
        qtmp = charge(i,icel)*vol_fact4*bs_fact3d

        do izz = 1, 8
          izs = vii(izz,3,i,icel)
          vzz = bsc(izz,3,i,icel)
          vzzd = bscd(izz,3,i,icel)
          do iyy = 1, 8
            iys = vii(iyy,2,i,icel)
            vyy = bsc(iyy,2,i,icel)
            vyyd = bscd(iyy,2,i,icel)
            do ixx = 1, 8
              ixs = vii(ixx,1,i,icel)
              force_local(1) = force_local(1) &
                             + bscd(ixx,1,i,icel)*vyy*vzz*qdf(ixs,iys,izs,1)
              force_local(2) = force_local(2)     &
                             + bsc(ixx,1,i,icel)*vyyd*vzz*qdf(ixs,iys,izs,1)
              force_local(3) = force_local(3)     &
                             + bsc(ixx,1,i,icel)*vyy*vzzd*qdf(ixs,iys,izs,1)
            end do
          end do
        end do

        force(1:3,i,icel) = - force_local(1:3)*qtmp*r_scale(1:3)

      end do
    end do

    deallocate(work1, work2)

    !$omp end parallel
 
    return

  end subroutine pme_recip_opt_1dalltoall_8

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_recip_opt_1dalltoall_lj_8
  !> @brief        Calculate PME reciprocal part with domain decomposition
  !                (B spline order 8)                                     
  !! @authors      JJ, NT
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_recip_opt_1dalltoall_lj_8(domain, enefunc, force, virial, &
                                           eelec, evdw)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wip),               intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw (nthread)

    ! local variables
    real(wp)                 :: elec_temp, evdw_temp, temp
    real(wp)                 :: vr(4), dv(4), grid(4), half_grid(4)
    real(wp)                 :: bscx(4), bscy(4), bscz(4)
    real(wp)                 :: bscdx(4), bscdy(4), bscdz(4)
    real(wp)                 :: vxx, vyy, vzz, vyyd, vzzd
    real(wp)                 :: efft_real, efft_imag
    real(wp)                 :: vfft_real, vfft_imag
    real(wp)                 :: tqq, tll, qtmp, ltmp
    real(wp)                 :: force_local(3), force_local_lj(3)
    real(wp)                 :: vol_factor
    real(wp)                 :: u_2(3), u_3(3), u_4(3), u_5(3)
    real(wp)                 :: u_6(3), u_7(3)
    real(wp)                 :: v_1(3), v_4(3), v_6(3), v_7(3)
    integer                  :: iatmcls
    integer                  :: i, k, icel, iproc, id
    integer                  :: ix, iv, iy, iz, ixs, iys, izs
    integer                  :: ixyz, nix, nix1, niy, niz, nizx, nizy
    integer                  :: ixx, iyy, izz, ii(4)
    integer                  :: iix(4), iiy(4), iiz(4)
    integer                  :: start(3), end(3)
    integer                  :: ncell, ncell_local, omp_get_thread_num
    integer                  :: kk
    integer                  :: iprocx, iprocy
    integer                  :: ix_start, ix_end, iz_start, iz_end

    complex(wp), allocatable :: work1(:), work2(:), work3(:), work4(:)
    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: nonb_lj6_factor(:)
    integer,         pointer :: natom(:)
    integer,         pointer :: atmcls(:,:)

    coord  => domain%coord
    charge => domain%charge
    natom  => domain%num_atom
    atmcls => domain%atom_cls_no
    nonb_lj6_factor => enefunc%nonb_lj6_factor

    ncell  = domain%num_cell_local + domain%num_cell_boundary
    ncell_local  = domain%num_cell_local

    ! Initializing the energy and force
    !
    qdf_real(1:nlocalx*nlocaly*nlocalz) = 0.0_wp
    ldf_real(1:nlocalx*nlocaly*nlocalz) = 0.0_wp

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, iv, ix, iy, iz, icel, k, ii, dv, izz, izs, iyy, iys,  &
    !$omp         iz_start, iz_end, ix_start, ix_end, ixx, ixs, vxx, vyy, vzz, &
    !$omp         iproc, tqq, vr, work1, work2, elec_temp, evdw_temp, temp,    &
    !$omp         force_local, force_local_lj, ixyz, nix, nix1, niy, niz,      &
    !$omp         bscx, bscy, bscz, bscdx, bscdy, bscdz, kk, iprocx, iprocy,   &
    !$omp         qtmp, grid, half_grid, start, end, nizx, nizy, iix, iiy,     &
    !$omp         iiz, u_2, u_3, vyyd, vzzd, iatmcls, ltmp, efft_real,         &
    !$omp         efft_imag, vfft_real, vfft_imag, u_4, u_5, u_6, u_7,         &
    !$omp         v_1, v_4, v_6, v_7, tll, work3, work4)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    ! initialization
    !
    allocate(work1(ngridmax), work2(2*ngridmax),  &
             work3(ngridmax), work4(2*ngridmax))

    start(1) = x_start
    start(2) = y_start
    start(3) = z_start
    end(1) = x_end
    end(2) = y_end
    end(3) = z_end

    ! Initialization of Q-fanction
    !
    do iz = 1, nlocalz+16
      do iy = 1, nlocaly+16
        do ix = 1, nlocalx+16
          qdf(ix,iy,iz,id+1) = 0.0_wp
          ldf(ix,iy,iz,id+1) = 0.0_wp
        end do
      end do
    end do

    ! Calculateing Q-fanction
    !
    grid(1:3) = real(ngrid(1:3),wp)
    half_grid(1:3) = real(ngrid(1:3)/2,wp)
    grid(4) = 2.0_wp
    half_grid(4) = 1.0_wp
    ngrid(4) = 100000

    do icel = id+1, ncell_local, nthread

      do i = 1, natom(icel)

        vr(1) = coord(1,i,icel) * r_scale(1)
        vr(2) = coord(2,i,icel) * r_scale(2)
        vr(3) = coord(3,i,icel) * r_scale(3)
        vr(1) = vr(1) + half_grid(1) - grid(1)*anint(vr(1)/grid(1))
        vr(2) = vr(2) + half_grid(2) - grid(2)*anint(vr(2)/grid(2))
        vr(3) = vr(3) + half_grid(3) - grid(3)*anint(vr(3)/grid(3))
        vr(1) = min(vr(1),grid(1)-0.0000001_wp)
        vr(2) = min(vr(2),grid(2)-0.0000001_wp)
        vr(3) = min(vr(3),grid(3)-0.0000001_wp)
        vi(1,i,icel) = min(int(vr(1)),ngrid(1)-1)
        vi(2,i,icel) = min(int(vr(2)),ngrid(2)-1)
        vi(3,i,icel) = min(int(vr(3)),ngrid(3)-1)
        dv(1) = vr(1) - real(vi(1,i,icel),wp)
        dv(2) = vr(2) - real(vi(2,i,icel),wp)
        dv(3) = vr(3) - real(vi(3,i,icel),wp)

        v_1(1) = 1.0_wp - dv(1)
        v_1(2) = 1.0_wp - dv(2)
        v_1(3) = 1.0_wp - dv(3)
        v_4(1) = v_1(1)*v_1(1)*v_1(1)*v_1(1)
        v_4(2) = v_1(2)*v_1(2)*v_1(2)*v_1(2)
        v_4(3) = v_1(3)*v_1(3)*v_1(3)*v_1(3)
        v_6(1) = v_4(1)*v_1(1)*v_1(1)
        v_6(2) = v_4(1)*v_1(2)*v_1(2)
        v_6(3) = v_4(1)*v_1(3)*v_1(3)
        v_7(1) = v_6(1)*v_1(1)
        v_7(2) = v_6(2)*v_1(2)
        v_7(3) = v_6(3)*v_1(3)
        u_2(1) = dv(1)*dv(1)
        u_2(2) = dv(2)*dv(2)
        u_2(3) = dv(3)*dv(3)
        u_3(1) = u_2(1)*dv(1)
        u_3(2) = u_2(2)*dv(2)
        u_3(3) = u_2(3)*dv(3)
        u_4(1) = u_2(1)*u_2(1)
        u_4(2) = u_2(2)*u_2(2)
        u_4(3) = u_2(3)*u_2(3)
        u_5(1) = u_2(1)*u_3(1)
        u_5(2) = u_2(2)*u_3(2)
        u_5(3) = u_2(3)*u_3(3)
        u_6(1) = u_3(1)*u_3(1)
        u_6(2) = u_3(2)*u_3(2)
        u_6(3) = u_3(3)*u_3(3)
        u_7(1) = u_3(1)*u_4(1)
        u_7(2) = u_3(2)*u_4(2)
        u_7(3) = u_3(3)*u_4(3)
        bsc(1,1,i,icel) = u_7(1)
        bsc(1,2,i,icel) = u_7(2)
        bsc(1,3,i,icel) = u_7(3)
        bsc(2,1,i,icel) = -7.0_wp*u_7(1)+7.0_wp*u_6(1)+21.0_wp*u_5(1)  &
                        + 35.0_wp*u_4(1)+35.0_wp*u_3(1)+21.0_wp*u_2(1) &
                        + 7.0_wp*dv(1)+1.0_wp
        bsc(2,2,i,icel) = -7.0_wp*u_7(2)+7.0_wp*u_6(2)+21.0_wp*u_5(2)  &
                        + 35.0_wp*u_4(2)+35.0_wp*u_3(2)+21.0_wp*u_2(2) &
                        + 7.0_wp*dv(2)+1.0_wp
        bsc(2,3,i,icel) = -7.0_wp*u_7(3)+7.0_wp*u_6(3)+21.0_wp*u_5(3)  &
                        + 35.0_wp*u_4(3)+35.0_wp*u_3(3)+21.0_wp*u_2(3) &
                        + 7.0_wp*dv(3)+1.0_wp
        bsc(3,1,i,icel) = 21.0_wp*u_7(1)-42.0_wp*u_6(1)-84.0_wp*u_5(1) &
                        + 280.0_wp*u_3(1)+504.0_wp*u_2(1)              &
                        + 392.0_wp*dv(1)+120.0_wp
        bsc(3,2,i,icel) = 21.0_wp*u_7(2)-42.0_wp*u_6(2)-84.0_wp*u_5(2) &
                        + 280.0_wp*u_3(2)+504.0_wp*u_2(2)              &
                        + 392.0_wp*dv(2)+120.0_wp
        bsc(3,3,i,icel) = 21.0_wp*u_7(3)-42.0_wp*u_6(3)-84.0_wp*u_5(3) &
                        + 280.0_wp*u_3(3)+504.0_wp*u_2(3)              &
                        + 392.0_wp*dv(3)+120.0_wp
        bsc(4,1,i,icel) = -35.0_wp*u_7(1)+105.0_wp*u_6(1)+105.0_wp*u_5(1) &
                        - 315.0_wp*u_4(1)-665.0_wp*u_3(1)+315.0_wp*u_2(1) &
                        + 1715.0_wp*dv(1)+1191.0_wp
        bsc(4,2,i,icel) = -35.0_wp*u_7(2)+105.0_wp*u_6(2)+105.0_wp*u_5(2) &
                        - 315.0_wp*u_4(2)-665.0_wp*u_3(2)+315.0_wp*u_2(2) &
                        + 1715.0_wp*dv(2)+1191.0_wp
        bsc(4,3,i,icel) = -35.0_wp*u_7(3)+105.0_wp*u_6(3)+105.0_wp*u_5(3) &
                        - 315.0_wp*u_4(3)-665.0_wp*u_3(3)+315.0_wp*u_2(3) &
                        + 1715.0_wp*dv(3)+1191.0_wp
        bsc(5,1,i,icel) = 35.0_wp*u_7(1)-140.0_wp*u_6(1)+560.0_wp*u_4(1) &
                        - 1680.0_wp*u_2(1)+2416.0_wp
        bsc(5,2,i,icel) = 35.0_wp*u_7(2)-140.0_wp*u_6(2)+560.0_wp*u_4(2) &
                        - 1680.0_wp*u_2(2)+2416.0_wp
        bsc(5,3,i,icel) = 35.0_wp*u_7(3)-140.0_wp*u_6(3)+560.0_wp*u_4(3) &
                        - 1680.0_wp*u_2(3)+2416.0_wp
        bsc(6,1,i,icel) = -21.0_wp*u_7(1)+105.0_wp*u_6(1)-105.0_wp*u_5(1) &
                        - 315.0_wp*u_4(1)+665.0_wp*u_3(1)+315.0_wp*u_2(1) &
                        - 1715.0_wp*dv(1)+1191.0_wp
        bsc(6,2,i,icel) = -21.0_wp*u_7(2)+105.0_wp*u_6(2)-105.0_wp*u_5(2) &
                        - 315.0_wp*u_4(2)+665.0_wp*u_3(2)+315.0_wp*u_2(2) &
                        - 1715.0_wp*dv(2)+1191.0_wp
        bsc(6,3,i,icel) = -21.0_wp*u_7(3)+105.0_wp*u_6(3)-105.0_wp*u_5(3) &
                        - 315.0_wp*u_4(3)+665.0_wp*u_3(3)+315.0_wp*u_2(3) &
                        - 1715.0_wp*dv(3)+1191.0_wp
        bsc(7,1,i,icel) = 7.0_wp*u_7(1)-42.0_wp*u_6(1)+84.0_wp*u_5(1) &
                        - 280.0_wp*u_3(1)+504.0_wp*u_2(1)             &
                        - 392.0_wp*dv(1)+120.0_wp
        bsc(7,2,i,icel) = 7.0_wp*u_7(2)-42.0_wp*u_6(2)+84.0_wp*u_5(2) &
                        - 280.0_wp*u_3(2)+504.0_wp*u_2(2)             &
                        - 392.0_wp*dv(2)+120.0_wp
        bsc(7,3,i,icel) = 7.0_wp*u_7(3)-42.0_wp*u_6(3)+84.0_wp*u_5(3) &
                        - 280.0_wp*u_3(3)+504.0_wp*u_2(3)             &
                        - 392.0_wp*dv(3)+120.0_wp
        bsc(8,1,i,icel) = v_7(1)
        bsc(8,2,i,icel) = v_7(2)
        bsc(8,3,i,icel) = v_7(3)
        bscd(1,1,i,icel) = u_6(1)
        bscd(1,2,i,icel) = u_6(2)
        bscd(1,3,i,icel) = u_6(3)
        bscd(2,1,i,icel) = -7.0_wp*u_6(1)+6.0_wp*u_5(1)+15.0_wp*u_4(1) &
                         + 20.0_wp*u_3(1)+15.0_wp*u_2(1)+6.0_wp*dv(1)  &
                         + 1.0_wp
        bscd(2,2,i,icel) = -7.0_wp*u_6(2)+6.0_wp*u_5(2)+15.0_wp*u_4(2) &
                         + 20.0_wp*u_3(2)+15.0_wp*u_2(2)+6.0_wp*dv(2)  &
                         + 1.0_wp
        bscd(2,3,i,icel) = -7.0_wp*u_6(3)+6.0_wp*u_5(3)+15.0_wp*u_4(3) &
                         + 20.0_wp*u_3(3)+15.0_wp*u_2(3)+6.0_wp*dv(3)  &
                         + 1.0_wp
        bscd(3,1,i,icel) = 21.0_wp*u_6(1)-36.0_wp*u_5(1)-60.0_wp*u_4(1) &
                         + 120.0_wp*u_2(1)+144.0_wp*dv(1)+56.0_wp
        bscd(3,2,i,icel) = 21.0_wp*u_6(2)-36.0_wp*u_5(2)-60.0_wp*u_4(2) &
                         + 120.0_wp*u_2(2)+144.0_wp*dv(2)+56.0_wp
        bscd(3,3,i,icel) = 21.0_wp*u_6(3)-36.0_wp*u_5(3)-60.0_wp*u_4(3) &
                         + 120.0_wp*u_2(3)+144.0_wp*dv(3)+56.0_wp
        bscd(4,1,i,icel) = -35.0_wp*u_6(1)+90.0_wp*u_5(1)+75.0_wp*u_4(1) &
                         - 180.0_wp*u_3(1)-285.0_wp*u_2(1)+90.0_wp*dv(1) &
                         + 245.0_wp
        bscd(4,2,i,icel) = -35.0_wp*u_6(2)+90.0_wp*u_5(2)+75.0_wp*u_4(2) &
                         - 180.0_wp*u_3(2)-285.0_wp*u_2(2)+90.0_wp*dv(2) &
                         + 245.0_wp
        bscd(4,3,i,icel) = -35.0_wp*u_6(3)+90.0_wp*u_5(3)+75.0_wp*u_4(3) &
                         - 180.0_wp*u_3(3)-285.0_wp*u_2(3)+90.0_wp*dv(3) &
                         + 245.0_wp
        bscd(5,1,i,icel) = 35.0_wp*u_6(1)-120.0_wp*u_5(1)+320.0_wp*u_3(1) &
                         - 480.0_wp*dv(1)
        bscd(5,2,i,icel) = 35.0_wp*u_6(2)-120.0_wp*u_5(2)+320.0_wp*u_3(2) &
                         - 480.0_wp*dv(2)
        bscd(5,3,i,icel) = 35.0_wp*u_6(3)-120.0_wp*u_5(3)+320.0_wp*u_3(3) &
                         - 480.0_wp*dv(3)
        bscd(6,1,i,icel) = -21.0_wp*u_6(1)+90.0_wp*u_5(1)-75.0_wp*u_4(1) &
                         - 180.0_wp*u_3(1)+285.0_wp*u_2(1)+90.0_wp*dv(1) &
                         -245.0_wp
        bscd(6,2,i,icel) = -21.0_wp*u_6(2)+90.0_wp*u_5(2)-75.0_wp*u_4(2) &
                         - 180.0_wp*u_3(2)+285.0_wp*u_2(2)+90.0_wp*dv(2) &
                         -245.0_wp
        bscd(6,3,i,icel) = -21.0_wp*u_6(3)+90.0_wp*u_5(3)-75.0_wp*u_4(3) &
                         - 180.0_wp*u_3(3)+285.0_wp*u_2(3)+90.0_wp*dv(3) &
                         -245.0_wp
        bscd(7,1,i,icel) = 7.0_wp*u_6(1)-36.0_wp*u_5(1)+60.0_wp*u_4(1) &
                         -120.0_wp*u_2(1)+144.0_wp*dv(1)-56.0_wp
        bscd(7,2,i,icel) = 7.0_wp*u_6(2)-36.0_wp*u_5(2)+60.0_wp*u_4(2) &
                         -120.0_wp*u_2(2)+144.0_wp*dv(2)-56.0_wp
        bscd(7,3,i,icel) = 7.0_wp*u_6(3)-36.0_wp*u_5(3)+60.0_wp*u_4(3) &
                         -120.0_wp*u_2(3)+144.0_wp*dv(3)-56.0_wp
        bscd(8,1,i,icel) = -v_6(1)
        bscd(8,2,i,icel) = -v_6(2)
        bscd(8,3,i,icel) = -v_6(3)

      end do

    end do

    do icel = id+1, ncell_local, nthread
      do i = 1, natom(icel)

        qtmp = charge(i,icel)*bs_fact3
        iatmcls = atmcls(i,icel)
        ltmp = nonb_lj6_factor(iatmcls)*bs_fact3
        ii(1:3) = vi(1:3,i,icel)

        do ixyz = 1, 8
          iz = ii(3) - ixyz + 2
          if (iz <= 0) iz = iz + ngrid(3)
          vii(ixyz,3,i,icel) = grid_g2lz(iz)
          iy = ii(2) - ixyz + 2
          if (iy <= 0) iy = iy + ngrid(2)
          vii(ixyz,2,i,icel) = grid_g2ly(iy)
          ix = ii(1) - ixyz + 2
          if (ix <= 0) ix = ix + ngrid(1)
          vii(ixyz,1,i,icel) = grid_g2lx(ix)
        end do

        do izz = 1, 8
          izs = vii(izz,3,i,icel)
          vzz = bsc(izz,3,i,icel)
          do iyy = 1, 8
            iys = vii(iyy,2,i,icel)
            vyy = bsc(iyy,2,i,icel)
!ocl norecurrence(qdf)
            do ixx = 1, 8
              ixs = vii(ixx,1,i,icel)
              qdf(ixs,iys,izs,id+1) = qdf(ixs,iys,izs,id+1)            &
                                    + bsc(ixx,1,i,icel)*vyy*vzz*qtmp
              ldf(ixs,iys,izs,id+1) = ldf(ixs,iys,izs,id+1)            &
                                    + bsc(ixx,1,i,icel)*vyy*vzz*ltmp
            end do
          end do
        end do
      end do
    end do

    !$omp barrier
    do iproc = 2, nthread
      do iz = id+1, nlocalz+16,nthread
        do iy = 1, nlocaly+16
          do ix = 1, nlocalx+16
            qdf(ix,iy,iz,1) = qdf(ix,iy,iz,1) + qdf(ix,iy,iz,iproc)
            ldf(ix,iy,iz,1) = ldf(ix,iy,iz,1) + ldf(ix,iy,iz,iproc)
          end do
        end do
      end do
    end do

    !$omp barrier
    !$omp master
    call communicate_pme_pre_lj
    !$omp end master
    !$omp barrier

    do iz = id+1, nlocalz, nthread
      izs = iz + 8
      do iy = 1, nlocaly
        k = (iy-1)*nlocalx + (iz-1)*nlocalx*nlocaly
        iys = iy + 8
        do ix = 1, nlocalx
          ixs = ix + 8
          qdf_real(k+ix) = qdf_real(k+ix) + qdf(ixs,iys,izs,1)
          ldf_real(k+ix) = ldf_real(k+ix) + ldf(ixs,iys,izs,1)
        end do
      end do
    end do

    !$omp barrier
    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdf_real, nlocalx*nlocaly*nlocalz/nprocx, mpi_wp_real, &
                      qdf_work, nlocalx*nlocaly*nlocalz/nprocx, mpi_wp_real, &
                      grid_commx, ierror)
    call mpi_alltoall(ldf_real, nlocalx*nlocaly*nlocalz/nprocx, mpi_wp_real, &
                      ldf_work, nlocalx*nlocaly*nlocalz/nprocx, mpi_wp_real, &
                      grid_commx, ierror)
#else
    qdf_work(1:nlocalx*nlocaly*nlocalz,1) = qdf_real(1:nlocalx*nlocaly*nlocalz)
    ldf_work(1:nlocalx*nlocaly*nlocalz,1) = ldf_real(1:nlocalx*nlocaly*nlocalz)
#endif
    !$omp end master
    !$omp barrier

    nizx = nlocalz / nprocx
    nizy = nlocalz / nprocy
    call fft3d_1d_alltoall_lj(qdf_work, ftqdf, ftqdf, ftqdf, ftqdf_work,       &
                              ftqdf_work, ftqdf_work, ftqdf_work, ftqdf_work2, &
                              ldf_work, ftldf, ftldf, ftldf, ftldf_work,       &
                              ftldf_work, ftldf_work, ftldf_work, ftldf_work2, &
                              work1, work2, work3, work4,                      &
                              nlocalx, nlocaly, nlocalz, nizx, nizy,           &
                              nlocalx1, x_local1, ngrid(1), ngrid(2),          &
                              ngrid(3), nprocx, nprocy, nprocz, id, nthread,   &
                              grid_commx, grid_commy, grid_commz)


    ! Energy calculation
    !
    iproc = my_z_rank + 1
    niy = nlocaly / nprocz

    do ix = 1, nlocalx/2
      do iy = id+1, niy, nthread
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp
        vxx = 0.0_wp
        vyy = 0.0_wp
        vzz = 0.0_wp

        do iz = 1, ngrid(3)
          efft_real = real(ftqdf_work2(iz,ix,iy),wp)
          efft_imag = imag(ftqdf_work2(iz,ix,iy))
          vfft_real = real(ftldf_work2(iz,ix,iy),wp)
          vfft_imag = imag(ftldf_work2(iz,ix,iy))
          tqq = efft_real*efft_real + efft_imag*efft_imag
          tqq = tqq * theta(iz,iy,ix) * vol_fact4
          elec_temp = elec_temp + tqq
          tll = vfft_real*vfft_real + vfft_imag*vfft_imag
          tll = tll * vol_fact4_lj
          evdw_temp = evdw_temp - tll*theta_lj(iz,iy,ix)

          ! virial
          !
          vxx = vxx + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gx(ix) * gx(ix))
          vyy = vyy + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gy(iy) * gy(iy))
          vzz = vzz + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gz(iz) * gz(iz))
          vxx = vxx  &
              - tll*(theta_lj(iz,iy,ix)+vir_fact_lj(iz,iy,ix)*gx(ix)*gx(ix))
          vyy = vyy  &
              - tll*(theta_lj(iz,iy,ix)+vir_fact_lj(iz,iy,ix)*gy(iy)*gy(iy))
          vzz = vzz  &
              - tll*(theta_lj(iz,iy,ix)+vir_fact_lj(iz,iy,ix)*gz(iz)*gz(iz))
        end do

        eelec(id+1) = eelec(id+1) + elec_temp
        evdw (id+1) = evdw (id+1) + evdw_temp
        virial(1,1,id+1) = virial(1,1,id+1) + vxx
        virial(2,2,id+1) = virial(2,2,id+1) + vyy
        virial(3,3,id+1) = virial(3,3,id+1) + vzz
      end do
    end do

    if (my_x_rank == 0) then

      do iy = id+1, niy, nthread

        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp
        vxx = 0.0_wp
        vyy = 0.0_wp
        vzz = 0.0_wp

        do iz = 1, ngrid(3)

          efft_real = real(ftqdf_work2(iz,1,iy),wp)
          efft_imag = imag(ftqdf_work2(iz,1,iy))
          vfft_real = real(ftldf_work2(iz,1,iy),wp)
          vfft_imag = imag(ftldf_work2(iz,1,iy))
          tqq = efft_real*efft_real + efft_imag*efft_imag
          tqq = tqq * theta(iz,iy,1) * vol_fact2
          elec_temp = elec_temp - tqq
          tll = vfft_real*vfft_real + vfft_imag*vfft_imag
          tll = tll * vol_fact2_lj
          evdw_temp = evdw_temp + tll*theta_lj(iz,iy,1)

          vxx = vxx - tqq * (1.0_wp + vir_fact(iz,iy,1) * gx(1) * gx(1))
          vyy = vyy - tqq * (1.0_wp + vir_fact(iz,iy,1) * gy(iy) * gy(iy))
          vzz = vzz - tqq * (1.0_wp + vir_fact(iz,iy,1) * gz(iz) * gz(iz))
          vxx = vxx  &
              + tll*(theta_lj(iz,iy,1)+vir_fact_lj(iz,iy,1)*gx(1 )*gx(1 ))
          vyy = vyy  &
              + tll*(theta_lj(iz,iy,1)+vir_fact_lj(iz,iy,1)*gy(iy)*gy(iy))
          vzz = vzz  &
              + tll*(theta_lj(iz,iy,1)+vir_fact_lj(iz,iy,1)*gz(iz)*gz(iz))

        end do
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw (id+1) = evdw (id+1) + evdw_temp
        virial(1,1,id+1) = virial(1,1,id+1) + vxx
        virial(2,2,id+1) = virial(2,2,id+1) + vyy
        virial(3,3,id+1) = virial(3,3,id+1) + vzz
      end do

      ix = nlocalx/2 + 1

      do iy = id+1, niy, nthread

        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp
        vxx = 0.0_wp
        vyy = 0.0_wp
        vzz = 0.0_wp

        do iz = 1, ngrid(3)
  
          efft_real = real(ftqdf_work2(iz,ix,iy),wp)
          efft_imag = imag(ftqdf_work2(iz,ix,iy))
          vfft_real = real(ftldf_work2(iz,ix,iy),wp)
          vfft_imag = imag(ftldf_work2(iz,ix,iy))
          tqq = efft_real*efft_real + efft_imag*efft_imag
          tqq = tqq * theta(iz,iy,ix) * vol_fact2
          elec_temp = elec_temp + tqq
          tll = vfft_real*vfft_real + vfft_imag*vfft_imag
          tll = tll * vol_fact2_lj
          evdw_temp = evdw_temp - tll*theta_lj(iz,iy,ix)

          vxx = vxx + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gx(ix) * gx(ix))
          vyy = vyy + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gy(iy) * gy(iy))
          vzz = vzz + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gz(iz) * gz(iz))
          vxx = vxx &
              - tll*(theta_lj(iz,iy,ix)+vir_fact_lj(iz,iy,ix)*gx(ix)*gx(ix))
          vyy = vyy &
              - tll*(theta_lj(iz,iy,ix)+vir_fact_lj(iz,iy,ix)*gy(iy)*gy(iy))
          vzz = vzz &
              - tll*(theta_lj(iz,iy,ix)+vir_fact_lj(iz,iy,ix)*gz(iz)*gz(iz))

        end do

        eelec(id+1) = eelec(id+1) + elec_temp
        evdw (id+1) = evdw (id+1) + evdw_temp
        virial(1,1,id+1) = virial(1,1,id+1) + vxx
        virial(2,2,id+1) = virial(2,2,id+1) + vyy
        virial(3,3,id+1) = virial(3,3,id+1) + vzz

      end do
    end if

    ! add self dispersion term
    !
    !$omp master
    if (my_city_rank == 0) then
      temp    = PI*sqrt(PI)*alpha3*inv_volume/(6.0_wp)*enefunc%pme_dispersion_self1
      evdw(1) = evdw(1) - temp
      evdw(1) = evdw(1) + alpha6/12.0_wp*enefunc%pme_dispersion_self2
      virial(1,1,1) = virial(1,1,1) - temp
      virial(2,2,1) = virial(2,2,1) - temp
      virial(3,3,1) = virial(3,3,1) - temp
    end if
    !$omp end master

    ! F^-1[Th]*F^-1[Q] (=X)
    !
    !$omp barrier
    do ix = 1, x_local1
      do iy = id+1, niy, nthread
        k = (iy-1)*ngrid(3) + (ix-1)*ngrid(3)*niy
        do iz = 1, ngrid(3)
          ftqdf_work(k+iz) = ftqdf_work2(iz,ix,iy) &
                            *cmplx(theta(iz,iy,ix),0.0_wp)
          ftldf_work(k+iz) = ftldf_work2(iz,ix,iy) &
                            *cmplx(theta_lj(iz,iy,ix),0.0_wp)
        end do
      end do
    end do

    !$omp barrier
    call bfft3d_1d_alltoall_lj(qdf_work, qdf_real, ftqdf, ftqdf, ftqdf_work, &
                               ftqdf_work, ftqdf_work2, ftqdf_work,          &
                               ldf_work, ldf_real, ftldf, ftldf, ftldf_work, &
                               ftldf_work, ftldf_work2, ftldf_work,          &
                               work1, work2, work3, work4,                   &
                               nlocalx, nlocaly, nlocalz, nlocalx1,          &
                               x_local1, ngrid(1), ngrid(2), ngrid(3),       &
                               nprocx, nprocy, nprocz, id, nthread,          &
                               my_x_rank, grid_commx, grid_commy,            &
                               grid_commz, niy, nizy, nizx)

    ! X is saved on qdf
    !
    do iz = id+1, nlocalz, nthread
      izs = iz + 8
      do iy = 1, nlocaly
        iys = iy + 8
        k = (iy-1)*nlocalx + (iz-1)*nlocalx*nlocaly
        do ix = 1, nlocalx
          ixs = ix + 8
          qdf(ixs,iys,izs,1) = qdf_real(k+ix)
          ldf(ixs,iys,izs,1) = ldf_real(k+ix)
        end do
      end do
    end do

    !$omp barrier
    !$omp master
    call communicate_pme_post_lj
    !$omp end master
    !$omp barrier

    do icel = id+1, ncell_local, nthread
      do i = 1, natom(icel)

        force_local(1:3) = 0.0_wp
        force_local_lj(1:3) = 0.0_wp
        qtmp = charge(i,icel)*vol_fact4*bs_fact3d
        iatmcls = atmcls(i,icel)
        ltmp = nonb_lj6_factor(iatmcls)*vol_fact4_lj*bs_fact3d

        do izz = 1, 8
          izs = vii(izz,3,i,icel)
          vzz = bsc(izz,3,i,icel)
          vzzd = bscd(izz,3,i,icel)
          do iyy = 1, 8
            iys = vii(iyy,2,i,icel)
            vyy = bsc(iyy,2,i,icel)
            vyyd = bscd(iyy,2,i,icel)
            do ixx = 1, 8
              ixs = vii(ixx,1,i,icel)
              force_local(1)    = force_local(1) &
                                + bscd(ixx,1,i,icel)*vyy*vzz*qdf(ixs,iys,izs,1)
              force_local(2)    = force_local(2)     &
                                + bsc(ixx,1,i,icel)*vyyd*vzz*qdf(ixs,iys,izs,1)
              force_local(3)    = force_local(3)     &
                                + bsc(ixx,1,i,icel)*vyy*vzzd*qdf(ixs,iys,izs,1)
              force_local_lj(1) = force_local_lj(1) &
                                + bscd(ixx,1,i,icel)*vyy*vzz*ldf(ixs,iys,izs,1)
              force_local_lj(2) = force_local_lj(2)     &
                                + bsc(ixx,1,i,icel)*vyyd*vzz*ldf(ixs,iys,izs,1)
              force_local_lj(3) = force_local_lj(3)     &
                                + bsc(ixx,1,i,icel)*vyy*vzzd*ldf(ixs,iys,izs,1)
            end do
          end do
        end do

        force(1:3,i,icel) = - force_local(1:3)*qtmp*r_scale(1:3)    &
                            + force_local_lj(1:3)*ltmp*r_scale(1:3)

      end do
    end do

    deallocate(work1, work2, work3, work4)

    !$omp end parallel
 
    return

  end subroutine pme_recip_opt_1dalltoall_lj_8

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    communicate_pme_pre
  !> @brief        communicate boundary grid data before FFT
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine communicate_pme_pre

    ! local variable
    integer                  :: k, ix, iy, iz
    integer                  :: irequest1, irequest2, irequest3, irequest4
    integer                  :: upper_rank, lower_rank
#ifdef HAVE_MPI_GENESIS
    integer                  :: istatus(mpi_status_size)
#endif


    ! z direction
    !
    k = 0
    do iy = 1, nlocaly+2*n_bs
      do ix = 1, nlocalx+2*n_bs
        do iz = 1, n_bs
          buf_send(iz+k,1) = qdf(ix,iy,iz,1)  ! send to lower
          buf_send(iz+k,2) = qdf(ix,iy,n_bs+nlocalz+iz,1) !send to upper
        end do
        k = k + n_bs
      end do
    end do

#ifdef HAVE_MPI_GENESIS
    upper_rank = mod(my_z_rank+1,nprocz)
    lower_rank = mod(my_z_rank-1+nprocz,nprocz)
    ! send to lower, receive from upper
    call mpi_irecv(buf_recv(1,1), k, mpi_wp_real, upper_rank, &
                   upper_rank+nprocz+my_z_rank,               &
                   grid_commz, irequest1, ierror)
    call mpi_isend(buf_send(1,1), k, mpi_wp_real, lower_rank, &
                   my_z_rank+nprocz+lower_rank,               &
                   grid_commz, irequest2, ierror)
    ! send to upper, receive from lower
    call mpi_irecv(buf_recv(1,2), k, mpi_wp_real, lower_rank, &
                   my_z_rank+2*nprocz+lower_rank,               &
                   grid_commz, irequest3, ierror)
    call mpi_isend(buf_send(1,2), k, mpi_wp_real, upper_rank, &
                   upper_rank+2*nprocz+my_z_rank,               &
                   grid_commz, irequest4, ierror)
    call mpi_wait(irequest1, istatus, ierror)
    call mpi_wait(irequest2, istatus, ierror)
    call mpi_wait(irequest3, istatus, ierror)
    call mpi_wait(irequest4, istatus, ierror)
#else
    buf_recv(1:k,1) = buf_send(1:k,1)
    buf_recv(1:k,2) = buf_send(1:k,2)
#endif

    k = 0
    do iy = 1, nlocaly+2*n_bs
      do ix = 1, nlocalx+2*n_bs
        do iz = 1, n_bs
          qdf(ix,iy,nlocalz+iz,1) = &     !receive from upper
            qdf(ix,iy,nlocalz+iz,1) + buf_recv(iz+k,1)
          qdf(ix,iy,n_bs+iz,1)       = &     !receive from lower
           qdf(ix,iy,n_bs+iz,1) + buf_recv(iz+k,2)
        end do
        k = k + n_bs
      end do
    end do

    ! y direction
    !
    k = 0
    do iz = n_bs+1, n_bs+nlocalz
      do ix = 1, nlocalx+2*n_bs
        do iy = 1, n_bs
          buf_send(iy+k,1) = qdf(ix,iy,iz,1)
          buf_send(iy+k,2) = qdf(ix,n_bs+nlocaly+iy,iz,1)
        end do
        k = k + n_bs
      end do
    end do

#ifdef HAVE_MPI_GENESIS
    upper_rank = mod(my_y_rank+1,nprocy)
    lower_rank = mod(my_y_rank-1+nprocy,nprocy)
    ! send to lower, receive from upper
    call mpi_irecv(buf_recv(1,1), k, mpi_wp_real, upper_rank, &
                   upper_rank+nprocy+my_y_rank,               &
                   grid_commy, irequest1, ierror)
    call mpi_isend(buf_send(1,1), k, mpi_wp_real, lower_rank, &
                   my_y_rank+nprocy+lower_rank,               &
                   grid_commy, irequest2, ierror)
    ! send to upper, receive from lower
    call mpi_irecv(buf_recv(1,2), k, mpi_wp_real, lower_rank, &
                   my_y_rank+2*nprocy+lower_rank,             &
                   grid_commy, irequest3, ierror)
    call mpi_isend(buf_send(1,2), k, mpi_wp_real, upper_rank, &
                   upper_rank+2*nprocy+my_y_rank,             &
                   grid_commy, irequest4, ierror)
    call mpi_wait(irequest1, istatus, ierror)
    call mpi_wait(irequest2, istatus, ierror)
    call mpi_wait(irequest3, istatus, ierror)
    call mpi_wait(irequest4, istatus, ierror)
#else
    buf_recv(1:k,1) = buf_send(1:k,1)
    buf_recv(1:k,2) = buf_send(1:k,2)
#endif

    k = 0
    do iz = n_bs+1, n_bs+nlocalz
      do ix = 1, nlocalx+2*n_bs
        do iy = 1, n_bs
          qdf(ix,nlocaly+iy,iz,1) = &
            qdf(ix,nlocaly+iy,iz,1) + buf_recv(iy+k,1)
          qdf(ix,n_bs+iy,iz,1)       = &
            qdf(ix,n_bs+iy,iz,1) + buf_recv(iy+k,2)
        end do
        k = k + n_bs
      end do
    end do

    ! x direction
    !
    k = 0
    do iz = n_bs+1, n_bs+nlocalz
      do iy = n_bs+1, n_bs+nlocaly
        do ix = 1, n_bs
          buf_send(ix+k,1) = qdf(ix,iy,iz,1)
          buf_send(ix+k,2) = qdf(n_bs+nlocalx+ix,iy,iz,1)
        end do
        k = k + n_bs
      end do
    end do

#ifdef HAVE_MPI_GENESIS
    upper_rank = mod(my_x_rank+1,nprocx)
    lower_rank = mod(my_x_rank-1+nprocx,nprocx)
    ! send to lower, receive from upper
    call mpi_irecv(buf_recv(1,1), k, mpi_wp_real, upper_rank, &
                   upper_rank+nprocx+my_x_rank,               &
                   grid_commx, irequest1, ierror)
    call mpi_isend(buf_send(1,1), k, mpi_wp_real, lower_rank, &
                   my_x_rank+nprocx+lower_rank,               &
                   grid_commx, irequest2, ierror)
    ! send to upper, receive from lower
    call mpi_irecv(buf_recv(1,2), k, mpi_wp_real, lower_rank, &
                   my_x_rank+2*nprocx+lower_rank,             &
                   grid_commx, irequest3, ierror)
    call mpi_isend(buf_send(1,2), k, mpi_wp_real, upper_rank, &
                   upper_rank+2*nprocx+my_x_rank,             &
                   grid_commx, irequest4, ierror)
    call mpi_wait(irequest1, istatus, ierror)
    call mpi_wait(irequest2, istatus, ierror)
    call mpi_wait(irequest3, istatus, ierror)
    call mpi_wait(irequest4, istatus, ierror)
#else
    buf_recv(1:k,1) = buf_send(1:k,1)
    buf_recv(1:k,2) = buf_send(1:k,2)
#endif

    k = 0
    do iz = n_bs+1, n_bs+nlocalz
      do iy = n_bs+1, n_bs+nlocaly
        do ix = 1, n_bs
          qdf(nlocalx+ix,iy,iz,1) = &
            qdf(nlocalx+ix,iy,iz,1) + buf_recv(ix+k,1)
          qdf(n_bs+ix,iy,iz,1) =       &
            qdf(n_bs+ix,iy,iz,1) + buf_recv(ix+k,2)
        end do
        k = k + n_bs
      end do
    end do

    return

  end subroutine communicate_pme_pre

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    communicate_pme_post
  !> @brief        communicate boundary grid data after bFFT
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine communicate_pme_post

    ! local variable
    integer                  :: k, ix, iy, iz
    integer                  :: irequest1, irequest2, irequest3, irequest4
    integer                  :: upper_rank, lower_rank
#ifdef HAVE_MPI_GENESIS
    integer                  :: istatus(mpi_status_size)
#endif


    ! x direction
    !
    k = 0
    do iz = n_bs+1, n_bs+nlocalz
      do iy = n_bs+1, n_bs+nlocaly
        do ix = 1, n_bs
          buf_send(ix+k,1) = qdf(n_bs+ix,iy,iz,1)     ! send to lower
          buf_send(ix+k,2) = qdf(nlocalx+ix,iy,iz,1)  ! send to upper
        end do
        k = k + n_bs
      end do
    end do

#ifdef HAVE_MPI_GENESIS
    upper_rank = mod(my_x_rank+1,nprocx)
    lower_rank = mod(my_x_rank-1+nprocx,nprocx)
    ! send to lower, receive from upper
    call mpi_irecv(buf_recv(1,1), k, mpi_wp_real, upper_rank, &
                   upper_rank+nprocx+my_x_rank,               &
                   grid_commx, irequest1, ierror)
    call mpi_isend(buf_send(1,1), k, mpi_wp_real, lower_rank, &
                   my_x_rank+nprocx+lower_rank,               &
                   grid_commx, irequest2, ierror)
    ! send to upper, receive from lower
    call mpi_irecv(buf_recv(1,2), k, mpi_wp_real, lower_rank, &
                   my_x_rank+2*nprocx+lower_rank,             &
                   grid_commx, irequest3, ierror)
    call mpi_isend(buf_send(1,2), k, mpi_wp_real, upper_rank, &
                   upper_rank+2*nprocx+my_x_rank,             &
                   grid_commx, irequest4, ierror)
    call mpi_wait(irequest1, istatus, ierror)
    call mpi_wait(irequest2, istatus, ierror)
    call mpi_wait(irequest3, istatus, ierror)
    call mpi_wait(irequest4, istatus, ierror)
#else
    buf_recv(1:k,1) = buf_send(1:k,1)
    buf_recv(1:k,2) = buf_send(1:k,2)
#endif

    k = 0
    do iz = n_bs+1, n_bs+nlocalz
      do iy = n_bs+1, n_bs+nlocaly
        do ix = 1, n_bs
          qdf(n_bs+nlocalx+ix,iy,iz,1) = buf_recv(ix+k,1) ! receive from upper
          qdf(ix,iy,iz,1)  = buf_recv(ix+k,2) ! receive from lower
        end do
        k = k + n_bs
      end do
    end do

    ! y direction
    !
    k = 0
    do iz = n_bs+1, n_bs+nlocalz
      do ix = 1, nlocalx+2*n_bs
        do iy = 1, n_bs
          buf_send(iy+k,1) = qdf(ix,n_bs+iy,iz,1)      ! send to lower
          buf_send(iy+k,2) = qdf(ix,nlocaly+iy,iz,1)   ! send to upper
        end do
        k = k + n_bs
      end do
    end do

#ifdef HAVE_MPI_GENESIS
    upper_rank = mod(my_y_rank+1,nprocy)
    lower_rank = mod(my_y_rank-1+nprocy,nprocy)
    ! send to lower, receive from upper
    call mpi_irecv(buf_recv(1,1), k, mpi_wp_real, upper_rank, &
                   upper_rank+nprocy+my_y_rank,               &
                   grid_commy, irequest1, ierror)
    call mpi_isend(buf_send(1,1), k, mpi_wp_real, lower_rank, &
                   my_y_rank+nprocy+lower_rank,               &
                   grid_commy, irequest2, ierror)
    ! send to upper, receive from lower
    call mpi_irecv(buf_recv(1,2), k, mpi_wp_real, lower_rank, &
                   my_y_rank+2*nprocy+lower_rank,             &
                   grid_commy, irequest3, ierror)
    call mpi_isend(buf_send(1,2), k, mpi_wp_real, upper_rank, &
                   upper_rank+2*nprocy+my_y_rank,             &
                   grid_commy, irequest4, ierror)
    call mpi_wait(irequest1, istatus, ierror)
    call mpi_wait(irequest2, istatus, ierror)
    call mpi_wait(irequest3, istatus, ierror)
    call mpi_wait(irequest4, istatus, ierror)
#else
    buf_recv(1:k,1) = buf_send(1:k,1)
    buf_recv(1:k,2) = buf_send(1:k,2)
#endif

    k = 0
    do iz = n_bs+1, n_bs+nlocalz
      do ix = 1, nlocalx+2*n_bs
        do iy = 1, n_bs
          qdf(ix,n_bs+nlocaly+iy,iz,1) = buf_recv(iy+k,1)
          qdf(ix,iy,iz,1)    = buf_recv(iy+k,2)
        end do
        k = k + n_bs
      end do
    end do

    ! z direction
    !
    k = 0
    do iy = 1, nlocaly+2*n_bs
      do ix = 1, nlocalx+2*n_bs
        do iz = 1, n_bs
          buf_send(iz+k,1) = qdf(ix,iy,n_bs+iz,1)  ! send to lower
          buf_send(iz+k,2) = qdf(ix,iy,nlocalz+iz,1) !send to upper
        end do
        k = k + n_bs
      end do
    end do

#ifdef HAVE_MPI_GENESIS
    upper_rank = mod(my_z_rank+1,nprocz)
    lower_rank = mod(my_z_rank-1+nprocz,nprocz)
    ! send to lower, receive from upper
    call mpi_irecv(buf_recv(1,1), k, mpi_wp_real, upper_rank, &
                   upper_rank+nprocz+my_z_rank,               &
                   grid_commz, irequest1, ierror)
    call mpi_isend(buf_send(1,1), k, mpi_wp_real, lower_rank, &
                   my_z_rank+nprocz+lower_rank,               &
                   grid_commz, irequest2, ierror)
    ! send to upper, receive from lower
    call mpi_irecv(buf_recv(1,2), k, mpi_wp_real, lower_rank, &
                   my_z_rank+2*nprocz+lower_rank,             &
                   grid_commz, irequest3, ierror)
    call mpi_isend(buf_send(1,2), k, mpi_wp_real, upper_rank, &
                   upper_rank+2*nprocz+my_z_rank,             &
                   grid_commz, irequest4, ierror)
    call mpi_wait(irequest1, istatus, ierror)
    call mpi_wait(irequest2, istatus, ierror)
    call mpi_wait(irequest3, istatus, ierror)
    call mpi_wait(irequest4, istatus, ierror)
#else
    buf_recv(1:k,1) = buf_send(1:k,1)
    buf_recv(1:k,2) = buf_send(1:k,2)
#endif

    k = 0
    do iy = 1, nlocaly+2*n_bs
      do ix = 1, nlocalx+2*n_bs
        do iz = 1, n_bs
          qdf(ix,iy,n_bs+nlocalz+iz,1) = buf_recv(iz+k,1)
          qdf(ix,iy,iz,1) = buf_recv(iz+k,2)
        end do
        k = k + n_bs
      end do
    end do

    return

  end subroutine communicate_pme_post

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    communicate_pme_pre_lj
  !> @brief        communicate boundary grid data before FFT
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine communicate_pme_pre_lj

    ! local variable
    integer                  :: k, ix, iy, iz
    integer                  :: irequest1, irequest2, irequest3, irequest4
    integer                  :: upper_rank, lower_rank
#ifdef HAVE_MPI_GENESIS
    integer                  :: istatus(mpi_status_size)
#endif


    ! z direction
    !
    k = 0
    do iy = 1, nlocaly+2*n_bs
      do ix = 1, nlocalx+2*n_bs
        do iz = 1, n_bs
          buf_send(iz+k,1) = qdf(ix,iy,iz,1)  ! send to lower
          buf_send(iz+k,2) = qdf(ix,iy,n_bs+nlocalz+iz,1) !send to upper
        end do
        k = k + n_bs
      end do
    end do
    do iy = 1, nlocaly+2*n_bs
      do ix = 1, nlocalx+2*n_bs
        do iz = 1, n_bs
          buf_send(iz+k,1) = ldf(ix,iy,iz,1)  ! send to lower
          buf_send(iz+k,2) = ldf(ix,iy,n_bs+nlocalz+iz,1) !send to upper
        end do
        k = k + n_bs
      end do
    end do

#ifdef HAVE_MPI_GENESIS
    upper_rank = mod(my_z_rank+1,nprocz)
    lower_rank = mod(my_z_rank-1+nprocz,nprocz)
    ! send to lower, receive from upper
    call mpi_irecv(buf_recv(1,1), k, mpi_wp_real, upper_rank, &
                   upper_rank+nprocz+my_z_rank,               &
                   grid_commz, irequest1, ierror)
    call mpi_isend(buf_send(1,1), k, mpi_wp_real, lower_rank, &
                   my_z_rank+nprocz+lower_rank,               &
                   grid_commz, irequest2, ierror)
    ! send to upper, receive from lower
    call mpi_irecv(buf_recv(1,2), k, mpi_wp_real, lower_rank, &
                   my_z_rank+2*nprocz+lower_rank,             &
                   grid_commz, irequest3, ierror)
    call mpi_isend(buf_send(1,2), k, mpi_wp_real, upper_rank, &
                   upper_rank+2*nprocz+my_z_rank,             &
                   grid_commz, irequest4, ierror)
    call mpi_wait(irequest1, istatus, ierror)
    call mpi_wait(irequest2, istatus, ierror)
    call mpi_wait(irequest3, istatus, ierror)
    call mpi_wait(irequest4, istatus, ierror)
#else
    buf_recv(1:k,1) = buf_send(1:k,1)
    buf_recv(1:k,2) = buf_send(1:k,2)
#endif

    k = 0
    do iy = 1, nlocaly+2*n_bs
      do ix = 1, nlocalx+2*n_bs
        do iz = 1, n_bs
          qdf(ix,iy,nlocalz+iz,1) = &     !receive from upper
            qdf(ix,iy,nlocalz+iz,1) + buf_recv(iz+k,1)
          qdf(ix,iy,n_bs+iz,1)       = &     !receive from lower
           qdf(ix,iy,n_bs+iz,1) + buf_recv(iz+k,2)
        end do
        k = k + n_bs
      end do
    end do
    do iy = 1, nlocaly+2*n_bs
      do ix = 1, nlocalx+2*n_bs
        do iz = 1, n_bs
          ldf(ix,iy,nlocalz+iz,1) = &     !receive from upper
            ldf(ix,iy,nlocalz+iz,1) + buf_recv(iz+k,1)
          ldf(ix,iy,n_bs+iz,1)       = &     !receive from lower
           ldf(ix,iy,n_bs+iz,1) + buf_recv(iz+k,2)
        end do
        k = k + n_bs
      end do
    end do

    ! y direction
    !
    k = 0
    do iz = n_bs+1, n_bs+nlocalz
      do ix = 1, nlocalx+2*n_bs
        do iy = 1, n_bs
          buf_send(iy+k,1) = qdf(ix,iy,iz,1)
          buf_send(iy+k,2) = qdf(ix,n_bs+nlocaly+iy,iz,1)
        end do
        k = k + n_bs
      end do
    end do
    do iz = n_bs+1, n_bs+nlocalz
      do ix = 1, nlocalx+2*n_bs
        do iy = 1, n_bs
          buf_send(iy+k,1) = ldf(ix,iy,iz,1)
          buf_send(iy+k,2) = ldf(ix,n_bs+nlocaly+iy,iz,1)
        end do
        k = k + n_bs
      end do
    end do


#ifdef HAVE_MPI_GENESIS
    upper_rank = mod(my_y_rank+1,nprocy)
    lower_rank = mod(my_y_rank-1+nprocy,nprocy)
    ! send to lower, receive from upper
    call mpi_irecv(buf_recv(1,1), k, mpi_wp_real, upper_rank, &
                   upper_rank+nprocy+my_y_rank,               &
                   grid_commy, irequest1, ierror)
    call mpi_isend(buf_send(1,1), k, mpi_wp_real, lower_rank, &
                   my_y_rank+nprocy+lower_rank,               &
                   grid_commy, irequest2, ierror)
    ! send to upper, receive from lower
    call mpi_irecv(buf_recv(1,2), k, mpi_wp_real, lower_rank, &
                   my_y_rank+2*nprocy+lower_rank,             &
                   grid_commy, irequest3, ierror)
    call mpi_isend(buf_send(1,2), k, mpi_wp_real, upper_rank, &
                   upper_rank+2*nprocy+my_y_rank,             &
                   grid_commy, irequest4, ierror)
    call mpi_wait(irequest1, istatus, ierror)
    call mpi_wait(irequest2, istatus, ierror)
    call mpi_wait(irequest3, istatus, ierror)
    call mpi_wait(irequest4, istatus, ierror)
#else
    buf_recv(1:k,1) = buf_send(1:k,1)
    buf_recv(1:k,2) = buf_send(1:k,2)
#endif

    k = 0
    do iz = n_bs+1, n_bs+nlocalz
      do ix = 1, nlocalx+2*n_bs
        do iy = 1, n_bs
          qdf(ix,nlocaly+iy,iz,1) = &
            qdf(ix,nlocaly+iy,iz,1) + buf_recv(iy+k,1)
          qdf(ix,n_bs+iy,iz,1)       = &
            qdf(ix,n_bs+iy,iz,1) + buf_recv(iy+k,2)
        end do
        k = k + n_bs
      end do
    end do
    do iz = n_bs+1, n_bs+nlocalz
      do ix = 1, nlocalx+2*n_bs
        do iy = 1, n_bs
          ldf(ix,nlocaly+iy,iz,1) = &
            ldf(ix,nlocaly+iy,iz,1) + buf_recv(iy+k,1)
          ldf(ix,n_bs+iy,iz,1)       = &
            ldf(ix,n_bs+iy,iz,1) + buf_recv(iy+k,2)
        end do
        k = k + n_bs
      end do
    end do

    ! x direction
    !
    k = 0
    do iz = n_bs+1, n_bs+nlocalz
      do iy = n_bs+1, n_bs+nlocaly
        do ix = 1, n_bs
          buf_send(ix+k,1) = qdf(ix,iy,iz,1)
          buf_send(ix+k,2) = qdf(n_bs+nlocalx+ix,iy,iz,1)
        end do
        k = k + n_bs
      end do
    end do
    do iz = n_bs+1, n_bs+nlocalz
      do iy = n_bs+1, n_bs+nlocaly
        do ix = 1, n_bs
          buf_send(ix+k,1) = ldf(ix,iy,iz,1)
          buf_send(ix+k,2) = ldf(n_bs+nlocalx+ix,iy,iz,1)
        end do
        k = k + n_bs
      end do
    end do

#ifdef HAVE_MPI_GENESIS
    upper_rank = mod(my_x_rank+1,nprocx)
    lower_rank = mod(my_x_rank-1+nprocx,nprocx)
    ! send to lower, receive from upper
    call mpi_irecv(buf_recv(1,1), k, mpi_wp_real, upper_rank, &
                   upper_rank+nprocx+my_x_rank,               &
                   grid_commx, irequest1, ierror)
    call mpi_isend(buf_send(1,1), k, mpi_wp_real, lower_rank, &
                   my_x_rank+nprocx+lower_rank,               &
                   grid_commx, irequest2, ierror)
    ! send to upper, receive from lower
    call mpi_irecv(buf_recv(1,2), k, mpi_wp_real, lower_rank, &
                   my_x_rank+2*nprocx+lower_rank,             &
                   grid_commx, irequest3, ierror)
    call mpi_isend(buf_send(1,2), k, mpi_wp_real, upper_rank, &
                   upper_rank+2*nprocx+my_x_rank,             &
                   grid_commx, irequest4, ierror)
    call mpi_wait(irequest1, istatus, ierror)
    call mpi_wait(irequest2, istatus, ierror)
    call mpi_wait(irequest3, istatus, ierror)
    call mpi_wait(irequest4, istatus, ierror)
#else
    buf_recv(1:k,1) = buf_send(1:k,1)
    buf_recv(1:k,2) = buf_send(1:k,2)
#endif

    k = 0
    do iz = n_bs+1, n_bs+nlocalz
      do iy = n_bs+1, n_bs+nlocaly
        do ix = 1, n_bs
          qdf(nlocalx+ix,iy,iz,1) = &
            qdf(nlocalx+ix,iy,iz,1) + buf_recv(ix+k,1)
          qdf(n_bs+ix,iy,iz,1) =       &
            qdf(n_bs+ix,iy,iz,1) + buf_recv(ix+k,2)
        end do
        k = k + n_bs
      end do
    end do
    do iz = n_bs+1, n_bs+nlocalz
      do iy = n_bs+1, n_bs+nlocaly
        do ix = 1, n_bs
          ldf(nlocalx+ix,iy,iz,1) = &
            ldf(nlocalx+ix,iy,iz,1) + buf_recv(ix+k,1)
          ldf(n_bs+ix,iy,iz,1) =       &
            ldf(n_bs+ix,iy,iz,1) + buf_recv(ix+k,2)
        end do
        k = k + n_bs
      end do
    end do

    return

  end subroutine communicate_pme_pre_lj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    communicate_pme_post_lj
  !> @brief        communicate boundary grid data after bFFT
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine communicate_pme_post_lj

    ! local variable
    integer                  :: k, ix, iy, iz
    integer                  :: irequest1, irequest2, irequest3, irequest4
    integer                  :: upper_rank, lower_rank
#ifdef HAVE_MPI_GENESIS
    integer                  :: istatus(mpi_status_size)
#endif


    ! x direction
    !
    k = 0
    do iz = n_bs+1, n_bs+nlocalz
      do iy = n_bs+1, n_bs+nlocaly
        do ix = 1, n_bs
          buf_send(ix+k,1) = qdf(n_bs+ix,iy,iz,1)     ! send to lower
          buf_send(ix+k,2) = qdf(nlocalx+ix,iy,iz,1)  ! send to upper
        end do
        k = k + n_bs
      end do
    end do
    do iz = n_bs+1, n_bs+nlocalz
      do iy = n_bs+1, n_bs+nlocaly
        do ix = 1, n_bs
          buf_send(ix+k,1) = ldf(n_bs+ix,iy,iz,1)     ! send to lower
          buf_send(ix+k,2) = ldf(nlocalx+ix,iy,iz,1)  ! send to upper
        end do
        k = k + n_bs
      end do
    end do

#ifdef HAVE_MPI_GENESIS
    upper_rank = mod(my_x_rank+1,nprocx)
    lower_rank = mod(my_x_rank-1+nprocx,nprocx)
    ! send to lower, receive from upper
    call mpi_irecv(buf_recv(1,1), k, mpi_wp_real, upper_rank, &
                   upper_rank+nprocx+my_x_rank,               &
                   grid_commx, irequest1, ierror)
    call mpi_isend(buf_send(1,1), k, mpi_wp_real, lower_rank, &
                   my_x_rank+nprocx+lower_rank,               &
                   grid_commx, irequest2, ierror)
    ! send to upper, receive from lower
    call mpi_irecv(buf_recv(1,2), k, mpi_wp_real, lower_rank, &
                   my_x_rank+2*nprocx+lower_rank,             &
                   grid_commx, irequest3, ierror)
    call mpi_isend(buf_send(1,2), k, mpi_wp_real, upper_rank, &
                   upper_rank+2*nprocx+my_x_rank,             &
                   grid_commx, irequest4, ierror)
    call mpi_wait(irequest1, istatus, ierror)
    call mpi_wait(irequest2, istatus, ierror)
    call mpi_wait(irequest3, istatus, ierror)
    call mpi_wait(irequest4, istatus, ierror)
#else
    buf_recv(1:k,1) = buf_send(1:k,1)
    buf_recv(1:k,2) = buf_send(1:k,2)
#endif

    k = 0
    do iz = n_bs+1, n_bs+nlocalz
      do iy = n_bs+1, n_bs+nlocaly
        do ix = 1, n_bs
          qdf(n_bs+nlocalx+ix,iy,iz,1) = buf_recv(ix+k,1) ! receive from upper
          qdf(ix,iy,iz,1)  = buf_recv(ix+k,2) ! receive from lower
        end do
        k = k + n_bs
      end do
    end do
    do iz = n_bs+1, n_bs+nlocalz
      do iy = n_bs+1, n_bs+nlocaly
        do ix = 1, n_bs
          ldf(n_bs+nlocalx+ix,iy,iz,1) = buf_recv(ix+k,1) ! receive from upper
          ldf(ix,iy,iz,1)  = buf_recv(ix+k,2) ! receive from lower
        end do
        k = k + n_bs
      end do
    end do

    ! y direction
    !
    k = 0
    do iz = n_bs+1, n_bs+nlocalz
      do ix = 1, nlocalx+2*n_bs
        do iy = 1, n_bs
          buf_send(iy+k,1) = qdf(ix,n_bs+iy,iz,1)      ! send to lower
          buf_send(iy+k,2) = qdf(ix,nlocaly+iy,iz,1)   ! send to upper
        end do
        k = k + n_bs
      end do
    end do
    do iz = n_bs+1, n_bs+nlocalz
      do ix = 1, nlocalx+2*n_bs
        do iy = 1, n_bs
          buf_send(iy+k,1) = ldf(ix,n_bs+iy,iz,1)      ! send to lower
          buf_send(iy+k,2) = ldf(ix,nlocaly+iy,iz,1)   ! send to upper
        end do
        k = k + n_bs
      end do
    end do

#ifdef HAVE_MPI_GENESIS
    upper_rank = mod(my_y_rank+1,nprocy)
    lower_rank = mod(my_y_rank-1+nprocy,nprocy)
    ! send to lower, receive from upper
    call mpi_irecv(buf_recv(1,1), k, mpi_wp_real, upper_rank, &
                   upper_rank+nprocy+my_y_rank,               &
                   grid_commy, irequest1, ierror)
    call mpi_isend(buf_send(1,1), k, mpi_wp_real, lower_rank, &
                   my_y_rank+nprocy+lower_rank,               &
                   grid_commy, irequest2, ierror)
    ! send to upper, receive from lower
    call mpi_irecv(buf_recv(1,2), k, mpi_wp_real, lower_rank, &
                   my_y_rank+2*nprocy+lower_rank,             &
                   grid_commy, irequest3, ierror)
    call mpi_isend(buf_send(1,2), k, mpi_wp_real, upper_rank, &
                   upper_rank+2*nprocy+my_y_rank,             &
                   grid_commy, irequest4, ierror)
    call mpi_wait(irequest1, istatus, ierror)
    call mpi_wait(irequest2, istatus, ierror)
    call mpi_wait(irequest3, istatus, ierror)
    call mpi_wait(irequest4, istatus, ierror)
#else
    buf_recv(1:k,1) = buf_send(1:k,1)
    buf_recv(1:k,2) = buf_send(1:k,2)
#endif

    k = 0
    do iz = n_bs+1, n_bs+nlocalz
      do ix = 1, nlocalx+2*n_bs
        do iy = 1, n_bs
          qdf(ix,n_bs+nlocaly+iy,iz,1) = buf_recv(iy+k,1)
          qdf(ix,iy,iz,1)    = buf_recv(iy+k,2)
        end do
        k = k + n_bs
      end do
    end do
    do iz = n_bs+1, n_bs+nlocalz
      do ix = 1, nlocalx+2*n_bs
        do iy = 1, n_bs
          ldf(ix,n_bs+nlocaly+iy,iz,1) = buf_recv(iy+k,1)
          ldf(ix,iy,iz,1)    = buf_recv(iy+k,2)
        end do
        k = k + n_bs
      end do
    end do

    ! z direction
    !
    k = 0
    do iy = 1, nlocaly+2*n_bs
      do ix = 1, nlocalx+2*n_bs
        do iz = 1, n_bs
          buf_send(iz+k,1) = qdf(ix,iy,n_bs+iz,1)  ! send to lower
          buf_send(iz+k,2) = qdf(ix,iy,nlocalz+iz,1) !send to upper
        end do
        k = k + n_bs
      end do
    end do
    do iy = 1, nlocaly+2*n_bs
      do ix = 1, nlocalx+2*n_bs
        do iz = 1, n_bs
          buf_send(iz+k,1) = ldf(ix,iy,n_bs+iz,1)  ! send to lower
          buf_send(iz+k,2) = ldf(ix,iy,nlocalz+iz,1) !send to upper
        end do
        k = k + n_bs
      end do
    end do


#ifdef HAVE_MPI_GENESIS
    upper_rank = mod(my_z_rank+1,nprocz)
    lower_rank = mod(my_z_rank-1+nprocz,nprocz)
    ! send to lower, receive from upper
    call mpi_irecv(buf_recv(1,1), k, mpi_wp_real, upper_rank, &
                   upper_rank+nprocz+my_z_rank,               &
                   grid_commz, irequest1, ierror)
    call mpi_isend(buf_send(1,1), k, mpi_wp_real, lower_rank, &
                   my_z_rank+nprocz+lower_rank,               &
                   grid_commz, irequest2, ierror)
    ! send to upper, receive from lower
    call mpi_irecv(buf_recv(1,2), k, mpi_wp_real, lower_rank, &
                   my_z_rank+2*nprocz+lower_rank,             &
                   grid_commz, irequest3, ierror)
    call mpi_isend(buf_send(1,2), k, mpi_wp_real, upper_rank, &
                   upper_rank+2*nprocz+my_z_rank,             &
                   grid_commz, irequest4, ierror)
    call mpi_wait(irequest1, istatus, ierror)
    call mpi_wait(irequest2, istatus, ierror)
    call mpi_wait(irequest3, istatus, ierror)
    call mpi_wait(irequest4, istatus, ierror)
#else
    buf_recv(1:k,1) = buf_send(1:k,1)
    buf_recv(1:k,2) = buf_send(1:k,2)
#endif

    k = 0
    do iy = 1, nlocaly+2*n_bs
      do ix = 1, nlocalx+2*n_bs
        do iz = 1, n_bs
          qdf(ix,iy,n_bs+nlocalz+iz,1) = buf_recv(iz+k,1)
          qdf(ix,iy,iz,1) = buf_recv(iz+k,2)
        end do
        k = k + n_bs
      end do
    end do
    do iy = 1, nlocaly+2*n_bs
      do ix = 1, nlocalx+2*n_bs
        do iz = 1, n_bs
          ldf(ix,iy,n_bs+nlocalz+iz,1) = buf_recv(iz+k,1)
          ldf(ix,iy,iz,1) = buf_recv(iz+k,2)
        end do
        k = k + n_bs
      end do
    end do

    return

  end subroutine communicate_pme_post_lj

end module sp_energy_pme_opt_1dalltoall_mod
