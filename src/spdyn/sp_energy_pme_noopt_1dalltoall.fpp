!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_pme_noopt_1dalltoall_mod
!> @brief   Smooth particle mesh ewald method
!! @authors Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_pme_noopt_1dalltoall_mod

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
  real(dp), public, save :: u_self_no1   ! Ewald self energy   
  real(wp),         save :: inv_volume   ! inverse of volume
  real(wp),         save :: box(3)       ! box size
  real(wp),         save :: box_inv(3)   ! Inverse of box size
  real(wp),         save :: bs_fact      ! B-spline factor (1/(n-1)...2)
  real(wp),         save :: bs_fact3     ! bs_fact^3
  real(wp),         save :: bs_fact3d    ! bs_fact^3*(n-1)
  real(wp),         save :: r_scale(3)   ! coordinate-scaling factor (I/L)
  integer,          save :: n_bs         ! Order of B-spline
  integer,          save :: ngrid(3)     ! Number of grid
  integer,          save :: nx           ! process number in x dimension
  integer,          save :: ny           ! process number in y dimension
  integer,          save :: nz           ! process number in z dimension
  integer,          save :: x_start
  integer,          save :: x_end
  integer,          save :: x_local  
  integer,          save :: x_start1
  integer,          save :: x_end1
  integer,          save :: x_local1 
  integer,          save :: y_start
  integer,          save :: y_end
  integer,          save :: y_local  
  integer,          save :: z_start
  integer,          save :: z_end
  integer,          save :: z_local  
  integer,          save :: nlocalx
  integer,          save :: nlocalx1
  integer,          save :: nlocaly
  integer,          save :: nlocalz
  integer,          save :: maxproc
  integer,          save :: ngridmax

  real(wp),         save, allocatable :: b2(:,:)        ! b^2(hx,hy,hz)
  real(wp),         save, allocatable :: gx(:)
  real(wp),         save, allocatable :: gy(:)
  real(wp),         save, allocatable :: gz(:)
  real(wp),         save, allocatable :: vir_fact(:,:,:)! -2*(1+G^2/4a)/G^2
  real(wp),         save, allocatable :: vir_fact_lj(:,:,:)
  real(wp),         save, allocatable :: theta(:,:,:)   ! prefactor in recip energy
  real(wp),         save, allocatable :: theta_lj(:,:,:) ! prefactor in recip energy (LJ case)

  real(wp),         save, allocatable :: f(:,:,:)
  integer,          save, allocatable :: vi(:,:,:)
  real(wp),         save, allocatable :: bsc(:,:,:,:)
  real(wp),         save, allocatable :: bscd(:,:,:,:)
  real(wp),         save, allocatable :: qdf(:,:,:,:)
  real(wp),         save, allocatable :: qdf_real(:)
  real(wp),         save, allocatable :: qdf_work(:,:)
  complex(wp),      save, allocatable :: ftqdf(:)
  complex(wp),      save, allocatable :: ftqdf_work(:,:)
  real(wp),         save, allocatable :: ldf(:,:,:,:)
  real(wp),         save, allocatable :: ldf_real(:)
  real(wp),         save, allocatable :: ldf_work(:,:)
  complex(wp),      save, allocatable :: ftldf(:)
  complex(wp),      save, allocatable :: ftldf_work(:,:)

  !$omp threadprivate(f)
  !$omp threadprivate(bsc)
  !$omp threadprivate(bscd)

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
  public  :: setup_pme_noopt_1dalltoall
  public  :: dealloc_pme_noopt_1dalltoall
  public  :: pme_pre_noopt_1dalltoall
  public  :: pme_pre_noopt_1dalltoall_lj
  public  :: pme_recip_noopt_1dalltoall
  public  :: pme_recip_noopt_1dalltoall_lj

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_pme_noopt_1dalltoall
  !> @brief        Setup for PME with domain decomposition
  !! @authors      JJ, NT
  !! @param[in]    domain   : domain information
  !! @param[in]    boundary : boundary information
  !! @param[inout] enefunc  : energy function
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_pme_noopt_1dalltoall(domain, boundary, enefunc)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_boundary),        intent(in)    :: boundary
    type(s_enefunc),         intent(in)    :: enefunc

    ! local variables
    real(wp)                 :: bcos, bsin, fact, grid_space
    integer                  :: i, j, k, js
    integer                  :: localmax, ncell, maxgrid
    integer                  :: iproc(3), index_x, index_y, index_z
    integer                  :: remainder, quotient, expo
    integer                  :: alloc_stat
    integer                  :: nix, niy, niz

    real(wp),    allocatable :: bs(:)


    ncell    =  domain%num_cell_local + domain%num_cell_boundary

    ! Setting the number of grids (even number)
    !
    ngrid(1) = enefunc%pme_ngrid_x
    ngrid(2) = enefunc%pme_ngrid_y
    ngrid(3) = enefunc%pme_ngrid_z
    ngridmax = max(ngrid(1),ngrid(2),ngrid(3))

    if (ngrid(1) == 0 .and. ngrid(2) == 0 .and. ngrid(3) ==0) then

      grid_space = enefunc%pme_max_spacing
      ngrid(1) = int(boundary%box_size_x/grid_space)
      ngrid(2) = int(boundary%box_size_y/grid_space)
      ngrid(3) = int(boundary%box_size_z/grid_space)
      ngridmax = max(ngrid(1),ngrid(2),ngrid(3))

    end if

    ! grid partition
    !
    if (reciprocal_calc) then

      ! Setting common parameters
      !
      el_fact   = ELECOEF/enefunc%dielec_const
      n_bs      = enefunc%pme_nspline
      alpha     = enefunc%pme_alpha
      alpha2m   = -alpha*alpha
      alpha3    = alpha*alpha*alpha
      alpha6    = alpha3*alpha3
      inv_alpha = 1.0_wp / alpha

      nx = boundary%num_domain(1)
      ny = boundary%num_domain(2)
      nz = boundary%num_domain(3)

      ! Check if charges in boundary cells are sufficient for charge
      ! grid data in each subdomain
      !
      grid_space = boundary%box_size_x / real(ngrid(1),wp)
      ngridMax   = int(boundary%cell_size_x / grid_space)
      grid_space = boundary%box_size_y / real(ngrid(2),wp)
      ngridMax   = min(ngridMax,int(boundary%cell_size_y/grid_space))
      grid_space = boundary%box_size_z / real(ngrid(3),wp)
      ngridMax   = min(ngridMax,int(boundary%cell_size_z/grid_space))
      if (ngridmax < n_bs) &
        call error_msg('Cell size should be larger to obtain the charge'// &
                       ' data from neighboring cells')

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

      x_local  = nlocalx
      y_local  = nlocaly
      z_local  = nlocalz

      localmax = max(nlocalx*nlocaly,nlocalx*nlocalz,nlocaly*nlocalz)
      index_x  = iproc(2)*nlocaly*nlocalz + iproc(3)
      index_y  = iproc(3)*nlocalx*nlocalz + iproc(1)
      index_z  = iproc(1)*nlocalx*nlocaly + iproc(2)

      maxproc = max(nprocx,nprocy,nprocz,1)

      nlocalx1 = nlocalx/2 + 1
      if (my_x_rank == 0) then
        x_local1 = nlocalx/2+1
        x_start1 = 1
        x_end1   = nlocalx/2+1
      else
        x_local1 = nlocalx/2
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
      bs_fact3  = bs_fact*bs_fact*bs_fact
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

          b2(j, k) = 1.0_wp/(bcos**2 + bsin**2)

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

      ! Allocation
      !
      if (enefunc%vdw == VDWPME) then
        allocate(gx(x_local1), gy(nlocaly), gz(ngrid(3)),       &
                 vir_fact(ngrid(3), nlocaly, x_local1),         &
                 vir_fact_lj(ngrid(3), nlocaly, x_local1),      &
                 vi(3,MaxAtom,ncell),                           &
                 qdf(nlocalx,nlocaly,nlocalz,nthread),          &
                 ldf(nlocalx,nlocaly,nlocalz,nthread),          &
                 qdf_real(nlocalx*nlocaly*nlocalz),             &
                 ldf_real(nlocalx*nlocaly*nlocalz),             &
                 theta(ngrid(3),nlocaly,x_local1),              &
                 theta_lj(ngrid(3),nlocaly,x_local1),           &
                 qdf_work(nlocalx*nlocaly*nlocalz,maxproc),     &
                 ldf_work(nlocalx*nlocaly*nlocalz,maxproc),     &
                 ftqdf(nlocalx1*nlocaly*nlocalz),               &
                 ftldf(nlocalx1*nlocaly*nlocalz),               &
                 ftqdf_work(nlocalx1*nlocaly*nlocalz*maxproc,1),&
                 ftldf_work(nlocalx1*nlocaly*nlocalz*maxproc,1),&
                 ftqdf_work2(ngrid(3),nlocalx1,nlocaly),        &
                 ftldf_work2(ngrid(3),nlocalx1,nlocaly),        &
                 stat = alloc_stat)
      else
        allocate(gx(x_local1), gy(nlocaly), gz(ngrid(3)),       &
                 vir_fact(ngrid(3), nlocaly, x_local1),         &
                 vi(3,MaxAtom,ncell),                           &
                 qdf(nlocalx,nlocaly,nlocalz,nthread),          &
                 qdf_real(nlocalx*nlocaly*nlocalz),             &
                 theta(ngrid(3),nlocaly,x_local1),              &
                 qdf_work(nlocalx*nlocaly*nlocalz,maxproc),     &
                 ftqdf(nlocalx1*nlocaly*nlocalz),               &
                 ftqdf_work(nlocalx1*nlocaly*nlocalz*maxproc,1),&
                 ftqdf_work2(ngrid(3),nlocalx1,nlocaly),        &
                 stat = alloc_stat)
      end if
      if (alloc_stat /= 0) call error_msg_alloc

      !$omp parallel private(alloc_stat)
      allocate(f(3,MaxAtom,ncell),                                       &
               bsc(n_bs,3,MaxAtom,ncell),bscd(n_bs,3,MaxAtom,ncell),     &
               stat = alloc_stat)
      if (alloc_stat /= 0)   call error_msg_alloc
      !$omp end parallel

      deallocate(bs)

    end if

    return

  end subroutine setup_pme_noopt_1dalltoall

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_pme_noopt_1dalltoall
  !> @brief        deallocate all arrays used in PME
  !! @authors      JJ, NT
  !! @param[in]    enefunc : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_pme_noopt_1dalltoall(enefunc)

    type(s_enefunc),           intent(in) :: enefunc

    ! local variables
    integer :: dealloc_stat, ierr


    if (.not. allocated(f)) &
      return 

    !$omp parallel private(dealloc_stat)
    deallocate(f, bsc, bscd, stat = dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    !$omp end parallel

    if (enefunc%vdw == VDWPME) then
      deallocate(gx, gy, gz, vir_fact, vir_fact_lj, vi,         &
                 qdf, qdf_real, ldf, ldf_real, theta, theta_lj, &
                 qdf_work, ftqdf, ftqdf_work, ftqdf_work2,      &
                 ldf_work, ftldf, ftldf_work, ftldf_work2,      &
                 stat = dealloc_stat)
    else
      deallocate(gx, gy, gz, vir_fact, vi, qdf, qdf_real, &
                 theta, qdf_work, ftqdf, ftqdf_work,      &
                 ftqdf_work2, stat = dealloc_stat)
    end if
    if (dealloc_stat /= 0) call error_msg_dealloc

    deallocate(b2, stat = dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_pme_noopt_1dalltoall

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_pre_noopt_1dalltoall
  !> @brief        Prepare functions for PME calculation
  !! @authors      JJ, NT
  !! @param[in]    domain   : domain information
  !! @param[in]    boundary : boundary information
  !! @note         Extracted from setup_pme for NPT calculation
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_pre_noopt_1dalltoall(domain, boundary)

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

    inv_volume = box_inv(1)*box_inv(2)*box_inv(3)
    vol_fact2 = 2.0_wp * PI * el_fact * inv_volume
    vol_fact4 = 2.0_wp * vol_fact2

    ! Prepareing theta (prefactor of reciprocal energy)
    ! Gx, Gy, Gz, and vir_fact=-2(1+G^2/4a)/G^2 are also prepared 
    ! for virial calculation
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
          if (g2 > EPS) then
            vir_fact(k,iy,ix) = -2.0_wp * (1.0_wp - g2 * fact)/g2
            if (g2*fact < -80.0_wp) then
              theta(k,iy,ix) = 0.0_wp
            else
              theta(k,iy,ix) = b2(i,1) * b2(j,2) * b2(k,3) * exp(g2 * fact) / g2
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
          g2 = gx(ix)*gx(ix) + gy(iy)*gy(iy) + gz(k)*gz(k)
          if (g2 > EPS) then
            vir_fact(k,iy,ix) = -2.0_wp * (1.0_wp - g2 * fact)/g2
            if (g2*fact < -80.0_wp) then
              theta(k,iy,ix) = 0.0_wp
            else
              theta(k,iy,ix) = b2(i,1) * b2(j,2) * b2(k,3) * exp(g2 * fact) / g2
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
    u_self_no1  = 0.0_dp

    do i = 1, domain%num_cell_local
      do ix = 1, domain%num_atom(i)
        u_self_no1 = u_self_no1 + domain%charge(ix,i)*domain%charge(ix,i)
      end do 
    end do

    u_self_no1 = - u_self_no1 * el_fact * alpha/sqrt(PI)

    return

  end subroutine pme_pre_noopt_1dalltoall

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_pre_noopt_1dalltoall_lj
  !> @brief        Prepare functions for PME calculation (Elec and LJ)
  !! @authors      JJ, NT
  !! @param[in]    domain   : domain information
  !! @param[in]    boundary : boundary information
  !! @note         Extracted from setup_pme for NPT calculation
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_pre_noopt_1dalltoall_lj(domain, boundary)

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

    ! factor that convert from m to G=2*pi*m/L
    !
    gfact(1:3) = 2.0_wp * PI * box_inv(1:3)

    ! -1/(4*alpha^2)
    !
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

          ! g, g^2, g^3
          !
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

          ! calculate theta (prefacor in the reciprocal energy) &
          ! calculate vir_fact (prefacor in the reciprocal force) 
          !
          g_lj = (4.0_wp*alpha3-2.0_wp*alpha*g2)*fact2 &
               + SQRT_PI*g3*fact4
          theta_lj(k,iy,ix) = fact5 * g_lj
          vir_fact_lj(k,iy,ix) = 3.0_wp*g1*SQRT_PI*fact4-6.0_wp*alpha*fact2
          vir_fact_lj(k,iy,ix) = vir_fact_lj(k,iy,ix) * fact5
          if (g2 > EPS) then
            vir_fact(k,iy,ix) = -2.0_wp * (1.0_wp-fact1)/g2
            theta(k,iy,ix) = fact5 * fact2/g2
            if (abs(theta(k,iy,ix)) < 1.0e-15) theta(k,iy,ix) = 0.0_wp
          else
            vir_fact(k,iy,ix)    = 0.0_wp
            theta(k,iy,ix)       = 0.0_wp
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

          ! calculate theta (prefacor in the reciprocal energy) &
          ! calculate vir_fact (prefacor in the reciprocal force) 
          !
          g_lj = (4.0_wp*alpha3-2.0_wp*alpha*g2)*fact2 &
               + SQRT_PI*g3*fact4
          theta_lj(k,iy,ix) = fact5 * g_lj
          vir_fact_lj(k,iy,ix) = 3.0_wp*g1*SQRT_PI*fact4-6.0_wp*alpha*fact2
          vir_fact_lj(k,iy,ix) = vir_fact_lj(k,iy,ix) * fact5
          if (g2 > EPS) then
            vir_fact(k,iy,ix) = -2.0_wp * (1.0_wp-fact1)/g2
            theta(k,iy,ix) = fact5 * fact2/g2
            if (abs(theta(k,iy,ix)) < 1.0e-15) theta(k,iy,ix) = 0.0_wp
          else
            vir_fact(k,iy,ix)    = 0.0_wp
            theta(k,iy,ix)       = 0.0_wp
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

    ! Calculating self energy (self energy of LJ => Later)
    !
    u_self_no1  = 0.0_dp

    do i = 1, domain%num_cell_local
      do ix = 1, domain%num_atom(i)
        u_self_no1 = u_self_no1 + domain%charge(ix,i)*domain%charge(ix,i)
      end do 
    end do

    u_self_no1  = - u_self_no1 * el_fact * alpha/sqrt(PI)

    return

  end subroutine pme_pre_noopt_1dalltoall_lj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_recip_noopt_1dalltoall
  !> @brief        Calculate PME reciprocal part with domain decomposition
  !! @authors      JJ, NT
  !! @param[in]    domain : domain information
  !! @param[inout] force  : forces of target systems
  !! @param[inout] virial : virial term of target systems
  !! @param[inout] eelec  : electrostatic energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_recip_noopt_1dalltoall(domain, force, virial, eelec)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    real(wip),               intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: elec_temp
    real(wp)                 :: vr(3), dv(3), grid(3), half_grid(3)
    real(wp)                 :: bsc_tmp(3,n_bs), bscd_tmp(3,n_bs)
    real(wp)                 :: vxx, vyy, vzz, temp
    real(wp)                 :: tqq, qtmp
    real(wp)                 :: force_local(3) 
    real(wp)                 :: vol_factor
    integer                  :: i, k, k1, icel, iproc, id
    integer                  :: ix, iv, iy, iz, ixs, iys, izs
    integer                  :: ixyz, nix, nix1, niy, niz, nizx, nizy
    integer                  :: ixx, iyy, izz, ii(3), is(n_bs,3)
    integer                  :: start(3), end(3)
    integer                  :: vx_tmp(n_bs), vy_tmp(n_bs), vz_tmp(n_bs)
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
    !$omp         ixx, ixs, vxx, vyy, vzz, iproc, tqq, k1, temp, vr, &
    !$omp         work1, work2, elec_temp, force_local, ixyz, vx_tmp, vy_tmp,  &
    !$omp         vz_tmp, nix, nix1, niy, niz, bsc_tmp, bscd_tmp, kk, iprocx,  &
    !$omp         iprocy, qtmp, grid, half_grid, is, start, end, nizx, nizy)
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
    do iz = 1, nlocalz
      do iy = 1, nlocaly
        do ix = 1, nlocalx
          qdf(ix,iy,iz,id+1) = 0.0_wp
        end do
      end do
    end do

    ! Calculateing Q-fanction
    !
    grid(1:3) = real(ngrid(1:3),wp)
    half_grid(1:3) = real(ngrid(1:3)/2,wp)

    do icel = id+1, ncell_local, nthread
      do i = 1, natom(icel)

        qtmp = charge(i,icel)
        vr(1:3) = coord(1:3,i,icel) * r_scale(1:3)
        vr(1:3) = vr(1:3) + half_grid(1:3) &
                - grid(1:3)*anint(vr(1:3)/grid(1:3))

        do k = 1, 3
          vr(k) = min(vr(k),grid(k)-0.000001_wp)
          ii(k) = min(int(vr(k)),ngrid(k)-1)
        end do

        nix = 0
        niy = 0
        niz = 0
!ocl nosimd
        do ixyz = 1, n_bs
          ixs = ii(1) - ixyz + 2
          if (ixs <= 0) ixs = ixs + ngrid(1)
          if (ixs >= x_start .and. ixs <= x_end) then
            nix = nix + 1
            vx_tmp(nix) = ixs - x_start + 1
            is(nix,1) = ixyz
          end if
        end do

!ocl nosimd
        do ixyz = 1, n_bs
          iys = ii(2) - ixyz + 2
          if (iys <= 0) iys = iys + ngrid(2)
          if (iys >= y_start .and. iys <= y_end) then
            niy = niy + 1
            vy_tmp(niy) = iys - y_start + 1
            is(niy,2) = ixyz
          end if
        end do

!ocl nosimd
        do ixyz = 1, n_bs
          izs = ii(3) - ixyz + 2
          if (izs <= 0) izs = izs + ngrid(3)
          if (izs >= z_start .and. izs <= z_end) then
            niz = niz + 1
            vz_tmp(niz) = izs - z_start + 1
            is(niz,3) = ixyz
          end if
        end do

        force(1,i,icel) = 1.0_wp

        do k = 1, 3
          vi(k,i,icel) = ii(k)
          dv(k) = vr(k) - real(ii(k),wp)
          call b_spline_dev_coef(n_bs, dv(k), bsc(1,k,i,icel), bscd(1,k,i,icel))
        end do

!ocl nosimd
        do ix = 1, nix
          ixyz = is(ix,1)
          bsc_tmp(1,ix) = bsc(ixyz,1,i,icel)
        end do
!ocl nosimd
        do iy = 1, niy
          ixyz = is(iy,2)
          bsc_tmp(2,iy) = bsc(ixyz,2,i,icel)
        end do
!ocl nosimd
        do iz = 1, niz
          ixyz = is(iz,3)
          bsc_tmp(3,iz) = bsc(ixyz,3,i,icel)
        end do

        do iz = 1, niz
          izs = vz_tmp(iz)
          do iy = 1, niy
            iys = vy_tmp(iy)
!ocl norecurrence(qdf)
!ocl nosimd
             do ix = 1, nix
               ixs = vx_tmp(ix)
              qdf(ixs,iys,izs,id+1) = qdf(ixs,iys,izs,id+1)           &
                                    + bsc_tmp(1,ix)*bsc_tmp(2,iy)     &
                                     * bsc_tmp(3,iz)*qtmp   &
                                     * bs_fact3
            end do
          end do
        end do
      end do
    end do

    do icel = id+ncell_local+1, ncell, nthread
      do i = 1, natom(icel)

        qtmp = charge(i,icel)
        vr(1:3) = coord(1:3,i,icel) * r_scale(1:3)
        vr(1:3) = vr(1:3) + half_grid(1:3) &
                - grid(1:3)*anint(vr(1:3)/grid(1:3))

!ocl nosimd
        do k = 1, 3
          vr(k) = min(vr(k),grid(k)-0.000001_wp)
          ii(k) = min(int(vr(k)),ngrid(k)-1)
        end do

        nix = 0
        niy = 0
        niz = 0
 
        ixs = ii(1) - n_bs + 2
        if (ixs <= 0) ixs = ixs + ngrid(1)
        if (ixs >= x_start .and. ixs <= x_end) then
          nix = nix + 1
          vx_tmp(nix) = ixs - x_start + 1
          is (nix,1) = n_bs
        end if
        ixs = ii(1) + 1
        if (ixs <= 0) ixs = ixs + ngrid(1)
        if (ixs >= x_start .and. ixs <= x_end) then
          nix = nix + 1
          vx_tmp(nix) = ixs - x_start + 1
          is (nix,1) = 1
        end if
        if (nix == 0) cycle

        iys = ii(2) - n_bs + 2
        if (iys <= 0) iys = iys + ngrid(2)
        if (iys >= y_start .and. iys <= y_end) then
          niy = niy + 1
          vy_tmp(niy) = iys - y_start + 1
          is (niy,2) = n_bs
        end if
        iys = ii(2) + 1
        if (iys <= 0) iys = iys + ngrid(2)
        if (iys >= y_start .and. iys <= y_end) then
          niy = niy + 1
          vy_tmp(niy) = iys - y_start + 1
          is (niy,2) = 1
        end if
        if (niy == 0) cycle

        izs = ii(3) - n_bs + 2
        if (izs <= 0) izs = izs + ngrid(3)
        if (izs >= z_start .and. izs <= z_end) then
          niz = niz + 1
          vz_tmp(niz) = izs - z_start + 1
          is (niz,3) = n_bs
        end if
        izs = ii(3) + 1
        if (izs <= 0) izs = izs + ngrid(3)
        if (izs >= z_start .and. izs <= z_end) then
          niz = niz + 1
          vz_tmp(niz) = izs - z_start + 1
          is (niz,3) = 1
        end if
        if (niz == 0) cycle

!ocl nosimd
        do ixyz = 2, n_bs-1
          ixs = ii(1) - ixyz + 2
          if (ixs <= 0) ixs = ixs + ngrid(1)
          if (ixs >= x_start .and. ixs <= x_end) then
            nix = nix + 1
            vx_tmp(nix) = ixs - x_start + 1
            is(nix,1) = ixyz
          end if
        end do

!ocl nosimd
        do ixyz = 2, n_bs-1
          iys = ii(2) - ixyz + 2
          if (iys <= 0) iys = iys + ngrid(2)
          if (iys >= y_start .and. iys <= y_end) then
            niy = niy + 1
            vy_tmp(niy) = iys - y_start + 1
            is(niy,2) = ixyz
          end if
        end do

!ocl nosimd
        do ixyz = 2, n_bs-1
          izs = ii(3) - ixyz + 2
          if (izs <= 0) izs = izs + ngrid(3)
          if (izs >= z_start .and. izs <= z_end) then
            niz = niz + 1
            vz_tmp(niz) = izs - z_start + 1
            is(niz,3) = ixyz
          end if
        end do

        force(1,i,icel) = 1.0_wp

        do k = 1, 3
          vi(k,i,icel) = ii(k)
          dv(k) = vr(k) - real(ii(k),wp)
          call b_spline_dev_coef(n_bs, dv(k), bsc(1,k,i,icel), bscd(1,k,i,icel))
        end do

!ocl nosimd
        do ix = 1, nix
          ixyz = is(ix,1)
          bsc_tmp(1,ix) = bsc(ixyz,1,i,icel)
        end do
!ocl nosimd
        do iy = 1, niy
          ixyz = is(iy,2)
          bsc_tmp(2,iy) = bsc(ixyz,2,i,icel)
        end do
!ocl nosimd
        do iz = 1, niz
          ixyz = is(iz,3)
          bsc_tmp(3,iz) = bsc(ixyz,3,i,icel)
        end do

        do iz = 1, niz
          izs = vz_tmp(iz)
          do iy = 1, niy
            iys = vy_tmp(iy)
!ocl norecurrence(qdf)
!ocl nosimd
             do ix = 1, nix
               ixs = vx_tmp(ix)
              qdf(ixs,iys,izs,id+1) = qdf(ixs,iys,izs,id+1)           &
                                    + bsc_tmp(1,ix)*bsc_tmp(2,iy)     &
                                     * bsc_tmp(3,iz)*qtmp   &
                                     * bs_fact3
            end do
          end do
        end do
      end do
    end do
    !$omp barrier
    do iproc = 1, nthread
      !$omp do
      do iz = 1, nlocalz
        do iy = 1, nlocaly
          k = (iy-1)*nlocalx + (iz-1)*nlocalx*nlocaly
          do ix = 1, nlocalx
            qdf_real(k+ix) = qdf_real(k+ix) + qdf(ix,iy,iz,iproc)
          end do
        end do
      end do
    end do

    ! FFT (forward)
    !
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
          tqq = real(ftqdf_work2(iz,ix,iy),wp)**2 &
              + imag(ftqdf_work2(iz,ix,iy))**2
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
          tqq = real(ftqdf_work2(iz,1,iy),wp)**2 + &
                imag(ftqdf_work2(iz,1,iy))**2
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
          tqq = real(ftqdf_work2(iz,ix,iy),wp)**2 + &
                imag(ftqdf_work2(iz,ix,iy))**2
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

    ! F^-1[Th]*F^-1[Q] (=X)
    !
    !$omp barrier
    do ix = 1, x_local1
      do iy = id+1, niy, nthread
        k = (iy-1)*ngrid(3) + (ix-1)*ngrid(3)*niy
        do iz = 1, ngrid(3)
          ftqdf_work(k+iz,1) = ftqdf_work2(iz,ix,iy)   &
                             * cmplx(theta(iz,iy,ix),0.0_wp)
        end do
      end do
    end do

    ! FFT(backward)
    !
    !$omp barrier
    call bfft3d_1d_alltoall(qdf_work, qdf_real, ftqdf, ftqdf, ftqdf_work,    &
                    ftqdf_work, ftqdf_work2, ftqdf_work, work1, work2,       &
                    nlocalx, nlocaly, nlocalz, nlocalx1, x_local1, ngrid(1), &
                    ngrid(2), ngrid(3), nprocx, nprocy, nprocz, id, nthread, &
                    my_x_rank, grid_commx, grid_commy, grid_commz, niy,      &
                    nizy, nizx)

    !$omp barrier

    do iz = id+1, nlocalz, nthread
      do iy = 1, nlocaly
        k = (iy-1)*nlocalx + (iz-1)*nlocalx*nlocaly
        do ix = 1, nlocalx
          qdf(ix,iy,iz,1) = qdf_real(k+ix)
        end do
      end do
    end do

    !$omp barrier

    do icel = id+1, ncell_local, nthread
      do i = 1, natom(icel)

        force_local(1:3) = 0.0_wp
        qtmp = charge(i,icel)

!ocl nosimd
        do k = 1, 3
          ii(k) = vi(k,i,icel)
        end do

        nix=0
        niy=0
        niz=0

!ocl nosimd
        do ixyz = 1, n_bs

          ixs = ii(1) - (ixyz-1) + 1
          iys = ii(2) - (ixyz-1) + 1
          izs = ii(3) - (ixyz-1) + 1
          if (ixs <= 0) ixs = ixs + ngrid(1)
          if (iys <= 0) iys = iys + ngrid(2)
          if (izs <= 0) izs = izs + ngrid(3)

          if (ixs >= x_start .and. ixs <= x_end) then
            nix = nix + 1
            vx_tmp(nix) = ixs - x_start + 1
            bsc_tmp(1,nix)  = bsc(ixyz,1,i,icel)
            bscd_tmp(1,nix) = bscd(ixyz,1,i,icel)
          end if

          if (iys >= y_start .and. iys <= y_end) then
            niy = niy + 1
            vy_tmp(niy) = iys - y_start + 1
            bsc_tmp(2,niy)  = bsc(ixyz,2,i,icel)
            bscd_tmp(2,niy) = bscd(ixyz,2,i,icel)
          end if

          if (izs >= z_start .and. izs <= z_end) then
            niz = niz + 1
            vz_tmp(niz) = izs - z_start + 1
            bsc_tmp(3,niz)  = bsc(ixyz,3,i,icel)
            bscd_tmp(3,niz) = bscd(ixyz,3,i,icel)
          end if
        end do
       
!ocl nosimd 
        do iz = 1, niz
          izs = vz_tmp(iz)
          do iy = 1, niy
            iys = vy_tmp(iy)
            do ix = 1, nix
              ixs = vx_tmp(ix)

              force_local(1) = force_local(1)            &
                  + bscd_tmp(1,ix)*bsc_tmp(2,iy)         &
                  * bsc_tmp(3,iz)*qdf(ixs,iys,izs,1)
              force_local(2) = force_local(2)            &
                  + bsc_tmp(1,ix)*bscd_tmp(2,iy)         &
                  * bsc_tmp(3,iz)*qdf(ixs,iys,izs,1)
              force_local(3) = force_local(3)            &
                  + bsc_tmp(1,ix)*bsc_tmp(2,iy)          &
                  * bscd_tmp(3,iz)*qdf(ixs,iys,izs,1)
            end do
          end do
        end do

        force(1:3,i,icel) = - force_local(1:3) * &
                              qtmp*r_scale(1:3)*vol_fact4*bs_fact3d

      end do
    end do

    do icel = id+ncell_local+1, ncell, nthread
      do i = 1, natom(icel)

        if (abs(force(1,i,icel)) < EPS) cycle

        force_local(1:3) = 0.0_wp
        qtmp = charge(i,icel)

        do k = 1, 3
          ii(k) = vi(k,i,icel)
        end do

        nix=0
        niy=0
        niz=0

!ocl nosimd
        do ixyz = 1, n_bs

          ixs = ii(1) - (ixyz-1) + 1
          iys = ii(2) - (ixyz-1) + 1
          izs = ii(3) - (ixyz-1) + 1
          if (ixs <= 0) ixs = ixs + ngrid(1)
          if (iys <= 0) iys = iys + ngrid(2)
          if (izs <= 0) izs = izs + ngrid(3)

          if (ixs >= x_start .and. ixs <= x_end) then
            nix = nix + 1
            vx_tmp(nix) = ixs - x_start + 1
            bsc_tmp(1,nix)  = bsc(ixyz,1,i,icel)
            bscd_tmp(1,nix) = bscd(ixyz,1,i,icel)
          end if

          if (iys >= y_start .and. iys <= y_end) then
            niy = niy + 1
            vy_tmp(niy) = iys - y_start + 1
            bsc_tmp(2,niy)  = bsc(ixyz,2,i,icel)
            bscd_tmp(2,niy) = bscd(ixyz,2,i,icel)
          end if

          if (izs >= z_start .and. izs <= z_end) then
            niz = niz + 1
            vz_tmp(niz) = izs - z_start + 1
            bsc_tmp(3,niz)  = bsc(ixyz,3,i,icel)
            bscd_tmp(3,niz) = bscd(ixyz,3,i,icel)
          end if
        end do
       
!ocl nosimd 
        do iz = 1, niz
          izs = vz_tmp(iz)
          do iy = 1, niy
            iys = vy_tmp(iy)
            do ix = 1, nix
              ixs = vx_tmp(ix)

              force_local(1) = force_local(1)            &
                  + bscd_tmp(1,ix)*bsc_tmp(2,iy)         &
                  * bsc_tmp(3,iz)*qdf(ixs,iys,izs,1)
              force_local(2) = force_local(2)            &
                  + bsc_tmp(1,ix)*bscd_tmp(2,iy)         &
                  * bsc_tmp(3,iz)*qdf(ixs,iys,izs,1)
              force_local(3) = force_local(3)            &
                  + bsc_tmp(1,ix)*bsc_tmp(2,iy)          &
                  * bscd_tmp(3,iz)*qdf(ixs,iys,izs,1)
            end do
          end do
        end do

        force(1:3,i,icel) = - force_local(1:3) * &
                              qtmp*r_scale(1:3)*vol_fact4*bs_fact3d

      end do
    end do

    deallocate(work1, work2)

    !$omp end parallel

    return

  end subroutine pme_recip_noopt_1dalltoall

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_recip_noopt_1dalltoall_lj
  !> @brief        Calculate PME reciprocal part (Elec+LJ)
  !! @authors      JJ, NT
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy function information
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_recip_noopt_1dalltoall_lj(domain, enefunc, force, virial, &
                                           eelec, evdw)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wip),               intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw (nthread)

    ! local variables
    real(wp)                 :: elec_temp, evdw_temp
    real(wp)                 :: vr(3), dv(3), grid(3), half_grid(3)
    real(wp)                 :: bsc_tmp(3,n_bs), bscd_tmp(3,n_bs)
    real(wp)                 :: vxx, vyy, vzz, temp
    real(wp)                 :: efft_real, efft_imag
    real(wp)                 :: vfft_real, vfft_imag
    real(wp)                 :: tqq, tll, qtmp, ltmp
    real(wp)                 :: force_local(3), force_local_lj(3)
    real(wp)                 :: vol_factor
    integer                  :: iatmcls
    integer                  :: i, k, k1, icel, iproc, id
    integer                  :: ix, iv, iy, iz, ixs, iys, izs
    integer                  :: ixyz, nix, nix1, niy, niz, nizx, nizy
    integer                  :: ixx, iyy, izz, ii(3), is(n_bs,3)
    integer                  :: start(3), end(3)
    integer                  :: vx_tmp(n_bs), vy_tmp(n_bs), vz_tmp(n_bs)
    integer                  :: ncell, ncell_local, omp_get_thread_num
    integer                  :: kk, iorg, k_f, k_t
    integer                  :: iprocx, iprocy
    integer                  :: ix_start, ix_end, iz_start, iz_end

    complex(wp), allocatable :: work1(:), work2(:), work3(:), work4(:)
    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: nonb_lj6_factor(:)
    integer,         pointer :: natom(:)
    integer,         pointer :: atmcls(:,:)


    coord           => domain%coord
    charge          => domain%charge
    natom           => domain%num_atom
    atmcls          => domain%atom_cls_no
    nonb_lj6_factor => enefunc%nonb_lj6_factor

    ncell  = domain%num_cell_local + domain%num_cell_boundary
    ncell_local  = domain%num_cell_local 

    ! Initializing the energy and force
    !
    qdf_real(1:nlocalx*nlocaly*nlocalz) = 0.0_wp
    ldf_real(1:nlocalx*nlocaly*nlocalz) = 0.0_wp

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, iv, ix, iy, iz, icel, k, ii, dv, izz, izs, iyy, iys,  &
    !$omp         iz_start, iz_end, ix_start, ix_end, iorg, k_f, k_t,          &
    !$omp         ixx, ixs, vxx, vyy, vzz, iproc, tqq, tll, k1, temp, vr,      &
    !$omp         work1, work2, work3, work4, elec_temp, force_local, ixyz,    &
    !$omp         vx_tmp, vy_tmp, vz_tmp, nix, nix1, niy, niz, bsc_tmp,        &
    !$omp         bscd_tmp, kk, iprocx, iprocy, qtmp, grid, half_grid, is,     &
    !$omp         start, end, nizx, nizy, iatmcls, ltmp, efft_real, efft_imag, &
    !$omp         vfft_real, vfft_imag, force_local_lj, evdw_temp)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    ! initialization
    !
    allocate(work1(ngridmax), work2(2*ngridmax), &
             work3(ngridmax), work4(2*ngridmax))

    start(1) = x_start
    start(2) = y_start
    start(3) = z_start
    end(1) = x_end
    end(2) = y_end
    end(3) = z_end

    ! Initialization of Q-fanction
    !
    do iz = 1, nlocalz
      do iy = 1, nlocaly
        do ix = 1, nlocalx
          qdf(ix,iy,iz,id+1) = 0.0_wp
          ldf(ix,iy,iz,id+1) = 0.0_wp
        end do
      end do
    end do

    ! Calculateing Q-fanction
    !
    grid(1:3) = real(ngrid(1:3),wp)
    half_grid(1:3) = real(ngrid(1:3)/2,wp)

    do icel = id+1, ncell_local, nthread
      do i = 1, natom(icel)

        qtmp = charge(i,icel)
        iatmcls = atmcls(i,icel)
        ltmp = nonb_lj6_factor(iatmcls)
        vr(1:3) = coord(1:3,i,icel) * r_scale(1:3)
        vr(1:3) = vr(1:3) + half_grid(1:3) - grid(1:3)*anint(vr(1:3)/grid(1:3))

        do k = 1, 3
          vr(k) = min(vr(k),grid(k)-0.000001_wp)
          ii(k) = min(int(vr(k)),ngrid(k)-1)
        end do

        nix = 0
        niy = 0
        niz = 0
        do ixyz = 1, n_bs
          ixs = ii(1) - ixyz + 2
          if (ixs <= 0) ixs = ixs + ngrid(1)
          if (ixs >= x_start .and. ixs <= x_end) then
            nix = nix + 1
            vx_tmp(nix) = ixs - x_start + 1
            is(nix,1) = ixyz
          end if
        end do

        do ixyz = 1, n_bs
          iys = ii(2) - ixyz + 2
          if (iys <= 0) iys = iys + ngrid(2)
          if (iys >= y_start .and. iys <= y_end) then
            niy = niy + 1
            vy_tmp(niy) = iys - y_start + 1
            is(niy,2) = ixyz
          end if
        end do

        do ixyz = 1, n_bs
          izs = ii(3) - ixyz + 2
          if (izs <= 0) izs = izs + ngrid(3)
          if (izs >= z_start .and. izs <= z_end) then
            niz = niz + 1
            vz_tmp(niz) = izs - z_start + 1
            is(niz,3) = ixyz
          end if
        end do

        force(1,i,icel) = 1.0_wp

        do k = 1, 3
          vi(k,i,icel) = ii(k)
          dv(k) = vr(k) - real(ii(k),wp)
          call b_spline_dev_coef(n_bs, dv(k), bsc(1,k,i,icel), bscd(1,k,i,icel))
        end do

        do ix = 1, nix
          ixyz = is(ix,1)
          bsc_tmp(1,ix) = bsc(ixyz,1,i,icel)
        end do
        do iy = 1, niy
          ixyz = is(iy,2)
          bsc_tmp(2,iy) = bsc(ixyz,2,i,icel)
        end do
        do iz = 1, niz
          ixyz = is(iz,3)
          bsc_tmp(3,iz) = bsc(ixyz,3,i,icel)
        end do

        do iz = 1, niz
          izs = vz_tmp(iz)
          do iy = 1, niy
            iys = vy_tmp(iy)
!ocl norecurrence(qdf)
!ocl norecurrence(ldf)
             do ix = 1, nix
               ixs = vx_tmp(ix)
              qdf(ixs,iys,izs,id+1) = qdf(ixs,iys,izs,id+1)           &
                                    + bsc_tmp(1,ix)*bsc_tmp(2,iy)     &
                                     * bsc_tmp(3,iz)*qtmp   &
                                     * bs_fact3
              ldf(ixs,iys,izs,id+1) = ldf(ixs,iys,izs,id+1)           &
                                    + bsc_tmp(1,ix)*bsc_tmp(2,iy)     &
                                     * bsc_tmp(3,iz)*ltmp   &
                                     * bs_fact3
            end do
          end do
        end do
      end do
    end do

    do icel = id+ncell_local+1, ncell, nthread
      do i = 1, natom(icel)

        qtmp = charge(i,icel)
        iatmcls = atmcls(i,icel)
        ltmp = nonb_lj6_factor(iatmcls)
        vr(1:3) = coord(1:3,i,icel) * r_scale(1:3)
        vr(1:3) = vr(1:3) + half_grid(1:3) - grid(1:3)*anint(vr(1:3)/grid(1:3))

        do k = 1, 3
          vr(k) = min(vr(k),grid(k)-0.000001_wp)
          ii(k) = min(int(vr(k)),ngrid(k)-1)
        end do

        nix = 0
        niy = 0
        niz = 0
 
        ixs = ii(1) - n_bs + 2
        if (ixs <= 0) ixs = ixs + ngrid(1)
        if (ixs >= x_start .and. ixs <= x_end) then
          nix = nix + 1
          vx_tmp(nix) = ixs - x_start + 1
          is (nix,1) = n_bs
        end if
        ixs = ii(1) + 1
        if (ixs <= 0) ixs = ixs + ngrid(1)
        if (ixs >= x_start .and. ixs <= x_end) then
          nix = nix + 1
          vx_tmp(nix) = ixs - x_start + 1
          is (nix,1) = 1
        end if
        if (nix == 0) cycle

        iys = ii(2) - n_bs + 2
        if (iys <= 0) iys = iys + ngrid(2)
        if (iys >= y_start .and. iys <= y_end) then
          niy = niy + 1
          vy_tmp(niy) = iys - y_start + 1
          is (niy,2) = n_bs
        end if
        iys = ii(2) + 1
        if (iys <= 0) iys = iys + ngrid(2)
        if (iys >= y_start .and. iys <= y_end) then
          niy = niy + 1
          vy_tmp(niy) = iys - y_start + 1
          is (niy,2) = 1
        end if
        if (niy == 0) cycle

        izs = ii(3) - n_bs + 2
        if (izs <= 0) izs = izs + ngrid(3)
        if (izs >= z_start .and. izs <= z_end) then
          niz = niz + 1
          vz_tmp(niz) = izs - z_start + 1
          is (niz,3) = n_bs
        end if
        izs = ii(3) + 1
        if (izs <= 0) izs = izs + ngrid(3)
        if (izs >= z_start .and. izs <= z_end) then
          niz = niz + 1
          vz_tmp(niz) = izs - z_start + 1
          is (niz,3) = 1
        end if
        if (niz == 0) cycle

        do ixyz = 2, n_bs-1
          ixs = ii(1) - ixyz + 2
          if (ixs <= 0) ixs = ixs + ngrid(1)
          if (ixs >= x_start .and. ixs <= x_end) then
            nix = nix + 1
            vx_tmp(nix) = ixs - x_start + 1
            is(nix,1) = ixyz
          end if
        end do

        do ixyz = 2, n_bs-1
          iys = ii(2) - ixyz + 2
          if (iys <= 0) iys = iys + ngrid(2)
          if (iys >= y_start .and. iys <= y_end) then
            niy = niy + 1
            vy_tmp(niy) = iys - y_start + 1
            is(niy,2) = ixyz
          end if
        end do

        do ixyz = 2, n_bs-1
          izs = ii(3) - ixyz + 2
          if (izs <= 0) izs = izs + ngrid(3)
          if (izs >= z_start .and. izs <= z_end) then
            niz = niz + 1
            vz_tmp(niz) = izs - z_start + 1
            is(niz,3) = ixyz
          end if
        end do

        force(1,i,icel) = 1.0_wp

        do k = 1, 3
          vi(k,i,icel) = ii(k)
          dv(k) = vr(k) - real(ii(k),wp)
          call b_spline_dev_coef(n_bs, dv(k), bsc(1,k,i,icel), bscd(1,k,i,icel))
        end do

        do ix = 1, nix
          ixyz = is(ix,1)
          bsc_tmp(1,ix) = bsc(ixyz,1,i,icel)
        end do
        do iy = 1, niy
          ixyz = is(iy,2)
          bsc_tmp(2,iy) = bsc(ixyz,2,i,icel)
        end do
        do iz = 1, niz
          ixyz = is(iz,3)
          bsc_tmp(3,iz) = bsc(ixyz,3,i,icel)
        end do

        do iz = 1, niz
          izs = vz_tmp(iz)
          do iy = 1, niy
            iys = vy_tmp(iy)
!ocl norecurrence(qdf)
!ocl norecurrence(ldf)
             do ix = 1, nix
               ixs = vx_tmp(ix)
              qdf(ixs,iys,izs,id+1) = qdf(ixs,iys,izs,id+1)           &
                                    + bsc_tmp(1,ix)*bsc_tmp(2,iy)     &
                                     * bsc_tmp(3,iz)*qtmp   &
                                     * bs_fact3
              ldf(ixs,iys,izs,id+1) = ldf(ixs,iys,izs,id+1)           &
                                    + bsc_tmp(1,ix)*bsc_tmp(2,iy)     &
                                     * bsc_tmp(3,iz)*ltmp   &
                                     * bs_fact3
            end do
          end do
        end do
      end do
    end do

    !$omp barrier
    do iproc = 1, nthread
      !$omp do
      do iz = 1, nlocalz
        do iy = 1, nlocaly
          k = (iy-1)*nlocalx + (iz-1)*nlocalx*nlocaly
          do ix = 1, nlocalx
            qdf_real(k+ix) = qdf_real(k+ix) + qdf(ix,iy,iz,iproc)
            ldf_real(k+ix) = ldf_real(k+ix) + ldf(ix,iy,iz,iproc)
          end do
        end do
      end do
    end do

    ! FFT (forward)
    !
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
          vxx = vxx - tll * (theta_lj(iz,iy,ix) + vir_fact_lj(iz,iy,ix)*gx(ix)*gx(ix))
          vyy = vyy - tll * (theta_lj(iz,iy,ix) + vir_fact_lj(iz,iy,ix)*gy(iy)*gy(iy))
          vzz = vzz - tll * (theta_lj(iz,iy,ix) + vir_fact_lj(iz,iy,ix)*gz(iz)*gz(iz))

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
          vxx = vxx + tll * (theta_lj(iz,iy,1) + vir_fact_lj(iz,iy,1)*gx(1)*gx(1))
          vyy = vyy + tll * (theta_lj(iz,iy,1) + vir_fact_lj(iz,iy,1)*gy(iy)*gy(iy))
          vzz = vzz + tll * (theta_lj(iz,iy,1) + vir_fact_lj(iz,iy,1)*gz(iz)*gz(iz))

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
          vxx = vxx - tll * (theta_lj(iz,iy,ix) + vir_fact_lj(iz,iy,ix)*gx(ix)*gx(ix))
          vyy = vyy - tll * (theta_lj(iz,iy,ix) + vir_fact_lj(iz,iy,ix)*gy(iy)*gy(iy))
          vzz = vzz - tll * (theta_lj(iz,iy,ix) + vir_fact_lj(iz,iy,ix)*gz(iz)*gz(iz))

        end do

        eelec(id+1) = eelec(id+1) + elec_temp
        evdw (id+1) = evdw (id+1) + evdw_temp
        virial(1,1,id+1) = virial(1,1,id+1) + vxx
        virial(2,2,id+1) = virial(2,2,id+1) + vyy
        virial(3,3,id+1) = virial(3,3,id+1) + vzz

      end do
    end if

    ! add self term
    !
    !$omp master
    if (my_x_rank == 0 .and. my_y_rank == 0 .and. my_z_rank == 0) then
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
          ftqdf_work(k+iz,1) = ftqdf_work2(iz,ix,iy)   &
                             * cmplx(theta(iz,iy,ix),0.0_wp)
          ftldf_work(k+iz,1) = ftldf_work2(iz,ix,iy)   &
                             * cmplx(theta_lj(iz,iy,ix),0.0_wp)
        end do
      end do
    end do

    ! FFT(backward)
    !
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

    !$omp barrier

    do iz = id+1, nlocalz, nthread
      do iy = 1, nlocaly
        k = (iy-1)*nlocalx + (iz-1)*nlocalx*nlocaly
        do ix = 1, nlocalx
          qdf(ix,iy,iz,1) = qdf_real(k+ix)
          ldf(ix,iy,iz,1) = ldf_real(k+ix)
        end do
      end do
    end do

    !$omp barrier

    do icel = id+1, ncell_local, nthread
      do i = 1, natom(icel)

        force_local(1:3) = 0.0_wp
        force_local_lj(1:3) = 0.0_wp
        qtmp = charge(i,icel)
        iatmcls = atmcls(i,icel)
        ltmp = nonb_lj6_factor(iatmcls)

        do k = 1, 3
          ii(k) = vi(k,i,icel)
        end do

        nix=0
        niy=0
        niz=0

        do ixyz = 1, n_bs

          ixs = ii(1) - (ixyz-1) + 1
          iys = ii(2) - (ixyz-1) + 1
          izs = ii(3) - (ixyz-1) + 1
          if (ixs <= 0) ixs = ixs + ngrid(1)
          if (iys <= 0) iys = iys + ngrid(2)
          if (izs <= 0) izs = izs + ngrid(3)

          if (ixs >= x_start .and. ixs <= x_end) then
            nix = nix + 1
            vx_tmp(nix) = ixs - x_start + 1
            bsc_tmp(1,nix)  = bsc(ixyz,1,i,icel)
            bscd_tmp(1,nix) = bscd(ixyz,1,i,icel)
          end if

          if (iys >= y_start .and. iys <= y_end) then
            niy = niy + 1
            vy_tmp(niy) = iys - y_start + 1
            bsc_tmp(2,niy)  = bsc(ixyz,2,i,icel)
            bscd_tmp(2,niy) = bscd(ixyz,2,i,icel)
          end if

          if (izs >= z_start .and. izs <= z_end) then
            niz = niz + 1
            vz_tmp(niz) = izs - z_start + 1
            bsc_tmp(3,niz)  = bsc(ixyz,3,i,icel)
            bscd_tmp(3,niz) = bscd(ixyz,3,i,icel)
          end if
        end do
        
        do iz = 1, niz
          izs = vz_tmp(iz)
          do iy = 1, niy
            iys = vy_tmp(iy)
            do ix = 1, nix
              ixs = vx_tmp(ix)

              force_local(1) = force_local(1)            &
                  + bscd_tmp(1,ix)*bsc_tmp(2,iy)         &
                  * bsc_tmp(3,iz)*qdf(ixs,iys,izs,1)
              force_local(2) = force_local(2)            &
                  + bsc_tmp(1,ix)*bscd_tmp(2,iy)         &
                  * bsc_tmp(3,iz)*qdf(ixs,iys,izs,1)
              force_local(3) = force_local(3)            &
                  + bsc_tmp(1,ix)*bsc_tmp(2,iy)          &
                  * bscd_tmp(3,iz)*qdf(ixs,iys,izs,1)
              force_local_lj(1) = force_local_lj(1)      &
                  + bscd_tmp(1,ix)*bsc_tmp(2,iy)         &
                  * bsc_tmp(3,iz)*ldf(ixs,iys,izs,1)
              force_local_lj(2) = force_local_lj(2)      &
                  + bsc_tmp(1,ix)*bscd_tmp(2,iy)         &
                  * bsc_tmp(3,iz)*ldf(ixs,iys,izs,1)
              force_local_lj(3) = force_local_lj(3)      &
                  + bsc_tmp(1,ix)*bsc_tmp(2,iy)          &
                  * bscd_tmp(3,iz)*ldf(ixs,iys,izs,1)
            end do
          end do
        end do

        force(1:3,i,icel) =                                          &
             -force_local(1:3)*qtmp*r_scale(1:3)*vol_fact4*bs_fact3d &
             +force_local_lj(1:3)*ltmp*r_scale(1:3)*vol_fact4_lj*bs_fact3d
                            

      end do
    end do

    do icel = id+ncell_local+1, ncell, nthread
      do i = 1, natom(icel)

        if (abs(force(1,i,icel)) < EPS) cycle

        force_local(1:3) = 0.0_wp
        force_local_lj(1:3) = 0.0_wp
        qtmp = charge(i,icel)
        iatmcls = atmcls(i,icel)
        ltmp = nonb_lj6_factor(iatmcls)

        do k = 1, 3
          ii(k) = vi(k,i,icel)
        end do

        nix=0
        niy=0
        niz=0

        do ixyz = 1, n_bs

          ixs = ii(1) - (ixyz-1) + 1
          iys = ii(2) - (ixyz-1) + 1
          izs = ii(3) - (ixyz-1) + 1
          if (ixs <= 0) ixs = ixs + ngrid(1)
          if (iys <= 0) iys = iys + ngrid(2)
          if (izs <= 0) izs = izs + ngrid(3)

          if (ixs >= x_start .and. ixs <= x_end) then
            nix = nix + 1
            vx_tmp(nix) = ixs - x_start + 1
            bsc_tmp(1,nix)  = bsc(ixyz,1,i,icel)
            bscd_tmp(1,nix) = bscd(ixyz,1,i,icel)
          end if

          if (iys >= y_start .and. iys <= y_end) then
            niy = niy + 1
            vy_tmp(niy) = iys - y_start + 1
            bsc_tmp(2,niy)  = bsc(ixyz,2,i,icel)
            bscd_tmp(2,niy) = bscd(ixyz,2,i,icel)
          end if

          if (izs >= z_start .and. izs <= z_end) then
            niz = niz + 1
            vz_tmp(niz) = izs - z_start + 1
            bsc_tmp(3,niz)  = bsc(ixyz,3,i,icel)
            bscd_tmp(3,niz) = bscd(ixyz,3,i,icel)
          end if
        end do
        
        do iz = 1, niz
          izs = vz_tmp(iz)
          do iy = 1, niy
            iys = vy_tmp(iy)
            do ix = 1, nix
              ixs = vx_tmp(ix)

              force_local(1) = force_local(1)            &
                  + bscd_tmp(1,ix)*bsc_tmp(2,iy)         &
                  * bsc_tmp(3,iz)*qdf(ixs,iys,izs,1)
              force_local(2) = force_local(2)            &
                  + bsc_tmp(1,ix)*bscd_tmp(2,iy)         &
                  * bsc_tmp(3,iz)*qdf(ixs,iys,izs,1)
              force_local(3) = force_local(3)            &
                  + bsc_tmp(1,ix)*bsc_tmp(2,iy)          &
                  * bscd_tmp(3,iz)*qdf(ixs,iys,izs,1)
              force_local_lj(1) = force_local_lj(1)      &
                  + bscd_tmp(1,ix)*bsc_tmp(2,iy)         &
                  * bsc_tmp(3,iz)*ldf(ixs,iys,izs,1)
              force_local_lj(2) = force_local_lj(2)      &
                  + bsc_tmp(1,ix)*bscd_tmp(2,iy)         &
                  * bsc_tmp(3,iz)*ldf(ixs,iys,izs,1)
              force_local_lj(3) = force_local_lj(3)      &
                  + bsc_tmp(1,ix)*bsc_tmp(2,iy)          &
                  * bscd_tmp(3,iz)*ldf(ixs,iys,izs,1)
            end do
          end do
        end do

        force(1:3,i,icel) = &
            -force_local(1:3)*qtmp*r_scale(1:3)*vol_fact4*bs_fact3d &
            +force_local_lj(1:3)*ltmp*r_scale(1:3)*vol_fact4_lj*bs_fact3d

      end do
    end do

    deallocate(work1, work2, work3, work4)

    !$omp end parallel

    return

  end subroutine pme_recip_noopt_1dalltoall_lj

end module sp_energy_pme_noopt_1dalltoall_mod
