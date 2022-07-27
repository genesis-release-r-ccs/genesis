!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_energy_pme_mod
!> @brief   Smooth particle mesh ewald method
!! @authors Takashi Imai (TI), Jaewoon Jung (JJ), Takaharu Mori (TM), 
!!          Chigusa Kobayashi (CK), Motoshi Kamiya (MK), Yuji Sugita (YS)
!! @note    TI modified for NPT (2010/11/23)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_energy_pme_mod

  use at_boundary_str_mod
  use at_pairlist_str_mod
  use at_enefunc_str_mod
  use at_boundary_mod
  use molecules_str_mod
  use math_libs_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use fft3d_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! variables
  integer, public, parameter :: PmeMaxNspline = 10

  real(wp),          allocatable :: myqdf(:,:,:,:)

  ! for FFTW
  integer,                  save :: nlocalx1
  integer,                  save :: ngridmax
  integer,                  save :: maxproc

#ifdef FFTW
  ! FFTW plans and work arrays
  integer*8,                save :: plan_fx
  integer*8,                save :: plan_fy
  integer*8,                save :: plan_fz
  integer*8,                save :: plan_bz
  integer*8,                save :: plan_by
  integer*8,                save :: plan_bx

  complex(wp), save, allocatable :: ftqdf(:)
  complex(wp), save, allocatable :: ftqdf_work(:,:)

  complex(wp), save, allocatable :: ftqdf2(:)
  complex(wp), save, allocatable :: ftqdf3(:)
  complex(wp), save, allocatable :: ftqdf_work2(:,:,:)
#endif

  ! subroutines
  public  :: pme_pre
  public  :: pme_recip
  public  :: setup_pme
  public  :: dealloc_pme

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_pre
  !> @brief        Prepare functions for PME calculation
  !! @authors      TI, JJ, TM
  !! @param[in]    boundary : boundary conditions information
  !! @note         Extracted from setup_pme for NPT calculation
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_pre(boundary, enefunc)

    ! formal arguments
    type(s_boundary),         intent(in)    :: boundary
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: fact, gfact(3), g2, alpha2m
    integer                  :: i, j, k, is, js, ks
    integer                  :: ngrid(3)

    real(wp), pointer        :: gx(:), gy(:), gz(:), b2(:,:)
    real(wp), pointer        :: vir_fact(:,:,:), theta(:,:,:)


    ! use pointers
    !
    gx         => enefunc%pme%gx
    gy         => enefunc%pme%gy
    gz         => enefunc%pme%gz
    vir_fact   => enefunc%pme%vir_fact
    theta      => enefunc%pme%theta
    b2         => enefunc%pme%b2
    ngrid(1:3) =  enefunc%pme%ngrid(1:3)
    alpha2m    =  enefunc%pme%alpha2m

    ! Prepareing theta = F^-1[theta](h), h shifted
    ! Gx, Gy, Gz, and vir_fact=-2(1+G^2/4a)/G^2 are also prepared 
    ! for virial calculation
    !
    gfact(1) = 2.0_wp*PI/boundary%box_size_x
    gfact(2) = 2.0_wp*PI/boundary%box_size_y
    gfact(3) = 2.0_wp*PI/boundary%box_size_z

    fact = 0.25_wp / alpha2m

    !$omp parallel do default(none)                                   &
    !$omp private(i, j, k, g2, ks, is, js)                            &
    !$omp shared(ngrid, gfact, gx, gy, gz, vir_fact, theta, fact, b2)
    !
    do k = 1, ngrid(3)/2 + 1
      ks    = k - 1
      gz(k) = gfact(3) * real(ks,wp)
      do j = 1, ngrid(2)/2 + 1
        js    = j - 1
        gy(j) = gfact(2) * real(js,wp)
        do i =1, ngrid(1)/2 + 1
          is    = i - 1
          gx(i) = gfact(1) * real(is,wp)
          g2    = gx(i)**2 + gy(j)**2 + gz(k)**2
          if (g2 > EPS) then
            vir_fact(i,j,k) = -2.0_wp * (1.0_wp - g2 * fact)/g2
            if (g2*fact < -80.0_wp) then
              theta(i,j,k)    = 0.0_wp
            else
              theta(i,j,k)    = b2(i,1) * b2(j,2) * b2(k,3) * exp(g2 * fact)/g2
            endif
            if (abs(theta(i,j,k)) < 1.0e-15_wp) theta(i,j,k) = 0.0_wp
          else
            vir_fact(i,j,k) = 0.0_wp
            theta(i,j,k)    = 0.0_wp
          end if
        end do
      end do

      do j = ngrid(2)/2 + 2, ngrid(2)
        js    = j - 1 - ngrid(2)
        gy(j) = gfact(2) * real(js,wp)
        do i = 1, ngrid(1)/2 + 1
          is    = i - 1
          gx(i) = gfact(1) * real(is,wp)
          g2    = gx(i)**2 + gy(j)**2 + gz(k)**2
          if (g2 > EPS) then
            vir_fact(i,j,k) = -2.0_wp * (1.0_wp - g2 * fact)/g2
            if (g2*fact < -80.0_wp) then
              theta(i,j,k)    = 0.0_wp
            else
              theta(i,j,k)    = b2(i,1) * b2(j,2) * b2(k,3) * exp(g2 * fact)/g2
            end if
            if (abs(theta(i,j,k)) < 1.0e-15_wp) theta(i,j,k) = 0.0_wp
          else
            vir_fact(i,j,k) = 0.0_wp
            theta(i,j,k)    = 0.0_wp
          end if
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do default(none)                                   &
    !$omp private(i, j, k, g2, ks, is, js)                            &
    !$omp shared(ngrid, gfact, gx, gy, gz, vir_fact, theta, fact, b2)
    !
    do k = ngrid(3)/2 + 2, ngrid(3)
      ks    = k - 1 - ngrid(3)
      gz(k) = gfact(3) * real(ks,wp)
      do j = 1, ngrid(2)/2 + 1
        js    = j - 1
        gy(j) = gfact(2) * real(js,wp)
        do i = 1, ngrid(1)/2 + 1
          is    = i - 1
          gx(i) = gfact(1) * real(is,wp)
          g2    = gx(i)**2 + gy(j)**2 + gz(k)**2
          if (g2 > EPS) then
            vir_fact(i,j,k) = -2.0_wp * (1.0_wp - g2 * fact)/g2
            if (g2*fact < -80.0_wp) then
              theta(i,j,k)    = 0.0_wp
            else
              theta(i,j,k)    = b2(i,1) * b2(j,2) * b2(k,3) * exp(g2 * fact)/g2
            end if
            if (abs(theta(i,j,k)) < 1.0e-15_wp) theta(i,j,k) = 0.0_wp
          else
            vir_fact(i,j,k) = 0.0_wp
            theta(i,j,k)    = 0.0_wp
          end if
        end do
      end do

      do j = ngrid(2)/2 + 2, ngrid(2)
        js    = j - 1 - ngrid(2)
        gy(j) = gfact(2) * real(js,wp)
        do i = 1, ngrid(1)/2 + 1
          is    = i - 1
          gx(i) = gfact(1) * real(is,wp)
          g2    = gx(i)**2 + gy(j)**2 + gz(k)**2
          if (g2 > EPS) then
            vir_fact(i,j,k) = -2.0_wp * (1.0_wp - g2 * fact)/g2
            if (g2*fact < -80.0_wp) then
              theta(i,j,k)    = 0.0_wp
            else
              theta(i,j,k)    = b2(i,1) * b2(j,2) * b2(k,3) * exp(g2 * fact)/g2
            end if
            if (abs(theta(i,j,k)) < 1.0e-15_wp) theta(i,j,k) = 0.0_wp
          else
            vir_fact(i,j,k) = 0.0_wp
            theta(i,j,k)    = 0.0_wp
          end if
        end do
      end do
    end do
    !$omp end parallel do

    do k = 1, ngrid(3)
      do j = 1, ngrid(2)
        do i = ngrid(1)/2+2, ngrid(1)
          theta(i,j,k)    = theta(ngrid(1)-i+2,j,k)
          vir_fact(i,j,k) = vir_fact(ngrid(1)-i+2,j,k)
        end do
      end do
    end do

    do i = ngrid(1)/2 + 2, ngrid(1)
      gx(i) = gx(ngrid(1)-i+2)
    end do

    vir_fact(1,1,1) = 0.0_wp
    theta(1,1,1)    = 0.0_wp

    return

  end subroutine pme_pre

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_recip
  !> @brief        Calculate PME reciprocal part
  !! @authors      TI, JJ, CK, TM, MK
  !! @param[in]    molecule : molecule information
  !! @param[in]    r        : coordinates 
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_recip(enefunc, boundary, molecule, r, force, virial, eelec)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary),         intent(in)    :: boundary
    type(s_molecule), target, intent(in)    :: molecule
    real(wp),                 intent(in)    :: r(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eelec

    ! local variables
    real(wp)                  :: u, el_fact, origin(3), box_inv(3), r_scale(3)
    real(wp)                  :: dv(3), vol_fact2, vol_fact4
    real(wp)                  :: vxx,vyx,vzx,vyy,vzy,vzz
    real(wp)                  :: tqq
    real(wp)                  :: bsc(PmeMaxNspline,3), bscd(PmeMaxNspline,3)
    real(wp)                  :: bs_fact3, bs_fact3d
    real(wp)                  :: vir_local(3,3)
    integer                   :: i, k
    integer                   :: ix, iy, iz, ixs, iys, izs, ixx, iyy, izz, ii(3)
    integer                   :: iproc, jproc
    integer                   :: istart_atom, iend_atom
    integer                   :: ngrid(3), n_bs
    integer                   :: niy, iyss
    integer                   :: id, omp_get_thread_num

#ifdef FFTW
    complex(wp), allocatable  :: work1(:), work2(:)
#endif

    real(wp),         pointer :: q(:), f(:,:), v(:,:), gx(:), gy(:), gz(:)
    real(wp),         pointer :: vir_fact(:,:,:), theta(:,:,:)
    real(wp),         pointer :: qdf(:,:,:), qdf_recv(:,:,:,:) 
    real(wp),         pointer :: qdf_send(:,:,:,:)
    integer,          pointer :: natm
    integer,          pointer :: x_start1(:), x_end1(:)
    integer,          pointer :: y_start(:), y_end(:)
    integer,          pointer :: y_local(:), y_start1(:), y_end1(:)
    integer,          pointer :: z_start(:), z_end(:), z_local(:)
    complex(wp),      pointer :: ftqdf_localA(:,:,:), ftqdf_localB(:,:,:)


    ! use pointers
    !
    natm         => molecule%num_atoms
    q            => molecule%charge
    f            => enefunc%pme%f
    v            => enefunc%pme%v
    gx           => enefunc%pme%gx
    gy           => enefunc%pme%gy
    gz           => enefunc%pme%gz
    qdf          => enefunc%pme%qdf
    x_start1     => enefunc%pme%x_start1
    x_end1       => enefunc%pme%x_end1
    y_start      => enefunc%pme%y_start
    y_end        => enefunc%pme%y_end
    y_local      => enefunc%pme%y_local
    y_start1     => enefunc%pme%y_start1
    y_end1       => enefunc%pme%y_end1
    z_start      => enefunc%pme%z_start
    z_end        => enefunc%pme%z_end
    z_local      => enefunc%pme%z_local
    qdf_recv     => enefunc%pme%qdf_recv
    qdf_send     => enefunc%pme%qdf_send
    ftqdf_localA => enefunc%pme%ftqdf_localA
    ftqdf_localB => enefunc%pme%ftqdf_localB
    theta        => enefunc%pme%theta
    vir_fact     => enefunc%pme%vir_fact

    istart_atom  =  enefunc%pme%istart_atom
    iend_atom    =  enefunc%pme%iend_atom
    ngrid(1:3)   =  enefunc%pme%ngrid(1:3)
    bs_fact3     =  enefunc%pme%bs_fact3
    bs_fact3d    =  enefunc%pme%bs_fact3d
    n_bs         =  enefunc%pme%n_bs
    el_fact      =  ELECOEF/enefunc%dielec_const

    origin(1)    = boundary%origin_x
    origin(2)    = boundary%origin_y
    origin(3)    = boundary%origin_z
    box_inv(1)   = 1.0_wp/boundary%box_size_x
    box_inv(2)   = 1.0_wp/boundary%box_size_y
    box_inv(3)   = 1.0_wp/boundary%box_size_z
    r_scale(1)   = real(ngrid(1),wp)*box_inv(1)
    r_scale(2)   = real(ngrid(2),wp)*box_inv(2)
    r_scale(3)   = real(ngrid(3),wp)*box_inv(3)

    vol_fact2    = 2.0_wp * PI * el_fact * box_inv(1)*box_inv(2)*box_inv(3) 
    vol_fact4    = 2.0_wp * vol_fact2

    ! Initializing the energy and force
    !
    u = 0.0_wp
    f(1:3,1:natm) = 0.0_wp
    vir_local(1:3,1:3) = 0.0_wp

    niy = y_local(my_city_rank) / nprocz

#ifdef FFTW
    !$omp parallel default(shared) &
    !$omp private(i, k, ix, iy, iz, ixs, iys, izs, ixx, iyy, izz, ii, tqq, &
    !$omp         dv, vxx, vyx, vzx, vyy, vzy, vzz, &
    !$omp         id, iyss, iproc, jproc, work1, work2, bsc, bscd)
    !
    allocate(work1(ngridmax), work2(2*ngridmax))
#else
    !$omp parallel default(shared) &
    !$omp private(i, k, ix, iy, iz, ixs, iys, izs, ixx, iyy, izz, ii, tqq, &
    !$omp         dv, id, iyss, iproc, jproc, bsc)
    !
#endif

#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif


    ! Translation to scaled fractional coordinates (v(x,y,z))
    !
    do i = istart_atom+id, iend_atom, nthread
      do k = 1, 3
        v(k,i) = (r(k,i) - origin(k)) * r_scale(k)
        v(k,i) = v(k,i) + real(ngrid(k) / 2,wp) &
                - real(ngrid(k),wp) * anint(v(k,i) / real(ngrid(k),wp))
      end do
    end do

    !$omp barrier
    !

    ! Initialization of Q-fanction
    !
    myqdf(1:ngrid(1),1:ngrid(2),1:ngrid(3),id+1) = 0.0_wp
    do iz = id+1, ngrid(3), nthread
      do iy = 1, ngrid(2)
        do ix = 1, ngrid(1)
          qdf(ix,iy,iz) = 0.0_wp
        end do
      end do
    end do

    !$omp barrier
    !

    ! Calculateing Q-fanction
    !
    !$omp do
    do i = istart_atom, iend_atom

      do k = 1, 3
        ii(k) = int(v(k,i))
        dv(k) = v(k,i) - real(ii(k),wp)
        call b_spline_coef(n_bs, dv(k), bsc(1:n_bs,k))
      end do

      do iz = 1, n_bs

        izz = ii(3) - iz + 1
        if (izz < 0) then
          izs = izz + ngrid(3) + 1
        else if (izz == ngrid(3)) then
          izs = izz
        else
          izs = izz + 1
        end if

        do iy = 1, n_bs

          iyy = ii(2) - iy + 1
          if (iyy < 0) then
            iys = iyy + ngrid(2) + 1
          else if (iyy == ngrid(2)) then
            iys = iyy
          else
            iys = iyy + 1
          end if

          do ix = 1, n_bs

            ixx = ii(1) - ix + 1
            if (ixx < 0) then
              ixs = ixx + ngrid(1) + 1
            else if (ixx == ngrid(1)) then
              ixs = ixx
            else
              ixs = ixx + 1
            end if
            myqdf(ixs,iys,izs,id+1) = myqdf(ixs,iys,izs,id+1) &
                             + q(i) * bsc(ix,1)*bsc(iy,2)*bsc(iz,3)*bs_fact3
          end do
        end do
      end do
    end do
    !$omp barrier

    ! gather thread distributed data
    do i = 1, nthread
      do iz = id+1, ngrid(3), nthread
        do iy = 1, ngrid(2)
          do ix = 1, ngrid(1)
            qdf(ix,iy,iz) = qdf(ix,iy,iz) + myqdf(ix,iy,iz,i)
          end do
        end do
      end do
    end do
    !$omp barrier

    ! Q-fanction is saved to qdf_send for communication
    !
    do jproc = 0, my_city_rank - 1
      do iz = z_start(jproc)+id, z_end(jproc), nthread
        izs = iz - z_start(jproc) + 1
        do iy = y_start(jproc), y_end(jproc)
          iys = iy - y_start(jproc) + 1
          do ix = 1, ngrid(1)
            qdf_send(ix,iys,izs,jproc+1) = qdf(ix,iy,iz)
          end do
        end do
      end do
    end do

    !$omp barrier

    do jproc = my_city_rank + 1, nproc_city - 1
      do iz = z_start(jproc)+id, z_end(jproc), nthread
        izs = iz - z_start(jproc) + 1
        do iy = y_start(jproc), y_end(jproc)
          iys = iy - y_start(jproc) + 1
          do ix = 1, ngrid(1)
            qdf_send(ix,iys,izs,jproc+1) = qdf(ix,iy,iz)
          end do
        end do
      end do
    end do


    !$omp barrier

    do iz = z_start(my_city_rank)+id, z_end(my_city_rank), nthread
      izs = iz - z_start(my_city_rank) + 1
      do iy = y_start(my_city_rank), y_end(my_city_rank)
        iys = iy - y_start(my_city_rank) + 1
        do ix = 1, ngrid(1)
          ftqdf_localA(ix,iys,izs) = qdf(ix,iy,iz)
        end do
      end do
    end do

#ifdef HAVE_MPI_GENESIS
    !$omp barrier
    !$omp master
    call mpi_alltoall &
         (qdf_send, ngrid(1)*y_local(0)*z_local(0), mpi_wp_real, &
          qdf_recv, ngrid(1)*y_local(0)*z_local(0), mpi_wp_real, &
          mpi_comm_city, ierror)
    !$omp end master
#endif

    !$omp barrier

    ! Q-fanction is accumulated to qdf_send for each processor
    !
    do iproc = 0, my_city_rank - 1
      do iz = id+1, z_local(my_city_rank), nthread
        do iy = 1, y_local(my_city_rank)
          do ix = 1, ngrid(1)
            ftqdf_localA(ix,iy,iz) = ftqdf_localA(ix,iy,iz) + &
                 qdf_recv(ix,iy,iz,iproc+1)
          end do
        end do
      end do
    end do

    !$omp barrier

    do iproc = my_city_rank + 1, nproc_city - 1
      do iz = id+1, z_local(my_city_rank), nthread
        do iy = 1, y_local(my_city_rank)
          do ix = 1, ngrid(1)
            ftqdf_localA(ix,iy,iz) = ftqdf_localA(ix,iy,iz) + &
                 qdf_recv(ix,iy,iz,iproc+1)
          end do
        end do
      end do
    end do

    ! Q -> F^-1[Q]
    !
    !$omp barrier
#ifdef FFTW
    call fft3d_1d_alltoall_x_serial(plan_fx, plan_fy, plan_fz,                 &
                                 ftqdf_localA, ftqdf, ftqdf2, ftqdf3,          &
                                 ftqdf_work, ftqdf_work, ftqdf_localB,         &
                                 work1, work2,                                 &
                                 y_local(my_city_rank), z_local(my_city_rank), &
                                 nlocalx1, ngrid(1), ngrid(2), ngrid(3),       &
                                 nprocy, nprocz, id, nthread,                  &
                                 my_y_rank, my_z_rank,                         &
                                 y_start(my_city_rank), y_end(my_city_rank),   &
                                 z_start(my_city_rank), z_end(my_city_rank),   &
                                 grid_commy, grid_commz)
#else
    !$omp end parallel
    !

    call pzfft3dv(ftqdf_localA, ftqdf_localB, ngrid(1), ngrid(2), ngrid(3), &
                  grid_commy, grid_commz, nprocy, nprocz, 0)
    call pzfft3dv(ftqdf_localA, ftqdf_localB, ngrid(1), ngrid(2), ngrid(3), &
                  grid_commy, grid_commz, nprocy, nprocz, -2)

    !$omp parallel default(shared) &
    !$omp private(i, k, ix, iy, iz, ixs, iys, izs, ixx, iyy, izz, ii, tqq, &
    !$omp         dv, vxx, vyx, vzx, vyy, vzy, vzz, &
    !$omp         id, iyss, iproc, jproc)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
#endif

    !$omp barrier
    !$omp do reduction(+:u) reduction(+:vir_local)
#ifdef FFTW
    do ix = 1, nlocalx1
      do iy = 1, niy
        iys = y_start(my_city_rank) + niy * my_z_rank - 1 + iy
        do iz = 1, ngrid(3)
          tqq = real(ftqdf_localB(iz,ix,iy),wp)**2 + &
                imag(ftqdf_localB(iz,ix,iy))**2
          tqq = tqq * theta(ix,iys,iz) * vol_fact4
#else
    do iz = 1, ngrid(3)
      do iy = y_start1(my_city_rank), y_end1(my_city_rank)

        iys = iy - y_start1(my_city_rank) + 1

        do ix = x_start1(my_city_rank), x_end1(my_city_rank)
          ixs = ix - x_start1(my_city_rank) + 1
          tqq = real(ftqdf_localB(ixs,iys,iz),wp)**2 + &
                imag(ftqdf_localB(ixs,iys,iz))**2
          tqq = tqq * theta(ix,iy,iz) * vol_fact2
#endif
          u = u + tqq

          ! virial
          !
#ifdef FFTW
          vxx = tqq * (1.0_wp + vir_fact(ix,iys,iz) * gx(ix) * gx(ix))
          vyx = tqq * vir_fact(ix,iys,iz) * gy(iys) * gx(ix)
          vzx = tqq * vir_fact(ix,iys,iz) * gz(iz) * gx(ix)
          vyy = tqq * (1.0_wp + vir_fact(ix,iys,iz) * gy(iys) * gy(iys))
          vzy = tqq * vir_fact(ix,iys,iz) * gz(iz) * gy(iys)
          vzz = tqq * (1.0_wp + vir_fact(ix,iys,iz) * gz(iz) * gz(iz))
#else
          vxx = tqq * (1.0_wp + vir_fact(ix,iy,iz) * gx(ix) * gx(ix))
          vyx = tqq * vir_fact(ix,iy,iz) * gy(iy) * gx(ix)
          vzx = tqq * vir_fact(ix,iy,iz) * gz(iz) * gx(ix)
          vyy = tqq * (1.0_wp + vir_fact(ix,iy,iz) * gy(iy) * gy(iy))
          vzy = tqq * vir_fact(ix,iy,iz) * gz(iz) * gy(iy)
          vzz = tqq * (1.0_wp + vir_fact(ix,iy,iz) * gz(iz) * gz(iz))
#endif

          vir_local(1,1) = vir_local(1,1) + vxx
          vir_local(2,1) = vir_local(2,1) + vyx
          vir_local(3,1) = vir_local(3,1) + vzx
          vir_local(1,2) = vir_local(1,2) + vyx
          vir_local(2,2) = vir_local(2,2) + vyy
          vir_local(3,2) = vir_local(3,2) + vzy
          vir_local(1,3) = vir_local(1,3) + vzx
          vir_local(2,3) = vir_local(2,3) + vzy
          vir_local(3,3) = vir_local(3,3) + vzz
        end do
      end do

#ifdef FFTW
      if (ix == 1 .or. ix == nlocalx1) then
        do iy = 1, niy
          iys = y_start(my_city_rank) + niy * my_z_rank - 1 + iy
          do iz = 1, ngrid(3)
            tqq = real(ftqdf_localB(iz,ix,iy),wp)**2 + &
                  imag(ftqdf_localB(iz,ix,iy))**2
            tqq = tqq * theta(ix,iys,iz) * vol_fact2
            u = u - tqq

            ! virial
            !
            vxx = tqq * (1.0_wp + vir_fact(ix,iys,iz) * gx(ix) * gx(ix))
            vyx = tqq * vir_fact(ix,iys,iz) * gy(iys) * gx(ix)
            vzx = tqq * vir_fact(ix,iys,iz) * gz(iz) * gx(ix)
            vyy = tqq * (1.0_wp + vir_fact(ix,iys,iz) * gy(iys) * gy(iys))
            vzy = tqq * vir_fact(ix,iys,iz) * gz(iz) * gy(iys)
            vzz = tqq * (1.0_wp + vir_fact(ix,iys,iz) * gz(iz) * gz(iz))

            vir_local(1,1) = vir_local(1,1) - vxx
            vir_local(2,1) = vir_local(2,1) - vyx
            vir_local(3,1) = vir_local(3,1) - vzx
            vir_local(1,2) = vir_local(1,2) - vyx
            vir_local(2,2) = vir_local(2,2) - vyy
            vir_local(3,2) = vir_local(3,2) - vzy
            vir_local(1,3) = vir_local(1,3) - vzx
            vir_local(2,3) = vir_local(2,3) - vzy
            vir_local(3,3) = vir_local(3,3) - vzz
          end do
        end do
      end if
#endif

    end do

    !$omp barrier

    ! F^-1[Th]*F^-1[Q] (=X)
    !
#ifdef FFTW
    do ix = 1, nlocalx1
      do iy = id+1, niy, nthread
        iys = y_start(my_city_rank) + niy * my_z_rank - 1 + iy
        do iz = 1, ngrid(3)
          ftqdf_localB(iz,ix,iy) = cmplx(theta(ix,iys,iz)) * &
               ftqdf_localB(iz,ix,iy)
        end do
      end do
    end do
#else
    do iz = 1, ngrid(3)
      do iy = y_start1(my_city_rank)+id, y_end1(my_city_rank), nthread
        iys = iy - y_start1(my_city_rank) + 1
        do ix = x_start1(my_city_rank), x_end1(my_city_rank)
          ixs = ix - x_start1(my_city_rank) + 1
          ftqdf_localB(ixs,iys,iz) = cmplx(theta(ix,iy,iz)) * &
               ftqdf_localB(ixs,iys,iz)
        end do
      end do
    end do
#endif

    ! X -> F[X] (=Conv[Th*Q])
    !
    !$omp barrier
#ifdef FFTW
    call bfft3d_1d_alltoall_x_serial(plan_bx, plan_by, plan_bz,              &
                         ftqdf_localA, ftqdf, ftqdf2,                        &
                         ftqdf_work, ftqdf_work, ftqdf_localB, work1, work2, &
                         y_local(my_city_rank), z_local(my_city_rank),       &
                         nlocalx1,                                           &
                         ngrid(1), ngrid(2), ngrid(3), nprocy, nprocz,       &
                         id, nthread, my_y_rank, my_z_rank,                  &
                         y_start(my_city_rank), y_end(my_city_rank),         &
                         z_start(my_city_rank), z_end(my_city_rank),         &
                         grid_commy, grid_commz, niy)
#else
    !$omp end parallel

    call pzfft3dv(ftqdf_localB, ftqdf_localA, ngrid(1), ngrid(2), ngrid(3),  &
                  grid_commy, grid_commz, nprocy, nprocz, 2)

    !$omp parallel default(shared) &
    !$omp private(i, k, ix, iy, iz, ixs, iys, izs, ixx, iyy, izz, ii, tqq,   &
    !$omp         dv, vxx, vyx, vzx, vyy, vzy, vzz, &
    !$omp         id, iyss, iproc, jproc, bsc, bscd)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
#endif

    do iz = id+1, z_local(my_city_rank), nthread
      do iy = 1, y_local(my_city_rank)
        do ix = 1, ngrid(1)
          qdf_send(ix,iy,iz,my_city_rank+1) = real(ftqdf_localA(ix,iy,iz),wp)
        end do
      end do
    end do

    ! X is saved on qdf_send and distributed
    !
#ifdef HAVE_MPI_GENESIS
    !$omp barrier
    !$omp master
    call mpi_allgather(qdf_send(1,1,1,my_city_rank+1), &
                       z_local(0)*y_local(0)*ngrid(1), &
                       mpi_wp_real, qdf_recv,          &
                       z_local(0)*y_local(0)*ngrid(1), &
                       mpi_wp_real, mpi_comm_city, ierror)
    !$omp end master
#else

    do iz = id+1, z_local(my_city_rank), nthread
      do iy = 1, y_local(my_city_rank)
        do ix = 1, ngrid(1)
          qdf_recv(ix,iy,iz,my_city_rank+1) = &
               qdf_send(ix,iy,iz,my_city_rank+1)
        end do
      end do
    end do
#endif
    !$omp barrier

    ! X is saved on qdf
    !
    do iproc = 0, nproc_city - 1
      do iz = z_start(iproc)+id, z_end(iproc), nthread
        izs = iz - z_start(iproc) + 1
        do iy = y_start(iproc), y_end(iproc)
          iys = iy - y_start(iproc) + 1
          do ix = 1, ngrid(1)
            qdf(ix,iy,iz) = qdf_recv(ix,iys,izs,iproc+1)
          end do
        end do
      end do
    end do

    !$omp barrier

    ! Calculating Fi=Sum{dQ/dri*Conv[Th*Q]}
    !
    !
    do i = istart_atom+id, iend_atom, nthread

      do k = 1, 3
        ii(k) = int(v(k,i))
        dv(k) = v(k,i) - real(ii(k),wp)
        call b_spline_dev_coef(n_bs, dv(k), bsc(1:n_bs,k), bscd(1:n_bs,k))
      end do

      do iz = 1, n_bs

        izz = ii(3) - iz + 1
        if (izz < 0) then
          izs = izz + ngrid(3) + 1
        else if (izz == ngrid(3)) then
          izs = izz
        else
          izs = izz + 1
        end if

        do iy = 1, n_bs

          iyy = ii(2) - iy + 1
          if (iyy < 0) then
            iys = iyy + ngrid(2) + 1
          else if (iyy == ngrid(2)) then
            iys = iyy
          else
            iys = iyy + 1
          end if

          do ix = 1, n_bs

            ixx = ii(1) - ix + 1
            if (ixx < 0) then
              ixs = ixx + ngrid(1) + 1
            else if (ixx == ngrid(1)) then
              ixs = ixx
            else
              ixs = ixx + 1
            end if

            f(1,i) = f(1,i) + bscd(ix,1)*bsc (iy,2)*bsc (iz,3)*qdf(ixs,iys,izs)
            f(2,i) = f(2,i) + bsc (ix,1)*bscd(iy,2)*bsc (iz,3)*qdf(ixs,iys,izs)
            f(3,i) = f(3,i) + bsc (ix,1)*bsc (iy,2)*bscd(iz,3)*qdf(ixs,iys,izs)

          end do
        end do
      end do

      do k = 1, 3
        f(k,i) = -f(k,i) * q(i) * r_scale(k) * vol_fact4 * bs_fact3d
      end do

    end do 

#ifdef FFTW
    deallocate(work1,work2)
#endif

    !$omp end parallel

    eelec = eelec + u
    virial(1:3,1:3) = virial(1:3,1:3) + vir_local(1:3,1:3)
    force(1:3,1:natm,1) = force(1:3,1:natm,1) + f(1:3,1:natm)

    return

  end subroutine pme_recip

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_pme
  !> @brief        Memory allocation and other preparations for FFTW/FFTE
  !! @note         Extracted from setup_pme of spdyn
  !! @authors      MK
  !! @param[inout] enefunc  : energy function
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_pme(enefunc)

#ifdef FFTW
    include "fftw3.f"
#endif

    ! formal arguments
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variables
    integer                  :: nlocaly, nlocalz, ngrid(3)
    complex(wp), allocatable :: work1(:), work2(:)
    real(wp),    allocatable :: work3(:)


    ngrid(1:3) = enefunc%pme%ngrid(1:3)

    ngridmax   = max(ngrid(1),ngrid(2),ngrid(3))
    nlocalx1   = ngrid(1)/2 + 1
    nlocaly    = enefunc%pme%y_local(0)
    nlocalz    = enefunc%pme%z_local(0)
    maxproc    = max(1,nprocy,nprocz)

#ifndef FFTW
    allocate(myqdf(ngrid(1),ngrid(2),ngrid(3),nthread))
#else
    allocate(ftqdf(nlocalx1*nlocaly*nlocalz), &
             ftqdf2(nlocalx1*nlocaly*nlocalz), &
             ftqdf3(nlocalx1*nlocaly*nlocalz), &
             ftqdf_work(nlocalx1*nlocaly*nlocalz,maxproc), &
             myqdf(ngrid(1),ngrid(2),ngrid(3),nthread))
 
    allocate(work1(ngridmax),work2(2*ngridmax),work3(ngrid(1)))

    ! FFTW_PATIENT (instead of FFTW_MEASURE) might be desirable
#ifdef _SINGLE
    call sfftw_plan_dft_1d(plan_fx,ngrid(1),work1,work2,FFTW_FORWARD,FFTW_MEASURE)
    call sfftw_plan_dft_1d(plan_fy,ngrid(2),work1,work2,FFTW_FORWARD,FFTW_MEASURE)
    call sfftw_plan_dft_1d(plan_fz,ngrid(3),work1,work2,FFTW_FORWARD,FFTW_MEASURE)
    call sfftw_plan_dft_1d(plan_bz,ngrid(3),work1,work2,FFTW_BACKWARD,FFTW_MEASURE)
    call sfftw_plan_dft_1d(plan_by,ngrid(2),work1,work2,FFTW_BACKWARD,FFTW_MEASURE)
    call sfftw_plan_dft_c2r_1d(plan_bx,ngrid(1),work1,work3,FFTW_MEASURE)
#else
    call dfftw_plan_dft_1d(plan_fx,ngrid(1),work1,work2,FFTW_FORWARD,FFTW_MEASURE)
    call dfftw_plan_dft_1d(plan_fy,ngrid(2),work1,work2,FFTW_FORWARD,FFTW_MEASURE)
    call dfftw_plan_dft_1d(plan_fz,ngrid(3),work1,work2,FFTW_FORWARD,FFTW_MEASURE)
    call dfftw_plan_dft_1d(plan_bz,ngrid(3),work1,work2,FFTW_BACKWARD,FFTW_MEASURE)
    call dfftw_plan_dft_1d(plan_by,ngrid(2),work1,work2,FFTW_BACKWARD,FFTW_MEASURE)
    call dfftw_plan_dft_c2r_1d(plan_bx,ngrid(1),work1,work3,FFTW_MEASURE)
#endif

    deallocate(work1,work2,work3)

#endif

    return

  end subroutine setup_pme

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_pme
  !> @brief        Memory deallocation
  !! @authors      MK
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_pme

    ! local variables
    integer                  :: dealloc_stat


#ifndef FFTW
    if (allocated(myqdf)) deallocate(myqdf)
#else
    if (.not. allocated(ftqdf)) return

    deallocate(ftqdf, &
               ftqdf2, &
               ftqdf3, &
               ftqdf_work, &
               myqdf, &
               stat = dealloc_stat)

    if (dealloc_stat /= 0) call error_msg_dealloc

#ifdef _SINGLE
    call sfftw_destroy_plan(plan_fx)
    call sfftw_destroy_plan(plan_fy)
    call sfftw_destroy_plan(plan_fz)
    call sfftw_destroy_plan(plan_bz)
    call sfftw_destroy_plan(plan_by)
    call sfftw_destroy_plan(plan_bx)
#else
    call dfftw_destroy_plan(plan_fx)
    call dfftw_destroy_plan(plan_fy)
    call dfftw_destroy_plan(plan_fz)
    call dfftw_destroy_plan(plan_bz)
    call dfftw_destroy_plan(plan_by)
    call dfftw_destroy_plan(plan_bx)
#endif
#endif

    return

  end subroutine dealloc_pme

end module at_energy_pme_mod
