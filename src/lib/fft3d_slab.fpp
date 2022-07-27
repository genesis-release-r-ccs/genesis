!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fft3d_slab_mod
!> @brief   routines for fast fourier transform with slab decomposition
!! @authors Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fft3d_slab_mod

  use constants_mod
  use timers_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    fft3d_slab
  !> @brief        calculate 3D fast fourier transform with 1D decomposition
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine fft3d_slab(qdf_work, fftqdfyzx, fftqdfzyx, fftqdfxyz_work,  &
                        fftqdfyzx_work, ftqdfyz, ftqdfzy, ftqdfzy_work,  &
                        work1, work2, nlocalx, nlocaly, nix, niy, niz,   &
                        nx, ny, nz, nprocx, nprocy, id, nthread, nproc,  &
                        comm_city)

#ifdef FFTW
    include 'fftw3.f'
#endif

    ! formal arguments
    real(wp),         intent(inout) :: qdf_work(nlocalx,nlocaly,niz,nprocx,*)
    complex(wp),      intent(inout) :: fftqdfyzx(ny,niz,*)
    complex(wp),      intent(inout) :: fftqdfzyx(nz,ny,*)
    complex(wp),      intent(inout) :: fftqdfxyz_work(nx/2+1,ny,*)
    complex(wp),      intent(inout) :: fftqdfyzx_work(ny,niz,nix/2,*)
    complex(wp),      intent(inout) :: ftqdfyz(niz,*), ftqdfzy(nz,*)
    complex(wp),      intent(inout) :: ftqdfzy_work(niz,niy,*)
    complex(wp),      intent(inout) :: work1(*), work2(*)
    integer,          intent(in)    :: nlocalx, nlocaly
    integer,          intent(in)    :: nix, niy, niz
    integer,          intent(in)    :: nx, ny, nz
    integer,          intent(in)    :: nprocx, nprocy
    integer,          intent(in)    :: id, nthread
    integer,          intent(in)    :: nproc
    integer,          intent(in)    :: comm_city

    ! local variables
    integer(kind=8)   :: plan
    integer           :: ix, iy, iz
    integer           :: iprocx, iprocy
    integer           :: kx, ky
    integer           :: g_ix, g_iy
    integer           :: k, iproc, ierror


    ! x direction
    !
    !$omp barrier
    call zfft1d(work1, nx, 0, work2)

    do iz = id+1, niz, nthread
      do iprocy = 1, nprocy
        ky = (iprocy-1)*nlocaly
        do iy = 1, nlocaly, 2
          g_iy = iy + ky
          do iprocx = 1, nprocx
            kx = (iprocx-1)*nlocalx
            do ix = 1, nlocalx
              g_ix = ix + kx
              work1(g_ix) = cmplx(qdf_work(ix,iy,  iz,iprocx,iprocy), &
                                  qdf_work(ix,iy+1,iz,iprocx,iprocy), &
                                  kind=wp)
            end do
          end do
          call zfft1d(work1, nx, -1, work2)
          fftqdfxyz_work(1,g_iy,iz)   = real(work1(1))
          fftqdfxyz_work(1,g_iy+1,iz) = imag(work1(1))
          do ix = 2, nx/2+1
            fftqdfxyz_work(ix,g_iy,iz) = 0.5_wp             &
                               *(work1(ix)+conjg(work1(nx-ix+2)))
            fftqdfxyz_work(ix,g_iy+1,iz) = (0.0_wp,-0.5_wp) &
                               *(work1(ix)-conjg(work1(nx-ix+2)))
          end do
        end do
      end do
    end do

    ! y direction
    !
    !$omp barrier
    call zfft1d(work1, ny, 0, work2)
    do iz = id+1, niz, nthread
      do ix = 1, nx/2
        do iy = 1, ny
          work1(iy) = fftqdfxyz_work(ix,iy,iz)
        end do
        call zfft1d(work1, ny, -1, work2)
        fftqdfyzx(1:ny,iz,ix) = work1(1:ny)
      end do
    end do
    ix = nx/2 + 1
    do iz = id+1, niz, nthread
      do iy = 1, ny
        work1(iy) = fftqdfxyz_work(ix,iy,iz)
      end do
      call zfft1d(work1, ny, -1, work2)
      do iy = 1, ny
        ftqdfyz(iz,iy) = work1(iy)
      end do
    end do

#ifdef HAVE_MPI_GENESIS

    !$omp barrier
    !$omp master
    call mpi_alltoall(fftqdfyzx,      ny*niz*(nix/2), mpi_wp_complex, &
                      fftqdfyzx_work, ny*niz*(nix/2), mpi_wp_complex, &
                      comm_city, ierror)
    call mpi_alltoall(ftqdfyz,      niy*niz, mpi_wp_complex, &
                      ftqdfzy_work, niy*niz, mpi_wp_complex, &
                      comm_city, ierror)  
    !$omp end master
    !$omp barrier

#endif

    call zfft1d(work1, nz, 0, work2)
    do ix = id+1, nix/2, nthread
      do iy = 1, ny
        do iproc = 1, nproc
          k = (iproc-1)*niz
          do iz = 1, niz
            work1(iz+k) = fftqdfyzx_work(iy,iz,ix,iproc)
          end do
        end do
        call zfft1d(work1, nz, -1, work2)
        fftqdfzyx(1:nz,iy,ix) = work1(1:nz)
      end do
    end do
    do iy = id+1, niy, nthread
      do iproc = 1, nproc
        k = (iproc-1)*niz
        do iz = 1, niz
          work1(iz+k) = ftqdfzy_work(iz,iy,iproc)
        end do
      end do
      call zfft1d(work1, nz, -1, work2)
      ftqdfzy(1:nz,iy) = work1(1:nz)
    end do

    return

  end subroutine fft3d_slab

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bfft3d_slab
  !> @brief        calculate inverse 3D fast fourier transform with 
  !                1D decomposition
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine bfft3d_slab(qdf_work, fftqdfyxz, fftqdfzyx, fftqdfzyx_work, &
                         fftqdfyxz_work, ftqdfyz, ftqdfzy, ftqdf_work,   &
                         work1, work2, nlocalx, nlocaly, nix, niy, niz,  &
                         nx, ny, nz, nprocx, nprocy, id, nthread, nproc, &
                         comm_city)

#ifdef FFTW
    include 'fftw3.f'
#endif

    ! formal arguments
    real(wp),         intent(inout) :: qdf_work(nlocalx,nlocaly,niz,nprocx,*)
    complex(wp),      intent(inout) :: fftqdfzyx(nz,ny,*)
    complex(wp),      intent(inout) :: fftqdfyxz(ny,nix/2,*)
    complex(wp),      intent(inout) :: fftqdfzyx_work(nz,ny,*)
    complex(wp),      intent(inout) :: fftqdfyxz_work(ny,nix/2,niz,*)
    complex(wp),      intent(inout) :: ftqdfyz(niy,niz,*), ftqdfzy(nz,*)
    complex(wp),      intent(inout) :: ftqdf_work(*)
    complex(wp),      intent(inout) :: work1(*), work2(*)
    integer,          intent(in)    :: nlocalx, nlocaly
    integer,          intent(in)    :: nix, niy, niz
    integer,          intent(in)    :: nx, ny, nz
    integer,          intent(in)    :: nprocx, nprocy
    integer,          intent(in)    :: id, nthread
    integer,          intent(in)    :: nproc
    integer,          intent(in)    :: comm_city

    ! local variables
    integer(kind=8)   :: plan
    integer           :: ix, iy, iz
    integer           :: iprocx, iprocy
    integer           :: kx, ky
    integer           :: g_ix, g_iy
    integer           :: k, iproc, ierror
    complex(wp)       :: temp


    ! z direction
    !
    !$omp barrier
    call zfft1d(work1, nz, 0, work2)
    do ix = id+1, nix/2, nthread
      do iy = 1, ny
        work1(1:nz) = fftqdfzyx(1:nz,iy,ix)
        call zfft1d(work1, nz, 1, work2)
        fftqdfzyx_work(1:nz,iy,ix) = work1(1:nz)
      end do
    end do 
    do iy = id+1, niy, nthread
      work1(1:nz) = ftqdfzy(1:nz,iy)
      call zfft1d(work1, nz, 1, work2)
      do iz = 1, nz
        ftqdf_work(iy+(iz-1)*niy) = work1(iz)
      end do
    end do

#ifdef HAVE_MPI_GENESIS

    !$omp barrier
    !$omp master
    call mpi_alltoall(ftqdf_work, niy*niz, mpi_wp_complex, &
                      ftqdfyz,    niy*niz, mpi_wp_complex, &
                      comm_city, ierror)
    !$omp end master
    !$omp barrier

#endif

    ! y direction
    !
    call zfft1d(work1, ny, 0, work2)
    do ix = id+1, nix/2, nthread
      do iz = 1, nz
        do iy = 1, ny
          work1(iy) = fftqdfzyx_work(iz,iy,ix)
        end do
        call zfft1d(work1, ny, 1, work2)
        fftqdfyxz(1:ny,ix,iz) = work1(1:ny)
      end do
    end do

    do iz = id+1, niz, nthread
      do iproc = 1, nproc
        k = (iproc-1)*niy
        do iy = 1, niy
          work1(iy+k) = ftqdfyz(iy,iz,iproc)
        end do
      end do
      call zfft1d(work1, ny, 1, work2)
      do iy = 1, ny
        ftqdf_work(iy+(iz-1)*ny) = work1(iy)
      end do
    end do

#ifdef HAVE_MPI_GENESIS

    !$omp barrier
    !$omp master
    call mpi_alltoall(fftqdfyxz,      ny*niz*(nix/2), mpi_wp_complex, &
                      fftqdfyxz_work, ny*niz*(nix/2), mpi_wp_complex, &
                      comm_city, ierror)
    !$omp end master
    !$omp barrier

#endif

    call zfft1d(work1, nx, 0, work2)
    do iz = id+1, niz, nthread
      do iprocy = 1, nprocy
        ky = (iprocy-1)*nlocaly
        do iy = 1, nlocaly, 2
          g_iy = iy + ky
          work1(1) = cmplx(real(fftqdfyxz_work(g_iy,  1,iz,1),wp), &
                           real(fftqdfyxz_work(g_iy+1,1,iz,1),wp), &
                           kind=wp)
          do ix = 2, nix/2
            temp = (0.0_wp,1.0_wp)*fftqdfyxz_work(g_iy+1,ix,iz,1) 
            work1(ix)      = fftqdfyxz_work(g_iy,ix,iz,1)+temp
            work1(nx-ix+2) = conjg(fftqdfyxz_work(g_iy,ix,iz,1)-temp)
          end do
          do iproc = 2, nproc
            kx = (iproc-1)*nix/2
            do ix = 1, nix/2
              temp = (0.0_wp,1.0_wp)*fftqdfyxz_work(g_iy+1,ix,iz,iproc)
              work1(ix+kx) = fftqdfyxz_work(g_iy,ix,iz,iproc)+temp
              work1(nx-ix-kx+2) = conjg(fftqdfyxz_work(g_iy,ix,iz,iproc)-temp)
            end do
          end do
          work1(nx/2+1) = cmplx(real(ftqdf_work(g_iy+(iz-1)*ny),wp),   &
                                real(ftqdf_work(g_iy+1+(iz-1)*ny),wp), &
                                kind=wp)
          call zfft1d(work1, nx, 1, work2)
          do iprocx = 1, nprocx
            kx = (iprocx-1)*nlocalx
            do ix = 1, nlocalx
              g_ix = ix + kx
              qdf_work(ix,iy,iz,iprocx,iprocy) = real(work1(g_ix),wp)
              qdf_work(ix,iy+1,iz,iprocx,iprocy) = imag(work1(g_ix))
            end do
          end do
        end do
      end do
    end do

    return

  end subroutine bfft3d_slab

end module fft3d_slab_mod
