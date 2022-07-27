!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fft3d_2dalltoall_mod
!> @brief   routines for fast fourier transform with 2d_alltoall
!! @authors Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fft3d_2dalltoall_mod

  use constants_mod
  use timers_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    fft3d_2d_alltoall
  !> @brief        calculate 3D fast fourier transform 
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine fft3d_2d_alltoall(qdf, qdfxyz, qdfyzx_procz, qdfyzx_procy, &
                       qdfxyz_work, qdfyzx, c_work, work1, work2,       &
                       nlocalx, nlocaly, nlocalz, nx, ny, nz,           &
                       nprocx, nprocy, nprocz, nix, nix1, niz, niy,     &
                       id, nthread, my_x_rank, grid_commz, grid_commxy)

#ifdef FFTW
    include 'fftw3.f'
#endif

    ! formal arguments
    real(wp),         intent(inout) :: qdf(nlocalx,nlocaly,niz,*)
    complex(wp),      intent(inout) :: qdfxyz(nix1,nlocaly,niz,*)
    complex(wp),      intent(inout) :: qdfxyz_work(nix1,nlocaly,niz,nprocx,*)
    complex(wp),      intent(inout) :: work1(*), work2(*)
    complex(wp),      intent(inout) :: c_work(nz,niy,*)
    complex(wp),      intent(inout) :: qdfyzx_procy(nlocaly,nlocalz,nix1,*)
    complex(wp),      intent(inout) :: qdfyzx_procz(niy,nlocalz,nix1,*)
    complex(wp),      intent(inout) :: qdfyzx(niy,nlocalz,nix1,*) 

    integer,          intent(in)    :: nlocalx, nlocaly, nlocalz
    integer,          intent(in)    :: nx, ny, nz
    integer,          intent(in)    :: nprocx, nprocy, nprocz
    integer,          intent(in)    :: nix, nix1, niz, niy
    integer,          intent(in)    :: id, nthread
    integer,          intent(in)    :: my_x_rank
    integer,          intent(in)    :: grid_commz
    integer,          intent(in)    :: grid_commxy

    ! local variables
    integer(kind=8)   :: plan
    integer           :: ix, iy, iz, izs
    integer           :: izy, izx, iyx
    integer           :: iz_start, iz_end
    integer           :: k, iproc, ierror
    integer           :: iprocx, iprocy


    !$omp barrier

    ! x direction
    !
    call zfft1d(work1, nx, 0, work2)

    !$omp do collapse(2)
    do iz = 1, niz
      do iy = 1, nlocaly, 2
       
        do iproc = 1, nprocx
          k = (iproc-1)*nlocalx
          do ix = 1, nlocalx
            work1(k+ix) = cmplx(qdf(ix,iy,  iz,iproc), &
                                qdf(ix,iy+1,iz,iproc), &
                                kind=wp) 
          end do
        end do
        call zfft1d(work1, nx, -1, work2)

        qdfxyz(1, iy,  iz,1) = real(work1(1),wp)
        qdfxyz(1, iy+1,iz,1) = imag(work1(1))
        do ix = 2, nix
          qdfxyz(ix,iy,  iz,1) = &
            0.5_wp          *(work1(ix)+conjg(work1(nx-ix+2)))
          qdfxyz(ix,iy+1,iz,1) = &
            (0.0_wp,-0.5_wp)*(work1(ix)-conjg(work1(nx-ix+2)))
        end do
        qdfxyz(nix1,iy,  iz,1) =  real(work1(nx/2+1))
        qdfxyz(nix1,iy+1,iz,1) =  imag(work1(nx/2+1))
        do iproc = 2, nprocx*nprocy
          k = (iproc-1)*nix
          do ix = 1, nix
            qdfxyz(ix,iy,  iz,iproc) = &
              0.5_wp          *(work1(ix+k)+conjg(work1(nx-ix-k+2)))
            qdfxyz(ix,iy+1,iz,iproc) = &
              (0.0_wp,-0.5_wp)*(work1(ix+k)-conjg(work1(nx-ix-k+2)))
          end do
        end do
      end do
    end do

    !$omp barrier

    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdfxyz,      nlocaly*nix1*niz, mpi_wp_complex, &
                      qdfxyz_work, nlocaly*nix1*niz, mpi_wp_complex, &
                      grid_commxy, ierror)
#endif

    !$omp end master
    !$omp barrier

    !$omp do collapse(2)
    do iprocy = 1, nprocy
      do iprocx = 1, nprocx
        iz_start = (iprocx-1)*niz + 1
        iz_end   = iz_start + niz - 1
        do ix = 1, nix1
          do iz = iz_start, iz_end
            izs = iz - iz_start + 1
            do iy = 1, nlocaly
              qdfyzx_procy(iy,iz,ix,iprocy) = qdfxyz_work(ix,iy,izs,iprocx,iprocy)
            end do
          end do
        end do
      end do
    end do

    !$omp barrier

    ! y direction
    !
    call zfft1d(work1, ny, 0, work2)

    !$omp do collapse(2)
    do ix = 1, nix1
      do iz = 1, nlocalz
        do iproc = 1, nprocy
          k = (iproc-1)*nlocaly
          do iy = 1, nlocaly
            work1(k+iy) = qdfyzx_procy(iy,iz,ix,iproc)
          end do
        end do
        call zfft1d(work1, ny, -1, work2)
        do iproc = 1, nprocz
          k = (iproc-1)*niy
          do iy = 1, niy
            qdfyzx(iy,iz,ix,iproc) = work1(k+iy)
          end do
        end do
      end do
    end do

    !$omp barrier

    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdfyzx, nix1*niy*nlocalz, mpi_wp_complex, &
                      qdfyzx_procz, nix1*niy*nlocalz, mpi_wp_complex, &
                      grid_commz, ierror)
#endif
    !$omp end master
    !$omp barrier

    ! z direction
    !
    call zfft1d(work1, nz, 0, work2)

    !$omp do collapse(2)
    do ix = 1, nix1
      do iy = 1, niy
        do iproc = 1, nprocz
          k = (iproc-1)*nlocalz
          do iz = 1, nlocalz
            work1(k+iz) = qdfyzx_procz(iy,iz,ix,iproc)
          end do
        end do
        call zfft1d(work1, nz, -1, work2)
        do iz = 1, nz
          c_work(iz,iy,ix) = work1(iz)
        end do
      end do
    end do

    !$omp barrier

    return

  end subroutine fft3d_2d_alltoall

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    fft3d_2d_alltoall_lj
  !> @brief        calculate 3D fast fourier transform 
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine fft3d_2d_alltoall_lj( &
                         qdf, qdfxyz, qdfyzx_procz, qdfyzx_procy,       &
                         qdfxyz_work, qdfyzx, c_work, work1, work2,     &
                         ldf, ldfxyz, ldfyzx_procz, ldfyzx_procy,       &
                         ldfxyz_work, ldfyzx, c_work_lj, work3, work4,  &
                         nlocalx, nlocaly, nlocalz, nx, ny, nz,         &
                         nprocx, nprocy, nprocz, nix, nix1, niz, niy,   &
                         id, nthread, my_x_rank, grid_commz, grid_commxy)

#ifdef FFTW
    include 'fftw3.f'
#endif

    ! formal arguments
    real(wp),         intent(inout) :: qdf(nlocalx,nlocaly,niz,*)
    complex(wp),      intent(inout) :: qdfxyz(nix1,nlocaly,niz,*)
    complex(wp),      intent(inout) :: qdfxyz_work(nix1,nlocaly,niz,nprocx,*)
    complex(wp),      intent(inout) :: work1(*), work2(*)
    complex(wp),      intent(inout) :: c_work(nz,niy,*)
    complex(wp),      intent(inout) :: qdfyzx_procy(nlocaly,nlocalz,nix1,*)
    complex(wp),      intent(inout) :: qdfyzx_procz(niy,nlocalz,nix1,*)
    complex(wp),      intent(inout) :: qdfyzx(niy,nlocalz,nix1,*) 
    real(wp),         intent(inout) :: ldf(nlocalx,nlocaly,niz,*)
    complex(wp),      intent(inout) :: ldfxyz(nix1,nlocaly,niz,*)
    complex(wp),      intent(inout) :: ldfxyz_work(nix1,nlocaly,niz,nprocx,*)
    complex(wp),      intent(inout) :: work3(*), work4(*)
    complex(wp),      intent(inout) :: c_work_lj(nz,niy,*)
    complex(wp),      intent(inout) :: ldfyzx_procy(nlocaly,nlocalz,nix1,*)
    complex(wp),      intent(inout) :: ldfyzx_procz(niy,nlocalz,nix1,*)
    complex(wp),      intent(inout) :: ldfyzx(niy,nlocalz,nix1,*)

    integer,          intent(in)    :: nlocalx, nlocaly, nlocalz
    integer,          intent(in)    :: nx, ny, nz
    integer,          intent(in)    :: nprocx, nprocy, nprocz
    integer,          intent(in)    :: nix, nix1, niz, niy
    integer,          intent(in)    :: id, nthread
    integer,          intent(in)    :: my_x_rank
    integer,          intent(in)    :: grid_commz
    integer,          intent(in)    :: grid_commxy

    ! local variables
    integer(kind=8)   :: plan1, plan2
    integer           :: ix, iy, iz, izs
    integer           :: izy, izx, iyx
    integer           :: iz_start, iz_end
    integer           :: k, iproc, ierror
    integer           :: iprocx, iprocy


    !$omp barrier

    ! x direction
    !
    call zfft1d(work1, nx, 0, work2)
    call zfft1d(work3, nx, 0, work4)

    !$omp do collapse(2)
    do iz = 1, niz
      do iy = 1, nlocaly, 2
       
        do iproc = 1, nprocx
          k = (iproc-1)*nlocalx
          do ix = 1, nlocalx
            work1(k+ix) = cmplx(qdf(ix,iy,  iz,iproc), &
                                qdf(ix,iy+1,iz,iproc))
            work3(k+ix) = cmplx(ldf(ix,iy,  iz,iproc), &
                                ldf(ix,iy+1,iz,iproc))
          end do
        end do
        call zfft1d(work1, nx, -1, work2)
        call zfft1d(work3, nx, -1, work4)

        qdfxyz(1, iy,  iz,1) = real(work1(1),wp)
        qdfxyz(1, iy+1,iz,1) = imag(work1(1))
        ldfxyz(1, iy,  iz,1) = real(work3(1),wp)
        ldfxyz(1, iy+1,iz,1) = imag(work3(1))
        do ix = 2, nix
          qdfxyz(ix,iy,  iz,1) = &
            0.5_wp          *(work1(ix)+conjg(work1(nx-ix+2)))
          qdfxyz(ix,iy+1,iz,1) = &
            (0.0_wp,-0.5_wp)*(work1(ix)-conjg(work1(nx-ix+2)))
          ldfxyz(ix,iy,  iz,1) = &
            0.5_wp          *(work3(ix)+conjg(work3(nx-ix+2)))
          ldfxyz(ix,iy+1,iz,1) = &
            (0.0_wp,-0.5_wp)*(work3(ix)-conjg(work3(nx-ix+2)))
        end do
        qdfxyz(nix1,iy,  iz,1) =  real(work1(nx/2+1))
        qdfxyz(nix1,iy+1,iz,1) =  imag(work1(nx/2+1))
        ldfxyz(nix1,iy,  iz,1) =  real(work3(nx/2+1))
        ldfxyz(nix1,iy+1,iz,1) =  imag(work3(nx/2+1))
        do iproc = 2, nprocx*nprocy
          k = (iproc-1)*nix
          do ix = 1, nix
            qdfxyz(ix,iy,  iz,iproc) = &
              0.5_wp          *(work1(ix+k)+conjg(work1(nx-ix-k+2)))
            qdfxyz(ix,iy+1,iz,iproc) = &
              (0.0_wp,-0.5_wp)*(work1(ix+k)-conjg(work1(nx-ix-k+2)))
            ldfxyz(ix,iy,  iz,iproc) = &
              0.5_wp          *(work3(ix+k)+conjg(work3(nx-ix-k+2)))
            ldfxyz(ix,iy+1,iz,iproc) = &
              (0.0_wp,-0.5_wp)*(work3(ix+k)-conjg(work3(nx-ix-k+2)))
          end do
        end do
      end do
    end do

    !$omp barrier

    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdfxyz,      nlocaly*nix1*niz, mpi_wp_complex, &
                      qdfxyz_work, nlocaly*nix1*niz, mpi_wp_complex, &
                      grid_commxy, ierror)
    call mpi_alltoall(ldfxyz,      nlocaly*nix1*niz, mpi_wp_complex, &
                      ldfxyz_work, nlocaly*nix1*niz, mpi_wp_complex, &
                      grid_commxy, ierror)
#endif

    !$omp end master
    !$omp barrier

    !$omp do collapse(2)
    do iprocy = 1, nprocy
      do iprocx = 1, nprocx
        iz_start = (iprocx-1)*niz + 1
        iz_end   = iz_start + niz - 1
        do ix = 1, nix1
          do iz = iz_start, iz_end
            izs = iz - iz_start + 1
            do iy = 1, nlocaly
              qdfyzx_procy(iy,iz,ix,iprocy) = qdfxyz_work(ix,iy,izs,iprocx,iprocy)
              ldfyzx_procy(iy,iz,ix,iprocy) = ldfxyz_work(ix,iy,izs,iprocx,iprocy)
            end do
          end do
        end do
      end do
    end do

    !$omp barrier

    ! y direction
    !
    call zfft1d(work1, ny, 0, work2)
    call zfft1d(work3, ny, 0, work4)

    !$omp do collapse(2)
    do ix = 1, nix1
      do iz = 1, nlocalz
        do iproc = 1, nprocy
          k = (iproc-1)*nlocaly
          do iy = 1, nlocaly
            work1(k+iy) = qdfyzx_procy(iy,iz,ix,iproc)
            work3(k+iy) = ldfyzx_procy(iy,iz,ix,iproc)
          end do
        end do
        call zfft1d(work1, ny, -1, work2)
        call zfft1d(work3, ny, -1, work4)
        do iproc = 1, nprocz
          k = (iproc-1)*niy
          do iy = 1, niy
            qdfyzx(iy,iz,ix,iproc) = work1(k+iy)
            ldfyzx(iy,iz,ix,iproc) = work3(k+iy)
          end do
        end do
      end do
    end do

    !$omp barrier

    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdfyzx, nix1*niy*nlocalz, mpi_wp_complex,       &
                      qdfyzx_procz, nix1*niy*nlocalz, mpi_wp_complex, &
                      grid_commz, ierror)
    call mpi_alltoall(ldfyzx, nix1*niy*nlocalz, mpi_wp_complex,       &
                      ldfyzx_procz, nix1*niy*nlocalz, mpi_wp_complex, &
                      grid_commz, ierror)
#endif
    !$omp end master
    !$omp barrier

    ! z direction
    !
    call zfft1d(work1, nz, 0, work2)
    call zfft1d(work3, nz, 0, work4)

    !$omp do collapse(2)
    do ix = 1, nix1
      do iy = 1, niy
        do iproc = 1, nprocz
          k = (iproc-1)*nlocalz
          do iz = 1, nlocalz
            work1(k+iz) = qdfyzx_procz(iy,iz,ix,iproc)
            work3(k+iz) = ldfyzx_procz(iy,iz,ix,iproc)
          end do
        end do
        call zfft1d(work1, nz, -1, work2)
        call zfft1d(work3, nz, -1, work4)
        do iz = 1, nz
          c_work(iz,iy,ix) = work1(iz)
          c_work_lj(iz,iy,ix) = work3(iz)
        end do
      end do
    end do

    !$omp barrier

    return

  end subroutine fft3d_2d_alltoall_lj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bfft3d_2d_alltoall
  !> @brief        calculate inverse 3D fast fourier transform 
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine bfft3d_2d_alltoall(qdf, qdf_real, qdfxyz, qdfzyx_work, &
                        qdfyxz_procy, qdfxyz_work, qdfzyx, c_work,  &
                        work1, work2,                               &
                        nlocalx, nlocaly, nlocalz,                  &
                        nx, ny, nz, nprocx, nprocy, nprocz,         &
                        nix, nix1, niy, niz, id, nthread,           &
                        grid_commx,                                 &
                        grid_commz, grid_commxy)

    ! formal arguments
    real(wp),         intent(inout) :: qdf(nlocalx,nlocaly,niz,*)
    real(wp),         intent(inout) :: qdf_real(*)
    complex(wp),      intent(inout) :: qdfxyz(nix1,nlocaly,niz,*)
    complex(wp),      intent(inout) :: qdfxyz_work(nix1,nlocaly,niz,nprocx,*)
    complex(wp),      intent(inout) :: work1(*), work2(*)
    complex(wp),      intent(inout) :: c_work(nz,niy,*)
    complex(wp),      intent(inout) :: qdfyxz_procy(nlocaly,nix1,nlocalz,*)
    complex(wp),      intent(inout) :: qdfzyx_work(nlocalz,niy,nix1,*)
    complex(wp),      intent(inout) :: qdfzyx(nlocalz,niy,nix1,*)
    integer,          intent(in)    :: nlocalx, nlocaly, nlocalz
    integer,          intent(in)    :: nx, ny, nz
    integer,          intent(in)    :: nprocx, nprocy, nprocz
    integer,          intent(in)    :: nix, nix1, niy, niz
    integer,          intent(in)    :: id, nthread
    integer,          intent(in)    :: grid_commx, grid_commz
    integer,          intent(in)    :: grid_commxy

    ! local variables
    integer           :: ix, iy, iz, izs
    integer           :: izx, iyx
    integer           :: iz_start, iz_end
    integer           :: k, iproc, ierror
    integer           :: iprocx, iprocy
    complex(wp)       :: temp


    !$omp barrier

    ! zdirection
    !
    call zfft1d(work1, nz, 0, work2)

    !$omp do collapse(2)
    do ix = 1, nix1
      do iy = 1, niy
        do iz = 1, nz
          work1(iz) = c_work(iz,iy,ix)
        end do
        call zfft1d(work1, nz, 1, work2)
        do iproc = 1, nprocz
          k = (iproc-1)*nlocalz
          do iz = 1, nlocalz
            qdfzyx_work(iz,iy,ix,iproc) = work1(k+iz)
          end do
        end do
      end do
    end do

    !$omp barrier

    !$omp master

#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdfzyx_work, nix1*niy*nlocalz, mpi_wp_complex, &
                      qdfzyx,      nix1*niy*nlocalz, mpi_wp_complex, &
                      grid_commz, ierror)
#endif
    !$omp end master
    !$omp barrier

    ! y direction
    !
    call zfft1d(work1, ny, 0, work2)
    !$omp do collapse(2)
    do ix = 1, nix1
      do iz = 1, nlocalz
        do iproc = 1, nprocz
          k = (iproc-1)*niy
          do iy = 1, niy
            work1(k+iy) = qdfzyx(iz,iy,ix,iproc)
          end do
        end do
        call zfft1d(work1, ny, 1, work2)
        do iproc = 1, nprocy
          k = (iproc-1)*nlocaly
          do iy = 1, nlocaly
            qdfyxz_procy(iy,ix,iz,iproc) = work1(k+iy)
          end do
        end do
      end do
    end do

    !$omp barrier

    !$omp do collapse(2)
    do iprocy = 1, nprocy
      do iprocx = 1, nprocx
        iz_start = (iprocx-1)*niz + 1
        iz_end   = iz_start + niz - 1
!       !$omp do
        do iz = iz_start, iz_end
          izs = iz - iz_start + 1
          do iy = 1, nlocaly
            do ix = 1, nix1
              qdfxyz_work(ix,iy,izs,iprocx,iprocy) = qdfyxz_procy(iy,ix,iz,iprocy)
            end do
          end do
        end do
      end do
    end do

    !$omp barrier 

    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdfxyz_work, nlocaly*nix1*niz, mpi_wp_complex, &
                      qdfxyz,      nlocaly*nix1*niz, mpi_wp_complex, &
                      grid_commxy, ierror)
#endif
    !$omp end master
    !$omp barrier

    ! x direction
    !
    call zfft1d(work1, nx, 0, work2)

    !$omp do collapse(2)
    do iz = 1, niz
!     do iy = 2*id+1, nlocaly, 2*nthread
      do iy = 1, nlocaly, 2
        work1(1) = cmplx(real(qdfxyz(1,iy,  iz,1),wp), &
                         real(qdfxyz(1,iy+1,iz,1),wp), &
                         kind=wp)
!ocl norecurrence(work1)
        do ix = 2, nix
          temp = (0.0_wp,1.0_wp)*qdfxyz(ix,iy+1,iz,1)
          work1(ix)      = qdfxyz(ix,iy,iz,1) + temp
          work1(nx-ix+2) = conjg(qdfxyz(ix,iy,iz,1)-temp)
        end do
        do iproc = 2, nprocx*nprocy
          k = (iproc-1)*nix
!ocl norecurrence(work1)
          do ix = 1, nix
            temp = (0.0_wp,1.0_wp)*qdfxyz(ix,iy+1,iz,iproc)
            work1(k+ix) = qdfxyz(ix,iy,iz,iproc) + temp
            work1(nx-ix-k+2) = conjg(qdfxyz(ix,iy,iz,iproc)-temp)
          end do
        end do
        work1(nx/2+1) = cmplx(real(qdfxyz(nix1,iy,  iz,1),wp), &
                              real(qdfxyz(nix1,iy+1,iz,1),wp), &
                              kind=wp)
        call zfft1d(work1, nx, 1, work2)
        do iproc = 1, nprocx
          k = (iproc-1)*nlocalx
          do ix = 1, nlocalx
            qdf(ix,iy,  iz,iproc) = real(work1(k+ix),wp)
            qdf(ix,iy+1,iz,iproc) = imag(work1(k+ix))
          end do
        end do
      end do
    end do

    !$omp barrier

    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdf,      nlocalx*nlocaly*niz, mpi_wp_real, &
                      qdf_real, nlocalx*nlocaly*niz, mpi_wp_real, &
                      grid_commx, ierror)
#endif
    !$omp end master
    !$omp barrier

    return

  end subroutine bfft3d_2d_alltoall

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bfft3d_2d_alltoall_lj
  !> @brief        calculate inverse 3D fast fourier transform 
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine bfft3d_2d_alltoall_lj( &
                         qdf, qdf_real, qdfxyz, qdfzyx_work, qdfyxz_procy, &
                         qdfxyz_work, qdfzyx, c_work, work1, work2,        &
                         ldf, ldf_real, ldfxyz, ldfzyx_work, ldfyxz_procy, &
                         ldfxyz_work, ldfzyx, c_work_lj, work3, work4,     &
                         nlocalx, nlocaly, nlocalz, nx, ny, nz,            &
                         nprocx, nprocy, nprocz, nix, nix1, niy, niz, id,  &
                         nthread, grid_commx, grid_commz, grid_commxy)

    ! formal arguments
    real(wp),         intent(inout) :: qdf(nlocalx,nlocaly,niz,*)
    real(wp),         intent(inout) :: qdf_real(*)
    complex(wp),      intent(inout) :: qdfxyz(nix1,nlocaly,niz,*)
    complex(wp),      intent(inout) :: qdfxyz_work(nix1,nlocaly,niz,nprocx,*)
    complex(wp),      intent(inout) :: work1(*), work2(*)
    complex(wp),      intent(inout) :: c_work(nz,niy,*)
    complex(wp),      intent(inout) :: qdfyxz_procy(nlocaly,nix1,nlocalz,*)
    complex(wp),      intent(inout) :: qdfzyx_work(nlocalz,niy,nix1,*)
    complex(wp),      intent(inout) :: qdfzyx(nlocalz,niy,nix1,*)
    real(wp),         intent(inout) :: ldf(nlocalx,nlocaly,niz,*)
    real(wp),         intent(inout) :: ldf_real(*)
    complex(wp),      intent(inout) :: ldfxyz(nix1,nlocaly,niz,*)
    complex(wp),      intent(inout) :: ldfxyz_work(nix1,nlocaly,niz,nprocx,*)
    complex(wp),      intent(inout) :: work3(*), work4(*)
    complex(wp),      intent(inout) :: c_work_lj(nz,niy,*)
    complex(wp),      intent(inout) :: ldfyxz_procy(nlocaly,nix1,nlocalz,*)
    complex(wp),      intent(inout) :: ldfzyx_work(nlocalz,niy,nix1,*)
    complex(wp),      intent(inout) :: ldfzyx(nlocalz,niy,nix1,*)
    integer,          intent(in)    :: nlocalx, nlocaly, nlocalz
    integer,          intent(in)    :: nx, ny, nz
    integer,          intent(in)    :: nprocx, nprocy, nprocz
    integer,          intent(in)    :: nix, nix1, niy, niz
    integer,          intent(in)    :: id, nthread
    integer,          intent(in)    :: grid_commx, grid_commz
    integer,          intent(in)    :: grid_commxy

    ! local variables
    integer           :: ix, iy, iz, izs
    integer           :: izx, iyx
    integer           :: iz_start, iz_end
    integer           :: k, iproc, ierror
    integer           :: iprocx, iprocy
    complex(wp)       :: temp


    !$omp barrier

    ! zdirection
    !
    call zfft1d(work1, nz, 0, work2)
    call zfft1d(work3, nz, 0, work4)

    !$omp do collapse(2)
    do ix = 1, nix1
      do iy = 1, niy
        do iz = 1, nz
          work1(iz) = c_work(iz,iy,ix)
          work3(iz) = c_work_lj(iz,iy,ix)
        end do
        call zfft1d(work1, nz, 1, work2)
        call zfft1d(work3, nz, 1, work4)
        do iproc = 1, nprocz
          k = (iproc-1)*nlocalz
          do iz = 1, nlocalz
            qdfzyx_work(iz,iy,ix,iproc) = work1(k+iz)
            ldfzyx_work(iz,iy,ix,iproc) = work3(k+iz)
          end do
        end do
      end do
    end do

    !$omp barrier
    !$omp master

#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdfzyx_work, nix1*niy*nlocalz, mpi_wp_complex, &
                      qdfzyx,      nix1*niy*nlocalz, mpi_wp_complex, &
                      grid_commz, ierror)
    call mpi_alltoall(ldfzyx_work, nix1*niy*nlocalz, mpi_wp_complex, &
                      ldfzyx,      nix1*niy*nlocalz, mpi_wp_complex, &
                      grid_commz, ierror)
#endif
    !$omp end master
    !$omp barrier

    ! y direction
    !
    call zfft1d(work1, ny, 0, work2)
    call zfft1d(work3, ny, 0, work4)
    !$omp do collapse(2)
    do ix = 1, nix1
      do iz = 1, nlocalz
        do iproc = 1, nprocz
          k = (iproc-1)*niy
          do iy = 1, niy
            work1(k+iy) = qdfzyx(iz,iy,ix,iproc)
            work3(k+iy) = ldfzyx(iz,iy,ix,iproc)
          end do
        end do
        call zfft1d(work1, ny, 1, work2)
        call zfft1d(work3, ny, 1, work4)
        do iproc = 1, nprocy
          k = (iproc-1)*nlocaly
          do iy = 1, nlocaly
            qdfyxz_procy(iy,ix,iz,iproc) = work1(k+iy)
            ldfyxz_procy(iy,ix,iz,iproc) = work3(k+iy)
          end do
        end do
      end do
    end do

    !$omp barrier

    !$omp do collapse(2)
    do iprocy = 1, nprocy
      do iprocx = 1, nprocx
        iz_start = (iprocx-1)*niz + 1
        iz_end   = iz_start + niz - 1
!       !$omp do
        do iz = iz_start, iz_end
          izs = iz - iz_start + 1
          do iy = 1, nlocaly
            do ix = 1, nix1
              qdfxyz_work(ix,iy,izs,iprocx,iprocy) = qdfyxz_procy(iy,ix,iz,iprocy)
              ldfxyz_work(ix,iy,izs,iprocx,iprocy) = ldfyxz_procy(iy,ix,iz,iprocy)
            end do
          end do
        end do
      end do
    end do

    !$omp barrier 
    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdfxyz_work, nlocaly*nix1*niz, mpi_wp_complex, &
                      qdfxyz,      nlocaly*nix1*niz, mpi_wp_complex, &
                      grid_commxy, ierror)
    call mpi_alltoall(ldfxyz_work, nlocaly*nix1*niz, mpi_wp_complex, &
                      ldfxyz,      nlocaly*nix1*niz, mpi_wp_complex, &
                      grid_commxy, ierror)
#endif
    !$omp end master
    !$omp barrier

    ! x direction
    !
    call zfft1d(work1, nx, 0, work2)
    call zfft1d(work3, nx, 0, work4)

    !$omp do collapse(2)
    do iz = 1, niz
!     do iy = 2*id+1, nlocaly, 2*nthread
      do iy = 1, nlocaly, 2
        work1(1) = cmplx(real(qdfxyz(1,iy,  iz,1),wp), &
                         real(qdfxyz(1,iy+1,iz,1),wp))
        work3(1) = cmplx(real(ldfxyz(1,iy,  iz,1),wp), &
                         real(ldfxyz(1,iy+1,iz,1),wp))
!ocl norecurrence(work1,work3)
        do ix = 2, nix
          temp = (0.0_wp,1.0_wp)*qdfxyz(ix,iy+1,iz,1)
          work1(ix)      = qdfxyz(ix,iy,iz,1) + temp
          work1(nx-ix+2) = conjg(qdfxyz(ix,iy,iz,1)-temp)
          temp = (0.0_wp,1.0_wp)*ldfxyz(ix,iy+1,iz,1)
          work3(ix)      = ldfxyz(ix,iy,iz,1) + temp
          work3(nx-ix+2) = conjg(ldfxyz(ix,iy,iz,1)-temp)
        end do
        do iproc = 2, nprocx*nprocy
          k = (iproc-1)*nix
!ocl norecurrence(work1,work3)
          do ix = 1, nix
            temp = (0.0_wp,1.0_wp)*qdfxyz(ix,iy+1,iz,iproc)
            work1(k+ix) = qdfxyz(ix,iy,iz,iproc) + temp
            work1(nx-ix-k+2) = conjg(qdfxyz(ix,iy,iz,iproc)-temp)
            temp = (0.0_wp,1.0_wp)*ldfxyz(ix,iy+1,iz,iproc)
            work3(k+ix) = ldfxyz(ix,iy,iz,iproc) + temp
            work3(nx-ix-k+2) = conjg(ldfxyz(ix,iy,iz,iproc)-temp)
          end do
        end do
        work1(nx/2+1) = cmplx(real(qdfxyz(nix1,iy,  iz,1),wp), &
                              real(qdfxyz(nix1,iy+1,iz,1),wp))
        work3(nx/2+1) = cmplx(real(ldfxyz(nix1,iy,  iz,1),wp), &
                              real(ldfxyz(nix1,iy+1,iz,1),wp))
        call zfft1d(work1, nx, 1, work2)
        call zfft1d(work3, nx, 1, work4)
        do iproc = 1, nprocx
          k = (iproc-1)*nlocalx
          do ix = 1, nlocalx
            qdf(ix,iy,  iz,iproc) = real(work1(k+ix),wp)
            qdf(ix,iy+1,iz,iproc) = imag(work1(k+ix))
            ldf(ix,iy,  iz,iproc) = real(work3(k+ix),wp)
            ldf(ix,iy+1,iz,iproc) = imag(work3(k+ix))
          end do
        end do
      end do
    end do

    !$omp barrier
    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdf,      nlocalx*nlocaly*niz, mpi_wp_real, &
                      qdf_real, nlocalx*nlocaly*niz, mpi_wp_real, &
                      grid_commx, ierror)
    call mpi_alltoall(ldf,      nlocalx*nlocaly*niz, mpi_wp_real, &
                      ldf_real, nlocalx*nlocaly*niz, mpi_wp_real, &
                      grid_commx, ierror)
#endif
    !$omp end master
    !$omp barrier

    return

  end subroutine bfft3d_2d_alltoall_lj

end module fft3d_2dalltoall_mod
