!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fft3d_1dalltoall_mod
!> @brief   routines for fast fourier transform
!! @authors Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fft3d_1dalltoall_mod

  use constants_mod
  use timers_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif
  use, intrinsic :: iso_c_binding

  implicit none

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    fft3d_1d_alltoall
  !> @brief        calculate 3D fast fourier transform 
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine fft3d_1d_alltoall(   &
                         qdf, qdfyxz, qdfzxy, qdfyxz2, qdfyxz2_work,           &
                         qdfyxz_work, qdfyxz1_work, qdfzxy_work, qdfzxy_work2, &
                         work1, work2, nlocalx, nlocaly, nlocalz, nizx, nizy,  &
                         nlocalx1, x_local1, nx, ny, nz, nprocx, nprocy,       &
                         nprocz, id, nthread, grid_commx, grid_commy,          &
                         grid_commz)

#ifdef FFTW
    include 'fftw3.f03'
#endif

    ! formal arguments
    real(wp),          intent(inout) :: qdf(nlocalx,nlocaly,nizx,*)
    complex(wp),       intent(inout) :: qdfyxz(nlocaly,nlocalx1,*)
    complex(wp),       intent(inout) :: qdfzxy(nlocalz,nlocalx1,*)
    complex(wp),       intent(inout) :: qdfyxz2(nlocaly,nlocalx1,nizy,*)
    complex(wp),       intent(inout) :: qdfyxz2_work(nlocaly,nlocalx1,*)
    complex(wp),       intent(inout) :: qdfyxz_work(nlocaly,nlocalx1,nizy, *)
    complex(wp),       intent(inout) :: qdfyxz1_work(nlocaly,nlocalx1,nizx,*)
    complex(wp),       intent(inout) :: qdfzxy_work(nlocalz,nlocalx1,nlocaly,*)
    complex(wp),       intent(inout) :: qdfzxy_work2(nz,nlocalx1,*)
    complex(wp),       intent(inout) :: work1(*), work2(*)
    integer,           intent(in)    :: nlocalx, nlocaly, nlocalz, nizx, nizy
    integer,           intent(in)    :: nlocalx1
    integer,           intent(in)    :: x_local1
    integer,           intent(in)    :: nx, ny, nz
    integer,           intent(in)    :: nprocx, nprocy, nprocz
    integer,           intent(in)    :: id, nthread
    integer,           intent(in)    :: grid_commx, grid_commy, grid_commz

    ! local variables
    integer(kind=8)    :: plan
    type(C_PTR)        :: p1, p2
    integer            :: ix, iy, iz, ixs, ixe, iys, izs, niz, niy, iyss
    integer            :: nlocalx_half
    integer            :: ix_start1, ix_end1
    integer            :: k, iproc, ierror
    integer            :: nmax
    integer, parameter :: csize = 8


    !$omp barrier

    nmax = max(nx, ny, nz)
#ifdef FFTW
#ifdef _SINGLE
    p1 = fftwf_alloc_complex(int(nmax,csize))
    p2 = fftwf_alloc_complex(int(nmax,csize))
#else
    p1 = fftw_alloc_complex(int(nmax,csize))
    p2 = fftw_alloc_complex(int(nmax,csize))
#endif
    call c_f_pointer(p1, work1, [nmax])
    call c_f_pointer(p2, work2, [nmax])
#endif

    ! x direction
    !
#ifdef FFTW

#ifdef _SINGLE
    call sfftw_plan_dft_1d(plan,nx,work1,work2,FFTW_FORWARD,FFTW_MEASURE)
#else
    call dfftw_plan_dft_1d(plan,nx,work1,work2,FFTW_FORWARD,FFTW_MEASURE)
#endif
#else
    call zfft1d(work1, nx, 0, work2)
#endif

    niz      = nlocalz/nprocx
    ixe = nlocalx/2 + 1
    nlocalx_half = nlocalx/2

    do iz = 1, niz
      do iy = 2*id+1, nlocaly, 2*nthread
        do iproc = 1, nprocx
          k = (iproc-1)*nlocalx
          do ix = 1, nlocalx
            work1(k+ix) = cmplx(qdf(ix,iy,  iz,iproc), &
                                qdf(ix,iy+1,iz,iproc), &
                                kind=wp) 
          end do
        end do
#ifdef FFTW

#ifdef _SINGLE
        call sfftw_execute_dft(plan, work1, work2)
#else
        call dfftw_execute_dft(plan, work1, work2)
#endif
        qdfyxz1_work(iy,1,  iz,1) = real(work2(1),wp)
        qdfyxz1_work(iy+1,1,iz,1) = imag(work2(1))

        do ix = 2, nlocalx_half
          qdfyxz1_work(iy,ix,iz,1)   = 0.5_wp          &
                             *(work2(ix)+conjg(work2(nx-ix+2)))
          qdfyxz1_work(iy+1,ix,iz,1) = (0.0_wp,-0.5_wp) &
                             *(work2(ix)-conjg(work2(nx-ix+2)))
        end do
        qdfyxz1_work(iy,  ixe,iz,1) = real(work2(nx/2+1))
        qdfyxz1_work(iy+1,ixe,iz,1) = imag(work2(nx/2+1))

        do iproc = 2, nprocx

          ix_start1 = (iproc-1)*nlocalx_half + 1
          ix_end1   = ix_start1 + nlocalx_half -1

          do ix = ix_start1, ix_end1
            ixs = ix - ix_start1 + 1
            qdfyxz1_work(iy,ixs,iz,iproc) = 0.5_wp            &
                             *(work2(ix)+conjg(work2(nx-ix+2)))
            qdfyxz1_work(iy+1,ixs,iz,iproc) = (0.0_wp,-0.5_wp) &
                             *(work2(ix)-conjg(work2(nx-ix+2)))
          end do
        end do
#else
        call zfft1d(work1, nx, -1, work2)

        qdfyxz1_work(iy,  1,iz,1) = real(work1(1),wp)
        qdfyxz1_work(iy+1,1,iz,1) = imag(work1(1))

        do ix = 2, nlocalx_half
          qdfyxz1_work(iy,  ix,iz,1)   = 0.5_wp          &
                             *(work1(ix)+conjg(work1(nx-ix+2)))
          qdfyxz1_work(iy+1,ix,iz,1) = (0.0_wp,-0.5_wp) &
                             *(work1(ix)-conjg(work1(nx-ix+2)))
        end do
        qdfyxz1_work(iy,  ixe,iz,1) = real(work1(nx/2+1))
        qdfyxz1_work(iy+1,ixe,iz,1) = imag(work1(nx/2+1))

        do iproc = 2, nprocx 

          ix_start1 = (iproc-1)*nlocalx_half + 1
          ix_end1   = ix_start1 + nlocalx_half -1

          do ix = ix_start1, ix_end1
            ixs = ix - ix_start1 + 1
            qdfyxz1_work(iy,ixs,iz,iproc) = 0.5_wp            &
                             *(work1(ix)+conjg(work1(nx-ix+2)))
            qdfyxz1_work(iy+1,ixs,iz,iproc) = (0.0_wp,-0.5_wp) &
                             *(work1(ix)-conjg(work1(nx-ix+2)))
          end do
        end do
#endif
      end do
    end do

#ifdef FFTW
#ifdef _SINGLE
    call sfftw_destroy_plan(plan)
#else
    call dfftw_destroy_plan(plan)
#endif
#endif

    !$omp barrier
    !$omp master
#ifdef HAVE_MPI_GENESIS

    call mpi_alltoall(qdfyxz1_work, &
                      nlocalx1*nlocaly*niz, mpi_wp_complex, &
                      qdfyxz, &
                      nlocalx1*nlocaly*niz, mpi_wp_complex, &
                      grid_commx, ierror)

#endif

#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdfyxz, &
                      nlocalx1*nlocaly*nlocalz/nprocy, mpi_wp_complex, &
                      qdfyxz_work, &
                      nlocalx1*nlocaly*nlocalz/nprocy, mpi_wp_complex, &
                      grid_commy, ierror)
#endif

    !$omp end master
    !$omp barrier

    ! y direction
    !

#ifdef FFTW
#ifdef _SINGLE
    call sfftw_plan_dft_1d(plan, ny, work1, work2, FFTW_FORWARD, FFTW_MEASURE)
#else
    call dfftw_plan_dft_1d(plan, ny, work1, work2, FFTW_FORWARD, FFTW_MEASURE)
#endif
#else
    call zfft1d(work1, ny, 0, work2)
#endif

    niz = nlocalz/nprocy

    do izs = id+1, niz*x_local1, nthread

      iz = mod(izs-1,niz) + 1
      ix = (izs-1)/niz + 1

      do iproc = 1, nprocy
        k = (iproc-1)*nlocaly
        do iy = 1, nlocaly
          work1(k+iy) = qdfyxz_work(iy,ix,iz,iproc)
        end do
      end do

#ifdef FFTW
#ifdef _SINGLE
      call sfftw_execute_dft(plan, work1, work2)
#else
      call dfftw_execute_dft(plan, work1, work2)
#endif
      do iproc = 1, nprocy
        k = (iproc-1)*nlocaly
        do iy = 1, nlocaly
          qdfyxz2(iy,ix,iz,iproc) = work2(k+iy)
        end do
      end do
#else
      call zfft1d(work1, ny, -1, work2)

      do iproc = 1, nprocy
        k = (iproc-1)*nlocaly
        do iy = 1, nlocaly
          qdfyxz2(iy,ix,iz,iproc) = work1(k+iy)
        end do
      end do
#endif
    end do

#ifdef FFTW
#ifdef _SINGLE
    call sfftw_destroy_plan(plan)
#else
    call dfftw_destroy_plan(plan)
#endif
#endif

    !$omp barrier
    !$omp master

#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdfyxz2, &
                      nlocaly*nlocalx1*niz, mpi_wp_complex, &
                      qdfyxz2_work,  &
                      nlocaly*nlocalx1*niz, mpi_wp_complex, &
                      grid_commy, ierror)
#endif
    !$omp end master
    !$omp barrier

    do iy = id+1, nlocaly, nthread
      do ix = 1, nlocalx1
        do iz = 1, nlocalz
          qdfzxy(iz,ix,iy) = qdfyxz2_work(iy,ix,iz)
        end do
      end do
    end do

    !$omp barrier
    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdfzxy, &
                      nlocalx1*nlocaly*nlocalz/nprocz,mpi_wp_complex, &
                      qdfzxy_work, &
                      nlocalx1*nlocaly*nlocalz/nprocz,mpi_wp_complex, &
                      grid_commz, ierror)
#endif
    !$omp end master
    !$omp barrier

    ! z direction
    !

#ifdef FFTW
#ifdef _SINGLE
    call sfftw_plan_dft_1d(plan, nz, work1, work2, FFTW_FORWARD, FFTW_MEASURE)
#else
    call dfftw_plan_dft_1d(plan, nz, work1, work2, FFTW_FORWARD, FFTW_MEASURE)
#endif
#else
    call zfft1d(work1, nz, 0, work2)
#endif

    niy = nlocaly/nprocz
    do iys = id+1, niy*x_local1, nthread

      iy = mod(iys-1,niy) + 1
      ix = (iys-1)/niy + 1
      do iproc = 1, nprocz
        k = (iproc-1)*nlocalz
        iyss = (iproc-1)*niy + iy
#ifdef FFTW
        do iz = 1, nlocalz
          work1(k+iz) = qdfzxy_work(iz,ix,iyss,1)
        end do
      end do
#ifdef _SINGLE
      call sfftw_execute_dft(plan, work1, work2)
#else
      call dfftw_execute_dft(plan, work1, work2)
#endif 
      do k = 1, nz
        qdfzxy_work2(k,ix,iy) = work2(k)
      end do
#else
        do iz = 1, nlocalz
          qdfzxy_work2(k+iz,ix,iy) = qdfzxy_work(iz,ix,iyss,1)
        end do
      end do
      call zfft1d(qdfzxy_work2(1,ix,iy), nz, -1, work2)
#endif
    end do

    !$omp barrier

#ifdef FFTW
#ifdef _SINGLE
    call sfftw_destroy_plan(plan)
#else
    call dfftw_destroy_plan(plan)
#endif

#ifdef _SINGLE
  call fftwf_free(p1)
  call fftwf_free(p2)
#else
  call fftw_free(p1)
  call fftw_free(p2)
#endif
#endif

    return

  end subroutine fft3d_1d_alltoall

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    fft3d_1d_alltoall_lj
  !> @brief        calculate 3D fast fourier transform 
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine fft3d_1d_alltoall_lj( &
                         qdf, qdfyxz, qdfzxy, qdfyxz2, qdfyxz2_work,           &
                         qdfyxz_work, qdfyxz1_work, qdfzxy_work, qdfzxy_work2, &
                         ldf, ldfyxz, ldfzxy, ldfyxz2, ldfyxz2_work,           &
                         ldfyxz_work, ldfyxz1_work, ldfzxy_work, ldfzxy_work2, &
                         work1, work2, work3, work4,                           &
                         nlocalx, nlocaly, nlocalz, nizx, nizy,                &
                         nlocalx1, x_local1, nx, ny, nz, nprocx, nprocy,       &
                         nprocz, id, nthread, grid_commx, grid_commy,          &
                         grid_commz)

#ifdef FFTW
    include 'fftw3.f03'
#endif

    ! formal arguments
    real(wp),         intent(inout) :: qdf(nlocalx,nlocaly,nizx,*)
    real(wp),         intent(inout) :: ldf(nlocalx,nlocaly,nizx,*)
    complex(wp),      intent(inout) :: qdfyxz(nlocaly,nlocalx1,*)
    complex(wp),      intent(inout) :: qdfzxy(nlocalz,nlocalx1,*)
    complex(wp),      intent(inout) :: qdfyxz2(nlocaly,nlocalx1,nizy,*)
    complex(wp),      intent(inout) :: qdfyxz2_work(nlocaly,nlocalx1,*)
    complex(wp),      intent(inout) :: qdfyxz_work(nlocaly,nlocalx1,nizy, *)
    complex(wp),      intent(inout) :: qdfyxz1_work(nlocaly,nlocalx1,nizx,*)
    complex(wp),      intent(inout) :: qdfzxy_work(nlocalz,nlocalx1,nlocaly,*)
    complex(wp),      intent(inout) :: qdfzxy_work2(nz,nlocalx1,*)
    complex(wp),      intent(inout) :: ldfyxz(nlocaly,nlocalx1,*)
    complex(wp),      intent(inout) :: ldfzxy(nlocalz,nlocalx1,*)
    complex(wp),      intent(inout) :: ldfyxz2(nlocaly,nlocalx1,nizy,*)
    complex(wp),      intent(inout) :: ldfyxz2_work(nlocaly,nlocalx1,*)
    complex(wp),      intent(inout) :: ldfyxz_work(nlocaly,nlocalx1,nizy, *)
    complex(wp),      intent(inout) :: ldfyxz1_work(nlocaly,nlocalx1,nizx,*)
    complex(wp),      intent(inout) :: ldfzxy_work(nlocalz,nlocalx1,nlocaly,*)
    complex(wp),      intent(inout) :: ldfzxy_work2(nz,nlocalx1,*)
    complex(wp),      intent(inout) :: work1(*), work2(*)
    complex(wp),      intent(inout) :: work3(*), work4(*)
    integer,          intent(in)    :: nlocalx, nlocaly, nlocalz, nizx, nizy
    integer,          intent(in)    :: nlocalx1
    integer,          intent(in)    :: x_local1
    integer,          intent(in)    :: nx, ny, nz
    integer,          intent(in)    :: nprocx, nprocy, nprocz
    integer,          intent(in)    :: id, nthread
    integer,          intent(in)    :: grid_commx, grid_commy, grid_commz

    ! local variables
    integer(kind=8)   :: plan1, plan2
    type(C_PTR)       :: p1, p2
    integer           :: ix, iy, iz, ixs, ixe, iys, izs, niz, niy, iyss
    integer           :: nlocalx_half
    integer           :: ix_start1, ix_end1
    integer           :: k, iproc, ierror
    integer           :: nmax
    integer, parameter :: csize = 8


    nmax = max(nx, ny, nz)
#ifdef FFTW
#ifdef _SINGLE
    p1 = fftwf_alloc_complex(int(nmax,csize))
    p2 = fftwf_alloc_complex(int(nmax,csize))
#else
    p1 = fftw_alloc_complex(int(nmax,csize))
    p2 = fftw_alloc_complex(int(nmax,csize))
#endif
    call c_f_pointer(p1, work1, [nmax])
    call c_f_pointer(p2, work2, [nmax])
#endif

    ! x direction
    !
#ifdef FFTW

#ifdef _SINGLE
    call sfftw_plan_dft_1d(plan1, nx, work1, work2, FFTW_FORWARD, FFTW_MEASURE)
    call sfftw_plan_dft_1d(plan2, nx, work3, work4, FFTW_FORWARD, FFTW_MEASURE)
#else
    call dfftw_plan_dft_1d(plan1, nx, work1, work2, FFTW_FORWARD, FFTW_MEASURE)
    call dfftw_plan_dft_1d(plan2, nx, work3, work4, FFTW_FORWARD, FFTW_MEASURE)
#endif
#else
    call zfft1d(work1, nx, 0, work2)
    call zfft1d(work3, nx, 0, work4)
#endif

    niz = nlocalz/nprocx
    ixe = nlocalx/2 + 1
    nlocalx_half = nlocalx/2

    do iz = 1, niz
      do iy = 2*id+1, nlocaly, 2*nthread

        ! Charge part
        !
        do iproc = 1, nprocx
          k = (iproc-1)*nlocalx
          do ix = 1, nlocalx
            work1(k+ix) = cmplx(qdf(ix,iy,  iz,iproc), &
                                qdf(ix,iy+1,iz,iproc), &
                                kind=wp) 
            work3(k+ix) = cmplx(ldf(ix,iy,  iz,iproc), &
                                ldf(ix,iy+1,iz,iproc), &
                                kind=wp) 
          end do
        end do
#ifdef FFTW
#ifdef _SINGLE
        call sfftw_execute_dft(plan1, work1, work2)
        call sfftw_execute_dft(plan2, work3, work4)
#else
        call dfftw_execute_dft(plan1, work1, work2)
        call dfftw_execute_dft(plan2, work3, work4)
#endif
        qdfyxz1_work(iy,1,  iz,1) = real(work2(1),wp)
        qdfyxz1_work(iy+1,1,iz,1) = imag(work2(1))
        ldfyxz1_work(iy,1,  iz,1) = real(work4(1),wp)
        ldfyxz1_work(iy+1,1,iz,1) = imag(work4(1))

        do ix = 2, nlocalx_half
          qdfyxz1_work(iy,ix,iz,1)   = 0.5_wp*(work2(ix)+conjg(work2(nx-ix+2)))
          qdfyxz1_work(iy+1,ix,iz,1) = (0.0_wp,-0.5_wp)*(work2(ix)-conjg(work2(nx-ix+2)))
          ldfyxz1_work(iy,ix,iz,1)   = 0.5_wp*(work4(ix)+conjg(work4(nx-ix+2)))
          ldfyxz1_work(iy+1,ix,iz,1) = (0.0_wp,-0.5_wp)*(work4(ix)-conjg(work4(nx-ix+2)))
        end do
        qdfyxz1_work(iy,  ixe,iz,1) = real(work2(nx/2+1))
        qdfyxz1_work(iy+1,ixe,iz,1) = imag(work2(nx/2+1))
        ldfyxz1_work(iy,  ixe,iz,1) = real(work4(nx/2+1))
        ldfyxz1_work(iy+1,ixe,iz,1) = imag(work4(nx/2+1))

        do iproc = 2, nprocx

          ix_start1 = (iproc-1)*nlocalx_half + 1
          ix_end1   = ix_start1 + nlocalx_half -1

          do ix = ix_start1, ix_end1
            ixs = ix - ix_start1 + 1
            qdfyxz1_work(iy,ixs,iz,iproc) = 0.5_wp*(work2(ix)+conjg(work2(nx-ix+2)))
            qdfyxz1_work(iy+1,ixs,iz,iproc) = (0.0_wp,-0.5_wp)*(work2(ix)-conjg(work2(nx-ix+2)))
            ldfyxz1_work(iy,ixs,iz,iproc) = 0.5_wp*(work2(ix)+conjg(work4(nx-ix+2)))
            ldfyxz1_work(iy+1,ixs,iz,iproc) = (0.0_wp,-0.5_wp)*(work4(ix)-conjg(work4(nx-ix+2)))
          end do
        end do
#else
        call zfft1d(work1, nx, -1, work2)
        call zfft1d(work3, nx, -1, work4)

        qdfyxz1_work(iy,  1,iz,1) = real(work1(1),wp)
        qdfyxz1_work(iy+1,1,iz,1) = imag(work1(1))
        ldfyxz1_work(iy,  1,iz,1) = real(work3(1),wp)
        ldfyxz1_work(iy+1,1,iz,1) = imag(work3(1))

        do ix = 2, nlocalx_half
          qdfyxz1_work(iy,  ix,iz,1)   = 0.5_wp*(work1(ix)+conjg(work1(nx-ix+2)))
          qdfyxz1_work(iy+1,ix,iz,1) = (0.0_wp,-0.5_wp)*(work1(ix)-conjg(work1(nx-ix+2)))
          ldfyxz1_work(iy,  ix,iz,1)   = 0.5_wp*(work3(ix)+conjg(work3(nx-ix+2)))
          ldfyxz1_work(iy+1,ix,iz,1) = (0.0_wp,-0.5_wp)*(work3(ix)-conjg(work3(nx-ix+2)))
        end do
        qdfyxz1_work(iy,  ixe,iz,1) = real(work1(nx/2+1))
        qdfyxz1_work(iy+1,ixe,iz,1) = imag(work1(nx/2+1))
        ldfyxz1_work(iy,  ixe,iz,1) = real(work3(nx/2+1))
        ldfyxz1_work(iy+1,ixe,iz,1) = imag(work3(nx/2+1))

        do iproc = 2, nprocx 

          ix_start1 = (iproc-1)*nlocalx_half + 1
          ix_end1   = ix_start1 + nlocalx_half -1

          do ix = ix_start1, ix_end1
            ixs = ix - ix_start1 + 1
            qdfyxz1_work(iy,ixs,iz,iproc) = 0.5_wp*(work1(ix)+conjg(work1(nx-ix+2)))
            qdfyxz1_work(iy+1,ixs,iz,iproc) = (0.0_wp,-0.5_wp)*(work1(ix)-conjg(work1(nx-ix+2)))
            ldfyxz1_work(iy,ixs,iz,iproc) = 0.5_wp*(work3(ix)+conjg(work3(nx-ix+2)))
            ldfyxz1_work(iy+1,ixs,iz,iproc) = (0.0_wp,-0.5_wp)*(work3(ix)-conjg(work3(nx-ix+2)))
          end do
        end do
#endif

      end do
    end do

#ifdef FFTW
#ifdef _SINGLE
    call sfftw_destroy_plan(plan1)
    call sfftw_destroy_plan(plan2)
#else
    call dfftw_destroy_plan(plan1)
    call dfftw_destroy_plan(plan2)
#endif
#endif

    !$omp barrier
    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdfyxz1_work, &
                      nlocalx1*nlocaly*niz, mpi_wp_complex, &
                      qdfyxz, &
                      nlocalx1*nlocaly*niz, mpi_wp_complex, &
                      grid_commx, ierror)
    call mpi_alltoall(ldfyxz1_work, &
                      nlocalx1*nlocaly*niz, mpi_wp_complex, &
                      ldfyxz, &
                      nlocalx1*nlocaly*niz, mpi_wp_complex, &
                      grid_commx, ierror)
    call mpi_alltoall(qdfyxz, &
                      nlocalx1*nlocaly*nlocalz/nprocy, mpi_wp_complex, &
                      qdfyxz_work, &
                      nlocalx1*nlocaly*nlocalz/nprocy, mpi_wp_complex, &
                      grid_commy, ierror)
    call mpi_alltoall(ldfyxz, &
                      nlocalx1*nlocaly*nlocalz/nprocy, mpi_wp_complex, &
                      ldfyxz_work, &
                      nlocalx1*nlocaly*nlocalz/nprocy, mpi_wp_complex, &
                      grid_commy, ierror)
#endif
    !$omp end master
    !$omp barrier

    ! y direction
    !
#ifdef FFTW
#ifdef _SINGLE
    call sfftw_plan_dft_1d(plan1, ny, work1, work2, FFTW_FORWARD, FFTW_MEASURE)
    call sfftw_plan_dft_1d(plan2, ny, work3, work4, FFTW_FORWARD, FFTW_MEASURE)
#else
    call dfftw_plan_dft_1d(plan1, ny, work1, work2, FFTW_FORWARD, FFTW_MEASURE)
    call dfftw_plan_dft_1d(plan2, ny, work3, work4, FFTW_FORWARD, FFTW_MEASURE)
#endif
#else
    call zfft1d(work1, ny, 0, work2)
    call zfft1d(work3, ny, 0, work4)
#endif

    niz = nlocalz/nprocy

    do izs = id+1, niz*x_local1, nthread

      iz = mod(izs-1,niz) + 1
      ix = (izs-1)/niz + 1

      do iproc = 1, nprocy
        k = (iproc-1)*nlocaly
        do iy = 1, nlocaly
          work1(k+iy) = qdfyxz_work(iy,ix,iz,iproc)
          work3(k+iy) = ldfyxz_work(iy,ix,iz,iproc)
        end do
      end do

#ifdef FFTW
#ifdef _SINGLE
      call sfftw_execute_dft(plan1, work1, work2)
      call sfftw_execute_dft(plan2, work3, work4)
#else
      call dfftw_execute_dft(plan1, work1, work2)
      call dfftw_execute_dft(plan2, work3, work4)
#endif
      do iproc = 1, nprocy
        k = (iproc-1)*nlocaly
        do iy = 1, nlocaly
          qdfyxz2(iy,ix,iz,iproc) = work2(k+iy)
          ldfyxz2(iy,ix,iz,iproc) = work4(k+iy)
        end do
      end do
#else
      call zfft1d(work1, ny, -1, work2)
      call zfft1d(work3, ny, -1, work4)

      do iproc = 1, nprocy
        k = (iproc-1)*nlocaly
        do iy = 1, nlocaly
          qdfyxz2(iy,ix,iz,iproc) = work1(k+iy)
          ldfyxz2(iy,ix,iz,iproc) = work3(k+iy)
        end do
      end do
#endif

    end do

#ifdef FFTW
#ifdef _SINGLE
    call sfftw_destroy_plan(plan1)
    call sfftw_destroy_plan(plan2)
#else
    call dfftw_destroy_plan(plan1)
    call dfftw_destroy_plan(plan2)
#endif
#endif

    !$omp barrier
    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdfyxz2, &
                      nlocaly*nlocalx1*niz, mpi_wp_complex, &
                      qdfyxz2_work,  &
                      nlocaly*nlocalx1*niz, mpi_wp_complex, &
                      grid_commy, ierror)
    call mpi_alltoall(ldfyxz2, &
                      nlocaly*nlocalx1*niz, mpi_wp_complex, &
                      ldfyxz2_work,  &
                      nlocaly*nlocalx1*niz, mpi_wp_complex, &
                      grid_commy, ierror)
#endif
    !$omp end master
    !$omp barrier

    do iy = id+1, nlocaly, nthread
      do ix = 1, nlocalx1
        do iz = 1, nlocalz
          qdfzxy(iz,ix,iy) = qdfyxz2_work(iy,ix,iz)
          ldfzxy(iz,ix,iy) = ldfyxz2_work(iy,ix,iz)
        end do
      end do
    end do

    !$omp barrier
    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdfzxy, &
                      nlocalx1*nlocaly*nlocalz/nprocz,mpi_wp_complex, &
                      qdfzxy_work, &
                      nlocalx1*nlocaly*nlocalz/nprocz,mpi_wp_complex, &
                      grid_commz, ierror)
    call mpi_alltoall(ldfzxy, &
                      nlocalx1*nlocaly*nlocalz/nprocz,mpi_wp_complex, &
                      ldfzxy_work, &
                      nlocalx1*nlocaly*nlocalz/nprocz,mpi_wp_complex, &
                      grid_commz, ierror)
#endif
    !$omp end master
    !$omp barrier

    ! z direction
    !
#ifdef FFTW
#ifdef _SINGLE
    call sfftw_plan_dft_1d(plan1, nz, work1, work2, FFTW_FORWARD, FFTW_MEASURE)
    call sfftw_plan_dft_1d(plan2, nz, work3, work4, FFTW_FORWARD, FFTW_MEASURE)
#else
    call dfftw_plan_dft_1d(plan1, nz, work1, work2, FFTW_FORWARD, FFTW_MEASURE)
    call dfftw_plan_dft_1d(plan2, nz, work3, work4, FFTW_FORWARD, FFTW_MEASURE)
#endif
#else
    call zfft1d(work1, nz, 0, work2)
    call zfft1d(work3, nz, 0, work4)
#endif

    niy   = nlocaly/nprocz
    do iys = id+1, niy*x_local1, nthread

      iy = mod(iys-1,niy) + 1
      ix = (iys-1)/niy + 1

      do iproc = 1, nprocz
        k = (iproc-1)*nlocalz
        iyss = (iproc-1)*niy + iy
#ifdef FFTW
        do iz = 1, nlocalz
          work1(k+iz) = qdfzxy_work(iz,ix,iyss,1)
          work3(k+iz) = ldfzxy_work(iz,ix,iyss,1)
        end do
      end do
#ifdef _SINGLE
      call sfftw_execute_dft(plan1, work1, work2)
      call sfftw_execute_dft(plan2, work3, work4)
#else
      call dfftw_execute_dft(plan1, work1, work2)
      call dfftw_execute_dft(plan2, work3, work4)
#endif 
      do k = 1, nz
        qdfzxy_work2(k,ix,iy) = work2(k)
        ldfzxy_work2(k,ix,iy) = work4(k)
      end do
#else
        do iz = 1, nlocalz
          qdfzxy_work2(k+iz,ix,iy) = qdfzxy_work(iz,ix,iyss,1)
          ldfzxy_work2(k+iz,ix,iy) = ldfzxy_work(iz,ix,iyss,1)
        end do
      end do
      call zfft1d(qdfzxy_work2(1,ix,iy), nz, -1, work2)
      call zfft1d(ldfzxy_work2(1,ix,iy), nz, -1, work4)
#endif

    end do

    !$omp barrier

#ifdef FFTW
#ifdef _SINGLE
    call sfftw_destroy_plan(plan1)
    call sfftw_destroy_plan(plan2)
#else
    call dfftw_destroy_plan(plan1)
    call dfftw_destroy_plan(plan2)
#endif
#ifdef _SINGLE
  call fftwf_free(p1)
  call fftwf_free(p2)
#else
  call fftw_free(p1)
  call fftw_free(p2)
#endif
#endif

    return

  end subroutine fft3d_1d_alltoall_lj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bfft3d_1d_alltoall
  !> @brief        calculate inverse 3D fast fourier transform 
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine bfft3d_1d_alltoall( &
                          qdf, qdf_real, qdfzxy, qdfyxz, qdfzxy_work,          &
                          qdfyxz_work, qdfyxz1_work, qdfzyx_work, work1,       &
                          work2, nlocalx, nlocaly, nlocalz, nlocalx1,          &
                          x_local1, nx, ny, nz, nprocx, nprocy, nprocz,        &
                          id, nthread, my_x_rank, grid_commx, grid_commy,      &
                          grid_commz, niy, nizy, nizx)

#ifdef FFTW
    include 'fftw3.f03'
#endif

    ! formal arguments
    real(wp),          intent(inout) :: qdf(nlocalx,nlocaly,nizx,*)
    real(wp),          intent(inout) :: qdf_real(*)
    complex(wp),       intent(inout) :: qdfzxy(nlocalz,nlocalx1,niy,*)
    complex(wp),       intent(inout) :: qdfyxz(nlocaly,nlocalx1,*)
    complex(wp),       intent(inout) :: qdfzxy_work(nlocalz,nlocalx1,*)
    complex(wp),       intent(inout) :: qdfyxz_work(nlocaly,nlocalx1,nizy, *)
    complex(wp),       intent(inout) :: qdfyxz1_work(nlocaly,nlocalx1,nizx,*)
    complex(wp),       intent(inout) :: qdfzyx_work(nz,niy,*)
    complex(wp),       intent(inout) :: work1(*), work2(*)
    integer,           intent(in)    :: nlocalx, nlocaly, nlocalz
    integer,           intent(in)    :: nlocalx1
    integer,           intent(in)    :: x_local1
    integer,           intent(in)    :: nx, ny, nz
    integer,           intent(in)    :: nprocx, nprocy, nprocz
    integer,           intent(in)    :: id, nthread
    integer,           intent(in)    :: my_x_rank
    integer,           intent(in)    :: grid_commx, grid_commy, grid_commz
    integer,           intent(in)    :: niy, nizy, nizx

    ! local variables
    complex(wp)        :: temp
    integer(kind=8)    :: plan
    type(C_PTR)        :: p1, p2
    integer            :: ix, iy, iz, iys, izs, k1
    integer            :: niz, iz_start, iz_end
    integer            :: k, iproc, ierror
    integer            :: nmax
    integer, parameter :: csize = 8


    nmax = max(nx, ny, nz)
#ifdef FFTW
#ifdef _SINGLE
    p1 = fftwf_alloc_complex(int(nmax,csize))
    p2 = fftwf_alloc_complex(int(nmax,csize))
#else
    p1 = fftw_alloc_complex(int(nmax,csize))
    p2 = fftw_alloc_complex(int(nmax,csize))
#endif
    call c_f_pointer(p1, work1, [nmax])
    call c_f_pointer(p2, work2, [nmax])
#endif

    !$omp barrier

    ! z direction
    !
#ifdef FFTW
#ifdef _SINGLE
    call sfftw_plan_dft_1d(plan, nz, work1, work2, FFTW_BACKWARD, FFTW_MEASURE)
#else
    call dfftw_plan_dft_1d(plan, nz, work1, work2, FFTW_BACKWARD, FFTW_MEASURE)
#endif
#else
    call zfft1d(work1, nz, 0, work2)
#endif

    do iys = id+1, niy*x_local1, nthread 

      iy = mod(iys-1,niy) + 1
      ix = (iys-1)/niy + 1

      work1(1:nz) = qdfzyx_work(1:nz,iy,ix)

#ifdef FFTW
#ifdef _SINGLE
      call sfftw_execute_dft(plan, work1, work2)
#else
      call dfftw_execute_dft(plan, work1, work2)
#endif
      do iproc = 1, nprocz
        k = (iproc-1)*nlocalz
        do iz = 1, nlocalz
          qdfzxy(iz,ix,iy,iproc) = work2(k+iz)
        end do
      end do
#else
      call zfft1d(work1, nz, 1, work2)
      do iproc = 1, nprocz
        k = (iproc-1)*nlocalz
        do iz = 1, nlocalz
          qdfzxy(iz,ix,iy,iproc) = work1(k+iz)
        end do
      end do
#endif
    end do

#ifdef FFTW
#ifdef _SINGLE
    call sfftw_destroy_plan(plan)
#else
    call dfftw_destroy_plan(plan)
#endif
#endif

    !$omp barrier

    !$omp barrier
    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdfzxy, nlocalz*nlocalx1*niy, mpi_wp_complex,      &
                      qdfzxy_work, nlocalz*nlocalx1*niy, mpi_wp_complex, &
                      grid_commz, ierror)
#endif
    !$omp end master
    !$omp barrier

    do iz = id+1, nlocalz, nthread
      do ix = 1, nlocalx1
        do iy = 1, nlocaly
          qdfyxz(iy,ix,iz) = qdfzxy_work(iz,ix,iy)
        end do
      end do
    end do

    niz      = nlocalz/nprocy

    !$omp barrier
    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdfyxz, nlocalx1*nlocaly*niz, mpi_wp_complex,      &
                      qdfyxz_work, nlocalx1*nlocaly*niz, mpi_wp_complex, &
                      grid_commy, ierror)

#endif
    !$omp end master
    !$omp barrier

    ! y direction
    !

#ifdef FFTW
#ifdef _SINGLE
    call sfftw_plan_dft_1d(plan, ny, work1, work2, FFTW_BACKWARD, FFTW_MEASURE)
#else
    call dfftw_plan_dft_1d(plan, ny, work1, work2, FFTW_BACKWARD, FFTW_MEASURE)
#endif
#else
    call zfft1d(work1, ny, 0, work2)
#endif

    do izs = id+1, niz*x_local1, nthread

      iz = mod(izs-1,niz) + 1
      ix = (izs-1)/niz + 1

      do iproc = 1, nprocy
        k = (iproc-1)*nlocaly
        do iy = 1, nlocaly
          work1(k+iy) = qdfyxz_work(iy,ix,iz,iproc)
        end do
      end do

#ifdef FFTW
#ifdef _SINGLE
      call sfftw_execute_dft(plan, work1, work2)
#else
      call dfftw_execute_dft(plan, work1, work2)
#endif
      do iproc = 1, nprocy
        k = (iproc-1)*nlocaly
        do iy = 1, nlocaly
          qdfyxz_work(iy,ix,iz,iproc) = work2(k+iy)
        end do
      end do
#else
      call zfft1d(work1, ny, 1, work2)
      do iproc = 1, nprocy
        k = (iproc-1)*nlocaly
        do iy = 1, nlocaly
          qdfyxz_work(iy,ix,iz,iproc) = work1(k+iy)
        end do
      end do
#endif
    end do

#ifdef FFTW
#ifdef _SINGLE
    call sfftw_destroy_plan(plan)
#else
    call dfftw_destroy_plan(plan)
#endif
#endif

    !$omp barrier
    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdfyxz_work, &
                      nlocaly*nlocalx1*niz, mpi_wp_complex, &
                      qdfyxz, &
                      nlocaly*nlocalx1*niz, mpi_wp_complex, &
                      grid_commy, ierror)

    call mpi_alltoall(qdfyxz, &
                      nlocalx1*nlocaly*nlocalz/nprocx, mpi_wp_complex, &
                      qdfyxz1_work, &
                      nlocalx1*nlocaly*nlocalz/nprocx, mpi_wp_complex, &
                      grid_commx, ierror)
#endif
    !$omp end master
    !$omp barrier

    ! x direction
    !
#ifdef FFTW
#ifdef _SINGLE
    call sfftw_plan_dft_1d(plan, nx, work1, work2, FFTW_BACKWARD, FFTW_MEASURE)
#else
    call dfftw_plan_dft_1d(plan, nx, work1, work2, FFTW_BACKWARD, FFTW_MEASURE)
#endif
#else
    call zfft1d(work1, nx, 0, work2)
#endif

    niz      = nlocalz/nprocx
    iz_start = niz * my_x_rank + 1
    iz_end   = iz_start + niz - 1

    do iz = 1, niz
      do iy = 2*id+1, nlocaly, 2*nthread

        work1(1) = cmplx(real(qdfyxz1_work(iy,1,iz,1),wp),   &
                         real(qdfyxz1_work(iy+1,1,iz,1),wp), &
                         kind=wp)
        do ix = 2, nlocalx1-1
          temp = (0.0_wp,1.0_wp)*qdfyxz1_work(iy+1,ix,iz,1)
          work1(ix) = qdfyxz1_work(iy,ix,iz,1)+temp
          work1(nx-ix+2) = conjg(qdfyxz1_work(iy,ix,iz,1)-temp)
        end do
        temp = (0.0_wp,1.0_wp)*qdfyxz1_work(iy+1,nlocalx1,iz,1)
        work1(nx/2+1) = cmplx(real(qdfyxz1_work(iy,  nlocalx1,iz,1)), &
                              real(qdfyxz1_work(iy+1,nlocalx1,iz,1)), &
                              kind=wp)

        do iproc = 2, nprocx
          k = (iproc-1)*(nlocalx1-1)
          do ix = 1, nlocalx1-1
            k1 = k + ix
            temp = (0.0_wp,1.0_wp)*qdfyxz1_work(iy+1,ix,iz,iproc)
            work1(k1) = qdfyxz1_work(iy,ix,iz,iproc)+temp
            work1(nx-k1+2) = conjg(qdfyxz1_work(iy,ix,iz,iproc)-temp)
          end do
        end do

#ifdef FFTW
#ifdef _SINGLE
        call sfftw_execute_dft(plan, work1, work2)
#else
        call dfftw_execute_dft(plan, work1, work2)
#endif
        do iproc = 1, nprocx
          k = (iproc-1)*nlocalx
          do ix = 1, nlocalx
            qdf(ix,iy,  iz,iproc) = real(work2(k+ix),wp)
            qdf(ix,iy+1,iz,iproc) = imag(work2(k+ix))
          end do
        end do
#else 
        call zfft1d(work1, nx, 1, work2)
        do iproc = 1, nprocx
          k = (iproc-1)*nlocalx
          do ix = 1, nlocalx
            qdf(ix,iy,  iz,iproc) = real(work1(k+ix),wp)
            qdf(ix,iy+1,iz,iproc) = imag(work1(k+ix))
          end do
        end do
#endif
      end do
    end do

#ifdef FFTW
#ifdef _SINGLE
    call sfftw_destroy_plan(plan)
#else
    call dfftw_destroy_plan(plan)
#endif
#endif

    !$omp barrier
    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdf, nlocalx*nlocaly*niz, mpi_wp_real,      &
                      qdf_real, nlocalx*nlocaly*niz, mpi_wp_real, &
                      grid_commx, ierror)

#endif
    !$omp end master
    !$omp barrier

#ifdef FFTW
#ifdef _SINGLE
  call fftwf_free(p1)
  call fftwf_free(p2)
#else
  call fftw_free(p1)
  call fftw_free(p2)
#endif
#endif

    return

  end subroutine bfft3d_1d_alltoall

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bfft3d_1d_alltoall_lj
  !> @brief        calculate inverse 3D fast fourier transform 
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine bfft3d_1d_alltoall_lj( &
                          qdf, qdf_real, qdfzxy, qdfyxz, qdfzxy_work,          &
                          qdfyxz_work, qdfyxz1_work, qdfzyx_work,              &
                          ldf, ldf_real, ldfzxy, ldfyxz, ldfzxy_work,          &
                          ldfyxz_work, ldfyxz1_work, ldfzyx_work,              &
                          work1, work2, work3, work4,                          &
                          nlocalx, nlocaly, nlocalz, nlocalx1, x_local1,       &
                          nx, ny, nz, nprocx, nprocy, nprocz, id, nthread,     &
                          my_x_rank, grid_commx, grid_commy, grid_commz,       &
                          niy, nizy, nizx)

#ifdef FFTW
    include 'fftw3.f03'
#endif

    ! formal arguments
    real(wp),         intent(inout) :: qdf(nlocalx,nlocaly,nizx,*)
    real(wp),         intent(inout) :: qdf_real(*)
    complex(wp),      intent(inout) :: qdfzxy(nlocalz,nlocalx1,niy,*)
    complex(wp),      intent(inout) :: qdfyxz(nlocaly,nlocalx1,*)
    complex(wp),      intent(inout) :: qdfzxy_work(nlocalz,nlocalx1,*)
    complex(wp),      intent(inout) :: qdfyxz_work(nlocaly,nlocalx1,nizy, *)
    complex(wp),      intent(inout) :: qdfyxz1_work(nlocaly,nlocalx1,nizx,*)
    complex(wp),      intent(inout) :: qdfzyx_work(nz,niy,*)
    real(wp),         intent(inout) :: ldf(nlocalx,nlocaly,nizx,*)
    real(wp),         intent(inout) :: ldf_real(*)
    complex(wp),      intent(inout) :: ldfzxy(nlocalz,nlocalx1,niy,*)
    complex(wp),      intent(inout) :: ldfyxz(nlocaly,nlocalx1,*)
    complex(wp),      intent(inout) :: ldfzxy_work(nlocalz,nlocalx1,*)
    complex(wp),      intent(inout) :: ldfyxz_work(nlocaly,nlocalx1,nizy, *)
    complex(wp),      intent(inout) :: ldfyxz1_work(nlocaly,nlocalx1,nizx,*)
    complex(wp),      intent(inout) :: ldfzyx_work(nz,niy,*)
    complex(wp),      intent(inout) :: work1(*), work2(*), work3(*), work4(*)
    integer,          intent(in)    :: nlocalx, nlocaly, nlocalz
    integer,          intent(in)    :: nlocalx1
    integer,          intent(in)    :: x_local1
    integer,          intent(in)    :: nx, ny, nz
    integer,          intent(in)    :: nprocx, nprocy, nprocz
    integer,          intent(in)    :: id, nthread
    integer,          intent(in)    :: my_x_rank
    integer,          intent(in)    :: grid_commx, grid_commy, grid_commz
    integer,          intent(in)    :: niy, nizy, nizx

    ! local variables
    complex(wp)       :: temp
    integer(kind=8)   :: plan1, plan2
    type(C_PTR)       :: p1, p2
    integer           :: ix, iy, iz, iys, izs, k1
    integer           :: niz, iz_start, iz_end
    integer           :: k, iproc, ierror
    integer           :: nmax
    integer, parameter :: csize = 8

    nmax = max(nx, ny, nz)
#ifdef FFTW
#ifdef _SINGLE
    p1 = fftwf_alloc_complex(int(nmax,csize))
    p2 = fftwf_alloc_complex(int(nmax,csize))
#else
    p1 = fftw_alloc_complex(int(nmax,csize))
    p2 = fftw_alloc_complex(int(nmax,csize))
#endif
    call c_f_pointer(p1, work1, [nmax])
    call c_f_pointer(p2, work2, [nmax])
#endif

    ! z direction
    !
#ifdef FFTW
#ifdef _SINGLE
    call sfftw_plan_dft_1d(plan1, nz, work1, work2, FFTW_BACKWARD, FFTW_MEASURE)
    call sfftw_plan_dft_1d(plan2, nz, work3, work4, FFTW_BACKWARD, FFTW_MEASURE)
#else
    call dfftw_plan_dft_1d(plan1, nz, work1, work2, FFTW_BACKWARD, FFTW_MEASURE)
    call dfftw_plan_dft_1d(plan2, nz, work3, work4, FFTW_BACKWARD, FFTW_MEASURE)
#endif
#else
    call zfft1d(work1, nz, 0, work2)
    call zfft1d(work3, nz, 0, work4)
#endif

    do iys = id+1, niy*x_local1, nthread 

      iy = mod(iys-1,niy) + 1
      ix = (iys-1)/niy + 1

      work1(1:nz) = qdfzyx_work(1:nz,iy,ix)
      work3(1:nz) = ldfzyx_work(1:nz,iy,ix)

#ifdef FFTW
#ifdef _SINGLE
      call sfftw_execute_dft(plan1, work1, work2)
      call sfftw_execute_dft(plan2, work3, work4)
#else
      call dfftw_execute_dft(plan1, work1, work2)
      call dfftw_execute_dft(plan2, work3, work4)
#endif
      do iproc = 1, nprocz
        k = (iproc-1)*nlocalz
        do iz = 1, nlocalz
          qdfzxy(iz,ix,iy,iproc) = work2(k+iz)
          ldfzxy(iz,ix,iy,iproc) = work4(k+iz)
        end do
      end do
#else
      call zfft1d(work1, nz, 1, work2)
      call zfft1d(work3, nz, 1, work4)
      do iproc = 1, nprocz
        k = (iproc-1)*nlocalz
        do iz = 1, nlocalz
          qdfzxy(iz,ix,iy,iproc) = work1(k+iz)
          ldfzxy(iz,ix,iy,iproc) = work3(k+iz)
        end do
      end do
#endif

    end do

#ifdef FFTW
#ifdef _SINGLE
    call sfftw_destroy_plan(plan1)
    call sfftw_destroy_plan(plan2)
#else
    call dfftw_destroy_plan(plan1)
    call dfftw_destroy_plan(plan2)
#endif
#endif

    !$omp barrier
    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdfzxy, nlocalz*nlocalx1*niy, mpi_wp_complex,      &
                      qdfzxy_work, nlocalz*nlocalx1*niy, mpi_wp_complex, &
                      grid_commz, ierror)
    call mpi_alltoall(ldfzxy, nlocalz*nlocalx1*niy, mpi_wp_complex,      &
                      ldfzxy_work, nlocalz*nlocalx1*niy, mpi_wp_complex, &
                      grid_commz, ierror)
#endif
    !$omp end master
    !$omp barrier

    do iz = id+1, nlocalz, nthread
      do ix = 1, nlocalx1
        do iy = 1, nlocaly
          qdfyxz(iy,ix,iz) = qdfzxy_work(iz,ix,iy)
          ldfyxz(iy,ix,iz) = ldfzxy_work(iz,ix,iy)
        end do
      end do
    end do

    niz      = nlocalz/nprocy

    !$omp barrier
    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdfyxz, nlocalx1*nlocaly*niz, mpi_wp_complex,      &
                      qdfyxz_work, nlocalx1*nlocaly*niz, mpi_wp_complex, &
                      grid_commy, ierror)
    call mpi_alltoall(ldfyxz, nlocalx1*nlocaly*niz, mpi_wp_complex,      &
                      ldfyxz_work, nlocalx1*nlocaly*niz, mpi_wp_complex, &
                      grid_commy, ierror)
#endif
    !$omp end master
    !$omp barrier

    ! y direction
    !
#ifdef FFTW
#ifdef _SINGLE
    call sfftw_plan_dft_1d(plan1, ny, work1, work2, FFTW_BACKWARD, FFTW_MEASURE)
    call sfftw_plan_dft_1d(plan2, ny, work1, work2, FFTW_BACKWARD, FFTW_MEASURE)
#else
    call dfftw_plan_dft_1d(plan1, ny, work1, work2, FFTW_BACKWARD, FFTW_MEASURE)
    call dfftw_plan_dft_1d(plan2, ny, work1, work2, FFTW_BACKWARD, FFTW_MEASURE)
#endif
#else
    call zfft1d(work1, ny, 0, work2)
    call zfft1d(work3, ny, 0, work4)
#endif

    do izs = id+1, niz*x_local1, nthread

      iz = mod(izs-1,niz) + 1
      ix = (izs-1)/niz + 1

      do iproc = 1, nprocy
        k = (iproc-1)*nlocaly
        do iy = 1, nlocaly
          work1(k+iy) = qdfyxz_work(iy,ix,iz,iproc)
          work3(k+iy) = ldfyxz_work(iy,ix,iz,iproc)
        end do
      end do

#ifdef FFTW
#ifdef _SINGLE
      call sfftw_execute_dft(plan1, work1, work2)
      call sfftw_execute_dft(plan2, work3, work4)
#else
      call dfftw_execute_dft(plan1, work1, work2)
      call dfftw_execute_dft(plan2, work3, work4)
#endif
      do iproc = 1, nprocy
        k = (iproc-1)*nlocaly
        do iy = 1, nlocaly
          qdfyxz_work(iy,ix,iz,iproc) = work2(k+iy)
          ldfyxz_work(iy,ix,iz,iproc) = work4(k+iy)
        end do
      end do
#else
      call zfft1d(work1, ny, 1, work2)
      call zfft1d(work3, ny, 1, work4)
      do iproc = 1, nprocy
        k = (iproc-1)*nlocaly
        do iy = 1, nlocaly
          qdfyxz_work(iy,ix,iz,iproc) = work1(k+iy)
          ldfyxz_work(iy,ix,iz,iproc) = work3(k+iy)
        end do
      end do
#endif
    end do

#ifdef FFTW
#ifdef _SINGLE
    call sfftw_destroy_plan(plan1)
    call sfftw_destroy_plan(plan2)
#else
    call dfftw_destroy_plan(plan1)
    call dfftw_destroy_plan(plan2)
#endif
#endif
    !$omp barrier
    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdfyxz_work, &
                      nlocaly*nlocalx1*niz, mpi_wp_complex, &
                      qdfyxz, &
                      nlocaly*nlocalx1*niz, mpi_wp_complex, &
                      grid_commy, ierror)
    call mpi_alltoall(ldfyxz_work, &
                      nlocaly*nlocalx1*niz, mpi_wp_complex, &
                      ldfyxz, &
                      nlocaly*nlocalx1*niz, mpi_wp_complex, &
                      grid_commy, ierror)
    call mpi_alltoall(qdfyxz, &
                      nlocalx1*nlocaly*nlocalz/nprocx, mpi_wp_complex, &
                      qdfyxz1_work, &
                      nlocalx1*nlocaly*nlocalz/nprocx, mpi_wp_complex, &
                      grid_commx, ierror)
    call mpi_alltoall(ldfyxz, &
                      nlocalx1*nlocaly*nlocalz/nprocx, mpi_wp_complex, &
                      ldfyxz1_work, &
                      nlocalx1*nlocaly*nlocalz/nprocx, mpi_wp_complex, &
                      grid_commx, ierror)
#endif
    !$omp end master
    !$omp barrier

    ! x direction
    !
#ifdef FFTW
#ifdef _SINGLE
    call sfftw_plan_dft_1d(plan1, nx, work1, work2, FFTW_BACKWARD, FFTW_MEASURE)
    call sfftw_plan_dft_1d(plan2, nx, work3, work4, FFTW_BACKWARD, FFTW_MEASURE)
#else
    call dfftw_plan_dft_1d(plan1, nx, work1, work2, FFTW_BACKWARD, FFTW_MEASURE)
    call dfftw_plan_dft_1d(plan2, nx, work3, work4, FFTW_BACKWARD, FFTW_MEASURE)
#endif
#else
    call zfft1d(work1, nx, 0, work2)
    call zfft1d(work3, nx, 0, work4)
#endif

    niz      = nlocalz/nprocx
    iz_start = niz * my_x_rank + 1
    iz_end   = iz_start + niz - 1

    do iz = 1, niz
      do iy = 2*id+1, nlocaly, 2*nthread

        work1(1) = cmplx(real(qdfyxz1_work(iy,1,iz,1),wp),   &
                         real(qdfyxz1_work(iy+1,1,iz,1),wp), &
                         kind=wp)
        work3(1) = cmplx(real(ldfyxz1_work(iy,1,iz,1),wp),   &
                         real(ldfyxz1_work(iy+1,1,iz,1),wp), &
                         kind=wp)
        do ix = 2, nlocalx1-1
          temp = (0.0_wp,1.0_wp)*qdfyxz1_work(iy+1,ix,iz,1)
          work1(ix) = qdfyxz1_work(iy,ix,iz,1)+temp
          work1(nx-ix+2) = conjg(qdfyxz1_work(iy,ix,iz,1)-temp)
          temp = (0.0_wp,1.0_wp)*ldfyxz1_work(iy+1,ix,iz,1)
          work3(ix) = ldfyxz1_work(iy,ix,iz,1)+temp
          work3(nx-ix+2) = conjg(ldfyxz1_work(iy,ix,iz,1)-temp)
        end do
        work1(nx/2+1) = cmplx(real(qdfyxz1_work(iy,  nlocalx1,iz,1)), &
                              real(qdfyxz1_work(iy+1,nlocalx1,iz,1)), &
                              kind=wp)
        work3(nx/2+1) = cmplx(real(ldfyxz1_work(iy,  nlocalx1,iz,1)), &
                              real(ldfyxz1_work(iy+1,nlocalx1,iz,1)), &
                              kind=wp)

        do iproc = 2, nprocx
          k = (iproc-1)*(nlocalx1-1)
          do ix = 1, nlocalx1-1
            k1 = k + ix
            temp = (0.0_wp,1.0_wp)*qdfyxz1_work(iy+1,ix,iz,iproc)
            work1(k1) = qdfyxz1_work(iy,ix,iz,iproc)+temp
            work1(nx-k1+2) = conjg(qdfyxz1_work(iy,ix,iz,iproc)-temp)
            temp = (0.0_wp,1.0_wp)*ldfyxz1_work(iy+1,ix,iz,iproc)
            work3(k1) = ldfyxz1_work(iy,ix,iz,iproc)+temp
            work3(nx-k1+2) = conjg(ldfyxz1_work(iy,ix,iz,iproc)-temp)
          end do
        end do

#ifdef FFTW
#ifdef _SINGLE
        call sfftw_execute_dft(plan1, work1, work2)
        call sfftw_execute_dft(plan2, work3, work4)
#else
        call dfftw_execute_dft(plan1, work1, work2)
        call dfftw_execute_dft(plan2, work3, work4)
#endif
        do iproc = 1, nprocx
          k = (iproc-1)*nlocalx
          do ix = 1, nlocalx
            qdf(ix,iy,  iz,iproc) = real(work2(k+ix),wp)
            qdf(ix,iy+1,iz,iproc) = imag(work2(k+ix))
            ldf(ix,iy,  iz,iproc) = real(work4(k+ix),wp)
            ldf(ix,iy+1,iz,iproc) = imag(work4(k+ix))
          end do
          end do
        end do
#else 
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
#endif
      end do
    end do

#ifdef FFTW
#ifdef _SINGLE
    call sfftw_destroy_plan(plan1)
    call sfftw_destroy_plan(plan2)
#else
    call dfftw_destroy_plan(plan1)
    call dfftw_destroy_plan(plan2)
#endif
#endif
    !$omp barrier
    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdf, nlocalx*nlocaly*niz, mpi_wp_real,      &
                      qdf_real, nlocalx*nlocaly*niz, mpi_wp_real, &
                      grid_commx, ierror)
    call mpi_alltoall(ldf, nlocalx*nlocaly*niz, mpi_wp_real,      &
                      ldf_real, nlocalx*nlocaly*niz, mpi_wp_real, &
                      grid_commx, ierror)
#endif
    !$omp end master
    !$omp barrier

#ifdef FFTW
#ifdef _SINGLE
  call fftwf_free(p1)
  call fftwf_free(p2)
#else
  call fftw_free(p1)
  call fftw_free(p2)
#endif
#endif

    return

  end subroutine bfft3d_1d_alltoall_lj

end module fft3d_1dalltoall_mod
