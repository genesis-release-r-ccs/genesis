!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fft3d_mod
!> @brief   routines for fast fourier transform
!! @authors Jaewoon Jung (JJ), Motoshi Kamiya (MK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fft3d_mod

  use constants_mod
  use timers_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    fft3d
  !> @brief        calculate 3D fast fourier transform 
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine fft3d(qdf, qdfxyz, qdfyxz, qdfzyx, qdfxyz_work,        &
                   qdfyxz_work, qdfzyx_work, work1, work2,          &
                   nlocalx, nlocaly, nlocalz, nlocalx1, x_local1,   &
                   nx, ny, nz, nprocx, nprocy, nprocz, id, nthread, &
                   my_x_rank, x_start1, x_end1, y_start, y_end,     &
                   z_start, z_end, grid_commx, grid_commy,          &
                   grid_commz)

#ifdef FFTW
    include 'fftw3.f'
#endif

    ! formal arguments
    real(wp),         intent(inout) :: qdf(nlocalx,nlocaly,nlocalz,*)
    complex(wp),      intent(inout) :: qdfxyz(nlocalx1,nlocaly,*)
    complex(wp),      intent(inout) :: qdfyxz(nlocaly,nlocalx1,*)
    complex(wp),      intent(inout) :: qdfzyx(nlocalz,nlocaly,*)
    complex(wp),      intent(inout) :: qdfxyz_work(nlocalx1,nlocaly,nlocalz,*)
    complex(wp),      intent(inout) :: qdfyxz_work(nlocaly,nlocalx1,nlocalz,*)
    complex(wp),      intent(inout) :: qdfzyx_work(nlocalz,nlocaly,nlocalx1,*)
    complex(wp),      intent(inout) :: work1(*), work2(*)
    integer,          intent(in)    :: nlocalx, nlocaly, nlocalz
    integer,          intent(in)    :: nlocalx1
    integer,          intent(in)    :: x_local1
    integer,          intent(in)    :: nx, ny, nz
    integer,          intent(in)    :: nprocx, nprocy, nprocz
    integer,          intent(in)    :: id, nthread
    integer,          intent(in)    :: my_x_rank
    integer,          intent(in)    :: x_start1, x_end1
    integer,          intent(in)    :: y_start, y_end
    integer,          intent(in)    :: z_start, z_end
    integer,          intent(in)    :: grid_commx, grid_commy, grid_commz

    ! local variables
    integer           :: ix, iy, iz, ixs, iys
    integer           :: i, k, iproc, ierror


    ! x direction
    !

#ifdef FFTW

#ifdef _SINGLE
    call sfftw_plan_dft_1d(plan, nx, work1, work2, FFTW_FORWARD, FFTW_MEASURE)
#else
    call dfftw_plan_dft_1d(plan, nx, work1, work2, FFTW_FORWARD, FFTW_MEASURE)
#endif

#else
    call zfft1d(work1, nx, 0, work2)
#endif

    do iz = id+1, nlocalz, nthread
      do iy = 1, nlocaly, 2

        do iproc = 1, nprocx
          k = (iproc-1)*nlocalx
          do i = 1, nlocalx, 16
            do ix = i, min0(i+15,nlocalx)
              work1(k+ix) = cmplx(qdf(ix,iy,iz,iproc),  &
                                  qdf(ix,iy+1,iz,iproc),&
                                  kind=wp) 
            end do
          end do
        end do

#ifdef FFTW

#ifdef _SINGLE
        call sfftw_execute_dft(plan, work1, work2)
#else
        call dfftw_execute_dft(plan, work1, work2)
#endif

        if (my_x_rank == 0) then
          qdfyxz(iy,1,iz) = real(work2(1),wp)
          qdfyxz(iy+1,1,iz) = imag(work2(1))
          do ix = 2, x_local1
            qdfyxz(iy,ix,iz) = 0.5_wp                 &
                             *(work2(ix)+conjg(work2(nx-ix+2)))
            qdfyxz(iy+1,ix,iz) = (0.0_wp,-0.5_wp)     &
                             *(work2(ix)-conjg(work2(nx-ix+2)))
          end do
        else
          do ix = x_start1, x_end1
            ixs = ix - x_start1 + 1
            qdfyxz(iy,ixs,iz) = 0.5_wp                &
                             *(work2(ix)+conjg(work2(nx-ix+2)))
            qdfyxz(iy+1,ixs,iz) = (0.0_wp,-0.5_wp)    &
                             *(work2(ix)-conjg(work2(nx-ix+2)))
          end do
        end if
#else
        call zfft1d(work1, nx, -1, work2)

        if (my_x_rank == 0) then
          qdfyxz(iy,1,iz) = real(work1(1),wp)
          qdfyxz(iy+1,1,iz) = imag(work1(1))
          do ix = 2, x_local1
            qdfyxz(iy,ix,iz) = 0.5_wp                 &
                             *(work1(ix)+conjg(work1(nx-ix+2)))
            qdfyxz(iy+1,ix,iz) = (0.0_wp,-0.5_wp)     &
                             *(work1(ix)-conjg(work1(nx-ix+2)))
          end do
        else
          do ix = x_start1, x_end1
            ixs = ix - x_start1 + 1
            qdfyxz(iy,ixs,iz) = 0.5_wp                &
                             *(work1(ix)+conjg(work1(nx-ix+2)))
            qdfyxz(iy+1,ixs,iz) = (0.0_wp,-0.5_wp)    &
                             *(work1(ix)-conjg(work1(nx-ix+2)))
          end do
        end if
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
    call mpi_allgather(qdfyxz, nlocalx1*nlocaly*nlocalz,         &
                       mpi_wp_complex, qdfyxz_work,              &
                       nlocalx1*nlocaly*nlocalz, mpi_wp_complex, &
                       grid_commy, ierror)
#else
    qdfyxz_work(1:nlocaly,1:nlocalx1,1:nlocalz,1)  &
     = qdfyxz(1:nlocaly,1:nlocalx1,1:nlocalz)
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

    do iz = id+1, nlocalz, nthread
      do ix = 1, x_local1

        do iproc = 1, nprocy
          k = (iproc-1)*nlocaly
          do i = 1, nlocaly, 16
            do iy = i, min0(i+15,nlocaly)
              work1(k+iy) = qdfyxz_work(iy,ix,iz,iproc)
            end do
          end do
        end do

#ifdef FFTW

#ifdef _SINGLE
        call sfftw_execute_dft(plan, work1, work2)
#else
        call dfftw_execute_dft(plan, work1, work2)
#endif

        do i = y_start, y_end, 16
          do iy = i, min0(i+15,y_end)
            iys = iy - y_start + 1
            qdfzyx(iz,iys,ix) = work2(iy)
          end do
        end do
#else
        call zfft1d(work1, ny, -1, work2)

        do i = y_start, y_end, 16
          do iy = i, min0(i+15,y_end)
            iys = iy - y_start + 1
            qdfzyx(iz,iys,ix) = work1(iy)
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
    call mpi_allgather(qdfzyx, nlocalx1*nlocaly*nlocalz,         &
                       mpi_wp_complex, qdfzyx_work,              &
                       nlocalx1*nlocaly*nlocalz, mpi_wp_complex, &
                       grid_commz, ierror)
#else
    qdfzyx_work(1:nlocalz,1:nlocaly,1:nlocalx1,1)  &
     = qdfzyx(1:nlocalz,1:nlocaly,1:nlocalx1)
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

    do ix = 1, x_local1
      do iy = id+1, nlocaly, nthread

        do iproc = 1, nprocz
          k = (iproc-1)*nlocalz
          do i = 1, nlocalz, 16
            do iz = 1, min0(i+15,nlocalz)
              work1(k+iz) = qdfzyx_work(iz,iy,ix,iproc)
            end do
          end do
        end do

#ifdef FFTW

#ifdef _SINGLE
        call sfftw_execute_dft(plan, work1, work2)
#else
        call dfftw_execute_dft(plan, work1, work2)
#endif

        do iproc = 1, nprocz
          k = (iproc-1)*nlocalz
          qdfzyx_work(1:nlocalz,iy,ix,iproc) = work2(k+1:k+nlocalz)
        end do
#else
        call zfft1d(work1, nz, -1, work2)

        do iproc = 1, nprocz
          k = (iproc-1)*nlocalz
          qdfzyx_work(1:nlocalz,iy,ix,iproc) = work1(k+1:k+nlocalz)
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

    return

  end subroutine fft3d

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bfft3d
  !> @brief        calculate inverse 3D fast fourier transform 
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine bfft3d(qdf, qdf_real, qdfxyz, qdfyxz, qdfzyx,        &
                    qdfxyz_work, qdfyxz_work, qdfzyx_work, work1, &
                    work2, nlocalx, nlocaly, nlocalz, nlocalx1,   &
                    x_local1, nx, ny, nz, nprocx, nprocy, nprocz, &
                    id, nthread, my_x_rank, x_start1, x_end1,     &
                    y_start, y_end, z_start, z_end, grid_commx,   &
                    grid_commy, grid_commz)

#ifdef FFTW
    include 'fftw3.f'
#endif

    ! formal arguments
    real(wp),         intent(inout) :: qdf(nlocalx,nlocaly,nlocalz,*)
    real(wp),         intent(inout) :: qdf_real(*)
    complex(wp),      intent(inout) :: qdfxyz(nlocalx1,nlocaly,*)
    complex(wp),      intent(inout) :: qdfyxz(nlocaly,nlocalx1,*)
    complex(wp),      intent(inout) :: qdfzyx(nlocalz,nlocaly,*)
    complex(wp),      intent(inout) :: qdfxyz_work(nlocalx1,nlocaly,nlocalz,*)
    complex(wp),      intent(inout) :: qdfyxz_work(nlocaly,nlocalx1,nlocalz,*)
    complex(wp),      intent(inout) :: qdfzyx_work(nlocalz,nlocaly,nlocalx1,*)
    complex(wp),      intent(inout) :: work1(*), work2(*)
    integer,          intent(in)    :: nlocalx, nlocaly, nlocalz
    integer,          intent(in)    :: nlocalx1
    integer,          intent(in)    :: x_local1
    integer,          intent(in)    :: nx, ny, nz
    integer,          intent(in)    :: nprocx, nprocy, nprocz
    integer,          intent(in)    :: id, nthread
    integer,          intent(in)    :: my_x_rank
    integer,          intent(in)    :: x_start1, x_end1
    integer,          intent(in)    :: y_start, y_end
    integer,          intent(in)    :: z_start, z_end
    integer,          intent(in)    :: grid_commx, grid_commy, grid_commz

    ! local variables
    integer           :: ix, iy, iz, iys, izs
    integer           :: i, k, k1, iproc, ierror


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

    do ix = 1, x_local1
      do iy = id+1, nlocaly, nthread

        do iproc = 1, nprocz
          k = (iproc-1)*nlocalz
          do i = 1, nlocalz, 16
            do iz = i, min0(i+15,nlocalz)
              work1(k+iz) = qdfzyx_work(iz,iy,ix,iproc)
            end do
          end do
        end do

#ifdef FFTW

#ifdef _SINGLE
        call sfftw_execute_dft(plan, work1, work2)
#else
        call dfftw_execute_dft(plan, work1, work2)
#endif

        do iz = z_start, z_end
          izs = iz - z_start + 1
          qdfzyx(izs,iy,ix) = work2(iz)
        end do
#else
        call zfft1d(work1, nz, 1, work2)

        do iz = z_start, z_end
          izs = iz - z_start + 1
          qdfzyx(izs,iy,ix) = work1(iz)
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
    call mpi_allgather(qdfzyx, nlocalx1*nlocaly*nlocalz,         &
                       mpi_wp_complex, qdfzyx_work,              &
                       nlocalx1*nlocaly*nlocalz, mpi_wp_complex, &
                       grid_commy, ierror)
#else
    qdfzyx_work(1:nlocalz,1:nlocaly,1:nlocalx1,1)  &
     = qdfzyx(1:nlocalz,1:nlocaly,1:nlocalx1)
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

    do iz = id+1, nlocalz, nthread
      do ix = 1, x_local1

        do iproc = 1, nprocy
          k = (iproc-1)*nlocaly
          do i = 1, nlocaly, 16
            do iy = i, min0(i+15,nlocaly)
              work1(k+iy) = qdfzyx_work(iz,iy,ix,iproc)
            end do
          end do
        end do

#ifdef FFTW

#ifdef _SINGLE
        call sfftw_execute_dft(plan, work1, work2)
#else
        call dfftw_execute_dft(plan, work1, work2)
#endif

        do iy = y_start, y_end
          iys = iy - y_start + 1
          qdfyxz(iys,ix,iz) = work2(iy)
        end do
#else
        call zfft1d(work1, ny, 1, work2)

        do iy = y_start, y_end
          iys = iy - y_start + 1
          qdfyxz(iys,ix,iz) = work1(iy)
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
    call mpi_allgather(qdfyxz, nlocalx1*nlocaly*nlocalz,         &
                       mpi_wp_complex, qdfyxz_work,              &
                       nlocalx1*nlocaly*nlocalz, mpi_wp_complex, &
                       grid_commx, ierror)
#else
    qdfyxz_work(1:nlocaly,1:nlocalx1,1:nlocalz,1)  &
     = qdfyxz(1:nlocaly,1:nlocalx1,1:nlocalz)
#endif
    !$omp end master
    !$omp barrier


    ! x direction
    !

#ifdef FFTW

#ifdef _SINGLE
    call sfftw_plan_dft_c2r_1d(plan, nx, work1, work3, FFTW_MEASURE)
#else
    call dfftw_plan_dft_c2r_1d(plan, nx, work1, work3, FFTW_MEASURE)
#endif

#else
    call zfft1d(work1, nx, 0, work2)
#endif

    do iz = id+1, nlocalz, nthread
      do iy = 1, nlocaly

        do i = 1, nlocalx1, 16
          do ix = i, min0(i+15,nlocalx1)
            work1(ix) = qdfyxz_work(iy,ix,iz,1)
          end do
        end do

        do iproc = 2, nprocx
          k = nlocalx1 + (iproc-2)*(nlocalx1-1)
          do i = 1, nlocalx1-1, 16
            do ix = 1, min0(i+15,nlocalx1-1)
              work1(k+ix) = qdfyxz_work(iy,ix,iz,iproc)
            end do
          end do
        end do

#ifdef FFTW

#ifdef _SINGLE
        call sfftw_execute_dft_c2r(plan, work1, work3)
#else
        call dfftw_execute_dft_c2r(plan, work1, work3)
#endif

        k = (iy-1)*nlocalx + (iz-1)*nlocalx*nlocaly
        k1 = my_x_rank*nlocalx
        qdf_real(k+1:k+nlocalx) = work3(k1+1:k1+nlocalx)
#else 
        call zfft1d(work1, nx, 2, work2)

        k = (iy-1)*nlocalx + (iz-1)*nlocalx*nlocaly
        k1 = my_x_rank*nlocalx
        qdf_real(k+1:k+nlocalx) = work1(k1+1:k1+nlocalx)
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

    return

  end subroutine bfft3d

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    fft3d_2d_alltoall
  !> @brief        calculate 3D fast fourier transform 
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine fft3d_2d_alltoall(qdf, qdfyxz_work, work1, work2,          &
                       nlocalx, nlocaly, nlocalz, nlocalx1, x_local1,   &
                       nx, ny, nz, nprocx, nprocy, nprocz, id, nthread, &
                       my_x_rank, my_y_rank, my_z_rank,                 &
                       x_start1, x_end1, y_start, y_end,                &
                       z_start, z_end, grid_commx, grid_commy,          &
                       grid_commz, grid_commxy, c_work,                 &
                       c_work1, c_work2, c_work3)

#ifdef FFTW
    include 'fftw3.f'
#endif

    ! formal arguments
    real(wp),         intent(inout) :: qdf(nlocalx,nlocaly,nlocalz,*)
    complex(wp),      intent(inout) :: qdfyxz_work(nlocaly,nlocalx1,nlocalz,*)
    complex(wp),      intent(inout) :: work1(*), work2(*)
    complex(wp),      intent(inout) :: c_work(nlocalz,nlocalx,nlocaly)
    complex(wp),      intent(inout) :: c_work1(nlocalx*nlocaly*nlocalz)
    complex(wp),      intent(inout) :: c_work2(nlocalx*nlocaly*nlocalz)
    complex(wp),      intent(inout) :: c_work3(nlocaly,nlocalx,nlocalz)

    integer,          intent(in)    :: nlocalx, nlocaly, nlocalz
    integer,          intent(in)    :: nlocalx1
    integer,          intent(in)    :: x_local1
    integer,          intent(in)    :: nx, ny, nz
    integer,          intent(in)    :: nprocx, nprocy, nprocz
    integer,          intent(in)    :: id, nthread
    integer,          intent(in)    :: my_x_rank, my_y_rank, my_z_rank
    integer,          intent(in)    :: x_start1, x_end1
    integer,          intent(in)    :: y_start, y_end
    integer,          intent(in)    :: z_start, z_end
    integer,          intent(in)    :: grid_commx, grid_commy, grid_commz
    integer,          intent(in)    :: grid_commxy

    ! local variables
    integer           :: ix, iy, iz, ixs, ixe, iys, izs, niz, niy, izss
    integer           :: ix_start1, ix_end1, iz_start, iz_end
    integer           :: k, ii, kk, iproc, ierror
    integer           :: nix, nix1, iorg
    integer           :: iprocx, iprocy
    integer           :: ix_start, ix_end
    complex(wp)       :: work1m(nprocx*nlocalx+nprocy*nlocaly+6,4)


    ! x direction
    !
#ifdef FFTW

#ifdef _SINGLE
    call sfftw_plan_dft_1d(plan, nx, work1, work2, FFTW_FORWARD, FFTW_MEASURE)
#else
    call dfftw_plan_dft_1d(plan, nx, work1, work2, FFTW_FORWARD, FFTW_MEASURE)
#endif

#else
    call zfft1d(work1, nx, 0, work2)
#endif

    niz      = nlocalz/nprocx
    iz_start = niz * my_x_rank + 1
    iz_end   = iz_start + niz - 1

    do iz = iz_start, iz_end

      izs = iz - iz_start + 1

      !$omp do
      do iy = 1, nlocaly, 2

        do iproc = 1, nprocx
          k = (iproc-1)*nlocalx
          izss = (iproc-1)*niz + izs
          do ix = 1, nlocalx
            work1(k+ix) = cmplx(qdf(ix,iy,izss,1),   &
                                qdf(ix,iy+1,izss,1), &
                                kind=wp) 
          end do
        end do
#ifdef FFTW

#ifdef _SINGLE
        call sfftw_execute_dft(plan, work1, work2)
#else
        call dfftw_execute_dft(plan, work1, work2)
#endif
        ixe = nlocalx/2 + 1
        qdfyxz_work(iy,1,izs,1) = real(work2(1),wp)
        qdfyxz_work(iy+1,1,izs,1) = imag(work2(1))

        do ix = 2, ixe
          qdfyxz_work(iy,ix,izs,1)   = 0.5_wp           &
                             *(work2(ix)+conjg(work2(nx-ix+2)))
          qdfyxz_work(iy+1,ix,izs,1) = (0.0_wp,-0.5_wp) &
                             *(work2(ix)-conjg(work2(nx-ix+2)))
        end do

        do iproc = 2, nprocx

          ix_start1 = (iproc-1)*nlocalx/2 + 2
          ix_end1   = ix_start1 + nlocalx/2 -1
          izss = (iproc-1)*niz + izs

          do ix = ix_start1, ix_end1
            ixs = ix - ix_start1 + 1
            qdfyxz_work(iy,ixs,izss,1) = 0.5_wp             &
                             *(work2(ix)+conjg(work2(nx-ix+2)))
            qdfyxz_work(iy+1,ixs,izss,1) = (0.0_wp,-0.5_wp) &
                             *(work2(ix)-conjg(work2(nx-ix+2)))
          end do
        end do

#else
        call zfft1d(work1, nx, -1, work2)

        ixe = nlocalx/2 + 1
        qdfyxz_work(iy,1,izs,1) = real(work1(1),wp)
        qdfyxz_work(iy+1,1,izs,1) = imag(work1(1))

        do ix = 2, ixe
          qdfyxz_work(iy,ix,izs,1)   = 0.5_wp           &
                             *(work1(ix)+conjg(work1(nx-ix+2)))
          qdfyxz_work(iy+1,ix,izs,1) = (0.0_wp,-0.5_wp) &
                             *(work1(ix)-conjg(work1(nx-ix+2)))
        end do

        do iproc = 2, nprocx 

          ix_start1 = (iproc-1)*nlocalx/2 + 2
          ix_end1   = ix_start1 + nlocalx/2 -1
          izss = (iproc-1)*niz + izs 

          do ix = ix_start1, ix_end1
            ixs = ix - ix_start1 + 1
            qdfyxz_work(iy,ixs,izss,1) = 0.5_wp             &
                             *(work1(ix)+conjg(work1(nx-ix+2)))
            qdfyxz_work(iy+1,ixs,izss,1) = (0.0_wp,-0.5_wp) &
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
    nix  = (nlocalx1-1)/nprocy
    nix1 = nix+1
    niz  = nlocalz/nprocx

    !$omp do
    do iprocx = 0, nprocx-1
      iz_start = iprocx*niz+1
      iz_end = iz_start+niz-1
      do iprocy = 0, nprocy-1
        ix_start = iprocy*nix+1
        ix_end   = ix_start+nix-1
        if(iprocy .eq. nprocy-1) ix_end = ix_start+nix1-1
        iorg = nlocaly*nix1*niz*(iprocx+nprocx*iprocy)
        do iz=iz_start, iz_end
          izs = iz - iz_start
          do ix=ix_start, ix_end
            ixs = ix - ix_start
            k = iorg+nlocaly*ixs+nlocaly*nix1*izs
            do iy = 1, nlocaly
              c_work1(k+iy) = qdfyxz_work(iy,ix,iz,1)
            end do
          end do
        end do
      end do
    end do

    !$omp barrier
    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(c_work1, nlocaly*nix1*niz, mpi_wp_complex, &
                      c_work2, nlocaly*nix1*niz, mpi_wp_complex, &
                      grid_commxy, ierror)
#endif
    !$omp end master
    !$omp barrier

    !$omp do
    do iprocx = 0, nprocx-1
      iz_start = iprocx*niz+1
      iz_end = iz_start+niz-1
      do iprocy = 0, nprocy-1
        ix_start = iprocy*nix1+1
        ix_end   = ix_start+nix1-1
        iorg = nlocaly*nix1*niz*(iprocx+nprocx*iprocy)
        do iz = iz_start, iz_end
          izs = iz - iz_start
          do ix = ix_start, ix_end
            ixs = ix - ix_start
            k = iorg+nlocaly*ixs+nlocaly*nix1*izs
            do iy = 1, nlocaly
              c_work3(iy,ix,iz) = c_work2(k+iy)
            end do
          end do
        end do
      end do
    end do
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

    if (nlocalz .ge. 4) then

      nix = (nlocalx1-1)/nprocy
      nix1 = nix+1

      !$omp do
      do iz = 1, nlocalz-mod(nlocalz,4), 4
        do ix = 1, nix1
          do ii = 1, 4
            do iproc = 1, nprocy
              k = (iproc-1)*nlocaly
              ixs = (iproc-1)*nix1 + ix
              do iy = 1, nlocaly
                work1m(k+iy,ii) = c_work3(iy,ixs,iz+ii-1)
              end do
            end do
            call zfft1d(work1m(1,ii), ny, -1, work2)
          end do
          do iproc = 1, nprocy
            k = (iproc-1)*nlocaly
            ixs = (iproc-1)*nix1 + ix
            do iy = 1, nlocaly
              kk = iz + nlocalz*(ixs-1) + nlocalz*nix1*nprocy*(iy-1)
              do ii = 1, 4
                c_work2(kk+ii-1) = work1m(k+iy,ii)
              end do
            end do
          end do
        end do
      end do

      if(mod(nlocalz,4).ne.0) then
        do iz = nlocalz-mod(nlocalz,4)+1, nlocalz

          !$omp do
          do ix = 1, nix1
            do iproc = 1, nprocy
              k = (iproc-1)*nlocaly
              ixs = (iproc-1)*nix1 + ix
              do iy = 1, nlocaly
                work1(k+iy) = c_work3(iy,ixs,iz)
              end do
            end do
            call zfft1d(work1, ny, -1, work2)
            do iproc = 1, nprocy
              k = (iproc-1)*nlocaly
              ixs = (iproc-1)*nix1 + ix
              do iy = 1, nlocaly
                kk = iz + nlocalz*(ixs-1) + nlocalz*nix1*nprocy*(iy-1)
                c_work2(kk) = work1(k+iy)
              end do
            end do
          end do
        end do
      endif

    else

      nix = (nlocalx1-1)/nprocy
      nix1 = nix+1

      !$omp do
      do iz = 1, nlocalz
        do ix = 1, nix1
          do iproc = 1, nprocy
            k = (iproc-1)*nlocaly
            ixs = (iproc-1)*nix1 + ix
            do iy = 1, nlocaly
              work1(k+iy) = c_work3(iy,ixs,iz)
            end do
          end do
          call zfft1d(work1, ny, -1, work2)
          do iproc = 1, nprocy
            k = (iproc-1)*nlocaly
            ixs = (iproc-1)*nix1 + ix
            do iy = 1, nlocaly
              kk = iz + nlocalz*(ixs-1) + nlocalz*nix1*nprocy*(iy-1)
              c_work2(kk) = work1(k+iy)
            end do
          end do
        end do
      end do

    endif

    !$omp barrier
    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(c_work2, &
                      nix1*nprocy*nlocaly*nlocalz/nprocz, &
                      mpi_wp_complex, &
                      c_work1, &
                      nix1*nprocy*nlocaly*nlocalz/nprocz, &
                      mpi_wp_complex, &
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
    do iy = 1, niy

      !$omp do
      do ix = 1, nix1*nprocy
        do iproc = 1, nprocz
          k = (iproc-1)*nlocalz
          iys = (iproc-1)*niy + iy
          do iz = 1, nlocalz
            kk = iz + nlocalz*(ix-1) + nlocalz*nix1*nprocy*(iys-1)
            work1(k+iz) = c_work1(kk)
          end do
        end do
        call zfft1d(work1, nz, -1, work2)
        do iproc = 1, nprocz
          k = (iproc-1)*nlocalz
          iys = (iproc-1)*niy + iy
          do iz = 1, nlocalz
            c_work(iz,ix,iys) = work1(k+iz)
          end do
        end do
      end do
    end do

    return

  end subroutine fft3d_2d_alltoall

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bfft3d_2d_alltoall
  !> @brief        calculate inverse 3D fast fourier transform 
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine bfft3d_2d_alltoall(qdf, qdf_real, qdfyxz_work, work1, work2, &
                        nlocalx, nlocaly, nlocalz, nlocalx1,              &
                        x_local1, nx, ny, nz, nprocx, nprocy, nprocz,     &
                        id, nthread, my_x_rank, my_y_rank, my_z_rank,     &
                        x_start1, x_end1,                                 &
                        y_start, y_end, z_start, z_end, grid_commx,       &
                        grid_commy, grid_commz, grid_commxy, c_work,      &
                        c_work1,c_work2,c_work3)

#ifdef FFTW
    include 'fftw3.f'
#endif

    ! formal arguments
    real(wp),         intent(inout) :: qdf(nlocalx,nlocaly,nlocalz,*)
    real(wp),         intent(inout) :: qdf_real(*)
    complex(wp),      intent(inout) :: qdfyxz_work(nlocaly,nlocalx1,nlocalz,*)
    complex(wp),      intent(inout) :: c_work(nlocalz,nlocalx,nlocaly)
    complex(wp),      intent(inout) :: c_work1(nlocalx*nlocaly*nlocalz)
    complex(wp),      intent(inout) :: c_work2(nlocalx*nlocaly*nlocalz)
    complex(wp),      intent(inout) :: c_work3(nlocaly,nlocalx,nlocalz)
    complex(wp),      intent(inout) :: work1(*), work2(*)
    integer,          intent(in)    :: nlocalx, nlocaly, nlocalz
    integer,          intent(in)    :: nlocalx1
    integer,          intent(in)    :: x_local1
    integer,          intent(in)    :: nx, ny, nz
    integer,          intent(in)    :: nprocx, nprocy, nprocz
    integer,          intent(in)    :: id, nthread
    integer,          intent(in)    :: my_x_rank, my_y_rank, my_z_rank
    integer,          intent(in)    :: x_start1, x_end1
    integer,          intent(in)    :: y_start, y_end
    integer,          intent(in)    :: z_start, z_end
    integer,          intent(in)    :: grid_commx, grid_commy, grid_commz
    integer,          intent(in)    :: grid_commxy

    ! local variables
    complex(wp)       :: work1m(nprocx*nlocalx+nprocy*nlocaly+6,4)
    integer           :: ix, iy, iz, iys, izs, izss
    integer           :: niy, niz, iz_start, iz_end
    integer           :: i, k, ii, kk, iproc, ierror
    integer           :: nix, nix1, iorg
    integer           :: iprocx, iprocy, ixs
    integer           :: ix_start, ix_end

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

    !$omp barrier
    nix = (nlocalx1-1)/nprocy
    nix1 = nix+1
    niy = nlocaly/nprocz
    do iy = 1, niy

      !$omp do
      do ix = 1, nix1*nprocy
        do iproc = 1, nprocz
          k = (iproc-1)*nlocalz
          iys = (iproc-1)*niy + iy
          do iz = 1, nlocalz
            work1(k+iz) = c_work(iz,ix,iys)
          end do
        end do
        call zfft1d(work1, nz, 1, work2)
        do iproc = 1, nprocz
          k = (iproc-1)*nlocalz
          iys = (iproc-1)*niy + iy
          do iz = 1, nlocalz
            kk = iz + nlocalz*(ix-1) + nlocalz*nix1*nprocy*(iys-1)
            c_work1(kk) = work1(k+iz)
          end do
        end do
      end do
    end do

#ifdef FFTW
    call dfftw_destroy_plan(plan)
#endif

    !$omp barrier
    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(c_work1, &
                      nix1*nprocy*nlocaly*nlocalz/nprocz, &
                      mpi_wp_complex, &
                      c_work2, &
                      nix1*nprocy*nlocaly*nlocalz/nprocz, &
                      mpi_wp_complex, &
                      grid_commz, ierror)
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

    if (nlocalz .ge. 4) then

      nix = (nlocalx1-1)/nprocy
      nix1 = nix+1

      !$omp do
      do iz = 1, nlocalz-mod(nlocalz,4), 4
        do ix = 1, nix1
          do iproc = 1, nprocy
            k = (iproc-1)*nlocaly
            ixs = (iproc-1)*nix1 + ix
            do iy = 1, nlocaly
              kk = iz + nlocalz*(ixs-1) + nlocalz*nix1*nprocy*(iy-1)
              do ii = 1, 4
                work1m(k+iy,ii) = c_work2(kk+ii-1)
              end do
            end do
          end do
          do ii = 1, 4
            call zfft1d(work1m(1,ii), ny, 1, work2)
            do iproc = 1, nprocy
              k = (iproc-1)*nlocaly
              ixs = (iproc-1)*nix1 + ix
              do iy = 1, nlocaly
                c_work3(iy,ixs,iz+ii-1) = work1m(k+iy,ii)
              end do
            end do
          end do 
        end do
      end do

      if (mod(nlocalz,4) .ne. 0) then

        do iz = nlocalz-mod(nlocalz,4)+1, nlocalz
 
          !$omp do
          do ix = 1, nix1
            do iproc = 1, nprocy
              k = (iproc-1)*nlocaly
              ixs = (iproc-1)*nix1 + ix
              do iy = 1, nlocaly
                kk = iz + nlocalz*(ixs-1) + nlocalz*nix1*nprocy*(iy-1)
                work1(k+iy) = c_work2(kk)
              end do
            end do
            call zfft1d(work1, ny, 1, work2)
            do iproc = 1, nprocy 
              k = (iproc-1)*nlocaly
              ixs = (iproc-1)*nix1 + ix
              do iy = 1, nlocaly
                c_work3(iy,ixs,iz) = work1(k+iy)
              end do 
            end do
          end do
        end do
      endif

    else

      nix = (nlocalx1-1)/nprocy
      nix1 = nix+1

      !$omp do
      do iz = 1, nlocalz
        do ix = 1, nix1
          do iproc = 1, nprocy
            k = (iproc-1)*nlocaly
            ixs = (iproc-1)*nix1 + ix
            do iy = 1, nlocaly
              kk = iz + nlocalz*(ixs-1) + nlocalz*nix1*nprocy*(iy-1)
              work1(k+iy) = c_work2(kk)
            end do
          end do
          call zfft1d(work1, ny, 1, work2)
          do iproc = 1, nprocy
            k = (iproc-1)*nlocaly
            ixs = (iproc-1)*nix1 + ix
            do iy = 1, nlocaly
              c_work3(iy,ixs,iz) = work1(k+iy)
            end do
          end do
        end do
      end do

    endif

    nix = (nlocalx1-1)/nprocy
    nix1 = nix+1
    niz = nlocalz/nprocx

    !$omp do
    do iprocx=0, nprocx-1
      iz_start = iprocx*niz+1
      iz_end = iz_start+niz-1
      do iprocy=0, nprocy-1
        ix_start = iprocy*nix1+1
        ix_end   = ix_start+nix1-1
        iorg = nlocaly*nix1*niz*(iprocx+nprocx*iprocy)
        do iz=iz_start, iz_end
          izs = iz - iz_start
          do ix=ix_start, ix_end
            ixs = ix - ix_start
            k = iorg+nlocaly*ixs+nlocaly*nix1*izs
            do iy=1, nlocaly
              c_work1(k+iy) = c_work3(iy,ix,iz)
            end do
          end do
        end do
      end do
    end do

    !$omp barrier 
    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(c_work1, nlocaly*nix1*niz, mpi_wp_complex, &
                      c_work2, nlocaly*nix1*niz, mpi_wp_complex, &
                      grid_commxy, ierror)
#endif

    !$omp end master
    !$omp barrier

    !$omp do
    do iprocx=0, nprocx-1
      iz_start = iprocx*niz+1
      iz_end = iz_start+niz-1
      do iprocy=0, nprocy-1
        ix_start = iprocy*nix+1
        ix_end   = ix_start+nix-1
        if(iprocy.eq.nprocy-1) ix_end = ix_start+nix1-1
        iorg = nlocaly*nix1*niz*(iprocx+nprocx*iprocy)
        do iz=iz_start, iz_end
          izs = iz - iz_start 
          do ix=ix_start, ix_end
            ixs = ix - ix_start
            k = iorg+nlocaly*ixs+nlocaly*nix1*izs
            do iy=1, nlocaly
              qdfyxz_work(iy,ix,iz,1) = c_work2(k+iy)
            end do
          end do
        end do
      end do
    end do


    ! x direction
    !

#ifdef FFTW

#ifdef _SINGLE
    call sfftw_plan_dft_c2r_1d(plan, nx, work1, work3, FFTW_MEASURE)
#else
    call dfftw_plan_dft_c2r_1d(plan, nx, work1, work3, FFTW_MEASURE)
#endif

#else
    call zfft1d(work1, nx, 0, work2)
#endif

    if(nlocaly.ge.4) then

      niz      = nlocalz/nprocx
      iz_start = niz * my_x_rank + 1
      iz_end   = iz_start + niz - 1
      do iz = iz_start, iz_end
        izs = iz - iz_start + 1

        !$omp do
        do iy = 1, nlocaly-mod(nlocaly,4), 4
          do ix = 1, nlocalx1
            do ii = 1, 4
              work1m(ix,ii) = qdfyxz_work(iy+ii-1,ix,izs,1)
            end do
          end do
          do iproc = 2, nprocx
            k = nlocalx1 + (iproc-2)*(nlocalx1-1)
            izss = (iproc-1)*niz + izs
            do ix = 1, nlocalx1
              do ii = 1, 4
                work1m(k+ix,ii) = qdfyxz_work(iy+ii-1,ix,izss,1)
              end do
            end do
          end do
          do ii = 1, 4
            call zfft1d(work1m(1,ii), nx, 2, work2)
            do iproc = 1, nprocx
              k = (iproc-1)*nlocalx
              izss = (iproc-1)*niz + izs
              do ix = 1, nlocalx
                qdf(ix,iy+ii-1,izss,1) = work1m(k+ix,ii)
              end do
            end do
          end do
        end do
      end do

      if(mod(nlocaly,4).ne.0) then

        do iz = iz_start, iz_end
        izs = iz - iz_start + 1

        !$omp do
          do iy = nlocaly-mod(nlocaly,4)+1, nlocaly
            do ix = 1, nlocalx1
              work1(ix) = qdfyxz_work(iy,ix,izs,1)
            end do
            do iproc = 2, nprocx
              k = nlocalx1 + (iproc-2)*(nlocalx1-1)
              izss = (iproc-1)*niz + izs
              do ix = 1, nlocalx1
                work1(k+ix) = qdfyxz_work(iy,ix,izss,1)
              end do
            end do
            call zfft1d(work1, nx, 2, work2)
            do iproc = 1, nprocx
              k = (iproc-1)*nlocalx
              izss = (iproc-1)*niz + izs
              do ix = 1, nlocalx
                qdf(ix,iy,izss,1) = work1(k+ix)
              end do
            end do
          end do
        end do
      endif

    else

      niz      = nlocalz/nprocx
      iz_start = niz * my_x_rank + 1
      iz_end   = iz_start + niz - 1
      do iz = iz_start+id, iz_end, nthread
        izs = iz - iz_start + 1
        do iy = 1, nlocaly
          do i = 1, nlocalx1, 16
            do ix = i, min0(i+15,nlocalx1)
              work1(ix) = qdfyxz_work(iy,ix,izs,1)
            end do
          end do
          do iproc = 2, nprocx
            k = nlocalx1 + (iproc-2)*(nlocalx1-1)
            izss = (iproc-1)*niz + izs
            do i = 1, nlocalx1-1, 16
              do ix = i, min0(i+15,nlocalx1-1)
                work1(k+ix) = qdfyxz_work(iy,ix,izss,1)
              end do
            end do
          end do
#ifdef FFTW
          call dfftw_execute_dft_c2r(plan, work1, work3)
          k = (iy-1)*nlocalx + (iz-1)*nlocalx*nlocaly
          k1 = my_x_rank*nlocalx
          qdf_real(k+1:k+nlocalx) = work3(k1+1:k1+nlocalx)
#else 
          call zfft1d(work1, nx, 2, work2)
          do iproc = 1, nprocx
            k = (iproc-1)*nlocalx
            izss = (iproc-1)*niz + izs
            do ix = 1, nlocalx
              qdf(ix,iy,izss,1) = work1(k+ix)
            end do
          end do
#endif
        end do
      end do

    endif

#ifdef FFTW
      call dfftw_destroy_plan(plan)
#endif
 
    !$omp barrier
    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdf,                             &
                     nlocalx*nlocaly*niz, mpi_wp_real, &
                     qdf_real,                         &
                     nlocalx*nlocaly*niz, mpi_wp_real, &
                     grid_commx, ierror)
#endif
    !$omp end master
    !$omp barrier

    return

  end subroutine bfft3d_2d_alltoall

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    fft3d_1d_alltoall_x_seiral
  !> @brief        3D fast fourier transform by FFTW for atdyn
  !> @authors      JJ, MK
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine fft3d_1d_alltoall_x_serial(plan_fx, plan_fy, plan_fz, &
                                     qdf, qdfyxz, qdfzxy, qdfyxz2, &
                                     qdfyxz_work, qdfzxy_work, qdfzxy_work2, &
                                     work1, work2, &
                                     nlocaly, nlocalz, nlocalx1, &
                                     nx, ny, nz, nprocy, nprocz, id, nthread, &
                                     my_y_rank, my_z_rank, &
                                     y_start, y_end, z_start, z_end, &
                                     grid_commy, grid_commz)

#ifdef FFTW
    include 'fftw3.f'
#endif

    ! formal arguments
    integer*8,        intent(in)    :: plan_fx, plan_fy, plan_fz
    complex(wp),      intent(inout) :: qdf(nx,nlocaly,nlocalz)
    complex(wp),      intent(inout) :: qdfyxz(nlocaly,nlocalx1,*)
    complex(wp),      intent(inout) :: qdfzxy(nlocalz,nlocalx1,*)
    complex(wp),      intent(inout) :: qdfyxz2(nlocaly,nlocalx1,*)
    complex(wp),      intent(inout) :: qdfyxz_work(nlocaly,nlocalx1,nlocalz,*)
    complex(wp),      intent(inout) :: qdfzxy_work(nlocalz,nlocalx1,nlocaly,*)
    complex(wp),      intent(inout) :: qdfzxy_work2(nz,nlocalx1,nlocaly)
    complex(wp),      intent(inout) :: work1(*), work2(*)
    integer,          intent(in)    :: nlocaly, nlocalz
    integer,          intent(in)    :: nlocalx1
    integer,          intent(in)    :: nx, ny, nz
    integer,          intent(in)    :: nprocy, nprocz
    integer,          intent(in)    :: id, nthread
    integer,          intent(in)    :: my_y_rank, my_z_rank
    integer,          intent(in)    :: y_start, y_end
    integer,          intent(in)    :: z_start, z_end
    integer,          intent(in)    :: grid_commy, grid_commz

    ! local variables
    integer           :: ix, iy, iz, iys, izs, izss, iyss
    integer           :: niy, niz
    integer           :: iy_start, iy_end, iz_start, iz_end
    integer           :: k, iproc, ierror


    ! x direction
    !
    !niz      = nlocalz/nprocx
    niz      = nlocalz

    do iz = id+1, nlocalz, nthread
      do iy = 1, nlocaly, 2
        do ix = 1, nx
          work1(ix) = cmplx(real(qdf(ix,iy,iz),wp),   &
                            real(qdf(ix,iy+1,iz),wp), &
                            kind=wp)
        end do

#ifdef FFTW

#ifdef SINGLE
        call sfftw_execute_dft(plan_fx, work1, work2)
#else
        call dfftw_execute_dft(plan_fx, work1, work2)
#endif

#endif

        qdfyxz(iy,1,iz) = real(work2(1),wp)
        qdfyxz(iy+1,1,iz) = imag(work2(1))

        do ix = 2, nlocalx1
          qdfyxz(iy,ix,iz) = 0.5_wp &
                             *(work2(ix)+conjg(work2(nx-ix+2)))
          qdfyxz(iy+1,ix,iz) = (0.0_wp,-0.5_wp) &
                             *(work2(ix)-conjg(work2(nx-ix+2)))
        end do

      end do
    end do

    !$omp barrier
    !$omp master

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
    niz      = nlocalz/nprocy
    iz_start = niz * my_y_rank + 1
    iz_end   = iz_start + niz - 1

    do iz = iz_start+id, iz_end, nthread
      izs = iz - iz_start + 1
      do ix = 1, nlocalx1

        do iproc = 1, nprocy
          k = (iproc-1)*nlocaly
          izss = (iproc-1)*niz + izs
          do iy = 1, nlocaly
            work1(k+iy) = qdfyxz_work(iy,ix,izss,1)
          end do
        end do

#ifdef FFTW

#ifdef SINGLE
       call sfftw_execute_dft(plan_fy, work1, work2)
#else
       call dfftw_execute_dft(plan_fy, work1, work2)
#endif

#endif

       do iproc = 1, nprocy
         k = (iproc-1)*nlocaly
         izss = (iproc-1)*niz + izs
         do iy = 1, nlocaly
           qdfyxz2(iy,ix,izss) = work2(k+iy)
         end do
       end do

      end do
    end do

    !$omp barrier
    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdfyxz2, &
                      nlocaly*nlocalx1*niz, mpi_wp_complex, &
                      qdfyxz, &
                      nlocaly*nlocalx1*niz, mpi_wp_complex, &
                      grid_commy, ierror)
#endif
    !$omp end master
    !$omp barrier

    !$omp do
    do iz = 1, nlocalz
      do ix = 1, nlocalx1
        do iy = 1, nlocaly
          qdfzxy(iz,ix,iy) = qdfyxz(iy,ix,iz)
        end do
      end do
    end do

    !$omp master
#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdfzxy, &
                      nlocalx1*nlocaly*nlocalz/nprocz, mpi_wp_complex, &
                      qdfzxy_work, &
                      nlocalx1*nlocaly*nlocalz/nprocz, mpi_wp_complex, &
                      grid_commz, ierror)
#endif
    !$omp end master
    !$omp barrier


    ! z direction
    !
    niy      = nlocaly/nprocz
    iy_start = niy * my_z_rank + 1
    iy_end   = iy_start + niy - 1

    do ix = 1, nlocalx1
      do iy = iy_start+id, iy_end, nthread
        iys = iy - iy_start + 1

        do iproc = 1, nprocz
          k = (iproc-1)*nlocalz
          iyss = (iproc-1)*niy + iys
          do iz = 1, nlocalz
            qdfzxy_work2(k+iz,ix,iy) = qdfzxy_work(iz,ix,iyss,1)
          end do
        end do

#ifdef FFTW

#ifdef SINGLE
        call sfftw_execute_dft(plan_fz, qdfzxy_work2(1,ix,iy), work2)
#else
        call dfftw_execute_dft(plan_fz, qdfzxy_work2(1,ix,iy), work2)
#endif

#endif

        do k = 1, nz
          qdfzxy_work2(k,ix,iys) = work2(k)
        end do

      end do
    end do

    return

  end subroutine fft3d_1d_alltoall_x_serial

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bfft3d_1d_alltoall_x_serial
  !> @brief        3D inverse fourier transform for atdyn by FFTW
  !> @authors      JJ, MK
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine bfft3d_1d_alltoall_x_serial(plan_bx, plan_by, plan_bz, &
                                      qdf, qdfyxz, qdfzxy, &
                                      qdfyxz_work, qdfzxy_work, qdfzxy_work2, &
                                      work1, work2, &
                                      nlocaly, nlocalz, nlocalx1, &
                                      nx, ny, nz, nprocy, nprocz, id, nthread, &
                                      my_y_rank, my_z_rank, &
                                      y_start, y_end, z_start, z_end, &
                                      grid_commy, grid_commz, niy)

#ifdef FFTW
    include "fftw3.f"
#endif

    ! formal arguments
    integer*8,        intent(in)    :: plan_bx, plan_by, plan_bz
    complex(wp),      intent(inout) :: qdf(nx,nlocaly,nlocalz,*)
    complex(wp),      intent(inout) :: qdfyxz(nlocaly,nlocalx1,*)
    complex(wp),      intent(inout) :: qdfzxy(nlocalz,nlocalx1,*)
    complex(wp),      intent(inout) :: qdfyxz_work(nlocaly,nlocalx1,nlocalz,*)
    complex(wp),      intent(inout) :: qdfzxy_work(nlocalz,nlocalx1,nlocaly,*)
    complex(wp),      intent(inout) :: qdfzxy_work2(nz,nlocalx1,nlocaly)
    complex(wp),      intent(inout) :: work1(*), work2(*)
    integer,          intent(in)    :: nlocaly, nlocalz
    integer,          intent(in)    :: nlocalx1
    integer,          intent(in)    :: nx, ny, nz
    integer,          intent(in)    :: nprocy, nprocz
    integer,          intent(in)    :: id, nthread
    integer,          intent(in)    :: my_y_rank, my_z_rank
    integer,          intent(in)    :: y_start, y_end
    integer,          intent(in)    :: z_start, z_end
    integer,          intent(in)    :: grid_commy, grid_commz
    integer,          intent(in)    :: niy

    ! local variables
    real(wp)          :: work3(nx)
    integer           :: ix, iy, iz, iys, iyss, izs, izss
    integer           :: niz, iy_start, iy_end, iz_start, iz_end
    integer           :: i, k, iproc, ierror


    ! z direction
    !

    iy_start = niy * my_z_rank + 1
    iy_end   = iy_start + niy - 1

    do ix = 1, nlocalx1
      do iy = iy_start+id, iy_end, nthread
        iys = iy - iy_start + 1

        do i = 1, nz, 16
          do iz = i, min0(i+15,nz)
            work1(iz) = qdfzxy_work2(iz,ix,iys)
          end do
        end do

#ifdef FFTW

#ifdef _SINGLE
        call sfftw_execute_dft(plan_bz, work1, work2)
#else
        call dfftw_execute_dft(plan_bz, work1, work2)
#endif

#endif

        do iproc = 1, nprocz
          k = (iproc-1)*nlocalz
          iyss = (iproc-1)*niy + iys
          do iz = 1, nlocalz
            qdfzxy_work(iz,ix,iyss,1) = work2(k+iz)
          end do
        end do

      end do
    end do

    !$omp barrier
    !$omp master

#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdfzxy_work, &
                      nlocalz*nlocalx1*niy, mpi_wp_complex, &
                      qdfzxy, &
                      nlocalz*nlocalx1*niy, mpi_wp_complex, &
                      grid_commz, ierror)
#endif

    !$omp end master
    !$omp barrier

    !$omp do
    do iz = 1, nlocalz
      do ix = 1, nlocalx1
        do iy = 1, nlocaly
          qdfyxz(iy,ix,iz) = qdfzxy(iz,ix,iy)
        end do
      end do
    end do

    !$omp master

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

    niz      = nlocalz / nprocy
    iz_start = niz * my_y_rank + 1
    iz_end   = iz_start + niz - 1

    do iz = iz_start+id, iz_end, nthread
      izs = iz - iz_start + 1
      do ix = 1, nlocalx1

        do iproc = 1, nprocy
          k = (iproc-1)*nlocaly
          izss = (iproc-1)*niz + izs
          do i = 1, nlocaly, 16
            do iy = i, min0(i+15,nlocaly)
              work1(k+iy) = qdfyxz_work(iy,ix,izss,1)
            end do
          end do
        end do

#ifdef FFTW

#ifdef _SINGLE
        call sfftw_execute_dft(plan_by, work1, work2)
#else
        call dfftw_execute_dft(plan_by, work1, work2)
#endif

#endif

        do iproc = 1, nprocy
          k = (iproc-1)*nlocaly
          izss = (iproc-1)*niz + izs
          do iy = 1, nlocaly
            qdfyxz_work(iy,ix,izss,1) = work2(k+iy)
          end do
        end do

      end do
    end do


    !$omp barrier
    !$omp master

#ifdef HAVE_MPI_GENESIS
    call mpi_alltoall(qdfyxz_work, &
                      nlocaly*nlocalx1*niz, mpi_wp_complex, &
                      qdfyxz, &
                      nlocaly*nlocalx1*niz, mpi_wp_complex, &
                      grid_commy, ierror)
#endif

    !$omp end master
    !$omp barrier

    ! x direction
    !
    do iz = id+1, nlocalz, nthread
      do iy = 1, nlocaly

        do i = 1, nlocalx1, 16
          do ix = i, min0(i+15,nlocalx1)
            work1(ix) = qdfyxz(iy,ix,iz)
          end do
        end do

#ifdef FFTW

#ifdef _SINGLE
        call sfftw_execute_dft_c2r(plan_bx, work1, work3)
#else
        call dfftw_execute_dft_c2r(plan_bx, work1, work3)
#endif

#endif

        do ix = 1, nx
          qdf(ix,iy,iz,1) = work3(ix)
        end do

      end do
    end do

    !$omp barrier

    return

  end subroutine bfft3d_1d_alltoall_x_serial

end module fft3d_mod
