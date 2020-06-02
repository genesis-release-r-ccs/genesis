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
#ifdef MPI
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

#ifdef PKTIMER
    use Ctim
#endif

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
    integer*8         :: plan
    integer           :: ix, iy, iz, izs
    integer           :: izy, izx, iyx
    integer           :: iz_start, iz_end
    integer           :: k, iproc, ierror
    integer           :: iprocx, iprocy

#ifdef PKTIMER
    real(dp)         :: st,et
#endif

    !$omp barrier
#ifdef PKTIMER
    !$omp master
    call timer_sta(23)
    !$omp end master
#ifdef FJ_TIMER_DETAIL
    !$omp master
    call timer_sta(139)
    !$omp end master
#ifdef FJ_PROF_FAPP
    call fapp_start("FFT3D_fft3d_2d_alltoall_1",139,0)
#endif
#else
#ifdef FJ_PROF_FAPP
    call fapp_start("FFT3D",23,0)
#endif
#endif
#endif

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

#ifdef PKTIMER
    !$omp master
    call timer_end(23)
    !$omp end master
#ifdef FJ_TIMER_DETAIL
    !$omp master
    call timer_end(139)
    !$omp end master
#ifdef FJ_PROF_FAPP
    call fapp_stop("FFT3D_fft3d_2d_alltoall_1",139,0)
#endif
#else
#ifdef FJ_PROF_FAPP
    call fapp_stop("FFT3D",23,0)
#endif
#endif
#endif

    !$omp master
#ifdef MPI

#ifdef PKTIMER
    call gettod(st)
    call mpi_barrier(grid_commxy,ierror)
    call gettod(et)
    mpi_bari(1)=mpi_bari(1)+(et-st)
    call gettod(st)
#endif

    call mpi_alltoall(qdfxyz,      nlocaly*nix1*niz, mpi_wp_complex, &
                      qdfxyz_work, nlocaly*nix1*niz, mpi_wp_complex, &
                      grid_commxy, ierror)

#ifdef PKTIMER
    call gettod(et)
    mpi_tran(1,2)=mpi_tran(1,2)+(et-st)
#endif

#endif
    !$omp end master
    !$omp barrier

#ifdef PKTIMER
    !$omp master
    call timer_sta(23)
    !$omp end master
#ifdef FJ_TIMER_DETAIL
    !$omp master
    call timer_sta(140)
    !$omp end master
#ifdef FJ_PROF_FAPP
    call fapp_start("FFT3D_fft3d_2d_alltoall_2",140,0)
#endif
#else
#ifdef FJ_PROF_FAPP
    call fapp_start("FFT3D",23,0)
#endif
#endif
#endif

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

#ifdef PKTIMER
    !$omp master
    call timer_end(23)
    !$omp end master
#ifdef FJ_TIMER_DETAIL
    !$omp master
    call timer_end(140)
    !$omp end master
#ifdef FJ_PROF_FAPP
    call fapp_stop("FFT3D_fft3d_2d_alltoall_2",140,0)
#endif
#else
#ifdef FJ_PROF_FAPP
    call fapp_stop("FFT3D",23,0)
#endif
#endif
#endif

    !$omp master
#ifdef MPI

#ifdef PKTIMER
    call gettod(st)
    call mpi_barrier(grid_commz,ierror)
    call gettod(et)
    mpi_bari(2)=mpi_bari(2)+(et-st)
#endif

    call mpi_alltoall(qdfyzx, nix1*niy*nlocalz, mpi_wp_complex, &
                      qdfyzx_procz, nix1*niy*nlocalz, mpi_wp_complex, &
                      grid_commz, ierror)

#ifdef PKTIMER
    call gettod(et)
    mpi_tran(2,2)=mpi_tran(2,2)+(et-st)
#endif

#endif
    !$omp end master
    !$omp barrier

#ifdef PKTIMER
    !$omp master
    call timer_sta(23)
    !$omp end master
#ifdef FJ_TIMER_DETAIL
    !$omp master
    call timer_sta(141)
    !$omp end master
#ifdef FJ_PROF_FAPP
    call fapp_start("FFT3D_fft3d_2d_alltoall_3",141,0)
#endif
#else
#ifdef FJ_PROF_FAPP
    call fapp_start("FFT3D",23,0)
#endif
#endif
#endif
 
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

#ifdef PKTIMER
    !$omp master
    call timer_end(23)
    !$omp end master
#ifdef FJ_TIMER_DETAIL
    !$omp master
    call timer_end(141)
    !$omp end master
#ifdef FJ_PROF_FAPP
    call fapp_stop("FFT3D_fft3d_2d_alltoall_3",141,0)
#endif
#else
#ifdef FJ_PROF_FAPP
    call fapp_stop("FFT3D",23,0)
#endif
#endif
#endif

    return

  end subroutine fft3d_2d_alltoall

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bfft3d_2d_alltoall
  !> @brief        calculate inverse 3D fast fourier transform 
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine bfft3d_2d_alltoall(qdf, qdf_real, qdfxyz, qdfzyx_work,       &
                        qdfyxz_procy, qdfxyz_work, qdfzyx, c_work,        &
                        work1, work2,                                     &
                        nlocalx, nlocaly, nlocalz,                        &
                        nx, ny, nz, nprocx, nprocy, nprocz,     &
                        nix, nix1, niy, niz, id, nthread,                 &
                        grid_commx,       &
                        grid_commz, grid_commxy)

#ifdef PKTIMER
    use Ctim
#endif

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

#ifdef PKTIMER
    real(dp)         :: st,et
#endif

    !$omp barrier
#ifdef PKTIMER
    !$omp master
    call timer_sta(23)
    !$omp end master
#ifdef FJ_TIMER_DETAIL
    !$omp master
    call timer_sta(142)
    !$omp end master
#ifdef FJ_PROF_FAPP
    call fapp_start("FFT3D_bfft3d_2d_alltoall_1",142,0)
#endif
#else
#ifdef FJ_PROF_FAPP
    call fapp_start("FFT3D",23,0)
#endif
#endif
#endif

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

#ifdef PKTIMER
    !$omp master
    call timer_end(23)
    !$omp end master
#ifdef FJ_TIMER_DETAIL
    !$omp master
    call timer_end(142)
    !$omp end master
#ifdef FJ_PROF_FAPP
    call fapp_stop("FFT3D_bfft3d_2d_alltoall_1",142,0)
#endif
#else
#ifdef FJ_PROF_FAPP
    call fapp_stop("FFT3D",23,0)
#endif
#endif
#endif

    !$omp master

#ifdef PKTIMER
    call gettod(st)
    call mpi_barrier(grid_commz,ierror)
    call gettod(et)
    mpi_bari(4)=mpi_bari(4)+(et-st)
    call gettod(st)
#endif

#ifdef MPI
    call mpi_alltoall(qdfzyx_work, nix1*niy*nlocalz, mpi_wp_complex, &
                      qdfzyx,      nix1*niy*nlocalz, mpi_wp_complex, &
                      grid_commz, ierror)

#ifdef PKTIMER
    call gettod(et)
    mpi_tran(5,2)=mpi_tran(5,2)+(et-st)
#endif

#endif
    !$omp end master
    !$omp barrier

    ! y direction
    !

#ifdef PKTIMER
    !$omp master
    call timer_sta(23)
    !$omp end master
#ifdef FJ_TIMER_DETAIL
    !$omp master
    call timer_sta(143)
    !$omp end master
#ifdef FJ_PROF_FAPP
    call fapp_start("FFT3D_bfft3d_2d_alltoall_2",143,0)
#endif
#else
#ifdef FJ_PROF_FAPP
    call fapp_start("FFT3D",23,0)
#endif
#endif
#endif

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

#ifdef PKTIMER
    !$omp master
    call timer_end(23)
    !$omp end master
#ifdef FJ_TIMER_DETAIL
    !$omp master
    call timer_end(143)
    !$omp end master
#ifdef FJ_PROF_FAPP
    call fapp_stop("FFT3D_bfft3d_2d_alltoall_2",143,0)
#endif
#else
#ifdef FJ_PROF_FAPP
    call fapp_stop("FFT3D",23,0)
#endif
#endif
#endif

    !$omp master
#ifdef MPI

#ifdef PKTIMER
    call gettod(st)
    call mpi_barrier(grid_commxy,ierror)
    call gettod(et)
    mpi_bari(5)=mpi_bari(5)+(et-st)
    call gettod(st)
#endif

    call mpi_alltoall(qdfxyz_work, nlocaly*nix1*niz, mpi_wp_complex, &
                      qdfxyz,      nlocaly*nix1*niz, mpi_wp_complex, &
                      grid_commxy, ierror)

#ifdef PKTIMER
    call gettod(et)
    mpi_tran(6,2)=mpi_tran(6,2)+(et-st)
#endif

#endif
    !$omp end master
    !$omp barrier

    ! x direction
    !
#ifdef PKTIMER
    !$omp master
    call timer_sta(23)
    !$omp end master
#ifdef FJ_TIMER_DETAIL
    !$omp master
    call timer_sta(144)
    !$omp end master
#ifdef FJ_PROF_FAPP
    call fapp_start("FFT3D_bfft3d_2d_alltoall_3",144,0)
#endif
#else
#ifdef FJ_PROF_FAPP
    call fapp_start("FFT3D",23,0)
#endif
#endif
#endif

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

#ifdef PKTIMER
    !$omp master
    call timer_end(23)
    !$omp end master
#ifdef FJ_TIMER_DETAIL
    !$omp master
    call timer_end(144)
    !$omp end master
#ifdef FJ_PROF_FAPP
    call fapp_stop("FFT3D_bfft3d_2d_alltoall_3",144,0)
#endif
#else
#ifdef FJ_PROF_FAPP
    call fapp_stop("FFT3D",23,0)
#endif
#endif
#endif

    !$omp master
#ifdef MPI

#ifdef PKTIMER
    call gettod(st)
    call mpi_barrier(grid_commx,ierror)
    call gettod(et)
    mpi_bari(6)=mpi_bari(6)+(et-st)
    call gettod(st)
#endif

    call mpi_alltoall(qdf,      nlocalx*nlocaly*niz, mpi_wp_real, &
                      qdf_real, nlocalx*nlocaly*niz, mpi_wp_real, &
                      grid_commx, ierror)

#ifdef PKTIMER
    call gettod(et)
    mpi_tran(7,2)=mpi_tran(7,2)+(et-st)
#endif

#endif
    !$omp end master
    !$omp barrier

    return

  end subroutine bfft3d_2d_alltoall


end module fft3d_2dalltoall_mod
