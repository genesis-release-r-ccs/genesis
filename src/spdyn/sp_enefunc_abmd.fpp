!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_enefunc_abmd_mod
!> @brief   define potential energy functions for ABMD
!           Refs: 
!! @authors Junichiro Koga (JK)
!  
!  (c) Copyright 2020 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_enefunc_abmd_mod

  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use sp_dynvars_str_mod
  use sp_dynamics_str_mod
  use mpi_parallel_mod
  use messages_mod
  use constants_mod
  use bspline_mod
  use biweight_mod
  use sp_remd_str_mod
  use sp_bias_str_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public :: update_enefunc_abmd

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_enefunc_abmd
  !> @brief        update B-Spline coeffs 
  !! @authors      JK
  !! @param[in]    dynvars    : dynamics variables information
  !! @param[in]    dynamics   : dynamics information
  !! @param[in]    kbT        : kbT in kcal/mol
  !! @param[inout] enefunc    : potential energy functions information
  !! @param[inout] remd       : REMD information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_enefunc_abmd(dynvars, dynamics, kbT, enefunc, remd, bias)
    type(s_dynvars),  intent(in)            :: dynvars
    type(s_dynamics), intent(in)            :: dynamics
    real(dp),         intent(in)            :: kbT
    type(s_enefunc),  intent(inout), target :: enefunc
    type(s_remd), optional, intent(in)      :: remd
    type(s_bias), optional, intent(in)      :: bias


    type(s_enefunc_reac_coord), pointer :: reac_coord
    type(s_enefunc_abmd),       pointer :: abmd
    real(dp),                   pointer :: uk(:)
    real(wp),                   pointer :: invdxi(:)
    integer,                    pointer :: indices(:)
    integer,                    pointer :: istarts(:), nnonzeros(:)
    real(wip)                           :: dt 
    real(dp)                            :: factor, utmp, wt_factor
    integer                             :: i, j, k, l, m, n, ista_xi, iend_xi
    real(wp)                            :: xi
    real(dp)                            :: avtemp
    integer                             :: nreplicas, nmax_nnonzero, ierr
    integer,                allocatable :: irep_nnonzero(:)           
    integer,                allocatable :: indicesk(:,:)
    integer,                allocatable :: indicesrep(:,:,:)
    real(dp),               allocatable :: ukperrep(:,:)
    logical                             :: out_of_range
    integer                             :: num_degree

    ! local variables
    integer                             :: nextk, indj, nloop, loop_count
    integer                             :: inum, jnum, replicaid, dockid
    integer,                allocatable :: update_ind(:), current_indices(:)
    integer,                allocatable :: current_indicesrep(:,:)
    logical                             :: multiple_walker, exit_periodicity
    logical                             :: remd_use,needs_periodic_copy
    real(dp)                            :: utmp_wt

    if (.not. enefunc%abmd%update_biasing_potential) then
      return
    endif

    out_of_range = .false.
    abmd        => enefunc%abmd
    reac_coord  => enefunc%reac_coord
    uk          => enefunc%abmd%uk 
    indices     => reac_coord%indices
    istarts     => reac_coord%istarts
    nnonzeros   => reac_coord%nnonzeros
    invdxi      => abmd%invdxi

    allocate(update_ind(abmd%nreac_coord))
    allocate(current_indices(abmd%nreac_coord))

    multiple_walker = .false.
    remd_use        = .false.
    
    if (present(remd)) then
      if (remd%total_nreplicas>=1) then
        if (remd%exchange_period==dynamics%nsteps) then 
          multiple_walker = .true.
        else
          remd_use = .true.
        end if
      end if
    end if

    factor   =  (dynamics%timestep) * kbT/enefunc%abmd%tauf

    !if (remd_use) then
    !  call mpi_allreduce(kbT, avtemp, 1, mpi_double_precision, &
    !                     mpi_sum, mpi_comm_country, ierr)
    !                    
    !  avtemp   =  avtemp/real(remd%total_nreplicas,wp)  
    !  factor   =  (dynamics%timestep) * avtemp/enefunc%abmd%tauf
    !else
    !  factor   =  (dynamics%timestep) * kbT/enefunc%abmd%tauf
    !  !factor   =  (dynamics%timestep) * kbT/enefunc%abmd%tauf
    !endif

    if(multiple_walker) then
      nreplicas  =  remd%total_nreplicas
      allocate(irep_nnonzero(nreplicas)); irep_nnonzero = 0

#ifdef HAVE_MPI_GENESIS
      call mpi_allgather(reac_coord%nnonzero, 1, mpi_integer, &
      	                 irep_nnonzero,       1, mpi_integer, &
      	                 mpi_comm_airplane, ierr)
#endif
      nmax_nnonzero = maxval(irep_nnonzero)
      allocate(indicesk(nmax_nnonzero, nreplicas)); indicesk = 0
      allocate(ukperrep(nmax_nnonzero, nreplicas)); ukperrep = 0.0_dp

      ! For periodic CVs with REPLICAs
      allocate(indicesrep(nmax_nnonzero, abmd%nreac_coord, &
                          nreplicas)); indicesrep = 0
      allocate(current_indicesrep(nmax_nnonzero, abmd%nreac_coord)); current_indicesrep = 0      
    endif

    if (present(bias)) dockid = my_country_no + 1

    wt_factor = 1.0_dp
    do i = 1, reac_coord%nnonzero
      call bspline_resolve_indices(i, abmd%nreac_coord,  &
                                   nnonzeros, indices, npad = 0)

      utmp = 1.0_dp
      do j = 1, abmd%nreac_coord
        if (present(bias)) then
          jnum = abmd%dockid2repsetid(dockid)
          inum = abmd%ifunc(bias%repsetid2paramid(jnum,1), j)
          xi  = reac_coord%xi(inum)
        else
          xi  = reac_coord%xi(j)
        end if

        out_of_range = xi>abmd%ximax(j) .or. xi<abmd%ximin(j)
        if (out_of_range) exit
        xi  = (xi-abmd%ximin(j)) * invdxi(j)
        indices (j) = istarts(j) + indices(j) - 1
        utmp        = utmp * biweight(xi-indices(j))
      enddo

      if(.not. out_of_range) then
        k = bspline_get_1d_index(abmd%nreac_coord, indices, abmd%nxi, &
                                 abmd%periodicity)
        if (multiple_walker) then
          indicesk  (i,    my_country_no+1) = k
          ukperrep  (i,    my_country_no+1) = utmp
          indicesrep(i, :, my_country_no+1) = indices
        else
          if (enefunc%abmd%wt) then
            utmp_wt = uk(k) + utmp * factor
            wt_factor =  exp(-utmp_wt / (KBOLTZ * enefunc%abmd%wt_temperature))
          end if
          uk(k) = uk(k) + utmp * factor * wt_factor
      
          !! Check periodicity of CVs for REPLICAs
          if (reac_coord%chk_periodicity) then
            do j = 1, abmd%nreac_coord
              do l = 1, abmd%nreac_coord
                update_ind(l) = l
              end do
        
              do
                current_indices = indices
        
                exit_periodicity = .false.
                do l = 1, j
                  if (.not. abmd%periodicity(update_ind(l))) then
                    exit_periodicity = .true.
                    exit
                  end if
                end do
        
                if (.not. exit_periodicity) then
                  needs_periodic_copy = has_periodic_boundary_index(indices, update_ind, j, abmd%nxi)
                
                  if (needs_periodic_copy) then
                    do l = 1, j
                      current_indices(update_ind(l)) = map_periodic_index( &
                           current_indices(update_ind(l)), abmd%nxi(update_ind(l)) )
                    end do
                
                    nextk = bspline_get_1d_index(abmd%nreac_coord, current_indices, &
                                                 abmd%nxi, abmd%periodicity)
                    uk(nextk) = uk(k)
                  end if
                end if
        
                if (.not. next_combination_prefix(update_ind, j, abmd%nreac_coord)) exit
              end do
            end do
          end if
          !! End of checking periodicity of CVs

        end if
      endif
    enddo

    if (multiple_walker) then
#ifdef HAVE_MPI_GENESIS      
      call mpi_allreduce(mpi_in_place, indicesk, nreplicas*nmax_nnonzero, &
                         mpi_integer, mpi_sum, mpi_comm_airplane, ierr)
      call mpi_allreduce(mpi_in_place, ukperrep, nreplicas*nmax_nnonzero, &
                         mpi_double_precision, mpi_sum, mpi_comm_airplane, ierr)
      call mpi_allreduce(mpi_in_place, indicesrep, nreplicas*nmax_nnonzero*abmd%nreac_coord, &
                         mpi_integer, mpi_sum, mpi_comm_airplane, ierr)
#endif                          
      do i = 1, nreplicas
        do j = 1, irep_nnonzero(i)
          if (enefunc%abmd%wt) then
            utmp_wt = uk(indicesk(j, i)) + ukperrep(j, i) * factor
            wt_factor =  exp(-utmp_wt / (KBOLTZ * enefunc%abmd%wt_temperature))
          end if
          uk(indicesk(j, i)) = uk(indicesk(j, i)) + ukperrep(j, i) * factor * wt_factor

          !! Check periodicity of CVs for REPLICAs
          if (reac_coord%chk_periodicity) then
            do n = 1, abmd%nreac_coord
              do l = 1, abmd%nreac_coord
                update_ind(l) = l
              end do

              do
                current_indicesrep(j,:) = indicesrep(j,:,i)

                exit_periodicity = .false.
                do l = 1, n
                  if (.not. abmd%periodicity(update_ind(l))) then
                    exit_periodicity = .true.
                    exit
                  end if
                end do

                if (.not. exit_periodicity) then
                  needs_periodic_copy = has_periodic_boundary_index( &
                       indicesrep(j,:,i), update_ind, n, abmd%nxi)

                  if (needs_periodic_copy) then
                    do l = 1, n
                      current_indicesrep(j,update_ind(l)) = map_periodic_index( &
                           current_indicesrep(j,update_ind(l)), abmd%nxi(update_ind(l)) )
                    end do

                    nextk = bspline_get_1d_index(abmd%nreac_coord, current_indicesrep(j,:), &
                                                 abmd%nxi, abmd%periodicity)
                    uk(nextk) = uk(indicesk(j, i))
                  end if
                end if

                if (.not. next_combination_prefix(update_ind, n, abmd%nreac_coord)) exit
              end do
            end do
          end if
          !! End of checking periodicity of CVs

        enddo
      enddo
 
      deallocate(irep_nnonzero)
      deallocate(indicesk)
      deallocate(ukperrep)
      deallocate(indicesrep)
      deallocate(current_indicesrep)
    endif

    deallocate(update_ind)
    deallocate(current_indices)

    return

  end subroutine update_enefunc_abmd

  integer function map_periodic_index(idx, nxi) result(new_idx)
    implicit none
    integer, intent(in) :: idx, nxi

    if (idx == -1) then
      new_idx = nxi - 2
    else if (idx == 0) then
      new_idx = nxi - 1
    else if (idx == 1) then
      new_idx = nxi
    else if (idx == nxi - 2) then
      new_idx = -1
    else if (idx == nxi - 1) then
      new_idx = 0
    else if (idx == nxi) then
      new_idx = 1
    else
      new_idx = idx
    end if
  end function map_periodic_index

  logical function next_combination_prefix(update_ind, j, n) result(has_next)
    implicit none
    integer, intent(inout) :: update_ind(:)
    integer, intent(in)    :: j, n
    integer :: l, m

    do l = j, 1, -1
      if (update_ind(l) < n - j + l) then
        update_ind(l) = update_ind(l) + 1
        do m = l + 1, j
          update_ind(m) = update_ind(m - 1) + 1
        end do
        has_next = .true.
        return
      end if
    end do

    has_next = .false.
  end function next_combination_prefix

  logical function has_periodic_boundary_index(indices, update_ind, j, nxi) result(found)
    implicit none
    integer, intent(in) :: indices(:)
    integer, intent(in) :: update_ind(:)
    integer, intent(in) :: j
    integer, intent(in) :: nxi(:)
    integer :: l, idx, idim  

    found = .false.  

    do l = 1, j
      idim = update_ind(l)
      idx  = indices(idim)  

      if (idx == -1 .or. idx == 0 .or. idx == 1 .or. &
          idx == nxi(idim) - 2 .or. idx == nxi(idim) - 1 .or. idx == nxi(idim)) then
        found = .true.
        return
      end if
    end do
  end function has_periodic_boundary_index

end module sp_enefunc_abmd_mod

