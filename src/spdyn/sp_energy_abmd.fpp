!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_abmd_mod
!> @brief   compute energy for ABMD
!           Refs: 
!! @authors Junichiro Koga (JK)
!  
!  (c) Copyright 2020 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_abmd_mod

  use sp_enefunc_str_mod
  use sp_energy_str_mod
  use sp_domain_str_mod
  use sp_boundary_str_mod
  use sp_bias_str_mod
  use fileio_control_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use math_libs_mod
  use string_mod
  use bspline_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public :: compute_bias_potential_abmd

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_bias_potential_abmd
  !> @brief        add the bias potential energy for the ABMD method
  !! @authors      JK
  !! @param[in]    domain     : domain information
  !! @param[inout] enefunc    : potential energy functions information
  !! @param[inout] energy     : energy information
  !! @param[inout] force      : forces of target systems
  !! @param[inout] virial     : virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_bias_potential_abmd(domain, enefunc, energy, force, &
                                         virial, sumu, xiex, update_u)

    ! formal arguments
    type(s_domain),  target,     intent(in)    :: domain
    type(s_enefunc), target,     intent(inout) :: enefunc
    type(s_energy), optional,    intent(inout) :: energy
    real(wip), optional,         intent(inout) :: force(:,:,:)
    real(dp),  optional,         intent(inout) :: virial(3,3)
    real(dp),  optional,         intent(out)   :: sumu
    real(wp),  optional, target, intent(in)    :: xiex(:)
    logical,   optional,         intent(in)    :: update_u
    !type(s_bias), optional,      intent(inout) :: bias

    type(s_enefunc_abmd),       pointer :: abmd
    type(s_enefunc_reac_coord), pointer :: reac_coord

    integer                     :: i, j, k, id, iatm, ix, icel, ncell_local, next, ifor
    integer                     :: repi, inum, jnum
    integer                     :: dockid 
    real(wp)                    :: xi
    real(wp),       pointer     :: xis(:)
    real(wp),       pointer     :: bk(:,:), bk2(:,:), invdxi(:), sumdbk(:)
    real(wp)                    :: pbk, pdbk, udiff
    real(dp)                    :: utmp, uoff
    real(dp),       pointer     :: uk(:)
    integer                     :: ista_xi, iend_xi
    integer(int2),  pointer     :: id_g2l(:,:)
    integer,        pointer     :: bonds_to_atom(:)
    integer,        pointer     :: indices(:)
    integer,        pointer     :: istarts(:), nnonzeros(:)
    integer                     :: ind
    integer                     :: iscomb, iecomb, iwork, ierr
    logical                     :: out_of_range, overlap

    abmd          => enefunc%abmd
    reac_coord    => enefunc%reac_coord
    abmd%u        =  0.0_dp
    if(present(sumu)) &
    sumu          =  0.0_dp
    utmp          =  0.0_dp
    id_g2l        => domain%id_g2l
    ncell_local   =  domain%num_cell_local
    indices       => reac_coord%indices
    istarts       => reac_coord%istarts
    nnonzeros     => reac_coord%nnonzeros
    invdxi        => abmd%invdxi
    bk            => abmd%bk
    bk2           => abmd%bk2
    sumdbk        => abmd%sumdbk

    if(present(xiex)) then
      xis         => xiex
      uk          => abmd%puk
    else
      xis         => reac_coord%xi
      uk          => abmd%uk
      abmd%puk    =  abmd%uk
    endif
    out_of_range  = .false.
    bk  = 0.0_wp
    bk2 = 0.0_wp

    if (enefunc%bias_use) dockid = my_country_no + 1

    do j = 1, abmd%nreac_coord
      inum = j
      if (enefunc%bias_use) then
        jnum = abmd%dockid2repsetid(dockid)
        inum = abmd%ifunc(abmd%repsetid2paramid(jnum,1), j)
      end if
      xi  = xis(inum)
      out_of_range = xi>abmd%ximax(j) .or. xi<abmd%ximin(j)
      if (out_of_range) exit
      xi  = (xi-abmd%ximin(j)) * invdxi(j)
      call bspline(-1, abmd%nxi(j),   xi, bk (:,j))
      call bspline(-1, abmd%nxi(j)-1, xi, bk2(:,j), 2)
    enddo

    reac_coord%nnonzero = 1

    if(.not.out_of_range) then
      do j = 1, abmd%nreac_coord
        inum = j
        if (enefunc%bias_use) then
          jnum = abmd%dockid2repsetid(dockid)
          inum = abmd%ifunc(abmd%repsetid2paramid(jnum,1), j)
        end if
        xi                  = xis(inum)
        xi                  = (xi-abmd%ximin(j)) * invdxi(j)
        ista_xi             = max(int(ceiling(xi - 2)), -1)
        iend_xi             = min(int(floor  (xi + 2)), abmd%nxi(j))
        istarts(j)          = ista_xi
        nnonzeros(j)        = (iend_xi - ista_xi + 1)
        reac_coord%nnonzero = reac_coord%nnonzero * nnonzeros(j)
      enddo

      iwork  = (reac_coord%nnonzero-1)/nproc_country + 1
      iscomb = min(my_country_rank*iwork+1, reac_coord%nnonzero+1)
      iecomb = min(iscomb+iwork-1,reac_coord%nnonzero)

        do i=1,reac_coord%nnonzero
          call bspline_resolve_indices(i, abmd%nreac_coord,  &
                                       nnonzeros, indices, npad = 0)
          pbk  = 1.0_wp
          uoff = 0.0_wp
          do j = 1, abmd%nreac_coord
            indices(j)  = istarts(j) + indices(j) - 1
            pbk         = pbk * bk(indices(j),j)
            !if (indices(j)==-1 .or. indices(j)==reac_coord%nxi(j)) then
            !  uoff = uoff + reac_coord%uboundary(j)
            !endif
          enddo
          k    = bspline_get_1d_index(abmd%nreac_coord, &
                                      indices, abmd%nxi, &
                                      abmd%periodicity)
          utmp = utmp + pbk * (uk(k)+uoff)
        enddo
#ifdef HAVE_MPI_GENESIS
      call mpi_bcast(utmp, 1, mpi_double_precision, 0, mpi_comm_country, ierr)
#endif
    endif

    if(present(update_u)) then
      if(update_u) then
        abmd%u = utmp
      endif
    else
      abmd%u = utmp
    endif

    if(present(sumu)) sumu = utmp

    if( .not. present(force) .or. .not. present(virial)) then
      return
    endif

    if(.not.out_of_range) then
      reac_coord%nnonzero = 1
      do j = 1, abmd%nreac_coord
        inum = j
        if (enefunc%bias_use) then
          jnum = abmd%dockid2repsetid(dockid)
          inum = abmd%ifunc(abmd%repsetid2paramid(jnum,1), j)
        end if
        xi                  = xis(inum)
        xi                  = (xi-abmd%ximin(j)) * invdxi(j)
        ista_xi             = max(int(ceiling(xi - 2)), -1)
        iend_xi             = min(int(floor  (xi + 2)), abmd%nxi(j)-1)
        istarts(j)          = ista_xi
        nnonzeros(j)        = (iend_xi - ista_xi + 1)
        reac_coord%nnonzero = reac_coord%nnonzero * nnonzeros(j)
      enddo

      iwork  = (reac_coord%nnonzero-1)/nproc_country + 1
      iscomb = min(my_country_rank*iwork+1, reac_coord%nnonzero+1)
      iecomb = min(iscomb+iwork-1,reac_coord%nnonzero)
      sumdbk = 0.0_wp

      do i = 1, reac_coord%nnonzero
        call bspline_resolve_indices(i, abmd%nreac_coord,  &
                                     nnonzeros, indices, npad = 0)
        do j = 1, abmd%nreac_coord
          indices (j) = istarts(j) + indices(j) - 1
        enddo
        ind  = bspline_get_1d_index(abmd%nreac_coord, indices, &
                                    abmd%nxi, abmd%periodicity)
        do j = 1, abmd%nreac_coord
          next  = ind + bspline_get_next_index_offset(j, abmd%nreac_coord, &
                                                      abmd%nxi)
          uoff = 0.0_wp
          if( indices(j) == -1 ) then
            !uoff = -reac_coord%uboundary(j)
          endif
          if( indices(j) == reac_coord%nxi(j)-1 ) then
            !uoff =  reac_coord%uboundary(j)
          endif
          udiff = uk(next) - uk(ind) + uoff
          pdbk  = udiff * bk2(indices(j),j) * invdxi(j)
          do k = 1, abmd%nreac_coord
            if(k == j) cycle
            pdbk = pdbk * bk(indices(k),k)
          enddo
          sumdbk(j) = sumdbk(j) + pdbk
        enddo
      enddo
#ifdef HAVE_MPI_GENESIS
      call mpi_bcast(sumdbk, abmd%nreac_coord, mpi_wp_real, 0, mpi_comm_country, ierr)
#endif
      do i = 1, abmd%nreac_coord
        inum = i
        if (enefunc%bias_use) then
          jnum = abmd%dockid2repsetid(dockid)
          inum = abmd%ifunc(abmd%repsetid2paramid(jnum,1), i)
        end if
        bonds_to_atom => abmd%bondslist_to_atomlist(:,i)
        do iatm = 1, reac_coord%num_atoms_bonds_reac_coord(inum)
          icel = id_g2l(1,bonds_to_atom(iatm))          
          overlap = .false.

          do j = iatm+1, reac_coord%num_atoms_bonds_reac_coord(inum)
            if (bonds_to_atom(iatm) == bonds_to_atom(j)) then
              overlap = .true.
              exit
            end if 
          end do
        
          if (icel > 0 .and. icel <= ncell_local) then
            if (.not. overlap) then
              ix = id_g2l(2,bonds_to_atom(iatm))
              force(1:3,ix,icel) = &
              force(1:3,ix,icel) + &
              abmd%force_xi(i,1:3,ix,icel) * sumdbk(i)
            end if
          end if
        enddo
        if(main_rank .or. replica_main_rank) then
          virial(:,:) = virial(:,:) + abmd%virial_xi(i,:,:) * sumdbk(i)
        endif
      enddo
    endif

    if(present(energy)) then
      energy%bias_abmd = abmd%u
    endif

    return

  end subroutine compute_bias_potential_abmd

end module sp_energy_abmd_mod
