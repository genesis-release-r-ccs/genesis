!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_energy_morph_mod
!> @brief   calculate morph energy
!! @authors Chigusa Kobayashi (CK)
!! @date    2015/02/06 (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
  
module at_energy_morph_mod

  use at_enefunc_str_mod
  use timers_mod
  use mpi_parallel_mod
  use messages_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public :: compute_energy_morph

contains
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine  compute_energy_morph
  !> @brief      calculate morph energy
  !! @authors    CK
  !! @param[in]  enefunc potential energy functions [str]
  !! @param[in]  coord coordinates of target systems [dble]
  !! @param[out] force forces of target systems [dble]
  !! @param[out] virial of target systems [dble]
  !! @param[out] morpgh energy of target systems [dble]
  !! @date       2015/02/06 (CK)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_morph(enefunc, coord, force, virial, emorph, cv)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: emorph
    real(wp),                intent(inout) :: cv(1:2)

    ! local variables
    integer                  :: i, j, k
    integer                  :: istart_bb, iend_bb
    integer                  :: istart_sc, iend_sc
    integer                  :: i1, i2
    integer                  :: num_morph_bb, num_morph_sc
    real(wp)                 :: d12(1:3), r12, vtmp
    real(wp)                 :: t, tmp, dn, coef, coef_t, tott
    real(wp)                 :: vxx, vyx, vzx, vyy, vzy, vzz
    real(wp)                 :: work(3)
    integer,         pointer :: list_bb(:,:), list_sc(:,:)
    real(wp),        pointer :: rmin_bb(:), rmin_sc(:)


    if (enefunc%morph_ene_flag == MorphEneFlagBB) then

      call compute_energy_morph_bb(enefunc, coord, force, virial, emorph, cv)

    else if (enefunc%morph_ene_flag == MorphEneFlagSC) then
    
      call compute_energy_morph_sc(enefunc, coord, force, virial, emorph, cv)

    else if (enefunc%morph_ene_flag == MorphEneFlagTot) then

      call compute_energy_morph_tot(enefunc, coord, force, virial, emorph, cv)

    end if

    return

  end subroutine compute_energy_morph

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine  compute_energy_morph_tot
  !> @brief      calculate morph energy
  !! @authors    CK
  !! @param[in]  enefunc potential energy functions [str]
  !! @param[in]  coord coordinates of target systems [dble]
  !! @param[out] force forces of target systems [dble]
  !! @param[out] virial of target systems [dble]
  !! @param[out] morpgh energy of target systems [dble]
  !! @date       2015/02/06 (CK)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_morph_tot(enefunc, coord, force, virial, emorph, cv)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: emorph
    real(wp),                intent(inout) :: cv(1:2)

    ! local variables
    integer                  :: i, j, k, id, omp_get_thread_num
    integer                  :: istart_bb, iend_bb
    integer                  :: istart_sc, iend_sc
    integer                  :: i1, i2
    integer                  :: num_morph_bb, num_morph_sc
    real(wp)                 :: d12(1:3), r12, vtmp
    real(wp)                 :: t, tt(2), tmp, dn, coef, coef_t, tott(2)
    real(wp)                 :: vxx, vyx, vzx, vyy, vzy, vzz
    real(wp)                 :: work(3)
    integer,         pointer :: list_bb(:,:), list_sc(:,:)
    real(wp),        pointer :: rmin_bb(:), rmin_sc(:)
    real(wp),        pointer :: rmin_bb_other(:), rmin_sc_other(:)


    call timer(TimerMorph, TimerOn)

    istart_bb = enefunc%istart_morph_bb
    iend_bb   = enefunc%iend_morph_bb

    list_bb   => enefunc%morph_list_bb
    rmin_bb   => enefunc%morph_dist_bb
    rmin_bb_other   => enefunc%morph_dist_bb_other

    istart_sc = enefunc%istart_morph_sc
    iend_sc   = enefunc%iend_morph_sc

    list_sc   => enefunc%morph_list_sc
    rmin_sc   => enefunc%morph_dist_sc
    rmin_sc_other   => enefunc%morph_dist_sc_other

    num_morph_bb = enefunc%num_morph_bb
    num_morph_sc = enefunc%num_morph_sc


    emorph  = 0.0_wp
    t       = 0.0_wp
    tt(1:2) = 0.0_wp
    dn = real(num_morph_bb+num_morph_sc, wp)

    !$omp parallel do default(none)                          &
    !$omp private(i,  d12, r12,  tmp, i1, i2)    &
    !$omp shared(istart_bb, iend_bb, coord, list_bb, rmin_bb, num_morph_bb, &
    !$omp        rmin_bb_other)   &
    !$omp reduction(+:tt)
    !

    do i = istart_bb, iend_bb

      i1 = list_bb(1,i)
      i2 = list_bb(2,i)
      d12(1:3) = coord(1:3,i1) - coord(1:3,i2)
      r12   = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )
      tmp = r12 - rmin_bb(i)
      tt(1) = tt(1) + tmp*tmp
      tmp = r12 - rmin_bb_other(i)
      tt(2) = tt(2) + tmp*tmp

    end do
    !$omp end parallel do

    !$omp parallel do default(none)                          &
    !$omp private(i,  d12, r12,  tmp, i1, i2)    &
    !$omp shared(istart_sc, iend_sc, coord, list_sc, rmin_sc, num_morph_sc,   &
    !$omp        rmin_sc_other)   &
    !$omp reduction(+:tt)
    !

    do i = istart_sc, iend_sc

      i1 = list_sc(1,i)
      i2 = list_sc(2,i)
      d12(1:3) = coord(1:3,i1) - coord(1:3,i2)
      r12   = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )
      tmp = r12 - rmin_sc(i)
      tt(1) = tt(1) + tmp*tmp
      tmp = r12 - rmin_sc_other(i)
      tt(2) = tt(2) + tmp*tmp

    end do
    !$omp end parallel do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(tt, tott, 2, &
                    mpi_wp_real, mpi_sum, mpi_comm_country, ierror)
#else
    tott(1:2) = tt(1:2)
#endif
    cv(1:2) = sqrt(tott(1:2)/dn)

    t  = sqrt(tott(1)/dn)
    if (main_rank) &
      emorph = enefunc%morph_spring*t*t / (enefunc%morph_linear+t)
    coef_t = (1.0_wp/dn)*enefunc%morph_spring*(2.0d0*enefunc%morph_linear+t)/  &
              ((enefunc%morph_linear+t)*(enefunc%morph_linear+t))


    !$omp parallel default(none)                          &
    !$omp private(i, j, k, d12, r12, tmp, vtmp, coef, work, i1, i2,id)    &
    !$omp shared(istart_bb, iend_bb, coord, list_bb,coef_t, rmin_bb, dn, t, &
    !$omp        force, nthread)     &
    !$omp reduction(+:virial)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = istart_bb+id, iend_bb, nthread
      i1 = list_bb(1,i)
      i2 = list_bb(2,i)
      d12(1:3) = coord(1:3,i1) - coord(1:3,i2)
      r12   = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )
      tmp = r12 - rmin_bb(i)
      coef = coef_t*tmp/r12

      work(1:3) = coef * d12(1:3)

      ! virial
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = d12(k) * work(j) 
          virial(k,j) = virial(k,j) - vtmp
          virial(j,k) = virial(j,k) - vtmp
        end do
        vtmp = d12(j) * work(j) 
        virial(j,j) = virial(j,j) - vtmp
      end do
      force(1:3,i1,id+1) = force(1:3,i1,id+1) - work(1:3)
      force(1:3,i2,id+1) = force(1:3,i2,id+1) + work(1:3)
    end do

    !$omp end parallel

    !$omp parallel default(none)                          &
    !$omp private(i, j, k, d12, r12, tmp, vtmp, coef, work, i1, i2, id)    &
    !$omp shared(istart_sc, iend_sc, coord, list_sc,coef_t, rmin_sc, dn, t,&
    !$omp        force, nthread)     &
    !$omp reduction(+:virial)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = istart_sc+id, iend_sc, nthread
      i1 = list_sc(1,i)
      i2 = list_sc(2,i)
      d12(1:3) = coord(1:3,i1) - coord(1:3,i2)
      r12   = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )
      tmp = r12 - rmin_sc(i)
      coef = coef_t*tmp/r12

      work(1:3) = coef * d12(1:3)

      ! virial
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = d12(k) * work(j) 
          virial(k,j) = virial(k,j) - vtmp
          virial(j,k) = virial(j,k) - vtmp
        end do
        vtmp = d12(j) * work(j) 
        virial(j,j) = virial(j,j) - vtmp
      end do
      force(1:3,i1,id+1) = force(1:3,i1,id+1) - work(1:3)
      force(1:3,i2,id+1) = force(1:3,i2,id+1) + work(1:3)
    end do

    !$omp end parallel

    call timer(TimerMorph, TimerOff)

    return

  end subroutine compute_energy_morph_tot

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine  compute_energy_morph_bb
  !> @brief      calculate morph energy
  !! @authors    CK
  !! @param[in]  enefunc potential energy functions [str]
  !! @param[in]  coord coordinates of target systems [dble]
  !! @param[out] force forces of target systems [dble]
  !! @param[out] virial of target systems [dble]
  !! @param[out] morpgh energy of target systems [dble]
  !! @date       2015/02/06 (CK)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_morph_bb(enefunc, coord, force, virial, emorph, cv)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: emorph
    real(wp),                intent(inout) :: cv(1:2)

    ! local variables
    integer                  :: i, j, k, id, omp_get_thread_num
    integer                  :: istart_bb, iend_bb
    integer                  :: i1, i2
    integer                  :: num_morph_bb
    real(wp)                 :: d12(1:3), r12, vtmp
    real(wp)                 :: t, tt(2), tmp, dn, coef, coef_t, tott(2)
    real(wp)                 :: vxx, vyx, vzx, vyy, vzy, vzz
    real(wp)                 :: work(3)
    integer,         pointer :: list_bb(:,:)
    real(wp),        pointer :: rmin_bb(:)
    real(wp),        pointer :: rmin_bb_other(:)


    call timer(TimerMorph, TimerOn)

    istart_bb = enefunc%istart_morph_bb
    iend_bb   = enefunc%iend_morph_bb
    list_bb   => enefunc%morph_list_bb
    rmin_bb   => enefunc%morph_dist_bb
    rmin_bb_other   => enefunc%morph_dist_bb_other
    num_morph_bb = enefunc%num_morph_bb

    emorph  = 0.0_wp
    t       = 0.0_wp
    tt(1:2) = 0.0_wp
    dn = real(num_morph_bb, wp)

    !$omp parallel do default(none)                          &
    !$omp private(i,  d12, r12,  tmp, i1, i2)    &
    !$omp shared(istart_bb, iend_bb, coord, list_bb, rmin_bb, num_morph_bb, &
    !$omp rmin_bb_other)   &
    !$omp reduction(+:tt) 
    !

    do i = istart_bb, iend_bb

      i1 = list_bb(1,i)
      i2 = list_bb(2,i)
      d12(1:3) = coord(1:3,i1) - coord(1:3,i2)
      r12   = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )
      tmp = r12 - rmin_bb(i)
      tt(1) = tt(1) + tmp*tmp
      tmp = r12 - rmin_bb_other(i)
      tt(2) = tt(2) + tmp*tmp

    end do
    !$omp end parallel do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(tt, tott, 2, &
                    mpi_wp_real, mpi_sum, mpi_comm_country, ierror)
#else
    tott(1:2) = tt(1:2)
#endif
    cv(1:2) = sqrt(tott(1:2)/dn)

    t  = sqrt(tott(1)/dn)
    if (main_rank) &
      emorph = enefunc%morph_spring*t*t / (enefunc%morph_linear+t)

    coef_t = (1.0_wp/dn)*enefunc%morph_spring*(2.0d0*enefunc%morph_linear+t)/  &
              ((enefunc%morph_linear+t)*(enefunc%morph_linear+t))


    !$omp parallel default(none)                          &
    !$omp private(i, j, k, d12, r12, tmp, vtmp, coef, work, i1, i2, id)    &
    !$omp shared(istart_bb, iend_bb, coord, list_bb,coef_t, rmin_bb, dn, t,&
    !$omp        force, nthread)     &
    !$omp reduction(+:virial)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = istart_bb+id, iend_bb, nthread
      i1 = list_bb(1,i)
      i2 = list_bb(2,i)
      d12(1:3) = coord(1:3,i1) - coord(1:3,i2)
      r12   = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )
      tmp = r12 - rmin_bb(i)
      coef = coef_t*tmp/r12

      work(1:3) = coef * d12(1:3)

      ! virial
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = d12(k) * work(j) 
          virial(k,j) = virial(k,j) - vtmp
          virial(j,k) = virial(j,k) - vtmp
        end do
        vtmp = d12(j) * work(j) 
        virial(j,j) = virial(j,j) - vtmp
      end do
      force(1:3,i1,id+1) = force(1:3,i1,id+1) - work(1:3)
      force(1:3,i2,id+1) = force(1:3,i2,id+1) + work(1:3)
    end do

    !$omp end parallel 

    call timer(TimerMorph, TimerOff)

    return

  end subroutine compute_energy_morph_bb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine  compute_energy_morph_sc
  !> @brief      calculate morph energy
  !! @authors    CK
  !! @param[in]  enefunc potential energy functions [str]
  !! @param[in]  coord coordinates of target systems [dble]
  !! @param[out] force forces of target systems [dble]
  !! @param[out] virial of target systems [dble]
  !! @param[out] morpgh energy of target systems [dble]
  !! @date       2015/02/06 (CK)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_morph_sc(enefunc, coord, force, virial, emorph, cv)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: emorph
    real(wp),                intent(inout) :: cv(1:2)

    ! local variables
    integer                  :: i, j, k, id, omp_get_thread_num
    integer                  :: istart_sc, iend_sc
    integer                  :: i1, i2
    integer                  :: num_morph_sc
    real(wp)                 :: d12(1:3), r12, vtmp
    real(wp)                 :: t, tt(2), tmp, dn, coef, coef_t, tott(2)
    real(wp)                 :: vxx, vyx, vzx, vyy, vzy, vzz
    real(wp)                 :: work(3)
    integer,         pointer :: list_sc(:,:)
    real(wp),        pointer :: rmin_sc(:)
    real(wp),        pointer :: rmin_sc_other(:)


    call timer(TimerMorph, TimerOn)

    istart_sc = enefunc%istart_morph_sc
    iend_sc   = enefunc%iend_morph_sc

    list_sc   => enefunc%morph_list_sc
    rmin_sc   => enefunc%morph_dist_sc
    rmin_sc_other   => enefunc%morph_dist_sc_other

    num_morph_sc = enefunc%num_morph_sc

    emorph  = 0.0_wp
    t       = 0.0_wp
    tt(1:2) = 0.0_wp
    dn = real(num_morph_sc, wp)

    !$omp parallel do default(none)            &
    !$omp private(i,  d12, r12,  tmp, i1, i2)  &
    !$omp shared(istart_sc, iend_sc, coord, list_sc, rmin_sc, num_morph_sc,  &
    !$omp        rmin_sc_other)                &
    !$omp reduction(+:tt)
    !

    do i = istart_sc, iend_sc

      i1 = list_sc(1,i)
      i2 = list_sc(2,i)
      d12(1:3) = coord(1:3,i1) - coord(1:3,i2)
      r12   = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )
      tmp = r12 - rmin_sc(i)
      tt(1) = tt(1) + tmp*tmp
      tmp = r12 - rmin_sc_other(i)
      tt(2) = tt(2) + tmp*tmp

    end do
    !$omp end parallel do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(t, tott, 2, &
                    mpi_wp_real, mpi_sum, mpi_comm_country, ierror)
#else
    tott(1:2) = tt(1:2)
#endif
    cv(1:2) = sqrt(tott(1:2)/dn)

    t  = sqrt(tott(1)/dn)
    if (main_rank) &
      emorph = enefunc%morph_spring*t*t / (enefunc%morph_linear+t)

    coef_t = (1.0_wp/dn)*enefunc%morph_spring*(2.0d0*enefunc%morph_linear+t)/  &
              ((enefunc%morph_linear+t)*(enefunc%morph_linear+t))

    !$omp parallel default(none)                          &
    !$omp private(i, j, k, d12, r12, tmp, vtmp, coef, work, i1, i2, id)    &
    !$omp shared(istart_sc, iend_sc, coord, list_sc,coef_t, rmin_sc, dn, t,&
    !$omp        force, nthread)     &
    !$omp reduction(+:virial)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = istart_sc+id, iend_sc, nthread
      i1 = list_sc(1,i)
      i2 = list_sc(2,i)
      d12(1:3) = coord(1:3,i1) - coord(1:3,i2)
      r12   = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )
      tmp = r12 - rmin_sc(i)
      coef = coef_t*tmp/r12

      work(1:3) = coef * d12(1:3)

      ! virial
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = d12(k) * work(j) 
          virial(k,j) = virial(k,j) - vtmp
          virial(j,k) = virial(j,k) - vtmp
        end do
        vtmp = d12(j) * work(j) 
        virial(j,j) = virial(j,j) - vtmp
      end do
      force(1:3,i1,id+1) = force(1:3,i1,id+1) - work(1:3)
      force(1:3,i2,id+1) = force(1:3,i2,id+1) + work(1:3)
    end do

    !$omp end parallel

    call timer(TimerMorph, TimerOff)

    return

  end subroutine compute_energy_morph_sc

end module at_energy_morph_mod
