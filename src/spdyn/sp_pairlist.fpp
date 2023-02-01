!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_pairlist_mod
!> @brief   set pairlist for nonbonded interactions
!! @authors Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_pairlist_mod

  use sp_pairlist_gpu_mod
  use sp_pairlist_fugaku_mod
  use sp_pairlist_generic_mod
  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public  :: setup_pairlist
  public  :: update_pairlist_pbc
  public  :: update_pairlist_pbc_check

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_pairlist
  !> @brief        initialize/allocate/setup pairlist for nonbonded interactions
  !! @authors      JJ
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_pairlist(enefunc, domain, pairlist)

    ! formal arguments
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_domain),          intent(inout) :: domain
    type(s_pairlist),        intent(inout) :: pairlist


    call init_pairlist(pairlist)

    pairlist%pairlistdist = enefunc%pairlistdist

    if (enefunc%pairlist_check .or. (.not.enefunc%pme_use) .or. &
        enefunc%nonb_limiter) then

      call alloc_pairlist(pairlist, PairListGeneric, domain%num_cell_local)

    else

      select case (domain%pairlist_kernel)

      case (PLK_Generic)
        call alloc_pairlist(pairlist, PairListGeneric, domain%num_cell_local)
  
      case (PLK_Fugaku)
        call alloc_pairlist(pairlist, PairListFugaku, &
                            domain%num_cell_local+domain%num_cell_boundary)

      case (PLK_Intel)
        call alloc_pairlist(pairlist, PairListFugaku, &
                            domain%num_cell_local+domain%num_cell_boundary)

      case (PLK_GPU)
        call alloc_pairlist(pairlist, PairListGPU, domain%num_cell_local)

      end select

    end if


    if (domain%fep_use) then
      ! FEP
      call alloc_pairlist(pairlist, PairListFEP, domain%num_cell_local)
    end if

    call timer(TimerPairList, TimerOn)

    if (enefunc%pairlist_check) then
      call update_pairlist_pbc_check(enefunc, domain, pairlist)
    else
      call update_pairlist_pbc(enefunc, domain, pairlist)
    end if

    call timer(TimerPairList, TimerOff)

    return

  end subroutine setup_pairlist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_pbc
  !> @brief        update pairlist in each domain with periodic boundary
  !!               condition
  !! @authors      JJ
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_pbc(enefunc, domain, pairlist)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(inout) :: domain
    type(s_pairlist), target, intent(inout) :: pairlist

    if (domain%fep_use) then
      ! FEP
      if (.not. enefunc%pme_use .or. enefunc%nonb_limiter) then
        call update_pairlist_pbc_generic_fep(enefunc, domain, pairlist)

      else

        select case (domain%pairlist_kernel)

        case (PLK_Generic)
          call update_pairlist_pbc_generic_fep(enefunc, domain, pairlist)
    
        case (PLK_Fugaku)
          call update_pairlist_pbc_fugaku_fep(enefunc, domain, pairlist)

        case (PLK_Intel)
          call update_pairlist_pbc_fugaku_fep(enefunc, domain, pairlist)

        case (PLK_GPU)
          call update_pairlist_pbc_gpu(enefunc, domain, pairlist)

        end select

      end if

    else

      if (.not. enefunc%pme_use .or. enefunc%nonb_limiter) then
        call update_pairlist_pbc_generic(enefunc, domain, pairlist)

      else

        select case (domain%pairlist_kernel)

        case (PLK_Generic)
          call update_pairlist_pbc_generic(enefunc, domain, pairlist)
    
        case (PLK_Fugaku)
          call update_pairlist_pbc_fugaku(enefunc, domain, pairlist)

        case (PLK_Intel)
          call update_pairlist_pbc_fugaku(enefunc, domain, pairlist)

        case (PLK_GPU)
          call update_pairlist_pbc_gpu(enefunc, domain, pairlist)

        end select

      end if

    end if

    return

  end subroutine update_pairlist_pbc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_pbc_check
  !> @brief        update pairlist in each domain with periodic boundary
  !!               condition
  !! @authors      JJ, CK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_pbc_check(enefunc, domain, pairlist)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(in)    :: domain
    type(s_pairlist), target, intent(inout) :: pairlist


    if (domain%fep_use) then

      ! FEP
      call update_pairlist_pbc_check_generic_fep(enefunc, domain, pairlist)

    else

      call update_pairlist_pbc_check_generic(enefunc, domain, pairlist)

    end if

    return

  end subroutine update_pairlist_pbc_check

end module sp_pairlist_mod
