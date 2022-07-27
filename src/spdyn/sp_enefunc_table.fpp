!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_enefunc_table_mod
!> @brief   define lookup table for energy calculation
!! @authors Jaewoon Jung (JJ), Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_enefunc_table_mod

  use sp_boundary_mod
  use sp_energy_pme_mod
  use sp_energy_mod
  use sp_energy_str_mod
  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use molecules_str_mod
  use table_libs_mod
  use mpi_parallel_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! variables
  real(wp),         save      :: el_fact    ! e^2/(4*pi*eps) A*kcal/mol
  real(wp),         save      :: alpha      ! alpha
  real(wp),         save      :: alpha2m    ! -alpha^2
  real(wp),         save      :: alpha2sp   ! 2*alpha/sqrt(pi)

  public  :: setup_enefunc_table
  private :: setup_table_general_pme_linear
  private :: setup_table_general_cutoff_cubic

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_table
  !> @brief        define lookup table of potential energy function
  !! @authors      JJ
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_table(ene_info, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: cutoff, cutoff2
    real(wp)                 :: cutoff2_water, cutoff2_water1, density
    real(wp)                 :: mind
    integer                  :: cutoff_int, cutoff_int2, cutoff_int1


    if (.not. ene_info%table) &
      return

    ! set cutoff distance etc.
    !
    enefunc%table%table       = ene_info%table
    enefunc%table%table_order = ene_info%table_order
    enefunc%table%water_model = ene_info%water_model
    enefunc%table%density     = ene_info%table_density

    ! setting common parameters
    !
    el_fact  = ELECOEF/ene_info%dielec_const
    alpha    = ene_info%pme_alpha
    alpha2m  = -alpha**2
    alpha2sp = 2.0_wp*alpha/sqrt(PI)

    ! set parameters for table
    !
    cutoff   = ene_info%cutoffdist
    density  = ene_info%table_density

    cutoff2        = cutoff * cutoff
    cutoff2_water  = (cutoff+4.5_wp) * (cutoff+4.5_wp)
    cutoff2_water1 = (cutoff+1.5_wp) * (cutoff+1.5_wp)
    cutoff_int     = int(cutoff2*density)
    cutoff_int2    = int(cutoff2_water*density)
    cutoff_int1    = int(cutoff2_water1*density)

    if (enefunc%table%table_order == 1) then

      cutoff_int = int(2.1_wp*cutoff2_water1*density)
      cutoff_int2 = cutoff_int+1

    end if

    enefunc%table%cutoff_int = cutoff_int2

    if (enefunc%contact_check .or. enefunc%nonb_limiter) then
      if (enefunc%table%table_order /= 1) then
        if (main_rank) &
          write(MsgOut,'(A)') 'Setup_Enefunc_Table> Warning:'//&
          'concact_check is only in table_order =1'
      else
        mind=cutoff2*density/real(cutoff_int2,wp)+0.001_wp
        enefunc%minimum_contact=max(mind, enefunc%minimum_contact)
      endif
    endif

    call alloc_enefunc(enefunc, EneFuncTableDomain, cutoff_int2, 1)

    if (enefunc%pme_use) then

        call setup_table_general_pme_linear(ene_info, cutoff2, cutoff_int,  &
                                            cutoff_int2, enefunc)

    else

        call setup_table_general_cutoff_cubic(ene_info, cutoff2, cutoff_int,&
                                              cutoff_int2, enefunc)

    end if

    return

  end subroutine setup_enefunc_table

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_table_general_pme_linear
  !> @brief        define lookup table of potential energy function
  !!               (linear)
  !! @authors      JJ, CK
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    cutoff2     : square cutoff distance
  !! @param[in]    cutoff_int  : integer cutoff distance
  !! @param[in]    cutoff_int1 : square integer cutoff distance
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_table_general_pme_linear(ene_info, cutoff2, &
                                            cutoff_int, cutoff_int1, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    real(wp),                intent(in)    :: cutoff2
    integer,                 intent(in)    :: cutoff_int
    integer,                 intent(in)    :: cutoff_int1
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: switchdist, cutoff
    real(wp)                 :: switchdist2
    real(wp)                 :: density
    integer                  :: switch_int

    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: table_ecor(:), table_decor(:)
    real(wp),        pointer :: table_vcor(:), table_dvcor(:)


    table_ene   => enefunc%table%table_ene
    table_grad  => enefunc%table%table_grad
    table_ecor  => enefunc%table%table_ecor
    table_decor => enefunc%table%table_decor
    table_vcor  => enefunc%table%table_vcor
    table_dvcor => enefunc%table%table_dvcor

    density     = enefunc%table%density

    ! set cutoff and switch distance etc.
    !
    switchdist  = ene_info%switchdist
    cutoff      = ene_info%cutoffdist

    ! set lennard-jones parameters
    !
    switchdist2 = switchdist * switchdist
    if (switchdist2 > 0.0_wp) then
      switch_int  = int(1.0_wp/switchdist2*cutoff2*density)
    end if

    table_ene  (1:3*cutoff_int1) = 0.0_wp
    table_ecor (1:cutoff_int1)   = 0.0_wp
    table_decor(1:cutoff_int1)   = 0.0_wp
    table_grad (1:3*cutoff_int1) = 0.0_wp

    if (enefunc%vdw == VDWPME) then

      call table_pme_linear_pme(cutoff_int, density, cutoff2,              &
                                el_fact, alpha, alpha2sp, alpha2m,         &
                                table_ene, table_grad, table_ecor,         &
                                table_decor, table_vcor, table_dvcor)


    else if (enefunc%vdw == VDWFSW) then

      call table_pme_linear_fswitch(switch_int, cutoff_int, density,         &
                                   cutoff2, switchdist2, el_fact, alpha,     &
                                   alpha2sp, alpha2m,                        &
                                   table_ene, table_grad, table_ecor,        &
                                   table_decor)

    else if (enefunc%vdw == VDWCutoff) then

      enefunc%vdw_no_switch = .true.
      call table_pme_linear_noswitch(cutoff_int, density, cutoff2, el_fact, &
                                     alpha, alpha2sp, alpha2m, table_ene,   &
                                     table_grad, table_ecor, table_decor)

    else if (enefunc%forcefield == ForcefieldGROAMBER .or.                 &
             enefunc%forcefield == ForcefieldGROMARTINI) then

      if (enefunc%vdw == VDWShift) then

        call table_pme_linear_groshift(switch_int, cutoff_int, density,      &
                                   cutoff, switchdist, el_fact, alpha,       &
                                   alpha2sp, alpha2m,                        &
                                   table_ene, table_grad, table_ecor,        &
                                   table_decor)

      else

        call table_pme_linear_groswitch(switch_int, cutoff_int, density,     &
                                   cutoff, switchdist, el_fact, alpha,       &
                                   alpha2sp, alpha2m,                        &
                                   table_ene, table_grad, table_ecor,        &
                                   table_decor, enefunc%eswitch,             &
                                   enefunc%vswitch)

      endif

    else if (enefunc%vdw == VDWSwitch) then

      call table_pme_linear_switch(switch_int, cutoff_int, density,          &
                                   cutoff2, switchdist2, el_fact, alpha,     &
                                   alpha2sp, alpha2m,                        &
                                   table_ene, table_grad, table_ecor,        &
                                   table_decor)

    end if

    return

  end subroutine setup_table_general_pme_linear

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_table_general_cutoff_cubic
  !> @brief        define lookup table for cutoff
  !! @authors      JJ, CK
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    cutoff2     : square cutoff distance
  !! @param[in]    cutoff_int  : integer cutoff distance
  !! @param[in]    cutoff_int1 : square integer cutoff distance
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_table_general_cutoff_cubic(ene_info, cutoff2,  &
                                              cutoff_int, cutoff_int1, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    real(wp),                intent(in)    :: cutoff2
    integer,                 intent(in)    :: cutoff_int
    integer,                 intent(in)    :: cutoff_int1
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: switchdist
    real(wp)                 :: switchdist2, cutoff2_inv, cutoff
    real(wp)                 :: density
    integer                  :: switch_int

    real(wp),        pointer :: table_ene(:)
    real(wp),        pointer :: table_grad(:)


    density     = enefunc%table%density
    table_ene  => enefunc%table%table_ene
    table_grad => enefunc%table%table_grad

    cutoff2_inv = 1.0_wp/cutoff2
    el_fact = ELECOEF/enefunc%dielec_const 

    ! set switch distance etc.
    !
    switchdist  = ene_info%switchdist
    cutoff      = ene_info%cutoffdist

    ! set lennard-jones parameters
    !
    switchdist2 = switchdist * switchdist
    switch_int  = int(switchdist2*density)

    table_ene(1:6*cutoff_int1)  = 0.0_wp
    table_grad(1:6*cutoff_int1) = 0.0_wp

    if (enefunc%vdw == VDWFSW) then

      call table_cutoff_cubic_fswitch(switch_int, cutoff_int, cutoff_int1,  &
                                     density, el_fact,                      &
                                     cutoff2, switchdist2, table_ene,       &
                                     table_grad)


    else if (enefunc%vdw == VDWCutoff) then

      call table_cutoff_cubic_noswitch(switch_int, cutoff_int, cutoff_int1, &
                                       density, el_fact, cutoff2,           &
                                       table_ene, table_grad)

    else if (enefunc%forcefield == ForcefieldGROAMBER .or.                &
            enefunc%forcefield == ForcefieldGROMARTINI) then

      if (enefunc%vdw == VDWShift) then

        if (enefunc%forcefield == ForcefieldGROAMBER) then
          call table_cutoff_cubic_groshift(switch_int, cutoff_int,          &
                                           cutoff_int1, density, el_fact,   &
                                           cutoff, switchdist, table_ene,   &
                                           table_grad)
        else
          call table_cutoff_cubic_grodoubleshift(switch_int, cutoff_int,    &
                                         cutoff_int1, density, el_fact,     &
                                         cutoff, switchdist, table_ene,     &
                                         table_grad)
        endif

      else

        call table_cutoff_cubic_groswitch(switch_int, cutoff_int, cutoff_int1,&
                                         density, el_fact,                   &
                                         cutoff, switchdist, table_ene,      &
                                         table_grad)
      endif

    else if (enefunc%vdw == VDWSwitch) then 

      call table_cutoff_cubic_switch(switch_int, cutoff_int, cutoff_int1,    &
                                     density, el_fact,                       &
                                     cutoff2, switchdist2, table_ene,        &
                                     table_grad)

    end if

    return

  end subroutine setup_table_general_cutoff_cubic

end module sp_enefunc_table_mod
