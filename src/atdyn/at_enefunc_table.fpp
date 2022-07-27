!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_enefunc_table_mod
!> @brief   define lookup table for energy calculation
!! @authors Jaewoon Jung (JJ), Yuji Sugita (YS), Takaharu Mori (TM), 
!!          Chigusa Kobayashi (CK), Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_enefunc_table_mod

  use at_energy_mod
  use at_energy_str_mod
  use at_enefunc_str_mod
  use table_libs_mod
  use math_libs_mod
  use nbond_list_mod
  use molecules_str_mod
  use fileio_table_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! constants
  integer,          parameter :: MAX_EXCL = 52 
                               ! = 4(1-2) + 4x3 (1-3) + 3*4*3 (1-4:CH19) 
  integer,          parameter :: MAX_NB_14      = 36 ! = 3x4x3 (1-4) 
  real(wp),         parameter :: FACT_NUM_NB15  = 1.5_wp
  real(wp),         parameter :: FACT_THRESHOLD = 0.95_wp

  ! subroutines
  public  :: setup_enefunc_table
  private :: setup_table_general_pme_linear
  private :: setup_table_general_cutoff_cubic
  private :: setup_table_water_pme_linear
  private :: setup_user_table_general_pme_linear
  private :: setup_user_table_general_cutoff_cubic
  private :: setup_user_table_water_pme_linear
  private :: setup_solute_water_lists
  private :: count_nonb_excl_water
  private :: count_nonb_excl_c19
  private :: count_nonb_ecqm14

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_table
  !> @brief        define lookup table of potential energy function
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[in]    table    : lookup table information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !! @authors      JJ, TM
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_table(ene_info, table, molecule, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_table),           intent(in)    :: table
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: cutoff, cutoff2_water, cutoff2
    real(wp)                 :: el_fact, alpha, alpha2m, alpha2sp
    real(wp)                 :: mind
    integer                  :: cutoff_int
    integer                  :: cutoff_int2


    if (ene_info%table) then
      enefunc%table%table = .true.
    else
      enefunc%table%table = .false.

!CK20180325
      if (ene_info%forcefield == ForcefieldKBGO   .or. &
          ene_info%forcefield == ForcefieldAAGO   .or. &
          ene_info%forcefield == ForcefieldCAGO) then
        call setup_solute_water_lists(molecule, enefunc)
      end if
!CK20180325
      return
    end if

    ! set variables
    !
    enefunc%table%table_order = ene_info%table_order
    enefunc%table%water_model = ene_info%water_model
    enefunc%table%density     = ene_info%table_density

    ! set common parameters
    !
    cutoff        = ene_info%cutoffdist
    cutoff2       = cutoff * cutoff
    cutoff2_water = (cutoff + 4.5_wp)*(cutoff + 4.5_wp)        !TODO 4.5?

    if (enefunc%table%table_order == 3) then
      cutoff_int  = int(enefunc%table%density*cutoff2       )
      cutoff_int2 = int(enefunc%table%density*cutoff2_water )
    else
      cutoff_int  = int(enefunc%table%density*cutoff2*2.1_wp)  !TODO 2.1?
      cutoff_int2 = cutoff_int + 1
    end if

    if (enefunc%contact_check .or. enefunc%nonb_limiter) then
      if (enefunc%table%table_order .ne. 1) then
        if (main_rank) &
          write(MsgOut,'(A)') 'Setup_Enefunc_Table> WARNING: concact_check/nonb_limiter is only in table_order = 1'
      else
        mind=cutoff2*enefunc%table%density/real(cutoff_int2,wp)+0.001_wp
        enefunc%minimum_contact=max(mind, enefunc%minimum_contact)
      end if
    end if

    call alloc_enefunc(enefunc, EneFuncTable, cutoff_int2, 1)

    ! make water and solute lists
    !
    call setup_solute_water_lists(molecule, enefunc)

    ! define lookup table of potential energy function
    !
    if (ene_info%electrostatic == ElectrostaticPME) then

      if (.not. ene_info%user_def_table) then

        call setup_table_general_pme_linear(ene_info,            &
                                            cutoff2, cutoff_int, &
                                            cutoff_int2, enefunc)
        call setup_table_water_pme_linear  (ene_info,            &
                                            cutoff2, cutoff_int, &
                                            cutoff_int2, enefunc)
      else

        call setup_user_table_general_pme_linear( &
                                            ene_info, table,     &
                                            cutoff2, cutoff_int, &
                                            cutoff_int2, enefunc)
        call setup_user_table_water_pme_linear  ( &
                                            ene_info,            &
                                            cutoff2, cutoff_int, &
                                            cutoff_int2, enefunc)
      end if

    else

      if (.not. ene_info%user_def_table) then

        call setup_table_general_cutoff_cubic( &
                                            ene_info,            &
                                            cutoff2, cutoff_int, &
                                            cutoff_int2, enefunc)

      else

        call setup_user_table_general_cutoff_cubic( &
                                            ene_info, table,     &
                                            cutoff2, cutoff_int, &
                                            cutoff_int2, enefunc)

      end if

    end if

    ! exclude 1-2, 1-3 interactions for solute-solute interaction
    !
    if (enefunc%forcefield == ForcefieldCHARMM19) then
      call count_nonb_excl_c19(molecule, enefunc)
    else
      call count_nonb_excl_water(molecule, enefunc)
    end if
    if(enefunc%qmmm%do_qmmm .and. enefunc%qmmm%num_qmmmbonds > 0)  &
                                                call count_nonb_ecqm14(enefunc)

    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Setup_Enefunc_Table> Setup Variables for LookupTable'
      write(MsgOut,'(A20,I10,A20,I10)')                       &
            '  num_solutes     = ', enefunc%table%num_solute, &
            '  num_waters      = ', enefunc%table%num_water
      write(MsgOut,'(A)') ''
    end if

    return

  end subroutine setup_enefunc_table

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_table_general_pme_linear
  !> @brief        define lookup table of potential energy function (linear)
  !! @authors      JJ
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    table       : lookup table information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_table_general_pme_linear(ene_info, cutoff2, &
                                            cutoff_int, cutoff_int2, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    real(wp),                intent(in)    :: cutoff2
    integer,                 intent(in)    :: cutoff_int
    integer,                 intent(in)    :: cutoff_int2
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: switchdist, cutoff
    real(wp)                 :: switchdist2
    real(wp)                 :: density, el_fact
    real(wp)                 :: alpha, alpha2m, alpha2sp
    integer                  :: i, switch_int

    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: table_ecor(:), table_decor(:)


    el_fact     = ELECOEF / enefunc%dielec_const
    density     = enefunc%table%density
    alpha       = ene_info%pme_alpha
    alpha2m     = -alpha**2
    alpha2sp    = 2.0_wp*alpha/sqrt(PI)

    table_ene   => enefunc%table%table_ene
    table_ecor  => enefunc%table%table_ecor
    table_decor => enefunc%table%table_decor
    table_grad  => enefunc%table%table_grad

    ! set switch distance etc.
    !
    switchdist  = ene_info%switchdist
    cutoff      = ene_info%cutoffdist

    ! set lennard-jones parameters
    !
    switchdist2 = switchdist * switchdist
    if (switchdist2 > 0.0_wp) then
      switch_int  = int(1.0_wp/switchdist2*cutoff2*density)
    endif

    table_ene  (1:3*cutoff_int2) = 0.0_wp
    table_ecor (1:cutoff_int2)   = 0.0_wp
    table_decor(1:cutoff_int2)   = 0.0_wp
    table_grad (1:3*cutoff_int2) = 0.0_wp

    if (ene_info%vdw_force_switch) then

      call table_pme_linear_fswitch(switch_int, cutoff_int, density,         &
                                     cutoff2, switchdist2, el_fact, alpha,   &
                                     alpha2sp, alpha2m,                      &
                                     table_ene, table_grad,                  &
                                     table_ecor, table_decor)

    else if (enefunc%forcefield == ForcefieldAMBER) then

      call table_pme_linear_noswitch_atdyn(cutoff_int, density, cutoff2,     &
                                     el_fact, alpha, alpha2sp, alpha2m,      &
                                     table_ene, table_grad, table_ecor,      &
                                     table_decor)

    else if (enefunc%forcefield .eq. ForcefieldGROAMBER .or.                 &
            enefunc%forcefield .eq. ForcefieldGROMARTINI) then
      if (ene_info%vdw_shift) then
        call table_pme_linear_groshift(switch_int, cutoff_int, density,      &
                                   cutoff, switchdist, el_fact, alpha,       &
                                   alpha2sp, alpha2m,                        &
                                   table_ene, table_grad, table_ecor,        &
                                   table_decor)
      else if (abs(switchdist-cutoff) < EPS) then
        call table_pme_linear_noswitch_atdyn(cutoff_int, density, cutoff2,   &
                                       el_fact, alpha, alpha2sp, alpha2m,    &
                                       table_ene, table_grad, table_ecor,    &
                                       table_decor)
      else
        call table_pme_linear_groswitch(switch_int, cutoff_int, density,     &
                                   cutoff, switchdist, el_fact, alpha,       &
                                   alpha2sp, alpha2m,                        &
                                   table_ene, table_grad, table_ecor,        &
                                   table_decor, enefunc%eswitch,             &
                                   enefunc%vswitch)

      end if
 
    else

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
  !! @authors      JJ
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    table       : lookup table information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_table_general_cutoff_cubic(ene_info, cutoff2, &
                                              cutoff_int, cutoff_int2, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    real(wp),                intent(in)    :: cutoff2
    integer,                 intent(in)    :: cutoff_int
    integer,                 intent(in)    :: cutoff_int2
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: switchdist, cutoff
    real(wp)                 :: switchdist2, cutoff2_inv
    real(wp)                 :: density
    real(wp)                 :: el_fact
    integer                  :: i, switch_int

    real(wp),        pointer :: table_ene(:)
    real(wp),        pointer :: table_grad(:)


    density    =  enefunc%table%density
    table_ene  => enefunc%table%table_ene
    table_grad => enefunc%table%table_grad

    el_fact = ELECOEF/enefunc%dielec_const

    ! set switch distance etc.
    !
    switchdist  = ene_info%switchdist
    cutoff      = ene_info%cutoffdist

    ! set lennard-jones parameters
    !
    switchdist2 = switchdist * switchdist
    switch_int  = int(switchdist2*density)

    table_ene(1:6*cutoff_int2)  = 0.0_wp
    table_grad(1:6*cutoff_int2) = 0.0_wp


    if (enefunc%forcefield == ForcefieldAMBER .and. &
        abs(switchdist-cutoff) < EPS) then

      call table_cutoff_cubic_noswitch(switch_int, cutoff_int, cutoff_int2,  &
                                       density, el_fact, cutoff2, table_ene, &
                                       table_grad)

    else if (ene_info%vdw_force_switch) then

      call table_cutoff_cubic_fswitch(switch_int, cutoff_int, cutoff_int2, &
                                     density, el_fact,                     &
                                   cutoff2, switchdist2, table_ene, table_grad)

    else if (enefunc%forcefield == ForcefieldGROAMBER .or.                 &
             enefunc%forcefield == ForcefieldGROMARTINI) then

      if (ene_info%vdw_shift) then

        if (enefunc%forcefield == ForcefieldGROAMBER) then
          call table_cutoff_cubic_groshift(switch_int, cutoff_int,         &
                                         cutoff_int2, density, el_fact,    &
                                         cutoff, switchdist, table_ene,    &
                                         table_grad)
        else
          call table_cutoff_cubic_grodoubleshift(switch_int, cutoff_int,   &
                                         cutoff_int2, density, el_fact,    &
                                         cutoff, switchdist, table_ene,    &
                                         table_grad)
        end if

      else if (abs(switchdist-cutoff) < EPS) then 

        call table_cutoff_cubic_noswitch(switch_int, cutoff_int, cutoff_int2,  &
                                         density, el_fact, cutoff2,            &
                                         table_ene, table_grad)
      else

        call table_cutoff_cubic_groswitch(switch_int, cutoff_int, cutoff_int2, &
                                      density, el_fact,                        &
                                      cutoff, switchdist, table_ene,           &
                                      table_grad)
      end if
     
    else

      call table_cutoff_cubic_switch(switch_int, cutoff_int, cutoff_int2,      &
                                     density, el_fact,                         &
                                     cutoff2, switchdist2, table_ene,          &
                                     table_grad)

    end if

    return

  end subroutine setup_table_general_cutoff_cubic

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_table_water_pme_linear
  !> @brief        define lookup table of water (linear interpolation)
  !! @authors      JJ
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    table       : lookup table information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_table_water_pme_linear(ene_info, cutoff2, &
                                          cutoff_int, cutoff_int2, enefunc)

    ! formal arguments
    type(s_ene_info),         intent(in)    :: ene_info
    real(wp),                 intent(in)    :: cutoff2
    integer,                  intent(in)    :: cutoff_int
    integer,                  intent(in)    :: cutoff_int2
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variables
    real(wp)                  :: lj6_OO, lj6_OH, lj6_HH
    real(wp)                  :: lj12_OO, lj12_OH, lj12_HH
    real(wp)                  :: cc_OO, cc_OH, cc_HH, charge_O, charge_H
    real(wp)                  :: switchdist, cutoff
    real(wp)                  :: switchdist2
    real(wp)                  :: density, el_fact
    real(wp)                  :: alpha, alpha2m, alpha2sp
    integer                   :: i, switch_int, atmcls_O, atmcls_H

    real(wp),         pointer :: table_de_WW(:,:)
    real(wp),         pointer :: table_ene_WW(:,:)


    if (enefunc%table%num_water == 0) return

    el_fact     = ELECOEF / enefunc%dielec_const
    density     = enefunc%table%density
    alpha       = ene_info%pme_alpha
    alpha2m     = -alpha**2
    alpha2sp    = 2.0_wp*alpha/sqrt(PI)
    table_ene_WW => enefunc%table%table_ene_WW
    table_de_WW  => enefunc%table%table_de_WW

    ! set switch distance etc.
    !
    switchdist  = ene_info%switchdist
    cutoff      = ene_info%cutoffdist

    ! set lennard-jones parameters
    !
    switchdist2 = switchdist * switchdist
    if (switchdist2 > 0.0_wp) then
      switch_int  = int(1.0_wp/switchdist2*cutoff2*density)
    endif

    table_ene_WW(1:6*cutoff_int2,1:3) = 0.0_wp
    table_de_WW(1:2*cutoff_int2,1:3)  = 0.0_wp

    ! atom class and charge of water oxygen and hydrogen
    !
    atmcls_O = enefunc%table%atom_cls_no_O
    atmcls_H = enefunc%table%atom_cls_no_H
    charge_O = enefunc%table%charge_O
    charge_H = enefunc%table%charge_H

    lj6_OO = enefunc%nonb_lj6(atmcls_O,atmcls_O)
    lj6_OH = enefunc%nonb_lj6(atmcls_O,atmcls_H)
    lj6_HH = enefunc%nonb_lj6(atmcls_H,atmcls_H)

    lj12_OO = enefunc%nonb_lj12(atmcls_O,atmcls_O)
    lj12_OH = enefunc%nonb_lj12(atmcls_O,atmcls_H)
    lj12_HH = enefunc%nonb_lj12(atmcls_H,atmcls_H)

    cc_OO = charge_O * charge_O * el_fact
    cc_OH = charge_O * charge_H * el_fact
    cc_HH = charge_H * charge_H * el_fact

    if (ene_info%vdw_force_switch) then

      call table_water_pme_linear_fswitch(switch_int, cutoff_int, density,     &
                                          cutoff2, switchdist2, alpha,         &
                                          lj6_OO,  lj6_OH,  lj6_HH,            &
                                          lj12_OO, lj12_OH, lj12_HH,           &
                                          cc_OO, cc_OH, cc_HH,                 &
                                          alpha2sp, alpha2m,                   &
                                          table_ene_WW, table_de_WW)

    else if (enefunc%forcefield == ForcefieldAMBER) then

      call table_water_pme_linear_noswitch(cutoff_int, density, cutoff2,       &
                                           alpha, lj6_OO,  lj6_OH,  lj6_HH,    &
                                           lj12_OO, lj12_OH, lj12_HH,          &
                                           cc_OO, cc_OH, cc_HH,                &
                                           alpha2sp, alpha2m,                  &
                                           table_ene_WW, table_de_WW)

    else if (enefunc%forcefield == ForcefieldGROAMBER .or.                     &
             enefunc%forcefield == ForcefieldGROMARTINI) then
      if (ene_info%vdw_shift) then
        call table_water_pme_linear_groshift(switch_int, cutoff_int, density,  &
                                             cutoff, switchdist, alpha,        &
                                             lj6_OO,  lj6_OH,  lj6_HH,         &
                                             lj12_OO, lj12_OH, lj12_HH,        &
                                             cc_OO, cc_OH, cc_HH,              &
                                             alpha2sp, alpha2m,                &
                                             table_ene_WW, table_de_WW)
     else if (abs(switchdist-cutoff) < EPS) then
      call table_water_pme_linear_noswitch(cutoff_int, density, cutoff2,       &
                                           alpha, lj6_OO,  lj6_OH,  lj6_HH,    &
                                           lj12_OO, lj12_OH, lj12_HH,          &
                                           cc_OO, cc_OH, cc_HH,                &
                                           alpha2sp, alpha2m,                  &
                                           table_ene_WW, table_de_WW)
     else
        call table_water_pme_linear_groswitch(switch_int, cutoff_int, density, &
                                              cutoff, switchdist, alpha,       &
                                              lj6_OO,  lj6_OH,  lj6_HH,        &
                                              lj12_OO, lj12_OH, lj12_HH,       &
                                              cc_OO, cc_OH, cc_HH,             &
                                              alpha2sp, alpha2m,               &
                                              table_ene_WW, table_de_WW)

     end if

    else

      call table_water_pme_linear_switch(switch_int, cutoff_int, density,      &
                                         cutoff2, switchdist2, alpha,          &
                                         lj6_OO,  lj6_OH,  lj6_HH,             &
                                         lj12_OO, lj12_OH, lj12_HH,            &
                                         cc_OO, cc_OH, cc_HH,                  &
                                         alpha2sp, alpha2m,                    &
                                         table_ene_WW, table_de_WW)
    end if

    return

  end subroutine setup_table_water_pme_linear

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_user_table_general_pme_linear
  !> @brief        define lookup table of potential energy function (linear)
  !! @authors      JJ, NT
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    table       : lookup table information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_user_table_general_pme_linear(ene_info, table,     &
                                                 cutoff2, cutoff_int, &
                                                 cutoff_int2, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_table),           intent(in)    :: table
    real(wp),                intent(in)    :: cutoff2
    integer,                 intent(in)    :: cutoff_int
    integer,                 intent(in)    :: cutoff_int2
    type(s_enefunc),         intent(inout) :: enefunc


    if (.not. check_table(table, &
                     ene_info%cutoffdist,    &
                     ene_info%switchdist,    &
                     ene_info%dielec_const,  &
                     ene_info%table_density, &
                     ene_info%table_order,   &
                     ene_info%pme_alpha)) then
      call error_msg('Setup_Table_General_Pme_Linear> '//&
             'table data has different value from GENESIS settings')
    end if

    call get_table(table, &
                   enefunc%table%table_ene,      &
                   enefunc%table%table_grad,     &
                   enefunc%table%table_ecor,     &
                   enefunc%table%table_decor)
    write(MsgOut,'(A)') 'Setup_Table_General_Pme_Linear> '//&
        'Lookup table was built from user defined table file.'
    write(MsgOut,'(A)') ''

    return

  end subroutine setup_user_table_general_pme_linear

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_user_table_general_cutoff_cubic
  !> @brief        define lookup table for cutoff
  !! @authors      JJ
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    table       : lookup table information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_user_table_general_cutoff_cubic(ene_info, table,     &
                                                   cutoff2, cutoff_int, &
                                                   cutoff_int2, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_table),           intent(in)    :: table
    real(wp),                intent(in)    :: cutoff2
    integer,                 intent(in)    :: cutoff_int
    integer,                 intent(in)    :: cutoff_int2
    type(s_enefunc), target, intent(inout) :: enefunc


    if (.not. check_table(table, &
                   ene_info%cutoffdist,    &
                   ene_info%switchdist,    &
                   ene_info%dielec_const,  &
                   ene_info%table_density, &
                   ene_info%table_order,   &
                   ene_info%pme_alpha)) then
      call error_msg('Setup_Table_General_Pme_Linear> '//&
           'table data has different value from GENESIS settings')
    end if

    call get_table(table, &
                   enefunc%table%table_ene,      &
                   enefunc%table%table_grad)
    write(MsgOut,'(A)') 'Setup_Table_General_Cutoff_Cubic> '//&
        'Lookup table was built from user defined table file.'
    write(MsgOut,'(A)') ''

    return

  end subroutine setup_user_table_general_cutoff_cubic

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_user_table_water_pme_linear
  !> @brief        define lookup table of water (linear interpolation)
  !! @authors      JJ, CK
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    table       : lookup table information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_user_table_water_pme_linear(ene_info, cutoff2, &
                                               cutoff_int, cutoff_int2, enefunc)

    ! formal arguments
    type(s_ene_info),         intent(in)    :: ene_info
    real(wp),                 intent(in)    :: cutoff2
    integer,                  intent(in)    :: cutoff_int
    integer,                  intent(in)    :: cutoff_int2
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variables
    real(wp)                  :: lj6_OO, lj6_OH, lj6_HH
    real(wp)                  :: lj12_OO, lj12_OH, lj12_HH
    real(wp)                  :: cc_OO, cc_OH, cc_HH, charge_O, charge_H
    real(wp)                  :: switchdist, cutoff
    real(wp)                  :: switchdist2
    real(wp)                  :: density, el_fact
    real(wp)                  :: alpha, alpha2m, alpha2sp
    integer                   :: i, switch_int, atmcls_O, atmcls_H

    real(wp),         pointer :: table_de_WW(:,:)
    real(wp),         pointer :: table_ene_WW(:,:)


    if (enefunc%table%num_water == 0) return

    density     = enefunc%table%density
    table_ene_WW => enefunc%table%table_ene_WW
    table_de_WW  => enefunc%table%table_de_WW

    ! set switch distance etc.
    !
    switchdist  = ene_info%switchdist
    cutoff      = ene_info%cutoffdist

    ! set lennard-jones parameters
    !
    switchdist2 = switchdist * switchdist
    switch_int = int(1.0_wp/switchdist2*cutoff2*density)

    table_ene_WW(1:6*cutoff_int2,1:3) = 0.0_wp
    table_de_WW (1:2*cutoff_int2,1:3) = 0.0_wp

    ! atom class and charge of water oxygen and hydrogen
    !
    atmcls_O = enefunc%table%atom_cls_no_O
    atmcls_H = enefunc%table%atom_cls_no_H
    charge_O = enefunc%table%charge_O
    charge_H = enefunc%table%charge_H

    lj6_OO = enefunc%nonb_lj6(atmcls_O,atmcls_O)
    lj6_OH = enefunc%nonb_lj6(atmcls_O,atmcls_H)
    lj6_HH = enefunc%nonb_lj6(atmcls_H,atmcls_H)

    lj12_OO = enefunc%nonb_lj12(atmcls_O,atmcls_O)
    lj12_OH = enefunc%nonb_lj12(atmcls_O,atmcls_H)
    lj12_HH = enefunc%nonb_lj12(atmcls_H,atmcls_H)

    cc_OO = charge_O * charge_O 
    cc_OH = charge_O * charge_H 
    cc_HH = charge_H * charge_H 

      call table_water_pme_linear_user_defined(enefunc%table%table_ene,       &
                                               enefunc%table%table_grad,      &
                                               lj6_OO,  lj6_OH,  lj6_HH,      &
                                               lj12_OO, lj12_OH, lj12_HH,     &
                                               cc_OO, cc_OH, cc_HH,           &
                                               table_ene_WW, table_de_WW)  

    return

  end subroutine setup_user_table_water_pme_linear

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_solute_water_lists
  !> @brief        setup solute and water lists
  !! @authors      CK, TM, JJ
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_solute_water_lists(molecule, enefunc)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, natom, nwater, nsolute
    integer                  :: io, ih
    character(10)            :: water_model
    character(4)             :: an1, an2, an3


    water_model = enefunc%table%water_model

    natom   = molecule%num_atoms
    nwater  = 0
    nsolute = 0

    do i = 1, natom
      if (molecule%residue_name(i) .eq. water_model) then
        nwater  = nwater + 1
      else
        nsolute = nsolute + 1
      end if
    end do

    if ((natom-nsolute)/4 == nwater .and. nwater > 0) &
      call error_msg('Setup_Solute_Water_Lists> TIP4P is not avaialble in ATDYN.')
    if (mod(nwater, 3) /= 0) &
      call error_msg('Setup_Solute_Water_Lists> # of water is incorrect.')

    nwater = nwater / 3
    enefunc%table%num_water  = nwater
    enefunc%table%num_solute = nsolute

    call alloc_enefunc(enefunc, EneFuncTableWater,  nwater,  1)
    call alloc_enefunc(enefunc, EneFuncTableSolute, nsolute, natom)

    i = 1
    nwater  = 0
    nsolute = 0
    do while(.true.)

      if (i > natom) &
        exit

      if (molecule%residue_name(i) .eq. water_model) then

        an1 = molecule%atom_name(i)
        an2 = molecule%atom_name(i+1)
        an3 = molecule%atom_name(i+2)

        if (an1(1:1) .eq. 'H' .and. &
            an2(1:1) .eq. 'O' .and. &
            an3(1:1) .eq. 'H' .or.  &
            an1(1:1) .eq. 'h' .and. &
            an2(1:1) .eq. 'o' .and. &
            an3(1:1) .eq. 'h') then

          nwater = nwater + 1
          enefunc%table%water_list(1,nwater) = i+1
          enefunc%table%water_list(2,nwater) = i
          enefunc%table%water_list(3,nwater) = i+2

        else if ( &
            an1(1:1) .eq. 'O' .and. &
            an2(1:1) .eq. 'H' .and. &
            an3(1:1) .eq. 'H' .or.  &
            an1(1:1) .eq. 'o' .and. &
            an2(1:1) .eq. 'h' .and. &
            an3(1:1) .eq. 'h') then

          nwater = nwater + 1
          enefunc%table%water_list(1,nwater) = i
          enefunc%table%water_list(2,nwater) = i+1
          enefunc%table%water_list(3,nwater) = i+2

        else if ( &
            an1(1:1) .eq. 'H' .and. &
            an2(1:1) .eq. 'H' .and. &
            an3(1:1) .eq. 'O' .or.  &
            an1(1:1) .eq. 'h' .and. &
            an2(1:1) .eq. 'h' .and. &
            an3(1:1) .eq. 'o') then

          nwater = nwater + 1
          enefunc%table%water_list(1,nwater) = i+2
          enefunc%table%water_list(2,nwater) = i
          enefunc%table%water_list(3,nwater) = i+1

        end if

        i = i + 3

      else

        nsolute = nsolute + 1
        enefunc%table%solute_list(nsolute) = i
        enefunc%table%solute_list_inv(i) = nsolute

        i = i + 1

      end if

    end do

    if (nwater /= enefunc%table%num_water) &
      call error_msg('Setup_Solute_Water_Lists> number of water is incorrect')

    if (nsolute /= enefunc%table%num_solute) &
      call error_msg('Setup_Solute_Water_Lists> number of solute is incorrect')

    if (nwater > 0) then
      io = enefunc%table%water_list(1,1)
      ih = enefunc%table%water_list(2,1)
     
      enefunc%table%atom_cls_no_O = molecule%atom_cls_no(io)
      enefunc%table%atom_cls_no_H = molecule%atom_cls_no(ih)
      enefunc%table%charge_O      = molecule%charge(io)
      enefunc%table%charge_H      = molecule%charge(ih)
      enefunc%table%mass_O        = molecule%mass(io)
      enefunc%table%mass_H        = molecule%mass(ih)
    end if

    return

  end subroutine setup_solute_water_lists

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    count_nonb_excl_water
  !> @brief        exclude 1-2, 1-3 interactions for solute-solute interaction
  !! @authors      JJ
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine count_nonb_excl_water(molecule, enefunc)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: i, k, k1, k2, k3, k4, ii, nsolute
    integer                  :: num_excl, num_excl_total, num_excl_total1
    integer                  :: num_nb14, num_nb14_total
    integer                  :: max_nonb_excl_num, max_nb14_num
    integer                  :: alloc_stat, dealloc_stat
    logical                  :: duplicate
    integer                  :: ndihedral

    integer,         pointer :: solute(:)
    integer,         pointer :: dihelist(:,:)

    integer,     allocatable :: nonb_excl_list(:,:), nonb_excl_list1(:,:)
    integer,     allocatable :: nb14_calc_list(:,:)
    real(wp),    allocatable :: nb14_qq_scale (:,:)
    real(wp),    allocatable :: nb14_lj_scale (:,:)

    integer                  :: l, l1, l2, mm_atom
    integer                  :: num_bonds, num_angles
    type(s_qmmm),    pointer :: qmmm
    integer                  :: qm_natoms, ec_natoms
    integer,         pointer :: qmatom_id(:), ecatom_id(:)
    integer                  :: MAX_EXCL_LOC


    ! Return if Go model under PBC condition
    if ( enefunc%forcefield == ForcefieldKBGO .or. &
         enefunc%forcefield == ForcefieldAAGO .or. &
         enefunc%forcefield == ForcefieldCAGO .or. &
         enefunc%forcefield == ForcefieldSOFT ) then
      return
    end if

    solute  => enefunc%table%solute_list_inv
    nsolute = enefunc%table%num_solute

    qmmm      => enefunc%qmmm
    qm_natoms =  enefunc%qmmm%qm_natoms
    qmatom_id => enefunc%qmmm%qmatom_id
    ec_natoms =  enefunc%qmmm%ec_natoms
    ecatom_id => enefunc%qmmm%ecatom_id

    MAX_EXCL_LOC = MAX_EXCL
    if(qmmm%do_qmmm) MAX_EXCL_LOC = (qm_natoms - 1) + MAX_EXCL

    ! allocate nonbonded exclusion list (solute only and all) and 
    ! 1-4 interaction list
    !
    allocate(enefunc%table%num_nonb_excl(nsolute), stat=alloc_stat)
    if (alloc_stat /=0) call error_msg_alloc
!ky    allocate(nonb_excl_list(MAX_EXCL,nsolute),stat=alloc_stat)
    allocate(nonb_excl_list(MAX_EXCL_LOC,nsolute),stat=alloc_stat)
    if (alloc_stat /=0) call error_msg_alloc

    allocate(enefunc%num_nonb_excl(molecule%num_atoms), stat=alloc_stat)
    if (alloc_stat /=0) call error_msg_alloc
!ky    allocate(nonb_excl_list1(MAX_EXCL,molecule%num_atoms),stat=alloc_stat)
    allocate(nonb_excl_list1(MAX_EXCL_LOC,molecule%num_atoms),stat=alloc_stat)
    if (alloc_stat /=0) call error_msg_alloc

    allocate(enefunc%num_nb14_calc(nsolute), stat=alloc_stat)
    if (alloc_stat /=0) call error_msg_alloc
    allocate(nb14_calc_list(MAX_NB_14,nsolute), &
             nb14_qq_scale (MAX_NB_14,nsolute), &
             nb14_lj_scale (MAX_NB_14,nsolute),stat=alloc_stat)
    if (alloc_stat /=0) call error_msg_alloc

    enefunc%table%num_nonb_excl(1:nsolute) = 0
    enefunc%num_nb14_calc(1:nsolute) = 0
    enefunc%num_nonb_excl(1:molecule%num_atoms) = 0

    ! exclude 1-2 interaction
    !
    num_excl_total  = 0
    num_excl_total1 = 0

    if (enefunc%excl_level > 0) then
      if(qmmm%do_qmmm) then
        num_bonds = molecule%num_bonds - qmmm%qm_nbonds
      else
        num_bonds = molecule%num_bonds
      end if

      do k = 1, num_bonds
     
        k1 = solute(molecule%bond_list(1,k))
        k2 = solute(molecule%bond_list(2,k))
     
        if (k1 /= 0 .and. k2 /= 0) then
     
          if (k1 < k2) then
            num_excl = enefunc%table%num_nonb_excl(k1) + 1
            enefunc%table%num_nonb_excl(k1) = num_excl
            nonb_excl_list(num_excl,k1) = k2
            num_excl_total = num_excl_total + 1
          else
            num_excl = enefunc%table%num_nonb_excl(k2) + 1
            enefunc%table%num_nonb_excl(k2) = num_excl
            nonb_excl_list(num_excl,k2) = k1
            num_excl_total = num_excl_total + 1
          end if
     
        end if
     
        k1 = molecule%bond_list(1,k)
        k2 = molecule%bond_list(2,k)

        if (k1 < k2) then
          num_excl = enefunc%num_nonb_excl(k1) + 1
          enefunc%num_nonb_excl(k1) = num_excl
          nonb_excl_list1(num_excl,k1) = k2
          num_excl_total1 = num_excl_total1 + 1
        else
          num_excl = enefunc%num_nonb_excl(k2) + 1
          enefunc%num_nonb_excl(k2) = num_excl
          nonb_excl_list1(num_excl,k2) = k1
          num_excl_total1 = num_excl_total1 + 1
        end if
     
      end do

      if(qmmm%do_qmmm) then

        do l1 = 1, qm_natoms
          do l2 = 1, l1-1

            k1 = solute(qmatom_id(l1))
            k2 = solute(qmatom_id(l2))

            if (k1 /= 0 .and. k2 /= 0) then

              if (k1 < k2) then
                num_excl = enefunc%table%num_nonb_excl(k1) + 1
                enefunc%table%num_nonb_excl(k1) = num_excl
                nonb_excl_list(num_excl,k1) = k2
                num_excl_total = num_excl_total + 1
              else
                num_excl = enefunc%table%num_nonb_excl(k2) + 1
                enefunc%table%num_nonb_excl(k2) = num_excl
                nonb_excl_list(num_excl,k2) = k1
                num_excl_total = num_excl_total + 1
              end if

            end if
 
            k1 = qmatom_id(l1)
            k2 = qmatom_id(l2)

            if (k1 < k2) then
              num_excl = enefunc%num_nonb_excl(k1) + 1
              enefunc%num_nonb_excl(k1) = num_excl
              nonb_excl_list1(num_excl,k1) = k2
              num_excl_total1 = num_excl_total1 + 1
            else
              num_excl = enefunc%num_nonb_excl(k2) + 1
              enefunc%num_nonb_excl(k2) = num_excl
              nonb_excl_list1(num_excl,k2) = k1
              num_excl_total1 = num_excl_total1 + 1
            end if
     
          end do
        end do

      end if

    end if

    ! exclude 1-3 interaction
    !
    if (enefunc%excl_level > 1) then
      if(qmmm%do_qmmm) then
        num_angles = molecule%num_angles - qmmm%qm_nangles
      else
        num_angles = molecule%num_angles
      end if

      do k = 1, num_angles
     
        k1 = solute(molecule%angl_list(1,k))
        k3 = solute(molecule%angl_list(3,k))
     
        if (k1 /= 0 .and. k3 /= 0) then 
     
          if (k1 < k3) then
     
            num_excl = enefunc%table%num_nonb_excl(k1)
            duplicate = .false.
     
            do i = 1, enefunc%table%num_nonb_excl(k1)
              if (k3 == nonb_excl_list(i,k1)) duplicate = .true.
            end do
     
            if (.not.duplicate) then
              num_excl = num_excl + 1
              enefunc%table%num_nonb_excl(k1) = num_excl
              nonb_excl_list(num_excl,k1) = k3
              num_excl_total = num_excl_total + 1
            end if
     
          else
     
            num_excl = enefunc%table%num_nonb_excl(k3)
            duplicate = .false.
     
            do i = 1, enefunc%table%num_nonb_excl(k3)
              if (k1 == nonb_excl_list(i,k3)) duplicate = .true.
            end do
     
            if (.not.duplicate) then
              num_excl = num_excl + 1
              enefunc%table%num_nonb_excl(k3) = num_excl
              nonb_excl_list(num_excl,k3) = k1
              num_excl_total = num_excl_total + 1
            end if
     
          end if
        end if
     
        k1 = molecule%angl_list(1,k)
        k3 = molecule%angl_list(3,k)
     
        if (k1 < k3) then
     
          num_excl = enefunc%num_nonb_excl(k1)
          duplicate = .false.
     
          do i = 1, enefunc%num_nonb_excl(k1)
            if (k3 == nonb_excl_list1(i,k1)) duplicate = .true.
          end do
     
          if (.not.duplicate) then
             num_excl = num_excl + 1
             enefunc%num_nonb_excl(k1) = num_excl
             nonb_excl_list1(num_excl,k1) = k3
             num_excl_total1 = num_excl_total1 + 1
          end if
     
        else
     
          num_excl = enefunc%num_nonb_excl(k3)
          duplicate = .false.
     
          do i = 1, enefunc%num_nonb_excl(k3)
            if (k1 == nonb_excl_list1(i,k3)) duplicate = .true.
          end do
     
          if (.not.duplicate) then
            num_excl = num_excl + 1
            enefunc%num_nonb_excl(k3) = num_excl
            nonb_excl_list1(num_excl,k3) = k1
            num_excl_total1 = num_excl_total1 + 1
          end if
        end if
     
      end do
    end if

    ! count 1-4 interaction
    !
    num_nb14_total = 0
    if (enefunc%excl_level > 2) then
      do ii = 1, 2
        ndihedral = enefunc%num_dihedrals 
        dihelist => enefunc%dihe_list
        if (ii == 2) then
          ndihedral = enefunc%num_rb_dihedrals
          dihelist  => enefunc%rb_dihe_list
        end if
        do k = 1, ndihedral
     
          k1 = solute(dihelist(1,k))
          k4 = solute(dihelist(4,k))
     
          if (k1 /= 0 .and. k4 /= 0) then
     
            if (k1 < k4) then
     
              num_nb14 = enefunc%num_nb14_calc(k1)
              duplicate = .false.
     
              do i = 1, enefunc%table%num_nonb_excl(k1)
                if (k4 == nonb_excl_list(i,k1)) duplicate = .true.
              end do
     
              do i = 1, num_nb14
                if (k4 == nb14_calc_list(i,k1)) duplicate = .true.
              end do
     
              if (.not. duplicate) then
                num_nb14 = num_nb14 + 1
                enefunc%num_nb14_calc(k1) = num_nb14
                nb14_calc_list(num_nb14,k1) = k4
                if (enefunc%forcefield == ForcefieldAMBER) then
                  nb14_qq_scale (num_nb14,k1) = enefunc%dihe_scee(k)
                  nb14_lj_scale (num_nb14,k1) = enefunc%dihe_scnb(k)
                else if (enefunc%forcefield == ForcefieldGROAMBER .or.      &
                         enefunc%forcefield == ForcefieldGROMARTINI) then
                  nb14_qq_scale (num_nb14,k1) = enefunc%fudge_qq
                  nb14_lj_scale (num_nb14,k1) = enefunc%fudge_lj
                end if
                num_nb14_total = num_nb14_total + 1
              end if
     
            else
     
              num_nb14 = enefunc%num_nb14_calc(k4)
              duplicate = .false.
     
              do i = 1, enefunc%table%num_nonb_excl(k4)
                if (k1 == nonb_excl_list(i,k4)) duplicate = .true.
              end do
     
              do i = 1, num_nb14
                if (k1 == nb14_calc_list(i,k4)) duplicate = .true.
              end do
     
              if (.not.duplicate) then
                num_nb14 = num_nb14 + 1
                enefunc%num_nb14_calc(k4) = num_nb14
                nb14_calc_list(num_nb14,k4) = k1
                if (enefunc%forcefield == ForcefieldAMBER) then
                  nb14_qq_scale (num_nb14,k4) = enefunc%dihe_scee(k)
                  nb14_lj_scale (num_nb14,k4) = enefunc%dihe_scnb(k)
                else if (enefunc%forcefield == ForcefieldGROAMBER .or.  &
                         enefunc%forcefield == ForcefieldGROMARTINI) then
                  nb14_qq_scale (num_nb14,k4) = enefunc%fudge_qq
                  nb14_lj_scale (num_nb14,k4) = enefunc%fudge_lj
                end if
                num_nb14_total = num_nb14_total + 1
              end if
     
            end if
          end if
        end do
      end do

      if(qmmm%num_qmmmbonds > 0) then
        do k = molecule%num_dihedrals - qmmm%qm_ndihedrals + 1,  &
               molecule%num_dihedrals
          do l = 1, qmmm%num_qmmmbonds

            mm_atom = qmmm%qmmmbond_list(2,l)
            if(molecule%dihe_list(1,k) == mm_atom .or. &
               molecule%dihe_list(4,k) == mm_atom) then

              k1 = solute(molecule%dihe_list(1,k))
              k4 = solute(molecule%dihe_list(4,k))

              if (k1 /= 0 .and. k4 /= 0) then
                if (k1 < k4) then
         
                  num_nb14 = enefunc%num_nb14_calc(k1)
                  duplicate = .false.
         
                  do i = 1, enefunc%table%num_nonb_excl(k1)
                    if (k4 == nonb_excl_list(i,k1)) duplicate = .true.
                  end do
         
                  do i = 1, num_nb14
                    if (k4 == nb14_calc_list(i,k1)) duplicate = .true.
                  end do
         
                  if (.not.duplicate) then
                    num_nb14 = num_nb14 + 1
                    enefunc%num_nb14_calc(k1) = num_nb14
                    nb14_calc_list(num_nb14,k1) = k4
                    if (enefunc%forcefield == ForcefieldAMBER) then
                      nb14_qq_scale (num_nb14,k1) = enefunc%dihe_scee(k)
                      nb14_lj_scale (num_nb14,k1) = enefunc%dihe_scnb(k)
                    else if (enefunc%forcefield == ForcefieldGROAMBER .or. &
                             enefunc%forcefield == ForcefieldGROMARTINI) then
                      nb14_qq_scale (num_nb14,k1) = enefunc%fudge_qq
                      nb14_lj_scale (num_nb14,k1) = enefunc%fudge_lj
                    end if
                    num_nb14_total = num_nb14_total + 1
                  end if
         
                else
         
                  num_nb14 = enefunc%num_nb14_calc(k4)
                  duplicate = .false.
         
                  do i = 1, enefunc%table%num_nonb_excl(k4)
                    if (k1 == nonb_excl_list(i,k4)) duplicate = .true.
                  end do
         
                  do i = 1, num_nb14
                    if (k1 == nb14_calc_list(i,k4)) duplicate = .true.
                  end do
         
                  if (.not.duplicate) then
                    num_nb14 = num_nb14 + 1
                    enefunc%num_nb14_calc(k4) = num_nb14
                    nb14_calc_list(num_nb14,k4) = k1
                    if (enefunc%forcefield == ForcefieldAMBER) then
                      nb14_qq_scale (num_nb14,k4) = enefunc%dihe_scee(k)
                      nb14_lj_scale (num_nb14,k4) = enefunc%dihe_scnb(k)
                    else if (enefunc%forcefield == ForcefieldGROAMBER .or.  &
                             enefunc%forcefield == ForcefieldGROMARTINI) then
                      nb14_qq_scale (num_nb14,k4) = enefunc%fudge_qq
                      nb14_lj_scale (num_nb14,k4) = enefunc%fudge_lj
                    end if
                    num_nb14_total = num_nb14_total + 1
                  end if
         
                end if
              end if

            end if

          end do
        end do

      end if

    end if

    ! pack 2D-array into 1D-array
    !
    allocate(enefunc%nonb_excl_list(num_excl_total1), stat = alloc_stat)
    if (alloc_stat /=0) call error_msg_alloc

    call pack_array_i(molecule%num_atoms, enefunc%num_nonb_excl,  &
                      nonb_excl_list1, enefunc%nonb_excl_list)
    deallocate(nonb_excl_list1, stat = dealloc_stat)


    max_nonb_excl_num = max(1,maxval(enefunc%num_nonb_excl(1:nsolute)))
    allocate(enefunc%table%nonb_excl_list(max_nonb_excl_num,nsolute), &
             stat = alloc_stat)
    if (alloc_stat /=0) call error_msg_alloc

    enefunc%table%nonb_excl_list(1:max_nonb_excl_num,1:nsolute) = 0
    do i = 1, nsolute
      enefunc%table%nonb_excl_list(1:enefunc%num_nonb_excl(i),i) = &
             nonb_excl_list(1:enefunc%num_nonb_excl(i),i)
    end do
    deallocate(nonb_excl_list, stat = dealloc_stat)


    max_nb14_num = max(1,maxval(enefunc%num_nb14_calc(1:nsolute)))
    allocate(enefunc%nb14_calc_list(max_nb14_num,nsolute), &
             enefunc%nb14_qq_scale (max_nb14_num,nsolute), &
             enefunc%nb14_lj_scale (max_nb14_num,nsolute), &
             stat = alloc_stat)
    if (alloc_stat /=0) call error_msg_alloc

    enefunc%nb14_calc_list(1:max_nb14_num,1:nsolute) = 0
    enefunc%nb14_qq_scale (1:max_nb14_num,1:nsolute) = 0.0_wp
    enefunc%nb14_lj_scale (1:max_nb14_num,1:nsolute) = 0.0_wp
    do i = 1, nsolute
      enefunc%nb14_calc_list(1:enefunc%num_nb14_calc(i),i) = &
              nb14_calc_list(1:enefunc%num_nb14_calc(i),i)
      enefunc%nb14_qq_scale (1:enefunc%num_nb14_calc(i),i) = &
              nb14_qq_scale (1:enefunc%num_nb14_calc(i),i)
      enefunc%nb14_lj_scale (1:enefunc%num_nb14_calc(i),i) = &
              nb14_lj_scale (1:enefunc%num_nb14_calc(i),i)
    end do


    deallocate(nb14_calc_list, &
               nb14_qq_scale,  &
               nb14_lj_scale, stat = dealloc_stat)

    return

  end subroutine count_nonb_excl_water

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    count_nonb_excl_c19
  !> @brief        exclude 1-2, 1-3 interactions for solute-solute interaction
  !! @authors      JJ, CK
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine count_nonb_excl_c19(molecule, enefunc)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: i, k, k1, k2, k3, k4, ii, nsolute
    integer                  :: ik1, ik4
    integer                  :: num_excl, num_excl_total, num_excl_total1
    integer                  :: num_nb14, num_nb14_total
    integer                  :: max_nonb_excl_num, max_nb14_num
    integer                  :: alloc_stat, dealloc_stat
    logical                  :: duplicate, improp

    integer,         pointer :: solute(:)

    integer,     allocatable :: nonb_excl_list(:,:), nonb_excl_list1(:,:)
    integer,     allocatable :: nb14_calc_list(:,:)
    real(wp),    allocatable :: nb14_qq_scale (:,:)
    real(wp),    allocatable :: nb14_lj_scale (:,:)
    integer,     allocatable :: temporary_list(:,:)


    solute  => enefunc%table%solute_list_inv
    nsolute = enefunc%table%num_solute

    ! allocate nonbonded exclusion list (solute only and all) and 
    ! 1-4 interaction list
    !
    allocate(enefunc%table%num_nonb_excl(nsolute), stat=alloc_stat)
    if (alloc_stat /=0) call error_msg_alloc
    allocate(nonb_excl_list(MAX_EXCL,nsolute),stat=alloc_stat)
    if (alloc_stat /=0) call error_msg_alloc

    allocate(enefunc%num_nonb_excl(molecule%num_atoms), stat=alloc_stat)
    if (alloc_stat /=0) call error_msg_alloc
    allocate(nonb_excl_list1(MAX_EXCL,molecule%num_atoms),stat=alloc_stat)
    if (alloc_stat /=0) call error_msg_alloc

    allocate(enefunc%num_nb14_calc(nsolute), stat=alloc_stat)
    if (alloc_stat /=0) call error_msg_alloc
    !
    ! Do not use in CHARMM19. nb14_qq parameter is defined in at_enefunc_charmm
    !
    allocate(nb14_calc_list(MAX_NB_14,nsolute), &
             nb14_qq_scale (MAX_NB_14,nsolute), &
             nb14_lj_scale (MAX_NB_14,nsolute),stat=alloc_stat)
    if (alloc_stat /=0) call error_msg_alloc

    enefunc%table%num_nonb_excl(1:nsolute) = 0
    enefunc%num_nb14_calc(1:nsolute) = 0
    enefunc%num_nonb_excl(1:molecule%num_atoms) = 0

    ! exclude 1-2 interaction
    !
    num_excl_total  = 0
    num_excl_total1 = 0

    if (enefunc%excl_level > 0) then
      do k = 1, molecule%num_bonds
     
        k1 = solute(molecule%bond_list(1,k))
        k2 = solute(molecule%bond_list(2,k))
     
        if (k1 /= 0 .and. k2 /= 0) then
     
          if (k1 < k2) then
            num_excl = enefunc%table%num_nonb_excl(k1) + 1
            enefunc%table%num_nonb_excl(k1) = num_excl
            nonb_excl_list(num_excl,k1) = k2
            num_excl_total = num_excl_total + 1
          else
            num_excl = enefunc%table%num_nonb_excl(k2) + 1
            enefunc%table%num_nonb_excl(k2) = num_excl
            nonb_excl_list(num_excl,k2) = k1
            num_excl_total = num_excl_total + 1
          end if
     
        end if
     
        k1 = molecule%bond_list(1,k)
        k2 = molecule%bond_list(2,k)
     
        if (k1 < k2) then
          num_excl = enefunc%num_nonb_excl(k1) + 1
          enefunc%num_nonb_excl(k1) = num_excl
          nonb_excl_list1(num_excl,k1) = k2
          num_excl_total1 = num_excl_total1 + 1
        else
          num_excl = enefunc%num_nonb_excl(k2) + 1
          enefunc%num_nonb_excl(k2) = num_excl
          nonb_excl_list1(num_excl,k2) = k1
          num_excl_total1 = num_excl_total1 + 1
        end if
     
      end do
    end if

    ! exclude 1-3 interaction
    !
    if (enefunc%excl_level > 1) then
      call create_nbond_list(molecule%bond_list, 3, temporary_list)
     
      do k = 1, size(temporary_list(1,:))
     
        k1 = solute(temporary_list(1,k))
        k3 = solute(temporary_list(3,k))
     
        if (k1 /= 0 .and. k3 /= 0) then 
          if (k1 < k3) then
            num_excl = enefunc%table%num_nonb_excl(k1)
            duplicate = .false.
            do i = 1, enefunc%table%num_nonb_excl(k1)
              if (k3 == nonb_excl_list(i,k1)) duplicate = .true.
            end do
            if (.not. duplicate) then
              num_excl = num_excl + 1
              enefunc%table%num_nonb_excl(k1) = num_excl
              nonb_excl_list(num_excl,k1) = k3
              num_excl_total = num_excl_total + 1
            end if
          else
            num_excl = enefunc%table%num_nonb_excl(k3)
            duplicate = .false.
            do i = 1, enefunc%table%num_nonb_excl(k3)
              if (k1 == nonb_excl_list(i,k3)) duplicate = .true.
            end do
            if (.not. duplicate) then
              num_excl = num_excl + 1
              enefunc%table%num_nonb_excl(k3) = num_excl
              nonb_excl_list(num_excl,k3) = k1
              num_excl_total = num_excl_total + 1
            end if
          end if
        end if

        k1 = temporary_list(1,k)
        k3 = temporary_list(3,k)
     
        if (k1 < k3) then
          num_excl = enefunc%num_nonb_excl(k1)
          duplicate = .false.
          do i = 1, enefunc%num_nonb_excl(k1)
            if (k3 == nonb_excl_list1(i,k1)) duplicate = .true.
          end do
          if (.not. duplicate) then
            num_excl = num_excl + 1
            enefunc%num_nonb_excl(k1) = num_excl
            nonb_excl_list1(num_excl,k1) = k3
            num_excl_total1 = num_excl_total1 + 1
          end if

        else

          num_excl = enefunc%num_nonb_excl(k3)
          duplicate = .false.
          do i = 1, enefunc%num_nonb_excl(k3)
            if (k1 == nonb_excl_list1(i,k3)) duplicate = .true.
          end do
          if (.not. duplicate) then
            num_excl = num_excl + 1
            enefunc%num_nonb_excl(k3) = num_excl
            nonb_excl_list1(num_excl,k3) = k1
            num_excl_total1 = num_excl_total1 + 1
          end if
        end if

      end do
      deallocate(temporary_list)
    end if

    ! count 1-4 interaction
    !
    num_nb14_total = 0
    if (enefunc%excl_level > 2) then

      call create_nbond_list(molecule%bond_list, 4, temporary_list)
      do k = 1, size(temporary_list(1,:))

        k1 = solute(temporary_list(1,k))
        k4 = solute(temporary_list(4,k))
        improp = .false.
        if (k1 /= 0 .and. k4 /= 0) then
          do i = 1, enefunc%num_impropers
            ik1 = enefunc%impr_list(1,i)
            ik4 = enefunc%impr_list(4,i)
            if (min(k1,k4) == min(ik1,ik4) .and. &
                max(k1,k4) == max(ik1,ik4)) then
              improp = .true.
              exit
            end if
          end do
          if (k1 < k4) then
         
            num_nb14 = enefunc%num_nb14_calc(k1)
            num_excl = enefunc%table%num_nonb_excl(k1)
            duplicate = .false.
         
            do i = 1, enefunc%table%num_nonb_excl(k1)
              if (k4 == nonb_excl_list(i,k1)) duplicate = .true.
            end do
         
            do i = 1, num_nb14
              if (k4 == nb14_calc_list(i,k1)) duplicate = .true.
            end do
         
            if (.not. duplicate) then
              if (.not. improp) then
                num_nb14 = num_nb14 + 1
                enefunc%num_nb14_calc(k1) = num_nb14
                nb14_calc_list(num_nb14,k1) = k4
                num_nb14_total = num_nb14_total + 1
              else
                num_excl = num_excl + 1
                enefunc%table%num_nonb_excl(k1) = num_excl
                nonb_excl_list(num_excl,k1) = k4
                num_excl_total = num_excl_total + 1
              end if
            end if
         
          else
         
            num_nb14 = enefunc%num_nb14_calc(k4)
            num_excl = enefunc%table%num_nonb_excl(k4)
            duplicate = .false.
         
            do i = 1, enefunc%table%num_nonb_excl(k4)
              if (k1 == nonb_excl_list(i,k4)) duplicate = .true.
            end do
         
            do i = 1, num_nb14
              if (k1 == nb14_calc_list(i,k4)) duplicate = .true.
            end do
         
            if (.not. duplicate) then
              if (.not. improp) then
                num_nb14 = num_nb14 + 1
                enefunc%num_nb14_calc(k4) = num_nb14
                nb14_calc_list(num_nb14,k4) = k1
                num_nb14_total = num_nb14_total + 1
              else
                num_excl = num_excl + 1
                enefunc%table%num_nonb_excl(k4) = num_excl
                nonb_excl_list(num_excl,k4) = k1
                num_excl_total = num_excl_total + 1
              end if
            end if
          end if
        end if

        k1 = temporary_list(1,k)
        k4 = temporary_list(4,k)
        improp = .false.
        do i = 1, enefunc%num_impropers
          ik1 = enefunc%impr_list(1,i)
          ik4 = enefunc%impr_list(4,i)
          if (min(k1,k4) == min(ik1,ik4) .and. &
              max(k1,k4) == max(ik1,ik4)) then
            improp = .true.
            exit
          end if
        end do
        if (k1 < k4) then
        
          num_excl = enefunc%num_nonb_excl(k1)
          duplicate = .false.
        
          do i = 1, enefunc%num_nonb_excl(k1)
            if (k4 == nonb_excl_list1(i,k1)) duplicate = .true.
          end do
        
          if (.not. duplicate) then
            if (improp) then
              num_excl = num_excl + 1
              enefunc%num_nonb_excl(k1) = num_excl
              nonb_excl_list1(num_excl,k1) = k4
              num_excl_total1 = num_excl_total1 + 1
            end if
          end if
        
        else
        
          num_excl = enefunc%num_nonb_excl(k4)
          duplicate = .false.
        
          do i = 1, enefunc%num_nonb_excl(k4)
            if (k1 == nonb_excl_list1(i,k4)) duplicate = .true.
          end do
        
          if (.not. duplicate) then
            if (improp) then
              num_excl = num_excl + 1
              enefunc%num_nonb_excl(k4) = num_excl
              nonb_excl_list1(num_excl,k4) = k1
              num_excl_total1 = num_excl_total1 + 1
            end if
          end if
        end if
      end do
      deallocate(temporary_list)
    end if

    ! pack 2D-array into 1D-array
    !
    allocate(enefunc%nonb_excl_list(num_excl_total1), stat = alloc_stat)
    if (alloc_stat /=0) call error_msg_alloc

    call pack_array_i(molecule%num_atoms, enefunc%num_nonb_excl,  &
                      nonb_excl_list1, enefunc%nonb_excl_list)
    deallocate(nonb_excl_list1, stat = dealloc_stat)


    max_nonb_excl_num = max(1,maxval(enefunc%num_nonb_excl(1:nsolute)))
    allocate(enefunc%table%nonb_excl_list(max_nonb_excl_num,nsolute), &
             stat = alloc_stat)
    if (alloc_stat /=0) call error_msg_alloc

    enefunc%table%nonb_excl_list(1:max_nonb_excl_num,1:nsolute) = 0
    do i = 1, nsolute
      enefunc%table%nonb_excl_list(1:enefunc%num_nonb_excl(i),i) = &
             nonb_excl_list(1:enefunc%num_nonb_excl(i),i)
    end do
    deallocate(nonb_excl_list, stat = dealloc_stat)

    max_nb14_num = max(1,maxval(enefunc%num_nb14_calc(1:nsolute)))
    allocate(enefunc%nb14_calc_list(max_nb14_num,nsolute), &
             enefunc%nb14_qq_scale (max_nb14_num,nsolute), &
             enefunc%nb14_lj_scale (max_nb14_num,nsolute), &
             stat = alloc_stat)
    if (alloc_stat /=0) call error_msg_alloc

    enefunc%nb14_calc_list(1:max_nb14_num,1:nsolute) = 0
    enefunc%nb14_qq_scale (1:max_nb14_num,1:nsolute) = 0.0_wp
    enefunc%nb14_lj_scale (1:max_nb14_num,1:nsolute) = 0.0_wp
    do i = 1, nsolute
      enefunc%nb14_calc_list(1:enefunc%num_nb14_calc(i),i) = &
              nb14_calc_list(1:enefunc%num_nb14_calc(i),i)
    end do

    deallocate(nb14_calc_list, &
               nb14_qq_scale,  &
               nb14_lj_scale, stat = dealloc_stat)

    return

  end subroutine count_nonb_excl_c19

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    count_nonb_ecqm14
  !> @brief        create ecqm_nb14_list
  !! @authors      KY
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine count_nonb_ecqm14(enefunc)

    ! formal arguments
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer :: i, j, k, i1, i2, k1, k2, k3, k4
    integer :: alloc_stat

    integer,         pointer :: solute(:)
    integer,         pointer :: ecqm_num_nb14
    integer,         pointer :: ecqm_nb14_list(:,:)

    integer                  :: qm_natoms, ec_natoms
    integer,         pointer :: qmatom_id(:), ecatom_id(:)


    qm_natoms =  enefunc%qmmm%qm_natoms
    qmatom_id => enefunc%qmmm%qmatom_id
    ec_natoms =  enefunc%qmmm%ec_natoms
    ecatom_id => enefunc%qmmm%ecatom_id

    enefunc%table%ecqm_num_nb14 = 0
    allocate(enefunc%table%ecqm_nb14_list(2, qm_natoms*ec_natoms), &
             stat=alloc_stat)
    if (alloc_stat /=0) call error_msg_alloc

    solute          => enefunc%table%solute_list_inv
    ecqm_num_nb14   => enefunc%table%ecqm_num_nb14
    ecqm_nb14_list  => enefunc%table%ecqm_nb14_list

    ecqm_num_nb14 = 0

    do i1 = 1, ec_natoms
      k1 = solute(ecatom_id(i1))

      if (enefunc%num_nb14_calc(k1) > 0) then
        do k = 1, enefunc%num_nb14_calc(k1)
          k4 = enefunc%nb14_calc_list(k,k1)
          do i2 = 1, qm_natoms
            k2 = solute(qmatom_id(i2))
            if (k2 == k4) then
              ! NOTE: k1 < k4
              ecqm_num_nb14 = ecqm_num_nb14 + 1
              ecqm_nb14_list(1,ecqm_num_nb14) = k1
              ecqm_nb14_list(2,ecqm_num_nb14) = k4
            end if
          end do
        end do
      end if

    end do

    do i1 = 1, qm_natoms
      k1 = solute(qmatom_id(i1))

      if (enefunc%num_nb14_calc(k1) > 0) then
        do k = 1, enefunc%num_nb14_calc(k1)
          k4 = enefunc%nb14_calc_list(k,k1)
          do i2 = 1, ec_natoms
            k2 = solute(ecatom_id(i2))
            if(k2 == k4) then
              ! NOTE: k1 < k4
              ecqm_num_nb14 = ecqm_num_nb14 + 1
              ecqm_nb14_list(1,ecqm_num_nb14) = k1
              ecqm_nb14_list(2,ecqm_num_nb14) = k4
            end if
          end do
        end do
      end if

    end do

    !dbg write(MsgOut,'("QM-EC_14 pair")')
    !dbg write(MsgOut,'("ecqm_num_nb14=",i8)') enefunc%table%ecqm_num_nb14
    !dbg do i = 1, ecqm_num_nb14
    !dbg   write(MsgOut,'("k1,k4=",2i4)') ecqm_nb14_list(:,i)
    !dbg end do

    return

  end subroutine count_nonb_ecqm14

end module at_enefunc_table_mod
