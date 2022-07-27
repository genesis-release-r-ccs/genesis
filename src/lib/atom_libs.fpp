!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   atom_lib_mod
!> @brief   utilities of atom libraries
!! @authors Daisuke Matsuoka (DM), Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module atom_libs_mod

  use messages_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public :: atomic_number
  public :: atomic_number_by_name

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_atomic_number
  !> @brief        define atomic number from molecule
  !! @authors      TM
  !! @param[in]    mass      : mass
  !! @param[in]    atom_type : atom class name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function atomic_number(mass, atom_type)

    ! formal arguments
    real(wp),              intent(in)    :: mass
    character(6),          intent(inout) :: atom_type

    ! local variables
    integer                              :: atomic_number
    integer                              :: z, int_mass

    z = 0
    int_mass = nint(mass)

    if      (int_mass == 1  ) then  ! H
      z = 1
    else if (int_mass == 7  ) then  ! Li
      z = 3
    else if (int_mass == 12 ) then  ! C
      z = 6
    else if (int_mass == 14 ) then  ! N
      z = 7
    else if (int_mass == 16 ) then  ! O
      z = 8
    else if (int_mass == 19 ) then  ! F
      z = 9
    else if (int_mass == 23 ) then  ! Na
      z = 11
    else if (int_mass == 24 ) then  ! Mg
      z = 12
    else if (int_mass == 27 ) then  ! Al
      z = 13
    else if (int_mass == 31 ) then  ! P
      z = 15
    else if (int_mass == 32 ) then  ! S
      z = 16
    else if (int_mass == 35 ) then  ! Cl
      z = 17
    else if (int_mass == 39 ) then  ! K
      z = 19
    else if (int_mass == 40 ) then  ! Ca
      z = 20
    else if (int_mass == 56 ) then  ! Fe
      z = 26
    else if (int_mass == 64 ) then  ! Cu
      z = 29
    else if (int_mass == 65 ) then  ! Zn
      z = 30
    else if (int_mass == 79 ) then  ! Se
      z = 34
    else if (int_mass == 80 ) then  ! Br
      z = 35
    else if (int_mass == 85 ) then  ! Rb
      z = 37
    else if (int_mass == 112) then  ! Cd
      z = 48
    else if (int_mass == 127) then  ! I
      z = 53
    else if (int_mass == 133) then  ! Cs
      z = 55
    else if (int_mass == 137) then  ! Ba
      z = 56
    else
      z = atomic_number_by_name(atom_type)
    end if

    atomic_number = z

    return

  end function atomic_number

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      atomic_number_by_name
  !> @brief        get atomic number of an atom
  !! @authors      DM
  !! @param[in]    atom_type   : atom type
  !! @return       atomic number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
              
  function atomic_number_by_name(atom_type)

    ! formal variables
    character(6), intent(in) :: atom_type

    ! local variable
    integer                  :: atomic_number_by_name
    integer                  :: z
    character(3)             :: atom_element


    atom_element = adjustl(atom_type(1:3))

    if      (atom_element(1:1) .eq. 'H'  .or. atom_element(1:1) .eq. 'h'  .or. &
             atom_element(1:2) .eq. 'TH') then       ! amber param19_ipq.dat
      z = 1
  
    else if (atom_element(1:2) .eq. 'Na' .or. atom_element(1:3) .eq. 'SOD') then
      z = 11

    else if (atom_element(1:2) .eq. 'MG') then
      z = 12

    else if (atom_element(1:1) .eq. 'K'  .or. atom_element(1:3) .eq. 'POT') then
      z = 19
  
    else if (atom_element(1:2) .eq. 'Cl' .or.   &     ! parm10.dat
             atom_element(1:2) .eq. 'cl' .or.   &     ! gaff.dat
             atom_element(1:3) .eq. 'CLA') then
      z = 17

    else if (atom_element(1:2) .eq. 'Cu'  .or.  &
             atom_element(1:2) .eq. 'CU') then
      z = 29

    else if (atom_element(1:2) .eq. 'C0'   .or. &     ! amber param10.dat
             atom_element(1:3) .eq. 'CAL') then       ! CHARMM
      z = 20

    else if (atom_element(1:2) .eq. 'FE') then
      z = 26

    else if (atom_element(1:2) .eq. 'ZN' .or. atom_element(1:2) .eq. 'Zn') then
      z = 30

    else if (atom_element(1:1) .eq. 'C'  .or. atom_element(1:1) .eq. 'c') then
      z = 6

    else if (atom_element(1:2) .eq. '2C' .or. atom_element(1:2) .eq. '3C' .or. &
             atom_element(1:2) .eq. 'TG' .or. atom_element(1:2) .eq. 'TJ' .or. &
             atom_element(1:2) .eq. 'TP' .or. atom_element(1:2) .eq. 'TM' .or. &
             atom_element(1:2) .eq. 'TA') then    ! amber param19_ipq.dat
      z = 6
  
    else if (atom_element(1:1) .eq. 'N'  .or. atom_element(1:1) .eq. 'n') then
      z = 7
  
    else if (atom_element(1:2) .eq. 'TN') then    ! amber param14ipq
      z = 7

    else if (atom_element(1:1) .eq. 'O' .or. atom_element(1:1) .eq. 'o') then
      z = 8

    else if (atom_element(1:1) .eq. 'F' .or. atom_element(1:1) .eq. 'f') then
      z = 9
  
    else if (atom_element(1:1) .eq. 'S' .or. atom_element(1:1) .eq. 's') then
      z = 16
  
    else if (atom_element(1:1) .eq. 'P' .or. atom_element(1:1) .eq. 'p') then
      z = 15

    ! lone pair of TIP4P/TiP5P
    else if (atom_element(1:2) .eq. 'LP' .or.  &   ! charmm
             atom_element(1:2) .eq. 'EP' .or.  &   ! amber
             atom_element(1:2) .eq. 'IW') then     ! gromos
      z = 0

    else
      call error_msg('Atomic_Number_By_Name> Unknown atom element ['// &
           trim(atom_type)//']')
  
    end if

    atomic_number_by_name = z

    return

  end function atomic_number_by_name

end module atom_libs_mod
