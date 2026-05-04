!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   biweight_mod
!> @brief   libralies for calculating the bi-weight function
!! @authors Junichiro Koga (JK)
!
!  (c) Copyright 2020 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module biweight_mod

  use constants_mod

  implicit none
  private

  public :: biweight

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      biweight
  !> @brief        calculate the biweight function, 
  !                G(x) = a*(1-x**2/c**2)**2 if |x| < c, 0 otherwise.
  !! @authors      JK
  !! @param[in]    x    : the independent variable
  !! @param[in]    a    : the parameter a (optional) defaults to 48/41
  !! @param[in]    c    : the parameter c (optional) defaults to 2
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8
  function biweight(x, a, c) result(y)
    real(wp), intent(in)           :: x
    real(wp), intent(in), optional :: a
    real(wp), intent(in), optional :: c
    real(wp) :: y
    real(wp) :: parama, paramc, paramc2
    parama = 48.0_wp/41.0_wp
    if (present(a)) parama = a
    paramc  = 2.0_wp
    paramc2 = 4.0_wp
    if (present(c)) then
      paramc  = c
      paramc2 = c * c
    endif

    y = 0.0_wp
    if (abs(x)>paramc) then
      return
    endif

    y = parama * (1.0_wp - x*x/paramc2)**2

    return
  end function biweight

end module biweight_mod
