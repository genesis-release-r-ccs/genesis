!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   math_libs_mod
!> @brief   utilities of math libraries in ATDYN
!! @authors Yuji Sugita (YS), Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module math_libs_mod

  use constants_mod
  use random_mod

  implicit none
  private

  ! constants

  ! subroutines
  public  :: pack_array_i
  public  :: pack_array_r
  public  :: erf_n
  public  :: erfc04
  public  :: b_spline_coef
  public  :: b_spline_dev_coef
  public  :: compute_inverse_matrix
  public  :: get_ewald_alpha
  public  :: factorization_235
  public  :: powersinh
  public  :: powersinh_double
  public  :: sum_gauss
  private :: factorization_prime
  public  :: cubic_spline
  public  :: cubic_spline_coeff
  public  :: lcm

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pack_array_i
  !> @brief        pack 2D-array into 1D-array
  !! @authors      YS
  !! @param[in]    num_ele  : number of rows in 2d list
  !! @param[in]    num_list : number of columns in 2d list
  !! @param[in]    list_2d  : 2d list
  !! @param[out]   list_1d  : 1d list
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pack_array_i(num_ele, num_list, list_2d, list_1d)

    ! formal arguments
    integer,                 intent(in)    :: num_ele
    integer,                 intent(in)    :: num_list(:)
    integer,                 intent(in)    :: list_2d(:,:)
    integer,                 intent(inout) :: list_1d(:)

    ! local variables
    integer                  :: i, j, k


    k = 0
    do i = 1, num_ele
      do j = 1, num_list(i)
        k = k + 1
        list_1d(k) = list_2d(j,i)
      end do
    end do

    return

  end subroutine pack_array_i

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pack_array_r
  !> @brief        pack 2D-array into 1D-array
  !! @authors      YS
  !! @param[in]    num_ele  : number of rows in 2d list
  !! @param[in]    num_list : number of columns in 2d list
  !! @param[in]    list_2d  : 2d list
  !! @param[out]   list_1d  : 1d list
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pack_array_r(num_ele, num_list, list_2d, list_1d)

    ! formal arguments
    integer,                 intent(in)    :: num_ele
    integer,                 intent(in)    :: num_list(:)
    real(wp),                intent(in)    :: list_2d(:,:)
    real(wp),                intent(inout) :: list_1d(:)

    ! local variables
    integer                  :: i, j, k


    k = 0
    do i = 1, num_ele
      do j = 1, num_list(i)
        k = k + 1
        list_1d(k) = list_2d(j,i)
      end do
    end do

    return

  end subroutine pack_array_r

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      erf_n
  !> @brief        error function
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function erf_n (x) result(erf)

    ! formal arguments
    real(wp),                intent(in)    :: x

    ! return value
    real(wp)                 :: erf


    erf = 1.0_wp - erfc04(x)

    return

  end function erf_n

!-------------------------------------------------------------------------------
! Complementary error function
! Taken from SUN's FDLIBM version 5.2 and translated from c to fortran.
!-------------------------------------------------------------------------------
! Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
!
! Developed at SunSoft, a Sun Microsystems, Inc. business.
! Permission to use, copy, modify, and distribute this
! software is freely granted, provided that this notice
! is preserved.
!-------------------------------------------------------------------------------
! Definition:
!------------
!                       x
!                2      |\
! erf(x)  =  ---------  | exp(-t*t)dt
!             sqrt(pi) \|
!                       0
!
! erfc(x) =  1 - erf(x)
!
! Note that  erf(-x) = -erf(x)
!            erfc(-x) = 2 - erfc(x)
!
!
! Method:
!--------
!
! 1. For |x| in [0, 0.84375]
!    erf(x)  = x + x*R(x^2)
!    erfc(x) = 1 - erf(x)           if x in [-.84375,0.25]
!            = 0.5 + ((0.5-x)-x*R)  if x in [0.25,0.84375]
!    where R = P/Q where P is an odd poly of degree 8 and
!                        Q is an odd poly of degree 10.
!                                                 -57.90
!                        | R - (erf(x)-x)/x | <= 2
!
!
!    Remark. The formula is derived by noting
!            erf(x) = (2/sqrt(pi))*(x - x^3/3 + x^5/10 - x^7/42 + ....)
!    and that
!            2/sqrt(pi) = 1.128379167095512573896158903121545171688
!    is close to one. The interval is chosen because the fix
!    point of erf(x) is near 0.6174 (i.e., erf(x)=x when x is
!    near 0.6174), and by some experiment, 0.84375 is chosen to
!    guarantee the error is less than one ulp for erf.
!
! 2. For |x| in [0.84375,1.25], let s = |x| - 1, and c = 0.84506291151
!    rounded to single (24 bits)
!    erf(x)  = sign(x) * (c  + P1(s)/Q1(s))
!    erfc(x) = (1-c)  - P1(s)/Q1(s) if x > 0
!              1+(c+P1(s)/Q1(s))    if x < 0
!              |P1/Q1 - (erf(|x|)-c)| <= 2**-59.06
!
!    Remark: here we use the taylor series expansion at x=1.
!            erf(1+s) = erf(1) + s*Poly(s)
!                     = 0.845.. + P1(s)/Q1(s)
!            That is, we use rational approximation to approximate
!            erf(1+s) - (c = (single)0.84506291151)
!    Note that |P1/Q1|< 0.078 for x in [0.84375,1.25] where
!    P1(s) = degree 6 poly in s
!    Q1(s) = degree 6 poly in s
!
! 3. For x in [1.25,1/0.35(~2.857143)],
!    erfc(x) = (1/x)*exp(-x*x-0.5625+R1/S1)
!    erf(x)  = 1 - erfc(x)
!    where
!    R1(z) = degree 7 poly in z, (z=1/x^2)
!    S1(z) = degree 8 poly in z
!
! 4. For x in [1/0.35,28]
!    erfc(x) = (1/x)*exp(-x*x-0.5625+R2/S2) if x > 0
!            = 2.0 - (1/x)*exp(-x*x-0.5625+R2/S2) if -6<x<0
!            = 2.0 - tiny(if x <= -6)
!    erf(x)  = sign(x)*(1.0 - erfc(x)) if x < 6, else
!    erf(x)  = sign(x)*(1.0 - tiny)
!    where
!    R2(z) = degree 6 poly in z, (z=1/x^2)
!    S2(z) = degree 7 poly in z
!
!    Note1: To compute exp(-x*x-0.5625+R/S), let s be a single
!           precision number and s := x; then
!           -x*x = -s*s + (s-x)*(s+x)
!            exp(-x*x-0.5626+R/S) = exp(-s*s-0.5625)*exp((s-x)*(s+x)+R/S)
!    Note2: Here 4 and 5 make use of the asymptotic series
!                      exp(-x*x)
!           erfc(x) ~ ---------- * ( 1 + Poly(1/x^2) )
!                     x*sqrt(pi)
!           We use rational approximation to approximate
!           g(s) = f(1/x^2) = log(erfc(x)*x) - x*x + 0.5625
!           Here is the error bound for R1/S1 and R2/S2
!           |R1/S1 - f(x)|  < 2**(-62.57)
!           |R2/S2 - f(x)|  < 2**(-61.52)
!
! 5. For inf > x >= 28
!    erf(x)  = sign(x) *(1 - tiny)  (raise inexact)
!    erfc(x) = tiny*tiny (raise underflow) if x > 0
!            = 2 - tiny if x<0
!
! 7. Special case:
!    erf(0)  = 0, erf(inf)  = 1, erf(-inf) = -1,
!    erfc(0) = 1, erfc(inf) = 0, erfc(-inf) = 2,
!    erfc/erf(NaN) is NaN
!-------------------------------------------------------------------------------

  function erfc04(x) result(erfc)

    real(wp),         intent(in)    :: x
    real(dp)                        :: erfc
    real(dp)                        :: ax,p,q,r,s,y,z

    real(dp),         parameter     :: zero1=  0.0_wp
    real(dp),         parameter     :: half1=  0.5_wp
    real(dp),         parameter     :: one1 =  1.0_wp
    real(dp),         parameter     :: two1 =  2.0_wp
    real(dp),         parameter     :: erx  =  8.45062911510467529297e-01_dp
    ! Coefficients for approximation to  erf on [0,0.84375]
    real(dp),         parameter     :: efx  =  1.28379167095512586316e-01_dp
    real(dp),         parameter     :: efx8 =  1.02703333676410069053e+00_dp
    real(dp),         parameter     :: pp0  =  1.28379167095512558561e-01_dp
    real(dp),         parameter     :: pp1  = -3.25042107247001499370e-01_dp
    real(dp),         parameter     :: pp2  = -2.84817495755985104766e-02_dp
    real(dp),         parameter     :: pp3  = -5.77027029648944159157e-03_dp
    real(dp),         parameter     :: pp4  = -2.37630166566501626084e-05_dp
    real(dp),         parameter     :: qq1  =  3.97917223959155352819e-01_dp
    real(dp),         parameter     :: qq2  =  6.50222499887672944485e-02_dp
    real(dp),         parameter     :: qq3  =  5.08130628187576562776e-03_dp
    real(dp),         parameter     :: qq4  =  1.32494738004321644526e-04_dp
    real(dp),         parameter     :: qq5  = -3.96022827877536812320e-06_dp
    ! Coefficients for approximation to  erf  in [0.84375,1.25]
    real(dp),         parameter     :: pa0  = -2.36211856075265944077e-03_dp
    real(dp),         parameter     :: pa1  =  4.14856118683748331666e-01_dp
    real(dp),         parameter     :: pa2  = -3.72207876035701323847e-01_dp
    real(dp),         parameter     :: pa3  =  3.18346619901161753674e-01_dp
    real(dp),         parameter     :: pa4  = -1.10894694282396677476e-01_dp
    real(dp),         parameter     :: pa5  =  3.54783043256182359371e-02_dp
    real(dp),         parameter     :: pa6  = -2.16637559486879084300e-03_dp
    real(dp),         parameter     :: qa1  =  1.06420880400844228286e-01_dp
    real(dp),         parameter     :: qa2  =  5.40397917702171048937e-01_dp
    real(dp),         parameter     :: qa3  =  7.18286544141962662868e-02_dp
    real(dp),         parameter     :: qa4  =  1.26171219808761642112e-01_dp
    real(dp),         parameter     :: qa5  =  1.36370839120290507362e-02_dp
    real(dp),         parameter     :: qa6  =  1.19844998467991074170e-02_dp
    ! Coefficients for approximation to  erfc in [1.25,1/0.35]
    real(dp),         parameter     :: ra0  = -9.86494403484714822705e-03_dp
    real(dp),         parameter     :: ra1  = -6.93858572707181764372e-01_dp
    real(dp),         parameter     :: ra2  = -1.05586262253232909814e+01_dp
    real(dp),         parameter     :: ra3  = -6.23753324503260060396e+01_dp
    real(dp),         parameter     :: ra4  = -1.62396669462573470355e+02_dp
    real(dp),         parameter     :: ra5  = -1.84605092906711035994e+02_dp
    real(dp),         parameter     :: ra6  = -8.12874355063065934246e+01_dp
    real(dp),         parameter     :: ra7  = -9.81432934416914548592e+00_dp
    real(dp),         parameter     :: sa1  =  1.96512716674392571292e+01_dp
    real(dp),         parameter     :: sa2  =  1.37657754143519042600e+02_dp
    real(dp),         parameter     :: sa3  =  4.34565877475229228821e+02_dp
    real(dp),         parameter     :: sa4  =  6.45387271733267880336e+02_dp
    real(dp),         parameter     :: sa5  =  4.29008140027567833386e+02_dp
    real(dp),         parameter     :: sa6  =  1.08635005541779435134e+02_dp
    real(dp),         parameter     :: sa7  =  6.57024977031928170135e+00_dp
    real(dp),         parameter     :: sa8  = -6.04244152148580987438e-02_dp
    ! Coefficients for approximation to  erfc in [1/.35,28]
    real(dp),         parameter     :: rb0  = -9.86494292470009928597e-03_dp
    real(dp),         parameter     :: rb1  = -7.99283237680523006574e-01_dp
    real(dp),         parameter     :: rb2  = -1.77579549177547519889e+01_dp
    real(dp),         parameter     :: rb3  = -1.60636384855821916062e+02_dp
    real(dp),         parameter     :: rb4  = -6.37566443368389627722e+02_dp
    real(dp),         parameter     :: rb5  = -1.02509513161107724954e+03_dp
    real(dp),         parameter     :: rb6  = -4.83519191608651397019e+02_dp
    real(dp),         parameter     :: sb1  =  3.03380607434824582924e+01_dp
    real(dp),         parameter     :: sb2  =  3.25792512996573918826e+02_dp
    real(dp),         parameter     :: sb3  =  1.53672958608443695994e+03_dp
    real(dp),         parameter     :: sb4  =  3.19985821950859553908e+03_dp
    real(dp),         parameter     :: sb5  =  2.55305040643316442583e+03_dp
    real(dp),         parameter     :: sb6  =  4.74528541206955367215e+02_dp
    real(dp),         parameter     :: sb7  = -2.24409524465858183362e+01_dp

    ax = abs(x)

    if (ax < 0.84375_dp) then

      if (ax < epsilon(x)) then
        erfc = one1 - x
      else

        z = x**2
        r = pp0 + z*(pp1 + z*(pp2 + z*(pp3 + z*pp4)))
        s = one1 + z*(qq1 + z*(qq2 + z*(qq3 + z*(qq4 + z*qq5))))
        y = r/s

        if (x < 0.25_dp) then
          erfc =  one1 - (x + x*y)
        else
          r = x*y
          r = r + (x - half1)
          erfc = half1 - r
        end if

      end if

    else if (ax < 1.25_dp) then

      s = ax - one1
      p = pa0 + s*(pa1 + s*(pa2 + s*(pa3 + s*(pa4 + s*(pa5 + s*pa6)))))
      q = one1 + s*(qa1 + s*(qa2 + s*(qa3 + s*(qa4 + s*(qa5 + s*qa6)))))

      if (x > zero1) then
        z  = one1 - erx
        erfc = z - p/q
      else
        z = erx + p/q
        erfc = one1 + z
      end if

    else if (ax < 28.0_dp) then

      s = one1/(ax**2)

      if (ax < 2.857143_dp) then
        p = ra0 + s*(ra1 + s*(ra2 + s*(ra3 + s*(ra4 +   &
             s*(ra5 + s*(ra6 + s*ra7))))))
        q = one1 + s*(sa1 + s*(sa2 + s*(sa3 + s*(sa4 +   &
             s*(sa5 + s*(sa6 + s*(sa7 + s*sa8)))))))
      else

        if (x < -6.0_dp) then

          erfc = two1
          return

        end if

        p = rb0 + s*(rb1 + s*(rb2 + s*(rb3 + s*(rb4 + s*(rb5 + s*rb6)))))
        q = one1 + s*(sb1 + s*(sb2 + s*(sb3 + s*(sb4 +   &
             s*(sb5 + s*(sb6 + s*sb7))))))

      end if

      z = ax
      r = exp(-z**2-0.5625_dp+(z-ax)*(z+ax)+p/q)

      if (x > zero1) then
        erfc = r/x
      else
        erfc = two1 + r/x
      end if

    else

      if (x > zero1) then
        erfc = zero1
      else
        erfc = two1
      end if

    end if

    return

  end function erfc04

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    b_spline_coef
  !> @brief        Calculate B-spline coefficients
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine b_spline_coef(n, u, M)

    ! formal arguments
    integer,                 intent(in)    :: n
    real(wp),                intent(in)    :: u
    real(wp),                intent(out)   :: M(1:n)

    ! local variables
    integer   :: i,j
    real(wp)  :: u_3, u_2


    if (n == 4) then

      u_3 = u*u*u
      u_2 = u*u
      M(1) = u_3
      M(2) = -3.0_wp*u_3 + 3.0_wp*u_2 + 3.0_wp*u + 1.0_wp
      M(3) = 3.0_wp*u_3 - 6.0_wp*u_2 + 4.0_wp
      M(4) = (1.0_wp-u)*(1.0_wp-u)*(1.0_wp-u)

    else

      M(1) = u
      M(2) = 1.0_wp - u

      do j = 3, n
        M(j) = M(j-1) * (1.0_wp-u)
        do i = 1, j - 2
          M(j - i) = (real(j-i-1,wp) + u) * M(j-i) &
                   + (real(i+1,wp) - u) * M(j-i-1)
        end do
        M(1) = M(1) * u
      end do

    end if

    return

  end subroutine b_spline_coef

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    b_spline_dev_coef
  !> @brief        Calculate B-spline coefficients and its derivatives
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine b_spline_dev_coef(n, u, M, dM)

    ! formal arguments
    integer,                 intent(in)    :: n
    real(wp),                intent(in)    :: u
    real(wp),                intent(out)   :: M(1:n)
    real(wp),                intent(out)   :: dM(1:n)

    ! local variables
    integer                  :: i,j
    real(wp)                 :: u_5, u_4, u_3, u_2, v_1, v_2, v_4, v_5


    if (n == 4) then

      u_3 = u*u*u
      u_2 = u*u
      M(1) = u_3
      M(2) = -3.0_wp*u_3 + 3.0_wp*u_2 + 3.0_wp*u + 1.0_wp
      M(3) = 3.0_wp*u_3 - 6.0_wp*u_2 + 4.0_wp
      M(4) = (1.0_wp-u)*(1.0_wp-u)*(1.0_wp-u)
      dM(1) = u_2
      dM(2) = -3.0_wp*u_2 + 2.0_wp*u + 1.0_wp
      dM(3) = 3.0_wp*u_2 - 4.0_wp*u
      dM(4) = -(1.0_wp-u)*(1.0_wp-u)

    else if (n == 6) then

      u_2 = u*u
      u_3 = u_2*u
      u_4 = u_3*u
      u_5 = u_4*u
      v_1 = 1.0_wp - u
      v_2 = v_1*v_1
      v_4 = v_2*v_2
      v_5 = v_4*v_1
      M(1) = u_5
      M(2) = -5.0_wp*u_5+5.0_wp*u_4+10.0_wp*u_3+10.0_wp*u_2+5.0_wp*u+1.0_wp
      M(3) = 10.0_wp*u_5-20.0_wp*u_4-20.0_wp*u_3+20.0_wp*u_2+50.0_wp*u+26.0_wp
      M(4) = -10.0_wp*u_5+30.0_wp*u_4-60.0_wp*u_2+66.0_wp
      M(5) = 5.0_wp*u_5-20.0_wp*u_4+20.0_wp*u_3+20.0_wp*u_2-50.0_wp*u+26.0_wp
      M(6) = v_5
      dM(1) = u_4
      dM(2) = -5.0_wp*u_4+4.0_wp*u_3+6.0_wp*u_2+4.0_wp*u+1.0_wp
      dM(3) = 10.0_wp*u_4-16.0_wp*u_3-12.0_wp*u_2+8.0_wp*u+10.0_wp
      dM(4) = -10.0_wp*u_4+24.0_wp*u_3-24.0_wp*u
      dM(5) = 5.0_wp*u_4-16.0_wp*u_3+12.0_wp*u_2+8.0_wp*u-10.0_wp
      dM(6) = -v_4

    else

      M(1) = u
      M(2) = 1.0_wp - u

      do j = 3, n-1
        M(j) = M(j-1) * (1.0_wp-u)
        do i = 1, j-2
          M(j-i) = (real(j-i-1,wp)+u) * M(j-i) &
                 + (real(i+1,wp)-u) * M(j-i-1)
        end do
        M(1) = M(1) * u
      end do

      dM(1) = M(1)
      do i = 2, n - 1
        dM(i) = M(i) - M(i-1)
      end do
      dM(n) = - M(n-1)

      M(n) = M(n-1) * (1.0_wp-u)
      do i = 1, n - 2
        M(n-i) = (real(n-i-1,wp)+u) * M(n-i) &
               + (real(i+1,wp)-u) * M(n-i-1)
      end do
      M(1) = M(1) * u

    end if

    return

  end subroutine b_spline_dev_coef

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_inverse_matrix
  !> @brief        Gauss-Jordan elimination method with partial pivoting
  !! @authors      TM
  !! @param[in]    n     : matrix dimension
  !! @param[in]    m     : original matrix
  !! @param[out]   inv_m : inverse matrix
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_inverse_matrix(n, m, inv_m)

    ! formal arguments
    integer,      intent(in)    :: n
    real(wp),     intent(in)    :: m(n,n)
    real(wp),     intent(out)   :: inv_m(n,n)

    ! local variables
    integer                     :: i, j, k, imax
    real(wp)                    :: a(n,2*n), amax, swap, coef, pivot


    do i = 1, n
      do j = 1, n
        a(i,j) = m(i,j)
        if (i == j) a(i,j+n) = 1.0_wp
        if (i /= j) a(i,j+n) = 0.0_wp
      end do
    end do

    do k = 1, n
      amax = a(k,k)
      imax = k
      do i = 1, n
        if (abs(a(i,k)) > amax) then
          amax = a(i,k)
          imax = i
        end if
      end do

      do j = k, 2*n
        swap = a(k,j)
        a(k,j) = a(imax,j)
        a(imax,j) = swap
      end do

      pivot = 1.0_wp/a(k,k)
      do j = k, 2*n
        a(k,j) = a(k,j)*pivot
      end do
      do i = 1, n
        if (i /= k) then
          coef = a(i,k)
          do j = k, 2*n
            a(i,j) = a(i,j) - a(k,j)*coef
          end do
        end if
      end do
    end do

    do i = 1, n
      do j = 1, n
        inv_m(i,j) = a(i,j+n)
      end do
    end do

    return

  end subroutine compute_inverse_matrix

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      get_ewald_alpha
  !> @brief        calculate ewald alpha value from cutoff and tolerance
  !! @authors      MK
  !! @param[in]    cutoff : cutoff distance
  !! @param[in]    tol    : tolerance
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_ewald_alpha(cutoff, tol) result(alpha)

    real(wp), intent(in)  :: cutoff
    real(wp), intent(in)  :: tol
    real(wp)              :: alpha

    integer, parameter :: maxit = 100
    integer :: i
    real(wp) :: ew_lo, ew_hi


    alpha = 1.0_wp
    do while( erfc04( alpha * cutoff ) / cutoff >= tol )
      alpha = 2.0_wp * alpha
    end do

    ew_lo = 0.0_wp
    ew_hi = alpha
    do i = 1, maxit
      alpha = 0.5_wp * ( ew_lo + ew_hi )
      if ( erfc04( alpha * cutoff ) / cutoff >= tol ) then
        ew_lo = alpha
      else
        ew_hi = alpha
      end if
    end do

    return

  end function get_ewald_alpha

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      factorization_235
  !> @brief        prime factroization of 235
  !! @authors      CK
  !! @param[in]    int_number : int_number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function factorization_235(int_number) result(prime_flag)

    integer, intent(in)   :: int_number
    logical               :: prime_flag

    integer               :: itmp


    itmp = int_number
    prime_flag=.false.
    call factorization_prime(5, itmp, prime_flag)
    if (prime_flag) return
    call factorization_prime(3, itmp, prime_flag)
    if (prime_flag) return
    call factorization_prime(2, itmp, prime_flag)

    return

  end function factorization_235

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    factorization_prime
  !> @brief        prime factroization of prime
  !! @authors      CK
  !! @param[in]    prime_number : prime_number
  !! @param[inout] int_number   : int_number
  !! @param[out]   end_flag     : end_flag
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine factorization_prime(prime_number, int_number, end_flag)

    ! formal arguments
    integer, intent(in)      :: prime_number
    integer, intent(inout)   :: int_number
    logical, intent(inout)   :: end_flag

    ! local variables
    integer                  :: mod_num
    logical                  :: run_factorization


    run_factorization = .true.

    do while(run_factorization)

     mod_num = mod(int_number,prime_number)

     if (mod_num == 0) then
       int_number = int(int_number/prime_number)
       if (int_number == 1) then
         end_flag=.true.
         run_factorization=.false.
       end if
     else
       run_factorization=.false.
     end if

    end do

    return

  end subroutine factorization_prime

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      powersinh
  !> @brief        Power seriesa expansion of sinh(x)/x up to tenth order
  !! @authors      TA, JJ
  !! @param[in]    x : value x, which is vector
  !! @param[in]    n : size of x
  !! @return       y : approx. value of sinh(x)/x
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function powersinh(x) result(y)

    ! return values
    real(wip)               :: y(3)

    ! formal arguments
    real(wip),   intent(in) :: x(3)

    ! local variables
    integer                 :: j, k
    real(wip)               :: a(0:5)


    a(0)  = 1.0_wip
    a(1)  = a(0) / (2.0_wip*3.0_wip)
    a(2)  = a(1) / (4.0_wip*5.0_wip)
    a(3)  = a(2) / (6.0_wip*7.0_wip)
    a(4)  = a(3) / (8.0_wip*9.0_wip)
    a(5)  = a(4) / (10.0_wip*11.0_wip)

    y(1:3) = 0.0_wip
    do j = 0, 5
      do k = 1, 3
        if (real(2*j,wp)*log(abs(x(k))) > -15.0_wp) &
        y(k) = a(j)*x(k)**(2*j) + y(k)
      end do
    end do

  end function powersinh

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      powersinh_double
  !> @brief        Power seriesa expansion of sinh(x)/x up to tenth order
  !! @authors      TA, JJ
  !! @param[in]    x : value x, which is vector
  !! @param[in]    n : size of x
  !! @return       y : approx. value of sinh(x)/x
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function powersinh_double(x) result(y)

    ! return values
    real(dp)                :: y(3)

    ! formal arguments
    real(dp),    intent(in) :: x(3)

    ! local variables
    integer                 :: j, k
    real(dp)                :: a(0:5)


    a(0)  = 1.0_wip
    a(1)  = a(0) / (2.0_wip*3.0_wip)
    a(2)  = a(1) / (4.0_wip*5.0_wip)
    a(3)  = a(2) / (6.0_wip*7.0_wip)
    a(4)  = a(3) / (8.0_wip*9.0_wip)
    a(5)  = a(4) / (10.0_wip*11.0_wip)

    y(1:3) = 0.0_wip
    do j = 0, 5
      do k = 1, 3
        y(k) = a(j)*x(k)**(2*j) + y(k)
      end do
    end do

  end function powersinh_double

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      sum_gauss
  !> @brief        returns the sum of n independent gaussian noises squared
  !! @authors      TA, JJ
  !! @param[in]    n : number of noises
  !! @return       the sum of n independent gaussian noises squared
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function sum_gauss(n)

    ! return values
    real(wip)                :: sum_gauss

    ! formal arguments
    integer(iintegers),  intent(in) :: n
    real(wip)               :: random_temp1, random_temp


    if( mod(n,2) == 0 ) then
      sum_gauss = 2.0_wip*gamma_distr(n/2)
    else
      random_temp1 = 2.0_wip*gamma_distr((n-1)/2)
      random_temp  = random_get_gauss()
      sum_gauss    = random_temp1 + random_temp*random_temp
    end if

  end function sum_gauss

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      gamma_distr
  !> @brief        gamma-distributed random number
  !! @authors      TA, JJ, CK
  !! @param[in]    n  : integer order
  !! @return       Returns a deviate distributed as a gamma distribution
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function gamma_distr(n)

    ! return values
    real(wip)                :: gamma_distr

    ! formal arguments
    integer(iintegers), intent(in)   :: n

    ! local variables
    integer                  :: j
    real(wip)                :: random1, random2, random3
    real(wip)                :: a, b, c, d, m


    do while(.true.)

      random1 = 2.0_wip*random_get() - 1.0_wip
      random2 = 2.0_wip*random_get() - 1.0_wip
      if (random1*random1 + random2*random2 > 1.0_wip) cycle

      a = random2 / random1
      m = real(n-1,wip)
      b = sqrt(2.0_wip*m + 1.0_wip)
      c = b*a + m

      if (c <= 0.0_wip) cycle

      d = m*log(c/m)-b*a
      if (d > 80.0_wp) exit
      if (d < -80.0_wp) then
        d = 0.0_wp
      else
        d = (1.0_wip+a*a)*exp(d)
      end if

      random3 = random_get()
      if (random3 <= d) exit

    end do

    gamma_distr = c

  end function gamma_distr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cubic_spline
  !> @brief        natrual cubic spline function
  !! @authors      KY
  !! @param[in]    nn      : size of data
  !! @param[in]    xx(nn)  : x values
  !! @param[in]    yy(nn)  : y values
  !! @param[inout] coeff(4,nn-1) : coefficients of cubic functions
  !! @param[in]    xin     : x value
  !! @param[out]   yout    : y value at xin
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cubic_spline(nn, xx, yy, coeff, xin, yout)

    ! formal arguments
    integer, intent(in)    :: nn
    real(8), intent(in)    :: xx(0:nn-1)
    real(8), intent(in)    :: yy(0:nn-1)
    real(8), intent(inout) :: coeff(0:3,0:nn-2)
    real(8), intent(in)    :: xin
    real(8), intent(out)   :: yout

    ! local arguments
    integer :: i, idx
    real(8) :: dx, dxi


    idx = nn-2
    do i = 1, nn-1
      if(xin < xx(i)) then
        idx = i-1
        exit
      end if
    end do

    dx   = xin - xx(idx)
    dxi  = 1.0_wp
    yout = 0.0_wp
    do i = 0, 3
       yout = yout + coeff(i,idx)*dxi
       dxi  = dxi * dx
    end do

  end subroutine cubic_spline

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cubic_spline_coeff
  !> @brief        calculate the coefficients of cubic functions
  !! @authors      KY
  !! @param[in]    nn      : size of data
  !! @param[in]    xx(nn)  : x values
  !! @param[in]    yy(nn)  : y values
  !! @param[out]   coeff(4,nn-1) : coefficients of cubic functions
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cubic_spline_coeff(nn, xx, yy, coeff)

    ! formal arguments
    integer, intent(in)  :: nn
    real(8), intent(in)  :: xx(0:nn-1)
    real(8), intent(in)  :: yy(0:nn-1)
    real(8), intent(out) :: coeff(0:3,0:nn-2)

    ! local arguments
    integer :: i
    real(8) :: hh(0:nn-2), gg(0:nn-2)
    real(8) :: vii(1:nn-2), aii(1:nn-2), aij(1:nn-3), aji(2:nn-2)
    real(8) :: uii(0:nn-1)


    do i = 0, nn-2
      hh(i) = xx(i+1) - xx(i)
      gg(i) = (yy(i+1) - yy(i))/hh(i)
    end do

    do i = 1, nn-2
      vii(i) = 6.0_wp*(gg(i) - gg(i-1))
      aii(i) = 2.0_wp*(hh(i-1) + hh(i))
    end do

    do i = 1, nn-3
      aij(i)   = hh(i)
      aji(i+1) = hh(i)
    end do

    ! forward
    do i = 1, nn-3
      aij(i) = aij(i) / aii(i)
      vii(i) = vii(i) / aii(i)
      aii(i+1) = aii(i+1) - aij(i)*aji(i+1)
      vii(i+1) = vii(i+1) - vii(i)*aji(i+1)
    end do
    vii(nn-2) = vii(nn-2)/aii(nn-2)

    ! backward
    uii(nn-2) = vii(nn-2)
    do i = nn-3, 1, -1
      uii(i) = vii(i+1) - aij(i)*uii(i+1)
    end do
    uii(0)    = 0.0_wp
    uii(nn-1) = 0.0_wp

    ! calc coefficients
    do i = 0, nn-2
      coeff(0,i) = yy(i)
      coeff(1,i) = (yy(i+1) - yy(i))/hh(i) &
                 - hh(i)*(2.0_wp*uii(i) + uii(i+1))/6.0_wp
      coeff(2,i) = uii(i)*0.5_wp
      coeff(3,i) = (uii(i+1) - uii(i))/hh(i)/6.0_wp
    end do

  end subroutine cubic_spline_coeff

  function lcm(a, b) ! Greatest Common Divisor

    integer :: lcm
    integer, intent(in) :: a, b
    integer :: gcd, num_small, num_big, i

    gcd = 1
    i = 2
    if (a <= b) then
      num_small = a
      num_big   = b
    else
      num_small = b
      num_big   = a
    end if

    do while (i <= num_small)
      if (mod(num_small,i) == 0 .and. mod(num_big,i) == 0) then
        gcd = gcd * i
        num_small = num_small / i
        num_big   = num_big   / i
      else
        i = i + 1
      end if
    end do
    lcm = gcd * num_small * num_big

    end function lcm

end module math_libs_mod
