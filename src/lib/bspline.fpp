!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   bspline_mod
!> @brief   libralies for calculating the b-spline function
!! @authors Junichiro Koga (JK)
!
!  (c) Copyright 2020 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module bspline_mod

  use constants_mod

  implicit none
  private

  public :: bspline
  public :: dbspline
  public :: bspline_single_point
  public :: dbspline_single_point
  public :: bspline_resolve_indices
  public :: bspline_get_total_num_elements
  public :: bspline_get_next_index_offset
  public :: bspline_get_1d_index
  public :: bspline_interpolation
  public :: bspline_solve_for_coeffs

  real(wp), private, parameter     :: onethird  = 1.0_wp/3.0_wp
  real(wp), private, parameter     :: twothirds = 2.0_wp/3.0_wp
  real(wp), private, parameter     :: onesixth  = 1.0_wp/6.0_wp

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bspline
  !> @brief        calculate the b-spline function
  !! @authors      JK
  !! @param[in]    istart : start index of the output array
  !! @param[in]    iend   : end   index of the output array
  !! @param[in]    xin    : the x-value on which the b-spline shall be evaluated
  !! @param[out]   yout   : array of the b-splines
  !! @param[in]    ndim   : the dimension of the b-spline 
  !                         (optinal. either 2 or 3, defaults to 3)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine bspline(istart, iend, xin, yout, ndim)
    integer,  intent(in)           :: istart, iend
    real(wp), intent(in)           :: xin
    real(wp), intent(out)          :: yout(istart:iend)
    integer,  intent(in), optional :: ndim
    real(wp) :: xx, axx
    
    if(present(ndim)) then
      if (ndim==2) then
        call bspline2(istart, iend, xin, yout)
      else
        call bspline3(istart, iend, xin, yout)
      endif
    else
      call bspline3(istart, iend, xin, yout)
    endif
    return
  end subroutine bspline

  subroutine bspline3(istart, iend, xin, yout)
    integer,  intent(in)  :: istart, iend
    real(wp), intent(in)  :: xin
    real(wp), intent(out) :: yout(istart:iend)
    real(wp) :: xx, axx
    integer  :: i
    integer :: is, ie
    is = int(ceiling(xin - 2))
    ie = int(floor  (xin + 2))
    
    ! When single precision is used, is can be equal to istart-1
    ! or ie can be equal to iend+1. At that time, is is forced to
    ! change to istart or ie is forced to change to iend.
    if (is == (istart-1)) is = istart
    if (ie == (iend+1)) ie = iend

    yout = 0.0_wp
    do i = is, ie
      xx = xin - i 
      axx = abs(xx)
      if (axx >= 0 .and. axx <  1) then
        yout(i) = xx * xx * (axx-2) * 0.5_wp + twothirds
      else
        yout(i) = (2 - axx) * (2 - axx) * (2 - axx) * onesixth
      endif
    enddo

    return

  end subroutine bspline3

  subroutine bspline2(istart, iend, xin, yout)
    integer,  intent(in)  :: istart, iend
    real(wp), intent(in)  :: xin
    real(wp), intent(out) :: yout(istart:iend)
    real(wp) :: xx
    integer  :: i
    integer  :: is, ie
    yout = 0.0_wp
   
    is = int(ceiling(xin - 2))
    ie = int(floor  (xin + 1))
    
    ! When single precision is used, is can be equal to istart-1
    ! or ie can be equal to iend+1. At that time, is is forced to
    ! change to istart or ie is forced to change to iend.
    if (is == (istart-1)) is = istart
    if (ie == (iend+1)) ie = iend
    
    do i = is, ie
      xx = xin - i
      if ( xx>=-1 .and. xx<0 ) then
        yout(i) = (xx+1)*(xx+1)*0.5_wp
      else if (xx>=0 .and. xx< 1) then
        yout(i) = -xx*(xx-1) + 0.5_wp
      else if (xx>=1 .and. xx<2) then
        yout(i) = (xx-2)*(xx-2)*0.5_wp
      endif
    enddo

    return

  end subroutine bspline2

  subroutine dbspline(istart, iend, xin, uin, dx, yout)
    integer,  intent(in)  :: istart, iend
    real(wp), intent(in)  :: xin
    real(wp), intent(in)  :: uin (istart:iend)
    real(wp), intent(in)  :: dx
    real(wp), intent(out) :: yout(istart:iend-1)
    real(wp), allocatable :: b2(:)
    integer  :: i
    integer  :: is, ie
    real(wp) :: invdx

!    allocate(b2(istart:iend-1))
    invdx = 1.0_wp/dx
    call bspline(istart, iend-1, xin, yout, 2)
    is = int(ceiling(xin - 2))
    ie = int(floor  (xin + 1))
!    do i = istart, iend-1
    do i = is, ie
      yout(i) = (uin(i+1) - uin(i)) * yout(i) * invdx
    enddo

!    deallocate(b2)
    return
  end subroutine dbspline

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bspline_single_point
  !> @brief        calculate the b-spline function for a single y-value
  !! @authors      JK
  !! @param[in]    i    : the index for which the b-spline shall be evaluated
  !! @param[in]    xin  : the x-value on which the b-spline shall be evaluated
  !! @param[out]   yout : the value of the b-spline
  !! @param[in]    ndim : the dimension of the b-spline 
  !                         (optinal. either 2 or 3, defaults to 3)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine bspline_single_point(i,xin, yout, ndim)
    integer,  intent(in)           :: i
    real(wp), intent(in)           :: xin
    real(wp), intent(out)          :: yout
    integer,  intent(in), optional :: ndim
    real(wp) :: xx, axx
    
    if(present(ndim)) then
      if (ndim==2) then
        call bspline2_single_point(i, xin, yout)
      else
        call bspline3_single_point(i, xin, yout)
      endif
    else
      call bspline3_single_point(i, xin, yout)
    endif
    return
  end subroutine bspline_single_point

  subroutine bspline3_single_point(i, xin, yout)
    integer,  intent(in)  :: i
    real(wp), intent(in)  :: xin
    real(wp), intent(out) :: yout
    real(wp) :: xx, axx
    integer :: is, ie
    yout = 0.0_wp
    xx = xin - i 
    axx = abs(xx)
    if(axx>=2) then
      return
    end if
    if (axx >= 0 .and. axx <  1) then
      yout = xx * xx * (axx-2) * 0.5_wp + twothirds
    else if (axx < 2) then
      yout = (2 - axx) * (2 - axx) * (2 - axx) * onesixth
    endif

    return

  end subroutine bspline3_single_point

  subroutine bspline2_single_point(i, xin, yout)
    integer,  intent(in)  :: i
    real(wp), intent(in)  :: xin
    real(wp), intent(out) :: yout
    real(wp) :: xx
    integer  :: is, ie
    yout = 0.0_wp
   
    xx = xin - i
    if ( xx>=-1 .and. xx<0 ) then
      yout = (xx+1)*(xx+1)*0.5_wp
    else if (xx>=0 .and. xx< 1) then
      yout = -xx*(xx-1) + 0.5_wp
    else if (xx>=1 .and. xx<2) then
      yout = (xx-2)*(xx-2)*0.5_wp
    endif

    return

  end subroutine bspline2_single_point

  subroutine dbspline_single_point(i, istart, iend, xin, uin, dx, yout)
    integer,  intent(in)  :: i, istart, iend
    real(wp), intent(in)  :: xin
    real(wp), intent(in)  :: uin (istart:iend)
    real(wp), intent(in)  :: dx
    real(wp), intent(out) :: yout
    integer  :: is, ie
    real(wp) :: invdx, ytmp

    invdx = 1.0_wp/dx
    call bspline_single_point(i, xin, ytmp, 2)
    yout  = (uin(i+1) - uin(i)) * ytmp * invdx

    return

  end subroutine dbspline_single_point

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bspline_resolve_indices
  !> @brief        obtain an array of indices containing indices for each
  !                dimension from one-dimensionalized index
  !! @authors      JK
  !! @param[in]  i               : the one-dimensional index
  !! @param[in]  ndep            : the dimension of the problem (note: this has nothing 
  !                                to do with the dimension of the b-spline function!)
  !! @param[in]  nn              : the array storing the number of elements for each
  !                                dimension
  !! @param[out] indices         : the array storing the indices.
  !! @param[in], npad (optional) : the amount of padding. defaults to 2
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine bspline_resolve_indices(i, ndep, nn, indices, npad)
    integer, intent(in)           :: i, ndep
    integer, intent(in)           :: nn(ndep)
    integer, intent(out)          :: indices(ndep)
    integer, intent(in), optional :: npad
    integer :: n
    integer :: j, k, nsho, namari
    integer :: ix, iy, iz, iz2, currn, nt, np
    np = 2
    if(present(npad)) np = npad

    if(ndep == 1) then
      indices(1) = i-np
      return
    endif
    currn = i-1
    do j = 1, ndep-1
      nt = 1
      do k = 1, ndep-j
        nt = nt * (nn(k)+np)
      enddo
      nsho   = currn/nt
      namari = mod(currn, nt)
      indices(ndep-j+1) = nsho + 1 - np
      currn = namari
    enddo
    indices(1) = namari + 1 - np

    return

  end subroutine bspline_resolve_indices

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bspline_get_ncombo_excluding
  !> @brief        get the number of combinations excluding the specified
  !                dimension
  !! @authors      JK
  !! @param[in]  i               : the dimension to be ignored
  !! @param[in]  ndim            : the dimension of the problem
  !! @param[in]  nn              : the array storing the number of elements for each
  !                                dimension
  !! @param[in], npad (optional) : the amount of padding. defaults to 2
  !! @return                     : the number of combinations excluding the i-th
  !                                dimension
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  function bspline_get_ncombo_excluding(i, ndim, nn, npad) result(res)
    implicit none
    integer, intent(in)           :: i, ndim
    integer, intent(in), optional :: npad
    integer :: res
    integer :: nn(ndim)
    integer :: j, np
    np = 2
    if(present(npad)) np = npad
    res = 1
    do j = 1, ndim
      if(j==i) cycle
      res = res * (nn(j)+np)
    enddo

    return

  end function bspline_get_ncombo_excluding

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bspline_get_total_num_elements
  !> @brief        get the total number of combinations for the problem
  !! @authors      JK
  !! @param[in]  ndim            : the dimension of the problem
  !! @param[in]  nmesh           : the array storing the number of elements for each
  !                                dimension
  !! @param[in], npad (optional) : the amount of padding. defaults to 2
  !! @return                     : the total number of combinations
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  function bspline_get_total_num_elements(ndim, nmesh, npad) result(res)
    implicit none
    integer, intent(in)           :: ndim
    integer, intent(in)           :: nmesh(ndim)
    integer, intent(in), optional :: npad
    integer :: i, res, np
    np  = 2
    if(present(npad)) np = npad
    res = 1
    do i = 1, ndim
      res = res * (nmesh(i)+np)
    enddo

    return

  end function bspline_get_total_num_elements

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bspline_get_next_index_offset
  !> @brief        get the amount of offset in order to refer to the next
  !                element in multi-dimensional space.
  !! @authors      JK
  !! @param[in]  idimension      : the dimension of interest
  !! @param[in]  ndim            : the dimension of the problem
  !! @param[in]  nn              : the array storing the number of elements for each
  !                                dimension
  !! @param[in], npad (optional) : the amount of padding. defaults to 2
  !! @return                     : the amount of offset needed in order to refer 
  !                                to the next element for the idimension-th 
  !                                dimension
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  function bspline_get_next_index_offset(idimension, ndim, nn, npad) result(res)
    implicit none
    integer, intent(in)           :: idimension, ndim
    integer, intent(in)           :: nn(ndim)
    integer, intent(in), optional :: npad
    integer :: j, np, res
    np = 2
    if(present(npad)) np = npad
    res =  1
    do j=1, idimension-1
      res = res * (nn(j)+np)
    enddo

    return

  end function bspline_get_next_index_offset

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bspline_get_1d_index
  !> @brief        one-dimensionalize the specified array of indices
  !! @authors      JK
  !! @param[in]  n       : the dimension of the problem
  !! @param[in]  indices : the array storing the indices.
  !! @param[in]  nn      : the array storing the number of elements for each
  !                        dimension
  !! @param[in], npad (optional) : the amount of padding. defaults to 2
  !! @return                     : the one-dimensionalized index
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  function bspline_get_1d_index(n, indices, nsizes, periodicity, &
                                npad) result(res)
    integer, intent(in)           :: n
    integer, intent(in)           :: indices(n)
    integer, intent(in)           :: nsizes(n)
    logical, intent(in)           :: periodicity(n)
    integer, intent(in), optional :: npad
    integer :: res
    integer :: i, j, theindex
    integer :: nprod, np
    np = 2
    if(present(npad)) np = npad
    res = 1
    do i = n, 1, -1
      nprod = 1
      do j = 1, i-1
        nprod = nprod * (nsizes(j)+np)
      enddo
      theindex = indices(i)+np-1
      if(.not. periodicity(i)) then
        if ( indices(i) == nsizes(i) )  then
          theindex = theindex - np
        else if ( indices(i)==-1 ) then
          theindex = theindex + np
        end if 
      end if  
      res = res + theindex * nprod
    enddo

    return
  end function bspline_get_1d_index

  subroutine bspline_solve_1d(nout, yout)
    implicit none
    integer, intent(in)          :: nout
    real(kind=wp), intent(inout) :: yout(-1:nout)

    integer,       allocatable   :: ipiv(:)
    real(kind=dp), allocatable   :: yi(:,:)
    real(kind=dp), allocatable   :: smat(:,:)

    integer :: i, j, ista, iend, INFO

    allocate(yi(1:nout,1))
    do i = 0, nout-1
      yi(i+1,1) = yout(i)
    enddo

    allocate(smat(1:nout,1:nout))
    smat = 0.0_dp
    smat(1,1)         = twothirds
    smat(1,2)         = onethird
    smat(nout,nout)   = twothirds
    smat(nout,nout-1) = onethird
    do i = 1, nout-2
      smat(i+1, i  )  = onesixth
      smat(i+1, i+1)  = twothirds
      smat(i+1, i+2)  = onesixth
    enddo

    allocate(ipiv(nout))
    call dgesv(nout, 1, smat, nout, ipiv, yi, nout, INFO)

    do i = 1, nout
      yout(i-1) = yi(i,1)
    enddo
    yout(-1)   = yout(1)
    yout(nout) = yout(nout-2)

    deallocate(yi)
    deallocate(smat)
    deallocate(ipiv)

  end subroutine bspline_solve_1d

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bspline_interpolation
  !> @brief        perform bspline interpolation on to the specified mesh
  !! @authors      JK
  !! @param[in]   ntot_spline  : number of elements for the spline coeffs
  !! @param[in]   ntot_intpl   : number of elements for the interpolating
  !                              function
  !! @param[in]   ndim         : the dimension of the problem
  !! @param[in]   nmesh_spline : number of elements in each dimension for the
  !                              spline coeffs
  !! @param[in]   nmesh_intpl  : number of elements in each dimension for the
  !                              interpolating function
  !! @param[in]   xmin         : array storing the minimum value for each
  !                              dimension
  !! @param[in]   dxs          : array storing the delta value for the spline
  !                              coeffs in each dimension
  !! @param[in]   dxi          : array storing the delta value for the
  !                              interpolating function in each dimension
  !! @param[in]   splinecoeffs : array of spline coefficients
  !! @param[out]  interpolated : array of the interpolating function
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine bspline_interpolation(ntot_spline, ntot_intpl, ndim, &
                                   nmesh_spline, nmesh_intpl,     &
                                   xmin, dxs, dxi,                &
                                   splinecoeffs, interpolated)
    integer,       intent(in)  :: ntot_spline, ntot_intpl, ndim
    integer,       intent(in)  :: nmesh_spline(ndim), nmesh_intpl(ndim)
    real(kind=wp), intent(in)  :: xmin(ndim), dxs(ndim), dxi(ndim)
    real(kind=dp), intent(in)  :: splinecoeffs(ntot_spline)
    real(kind=dp), intent(out) :: interpolated(ntot_intpl)

    integer                    :: i, j, k, ncombo, indcombo, ind
    real(kind=wp)              :: x, ytmp, ytmp2
    integer,       allocatable :: indices(:), nmeshtmp(:), istas(:), iends(:)
    real(kind=wp), allocatable :: xs(:)

    allocate(indices(ndim))
    allocate(nmeshtmp(ndim))
    allocate(istas(ndim))
    allocate(iends(ndim))
    allocate(xs(ndim))
    interpolated = 0.0_wp
    do i = 1, ntot_intpl
      ncombo = 1
      call bspline_resolve_indices(i, ndim, nmesh_intpl, indices)
      do j = 1, ndim
        x           = xmin(j) + dxi(j)*indices(j)
        xs(j)       = (x-xmin(j))/dxs(j)
        istas(j)    = max(int(ceiling(xs(j) - 2)), -1)
        iends(j)    = min(int(floor  (xs(j) + 2)), nmesh_spline(j))
        nmeshtmp(j) = iends(j)-istas(j)+1
        ncombo      = ncombo * nmeshtmp(j)
      enddo
      do j = 1, ncombo
        call bspline_resolve_indices(j, ndim, nmeshtmp, indices, npad=0)
        ytmp = 1.0_wp
        indcombo = 1
        do k = 1, ndim
          ind = istas(k) + indices(k) - 1
          call bspline_single_point(ind, xs(k), ytmp2)
          ytmp = ytmp * ytmp2
          indcombo = indcombo + (ind+1) * bspline_get_next_index_offset( &
                                          k, ndim, nmesh_spline)
        enddo
        interpolated(i) = interpolated(i) + ytmp * splinecoeffs(indcombo)
      enddo
    enddo

    deallocate(indices)
    deallocate(nmeshtmp)
    deallocate(istas)
    deallocate(iends)
    deallocate(xs)

    return

  end subroutine bspline_interpolation

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bspline_solve_for_coeffs
  !> @brief        obtain spline coefficients which can reproduce the specified
  !                function
  !! @authors      JK
  !! @param[in]    ntot             : number of elements
  !! @param[in]    ndim             : dimension of the problem
  !! @param[in]    nmesh_spline_new : the number of elements in each dimension
  !! @param[inout] newsplinecoeffs  : on input, the function to be resolved. on
  !                                   output, the spline coefficients
  !! @param[in], optional positive  : specify whether the spline-coefficients
  !                                   shall be always positive or not.
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine bspline_solve_for_coeffs(ntot, ndim, nmesh_spline_new, &
                                      newsplinecoeffs, positive)
    integer,           intent(in)    :: ntot, ndim
    integer,           intent(in)    :: nmesh_spline_new(ndim)
    real(kind=dp),     intent(inout) :: newsplinecoeffs(ntot)
    logical, optional, intent(in) :: positive

    integer :: i, j, k, ncombo, noff, noffnext, amari, sho, ind
    integer,       allocatable :: indices(:)
    real(kind=wp), allocatable :: spltmp(:)

    allocate(indices(ndim))
    do i = 1, ndim
      ncombo = bspline_get_ncombo_excluding (i, ndim, nmesh_spline_new)
      noff   = bspline_get_next_index_offset(i, ndim, nmesh_spline_new)
      if(i == ndim) then
        noffnext = 0
      else
        noffnext = bspline_get_next_index_offset(i+1, ndim, nmesh_spline_new)
      endif
      allocate(spltmp(-1:nmesh_spline_new(i)))
      do j = 1, ncombo
        amari  = mod(j-1, noff)
        sho    = (j-1)/noff
        sho    = sho * noffnext
        spltmp = 0.0_wp
        do k = -1, nmesh_spline_new(i)
          ind       = (k+1)*noff + 1 + sho + amari
          spltmp(k) = newsplinecoeffs(ind)
        enddo
        call bspline_solve_1d(nmesh_spline_new(i),spltmp)
        do k = -1, nmesh_spline_new(i)
          ind                  = (k+1)*noff + 1 + sho + amari
          newsplinecoeffs(ind) = spltmp(k)
        enddo
      enddo
      deallocate(spltmp)
    enddo
    deallocate(indices)

    if(present(positive)) then
      if(positive) then
        ncombo = bspline_get_total_num_elements(ndim, nmesh_spline_new)
        do j = 1, ncombo
          newsplinecoeffs(j) = abs(newsplinecoeffs(j))
        enddo
      endif
    endif

    return
  
  end subroutine bspline_solve_for_coeffs

end module bspline_mod

