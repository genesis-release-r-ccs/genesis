!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   dihedral_libs_mod
!> @brief   libralies for dihedral angles (CMAP and proper)
!! @authors Yasuhiro Matsunaga (YM), Chigusa Kobayashi (CK), Takao Yoda (TY) 
!! @date    2010/12/27-2011/2/3 (TY)
!! @date    2014/4/16-2014/4/20 (YM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module dihedral_libs_mod

  use fileio_par_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  public  :: calculate_dihedral
  public  :: derive_cmap_coefficients_p   ! derive cmap parameters
                                          ! to calculate c_ij (periodic)
  public  :: derive_cmap_coefficients_np  ! derive cmap parameters
                                          ! to calculate c_ij (natural)

  private :: prep_splin_matrix_p ! (for periodic spline)
  private :: prep_splin_vector_p ! (for periodic spline)
  private :: simultaneous_eq     ! solve simultaneous linear equations
  private :: cal_y1              ! calc dEcmap/d(dih) by spline
  private :: cmap_cal_coefs      ! evaluate c_ij by bicubic interpolation
  private :: bicubic_coefs       ! calc coefs for bicubic interpolation 
  private :: solve_nspline       ! evaluate d2y/dx2 by natural cubic spline
  private :: cal_ecmap_only      ! for chck Ecmap
  private :: cal_decmap_ddih     ! for chck dEcmap/d(dih)
  private :: cal_d2ecmap         ! for chck d2Ecmap/d(dih)2

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calculate_dihedral
  !> @brief        calculate dihedral angle
  !! @authors      CK
  !! @param[in]    aindex  : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] cos_dih : cosin of dihedral angles
  !! @param[inout] sin_dih : sin of dihedral angles
  !! @param[inout] grad    : gradient of dihedral angles
  !! @note         Blondel and Karplus, J. Comput. Chem., 17, 1132-1141 (1996)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calculate_dihedral(aindex, coord, cos_dih, sin_dih, grad, v)

    ! formal arguments
    integer,   intent(in)    :: aindex(:)
    real(wp),  intent(in)    :: coord(:,:)
    real(wp),  intent(inout) :: cos_dih
    real(wp),  intent(inout) :: sin_dih
    real(wp),  intent(inout) :: grad(1:9)
    real(wp),  intent(inout) :: v(1:3,1:3)

    ! local variable
    integer                  :: i, j
    real(wp)                 :: dij(1:3), djk(1:3), dlk(1:3)
    real(wp)                 :: aijk(1:3), ajkl(1:3)
    real(wp)                 :: tmp(1:4)
    real(wp)                 :: raijk2, rajkl2, vtmp
    real(wp)                 :: inv_raijk2, inv_rajkl2, inv_raijkl
    real(wp)                 :: rjk, inv_rjk, dotpro_ijk, dotpro_jkl
      
    dij(1:3) = coord(1:3,aindex(1)) - coord(1:3,aindex(2))
    djk(1:3) = coord(1:3,aindex(2)) - coord(1:3,aindex(3))
    dlk(1:3) = coord(1:3,aindex(4)) - coord(1:3,aindex(3))

    aijk(1) = dij(2)*djk(3) - dij(3)*djk(2)
    aijk(2) = dij(3)*djk(1) - dij(1)*djk(3)
    aijk(3) = dij(1)*djk(2) - dij(2)*djk(1)

    ajkl(1) = dlk(2)*djk(3) - dlk(3)*djk(2)
    ajkl(2) = dlk(3)*djk(1) - dlk(1)*djk(3)
    ajkl(3) = dlk(1)*djk(2) - dlk(2)*djk(1)

    raijk2     = aijk(1)*aijk(1) + aijk(2)*aijk(2) + aijk(3)*aijk(3)
    rajkl2     = ajkl(1)*ajkl(1) + ajkl(2)*ajkl(2) + ajkl(3)*ajkl(3)

    inv_raijk2  = 1.0_wp / raijk2
    inv_rajkl2  = 1.0_wp / rajkl2

    inv_raijkl = sqrt(inv_raijk2*inv_rajkl2)

    cos_dih = (aijk(1)*ajkl(1) + aijk(2)*ajkl(2) + aijk(3)*ajkl(3))*inv_raijkl

    rjk     = sqrt(djk(1)*djk(1) + djk(2)*djk(2) + djk(3)*djk(3))
    inv_rjk = 1.0_wp/rjk

    tmp(1)  = aijk(1)*dlk(1) + aijk(2)*dlk(2) + aijk(3)*dlk(3)
    sin_dih = tmp(1) * rjk * inv_raijkl

    dotpro_ijk = dij(1)*djk(1) + dij(2)*djk(2) + dij(3)*djk(3)
    dotpro_jkl = djk(1)*dlk(1) + djk(2)*dlk(2) + djk(3)*dlk(3)

    tmp(1) = rjk*inv_raijk2
    tmp(2) = rjk*inv_rajkl2

    tmp(3) =  dotpro_ijk*inv_raijk2*inv_rjk
    tmp(4) =  dotpro_jkl*inv_rajkl2*inv_rjk

    grad(1:3) =  tmp(1)*aijk(1:3)
    grad(4:6) = -tmp(3)*aijk(1:3) + tmp(4)*ajkl(1:3)
    grad(7:9) =                   - tmp(2)*ajkl(1:3)

    do i = 1, 3
      do j = i+1, 3
      vtmp = - (grad(j)*dij(i) + grad(j+3)*djk(i) + grad(j+6)*dlk(i))
      v(j, i) = vtmp
      v(i, j) = vtmp
      end do
      vtmp = - (grad(i)*dij(i) + grad(i+3)*djk(i) + grad(i+6)*dlk(i))
      v(i, i) = vtmp
    end do

    return
  end subroutine calculate_dihedral

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    derive_cmap_coefficients_p
  !> @brief        evaluate coefficients (c_ij) values for cmap term (periodic)
  !! @authors      TY
  !! @param[in]    cmap_type_id : 1 <= cmap_type_id <= par%num_cmaps = 6
  !! @param[in]    par          : CHARMM par information
  !! @param[out]   c_ij         : cmap coefficients for each type of cmap
  !! @date         2010/12/27-2011/2/3 (TY)
  !
  ! (note on 2010/12/27 by TY)
  !      This routne calculates c_ij values by use of bicubic interpolation.
  !      This version uses 1-d cubic spline with PERIODIC boundary condition,
  !      instead of the "NATURAL" spline that seems to be used in the original
  !      CHARMM program.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine derive_cmap_coefficients_p(cmap_type_id, par, c_ij )
   
    !  formal arguments
    integer,      intent(in)    :: cmap_type_id
    type(s_par),  intent(in)    :: par
    real(wp),     intent(inout) :: c_ij(:,:,:,:)
   
    !  local variables
    real(wp)      :: dlt_x    ! distance between grid lines (degree)
    real(wp)      :: rdum,rdum1,rdum2,rdum3
    integer       :: ngrid0, ngrid
    integer       :: alloc_stat, dealloc_stat, i, j,  idih1, idih2

    real(wp),     allocatable :: y_data(:,:)    ! extended on-grid Ecmap tab
    real(wp),     allocatable :: deriv_1(:,:)   ! on-grid dEcmap/d(dih1) 
    real(wp),     allocatable :: deriv_2(:,:)   ! on-grid dEcmap/d(dih2) 
    real(wp),     allocatable :: deriv_12(:,:)  ! on-grid d2Ecmap/d(dih1)d(dih2)
    real(wp),     allocatable :: deriv_tmp(:,:) ! on-grid d2Ecmap/d(dih2)d(dih1)
    real(wp),     allocatable :: rmatrix(:,:)   ! matrix for simultaneous eq
    real(wp),     allocatable :: rvector(:)     ! vector for simultaneous eq
    real(wp),     allocatable :: y_1d_vec(:)    


    alloc_stat   = 0
    dealloc_stat = 0

    ! preparation for spline
    !

    ! (size of the cmap table)
    !
    ngrid0 = par%cmap_resolution(cmap_type_id)
    ngrid  = ngrid0 + 1
    dlt_x  = 360.0_wp/real(ngrid0,wp)

    ! (allocate an array for the extended on-grid-points Ecmap table)
    !
    allocate(y_data(ngrid,ngrid), stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    ! (copy and extend the on-grid-points Ecmap table)
    !
    y_data(1:ngrid0, 1:ngrid0) = par%cmap_data(1:ngrid0, 1:ngrid0, cmap_type_id)
    y_data(  ngrid , 1:ngrid0) = y_data(    1   , 1:ngrid0)
    y_data(1:ngrid ,   ngrid ) = y_data(1:ngrid ,     1   )

    ! (allocate the other local arrays )
    !
    allocate(rmatrix  (ngrid, ngrid), &
             deriv_1  (ngrid, ngrid), &
             deriv_2  (ngrid, ngrid), &
             deriv_12 (ngrid,ngrid),  &
             deriv_tmp(ngrid,ngrid),  &
             rvector  (ngrid),        &
             y_1d_vec (ngrid), stat = alloc_stat)

    if (alloc_stat /= 0) then
      write(MsgOut,*) 'Derive_Cmap_Coefficients_P> ERROR: allocation error in derive_cmap_coefficients_p.'
      call error_msg_alloc
    end if

    ! (prepare the matrix for the simultaneous equations)
    !
    call prep_splin_matrix_p(ngrid, rmatrix)

    ! calculate dy/d(dih1) on grid points by cubic spline
    !
    do idih2 = 1, ngrid
      y_1d_vec(1:ngrid) = y_data(1:ngrid, idih2)
      call prep_splin_vector_p(ngrid, dlt_x, y_1d_vec, rvector)

      ! (solve simultaneous equations for spline.  rvector will be d2y/d(dih1)2)
      !
      call simultaneous_eq(ngrid, rmatrix, rvector)

      ! (calculate dy/d(dih1) using y and precalculated d2y/d(dih1)2 values)
      !
      call cal_y1(ngrid, dlt_x, y_1d_vec, rvector, 0)

      ! (set the resultant dy/d(dih1) values to array)
      !
      deriv_1(1:ngrid, idih2) = rvector(1:ngrid)
    end do

    ! calculate dy/d(dih2) on grid points by cubic spline
    !
    do idih1 = 1, ngrid
      y_1d_vec(1:ngrid) = y_data(idih1, 1:ngrid)
      call prep_splin_vector_p(ngrid, dlt_x, y_1d_vec, rvector)
      call simultaneous_eq(ngrid, rmatrix, rvector)
      call cal_y1(ngrid, dlt_x, y_1d_vec, rvector, 0)
      deriv_2(idih1, 1:ngrid) = rvector(1:ngrid)
    end do

    ! calculate d2y/d(dih1)d(dih2) on grid points by cubic spline
    !
    do idih1 = 1, ngrid
      y_1d_vec(1:ngrid) = deriv_1(idih1, 1:ngrid)
      call prep_splin_vector_p(ngrid, dlt_x, y_1d_vec, rvector)
      call simultaneous_eq(ngrid, rmatrix, rvector)
      call cal_y1(ngrid, dlt_x, y_1d_vec, rvector, 0)
      deriv_12(idih1, 1:ngrid) = rvector(1:ngrid)
    end do

    ! calculate d2y/d(dih2)d(dih1) on grid points by cubic spline
    !
    do idih2 = 1, ngrid
      y_1d_vec(1:ngrid) = deriv_2(1:ngrid, idih2)
      call prep_splin_vector_p(ngrid, dlt_x, y_1d_vec, rvector)
      call simultaneous_eq(ngrid, rmatrix, rvector)
      call cal_y1(ngrid, dlt_x, y_1d_vec, rvector, 0)
      deriv_tmp(1:ngrid, idih2) = rvector(1:ngrid)
    end do
   
    ! examine the identity of cross derivative values 
    !
    do i = 1, ngrid
      do j = 1, ngrid
        rdum1 = deriv_12(i,j)
        rdum2 = deriv_tmp(i,j)
        rdum  = 0.5_wp*(rdum1 + rdum2)

        deriv_12(i,j) = rdum

        if (rdum /= 0.0_wp) then
          rdum3 = abs((rdum1 - rdum2)/rdum)
        else
          rdum3 = 0.0_wp
        end if

#ifdef DEBUG
        if (rdum3 > 0.00010_wp) then
          write(MsgOut,*) 'Derive_Cmap_Coefficients_P> WARN: d2y/dphidpsi seems to be inaccurate.'
          write(MsgOut,*) i,j,rdum1, rdum2
          write(MsgOut,*) ''
        end if
#endif
      end do
    end do

    ! preparation for bicubic interpolation has been finished
    !

    ! execute bicubic interpolation
    !
    call cmap_cal_coefs(ngrid0, dlt_x, y_data, deriv_1, deriv_2, deriv_12, c_ij)
   
    ! deallocate local arrays
    !
    deallocate(y_data,   &
               rmatrix,  &
               rvector,  &
               y_1d_vec, &
               deriv_1,  &
               deriv_2,  &
               deriv_12, &
               deriv_tmp, stat = dealloc_stat)
    if (dealloc_stat /= 0) &
      call error_msg_dealloc

    return

  end subroutine derive_cmap_coefficients_p

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    derive_cmap_coefficients_np
  !> @brief        evaluate coefficients (c_ij)values for cmap term(no-periodic)
  !! @authors      TY
  !! @param[in]    cmap_type_id : 1 <= cmap_type_id <= par%num_cmaps = 6
  !! @param[in]    par          : CHARMM par information
  !! @param[out]   c_ij         : cmap coefficients for each type of cmap
  !! @date         2011/2/18- (TY)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine derive_cmap_coefficients_np(cmap_type_id, par, c_ij)

    ! formal arguments
    integer,      intent(in)    :: cmap_type_id
    type(s_par),  intent(in)    :: par
    real(wp),     intent(inout) :: c_ij(:,:,:,:)

    ! local variables
    real(wp)      :: dlt_x    ! distance between grid lines (degree)
    real(wp)      :: rdum,rdum1,rdum2,rdum3
    integer       :: ngrid0, ngrid
    integer       :: alloc_stat, dealloc_stat, i, j,  idih1, idih2, i0, j0
    integer       :: XM

    real(wp),     allocatable :: y_data(:,:)    ! extended on-grid Ecmap tab
    real(wp),     allocatable :: deriv_1(:,:)   ! on-grid dEcmap/d(dih1) 
    real(wp),     allocatable :: deriv_2(:,:)   ! on-grid dEcmap/d(dih2) 
    real(wp),     allocatable :: deriv_12(:,:)  ! on-grid d2Ecmap/d(dih1)d(dih2)
    real(wp),     allocatable :: deriv_tmp(:,:) ! on-grid d2Ecmap/d(dih2)d(dih1)
    real(wp),     allocatable :: rmatrix(:,:)   ! matrix for simultaneous eq
    real(wp),     allocatable :: rvector(:)     ! vector for simultaneous eq
    real(wp),     allocatable :: y_1d_vec(:)

    real(wp),     allocatable :: y_data0(:,:)
    real(wp),     allocatable :: deriv_10(:,:)
    real(wp),     allocatable :: deriv_20(:,:)
    real(wp),     allocatable :: deriv_120(:,:)


    alloc_stat   = 0
    dealloc_stat = 0

    ! preparation for spline
    !

    ! size of the cmap table
    ngrid0 = par%cmap_resolution(cmap_type_id)

    if (mod(ngrid0,2) == 0) then
      XM = ngrid0/2
    else
      XM = ngrid0/2 + 1
    end if

    ngrid  = ngrid0 + 2*XM
    dlt_x  = 360.0_wp/real(ngrid0,wp)

    ! allocate an array for the extended on-grid-points Ecmap table
    !
    allocate(y_data(ngrid,ngrid), stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    ! copy and extend the on-grid-points Ecmap table
    !
    do i = 1, ngrid
      i0 = mod(i+ngrid0-XM-1, ngrid0) + 1
      do j = 1, ngrid
        j0 = mod(j+ngrid0-XM-1, ngrid0) + 1
        y_data(i,j) = par%cmap_data(i0, j0, cmap_type_id)
      end do
    end do

    ! allocate the other local arrays
    !
    allocate(rmatrix  (ngrid, ngrid), &
             deriv_1  (ngrid, ngrid), &
             deriv_2  (ngrid, ngrid), &
             deriv_12 (ngrid,ngrid),  &
             deriv_tmp(ngrid,ngrid),  &
             rvector  (ngrid),        &
             y_1d_vec (ngrid), stat = alloc_stat)
    if (alloc_stat /= 0) then
      write(MsgOut,*) 'Derive_Cmap_Coefficients_Np> ERROR: allocation error.'
      call error_msg_alloc
    end if

    ! calculate dy/d(dih1) on grid points by cubic spline
    !
    do idih2 = 1, ngrid
      y_1d_vec(1:ngrid) = y_data(1:ngrid, idih2)

      ! solve simultaneous equations for spline.  rvector will be d2y/d(dih1)2
      !
      call solve_nspline(dlt_x, y_1d_vec, ngrid, rvector)

      ! calculate dy/d(dih1) using y and precalculated d2y/d(dih1)2 values
      !
      call cal_y1(ngrid, dlt_x, y_1d_vec, rvector, 1)

      ! set the resultant dy/d(dih1) values to array
      !
      deriv_1(1:ngrid, idih2) = rvector(1:ngrid)
    end do

    ! calculate dy/d(dih2) on grid points by cubic spline
    !
    do idih1 = 1, ngrid
      y_1d_vec(1:ngrid) = y_data(idih1, 1:ngrid)
      call solve_nspline(dlt_x, y_1d_vec, ngrid, rvector)
      call cal_y1(ngrid, dlt_x, y_1d_vec, rvector, 1)
      deriv_2(idih1, 1:ngrid) = rvector(1:ngrid)
    end do

    ! calculate d2y/d(dih1)d(dih2) on grid points by cubic spline
    !
    do idih1 = 1, ngrid
      y_1d_vec(1:ngrid) = deriv_1(idih1, 1:ngrid)
      call solve_nspline(dlt_x, y_1d_vec, ngrid, rvector)
      call cal_y1(ngrid, dlt_x, y_1d_vec, rvector, 1)
      deriv_12(idih1, 1:ngrid) = rvector(1:ngrid)
    end do

    ! calculate d2y/d(dih2)d(dih1) on grid points by cubic spline
    !
    do idih2 = 1, ngrid
      y_1d_vec(1:ngrid) = deriv_2(1:ngrid, idih2)
      call solve_nspline(dlt_x, y_1d_vec, ngrid, rvector)
      call cal_y1(ngrid, dlt_x, y_1d_vec, rvector, 1)
      deriv_tmp(1:ngrid, idih2) = rvector(1:ngrid)
    end do

    ! examine the identity of cross derivative values 
    !
    do i = 1, ngrid
      do j = 1, ngrid
        rdum1 = deriv_12(i,j)
        rdum2 = deriv_tmp(i,j)
        rdum  = 0.5_wp*(rdum1 + rdum2)

        deriv_12(i,j) = rdum

        if (rdum /= 0.0_wp) then
          rdum3 = abs((rdum1 - rdum2)/rdum)
        else
          rdum3 = 0.0_wp
        end if

#ifdef DEBUG
        if (rdum3 > EPS_CMAP) then
          write(MsgOut,*) 'Derive_Cmap_Coefficients_Np> WARN: d2y/dphidpsi seems to be inaccurate.'
          write(MsgOut,*) i,j,rdum1, rdum2
          write(MsgOut,*) ''
        end if
#endif
      end do
    end do

    !  copy back on-grid-point y and derivatives values to
    !  (1+ngrid0)x(1+ngrid0) arrays
    !
    ngrid = ngrid0 + 1

    allocate(deriv_10 (ngrid, ngrid), &
             deriv_20 (ngrid, ngrid), &
             deriv_120(ngrid,ngrid),  &
             y_data0  (ngrid,ngrid), stat = alloc_stat)
    if (alloc_stat /= 0) then
      write(MsgOut,*) 'Derive_Cmap_Coefficients_Np> ERROR: allocation error.'
      call error_msg_alloc
    end if

    y_data0  (1:ngrid, 1:ngrid) = y_data  (XM+1:XM+ngrid, XM+1:XM+ngrid)
    deriv_10 (1:ngrid, 1:ngrid) = deriv_1 (XM+1:XM+ngrid, XM+1:XM+ngrid)
    deriv_20 (1:ngrid, 1:ngrid) = deriv_2 (XM+1:XM+ngrid, XM+1:XM+ngrid)
    deriv_120(1:ngrid, 1:ngrid) = deriv_12(XM+1:XM+ngrid, XM+1:XM+ngrid)

    ! preparation for bicubic interpolation has been finished
    !

    ! execute bicubic interpolation
    !
    call cmap_cal_coefs(ngrid0, dlt_x, y_data0, &
                        deriv_10, deriv_20, deriv_120, c_ij)

    ! deallocate local arrays
    !
    deallocate(y_data,    &
               rmatrix,   &
               rvector,   &
               y_1d_vec,  &
               deriv_1,   &
               deriv_2,   &
               deriv_12,  &
               deriv_tmp, &
               y_data0,   &
               deriv_10,  &
               deriv_20,  &
               deriv_120, stat = dealloc_stat)
    if (dealloc_stat /= 0) &
      call error_msg_dealloc

    return

  end subroutine derive_cmap_coefficients_np

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    prep_splin_matrix_p
  !> @brief        prepare the matrix that will be used for the simultaneous
  !!               equations for 1-d cubic spline interpolation with PERIODIC 
  !!               boundary condition.
  !! @authors      TY
  !! @param[in]    ngrid   : the number of equations to be solved simultaneously
  !! @param[out]   rmatrix : matrix for the simultaneous equations
  !! @date         2010/12/27 (TY)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine prep_splin_matrix_p(ngrid,  rmatrix)

    ! formal arguments
    integer,                 intent(in)    :: ngrid
    real(wp),                intent(out)   :: rmatrix(:,:)
   
    ! local variables
    integer                  :: i,  nerr


    ! examine the size of the matrices
    !    

    nerr = 0
    if (size(rmatrix, 1) /= ngrid) nerr = nerr + 1
    if (size(rmatrix, 2) /= ngrid) nerr = nerr + 1
    if (nerr > 0) then
      write(MsgOut) " the number of disagreements = ", nerr
      call error_msg('Enefunc> matrix size looks wrong.') 
    end if

    ! setup the matrix
    !
    rmatrix(:,:)           = 0.0_wp
    rmatrix(1,1)           = 4.0_wp
    rmatrix(1,2)           = 1.0_wp
    rmatrix(1,ngrid-1)     = 1.0_wp
    rmatrix(ngrid,2)       = 1.0_wp
    rmatrix(ngrid,ngrid-1) = 1.0_wp
    rmatrix(ngrid,ngrid  ) = 4.0_wp
   
    do i = 2, ngrid - 1
      rmatrix(i,i-1) = 1.0_wp
      rmatrix(i,i  ) = 4.0_wp
      rmatrix(i,i+1) = 1.0_wp
    end do

    return

  end subroutine prep_splin_matrix_p

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    prep_splin_vector_p
  !> @brief        prepare the vector that will be used for the simultaneous
  !!               equations for 1-d cubic spline interpolation with PERIODIC 
  !!               boundary condition.
  !! @authors      TY
  !! @param[in]    dlt_x    : distance between adjacent grid points
  !! @param[in]    ngrid    : the number of equations to be solved 
  !!                          simultaneously
  !! @param[in]    y_1d_vec : array carrying "y" values on grid points
  !! @param[out]   rvector  : vector for the simultaneous equations
  !! @date         2010/12/27 (TY)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine prep_splin_vector_p(ngrid, dlt_x, y_1d_vec, rvector)

    ! formal arguments
    integer,                 intent(in)    :: ngrid
    real(wp),                intent(in)    :: dlt_x
    real(wp),                intent(in)    :: y_1d_vec(:)
    real(wp),                intent(out)   :: rvector(:)

    ! local variables
    real(wp)                 :: rcoe
    integer                  :: i


    ! examine the size of vectors
    !
    if (size(y_1d_vec) /= ngrid) &
      call error_msg('Enefunc> size of y_1d_vec looks wrong.')

    if (size(rvector) /= ngrid) &
      call error_msg('Enefunc> size of rvector looks wrong.')


    ! setup the vector
    !
    rcoe = dlt_x**2
    rcoe = 6.0_wp/rcoe

    rvector(1)     = y_1d_vec(2)-y_1d_vec(1)+y_1d_vec(ngrid - 1)-y_1d_vec(ngrid)
    rvector(ngrid) = rvector(1)

    do  i = 2, ngrid - 1
      rvector(i) = y_1d_vec(i-1) - 2.0_wp*y_1d_vec(i) + y_1d_vec(i+1)
    end do

    rvector(1:ngrid) = rcoe * rvector(1:ngrid)

    return

  end subroutine prep_splin_vector_p

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    simultaneous_eq
  !> @brief        solve linear simultaneous equations for spline
  !! @authors      TY
  !! @param[in]    ngrid   : the number of equations to be solved simultaneously
  !! @param[in]    rmatrix : matrix for the simultaneous equations
  !! @param[inout] rvector : (i) vector for the simultaneous equations
  !!                         (o) solution of the equations
  !!
  !! @date       2010/12/29 (TY)
  !
  !  note by TY
  !
  !       This subroutine requires LAPACK and BLOS.
  !       A driver routine (dgesv.f) is used.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine simultaneous_eq(ngrid, rmatrix, rvector)

    ! local constant
    integer,                 parameter     :: MAXGRID = 48

    ! formal arguments
    integer,                 intent(in)    :: ngrid
    real(wp),                intent(in)    :: rmatrix(:,:)
    real(wp),                intent(inout) :: rvector(:)

    ! local variables
    real(wp)                 :: rmatrix2(MAXGRID,MAXGRID)
    real(wp)                 :: rvector2(MAXGRID)
    integer                  :: ipiv(MAXGRID) , info


    ! examine the size of arrays
    !
    if (ngrid > MAXGRID) then
      write(MsgOut,*) ''
      write(MsgOut,*) '  ngrid   =', ngrid
      write(MsgOut,*) '  maxgrid =', MAXGRID
      call error_msg('Setup_Enefunc> maxgrid is too small.')       
    end if

    ! copy the matrix and the vector to local arrays
    !
    rmatrix2(1:ngrid, 1:ngrid) = rmatrix(1:ngrid, 1:ngrid)
    rvector2(1:ngrid)          = rvector(1:ngrid)

    ! call lapack routine  
    !

#ifdef LAPACK

    call dgesv(ngrid, 1, rmatrix2, MAXGRID, ipiv, rvector2, MAXGRID, info)

    if (info > 0) then

      write(ErrOut,*) 'Simultaneous_Eq> ERROR. '
      write(ErrOut,*)  '       i =', info
      write(ErrOut,*) ' ---LAPACK says----'
      write(ErrOut,*)' U(i,i) is exactly zero.  The factorization'
      write(ErrOut,*)' has been completed, but the factor U is exactly'
      write(ErrOut,*)' singular, so the solution could not be computed.'
      call error_msg(' ----stop.----')

    elseif (info < 0) then
      write(ErrOut,*) 'Simultaneous_Eq> ERROR.'
      write(ErrOut,*) '       i = ', info, '* -1'
      write(ErrOut,*) ' ----LAPACK says----'
      write(ErrOut,*) 'The i-th argument had an illegal value'
      call error_msg(' ----stop.----')

    end if

#else
    call error_msg('Simultaneous_Eq> ERROR: This subroutine needs LAPACK.')
#endif

    ! copy back the solution to the array for output
    !
    rvector(1:ngrid) = rvector2(1:ngrid)

    return

  end subroutine simultaneous_eq

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cal_y1
  !> @brief        calculate dy/dx values by cubic spline,based on y and d2y/dx2
  !! @authors      TY
  !! @param[in]    ngrid    : the number of grids in one direction
  !! @param[in]    dlt_x    : distance between adjacent grid points
  !! @param[in]    y_1d_vec : on-grid-points y values
  !! @param[inout] rvector  : (i) on-grid-points d2y/dx2 values
  !!                          (o) on-grid-points dy/dx values
  !! @param[in]    mode     : mode for periodic boundary check, 
  !!                          mode==0: do it, mode==1: do not it
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cal_y1(ngrid, dlt_x, y_1d_vec, rvector, mode)

    ! formal arguments
    integer,                 intent(in)    :: ngrid
    real(wp),                intent(in)    :: dlt_x
    real(wp),                intent(in)    :: y_1d_vec(:)
    real(wp),                intent(inout) :: rvector(:)
    integer,                 intent(in)    :: mode

    ! local variables
    real(wp)                 :: rdum, rdum1, rdum2, rdum3, fact1, fact2
    integer                  :: i, alloc_stat, dealloc_stat

    real(wp),   allocatable  :: rrtmp(:)


    alloc_stat   = 0
    dealloc_stat = 0

    allocate(rrtmp(ngrid), stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    fact1 = 1.0_wp/dlt_x
    fact2 = dlt_x/6.0_wp

    ! evaluate dy/dx for not boundary points
    !
    do i = 2, ngrid - 1

      ! estimate dy/dx by two means
      rdum1 =   (y_1d_vec(i) - y_1d_vec(i-1))*fact1 &
              + (rvector(i-1) + 2.0_wp*rvector(i))*fact2

      rdum2 =   (y_1d_vec(i+1) - y_1d_vec(i))*fact1 &
              - (2.0_wp*rvector(i) + rvector(i+1))*fact2
      rdum  = 0.5_wp * (rdum1 + rdum2)
      rrtmp(i) = rdum

    end do

    ! evaluate dy/dx on boundary grid points
    !

    ! ( estimate dy/dx by two means )
    !
    rdum1 = (y_1d_vec(2)-y_1d_vec(1))*fact1 &
          - (2.0_wp*rvector(1)+rvector(2))*fact2

    rdum2 = (y_1d_vec(ngrid)-y_1d_vec(ngrid - 1))*fact1 &
          + (rvector(ngrid - 1)+2.0_wp*rvector(ngrid))*fact2

    ! for periodic boundary condition
    !
    if (mode == 0) then

      rdum = 0.5_wp*(rdum1 + rdum2)
      rrtmp(1)     = rdum
      rrtmp(ngrid) = rdum

      ! examine error
      !
      rdum3 = rdum1 - rdum2
      rdum3 = abs(rdum3)

      if (rdum1*rdum2 > 0) then

        rdum3 = rdum3 / rdum
        if (rdum3 > 0.00000010_wp) then
          if (main_rank) then
            write(MsgOut,'(a)') &
                 'Cal_Y1> INFO: disagreement in dy/dx values at boundary(1).'
            write(MsgOut,*) &
                 '    ',rdum1,rdum2,rdum3
          end if
        end if

      else if (rdum3 > 2.0e-17_wp) then
        if (main_rank) then
          write(MsgOut,'(a)') &
               'Cal_Y1> INFO: disagreement in dy/dx values at boundary(2).'
          write(MsgOut,*) &
               '    ', rdum1, rdum2, rdum3
        end if
      end if
      ! (examination over)

    else

      ! for not periodic boundaries, skip the examination
      !
      rrtmp(1)     = rdum1
      rrtmp(ngrid) = rdum2

    endif

    ! set the result back to rvector for output
    !
    rvector(1:ngrid) = rrtmp(1:ngrid)

    ! deallocation
    !
    deallocate(rrtmp, stat = dealloc_stat)
    if (dealloc_stat /= 0) &
      call error_msg_dealloc

    return

  end subroutine cal_y1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    camp_cal_coefs
  !> @brief        calculate coefficients c_ij for cmap energy term
  !! @authors      TY, YM
  !! @param[in]    ngrid0   : the number of cells in one direction
  !! @param[in]    dlt_x    : distance between adjacent grid points in degree
  !! @param[in]    y_data   : extended on-grid Ecmap tab
  !! @param[in]    deriv_1  : on-grid dEcmap/d(dih1) 
  !! @param[in]    deriv_2  : on-grid dEcmap/d(dih2) 
  !! @param[in]    deriv_12 : on-grid d2Ecmap/d(dih1)d(dih2)
  !! @param[inout] c_ij     : cmap coefficients 
  !!
  !! @note
  !!   c_ij(i,j,icell1,icell2)  : coefficients for the unit cell that has 
  !!                              four grid points:
  !!    where 1 <= icell? <= 24 (= ngrid0)
  !!
  !!     c_ij(i,j,icell1,icell2), where 1 <= i,j <= 4 represent
  !!     coefficients for the cell that has four edge points:
  !!
  !!            (dih1(icell1),     dih2(icell2)    ) <---- [1],
  !!            (dih1(icell1 + 1),     dih2(icell2)) <---- [2],
  !!            (dih1(icell1 + 1), dih2(icell2 + 1)) <---- [3], and
  !!            (dih1(icell1),     dih2(icell2 + 1)) <---- [4].
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cmap_cal_coefs(ngrid0, dlt_x, y_data, deriv_1, deriv_2, &
                            deriv_12, c_ij)

    ! formal arguments
    integer,      intent(in)    :: ngrid0
    real(wp),     intent(in)    :: dlt_x
    real(wp),     intent(in)    :: y_data(:,:)  ! extended on-grid Ecmap tab
    real(wp),     intent(in)    :: deriv_1(:,:) ! on-grid dEcmap/d(dih1) 
    real(wp),     intent(in)    :: deriv_2(:,:) ! on-grid dEcmap/d(dih2) 
    real(wp),     intent(in)    :: deriv_12(:,:)! on-grid d2Ecmap/d(dih1)d(dih2)
    real(wp),     intent(inout) :: c_ij(:,:,:,:)! cmap coefficients

    ! local variables
    real(wp)      :: c(1:4,1:4)
    integer       :: idih1, idih2, i, j


    ! examine the size of the array, and initialize c_ij
    !
    if (size(c_ij,1) /= 4      .or. &
        size(c_ij,2) /= 4      .or. &
        size(c_ij,3) /= ngrid0 .or. &
        size(c_ij,4) /= ngrid0) then
      write(MsgOut,*) size(c_ij,1),size(c_ij,2),size(c_ij,3),size(c_ij,4)
      call error_msg('Cmap_Cal_Coefs> ERROR: size of c_ij looks incorrect.')
    end if

    c_ij(1:4,1:4,1:ngrid0,1:ngrid0) = 0.0_wp

    ! set Ecmap and derivative values to local arrays
    !
    do idih2 = 1, ngrid0
      do idih1 = 1, ngrid0

        call bicubic_coefs(y_data(idih1:idih1+1,idih2:idih2+1),       &
                           deriv_1(idih1:idih1+1,idih2:idih2+1),      &
                           deriv_2(idih1:idih1+1,idih2:idih2+1),      &
                           deriv_12(idih1:idih1+1,idih2:idih2+1),     &
                           dlt_x, dlt_x, c) 

        do i = 1, 4
          do j = 1, 4
            c_ij(i,j,idih1,idih2) = c(i,j)
          end do
        end do

      end do
    end do

    return

  end subroutine cmap_cal_coefs

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bicubic_coefs
  !> @brief        calculate bicubic coefficients c_ij for one cell
  !! @authors      TY, YM
  !! @param[in]    y_data   : cmap energy on four edge points of the cell
  !! @param[in]    deriv_1  : dEcmap/d(dih1) on edge points
  !! @param[in]    deriv_2  : dEcmap/d(dih2) on edge points
  !! @param[in]    deriv_12 : d2Ecmap/d(dih1)d(dih2) on edge points
  !! @param[in]    dlt_x1   : size (in degree) of the cell along dih1 axis
  !! @param[in]    dlt_x2   : size (in degree) of the cell along dih2 axis
  !! @param[out]   c        : cmap coefficients 
  !!
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine bicubic_coefs(y_data, deriv_1, deriv_2, deriv_12, dlt_x1, dlt_x2, c)

    ! formal arguments
    real(wp),                intent(in)    :: y_data  (1:4)
    real(wp),                intent(in)    :: deriv_1 (1:4)
    real(wp),                intent(in)    :: deriv_2 (1:4)
    real(wp),                intent(in)    :: deriv_12(1:4)
    real(wp),                intent(in)    :: dlt_x1
    real(wp),                intent(in)    :: dlt_x2
    real(wp),                intent(out)   :: c(1:4,1:4)

    ! local variables
    real(wp)                 :: dlt_x1x2, cl(1:16), x(1:16)
    real(wp),         save   :: ainv(1:16,1:16)
    integer                  :: ainv_integer(1:16,1:16)
    integer                  :: i,j,k


    ainv_integer( 1,1:16) = (/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    ainv_integer( 2,1:16) = (/0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0/)
    ainv_integer( 3,1:16) = (/-3,0,3,0,0,0,0,0,-2,0,-1,0,0,0,0,0/)
    ainv_integer( 4,1:16) = (/2,0,-2,0,0,0,0,0,1,0,1,0,0,0,0,0/)
    ainv_integer( 5,1:16) = (/0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0/)
    ainv_integer( 6,1:16) = (/0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0/)
    ainv_integer( 7,1:16) = (/0,0,0,0,-3,0,3,0,0,0,0,0,-2,0,-1,0/)
    ainv_integer( 8,1:16) = (/0,0,0,0,2,0,-2,0,0,0,0,0,1,0,1,0/)
    ainv_integer( 9,1:16) = (/-3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0/)
    ainv_integer(10,1:16) = (/0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0/)
    ainv_integer(11,1:16) = (/9,-9,-9,9,6,3,-6,-3,6,-6,3,-3,4,2,2,1/)
    ainv_integer(12,1:16) = (/-6,6,6,-6,-4,-2,4,2,-3,3,-3,3,-2,-1,-2,-1/)
    ainv_integer(13,1:16) = (/2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0/)
    ainv_integer(14,1:16) = (/0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0/)
    ainv_integer(15,1:16) = (/-6,6,6,-6,-3,-3,3,3,-4,4,-2,2,-2,-2,-1,-1/)
    ainv_integer(16,1:16) = (/4,-4,-4,4,2,2,-2,-2,2,-2,2,-2,1,1,1,1/)

    ainv(1:16,1:16) = real(ainv_integer(1:16,1:16),wp)

    dlt_x1x2 = dlt_x1*dlt_x2

    x(1:4)   = y_data(1:4)
    x(5:8)   = deriv_1(1:4)*dlt_x1
    x(9:12)  = deriv_2(1:4)*dlt_x2
    x(13:16) = deriv_12(1:4)*dlt_x1x2

    cl = matmul(ainv,x)

    k = 0
    do i = 1, 4
      do j = 1, 4
        k = k + 1
        c(i,j) = cl(k)
      end do
    end do

    return

  end subroutine bicubic_coefs

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cal_ecmap_only
  !> @brief        calculate Ecmap using given dih1 (=psi) and dih2(=phi) values
  !! @authors      TY
  !! @param[in]    ngrid0 : the number of grid points in one direction
  !! @param[in]    dih1   : 1st dihedral angle value ( = psi) in degree
  !! @param[in]    dih2   : 2nd dihedral angle value ( = phi) in degree
  !! @param[in]    c_ij   : cmap coefficients to be used.
  !! @param[out]   ecmap  : resultant Ecmap value
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cal_ecmap_only(ngrid0, dih1, dih2, c_ij, ecmap)

    ! formal arguments
    integer,                 intent(in)    :: ngrid0
    real(wp),                intent(in)    :: dih1
    real(wp),                intent(in)    :: dih2
    real(wp),                intent(in)    :: c_ij(:,:,:,:)  ! cmap coeffs
    real(wp),                intent(out)   :: ecmap

    ! local variables
    real(wp)                 :: rdum10, rdum20, t, u, dlt
    integer                  :: i, j
    integer                  :: igrid(2)


    dlt      = 360.0_wp/real(ngrid0,wp)
    igrid(1) = int((dih1 + 180.0_wp)/dlt)
    igrid(2) = int((dih2 + 180.0_wp)/dlt)

    if (dih1 < 180.0_wp) then
      igrid(1) = igrid(1) + 1
    end if

    if (dih2 < 180.0_wp) then
      igrid(2) = igrid(2) + 1
    end if

    t = (dih1 - (dlt*real(igrid(1)-1,wp) - 180.0_wp)) / dlt
    u = (dih2 - (dlt*real(igrid(2)-1,wp) - 180.0_wp)) / dlt

    ecmap = 0.0_wp

    do i = 1, 4
      do j = 1, 4
        rdum10 = t**(real(i-1,wp))
        rdum20 = u**(real(j-1,wp))
        ecmap = ecmap + rdum10*rdum20*c_ij(i,j,igrid(1),igrid(2))
      end do
    end do

    return

  end subroutine cal_ecmap_only

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cal_decmap_ddih
  !> @brief        calc. dEcmap/d(dih) using given dih1 (=psi) and dih2 (=phi)
  !!               values
  !! @authors      TY
  !! @param[in]    ngrid0 : the number of grid points in one direction
  !! @param[in]    dih1   : 1st dihedral angle value ( = psi) in degree
  !! @param[in]    dih2   : 2nd dihedral angle value ( = phi) in degree
  !! @param[in]    c_ij   : cmap coefficients to be used.
  !! @param[out]   ddih1  : resultant dEcmap/d(dih1)
  !! @param[out]   ddih2  : resultant dEcmap/d(dih2)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cal_decmap_ddih(ngrid0, dih1, dih2, c_ij, ddih1, ddih2)

    ! formal arguments
    integer,                 intent(in)    :: ngrid0
    real(wp),                intent(in)    :: dih1
    real(wp),                intent(in)    :: dih2
    real(wp),                intent(in)    :: c_ij(:,:,:,:)  ! cmap coeffs
    real(wp),                intent(out)   :: ddih1
    real(wp),                intent(out)   :: ddih2

    ! local variables
    real(wp)                 :: rdum10, rdum20, t, u, dlt
    integer                  :: i, j
    integer                  :: igrid(2)


    dlt      = 360.0_wp/real(ngrid0,wp)
    igrid(1) = int((dih1 + 180.0_wp)/dlt)
    igrid(2) = int((dih2 + 180.0_wp)/dlt)

    if (dih1 < 180.0_wp) then
      igrid(1) = igrid(1) + 1
    end if

    if (dih2 < 180.0_wp) then
      igrid(2) = igrid(2) + 1
    end if

    t = (dih1 - (dlt*real(igrid(1)-1,wp) - 180.0_wp)) / dlt
    u = (dih2 - (dlt*real(igrid(2)-1,wp) - 180.0_wp)) / dlt

    ddih1 = 0.0_wp
    ddih2 = 0.0_wp

    do i = 2, 4
      do j = 1, 4
        rdum10 = t**(real(i-2,wp))
        rdum10 = rdum10*real(i-1,wp)
        rdum20 = u**(real(j-1,wp))
        ddih1 = ddih1 + rdum10*rdum20*c_ij(i,j,igrid(1),igrid(2))
      end do
    end do
    ddih1 = ddih1/dlt

    do i = 1, 4
      do j = 2, 4
        rdum10 = t**(real(i-1,wp))
        rdum20 = u**(real(j-2,wp))
        rdum20 = rdum20*real(j-1,wp)
        ddih2 = ddih2 + rdum10*rdum20*c_ij(i,j,igrid(1),igrid(2))
      end do
    end do

    ddih2 = ddih2/dlt

    return

  end subroutine cal_decmap_ddih

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cal_d2ecmap
  !> @brief        calc. d2Ecmap/d(dih)2 using given dih1 (=psi) and dih2 (=phi)
  !                values
  !! @authors      TY
  !! @param[in]    ngrid0 : the number of grid points in one direction
  !! @param[in]    dih1   : 1st dihedral angle value ( = psi) in degree
  !! @param[in]    dih2   : 2nd dihedral angle value ( = phi) in degree
  !! @param[in]    c_ij   : cmap coefficients to be used.
  !! @param[out]   d2dih1 : resultant d2Ecmap/d(dih1)2
  !! @param[out]   d2dih2 : resultant d2Ecmap/d(dih2)2
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cal_d2ecmap(ngrid0, dih1, dih2, c_ij, d2dih1, d2dih2)

    ! formal arguments
    integer,                 intent(in)    :: ngrid0
    real(wp),                intent(in)    :: dih1
    real(wp),                intent(in)    :: dih2
    real(wp),                intent(in)    :: c_ij(:,:,:,:) ! cmap coeffs
    real(wp),                intent(out)   :: d2dih1
    real(wp),                intent(out)   :: d2dih2

    ! local variables
    real(wp)                 :: rdum10, rdum20, t, u, dlt
    integer                  :: i, j
    integer                  :: igrid(2)


    dlt      = 360.0_wp/real(ngrid0,wp)
    igrid(1) = int((dih1 + 180.0_wp)/dlt)
    igrid(2) = int((dih2 + 180.0_wp)/dlt)

    if (dih1 < 180.0_wp) then
      igrid(1) = igrid(1) + 1
    end if

    if (dih2 < 180.0_wp) then
      igrid(2) = igrid(2) + 1
    end if

    t = (dih1 - (dlt*real(igrid(1)-1,wp) - 180.0_wp)) / dlt
    u = (dih2 - (dlt*real(igrid(2)-1,wp) - 180.0_wp)) / dlt

    d2dih1 = 0.0_wp
    d2dih2 = 0.0_wp

    do i = 3, 4
      do j = 1, 4
        rdum10 = t**(real(i-3,wp))
        rdum10 = rdum10*real(i-1,wp)*real(i-2,wp)
        rdum20 = u**(real(j-1,wp))
        d2dih1 = d2dih1 + rdum10*rdum20*c_ij(i,j,igrid(1),igrid(2))
      end do
    end do

    d2dih1 = d2dih1/(dlt*dlt)

    do i = 1, 4
      do j = 3, 4
        rdum10 = t**(real(i-1,wp))
        rdum20 = u**(real(j-3,wp))
        rdum20 = rdum20*real(j-1,wp)*real(j-2,wp)
        d2dih2 = d2dih2 + rdum10*rdum20*c_ij(i,j,igrid(1),igrid(2))
      end do
    end do

    d2dih2 = d2dih2/(dlt*dlt)

    return

  end subroutine cal_d2ecmap

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    solve_nspline
  !> @brief        calc d2y/dx2 of natural cubic spline
  !> @brief        by using the tridiagonal matrix algorithm
  !! @authors      TY, YM
  !! @param[in]    dlt_x   : distance between neighboring grid points
  !! @param[in]    y_vec   : input y value on grid points (array)
  !! @param[in]    ngrid   : the number of grid points in one direction
  !! @param[out]   rvector : resultant d2y/dx2 values on grid points
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine solve_nspline(dlt_x, y_vec, ngrid, rvector)

    ! formal arguments
    real(wp),                intent(in)    :: dlt_x
    real(wp),                intent(in)    :: y_vec(:)
    integer,                 intent(in)    :: ngrid
    real(wp),                intent(inout) :: rvector(:)

    ! local variables
    real(wp)                 :: p
    integer                  :: i, alloc_stat

    real(wp),                allocatable :: dp(:)


    ! allocate local array
    !
    allocate (dp(ngrid), stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc


    ! examine the size of arrays
    !
    if (size(y_vec) /= ngrid) then
      write(MsgOut,*)  'Solve_Nspline> ERROR(1).  ', size(y_vec), ngrid
      call error_msg('Solve_Nspline> unexpected array size. ')
    end if

    if (size(rvector) /= ngrid) then
      write(MsgOut,*)  'Solve_Nspline> ERROR(2).  ', size(rvector), ngrid
      call error_msg('Solve_Nspline> unexpected array size. ')
    end if

    ! boundary condition at the lower limit
    !
    rvector(1) = 0.0_wp
    dp(1) = 0.0_wp

    ! first sweep
    !
    do i = 2, ngrid - 1
      p = 0.5_wp*rvector(i-1) + 2.0_wp
      rvector(i) = -0.5_wp/p
      dp(i) = (y_vec(i+1) - 2.0_wp*y_vec(i) + y_vec(i-1))/dlt_x
      dp(i) = (3.0_wp*dp(i)/dlt_x - 0.5_wp*dp(i-1))/p
    end do

    ! boundary condition at the higher limit
    !
    rvector(ngrid) = 0.0_wp

    ! backward substitution
    ! 
    do i = ngrid-1, 1, -1
      rvector(i) = rvector(i)*rvector(i+1) + dp(i)
    end do

    return

  end subroutine solve_nspline

end module dihedral_libs_mod
