!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ma_matrix_mod
!> @brief   matrix routines
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ma_matrix_mod

  use messages_mod
  use constants_mod

  implicit none

  public  :: inv_matrix
  public  :: mul_matrix

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    inv_matrix
  !> @brief        get inverse matrix
  !! @authors      NT
  !! @param[in]    A  : input matrix
  !! @param[in]    AI : inverse matrix
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine inv_matrix(A, AI)

    ! parameters
    integer,                 parameter     :: NB = 64

    ! formal arguments
    real(wp),                intent(in)    :: A (:,:)
    real(wp),                intent(inout) :: AI(:,:)

    ! local variables
    integer                  :: info, m, n
    integer                  :: LWORK
    real(8),     allocatable :: WORK(:) ! dimension LWORK
    integer,     allocatable :: IPIV(:) ! dimension N


    m = size(A(:,1))
    n = size(A(1,:))

    lwork = n*nb
    allocate(ipiv(n))
    allocate(work(lwork))

    ai = a
    call dgetrf(m, n, ai, m, ipiv(1:min(m,n)), info)
    call dgetri(n, ai, m, ipiv, work, lwork, info)

    deallocate(ipiv)
    deallocate(work)

    if (info /= 0) then
      write(*,*) 'Inv_Matrix> ERROR.'
    end if

    return

  end subroutine inv_matrix

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    mul_matrix
  !> @brief        compute V0 = A x VI
  !! @authors      NT
  !! @param[in]    A  : input matrix
  !! @param[in]    VI : input vector
  !! @param[out]   VO : output vector
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine mul_matrix(A, VI, VO)

    ! formal arguments
    real(wp),                intent(in)    :: A (:,:)
    real(wp),                intent(in)    :: VI(:)
    real(wp),                intent(inout) :: VO(:)

    ! local variables
    integer                  :: i, j, m, n


    m = size(A(:,1))
    n = size(A(1,:))

    VO(:) = 0.0_wp

    do j = 1, n
      do i = 1, m
        VO(j) = VO(j) + A(i, j) * VI(i)
      end do
    end do

    return

  end subroutine mul_matrix

end module ma_matrix_mod
