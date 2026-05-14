!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_assign_velocity_mod
!> @brief   generate initial velocities and remove trans- and rotational motions
!! @authors Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_assign_velocity_mod

  use random_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use sp_dynvars_mod
#ifdef MPI
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: initial_velocity
  public  :: stop_trans_rotation
  private :: reduce_com

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    initial_velocity
  !> @brief        generate initial velocities
  !! @authors      JJ
  !! @param[in]    firstt   : initial temperature
  !! @param[in]    natom    : number of atom
  !! @param[in]    ncell    : total number of cells in each domain
  !! @param[in]    id_g2l   : map of atom index (global -> local)
  !! @param[in]    mass     : mass in each domain
  !! @param[inout] iseed    : random number seed
  !! @param[inout] velocity : velocities generated
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine initial_velocity(firstt, natom, ncell, id_g2l, mass, &
                              iseed, velocity)

    ! formal arguments
    real(dp),                intent(in)    :: firstt
    integer,                 intent(in)    :: natom
    integer,                 intent(in)    :: ncell
    integer,                 intent(in)    :: id_g2l(:,:)
    real(dp),                intent(in)    :: mass(:,:)
    integer,                 intent(inout) :: iseed
    real(dp),                intent(inout) :: velocity(:,:,:)

    ! local variables
    real(dp)                 :: alpha1, beta1, alpha2, beta2, alpha3, beta3
    real(dp)                 :: sigma, myu, kBT
    integer                  :: i, ix, icel


    if (main_rank .or. replica_main_rank) then
      write(DynvarsOut,'(A)') 'Initial_Velocity> Generate initial velocities'
      write(DynvarsOut,'(A20,I10)')   '  iseed           = ', iseed
      write(DynvarsOut,'(A20,F10.3)') '  temperature     = ', firstt
      write(DynvarsOut,'(A)') ' '
    end if

    ! generate initial velocities 
    !
    kBT = KBOLTZ * firstt
    myu = 0.0_dp

    do i = 1, natom

       alpha1 = sqrt( -2.0_dp * log(random_get_legacy(iseed)) )
       beta1  = cos ( 2.0_dp * PI * random_get_legacy(iseed) )

       alpha2 = sqrt( -2.0_dp * log(random_get_legacy(iseed)) )
       beta2  = cos ( 2.0_dp * PI * random_get_legacy(iseed) )

       alpha3 = sqrt( -2.0_dp * log(random_get_legacy(iseed)) )
       beta3  = cos ( 2.0_dp * PI * random_get_legacy(iseed) )

       if (id_g2l(1,i) <= ncell .and. id_g2l(2,i) /=0) then
          icel  = id_g2l(1,i)
          ix    = id_g2l(2,i)
          if (abs(mass(ix,icel)) > EPS) then
            sigma = sqrt( kBT / mass(ix,icel) )
          else
            sigma = 0.0_dp
          end if
          velocity(1,ix,icel) = sigma * alpha1 * beta1 + myu
          velocity(2,ix,icel) = sigma * alpha2 * beta2 + myu
          velocity(3,ix,icel) = sigma * alpha3 * beta3 + myu
       end if

    end do
   
    return 

  end subroutine initial_velocity

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    stop_trans_rotation
  !> @brief        remove trans and rotational motion about the center of mass
  !! @authors      JJ
  !! @param[in]    ncell      : number of cell
  !! @param[in]    natom      : number of atom in each cell
  !! @param[in]    stop_trans : flag for stop translational motion
  !! @param[in]    stop_rot   : flag for stop rotational motion
  !! @param[in]    mass       : mass in each domain
  !! @param[in]    coord      : coordinates in each domain
  !! @param[inout] velocity   : velocities in each domain
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine stop_trans_rotation(ncell, natom, stop_trans, stop_rot, &
                                 mass, coord, velocity)

    ! formal arguments
    integer,                 intent(in)    :: ncell
    integer,                 intent(in)    :: natom(:)
    logical,                 intent(in)    :: stop_trans
    logical,                 intent(in)    :: stop_rot
    real(dp),                intent(in)    :: mass(:,:)
    real(dp),                intent(in)    :: coord(:,:,:)
    real(dp),                intent(inout) :: velocity(:,:,:)

    ! local variables
    real(dp)                 :: inertia(3,3), inv_inertia(3,3)
    real(dp)                 :: total_mass, ccom(3), vcom(3)
    real(dp)                 :: l_angular(1:3)
    real(dp)                 :: c(3), v(3), omega(3)
    real(dp)                 :: xx, xy, xz, yy, yz, zz, x, y, z
    integer                  :: i, j, k, ix

    ! calculate the center of mass of coordinates and velocities
    !
    ccom(1:3)  = 0.0_dp
    vcom(1:3)  = 0.0_dp
    total_mass = 0.0_dp
    do i = 1, ncell
      do ix = 1, natom(i)
        ccom(1:3)  = ccom(1:3) + coord(1:3,ix,i)   *mass(ix,i)
        vcom(1:3)  = vcom(1:3) + velocity(1:3,ix,i)*mass(ix,i)
        total_mass = total_mass + mass(ix,i)
      end do
    end do

    call reduce_com(ccom, vcom, total_mass)
    
    ccom(1:3) = ccom(1:3) / total_mass
    vcom(1:3) = vcom(1:3) / total_mass

    ! calculate the angular momentum
    !
    l_angular(1:3) = 0.0_dp
    do i = 1, ncell
      do ix = 1, natom(i)
        c(1:3) = coord   (1:3,ix,i) - ccom(1:3)
        v(1:3) = velocity(1:3,ix,i) - vcom(1:3)
        l_angular(1) = l_angular(1) + (c(2)*v(3) - c(3)*v(2))*mass(ix,i)
        l_angular(2) = l_angular(2) + (c(3)*v(1) - c(1)*v(3))*mass(ix,i)
        l_angular(3) = l_angular(3) + (c(1)*v(2) - c(2)*v(1))*mass(ix,i)
      end do
    end do
 
#ifdef MPI
    call mpi_allreduce(mpi_in_place, l_angular, 3, mpi_real8, &
                       mpi_sum, mpi_comm_country, ierror)
#endif

    ! calculate the inertia tensor (I)
    !
    xx = 0.0_dp
    xy = 0.0_dp
    xz = 0.0_dp
    yy = 0.0_dp
    yz = 0.0_dp
    zz = 0.0_dp

    do i = 1, ncell
      do ix = 1, natom(i)
        c(1:3) = coord(1:3,ix,i) - ccom(1:3)
        xx = xx + c(1)*c(1)*mass(ix,i)
        xy = xy + c(1)*c(2)*mass(ix,i)
        xz = xz + c(1)*c(3)*mass(ix,i)
        yy = yy + c(2)*c(2)*mass(ix,i)
        yz = yz + c(2)*c(3)*mass(ix,i)
        zz = zz + c(3)*c(3)*mass(ix,i)
      end do
    end do

    inertia(1,1) =   yy + zz
    inertia(1,2) = - xy
    inertia(1,3) = - xz
    inertia(2,1) = - xy
    inertia(2,2) =   xx + zz
    inertia(2,3) = - yz
    inertia(3,1) = - xz
    inertia(3,2) = - yz
    inertia(3,3) =   xx + yy

#ifdef MPI
    call mpi_allreduce(mpi_in_place, inertia, 9, mpi_real8, &
                       mpi_sum, mpi_comm_country, ierror)
#endif

    ! calculate the inverse matrix of the inertia tensor
    !
    call compute_inverse_matrix(3, inertia, inv_inertia)
 
    ! calculate the angular velocity (OMG)
    !
    omega(1) = inv_inertia(1,1)*l_angular(1) + inv_inertia(1,2)*l_angular(2) &
             + inv_inertia(1,3)*l_angular(3)
    omega(2) = inv_inertia(2,1)*l_angular(1) + inv_inertia(2,2)*l_angular(2) &
             + inv_inertia(2,3)*l_angular(3)
    omega(3) = inv_inertia(3,1)*l_angular(1) + inv_inertia(3,2)*l_angular(2) &
             + inv_inertia(3,3)*l_angular(3)

    ! remove translational motion
    !
    if (stop_trans) then
      do i = 1, ncell
        do ix = 1, natom(i)
          velocity(1:3,ix,i) = velocity(1:3,ix,i) - vcom(1:3)
        end do
      end do
    end if

    ! remove rotational motion
    !
    if (stop_rot) then
      do i = 1, ncell
        do ix = 1, natom(i)
          c(1:3) = coord(1:3,ix,i) - ccom(1:3)
          velocity(1,ix,i) = velocity(1,ix,i) - (omega(2)*c(3) - omega(3)*c(2))
          velocity(2,ix,i) = velocity(2,ix,i) - (omega(3)*c(1) - omega(1)*c(3))
          velocity(3,ix,i) = velocity(3,ix,i) - (omega(1)*c(2) - omega(2)*c(1))
        end do
      end do
    end if

    return

  end subroutine stop_trans_rotation

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    reduce_com
  !> @brief
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine reduce_com(val1, val2, val3)

    ! formal arguments
    real(dp),                intent(inout) :: val1(:), val2(:), val3

    ! local variables
    real(dp)                      :: before_reduce(7), after_reduce(7)


    before_reduce(1:3) = val1(1:3)
    before_reduce(4:6) = val2(1:3)
    before_reduce(7)   = val3

#ifdef MPI
    call mpi_allreduce(before_reduce, after_reduce, 7, mpi_real8, &
                       mpi_sum, mpi_comm_city, ierror)
#else
    after_reduce(1:7) = before_reduce(1:7)
#endif

    val1(1:3)    = after_reduce(1:3)
    val2(1:3)    = after_reduce(4:6)
    val3         = after_reduce(7)

    return

  end subroutine reduce_com


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
    integer,      intent(in)  :: n
    real(dp),     intent(in)  :: m(n,n)
    real(dp),     intent(out) :: inv_m(n,n)

    ! local variables
    integer                   :: i, j, k, imax
    real(dp)                  :: a(n,2*n), amax, swap, coef, pivot

    do i = 1, n
      do j = 1, n
        a(i,j) = m(i,j)
        if (i == j) a(i,j+n) = 1.0_dp
        if (i /= j) a(i,j+n) = 0.0_dp
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

      pivot = 1.0_dp/a(k,k)
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

end module sp_assign_velocity_mod
