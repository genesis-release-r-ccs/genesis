!--------1---------2---------3---------4---------5---------6---------7---------8
! 
!  Module   sp_migration_mod
!> @brief   migration of atoms and energy functions in each domain
!! @authors Jaewoon Jung (JJ) , Chigusa Kobayshi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
! 
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_migration_mod

  use sp_constraints_str_mod
  use sp_boundary_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public :: update_outgoing_atom
  public :: update_outgoing_water
  public :: update_outgoing_HGr
  public :: update_incoming_atom
  public :: update_incoming_water
  public :: update_incoming_HGr
  public :: update_enefunc_enm
  public :: update_outgoing_enefunc_bond
  public :: update_incoming_enefunc_bond
  public :: update_outgoing_enefunc_angl
  public :: update_incoming_enefunc_angl
  public :: update_outgoing_enefunc_dihe
  public :: update_incoming_enefunc_dihe
  public :: update_outgoing_enefunc_rb_dihe
  public :: update_incoming_enefunc_rb_dihe
  public :: update_outgoing_enefunc_impr
  public :: update_incoming_enefunc_impr
  public :: update_outgoing_enefunc_cmap
  public :: update_incoming_enefunc_cmap
  public :: update_outgoing_enefunc_restraint
  public :: update_incoming_enefunc_restraint
  public :: update_outgoing_enefunc_fitting
  public :: update_incoming_enefunc_fitting
  ! FEP
  public :: update_outgoing_atom_fep
  public :: update_outgoing_water_fep
  public :: update_outgoing_HGr_fep
  public :: update_incoming_atom_fep
  public :: update_incoming_water_fep
  public :: update_incoming_HGr_fep
  public :: update_outgoing_enefunc_bond_fep
  public :: update_incoming_enefunc_bond_fep
  public :: update_outgoing_enefunc_angl_fep
  public :: update_incoming_enefunc_angl_fep
  public :: update_outgoing_enefunc_dihe_fep
  public :: update_incoming_enefunc_dihe_fep
  public :: update_outgoing_enefunc_rb_dihe_fep
  public :: update_incoming_enefunc_rb_dihe_fep
  public :: update_outgoing_enefunc_impr_fep
  public :: update_incoming_enefunc_impr_fep
  public :: update_outgoing_enefunc_cmap_fep
  public :: update_incoming_enefunc_cmap_fep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_atom
  !> @brief        check particles (not bonded to hydrogen) going other cells
  !! @authors      JJ
  !! @param[in]    boundary    : boundary condition information
  !! @param[in]    constraints : constraints information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_atom(boundary, constraints, domain)

    ! formal arguments
    type(s_boundary),    target, intent(in)    :: boundary
    type(s_constraints), target, intent(in)    :: constraints
    type(s_domain),      target, intent(inout) :: domain

    ! local variable
    real(wip)                    :: x_shift, y_shift, z_shift
    real(wip)                    :: move(3)
    integer                      :: i, k, ix, icx, icy, icz, icel, ncel
    integer                      :: icel_local, icel_bd

    real(wip),           pointer :: bsize_x, bsize_y, bsize_z
    real(wip),           pointer :: csize_x, csize_y, csize_z
    real(wip),           pointer :: coord(:,:,:), velocity(:,:,:)
    real(wp),            pointer :: charge(:,:)
    real(wip),           pointer :: mass(:,:)
    real(wip),           pointer :: buf_real(:,:,:)
    integer,             pointer :: ncel_x, ncel_y, ncel_z, ncel_local, ncel_bd
    integer,             pointer :: nsolute(:)
    integer(int2),       pointer :: cell_g2l(:), cell_g2b(:)
    integer,             pointer :: atmcls(:,:), id_l2g(:,:), id_l2g_sol(:,:)
    integer,             pointer :: ptl_add(:), ptl_exit(:), ptlindex(:,:)
    integer,             pointer :: buf_int(:,:,:)


    bsize_x    => boundary%box_size_x
    bsize_y    => boundary%box_size_y
    bsize_z    => boundary%box_size_z
    ncel_x     => boundary%num_cells_x
    ncel_y     => boundary%num_cells_y
    ncel_z     => boundary%num_cells_z
    csize_x    => boundary%cell_size_x
    csize_y    => boundary%cell_size_y
    csize_z    => boundary%cell_size_z

    nsolute    => constraints%No_HGr

    ncel_local => domain%num_cell_local
    ncel_bd    => domain%num_cell_boundary
    cell_g2l   => domain%cell_g2l
    cell_g2b   => domain%cell_g2b
    coord      => domain%coord
    velocity   => domain%velocity
    charge     => domain%charge
    mass       => domain%mass
    atmcls     => domain%atom_cls_no
    id_l2g     => domain%id_l2g
    id_l2g_sol => domain%id_l2g_solute
    ptl_add    => domain%ptl_add
    ptl_exit   => domain%ptl_exit
    ptlindex   => domain%ptl_exit_index
    buf_int    => domain%buf_integer
    buf_real   => domain%buf_real


    ! initializaiton
    !
    ncel = ncel_local + ncel_bd
    ptl_exit(1:ncel_local) = 0
    ptl_add(1:ncel) = 0

    ! Check outgoing particles
    !
    do i = 1, ncel_local

      k = 0
      do ix = 1, nsolute(i)

        x_shift = coord(1,ix,i) - boundary%origin_x
        y_shift = coord(2,ix,i) - boundary%origin_y
        z_shift = coord(3,ix,i) - boundary%origin_z

        !coordinate shifted to the first quadrant and set into the boundary box
        move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
        move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
        move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)

        x_shift = x_shift + move(1)
        y_shift = y_shift + move(2)
        z_shift = z_shift + move(3)

        !assign which cell
        icx = int(x_shift/csize_x)
        icy = int(y_shift/csize_y)
        icz = int(z_shift/csize_z)
        if (icx == ncel_x) icx = icx - 1
        if (icy == ncel_y) icy = icy - 1
        if (icz == ncel_z) icz = icz - 1
        icel = 1 + icx + icy*ncel_x + icz*ncel_x*ncel_y
        icel_local = cell_g2l(icel)
        icel_bd    = cell_g2b(icel)

        if (icel_local /= i) then

          ptl_exit(i) = ptl_exit(i) + 1
          ptlindex(ptl_exit(i),i) = ix

          if (icel_local /= 0) then

            ptl_add (icel_local) = ptl_add(icel_local) + 1
            buf_real(1,ptl_add(icel_local),icel_local) = coord(1,ix,i)
            buf_real(2,ptl_add(icel_local),icel_local) = coord(2,ix,i)
            buf_real(3,ptl_add(icel_local),icel_local) = coord(3,ix,i)
            buf_real(4,ptl_add(icel_local),icel_local) = velocity(1,ix,i)
            buf_real(5,ptl_add(icel_local),icel_local) = velocity(2,ix,i)
            buf_real(6,ptl_add(icel_local),icel_local) = velocity(3,ix,i)
            buf_real(7,ptl_add(icel_local),icel_local) = charge(ix,i)
            buf_real(8,ptl_add(icel_local),icel_local) = mass(ix,i)
            buf_int (1,ptl_add(icel_local),icel_local) = atmcls(ix,i)
            buf_int (2,ptl_add(icel_local),icel_local) = id_l2g(ix,i)
            buf_int (3,ptl_add(icel_local),icel_local) = id_l2g_sol(ix,i)

          else if (icel_bd /= 0) then

            icel_bd = icel_bd + ncel_local
            ptl_add (icel_bd) = ptl_add(icel_bd) + 1
            buf_real(1,ptl_add(icel_bd),icel_bd) = coord(1,ix,i)
            buf_real(2,ptl_add(icel_bd),icel_bd) = coord(2,ix,i)
            buf_real(3,ptl_add(icel_bd),icel_bd) = coord(3,ix,i)
            buf_real(4,ptl_add(icel_bd),icel_bd) = velocity(1,ix,i)
            buf_real(5,ptl_add(icel_bd),icel_bd) = velocity(2,ix,i)
            buf_real(6,ptl_add(icel_bd),icel_bd) = velocity(3,ix,i)
            buf_real(7,ptl_add(icel_bd),icel_bd) = charge(ix,i)
            buf_real(8,ptl_add(icel_bd),icel_bd) = mass(ix,i)
            buf_int (1,ptl_add(icel_bd),icel_bd) = atmcls(ix,i)
            buf_int (2,ptl_add(icel_bd),icel_bd) = id_l2g(ix,i)
            buf_int (3,ptl_add(icel_bd),icel_bd) = id_l2g_sol(ix,i)

          end if

        end if

       end do
    end do

    return

  end subroutine update_outgoing_atom

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_water
  !> @brief        check water particles going/remaining cells
  !! @authors      JJ
  !! @param[in]    water_atom : number of water atoms
  !! @param[in]    boundary   : boundary condition information
  !! @param[inout] domain     : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_water(water_atom, boundary, domain)

    ! formal arguments
    integer,                  intent(in)    :: water_atom
    type(s_boundary), target, intent(in)    :: boundary
    type(s_domain),   target, intent(inout) :: domain

    ! local variable
    real(wip)                 :: x_shift, y_shift, z_shift
    real(wip)                 :: move(3)
    integer                   :: i, k, iwater, list, ix
    integer                   :: icx, icy, icz, icel, ncel
    integer                   :: icel_local, icel_bd

    real(wip),        pointer :: bsize_x, bsize_y, bsize_z
    real(wip),        pointer :: csize_x, csize_y, csize_z
    real(wip),        pointer :: coord(:,:,:), velocity(:,:,:)
    real(wip),        pointer :: water_move_real(:,:,:), water_stay_real(:,:,:)
    integer,          pointer :: ncel_x, ncel_y, ncel_z, ncel_local, ncel_bd
    integer,          pointer :: nwater(:), water_list(:,:,:)
    integer,          pointer :: id_l2g(:,:)
    integer(int2),    pointer :: cell_g2l(:), cell_g2b(:)
    integer,          pointer :: water_move(:), water_stay(:)
    integer,          pointer :: water_move_int(:,:,:), water_stay_int(:,:,:)


    bsize_x         => boundary%box_size_x
    bsize_y         => boundary%box_size_y
    bsize_z         => boundary%box_size_z
    ncel_x          => boundary%num_cells_x
    ncel_y          => boundary%num_cells_y
    ncel_z          => boundary%num_cells_z
    csize_x         => boundary%cell_size_x
    csize_y         => boundary%cell_size_y
    csize_z         => boundary%cell_size_z

    ncel_local      => domain%num_cell_local
    ncel_bd         => domain%num_cell_boundary
    nwater          => domain%num_water
    water_list      => domain%water_list
    cell_g2l        => domain%cell_g2l
    cell_g2b        => domain%cell_g2b
    coord           => domain%coord
    velocity        => domain%velocity
    id_l2g          => domain%id_l2g
    water_move      => domain%water%move
    water_stay      => domain%water%stay
    water_move_real => domain%water%move_real
    water_stay_real => domain%water%stay_real
    water_move_int  => domain%water%move_integer
    water_stay_int  => domain%water%stay_integer


    ! initializaiton
    !
    ncel = ncel_local + ncel_bd
    water_stay(1:ncel_local) = 0
    water_move(1:ncel) = 0

    ! Check outgoing particles
    !
    do i = 1, ncel_local

      k = 0
      do iwater = 1, nwater(i)

        ix = water_list(1,iwater,i)
        x_shift = coord(1,ix,i) - boundary%origin_x
        y_shift = coord(2,ix,i) - boundary%origin_y
        z_shift = coord(3,ix,i) - boundary%origin_z

        !coordinate shifted to the first quadrant and set into the boundary box
        move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
        move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
        move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)

        x_shift = x_shift + move(1)
        y_shift = y_shift + move(2)
        z_shift = z_shift + move(3)

        !assign which cell
        icx = int(x_shift/csize_x)
        icy = int(y_shift/csize_y)
        icz = int(z_shift/csize_z)
        if (icx == ncel_x) icx = icx - 1
        if (icy == ncel_y) icy = icy - 1
        if (icz == ncel_z) icz = icz - 1
        icel = 1 + icx + icy*ncel_x + icz*ncel_x*ncel_y
        icel_local = cell_g2l(icel)
        icel_bd    = cell_g2b(icel)

        ! atom index
        if (icel_local /= i) then

          if (icel_local /= 0) then

            water_move(icel_local) = water_move(icel_local) + 1

            do list = 1, water_atom
              ix = water_list(list,iwater,i)
              water_move_real(6*(list-1)+1,water_move(icel_local),icel_local) &
                 = coord(1,ix,i)
              water_move_real(6*(list-1)+2,water_move(icel_local),icel_local) &
                 = coord(2,ix,i)
              water_move_real(6*(list-1)+3,water_move(icel_local),icel_local) &
                 = coord(3,ix,i)
              water_move_real(6*(list-1)+4,water_move(icel_local),icel_local) &
                 = velocity(1,ix,i)
              water_move_real(6*(list-1)+5,water_move(icel_local),icel_local) &
                 = velocity(2,ix,i)
              water_move_real(6*(list-1)+6,water_move(icel_local),icel_local) &
                 = velocity(3,ix,i)
              water_move_int(list,water_move(icel_local),icel_local) &
                 = id_l2g(ix,i)
            end do

          else if (icel_bd /= 0) then

            icel_bd = icel_bd + ncel_local
            water_move(icel_bd) = water_move(icel_bd) + 1

            do list = 1, water_atom
              ix = water_list(list,iwater,i)
              water_move_real(6*(list-1)+1,water_move(icel_bd),icel_bd)       &
                 = coord(1,ix,i)
              water_move_real(6*(list-1)+2,water_move(icel_bd),icel_bd)       &
                 = coord(2,ix,i)
              water_move_real(6*(list-1)+3,water_move(icel_bd),icel_bd)       &
                 = coord(3,ix,i)
              water_move_real(6*(list-1)+4,water_move(icel_bd),icel_bd)       &
                 = velocity(1,ix,i)
              water_move_real(6*(list-1)+5,water_move(icel_bd),icel_bd)       &
                 = velocity(2,ix,i)
              water_move_real(6*(list-1)+6,water_move(icel_bd),icel_bd)       &
                 = velocity(3,ix,i)
              water_move_int(list,water_move(icel_bd),icel_bd) &
                 = id_l2g(ix,i)
            end do
          end if

        else

          water_stay(i) = water_stay(i) + 1

          do list = 1, water_atom
            ix = water_list(list,iwater,i)
            water_stay_real(6*(list-1)+1,water_stay(i),i) = coord(1,ix,i)
            water_stay_real(6*(list-1)+2,water_stay(i),i) = coord(2,ix,i)
            water_stay_real(6*(list-1)+3,water_stay(i),i) = coord(3,ix,i)
            water_stay_real(6*(list-1)+4,water_stay(i),i) = velocity(1,ix,i)
            water_stay_real(6*(list-1)+5,water_stay(i),i) = velocity(2,ix,i)
            water_stay_real(6*(list-1)+6,water_stay(i),i) = velocity(3,ix,i)
            water_stay_int(list,water_stay(i),i) = id_l2g(ix,i)
          end do

        end if

      end do

    end do

    return

  end subroutine update_outgoing_water

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_HGr
  !> @brief        check hydrogen-bonded particles going/remaining cells
  !! @authors      JJ
  !! @param[in]    boundary    : boundary condition information
  !! @param[inout] constraints : constraints information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_HGr(boundary, constraints, domain)

    ! formal arguments
    type(s_boundary),    target, intent(in)    :: boundary
    type(s_constraints), target, intent(inout) :: constraints
    type(s_domain),      target, intent(inout) :: domain

    ! local variable
    real(wip)                    :: x_shift, y_shift, z_shift
    real(wip)                    :: move(3)
    integer                      :: i, j, k, list, ix
    integer                      :: num_move, num_stay
    integer                      :: icx, icy, icz, icel, ncel
    integer                      :: icel_local, icel_bd

    real(wip),           pointer :: bsize_x, bsize_y, bsize_z
    real(wip),           pointer :: csize_x, csize_y, csize_z
    real(wip),           pointer :: coord(:,:,:), vel(:,:,:)
    real(wp),            pointer :: charge(:,:)
    real(wip),           pointer :: mass(:,:)
    real(wip),           pointer :: HGr_bond_dist(:,:,:,:)
    real(wip),           pointer :: HGr_move_real(:,:,:,:)
    real(wip),           pointer :: HGr_stay_real(:,:,:,:)
    integer,             pointer :: ncel_x, ncel_y, ncel_z
    integer,             pointer :: ncel_local, ncel_bd
    integer,             pointer :: connect
    integer,             pointer :: atmcls(:,:)
    integer,             pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)
    integer,             pointer :: id_l2g(:,:), id_l2g_sol(:,:)
    integer(int2),       pointer :: cell_g2l(:), cell_g2b(:)
    integer,             pointer :: HGr_move(:,:), HGr_stay(:,:)
    integer,             pointer :: HGr_move_int(:,:,:,:), HGr_stay_int(:,:,:,:)


    bsize_x       => boundary%box_size_x
    bsize_y       => boundary%box_size_y
    bsize_z       => boundary%box_size_z
    ncel_x        => boundary%num_cells_x
    ncel_y        => boundary%num_cells_y
    ncel_z        => boundary%num_cells_z
    csize_x       => boundary%cell_size_x
    csize_y       => boundary%cell_size_y
    csize_z       => boundary%cell_size_z

    HGr_local     => constraints%HGr_local
    HGr_bond_list => constraints%HGr_bond_list
    HGr_move      => constraints%HGr_move
    HGr_stay      => constraints%HGr_stay
    HGr_bond_dist => constraints%HGr_bond_dist
    HGr_move_real => constraints%HGr_move_real
    HGr_stay_real => constraints%HGr_stay_real
    HGr_move_int  => constraints%HGr_move_int
    HGr_stay_int  => constraints%HGr_stay_int
    connect       => constraints%connect

    ncel_local    => domain%num_cell_local
    ncel_bd       => domain%num_cell_boundary
    cell_g2l      => domain%cell_g2l
    cell_g2b      => domain%cell_g2b
    coord         => domain%coord
    vel           => domain%velocity
    charge        => domain%charge
    mass          => domain%mass
    id_l2g        => domain%id_l2g
    id_l2g_sol    => domain%id_l2g_solute
    atmcls        => domain%atom_cls_no


    ! initializaiton
    !
    ncel = ncel_local + ncel_bd
    HGr_stay(1:connect,1:ncel_local) = 0
    HGr_move(1:connect,1:ncel) = 0

    ! Check outgoing particles
    !
    do i = 1, ncel_local
      do j = 1, connect
        do k = 1, HGr_local(j,i)

          ix = HGr_bond_list(1,k,j,i)
          x_shift = coord(1,ix,i) - boundary%origin_x
          y_shift = coord(2,ix,i) - boundary%origin_y
          z_shift = coord(3,ix,i) - boundary%origin_z

          move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
          move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
          move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)

          x_shift = x_shift + move(1)
          y_shift = y_shift + move(2)
          z_shift = z_shift + move(3)

          !assign which cell
          icx = int(x_shift/csize_x)
          icy = int(y_shift/csize_y)
          icz = int(z_shift/csize_z)
          if (icx == ncel_x) icx = icx - 1
          if (icy == ncel_y) icy = icy - 1
          if (icz == ncel_z) icz = icz - 1
          icel = 1 + icx + icy*ncel_x + icz*ncel_x*ncel_y
          icel_local = cell_g2l(icel)
          icel_bd    = cell_g2b(icel)

          ! atom index
          if (icel_local /= i) then
            if (icel_local /= 0) then

              HGr_move(j,icel_local) = HGr_move(j,icel_local) + 1
              num_move = HGr_move(j,icel_local)

              do list = 1, j+1
                ix = HGr_bond_list(list,k,j,i)
                HGr_move_real(9*(list-1)+1,num_move,j,icel_local) = &
                     coord(1,ix,i)
                HGr_move_real(9*(list-1)+2,num_move,j,icel_local) = &
                     coord(2,ix,i)
                HGr_move_real(9*(list-1)+3,num_move,j,icel_local) = &
                     coord(3,ix,i)
                HGr_move_real(9*(list-1)+4,num_move,j,icel_local) = vel(1,ix,i)
                HGr_move_real(9*(list-1)+5,num_move,j,icel_local) = vel(2,ix,i)
                HGr_move_real(9*(list-1)+6,num_move,j,icel_local) = vel(3,ix,i)
                HGr_move_real(9*(list-1)+7,num_move,j,icel_local) = charge(ix,i)
                HGr_move_real(9*(list-1)+8,num_move,j,icel_local) = mass(ix,i)
                HGr_move_real(9*(list-1)+9,num_move,j,icel_local) = &
                     HGr_bond_dist(list,k,j,i)
                HGr_move_int(3*(list-1)+1,num_move,j,icel_local)  = atmcls(ix,i)
                HGr_move_int(3*(list-1)+2,num_move,j,icel_local)  = id_l2g(ix,i)
                HGr_move_int(3*(list-1)+3,num_move,j,icel_local)  = id_l2g_sol(ix,i)
              end do

            else if (icel_bd /= 0) then

              icel_bd = icel_bd + ncel_local
              HGr_move(j,icel_bd) = HGr_move(j,icel_bd) + 1
              num_move = HGr_move(j,icel_bd)

              do list = 1, j+1
                ix = HGr_bond_list(list,k,j,i)
                hgr_move_real(9*(list-1)+1,num_move,j,icel_bd) = coord(1,ix,i)
                hgr_move_real(9*(list-1)+2,num_move,j,icel_bd) = coord(2,ix,i)
                hgr_move_real(9*(list-1)+3,num_move,j,icel_bd) = coord(3,ix,i)
                hgr_move_real(9*(list-1)+4,num_move,j,icel_bd) = vel(1,ix,i)
                hgr_move_real(9*(list-1)+5,num_move,j,icel_bd) = vel(2,ix,i)
                hgr_move_real(9*(list-1)+6,num_move,j,icel_bd) = vel(3,ix,i)
                hgr_move_real(9*(list-1)+7,num_move,j,icel_bd) = charge(ix,i)
                hgr_move_real(9*(list-1)+8,num_move,j,icel_bd) = mass(ix,i)
                HGr_move_real(9*(list-1)+9,num_move,j,icel_bd) = &
                     HGr_bond_dist(list,k,j,i)
                hgr_move_int(3*(list-1)+1,num_move,j,icel_bd)  = atmcls(ix,i)
                hgr_move_int(3*(list-1)+2,num_move,j,icel_bd)  = id_l2g(ix,i)
                hgr_move_int(3*(list-1)+3,num_move,j,icel_bd)  = id_l2g_sol(ix,i)
              end do
            end if

          else

            hgr_stay(j,i) = hgr_stay(j,i) + 1
            num_stay = hgr_stay(j,i)

            do list = 1, j+1
              ix = HGr_bond_list(list,k,j,i)
              hgr_stay_real(9*(list-1)+1,num_stay,j,i) = coord(1,ix,i)
              hgr_stay_real(9*(list-1)+2,num_stay,j,i) = coord(2,ix,i)
              hgr_stay_real(9*(list-1)+3,num_stay,j,i) = coord(3,ix,i)
              hgr_stay_real(9*(list-1)+4,num_stay,j,i) = vel(1,ix,i)
              hgr_stay_real(9*(list-1)+5,num_stay,j,i) = vel(2,ix,i)
              hgr_stay_real(9*(list-1)+6,num_stay,j,i) = vel(3,ix,i)
              hgr_stay_real(9*(list-1)+7,num_stay,j,i) = charge(ix,i)
              hgr_stay_real(9*(list-1)+8,num_stay,j,i) = mass(ix,i)
              hgr_stay_real(9*(list-1)+9,num_stay,j,i) = &
                   HGr_bond_dist(list,k,j,i)
              hgr_stay_int (3*(list-1)+1,num_stay,j,i) = atmcls(ix,i)
              hgr_stay_int (3*(list-1)+2,num_stay,j,i) = id_l2g(ix,i)
              hgr_stay_int (3*(list-1)+3,num_stay,j,i) = id_l2g_sol(ix,i)
            end do

          end if

        end do
      end do
    end do

    return

  end subroutine update_outgoing_HGr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_atom
  !> @brief        check particles (not bonded to hydrogen) incoming to each
  !!               cell
  !! @authors      JJ
  !! @param[inout] constraints : constraints information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_incoming_atom(constraints, domain)

    ! formal arguments
    type(s_constraints), target, intent(inout) :: constraints
    type(s_domain),      target, intent(inout) :: domain

    ! local variables
    integer                      :: i, j, k, ix, kx
    logical                      :: insert

    real(wip),        pointer    :: coord(:,:,:), velocity(:,:,:)
    real(wp),         pointer    :: charge(:,:)
    real(wip),        pointer    :: mass(:,:)
    real(wip),        pointer    :: buf_real(:,:,:)
    integer,          pointer    :: ncel_local, nsolute(:)
    integer,          pointer    :: atmcls(:,:), id_l2g(:,:), id_l2g_sol(:,:)
    integer(int2),    pointer    :: id_g2l(:,:)
    integer,          pointer    :: ptl_add(:), ptl_exit(:), ptlindex(:,:)
    integer,          pointer    :: buf_int(:,:,:)


    nsolute    => constraints%No_HGr

    ncel_local => domain%num_cell_local
    coord      => domain%coord
    velocity   => domain%velocity
    charge     => domain%charge
    mass       => domain%mass
    atmcls     => domain%atom_cls_no
    id_l2g     => domain%id_l2g
    id_l2g_sol => domain%id_l2g_solute
    id_g2l     => domain%id_g2l
    ptl_add    => domain%ptl_add
    ptl_exit   => domain%ptl_exit
    ptlindex   => domain%ptl_exit_index
    buf_int    => domain%buf_integer
    buf_real   => domain%buf_real


    ! incoming particles
    !
    do i = 1, ncel_local

#ifdef DEBUG
      if (ptl_add(i)+nsolute(i)-ptl_exit(i) > MaxAtom) &
        call error_msg('Debug: Update_Incoming_Atom> atom is exceed MaxAtom')
#endif

      ! when the number of coming particles is larger than that of outgoing ones
      !
      if (ptl_add(i) >= ptl_exit(i)) then

        do k = 1, ptl_exit(i)
          coord(1,ptlindex(k,i),i)    = buf_real(1,k,i)
          coord(2,ptlindex(k,i),i)    = buf_real(2,k,i)
          coord(3,ptlindex(k,i),i)    = buf_real(3,k,i)
          velocity(1,ptlindex(k,i),i) = buf_real(4,k,i)
          velocity(2,ptlindex(k,i),i) = buf_real(5,k,i)
          velocity(3,ptlindex(k,i),i) = buf_real(6,k,i)
          charge(ptlindex(k,i),i)     = buf_real(7,k,i)
          mass(ptlindex(k,i),i)       = buf_real(8,k,i)
          atmcls(ptlindex(k,i),i)     = buf_int(1,k,i)
          id_l2g(ptlindex(k,i),i)     = buf_int(2,k,i)
          id_l2g_sol(ptlindex(k,i),i) = buf_int(3,k,i)
          id_g2l(1,buf_int(3,k,i))    = i
          id_g2l(2,buf_int(3,k,i))    = ptlindex(k,i)
        end do

        do k = ptl_exit(i)+1, ptl_add(i)
          ix = k + nsolute(i) - ptl_exit(i)
          coord(1,ix,i)               = buf_real(1,k,i)
          coord(2,ix,i)               = buf_real(2,k,i)
          coord(3,ix,i)               = buf_real(3,k,i)
          velocity(1,ix,i)            = buf_real(4,k,i)
          velocity(2,ix,i)            = buf_real(5,k,i)
          velocity(3,ix,i)            = buf_real(6,k,i)
          charge(ix,i)                = buf_real(7,k,i)
          mass(ix,i)                  = buf_real(8,k,i)
          atmcls(ix,i)                = buf_int (1,k,i)
          id_l2g(ix,i)                = buf_int (2,k,i)
          id_l2g_sol(ix,i)            = buf_int (3,k,i)
          id_g2l(1,buf_int(3,k,i))    = i
          id_g2l(2,buf_int(3,k,i))    = ix
        end do
        nsolute(i) = nsolute(i) + ptl_add(i) - ptl_exit(i)

      ! when the number of coming particles is less than that of outgoing ones
      !
      else

        do k = 1, ptl_add(i)
          coord(1,ptlindex(k,i),i)    = buf_real(1,k,i)
          coord(2,ptlindex(k,i),i)    = buf_real(2,k,i)
          coord(3,ptlindex(k,i),i)    = buf_real(3,k,i)
          velocity(1,ptlindex(k,i),i) = buf_real(4,k,i)
          velocity(2,ptlindex(k,i),i) = buf_real(5,k,i)
          velocity(3,ptlindex(k,i),i) = buf_real(6,k,i)
          charge(ptlindex(k,i),i)     = buf_real(7,k,i)
          mass(ptlindex(k,i),i)       = buf_real(8,k,i)
          atmcls(ptlindex(k,i),i)     = buf_int (1,k,i)
          id_l2g(ptlindex(k,i),i)     = buf_int (2,k,i)
          id_l2g_sol(ptlindex(k,i),i) = buf_int (3,k,i)
          id_g2l(1,buf_int(3,k,i))    = i
          id_g2l(2,buf_int(3,k,i))    = ptlindex(k,i)
        end do

        j  = 0
        ix = nsolute(i)
        k  = ptl_add(i) + 1

        do while (j < (ptl_exit(i)-ptl_add(i)))

          insert = .true.

          do kx = k, ptl_exit(i)
            if (ix == ptlindex(kx,i)) then
              insert = .false.
              j = j + 1
              exit
            end if
          end do

          if (insert) then
            kx = ptlindex(k,i)
            coord(1,kx,i)          = coord(1,ix,i)
            coord(2,kx,i)          = coord(2,ix,i)
            coord(3,kx,i)          = coord(3,ix,i)
            velocity(1,kx,i)       = velocity(1,ix,i)
            velocity(2,kx,i)       = velocity(2,ix,i)
            velocity(3,kx,i)       = velocity(3,ix,i)
            charge(kx,i)           = charge(ix,i)
            mass(kx,i)             = mass(ix,i)
            atmcls(kx,i)           = atmcls(ix,i)
            id_l2g(kx,i)           = id_l2g(ix,i)
            id_l2g_sol(kx,i)       = id_l2g_sol(ix,i)
            id_g2l(1,id_l2g_sol(kx,i)) = i
            id_g2l(2,id_l2g_sol(kx,i)) = kx

            j = j + 1
            k = k + 1

          end if
          ix = ix - 1

        end do
        nsolute(i) = nsolute(i) + ptl_add(i) - ptl_exit(i)

      end if
    end do

    return

  end subroutine update_incoming_atom

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_water
  !> @brief        check water particles incoming cells
  !! @authors      JJ
  !! @param[in]    water_atom : number of water atoms
  !! @param[inout] domain     : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_incoming_water(water_atom, domain)

    ! formal arguments
    integer,                 intent(in)    :: water_atom
    type(s_domain),  target, intent(inout) :: domain

    ! local variable
    integer                  :: i, k, iwater, jwater, ix

    real(wip),       pointer :: coord(:,:,:), velocity(:,:,:)
    real(wip),       pointer :: water_move_real(:,:,:), water_stay_real(:,:,:)
    real(wp),        pointer :: water_charge(:), water_mass(:)
    real(wp),        pointer :: charge(:,:)
    real(wip),       pointer :: mass(:,:)
    integer,         pointer :: ncel_local
    integer,         pointer :: nsolute(:), natom(:), nwater(:)
    integer,         pointer :: water_list(:,:,:)
    integer,         pointer :: id_l2g(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: water_move(:), water_stay(:)
    integer,         pointer :: water_move_int(:,:,:), water_stay_int(:,:,:)
    integer,         pointer :: water_atmcls(:), atmcls(:,:)


    ncel_local      => domain%num_cell_local
    nsolute         => domain%num_solute
    nwater          => domain%num_water
    natom           => domain%num_atom
    water_list      => domain%water_list
    coord           => domain%coord
    velocity        => domain%velocity
    charge          => domain%charge
    mass            => domain%mass
    atmcls          => domain%atom_cls_no
    id_l2g          => domain%id_l2g
    id_g2l          => domain%id_g2l
    water_atmcls    => domain%water%atom_cls_no
    water_charge    => domain%water%charge
    water_mass      => domain%water%mass
    water_move      => domain%water%move
    water_stay      => domain%water%stay
    water_move_real => domain%water%move_real
    water_stay_real => domain%water%stay_real
    water_move_int  => domain%water%move_integer
    water_stay_int  => domain%water%stay_integer


    do i = 1, ncel_local

#ifdef DEBUG
      if (water_stay(i)+water_move(i) > MaxWater) then
        call error_msg('Debug: Update_Incoming_Water> water is exceed MaxWater')
      end if
#endif

      do iwater = 1, water_stay(i)
        do k = 1, water_atom
          water_list(k,iwater,i) = nsolute(i) + water_atom*(iwater-1) + k
          ix = water_list(k,iwater,i)

          coord(1,ix,i)          = water_stay_real(6*(k-1)+1,iwater,i)
          coord(2,ix,i)          = water_stay_real(6*(k-1)+2,iwater,i)
          coord(3,ix,i)          = water_stay_real(6*(k-1)+3,iwater,i)
          velocity(1,ix,i)       = water_stay_real(6*(k-1)+4,iwater,i)
          velocity(2,ix,i)       = water_stay_real(6*(k-1)+5,iwater,i)
          velocity(3,ix,i)       = water_stay_real(6*(k-1)+6,iwater,i)
          id_l2g(ix,i)           = water_stay_int (k,iwater,i)
          charge(ix,i)           = water_charge(k)
          mass(ix,i)             = water_mass(k)
          atmcls(ix,i)           = water_atmcls(k)
        end do
      end do

      do iwater = 1, water_move(i)
        jwater = iwater + water_stay(i)
        do k = 1, water_atom
          water_list(k,jwater,i) = nsolute(i) + water_atom*(jwater-1) + k
          ix = water_list(k,jwater,i)
          coord(1,ix,i)          = water_move_real(6*(k-1)+1,iwater,i)
          coord(2,ix,i)          = water_move_real(6*(k-1)+2,iwater,i)
          coord(3,ix,i)          = water_move_real(6*(k-1)+3,iwater,i)
          velocity(1,ix,i)       = water_move_real(6*(k-1)+4,iwater,i)
          velocity(2,ix,i)       = water_move_real(6*(k-1)+5,iwater,i)
          velocity(3,ix,i)       = water_move_real(6*(k-1)+6,iwater,i)
          id_l2g(ix,i)           = water_move_int (k,iwater,i)
          charge(ix,i)           = water_charge(k)
          mass(ix,i)             = water_mass(k)
          atmcls(ix,i)           = water_atmcls(k)
        end do
      end do

      nwater(i) = water_stay(i) + water_move(i)
      natom(i)  = nsolute(i) + water_atom*nwater(i)
    end do

    return

  end subroutine update_incoming_water

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_HGr
  !> @brief        check hydrogen-bonded particles incoming cells
  !! @authors      JJ
  !! @param[inout] constraints : constraints information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_incoming_HGr(constraints, domain)

    ! formal variable
    type(s_constraints), target, intent(inout)  :: constraints
    type(s_domain),      target, intent(inout)  :: domain

    ! local variable
    integer                      :: i, j, k, k1, list, ix, connect
    integer                      :: num_atom

    real(wip),        pointer    :: coord(:,:,:), vel(:,:,:)
    real(wip),        pointer    :: hgr_move_real(:,:,:,:)
    real(wip),        pointer    :: hgr_stay_real(:,:,:,:)
    real(wp),         pointer    :: charge(:,:)
    real(wip),        pointer    :: mass(:,:)
    real(wip),        pointer    :: HGr_bond_dist(:,:,:,:)
    integer,          pointer    :: ncel_local
    integer,          pointer    :: nsolute(:)
    integer,          pointer    :: id_l2g(:,:), id_l2g_sol(:,:), atmcls(:,:)
    integer(int2),    pointer    :: id_g2l(:,:)
    integer,          pointer    :: No_HGr(:), HGr_local(:,:)
    integer,          pointer    :: HGr_bond_list(:,:,:,:)
    integer,          pointer    :: HGr_move(:,:), HGr_stay(:,:)
    integer,          pointer    :: HGr_move_int(:,:,:,:), HGr_stay_int(:,:,:,:)


    No_HGr        => constraints%no_hgr
    HGr_local     => constraints%HGr_local
    HGr_bond_list => constraints%HGr_bond_list
    HGr_move      => constraints%HGr_move
    HGr_stay      => constraints%HGr_stay
    HGr_bond_dist => constraints%HGr_bond_dist
    HGr_move_real => constraints%HGr_move_real
    HGr_stay_real => constraints%HGr_stay_real
    HGr_move_int  => constraints%HGr_move_int
    HGr_stay_int  => constraints%HGr_stay_int
    connect       =  constraints%connect

    ncel_local    => domain%num_cell_local
    nsolute       => domain%num_solute
    coord         => domain%coord
    vel           => domain%velocity
    charge        => domain%charge
    mass          => domain%mass
    atmcls        => domain%atom_cls_no
    id_l2g        => domain%id_l2g
    id_l2g_sol    => domain%id_l2g_solute
    id_g2l        => domain%id_g2l


    do i = 1, ncel_local

      num_atom = No_HGr(i)

      do j = 1, connect

#ifdef DEBUG
        if (HGr_stay(j,i)+Hgr_move(j,i) > HGroupMax) then
          call error_msg('Debug: Update_Incoming_HGr> HGr is exceed HGroupMax')
        end if
#endif

        do k = 1, HGr_stay(j,i)
          do list = 1, j+1
            num_atom = num_atom + 1
            HGr_bond_list(list,k,j,i) = num_atom
            ix = num_atom
            coord(1,ix,i)             = HGr_stay_real(9*(list-1)+1,k,j,i)
            coord(2,ix,i)             = HGr_stay_real(9*(list-1)+2,k,j,i)
            coord(3,ix,i)             = HGr_stay_real(9*(list-1)+3,k,j,i)
            vel(1,ix,i)               = HGr_stay_real(9*(list-1)+4,k,j,i)
            vel(2,ix,i)               = HGr_stay_real(9*(list-1)+5,k,j,i)
            vel(3,ix,i)               = HGr_stay_real(9*(list-1)+6,k,j,i)
            charge(ix,i)              = HGr_stay_real(9*(list-1)+7,k,j,i)
            mass(ix,i)                = HGr_stay_real(9*(list-1)+8,k,j,i)
            HGr_bond_dist(list,k,j,i) = HGr_stay_real(9*(list-1)+9,k,j,i)
            atmcls(ix,i)              = HGr_stay_int (3*(list-1)+1,k,j,i)
            id_l2g(ix,i)              = HGr_stay_int (3*(list-1)+2,k,j,i)
            id_l2g_sol(ix,i)          = HGr_stay_int (3*(list-1)+3,k,j,i)
            id_g2l(1,id_l2g_sol(ix,i))= i
            id_g2l(2,id_l2g_sol(ix,i))= ix
          end do
        end do

        do k = 1, HGr_move(j,i)
          k1 = k + HGr_stay(j,i)
          do list = 1, j+1
            num_atom = num_atom + 1
            HGr_bond_list(list,k1,j,i) = num_atom
            ix = num_atom
            coord(1,ix,i)             = HGr_move_real(9*(list-1)+1,k,j,i)
            coord(2,ix,i)             = HGr_move_real(9*(list-1)+2,k,j,i)
            coord(3,ix,i)             = HGr_move_real(9*(list-1)+3,k,j,i)
            vel(1,ix,i)               = HGr_move_real(9*(list-1)+4,k,j,i)
            vel(2,ix,i)               = HGr_move_real(9*(list-1)+5,k,j,i)
            vel(3,ix,i)               = HGr_move_real(9*(list-1)+6,k,j,i)
            charge(ix,i)              = HGr_move_real(9*(list-1)+7,k,j,i)
            mass(ix,i)                = HGr_move_real(9*(list-1)+8,k,j,i)
            HGr_bond_dist(list,k1,j,i)= HGr_move_real(9*(list-1)+9,k,j,i)
            atmcls(ix,i)              = HGr_move_int (3*(list-1)+1,k,j,i)
            id_l2g(ix,i)              = HGr_move_int (3*(list-1)+2,k,j,i)
            id_l2g_sol(ix,i)          = HGr_move_int (3*(list-1)+3,k,j,i)
            id_g2l(1,id_l2g_sol(ix,i))= i
            id_g2l(2,id_l2g_sol(ix,i))= ix
          end do
        end do

        HGr_local(j,i) = HGr_stay(j,i) + HGr_move(j,i)

      end do

      nsolute(i) = num_atom

    end do

    return

  end subroutine update_incoming_HGr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_enefunc_enm
  !> @brief        update ENM term for each cell in potential energy
  !function
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_enefunc_enm(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    integer                  :: i, ix, icel_local, k
    integer                  :: icel1, icel2
    integer                  :: list1, list2

    real(wp),        pointer :: dist(:,:), force(:,:)
    real(wp),        pointer :: buf_real(:,:,:)
    integer,         pointer :: ncel_local, nboundary
    integer(int2),   pointer :: cell_pair(:,:), id_g2l(:,:)
    integer,         pointer :: enm(:), list(:,:,:)
    integer,         pointer :: enm_add(:), buf_int(:,:,:)


    ncel_local  => domain%num_cell_local
    nboundary   => domain%num_cell_boundary
    cell_pair   => domain%cell_pair
    id_g2l      => domain%id_g2l

    enm         => enefunc%num_enm
    list        => enefunc%enm_list
    force       => enefunc%enm_force_const
    dist        => enefunc%enm_dist_min
    enm_add     => enefunc%enm_add
    buf_int     => enefunc%buf_enm_integer
    buf_real    => enefunc%buf_enm_real

    !$omp parallel do private(i, ix)
    do i = 1, ncel_local+nboundary
      do ix = 1, enm(i)
        buf_int (1,ix,i) = list(1,ix,i)
        buf_int (2,ix,i) = list(2,ix,i)
        buf_real(1,ix,i) = force(ix,i)
        buf_real(2,ix,i) = dist(ix,i)
      end do
      enm_add(i) = enm(i)
    end do
    !$omp end parallel do

    enm(1:ncel_local+nboundary) = 0

    do i = 1, ncel_local+nboundary

      do ix = 1, enm_add(i)

        icel1 = id_g2l(1,buf_int(1,ix,i))
        icel2 = id_g2l(1,buf_int(2,ix,i))

        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1,icel2)

          if (icel_local > 0 .and. icel_local <= ncel_local) then
            enm(icel_local) = enm(icel_local) + 1
            k = enm(icel_local)
            force(k,icel_local)  = buf_real(1,ix,i)
            dist (k,icel_local)  = buf_real(2,ix,i)
            list(1,k,icel_local) = buf_int(1,ix,i)
            list(2,k,icel_local) = buf_int(2,ix,i)
          end if

        end if

      end do
    end do

    return

  end subroutine update_enefunc_enm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_enefunc_bond
  !> @brief        update BOND term for each cell in potential energy function
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_enefunc_bond(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    integer                  :: i, ix, icel_local
    integer                  :: i1, i2, icel1, icel2

    real(wp),        pointer :: fc(:,:), r0(:,:)
    real(wp),        pointer :: buf_real(:,:,:)
    integer,         pointer :: ncel_local, nboundary
    integer(int2),   pointer :: cell_pair(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: nbond(:), bondlist(:,:,:)
    integer(1),      pointer :: bondkind(:,:)
    integer,         pointer :: bond_pbc(:,:)
    integer,         pointer :: bond_exit(:), index(:,:)
    integer,         pointer :: bond_add(:), buf_int(:,:,:)


    ncel_local  => domain%num_cell_local
    nboundary   => domain%num_cell_boundary
    cell_pair   => domain%cell_pair
    id_g2l      => domain%id_g2l

    nbond       => enefunc%num_bond
    bondlist    => enefunc%bond_list
    bondkind    => enefunc%bond_kind
    bond_pbc    => enefunc%bond_pbc
    fc          => enefunc%bond_force_const
    r0          => enefunc%bond_dist_min
    bond_exit   => enefunc%bond_exit
    index       => enefunc%bond_exit_index
    bond_add    => enefunc%bond_add
    buf_int     => enefunc%buf_bond_integer
    buf_real    => enefunc%buf_bond_real

    bond_exit(1:ncel_local) = 0
    bond_add(1:ncel_local+nboundary) = 0

    do i = 1, ncel_local
      do ix = 1, nbond(i)

        icel1 = id_g2l(1,bondlist(1,ix,i))
        i1    = id_g2l(2,bondlist(1,ix,i))
        icel2 = id_g2l(1,bondlist(2,ix,i))
        i2    = id_g2l(2,bondlist(2,ix,i))

        icel_local = cell_pair(icel1,icel2)

        if (icel_local /= i) then

          bond_exit(i) = bond_exit(i) + 1
          index(bond_exit(i),i) = ix

          bond_add(icel_local) = bond_add(icel_local) + 1
          buf_real(1,bond_add(icel_local),icel_local) = r0(ix,i)
          buf_real(2,bond_add(icel_local),icel_local) = fc(ix,i)
          buf_int (1,bond_add(icel_local),icel_local) = bondlist(1,ix,i)
          buf_int (2,bond_add(icel_local),icel_local) = bondlist(2,ix,i)
          buf_int (3,bond_add(icel_local),icel_local) = bondkind(ix,i)
          buf_int (4,bond_add(icel_local),icel_local) = bond_pbc(ix,i)

        end if

      end do
    end do

    return

  end subroutine update_outgoing_enefunc_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_enefunc_bond
  !> @brief        update BOND term for each cell in potential energy function
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_incoming_enefunc_bond(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    integer                  :: i, j, k, ix, kx, found
    logical                  :: insert

    real(wp),        pointer :: buf_real(:,:,:), bonddist(:,:), bondforce(:,:)
    integer,         pointer :: ncel_local, nboundary
    integer(int2),   pointer :: cell_g2l(:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: bond_add(:), bond_exit(:), bond_exit_index(:,:)
    integer,         pointer :: buf_int(:,:,:), nbond(:)
    integer,         pointer :: bondlist(:,:,:)
    integer(1),      pointer :: bondkind(:,:)
    integer,         pointer :: bond_pbc(:,:)


    ncel_local      => domain%num_cell_local
    nboundary       => domain%num_cell_boundary
    cell_g2l        => domain%cell_g2l
    id_g2l          => domain%id_g2l

    bond_add        => enefunc%bond_add
    bond_exit       => enefunc%bond_exit
    bond_exit_index => enefunc%bond_exit_index
    buf_int         => enefunc%buf_bond_integer
    buf_real        => enefunc%buf_bond_real
    bondlist        => enefunc%bond_list
    bondforce       => enefunc%bond_force_const
    bonddist        => enefunc%bond_dist_min
    nbond           => enefunc%num_bond
    bondkind        => enefunc%bond_kind
    bond_pbc        => enefunc%bond_pbc

    found = 0
    do i = 1, ncel_local

#ifdef DEBUG
      if (bond_add(i)+nbond(i)-bond_exit(i) > MaxBond) then
        call error_msg( &
            'Debug: update_incoming_enefunc_bond> bond is exceed MaxBond')
      end if
#endif

      if (bond_add(i) > BondMove) &
        call error_msg( &
        'Debug: update_incoming_enefunc_bond> bond migration is exceed maximum')

      if (bond_add(i) >= bond_exit(i)) then

        do k = 1, bond_exit(i)
          bondlist (1,bond_exit_index(k,i),i) = buf_int (1,k,i)
          bondlist (2,bond_exit_index(k,i),i) = buf_int (2,k,i)
          bondkind(bond_exit_index(k,i),i)    = buf_int (3,k,i)
          bond_pbc(bond_exit_index(k,i),i)    = buf_int (4,k,i)
          bonddist (bond_exit_index(k,i),i)   = buf_real(1,k,i)
          bondforce(bond_exit_index(k,i),i)   = buf_real(2,k,i)
        end do

        do k = bond_exit(i)+1, bond_add(i)
          ix = k + nbond(i) - bond_exit(i)
          bondlist (1,ix,i) = buf_int(1,k,i)
          bondlist (2,ix,i) = buf_int(2,k,i)
          bondkind (ix,i)   = buf_int(3,k,i)
          bond_pbc (ix,i)   = buf_int(4,k,i)
          bonddist (ix,i)   = buf_real(1,k,i)
          bondforce(ix,i)   = buf_real(2,k,i)
        end do

        nbond(i) = nbond(i) + bond_add(i) - bond_exit(i)

      else

        do k = 1, bond_add(i)
          bondlist (1,bond_exit_index(k,i),i) = buf_int (1,k,i)
          bondlist (2,bond_exit_index(k,i),i) = buf_int (2,k,i)
          bondkind(bond_exit_index(k,i),i)    = buf_int (3,k,i)
          bond_pbc(bond_exit_index(k,i),i)    = buf_int (4,k,i)
          bonddist (bond_exit_index(k,i),i)   = buf_real(1,k,i)
          bondforce(bond_exit_index(k,i),i)   = buf_real(2,k,i)
        end do

        j  = 0
        ix = nbond(i)
        k  = bond_add(i) + 1

        do while (j < (bond_exit(i)-bond_add(i)))

          insert = .true.
          do kx = k, bond_exit(i)
            if (ix == bond_exit_index(kx,i)) then
              insert = .false.
              j = j + 1
              exit
            end if
          end do

          if (insert) then
            kx = bond_exit_index(k,i)
            bondlist (1,kx,i) = bondlist(1,ix,i)
            bondlist (2,kx,i) = bondlist(2,ix,i)
            bondkind (kx,i)   = bondkind(ix,i)
            bond_pbc (kx,i)   = bond_pbc(ix,i)
            bonddist (kx,i)   = bonddist(ix,i)
            bondforce(kx,i)   = bondforce(ix,i)
            j = j + 1
            k = k + 1
          end if
          ix = ix - 1

        end do

        nbond(i) = nbond(i) + bond_add(i) - bond_exit(i)

      end if

    end do

    return

  end subroutine update_incoming_enefunc_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_enefunc_angl
  !> @brief        update ANGLE term for each cell in potential energy function
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_enefunc_angl(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: icel1, icel2, icel3, i1, i2, i3
    integer                  :: i, ix, icel_local

    real(wp),        pointer :: fc(:,:), theta0(:,:), fc_ub(:,:), r0_ub(:,:)
    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: buf_real(:,:,:)
    integer,         pointer :: ncel_local, nboundary
    integer(int2),   pointer :: cell_g2l(:), cell_g2b(:)
    integer(int2),   pointer :: cell_pair(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: nangle(:), anglelist(:,:,:), angl_pbc(:,:,:)
    integer(1),      pointer :: anglekind(:,:)
    integer,         pointer :: angle_add(:), angle_exit(:)
    integer,         pointer :: angle_exit_index(:,:), buf_int(:,:,:)


    ncel_local       => domain%num_cell_local
    nboundary        => domain%num_cell_boundary
    cell_g2l         => domain%cell_g2l
    cell_g2b         => domain%cell_g2b
    cell_pair        => domain%cell_pair
    id_g2l           => domain%id_g2l
    coord            => domain%coord

    nangle           => enefunc%num_angle
    anglelist        => enefunc%angle_list
    anglekind        => enefunc%angle_kind
    angl_pbc         => enefunc%angle_pbc
    fc               => enefunc%angle_force_const
    theta0           => enefunc%angle_theta_min
    fc_ub            => enefunc%urey_force_const
    r0_ub            => enefunc%urey_rmin
    angle_add        => enefunc%angle_add
    angle_exit       => enefunc%angle_exit
    angle_exit_index => enefunc%angle_exit_index
    buf_int          => enefunc%buf_angle_integer
    buf_real         => enefunc%buf_angle_real

    angle_add (1:ncel_local+nboundary) = 0
    angle_exit(1:ncel_local) = 0

    do i = 1, ncel_local
      do ix = 1, nangle(i)

        icel1 = id_g2l(1,anglelist(1,ix,i))
        i1    = id_g2l(2,anglelist(1,ix,i))
        icel2 = id_g2l(1,anglelist(2,ix,i))
        i2    = id_g2l(2,anglelist(2,ix,i))
        icel3 = id_g2l(1,anglelist(3,ix,i))
        i3    = id_g2l(2,anglelist(3,ix,i))

        icel_local = cell_pair(icel1,icel3)

        if (icel_local /= i) then

          angle_exit(i) = angle_exit(i) + 1
          angle_exit_index(angle_exit(i),i) = ix
          angle_add(icel_local) = angle_add(icel_local) + 1
          buf_real(1,angle_add(icel_local),icel_local) = theta0(ix,i)
          buf_real(2,angle_add(icel_local),icel_local) = fc(ix,i)
          buf_real(3,angle_add(icel_local),icel_local) = r0_ub(ix,i)
          buf_real(4,angle_add(icel_local),icel_local) = fc_ub(ix,i)
          buf_int (1,angle_add(icel_local),icel_local) = anglelist(1,ix,i)
          buf_int (2,angle_add(icel_local),icel_local) = anglelist(2,ix,i)
          buf_int (3,angle_add(icel_local),icel_local) = anglelist(3,ix,i)
          buf_int (4,angle_add(icel_local),icel_local) = anglekind(ix,i)
          buf_int (5,angle_add(icel_local),icel_local) = angl_pbc(1,ix,i)
          buf_int (6,angle_add(icel_local),icel_local) = angl_pbc(2,ix,i)
          buf_int (7,angle_add(icel_local),icel_local) = angl_pbc(3,ix,i)
         

        end if

      end  do
    end do

    return

  end subroutine update_outgoing_enefunc_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_enefunc_angl
  !> @brief        update ANGLE term for each cell in potential energy function
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_incoming_enefunc_angl(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    integer                  :: i, j, k, ix, kx
    logical                  :: insert

    real(wp),        pointer :: buf_real(:,:,:)
    real(wp),        pointer :: angletheta(:,:), angleforce(:,:)
    real(wp),        pointer :: ubrmin(:,:), ubforce(:,:)
    integer,         pointer :: ncel_local, nboundary
    integer(int2),   pointer :: cell_g2l(:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: angle_add(:), angle_exit(:)
    integer,         pointer :: angle_exit_index(:,:)
    integer(1),      pointer :: anglekind(:,:)
    integer,         pointer :: buf_int(:,:,:), nangle(:), anglelist(:,:,:)
    integer,         pointer :: angl_pbc(:,:,:)

    ncel_local       => domain%num_cell_local
    nboundary        => domain%num_cell_boundary
    cell_g2l         => domain%cell_g2l
    id_g2l           => domain%id_g2l

    nangle           => enefunc%num_angle
    angle_add        => enefunc%angle_add
    angle_exit       => enefunc%angle_exit
    angle_exit_index => enefunc%angle_exit_index
    buf_int          => enefunc%buf_angle_integer
    buf_real         => enefunc%buf_angle_real
    anglelist        => enefunc%angle_list
    anglekind        => enefunc%angle_kind
    angl_pbc         => enefunc%angle_pbc 
    angleforce       => enefunc%angle_force_const
    angletheta       => enefunc%angle_theta_min
    ubforce          => enefunc%urey_force_const
    ubrmin           => enefunc%urey_rmin

    do i = 1, ncel_local

#ifdef DEBUG
      if (angle_add(i)+nangle(i)-angle_exit(i) > MaxAngle) &
        call error_msg( &
            'Debug: update_incoming_enefunc_angle> angle is exceed MaxAngle')
#endif

      if (angle_add(i) > AngleMove) &
        call error_msg( &
     'Debug: update_incoming_enefunc_angle> angle migration is exceed maximum')

      if (angle_add(i) >= angle_exit(i)) then

        do k = 1, angle_exit(i)
          anglelist (1,angle_exit_index(k,i),i) = buf_int (1,k,i)
          anglelist (2,angle_exit_index(k,i),i) = buf_int (2,k,i)
          anglelist (3,angle_exit_index(k,i),i) = buf_int (3,k,i)
          anglekind (angle_exit_index(k,i),i)   = buf_int (4,k,i)
          angl_pbc  (1,angle_exit_index(k,i),i) = buf_int (5,k,i)
          angl_pbc  (2,angle_exit_index(k,i),i) = buf_int (6,k,i)
          angl_pbc  (3,angle_exit_index(k,i),i) = buf_int (7,k,i)
          angletheta(angle_exit_index(k,i),i)   = buf_real(1,k,i)
          angleforce(angle_exit_index(k,i),i)   = buf_real(2,k,i)
          ubrmin    (angle_exit_index(k,i),i)   = buf_real(3,k,i)
          ubforce   (angle_exit_index(k,i),i)   = buf_real(4,k,i)
        end do

        do k = angle_exit(i)+1, angle_add(i)
          ix = k + nangle(i) - angle_exit(i)
          anglelist (1,ix,i) = buf_int(1,k,i)
          anglelist (2,ix,i) = buf_int(2,k,i)
          anglelist (3,ix,i) = buf_int(3,k,i)
          anglekind (ix,i)   = buf_int(4,k,i)
          angl_pbc  (1,ix,i) = buf_int(5,k,i)
          angl_pbc  (2,ix,i) = buf_int(6,k,i)
          angl_pbc  (3,ix,i) = buf_int(7,k,i)
          angletheta(ix,i)   = buf_real(1,k,i)
          angleforce(ix,i)   = buf_real(2,k,i)
          ubrmin    (ix,i)   = buf_real(3,k,i)
          ubforce   (ix,i)   = buf_real(4,k,i)
        end do
        nangle(i) = nangle(i) + angle_add(i) - angle_exit(i)

      else
        do k = 1, angle_add(i)
          anglelist (1,angle_exit_index(k,i),i) = buf_int (1,k,i)
          anglelist (2,angle_exit_index(k,i),i) = buf_int (2,k,i)
          anglelist (3,angle_exit_index(k,i),i) = buf_int (3,k,i)
          anglekind (angle_exit_index(k,i),i)   = buf_int (4,k,i)
          angl_pbc  (1,angle_exit_index(k,i),i) = buf_int (5,k,i)
          angl_pbc  (2,angle_exit_index(k,i),i) = buf_int (6,k,i)
          angl_pbc  (3,angle_exit_index(k,i),i) = buf_int (7,k,i)
          angletheta(angle_exit_index(k,i),i)   = buf_real(1,k,i)
          angleforce(angle_exit_index(k,i),i)   = buf_real(2,k,i)
          ubrmin    (angle_exit_index(k,i),i)   = buf_real(3,k,i)
          ubforce   (angle_exit_index(k,i),i)   = buf_real(4,k,i)
        end do

        j  = 0
        ix = nangle(i)
        k  = angle_add(i) + 1

        do while (j < (angle_exit(i)-angle_add(i)))

          insert = .true.

          do kx = k, angle_exit(i)
            if (ix == angle_exit_index(kx,i)) then
              insert = .false.
              j = j + 1
              exit
            end if
          end do

          if (insert) then
            kx = angle_exit_index(k,i)
            anglelist (1,kx,i) = anglelist(1,ix,i)
            anglelist (2,kx,i) = anglelist(2,ix,i)
            anglelist (3,kx,i) = anglelist(3,ix,i)
            anglekind (kx,i)   = anglekind(ix,i)
            angl_pbc  (1,kx,i) = angl_pbc (1,ix,i)
            angl_pbc  (2,kx,i) = angl_pbc (2,ix,i)
            angl_pbc  (3,kx,i) = angl_pbc (3,ix,i)
            angletheta(kx,i)   = angletheta(ix,i)
            angleforce(kx,i)   = angleforce(ix,i)
            ubrmin    (kx,i)   = ubrmin(ix,i)
            ubforce   (kx,i)   = ubforce(ix,i)
            j = j + 1
            k = k + 1

          end if
          ix = ix - 1

        end do

        nangle(i) = nangle(i) + angle_add(i) - angle_exit(i)

      end if

    end do

    return

  end subroutine update_incoming_enefunc_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_enefunc_dihe
  !> @brief        update DIHEDRAL term for each cell in potential energy
  !!               function
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_enefunc_dihe(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: icel1, icel2, icel3, icel4, i1, i2, i3, i4
    integer                  :: i, ix, icel_local

    real(wp),        pointer :: fc(:,:), phase(:,:)
    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: buf_real(:,:,:)
    integer,         pointer :: ncel_local, nboundary
    integer(int2),   pointer :: cell_g2l(:), cell_g2b(:)
    integer(int2),   pointer :: cell_pair(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: dihed_add(:), dihed_exit(:)
    integer,         pointer :: dihed_exit_index(:,:)
    integer,         pointer :: ndihe(:), dihelist(:,:,:), dihe_pbc(:,:,:)
    integer(1),      pointer :: dihekind(:,:)
    integer,         pointer :: nperiod(:,:), buf_int(:,:,:)


    ncel_local       => domain%num_cell_local
    nboundary        => domain%num_cell_boundary
    cell_g2l         => domain%cell_g2l
    cell_g2b         => domain%cell_g2b
    cell_pair        => domain%cell_pair
    id_g2l           => domain%id_g2l
    coord            => domain%coord

    ndihe            => enefunc%num_dihedral
    dihelist         => enefunc%dihe_list
    dihekind         => enefunc%dihe_kind
    dihe_pbc         => enefunc%dihe_pbc 
    dihed_add        => enefunc%dihed_add
    dihed_exit       => enefunc%dihed_exit
    dihed_exit_index => enefunc%dihed_exit_index
    fc               => enefunc%dihe_force_const
    nperiod          => enefunc%dihe_periodicity
    phase            => enefunc%dihe_phase
    buf_int          => enefunc%buf_dihed_integer
    buf_real         => enefunc%buf_dihed_real

    dihed_add (1:ncel_local+nboundary) = 0
    dihed_exit(1:ncel_local) = 0


    do i = 1, ncel_local
      do ix = 1, ndihe(i)

        icel1 = id_g2l(1,dihelist(1,ix,i))
        i1    = id_g2l(2,dihelist(1,ix,i))
        icel2 = id_g2l(1,dihelist(2,ix,i))
        i2    = id_g2l(2,dihelist(2,ix,i))
        icel3 = id_g2l(1,dihelist(3,ix,i))
        i3    = id_g2l(2,dihelist(3,ix,i))
        icel4 = id_g2l(1,dihelist(4,ix,i))
        i4    = id_g2l(2,dihelist(4,ix,i))

        icel_local = cell_pair(icel1,icel4)

        if (icel_local /= i) then

          dihed_exit(i) = dihed_exit(i) + 1
          dihed_exit_index(dihed_exit(i),i) = ix
          dihed_add(icel_local) = dihed_add(icel_local) + 1
          buf_real(1,dihed_add(icel_local),icel_local) = phase(ix,i)
          buf_real(2,dihed_add(icel_local),icel_local) = fc(ix,i)
          buf_int (1,dihed_add(icel_local),icel_local) = dihelist(1,ix,i)
          buf_int (2,dihed_add(icel_local),icel_local) = dihelist(2,ix,i)
          buf_int (3,dihed_add(icel_local),icel_local) = dihelist(3,ix,i)
          buf_int (4,dihed_add(icel_local),icel_local) = dihelist(4,ix,i)
          buf_int (5,dihed_add(icel_local),icel_local) = nperiod(ix,i)
          buf_int (6,dihed_add(icel_local),icel_local) = dihekind(ix,i)
          buf_int (7,dihed_add(icel_local),icel_local) = dihe_pbc(1,ix,i)
          buf_int (8,dihed_add(icel_local),icel_local) = dihe_pbc(2,ix,i)
          buf_int (9,dihed_add(icel_local),icel_local) = dihe_pbc(3,ix,i)

        end if

      end  do
    end do

    return

  end subroutine update_outgoing_enefunc_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_enefunc_dihe
  !> @brief        update DIHEDRAL term for each cell in potential energy
  !!               function
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_incoming_enefunc_dihe(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    integer                  :: i, j, k, ix, kx, ndihe
    logical                  :: insert

    real(wp),        pointer :: dihephase(:,:), diheforce(:,:), buf_real(:,:,:)
    integer,         pointer :: ncel_local, nboundary
    integer(int2),   pointer :: cell_g2l(:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: dihed_add(:), dihed_exit(:)
    integer,         pointer :: dihed_exit_index(:,:), buf_int(:,:,:)
    integer,         pointer :: ndihedral(:), dihelist(:,:,:), diheperio(:,:)
    integer(1),      pointer :: dihekind(:,:)
    integer,         pointer :: dihe_pbc(:,:,:)


    ncel_local       => domain%num_cell_local
    nboundary        => domain%num_cell_boundary
    cell_g2l         => domain%cell_g2l
    id_g2l           => domain%id_g2l

    ndihedral        => enefunc%num_dihedral
    dihelist         => enefunc%dihe_list
    dihekind         => enefunc%dihe_kind
    dihe_pbc         => enefunc%dihe_pbc 
    diheforce        => enefunc%dihe_force_const
    diheperio        => enefunc%dihe_periodicity
    dihephase        => enefunc%dihe_phase
    dihed_add        => enefunc%dihed_add
    dihed_exit       => enefunc%dihed_exit
    dihed_exit_index => enefunc%dihed_exit_index
    buf_int          => enefunc%buf_dihed_integer
    buf_real         => enefunc%buf_dihed_real

    ndihe = enefunc%num_dihe_all

    do i = 1, ncel_local

#ifdef DEBUG
      if (dihed_add(i)+ndihedral(i)-dihed_exit(i) > MaxDihe) &
        call error_msg( &
            'Debug: update_incoming_enefunc_dihed> dihed is exceed MaxDihe')
#endif

      if (dihed_add(i) > DiheMove) &
        call error_msg( &
    'Debug: update_incoming_enefunc_dihed> dihed migration is exceed maximum')

      if (dihed_add(i) >= dihed_exit(i)) then

        do k = 1, dihed_exit(i)
          dihelist (1,dihed_exit_index(k,i),i) = buf_int (1,k,i)
          dihelist (2,dihed_exit_index(k,i),i) = buf_int (2,k,i)
          dihelist (3,dihed_exit_index(k,i),i) = buf_int (3,k,i)
          dihelist (4,dihed_exit_index(k,i),i) = buf_int (4,k,i)
          diheperio(dihed_exit_index(k,i),i)   = buf_int (5,k,i)
          dihekind (dihed_exit_index(k,i),i)   = buf_int (6,k,i)
          dihe_pbc (1,dihed_exit_index(k,i),i) = buf_int (7,k,i)
          dihe_pbc (2,dihed_exit_index(k,i),i) = buf_int (8,k,i)
          dihe_pbc (3,dihed_exit_index(k,i),i) = buf_int (9,k,i)
          dihephase(dihed_exit_index(k,i),i)   = buf_real(1,k,i)
          diheforce(dihed_exit_index(k,i),i)   = buf_real(2,k,i)
        end do

        do k = dihed_exit(i)+1, dihed_add(i)
          ix = k + ndihedral(i) - dihed_exit(i)
          dihelist (1,ix,i) = buf_int (1,k,i)
          dihelist (2,ix,i) = buf_int (2,k,i)
          dihelist (3,ix,i) = buf_int (3,k,i)
          dihelist (4,ix,i) = buf_int (4,k,i)
          diheperio(ix,i)   = buf_int (5,k,i)
          dihekind (ix,i)   = buf_int (6,k,i)
          dihe_pbc (1,ix,i) = buf_int (7,k,i)
          dihe_pbc (2,ix,i) = buf_int (8,k,i)
          dihe_pbc (3,ix,i) = buf_int (9,k,i)
          dihephase(ix,i)   = buf_real(1,k,i)
          diheforce(ix,i)   = buf_real(2,k,i)
        end do

        ndihedral(i) = ndihedral(i) + dihed_add(i) - dihed_exit(i)

      else
        do k = 1, dihed_add(i)
          dihelist (1,dihed_exit_index(k,i),i) = buf_int (1,k,i)
          dihelist (2,dihed_exit_index(k,i),i) = buf_int (2,k,i)
          dihelist (3,dihed_exit_index(k,i),i) = buf_int (3,k,i)
          dihelist (4,dihed_exit_index(k,i),i) = buf_int (4,k,i)
          diheperio(dihed_exit_index(k,i),i)   = buf_int (5,k,i)
          dihekind (dihed_exit_index(k,i),i)   = buf_int (6,k,i)
          dihe_pbc (1,dihed_exit_index(k,i),i) = buf_int (7,k,i)
          dihe_pbc (2,dihed_exit_index(k,i),i) = buf_int (8,k,i)
          dihe_pbc (3,dihed_exit_index(k,i),i) = buf_int (9,k,i)
          dihephase(dihed_exit_index(k,i),i)   = buf_real(1,k,i)
          diheforce(dihed_exit_index(k,i),i)   = buf_real(2,k,i)
        end do

        j  = 0
        ix = ndihedral(i)
        k  = dihed_add(i) + 1

        do while (j < (dihed_exit(i)-dihed_add(i)))

          insert = .true.
          do kx = k, dihed_exit(i)
            if (ix == dihed_exit_index(kx,i)) then
              insert = .false.
              j = j + 1
              exit
            end if
          end do

          if (insert) then
            kx = dihed_exit_index(k,i)
            dihelist (1,kx,i) = dihelist (1,ix,i)
            dihelist (2,kx,i) = dihelist (2,ix,i)
            dihelist (3,kx,i) = dihelist (3,ix,i)
            dihelist (4,kx,i) = dihelist (4,ix,i)
            diheperio(kx,i)   = diheperio(ix,i)
            dihekind (kx,i)   = dihekind (ix,i)
            dihe_pbc (1,kx,i) = dihe_pbc (1,ix,i)
            dihe_pbc (2,kx,i) = dihe_pbc (2,ix,i)
            dihe_pbc (3,kx,i) = dihe_pbc (3,ix,i)
            dihephase(kx,i)   = dihephase(ix,i)
            diheforce(kx,i)   = diheforce(ix,i)
            j = j + 1
            k = k + 1
          end if

          ix = ix - 1

        end do
        ndihedral(i) = ndihedral(i) + dihed_add(i) - dihed_exit(i)

      end if

    end do

    return

  end subroutine update_incoming_enefunc_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_enefunc_rb_dihe
  !> @brief        update Ryckaert-Bellemans DIHEDRAL term for each cell in
  !!               potential energy function
  !! @authors      NT
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_enefunc_rb_dihe(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: icel1, icel2, icel3, icel4, i1, i2, i3, i4
    integer                  :: i, ix, icel_local

    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: buf_real(:,:,:), dihec(:,:,:)
    integer,         pointer :: ncel_local, nboundary
    integer(int2),   pointer :: cell_g2l(:), cell_g2b(:)
    integer(int2),   pointer :: cell_pair(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: dihed_add(:), dihed_exit(:)
    integer,         pointer :: dihed_exit_index(:,:)
    integer,         pointer :: ndihe(:), dihelist(:,:,:)
    integer,         pointer :: dihe_pbc(:,:,:)
    integer,         pointer :: buf_int(:,:,:)


    ncel_local       => domain%num_cell_local
    nboundary        => domain%num_cell_boundary
    cell_g2l         => domain%cell_g2l
    cell_g2b         => domain%cell_g2b
    cell_pair        => domain%cell_pair
    id_g2l           => domain%id_g2l
    coord            => domain%coord

    ndihe            => enefunc%num_rb_dihedral
    dihelist         => enefunc%rb_dihe_list
    dihec            => enefunc%rb_dihe_c
    dihed_add        => enefunc%rb_dihed_add
    dihed_exit       => enefunc%rb_dihed_exit
    dihed_exit_index => enefunc%rb_dihed_exit_index
    dihe_pbc         => enefunc%rb_dihe_pbc
    buf_int          => enefunc%buf_rb_dihed_integer
    buf_real         => enefunc%buf_rb_dihed_real

    dihed_add (1:ncel_local+nboundary) = 0
    dihed_exit(1:ncel_local) = 0


    do i = 1, ncel_local
      do ix = 1, ndihe(i)

        icel1 = id_g2l(1,dihelist(1,ix,i))
        i1    = id_g2l(2,dihelist(1,ix,i))
        icel2 = id_g2l(1,dihelist(2,ix,i))
        i2    = id_g2l(2,dihelist(2,ix,i))
        icel3 = id_g2l(1,dihelist(3,ix,i))
        i3    = id_g2l(2,dihelist(3,ix,i))
        icel4 = id_g2l(1,dihelist(4,ix,i))
        i4    = id_g2l(2,dihelist(4,ix,i))

        icel_local = cell_pair(icel1,icel4)

        if (icel_local /= i) then

          dihed_exit(i) = dihed_exit(i) + 1
          dihed_exit_index(dihed_exit(i),i) = ix
          dihed_add(icel_local) = dihed_add(icel_local) + 1
          buf_real(1,dihed_add(icel_local),icel_local) = dihec(1,ix,i)
          buf_real(2,dihed_add(icel_local),icel_local) = dihec(2,ix,i)
          buf_real(3,dihed_add(icel_local),icel_local) = dihec(3,ix,i)
          buf_real(4,dihed_add(icel_local),icel_local) = dihec(4,ix,i)
          buf_real(5,dihed_add(icel_local),icel_local) = dihec(5,ix,i)
          buf_real(6,dihed_add(icel_local),icel_local) = dihec(6,ix,i)
          buf_int (1,dihed_add(icel_local),icel_local) = dihelist(1,ix,i)
          buf_int (2,dihed_add(icel_local),icel_local) = dihelist(2,ix,i)
          buf_int (3,dihed_add(icel_local),icel_local) = dihelist(3,ix,i)
          buf_int (4,dihed_add(icel_local),icel_local) = dihelist(4,ix,i)
          buf_int (5,dihed_add(icel_local),icel_local) = dihe_pbc(1,ix,i)
          buf_int (6,dihed_add(icel_local),icel_local) = dihe_pbc(2,ix,i)
          buf_int (7,dihed_add(icel_local),icel_local) = dihe_pbc(3,ix,i)

        end if

      end  do
    end do

    return

  end subroutine update_outgoing_enefunc_rb_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_enefunc_rb_dihe
  !> @brief        update Ryckaert-Bellemans DIHEDRAL term for each cell in
  !!               potential energy function
  !! @authors      NT
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_incoming_enefunc_rb_dihe(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    integer                  :: i, j, k, ix, kx, ndihe
    logical                  :: insert

    real(wp),        pointer :: dihec(:,:,:), buf_real(:,:,:)
    integer,         pointer :: ncel_local, nboundary
    integer(int2),   pointer :: cell_g2l(:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: dihed_add(:), dihed_exit(:)
    integer,         pointer :: dihed_exit_index(:,:), buf_int(:,:,:)
    integer,         pointer :: ndihedral(:), dihelist(:,:,:)
    integer,         pointer :: dihe_pbc(:,:,:)



    ncel_local       => domain%num_cell_local
    nboundary        => domain%num_cell_boundary
    cell_g2l         => domain%cell_g2l
    id_g2l           => domain%id_g2l

    ndihedral        => enefunc%num_rb_dihedral
    dihelist         => enefunc%rb_dihe_list
    dihec            => enefunc%rb_dihe_c
    dihed_add        => enefunc%rb_dihed_add
    dihed_exit       => enefunc%rb_dihed_exit
    dihed_exit_index => enefunc%rb_dihed_exit_index
    dihe_pbc         => enefunc%rb_dihe_pbc
    buf_int          => enefunc%buf_rb_dihed_integer
    buf_real         => enefunc%buf_rb_dihed_real

    ndihe = enefunc%num_rb_dihe_all

    do i = 1, ncel_local

#ifdef DEBUG
      if (dihed_add(i)+ndihedral(i)-dihed_exit(i) > MaxDihe) &
        call error_msg( &
            'Debug: update_incoming_enefunc_rb_dihed> dihed is exceed MaxDihe')
#endif

      if (dihed_add(i) > DiheMove) &
        call error_msg( &
   'Debug: update_incoming_enefunc_rb_dihed> dihed migration is exceed maximum')

      if (dihed_add(i) >= dihed_exit(i)) then

        do k = 1, dihed_exit(i)
          dihelist (1,dihed_exit_index(k,i),i) = buf_int (1,k,i)
          dihelist (2,dihed_exit_index(k,i),i) = buf_int (2,k,i)
          dihelist (3,dihed_exit_index(k,i),i) = buf_int (3,k,i)
          dihelist (4,dihed_exit_index(k,i),i) = buf_int (4,k,i)
          dihe_pbc (1,dihed_exit_index(k,i),i) = buf_int (5,k,i)
          dihe_pbc (2,dihed_exit_index(k,i),i) = buf_int (6,k,i)
          dihe_pbc (3,dihed_exit_index(k,i),i) = buf_int (7,k,i)
          dihec    (1,dihed_exit_index(k,i),i) = buf_real(1,k,i)
          dihec    (2,dihed_exit_index(k,i),i) = buf_real(2,k,i)
          dihec    (3,dihed_exit_index(k,i),i) = buf_real(3,k,i)
          dihec    (4,dihed_exit_index(k,i),i) = buf_real(4,k,i)
          dihec    (5,dihed_exit_index(k,i),i) = buf_real(5,k,i)
          dihec    (6,dihed_exit_index(k,i),i) = buf_real(6,k,i)
        end do

        do k = dihed_exit(i)+1, dihed_add(i)
          ix = k + ndihedral(i) - dihed_exit(i)
          dihelist (1,ix,i) = buf_int (1,k,i)
          dihelist (2,ix,i) = buf_int (2,k,i)
          dihelist (3,ix,i) = buf_int (3,k,i)
          dihelist (4,ix,i) = buf_int (4,k,i)
          dihe_pbc (1,ix,i) = buf_int (5,k,i)
          dihe_pbc (2,ix,i) = buf_int (6,k,i)
          dihe_pbc (3,ix,i) = buf_int (7,k,i)
          dihec    (1,ix,i) = buf_real(1,k,i)
          dihec    (2,ix,i) = buf_real(2,k,i)
          dihec    (3,ix,i) = buf_real(3,k,i)
          dihec    (4,ix,i) = buf_real(4,k,i)
          dihec    (5,ix,i) = buf_real(5,k,i)
          dihec    (6,ix,i) = buf_real(6,k,i)
        end do

        ndihedral(i) = ndihedral(i) + dihed_add(i) - dihed_exit(i)

      else
        do k = 1, dihed_add(i)
          dihelist (1,dihed_exit_index(k,i),i) = buf_int (1,k,i)
          dihelist (2,dihed_exit_index(k,i),i) = buf_int (2,k,i)
          dihelist (3,dihed_exit_index(k,i),i) = buf_int (3,k,i)
          dihelist (4,dihed_exit_index(k,i),i) = buf_int (4,k,i)
          dihe_pbc (1,dihed_exit_index(k,i),i) = buf_int (5,k,i)
          dihe_pbc (2,dihed_exit_index(k,i),i) = buf_int (6,k,i)
          dihe_pbc (3,dihed_exit_index(k,i),i) = buf_int (7,k,i)
          dihec    (1,dihed_exit_index(k,i),i) = buf_real(1,k,i)
          dihec    (2,dihed_exit_index(k,i),i) = buf_real(2,k,i)
          dihec    (3,dihed_exit_index(k,i),i) = buf_real(3,k,i)
          dihec    (4,dihed_exit_index(k,i),i) = buf_real(4,k,i)
          dihec    (5,dihed_exit_index(k,i),i) = buf_real(5,k,i)
          dihec    (6,dihed_exit_index(k,i),i) = buf_real(6,k,i)
        end do

        j  = 0
        ix = ndihedral(i)
        k  = dihed_add(i) + 1

        do while (j < (dihed_exit(i)-dihed_add(i)))

          insert = .true.
          do kx = k, dihed_exit(i)
            if (ix == dihed_exit_index(kx,i)) then
              insert = .false.
              j = j + 1
              exit
            end if
          end do

          if (insert) then
            kx = dihed_exit_index(k,i)
            dihelist (1,kx,i) = dihelist (1,ix,i)
            dihelist (2,kx,i) = dihelist (2,ix,i)
            dihelist (3,kx,i) = dihelist (3,ix,i)
            dihelist (4,kx,i) = dihelist (4,ix,i)
            dihe_pbc (1,kx,i) = dihe_pbc (1,ix,i)
            dihe_pbc (2,kx,i) = dihe_pbc (2,ix,i)
            dihe_pbc (3,kx,i) = dihe_pbc (3,ix,i)
            dihec    (1,kx,i) = dihec    (1,ix,i)
            dihec    (2,kx,i) = dihec    (2,ix,i)
            dihec    (3,kx,i) = dihec    (3,ix,i)
            dihec    (4,kx,i) = dihec    (4,ix,i)
            dihec    (5,kx,i) = dihec    (5,ix,i)
            dihec    (6,kx,i) = dihec    (6,ix,i)
            j = j + 1
            k = k + 1
          end if

          ix = ix - 1

        end do
        ndihedral(i) = ndihedral(i) + dihed_add(i) - dihed_exit(i)

      end if

    end do

    return

  end subroutine update_incoming_enefunc_rb_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_enefunc_impr
  !> @brief        update IMPROPER DIHEDRAL term for each cell in potential
  !                energy function
  !! @authors      JJ, CK
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_enefunc_impr(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: icel1, icel2, icel3, icel4, i1, i2, i3, i4
    integer                  :: i, ix, icel_local

    real(wp),        pointer :: fc(:,:), phase(:,:)
    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: buf_real(:,:,:)
    integer,         pointer :: ncel_local, nboundary
    integer(int2),   pointer :: cell_g2l(:), cell_g2b(:)
    integer(int2),   pointer :: cell_pair(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: impr_add(:), impr_exit(:), impr_exit_index(:,:)
    integer,         pointer :: nperiod(:,:)
    integer,         pointer :: buf_int(:,:,:), nimproper(:), imprlist(:,:,:)
    integer,         pointer :: impr_pbc(:,:,:)



    ncel_local      => domain%num_cell_local
    nboundary       => domain%num_cell_boundary
    cell_g2l        => domain%cell_g2l
    cell_g2b        => domain%cell_g2b
    cell_pair       => domain%cell_pair
    id_g2l          => domain%id_g2l
    coord           => domain%coord

    nimproper       => enefunc%num_improper
    imprlist        => enefunc%impr_list
    impr_add        => enefunc%impr_add
    impr_exit       => enefunc%impr_exit
    impr_exit_index => enefunc%impr_exit_index
    fc              => enefunc%impr_force_const
    phase           => enefunc%impr_phase
    nperiod         => enefunc%impr_periodicity
    buf_int         => enefunc%buf_impr_integer
    buf_real        => enefunc%buf_impr_real
    impr_pbc        => enefunc%impr_pbc

    impr_exit(1:ncel_local) = 0
    impr_add (1:ncel_local+nboundary) = 0

    do i = 1, ncel_local
      do ix = 1, nimproper(i)

        icel1 = id_g2l(1,imprlist(1,ix,i))
        i1    = id_g2l(2,imprlist(1,ix,i))
        icel2 = id_g2l(1,imprlist(2,ix,i))
        i2    = id_g2l(2,imprlist(2,ix,i))
        icel3 = id_g2l(1,imprlist(3,ix,i))
        i3    = id_g2l(2,imprlist(3,ix,i))
        icel4 = id_g2l(1,imprlist(4,ix,i))
        i4    = id_g2l(2,imprlist(4,ix,i))

        icel_local = cell_pair(icel1,icel4)

        if (icel_local /= i) then

          impr_exit(i) = impr_exit(i) + 1
          impr_exit_index(impr_exit(i),i) = ix
          impr_add(icel_local) = impr_add(icel_local) + 1
          buf_real(1,impr_add(icel_local),icel_local) = phase(ix,i)
          buf_real(2,impr_add(icel_local),icel_local) = fc(ix,i)
          buf_int (1,impr_add(icel_local),icel_local) = imprlist(1,ix,i)
          buf_int (2,impr_add(icel_local),icel_local) = imprlist(2,ix,i)
          buf_int (3,impr_add(icel_local),icel_local) = imprlist(3,ix,i)
          buf_int (4,impr_add(icel_local),icel_local) = imprlist(4,ix,i)
          buf_int (5,impr_add(icel_local),icel_local) = nperiod(ix,i)
          buf_int (6,impr_add(icel_local),icel_local) = impr_pbc(1,ix,i)
          buf_int (7,impr_add(icel_local),icel_local) = impr_pbc(2,ix,i)
          buf_int (8,impr_add(icel_local),icel_local) = impr_pbc(3,ix,i)

        end if

      end  do
    end do

    return

  end subroutine update_outgoing_enefunc_impr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_enefunc_impr
  !> @brief        update IMPROPER DIHEDRAL term for each cell in potential
  !!               energy function
  !! @authors      JJ, CK
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_incoming_enefunc_impr(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    integer                  :: i, j, k, ix, kx
    logical                  :: insert

    real(wp),        pointer :: buf_real(:,:,:)
    real(wp),        pointer :: impr_force(:,:), impr_phase(:,:)
    integer,         pointer :: ncel_local, nboundary
    integer(2),      pointer :: cell_g2l(:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: impr_add(:), impr_exit(:), impr_exit_index(:,:)
    integer,         pointer :: impr_nperiod(:,:)
    integer,         pointer :: buf_int(:,:,:), nimproper(:), imprlist(:,:,:)
    integer,         pointer :: impr_pbc(:,:,:)


    ncel_local       => domain%num_cell_local
    nboundary        => domain%num_cell_boundary
    cell_g2l         => domain%cell_g2l
    id_g2l           => domain%id_g2l

    nimproper        => enefunc%num_improper
    imprlist         => enefunc%impr_list
    impr_add         => enefunc%impr_add
    impr_exit        => enefunc%impr_exit
    impr_exit_index  => enefunc%impr_exit_index
    impr_force       => enefunc%impr_force_const
    impr_phase       => enefunc%impr_phase
    impr_nperiod     => enefunc%impr_periodicity
    buf_int          => enefunc%buf_impr_integer
    buf_real         => enefunc%buf_impr_real
    impr_pbc         => enefunc%impr_pbc

    do i = 1, ncel_local
#ifdef DEBUG
      if (impr_add(i)+nimproper(i)-impr_exit(i) > MaxImpr) &
        call error_msg( &
            'Debug: update_incoming_enefunc_impr> improper is exceed MaxImpr')
#endif

      if (impr_add(i) > ImprMove) &
        call error_msg( &
   'Debug: update_incoming_enefunc_impr> improper migration is exceed maximum')

      if (impr_add(i) >= impr_exit(i)) then

        do k = 1, impr_exit(i)
          imprlist(1,impr_exit_index(k,i),i)   = buf_int (1,k,i)
          imprlist(2,impr_exit_index(k,i),i)   = buf_int (2,k,i)
          imprlist(3,impr_exit_index(k,i),i)   = buf_int (3,k,i)
          imprlist(4,impr_exit_index(k,i),i)   = buf_int (4,k,i)
          impr_nperiod(impr_exit_index(k,i),i) = buf_int (5,k,i)
          impr_pbc(1,impr_exit_index(k,i),i)   = buf_int (6,k,i)
          impr_pbc(2,impr_exit_index(k,i),i)   = buf_int (7,k,i)
          impr_pbc(3,impr_exit_index(k,i),i)   = buf_int (8,k,i)
          impr_phase(impr_exit_index(k,i),i)   = buf_real(1,k,i)
          impr_force(impr_exit_index(k,i),i)   = buf_real(2,k,i)
        end do

        do k = impr_exit(i)+1, impr_add(i)
          ix = k + nimproper(i) - impr_exit(i)
          imprlist(1,ix,i)   = buf_int (1,k,i)
          imprlist(2,ix,i)   = buf_int (2,k,i)
          imprlist(3,ix,i)   = buf_int (3,k,i)
          imprlist(4,ix,i)   = buf_int (4,k,i)
          impr_nperiod(ix,i) = buf_int (5,k,i)
          impr_pbc(1,ix,i)   = buf_int (6,k,i)
          impr_pbc(2,ix,i)   = buf_int (7,k,i)
          impr_pbc(3,ix,i)   = buf_int (8,k,i)
          impr_phase(ix,i)   = buf_real(1,k,i)
          impr_force(ix,i)   = buf_real(2,k,i)
        end do

        nimproper(i) = nimproper(i) + impr_add(i) - impr_exit(i)

      else
        do k = 1, impr_add(i)
          imprlist(1,impr_exit_index(k,i),i)   = buf_int (1,k,i)
          imprlist(2,impr_exit_index(k,i),i)   = buf_int (2,k,i)
          imprlist(3,impr_exit_index(k,i),i)   = buf_int (3,k,i)
          imprlist(4,impr_exit_index(k,i),i)   = buf_int (4,k,i)
          impr_nperiod(impr_exit_index(k,i),i) = buf_int (5,k,i)
          impr_pbc(1,impr_exit_index(k,i),i)   = buf_int (6,k,i)
          impr_pbc(2,impr_exit_index(k,i),i)   = buf_int (7,k,i)
          impr_pbc(3,impr_exit_index(k,i),i)   = buf_int (8,k,i)
          impr_phase(impr_exit_index(k,i),i)   = buf_real(1,k,i)
          impr_force(impr_exit_index(k,i),i)   = buf_real(2,k,i)
        end do

        j  = 0
        ix = nimproper(i)
        k  = impr_add(i) + 1

        do while (j < (impr_exit(i)-impr_add(i)))

          insert = .true.
          do kx = k, impr_exit(i)
            if (ix == impr_exit_index(kx,i)) then
              insert = .false.
              j = j + 1
              exit
            end if
          end do

          if (insert) then
            kx = impr_exit_index(k,i)
            imprlist(1,kx,i)   = imprlist(1,ix,i)
            imprlist(2,kx,i)   = imprlist(2,ix,i)
            imprlist(3,kx,i)   = imprlist(3,ix,i)
            imprlist(4,kx,i)   = imprlist(4,ix,i)
            impr_nperiod(kx,i) = impr_nperiod(ix,i)
            impr_pbc(1,kx,i)   = impr_pbc(1,ix,i)
            impr_pbc(2,kx,i)   = impr_pbc(2,ix,i)
            impr_pbc(3,kx,i)   = impr_pbc(3,ix,i)
            impr_phase(kx,i)   = impr_phase(ix,i)
            impr_force(kx,i)   = impr_force(ix,i)
            j = j + 1
            k = k + 1
          end if
          ix = ix - 1
        end do

        nimproper(i) = nimproper(i) + impr_add(i) - impr_exit(i)

      end if
    end do

    return

  end subroutine update_incoming_enefunc_impr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_enefunc_cmap
  !> @brief        update CMAP term for each cell in potential
  !                energy function
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_enefunc_cmap(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: icel1, icel2
    integer                  :: i, ix, k, icel_local

    real(wip),       pointer :: coord(:,:,:)
    integer,         pointer :: ncel_local, nboundary
    integer(int2),   pointer :: cell_g2l(:), cell_g2b(:)
    integer(int2),   pointer :: cell_pair(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: cmap_exit(:), cmap_add(:), cmap_exit_index(:,:)
    integer,         pointer :: ncmap(:), cmaplist(:,:,:), cmaptype(:,:)
    integer,         pointer :: buf_int(:,:,:)
    integer,         pointer :: cmap_pbc(:,:,:)


    ncel_local      => domain%num_cell_local
    nboundary       => domain%num_cell_boundary
    cell_g2l        => domain%cell_g2l
    cell_g2b        => domain%cell_g2b
    cell_pair       => domain%cell_pair
    id_g2l          => domain%id_g2l
    coord           => domain%coord

    ncmap           => enefunc%num_cmap
    cmaplist        => enefunc%cmap_list
    cmaptype        => enefunc%cmap_type
    cmap_add        => enefunc%cmap_add
    cmap_exit       => enefunc%cmap_exit
    cmap_pbc        => enefunc%cmap_pbc
    cmap_exit_index => enefunc%cmap_exit_index
    buf_int         => enefunc%buf_cmap_integer

    cmap_exit(1:ncel_local) = 0
    cmap_add (1:ncel_local+nboundary) = 0

    do i = 1, ncel_local
      do ix = 1, ncmap(i)

        icel1 = id_g2l(1,cmaplist(1,ix,i))
        icel2 = id_g2l(1,cmaplist(8,ix,i))
        icel_local = cell_pair(icel1,icel2)

        if (icel_local /= i) then

          cmap_exit(i) = cmap_exit(i) + 1
          cmap_exit_index(cmap_exit(i),i) = ix
          cmap_add(icel_local) = cmap_add(icel_local) + 1
          do k = 1, 8
            buf_int(k,cmap_add(icel_local),icel_local) =cmaplist(k,ix,i)
          end do
          buf_int(9,cmap_add(icel_local),icel_local) = cmaptype(ix,i)
          buf_int(10,cmap_add(icel_local),icel_local) = cmap_pbc(1,ix,i)
          buf_int(11,cmap_add(icel_local),icel_local) = cmap_pbc(2,ix,i)
          buf_int(12,cmap_add(icel_local),icel_local) = cmap_pbc(3,ix,i)
          buf_int(13,cmap_add(icel_local),icel_local) = cmap_pbc(4,ix,i)
          buf_int(14,cmap_add(icel_local),icel_local) = cmap_pbc(5,ix,i)
          buf_int(15,cmap_add(icel_local),icel_local) = cmap_pbc(6,ix,i)

        end if

      end  do
    end do

    return

  end subroutine update_outgoing_enefunc_cmap

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_enefunc_cmap
  !> @brief        update CMAP term for each cell in potential energy function
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_incoming_enefunc_cmap(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    integer                  :: i, j, k, ix, kx, ia
    logical                  :: insert

    integer,         pointer :: ncel_local, nboundary
    integer(int2),   pointer :: cell_g2l(:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: cmap_add(:), cmap_exit(:), cmap_exit_index(:,:)
    integer,         pointer :: buf_int(:,:,:), ncmap(:), cmap_list(:,:,:)
    integer,         pointer :: cmap_type(:,:)
    integer,         pointer :: cmap_pbc(:,:,:)



    ncel_local      => domain%num_cell_local
    nboundary       => domain%num_cell_boundary
    cell_g2l        => domain%cell_g2l
    id_g2l          => domain%id_g2l

    ncmap           => enefunc%num_cmap
    cmap_list       => enefunc%cmap_list
    cmap_type       => enefunc%cmap_type
    cmap_add        => enefunc%cmap_add
    cmap_exit       => enefunc%cmap_exit
    cmap_exit_index => enefunc%cmap_exit_index
    cmap_pbc        => enefunc%cmap_pbc
    buf_int         => enefunc%buf_cmap_integer


    do i = 1, ncel_local

#ifdef DEBUG
      if (cmap_add(i)+ncmap(i)-cmap_exit(i) > MaxCmap) &
        call error_msg( &
          'Debug: update_incoming_enefunc_cmap> cmap is exceed MaxCmap')
#endif

      if (cmap_add(i) > CmapMove) &
        call error_msg( &
      'Debug: update_incoming_enefunc_cmap> cmap migration is exceed maximum')

      if (cmap_add(i) >= cmap_exit(i)) then

        do k = 1, cmap_exit(i)
          do ia = 1, 8
            cmap_list(ia,cmap_exit_index(k,i),i) = buf_int(ia,k,i)
          end do
          cmap_type(cmap_exit_index(k,i),i) = buf_int(9,k,i)
          cmap_pbc (1,cmap_exit_index(k,i),i) = buf_int(10,k,i)
          cmap_pbc (2,cmap_exit_index(k,i),i) = buf_int(11,k,i)
          cmap_pbc (3,cmap_exit_index(k,i),i) = buf_int(12,k,i)
          cmap_pbc (4,cmap_exit_index(k,i),i) = buf_int(13,k,i)
          cmap_pbc (5,cmap_exit_index(k,i),i) = buf_int(14,k,i)
          cmap_pbc (6,cmap_exit_index(k,i),i) = buf_int(15,k,i)
        end do

        do k = cmap_exit(i)+1, cmap_add(i)

          ix = k + ncmap(i) - cmap_exit(i)
          do ia = 1, 8
            cmap_list(ia,ix,i) = buf_int(ia,k,i)
          end do
          cmap_type(ix,i) = buf_int(9,k,i)
          cmap_pbc (1,ix,i) = buf_int(10,k,i)
          cmap_pbc (2,ix,i) = buf_int(11,k,i)
          cmap_pbc (3,ix,i) = buf_int(12,k,i)
          cmap_pbc (4,ix,i) = buf_int(13,k,i)
          cmap_pbc (5,ix,i) = buf_int(14,k,i)
          cmap_pbc (6,ix,i) = buf_int(15,k,i)

        end do

        ncmap(i) = ncmap(i) + cmap_add(i) - cmap_exit(i)

      else

        do k = 1, cmap_add(i)
          do ia = 1, 8
            cmap_list(ia,cmap_exit_index(k,i),i) = buf_int(ia,k,i)
          end do
          cmap_type(cmap_exit_index(k,i),i)   = buf_int(9,k,i)
          cmap_pbc (1,cmap_exit_index(k,i),i) = buf_int(10,k,i)
          cmap_pbc (2,cmap_exit_index(k,i),i) = buf_int(11,k,i)
          cmap_pbc (3,cmap_exit_index(k,i),i) = buf_int(12,k,i)
          cmap_pbc (4,cmap_exit_index(k,i),i) = buf_int(13,k,i)
          cmap_pbc (5,cmap_exit_index(k,i),i) = buf_int(14,k,i)
          cmap_pbc (6,cmap_exit_index(k,i),i) = buf_int(15,k,i)
        end do

        j  = 0
        ix = ncmap(i)
        k  = cmap_add(i) + 1

        do while (j < (cmap_exit(i)-cmap_add(i)))

          insert = .true.
          do kx = k, cmap_exit(i)
            if (ix == cmap_exit_index(kx,i)) then
              insert = .false.
              j = j + 1
              exit
            end if
          end do

          if (insert) then

            kx = cmap_exit_index(k,i)
            do ia = 1, 8
              cmap_list(ia,kx,i) = cmap_list(ia,ix,i)
            end do

            cmap_type(kx,i) = cmap_type(ix,i)
            cmap_pbc (1,kx,i) = cmap_pbc(1,ix,i)
            cmap_pbc (2,kx,i) = cmap_pbc(2,ix,i)
            cmap_pbc (3,kx,i) = cmap_pbc(3,ix,i)
            cmap_pbc (4,kx,i) = cmap_pbc(4,ix,i)
            cmap_pbc (5,kx,i) = cmap_pbc(5,ix,i)
            cmap_pbc (6,kx,i) = cmap_pbc(6,ix,i)
            j = j + 1
            k = k + 1

          end if

          ix = ix - 1

        end do

        ncmap(i) = ncmap(i) + cmap_add(i) - cmap_exit(i)

      end if

    end do

    return

  end subroutine update_incoming_enefunc_cmap

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_enefunc_restraint
  !> @brief        update restraint term for each cell
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_enefunc_restraint(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    integer                  :: i, ix, icel_local

    real(wp),        pointer :: coord(:,:,:)
    real(wp),        pointer :: force(:,:,:), buf_real(:,:,:)
    integer,         pointer :: ncel_local, nboundary
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: nrestraint(:), restraint_atom(:,:)
    integer,         pointer :: restraint_exit(:), index(:,:)
    integer,         pointer :: restraint_add(:), buf_int(:,:)


    if ((.not. enefunc%restraint_posi) .and. (.not. enefunc%restraint_rmsd)) &
      return

    ncel_local     => domain%num_cell_local
    nboundary      => domain%num_cell_boundary
    id_g2l         => domain%id_g2l

    nrestraint     => enefunc%num_restraint
    restraint_atom => enefunc%restraint_atom
    force          => enefunc%restraint_force
    coord          => enefunc%restraint_coord
    restraint_exit => enefunc%restraint_exit
    index          => enefunc%restraint_exit_index
    restraint_add  => enefunc%restraint_add
    buf_int        => enefunc%buf_restraint_integer
    buf_real       => enefunc%buf_restraint_real

    restraint_exit(1:ncel_local) = 0
    restraint_add (1:ncel_local+nboundary) = 0

    do i = 1, ncel_local
      do ix = 1, nrestraint(i)

        icel_local = id_g2l(1,restraint_atom(ix,i))

        if (icel_local /= i) then

          restraint_exit(i) = restraint_exit(i) + 1
          index(restraint_exit(i),i) = ix

          restraint_add(icel_local) = restraint_add(icel_local) + 1
          buf_real(1:4,restraint_add(icel_local),icel_local) = force(1:4,ix,i)
          buf_real(5:7,restraint_add(icel_local),icel_local) = coord(1:3,ix,i)
          buf_int (restraint_add(icel_local),icel_local) = restraint_atom(ix,i)

        end if

      end do
    end do

    return

  end subroutine update_outgoing_enefunc_restraint

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_enefunc_restraint
  !> @brief        update incoing restraint term for each cell
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_incoming_enefunc_restraint(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    integer                  :: i, j, k, ix, kx, found
    logical                  :: insert

    real(wp),        pointer :: buf_real(:,:,:)
    real(wp),        pointer :: force(:,:,:), coord(:,:,:)
    integer,         pointer :: ncel_local, nboundary
    integer(int2),   pointer :: cell_g2l(:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: restraint_add(:), restraint_exit(:)
    integer,         pointer :: restraint_exit_index(:,:)
    integer,         pointer :: buf_int(:,:), nrestraint(:)
    integer,         pointer :: restraint_atom(:,:)


    if ((.not. enefunc%restraint_posi) .and. (.not. enefunc%restraint_rmsd)) &
      return

    ncel_local           => domain%num_cell_local
    nboundary            => domain%num_cell_boundary
    cell_g2l             => domain%cell_g2l
    id_g2l               => domain%id_g2l

    restraint_add        => enefunc%restraint_add
    restraint_exit       => enefunc%restraint_exit
    restraint_exit_index => enefunc%restraint_exit_index
    buf_int              => enefunc%buf_restraint_integer
    buf_real             => enefunc%buf_restraint_real
    restraint_atom       => enefunc%restraint_atom
    force                => enefunc%restraint_force
    coord                => enefunc%restraint_coord
    nrestraint           => enefunc%num_restraint

    found = 0
    do i = 1, ncel_local

#ifdef DEBUG
      if (restraint_add(i)+nrestraint(i)-restraint_exit(i) > MaxRest) &
        call error_msg( &
          'Debug: update_incoming_enefunc_rest> restraint is exceed MaxRest')
#endif

      if (restraint_add(i) > RestMove) &
        call error_msg( &
  'Debug: update_incoming_enefunc_rest> restraint migration is exceed maximum')

      if (restraint_add(i) >= restraint_exit(i)) then

        do k = 1, restraint_exit(i)
          restraint_atom(restraint_exit_index(k,i),i) = buf_int(k,i)
          force(1:4,restraint_exit_index(k,i),i) = buf_real(1:4,k,i)
          coord(1:3,restraint_exit_index(k,i),i) = buf_real(5:7,k,i)
        end do

        do k = restraint_exit(i)+1, restraint_add(i)
          ix = k + nrestraint(i) - restraint_exit(i)
          restraint_atom(ix,i) = buf_int(k,i)
          force(1:4,ix,i) = buf_real(1:4,k,i)
          coord(1:3,ix,i) = buf_real(5:7,k,i)
        end do

        nrestraint(i) = nrestraint(i) + restraint_add(i) - restraint_exit(i)

      else

        do k = 1, restraint_add(i)
          restraint_atom(restraint_exit_index(k,i),i) = buf_int(k,i)
          force(1:4,restraint_exit_index(k,i),i)  = buf_real(1:4,k,i)
          coord(1:3,restraint_exit_index(k,i),i)  = buf_real(5:7,k,i)
        end do

        j  = 0
        ix = nrestraint(i)
        k  = restraint_add(i) + 1

        do while (j < (restraint_exit(i)-restraint_add(i)))

          insert = .true.
          do kx = k, restraint_exit(i)
            if (ix == restraint_exit_index(kx,i)) then
              insert = .false.
              j = j + 1
              exit
            end if
          end do

          if (insert) then

            kx = restraint_exit_index(k,i)
            restraint_atom(kx,i) = restraint_atom(ix,i)
            force(1:4,kx,i) = force(1:4,ix,i)
            coord(1:3,kx,i) = coord(1:3,ix,i)

            j = j + 1
            k = k + 1

          end if

          ix = ix - 1

        end do

        nrestraint(i) = nrestraint(i) + restraint_add(i) - restraint_exit(i)

      end if
    end do

    return

  end subroutine update_incoming_enefunc_restraint

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_enefunc_fitting
  !> @brief        update fitting term for each cell
  !! @authors      KT
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_enefunc_fitting(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    integer                  :: i, ix, icel_local

    real(wp),        pointer :: fit_coord(:,:,:)
    real(wp),        pointer :: force(:,:,:), buf_real(:,:,:)
    integer,         pointer :: ncel_local, nboundary
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: nfitting(:), fitting_atom(:,:)
    integer,         pointer :: fitting_exit(:), index(:,:)
    integer,         pointer :: fitting_add(:), buf_int(:,:)


    if (.not. enefunc%do_fitting) return

    ncel_local     => domain%num_cell_local
    nboundary      => domain%num_cell_boundary
    id_g2l         => domain%id_g2l

    nfitting       => enefunc%nfitting
    fitting_atom   => enefunc%fitting_atom
    fit_coord      => enefunc%fit_coord
    fitting_exit   => enefunc%fitting_exit
    index          => enefunc%fitting_exit_index
    fitting_add    => enefunc%fitting_add
    buf_int        => enefunc%buf_fitting_integer
    buf_real       => enefunc%buf_fitting_real

    fitting_exit(1:ncel_local) = 0
    fitting_add (1:ncel_local+nboundary) = 0

    do i = 1, ncel_local
      do ix = 1, nfitting(i)

        icel_local = id_g2l(1,fitting_atom(ix,i))

        if (icel_local /= i) then

          fitting_exit(i) = fitting_exit(i) + 1
          index(fitting_exit(i),i) = ix

          fitting_add(icel_local) = fitting_add(icel_local) + 1
          buf_real(1:3,fitting_add(icel_local),icel_local) = fit_coord(1:3,ix,i)
          buf_int (fitting_add(icel_local),icel_local) = fitting_atom(ix,i)

        end if

      end do
    end do

    return

  end subroutine update_outgoing_enefunc_fitting

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_enefunc_fitting
  !> @brief        update incoming fitting term for each cell
  !! @authors      KT
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_incoming_enefunc_fitting(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    integer                  :: i, j, k, ix, kx, found
    logical                  :: insert

    real(wp),        pointer :: buf_real(:,:,:)
    real(wp),        pointer :: fit_coord(:,:,:)
    integer,         pointer :: ncel_local, nboundary
    integer(int2),   pointer :: cell_g2l(:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: fitting_add(:), fitting_exit(:)
    integer,         pointer :: fitting_exit_index(:,:)
    integer,         pointer :: buf_int(:,:), nfitting(:)
    integer,         pointer :: fitting_atom(:,:)


    if (.not. enefunc%do_fitting) return

    ncel_local           => domain%num_cell_local
    nboundary            => domain%num_cell_boundary
    cell_g2l             => domain%cell_g2l
    id_g2l               => domain%id_g2l

    fitting_add          => enefunc%fitting_add
    fitting_exit         => enefunc%fitting_exit
    fitting_exit_index   => enefunc%fitting_exit_index
    buf_int              => enefunc%buf_fitting_integer
    buf_real             => enefunc%buf_fitting_real
    fitting_atom         => enefunc%fitting_atom
    fit_coord            => enefunc%fit_coord
    nfitting             => enefunc%nfitting

    found = 0
    do i = 1, ncel_local

#ifdef DEBUG
      if (fitting_add(i)+nfitting(i)-fitting_exit(i) > MaxRest) &
        call error_msg( &
          'Debug: update_incoming_enefunc_fit> fitting atoms exceed MaxRest')
#endif

      if (fitting_add(i) > RestMove) &
        call error_msg( &
  'Debug: update_incoming_enefunc_fit> fitting migration exceeds maximum')

      if (fitting_add(i) >= fitting_exit(i)) then

        do k = 1, fitting_exit(i)
          fitting_atom(fitting_exit_index(k,i),i)  = buf_int(k,i)
          fit_coord(1:3,fitting_exit_index(k,i),i) = buf_real(1:3,k,i)
        end do

        do k = fitting_exit(i)+1, fitting_add(i)
          ix = k + nfitting(i) - fitting_exit(i)
          fitting_atom(ix,i)  = buf_int(k,i)
          fit_coord(1:3,ix,i) = buf_real(1:3,k,i)
        end do

        nfitting(i) = nfitting(i) + fitting_add(i) - fitting_exit(i)

      else

        do k = 1, fitting_add(i)
          fitting_atom(fitting_exit_index(k,i),i)  = buf_int(k,i)
          fit_coord(1:3,fitting_exit_index(k,i),i) = buf_real(1:3,k,i)
        end do

        j  = 0
        ix = nfitting(i)
        k  = fitting_add(i) + 1

        do while (j < (fitting_exit(i)-fitting_add(i)))

          insert = .true.
          do kx = k, fitting_exit(i)
            if (ix == fitting_exit_index(kx,i)) then
              insert = .false.
              j = j + 1
              exit
            end if
          end do

          if (insert) then

            kx = fitting_exit_index(k,i)
            fitting_atom(kx,i)  = fitting_atom(ix,i)
            fit_coord(1:3,kx,i) = fit_coord(1:3,ix,i)

            j = j + 1
            k = k + 1

          end if

          ix = ix - 1

        end do

        nfitting(i) = nfitting(i) + fitting_add(i) - fitting_exit(i)

      end if
    end do

    return

  end subroutine update_incoming_enefunc_fitting

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_atom_fep
  !> @brief        check particles (not bonded to hydrogen) going other cells
  !                for FEP
  !! @authors      HO
  !! @param[in]    boundary    : boundary condition information
  !! @param[in]    constraints : constraints information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_atom_fep(boundary, constraints, domain)

    ! formal arguments
    type(s_boundary),    target, intent(in)    :: boundary
    type(s_constraints), target, intent(in)    :: constraints
    type(s_domain),      target, intent(inout) :: domain

    ! local variable
    real(wip)                    :: x_shift, y_shift, z_shift
    real(wip)                    :: move(3)
    integer                      :: i, k, ix, icx, icy, icz, icel, ncel
    integer                      :: icel_local, icel_bd

    real(wip),           pointer :: bsize_x, bsize_y, bsize_z
    real(wip),           pointer :: csize_x, csize_y, csize_z
    real(wip),           pointer :: coord(:,:,:), velocity(:,:,:)
    real(wp),            pointer :: charge(:,:)
    real(wip),           pointer :: mass(:,:)
    real(wip),           pointer :: buf_real(:,:,:)
    integer,             pointer :: ncel_x, ncel_y, ncel_z, ncel_local, ncel_bd
    integer,             pointer :: nsolute(:)
    integer(int2),       pointer :: cell_g2l(:), cell_g2b(:)
    integer,             pointer :: atmcls(:,:), id_l2g(:,:), id_l2g_sol(:,:)
    integer,             pointer :: ptl_add(:), ptl_exit(:), ptlindex(:,:)
    integer,             pointer :: buf_int(:,:,:)

    ! FEP
    integer,             pointer :: fepgrp(:,:)
    integer,             pointer :: fep_atmcls_singleB(:,:)
    real(wp),            pointer :: fep_chargeA(:,:)
    real(wp),            pointer :: fep_chargeB(:,:)

    bsize_x    => boundary%box_size_x
    bsize_y    => boundary%box_size_y
    bsize_z    => boundary%box_size_z
    ncel_x     => boundary%num_cells_x
    ncel_y     => boundary%num_cells_y
    ncel_z     => boundary%num_cells_z
    csize_x    => boundary%cell_size_x
    csize_y    => boundary%cell_size_y
    csize_z    => boundary%cell_size_z

    nsolute    => constraints%No_HGr

    ncel_local => domain%num_cell_local
    ncel_bd    => domain%num_cell_boundary
    cell_g2l   => domain%cell_g2l
    cell_g2b   => domain%cell_g2b
    coord      => domain%coord
    velocity   => domain%velocity
    charge     => domain%charge
    mass       => domain%mass
    atmcls     => domain%atom_cls_no
    id_l2g     => domain%id_l2g
    id_l2g_sol => domain%id_l2g_solute
    ptl_add    => domain%ptl_add
    ptl_exit   => domain%ptl_exit
    ptlindex   => domain%ptl_exit_index
    buf_int    => domain%buf_integer
    buf_real   => domain%buf_real

    ! FEP
    fepgrp     => domain%fepgrp
    fep_chargeA => domain%fep_chargeA
    fep_chargeB => domain%fep_chargeB
    fep_atmcls_singleB => domain%fep_atmcls_singleB

    ! initializaiton
    !
    ncel = ncel_local + ncel_bd
    ptl_exit(1:ncel_local) = 0
    ptl_add(1:ncel) = 0

    ! Check outgoing particles
    !
    do i = 1, ncel_local

      k = 0
      do ix = 1, nsolute(i)

        x_shift = coord(1,ix,i) - boundary%origin_x
        y_shift = coord(2,ix,i) - boundary%origin_y
        z_shift = coord(3,ix,i) - boundary%origin_z

        !coordinate shifted to the first quadrant and set into the boundary box
        move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
        move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
        move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)

        x_shift = x_shift + move(1)
        y_shift = y_shift + move(2)
        z_shift = z_shift + move(3)

        !assign which cell
        icx = int(x_shift/csize_x)
        icy = int(y_shift/csize_y)
        icz = int(z_shift/csize_z)
        if (icx == ncel_x) icx = icx - 1
        if (icy == ncel_y) icy = icy - 1
        if (icz == ncel_z) icz = icz - 1
        icel = 1 + icx + icy*ncel_x + icz*ncel_x*ncel_y
        icel_local = cell_g2l(icel)
        icel_bd    = cell_g2b(icel)

        if (icel_local /= i) then

          ptl_exit(i) = ptl_exit(i) + 1
          ptlindex(ptl_exit(i),i) = ix

          if (icel_local /= 0) then

            ptl_add (icel_local) = ptl_add(icel_local) + 1
            buf_real( 1,ptl_add(icel_local),icel_local) = coord(1,ix,i)
            buf_real( 2,ptl_add(icel_local),icel_local) = coord(2,ix,i)
            buf_real( 3,ptl_add(icel_local),icel_local) = coord(3,ix,i)
            buf_real( 4,ptl_add(icel_local),icel_local) = velocity(1,ix,i)
            buf_real( 5,ptl_add(icel_local),icel_local) = velocity(2,ix,i)
            buf_real( 6,ptl_add(icel_local),icel_local) = velocity(3,ix,i)
            buf_real( 7,ptl_add(icel_local),icel_local) = charge(ix,i)
            buf_real( 8,ptl_add(icel_local),icel_local) = mass(ix,i)
            buf_real( 9,ptl_add(icel_local),icel_local) = fep_chargeA(ix,i)
            buf_real(10,ptl_add(icel_local),icel_local) = fep_chargeB(ix,i)
            buf_int ( 1,ptl_add(icel_local),icel_local) = atmcls(ix,i)
            buf_int ( 2,ptl_add(icel_local),icel_local) = id_l2g(ix,i)
            buf_int ( 3,ptl_add(icel_local),icel_local) = id_l2g_sol(ix,i)
            buf_int ( 4,ptl_add(icel_local),icel_local) = fepgrp(ix,i)
            buf_int ( 5,ptl_add(icel_local),icel_local) = fep_atmcls_singleB(ix,i)

          else if (icel_bd /= 0) then

            icel_bd = icel_bd + ncel_local
            ptl_add (icel_bd) = ptl_add(icel_bd) + 1
            buf_real( 1,ptl_add(icel_bd),icel_bd) = coord(1,ix,i)
            buf_real( 2,ptl_add(icel_bd),icel_bd) = coord(2,ix,i)
            buf_real( 3,ptl_add(icel_bd),icel_bd) = coord(3,ix,i)
            buf_real( 4,ptl_add(icel_bd),icel_bd) = velocity(1,ix,i)
            buf_real( 5,ptl_add(icel_bd),icel_bd) = velocity(2,ix,i)
            buf_real( 6,ptl_add(icel_bd),icel_bd) = velocity(3,ix,i)
            buf_real( 7,ptl_add(icel_bd),icel_bd) = charge(ix,i)
            buf_real( 8,ptl_add(icel_bd),icel_bd) = mass(ix,i)
            buf_real( 9,ptl_add(icel_bd),icel_bd) = fep_chargeA(ix,i)
            buf_real(10,ptl_add(icel_bd),icel_bd) = fep_chargeB(ix,i)
            buf_int ( 1,ptl_add(icel_bd),icel_bd) = atmcls(ix,i)
            buf_int ( 2,ptl_add(icel_bd),icel_bd) = id_l2g(ix,i)
            buf_int ( 3,ptl_add(icel_bd),icel_bd) = id_l2g_sol(ix,i)
            buf_int ( 4,ptl_add(icel_bd),icel_bd) = fepgrp(ix,i)
            buf_int ( 5,ptl_add(icel_bd),icel_bd) = fep_atmcls_singleB(ix,i)

          end if

        end if

       end do
    end do

    return

  end subroutine update_outgoing_atom_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_water_fep
  !> @brief        check water particles going/remaining cells for FEP
  !! @authors      HO
  !! @param[in]    water_atom : number of water atoms
  !! @param[in]    boundary   : boundary condition information
  !! @param[inout] domain     : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_water_fep(water_atom, boundary, domain)

    ! formal arguments
    integer,                  intent(in)    :: water_atom
    type(s_boundary), target, intent(in)    :: boundary
    type(s_domain),   target, intent(inout) :: domain

    ! local variable
    real(wip)                 :: x_shift, y_shift, z_shift
    real(wip)                 :: move(3)
    integer                   :: i, k, iwater, list, ix
    integer                   :: icx, icy, icz, icel, ncel
    integer                   :: icel_local, icel_bd

    real(wip),        pointer :: bsize_x, bsize_y, bsize_z
    real(wip),        pointer :: csize_x, csize_y, csize_z
    real(wip),        pointer :: coord(:,:,:), velocity(:,:,:)
    real(wip),        pointer :: water_move_real(:,:,:), water_stay_real(:,:,:)
    integer,          pointer :: ncel_x, ncel_y, ncel_z, ncel_local, ncel_bd
    integer,          pointer :: nwater(:), water_list(:,:,:)
    integer,          pointer :: id_l2g(:,:)
    integer(int2),    pointer :: cell_g2l(:), cell_g2b(:)
    integer,          pointer :: water_move(:), water_stay(:)
    integer,          pointer :: water_move_int(:,:,:), water_stay_int(:,:,:)

    ! FEP
    integer,          pointer :: fepgrp(:,:)

    bsize_x         => boundary%box_size_x
    bsize_y         => boundary%box_size_y
    bsize_z         => boundary%box_size_z
    ncel_x          => boundary%num_cells_x
    ncel_y          => boundary%num_cells_y
    ncel_z          => boundary%num_cells_z
    csize_x         => boundary%cell_size_x
    csize_y         => boundary%cell_size_y
    csize_z         => boundary%cell_size_z

    ncel_local      => domain%num_cell_local
    ncel_bd         => domain%num_cell_boundary
    nwater          => domain%num_water
    water_list      => domain%water_list
    cell_g2l        => domain%cell_g2l
    cell_g2b        => domain%cell_g2b
    coord           => domain%coord
    velocity        => domain%velocity
    id_l2g          => domain%id_l2g
    water_move      => domain%water%move
    water_stay      => domain%water%stay
    water_move_real => domain%water%move_real
    water_stay_real => domain%water%stay_real
    water_move_int  => domain%water%move_integer
    water_stay_int  => domain%water%stay_integer

    ! FEP
    fepgrp          => domain%fepgrp


    ! initializaiton
    !
    ncel = ncel_local + ncel_bd
    water_stay(1:ncel_local) = 0
    water_move(1:ncel) = 0

    ! Check outgoing particles
    !
    do i = 1, ncel_local

      k = 0
      do iwater = 1, nwater(i)

        ix = water_list(1,iwater,i)
        x_shift = coord(1,ix,i) - boundary%origin_x
        y_shift = coord(2,ix,i) - boundary%origin_y
        z_shift = coord(3,ix,i) - boundary%origin_z

        !coordinate shifted to the first quadrant and set into the boundary box
        move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
        move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
        move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)

        x_shift = x_shift + move(1)
        y_shift = y_shift + move(2)
        z_shift = z_shift + move(3)

        !assign which cell
        icx = int(x_shift/csize_x)
        icy = int(y_shift/csize_y)
        icz = int(z_shift/csize_z)
        if (icx == ncel_x) icx = icx - 1
        if (icy == ncel_y) icy = icy - 1
        if (icz == ncel_z) icz = icz - 1
        icel = 1 + icx + icy*ncel_x + icz*ncel_x*ncel_y
        icel_local = cell_g2l(icel)
        icel_bd    = cell_g2b(icel)

        ! atom index
        if (icel_local /= i) then

          if (icel_local /= 0) then

            water_move(icel_local) = water_move(icel_local) + 1

            do list = 1, water_atom
              ix = water_list(list,iwater,i)
              water_move_real(6*(list-1)+1,water_move(icel_local),icel_local) &
                 = coord(1,ix,i)
              water_move_real(6*(list-1)+2,water_move(icel_local),icel_local) &
                 = coord(2,ix,i)
              water_move_real(6*(list-1)+3,water_move(icel_local),icel_local) &
                 = coord(3,ix,i)
              water_move_real(6*(list-1)+4,water_move(icel_local),icel_local) &
                 = velocity(1,ix,i)
              water_move_real(6*(list-1)+5,water_move(icel_local),icel_local) &
                 = velocity(2,ix,i)
              water_move_real(6*(list-1)+6,water_move(icel_local),icel_local) &
                 = velocity(3,ix,i)
              water_move_int(2*(list-1)+1,water_move(icel_local),icel_local) &
                 = id_l2g(ix,i)
              water_move_int(2*(list-1)+2,water_move(icel_local),icel_local) &
                 = fepgrp(ix,i)
            end do

          else if (icel_bd /= 0) then

            icel_bd = icel_bd + ncel_local
            water_move(icel_bd) = water_move(icel_bd) + 1

            do list = 1, water_atom
              ix = water_list(list,iwater,i)
              water_move_real(6*(list-1)+1,water_move(icel_bd),icel_bd)       &
                 = coord(1,ix,i)
              water_move_real(6*(list-1)+2,water_move(icel_bd),icel_bd)       &
                 = coord(2,ix,i)
              water_move_real(6*(list-1)+3,water_move(icel_bd),icel_bd)       &
                 = coord(3,ix,i)
              water_move_real(6*(list-1)+4,water_move(icel_bd),icel_bd)       &
                 = velocity(1,ix,i)
              water_move_real(6*(list-1)+5,water_move(icel_bd),icel_bd)       &
                 = velocity(2,ix,i)
              water_move_real(6*(list-1)+6,water_move(icel_bd),icel_bd)       &
                 = velocity(3,ix,i)
              water_move_int(2*(list-1)+1,water_move(icel_bd),icel_bd) &
                 = id_l2g(ix,i)
              water_move_int(2*(list-1)+2,water_move(icel_bd),icel_bd) &
                 = fepgrp(ix,i)
            end do
          end if

        else

          water_stay(i) = water_stay(i) + 1

          do list = 1, water_atom
            ix = water_list(list,iwater,i)
            water_stay_real(6*(list-1)+1,water_stay(i),i) = coord(1,ix,i)
            water_stay_real(6*(list-1)+2,water_stay(i),i) = coord(2,ix,i)
            water_stay_real(6*(list-1)+3,water_stay(i),i) = coord(3,ix,i)
            water_stay_real(6*(list-1)+4,water_stay(i),i) = velocity(1,ix,i)
            water_stay_real(6*(list-1)+5,water_stay(i),i) = velocity(2,ix,i)
            water_stay_real(6*(list-1)+6,water_stay(i),i) = velocity(3,ix,i)
            water_stay_int(2*(list-1)+1,water_stay(i),i) = id_l2g(ix,i)
            water_stay_int(2*(list-1)+2,water_stay(i),i) = fepgrp(ix,i)
          end do

        end if

      end do

    end do

    return

  end subroutine update_outgoing_water_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_HGr_fep
  !> @brief        check hydrogen-bonded particles going/remaining cells
  !                for FEP
  !! @authors      HO
  !! @param[in]    boundary    : boundary condition information
  !! @param[inout] constraints : constraints information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_HGr_fep(boundary, constraints, domain)

    ! formal arguments
    type(s_boundary),    target, intent(in)    :: boundary
    type(s_constraints), target, intent(inout) :: constraints
    type(s_domain),      target, intent(inout) :: domain

    ! local variable
    real(wip)                    :: x_shift, y_shift, z_shift
    real(wip)                    :: move(3)
    integer                      :: i, j, k, list, ix
    integer                      :: num_move, num_stay
    integer                      :: icx, icy, icz, icel, ncel
    integer                      :: icel_local, icel_bd

    real(wip),           pointer :: bsize_x, bsize_y, bsize_z
    real(wip),           pointer :: csize_x, csize_y, csize_z
    real(wip),           pointer :: coord(:,:,:), vel(:,:,:)
    real(wp),            pointer :: charge(:,:)
    real(wip),           pointer :: mass(:,:)
    real(wip),           pointer :: HGr_bond_dist(:,:,:,:)
    real(wip),           pointer :: HGr_move_real(:,:,:,:)
    real(wip),           pointer :: HGr_stay_real(:,:,:,:)
    integer,             pointer :: ncel_x, ncel_y, ncel_z
    integer,             pointer :: ncel_local, ncel_bd
    integer,             pointer :: connect
    integer,             pointer :: atmcls(:,:)
    integer,             pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)
    integer,             pointer :: id_l2g(:,:), id_l2g_sol(:,:)
    integer(int2),       pointer :: cell_g2l(:), cell_g2b(:)
    integer,             pointer :: HGr_move(:,:), HGr_stay(:,:)
    integer,             pointer :: HGr_move_int(:,:,:,:), HGr_stay_int(:,:,:,:)

    ! FEP
    integer,             pointer :: fepgrp(:,:)
    integer,             pointer :: fep_atmcls_singleB(:,:)
    real(wp),            pointer :: fep_chargeA(:,:)
    real(wp),            pointer :: fep_chargeB(:,:)

    bsize_x       => boundary%box_size_x
    bsize_y       => boundary%box_size_y
    bsize_z       => boundary%box_size_z
    ncel_x        => boundary%num_cells_x
    ncel_y        => boundary%num_cells_y
    ncel_z        => boundary%num_cells_z
    csize_x       => boundary%cell_size_x
    csize_y       => boundary%cell_size_y
    csize_z       => boundary%cell_size_z

    HGr_local     => constraints%HGr_local
    HGr_bond_list => constraints%HGr_bond_list
    HGr_move      => constraints%HGr_move
    HGr_stay      => constraints%HGr_stay
    HGr_bond_dist => constraints%HGr_bond_dist
    HGr_move_real => constraints%HGr_move_real
    HGr_stay_real => constraints%HGr_stay_real
    HGr_move_int  => constraints%HGr_move_int
    HGr_stay_int  => constraints%HGr_stay_int
    connect       => constraints%connect

    ncel_local    => domain%num_cell_local
    ncel_bd       => domain%num_cell_boundary
    cell_g2l      => domain%cell_g2l
    cell_g2b      => domain%cell_g2b
    coord         => domain%coord
    vel           => domain%velocity
    charge        => domain%charge
    mass          => domain%mass
    id_l2g        => domain%id_l2g
    id_l2g_sol    => domain%id_l2g_solute
    atmcls        => domain%atom_cls_no

    ! FEP
    fepgrp        => domain%fepgrp
    fep_chargeA => domain%fep_chargeA
    fep_chargeB => domain%fep_chargeB
    fep_atmcls_singleB => domain%fep_atmcls_singleB

    ! initializaiton
    !
    ncel = ncel_local + ncel_bd
    HGr_stay(1:connect,1:ncel_local) = 0
    HGr_move(1:connect,1:ncel) = 0

    ! Check outgoing particles
    !
    do i = 1, ncel_local
      do j = 1, connect
        do k = 1, HGr_local(j,i)

          ix = HGr_bond_list(1,k,j,i)
          x_shift = coord(1,ix,i) - boundary%origin_x
          y_shift = coord(2,ix,i) - boundary%origin_y
          z_shift = coord(3,ix,i) - boundary%origin_z

          move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
          move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
          move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)

          x_shift = x_shift + move(1)
          y_shift = y_shift + move(2)
          z_shift = z_shift + move(3)

          !assign which cell
          icx = int(x_shift/csize_x)
          icy = int(y_shift/csize_y)
          icz = int(z_shift/csize_z)
          if (icx == ncel_x) icx = icx - 1
          if (icy == ncel_y) icy = icy - 1
          if (icz == ncel_z) icz = icz - 1
          icel = 1 + icx + icy*ncel_x + icz*ncel_x*ncel_y
          icel_local = cell_g2l(icel)
          icel_bd    = cell_g2b(icel)

          ! atom index
          if (icel_local /= i) then
            if (icel_local /= 0) then

              HGr_move(j,icel_local) = HGr_move(j,icel_local) + 1
              num_move = HGr_move(j,icel_local)

              do list = 1, j+1
                ix = HGr_bond_list(list,k,j,i)
                HGr_move_real(11*(list-1)+ 1,num_move,j,icel_local) = coord(1,ix,i)
                HGr_move_real(11*(list-1)+ 2,num_move,j,icel_local) = coord(2,ix,i)
                HGr_move_real(11*(list-1)+ 3,num_move,j,icel_local) = coord(3,ix,i)
                HGr_move_real(11*(list-1)+ 4,num_move,j,icel_local) = vel(1,ix,i)
                HGr_move_real(11*(list-1)+ 5,num_move,j,icel_local) = vel(2,ix,i)
                HGr_move_real(11*(list-1)+ 6,num_move,j,icel_local) = vel(3,ix,i)
                HGr_move_real(11*(list-1)+ 7,num_move,j,icel_local) = charge(ix,i)
                HGr_move_real(11*(list-1)+ 8,num_move,j,icel_local) = mass(ix,i)
                HGr_move_real(11*(list-1)+ 9,num_move,j,icel_local) = HGr_bond_dist(list,k,j,i)
                HGr_move_real(11*(list-1)+10,num_move,j,icel_local) = fep_chargeA(ix,i)
                HGr_move_real(11*(list-1)+11,num_move,j,icel_local) = fep_chargeB(ix,i)
                HGr_move_int ( 5*(list-1)+ 1,num_move,j,icel_local) = atmcls(ix,i)
                HGr_move_int ( 5*(list-1)+ 2,num_move,j,icel_local) = id_l2g(ix,i)
                HGr_move_int ( 5*(list-1)+ 3,num_move,j,icel_local) = id_l2g_sol(ix,i)
                HGr_move_int ( 5*(list-1)+ 4,num_move,j,icel_local) = fepgrp(ix,i)
                HGr_move_int ( 5*(list-1)+ 5,num_move,j,icel_local) = fep_atmcls_singleB(ix,i)
              end do

            else if (icel_bd /= 0) then

              icel_bd = icel_bd + ncel_local
              HGr_move(j,icel_bd) = HGr_move(j,icel_bd) + 1
              num_move = HGr_move(j,icel_bd)

              do list = 1, j+1
                ix = HGr_bond_list(list,k,j,i)
                hgr_move_real(11*(list-1)+ 1,num_move,j,icel_bd) = coord(1,ix,i)
                hgr_move_real(11*(list-1)+ 2,num_move,j,icel_bd) = coord(2,ix,i)
                hgr_move_real(11*(list-1)+ 3,num_move,j,icel_bd) = coord(3,ix,i)
                hgr_move_real(11*(list-1)+ 4,num_move,j,icel_bd) = vel(1,ix,i)
                hgr_move_real(11*(list-1)+ 5,num_move,j,icel_bd) = vel(2,ix,i)
                hgr_move_real(11*(list-1)+ 6,num_move,j,icel_bd) = vel(3,ix,i)
                hgr_move_real(11*(list-1)+ 7,num_move,j,icel_bd) = charge(ix,i)
                hgr_move_real(11*(list-1)+ 8,num_move,j,icel_bd) = mass(ix,i)
                HGr_move_real(11*(list-1)+ 9,num_move,j,icel_bd) = HGr_bond_dist(list,k,j,i)
                hgr_move_real(11*(list-1)+10,num_move,j,icel_bd) = fep_chargeA(ix,i)
                hgr_move_real(11*(list-1)+11,num_move,j,icel_bd) = fep_chargeB(ix,i)
                hgr_move_int ( 5*(list-1)+ 1,num_move,j,icel_bd) = atmcls(ix,i)
                hgr_move_int ( 5*(list-1)+ 2,num_move,j,icel_bd) = id_l2g(ix,i)
                hgr_move_int ( 5*(list-1)+ 3,num_move,j,icel_bd) = id_l2g_sol(ix,i)
                hgr_move_int ( 5*(list-1)+ 4,num_move,j,icel_bd) = fepgrp(ix,i)
                hgr_move_int ( 5*(list-1)+ 5,num_move,j,icel_bd) = fep_atmcls_singleB(ix,i)
              end do
            end if

          else

            hgr_stay(j,i) = hgr_stay(j,i) + 1
            num_stay = hgr_stay(j,i)

            do list = 1, j+1
              ix = HGr_bond_list(list,k,j,i)
              hgr_stay_real(11*(list-1)+ 1,num_stay,j,i) = coord(1,ix,i)
              hgr_stay_real(11*(list-1)+ 2,num_stay,j,i) = coord(2,ix,i)
              hgr_stay_real(11*(list-1)+ 3,num_stay,j,i) = coord(3,ix,i)
              hgr_stay_real(11*(list-1)+ 4,num_stay,j,i) = vel(1,ix,i)
              hgr_stay_real(11*(list-1)+ 5,num_stay,j,i) = vel(2,ix,i)
              hgr_stay_real(11*(list-1)+ 6,num_stay,j,i) = vel(3,ix,i)
              hgr_stay_real(11*(list-1)+ 7,num_stay,j,i) = charge(ix,i)
              hgr_stay_real(11*(list-1)+ 8,num_stay,j,i) = mass(ix,i)
              hgr_stay_real(11*(list-1)+ 9,num_stay,j,i) = HGr_bond_dist(list,k,j,i)
              hgr_stay_real(11*(list-1)+10,num_stay,j,i) = fep_chargeA(ix,i)
              hgr_stay_real(11*(list-1)+11,num_stay,j,i) = fep_chargeB(ix,i)
              hgr_stay_int ( 5*(list-1)+ 1,num_stay,j,i) = atmcls(ix,i)
              hgr_stay_int ( 5*(list-1)+ 2,num_stay,j,i) = id_l2g(ix,i)
              hgr_stay_int ( 5*(list-1)+ 3,num_stay,j,i) = id_l2g_sol(ix,i)
              hgr_stay_int ( 5*(list-1)+ 4,num_stay,j,i) = fepgrp(ix,i)
              hgr_stay_int ( 5*(list-1)+ 5,num_stay,j,i) = fep_atmcls_singleB(ix,i)
            end do

          end if

        end do
      end do
    end do

    return

  end subroutine update_outgoing_HGr_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_atom_fep
  !> @brief        check particles (not bonded to hydrogen) incoming to each
  !!               cell for FEP
  !! @authors      HO
  !! @param[inout] constraints : constraints information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_incoming_atom_fep(constraints, domain)

    ! formal arguments
    type(s_constraints), target, intent(inout) :: constraints
    type(s_domain),      target, intent(inout) :: domain

    ! local variables
    integer                      :: i, j, k, ix, kx
    logical                      :: insert

    real(wip),        pointer    :: coord(:,:,:), velocity(:,:,:)
    real(wp),         pointer    :: charge(:,:)
    real(wip),        pointer    :: mass(:,:)
    real(wip),        pointer    :: buf_real(:,:,:)
    integer,          pointer    :: ncel_local, nsolute(:)
    integer,          pointer    :: atmcls(:,:), id_l2g(:,:), id_l2g_sol(:,:)
    integer(int2),    pointer    :: id_g2l(:,:)
    integer,          pointer    :: ptl_add(:), ptl_exit(:), ptlindex(:,:)
    integer,          pointer    :: buf_int(:,:,:)

    ! FEP
    integer,          pointer    :: fepgrp(:,:)
    integer,          pointer    :: fep_atmcls_singleB(:,:)
    real(wp),         pointer    :: fep_chargeA(:,:)
    real(wp),         pointer    :: fep_chargeB(:,:)

    nsolute    => constraints%No_HGr

    ncel_local => domain%num_cell_local
    coord      => domain%coord
    velocity   => domain%velocity
    charge     => domain%charge
    mass       => domain%mass
    atmcls     => domain%atom_cls_no
    id_l2g     => domain%id_l2g
    id_l2g_sol => domain%id_l2g_solute
    id_g2l     => domain%id_g2l
    ptl_add    => domain%ptl_add
    ptl_exit   => domain%ptl_exit
    ptlindex   => domain%ptl_exit_index
    buf_int    => domain%buf_integer
    buf_real   => domain%buf_real

    ! FEP
    fepgrp     => domain%fepgrp
    fep_chargeA => domain%fep_chargeA
    fep_chargeB => domain%fep_chargeB
    fep_atmcls_singleB => domain%fep_atmcls_singleB

    ! incoming particles
    !
    do i = 1, ncel_local

#ifdef DEBUG
      if (ptl_add(i)+nsolute(i)-ptl_exit(i) > MaxAtom) &
        call error_msg('Debug: Update_Incoming_Atom> atom is exceed MaxAtom')
#endif

      ! when the number of coming particles is larger than that of outgoing ones
      !
      if (ptl_add(i) >= ptl_exit(i)) then

        do k = 1, ptl_exit(i)
          coord(1,ptlindex(k,i),i)     = buf_real(1,k,i)
          coord(2,ptlindex(k,i),i)     = buf_real(2,k,i)
          coord(3,ptlindex(k,i),i)     = buf_real(3,k,i)
          velocity(1,ptlindex(k,i),i)  = buf_real(4,k,i)
          velocity(2,ptlindex(k,i),i)  = buf_real(5,k,i)
          velocity(3,ptlindex(k,i),i)  = buf_real(6,k,i)
          charge(ptlindex(k,i),i)      = buf_real(7,k,i)
          mass(ptlindex(k,i),i)        = buf_real(8,k,i)
          fep_chargeA(ptlindex(k,i),i) = buf_real(9,k,i)
          fep_chargeB(ptlindex(k,i),i) = buf_real(10,k,i)
          atmcls(ptlindex(k,i),i)      = buf_int(1,k,i)
          id_l2g(ptlindex(k,i),i)      = buf_int(2,k,i)
          id_l2g_sol(ptlindex(k,i),i)  = buf_int(3,k,i)
          id_g2l(1,buf_int(3,k,i))     = i
          id_g2l(2,buf_int(3,k,i))     = ptlindex(k,i)
          fepgrp(ptlindex(k,i),i)      = buf_int(4,k,i)
          fep_atmcls_singleB(ptlindex(k,i),i) = buf_int(5,k,i)
        end do

        do k = ptl_exit(i)+1, ptl_add(i)
          ix = k + nsolute(i) - ptl_exit(i)
          coord(1,ix,i)               = buf_real(1,k,i)
          coord(2,ix,i)               = buf_real(2,k,i)
          coord(3,ix,i)               = buf_real(3,k,i)
          velocity(1,ix,i)            = buf_real(4,k,i)
          velocity(2,ix,i)            = buf_real(5,k,i)
          velocity(3,ix,i)            = buf_real(6,k,i)
          charge(ix,i)                = buf_real(7,k,i)
          mass(ix,i)                  = buf_real(8,k,i)
          fep_chargeA(ix,i)           = buf_real(9,k,i)
          fep_chargeB(ix,i)           = buf_real(10,k,i)
          atmcls(ix,i)                = buf_int (1,k,i)
          id_l2g(ix,i)                = buf_int (2,k,i)
          id_l2g_sol(ix,i)            = buf_int (3,k,i)
          id_g2l(1,buf_int(3,k,i))    = i
          id_g2l(2,buf_int(3,k,i))    = ix
          fepgrp(ix,i)                = buf_int (4,k,i)
          fep_atmcls_singleB(ix,i)    = buf_int (5,k,i)
        end do
        nsolute(i) = nsolute(i) + ptl_add(i) - ptl_exit(i)

      ! when the number of coming particles is less than that of outgoing ones
      !
      else

        do k = 1, ptl_add(i)
          coord(1,ptlindex(k,i),i)     = buf_real(1,k,i)
          coord(2,ptlindex(k,i),i)     = buf_real(2,k,i)
          coord(3,ptlindex(k,i),i)     = buf_real(3,k,i)
          velocity(1,ptlindex(k,i),i)  = buf_real(4,k,i)
          velocity(2,ptlindex(k,i),i)  = buf_real(5,k,i)
          velocity(3,ptlindex(k,i),i)  = buf_real(6,k,i)
          charge(ptlindex(k,i),i)      = buf_real(7,k,i)
          mass(ptlindex(k,i),i)        = buf_real(8,k,i)
          fep_chargeA(ptlindex(k,i),i) = buf_real(9,k,i)
          fep_chargeB(ptlindex(k,i),i) = buf_real(10,k,i)
          atmcls(ptlindex(k,i),i)      = buf_int (1,k,i)
          id_l2g(ptlindex(k,i),i)      = buf_int (2,k,i)
          id_l2g_sol(ptlindex(k,i),i)  = buf_int (3,k,i)
          id_g2l(1,buf_int(3,k,i))     = i
          id_g2l(2,buf_int(3,k,i))     = ptlindex(k,i)
          fepgrp(ptlindex(k,i),i)      = buf_int (4,k,i)
          fep_atmcls_singleB(ptlindex(k,i),i) = buf_int (5,k,i)
        end do

        j  = 0
        ix = nsolute(i)
        k  = ptl_add(i) + 1

        do while (j < (ptl_exit(i)-ptl_add(i)))

          insert = .true.

          do kx = k, ptl_exit(i)
            if (ix == ptlindex(kx,i)) then
              insert = .false.
              j = j + 1
              exit
            end if
          end do

          if (insert) then
            kx = ptlindex(k,i)
            coord(1,kx,i)              = coord(1,ix,i)
            coord(2,kx,i)              = coord(2,ix,i)
            coord(3,kx,i)              = coord(3,ix,i)
            velocity(1,kx,i)           = velocity(1,ix,i)
            velocity(2,kx,i)           = velocity(2,ix,i)
            velocity(3,kx,i)           = velocity(3,ix,i)
            charge(kx,i)               = charge(ix,i)
            mass(kx,i)                 = mass(ix,i)
            fep_chargeA(kx,i)          = fep_chargeA(ix,i)
            fep_chargeB(kx,i)          = fep_chargeB(ix,i)
            atmcls(kx,i)               = atmcls(ix,i)
            id_l2g(kx,i)               = id_l2g(ix,i)
            id_l2g_sol(kx,i)           = id_l2g_sol(ix,i)
            id_g2l(1,id_l2g_sol(kx,i)) = i
            id_g2l(2,id_l2g_sol(kx,i)) = kx
            fepgrp(kx,i)               = fepgrp(ix,i)
            fep_atmcls_singleB(kx,i)   = fep_atmcls_singleB(ix,i)

            j = j + 1
            k = k + 1

          end if
          ix = ix - 1

        end do
        nsolute(i) = nsolute(i) + ptl_add(i) - ptl_exit(i)

      end if
    end do

    return

  end subroutine update_incoming_atom_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_water_fep
  !> @brief        check water particles incoming cells for FEP
  !! @authors      HO
  !! @param[in]    water_atom : number of water atoms
  !! @param[inout] domain     : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_incoming_water_fep(water_atom, domain)

    ! formal arguments
    integer,                 intent(in)    :: water_atom
    type(s_domain),  target, intent(inout) :: domain

    ! local variable
    integer                  :: i, k, iwater, jwater, ix

    real(wip),       pointer :: coord(:,:,:), velocity(:,:,:)
    real(wip),       pointer :: water_move_real(:,:,:), water_stay_real(:,:,:)
    real(wp),        pointer :: water_charge(:), water_mass(:)
    real(wp),        pointer :: charge(:,:)
    real(wip),       pointer :: mass(:,:)
    integer,         pointer :: ncel_local
    integer,         pointer :: nsolute(:), natom(:), nwater(:)
    integer,         pointer :: water_list(:,:,:)
    integer,         pointer :: id_l2g(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: water_move(:), water_stay(:)
    integer,         pointer :: water_move_int(:,:,:), water_stay_int(:,:,:)
    integer,         pointer :: water_atmcls(:), atmcls(:,:)

    ! FEP
    integer,         pointer :: fepgrp(:,:)


    ncel_local      => domain%num_cell_local
    nsolute         => domain%num_solute
    nwater          => domain%num_water
    natom           => domain%num_atom
    water_list      => domain%water_list
    coord           => domain%coord
    velocity        => domain%velocity
    charge          => domain%charge
    mass            => domain%mass
    atmcls          => domain%atom_cls_no
    id_l2g          => domain%id_l2g
    id_g2l          => domain%id_g2l
    water_atmcls    => domain%water%atom_cls_no
    water_charge    => domain%water%charge
    water_mass      => domain%water%mass
    water_move      => domain%water%move
    water_stay      => domain%water%stay
    water_move_real => domain%water%move_real
    water_stay_real => domain%water%stay_real
    water_move_int  => domain%water%move_integer
    water_stay_int  => domain%water%stay_integer

    ! FEP
    fepgrp          => domain%fepgrp


    do i = 1, ncel_local

#ifdef DEBUG
      if (water_stay(i)+water_move(i) > MaxWater) then
        call error_msg('Debug: Update_Incoming_Water> water is exceed MaxWater')
      end if
#endif

      do iwater = 1, water_stay(i)
        do k = 1, water_atom
          water_list(k,iwater,i) = nsolute(i) + water_atom*(iwater-1) + k
          ix = water_list(k,iwater,i)

          coord(1,ix,i)          = water_stay_real(6*(k-1)+1,iwater,i)
          coord(2,ix,i)          = water_stay_real(6*(k-1)+2,iwater,i)
          coord(3,ix,i)          = water_stay_real(6*(k-1)+3,iwater,i)
          velocity(1,ix,i)       = water_stay_real(6*(k-1)+4,iwater,i)
          velocity(2,ix,i)       = water_stay_real(6*(k-1)+5,iwater,i)
          velocity(3,ix,i)       = water_stay_real(6*(k-1)+6,iwater,i)
          id_l2g(ix,i)           = water_stay_int (2*(k-1)+1,iwater,i)
          fepgrp(ix,i)           = water_stay_int (2*(k-1)+2,iwater,i)
          charge(ix,i)           = water_charge(k)
          mass(ix,i)             = water_mass(k)
          atmcls(ix,i)           = water_atmcls(k)
        end do
      end do

      do iwater = 1, water_move(i)
        jwater = iwater + water_stay(i)
        do k = 1, water_atom
          water_list(k,jwater,i) = nsolute(i) + water_atom*(jwater-1) + k
          ix = water_list(k,jwater,i)
          coord(1,ix,i)          = water_move_real(6*(k-1)+1,iwater,i)
          coord(2,ix,i)          = water_move_real(6*(k-1)+2,iwater,i)
          coord(3,ix,i)          = water_move_real(6*(k-1)+3,iwater,i)
          velocity(1,ix,i)       = water_move_real(6*(k-1)+4,iwater,i)
          velocity(2,ix,i)       = water_move_real(6*(k-1)+5,iwater,i)
          velocity(3,ix,i)       = water_move_real(6*(k-1)+6,iwater,i)
          id_l2g(ix,i)           = water_move_int (2*(k-1)+1,iwater,i)
          fepgrp(ix,i)           = water_move_int (2*(k-1)+2,iwater,i)
          charge(ix,i)           = water_charge(k)
          mass(ix,i)             = water_mass(k)
          atmcls(ix,i)           = water_atmcls(k)
        end do
      end do

      nwater(i) = water_stay(i) + water_move(i)
      natom(i)  = nsolute(i) + water_atom*nwater(i)
    end do

    return

  end subroutine update_incoming_water_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_HGr_fep
  !> @brief        check hydrogen-bonded particles incoming cells for FEP
  !! @authors      HO
  !! @param[inout] constraints : constraints information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_incoming_HGr_fep(constraints, domain)

    ! formal variable
    type(s_constraints), target, intent(inout)  :: constraints
    type(s_domain),      target, intent(inout)  :: domain

    ! local variable
    integer                      :: i, j, k, k1, list, ix, connect
    integer                      :: num_atom

    real(wip),        pointer    :: coord(:,:,:), vel(:,:,:)
    real(wip),        pointer    :: hgr_move_real(:,:,:,:)
    real(wip),        pointer    :: hgr_stay_real(:,:,:,:)
    real(wp),         pointer    :: charge(:,:)
    real(wip),        pointer    :: mass(:,:)
    real(wip),        pointer    :: HGr_bond_dist(:,:,:,:)
    integer,          pointer    :: ncel_local
    integer,          pointer    :: nsolute(:)
    integer,          pointer    :: id_l2g(:,:), id_l2g_sol(:,:), atmcls(:,:)
    integer(int2),    pointer    :: id_g2l(:,:)
    integer,          pointer    :: No_HGr(:), HGr_local(:,:)
    integer,          pointer    :: HGr_bond_list(:,:,:,:)
    integer,          pointer    :: HGr_move(:,:), HGr_stay(:,:)
    integer,          pointer    :: HGr_move_int(:,:,:,:), HGr_stay_int(:,:,:,:)

    ! FEP
    integer,          pointer    :: fepgrp(:,:)
    integer,          pointer    :: fep_atmcls_singleB(:,:)
    real(wp),         pointer    :: fep_chargeA(:,:)
    real(wp),         pointer    :: fep_chargeB(:,:)

    No_HGr        => constraints%no_hgr
    HGr_local     => constraints%HGr_local
    HGr_bond_list => constraints%HGr_bond_list
    HGr_move      => constraints%HGr_move
    HGr_stay      => constraints%HGr_stay
    HGr_bond_dist => constraints%HGr_bond_dist
    HGr_move_real => constraints%HGr_move_real
    HGr_stay_real => constraints%HGr_stay_real
    HGr_move_int  => constraints%HGr_move_int
    HGr_stay_int  => constraints%HGr_stay_int
    connect       =  constraints%connect

    ncel_local    => domain%num_cell_local
    nsolute       => domain%num_solute
    coord         => domain%coord
    vel           => domain%velocity
    charge        => domain%charge
    mass          => domain%mass
    atmcls        => domain%atom_cls_no
    id_l2g        => domain%id_l2g
    id_l2g_sol    => domain%id_l2g_solute
    id_g2l        => domain%id_g2l

    ! FEP
    fepgrp        => domain%fepgrp
    fep_chargeA => domain%fep_chargeA
    fep_chargeB => domain%fep_chargeB
    fep_atmcls_singleB => domain%fep_atmcls_singleB

    do i = 1, ncel_local

      num_atom = No_HGr(i)

      do j = 1, connect

#ifdef DEBUG
        if (HGr_stay(j,i)+Hgr_move(j,i) > HGroupMax) then
          call error_msg('Debug: Update_Incoming_HGr> HGr is exceed HGroupMax')
        end if
#endif

        do k = 1, HGr_stay(j,i)
          do list = 1, j+1
            num_atom = num_atom + 1
            HGr_bond_list(list,k,j,i) = num_atom
            ix = num_atom
            coord(1,ix,i)             = HGr_stay_real(11*(list-1)+1,k,j,i)
            coord(2,ix,i)             = HGr_stay_real(11*(list-1)+2,k,j,i)
            coord(3,ix,i)             = HGr_stay_real(11*(list-1)+3,k,j,i)
            vel(1,ix,i)               = HGr_stay_real(11*(list-1)+4,k,j,i)
            vel(2,ix,i)               = HGr_stay_real(11*(list-1)+5,k,j,i)
            vel(3,ix,i)               = HGr_stay_real(11*(list-1)+6,k,j,i)
            charge(ix,i)              = HGr_stay_real(11*(list-1)+7,k,j,i)
            mass(ix,i)                = HGr_stay_real(11*(list-1)+8,k,j,i)
            HGr_bond_dist(list,k,j,i) = HGr_stay_real(11*(list-1)+9,k,j,i)
            fep_chargeA(ix,i)         = HGr_stay_real(11*(list-1)+10,k,j,i)
            fep_chargeB(ix,i)         = HGr_stay_real(11*(list-1)+11,k,j,i)
            atmcls(ix,i)              = HGr_stay_int (5*(list-1)+1,k,j,i)
            id_l2g(ix,i)              = HGr_stay_int (5*(list-1)+2,k,j,i)
            id_l2g_sol(ix,i)          = HGr_stay_int (5*(list-1)+3,k,j,i)
            id_g2l(1,id_l2g_sol(ix,i))= i
            id_g2l(2,id_l2g_sol(ix,i))= ix
            fepgrp(ix,i)              = HGr_stay_int (5*(list-1)+4,k,j,i)
            fep_atmcls_singleB(ix,i)  = HGr_stay_int (5*(list-1)+5,k,j,i)
          end do
        end do

        do k = 1, HGr_move(j,i)
          k1 = k + HGr_stay(j,i)
          do list = 1, j+1
            num_atom = num_atom + 1
            HGr_bond_list(list,k1,j,i) = num_atom
            ix = num_atom
            coord(1,ix,i)             = HGr_move_real(11*(list-1)+1,k,j,i)
            coord(2,ix,i)             = HGr_move_real(11*(list-1)+2,k,j,i)
            coord(3,ix,i)             = HGr_move_real(11*(list-1)+3,k,j,i)
            vel(1,ix,i)               = HGr_move_real(11*(list-1)+4,k,j,i)
            vel(2,ix,i)               = HGr_move_real(11*(list-1)+5,k,j,i)
            vel(3,ix,i)               = HGr_move_real(11*(list-1)+6,k,j,i)
            charge(ix,i)              = HGr_move_real(11*(list-1)+7,k,j,i)
            mass(ix,i)                = HGr_move_real(11*(list-1)+8,k,j,i)
            HGr_bond_dist(list,k1,j,i)= HGr_move_real(11*(list-1)+9,k,j,i)
            fep_chargeA(ix,i)         = HGr_move_real(11*(list-1)+10,k,j,i)
            fep_chargeB(ix,i)         = HGr_move_real(11*(list-1)+11,k,j,i)
            atmcls(ix,i)              = HGr_move_int (5*(list-1)+1,k,j,i)
            id_l2g(ix,i)              = HGr_move_int (5*(list-1)+2,k,j,i)
            id_l2g_sol(ix,i)          = HGr_move_int (5*(list-1)+3,k,j,i)
            id_g2l(1,id_l2g_sol(ix,i))= i
            id_g2l(2,id_l2g_sol(ix,i))= ix
            fepgrp(ix,i)              = HGr_move_int (5*(list-1)+4,k,j,i)
            fep_atmcls_singleB(ix,i)  = HGr_move_int (5*(list-1)+5,k,j,i)
          end do
        end do

        HGr_local(j,i) = HGr_stay(j,i) + HGr_move(j,i)

      end do

      nsolute(i) = num_atom

    end do

    return

  end subroutine update_incoming_HGr_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_enefunc_bond_fep
  !> @brief        update BOND term for each cell in potential energy function
  !                for FEP
  !! @authors      HO
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_enefunc_bond_fep(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    integer                  :: i, ix, icel_local
    integer                  :: i1, i2, icel1, icel2

    real(wp),        pointer :: fc(:,:), r0(:,:)
    real(wp),        pointer :: buf_real(:,:,:)
    integer,         pointer :: ncel_local, nboundary
    integer(int2),   pointer :: cell_pair(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: nbond(:), bondlist(:,:,:)
    integer(1),      pointer :: bondkind(:,:)
    integer,         pointer :: bond_pbc(:,:)
    integer,         pointer :: bond_exit(:), index(:,:)
    integer,         pointer :: bond_add(:), buf_int(:,:,:)

    ! FEP
    integer,         pointer :: bond_singleB(:,:)

    ncel_local  => domain%num_cell_local
    nboundary   => domain%num_cell_boundary
    cell_pair   => domain%cell_pair
    id_g2l      => domain%id_g2l

    nbond       => enefunc%num_bond
    bondlist    => enefunc%bond_list
    bondkind    => enefunc%bond_kind
    bond_pbc    => enefunc%bond_pbc
    fc          => enefunc%bond_force_const
    r0          => enefunc%bond_dist_min
    bond_exit   => enefunc%bond_exit
    index       => enefunc%bond_exit_index
    bond_add    => enefunc%bond_add
    buf_int     => enefunc%buf_bond_integer
    buf_real    => enefunc%buf_bond_real

    ! FEP
    bond_singleB => enefunc%bond_singleB

    bond_exit(1:ncel_local) = 0
    bond_add(1:ncel_local+nboundary) = 0

    do i = 1, ncel_local
      do ix = 1, nbond(i)

        icel1 = id_g2l(1,bondlist(1,ix,i))
        i1    = id_g2l(2,bondlist(1,ix,i))
        icel2 = id_g2l(1,bondlist(2,ix,i))
        i2    = id_g2l(2,bondlist(2,ix,i))

        icel_local = cell_pair(icel1,icel2)

        if (icel_local /= i) then

          bond_exit(i) = bond_exit(i) + 1
          index(bond_exit(i),i) = ix

          bond_add(icel_local) = bond_add(icel_local) + 1
          buf_real(1,bond_add(icel_local),icel_local) = r0(ix,i)
          buf_real(2,bond_add(icel_local),icel_local) = fc(ix,i)
          buf_int (1,bond_add(icel_local),icel_local) = bondlist(1,ix,i)
          buf_int (2,bond_add(icel_local),icel_local) = bondlist(2,ix,i)
          buf_int (3,bond_add(icel_local),icel_local) = bondkind(ix,i)
          buf_int (4,bond_add(icel_local),icel_local) = bond_pbc(ix,i)
          buf_int (5,bond_add(icel_local),icel_local) = bond_singleB(ix,i)

        end if

      end do
    end do

    return

  end subroutine update_outgoing_enefunc_bond_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_enefunc_bond_fep
  !> @brief        update BOND term for each cell in potential energy function
  !                for FEP
  !! @authors      HO
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_incoming_enefunc_bond_fep(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    integer                  :: i, j, k, ix, kx, found
    logical                  :: insert

    real(wp),        pointer :: buf_real(:,:,:), bonddist(:,:), bondforce(:,:)
    integer,         pointer :: ncel_local, nboundary
    integer(int2),   pointer :: cell_g2l(:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: bond_add(:), bond_exit(:), bond_exit_index(:,:)
    integer,         pointer :: buf_int(:,:,:), nbond(:)
    integer,         pointer :: bondlist(:,:,:)
    integer(1),      pointer :: bondkind(:,:)
    integer,         pointer :: bond_pbc(:,:)

    ! FEP
    integer,         pointer :: bond_singleB(:,:)

    ncel_local      => domain%num_cell_local
    nboundary       => domain%num_cell_boundary
    cell_g2l        => domain%cell_g2l
    id_g2l          => domain%id_g2l

    bond_add        => enefunc%bond_add
    bond_exit       => enefunc%bond_exit
    bond_exit_index => enefunc%bond_exit_index
    buf_int         => enefunc%buf_bond_integer
    buf_real        => enefunc%buf_bond_real
    bondlist        => enefunc%bond_list
    bondforce       => enefunc%bond_force_const
    bonddist        => enefunc%bond_dist_min
    nbond           => enefunc%num_bond
    bondkind        => enefunc%bond_kind
    bond_pbc        => enefunc%bond_pbc

    ! FEP
    bond_singleB => enefunc%bond_singleB

    found = 0
    do i = 1, ncel_local

#ifdef DEBUG
      if (bond_add(i)+nbond(i)-bond_exit(i) > MaxBond) then
        call error_msg( &
            'Debug: update_incoming_enefunc_bond> bond is exceed MaxBond')
      end if
#endif

      if (bond_add(i) > BondMove) &
        call error_msg( &
        'Debug: update_incoming_enefunc_bond> bond migration is exceed maximum')

      if (bond_add(i) >= bond_exit(i)) then

        do k = 1, bond_exit(i)
          bondlist (1,bond_exit_index(k,i),i) = buf_int (1,k,i)
          bondlist (2,bond_exit_index(k,i),i) = buf_int (2,k,i)
          bondkind(bond_exit_index(k,i),i)    = buf_int (3,k,i)
          bond_pbc(bond_exit_index(k,i),i)    = buf_int (4,k,i)
          bond_singleB(bond_exit_index(k,i),i) = buf_int (5,k,i)
          bonddist (bond_exit_index(k,i),i)   = buf_real(1,k,i)
          bondforce(bond_exit_index(k,i),i)   = buf_real(2,k,i)
        end do

        do k = bond_exit(i)+1, bond_add(i)
          ix = k + nbond(i) - bond_exit(i)
          bondlist (1,ix,i) = buf_int(1,k,i)
          bondlist (2,ix,i) = buf_int(2,k,i)
          bondkind (ix,i)   = buf_int(3,k,i)
          bond_pbc (ix,i)   = buf_int(4,k,i)
          bond_singleB(ix,i) = buf_int(5,k,i)
          bonddist (ix,i)   = buf_real(1,k,i)
          bondforce(ix,i)   = buf_real(2,k,i)
        end do

        nbond(i) = nbond(i) + bond_add(i) - bond_exit(i)

      else

        do k = 1, bond_add(i)
          bondlist (1,bond_exit_index(k,i),i) = buf_int (1,k,i)
          bondlist (2,bond_exit_index(k,i),i) = buf_int (2,k,i)
          bondkind(bond_exit_index(k,i),i)    = buf_int (3,k,i)
          bond_pbc(bond_exit_index(k,i),i)    = buf_int (4,k,i)
          bond_singleB(bond_exit_index(k,i),i) = buf_int (5,k,i)
          bonddist (bond_exit_index(k,i),i)   = buf_real(1,k,i)
          bondforce(bond_exit_index(k,i),i)   = buf_real(2,k,i)
        end do

        j  = 0
        ix = nbond(i)
        k  = bond_add(i) + 1

        do while (j < (bond_exit(i)-bond_add(i)))

          insert = .true.
          do kx = k, bond_exit(i)
            if (ix == bond_exit_index(kx,i)) then
              insert = .false.
              j = j + 1
              exit
            end if
          end do

          if (insert) then
            kx = bond_exit_index(k,i)
            bondlist (1,kx,i) = bondlist(1,ix,i)
            bondlist (2,kx,i) = bondlist(2,ix,i)
            bondkind (kx,i)   = bondkind(ix,i)
            bond_pbc (kx,i)   = bond_pbc(ix,i)
            bond_singleB(kx,i) = bond_singleB(ix,i)
            bonddist (kx,i)   = bonddist(ix,i)
            bondforce(kx,i)   = bondforce(ix,i)
            j = j + 1
            k = k + 1
          end if
          ix = ix - 1

        end do

        nbond(i) = nbond(i) + bond_add(i) - bond_exit(i)

      end if

    end do

    return

  end subroutine update_incoming_enefunc_bond_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_enefunc_angl_fep
  !> @brief        update ANGLE term for each cell in potential energy function
  !                for FEP
  !! @authors      HO
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_enefunc_angl_fep(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: icel1, icel2, icel3, i1, i2, i3
    integer                  :: i, ix, icel_local

    real(wp),        pointer :: fc(:,:), theta0(:,:), fc_ub(:,:), r0_ub(:,:)
    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: buf_real(:,:,:)
    integer,         pointer :: ncel_local, nboundary
    integer(int2),   pointer :: cell_g2l(:), cell_g2b(:)
    integer(int2),   pointer :: cell_pair(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: nangle(:), anglelist(:,:,:), angl_pbc(:,:,:)
    integer(1),      pointer :: anglekind(:,:)
    integer,         pointer :: angle_add(:), angle_exit(:)
    integer,         pointer :: angle_exit_index(:,:), buf_int(:,:,:)

    ! FEP
    integer,         pointer :: angl_singleB(:,:)

    ncel_local       => domain%num_cell_local
    nboundary        => domain%num_cell_boundary
    cell_g2l         => domain%cell_g2l
    cell_g2b         => domain%cell_g2b
    cell_pair        => domain%cell_pair
    id_g2l           => domain%id_g2l
    coord            => domain%coord

    nangle           => enefunc%num_angle
    anglelist        => enefunc%angle_list
    anglekind        => enefunc%angle_kind
    angl_pbc         => enefunc%angle_pbc
    fc               => enefunc%angle_force_const
    theta0           => enefunc%angle_theta_min
    fc_ub            => enefunc%urey_force_const
    r0_ub            => enefunc%urey_rmin
    angle_add        => enefunc%angle_add
    angle_exit       => enefunc%angle_exit
    angle_exit_index => enefunc%angle_exit_index
    buf_int          => enefunc%buf_angle_integer
    buf_real         => enefunc%buf_angle_real

    ! FEP
    angl_singleB     => enefunc%angl_singleB

    angle_add (1:ncel_local+nboundary) = 0
    angle_exit(1:ncel_local) = 0

    do i = 1, ncel_local
      do ix = 1, nangle(i)

        icel1 = id_g2l(1,anglelist(1,ix,i))
        i1    = id_g2l(2,anglelist(1,ix,i))
        icel2 = id_g2l(1,anglelist(2,ix,i))
        i2    = id_g2l(2,anglelist(2,ix,i))
        icel3 = id_g2l(1,anglelist(3,ix,i))
        i3    = id_g2l(2,anglelist(3,ix,i))

        icel_local = cell_pair(icel1,icel3)

        if (icel_local /= i) then

          angle_exit(i) = angle_exit(i) + 1
          angle_exit_index(angle_exit(i),i) = ix
          angle_add(icel_local) = angle_add(icel_local) + 1
          buf_real(1,angle_add(icel_local),icel_local) = theta0(ix,i)
          buf_real(2,angle_add(icel_local),icel_local) = fc(ix,i)
          buf_real(3,angle_add(icel_local),icel_local) = r0_ub(ix,i)
          buf_real(4,angle_add(icel_local),icel_local) = fc_ub(ix,i)
          buf_int (1,angle_add(icel_local),icel_local) = anglelist(1,ix,i)
          buf_int (2,angle_add(icel_local),icel_local) = anglelist(2,ix,i)
          buf_int (3,angle_add(icel_local),icel_local) = anglelist(3,ix,i)
          buf_int (4,angle_add(icel_local),icel_local) = anglekind(ix,i)
          buf_int (5,angle_add(icel_local),icel_local) = angl_pbc(1,ix,i)
          buf_int (6,angle_add(icel_local),icel_local) = angl_pbc(2,ix,i)
          buf_int (7,angle_add(icel_local),icel_local) = angl_pbc(3,ix,i)
          buf_int (8,angle_add(icel_local),icel_local) = angl_singleB(ix,i)
         
        end if

      end  do
    end do

    return

  end subroutine update_outgoing_enefunc_angl_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_enefunc_angl_fep
  !> @brief        update ANGLE term for each cell in potential energy function
  !                for FEP
  !! @authors      HO
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_incoming_enefunc_angl_fep(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    integer                  :: i, j, k, ix, kx
    logical                  :: insert

    real(wp),        pointer :: buf_real(:,:,:)
    real(wp),        pointer :: angletheta(:,:), angleforce(:,:)
    real(wp),        pointer :: ubrmin(:,:), ubforce(:,:)
    integer,         pointer :: ncel_local, nboundary
    integer(int2),   pointer :: cell_g2l(:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: angle_add(:), angle_exit(:)
    integer,         pointer :: angle_exit_index(:,:)
    integer(1),      pointer :: anglekind(:,:)
    integer,         pointer :: buf_int(:,:,:), nangle(:), anglelist(:,:,:)
    integer,         pointer :: angl_pbc(:,:,:)

    ! FEP
    integer,         pointer :: angl_singleB(:,:)

    ncel_local       => domain%num_cell_local
    nboundary        => domain%num_cell_boundary
    cell_g2l         => domain%cell_g2l
    id_g2l           => domain%id_g2l

    nangle           => enefunc%num_angle
    angle_add        => enefunc%angle_add
    angle_exit       => enefunc%angle_exit
    angle_exit_index => enefunc%angle_exit_index
    buf_int          => enefunc%buf_angle_integer
    buf_real         => enefunc%buf_angle_real
    anglelist        => enefunc%angle_list
    anglekind        => enefunc%angle_kind
    angl_pbc         => enefunc%angle_pbc 
    angleforce       => enefunc%angle_force_const
    angletheta       => enefunc%angle_theta_min
    ubforce          => enefunc%urey_force_const
    ubrmin           => enefunc%urey_rmin

    ! FEP
    angl_singleB     => enefunc%angl_singleB

    do i = 1, ncel_local

#ifdef DEBUG
      if (angle_add(i)+nangle(i)-angle_exit(i) > MaxAngle) &
        call error_msg( &
            'Debug: update_incoming_enefunc_angle> angle is exceed MaxAngle')
#endif

      if (angle_add(i) > AngleMove) &
        call error_msg( &
     'Debug: update_incoming_enefunc_angle> angle migration is exceed maximum')

      if (angle_add(i) >= angle_exit(i)) then

        do k = 1, angle_exit(i)
          anglelist (1,angle_exit_index(k,i),i) = buf_int (1,k,i)
          anglelist (2,angle_exit_index(k,i),i) = buf_int (2,k,i)
          anglelist (3,angle_exit_index(k,i),i) = buf_int (3,k,i)
          anglekind (angle_exit_index(k,i),i)   = buf_int (4,k,i)
          angl_pbc  (1,angle_exit_index(k,i),i) = buf_int (5,k,i)
          angl_pbc  (2,angle_exit_index(k,i),i) = buf_int (6,k,i)
          angl_pbc  (3,angle_exit_index(k,i),i) = buf_int (7,k,i)
          angl_singleB(angle_exit_index(k,i),i) = buf_int (8,k,i)
          angletheta(angle_exit_index(k,i),i)   = buf_real(1,k,i)
          angleforce(angle_exit_index(k,i),i)   = buf_real(2,k,i)
          ubrmin    (angle_exit_index(k,i),i)   = buf_real(3,k,i)
          ubforce   (angle_exit_index(k,i),i)   = buf_real(4,k,i)
        end do

        do k = angle_exit(i)+1, angle_add(i)
          ix = k + nangle(i) - angle_exit(i)
          anglelist (1,ix,i) = buf_int(1,k,i)
          anglelist (2,ix,i) = buf_int(2,k,i)
          anglelist (3,ix,i) = buf_int(3,k,i)
          anglekind (ix,i)   = buf_int(4,k,i)
          angl_pbc  (1,ix,i) = buf_int(5,k,i)
          angl_pbc  (2,ix,i) = buf_int(6,k,i)
          angl_pbc  (3,ix,i) = buf_int(7,k,i)
          angl_singleB(ix,i) = buf_int(8,k,i)
          angletheta(ix,i)   = buf_real(1,k,i)
          angleforce(ix,i)   = buf_real(2,k,i)
          ubrmin    (ix,i)   = buf_real(3,k,i)
          ubforce   (ix,i)   = buf_real(4,k,i)
        end do
        nangle(i) = nangle(i) + angle_add(i) - angle_exit(i)

      else
        do k = 1, angle_add(i)
          anglelist (1,angle_exit_index(k,i),i) = buf_int (1,k,i)
          anglelist (2,angle_exit_index(k,i),i) = buf_int (2,k,i)
          anglelist (3,angle_exit_index(k,i),i) = buf_int (3,k,i)
          anglekind (angle_exit_index(k,i),i)   = buf_int (4,k,i)
          angl_pbc  (1,angle_exit_index(k,i),i) = buf_int (5,k,i)
          angl_pbc  (2,angle_exit_index(k,i),i) = buf_int (6,k,i)
          angl_pbc  (3,angle_exit_index(k,i),i) = buf_int (7,k,i)
          angl_singleB(angle_exit_index(k,i),i) = buf_int (8,k,i)
          angletheta(angle_exit_index(k,i),i)   = buf_real(1,k,i)
          angleforce(angle_exit_index(k,i),i)   = buf_real(2,k,i)
          ubrmin    (angle_exit_index(k,i),i)   = buf_real(3,k,i)
          ubforce   (angle_exit_index(k,i),i)   = buf_real(4,k,i)
        end do

        j  = 0
        ix = nangle(i)
        k  = angle_add(i) + 1

        do while (j < (angle_exit(i)-angle_add(i)))

          insert = .true.

          do kx = k, angle_exit(i)
            if (ix == angle_exit_index(kx,i)) then
              insert = .false.
              j = j + 1
              exit
            end if
          end do

          if (insert) then
            kx = angle_exit_index(k,i)
            anglelist (1,kx,i) = anglelist(1,ix,i)
            anglelist (2,kx,i) = anglelist(2,ix,i)
            anglelist (3,kx,i) = anglelist(3,ix,i)
            anglekind (kx,i)   = anglekind(ix,i)
            angl_pbc  (1,kx,i) = angl_pbc (1,ix,i)
            angl_pbc  (2,kx,i) = angl_pbc (2,ix,i)
            angl_pbc  (3,kx,i) = angl_pbc (3,ix,i)
            angl_singleB(kx,i) = angl_singleB(ix,i)
            angletheta(kx,i)   = angletheta(ix,i)
            angleforce(kx,i)   = angleforce(ix,i)
            ubrmin    (kx,i)   = ubrmin(ix,i)
            ubforce   (kx,i)   = ubforce(ix,i)
            j = j + 1
            k = k + 1

          end if
          ix = ix - 1

        end do

        nangle(i) = nangle(i) + angle_add(i) - angle_exit(i)

      end if

    end do

    return

  end subroutine update_incoming_enefunc_angl_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_enefunc_dihe_fep
  !> @brief        update DIHEDRAL term for each cell in potential energy
  !!               function for FEP
  !! @authors      HO
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_enefunc_dihe_fep(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: icel1, icel2, icel3, icel4, i1, i2, i3, i4
    integer                  :: i, ix, icel_local

    real(wp),        pointer :: fc(:,:), phase(:,:)
    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: buf_real(:,:,:)
    integer,         pointer :: ncel_local, nboundary
    integer(int2),   pointer :: cell_g2l(:), cell_g2b(:)
    integer(int2),   pointer :: cell_pair(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: dihed_add(:), dihed_exit(:)
    integer,         pointer :: dihed_exit_index(:,:)
    integer,         pointer :: ndihe(:), dihelist(:,:,:), dihe_pbc(:,:,:)
    integer(1),      pointer :: dihekind(:,:)
    integer,         pointer :: nperiod(:,:), buf_int(:,:,:)

    ! FEP
    integer,         pointer :: dihe_singleB(:,:)

    ncel_local       => domain%num_cell_local
    nboundary        => domain%num_cell_boundary
    cell_g2l         => domain%cell_g2l
    cell_g2b         => domain%cell_g2b
    cell_pair        => domain%cell_pair
    id_g2l           => domain%id_g2l
    coord            => domain%coord

    ndihe            => enefunc%num_dihedral
    dihelist         => enefunc%dihe_list
    dihekind         => enefunc%dihe_kind
    dihe_pbc         => enefunc%dihe_pbc 
    dihed_add        => enefunc%dihed_add
    dihed_exit       => enefunc%dihed_exit
    dihed_exit_index => enefunc%dihed_exit_index
    fc               => enefunc%dihe_force_const
    nperiod          => enefunc%dihe_periodicity
    phase            => enefunc%dihe_phase
    buf_int          => enefunc%buf_dihed_integer
    buf_real         => enefunc%buf_dihed_real

    ! FEP
    dihe_singleB     => enefunc%dihe_singleB

    dihed_add (1:ncel_local+nboundary) = 0
    dihed_exit(1:ncel_local) = 0


    do i = 1, ncel_local
      do ix = 1, ndihe(i)

        icel1 = id_g2l(1,dihelist(1,ix,i))
        i1    = id_g2l(2,dihelist(1,ix,i))
        icel2 = id_g2l(1,dihelist(2,ix,i))
        i2    = id_g2l(2,dihelist(2,ix,i))
        icel3 = id_g2l(1,dihelist(3,ix,i))
        i3    = id_g2l(2,dihelist(3,ix,i))
        icel4 = id_g2l(1,dihelist(4,ix,i))
        i4    = id_g2l(2,dihelist(4,ix,i))

        icel_local = cell_pair(icel1,icel4)

        if (icel_local /= i) then

          dihed_exit(i) = dihed_exit(i) + 1
          dihed_exit_index(dihed_exit(i),i) = ix
          dihed_add(icel_local) = dihed_add(icel_local) + 1
          buf_real(1,dihed_add(icel_local),icel_local) = phase(ix,i)
          buf_real(2,dihed_add(icel_local),icel_local) = fc(ix,i)
          buf_int (1,dihed_add(icel_local),icel_local) = dihelist(1,ix,i)
          buf_int (2,dihed_add(icel_local),icel_local) = dihelist(2,ix,i)
          buf_int (3,dihed_add(icel_local),icel_local) = dihelist(3,ix,i)
          buf_int (4,dihed_add(icel_local),icel_local) = dihelist(4,ix,i)
          buf_int (5,dihed_add(icel_local),icel_local) = nperiod(ix,i)
          buf_int (6,dihed_add(icel_local),icel_local) = dihekind(ix,i)
          buf_int (7,dihed_add(icel_local),icel_local) = dihe_pbc(1,ix,i)
          buf_int (8,dihed_add(icel_local),icel_local) = dihe_pbc(2,ix,i)
          buf_int (9,dihed_add(icel_local),icel_local) = dihe_pbc(3,ix,i)
          buf_int (10,dihed_add(icel_local),icel_local) = dihe_singleB(ix,i)

        end if

      end  do
    end do

    return

  end subroutine update_outgoing_enefunc_dihe_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_enefunc_dihe_fep
  !> @brief        update DIHEDRAL term for each cell in potential energy
  !!               function for FEP
  !! @authors      HO
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_incoming_enefunc_dihe_fep(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    integer                  :: i, j, k, ix, kx, ndihe
    logical                  :: insert

    real(wp),        pointer :: dihephase(:,:), diheforce(:,:), buf_real(:,:,:)
    integer,         pointer :: ncel_local, nboundary
    integer(int2),   pointer :: cell_g2l(:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: dihed_add(:), dihed_exit(:)
    integer,         pointer :: dihed_exit_index(:,:), buf_int(:,:,:)
    integer,         pointer :: ndihedral(:), dihelist(:,:,:), diheperio(:,:)
    integer(1),      pointer :: dihekind(:,:)
    integer,         pointer :: dihe_pbc(:,:,:)

    ! FEP
    integer,         pointer :: dihe_singleB(:,:)

    ncel_local       => domain%num_cell_local
    nboundary        => domain%num_cell_boundary
    cell_g2l         => domain%cell_g2l
    id_g2l           => domain%id_g2l

    ndihedral        => enefunc%num_dihedral
    dihelist         => enefunc%dihe_list
    dihekind         => enefunc%dihe_kind
    dihe_pbc         => enefunc%dihe_pbc 
    diheforce        => enefunc%dihe_force_const
    diheperio        => enefunc%dihe_periodicity
    dihephase        => enefunc%dihe_phase
    dihed_add        => enefunc%dihed_add
    dihed_exit       => enefunc%dihed_exit
    dihed_exit_index => enefunc%dihed_exit_index
    buf_int          => enefunc%buf_dihed_integer
    buf_real         => enefunc%buf_dihed_real

    ! FEP
    dihe_singleB     => enefunc%dihe_singleB

    ndihe = enefunc%num_dihe_all

    do i = 1, ncel_local

#ifdef DEBUG
      if (dihed_add(i)+ndihedral(i)-dihed_exit(i) > MaxDihe) &
        call error_msg( &
            'Debug: update_incoming_enefunc_dihed> dihed is exceed MaxDihe')
#endif

      if (dihed_add(i) > DiheMove) &
        call error_msg( &
    'Debug: update_incoming_enefunc_dihed> dihed migration is exceed maximum')

      if (dihed_add(i) >= dihed_exit(i)) then

        do k = 1, dihed_exit(i)
          dihelist (1,dihed_exit_index(k,i),i) = buf_int (1,k,i)
          dihelist (2,dihed_exit_index(k,i),i) = buf_int (2,k,i)
          dihelist (3,dihed_exit_index(k,i),i) = buf_int (3,k,i)
          dihelist (4,dihed_exit_index(k,i),i) = buf_int (4,k,i)
          diheperio(dihed_exit_index(k,i),i)   = buf_int (5,k,i)
          dihekind (dihed_exit_index(k,i),i)   = buf_int (6,k,i)
          dihe_pbc (1,dihed_exit_index(k,i),i) = buf_int (7,k,i)
          dihe_pbc (2,dihed_exit_index(k,i),i) = buf_int (8,k,i)
          dihe_pbc (3,dihed_exit_index(k,i),i) = buf_int (9,k,i)
          dihe_singleB(dihed_exit_index(k,i),i) = buf_int(10,k,i)
          dihephase(dihed_exit_index(k,i),i)   = buf_real(1,k,i)
          diheforce(dihed_exit_index(k,i),i)   = buf_real(2,k,i)
        end do

        do k = dihed_exit(i)+1, dihed_add(i)
          ix = k + ndihedral(i) - dihed_exit(i)
          dihelist (1,ix,i) = buf_int (1,k,i)
          dihelist (2,ix,i) = buf_int (2,k,i)
          dihelist (3,ix,i) = buf_int (3,k,i)
          dihelist (4,ix,i) = buf_int (4,k,i)
          diheperio(ix,i)   = buf_int (5,k,i)
          dihekind (ix,i)   = buf_int (6,k,i)
          dihe_pbc (1,ix,i) = buf_int (7,k,i)
          dihe_pbc (2,ix,i) = buf_int (8,k,i)
          dihe_pbc (3,ix,i) = buf_int (9,k,i)
          dihe_singleB(ix,i) = buf_int(10,k,i)
          dihephase(ix,i)   = buf_real(1,k,i)
          diheforce(ix,i)   = buf_real(2,k,i)
        end do

        ndihedral(i) = ndihedral(i) + dihed_add(i) - dihed_exit(i)

      else
        do k = 1, dihed_add(i)
          dihelist (1,dihed_exit_index(k,i),i) = buf_int (1,k,i)
          dihelist (2,dihed_exit_index(k,i),i) = buf_int (2,k,i)
          dihelist (3,dihed_exit_index(k,i),i) = buf_int (3,k,i)
          dihelist (4,dihed_exit_index(k,i),i) = buf_int (4,k,i)
          diheperio(dihed_exit_index(k,i),i)   = buf_int (5,k,i)
          dihekind (dihed_exit_index(k,i),i)   = buf_int (6,k,i)
          dihe_pbc (1,dihed_exit_index(k,i),i) = buf_int (7,k,i)
          dihe_pbc (2,dihed_exit_index(k,i),i) = buf_int (8,k,i)
          dihe_pbc (3,dihed_exit_index(k,i),i) = buf_int (9,k,i)
          dihe_singleB(dihed_exit_index(k,i),i) = buf_int (10,k,i)
          dihephase(dihed_exit_index(k,i),i)   = buf_real(1,k,i)
          diheforce(dihed_exit_index(k,i),i)   = buf_real(2,k,i)
        end do

        j  = 0
        ix = ndihedral(i)
        k  = dihed_add(i) + 1

        do while (j < (dihed_exit(i)-dihed_add(i)))

          insert = .true.
          do kx = k, dihed_exit(i)
            if (ix == dihed_exit_index(kx,i)) then
              insert = .false.
              j = j + 1
              exit
            end if
          end do

          if (insert) then
            kx = dihed_exit_index(k,i)
            dihelist (1,kx,i) = dihelist (1,ix,i)
            dihelist (2,kx,i) = dihelist (2,ix,i)
            dihelist (3,kx,i) = dihelist (3,ix,i)
            dihelist (4,kx,i) = dihelist (4,ix,i)
            diheperio(kx,i)   = diheperio(ix,i)
            dihekind (kx,i)   = dihekind (ix,i)
            dihe_pbc (1,kx,i) = dihe_pbc (1,ix,i)
            dihe_pbc (2,kx,i) = dihe_pbc (2,ix,i)
            dihe_pbc (3,kx,i) = dihe_pbc (3,ix,i)
            dihe_singleB(kx,i) = dihe_singleB(ix,i)
            dihephase(kx,i)   = dihephase(ix,i)
            diheforce(kx,i)   = diheforce(ix,i)
            j = j + 1
            k = k + 1
          end if

          ix = ix - 1

        end do
        ndihedral(i) = ndihedral(i) + dihed_add(i) - dihed_exit(i)

      end if

    end do

    return

  end subroutine update_incoming_enefunc_dihe_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_enefunc_rb_dihe_fep
  !> @brief        update Ryckaert-Bellemans DIHEDRAL term for each cell in
  !!               potential energy function for FEP
  !! @authors      HO
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_enefunc_rb_dihe_fep(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: icel1, icel2, icel3, icel4, i1, i2, i3, i4
    integer                  :: i, ix, icel_local

    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: buf_real(:,:,:), dihec(:,:,:)
    integer,         pointer :: ncel_local, nboundary
    integer(int2),   pointer :: cell_g2l(:), cell_g2b(:)
    integer(int2),   pointer :: cell_pair(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: dihed_add(:), dihed_exit(:)
    integer,         pointer :: dihed_exit_index(:,:)
    integer,         pointer :: ndihe(:), dihelist(:,:,:)
    integer,         pointer :: dihe_pbc(:,:,:)
    integer,         pointer :: buf_int(:,:,:)

    ! FEP
    integer,         pointer :: dihe_singleB(:,:)

    ncel_local       => domain%num_cell_local
    nboundary        => domain%num_cell_boundary
    cell_g2l         => domain%cell_g2l
    cell_g2b         => domain%cell_g2b
    cell_pair        => domain%cell_pair
    id_g2l           => domain%id_g2l
    coord            => domain%coord

    ndihe            => enefunc%num_rb_dihedral
    dihelist         => enefunc%rb_dihe_list
    dihec            => enefunc%rb_dihe_c
    dihed_add        => enefunc%rb_dihed_add
    dihed_exit       => enefunc%rb_dihed_exit
    dihed_exit_index => enefunc%rb_dihed_exit_index
    dihe_pbc         => enefunc%rb_dihe_pbc
    buf_int          => enefunc%buf_rb_dihed_integer
    buf_real         => enefunc%buf_rb_dihed_real

    ! FEP
    dihe_singleB     => enefunc%rb_dihe_singleB

    dihed_add (1:ncel_local+nboundary) = 0
    dihed_exit(1:ncel_local) = 0


    do i = 1, ncel_local
      do ix = 1, ndihe(i)

        icel1 = id_g2l(1,dihelist(1,ix,i))
        i1    = id_g2l(2,dihelist(1,ix,i))
        icel2 = id_g2l(1,dihelist(2,ix,i))
        i2    = id_g2l(2,dihelist(2,ix,i))
        icel3 = id_g2l(1,dihelist(3,ix,i))
        i3    = id_g2l(2,dihelist(3,ix,i))
        icel4 = id_g2l(1,dihelist(4,ix,i))
        i4    = id_g2l(2,dihelist(4,ix,i))

        icel_local = cell_pair(icel1,icel4)

        if (icel_local /= i) then

          dihed_exit(i) = dihed_exit(i) + 1
          dihed_exit_index(dihed_exit(i),i) = ix
          dihed_add(icel_local) = dihed_add(icel_local) + 1
          buf_real(1,dihed_add(icel_local),icel_local) = dihec(1,ix,i)
          buf_real(2,dihed_add(icel_local),icel_local) = dihec(2,ix,i)
          buf_real(3,dihed_add(icel_local),icel_local) = dihec(3,ix,i)
          buf_real(4,dihed_add(icel_local),icel_local) = dihec(4,ix,i)
          buf_real(5,dihed_add(icel_local),icel_local) = dihec(5,ix,i)
          buf_real(6,dihed_add(icel_local),icel_local) = dihec(6,ix,i)
          buf_int (1,dihed_add(icel_local),icel_local) = dihelist(1,ix,i)
          buf_int (2,dihed_add(icel_local),icel_local) = dihelist(2,ix,i)
          buf_int (3,dihed_add(icel_local),icel_local) = dihelist(3,ix,i)
          buf_int (4,dihed_add(icel_local),icel_local) = dihelist(4,ix,i)
          buf_int (5,dihed_add(icel_local),icel_local) = dihe_pbc(1,ix,i)
          buf_int (6,dihed_add(icel_local),icel_local) = dihe_pbc(2,ix,i)
          buf_int (7,dihed_add(icel_local),icel_local) = dihe_pbc(3,ix,i)
          buf_int (8,dihed_add(icel_local),icel_local) = dihe_singleB(ix,i)

        end if

      end  do
    end do

    return

  end subroutine update_outgoing_enefunc_rb_dihe_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_enefunc_rb_dihe_fep
  !> @brief        update Ryckaert-Bellemans DIHEDRAL term for each cell in
  !!               potential energy function for FEP
  !! @authors      HO
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_incoming_enefunc_rb_dihe_fep(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    integer                  :: i, j, k, ix, kx, ndihe
    logical                  :: insert

    real(wp),        pointer :: dihec(:,:,:), buf_real(:,:,:)
    integer,         pointer :: ncel_local, nboundary
    integer(int2),   pointer :: cell_g2l(:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: dihed_add(:), dihed_exit(:)
    integer,         pointer :: dihed_exit_index(:,:), buf_int(:,:,:)
    integer,         pointer :: ndihedral(:), dihelist(:,:,:)
    integer,         pointer :: dihe_pbc(:,:,:)

    ! FEP
    integer,         pointer :: dihe_singleB(:,:)

    ncel_local       => domain%num_cell_local
    nboundary        => domain%num_cell_boundary
    cell_g2l         => domain%cell_g2l
    id_g2l           => domain%id_g2l

    ndihedral        => enefunc%num_rb_dihedral
    dihelist         => enefunc%rb_dihe_list
    dihec            => enefunc%rb_dihe_c
    dihed_add        => enefunc%rb_dihed_add
    dihed_exit       => enefunc%rb_dihed_exit
    dihed_exit_index => enefunc%rb_dihed_exit_index
    dihe_pbc         => enefunc%rb_dihe_pbc
    buf_int          => enefunc%buf_rb_dihed_integer
    buf_real         => enefunc%buf_rb_dihed_real

    ! FEP
    dihe_singleB     => enefunc%rb_dihe_singleB

    ndihe = enefunc%num_rb_dihe_all

    do i = 1, ncel_local

#ifdef DEBUG
      if (dihed_add(i)+ndihedral(i)-dihed_exit(i) > MaxDihe) &
        call error_msg( &
            'Debug: update_incoming_enefunc_rb_dihed> dihed is exceed MaxDihe')
#endif

      if (dihed_add(i) > DiheMove) &
        call error_msg( &
   'Debug: update_incoming_enefunc_rb_dihed> dihed migration is exceed maximum')

      if (dihed_add(i) >= dihed_exit(i)) then

        do k = 1, dihed_exit(i)
          dihelist (1,dihed_exit_index(k,i),i) = buf_int (1,k,i)
          dihelist (2,dihed_exit_index(k,i),i) = buf_int (2,k,i)
          dihelist (3,dihed_exit_index(k,i),i) = buf_int (3,k,i)
          dihelist (4,dihed_exit_index(k,i),i) = buf_int (4,k,i)
          dihe_pbc (1,dihed_exit_index(k,i),i) = buf_int (5,k,i)
          dihe_pbc (2,dihed_exit_index(k,i),i) = buf_int (6,k,i)
          dihe_pbc (3,dihed_exit_index(k,i),i) = buf_int (7,k,i)
          dihe_singleB(dihed_exit_index(k,i),i) = buf_int (8,k,i)
          dihec    (1,dihed_exit_index(k,i),i) = buf_real(1,k,i)
          dihec    (2,dihed_exit_index(k,i),i) = buf_real(2,k,i)
          dihec    (3,dihed_exit_index(k,i),i) = buf_real(3,k,i)
          dihec    (4,dihed_exit_index(k,i),i) = buf_real(4,k,i)
          dihec    (5,dihed_exit_index(k,i),i) = buf_real(5,k,i)
          dihec    (6,dihed_exit_index(k,i),i) = buf_real(6,k,i)
        end do

        do k = dihed_exit(i)+1, dihed_add(i)
          ix = k + ndihedral(i) - dihed_exit(i)
          dihelist (1,ix,i) = buf_int (1,k,i)
          dihelist (2,ix,i) = buf_int (2,k,i)
          dihelist (3,ix,i) = buf_int (3,k,i)
          dihelist (4,ix,i) = buf_int (4,k,i)
          dihe_pbc (1,ix,i) = buf_int (5,k,i)
          dihe_pbc (2,ix,i) = buf_int (6,k,i)
          dihe_pbc (3,ix,i) = buf_int (7,k,i)
          dihe_singleB(ix,i) = buf_int (8,k,i)
          dihec    (1,ix,i) = buf_real(1,k,i)
          dihec    (2,ix,i) = buf_real(2,k,i)
          dihec    (3,ix,i) = buf_real(3,k,i)
          dihec    (4,ix,i) = buf_real(4,k,i)
          dihec    (5,ix,i) = buf_real(5,k,i)
          dihec    (6,ix,i) = buf_real(6,k,i)
        end do

        ndihedral(i) = ndihedral(i) + dihed_add(i) - dihed_exit(i)

      else
        do k = 1, dihed_add(i)
          dihelist (1,dihed_exit_index(k,i),i) = buf_int (1,k,i)
          dihelist (2,dihed_exit_index(k,i),i) = buf_int (2,k,i)
          dihelist (3,dihed_exit_index(k,i),i) = buf_int (3,k,i)
          dihelist (4,dihed_exit_index(k,i),i) = buf_int (4,k,i)
          dihe_pbc (1,dihed_exit_index(k,i),i) = buf_int (5,k,i)
          dihe_pbc (2,dihed_exit_index(k,i),i) = buf_int (6,k,i)
          dihe_pbc (3,dihed_exit_index(k,i),i) = buf_int (7,k,i)
          dihe_singleB(dihed_exit_index(k,i),i) = buf_int (8,k,i)
          dihec    (1,dihed_exit_index(k,i),i) = buf_real(1,k,i)
          dihec    (2,dihed_exit_index(k,i),i) = buf_real(2,k,i)
          dihec    (3,dihed_exit_index(k,i),i) = buf_real(3,k,i)
          dihec    (4,dihed_exit_index(k,i),i) = buf_real(4,k,i)
          dihec    (5,dihed_exit_index(k,i),i) = buf_real(5,k,i)
          dihec    (6,dihed_exit_index(k,i),i) = buf_real(6,k,i)
        end do

        j  = 0
        ix = ndihedral(i)
        k  = dihed_add(i) + 1

        do while (j < (dihed_exit(i)-dihed_add(i)))

          insert = .true.
          do kx = k, dihed_exit(i)
            if (ix == dihed_exit_index(kx,i)) then
              insert = .false.
              j = j + 1
              exit
            end if
          end do

          if (insert) then
            kx = dihed_exit_index(k,i)
            dihelist (1,kx,i) = dihelist (1,ix,i)
            dihelist (2,kx,i) = dihelist (2,ix,i)
            dihelist (3,kx,i) = dihelist (3,ix,i)
            dihelist (4,kx,i) = dihelist (4,ix,i)
            dihe_pbc (1,kx,i) = dihe_pbc (1,ix,i)
            dihe_pbc (2,kx,i) = dihe_pbc (2,ix,i)
            dihe_pbc (3,kx,i) = dihe_pbc (3,ix,i)
            dihe_singleB(kx,i) = dihe_singleB(ix,i)
            dihec    (1,kx,i) = dihec    (1,ix,i)
            dihec    (2,kx,i) = dihec    (2,ix,i)
            dihec    (3,kx,i) = dihec    (3,ix,i)
            dihec    (4,kx,i) = dihec    (4,ix,i)
            dihec    (5,kx,i) = dihec    (5,ix,i)
            dihec    (6,kx,i) = dihec    (6,ix,i)
            j = j + 1
            k = k + 1
          end if

          ix = ix - 1

        end do
        ndihedral(i) = ndihedral(i) + dihed_add(i) - dihed_exit(i)

      end if

    end do

    return

  end subroutine update_incoming_enefunc_rb_dihe_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_enefunc_impr_fep
  !> @brief        update IMPROPER DIHEDRAL term for each cell in potential
  !                energy function for FEP
  !! @authors      HO
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_enefunc_impr_fep(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: icel1, icel2, icel3, icel4, i1, i2, i3, i4
    integer                  :: i, ix, icel_local

    real(wp),        pointer :: fc(:,:), phase(:,:)
    real(wip),       pointer :: coord(:,:,:)
    real(wp),        pointer :: buf_real(:,:,:)
    integer,         pointer :: ncel_local, nboundary
    integer(int2),   pointer :: cell_g2l(:), cell_g2b(:)
    integer(int2),   pointer :: cell_pair(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: impr_add(:), impr_exit(:), impr_exit_index(:,:)
    integer,         pointer :: nperiod(:,:)
    integer,         pointer :: buf_int(:,:,:), nimproper(:), imprlist(:,:,:)
    integer,         pointer :: impr_pbc(:,:,:)

    ! FEP
    integer,         pointer :: impr_singleB(:,:)

    ncel_local      => domain%num_cell_local
    nboundary       => domain%num_cell_boundary
    cell_g2l        => domain%cell_g2l
    cell_g2b        => domain%cell_g2b
    cell_pair       => domain%cell_pair
    id_g2l          => domain%id_g2l
    coord           => domain%coord

    nimproper       => enefunc%num_improper
    imprlist        => enefunc%impr_list
    impr_add        => enefunc%impr_add
    impr_exit       => enefunc%impr_exit
    impr_exit_index => enefunc%impr_exit_index
    fc              => enefunc%impr_force_const
    phase           => enefunc%impr_phase
    nperiod         => enefunc%impr_periodicity
    buf_int         => enefunc%buf_impr_integer
    buf_real        => enefunc%buf_impr_real
    impr_pbc        => enefunc%impr_pbc

    ! FEP
    impr_singleB    => enefunc%impr_singleB

    impr_exit(1:ncel_local) = 0
    impr_add (1:ncel_local+nboundary) = 0

    do i = 1, ncel_local
      do ix = 1, nimproper(i)

        icel1 = id_g2l(1,imprlist(1,ix,i))
        i1    = id_g2l(2,imprlist(1,ix,i))
        icel2 = id_g2l(1,imprlist(2,ix,i))
        i2    = id_g2l(2,imprlist(2,ix,i))
        icel3 = id_g2l(1,imprlist(3,ix,i))
        i3    = id_g2l(2,imprlist(3,ix,i))
        icel4 = id_g2l(1,imprlist(4,ix,i))
        i4    = id_g2l(2,imprlist(4,ix,i))

        icel_local = cell_pair(icel1,icel4)

        if (icel_local /= i) then

          impr_exit(i) = impr_exit(i) + 1
          impr_exit_index(impr_exit(i),i) = ix
          impr_add(icel_local) = impr_add(icel_local) + 1
          buf_real(1,impr_add(icel_local),icel_local) = phase(ix,i)
          buf_real(2,impr_add(icel_local),icel_local) = fc(ix,i)
          buf_int (1,impr_add(icel_local),icel_local) = imprlist(1,ix,i)
          buf_int (2,impr_add(icel_local),icel_local) = imprlist(2,ix,i)
          buf_int (3,impr_add(icel_local),icel_local) = imprlist(3,ix,i)
          buf_int (4,impr_add(icel_local),icel_local) = imprlist(4,ix,i)
          buf_int (5,impr_add(icel_local),icel_local) = nperiod(ix,i)
          buf_int (6,impr_add(icel_local),icel_local) = impr_pbc(1,ix,i)
          buf_int (7,impr_add(icel_local),icel_local) = impr_pbc(2,ix,i)
          buf_int (8,impr_add(icel_local),icel_local) = impr_pbc(3,ix,i)
          buf_int (9,impr_add(icel_local),icel_local) = impr_singleB(ix,i)

        end if

      end  do
    end do

    return

  end subroutine update_outgoing_enefunc_impr_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_enefunc_impr_fep
  !> @brief        update IMPROPER DIHEDRAL term for each cell in potential
  !!               energy function for FEP
  !! @authors      HO
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_incoming_enefunc_impr_fep(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    integer                  :: i, j, k, ix, kx
    logical                  :: insert

    real(wp),        pointer :: buf_real(:,:,:)
    real(wp),        pointer :: impr_force(:,:), impr_phase(:,:)
    integer,         pointer :: ncel_local, nboundary
    integer(2),      pointer :: cell_g2l(:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: impr_add(:), impr_exit(:), impr_exit_index(:,:)
    integer,         pointer :: impr_nperiod(:,:)
    integer,         pointer :: buf_int(:,:,:), nimproper(:), imprlist(:,:,:)
    integer,         pointer :: impr_pbc(:,:,:)

    ! FEP
    integer,         pointer :: impr_singleB(:,:)

    ncel_local       => domain%num_cell_local
    nboundary        => domain%num_cell_boundary
    cell_g2l         => domain%cell_g2l
    id_g2l           => domain%id_g2l

    nimproper        => enefunc%num_improper
    imprlist         => enefunc%impr_list
    impr_add         => enefunc%impr_add
    impr_exit        => enefunc%impr_exit
    impr_exit_index  => enefunc%impr_exit_index
    impr_force       => enefunc%impr_force_const
    impr_phase       => enefunc%impr_phase
    impr_nperiod     => enefunc%impr_periodicity
    buf_int          => enefunc%buf_impr_integer
    buf_real         => enefunc%buf_impr_real
    impr_pbc         => enefunc%impr_pbc

    ! FEP
    impr_singleB     => enefunc%impr_singleB

    do i = 1, ncel_local
#ifdef DEBUG
      if (impr_add(i)+nimproper(i)-impr_exit(i) > MaxImpr) &
        call error_msg( &
            'Debug: update_incoming_enefunc_impr> improper is exceed MaxImpr')
#endif

      if (impr_add(i) > ImprMove) &
        call error_msg( &
   'Debug: update_incoming_enefunc_impr> improper migration is exceed maximum')

      if (impr_add(i) >= impr_exit(i)) then

        do k = 1, impr_exit(i)
          imprlist(1,impr_exit_index(k,i),i)   = buf_int (1,k,i)
          imprlist(2,impr_exit_index(k,i),i)   = buf_int (2,k,i)
          imprlist(3,impr_exit_index(k,i),i)   = buf_int (3,k,i)
          imprlist(4,impr_exit_index(k,i),i)   = buf_int (4,k,i)
          impr_nperiod(impr_exit_index(k,i),i) = buf_int (5,k,i)
          impr_pbc(1,impr_exit_index(k,i),i)   = buf_int (6,k,i)
          impr_pbc(2,impr_exit_index(k,i),i)   = buf_int (7,k,i)
          impr_pbc(3,impr_exit_index(k,i),i)   = buf_int (8,k,i)
          impr_singleB(impr_exit_index(k,i),i) = buf_int (9,k,i)
          impr_phase(impr_exit_index(k,i),i)   = buf_real(1,k,i)
          impr_force(impr_exit_index(k,i),i)   = buf_real(2,k,i)
        end do

        do k = impr_exit(i)+1, impr_add(i)
          ix = k + nimproper(i) - impr_exit(i)
          imprlist(1,ix,i)   = buf_int (1,k,i)
          imprlist(2,ix,i)   = buf_int (2,k,i)
          imprlist(3,ix,i)   = buf_int (3,k,i)
          imprlist(4,ix,i)   = buf_int (4,k,i)
          impr_nperiod(ix,i) = buf_int (5,k,i)
          impr_pbc(1,ix,i)   = buf_int (6,k,i)
          impr_pbc(2,ix,i)   = buf_int (7,k,i)
          impr_pbc(3,ix,i)   = buf_int (8,k,i)
          impr_singleB(ix,i) = buf_int (9,k,i)
          impr_phase(ix,i)   = buf_real(1,k,i)
          impr_force(ix,i)   = buf_real(2,k,i)
        end do

        nimproper(i) = nimproper(i) + impr_add(i) - impr_exit(i)

      else
        do k = 1, impr_add(i)
          imprlist(1,impr_exit_index(k,i),i)   = buf_int (1,k,i)
          imprlist(2,impr_exit_index(k,i),i)   = buf_int (2,k,i)
          imprlist(3,impr_exit_index(k,i),i)   = buf_int (3,k,i)
          imprlist(4,impr_exit_index(k,i),i)   = buf_int (4,k,i)
          impr_nperiod(impr_exit_index(k,i),i) = buf_int (5,k,i)
          impr_pbc(1,impr_exit_index(k,i),i)   = buf_int (6,k,i)
          impr_pbc(2,impr_exit_index(k,i),i)   = buf_int (7,k,i)
          impr_pbc(3,impr_exit_index(k,i),i)   = buf_int (8,k,i)
          impr_singleB(impr_exit_index(k,i),i) = buf_int (9,k,i)
          impr_phase(impr_exit_index(k,i),i)   = buf_real(1,k,i)
          impr_force(impr_exit_index(k,i),i)   = buf_real(2,k,i)
        end do

        j  = 0
        ix = nimproper(i)
        k  = impr_add(i) + 1

        do while (j < (impr_exit(i)-impr_add(i)))

          insert = .true.
          do kx = k, impr_exit(i)
            if (ix == impr_exit_index(kx,i)) then
              insert = .false.
              j = j + 1
              exit
            end if
          end do

          if (insert) then
            kx = impr_exit_index(k,i)
            imprlist(1,kx,i)   = imprlist(1,ix,i)
            imprlist(2,kx,i)   = imprlist(2,ix,i)
            imprlist(3,kx,i)   = imprlist(3,ix,i)
            imprlist(4,kx,i)   = imprlist(4,ix,i)
            impr_nperiod(kx,i) = impr_nperiod(ix,i)
            impr_pbc(1,kx,i)   = impr_pbc(1,ix,i)
            impr_pbc(2,kx,i)   = impr_pbc(2,ix,i)
            impr_pbc(3,kx,i)   = impr_pbc(3,ix,i)
            impr_singleB(kx,i) = impr_singleB(ix,i)
            impr_phase(kx,i)   = impr_phase(ix,i)
            impr_force(kx,i)   = impr_force(ix,i)
            j = j + 1
            k = k + 1
          end if
          ix = ix - 1
        end do

        nimproper(i) = nimproper(i) + impr_add(i) - impr_exit(i)

      end if
    end do

    return

  end subroutine update_incoming_enefunc_impr_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_enefunc_cmap_fep
  !> @brief        update CMAP term for each cell in potential
  !                energy function for FEP
  !! @authors      HO
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_enefunc_cmap_fep(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: icel1, icel2
    integer                  :: i, ix, k, icel_local

    real(wip),       pointer :: coord(:,:,:)
    integer,         pointer :: ncel_local, nboundary
    integer(int2),   pointer :: cell_g2l(:), cell_g2b(:)
    integer(int2),   pointer :: cell_pair(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: cmap_exit(:), cmap_add(:), cmap_exit_index(:,:)
    integer,         pointer :: ncmap(:), cmaplist(:,:,:), cmaptype(:,:)
    integer,         pointer :: buf_int(:,:,:)
    integer,         pointer :: cmap_pbc(:,:,:)

    ! FEP
    integer,         pointer :: cmap_singleB(:,:)

    ncel_local      => domain%num_cell_local
    nboundary       => domain%num_cell_boundary
    cell_g2l        => domain%cell_g2l
    cell_g2b        => domain%cell_g2b
    cell_pair       => domain%cell_pair
    id_g2l          => domain%id_g2l
    coord           => domain%coord

    ncmap           => enefunc%num_cmap
    cmaplist        => enefunc%cmap_list
    cmaptype        => enefunc%cmap_type
    cmap_add        => enefunc%cmap_add
    cmap_exit       => enefunc%cmap_exit
    cmap_pbc        => enefunc%cmap_pbc
    cmap_exit_index => enefunc%cmap_exit_index
    buf_int         => enefunc%buf_cmap_integer

    ! FEP
    cmap_singleB    => enefunc%cmap_singleB

    cmap_exit(1:ncel_local) = 0
    cmap_add (1:ncel_local+nboundary) = 0

    do i = 1, ncel_local
      do ix = 1, ncmap(i)

        icel1 = id_g2l(1,cmaplist(1,ix,i))
        icel2 = id_g2l(1,cmaplist(8,ix,i))
        icel_local = cell_pair(icel1,icel2)

        if (icel_local /= i) then

          cmap_exit(i) = cmap_exit(i) + 1
          cmap_exit_index(cmap_exit(i),i) = ix
          cmap_add(icel_local) = cmap_add(icel_local) + 1
          do k = 1, 8
            buf_int(k,cmap_add(icel_local),icel_local) =cmaplist(k,ix,i)
          end do
          buf_int(9,cmap_add(icel_local),icel_local) = cmaptype(ix,i)
          buf_int(10,cmap_add(icel_local),icel_local) = cmap_pbc(1,ix,i)
          buf_int(11,cmap_add(icel_local),icel_local) = cmap_pbc(2,ix,i)
          buf_int(12,cmap_add(icel_local),icel_local) = cmap_pbc(3,ix,i)
          buf_int(13,cmap_add(icel_local),icel_local) = cmap_pbc(4,ix,i)
          buf_int(14,cmap_add(icel_local),icel_local) = cmap_pbc(5,ix,i)
          buf_int(15,cmap_add(icel_local),icel_local) = cmap_pbc(6,ix,i)
          buf_int(16,cmap_add(icel_local),icel_local) = cmap_singleB(ix,i)

        end if

      end  do
    end do

    return

  end subroutine update_outgoing_enefunc_cmap_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_enefunc_cmap_fep
  !> @brief        update CMAP term for each cell in potential energy function
  !                for FEP
  !! @authors      HO
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_incoming_enefunc_cmap_fep(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    integer                  :: i, j, k, ix, kx, ia
    logical                  :: insert

    integer,         pointer :: ncel_local, nboundary
    integer(int2),   pointer :: cell_g2l(:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: cmap_add(:), cmap_exit(:), cmap_exit_index(:,:)
    integer,         pointer :: buf_int(:,:,:), ncmap(:), cmap_list(:,:,:)
    integer,         pointer :: cmap_type(:,:)
    integer,         pointer :: cmap_pbc(:,:,:)

    ! FEP
    integer,         pointer :: cmap_singleB(:,:)

    ncel_local      => domain%num_cell_local
    nboundary       => domain%num_cell_boundary
    cell_g2l        => domain%cell_g2l
    id_g2l          => domain%id_g2l

    ncmap           => enefunc%num_cmap
    cmap_list       => enefunc%cmap_list
    cmap_type       => enefunc%cmap_type
    cmap_add        => enefunc%cmap_add
    cmap_exit       => enefunc%cmap_exit
    cmap_exit_index => enefunc%cmap_exit_index
    cmap_pbc        => enefunc%cmap_pbc
    buf_int         => enefunc%buf_cmap_integer

    ! FEP
    cmap_singleB    => enefunc%cmap_singleB

    do i = 1, ncel_local

#ifdef DEBUG
      if (cmap_add(i)+ncmap(i)-cmap_exit(i) > MaxCmap) &
        call error_msg( &
          'Debug: update_incoming_enefunc_cmap> cmap is exceed MaxCmap')
#endif

      if (cmap_add(i) > CmapMove) &
        call error_msg( &
      'Debug: update_incoming_enefunc_cmap> cmap migration is exceed maximum')

      if (cmap_add(i) >= cmap_exit(i)) then

        do k = 1, cmap_exit(i)
          do ia = 1, 8
            cmap_list(ia,cmap_exit_index(k,i),i) = buf_int(ia,k,i)
          end do
          cmap_type(cmap_exit_index(k,i),i) = buf_int(9,k,i)
          cmap_pbc (1,cmap_exit_index(k,i),i) = buf_int(10,k,i)
          cmap_pbc (2,cmap_exit_index(k,i),i) = buf_int(11,k,i)
          cmap_pbc (3,cmap_exit_index(k,i),i) = buf_int(12,k,i)
          cmap_pbc (4,cmap_exit_index(k,i),i) = buf_int(13,k,i)
          cmap_pbc (5,cmap_exit_index(k,i),i) = buf_int(14,k,i)
          cmap_pbc (6,cmap_exit_index(k,i),i) = buf_int(15,k,i)
          cmap_singleB(cmap_exit_index(k,i),i) = buf_int(16,k,i)
        end do

        do k = cmap_exit(i)+1, cmap_add(i)

          ix = k + ncmap(i) - cmap_exit(i)
          do ia = 1, 8
            cmap_list(ia,ix,i) = buf_int(ia,k,i)
          end do
          cmap_type(ix,i) = buf_int(9,k,i)
          cmap_pbc (1,ix,i) = buf_int(10,k,i)
          cmap_pbc (2,ix,i) = buf_int(11,k,i)
          cmap_pbc (3,ix,i) = buf_int(12,k,i)
          cmap_pbc (4,ix,i) = buf_int(13,k,i)
          cmap_pbc (5,ix,i) = buf_int(14,k,i)
          cmap_pbc (6,ix,i) = buf_int(15,k,i)
          cmap_singleB(ix,i) = buf_int(16,k,i)

        end do

        ncmap(i) = ncmap(i) + cmap_add(i) - cmap_exit(i)

      else

        do k = 1, cmap_add(i)
          do ia = 1, 8
            cmap_list(ia,cmap_exit_index(k,i),i) = buf_int(ia,k,i)
          end do
          cmap_type(cmap_exit_index(k,i),i)   = buf_int(9,k,i)
          cmap_pbc (1,cmap_exit_index(k,i),i) = buf_int(10,k,i)
          cmap_pbc (2,cmap_exit_index(k,i),i) = buf_int(11,k,i)
          cmap_pbc (3,cmap_exit_index(k,i),i) = buf_int(12,k,i)
          cmap_pbc (4,cmap_exit_index(k,i),i) = buf_int(13,k,i)
          cmap_pbc (5,cmap_exit_index(k,i),i) = buf_int(14,k,i)
          cmap_pbc (6,cmap_exit_index(k,i),i) = buf_int(15,k,i)
          cmap_singleB(cmap_exit_index(k,i),i) = buf_int(16,k,i)
        end do

        j  = 0
        ix = ncmap(i)
        k  = cmap_add(i) + 1

        do while (j < (cmap_exit(i)-cmap_add(i)))

          insert = .true.
          do kx = k, cmap_exit(i)
            if (ix == cmap_exit_index(kx,i)) then
              insert = .false.
              j = j + 1
              exit
            end if
          end do

          if (insert) then

            kx = cmap_exit_index(k,i)
            do ia = 1, 8
              cmap_list(ia,kx,i) = cmap_list(ia,ix,i)
            end do

            cmap_type(kx,i) = cmap_type(ix,i)
            cmap_pbc (1,kx,i) = cmap_pbc(1,ix,i)
            cmap_pbc (2,kx,i) = cmap_pbc(2,ix,i)
            cmap_pbc (3,kx,i) = cmap_pbc(3,ix,i)
            cmap_pbc (4,kx,i) = cmap_pbc(4,ix,i)
            cmap_pbc (5,kx,i) = cmap_pbc(5,ix,i)
            cmap_pbc (6,kx,i) = cmap_pbc(6,ix,i)
            cmap_singleB(kx,i) = cmap_singleB(ix,i)
            j = j + 1
            k = k + 1

          end if

          ix = ix - 1

        end do

        ncmap(i) = ncmap(i) + cmap_add(i) - cmap_exit(i)

      end if

    end do

    return

  end subroutine update_incoming_enefunc_cmap_fep

end module sp_migration_mod
