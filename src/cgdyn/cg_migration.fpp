!--------1---------2---------3---------4---------5---------6---------7---------8
! 
!  Module   cg_migration_mod
!> @brief   Manages the migration of atoms and energy functions within each domain.
!! @authors Jaewoon Jung (JJ) , Chigusa Kobayshi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_migration_mod

  use cg_boundary_str_mod
  use cg_enefunc_str_mod
  use cg_domain_str_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef OMP
    use omp_lib, only: omp_get_thread_num
#endif

  implicit none
  private

  ! subroutines
  public :: update_outgoing_charge
  public :: update_outgoing_nocharge
  public :: update_incoming_ptl
  public :: update_outgoing_charge_martini
  public :: update_outgoing_nocharge_martini
  public :: update_incoming_ptl_martini
  public :: update_outgoing_enefunc_bond
  public :: update_incoming_enefunc_bond
  public :: update_outgoing_enefunc_bondsq
  public :: update_outgoing_enefunc_angl
  public :: update_outgoing_enefunc_anglflex
  public :: update_outgoing_enefunc_angllocal
  public :: update_incoming_enefunc_angl
  public :: update_outgoing_enefunc_dihe
  public :: update_outgoing_enefunc_diheflex
  public :: update_outgoing_enefunc_dihelocal
  public :: update_incoming_enefunc_dihe
  public :: update_outgoing_enefunc_stack
  public :: update_incoming_enefunc_stack
  public :: update_outgoing_enefunc_pwmcos
  public :: update_incoming_enefunc_pwmcos
  public :: update_outgoing_enefunc_pwmcosns
  public :: update_incoming_enefunc_pwmcosns
  public :: update_enefunc_contact
  public :: update_outgoing_enefunc_restraint
  public :: update_incoming_enefunc_restraint

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_charge
  !> @brief        Checks for charged particles moving to other cells.
  !! @authors      JJ
  !! @param[in]    boundary    : boundary condition information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_charge(boundary, domain)

    ! formal arguments
    type(s_boundary),    target, intent(in)    :: boundary
    type(s_domain),      target, intent(inout) :: domain

    ! local variable
    real(wip)                    :: x_shift, y_shift, z_shift
    real(wip)                    :: move(3)
    integer                      :: i, icx, icy, icz, icel, ncel
    integer                      :: ilist, ig, n_stay, num
    integer                      :: icel_local, icel_bd, ip
    integer                      :: ki, ko

    real(wip),           pointer :: bsize_x, bsize_y, bsize_z
    real(wip),           pointer :: csize_x, csize_y, csize_z
    real(wip),           pointer :: coord(:,:), velocity(:,:)
    real(wp),            pointer :: charge(:)
    real(wip),           pointer :: mass(:)
    real(wip),           pointer :: buf_stay_real(:)
    real(wip),           pointer :: buf_move_real(:)
    real(wip),           pointer :: buf_comm_real(:)
    integer,             pointer :: ncel_x, ncel_y, ncel_z, ncel_local, ncel_bd
    integer,             pointer :: natom(:), ncharge(:)
    integer,             pointer :: cell_g2l(:), cell_g2b(:)
    integer,             pointer :: cell_b2g(:)
    integer,             pointer :: atmcls(:), chain_id(:), atom_2_cell(:)
    integer,             pointer :: id_l2g(:), id_g2l(:)
    integer,             pointer :: atom_type(:)
    integer,             pointer :: charge_move(:), charge_comm(:)
    integer,             pointer :: charge_stay(:), charge_comm1(:)
    integer,             pointer :: cell_rank(:)
    integer,             pointer :: buf_stay_int(:)
    integer,             pointer :: buf_move_int(:)
    integer,             pointer :: buf_comm_int(:)

    bsize_x       => boundary%box_size_x
    bsize_y       => boundary%box_size_y
    bsize_z       => boundary%box_size_z
    ncel_x        => boundary%num_cells_x
    ncel_y        => boundary%num_cells_y
    ncel_z        => boundary%num_cells_z
    csize_x       => boundary%cell_size_x
    csize_y       => boundary%cell_size_y
    csize_z       => boundary%cell_size_z

    ncel_local    => domain%num_cell_local
    natom         => domain%num_atom
    ncharge       => domain%num_charge
    ncel_bd       => domain%num_cell_boundary
    cell_g2l      => domain%cell_g2l
    cell_g2b      => domain%cell_g2b
    cell_b2g      => domain%cell_b2g
    coord         => domain%coord
    velocity      => domain%velocity
    charge        => domain%charge
    mass          => domain%mass
    atmcls        => domain%atom_cls_no
    atom_2_cell   => domain%atom_2_cell
    id_l2g        => domain%id_l2g
    id_g2l        => domain%id_g2l
    chain_id      => domain%mol_chain_id
    atom_type     => domain%NA_base_type
    cell_rank     => domain%domain_cell_rank
    charge_move   => domain%type1_move
    charge_comm   => domain%type1_comm
    charge_comm1  => domain%type1_comm_move
    charge_stay   => domain%type1_stay
    buf_stay_real => domain%buf_var0_stay_real
    buf_move_real => domain%buf_var0_move_real
    buf_comm_real => domain%buf_var0_comm_real
    buf_stay_int  => domain%buf_var0_stay_int 
    buf_move_int  => domain%buf_var0_move_int 
    buf_comm_int  => domain%buf_var0_comm_int 

    ! initializaiton
    !
    ncel = ncel_local + ncel_bd
    num  = domain%num_comm_proc
    charge_move(1:ncel) = 0
    charge_comm(1:num) = 0
    charge_comm1(1:num) = 0
    charge_stay(1:ncel_local) = 0

    ! initialized the global index
    !
    !$omp parallel do private(ig,i)
    do i = 1, domain%num_atom_domain+domain%num_atom_boundary
      ig = id_l2g(i)
      id_g2l(ig) = 0
    end do
    !$omp end parallel do

    ! Check outgoing particles
    !
    n_stay = 0
    ki     = 0
    ko     = 0

    do ilist = 1, domain%num_atom_domain

      if (abs(charge(ilist)) > EPS) then

        i = atom_2_cell(ilist)
        x_shift = coord(ilist,1) - boundary%origin_x
        y_shift = coord(ilist,2) - boundary%origin_y
        z_shift = coord(ilist,3) - boundary%origin_z

        !coordinate shifted to the first quadrant and set into the boundary box
        if (boundary%type == BoundaryTypePBC) then
          move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
          move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
          move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)
          x_shift = x_shift + move(1)
          y_shift = y_shift + move(2)
          z_shift = z_shift + move(3)
        else
          x_shift = x_shift + bsize_x*0.5_wip
          y_shift = y_shift + bsize_y*0.5_wip
          z_shift = z_shift + bsize_z*0.5_wip
        end if

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

          if (icel_local /= 0) then

            ki = ki + 1
            charge_move(icel_local) = charge_move(icel_local) + 1
            buf_move_real(8*ki-7) = coord    (ilist,1)
            buf_move_real(8*ki-6) = coord    (ilist,2)
            buf_move_real(8*ki-5) = coord    (ilist,3)
            buf_move_real(8*ki-4) = velocity (ilist,1)
            buf_move_real(8*ki-3) = velocity (ilist,2)
            buf_move_real(8*ki-2) = velocity (ilist,3)
            buf_move_real(8*ki-1) = mass     (ilist  )
            buf_move_real(8*ki  ) = charge   (ilist  )
            buf_move_int (5*ki-4) = icel_local
            buf_move_int (5*ki-3) = atmcls   (ilist  )
            buf_move_int (5*ki-2) = id_l2g   (ilist  )
            buf_move_int (5*ki-1) = chain_id (ilist  ) 
            buf_move_int (5*ki  ) = atom_type(ilist  ) 

          else if (icel_bd /= 0) then

            ig = cell_b2g(icel_bd)
            icel_bd = icel_bd + ncel_local
            ip = cell_rank(icel_bd)
            charge_comm(ip) = charge_comm(ip) + 1
            ko = ko + 1
            buf_comm_real(8*ko-7) = coord    (ilist,1)
            buf_comm_real(8*ko-6) = coord    (ilist,2)
            buf_comm_real(8*ko-5) = coord    (ilist,3)
            buf_comm_real(8*ko-4) = velocity (ilist,1)
            buf_comm_real(8*ko-3) = velocity (ilist,2)
            buf_comm_real(8*ko-2) = velocity (ilist,3)
            buf_comm_real(8*ko-1) = mass     (ilist  )
            buf_comm_real(8*ko  ) = charge   (ilist  )
            buf_comm_int (6*ko-5) = ip
            buf_comm_int (6*ko-4) = ig
            buf_comm_int (6*ko-3) = atmcls   (ilist  )
            buf_comm_int (6*ko-2) = id_l2g   (ilist  )
            buf_comm_int (6*ko-1) = chain_id (ilist  )
            buf_comm_int (6*ko  ) = atom_type(ilist  )

          end if

        else

          n_stay = n_stay + 1
          charge_stay(i) = charge_stay(i) + 1
          buf_stay_real(8*n_stay-7) = coord    (ilist,1)
          buf_stay_real(8*n_stay-6) = coord    (ilist,2)
          buf_stay_real(8*n_stay-5) = coord    (ilist,3)
          buf_stay_real(8*n_stay-4) = velocity (ilist,1)
          buf_stay_real(8*n_stay-3) = velocity (ilist,2)
          buf_stay_real(8*n_stay-2) = velocity (ilist,3)
          buf_stay_real(8*n_stay-1) = mass     (ilist  )
          buf_stay_real(8*n_stay  ) = charge   (ilist  )
          buf_stay_int (5*n_stay-4) = i
          buf_stay_int (5*n_stay-3) = atmcls   (ilist  )
          buf_stay_int (5*n_stay-2) = id_l2g   (ilist  )
          buf_stay_int (5*n_stay-1) = chain_id (ilist  )
          buf_stay_int (5*n_stay  ) = atom_type(ilist  )

        end if

      end if

    end do

    domain%n_stay1 = n_stay
    domain%charge_move_domain = ki
    domain%charge_comm_domain = ko

    return

  end subroutine update_outgoing_charge

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_nocharge
  !> @brief        Checks for non-charged particles moving to other cells.
  !! @authors      JJ
  !! @param[in]    boundary    : boundary condition information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_nocharge(boundary, domain)

    ! formal arguments
    type(s_boundary),    target, intent(in)    :: boundary
    type(s_domain),      target, intent(inout) :: domain

    ! local variable
    real(wip)                    :: x_shift, y_shift, z_shift
    real(wip)                    :: move(3)
    integer                      :: i, icx, icy, icz, icel, ncel
    integer                      :: ilist, ig, n_stay, n_stay_int, n_stay_real
    integer                      :: icel_local, icel_bd, ip, num
    integer                      :: ki_int, ki_real, ko_int, ko_real, ki, ko

    real(wip),           pointer :: bsize_x, bsize_y, bsize_z
    real(wip),           pointer :: csize_x, csize_y, csize_z
    real(wip),           pointer :: coord(:,:), velocity(:,:)
    real(wp),            pointer :: charge(:)
    real(wip),           pointer :: mass(:)
    real(wip),           pointer :: buf_stay_real(:)
    real(wip),           pointer :: buf_move_real(:)
    real(wip),           pointer :: buf_comm_real(:)
    integer,             pointer :: ncel_x, ncel_y, ncel_z, ncel_local, ncel_bd
    integer,             pointer :: natom(:), ncharge(:) 
    integer,             pointer :: cell_g2l(:), cell_g2b(:), cell_b2g(:)
    integer,             pointer :: atmcls(:), chain_id(:), atom_2_cell(:)
    integer,             pointer :: id_l2g(:), id_g2l(:)
    integer,             pointer :: atom_type(:)
    integer,             pointer :: nocharge_move(:), nocharge_stay(:)
    integer,             pointer :: nocharge_comm(:), nocharge_comm1(:)
    integer,             pointer :: cell_rank(:)
    integer,             pointer :: buf_stay_int(:)
    integer,             pointer :: buf_move_int(:)
    integer,             pointer :: buf_comm_int(:)


    bsize_x       => boundary%box_size_x
    bsize_y       => boundary%box_size_y
    bsize_z       => boundary%box_size_z
    ncel_x        => boundary%num_cells_x
    ncel_y        => boundary%num_cells_y
    ncel_z        => boundary%num_cells_z
    csize_x       => boundary%cell_size_x
    csize_y       => boundary%cell_size_y
    csize_z       => boundary%cell_size_z

    ncel_local    => domain%num_cell_local
    natom         => domain%num_atom
    ncharge       => domain%num_charge
    ncel_bd       => domain%num_cell_boundary
    cell_g2l      => domain%cell_g2l
    cell_g2b      => domain%cell_g2b
    cell_b2g      => domain%cell_b2g
    coord         => domain%coord
    velocity      => domain%velocity
    charge        => domain%charge
    mass          => domain%mass
    atmcls        => domain%atom_cls_no
    atom_2_cell   => domain%atom_2_cell
    id_l2g        => domain%id_l2g
    id_g2l        => domain%id_g2l
    chain_id      => domain%mol_chain_id
    atom_type     => domain%NA_base_type
    cell_rank     => domain%domain_cell_rank
    nocharge_move => domain%type2_move
    nocharge_stay => domain%type2_stay
    nocharge_comm => domain%type2_comm
    nocharge_comm1=> domain%type2_comm_move
    buf_stay_real => domain%buf_var0_stay_real
    buf_move_real => domain%buf_var0_move_real
    buf_comm_real => domain%buf_var0_comm_real
    buf_stay_int  => domain%buf_var0_stay_int 
    buf_move_int  => domain%buf_var0_move_int 
    buf_comm_int  => domain%buf_var0_comm_int 

    ! initializaiton
    !
    ncel = ncel_local + ncel_bd
    num  = domain%num_comm_proc

    nocharge_move(1:ncel) = 0
    nocharge_comm(1:num) = 0
    nocharge_stay(1:ncel_local) = 0

    n_stay_int  = domain%n_stay1 * 5
    n_stay_real = domain%n_stay1 * 8
    ki_int      = domain%charge_move_domain * 5
    ki_real     = domain%charge_move_domain * 8
    ko_int      = domain%charge_comm_domain * 6
    ko_real     = domain%charge_comm_domain * 8

    ! Check outgoing particles
    !
    n_stay = 0
    ki = 0
    ko = 0

    do ilist = 1, domain%num_atom_domain

      if (abs(charge(ilist)) < EPS) then

        i = atom_2_cell(ilist)
        x_shift = coord(ilist,1) - boundary%origin_x
        y_shift = coord(ilist,2) - boundary%origin_y
        z_shift = coord(ilist,3) - boundary%origin_z

        !coordinate shifted to the first quadrant and set into the boundary box
        if (boundary%type == BoundaryTypePBC) then
          move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
          move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
          move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)
          x_shift = x_shift + move(1)
          y_shift = y_shift + move(2)
          z_shift = z_shift + move(3)
        else
          x_shift = x_shift + bsize_x*0.5_wip
          y_shift = y_shift + bsize_y*0.5_wip
          z_shift = z_shift + bsize_z*0.5_wip
        end if

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

          if (icel_local /= 0) then

            ki = ki + 1
            nocharge_move(icel_local) = nocharge_move(icel_local) + 1
            buf_move_real(ki_real+8*ki-7) = coord    (ilist,1)
            buf_move_real(ki_real+8*ki-6) = coord    (ilist,2)
            buf_move_real(ki_real+8*ki-5) = coord    (ilist,3)
            buf_move_real(ki_real+8*ki-4) = velocity (ilist,1)
            buf_move_real(ki_real+8*ki-3) = velocity (ilist,2)
            buf_move_real(ki_real+8*ki-2) = velocity (ilist,3)
            buf_move_real(ki_real+8*ki-1) = mass     (ilist  )
            buf_move_real(ki_real+8*ki  ) = charge   (ilist  )
            buf_move_int (ki_int +5*ki-4) = icel_local
            buf_move_int (ki_int +5*ki-3) = atmcls   (ilist  )
            buf_move_int (ki_int +5*ki-2) = id_l2g   (ilist  )
            buf_move_int (ki_int +5*ki-1) = chain_id (ilist  ) 
            buf_move_int (ki_int +5*ki  ) = atom_type(ilist  ) 

          else if (icel_bd /= 0) then

            ig = cell_b2g(icel_bd)
            icel_bd = icel_bd + ncel_local
            ip = cell_rank(icel_bd)
            nocharge_comm(ip) = nocharge_comm(ip) + 1
            ko = ko + 1
            buf_comm_real(ko_real+8*ko-7) = coord    (ilist,1)
            buf_comm_real(ko_real+8*ko-6) = coord    (ilist,2)
            buf_comm_real(ko_real+8*ko-5) = coord    (ilist,3)
            buf_comm_real(ko_real+8*ko-4) = velocity (ilist,1)
            buf_comm_real(ko_real+8*ko-3) = velocity (ilist,2)
            buf_comm_real(ko_real+8*ko-2) = velocity (ilist,3)
            buf_comm_real(ko_real+8*ko-1) = mass     (ilist  )
            buf_comm_real(ko_real+8*ko  ) = charge   (ilist  )
            buf_comm_int (ko_int +6*ko-5) = ip
            buf_comm_int (ko_int +6*ko-4) = ig
            buf_comm_int (ko_int +6*ko-3) = atmcls   (ilist  )
            buf_comm_int (ko_int +6*ko-2) = id_l2g   (ilist  )
            buf_comm_int (ko_int +6*ko-1) = chain_id (ilist  )
            buf_comm_int (ko_int +6*ko  ) = atom_type(ilist  )

          end if

        else

          n_stay = n_stay + 1
          nocharge_stay(i) = nocharge_stay(i) + 1
          buf_stay_real(n_stay_real+8*n_stay-7) = coord    (ilist,1)
          buf_stay_real(n_stay_real+8*n_stay-6) = coord    (ilist,2)
          buf_stay_real(n_stay_real+8*n_stay-5) = coord    (ilist,3)
          buf_stay_real(n_stay_real+8*n_stay-4) = velocity (ilist,1)
          buf_stay_real(n_stay_real+8*n_stay-3) = velocity (ilist,2)
          buf_stay_real(n_stay_real+8*n_stay-2) = velocity (ilist,3)
          buf_stay_real(n_stay_real+8*n_stay-1) = mass     (ilist  )
          buf_stay_real(n_stay_real+8*n_stay  ) = charge   (ilist  )
          buf_stay_int (n_stay_int +5*n_stay-4) = i
          buf_stay_int (n_stay_int +5*n_stay-3) = atmcls   (ilist  )
          buf_stay_int (n_stay_int +5*n_stay-2) = id_l2g   (ilist  )
          buf_stay_int (n_stay_int +5*n_stay-1) = chain_id (ilist  )
          buf_stay_int (n_stay_int +5*n_stay  ) = atom_type(ilist  )

        end if

      end if

    end do

    domain%n_stay2 = n_stay
    domain%nocharge_move_domain = ki
    domain%nocharge_comm_domain = ko

    return

  end subroutine update_outgoing_nocharge

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_ptl
  !> @brief        Checks for charged particles incoming to each cell.
  !! @authors      JJ
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_incoming_ptl(domain)

    ! formal arguments
    type(s_domain),      target, intent(inout) :: domain

    ! local variables
    integer                      :: i, j, k, ig, start_i, ilist, n_stay
    integer                      :: n_stay_int, n_stay_real
    integer                      :: charge_move_domain, nocharge_move_domain
    integer                      :: id, ip, num

    real(wip),        pointer    :: coord(:,:), velocity(:,:)
    real(wp),         pointer    :: charge(:)
    real(wip),        pointer    :: mass(:)
    real(wip),        pointer    :: buf_move_real(:), buf_comm_real(:)
    real(wip),        pointer    :: buf_stay_real(:)
    integer,          pointer    :: ncel_local, ncel_bd 
    integer,          pointer    :: natom(:), ncharge(:), start_atom(:)
    integer,          pointer    :: atmcls(:), id_l2g(:), id_g2l(:)
    integer,          pointer    :: chain_id(:)
    integer,          pointer    :: atom_type(:)
    integer,          pointer    :: cell_g2l(:)
    integer,          pointer    :: charge_move(:), charge_stay(:)
    integer,          pointer    :: nocharge_move(:), nocharge_stay(:)
    integer,          pointer    :: charge_comm(:), nocharge_comm(:)
    integer,          pointer    :: charge_comm1(:), nocharge_comm1(:)
    integer,          pointer    :: buf_move_int(:), buf_comm_int(:)
    integer,          pointer    :: buf_stay_int(:)

    ncel_local             => domain%num_cell_local
    natom                  => domain%num_atom
    ncharge                => domain%num_charge
    start_atom             => domain%start_atom
    ncel_bd                => domain%num_cell_boundary
    coord                  => domain%coord
    velocity               => domain%velocity
    charge                 => domain%charge
    mass                   => domain%mass
    atmcls                 => domain%atom_cls_no
    id_l2g                 => domain%id_l2g
    id_g2l                 => domain%id_g2l
    atom_type              => domain%NA_base_type
    chain_id               => domain%mol_chain_id
    cell_g2l               => domain%cell_g2l
    charge_move            => domain%type1_move
    charge_comm1           => domain%type1_comm_move
    charge_comm            => domain%type1_comm
    charge_stay            => domain%type1_stay
    nocharge_move          => domain%type2_move
    nocharge_comm1         => domain%type2_comm_move
    nocharge_comm          => domain%type2_comm
    nocharge_stay          => domain%type2_stay
    buf_stay_int           => domain%buf_var0_stay_int
    buf_stay_real          => domain%buf_var0_stay_real
    buf_move_real          => domain%buf_var0_move_real
    buf_move_int           => domain%buf_var0_move_int
    buf_comm_real          => domain%buf_var0_comm_real
    buf_comm_int           => domain%buf_var0_comm_int

    charge_move_domain   = domain%charge_move_domain
    nocharge_move_domain = domain%nocharge_move_domain
    num = domain%num_comm_proc

    do i = 1, ncel_local
      ncharge(i) = charge_stay(i) + charge_move(i)
      natom(i)   = ncharge(i) + nocharge_stay(i) + nocharge_move(i)
    end do
    start_atom(1:ncel_local+ncel_bd) = 0
    k = 0
    do i = 1, ncel_local-1
      k = k + natom(i)
      start_atom(i+1) = k
    end do
    k = k + natom(ncel_local)
    domain%num_atom_domain = k 

    !$omp parallel private(id)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    call get_para_range(1, k, nthread, id, domain%start_index(id), &
                        domain%end_index(id))
    !$omp end parallel

#ifdef DEBUG
    do i = 1, ncel_local
      if (start_atom(i)+natom(i) > MaxAtom_domain) &
        call error_msg('Debug: Update_Incoming_Atom: Atom number exceeds '//&
                       'MaxAtom_domain.')
    end do
#endif

    ! incoming charge
    !
!   !$omp parallel do private(n_stay, i, k, ilist, ig)
    charge_stay(1:ncel_local) = 0
    do n_stay = 1, domain%n_stay1
      i = buf_stay_int(5*n_stay-4)
      charge_stay(i) = charge_stay(i) + 1
      ilist = charge_stay(i) + start_atom(i)
      coord    (ilist,1)  = buf_stay_real(8*n_stay-7)
      coord    (ilist,2)  = buf_stay_real(8*n_stay-6)
      coord    (ilist,3)  = buf_stay_real(8*n_stay-5)
      velocity (ilist,1)  = buf_stay_real(8*n_stay-4)
      velocity (ilist,2)  = buf_stay_real(8*n_stay-3)
      velocity (ilist,3)  = buf_stay_real(8*n_stay-2)
      mass     (ilist  )  = buf_stay_real(8*n_stay-1)
      charge   (ilist  )  = buf_stay_real(8*n_stay  )
      atmcls   (ilist  )  = buf_stay_int (5*n_stay-3)
      id_l2g   (ilist  )  = buf_stay_int (5*n_stay-2)
      chain_id (ilist  )  = buf_stay_int (5*n_stay-1)
      atom_type(ilist  )  = buf_stay_int (5*n_stay  )
      ig = id_l2g(ilist)
      id_g2l(ig) = ilist
    end do
!   !$omp end parallel do

    n_stay_real = domain%n_stay1 * 8
    n_stay_int  = domain%n_stay1 * 5
    nocharge_stay(1:ncel_local) = 0
!   !$omp parallel do private(n_stay, i, k, ilist, ig)
    do n_stay = 1, domain%n_stay2
      i = buf_stay_int(n_stay_int+5*n_stay-4)
      nocharge_stay(i) = nocharge_stay(i) + 1
      ilist = ncharge(i) + nocharge_stay(i) + start_atom(i)
      coord    (ilist,1)  = buf_stay_real(n_stay_real+8*n_stay-7)
      coord    (ilist,2)  = buf_stay_real(n_stay_real+8*n_stay-6)
      coord    (ilist,3)  = buf_stay_real(n_stay_real+8*n_stay-5)
      velocity (ilist,1)  = buf_stay_real(n_stay_real+8*n_stay-4)
      velocity (ilist,2)  = buf_stay_real(n_stay_real+8*n_stay-3)
      velocity (ilist,3)  = buf_stay_real(n_stay_real+8*n_stay-2)
      mass     (ilist  )  = buf_stay_real(n_stay_real+8*n_stay-1)
      charge   (ilist  )  = buf_stay_real(n_stay_real+8*n_stay  )
      atmcls   (ilist  )  = buf_stay_int (n_stay_int +5*n_stay-3)
      id_l2g   (ilist  )  = buf_stay_int (n_stay_int +5*n_stay-2)
      chain_id (ilist  )  = buf_stay_int (n_stay_int +5*n_stay-1)
      atom_type(ilist  )  = buf_stay_int (n_stay_int +5*n_stay  )
      ig = id_l2g(ilist)
      id_g2l(ig) = ilist
    end do
!   !$omp end parallel do

    charge_comm1(1:ncel_local+ncel_bd) = 0
    nocharge_comm1(1:ncel_local+ncel_bd) = 0

    do k = 1, charge_move_domain
      i = buf_move_int(5*k-4)
      start_i = start_atom(i)
      charge_comm1(i) = charge_comm1(i) + 1
      ilist = charge_comm1(i) + charge_stay(i) + start_i
      coord    (ilist,1)  = buf_move_real(8*k-7)
      coord    (ilist,2)  = buf_move_real(8*k-6)
      coord    (ilist,3)  = buf_move_real(8*k-5)
      velocity (ilist,1)  = buf_move_real(8*k-4)
      velocity (ilist,2)  = buf_move_real(8*k-3)
      velocity (ilist,3)  = buf_move_real(8*k-2)
      mass     (ilist  )  = buf_move_real(8*k-1)
      charge   (ilist  )  = buf_move_real(8*k  )
      atmcls   (ilist  )  = buf_move_int (5*k-3)
      id_l2g   (ilist  )  = buf_move_int (5*k-2)
      chain_id (ilist  )  = buf_move_int (5*k-1)
      atom_type(ilist  )  = buf_move_int (5*k  )
      ig = id_l2g(ilist)
      id_g2l(ig) = ilist
    end do

    do k = charge_move_domain+1, charge_move_domain+nocharge_move_domain
      i = buf_move_int(5*k-4)
      start_i = start_atom(i)
      nocharge_comm1(i) = nocharge_comm1(i) + 1
      ilist = nocharge_comm1(i) + nocharge_stay(i) + ncharge(i) + start_i
      coord    (ilist,1)  = buf_move_real(8*k-7)
      coord    (ilist,2)  = buf_move_real(8*k-6)
      coord    (ilist,3)  = buf_move_real(8*k-5)
      velocity (ilist,1)  = buf_move_real(8*k-4)
      velocity (ilist,2)  = buf_move_real(8*k-3)
      velocity (ilist,3)  = buf_move_real(8*k-2)
      mass     (ilist  )  = buf_move_real(8*k-1)
      charge   (ilist  )  = buf_move_real(8*k  )
      atmcls   (ilist  )  = buf_move_int (5*k-3)
      id_l2g   (ilist  )  = buf_move_int (5*k-2)
      chain_id (ilist  )  = buf_move_int (5*k-1)
      atom_type(ilist  )  = buf_move_int (5*k  )
      ig = id_l2g(ilist)
      id_g2l(ig) = ilist
    end do

    j = 0
    do ip = 1, num
      do k = 1, charge_comm(ip)
        j = j + 1
        i = cell_g2l(buf_comm_int(5*j-4))
        start_i = start_atom(i)
        charge_comm1(i) = charge_comm1(i) + 1
        ilist = charge_comm1(i) + charge_stay(i) + start_i
        coord    (ilist,1)  = buf_comm_real(8*j-7)
        coord    (ilist,2)  = buf_comm_real(8*j-6)
        coord    (ilist,3)  = buf_comm_real(8*j-5)
        velocity (ilist,1)  = buf_comm_real(8*j-4)
        velocity (ilist,2)  = buf_comm_real(8*j-3)
        velocity (ilist,3)  = buf_comm_real(8*j-2)
        mass     (ilist  )  = buf_comm_real(8*j-1)
        charge   (ilist  )  = buf_comm_real(8*j  )
        atmcls   (ilist  )  = buf_comm_int (5*j-3)
        id_l2g   (ilist  )  = buf_comm_int (5*j-2)
        chain_id (ilist  )  = buf_comm_int (5*j-1)
        atom_type(ilist  )  = buf_comm_int (5*j  )
        ig = id_l2g(ilist)
        id_g2l(ig) = ilist
      end do
  
      do k = charge_comm(ip)+1,charge_comm(ip)+nocharge_comm(ip)
        j = j + 1 
        i = cell_g2l(buf_comm_int(5*j-4))
        start_i = start_atom(i)
        nocharge_comm1(i) = nocharge_comm1(i) + 1
        ilist = nocharge_comm1(i) + nocharge_stay(i) + ncharge(i) + start_i
        coord    (ilist,1)  = buf_comm_real(8*j-7)    
        coord    (ilist,2)  = buf_comm_real(8*j-6)    
        coord    (ilist,3)  = buf_comm_real(8*j-5)    
        velocity (ilist,1)  = buf_comm_real(8*j-4)
        velocity (ilist,2)  = buf_comm_real(8*j-3)
        velocity (ilist,3)  = buf_comm_real(8*j-2)
        mass     (ilist  )  = buf_comm_real(8*j-1)
        charge   (ilist  )  = buf_comm_real(8*j  )
        atmcls   (ilist  )  = buf_comm_int (5*j-3)
        id_l2g   (ilist  )  = buf_comm_int (5*j-2)
        chain_id (ilist  )  = buf_comm_int (5*j-1)
        atom_type(ilist  )  = buf_comm_int (5*j  )
        ig = id_l2g(ilist)
        id_g2l(ig) = ilist
      end do
    end do

    return

  end subroutine update_incoming_ptl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_charge_martini
  !> @brief        Checks for charged particles moving to other cells (Martini case).
  !! @authors      JJ
  !! @param[in]    boundary    : boundary condition information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_charge_martini(boundary, domain)

    ! formal arguments
    type(s_boundary),    target, intent(in)    :: boundary
    type(s_domain),      target, intent(inout) :: domain

    ! local variable
    real(wip)                    :: x_shift, y_shift, z_shift
    real(wip)                    :: move(3)
    integer                      :: i, icx, icy, icz, icel, ncel
    integer                      :: ilist, ig, n_stay, num
    integer                      :: icel_local, icel_bd, ip
    integer                      :: ki, ko

    real(wip),           pointer :: bsize_x, bsize_y, bsize_z
    real(wip),           pointer :: csize_x, csize_y, csize_z
    real(wip),           pointer :: coord(:,:), velocity(:,:)
    real(wp),            pointer :: charge(:)
    real(wip),           pointer :: mass(:)
    real(wip),           pointer :: buf_stay_real(:)
    real(wip),           pointer :: buf_move_real(:)
    real(wip),           pointer :: buf_comm_real(:)
    integer,             pointer :: ncel_x, ncel_y, ncel_z, ncel_local, ncel_bd
    integer,             pointer :: natom(:), ncharge(:)
    integer,             pointer :: cell_g2l(:), cell_g2b(:)
    integer,             pointer :: cell_b2g(:)
    integer,             pointer :: atmcls(:), atom_2_cell(:)
    integer,             pointer :: id_l2g(:), id_g2l(:)
    integer,             pointer :: charge_move(:), charge_comm(:)
    integer,             pointer :: charge_stay(:), charge_comm1(:)
    integer,             pointer :: cell_rank(:)
    integer,             pointer :: buf_stay_int(:)
    integer,             pointer :: buf_move_int(:)
    integer,             pointer :: buf_comm_int(:)

    bsize_x       => boundary%box_size_x
    bsize_y       => boundary%box_size_y
    bsize_z       => boundary%box_size_z
    ncel_x        => boundary%num_cells_x
    ncel_y        => boundary%num_cells_y
    ncel_z        => boundary%num_cells_z
    csize_x       => boundary%cell_size_x
    csize_y       => boundary%cell_size_y
    csize_z       => boundary%cell_size_z

    ncel_local    => domain%num_cell_local
    natom         => domain%num_atom
    ncharge       => domain%num_charge
    ncel_bd       => domain%num_cell_boundary
    cell_g2l      => domain%cell_g2l
    cell_g2b      => domain%cell_g2b
    cell_b2g      => domain%cell_b2g
    coord         => domain%coord
    velocity      => domain%velocity
    charge        => domain%charge
    mass          => domain%mass
    atmcls        => domain%atom_cls_no
    atom_2_cell   => domain%atom_2_cell
    id_l2g        => domain%id_l2g
    id_g2l        => domain%id_g2l
    cell_rank     => domain%domain_cell_rank
    charge_move   => domain%type1_move
    charge_comm   => domain%type1_comm
    charge_comm1  => domain%type1_comm_move
    charge_stay   => domain%type1_stay
    buf_stay_real => domain%buf_var0_stay_real
    buf_move_real => domain%buf_var0_move_real
    buf_comm_real => domain%buf_var0_comm_real
    buf_stay_int  => domain%buf_var0_stay_int 
    buf_move_int  => domain%buf_var0_move_int 
    buf_comm_int  => domain%buf_var0_comm_int 

    ! initializaiton
    !
    ncel = ncel_local + ncel_bd
    num  = domain%num_comm_proc
    charge_move(1:ncel) = 0
    charge_comm(1:num) = 0
    charge_comm1(1:num) = 0
    charge_stay(1:ncel_local) = 0

    ! initialized the global index
    !
    !$omp parallel do private(ig,i)
    do i = 1, domain%num_atom_domain+domain%num_atom_boundary
      ig = id_l2g(i)
      id_g2l(ig) = 0
    end do
    !$omp end parallel do

    ! Check outgoing particles
    !
    n_stay = 0
    ki     = 0
    ko     = 0

    do ilist = 1, domain%num_atom_domain

      if (abs(charge(ilist)) > EPS) then

        i = atom_2_cell(ilist)
        x_shift = coord(ilist,1) - boundary%origin_x
        y_shift = coord(ilist,2) - boundary%origin_y
        z_shift = coord(ilist,3) - boundary%origin_z

        !coordinate shifted to the first quadrant and set into the boundary box
        if (boundary%type == BoundaryTypePBC) then
          move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
          move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
          move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)
          x_shift = x_shift + move(1)
          y_shift = y_shift + move(2)
          z_shift = z_shift + move(3)
        else
          x_shift = x_shift + bsize_x*0.5_wip
          y_shift = y_shift + bsize_y*0.5_wip
          z_shift = z_shift + bsize_z*0.5_wip
        end if

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

          if (icel_local /= 0) then

            ki = ki + 1
            charge_move(icel_local) = charge_move(icel_local) + 1
            buf_move_real(8*ki-7) = coord    (ilist,1)
            buf_move_real(8*ki-6) = coord    (ilist,2)
            buf_move_real(8*ki-5) = coord    (ilist,3)
            buf_move_real(8*ki-4) = velocity (ilist,1)
            buf_move_real(8*ki-3) = velocity (ilist,2)
            buf_move_real(8*ki-2) = velocity (ilist,3)
            buf_move_real(8*ki-1) = mass     (ilist  )
            buf_move_real(8*ki  ) = charge   (ilist  )
            buf_move_int (3*ki-2) = icel_local
            buf_move_int (3*ki-1) = atmcls   (ilist  )
            buf_move_int (3*ki-0) = id_l2g   (ilist  )

          else if (icel_bd /= 0) then

            ig = cell_b2g(icel_bd)
            icel_bd = icel_bd + ncel_local
            ip = cell_rank(icel_bd)
            charge_comm(ip) = charge_comm(ip) + 1
            ko = ko + 1
            buf_comm_real(8*ko-7) = coord    (ilist,1)
            buf_comm_real(8*ko-6) = coord    (ilist,2)
            buf_comm_real(8*ko-5) = coord    (ilist,3)
            buf_comm_real(8*ko-4) = velocity (ilist,1)
            buf_comm_real(8*ko-3) = velocity (ilist,2)
            buf_comm_real(8*ko-2) = velocity (ilist,3)
            buf_comm_real(8*ko-1) = mass     (ilist  )
            buf_comm_real(8*ko  ) = charge   (ilist  )
            buf_comm_int (4*ko-3) = ip
            buf_comm_int (4*ko-2) = ig
            buf_comm_int (4*ko-1) = atmcls   (ilist  )
            buf_comm_int (4*ko-0) = id_l2g   (ilist  )

          end if

        else

          n_stay = n_stay + 1
          charge_stay(i) = charge_stay(i) + 1
          buf_stay_real(8*n_stay-7) = coord    (ilist,1)
          buf_stay_real(8*n_stay-6) = coord    (ilist,2)
          buf_stay_real(8*n_stay-5) = coord    (ilist,3)
          buf_stay_real(8*n_stay-4) = velocity (ilist,1)
          buf_stay_real(8*n_stay-3) = velocity (ilist,2)
          buf_stay_real(8*n_stay-2) = velocity (ilist,3)
          buf_stay_real(8*n_stay-1) = mass     (ilist  )
          buf_stay_real(8*n_stay  ) = charge   (ilist  )
          buf_stay_int (3*n_stay-2) = i
          buf_stay_int (3*n_stay-1) = atmcls   (ilist  )
          buf_stay_int (3*n_stay-0) = id_l2g   (ilist  )

        end if

      end if

    end do

    domain%n_stay1 = n_stay
    domain%charge_move_domain = ki
    domain%charge_comm_domain = ko

    return

  end subroutine update_outgoing_charge_martini

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_nocharge_martini
  !> @brief        Checks for non-charged particles moving to other cells (Martini).
  !! @authors      JJ
  !! @param[in]    boundary    : boundary condition information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_nocharge_martini(boundary, domain)

    ! formal arguments
    type(s_boundary),    target, intent(in)    :: boundary
    type(s_domain),      target, intent(inout) :: domain

    ! local variable
    real(wip)                    :: x_shift, y_shift, z_shift
    real(wip)                    :: move(3)
    integer                      :: i, icx, icy, icz, icel, ncel
    integer                      :: ilist, ig, n_stay, n_stay_int, n_stay_real
    integer                      :: icel_local, icel_bd, ip, num
    integer                      :: ki_int, ki_real, ko_int, ko_real, ki, ko

    real(wip),           pointer :: bsize_x, bsize_y, bsize_z
    real(wip),           pointer :: csize_x, csize_y, csize_z
    real(wip),           pointer :: coord(:,:), velocity(:,:)
    real(wp),            pointer :: charge(:)
    real(wip),           pointer :: mass(:)
    real(wip),           pointer :: buf_stay_real(:)
    real(wip),           pointer :: buf_move_real(:)
    real(wip),           pointer :: buf_comm_real(:)
    integer,             pointer :: ncel_x, ncel_y, ncel_z, ncel_local, ncel_bd
    integer,             pointer :: natom(:), ncharge(:) 
    integer,             pointer :: cell_g2l(:), cell_g2b(:), cell_b2g(:)
    integer,             pointer :: atmcls(:), atom_2_cell(:)
    integer,             pointer :: id_l2g(:), id_g2l(:)
    integer,             pointer :: nocharge_move(:), nocharge_stay(:)
    integer,             pointer :: nocharge_comm(:), nocharge_comm1(:)
    integer,             pointer :: cell_rank(:)
    integer,             pointer :: buf_stay_int(:)
    integer,             pointer :: buf_move_int(:)
    integer,             pointer :: buf_comm_int(:)


    bsize_x       => boundary%box_size_x
    bsize_y       => boundary%box_size_y
    bsize_z       => boundary%box_size_z
    ncel_x        => boundary%num_cells_x
    ncel_y        => boundary%num_cells_y
    ncel_z        => boundary%num_cells_z
    csize_x       => boundary%cell_size_x
    csize_y       => boundary%cell_size_y
    csize_z       => boundary%cell_size_z

    ncel_local    => domain%num_cell_local
    natom         => domain%num_atom
    ncharge       => domain%num_charge
    ncel_bd       => domain%num_cell_boundary
    cell_g2l      => domain%cell_g2l
    cell_g2b      => domain%cell_g2b
    cell_b2g      => domain%cell_b2g
    coord         => domain%coord
    velocity      => domain%velocity
    charge        => domain%charge
    mass          => domain%mass
    atmcls        => domain%atom_cls_no
    atom_2_cell   => domain%atom_2_cell
    id_l2g        => domain%id_l2g
    id_g2l        => domain%id_g2l
    cell_rank     => domain%domain_cell_rank
    nocharge_move => domain%type2_move
    nocharge_stay => domain%type2_stay
    nocharge_comm => domain%type2_comm
    nocharge_comm1=> domain%type2_comm_move
    buf_stay_real => domain%buf_var0_stay_real
    buf_move_real => domain%buf_var0_move_real
    buf_comm_real => domain%buf_var0_comm_real
    buf_stay_int  => domain%buf_var0_stay_int 
    buf_move_int  => domain%buf_var0_move_int 
    buf_comm_int  => domain%buf_var0_comm_int 

    ! initializaiton
    !
    ncel = ncel_local + ncel_bd
    num  = domain%num_comm_proc

    nocharge_move(1:ncel) = 0
    nocharge_comm(1:num) = 0
    nocharge_stay(1:ncel_local) = 0

    n_stay_int  = domain%n_stay1 * 3
    n_stay_real = domain%n_stay1 * 8
    ki_int      = domain%charge_move_domain * 3
    ki_real     = domain%charge_move_domain * 8
    ko_int      = domain%charge_comm_domain * 4
    ko_real     = domain%charge_comm_domain * 8

    ! Check outgoing particles
    !
    n_stay = 0
    ki = 0
    ko = 0

    do ilist = 1, domain%num_atom_domain

      if (abs(charge(ilist)) < EPS) then

        i = atom_2_cell(ilist)
        x_shift = coord(ilist,1) - boundary%origin_x
        y_shift = coord(ilist,2) - boundary%origin_y
        z_shift = coord(ilist,3) - boundary%origin_z

        !coordinate shifted to the first quadrant and set into the boundary box
        if (boundary%type == BoundaryTypePBC) then
          move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
          move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
          move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)
          x_shift = x_shift + move(1)
          y_shift = y_shift + move(2)
          z_shift = z_shift + move(3)
        else
          x_shift = x_shift + bsize_x*0.5_wip
          y_shift = y_shift + bsize_y*0.5_wip
          z_shift = z_shift + bsize_z*0.5_wip
        end if

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

          if (icel_local /= 0) then

            ki = ki + 1
            nocharge_move(icel_local) = nocharge_move(icel_local) + 1
            buf_move_real(ki_real+8*ki-7) = coord    (ilist,1)
            buf_move_real(ki_real+8*ki-6) = coord    (ilist,2)
            buf_move_real(ki_real+8*ki-5) = coord    (ilist,3)
            buf_move_real(ki_real+8*ki-4) = velocity (ilist,1)
            buf_move_real(ki_real+8*ki-3) = velocity (ilist,2)
            buf_move_real(ki_real+8*ki-2) = velocity (ilist,3)
            buf_move_real(ki_real+8*ki-1) = mass     (ilist  )
            buf_move_real(ki_real+8*ki  ) = charge   (ilist  )
            buf_move_int (ki_int +3*ki-2) = icel_local
            buf_move_int (ki_int +3*ki-1) = atmcls   (ilist  )
            buf_move_int (ki_int +3*ki-0) = id_l2g   (ilist  )

          else if (icel_bd /= 0) then

            ig = cell_b2g(icel_bd)
            icel_bd = icel_bd + ncel_local
            ip = cell_rank(icel_bd)
            nocharge_comm(ip) = nocharge_comm(ip) + 1
            ko = ko + 1
            buf_comm_real(ko_real+8*ko-7) = coord    (ilist,1)
            buf_comm_real(ko_real+8*ko-6) = coord    (ilist,2)
            buf_comm_real(ko_real+8*ko-5) = coord    (ilist,3)
            buf_comm_real(ko_real+8*ko-4) = velocity (ilist,1)
            buf_comm_real(ko_real+8*ko-3) = velocity (ilist,2)
            buf_comm_real(ko_real+8*ko-2) = velocity (ilist,3)
            buf_comm_real(ko_real+8*ko-1) = mass     (ilist  )
            buf_comm_real(ko_real+8*ko  ) = charge   (ilist  )
            buf_comm_int (ko_int +4*ko-3) = ip
            buf_comm_int (ko_int +4*ko-2) = ig
            buf_comm_int (ko_int +4*ko-1) = atmcls   (ilist  )
            buf_comm_int (ko_int +4*ko-0) = id_l2g   (ilist  )

          end if

        else

          n_stay = n_stay + 1
          nocharge_stay(i) = nocharge_stay(i) + 1
          buf_stay_real(n_stay_real+8*n_stay-7) = coord    (ilist,1)
          buf_stay_real(n_stay_real+8*n_stay-6) = coord    (ilist,2)
          buf_stay_real(n_stay_real+8*n_stay-5) = coord    (ilist,3)
          buf_stay_real(n_stay_real+8*n_stay-4) = velocity (ilist,1)
          buf_stay_real(n_stay_real+8*n_stay-3) = velocity (ilist,2)
          buf_stay_real(n_stay_real+8*n_stay-2) = velocity (ilist,3)
          buf_stay_real(n_stay_real+8*n_stay-1) = mass     (ilist  )
          buf_stay_real(n_stay_real+8*n_stay  ) = charge   (ilist  )
          buf_stay_int (n_stay_int +3*n_stay-2) = i
          buf_stay_int (n_stay_int +3*n_stay-1) = atmcls   (ilist  )
          buf_stay_int (n_stay_int +3*n_stay-0) = id_l2g   (ilist  )

        end if

      end if

    end do

    domain%n_stay2 = n_stay
    domain%nocharge_move_domain = ki
    domain%nocharge_comm_domain = ko

    return

  end subroutine update_outgoing_nocharge_martini

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_ptl_martini
  !> @brief        Checks for charged particles incoming to each cell (Martini).
  !! @authors      JJ
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_incoming_ptl_martini(domain)

    ! formal arguments
    type(s_domain),      target, intent(inout) :: domain

    ! local variables
    integer                      :: i, j, k, ig, start_i, ilist, n_stay
    integer                      :: n_stay_int, n_stay_real
    integer                      :: charge_move_domain, nocharge_move_domain
    integer                      :: id, ip, num

    real(wip),        pointer    :: coord(:,:), velocity(:,:)
    real(wp),         pointer    :: charge(:)
    real(wip),        pointer    :: mass(:)
    real(wip),        pointer    :: buf_move_real(:), buf_comm_real(:)
    real(wip),        pointer    :: buf_stay_real(:)
    integer,          pointer    :: ncel_local, ncel_bd 
    integer,          pointer    :: natom(:), ncharge(:), start_atom(:)
    integer,          pointer    :: atmcls(:), id_l2g(:), id_g2l(:)
    integer,          pointer    :: cell_g2l(:)
    integer,          pointer    :: charge_move(:), charge_stay(:)
    integer,          pointer    :: nocharge_move(:), nocharge_stay(:)
    integer,          pointer    :: charge_comm(:), nocharge_comm(:)
    integer,          pointer    :: charge_comm1(:), nocharge_comm1(:)
    integer,          pointer    :: buf_move_int(:), buf_comm_int(:)
    integer,          pointer    :: buf_stay_int(:)

    ncel_local             => domain%num_cell_local
    natom                  => domain%num_atom
    ncharge                => domain%num_charge
    start_atom             => domain%start_atom
    ncel_bd                => domain%num_cell_boundary
    coord                  => domain%coord
    velocity               => domain%velocity
    charge                 => domain%charge
    mass                   => domain%mass
    atmcls                 => domain%atom_cls_no
    id_l2g                 => domain%id_l2g
    id_g2l                 => domain%id_g2l
    cell_g2l               => domain%cell_g2l
    charge_move            => domain%type1_move
    charge_comm1           => domain%type1_comm_move
    charge_comm            => domain%type1_comm
    charge_stay            => domain%type1_stay
    nocharge_move          => domain%type2_move
    nocharge_comm1         => domain%type2_comm_move
    nocharge_comm          => domain%type2_comm
    nocharge_stay          => domain%type2_stay
    buf_stay_int           => domain%buf_var0_stay_int
    buf_stay_real          => domain%buf_var0_stay_real
    buf_move_real          => domain%buf_var0_move_real
    buf_move_int           => domain%buf_var0_move_int
    buf_comm_real          => domain%buf_var0_comm_real
    buf_comm_int           => domain%buf_var0_comm_int

    charge_move_domain   = domain%charge_move_domain
    nocharge_move_domain = domain%nocharge_move_domain
    num = domain%num_comm_proc

    do i = 1, ncel_local
      ncharge(i) = charge_stay(i) + charge_move(i)
      natom(i)   = ncharge(i) + nocharge_stay(i) + nocharge_move(i)
    end do
    start_atom(1:ncel_local+ncel_bd) = 0
    k = 0
    do i = 1, ncel_local-1
      k = k + natom(i)
      start_atom(i+1) = k
    end do
    k = k + natom(ncel_local)
    domain%num_atom_domain = k 

    !$omp parallel private(id)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    call get_para_range(1, k, nthread, id, domain%start_index(id), &
                        domain%end_index(id))
    !$omp end parallel

#ifdef DEBUG
    do i = 1, ncel_local
      if (start_atom(i)+natom(i) > MaxAtom_domain) &
        call error_msg('Debug: Update_Incoming_Atom: Atom number exceeds '//&
                       'MaxAtom_domain.')
    end do
#endif

    ! incoming charge
    !
!   !$omp parallel do private(n_stay, i, k, ilist, ig)
    charge_stay(1:ncel_local) = 0
    do n_stay = 1, domain%n_stay1
      i = buf_stay_int(3*n_stay-2)
      charge_stay(i) = charge_stay(i) + 1
      ilist = charge_stay(i) + start_atom(i)
      coord    (ilist,1)  = buf_stay_real(8*n_stay-7)
      coord    (ilist,2)  = buf_stay_real(8*n_stay-6)
      coord    (ilist,3)  = buf_stay_real(8*n_stay-5)
      velocity (ilist,1)  = buf_stay_real(8*n_stay-4)
      velocity (ilist,2)  = buf_stay_real(8*n_stay-3)
      velocity (ilist,3)  = buf_stay_real(8*n_stay-2)
      mass     (ilist  )  = buf_stay_real(8*n_stay-1)
      charge   (ilist  )  = buf_stay_real(8*n_stay  )
      atmcls   (ilist  )  = buf_stay_int (3*n_stay-1)
      id_l2g   (ilist  )  = buf_stay_int (3*n_stay-0)
      ig = id_l2g(ilist)
      id_g2l(ig) = ilist
    end do
!   !$omp end parallel do

    n_stay_real = domain%n_stay1 * 8
    n_stay_int  = domain%n_stay1 * 3
    nocharge_stay(1:ncel_local) = 0
!   !$omp parallel do private(n_stay, i, k, ilist, ig)
    do n_stay = 1, domain%n_stay2
      i = buf_stay_int(n_stay_int+3*n_stay-2)
      nocharge_stay(i) = nocharge_stay(i) + 1
      ilist = ncharge(i) + nocharge_stay(i) + start_atom(i)
      coord    (ilist,1)  = buf_stay_real(n_stay_real+8*n_stay-7)
      coord    (ilist,2)  = buf_stay_real(n_stay_real+8*n_stay-6)
      coord    (ilist,3)  = buf_stay_real(n_stay_real+8*n_stay-5)
      velocity (ilist,1)  = buf_stay_real(n_stay_real+8*n_stay-4)
      velocity (ilist,2)  = buf_stay_real(n_stay_real+8*n_stay-3)
      velocity (ilist,3)  = buf_stay_real(n_stay_real+8*n_stay-2)
      mass     (ilist  )  = buf_stay_real(n_stay_real+8*n_stay-1)
      charge   (ilist  )  = buf_stay_real(n_stay_real+8*n_stay  )
      atmcls   (ilist  )  = buf_stay_int (n_stay_int +3*n_stay-1)
      id_l2g   (ilist  )  = buf_stay_int (n_stay_int +3*n_stay-0)
      ig = id_l2g(ilist)
      id_g2l(ig) = ilist
    end do
!   !$omp end parallel do

    charge_comm1(1:ncel_local+ncel_bd) = 0
    nocharge_comm1(1:ncel_local+ncel_bd) = 0

    do k = 1, charge_move_domain
      i = buf_move_int(3*k-2)
      start_i = start_atom(i)
      charge_comm1(i) = charge_comm1(i) + 1
      ilist = charge_comm1(i) + charge_stay(i) + start_i
      coord    (ilist,1)  = buf_move_real(8*k-7)
      coord    (ilist,2)  = buf_move_real(8*k-6)
      coord    (ilist,3)  = buf_move_real(8*k-5)
      velocity (ilist,1)  = buf_move_real(8*k-4)
      velocity (ilist,2)  = buf_move_real(8*k-3)
      velocity (ilist,3)  = buf_move_real(8*k-2)
      mass     (ilist  )  = buf_move_real(8*k-1)
      charge   (ilist  )  = buf_move_real(8*k  )
      atmcls   (ilist  )  = buf_move_int (3*k-1)
      id_l2g   (ilist  )  = buf_move_int (3*k-0)
      ig = id_l2g(ilist)
      id_g2l(ig) = ilist
    end do

    do k = charge_move_domain+1, charge_move_domain+nocharge_move_domain
      i = buf_move_int(3*k-2)
      start_i = start_atom(i)
      nocharge_comm1(i) = nocharge_comm1(i) + 1
      ilist = nocharge_comm1(i) + nocharge_stay(i) + ncharge(i) + start_i
      coord    (ilist,1)  = buf_move_real(8*k-7)
      coord    (ilist,2)  = buf_move_real(8*k-6)
      coord    (ilist,3)  = buf_move_real(8*k-5)
      velocity (ilist,1)  = buf_move_real(8*k-4)
      velocity (ilist,2)  = buf_move_real(8*k-3)
      velocity (ilist,3)  = buf_move_real(8*k-2)
      mass     (ilist  )  = buf_move_real(8*k-1)
      charge   (ilist  )  = buf_move_real(8*k  )
      atmcls   (ilist  )  = buf_move_int (3*k-1)
      id_l2g   (ilist  )  = buf_move_int (3*k-0)
      ig = id_l2g(ilist)
      id_g2l(ig) = ilist
    end do

    j = 0
    do ip = 1, num
      do k = 1, charge_comm(ip)
        j = j + 1
        i = cell_g2l(buf_comm_int(3*j-2))
        start_i = start_atom(i)
        charge_comm1(i) = charge_comm1(i) + 1
        ilist = charge_comm1(i) + charge_stay(i) + start_i
        coord    (ilist,1)  = buf_comm_real(8*j-7)
        coord    (ilist,2)  = buf_comm_real(8*j-6)
        coord    (ilist,3)  = buf_comm_real(8*j-5)
        velocity (ilist,1)  = buf_comm_real(8*j-4)
        velocity (ilist,2)  = buf_comm_real(8*j-3)
        velocity (ilist,3)  = buf_comm_real(8*j-2)
        mass     (ilist  )  = buf_comm_real(8*j-1)
        charge   (ilist  )  = buf_comm_real(8*j  )
        atmcls   (ilist  )  = buf_comm_int (3*j-1)
        id_l2g   (ilist  )  = buf_comm_int (3*j-0)
        ig = id_l2g(ilist)
        id_g2l(ig) = ilist
      end do
  
      do k = charge_comm(ip)+1,charge_comm(ip)+nocharge_comm(ip)
        j = j + 1 
        i = cell_g2l(buf_comm_int(3*j-2))
        start_i = start_atom(i)
        nocharge_comm1(i) = nocharge_comm1(i) + 1
        ilist = nocharge_comm1(i) + nocharge_stay(i) + ncharge(i) + start_i
        coord    (ilist,1)  = buf_comm_real(8*j-7)    
        coord    (ilist,2)  = buf_comm_real(8*j-6)    
        coord    (ilist,3)  = buf_comm_real(8*j-5)    
        velocity (ilist,1)  = buf_comm_real(8*j-4)
        velocity (ilist,2)  = buf_comm_real(8*j-3)
        velocity (ilist,3)  = buf_comm_real(8*j-2)
        mass     (ilist  )  = buf_comm_real(8*j-1)
        charge   (ilist  )  = buf_comm_real(8*j  )
        atmcls   (ilist  )  = buf_comm_int (3*j-1)
        id_l2g   (ilist  )  = buf_comm_int (3*j-0)
        ig = id_l2g(ilist)
        id_g2l(ig) = ilist
      end do
    end do

    return

  end subroutine update_incoming_ptl_martini

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_enefunc_bondsq
  !> @brief        Updates BOND terms for each cell in the potential energy function.
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_enefunc_bondsq(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    integer                  :: ix, icel_local, k, n_stay, ko
    integer                  :: i1, i2, icel1, icel2, ip, num_proc

    real(wp),        pointer :: fc(:), r0(:), natom(:)
    real(wip),       pointer :: buf_move_real(:), buf_stay_real(:)
    integer,         pointer :: ncel_local, nboundary
    integer,         pointer :: id_g2l(:)
    integer,         pointer :: atom_2_cell(:)
    integer,         pointer :: bondlist(:,:), bondkind(:)
    integer,         pointer :: cell_rank(:)
    integer,         pointer :: bond_move(:)
    integer,         pointer :: buf_move_int(:), buf_stay_int(:)


    ncel_local    => domain%num_cell_local
    nboundary     => domain%num_cell_boundary
    id_g2l        => domain%id_g2l
    atom_2_cell   => domain%atom_2_cell
    natom         => domain%num_atom_t0
    cell_rank     => domain%domain_cell_rank

    bondlist      => enefunc%bond_list
    bondkind      => enefunc%bond_kind
    fc            => enefunc%bond_force_const
    r0            => enefunc%bond_dist_min

    bond_move     => domain%type1_comm
    buf_move_int  => domain%buf_var0_comm_int
    buf_stay_int  => domain%buf_var0_stay_int
    buf_move_real => domain%buf_var0_comm_real
    buf_stay_real => domain%buf_var0_stay_real

    num_proc      = domain%num_comm_proc
    bond_move(1:num_proc) = 0

    n_stay = 0
    ko = 0

    do ix = 1, enefunc%num_bondsq_domain

      i1 = id_g2l(bondlist(1,ix))
      i2 = id_g2l(bondlist(2,ix)) 
      icel1 = atom_2_cell(i1)
      icel2 = atom_2_cell(i2)
      if (icel1 <= ncel_local .and. icel2 <= ncel_local) then
        icel_local = min(icel1,icel2)
      else
        if (natom(icel1) >= natom(icel2)) then
          icel_local = icel2
        else if (natom(icel2) > natom(icel1)) then
          icel_local = icel1
        end if
      end if

      if (icel_local > ncel_local) then

        ip = cell_rank(icel_local)
        bond_move(ip) = bond_move(ip) + 1
        k = bond_move(ip)
        ko = ko + 1
        buf_move_real(2*ko-1) = r0(ix)
        buf_move_real(2*ko  ) = fc(ix)
        buf_move_int (4*ko-3) = ip
        buf_move_int (4*ko-2) = bondlist(1,ix)
        buf_move_int (4*ko-1) = bondlist(2,ix)
        buf_move_int (4*ko  ) = bondkind(  ix)

      else

        n_stay = n_stay + 1
        buf_stay_real(2*n_stay-1) = r0(ix)
        buf_stay_real(2*n_stay  ) = fc(ix)
        buf_stay_int (3*n_stay-2) = bondlist(1,ix)
        buf_stay_int (3*n_stay-1) = bondlist(2,ix)
        buf_stay_int (3*n_stay  ) = bondkind(  ix)
         
      end if

    end do 

    enefunc%bond_n_stay1 = n_stay
    enefunc%bonds_comm_domain = ko

    return

  end subroutine update_outgoing_enefunc_bondsq

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_enefunc_bond
  !> @brief        Updates BOND terms for each cell in the potential energy function.
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
    integer                  :: i, ix, icel_local, n_stay, n_stay1
    integer                  :: i1, i2, icel1, icel2, num_proc, ip
    integer                  :: ko_real, ko_int, ko

    real(wp),        pointer :: fc(:), r0(:), natom(:)
    real(wip),       pointer :: buf_move_real(:), buf_stay_real(:)
    integer,         pointer :: ncel_local, nboundary
    integer,         pointer :: id_g2l(:)
    integer,         pointer :: atom_2_cell(:)
    integer,         pointer :: bondlist(:,:), bondkind(:)
    integer,         pointer :: cell_rank(:)
    integer,         pointer :: bond_move(:)
    integer,         pointer :: buf_move_int(:), buf_stay_int(:)

    ncel_local    => domain%num_cell_local
    nboundary     => domain%num_cell_boundary
    id_g2l        => domain%id_g2l
    atom_2_cell   => domain%atom_2_cell
    natom         => domain%num_atom_t0
    cell_rank     => domain%domain_cell_rank

    bondlist      => enefunc%bond_list
    bondkind      => enefunc%bond_kind
    fc            => enefunc%bond_force_const
    r0            => enefunc%bond_dist_min

    bond_move     => domain%type2_comm
    buf_move_int  => domain%buf_var0_comm_int
    buf_stay_int  => domain%buf_var0_stay_int
    buf_move_real => domain%buf_var0_comm_real
    buf_stay_real => domain%buf_var0_stay_real

    num_proc      = domain%num_comm_proc
    bond_move(1:num_proc) = 0

    n_stay = 0
    ko_real = enefunc%bonds_comm_domain * 2
    ko_int  = enefunc%bonds_comm_domain * 4
    ko = 0

    do ix = enefunc%num_bondsq_domain+1,  &
            enefunc%num_bondsq_domain+enefunc%num_bond_domain

      i1    = id_g2l(bondlist(1,ix))
      i2    = id_g2l(bondlist(2,ix))
      i     = bondlist(3,ix)
      icel1 = atom_2_cell(i1)
      icel2 = atom_2_cell(i2)
      if (icel1 <= ncel_local .and. icel2 <= ncel_local) then
        icel_local = min(icel1,icel2)
      else
        if (natom(icel1) >= natom(icel2)) then
          icel_local = icel2
        else if (natom(icel2) > natom(icel1)) then
          icel_local = icel1
        end if
      end if

      if (icel_local > ncel_local) then

        ip = cell_rank(icel_local)
        bond_move(ip) = bond_move(ip) + 1
        ko = ko + 1
        buf_move_real(ko_real+2*ko-1) = r0(ix)
        buf_move_real(ko_real+2*ko  ) = fc(ix)
        buf_move_int (ko_int +4*ko-3) = ip
        buf_move_int (ko_int +4*ko-2) = bondlist(1,ix)
        buf_move_int (ko_int +4*ko-1) = bondlist(2,ix)
        buf_move_int (ko_int +4*ko  ) = bondkind(  ix)

      else

        n_stay = n_stay + 1
        n_stay1 = n_stay + enefunc%bond_n_stay1
        buf_stay_real(2*n_stay1-1) = r0(ix)
        buf_stay_real(2*n_stay1  ) = fc(ix)
        buf_stay_int (3*n_stay1-2) = bondlist(1,ix)
        buf_stay_int (3*n_stay1-1) = bondlist(2,ix)
        buf_stay_int (3*n_stay1  ) = bondkind(  ix)
         
      end if

    end do 
    enefunc%bond_n_stay2 = n_stay
    enefunc%bondq_comm_domain = ko

    return

  end subroutine update_outgoing_enefunc_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_enefunc_bond
  !> @brief        Updates BOND terms for each cell in the potential energy function.
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
    integer                  :: i, k, kx, n_stay, n_stay1, ix, ip, num
    integer                  :: list, lists, listq

    real(wip),       pointer :: buf_move_real(:), buf_stay_real(:)
    real(wp),        pointer :: fc(:), r0(:)
    integer,         pointer :: ncel_local, nboundary, id_g2l(:)
    integer,         pointer :: bondsq_move(:), bond_move(:)
    integer,         pointer :: buf_move_int(:), buf_stay_int(:)
    integer,         pointer :: bondlist(:,:), bondkind(:)

    ncel_local    => domain%num_cell_local
    nboundary     => domain%num_cell_boundary
    id_g2l        => domain%id_g2l

    bondlist      => enefunc%bond_list
    bondkind      => enefunc%bond_kind
    fc            => enefunc%bond_force_const
    r0            => enefunc%bond_dist_min

    bondsq_move   => domain%type1_comm
    bond_move     => domain%type2_comm
    buf_move_int  => domain%buf_var0_comm_int
    buf_stay_int  => domain%buf_var0_stay_int
    buf_move_real => domain%buf_var0_comm_real
    buf_stay_real => domain%buf_var0_stay_real

    num = domain%num_comm_proc

    k = enefunc%bond_n_stay1
    do ip = 1, domain%num_comm_proc
      k = k + bondsq_move(ip)
    end do
    enefunc%num_bondsq_domain = k

    k = enefunc%bond_n_stay2
    do ip = 1, domain%num_comm_proc
      k = k + bond_move(ip)
    end do
    enefunc%num_bond_domain = k

#ifdef DEBUG
    if (enefunc%num_bondsq_domain + enefunc%num_bond_domain > MaxBond) &
        call error_msg( &
            'Debug: update_incoming_enefunc_bond: Bond number exceeds MaxBond.')
#endif

    !$omp parallel do private(n_stay, i, ix)
    do n_stay = 1, enefunc%bond_n_stay1
      bondlist(1,n_stay) = buf_stay_int (3*n_stay-2)
      bondlist(2,n_stay) = buf_stay_int (3*n_stay-1)
      bondkind(  n_stay) = buf_stay_int (3*n_stay  )
      r0      (  n_stay) = buf_stay_real(2*n_stay-1)
      fc      (  n_stay) = buf_stay_real(2*n_stay  )
    end do
    !$omp end parallel do
    !$omp parallel do private(n_stay, k, n_stay1)
    do n_stay = 1, enefunc%bond_n_stay2
      k = enefunc%num_bondsq_domain+n_stay
      n_stay1 = enefunc%bond_n_stay1 + n_stay
      bondlist(1,k) = buf_stay_int (3*n_stay1-2)
      bondlist(2,k) = buf_stay_int (3*n_stay1-1)
      bondkind(  k) = buf_stay_int (3*n_stay1  )
      r0      (  k) = buf_stay_real(2*n_stay1-1)
      fc      (  k) = buf_stay_real(2*n_stay1  )
    end do
    !$omp end parallel do

    n_stay = enefunc%bond_n_stay1
    n_stay1 = enefunc%bond_n_stay2
    list = 0
    lists = 0
    listq = 0
    do ip = 1, num
      do k = 1, bondsq_move(ip)
        list = list + 1
        lists = lists + 1
        kx = lists + n_stay
        bondlist(1,kx) = buf_move_int (3*list-2)
        bondlist(2,kx) = buf_move_int (3*list-1)
        bondkind(  kx) = buf_move_int (3*list  )
        r0      (  kx) = buf_move_real(2*list-1)
        fc      (  kx) = buf_move_real(2*list  )
      end do
      do k = 1, bond_move(ip)
        list = list + 1
        listq = listq + 1
        kx = enefunc%num_bondsq_domain + n_stay1 + listq 
        bondlist(1,kx) = buf_move_int (3*list-2)
        bondlist(2,kx) = buf_move_int (3*list-1)
        bondkind(  kx) = buf_move_int (3*list  )
        r0      (  kx) = buf_move_real(2*list-1)
        fc      (  kx) = buf_move_real(2*list  )
      end do
    end do

    return

  end subroutine update_incoming_enefunc_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_enefunc_anglflex
  !> @brief        Updates flexible ANGLE terms.
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_enefunc_anglflex(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    integer                  :: ix, icel_local, n_stay, ko
    integer                  :: i1, i2, icel1, icel2, ip, num_proc

    real(wp),        pointer :: natom(:)
    integer,         pointer :: ncel_local, nboundary
    integer,         pointer :: id_g2l(:)
    integer,         pointer :: atom_2_cell(:)
    integer,         pointer :: angllist(:,:), anglkind(:)
    integer,         pointer :: cell_rank(:)
    integer,         pointer :: angl_move(:)
    integer,         pointer :: buf_move_int(:), buf_stay_int(:)

    ncel_local      => domain%num_cell_local
    nboundary       => domain%num_cell_boundary
    id_g2l          => domain%id_g2l
    atom_2_cell     => domain%atom_2_cell
    natom           => domain%num_atom_t0
    cell_rank       => domain%domain_cell_rank

    angllist        => enefunc%angl_list
    anglkind        => enefunc%angl_kind
    angl_move       => domain%type1_comm
    buf_move_int    => domain%buf_var0_comm_int
    buf_stay_int    => domain%buf_var0_stay_int

    num_proc      = domain%num_comm_proc
    angl_move(1:num_proc) = 0

    n_stay = 0
    ko = 0

    do ix = 1, enefunc%num_angle_flexible_domain

      i1    = id_g2l(angllist(1,ix))
      i2    = id_g2l(angllist(3,ix))
      icel1 = atom_2_cell(i1)
      icel2 = atom_2_cell(i2)
      if (icel1 <= ncel_local .and. icel2 <= ncel_local) then
        icel_local = min(icel1,icel2)
      else
        if (natom(icel1) >= natom(icel2)) then
          icel_local = icel2
        else if (natom(icel2) > natom(icel1)) then
          icel_local = icel1
        end if
      end if

      if (icel_local > ncel_local) then

        ip = cell_rank(icel_local)
        angl_move(ip) = angl_move(ip) + 1
        ko = ko + 1
        buf_move_int (5*ko-4) = ip
        buf_move_int (5*ko-3) = angllist(1,ix)
        buf_move_int (5*ko-2) = angllist(2,ix)
        buf_move_int (5*ko-1) = angllist(3,ix)
        buf_move_int (5*ko  ) = anglkind(  ix)

      else

        n_stay = n_stay + 1
        buf_stay_int(4*n_stay-3) = angllist(1,ix)
        buf_stay_int(4*n_stay-2) = angllist(2,ix)
        buf_stay_int(4*n_stay-1) = angllist(3,ix)
        buf_stay_int(4*n_stay  ) = anglkind(  ix)

      end if

    end do

    enefunc%angl_n_stay1 = n_stay
    enefunc%anglf_comm_domain = ko

    return

  end subroutine update_outgoing_enefunc_anglflex

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_enefunc_angllocal
  !> @brief        Updates local ANGLE terms.
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_enefunc_angllocal(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    integer                  :: ix, icel_local, ko, ko_real, ko_int
    integer                  :: n_stay_int, n_stay
    integer                  :: i1, i2, icel1, icel2, num_proc, ip

    real(wp),        pointer :: fc(:), r0(:), width(:), natom(:)
    real(wip),       pointer :: buf_move_real(:), buf_stay_real(:)
    integer,         pointer :: ncel_local, nboundary
    integer,         pointer :: id_g2l(:)
    integer,         pointer :: atom_2_cell(:)
    integer,         pointer :: cell_rank(:)
    integer,         pointer :: angllist(:,:)
    integer,         pointer :: angl_move(:)
    integer,         pointer :: buf_move_int(:), buf_stay_int(:)

    ncel_local       => domain%num_cell_local
    nboundary        => domain%num_cell_boundary
    id_g2l           => domain%id_g2l
    atom_2_cell      => domain%atom_2_cell
    natom            => domain%num_atom_t0
    cell_rank        => domain%domain_cell_rank

    angllist         => enefunc%angl_list
    fc               => enefunc%angl_force_const
    r0               => enefunc%angl_theta_min
    width            => enefunc%angl_width
    angl_move        => domain%type2_comm
    buf_move_int     => domain%buf_var0_comm_int
    buf_stay_int     => domain%buf_var0_stay_int
    buf_move_real    => domain%buf_var0_comm_real
    buf_stay_real    => domain%buf_var0_stay_real

    num_proc         = domain%num_comm_proc
    angl_move(1:num_proc) = 0

    n_stay_int  = enefunc%angl_n_stay1 * 4
    n_stay = 0
    ko_real = 0
    ko_int  = enefunc%anglf_comm_domain * 5
    ko = 0

    do ix = enefunc%num_angle_flexible_domain+1, &
            enefunc%num_angle_flexible_domain+enefunc%num_angle_local_domain

      i1    = id_g2l(angllist(1,ix))
      i2    = id_g2l(angllist(3,ix))
      icel1 = atom_2_cell(i1)
      icel2 = atom_2_cell(i2)
      if (icel1 <= ncel_local .and. icel2 <= ncel_local) then
        icel_local = min(icel1,icel2)
      else
        if (natom(icel1) >= natom(icel2)) then
          icel_local = icel2
        else if (natom(icel2) > natom(icel1)) then
          icel_local = icel1
        end if
      end if

      if (icel_local > ncel_local) then

        ip = cell_rank(icel_local)
        angl_move(ip) = angl_move(ip) + 1
        ko = ko + 1
        buf_move_int (ko_int +4*ko-3) = ip
        buf_move_int (ko_int +4*ko-2) = angllist(1,ix)
        buf_move_int (ko_int +4*ko-1) = angllist(2,ix)
        buf_move_int (ko_int +4*ko  ) = angllist(3,ix)
        buf_move_real(ko_real+3*ko-2) = fc      (  ix)
        buf_move_real(ko_real+3*ko-1) = r0      (  ix)
        buf_move_real(ko_real+3*ko  ) = width   (  ix)

      else

        n_stay = n_stay + 1
        buf_stay_int (n_stay_int+3*n_stay-2) = angllist(1,ix)
        buf_stay_int (n_stay_int+3*n_stay-1) = angllist(2,ix)
        buf_stay_int (n_stay_int+3*n_stay  ) = angllist(3,ix)
        buf_stay_real(           3*n_stay-2) = fc      (  ix)
        buf_stay_real(           3*n_stay-1) = r0      (  ix)
        buf_stay_real(           3*n_stay  ) = width   (  ix)

      end if

    end do

    enefunc%angl_n_stay2  = n_stay
    enefunc%angll_comm_domain = ko

    return

  end subroutine update_outgoing_enefunc_angllocal

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_enefunc_angl
  !> @brief        Updates ANGLE terms for each cell in the potential energy function.
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
    integer                  :: icel1, icel2, icel3, i1, i2, i3, ko
    integer                  :: ix, icel_local, num_proc, ip
    integer                  :: n_stay_int, n_stay_real, n_stay, ko_real, ko_int

    real(wp),        pointer :: fc(:), theta0(:), fc_ub(:), r0_ub(:)
    real(wip),       pointer :: buf_move_real(:), buf_stay_real(:)
    real(wp),        pointer :: natom(:)
    integer,         pointer :: ncel_local, nboundary
    integer,         pointer :: cell_g2l(:), cell_g2b(:)
    integer,         pointer :: id_g2l(:)
    integer,         pointer :: atom_2_cell(:)
    integer,         pointer :: angllist(:,:), anglkind(:)
    integer,         pointer :: cell_rank(:)
    integer,         pointer :: angl_move(:)
    integer,         pointer :: buf_move_int(:), buf_stay_int(:)

    ncel_local       => domain%num_cell_local
    nboundary        => domain%num_cell_boundary
    cell_g2l         => domain%cell_g2l
    cell_g2b         => domain%cell_g2b
    id_g2l           => domain%id_g2l
    atom_2_cell      => domain%atom_2_cell
    natom            => domain%num_atom_t0
    cell_rank        => domain%domain_cell_rank

    angllist         => enefunc%angl_list
    anglkind         => enefunc%angl_kind
    fc               => enefunc%angl_force_const
    theta0           => enefunc%angl_theta_min
    fc_ub            => enefunc%urey_force_const
    r0_ub            => enefunc%urey_rmin
    angl_move        => domain%type3_comm
    buf_move_int     => domain%buf_var0_comm_int
    buf_stay_int     => domain%buf_var0_stay_int
    buf_move_real    => domain%buf_var0_comm_real
    buf_stay_real    => domain%buf_var0_stay_real

    num_proc      = domain%num_comm_proc
    angl_move(1:num_proc) = 0

    n_stay_int  = enefunc%angl_n_stay1*4+enefunc%angl_n_stay2*3
    n_stay_real = enefunc%angl_n_stay2*3
    n_stay      = 0
    ko_real     = enefunc%angll_comm_domain*3
    ko_int      = enefunc%anglf_comm_domain*5+enefunc%angll_comm_domain*4
    ko = 0

    do ix = enefunc%num_angle_flexible_domain+enefunc%num_angle_local_domain+1,&
            enefunc%num_angle_flexible_domain+enefunc%num_angle_local_domain+enefunc%num_angle_domain

      i1    = id_g2l(angllist(1,ix))
      i2    = id_g2l(angllist(2,ix))
      i3    = id_g2l(angllist(3,ix))
      icel1 = atom_2_cell(i1)
      icel2 = atom_2_cell(i2)
      icel3 = atom_2_cell(i3)

      if (icel1 <= ncel_local .and. icel3 <= ncel_local) then
        icel_local = min(icel1,icel3)
      else
        if (natom(icel1) >= natom(icel3)) then
          icel_local = icel3
        else if (natom(icel3) > natom(icel1)) then
          icel_local = icel1
        end if
      end if

      if (icel_local > ncel_local) then

        ip = cell_rank(icel_local)
        angl_move(ip) = angl_move(ip) + 1
        ko = ko + 1
        buf_move_int (ko_int +5*ko-4) = ip
        buf_move_int (ko_int +5*ko-3) = angllist(1,ix)
        buf_move_int (ko_int +5*ko-2) = angllist(2,ix)
        buf_move_int (ko_int +5*ko-1) = angllist(3,ix)
        buf_move_int (ko_int +5*ko  ) = anglkind(  ix)
        buf_move_real(ko_real+4*ko-3) = fc      (  ix)
        buf_move_real(ko_real+4*ko-2) = theta0  (  ix)
        buf_move_real(ko_real+4*ko-1) = fc_ub   (  ix)
        buf_move_real(ko_real+4*ko  ) = r0_ub   (  ix)

      else

        n_stay = n_stay + 1
        buf_stay_int (n_stay_int +4*n_stay-3) = angllist(1,ix)
        buf_stay_int (n_stay_int +4*n_stay-2) = angllist(2,ix)
        buf_stay_int (n_stay_int +4*n_stay-1) = angllist(3,ix)
        buf_stay_int (n_stay_int +4*n_stay  ) = anglkind(  ix)
        buf_stay_real(n_stay_real+4*n_stay-3) = fc      (  ix)
        buf_stay_real(n_stay_real+4*n_stay-2) = theta0  (  ix)
        buf_stay_real(n_stay_real+4*n_stay-1) = fc_ub   (  ix)
        buf_stay_real(n_stay_real+4*n_stay  ) = r0_ub   (  ix)

      end if

    end do

    enefunc%angl_n_stay3  = n_stay
    enefunc%angl_comm_domain = ko

    return

  end subroutine update_outgoing_enefunc_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_enefunc_angl
  !> @brief        Updates ANGLE terms.
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
    integer                  :: i, k, n_stay, ix
    integer                  :: num1, num2, num3
    integer                  :: n_stay_int, n_stay_real

    integer,         pointer :: ncel_local, nboundary
    integer,         pointer :: id_g2l(:)
    integer,         pointer :: angllist(:,:), anglkind(:)
    integer,         pointer :: buf_move_int(:), buf_stay_int(:)
    real(wp),        pointer :: fc(:), r0(:), width(:), theta0(:)
    real(wp),        pointer :: fc_ub(:), r0_ub(:)
    real(wip),       pointer :: buf_move_real(:), buf_stay_real(:)

    ncel_local      => domain%num_cell_local
    nboundary       => domain%num_cell_boundary
    id_g2l          => domain%id_g2l

    angllist        => enefunc%angl_list
    anglkind        => enefunc%angl_kind
    fc              => enefunc%angl_force_const
    r0              => enefunc%angl_theta_min
    width           => enefunc%angl_width
    theta0          => enefunc%angl_theta_min
    fc_ub           => enefunc%urey_force_const
    r0_ub           => enefunc%urey_rmin
    buf_move_int    => domain%buf_var0_comm_int
    buf_stay_int    => domain%buf_var0_stay_int
    buf_move_real   => domain%buf_var0_comm_real
    buf_stay_real   => domain%buf_var0_stay_real

    k = enefunc%angl_n_stay1 + enefunc%angl_n_stay2 + enefunc%angl_n_stay3
    k = k + enefunc%anglf_comm_domain + enefunc%angll_comm_domain &
          + enefunc%angl_comm_domain
#ifdef DEBUG
    if (k > MaxAngl)  &
      call error_msg( &
            'Debug: update_incoming_enefunc_angl: Angle number exceeds MaxAngl.')
#endif
    ! flexible angle
    !
    enefunc%num_angle_flexible_domain = enefunc%angl_n_stay1 &
                                      + enefunc%anglf_comm_domain

    !$omp parallel do private(n_stay)
    do n_stay = 1, enefunc%angl_n_stay1
      angllist(1,n_stay) = buf_stay_int(4*n_stay-3)
      angllist(2,n_stay) = buf_stay_int(4*n_stay-2)
      angllist(3,n_stay) = buf_stay_int(4*n_stay-1)
      anglkind(  n_stay) = buf_stay_int(4*n_stay  )
    end do
    !$omp end parallel do

    num1 = enefunc%angl_n_stay1
    do i = 1, enefunc%anglf_comm_domain
      angllist(1,i+num1) = buf_move_int(4*i-3)
      angllist(2,i+num1) = buf_move_int(4*i-2)
      angllist(3,i+num1) = buf_move_int(4*i-1)
      anglkind(  i+num1) = buf_move_int(4*i  )
    end do

    ! local angle
    !
    enefunc%num_angle_local_domain = enefunc%angl_n_stay2 &
                                   + enefunc%angll_comm_domain
    n_stay_int  = 4*enefunc%angl_n_stay1
    num1 = enefunc%num_angle_flexible_domain
   
    !$omp parallel do private(n_stay, i, ix)
    do n_stay = 1, enefunc%angl_n_stay2
      angllist(1,n_stay+num1) = buf_stay_int (n_stay_int+3*n_stay-2)
      angllist(2,n_stay+num1) = buf_stay_int (n_stay_int+3*n_stay-1)
      angllist(3,n_stay+num1) = buf_stay_int (n_stay_int+3*n_stay  )
      fc      (  n_stay+num1) = buf_stay_real(           3*n_stay-2)
      r0      (  n_stay+num1) = buf_stay_real(           3*n_stay-1)
      width   (  n_stay+num1) = buf_stay_real(           3*n_stay  )
    end do
    !$omp end parallel do

    num1 = enefunc%num_angle_flexible_domain + enefunc%angl_n_stay2
    num2 = enefunc%anglf_comm_domain*4
    do i = 1, enefunc%angll_comm_domain
      angllist(1,i+num1) = buf_move_int (3*i-2+num2)
      angllist(2,i+num1) = buf_move_int (3*i-1+num2)
      angllist(3,i+num1) = buf_move_int (3*i  +num2)
      fc      (  i+num1) = buf_move_real(3*i-2     )
      r0      (  i+num1) = buf_move_real(3*i-1     )
      width   (  i+num1) = buf_move_real(3*i       )
    end do

    ! angle
    !
    enefunc%num_angle_domain = enefunc%angl_n_stay3 &
                             + enefunc%angl_comm_domain
    n_stay_int  = 4*enefunc%angl_n_stay1 + 3*enefunc%angl_n_stay2
    n_stay_real = 3*enefunc%angl_n_stay2
    num1 = enefunc%num_angle_flexible_domain + enefunc%num_angle_local_domain 

    !$omp parallel do private(n_stay, i, ix)
    do n_stay = 1, enefunc%angl_n_stay3
      angllist(1,n_stay+num1) = buf_stay_int (n_stay_int +4*n_stay-3)
      angllist(2,n_stay+num1) = buf_stay_int (n_stay_int +4*n_stay-2)
      angllist(3,n_stay+num1) = buf_stay_int (n_stay_int +4*n_stay-1)
      anglkind(  n_stay+num1) = buf_stay_int (n_stay_int +4*n_stay  )
      fc      (  n_stay+num1) = buf_stay_real(n_stay_real+4*n_stay-3)
      theta0  (  n_stay+num1) = buf_stay_real(n_stay_real+4*n_stay-2)
      fc_ub   (  n_stay+num1) = buf_stay_real(n_stay_real+4*n_stay-1)
      r0_ub   (  n_stay+num1) = buf_stay_real(n_stay_real+4*n_stay  )
    end do
    !$omp end parallel do

    num1 = enefunc%num_angle_flexible_domain + enefunc%num_angle_local_domain &
         + enefunc%angl_n_stay3
    num2 = enefunc%anglf_comm_domain*4 + enefunc%angll_comm_domain*3
    num3 = enefunc%angll_comm_domain*3
    do i = 1, enefunc%angl_comm_domain
      angllist(1,i+num1) = buf_move_int (4*i-3+num2)
      angllist(2,i+num1) = buf_move_int (4*i-2+num2)
      angllist(3,i+num1) = buf_move_int (4*i-1+num2)
      anglkind(  i+num1) = buf_move_int (4*i  +num2)
      fc      (  i+num1) = buf_move_real(4*i-3+num3)
      theta0  (  i+num1) = buf_move_real(4*i-2+num3)
      fc_ub   (  i+num1) = buf_move_real(4*i-1+num3)
      r0_ub   (  i+num1) = buf_move_real(4*i  +num3)
    end do

    return

  end subroutine update_incoming_enefunc_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_enefunc_diheflex
  !> @brief        Updates flexible DIHEDRAL terms.
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_enefunc_diheflex(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    integer                  :: i, ix, icel_local, n_stay, ko
    integer                  :: i1, i2, icel1, icel2, ip, num_proc

    real(wp),        pointer :: natom(:)
    integer,         pointer :: ncel_local, nboundary
    integer,         pointer :: id_g2l(:)
    integer,         pointer :: atom_2_cell(:)
    integer,         pointer :: dihelist(:,:), dihekind(:)
    integer,         pointer :: cell_rank(:)
    integer,         pointer :: dihefunc(:)
    integer,         pointer :: dihe_move(:)
    integer,         pointer :: buf_move_int(:), buf_stay_int(:)

    ncel_local      => domain%num_cell_local
    nboundary       => domain%num_cell_boundary
    id_g2l          => domain%id_g2l
    atom_2_cell     => domain%atom_2_cell
    natom           => domain%num_atom_t0
    cell_rank       => domain%domain_cell_rank

    dihelist        => enefunc%dihe_list
    dihekind        => enefunc%dihe_kind
    dihefunc        => enefunc%dihe_periodicity
    dihe_move       => domain%type1_comm
    buf_move_int    => domain%buf_var0_comm_int
    buf_stay_int    => domain%buf_var0_stay_int

    num_proc      = domain%num_comm_proc
    dihe_move(1:num_proc) = 0

    n_stay = 0
    ko = 0

    do ix = 1, enefunc%num_dihe_flexible_domain

      i1    = id_g2l(dihelist(1,ix))
      i2    = id_g2l(dihelist(4,ix))
      i     = dihelist(5,ix)
      icel1 = atom_2_cell(i1)
      icel2 = atom_2_cell(i2)
      if (icel1 <= ncel_local .and. icel2 <= ncel_local) then
        icel_local = min(icel1,icel2)
      else
        if (natom(icel1) >= natom(icel2)) then
          icel_local = icel2
        else if (natom(icel2) > natom(icel1)) then
          icel_local = icel1
        end if
      end if

      if (icel_local > ncel_local) then

        ip = cell_rank(icel_local)
        dihe_move(ip) = dihe_move(ip) + 1
        ko = ko + 1
        buf_move_int (7*ko-6) = ip
        buf_move_int (7*ko-5) = dihelist(1,ix)
        buf_move_int (7*ko-4) = dihelist(2,ix)
        buf_move_int (7*ko-3) = dihelist(3,ix)
        buf_move_int (7*ko-2) = dihelist(4,ix)
        buf_move_int (7*ko-1) = dihekind(  ix)
        buf_move_int (7*ko  ) = dihefunc(  ix)

      else

        n_stay = n_stay + 1
        buf_stay_int(6*n_stay-5) = dihelist(1,ix)
        buf_stay_int(6*n_stay-4) = dihelist(2,ix)
        buf_stay_int(6*n_stay-3) = dihelist(3,ix)
        buf_stay_int(6*n_stay-2) = dihelist(4,ix)
        buf_stay_int(6*n_stay-1) = dihekind(  ix)
        buf_stay_int(6*n_stay  ) = dihefunc(  ix)

      end if

    end do

    enefunc%dihe_n_stay1 = n_stay
    enefunc%dihef_comm_domain = ko

    return

  end subroutine update_outgoing_enefunc_diheflex

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_enefunc_dihelocal
  !> @brief        Updates local DIHEDRAL terms.
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_enefunc_dihelocal(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    integer                  :: i, ix, icel_local, ko, ko_real, ko_int
    integer                  :: n_stay_int, n_stay
    integer                  :: i1, i2, icel1, icel2, num_proc, ip

    real(wp),        pointer :: fc(:), phase(:), width(:), natom(:)
    real(wip),       pointer :: buf_move_real(:), buf_stay_real(:)
    integer,         pointer :: ncel_local, nboundary
    integer,         pointer :: id_g2l(:)
    integer,         pointer :: atom_2_cell(:)
    integer,         pointer :: cell_rank(:)
    integer,         pointer :: dihelist(:,:), func(:)
    integer,         pointer :: dihe_move(:)
    integer,         pointer :: buf_move_int(:), buf_stay_int(:)

    ncel_local       => domain%num_cell_local
    nboundary        => domain%num_cell_boundary
    id_g2l           => domain%id_g2l
    atom_2_cell      => domain%atom_2_cell
    natom            => domain%num_atom_t0
    cell_rank        => domain%domain_cell_rank

    dihelist         => enefunc%dihe_list
    fc               => enefunc%dihe_force_const
    phase            => enefunc%dihe_phase
    width            => enefunc%dihe_width
    func             => enefunc%dihe_kind
    dihe_move        => domain%type2_comm
    buf_move_int     => domain%buf_var0_comm_int
    buf_stay_int     => domain%buf_var0_stay_int
    buf_move_real    => domain%buf_var0_comm_real
    buf_stay_real    => domain%buf_var0_stay_real

    num_proc         = domain%num_comm_proc
    dihe_move(1:num_proc) = 0

    n_stay_int  = enefunc%dihe_n_stay1 * 6
    n_stay = 0
    ko_real = 0
    ko_int  = enefunc%dihef_comm_domain * 7
    ko = 0

    do ix = enefunc%num_dihe_flexible_domain+1, &
            enefunc%num_dihe_flexible_domain+enefunc%num_dihe_local_domain

      i1    = id_g2l(dihelist(1,ix))
      i2    = id_g2l(dihelist(4,ix))
      i     = dihelist(5,ix)
      icel1 = atom_2_cell(i1)
      icel2 = atom_2_cell(i2)
      if (icel1 <= ncel_local .and. icel2 <= ncel_local) then
        icel_local = min(icel1,icel2)
      else
        if (natom(icel1) >= natom(icel2)) then
          icel_local = icel2
        else if (natom(icel2) > natom(icel1)) then
          icel_local = icel1
        end if
      end if

      if (icel_local > ncel_local) then

        ip = cell_rank(icel_local)
        dihe_move(ip) = dihe_move(ip) + 1
        ko = ko + 1
        buf_move_int (ko_int +6*ko-5) = ip
        buf_move_int (ko_int +6*ko-4) = dihelist(1,ix)
        buf_move_int (ko_int +6*ko-3) = dihelist(2,ix)
        buf_move_int (ko_int +6*ko-2) = dihelist(3,ix)
        buf_move_int (ko_int +6*ko-1) = dihelist(4,ix)
        buf_move_int (ko_int +6*ko  ) = func    (  ix)
        buf_move_real(        3*ko-2) = fc      (  ix)
        buf_move_real(        3*ko-1) = phase   (  ix)
        buf_move_real(        3*ko  ) = width   (  ix)

      else

        n_stay = n_stay + 1 
        buf_stay_int (n_stay_int+5*n_stay-4) = dihelist(1,ix)
        buf_stay_int (n_stay_int+5*n_stay-3) = dihelist(2,ix)
        buf_stay_int (n_stay_int+5*n_stay-2) = dihelist(3,ix)
        buf_stay_int (n_stay_int+5*n_stay-1) = dihelist(4,ix)
        buf_stay_int (n_stay_int+5*n_stay  ) = func    (  ix)
        buf_stay_real(           3*n_stay-2) = fc      (  ix)
        buf_stay_real(           3*n_stay-1) = phase   (  ix)
        buf_stay_real(           3*n_stay  ) = width   (  ix)

      end if

    end do

    enefunc%dihe_n_stay2  = n_stay
    enefunc%dihel_comm_domain = ko


    return

  end subroutine update_outgoing_enefunc_dihelocal

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_enefunc_dihe
  !> @brief        Updates DIHEDRAL terms for each cell in the potential energy function.
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
    integer                  :: icel1, icel2, i1, i2, ko
    integer                  :: ix, icel_local, num_proc, ip 
    integer                  :: n_stay_int, n_stay_real, n_stay, ko_real, ko_int

    real(wp),        pointer :: fc(:), phase(:), natom(:)
    real(wip),       pointer :: buf_move_real(:), buf_stay_real(:)
    integer,         pointer :: ncel_local, nboundary
    integer,         pointer :: cell_g2l(:), cell_g2b(:)
    integer,         pointer :: id_g2l(:)
    integer,         pointer :: atom_2_cell(:)
    integer,         pointer :: dihelist(:,:), dihekind(:)
    integer,         pointer :: cell_rank(:)
    integer,         pointer :: nperiod(:)
    integer,         pointer :: dihe_move(:)
    integer,         pointer :: buf_move_int(:), buf_stay_int(:)

    ncel_local       => domain%num_cell_local
    nboundary        => domain%num_cell_boundary
    cell_g2l         => domain%cell_g2l
    cell_g2b         => domain%cell_g2b
    id_g2l           => domain%id_g2l
    atom_2_cell      => domain%atom_2_cell
    natom            => domain%num_atom_t0
    cell_rank        => domain%domain_cell_rank

    dihelist         => enefunc%dihe_list
    dihekind         => enefunc%dihe_kind
    fc               => enefunc%dihe_force_const
    phase            => enefunc%dihe_phase
    nperiod          => enefunc%dihe_periodicity
    dihe_move        => domain%type3_comm
    buf_move_int     => domain%buf_var0_comm_int
    buf_stay_int     => domain%buf_var0_stay_int
    buf_move_real    => domain%buf_var0_comm_real
    buf_stay_real    => domain%buf_var0_stay_real

    num_proc      = domain%num_comm_proc
    dihe_move(1:num_proc) = 0

    n_stay_int  = enefunc%dihe_n_stay1*6 + enefunc%dihe_n_stay2*5
    n_stay_real = enefunc%dihe_n_stay2*3
    n_stay = 0
    ko_int      = enefunc%dihef_comm_domain*7 + enefunc%dihel_comm_domain*6
    ko_real     = enefunc%dihel_comm_domain*3
    ko          = 0

    do ix = enefunc%num_dihe_flexible_domain &
              + enefunc%num_dihe_local_domain+1, &
            enefunc%num_dihe_flexible_domain &
              + enefunc%num_dihe_local_domain+enefunc%num_dihe_domain

      i1    = id_g2l(dihelist(1,ix))
      i2    = id_g2l(dihelist(4,ix))
      icel1 = atom_2_cell(i1)
      icel2 = atom_2_cell(i2)
      if (icel1 <= ncel_local .and. icel2 <= ncel_local) then
        icel_local = min(icel1,icel2)
      else
        if (natom(icel1) >= natom(icel2)) then
          icel_local = icel2
        else if (natom(icel2) > natom(icel1)) then
          icel_local = icel1
        end if
      end if

      if (icel_local > ncel_local) then

        ip = cell_rank(icel_local)
        dihe_move(ip) = dihe_move(ip) + 1
        ko = ko + 1
        buf_move_int (ko_int +7*ko-6) = ip
        buf_move_int (ko_int +7*ko-5) = dihelist(1,ix)
        buf_move_int (ko_int +7*ko-4) = dihelist(2,ix)
        buf_move_int (ko_int +7*ko-3) = dihelist(3,ix)
        buf_move_int (ko_int +7*ko-2) = dihelist(4,ix)
        buf_move_int (ko_int +7*ko-1) = nperiod (  ix)
        buf_move_int (ko_int +7*ko  ) = dihekind(  ix)
        buf_move_real(ko_real+2*ko-1) = fc      (  ix)
        buf_move_real(ko_real+2*ko  ) = phase   (  ix)

      else

        n_stay = n_stay + 1
        buf_stay_int (n_stay_int +6*n_stay-5) = dihelist(1,ix)
        buf_stay_int (n_stay_int +6*n_stay-4) = dihelist(2,ix)
        buf_stay_int (n_stay_int +6*n_stay-3) = dihelist(3,ix)
        buf_stay_int (n_stay_int +6*n_stay-2) = dihelist(4,ix)
        buf_stay_int (n_stay_int +6*n_stay-1) = nperiod (  ix)
        buf_stay_int (n_stay_int +6*n_stay  ) = dihekind(  ix)
        buf_stay_real(n_stay_real+2*n_stay-1) = fc      (  ix)
        buf_stay_real(n_stay_real+2*n_stay  ) = phase   (  ix)

      end if

    end do

    enefunc%dihe_n_stay3  = n_stay
    enefunc%dihe_comm_domain = ko

    return

  end subroutine update_outgoing_enefunc_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_enefunc_dihe
  !> @brief        Updates DIHEDRAL terms.
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
    integer                  :: i, k, n_stay, ix
    integer                  :: num1, num2, num3
    integer                  :: n_stay_int, n_stay_real

    integer,         pointer :: ncel_local, nboundary
    integer,         pointer :: id_g2l(:)
    integer,         pointer :: dihelist(:,:), dihekind(:)
    integer,         pointer :: nperiod(:)
    integer,         pointer :: buf_move_int(:), buf_stay_int(:)
    real(wp),        pointer :: fc(:), phase(:), width(:)
    real(wip),       pointer :: buf_move_real(:), buf_stay_real(:)

    ncel_local       => domain%num_cell_local
    nboundary        => domain%num_cell_boundary
    id_g2l           => domain%id_g2l

    dihelist         => enefunc%dihe_list
    fc               => enefunc%dihe_force_const
    nperiod          => enefunc%dihe_periodicity
    phase            => enefunc%dihe_phase
    width            => enefunc%dihe_width
    dihekind         => enefunc%dihe_kind
    buf_move_int     => domain%buf_var0_comm_int
    buf_stay_int     => domain%buf_var0_stay_int
    buf_move_real    => domain%buf_var0_comm_real
    buf_stay_real    => domain%buf_var0_stay_real

    k = enefunc%dihe_n_stay1 + enefunc%dihe_n_stay2 + enefunc%dihe_n_stay3
    k = k + enefunc%dihef_comm_domain + enefunc%dihel_comm_domain &
          + enefunc%dihe_comm_domain

#ifdef DEBUG
    if (k > MaxDihe)  &
      call error_msg( &
            'Debug: update_incoming_enefunc_dihe: Dihedral number exceeds MaxDihe.')
#endif

    ! flexible dihedral
    !
    enefunc%num_dihe_flexible_domain = enefunc%dihe_n_stay1 &
                                     + enefunc%dihef_comm_domain

    !$omp parallel do private(n_stay, i, ix)
    do n_stay = 1, enefunc%dihe_n_stay1
      dihelist(1,n_stay) = buf_stay_int(6*n_stay-5)
      dihelist(2,n_stay) = buf_stay_int(6*n_stay-4)
      dihelist(3,n_stay) = buf_stay_int(6*n_stay-3)
      dihelist(4,n_stay) = buf_stay_int(6*n_stay-2)
      dihekind(  n_stay) = buf_stay_int(6*n_stay-1)
      nperiod (  n_stay) = buf_stay_int(6*n_stay  )
    end do
    !$omp end parallel do

    num1 = enefunc%dihe_n_stay1
    do i = 1, enefunc%dihef_comm_domain
      dihelist(1,i+num1) = buf_move_int(6*i-5) 
      dihelist(2,i+num1) = buf_move_int(6*i-4) 
      dihelist(3,i+num1) = buf_move_int(6*i-3) 
      dihelist(4,i+num1) = buf_move_int(6*i-2) 
      dihekind(  i+num1) = buf_move_int(6*i-1)
      nperiod (  i+num1) = buf_move_int(6*i  )
    end do

    ! local dihedral
    !
    enefunc%num_dihe_local_domain = enefunc%dihe_n_stay2 &
                                  + enefunc%dihel_comm_domain
    n_stay_int = 6*enefunc%dihe_n_stay1
    num1 = enefunc%num_dihe_flexible_domain

    !$omp parallel do private(n_stay, i, ix)
    do n_stay = 1, enefunc%dihe_n_stay2
      dihelist(1,n_stay+num1) = buf_stay_int (n_stay_int+5*n_stay-4) 
      dihelist(2,n_stay+num1) = buf_stay_int (n_stay_int+5*n_stay-3) 
      dihelist(3,n_stay+num1) = buf_stay_int (n_stay_int+5*n_stay-2) 
      dihelist(4,n_stay+num1) = buf_stay_int (n_stay_int+5*n_stay-1) 
      dihekind(  n_stay+num1) = buf_stay_int (n_stay_int+5*n_stay  )
      fc      (  n_stay+num1) = buf_stay_real(           3*n_stay-2)
      phase   (  n_stay+num1) = buf_stay_real(           3*n_stay-1)
      width   (  n_stay+num1) = buf_stay_real(           3*n_stay  )
    end do
    !$omp end parallel do

    num1 = enefunc%num_dihe_flexible_domain + enefunc%dihe_n_stay2
    num2 = enefunc%dihef_comm_domain*6
    do i = 1, enefunc%dihel_comm_domain
      dihelist(1,i+num1) = buf_move_int (5*i-4+num2)
      dihelist(2,i+num1) = buf_move_int (5*i-3+num2)
      dihelist(3,i+num1) = buf_move_int (5*i-2+num2)
      dihelist(4,i+num1) = buf_move_int (5*i-1+num2)
      dihekind(  i+num1) = buf_move_int (5*i  +num2)
      fc      (  i+num1) = buf_move_real(3*i-2     )
      phase   (  i+num1) = buf_move_real(3*i-1     )
      width   (  i+num1) = buf_move_real(3*i       )
    end do

    ! dihedral
    !
    enefunc%num_dihe_domain = enefunc%dihe_n_stay3 &
                            + enefunc%dihe_comm_domain
    n_stay_int  = 6*enefunc%dihe_n_stay1 + 5*enefunc%dihe_n_stay2
    n_stay_real = 3*enefunc%dihe_n_stay2
    num1 = enefunc%num_dihe_flexible_domain + enefunc%num_dihe_local_domain

    !$omp parallel do private(n_stay, i, ix)
    do n_stay = 1, enefunc%dihe_n_stay3
      dihelist(1,n_stay+num1) = buf_stay_int (n_stay_int +6*n_stay-5)
      dihelist(2,n_stay+num1) = buf_stay_int (n_stay_int +6*n_stay-4)
      dihelist(3,n_stay+num1) = buf_stay_int (n_stay_int +6*n_stay-3)
      dihelist(4,n_stay+num1) = buf_stay_int (n_stay_int +6*n_stay-2)
      nperiod (  n_stay+num1) = buf_stay_int (n_stay_int +6*n_stay-1)
      dihekind(  n_stay+num1) = buf_stay_int (n_stay_int +6*n_stay  )
      fc      (  n_stay+num1) = buf_stay_real(n_stay_real+2*n_stay-1)
      phase   (  n_stay+num1) = buf_stay_real(n_stay_real+2*n_stay  )
    end do
    !$omp end parallel do

    num1 = enefunc%num_dihe_flexible_domain + enefunc%num_dihe_local_domain &
         + enefunc%dihe_n_stay3
    num2 = enefunc%dihef_comm_domain*6 + enefunc%dihel_comm_domain*5
    num3 = enefunc%dihel_comm_domain*3

    do i = 1, enefunc%dihe_comm_domain
      dihelist(1,i+num1) = buf_move_int (6*i-5+num2)
      dihelist(2,i+num1) = buf_move_int (6*i-4+num2)
      dihelist(3,i+num1) = buf_move_int (6*i-3+num2)
      dihelist(4,i+num1) = buf_move_int (6*i-2+num2)
      nperiod (  i+num1) = buf_move_int (6*i-1+num2)
      dihekind(  i+num1) = buf_move_int (6*i  +num2)
      fc      (  i+num1) = buf_move_real(2*i-1+num3)
      phase   (  i+num1) = buf_move_real(2*i  +num3)
    end do

    return

  end subroutine update_incoming_enefunc_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_enefunc_stack
  !> @brief        Updates BASE Stack terms for each cell.
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_enefunc_stack(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: icel1, icel2, icel3, i1, i2, i3, ko
    integer                  :: i, icel_local, n_stay
    integer                  :: ip, num_proc

    real(wp),        pointer :: epsilon(:), sigma(:), theta(:), natom(:)
    real(wip),       pointer :: buf_move_real(:), buf_stay_real(:)
    integer,         pointer :: ncel_local, nboundary
    integer,         pointer :: cell_g2l(:), cell_g2b(:)
    integer,         pointer :: id_g2l(:)
    integer,         pointer :: atom_2_cell(:)
    integer,         pointer :: stacklist(:,:)
    integer,         pointer :: cell_rank(:)
    integer,         pointer :: stack_move(:)
    integer,         pointer :: buf_move_int(:), buf_stay_int(:)

    ncel_local       => domain%num_cell_local
    nboundary        => domain%num_cell_boundary
    cell_g2l         => domain%cell_g2l
    cell_g2b         => domain%cell_g2b
    id_g2l           => domain%id_g2l
    atom_2_cell      => domain%atom_2_cell
    natom            => domain%num_atom_t0
    cell_rank       => domain%domain_cell_rank

    stacklist        => enefunc%base_stack_list
    epsilon          => enefunc%base_stack_epsilon
    sigma            => enefunc%base_stack_sigma
    theta            => enefunc%base_stack_theta_bs

    stack_move       => domain%type1_comm
    buf_move_int     => domain%buf_var0_comm_int
    buf_stay_int     => domain%buf_var0_stay_int
    buf_move_real    => domain%buf_var0_comm_real
    buf_stay_real    => domain%buf_var0_stay_real

    num_proc      = domain%num_comm_proc
    stack_move(1:num_proc) = 0

    n_stay = 0
    ko = 0

    do i = 1, enefunc%num_stack_domain

      i1    = id_g2l(stacklist(1,i))
      i2    = id_g2l(stacklist(2,i))
      i3    = id_g2l(stacklist(3,i))
      icel1 = atom_2_cell(i1)
      icel2 = atom_2_cell(i2)
      icel3 = atom_2_cell(i3)

      if (icel2 <= ncel_local .and. icel3 <= ncel_local) then
        icel_local = min(icel2,icel3)
      else
        if (natom(icel2) >= natom(icel3)) then
          icel_local = icel3
        else if (natom(icel3) > natom(icel2)) then
          icel_local = icel2
        end if
      end if

      if (icel_local > ncel_local) then

        ip = cell_rank(icel_local)
        stack_move(ip) = stack_move(ip) + 1
        ko = ko + 1
        buf_move_int (4*ko-3) = ip
        buf_move_int (4*ko-2) = stacklist(1,i)
        buf_move_int (4*ko-1) = stacklist(2,i)
        buf_move_int (4*ko  ) = stacklist(3,i)
        buf_move_real(3*ko-2) = epsilon  (  i)
        buf_move_real(3*ko-1) = sigma    (  i)
        buf_move_real(3*ko  ) = theta    (  i)

      else

        n_stay = n_stay + 1
        buf_stay_int (3*n_stay-2) = stacklist(1,i)
        buf_stay_int (3*n_stay-1) = stacklist(2,i)
        buf_stay_int (3*n_stay  ) = stacklist(3,i)
        buf_stay_real(3*n_stay-2) = epsilon  (  i)
        buf_stay_real(3*n_stay-1) = sigma    (  i)
        buf_stay_real(3*n_stay  ) = theta    (  i)

      end if

    end do

    enefunc%stack_n_stay = n_stay
    enefunc%stack_comm_domain = ko

    return

  end subroutine update_outgoing_enefunc_stack

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_enefunc_stack
  !> @brief        Updates Base Stack terms for each cell.
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_incoming_enefunc_stack(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: i, ix, k, n_stay, num

    real(wp),        pointer :: epsilon(:), sigma(:), theta(:)
    real(wip),       pointer :: buf_move_real(:), buf_stay_real(:)
    integer,         pointer :: ncel_local, nboundary
    integer,         pointer :: cell_g2l(:), cell_g2b(:)
    integer,         pointer :: id_g2l(:)
    integer,         pointer :: stacklist(:,:)
    integer,         pointer :: stack_move(:)
    integer,         pointer :: buf_move_int(:), buf_stay_int(:)

    ncel_local       => domain%num_cell_local
    nboundary        => domain%num_cell_boundary
    cell_g2l         => domain%cell_g2l
    cell_g2b         => domain%cell_g2b
    id_g2l           => domain%id_g2l

    stacklist        => enefunc%base_stack_list
    epsilon          => enefunc%base_stack_epsilon
    sigma            => enefunc%base_stack_sigma  
    theta            => enefunc%base_stack_theta_bs

    stack_move       => domain%type1_comm
    buf_move_int     => domain%buf_var0_comm_int
    buf_stay_int     => domain%buf_var0_stay_int
    buf_move_real    => domain%buf_var0_comm_real
    buf_stay_real    => domain%buf_var0_stay_real

    k = enefunc%stack_n_stay + enefunc%stack_comm_domain
#ifdef DEBUG
    if (k > MaxStack) &
      call error_msg( &
        'Debug: update_incoming_enefunc_stack: Stack number exceeds MaxStack.')
#endif
    enefunc%num_stack_domain = k

    !$omp parallel do private(n_stay, i, ix)
    do n_stay = 1, enefunc%stack_n_stay
      stacklist(1,n_stay) = buf_stay_int (3*n_stay-2)
      stacklist(2,n_stay) = buf_stay_int (3*n_stay-1)
      stacklist(3,n_stay) = buf_stay_int (3*n_stay  )
      epsilon  (  n_stay) = buf_stay_real(3*n_stay-2)
      sigma    (  n_stay) = buf_stay_real(3*n_stay-1)
      theta    (  n_stay) = buf_stay_real(3*n_stay  )
    end do 
    !$omp end parallel do

    num = enefunc%stack_n_stay
    do i = 1, enefunc%stack_comm_domain
      stacklist(1,i+num) = buf_move_int (3*i-2)
      stacklist(2,i+num) = buf_move_int (3*i-1)
      stacklist(3,i+num) = buf_move_int (3*i  ) 
      epsilon  (  i+num) = buf_move_real(3*i-2)
      sigma    (  i+num) = buf_move_real(3*i-1)
      theta    (  i+num) = buf_move_real(3*i  )
    end do

    return

  end subroutine update_incoming_enefunc_stack

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_enefunc_pwmcos
  !> @brief        Updates PWMCOS terms for each cell.
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_enefunc_pwmcos(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: i1, i2, k1, ko, ip, num_proc
    integer                  :: i, icel_local, n_stay

    real(wp),        pointer :: r0(:,:)
    real(wp),        pointer :: theta1(:,:), theta2(:,:), theta3(:,:)
    real(wp),        pointer :: ene_A(:,:), ene_C(:,:), ene_G(:,:)
    real(wp),        pointer :: ene_T(:,:), gamma(:,:), eps(:,:)
    real(wip),       pointer :: buf_move_real(:), buf_stay_real(:)
    integer,         pointer :: ncel_local, nboundary
    integer,         pointer :: cell_g2l(:), cell_g2b(:)
    integer,         pointer :: id_g2l(:)
    integer,         pointer :: atom_2_cell(:)
    integer,         pointer :: cell_rank(:)
    integer,         pointer :: count(:)
    integer,         pointer :: id(:), id_N(:), id_C(:)
    integer,         pointer :: pwmcos_move(:)
    integer,         pointer :: buf_move_int(:), buf_stay_int(:)

    ncel_local       => domain%num_cell_local
    nboundary        => domain%num_cell_boundary
    cell_g2l         => domain%cell_g2l
    cell_g2b         => domain%cell_g2b
    id_g2l           => domain%id_g2l
    atom_2_cell      => domain%atom_2_cell
    cell_rank        => domain%domain_cell_rank

    count            => enefunc%pwmcos_count
    id               => enefunc%pwmcos_protein_id
    id_N             => enefunc%pwmcos_protein_id_N
    id_C             => enefunc%pwmcos_protein_id_C
    r0               => enefunc%pwmcos_r0
    theta1           => enefunc%pwmcos_theta1
    theta2           => enefunc%pwmcos_theta2
    theta3           => enefunc%pwmcos_theta3
    ene_A            => enefunc%pwmcos_ene_A 
    ene_C            => enefunc%pwmcos_ene_C 
    ene_G            => enefunc%pwmcos_ene_G 
    ene_T            => enefunc%pwmcos_ene_T 
    gamma            => enefunc%pwmcos_gamma 
    eps              => enefunc%pwmcos_eps   

    pwmcos_move      => domain%type1_comm
    buf_move_int     => domain%buf_var0_comm_int
    buf_stay_int     => domain%buf_var0_stay_int
    buf_move_real    => domain%buf_var0_comm_real
    buf_stay_real    => domain%buf_var0_stay_real

    num_proc = domain%num_comm_proc
    pwmcos_move(1:num_proc) = 0

    n_stay = 0
    ko = 0

    do i = 1, enefunc%num_pwmcos_domain

      i1 = id_g2l(id(i))
      icel_local = atom_2_cell(i1)

      if (icel_local > ncel_local) then
    
        ip = cell_rank(icel_local)  
        pwmcos_move(ip) = pwmcos_move(ip) + 1
        ko = ko + 1
        buf_move_int(5*ko-4) = ip
        buf_move_int(5*ko-3) = id   (i)  
        buf_move_int(5*ko-2) = id_N (i)  
        buf_move_int(5*ko-1) = id_C (i)  
        buf_move_int(5*ko  ) = count(i)

        do i2 = 1, count(i)
          k1 = 60*(ko-1) + 10*(i2-1)
          buf_move_real(k1+ 1) = r0    (i2,i)
          buf_move_real(k1+ 2) = theta1(i2,i)
          buf_move_real(k1+ 3) = theta2(i2,i)
          buf_move_real(k1+ 4) = theta3(i2,i)
          buf_move_real(k1+ 5) = ene_A (i2,i)
          buf_move_real(k1+ 6) = ene_C (i2,i)
          buf_move_real(k1+ 7) = ene_G (i2,i)
          buf_move_real(k1+ 8) = ene_T (i2,i)
          buf_move_real(k1+ 9) = gamma (i2,i)
          buf_move_real(k1+10) = eps   (i2,i)
        end do

      else

        n_stay = n_stay + 1
        buf_stay_int (4*n_stay-3) = id   (i)
        buf_stay_int (4*n_stay-2) = id_N (i)
        buf_stay_int (4*n_stay-1) = id_C (i)
        buf_stay_int (4*n_stay  ) = count(i)

        do i2 = 1, count(i)
          k1 = 60*(n_stay-1) + 10*(i2-1)
          buf_stay_real(k1+ 1) = r0    (i2,i)
          buf_stay_real(k1+ 2) = theta1(i2,i)
          buf_stay_real(k1+ 3) = theta2(i2,i)
          buf_stay_real(k1+ 4) = theta3(i2,i)
          buf_stay_real(k1+ 5) = ene_A (i2,i)
          buf_stay_real(k1+ 6) = ene_C (i2,i)
          buf_stay_real(k1+ 7) = ene_G (i2,i)
          buf_stay_real(k1+ 8) = ene_T (i2,i)
          buf_stay_real(k1+ 9) = gamma (i2,i)
          buf_stay_real(k1+10) = eps   (i2,i)
        end do

      end if

    end do
   
    enefunc%pwmcos_n_stay = n_stay
    enefunc%pwmcos_comm_domain = ko

    return

  end subroutine update_outgoing_enefunc_pwmcos

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_enefunc_pwmcos
  !> @brief        Updates PWMCOS terms for each cell.
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine update_incoming_enefunc_pwmcos(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: i2, k, k1, n_stay, num
    integer                  :: i

    real(wp),        pointer :: r0(:,:)
    real(wp),        pointer :: theta1(:,:), theta2(:,:), theta3(:,:)
    real(wp),        pointer :: ene_A(:,:), ene_C(:,:), ene_G(:,:)
    real(wp),        pointer :: ene_T(:,:), gamma(:,:), eps(:,:)
    real(wip),       pointer :: buf_move_real(:), buf_stay_real(:)
    integer,         pointer :: ncel_local
    integer,         pointer :: count(:)
    integer,         pointer :: id(:), id_N(:), id_C(:)
    integer,         pointer :: pwmcos_move(:)
    integer,         pointer :: buf_move_int(:), buf_stay_int(:)

    ncel_local       => domain%num_cell_local

    count            => enefunc%pwmcos_count
    id               => enefunc%pwmcos_protein_id
    id_N             => enefunc%pwmcos_protein_id_N
    id_C             => enefunc%pwmcos_protein_id_C
    r0               => enefunc%pwmcos_r0
    theta1           => enefunc%pwmcos_theta1
    theta2           => enefunc%pwmcos_theta2
    theta3           => enefunc%pwmcos_theta3
    ene_A            => enefunc%pwmcos_ene_A
    ene_C            => enefunc%pwmcos_ene_C
    ene_G            => enefunc%pwmcos_ene_G
    ene_T            => enefunc%pwmcos_ene_T
    gamma            => enefunc%pwmcos_gamma
    eps              => enefunc%pwmcos_eps

    pwmcos_move      => domain%type1_comm
    buf_move_int     => domain%buf_var0_comm_int
    buf_stay_int     => domain%buf_var0_stay_int
    buf_move_real    => domain%buf_var0_comm_real
    buf_stay_real    => domain%buf_var0_stay_real

    k = enefunc%pwmcos_n_stay + enefunc%pwmcos_comm_domain
#ifdef DEBUG
    if (k > MaxPwmCos) then
      call error_msg( &
       'Debug: update_incoming_enefunc_pwmcos: PWMCOS number exceeds MaxPwmCos.')
    end if
#endif
    enefunc%num_pwmcos_domain = k

    do n_stay = 1, enefunc%pwmcos_n_stay

      id   (n_stay) = buf_stay_int(4*n_stay-3)
      id_N (n_stay) = buf_stay_int(4*n_stay-2)
      id_C (n_stay) = buf_stay_int(4*n_stay-1)
      count(n_stay) = buf_stay_int(4*n_stay  )
      do i2 = 1, count(n_stay)
        k1 = 60*(n_stay-1) + 10*(i2-1)
        r0    (i2,n_stay) = buf_stay_real(k1+ 1)
        theta1(i2,n_stay) = buf_stay_real(k1+ 2)
        theta2(i2,n_stay) = buf_stay_real(k1+ 3)
        theta3(i2,n_stay) = buf_stay_real(k1+ 4)
        ene_A (i2,n_stay) = buf_stay_real(k1+ 5)
        ene_C (i2,n_stay) = buf_stay_real(k1+ 6)
        ene_G (i2,n_stay) = buf_stay_real(k1+ 7)
        ene_T (i2,n_stay) = buf_stay_real(k1+ 8)
        gamma (i2,n_stay) = buf_stay_real(k1+ 9)
        eps   (i2,n_stay) = buf_stay_real(k1+10)
      end do
    end do

    num = enefunc%pwmcos_n_stay

    do i = 1, enefunc%pwmcos_comm_domain
      id   (i+num) = buf_move_int(4*i-3)
      id_N (i+num) = buf_move_int(4*i-2)
      id_C (i+num) = buf_move_int(4*i-1)
      count(i+num) = buf_move_int(4*i  )
      do i2 = 1, count(i+num)
        k1 = 60*(i-1) + 10*(i2-1)
        r0    (i2,i+num) = buf_move_real(k1+ 1)
        theta1(i2,i+num) = buf_move_real(k1+ 2)
        theta2(i2,i+num) = buf_move_real(k1+ 3)
        theta3(i2,i+num) = buf_move_real(k1+ 4)
        ene_A (i2,i+num) = buf_move_real(k1+ 5)
        ene_C (i2,i+num) = buf_move_real(k1+ 6)
        ene_G (i2,i+num) = buf_move_real(k1+ 7)
        ene_T (i2,i+num) = buf_move_real(k1+ 8)
        gamma (i2,i+num) = buf_move_real(k1+ 9)
        eps   (i2,i+num) = buf_move_real(k1+10)
      end do
    end  do

    return

  end subroutine update_incoming_enefunc_pwmcos

!======1=========2=========3=========4=========5=========6=========7=========8
!
!  Subroutine    update_outgoing_enefunc_pwmcosns
!> @brief        Updates PWMCOSns terms for each cell.
!! @authors      JJ
!! @param[in]    domain  : domain information
!! @param[inout] enefunc : energy potential functions information
!
!======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_enefunc_pwmcosns(domain, enefunc)

! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

! local variables
    integer                  :: i1, i2, k1, ko, ip, num_proc
    integer                  :: i, icel_local, n_stay

    real(wp),        pointer :: r0(:,:)
    real(wp),        pointer :: theta1(:,:), theta2(:,:)
    real(wp),        pointer :: ene(:,:)
    real(wip),       pointer :: buf_move_real(:), buf_stay_real(:)
    integer,         pointer :: ncel_local, nboundary
    integer,         pointer :: cell_g2l(:), cell_g2b(:)
    integer,         pointer :: id_g2l(:)
    integer,         pointer :: atom_2_cell(:)
    integer,         pointer :: cell_rank(:)
    integer,         pointer :: count(:)
    integer,         pointer :: id(:), id_N(:), id_C(:)
    integer,         pointer :: specificity(:,:)
    integer,         pointer :: pwmcosns_move(:)
    integer,         pointer :: buf_move_int(:), buf_stay_int(:)

    ncel_local       => domain%num_cell_local
    nboundary        => domain%num_cell_boundary
    cell_g2l         => domain%cell_g2l
    cell_g2b         => domain%cell_g2b
    id_g2l           => domain%id_g2l
    atom_2_cell      => domain%atom_2_cell
    cell_rank        => domain%domain_cell_rank

    count            => enefunc%pwmcosns_count
    id               => enefunc%pwmcosns_protein_id
    id_N             => enefunc%pwmcosns_protein_id_N
    id_C             => enefunc%pwmcosns_protein_id_C
    r0               => enefunc%pwmcosns_r0
    theta1           => enefunc%pwmcosns_theta1
    theta2           => enefunc%pwmcosns_theta2
    ene              => enefunc%pwmcosns_ene
    specificity      => enefunc%pwmcosns_specificity

    pwmcosns_move    => domain%type1_comm
    buf_move_int     => domain%buf_var0_comm_int
    buf_stay_int     => domain%buf_var0_stay_int
    buf_move_real    => domain%buf_var0_comm_real
    buf_stay_real    => domain%buf_var0_stay_real

    num_proc = domain%num_comm_proc
    pwmcosns_move(1:num_proc) = 0

    n_stay = 0
    ko = 0

    do i = 1, enefunc%num_pwmcosns_domain

      i1 = id_g2l(id(i))
      icel_local = atom_2_cell(i1)

      if (icel_local > ncel_local) then

        ip = cell_rank(icel_local)
        pwmcosns_move(ip) = pwmcosns_move(ip) + 1
        ko = ko + 1
        buf_move_int(11*ko-10) = ip
        buf_move_int(11*ko- 9) = id   (i)
        buf_move_int(11*ko- 8) = id_N (i)
        buf_move_int(11*ko- 7) = id_C (i)
        buf_move_int(11*ko- 6) = count(i)

        do i2 = 1, count(i)
          k1 = 24*(ko-1) + 4*(i2-1)
          buf_move_int (11*ko-6+i2) = specificity(i2,i)
          buf_move_real(k1+1)       = r0         (i2,i)
          buf_move_real(k1+2)       = theta1     (i2,i)
          buf_move_real(k1+3)       = theta2     (i2,i)
          buf_move_real(k1+4)       = ene        (i2,i)
        end do

      else

        n_stay = n_stay + 1
        buf_stay_int (10*n_stay-9) = id   (i)
        buf_stay_int (10*n_stay-8) = id_N (i)
        buf_stay_int (10*n_stay-7) = id_C (i)
        buf_stay_int (10*n_stay-6) = count(i)

        do i2 = 1, count(i)
          k1 = 24*(n_stay-1) + 4*(i2-1)
          buf_stay_int (10*n_stay-6+i2) = specificity(i2,i)
          buf_stay_real(k1+1) = r0    (i2,i)
          buf_stay_real(k1+2) = theta1(i2,i)
          buf_stay_real(k1+3) = theta2(i2,i)
          buf_stay_real(k1+4) = ene   (i2,i)
        end do

      end if

    end do

    enefunc%pwmcosns_n_stay = n_stay
    enefunc%pwmcosns_comm_domain = ko

    return

  end subroutine update_outgoing_enefunc_pwmcosns

!======1=========2=========3=========4=========5=========6=========7=========8
!
!  Subroutine    update_incoming_enefunc_pwmcosns
!> @brief        Updates PWMCOSns terms for each cell.
!! @authors      JJ
!! @param[in]    domain  : domain information
!! @param[inout] enefunc : energy potential functions information
!
!======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_incoming_enefunc_pwmcosns(domain, enefunc)

! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

! local variables
    integer                  :: i2, k, k1, n_stay, num
    integer                  :: i

    real(wp),        pointer :: r0(:,:)
    real(wp),        pointer :: theta1(:,:), theta2(:,:)
    real(wp),        pointer :: ene(:,:)
    real(wip),       pointer :: buf_move_real(:), buf_stay_real(:)
    integer,         pointer :: ncel_local
    integer,         pointer :: specificity(:,:)
    integer,         pointer :: count(:)
    integer,         pointer :: id(:), id_N(:), id_C(:)
    integer,         pointer :: pwmcosns_move(:)
    integer,         pointer :: buf_move_int(:), buf_stay_int(:)

    ncel_local       => domain%num_cell_local

    count            => enefunc%pwmcosns_count
    id               => enefunc%pwmcosns_protein_id
    id_N             => enefunc%pwmcosns_protein_id_N
    id_C             => enefunc%pwmcosns_protein_id_C
    r0               => enefunc%pwmcosns_r0
    theta1           => enefunc%pwmcosns_theta1
    theta2           => enefunc%pwmcosns_theta2
    ene              => enefunc%pwmcosns_ene
    specificity      => enefunc%pwmcosns_specificity

    pwmcosns_move    => domain%type1_comm
    buf_move_int     => domain%buf_var0_comm_int
    buf_stay_int     => domain%buf_var0_stay_int
    buf_move_real    => domain%buf_var0_comm_real
    buf_stay_real    => domain%buf_var0_stay_real

    k = enefunc%pwmcosns_n_stay + enefunc%pwmcosns_comm_domain
#ifdef DEBUG
    if (k > MaxPwmCosns) then
      call error_msg( &
          'Debug: update_incoming_enefunc_pwmcosns: PWMCOSNS number exceeds MaxPwmCosns.')
    end if
#endif
    enefunc%num_pwmcosns_domain = k     

    do n_stay = 1, enefunc%pwmcosns_n_stay

      id   (n_stay) = buf_stay_int(10*n_stay-9)
      id_N (n_stay) = buf_stay_int(10*n_stay-8)
      id_C (n_stay) = buf_stay_int(10*n_stay-7)
      count(n_stay) = buf_stay_int(10*n_stay-6)
      do i2 = 1, count(n_stay)
        k1 = 24*(n_stay-1) + 4*(i2-1)
        specificity(i2,n_stay) = buf_stay_int(10*n_stay-6+i2)
        r0         (i2,n_stay) = buf_stay_real(k1+1)
        theta1     (i2,n_stay) = buf_stay_real(k1+2)
        theta2     (i2,n_stay) = buf_stay_real(k1+3)
        ene        (i2,n_stay) = buf_stay_real(k1+4)
      end do
    end do

    num = enefunc%pwmcosns_n_stay

    do i = 1, enefunc%pwmcosns_comm_domain
      id   (i+num) = buf_move_int(10*i-9)
      id_N (i+num) = buf_move_int(10*i-8)
      id_C (i+num) = buf_move_int(10*i-7)
      count(i+num) = buf_move_int(10*i-6)
      do i2 = 1, count(i+num)
        specificity(i2,i+num) = buf_move_int(10*i-6+i2)
        k1 = 24*(i-1) + 4*(i2-1)
        r0    (i2,i+num) = buf_move_real(k1+1)
        theta1(i2,i+num) = buf_move_real(k1+2)
        theta2(i2,i+num) = buf_move_real(k1+3)
        ene   (i2,i+num) = buf_move_real(k1+4)
      end do

    end do

    return

  end subroutine update_incoming_enefunc_pwmcosns

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_enefunc_contact
  !> @brief        Updates CONTACT terms for each cell in the potential energy function.
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_enefunc_contact(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    integer                  :: i, icel_local, k, j, ko, ip, num_contact
    integer                  :: icel1, icel2
    integer                  :: list1, list2
    integer                  :: num

    real(wp),        pointer :: lj12(:), lj10(:), lj6(:), natom(:)
    real(wip),       pointer :: buf_move_real(:)
    integer,         pointer :: ncel_local, nboundary
    integer,         pointer :: id_g2l(:)
    integer,         pointer :: atom_2_cell(:)
    integer,         pointer :: num_proc(:), proc_list(:,:)
    integer,         pointer :: contactlist(:,:)
    integer,         pointer :: contactfunc(:)
    integer,         pointer :: contact_move(:), buf_move_int(:)

    ncel_local    => domain%num_cell_local
    nboundary     => domain%num_cell_boundary
    id_g2l        => domain%id_g2l
    atom_2_cell   => domain%atom_2_cell
    natom         => domain%num_atom_t0
    num_proc      => domain%num_proc
    proc_list     => domain%proc_list

    contactlist   => enefunc%contact_list
    contactfunc   => enefunc%contact_func
    lj12          => enefunc%contact_lj12    
    lj10          => enefunc%contact_lj10   
    lj6           => enefunc%contact_lj6    

    contact_move  => domain%type1_comm
    buf_move_int  => domain%buf_var0_comm_int
    buf_move_real => domain%buf_var0_comm_real

    num_contact = enefunc%num_contact_domain + enefunc%num_contact_boundary
    num = domain%num_comm_proc

    k = 0
    ko = 0
    contact_move(1:num) = 0

    do i = 1, num_contact
      list1 = contactlist(1,i)
      list2 = contactlist(2,i)
      list1 = id_g2l(list1)
      list2 = id_g2l(list2)
      icel_local = 0
      if (list1 > 0 .and. list2 > 0) then
        icel1 = atom_2_cell(list1)
        icel2 = atom_2_cell(list2)
        if (icel1 <= ncel_local .and. icel2 <= ncel_local) then
          icel_local = icel1
        else
          if (natom(icel1) >= natom(icel2)) then
            icel_local = icel2
          else if (natom(icel2) > natom(icel1)) then
            icel_local = icel1
          end if
        end if
      end if
      if (icel_local <= ncel_local .and. icel_local > 0) then
        k = k + 1
        contactlist(1,k) = contactlist(1,i)
        contactlist(2,k) = contactlist(2,i)
        contactfunc(  k) = contactfunc(  i)
        lj12       (  k) = lj12       (  i)
        lj10       (  k) = lj10       (  i)
        lj6        (  k) = lj6        (  i)
        if (num_proc(icel_local) > 0) then
          do j = 1, num_proc(icel_local)
            ip = proc_list(j,icel_local)
            contact_move(ip) = contact_move(ip) + 1
            ko = ko + 1
            buf_move_int(4*ko-3)  = ip
            buf_move_int(4*ko-2)  = contactlist(1,k)
            buf_move_int(4*ko-1)  = contactlist(2,k)
            buf_move_int(4*ko  )  = contactfunc(  k)
            buf_move_real(3*ko-2) = lj12(k)
            buf_move_real(3*ko-1) = lj10(k)
            buf_move_real(3*ko  ) = lj6 (k)
          end do
        end if
      end if
    end do

    enefunc%num_contact_domain = k
    enefunc%num_contact_boundary = ko

    return

  end subroutine update_enefunc_contact

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_enefunc_restraint
  !> @brief        Updates restraint terms.
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_enefunc_restraint(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: icel1, i1
    integer                  :: ko, ko_real, ko_int, n_stay_int, n_stay_real
    integer                  :: i, n_stay
    integer                  :: ip, num_proc

    real(wp),        pointer :: restraint_coord(:,:), restraint_force(:,:)
    real(wp),        pointer :: fit_coord(:,:)
    real(wip),       pointer :: buf_move_real(:), buf_stay_real(:)
    integer,         pointer :: ncel_local, nboundary
    integer,         pointer :: cell_g2l(:), cell_g2b(:)
    integer,         pointer :: id_g2l(:)
    integer,         pointer :: atom_2_cell(:)
    integer,         pointer :: restraint_atom(:), fitting_atom(:)
    integer,         pointer :: cell_rank(:)
    integer,         pointer :: rest_move(:), fit_move(:)
    integer,         pointer :: buf_move_int(:), buf_stay_int(:)


    ncel_local       => domain%num_cell_local
    nboundary        => domain%num_cell_boundary
    cell_g2l         => domain%cell_g2l
    cell_g2b         => domain%cell_g2b
    id_g2l           => domain%id_g2l
    atom_2_cell      => domain%atom_2_cell
    cell_rank        => domain%domain_cell_rank

    restraint_atom   => enefunc%restraint_atom 
    restraint_coord  => enefunc%restraint_coord   
    restraint_force  => enefunc%restraint_force   
    fitting_atom     => enefunc%fitting_atom
    fit_coord        => enefunc%fit_coord

    rest_move        => domain%type1_comm
    fit_move         => domain%type2_comm
    buf_move_int     => domain%buf_var0_comm_int
    buf_stay_int     => domain%buf_var0_stay_int
    buf_move_real    => domain%buf_var0_comm_real
    buf_stay_real    => domain%buf_var0_stay_real

    num_proc      = domain%num_comm_proc
    rest_move(1:num_proc) = 0
    fit_move(1:num_proc) = 0

    n_stay = 0
    ko = 0

    do i = 1, enefunc%num_rest_domain

      i1    = id_g2l(restraint_atom(i))
      icel1 = atom_2_cell(i1)

      if (icel1 > ncel_local) then

        ip = cell_rank(icel1)
        rest_move(ip) = rest_move(ip) + 1
        ko = ko + 1
        buf_move_int (2*ko-1) = ip
        buf_move_int (2*ko  ) = restraint_atom(i)
        buf_move_real(7*ko-6) = restraint_coord(1,i)
        buf_move_real(7*ko-5) = restraint_coord(2,i)
        buf_move_real(7*ko-4) = restraint_coord(3,i)
        buf_move_real(7*ko-3) = restraint_force(1,i)
        buf_move_real(7*ko-2) = restraint_force(2,i)
        buf_move_real(7*ko-1) = restraint_force(3,i)
        buf_move_real(7*ko  ) = restraint_force(4,i)

      else if (icel1 > 0 .and. icel1 <= ncel_local) then

        n_stay = n_stay + 1
        buf_stay_int (n_stay    ) = restraint_atom(i)
        buf_stay_real(7*n_stay-6) = restraint_coord(1,i)
        buf_stay_real(7*n_stay-5) = restraint_coord(2,i)
        buf_stay_real(7*n_stay-4) = restraint_coord(3,i)
        buf_stay_real(7*n_stay-3) = restraint_force(1,i)
        buf_stay_real(7*n_stay-2) = restraint_force(2,i)
        buf_stay_real(7*n_stay-1) = restraint_force(3,i)
        buf_stay_real(7*n_stay  ) = restraint_force(4,i)

      end if

    end do

    enefunc%rest_n_stay = n_stay
    enefunc%rest_comm_domain = ko

    ko_real = ko * 6
    ko_int  = ko * 2
    ko      = 0
    n_stay_int  = n_stay 
    n_stay_real = n_stay * 6
    n_stay  = 0

    do i = 1, enefunc%num_fit_domain

      i1    = id_g2l(fitting_atom(i))
      icel1 = atom_2_cell(i1)

      if (icel1 > ncel_local) then

        ip = cell_rank(icel1)
        fit_move(ip) = fit_move(ip) + 1
        ko = ko + 1
        buf_move_int (ko_int +2*ko-1) = ip
        buf_move_int (ko_int +2*ko  ) = fitting_atom(i)
        buf_move_real(ko_real+3*ko-2) = fit_coord(1,i)
        buf_move_real(ko_real+3*ko-1) = fit_coord(2,i)
        buf_move_real(ko_real+3*ko  ) = fit_coord(3,i)
       
      else

        n_stay = n_stay + 1
        buf_stay_int (n_stay_int +  n_stay  ) = fitting_atom(i)
        buf_stay_real(n_stay_real+3*n_stay-2) = fit_coord(1,i) 
        buf_stay_real(n_stay_real+3*n_stay-1) = fit_coord(2,i) 
        buf_stay_real(n_stay_real+3*n_stay  ) = fit_coord(3,i) 

      end if

    end do

    enefunc%fit_n_stay = n_stay
    enefunc%fit_comm_domain = ko 

    return

  end subroutine update_outgoing_enefunc_restraint

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_enefunc_restraint
  !> @brief        Updates restraint terms.
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_incoming_enefunc_restraint(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: k
    integer                  :: ko_real, ko_int
    integer                  :: i, ix, n_stay
    integer                  :: num
    integer                  :: n_stay_real, n_stay_int

    real(wp),        pointer :: restraint_coord(:,:), restraint_force(:,:)
    real(wp),        pointer :: fit_coord(:,:)
    real(wip),       pointer :: buf_move_real(:), buf_stay_real(:)
    integer,         pointer :: ncel_local, nboundary
    integer,         pointer :: cell_g2l(:), cell_g2b(:)
    integer,         pointer :: id_g2l(:)
    integer,         pointer :: atom_2_cell(:)
    integer,         pointer :: restraint_atom(:), fitting_atom(:)
    integer,         pointer :: cell_rank(:)
    integer,         pointer :: rest_move(:), fit_move(:)
    integer,         pointer :: buf_move_int(:), buf_stay_int(:)

    ncel_local       => domain%num_cell_local
    nboundary        => domain%num_cell_boundary
    cell_g2l         => domain%cell_g2l
    cell_g2b         => domain%cell_g2b
    id_g2l           => domain%id_g2l
    atom_2_cell      => domain%atom_2_cell
    cell_rank        => domain%domain_cell_rank

    restraint_atom   => enefunc%restraint_atom
    restraint_coord  => enefunc%restraint_coord
    restraint_force  => enefunc%restraint_force
    fitting_atom     => enefunc%fitting_atom
    fit_coord        => enefunc%fit_coord

    rest_move        => domain%type1_comm
    fit_move         => domain%type2_comm
    buf_move_int     => domain%buf_var0_comm_int
    buf_stay_int     => domain%buf_var0_stay_int
    buf_move_real    => domain%buf_var0_comm_real
    buf_stay_real    => domain%buf_var0_stay_real

    k = enefunc%rest_n_stay + enefunc%rest_comm_domain
    enefunc%num_rest_domain = k
#ifdef DEBUG
    if (k > MaxRest) &
      call error_msg( &
        'Debug: update_incoming_enefunc_rest: Restraint number exceeds MaxRest.')
#endif
    k = enefunc%fit_n_stay + enefunc%fit_comm_domain
    enefunc%num_fit_domain = k
#ifdef DEBUG
    if (k > MaxFit) &
      call error_msg( &
        'Debug: update_incoming_enefunc_fit: Fit number exceeds MaxFit.')
#endif

    !$omp parallel do private(n_stay, i, ix)
    do n_stay = 1, enefunc%rest_n_stay
      restraint_atom (n_stay) = buf_stay_int (n_stay)
      restraint_coord(1,n_stay) = buf_stay_real(7*n_stay-6)
      restraint_coord(2,n_stay) = buf_stay_real(7*n_stay-5)
      restraint_coord(3,n_stay) = buf_stay_real(7*n_stay-4)
      restraint_force(1,n_stay) = buf_stay_real(7*n_stay-3)
      restraint_force(2,n_stay) = buf_stay_real(7*n_stay-2)
      restraint_force(3,n_stay) = buf_stay_real(7*n_stay-1)
      restraint_force(4,n_stay) = buf_stay_real(7*n_stay  )
    end do 
    !$omp end parallel do

    num = enefunc%rest_n_stay
    do i = 1, enefunc%rest_comm_domain
      restraint_atom (  i+num) = buf_move_int (i)
      restraint_coord(1,i+num) = buf_move_real(7*i-6)
      restraint_coord(2,i+num) = buf_move_real(7*i-5)
      restraint_coord(3,i+num) = buf_move_real(7*i-4)
      restraint_force(1,i+num) = buf_move_real(7*i-3)
      restraint_force(2,i+num) = buf_move_real(7*i-2)
      restraint_force(3,i+num) = buf_move_real(7*i-1)
      restraint_force(4,i+num) = buf_move_real(7*i  )
    end do

    n_stay_real = 6*n_stay
    n_stay_int  = n_stay

    !$omp parallel do private(n_stay, i, ix)
    do n_stay = 1, enefunc%fit_n_stay
      fitting_atom(n_stay) = buf_stay_int (n_stay+n_stay_int)
      fit_coord (1,n_stay) = buf_stay_real(n_stay_real+3*n_stay-2)
      fit_coord (2,n_stay) = buf_stay_real(n_stay_real+3*n_stay-1)
      fit_coord (3,n_stay) = buf_stay_real(n_stay_real+3*n_stay  )
    end do
    !$omp end parallel do

    ko_int  = enefunc%rest_comm_domain
    ko_real = enefunc%rest_comm_domain * 6
    num = enefunc%fit_comm_domain
    do i = 1, enefunc%fit_comm_domain
      fitting_atom(i+num) = buf_move_int(ko_int  +  i  )
      fit_coord (1,i+num) = buf_move_real(ko_real+3*i-2)
      fit_coord (2,i+num) = buf_move_real(ko_real+3*i-1)
      fit_coord (3,i+num) = buf_move_real(ko_real+3*i  )
    end do

    return

  end subroutine update_incoming_enefunc_restraint

end module cg_migration_mod
