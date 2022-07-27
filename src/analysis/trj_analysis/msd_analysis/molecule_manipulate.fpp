!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   molecule_manipulate_mod
!> @brief   subroutines to manipulate s_one_molecule structures
!! @authors Donatas Surblys (DS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

! in older gfortran versions move_alloc is not pure
! delete pure keyword in such cases to allow compilation
#ifdef MOVE_ALLOC_IS_NOT_PURE
#define pure
#endif /* MOVE_ALLOC_IS_NOT_PURE */

module molecule_manipulate_mod

  use one_molecule_str_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public :: unwrap_molecules
  public :: unwrap_according_to_previous
  public :: cluster_bonds

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    unwrap_molecules
  !> @brief        unwrap s_one_molecule in an array
  !! @authors      DS
  !! @param[inout] molecules : a 1d s_one_molecule array
  !! @param[in]    pbc_box   : a real array(3, 3) with pbc information
  !! @param[inout] coord     : coordinates for the whole system
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  pure subroutine unwrap_molecules(pbc_box, molecules, coord)

    real(wp),             dimension(3, 3), intent(in)              :: pbc_box
    type(s_one_molecule), dimension(:),    intent(inout)           :: molecules
    real(wp),             dimension(:, :), intent(inout)           :: coord

    integer :: i

    do i=1,size(molecules)
      call unwrap_one_molecule(pbc_box, molecules(i), coord)
      call align_substructures_in_one_molecule(pbc_box, molecules(i), coord)
    end do

    return

  end subroutine unwrap_molecules

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    unwrap_one_molecules
  !> @brief        unwrap single s_one_molecule by traversing bonds
  !! @authors      DS
  !! @param[inout] one_molecule : a s_one_molecule structure
  !! @param[in]    pbc_box      : a real array(3, 3) with pbc information
  !! @param[inout] coord        : coordinates for the whole system
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  pure subroutine unwrap_one_molecule(pbc_box, one_molecule, coord)

    real(wp), dimension(3, 3), intent(in)                      :: pbc_box
    type(s_one_molecule),      intent(inout)                   :: one_molecule
    real(wp), dimension(:, :), intent(inout)                   :: coord

    integer                                :: ibond, iatm1, iatm2

    integer,  dimension(:, :), allocatable :: unwrap_order


    call set_unwrap_order_from_bonds(one_molecule%bond_list, unwrap_order)

    do ibond = 1, size(unwrap_order, 2)

      iatm1 = unwrap_order(1, ibond)
      iatm2 = unwrap_order(2, ibond)

      ! wrap pos1 according to pos2

      call apply_pdb_condition_to_position(pbc_box, coord(:, iatm2), &
                                                    coord(:, iatm1))

    end do

    return

  end subroutine unwrap_one_molecule

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    align_substructures_in_one_moleucle
  !> @brief        align substructures inside s_one_molecule
  !! @authors      DS
  !! @param[in]    pbc_box      : a real array(3, 3) with pbc information
  !! @param[inout] one_molecule : a s_one_molecule structure
  !! @param[inout] coord        : coordinates for the whole system
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  pure subroutine align_substructures_in_one_molecule(pbc_box, one_molecule, &
    coord)

    real(wp), dimension(3, 3), intent(in)                      :: pbc_box
    type(s_one_molecule),      intent(inout)                   :: one_molecule
    real(wp), dimension(:, :), intent(inout)                   :: coord

    integer                                :: nmaset, imaset, iatm, idx
    real(wp), dimension(:, :), allocatable :: pos_centers, dpos_centers

    nmaset = one_molecule%num_sub
    if (nmaset < 2 ) then
      return
    end if

    allocate(pos_centers(3, nmaset), dpos_centers(3, nmaset))

    pos_centers = 0
    dpos_centers = 0

    do iatm = 1, one_molecule%num_atoms
      imaset = one_molecule%sub_idx(iatm)
      idx = one_molecule%atom_idx(iatm)
      pos_centers(:, imaset) = pos_centers(:, imaset) + coord(:, idx)
    end do

    ! use centroids as centers
    forall (imaset=1:nmaset)
      pos_centers(:, imaset) = pos_centers(:, imaset) &
        / count(one_molecule%sub_idx == imaset)
    end forall

    do imaset = 2, nmaset
      dpos_centers(:, imaset) = pos_centers(:, imaset) - pos_centers(:, 1)
      call apply_pdb_condition_to_displacement(pbc_box, dpos_centers(:, imaset))
    end do

    do iatm = 1, one_molecule%num_atoms
      imaset = one_molecule%sub_idx(iatm)
      idx = one_molecule%atom_idx(iatm)
      coord(:, idx) = coord(:, idx) &
        - pos_centers(:, imaset) + pos_centers(:, 1) + dpos_centers(:, imaset)
    end do

    return

  end subroutine align_substructures_in_one_molecule

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    unwrap_according_to_previous
  !> @brief        unwrap current coordinates according to previous coordinates
  !! @authors      DS
  !! @param[in]    pbc_box    : PBC array(3,3)
  !! @param[in]    coord_prev : coordinate array(3, :)
  !! @param[inout] coord      : coordinate array(3, :)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  pure subroutine unwrap_according_to_previous(pbc_box, coord_prev, coord)

    real(wp), dimension(3, 3), intent(in)    :: pbc_box
    real(wp), dimension(:, :), intent(in)    :: coord_prev
    real(wp), dimension(:, :), intent(inout) :: coord

    integer :: iatm

    do iatm = 1, size(coord, 2)
      call apply_pdb_condition_to_position(pbc_box, coord_prev(:, iatm), &
        coord(:, iatm))
    end do

    return

  end subroutine unwrap_according_to_previous

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    apply_pdb_condition_to_position
  !> @brief        apply periodic boundary condition to pos1
  !! @authors      DS
  !! @param[in]    pbc_box : a real array(3, 3) with pbc information
  !! @param[in]    pos2    : a real array(3) with particle position
  !! @param[inout] pos1    : a real array(3) with particle position
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  pure subroutine apply_pdb_condition_to_position(pbc_box, pos2, pos1)

    real(wp), dimension(3, 3), intent(in)    :: pbc_box
    real(wp), dimension(3),    intent(in)    :: pos2
    real(wp), dimension(3),    intent(inout) :: pos1

    real(wp), dimension(3)                   :: dpos


    dpos = pos1 - pos2
    call apply_pdb_condition_to_displacement(pbc_box, dpos)
    pos1 = dpos + pos2

    return

  end subroutine apply_pdb_condition_to_position

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    apply_pdb_condition_to_displacement
  !> @brief        apply periodic boundary condition to displacement
  !! @authors      DS
  !! @param[inout] dpos    : a real array(3) with particle displacement
  !! @param[in]    pbc_box : a real array(3, 3) with pbc information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  pure subroutine apply_pdb_condition_to_displacement(pbc_box, dpos)

    real(wp), dimension(3, 3), intent(in)    :: pbc_box
    real(wp), dimension(3),    intent(inout) :: dpos

    integer                                  :: iaxis

    forall (iaxis=1:3)
      dpos(iaxis) = dpos(iaxis) &
        - nint(dpos(iaxis)/pbc_box(iaxis, iaxis))*pbc_box(iaxis, iaxis)
    end forall

    return

  end subroutine apply_pdb_condition_to_displacement

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    set_unwrap_order_from_bonds
  !> @brief        set unwrap order for molecule
  !! @authors      DS
  !! @param[in]    bond_list   : a 2d integer array containing bond list
  !! @param[out]   unwrap_order: a 2d integer array containing unwrap order
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  pure subroutine set_unwrap_order_from_bonds(bond_list, unwrap_order)

    integer, dimension(:, :),              intent(in)   :: bond_list
    integer, dimension(:, :), allocatable, intent(out)  :: unwrap_order


    integer                                :: ibond, iatm1, iatm2, &
                                              bond_start, bond_end, &
                                              iclust, nclust, &
                                              cnt, idx

    logical, dimension(size(bond_list, 2)) :: bond_done
    integer, dimension(2, size(bond_list)) :: work
    logical, dimension(:),    allocatable  :: atom_done
    integer, dimension(:),    allocatable  :: bond_counts
    integer, dimension(:, :), allocatable  :: bonds


    ! do nothing if there are no bonds
    if (size(bond_list, 2) < 1) then

      allocate(unwrap_order(2,0))
      return

    end if

    call cluster_bonds(bond_list, bonds=bonds, bond_counts=bond_counts)

    bond_done = .false.

    allocate(atom_done(minval(bonds):maxval(bonds)))
    atom_done = .false.

    bond_end = 0

    cnt = 0
    nclust = size(bond_counts)
    do iclust = 1, nclust

      bond_start = bond_end + 1
      bond_end = bond_end + bond_counts(iclust)

      idx = bond_start
      cnt = cnt + 1
      work(:, cnt) = bonds(:, idx)
      bond_done(idx) = .true.

      atom_done(bonds(1, idx)) = .true.
      atom_done(bonds(2, idx)) = .true.

      do while (.not. all(bond_done(bond_start:bond_end)))

        do ibond = bond_start+1, bond_end

          if (bond_done(ibond)) cycle

          iatm1 = bonds(1, ibond)
          iatm2 = bonds(2, ibond)

          ! atom_done has iatm1
          if (atom_done(iatm1)) then
            ! atom_done doesn't have iatm2
            if (.not. atom_done(iatm2)) then
              cnt = cnt + 1
              work(:, cnt) = [ iatm2, iatm1 ]
              atom_done(iatm2) = .true.
            end if
            bond_done(ibond) = .true.
          ! atom_done has iatm2, but not iatm1
          else if (atom_done(iatm2)) then
            cnt = cnt + 1
            work(:, cnt) = [ iatm1, iatm2 ]
            atom_done(iatm1) = .true.
            bond_done(ibond) = .true.
          end if

        end do

      end do

    end do

    allocate(unwrap_order(2, cnt))
    unwrap_order(:, :) = work(1:2, 1:cnt)

    return

  end subroutine set_unwrap_order_from_bonds

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cluster_bonds
  !> @brief        clustering based on bond information
  !! @authors      DS
  !! @param[in]    bond_list   : integer array(2,:) containing bond information
  !! @param[out]   atoms       : array of atoms sorted by cluster
  !! @param[out]   atom_counts : count of atoms in each cluster
  !! @param[out]   bonds       : array of bonds sorted by cluster
  !! @param[out]   bond_counts : count of bonds in each cluster
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  pure subroutine cluster_bonds(bond_list, atoms, atom_counts, &
    bonds, bond_counts)

    integer, dimension(:, :),                        intent(in)  :: bond_list

    integer, dimension(:),    allocatable, optional, intent(out) :: atoms
    integer, dimension(:),    allocatable, optional, intent(out) :: atom_counts
    integer, dimension(:, :), allocatable, optional, intent(out) :: bonds
    integer, dimension(:),    allocatable, optional, intent(out) :: bond_counts

    integer, dimension(:),    allocatable :: work, &
                                             iatm2iclust, &
                                             atom_counts_tmp, bond_counts_tmp, &
                                             atoms_tmp
    integer, dimension(:, :), allocatable :: bonds_tmp

    integer                               :: maxatm, minatm, maxnatm, &
                                             iatm1, iatm2, &
                                             natm_iclust1, natm_iclust2, &
                                             nbond_iclust1, nbond_iclust2, &
                                             iclust, iclust1, iclust2, nclust, &
                                             ibond, nbond, &
                                             nbond_tmp, natm_tmp, &
                                             atom_start, atom_end, &
                                             bond_start, bond_end

    logical, dimension(:),    allocatable  :: empty_cluster

    type :: s_integer_array
      integer, dimension(:), allocatable :: i
    end type

    type(s_integer_array), dimension(:), allocatable &
                                         :: iclust2iatms, iclust2ibonds

    minatm = minval(bond_list)
    maxatm = maxval(bond_list)
    maxnatm = maxatm - minatm + 1
    nbond = size(bond_list, 2)

    allocate(iatm2iclust(minatm:maxatm), &
      iclust2iatms(minatm:maxatm), iclust2ibonds(minatm:maxatm), &
      empty_cluster(minatm:maxatm))

    do iclust = minatm, maxatm

      allocate(iclust2ibonds(iclust)%i(0))

    end do

    iatm2iclust = -1
    empty_cluster = .true.

    do ibond = 1, nbond

      iatm1 = bond_list(1, ibond)
      iatm2 = bond_list(2, ibond)

      if (.not.  allocated(iclust2iatms(iatm1)%i)) then
        allocate(iclust2iatms(iatm1)%i(1))
        iclust2iatms(iatm1)%i(1) = iatm1
        iatm2iclust(iatm1) = iatm1
        empty_cluster(iatm1) = .false.
      end if

      if (.not.  allocated(iclust2iatms(iatm2)%i)) then
        allocate(iclust2iatms(iatm2)%i(1))
        iclust2iatms(iatm2)%i(1) = iatm2
        iatm2iclust(iatm2) = iatm2
        empty_cluster(iatm2) = .false.
      end if


    end do

    do ibond = 1, nbond

      iatm1 = bond_list(1, ibond)
      iatm2 = bond_list(2, ibond)

      iclust1 = iatm2iclust(iatm1)
      iclust2 = iatm2iclust(iatm2)

      if (iclust1 /= iclust2) then

        ! transfer atoms
        natm_iclust1 = size(iclust2iatms(iclust1)%i)
        natm_iclust2 = size(iclust2iatms(iclust2)%i)

        iatm2iclust(iclust2iatms(iclust2)%i) = iclust1

        allocate(work(natm_iclust1+natm_iclust2))
        work(1:natm_iclust1) = iclust2iatms(iclust1)%i
        work(natm_iclust1+1:size(work)) = iclust2iatms(iclust2)%i

        call move_alloc(work, iclust2iatms(iclust1)%i)
        deallocate(iclust2iatms(iclust2)%i)


        ! transfer bonds
        nbond_iclust1 = size(iclust2ibonds(iclust1)%i)
        nbond_iclust2 = size(iclust2ibonds(iclust2)%i)

        allocate(work(nbond_iclust1+nbond_iclust2+1))
        work(1:nbond_iclust1) = iclust2ibonds(iclust1)%i
        work(nbond_iclust1+1:nbond_iclust1+nbond_iclust2) &
          = iclust2ibonds(iclust2)%i
        work(size(work)) = ibond

        call move_alloc(work, iclust2ibonds(iclust1)%i)
        deallocate(iclust2ibonds(iclust2)%i)

        empty_cluster(iclust2) = .true.

      else

        ! append 1 bond

        nbond_iclust1 = size(iclust2ibonds(iclust1)%i)

        allocate(work(nbond_iclust1+1))
        work(1:nbond_iclust1) = iclust2ibonds(iclust1)%i
        work(size(work)) = ibond

        call move_alloc(work, iclust2ibonds(iclust1)%i)

      end if

    end do


    allocate(atoms_tmp(maxnatm), atom_counts_tmp(maxnatm), &
      bonds_tmp(2, nbond), bond_counts_tmp(maxnatm))
    atom_end = 0
    bond_end = 0

    nclust = 0
    do iclust = minatm, maxatm


      if (empty_cluster(iclust)) cycle

      nclust = nclust + 1

      natm_tmp = size(iclust2iatms(iclust)%i)
      nbond_tmp = size(iclust2ibonds(iclust)%i)

      atom_counts_tmp(nclust) = natm_tmp
      bond_counts_tmp(nclust) = nbond_tmp

      atom_start = atom_end + 1
      atom_end = atom_end + natm_tmp

      bond_start = bond_end + 1
      bond_end = bond_end + nbond_tmp

      ! sort by atom number before assigning
      call quicksort(iclust2iatms(iclust)%i)
      atoms_tmp(atom_start:atom_end) = iclust2iatms(iclust)%i

      bonds_tmp(:, bond_start:bond_end) = bond_list(:, iclust2ibonds(iclust)%i)


    end do

    if (present(atoms)) then
      allocate(atoms(atom_end))
      atoms(:) = atoms_tmp(:atom_end)
    end if

    if (present(atom_counts)) then
      allocate(atom_counts(nclust))
      atom_counts(:) = atom_counts_tmp(1:nclust)
    end if

    if (present(bonds)) then
      allocate(bonds(2, nbond))
      bonds(:, :) = bonds_tmp(:, 1:nbond)
    end if

    if (present(bond_counts)) then
      allocate(bond_counts(nclust))
      bond_counts(:) = bond_counts_tmp(1:nclust)
    end if

    return

  end subroutine cluster_bonds

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    quicksort
  !> @brief        sort integer array using quicksort algorithm
  !! @authors      DS
  !! @param[inout] array : integer array
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  pure recursive subroutine quicksort(array)

    integer, dimension(:), intent(inout) :: array

    integer                              :: nelem, pivot, ileft, iright

    nelem = size(array)

    if (nelem < 2) return

    pivot = array(nelem/2)

    ileft = 1
    iright = nelem

    do

      do while (array(ileft) < pivot)
        ileft = ileft + 1
      end do

      do while (pivot < array(iright))
        iright = iright - 1
      end do

      if (iright <= ileft) exit

      array([ileft, iright]) = array([iright, ileft])
      ileft = ileft + 1
      iright = iright - 1

    end do

    call quicksort(array(:ileft-1))
    call quicksort(array(iright+1:))

    return

  end subroutine quicksort

end module molecule_manipulate_mod
