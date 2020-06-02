!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   nbond_list_mod
!> @brief   create N-bond list from 2-bond list
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module nbond_list_mod

  use messages_mod

  implicit none
  private

  ! subroutine
  public  :: create_nbond_list
  public  :: create_bond_map
  public  :: create_nbond_list_mol
  public  :: create_nbond_list_atom

  private :: search_nbond
  private :: add_nbond
  private :: bind_list
  private :: sort_list
  private :: r_side_greator

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    create_nbond_list
  !> @brief        create N-bond list from 2-bond list
  !! @authors      NT
  !! @param[in]    tbond_list : 2-bond list
  !! @param[in]    n          : level of N-bond
  !! @param[inout] nbond_list : N-bond list
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine create_nbond_list(tbond_list, n, nbond_list)

    ! formal arguments
    integer,                 intent(in)    :: tbond_list(:,:)
    integer,                 intent(in)    :: n
    integer,    allocatable, intent(inout) :: nbond_list(:,:)

    ! local variables
    integer,    allocatable  :: bond_map(:,:)


    ! create bond map
    !
    call create_bond_map(tbond_list, bond_map)


    ! create nbond list 
    !
    call create_nbond_list_mol(bond_map, n, nbond_list)


    deallocate(bond_map)

    return

  end subroutine create_nbond_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    create_bond_map
  !> @brief        create bond map
  !! @param[in]    tbond_list : 2-bond list
  !! @param[inout] bond_map   : status map of bond list
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine create_bond_map(tbond_list, bond_map)

    ! formal arguments
    integer,                 intent(in)    :: tbond_list(:,:)
    integer,    allocatable, intent(inout) :: bond_map(:,:)

    ! local variables
    integer                  :: i, i1, i2, ii, natom, nbond
    integer                  :: alloc_stat


    nbond = size(tbond_list(1,:))
    natom = 0
    do i = 1, nbond
      natom = max(natom, tbond_list(1,i))
      natom = max(natom, tbond_list(2,i))
    end do

    allocate(bond_map(5,natom), stat=alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    do i = 1, natom
      bond_map(1:5,i) = 0
    end do

    do i = 1, nbond

      i1 = tbond_list(1,i)
      i2 = tbond_list(2,i)

      do ii = 1, 5
        if (bond_map(ii,i1) == 0) &
          exit
      end do
      if (ii > 5) &
        call error_msg('Create_Bond_Map> Duplicate was overflow. > 5')
      bond_map(ii,i1) = i2

      do ii = 1, 5
        if (bond_map(ii,i2) == 0) &
          exit
      end do
      if (ii > 5) &
        call error_msg('Create_Bond_Map> Duplicate was overflow. > 5')
      bond_map(ii,i2) = i1

    end do

    return

  end subroutine create_bond_map

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    create_nbond_list_mol
  !> @brief        create molecule N-bond list from bond map
  !! @authors      NT
  !! @param[in]    bond_map   : status map of bond list
  !! @param[in]    n          : level of N-bond
  !! @param[inout] nbond_list : N-bond list
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine create_nbond_list_mol(bond_map, n, nbond_list)

    ! formal arguments
    integer,                 intent(in)    :: bond_map(:,:)
    integer,                 intent(in)    :: n
    integer,    allocatable, intent(inout) :: nbond_list(:,:)

    ! local variables
    integer                  :: i, natom, nlist
    integer,    allocatable  :: nbond(:)


    if (allocated(nbond_list)) &
      deallocate(nbond_list)


    ! make nbond-list about all atoms
    !
    natom = size(bond_map(1,:))

    nlist = 0

    allocate(nbond(n))

    do i = 1, natom
      call search_nbond(bond_map, i, 0, n, nbond, nbond_list, nlist)
    end do

    deallocate(nbond)


    ! bind nbond-list
    !
    call bind_list(nbond_list, nlist)

    return

  end subroutine create_nbond_list_mol

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    create_nbond_list_atom
  !> @brief        create atom N-bond list from bond map
  !! @note         N-bond list is not binded (list size was returned as nlist)
  !! @authors      NT
  !! @param[in]    bond_map   : status map of bond list
  !! @param[in]    n          : level of N-bond
  !! @param[in]    iatom      : atom index
  !! @param[inout] nbond_list : N-bond list
  !! @param[inout] nlist      : number of lists
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine create_nbond_list_atom(bond_map, n, iatom, nbond_list, nlist)

    ! formal arguments
    integer,                 intent(in)    :: bond_map(:,:)
    integer,                 intent(in)    :: n
    integer,                 intent(in)    :: iatom
    integer,    allocatable, intent(inout) :: nbond_list(:,:)
    integer,                 intent(inout) :: nlist

    ! local variables
    integer,    allocatable  :: nbond(:)


    if (allocated(nbond_list)) then
      if (n > size(nbond_list(:,1))) &
        call error_msg( &
        'Create_Nbond_List_Atom> n and size(nbond_list(:,1) is different.')
    end if


    ! make nbond-list about indicated atom
    !
    nlist = 0

    if (iatom > size(bond_map(1,:))) &
      return

    allocate(nbond(n))

    call search_nbond(bond_map, iatom, 0, n, nbond, nbond_list, nlist)

    deallocate(nbond)

    return

  end subroutine create_nbond_list_atom

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    search_nbond
  !> @brief        search N-bond around iatom
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine search_nbond(bond_map, iatom, ifrom, ilevel, nbond, &
                                    nbond_list, ilist)

    ! formal arguments
    integer,                 intent(in)    :: bond_map(:,:)
    integer,                 intent(in)    :: iatom
    integer,                 intent(in)    :: ifrom
    integer,                 intent(in)    :: ilevel
    integer,                 intent(inout) :: nbond(:)
    integer,    allocatable, intent(inout) :: nbond_list(:,:)
    integer,                 intent(inout) :: ilist

    ! local variables
    integer                  :: i, ito


    nbond(ilevel) = iatom

    if (ilevel == 1) then
      call add_nbond(nbond, nbond_list, ilist)
      return
    end if

    do i = 1, 5
      ito = bond_map(i,iatom)
      if (ito == 0) &
        exit
      if (ito == ifrom) &
        cycle
      call search_nbond(bond_map, ito, iatom, ilevel-1, &
                        nbond, nbond_list, ilist)
    end do

    return

  end subroutine search_nbond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    add_nbond
  !> @brief        add N-list to list
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine add_nbond(nbond, nbond_list, ilist)

    ! parameters
    integer,                 parameter     :: GrowSize = 1024

    ! formal arguments
    integer,                 intent(in)    :: nbond(:)
    integer,    allocatable, intent(inout) :: nbond_list(:,:)
    integer,                 intent(inout) :: ilist

    ! local variables
    integer                  :: i, nl, nt, alloc_stat
    integer,    allocatable  :: tmp(:,:)


    nl = size(nbond)

    if (allocated(nbond_list)) then
      nt = size(nbond_list(1,:))
    else
      nt = 0
    end if

    ilist = ilist + 1

    if (nt == 0) then

      allocate(nbond_list(nl,GrowSize),stat=alloc_stat)
      if (alloc_stat /= 0) &
        call error_msg_alloc

    else if (ilist > nt) then

      allocate(tmp(nl,nt),stat=alloc_stat)
      if (alloc_stat /= 0) &
        call error_msg_alloc

      do i = 1, nt
        tmp(1:nl,i) = nbond_list(1:nl,i)
      end do

      deallocate(nbond_list)

      allocate(nbond_list(nl,nt+GrowSize),stat=alloc_stat)
      if (alloc_stat /= 0) &
        call error_msg_alloc

      do i = 1, nt
        nbond_list(1:nl,i) = tmp(1:nl,i)
      end do

      nt = nt + GrowSize

      deallocate(tmp)

    end if

    nbond_list(1:nl,ilist) = nbond(1:nl)

    return

  end subroutine add_nbond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bind_list
  !> @brief        bind N-bond list
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine bind_list(nbond_list, nlist)

    ! formal arguments
    integer,    allocatable, intent(inout) :: nbond_list(:,:)
    integer,                 intent(inout) :: nlist

    ! local variables
    integer                  :: i, j, nl, nt, alloc_stat
    integer,    allocatable  :: nbond(:)
    integer,    allocatable  :: tmp(:,:)


    if (.not. allocated(nbond_list)) &
      return

    nl = size(nbond_list(:,1))
    nt = nlist


    ! unify direction
    !
    allocate(nbond(nl))

    do i = 1, nt
      if (nbond_list(1,i) > nbond_list(nl,i)) then
        do j = 1, nl
          nbond(j) = nbond_list(nl-j+1,i)
        end do
        nbond_list(1:nl,i) = nbond(1:nl)
      end if
    end do

    deallocate(nbond)


    ! sort less
    !
    call sort_list(nbond_list, nl, 1, nt)


    ! remove overlapped
    !
    allocate(tmp(nl,nt),stat=alloc_stat)

    nt = 0

    do i = 1, nlist

      if (nt /= 0) then
        do j = 1, nl
          if (tmp(j,nt) /= nbond_list(j,i)) &
            exit
        end do
        if (j > nl) &
          cycle
      end if

      nt = nt + 1
      tmp(1:nl,nt) = nbond_list(1:nl,i)

    end do


    ! trim array
    !
    deallocate(nbond_list)
    allocate(nbond_list(nl,nt),stat=alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    do i = 1, nt
      nbond_list(1:nl,i) = tmp(1:nl,i)
    end do

    deallocate(tmp)

    return

  end subroutine bind_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    sort_list
  !> @brief        sort N-bond list
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine sort_list(nbond_list, nl, start, end)

    ! formal arguments
    integer,                 intent(inout) :: nbond_list(:,:)
    integer,                 intent(in)    :: nl
    integer,                 intent(in)    :: start
    integer,                 intent(in)    :: end

    ! local variables
    integer,     allocatable :: x(:), t(:)
    integer                  :: i, j, ii


    allocate(x(nl), t(nl))

    x(:) = nbond_list(:,(start + end) / 2)
    i    = start
    j    = end

    do
      do while (r_side_greator(nl, nbond_list(:,i), x))
        i = i + 1
      end do

      do while (r_side_greator(nl, x, nbond_list(:,j)))
        j = j - 1
      end do

      if (i >= j) exit

      t(:) = nbond_list(:,i)
      nbond_list(:,i) = nbond_list(:,j)
      nbond_list(:,j) = t
      i = i + 1
      j = j - 1

    end do

    if (start < i - 1) &
      call sort_list(nbond_list, nl, start, i - 1)

    if (j + 1 < end) &
      call sort_list(nbond_list, nl, j + 1, end)

    deallocate(t, x)

    return

  end subroutine sort_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      r_side_greator
  !> @brief        check right side argument is grator or not
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function r_side_greator(nl, lv, rv)

    ! return value
    logical                  :: r_side_greator

    ! formal arguments
    integer,                 intent(in)    :: nl
    integer,                 intent(in)    :: lv(:)
    integer,                 intent(in)    :: rv(:)

    ! local varabiles
    integer                  :: i


    do i = 1, nl
      if (lv(i) > rv(i)) then
        r_side_greator = .false.
        return
      else if (lv(i) < rv(i)) then
        r_side_greator = .true.
        return
      end if
    end do

    r_side_greator = .false.

    return

  end function r_side_greator

end module nbond_list_mod
