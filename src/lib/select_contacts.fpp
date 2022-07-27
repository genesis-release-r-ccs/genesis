!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   select_contacts_mod
!> @brief   contact atom selection module
!! @authors Norio Takase (NT), Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module select_contacts_mod

  use molecules_str_mod
  use select_atoms_str_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type s_cell_linked_list

     ! origin      : simulation box origin
     ! region      : simulation box region
     ! cell_count  : cell count for each axis
     ! cell_length : edge length of each cell
     !
     ! head        : head[c] holds the index of the first atom in the c-th cell
     ! lscl        : An array implementation of the linked lists

     real(wp)                  :: origin(3)
     real(wp)                  :: region(3)
     integer                   :: cell_count(3)
     real(wp)                  :: cell_length(3)

     integer, allocatable      :: head(:)
     integer, allocatable      :: lscl(:)

     real(wp),         pointer :: coord(:,:)

  end type s_cell_linked_list

  type s_sdata_atom
     real(wp)                  :: sq_dis
     integer                   :: idx
  end type s_sdata_atom

  type s_sdata_mole
     real(wp)                  :: sq_dis
     integer                   :: idx1
     integer                   :: idx2
  end type s_sdata_mole

  ! parametes
  integer,  public,  parameter :: SelectModeAtom = 1
  integer,  public,  parameter :: SelectModeResi = 2
  integer,  public,  parameter :: SelectModeMole = 3
  integer,  private, parameter :: Empty = 0

  ! subroutines
  public  :: select_contacts
  private :: init_select
  private :: select_atom_contact
  private :: select_mole_contact
  private :: make_resi_mole_no
  private :: quicksort_atom
  private :: quicksort_mole
  private :: build_cell_linked_list
  private :: dealloc_cell_linked_list
  private :: search_neighbour_coord
  private :: compute_cell_by_coords
  private :: compute_cell_by_region

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    select_contacts
  !> @brief        select atoms in less than distance of atoms-A from atoms-B
  !! @authors      NT
  !! @param[in]    mode        : selection mode
  !! @param[in]    molecule_a  : array of atom coordinate A
  !! @param[in]    molecule_b  : array of atom coordinate B
  !! @param[in]    selatoms_a  : atom seleciton A
  !! @param[in]    selatoms_b  : atom seleciton B
  !! @param[in]    distance    : contact distance 
  !! @param[inout] selatoms    : result atom selection
  !! @param[in]    o_max_atom  : maximum atom count    (optional)
  !! @param[in]    o_sort_list : sort list by distance (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine select_contacts(mode,       &
                             molecule_a, &
                             molecule_b, &
                             selatoms_a, &
                             selatoms_b, &
                             distance,   &
                             selatoms,   &
                             o_max_atom, &
                             o_sort_list)

    ! formal arguments
    integer,                 intent(in)    :: mode
    type(s_molecule),        intent(in)    :: molecule_a
    type(s_molecule),        intent(in)    :: molecule_b
    type(s_selatoms),        intent(in)    :: selatoms_a
    type(s_selatoms),        intent(in)    :: selatoms_b
    real(wp),                intent(in)    :: distance
    type(s_selatoms),        intent(inout) :: selatoms
    integer,       optional, intent(in)    :: o_max_atom
    logical,       optional, intent(in)    :: o_sort_list
                                 
    ! local variables
    type(s_cell_linked_list) :: list_a, list_b
    integer                  :: max_atom, alloc_stat
    logical                  :: sort_list
    
    integer,     allocatable :: resi_mole_no(:)
    integer,     allocatable :: resi_mole_atom_no(:)

    character(*),  parameter :: ModeStr(3) = (/'Atom','Resi','Mole'/)

    
    if (present(o_max_atom)) then
      max_atom = o_max_atom
    else
      max_atom = 0
    end if

    if (present(o_sort_list)) then
      sort_list = o_sort_list
    else
      sort_list = .false.
    end if


    if (main_rank) then

      write(MsgOut, '(A,A4)') &
           'Select_Contacts> mode        : ', ModeStr(mode)
      write(MsgOut, '(A,I7)') &
           '                 sel atoms A : ', size(selatoms_a%idx)
      write(MsgOut, '(A,I7)') &
           '                 sel atoms B : ', size(selatoms_b%idx)
      write(MsgOut, '(A,F7.2)') &
           '                 distance    : ', distance
    end if

    ! initialize selection
    !
    call init_select(molecule_a, &
                     molecule_b, &
                     selatoms_a, &
                     selatoms_b, &
                     distance,   &
                     list_a,     &
                     list_b)

    select case(mode) 

    case (SelectModeAtom)

      call select_atom_contact(distance,   &
                               max_atom,   &
                               sort_list,  &
                               list_a,     &
                               list_b,     &
                               selatoms_b, &
                               selatoms)

    case (SelectModeResi)

      call make_resi_mole_no(molecule_b,   &
                             resi_mole_no, &
                             resi_mole_atom_no)

      call select_mole_contact(distance,   &
                               max_atom,   &
                               sort_list,  &
                               list_a,     &
                               list_b,     &
                               selatoms_b, &
                               resi_mole_no,      &
                               resi_mole_atom_no, &
                               selatoms)

      deallocate(resi_mole_no, &
                 resi_mole_atom_no, stat = alloc_stat)
      if (alloc_stat /= 0) &
        call error_msg_dealloc

    case (SelectModeMole)

      call select_mole_contact(distance,   &
                               max_atom,   &
                               sort_list,  &
                               list_a,     &
                               list_b,     &
                               selatoms_b, &
                               molecule_b%molecule_no,      &
                               molecule_b%molecule_atom_no, &
                               selatoms)

    end select

    if (main_rank) then
      write(MsgOut, '(A,I7)') &
           '                contact atom list : ', size(selatoms%idx)
      write(MsgOut, *) ''
    end if


    ! deallocate
    !
    deallocate(list_a%coord, list_b%coord, stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_dealloc

    call dealloc_cell_linked_list(list_a)
    call dealloc_cell_linked_list(list_b)

    return

  end subroutine select_contacts

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_select
  !> @brief        initialize cell-linked list
  !! @authors      NT
  !! @param[in]    molecule_a : array of atom coordinate A
  !! @param[in]    molecule_b : array of atom coordinate B
  !! @param[in]    selatoms_a : atom seleciton A
  !! @param[in]    selatoms_b : atom seleciton B
  !! @param[in]    distance   : contact distance 
  !! @param[inout] list_a     : cell-linked list A
  !! @param[inout] list_b     : cell-linked list B
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_select(molecule_a, &
                         molecule_b, &
                         selatoms_a, &
                         selatoms_b, &
                         distance,   &
                         list_a,     &
                         list_b)

    ! formal arguments
    type(s_molecule),         intent(in)    :: molecule_a
    type(s_molecule),         intent(in)    :: molecule_b
    type(s_selatoms),         intent(in)    :: selatoms_a
    type(s_selatoms),         intent(in)    :: selatoms_b
    real(wp),                 intent(in)    :: distance
    type(s_cell_linked_list), intent(inout) :: list_a
    type(s_cell_linked_list), intent(inout) :: list_b

    ! local variables
    real(wp)                  :: minbou(3), maxbou(3), mincl(3), region(3)
    integer                   :: i, alloc_stat
    integer                   :: num_a, num_b

    real(wp)                  :: cell_length(3)
    integer                   :: cell_count(3)

    real(wp),         pointer :: coord_a(:,:)
    real(wp),         pointer :: coord_b(:,:)


    ! setup coordinates
    !
    num_a = size(selatoms_a%idx)
    num_b = size(selatoms_b%idx)

    allocate(coord_a(3, num_a), &
             coord_b(3, num_b), &
             stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc


    minbou(1:3) =  1000000.0_wp
    maxbou(1:3) = -1000000.0_wp

    if (allocated(molecule_a%atom_refcoord)) then

      ! select from reference coordinates
      do i = 1, num_a
        coord_a(1:3, i) = molecule_a%atom_refcoord(1:3, selatoms_a%idx(i))
        minbou(1:3) = min(minbou(1:3), coord_a(1:3, i))
        maxbou(1:3) = max(maxbou(1:3), coord_a(1:3, i))
      end do
    
      do i = 1, num_b
        coord_b(1:3, i) = molecule_b%atom_refcoord(1:3, selatoms_b%idx(i))
        minbou(1:3) = min(minbou(1:3), coord_b(1:3, i))
        maxbou(1:3) = max(maxbou(1:3), coord_b(1:3, i))
      end do
    
    else
    
     ! select from input coordinates
     do i = 1, num_a
       coord_a(1:3, i) = molecule_a%atom_coord(1:3, selatoms_a%idx(i))
       minbou(1:3) = min(minbou(1:3), coord_a(1:3, i))
       maxbou(1:3) = max(maxbou(1:3), coord_a(1:3, i))
     end do
    
     do i = 1, num_b
       coord_b(1:3, i) = molecule_b%atom_coord(1:3, selatoms_b%idx(i))
       minbou(1:3) = min(minbou(1:3), coord_b(1:3, i))
       maxbou(1:3) = max(maxbou(1:3), coord_b(1:3, i))
     end do
    
    end if

    !ky uncomment to use coordinates in reffile 
    !ky 
    !ky if (allocated(molecule_a%atom_refcoord)) then
    !ky   ! select from reference coordinates
    !ky   do i = 1, num_a
    !ky     coord_a(1:3, i) = molecule_a%atom_refcoord(1:3, selatoms_a%idx(i))
    !ky     minbou(1:3) = min(minbou(1:3), coord_a(1:3, i))
    !ky     maxbou(1:3) = max(maxbou(1:3), coord_a(1:3, i))
    !ky   end do
    !ky
    !ky   do i = 1, num_b
    !ky     coord_b(1:3, i) = molecule_b%atom_refcoord(1:3, selatoms_b%idx(i))
    !ky     minbou(1:3) = min(minbou(1:3), coord_b(1:3, i))
    !ky     maxbou(1:3) = max(maxbou(1:3), coord_b(1:3, i))
    !ky   end do
    !ky
    !ky else
    !ky
    !ky  ! select from input coordinates
    !ky  do i = 1, num_a
    !ky    coord_a(1:3, i) = molecule_a%atom_coord(1:3, selatoms_a%idx(i))
    !ky    minbou(1:3) = min(minbou(1:3), coord_a(1:3, i))
    !ky    maxbou(1:3) = max(maxbou(1:3), coord_a(1:3, i))
    !ky  end do
    !ky
    !ky  do i = 1, num_b
    !ky    coord_b(1:3, i) = molecule_b%atom_coord(1:3, selatoms_b%idx(i))
    !ky    minbou(1:3) = min(minbou(1:3), coord_b(1:3, i))
    !ky    maxbou(1:3) = max(maxbou(1:3), coord_b(1:3, i))
    !ky  end do
    !ky
    !ky end if

    minbou(1:3) = minbou(1:3) - EPS
    maxbou(1:3) = maxbou(1:3) + EPS
    

    ! compute cell parameter
    !
    mincl(1:3)  = distance
    region(1:3) = maxbou(1:3) - minbou(1:3)

    call compute_cell_by_region(region,     &
                                mincl,      &
                                cell_count, &
                                cell_length)


    ! build cell linked-list
    !
    call build_cell_linked_list(minbou(1:3), &
                                region,      &
                                cell_count,  &
                                coord_a,     &
                                list_a)

    call build_cell_linked_list(minbou(1:3), &
                                region,      &
                                cell_count,  &
                                coord_b,     &
                                list_b)

    return

  end subroutine init_select

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    select_atom_contact
  !> @brief        select atoms in less than distance of atoms-A from atoms-B
  !! @authors      NT
  !! @param[in]    distance   : contact distance 
  !! @param[in]    max_atom   : maximum selection atom count
  !! @param[in]    sort_list  : sort list
  !! @param[in]    list_a     : cell-linked list A
  !! @param[in]    list_b     : cell-linked list B
  !! @param[in]    selatoms_b : atom seleciton B
  !! @param[inout] selatoms   : contact atom seleciton
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine select_atom_contact(distance,   &
                                 max_atom,   &
                                 sort_list,  &
                                 list_a,     &
                                 list_b,     &
                                 selatoms_b, &
                                 selatoms)

    ! formal arguments
    real(wp),                 intent(in)    :: distance
    integer,                  intent(in)    :: max_atom
    logical,                  intent(in)    :: sort_list
    type(s_cell_linked_list), intent(in)    :: list_a
    type(s_cell_linked_list), intent(in)    :: list_b
    type(s_selatoms),         intent(in)    :: selatoms_b
    type(s_selatoms),         intent(inout) :: selatoms

    ! local variables
    real(wp)                  :: sqr_d, dx, dy, dz, dxyz
    integer                   :: ccx, ccy, ccz, ccyz
    integer                   :: cxa, cya, cza, cxb, cyb, czb, ca, cb
    integer                   :: i, natom, natom_b, nsel, a, b, alloc_stat

    type(s_sdata_atom),       allocatable :: sort_array(:)
    real(wp),                 allocatable :: atom_dis(:)
    logical,                  allocatable :: atom_flg(:)

    real(wp),                 pointer     :: crda(:,:), crdb(:,:)


    crda => list_a%coord
    crdb => list_b%coord
    sqr_d = distance * distance

    natom_b = size(selatoms_b%idx)


    ! allocate
    !
    allocate(atom_flg(natom_b), &
             atom_dis(natom_b), stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    atom_flg(1:natom_b) = .false.


    ! search
    !
    ccx  = list_a%cell_count(1)
    ccy  = list_a%cell_count(2)
    ccz  = list_a%cell_count(3)
    ccyz = ccy * ccz

    do cxa = 0, ccx - 1
    do cya = 0, ccy - 1
    do cza = 0, ccz - 1
      
      ca = cxa * ccyz + &
           cya * ccz +  &
           cza + 1

      do cxb = cxa - 1, cxa + 1
        if (cxb < 0 .or. cxb >= ccx) cycle
      do cyb = cya - 1, cya + 1
        if (cyb < 0 .or. cyb >= ccy) cycle
      do czb = cza - 1, cza + 1
        if (czb < 0 .or. czb >= ccz) cycle

        cb = cxb * ccyz + &
             cyb * ccz +  &
             czb + 1

        a = list_a%head(ca)

        do while(a /= Empty)

          b = list_b%head(cb)

          do while(b /= Empty)

            dx = crda(1, a) - crdb(1, b)
            dy = crda(2, a) - crdb(2, b)
            dz = crda(3, a) - crdb(3, b)
            dxyz = dx*dx + dy*dy + dz*dz

            if (sqr_d > dxyz) then
              atom_flg(b) = .true.
              atom_dis(b) = dxyz
            end if

            b = list_b%lscl(b)

          end do
          a = list_a%lscl(a)

        end do

      end do
      end do
      end do

    end do
    end do
    end do


    ! setup result
    !

    nsel = 0
    do i = 1, natom_b
      if (atom_flg(i)) &
        nsel = nsel + 1
    end do

    if (max_atom /= 0) then
      natom = min(max_atom, nsel)
    else
      natom = nsel
    end if

    call alloc_selatoms(selatoms, natom)

    if (nsel /= 0) then

      ! sort
      if (sort_list) then

        allocate(sort_array(nsel), stat = alloc_stat)
        if (alloc_stat /= 0) &
          call error_msg_alloc

        nsel = 0
        do i = 1, natom_b
          if (atom_flg(i)) then
            nsel = nsel + 1
            sort_array(nsel)%sq_dis = atom_dis(i)
            sort_array(nsel)%idx    = selatoms_b%idx(i)
          end if
        end do

        call quicksort_atom(sort_array, 1, nsel)

        do i = 1, natom
          selatoms%idx(i) = sort_array(i)%idx
        end do

        deallocate(sort_array, stat = alloc_stat)

        ! not sort
      else

        nsel = 0
        do i = 1, natom_b
          if (atom_flg(i)) then
            nsel = nsel + 1
            if (nsel > natom) &
              exit
            selatoms%idx(nsel) = selatoms_b%idx(i)
          end if
        end do

      end if

    end if


    ! dealloc
    !
    deallocate(atom_flg, atom_dis, stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_dealloc

    return

  end subroutine select_atom_contact

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    select_mole_contact
  !> @brief        select molecule in less than distance of atoms A from atoms B
  !! @authors      NT
  !! @param[in]    distance       : contact distance 
  !! @param[in]    max_atom       : maximum selection atom count
  !! @param[in]    sort_list      : sort list
  !! @param[in]    list_a         : cell-linked list A
  !! @param[in]    list_b         : cell-linked list B
  !! @param[in]    selatoms_b     : atom seleciton B
  !! @param[in]    mole_no_b      : molecule number list
  !! @param[in]    mole_atom_no_b : atom number list
  !! @param[inout] selatoms       : contact atom seleciton
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine select_mole_contact(distance,       &
                                 max_atom,       &
                                 sort_list,      &
                                 list_a,         &
                                 list_b,         &
                                 selatoms_b,     &
                                 mole_no_b,      &
                                 mole_atom_no_b, &
                                 selatoms)

    ! formal arguments
    real(wp),                 intent(in)    :: distance
    integer,                  intent(in)    :: max_atom
    logical,                  intent(in)    :: sort_list
    type(s_cell_linked_list), intent(in)    :: list_a
    type(s_cell_linked_list), intent(in)    :: list_b
    type(s_selatoms),         intent(in)    :: selatoms_b
    integer,                  intent(in)    :: mole_no_b(:)
    integer,                  intent(in)    :: mole_atom_no_b(:)
    type(s_selatoms),         intent(inout) :: selatoms

    ! local variables
    real(wp)                  :: sqr_d, dx, dy, dz, dxyz
    integer                   :: ccx, ccy, ccz, ccyz
    integer                   :: cxa, cya, cza, cxb, cyb, czb, ca, cb
    integer                   :: natom, natom0, nmole, nmole_b, midx, a, b
    integer                   :: i, j, alloc_stat
    integer                   :: first, last

    type(s_sdata_mole),       allocatable :: sort_array(:)
    real(wp),                 allocatable :: mole_dis(:)
    logical,                  allocatable :: mole_flg(:)
    real(wp),                 pointer     :: crda(:,:), crdb(:,:)


    crda => list_a%coord
    crdb => list_b%coord
    sqr_d = distance * distance

    nmole_b = size(mole_atom_no_b)


    ! check molecule atom count
    !
    if (nmole_b == 0) then
      call alloc_selatoms(selatoms, 0)
      return 
    end if

    ! allocate
    !
    allocate(mole_flg(nmole_b), &
             mole_dis(nmole_b), stat = alloc_stat)
    if (alloc_stat /= 0) &
       call error_msg_alloc


    mole_flg(1:nmole_b) = .false.
    mole_dis(1:nmole_b) = sqr_d

    ccx  = list_a%cell_count(1)
    ccy  = list_a%cell_count(2)
    ccz  = list_a%cell_count(3)
    ccyz = ccy * ccz

    do cxa = 0, ccx - 1
    do cya = 0, ccy - 1
    do cza = 0, ccz - 1

      ca = cxa * ccyz + &
           cya * ccz +  &
           cza + 1

      do cxb = cxa - 1, cxa + 1
        if (cxb < 0 .or. cxb >= ccx) cycle
      do cyb = cya - 1, cya + 1
        if (cyb < 0 .or. cyb >= ccy) cycle
      do czb = cza - 1, cza + 1
        if (czb < 0 .or. czb >= ccz) cycle

        cb = cxb * ccyz + &
             cyb * ccz +  &
             czb + 1

        a = list_a%head(ca)

        do while(a /= Empty)

          b = list_b%head(cb)

          do while(b /= Empty)

            dx = crda(1, a) - crdb(1, b)
            dy = crda(2, a) - crdb(2, b)
            dz = crda(3, a) - crdb(3, b)
            dxyz = dx*dx + dy*dy + dz*dz

            if (sqr_d > dxyz) then

              midx = mole_no_b(selatoms_b%idx(b))
              mole_flg(midx) = .true.
              mole_dis(midx) = min(mole_dis(midx), dxyz)

            end if

            b = list_b%lscl(b)

          end do
          a = list_a%lscl(a)

        end do

      end do
      end do
      end do

    end do
    end do
    end do


    ! setup result
    !

    ! sort
    if (sort_list) then

      nmole = 0
      do i = 1, nmole_b
        if (mole_flg(i)) &
          nmole = nmole + 1
      end do

      allocate(sort_array(nmole), stat = alloc_stat)
      if (alloc_stat /= 0) &
        call error_msg_alloc

      nmole = 0
      do i = 1, nmole_b
        if (mole_flg(i)) then

          nmole = nmole + 1

          first = mole_atom_no_b(i)
          if (i /= nmole_b) then
            last = mole_atom_no_b(i+1) - 1
          else
            last = size(mole_no_b)
          end if

          sort_array(nmole)%sq_dis = mole_dis(i)
          sort_array(nmole)%idx1 = first
          sort_array(nmole)%idx2 = last

        end if

      end do

      if (nmole /= 0) &
        call quicksort_mole(sort_array, 1, nmole)
        
      natom = 0
      do i = 1, nmole

        natom0 = sort_array(i)%idx2 - sort_array(i)%idx1 + 1

        if (max_atom /= 0) then
          if (natom == max_atom) &
            exit
          if (natom + natom0 > max_atom) &
            cycle
        end if

        natom = natom + natom0

      end do

      call alloc_selatoms(selatoms, natom)

      natom = 0

      do i = 1, nmole

        natom0 = sort_array(i)%idx2 - sort_array(i)%idx1 + 1

        if (max_atom /= 0) then
          if (natom == max_atom) &
            exit
          if (natom + natom0 > max_atom) &
            cycle
        end if

        do j = sort_array(i)%idx1, sort_array(i)%idx2
          natom = natom + 1
          selatoms%idx(natom) = j
        end do

      end do

      deallocate(sort_array, stat = alloc_stat)

    ! not sort
    else

      natom = 0

      do i = 1, nmole_b
        if (mole_flg(i)) then

          first = mole_atom_no_b(i)

          if (i /= nmole_b) then
            last = mole_atom_no_b(i+1) - 1
          else
            last = size(mole_no_b)
          end if

          natom0 = last - first + 1

          if (max_atom /= 0 .and. natom + natom0 > max_atom) &
            exit

          natom = natom + natom0

        end if
      end do

      call alloc_selatoms(selatoms, natom)

      natom = 0

      do i = 1, nmole_b
        if (mole_flg(i)) then

          first = mole_atom_no_b(i)

          if (i /= nmole_b) then
            last = mole_atom_no_b(i+1) - 1
          else
            last = size(mole_no_b)
          end if

          natom0 = last - first + 1

          if (max_atom /= 0 .and. natom + natom0 > max_atom) &
            exit
          do j = first, last
            natom = natom + 1
            selatoms%idx(natom) = j
          end do

        end if
      end do

    end if


    ! dealloc
    !
    deallocate(mole_dis, mole_flg, stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_dealloc

    return

  end subroutine select_mole_contact

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_resi_mole_no
  !> @brief        make_resi_mole_no
  !! @authors      NT
  !! @param[in]    molecule          : information of molecule
  !! @param[out]   resi_mole_no      : array of residue molecule number
  !! @param[out]   resi_mole_atom_no : array of residue molecule atom number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine make_resi_mole_no(molecule,     &
                               resi_mole_no, &
                               resi_mole_atom_no)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    integer,     allocatable, intent(inout) :: resi_mole_no(:)
    integer,     allocatable, intent(inout) :: resi_mole_atom_no(:)

    ! local variables
    integer                   :: i, natom, nresi
    integer                   :: alloc_stat

    integer,          pointer :: presno(:)
    character(4),     pointer :: psegname(:)


    psegname => molecule%segment_name
    presno   => molecule%residue_no

    natom = molecule%num_atoms
    if (natom == 0) &
      return

    nresi = 1
    do i = 2, natom
      if (psegname(i-1) /= psegname(i) .or. &
          presno  (i-1) /= presno  (i)) &
        nresi = nresi + 1
    end do

    allocate(resi_mole_no(natom), &
             resi_mole_atom_no(nresi), &
             stat = alloc_stat)

    nresi = 1
    resi_mole_atom_no(nresi) = 1
    resi_mole_no(1) = nresi

    do i = 2, natom
      if (psegname(i-1) /= psegname(i) .or. &
          presno  (i-1) /= presno  (i)) then

        nresi = nresi + 1
        resi_mole_atom_no(nresi) = i
      end if
      resi_mole_no(i) = nresi
    end do

    return

  end subroutine make_resi_mole_no

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    quicksort_atom
  !> @brief        sorting function
  !! @authors      NT
  !! @param[inout] a     : array
  !! @param[in]    start : start index
  !! @param[in]    end   : end index
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine quicksort_atom(a, start, end)

    ! formal arguments
    type(s_sdata_atom),      intent(inout) :: a(:)
    integer,                 intent(in)    :: start
    integer,                 intent(in)    :: end

    ! local variables
    type(s_sdata_atom)       :: t, x
    integer                  :: i, j


    x = a((start + end) / 2)
    i = start
    j = end

    do

      do while (a(i)%sq_dis < x%sq_dis)
        i = i + 1
      end do

      do while (x%sq_dis < a(j)%sq_dis)
        j = j - 1
      end do

      if (i >= j) &
        exit

      t = a(i)
      a(i) = a(j)
      a(j) = t
      i = i + 1
      j = j - 1

    end do

    if (start < i - 1) &
      call quicksort_atom(a, start, i - 1)
    if (j + 1 < end) &
      call quicksort_atom(a, j + 1, end)

    return

  end subroutine quicksort_atom

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    quicksort_mole
  !> @brief        sorting function
  !! @authors      NT
  !! @param[inout] a     : array
  !! @param[in]    start : start index
  !! @param[in]    end   : end index
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine quicksort_mole(a, start, end)

    ! formal arguments
    type(s_sdata_mole),      intent(inout) :: a(:)
    integer,                 intent(in)    :: start
    integer,                 intent(in)    :: end

    ! local variables
    type(s_sdata_mole)       :: t, x
    integer                  :: i, j


    x = a((start + end) / 2)
    i = start
    j = end

    do

      do while (a(i)%sq_dis < x%sq_dis)
        i = i + 1
      end do

      do while (x%sq_dis < a(j)%sq_dis)
        j = j - 1
      end do

      if (i >= j) &
        exit

      t = a(i)
      a(i) = a(j)
      a(j) = t
      i = i + 1
      j = j - 1

    end do

    if (start < i - 1) &
      call quicksort_mole(a, start, i - 1)
    if (j + 1 < end) &
      call quicksort_mole(a, j + 1, end)

    return

  end subroutine quicksort_mole

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    build_cell_linked_list
  !> @brief        build cell-linked list
  !! @authors      NT
  !! @param[in]    origin     : simulation box origin
  !! @param[in]    region     : simulation box region
  !! @param[in]    cell_count : cell count for each axis
  !! @param[in]    coord      : array of 3d-coordinate
  !! @param[out]   list       : cell linked-list
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine build_cell_linked_list(origin,     &
                                    region,     &
                                    cell_count, &
                                    coord,      &
                                    list)

    ! formal arguments
    real(wp),                 intent(in)    :: origin(3)
    real(wp),                 intent(in)    :: region(3)
    integer,                  intent(in)    :: cell_count(3)
    real(wp),         target, intent(in)    :: coord(:,:)
    type(s_cell_linked_list), intent(out)   :: list

    ! local variables
    integer                   :: total_cell, total_cell_yz, total_coord
    integer                   :: i, alloc_stat, dealloc_stat
    integer                   :: ci, cix, ciy, ciz


    ! check parameter
    !
    if (region(1) <= 0.0_wp .or. &
        region(2) <= 0.0_wp .or. &
        region(3) <= 0.0_wp .or. &
        cell_count(1) <= 0 .or. &
        cell_count(2) <= 0 .or. &
        cell_count(3) <= 0) &
      call error_msg('build_cell_linked_list> ERROR : cell linked-list'//&
                     ' parameter.')
       
    total_cell_yz = cell_count(2) * cell_count(3)
    total_cell    = cell_count(1) * total_cell_yz
    total_coord   = size(coord(1,:))

    if (total_coord == 0) &
      call error_msg('build_cell_linked_list> ERROR : There are no'//&
                     ' coordinates.')

    list%origin(1:3)      = origin(1:3)
    list%region(1:3)      = region(1:3)
    list%cell_count(1:3)  = cell_count(1:3)
    list%cell_length(1:3) = region(1:3) / real(cell_count(1:3),wp)

    list%coord => coord

    !write(MsgOut, *) ' '
    !write(MsgOut, fmt='(" build_cell_linked_list> origin      : '//&
    !                  '(",3f8.3,")")'), list%origin(1:3)
    !write(MsgOut, fmt='("                         region      : '//&
    !                  '(",3f8.3,")")'), list%region(1:3)
    !write(MsgOut, fmt='("                         cell count  : '//&
    !                  '(",3i3,")")'), list%cell_count(1:3)
    !write(MsgOut, fmt='("                         cell length : '//&
    !                  '(",3f8.3,")")'), list%cell_length(1:3)
    !write(MsgOut, fmt='("                         total cell  : ",i8)'), &
    !     total_cell
    !write(MsgOut, fmt='("                         total coord : ",i8)'), &
    !     total_coord


    ! allocate list
    !
    if (allocated(list%head)) then
      deallocate (list%head, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(list%lscl)) then
      deallocate (list%lscl, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    allocate(list%head(total_cell),  &
             list%lscl(total_coord), &
             stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc


    ! build
    !

    do i = 1, total_cell
       list%head(i) = EMPTY
    end do

    do i = 1, total_coord
      cix = floor((coord(1, i) - list%origin(1)) / list%cell_length(1))
      ciy = floor((coord(2, i) - list%origin(2)) / list%cell_length(2))
      ciz = floor((coord(3, i) - list%origin(3)) / list%cell_length(3))

      ci = cix * total_cell_yz +      &
           ciy * list%cell_count(3) + &
           ciz + 1

      if (ci < 1 .or. ci > total_cell) &
        call error_msg('build_cell_linked_list> coordinate is out of region.')

      list%lscl(i)  = list%head(ci)
      list%head(ci) = i

    end do

    return

  end subroutine build_cell_linked_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_cell_linked_list
  !> @brief        dealloc cell-linked list
  !! @authors      NT
  !! @param[inout] list : cell linked-list
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_cell_linked_list(list)

    ! formal arguments
    type(s_cell_linked_list), intent(inout) :: list

    ! local variables
    integer                   :: dealloc_stat


    if (allocated(list%head)) then
      deallocate (list%head, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(list%lscl)) then
      deallocate (list%lscl, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    nullify(list%coord)

    return

  end subroutine dealloc_cell_linked_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    search_neighbour_coord
  !> @brief        search neighbour coordinates of target coordinate
  !! @authors      NT
  !! @param[in]    list       : cell linked-list
  !! @param[in]    distance   : distance to judge as a neighbour coord
  !! @param[in]    coord_tgt  : target coordinate
  !! @param[out]   coord_list : coordinate index list
  !! @param[out]   num_list   : number of index list
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine search_neighbour_coord(list,       &
                                    distance,   &
                                    coord_tgt,  &
                                    coord_list, &
                                    num_list)

    ! formal arguments
    type(s_cell_linked_list), intent(in)    :: list
    real(wp),                 intent(in)    :: distance
    real(wp),                 intent(in)    :: coord_tgt(3)
    integer,                  intent(inout) :: coord_list(:)
    integer,                  intent(inout) :: num_list

    ! local variables
    real(wp)                  :: dis2, d2
    integer                   :: i, cx, cy, cz, dcx, dcy, dcz, dc
    integer                   :: ccx, ccy, ccz, ccyz

    real(wp),         pointer :: coord(:,:)


    num_list = 0

    ccx  = list%cell_count(1)
    ccy  = list%cell_count(2)
    ccz  = list%cell_count(3)
    ccyz = ccy * ccz

    cx = floor((coord_tgt(1) - list%origin(1)) / list%cell_length(1))
    cy = floor((coord_tgt(2) - list%origin(2)) / list%cell_length(2))
    cz = floor((coord_tgt(3) - list%origin(3)) / list%cell_length(3))

    dis2 = distance**2

    coord => list%coord

    if (cx < 0 .or. cx >= ccx .or. &
        cy < 0 .or. cy >= ccy .or. &
        cz < 0 .or. cz >= ccz)     &
      call error_msg('search_cell_linked_list> coordinate is out of region.')

    do dcx = cx - 1, cx + 1
      if (dcx < 0 .or. dcx >= ccx) &
        cycle
    do dcy = cy - 1, cy + 1
      if (dcy < 0 .or. dcy >= ccy) &
        cycle
    do dcz = cz - 1, cz + 1
      if (dcz < 0 .or. dcz >= ccz) &
        cycle

       dc = dcx * ccyz + &
            dcy * ccz +  &
            dcz + 1
             
       i = list%head(dc)

       do while(i /= EMPTY)
          d2 = (coord_tgt(1) - coord(1,i))**2 + &
               (coord_tgt(2) - coord(2,i))**2 + &
               (coord_tgt(3) - coord(3,i))**2
          if (d2 > EPS .and. d2 <= dis2) then
            num_list = num_list + 1
            coord_list(num_list) = i
          end if
          i = list%lscl(i)
       end do

    end do
    end do
    end do

    return

  end subroutine search_neighbour_coord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_cell_by_coords
  !> @brief        compute the most suitable cell parameter 
  !! @authors      NT
  !! @param[in]    coord        : array of 3d-coordinate
  !! @param[in]    min_cell_len : minimum cell length
  !! @param[out]   origin       : simulation box origin
  !! @param[out]   region       : simulation box region
  !! @param[out]   cell_count   : computed cell count
  !! @param[out]   cell_length  : computed cell length
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_cell_by_coords(coord,        &
                                    min_cell_len, &
                                    origin,       &
                                    region,       &
                                    cell_count,   &
                                    cell_length)

    ! formal arguments
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(in)    :: min_cell_len(3)
    real(wp),                intent(out)   :: origin(3)
    real(wp),                intent(out)   :: region(3)
    integer,                 intent(out)   :: cell_count(3) 
    real(wp),                intent(out)   :: cell_length(3)

    ! local variables
    real(wp)                 :: minb(3), maxb(3)
    integer                  :: i, num


    ! compute boundary
    num = size(coord(1,:))
    if (num <= 0) &
      call error_msg('compute_cell_by_coords')

    minb(1:3) = coord(1:3, 1)
    maxb(1:3) = coord(1:3, 1)
    
    do i = 2, num
      minb(1:3) = min(minb(1:3), coord(1:3, i))
      maxb(1:3) = max(maxb(1:3), coord(1:3, i))
    end do
    
    minb(1:3) = minb(1:3) - EPS
    maxb(1:3) = maxb(1:3) + EPS

    origin(1:3) = minb(1:3)
    region(1:3) = maxb(1:3) - minb(1:3)

    !write(MsgOut, fmt='(" compute_cell_by_coords> origin : '//&
    !                  '(",3f8.3,")")'), origin(1:3)
    !write(MsgOut, fmt='("                         region : '//&
    !                  '(",3f8.3,")")'), region(1:3)

    call compute_cell_by_region(region,       &
                                min_cell_len, &
                                cell_count,   &
                                cell_length)

    return

  end subroutine compute_cell_by_coords

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_cell_by_region
  !> @brief        compute the most suitable cell parameter 
  !! @authors      NT
  !! @param[in]    region       : simulation box region
  !! @param[in]    min_cell_len : minimum cell length
  !! @param[out]   cell_count   : computed cell count
  !! @param[out]   cell_length  : computed cell length
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_cell_by_region(region,       &
                                    min_cell_len, &
                                    cell_count,   &
                                    cell_length)

    ! formal arguments
    real(wp),                intent(in)    :: region(3)
    real(wp),                intent(in)    :: min_cell_len(3)
    integer,                 intent(out)   :: cell_count(3) 
    real(wp),                intent(out)   :: cell_length(3)

    ! local variables
    real(wp)                 :: r_cell_count(3)


    if (min_cell_len(1) == 0.0_wp .or. &
        min_cell_len(2) == 0.0_wp .or. &
        min_cell_len(3) == 0.0_wp) &
      call error_msg('compute_cell_by_region> ERROR : min_cell_len is zero.')

    r_cell_count(1:3) = region(1:3) / min_cell_len(1:3)

    cell_count(1:3) = int(floor(r_cell_count(1:3)))
    if (cell_count(1) < 1) cell_count(1) = 1
    if (cell_count(2) < 1) cell_count(2) = 1
    if (cell_count(3) < 1) cell_count(3) = 1

    cell_length(1:3) = region(1:3) / real(cell_count(1:3),wp)

    !write(MsgOut, fmt='(" compute_cell_by_region> cell count  : '//&
    !                  '(",3i3,")")'), cell_count(1:3)
    !write(MsgOut, fmt='("                         cell length : '//&
    !                  '(",3f8.3,")")'), cell_length(1:3)

    return

  end subroutine compute_cell_by_region

end module select_contacts_mod
