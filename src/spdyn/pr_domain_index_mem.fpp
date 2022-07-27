!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pr_domain_index_mod
!> @brief   handle domain index file (memory version)
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

module pr_domain_index_mod

  use pr_huge_molecule_mod
  use sp_boundary_str_mod
  use fileio_mod
  use messages_mod
  use constants_mod
  use mpi_parallel_mod

  implicit none
  private

  ! public
  !

  ! structures
  type, public :: s_selatoms_list
    integer                    :: n
    integer, allocatable       :: i(:)
  end type s_selatoms_list

  type, public :: s_domain_index
    integer                    :: natom_list = 0
    integer                    :: nbond_list = 0
    integer                    :: nangl_list = 0
    integer                    :: ndihe_list = 0
    integer                    :: nimpr_list = 0
    integer                    :: ncmap_list = 0
    integer                    :: nwater_list = 0
    integer                    :: nsolute_list = 0
    integer, allocatable       :: iatom_list(:)
    integer, allocatable       :: ibond_list(:)
    integer, allocatable       :: iangl_list(:)
    integer, allocatable       :: idihe_list(:)
    integer, allocatable       :: iimpr_list(:)
    integer, allocatable       :: icmap_list(:)
    integer, allocatable       :: iwater_list(:)
    integer, allocatable       :: isolute_list(:)

    type(s_selatoms_list), allocatable :: selatoms_list(:)

  end type s_domain_index

  ! variables
  character(200), public       :: di_file_dir = './'
  character(200), public       :: di_filename = ''

  ! subroutine
  public  :: di_create
  public  :: di_delete
  public  :: di_get_index


  ! private
  !

  ! constants
  integer,           parameter :: GrowSize = 4096

  ! variables
  type(s_domain_index), target, allocatable :: domain_index_list(:)

  real(wp)                     :: di_bsize_x
  real(wp)                     :: di_bsize_y
  real(wp)                     :: di_bsize_z
  real(wp)                     :: di_borg_x
  real(wp)                     :: di_borg_y
  real(wp)                     :: di_borg_z
  real(wp)                     :: di_dsize_x
  real(wp)                     :: di_dsize_y
  real(wp)                     :: di_dsize_z
  integer                      :: di_dnum_x
  integer                      :: di_dnum_y
  integer                      :: di_dnum_z

  ! subroutines
  private :: di_get_domain_contains
  private :: di_write_index
  private :: di_realloc_memory
  private :: di_quicksort
  private :: di_show_progress

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    di_create
  !> @brief        create domain index file
  !! @authors      NT
  !! @param[in]    boundary : boundary condition information
  !! @param[in]    domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine di_create(boundary)

    ! formal arguments
    type(s_boundary),        intent(in)    :: boundary

    ! local variables
    integer                  :: i, j, iatm, idom, ndom
    type(s_domain_index), pointer :: di


    ! setup global variables
    !
    di_bsize_x = boundary%box_size_x
    di_bsize_y = boundary%box_size_y
    di_bsize_z = boundary%box_size_z

    di_borg_x  = boundary%origin_x
    di_borg_y  = boundary%origin_y
    di_borg_z  = boundary%origin_z

    di_dnum_x  = boundary%num_domain(1)
    di_dnum_y  = boundary%num_domain(2)
    di_dnum_z  = boundary%num_domain(3)

    di_dsize_x = di_bsize_x / real(di_dnum_x, wp)
    di_dsize_y = di_bsize_y / real(di_dnum_y, wp)
    di_dsize_z = di_bsize_z / real(di_dnum_z, wp)

    di_filename = hm_work_name

    if (main_rank) then

      write(MsgOut,'(a)')      'Di_Create> '
      write(MsgOut,'(a,f8.2)') '     boundary size x : ', di_bsize_x
      write(MsgOut,'(a,f8.2)') '                 , y : ', di_bsize_y
      write(MsgOut,'(a,f8.2)') '                 , z : ', di_bsize_z
      write(MsgOut,'(a,f8.2)') '   boundary origin x : ', di_borg_x
      write(MsgOut,'(a,f8.2)') '                 , y : ', di_borg_y
      write(MsgOut,'(a,f8.2)') '                 , z : ', di_borg_z
      write(MsgOut,'(a)')      ' '
      write(MsgOut,'(a,i8)')   '  number of domain x : ', di_dnum_x
      write(MsgOut,'(a,i8)')   '                 , y : ', di_dnum_y
      write(MsgOut,'(a,i8)')   '                 , z : ', di_dnum_z
      write(MsgOut,'(a,f8.2)') '       domain size x : ', di_dsize_x
      write(MsgOut,'(a,f8.2)') '                 , y : ', di_dsize_y
      write(MsgOut,'(a,f8.2)') '                 , z : ', di_dsize_z
      write(MsgOut,'(a)')      ' '

    end if

    ! allocate domain index list
    !

    ndom = di_dnum_x*di_dnum_y*di_dnum_z

    allocate(domain_index_list(ndom))

    do idom = 1, ndom
      di => domain_index_list(idom)
      allocate(di%iatom_list(GrowSize))
      allocate(di%ibond_list(GrowSize))
      allocate(di%iangl_list(GrowSize))
      allocate(di%idihe_list(GrowSize))
      allocate(di%iimpr_list(GrowSize))
      allocate(di%icmap_list(GrowSize))
      allocate(di%iwater_list(GrowSize))
      allocate(di%isolute_list(GrowSize))
      allocate(di%selatoms_list(hm_num_selatoms_list))
      do j = 1, hm_num_selatoms_list
        allocate(di%selatoms_list(j)%i(GrowSize))
        di%selatoms_list(j)%n=0
      end do
    end do

    ! atom_list
    !
    if (main_rank) &
      write(MsgOut,'(a)') 'Di_Create> Atom-list: '

    do i = 1, hm_num_atoms
      idom = di_get_domain_contains(i)
      di => domain_index_list(idom)
      call di_write_index(di%iatom_list, di%natom_list, i)
      call di_show_progress(i, hm_num_atoms, 'atom')
    end do

    if (main_rank) &
      write(MsgOut,'(a)') '  done.'

    ! bond_list
    !
    if (main_rank) &
      write(MsgOut,'(a)') 'Di_Create> Bond-list: '

    do i = 1, hm_num_bonds
      iatm = hm_bond_list(1, i)
      idom = di_get_domain_contains(iatm)
      di => domain_index_list(idom)
      call di_write_index(di%ibond_list, di%nbond_list, i)
      call di_show_progress(i, hm_num_bonds, 'bond')
    end do

    if (main_rank) &
      write(MsgOut,'(a)') '  done.'

    ! angl_list
    !
    if (main_rank) &
      write(MsgOut,'(a)') 'Di_Create> Angl-list: '

    do i = 1, hm_num_angles
      iatm = hm_angl_list(1, i)
      idom = di_get_domain_contains(iatm)
      di => domain_index_list(idom)
      call di_write_index(di%iangl_list, di%nangl_list, i)
      call di_show_progress(i, hm_num_angles, 'angl')
    end do

    if (main_rank) &
      write(MsgOut,'(a)') '  done.'


    ! dihe_list
    !
    if (main_rank) &
      write(MsgOut,'(a)') 'Di_Create> Dihe-list: '

    do i = 1, hm_num_dihedrals
      iatm = hm_dihe_list(1, i)
      idom = di_get_domain_contains(iatm)
      di => domain_index_list(idom)
      call di_write_index(di%idihe_list, di%ndihe_list, i)
      call di_show_progress(i, hm_num_dihedrals, 'dihe')
    end do

    if (main_rank) &
      write(MsgOut,'(a)') '  done.'


    ! impr_list
    !
    if (main_rank) &
      write(MsgOut,'(a)') 'Di_Create> Impr-list: '

    do i = 1, hm_num_impropers
      iatm = hm_impr_list(1, i)
      idom = di_get_domain_contains(iatm)
      di => domain_index_list(idom)
      call di_write_index(di%iimpr_list, di%nimpr_list, i)
      call di_show_progress(i, hm_num_impropers, 'impr')
    end do

    if (main_rank) &
      write(MsgOut,'(a)') '  done.'


    ! cmap_list
    !
    if (main_rank) &
      write(MsgOut,'(a)') 'Di_Create> Cmap-list: '

    do i = 1, hm_num_cmaps
      iatm = hm_cmap_list(1, i)
      idom = di_get_domain_contains(iatm)
      di => domain_index_list(idom)
      call di_write_index(di%icmap_list, di%ncmap_list, i)
      call di_show_progress(i, hm_num_cmaps, 'cmap')
    end do

    if (main_rank) &
      write(MsgOut,'(a)') '  done.'


    ! water_list
    !
    if (main_rank) &
      write(MsgOut,'(a)') 'Di_Create> Water-list: '

    do i = 1, hm_num_water_list
      iatm = hm_water_list(1, i)
      idom = di_get_domain_contains(iatm)
      di => domain_index_list(idom)
      call di_write_index(di%iwater_list, di%nwater_list, i)
      call di_show_progress(i, hm_num_water_list, 'water_list')
    end do

    if (main_rank) &
      write(MsgOut,'(a)') '  done.'


    ! solute_list
    !
    if (main_rank) &
      write(MsgOut,'(a)') 'Di_Create> Solute-list: '

    do i = 1, hm_num_solute_list
      iatm = hm_solute_list(i)
      idom = di_get_domain_contains(iatm)
      di => domain_index_list(idom)
      call di_write_index(di%isolute_list, di%nsolute_list, i)
      call di_show_progress(i, hm_num_solute_list, 'solute_list')
    end do

    if (main_rank) &
      write(MsgOut,'(a)') '  done.'


    !  selatoms_list
    if (main_rank) &
      write(MsgOut,'(a)') 'Di_Create> SelAtoms_list: '

    do j = 1, hm_num_selatoms_list
      do i = 1, hm_num_selatoms(j)
        iatm = hm_selatoms(i, j)
        idom = di_get_domain_contains(iatm)
        di => domain_index_list(idom)
        call di_write_index(di%selatoms_list(j)%i, di%selatoms_list(j)%n, i)
        call di_show_progress(i, hm_num_solute_list, 'selatoms')
      end do
    end do

    if (main_rank) &
      write(MsgOut,'(a)') '  done.'

    return

  end subroutine di_create

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    di_delete
  !> @brief        delete domain index file
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine di_delete

    ! local variables
    integer                  :: idom, ndom, j
    type(s_domain_index), pointer :: di


    if (main_rank) then

      write(MsgOut,'(a)') 'Di_Delete> cleanup memory.'
      write(MsgOut,'(a)') ' '

    end if

    ndom = di_dnum_x*di_dnum_y*di_dnum_z

    do idom = 1, ndom
      di => domain_index_list(idom)
      deallocate(di%iatom_list)
      deallocate(di%ibond_list)
      deallocate(di%iangl_list)
      deallocate(di%idihe_list)
      deallocate(di%iimpr_list)
      deallocate(di%icmap_list)
      deallocate(di%iwater_list)
      deallocate(di%isolute_list)
      do j = 1, hm_num_selatoms_list
        deallocate(di%selatoms_list(j)%i)
      end do
    end do

    return

  end subroutine di_delete

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    di_get_index
  !> @brief        get index information for domain
  !! @authors      NT
  !! @param[in]    idom         : index of domain
  !! @param[inout] domain_index : domain index information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine di_get_index(idom, domain_index)

    ! parameters
    integer, parameter    :: X=1, Y=2, Z=3, Lo=1, Cn=2, Up=3

    ! formal arguments
    integer,                      intent(in)    :: idom
    type(s_domain_index), target, intent(inout) :: domain_index

    ! local variables
    integer               :: id(X:Z)          ! domain x, y, z index
    integer               :: nd(X:Z)          ! number of domain x, y, z
    integer               :: idLU(Lo:Up, X:Z) ! domain x, y, z lo-up index
    integer               :: axis, ix, iy, iz, jdom, nlist, ilist, n, i
    integer               :: list(27)

    type(s_domain_index), pointer :: di
    type(s_domain_index), pointer :: di0

    
    di => domain_index

    ! compute neighbour domain
    !

    nd(X) = di_dnum_x
    nd(Y) = di_dnum_y
    nd(Z) = di_dnum_z

    id(X) = mod( idom - 1,                nd(X)) + 1
    id(Y) = mod((idom - 1)/ nd(X),        nd(Y)) + 1
    id(Z) = mod((idom - 1)/(nd(X)*nd(Y)), nd(Z)) + 1

    idLU(Lo:Up, X:Z) = 0

    do axis = X, Z
    
      if (id(axis) == 1) then

        if (nd(axis) > 2) idLU(Lo, axis) = nd(axis)
        if (nd(axis) > 1) idLU(Up, axis) = id(axis) + 1

      else if (id(axis) == nd(axis)) then

        if (nd(axis) > 1) idLU(Lo, axis) = id(axis) - 1
        if (nd(axis) > 2) idLU(Up, axis) = 1

      else

        idLU(Lo, axis) = id(axis) - 1
        idLU(Up, axis) = id(axis) + 1

      end if

    end do

    idLU(Cn, X:Z) = id(X:Z)


    ! read domain index file
    !
    nlist = 0

    do ix = Lo, Up
    do iy = Lo, Up
    do iz = Lo, Up

      if (idLU(ix,X) == 0 .or. &
          idLU(iy,Y) == 0 .or. &
          idLU(iz,Z) == 0) &
        cycle

      jdom = (idLU(ix,X)-1) + &
             (idLU(iy,Y)-1) * nd(X) + &
             (idLU(iz,Z)-1) * nd(X) * nd(Y) + 1

      !id(X) = mod( jdoms - 1,                nd(X)) + 1
      !id(Y) = mod((jdoms - 1)/ nd(X),        nd(Y)) + 1
      !id(Z) = mod((jdoms - 1)/(nd(X)*nd(Y)), nd(Z)) + 1

      nlist = nlist + 1
      list(nlist) = jdom

    end do
    end do
    end do


    ! merge domain index
    !

    if (hm_num_selatoms_list /= 0 .and. .not. allocated(di%selatoms_list)) &
      allocate(di%selatoms_list(hm_num_selatoms_list))

    di%natom_list   = 0
    di%nbond_list   = 0
    di%nangl_list   = 0
    di%ndihe_list   = 0
    di%nimpr_list   = 0
    di%ncmap_list   = 0
    di%nwater_list  = 0
    di%nsolute_list = 0
    di%selatoms_list(1:hm_num_selatoms_list)%n = 0
    do ilist = 1, nlist
      di0 => domain_index_list(list(ilist))
      di%natom_list   = di%natom_list   + di0%natom_list
      di%nbond_list   = di%nbond_list   + di0%nbond_list
      di%nangl_list   = di%nangl_list   + di0%nangl_list
      di%ndihe_list   = di%ndihe_list   + di0%ndihe_list
      di%nimpr_list   = di%nimpr_list   + di0%nimpr_list
      di%ncmap_list   = di%ncmap_list   + di0%ncmap_list
      di%nwater_list  = di%nwater_list  + di0%nwater_list
      di%nsolute_list = di%nsolute_list + di0%nsolute_list
      do i = 1, hm_num_selatoms_list
        di%selatoms_list(i)%n = di%selatoms_list(i)%n + &
                                          di0%selatoms_list(i)%n
      end do
    end do

    call di_realloc_memory(di%iatom_list,   di%natom_list  )
    call di_realloc_memory(di%ibond_list,   di%nbond_list  )
    call di_realloc_memory(di%iangl_list,   di%nangl_list  )
    call di_realloc_memory(di%idihe_list,   di%ndihe_list  )
    call di_realloc_memory(di%iimpr_list,   di%nimpr_list  )
    call di_realloc_memory(di%icmap_list,   di%ncmap_list  )
    call di_realloc_memory(di%iwater_list,  di%nwater_list )
    call di_realloc_memory(di%isolute_list, di%nsolute_list)
    do i = 1, hm_num_selatoms_list
      call di_realloc_memory(di%selatoms_list(i)%i, di%selatoms_list(i)%n)
    end do

    di%natom_list   = 0
    di%nbond_list   = 0
    di%nangl_list   = 0
    di%ndihe_list   = 0
    di%nimpr_list   = 0
    di%ncmap_list   = 0
    di%nwater_list  = 0
    di%nsolute_list = 0
    di%selatoms_list(1:hm_num_selatoms_list)%n = 0
    do ilist = 1, nlist
      di0 => domain_index_list(list(ilist))

      n = di%natom_list
      di%iatom_list(n+1:n+di0%natom_list) = di0%iatom_list(1:di0%natom_list)
      di%natom_list   = di%natom_list   + di0%natom_list

      n = di%nbond_list
      di%ibond_list(n+1:n+di0%nbond_list) = di0%ibond_list(1:di0%nbond_list)
      di%nbond_list   = di%nbond_list   + di0%nbond_list

      n = di%nangl_list
      di%iangl_list(n+1:n+di0%nangl_list) = di0%iangl_list(1:di0%nangl_list)
      di%nangl_list   = di%nangl_list   + di0%nangl_list

      n = di%ndihe_list
      di%idihe_list(n+1:n+di0%ndihe_list) = di0%idihe_list(1:di0%ndihe_list)
      di%ndihe_list   = di%ndihe_list   + di0%ndihe_list

      n = di%nimpr_list
      di%iimpr_list(n+1:n+di0%nimpr_list) = di0%iimpr_list(1:di0%nimpr_list)
      di%nimpr_list   = di%nimpr_list   + di0%nimpr_list

      n = di%ncmap_list
      di%icmap_list(n+1:n+di0%ncmap_list) = di0%icmap_list(1:di0%ncmap_list)
      di%ncmap_list   = di%ncmap_list   + di0%ncmap_list

      n = di%nwater_list
      di%iwater_list(n+1:n+di0%nwater_list) = &
           di0%iwater_list(1:di0%nwater_list)
      di%nwater_list  = di%nwater_list  + di0%nwater_list

      n = di%nsolute_list
      di%isolute_list(n+1:n+di0%nsolute_list) = &
           di0%isolute_list(1:di0%nsolute_list)
      di%nsolute_list = di%nsolute_list + di0%nsolute_list

      do i = 1, hm_num_selatoms_list
        n = di%selatoms_list(i)%n
        di%selatoms_list(i)%i(n+1:n+di0%selatoms_list(i)%n) = &
           di0%selatoms_list(i)%i(1:di0%selatoms_list(i)%n)
        di%selatoms_list(i)%n = di%selatoms_list(i)%n + di0%selatoms_list(i)%n
      end do

    end do


    ! sort index array
    !
    if (di%natom_list > 0) &
      call di_quicksort(di%iatom_list,   1, di%natom_list)
    if (di%nbond_list > 0) &
      call di_quicksort(di%ibond_list,   1, di%nbond_list)
    if (di%nangl_list > 0) &
      call di_quicksort(di%iangl_list,   1, di%nangl_list)
    if (di%ndihe_list > 0) &
      call di_quicksort(di%idihe_list,   1, di%ndihe_list)
    if (di%nimpr_list > 0) &
      call di_quicksort(di%iimpr_list,   1, di%nimpr_list)
    if (di%ncmap_list > 0) &
      call di_quicksort(di%icmap_list,   1, di%ncmap_list)
    if (di%nwater_list > 0) &
      call di_quicksort(di%iwater_list,  1, di%nwater_list)
    if (di%nsolute_list > 0) &
      call di_quicksort(di%isolute_list, 1, di%nsolute_list)
    do i = 1, hm_num_selatoms_list
      if (di%selatoms_list(i)%n > 0) &
        call di_quicksort(di%selatoms_list(i)%i, 1, di%selatoms_list(i)%n)
    end do

    return

  end subroutine di_get_index

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      di_get_domain_contains
  !> @brief        check that domain contains the indicated atom
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function di_get_domain_contains(atom_index)

    ! return value
    integer                  :: di_get_domain_contains

    ! formal arguments
    integer,                 intent(in)    :: atom_index

    ! local variables
    real(wp)                 :: shift_x, shift_y, shift_z
    real(wp)                 :: move_x, move_y, move_z
    integer                  :: idx, idy, idz


    ! coordinate shifted against the origin
    !
    shift_x = hm_atom_coord(1,atom_index) - di_borg_x
    shift_y = hm_atom_coord(2,atom_index) - di_borg_y
    shift_z = hm_atom_coord(3,atom_index) - di_borg_z

    ! coordinate shifted to the first quadrant and set into the boundary box
    move_x  = di_bsize_x*0.5_wp - di_bsize_x*anint(shift_x/di_bsize_x)
    move_y  = di_bsize_y*0.5_wp - di_bsize_y*anint(shift_y/di_bsize_y)
    move_z  = di_bsize_z*0.5_wp - di_bsize_z*anint(shift_z/di_bsize_z)
    shift_x = shift_x + move_x
    shift_y = shift_y + move_y
    shift_z = shift_z + move_z

    ! assign which domain
    !
    idx = int(shift_x/di_dsize_x)
    idy = int(shift_y/di_dsize_y)
    idz = int(shift_z/di_dsize_z)
    if (idx == di_dnum_x) idx = idx - 1
    if (idy == di_dnum_y) idy = idy - 1
    if (idz == di_dnum_z) idz = idz - 1
    di_get_domain_contains = 1 + idx + &
                                 idy*di_dnum_x + &
                                 idz*di_dnum_x*di_dnum_y

    return

  end function di_get_domain_contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    di_write_index
  !> @brief        write the atom index to domain index list
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine di_write_index(list, nlist, index)

    ! formal arguments
    integer,                 allocatable   :: list(:)
    integer,                 intent(inout) :: nlist
    integer,                 intent(in)    :: index

    integer                  :: i, nlist_new
    integer,     allocatable :: work(:)


    nlist = nlist + 1

    if (nlist > size(list)) then
      nlist_new = size(list)+GrowSize
      allocate(work(nlist_new))
      do i = 1, size(list)
        work(i) = list(i)
      end do
      deallocate(list)
      allocate(list(nlist_new))
      do i = 1, nlist_new
        list(i) = work(i)
      end do
      deallocate(work)
    end if

    list(nlist) = index

    return

  end subroutine di_write_index

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    di_realloc_memory
  !> @brief        reallocate memory
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine di_realloc_memory(var, newsize)

    ! formal arguments
    integer,    allocatable, intent(inout) :: var(:)
    integer,                 intent(in)    :: newsize


    if (.not. allocated(var)) then
      allocate(var(newsize))

    else if (newsize > size(var)) then
      if (allocated(var)) &
        deallocate(var)
      allocate(var(newsize))

    end if

    return

  end subroutine di_realloc_memory

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    di_quicksort
  !> @brief        sorting function
  !! @authors      NT
  !! @param[inout] a   : list of integers
  !! @param[in]    numa: number of array a
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine di_quicksort(a, start, end)

    ! formal arguments
    integer,                 intent(inout) ::  a(*)
    integer,                 intent(in)    :: start, end

    ! local variables
    integer                  :: i, j, t, x


    x = a((start + end) / 2)
    i = start
    j = end
    do
      do while (a(i) < x)
        i = i + 1
      end do
      do while (x < a(j))
        j = j - 1
      end do
      if (i >= j) exit
      t = a(i)
      a(i) = a(j)
      a(j) = t
      i = i + 1
      j = j - 1
    end do

    if (start < i - 1) call di_quicksort(a, start, i - 1)
    if (j + 1 < end)  call di_quicksort(a, j + 1, end)

    return

  end subroutine di_quicksort

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    di_show_progress
  !> @brief        show progrss
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine di_show_progress(ip, nt, str)

    ! parameter
    integer,                 parameter     :: Unit = 1000000

    ! formal arguments
    integer,                 intent(in)    :: ip
    integer,                 intent(in)    :: nt
    character(*),            intent(in)    :: str


    if (mod(ip, Unit) /= 0 .or. .not. main_rank) &
      return

    write(MsgOut,'(a,i9,a,i9,a)') '   '//trim(str)//': ', &
                           ip / Unit, &
                           '/',       &
                           nt / Unit, ' million'

    return

  end subroutine di_show_progress

end module pr_domain_index_mod
