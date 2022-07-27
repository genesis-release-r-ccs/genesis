!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pr_domain_index_mod
!> @brief   handle domain index (file version)
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
    integer                    :: natom_list
    integer                    :: nbond_list
    integer                    :: nangl_list
    integer                    :: ndihe_list
    integer                    :: nimpr_list
    integer                    :: ncmap_list
    integer                    :: nwater_list
    integer                    :: nsolute_list
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

  ! parameters
  integer,           parameter :: DINumOpened = 100

  ! structures
  type s_di_opened
    integer                    :: file   = 0
    integer                    :: domain = 0
    integer(8)                 :: pos    = 0
    type(s_di_opened), pointer :: next   => null()
  end type s_di_opened

  ! variables
  type(s_domain_index), target :: domain_index_list(27)

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

  integer(8)                   :: di_fast  = 0
  integer(8)                   :: di_slow  = 0
  integer(8)                   :: di_open  = 0
  integer(8)                   :: di_close = 0

  type(s_di_opened), target,save  :: di_opened(DINumOpened+1)

  ! subroutines
  private :: di_get_domain_contains
  private :: di_check_natom_per_domain
  private :: di_read_index
  private :: di_write_index
  private :: di_write_head
  private :: di_write_nindex
  private :: di_write_nindex_
  private :: di_open_file
  private :: di_close_file
  private :: di_get_filename
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
    integer                  :: i, j, iatm, idom


    ! setup opend files linked-list
    !
    do i = 1, DINumOpened + 1
      di_opened(i)%file   = 0
      di_opened(i)%domain = 0
      if (i == DINumOpened + 1) then
        di_opened(i)%next => null()
      else
        di_opened(i)%next => di_opened(i+1)
      end if
    end do


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

    if (.not. main_rank) &
      return

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


    ! Check for number of atom inside domain (for DEBUG)
    !
    !call di_check_natom_per_domain


    ! atom_list
    !
    write(MsgOut,'(a)') 'Di_Create> Atom-list: '

    call di_write_head(.true.)
    do i = 1, hm_num_atoms
      idom = di_get_domain_contains(i)
      call di_write_index(idom, i)
      call di_show_progress(i, hm_num_atoms, 'atom')
    end do

    write(MsgOut,'(a)') '  done.'

    ! bond_list
    !
    write(MsgOut,'(a)') 'Di_Create> Bond-list: '

    call di_write_head(.false.)
    do i = 1, hm_num_bonds
      iatm = hm_bond_list(1, i)
      idom = di_get_domain_contains(iatm)
      call di_write_index(idom, i)
      call di_show_progress(i, hm_num_bonds, 'bond')
    end do

    write(MsgOut,'(a)') '  done.'

    ! angl_list
    !
    write(MsgOut,'(a)') 'Di_Create> Angl-list: '

    call di_write_head(.false.)
    do i = 1, hm_num_angles
      iatm = hm_angl_list(1, i)
      idom = di_get_domain_contains(iatm)
      call di_write_index(idom, i)
      call di_show_progress(i, hm_num_angles, 'angl')
    end do

    write(MsgOut,'(a)') '  done.'


    ! dihe_list
    !
    write(MsgOut,'(a)') 'Di_Create> Dihe-list: '

    call di_write_head(.false.)
    do i = 1, hm_num_dihedrals
      iatm = hm_dihe_list(1, i)
      idom = di_get_domain_contains(iatm)
      call di_write_index(idom, i)
      call di_show_progress(i, hm_num_dihedrals, 'dihe')
    end do

    write(MsgOut,'(a)') '  done.'


    ! impr_list
    !
    write(MsgOut,'(a)') 'Di_Create> Impr-list: '

    call di_write_head(.false.)
    do i = 1, hm_num_impropers
      iatm = hm_impr_list(1, i)
      idom = di_get_domain_contains(iatm)
      call di_write_index(idom, i)
      call di_show_progress(i, hm_num_impropers, 'impr')
    end do

    write(MsgOut,'(a)') '  done.'


    ! cmap_list
    !
    write(MsgOut,'(a)') 'Di_Create> Cmap-list: '

    call di_write_head(.false.)
    do i = 1, hm_num_cmaps
      iatm = hm_cmap_list(1, i)
      idom = di_get_domain_contains(iatm)
      call di_write_index(idom, i)
      call di_show_progress(i, hm_num_cmaps, 'cmap')
    end do

    write(MsgOut,'(a)') '  done.'


    ! water_list
    !
    write(MsgOut,'(a)') 'Di_Create> Water-list: '

    call di_write_head(.false.)
    do i = 1, hm_num_water_list
      iatm = hm_water_list(1, i)
      idom = di_get_domain_contains(iatm)
      call di_write_index(idom, i)
      call di_show_progress(i, hm_num_water_list, 'water_list')
    end do

    write(MsgOut,'(a)') '  done.'


    ! solute_list
    !
    write(MsgOut,'(a)') 'Di_Create> Solute-list: '

    call di_write_head(.false.)
    do i = 1, hm_num_solute_list
      iatm = hm_solute_list(i)
      idom = di_get_domain_contains(iatm)
      call di_write_index(idom, i)
      call di_show_progress(i, hm_num_solute_list, 'solute_list')
    end do

    write(MsgOut,'(a)') '  done.'


    !  selatoms_list
    write(MsgOut,'(a)') 'Di_Create> SelAtoms_list: '

    do j = 1, hm_num_selatoms_list
      call di_write_head(.false.)
      do i = 1, hm_num_selatoms(j)
        iatm = hm_selatoms(i, j)
        idom = di_get_domain_contains(iatm)
        call di_write_index(idom, i)
        call di_show_progress(i, hm_num_solute_list, 'selatoms')
      end do
    end do

    write(MsgOut,'(a)') '  done.'

    call di_write_nindex


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
    integer                  :: i, j
    character(200)           :: filename


    do i = 1, size(domain_index_list)
      
      if (allocated(domain_index_list(i)%iatom_list))   &
         deallocate(domain_index_list(i)%iatom_list)

      if (allocated(domain_index_list(i)%ibond_list))   &
         deallocate(domain_index_list(i)%ibond_list)

      if (allocated(domain_index_list(i)%iangl_list))   &
         deallocate(domain_index_list(i)%iangl_list)

      if (allocated(domain_index_list(i)%idihe_list))   &
         deallocate(domain_index_list(i)%idihe_list)

      if (allocated(domain_index_list(i)%iimpr_list))   &
         deallocate(domain_index_list(i)%iimpr_list)

      if (allocated(domain_index_list(i)%icmap_list))   &
         deallocate(domain_index_list(i)%icmap_list)

      if (allocated(domain_index_list(i)%iwater_list))  &
         deallocate(domain_index_list(i)%iwater_list)

      if (allocated(domain_index_list(i)%isolute_list)) &
         deallocate(domain_index_list(i)%isolute_list)

      if (allocated(domain_index_list(i)%selatoms_list)) then
        do j = 1, size(domain_index_list(i)%selatoms_list)
          if(allocated(domain_index_list(i)%selatoms_list(j)%i)) &
           deallocate (domain_index_list(i)%selatoms_list(j)%i)
        end do
      end if

    end do


    if (di_dnum_x == 0 .or. di_dnum_y == 0 .or. di_dnum_z == 0) &
      call error_msg('Di_Delete> Domain index file is not created.')
    
    if (.not. main_rank) &
      return

    write(filename, '(a,a,a,a)') &
          trim(di_file_dir), '/', &
          trim(di_filename), '_domain_index_file_*'

    call system('rm -f '//trim(filename))


    ! print debug information
    !
    write(MsgOut,'(a,i14)') 'Di_Delete> Fast search: ', di_fast
    write(MsgOut,'(a,i14)') '           Slow search: ', di_slow
    write(MsgOut,'(a,i14)') '           Open file  : ', di_open
    write(MsgOut,'(a,i14)') '           Close file : ', di_close
    write(MsgOut,'(a)')   ' '

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
      call di_read_index(nlist, jdom)

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
      di0 => domain_index_list(ilist)
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
      di0 => domain_index_list(ilist)

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
  !  Subroutine    di_check_natom_per_domain
  !> @brief        check number of atom of per domain
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine di_check_natom_per_domain

    ! local variables
    integer                  :: i, iatm, idom
    integer, allocatable     :: ndom(:)


    allocate(ndom(di_dnum_x * di_dnum_y * di_dnum_z))
    ndom(:) = 0

    do i = 1, hm_num_water_list
      iatm = hm_water_list(1, i)
      idom = di_get_domain_contains(iatm)
      ndom(idom) = ndom(idom) + 3
    end do

    do i = 1, hm_num_solute_list
      iatm = hm_solute_list(i)
      idom = di_get_domain_contains(iatm)
      ndom(idom) = ndom(idom) + 1
    end do

    write(MsgOut,'(a)') 'Di_Check_Natom_Per_Domain> (DEBUG)'
    do i = 1, size(ndom)
      write(MsgOut,'(a,i5,a,i9)') '   Domain: ', i, '  = ', ndom(i)
    end do
    write(MsgOut,'(a)') ''

    deallocate(ndom)

    return

  end subroutine di_check_natom_per_domain

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    di_read_index
  !> @brief        read the domain index list
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine di_read_index(nlist, domain)

    ! formal arguments
    integer,                 intent(in)    :: nlist
    integer,                 intent(in)    :: domain

    ! local variables
    integer                  :: file, i
    character(200)           :: filename

    type(s_domain_index),    pointer :: di


    di => domain_index_list(nlist)

    file = get_unit_no()

    filename = di_get_filename(domain)

#ifndef RICC
    open(file,                   &
         file   = filename,      &
         status = 'old',         &
         form   = 'unformatted', &
         access = 'stream')
#else
#warning
#warning "WARNING: Cannot run prst_setup on RICC platform."
#warning
    call error_msg( &
         'Di_Read_Index> ERROR: Cannot run prst_setup on RICC platform.')
#endif

    ! atom-list
    read(file) di%natom_list
    call di_realloc_memory(di%iatom_list, di%natom_list)
    if (di%natom_list > 0) &
      read(file) di%iatom_list(1:di%natom_list)

    ! bond-list
    read(file) di%nbond_list
    call di_realloc_memory(di%ibond_list, di%nbond_list)
    if (di%nbond_list > 0) &
      read(file) di%ibond_list(1:di%nbond_list)

    ! angl-list
    read(file) di%nangl_list
    call di_realloc_memory(di%iangl_list, di%nangl_list)
    if (di%nangl_list > 0) &
      read(file) di%iangl_list(1:di%nangl_list)

    ! dihe-list
    read(file) di%ndihe_list
    call di_realloc_memory(di%idihe_list, di%ndihe_list)
    if (di%ndihe_list > 0) &
      read(file) di%idihe_list(1:di%ndihe_list)

    ! impr-list
    read(file) di%nimpr_list
    call di_realloc_memory(di%iimpr_list, di%nimpr_list)
    if (di%nimpr_list > 0) &
      read(file) di%iimpr_list(1:di%nimpr_list)

    ! cmap-list
    read(file) di%ncmap_list
    call di_realloc_memory(di%icmap_list, di%ncmap_list)
    if (di%ncmap_list > 0) &
      read(file) di%icmap_list(1:di%ncmap_list)

    ! water-list
    read(file) di%nwater_list
    call di_realloc_memory(di%iwater_list, di%nwater_list)
    if (di%nwater_list > 0) &
      read(file) di%iwater_list(1:di%nwater_list)

    ! solute-list
    read(file) di%nsolute_list
    call di_realloc_memory(di%isolute_list, di%nsolute_list)
    if (di%nsolute_list > 0) &
      read(file,err=10) di%isolute_list(1:di%nsolute_list)

    ! selatoms-list
10  if (hm_num_selatoms_list /= 0 .and. .not. allocated(di%selatoms_list)) &
      allocate(di%selatoms_list(hm_num_selatoms_list))
    do i = 1, hm_num_selatoms_list
      read(file) di%selatoms_list(i)%n
      call di_realloc_memory(di%selatoms_list(i)%i, di%selatoms_list(i)%n)
      if (di%selatoms_list(i)%n > 0) &
        read(file,err=20) di%selatoms_list(i)%i(1:di%selatoms_list(i)%n)
    end do
20  continue

    call close_file(file)

    return

  end subroutine di_read_index

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    di_write_index
  !> @brief        write the domain index list
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine di_write_index(domain, index)

    ! formal arguments
    integer,                 intent(in)    :: domain
    integer,                 intent(in)    :: index

    ! local arguments
    type(s_di_opened),       pointer :: opened, opened_prv


    opened => di_opened(1)%next

    ! check head
    !
    if (opened%domain == domain) then

      write(opened%file,rec=opened%pos) index
      opened%pos = opened%pos + 1

      di_fast = di_fast + 1
      return

    else if (opened%file == 0) then

      call di_open_file(domain, opened)

      write(opened%file,rec=opened%pos) index
      opened%pos = opened%pos + 1

      di_fast = di_fast + 1
      return

    end if

    ! check tails
    !

    do while(.true.)

      opened_prv => opened
      opened     => opened%next

      if (opened%domain == domain) then

        write(opened%file,rec=opened%pos) index
        opened%pos = opened%pos + 1

        di_slow = di_slow + 1

        opened_prv%next   => opened%next
        opened%next       => di_opened(1)%next
        di_opened(1)%next => opened

        return

      else if (opened%file == 0) then

        call di_open_file(domain, opened)

        write(opened%file,rec=opened%pos) index
        opened%pos = opened%pos + 1

        di_slow = di_slow + 1

        opened_prv%next   => opened%next
        opened%next       => di_opened(1)%next
        di_opened(1)%next => opened

        return

      end if

      if (.not. associated(opened%next)) &
        exit

    end do

    ! check tip
    !

    call di_close_file(opened)

    call di_open_file (domain, opened)
    write(opened%file,rec=opened%pos) index
    opened%pos = opened%pos + 1

    di_slow = di_slow + 1

    opened_prv%next   => opened%next
    opened%next       => di_opened(1)%next
    di_opened(1)%next => opened

    return

  end subroutine di_write_index

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    di_write_head
  !> @brief        write header informaiton
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine di_write_head(at_first)

    ! formal arguments
    logical,                   intent(in)    :: at_first

    ! logical variables
    integer(8)                 :: ftell
    integer                    :: idom, ndom, file
    character(200)             :: filename
    logical                    :: lexist

#ifdef KCOMP
    integer                    :: flen
#endif

    type(s_di_opened), pointer :: opened


    !write(MsgOut,'(a)') '  Di_write_head> '

    ! close opened files
    !

    opened => di_opened(1)%next

    do while(associated(opened))
      if (opened%file /= 0) &
        call di_close_file(opened)
      opened => opened%next
    end do


    ! write head sign
    !
    ndom = di_dnum_x * di_dnum_y * di_dnum_z

    do idom = 1, ndom

      filename = di_get_filename(idom)

      ! check file exist
      if (at_first) then
        inquire(file=filename, exist=lexist)
        if (lexist) goto 910
      end if

      ! write head
      file = get_unit_no()
      open(file,                   &
           file   = filename,      &
           form   = 'unformatted', &
           access = 'direct',      &
           recl   = 4,             &
           err    = 900)

#ifdef KCOMP
      inquire(file,flen=flen)
      write(file,rec=(flen / 4 + 1)) 0
#else
      call fseek(file, 0, 2)
      write(file,rec=(ftell(file) / 4 + 1)) 0
#endif

      close(file)
      call free_unit_no(file)

    end do

    return

900 call error_msg('Di_Write_Head> File open error. '//trim(filename))
910 call error_msg('Di_Write_Head> File is already exist. '//trim(filename))

  end subroutine di_write_head

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    di_write_nindex
  !> @brief        write number of index
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine di_write_nindex

    ! logical variables
    integer                    :: idom, ndom

    type(s_di_opened), pointer :: opened


    write(MsgOut,'(a)') '  Di_write_nindex> '

    ! close opened files
    !

    opened => di_opened(1)%next

    do while(associated(opened))
      if (opened%file /= 0) &
        call di_close_file(opened)
      opened => opened%next
    end do


    ! write nindex
    !
    ndom = di_dnum_x * di_dnum_y * di_dnum_z

    do idom = 1, ndom

      call di_write_nindex_(idom)

    end do

    return

  end subroutine di_write_nindex

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    di_write_nindex_
  !> @brief        
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine di_write_nindex_(idom)

    ! formal arguments
    integer,                 intent(in)    :: idom

    ! logical variables
    integer                  :: file, fpos, npos, ierr, ival, nval
    character(200)           :: filename


    file = get_unit_no()

    filename = di_get_filename(idom)

    open(file,                   &
         file   = filename,      &
         status = 'old',         &
         form   = 'unformatted', &
         access = 'direct',      &
         recl   = 4)

    nval = 0
    fpos = 1

    do while(.true.)

      read(file,rec=fpos,iostat=ierr) ival
      if (ierr /= 0) &
           exit

      if (ival == 0) then

        if (nval /= 0) then
          npos = fpos - nval - 1
          write(file,rec=npos) nval
        end if

        nval = 0

      else
        nval = nval + 1

      end if

      fpos = fpos + 1

    end do

    if (nval /= 0) then
      npos = fpos - nval - 1
      write(file,rec=npos) nval
    end if

    close(file)

    call free_unit_no(file)

    return

900 call error_msg('Di_Write_NIndex_> File open error. '//trim(filename))

  end subroutine di_write_nindex_

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    di_open_file
  !> @brief        open the domain index file
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine di_open_file(domain, opened)

    ! formal arguments
    integer,                 intent(in)    :: domain
    type(s_di_opened),       intent(inout) :: opened

    ! local variables
    integer(8)               :: ftell
    character(200)           :: filename

#ifdef KCOMP
    integer                  :: flen
#endif


    filename = di_get_filename(domain)

    opened%file   = get_unit_no()
    opened%domain = domain

    open(opened%file,            &
         file   = filename,      &
         status = 'old',         &
         form   = 'unformatted', &
         access = 'direct',      &
         recl   = 4,             &
         err    = 900)

#ifdef KCOMP
    inquire(opened%file,flen=flen)
    opened%pos = flen / 4 + 1
#else
    call fseek(opened%file, 0, 2)
    opened%pos = ftell(opened%file) / 4 + 1
#endif

    di_open = di_open + 1

    return

900 call error_msg('Di_Open_File> File open error. '//trim(filename))

  end subroutine di_open_file

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    di_close_file
  !> @brief        close the domain index file
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine di_close_file(opened)

    ! formal arguments
    type(s_di_opened),       intent(inout) :: opened


    if (opened%file == 0) &
      call error_msg('Di_Close_File> Unexpected error.')

    close(opened%file)
    call free_unit_no(opened%file)

    di_close = di_close + 1

    opened%file   = 0
    opened%domain = 0

    return

  end subroutine di_close_file

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      di_get_filename
  !> @brief        get the domain index filename
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function di_get_filename(domain)

    ! return value
    character(200)           :: di_get_filename

    ! formal arguments
    integer,                 intent(in)    :: domain


    write(di_get_filename, '(a,a,a,a,i6.6)') &
          trim(di_file_dir), '/', &
          trim(di_filename), '_domain_index_file_', domain

    return

  end function di_get_filename

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


    if (mod(ip, Unit) /= 0) &
      return

    write(MsgOut,'(a,i9,a,i9,a)') '   '//trim(str)//': ', &
                           ip / Unit, &
                           '/',       &
                           nt / Unit, ' million'

    return

  end subroutine di_show_progress

end module pr_domain_index_mod
