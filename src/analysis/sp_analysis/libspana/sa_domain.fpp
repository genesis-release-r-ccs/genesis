!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sa_domain_mod
!> @brief   utilities for domain decomposition
!! @authors Jaewoon Jung (JJ), Takaharu Mori (TM), Isseki Yu (IY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module sa_domain_mod

  use sa_boundary_mod
  use sa_domain_str_mod
  use sa_boundary_str_mod
  use sa_ensemble_str_mod
  use sa_option_str_mod
  use trajectory_str_mod
  use select_atoms_mod
  use select_atoms_str_mod
  use molecules_str_mod
  use timers_mod
  use random_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
#ifdef HAVE_MPI_GENESIS
#ifdef MSMPI
!GCC$ ATTRIBUTES DLLIMPORT :: MPI_BOTTOM, MPI_IN_PLACE
#endif
#endif
  private

  ! parameters for setup cell capacity
  real(wp),         private, parameter :: VolumeBox8 = 512.0_wp
  real(wp),         private, parameter :: ShrinkRate = 1.3_wp
! integer,          private, parameter :: NAtomBox8  = 100
  integer,          private, parameter :: NAtomBox8  = 80
  integer,          private, parameter :: NBondBox8  = 90
  integer,          private, parameter :: NAnglBox8  = 215
  integer,          private, parameter :: NDiheBox8  = 600
  integer,          private, parameter :: NImprBox8  = 200
  integer,          private, parameter :: NCmapBox8  = 12
  integer,          private, parameter :: NHGrpBox8  = 30
  integer,          private, parameter :: NHMovBox8  = 5
  integer,          private, parameter :: NContBox8  = 300

  ! subroutines
  public  :: setup_domain
  public  :: setup_atom
  private :: setup_processor_rank
  private :: setup_cell_capacity
  private :: setup_cell_boundary
  private :: molecule_to_domain

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_domain
  !> @brief        setup domain information
  !! @authors      JJ, IY
  !! @param[in]    boundary    : boundary condition information
  !! @param[in]    molecule    : molecule information
  !! @param[inout] constraints : constraints information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_domain(boundary, molecule, domain, option)

    type(s_boundary),        intent(inout) :: boundary
    type(s_molecule),        intent(in)    :: molecule
    type(s_domain),          intent(inout) :: domain
    type(s_option),          intent(in)    :: option

    ! local variables
    integer                  :: i, j, k, cell(3)
    integer                  :: icel_local, icel
    integer                  :: ncel_local, ncel_bound, ncel_all


    ! initialize structure informations
    !

    call init_domain(boundary, domain, option) !defined in sa_domain_str

    domain%num_atom_all       = molecule%num_atoms

    call setup_processor_rank(option, boundary, domain, cell)

    ! decide cell capacity (max**) for memory allocation
    !
    call setup_cell_capacity(boundary, domain)

    ! memory allocaltion of maps connecting local to global cell indices
    !
    ncel_local = domain%num_cell_local
    ncel_bound = domain%num_cell_boundary
    ncel_all   = ncel_local + ncel_bound


    call alloc_domain(domain, option, DomainCellGlobal, cell(1),cell(2),cell(3))
    call alloc_domain(domain, option, DomainCellLocal,    ncel_local, 1, 1)
    call alloc_domain(domain, option, DomainCellLocBou,   ncel_all,   1, 1)
    call alloc_domain(domain, option, DomainCellBoundary, ncel_bound, 1, 1)

    ! assign global<->local mapping of cell indexa
    ! i, j, k give global position of the cell  (IY)
    icel_local = 0
    do i = domain%cell_start(3), domain%cell_end(3)
      do j = domain%cell_start(2), domain%cell_end(2)
        do k = domain%cell_start(1), domain%cell_end(1)
          icel_local = icel_local + 1
          icel = k + (j-1)*cell(1) + (i-1)*cell(1)*cell(2)
          domain%cell_g2l(icel) = icel_local
          domain%cell_l2g(icel_local) = icel
          domain%cell_l2gx(icel_local) = k
          domain%cell_l2gy(icel_local) = j
          domain%cell_l2gz(icel_local) = i
          domain%cell_l2gx_orig(icel_local) = k
          domain%cell_l2gy_orig(icel_local) = j
          domain%cell_l2gz_orig(icel_local) = i
          domain%cell_gxyz2l(k,j,i) = icel_local
        end do
      end do
    end do

    ! assigin each boundary cell
    ! determin the boundary cells around each domain (IY)
    call setup_cell_boundary(option, cell, boundary%num_domain, domain)

#ifdef DEBUG
    ! debug
    !

    if (main_rank) then
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_domain_str     ::MaxAtom      : ', MaxAtom
      write(MsgOut,*) 'sp_domain_str     ::MaxWater     : ', MaxWater
      write(MsgOut,*) 'sp_domain_str     ::MaxMove      : ', MaxMove
      write(MsgOut,*) 'sp_domain_str     ::MaxWaterMove : ', MaxWaterMove
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_enefunc_str    ::MaxNb15      : ', MaxNb15
      write(MsgOut,*) 'sp_enefunc_str    ::MaxNb15Water : ', MaxNb15Water
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_constraints_str::HGroupMax    : ', HGroupMax
      write(MsgOut,*) 'sp_constraints_str::HGrpMaxMove  : ', HGrpMaxMove
      write(MsgOut,*) ''

    end if
#endif

    return

  end subroutine setup_domain

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_atom
  !> @brief        setup atom maps with whole atoms
  !! @authors      JJ, IY
  !! @param[in]    option     : option section information
  !! @param[in]    ensemble   : ensemble information
  !! @param[in]    trajectory : trajectory information
  !! @param[in]    selec_atom : atom selection information
  !! @param[inout] molecule   : molecule information
  !! @param[inout] boundary   : boundary condition information
  !! @param[inout] domain     : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_atom(option, ensemble, trajectory, selec_atom, &
                        molecule, boundary, domain)

    ! formal arguments
    type(s_option),           intent(in)    :: option
    type(s_ensemble),         intent(in)    :: ensemble
    type(s_trajectory),       intent(in)    :: trajectory
    type(s_parray),           intent(in)    :: selec_atom(:)
    type(s_molecule), target, intent(inout) :: molecule
    type(s_boundary), target, intent(inout) :: boundary
    type(s_domain),   target, intent(inout) :: domain

    ! local variable
    real(wp)                  :: x_shift, y_shift, z_shift
    real(wp)                  :: move(3), origin(3)
    integer                   :: step
    integer                   :: i, ig, icx, icy, icz, icel, is
    integer                   :: igroup
    integer                   :: icel_local
    integer                   :: ncel_local, ncel
    integer                   :: natom_all, natom_selec
    integer                   :: ngroup
    integer                   :: alloc_stat
    integer                   :: dealloc_stat
    integer                   :: inside_x, inside_y, inside_z
    integer                   :: insidebox

    integer,      allocatable :: natom_group(:)
    real(wp),         pointer :: bsize_x, bsize_y, bsize_z
    real(wp),         pointer :: csize_x, csize_y, csize_z
    integer,          pointer :: celnatom_group(:,:)
    integer,          pointer :: ncel_x, ncel_y, ncel_z
    integer,          pointer :: cell_g2l(:), cell_g2b(:)
    integer,          pointer :: natom(:)


    bsize_x      => boundary%box_size_x
    bsize_y      => boundary%box_size_y
    bsize_z      => boundary%box_size_z
    csize_x      => boundary%cell_size_x
    csize_y      => boundary%cell_size_y
    csize_z      => boundary%cell_size_z
    ncel_x       => boundary%num_cells_x
    ncel_y       => boundary%num_cells_y
    ncel_z       => boundary%num_cells_z

    cell_g2l     => domain%cell_g2l
    cell_g2b     => domain%cell_g2b

    origin(1)    = boundary%origin_x
    origin(2)    = boundary%origin_y
    origin(3)    = boundary%origin_z

    ncel         = domain%num_cell_local + domain%num_cell_boundary
    ncel_local   = domain%num_cell_local
    natom_all    = domain%num_atom_all
    natom_selec  = domain%num_atom_selec
    ngroup       = domain%num_group

    alloc_stat   = 0
    dealloc_stat = 0


    call alloc_domain(domain, option, DomainDynvar, ncel, ngroup, 1)

    ! these should be defined after allocation
    natom           => domain%num_atom
    celnatom_group  => domain%num_atom_group


    ! (IY)
    ! allocation size corresponds  number of selected atom (not for all atoms)
    ! id_s2l, selec_atom2cell are alocasted and initialized with zero
    call alloc_domain(domain, option, DomainGlobal, domain%num_atom_all, 1, 1)
    call alloc_domain(domain, option, DomainSelec, domain%num_atom_selec, 1, 1)

    !determin the natom for each selec_atom group
    !
    if( .not. allocated(natom_group)) then
      allocate(natom_group(ngroup), stat = alloc_stat)
    end if

    do igroup  = 1, ngroup
      natom_group(igroup) = size(selec_atom(igroup)%idx)
    end do

    !mapping serial id of selec_atom to global atom id (only once)
    is = 0
    if (domain%id_selec2global(1) == 0) then
      do igroup = 1, ngroup
        do i = 1, natom_group(igroup)
          is = is +1
          ig = selec_atom(igroup)%idx(i)
          domain%id_selec2global(is)=ig
          domain%id_global2selec(ig)=is
        end do
      end do
    end if

    do step = 1, 2

      natom(1:ncel) = 0
      celnatom_group(1:ngroup, 1:ncel) = 0
      domain%num_atom_t0(1:ncel)    = 0

      ! id_l2s, id_g2l,  should be initialized
      ! in every timestep (IY)
      if (step == 2) then
        domain%id_l2s(1:MaxAtom, 1:ncel)  = 0
        domain%id_s2l(1:2,1:domain%num_atom_selec) = 0
      end if

      is = 0
      do igroup = 1, ngroup
      do i = 1, natom_group(igroup)

        insidebox = 1
        inside_x  = 1
        inside_y  = 1
        inside_z  = 1

        is = is +1

        ! global id of atom i
        ig = selec_atom(igroup)%idx(i)

        x_shift = molecule%atom_coord(1,ig)
        y_shift = molecule%atom_coord(2,ig)
        z_shift = molecule%atom_coord(3,ig)

        ! judge if the atom is inside the box or not
        ! this is needed 1: when the trajectory is not wrapped
        ! or 2: when trajectory is wrapped by PBC box, but
        ! box size is manually determined by ctrl, and is smaller than
        ! PBC box size
        if (x_shift > bsize_x .or. x_shift < 0) &
          inside_x = 0
        if (y_shift > bsize_y .or. y_shift < 0) &
          inside_y = 0
        if (z_shift > bsize_z .or. z_shift < 0) &
          inside_z = 0
        if(inside_x*inside_y*inside_z == 0) &
          insidebox = 0

        !assign which cell
        !
        if (insidebox == 1) then
          icx = int(x_shift/csize_x)
          icy = int(y_shift/csize_y)
          icz = int(z_shift/csize_z)
          if (icx == ncel_x) icx = icx - 1
          if (icy == ncel_y) icy = icy - 1
          if (icz == ncel_z) icz = icz - 1
          icel = 1 + icx + icy*ncel_x + icz*ncel_x*ncel_y
          domain%selec_atom2cell(is) = icel ! isseki

          ! atoms inside the domain
          !
          if (cell_g2l(icel) /= 0 ) then

            ! local cell index
            !
            icel_local = cell_g2l(icel)
            natom(icel_local)   = natom(icel_local)   + 1
            celnatom_group(igroup, icel_local) = &
            celnatom_group(igroup, icel_local) + 1

            if (step == 2) &
              call molecule_to_domain(molecule, move, origin, is, &
                                      domain, icel_local, natom(icel_local))

          ! atoms in a boundary
          !
          else if (cell_g2b(icel) /= 0) then

            ! local cell index
            !
            icel_local = cell_g2b(icel) + ncel_local
            natom(icel_local)   = natom(icel_local)   + 1
            celnatom_group(igroup, icel_local) = &
            celnatom_group(igroup, icel_local) + 1

            if (step == 2) &
              call molecule_to_domain(molecule, move, origin, is, &
                                      domain, icel_local, natom(icel_local))
          end if
        end if ! insidebox

      end do ! i(atom)
      end do ! igroup

      if (step == 1) then

        MaxAtom = 0
        do i = 1, ncel
          MaxAtom = max(MaxAtom, natom(i))
        end do

#ifdef HAVE_MPI_GENESIS
        call mpi_allreduce(mpi_in_place, MaxAtom, 1, mpi_integer, mpi_max, &
                           mpi_comm_country, ierror)
#endif
        !MaxAtom   = MaxAtom   * 3/2
        call alloc_domain(domain, option, DomainDynvar_Atom, ncel, 1, 1)

      end if

    end do ! step

    domain%num_atom_t0(1:ncel) = natom(1:ncel)

    return

  end subroutine setup_atom

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_processor_rank
  !> @brief        define the processor rank in each dimension and decide the
  !!               number of cells in each domain
  !! @authors      JJ
  !! @param[in]    option   : option information
  !! @param[in]    boundary : boundary condition information
  !! @param[inout] domain   : domain information
  !! @param[out]   cell     : cells in boundary
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_processor_rank(option, boundary, domain, cell)

    ! formal arguments
    type(s_option),          intent(in)    :: option
    type(s_boundary),target, intent(in)    :: boundary
    type(s_domain),  target, intent(inout) :: domain
    integer,                 intent(inout) :: cell(3)

    ! local variables
    integer                  :: proc_rank, iproc(3)
    integer                  :: i, j, k, quotient, remainder
    integer                  :: mbl

    integer, pointer         :: ncell_local, ncell_boundary, ndomain(:)
    integer, pointer         :: cell_start(:), cell_end(:), cell_length(:)
    integer, pointer         :: iproc_lower(:), iproc_upper(:), neighbor(:,:,:)


    ndomain        => boundary%num_domain
    ncell_local    => domain%num_cell_local
    ncell_boundary => domain%num_cell_boundary
    cell_start     => domain%cell_start
    cell_end       => domain%cell_end
    cell_length    => domain%cell_length
    iproc_lower    => domain%iproc_lower
    iproc_upper    => domain%iproc_upper
    neighbor       => domain%neighbor

    mbl            = domain%max_boundary_layer

#ifdef DEBUG
    if (main_rank) write(MsgOut,'(A),1x,i0') "max_boundary_layer=", mbl
#endif

    ! Assign the rank of each dimension from my_rank
    ! proc_rank = rank of process.
    !
    proc_rank = my_city_rank

    iproc(1) = mod(proc_rank, ndomain(1))
    iproc(2) = mod(proc_rank/ndomain(1), ndomain(2))
    iproc(3) = proc_rank/(ndomain(1)*ndomain(2))

    ! Cell number of each dimension
    !
    cell(1) = boundary%num_cells_x
    cell(2) = boundary%num_cells_y
    cell(3) = boundary%num_cells_z

    ! Default value of the number of cell in each domain
    !
    ncell_local = 1

    ! Assign the cell index for each processor and the total number of cells
    !
    do i = 1, 3

      quotient = cell(i) / ndomain(i)
      remainder = mod(cell(i), ndomain(i))

      if (iproc(i) <= (remainder -1)) then
        quotient       = quotient + 1
        cell_start(i)  = quotient * iproc(i) + 1
        cell_end(i)    = cell_start(i) + quotient - 1
        cell_length(i) = quotient
        ncell_local    = ncell_local * quotient
      else
        cell_start(i)  = (quotient+1)*remainder + quotient*(iproc(i)-remainder)
        cell_start(i)  = cell_start(i) + 1
        cell_end(i)    = cell_start(i) + quotient - 1
        cell_length(i) = quotient
        ncell_local    = ncell_local * quotient
      end if

    end do

    ! Assign the lower processor index
    !
    do i = 1, 3

      iproc(i) = iproc(i) - 1
      if (iproc(i) == -1) &
        iproc(i) = ndomain(i) - 1

      iproc_lower(i) = iproc(1) + &
                       iproc(2)*ndomain(1) + &
                       iproc(3)*ndomain(1)*ndomain(2)

      iproc(i) = iproc(i) + 1
      if (iproc(i) == ndomain(i)) &
        iproc(i) = 0

    end do

    ! Assign the upper processor index
    !
    do i = 1, 3

      iproc(i) = iproc(i) + 1
      if (iproc(i) == ndomain(i)) &
        iproc(i) = 0

      iproc_upper(i) = iproc(1) + &
                       iproc(2)*ndomain(1) + &
                       iproc(3)*ndomain(1)*ndomain(2)

      iproc(i) = iproc(i) - 1
      if (iproc(i) == -1) &
        iproc(i) = ndomain(i) - 1
    end do

    ! Assign the neighboring process index
    !
    do i = -1, 1

      if (i == -1) then

        iproc(1) = iproc(1) - 1
        if (iproc(1) == -1) &
          iproc(1) = ndomain(1) - 1

      else if (i == 1) then

        iproc(1) = iproc(1) + 1
        if (iproc(1) == ndomain(1)) &
          iproc(1) = 0

      end if

      do j = -1, 1

        if (j == -1) then

          iproc(2) = iproc(2) - 1
          if (iproc(2) == -1) &
            iproc(2) = ndomain(2) - 1

        else if (j == 1) then

          iproc(2) = iproc(2) + 1
          if (iproc(2) == ndomain(2)) &
            iproc(2) = 0

        end if

        do k = -1, 1

          if (k == -1) then

            iproc(3) = iproc(3) - 1
            if (iproc(3) == -1) &
              iproc(3) = ndomain(3) - 1

          else if (k == 1) then

            iproc(3) = iproc(3) + 1
            if (iproc(3) == ndomain(3)) &
              iproc(3) = 0
          end if

          neighbor(i,j,k) = iproc(1) + &
                            iproc(2)*ndomain(1) + &
                            iproc(3)*ndomain(1)*ndomain(2)
          if (real_calc) &
           neighbor(i,j,k) = neighbor(i,j,k)
          if (k == -1) then

            iproc(3) = iproc(3) + 1
            if (iproc(3) == ndomain(3)) &
              iproc(3) = 0

          else if (k == 1) then

            iproc(3) = iproc(3) - 1
            if (iproc(3) == -1) &
              iproc(3) = ndomain(3) - 1

          end if
        end do

        if (j == -1) then

          iproc(2) = iproc(2) + 1
          if (iproc(2) == ndomain(2)) &
            iproc(2) = 0

        else if (j == 1) then

          iproc(2) = iproc(2) - 1
          if (iproc(2) == -1) &
            iproc(2) = ndomain(2) - 1
        end if

      end do

      if (i == -1) then

        iproc(1) = iproc(1) + 1
        if (iproc(1) == ndomain(1)) &
          iproc(1) = 0

      else if (i == 1) then

        iproc(1) = iproc(1) - 1
        if (iproc(1) == -1) &
          iproc(1) = ndomain(1) - 1

      end if
    end do

    ! chenged calculation formula (IY)
    ncell_boundary = (cell_length(1)+2*mbl)*(cell_length(2)+2*mbl)*&
                     (cell_length(3)+2*mbl)-cell_length(1)*cell_length(2)&
                                                          *cell_length(3)

    return

  end subroutine setup_processor_rank

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_cell_capacity
  !> @brief        setup cell capacity for memory allocation
  !! @authors      JJ
  !! @param[in]    boundary : boundary information
  !! @param[in]    domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_cell_capacity(boundary, domain)

    ! formal arguments
    type(s_boundary),        intent(in)    :: boundary
    type(s_domain),          intent(in)    :: domain

    ! local variables
    real(wp)                 :: v, v_rate


    v = boundary%box_size_x * &
        boundary%box_size_y * &
        boundary%box_size_z

    v = v / real(boundary%num_domain(1)* &
                 boundary%num_domain(2)* &
                 boundary%num_domain(3),wp)
    v = v / real(domain%num_cell_local,wp)

    v_rate = v / VolumeBox8


    MaxContact   = int(v_rate * real(NContBox8,wp) * ShrinkRate)

    ContactMove  = MaxContact / 2

    ! sp_domain_str
    !

    MaxAtom      = int(v_rate * real(NAtomBox8,wp) * ShrinkRate)
    MaxSolute    = MaxAtom
    MaxWater     = int(v_rate * real(NAtomBox8,wp) / 3.0_wp * ShrinkRate)
    MaxMove      = int(v_rate * real(NAtomBox8,wp) / 5.0_wp * ShrinkRate)
    MaxWaterMove = int(v_rate * real(NAtomBox8,wp) / 7.0_wp * ShrinkRate)


    ! sp_enefunc_str
    !

    MaxBond      = int(v_rate * real(NBondBox8,wp) * ShrinkRate)
    MaxAngle     = int(v_rate * real(NAnglBox8,wp) * ShrinkRate)
    MaxDihe      = int(v_rate * real(NDiheBox8,wp) * ShrinkRate)
    MaxImpr      = int(v_rate * real(NImprBox8,wp) * ShrinkRate)
    MaxCmap      = int(v_rate * real(NCmapBox8,wp) * ShrinkRate)

    BondMove     = MaxBond  / 2
    AngleMove    = MaxAngle / 2
    DiheMove     = MaxDihe  / 2
    ImprMove     = MaxImpr  / 2

    ! sp_pairlist_str
    !

    MaxNb15      = int(v_rate * real(NAtomBox8,wp) * ShrinkRate)
    MaxNb15      = MaxNb15 ** 2

    MaxNb15Water = int(v_rate * real(NAtomBox8,wp) * ShrinkRate)
    MaxNb15Water = MaxNb15Water ** 2
    MaxNb15water = MaxNb15Water / 10

    ! sp_constraints_str
    !

    HGroupMax    = int(v_rate * real(NHGrpBox8,wp) * ShrinkRate)
    HGrpMaxMove  = int(v_rate * real(NHMovBox8,wp) * ShrinkRate)

#ifdef DEBUG
    ! debug
    !

    if (main_rank) then
      write(MsgOut,*) 'Cell volume                      : ', v
      write(MsgOut,*) ''

    end if
#endif

    return

  end subroutine setup_cell_capacity

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_cell_boundary
  !> @brief        define boundary cells in each domain
  !! @authors      JJ, IY
  !! @param[in]    option : option information
  !! @param[in]    cell   : cell count for each axis
  !! @param[in]    ndom   : domain count for each axis
  !! @param[inout] domain : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_cell_boundary(option, cell, ndom, domain)

    ! formal arguments
    type(s_option),          intent(in)    :: option
    integer,                 intent(in)    :: cell(3)
    integer,                 intent(in)    :: ndom(3)
    type(s_domain),  target, intent(inout) :: domain

    ! local variables
    integer                  :: i, j, k, ic, jc, kc, icel, icel_local
    integer                  :: iw, jw, kw
    integer                  :: ibl, mbl

    integer,         pointer :: ncell, ncell_boundary
    integer,         pointer :: cell_start(:), cell_end(:)
    integer,         pointer :: cell_g2b(:), cell_b2g(:)
    integer,         pointer :: cell_gxyz2l(:,:,:)
    integer,         pointer :: cell_l2gx(:), cell_l2gy(:), cell_l2gz(:)
    integer,         pointer :: cell_l2gx_orig(:)
    integer,         pointer :: cell_l2gy_orig(:)
    integer,         pointer :: cell_l2gz_orig(:)


    ncell          => domain%num_cell_local
    ncell_boundary => domain%num_cell_boundary
    cell_start     => domain%cell_start
    cell_end       => domain%cell_end
    cell_g2b       => domain%cell_g2b
    cell_b2g       => domain%cell_b2g
    cell_gxyz2l    => domain%cell_gxyz2l
    cell_l2gx      => domain%cell_l2gx
    cell_l2gy      => domain%cell_l2gy
    cell_l2gz      => domain%cell_l2gz
    cell_l2gx_orig => domain%cell_l2gx_orig
    cell_l2gy_orig => domain%cell_l2gy_orig
    cell_l2gz_orig => domain%cell_l2gz_orig

    mbl            =  domain%max_boundary_layer

    icel_local = 0

    !-------- (IY)----------
    ! ibl means ith layer of boundary
    ! iw, jw, kw are wrapped position
    ! ic, jc, kc are absolute position
    do ibl = 1, mbl

    if (ndom(1) > 1) then
      ! face boundary (x direction, upper)
      !
      do j = cell_start(2)-ibl+1, cell_end(2)+ibl-1
        ! absolute position
        jc = j
        ! wrapped position
        jw = j

        if (jc < 1) then
          jw  = cell(2) + jc
        else if ( jc > cell(2)) then
          jw  = jc - cell (2)
        end if

        do k = cell_start(3)-ibl+1, cell_end(3)+ibl-1
          ! absolute position for k
          kc = k
          ! wrapped position
          kw = k

          if (kc < 1) then
            kw  = cell(3) + kc
          else if ( kc > cell(3)) then
            kw  = kc - cell (3)
          end if


          icel_local = icel_local + 1 ! count ncell

          ic = cell_end(1)+ibl ! absolute position for i
          iw = ic
          if (ic >  cell(1)) then
            iw = ic - cell(1)
          end if

          icel = iw + (jw-1)*cell(1) + (kw-1)*cell(1)*cell(2)
          cell_g2b(icel) = icel_local
          cell_b2g(icel_local) = icel
          cell_l2gx(icel_local+ncell) = ic
          cell_l2gy(icel_local+ncell) = jc
          cell_l2gz(icel_local+ncell) = kc
          cell_l2gx_orig(icel_local+ncell) = iw
          cell_l2gy_orig(icel_local+ncell) = jw
          cell_l2gz_orig(icel_local+ncell) = kw
          cell_gxyz2l(ic,jc,kc) = icel_local+ncell

        end do
      end do


      ! face boundary (x direction, lower)
      !
      do j = cell_start(2)-ibl+1, cell_end(2)+ibl-1
        ! absolute position
        jc = j
        ! wrapped position
        jw = j

        if (jc < 1) then
          jw  = cell(2) + jc
        else if ( jc > cell(2)) then
          jw  = jc - cell (2)
        end if

        do k = cell_start(3)-ibl+1, cell_end(3)+ibl-1
          ! absolute position for k
          kc = k
          ! wrapped position
          kw = kc

          if (kc < 1) then
            kw  = cell(3) + kc
          else if ( kc > cell(3)) then
            kw  = kc - cell (3)
          end if

          icel_local = icel_local + 1

          ic = cell_start(1)-ibl
          iw = ic
          if (ic < 1) then
            iw = cell(1) + ic
          end if

          icel = iw + (jw-1)*cell(1) + (kw-1)*cell(1)*cell(2)
          cell_g2b(icel) = icel_local
          cell_b2g(icel_local) = icel
          cell_l2gx(icel_local+ncell) = ic
          cell_l2gy(icel_local+ncell) = jc
          cell_l2gz(icel_local+ncell) = kc
          cell_l2gx_orig(icel_local+ncell) = iw
          cell_l2gy_orig(icel_local+ncell) = jw
          cell_l2gz_orig(icel_local+ncell) = kw
          cell_gxyz2l(ic,jc,kc) = icel_local+ncell

        end do
      end do
    end if

    if (ndom(2) > 1) then
      ! face boundary (y direction, upper)
      !
      do i = cell_start(1)-ibl, cell_end(1)+ibl
        ic = i
        iw = i
        if (ic < 1) then
          iw = cell(1) + ic
        else if (ic > cell(1)) then
          iw = ic - cell(1)
        end if

        do k = cell_start(3)-ibl+1, cell_end(3)+ibl-1
          kc = k
          kw = k
          if (kc < 1) then
            kw = cell(3) + kc
          else if (kc > cell(3)) then
            kw = kc - cell(3)
          end if

          icel_local = icel_local + 1

          jc = cell_end(2)+ibl
          jw = jc
          if (jc > cell(2)) then
            jw = jc - cell(2)
          end if

          icel = iw + (jw-1)*cell(1) + (kw-1)*cell(1)*cell(2)
          cell_g2b(icel) = icel_local
          cell_b2g(icel_local) = icel
          cell_l2gx(icel_local+ncell) = ic
          cell_l2gy(icel_local+ncell) = jc
          cell_l2gz(icel_local+ncell) = kc
          cell_l2gx_orig(icel_local+ncell) = iw
          cell_l2gy_orig(icel_local+ncell) = jw
          cell_l2gz_orig(icel_local+ncell) = kw
          cell_gxyz2l(ic,jc,kc) = icel_local+ncell

        end do
      end do


      ! face boundary (y direction, lower)
      !
      do i = cell_start(1)-ibl, cell_end(1)+ibl
        ic = i
        iw = i
        if (ic < 1) then
          iw = cell(1) + ic
        else if (ic > cell(1)) then
          iw = ic - cell(1)
        end if

        do k = cell_start(3)-ibl+1, cell_end(3)+ibl-1
          kc = k
          kw = k
          if (kc < 1) then
            kw = cell(3) + kc
          else if (kc > cell(3)) then
            kw = kc - cell(3)
          end if

          icel_local = icel_local + 1

          jc = cell_start(2)-ibl
          jw = jc
          if (jc < 1) then
            jw = cell(2) +jc
          end if

          icel = iw + (jw-1)*cell(1) + (kw-1)*cell(1)*cell(2)
          cell_g2b(icel) = icel_local
          cell_b2g(icel_local) = icel
          cell_l2gx(icel_local+ncell) = ic
          cell_l2gy(icel_local+ncell) = jc
          cell_l2gz(icel_local+ncell) = kc
          cell_l2gx_orig(icel_local+ncell) = iw
          cell_l2gy_orig(icel_local+ncell) = jw
          cell_l2gz_orig(icel_local+ncell) = kw
          cell_gxyz2l(ic,jc,kc) = icel_local+ncell

        end do
      end do
    end if

    if (ndom(3) > 1) then

      ! face boundary (z direction, upper)
      !
      do i = cell_start(1)-ibl, cell_end(1)+ibl
        ic = i
        iw = i
        if (ic < 1) then
          iw = cell(1) + ic
        else if (ic > cell(1)) then
          iw = ic - cell (1)
        end if

        do j = cell_start(2)-ibl, cell_end(2)+ibl
          jc = j
          jw = j
          if (jc < 1) then
            jw = cell(2) + jc
          else if (jc > cell(2)) then
            jw = jc - cell(2)
          end if

          icel_local = icel_local + 1

          kc = cell_end(3) + ibl
          kw = kc
          if (kc > cell(3)) then
            kw = kc - cell(3)
          end if

          icel = iw + (jw-1)*cell(1) + (kw-1)*cell(1)*cell(2)
          cell_g2b(icel) = icel_local
          cell_b2g(icel_local) = icel
          cell_l2gx(icel_local+ncell) = ic
          cell_l2gy(icel_local+ncell) = jc
          cell_l2gz(icel_local+ncell) = kc
          cell_l2gx_orig(icel_local+ncell) = iw
          cell_l2gy_orig(icel_local+ncell) = jw
          cell_l2gz_orig(icel_local+ncell) = kw
          cell_gxyz2l(ic,jc,kc) = icel_local+ncell

        end do
      end do


      ! face boundary (z direction, lower)
      !
      do i = cell_start(1)-ibl, cell_end(1)+ibl
        ic = i
        iw = i
        if (ic < 1) then
          iw = cell(1) + ic
        else if (ic > cell(1)) then
          iw = ic - cell (1)
        end if

        do j = cell_start(2)-ibl, cell_end(2)+ibl
          jc = j
          jw = j
          if (jc < 1) then
            jw = cell(2) + jc
          else if (jc > cell(2)) then
            jw = jc - cell(2)
          end if

          icel_local = icel_local + 1

          kc = cell_start(3)-ibl
          kw = kc
          if (kc < 1) then
            kw = cell(3) + kc
          end if

          icel = iw + (jw-1)*cell(1) + (kw-1)*cell(1)*cell(2)
          cell_g2b(icel) = icel_local
          cell_b2g(icel_local) = icel
          cell_l2gx(icel_local+ncell) = ic
          cell_l2gy(icel_local+ncell) = jc
          cell_l2gz(icel_local+ncell) = kc
          cell_l2gx_orig(icel_local+ncell) = iw
          cell_l2gy_orig(icel_local+ncell) = jw
          cell_l2gz_orig(icel_local+ncell) = kw
          cell_gxyz2l(ic,jc,kc) = icel_local+ncell

        end do
      end do
    end if

    end do   ! roop for boundary layer


    ! total number of boundary cells
    !
    ncell_boundary = icel_local

    return

  end subroutine setup_cell_boundary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    molecule_to_domain
  !> @brief        copy molecule information to domain
  !! @authors      JJ, IY
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine molecule_to_domain(molecule, move, origin, is, &
                                domain, icel, icel_atom)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    real(wp),                 intent(in)    :: move(3)
    real(wp),                 intent(in)    :: origin(3)
    integer,                  intent(in)    :: is
    type(s_domain),   target, intent(inout) :: domain
    integer,                  intent(in)    :: icel
    integer,                  intent(in)    :: icel_atom

    ! local variables
    integer                   :: ig
    real(wp),         pointer :: coord(:,:) 
    real(wp),         pointer :: coord_local (:,:,:)
    integer,          pointer :: id_l2s      (:,:)
    integer,          pointer :: id_s2l      (:,:)


    coord        => molecule%atom_coord
    id_s2l       => domain%id_s2l
    id_l2s       => domain%id_l2s
    coord_local  => domain%coord

    ! ig means  is-th selected atom
    ig = domain%id_selec2global(is)

    id_l2s(icel_atom,icel) = is
    id_s2l(1,is) = icel
    id_s2l(2,is) = icel_atom

    coord_local(1:3,icel_atom,icel) = coord (1:3,ig) - origin(1:3)

    return

  end subroutine molecule_to_domain

end module sa_domain_mod
