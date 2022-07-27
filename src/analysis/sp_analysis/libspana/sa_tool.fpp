!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sa_tool_mod
!> @brief
!! @authors Isseki Yu (IY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module sa_tool_mod

  use sa_boundary_mod
  use sa_boundary_str_mod
  use sa_option_str_mod
  use sa_ensemble_str_mod
  use measure_mod
  use trajectory_str_mod
  use output_str_mod
  use molecules_str_mod
  use fileio_grocrd_mod
  use fileio_grotop_mod
  use fileio_ambcrd_mod
  use fileio_prmtop_mod
  use fileio_localres_mod
  use fileio_psf_mod
  use fileio_rst_mod
  use fileio_str_mod
  use fileio_par_mod
  use fileio_top_mod
  use fileio_crd_mod
  use fileio_pdb_mod
  use fileio_mode_mod
  use fileio_control_mod
  use fileio_mod
  use messages_mod
  use string_mod
  use constants_mod
  use mpi_parallel_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif
#ifdef OMP
  use omp_lib
#endif

  implicit none
  private

  ! structures
  type, public :: s_spana_inp_info
    character(MaxFilename) :: topfile    = ''
    character(MaxFilename) :: parfile    = ''
    character(MaxFilename) :: strfile    = ''
    character(MaxFilename) :: psffile    = ''
    character(MaxFilename) :: prmtopfile = ''
    character(MaxFilename) :: grotopfile = ''
    character(MaxFilename) :: pdbfile    = ''
    character(MaxFilename) :: crdfile    = ''
    character(MaxFilename) :: ambcrdfile = ''
    character(MaxFilename) :: grocrdfile = ''
    character(MaxFilename) :: rstfile    = ''
    character(MaxFilename) :: reffile    = ''
    character(MaxFilename) :: ambreffile = ''
    character(MaxFilename) :: groreffile = ''
    character(MaxFilename) :: local_resfile = ''
    character(MaxFilename) :: modefile   = ''
  end type s_spana_inp_info

  ! subroutines
  public  :: mindist_cell2atom
  public  :: get_cofm
  public  :: shift_molecule
  public  :: wrap_molecule
  public  :: input_parameter
  private :: wrap_atom
  private :: shift_firstq

contains
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    mindist_grid2atom
  !> @brief        compute minimum distance of atoms inside bounds
  !! @authors      IY
  !! @param[in]    molecule   : molecule information
  !! @param[in]    boundary   : boundary information
  !! @param[in]    cell_iposi : cell position in index
  !! @param[in]    atoms_in_bound : atoms in boundary
  !! @param[out]   mindist    : minimum distance
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine mindist_cell2atom(molecule, boundary, &
                               cell_iposi, atoms_in_bound, mindist)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_boundary),        intent(in)    :: boundary
    integer,                 intent(in)    :: cell_iposi(3)
    integer,                 intent(in)    :: atoms_in_bound(:)
    real(wp),                intent(out)   :: mindist

    ! real position of the cell
    real(wp)                 :: cell_px, cell_py, cell_pz
    real(wp)                 :: atom_px, atom_py, atom_pz

    ! size of a cell for each direction
    real(wp)                 :: cell_sx, cell_sy, cell_sz
    real(wp)                 :: bsize_x, bsize_y, bsize_z
    real(wp)                 :: dist
    integer                  :: iatom, iatom_global


    bsize_x = boundary%box_size_x
    bsize_y = boundary%box_size_y
    bsize_z = boundary%box_size_z

    cell_sx = bsize_x/boundary%num_cells_x
    cell_sy = bsize_y/boundary%num_cells_y
    cell_sz = bsize_z/boundary%num_cells_z

    ! position of the cell center
    cell_px = (real(cell_iposi(1),wp)-0.5_wp)*cell_sx
    cell_py = (real(cell_iposi(2),wp)-0.5_wp)*cell_sy
    cell_pz = (real(cell_iposi(3),wp)-0.5_wp)*cell_sz

    mindist = 100000.0_wp
    do iatom = 1, size(atoms_in_bound)
      iatom_global = atoms_in_bound(iatom)
      atom_px      = molecule%atom_coord(1,iatom_global)
      atom_py      = molecule%atom_coord(2,iatom_global)
      atom_pz      = molecule%atom_coord(3,iatom_global)

      atom_px  =  atom_px  - bsize_x*anint((atom_px-cell_px)/bsize_x)
      atom_py  =  atom_py  - bsize_y*anint((atom_py-cell_py)/bsize_y)
      atom_pz  =  atom_pz  - bsize_z*anint((atom_pz-cell_pz)/bsize_z)

      dist = (atom_px - cell_px)* (atom_px - cell_px)+&
             (atom_py - cell_py)* (atom_py - cell_py)+&
             (atom_pz - cell_pz)* (atom_pz - cell_pz)

       if(dist < mindist)then
         mindist = dist
       end if
    end do

    mindist = sqrt(mindist)

    return

  end subroutine mindist_cell2atom

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_cofm
  !> @brief        get center of mass
  !! @authors      IY
  !! @param[in]    molecule     : molecule information
  !! @param[in]    select_atom  : selected atom
  !! @param[in]    target_group : target group
  !! @param[out]   cofm         : center of mass
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine get_cofm(molecule, selec_atom, target_group, cofm)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_parray),          intent(in)    :: selec_atom(:)
    integer,                 intent(in)    :: target_group
    real(wp),                intent(out)   :: cofm(3)

    ! local variables
    integer                  :: i, ig
    integer                  :: natom
    real(wp)                 :: p_x, p_y, p_z
    real(wp)                 :: mass, totmass


    natom = size(selec_atom(target_group)%idx)

    cofm(1:3) = 0
    totmass   = 0

    do i = 1, natom

      ! global id of atom i
      ig = selec_atom(target_group)%idx(i)

      p_x = molecule%atom_coord(1,ig)
      p_y = molecule%atom_coord(2,ig)
      p_z = molecule%atom_coord(3,ig)

      mass = molecule%mass(ig)

      totmass = totmass + mass

      cofm(1) = cofm(1) + mass*p_x
      cofm(2) = cofm(2) + mass*p_y
      cofm(3) = cofm(3) + mass*p_z

    end do

    cofm(1) = cofm(1)/totmass
    cofm(2) = cofm(2)/totmass
    cofm(3) = cofm(3)/totmass

    return

  end subroutine get_cofm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    shift_molecule
  !> @brief        shift molecule coordinates
  !! @authors      IY
  !! @param[in]    move     : move amount
  !! @param[inout] molecule : molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine shift_molecule(move, molecule)

    ! formal arguments
    real(wp),                intent(in)    :: move(3)
    type(s_molecule),        intent(inout) :: molecule

    ! local variables
    integer                  :: i


    do i = 1, molecule%num_atoms
      molecule%atom_coord(1,i) = molecule%atom_coord(1,i) + move(1)
      molecule%atom_coord(2,i) = molecule%atom_coord(2,i) + move(2)
      molecule%atom_coord(3,i) = molecule%atom_coord(3,i) + move(3)
    end do

    return

  end subroutine shift_molecule

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    wrap_molecule
  !> @brief        wrap molecule coordinates
  !! @authors      NT
  !! @param[in]    option     : option information
  !! @param[in]    trajectory : trajectory information
  !! @param[in]    ensemble   : ensemble information
  !! @param[inout] boundary   : boundary information
  !! @param[inout] molecule   : molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine wrap_molecule(option, trajectory, ensemble, boundary, molecule)

    ! formal arguments
    type(s_option),           intent(in)    :: option
    type(s_trajectory),       intent(in)    :: trajectory
    type(s_ensemble),         intent(in)    :: ensemble
    type(s_boundary),         intent(inout) :: boundary
    type(s_molecule),         intent(inout) :: molecule


    if (boundary%type == BoundaryTypePBC) then
      if(option%determine_box == DetermineBoxTrajectory) then
        if (ensemble%ensemble == EnsembleNPT  .or. &
            ensemble%ensemble == EnsembleNPAT .or. &
            ensemble%ensemble == EnsembleNPgT ) then
          call refresh_boundary(trajectory, boundary)
        end if
      end if

      if (option%wrap) then
        call wrap_atom(trajectory, boundary, molecule)
      else
        call shift_firstq(trajectory, boundary, molecule)
      end if
    end if

    return

  end subroutine wrap_molecule

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_parameter
  !> @brief        input data
  !! @authors      IY
  !! @param[in]    spana_inp_info : spana input information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_parameter(spana_inp_info, top, par, psf, prmtop, grotop, &
                             pdb, crd, ambcrd, grocrd, rst, ref, ambref,    &
                             groref, localres, mode)

    ! formal arguments
    type(s_spana_inp_info),  intent(in)    :: spana_inp_info
    type(s_top),             intent(inout) :: top
    type(s_par),             intent(inout) :: par
    type(s_psf),             intent(inout) :: psf
    type(s_prmtop),          intent(inout) :: prmtop
    type(s_grotop),          intent(inout) :: grotop
    type(s_pdb),             intent(inout) :: pdb
    type(s_crd),             intent(inout) :: crd
    type(s_ambcrd),          intent(inout) :: ambcrd
    type(s_grocrd),          intent(inout) :: grocrd
    type(s_rst),             intent(inout) :: rst
    type(s_pdb),             intent(inout) :: ref
    type(s_ambcrd),          intent(inout) :: ambref
    type(s_grocrd),          intent(inout) :: groref
    type(s_localres),        intent(inout) :: localres
    type(s_mode),            intent(inout) :: mode


    call init_top(top)
    call init_par(par)
    call init_psf(psf)
    call init_prmtop(prmtop)
    call init_grotop(grotop)
    call init_pdb(pdb)
    call init_crd(crd)
    call init_ambcrd(ambcrd)
    call init_grocrd(grocrd)
    call init_pdb(ref)
    call init_ambcrd(ambref)
    call init_grocrd(groref)

    if (spana_inp_info%topfile .ne. '') then
      call input_top(spana_inp_info%topfile, top)
    end if

    if (spana_inp_info%parfile .ne. '') then
      call input_par(spana_inp_info%parfile, par)
    end if

    if (spana_inp_info%strfile .ne. '') then
      call input_str(spana_inp_info%strfile, top, par)
    end if

    if (spana_inp_info%psffile .ne. '') then
      call input_psf(spana_inp_info%psffile, psf)
    end if

    if (spana_inp_info%prmtopfile .ne. '') then
      call input_prmtop(spana_inp_info%prmtopfile, prmtop)
    end if

    if (spana_inp_info%grotopfile .ne. '') then
      call input_grotop(spana_inp_info%grotopfile, grotop)
    end if

    if (spana_inp_info%pdbfile .ne. '') then
      call input_pdb(spana_inp_info%pdbfile, pdb)
    end if

    if (spana_inp_info%crdfile .ne. '') then
      call input_crd(spana_inp_info%crdfile, crd)
    end if

    if (spana_inp_info%ambcrdfile .ne. '') then
      call input_ambcrd(spana_inp_info%ambcrdfile, ambcrd)
    end if

    if (spana_inp_info%grocrdfile .ne. '') then
      call input_grocrd(spana_inp_info%grocrdfile, grocrd)
    end if

    if (spana_inp_info%rstfile .ne. '') then
      call input_rst(spana_inp_info%rstfile, rst)
    end if

    if (spana_inp_info%reffile .ne. '') then
      call input_pdb(spana_inp_info%reffile, ref)
    end if

    if (spana_inp_info%ambreffile .ne. '') then
      call input_ambcrd(spana_inp_info%ambreffile, ambref)
    end if

    if (spana_inp_info%groreffile .ne. '') then
      call input_grocrd(spana_inp_info%groreffile, groref)
    end if

    if (spana_inp_info%local_resfile .ne. '') then
      call input_localres(spana_inp_info%local_resfile, localres)
    end if

    if (spana_inp_info%modefile .ne. '') then
      call input_mode(spana_inp_info%modefile, mode)
    end if

    return

  end subroutine input_parameter

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    wrap_atom
  !> @brief        wrap atom by boundary
  !! @authors      IY
  !! @param[in]    trajectory : trajectory information
  !! @param[in]    boundary   : boundary information
  !! @param[inout] molecule   : molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine wrap_atom(trajectory, boundary, molecule)

    ! formal arguments
    type(s_trajectory),      intent(in)    :: trajectory
    type(s_boundary),        intent(in)    :: boundary
    type(s_molecule),        intent(inout) :: molecule

    ! local variables
    real(wp)                 :: x_shift, y_shift, z_shift
    real(wp)                 :: move(3)
    real(wp)                 :: bsize_x, bsize_y, bsize_z
    integer                  :: i


    bsize_x   = trajectory%pbc_box(1,1)
    bsize_y   = trajectory%pbc_box(2,2)
    bsize_z   = trajectory%pbc_box(3,3)

    do i = 1, molecule%num_atoms

      x_shift = molecule%atom_coord(1,i) - boundary%origin_x
      y_shift = molecule%atom_coord(2,i) - boundary%origin_y
      z_shift = molecule%atom_coord(3,i) - boundary%origin_z

      !coordinate shifted to the first quadrant and set into the boundary box
      !
      move(1) = bsize_x*0.5_wp - bsize_x*anint(x_shift/bsize_x)
      move(2) = bsize_y*0.5_wp - bsize_y*anint(y_shift/bsize_y)
      move(3) = bsize_z*0.5_wp - bsize_z*anint(z_shift/bsize_z)

      x_shift = x_shift + move(1)
      y_shift = y_shift + move(2)
      z_shift = z_shift + move(3)

      molecule%atom_coord(1,i) = x_shift
      molecule%atom_coord(2,i) = y_shift
      molecule%atom_coord(3,i) = z_shift

    end do

    return

  end subroutine wrap_atom

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    shift_firstq
  !> @brief        shift atom coordinate
  !! @authors      IY
  !! @param[in]    trajectory : trajectory information
  !! @param[in]    boundary   : boundary information
  !! @param[inout] molecule   : molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine shift_firstq(trajectory, boundary, molecule)

    ! formal arguments
    type(s_trajectory),      intent(in)    :: trajectory
    type(s_boundary),        intent(in)    :: boundary
    type(s_molecule),        intent(inout) :: molecule

    real(wp)                 :: x_shift, y_shift, z_shift
    real(wp)                 :: bsize_x, bsize_y, bsize_z
    integer                  :: i


    bsize_x   = trajectory%pbc_box(1,1)
    bsize_y   = trajectory%pbc_box(2,2)
    bsize_z   = trajectory%pbc_box(3,3)

    do i = 1, molecule%num_atoms

      x_shift = molecule%atom_coord(1,i)
      y_shift = molecule%atom_coord(2,i)
      z_shift = molecule%atom_coord(3,i)

      x_shift = x_shift + bsize_x*0.5_wp - boundary%origin_x
      y_shift = y_shift + bsize_y*0.5_wp - boundary%origin_y
      z_shift = z_shift + bsize_z*0.5_wp - boundary%origin_z

      molecule%atom_coord(1,i) = x_shift
      molecule%atom_coord(2,i) = y_shift
      molecule%atom_coord(3,i) = z_shift

    end do

    return

  end subroutine shift_firstq

end module sa_tool_mod
