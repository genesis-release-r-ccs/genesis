!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_enefunc_spica_mod
!> @brief   define potential energy functions
!! @authors Jaewoon Jung (JJ), Takeshi Imai (TI), Chigusa Kobayashi (CK),
!!          Takaharu Mori (TM), Yuji Sugita (YS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_enefunc_spica_mod

  use sp_enefunc_localres_mod
  use sp_enefunc_restraints_mod
  use sp_enefunc_fit_mod
  use sp_enefunc_table_mod
  use sp_energy_mod
  use sp_restraints_str_mod
  use sp_constraints_str_mod
  use sp_enefunc_str_mod
  use sp_energy_str_mod
  use sp_domain_str_mod
  use sp_communicate_mod
  use dihedral_libs_mod
  use nbond_list_mod
  use molecules_str_mod
  use fileio_par_mod
  use fileio_localres_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif
#ifdef OMP
  use omp_lib, only: omp_get_thread_num
#endif

  implicit none
  private

  ! subroutines
  public  :: define_enefunc_spica
  private :: setup_enefunc_bond
  private :: setup_enefunc_enm 
  private :: setup_enefunc_angl
  private :: setup_enefunc_dihe
  private :: setup_enefunc_impr ! yyzhang def spica improper
  private :: setup_enefunc_nonb
  public  :: count_nonb_excl_spica

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc_spica
  !> @brief        a driver subroutine for defining potential energy functions
  !! @authors      YS, TI, JJ, CK
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    par         : SPICA PAR information
  !! @param[in]    localres    : local restraint information
  !! @param[in]    molecule    : molecule information
  !! @param[inout] constraints : constraints information
  !! @param[in]    restraints  : restraints information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc_spica(ene_info, par, localres, molecule,        &
                                  constraints, restraints, domain, enefunc, &
                                  comm)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_par),             intent(in)    :: par
    type(s_localres),        intent(in)    :: localres
    type(s_molecule),        intent(in)    :: molecule
    type(s_constraints),     intent(inout) :: constraints
    type(s_restraints),      intent(in)    :: restraints
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_comm),            intent(inout) :: comm   

    ! local variables
    integer                  :: ncel, ncelb


    ! base
    !
    ncel  = domain%num_cell_local
    ncelb = domain%num_cell_local + domain%num_cell_boundary

    call alloc_enefunc(enefunc, EneFuncBase, ncel, ncelb)
    if (molecule%num_enm_bonds > 0) then
      enefunc%enm_use = .true.
    end if
    call alloc_enefunc(enefunc, EneFuncExcl, ncel, ncel)

    ! bond
    !
    call setup_enefunc_bond(par, molecule, domain, constraints, enefunc)

    ! enm
    !
    if (enefunc%enm_use) &
      call setup_enefunc_enm(par, molecule, domain, enefunc, comm)

    ! angle
    !
    call setup_enefunc_angl(par, molecule, domain, constraints, enefunc)

    ! dihedral
    !
    call setup_enefunc_dihe(par, molecule, domain, enefunc)

    ! improper yyzhang def spica improper coplanar
    !
    call setup_enefunc_impr(par, molecule, domain, enefunc)

    ! nonbonded
    !
    call setup_enefunc_nonb(par, molecule, constraints, domain, enefunc)

    ! lookup table
    !
    if (ene_info%table) &
    call setup_enefunc_table(ene_info, enefunc)

    ! restraint
    !
    call setup_enefunc_restraints(molecule, restraints, domain, enefunc)

    call setup_enefunc_localres(localres, domain, enefunc)

    call alloc_enefunc(enefunc, EneFuncBondCell, ncel, ncelb)

    ! write summary of energy function
    !
    if (main_rank) then
      write(MsgOut,'(A)') &
           'Define_Enefunc_Spica> Number of Interactions in Each Term'
      write(MsgOut,'(A20,I10,A20,I10)')                  &
           '  bond_ene        = ', enefunc%num_bond_all, &
           '  angle_ene       = ', enefunc%num_angl_all
      write(MsgOut,'(A20,I10,A20,I10)')                  &
           '  torsion_ene     = ', enefunc%num_dihe_all, &
           '  improper_ene    = ', enefunc%num_impr_all
      write(MsgOut,'(A20,I10)')                          &
           '  cmap_ene        = ', enefunc%num_cmap_all
      write(MsgOut,'(A20,I10,A20,I10)')                  &
           '  nb_exclusions   = ', enefunc%num_excl_all, &
           '  nb14_calc       = ', enefunc%num_nb14_all
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine define_enefunc_spica

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_bond
  !> @brief        define BOND term for each cell in potential energy function
  !! @authors      YS, JJ, TM
  !! @param[in]    par         : CHARMM PAR information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    domain      : domain information
  !! @param[in]    constraints : constraints information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond(par, molecule, domain, constraints, enefunc)

    ! formal arguments
    type(s_par),                 intent(in)    :: par
    type(s_molecule),    target, intent(in)    :: molecule
    type(s_domain),      target, intent(in)    :: domain
    type(s_constraints), target, intent(inout) :: constraints
    type(s_enefunc),     target, intent(inout) :: enefunc

    ! local variable
    integer                  :: dupl, ioffset
    integer                  :: i, j, icel_local, k
    integer                  :: icel1, icel2, ix1, ix2
    integer                  :: nbond, nbond_p, i1, i2, i3, found, found1
    integer                  :: nenmt_p, nenmp_p, key, ndat, num
    character(6)             :: ci1, ci2, ci3, ri1, ri2
    real(wp)                 :: dij(3), rij2, pairdist2
    integer,     allocatable :: temporary_list(:,:)

    real(wp),        pointer :: force(:,:), dist(:,:)
    integer,         pointer :: nwater(:), bond(:), list(:,:,:)
    integer,         pointer :: nexcl(:), n14(:)
    integer,         pointer :: excl_list(:,:,:), n14_list(:,:,:)
    integer,         pointer :: ncel
    integer(int2),   pointer :: cell_pair(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: mol_bond_list(:,:)
    character(6),    pointer :: mol_cls_name(:)
    real(wp),        pointer :: coord_pbc(:,:,:)
    integer,         pointer :: bond_pbc(:,:)


    mol_bond_list => molecule%bond_list
    mol_cls_name  => molecule%atom_cls_name

    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l
    coord_pbc => domain%coord_pbc 
    nwater    => domain%num_water

    bond      => enefunc%num_bond

    nbond     = molecule%num_bonds
    nbond_p   = par%num_bonds
    nenmt_p   = par%num_enmt
    nenmp_p   = par%num_enmp

    pairdist2 = enefunc%pairlistdist * enefunc%pairlistdist
    enefunc%table%TIP2_bond_calc = .false.

    do dupl = 1, domain%num_duplicate

      ioffset = (dupl-1) * enefunc%table%num_all

      do i = 1, nbond

        i1  = molecule%bond_list(1,i)
        i2  = molecule%bond_list(2,i)
        ci1 = mol_cls_name(i1)
        ci2 = mol_cls_name(i2)
        ri1 = molecule%residue_name(i1)
        ri2 = molecule%residue_name(i2)

        if (ri1(1:3) .ne. 'TIP' .and. ri1(1:3) .ne. 'WAT' .and. &
            ri1(1:3) .ne. 'SOL' .and. ri1(1:3) .ne. 'MIZ' .and. &
            ri1(1:3) .ne. 'PWA' .and. ri2(1:3) .ne. 'PWA' .and. &
            ri2(1:3) .ne. 'TIP' .and. ri2(1:3) .ne. 'WAT' .and. &
            ri2(1:3) .ne. 'SOL' .and. ri2(1:3) .ne. 'MIZ') then

          i1  = i1 + ioffset
          i2  = i2 + ioffset
          icel1 = id_g2l(1,i1)
          icel2 = id_g2l(1,i2)

          ! Check if it is in my domain
          !
          if (icel1 /= 0 .and. icel2 /= 0) then

            icel_local = cell_pair(icel1,icel2)

            if (icel_local > 0 .and. icel_local <= ncel) then

              do j = 1, nbond_p
                if ((ci1 .eq. par%bond_atom_cls(1,j) .and.  &
                     ci2 .eq. par%bond_atom_cls(2,j)) .or.  &
                    (ci1 .eq. par%bond_atom_cls(2,j) .and.  &
                     ci2 .eq. par%bond_atom_cls(1,j))) then
 
                  bond (icel_local) = bond(icel_local) + 1
                          
                  exit
                end if
              end do

            end if
          end if

        else if (.not.constraints%rigid_bond) then

          do j = 1, nbond_p
            if ((ci1 .eq. par%bond_atom_cls(1,j) .and.  &
                 ci2 .eq. par%bond_atom_cls(2,j)) .or.  &
                (ci1 .eq. par%bond_atom_cls(2,j) .and.  &
                 ci2 .eq. par%bond_atom_cls(1,j))) then

              enefunc%table%TIP2_bond_calc = .true.
              enefunc%table%TIP2_bond = par%bond_dist_min(j)
              enefunc%table%TIP2_force = par%bond_force_const(j)

              exit
            end if
          end do

        end if

      end do
    end do

    k = 0
    do i = 1, ncel
      k = max(k,bond(i))
    end do
#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, k, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
#endif
    MaxBond  = 2*k
    BondMove = MaxBond / 3 ! yyzhang for small testing system, just set BondMove = MaxBond
    bond(1:ncel) = 0
    if (main_rank) then
      write(MsgOut,*) 'sp_enefunc_str    ::MaxBond      : ', MaxBond
      write(MsgOut,*) 'sp_enefunc_str    ::BondMove     : ', BondMove
      write(MsgOut,*) ' '
    end if
    call alloc_enefunc(enefunc, EneFuncBond, ncel, ncel)
    list      => enefunc%bond_list
    force     => enefunc%bond_force_const
    dist      => enefunc%bond_dist_min
    bond_pbc  => enefunc%bond_pbc

    do dupl = 1, domain%num_duplicate

      ioffset = (dupl-1) * enefunc%table%num_all

      do i = 1, nbond

        i1  = molecule%bond_list(1,i)
        i2  = molecule%bond_list(2,i)
        ci1 = mol_cls_name(i1)
        ci2 = mol_cls_name(i2)
        ri1 = molecule%residue_name(i1)
        ri2 = molecule%residue_name(i2)

        if (ri1(1:3) .ne. 'TIP' .and. ri1(1:3) .ne. 'WAT' .and. &
            ri1(1:3) .ne. 'SOL' .and. ri1(1:3) .ne. 'MIZ' .and. &
            ri1(1:3) .ne. 'PWA' .and. ri2(1:3) .ne. 'PWA' .and. &
            ri2(1:3) .ne. 'TIP' .and. ri2(1:3) .ne. 'WAT' .and. &
            ri2(1:3) .ne. 'SOL' .and. ri2(1:3) .ne. 'MIZ') then

          i1  = i1 + ioffset
          i2  = i2 + ioffset
          icel1 = id_g2l(1,i1)
          icel2 = id_g2l(1,i2)

          ! Check if it is in my domain
          !
          if (icel1 /= 0 .and. icel2 /= 0) then

            icel_local = cell_pair(icel1,icel2)

            if (icel_local > 0 .and. icel_local <= ncel) then

              do j = 1, nbond_p
                if ((ci1 .eq. par%bond_atom_cls(1,j) .and.  &
                     ci2 .eq. par%bond_atom_cls(2,j)) .or.  &
                    (ci1 .eq. par%bond_atom_cls(2,j) .and.  &
                     ci2 .eq. par%bond_atom_cls(1,j))) then

                  bond (icel_local) = bond(icel_local) + 1
                  list (1,bond(icel_local),icel_local) = i1
                  list (2,bond(icel_local),icel_local) = i2
                  force(bond(icel_local),icel_local)   = &
                        par%bond_force_const(j)
                  dist (bond(icel_local),icel_local)   = &
                        par%bond_dist_min(j)
                  bond_pbc(bond(icel_local),icel_local) = 13
                          
                  exit
                end if
              end do

              if (j == nbond_p + 1) &
                write(MsgOut,*) &
                  'Setup_Enefunc_Bond> not found BOND: [', &
                  ci1, ']-[', ci2, '] in parameter file. (ERROR)'

            end if
 
          end if

        else if (.not.constraints%rigid_bond) then

          do j = 1, nbond_p
            if ((ci1 .eq. par%bond_atom_cls(1,j) .and.  &
                 ci2 .eq. par%bond_atom_cls(2,j)) .or.  &
                (ci1 .eq. par%bond_atom_cls(2,j) .and.  &
                 ci2 .eq. par%bond_atom_cls(1,j))) then

              enefunc%table%TIP2_bond_calc = .true.
              enefunc%table%TIP2_bond = par%bond_dist_min(j)
              enefunc%table%TIP2_force = par%bond_force_const(j)

              exit
            end if
          end do

        end if

      end do
    end do

    found  = 0
    found1 = 0
    do i = 1, ncel
      found = found + bond(i)
      if (constraints%water_type == TIP2) found = found + nwater(i)
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_bond_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_bond_all = found
#endif

!   if (enefunc%num_bond_all /= nbond*domain%num_duplicate) &
!     call error_msg('Setup_Enefunc_Bond> Some bond paremeters are missing.')

    return

  end subroutine setup_enefunc_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_enm 
  !> @brief        define BOND term for each cell in potential energy function
  !! @authors      YS, JJ, TM
  !! @param[in]    par      : CHARMM PAR information
  !! @param[in]    molecule : molecule information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_enm(par, molecule, domain, enefunc, comm)

    ! formal arguments
    type(s_par),                 intent(in)    :: par
    type(s_molecule),    target, intent(in)    :: molecule
    type(s_domain),      target, intent(in)    :: domain
    type(s_enefunc),     target, intent(inout) :: enefunc
    type(s_comm),        target, intent(inout) :: comm

    ! local variable
    integer                  :: dupl, ioffset
    integer                  :: i, j, icel_local, k, ncelb
    integer                  :: icel1, icel2, ix1, ix2
    integer                  :: nbond, nbond_p, i1, i2, i3, found, found1
    integer                  :: nenm, nenmt_p, nenmp_p, key, ndat, num
    character(6)             :: ci1, ci2, ci3, ri1, ri2
    real(wp)                 :: dij(3), rij2, pairdist2
    integer,     allocatable :: temporary_list(:,:)

    real(wp),        pointer :: force(:,:), dist(:,:)
    integer,         pointer :: enm(:), list(:,:,:)
    integer,         pointer :: nexcl(:), n14(:)
    integer,         pointer :: excl_list(:,:,:), n14_list(:,:,:)
    integer,         pointer :: ncel
    integer(int2),   pointer :: cell_pair(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: mol_bond_list(:,:)
    character(6),    pointer :: mol_cls_name(:)
    real(wp),        pointer :: coord_pbc(:,:,:)

    mol_bond_list => molecule%enm_list
    mol_cls_name  => molecule%atom_cls_name

    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l
    coord_pbc => domain%coord_pbc 

    enm       => enefunc%num_enm

    nenm      = molecule%num_enm_bonds
    nenmt_p   = par%num_enmt
    nenmp_p   = par%num_enmp

    ncelb     = ncel + domain%num_cell_boundary

    pairdist2 = enefunc%pairlistdist * enefunc%pairlistdist

    do dupl = 1, domain%num_duplicate

      ioffset = (dupl-1) * enefunc%table%num_all

      do i = 1, nenm

        i1  = molecule%enm_list(1,i)
        i2  = molecule%enm_list(2,i)
        ci1 = mol_cls_name(i1)
        ci2 = mol_cls_name(i2)
        ri1 = molecule%residue_name(i1)
        ri2 = molecule%residue_name(i2)

        i1  = i1 + ioffset
        i2  = i2 + ioffset

        icel1 = id_g2l(1,i1)
        icel2 = id_g2l(1,i2)

        ! Check if it is in my domain
        !
        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1,icel2)

          if (icel_local > 0 .and. icel_local <= ncel) then

            key = 0
            do j = 1, nenmt_p
              if ((par%bond_enm_atom_id(1,j) == i1 .and. &
                   par%bond_enm_atom_id(2,j) == i2) .or. &
                  (par%bond_enm_atom_id(1,j) == i2 .and. &
                   par%bond_enm_atom_id(2,j) == i1)) then
                do k = 1, nenmp_p
                  if (par%bond_enm_param_type(k) == &
                      par%bond_enm_ij_type(j)) then
                    enm(icel_local) = enm(icel_local) + 1
                    key = 1
                    exit
                  end if
                end do
              end if
              if (key == 1) exit
            end do

            if (j == nenmt_p + 1) &
                write(MsgOut,*) &
                  'Setup_Enefunc_Bond> not found BOND: [', &
                  ci1, ']-[', ci2, '] in parameter file. (ERROR)'
 
          end if
        end if
      end do
    end do

    found = 0
    do i = 1, ncel
      found = max(found,enm(i))
    end do
#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, found, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
#endif
    MaxEnm  = 3*found
    enm(1:ncelb) = 0
    if (main_rank) then
      write(MsgOut,*) 'sp_enefunc_str    ::MaxEnm       : ', MaxEnm
      write(MsgOut,*) ' '
    end if
    call alloc_enefunc(enefunc, EneFuncEnm, ncelb, ncelb)
    list      => enefunc%enm_list
    force     => enefunc%enm_force_const
    dist      => enefunc%enm_dist_min

    do dupl = 1, domain%num_duplicate

      ioffset = (dupl-1) * enefunc%table%num_all

      do i = 1, nenm

        i1  = molecule%enm_list(1,i)
        i2  = molecule%enm_list(2,i)
        ci1 = mol_cls_name(i1)
        ci2 = mol_cls_name(i2)
        ri1 = molecule%residue_name(i1)
        ri2 = molecule%residue_name(i2)

        i1  = i1 + ioffset
        i2  = i2 + ioffset

        icel1 = id_g2l(1,i1)
        icel2 = id_g2l(1,i2)

        ! Check if it is in my domain
        !
        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1,icel2)

!         if (icel_local > 0 .and. icel_local <= ncel) then
          if (icel_local > 0 .and. icel_local <= ncelb) then

            key = 0
            do j = 1, nenmt_p
              if ((par%bond_enm_atom_id(1,j) == i1 .and. &
                   par%bond_enm_atom_id(2,j) == i2) .or. &
                  (par%bond_enm_atom_id(1,j) == i2 .and. &
                   par%bond_enm_atom_id(2,j) == i1)) then
                do k = 1, nenmp_p
                  if (par%bond_enm_param_type(k) == &
                      par%bond_enm_ij_type(j)) then
                    enm(icel_local) = enm(icel_local) + 1
                    list(1,enm(icel_local),icel_local) = i1
                    list(2,enm(icel_local),icel_local) = i2
                    force(enm(icel_local),icel_local) = &
                          par%bond_enm_force_const(k)
                    dist (enm(icel_local),icel_local) = &
                          par%bond_enm_dist_min(k)
                    key = 1
                    exit
                  end if
                end do
              end if
              if (key == 1) exit
            end do

            if (j == nenmt_p + 1) &
                write(MsgOut,*) &
                  'Setup_Enefunc_Bond> not found ENM BOND: [', &
                  ci1, ']-[', ci2, '] in parameter file. (ERROR)'
 
          end if
        end if
      end do
    end do

    found  = 0
    found1 = 0
    do i = 1, ncel
      found = found + enm(i)
    end do
#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_enm_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_enm_all = found
#endif

    if (enefunc%num_enm_all /= nenm*domain%num_duplicate) &
      call error_msg('Setup_Enefunc_Bond> Some bond paremeters are '//&
                     'missing.')

!   call update_cell_size_enm(domain, enefunc, comm)
!   call update_cell_boundary_enm(domain, enefunc, comm)

    return

  end subroutine setup_enefunc_enm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_angl
  !> @brief        define ANGLE term for each cell in potential energy function
  !! @authors      YS, JJ, TM
  !! @param[in]    par         : CHARMM PAR information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    domain      : domain information
  !! @param[in]    constraints : constraints information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_angl(par, molecule, domain, constraints, enefunc)

    ! formal arguments
    type(s_par),         target, intent(in)    :: par
    type(s_molecule),    target, intent(in)    :: molecule
    type(s_domain),      target, intent(in)    :: domain
    type(s_constraints), target, intent(inout) :: constraints
    type(s_enefunc),     target, intent(inout) :: enefunc

    ! local variables
    integer                  :: dupl, ioffset
    integer                  :: i, j, icel_local, k
    integer                  :: icel1, icel2
    integer                  :: nangl, nangl_p, found
    integer                  :: nenmu_p, nenmq_p, key
    integer                  :: list(3)
    character(6)             :: ci1, ci2, ci3, ri1, ri2, ri3

    real(wp),        pointer :: force(:,:), theta(:,:)
    real(wp),        pointer :: ubforce(:,:), ubrmin(:,:)
    integer,         pointer :: angle(:), alist(:,:,:)
    integer,         pointer :: ncel
    integer(int2),   pointer :: cell_pair(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: mol_angl_list(:,:)
    integer,         pointer :: nwater(:)
    character(6),    pointer :: mol_cls_name(:), mol_res_name(:)
    integer,         pointer :: angl_pbc(:,:,:)


    mol_angl_list => molecule%angl_list
    mol_cls_name  => molecule%atom_cls_name
    mol_res_name  => molecule%residue_name

    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l
    nwater    => domain%num_water

    angle     => enefunc%num_angle

    nangl     = molecule%num_angles
    nangl_p   = par%num_angles
    nenmu_p = par%num_enmu
    nenmq_p = par%num_enmq

    do dupl = 1, domain%num_duplicate

      ioffset = (dupl-1)*enefunc%table%num_all

      do i = 1, nangl

        list(1:3) = mol_angl_list(1:3,i)
        ci1 = mol_cls_name(list(1))
        ci2 = mol_cls_name(list(2))
        ci3 = mol_cls_name(list(3))
        ri1 = mol_res_name(list(1))
        ri2 = mol_res_name(list(2))
        ri3 = mol_res_name(list(3))

        list(1:3) = list(1:3) + ioffset
        icel1 = id_g2l(1,list(1))
        icel2 = id_g2l(1,list(3))

        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1,icel2)

          if (icel_local >= 1 .and. icel_local <= ncel) then

            if (nenmu_p > 0) then
              key = 0
              do j = 1, nenmu_p
                if ((par%angle_enm_atom_id(1,j) == list(1) .and. &
                     par%angle_enm_atom_id(2,j) == list(2) .and. &
                     par%angle_enm_atom_id(3,j) == list(3)) .or. &
                    (par%angle_enm_atom_id(3,j) == list(1) .and. &
                     par%angle_enm_atom_id(2,j) == list(2) .and. &
                     par%angle_enm_atom_id(1,j) == list(3))) then
                  do k = 1, nenmq_p
                    if (par%angle_enm_param_type(k) == &
                        par%angle_enm_ij_type(j)) then
                      angle(icel_local) = angle(icel_local) + 1
                      key = 1
                      exit
                    end if
                  end do
                end if
                if (key == 1) exit
              end do
            end if

            if (nenmu_p == 0 .or. nenmq_p == 0 .or. &
                (nenmu_p > 0 .and. key == 0)) then
              do j = 1, nangl_p
                if ((ci1 .eq. par%angl_atom_cls(1,j) .and. &
                     ci2 .eq. par%angl_atom_cls(2,j) .and. &
                     ci3 .eq. par%angl_atom_cls(3,j)) .or. &
                    (ci1 .eq. par%angl_atom_cls(3,j) .and. &
                     ci2 .eq. par%angl_atom_cls(2,j) .and. &
                     ci3 .eq. par%angl_atom_cls(1,j))) then
                  angle(icel_local) = angle(icel_local) + 1
                  exit
                end if
              end do
  
              if (j == nangl_p + 1) &
                write(MsgOut,*) &
                  'Setup_Enefunc_Angl> not found ANGL: [',ci1, &
                  ']-[', ci2, ']-[', ci3, '] in parameter file. (ERROR)'
    
            end if
          end if 
        end if
      end do
    end do

    found = 0
    do i = 1, ncel
      found = max(found,angle(i))
    end do
#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, found, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
#endif
    MaxAngle  = 2*found
    AngleMove = MaxAngle / 3
    angle(1:ncel) = 0
    if (main_rank) then
      write(MsgOut,*) 'sp_enefunc_str    ::MaxAngle     : ', MaxAngle
      write(MsgOut,*) 'sp_enefunc_str    ::AngleMove    : ', AngleMove
      write(MsgOut,*) ' '
    end if
    call alloc_enefunc(enefunc, EneFuncAngl, ncel, ncel)
    alist     => enefunc%angle_list
    force     => enefunc%angle_force_const
    theta     => enefunc%angle_theta_min
    ubforce   => enefunc%urey_force_const
    ubrmin    => enefunc%urey_rmin
    angl_pbc  => enefunc%angle_pbc
      
    do dupl = 1, domain%num_duplicate

      ioffset = (dupl-1)*enefunc%table%num_all

      do i = 1, nangl

        list(1:3) = mol_angl_list(1:3,i)
        ci1 = mol_cls_name(list(1))
        ci2 = mol_cls_name(list(2))
        ci3 = mol_cls_name(list(3))
        ri1 = mol_res_name(list(1))
        ri2 = mol_res_name(list(2))
        ri3 = mol_res_name(list(3))

        list(1:3) = list(1:3) + ioffset

        icel1 = id_g2l(1,list(1))
        icel2 = id_g2l(1,list(3))

        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1,icel2)

          if (icel_local >= 1 .and. icel_local <= ncel) then

            if (nenmu_p > 0) then
              key = 0
              do j = 1, nenmu_p
                if ((par%angle_enm_atom_id(1,j) == list(1) .and. &
                     par%angle_enm_atom_id(2,j) == list(2) .and. &
                     par%angle_enm_atom_id(3,j) == list(3)) .or. &
                    (par%angle_enm_atom_id(3,j) == list(1) .and. &
                     par%angle_enm_atom_id(2,j) == list(2) .and. &
                     par%angle_enm_atom_id(1,j) == list(3))) then
                  do k = 1, nenmq_p
                    if (par%angle_enm_param_type(k) == &
                        par%angle_enm_ij_type(j)) then
                      angle(icel_local) = angle(icel_local) + 1
                      alist(1:3,angle(icel_local),icel_local) = list(1:3)
                      force(angle(icel_local),icel_local) = &
                            par%angle_enm_force_const(k)
                      theta(angle(icel_local),icel_local) = &
                            par%angle_enm_dist_min(k) * RAD
                      ubforce(angle(icel_local),icel_local) = &
                            par%angle_enm_lj_eps(k)
                      ubrmin(angle(icel_local),icel_local) = &
                            par%angle_enm_lj_sigma(k)
                      angl_pbc(1:3,angle(icel_local),icel_local) = 13
                      key = 1
                      exit
                    end if
                  end do
                end if
                if (key == 1) exit
              end do
            end if

            if (nenmu_p == 0 .or. nenmq_p == 0 .or. &
                (nenmu_p > 0 .and. key == 0)) then
              do j = 1, nangl_p
                if ((ci1 .eq. par%angl_atom_cls(1,j) .and. &
                     ci2 .eq. par%angl_atom_cls(2,j) .and. &
                     ci3 .eq. par%angl_atom_cls(3,j)) .or. &
                    (ci1 .eq. par%angl_atom_cls(3,j) .and. &
                     ci2 .eq. par%angl_atom_cls(2,j) .and. &
                     ci3 .eq. par%angl_atom_cls(1,j))) then
    
                  angle(icel_local) = angle(icel_local) + 1
                  alist(1:3,angle(icel_local),icel_local) = list(1:3)
      
                  force(angle(icel_local),icel_local)   = &
                       par%angl_force_const(j)
                  theta(angle(icel_local),icel_local)   = &
                       par%angl_theta_min(j)*RAD
                  ubforce(angle(icel_local),icel_local) = &
                       par%angl_ub_force_const(j)
                  ubrmin(angle(icel_local),icel_local)  = &
                     par%angl_ub_rmin(j)
                  angl_pbc(1:3,angle(icel_local),icel_local) = 13
                  exit
    
                end if
              end do
  
              if (j == nangl_p + 1) &
                write(MsgOut,*) &
                  'Setup_Enefunc_Angl> not found ANGL: [',ci1, &
                  ']-[', ci2, ']-[', ci3, '] in parameter file. (ERROR)'
    
            end if
          end if 
        end if
      end do
    end do

    found = 0
    do i = 1, ncel
      found = found + angle(i)
      if (enefunc%table%water_bond_calc_OH) found = found + nwater(i)
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_angl_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_angl_all = found
#endif

    if (enefunc%num_angl_all /= nangl*domain%num_duplicate) &
      call error_msg('Setup_Enefunc_Angl> Some angle paremeters are '//&
                     'missing.')

    return

  end subroutine setup_enefunc_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_dihe
  !> @brief        define DIHEDRAL term in potential energy function
  !! @authors      YS, JJ, TM
  !! @param[in]    par      : CHARMM PAR information
  !! @param[in]    molecule : molecule information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_dihe(par, molecule, domain, enefunc)

    ! formal arguments
    type(s_par),      target, intent(in)    :: par
    type(s_molecule), target, intent(in)    :: molecule
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variables
    integer                   :: dupl, ioffset
    integer                   :: ndihe, ndihe_p
    integer                   :: i, j, icel_local
    integer                   :: icel1, icel2
    integer                   :: found, nw_found
    integer                   :: list(4)
    character(6)              :: ci1, ci2, ci3, ci4

    real(wp),         pointer :: force(:,:), phase(:,:)
    integer,          pointer :: dihedral(:), dlist(:,:,:), period(:,:)
    integer,          pointer :: ncel
    integer(int2),    pointer :: cell_pair(:,:)
    integer(int2),    pointer :: id_g2l(:,:)
    logical,      allocatable :: no_wild(:)
    integer,          pointer :: notation
    integer,          pointer :: mol_dihe_list(:,:)
    character(6),     pointer :: mol_cls_name(:)
    integer,          pointer :: dihe_pbc(:,:,:)


    mol_dihe_list  => molecule%dihe_list
    mol_cls_name   => molecule%atom_cls_name

    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l

    dihedral  => enefunc%num_dihedral
    notation  => enefunc%notation_14types
    notation = 100

    ndihe     = molecule%num_dihedrals
    ndihe_p   = par%num_dihedrals

    ! check usage of wild card
    !
    allocate(no_wild(ndihe_p))

    do i = 1, ndihe_p

      if ((par%dihe_atom_cls(1,i) .ne. WildCard) .and. &
          (par%dihe_atom_cls(4,i) .ne. WildCard)) then
        ! A-B-C-D type
        no_wild(i) = .true.
      else
        ! X-B-C-D type
        no_wild(i) = .false.
      end if

    end do

    ! find number of interactions
    !
    do dupl = 1, domain%num_duplicate

      ioffset = (dupl-1) * enefunc%table%num_all

      do i = 1, ndihe

        list(1:4) = mol_dihe_list(1:4,i)
        ci1 = mol_cls_name(list(1))
        ci2 = mol_cls_name(list(2))
        ci3 = mol_cls_name(list(3))
        ci4 = mol_cls_name(list(4))

        list(1:4) = list(1:4) + ioffset
        icel1 = id_g2l(1,list(1))
        icel2 = id_g2l(1,list(4))

        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1,icel2)

          if (icel_local >= 1 .and. icel_local <= ncel) then

            nw_found = 0
            do j = 1, ndihe_p
              if (no_wild(j)) then
                if (((ci1 .eq. par%dihe_atom_cls(1,j)) .and. &
                     (ci2 .eq. par%dihe_atom_cls(2,j)) .and. &
                     (ci3 .eq. par%dihe_atom_cls(3,j)) .and. &
                     (ci4 .eq. par%dihe_atom_cls(4,j))) .or. &
                    ((ci1 .eq. par%dihe_atom_cls(4,j)) .and. &
                     (ci2 .eq. par%dihe_atom_cls(3,j)) .and. &
                     (ci3 .eq. par%dihe_atom_cls(2,j)) .and. &
                     (ci4 .eq. par%dihe_atom_cls(1,j)))) then
  
                  nw_found = nw_found + 1
                  dihedral(icel_local) = dihedral(icel_local) + 1
                  exit
  
                end if
              end if
            end do
  
            if (nw_found == 0) then
              do j = 1, ndihe_p
                if (.not.no_wild(j)) then
                  if (((ci2 .eq. par%dihe_atom_cls(2,j)) .and. &
                       (ci3 .eq. par%dihe_atom_cls(3,j))) .or. &
                      ((ci2 .eq. par%dihe_atom_cls(3,j)) .and. &
                       (ci3 .eq. par%dihe_atom_cls(2,j)))) then
  
                    dihedral(icel_local) = dihedral(icel_local) + 1
                    exit  
                  end if
                end if
              end do
            end if
          end if
        end if
      end do
    end do

    found = 0
    do i = 1, ncel
      found = max(found, dihedral(i))
    end do
#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, found, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
#endif
    MaxDihe = MaxDihe * 2
    DiheMove = MaxDihe / 3
    if (main_rank) then
      write(MsgOut,*) 'sp_enefunc_str    ::MaxDihe      : ', MaxDihe
      write(MsgOut,*) 'sp_enefunc_str    ::DiheMove     : ', DiheMove
      write(MsgOut,*) ' '
    end if
    call alloc_enefunc(enefunc, EneFuncDihe, ncel, ncel)
    dlist     => enefunc%dihe_list
    force     => enefunc%dihe_force_const
    period    => enefunc%dihe_periodicity
    phase     => enefunc%dihe_phase
    dihe_pbc  => enefunc%dihe_pbc
    dihedral(1:ncel) = 0

    do dupl = 1, domain%num_duplicate

      ioffset = (dupl-1) * enefunc%table%num_all

      do i = 1, ndihe

        list(1:4) = mol_dihe_list(1:4,i)
        ci1 = mol_cls_name(list(1))
        ci2 = mol_cls_name(list(2))
        ci3 = mol_cls_name(list(3))
        ci4 = mol_cls_name(list(4))

        list(1:4) = list(1:4) + ioffset
        icel1 = id_g2l(1,list(1))
        icel2 = id_g2l(1,list(4))

        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1,icel2)

          if (icel_local >= 1 .and. icel_local <= ncel) then

            nw_found = 0
            do j = 1, ndihe_p
              if (no_wild(j)) then
                if (((ci1 .eq. par%dihe_atom_cls(1,j)) .and. &
                     (ci2 .eq. par%dihe_atom_cls(2,j)) .and. &
                     (ci3 .eq. par%dihe_atom_cls(3,j)) .and. &
                     (ci4 .eq. par%dihe_atom_cls(4,j))) .or. &
                    ((ci1 .eq. par%dihe_atom_cls(4,j)) .and. &
                     (ci2 .eq. par%dihe_atom_cls(3,j)) .and. &
                     (ci3 .eq. par%dihe_atom_cls(2,j)) .and. &
                     (ci4 .eq. par%dihe_atom_cls(1,j)))) then
  
                  nw_found = nw_found + 1
                  dihedral(icel_local) = dihedral(icel_local) + 1
                  dlist(1:4,dihedral(icel_local),icel_local) = list(1:4)
  
                  force (dihedral(icel_local),icel_local) = &
                       par%dihe_force_const(j)
                  period(dihedral(icel_local),icel_local) = &
                       par%dihe_periodicity(j)
                  phase (dihedral(icel_local),icel_local) = &
                       par%dihe_phase(j) * RAD
                  dihe_pbc(1:3,dihedral(icel_local),icel_local) = 13 
                  if (period(dihedral(icel_local),icel_local) >  &
                      enefunc%notation_14types) &
                    call error_msg('Setup_Enefunc_Dihe> '//&
                                   'Too many periodicity.')

                  exit
  
                end if
              end if
            end do
  
            if (nw_found == 0) then
              do j = 1, ndihe_p
                if (.not.no_wild(j)) then
                  if (((ci2 .eq. par%dihe_atom_cls(2,j)) .and. &
                       (ci3 .eq. par%dihe_atom_cls(3,j))) .or. &
                      ((ci2 .eq. par%dihe_atom_cls(3,j)) .and. &
                       (ci3 .eq. par%dihe_atom_cls(2,j)))) then
  
                    dihedral(icel_local) = dihedral(icel_local) + 1
                    dlist(1:4,dihedral(icel_local),icel_local) = list(1:4)
  
                    force (dihedral(icel_local),icel_local) = &
                         par%dihe_force_const(j)
                    period(dihedral(icel_local),icel_local) = &
                         par%dihe_periodicity(j)
                    phase (dihedral(icel_local),icel_local) = &
                         par%dihe_phase(j) * RAD
                    dihe_pbc(1:3,dihedral(icel_local),icel_local) = 13 
                    if (period(dihedral(icel_local),icel_local) >  &
                        enefunc%notation_14types) &
                      call error_msg('Setup_Enefunc_Dihe> '//&
                                     'Too many periodicity.')
                    exit  
                  end if
                end if
              end do
            end if
  
          end if
        end if
      end do
    end do

    deallocate(no_wild)

    found = 0
    do i = 1, ncel
      found = found + dihedral(i)
    end do
#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_dihe_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_dihe_all = found
#endif

    if (enefunc%num_dihe_all < ndihe*domain%num_duplicate) &
      call error_msg( &
         'Setup_Enefunc_Dihe> Some dihedral paremeters are missing.')

    return

  end subroutine setup_enefunc_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_impr
  !> @brief        define IMPROPER term in potential energy function
  !! @authors      YS,JJ,HO
  !! @param[in]    par      : CHARMM PAR information
  !! @param[in]    molecule : molecule information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_impr(par, molecule, domain, enefunc)

    ! formal variables
    type(s_par),      target, intent(in)    :: par
    type(s_molecule), target, intent(in)    :: molecule
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variables
    integer                   :: dupl, ioffset
    integer                   :: nimpr, nimpr_p
    integer                   :: i, j, icel_local
    integer                   :: icel1, icel2
    integer                   :: found, nw_found
    integer                   :: list(4), lists(4)
    character(6)              :: ci1, ci2, ci3, ci4

    real(wp),         pointer :: force(:,:), phase(:,:), coord(:,:)
    real(wp),         pointer :: box_size(:)
    integer,          pointer :: improper(:), ilist(:,:,:)
    integer,          pointer :: ncel
    integer(int2),    pointer :: cell_pair(:,:)
    integer(int2),    pointer :: id_g2l(:,:)
    integer,      allocatable :: wc_type(:)
    logical,      allocatable :: no_wild(:)
    integer,          pointer :: mol_impr_list(:,:)
    character(6),     pointer :: mol_cls_name(:)
    integer,          pointer :: impr_pbc(:,:,:)


    mol_impr_list  => molecule%impr_list
    mol_cls_name   => molecule%atom_cls_name

    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l

    improper  => enefunc%num_improper

    nimpr     = molecule%num_impropers
    nimpr_p   = par%num_impropers


    ! check usage of wild card
    ! but spica has no wild card, so this can be simplified later
    allocate(wc_type(nimpr_p), no_wild(nimpr))

    do i = 1, nimpr_p

      if ((par%impr_atom_cls(1,i) .ne. WildCard) .and. &
          (par%impr_atom_cls(2,i) .ne. WildCard) .and. &
          (par%impr_atom_cls(3,i) .ne. WildCard) .and. &
          (par%impr_atom_cls(4,i) .ne. WildCard)) then

        ! A-B-C-D type
        wc_type(i) = 0

      else if ((par%impr_atom_cls(1,i) .eq. WildCard) .and. &
               (par%impr_atom_cls(2,i) .ne. WildCard) .and. &
               (par%impr_atom_cls(3,i) .ne. WildCard) .and. &
               (par%impr_atom_cls(4,i) .ne. WildCard)) then

        ! X-B-C-D type
        wc_type(i) = 1

      else if ((par%impr_atom_cls(1,i) .eq. WildCard) .and. &
               (par%impr_atom_cls(2,i) .eq. WildCard) .and. &
               (par%impr_atom_cls(3,i) .ne. WildCard) .and. &
               (par%impr_atom_cls(4,i) .ne. WildCard)) then

        ! X-X-C-D type
        wc_type(i) = 2

      else if ((par%impr_atom_cls(1,i) .ne. WildCard) .and. &
               (par%impr_atom_cls(2,i) .eq. WildCard) .and. &
               (par%impr_atom_cls(3,i) .eq. WildCard) .and. &
               (par%impr_atom_cls(4,i) .ne. WildCard)) then

        ! A-X-X-D type
        wc_type(i) = 3

      else
        call error_msg('Setup_Enefunc_Impr> Undefined Wild Card')

      end if

    end do

    ! find number of interactions
    !
    do dupl = 1, domain%num_duplicate

      ioffset = (dupl-1) * enefunc%table%num_all

      do i = 1, nimpr

        no_wild(i) = .false.
  
        list(1:4) = mol_impr_list(1:4,i)
        ci1 = mol_cls_name(list(1))
        ci2 = mol_cls_name(list(2))
        ci3 = mol_cls_name(list(3))
        ci4 = mol_cls_name(list(4))

        lists(1:4) = list(1:4) + ioffset

        icel1 = id_g2l(1,lists(1))
        icel2 = id_g2l(1,lists(4))

        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1,icel2)

          if (icel_local >= 1 .and. icel_local <= ncel) then

            ! A-B-C-D type
            !
            do j = 1, nimpr_p
              if (wc_type(j) == 0) then
                if (((ci1 .eq. par%impr_atom_cls(1,j)) .and. &
                     (ci2 .eq. par%impr_atom_cls(2,j)) .and. &
                     (ci3 .eq. par%impr_atom_cls(3,j)) .and. &
                     (ci4 .eq. par%impr_atom_cls(4,j))) .or. &
                    ((ci1 .eq. par%impr_atom_cls(4,j)) .and. &
                     (ci2 .eq. par%impr_atom_cls(3,j)) .and. &
                     (ci3 .eq. par%impr_atom_cls(2,j)) .and. &
                     (ci4 .eq. par%impr_atom_cls(1,j)))) then
  
                  improper(icel_local) = improper(icel_local) + 1
                  no_wild(i) = .true.

                  exit

                end if
              end if
            end do
  
            ! X-B-C-D type
            !
            if (.not.no_wild(i)) then
              do j = 1, nimpr_p
                if (wc_type(j) == 1) then
                  if (((ci2 .eq. par%impr_atom_cls(2,j)) .and. &
                       (ci3 .eq. par%impr_atom_cls(3,j)) .and. &
                       (ci4 .eq. par%impr_atom_cls(4,j))) .or. &
                      ((ci2 .eq. par%impr_atom_cls(4,j)) .and. &
                       (ci3 .eq. par%impr_atom_cls(3,j)) .and. &
                       (ci4 .eq. par%impr_atom_cls(2,j)))) then
  
                    improper(icel_local) = improper(icel_local) + 1
                    no_wild(i) = .true.
                    
                    exit

                  end if
                end if
              end do
            end if
    
            ! X-X-C-D type
            !
            if (.not.no_wild(i)) then
              do j = 1, nimpr_p
                if (wc_type(j) == 2) then
                  if (((ci3 .eq. par%impr_atom_cls(3,j)) .and. &
                       (ci4 .eq. par%impr_atom_cls(4,j))) .or. &
                      ((ci3 .eq. par%impr_atom_cls(4,j)) .and. &
                       (ci4 .eq. par%impr_atom_cls(3,j)))) then
  
                    improper(icel_local) = improper(icel_local) + 1
                    no_wild(i) = .true.

                    exit
  
                  end if
                end if
              end do
            end if
  
            ! A-X-X-D type
            !
            if (.not.no_wild(i)) then
              do j = 1, nimpr_p
                if (wc_type(j) == 3) then
                  if (((ci1 .eq. par%impr_atom_cls(1,j)) .and. &
                       (ci4 .eq. par%impr_atom_cls(4,j))) .or. &
                      ((ci1 .eq. par%impr_atom_cls(4,j)) .and. &
                       (ci4 .eq. par%impr_atom_cls(1,j)))) then
      
                    improper(icel_local) = improper(icel_local) + 1
                    no_wild(i) = .true.

                    exit
  
                  end if
                end if
              end do
            end if

            if (.not.no_wild(i)) &
              write(MsgOut,*) &
                'Setup_Enefunc_Impr> Unknown IMPR type. [', &
                ci1, ']-[', ci2, ']-[', ci3, ']-[', ci4, '] (ERROR)'
  
          end if
        end if

      end do

    end do

    found = 0
    do i = 1, ncel
      found = max(found, improper(i))
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, found, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
#endif
    MaxImpr = MaxImpr * 2
    ImprMove = MaxImpr / 3 ! yyzhang for small testing system, just set ImprMove = MaxImpr
    if (main_rank) then
      write(MsgOut,*) 'sp_enefunc_str    ::MaxImpr      : ', MaxImpr
      write(MsgOut,*) 'sp_enefunc_str    ::ImprMove     : ', ImprMove
      write(MsgOut,*) ' '
    end if
    call alloc_enefunc(enefunc, EneFuncImpr, ncel, ncel)
    ilist     => enefunc%impr_list
    force     => enefunc%impr_force_const
    phase     => enefunc%impr_phase
    impr_pbc  => enefunc%impr_pbc
    improper(1:ncel) = 0

    ! set parameter
    !
    do dupl = 1, domain%num_duplicate

      ioffset = (dupl-1) * enefunc%table%num_all

      do i = 1, nimpr

        no_wild(i) = .false.
  
        list(1:4) = mol_impr_list(1:4,i)
        ci1 = mol_cls_name(list(1))
        ci2 = mol_cls_name(list(2))
        ci3 = mol_cls_name(list(3))
        ci4 = mol_cls_name(list(4))

        lists(1:4) = list(1:4) + ioffset

        icel1 = id_g2l(1,lists(1))
        icel2 = id_g2l(1,lists(4))

        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1,icel2)

          if (icel_local >= 1 .and. icel_local <= ncel) then

            ! A-B-C-D type
            !
            do j = 1, nimpr_p
              if (wc_type(j) == 0) then
                if (((ci1 .eq. par%impr_atom_cls(1,j)) .and. &
                     (ci2 .eq. par%impr_atom_cls(2,j)) .and. &
                     (ci3 .eq. par%impr_atom_cls(3,j)) .and. &
                     (ci4 .eq. par%impr_atom_cls(4,j))) .or. &
                    ((ci1 .eq. par%impr_atom_cls(4,j)) .and. &
                     (ci2 .eq. par%impr_atom_cls(3,j)) .and. &
                     (ci3 .eq. par%impr_atom_cls(2,j)) .and. &
                     (ci4 .eq. par%impr_atom_cls(1,j)))) then
  
                  improper(icel_local) = improper(icel_local) + 1

                  ilist(1:4,improper(icel_local),icel_local) = list(1:4)
                  force (improper(icel_local),icel_local) = &
                       par%impr_force_const(j)
                  phase (improper(icel_local),icel_local) = &
                       par%impr_phase(j) * RAD
                  impr_pbc(1:3,improper(icel_local),icel_local) = 13 

                  no_wild(i) = .true.

                  exit

                end if
              end if
            end do
  
            ! X-B-C-D type
            !
            if (.not.no_wild(i)) then
              do j = 1, nimpr_p
                if (wc_type(j) == 1) then
                  if (((ci2 .eq. par%impr_atom_cls(2,j)) .and. &
                       (ci3 .eq. par%impr_atom_cls(3,j)) .and. &
                       (ci4 .eq. par%impr_atom_cls(4,j))) .or. &
                      ((ci2 .eq. par%impr_atom_cls(4,j)) .and. &
                       (ci3 .eq. par%impr_atom_cls(3,j)) .and. &
                       (ci4 .eq. par%impr_atom_cls(2,j)))) then
  
                    improper(icel_local) = improper(icel_local) + 1

                    ilist(1:4,improper(icel_local),icel_local) = list(1:4)
                    force (improper(icel_local),icel_local) = &
                         par%impr_force_const(j)
                    phase (improper(icel_local),icel_local) = &
                         par%impr_phase(j) * RAD
                    impr_pbc(1:3,improper(icel_local),icel_local) = 13 

                    no_wild(i) = .true.
                    
                    exit

                  end if
                end if
              end do
            end if
    
            ! X-X-C-D type
            !
            if (.not.no_wild(i)) then
              do j = 1, nimpr_p
                if (wc_type(j) == 2) then
                  if (((ci3 .eq. par%impr_atom_cls(3,j)) .and. &
                       (ci4 .eq. par%impr_atom_cls(4,j))) .or. &
                      ((ci3 .eq. par%impr_atom_cls(4,j)) .and. &
                       (ci4 .eq. par%impr_atom_cls(3,j)))) then
  
                    improper(icel_local) = improper(icel_local) + 1

                    ilist(1:4,improper(icel_local),icel_local) = list(1:4)
                    force (improper(icel_local),icel_local) = &
                         par%impr_force_const(j)
                    phase (improper(icel_local),icel_local) = &
                         par%impr_phase(j) * RAD
                    impr_pbc(1:3,improper(icel_local),icel_local) = 13 

                    no_wild(i) = .true.

                    exit
  
                  end if
                end if
              end do
            end if
  
            ! A-X-X-D type
            !
            if (.not.no_wild(i)) then
              do j = 1, nimpr_p
                if (wc_type(j) == 3) then
                  if (((ci1 .eq. par%impr_atom_cls(1,j)) .and. &
                       (ci4 .eq. par%impr_atom_cls(4,j))) .or. &
                      ((ci1 .eq. par%impr_atom_cls(4,j)) .and. &
                       (ci4 .eq. par%impr_atom_cls(1,j)))) then
      
                    improper(icel_local) = improper(icel_local) + 1

                    ilist(1:4,improper(icel_local),icel_local) = list(1:4)
                    force (improper(icel_local),icel_local) = &
                         par%impr_force_const(j)
                    phase (improper(icel_local),icel_local) = &
                         par%impr_phase(j) * RAD
                    impr_pbc(1:3,improper(icel_local),icel_local) = 13 

                    no_wild(i) = .true.

                    exit
  
                  end if
                end if
              end do
            end if

  
          end if
        end if

      end do

    end do


    deallocate(wc_type, no_wild)

    found = 0
    do i = 1, ncel
      found = found + improper(i)
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_impr_all, 1, mpi_integer, mpi_sum, &
                       mpi_comm_country, ierror)
#else
    enefunc%num_impr_all = found
#endif

    if (enefunc%num_impr_all < nimpr*domain%num_duplicate) &
        call error_msg( &
        'Setup_Enefunc_Impr> Some improper paremeters are missing.')

    return

  end subroutine setup_enefunc_impr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_nonb
  !> @brief        define NON-BOND term in potential energy function
  !! @authors      YS, JJ, TI
  !! @param[in]    par         : CHARMM PAR information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    constraints : constraints information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_nonb(par, molecule, constraints, domain, enefunc)

    ! formal arguments
    type(s_par),             intent(in)    :: par
    type(s_molecule),        intent(in)    :: molecule
    type(s_constraints),     intent(in)    :: constraints
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: eps14, rmin14, eps, rmin, lamda_i, lamda_j
    real(dp)                 :: vdw_self1, vdw_self2
    integer                  :: i, j, k, ix, jx, kk
    integer                  :: nonb_p, nbfx_p, cls_local, ncel
    character(6)             :: ci1, ci2

    integer,  allocatable    :: nonb_atom_cls(:), check_cls(:)
    integer,  allocatable    :: atmcls_map_g2l(:), atmcls_map_l2g(:)
    integer,  allocatable    :: nonb_lj(:,:)
    real(wp), allocatable    :: nb14_lj6(:,:), nb14_lj12(:,:)
    real(wp), allocatable    :: nonb_lj6(:,:), nonb_lj12(:,:)
    real(wp), allocatable    :: lj_coef(:,:)


    enefunc%num_atom_cls = par%num_atom_cls

    ELECOEF          = ELECOEF_CHARMM

    ! set lennard-jones parameters
    !
    nonb_p = enefunc%num_atom_cls

    allocate(nonb_atom_cls(nonb_p),     &
             check_cls(nonb_p),         &
             atmcls_map_g2l(nonb_p),    &
             atmcls_map_l2g(nonb_p),    &
             nb14_lj6 (nonb_p, nonb_p), &
             nb14_lj12(nonb_p, nonb_p), &
             nonb_lj6 (nonb_p, nonb_p), &
             nonb_lj12(nonb_p, nonb_p), &
             nonb_lj  (nonb_p, nonb_p), &
             lj_coef(2,nonb_p))

    nonb_atom_cls(1:nonb_p)       = 0
    check_cls(1:nonb_p)           = 0
    nonb_lj  (1:nonb_p, 1:nonb_p) = 0
    nb14_lj6 (1:nonb_p, 1:nonb_p) = 0.0_wp
    nb14_lj12(1:nonb_p, 1:nonb_p) = 0.0_wp
    nonb_lj6 (1:nonb_p, 1:nonb_p) = 0.0_wp
    nonb_lj12(1:nonb_p, 1:nonb_p) = 0.0_wp
    lj_coef  (1:2,1:nonb_p)       = 0.0_wp 

    do i = 1, nonb_p

      lj_coef(1,i) = abs(par%nonb_eps(i))
      lj_coef(2,i) = par%nonb_rmin(i)

      do j = 1, nonb_p

        ! combination rule
        eps14  = sqrt(par%nonb_eps_14(i) * par%nonb_eps_14(j))
        rmin14 = par%nonb_rmin_14(i) + par%nonb_rmin_14(j)
        eps    = sqrt(par%nonb_eps(i) * par%nonb_eps(j))
        rmin   = par%nonb_rmin(i) + par%nonb_rmin(j)

        ! set parameters
        nb14_lj12(i,j) = eps14 * (rmin14 ** 12)
        nb14_lj6(i,j)  = 2.0_wp * eps14 * (rmin14 ** 6)
        nonb_lj12(i,j) = eps * (rmin ** 12)
        nonb_lj6(i,j)  = 2.0_wp * eps * (rmin ** 6)

      end do
    end do

    ! overwrite lennard-jones parameters by NBFIX parameters
    !
    nbfx_p = par%num_nbfix

    do k = 1, nbfx_p
      ci1 = par%nbfi_atom_cls(1,k)
      ci2 = par%nbfi_atom_cls(2,k)
      do i = 1, nonb_p
        do j = 1, nonb_p
          if ((ci1 == par%nonb_atom_cls(i)  .and. &
               ci2 == par%nonb_atom_cls(j)) .or.  &
              (ci2 == par%nonb_atom_cls(i) .and. &
               ci1 == par%nonb_atom_cls(j))) then

            ! combination rule
            !
            eps14  = abs(par%nbfi_eps_14 (k)) !TODO CHECK
            rmin14 = par%nbfi_rmin_14(k)
            eps    = abs(par%nbfi_eps    (k))
            rmin   = par%nbfi_rmin   (k)

            ! nonb_lj(i,j) = par%nbfi_repuls(k)

            if (par%nbfi_repuls(k)==12 .and. par%nbfi_attract(k)==4) then
              nonb_lj(i,j) = 12
              nb14_lj12(i,j) = 1.5_wp * sqrt(3.0_wp)*eps14 * (rmin14**12)
              nb14_lj6 (i,j) = 1.5_wp * sqrt(3.0_wp)*eps14 * (rmin14**4)
              nonb_lj12(i,j) = 1.5_wp * sqrt(3.0_wp)*eps   * (rmin**12)
              nonb_lj6 (i,j) = 1.5_wp * sqrt(3.0_wp)*eps   * (rmin**4)
            else if (par%nbfi_repuls(k)==9 .and. par%nbfi_attract(k)==6) then
              nonb_lj(i,j) = 9
              nb14_lj12(i,j) = 6.75_wp * eps14 * (rmin14**9)
              nb14_lj6 (i,j) = 6.75_wp * eps14 * (rmin14**6)
              nonb_lj12(i,j) = 6.75_wp * eps   * (rmin**9)
              nonb_lj6 (i,j) = 6.75_wp * eps   * (rmin**6)
            else  if(par%nbfi_repuls(k)==12  .and. par%nbfi_attract(k)==5 ) then
              nonb_lj(i,j) = 5
              nb14_lj12(i,j) = 3.203779841*eps14 * (rmin14 ** 12)
              nb14_lj6 (i,j) = 3.203779841*eps14 * (rmin14 ** 5)
              nonb_lj12(i,j) = 3.203779841*eps * (rmin ** 12)
              nonb_lj6 (i,j) = 3.203779841*eps * (rmin ** 5)
            else if (par%nbfi_repuls(k)==-1 .and. par%nbfi_attract(k)==-1) then
              nonb_lj(i,j) = -1
              nb14_lj12(i,j) = eps14 * (rmin14**12)
              nb14_lj6 (i,j) = 2.0_wp * eps14 * (rmin14**6)
              nonb_lj12(i,j) = eps * (rmin**12)
              nonb_lj6 (i,j) = 2.0_wp * eps * (rmin**6)
            else
              call error_msg( &
                 'Setup_Enefunc_NBFIX> unavilable functional type is '//&
                 'defined: "'//trim(molecule%atom_cls_name(i))//'"')
            end if

          end if
        end do
      end do
    end do

    ! check the usage of atom class
    !
    do i = 1, molecule%num_atoms
      k = molecule%atom_cls_no(i)
      if (k < 1) then
        call error_msg( &
        'Setup_Enefunc_Nonb> atom class is not defined: "'&
        //trim(molecule%atom_cls_name(i))//'"')
      end if
      check_cls(k) = check_cls(k) + 1
    end do

    k = 0
    do i = 1, nonb_p
      if (check_cls(i) > 0) then
        k = k + 1
        atmcls_map_g2l(i) = k
        atmcls_map_l2g(k) = i
        ! else if (check_cls(i) < 0) then
        !    write(*,*) "debug check_cls(i)", i
        ! write(*,*)"classnum k",k
      end if
    end do
    cls_local = k
    max_class = cls_local

    call alloc_enefunc(enefunc, EneFuncNbon, cls_local)

    do i = 1, cls_local
      ix = atmcls_map_l2g(i)
      enefunc%nonb_lj6_factor(i)  = sqrt(2.0_wp)*8.0_wp*sqrt(lj_coef(1,ix))  &
                                   *lj_coef(2,ix)**3
      enefunc%nonb_atom_cls_no(i) = check_cls(ix)
    end do

    do i = 1, cls_local
      ix = atmcls_map_l2g(i)
      do j = 1, cls_local
        jx = atmcls_map_l2g(j)
        enefunc%nonb_ljtype(i,j) = nonb_lj  (ix,jx)
        enefunc%nb14_lj12  (i,j) = nb14_lj12(ix,jx)
        enefunc%nb14_lj6   (i,j) = nb14_lj6 (ix,jx)
        enefunc%nonb_lj12  (i,j) = nonb_lj12(ix,jx)
        enefunc%nonb_lj6   (i,j) = nonb_lj6 (ix,jx)
      end do
    end do

    ! Sum of lamda_i*lamda_j and lamda_i*lamda_i
    !
    vdw_self1 = 0.0_dp
    vdw_self2 = 0.0_dp
    k = 0
    do i = 1, cls_local
      ix = enefunc%nonb_atom_cls_no(i)
      lamda_i = enefunc%nonb_lj6_factor(i)
      vdw_self2 = vdw_self2 + ix*lamda_i*lamda_i
      do j = 1, cls_local 
        jx = enefunc%nonb_atom_cls_no(j)
        lamda_j = enefunc%nonb_lj6_factor(j)
        vdw_self1 = vdw_self1 + ix*jx*lamda_i*lamda_j
      end do
    end do
    enefunc%pme_dispersion_self1 = real(vdw_self1,wp)
    enefunc%pme_dispersion_self2 = real(vdw_self2,wp)

    ! update domain information
    !
    do i = 1, domain%num_cell_local+domain%num_cell_boundary
      do ix = 1, domain%num_atom(i)
        domain%atom_cls_no(ix,i) = &
           atmcls_map_g2l(domain%atom_cls_no(ix,i))
      end do
    end do
    if (constraints%water_type == TIP2) then
!     domain%water%atom_cls_no(1:2)  &
!       = atmcls_map_g2l(domain%water%atom_cls_no(1:2))
      do i = 1, 3
        if (enefunc%table%atom_cls_no_1(i) /= 0) then
          enefunc%table%atom_cls_no_1(i)  &
            = atmcls_map_g2l(enefunc%table%atom_cls_no_1(i))
          enefunc%table%atom_cls_no_2(i)  &
            = atmcls_map_g2l(enefunc%table%atom_cls_no_2(i))
        end if
      end do
    else if (constraints%water_type == TIP1) then
      domain%water%atom_cls_no(1) &
        = atmcls_map_g2l(domain%water%atom_cls_no(1))
    end if

    deallocate(nonb_atom_cls,  &
               check_cls,      &
               atmcls_map_g2l, &
               atmcls_map_l2g, &
               nb14_lj6,       &
               nb14_lj12,      &
               nonb_lj6,       &
               nonb_lj12,      & 
               lj_coef)

    enefunc%num_atom_cls = cls_local

    ! treatment for 1-2, 1-3, 1-4 interactions
    !
    ncel   = domain%num_cell_local

    call alloc_enefunc(enefunc, EneFuncNonb,     ncel, maxcell_near)
    call alloc_enefunc(enefunc, EneFuncNonbList, ncel, maxcell_near)

    call count_nonb_excl_spica(.true., constraints, domain, enefunc)

    return

  end subroutine setup_enefunc_nonb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    count_nonb_excl_spica
  !> @brief        exclude 1-2, 1-3 interactions and constraints
  !! @authors      JJ
  !! @param[in]    first       : flag for first call or not
  !! @param[in]    constraints : constraints information   
  !! @param[inout] domain      : structure of domain
  !! @param[inout] enefunc     : structure of enefunc
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine count_nonb_excl_spica(first, constraints, domain, enefunc)

    ! formal arguments
    logical,                     intent(in)    :: first
    type(s_constraints), target, intent(in)    :: constraints
    type(s_domain),      target, intent(inout) :: domain
    type(s_enefunc),     target, intent(inout) :: enefunc

    ! local variables
    integer                  :: ncell, ncell_local, i, ii, ix, k, i1, i2, i3
    integer                  :: icel, icel1, icel2
    integer                  :: ic, j, ih, ij, index(4)
    integer                  :: num_excl, num_nb14, id
    integer                  :: found1, found2
    integer                  :: fkind
    integer                  :: list1, list2, start_i, ncharge

    integer,         pointer :: natom(:), start_atom(:), nwater(:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: water_list(:,:,:)
    integer,         pointer :: nbond(:), bond_list(:,:,:)
    integer,         pointer :: nenm(:), enm_list(:,:,:)
    integer,         pointer :: nangle(:), angl_list(:,:,:)
    integer,         pointer :: ndihedral(:), dihe_list(:,:,:)
    integer,         pointer :: nimproper(:), impr_list(:,:,:) ! yyzhang def spica improper coplanar
    integer,         pointer :: nexcl(:), excl_list(:,:,:)
    integer,         pointer :: n14(:), n14_list(:,:,:)
    integer,         pointer :: cell_pairlist2(:,:)
    integer(int2),   pointer :: cell_pair(:,:)
    integer,         pointer :: nonb_excl_list(:,:,:)
    integer,         pointer :: nb14_calc_list(:,:,:)
    integer(1),      pointer :: exclusion_mask(:,:,:), exclusion_mask1(:,:,:)
    integer,         pointer :: sc_calc_list(:,:)
    integer,         pointer :: num_nonb_excl(:), num_nb14_calc(:)
    integer,         pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)
    real(wp),        pointer :: nb14_qq_scale(:,:), nb14_lj_scale(:,:)
    real(wp),        pointer :: dihe_scnb(:), dihe_scee(:) ! yyzhang readcode not used ?
    real(wp),        pointer :: impr_scnb(:), impr_scee(:) ! yyzhang def improper coplanar not used
    real(wp),        pointer :: charge(:,:)


    natom           => domain%num_atom
    nwater          => domain%num_water
    start_atom      => domain%start_atom
    water_list      => domain%water_list
    id_g2l          => domain%id_g2l
    cell_pair       => domain%cell_pairlist1
    cell_pairlist2  => domain%cell_pairlist2
    charge          => domain%charge

    nbond           => enefunc%num_bond
    nenm            => enefunc%num_enm 
    nangle          => enefunc%num_angle
    ndihedral       => enefunc%num_dihedral
    nimproper       => enefunc%num_improper ! yyzhang def improper coplanar
    bond_list       => enefunc%bond_list
    enm_list        => enefunc%enm_list
    angl_list       => enefunc%angle_list
    dihe_list       => enefunc%dihe_list
    impr_list       => enefunc%impr_list ! yyzhang def improper coplanar
    nexcl           => enefunc%num_excl  
    n14             => enefunc%num_excl14
    excl_list       => enefunc%excl_list
    n14_list        => enefunc%excl14_list
    nonb_excl_list  => enefunc%nonb_excl_list
    nb14_calc_list  => enefunc%nb14_calc_list
    num_nonb_excl   => enefunc%num_nonb_excl
    num_nb14_calc   => enefunc%num_nb14_calc
    exclusion_mask  => enefunc%exclusion_mask
    exclusion_mask1 => enefunc%exclusion_mask1

    ncell_local = domain%num_cell_local
    ncell       = domain%num_cell_local + domain%num_cell_boundary

    domain%max_num_atom = 0
    do i = 1, ncell
      domain%max_num_atom = max(domain%max_num_atom,domain%num_atom(i))
    end do

    ! initialization
    !
    num_nonb_excl(1:ncell_local)  = 0
    num_nb14_calc(1:ncell_local)  = 0

    start_atom(1:ncell) = 0
    k = 0
    do i = 1, ncell-1
      k = k + natom(i)
      start_atom(i+1) = k
    end do
    k = k + natom(ncell)
    domain%num_atom_domain = k

    ! exclude 1-2 interaction
    !
    !$omp parallel default(shared)                                     &
    !$omp private(id, i, ix, icel1, icel2, icel, i1, i2, i3, num_excl, &
    !$omp         k, num_nb14, ic, j, ih, list1, list2, ii, ij, fkind, &
    !$omp         index)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell_local, nthread
      k = natom(i)
      exclusion_mask1(1:k,1:k,i) = 1
      do i1 = 1, k
        exclusion_mask1(i1,i1,i) = 0
      end do
    end do
    do ij = id+1, maxcell_near, nthread
      i = cell_pair(1,ij)
      j = cell_pair(2,ij)
      i1 = max(natom(i),natom(j))
      exclusion_mask(1:i1,1:i1,ij) = 1
    end do
    !$omp barrier

    if (enefunc%excl_level > 0) then

      do i = id+1, ncell_local, nthread

        if (constraints%water_type == TIP2) then
          do ix = 1, nwater(i)
            list1 = water_list(1,ix,i)
            list2 = water_list(2,ix,i)
            if (abs(charge(list1,i)) > EPS .and. &
                abs(charge(list2,i)) > EPS) then
              num_excl = num_nonb_excl(i) + 1
              num_nonb_excl(i) = num_excl
              nonb_excl_list(1,num_excl,i) = start_atom(i) + list1
              nonb_excl_list(2,num_excl,i) = start_atom(i) + list2
            end if
            exclusion_mask1(list1,list2,i) = 0
            exclusion_mask1(list2,list1,i) = 0
          end do
        end if

        do ix = 1, nbond(i)

          list1 = bond_list(1,ix,i)
          list2 = bond_list(2,ix,i)
          icel1 = id_g2l(1,list1)
          icel2 = id_g2l(1,list2)
          i1    = id_g2l(2,list1)
          i2    = id_g2l(2,list2)
          if (abs(charge(i1,icel1)) > EPS .and. &
              abs(charge(i2,icel2)) > EPS) then
            num_excl = num_nonb_excl(i) + 1
            num_nonb_excl(i) = num_excl
            nonb_excl_list(1,num_excl,i) = start_atom(icel1) + i1
            nonb_excl_list(2,num_excl,i) = start_atom(icel2) + i2
          end if

          if (icel1 /= icel2) then
            icel  = cell_pairlist2(icel1,icel2)
            if (icel1 < icel2) then
              exclusion_mask(i2,i1,icel) = 0
            else if (icel1 > icel2) then
              exclusion_mask(i1,i2,icel) = 0
            end if
          else
            if (i1 < i2) then
              exclusion_mask1(i2,i1,i) = 0
            else if (i1 > i2) then
              exclusion_mask1(i1,i2,i) = 0
            end if
          end if

        end do

        do ix = 1, nenm(i)

          list1 = enm_list(1,ix,i)
          list2 = enm_list(2,ix,i)
          icel1 = id_g2l(1,list1)
          icel2 = id_g2l(1,list2)
          i1    = id_g2l(2,list1)
          i2    = id_g2l(2,list2)
          icel  = cell_pairlist2(icel1,icel2)

          if (icel < maxcell_near) then

            if (abs(charge(i1,icel1)) > EPS .and. &
                abs(charge(i2,icel2)) > EPS) then
              num_excl = num_nonb_excl(i) + 1
              num_nonb_excl(i) = num_excl
              nonb_excl_list(1,num_excl,i) = start_atom(icel1) + i1
              nonb_excl_list(2,num_excl,i) = start_atom(icel2) + i2
            end if

            if (icel1 /= icel2) then
              icel  = cell_pairlist2(icel1,icel2)
              if (icel1 < icel2) then
                exclusion_mask(i2,i1,icel) = 0
              else if (icel1 > icel2) then
                exclusion_mask(i1,i2,icel) = 0
              end if
            else
              if (i1 < i2) then
                exclusion_mask1(i2,i1,i) = 0
              else if (i1 > i2) then
                exclusion_mask1(i1,i2,i) = 0
              end if
            end if
          end if
        end do

      end do
    end if

    !$omp barrier

    if (enefunc%excl_level > 1) then

      ! exclude 1-3 interaction
      !
      do i = id+1, ncell_local, nthread
        do ix = 1, nangle(i)

          list1 = angl_list(1,ix,i)
          list2 = angl_list(3,ix,i)
          icel1 = id_g2l(1,list1)
          icel2 = id_g2l(1,list2)
          i1    = id_g2l(2,list1)
          i2    = id_g2l(2,list2)

          if (icel1 /= icel2) then
            icel  = cell_pairlist2(icel1,icel2)
            if (icel1 < icel2) then
              if (exclusion_mask(i2,i1,icel)==1) then
                if (abs(charge(i1,icel1)) > EPS .and. &
                    abs(charge(i2,icel2)) > EPS) then
                  num_excl = num_nonb_excl(i) + 1
                  num_nonb_excl(i) = num_excl
                  nonb_excl_list(1,num_excl,i) = start_atom(icel1) + i1
                  nonb_excl_list(2,num_excl,i) = start_atom(icel2) + i2
                end if
                exclusion_mask(i2,i1,icel) = 0
              end if
            else if (icel1 > icel2) then
              if (exclusion_mask(i1,i2,icel)==1) then
                if (abs(charge(i1,icel1)) > EPS .and. &
                    abs(charge(i2,icel2)) > EPS) then
                  num_excl = num_nonb_excl(i) + 1
                  num_nonb_excl(i) = num_excl
                  nonb_excl_list(1,num_excl,i) = start_atom(icel1) + i1
                  nonb_excl_list(2,num_excl,i) = start_atom(icel2) + i2
                end if
                exclusion_mask(i1,i2,icel) = 0
              end if
            end if
          else
            if (i1 < i2) then
              if (exclusion_mask1(i2,i1,i)==1) then
                if (abs(charge(i1,icel1)) > EPS .and. &
                    abs(charge(i2,icel2)) > EPS) then
                  num_excl = num_nonb_excl(i) + 1
                  num_nonb_excl(i) = num_excl
                  nonb_excl_list(1,num_excl,i) = start_atom(icel1) + i1
                  nonb_excl_list(2,num_excl,i) = start_atom(icel2) + i2
                end if
                exclusion_mask1(i2,i1,i) = 0
              end if
            else if (i1 > i2) then
              if (exclusion_mask1(i1,i2,i)==1) then
                if (abs(charge(i1,icel1)) > EPS .and. &
                    abs(charge(i2,icel2)) > EPS) then
                  num_excl = num_nonb_excl(i) + 1
                  num_nonb_excl(i) = num_excl
                  nonb_excl_list(1,num_excl,i) = start_atom(icel1) + i1
                  nonb_excl_list(2,num_excl,i) = start_atom(icel2) + i2
                end if
                exclusion_mask1(i1,i2,i) = 0
              end if
            end if
          end if

        end do
      end do

      do i = id+1, ncell_local, nthread
        do ix = 1, nexcl(i)

          list1 = excl_list(1,ix,i)
          list2 = excl_list(2,ix,i)
          icel1 = id_g2l(1,list1)
          icel2 = id_g2l(1,list2)
          i1    = id_g2l(2,list1)
          i2    = id_g2l(2,list2)

          if (icel1 /= icel2) then
            icel  = cell_pairlist2(icel1,icel2)
            if (icel1 < icel2) then
              if (exclusion_mask(i2,i1,icel)==1) then
                if (abs(charge(i1,icel1)) > EPS .and. &
                    abs(charge(i2,icel2)) > EPS) then
                  num_excl = num_nonb_excl(i) + 1
                  num_nonb_excl(i) = num_excl
                  nonb_excl_list(1,num_excl,i) = start_atom(icel1) + i1
                  nonb_excl_list(2,num_excl,i) = start_atom(icel2) + i2
                end if
                exclusion_mask(i2,i1,icel) = 0
              end if
            else if (icel1 > icel2) then
              if (exclusion_mask(i1,i2,icel)==1) then
                if (abs(charge(i1,icel1)) > EPS .and. &
                    abs(charge(i2,icel2)) > EPS) then
                  num_excl = num_nonb_excl(i) + 1
                  num_nonb_excl(i) = num_excl
                  nonb_excl_list(1,num_excl,i) = start_atom(icel1) + i1
                  nonb_excl_list(2,num_excl,i) = start_atom(icel2) + i2
                end if
                exclusion_mask(i1,i2,icel) = 0
              end if
            end if
          else
            if (i1 < i2) then
              if (exclusion_mask1(i2,i1,i)==1) then
                if (abs(charge(i1,icel1)) > EPS .and. &
                    abs(charge(i2,icel2)) > EPS) then
                  num_excl = num_nonb_excl(i) + 1
                  num_nonb_excl(i) = num_excl
                  nonb_excl_list(1,num_excl,i) = start_atom(icel1) + i1
                  nonb_excl_list(2,num_excl,i) = start_atom(icel2) + i2
                end if
                exclusion_mask1(i2,i1,i) = 0
              end if
            else if (i1 > i2) then
              if (exclusion_mask1(i1,i2,i)==1) then
                if (abs(charge(i1,icel1)) > EPS .and. &
                    abs(charge(i2,icel2)) > EPS) then
                  num_excl = num_nonb_excl(i) + 1
                  num_nonb_excl(i) = num_excl
                  nonb_excl_list(1,num_excl,i) = start_atom(icel1) + i1
                  nonb_excl_list(2,num_excl,i) = start_atom(icel2) + i2
                end if
                exclusion_mask1(i1,i2,i) = 0
              end if
            end if
          end if

        end do
      end do

    end if

    !$omp end parallel

    ! count 1-4 interaction
    !
    if (enefunc%excl_level > 2) then

      !$omp parallel default(shared)                                     &
      !$omp private(id, i, ix, icel1, icel2, icel, i1, i2, list1, list2, &
      !$omp         num_nb14, fkind)
      !
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      ! yyzhang readcode exclude dihedral start
      do i = id+1, ncell_local, nthread
        do ix = 1, ndihedral(i)

          list1 = dihe_list(1,ix,i)
          list2 = dihe_list(4,ix,i)
          icel1 = id_g2l(1,list1)
          icel2 = id_g2l(1,list2)
          i1    = id_g2l(2,list1)
          i2    = id_g2l(2,list2)

          if (icel1 /= icel2) then
            icel  = cell_pairlist2(icel1,icel2)
            if (icel1 < icel2) then
              if (exclusion_mask(i2,i1,icel) == 1) then
                if (abs(charge(i1,icel1)) > EPS .and. &
                    abs(charge(i2,icel2)) > EPS) then
                  num_nb14 = num_nb14_calc(i) + 1
                  num_nb14_calc(i) = num_nb14
                  nb14_calc_list(1,num_nb14,i) = start_atom(icel1) + i1
                  nb14_calc_list(2,num_nb14,i) = start_atom(icel2) + i2
                end if
                exclusion_mask(i2,i1,icel) = 0
              end if
            else if (icel1 > icel2) then
              if (exclusion_mask(i1,i2,icel) == 1) then
                if (abs(charge(i1,icel1)) > EPS .and. &
                    abs(charge(i2,icel2)) > EPS) then
                  num_nb14 = num_nb14_calc(i) + 1
                  num_nb14_calc(i) = num_nb14
                  nb14_calc_list(1,num_nb14,i) = start_atom(icel1) + i1
                  nb14_calc_list(2,num_nb14,i) = start_atom(icel2) + i2
                end if
                exclusion_mask(i1,i2,icel) = 0
              end if
            end if
          else
            if (i1 < i2) then
              if (exclusion_mask1(i2,i1,i) == 1) then
                if (abs(charge(i1,icel1)) > EPS .and. &
                    abs(charge(i2,icel2)) > EPS) then
                  num_nb14 = num_nb14_calc(i) + 1
                  num_nb14_calc(i) = num_nb14
                  nb14_calc_list(1,num_nb14,i) = start_atom(icel1) + i1
                  nb14_calc_list(2,num_nb14,i) = start_atom(icel2) + i2
                end if
                exclusion_mask1(i2,i1,i) = 0
              end if
            else if (i1 > i2) then
              if (exclusion_mask1(i1,i2,i) == 1) then
                if (abs(charge(i1,icel1)) > EPS .and. &
                    abs(charge(i2,icel2)) > EPS) then
                  num_nb14 = num_nb14_calc(i) + 1
                  num_nb14_calc(i) = num_nb14
                  nb14_calc_list(1,num_nb14,i) = start_atom(icel1) + i1
                  nb14_calc_list(2,num_nb14,i) = start_atom(icel2) + i2
                end if
                exclusion_mask1(i1,i2,i) = 0
              end if
            end if
          end if

        end do
    end do
      ! yyzhang readcode exclude dihedral end

      ! yyzhang exclude improper start
      do i = id+1, ncell_local, nthread
        do ix = 1, nimproper(i)

          list1 = impr_list(1,ix,i)
          list2 = impr_list(4,ix,i)
          icel1 = id_g2l(1,list1)
          icel2 = id_g2l(1,list2)
          i1    = id_g2l(2,list1)
          i2    = id_g2l(2,list2)

          if (icel1 /= icel2) then
            icel  = cell_pairlist2(icel1,icel2)
            if (icel1 < icel2) then
              if (exclusion_mask(i2,i1,icel) == 1) then
                if (abs(charge(i1,icel1)) > EPS .and. &
                    abs(charge(i2,icel2)) > EPS) then
                  num_nb14 = num_nb14_calc(i) + 1
                  num_nb14_calc(i) = num_nb14
                  nb14_calc_list(1,num_nb14,i) = start_atom(icel1) + i1
                  nb14_calc_list(2,num_nb14,i) = start_atom(icel2) + i2
                end if
                exclusion_mask(i2,i1,icel) = 0
              end if
            else if (icel1 > icel2) then
              if (exclusion_mask(i1,i2,icel) == 1) then
                if (abs(charge(i1,icel1)) > EPS .and. &
                    abs(charge(i2,icel2)) > EPS) then
                  num_nb14 = num_nb14_calc(i) + 1
                  num_nb14_calc(i) = num_nb14
                  nb14_calc_list(1,num_nb14,i) = start_atom(icel1) + i1
                  nb14_calc_list(2,num_nb14,i) = start_atom(icel2) + i2
                end if
                exclusion_mask(i1,i2,icel) = 0
              end if
            end if
          else
            if (i1 < i2) then
              if (exclusion_mask1(i2,i1,i) == 1) then
                if (abs(charge(i1,icel1)) > EPS .and. &
                    abs(charge(i2,icel2)) > EPS) then
                  num_nb14 = num_nb14_calc(i) + 1
                  num_nb14_calc(i) = num_nb14
                  nb14_calc_list(1,num_nb14,i) = start_atom(icel1) + i1
                  nb14_calc_list(2,num_nb14,i) = start_atom(icel2) + i2
                end if
                exclusion_mask1(i2,i1,i) = 0
              end if
            else if (i1 > i2) then
              if (exclusion_mask1(i1,i2,i) == 1) then
                if (abs(charge(i1,icel1)) > EPS .and. &
                    abs(charge(i2,icel2)) > EPS) then
                  num_nb14 = num_nb14_calc(i) + 1
                  num_nb14_calc(i) = num_nb14
                  nb14_calc_list(1,num_nb14,i) = start_atom(icel1) + i1
                  nb14_calc_list(2,num_nb14,i) = start_atom(icel2) + i2
                end if
                exclusion_mask1(i1,i2,i) = 0
              end if
            end if
          end if

        end do
      end do
      ! yyzhang exclude improper end

      do i = id+1, ncell_local, nthread
        do ix = 1, n14(i)

          list1 = n14_list(1,ix,i)
          list2 = n14_list(2,ix,i)
          icel1 = id_g2l(1,list1)
          icel2 = id_g2l(1,list2)
          i1    = id_g2l(2,list1)
          i2    = id_g2l(2,list2)

          if (icel1 /= icel2) then
            icel  = cell_pairlist2(icel1,icel2)
            if (icel1 < icel2) then
              if (exclusion_mask(i2,i1,icel) == 1) then
                if (abs(charge(i1,icel1)) > EPS .and. &
                    abs(charge(i2,icel2)) > EPS) then
                  num_nb14 = num_nb14_calc(i) + 1
                  num_nb14_calc(i) = num_nb14
                  nb14_calc_list(1,num_nb14,i) = start_atom(icel1) + i1
                  nb14_calc_list(2,num_nb14,i) = start_atom(icel2) + i2
                end if
                exclusion_mask(i2,i1,icel) = 0
              end if
            else if (icel1 > icel2) then
              if (exclusion_mask(i1,i2,icel) == 1) then
                if (abs(charge(i1,icel1)) > EPS .and. &
                    abs(charge(i2,icel2)) > EPS) then
                  num_nb14 = num_nb14_calc(i) + 1
                  num_nb14_calc(i) = num_nb14
                  nb14_calc_list(1,num_nb14,i) = start_atom(icel1) + i1
                  nb14_calc_list(2,num_nb14,i) = start_atom(icel2) + i2
                end if
                exclusion_mask(i1,i2,icel) = 0
              end if
            end if
          else
            if (i1 < i2) then
              if (exclusion_mask1(i2,i1,i) == 1) then
                if (abs(charge(i1,icel1)) > EPS .and. &
                    abs(charge(i2,icel2)) > EPS) then
                  num_nb14 = num_nb14_calc(i) + 1
                  num_nb14_calc(i) = num_nb14
                  nb14_calc_list(1,num_nb14,i) = start_atom(icel1) + i1
                  nb14_calc_list(2,num_nb14,i) = start_atom(icel2) + i2
                end if
                exclusion_mask1(i2,i1,i) = 0
              end if
            else if (i1 > i2) then
              if (exclusion_mask1(i1,i2,i) == 1) then
                if (abs(charge(i1,icel1)) > EPS .and. &
                    abs(charge(i2,icel2)) > EPS) then
                  num_nb14 = num_nb14_calc(i) + 1
                  num_nb14_calc(i) = num_nb14
                  nb14_calc_list(1,num_nb14,i) = start_atom(icel1) + i1
                  nb14_calc_list(2,num_nb14,i) = start_atom(icel2) + i2
                end if
                exclusion_mask1(i1,i2,i) = 0
              end if
            end if
          end if

        end do
      end do
      !$omp end parallel

    end if

    ! Check the total number of exclusion list
    !
    if (first) then

      found1 = 0
      found2 = 0

      do icel = 1, ncell_local
        found1 = found1 + num_nonb_excl(icel)
        found2 = found2 + num_nb14_calc(icel)
      end do

#ifdef HAVE_MPI_GENESIS
      call mpi_reduce(found1, enefunc%num_excl_all, 1, mpi_integer, mpi_sum, &
                      0, mpi_comm_country, ierror)
      call mpi_reduce(found2, enefunc%num_nb14_all, 1, mpi_integer, mpi_sum, &
                      0, mpi_comm_country, ierror)
#else
    enefunc%num_excl_all = found1
    enefunc%num_nb14_all = found2
#endif
    end if

    return

  end subroutine count_nonb_excl_spica

end module sp_enefunc_spica_mod
