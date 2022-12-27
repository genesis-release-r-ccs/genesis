!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_enefunc_gromacs_mod
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

module sp_enefunc_gromacs_mod

  use sp_enefunc_charmm_mod
  use sp_enefunc_restraints_mod
  use sp_enefunc_table_mod
  use sp_energy_mod
  use sp_restraints_str_mod
  use sp_constraints_str_mod
  use sp_enefunc_str_mod
  use sp_energy_str_mod
  use sp_domain_str_mod
  use molecules_str_mod
  use fileio_grotop_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: define_enefunc_gromacs
  private :: setup_enefunc_bond
  private :: setup_enefunc_bond_constraint
  private :: setup_enefunc_angl
  private :: setup_enefunc_angl_constraint
  private :: setup_enefunc_dihe
  private :: setup_enefunc_rb_dihe
  private :: setup_enefunc_impr
  private :: setup_enefunc_nonb
  private :: setup_enefunc_gro_restraints

  ! paramters
  integer, parameter :: WaterIdx(2,3) = reshape((/1,2,1,3,2,3/),shape=(/2,3/))

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc_gromacs
  !> @brief        a driver subroutine for defining potential energy functions
  !! @authors      NT
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    grotop      : GROMACS TOP information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    constraints : constraints information
  !! @param[in]    restraints  : restraints information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc_gromacs(ene_info, grotop, molecule, &
                                    constraints, restraints, domain, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info 
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_constraints),     intent(inout) :: constraints
    type(s_restraints),      intent(in)    :: restraints
    type(s_domain),          intent(inout) :: domain 
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: ncel, ncelb


    ! base
    !
    ncel  = domain%num_cell_local
    ncelb = domain%num_cell_local + domain%num_cell_boundary

    call alloc_enefunc(enefunc, EneFuncBase,     ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncBond,     ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncAngl,     ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncDihe,     ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncRBDihe,   ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncImpr,     ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncBondCell, ncel, ncelb)

    if (.not. constraints%rigid_bond) then

      ! bond
      !
      call setup_enefunc_bond(molecule, grotop, domain, constraints, enefunc)

      ! angle
      !
      call setup_enefunc_angl(molecule, grotop, domain, enefunc)

    else

      ! bond
      !
      call setup_enefunc_bond_constraint(molecule, grotop, domain, &
                                         constraints, enefunc)

      ! angle
      !
      call setup_enefunc_angl_constraint(molecule, grotop, domain, &
                                         constraints, enefunc)

    end if

    ! dihedral
    !
    call setup_enefunc_dihe(molecule, grotop, domain, enefunc)

    ! Ryckaert-Bellemans dihedral
    !
    call setup_enefunc_rb_dihe(molecule, grotop, domain, enefunc)

    ! improper
    !
    call setup_enefunc_impr(molecule, grotop, domain, enefunc)

    ! nonbonded
    !
    call setup_enefunc_nonb(grotop, molecule, constraints, domain, enefunc)

    ! lookup table
    !
    call setup_enefunc_table(ene_info, enefunc)

    ! restraints
    !
    call setup_enefunc_restraints(molecule, restraints, domain, enefunc)

    ! restraints (gromacs)
    !
    call setup_enefunc_gro_restraints(molecule, grotop, domain, enefunc)


    ! write summary of energy function
    !
    if (main_rank) then
      write(MsgOut,'(A)') &
           'Define_Enefunc_Gromacs> Number of Interactions in Each Term'
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           '  bond_ene        = ', enefunc%num_bond_all,        &
           '  angle_ene       = ', enefunc%num_angl_all
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           '  torsion_ene     = ', enefunc%num_dihe_all,        &
           '  rb_torsion_ene  = ', enefunc%num_rb_dihe_all
      write(MsgOut,'(A20,I10)')                                 &
           '  improper_ene    = ', enefunc%num_impr_all
      if (.not. ene_info%table) then
        write(MsgOut,'(A20,I10,A20,I10)')                       &
           '  nb_exclusions   = ', enefunc%num_excl_all,        &
           '  nb14_calc       = ', enefunc%num_nb14_all
      end if
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           ' restraint_groups = ', enefunc%num_restraintgroups, &
           ' restraint_funcs  = ', enefunc%num_restraintfuncs
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine define_enefunc_gromacs

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_bond
  !> @brief        define BOND term for each cell in potential energy function
  !! @authors      NT
  !! @param[in]    molecule    : molecule information
  !! @param[in]    grotop      : GROMACS TOP information
  !! @param[in]    domain      : domain information
  !! @param[in]    constraints : constraints information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond(molecule, grotop, domain, constraints, enefunc)

    ! formal arguments
    type(s_molecule),target, intent(in)    :: molecule
    type(s_grotop),          intent(in)    :: grotop
    type(s_domain),  target, intent(in)    :: domain
    type(s_constraints),     intent(in)    :: constraints
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: water_dist(3), m1, m2
    integer                  :: dupl, ioffset_dupl
    integer                  :: i, j, k, i1, i2, ia, ib, pbc_int
    integer                  :: ioffset, nbond_a
    integer                  :: idx1, idx2, icel1, icel2, icel_local
    character(6)             :: res1, res2
    real(wp)                 :: cwork(3,2), dij(3)

    type(s_grotop_mol), pointer :: gromol
    real(wp),           pointer :: force(:,:), dist(:,:), coord(:,:)
    real(wp),           pointer :: box_size(:)
    integer,            pointer :: bond(:), list(:,:,:)
    integer,            pointer :: ncel
    integer(int2),      pointer :: cell_pair(:,:)
    integer(int2),      pointer :: id_g2l(:,:)
    integer,            pointer :: bond_pbc(:,:)


    coord     => molecule%atom_coord

    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l
    box_size  => domain%system_size

    bond      => enefunc%num_bond
    list      => enefunc%bond_list
    force     => enefunc%bond_force_const
    dist      => enefunc%bond_dist_min
    bond_pbc  => enefunc%bond_pbc

    nbond_a   = 0

    do dupl = 1, domain%num_duplicate

      ioffset_dupl = (dupl-1) * enefunc%table%num_all
      ioffset   = 0

      do i = 1, grotop%num_molss
        gromol => grotop%molss(i)%moltype%mol
        do j = 1, grotop%molss(i)%count

          if (gromol%settles%func == 0) then

            do k = 1, gromol%num_bonds

              idx1 = gromol%bonds(k)%atom_idx1 
              idx2 = gromol%bonds(k)%atom_idx2
              res1 = gromol%atoms(idx1)%residue_name
              res2 = gromol%atoms(idx2)%residue_name
  
              if (res1(1:3) .ne. 'TIP' .and. res1(1:3) .ne. 'WAT' .and. &
                  res1(1:3) .ne. 'SOL' .and. res2(1:3) .ne. 'TIP' .and. &
                  res2(1:3) .ne. 'WAT' .and. res2(1:3) .ne. 'SOL') then

                idx1 = idx1 + ioffset  
                idx2 = idx2 + ioffset  
                ia   = idx1 + ioffset_dupl
                ib   = idx2 + ioffset_dupl
  
                icel1 = id_g2l(1,ia)
                icel2 = id_g2l(1,ib)
  
                if (icel1 /= 0 .and. icel2 /= 0) then
  
                  icel_local = cell_pair(icel1,icel2)
  
                  if (icel_local > 0 .and. icel_local <= ncel) then
  
                    if (bond(icel_local)+1 > MaxBond) &
                      call error_msg('Setup_Enefunc_Bond> Too many bonds.') 
  
                    nbond_a = nbond_a + 1
  
                    bond(icel_local) = bond(icel_local) + 1
                    list (1,bond(icel_local),icel_local) = ia
                    list (2,bond(icel_local),icel_local) = ib
                    force(  bond(icel_local),icel_local) = &
                            gromol%bonds(k)%kb * 0.01_wp * JOU2CAL * 0.5_wp
                    dist (  bond(icel_local),icel_local) = &
                            gromol%bonds(k)%b0 * 10.0_wp
                    cwork(1:3,1) = coord(1:3,idx1)
                    cwork(1:3,2) = coord(1:3,idx2)
                    dij(1:3) = cwork(1:3,1) - cwork(1:3,2)
                    call check_pbc(box_size, dij, pbc_int)
                    bond_pbc(bond(icel_local),icel_local) = pbc_int

                  end if
  
                end if
  
              end if

            end do

          end if

          ioffset = ioffset + gromol%num_atoms

        end do
      end do

    end do

    ! water
    !
    ioffset = 0
    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        if (gromol%settles%func == 0) then

          do k = 1, gromol%num_bonds

            idx1 = gromol%bonds(k)%atom_idx1 
            idx2 = gromol%bonds(k)%atom_idx2
            res1 = gromol%atoms(idx1)%residue_name
            res2 = gromol%atoms(idx2)%residue_name
            m1   = gromol%atoms(idx1)%mass
            m2   = gromol%atoms(idx2)%mass

            if (res1(1:3) .eq. 'TIP' .or. res1(1:3) .eq. 'WAT' .or. &
                res1(1:3) .eq. 'SOL') then
 
              if (m1 /= m2) then
                nbond_a = nbond_a + domain%num_duplicate
                enefunc%table%water_bond_calc_OH = .true.
                enefunc%table%water_bond_calc = .true.
                enefunc%table%OH_bond = gromol%bonds(k)%b0 * 10.0_wp
                enefunc%table%OH_force = &
                    gromol%bonds(k)%kb * 0.01_wp * JOU2CAL * 0.5_wp
              else if (m1 == m2 .and. m1 < LIGHT_ATOM_MASS_LIMIT) then
                nbond_a = nbond_a + domain%num_duplicate
                enefunc%table%HH_bond = gromol%bonds(k)%b0 * 10.0_wp
                enefunc%table%HH_force = &
                    gromol%bonds(k)%kb * 0.01_wp * JOU2CAL * 0.5_wp
                if (enefunc%table%HH_force > EPS) &
                  enefunc%table%water_bond_calc_HH = .true.
              end if
            end if

          end do
        end if
      end do
    end do

    if (.not. enefunc%table%water_bond_calc) then
      do i = 1, grotop%num_molss
        gromol => grotop%molss(i)%moltype%mol
        do j = 1, grotop%molss(i)%count
          if (gromol%settles%func == 1) then
            water_dist = (/gromol%settles%doh * 10.0_wp, &
                           gromol%settles%doh * 10.0_wp, &
                           gromol%settles%dhh * 10.0_wp/)
            enefunc%table%water_bond_calc_OH = .true.
            enefunc%table%OH_bond  = water_dist(1)
            enefunc%table%OH_force = 0.0_wp
            enefunc%table%HH_bond  = water_dist(3)
            enefunc%table%HH_force = 0.0_wp
            enefunc%table%water_bond_calc = .true.
            enefunc%table%water_bond_calc_OH = .true.
            enefunc%table%water_bond_calc_HH = .true.
            nbond_a = nbond_a + 3*domain%num_duplicate
          end if
        end do
      end do
    end if
    
#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(nbond_a, enefunc%num_bond_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_bond_all = nbond_a
#endif

    return
    
  end subroutine setup_enefunc_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_bond_constraint
  !> @brief        define BOND term between heavy atoms
  !! @authors      NT
  !! @param[in]    molecule    : molecule information
  !! @param[in]    grotop      : CHARMM grotop information
  !! @param[in]    domain      : domain information
  !! @param[in]    constraints : constraints information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond_constraint(molecule, grotop, domain, &
                                           constraints, enefunc)

    ! formal arguments
    type(s_molecule),    target, intent(in)    :: molecule
    type(s_grotop),              intent(in)    :: grotop
    type(s_domain),      target, intent(in)    :: domain
    type(s_constraints), target, intent(inout) :: constraints
    type(s_enefunc),     target, intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: water_dist(3)
    integer                  :: dupl, ioffset_dupl, pbc_int
    integer                  :: i, j, k, m, n
    integer                  :: nwat, ioffset, nbond_a, nbond_c
    integer                  :: icel, connect, ih, ih1, ih2, ia, ib
    integer                  :: i1, i2, idx1, idx2, icel1, icel2, icel_local
    integer                  :: wat_bonds, wat_found, nbond_sys
    character(6)             :: res1, res2
    character(6)             :: atm1, atm2
    logical                  :: cl1, cl2
    real(wp)                 :: cwork(3,2), dij(3)

    type(s_grotop_mol), pointer :: gromol
    real(wp),           pointer :: force(:,:), dist(:,:)
    real(wp),           pointer :: coord(:,:), box_size(:)
    real(wip),          pointer :: HGr_bond_dist(:,:,:,:)
    integer,            pointer :: bond(:), list(:,:,:)
    integer,            pointer :: ncel
    integer(int2),      pointer :: cell_pair(:,:)
    integer,            pointer :: id_l2g(:,:)
    integer(int2),      pointer :: id_g2l(:,:)
    integer,            pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)
    integer,            pointer :: bond_pbc(:,:)


    coord         => molecule%atom_coord

    ncel          => domain%num_cell_local
    cell_pair     => domain%cell_pair
    id_g2l        => domain%id_g2l
    id_l2g        => domain%id_l2g
    box_size      => domain%system_size

    HGr_local     => constraints%HGr_local
    HGr_bond_list => constraints%HGr_bond_list
    HGr_bond_dist => constraints%HGr_bond_dist
    connect       =  constraints%connect

    bond          => enefunc%num_bond
    list          => enefunc%bond_list
    force         => enefunc%bond_force_const
    dist          => enefunc%bond_dist_min
    bond_pbc      => enefunc%bond_pbc
    nwat          =  enefunc%table%num_water

    nbond_sys = 0

    nbond_a   = 0
    nbond_c   = 0

    do dupl = 1, domain%num_duplicate

      ioffset_dupl = (dupl-1) * enefunc%table%num_all
      ioffset   = 0

      do i = 1, grotop%num_molss
        gromol => grotop%molss(i)%moltype%mol
        do j = 1, grotop%molss(i)%count

          if (gromol%settles%func == 0) then

            do k = 1, gromol%num_bonds

              i1 = gromol%bonds(k)%atom_idx1
              i2 = gromol%bonds(k)%atom_idx2

              atm1 = gromol%atoms(i1)%atom_type
              atm2 = gromol%atoms(i2)%atom_type

              cl1 = (atm1(1:1) .ne. 'H' .and. atm1(1:1) .ne. 'h')
              cl2 = (atm2(1:1) .ne. 'H' .and. atm2(1:1) .ne. 'h')
              if (constraints%hydrogen_type == ConstraintAtomMass) then
                cl1 = (gromol%atoms(i1)%mass > LIGHT_ATOM_MASS_LIMIT) 
                cl2 = (gromol%atoms(i2)%mass > LIGHT_ATOM_MASS_LIMIT) 
              else if (constraints%hydrogen_type == ConstraintAtomBoth) then
                cl1 = (cl1 .and. &
                    gromol%atoms(i1)%mass > LIGHT_ATOM_MASS_LIMIT) 
                cl2 = (cl2 .and. &
                   gromol%atoms(i2)%mass > LIGHT_ATOM_MASS_LIMIT) 
              end if
              idx1 = i1 + ioffset
              idx2 = i2 + ioffset
              ia   = idx1 + ioffset_dupl
              ib   = idx2 + ioffset_dupl

              if (cl1 .and. cl2) then

                icel1 = id_g2l(1,ia)
                icel2 = id_g2l(1,ib)

                if (icel1 /= 0 .and. icel2 /= 0) then

                  icel_local = cell_pair(icel1,icel2)
  
                  if (icel_local > 0 .and. icel_local <= ncel) then
  
                    if (bond(icel_local)+1 > MaxBond) &
                 call error_msg('Setup_Enefunc_Bond_Constraint> Too many bonds.') 
  
                    nbond_a = nbond_a + 1

                    bond(icel_local) = bond(icel_local) + 1
                    list (1,bond(icel_local),icel_local) = ia
                    list (2,bond(icel_local),icel_local) = ib
                    force(  bond(icel_local),icel_local) = &
                            gromol%bonds(k)%kb * 0.01_wp * JOU2CAL * 0.5_wp
                    dist (  bond(icel_local),icel_local) = &
                            gromol%bonds(k)%b0 * 10.0_wp
                    cwork(1:3,1) = coord(1:3,idx1)
                    cwork(1:3,2) = coord(1:3,idx2)
                    dij(1:3) = cwork(1:3,1) - cwork(1:3,2)
                    call check_pbc(box_size, dij, pbc_int)
                    bond_pbc(bond(icel_local),icel_local) = pbc_int

                  end if

                end if

              else

                icel1 = id_g2l(1,idx1)
                icel2 = id_g2l(1,idx2)

                if (icel1 /= 0 .and. icel2 /= 0) then

                  icel = cell_pair(icel1,icel2)

                  if (icel > 0 .and. icel <= ncel) then
                  
                    do m = 1, connect

                      do n = 1, HGr_local(m,icel)
                        ih1 = id_l2g(HGr_bond_list(1,n,m,icel),icel)
                        do ih = 1, m
                          ih2 = id_l2g(HGr_bond_list(ih+1,n,m,icel),icel)

                          if (ih1 == idx1 .and. ih2 == idx2 .or. &
                            ih2 == idx1 .and. ih1 == idx2) then

                            nbond_c = nbond_c + 1
                            HGr_bond_dist(ih+1,n,m,icel) = &
                                       gromol%bonds(k)%b0 * 10.0_wp
  
                            goto 1

                          end if

                        end do
                      end do

                    end do
                  end if
                end if
1             continue

              end if

            end do

            nbond_sys = nbond_sys + gromol%num_bonds

          else

            water_dist = (/gromol%settles%doh * 10.0_wp, &
                           gromol%settles%doh * 10.0_wp, &
                           gromol%settles%dhh * 10.0_wp/)

            nbond_sys = nbond_sys + 3

          end if

          ioffset   = ioffset   + gromol%num_atoms

        end do
      end do

    end do

    ! for water molecule
    !
    wat_bonds = 0
    wat_found = 0

    if (constraints%fast_water) then

      do i = 1, grotop%num_molss
        gromol => grotop%molss(i)%moltype%mol
        do j = 1, grotop%molss(i)%count

          if (gromol%settles%func == 0) then

            do k = 1, gromol%num_bonds

              i1 = gromol%bonds(k)%atom_idx1
              i2 = gromol%bonds(k)%atom_idx2

              res1 = gromol%atoms(i1)%residue_name
              res2 = gromol%atoms(i2)%residue_name

              if (res1 .eq. constraints%water_model .and. &
                  res2 .eq. constraints%water_model) then

                wat_bonds = wat_bonds+1

                if (k == 1) &
                wat_found = wat_found + 1

              end if

            end do

          else

            wat_bonds = wat_bonds + 3
            wat_found = wat_found + 1

          end if

        end do
      end do

      if (wat_found /= nwat) &
           call error_msg( &
           'Setup_Enefunc_Bond_Constraint> # of water is incorrect')

    end if

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(nbond_a, enefunc%num_bond_all,  1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)

    call mpi_allreduce(nbond_c, constraints%num_bonds, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_bond_all  = nbond_a
    constraints%num_bonds = nbond_c
#endif

    if (constraints%fast_water) then

      if (enefunc%num_bond_all /= &
          (nbond_sys-constraints%num_bonds-wat_bonds*domain%num_duplicate)) then
        call error_msg( &
          'Setup_Enefunc_Bond_Constraint> Some bond paremeters are missing.')
      end if

    else

      if (enefunc%num_bond_all /= &
          (nbond_sys-constraints%num_bonds-3*nwat*domain%num_duplicate))then
        call error_msg( &
          'Setup_Enefunc_Bond_Constraint> Some bond paremeters are missing.')
      end if

    end if

    return
    
  end subroutine setup_enefunc_bond_constraint

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_angl
  !> @brief        define ANGLE term for each cell in potential energy function
  !! @authors      NT
  !! @param[in]    molecule : molecule information
  !! @param[in]    grotop   : GROMACS topology informaiton
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_angl(molecule, grotop, domain, enefunc)

    ! formal arguments
    type(s_molecule),target, intent(in)    :: molecule
    type(s_grotop),  target, intent(in)    :: grotop
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: dupl, ioffset_dupl, pbc_int
    integer                  :: i, j, k, ia, ib, ic
    integer                  :: ioffset, nangl_a
    integer                  :: idx1, idx2, idx3, icel1, icel2, icel_local
    character(6)             :: res1, res2, res3
    real(wp)                 :: cwork(3,3), dij(3)

    type(s_grotop_mol), pointer :: gromol
    real(wp),           pointer :: force(:,:), theta(:,:)
    real(wp),           pointer :: coord(:,:), box_size(:)
    integer,            pointer :: angle(:), list(:,:,:)
    integer,            pointer :: ncel
    integer(int2),      pointer :: cell_pair(:,:)
    integer(int2),      pointer :: id_g2l(:,:)
    integer,            pointer :: angl_pbc(:,:,:)


    coord     => molecule%atom_coord

    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l
    box_size  => domain%system_size

    angle     => enefunc%num_angle
    list      => enefunc%angle_list
    force     => enefunc%angle_force_const
    theta     => enefunc%angle_theta_min
    angl_pbc  => enefunc%angle_pbc

    nangl_a   = 0

    do dupl = 1, domain%num_duplicate

      ioffset_dupl = (dupl-1) * enefunc%table%num_all
      ioffset = 0

      do i = 1, grotop%num_molss
        gromol => grotop%molss(i)%moltype%mol
        do j = 1, grotop%molss(i)%count

          if (gromol%settles%func == 0) then

            do k = 1, gromol%num_angls

              idx1 = gromol%angls(k)%atom_idx1 
              idx2 = gromol%angls(k)%atom_idx2 
              idx3 = gromol%angls(k)%atom_idx3 

              res1 = gromol%atoms(idx1)%residue_name
              res2 = gromol%atoms(idx2)%residue_name
              res3 = gromol%atoms(idx3)%residue_name

              if (res1(1:3) .ne. 'TIP' .and. res1(1:3) .ne. 'WAT' .and. &
                  res1(1:3) .ne. 'SOL') then

                idx1 = idx1 + ioffset 
                idx2 = idx2 + ioffset 
                idx3 = idx3 + ioffset 
                ia   = idx1 + ioffset_dupl
                ib   = idx2 + ioffset_dupl
                ic   = idx3 + ioffset_dupl
                
                icel1 = id_g2l(1,ia)
                icel2 = id_g2l(1,ic)

                if (icel1 /= 0 .and. icel2 /= 0) then

                  icel_local = cell_pair(icel1,icel2)

                  if (icel_local > 0 .and. icel_local <= ncel) then

                    if (angle(icel_local)+1 > MaxAngle) &
                      call error_msg('Setup_Enefunc_Angl> Too many angles.') 

                    nangl_a = nangl_a + 1

                    angle(icel_local) = angle(icel_local) + 1
                    list (1:3,angle(icel_local),icel_local) = (/ia, ib, ic/)
                    force(    angle(icel_local),icel_local) = &
                                      gromol%angls(k)%kt * JOU2CAL * 0.5_wp
                    theta(    angle(icel_local),icel_local) = &
                                      gromol%angls(k)%theta_0 * RAD
                    cwork(1:3,1) = coord(1:3,idx1)
                    cwork(1:3,2) = coord(1:3,idx2)
                    cwork(1:3,3) = coord(1:3,idx3)
                    dij(1:3) = cwork(1:3,1) - cwork(1:3,2)
                    call check_pbc(box_size, dij, pbc_int)
                    angl_pbc(1,angle(icel_local),icel_local) = pbc_int
                    dij(1:3) = cwork(1:3,3) - cwork(1:3,2)
                    call check_pbc(box_size, dij, pbc_int)
                    angl_pbc(2,angle(icel_local),icel_local) = pbc_int
                    dij(1:3) = cwork(1:3,1) - cwork(1:3,3)
                    call check_pbc(box_size, dij, pbc_int)
                    angl_pbc(3,angle(icel_local),icel_local) = pbc_int

                  end if
  
                end if

              end if

            end do

          end if

          ioffset = ioffset + gromol%num_atoms

        end do
      end do

    end do

    ! water
    !
    ioffset = 0
    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        if (gromol%settles%func == 0) then

          do k = 1, gromol%num_angls

            idx1 = gromol%angls(k)%atom_idx1
            idx2 = gromol%angls(k)%atom_idx2
            idx3 = gromol%angls(k)%atom_idx2
            res1 = gromol%atoms(idx1)%residue_name
            res2 = gromol%atoms(idx2)%residue_name
            res3 = gromol%atoms(idx3)%residue_name

            if (res1(1:3) .eq. 'TIP' .or. res1(1:3) .eq. 'WAT' .or. &
                res1(1:3) .eq. 'SOL') then

              nangl_a = nangl_a + domain%num_duplicate
              enefunc%table%water_angle_calc = .true.
              enefunc%table%HOH_angle = gromol%angls(k)%theta_0 * RAD
              enefunc%table%HOH_force = &
                    gromol%angls(k)%kt * JOU2CAL * 0.5_wp
            end if

          end do
        end if
      end do
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(nangl_a, enefunc%num_angl_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_angl_all = nangl_a
#endif

    return

  end subroutine setup_enefunc_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_angl_constraint
  !> @brief        define ANGLE term for each cell in potential energy function
  !                with SETTLE constraint
  !! @authors      NT
  !! @param[in]    molecule    : molecule information
  !! @param[in]    grotop      : GROMACS topology information
  !! @param[in]    domain      : domain information
  !! @param[in]    constraints : constraints information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_angl_constraint(molecule, grotop, domain, &
                                           constraints, enefunc)

    ! formal arguments
    type(s_molecule),    target, intent(in)    :: molecule
    type(s_grotop),              intent(in)    :: grotop
    type(s_domain),      target, intent(in)    :: domain
    type(s_constraints), target, intent(in)    :: constraints
    type(s_enefunc),     target, intent(inout) :: enefunc

    ! local variables
    integer                  :: dupl, ioffset_dupl, pbc_int
    integer                  :: i, j, k, ia, ib, ic
    integer                  :: nwat, ioffset, icel_local, nangl_a
    integer                  :: i1, i2, i3, idx1, idx2, idx3, icel1, icel2
    integer                  :: nangl_sys
    character(6)             :: res1, res2, res3
    real(wp)                 :: cwork(3,3), dij(3)

    type(s_grotop_mol), pointer :: gromol
    real(wp),           pointer :: force(:,:), theta(:,:)
    real(wp),           pointer :: coord(:,:), box_size(:)
    integer,            pointer :: angle(:), list(:,:,:)
    integer,            pointer :: ncel
    integer(int2),      pointer :: cell_pair(:,:)
    integer(int2),      pointer :: id_g2l(:,:)
    integer,            pointer :: angl_pbc(:,:,:)


    coord     => molecule%atom_coord

    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l
    box_size  => domain%system_size

    angle     => enefunc%num_angle
    list      => enefunc%angle_list
    force     => enefunc%angle_force_const
    theta     => enefunc%angle_theta_min
    angl_pbc  => enefunc%angle_pbc
    nwat      =  enefunc%table%num_water

    nangl_sys = 0
    nangl_a   = 0

    do dupl = 1, domain%num_duplicate

      ioffset_dupl = (dupl-1) * enefunc%table%num_all
      ioffset = 0

      do i = 1, grotop%num_molss
        gromol => grotop%molss(i)%moltype%mol
        do j = 1, grotop%molss(i)%count

          if (gromol%settles%func == 0) then

            do k = 1, gromol%num_angls

              i1 = gromol%angls(k)%atom_idx1
              i2 = gromol%angls(k)%atom_idx2
              i3 = gromol%angls(k)%atom_idx3

              res1 = gromol%atoms(i1)%residue_name
              res2 = gromol%atoms(i2)%residue_name
              res3 = gromol%atoms(i3)%residue_name

              if (res1(1:4) .ne. constraints%water_model .and. &
                  res2(1:4) .ne. constraints%water_model .and. &
                  res3(1:4) .ne. constraints%water_model) then

                idx1 = i1 + ioffset
                idx2 = i2 + ioffset
                idx3 = i3 + ioffset
                ia   = idx1 + ioffset_dupl
                ib   = idx2 + ioffset_dupl
                ic   = idx3 + ioffset_dupl

                icel1 = id_g2l(1,ia)
                icel2 = id_g2l(1,ic)

                if (icel1 /= 0 .and. icel2 /= 0) then

                  icel_local = cell_pair(icel1,icel2)

                  if (icel_local > 0 .and. icel_local <= ncel) then

                    if (angle(icel_local)+1 > MaxAngle) &
                      call error_msg( &
                              'Setup_Enefunc_Angl_Constraint> Too many angles.')

                    nangl_a = nangl_a + 1

                    angle(icel_local) = angle(icel_local) + 1
                    list (1:3,angle(icel_local),icel_local) = (/ia,ib,ic/)
                    force(    angle(icel_local),icel_local) = &
                                      gromol%angls(k)%kt * JOU2CAL * 0.5_wp
                    theta(    angle(icel_local),icel_local) = &
                                      gromol%angls(k)%theta_0 * RAD
                    cwork(1:3,1) = coord(1:3,idx1)
                    cwork(1:3,2) = coord(1:3,idx2)
                    cwork(1:3,3) = coord(1:3,idx3)
                    dij(1:3) = cwork(1:3,1) - cwork(1:3,2)
                    call check_pbc(box_size, dij, pbc_int)
                    angl_pbc(1,angle(icel_local),icel_local) = pbc_int
                    dij(1:3) = cwork(1:3,3) - cwork(1:3,2)
                    call check_pbc(box_size, dij, pbc_int)
                    angl_pbc(2,angle(icel_local),icel_local) = pbc_int
                    dij(1:3) = cwork(1:3,1) - cwork(1:3,3)
                    call check_pbc(box_size, dij, pbc_int)
                    angl_pbc(3,angle(icel_local),icel_local) = pbc_int

                  end if

                end if

              end if

            end do

            nangl_sys = nangl_sys + gromol%num_angls

          else

            nangl_sys = nangl_sys + 1

          end if

          ioffset   = ioffset   + gromol%num_atoms
        
        end do
      end do

    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(nangl_a, enefunc%num_angl_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_angl_all = nangl_a
#endif

    if (enefunc%num_angl_all /= (nangl_sys - nwat*domain%num_duplicate)) &
      call error_msg( &
        'Setup_Enefunc_Angl_Constraint> Some angle paremeters are missing.')

    return

  end subroutine setup_enefunc_angl_constraint

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_dihe
  !> @brief        define DIHEDRAL term in potential energy function
  !! @authors      NT
  !! @param[in]    molecule : molecule information
  !! @param[in]    grotop   : GROMACS TOP information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_dihe(molecule, grotop, domain, enefunc)

    ! formal arguments
    type(s_molecule),target, intent(in)    :: molecule
    type(s_grotop),  target, intent(in)    :: grotop
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: dupl, ioffset_dupl, pbc_int
    integer                  :: i, j, k, ia, ib, ic, id
    integer                  :: ioffset, found
    integer                  :: idx1, idx2, idx3, idx4, icel1, icel2, icel_local
    real(wp)                 :: cwork(3,4), dij(3)

    type(s_grotop_mol), pointer :: gromol
    real(wp),           pointer :: force(:,:), phase(:,:)
    real(wp),           pointer :: coord(:,:), box_size(:)
    integer,            pointer :: dihe(:), list(:,:,:), period(:,:)
    integer,            pointer :: ncel
    integer(int2),      pointer :: cell_pair(:,:)
    integer(int2),      pointer :: id_g2l(:,:)
    integer,            pointer :: ndihe
    integer,            pointer :: notation
    integer,            pointer :: dihe_pbc(:,:,:)


    coord     => molecule%atom_coord

    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l
    box_size  => domain%system_size

    dihe      => enefunc%num_dihedral
    list      => enefunc%dihe_list
    force     => enefunc%dihe_force_const
    phase     => enefunc%dihe_phase
    period    => enefunc%dihe_periodicity
    notation  => enefunc%notation_14types
    dihe_pbc  => enefunc%dihe_pbc
    notation  = 100

    do dupl = 1, domain%num_duplicate

      ioffset_dupl = (dupl-1)*enefunc%table%num_all
      ioffset = 0

      do i = 1, grotop%num_molss
        gromol => grotop%molss(i)%moltype%mol
        do j = 1, grotop%molss(i)%count

          do k = 1, gromol%num_dihes

            if (gromol%dihes(k)%func /= 1 .and. gromol%dihes(k)%func /= 4 &
                .and. gromol%dihes(k)%func /= 9) &
              cycle

            idx1 = gromol%dihes(k)%atom_idx1 + ioffset
            idx2 = gromol%dihes(k)%atom_idx2 + ioffset
            idx3 = gromol%dihes(k)%atom_idx3 + ioffset
            idx4 = gromol%dihes(k)%atom_idx4 + ioffset
            ia   = idx1 + ioffset_dupl
            ib   = idx2 + ioffset_dupl
            ic   = idx3 + ioffset_dupl
            id   = idx4 + ioffset_dupl

            icel1 = id_g2l(1,ia)
            icel2 = id_g2l(1,id)

            if (icel1 /= 0 .and. icel2 /= 0) then

              icel_local = cell_pair(icel1,icel2)

              if (icel_local > 0 .and. icel_local <= ncel) then

                ndihe => dihe(icel_local)
                ndihe = ndihe + 1

                if (ndihe > MaxDihe) &
                  call error_msg('Setup_Enefunc_Dihe> Too many dihedrals.') 

                list (1:4,ndihe,icel_local) = (/ia, ib, ic, id/)
                force (   ndihe,icel_local) = gromol%dihes(k)%kp * JOU2CAL
                phase (   ndihe,icel_local) = gromol%dihes(k)%ps * RAD
                period(   ndihe,icel_local) = gromol%dihes(k)%multiplicity
                cwork(1:3,1) = coord(1:3,idx1)
                cwork(1:3,2) = coord(1:3,idx2)
                cwork(1:3,3) = coord(1:3,idx3)
                cwork(1:3,4) = coord(1:3,idx4)
                dij(1:3) = cwork(1:3,1) - cwork(1:3,2)
                call check_pbc(box_size, dij, pbc_int)
                dihe_pbc(1,ndihe,icel_local) = pbc_int
                dij(1:3) = cwork(1:3,2) - cwork(1:3,3)
                call check_pbc(box_size, dij, pbc_int)
                dihe_pbc(2,ndihe,icel_local) = pbc_int
                dij(1:3) = cwork(1:3,4) - cwork(1:3,3)
                call check_pbc(box_size, dij, pbc_int)
                dihe_pbc(3,ndihe,icel_local) = pbc_int

                if (period(ndihe,icel_local) >  enefunc%notation_14types) &
              call error_msg('Setup_Enefunc_Dihe> Too many periodicity.')

              end if

            end if

          end do

          ioffset = ioffset + gromol%num_atoms

        end do
      end do

    end do

    found = 0 
    do i = 1, ncel
      found = found + dihe(i)
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_dihe_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_dihe_all = found
#endif

    return

  end subroutine setup_enefunc_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_rb_dihe
  !> @brief        define Ryckaert-Bellemans DIHEDRAL term in potential energy
  !                function
  !! @authors      NT
  !! @param[in]    molecule : molecule information
  !! @param[in]    grotop   : GROMACS TOP information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_rb_dihe(molecule, grotop, domain, enefunc)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_grotop),   target, intent(in)    :: grotop
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variables
    integer                  :: dupl, ioffset_dupl, pbc_int
    integer                  :: i, j, k
    integer                  :: ioffset, found
    integer                  :: idx1, idx2, idx3, idx4, icel1, icel2, icel_local
    integer                  :: ia, ib, ic, id
    real(wp)                 :: cwork(3,4), dij(3)

    type(s_grotop_mol), pointer :: gromol
    real(wp),           pointer :: c(:,:,:), box_size(:), coord(:,:)
    integer,            pointer :: dihe(:), list(:,:,:)
    integer,            pointer :: ncel
    integer(int2),      pointer :: cell_pair(:,:)
    integer(int2),      pointer :: id_g2l(:,:)
    integer,            pointer :: ndihe
    integer,            pointer :: dihe_pbc(:,:,:)


    coord     => molecule%atom_coord

    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l
    box_size  => domain%system_size

    dihe      => enefunc%num_rb_dihedral
    list      => enefunc%rb_dihe_list
    c         => enefunc%rb_dihe_c
    dihe_pbc  => enefunc%rb_dihe_pbc

    do dupl = 1, domain%num_duplicate

      ioffset_dupl = (dupl-1)*enefunc%table%num_all
      ioffset = 0

      do i = 1, grotop%num_molss
        gromol => grotop%molss(i)%moltype%mol
        do j = 1, grotop%molss(i)%count

          do k = 1, gromol%num_dihes

            if (gromol%dihes(k)%func /= 3) &
              cycle

            idx1 = gromol%dihes(k)%atom_idx1 + ioffset
            idx2 = gromol%dihes(k)%atom_idx2 + ioffset
            idx3 = gromol%dihes(k)%atom_idx3 + ioffset
            idx4 = gromol%dihes(k)%atom_idx4 + ioffset
            ia   = idx1 + ioffset_dupl
            ib   = idx2 + ioffset_dupl
            ic   = idx3 + ioffset_dupl
            id   = idx4 + ioffset_dupl

            icel1 = id_g2l(1,ia)
            icel2 = id_g2l(1,id)

            if (icel1 /= 0 .and. icel2 /= 0) then

              icel_local = cell_pair(icel1,icel2)

              if (icel_local > 0 .and. icel_local <= ncel) then

                ndihe => dihe(icel_local)
                ndihe = ndihe + 1

                if (ndihe > MaxDihe) &
                call error_msg('Setup_Enefunc_RB_Dihe> Too many dihedrals.') 
  
                list (1:4,ndihe,icel_local) = (/ia, ib, ic, id/)
                c    (1:6,ndihe,icel_local) = gromol%dihes(k)%c(1:6) * JOU2CAL
                cwork(1:3,1) = coord(1:3,idx1)
                cwork(1:3,2) = coord(1:3,idx2)
                cwork(1:3,3) = coord(1:3,idx3)
                cwork(1:3,4) = coord(1:3,idx4)
                dij(1:3) = cwork(1:3,1) - cwork(1:3,2)
                call check_pbc(box_size, dij, pbc_int)
                dihe_pbc(1,ndihe,icel_local) = pbc_int
                dij(1:3) = cwork(1:3,2) - cwork(1:3,3)
                call check_pbc(box_size, dij, pbc_int)
                dihe_pbc(2,ndihe,icel_local) = pbc_int
                dij(1:3) = cwork(1:3,4) - cwork(1:3,3)
                call check_pbc(box_size, dij, pbc_int)
                dihe_pbc(3,ndihe,icel_local) = pbc_int
              end if

            end if

          end do

          ioffset = ioffset + gromol%num_atoms

        end do
      end do

    end do

    found = 0 
    do i = 1, ncel
      found = found + dihe(i)
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_rb_dihe_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_rb_dihe_all = found
#endif

    return

  end subroutine setup_enefunc_rb_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_impr
  !> @brief        define improper DIHEDRAL term in potential energy function
  !! @authors      NT
  !! @param[in]    molecule : molecule information
  !! @param[in]    grotop   : GROMACS TOP information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_impr(molecule, grotop, domain, enefunc)

    ! formal arguments
    type(s_molecule),target, intent(in)    :: molecule
    type(s_grotop),  target, intent(in)    :: grotop
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables 
    integer                  :: dupl, ioffset_dupl
    integer                  :: i, j, k, pbc_int
    integer                  :: ioffset, found
    integer                  :: idx1, idx2, idx3, idx4, icel1, icel2, icel_local
    integer                  :: ia, ib, ic, id
    real(wp)                 :: cwork(3,4), dij(3)

    type(s_grotop_mol), pointer :: gromol
    real(wp),           pointer :: force(:,:), phase(:,:)
    real(wp),           pointer :: coord(:,:), box_size(:)
    integer,            pointer :: impr(:), list(:,:,:), period(:,:)
    integer,            pointer :: ncel
    integer(int2),      pointer :: cell_pair(:,:)
    integer(int2),      pointer :: id_g2l(:,:)
    integer,            pointer :: nimpr
    integer,            pointer :: impr_pbc(:,:,:)


    coord     => molecule%atom_coord

    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l
    box_size  => domain%system_size

    impr      => enefunc%num_improper
    list      => enefunc%impr_list
    force     => enefunc%impr_force_const
    phase     => enefunc%impr_phase
    period    => enefunc%impr_periodicity
    impr_pbc  => enefunc%impr_pbc

    do dupl = 1, domain%num_duplicate

      ioffset_dupl = (dupl-1)*enefunc%table%num_all
      ioffset = 0

      do i = 1, grotop%num_molss
        gromol => grotop%molss(i)%moltype%mol
        do j = 1, grotop%molss(i)%count

          do k = 1, gromol%num_dihes

            if (gromol%dihes(k)%func /= 2) &
              cycle

            idx1 = gromol%dihes(k)%atom_idx1 + ioffset 
            idx2 = gromol%dihes(k)%atom_idx2 + ioffset 
            idx3 = gromol%dihes(k)%atom_idx3 + ioffset 
            idx4 = gromol%dihes(k)%atom_idx4 + ioffset 
            ia   = idx1 + ioffset_dupl
            ib   = idx2 + ioffset_dupl
            ic   = idx3 + ioffset_dupl
            id   = idx4 + ioffset_dupl

            icel1 = id_g2l(1,ia)
            icel2 = id_g2l(1,id)

            if (icel1 /= 0 .and. icel2 /= 0) then

              icel_local = cell_pair(icel1,icel2)

              if (icel_local > 0 .and. icel_local <= ncel) then

                nimpr => impr(icel_local)
                nimpr = nimpr + 1

                if (nimpr > MaxImpr) &
                  call error_msg('Setup_Enefunc_Impr> Too many impropers.') 

                list (1:4,nimpr,icel_local) = (/ia, ib, ic, id/)
                force (   nimpr,icel_local) = gromol%dihes(k)%kp & 
                                            * JOU2CAL * 0.5_wp
                phase (   nimpr,icel_local) = gromol%dihes(k)%ps * RAD
                period(   nimpr,icel_local) = gromol%dihes(k)%multiplicity
                cwork(1:3,1) = coord(1:3,idx1)
                cwork(1:3,2) = coord(1:3,idx2)
                cwork(1:3,3) = coord(1:3,idx3)
                cwork(1:3,4) = coord(1:3,idx4)
                dij(1:3) = cwork(1:3,1) - cwork(1:3,2)
                call check_pbc(box_size, dij, pbc_int)
                impr_pbc(1,nimpr,icel_local) = pbc_int
                dij(1:3) = cwork(1:3,2) - cwork(1:3,3)
                call check_pbc(box_size, dij, pbc_int)
                impr_pbc(2,nimpr,icel_local) = pbc_int
                dij(1:3) = cwork(1:3,4) - cwork(1:3,3)
                call check_pbc(box_size, dij, pbc_int)
                impr_pbc(3,nimpr,icel_local) = pbc_int

              end if

            end if

          end do

          ioffset = ioffset + gromol%num_atoms

        end do
      end do

    end do

    found = 0 
    do i = 1, ncel
      found = found + impr(i)
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_impr_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_impr_all = found
#endif

    return

  end subroutine setup_enefunc_impr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_nonb
  !> @brief        define NON-BOND term in potential energy function
  !! @authors      NT
  !! @param[in]    molecule    : molecule information
  !! @param[in]    grotop      : GROMACS information
  !! @param[in]    constraints : constraints information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_nonb(grotop, molecule, constraints, domain, enefunc)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_constraints),     intent(in)    :: constraints
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: eps1, sig, ei, ej, si, sj
    real(wp)                 :: c6i, c6j, c12i, c12j, c6, c12
    real(wp)                 :: vi, vj, wi, wj, vij, wij
    integer                  :: nnonb, ncel, i, j, k, excl_level
    integer                  :: ix, jx
    integer                  :: cls_local

    integer,    allocatable  :: check_cls(:)
    integer,    allocatable  :: atmcls_map_g2l(:), atmcls_map_l2g(:)
    real(wp),   allocatable  :: nb14_lj6(:,:), nb14_lj12(:,:)
    real(wp),   allocatable  :: nonb_lj6(:,:), nonb_lj12(:,:)


    enefunc%num_atom_cls = grotop%num_atomtypes
    enefunc%fudge_lj     = grotop%defaults%fudge_lj
    enefunc%fudge_qq     = grotop%defaults%fudge_qq

    ncel                 = domain%num_cell_local
    ELECOEF              = ELECOEF_GROMACS

    ! set lennard-jones parameters
    !
    nnonb = enefunc%num_atom_cls

    allocate(check_cls(nnonb),        &
             atmcls_map_g2l(nnonb),   &
             atmcls_map_l2g(nnonb),   &
             nb14_lj6 (nnonb, nnonb), &
             nb14_lj12(nnonb, nnonb), &
             nonb_lj6 (nnonb, nnonb), &
             nonb_lj12(nnonb, nnonb))

    check_cls(1:nnonb)          = 0
    nb14_lj6 (1:nnonb, 1:nnonb) = 0.0_wp
    nb14_lj12(1:nnonb, 1:nnonb) = 0.0_wp
    nonb_lj6 (1:nnonb, 1:nnonb) = 0.0_wp
    nonb_lj12(1:nnonb, 1:nnonb) = 0.0_wp

    do i = 1, nnonb

!     lj_coef(1,i) = grotop%atomtypes(i)%v * 10.0_wp
!     lj_coef(2,i) = grotop%atomtypes(i)%w * JOU2CAL

      do j = 1, nnonb

        vi = grotop%atomtypes(i)%v
        vj = grotop%atomtypes(j)%v
        wi = grotop%atomtypes(i)%w
        wj = grotop%atomtypes(j)%w

        if (grotop%defaults%combi_rule == 2) then

          si = vi * 10.0_wp
          sj = vj * 10.0_wp

          ei = wi * JOU2CAL
          ej = wj * JOU2CAL

          sig = (si + sj) * 0.5_wp
          eps1 = sqrt(ei * ej)

          c6  = 4.0_wp * eps1 * (sig ** 6)
          c12 = 4.0_wp * eps1 * (sig ** 12)

        else ! combi_rule == 1 or 3

          c6i  = vi * 1000000.0_wp * JOU2CAL
          c6j  = vj * 1000000.0_wp * JOU2CAL

          c12i = wi * 1000000.0_wp * 1000000.0_wp * JOU2CAL
          c12j = wj * 1000000.0_wp * 1000000.0_wp * JOU2CAL

          c6  = sqrt(c6i  * c6j)
          c12 = sqrt(c12i * c12j)

        end if

        if (grotop%num_nbonparms > 0) then

          vij = 0.0_wp
          wij = 0.0_wp

          do k = 1, grotop%num_nbonparms
            if (grotop%atomtypes(i)%type_name .eq. &
                  grotop%nbonparms(k)%atom_type1 .and. &
                grotop%atomtypes(j)%type_name .eq. &
                  grotop%nbonparms(k)%atom_type2 .or.  &
                grotop%atomtypes(j)%type_name .eq. &
                  grotop%nbonparms(k)%atom_type1 .and. &
                grotop%atomtypes(i)%type_name .eq. &
                  grotop%nbonparms(k)%atom_type2) then

              vij = grotop%nbonparms(k)%v
              wij = grotop%nbonparms(k)%w

              exit
            end if
          end do

          if (abs(vij) > eps .and. abs(wij) > eps) then

            if (grotop%defaults%combi_rule == 2) then

              sig = vij * 10.0_wp
              eps1 = wij * JOU2CAL
  
              c6  = 4.0_wp * eps1 * (sig ** 6)
              c12 = 4.0_wp * eps1 * (sig ** 12)

            else ! combi_rule = 1 or 3

              c6  = vij * 1000000.0_wp * JOU2CAL
              c12 = wij * 1000000.0_wp * 1000000.0_wp * JOU2CAL

            end if

          end if

        end if

        ! set parameters
        !
        nb14_lj12(i,j) = c12
        nb14_lj6 (i,j) = c6

        nonb_lj12(i,j) = c12
        nonb_lj6 (i,j) = c6

      end do
    end do

    ! check # of exclusion level
    !
    if (enefunc%forcefield == ForcefieldGROMARTINI) then

      !TODO

      enefunc%excl_level = -1

      do i = 1, grotop%num_molss
        excl_level = grotop%molss(i)%moltype%exclude_nbon
        if (enefunc%excl_level == -1) then
          enefunc%excl_level = excl_level
        else if (enefunc%excl_level /= excl_level) then
          call error_msg( &
               'Setup_Enefunc_Nonb> multiple "exclude_nbon" is not supported.')
        end if

      end do

    end if

    ! check the usage of atom class
    !
    do i = 1, molecule%num_atoms
      k = molecule%atom_cls_no(i)
      if (k < 1) then
        call error_msg( &
        'Setup_Enefunc_Nonb> atom class is not defined: "'&
        //trim(molecule%atom_cls_name(i))//'"')
      end if
      check_cls(k) = 1
    end do

    k = 0
    do i = 1, nnonb
      if (check_cls(i) == 1) then
        k = k + 1
        atmcls_map_g2l(i) = k
        atmcls_map_l2g(k) = i
      end if
    end do
    cls_local = k
    max_class = cls_local

    call alloc_enefunc(enefunc, EneFuncNbon, cls_local)

    do i = 1, cls_local
      ix = atmcls_map_l2g(i)
      do j = 1, cls_local
        jx = atmcls_map_l2g(j)
        enefunc%nb14_lj12(i,j) = nb14_lj12(ix,jx)
        enefunc%nb14_lj6 (i,j) = nb14_lj6 (ix,jx)
        enefunc%nonb_lj12(i,j) = nonb_lj12(ix,jx)
        enefunc%nonb_lj6 (i,j) = nonb_lj6 (ix,jx)
      end do
    end do

    ! update domain information
    !
    do i = 1, domain%num_cell_local+domain%num_cell_boundary
      do ix = 1, domain%num_atom(i)
        domain%atom_cls_no(ix,i) = atmcls_map_g2l(domain%atom_cls_no(ix,i))
      end do
    end do

    if (enefunc%table%num_water > 0) then
      domain%water%atom_cls_no(1:3)  &
        = atmcls_map_g2l(domain%water%atom_cls_no(1:3))
      enefunc%table%atom_cls_no_O = atmcls_map_g2l(enefunc%table%atom_cls_no_O)
      enefunc%table%atom_cls_no_H = atmcls_map_g2l(enefunc%table%atom_cls_no_H)
      if (constraints%water_type == TIP4) then
        enefunc%table%atom_cls_no_D =atmcls_map_g2l(enefunc%table%atom_cls_no_D)
        domain%water%atom_cls_no(4) =atmcls_map_g2l(domain%water%atom_cls_no(4))
      end if
    end if
 
    deallocate(check_cls,      &
               atmcls_map_g2l, &
               atmcls_map_l2g, &
               nb14_lj6,       &
               nb14_lj12,      &
               nonb_lj6,       &
               nonb_lj12)

    enefunc%num_atom_cls = cls_local

    ! treatment for 1-2, 1-3, 1-4 interactions
    !
    ncel   = domain%num_cell_local 

    call alloc_enefunc(enefunc, EneFuncNonb,     ncel, maxcell_near)
    call alloc_enefunc(enefunc, EneFuncNonbList, ncel, maxcell_near)

    if (constraints%rigid_bond) then

      call count_nonb_excl(.true., .true., constraints, domain, enefunc)

    else

      call count_nonb_excl(.true., .false., constraints, domain, enefunc)

    end if

    return

  end subroutine setup_enefunc_nonb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_gro_restraints
  !> @brief        setup restraints from GROTOP information
  !! @authors      NT
  !! @param[in]    molecule : molecule information
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_gro_restraints(molecule, grotop, domain, enefunc)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_grotop),          intent(in)    :: grotop
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    type(s_enefunc)          :: ef0
    real(wp)                 :: kx, ky, kz
    integer                  :: nposres, npr_atom, max_pr_atom, natom, ioffset
    integer                  :: i, j, k, n, n2, n3, group0, func0
    integer                  :: ix, iatm, icel

    type(s_grotop_mol), pointer :: gromol


    ! count # of position restraints
    !
    nposres     = 0
    max_pr_atom = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol

      if (gromol%num_posress > 0) &
        nposres = nposres + 1

      npr_atom = grotop%molss(i)%count * gromol%num_posress
      max_pr_atom = max(max_pr_atom,npr_atom)
    end do

    if (nposres == 0) &
      return

    enefunc%restraint = .true.


    ! setup EneFuncRefc
    !
    if (.not. allocated(enefunc%restraint_refcoord)) then
      n = molecule%num_atoms
      call alloc_enefunc(enefunc, EneFuncRefc, n)
      enefunc%restraint_refcoord(1:3,1:n) = molecule%atom_coord(1:3,1:n)
    end if


    ! setup EneFuncRefg
    !
    if (allocated(enefunc%restraint_numatoms)) then

      n  = size(enefunc%restraint_numatoms)
      n2 = size(enefunc%restraint_atomlist(:,1))
      call alloc_enefunc(ef0, EneFuncRefg, n, n2)
      ef0%restraint_numatoms     (1:n) = enefunc%restraint_numatoms     (1:n)
      ef0%restraint_atomlist(1:n2,1:n) = enefunc%restraint_atomlist(1:n2,1:n)
      ef0%restraint_masscoef(1:n2,1:n) = enefunc%restraint_masscoef(1:n2,1:n)

      call alloc_enefunc(enefunc, EneFuncRefg, n+nposres, max(n2,max_pr_atom))
      enefunc%restraint_numatoms     (1:n) = ef0%restraint_numatoms     (1:n)
      enefunc%restraint_atomlist(1:n2,1:n) = ef0%restraint_atomlist(1:n2,1:n)
      enefunc%restraint_masscoef(1:n2,1:n) = ef0%restraint_masscoef(1:n2,1:n)

      call dealloc_enefunc(ef0, EneFuncRefg)

    else

      n = 0
      call alloc_enefunc(enefunc, EneFuncRefg, nposres, max_pr_atom)

    end if


    group0 = n

    natom   = 0
    nposres = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol

      if (gromol%num_posress > 0) &
        nposres = nposres + 1

      npr_atom = 0
      do j = 1, grotop%molss(i)%count
        ioffset = natom
        natom   = natom + gromol%num_atoms

        do k = 1, gromol%num_posress
          npr_atom = npr_atom + 1
          enefunc%restraint_atomlist(npr_atom, group0+nposres) = &
               gromol%posress(k)%atom_idx + ioffset

          if (k > 1 .and. main_rank) then
            if (gromol%posress(k-1)%kx /= gromol%posress(k)%kx .or. &
                gromol%posress(k-1)%ky /= gromol%posress(k)%ky .or. &
                gromol%posress(k-1)%kz /= gromol%posress(k)%kz) then
              write(MsgOut,'(a)') 'Setup_Enefunc_Gro_Restraints> WARNING:'
              write(MsgOut,'(a,a,a)') &
  '   different restraint constant between foreach atoms is not supported. [', &
  trim(grotop%molss(i)%moltype%name), ']'
              write(MsgOut,'(a)') ' '
            end if
          end if

        end do

      end do

      if (npr_atom > 0) &
        enefunc%restraint_numatoms(group0+nposres) = npr_atom

    end do


    ! setup EneFuncReff
    !
    if (allocated(enefunc%restraint_kind)) then

      n  = size(enefunc%restraint_kind)
      n2 = size(enefunc%restraint_grouplist(:,1))
      n3 = max(int(n2/2),1)
      call alloc_enefunc(ef0, EneFuncReff, n, n2)
      ef0%restraint_kind          (1:n) = enefunc%restraint_kind          (1:n)
      ef0%restraint_grouplist(1:n2,1:n) = enefunc%restraint_grouplist(1:n2,1:n)
      ef0%restraint_const    (1:4, 1:n) = enefunc%restraint_const    (1:4, 1:n)
      ef0%restraint_ref      (1:2, 1:n) = enefunc%restraint_ref      (1:2, 1:n)
      ef0%restraint_funcgrp       (1:n) = enefunc%restraint_funcgrp       (1:n)
      ef0%restraint_exponent_func (1:n) = enefunc%restraint_exponent_func (1:n)
      ef0%restraint_exponent_dist (1:n3,1:n) &
                                   = enefunc%restraint_exponent_dist (1:n3,1:n)
      ef0%restraint_weight_dist   (1:n3,1:n) &
                                   = enefunc%restraint_weight_dist   (1:n3,1:n)

      call alloc_enefunc(enefunc, EneFuncReff, n+nposres, max(n2,1))
      enefunc%restraint_kind          (1:n) = ef0%restraint_kind          (1:n)
      enefunc%restraint_grouplist(1:n2,1:n) = ef0%restraint_grouplist(1:n2,1:n)
      enefunc%restraint_const    (1:4, 1:n) = ef0%restraint_const    (1:4, 1:n)
      enefunc%restraint_ref      (1:2, 1:n) = ef0%restraint_ref      (1:2, 1:n)
      enefunc%restraint_funcgrp       (1:n) = ef0%restraint_funcgrp       (1:n)
      enefunc%restraint_exponent_func (1:n) = ef0%restraint_exponent_func (1:n)
      enefunc%restraint_exponent_dist (1:n3,1:n) &
                                       = ef0%restraint_exponent_dist (1:n3,1:n)
      enefunc%restraint_weight_dist   (1:n3,1:n) &
                                       = ef0%restraint_weight_dist   (1:n3,1:n)

      call dealloc_enefunc(ef0, EneFuncReff)

    else

      n = 0
      call alloc_enefunc(enefunc, EneFuncReff, nposres, 1)

    end if


    func0 = n

    nposres = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol

      if (gromol%num_posress == 0) &
        cycle

      nposres = nposres + 1

      kx = gromol%posress(1)%kx * 0.01_wp * JOU2CAL * 0.5_wp
      ky = gromol%posress(1)%ky * 0.01_wp * JOU2CAL * 0.5_wp
      kz = gromol%posress(1)%kz * 0.01_wp * JOU2CAL * 0.5_wp

      enefunc%restraint_kind         (func0+nposres) = RestraintsFuncPOSI
      enefunc%restraint_funcgrp      (func0+nposres) = 1
      enefunc%restraint_grouplist  (1,func0+nposres) = group0 + nposres
      enefunc%restraint_const      (1,func0+nposres) = 1.0_wp
      enefunc%restraint_const      (2,func0+nposres) = kx
      enefunc%restraint_const      (3,func0+nposres) = ky
      enefunc%restraint_const      (4,func0+nposres) = kz
      enefunc%restraint_exponent_func(func0+nposres) = 2

    end do

    enefunc%num_restraintgroups = group0+nposres
    enefunc%num_restraintfuncs  = func0+nposres


    ! setup domain restraints
    !
    call alloc_enefunc(enefunc, EneFuncRest, domain%num_cell_local, 0)

    do i = 1, nposres

      do ix = 1, enefunc%restraint_numatoms( &
                    enefunc%restraint_grouplist(1,func0+i))

        iatm = enefunc%restraint_atomlist( &
                ix, enefunc%restraint_grouplist(1,func0+i))

        icel = domain%id_g2l(1,iatm)

        if (icel > 0 .and. icel <= domain%num_cell_local) then
          enefunc%num_restraint(icel) = enefunc%num_restraint(icel) + 1

          n = enefunc%num_restraint(icel)
          enefunc%restraint_atom (    n,icel) = iatm
          enefunc%restraint_force(1:4,n,icel) = enefunc%restraint_const(1:4,i)
          enefunc%restraint_coord(1:3,n,icel) = &
                                          enefunc%restraint_refcoord(1:3,iatm)
        end if

      end do

    end do


    ! summary of setup enefunc_gro_restraints
    !
    if (main_rank) then

      write(MsgOut,'(A)')&
           'Setup_Enefunc_Gro_Restraints> Setup restraint functions'

      do i = func0+1, enefunc%num_restraintfuncs

        write(MsgOut,'(A,I5,A,I5)') &
          ' func  = ', i, ' kind  = ', enefunc%restraint_kind(i)

        ! summary for positional restraint
        !
        write(MsgOut,'(A,4F8.3)') &
             ' const(total, x, y, z) = ', enefunc%restraint_const(1:4,i) 
        write(MsgOut,'(A,I5)') ' exponend of function = ', &
                                          enefunc%restraint_exponent_func(i)

        write(MsgOut,'(A,I5)') ' # of groups  = ', enefunc%restraint_funcgrp(i)
        write(MsgOut,'(" grouplist: ",$)') 
        do j = 1, enefunc%restraint_funcgrp(i) 
          write(MsgOut,'(i3,$)') enefunc%restraint_grouplist(j,i)
          if (mod(j,20) == 0 .and. j /= enefunc%restraint_funcgrp(i) )  &
            write(MsgOut,'(A)') ''
        end do
        write(MsgOut,'(A)') ''

      end do

      write(MsgOut,*) ''

    end if

    return

  end subroutine setup_enefunc_gro_restraints

end module sp_enefunc_gromacs_mod
