!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_enefunc_amber_mod
!> @brief   define potential energy functions
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_enefunc_amber_mod

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
  use fileio_prmtop_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef MPI
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: define_enefunc_amber
  private :: setup_enefunc_bond
  private :: setup_enefunc_bond_constraint
  private :: setup_enefunc_angl
  private :: setup_enefunc_angl_constraint
  private :: setup_enefunc_dihe
  private :: setup_enefunc_impr
  private :: setup_enefunc_nonb

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc_amber
  !> @brief        a driver subroutine for defining potential energy
  !! @authors      NT
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    prmtop      : AMBER parameter topology information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    constraints : constraints information
  !! @param[in]    restraints  : restraints information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc_amber(ene_info, prmtop, molecule, &
                                  constraints, restraints, domain, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info 
    type(s_prmtop),          intent(in)    :: prmtop
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
      call setup_enefunc_bond(molecule, prmtop, domain, enefunc)

      ! angle
      !
      call setup_enefunc_angl(molecule, prmtop, domain, enefunc)

    else

      ! bond
      !
      call setup_enefunc_bond_constraint(prmtop, molecule, &
                                         domain, constraints, enefunc)

      ! angle
      !
      call setup_enefunc_angl_constraint(prmtop, molecule, &
                                         domain, constraints, enefunc)

    end if

    ! dihedral
    !
    call setup_enefunc_dihe(molecule, prmtop, domain, enefunc)

    ! improper
    !
    call setup_enefunc_impr(molecule, prmtop, domain, enefunc)

    ! nonbonded
    !
    call setup_enefunc_nonb(prmtop, molecule, constraints, domain, enefunc)

    ! lookup table
    !
    call setup_enefunc_table(ene_info, enefunc)

    ! restraints
    !
    call setup_enefunc_restraints(molecule, restraints, domain, enefunc)

    ! write summary of energy function
    !
    if (main_rank) then

      write(MsgOut,'(A)') &
           'Define_Enefunc_Amber> Number of Interactions in Each Term'
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           '  bond_ene        = ', enefunc%num_bond_all,        &
           '  angle_ene       = ', enefunc%num_angl_all
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           '  torsion_ene     = ', enefunc%num_dihe_all,        &
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

  end subroutine define_enefunc_amber

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_bond
  !> @brief        define BOND term in potential energy function
  !! @authors      NT
  !! @param[in]    prmtop   : AMBER parameter topology information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond(molecule, prmtop, domain, enefunc)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: m1, m2
    integer                  :: dupl, ioffset
    integer                  :: i, found, i1, i2, icel1, icel2
    integer                  :: icel_local
    character(6)             :: ri1, ri2

    real(wp),        pointer :: force(:,:), dist(:,:)
    integer,         pointer :: nwater(:), bond(:), list(:,:,:)
    integer,         pointer :: ncel
    integer(int2),   pointer :: cell_pair(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: nbond, sollist(:)


    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l
    nwater    => domain%num_water

    bond      => enefunc%num_bond
    list      => enefunc%bond_list
    force     => enefunc%bond_force_const
    dist      => enefunc%bond_dist_min
    sollist   => enefunc%table%solute_list_inv

    do dupl = 1, domain%num_duplicate

      ioffset = (dupl-1) * enefunc%table%num_solute

      do i = 1, prmtop%num_bondh

        i1 = prmtop%bond_inc_hy(1,i) / 3 + 1
        i2 = prmtop%bond_inc_hy(2,i) / 3 + 1
        ri1 = molecule%residue_name(i1)
        ri2 = molecule%residue_name(i2)

        if (ri1(1:3) /= 'TIP' .and. ri1(1:3) /= 'WAT' .and. &
            ri1(1:3) /= 'SOL' .and. ri2(1:3) /= 'TIP' .and. &
            ri2(1:3) /= 'WAT' .and. ri2(1:3) /= 'SOL') then

          i1 = sollist(i1) + ioffset
          i2 = sollist(i2) + ioffset

          icel1 = id_g2l(1,i1)
          icel2 = id_g2l(1,i2)

          if (icel1 /= 0 .and. icel2 /= 0) then

            icel_local = cell_pair(icel1,icel2)

            if (icel_local > 0 .and. icel_local <= ncel) then

              nbond => bond(icel_local)
              nbond = nbond + 1

              if (nbond > MaxBond) &
                call error_msg('Setup_Enefunc_Bond> Too many bonds.') 

              list (1,nbond,icel_local) = i1
              list (2,nbond,icel_local) = i2
              force(  nbond,icel_local) = &
                             prmtop%bond_fcons_uniq(prmtop%bond_inc_hy(3,i))
              dist (  nbond,icel_local) = &
                             prmtop%bond_equil_uniq(prmtop%bond_inc_hy(3,i))
            end if

          end if

        end if

      end do

      do i = 1, prmtop%num_mbonda

        i1 = prmtop%bond_wo_hy(1,i) / 3 + 1
        i2 = prmtop%bond_wo_hy(2,i) / 3 + 1
 
        i1 = sollist(i1) + ioffset
        i2 = sollist(i2) + ioffset

        icel1 = id_g2l(1,i1)
        icel2 = id_g2l(1,i2)

        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1,icel2)

          if (icel_local > 0 .and. icel_local <= ncel) then

            nbond => bond(icel_local)
            nbond = nbond + 1

            if (nbond > MaxBond) &
              call error_msg('Setup_Enefunc_Bond> Too many bonds.') 

            list (1,nbond,icel_local) = i1
            list (2,nbond,icel_local) = i2
            force(  nbond,icel_local) = &
                             prmtop%bond_fcons_uniq(prmtop%bond_wo_hy(3,i))
            dist (  nbond,icel_local) = &
                             prmtop%bond_equil_uniq(prmtop%bond_wo_hy(3,i))
          end if

        end if

      end do

    end do

    ! water molecules
    !
    do i = 1, prmtop%num_bondh

      i1 = prmtop%bond_inc_hy(1,i) / 3 + 1
      i2 = prmtop%bond_inc_hy(2,i) / 3 + 1
      ri1 = molecule%residue_name(i1)
      ri2 = molecule%residue_name(i2)
      m1 = molecule%mass(i1)
      m2 = molecule%mass(i2)

      if (ri1(1:3) == 'TIP' .or. ri1(1:3) == 'WAT' .or. &
          ri1(1:3) == 'SOL') then

        if (abs(m1-m2) > EPS) then
          enefunc%table%water_bond_calc = .true.
          enefunc%table%water_bond_calc_OH = .true.
          enefunc%table%OH_bond = &
               prmtop%bond_equil_uniq(prmtop%bond_inc_hy(3,i))
          enefunc%table%OH_force = &
               prmtop%bond_fcons_uniq(prmtop%bond_inc_hy(3,i))
        else if (m1 == m2 .and. m1 < LIGHT_ATOM_MASS_LIMIT) then
          enefunc%table%water_bond_calc = .true.
          enefunc%table%water_bond_calc_HH = .true.
          enefunc%table%HH_bond = &
               prmtop%bond_equil_uniq(prmtop%bond_inc_hy(3,i))
          enefunc%table%HH_force = &
               prmtop%bond_fcons_uniq(prmtop%bond_inc_hy(3,i))
        end if
      end if

    end do

    found = 0 
    do i = 1, ncel
      found = found + bond(i)
      if (enefunc%table%water_bond_calc_OH) found = found + 2*nwater(i)
      if (enefunc%table%water_bond_calc_HH) found = found + nwater(i)
    end do

#ifdef MPI
    call mpi_allreduce(found, enefunc%num_bond_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_bond_all = found
#endif

    return

  end subroutine setup_enefunc_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_bond_constraint
  !> @brief        define BOND term in potential energy function
  !! @authors      NT
  !! @param[in]    prmtop      : AMBER parameter topology information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    domain      : domain information
  !! @param[in]    constraints : constraints information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond_constraint(prmtop, molecule, domain, &
                                           constraints, enefunc)

    ! formal arguments
    type(s_prmtop),              intent(in)    :: prmtop
    type(s_molecule),            intent(in)    :: molecule
    type(s_domain),      target, intent(in)    :: domain
    type(s_constraints), target, intent(inout) :: constraints
    type(s_enefunc),     target, intent(inout) :: enefunc

    ! local variables
    integer                      :: dupl, ioffset
    integer                      :: i, j, k, ih, icel_local, connect
    integer                      :: i1, i2, ih1, ih2, icel1, icel2, icel
    integer                      :: nbond_a, nbond_c
    integer                      :: wat_bonds, wat_found
    integer                      :: tmp_mole_no, mole_no
    character(6)                 :: ri1, ri2

    real(wp),            pointer :: force(:,:), dist(:,:)
    real(wip),           pointer :: HGr_bond_dist(:,:,:,:)
    integer,             pointer :: bond(:), list(:,:,:), ncel
    integer,             pointer :: id_l2g(:,:), sollist(:)
    integer(int2),       pointer :: cell_pair(:,:)
    integer(int2),       pointer :: id_g2l(:,:)
    integer,             pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)


    ncel          => domain%num_cell_local
    cell_pair     => domain%cell_pair
    id_g2l        => domain%id_g2l
    id_l2g        => domain%id_l2g

    HGr_local     => constraints%HGr_local
    HGr_bond_list => constraints%HGr_bond_list
    HGr_bond_dist => constraints%HGr_bond_dist

    bond          => enefunc%num_bond
    list          => enefunc%bond_list
    force         => enefunc%bond_force_const
    dist          => enefunc%bond_dist_min
    sollist       => enefunc%table%solute_list_inv

    connect       =  constraints%connect

    nbond_a       =  0
    nbond_c       =  0

    do dupl = 1, domain%num_duplicate

      ioffset = (dupl-1) * enefunc%table%num_solute

      do i = 1, prmtop%num_mbonda

        i1 = prmtop%bond_wo_hy(1,i) / 3 + 1
        i2 = prmtop%bond_wo_hy(2,i) / 3 + 1

        ri1 = molecule%residue_name(i1)
        ri2 = molecule%residue_name(i2)

        i1  = sollist(i1) + ioffset
        i2  = sollist(i2) + ioffset

        if (ri1 /= constraints%water_model .and. &
            ri2 /= constraints%water_model) then

          icel1 = id_g2l(1,i1)
          icel2 = id_g2l(1,i2)

          if (icel1 /= 0 .and. icel2 /= 0) then

            icel_local = cell_pair(icel1,icel2)

            if (icel_local > 0 .and. icel_local <= ncel) then

              nbond_a = nbond_a + 1
              bond(icel_local) = bond(icel_local) + 1

              if (bond(icel_local) > MaxBond) &
                call error_msg('Setup_Enefunc_Bond_Constraint> Too many bonds.') 

              list (1,bond(icel_local),icel_local) = i1
              list (2,bond(icel_local),icel_local) = i2
              force(  bond(icel_local),icel_local) = &
                               prmtop%bond_fcons_uniq(prmtop%bond_wo_hy(3,i))
              dist (  bond(icel_local),icel_local) = &
                               prmtop%bond_equil_uniq(prmtop%bond_wo_hy(3,i))
            end if
          end if

        end if

      end do

      do i = 1, prmtop%num_bondh

        i1 = prmtop%bond_inc_hy(1,i) / 3 + 1
        i2 = prmtop%bond_inc_hy(2,i) / 3 + 1

        ri1 = molecule%residue_name(i1)
        ri2 = molecule%residue_name(i2)

        i1 = sollist(i1) + ioffset
        i2 = sollist(i2) + ioffset

        if (ri1 /= constraints%water_model .and. &
            ri2 /= constraints%water_model) then

          icel1 = id_g2l(1,i1)
          icel2 = id_g2l(1,i2)

          if (icel1 /= 0 .and. icel2 /= 0) then

            icel = cell_pair(icel1,icel2)

            if (icel > 0 .and. icel <= ncel) then

              do j = 1, connect

                do k = 1, HGr_local(j,icel)
                  ih1 = id_l2g(HGr_bond_list(1,k,j,icel),icel)
                  do ih = 1, j
                    ih2 = id_l2g(HGr_bond_list(ih+1,k,j,icel),icel)
  
                    if (ih1 == i1 .and. ih2 == i2 .or. &
                        ih2 == i1 .and. ih1 == i2) then
  
                      nbond_c = nbond_c + 1
                      HGr_bond_dist(ih+1,k,j,icel) = &
                              prmtop%bond_equil_uniq(prmtop%bond_inc_hy(3,i))
  
                    end if
  
                  end do
                end do
              end do
            end if
          end if

        end if

      end do

    end do

    ! for water molecule
    !
    wat_bonds = 0
    wat_found = 0

    if (constraints%fast_water) then

      tmp_mole_no = -1

      do i = 1, prmtop%num_mbonda

        i1 = prmtop%bond_wo_hy(1,i) / 3 + 1
        i2 = prmtop%bond_wo_hy(2,i) / 3 + 1

        ri1 = molecule%residue_name(i1)
        ri2 = molecule%residue_name(i2)

        if (ri1 == constraints%water_model .and. &
            ri2 == constraints%water_model) then

          wat_bonds = wat_bonds+1

        end if

      end do

      do i = 1, prmtop%num_bondh
  
        i1 = prmtop%bond_inc_hy(1,i) / 3 + 1
        i2 = prmtop%bond_inc_hy(2,i) / 3 + 1
  
        ri1 = molecule%residue_name(i1)
        ri2 = molecule%residue_name(i2)

        if (ri1 == constraints%water_model .and. &
            ri2 == constraints%water_model) then

          wat_bonds = wat_bonds+1
          mole_no = molecule%molecule_no(molecule%bond_list(1,i))

          if (mole_no /= tmp_mole_no) then
            wat_found = wat_found +1
            tmp_mole_no = mole_no
          end if

        end if

      end do

      if (wat_found /= enefunc%table%num_water) &
        call error_msg( &
          'Setup_Enefunc_Bond_Constraint> # of water is incorrect')

    end if

#ifdef MPI
    call mpi_allreduce(nbond_a, enefunc%num_bond_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
    call mpi_allreduce(nbond_c, constraints%num_bonds, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_bond_all  = nbond_a
    constraints%num_bonds = nbond_c
#endif

    if (constraints%fast_water) then

      if (enefunc%num_bond_all /= prmtop%num_bondh *domain%num_duplicate &
                                  + prmtop%num_mbonda*domain%num_duplicate &
                                  - constraints%num_bonds &
                                  - wat_bonds*domain%num_duplicate) then
        call error_msg( &
          'Setup_Enefunc_Bond_Constraint> Some bond paremeters are missing.')
      end if

    else

      if (enefunc%num_bond_all /=   prmtop%num_bondh *domain%num_duplicate &
                                    + prmtop%num_mbonda*domain%num_duplicate &
                                    - constraints%num_bonds &
                                    -3*enefunc%table%num_water*domain%num_duplicate) then
        call error_msg( &
          'Setup_Enefunc_Bond_Constraint> Some bond paremeters are missing.')
      end if

    end if

    return

  end subroutine setup_enefunc_bond_constraint

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_angl
  !> @brief        define ANGLE term in potential energy function
  !! @authors      NT
  !! @param[in]    prmtop   : AMBER parameter topology information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_angl(molecule, prmtop, domain, enefunc)

    ! formal arguments
    type(s_molecule),target, intent(in)    :: molecule
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: dupl, ioffset
    integer                  :: i, found, i1, i2, i3, icel1, icel2
    integer                  :: icel_local
    character(6)             :: ri1, ri2, ri3

    real(wp),        pointer :: force(:,:), theta(:,:)
    integer,         pointer :: angle(:), alist(:,:,:)
    integer,         pointer :: ncel
    integer(int2),   pointer :: cell_pair(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: nangl, sollist(:), nwater(:)
    character(6),    pointer :: mol_res_name(:)

    mol_res_name => molecule%residue_name

    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l
    nwater    => domain%num_water

    angle     => enefunc%num_angle
    alist     => enefunc%angle_list
    force     => enefunc%angle_force_const
    theta     => enefunc%angle_theta_min
    sollist   => enefunc%table%solute_list_inv

    do dupl = 1, domain%num_duplicate

      ioffset = (dupl-1)*enefunc%table%num_solute

      do i = 1, prmtop%num_anglh

        i1 = prmtop%angl_inc_hy(1,i) / 3 + 1
        i2 = prmtop%angl_inc_hy(2,i) / 3 + 1
        i3 = prmtop%angl_inc_hy(3,i) / 3 + 1

        ri1 = mol_res_name(i1)
        ri2 = mol_res_name(i2)
        ri3 = mol_res_name(i3)

        if (ri1(1:3) /= 'TIP' .and. ri1(1:3) /= 'WAT' .and. &
            ri1(1:3) /= 'SOL' .and. ri3(1:3) /= 'TIP' .and. &
            ri3(1:3) /= 'WAT' .and. ri3(1:3) /= 'SOL') then

          i1 = sollist(i1) + ioffset
          i2 = sollist(i2) + ioffset
          i3 = sollist(i3) + ioffset

          icel1 = id_g2l(1,i1)
          icel2 = id_g2l(1,i3)

          if (icel1 /= 0 .and. icel2 /= 0) then

            icel_local = cell_pair(icel1,icel2)

            if (icel_local > 0 .and. icel_local <= ncel) then

              nangl => angle(icel_local)
              nangl = nangl + 1

              if (nangl > MaxAngle) &
                call error_msg('Setup_Enefunc_Angl> Too many angles.') 

              alist(1:3,nangl,icel_local) = (/i1, i2, i3/)
              force(    nangl,icel_local) = &
                               prmtop%angl_fcons_uniq(prmtop%angl_inc_hy(4,i))
              theta(    nangl,icel_local) = &
                               prmtop%angl_equil_uniq(prmtop%angl_inc_hy(4,i))
            end if

          end if

        end if

      end do

      do i = 1, prmtop%num_mangla

        i1 = prmtop%angl_wo_hy(1,i) / 3 + 1
        i2 = prmtop%angl_wo_hy(2,i) / 3 + 1
        i3 = prmtop%angl_wo_hy(3,i) / 3 + 1

        i1 = sollist(i1) + ioffset
        i2 = sollist(i2) + ioffset
        i3 = sollist(i3) + ioffset

        icel1 = id_g2l(1,i1)
        icel2 = id_g2l(1,i3)

        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1,icel2)

          if (icel_local > 0 .and. icel_local <= ncel) then

            nangl => angle(icel_local)
            nangl = nangl + 1

            if (nangl > MaxAngle) &
              call error_msg('Setup_Enefunc_Angl> Too many angles.') 

            alist(1:3,nangl,icel_local) = (/i1, i2, i3/)
            force(    nangl,icel_local) = &
                             prmtop%angl_fcons_uniq(prmtop%angl_wo_hy(4,i))
            theta(    nangl,icel_local) = &
                             prmtop%angl_equil_uniq(prmtop%angl_wo_hy(4,i))
          end if

        end if

      end do

    end do

    ! angle energy for water
    !
    do i = 1, prmtop%num_anglh

      i1 = prmtop%angl_inc_hy(1,i) / 3 + 1
      i2 = prmtop%angl_inc_hy(2,i) / 3 + 1
      i3 = prmtop%angl_inc_hy(3,i) / 3 + 1

      ri1 = mol_res_name(i1)
      ri2 = mol_res_name(i2)
      ri3 = mol_res_name(i3)

      if (ri1(1:3) == 'TIP' .or. ri1(1:3) == 'WAT' .or. &
          ri1(1:3) == 'SOL') then

        enefunc%table%HOH_angle = &
            prmtop%angl_equil_uniq(prmtop%angl_inc_hy(4,i))
        enefunc%table%HOH_force = &
            prmtop%angl_fcons_uniq(prmtop%angl_inc_hy(4,i))
        enefunc%table%water_angle_calc = .true.
        exit
      end if
    end do

    found = 0 
    do i = 1, ncel
      found = found + angle(i)
      if (enefunc%table%water_bond_calc_OH) found = found + nwater(i)
    end do

#ifdef MPI
    call mpi_allreduce(found, enefunc%num_angl_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_angl_all = found
#endif

    return

  end subroutine setup_enefunc_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_angl_constraint
  !> @brief        define ANGLE term in potential energy function
  !! @authors      NT
  !! @param[in]    prmtop      : AMBER parameter topology information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    domain      : domain information
  !! @param[in]    constraints : constraints information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_angl_constraint(prmtop, molecule, domain, &
                                           constraints, enefunc)

    ! formal arguments
    type(s_prmtop),              intent(in)    :: prmtop
    type(s_molecule),            intent(in)    :: molecule
    type(s_domain),      target, intent(in)    :: domain
    type(s_constraints), target, intent(in)    :: constraints
    type(s_enefunc),     target, intent(inout) :: enefunc

    ! local variables
    integer                  :: dupl, ioffset
    integer                  :: i, found, i1, i2, i3, icel1, icel2
    integer                  :: icel_local
    character(6)             :: res1, res2, res3

    real(wp),        pointer :: force(:,:), theta(:,:)
    integer,         pointer :: angle(:), alist(:,:,:)
    integer,         pointer :: ncel
    integer(int2),   pointer :: cell_pair(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: nangl, sollist(:)


    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l

    angle     => enefunc%num_angle
    alist     => enefunc%angle_list
    force     => enefunc%angle_force_const
    theta     => enefunc%angle_theta_min
    sollist   => enefunc%table%solute_list_inv

    do dupl = 1, domain%num_duplicate

      ioffset = (dupl-1) * enefunc%table%num_solute

      do i = 1, prmtop%num_anglh

        i1 = prmtop%angl_inc_hy(1,i) / 3 + 1
        i2 = prmtop%angl_inc_hy(2,i) / 3 + 1
        i3 = prmtop%angl_inc_hy(3,i) / 3 + 1

        res1 = molecule%residue_name(i1)
        res2 = molecule%residue_name(i2)
        res3 = molecule%residue_name(i3)

        i1 = sollist(i1) + ioffset
        i2 = sollist(i2) + ioffset
        i3 = sollist(i3) + ioffset

        if (res1 /= constraints%water_model .and. &
            res2 /= constraints%water_model .and. &
            res3 /= constraints%water_model) then

          icel1 = id_g2l(1,i1)
          icel2 = id_g2l(1,i3)
  
          if (icel1 /= 0 .and. icel2 /= 0) then
  
            icel_local = cell_pair(icel1,icel2)
  
            if (icel_local > 0 .and. icel_local <= ncel) then
  
              nangl => angle(icel_local)
              nangl = nangl + 1
  
              if (nangl > MaxAngle) &
                call error_msg('Setup_Enefunc_Angl_Constraint> Too many angles.') 
  
              alist(1:3,nangl,icel_local) = (/i1, i2, i3/)
              force(    nangl,icel_local) = &
                               prmtop%angl_fcons_uniq(prmtop%angl_inc_hy(4,i))
              theta(    nangl,icel_local) = &
                               prmtop%angl_equil_uniq(prmtop%angl_inc_hy(4,i))
            end if
  
          end if

        end if

      end do

      do i = 1, prmtop%num_mangla

        i1 = prmtop%angl_wo_hy(1,i) / 3 + 1
        i2 = prmtop%angl_wo_hy(2,i) / 3 + 1
        i3 = prmtop%angl_wo_hy(3,i) / 3 + 1

        i1 = sollist(i1) + ioffset
        i2 = sollist(i2) + ioffset
        i3 = sollist(i3) + ioffset

        icel1 = id_g2l(1,i1)
        icel2 = id_g2l(1,i3)

        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1,icel2)

          if (icel_local > 0 .and. icel_local <= ncel) then

            nangl => angle(icel_local)
            nangl = nangl + 1

            if (nangl > MaxAngle) &
              call error_msg('Setup_Enefunc_Angl_Constraint> Too many angles.') 

            alist(1:3,nangl,icel_local) = (/i1, i2, i3/)
            force(    nangl,icel_local) = &
                             prmtop%angl_fcons_uniq(prmtop%angl_wo_hy(4,i))
            theta(    nangl,icel_local) = &
                             prmtop%angl_equil_uniq(prmtop%angl_wo_hy(4,i))
          end if

        end if

      end do

    end do

    found = 0 
    do i = 1, ncel
      found = found + angle(i)
    end do

#ifdef MPI
    call mpi_allreduce(found, enefunc%num_angl_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_angl_all = found
#endif

    if (enefunc%num_angl_all /= (prmtop%num_anglh + prmtop%num_mangla)       &
                                 *domain%num_duplicate                 .and. &
        enefunc%num_angl_all /= (prmtop%num_anglh + prmtop%num_mangla - enefunc%table%num_water)  &
                                 *domain%num_duplicate ) &
      call error_msg( &
        'Setup_Enefunc_Angl_Constraint> Some angle paremeters are missing.')

    return

  end subroutine setup_enefunc_angl_constraint

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_dihe
  !> @brief        define DIHEDRAL term in potential energy function
  !! @authors      NT
  !! @param[in]    prmtop   : AMBER parameter topology information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_dihe(molecule, prmtop, domain, enefunc)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: dupl, ioffset
    integer                  :: i, found
    integer                  :: i1, i2, i3, i4, icel1, icel2, icel_local

    real(wp),           pointer :: force(:,:), phase(:,:)
    integer,            pointer :: dihe(:), list(:,:,:), period(:,:)
    integer,            pointer :: ncel
    integer(int2),      pointer :: cell_pair(:,:)
    integer(int2),      pointer :: id_g2l(:,:)
    integer,            pointer :: ndihe
    integer,            pointer :: notation


    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l

    dihe      => enefunc%num_dihedral
    list      => enefunc%dihe_list
    force     => enefunc%dihe_force_const
    phase     => enefunc%dihe_phase
    period    => enefunc%dihe_periodicity
    notation  => enefunc%notation_14types

    notation = 100
    if (prmtop%num_uniqdihe > 100) then
      notation = 1000
      if (prmtop%num_uniqdihe > 1000) then
        call error_msg('Setup_Enefunc_Dihe> Too many dihedrals.') 
      endif
    endif

    do dupl = 1, domain%num_duplicate

      ioffset = (dupl-1)*molecule%num_atoms

      do i = 1, prmtop%num_diheh
      
        if (prmtop%dihe_inc_hy(4,i) < 0) &
          cycle

        i1 =      prmtop%dihe_inc_hy(1,i)  / 3 + 1
        i2 =      prmtop%dihe_inc_hy(2,i)  / 3 + 1
        i3 = iabs(prmtop%dihe_inc_hy(3,i)) / 3 + 1
        i4 =      prmtop%dihe_inc_hy(4,i)  / 3 + 1

        i1 = i1 + ioffset
        i2 = i2 + ioffset
        i3 = i3 + ioffset
        i4 = i4 + ioffset

        icel1 = id_g2l(1,i1)
        icel2 = id_g2l(1,i4)

        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1,icel2)

          if (icel_local > 0 .and. icel_local <= ncel) then

            ndihe => dihe(icel_local)
            ndihe = ndihe + 1

            if (ndihe > MaxDihe) &
              call error_msg('Setup_Enefunc_Dihe> Too many dihedrals.') 

            list(1:4,ndihe,icel_local) = (/i1,i2,i3,i4/)
            force (  ndihe,icel_local) = &
                              prmtop%dihe_fcons_uniq(prmtop%dihe_inc_hy(5,i))
            phase (  ndihe,icel_local) = &
                              prmtop%dihe_phase_uniq(prmtop%dihe_inc_hy(5,i))
            period(  ndihe,icel_local) = &
                          int(prmtop%dihe_perio_uniq(prmtop%dihe_inc_hy(5,i)))
            if (prmtop%lscee_scale_factor .or. prmtop%lscnb_scale_factor ) then
              period(ndihe,icel_local) = period(  ndihe,icel_local) &
                              + prmtop%dihe_inc_hy(5,i)*notation
            endif

          end if

        end if

      end do

      do i = 1, prmtop%num_mdihea

        if (prmtop%dihe_wo_hy(4,i) < 0) &
          cycle

        i1 =      prmtop%dihe_wo_hy(1,i)  / 3 + 1
        i2 =      prmtop%dihe_wo_hy(2,i)  / 3 + 1
        i3 = iabs(prmtop%dihe_wo_hy(3,i)) / 3 + 1
        i4 =      prmtop%dihe_wo_hy(4,i)  / 3 + 1

        i1 = i1 + ioffset
        i2 = i2 + ioffset
        i3 = i3 + ioffset
        i4 = i4 + ioffset

        icel1 = id_g2l(1,i1)
        icel2 = id_g2l(1,i4)

        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1,icel2)

          if (icel_local > 0 .and. icel_local <= ncel) then

            ndihe => dihe(icel_local)
            ndihe = ndihe + 1

            if (ndihe > MaxDihe) &
              call error_msg('Setup_Enefunc_Dihe> Too many dihedrals.') 

            list(1:4,ndihe,icel_local) = (/i1,i2,i3,i4/)
            force (  ndihe,icel_local) = &
                              prmtop%dihe_fcons_uniq(prmtop%dihe_wo_hy(5,i))
            phase (  ndihe,icel_local) = &
                              prmtop%dihe_phase_uniq(prmtop%dihe_wo_hy(5,i))
            period(  ndihe,icel_local) = &
                          int(prmtop%dihe_perio_uniq(prmtop%dihe_wo_hy(5,i)))
            if (prmtop%lscee_scale_factor .or. prmtop%lscnb_scale_factor ) then
              period(ndihe,icel_local) = period(  ndihe,icel_local) &
                              + prmtop%dihe_wo_hy(5,i)*notation
            endif

          end if

        end if

      end do

    end do

    found = 0 
    do i = 1, ncel
      found = found + dihe(i)
    end do

#ifdef MPI
    call mpi_allreduce(found, enefunc%num_dihe_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_dihe_all = found
#endif

    return

  end subroutine setup_enefunc_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_impr
  !> @brief        define IMPROPER dihedral term in potential energy function
  !! @authors      NT
  !! @param[in]    prmtop   : AMBER parameter topology information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_impr(molecule, prmtop, domain, enefunc)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: dupl, ioffset
    integer                  :: i, found
    integer                  :: i1, i2, i3, i4, icel1, icel2, icel_local

    real(wp),           pointer :: force(:,:), phase(:,:)
    integer,            pointer :: impr(:), list(:,:,:), period(:,:)
    integer,            pointer :: ncel
    integer(int2),      pointer :: cell_pair(:,:)
    integer(int2),      pointer :: id_g2l(:,:)
    integer,            pointer :: nimpr
    integer,            pointer :: notation


    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l

    impr      => enefunc%num_improper
    list      => enefunc%impr_list
    force     => enefunc%impr_force_const
    phase     => enefunc%impr_phase
    period    => enefunc%impr_periodicity
    notation  => enefunc%notation_14types

    do dupl = 1, domain%num_duplicate

      ioffset = molecule%num_atoms*(dupl-1)

      do i = 1, prmtop%num_diheh
      
        if (prmtop%dihe_inc_hy(4,i) >= 0) &
          cycle

        i1 =      prmtop%dihe_inc_hy(1,i)  / 3 + 1
        i2 =      prmtop%dihe_inc_hy(2,i)  / 3 + 1
        i3 = iabs(prmtop%dihe_inc_hy(3,i)) / 3 + 1
        i4 = iabs(prmtop%dihe_inc_hy(4,i)) / 3 + 1

        i1 = i1 + ioffset
        i2 = i2 + ioffset
        i3 = i3 + ioffset
        i4 = i4 + ioffset

        icel1 = id_g2l(1,i1)
        icel2 = id_g2l(1,i4)

        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1,icel2)

          if (icel_local > 0 .and. icel_local <= ncel) then

            nimpr => impr(icel_local)
            nimpr = nimpr + 1

            if (nimpr > MaxImpr) &
              call error_msg('Setup_Enefunc_Impr> Too many impropers.') 

            list(1:4,nimpr,icel_local) = (/i1,i2,i3,i4/)
            force (  nimpr,icel_local) = &
                              prmtop%dihe_fcons_uniq(prmtop%dihe_inc_hy(5,i))
            phase (  nimpr,icel_local) = &
                              prmtop%dihe_phase_uniq(prmtop%dihe_inc_hy(5,i))
            period(  nimpr,icel_local) = &
                         int(prmtop%dihe_perio_uniq(prmtop%dihe_inc_hy(5,i)))
            if (prmtop%lscee_scale_factor .or. prmtop%lscnb_scale_factor ) then
              period(nimpr,icel_local) = period(  nimpr,icel_local) &
                              + prmtop%dihe_inc_hy(5,i)*notation
            endif

          end if

        end if

      end do

      do i = 1, prmtop%num_mdihea

        if (prmtop%dihe_wo_hy(4,i) >= 0) &
          cycle

        i1 =      prmtop%dihe_wo_hy(1,i)  / 3 + 1
        i2 =      prmtop%dihe_wo_hy(2,i)  / 3 + 1
        i3 = iabs(prmtop%dihe_wo_hy(3,i)) / 3 + 1
        i4 = iabs(prmtop%dihe_wo_hy(4,i)) / 3 + 1

        i1 = i1 + ioffset
        i2 = i2 + ioffset
        i3 = i3 + ioffset
        i4 = i4 + ioffset

        icel1 = id_g2l(1,i1)
        icel2 = id_g2l(1,i4)

        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1,icel2)

          if (icel_local > 0 .and. icel_local <= ncel) then

            nimpr => impr(icel_local)
            nimpr = nimpr + 1

            if (nimpr > MaxImpr) &
              call error_msg('Setup_Enefunc_Impr> Too many impropers.') 

            list(1:4,nimpr,icel_local) = (/i1,i2,i3,i4/)
            force (  nimpr,icel_local) = &
                              prmtop%dihe_fcons_uniq(prmtop%dihe_wo_hy(5,i))
            phase (  nimpr,icel_local) = &
                              prmtop%dihe_phase_uniq(prmtop%dihe_wo_hy(5,i))
            period(  nimpr,icel_local) = &
                            int(prmtop%dihe_perio_uniq(prmtop%dihe_wo_hy(5,i)))

            if (prmtop%lscee_scale_factor .or. prmtop%lscnb_scale_factor ) then
              period(nimpr,icel_local) = period(  nimpr,icel_local) &
                              + prmtop%dihe_wo_hy(5,i)*notation
            endif

          end if

        end if

      end do

    end do

    found = 0 
    do i = 1, ncel
      found = found + impr(i)
    end do

#ifdef MPI
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
  !! @param[in]    prmtop      : AMBER parameter topology information
  !! @param[in]    constraints : constraints information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_nonb(prmtop, molecule, constraints, domain, enefunc)

    ! formal arguments
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_molecule),        intent(in)    :: molecule
    type(s_constraints),     intent(in)    :: constraints
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, nnonb, ncel, ndihetype, ndihe
    integer                  :: nb_par_idx
    integer                  :: k, ix, jx, cls_local

    integer,  allocatable    :: check_cls(:)
    integer,  allocatable    :: atmcls_map_g2l(:), atmcls_map_l2g(:)
    real(wp), allocatable    :: nb14_lj6(:,:), nb14_lj12(:,:)
    real(wp), allocatable    :: nonb_lj6(:,:), nonb_lj12(:,:)

    ndihetype = size(prmtop%scee_scale_fact)
    call alloc_enefunc(enefunc, EneFuncAMBERScale, ndihetype, 0)

    ELECOEF          = ELECOEF_AMBER

    ndihe = 0

    do i = 1, prmtop%num_diheh
      if (prmtop%dihe_inc_hy(4,i) < 0) &
        cycle

      ndihe = ndihe + 1

      if (prmtop%lscee_scale_factor) then
        enefunc%dihe_scee(prmtop%dihe_inc_hy(5,i)) = &
         1.0_wp/prmtop%scee_scale_fact(prmtop%dihe_inc_hy(5,i))
      else
        enefunc%dihe_scee(0) = 1.0_wp/1.2_wp
      end if

      if (prmtop%lscnb_scale_factor) then
        enefunc%dihe_scnb(prmtop%dihe_inc_hy(5,i)) = &
         1.0_wp/prmtop%scnb_scale_fact(prmtop%dihe_inc_hy(5,i))
      else
        enefunc%dihe_scnb(0) = 1.0_wp/2.0_wp
      end if

    end do

    do i = 1, prmtop%num_mdihea
      if (prmtop%dihe_wo_hy(4,i) < 0) &
        cycle

      ndihe = ndihe + 1

      if (prmtop%lscee_scale_factor) then
        enefunc%dihe_scee(prmtop%dihe_wo_hy(5,i)) =  &
          1.0_wp/prmtop%scee_scale_fact(prmtop%dihe_wo_hy(5,i))
      else
        enefunc%dihe_scee(0) = 1.0_wp/1.2_wp
      end if

      if (prmtop%lscnb_scale_factor) then
        enefunc%dihe_scnb(prmtop%dihe_wo_hy(5,i)) = &
          1.0_wp/prmtop%scnb_scale_fact(prmtop%dihe_wo_hy(5,i))
      else
        enefunc%dihe_scnb(0) = 1.0_wp/2.0_wp
      end if

    end do

    nnonb = prmtop%num_types

    allocate(check_cls(nnonb),             &
             atmcls_map_g2l(nnonb),        &
             atmcls_map_l2g(nnonb),        &
             nb14_lj6      (nnonb, nnonb), &
             nb14_lj12     (nnonb, nnonb), &
             nonb_lj6      (nnonb, nnonb), &
             nonb_lj12     (nnonb, nnonb))

    check_cls(1:nnonb)          = 0
    nb14_lj6 (1:nnonb, 1:nnonb) = 0.0_wp
    nb14_lj12(1:nnonb, 1:nnonb) = 0.0_wp
    nonb_lj6 (1:nnonb, 1:nnonb) = 0.0_wp
    nonb_lj12(1:nnonb, 1:nnonb) = 0.0_wp

    do i = 1, nnonb
      do j = 1, nnonb

        nb_par_idx = prmtop%nb_par_idx(nnonb*(i-1)+j)

        if (nb_par_idx < 0) then
          nb14_lj12(i,j) = 0.0_wp
          nb14_lj6 (i,j) = 0.0_wp
          nonb_lj12(i,j) = 0.0_wp
          nonb_lj6 (i,j) = 0.0_wp
        else
          nb14_lj12(i,j) = prmtop%lennarda(nb_par_idx)
          nb14_lj6 (i,j) = prmtop%lennardb(nb_par_idx)
          nonb_lj12(i,j) = prmtop%lennarda(nb_par_idx)
          nonb_lj6 (i,j) = prmtop%lennardb(nb_par_idx)
        end if

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
      endif
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
    domain%water%atom_cls_no(1:3)  &
      = atmcls_map_g2l(domain%water%atom_cls_no(1:3))
    enefunc%table%atom_cls_no_O = atmcls_map_g2l(enefunc%table%atom_cls_no_O)
    enefunc%table%atom_cls_no_H = atmcls_map_g2l(enefunc%table%atom_cls_no_H)
    if (constraints%tip4 .or. enefunc%table%tip4) then
      enefunc%table%atom_cls_no_D = atmcls_map_g2l(enefunc%table%atom_cls_no_D)
      domain%water%atom_cls_no(4) = atmcls_map_g2l(domain%water%atom_cls_no(4))
    end if

    deallocate(check_cls,     &
               atmcls_map_g2l,&
               atmcls_map_l2g,&
               nb14_lj6,      &
               nb14_lj12,     &
               nonb_lj6,      &
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

end module sp_enefunc_amber_mod
