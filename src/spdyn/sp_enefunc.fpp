!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_enefunc_mod
!> @brief   define potential energy functions in each domain
!! @authors Jaewoon Jung (JJ), Yuji Sugita (YS), Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_enefunc_mod

  use sp_enefunc_gromacs_mod
  use sp_enefunc_amber_mod
  use sp_enefunc_charmm_mod
  use sp_enefunc_restraints_mod
  use sp_enefunc_localres_mod
  use sp_enefunc_table_mod
  use sp_communicate_mod
  use sp_migration_mod
  use sp_energy_mod
  use sp_restraints_str_mod
  use sp_constraints_str_mod
  use sp_enefunc_str_mod
  use sp_energy_str_mod
  use sp_domain_str_mod
  use molecules_str_mod
  use fileio_localres_mod
  use fileio_grotop_mod
  use fileio_prmtop_mod
  use fileio_par_mod
  use timers_mod
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

  ! subroutines
  public  :: define_enefunc
  public  :: define_enefunc_pio
  public  :: update_enefunc
  private :: setup_enefunc_bond_pio
  private :: setup_enefunc_bond_constraint_pio
  private :: setup_enefunc_angl_pio
  private :: setup_enefunc_dihe_pio
  private :: setup_enefunc_rb_dihe_pio
  private :: setup_enefunc_impr_pio
  private :: setup_enefunc_cmap_pio
  private :: setup_enefunc_nonb_pio
  private :: setup_enefunc_dispcorr
  private :: check_bonding
  ! FEP
  public  :: update_enefunc_fep
  private :: setup_enefunc_dispcorr_fep
  private :: set_fepgrp_nonb

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc
  !> @brief        a driver subroutine for defining potential energy functions
  !! @authors      YS, JJ, CK
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    par         : CHARMM PAR information
  !! @param[in]    prmtop      : AMBER parameter topology information
  !! @param[in]    grotop      : GROMACS parameter topology information
  !! @param[in]    localres    : local restraint information
  !! @param[in]    molecule    : molecule information
  !! @param[inout] constraints : constraints information
  !! @param[in]    restraints  : restraints information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : energy potential functions information
  !! @param[inout] comm        : communication information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc(ene_info, par, prmtop, grotop, localres, molecule, &
                            constraints, restraints, domain, enefunc, comm)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_par),             intent(in)    :: par
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_grotop),          intent(in)    :: grotop
    type(s_localres),        intent(in)    :: localres
    type(s_molecule),        intent(in)    :: molecule
    type(s_constraints),     intent(inout) :: constraints
    type(s_restraints),      intent(in)    :: restraints
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_comm),            intent(inout) :: comm   


    enefunc%forcefield        = ene_info%forcefield
    enefunc%output_style      = ene_info%output_style
    enefunc%table%water_table = ene_info%table .and. &
                                (ene_info%water_model(1:4) .ne. 'NONE' .and. &
                                (.not. ene_info%nonb_limiter))

    enefunc%switchdist        = ene_info%switchdist
    enefunc%cutoffdist        = ene_info%cutoffdist
    enefunc%pairlistdist      = ene_info%pairlistdist
    enefunc%dielec_const      = ene_info%dielec_const
    enefunc%force_switch      = ene_info%vdw_force_switch
    enefunc%vdw_shift         = ene_info%vdw_shift

    enefunc%vdw               = ene_info%vdw
    enefunc%pme_use           = ene_info%electrostatic == ElectrostaticPME
    enefunc%pme_alpha         = ene_info%pme_alpha
    enefunc%pme_ngrid_x       = ene_info%pme_ngrid_x
    enefunc%pme_ngrid_y       = ene_info%pme_ngrid_y
    enefunc%pme_ngrid_z       = ene_info%pme_ngrid_z
    enefunc%pme_nspline       = ene_info%pme_nspline
    enefunc%pme_scheme        = ene_info%pme_scheme
    enefunc%pme_max_spacing   = ene_info%pme_max_spacing
    enefunc%dispersion_corr   = ene_info%dispersion_corr
    enefunc%contact_check     = ene_info%contact_check
    enefunc%nonb_limiter      = ene_info%nonb_limiter
    enefunc%minimum_contact   = ene_info%minimum_contact
    enefunc%err_minimum_contact = ene_info%err_minimum_contact

    enefunc%efield(1)         = ene_info%efield_x           
    enefunc%efield(2)         = ene_info%efield_y           
    enefunc%efield(3)         = ene_info%efield_z           
    enefunc%efield_virial     = ene_info%efield_virial
    enefunc%efield_normal     = ene_info%efield_normal

    ! For Vacuum
    enefunc%vacuum            = ene_info%vacuum

    ! FEP
    if (domain%fep_use) then
      enefunc%fep_topology = molecule%fep_topology
      call set_fepgrp_nonb(enefunc)
    end if

    if (abs(enefunc%efield(1)) > EPS .or. abs(enefunc%efield(2)) > EPS .or. &
        abs(enefunc%efield(3)) > EPS) then
      enefunc%use_efield = .true.
    end if

    if (ene_info%structure_check == StructureCheckDomain) then
      enefunc%pairlist_check = .true.
      enefunc%bonding_check  = .true.
    end if

    ! charmm
    !
    if (par%num_bonds > 0) then

      call define_enefunc_charmm(ene_info, par, localres, molecule, &
                                 constraints, restraints, domain, enefunc)

    ! amber
    !
    else if (prmtop%num_atoms > 0) then
      if (localres%num_funcs > 0) &
         call error_msg('Define_Enefunc> Localres is not available in this FF')

      call define_enefunc_amber (ene_info, prmtop, molecule, &
                                 constraints, restraints, domain, enefunc)

    ! gromacs
    !
    else if (grotop%num_atomtypes > 0) then
      if (localres%num_funcs > 0) &
         call error_msg('Define_Enefunc> Localres is not available in this FF')

      call define_enefunc_gromacs(ene_info, grotop, molecule, &
                                 constraints, restraints, domain, enefunc)

    end if

    ! dispersion correction
    !
    if (domain%fep_use) then
      ! FEP
      call setup_enefunc_dispcorr_fep(ene_info, domain, enefunc, constraints)
    else
      call setup_enefunc_dispcorr(ene_info, domain, enefunc, constraints)
    end if

    ! bonding_checker
    !
    if (ene_info%structure_check /= StructureCheckNone)  &
      call check_bonding(enefunc, domain)

    return

  end subroutine define_enefunc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc_pio
  !> @brief        a driver subroutine for defining potential energy functions
  !! @authors      JJ
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    localres    : local restraint information
  !! @param[inout] comm        : communication information
  !! @param[inout] constraints : constraints information
  !! @param[inout] restraints  : restraints information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc_pio(ene_info, localres, comm, constraints, &
                                restraints, domain, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_localres),        intent(in)    :: localres
    type(s_comm),            intent(inout) :: comm
    type(s_constraints),     intent(inout) :: constraints
    type(s_restraints),      intent(inout) :: restraints 
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: ncel, ncelb


    enefunc%forcefield        = ene_info%forcefield
    enefunc%output_style      = ene_info%output_style
    enefunc%table%water_table = ene_info%table .and. &
                                (ene_info%water_model(1:4) .ne. 'NONE')

    enefunc%switchdist        = ene_info%switchdist
    enefunc%cutoffdist        = ene_info%cutoffdist
    enefunc%pairlistdist      = ene_info%pairlistdist
    enefunc%dielec_const      = ene_info%dielec_const

    enefunc%pme_use           = ene_info%electrostatic == ElectrostaticPME
    enefunc%vdw               = ene_info%vdw
    enefunc%pme_alpha         = ene_info%pme_alpha
    enefunc%pme_ngrid_x       = ene_info%pme_ngrid_x
    enefunc%pme_ngrid_y       = ene_info%pme_ngrid_y
    enefunc%pme_ngrid_z       = ene_info%pme_ngrid_z
    enefunc%pme_nspline       = ene_info%pme_nspline
    enefunc%pme_max_spacing   = ene_info%pme_max_spacing
    enefunc%pme_scheme        = ene_info%pme_scheme
    enefunc%dispersion_corr   = ene_info%dispersion_corr
    enefunc%contact_check     = ene_info%contact_check
    enefunc%nonb_limiter      = ene_info%nonb_limiter
    enefunc%minimum_contact   = ene_info%minimum_contact
    enefunc%err_minimum_contact = ene_info%err_minimum_contact

    ! base
    !
    ncel  = domain%num_cell_local
    ncelb = domain%num_cell_local + domain%num_cell_boundary

    call alloc_enefunc(enefunc, EneFuncBase, ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncBond, ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncAngl, ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncDihe, ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncImpr, ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncBondCell, ncel, ncelb)

    if (constraints%rigid_bond) then

      ! bond
      !
      call setup_enefunc_bond_constraint_pio(domain, constraints, enefunc)

    else

      if (.not. constraints%fast_water) then
        enefunc%table%water_bond_calc = .true.
        if (enefunc%table%HOH_force > 0.0_wp) &
          enefunc%table%water_angle_calc = .true.
      end if

      ! bond
      !
      call setup_enefunc_bond_pio(domain, enefunc)

    end if

    ! angle
    !
    call setup_enefunc_angl_pio(domain, enefunc)

    ! dihedral
    !
    call setup_enefunc_dihe_pio(domain, enefunc)

    ! Ryckaert-Bellemans dihedral
    !
    call setup_enefunc_rb_dihe_pio(domain, enefunc)

    ! improper
    !
    call setup_enefunc_impr_pio(domain, enefunc)

    ! cmap
    !
    call setup_enefunc_cmap_pio(ene_info, domain, enefunc)

    ! restraint
    !
    call setup_enefunc_restraints_pio(restraints, domain, enefunc)

    ! reassign bond information
    !
    call update_enefunc_pio(domain, comm, enefunc, constraints)
    call dealloc_enefunc(enefunc, EneFuncBondCell)

    ! nonbonded
    !
    call setup_enefunc_nonb_pio(ene_info, constraints, domain, enefunc)

    ! dispersion correction
    !
    call setup_enefunc_dispcorr(ene_info, domain, enefunc, constraints)

    ! lookup table
    !
    if(ene_info%table) &
    call setup_enefunc_table(ene_info, enefunc)

    ! restraint
    !
    call setup_enefunc_localres(localres, domain, enefunc)

    if (ene_info%structure_check /= StructureCheckNone)  &
      call check_bonding(enefunc, domain)

    call dealloc_domain(domain, DomainDynvar_pio)
    call dealloc_enefunc(enefunc, EneFuncBase_pio)
    call dealloc_enefunc(enefunc, EneFuncBond_pio)
    call dealloc_enefunc(enefunc, EneFuncAngl_pio)
    call dealloc_enefunc(enefunc, EneFuncDihe_pio)
    call dealloc_enefunc(enefunc, EneFuncRBDihe_pio)
    call dealloc_enefunc(enefunc, EneFuncImpr_pio)
    call dealloc_enefunc(enefunc, EneFuncCmap_pio)

    ! write summary of energy function
    !
    if (main_rank) then
      write(MsgOut,'(A)') &
           'Define_Enefunc_Pio> Number of Interactions in Each Term'
      write(MsgOut,'(A20,I10,A20,I10)')                  &
           '  bond_ene        = ', enefunc%num_bond_all, &
           '  angle_ene       = ', enefunc%num_angl_all
      write(MsgOut,'(A20,I10,A20,I10)')                  &
           '  torsion_ene     = ', enefunc%num_dihe_all, &
           '  rb_torsion_ene  = ', enefunc%num_rb_dihe_all
      write(MsgOut,'(A20,I10,A20,I10)')                  &
           '  improper_ene    = ', enefunc%num_impr_all, &
           '  cmap_ene        = ', enefunc%num_cmap_all
      write(MsgOut,'(A20,I10,A20,I10)')                  &
           '  nb_exclusions   = ', enefunc%num_excl_all, &
           '  nb14_calc       = ', enefunc%num_nb14_all
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine define_enefunc_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_enefunc
  !> @brief        a driver subroutine for updating potential energy functions
  !! @authors      JJ
  !! @param[inout] domain      : domain information
  !! @param[inout] comm        : communication information
  !! @param[inout] enefunc     : energy potential functions information
  !! @param[inout] constraints : constraints information [optional]
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_enefunc(domain, comm, enefunc, constraints)

    ! formal arguments
    type(s_domain),                intent(inout) :: domain
    type(s_comm),                  intent(inout) :: comm
    type(s_enefunc),               intent(inout) :: enefunc
    type(s_constraints), optional, intent(inout) :: constraints

    ! local variables
    logical                        :: first


    ! sending the bonding information to other domain
    !

    ! bond
    !
    call update_outgoing_enefunc_bond(domain, enefunc)

    ! enm
    !
    if (enefunc%enm_use) then
      call update_enefunc_enm(domain, enefunc)
      call update_cell_size_enm(domain, enefunc, comm)
      call update_cell_boundary_enm(domain, enefunc, comm)
    end if

    ! angle
    !
    call update_outgoing_enefunc_angl(domain, enefunc)

    ! dihedral
    !
    call update_outgoing_enefunc_dihe(domain, enefunc)

    ! Ryckaert-Bellemans dihedral
    !
    call update_outgoing_enefunc_rb_dihe(domain, enefunc)

    ! improper dihedral
    !
    call update_outgoing_enefunc_impr(domain, enefunc)

    ! cmap
    !
    call update_outgoing_enefunc_cmap(domain, enefunc)

    ! restraint
    !
    call update_outgoing_enefunc_restraint(domain, enefunc)

    ! fitting
    !
    call update_outgoing_enefunc_fitting(domain, enefunc)

    ! communicate neighbour domain
    !
    call communicate_bond(domain, comm, enefunc)


    ! bond
    !
    call update_incoming_enefunc_bond(domain, enefunc)

    ! angle
    !
    call update_incoming_enefunc_angl(domain, enefunc)

    ! dihedral
    !
    call update_incoming_enefunc_dihe(domain, enefunc)

    ! Ryckaert-Bellemans dihedral
    !
    call update_incoming_enefunc_rb_dihe(domain, enefunc)

    ! improper dihedral
    !
    call update_incoming_enefunc_impr(domain, enefunc)

    ! cmap
    !
    call update_incoming_enefunc_cmap(domain, enefunc)

    ! restraint
    !
    call update_incoming_enefunc_restraint(domain, enefunc)

    ! fitting
    !
    call update_incoming_enefunc_fitting(domain, enefunc)

    ! re-count nonbond exclusion list
    !

    first = .false.

    if (constraints%rigid_bond) then

      call count_nonb_excl(first, .true., constraints, domain, enefunc)

    else

      call count_nonb_excl(first, .false., constraints, domain, enefunc)

    end if

    if (enefunc%bonding_check) call check_bonding(enefunc, domain)

    return

  end subroutine update_enefunc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_enefunc_pio
  !> @brief        a driver subroutine for updating potential energy functions
  !! @authors      JJ
  !! @param[inout] domain      : domain information
  !! @param[inout] comm        : communication information
  !! @param[inout] enefunc     : energy potential functions information
  !! @param[inout] constraints : constraints information [optional]
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_enefunc_pio(domain, comm, enefunc, constraints)

    ! formal arguments
    type(s_domain),                intent(inout) :: domain
    type(s_comm),                  intent(inout) :: comm
    type(s_enefunc),               intent(inout) :: enefunc
    type(s_constraints), optional, intent(inout) :: constraints

    ! local variables
    logical                        :: first


    ! sending the bonding information to other domain
    !

    ! bond
    !
    call update_outgoing_enefunc_bond(domain, enefunc)

    ! angle
    !
    call update_outgoing_enefunc_angl(domain, enefunc)

    ! dihedral
    !
    call update_outgoing_enefunc_dihe(domain, enefunc)

    ! Ryckaert-Bellemans dihedral
    !
    call update_outgoing_enefunc_rb_dihe(domain, enefunc)

    ! improper dihedral
    !
    call update_outgoing_enefunc_impr(domain, enefunc)

    ! cmap
    !
    call update_outgoing_enefunc_cmap(domain, enefunc)

    ! restraint
    !
    call update_outgoing_enefunc_restraint(domain, enefunc)

    ! fitting
    !
    call update_outgoing_enefunc_fitting(domain, enefunc)

    ! communicate neighbour domain
    !
    call communicate_bond(domain, comm, enefunc)


    ! bond
    !
    call update_incoming_enefunc_bond(domain, enefunc)

    ! angle
    !
    call update_incoming_enefunc_angl(domain, enefunc)

    ! dihedral
    !
    call update_incoming_enefunc_dihe(domain, enefunc)

    ! Ryckaert-Bellemans dihedral
    !
    call update_incoming_enefunc_rb_dihe(domain, enefunc)

    ! improper dihedral
    !
    call update_incoming_enefunc_impr(domain, enefunc)

    ! cmap
    !
    call update_incoming_enefunc_cmap(domain, enefunc)

    ! restraint
    !
    call update_incoming_enefunc_restraint(domain, enefunc)

    ! fitting
    !
    call update_incoming_enefunc_fitting(domain, enefunc)

    return

  end subroutine update_enefunc_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_bond_pio
  !> @brief        define BOND term for each cell in potential energy function
  !! @authors      NT
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond_pio(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    integer                  :: i, ix, ic, icel, found, ncell, ncell_pio
    integer                  :: file_num, file_tot_num
    integer,         pointer :: nwater(:)
    integer,         pointer :: bond(:), bond_list(:,:,:)
    real(wp),        pointer :: bond_force(:,:), bond_dist(:,:)
    integer,         pointer :: bond_pio(:,:), bond_list_pio(:,:,:,:)
    integer,         pointer :: bond_pbc_pio(:,:,:)
    real(wp),        pointer :: bond_force_pio(:,:,:), bond_dist_pio(:,:,:)
    integer,         pointer :: cell_l2g_pio(:,:)
    integer,         pointer :: bond_pbc(:,:)
    integer(int2),   pointer :: cell_g2l(:)


    cell_l2g_pio   => domain%cell_l2g_pio
    cell_g2l       => domain%cell_g2l
    nwater         => domain%num_water
    bond           => enefunc%num_bond
    bond_list      => enefunc%bond_list
    bond_force     => enefunc%bond_force_const
    bond_dist      => enefunc%bond_dist_min
    bond_pbc       => enefunc%bond_pbc
    bond_pio       => enefunc%num_bond_pio
    bond_list_pio  => enefunc%bond_list_pio
    bond_pbc_pio   => enefunc%bond_pbc_pio
    bond_force_pio => enefunc%bond_force_const_pio
    bond_dist_pio  => enefunc%bond_dist_min_pio
    ncell          =  domain%num_cell_local
    ncell_pio      =  domain%ncell_local_pio

    file_tot_num = domain%file_tot_num
    found = 0

    do file_num = 1, file_tot_num

      do icel = 1, ncell_pio

        ic = cell_l2g_pio(icel,file_num)
        i = cell_g2l(ic)

        if (i /= 0) then

          bond(i) = bond_pio(icel,file_num)

          do ix = 1, bond(i)
            bond_list (1,ix,i) = bond_list_pio (1,ix,icel,file_num)
            bond_list (2,ix,i) = bond_list_pio (2,ix,icel,file_num)
            bond_force(  ix,i) = bond_force_pio(  ix,icel,file_num)
            bond_dist (  ix,i) = bond_dist_pio (  ix,icel,file_num)
            bond_pbc  (  ix,i) = bond_pbc_pio  (  ix,icel,file_num)
          end do  
          found = found + bond(i) + 2*nwater(i)  !! 2 is two O-H bond from water molecules
          if (bond(i) > MaxBond) &
            call error_msg('Setup_Enefunc_Bond_Pio> Too many bonds.')

        end if

      end do
    end do

    enefunc%table%water_bond_calc = .true.

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_bond_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_bond_all = found
#endif

    return

  end subroutine setup_enefunc_bond_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_bond_constraint_pio
  !> @brief        define BOND term for each cell in potential energy function
  !! @authors      NT
  !! @param[in]    domain      : domain information
  !! @param[in]    constraints : constraints information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond_constraint_pio(domain, constraints, enefunc)

    ! formal arguments
    type(s_domain),      target, intent(in)    :: domain
    type(s_constraints), target, intent(inout) :: constraints
    type(s_enefunc),     target, intent(inout) :: enefunc

    ! local variable
    logical                  :: calc_bond
    integer                  :: i, ix, j, k, ih, ih1, ih2, i1, i2, ia, ib, ig
    integer                  :: file_num, file_tot_num
    integer                  :: found, ncel, connect, pbond
    integer                  :: icel, ncell, ncell_pio
    integer,         pointer :: bond(:), bond_pio(:,:)
    integer,         pointer :: bond_list(:,:,:), bond_list_pio(:,:,:,:)
    integer,         pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)
    integer,         pointer :: id_l2g_sol(:,:)
    integer,         pointer :: cell_l2g_pio(:,:)
    integer(int2),   pointer :: cell_g2l(:)
    real(wp),        pointer :: bond_force(:,:), bond_force_pio(:,:,:)
    real(wp),        pointer :: bond_dist(:,:), bond_dist_pio(:,:,:)


    cell_l2g_pio   => domain%cell_l2g_pio
    cell_g2l       => domain%cell_g2l
    bond           => enefunc%num_bond
    bond_list      => enefunc%bond_list
    bond_force     => enefunc%bond_force_const
    bond_dist      => enefunc%bond_dist_min
    bond_pio       => enefunc%num_bond_pio
    bond_list_pio  => enefunc%bond_list_pio
    bond_force_pio => enefunc%bond_force_const_pio
    bond_dist_pio  => enefunc%bond_dist_min_pio

    id_l2g_sol     => domain%id_l2g_solute

    HGr_local      => constraints%HGr_local
    HGr_bond_list  => constraints%HGr_bond_list

    found = 0
    file_tot_num = domain%file_tot_num

    ncell     = domain%num_cell_local
    ncell_pio = domain%ncell_local_pio

    bond(1:ncell) = 0

    do file_num = 1, file_tot_num

      do icel = 1, ncell_pio

        ig = cell_l2g_pio(icel,file_num)
        i  = cell_g2l(ig)

        if (i /= 0) then

          do ix = 1, bond_pio(icel,file_num)

            calc_bond = .true.
            i1 = bond_list_pio(1,ix,icel,file_num)
            i2 = bond_list_pio(2,ix,icel,file_num)
        
            do j = 1, constraints%connect
              do k = 1, HGr_local(j,i)
                ia = HGr_bond_list(1,k,j,i)
                ih1 = id_l2g_sol(ia,i)
                do ih = 1, j
                  ib = HGr_bond_list(ih+1,k,j,i)
                  ih2 = id_l2g_sol(ib,i)
                  if ((ih1 == i1 .and. ih2 == i2) .or. &
                      (ih1 == i2 .and. ih2 == i1)) then
                    calc_bond = .false.
                    exit
                  end if
                end do
                if (.not.calc_bond) exit
              end do
              if (.not.calc_bond) exit
            end do
        
            if (calc_bond) then
              bond(i) = bond(i) + 1
              found = found + 1
              pbond = bond(i)
              bond_list(1,pbond,i) = i1
              bond_list(2,pbond,i) = i2
              bond_force(pbond,i) = bond_force_pio(ix,icel,file_num)
              bond_dist (pbond,i) = bond_dist_pio(ix,icel,file_num)
            end if

!           found = found + bond(i)
            if (bond(i) > MaxBond) call error_msg &
              ('Setup_Enefunc_Bond_Constraint_Pio> Too many bonds.')

          end do

        end if

      end do
    end do
#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_bond_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_bond_all = found
#endif

    ! count # of constraints bonds
    HGr_local => constraints%HGr_local

    ncel    = domain%num_cell_local
    connect = constraints%connect

    found = 0

    do i = 1, ncel
      do j = 1, connect
        do k = 1, HGr_local(j,i)
          do ih = 1, j
            found = found + 1
          end do
        end do
      end do
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, constraints%num_bonds, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    constraints%num_bonds = found
#endif

    return

  end subroutine setup_enefunc_bond_constraint_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_angl_pio
  !> @brief        define ANGLE term for each cell in potential energy function
  !! @authors      NT
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_angl_pio(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: i, ix, icel, ic, found
    integer                  :: file_num, file_tot_num
    integer,         pointer :: angle(:), angle_pio(:,:)
    integer,         pointer :: list(:,:,:), list_pio(:,:,:,:)
    integer,         pointer :: angl_pbc(:,:,:), angl_pbc_pio(:,:,:,:)
    integer,         pointer :: nwater(:)
    integer,         pointer :: cell_l2g_pio(:,:)
    integer(int2),   pointer :: cell_g2l(:)
    real(wp),        pointer :: force(:,:), theta(:,:)
    real(wp),        pointer :: force_pio(:,:,:), theta_pio(:,:,:)
    real(wp),        pointer :: ubforce(:,:), ubrmin(:,:)
    real(wp),        pointer :: ubforce_pio(:,:,:), ubrmin_pio(:,:,:)


    angle        => enefunc%num_angle
    list         => enefunc%angle_list
    force        => enefunc%angle_force_const
    theta        => enefunc%angle_theta_min
    ubforce      => enefunc%urey_force_const
    ubrmin       => enefunc%urey_rmin
    angl_pbc     => enefunc%angle_pbc
    angl_pbc_pio => enefunc%angle_pbc_pio
    angle_pio    => enefunc%num_angle_pio
    list_pio     => enefunc%angle_list_pio
    force_pio    => enefunc%angle_force_const_pio
    theta_pio    => enefunc%angle_theta_min_pio
    ubforce_pio  => enefunc%urey_force_const_pio
    ubrmin_pio   => enefunc%urey_rmin_pio
    nwater       => domain%num_water
    cell_l2g_pio => domain%cell_l2g_pio
    cell_g2l     => domain%cell_g2l

    found = 0
    file_tot_num = domain%file_tot_num

    do file_num = 1, domain%file_tot_num

      do icel = 1, domain%ncell_local_pio

        ic = cell_l2g_pio(icel,file_num)
        i  = cell_g2l(ic)

        if (i /= 0) then

          angle(i)   = angle_pio(icel,file_num)

          do ix = 1, angle(i)
            list(1:3,ix,i) = list_pio(1:3,ix,icel,file_num)
            force(ix,i)   = force_pio(ix,icel,file_num)
            theta(ix,i)   = theta_pio(ix,icel,file_num)
            ubforce(ix,i) = ubforce_pio(ix,icel,file_num)
            ubrmin(ix,i)  = ubrmin_pio(ix,icel,file_num)
          end do
    
          if (enefunc%table%water_angle_calc) then
            found = found + angle(i) + nwater(i)  ! angle from water
          else
            found = found + angle(i)
          end if
 
          if (angle(i) > MaxAngle) &
            call error_msg('Setup_Enefunc_Angl_Pio> Too many angles.')

        end if

      end do
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_angl_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_angl_all = found
#endif

    return

  end subroutine setup_enefunc_angl_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_dihe_pio
  !> @brief        define DIHEDRAL term in potential energy function
  !! @authors      NT
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_dihe_pio(domain, enefunc)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variables
    integer                   :: i, ix, ic, icel, found, ncell, ncell_pio
    integer                   :: file_num, file_tot_num
    integer,          pointer :: dihedral(:), dihedral_pio(:,:)
    integer,          pointer :: list(:,:,:), list_pio(:,:,:,:)
    integer,          pointer :: period(:,:), period_pio(:,:,:)
    integer,          pointer :: cell_l2g_pio(:,:)
    integer,          pointer :: dihe_pbc(:,:,:), dihe_pbc_pio(:,:,:,:)
    integer(int2),    pointer :: cell_g2l(:)
    real(wp),         pointer :: force(:,:), force_pio(:,:,:)
    real(wp),         pointer :: phase(:,:), phase_pio(:,:,:)


    cell_l2g_pio  => domain%cell_l2g_pio
    cell_g2l      => domain%cell_g2l
    dihedral      => enefunc%num_dihedral
    list          => enefunc%dihe_list
    force         => enefunc%dihe_force_const
    period        => enefunc%dihe_periodicity
    phase         => enefunc%dihe_phase
    dihe_pbc      => enefunc%dihe_pbc
    dihedral_pio  => enefunc%num_dihedral_pio
    list_pio      => enefunc%dihe_list_pio
    force_pio     => enefunc%dihe_force_const_pio
    period_pio    => enefunc%dihe_periodicity_pio
    phase_pio     => enefunc%dihe_phase_pio
    dihe_pbc_pio  => enefunc%dihe_pbc_pio

    ncell         = domain%num_cell_local
    ncell_pio     = domain%ncell_local_pio

    found = 0

    do file_num = 1, domain%file_tot_num

      do icel = 1, ncell_pio

        ic = cell_l2g_pio(icel,file_num)
        i  = cell_g2l(ic)

        if (i /= 0) then

          dihedral(i) = dihedral_pio(icel,file_num)
          do ix = 1, dihedral(i)
            list(1:4,ix,i) = list_pio(1:4,ix,icel,file_num)
            force(ix,i) = force_pio(ix,icel,file_num)
            period(ix,i) = period_pio(ix,icel,file_num)
            phase(ix,i) = phase_pio(ix,icel,file_num)
            dihe_pbc(1:3,ix,i) = dihe_pbc_pio(1:3,ix,icel,file_num)
          end do
          found = found + dihedral(i)
          if (dihedral(i) > MaxDihe) &
            call error_msg('Setup_Enefunc_Dihe_Pio> Too many dihedral angles.')
  
        end if

      end do
    end do
#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_dihe_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_dihe_all = found
#endif

    return

  end subroutine setup_enefunc_dihe_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_rb_dihe_pio
  !> @brief        define RB_DIHEDRAL term in potential energy function
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_rb_dihe_pio(domain, enefunc)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variables
    integer                   :: i, ix, ic, icel, found, ncell, ncell_pio
    integer                   :: file_num, file_tot_num
    integer,          pointer :: dihedral(:), dihedral_pio(:,:)
    integer,          pointer :: list(:,:,:), list_pio(:,:,:,:)
    integer,          pointer :: cell_l2g_pio(:,:)
    integer,          pointer :: dihe_pbc(:,:,:), dihe_pbc_pio(:,:,:,:)
    integer(int2),    pointer :: cell_g2l(:)
    real(wp),         pointer :: c(:,:,:), c_pio(:,:,:,:)

    cell_l2g_pio  => domain%cell_l2g_pio
    cell_g2l      => domain%cell_g2l
    dihedral      => enefunc%num_rb_dihedral
    list          => enefunc%rb_dihe_list
    c             => enefunc%rb_dihe_c
    dihe_pbc      => enefunc%rb_dihe_pbc
    dihedral_pio  => enefunc%num_rb_dihedral_pio
    list_pio      => enefunc%rb_dihe_list_pio
    c_pio         => enefunc%rb_dihe_c_pio
    dihe_pbc_pio  => enefunc%rb_dihe_pbc_pio

    ncell         = domain%num_cell_local
    ncell_pio     = domain%ncell_local_pio

    found = 0

    do file_num = 1, domain%file_tot_num

      do icel = 1, ncell_pio

        ic = cell_l2g_pio(icel,file_num)
        i  = cell_g2l(ic)

        if (i /= 0) then

          dihedral(i) = dihedral_pio(icel,file_num)
          do ix = 1, dihedral(i)
            list(1:4,ix,i) = list_pio(1:4,ix,icel,file_num)
            c(1:6,ix,i) = c_pio(1:6,ix,icel,file_num)
            dihe_pbc(1:3,ix,i) = dihe_pbc_pio(1:3,ix,icel,file_num)
          end do
          found = found + dihedral(i)
          if (dihedral(i) > MaxDihe) &
            call error_msg('Setup_Enefunc_Dihe_Pio> Too many dihedral angles.')
  
        end if

      end do
    end do
#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_rb_dihe_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_rb_dihe_all = found
#endif

    return

  end subroutine setup_enefunc_rb_dihe_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_impr_pio
  !> @brief        define IMPROPER term in potential energy function
  !! @authors      NT
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_impr_pio(domain, enefunc)

    ! formal variables
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variables
    integer                   :: i, ix, ic, icel, found, ncell, ncell_pio
    integer                   :: file_num, file_tot_num
    integer,          pointer :: cell_l2g_pio(:,:)
    integer(int2),    pointer :: cell_g2l(:)
    integer,          pointer :: improper(:), improper_pio(:,:)
    integer,          pointer :: list(:,:,:), list_pio(:,:,:,:)
    integer,          pointer :: impr_pbc(:,:,:), impr_pbc_pio(:,:,:,:)
    integer,          pointer :: period(:,:), period_pio(:,:,:)
    real(wp),         pointer :: force(:,:), force_pio(:,:,:)
    real(wp),         pointer :: phase(:,:), phase_pio(:,:,:)


    cell_l2g_pio  => domain%cell_l2g_pio
    cell_g2l      => domain%cell_g2l
    improper      => enefunc%num_improper
    list          => enefunc%impr_list
    force         => enefunc%impr_force_const
    phase         => enefunc%impr_phase
    period        => enefunc%impr_periodicity
    impr_pbc      => enefunc%impr_pbc
    improper_pio  => enefunc%num_improper_pio
    list_pio      => enefunc%impr_list_pio
    force_pio     => enefunc%impr_force_const_pio
    phase_pio     => enefunc%impr_phase_pio
    impr_pbc_pio  => enefunc%impr_pbc_pio
    period_pio    => enefunc%impr_periodicity_pio

    ncell     = domain%num_cell_local
    ncell_pio = domain%ncell_local_pio

    found = 0
    file_tot_num = domain%file_tot_num

    do file_num = 1, file_tot_num

      do icel = 1, ncell_pio

        ic = cell_l2g_pio(icel,file_num)
        i  = cell_g2l(ic)

        if (i /= 0) then

          improper(i) = improper_pio(icel,file_num)
          do ix = 1, improper(i)
            list(1:4,ix,i) = list_pio(1:4,ix,icel,file_num)
            force(ix,i) = force_pio(ix,icel,file_num)
            phase(ix,i) = phase_pio(ix,icel,file_num)
            period(ix,i) = period_pio(ix,icel,file_num)
            impr_pbc(1:3,ix,i) = impr_pbc_pio(1:3,ix,icel,file_num)
          end do
          found = found + improper(i)
          if (improper(i) > MaxImpr) &
            call error_msg( &
              'Setup_Enefunc_Impr_Pio> Too many improper dihedral angles')
 
        end if

      end do
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_impr_all, 1, mpi_integer, mpi_sum, &
                       mpi_comm_country, ierror)
#else
    enefunc%num_impr_all = found
#endif

    return

  end subroutine setup_enefunc_impr_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cmap_pio
  !> @brief        define cmap term in potential energy function with DD
  !! @authors      NT
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : energy potential functions informationn
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cmap_pio(ene_info, domain, enefunc)

    ! formal variables
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                   :: i, j, k, l, m, ix, found, ncell
    integer                   :: file_num
    integer                   :: icel, ic, ncell_pio
    integer                   :: ngrid0, ncmap_p
    integer,          pointer :: cmap(:), cmap_pio(:,:)
    integer,          pointer :: list_pio(:,:,:,:)
    integer,          pointer :: ctype_pio(:,:,:) 
    integer,          pointer :: cmap_pbc_pio(:,:,:,:) 
    integer,          pointer :: cell_l2g_pio(:,:)
    integer(int2),    pointer :: cell_g2l(:)


    cell_l2g_pio => domain%cell_l2g_pio
    cell_g2l     => domain%cell_g2l
    cmap         => enefunc%num_cmap
    cmap_pio     => enefunc%num_cmap_pio
    list_pio     => enefunc%cmap_list_pio
    ctype_pio    => enefunc%cmap_type_pio
    cmap_pbc_pio => enefunc%cmap_pbc_pio

    ncell = domain%num_cell_local
    ncell_pio = domain%ncell_local_pio

    ngrid0  = enefunc%cmap_ngrid0
    ncmap_p = enefunc%cmap_ncmap_p

    call alloc_enefunc(enefunc, EneFuncCmap, ncell, ngrid0, ncmap_p)

    do i = 1, ncmap_p
      enefunc%cmap_resolution(i) = enefunc%cmap_resolution_pio(i)
      do m = 1, ngrid0
        do l = 1, ngrid0
          do k = 1, 4
            do j = 1, 4
              enefunc%cmap_coef(j,k,l,m,i) = enefunc%cmap_coef_pio(j,k,l,m,i)
            end do
          end do
        end do
      end do
    end do
 
    found = 0
    do file_num = 1, domain%file_tot_num

      do icel = 1, ncell_pio

        ic = cell_l2g_pio(icel,file_num)
        i  = cell_g2l(ic)

        if (i /= 0) then

          cmap(i) = cmap_pio(icel,file_num)
          do ix = 1, cmap(i)
            enefunc%cmap_list(1:8,ix,i) = list_pio(1:8,ix,icel,file_num)
            enefunc%cmap_type(ix,i) = ctype_pio(ix,icel,file_num)
            enefunc%cmap_pbc (1:6,ix,i) = cmap_pbc_pio(1:6,ix,icel,file_num)
          end do
          found = found + enefunc%num_cmap(i)
          if (enefunc%num_cmap(i) > MaxCmap) &
            call error_msg('Setup_Enefunc_Cmap_Pio> Too many cmaps.')

        end if

      end do
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_cmap_all, 1, mpi_integer, mpi_sum, &
                       mpi_comm_country, ierror)
#else
    enefunc%num_cmap_all = found
#endif

    return

  end subroutine setup_enefunc_cmap_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_nonb_pio
  !> @brief        define NON-BOND term in potential energy function
  !! @authors      NT
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    constraints : constraints information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_nonb_pio(ene_info, constraints, domain, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_constraints),     intent(in)    :: constraints
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: ncel


    ! treatment for 1-2, 1-3, 1-4 interactions
    !
    ncel   = domain%num_cell_local

    call alloc_enefunc(enefunc, EneFuncNonb,     ncel, maxcell)
    call alloc_enefunc(enefunc, EneFuncNonbList, ncel, maxcell)

    if (constraints%rigid_bond) then

      call count_nonb_excl(.true., .true., constraints, domain, enefunc)

    else

      call count_nonb_excl(.true., .false., constraints, domain, enefunc)

    end if

    return

  end subroutine setup_enefunc_nonb_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_dispcorr
  !> @brief        define dispersion correction term
  !! @authors      CK
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : energy potential functions information
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine setup_enefunc_dispcorr(ene_info, domain, enefunc, constraints)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_domain),  target, intent(inout) :: domain
    type(s_enefunc), target, intent(inout) :: enefunc
    type(s_constraints),     intent(inout) :: constraints

    ! local variables
    integer                  :: i, j, iatmcls, jatmcls, ntypes
    integer                  :: icel1, icel2, icel3, icel4
    integer                  :: i1, i2, i3, i4, k, ih
    integer                  :: num_all_atoms, natom2, nexpair
    real(wip)                :: lj6_tot, lj6_diff, lj6_ex
    real(wip)                :: factor, rpair
    real(wip)                :: diff_cs, diff_cs2, diff_cs3, diff_cs4
    real(wip)                :: cutoff , cutoff2, cutoff3, cutoff4
    real(wip)                :: cutoff5, cutoff6, cutoff7, cutoff8
    real(wip)                :: cutoff14
    real(wip)                :: inv_cutoff3, inv_cutoff6, inv_cutoff12
    real(wip)                :: switchdist , switchdist2, switchdist3
    real(wip)                :: switchdist4, switchdist5
    real(wip)                :: switchdist6, switchdist7, switchdist8
    real(wip)                :: shift_a, shift_b, shift_c
    real(wip)                :: vswitch, eswitch, vlong

    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: bondlist(:,:,:),anglelist(:,:,:)
    integer,         pointer :: dihelist(:,:,:),rb_dihelist(:,:,:)
    integer,         pointer :: atmcls(:,:),imprlist(:,:,:)
    integer,     allocatable :: atype(:)


    if (ene_info%dispersion_corr == Disp_corr_NONE) return

    bondlist    => enefunc%bond_list
    anglelist   => enefunc%angle_list
    dihelist    => enefunc%dihe_list
    rb_dihelist => enefunc%rb_dihe_list
    imprlist    => enefunc%impr_list
    atmcls      => domain%atom_cls_no
    id_g2l      => domain%id_g2l

    ntypes = enefunc%num_atom_cls
    allocate(atype(1:ntypes))

    atype(1:ntypes) = 0
    num_all_atoms   = 0

    do i = 1, domain%num_cell_local
      do j = 1, domain%num_atom(i)
        iatmcls = atmcls(j,i)
        atype(iatmcls) = atype(iatmcls)+1
      end do
      num_all_atoms = num_all_atoms + domain%num_atom(i)
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, atype, ntypes, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#endif

    lj6_tot = 0.0_wip
    do i = 1, ntypes
      do j = 1, ntypes
        lj6_tot = lj6_tot + enefunc%nonb_lj6(i,j)*atype(i)*atype(j)
      end do
    end do
    deallocate(atype)

    cutoff       = enefunc%cutoffdist
    cutoff2      = cutoff*cutoff
    cutoff3      = cutoff2*cutoff
    inv_cutoff3  = 1.0_wip/cutoff3

    eswitch = 0.0_wip
    vswitch = 0.0_wip
    vlong   = inv_cutoff3/3.0_wip

    if (enefunc%forcefield == ForcefieldAMBER ) then

      factor       = 2.0_wip*PI*lj6_tot
      enefunc%dispersion_energy = -factor*vlong
      enefunc%dispersion_virial = -2.0_wip*factor*vlong

    else if (enefunc%forcefield == ForcefieldGROAMBER .or.  &
             enefunc%forcefield == ForcefieldGROMARTINI) then
      !
      ! remove exclusion
      !
      lj6_ex = 0.0_wip
      nexpair = 0
      do i = 1, domain%num_cell_local

        !constraint
        !
        do j = 1, constraints%connect
          do k = 1, constraints%HGr_local(j,i)
            iatmcls = atmcls(constraints%HGr_bond_list(1,k,j,i),i)
            do ih = 1, j
              jatmcls = atmcls(constraints%HGr_bond_list(ih+1,k,j,i),i)
              lj6_ex = lj6_ex + 2.0_wp*enefunc%nonb_lj6(iatmcls,jatmcls)
              nexpair = nexpair + 2
            end do
          end do
        end do

        do j = 1, domain%num_water(i)
          iatmcls = atmcls(domain%water_list(1,j,i),i)
          jatmcls = atmcls(domain%water_list(2,j,i),i)
          lj6_ex  = lj6_ex + 2.0_wp*enefunc%nonb_lj6(iatmcls,jatmcls)
          nexpair = nexpair + 2
          iatmcls = atmcls(domain%water_list(1,j,i),i)
          jatmcls = atmcls(domain%water_list(3,j,i),i)
          lj6_ex  = lj6_ex + 2.0_wp*enefunc%nonb_lj6(iatmcls,jatmcls)
          nexpair = nexpair + 2
          iatmcls = atmcls(domain%water_list(2,j,i),i)
          jatmcls = atmcls(domain%water_list(3,j,i),i)
          lj6_ex  = lj6_ex + 2.0_wp*enefunc%nonb_lj6(iatmcls,jatmcls)
          nexpair = nexpair + 2
        end do

        ! self
        do j = 1, domain%num_atom(i)
          iatmcls = atmcls(j,i)
          lj6_ex  = lj6_ex + enefunc%nonb_lj6(iatmcls,iatmcls)
        end do

        ! bonds
        do j = 1, enefunc%num_bond(i)
          icel1 = id_g2l(1,bondlist(1,j,i))
          i1    = id_g2l(2,bondlist(1,j,i))
          icel2 = id_g2l(1,bondlist(2,j,i))
          i2    = id_g2l(2,bondlist(2,j,i))
          lj6_ex= lj6_ex + 2.0_wp*enefunc%nb14_lj6(atmcls(i1,icel1),atmcls(i2,icel2))
        end do

        ! angles
        do j = 1, enefunc%num_angle(i)
          icel1 = id_g2l(1,anglelist(1,j,i))
          i1    = id_g2l(2,anglelist(1,j,i))
          icel3 = id_g2l(1,anglelist(3,j,i))
          i3    = id_g2l(2,anglelist(3,j,i))
          lj6_ex= lj6_ex + 2.0_wp*enefunc%nb14_lj6(atmcls(i1,icel1),atmcls(i3,icel3))
        end do

        ! dihedral
        do j = 1, enefunc%num_dihedral(i)
          icel1 = id_g2l(1,dihelist(1,j,i))
          i1    = id_g2l(2,dihelist(1,j,i))
          icel4 = id_g2l(1,dihelist(4,j,i))
          i4    = id_g2l(2,dihelist(4,j,i))
          lj6_ex= lj6_ex + 2.0_wp*enefunc%nb14_lj6(atmcls(i1,icel1),atmcls(i4,icel4))
        end do

        ! RB dihedral
        do j = 1, enefunc%num_rb_dihedral(i)
          icel1 = id_g2l(1,rb_dihelist(1,j,i))
          i1    = id_g2l(2,rb_dihelist(1,j,i))
          icel4 = id_g2l(1,rb_dihelist(4,j,i))
          i4    = id_g2l(2,rb_dihelist(4,j,i))
          lj6_ex= lj6_ex + 2.0_wp*enefunc%nb14_lj6(atmcls(i1,icel1),atmcls(i4,icel4))
        end do

        ! improper
        do j = 1, enefunc%num_improper(i)
          icel1 = id_g2l(1,imprlist(1,j,i))
          i1    = id_g2l(2,imprlist(1,j,i))
          icel4 = id_g2l(1,imprlist(4,j,i))
          i4    = id_g2l(2,imprlist(4,j,i))
          lj6_ex= lj6_ex + 2.0_wp*enefunc%nb14_lj6(atmcls(i1,icel1),atmcls(i4,icel4))
        end do

        nexpair = nexpair + domain%num_atom(i)           &
                          + 2*enefunc%num_bond(i)        &
                          + 2*enefunc%num_angle(i)       &
                          + 2*enefunc%num_dihedral(i)    &
                          + 2*enefunc%num_rb_dihedral(i) &
                          + 2*enefunc%num_improper(i)
      end do
#ifdef HAVE_MPI_GENESIS
      call mpi_allreduce(mpi_in_place, num_all_atoms, 1, mpi_integer, &
                         mpi_sum, mpi_comm_country, ierror)
      call mpi_allreduce(mpi_in_place, nexpair, 1, mpi_integer, &
                         mpi_sum, mpi_comm_country, ierror)
      call mpi_allreduce(mpi_in_place, lj6_ex, 1, mpi_wip_real, &
                         mpi_sum, mpi_comm_country, ierror)
#endif
      lj6_diff = (lj6_tot - lj6_ex)

      natom2 = num_all_atoms*num_all_atoms
      rpair  = real(natom2,wip)/real(natom2-nexpair,wip)
      factor       = 2.0_wip*PI*rpair*lj6_diff

      switchdist   = enefunc%switchdist
      diff_cs      = (cutoff - switchdist)

      if (diff_cs > EPS) then

        if (enefunc%vdw_shift) then
          cutoff4      = cutoff3*cutoff
          cutoff5      = cutoff4*cutoff
          cutoff6      = cutoff5*cutoff
          cutoff7      = cutoff6*cutoff
          cutoff8      = cutoff7*cutoff
          cutoff14     = cutoff7*cutoff7
          inv_cutoff6  = inv_cutoff3*inv_cutoff3
          inv_cutoff12 = inv_cutoff6*inv_cutoff6
  
          diff_cs2     = diff_cs*diff_cs
          diff_cs3     = diff_cs2*diff_cs
          diff_cs4     = diff_cs3*diff_cs
  
          switchdist2  = switchdist*switchdist
          switchdist3  = switchdist2*switchdist
          switchdist4  = switchdist3*switchdist
          switchdist5  = switchdist4*switchdist
          switchdist6  = switchdist5*switchdist
          switchdist7  = switchdist6*switchdist
          switchdist8  = switchdist7*switchdist
  
          ! LJ6
          !
          shift_a = -(10.0_wip*cutoff - 7.0_wip*switchdist)/(cutoff8*diff_cs2)
          shift_b =  ( 9.0_wip*cutoff - 7.0_wip*switchdist)/(cutoff8*diff_cs3)
  
          shift_c = inv_cutoff6 - 2.0_wip * shift_a * diff_cs3  &
                    - 1.5_wip * shift_b * diff_cs4
  
          eswitch = -2.0_wip * shift_a * ((1.0_wip/6.0_wip)*cutoff6            &
                                        -(3.0_wip/5.0_wip)*cutoff5*switchdist  &
                                        +(3.0_wip/4.0_wip)*cutoff4*switchdist2 &
                                        -(1.0_wip/3.0_wip)*cutoff3*switchdist3 &
                                        +(1.0_wip/6.0e1_wip)*switchdist6)      &
                    -1.5_wip * shift_b * ((1.0_wip/7.0_wip)*cutoff7            &
                                        -(2.0_wip/3.0_wip)*cutoff6*switchdist  &
                                        +(6.0_wip/5.0_wip)*cutoff5*switchdist2 &
                                        -                cutoff4*switchdist3   &
                                        +(1.0_wip/3.0_wip)*cutoff3*switchdist4 &
                                        -(1.0_wip/1.05e2_wip)*switchdist7)     &
                    -(1.0_wip/3.0_wip) * shift_c * (cutoff3)
    
          ! LJ12
          !
          shift_a = -(16.0_wip*cutoff - 13.0_wip*switchdist)/(cutoff14*diff_cs2)
          shift_b =  (15.0_wip*cutoff - 13.0_wip*switchdist)/(cutoff14*diff_cs3)
          shift_c = inv_cutoff12 - 2.0_wip * shift_a * diff_cs3  &
                    - 1.5_wip * shift_b * diff_cs4
  
 
          shift_a = -(10.0_wip*cutoff - 7.0_wip*switchdist)/(cutoff8*diff_cs2)
          shift_b =  ( 9.0_wip*cutoff - 7.0_wip*switchdist)/(cutoff8*diff_cs3)
 
          vswitch = shift_a * ( (1.0_wip/6.0_wip)*cutoff6                      &
                               -(2.0_wip/5.0_wip)*cutoff5*switchdist           &
                               +(1.0_wip/4.0_wip)*cutoff4*switchdist2          &
                               -(1.0_wip/6.0e1_wip)*switchdist6)               &
                   +shift_b * ( (1.0_wip/7.0_wip)*cutoff7                      &
                               -(1.0_wip/2.0_wip)*cutoff6*switchdist           &
                               +(3.0_wip/5.0_wip)*cutoff5*switchdist2          &
                               -(1.0_wip/4.0_wip)*cutoff4*switchdist3          &
                               +(1.0_wip/1.4e2_wip)*switchdist7)
        enefunc%dispersion_energy = factor*(eswitch-vlong)
        enefunc%dispersion_virial = -2.0_wip*factor*(-vswitch+vlong)

        else

          eswitch = enefunc%eswitch
          vswitch = enefunc%vswitch
          enefunc%dispersion_energy = factor*(eswitch-vlong)
          enefunc%dispersion_virial = -factor*(vswitch+vlong)

        end if

      else 

        enefunc%dispersion_energy = factor*(eswitch-vlong)
        enefunc%dispersion_virial = -2.0_wip*factor*(-vswitch+vlong)

      end if

    else
      call error_msg('Setup_Enefunc_DispCorr> This force field is not allowed')
    end if

!   enefunc%dispersion_energy = factor*(eswitch-vlong)
!   enefunc%dispersion_virial = -2.0_wp*factor*(-vswitch+vlong)

    return

  end subroutine setup_enefunc_dispcorr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_bonding
  !> @brief        check bonds
  !! @authors      CK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    domain   : domain information
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_bonding(enefunc, domain)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    type(s_domain),  target, intent(in)    :: domain

    ! local variables
    real(wp)                 :: d12(1:3), r12, r_dif
    integer                  :: i, j, ix, icel1, icel2, i1, i2
    integer                  :: icel3, i3,icel4, i4
    integer                  :: id, my_id, omp_get_thread_num
    integer                  :: pbc_int, k1, k2, k3, kk1, kk2, kk3
    real(wp), parameter      :: maxdistance = 0.5_wp
    real(wp)                 :: maxcell_size 

    real(wp),        pointer :: r0(:,:)
    integer,         pointer :: nbond(:), bondlist(:,:,:)
    integer,         pointer :: nangle(:), anglelist(:,:,:)
    integer,         pointer :: ndihe(:),  dihelist(:,:,:)
    integer,         pointer :: nrbdihe(:),  rb_dihelist(:,:,:)
    integer,         pointer :: nimpr(:),  imprlist(:,:,:)
    integer,         pointer :: ncell_local
    integer(int2),   pointer :: id_g2l(:,:)
    real(wip),       pointer :: coord(:,:,:)
    integer,         pointer :: bond_pbc(:,:), angl_pbc(:,:,:), dihe_pbc(:,:,:)


    ncell_local => domain%num_cell_local
    id_g2l      => domain%id_g2l
    coord       => domain%coord

    maxcell_size = max(domain%cell_size(1),  &
                       domain%cell_size(2),  &
                       domain%cell_size(3))

    nbond       => enefunc%num_bond
    bondlist    => enefunc%bond_list
    r0          => enefunc%bond_dist_min
    bond_pbc    => enefunc%bond_pbc

    nangle      => enefunc%num_angle
    anglelist   => enefunc%angle_list
    angl_pbc    => enefunc%angle_pbc

    ndihe       => enefunc%num_dihedral
    dihelist    => enefunc%dihe_list
    dihe_pbc    => enefunc%dihe_pbc

    nrbdihe     => enefunc%num_rb_dihedral
    rb_dihelist => enefunc%rb_dihe_list

    nimpr       => enefunc%num_improper
    imprlist    => enefunc%impr_list

    !$omp parallel default(shared)                                     &
    !$omp private(id, i, j, ix, icel1, i1, icel2, i2, d12, r12, r_dif, &
    !$omp         my_id, icel3, i3, icel4, i4, pbc_int, k1, k2, k3,    &
    !$omp         kk1, kk2, kk3)
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif
    my_id = id

    do i = my_id+1, ncell_local, nthread

      do ix = 1, nbond(i)

        icel1 = id_g2l(1,bondlist(1,ix,i))
        i1    = id_g2l(2,bondlist(1,ix,i))
        icel2 = id_g2l(1,bondlist(2,ix,i))
        i2    = id_g2l(2,bondlist(2,ix,i))

        pbc_int = bond_pbc(ix,i)
        k3 = pbc_int / 9
        pbc_int = pbc_int - k3*9
        k2 = pbc_int / 3
        k1 = pbc_int - k2*3
        k1 = k1 - 1
        k2 = k2 - 1
        k3 = k3 - 1

        ! bond energy: E=K[b-b0]^2
        !
        d12(1) = coord(1,i1,icel1) - coord(1,i2,icel2) &
               + domain%system_size(1)*real(k1,wp)
        d12(2) = coord(2,i1,icel1) - coord(2,i2,icel2) &
               + domain%system_size(2)*real(k2,wp)
        d12(3) = coord(3,i1,icel1) - coord(3,i2,icel2) &
               + domain%system_size(3)*real(k3,wp)
        r12   = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )
        r_dif = r12 - r0(ix,i)
        if (r_dif > maxdistance) &
           write(MsgOut,'(A,I10,I10,F10.5)') &
          'WARNING: too long bond:',bondlist(1,ix,i),bondlist(2,ix,i),r12
        if (r_dif < -maxdistance) &
           write(MsgOut,'(A,I10,I10,F10.5)') &
          'WARNING: too short bond:',bondlist(1,ix,i),bondlist(2,ix,i),r12
        if (r12 > maxcell_size) then
           write(MsgOut,'(A,2I10,F10.5)') &
          'Check_bonding> distance is grater than cellsize:', &
           bondlist(1,ix,i),bondlist(2,ix,i),r12
           call error_msg('')
        end if

      end do

      do ix = 1, nangle(i)

        icel1 = id_g2l(1,anglelist(1,ix,i))
        i1    = id_g2l(2,anglelist(1,ix,i))
        icel3 = id_g2l(1,anglelist(3,ix,i))
        i3    = id_g2l(2,anglelist(3,ix,i))

        pbc_int = angl_pbc(3,ix,i)
        k3 = pbc_int / 9
        pbc_int = pbc_int - k3*9
        k2 = pbc_int / 3
        k1 = pbc_int - k2*3
        k1 = k1 - 1
        k2 = k2 - 1
        k3 = k3 - 1
        d12(1) = coord(1,i1,icel1) - coord(1,i3,icel3) &
               + domain%system_size(1)*real(k1,wp)
        d12(2) = coord(2,i1,icel1) - coord(2,i3,icel3) &
               + domain%system_size(2)*real(k2,wp)
        d12(3) = coord(3,i1,icel1) - coord(3,i3,icel3) &
               + domain%system_size(3)*real(k3,wp)
        r12   = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )

        if (r12 > maxcell_size) then
           write(MsgOut,'(A,2I10,F10.5)') &
           'Check_bonding> distance in angle is grater than cellsize:', &
           anglelist(1,ix,i),anglelist(3,ix,i),r12
           call error_msg('')
        end if

      end do

      do ix = 1, ndihe(i)

        icel1 = id_g2l(1,dihelist(1,ix,i))
        i1    = id_g2l(2,dihelist(1,ix,i))
        icel4 = id_g2l(1,dihelist(4,ix,i))
        i4    = id_g2l(2,dihelist(4,ix,i))

        pbc_int = dihe_pbc(1,ix,i)
        k3 = pbc_int / 9
        pbc_int = pbc_int - k3*9
        k2 = pbc_int / 3
        k1 = pbc_int - k2*3
        kk1 = k1 - 1
        kk2 = k2 - 1
        kk3 = k3 - 1

        pbc_int = dihe_pbc(2,ix,i)
        k3 = pbc_int / 9
        pbc_int = pbc_int - k3*9
        k2 = pbc_int / 3
        k1 = pbc_int - k2*3
        kk1 = kk1 + k1 - 1
        kk2 = kk2 + k2 - 1
        kk3 = kk3 + k3 - 1

        pbc_int = dihe_pbc(3,ix,i)
        k3 = pbc_int / 9
        pbc_int = pbc_int - k3*9
        k2 = pbc_int / 3
        k1 = pbc_int - k2*3
        kk1 = kk1 - k1 + 1
        kk2 = kk2 - k2 + 1
        kk3 = kk3 - k3 + 1

        d12(1) = coord(1,i1,icel1) - coord(1,i4,icel4) &
               + domain%system_size(1)*real(kk1,wp)
        d12(2) = coord(2,i1,icel1) - coord(2,i4,icel4) &
               + domain%system_size(2)*real(kk2,wp)
        d12(3) = coord(3,i1,icel1) - coord(3,i4,icel4) &
               + domain%system_size(3)*real(kk3,wp)

        r12   = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )

        if (r12 > maxcell_size) then
           write(MsgOut,'(A,2I10,F10.5)') &
           'Check_bonding> distance in dihedral is grater than cellsize:', &
           dihelist(1,ix,i),dihelist(4,ix,i),r12
           call error_msg('')
        end if

      end do

      if (nrbdihe(i) > 0) then
        do ix = 1, nrbdihe(i)
       
          icel1 = id_g2l(1,rb_dihelist(1,ix,i))
          i1    = id_g2l(2,rb_dihelist(1,ix,i))
          icel4 = id_g2l(1,rb_dihelist(4,ix,i))
          i4    = id_g2l(2,rb_dihelist(4,ix,i))
       
          d12(1:3) = coord(1:3,i1,icel1) - coord(1:3,i4,icel4)
          r12   = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )
          if (r12 > maxcell_size) then
            write(MsgOut,'(A,2I10,F10.5)') &
           'Check_bonding> distance in rb dihedral is grater than cellsize:', &
            rb_dihelist(1,ix,i), rb_dihelist(4,ix,i),r12
            call error_msg('')
          end if
       
        end do
      end if

      do ix = 1, nimpr(i)

        icel1 = id_g2l(1,imprlist(1,ix,i))
        i1    = id_g2l(2,imprlist(1,ix,i))
        icel4 = id_g2l(1,imprlist(4,ix,i))
        i4    = id_g2l(2,imprlist(4,ix,i))

        d12(1:3) = coord(1:3,i1,icel1) - coord(1:3,i4,icel4)
        r12   = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )

        if (r12 > maxcell_size) then
          write(MsgOut,'(A,2I10,F10.5)') &
      'Check_bonding> distance in improper dihedral is grater than cellsize:', &
          imprlist(1,ix,i), imprlist(4,ix,i),r12
          call error_msg('')
        end if
      end do

    end do

    !$omp end parallel 

    return

  end subroutine check_bonding

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_enefunc_fep
  !> @brief        a driver subroutine for updating potential energy functions
  !                for FEP
  !! @authors      HO
  !! @param[inout] domain      : domain information
  !! @param[inout] comm        : communication information
  !! @param[inout] enefunc     : energy potential functions information
  !! @param[inout] constraints : constraints information [optional]
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_enefunc_fep(domain, comm, enefunc, constraints)

    ! formal arguments
    type(s_domain),                intent(inout) :: domain
    type(s_comm),                  intent(inout) :: comm
    type(s_enefunc),               intent(inout) :: enefunc
    type(s_constraints), optional, intent(inout) :: constraints

    ! local variables
    logical                        :: first


    ! sending the bonding information to other domain
    !

    ! bond
    !
    call update_outgoing_enefunc_bond_fep(domain, enefunc)

    ! enm
    !
    if (enefunc%enm_use) then
      call update_enefunc_enm(domain, enefunc)
      call update_cell_size_enm(domain, enefunc, comm)
      call update_cell_boundary_enm(domain, enefunc, comm)
    end if

    ! angle
    !
    call update_outgoing_enefunc_angl_fep(domain, enefunc)

    ! dihedral
    !
    call update_outgoing_enefunc_dihe_fep(domain, enefunc)

    ! Ryckaert-Bellemans dihedral
    !
    call update_outgoing_enefunc_rb_dihe_fep(domain, enefunc)

    ! improper dihedral
    !
    call update_outgoing_enefunc_impr_fep(domain, enefunc)

    ! cmap
    !
    call update_outgoing_enefunc_cmap_fep(domain, enefunc)

    ! restraint
    !
    call update_outgoing_enefunc_restraint(domain, enefunc)

    ! fitting
    !
    call update_outgoing_enefunc_fitting(domain, enefunc)

    ! communicate neighbour domain
    !
    call communicate_bond_fep(domain, comm, enefunc)

    ! bond
    !
    call update_incoming_enefunc_bond_fep(domain, enefunc)

    ! angle
    !
    call update_incoming_enefunc_angl_fep(domain, enefunc)

    ! dihedral
    !
    call update_incoming_enefunc_dihe_fep(domain, enefunc)

    ! Ryckaert-Bellemans dihedral
    !
    call update_incoming_enefunc_rb_dihe_fep(domain, enefunc)

    ! improper dihedral
    !
    call update_incoming_enefunc_impr_fep(domain, enefunc)

    ! cmap
    !
    call update_incoming_enefunc_cmap_fep(domain, enefunc)

    ! restraint
    !
    call update_incoming_enefunc_restraint(domain, enefunc)

    ! fitting
    !
    call update_incoming_enefunc_fitting(domain, enefunc)

    ! re-count nonbond exclusion list
    !

    first = .false.

    if (constraints%rigid_bond) then

      call count_nonb_excl_fep(first, .true., constraints, domain, enefunc)

    else

      call count_nonb_excl_fep(first, .false., constraints, domain, enefunc)

    end if

    if (enefunc%bonding_check) call check_bonding(enefunc, domain)

    return

  end subroutine update_enefunc_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_dispcorr_fep
  !> @brief        define dispersion correction term for FEP
  !! @authors      HO
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : energy potential functions information
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_dispcorr_fep(ene_info, domain, enefunc, constraints)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_domain),  target, intent(inout) :: domain
    type(s_enefunc), target, intent(inout) :: enefunc
    type(s_constraints),     intent(inout) :: constraints

    ! local variables
    integer                  :: i, j, iatmcls, jatmcls, ntypes
    integer                  :: icel1, icel2, icel3, icel4
    integer                  :: i1, i2, i3, i4, k, ih
    real(wip)                :: diff_cs, diff_cs2, diff_cs3, diff_cs4
    real(wip)                :: cutoff , cutoff2, cutoff3, cutoff4
    real(wip)                :: cutoff5, cutoff6, cutoff7, cutoff8
    real(wip)                :: cutoff14
    real(wip)                :: inv_cutoff3, inv_cutoff6, inv_cutoff12
    real(wip)                :: switchdist , switchdist2, switchdist3
    real(wip)                :: switchdist4, switchdist5
    real(wip)                :: switchdist6, switchdist7, switchdist8
    real(wip)                :: shift_a, shift_b, shift_c
    real(wip)                :: vswitch, eswitch, vlong

    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: bondlist(:,:,:),anglelist(:,:,:)
    integer,         pointer :: dihelist(:,:,:),rb_dihelist(:,:,:)
    integer,         pointer :: atmcls(:,:),imprlist(:,:,:)

    ! FEP
    integer                  :: num_all_atoms_C
    integer                  :: num_all_atoms_A
    integer                  :: num_all_atoms_B
    integer                  :: fg1, fg2
    real(wip)                :: lj6_tot_CC, lj6_tot_AC, lj6_tot_BC
    real(wip)                :: lj6_tot_AA, lj6_tot_BB, lj6_tot_AB
    real(wip)                :: factor_CC, factor_AC, factor_BC
    real(wip)                :: factor_AA, factor_BB, factor_AB
    integer,     allocatable :: atype_C(:), atype_A(:), atype_B(:)

    if (ene_info%dispersion_corr == Disp_corr_NONE) return

    bondlist    => enefunc%bond_list
    anglelist   => enefunc%angle_list
    dihelist    => enefunc%dihe_list
    rb_dihelist => enefunc%rb_dihe_list
    imprlist    => enefunc%impr_list
    atmcls      => domain%atom_cls_no
    id_g2l      => domain%id_g2l

    ntypes = enefunc%num_atom_cls
    allocate(atype_C(1:ntypes), atype_A(1:ntypes), atype_B(1:ntypes))

    atype_C(1:ntypes) = 0
    atype_A(1:ntypes) = 0
    atype_B(1:ntypes) = 0
    num_all_atoms_C   = 0
    num_all_atoms_A   = 0
    num_all_atoms_B   = 0

    do i = 1, domain%num_cell_local
      do j = 1, domain%num_atom(i)
        iatmcls = atmcls(j,i)
        if(enefunc%fep_topology == 2) then
          ! Dual
          if (domain%fepgrp(j,i) == 5) then
            atype_C(iatmcls) = atype_C(iatmcls)+1
            num_all_atoms_C = num_all_atoms_C + 1
          else if (domain%fepgrp(j,i) == 3) then
            atype_A(iatmcls) = atype_A(iatmcls)+1
            num_all_atoms_A = num_all_atoms_A + 1
          else if (domain%fepgrp(j,i) == 4) then
            atype_B(iatmcls) = atype_B(iatmcls)+1
            num_all_atoms_B = num_all_atoms_B + 1
          end if
        else
          ! Hybrid
          if (domain%fepgrp(j,i) == 5) then
            atype_C(iatmcls) = atype_C(iatmcls)+1
            num_all_atoms_C = num_all_atoms_C + 1
          else if (domain%fepgrp(j,i) == 1) then
            atype_A(iatmcls) = atype_A(iatmcls)+1
            atype_B(iatmcls) = atype_B(iatmcls)+1
            num_all_atoms_A = num_all_atoms_A + 1
            num_all_atoms_B = num_all_atoms_B + 1
          else if (domain%fepgrp(j,i) == 3) then
            atype_A(iatmcls) = atype_A(iatmcls)+1
            num_all_atoms_A = num_all_atoms_A + 1
          else if (domain%fepgrp(j,i) == 4) then
            atype_B(iatmcls) = atype_B(iatmcls)+1
            num_all_atoms_B = num_all_atoms_B + 1
          end if
        end if
      end do
    end do


#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, atype_C, ntypes, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
    call mpi_allreduce(mpi_in_place, atype_A, ntypes, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
    call mpi_allreduce(mpi_in_place, atype_B, ntypes, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#endif

    lj6_tot_CC = 0.0_wip
    lj6_tot_AC = 0.0_wip
    lj6_tot_BC = 0.0_wip
    lj6_tot_AA = 0.0_wip
    lj6_tot_BB = 0.0_wip
    lj6_tot_AB = 0.0_wip
    do i = 1, ntypes
      do j = 1, ntypes
        lj6_tot_CC = lj6_tot_CC + enefunc%nonb_lj6(i,j)*atype_C(i)*atype_C(j)
        lj6_tot_AC = lj6_tot_AC + enefunc%nonb_lj6(i,j)*atype_A(i)*atype_C(j)
        lj6_tot_AC = lj6_tot_AC + enefunc%nonb_lj6(i,j)*atype_C(i)*atype_A(j)
        lj6_tot_AA = lj6_tot_AA + enefunc%nonb_lj6(i,j)*atype_A(i)*atype_A(j)
        lj6_tot_BC = lj6_tot_BC + enefunc%nonb_lj6(i,j)*atype_B(i)*atype_C(j)
        lj6_tot_BC = lj6_tot_BC + enefunc%nonb_lj6(i,j)*atype_C(i)*atype_B(j)
        lj6_tot_BB = lj6_tot_BB + enefunc%nonb_lj6(i,j)*atype_B(i)*atype_B(j)
        lj6_tot_AB = lj6_tot_AB + enefunc%nonb_lj6(i,j)*atype_A(i)*atype_B(j)
        lj6_tot_AB = lj6_tot_AB + enefunc%nonb_lj6(i,j)*atype_B(i)*atype_A(j)
      end do
    end do

    deallocate(atype_C, atype_A, atype_B)

    cutoff       = enefunc%cutoffdist
    cutoff2      = cutoff*cutoff
    cutoff3      = cutoff2*cutoff
    inv_cutoff3  = 1.0_wip/cutoff3

    eswitch = 0.0_wip
    vswitch = 0.0_wip
    vlong   = inv_cutoff3/3.0_wip

    if (enefunc%forcefield == ForcefieldAMBER ) then

      factor_CC = 2.0_wip*PI*lj6_tot_CC
      factor_AC = 2.0_wip*PI*lj6_tot_AC
      factor_BC = 2.0_wip*PI*lj6_tot_BC
      factor_AA = 2.0_wip*PI*lj6_tot_AA
      factor_BB = 2.0_wip*PI*lj6_tot_BB
      factor_AB = 2.0_wip*PI*lj6_tot_AB
      enefunc%dispersion_energy_CC = -factor_CC*vlong
      enefunc%dispersion_energy_AC = -factor_AC*vlong
      enefunc%dispersion_energy_BC = -factor_BC*vlong
      enefunc%dispersion_energy_AA = -factor_AA*vlong
      enefunc%dispersion_energy_BB = -factor_BB*vlong
      enefunc%dispersion_energy_AB = -factor_AB*vlong
      enefunc%dispersion_virial_CC = -2.0_wip*factor_CC*vlong
      enefunc%dispersion_virial_AC = -2.0_wip*factor_AC*vlong
      enefunc%dispersion_virial_BC = -2.0_wip*factor_BC*vlong
      enefunc%dispersion_virial_AA = -2.0_wip*factor_AA*vlong
      enefunc%dispersion_virial_BB = -2.0_wip*factor_BB*vlong
      enefunc%dispersion_virial_AB = -2.0_wip*factor_AB*vlong

    else
      call error_msg('Setup_Enefunc_DispCorr_Fep> This force field is not available in FEP')
    end if

    return

  end subroutine setup_enefunc_dispcorr_fep

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine set_fepgrp_nonb(enefunc)

    ! formal arguments
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: fg1, fg2
    integer                  :: i1, i2, i3, i4
    integer                  :: i5, i6, i7, i8
    integer                  :: idx
    integer,         pointer :: fepgrp_nonb(:,:)

    fepgrp_nonb       => enefunc%fepgrp_nonb

    ! Make lambda table for calculation of nonbonded energy
    !
    fepgrp_nonb(:,:) = 0

    do fg1 = 1, 5
      do fg2 = 1, 5

        ! Make table of lambda and softcore in FEP
        !
        if (((fg1==1).or.(fg1==2).or.(fg1==5)) .and. &
          ((fg2==1).or.(fg2==2).or.(fg2==5))) then
          fepgrp_nonb(fg1,fg2) = 5
        else if ((fg1==3).or.(fg2==3)) then
          fepgrp_nonb(fg1,fg2) = 3
        else if ((fg1==4).or.(fg2==4)) then
          fepgrp_nonb(fg1,fg2) = 4
        end if

        ! Interactons between partA and partB are set to 0,
        ! while others are set to 1.
        !
        if (((fg1==3).and.(fg2==4)) .or. ((fg1==4).and.(fg2==3))) then
          fepgrp_nonb(fg1,fg2) = 0
        end if

      end do
    end do

    return

  end subroutine set_fepgrp_nonb

end module sp_enefunc_mod
