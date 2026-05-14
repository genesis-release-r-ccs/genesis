!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_enefunc_gromacs_mod
!> @brief   Defines potential energy functions.
!! @authors Jaewoon Jung (JJ), Takeshi Imai (TI), Chigusa Kobayashi (CK),
!!          Takaharu Mori (TM), Yuji Sugita (YS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_enefunc_gromacs_mod

  use cg_enefunc_restraints_mod
  use cg_energy_mod
  use cg_restraints_str_mod
  use cg_enefunc_str_mod
  use cg_energy_str_mod
  use cg_domain_str_mod
  use cg_communicate_str_mod
  use molecules_str_mod
  use table_libs_mod
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
  public  :: define_enefunc_gromacs_lb
  public  :: count_nonb_excl_go
  public  :: count_nonb_excl
  private :: setup_enefunc_bond
  private :: setup_enefunc_angl
  private :: setup_enefunc_dihe
  private :: setup_enefunc_cgDNA_nonb
  private :: setup_enefunc_cg_basetype
  private :: setup_enefunc_nonb
  private :: setup_enefunc_contact

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc_gromacs
  !> @brief        A driver subroutine for defining potential energy functions.
  !! @authors      NT
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    grotop      : GROMACS TOP information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    restraints  : restraints information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc_gromacs(ene_info, grotop, molecule, &
                                    restraints, domain, enefunc, comm)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
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

    ! bond
    !
    call setup_enefunc_bond_alloc(grotop, domain, enefunc)
    call alloc_enefunc(enefunc, EneFuncBond, MaxBond, 1)
    call setup_enefunc_bond(grotop, domain, enefunc)

    ! angle
    !
    call setup_enefunc_angl_alloc(grotop, domain, enefunc)
    call alloc_enefunc(enefunc, EneFuncAngl, MaxAngl, 1)
    call setup_enefunc_angl(grotop, domain, enefunc)

    ! dihedral
    !
    call setup_enefunc_dihe_alloc(grotop, domain, enefunc)
    call alloc_enefunc(enefunc, EneFuncDihe, MaxDihe, 1)
    call setup_enefunc_dihe(grotop, domain, enefunc)

    if (enefunc%forcefield == ForcefieldGroMartini) then
      MaxExcl = MaxBond+MaxAngl+MaxDihe
      call alloc_enefunc(enefunc, EneFuncExcl, MaxExcl, 1)
    end if

    if (enefunc%forcefield == ForcefieldRESIDCG) then

      !CG protein-protein KH model
      !
      call setup_enefunc_cg_KH(grotop, molecule, domain, enefunc)

      ! CG IDR: KH model
      !
      call setup_enefunc_cg_IDR_KH(grotop, molecule, domain, enefunc)

      ! CG IDR: HPS model
      !
      enefunc%cg_IDR_HPS_epsilon = ene_info%cg_IDR_HPS_epsilon
      call setup_enefunc_cg_IDR_HPS(grotop, molecule, domain, enefunc)

      ! CG IDR: MPIPI model
      !
      enefunc%cg_WangFrenkel_cutoff_scale = ene_info%cg_IDR_MPIPI_cutoff_scale
      call setup_enefunc_cg_IDR_MPIPI(grotop, molecule, domain, enefunc)

      ! base stacking
      !
      enefunc%cg_infinite_DNA  = ene_info%cg_infinite_DNA
      call setup_enefunc_cgDNA_base_stack_alloc(grotop, molecule, domain, &
                                                enefunc)
      call setup_enefunc_cgDNA_base_stack(grotop, molecule, domain, enefunc)

    end if

    ! nonbonded
    !
    call setup_enefunc_nonb(ene_info, grotop, molecule, domain, enefunc, comm)

    if (enefunc%forcefield == ForcefieldRESIDCG) then

      ! short range non-bonded
      !
      call setup_enefunc_cgDNA_nonb(ene_info, grotop, enefunc)

      ! CG ele: Debye-Huckel
      !
      call setup_enefunc_cg_ele(ene_info, grotop, molecule, domain, enefunc)

      ! CG PWMcos: protein-DNA interaction
      !
      call setup_enefunc_cg_PWMcos(ene_info, grotop, molecule, domain, enefunc)

      ! CG PWMcosns: protein-DNA interaction
      !
      enefunc%pwmcosns_sigma   = ene_info%cg_PWMcosns_sigma
      enefunc%pwmcosns_phi     = ene_info%cg_PWMcosns_phi * RAD
      call setup_enefunc_cg_PWMcosns(grotop, molecule, domain, enefunc)

    end if

    ! restraints
    !
    call setup_enefunc_restraints(molecule, restraints, domain, enefunc)


    ! write summary of energy function
    !
    if (main_rank) then
      write(MsgOut,'(A)') &
           'Define_Enefunc_Gromacs> Number of interactions in each term.'
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           '  bond_ene        = ', enefunc%num_bonds,           &
           '  angle_ene       = ', enefunc%num_angles+enefunc%num_anglocal
      if (enefunc%num_bonds_quartic > 0)  then
         write(MsgOut,'(A20,I10)')                         &
              '  bond_ene_cgDNA  = ', enefunc%num_bonds_quartic
      endif
      if (enefunc%num_angflex > 0)  then
        write(MsgOut,'(A20,I10)')                         &
           '  flex_angle_ene  = ', enefunc%num_angflex
      endif
      write(MsgOut,'(A20,I10)')                          &
           '  torsion_ene     = ', enefunc%num_dihe_all
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           ' restraint_groups = ', enefunc%num_restraintgroups, &
           ' restraint_funcs  = ', enefunc%num_restraintfuncs
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine define_enefunc_gromacs

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc_gromacs_lb
  !> @brief        Defines potential energy functions after load balancing.
  !! @authors      JJ
  !! @param[in]    grotop      : GROMACS TOP information
  !! @param[in]    molecule    : molecule information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc_gromacs_lb(grotop, molecule, domain, enefunc, comm)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_comm),            intent(inout) :: comm    

    ! local variables
    integer                  :: ncel, ncelb
    integer                  :: i, j, k, l

    ! base
    !
    ncel  = domain%num_cell_local
    ncelb = domain%num_cell_local + domain%num_cell_boundary

    ! bond
    !
    call setup_enefunc_bond_alloc(grotop, domain, enefunc)
    call alloc_enefunc(enefunc, EneFuncBond, MaxBond, 1)
    call setup_enefunc_bond(grotop, domain, enefunc)

    ! angle
    !
    call setup_enefunc_angl_alloc(grotop, domain, enefunc)
    call alloc_enefunc(enefunc, EneFuncAngl, MaxAngl, 1)
    call setup_enefunc_angl(grotop, domain, enefunc)

    ! dihedral
    !
    call setup_enefunc_dihe_alloc(grotop, domain, enefunc)
    call alloc_enefunc(enefunc, EneFuncDihe, MaxDihe, 1)
    call setup_enefunc_dihe(grotop, domain, enefunc)

    if (enefunc%forcefield == ForcefieldRESIDCG) then

      !CG protein-protein KH model
      !
      call setup_enefunc_cg_KH(grotop, molecule, domain, enefunc)

      ! CG IDR: KH model
      !
      call setup_enefunc_cg_IDR_KH(grotop, molecule, domain, enefunc)

      ! CG IDR: HPS model
      !
      call setup_enefunc_cg_IDR_HPS(grotop, molecule, domain, enefunc)

      ! CG IDR: MPIPI model
      !
      call setup_enefunc_cg_IDR_MPIPI(grotop, molecule, domain, enefunc)

      ! base stacking
      !
      call setup_enefunc_cgDNA_base_stack_alloc(grotop, molecule, domain, &
                                                enefunc)
      call setup_enefunc_cgDNA_base_stack(grotop, molecule, domain, enefunc)

    end if

    ! nonbonded
    !
    call setup_enefunc_nonb_lb(grotop, molecule, domain, enefunc, comm)

    if (enefunc%forcefield == ForcefieldRESIDCG) then

      ! CG PWMcos: protein-DNA interaction
      !
      call setup_enefunc_cg_PWMcos_lb(grotop, molecule, domain, enefunc)

      ! CG PWMcosns: protein-DNA interaction
      !
      call setup_enefunc_cg_PWMcosns(grotop, molecule, domain, enefunc)

    end if

    return

  end subroutine define_enefunc_gromacs_lb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_bond_alloc
  !> @brief        Defines BOND term size.
  !! @authors      NT
  !! @param[in]    grotop   : GROMACS TOP information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond_alloc(grotop, domain, enefunc)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: ioffset, nbond_s, nbond_q
    integer                  :: idx1, idx2, icel1, icel2, icel_local

    type(s_grotop_mol), pointer :: gromol
    integer,            pointer :: ncel
    integer,            pointer :: id_g2l(:)
    integer,            pointer :: atom_2_cell(:)
    real(wp),           pointer :: natom(:)

    ncel        => domain%num_cell_local
    id_g2l      => domain%id_g2l
    atom_2_cell => domain%atom_2_cell
    natom       => domain%num_atom_t0

    nbond_s      =  0
    nbond_q      =  0

    ! first we make a bond list of square term
    !
    ioffset   = 0
    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        do k = 1, gromol%num_bonds

          if (gromol%bonds(k)%func /= 21) then

            idx1 = gromol%bonds(k)%atom_idx1 + ioffset
            idx2 = gromol%bonds(k)%atom_idx2 + ioffset
            icel1 = id_g2l(idx1)
            icel2 = id_g2l(idx2)
  
            if (icel1 /= 0 .and. icel2 /= 0) then

              icel1 = atom_2_cell(icel1)
              icel2 = atom_2_cell(icel2)
              if (icel1 <= ncel .and. icel2 <= ncel) then
                icel_local = min(icel1,icel2)
              else
                if (natom(icel1) >= natom(icel2)) then
                  icel_local = icel2
                else if (natom(icel2) > natom(icel1)) then
                  icel_local = icel1
                end if
              end if
              if (icel_local > 0 .and. icel_local <= ncel) &
                nbond_s = nbond_s + 1
        
            end if

          else if (gromol%bonds(k)%func == 21) then

            idx1 = gromol%bonds(k)%atom_idx1 + ioffset
            idx2 = gromol%bonds(k)%atom_idx2 + ioffset
            icel1 = id_g2l(idx1)
            icel2 = id_g2l(idx2)

            if (icel1 /= 0 .and. icel2 /= 0) then

              icel1 = atom_2_cell(icel1)
              icel2 = atom_2_cell(icel2)
              if (icel1 <= ncel .and. icel2 <= ncel) then
                icel_local = min(icel1,icel2)
              else
                if (natom(icel1) >= natom(icel2)) then
                  icel_local = icel2
                else if (natom(icel2) > natom(icel1)) then
                  icel_local = icel1
                end if
              end if

              if (icel_local > 0 .and. icel_local <= ncel) &
                nbond_q = nbond_q + 1

            end if

          end if

        end do

        ioffset = ioffset + gromol%num_atoms

      end do
    end do

    enefunc%num_bondsq_domain = nbond_s
    k = nbond_s
    call mpi_allreduce(k, BondMoveS, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
    BondMoveS = BondMoveS / 2
    enefunc%num_bond_domain   = nbond_q
    k = nbond_q
    call mpi_allreduce(k, BondMoveQ, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
    BondMoveQ = BondMoveQ / 2
    k = nbond_s + nbond_q
    call mpi_allreduce(k, MaxBond, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
    if (MaxBond > 0) MaxBond = max(MaxBond*2, min_atom)

    return

  end subroutine setup_enefunc_bond_alloc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_bond
  !> @brief        Defines BOND term for each cell in potential energy function
  !! @authors      NT
  !! @param[in]    grotop   : GROMACS TOP information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond(grotop, domain, enefunc)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k, b_idx
    integer                  :: ioffset, nbond_s, nbond_q
    integer                  :: idx1, idx2, icel1, icel2, icel_local

    type(s_grotop_mol), pointer :: gromol
    real(wp),           pointer :: force(:), dist(:)
    integer,            pointer :: list(:,:)
    integer,            pointer :: ncel
    integer,            pointer :: id_g2l(:)
    integer,            pointer :: atom_2_cell(:)
    real(wp),           pointer :: natom(:)

    ncel          => domain%num_cell_local
    id_g2l        => domain%id_g2l
    atom_2_cell   => domain%atom_2_cell
    natom         => domain%num_atom_t0

    list          => enefunc%bond_list
    force         => enefunc%bond_force_const
    dist          => enefunc%bond_dist_min

    nbond_s = 0
    nbond_q = 0

    ! first we make a bond list of square term
    !
    ioffset   = 0
    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count

        do k = 1, gromol%num_bonds

          if (gromol%bonds(k)%func /= 21) then

            idx1 = gromol%bonds(k)%atom_idx1 + ioffset
            idx2 = gromol%bonds(k)%atom_idx2 + ioffset
            icel1 = id_g2l(idx1)
            icel2 = id_g2l(idx2)

            if (icel1 /= 0 .and. icel2 /= 0) then

              icel1 = atom_2_cell(icel1)
              icel2 = atom_2_cell(icel2)
              if (icel1 <= ncel .and. icel2 <= ncel) then
                icel_local = min(icel1,icel2)
              else
                if (natom(icel1) >= natom(icel2)) then
                  icel_local = icel2
                else if (natom(icel2) > natom(icel1)) then
                  icel_local = icel1
                end if
              end if

              if (icel_local > 0 .and. icel_local <= ncel) then
                nbond_s = nbond_s + 1
                b_idx = nbond_s
                list (1,b_idx) = idx1
                list (2,b_idx) = idx2
                list (3,b_idx) = icel_local
                if (gromol%bonds(k)%func == 1) then
                  force(b_idx) = gromol%bonds(k)%kb*0.01_wp*JOU2CAL*0.5_wp
                else
                  force(b_idx) = gromol%bonds(k)%kb*0.0001_wp*JOU2CAL*0.25_wp
                end if
                dist(b_idx) = gromol%bonds(k)%b0 * 10.0_wp
              end if

            end if
          end if
        end do

        ioffset = ioffset + gromol%num_atoms

      end do
    end do

    ! second we make a bond list of quartic term
    !
    ioffset   = 0
    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count

        do k = 1, gromol%num_bonds

          if (gromol%bonds(k)%func == 21) then

            idx1 = gromol%bonds(k)%atom_idx1 + ioffset
            idx2 = gromol%bonds(k)%atom_idx2 + ioffset
            icel1 = id_g2l(idx1)
            icel2 = id_g2l(idx2)

            if (icel1 /= 0 .and. icel2 /= 0) then

              icel1 = atom_2_cell(icel1)
              icel2 = atom_2_cell(icel2)
              if (icel1 <= ncel .and. icel2 <= ncel) then
                icel_local = min(icel1,icel2)
              else
                if (natom(icel1) >= natom(icel2)) then
                  icel_local = icel2
                else if (natom(icel2) > natom(icel1)) then
                  icel_local = icel1
                end if
              end if

              if (icel_local > 0 .and. icel_local <= ncel) then

                nbond_q = nbond_q + 1
                b_idx = nbond_s + nbond_q
                list (1,b_idx) = idx1
                list (2,b_idx) = idx2
                list (3,b_idx) = icel_local

                force(b_idx) = gromol%bonds(k)%kb*0.01_wp*JOU2CAL*0.5_wp
                dist (b_idx) = gromol%bonds(k)%b0*10.0_wp
              end if
            end if
          end if
        end do

        ioffset = ioffset + gromol%num_atoms

      end do
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(nbond_s, enefunc%num_bonds, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
    call mpi_allreduce(nbond_q, enefunc%num_bonds_quartic, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_bonds = nbond_s
    enefunc%num_bonds_quartic = nbond_q
#endif

    return

  end subroutine setup_enefunc_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_angl_alloc
  !> @brief        Defines ANGLE size.
  !! @authors      NT
  !! @param[in]    grotop   : GROMACS topology informaiton
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_angl_alloc(grotop, domain, enefunc)

    ! formal arguments
    type(s_grotop),  target, intent(in)    :: grotop
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                    :: i, j, k, nangl_f, nangl_l, nangl
    integer                    :: ioffset
    integer                    :: idx1, idx2, idx3, icel1, icel2, icel_local

    type(s_grotop_mol), pointer :: gromol
    integer,            pointer :: ncel
    integer,            pointer :: id_g2l(:)
    integer,            pointer :: atom_2_cell(:)
    real(wp),           pointer :: natom(:)

    ncel           => domain%num_cell_local
    id_g2l         => domain%id_g2l
    atom_2_cell    => domain%atom_2_cell
    natom          => domain%num_atom_t0

    nangl_f = 0
    nangl_l = 0
    nangl   = 0
    ioffset = 0

    do i = 1, grotop%num_molss

      gromol => grotop%molss(i)%moltype%mol

      do j = 1, grotop%molss(i)%count

        do k = 1, gromol%num_angls

          if (gromol%angls(k)%func == 22 .or. &
              gromol%angls(k)%func == 52) then

            idx1 = gromol%angls(k)%atom_idx1 + ioffset
            idx2 = gromol%angls(k)%atom_idx2 + ioffset
            idx3 = gromol%angls(k)%atom_idx3 + ioffset
            icel1 = id_g2l(idx1)
            icel2 = id_g2l(idx3)

            if (icel1 /= 0 .and. icel2 /= 0) then

              icel1 = atom_2_cell(icel1)
              icel2 = atom_2_cell(icel2)
              if (icel1 <= ncel .and. icel2 <= ncel) then
                icel_local = min(icel1,icel2)
              else 
                if (natom(icel1) >= natom(icel2)) then
                  icel_local = icel2
                else if (natom(icel2) > natom(icel1)) then
                  icel_local = icel1
                end if
              end if
  
              if (icel_local > 0 .and. icel_local <= ncel) &
                nangl_f = nangl_f + 1

            end if

          else if (gromol%angls(k)%func == 21) then

            idx1 = gromol%angls(k)%atom_idx1 + ioffset
            idx2 = gromol%angls(k)%atom_idx2 + ioffset
            idx3 = gromol%angls(k)%atom_idx3 + ioffset
            icel1 = id_g2l(idx1)
            icel2 = id_g2l(idx3)

            if (icel1 /= 0 .and. icel2 /= 0) then

              icel1 = atom_2_cell(icel1)
              icel2 = atom_2_cell(icel2)
              if (icel1 <= ncel .and. icel2 <= ncel) then
                icel_local = min(icel1,icel2)
              else
                if (natom(icel1) >= natom(icel2)) then
                  icel_local = icel2
                else if (natom(icel2) > natom(icel1)) then
                  icel_local = icel1
                end if
              end if

              if (icel_local > 0 .and. icel_local <= ncel) &
                nangl_l = nangl_l + 1

            end if

          else if (gromol%angls(k)%func == 1) then

            idx1 = gromol%angls(k)%atom_idx1 + ioffset
            idx2 = gromol%angls(k)%atom_idx2 + ioffset
            idx3 = gromol%angls(k)%atom_idx3 + ioffset
            icel1 = id_g2l(idx1)
            icel2 = id_g2l(idx3)

            if (icel1 /= 0 .and. icel2 /= 0) then

              icel1 = atom_2_cell(icel1)
              icel2 = atom_2_cell(icel2)
              if (icel1 <= ncel .and. icel2 <= ncel) then
                icel_local = min(icel1,icel2)
              else
                if (natom(icel1) >= natom(icel2)) then
                  icel_local = icel2
                else if (natom(icel2) > natom(icel1)) then
                  icel_local = icel1
                end if
              end if

              if (icel_local > 0 .and. icel_local <= ncel) &
                nangl = nangl + 1

            end if

          else if (gromol%angls(k)%func == 2) then

            idx1 = gromol%angls(k)%atom_idx1 + ioffset
            idx2 = gromol%angls(k)%atom_idx2 + ioffset
            idx3 = gromol%angls(k)%atom_idx3 + ioffset
            icel1 = id_g2l(idx1)
            icel2 = id_g2l(idx3)

            if (icel1 /= 0 .and. icel2 /= 0) then

              icel1 = atom_2_cell(icel1)
              icel2 = atom_2_cell(icel2)
              if (icel1 <= ncel .and. icel2 <= ncel) then
                icel_local = min(icel1,icel2)
              else
                if (natom(icel1) >= natom(icel2)) then
                  icel_local = icel2
                else if (natom(icel2) > natom(icel1)) then
                  icel_local = icel1
                end if
              end if

              if (icel_local > 0 .and. icel_local <= ncel) &
                nangl = nangl + 1

            end if

          end if

        end do

        ioffset = ioffset + gromol%num_atoms

      end do
    end do

    enefunc%num_angle_flexible_domain = nangl_f
    k = nangl_f
    call mpi_allreduce(k, AnglMoveF, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
    AnglMoveF = AnglMoveF / 2

    enefunc%num_angle_local_domain = nangl_l
    k = nangl_l
    call mpi_allreduce(k, AnglMoveL, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
    AnglMoveL = AnglMoveL / 2

    enefunc%num_angle_domain = nangl
    k = nangl
    call mpi_allreduce(k, AnglMove, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
    AnglMoveL = AnglMove / 2

    k = enefunc%num_angle_flexible_domain + enefunc%num_angle_local_domain &
      + enefunc%num_angle_domain 

    call mpi_allreduce(k, MaxAngl, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
    if (min_atom > 0) MaxAngl = max(MaxAngl*2, min_atom)

    return

  end subroutine setup_enefunc_angl_alloc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_angl
  !> @brief        Defines ANGLE term for each cell in potential energy function
  !! @authors      NT
  !! @param[in]    grotop   : GROMACS topology informaiton
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_angl(grotop, domain, enefunc)

    ! formal arguments
    type(s_grotop),  target, intent(in)    :: grotop
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    real(wp)                    :: min_th, max_th, min_th_ener, max_th_ener
    real(wp)                    :: min_energy, etmp, center, t123, gradient
    integer                     :: i, j, k, a_idx, nangl_f, nangl_l, nangl
    integer                     :: ioffset
    integer                     :: idx1, idx2, idx3, icel1, icel2, icel_local
    integer                     :: nanglflextypes, ntable

    type(s_grotop_mol), pointer :: gromol
    real(wp),           pointer :: force(:), theta(:), width(:)
    integer,            pointer :: list(:,:), atype(:)
    integer,            pointer :: ncel
    integer,            pointer :: id_g2l(:)
    integer,            pointer :: atom_2_cell(:)
    real(wp),           pointer :: natom(:)

    ncel                 => domain%num_cell_local
    id_g2l               => domain%id_g2l
    atom_2_cell          => domain%atom_2_cell
    natom                => domain%num_atom_t0

    list                 => enefunc%angl_list
    force                => enefunc%angl_force_const
    theta                => enefunc%angl_theta_min
    width                => enefunc%angl_width
    atype                => enefunc%angl_kind 

    nangl_f = 0
    nangl_l = 0
    nangl   = 0

    ! table for flexible angle
    !
    nanglflextypes = grotop%num_flangltypes

    if (nanglflextypes > 0) then

      ntable = size(grotop%flangltypes(1)%theta(:))
      call alloc_enefunc(enefunc, EneFuncAngFlexTbl, nanglflextypes, ntable)
      center = ( AICG2P_FBA_MAX_ANG - AICG2P_FBA_MIN_ANG ) * 0.5_wp

      do i = 1, nanglflextypes

        enefunc%anglflex_theta (1:ntable, i) = &
                                grotop%flangltypes(i)%theta(1:ntable)
        enefunc%anglflex_efunc (1:ntable, i) = &
                                grotop%flangltypes(i)%efunc(1:ntable) * JOU2CAL
        enefunc%anglflex_d2func(1:ntable, i) = &
                                grotop%flangltypes(i)%d2func(1:ntable) * JOU2CAL

        t123   = AICG2P_FBA_MIN_ANG
        min_th = AICG2P_FBA_MIN_ANG
        max_th = AICG2P_FBA_MIN_ANG

        call table_flexibleangle(i, t123, enefunc%anglflex_theta,   &
                                 enefunc%anglflex_efunc,            &
                                 enefunc%anglflex_d2func,           &
                                 etmp, gradient)
        min_th_ener = etmp
        min_energy  = etmp

        do while(t123 <= AICG2P_FBA_MAX_ANG)

          t123 = t123 + AICG2P_FBA_DTHEATA
          call table_flexibleangle(i, t123, enefunc%anglflex_theta,   &
                                   enefunc%anglflex_efunc,             &
                                   enefunc%anglflex_d2func,            &
                                   etmp, gradient)
          min_energy = min(min_energy, etmp)
          if (gradient < AICG2P_FBA_MIN_ANG_FORCE) then
            min_th      = t123
            min_th_ener =  etmp
          endif
          if (t123 > center .and.                        &
              abs(max_th-AICG2P_FBA_MIN_ANG) < EPS .and. &
              gradient > AICG2P_FBA_MAX_ANG_FORCE) then
            max_th      = t123
            max_th_ener =  etmp
          endif

        end do

        enefunc%anglflex_min_th(1, i) = min_th
        enefunc%anglflex_max_th(1, i) = max_th
        enefunc%anglflex_min_th(2, i) = min_th_ener
        enefunc%anglflex_max_th(2, i) = max_th_ener
        enefunc%anglflex_ener_corr(i) = min_energy

      end do

    endif

    ! flexible angle setup
    !
    ioffset   = 0
    do i = 1, grotop%num_molss

      gromol => grotop%molss(i)%moltype%mol

      do j = 1, grotop%molss(i)%count

        if (gromol%settles%func == 0) then

          do k = 1, gromol%num_angls

            if (gromol%angls(k)%func == 22 .or. &
                gromol%angls(k)%func == 52) then

              idx1 = gromol%angls(k)%atom_idx1 + ioffset
              idx2 = gromol%angls(k)%atom_idx2 + ioffset
              idx3 = gromol%angls(k)%atom_idx3 + ioffset

              icel1 = id_g2l(idx1)
              icel2 = id_g2l(idx3)

              if (icel1 /= 0 .and. icel2 /= 0) then

                icel1 = atom_2_cell(icel1)
                icel2 = atom_2_cell(icel2)
                if (icel1 <= ncel .and. icel2 <= ncel) then
                  icel_local = min(icel1,icel2)
                else
                  if (natom(icel1) >= natom(icel2)) then
                    icel_local = icel2
                  else if (natom(icel2) > natom(icel1)) then
                    icel_local = icel1
                  end if
                end if

                if (icel_local > 0 .and. icel_local <= ncel) then

                  nangl_f = nangl_f + 1
                  a_idx = nangl_f
                  list (1:3,a_idx) = (/idx1, idx2, idx3/)
                  list (4,a_idx) = icel_local
                  atype (a_idx) = gromol%angls(k)%types

                end if
              end if
            end if
          end do

        end if

        ioffset = ioffset + gromol%num_atoms

      end do
    end do

    if (enefunc%num_angflex > 0 .and. nanglflextypes <= 0 .and. main_rank) &
      call error_msg('Error: Setup_Enefunc_Angl - Flexible angle type is '//&
                     'not defined.')

    ! setup for angle_local (gaussian type angle energy)
    !
    ioffset = 0
    do i = 1, grotop%num_molss

      gromol => grotop%molss(i)%moltype%mol

      do j = 1, grotop%molss(i)%count

        if (gromol%settles%func == 0) then

          do k = 1, gromol%num_angls

            if (gromol%angls(k)%func == 21) then

              idx1 = gromol%angls(k)%atom_idx1 + ioffset
              idx2 = gromol%angls(k)%atom_idx2 + ioffset
              idx3 = gromol%angls(k)%atom_idx3 + ioffset

              icel1 = id_g2l(idx1)
              icel2 = id_g2l(idx3)

              if (icel1 /= 0 .and. icel2 /= 0) then

                icel1 = atom_2_cell(icel1)
                icel2 = atom_2_cell(icel2)
                if (icel1 <= ncel .and. icel2 <= ncel) then
                  icel_local = min(icel1,icel2)
                else
                  if (natom(icel1) >= natom(icel2)) then
                    icel_local = icel2
                  else if (natom(icel2) > natom(icel1)) then
                    icel_local = icel1
                  end if
                end if

                if (icel_local > 0 .and. icel_local <= ncel) then

                  nangl_l = nangl_l + 1
                  a_idx = nangl_l + nangl_f
                  list (1:3,a_idx) = (/idx1, idx2, idx3/)
                  list (4,a_idx) = icel_local
                  force(    a_idx) = -gromol%angls(k)%kt * JOU2CAL
                  theta(    a_idx) = gromol%angls(k)%theta_0 * 10.0_wp
                  width(    a_idx) = gromol%angls(k)%w * 10.0_wp
                end if

              end if
            end if
          end do

        end if

        ioffset = ioffset + gromol%num_atoms

      end do
    end do

    ! normal angle term
    !
    ioffset = 0
    do i = 1, grotop%num_molss

      gromol => grotop%molss(i)%moltype%mol

      do j = 1, grotop%molss(i)%count

        if (gromol%settles%func == 0) then

          do k = 1, gromol%num_angls

            if (gromol%angls(k)%func == 1) then

              idx1 = gromol%angls(k)%atom_idx1 + ioffset
              idx2 = gromol%angls(k)%atom_idx2 + ioffset
              idx3 = gromol%angls(k)%atom_idx3 + ioffset

              icel1 = id_g2l(idx1)
              icel2 = id_g2l(idx3)

              if (icel1 /= 0 .and. icel2 /= 0) then

                icel1 = atom_2_cell(icel1)
                icel2 = atom_2_cell(icel2)
                if (icel1 <= ncel .and. icel2 <= ncel) then
                  icel_local = min(icel1,icel2)
                else
                  if (natom(icel1) >= natom(icel2)) then
                    icel_local = icel2
                  else if (natom(icel2) > natom(icel1)) then
                    icel_local = icel1
                  end if
                end if

                if (icel_local > 0 .and. icel_local <= ncel) then

                  if (nangl > MaxAngl) &
                  call error_msg('Error: Setup_Enefunc_Angl - Too many angles.')

                  nangl = nangl + 1
                  a_idx = nangl + nangl_f + nangl_l
                  list (1:3,a_idx) = (/idx1, idx2, idx3/)
                  list (4,a_idx) = icel_local
                  force(    a_idx) = gromol%angls(k)%kt * JOU2CAL * 0.5_wp
                  theta(    a_idx) = gromol%angls(k)%theta_0 * RAD
                end if

              end if

            else if (gromol%angls(k)%func == 2) then

              idx1 = gromol%angls(k)%atom_idx1 + ioffset
              idx2 = gromol%angls(k)%atom_idx2 + ioffset
              idx3 = gromol%angls(k)%atom_idx3 + ioffset
              icel1 = id_g2l(idx1)
              icel2 = id_g2l(idx3)

              if (icel1 /= 0 .and. icel2 /= 0) then

                icel1 = atom_2_cell(icel1)
                icel2 = atom_2_cell(icel2)
                if (icel1 <= ncel .and. icel2 <= ncel) then
                  icel_local = min(icel1,icel2)
                else
                  if (natom(icel1) >= natom(icel2)) then
                    icel_local = icel2
                  else if (natom(icel2) > natom(icel1)) then
                    icel_local = icel1
                  end if
                end if

                if (icel_local > 0 .and. icel_local <= ncel) then

                  if (nangl > MaxAngl) &
                  call error_msg('Error: Setup_Enefunc_Angl - Too many angles.')

                  nangl = nangl + 1
                  list (1:3,nangl) = (/idx1, idx2, idx3/)
                  list (4  ,nangl) = icel_local
                  force(    nangl) = gromol%angls(k)%kt * JOU2CAL * 0.5_wp
                  theta(    nangl) = gromol%angls(k)%theta_0 * RAD
                end if

              end if

            end if
          end do
        end if

        ioffset = ioffset + gromol%num_atoms

      end do
    end do

#ifdef HAVE_MPI_GENESIS
    nangl   = enefunc%num_angle_domain
    nangl_f = enefunc%num_angle_flexible_domain
    nangl_l = enefunc%num_angle_local_domain

    call mpi_allreduce(nangl, enefunc%num_angles, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
    call mpi_allreduce(nangl_f, enefunc%num_angflex, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
    call mpi_allreduce(nangl_l, enefunc%num_anglocal, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#endif

    return

  end subroutine setup_enefunc_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_dihe_alloc
  !> @brief        Defines DIHEDRAL size.
  !! @authors      NT
  !! @param[in]    grotop   : GROMACS TOP information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_dihe_alloc(grotop, domain, enefunc)

    ! formal arguments
    type(s_grotop),  target, intent(in)    :: grotop
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k, l
    integer                  :: ioffset
    integer                  :: idx(4), icel1, icel2, icel
    integer                  :: num_dihe, num_flex, num_local

    type(s_grotop_mol), pointer :: gromol
    integer,            pointer :: ncel
    integer,            pointer :: id_g2l(:)
    integer,            pointer :: atom_2_cell(:)
    real(wp),           pointer :: natom(:)

    ncel          => domain%num_cell_local
    id_g2l        => domain%id_g2l
    atom_2_cell   => domain%atom_2_cell
    natom         => domain%num_atom_t0

    ioffset   = 0
    num_dihe  = 0
    num_flex  = 0
    num_local = 0
    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count

        do k = 1, gromol%num_dihes

          if (gromol%dihes(k)%func == 22 .or. &
              gromol%dihes(k)%func == 52) then

            idx(1) = gromol%dihes(k)%atom_idx1 + ioffset
            idx(2) = gromol%dihes(k)%atom_idx2 + ioffset
            idx(3) = gromol%dihes(k)%atom_idx3 + ioffset
            idx(4) = gromol%dihes(k)%atom_idx4 + ioffset
  
            icel1 = id_g2l(idx(1))
            icel2 = id_g2l(idx(4))

            if (icel1 /= 0 .and. icel2 /= 0) then
  
              icel1 = atom_2_cell(icel1)
              icel2 = atom_2_cell(icel2)
              if (icel1 <= ncel .and. icel2 <= ncel) then
                icel = min(icel1,icel2)
              else
                if (natom(icel1) >= natom(icel2)) then
                  icel = icel2
                else if (natom(icel2) > natom(icel1)) then
                  icel = icel1
                end if
              end if

              if (icel > 0 .and. icel <= ncel) then
                num_flex = num_flex + 1
              end if
            end if

          else if (gromol%dihes(k)%func == 21 .or. &
                   gromol%dihes(k)%func == 41) then

            idx(1) = gromol%dihes(k)%atom_idx1 + ioffset
            idx(2) = gromol%dihes(k)%atom_idx2 + ioffset
            idx(3) = gromol%dihes(k)%atom_idx3 + ioffset
            idx(4) = gromol%dihes(k)%atom_idx4 + ioffset

            icel1 = id_g2l(idx(1))
            icel2 = id_g2l(idx(4))

            if (icel1 /= 0 .and. icel2 /= 0) then

              icel1 = atom_2_cell(icel1)
              icel2 = atom_2_cell(icel2)
              if (icel1 <= ncel .and. icel2 <= ncel) then
                icel = min(icel1,icel2)
              else
                if (natom(icel1) >= natom(icel2)) then
                  icel = icel2
                else if (natom(icel2) > natom(icel1)) then
                  icel = icel1
                end if
              end if

              if (icel > 0 .and. icel <= ncel) then
                num_local = num_local + 1
              end if
            end if

          else if (gromol%dihes(k)%func == 1  .or. &
                   gromol%dihes(k)%func == 31 .or. &
                   gromol%dihes(k)%func == 32) then

            idx(1) = gromol%dihes(k)%atom_idx1 + ioffset
            idx(2) = gromol%dihes(k)%atom_idx2 + ioffset
            idx(3) = gromol%dihes(k)%atom_idx3 + ioffset
            idx(4) = gromol%dihes(k)%atom_idx4 + ioffset

            icel1 = id_g2l(idx(1))
            icel2 = id_g2l(idx(4))

            if (icel1 /= 0 .and. icel2 /= 0) then

              icel1 = atom_2_cell(icel1)
              icel2 = atom_2_cell(icel2)
              if (icel1 <= ncel .and. icel2 <= ncel) then
                icel = min(icel1,icel2)
              else
                if (natom(icel1) >= natom(icel2)) then
                  icel = icel2
                else if (natom(icel2) > natom(icel1)) then
                  icel = icel1
                end if
              end if

              if (icel > 0 .and. icel <= ncel) then
                num_dihe = num_dihe + 1
              end if
            end if

          end if

        end do
        ioffset = ioffset + gromol%num_atoms
      end do
    end do
 
    enefunc%num_dihe_flexible_domain = num_flex
    call mpi_allreduce(l, DiheMoveF, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)

    enefunc%num_dihe_local_domain = num_local
    call mpi_allreduce(l, DiheMoveL, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)

    enefunc%num_dihe_domain = num_dihe
    k = enefunc%num_dihe_flexible_domain + enefunc%num_dihe_local_domain &
      + enefunc%num_dihe_domain

    call mpi_allreduce(k, MaxDihe, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
    if (MaxDihe > 0) MaxDihe = max(MaxDihe*2, min_atom)

    return

  end subroutine setup_enefunc_dihe_alloc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_dihe
  !> @brief        Defines DIHEDRAL term in potential energy function
  !! @authors      NT
  !! @param[in]    grotop   : GROMACS TOP information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_dihe(grotop, domain, enefunc)

    ! formal arguments
    type(s_grotop),  target, intent(in)    :: grotop
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: ioffset
    integer                  :: idx(4), icel1, icel2, icel, d_idx
    integer                  :: ndiheflex, ndiheflextypes, ntable
    integer                  :: ndihe_local, ndihe
    integer                  :: num_flex, num_local, num_dihe
    real(wp)                 :: th, cos_dih, sin_dih, min_energy, etmp

    type(s_grotop_mol), pointer :: gromol
    real(wp),           pointer :: force(:), phase(:), width(:)
    integer,            pointer :: list(:,:), dtype(:), period(:)
    integer,            pointer :: ncel
    integer,            pointer :: id_g2l(:)
    integer,            pointer :: atom_2_cell(:)
    integer,            pointer :: notation
    real(wp),           pointer :: natom(:)

    ncel             => domain%num_cell_local
    id_g2l           => domain%id_g2l
    atom_2_cell      => domain%atom_2_cell
    natom            => domain%num_atom_t0

    list             => enefunc%dihe_list
    force            => enefunc%dihe_force_const
    phase            => enefunc%dihe_phase
    width            => enefunc%dihe_width
    dtype            => enefunc%dihe_kind
    period           => enefunc%dihe_periodicity
    notation         => enefunc%notation_14types
    notation         = 100

    ndiheflex = enefunc%num_dihe_flexible_domain
    ndihe_local = enefunc%num_dihe_local_domain
    ndihe = enefunc%num_dihe_domain

    ! parameter for each flexible dihderal angle type
    !
    ndiheflextypes = grotop%num_fldihetypes

    ! parameter for each flexible dihedral angle type
    !
    if (ndiheflextypes > 0) then
      ntable = size(grotop%fldihetypes(1)%coef(:))
      call alloc_enefunc(enefunc, EneFuncDiheFlexTbl, ndiheflextypes, ntable)
      do i = 1, ndiheflextypes
        enefunc%diheflex_coef(1:ntable,i) = &
            grotop%fldihetypes(i)%coef(1:ntable)*JOU2CAL
        th = -PI
        cos_dih = cos(th)
        sin_dih = sin(th)
        min_energy = enefunc%diheflex_coef(1,i)                      &
                   + enefunc%diheflex_coef(2,i)*cos_dih              &
                   + enefunc%diheflex_coef(3,i)*sin_dih              &
                   + enefunc%diheflex_coef(4,i)                      &
                    *(2.0_wp*cos_dih*cos_dih-1.0_wp)                 &
                   + enefunc%diheflex_coef(5,i)                      &
                    *2.0_wp*cos_dih*sin_dih                          &
                   + enefunc%diheflex_coef(6,i)                      &
                    *(4.0_wp*cos_dih*cos_dih*cos_dih-3.0_wp*cos_dih) &
                   + enefunc%diheflex_coef(7,i)                      &
                    *(-4.0_wp*sin_dih*sin_dih*sin_dih+3.0_wp*sin_dih)
        do while (th <= PI)
          th = th + AICG2P_FBA_DTHEATA
          cos_dih = cos(th)
          sin_dih = sin(th)
          etmp = enefunc%diheflex_coef(1,i) &
                     + enefunc%diheflex_coef(2,i)*cos_dih              &
                     + enefunc%diheflex_coef(3,i)*sin_dih              &
                     + enefunc%diheflex_coef(4,i)                      &
                      *(2.0_wp*cos_dih*cos_dih-1.0_wp)                 &
                     + enefunc%diheflex_coef(5,i)                      &
                      *2.0_wp*cos_dih*sin_dih                          &
                     + enefunc%diheflex_coef(6,i)                      &
                      *(4.0_wp*cos_dih*cos_dih*cos_dih-3.0_wp*cos_dih) &
                     + enefunc%diheflex_coef(7,i)                      &
                      *(-4.0_wp*sin_dih*sin_dih*sin_dih+3.0_wp*sin_dih)
          min_energy = min(min_energy, etmp)
        end do
        enefunc%diheflex_ener_corr(i) = min_energy
      end do
    end if

    ! flexible dihedral angle setup
    !
    ioffset   = 0
    num_flex  = 0 
    num_local = 0
    num_dihe  = 0
    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count

        do k = 1, gromol%num_dihes

          if (gromol%dihes(k)%func == 32 .or. &
              gromol%dihes(k)%func == 41 .or. &
              gromol%dihes(k)%func == 52) &
             enefunc%cg_safe_dihedral_calc   = .true.

          idx(1) = gromol%dihes(k)%atom_idx1 + ioffset
          idx(2) = gromol%dihes(k)%atom_idx2 + ioffset
          idx(3) = gromol%dihes(k)%atom_idx3 + ioffset
          idx(4) = gromol%dihes(k)%atom_idx4 + ioffset

          icel1 = id_g2l(idx(1))
          icel2 = id_g2l(idx(4))

          if (icel1 /= 0 .and. icel2 /= 0) then

            icel1 = atom_2_cell(icel1)
            icel2 = atom_2_cell(icel2)
            if (icel1 <= ncel .and. icel2 <= ncel) then
              icel = min(icel1,icel2)
            else
              if (natom(icel1) >= natom(icel2)) then
                icel = icel2
              else if (natom(icel2) > natom(icel1)) then
                icel = icel1
              end if
            end if

            if (icel > 0 .and. icel <= ncel) then

              if (gromol%dihes(k)%func == 22 .or. &
                  gromol%dihes(k)%func == 52) then
                num_flex = num_flex + 1
                d_idx = num_flex
                list (1:4,d_idx) = idx(1:4)
                list (5,d_idx) = icel
                dtype(    d_idx) = gromol%dihes(k)%types
                period(   d_idx) = gromol%dihes(k)%func
              else if (gromol%dihes(k)%func == 21 .or. &
                       gromol%dihes(k)%func == 41) then
                num_local = num_local + 1
                d_idx = ndiheflex + num_local
                list (1:4,d_idx) = idx(1:4)
                list (5,d_idx) = icel
                force(    d_idx) = - gromol%dihes(k)%kp * JOU2CAL
                phase(    d_idx) = gromol%dihes(k)%theta_0 * RAD
                width(    d_idx) = gromol%dihes(k)%w
                dtype(    d_idx) = gromol%dihes(k)%func
              else if (gromol%dihes(k)%func == 1  .or. &
                       gromol%dihes(k)%func == 31 .or. &
                       gromol%dihes(k)%func == 32) then
                num_dihe = num_dihe + 1
                d_idx = ndiheflex + ndihe_local + num_dihe
                list (1:4,d_idx) = idx(1:4)
                list (5,d_idx) = icel
                force (   d_idx) = gromol%dihes(k)%kp * JOU2CAL
                phase (   d_idx) = gromol%dihes(k)%ps * RAD
                period(   d_idx) = gromol%dihes(k)%multiplicity
                dtype(    d_idx) = gromol%dihes(k)%func
              end if
            end if
          end if
        end do
        ioffset = ioffset + gromol%num_atoms
      end do
    end do

#ifdef HAVE_MPI_GENESIS
    enefunc%num_dihe_flexible_domain = ndiheflex
    enefunc%num_dihe_local_domain = ndihe_local
    enefunc%num_dihe_domain = ndihe
    call mpi_allreduce(ndihe+ndiheflex+ndihe_local, enefunc%num_dihe_all, 1, &
                       mpi_integer, mpi_sum, mpi_comm_country, ierror)
#endif
 
    return

  end subroutine setup_enefunc_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cg_basetype
  !> @brief        Defines parameters in CG potential energy function
  !! @authors      JJ
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[in]    molecule : molecule including molecular information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cg_basetype(grotop, molecule, atom_cls_2_base_type, &
                                       enefunc)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),target, intent(in)    :: molecule
    integer,                 intent(inout) :: atom_cls_2_base_type(:)
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer              :: i, l, k, i_B6, i_B5, natom
    logical              :: chain_check
    integer,     pointer :: chain_id(:), DNA_end(:,:)

    chain_id => molecule%molecule_no
    natom    = molecule%num_atoms
    DNA_end  => enefunc%DNA_end
!   integer, allocatable :: atom_cls_2_base_type(:)

    ! ------------
    ! NA_base_type
    ! ------------
    !
    do i = 1, grotop%num_atomtypes
      if ( grotop%atomtypes(i)%type_name == 'DA' ) then
        atom_cls_2_base_type(i) = NABaseTypeDBA
      else if ( grotop%atomtypes(i)%type_name == 'DC' ) then
        atom_cls_2_base_type(i) = NABaseTypeDBC
      else if ( grotop%atomtypes(i)%type_name == 'DG' ) then
        atom_cls_2_base_type(i) = NABaseTypeDBG
      else if ( grotop%atomtypes(i)%type_name == 'DT' ) then
        atom_cls_2_base_type(i) = NABaseTypeDBT
      else if ( grotop%atomtypes(i)%type_name == 'DP' ) then
        atom_cls_2_base_type(i) = NABaseTypeDP
      else if ( grotop%atomtypes(i)%type_name == 'DS' ) then
        atom_cls_2_base_type(i) = NABaseTypeDS
      else if ( grotop%atomtypes(i)%type_name == 'RA' ) then
        atom_cls_2_base_type(i) = NABaseTypeRBA
      else if ( grotop%atomtypes(i)%type_name == 'RC' ) then
        atom_cls_2_base_type(i) = NABaseTypeRBC
      else if ( grotop%atomtypes(i)%type_name == 'RG' ) then
        atom_cls_2_base_type(i) = NABaseTypeRBG
      else if ( grotop%atomtypes(i)%type_name == 'RU' ) then
        atom_cls_2_base_type(i) = NABaseTypeRBU
      else if ( grotop%atomtypes(i)%type_name == 'RP' ) then
        atom_cls_2_base_type(i) = NABaseTypeRP
      else if ( grotop%atomtypes(i)%type_name == 'RS' ) then
        atom_cls_2_base_type(i) = NABaseTypeRS
      else
        atom_cls_2_base_type(i) = NABaseTypeProtein
      end if
    end do

    do i = 1, natom
      l = molecule%atom_cls_no(i)
      k = atom_cls_2_base_type(l)
      if (k <= NABaseTypeDBMAX) then
        i_B6 = i + 3
        i_B5 = i - 3
        chain_check = .false.
        if (i_B6 > natom) then
          i_B6 = i - 3
          chain_check = .true.
        else if (chain_id(i_B6) /= chain_id(i)) then
          i_B6 = i - 3
          chain_check = .true.
        end if
        if (chain_check) then
          do while(.true.)
            i_B6 = i_B6 - 3
            if (i_B6 <= 0) exit
            if (chain_id(i_B6) /= chain_id(i)) exit
          end do
          i_B6 = i_B6 + 3
          DNA_end(1,chain_id(i)) = i_B6
        end if
        chain_check = .false.
        if ( i_B5 <= 0) then
          i_B5 = i + 3
          chain_check = .true.
        else if (chain_id(i_B5) /= chain_id(i)) then
          i_B5 = i + 3
          chain_check = .true.
        end if
        if (chain_check) then
          do while(.true.)
            i_B5 = i_B5 + 3
            if (i_B5 > natom) exit
            if (chain_id(i_B5) /= chain_id(i)) exit
          end do 
          i_B5 = i_B5 - 3
          DNA_end(2,chain_id(i)) = i_B5
        end if
      end if
    end do
        
    return

  end subroutine setup_enefunc_cg_basetype

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cg_KH
  !> @brief        Defines parameters for CG protein KH model.
  !! @authors      CT, JJ
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========2=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cg_KH(grotop, molecule, domain, enefunc)

    ! formal arguments
    type(s_grotop),  target, intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_domain),  target, intent(inout) :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                    :: i, j, k, ncel, n_mols, ig, num_kh

    ncel = domain%num_cell_local + domain%num_cell_boundary

    if (grotop%num_cgkhmolpairs > 0) enefunc%cg_KH_calc = .true.

    n_mols = molecule%num_molecules
    call alloc_enefunc(enefunc, EneFuncCGKHmol, n_mols)

    do i = 1, grotop%num_cgkhmolpairs

      if (grotop%cg_KH_mol_pairs(i)%is_intermol) then
        do j = grotop%cg_KH_mol_pairs(i)%grp1_start, &
               grotop%cg_KH_mol_pairs(i)%grp1_end
          do k = grotop%cg_KH_mol_pairs(i)%grp2_start, &
                 grotop%cg_KH_mol_pairs(i)%grp2_end
            enefunc%cg_KH_mol_pair(j,k) = int(grotop%cg_KH_mol_pairs(i)%func,kind=1)
            enefunc%cg_KH_mol_pair(k,j) = int(grotop%cg_KH_mol_pairs(i)%func,kind=1)
          end do
        end do
      else
        do j = grotop%cg_KH_mol_pairs(i)%grp1_start, &
               grotop%cg_KH_mol_pairs(i)%grp1_end
          enefunc%cg_KH_mol_pair(j,j) = int(grotop%cg_KH_mol_pairs(i)%func,kind=1)
        end do
      end if

    end do

    num_kh = 0
    do i = 1, domain%num_atom_domain+domain%num_atom_boundary
      ig = domain%id_l2g(i)
      do j = 1, grotop%num_cg_KH_atomtypes
        if (grotop%cg_KH_atomtypes(j)%type_name &
            == molecule%atom_cls_name(ig)) then
          num_kh = num_kh + 1
          domain%cg_pro_use_KH(i) = 1
          domain%NA_base_type(i) = NABaseTypeKH
          exit
        end if
      end do
    end do
    enefunc%num_cg_KH = num_kh
#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, num_kh, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
#endif
    if (num_kh > 0) then
      Max_cg_KH = max(num_kh*2,min_atom)
    else
      Max_cg_KH = 0
    end if
    num_kh = Max_cg_KH

    call alloc_enefunc(enefunc, EneFuncCGKHList, num_kh)
    call alloc_enefunc(enefunc, EneFuncCGKHInvList, MaxAtom_domain)

    num_kh = 0
    do i = 1, domain%num_atom_domain+domain%num_atom_boundary
      ig = domain%id_l2g(i)
      do j = 1, grotop%num_cg_KH_atomtypes
        if (grotop%cg_KH_atomtypes(j)%type_name &
            == molecule%atom_cls_name(ig)) then
          num_kh = num_kh + 1
          enefunc%cg_KH_list(num_kh) = i
          enefunc%cg_KH_list_inv(i) = num_kh
          exit
        end if
      end do
    end do

    return

  end subroutine setup_enefunc_cg_KH

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cg_IDR_KH
  !> @brief        Defines parameters for CG protein KH model.
  !! @authors      CT, JJ
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cg_IDR_KH(grotop, molecule, domain, enefunc)

    ! formal arguments
    type(s_grotop),  target, intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_domain),  target, intent(inout) :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    ! local variables
    integer              :: n_atoms, n_idr_region, num_kh
    integer              :: ioffset
    integer              :: i, j, k, l, ix, ig, ncel

    type(s_grotop_mol), pointer :: gromol

    ! ------------------------------------
    ! count IDR KH "regions" in itp files
    ! ------------------------------------
    !
    n_idr_region = 0
    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        do k = 1, gromol%num_idr_kh
          n_idr_region = n_idr_region + 1
        end do
      end do
    end do

    ! switch on
    if ( n_idr_region > 0 ) &
      enefunc%cg_IDR_KH_calc = .true.

    ! -------------------------------------
    ! Set IDR regions if specified in input
    ! -------------------------------------
    !
    n_atoms = 0
    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        ioffset = n_atoms
        n_atoms = n_atoms + gromol%num_atoms
        do k = 1, gromol%num_idr_kh
          do l = gromol%idr_kh(k)%grp_start, gromol%idr_kh(k)%grp_end
            ix = domain%id_g2l(l+ioffset)
            if (ix > 0 .and. &
                ix <= domain%num_atom_domain+domain%num_atom_boundary) &
              domain%cg_IDR_KH(ix) = 1
          end do
        end do
      end do
    end do

    num_kh = 0
    do i = 1, domain%num_atom_domain+domain%num_atom_boundary
      if (domain%cg_IDR_KH(i) == 0) cycle
      ig = domain%id_l2g(i)
      do j = 1, grotop%num_cg_KH_atomtypes
        if (grotop%cg_KH_atomtypes(j)%type_name &
            == molecule%atom_cls_name(ig)) then
          num_kh = num_kh + 1
          domain%charge(i) = grotop%cg_KH_atomtypes(j)%charge
          if (domain%cg_pro_use_KH(i) == 1) then
            domain%NA_base_type(i) = NABaseTypeBothKH
          else
            domain%NA_base_type(i) = NABaseTypeIDRKH
          end if
          exit
        end if
      end do
    end do
    enefunc%num_cg_IDR_KH = num_kh
#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, num_kh, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
#endif
    if (num_kh > 0) then
      Max_cg_IDR_KH = max(num_kh*2,min_atom)
    else
      Max_cg_IDR_KH = 0
    end if
    num_kh = Max_cg_IDR_KH

    call alloc_enefunc(enefunc, EneFuncCGIDRKHList, num_kh)
    call alloc_enefunc(enefunc, EneFuncCGIDRKHInvList, MaxAtom_domain)

    num_kh = 0
    do i = 1, domain%num_atom_domain+domain%num_atom_boundary
      if (domain%cg_IDR_KH(i) == 0) cycle
      ig = domain%id_l2g(i)
      do j = 1, grotop%num_cg_KH_atomtypes
        if (grotop%cg_KH_atomtypes(j)%type_name &
            == molecule%atom_cls_name(ig)) then
          num_kh = num_kh + 1
          enefunc%cg_IDR_KH_list(num_kh) = i
          enefunc%cg_IDR_KH_list_inv(i) = num_kh
        end if
      end do
    end do

    return

  end subroutine setup_enefunc_cg_IDR_KH

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cg_IDR_HPS
  !> @brief        Defines parameters for CG protein KH model.
  !! @authors      CT, JJ
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cg_IDR_HPS(grotop, molecule, domain, enefunc)

    ! formal arguments
    type(s_grotop),  target, intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_domain),  target, intent(inout) :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    ! local variables
    integer              :: n_atoms, n_idr_region, num_hps
    integer              :: ioffset
    integer              :: i, j, k, l, ix, ig, ncel

    type(s_grotop_mol), pointer :: gromol

    ! ------------------------------------
    ! count IDR HPS "regions" in itp files
    ! ------------------------------------
    !
    n_idr_region = 0
    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        do k = 1, gromol%num_idr_hps
          n_idr_region = n_idr_region + 1
        end do
      end do
    end do

    ! switch on
    if ( n_idr_region > 0 ) &
      enefunc%cg_IDR_HPS_calc = .true.

    ncel = domain%num_cell_local + domain%num_cell_boundary

    ! -------------------------------------
    ! Set IDR regions if specified in input
    ! -------------------------------------
    !
    n_atoms = 0
    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        ioffset = n_atoms
        n_atoms = n_atoms + gromol%num_atoms
        do k = 1, gromol%num_idr_hps
          do l = gromol%idr_hps(k)%grp_start, gromol%idr_hps(k)%grp_end
            ix = domain%id_g2l(l+ioffset)
            if (ix > 0 .and. &
                ix <= domain%num_atom_domain+domain%num_atom_boundary) &
              domain%cg_IDR_HPS(ix) = 1
          end do
        end do
      end do
    end do

    num_hps = 0
    do i = 1, domain%num_atom_domain+domain%num_atom_boundary
      if (domain%cg_IDR_HPS(i) == 0) cycle
      ig = domain%id_l2g(i)
      do j = 1, grotop%num_cg_IDR_HPS_atomtypes
        if (grotop%cg_IDR_HPS_atomtypes(j)%type_name &
            == molecule%atom_cls_name(ig)) then
          num_hps = num_hps + 1
          domain%charge(i) = grotop%cg_IDR_HPS_atomtypes(j)%charge
          domain%NA_base_type(i) = NABaseTypeIDRHPS
          exit
        end if
      end do
    end do
    enefunc%num_cg_IDR_HPS = num_hps
#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, num_hps, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
#endif
    if (num_hps > 0) then
      Max_cg_IDR_HPS = max(num_hps*2,min_atom)
    else
      Max_cg_IDR_HPS = 0
    end if
    num_hps = Max_cg_IDR_HPS

    call alloc_enefunc(enefunc, EneFuncCGIDRHPSList, num_hps)
    call alloc_enefunc(enefunc, EneFuncCGIDRHPSInvList, MaxAtom_domain)

    num_hps = 0
    do i = 1, domain%num_atom_domain+domain%num_atom_boundary
      if (domain%cg_IDR_HPS(i) == 0) cycle
      ig = domain%id_l2g(i)
      do j = 1, grotop%num_cg_IDR_HPS_atomtypes
        if (grotop%cg_IDR_HPS_atomtypes(j)%type_name &
            == molecule%atom_cls_name(ig)) then
          num_hps = num_hps + 1
          enefunc%cg_IDR_HPS_list(num_hps) = i
          enefunc%cg_IDR_HPS_list_inv(i) = num_hps
          exit
        end if
      end do
    end do

    return

  end subroutine setup_enefunc_cg_IDR_HPS


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cg_IDR_MPIPI
  !> @brief        Defines parameters for CG protein KH model.
  !! @authors      CT
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cg_IDR_MPIPI(grotop, molecule, domain, enefunc)

    ! formal arguments
    type(s_grotop),  target, intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_domain),  target, intent(inout) :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer              :: n_atoms, n_idr_region, num_mpipi
    integer              :: ioffset
    integer              :: i, j, k, l, ix, ig, ncel

    type(s_grotop_mol), pointer :: gromol

    ! ------------------------------------
    ! count IDR MPIPI "regions" in itp files
    ! ------------------------------------
    !
    n_idr_region = 0
    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        do k = 1, gromol%num_idr_mpipi
          n_idr_region = n_idr_region + 1
        end do
      end do
    end do

    ! switch on
    if ( n_idr_region > 0 ) &
      enefunc%cg_IDR_MPIPI_calc = .true.

    ncel = domain%num_cell_local + domain%num_cell_boundary

    ! -------------------------------------
    ! Set IDR regions if specified in input
    ! -------------------------------------
    !
    n_atoms = 0
    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        ioffset = n_atoms
        n_atoms = n_atoms + gromol%num_atoms
        do k = 1, gromol%num_idr_mpipi
          do l = gromol%idr_mpipi(k)%grp_start, gromol%idr_mpipi(k)%grp_end
            ix = domain%id_g2l(l+ioffset)
            if (ix > 0 .and. &
                ix <= domain%num_atom_domain+domain%num_atom_boundary) &
              domain%cg_IDR_MPIPI(ix) = 1
          end do
        end do
      end do
    end do

    num_mpipi = 0
    do i = 1, domain%num_atom_domain+domain%num_atom_boundary
      if (domain%cg_IDR_MPIPI(i) == 0) cycle
      ig = domain%id_l2g(i)
      do j = 1, grotop%num_cg_IDR_MPIPI_atomtypes
        if (grotop%cg_IDR_MPIPI_atomtypes(j)%type_name &
            == molecule%atom_cls_name(ig)) then
          num_mpipi = num_mpipi + 1
          domain%charge(i) = grotop%cg_IDR_MPIPI_atomtypes(j)%charge
          domain%NA_base_type(i) = NABaseTypeIDRMPIPI
          exit
        end if
      end do
    end do
    enefunc%num_cg_IDR_MPIPI = num_mpipi
#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, num_mpipi, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
#endif
    Max_cg_IDR_MPIPI = num_mpipi * 2
    num_mpipi = Max_cg_IDR_MPIPI

    call alloc_enefunc(enefunc, EneFuncCGIDRMPIPIList, num_mpipi)
    call alloc_enefunc(enefunc, EneFuncCGIDRMPIPIInvList, MaxAtom_domain)

    num_mpipi = 0
    do i = 1, domain%num_atom_domain+domain%num_atom_boundary
      if (domain%cg_IDR_MPIPI(i) == 0) cycle
      ig = domain%id_l2g(i)
      do j = 1, grotop%num_cg_IDR_MPIPI_atomtypes
        if (grotop%cg_IDR_MPIPI_atomtypes(j)%type_name &
            == molecule%atom_cls_name(ig)) then
          num_mpipi = num_mpipi + 1
          enefunc%cg_IDR_MPIPI_list(num_mpipi) = i
          enefunc%cg_IDR_MPIPI_list_inv(i) = num_mpipi
          exit
        end if
      end do
    end do

    return

  end subroutine setup_enefunc_cg_IDR_MPIPI

 
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cgDNA_base_stack_alloc
  !> @brief        Defines BASE_STACK term in potential energy function.
  !! @authors      JJ
  !! @param[in]    grotop   : GROMACS topology informaiton
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cgDNA_base_stack_alloc(grotop, molecule, domain, &
                                                  enefunc)

    ! formal arguments
    type(s_grotop),  target, intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                    :: i, j, k, l
    integer                    :: ioffset, natom, nstack
    integer                    :: idx1, idx2, idx3, icel1, icel2, icel
    integer                    :: i_S5, i_B5, i_B3, i_B_last

    type(s_grotop_mol),    pointer :: gromol
    type(s_basestacktype), pointer :: base(:)
    real(wp),              pointer :: eps(:), sigma(:), theta_bs(:)
    integer,               pointer :: list(:,:), func(:)
    integer,               pointer :: ncel
    integer,               pointer :: id_g2l(:)
    integer,               pointer :: atom_2_cell(:)
    real(wp),              pointer :: num_atom(:)

    ncel           => domain%num_cell_local
    id_g2l         => domain%id_g2l
    atom_2_cell    => domain%atom_2_cell
    num_atom       => domain%num_atom_t0

    list           => enefunc%base_stack_list
    func           => enefunc%base_stack_func
    eps            => enefunc%base_stack_epsilon
    sigma          => enefunc%base_stack_sigma
    theta_bs       => enefunc%base_stack_theta_bs

    base           => grotop%basestacktypes

    enefunc%base_stack_alpha = 3.0_wp
    enefunc%base_stack_K     = 6.0_wp

    natom = 0
    nstack = 0

    ! check size for allocation
    !
    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        ioffset = natom
        natom   = natom + gromol%num_atoms
        do k = 1, gromol%num_atoms
          if (gromol%atoms(k)%atom_name == 'DB') then
            if (.not. enefunc%cg_infinite_DNA) then
              if (k < 4) cycle
              if (molecule%molecule_no(ioffset+k) /=        &
                  molecule%molecule_no(ioffset+k-3)) cycle
            end if
            if (k > 3) then
              if (molecule%molecule_no(ioffset+k) == &
                  molecule%molecule_no(ioffset+k-3)) then
                i_B5 = k - 3
              else
                do l = k+1, gromol%num_atoms
                  if (molecule%molecule_no(ioffset + k) == &
                      molecule%molecule_no(ioffset + l) .and. &
                      gromol%atoms(l)%atom_name == "DB" ) then
                    i_B_last = l
                  end if
                end do
                i_B5 = i_B_last
              end if
            else
              do l = k+1, gromol%num_atoms
                if (molecule%molecule_no(ioffset + k) == &
                    molecule%molecule_no(ioffset + l) .and. &
                    gromol%atoms(l)%atom_name == "DB" ) then
                  i_B_last = l
                end if
              end do
              i_B5 = i_B_last
            end if
            i_S5 = i_B5 - 1
            i_B3 = k
            idx1 = i_S5 + ioffset
            idx2 = i_B5 + ioffset
            idx3 = i_B3 + ioffset

            icel1 = id_g2l(idx2)
            icel2 = id_g2l(idx3)

            if (icel1 /= 0 .and. icel2 /= 0) then

              icel1 = atom_2_cell(icel1)
              icel2 = atom_2_cell(icel2)
              if (icel1 <= ncel .and. icel2 <= ncel) then
                icel = min(icel1,icel2)
              else
                if (num_atom(icel1) >= num_atom(icel2)) then
                  icel = icel2
                else if (num_atom(icel2) > num_atom(icel1)) then
                  icel = icel1
                end if
              end if

              if (icel > 0 .and. icel <= ncel) nstack = nstack + 1

            end if
          end if

        end do

      end do
    end do

    enefunc%num_stack_domain = nstack

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(k, MaxStack, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
#endif
    if (MaxStack > 0) MaxStack = max(MaxStack*2,min_atom)

    call alloc_enefunc(enefunc, EneFuncBaseStack, MaxStack, 1)

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(enefunc%num_stack_domain, enefunc%num_base_stack_all, &
                       1, mpi_integer, mpi_sum, mpi_comm_country, ierror)
#endif

    return

  end subroutine setup_enefunc_cgDNA_base_stack_alloc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cgDNA_base_stack
  !> @brief        Defines BASE_STACK term in potential energy function.
  !! @authors      JJ
  !! @param[in]    grotop   : GROMACS topology informaiton
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cgDNA_base_stack(grotop, molecule, domain, enefunc)

    ! formal arguments
    type(s_grotop),  target, intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                    :: i, j, k, l, s_idx
    integer                    :: ioffset, natom, nstack
    integer                    :: idx1, idx2, idx3, icel1, icel2, icel
    integer                    :: i_S5, i_B5, i_B3, i_B_last
    character(6)               :: basetype_5, basetype_3

    type(s_grotop_mol),    pointer :: gromol
    type(s_basestacktype), pointer :: base(:)
    real(wp),              pointer :: eps(:), sigma(:), theta_bs(:)
    integer,               pointer :: list(:,:), func(:)
    integer,               pointer :: ncel
    integer,               pointer :: id_g2l(:)
    integer,               pointer :: atom_2_cell(:)
    real(wp),              pointer :: num_atom(:)

    ncel           => domain%num_cell_local
    id_g2l         => domain%id_g2l
    atom_2_cell    => domain%atom_2_cell
    num_atom       => domain%num_atom_t0

    list           => enefunc%base_stack_list
    func           => enefunc%base_stack_func
    eps            => enefunc%base_stack_epsilon
    sigma          => enefunc%base_stack_sigma
    theta_bs       => enefunc%base_stack_theta_bs

    base           => grotop%basestacktypes

    enefunc%base_stack_alpha = 3.0_wp
    enefunc%base_stack_K     = 6.0_wp

    natom = 0
    nstack = 0
    do i = 1, grotop%num_molss

      gromol => grotop%molss(i)%moltype%mol

      do j = 1, grotop%molss(i)%count

        ioffset = natom
        natom   = natom + gromol%num_atoms

        do k = 1, gromol%num_atoms

          if (gromol%atoms(k)%atom_name == 'DB') then

            if (.not. enefunc%cg_infinite_DNA) then
              if (k < 4) cycle
              if (molecule%molecule_no(ioffset+k) /=        &
                  molecule%molecule_no(ioffset+k-3)) cycle
            end if
            if (k > 3) then
              if(molecule%molecule_no(ioffset+k) == &
                 molecule%molecule_no(ioffset+k-3)) then
                i_B5 = k - 3
              else
                do l = k+1, gromol%num_atoms
                  if (molecule%molecule_no(ioffset + k) == &
                      molecule%molecule_no(ioffset + l) .and. &
                      gromol%atoms(l)%atom_name == "DB" ) then
                    i_B_last = l
                  end if
                end do
                i_B5 = i_B_last
              end if
            else
              do l = k+1, gromol%num_atoms
                if (molecule%molecule_no(ioffset + k) == &
                    molecule%molecule_no(ioffset + l) .and. &
                    gromol%atoms(l)%atom_name == "DB" ) then
                  i_B_last = l
                end if
              end do
              i_B5 = i_B_last
            end if
            i_S5 = i_B5 - 1
            i_B3 = k
            idx1 = i_S5 + ioffset
            idx2 = i_B5 + ioffset
            idx3 = i_B3 + ioffset

            icel1 = id_g2l(idx2)
            icel2 = id_g2l(idx3)

            if (icel1 /= 0 .and. icel2 /= 0) then

              icel1 = atom_2_cell(icel1)
              icel2 = atom_2_cell(icel2)
              if (icel1 <= ncel .and. icel2 <= ncel) then
                icel = min(icel1,icel2)
              else
                if (num_atom(icel1) >= num_atom(icel2)) then
                  icel = icel2
                else if (num_atom(icel2) > num_atom(icel1)) then
                  icel = icel1
                end if
              end if

              if (icel > 0 .and. icel <= ncel) then

                basetype_5 = gromol%atoms(i_B5)%atom_type
                basetype_3 = gromol%atoms(i_B3)%atom_type

                nstack = nstack + 1
                s_idx = nstack
                list(1:3,s_idx) = (/idx1, idx2, idx3/)
                do l = 1, grotop%num_basestacktypes
                  if (base(l)%base_type5 == basetype_5 .and. &
                      base(l)%base_type3 == basetype_3) then
                    func (s_idx) = base(l)%func 
                    eps  (s_idx) = base(l)%epsilon * JOU2CAL
                    sigma(s_idx) = base(l)%sigma * 10.0_wp
                    theta_bs(s_idx) = base(l)%theta_bs * RAD
                  end if
                end do

              end if

            end if
          end if

        end do

      end do
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(enefunc%num_stack_domain, enefunc%num_base_stack_all, &
                       1, mpi_integer, mpi_sum, mpi_comm_country, ierror)
#endif

    return

  end subroutine setup_enefunc_cgDNA_base_stack

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cgDNA_nonb
  !> @brief        Defines nonbonded term for 3SPN.2C DNA model.
  !! @authors      CT
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cgDNA_nonb(ene_info, grotop, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_grotop),          intent(in)    :: grotop
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer              :: i, j, k
    integer              :: itype, jtype
    integer              :: alloc_stat
    logical              :: be_matched
    integer, allocatable :: base_type_2_atom_cls(:)
    real(wp)             :: sigma_i, sigma_j


    ! ----------------
    ! Model Parameters
    ! ----------------
    !
    enefunc%base_pair_alpha      = 2.0_wp
    enefunc%base_pair_K          = 12.0_wp
    enefunc%base_cross_alpha     = 4.0_wp
    enefunc%base_cross_K         = 8.0_wp
    enefunc%cgDNA_exv_epsilon    = 1.0_wp * JOU2CAL

    ! ===========================
    ! set cutoff and pairlistdist
    ! ===========================
    !
    enefunc%cg_cutoffdist_ele      = ene_info%cg_cutoffdist_ele
    enefunc%cg_cutoffdist_126      = ene_info%cg_cutoffdist_126
    enefunc%cg_cutoffdist_DNAbp    = ene_info%cg_cutoffdist_DNAbp
    enefunc%cg_pairlistdist_ele    = ene_info%cg_pairlistdist_ele
    enefunc%cg_pairlistdist_126    = ene_info%cg_pairlistdist_126
    enefunc%cg_pairlistdist_WF     = ene_info%cg_pairlistdist_WF 
    enefunc%cg_pairlistdist_PWMcos = ene_info%cg_pairlistdist_PWMcos
    enefunc%cg_pairlistdist_DNAbp  = ene_info%cg_pairlistdist_DNAbp
    enefunc%cg_pairlistdist_exv    = ene_info%cg_pairlistdist_exv
    enefunc%cg_ele_sol_T           = ene_info%cg_sol_temperature
    enefunc%cg_ele_sol_IC          = ene_info%cg_sol_ionic_strength
    enefunc%buffer_min = &
        min(enefunc%cg_pairlistdist_ele   - enefunc%cg_cutoffdist_ele,   &
            enefunc%cg_pairlistdist_DNAbp - enefunc%cg_cutoffdist_DNAbp, &
            enefunc%cg_pairlistdist_126   - enefunc%cg_cutoffdist_126)
    enefunc%buffer_min = enefunc%buffer_min / 2.0_wp

    ! =================================
    ! MAPPING base-pair <==> atom-class
    ! =================================
    ! count number of base-pair types
    ! and make the local mapping arrays
    !
    allocate(base_type_2_atom_cls(NABaseTypeNAMAX), &
        stat = alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    base_type_2_atom_cls(1:NABaseTypeNAMAX) = 0
    !
    do i = 1, grotop%num_atomtypes
      if (grotop%atomtypes(i)%type_name == 'DA' .or. &
          grotop%atomtypes(i)%type_name == 'DC' .or. &
          grotop%atomtypes(i)%type_name == 'DG' .or. &
          grotop%atomtypes(i)%type_name == 'DT' &
          ! Add something here if ~CG~ RNA will be used
          ) then
        if ( grotop%atomtypes(i)%type_name == 'DA' ) then
          base_type_2_atom_cls(NABaseTypeDBA) = i
        end if
        if ( grotop%atomtypes(i)%type_name == 'DC' ) then
          base_type_2_atom_cls(NABaseTypeDBC) = i
        end if
        if ( grotop%atomtypes(i)%type_name == 'DG' ) then
          base_type_2_atom_cls(NABaseTypeDBG) = i
        end if
        if ( grotop%atomtypes(i)%type_name == 'DT' ) then
          base_type_2_atom_cls(NABaseTypeDBT) = i
        end if
      else if ( grotop%atomtypes(i)%type_name == 'DP' ) then
        base_type_2_atom_cls(NABaseTypeDP) = i
      else if ( grotop%atomtypes(i)%type_name == 'DS' ) then
        base_type_2_atom_cls(NABaseTypeDS) = i
      end if
    end do

    ! --------------------------------------
    ! transfer params from grotop to enefunc
    ! --------------------------------------
    !
    call alloc_enefunc(enefunc, EneFuncBasePair, NABaseTypeBMAX)
    call alloc_enefunc(enefunc, EneFuncCGDNAExv, NABaseTypeNAMAX)

    ! enefunc%base_pair_xxx(:)
    do k = 1, grotop%num_basepairtypes
      do i = 1, NABaseTypeBMAX
        itype = base_type_2_atom_cls(i)
        if ( itype == 0 ) &
            cycle

        if (grotop%atomtypes(itype)%type_name == &
            grotop%basepairtypes(k)%base_type_a) then

          enefunc%base_pair_theta_1(i) = &
              grotop%basepairtypes(k)%theta_1 * RAD
          enefunc%base_pair_theta_2(i) = &
              grotop%basepairtypes(k)%theta_2 * RAD
          enefunc%base_pair_theta_3(i) = &
              grotop%basepairtypes(k)%theta_3 * RAD
          enefunc%base_pair_phi_1(i)   = &
              grotop%basepairtypes(k)%phi_1 * RAD
          enefunc%base_pair_sigma(i)   = &
              grotop%basepairtypes(k)%sigma * 10.0_wp
          enefunc%base_pair_epsilon(i) = &
              grotop%basepairtypes(k)%epsilon * JOU2CAL
          exit

        end if
      end do
    end do

    ! enefunc%base_cross_xxx(:, :)
    do k = 1, grotop%num_basecrosstypes
      be_matched = .false.
      do i = 1, NABaseTypeBMAX
        itype = base_type_2_atom_cls(i)
        if ( itype == 0 ) &
            cycle
        do j = 1, NABaseTypeBMAX
          jtype = base_type_2_atom_cls(j)
          if ( jtype == 0 ) &
              cycle
          if (grotop%atomtypes(itype)%type_name ==    &
              grotop%basecrosstypes(k)%base_type_a .and. &
              grotop%atomtypes(jtype)%type_name ==    &
              grotop%basecrosstypes(k)%base_type_b) then

            if ( grotop%basecrosstypes(k)%func == 1 ) then

              enefunc%base_cross_1_epsilon(i, j)  = &
                  grotop%basecrosstypes(k)%epsilon * JOU2CAL
              enefunc%base_cross_1_sigma(i, j)    = &
                  grotop%basecrosstypes(k)%sigma * 10.0_wp
              enefunc%base_cross_1_theta_cs(i, j) = &
                  grotop%basecrosstypes(k)%theta_cs * RAD
            else
              enefunc%base_cross_2_epsilon(i, j)  = &
                  grotop%basecrosstypes(k)%epsilon * JOU2CAL
              enefunc%base_cross_2_sigma(i, j)    = &
                  grotop%basecrosstypes(k)%sigma * 10.0_wp
              enefunc%base_cross_2_theta_cs(i, j) = &
                  grotop%basecrosstypes(k)%theta_cs * RAD

            end if
            be_matched = .true.
            exit
          end if
        end do                ! loop j
        if ( be_matched ) exit
      end do                  ! loop i
    end do

    ! --------------------------
    ! Find out the base pairs!!!
    ! --------------------------
    !
    enefunc%base_pair_is_WC(1:NABaseTypeBMAX, 1:NABaseTypeBMAX) = .false.
    do i = 1, NABaseTypeBMAX
      itype = base_type_2_atom_cls(i)
      if ( itype == 0 ) &
          cycle
      do j = 1, NABaseTypeBMAX
        jtype = base_type_2_atom_cls(j)
        if ( jtype == 0 ) &
            cycle
        if (grotop%atomtypes(itype)%type_name == "DA" .and. &
            grotop%atomtypes(jtype)%type_name == "DT" ) then
          enefunc%base_pair_is_WC(i, j) = .true.
        else if (grotop%atomtypes(itype)%type_name == "DC" .and. &
            grotop%atomtypes(jtype)%type_name == "DG" ) then
          enefunc%base_pair_is_WC(i, j) = .true.
        else if (grotop%atomtypes(itype)%type_name == "DG" .and. &
            grotop%atomtypes(jtype)%type_name == "DC" ) then
          enefunc%base_pair_is_WC(i, j) = .true.
        else if (grotop%atomtypes(itype)%type_name == "DT" .and. &
            grotop%atomtypes(jtype)%type_name == "DA" ) then
          enefunc%base_pair_is_WC(i, j) = .true.
        end if
      end do                ! loop j
    end do                  ! loop i

    ! enefunc%cgDNA_exv_sigma(:)
    !
    do i = 1, NABaseTypeNAMAX
      itype = base_type_2_atom_cls(i)
      if ( itype == 0 ) &
        cycle

      do k = 1, grotop%num_cgdnaexvtypes
        if (grotop%atomtypes(itype)%type_name == &
            grotop%cgdnaexvtypes(k)%base_type) then
          sigma_i = grotop%cgdnaexvtypes(k)%sigma * 10.0_wp
          exit
        end if
      end do

      do j = 1, NABaseTypeNAMAX
        jtype = base_type_2_atom_cls(j)
        if ( jtype == 0 ) &
          cycle

        do k = 1, grotop%num_cgdnaexvtypes
          if (grotop%atomtypes(jtype)%type_name == &
              grotop%cgdnaexvtypes(k)%base_type) then
            sigma_j = grotop%cgdnaexvtypes(k)%sigma * 10.0_wp
            exit
          end if
        end do

        enefunc%cgDNA_exv_sigma(i, j) = 0.5_wp * (sigma_i + sigma_j)

      end do                    ! loop j
    end do                      ! loop i

    deallocate(base_type_2_atom_cls, stat = alloc_stat)

    return

  end subroutine setup_enefunc_cgDNA_nonb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_nonb_lb
  !> @brief        Defines NON-BOND term in potential energy function.
  !! @authors      JJ
  !! @param[in]    grotop      : GROMACS information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_nonb_lb(grotop, molecule, domain, enefunc, comm)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_comm),            intent(inout) :: comm  

    integer                  :: i, ix, ixx, jx, jxx, k, kx, kxx, l
    integer                  :: num_base, num_dna, start_i 

    ! create native contact list
    if (enefunc%forcefield == ForcefieldKBGO .or. &
        enefunc%forcefield == ForcefieldAAGO .or. &
        enefunc%forcefield == ForcefieldCAGO .or. &
        enefunc%forcefield == ForcefieldSOFT .or. &
        enefunc%forcefield == ForcefieldRESIDCG) then
      call setup_enefunc_contact(grotop, domain, enefunc, comm)
    end if

    ! base type and molecular type
    !
    if (enefunc%forcefield == ForcefieldRESIDCG) then

      ! allocation of base
      !
      num_dna  = 0
      num_base = 0
      do i = 1, domain%num_atom_domain+domain%num_atom_boundary
        l = domain%atom_cls_no(i)
        k = enefunc%atom_cls_2_base_type(l)
        if (k <= NABaseTypeDBMax .or. k == NABaseTypeDP .or. &
            k == NABaseTypeDS) then
          num_dna = num_dna + 1
          if (k <= NABaseTypeDBMax) num_base = num_base + 1
        end if
      end do
      enefunc%num_cg_base = num_base
      enefunc%num_cg_DNA  = num_dna

#ifdef HAVE_MPI_GENESIS
      call mpi_allreduce(mpi_in_place, num_base, 1, mpi_integer, &
                         mpi_max, mpi_comm_country, ierror)
      call mpi_allreduce(mpi_in_place, num_dna , 1, mpi_integer, &
                         mpi_max, mpi_comm_country, ierror)
#endif
      if (num_dna > 0) then
        enefunc%cg_DNA_base_pair_calc = .true.
        enefunc%cg_DNA_exv_calc       = .true.
      end if

      if (num_base > 0) then
        Max_cg_base = max(num_base*2,min_atom)
      else
        Max_cg_base = 0
      end if
      if (num_dna > 0) then
        Max_cg_dna  = max(num_dna*2,min_atom)
      else
        Max_cg_dna  = 0
      end if
      num_dna     = Max_cg_dna
      num_base    = Max_cg_base

      call alloc_enefunc(enefunc, EneFuncCGDNAList,  num_dna)
      call alloc_enefunc(enefunc, EneFuncCGDNAInvList, MaxAtom_domain)
      call alloc_enefunc(enefunc, EneFuncCGBaseList, num_base)
      call alloc_enefunc(enefunc, EneFuncCGBaseInvList, MaxAtom_domain)

      num_dna  = 0
      num_base = 0
      do i = 1, domain%num_cell_local+domain%num_cell_boundary
        kx = 0
        jx = 0
        start_i = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)
          ixx = start_i + ix
          l = domain%atom_cls_no(ixx)
          k = enefunc%atom_cls_2_base_type(l)
          if (domain%NA_base_type(ixx) == 0) &
            domain%NA_base_type(ixx) = k
          if (k <= NABaseTypeDBMax .or. k == NABaseTypeDP .or. &
              k == NABaseTypeDS) then
            num_dna = num_dna + 1
            enefunc%cg_dna_list(num_dna) = ixx
            enefunc%cg_dna_list_inv(ixx) = num_dna
            if (k <= NABaseTypeDBMax) then
              num_base = num_base + 1
              enefunc%cg_base_list(num_base) = ixx
              enefunc%cg_base_list_inv(ixx) = num_base
              domain%dna_check(ixx) = 1
              kx = kx + 1
              kxx = kx + start_i
              domain%base_list(kxx) = ixx
            else if (k == NABaseTypeDP) then
              domain%dna_check(ixx) = 2
              jx = jx + 1
              jxx = jx + start_i
              domain%phos_list(jxx) = ixx
            else if (k == NABaseTypeDS) then
              domain%dna_check(ixx) = 3
            end if
          else
            domain%dna_check(ixx) = 0
          end if
        end do
        domain%num_base(i) = kx
        domain%num_phos(i) = jx
      end do

    end if


    return

  end subroutine setup_enefunc_nonb_lb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_nonb
  !> @brief        Defines NON-BOND term in potential energy function.
  !! @authors      JJ
  !! @param[in]    grotop      : GROMACS information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_nonb(ene_info, grotop, molecule, domain, enefunc, &
                                comm)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_comm),            intent(inout) :: comm   

    ! local variables
    real(wp)                 :: eps, sig, kij, ei, ej, si, sj, ki, kj
    real(wp)                 :: c6i, c6j, c12i, c12j, c6, c10, c12
    real(wp)                 :: vi, vj, wi, wj, vij, wij
    integer                  :: nnonb, ncel, i, j, k, l, excl_level
    integer                  :: ix, ixx, jx, jxx, kx, kxx, start_i
    integer                  :: cls_local, num_dna, num_base

    integer,    allocatable  :: check_cls(:)
    integer,    allocatable  :: atmcls_map_g2l(:), atmcls_map_l2g(:)
    integer,    allocatable  :: atom_cls_2_base_type(:)
    real(wp),   allocatable  :: nb14_lj6(:,:), nb14_lj10(:,:), nb14_lj12(:,:)
    real(wp),   allocatable  :: nonb_lj6(:,:), nonb_lj10(:,:), nonb_lj12(:,:)
    real(wp),   allocatable  :: nonb_aicg_eps(:,:), nonb_aicg_sig(:,:)
    real(wp),   allocatable  :: cg_KH_epsilon(:,:), cg_KH_sigma(:,:)
    real(wp),   allocatable  :: cg_IDR_KH_epsilon(:,:), cg_IDR_KH_sigma(:,:)
    real(wp),   allocatable  :: cg_IDR_HPS_lambda(:,:), cg_IDR_HPS_sigma(:,:)
    real(wp),   allocatable  :: cg_WangFrenkel_sigma(:,:), cg_WangFrenkel_epsilon(:,:)
    real(wp),   allocatable  :: cg_WangFrenkel_mu(:,:),    cg_WangFrenkel_nu(:,:)

    enefunc%num_atom_cls = grotop%num_atomtypes
    enefunc%fudge_lj     = grotop%defaults%fudge_lj
    enefunc%fudge_qq     = grotop%defaults%fudge_qq
    enefunc%cg_cutoffdist_ele      = ene_info%cg_cutoffdist_ele
    enefunc%cg_cutoffdist_vdw      = ene_info%cg_cutoffdist_vdw
    enefunc%cg_pairlistdist_ele    = ene_info%cg_pairlistdist_ele
    enefunc%cg_pairlistdist_vdw    = ene_info%cg_pairlistdist_vdw
    enefunc%cg_switchdist_vdw      = ene_info%cg_switchdist_vdw

    ncel                 = domain%num_cell_local
    ELECOEF              = ELECOEF_GROMACS

    ! set lennard-jones parameters
    !
    nnonb = enefunc%num_atom_cls

    allocate(check_cls(nnonb),        &
             atmcls_map_g2l(nnonb),   &
             atmcls_map_l2g(nnonb),   &
             nb14_lj6 (nnonb, nnonb), &
             nb14_lj10(nnonb, nnonb), &
             nb14_lj12(nnonb, nnonb), &
             nonb_lj6 (nnonb, nnonb), &
             nonb_lj10(nnonb, nnonb), &
             nonb_lj12(nnonb, nnonb))

    if (enefunc%forcefield == ForcefieldRESIDCG) &
      allocate(nonb_aicg_eps(nnonb, nnonb),     &
               nonb_aicg_sig(nnonb, nnonb),     &
               atom_cls_2_base_type(nnonb))
    if (enefunc%forcefield == ForcefieldRESIDCG) then
      if (enefunc%cg_KH_calc) &
        allocate(cg_KH_epsilon(nnonb, nnonb),   &
                 cg_KH_sigma  (nnonb, nnonb))
      if (enefunc%cg_IDR_KH_calc) &
        allocate(cg_IDR_KH_epsilon(nnonb, nnonb),&
                 cg_IDR_KH_sigma  (nnonb, nnonb))
      if (enefunc%cg_IDR_HPS_calc) &
        allocate(cg_IDR_HPS_lambda(nnonb, nnonb),&
                 cg_IDR_HPS_sigma (nnonb, nnonb))
      if (enefunc%cg_IDR_MPIPI_calc) &
        allocate(cg_WangFrenkel_epsilon(nnonb, nnonb),&
                 cg_WangFrenkel_sigma  (nnonb, nnonb),&
                 cg_WangFrenkel_mu     (nnonb, nnonb),&
                 cg_WangFrenkel_nu     (nnonb, nnonb))
    end if

    check_cls(1:nnonb)          = 0

    do i = 1, nnonb
      do j = 1, nnonb

        ! combination rule
        !
        vij = 0.0_wp
        wij = 0.0_wp
        c12 = 0.0_wp
        c10 = 0.0_wp
        c6  = 0.0_wp

        if (grotop%num_nbonparms > 0) then

          do k = 1, grotop%num_nbonparms
            if (grotop%atomtypes(i)%type_name == &
                  grotop%nbonparms(k)%atom_type1 .and. &
                grotop%atomtypes(j)%type_name == &
                  grotop%nbonparms(k)%atom_type2 .or.  &
                grotop%atomtypes(j)%type_name == &
                  grotop%nbonparms(k)%atom_type1 .and. &
                grotop%atomtypes(i)%type_name == &
                  grotop%nbonparms(k)%atom_type2) then

              vij = grotop%nbonparms(k)%v
              wij = grotop%nbonparms(k)%w

              exit
            end if
          end do

          if (grotop%defaults%combi_rule == 2) then

            sig = vij * 10.0_wp
            eps = wij * JOU2CAL

            c6  = 4.0_wp * eps * (sig ** 6)
            c10 = 0.0_wp
            c12 = 4.0_wp * eps * (sig ** 12)

          else ! combi_rule = 1 or 3

            c6  = vij * 1000000.0_wp * JOU2CAL
            c10 = 0.0_wp
            c12 = wij * 1000000.0_wp * 1000000.0_wp * JOU2CAL

          end if

        end if

        if (c6 == 0.0_wp .and. c10 == 0.0_wp .and. c12 == 0.0_wp) then

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
            eps = sqrt(ei * ej)

            if (enefunc%forcefield == ForcefieldRESIDCG) then

              nonb_aicg_eps(i,j) = eps
              nonb_aicg_sig(i,j) = sig * ene_info%cg_exv_sigma_scaling
              c6  = 4.0_wp * eps * (sig ** 6)
              c10 = 0.0_wp
              c12 = 4.0_wp * eps * (sig ** 12)

            else if (enefunc%forcefield /= ForcefieldSOFT) then

              c6  = 4.0_wp * eps * (sig ** 6)
              c10 = 0.0_wp
              c12 = 4.0_wp * eps * (sig ** 12)

            else

              ki = grotop%atomtypes(i)%khh * JOU2CAL / 100.0_wp
              kj = grotop%atomtypes(j)%khh * JOU2CAL / 100.0_wp

              kij = sqrt(ki*kj)

              c6  = kij
              c10 = sig
              c12 = eps

            end if

          else ! combi_rule == 1 or 3

            c6i  = vi * 1.0E6_wp * JOU2CAL
            c6j  = vj * 1.0E6_wp * JOU2CAL

            c12i = wi * 1.0E12_wp * JOU2CAL
            c12j = wj * 1.0E12_wp * JOU2CAL

            c6  = sqrt(c6i  * c6j)
            c10 = 0.0_wp
            c12 = sqrt(c12i * c12j)

            if (enefunc%forcefield == ForcefieldSOFT) then
              call error_msg('Error: Setup_Enefunc_Nonb - Combination rule '//&
		             'should be 2 for SOFT.')
            end if

          end if

        end if

        if (main_rank) then
          if (c6 == 0.0_wp .and. c10 == 0.0_wp .and. c12 == 0.0_wp) &
            write(MsgOut,'(A,A,A)') &
              'Setup_Enefunc_Nonb> WARNING, combination is not found.',&
               grotop%atomtypes(i)%type_name, grotop%atomtypes(j)%type_name
        endif

        ! set parameters
        !
        nb14_lj12(i,j) = c12
        nb14_lj10(i,j) = c10
        nb14_lj6 (i,j) = c6

        nonb_lj12(i,j) = c12
        nonb_lj10(i,j) = c10
        nonb_lj6 (i,j) = c6

      end do
    end do

    ! base type
    if (enefunc%forcefield == ForcefieldRESIDCG) then

      call alloc_enefunc(enefunc, EneFuncCGDNAmol, molecule%num_molecules)
      call setup_enefunc_cg_basetype(grotop, molecule, atom_cls_2_base_type, &
                                     enefunc)

      if (enefunc%cg_KH_calc) then
        do i = 1, nnonb
          do j = 1, nnonb

            do k = 1, grotop%num_cg_pair_MJ_eps
              if (grotop%atomtypes(i)%type_name ==           &
                  grotop%cg_pair_MJ_eps(k)%type_name_1 .and. &
                  grotop%atomtypes(j)%type_name ==           &
                  grotop%cg_pair_MJ_eps(k)%type_name_2 ) then
                cg_KH_epsilon(i,j) = grotop%cg_pair_MJ_eps(k)%epsilon
              end if
            end do

            do k = 1, grotop%num_cg_KH_atomtypes
              if (grotop%atomtypes(i)%type_name ==          &
                  grotop%cg_KH_atomtypes(k)%type_name)      &
                si = grotop%cg_KH_atomtypes(k)%sigma * 5.0_wp
            end do
            do k = 1, grotop%num_cg_KH_atomtypes 
              if (grotop%atomtypes(j)%type_name ==          &
                  grotop%cg_KH_atomtypes(k)%type_name)      &
                sj = grotop%cg_KH_atomtypes(k)%sigma * 5.0_wp
            end do
            cg_KH_sigma(i,j) = (si+sj)**6            

          end do
        end do
      end if

      if (enefunc%cg_IDR_KH_calc) then

        wi = enefunc%cg_KH_mod_D_lambda
        wj = enefunc%cg_KH_mod_D_eps_0         
        do i = 1, nnonb
          do j = 1, nnonb

            do k = 1, grotop%num_cg_pair_MJ_eps
              if (grotop%atomtypes(i)%type_name ==           &
                  grotop%cg_pair_MJ_eps(k)%type_name_1 .and. &
                  grotop%atomtypes(j)%type_name ==           &
                  grotop%cg_pair_MJ_eps(k)%type_name_2 ) then
                cg_IDR_KH_epsilon(i,j) = wi*(grotop%cg_pair_MJ_eps(k)%epsilon-wj)*0.593_wp
              end if
            end do

            do k = 1, grotop%num_cg_KH_atomtypes
              if (grotop%atomtypes(i)%type_name ==          &
                  grotop%cg_KH_atomtypes(k)%type_name)      &
                si = grotop%cg_KH_atomtypes(k)%sigma * 5.0_wp
            end do
            do k = 1, grotop%num_cg_KH_atomtypes
              if (grotop%atomtypes(j)%type_name ==          &
                  grotop%cg_KH_atomtypes(k)%type_name)      &
                sj = grotop%cg_KH_atomtypes(k)%sigma * 5.0_wp
            end do
            cg_IDR_KH_sigma(i,j) = (si+sj)**6

          end do
        end do
      end if

      if (enefunc%cg_IDR_HPS_calc) then

        do i = 1, nnonb
          do j = 1, nnonb

            do k = 1, grotop%num_cg_IDR_HPS_atomtypes
              if (grotop%atomtypes(i)%type_name ==  &
                  grotop%cg_IDR_HPS_atomtypes(k)%type_name) then
                si = grotop%cg_IDR_HPS_atomtypes(k)%sigma * 5.0_wp
                ei = grotop%CG_IDR_HPS_atomtypes(k)%lambda * 0.5_wp
              end if
            end do
            do k = 1, grotop%num_cg_IDR_HPS_atomtypes
              if (grotop%atomtypes(j)%type_name ==  &
                  grotop%cg_IDR_HPS_atomtypes(k)%type_name) then
                sj = grotop%cg_IDR_HPS_atomtypes(k)%sigma * 5.0_wp
                ej = grotop%CG_IDR_HPS_atomtypes(k)%lambda * 0.5_wp
              end if
            end do
            cg_IDR_HPS_sigma(i,j) = (si+sj)**6
            cg_IDR_HPS_lambda(i,j) = ei+ej

          end do
        end do
      end if

      if (enefunc%cg_IDR_MPIPI_calc) then

        do i = 1, nnonb
          do j = 1, nnonb

            do k = 1, grotop%num_cg_IDR_MPIPI_pairs
              if (grotop%atomtypes(i)%type_name == &
                  grotop%cg_IDR_MPIPI_pairs(k)%type_name_1 .and. &
                  grotop%atomtypes(j)%type_name == &
                  grotop%cg_IDR_MPIPI_pairs(k)%type_name_2 ) then
                ! i - j
                cg_WangFrenkel_epsilon(i, j) = grotop%cg_IDR_MPIPI_pairs(k)%epsilon
                cg_WangFrenkel_sigma(i, j) = grotop%cg_IDR_MPIPI_pairs(k)%sigma * 10.0_wp
                cg_WangFrenkel_mu(i, j) = grotop%cg_IDR_MPIPI_pairs(k)%mu
                cg_WangFrenkel_nu(i, j) = 1

                cg_WangFrenkel_epsilon(j, i) = grotop%cg_IDR_MPIPI_pairs(k)%epsilon
                cg_WangFrenkel_sigma(j, i) = grotop%cg_IDR_MPIPI_pairs(k)%sigma * 10.0_wp
                cg_WangFrenkel_mu(j, i) = grotop%cg_IDR_MPIPI_pairs(k)%mu
                cg_WangFrenkel_nu(j, i) = 1
              end if
            end do                  !k

          end do                    !j
        end do                      !i
      end if

    end if

    ! create native contact list
    if (enefunc%forcefield == ForcefieldKBGO .or. &
        enefunc%forcefield == ForcefieldAAGO .or. &
        enefunc%forcefield == ForcefieldCAGO .or. &
        enefunc%forcefield == ForcefieldSOFT .or. &
        enefunc%forcefield == ForcefieldRESIDCG) then
      call setup_enefunc_contact(grotop, domain, enefunc, comm)
    end if

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
          call error_msg('Error: Setup_Enefunc_Nonb - Multiple '//&
		         '"exclude_nbon" are not supported.')
        end if

      end do

    end if

    ! check the usage of atom class
    !
    do i = 1, molecule%num_atoms
      k = molecule%atom_cls_no(i)
      if (k < 1) then
        call error_msg('Error: Setup_Enefunc_Nonb - Atom class is '//&
		       'not defined: "'//trim(molecule%atom_cls_name(i))//'".')
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
        enefunc%nb14_lj10(i,j) = nb14_lj10(ix,jx)
        enefunc%nb14_lj6 (i,j) = nb14_lj6 (ix,jx)
        enefunc%nonb_lj12(i,j) = nonb_lj12(ix,jx)
        enefunc%nonb_lj10(i,j) = nonb_lj10(ix,jx)
        enefunc%nonb_lj6 (i,j) = nonb_lj6 (ix,jx)
        if (enefunc%forcefield == ForcefieldRESIDCG) then
          enefunc%nonb_aicg_eps(i,j) = nonb_aicg_eps(ix,jx)
          enefunc%nonb_aicg_sig(i,j) = nonb_aicg_sig(ix,jx)
          if(enefunc%cg_KH_calc) then
            enefunc%cg_KH_epsilon(i,j) = cg_KH_epsilon(ix,jx)
            enefunc%cg_KH_sigma(i,j)   = cg_KH_sigma  (ix,jx) 
          end if
          if(enefunc%cg_IDR_KH_calc) then
            enefunc%cg_IDR_KH_epsilon(i,j) = cg_IDR_KH_epsilon(ix,jx)
            enefunc%cg_IDR_KH_sigma  (i,j) = cg_IDR_KH_sigma  (ix,jx)
          end if
          if(enefunc%cg_IDR_HPS_calc) then
            enefunc%cg_IDR_HPS_lambda(i,j) = cg_IDR_HPS_lambda(ix,jx)
            enefunc%cg_IDR_HPS_sigma (i,j) = cg_IDR_HPS_sigma (ix,jx)
          end if
          if(enefunc%cg_IDR_MPIPI_calc) then
            enefunc%cg_WangFrenkel_epsilon(i, j) = cg_WangFrenkel_epsilon(ix, jx)
            enefunc%cg_WangFrenkel_sigma(i, j)   = cg_WangFrenkel_sigma(ix, jx)
            enefunc%cg_WangFrenkel_mu(i, j)      = cg_WangFrenkel_mu(ix, jx)
            enefunc%cg_WangFrenkel_nu(i, j)      = cg_WangFrenkel_nu(ix, jx)
          end if
        end if
      end do
      if (enefunc%forcefield == ForcefieldRESIDCG) &
        enefunc%atom_cls_2_base_type(i) = atom_cls_2_base_type(ix)
    end do

    ! update domain information
    !
    do i = 1, domain%num_atom_domain+domain%num_atom_boundary
      domain%atom_cls_no(i) = atmcls_map_g2l(domain%atom_cls_no(i))
    end do

    ! base type and molecular type
    !
    if (enefunc%forcefield == ForcefieldRESIDCG) then

      ! allocation of base
      !
      num_dna  = 0
      num_base = 0
      do i = 1, domain%num_atom_domain+domain%num_atom_boundary
        l = domain%atom_cls_no(i)
        k = enefunc%atom_cls_2_base_type(l)
        if (k <= NABaseTypeDBMax .or. k == NABaseTypeDP .or. &
            k == NABaseTypeDS) then
          num_dna = num_dna + 1
          if (k <= NABaseTypeDBMax) num_base = num_base + 1
        end if
      end do
      enefunc%num_cg_base = num_base
      enefunc%num_cg_DNA  = num_dna

#ifdef HAVE_MPI_GENESIS
      call mpi_allreduce(mpi_in_place, num_base, 1, mpi_integer, &
                         mpi_max, mpi_comm_country, ierror)
      call mpi_allreduce(mpi_in_place, num_dna , 1, mpi_integer, &
                         mpi_max, mpi_comm_country, ierror)
#endif
      if (num_dna > 0) then
        enefunc%cg_DNA_base_pair_calc = .true.
        enefunc%cg_DNA_exv_calc       = .true.
      end if

      if (num_base > 0) then
        Max_cg_base = max(num_base*2, min_atom)
      else
        Max_cg_base = 0
      end if
      if (num_dna > 0) then
        Max_cg_dna  = max(num_dna*2, min_atom)
      else
        Max_cg_dna  = 0
      end if
      num_dna     = Max_cg_dna
      num_base    = Max_cg_base
     
      call alloc_enefunc(enefunc, EneFuncCGDNAList,  num_dna) 
      call alloc_enefunc(enefunc, EneFuncCGDNAInvList, MaxAtom_domain) 
      call alloc_enefunc(enefunc, EneFuncCGBaseList, num_base) 
      call alloc_enefunc(enefunc, EneFuncCGBaseInvList, MaxAtom_domain) 

      num_dna  = 0
      num_base = 0
      do i = 1, domain%num_cell_local+domain%num_cell_boundary
        kx = 0
        jx = 0
        start_i = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)
          ixx = start_i + ix
          l = domain%atom_cls_no(ixx)
          k = enefunc%atom_cls_2_base_type(l)
          if (domain%NA_base_type(ixx) == 0) &
            domain%NA_base_type(ixx) = k
          if (k <= NABaseTypeDBMax .or. k == NABaseTypeDP .or. &
              k == NABaseTypeDS) then
            num_dna = num_dna + 1
            enefunc%cg_dna_list(num_dna) = ixx
            enefunc%cg_dna_list_inv(ixx) = num_dna
            if (k <= NABaseTypeDBMax) then
              num_base = num_base + 1
              enefunc%cg_base_list(num_base) = ixx
              enefunc%cg_base_list_inv(ixx) = num_base
              domain%dna_check(ixx) = 1
              kx = kx + 1
              kxx = kx + start_i
              domain%base_list(kxx) = ixx
            else if (k == NABaseTypeDP) then
              domain%dna_check(ixx) = 2
              jx = jx + 1
              jxx = jx + start_i
              domain%phos_list(jxx) = ixx
            else if (k == NABaseTypeDS) then
              domain%dna_check(ixx) = 3
            end if
          else
            domain%dna_check(ixx) = 0
          end if
        end do
        domain%num_base(i) = kx
        domain%num_phos(i) = jx
      end do

    end if

    deallocate(check_cls,      &
               atmcls_map_g2l, &
               atmcls_map_l2g, &
               nb14_lj6,       &
               nb14_lj10,      &
               nb14_lj12,      &
               nonb_lj6,       &
               nonb_lj10,      &
               nonb_lj12)
    if (enefunc%forcefield == ForcefieldRESIDCG) then
      deallocate(nonb_aicg_eps,  &
                 nonb_aicg_sig,  &
                 atom_cls_2_base_type)
      if (enefunc%cg_KH_calc) &
        deallocate(cg_KH_epsilon, cg_KH_sigma)
      if (enefunc%cg_IDR_KH_calc) &
        deallocate(cg_IDR_KH_epsilon, cg_IDR_KH_sigma)
      if (enefunc%cg_IDR_HPS_calc) &
        deallocate(cg_IDR_HPS_lambda, cg_IDR_HPS_sigma)
      if (enefunc%cg_IDR_MPIPI_calc)       &
        deallocate(cg_WangFrenkel_epsilon, &
                   cg_WangFrenkel_sigma,   &
                   cg_WangFrenkel_mu,      &
                   cg_WangFrenkel_nu)
    end if

    enefunc%num_atom_cls = cls_local

    ! treatment for 1-2, 1-3, 1-4 interactions
    !
    ncel        = domain%num_cell_local

    call alloc_enefunc(enefunc, EneFuncNonb, MaxAtom_domain)

    if (enefunc%forcefield == ForcefieldAAGO .or. &
        enefunc%forcefield == ForcefieldCAGO .or. &
        enefunc%forcefield == ForcefieldRESIDCG) then

      call count_nonb_excl_go(domain, enefunc)

    else if (enefunc%forcefield == ForcefieldGroMartini) then

      call count_nonb_excl(.true., domain, enefunc)

    end if

    return

  end subroutine setup_enefunc_nonb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cg_ele
  !> @brief        Defines electrostatic parameters.
  !! @authors      CT
  !! @param[in]    ene_info : ENERGY section control parameters
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cg_ele(ene_info, grotop, molecule, domain, &
                                  enefunc)

    ! formal arguments
    type(s_ene_info),         intent(in)    :: ene_info
    type(s_grotop),           intent(in)    :: grotop
    type(s_molecule),         intent(in)    :: molecule
    type(s_domain),   target, intent(inout) :: domain
    type(s_enefunc),          intent(inout) :: enefunc

    ! local variables
    integer              :: n_mols
    real(wp)             :: e_T, a_C
    real(wp)             :: sol_T, sol_C
    real(wp)             :: debye_length
    integer              :: i, j, k, ix, base_type, start_i

    integer,          pointer :: natom(:), ncharge(:)
    integer,          pointer :: id_l2g(:)
    integer,          pointer :: chain_id(:)
    integer,          pointer :: NA_base_type(:)
    integer(1),       pointer :: charge_type(:)
    integer,          pointer :: ncell, nboundary
    real(wip),        pointer :: coord(:,:)
    real(wp),         pointer :: charge(:)

    coord          => domain%coord
    charge         => domain%charge
    natom          => domain%num_atom
    ncharge        => domain%num_charge
    ncell          => domain%num_cell_local
    nboundary      => domain%num_cell_boundary
    chain_id       => domain%mol_chain_id
    NA_base_type   => domain%NA_base_type
    charge_type    => domain%charge_type
    ncell          => domain%num_cell_local
    id_l2g         => domain%id_l2g

    ! --------------------
    ! electrostatic params
    ! --------------------
    !
    if ( grotop%num_cgelemolpairs > 0 ) then
      enefunc%cg_ele_calc       = .true.
    end if

    enefunc%cg_ele_coef = ELEMENT_CHARGE1 * ELEMENT_CHARGE1 /   &
        (4.0_wp * PI * ELECTRIC_CONST1) * AVOGADRO1 * JOU2CAL * 1e4_wp

    sol_T = enefunc%cg_ele_sol_T
    sol_C = enefunc%cg_ele_sol_IC
    e_T   = 2.494e2_wp - 7.88e-1_wp * sol_T &
        + 7.2e-4_wp * sol_T * sol_T
    a_C   = 1.0e0_wp - 2.551e-1_wp * sol_C  &
        + 5.151e-2_wp * sol_C * sol_C       &
        - 6.889e-3_wp * sol_C * sol_C * sol_C
    enefunc%cg_dielec_const = e_T * a_C

!   enefunc%cg_debye_length = 1.0e10_wp             &
    debye_length =              &
         sqrt(                                     &
        (CAL2JOU * ELECTRIC_CONST1                  &
        * enefunc%cg_dielec_const * KBOLTZ * sol_T) &
        /                                           &
        (2.0_dp * AVOGADRO1 * AVOGADRO1               &
        * ELEMENT_CHARGE1 * ELEMENT_CHARGE1 * sol_C)  &
        )
    enefunc%cg_debye_length = debye_length

    ! --------------------------------------
    ! Set charges of CG particles in enefunc
    ! --------------------------------------
    !
    n_mols = molecule%num_molecules
    call alloc_enefunc(enefunc, EneFuncCGele, n_mols)

    enefunc%cg_pro_DNA_ele_scale_Q = ene_info%cg_pro_DNA_ele_scale_Q

    ! ---------------------------------
    ! set CG debye-huckel mol-mol pairs
    ! ---------------------------------
    do i = 1, grotop%num_cgelemolpairs

      if (grotop%cg_ele_mol_pairs(i)%is_intermol) then

        do j = grotop%cg_ele_mol_pairs(i)%grp1_start, grotop%cg_ele_mol_pairs(i)%grp1_end
          do k = grotop%cg_ele_mol_pairs(i)%grp2_start, grotop%cg_ele_mol_pairs(i)%grp2_end
            enefunc%cg_ele_mol_pair(j,k) = int(grotop%cg_ele_mol_pairs(i)%func,kind=1) 
            enefunc%cg_ele_mol_pair(k,j) = int(grotop%cg_ele_mol_pairs(i)%func,kind=1)
          end do
        end do

      else

        do j = grotop%cg_ele_mol_pairs(i)%grp1_start, grotop%cg_ele_mol_pairs(i)%grp1_end
          enefunc%cg_ele_mol_pair(j,j) = int(grotop%cg_ele_mol_pairs(i)%func,kind=1)
        end do

      end if

    end do

    ! molecule type, coordinate, charge in elec cell
    !
    do i = 1, ncell+nboundary
      start_i = domain%start_atom(i)
      do ix = 1, ncharge(i)
        base_type = NA_base_type(ix+start_i)
        if (base_type == NABaseTypeDP) then
          charge_type(ix+start_i) = 1
        else if (base_type > NABaseTypeNAMax) then
          charge_type(ix+start_i) = 2
        else
          charge_type(ix+start_i) = 3
        end if
      end do
    end do

    return

  end subroutine setup_enefunc_cg_ele

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cg_PWMcos
  !> @brief        Defines PWMcos parameters for protein-DNA interactions.
  !! @authors      JJ
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cg_PWMcos(ene_info, grotop, molecule, domain, enefunc)

    ! formal arguments
    type(s_ene_info),         intent(in)    :: ene_info
    type(s_grotop),           intent(in)    :: grotop
    type(s_molecule),         intent(in)    :: molecule
    type(s_domain),   target, intent(inout) :: domain
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variables
    integer                     :: n_atoms
    integer                     :: n_pwmcos, pwm_count
    integer                     :: n_mols
    integer                     :: i, j, k, icel, lold, lnew, ioffset
    integer                     :: ix, idx
    integer                     :: alloc_stat
    integer                     :: ncel

    integer,            pointer :: id_g2l(:)
    integer,            pointer :: atom_2_cell(:)
    integer(1),         pointer :: pwmcos_protein(:)
    type(s_grotop_mol), pointer :: gromol

    ncel           =  domain%num_cell_local + domain%num_cell_boundary
    id_g2l         => domain%id_g2l
    atom_2_cell    => domain%atom_2_cell
    pwmcos_protein => domain%pwmcos_protein

    ! -----------------------
    ! setup PWMcos parameters
    ! -----------------------
    !
    enefunc%pwmcos_sigma   = ene_info%cg_PWMcos_sigma
    enefunc%pwmcos_phi     = ene_info%cg_PWMcos_phi * RAD

    ! --------------------------------------------------
    ! setup pwmcos interaction for every protein residue
    ! --------------------------------------------------
    ! count number of pwmcos
    n_pwmcos    = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        do k = 1, gromol%num_pwmcos
          n_pwmcos = n_pwmcos + 1
        end do
      end do
    end do

    if (n_pwmcos > 0) enefunc%cg_pwmcos_calc = .true.

    ! -----------------------------
    ! setup CG PWMcos mol-mol pairs
    ! -----------------------------
    if (n_pwmcos > 0) then

      n_mols = molecule%num_molecules
      allocate(enefunc%pwmcos_mol_pair(n_mols, n_mols), stat=alloc_stat)
      enefunc%pwmcos_mol_pair(:, :) = 0

      do i = 1, grotop%num_pwmcosmolpairs

        do j = grotop%pwmcos_mol_pairs(i)%grp1_start, &
               grotop%pwmcos_mol_pairs(i)%grp1_end
          do k = grotop%pwmcos_mol_pairs(i)%grp2_start, &
                 grotop%pwmcos_mol_pairs(i)%grp2_end
            enefunc%pwmcos_mol_pair(j, k) = &
                 int(grotop%pwmcos_mol_pairs(i)%func,kind=1)
            enefunc%pwmcos_mol_pair(k, j) = &
                 int(grotop%pwmcos_mol_pairs(i)%func,kind=1)
          end do
        end do

      end do
   
    end if

    if (enefunc%cg_pwmcos_calc) then

      enefunc%cg_pwmcos_calc = .true.


      ! set parameters for every pwmcos term...
      n_atoms     = 0
      n_pwmcos    = 0
      ioffset     = 0

      do i = 1, grotop%num_molss
        gromol => grotop%molss(i)%moltype%mol
        do j = 1, grotop%molss(i)%count
          lold = 0
          do k = 1, gromol%num_pwmcos
            lnew = gromol%pwmcos(k)%protein_idx
            idx  = lnew + ioffset
            ix   = id_g2l(idx)
            if (ix > 0 .and. ix <= domain%num_atom_domain) then
              icel = atom_2_cell(ix)
              if ( lnew /= lold ) then
                lold = lnew
                n_pwmcos = n_pwmcos + 1
              end if
            end if
          end do
          ioffset = ioffset + gromol%num_atoms
        end do
      end do

      enefunc%num_pwmcos_domain = n_pwmcos 
#ifdef HAVE_MPI_GENESIS
      call mpi_allreduce(n_pwmcos, MaxPwmCos, 1, mpi_integer, &
                         mpi_max, mpi_comm_country, ierror)
#endif
      if (MaxPwmCos > 0) MaxPwmCos = max(MaxPwmCos*2, min_atom)
      call alloc_enefunc(enefunc, EneFuncPWMcos, MaxPwmCos)

      ioffset     = 0
      n_pwmcos = 0
      do i = 1, grotop%num_molss

        gromol => grotop%molss(i)%moltype%mol

        do j = 1, grotop%molss(i)%count

          lold = 0

          do k = 1, gromol%num_pwmcos

            lnew = gromol%pwmcos(k)%protein_idx
            idx  = lnew + ioffset
            ix   = id_g2l(idx)

            if (ix > 0 .and. ix <= domain%num_atom_domain) then

              icel = atom_2_cell(ix)
              if ( lnew /= lold ) then

                pwm_count = 0
                lold = lnew
                n_pwmcos = n_pwmcos + 1
                enefunc%pwmcos_protein_id(n_pwmcos) = idx
                if (lnew <= 1) then
                  enefunc%pwmcos_protein_id_N(n_pwmcos) = idx
                else
                  enefunc%pwmcos_protein_id_N(n_pwmcos) = idx - 1
                end if
                if (lnew >= gromol%num_atoms) then
                  enefunc%pwmcos_protein_id_C(n_pwmcos) = idx
                else
                  enefunc%pwmcos_protein_id_C(n_pwmcos) = idx + 1
                end if
  
              end if
  
              pwm_count = pwm_count + 1 
              enefunc%pwmcos_r0    (pwm_count, n_pwmcos)  &
                                      = gromol%pwmcos(k)%r0 * 10.0_wp
              enefunc%pwmcos_theta1(pwm_count, n_pwmcos)  & 
                                      = gromol%pwmcos(k)%theta1 * RAD
              enefunc%pwmcos_theta2(pwm_count, n_pwmcos)  &
                                      = gromol%pwmcos(k)%theta2 * RAD
              enefunc%pwmcos_theta3(pwm_count, n_pwmcos)  &
                                    = gromol%pwmcos(k)%theta3 * RAD
              enefunc%pwmcos_ene_A (pwm_count, n_pwmcos)  & 
                                      = gromol%pwmcos(k)%ene_A*KBOLTZ*300.0_wp
              enefunc%pwmcos_ene_C (pwm_count, n_pwmcos)  &
                                      = gromol%pwmcos(k)%ene_C*KBOLTZ*300.0_wp
              enefunc%pwmcos_ene_G (pwm_count, n_pwmcos)  &
                                      = gromol%pwmcos(k)%ene_G*KBOLTZ*300.0_wp
              enefunc%pwmcos_ene_T (pwm_count, n_pwmcos)  &
                                      = gromol%pwmcos(k)%ene_T*KBOLTZ*300.0_wp
              enefunc%pwmcos_gamma (pwm_count, n_pwmcos)  &
                                      = gromol%pwmcos(k)%gamma
              enefunc%pwmcos_eps   (pwm_count, n_pwmcos)  &
                                      = gromol%pwmcos(k)%eps_shift
              enefunc%pwmcos_specificity(pwm_count, n_pwmcos)  &
                                      = gromol%pwmcos(k)%func

              enefunc%pwmcos_count(n_pwmcos) = pwm_count
  
            end if
  
          end do
  
          ioffset = ioffset + gromol%num_atoms
  
        end do
      end do

    end if

    return

  end subroutine setup_enefunc_cg_pwmcos

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cg_PWMcos_lb
  !> @brief        Defines PWMcos parameters for protein-DNA interactions.
  !! @authors      JJ
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cg_PWMcos_lb(grotop, molecule, domain, enefunc)

    ! formal arguments
    type(s_grotop),           intent(in)    :: grotop
    type(s_molecule),         intent(in)    :: molecule
    type(s_domain),   target, intent(inout) :: domain
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variables
    integer                     :: n_atoms
    integer                     :: n_pwmcos, pwm_count
    integer                     :: n_mols
    integer                     :: i, j, k, icel, lold, lnew, ioffset
    integer                     :: ix, idx
    integer                     :: alloc_stat
    integer                     :: ncel

    integer,            pointer :: id_g2l(:)
    integer,            pointer :: atom_2_cell(:)
    integer(1),         pointer :: pwmcos_protein(:)
    type(s_grotop_mol), pointer :: gromol

    ncel           =  domain%num_cell_local + domain%num_cell_boundary
    id_g2l         => domain%id_g2l
    atom_2_cell    => domain%atom_2_cell
    pwmcos_protein => domain%pwmcos_protein

    ! --------------------------------------------------
    ! setup pwmcos interaction for every protein residue
    ! --------------------------------------------------
    ! count number of pwmcos
    n_pwmcos    = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        do k = 1, gromol%num_pwmcos
          n_pwmcos = n_pwmcos + 1
        end do
      end do
    end do

    if (n_pwmcos > 0) enefunc%cg_pwmcos_calc = .true.

    if (enefunc%cg_pwmcos_calc) then

      ! set parameters for every pwmcos term...
      n_atoms     = 0
      n_pwmcos    = 0
      ioffset     = 0

      do i = 1, grotop%num_molss
        gromol => grotop%molss(i)%moltype%mol
        do j = 1, grotop%molss(i)%count
          lold = 0
          do k = 1, gromol%num_pwmcos
            lnew = gromol%pwmcos(k)%protein_idx
            idx  = lnew + ioffset
            ix   = id_g2l(idx)
            if (ix > 0 .and. ix <= domain%num_atom_domain) then
              icel = atom_2_cell(ix)
              if ( lnew /= lold ) then
                lold = lnew
                n_pwmcos = n_pwmcos + 1
              end if
            end if
          end do
          ioffset = ioffset + gromol%num_atoms
        end do
      end do

      enefunc%num_pwmcos_domain = n_pwmcos 
#ifdef HAVE_MPI_GENESIS
      call mpi_allreduce(n_pwmcos, MaxPwmCos, 1, mpi_integer, &
                         mpi_max, mpi_comm_country, ierror)
#endif
      if (MaxPwmcos > 0) MaxPwmCos = max(MaxPwmCos*2,min_atom)
      call alloc_enefunc(enefunc, EneFuncPWMcos, MaxPwmCos)

      ioffset     = 0
      n_pwmcos = 0
      do i = 1, grotop%num_molss

        gromol => grotop%molss(i)%moltype%mol

        do j = 1, grotop%molss(i)%count

          lold = 0

          do k = 1, gromol%num_pwmcos

            lnew = gromol%pwmcos(k)%protein_idx
            idx  = lnew + ioffset
            ix   = id_g2l(idx)

            if (ix > 0 .and. ix <= domain%num_atom_domain) then

              icel = atom_2_cell(ix)
              if ( lnew /= lold ) then

                pwm_count = 0
                lold = lnew
                n_pwmcos = n_pwmcos + 1
                enefunc%pwmcos_protein_id(n_pwmcos) = idx
                if (lnew <= 1) then
                  enefunc%pwmcos_protein_id_N(n_pwmcos) = idx
                else
                  enefunc%pwmcos_protein_id_N(n_pwmcos) = idx - 1
                end if
                if (lnew >= gromol%num_atoms) then
                  enefunc%pwmcos_protein_id_C(n_pwmcos) = idx
                else
                  enefunc%pwmcos_protein_id_C(n_pwmcos) = idx + 1
                end if
  
              end if
  
              pwm_count = pwm_count + 1 
              enefunc%pwmcos_r0    (pwm_count, n_pwmcos)  &
                                      = gromol%pwmcos(k)%r0 * 10.0_wp
              enefunc%pwmcos_theta1(pwm_count, n_pwmcos)  & 
                                      = gromol%pwmcos(k)%theta1 * RAD
              enefunc%pwmcos_theta2(pwm_count, n_pwmcos)  &
                                      = gromol%pwmcos(k)%theta2 * RAD
              enefunc%pwmcos_theta3(pwm_count, n_pwmcos)  &
                                    = gromol%pwmcos(k)%theta3 * RAD
              enefunc%pwmcos_ene_A (pwm_count, n_pwmcos)  & 
                                      = gromol%pwmcos(k)%ene_A*KBOLTZ*300.0_wp
              enefunc%pwmcos_ene_C (pwm_count, n_pwmcos)  &
                                      = gromol%pwmcos(k)%ene_C*KBOLTZ*300.0_wp
              enefunc%pwmcos_ene_G (pwm_count, n_pwmcos)  &
                                      = gromol%pwmcos(k)%ene_G*KBOLTZ*300.0_wp
              enefunc%pwmcos_ene_T (pwm_count, n_pwmcos)  &
                                      = gromol%pwmcos(k)%ene_T*KBOLTZ*300.0_wp
              enefunc%pwmcos_gamma (pwm_count, n_pwmcos)  &
                                      = gromol%pwmcos(k)%gamma
              enefunc%pwmcos_eps   (pwm_count, n_pwmcos)  &
                                      = gromol%pwmcos(k)%eps_shift
              enefunc%pwmcos_specificity(pwm_count, n_pwmcos)  &
                                      = gromol%pwmcos(k)%func

              enefunc%pwmcos_count(n_pwmcos) = pwm_count
  
            end if
  
          end do
  
          ioffset = ioffset + gromol%num_atoms
  
        end do
      end do

    end if

    return

  end subroutine setup_enefunc_cg_PWMcos_lb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cg_PWMcosns
  !> @brief        Defines PWMcosns parameters for protein-DNA interactions.
  !! @authors      JJ
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cg_PWMcosns(grotop, molecule, domain, enefunc)

    ! formal arguments
    type(s_grotop),           intent(in)    :: grotop
    type(s_molecule),         intent(in)    :: molecule
    type(s_domain),   target, intent(inout) :: domain
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variables
    integer                     :: n_atoms
    integer                     :: n_pwmcosns, pwm_count
    integer                     :: n_mols
    integer                     :: i, j, k, icel, lold, lnew, ioffset
    integer                     :: ix, idx
    integer                     :: alloc_stat
    integer                     :: ncel

    integer,            pointer :: id_g2l(:)
    integer,            pointer :: atom_2_cell(:)
    integer(1),         pointer :: pwmcosns_protein(:)
    type(s_grotop_mol), pointer :: gromol

    ncel             =  domain%num_cell_local + domain%num_cell_boundary
    id_g2l           => domain%id_g2l
    atom_2_cell      => domain%atom_2_cell
    pwmcosns_protein => domain%pwmcosns_protein

    ! --------------------------------------------------
    ! setup pwmcosns interaction for every protein residue
    ! --------------------------------------------------
    ! count number of pwmcosns
    n_pwmcosns    = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        do k = 1, gromol%num_pwmcosns
          n_pwmcosns = n_pwmcosns + 1
        end do
      end do
    end do

    ! -----------------------------
    ! setup CG PWMcosns mol-mol pairs
    ! -----------------------------
    if (n_pwmcosns > 0) then

      n_mols = molecule%num_molecules
      allocate(enefunc%pwmcosns_mol_pair(n_mols, n_mols), stat=alloc_stat)
      enefunc%pwmcosns_mol_pair(:, :) = 0

      do i = 1, grotop%num_pwmcosnsmolpairs

        do j = grotop%pwmcosns_mol_pairs(i)%grp1_start, &
               grotop%pwmcosns_mol_pairs(i)%grp1_end
          do k = grotop%pwmcosns_mol_pairs(i)%grp2_start, &
                 grotop%pwmcosns_mol_pairs(i)%grp2_end
            enefunc%pwmcosns_mol_pair(j, k) = &
               int(grotop%pwmcosns_mol_pairs(i)%func,kind=1)
            enefunc%pwmcosns_mol_pair(k, j) = &
               int(grotop%pwmcosns_mol_pairs(i)%func,kind=1)
          end do
        end do

      end do

    end if

    if (n_pwmcosns > 0) enefunc%cg_pwmcosns_calc = .true.

    if (enefunc%cg_pwmcosns_calc) then

      ! set parameters for every pwmcos term...
      n_atoms     = 0
      n_pwmcosns  = 0
      ioffset     = 0

      do i = 1, grotop%num_molss
        gromol => grotop%molss(i)%moltype%mol
        do j = 1, grotop%molss(i)%count
          lold = 0
          do k = 1, gromol%num_pwmcosns
            lnew = gromol%pwmcosns(k)%protein_idx
            idx  = lnew + ioffset
            ix   = id_g2l(idx)
            if (ix > 0 .and. ix <= domain%num_atom_domain) then
              icel = atom_2_cell(ix)
              if ( lnew /= lold ) then
                lold = lnew
                n_pwmcosns = n_pwmcosns + 1
              end if
            end if
          end do
          ioffset = ioffset + gromol%num_atoms
        end do
      end do

      enefunc%num_pwmcosns_domain = n_pwmcosns

#ifdef HAVE_MPI_GENESIS
      call mpi_allreduce(n_pwmcosns, MaxPwmCosns, 1, mpi_integer, &
                         mpi_max, mpi_comm_country, ierror)
#endif
      if (maxPwmCosns > 0) MaxPwmCosns = max(MaxPwmCosns*2, min_atom)
      call alloc_enefunc(enefunc, EneFuncPWMcosns, MaxPwmCosns)

      ioffset = 0
      n_pwmcosns = 0

      do i = 1, grotop%num_molss

        gromol => grotop%molss(i)%moltype%mol

        do j = 1, grotop%molss(i)%count

          lold = 0

          do k = 1, gromol%num_pwmcosns

            lnew = gromol%pwmcosns(k)%protein_idx
            idx  = lnew + ioffset
            ix   = id_g2l(idx)

            if (ix > 0 .and. ix <= domain%num_atom_domain) then

              icel = atom_2_cell(ix)
              if ( lnew /= lold ) then

                pwm_count = 0
                lold = lnew
                n_pwmcosns = n_pwmcosns + 1
                enefunc%pwmcosns_protein_id(n_pwmcosns) = idx
                if (lnew <= 1) then
                  enefunc%pwmcosns_protein_id_N(n_pwmcosns) = idx
                else
                  enefunc%pwmcosns_protein_id_N(n_pwmcosns) = idx - 1
                end if
                if (lnew >= gromol%num_atoms) then
                  enefunc%pwmcosns_protein_id_C(n_pwmcosns) = idx
                else
                  enefunc%pwmcosns_protein_id_C(n_pwmcosns) = idx + 1
                end if
  
              end if
  
              pwm_count = pwm_count + 1 
              enefunc%pwmcosns_r0    (pwm_count, n_pwmcosns)  &
                                      = gromol%pwmcosns(k)%r0 * 10.0_wp
              enefunc%pwmcosns_theta1(pwm_count, n_pwmcosns)  & 
                                      = gromol%pwmcosns(k)%theta1 * RAD
              enefunc%pwmcosns_theta2(pwm_count, n_pwmcosns)  & 
                                      = gromol%pwmcosns(k)%theta2 * RAD
              enefunc%pwmcosns_ene   (pwm_count, n_pwmcosns)  & 
                                      = gromol%pwmcosns(k)%ene
              enefunc%pwmcosns_specificity(pwm_count, n_pwmcosns)  &
                                      = gromol%pwmcosns(k)%func
  
              enefunc%pwmcosns_count(n_pwmcosns) = pwm_count
  
            end if
  
          end do
  
          ioffset = ioffset + gromol%num_atoms
  
        end do
      end do

    end if

    return

  end subroutine setup_enefunc_cg_pwmcosns

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    count_nonb_excl_go
  !> @brief        Excludes 1-2, 1-3 interactions with Go potential.
  !! @authors      JJ
  !! @param[in]    first   : flag for first call or not
  !! @param[inout] domain  : structure of domain
  !! @param[inout] enefunc : structure of enefunc
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine count_nonb_excl_go(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(inout) :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: ncell, i, j, num_atom
    integer                  :: list1, list2
    integer                  :: start_i, start_j, ix, ixx, iy, iyy, inbc

    real(wp),        pointer :: charge(:)
    integer,         pointer :: ncell_local, ncell_boundary, max_atom
    integer,         pointer :: natom(:), start_atom(:)
    integer,         pointer :: near_cells_count(:), near_cells(:,:)
    integer,         pointer :: atom_2_cell(:)
    integer,         pointer :: id_g2l(:)
    integer,         pointer :: num_nonb_excl(:,:)
    integer,         pointer :: num_nonb_excl1(:,:)
    integer,         pointer :: bondlist(:,:)
    integer,         pointer :: contactlist(:,:)
    integer,         pointer :: anglelist(:,:)
    integer,         pointer :: dihelist(:,:)
    integer,         pointer :: stacklist(:,:)
    integer(1),      pointer :: exclusion_mask(:,:)

    ncell_local          => domain%num_cell_local
    ncell_boundary       => domain%num_cell_boundary
    max_atom             => domain%max_num_atom
    natom                => domain%num_atom
    start_atom           => domain%start_atom
    near_cells_count     => domain%near_cells_count
    near_cells           => domain%near_cells      
    atom_2_cell          => domain%atom_2_cell
    atom_2_cell          => domain%atom_2_cell
    id_g2l               => domain%id_g2l
    charge               => domain%charge
    
    bondlist             => enefunc%bond_list
    contactlist          => enefunc%contact_list
    anglelist            => enefunc%angl_list
    dihelist             => enefunc%dihe_list
    stacklist            => enefunc%base_stack_list
    num_nonb_excl        => enefunc%num_nonb_excl
    num_nonb_excl1       => enefunc%num_nonb_excl1
    exclusion_mask       => enefunc%exclusion_mask

    ncell = domain%num_cell_local + domain%num_cell_boundary

    num_atom = domain%num_atom_domain+domain%num_atom_boundary

    call timer(TimerMnonb, TimerOn)
    !$omp parallel do private(i, start_i, ix, ixx, iy, iyy, inbc, j, start_j)
    do i = 1, ncell
      start_i = start_atom(i)
      do ix = 1, natom(i)
        ixx = ix + start_i
        do iy = 1, natom(i)
          iyy = iy + start_i
          exclusion_mask(iyy,ixx) = 0
          exclusion_mask(ixx,iyy) = 0
        end do
      end do
      do inbc = 1, near_cells_count(i)
        j = near_cells(inbc,i)
        start_j = start_atom(j)
        do ix = 1, natom(i)
          ixx = ix + start_i
          do iy = 1, natom(j)
            iyy = iy + start_j
            exclusion_mask(iyy,ixx) = 0
            exclusion_mask(ixx,iyy) = 0
          end do
        end do
      end do
    end do
!   do i = 1, num_atom
!     do j = 1, num_atom
!       exclusion_mask(j,i) = 0
!     end do
!   end do
    !$omp end parallel do
    call timer(TimerMnonb, TimerOff)

    ! initialization
    !
    !$omp parallel do
    do i = 1, num_atom
      exclusion_mask(i,i) = 1
    end do
    !$omp end parallel do

    ! exclude 1-2 interaction
    !
    if (enefunc%excl_level > 0) then

      !$omp parallel do private(list1,list2)
      do i = 1, enefunc%num_bondsq_domain + enefunc%num_bond_domain

        list1      = bondlist(1,i)
        list2      = bondlist(2,i)
        list1      = id_g2l(list1)
        list2      = id_g2l(list2)
        exclusion_mask(list1,list2) = 1 
        exclusion_mask(list2,list1) = 1 

      end do
      !$omp end parallel do
    end if

    ! exclusion lists
    !
    !$omp parallel do private(list1,list2)
    do i = 1, enefunc%num_contact_domain

      list1      = contactlist(1,i)
      list2      = contactlist(2,i)
      list1      = id_g2l(list1)
      list2      = id_g2l(list2)
      exclusion_mask(list1,list2) = 1 
      exclusion_mask(list2,list1) = 1 
    end do
    !$omp end parallel do

    ! exclusion of base stack
    !
    !$omp parallel do private(list1,list2)
    do i = 1, enefunc%num_stack_domain

      list1  = stacklist(2,i)
      list2  = stacklist(3,i)
      list1  = id_g2l(list1)
      list2  = id_g2l(list2)
      exclusion_mask(list1,list2) = 1 
      exclusion_mask(list2,list1) = 1 

    end do
    !$omp end parallel do

    ! exclude 1-3 interaction
    !
    if (enefunc%excl_level > 1) then

      !$omp parallel do private(list1,list2)
      do i = 1, enefunc%num_angle_flexible_domain + &
                enefunc%num_angle_local_domain + enefunc%num_angle_domain

        list1  = anglelist(1,i)
        list2  = anglelist(3,i)
        list1  = id_g2l(list1)
        list2  = id_g2l(list2)
        exclusion_mask(list1,list2) = 1 
        exclusion_mask(list2,list1) = 1 

      end do
      !$omp end parallel do

    end if

    ! count 1-4 interaction
    !
    if (enefunc%excl_level > 2) then

      !$omp parallel do private(list1,list2)
      do i = 1, enefunc%num_dihe_flexible_domain + &
                enefunc%num_dihe_local_domain + enefunc%num_dihe_domain
        list1      = dihelist(1,i)
        list2      = dihelist(4,i)
        list1      = id_g2l(list1)
        list2      = id_g2l(list2)
        exclusion_mask(list1,list2) = 1 
        exclusion_mask(list2,list1) = 1 

      end do
      !$omp end parallel do

    end if

    call mpi_barrier(mpi_comm_country, ierror)

  end subroutine count_nonb_excl_go

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    count_nonb_excl
  !> @brief        Excludes 1-2, 1-3 interactions with Martini potential.
  !! @authors      JJ
  !! @param[in]    first   : flag for first call or not
  !! @param[inout] domain  : structure of domain
  !! @param[inout] enefunc : structure of enefunc
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine count_nonb_excl(first, domain, enefunc)

    ! formal arguments
    logical,                 intent(in   ) :: first
    type(s_domain),  target, intent(inout) :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: ncell, i, j, num_atom
    integer                  :: list1, list2
    integer                  :: start_i, start_j, ix, ixx, iy, iyy, inbc

    real(wp),        pointer :: charge(:)
    integer,         pointer :: ncell_local, ncell_boundary, max_atom
    integer,         pointer :: natom(:), start_atom(:)
    integer,         pointer :: near_cells_count(:), near_cells(:,:)
    integer,         pointer :: atom_2_cell(:)
    integer,         pointer :: id_g2l(:)
    integer,         pointer :: num_nonb_excl(:,:)
    integer,         pointer :: num_nonb_excl1(:,:)
    integer,         pointer :: bondlist(:,:)
    integer,         pointer :: anglelist(:,:)
    integer,         pointer :: dihelist(:,:)
    integer,         pointer :: excllist(:,:)
    integer(1),      pointer :: exclusion_mask(:,:)

    ncell_local          => domain%num_cell_local
    ncell_boundary       => domain%num_cell_boundary
    max_atom             => domain%max_num_atom
    natom                => domain%num_atom
    start_atom           => domain%start_atom
    near_cells_count     => domain%near_cells_count
    near_cells           => domain%near_cells      
    atom_2_cell          => domain%atom_2_cell
    id_g2l               => domain%id_g2l
    charge               => domain%charge
    
    bondlist             => enefunc%bond_list
    anglelist            => enefunc%angl_list
    dihelist             => enefunc%dihe_list
    excllist             => enefunc%excl_list
    num_nonb_excl        => enefunc%num_nonb_excl
    num_nonb_excl1       => enefunc%num_nonb_excl1
    exclusion_mask       => enefunc%exclusion_mask

    ncell = domain%num_cell_local + domain%num_cell_boundary

    num_atom = domain%num_atom_domain+domain%num_atom_boundary

    call timer(TimerMnonb, TimerOn)
    !$omp parallel do private(i, start_i, ix, ixx, iy, iyy, inbc, j, start_j)
    do i = 1, ncell
      start_i = start_atom(i)
      do ix = 1, natom(i)
        ixx = ix + start_i
        do iy = 1, natom(i)
          iyy = iy + start_i
          exclusion_mask(iyy,ixx) = 0
          exclusion_mask(ixx,iyy) = 0
        end do
      end do
      do inbc = 1, near_cells_count(i)
        j = near_cells(inbc,i)
        start_j = start_atom(j)
        do ix = 1, natom(i)
          ixx = ix + start_i
          do iy = 1, natom(j)
            iyy = iy + start_j
            exclusion_mask(iyy,ixx) = 0
            exclusion_mask(ixx,iyy) = 0
          end do
        end do
      end do
    end do
    !$omp end parallel do
    call timer(TimerMnonb, TimerOff)

    ! initialization
    !
    !$omp parallel do
    do i = 1, num_atom
      exclusion_mask(i,i) = 1
    end do
    !$omp end parallel do

    ! exclude 1-2 interaction
    !
    if (enefunc%excl_level > 0) then

      !$omp parallel do private(list1,list2)
      do i = 1, enefunc%num_bondsq_domain
        list1      = bondlist(1,i)
        list2      = bondlist(2,i)
        list1      = id_g2l(list1)
        list2      = id_g2l(list2)
        exclusion_mask(list1,list2) = 1 
        exclusion_mask(list2,list1) = 1 
        excllist(1,i) = list1
        excllist(2,i) = list2
      end do
      !$omp end parallel do
      ix = enefunc%num_bondsq_domain
    end if

    ! exclude 1-3 interaction
    !
    if (enefunc%excl_level > 1) then

      !$omp parallel do private(list1,list2)
      do i = 1, enefunc%num_angle_domain
        list1  = anglelist(1,i)
        list2  = anglelist(3,i)
        list1  = id_g2l(list1)
        list2  = id_g2l(list2)
        exclusion_mask(list1,list2) = 1 
        exclusion_mask(list2,list1) = 1 
        excllist(1,i+ix) = list1
        excllist(2,i+ix) = list2
      end do
      !$omp end parallel do
      ix = ix+enefunc%num_angle_domain

    end if

    ! count 1-4 interaction
    !
    if (enefunc%excl_level > 2) then

      ix = ix+enefunc%num_angle_domain
      !$omp parallel do private(list1,list2)
      do i = 1, enefunc%num_dihe_domain
        list1      = dihelist(1,i)
        list2      = dihelist(4,i)
        list1      = id_g2l(list1)
        list2      = id_g2l(list2)
        exclusion_mask(list1,list2) = 1 
        exclusion_mask(list2,list1) = 1 
        excllist(1,i+ix) = list1
        excllist(2,i+ix) = list2
      end do
      !$omp end parallel do
      ix = ix + enefunc%num_dihe_domain

    end if

    enefunc%num_excl = ix 

    call mpi_barrier(mpi_comm_country, ierror)

  end subroutine count_nonb_excl


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_contact
  !> @brief        Defines Contact term for each cell.
  !! @authors      JJ
  !! @param[in]    grotop   : GROMACS TOP information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_contact(grotop, domain, enefunc, comm)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc
    type(s_comm),    target, intent(inout) :: comm   

    ! local variables
    real(wp)                 :: sig, eps
    integer                  :: step
    integer                  :: i, j, k, c_idx
    integer                  :: ioffset, ncontact
    integer                  :: idx1, idx2, icel1, icel2, icel_local
    integer                  :: ncel, ncelb

    type(s_grotop_mol), pointer :: gromol
    real(wp),           pointer :: lj12(:), lj10(:), lj6(:)
    integer,            pointer :: list(:,:), func(:)
    integer,            pointer :: id_g2l(:), atom_2_cell(:)
    real(wp),           pointer :: natom(:)

    id_g2l        => domain%id_g2l
    atom_2_cell   => domain%atom_2_cell
    natom         => domain%num_atom_t0

    list          => enefunc%contact_list
    func          => enefunc%contact_func
    lj12          => enefunc%contact_lj12
    lj10          => enefunc%contact_lj10
    lj6           => enefunc%contact_lj6

    ncel          =  domain%num_cell_local
    ncelb         =  domain%num_cell_local + domain%num_cell_boundary

    do step = 1, 2

      ioffset   = 0
      ncontact  = 0

      do i = 1, grotop%num_molss
        gromol => grotop%molss(i)%moltype%mol
        do j = 1, grotop%molss(i)%count

          do k = 1, gromol%num_pairs

            idx1 = gromol%pairs(k)%atom_idx1 + ioffset
            idx2 = gromol%pairs(k)%atom_idx2 + ioffset

            icel1 = id_g2l(idx1)
            icel2 = id_g2l(idx2)

            if (icel1 /= 0 .and. icel2 /= 0) then

              icel1 = atom_2_cell(icel1)
              icel2 = atom_2_cell(icel2)
              if (icel1 == icel2) then
                icel_local = icel1
              else if (natom(icel1) >= natom(icel2)) then
                icel_local = icel2
              else if (natom(icel1) < natom(icel2)) then
                icel_local = icel1
              end if

              if (icel_local > 0 .and. icel_local <= ncel) then

                ncontact = ncontact + 1

                if (step == 2) then

                  c_idx = ncontact
                  enefunc%contact_list(1,c_idx) = idx1
                  enefunc%contact_list(2,c_idx) = idx2

                  if (grotop%defaults%combi_rule == 2) then
                    sig = 10.0_wp * gromol%pairs(k)%v
                    eps = JOU2CAL * gromol%pairs(k)%w

                    if (enefunc%forcefield == ForcefieldKBGO) then
                      enefunc%contact_lj12(c_idx) = 13.0_wp*eps*(sig**12)
                      enefunc%contact_lj10(c_idx) = 18.0_wp*eps*(sig**10)
                      enefunc%contact_lj6 (c_idx) =  4.0_wp*eps*(sig**6)
                    else if (enefunc%forcefield == forcefieldCAGO .or. &
                             enefunc%forcefield == forcefieldRESIDCG) then
                      enefunc%contact_lj12(c_idx) = 5.0_wp*eps*(sig**12)
                      enefunc%contact_lj10(c_idx) = 6.0_wp*eps*(sig**10)
                      enefunc%contact_lj6 (c_idx) = 0.0_wp
                    else
                      enefunc%contact_lj12(c_idx) = eps * (sig**12)
                      enefunc%contact_lj10(c_idx) = 0.0_wp
                      enefunc%contact_lj6 (c_idx) = 2.0_wp * eps * (sig** 6)
                    end if

                  else if ((gromol%pairs(k)%func == 10) .or. &
                           (gromol%pairs(k)%func == 11) .or. &
                           (gromol%pairs(k)%func == 14) .or. &
                           (gromol%pairs(k)%func == 15)) then
                    enefunc%contact_func(c_idx) = gromol%pairs(k)%func
                    enefunc%contact_lj12(c_idx) = JOU2CAL * gromol%pairs(k)%w 
                    enefunc%contact_lj10(c_idx) = 10.0_wp * gromol%pairs(k)%v 
                    enefunc%contact_lj6 (c_idx) = JOU2CAL * gromol%pairs(k)%khh/100.0_wp 

                  else if ((gromol%pairs(k)%func == 10) .or. &
                           (gromol%pairs(k)%func == 15)) then
                    call error_msg('Error: Setup_Enefunc_Nonb - Functions 12 '//&
			           'and 13 are not available.')

                  else if (gromol%pairs(k)%func == 21) then
                    sig = 10.0_wp * gromol%pairs(k)%r0
                    eps = JOU2CAL * gromol%pairs(k)%khh
                    enefunc%contact_lj12(c_idx) = eps*(sig**12)
                    enefunc%contact_lj10(c_idx) = eps*(sig**10)

                  else
                    
                    if (enefunc%forcefield == forcefieldCAGO) then
                      enefunc%contact_lj12(c_idx) = &
                                         gromol%pairs(k)%w*1.0E12_wp*JOU2CAL
                      enefunc%contact_lj10(c_idx) = &
                                         gromol%pairs(k)%v*1.0E10_wp*JOU2CAL
                      enefunc%contact_lj6 (c_idx) = 0.0_wp
                    else
                      enefunc%contact_lj12(c_idx) = &
                                         gromol%pairs(k)%w*1.0E12_wp*JOU2CAL
                      enefunc%contact_lj10(c_idx) = 0.0_wp
                      enefunc%contact_lj6 (c_idx) = &
                                         gromol%pairs(k)%v*1.0E6_wp*JOU2CAL
                    end if
                  end if

                end if

              end if
            end if

          end do

          ioffset = ioffset + gromol%num_atoms

        end do
      end do

      if (step == 1) then

        enefunc%num_contact_domain = ncontact
        enefunc%num_contact_boundary = 0

        k = enefunc%num_contact_domain * domain%num_comm_proc

        call mpi_allreduce(k, MaxContact, 1, mpi_integer, &
                           mpi_max, mpi_comm_country, ierror)
        if (MaxContact > 0) MaxContact = max(MaxContact*2,min_atom)

        call alloc_enefunc(enefunc, EneFuncContact, MaxContact, 1)

      end if

    end do
    call mpi_allreduce(enefunc%num_contact_domain, enefunc%num_contact_all, &
                       1, mpi_integer, mpi_sum, mpi_comm_country, ierror)

    return

  end subroutine setup_enefunc_contact

end module cg_enefunc_gromacs_mod
