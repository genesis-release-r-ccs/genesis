!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_enefunc_restraints_mod
!> @brief   restraint energy functions
!! @authors Chigusa Kobayashi (CK), Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8
  
#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_enefunc_restraints_mod

  use sp_restraints_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use molecules_str_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: setup_enefunc_restraints
  public  :: setup_enefunc_restraints_pio
  public  :: setup_restraint_mode
  private :: setup_enefunc_rest_group
  private :: setup_enefunc_rest_func
  private :: setup_enefunc_rest_domain
  public  :: setup_restraints_grouplist
  public  :: setup_restraints_exponent_func
  public  :: setup_restraints_constants
  public  :: setup_restraints_reference
  public  :: setup_restraints_exponent_dist
  public  :: setup_restraints_weight_dist
  private :: setup_enefunc_rest_refcoord
  private :: setup_enefunc_rest_mode

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_restraints
  !> @brief        define restraint potential
  !! @authors      CK, TM
  !! @param[in]    molecule   : molecule information
  !! @param[in]    restraints : restraints information
  !! @param[in]    domain     : domain information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_restraints(molecule, restraints, domain, enefunc)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_restraints),      intent(in)    :: restraints
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: iatm, ires, oldres
    character(4)             :: segname, oldsegnm
    character(6)             :: resname, oldresnm
    character(4)             :: atname


    enefunc%restraint         = restraints%restraint_flag
    enefunc%pressure_position = restraints%pressure_position
    enefunc%pressure_rmsd     = restraints%pressure_rmsd  

    if (.not. restraints%restraint_flag) &
      return

    ! setup group
    !
    call setup_enefunc_rest_group(molecule, restraints, enefunc)

    ! setup func
    !
    call setup_enefunc_rest_func(molecule, restraints, enefunc)

    ! setup domain
    !
    call setup_enefunc_rest_domain(domain, enefunc)

    ! we do not allow duplication option with restraint
    !
    if (domain%num_duplicate > 1) &
      call error_msg('Setup_Enefunc_Restraint> duplication is not allowed with restraint')

    ! summary of setup enefunc_restraint
    !
    if (main_rank) then

      write(MsgOut,'(A)') 'Setup_Enefunc_Restraints> Setup restraint groups'

      do i = 1, enefunc%num_restraintgroups
        write(MsgOut,'(A,I5,3A)') ' group = ', i, ', "', trim(restraints%group(i)),'"'
        write(MsgOut,'(A,I5)')   ' # of atoms = ', enefunc%restraint_numatoms(i)
        write(MsgOut,'(A)')      ' atomlist: ' 
        do j = 1, enefunc%restraint_numatoms(i) 
          write(MsgOut,'(i7,$)') enefunc%restraint_atomlist(j,i)
          if (mod(j,10) == 0 .and. j /= enefunc%restraint_numatoms(i)) then
            write(MsgOut,'(A)') ''
          end if
        end do
        write(MsgOut,'(A)') ''
      end do
      write(MsgOut,'(A)') ''

      write(MsgOut,'(A)') 'Setup_Enefunc_Restraints> Setup restraint functions'

      do i = 1, enefunc%num_restraintfuncs

        write(MsgOut,'(A,I5,A,I5)') &
          ' func  = ', i, ' kind  = ', enefunc%restraint_kind(i)

        ! summary for positional restraint
        !
        if (enefunc%restraint_kind(i) == RestraintsFuncPOSI) then
          write(MsgOut,'(A,4F8.3)') &
            ' const(total, x, y, z) = ', enefunc%restraint_const(1:4,i) 
          write(MsgOut,'(A,I5)') ' exponent of function = ',              &
                                 enefunc%restraint_exponent_func(i)
        else if (enefunc%restraint_kind(i) == RestraintsFuncRMSD .or.     &
                 enefunc%restraint_kind(i) == RestraintsFuncRMSDCOM) then
          write(MsgOut,'(A,F8.3)') &
            ' force constant        = ', enefunc%restraint_const(1,i)
          write(MsgOut,'(A,F8.3)') &
            ' target RMSD           = ', enefunc%restraint_ref(1,i)
          write(MsgOut,'(A,I5)') ' exponent of function = ',              &
                                 enefunc%restraint_exponent_func(i)
        else
          write(MsgOut,'(A,F8.3)') &
            ' const                 = ',   enefunc%restraint_const(1,i)
          write(MsgOut,'(A,F8.3)') &
            ' ref                   = ',   enefunc%restraint_ref(1,i)
          write(MsgOut,'(A,I5)') ' exponent of function = ',              &
                                 enefunc%restraint_exponent_func(i)
        end if
 
        ! summary for distance restraint
        !
        if (enefunc%restraint_kind(i) == RestraintsFuncDIST .or.          &
            enefunc%restraint_kind(i) == RestraintsFuncDISTCOM) then

          write(MsgOut,'(A,I5)') ' # of distances  = ', enefunc%restraint_funcgrp(i)/2
          write(MsgOut,'(A)')    ' distancelist: '
          do j = 1, enefunc%restraint_funcgrp(i)/2
            write(MsgOut,'(" group = (",i3,",",i3,")"$)')                 &
                           enefunc%restraint_grouplist(2*j-1,i),          &
                           enefunc%restraint_grouplist(2*j,i)
            write(MsgOut,'(" weight = ",f8.3," exponent = ",i5,$)')       &
                           enefunc%restraint_weight_dist(j,i),            &
                           enefunc%restraint_exponent_dist(j,i)
            write(MsgOut,'(A)') ''
          end do

        ! summary for angle restraint
        !
        else if (enefunc%restraint_kind(i) == RestraintsFuncANGLE .or.    &
                 enefunc%restraint_kind(i) == RestraintsFuncANGLECOM ) then

          write(MsgOut,'(A,I5)') ' # of angles     = ', enefunc%restraint_funcgrp(i)/3
          write(MsgOut,'(A)')    ' anglelist   : '
          do j = 1, enefunc%restraint_funcgrp(i)/3
            write(MsgOut,'(" group = (",i3,",",i3,",",i3,")"$)')          &
                           enefunc%restraint_grouplist(3*j-2,i),          &
                           enefunc%restraint_grouplist(3*j-1,i),          &
                           enefunc%restraint_grouplist(3*j,i)
            write(MsgOut,'(A)') ''
          end do

        ! summary for dihedral angle restraint
        !
        else if (enefunc%restraint_kind(i) == RestraintsFuncDIHED  .or.   &
                 enefunc%restraint_kind(i) == RestraintsFuncDIHEDCOM ) then

          write(MsgOut,'(A,I5)') ' # of dihedrals  = ', enefunc%restraint_funcgrp(i)/4
          write(MsgOut,'(A)')    ' dihedrallist: ' 
          do j = 1, enefunc%restraint_funcgrp(i)/4
            write(MsgOut,'(" group = (",i3,",",i3,",",i3,",",i3,")"$)')   &
                           enefunc%restraint_grouplist(4*j-3,i),          &
                           enefunc%restraint_grouplist(4*j-2,i),          &
                           enefunc%restraint_grouplist(4*j-1,i),          &
                           enefunc%restraint_grouplist(4*j,i)
            write(MsgOut,'(A)') ''
          end do

        ! summary for repulsive restraint
        !
        else if (enefunc%restraint_kind(i) == RestraintsFuncREPUL  .or.   &
                 enefunc%restraint_kind(i) == RestraintsFuncREPULCOM ) then

          write(MsgOut,'(A,I5)') ' # of distances  = ',                   &
            enefunc%restraint_funcgrp(i)*(enefunc%restraint_funcgrp(i)-1)/2
          write(MsgOut,'(A)')    ' distancelist: '
          do j = 1, enefunc%restraint_funcgrp(i)
            do k = j+1, enefunc%restraint_funcgrp(i)
              write(MsgOut,'(" group = (",i3,",",i3,")"$)')               &
                             enefunc%restraint_grouplist(j,i),            &
                             enefunc%restraint_grouplist(k,i)
              write(MsgOut,'(A)') ''
            end do
          end do

        ! summary for flat bottom restraint
        !
        else if (enefunc%restraint_kind(i) == RestraintsFuncFB  .or.      &
                 enefunc%restraint_kind(i) == RestraintsFuncFBCOM ) then

          write(MsgOut,'(A,I5)') ' # of distances  = ',                   &
            enefunc%restraint_funcgrp(i)*(enefunc%restraint_funcgrp(i)-1)/2
          write(MsgOut,'(A)')    ' distancelist: '
          do j = 1, enefunc%restraint_funcgrp(i)
            do k = j+1, enefunc%restraint_funcgrp(i)
              write(MsgOut,'(" group = (",i3,",",i3,")"$)')               &
                             enefunc%restraint_grouplist(j,i),            &
                             enefunc%restraint_grouplist(k,i)
              write(MsgOut,'(A)') ''
            end do
          end do

        ! summary for positinal and RMSD restraint
        !
        else 

          write(MsgOut,'(A,I5)') ' # of groups  = ', enefunc%restraint_funcgrp(i)
          write(MsgOut,'(" grouplist: ",$)') 
          do j = 1, enefunc%restraint_funcgrp(i) 
            write(MsgOut,'(i3,$)') enefunc%restraint_grouplist(j,i)
            if (mod(j,20) == 0 .and. j /= enefunc%restraint_funcgrp(i) )  &
              write(MsgOut,'(A)') ''
          end do
          write(MsgOut,'(A)') ''

        end if
        write(MsgOut,'(A)') ''

      end do
    end if

    ! error message when both position and rmsd restraints are assigned
    !
    if (enefunc%restraint_posi .and. enefunc%restraint_rmsd) &
      call error_msg('Setup_Enefunc_Restraint> Choose either position or rmsd restraint')

    return

  end subroutine setup_enefunc_restraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_restraints_pio
  !> @brief        setup restraint information about domain
  !! @authors      JJ
  !! @param[in]    restraints : restraints information
  !! @param[in]    domain     : domain information
  !! @param[inout] enefunc    : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_restraints_pio(restraints, domain, enefunc)

    ! formal arguments
    type(s_restraints),  target, intent(in)    :: restraints
    type(s_domain),      target, intent(in)    :: domain
    type(s_enefunc),     target, intent(inout) :: enefunc

    ! local variable
    integer                  :: i, ig, j, k, ix, ic, iatm, icel
    integer                  :: file_num
    integer                  :: max_grp
    logical                  :: position_flag
    logical                  :: rmsd_flag
    integer                  :: num_indexes(1)

    real(wp),        pointer :: restraint_force_pio(:,:,:,:)
    real(wp),        pointer :: restraint_coord_pio(:,:,:,:)
    integer,         pointer :: ncell, ncell_pio
    integer,         pointer :: cell_l2g_pio(:,:)
    integer(int2),   pointer :: cell_g2l(:)
    integer,         pointer :: nrestraint(:), nrestraint_pio(:,:)
    integer,         pointer :: restraint_atom_pio(:,:,:)


    ncell               => domain%num_cell_local
    ncell_pio           => domain%ncell_local_pio
    cell_l2g_pio        => domain%cell_l2g_pio
    cell_g2l            => domain%cell_g2l

    nrestraint          => enefunc%num_restraint
    nrestraint_pio      => enefunc%num_restraint_pio
    restraint_atom_pio  => enefunc%restraint_atom_pio
    restraint_force_pio => enefunc%restraint_force_pio
    restraint_coord_pio => enefunc%restraint_coord_pio

    k = 0
    do i = 1, restraints%num_groups
      k = k + restraints%num_atoms(i)
    end do
    enefunc%num_atoms_bonds_restraint = k

    ! memory allocation for restraint
    !
    if (restraints%nfunctions > 0) then
      call alloc_enefunc(enefunc, EneFuncReff, 1, 0)
      call alloc_enefunc(enefunc, EneFuncRest, ncell, k, 0)
    end if

    position_flag=.false.

    ! position/rmsd/pc restraint
    !
    if (restraints%nfunctions > 0) then
      enefunc%num_restraintfuncs = restraints%nfunctions
      do i = 1, enefunc%num_restraintfuncs
        enefunc%restraint_kind(i) = restraints%function(i)
      end do
    end if

    if (restraints%nfunctions > 0) then
      if (enefunc%restraint_kind(1) == RestraintsFuncPOSI) then

        enefunc%restraint_posi = .true.
        position_flag=.true.

        call setup_restraints_constants(1, restraints, enefunc)

        do file_num = 1, domain%file_tot_num
  
          do icel = 1, ncell_pio
  
            ic = cell_l2g_pio(icel,file_num)
            i  = cell_g2l(ic)
  
            if (i /= 0) then
  
              nrestraint(i) = nrestraint_pio(icel,file_num)
  
              do ix = 1, nrestraint(i)
  
                enefunc%restraint_atom(ix,i)      = &
                        restraint_atom_pio(ix,icel,file_num)
                enefunc%restraint_force(1:4,ix,i) = &
                        enefunc%restraint_const(1:4,1)
                enefunc%restraint_coord(1:3,ix,i) = &
                        restraint_coord_pio(1:3,ix,icel,file_num)
    
              end do
    
            end if
  
          end do
        end do

      end if
    end if

    if (position_flag) then
      call alloc_enefunc(enefunc, EneFuncRestDomain, ncell)
    end if

    return

  end subroutine setup_enefunc_restraints_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_rest_group
  !> @brief        define restraint atom group
  !! @authors      CK, TM
  !! @param[in]    molecule   : molecule information
  !! @param[in]    restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_rest_group(molecule, restraints, enefunc)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_restraints),      intent(in)    :: restraints
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    real(wip)                :: totmass, coef
    integer                  :: num_groups, max_atoms
    integer                  :: i, j, jatom


    enefunc%num_restraintgroups    = restraints%num_groups
    enefunc%max_restraint_numatoms = restraints%max_atoms
    num_groups = enefunc%num_restraintgroups
    max_atoms  = enefunc%max_restraint_numatoms

    call alloc_enefunc(enefunc, EneFuncRefg, num_groups, max_atoms)

    ! setup atom lists and num atoms
    !
    do j = 1, num_groups
      enefunc%restraint_numatoms(j) = restraints%num_atoms(j)
      do i = 1, max_atoms
        enefunc%restraint_atomlist(i,j) = restraints%atomlist(i,j)
      end do
    end do

    ! check number of atoms in groups
    !
    do i = 1, enefunc%num_restraintgroups
      if (enefunc%restraint_numatoms(i) == 0) &
        call error_msg('Setup_Enefunc_Rest_Group> no atom in group')
    end do

    ! setup mass_coefficient
    !
    enefunc%restraint_masscoef(1:max_atoms,1:num_groups) = 0.0_wip
    do i = 1, enefunc%num_restraintgroups
      totmass = 0.0_wip
      do j = 1, enefunc%restraint_numatoms(i) 
        jatom = enefunc%restraint_atomlist(j,i)
        totmass = totmass + molecule%mass(jatom)
      end do
      coef = 1.0_wip/totmass

      do j = 1, enefunc%restraint_numatoms(i) 
        jatom = enefunc%restraint_atomlist(j,i)
        enefunc%restraint_masscoef(j,i) = molecule%mass(jatom) * coef
      end do
    end do

    return

  end subroutine setup_enefunc_rest_group

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_rest_func
  !> @brief        define restraint function
  !! @authors      CK, TM
  !! @param[in]    molecule   : molecule information
  !! @param[in]    restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_rest_func(molecule, restraints, enefunc)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_restraints),      intent(in)    :: restraints
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, ndata
    integer                  :: num_funcs, max_grp
    integer,     allocatable :: num_indexes(:)


    ! get total number of groups in index
    !
    enefunc%num_restraintfuncs = restraints%nfunctions
    num_funcs = enefunc%num_restraintfuncs
    allocate(num_indexes(1:num_funcs))

    do i = 1, num_funcs
      ndata = split_num(restraints%select_index(i))

      select case(restraints%function(i))
      case(RestraintsFuncDIST:RestraintsFuncDISTCOM)
        if (ndata <  2 .or. mod(ndata,2) /= 0) &
          call error_msg('Setup_Enefunc_Rest_Func> Error in index for DIST')

      case(RestraintsFuncANGLE:RestraintsFuncANGLECOM)
        if (ndata /= 3) &
          call error_msg('Setup_Enefunc_Rest_Func> Error in index for ANGLE')

      case(RestraintsFuncDIHED:RestraintsFuncDIHEDCOM)
        if (ndata /= 4) &
          call error_msg('Setup_Enefunc_Rest_Func> Error in index for DIHE')

      case(RestraintsFuncREPUL:RestraintsFuncREPULCOM)
        if (ndata <  2) &
          call error_msg('Setup_Enefunc_Rest_Func> Error in index for REPUL')

      case(RestraintsFuncFB:RestraintsFuncFBCOM)
        if (ndata <  2) &
          call error_msg('Setup_Enefunc_Rest_Func> Error in index for FB')

      case(RestraintsFuncPC:RestraintsFuncPCCOM)
        if (ndata == 0) &
          call error_msg('Setup_Enefunc_Rest_Func> Error in index for PC')

      case default
        if (ndata /= 1) & 
          call error_msg('Setup_Enefunc_Rest_Func> Error in index for POSI/RMSD')

      end select

      num_indexes(i) = ndata
    end do

    ! get maximum number of groups among all 'index's
    !
    max_grp = maxval(num_indexes(1:num_funcs))
    enefunc%max_restraint_numgrps = max_grp

    call alloc_enefunc(enefunc, EneFuncReff, num_funcs, max_grp)
 
    ! setup parameters
    !
    do i = 1, num_funcs

      enefunc%restraint_funcgrp(i) = num_indexes(i)
      enefunc%restraint_kind(i)    = restraints%function(i)

      ! setup group list
      !
      call setup_restraints_grouplist(i, restraints, enefunc)

      ! setup exponent_func
      !
      call setup_restraints_exponent_func(i, restraints, enefunc)

      ! setup force constants
      !
      call setup_restraints_constants(i, restraints, enefunc)

      ! setup reference   (except for POSI)
      !
      call setup_restraints_reference(i, restraints, enefunc)

      ! setup exponent_dist (for DIST only)
      !
      call setup_restraints_exponent_dist(i, restraints, enefunc)

      ! setup weight_dist   (for DIST only)
      !
      call setup_restraints_weight_dist(i, restraints, enefunc)

      ! setup mode   (for PC only)
      !
      call setup_restraint_mode(i, restraints, enefunc)

    end do

    ! setup reference coordinates (for POSI and RMSD)
    !
    do i = 1, num_funcs
      if (enefunc%restraint_kind(i) == RestraintsFuncPOSI .or.    &
          enefunc%restraint_kind(i) == RestraintsFuncRMSD .or.    &
          enefunc%restraint_kind(i) == RestraintsFuncRMSDCOM .or. &
          enefunc%restraint_kind(i) == RestraintsFuncPC .or.      &
          enefunc%restraint_kind(i) == RestraintsFuncPCCOM) then
        call setup_enefunc_rest_refcoord(molecule, enefunc)
        if (enefunc%restraint_kind(i) == RestraintsFuncPOSI) &
          enefunc%restraint_posi = .true.
        if (enefunc%restraint_kind(i) == RestraintsFuncRMSD .or.  &
            enefunc%restraint_kind(i) == RestraintsFuncRMSDCOM)   &
          enefunc%restraint_rmsd = .true.
        if (enefunc%restraint_kind(i) == RestraintsFuncPC)   &
          enefunc%restraint_pc = .true.
      else if (enefunc%restraint_kind(i) == RestraintsFuncEM) then
        enefunc%restraint_emfit = .true.
      end if
    end do

    ! setup principal component mode (for PC)
    !
    do i = 1, num_funcs
      if (enefunc%restraint_kind(i) == RestraintsFuncPC     .or. &
          enefunc%restraint_kind(i) == RestraintsFuncPCCOM) then
        call setup_enefunc_rest_mode(molecule, enefunc)
        !exit
      end if
    end do

    ! deallocate local array
    !
    deallocate(num_indexes)

    return

  end subroutine setup_enefunc_rest_func

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_rest_domain
  !> @brief        setup restraint information about domain
  !! @authors      CK
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_rest_domain(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    integer                  :: i, ig, j, k, ix, iatm, icel
    logical                  :: position_flag
    logical                  :: rmsd_flag

    real(wp),        pointer :: const(:,:), ref(:,:), coord(:,:)
    real(wp),        pointer :: restraint_force(:,:,:)
    real(wp),        pointer :: restraint_coord(:,:,:)
    integer,         pointer :: ncell
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: num_funcs, num_atoms(:)
    integer,         pointer :: grouplist(:,:), atomlist(:,:), bondslist(:,:)
    integer,         pointer :: restraint_g2pc(:)
    integer,         pointer :: nrestraint(:), restraint_atom(:,:)


    ncell           => domain%num_cell_local
    id_g2l          => domain%id_g2l

    num_funcs       => enefunc%num_restraintfuncs
    num_atoms       => enefunc%restraint_numatoms
    grouplist       => enefunc%restraint_grouplist
    atomlist        => enefunc%restraint_atomlist
    bondslist       => enefunc%restraint_bondslist
    const           => enefunc%restraint_const
    ref             => enefunc%restraint_ref
    coord           => enefunc%restraint_refcoord
    nrestraint      => enefunc%num_restraint
    restraint_g2pc  => enefunc%restraint_g2pc

    position_flag   = .false.
    rmsd_flag       = .false.

    ! check the size for distance restraint
    !
    k = 0
    do i = 1, num_funcs
      if (enefunc%restraint_kind(i) == RestraintsFuncDIST     .or. &
          enefunc%restraint_kind(i) == RestraintsFuncDISTCOM  .or. &
          enefunc%restraint_kind(i) == RestraintsFuncANGLE    .or. &
          enefunc%restraint_kind(i) == RestraintsFuncANGLECOM .or. &
          enefunc%restraint_kind(i) == RestraintsFuncDIHED    .or. &
          enefunc%restraint_kind(i) == RestraintsFuncDIHEDCOM .or. &
          enefunc%restraint_kind(i) == RestraintsFuncREPUL    .or. &
          enefunc%restraint_kind(i) == RestraintsFuncREPULCOM .or. &
          enefunc%restraint_kind(i) == RestraintsFuncFB       .or. &
          enefunc%restraint_kind(i) == RestraintsFuncFBCOM) then
        do ig = 1, enefunc%restraint_funcgrp(i)
          do ix = 1, num_atoms(grouplist(ig,i))
            k = k + 1
          end do
        end do
      end if
    end do
    enefunc%num_atoms_bonds_restraint = k

    ! memory allocation for restraint
    !
    call alloc_enefunc(enefunc, EneFuncRest, ncell, k, j)

    restraint_atom  => enefunc%restraint_atom
    restraint_force => enefunc%restraint_force
    restraint_coord => enefunc%restraint_coord

    ! position/rmsd/pc restraint
    !
    do i = 1, num_funcs

      if (enefunc%restraint_kind(i) == RestraintsFuncPOSI) then
         position_flag=.true.

        do ix = 1, num_atoms(grouplist(1,i))

          iatm = atomlist(ix,grouplist(1,i))
          icel = id_g2l(1,iatm)

          if (icel > 0 .and. icel <= ncell) then
            nrestraint(icel) = nrestraint(icel) + 1
            restraint_atom(nrestraint(icel),icel) = iatm
            restraint_force(1:4,nrestraint(icel),icel) = const(1:4,i)
            restraint_coord(1:3,nrestraint(icel),icel) = coord(1:3,iatm)
          end if

        end do

      else if (enefunc%restraint_kind(i) == REstraintsFuncRMSD .or.     &
               enefunc%restraint_kind(i) == REstraintsFuncRMSDCOM) then
         rmsd_flag=.true.

        k = 0
        do ix = 1, num_atoms(grouplist(1,i))

          iatm = atomlist(ix,grouplist(1,i))
          icel = id_g2l(1,iatm)

          if (icel > 0 .and. icel <= ncell) then
            nrestraint(icel) = nrestraint(icel) + 1
            k = k + 1
            restraint_atom(nrestraint(icel),icel) = iatm
            restraint_coord(1:3,nrestraint(icel),icel) = coord(1:3,iatm)
          end if

        end do
#ifdef HAVE_MPI_GENESIS
        call mpi_allreduce(k, enefunc%nrmsd, 1, mpi_integer, mpi_sum, &
                           mpi_comm_country, ierror)
#else
        enefunc%nrmsd = k
#endif
        enefunc%rmsd_force  = const(1,i)
        enefunc%rmsd_target = ref(1,i)
        enefunc%rmsd_withmass =  &
          (enefunc%restraint_kind(i) == REstraintsFuncRMSDCOM)

      end if

    end do

    if (position_flag .and. rmsd_flag) &
     call error_msg('Setup_Enefunc_Rest_Domain> RMSD and POSI'//&
                    'are not allowed at once')

    j = 0
    do i = 1, num_funcs

      if (enefunc%restraint_kind(i) == RestraintsFuncPC) then

        j = j + 1
        k = 0
        if (j == 1) then
          do ix = 1, num_atoms(grouplist(1,i))
            iatm = atomlist(ix,grouplist(1,i))
            icel = id_g2l(1,iatm)
            restraint_g2pc(iatm) = ix

            if (icel > 0 .and. icel <= ncell) then
              nrestraint(icel) = nrestraint(icel) + 1
              k = k + 1
              restraint_atom(nrestraint(icel),icel) = iatm
              restraint_coord(1:3,nrestraint(icel),icel) = coord(1:3,iatm)
            end if
          end do
        end if
        enefunc%num_atoms_pc_restraint = num_atoms(grouplist(1,i))
        enefunc%pc_force(j)  = const(1,i)
        enefunc%pc_target(j) = ref(1,i)

      end if

    end do
    enefunc%num_pc_modes = j

    ! distance restraint
    !
    k = 0
    do i = 1, num_funcs

      if (enefunc%restraint_kind(i) == RestraintsFuncDIST     .or. &
          enefunc%restraint_kind(i) == RestraintsFuncDISTCOM  .or. &
          enefunc%restraint_kind(i) == RestraintsFuncANGLE    .or. &
          enefunc%restraint_kind(i) == RestraintsFuncANGLECOM .or. &
          enefunc%restraint_kind(i) == RestraintsFuncDIHED    .or. &
          enefunc%restraint_kind(i) == RestraintsFuncDIHEDCOM .or. &
          enefunc%restraint_kind(i) == RestraintsFuncREPUL    .or. &
          enefunc%restraint_kind(i) == RestraintsFuncREPULCOM .or. &
          enefunc%restraint_kind(i) == RestraintsFuncFB       .or. &
          enefunc%restraint_kind(i) == RestraintsFuncFBCOM) then
 
        do ig = 1, enefunc%restraint_funcgrp(i)
          do ix = 1, num_atoms(grouplist(ig,i))
            k = k + 1
            enefunc%restraint_bondslist_to_atomlist(k) =  &
                 atomlist(ix,grouplist(ig,i))
            bondslist(ix,grouplist(ig,i)) = k
          end do
        end do

      end if

    end do

    if (position_flag .or. rmsd_flag) then
      call alloc_enefunc(enefunc, EneFuncRestDomain, ncell)
    end if

    return

  end subroutine setup_enefunc_rest_domain

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_restraints_grouplist
  !> @brief        setup restraint grouplist
  !! @authors      CK, TM
  !! @param[in]    ifunc      : number for function
  !! @param[in]    restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_restraints_grouplist(ifunc, restraints, enefunc)

    ! formal arguments
    integer,                 intent(in)    :: ifunc
    type(s_restraints),      intent(in)    :: restraints
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variable
    integer                       :: j, max_grp
    character(MaxLine)            :: string
    integer,          allocatable :: idata(:)


    max_grp = enefunc%max_restraint_numgrps
    allocate(idata(max_grp))

    string = restraints%select_index(ifunc)

    call split(split_num(string), max_grp, string, idata)
    do j = 1, max_grp 
      enefunc%restraint_grouplist(j,ifunc) = idata(j)
    end do

    deallocate(idata)

    return

  end subroutine setup_restraints_grouplist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_restraints_exponent_func
  !> @brief        setup restraint exponent func
  !! @authors      CK, TM
  !! @param[in]    ifunc      : number for function
  !! @param[in]    restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_restraints_exponent_func(ifunc, restraints, enefunc)

    ! formal arguments
    integer,                 intent(in)    :: ifunc
    type(s_restraints),      intent(in)    :: restraints
    type(s_enefunc),         intent(inout) :: enefunc


    if (restraints%exponent(ifunc) > 0) then
      enefunc%restraint_exponent_func(ifunc) = restraints%exponent(ifunc)
    else
      enefunc%restraint_exponent_func(ifunc) = 2
    end if

    return

  end subroutine setup_restraints_exponent_func

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_restraints_constants
  !> @brief        setup restraint constants
  !! @authors      CK, TM
  !! @param[in]    ifunc      : number for function
  !! @param[in]    restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_restraints_constants(ifunc, restraints, enefunc)

    ! formal arguments
    integer,                 intent(in)    :: ifunc
    type(s_restraints),      intent(in)    :: restraints
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variable
    integer                  :: i
    integer                  :: ndata
    integer                  :: ndata_max
    real(wp),    allocatable :: ddata(:)
    character(MaxLine)       :: string
    integer                  :: num_funcs


    num_funcs = enefunc%num_restraintfuncs
    string = restraints%constant(ifunc)
    ndata = split_num(string)

    if (ndata == 0) then
      call error_msg('Setup_Restraint_Constants> constant is not given')
    end if

    if (allocated(enefunc%restraint_const_replica)) then
      ndata_max = size(enefunc%restraint_const_replica(:,1))
      if (ndata > ndata_max) then
        call error_msg('Setup_Restraint_Constants> the number of values for constant shoule be same in all functions')
      end if
    else
      ndata_max = ndata;
      call alloc_enefunc(enefunc, EneFuncRefr, num_funcs, ndata)
    end if

    allocate(ddata(ndata_max))

    call split(ndata, ndata_max, string, ddata)

    do i = 1, ndata
      enefunc%restraint_const_replica(i,ifunc) = ddata(i)
    end do

    deallocate(ddata)

    if (enefunc%restraint_kind(ifunc) == RestraintsFuncPOSI) then

      ! if (ndata == 1) then
      !   call split(ndata, 4, string, ddata)
      !   enefunc%restraint_const(1,ifunc) = ddata(1)
      !   enefunc%restraint_const(2,ifunc) = 1.0_wp
      !   enefunc%restraint_const(3,ifunc) = 1.0_wp
      !   enefunc%restraint_const(4,ifunc) = 1.0_wp
      ! else if (ndata == 4) then
      !   call split(ndata, 4, string, ddata)
      !   enefunc%restraint_const(1,ifunc) = ddata(1)
      !   enefunc%restraint_const(2,ifunc) = ddata(2)
      !   enefunc%restraint_const(3,ifunc) = ddata(3)
      !   enefunc%restraint_const(4,ifunc) = ddata(4)
      ! else
      !   call error_msg('Setup_Restraints_Constants> Error in const')
      ! end if

      enefunc%restraint_const(1,ifunc) = enefunc%restraint_const_replica(1,ifunc)
      if (restraints%direction(ifunc) == RestraintsDirALL) then
        enefunc%restraint_const(2:4,ifunc) = 1.0_wp
      else
        enefunc%restraint_const(2:4,ifunc) = 0.0_wp
        enefunc%restraint_const(restraints%direction(ifunc),ifunc) = 1.0_wp
      end if

    else

      ! if (ndata == 1) then
      !   call split(ndata, 4, string, ddata)
      !   enefunc%restraint_const(1,ifunc) = ddata(1)
      !   enefunc%restraint_const(2,ifunc) = ddata(1)
      !   enefunc%restraint_const(3,ifunc) = ddata(3) ! dummy
      !   enefunc%restraint_const(4,ifunc) = ddata(4) ! dummy
      ! else if (ndata == 2) then
      !   call split(ndata, 4, string, ddata)
      !   enefunc%restraint_const(1,ifunc) = ddata(1)
      !   enefunc%restraint_const(2,ifunc) = ddata(2)
      !   enefunc%restraint_const(3,ifunc) = ddata(3) ! dummy
      !   enefunc%restraint_const(4,ifunc) = ddata(4) ! dummy
      ! else
      !   call error_msg('Setup_Restraints_Constants> Error in const')
      ! end if

      enefunc%restraint_const(1,ifunc) = enefunc%restraint_const_replica(1,ifunc)
      enefunc%restraint_const(2,ifunc) = enefunc%restraint_const_replica(1,ifunc)
      enefunc%restraint_const(3,ifunc) = 1.0_wp
      enefunc%restraint_const(4,ifunc) = 1.0_wp

    end if

    return

  end subroutine setup_restraints_constants

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_restraints_reference
  !> @brief        setup restraint reference
  !! @authors      CK, TM
  !! @param[in]    ifunc      : number for function
  !! @param[in]    restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_restraints_reference(ifunc, restraints, enefunc)

    ! formal arguments
    integer,                 intent(in)    :: ifunc
    type(s_restraints),      intent(in)    :: restraints
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variable
    integer                  :: i
    integer                  :: ndata
    integer                  :: ndata_max
    real(wp),    allocatable :: ddata(:)
    character(MaxLine)       :: string
    integer                  :: num_funcs


    num_funcs = enefunc%num_restraintfuncs
    string = restraints%reference(ifunc)
    ndata = split_num(string)

    if (ndata > 0) then
      if (enefunc%restraint_kind(ifunc) ==  RestraintsFuncPOSI) then
      call error_msg('Setup_Restraint_Reference> reference of POSI is not needed')
      end if
      if (enefunc%restraint_kind(ifunc) ==  RestraintsFuncEM) then
        call error_msg('Setup_Restraint_Reference> reference of EM is not needed')
      end if

    end if

    if (ndata == 0) then
      if  (enefunc%restraint_kind(ifunc) /= RestraintsFuncPOSI    .and. &
           enefunc%restraint_kind(ifunc) /= RestraintsFuncEM) then
        call error_msg('Setup_Restraint_Reference> reference is not given')
      end if
    end if

    if (allocated(enefunc%restraint_ref_replica)) then
      ndata_max = size(enefunc%restraint_const_replica(:,1))
      if (ndata > ndata_max) then
        call error_msg('Setup_Restraint_Reference> the number of values for constant shoule be same in all functions')
      end if
    else
      ndata_max = ndata;
      call alloc_enefunc(enefunc, EneFuncRefr, num_funcs, ndata)
    end if

    allocate(ddata(ndata_max))

    call split(ndata, ndata_max, string, ddata)

    do i = 1, ndata
      enefunc%restraint_ref_replica(i,ifunc) = ddata(i)
    end do

    deallocate(ddata)

    if (ndata > 0) then
      enefunc%restraint_ref(1,ifunc) = enefunc%restraint_ref_replica(1,ifunc)
      enefunc%restraint_ref(2,ifunc) = enefunc%restraint_ref_replica(1,ifunc)
    end if

    return

  end subroutine setup_restraints_reference

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_restraints_exponent_dist
  !> @brief        setup restraint exponent dist
  !! @authors      CK, TM
  !! @param[in]    ifunc      : number for function
  !! @param[in]    restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_restraints_exponent_dist(ifunc, restraints, enefunc)

    ! formal arguments
    integer,                 intent(in)    :: ifunc
    type(s_restraints),      intent(in)    :: restraints
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variable
    integer                       :: ndata, factor, maxfunc, max_grp
    character(MaxLine)            :: string
    integer,          allocatable :: idata(:)


    select case (enefunc%restraint_kind(ifunc))
    case(RestraintsFuncDIST:RestraintsFuncDISTCOM)
      factor = 2
    case default
      return
    end select

    max_grp = enefunc%max_restraint_numgrps
    maxfunc = max(int(max_grp/factor), 1)

    ! initialize
    !
    enefunc%restraint_exponent_dist(1:maxfunc,ifunc) = 1

    ! set exponent_dist
    !
    string = restraints%exponent_dist(ifunc)
    ndata = split_num(string)

    allocate(idata(1:maxfunc))
    if (ndata > 0) then
      if (ndata /= enefunc%restraint_funcgrp(ifunc)/factor) &
        call error_msg('Setup_Restraints_Exponent_Dist> Error in exponent_dist')

      call split(ndata, maxfunc, string, idata)
      enefunc%restraint_exponent_dist(1:maxfunc,ifunc) = idata(1:maxfunc)
    end if
    deallocate(idata)

    return

  end subroutine setup_restraints_exponent_dist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_restraints_weight_dist
  !> @brief        setup restraint weight dist
  !! @authors      CK, TM
  !! @param[in]    ifunc      : number for function
  !! @param[in]    restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_restraints_weight_dist(ifunc, restraints, enefunc)

    ! formal arguments
    integer,                 intent(in)    :: ifunc
    type(s_restraints),      intent(in)    :: restraints
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variable
    integer                  :: ndata, factor, maxfunc, max_grp
    character(MaxLine)       :: string
    real(wp),    allocatable :: ddata(:)


    select case (enefunc%restraint_kind(ifunc))
    case(RestraintsFuncDIST:RestraintsFuncDISTCOM)
      factor = 2
    case default
      return
    end select

    max_grp = enefunc%max_restraint_numgrps
    maxfunc = max(int(max_grp/factor), 1)

    ! initialize
    !
    enefunc%restraint_weight_dist(1:maxfunc,ifunc) = 1.0_wp

    ! set weight_dist
    !
    string = restraints%weight_dist(ifunc)
    ndata = split_num(string)

    allocate(ddata(1:maxfunc))
    if (ndata > 0) then
      if (ndata /= enefunc%restraint_funcgrp(ifunc)/factor) &
        call error_msg('Setup_Restraints_Weight_Dist> Error in weight_dist')

      call split(ndata, maxfunc, string, ddata)
      enefunc%restraint_weight_dist(1:maxfunc,ifunc) = ddata(1:maxfunc)
    end if
    deallocate(ddata)

    return

  end subroutine setup_restraints_weight_dist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_rest_refcoord
  !> @brief        setup enefunc restraint reference coords
  !! @authors      CK, TM
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_rest_refcoord(molecule, enefunc)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: j


    if (size(molecule%atom_refcoord(1,:)) /= molecule%num_atoms) &
      call error_msg('Setup_Enefunc_Rest_Refcoord> bad reffile in [INPUT]')

    enefunc%num_atoms_ref = molecule%num_atoms

    call alloc_enefunc(enefunc, EneFuncRefc, molecule%num_atoms)

    do j = 1, molecule%num_atoms
      enefunc%restraint_refcoord(1,j) = molecule%atom_refcoord(1,j)
      enefunc%restraint_refcoord(2,j) = molecule%atom_refcoord(2,j)
      enefunc%restraint_refcoord(3,j) = molecule%atom_refcoord(3,j)
    end do

    return

  end subroutine setup_enefunc_rest_refcoord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_restraint_mode
  !> @brief        setup restraint mode (for PC)
  !! @authors      YK, JJ
  !! @param[in]    ifunc      : number for function
  !! @param[in]    restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_restraint_mode(ifunc, restraints, enefunc)

    ! formal arguments
    integer,                 intent(in)    :: ifunc
    type(s_restraints),      intent(in)    :: restraints
    type(s_enefunc),         intent(inout) :: enefunc


    if (restraints%mode(ifunc) > 0) then
      enefunc%restraint_mode(ifunc) = restraints%mode(ifunc)
    else
      enefunc%restraint_mode(ifunc) = 1
    end if

    return

  end subroutine setup_restraint_mode

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_rest_mode
  !> @brief        setup enefunc principal component mode
  !! @authors      YK, JJ
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_rest_mode(molecule, enefunc)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc


    if (size(molecule%pc_mode(:)) /= molecule%num_pc_modes) &
      call error_msg('Setup_Enefunc_Rest_Modes> bad modefile in [INPUT]')

    call alloc_enefunc(enefunc, EneFuncMode, molecule%num_pc_modes, &
                       molecule%num_atoms)

    enefunc%pc_mode(:) = molecule%pc_mode(:)

    return

  end subroutine setup_enefunc_rest_mode

end module sp_enefunc_restraints_mod
