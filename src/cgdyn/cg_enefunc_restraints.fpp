!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_enefunc_restraints_mod
!> @brief   Provides restraint energy functions.
!! @authors Chigusa Kobayashi (CK), Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8
  
#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_enefunc_restraints_mod

  use cg_restraints_str_mod
  use cg_enefunc_str_mod
  use cg_domain_str_mod
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
  public   :: setup_enefunc_restraints
  public   :: setup_restraint_mode
  private  :: setup_enefunc_rest_group
  private  :: setup_enefunc_rest_func
  public   :: setup_restraints_grouplist
  public   :: setup_restraints_exponent_func
  public   :: setup_restraints_constants
  public   :: setup_restraints_reference
  public   :: setup_restraints_exponent_dist
  public   :: setup_restraints_weight_dist
  private  :: setup_enefunc_rest_refcoord
  private  :: setup_enefunc_rest_mode

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_restraints
  !> @brief        Defines restraint potential.
  !! @authors      CK, TM
  !! @param[in]    molecule   : molecule information
  !! @param[in]    restraints : restraints information
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
    integer                  :: i, j
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

    ! summary of setup enefunc_restraint
    !
    if (main_rank) then

      write(MsgOut,'(A)') 'Setup_Enefunc_Restraints> Set up restraint groups.'

      do i = 1, enefunc%num_restraintgroups
        write(MsgOut,'(A,I5,3A)') ' Group = ', i, ', "', trim(restraints%group(i)),'".'
        write(MsgOut,'(A,I5)')   ' Number of atoms = ', enefunc%restraint_numatoms(i)
        write(MsgOut,'(A)')      ' Atom list: ' 
        if (restraints%verbose) then
          oldres   =-10000
          oldsegnm = "    "
          oldresnm = "      "
          do j = 1, enefunc%restraint_numatoms(i) 
            iatm=enefunc%restraint_atomlist(j,i)
            if (iatm <= 0 .or. iatm > molecule%num_atoms) &
              call error_msg('Error: Setup_Enefunc_Restraint - Atom index '//&
	                     'is not correct.')
            segname = molecule%segment_name(iatm)
            resname = molecule%residue_name(iatm)
            atname  = molecule%atom_name(iatm)
            ires    = molecule%residue_no(iatm)
            if (ires /= oldres .or. resname /= oldresnm .or. &
               segname /= oldsegnm) then
              if (oldres /= -10000) then
                write(MsgOut,'(A)') ''
              endif
              write(MsgOut,'(5X,A4,A,A6,A,i6,A,$)') & 
                  trim(segname),':',trim(resname),':',ires,':'
              write(MsgOut,'(1X,A,$)') trim(atname) 
              oldres = ires
              oldresnm = resname
              oldsegnm = segname
            else
              write(MsgOut,'(A,1X,A,$)') ',',trim(atname) 
            endif
          end do
        else
          do j = 1, enefunc%restraint_numatoms(i) 
            write(MsgOut,'(i7,$)') enefunc%restraint_atomlist(j,i)
            if (mod(j,10) == 0 .and. j /= enefunc%restraint_numatoms(i)) then
              write(MsgOut,'(A)') ''
            end if
          end do
        endif
        write(MsgOut,'(A)') ''
      end do

      write(MsgOut,'(A)') 'Setup_Enefunc_Restraints> Set up restraint functions.'

      do i = 1, enefunc%num_restraintfuncs

        write(MsgOut,'(A,I5,A,I5)') &
          ' Function = ', i, ' Kind = ', enefunc%restraint_kind(i)

        ! summary for positional restraint
        !
        if (enefunc%restraint_kind(i) == RestraintsFuncPOSI .or. &
            enefunc%restraint_kind(i) == RestraintsFuncEM) then
          write(MsgOut,'(A,4F8.3)') &
            ' Constant (total, x, y, z) = ', enefunc%restraint_const(1:4,i) 
        else if (enefunc%restraint_kind(i) == RestraintsFuncRMSD .or. &
                 enefunc%restraint_kind(i) == RestraintsFuncRMSDCOM) then
          write(MsgOut,'(A,F8.3)') &
            ' Force constant        = ', enefunc%restraint_const(1,i)
          write(MsgOut,'(A,F8.3)') &
            ' Target RMSD           = ', enefunc%restraint_ref(1,i)
        else
          write(MsgOut,'(A,F8.3)') &
            ' Constant                 = ',   enefunc%restraint_const(1,i)
          write(MsgOut,'(A,F8.3)') &
            ' Reference                   = ',   enefunc%restraint_ref(1,i)
        end if
 
        write(MsgOut,'(A,I5)') ' Exponent of function = ',                &
                               enefunc%restraint_exponent_func(i)

        ! summary for distance restraint
        !
        if (enefunc%restraint_kind(i) == RestraintsFuncDIST .or.          &
            enefunc%restraint_kind(i) == RestraintsFuncDISTCOM) then

          write(MsgOut,'(A,I5)') ' Number of distances  = ', enefunc%restraint_funcgrp(i)/2
          write(MsgOut,'(A)')    ' Distance list: '
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

          write(MsgOut,'(A,I5)') ' Number of angles     = ', enefunc%restraint_funcgrp(i)/3
          write(MsgOut,'(A)')    ' Angle list: '
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

          write(MsgOut,'(A,I5)') ' Number of dihedrals  = ', enefunc%restraint_funcgrp(i)/4
          write(MsgOut,'(A)')    ' Dihedral list: ' 
          do j = 1, enefunc%restraint_funcgrp(i)/4
            write(MsgOut,'(" group = (",i3,",",i3,",",i3,",",i3,")"$)')   &
                           enefunc%restraint_grouplist(4*j-3,i),          &
                           enefunc%restraint_grouplist(4*j-2,i),          &
                           enefunc%restraint_grouplist(4*j-1,i),          &
                           enefunc%restraint_grouplist(4*j,i)
            write(MsgOut,'(A)') ''
          end do

        ! summary for positinal and RMSD restraint
        !
        else 

          write(MsgOut,'(A,I5)') ' Number of groups  = ', enefunc%restraint_funcgrp(i)
          write(MsgOut,'(A)') "Group list:" 
          do j = 1, enefunc%restraint_funcgrp(i) 
            write(MsgOut,'(i3,$)') enefunc%restraint_grouplist(j,i)
            if (mod(j,20) == 0 .and. j /= enefunc%restraint_funcgrp(i) )  &
              write(MsgOut,'(A)') ''
          end do
          write(MsgOut,'(A)') ''

        end if
        write(MsgOut,'(A)') ''

      end do
    endif

    ! error message when both position and rmsd restraints are assigned
    !
    if (enefunc%restraint_posi .and. enefunc%restraint_rmsd) &
      call error_msg('Error: Setup_Enefunc_Restraint - Choose either '//&
                     'position or RMSD restraint.')

    return

  end subroutine setup_enefunc_restraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_rest_group
  !> @brief        Defines restraint atom groups.
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
    real(dp)                 :: totmass, coef
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
        call error_msg('Error: Setup_Enefunc_Rest_Group - No atom in group.')
    end do

    ! setup mass_coefficient
    !
    enefunc%restraint_masscoef(1:max_atoms,1:num_groups) = 0.0_dp
    do i = 1, enefunc%num_restraintgroups
      totmass = 0.0_dp
      do j = 1, enefunc%restraint_numatoms(i) 
        jatom = enefunc%restraint_atomlist(j,i)
        totmass = totmass + molecule%mass(jatom)
      end do
      coef = 1.0_dp/totmass

      do j = 1, enefunc%restraint_numatoms(i) 
        jatom = enefunc%restraint_atomlist(j,i)
        enefunc%restraint_masscoef(j,i) = molecule%mass(jatom) * real(coef,wp)
      end do
    end do

    return

  end subroutine setup_enefunc_rest_group

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_rest_func
  !> @brief        Defines restraint functions.
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
          call error_msg('Error: Setup_Enefunc_Rest_Func - Invalid index '//&
	                 'for DIST.')

      case(RestraintsFuncANGLE:RestraintsFuncANGLECOM)
        if (ndata /= 3) &
          call error_msg('Error: Setup_Enefunc_Rest_Func - Invalid index '//&
	                 'for ANGLE.')

      case(RestraintsFuncDIHED:RestraintsFuncDIHEDCOM)
        if (ndata /= 4) &
          call error_msg('Error: Setup_Enefunc_Rest_Func - Invalid index '//&
	                 'for DIHE.')

      case(RestraintsFuncPC:RestraintsFuncPCCOM)
        if (ndata == 0) &
          call error_msg('Error: Setup_Enefunc_Rest_Func - Invalid index '//&
	                 'for PC.')

      case default
        if (ndata /= 1) & 
          call error_msg('Error: Setup_Enefunc_Rest_Func - Invalid index '//&
	                 'for POSI, RMSD, or RG.')

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
  !> @brief        Sets up restraint information about the domain.
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
    integer                  :: i, ig, j, k, ix, iatm, icel, step
    integer                  :: num
    logical                  :: position_flag
    logical                  :: rmsd_flag

    real(wp),        pointer :: const(:,:), ref(:,:), coord(:,:)
    real(wp),        pointer :: restraint_force(:,:)
    real(wp),        pointer :: restraint_coord(:,:)
    integer,         pointer :: id_g2l(:), ncell
    integer,         pointer :: num_funcs, num_atoms(:)
    integer,         pointer :: grouplist(:,:), atomlist(:,:), bondslist(:,:)
    integer,         pointer :: restraint_g2pc(:)
    integer,         pointer :: restraint_atom(:), atom_2_cell(:)


    ncell           => domain%num_cell_local
    id_g2l          => domain%id_g2l
    atom_2_cell     => domain%atom_2_cell

    num_funcs       => enefunc%num_restraintfuncs
    num_atoms       => enefunc%restraint_numatoms
    grouplist       => enefunc%restraint_grouplist
    atomlist        => enefunc%restraint_atomlist
    bondslist       => enefunc%restraint_bondslist
    const           => enefunc%restraint_const
    ref             => enefunc%restraint_ref
    coord           => enefunc%restraint_refcoord
    restraint_g2pc  => enefunc%restraint_g2pc

    position_flag=.false.
    rmsd_flag=.false.

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
          enefunc%restraint_kind(i) == RestraintsFuncFBCOM    .or. &
          enefunc%restraint_kind(i) == RestraintsFuncCONF     .or. &
          enefunc%restraint_kind(i) == RestraintsFuncCONFCOM) then
        do ig = 1, enefunc%restraint_funcgrp(i)
          do ix = 1, num_atoms(grouplist(ig,i))
            k = k + 1
          end do
        end do
      end if
    end do
    enefunc%num_atoms_bonds_restraint = k

    ! position/rmsd/pc restraint
    !
    do step = 1, 2

      num = 0
      j = 0

      do i = 1, num_funcs

        if (enefunc%restraint_kind(i) == RestraintsFuncPOSI) then

          position_flag=.true.
          do ix = 1, num_atoms(grouplist(1,i))
            iatm = atomlist(ix,grouplist(1,i))
            icel = id_g2l(iatm)
            if (icel > 0) then
              icel = atom_2_cell(icel)
              if (icel > 0 .and. icel <= ncell) then
                num = num + 1
                if (step == 2) then
                  restraint_atom(num) = iatm
                  restraint_force(1:4,num) = const(1:4,i)
                  restraint_coord(1:3,num) = coord(1:3,iatm)
                end if
              end if
            end if
          end do

        else if (enefunc%restraint_kind(i) == REstraintsFuncRMSD .or.     &
                 enefunc%restraint_kind(i) == REstraintsFuncRMSDCOM) then

          rmsd_flag=.true.
          do ix = 1, num_atoms(grouplist(1,i))
            iatm = atomlist(ix,grouplist(1,i))
            icel = id_g2l(iatm)
            if (icel > 0) then
              icel = atom_2_cell(icel)
              if (icel > 0 .and. icel <= ncell) then
                num = num + 1
                if (step == 2) then
                  restraint_atom(num) = iatm
                  restraint_coord(1:3,num) = coord(1:3,iatm)
                end if
              end if
            end if
          end do
          if (step == 1) then
            call mpi_allreduce(num, enefunc%nrmsd, 1, mpi_integer, &
                               mpi_sum, mpi_comm_country, ierror)
            enefunc%rmsd_force  = const(1,i)
            enefunc%rmsd_target = ref(1,i)
            enefunc%rmsd_withmass =  &
              (enefunc%restraint_kind(i) == REstraintsFuncRMSDCOM )
          end if

        else if (enefunc%restraint_kind(i) == RestraintsFuncPC) then

          j = j + 1
          if (j == 1) then
            do ix = 1, num_atoms(grouplist(1,i))
              iatm = atomlist(ix,grouplist(1,i))
              icel = id_g2l(iatm)
              if (icel > 0) then
                icel = atom_2_cell(icel)
                restraint_g2pc(iatm) = ix
                if (icel > 0 .and. icel <= ncell) then
                  num = num + 1
                  if (step == 2) then
                    restraint_atom(num) = iatm
                    restraint_coord(1:3,num) = coord(1:3,iatm)
                  end if
                end if
              end if
            end do
          end if
          if (step == 1) then
            enefunc%num_atoms_pc_restraint = num_atoms(grouplist(1,i))
            enefunc%pc_force(j)  = const(1,i)
            enefunc%pc_target(j) = ref(1,i)
          end if

        end if

      end do

      if (step == 1) then
        if (position_flag .and. rmsd_flag) &
          call error_msg('Setup_Enefunc_Rest_Domain> RMSD and POSI'//&
                         'are not allowed at once')
        if (position_flag .and. j > 0) &
          call error_msg('Setup_Enefunc_Rest_Domain> POSI and PC'//&
                         'are not allowed at once')
        if (rmsd_flag .and. j > 0) &
          call error_msg('Setup_Enefunc_Rest_Domain> RMSD and PC'//&
                         'are not allowed at once')
        enefunc%num_pc_modes = j
        enefunc%num_rest_domain = num
        k = num
        call mpi_allreduce(k, MaxRest, 1, mpi_integer, mpi_max, &
                           mpi_comm_country, ierror)
        MaxRest = MaxRest * 2
        k = enefunc%num_atoms_bonds_restraint
        call alloc_enefunc(enefunc, EneFuncRest, k)
        if (position_flag .or. rmsd_flag) &
          call alloc_enefunc(enefunc, EneFuncRestDomain, MaxRest)

        restraint_atom  => enefunc%restraint_atom
        restraint_force => enefunc%restraint_force
        restraint_coord => enefunc%restraint_coord
      end if

    end do

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
          enefunc%restraint_kind(i) == RestraintsFuncFBCOM    .or. &
          enefunc%restraint_kind(i) == RestraintsFuncCONF     .or. &
          enefunc%restraint_kind(i) == RestraintsFuncCONFCOM) then

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

    return

  end subroutine setup_enefunc_rest_domain

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_restraints_grouplist
  !> @brief        Sets up restraint grouplists.
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
  !> @brief        Sets up restraint exponent functions.
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
  !> @brief        Sets up restraint constants.
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
      call error_msg('Error: Setup_Restraint_Constants - Constant is not '//&
                     'given.')
    end if

    if (allocated(enefunc%restraint_const_replica)) then
      ndata_max = size(enefunc%restraint_const_replica(:,1))
      if (ndata > ndata_max) then
        call error_msg('Error: Setup_Restraint_Constants - The number of '//&
		       'values for constant should be the same in all    '//&
                       'functions.')
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
      !   call error_msg('Error: Setup_Restraints_Constants - Invalid constant.')
      ! end if

      enefunc%restraint_const(1,ifunc) = &
              enefunc%restraint_const_replica(1,ifunc)
      if (restraints%direction(ifunc) == RestraintsDirALL) then
        enefunc%restraint_const(2:4,ifunc) = 1.0_wp
      else
        enefunc%restraint_const(2:4,ifunc) = 0.0_wp
        enefunc%restraint_const(restraints%direction(ifunc),ifunc) = 1.0_wp
      endif

    else if (enefunc%restraint_kind(ifunc) == RestraintsFuncEM) then

      enefunc%restraint_const(1,  ifunc) = &
              enefunc%restraint_const_replica(1,ifunc)
      enefunc%restraint_const(2:4,ifunc) = 1.0_wp

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
      !   call error_msg('Error: Setup_Restraints_Constants - Invalid constant.')
      ! end if

      enefunc%restraint_const(1,ifunc) = &
              enefunc%restraint_const_replica(1,ifunc)
      enefunc%restraint_const(2,ifunc) = &
              enefunc%restraint_const_replica(1,ifunc)
      enefunc%restraint_const(3,ifunc) = 1.0_wp
      enefunc%restraint_const(4,ifunc) = 1.0_wp

    end if

    return

  end subroutine setup_restraints_constants

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_restraints_reference
  !> @brief        Sets up restraint references.
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
        call error_msg('Error: Setup_Restraint_Reference - Reference for '//&
		       'POSI is not needed.')
      end if
      if (enefunc%restraint_kind(ifunc) ==  RestraintsFuncEM) then
        call error_msg('Error: Setup_Restraint_Reference - Reference for '//&
		       'EM is not needed.')
      end if
    end if

    if (ndata == 0) then
      if (enefunc%restraint_kind(ifunc) /=  RestraintsFuncPOSI .and. &
         (enefunc%restraint_kind(ifunc) /=  RestraintsFuncEM  )) then
        call error_msg('Error: Setup_Restraint_Reference - Reference is '//&
		       'not given.')
      end if
    end if

    if (allocated(enefunc%restraint_ref_replica)) then
      ndata_max = size(enefunc%restraint_const_replica(:,1))
      if (ndata > ndata_max) then
        call error_msg('Setup_Restraint_Reference> the number of values '//&
		       'for constant shoule be same in all functions')
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
  !> @brief        Sets up restraint exponent distances.
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
              call error_msg('Error: Setup_Restraints_Exponent_Dist - '//&
	             'Invalid exponent_dist.')
      call split(ndata, maxfunc, string, idata)
      enefunc%restraint_exponent_dist(1:maxfunc,ifunc) = idata(1:maxfunc)
    end if
    deallocate(idata)

    return

  end subroutine setup_restraints_exponent_dist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_restraints_weight_dist
  !> @brief        Sets up restraint weight distances.
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
              call error_msg('Error: Setup_Restraints_Weight_Dist - '//&
	             'Invalid weight_dist.')
      call split(ndata, maxfunc, string, ddata)
      enefunc%restraint_weight_dist(1:maxfunc,ifunc) = ddata(1:maxfunc)
    end if
    deallocate(ddata)

    return

  end subroutine setup_restraints_weight_dist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_rest_refcoord
  !> @brief        Sets up energy function restraint reference coordinates.
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
      call error_msg('Error: Setup_Enefunc_Rest_Refcoord - '//&
             'Invalid reffile in [INPUT].')

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
  !> @brief        Sets up restraint mode (for PC).
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
  !> @brief        Sets up energy function principal component mode.
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
      call error_msg('Error: Setup_Enefunc_Rest_Modes - '//&
             'Invalid modefile in [INPUT].')

    call alloc_enefunc(enefunc, EneFuncMode, molecule%num_pc_modes, &
                       molecule%num_atoms)

    enefunc%pc_mode(:) = molecule%pc_mode(:)

    return

  end subroutine setup_enefunc_rest_mode

end module cg_enefunc_restraints_mod
