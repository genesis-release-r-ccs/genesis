!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_enefunc_restraint_mod
!> @brief   restraint energy functions
!! @authors Chigusa Kobayashi (CK), Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_enefunc_restraints_mod

  use at_restraints_str_mod
  use at_enefunc_str_mod
  use molecules_str_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public   :: setup_enefunc_restraints
  private  :: setup_enefunc_rest_group
  private  :: setup_enefunc_rest_func
  private  :: setup_restraint_grouplist
  private  :: setup_restraint_exponent_func
  public   :: setup_restraint_constants      ! this is also called from remd
  public   :: setup_restraint_reference      ! this is also called from remd
  private  :: setup_restraint_exponent_dist
  private  :: setup_restraint_weight_dist
  private  :: setup_restraint_mode
  private  :: setup_enefunc_rest_refcoord
  private  :: setup_enefunc_rest_mode
  private  :: setup_restraint_caging
  private  :: setup_restraint_flat_radius
  private  :: setup_enefunc_rest_rg

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_restraints
  !> @brief        define restraint potential
  !! @authors      CK, TM
  !! @param[in]    molecule   : molecule information
  !! @param[in]    restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_restraints(molecule, restraints, enefunc)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_restraints),      intent(in)    :: restraints
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: istart, iend
    integer                  :: i, j, k
    integer                  :: iatm, ires, oldres
    character(4)             :: segname, oldsegnm
    character(6)             :: resname, oldresnm
    character(4)             :: atname


    enefunc%restraint_flag     = restraints%restraint_flag
    enefunc%pressure_position  = restraints%pressure_position
    enefunc%pressure_rmsd      = restraints%pressure_rmsd  

    if (.not. restraints%restraint_flag) return

    ! setup group
    !
    call setup_enefunc_rest_group(molecule, restraints, enefunc)

    ! setup func
    !
    call setup_enefunc_rest_func(molecule, restraints, enefunc)

    ! define restraint functions for each processor
    !
    ! call get_loop_index(enefunc%num_restraintfuncs, istart, iend)
    !
    ! enefunc%istart_restraint = istart
    ! enefunc%iend_restraint   = iend
    if (main_rank .or. replica_main_rank) then
      enefunc%istart_restraint = 1
      enefunc%iend_restraint   = enefunc%num_restraintfuncs
    else
      enefunc%istart_restraint = 0
      enefunc%iend_restraint   = -1
    end if

    ! summary of setup enefunc_restraint
    !
    if (main_rank) then

      write(MsgOut,'(A)') 'Setup_Enefunc_Restraint> Setup restraint groups'

      do i = 1, enefunc%num_restraintgroups
        write(MsgOut,'(A,I5,3A)') ' group = ', i, ', "', trim(restraints%group(i)),'"'
        write(MsgOut,'(A,I5)')   ' # of atoms = ', enefunc%restraint_numatoms(i)
        write(MsgOut,'(A)')      ' atomlist: ' 
        if (restraints%verbose) then
          oldres   =-10000
          oldsegnm = "    "
          oldresnm = "      "
          do j = 1, enefunc%restraint_numatoms(i) 
            iatm=enefunc%restraint_atomlist(j,i)
            if (iatm <= 0 .or. iatm > molecule%num_atoms) &
              call error_msg('Setup_Enefunc_Restraint> atom idx is not correct.')
            segname = molecule%segment_name(iatm)
            resname = molecule%residue_name(iatm)
            atname  = molecule%atom_name(iatm)
            ires    = molecule%residue_no(iatm)
            if (ires .ne. oldres .or. resname .ne. oldresnm .or. &
               segname .ne. oldsegnm) then
              if (oldres .ne. -10000) then
                write(MsgOut,'(A)') ''
              end if
              write(MsgOut,'(5X,A4,A,A6,A,i6,A,$)') & 
                  trim(segname),'.',trim(resname),'.',ires,':'
              write(MsgOut,'(1X,A,$)') trim(atname) 
              oldres = ires
              oldresnm = resname
              oldsegnm = segname
            else
              write(MsgOut,'(A,1X,A,$)') ',',trim(atname) 
            end if
          end do
        else
          do j = 1, enefunc%restraint_numatoms(i) 
            write(MsgOut,'(i7,$)') enefunc%restraint_atomlist(j,i)
            if (mod(j,10) == 0 .and. j /= enefunc%restraint_numatoms(i)) then
              write(MsgOut,'(A)') ''
            end if
          end do
        end if
        write(MsgOut,'(A)') ''
      end do

      write(MsgOut,'(A)') 'Setup_Enefunc_Restraint> Setup restraint functions'

      do i = 1, enefunc%num_restraintfuncs

        write(MsgOut,'(A,I5,A,I5)') &
          ' func  = ', i, ' kind  = ', enefunc%restraint_kind(i)

        ! summary for positional restraint
        !
        if (enefunc%restraint_kind(i) == RestraintsFuncPOSI .or. &
            enefunc%restraint_kind(i) == RestraintsFuncEM) then
          write(MsgOut,'(A,4F8.3)') &
            ' const(total, x, y, z) = ', enefunc%restraint_const(1:4,i) 
        else
          write(MsgOut,'(A,F8.3)') &
            ' const             = ',   enefunc%restraint_const(1,i)
          write(MsgOut,'(A,F8.3)') &
            ' ref               = ',   enefunc%restraint_ref(1,i)
        end if
 
        write(MsgOut,'(A,I5)') ' exponend of function = ',                &
                               enefunc%restraint_exponent_func(i)

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

          write(MsgOut,'(A,I5)') ' # of distances  = ', &
              enefunc%restraint_funcgrp(i)*(enefunc%restraint_funcgrp(i)-1)/2
          write(MsgOut,'(A)')    ' distancelist: '
          do j = 1, enefunc%restraint_funcgrp(i)
              do k = j + 1, enefunc%restraint_funcgrp(i)
                write(MsgOut,'(" group = (",i3,",",i3,")"$)')             &
                               enefunc%restraint_grouplist(j,i),          &
                               enefunc%restraint_grouplist(k,i)
                write(MsgOut,'(A)') ''
            end do
          end do

        ! summary for flat-bottom restraint
        !
        else if (enefunc%restraint_kind(i) == RestraintsFuncFB  .or.   &
                 enefunc%restraint_kind(i) == RestraintsFuncFBCOM ) then

          write(MsgOut,'(A,I5)') ' # of distances  = ', &
              enefunc%restraint_funcgrp(i)*(enefunc%restraint_funcgrp(i)-1)/2
          write(MsgOut,'(A)')    ' distancelist: '
          do j = 1, enefunc%restraint_funcgrp(i)
              do k = j + 1, enefunc%restraint_funcgrp(i)
                write(MsgOut,'(" group = (",i3,",",i3,")"$)')             &
                               enefunc%restraint_grouplist(j,i),          &
                               enefunc%restraint_grouplist(k,i)
                write(MsgOut,'(A)') ''
            end do
          end do

        ! summary for positinal, RMSD and principal component restraint
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

    return

  end subroutine setup_enefunc_restraints

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
    real(wp)                 :: totmass, coef
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
    enefunc%restraint_masscoef(1:max_atoms,1:num_groups) = 0.0_wp
    do i = 1, enefunc%num_restraintgroups
      totmass = 0.0_wp
      do j = 1, enefunc%restraint_numatoms(i) 
        jatom = enefunc%restraint_atomlist(j,i)
        totmass = totmass + molecule%mass(jatom)
      end do
      enefunc%restraint_totmass_group(i) = totmass
      coef = 1.0_wp/totmass

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
    integer                  :: i, ndata, itmd, ndata_max
    integer                  :: num_funcs, max_grp
    character(MaxLine)       :: string

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
      case(RestraintsFuncPC:RestraintsFuncPCCOM)
        if (ndata == 0) &
          call error_msg('Setup_Enefunc_Rest_Func> Error in index for PC')
      case(RestraintsFuncREPUL:RestraintsFuncREPULCOM)
        if (ndata <  2) &
          call error_msg('Setup_Enefunc_Rest_Func> Error in index for REPUL')
      case(RestraintsFuncFB:RestraintsFuncFBCOM)
        if (ndata <  2) &
          call error_msg('Setup_Enefunc_Rest_Func> Error in index for FB')
      case default
        if (ndata /= 1) & 
          call error_msg('Setup_Enefunc_Rest_Func> Error in index for '//  &
                         'POSI or RMSD or RG')
      end select

      num_indexes(i) = ndata
    end do

    ! get maximum number of groups among all 'index's
    !
    max_grp = maxval(num_indexes(1:num_funcs))
    enefunc%max_restraint_numgrps = max_grp
    !write(*,*)'enefunc%max_restraint_numgrps = ', enefunc%max_restraint_numgrps

    call alloc_enefunc(enefunc, EneFuncReff, num_funcs, max_grp)
    call alloc_enefunc(enefunc, EneFuncRefw, molecule%num_atoms)

    ndata_max = 0
    do i = 1, num_funcs
      string = restraints%constant(i)
      ndata = split_num(string)
      if (ndata_max < ndata) ndata_max = ndata
    end do
    call alloc_enefunc(enefunc, EneFuncRefr, num_funcs, ndata_max)
 
    ! setup parameters
    !
    do i = 1, num_funcs

      enefunc%restraint_funcgrp(i) = num_indexes(i)
      enefunc%restraint_kind(i)    = restraints%function(i)

      ! setup group list
      !
      call setup_restraint_grouplist(i, restraints, enefunc)

      ! setup exponent_func
      !
      call setup_restraint_exponent_func(i, restraints, enefunc)

      ! setup force constants
      !
      call setup_restraint_constants(i, restraints, enefunc)

      ! setup reference   (except for POSI)
      !
      call setup_restraint_reference(i, restraints, enefunc)

      ! setup exponent_dist (for DIST only)
      !
      call setup_restraint_exponent_dist(i, restraints, enefunc)

      ! setup weight_dist   (for DIST only)
      !
      call setup_restraint_weight_dist(i, restraints, enefunc)

      ! setup mode   (for PC only)
      !
      call setup_restraint_mode(i, restraints, enefunc)

      ! caging   (for RG only)
      !
      call setup_restraint_caging(i, restraints, enefunc)

      ! flat_radius   (for RG only)
      !
      call setup_restraint_flat_radius(i, restraints, enefunc)

    end do

    ! setup reference coordinates (for POSI, RMSD, and PC)
    !
    do i = 1, num_funcs
      if (enefunc%restraint_kind(i) == RestraintsFuncPOSI    .or. &
          enefunc%restraint_kind(i) == RestraintsFuncRMSD    .or. &
          enefunc%restraint_kind(i) == RestraintsFuncRMSDCOM .or. &
          enefunc%restraint_kind(i) == RestraintsFuncPC      .or. &
          enefunc%restraint_kind(i) == RestraintsFuncPCCOM) then
        call setup_enefunc_rest_refcoord(molecule, enefunc)
        exit
      end if
    end do

    ! setup principal component mode (for PC)
    !
    do i = 1, num_funcs
      if (enefunc%restraint_kind(i) == RestraintsFuncPC     .or. &
          enefunc%restraint_kind(i) == RestraintsFuncPCCOM) then
        call setup_enefunc_rest_mode(molecule, enefunc)
        exit
      end if
    end do

    ! setup work arrays (for RG)
    do i = 1, num_funcs
      if (enefunc%restraint_kind(i) == RestraintsFuncRG       .or. &
          enefunc%restraint_kind(i) == RestraintsFuncRGWOMASS) then
        if (.not. allocated(enefunc%restraint_masstmp)) then
          call setup_enefunc_rest_rg(molecule, enefunc)
        end if
        exit
      end if
    end do

    ! for steered & targeted MD
    !
    itmd = restraints%target_function
    if (itmd == 0) then
      enefunc%target_function = 0
      enefunc%steered_function = 0
    else
      if (itmd < 0 .or. itmd > num_funcs) &
        call error_msg('Setup_Enefunc_Rest_Func> Error in Targeted function')
      if (itmd > 0                                               .and. &
          enefunc%restraint_kind(itmd) /= RestraintsFuncDIST     .and. &
          enefunc%restraint_kind(itmd) /= RestraintsFuncDISTCOM  .and. &
          enefunc%restraint_kind(itmd) /= RestraintsFuncANGLE    .and. &
          enefunc%restraint_kind(itmd) /= RestraintsFuncANGLECOM .and. &
          enefunc%restraint_kind(itmd) /= RestraintsFuncDIHED    .and. &
          enefunc%restraint_kind(itmd) /= RestraintsFuncDIHEDCOM .and. &
          enefunc%restraint_kind(itmd) /= RestraintsFuncRMSD     .and. &
          enefunc%restraint_kind(itmd) /= RestraintsFuncRMSDCOM  .and. &
          enefunc%restraint_kind(itmd) /= RestraintsFuncREPUL    .and. &
          enefunc%restraint_kind(itmd) /= RestraintsFuncREPULCOM .and. &
          enefunc%restraint_kind(itmd) /= RestraintsFuncFB       .and. &
          enefunc%restraint_kind(itmd) /= RestraintsFuncFBCOM) then
        call error_msg('Setup_Enefunc_Rest_Func> Targeted function is not allowed')
      end if
      enefunc%target_function  = restraints%target_function
      enefunc%steered_function = restraints%target_function
      enefunc%target_value     = enefunc%restraint_ref(1,itmd)
    end if

    ! deallocate local array
    !
    deallocate(num_indexes)

    return

  end subroutine setup_enefunc_rest_func

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_restraint_grouplist
  !> @brief        setup restraint grouplist
  !! @authors      CK, TM
  !! @param[in]    ifunc      : number for function
  !! @param[in]    restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_restraint_grouplist(ifunc, restraints, enefunc)

    ! formal arguments
    integer,                 intent(in)    :: ifunc
    type(s_restraints),      intent(in)    :: restraints
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variable
    integer                  :: j, max_grp
    character(MaxLine)       :: string
    integer,     allocatable :: idata(:)


    max_grp = enefunc%max_restraint_numgrps
    allocate(idata(max_grp))

    string = restraints%select_index(ifunc)

    call split(split_num(string), max_grp, string, idata)
    do j = 1, max_grp 
      enefunc%restraint_grouplist(j,ifunc) = idata(j)
    end do

    deallocate(idata)

    return

  end subroutine setup_restraint_grouplist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_restraint_exponent_func
  !> @brief        setup restraint exponent func
  !! @authors      CK, TM
  !! @param[in]    ifunc      : number for function
  !! @param[in]    restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_restraint_exponent_func(ifunc, restraints, enefunc)

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

  end subroutine setup_restraint_exponent_func

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_restraint_constants
  !> @brief        setup restraint constants
  !! @authors      CK, TM
  !! @param[in]    ifunc      : number for function
  !! @param[in]    restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_restraint_constants(ifunc, restraints, enefunc)

    ! formal arguments
    integer,                 intent(in)    :: ifunc
    type(s_restraints),      intent(in)    :: restraints
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variable
    integer                  :: i
    integer                  :: ndata
    integer                  :: ndata_max
    integer                  :: num_funcs
    character(MaxLine)       :: string

    real(wp),    allocatable :: ddata(:)


    num_funcs = enefunc%num_restraintfuncs
    string = restraints%constant(ifunc)
    ndata = split_num(string)

    if (ndata == 0) then
      call error_msg('Setup_Restraint_Constants> constant is not given')
    end if

    ndata_max = size(enefunc%restraint_const_replica(:,1))
    allocate(ddata(ndata_max))

    call split(ndata, ndata_max, string, ddata)

    do i = 1, ndata
      enefunc%restraint_const_replica(i,ifunc) = ddata(i)
    end do

    deallocate(ddata)

    if (enefunc%restraint_kind(ifunc) == RestraintsFuncPOSI) then
      enefunc%restraint_const(1,ifunc) = enefunc%restraint_const_replica(1,ifunc)
      if (restraints%direction(ifunc) == RestraintsDirALL) then
        enefunc%restraint_const(2:4,ifunc) = 1.0_wp
      else
        enefunc%restraint_const(2:4,ifunc) = 0.0_wp
        enefunc%restraint_const(restraints%direction(ifunc),ifunc) = 1.0_wp
      end if
    else
      enefunc%restraint_const(1,ifunc) = enefunc%restraint_const_replica(1,ifunc)
      enefunc%restraint_const(2,ifunc) = enefunc%restraint_const_replica(1,ifunc)
      enefunc%restraint_const(3,ifunc) = 1.0_wp
      enefunc%restraint_const(4,ifunc) = 1.0_wp
      if (enefunc%restraint_kind(ifunc) == RestraintsFuncEM) then
        enefunc%restraint_const(2:4,ifunc) = 1.0_wp
      end if
    end if

    return

  end subroutine setup_restraint_constants

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_restraint_reference
  !> @brief        setup restraint reference
  !! @authors      CK, TM
  !! @param[in]    ifunc      : number for function
  !! @param[in]    restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_restraint_reference(ifunc, restraints, enefunc)

    ! formal arguments
    integer,                 intent(in)    :: ifunc
    type(s_restraints),      intent(in)    :: restraints
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variable
    integer                  :: i
    integer                  :: ndata
    integer                  :: ndata_max
    integer                  :: num_funcs
    character(MaxLine)       :: string

    real(wp),    allocatable :: ddata(:)


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

    if ((ndata > 0) .and. &
        (enefunc%restraint_kind(ifunc) ==  RestraintsFuncPOSI)) then
      call error_msg('Setup_Restraint_Reference> reference of POSI is not needed')
    end if

    if ((ndata == 0) .and. &
        (enefunc%restraint_kind(ifunc) /=  RestraintsFuncPOSI) .and. &
        (enefunc%restraint_kind(ifunc) /=  RestraintsFuncEM  )) then
      call error_msg('Setup_Restraint_Reference> reference is not given')
    end if

    ndata_max = size(enefunc%restraint_const_replica(:,1))

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

  end subroutine setup_restraint_reference

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_restraint_exponent_dist
  !> @brief        setup restraint exponent dist
  !! @authors      CK, TM
  !! @param[in]    ifunc      : number for function
  !! @param[in]    restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_restraint_exponent_dist(ifunc, restraints, enefunc)

    ! formal arguments
    integer,                 intent(in)    :: ifunc
    type(s_restraints),      intent(in)    :: restraints
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variable
    integer                  :: ndata, factor, maxfunc, max_grp
    character(MaxLine)       :: string

    integer,     allocatable :: idata(:)


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
        call error_msg('Setup_Restraint_Exponent_Dist> Error in exponent_dist')

      call split(ndata, maxfunc, string, idata)
      enefunc%restraint_exponent_dist(1:maxfunc,ifunc) = idata(1:maxfunc)
    end if
    deallocate(idata)

    return

  end subroutine setup_restraint_exponent_dist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_restraint_weight_dist
  !> @brief        setup restraint weight dist
  !! @authors      CK, TM
  !! @param[in]    ifunc      : number for function
  !! @param[in]    restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_restraint_weight_dist(ifunc, restraints, enefunc)

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
        call error_msg('Setup_Restraint_Weight_Dist> error in weight_dist')

      call split(ndata, maxfunc, string, ddata)
      enefunc%restraint_weight_dist(1:maxfunc,ifunc) = ddata(1:maxfunc)
    end if
    deallocate(ddata)

    return

  end subroutine setup_restraint_weight_dist


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_restraint_mode
  !> @brief        setup restraint mode (for PC)
  !! @authors      YK
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
  !  Subroutine    setup_restraint_caging
  !> @brief        setup restraint caging (for RG)
  !! @authors      MK
  !! @param[in]    ifunc      : number for function
  !! @param[in]    restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_restraint_caging(ifunc, restraints, enefunc)

    ! formal arguments
    integer,                 intent(in)    :: ifunc
    type(s_restraints),      intent(in)    :: restraints
    type(s_enefunc),         intent(inout) :: enefunc

    integer      :: ndata
    character(6) :: cdata(1)


    cdata(1) = ''
    if ( enefunc%restraint_kind(ifunc) /= RestraintsFuncRG .and. &
         enefunc%restraint_kind(ifunc) /= RestraintsFuncRGWOMASS ) then
      return
    end if

    ndata = split_num(restraints%caging(ifunc))
    if (ndata > 0) then
      call split(1, 1, restraints%caging(ifunc), cdata )
    end if

    call tolower(cdata(1))
    enefunc%restraint_caging(ifunc) = (cdata(1) .eq. 'yes' .or. &
                                       cdata(1) .eq. 'true' )

    return

  end subroutine setup_restraint_caging

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_restraint_flat_radius
  !> @brief        setup restraint flat radius (for RG)
  !! @authors      MK
  !! @param[in]    ifunc      : number for function
  !! @param[in]    restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_restraint_flat_radius(ifunc, restraints, enefunc)

    ! formal arguments
    integer,                 intent(in)    :: ifunc
    type(s_restraints),      intent(in)    :: restraints
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: ndata


    if (enefunc%restraint_kind(ifunc) /= RestraintsFuncRG .and. &
        enefunc%restraint_kind(ifunc) /= RestraintsFuncRGWOMASS) then
      return
    end if
    if (len_trim(restraints%flat_radius(ifunc)) == 0) then
      enefunc%restraint_flat_radius(ifunc) = 0.0_wp
      return
    end if
    read(restraints%flat_radius(ifunc),*) enefunc%restraint_flat_radius(ifunc)

    return

  end subroutine setup_restraint_flat_radius

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
  !  Subroutine    setup_enefunc_rest_mode
  !> @brief        setup enefunc principal component mode
  !! @authors      YK
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

    enefunc%num_pc_modes = molecule%num_pc_modes

    call alloc_enefunc(enefunc, EneFuncMode, molecule%num_pc_modes)
    call alloc_enefunc(enefunc, EneFuncModeRef, molecule%num_atoms)

    enefunc%pc_mode(:) = molecule%pc_mode(:)

    return

  end subroutine setup_enefunc_rest_mode

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_rest_rg
  !> @brief        setup enefunc for Rg calculation
  !! @authors      MK
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_rest_rg(molecule, enefunc)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc


    enefunc%num_atoms_ref = molecule%num_atoms
    call alloc_enefunc(enefunc, EneFuncRefc, molecule%num_atoms)

    return

  end subroutine setup_enefunc_rest_rg

end module at_enefunc_restraints_mod
