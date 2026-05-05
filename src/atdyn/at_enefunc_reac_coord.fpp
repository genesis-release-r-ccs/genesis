!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_enefunc_reac_coord_mod
!> @brief   reaction coordinate energy functions
!! @authors Junichiro Koga (JK), Shingo Ito (SI)
!
!  (c) Copyright 2020 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8
  
#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_enefunc_reac_coord_mod

  use at_restraints_str_mod
  use at_reac_coord_str_mod
  use at_enefunc_str_mod
  use at_remd_str_mod
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

  public   :: setup_enefunc_reac_coord
  public   :: reac_coords_are_out_of_range
  private  :: setup_enefunc_reac_group
  private  :: setup_enefunc_reac_func
  private  :: setup_reac_coord_grouplist
  private  :: setup_reac_coord_exponent_dist
  private  :: setup_reac_coord_weight_dist
  private  :: setup_reac_coord_unitv
  private  :: setup_enefunc_reac_refcoord
  private  :: setup_enefunc_reac_reac_coord
  private  :: alloc_enefunc_reac_coord

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_reac_coord
  !> @brief        define reaction coordinate potential
  !! @authors      JK
  !! @param[in]    molecule   : molecule information
  !! @param[in]    reac_coord : reaction coordinate information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_reac_coord(molecule, reac_coord, enefunc)

    ! formal arguments
    type(s_molecule),   intent(in)    :: molecule
    type(s_reac_coord), intent(in)    :: reac_coord
    type(s_enefunc),    intent(inout) :: enefunc
    integer :: i, j, k
    
    enefunc%reac_coord%reac_coord_flag = reac_coord%reac_coord_flag
    if( .not. reac_coord%reac_coord_flag ) return

    enefunc%reac_coord%nreac_coord = reac_coord%nreac_coord
    enefunc%reac_coord%out_of_range = .false.

    ! setup group
    !
    call setup_enefunc_reac_group(molecule, reac_coord, enefunc)

    ! setup func
    !
    call setup_enefunc_reac_func(molecule, reac_coord, enefunc)

    ! setup domain
    !
    call setup_enefunc_reac_domain(enefunc)

    ! setup reaction coordinate mesh
    call setup_enefunc_reac_reac_coord(molecule, reac_coord, enefunc)

    ! we do not allow duplication option with reaction coordinate
    !
    !if (domain%num_duplicate > 1) &
    !  call error_msg('Setup_Enefunc_ReacCoord> duplication is not allowed with reaction coordinate')

    ! summary of setup enefunc_reac_coord
    !
    if (main_rank) then

      write(MsgOut,'(A)') 'Setup_Enefunc_ReacCoord> Setup reaction coordinate groups'

      do i = 1, enefunc%reac_coord%num_groups
        write(MsgOut,'(A,I5,3A)') ' group = ', i, ', "', trim(reac_coord%group(i)),'"'
        write(MsgOut,'(A,I5)')   ' # of atoms = ', enefunc%reac_coord%numatoms(i)
        write(MsgOut,'(A)')      ' atomlist: ' 
        do j = 1, enefunc%reac_coord%numatoms(i) 
          write(MsgOut,'(i7,$)') enefunc%reac_coord%atomlist(j,i)
          if (mod(j,10) == 0 .and. j /= enefunc%reac_coord%numatoms(i)) then
            write(MsgOut,'(A)') ''
          end if
        end do
        write(MsgOut,'(A)') ''
      end do

      write(MsgOut,'(A)') 'Setup_Enefunc_ReacCoord> Setup reaction coordinate functions'

      do i = 1, enefunc%reac_coord%nreac_coord

        write(MsgOut,'(A,I5,A,I5)') &
          ' func  = ', i, ' kind  = ', enefunc%reac_coord%kind(i)

        ! summary for distance reaction coordinate
        !
        if (enefunc%reac_coord%kind(i) == ReacCoordsFuncDIST .or.          &
            enefunc%reac_coord%kind(i) == ReacCoordsFuncDISTCOM) then

          write(MsgOut,'(A,I5)') ' # of distances  = ', enefunc%reac_coord%funcgrp(i)/2
          write(MsgOut,'(A)')    ' distancelist: '
          do j = 1, enefunc%reac_coord%funcgrp(i)/2
            write(MsgOut,'(" group = (",i3,",",i3,")"$)')                 &
                           enefunc%reac_coord%grouplist(2*j-1,i),          &
                           enefunc%reac_coord%grouplist(2*j,i)
            write(MsgOut,'(" weight = ",f8.3," exponent = ",i5,$)')       &
                           enefunc%reac_coord%weight_dist(j,i),            &
                           enefunc%reac_coord%exponent_dist(j,i)
            write(MsgOut,'(A)') ''
          end do

        ! summary for angle reaction coordinate
        !
        else if (enefunc%reac_coord%kind(i) == ReacCoordsFuncANGLE .or.    &
                 enefunc%reac_coord%kind(i) == ReacCoordsFuncANGLECOM ) then

          write(MsgOut,'(A,I5)') ' # of angles     = ', enefunc%reac_coord%funcgrp(i)/3
          write(MsgOut,'(A)')    ' anglelist   : '
          do j = 1, enefunc%reac_coord%funcgrp(i)/3
            write(MsgOut,'(" group = (",i3,",",i3,",",i3,")"$)')          &
                           enefunc%reac_coord%grouplist(3*j-2,i),          &
                           enefunc%reac_coord%grouplist(3*j-1,i),          &
                           enefunc%reac_coord%grouplist(3*j,i)
            write(MsgOut,'(A)') ''
          end do

        ! summary for dihedral angle reaction coordinate
        !
        else if (enefunc%reac_coord%kind(i) == ReacCoordsFuncDIHED  .or.   &
                 enefunc%reac_coord%kind(i) == ReacCoordsFuncDIHEDCOM ) then

          write(MsgOut,'(A,I5)') ' # of dihedrals  = ', enefunc%reac_coord%funcgrp(i)/4
          write(MsgOut,'(A)')    ' dihedrallist: ' 
          do j = 1, enefunc%reac_coord%funcgrp(i)/4
            write(MsgOut,'(" group = (",i3,",",i3,",",i3,",",i3,")"$)')   &
                           enefunc%reac_coord%grouplist(4*j-3,i),          &
                           enefunc%reac_coord%grouplist(4*j-2,i),          &
                           enefunc%reac_coord%grouplist(4*j-1,i),          &
                           enefunc%reac_coord%grouplist(4*j,i)
            write(MsgOut,'(A)') ''
          end do

        ! summary for positinal and RMSD reaction coordinate
        !
        else 

          write(MsgOut,'(A,I5)') ' # of groups  = ', enefunc%reac_coord%funcgrp(i)
          write(MsgOut,'(" grouplist: ",$)') 
          do j = 1, enefunc%reac_coord%funcgrp(i) 
            write(MsgOut,'(i3,$)') enefunc%reac_coord%grouplist(j,i)
            if (mod(j,20) == 0 .and. j /= enefunc%reac_coord%funcgrp(i) )  &
              write(MsgOut,'(A)') ''
          end do
          write(MsgOut,'(A)') ''

        end if
        write(MsgOut,'(A)') ''

        write(MsgOut,'(A,I5)')    ' nxi          = ',enefunc%reac_coord%nxi      (i)
        write(MsgOut,'(A,F10.5)') ' ximin        = ',enefunc%reac_coord%ximin    (i)
        write(MsgOut,'(A,F10.5)') ' ximax        = ',enefunc%reac_coord%ximax    (i)
        write(MsgOut,'(A,F10.5)') ' dxi          = ',enefunc%reac_coord%dxi      (i)
        write(MsgOut,'(A,F10.5)') ' invdxi       = ',enefunc%reac_coord%invdxi   (i)
        write(MsgOut,'(A,F10.5)') ' uboundary    = ',enefunc%reac_coord%uboundary(i)
        write(MsgOut,'(A,L)')     ' periodicity  = ',enefunc%reac_coord%periodicity(i)

      end do
    endif

    return
  end subroutine setup_enefunc_reac_coord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_reac_group
  !> @brief        define reaction coordinate atom group
  !! @authors      JK
  !! @param[in]    molecule : molecule information
  !! @param[in]    reac     : reaction coordinate information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_reac_group(molecule, reac, enefunc)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_reac_coord),      intent(in)    :: reac
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    real(dp)                 :: totmass, coef
    integer                  :: num_groups, max_atoms
    integer                  :: i, j, jatom


    enefunc%reac_coord%monitor_period = reac%monitor_period
    enefunc%reac_coord%num_groups     = reac%num_groups
    enefunc%reac_coord%max_atoms      = reac%max_atoms
    num_groups                        = enefunc%reac_coord%num_groups
    max_atoms                         = enefunc%reac_coord%max_atoms
    call alloc_enefunc_reac_coord(enefunc%reac_coord, EneFuncReacg, &
                                  num_groups, max_atoms)

    ! setup atom lists and num atoms
    !
    do j = 1, num_groups
      enefunc%reac_coord%numatoms(j) = reac%num_atoms(j)
      do i = 1, max_atoms
        enefunc%reac_coord%atomlist(i,j) = reac%atomlist(i,j)
      end do
    end do

    ! check number of atoms in groups
    !
    do i = 1, enefunc%reac_coord%num_groups
      if (enefunc%reac_coord%numatoms(i) == 0) &
        call error_msg('Setup_Enefunc_Reac_Group> The specified atom group '//&
	               'is empty')
    end do

    ! setup mass_coefficient
    !
    enefunc%reac_coord%masscoef(1:max_atoms,1:num_groups) = 0.0_dp
    do i = 1, enefunc%reac_coord%num_groups
      totmass = 0.0_dp
      do j = 1, enefunc%reac_coord%numatoms(i) 
        jatom = enefunc%reac_coord%atomlist(j,i)
        totmass = totmass + molecule%mass(jatom)
      end do
      coef = 1.0_dp/totmass

      do j = 1, enefunc%reac_coord%numatoms(i) 
        jatom = enefunc%reac_coord%atomlist(j,i)
        enefunc%reac_coord%masscoef(j,i) = molecule%mass(jatom) * coef
      end do
    end do

    return

  end subroutine setup_enefunc_reac_group

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_reac_func
  !> @brief        define reaction coordinate function
  !! @authors      CK, TM
  !! @param[in]    molecule   : molecule information
  !! @param[in]    reac_coord : reaction coordinate information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_reac_func(molecule, reac_coord, enefunc)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_reac_coord),      intent(in)    :: reac_coord
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, ndata
    integer                  :: num_funcs, max_grp
    integer,     allocatable :: num_indexes(:)


    ! get total number of groups in index
    !

    num_funcs = enefunc%reac_coord%nreac_coord
    allocate(num_indexes(1:num_funcs))

    do i = 1, num_funcs
      ndata = split_num(reac_coord%select_index(i))

      select case(reac_coord%function(i))
      case(ReacCoordsFuncDIST:ReacCoordsFuncDISTCOM)
        if (ndata <  2 .or. mod(ndata,2) /= 0) &
          call error_msg('Setup_Enefunc_ReacCoord_Func> Invalid atom index '//&
	                 'specified for DIST')

      case(ReacCoordsFuncANGLE:ReacCoordsFuncANGLECOM)
        if (ndata /= 3) &
          call error_msg('Setup_Enefunc_ReacCoord_Func> Invalid atom index '//&
	                 'specified for ANGLE')

      case(ReacCoordsFuncDIHED:ReacCoordsFuncDIHEDCOM)
        if (ndata /= 4) &
          call error_msg('Setup_Enefunc_ReacCoord_Func> Invalid atom index '//&
	                 'specified for DIHE')
      case(ReacCoordsFuncPOSI, ReacCoordsFuncPOSICOM)
        if (ndata == 0) &
          call error_msg('Setup_Enefunc_ReacCoord_Func> Invalid atom index '//&
	                 'specified for POSI')
      case default
        if (ndata /= 1) & 
          call error_msg('Setup_Enefunc_ReacCoord_Func> Invalid atom index '//&
	                 'specified for POSI or RMSD or RG')

      end select

      num_indexes(i) = ndata
    end do

    ! get maximum number of groups among all 'index's
    !
    max_grp = maxval(num_indexes(1:num_funcs))
    enefunc%reac_coord%max_num_groups = max_grp

    call alloc_enefunc_reac_coord(enefunc%reac_coord, EneFuncReacf, num_funcs, max_grp)
 
    ! setup parameters
    !
    do i = 1, num_funcs

      enefunc%reac_coord%funcgrp(i) = num_indexes(i)
      enefunc%reac_coord%kind(i)    = reac_coord%function(i)

      ! setup group list
      !
      call setup_reac_coord_grouplist(i, reac_coord, enefunc)

      ! setup exponent_dist (for DIST only)
      !
      call setup_reac_coord_exponent_dist(i, reac_coord, enefunc)

      ! setup weight_dist   (for DIST only)
      !
      call setup_reac_coord_weight_dist(i, reac_coord, enefunc)

      ! setup direction (for POSI only)
      !
      call setup_reac_coord_direction(i, reac_coord, enefunc)

      ! setup weight_dist   (for POSI only)
      !
      call setup_reac_coord_unitv(i, reac_coord, enefunc)

    end do

    ! deallocate local array
    !
    deallocate(num_indexes)

    return

  end subroutine setup_enefunc_reac_func

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_reac_domain
  !> @brief        setup reaction coordinate information about domain
  !! @authors      JK
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_reac_domain(enefunc)

    ! formal arguments
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    integer                  :: i, ig, j, k, ix, iatm, icel
    logical                  :: position_flag
    logical                  :: rmsd_flag

    real(wp),        pointer :: const(:,:), ref(:,:)
    real(wp),        pointer :: reac_coord_force(:,:)
    real(wp),        pointer :: reac_coord_coord(:,:)
    integer,         pointer :: num_funcs, num_atoms(:)
    integer,         pointer :: grouplist(:,:), atomlist(:,:), bondslist(:,:,:)
    integer,         pointer :: bonds_to_atom(:,:)
    integer,         pointer :: reac_coord_atom(:,:)
    integer                  :: iar(4)
    integer                  :: dummy 


    num_funcs       => enefunc%reac_coord%nreac_coord
    num_atoms       => enefunc%reac_coord%numatoms
    grouplist       => enefunc%reac_coord%grouplist
    atomlist        => enefunc%reac_coord%atomlist
    position_flag=.false.
    rmsd_flag=.false.

    ! check the size for distance reaction coordinate
    !
    do i = 1, num_funcs
      if (enefunc%reac_coord%kind(i) == ReacCoordsFuncDIST     .or. &
          enefunc%reac_coord%kind(i) == ReacCoordsFuncDISTCOM  .or. &
          enefunc%reac_coord%kind(i) == ReacCoordsFuncANGLE    .or. &
          enefunc%reac_coord%kind(i) == ReacCoordsFuncANGLECOM .or. &
          enefunc%reac_coord%kind(i) == ReacCoordsFuncDIHED    .or. &
          enefunc%reac_coord%kind(i) == ReacCoordsFuncDIHEDCOM .or. &
          enefunc%reac_coord%kind(i) == ReacCoordsFuncPOSI     .or. &
          enefunc%reac_coord%kind(i) == ReacCoordsFuncPOSICOM ) then
        k = 0
        do ig = 1, enefunc%reac_coord%funcgrp(i)
          do ix = 1, num_atoms(grouplist(ig,i))
            k = k + 1
          end do
        end do
        enefunc%reac_coord%num_atoms_bonds_reac_coord(i) = k
      end if
    end do

    ! memory allocation for reaction coordinate
    !
    !            maxval(enefunc%reac_coord%num_atoms_bonds_reac_coord), &
    !            enefunc%reac_coord%nreac_coord
    dummy = 0
    call alloc_enefunc_reac_coord(enefunc%reac_coord, EneFuncReacatm, &
         dummy, maxval(enefunc%reac_coord%num_atoms_bonds_reac_coord))

    reac_coord_atom  => enefunc%reac_coord%atom
    reac_coord_force => enefunc%reac_coord%force
    reac_coord_coord => enefunc%reac_coord%coord

    if (position_flag .and. rmsd_flag) &
     call error_msg('Setup_Enefunc_Reac_Domain> RMSD and POSI'//&
                    'cannot be used simultaneously.')

    bondslist       => enefunc%reac_coord%bondslist
    bonds_to_atom   => enefunc%reac_coord%bondslist_to_atomlist
    ! distance reaction coordinate
    !
    do i = 1, num_funcs
      if (enefunc%reac_coord%kind(i) == ReacCoordsFuncDIST     .or. &
          enefunc%reac_coord%kind(i) == ReacCoordsFuncDISTCOM  .or. &
          enefunc%reac_coord%kind(i) == ReacCoordsFuncANGLE    .or. &
          enefunc%reac_coord%kind(i) == ReacCoordsFuncANGLECOM .or. &
          enefunc%reac_coord%kind(i) == ReacCoordsFuncDIHED    .or. &
          enefunc%reac_coord%kind(i) == ReacCoordsFuncDIHEDCOM .or. &
          enefunc%reac_coord%kind(i) == ReacCoordsFuncPOSI     .or. &
          enefunc%reac_coord%kind(i) == ReacCoordsFuncPOSICOM ) then

        k = 0
        do ig = 1, enefunc%reac_coord%funcgrp(i)
          do ix = 1, num_atoms(grouplist(ig,i))
            k = k + 1
            bonds_to_atom(k,i)              = atomlist(ix,grouplist(ig,i))
            bondslist(ix,grouplist(ig,i),i) = k
          end do
        end do

      end if
    end do

    return

  end subroutine setup_enefunc_reac_domain

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_reac_coord_grouplist
  !> @brief        setup reaction coordinate grouplist
  !! @authors      JK
  !! @param[in]    ifunc      : number for function
  !! @param[in]    reac_coord : reaction coordinate information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_reac_coord_grouplist(ifunc, reac_coord, enefunc)

    ! formal arguments
    integer,                 intent(in)    :: ifunc
    type(s_reac_coord),      intent(in)    :: reac_coord
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variable
    integer                       :: j, max_grp
    character(MaxLine)            :: string
    integer,          allocatable :: idata(:)


    max_grp = enefunc%reac_coord%max_num_groups
    allocate(idata(max_grp))

    string = reac_coord%select_index(ifunc)

    call split(split_num(string), max_grp, string, idata)
    do j = 1, max_grp 
      enefunc%reac_coord%grouplist(j,ifunc) = idata(j)
    end do

    deallocate(idata)

    return

  end subroutine setup_reac_coord_grouplist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_reac_coord_exponent_dist
  !> @brief        setup reaction coordinate exponent dist
  !! @authors      CK, TM
  !! @param[in]    ifunc      : number for function
  !! @param[in]    reac_coord : reaction coordinate information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_reac_coord_exponent_dist(ifunc, reac_coord, enefunc)

    ! formal arguments
    integer,                 intent(in)    :: ifunc
    type(s_reac_coord),      intent(in)    :: reac_coord
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variable
    integer                       :: ndata, factor, maxfunc, max_grp
    character(MaxLine)            :: string
    integer,          allocatable :: idata(:)


    select case (enefunc%reac_coord%kind(ifunc))
    case(ReacCoordsFuncDIST:ReacCoordsFuncDISTCOM)
      factor = 2
    case default
      return
    end select

    max_grp = enefunc%reac_coord%max_num_groups
    maxfunc = max(int(max_grp/factor), 1)

    ! initialize
    !
    enefunc%reac_coord%exponent_dist(1:maxfunc,ifunc) = 1

    ! set exponent_dist
    !
    string = reac_coord%exponent_dist(ifunc)
    ndata = split_num(string)

    allocate(idata(1:maxfunc))
    if (ndata > 0) then
      if (ndata /= enefunc%reac_coord%funcgrp(ifunc)/factor) &
        call error_msg('Setup_ReacCoord_Exponent_Dist> Invalid value of '//&
	               'exponent_dist')

      call split(ndata, maxfunc, string, idata)
      enefunc%reac_coord%exponent_dist(1:maxfunc,ifunc) = idata(1:maxfunc)
    end if
    deallocate(idata)

    return

  end subroutine setup_reac_coord_exponent_dist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_reac_coord_weight_dist
  !> @brief        setup reaction coordinate weight dist
  !! @authors      JK
  !! @param[in]    ifunc      : number for function
  !! @param[in]    reac_coord : reaction coordinate information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_reac_coord_weight_dist(ifunc, reac_coord, enefunc)

    ! formal arguments
    integer,                 intent(in)    :: ifunc
    type(s_reac_coord),      intent(in)    :: reac_coord
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variable
    integer                  :: ndata, factor, maxfunc, max_grp
    character(MaxLine)       :: string
    real(wp),    allocatable :: ddata(:)


    select case (enefunc%reac_coord%kind(ifunc))
    case(ReacCoordsFuncDIST:ReacCoordsFuncDISTCOM)
      factor = 2
    case default
      return
    end select

    max_grp = enefunc%reac_coord%max_num_groups
    maxfunc = max(int(max_grp/factor), 1)

    ! initialize
    !
    enefunc%reac_coord%weight_dist(1:maxfunc,ifunc) = 1.0_wp

    ! set weight_dist
    !
    string = reac_coord%weight_dist(ifunc)
    ndata = split_num(string)

    allocate(ddata(1:maxfunc))
    if (ndata > 0) then
      if (ndata /= enefunc%reac_coord%funcgrp(ifunc)/factor) &
        call error_msg('Setup_ReacCoord_Weight_Dist> Invalid value of weight_dist')

      call split(ndata, maxfunc, string, ddata)
      enefunc%reac_coord%weight_dist(1:maxfunc,ifunc) = ddata(1:maxfunc)
    end if
    deallocate(ddata)

    return

  end subroutine setup_reac_coord_weight_dist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_reac_coord_unitv
  !> @brief        setup reaction coordinate unit vector
  !! @authors      JK
  !! @param[in]    ifunc      : number for function
  !! @param[in]    reac_coord : reaction coordinate information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_reac_coord_unitv(ifunc, reac_coord, enefunc)

    ! formal arguments
    integer,                 intent(in)    :: ifunc
    type(s_reac_coord),      intent(in)    :: reac_coord
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variable
    integer                  :: ndim, ndata, factor, maxfunc, max_grp
    character(MaxLine)       :: string
    real(wp),    allocatable :: ddata(:)
    real(wp)                 :: norm


    select case (enefunc%reac_coord%kind(ifunc))
    case(ReacCoordsFuncPOSI:ReacCoordsFuncPOSICOM)
      ndim = 3
    case default
      return
    end select

    max_grp = enefunc%reac_coord%max_num_groups
    maxfunc = max_grp

    ! initialize
    !
    enefunc%reac_coord%unitv = 1.0_wp
    enefunc%reac_coord%unitv(1,ifunc) = 1.0_wp

    string = reac_coord%unitv(ifunc)
    ndata = split_num(string)
    if(ndata /= ndim) then
      call error_msg('Setup_ReacCoordUnitV> The unit vector must be defined '//&
	             'by exactly three real numbers (x, y, z).')
    endif

    allocate(ddata(1:3))
    call split(ndata, 3, string, ddata)
    enefunc%reac_coord%unitv(1:3, ifunc) = ddata(1:3)
    norm = sqrt(dot_product(enefunc%reac_coord%unitv(:,ifunc), enefunc%reac_coord%unitv(:,ifunc)))
    if(norm.le.EPS) then
      call error_msg('Setup_ReacCoordUnitV> The norm of the provided unit '//&
	             'vector is near zero.')
    endif

    enefunc%reac_coord%unitv(1:3, ifunc) = enefunc%reac_coord%unitv(1:3, ifunc)/norm
    deallocate(ddata)

    return

  end subroutine setup_reac_coord_unitv

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_reac_coord_direction
  !> @brief        setup reaction coordinate direction
  !! @authors      JK
  !! @param[in]    ifunc      : number for function
  !! @param[in]    reac_coord : reaction coordinate information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_reac_coord_direction(ifunc, reac_coord, enefunc)

    ! formal arguments
    integer,                 intent(in)    :: ifunc
    type(s_reac_coord),      intent(in)    :: reac_coord
    type(s_enefunc),         intent(inout) :: enefunc

    if (enefunc%reac_coord%kind(ifunc) == ReacCoordsFuncPOSI     .or. &
        enefunc%reac_coord%kind(ifunc) == ReacCoordsFuncPOSICOM  .or. &
        enefunc%reac_coord%kind(ifunc) == ReacCoordsFuncDIST     .or. &
        enefunc%reac_coord%kind(ifunc) == ReacCoordsFuncDISTCOM) then
      enefunc%reac_coord%direction(ifunc) = reac_coord%direction(ifunc)
    endif

    select case (enefunc%reac_coord%direction(ifunc))

    case(1) !! ALL
      enefunc%reac_coord%const(1:3,ifunc) = 1.0_wp

    case(2, 5) !! X or X2
      enefunc%reac_coord%const(1,ifunc)   = 1.0_wp
      enefunc%reac_coord%const(2:3,ifunc) = 0.0_wp

    case(3, 6) !! Y or Y2
      enefunc%reac_coord%const(1,ifunc)   = 0.0_wp
      enefunc%reac_coord%const(2,ifunc)   = 1.0_wp
      enefunc%reac_coord%const(3,ifunc)   = 0.0_wp

    case(4, 7) !! Z or Z2
      enefunc%reac_coord%const(1:2,ifunc) = 0.0_wp
      enefunc%reac_coord%const(3,ifunc)   = 1.0_wp

    end select
      
    return

  end subroutine setup_reac_coord_direction

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_reac_refcoord
  !> @brief        setup enefunc reaction coordinate reference coords
  !! @authors      JK
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_reac_refcoord(molecule, enefunc)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: j


!    if (size(molecule%atom_refcoord(1,:)) /= molecule%num_atoms) &
!      call error_msg('Setup_Enefunc_Reac_Refcoord> bad reffile in [INPUT]')

!    enefunc%reac_coord%num_atoms_ref = molecule%num_atoms

    call alloc_enefunc_reac_coord(enefunc%reac_coord, EneFuncRefc, molecule%num_atoms)

!    do j = 1, molecule%num_atoms
!      enefunc%reac_coord%refcoord(1,j) = molecule%atom_refcoord(1,j)
!      enefunc%reac_coord%refcoord(2,j) = molecule%atom_refcoord(2,j)
!      enefunc%reac_coord%refcoord(3,j) = molecule%atom_refcoord(3,j)
!    end do

    return

  end subroutine setup_enefunc_reac_refcoord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_reac_reac_coord
  !> @brief        define reaction coordinate 
  !! @authors      JK
  !! @param[in]    reac     : reaction coordinate information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_reac_reac_coord(molecule, reac, enefunc)
    type(s_molecule),         intent(in)    :: molecule
    type(s_reac_coord),       intent(in)    :: reac
    type(s_enefunc), target,  intent(inout) :: enefunc

    integer,  pointer  :: nxi      (:)
    integer,  pointer  :: direct   (:)
    real(wp), pointer  :: ximin    (:)
    real(wp), pointer  :: ximax    (:)
    real(wp), pointer  :: dxi      (:)
    real(wp), pointer  :: invdxi   (:)
    real(wp), pointer  :: uboundary(:)
    logical,  pointer  :: periodicity(:)
    logical,  pointer  :: chk_periodicity

    integer            :: ireac, i

    call alloc_enefunc_reac_coord(enefunc%reac_coord, EneFuncReacCoord, &
                                  reac%nreac_coord, &
                          maxval(enefunc%reac_coord%num_atoms_bonds_reac_coord), &
                                  molecule%num_atoms)

    enefunc%reac_coord%mol_num_atoms = molecule%num_atoms                              
    nxi         => enefunc%reac_coord%nxi 
    ximin       => enefunc%reac_coord%ximin
    ximax       => enefunc%reac_coord%ximax 
    dxi         => enefunc%reac_coord%dxi
    invdxi      => enefunc%reac_coord%invdxi
    uboundary   => enefunc%reac_coord%uboundary
    periodicity => enefunc%reac_coord%periodicity
    chk_periodicity => enefunc%reac_coord%chk_periodicity

    chk_periodicity = .false.
    do ireac = 1, reac%nreac_coord
      ximin    (ireac)   = reac%ximin    (ireac)
      ximax    (ireac)   = reac%ximax    (ireac)
      nxi      (ireac)   = reac%nxi      (ireac)
      dxi      (ireac)   = reac%dxi      (ireac)
      uboundary(ireac)   = reac%uboundary(ireac)
      periodicity(ireac) = reac%periodicity(ireac)
  
      if(periodicity(ireac)) chk_periodicity = .true.

      if(ximax(ireac) < ximin(ireac)) &
      & call error_msg('Setup_Enefunc_ReacCoord> ximax must be greater than ximin.')
      if(nxi(ireac) > 1) then
        dxi   (ireac) = (ximax(ireac)-ximin(ireac))/real(nxi(ireac)-1,wp)
        invdxi(ireac) = 1.0_wp/dxi(ireac)
      else if (dxi(ireac) > 1.e-5) then
        nxi   (ireac) = nint((ximax(ireac)-ximin(ireac))/dxi(ireac))+1
        dxi   (ireac) = (ximax(ireac)-ximin(ireac))/real(nxi(ireac)-1,wp)
        invdxi(ireac) = 1.0_wp/dxi(ireac)
      else
        call error_msg('Setup_Enefunc_ReacCoord> Either nxi or dxi must be '//&
		       'specified.')
      endif
      if(uboundary(ireac)<0) &
       call error_msg('uboundary must be a non-negative value.')

      invdxi(ireac) =  1.0_wp/dxi(ireac)
    enddo

    enefunc%reac_coord%nreac_coord_comb = 1
    do i = 1, enefunc%reac_coord%nreac_coord
      enefunc%reac_coord%nreac_coord_comb = &
      enefunc%reac_coord%nreac_coord_comb * (enefunc%reac_coord%nxi(i)+2)
    enddo

  end subroutine setup_enefunc_reac_reac_coord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_enefunc_reac_coord
  !> @brief        allocate energy functions information for the reaction
  !                coordinate
  !! @authors      JK, SI
  !! @param[out]   enefunc   : potential energy functions information
  !! @param[in]    variable  : selected variable
  !! @param[in]    var_size  : size of the selected variable
  !! @param[in]    var_size1 : 2nd size of the selected variable (optional)
  !! @param[in]    var_size2 : 3rd size of the selected variable (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_enefunc_reac_coord(reac_coord, variable, var_size, &
                                      var_size1, var_size2)
    ! formal arguments
    type(s_enefunc_reac_coord), intent(inout) :: reac_coord
    integer,                    intent(in)    :: variable
    integer,                    intent(in)    :: var_size
    integer,          optional, intent(in)    :: var_size1
    integer,          optional, intent(in)    :: var_size2

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat
    integer                  :: var_size3, var_size4
    integer                  :: nreac

    alloc_stat   = 0
    dealloc_stat = 0
    nreac        = reac_coord%nreac_coord

    ! allocate selected variables
    !
    ! var_size  = num_groups
    ! var_size1 = max_atoms
    !
    select case (variable)
    case (EneFuncReacg)

      if (allocated(reac_coord%numatoms)) then
        if (size(reac_coord%numatoms) /= var_size) &
          deallocate(reac_coord%numatoms,  &
                     reac_coord%atomlist,  &
                     reac_coord%bondslist, &
                     reac_coord%masscoef,  &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(reac_coord%numatoms)) &
        allocate(reac_coord%numatoms (var_size),                        &
                 reac_coord%atomlist (1:var_size1, 1:var_size),         &
                 reac_coord%bondslist(1:var_size1, 1:var_size, 1:nreac),&
                 reac_coord%masscoef (1:var_size1, 1:var_size),         &
                 stat = alloc_stat)

      reac_coord%numatoms (1:var_size)                       = 0
      reac_coord%atomlist (1:var_size1, 1:var_size)          = 0
      reac_coord%bondslist(1:var_size1, 1:var_size, 1:nreac) = 0
      reac_coord%masscoef (1:var_size1, 1:var_size)          = 0.0_wp

    case (EneFuncReacf)
      ! var_size  = num_funcs
      ! var_size1 = max_grp
      !
      if (allocated(reac_coord%kind)) then
        if (size(reac_coord%kind) /= var_size)              &
          deallocate(reac_coord%kind,                       &
                     reac_coord%grouplist,                  &
                     reac_coord%funcgrp,                    &
                     reac_coord%exponent_dist,              &
                     reac_coord%unitv,                      &
                     reac_coord%weight_dist,                &
                     reac_coord%direction,                  &
                     reac_coord%const,                      &
                     reac_coord%num_atoms_bonds_reac_coord, &
                     reac_coord%indices,                    &
                     reac_coord%indices2,                   &
                     reac_coord%istarts,                    &
                     reac_coord%nnonzeros,                  &
                     stat = dealloc_stat)
      end if

      var_size3 = max(int(var_size1/2),1)
      var_size4 = int(var_size1*(var_size1-1)/2)

      if (.not. allocated(reac_coord%kind)) &
        allocate(reac_coord%kind         (var_size),                &
                 reac_coord%grouplist    (1:var_size1, 1:var_size), &
                 reac_coord%funcgrp      (1:var_size),              &
                 reac_coord%exponent_dist(1:var_size3, 1:var_size), &
                 reac_coord%weight_dist  (1:var_size3, 1:var_size), &
                 reac_coord%unitv        (1:3,         1:var_size), &
                 reac_coord%direction    (1:var_size),              &
                 reac_coord%const        (1:3,         1:var_size), &
                 reac_coord%num_atoms_bonds_reac_coord(1:var_size), &
                 reac_coord%indices      (1:var_size),              &
                 reac_coord%indices2     (1:var_size),              &
                 reac_coord%istarts      (1:var_size),              &
                 reac_coord%nnonzeros    (1:var_size),              &
                 stat = alloc_stat)

      reac_coord%kind         (1:var_size)              = 0
      reac_coord%grouplist    (1:var_size1, 1:var_size) = 0
      reac_coord%funcgrp      (1:var_size)              = 0
      reac_coord%exponent_dist(1:var_size3, 1:var_size) = 0
      reac_coord%unitv        (1:3,         1:var_size) = 0.0_wp
      reac_coord%unitv        (1:1,         1:var_size) = 1.0_wp
      reac_coord%weight_dist  (1:var_size3, 1:var_size) = 0.0_wp
      reac_coord%direction    (1:var_size)              = 1
      reac_coord%const        (1:3,         1:var_size) = 1.0_wp
      reac_coord%indices      (1:var_size)              = 0
      reac_coord%indices2     (1:var_size)              = 0
      reac_coord%istarts      (1:var_size)              = 0
      reac_coord%nnonzeros    (1:var_size)              = 0

    case (EneFuncReacatm)
      ! var_size = maxval(enefunc%reac_coord%num_atoms_bonds_reac_coord)
      !
      if (.not. allocated(reac_coord%bondslist_to_atomlist)) &
        allocate(reac_coord%bondslist_to_atomlist(var_size1, nreac), &
                 reac_coord%bonds_coord(3, var_size1),               &
                 reac_coord%bonds_force(3, var_size1),               &
                 stat = alloc_stat)

      reac_coord%bondslist_to_atomlist(1:var_size1, 1:nreac) = 0
      reac_coord%bonds_coord(1:3, 1:var_size1)               = 0.0_wp
      reac_coord%bonds_force(1:3, 1:var_size1)               = 0.0_wp

    case (EneFuncReacCoord)
      ! var_size  = reac%nreac_coord
      ! var_size1 = maxval(enefunc%reac_coord%num_atoms_bonds_reac_coord)
      !
      if(allocated(reac_coord%nxi)) then
        if(size(reac_coord%nxi) /= var_size) then
          deallocate(reac_coord%nxi,       &
                     reac_coord%ximin,     &
                     reac_coord%ximax,     &
                     reac_coord%dxi,       &
                     reac_coord%invdxi,    &
                     reac_coord%xi,        &
                     reac_coord%force_xi,  &
                     reac_coord%virial_xi, &
                     reac_coord%uboundary, &
                     reac_coord%periodicity, &
                     stat = dealloc_stat)
        endif
      endif

      if (.not. allocated(reac_coord%nxi))        &
        allocate(reac_coord%nxi      (var_size),  &
                 reac_coord%ximin    (var_size),  &
                 reac_coord%ximax    (var_size),  &
                 reac_coord%dxi      (var_size),  &
                 reac_coord%invdxi   (var_size),  &
                 reac_coord%xi       (var_size),  &
                 reac_coord%uboundary(var_size),  &
                 reac_coord%force_xi (var_size,   &
                 	                    3,          &
                                      var_size2), &  
                 reac_coord%virial_xi(var_size,   &
                                       3, 3),     &
                 reac_coord%periodicity(var_size), &
                 stat = alloc_stat)

    case default

      call error_msg('Alloc_Enefunc_Reac_Coord> Invalid variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc    

    return

  end subroutine alloc_enefunc_reac_coord

  function reac_coords_are_out_of_range(enefunc, remd) result(res)
    type(s_enefunc),  intent(inout), target   :: enefunc
    type(s_remd),     intent(in),    optional :: remd

    logical                             :: res
    type(s_enefunc_reac_coord), pointer :: reac_coord
    integer                             :: i, ierr
    reac_coord => enefunc%reac_coord
    res = .true.
    do i=1, reac_coord%nreac_coord
      if(reac_coord%xi(i) >= reac_coord%ximin(i) .and. &
         reac_coord%xi(i) <= reac_coord%ximax(i)) then
         res = .false.
         exit
      endif
    enddo

    if (present(remd)) then
      if(remd%total_nreplicas<1) then
        enefunc%reac_coord%out_of_range = res
        if(main_rank .and. res) then
          write(MsgOut,'(A)') 'All reaction coordinates are out of range'
        endif
        return
      endif
    else
      enefunc%reac_coord%out_of_range = res
      if(main_rank .and. res) then
        write(MsgOut,'(A)') 'All reaction coordinates are out of range'
      endif
      return
    end if  

!    if(replica_main_rank .and. res) then
!      write(MsgOut,'(A,i8,a)') 'All reaction coordinates for replica ', my_country_no+1, ' are out of range'
!    endif
#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, res, 1, mpi_logical, &
                       mpi_land, mpi_comm_airplane, ierr)
#endif
    enefunc%reac_coord%out_of_range = res
    if(main_rank .and. res) then
      write(MsgOut,'(A)') 'All reaction coordinates are out of range'
    endif

    return
  end function reac_coords_are_out_of_range

end module at_enefunc_reac_coord_mod
