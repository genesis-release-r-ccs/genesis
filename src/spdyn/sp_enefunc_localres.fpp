!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_enefunc_localres_mod
!> @brief   restraint energy functions
!! @authors Chigusa Kobayashi (CK), Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8
  
#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_enefunc_localres_mod

  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use fileio_localres_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef MPI
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: setup_enefunc_localres
  private :: setup_enefunc_localres_bond
  private :: setup_enefunc_localres_angle
  private :: setup_enefunc_localres_dihed

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_localres
  !> @brief        define local restraint for each cell in 
  !                potential energy function
  !! @authors      CK
  !! @param[in]    localres : local restraint information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_localres(localres, domain, enefunc)

    ! formal arguments
    type(s_localres), target, intent(in)    :: localres
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(inout) :: enefunc


    ! bond
    !
    call setup_enefunc_localres_bond(localres, domain, enefunc)

    ! angle
    !
    call setup_enefunc_localres_angle(localres, domain, enefunc)

    ! dihedral
    !
    call setup_enefunc_localres_dihed(localres, domain, enefunc)


    return

  end subroutine setup_enefunc_localres

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_localres_bond
  !> @brief        define local restraint for each cell in 
  !                potential energy function
  !! @authors      CK
  !! @param[in]    localres : local restraint information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_localres_bond(localres, domain, enefunc)

    ! formal arguments
    type(s_localres), target, intent(in)    :: localres
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variable
    integer                   :: i, icel_local, i1, i2
    integer                   :: irest, icel1, icel2, found, oldbonds

    real(wp),         pointer :: const(:),ref(:)
    real(wp),         pointer :: force(:,:), dist(:,:)
    integer,          pointer :: ncell
    integer(int2),    pointer :: id_g2l(:,:)
    integer,          pointer :: num_funcs
    integer,          pointer :: index_atoms(:,:)
    integer,          pointer :: func(:)
    integer,          pointer :: bond(:), list(:,:,:)
    integer(int2),    pointer :: cell_pair(:,:)
    integer(1),       pointer :: bkind(:,:)


    num_funcs   => localres%num_funcs
    func        => localres%func
    index_atoms => localres%index_atoms
    const       => localres%const
    ref         => localres%ref

    ncell       => domain%num_cell_local
    cell_pair   => domain%cell_pair
    id_g2l      => domain%id_g2l

    bond        => enefunc%num_bond
    list        => enefunc%bond_list
    force       => enefunc%bond_force_const
    dist        => enefunc%bond_dist_min
    bkind       => enefunc%bond_kind


    if (num_funcs > 0) then
      enefunc%local_restraint = .true.
    else
      return
    end if

    ! locan bond restraint
    !
    irest = 0
    do i = 1, num_funcs

      if (localres%func(i) == LocalResBonds) then

        i1 = index_atoms(1,i)
        i2 = index_atoms(2,i)

        icel1 = domain%id_g2l(1,i1)
        icel2 = domain%id_g2l(1,i2)
        irest = irest+1

        if (icel1 /= 0 .and. icel2 /= 0) then
          icel_local = cell_pair(icel1,icel2)

          if (icel_local > 0 .and. icel_local <= ncell) then
            bond (icel_local) = bond(icel_local) + 1
            list (1:2,bond(icel_local),icel_local) = index_atoms(1:2,i)
            force(bond(icel_local),icel_local)     = const(i)
            dist (bond(icel_local),icel_local)     = ref(i)
            bkind(bond(icel_local),icel_local)     = 1
          end if
        end if

      end if

    end do

    found = 0
    oldbonds = enefunc%num_bond_all

    do i = 1,ncell
      found = found + bond(i)
      if (bond(i) > MaxBond) &
        call error_msg('Setup_Enefunc_Localres_Bond> Too many bonds.') 
    end do

#ifdef MPI
    call mpi_allreduce(found, enefunc%num_bond_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_bond_all = found
#endif

    if (enefunc%num_bond_all-oldbonds /= irest .and.                           &
        enefunc%num_bond_all+2*enefunc%table%num_water-oldbonds /= irest .and. &
        enefunc%num_bond_all+3*enefunc%table%num_water-oldbonds /= irest)      &
      call error_msg( &
           'Setup_Enefunc_Localres_Bond> Some Bond Paremeters are missing.')

    return

  end subroutine setup_enefunc_localres_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_localres_angle
  !> @brief        define local restraint for each cell in 
  !                potential energy function
  !! @authors      CK
  !! @param[in]    localres : local restraint information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_localres_angle(localres, domain, enefunc)

    ! formal arguments
    type(s_localres), target, intent(in)    :: localres
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variable
    integer                   :: i, icel_local
    integer                   :: i1, i2, i3
    integer                   :: irest, icel1, icel2, found, oldangle

    real(wp),         pointer :: const(:),ref(:)
    real(wp),         pointer :: ubforce(:,:), ubrmin(:,:)
    real(wp),         pointer :: force(:,:), theta(:,:)
    integer,          pointer :: ncell
    integer(int2),    pointer :: id_g2l(:,:)
    integer,          pointer :: num_funcs
    integer,          pointer :: index_atoms(:,:)
    integer,          pointer :: func(:)
    integer,          pointer :: angle(:), list(:,:,:)
    integer(int2),    pointer :: cell_pair(:,:)
    integer(1),       pointer :: akind(:,:)


    num_funcs   => localres%num_funcs
    func        => localres%func
    index_atoms => localres%index_atoms
    const       => localres%const
    ref         => localres%ref

    ncell       => domain%num_cell_local
    cell_pair   => domain%cell_pair
    id_g2l      => domain%id_g2l

    angle       => enefunc%num_angle
    list        => enefunc%angle_list
    force       => enefunc%angle_force_const
    theta       => enefunc%angle_theta_min
    ubforce     => enefunc%urey_force_const
    ubrmin      => enefunc%urey_rmin
    akind       => enefunc%angle_kind

    if (num_funcs > 0) then
      enefunc%local_restraint = .true.
    else
      return
    end if

    ! position restraint
    !
    irest = 0
    do i = 1, num_funcs

      if (localres%func(i) == LocalResAngles) then

        i1 = index_atoms(1,i)
        i2 = index_atoms(2,i)
        i3 = index_atoms(3,i)

        icel1 = domain%id_g2l(1,i1)
        icel2 = domain%id_g2l(1,i3)
        irest = irest+1

        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1,icel2)

          if (icel_local > 0 .and. icel_local <= ncell) then
            angle  (icel_local) = angle(icel_local) + 1
            list   (1:3,angle(icel_local),icel_local) = index_atoms(1:3,i)
            force  (angle(icel_local),icel_local)     = const(i)
            theta  (angle(icel_local),icel_local)     = ref(i)*RAD
            ubforce(angle(icel_local),icel_local)     = 0.0_wp
            ubrmin (angle(icel_local),icel_local)     = 0.0_wp
            akind  (angle(icel_local),icel_local)     = 1
          end if
        end if

      end if

    end do

    found = 0
    oldangle = enefunc%num_angl_all

    do i = 1,ncell
      found = found + angle(i)
      if (angle(i) > MaxAngle) &
        call error_msg('Setup_Enefunc_Localres_Angle> Too many angles.') 
    end do

#ifdef MPI
    call mpi_allreduce(found, enefunc%num_angl_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_angl_all = found
#endif

    if (enefunc%num_angl_all-oldangle /= irest .and.                    &
        enefunc%num_angl_all-oldangle+enefunc%table%num_water /= irest) &
      call error_msg( &
           'Setup_Enefunc_Localres_Angle> Some Angle Paremeters are missing.')

    return

  end subroutine setup_enefunc_localres_angle

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_localres_dihed
  !> @brief        define local restraint for each cell in 
  !                potential energy function
  !! @authors      CK
  !! @param[in]    localres : local restraint information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_localres_dihed(localres, domain, enefunc)

    ! formal arguments
    type(s_localres), target, intent(in)    :: localres
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variable
    integer                   :: i, icel_local 
    integer                   :: i1, i2, i3, i4
    integer                   :: irest, icel1, icel2, found, olddihed

    real(wp),         pointer :: const(:),ref(:)
    real(wp),         pointer :: force(:,:), theta(:,:)
    integer,          pointer :: ncell
    integer(int2),    pointer :: id_g2l(:,:)
    integer,          pointer :: num_funcs
    integer,          pointer :: index_atoms(:,:)
    integer,          pointer :: func(:)
    integer,          pointer :: dihed(:), list(:,:,:)
    integer(int2),    pointer :: cell_pair(:,:)
    integer(1),       pointer :: dkind(:,:)


    num_funcs   => localres%num_funcs
    func        => localres%func
    index_atoms => localres%index_atoms
    const       => localres%const
    ref         => localres%ref

    ncell       => domain%num_cell_local
    cell_pair   => domain%cell_pair
    id_g2l      => domain%id_g2l

    dihed       => enefunc%num_dihedral
    list        => enefunc%dihe_list
    force       => enefunc%dihe_force_const
    theta       => enefunc%dihe_phase
    dkind       => enefunc%dihe_kind

    if (num_funcs > 0) then
      enefunc%local_restraint = .true.
    else
      return
    end if

    ! position restraint
    !
    irest = 0
    do i = 1, num_funcs

      if (localres%func(i) == LocalResDihedrals) then

        i1 = index_atoms(1,i)
        i2 = index_atoms(2,i)
        i3 = index_atoms(3,i)
        i4 = index_atoms(4,i)

        icel1 = domain%id_g2l(1,i1)
        icel2 = domain%id_g2l(1,i4)
        irest = irest+1

        if (icel1 /= 0 .and. icel2 /= 0) then
          icel_local = cell_pair(icel1,icel2)

          if (icel_local > 0 .and. icel_local <= ncell) then
            dihed(icel_local) = dihed(icel_local) + 1
            list (1:4,dihed(icel_local),icel_local) = index_atoms(1:4,i)
            force(dihed(icel_local),icel_local)     = const(i)
            theta(dihed(icel_local),icel_local)     = ref(i)*RAD
            dkind(dihed(icel_local),icel_local)     = 1
          end if
        end if

      end if

    end do

    found = 0
    olddihed = enefunc%num_dihe_all

    do i = 1,ncell
      found = found + dihed(i)
      if (dihed(i) > maxdihe) &
        call error_msg( &
             'Setup_Enefunc_Localres_Dihed> Too many dihedral angles.') 
    end do

#ifdef MPI
    call mpi_allreduce(found, enefunc%num_dihe_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_dihe_all = found
#endif

    if (enefunc%num_dihe_all-olddihed /= irest) &
      call error_msg( &
           'Setup_Enefunc_Localres_Dihed> Some Dihed Paremeters are missing.')

    return

  end subroutine setup_enefunc_localres_dihed

end module sp_enefunc_localres_mod
