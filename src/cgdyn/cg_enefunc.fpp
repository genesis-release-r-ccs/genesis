!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_enefunc_mod
!> @brief   define potential energy functions in each domain
!! @authors Jaewoon Jung (JJ), Yuji Sugita (YS), Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_enefunc_mod

  use cg_enefunc_gromacs_mod
  use cg_enefunc_localres_mod
  use cg_communicate_str_mod
  use cg_communicate_mod
  use cg_migration_mod
  use cg_energy_mod
  use cg_restraints_str_mod
  use cg_enefunc_str_mod
  use cg_energy_str_mod
  use cg_domain_str_mod
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
#ifdef OMP
    use omp_lib, only: omp_get_thread_num
#endif

  implicit none
  private

  ! subroutines
  public  :: define_enefunc
  public  :: define_enefunc_lb
  public  :: update_enefunc_aicg
  public  :: update_enefunc_martini
  public  :: copy_bond_information
  public  :: copy_angl_information
  public  :: copy_dihe_information
  public  :: copy_stack_information
  public  :: copy_contact_information
  public  :: copy_pwmcos_information
  public  :: copy_pwmcosns_information
  public  :: copy_rest_information
  private :: check_bonding

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
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc(ene_info, grotop, molecule, restraints, &
                            domain, enefunc, comm)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(inout) :: molecule
    type(s_restraints),      intent(in)    :: restraints
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_comm),            intent(inout) :: comm   


    enefunc%forcefield        = ene_info%forcefield
    enefunc%output_style      = ene_info%output_style
    enefunc%dielec_const      = ene_info%dielec_const
    enefunc%electrostatic     = ene_info%electrostatic
    enefunc%epsilon_rf        = ene_info%epsilon_rf  
    enefunc%assign_force_max  = ene_info%assign_force_max
    enefunc%upper_force_value = ene_info%upper_force_value

    ! gromacs
    !
    if (grotop%num_atomtypes > 0) then

      call define_enefunc_gromacs(ene_info, grotop, molecule, &
                                  restraints, domain, enefunc, comm)

    end if

    return

  end subroutine define_enefunc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc_lb
  !> @brief        update potential energy functions after load balance
  !! @authors      JJ
  !! @param[in]    grotop      : GROMACS parameter topology information
  !! @param[in]    molecule    : molecule information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : energy potential functions information
  !! @param[inout] comm        : communication information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc_lb(grotop, molecule, domain, enefunc, comm)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(inout) :: molecule
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_comm),            intent(inout) :: comm

    ! gromacs
    !
    if (grotop%num_atomtypes > 0) then

      call define_enefunc_gromacs_lb(grotop, molecule, domain, enefunc, comm)

    end if

    return

  end subroutine define_enefunc_lb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_enefunc_aicg
  !> @brief        a driver subroutine for updating potential energy functions
  !! @authors      JJ
  !! @param[in]    table       : flag for table or not
  !! @param[inout] domain      : domain information
  !! @param[inout] comm        : communication information
  !! @param[inout] enefunc     : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_enefunc_aicg(domain, comm, enefunc)

    ! formal arguments
    type(s_domain),                intent(inout) :: domain
    type(s_comm),                  intent(inout) :: comm
    type(s_enefunc),               intent(inout) :: enefunc

    ! local variables
    logical                        :: first

    ! sending the bonding information to other domain
    !

    ! bond
    !
    call timer(TimerMBond, TimerOn)
    call update_outgoing_enefunc_bondsq(domain, enefunc)
    call update_outgoing_enefunc_bond(domain, enefunc)
    call timer(TimerMBond, TimerOff)
    call timer(TimerComm3, TimerOn)
    call communicate_bond(domain, comm, enefunc)
    call timer(TimerComm3, TimerOff)
    call timer(TimerMBond, TimerOn)
    call update_incoming_enefunc_bond(domain, enefunc)
    call timer(TimerMBond, TimerOff)

    ! contact
    !
    call timer(TimerMContact, TimerOn)
    call update_enefunc_contact(domain, enefunc)
    call timer(TimerMContact, TimerOff)
    call timer(TimerComm3, TimerOn)
    call communicate_contact(domain, comm, enefunc)
    call timer(TimerComm3, TimerOff)

    ! angle
    !
    call timer(TimerMAngl, TimerOn)
    call update_outgoing_enefunc_anglflex(domain, enefunc)
    call update_outgoing_enefunc_angllocal(domain, enefunc)
    call update_outgoing_enefunc_angl(domain, enefunc)
    call timer(TimerMAngl, TimerOff)
    call timer(TimerComm3, TimerOn)
    call communicate_angl(domain, comm, enefunc)
    call timer(TimerComm3, TimerOff)
    call timer(TimerMAngl, TimerOn)
    call update_incoming_enefunc_angl(domain, enefunc)
    call timer(TimerMAngl, TimerOff)

    ! dihedral
    !
    call timer(TimerMDihe, TimerOn)
    call update_outgoing_enefunc_diheflex(domain, enefunc)
    call update_outgoing_enefunc_dihelocal(domain, enefunc)
    call update_outgoing_enefunc_dihe(domain, enefunc)
    call timer(TimerMDihe, TimerOff)
    call timer(TimerComm3, TimerOn)
    call communicate_dihe(domain, comm, enefunc)
    call timer(TimerComm3, TimerOff)
    call timer(TimerMDihe, TimerOn)
    call update_incoming_enefunc_dihe(domain, enefunc)
    call timer(TimerMDihe, TimerOff)


    ! restraint
    ! 
    if (enefunc%restraint) then
      call update_outgoing_enefunc_restraint(domain, enefunc)
      call timer(TimerComm3, TimerOn)
      call communicate_restraint(domain, comm, enefunc)
      call timer(TimerComm3, TimerOff)
      call update_incoming_enefunc_restraint(domain, enefunc)
    end if

    ! stack, pwmcos, pwmcosns
    !
    if (enefunc%forcefield == ForceFieldRESIDCG) then
      if (enefunc%num_base_stack_all > 0) then
        call timer(TimerMStack, TimerOn)
        call update_outgoing_enefunc_stack(domain, enefunc)
        call timer(TimerMStack, TimerOff)
        call timer(TimerComm3, TimerOn)
        call communicate_stack(domain, comm, enefunc)
        call timer(TimerComm3, TimerOff)
        call timer(TimerMStack, TimerOn)
        call update_incoming_enefunc_stack(domain, enefunc)
        call timer(TimerMStack, TimerOff)
      end if
      if (enefunc%cg_pwmcos_calc) then
        call timer(TimerMPWMcos, TimerOn)
        call update_outgoing_enefunc_pwmcos(domain, enefunc)
        call timer(TimerMPWMcos, TimerOff)
        call timer(TimerComm3, TimerOn)
        call communicate_pwmcos(domain, comm, enefunc)
        call timer(TimerComm3, TimerOff)
        call timer(TimerMPWMcos, TimerOn)
        call update_incoming_enefunc_pwmcos(domain, enefunc)
        call timer(TimerMPWMcos, TimerOff)
      end if
      if (enefunc%cg_pwmcosns_calc) then
        call timer(TimerMPWMcosns, TimerOn)
        call update_outgoing_enefunc_pwmcosns(domain, enefunc)
        call timer(TimerMPWMcosns, TimerOff)
        call timer(TimerComm3, TimerOn)
        call communicate_pwmcosns(domain, comm, enefunc)
        call timer(TimerComm3, TimerOff)
        call timer(TimerMPWMcosns, TimerOn)
        call update_incoming_enefunc_pwmcosns(domain, enefunc)
        call timer(TimerMPWMcosns, TimerOff)
      end if
    end if

    ! re-count nonbond exclusion list
    !
    first = .false.
!   call timer(TimerMNonb, TimerOn)
    call count_nonb_excl_go(domain, enefunc)
!   call timer(TimerMNonb, TimerOff)

    return

  end subroutine update_enefunc_aicg

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_enefunc_martini
  !> @brief        a driver subroutine for updating potential energy functions
  !! @authors      JJ
  !! @param[in]    table       : flag for table or not
  !! @param[inout] domain      : domain information
  !! @param[inout] comm        : communication information
  !! @param[inout] enefunc     : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_enefunc_martini(domain, comm, enefunc)

    ! formal arguments
    type(s_domain),                intent(inout) :: domain
    type(s_comm),                  intent(inout) :: comm
    type(s_enefunc),               intent(inout) :: enefunc

    ! local variables
    logical                        :: first

    ! sending the bonding information to other domain
    !

    ! bond
    !
    call timer(TimerMBond, TimerOn)
    call update_outgoing_enefunc_bondsq(domain, enefunc)
    call update_outgoing_enefunc_bond(domain, enefunc)
    call timer(TimerMBond, TimerOff)
    call timer(TimerComm3, TimerOn)
    call communicate_bond(domain, comm, enefunc)
    call timer(TimerComm3, TimerOff)
    call timer(TimerMBond, TimerOn)
    call update_incoming_enefunc_bond(domain, enefunc)
    call timer(TimerMBond, TimerOff)

    ! angle
    !
    call timer(TimerMAngl, TimerOn)
    call update_outgoing_enefunc_anglflex(domain, enefunc)
    call update_outgoing_enefunc_angllocal(domain, enefunc)
    call update_outgoing_enefunc_angl(domain, enefunc)
    call timer(TimerMAngl, TimerOff)
    call timer(TimerComm3, TimerOn)
    call communicate_angl(domain, comm, enefunc)
    call timer(TimerComm3, TimerOff)
    call timer(TimerMAngl, TimerOn)
    call update_incoming_enefunc_angl(domain, enefunc)
    call timer(TimerMAngl, TimerOff)

    ! dihedral
    !
    call timer(TimerMDihe, TimerOn)
    call update_outgoing_enefunc_diheflex(domain, enefunc)
    call update_outgoing_enefunc_dihelocal(domain, enefunc)
    call update_outgoing_enefunc_dihe(domain, enefunc)
    call timer(TimerMDihe, TimerOff)
    call timer(TimerComm3, TimerOn)
    call communicate_dihe(domain, comm, enefunc)
    call timer(TimerComm3, TimerOff)
    call timer(TimerMDihe, TimerOn)
    call update_incoming_enefunc_dihe(domain, enefunc)
    call timer(TimerMDihe, TimerOff)

    ! restraint
    ! 
    if (enefunc%restraint) then
      call update_outgoing_enefunc_restraint(domain, enefunc)
      call timer(TimerComm3, TimerOn)
      call communicate_restraint(domain, comm, enefunc)
      call timer(TimerComm3, TimerOff)
      call update_incoming_enefunc_restraint(domain, enefunc)
    end if

    ! re-count nonbond exclusion list
    !
    first = .false.
!   call timer(TimerMNonb, TimerOn)
    call count_nonb_excl(.false., domain, enefunc)
!   call timer(TimerMNonb, TimerOff)

    return

  end subroutine update_enefunc_martini

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
    integer                  :: i, j, ix, icel1, icel2, i1, i2, start_i
    integer                  :: icel3, i3,icel4, i4
    integer                  :: id, my_id
    real(wp), parameter      :: maxdistance = 0.5_wp
    real(wp)                 :: maxcell_size 

    real(wp),        pointer :: r0(:)
    integer,         pointer :: bondlist(:,:)
    integer,         pointer :: angllist(:,:)
    integer,         pointer :: dihelist(:,:)
    integer,         pointer :: ncell_local
    integer,         pointer :: id_g2l(:)
    real(wip),       pointer :: coord(:,:)

    ncell_local => domain%num_cell_local
    id_g2l      => domain%id_g2l
    coord       => domain%coord

    maxcell_size = max(domain%cell_size(1),  &
                       domain%cell_size(2),  &
                       domain%cell_size(3))

    bondlist    => enefunc%bond_list
    r0          => enefunc%bond_dist_min

    angllist    => enefunc%angl_list

    dihelist    => enefunc%dihe_list

    !$omp parallel default(shared)                                     &
    !$omp private(id, i, j, ix, icel1, i1, icel2, i2, d12, r12, r_dif, &
    !$omp         my_id, icel3, i3, icel4, i4, start_i)
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif
    my_id = id

    do i = my_id+1, enefunc%num_bondsq_domain+enefunc%num_bond_domain, nthread

      i1    = id_g2l(bondlist(1,i))
      i2    = id_g2l(bondlist(2,i))

      d12(1:3) = real(coord(i1,1:3),wp) - real(coord(i2,1:3),wp)
      r12   = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )
      r_dif = r12 - r0(i)
      if (r_dif > maxdistance) &
         write(MsgOut,'(A,I10,I10,F10.5)') &
        'WARNING: too long bond:',bondlist(1,i),bondlist(2,i),r12
      if (r_dif < -maxdistance) &
         write(MsgOut,'(A,I10,I10,F10.5)') &
        'WARNING: too short bond:',bondlist(1,i),bondlist(2,i),r12
      if (r12 > maxcell_size) then
         write(MsgOut,'(A,2I10,F10.5)') &
        'Check_bonding> distance is grater than cellsize:', &
         bondlist(1,i),bondlist(2,i),r12
         call error_msg('')
      endif

    end do

    do i = my_id+1, enefunc%num_angle_flexible_domain + &
                    enefunc%num_angle_local_domain    + &
                    enefunc%num_angle_domain, nthread

      i1    = id_g2l(angllist(1,i))
      i3    = id_g2l(angllist(3,i))

      d12(1:3) = real(coord(i1,1:3),wp) - real(coord(i3,1:3),wp)
      r12   = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )

      if (r12 > maxcell_size) then
         write(MsgOut,'(A,2I10,F10.5)') &
         'Check_bonding> distance in angle is grater than cellsize:', &
         angllist(1,i),angllist(3,i),r12
         call error_msg('')
      endif

    end do

    do i = my_id+1, enefunc%num_dihe_flexible_domain + &
                    enefunc%num_dihe_local_domain    + &
                    enefunc%num_dihe_domain, nthread

      i1    = id_g2l(dihelist(1,i))
      i4    = id_g2l(dihelist(4,i))

      d12(1:3) = real(coord(i1,1:3),wp) - real(coord(i4,1:3),wp)
      r12   = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )

      if (r12 > maxcell_size) then
         write(MsgOut,'(A,2I10,F10.5)') &
         'Check_bonding> distance in dihedral is grater than cellsize:', &
         dihelist(1,i),dihelist(4,i),r12
         call error_msg('')
      endif

    end do

    !$omp end parallel 

    return

  end subroutine check_bonding

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    copy_bond_information
  !> @brief        save/load bond data using temporary data 
  !! @authors      JJ
  !! @param[in]    direction   : 1 to temporary, 2 from temporary
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : enefunc information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine copy_bond_information(direction, domain, enefunc)

    ! formaal arguments
    integer,                 intent(in   )  :: direction
    type(s_domain),  target, intent(inout)  :: domain
    type(s_enefunc), target, intent(inout)  :: enefunc

    integer                  :: i
    real(wp),        pointer :: fc(:), r0(:)
    real(wip),       pointer :: buf_real(:)
    integer,         pointer :: bondlist(:,:), bondkind(:)
    integer,         pointer :: buf_int(:)

    bondlist      => enefunc%bond_list
    bondkind      => enefunc%bond_kind
    fc            => enefunc%bond_force_const
    r0            => enefunc%bond_dist_min

    buf_real      => domain%buf_var0_stay_real
    buf_int       => domain%buf_var0_stay_int

    if (direction == 1) then

      do i = 1, enefunc%num_bondsq_domain+enefunc%num_bond_domain
        buf_int (3*i-2) = bondlist(1,i)
        buf_int (3*i-1) = bondlist(2,i)
        buf_int (3*i  ) = bondkind(  i)
        buf_real(2*i-1) = fc      (  i)
        buf_real(2*i  ) = r0      (  i)
      end do

    else if (direction == 2) then

      do i = 1, enefunc%num_bondsq_domain+enefunc%num_bond_domain
        bondlist(1,i) = buf_int (3*i-2)
        bondlist(2,i) = buf_int (3*i-1)
        bondkind(  i) = buf_int (3*i  )
        fc      (  i) = buf_real(2*i-1)
        r0      (  i) = buf_real(2*i  )
      end do

    end if

    return

  end subroutine copy_bond_information

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    copy_angl_information
  !> @brief        save/load angle data using temporary data 
  !! @authors      JJ
  !! @param[in]    direction   : 1 to temporary, 2 from temporary
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : enefunc information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine copy_angl_information(direction, domain, enefunc)

    ! formaal arguments
    integer,                 intent(in   )  :: direction
    type(s_domain),  target, intent(inout)  :: domain
    type(s_enefunc), target, intent(inout)  :: enefunc

    integer                  :: i, num1, num2, num3
    integer                  :: nangl, nangl_flex, nangl_local
    real(wp),        pointer :: fc(:), r0(:), width(:), theta0(:)
    real(wp),        pointer :: fc_ub(:), r0_ub(:)
    real(wip),       pointer :: buf_real(:)
    integer,         pointer :: angllist(:,:), anglkind(:)
    integer,         pointer :: buf_int(:)

    angllist      => enefunc%angl_list
    anglkind      => enefunc%angl_kind
    fc            => enefunc%angl_force_const
    r0            => enefunc%angl_theta_min
    width         => enefunc%angl_width
    theta0        => enefunc%angl_theta_min
    fc_ub         => enefunc%urey_force_const
    r0_ub         => enefunc%urey_rmin

    buf_real      => domain%buf_var0_stay_real
    buf_int       => domain%buf_var0_stay_int

    nangl_flex    = enefunc%num_angle_flexible_domain
    nangl_local   = enefunc%num_angle_local_domain
    nangl         = enefunc%num_angle_domain

    if (direction == 1) then

      do i = 1, nangl_flex
        buf_int(4*i-3) = angllist(1,i)
        buf_int(4*i-2) = angllist(2,i)
        buf_int(4*i-1) = angllist(3,i)
        buf_int(4*i  ) = anglkind(  i)
      end do
      num1 = nangl_flex
      num2 = 4*num1
      do i = 1, nangl_local
        buf_int (num2+3*i-2) = angllist(1,num1+i)
        buf_int (num2+3*i-1) = angllist(2,num1+i)
        buf_int (num2+3*i  ) = angllist(3,num1+i)
        buf_real(     3*i-2) = fc      (  num1+i)
        buf_real(     3*i-1) = r0      (  num1+i)
        buf_real(     3*i  ) = width   (  num1+i)
      end do
      num1 = nangl_flex + nangl_local
      num2 = 4*nangl_flex + 3*nangl_local
      num3 = 3*nangl_local
      do i = 1, nangl
        buf_int (num2+4*i-3) = angllist(1,num1+i)
        buf_int (num2+4*i-2) = angllist(2,num1+i)
        buf_int (num2+4*i-1) = angllist(3,num1+i)
        buf_int (num2+4*i  ) = anglkind(  num1+i)
        buf_real(num3+4*i-3) = fc      (  num1+i)
        buf_real(num3+4*i-2) = theta0  (  num1+i)
        buf_real(num3+4*i-1) = fc_ub   (  num1+i)
        buf_real(num3+4*i  ) = r0_ub   (  num1+i)
      end do 

    else if (direction == 2) then

      do i = 1, nangl_flex
        angllist(1,i) = buf_int(4*i-3)
        angllist(2,i) = buf_int(4*i-2)
        angllist(3,i) = buf_int(4*i-1)
        anglkind(  i) = buf_int(4*i  )
      end do
      num1 = nangl_flex
      num2 = 4*num1
      do i = 1, nangl_local
        angllist(1,num1+i) = buf_int (num2+3*i-2)
        angllist(2,num1+i) = buf_int (num2+3*i-1)
        angllist(3,num1+i) = buf_int (num2+3*i  )
        fc      (  num1+i) = buf_real(     3*i-2)
        r0      (  num1+i) = buf_real(     3*i-1)
        width   (  num1+i) = buf_real(     3*i  )
      end do
      num1 = nangl_flex + nangl_local
      num2 = 4*nangl_flex + 3*nangl_local
      num3 = 3*nangl_local
      do i = 1, nangl
        angllist(1,num1+i) = buf_int (num2+4*i-3)
        angllist(2,num1+i) = buf_int (num2+4*i-2)
        angllist(3,num1+i) = buf_int (num2+4*i-1)
        anglkind(  num1+i) = buf_int (num2+4*i  )
        fc      (  num1+i) = buf_real(num3+4*i-3)
        theta0  (  num1+i) = buf_real(num3+4*i-2)
        fc_ub   (  num1+i) = buf_real(num3+4*i-1)
        r0_ub   (  num1+i) = buf_real(num3+4*i  )
      end do

    end if

    return

  end subroutine copy_angl_information

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    copy_dihe_information
  !> @brief        save/load dihedral data using temporary data 
  !! @authors      JJ
  !! @param[in]    direction   : 1 to temporary, 2 from temporary
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : enefunc information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine copy_dihe_information(direction, domain, enefunc)

    ! formaal arguments
    integer,                 intent(in   )  :: direction
    type(s_domain),  target, intent(inout)  :: domain
    type(s_enefunc), target, intent(inout)  :: enefunc

    integer                  :: i, num1, num2, num3
    integer                  :: ndihe, ndihe_flex, ndihe_local
    real(wp),        pointer :: fc(:), phase(:), width(:)
    real(wip),       pointer :: buf_real(:)
    integer,         pointer :: dihelist(:,:), dihekind(:), nperiod(:)
    integer,         pointer :: buf_int(:)

    dihelist      => enefunc%dihe_list
    dihekind      => enefunc%dihe_kind
    nperiod       => enefunc%dihe_periodicity
    fc            => enefunc%dihe_force_const
    phase         => enefunc%dihe_phase
    width         => enefunc%dihe_width

    buf_real      => domain%buf_var0_stay_real
    buf_int       => domain%buf_var0_stay_int

    ndihe_flex    = enefunc%num_dihe_flexible_domain
    ndihe_local   = enefunc%num_dihe_local_domain
    ndihe         = enefunc%num_dihe_domain

    if (direction == 1) then

      do i = 1, ndihe_flex
        buf_int(6*i-5) = dihelist(1,i)
        buf_int(6*i-4) = dihelist(2,i)
        buf_int(6*i-3) = dihelist(3,i)
        buf_int(6*i-2) = dihelist(4,i)
        buf_int(6*i-1) = dihekind(  i)
        buf_int(6*i  ) = nperiod (  i)
      end do
      num1 = ndihe_flex
      num2 = 6*num1
      do i = 1, ndihe_local
        buf_int (num2+5*i-4) = dihelist(1,num1+i)
        buf_int (num2+5*i-3) = dihelist(2,num1+i)
        buf_int (num2+5*i-2) = dihelist(3,num1+i)
        buf_int (num2+5*i-1) = dihelist(4,num1+i)
        buf_int (num2+5*i  ) = dihekind(  num1+i)
        buf_real(     3*i-2) = fc      (  num1+i)
        buf_real(     3*i-1) = phase   (  num1+i)
        buf_real(     3*i  ) = width   (  num1+i)
      end do
      num1 = ndihe_flex + ndihe_local
      num2 = 6*ndihe_flex + 5*ndihe_local
      num3 = 3*ndihe_local
      do i = 1, ndihe
        buf_int (num2+6*i-5) = dihelist(1,num1+i)
        buf_int (num2+6*i-4) = dihelist(2,num1+i)
        buf_int (num2+6*i-3) = dihelist(3,num1+i)
        buf_int (num2+6*i-2) = dihelist(4,num1+i)
        buf_int (num2+6*i-1) = nperiod (  num1+i)
        buf_int (num2+6*i  ) = dihekind(  num1+i)
        buf_real(num3+2*i-3) = fc      (  num1+i)
        buf_real(num3+2*i-2) = phase   (  num1+i)
      end do 

    else if (direction == 2) then

      do i = 1, ndihe_flex
        dihelist(1,i) = buf_int(6*i-5) 
        dihelist(2,i) = buf_int(6*i-4) 
        dihelist(3,i) = buf_int(6*i-3) 
        dihelist(4,i) = buf_int(6*i-2) 
        dihekind(  i) = buf_int(6*i-1) 
        nperiod (  i) = buf_int(6*i  ) 
      end do
      num1 = ndihe_flex
      num2 = 6*num1
      do i = 1, ndihe_local
        dihelist(1,num1+i) = buf_int (num2+5*i-4)
        dihelist(2,num1+i) = buf_int (num2+5*i-3)
        dihelist(3,num1+i) = buf_int (num2+5*i-2)
        dihelist(4,num1+i) = buf_int (num2+5*i-1)
        dihekind(  num1+i) = buf_int (num2+5*i  )
        fc      (  num1+i) = buf_real(     3*i-2)
        phase   (  num1+i) = buf_real(     3*i-1)
        width   (  num1+i) = buf_real(     3*i  )
      end do
      num1 = ndihe_flex + ndihe_local
      num2 = 6*ndihe_flex + 5*ndihe_local
      num3 = 3*ndihe_local
      do i = 1, ndihe
        dihelist(1,num1+i) = buf_int (num2+6*i-5) 
        dihelist(2,num1+i) = buf_int (num2+6*i-4) 
        dihelist(3,num1+i) = buf_int (num2+6*i-3) 
        dihelist(4,num1+i) = buf_int (num2+6*i-2) 
        nperiod (  num1+i) = buf_int (num2+6*i-1) 
        dihekind(  num1+i) = buf_int (num2+6*i  ) 
        fc      (  num1+i) = buf_real(num3+2*i-3) 
        phase   (  num1+i) = buf_real(num3+2*i-2) 
      end do

    end if

    return

  end subroutine copy_dihe_information

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    copy_stack_information
  !> @brief        save/load stack data using temporary data 
  !! @authors      JJ
  !! @param[in]    direction   : 1 to temporary, 2 from temporary
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : enefunc information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine copy_stack_information(direction, domain, enefunc)

    ! formaal arguments
    integer,                 intent(in   )  :: direction
    type(s_domain),  target, intent(inout)  :: domain
    type(s_enefunc), target, intent(inout)  :: enefunc

    integer                  :: i
    real(wp),        pointer :: epsilon(:), sigma(:), theta(:)
    real(wip),       pointer :: buf_real(:)
    integer,         pointer :: stacklist(:,:)
    integer,         pointer :: buf_int(:)

    stacklist      => enefunc%base_stack_list
    epsilon        => enefunc%base_stack_epsilon
    sigma          => enefunc%base_stack_sigma
    theta          => enefunc%base_stack_theta_bs

    buf_real      => domain%buf_var0_stay_real
    buf_int       => domain%buf_var0_stay_int

    if (direction == 1) then

      do i = 1, enefunc%num_stack_domain
        buf_int (3*i-2) = stacklist(1,i)
        buf_int (3*i-1) = stacklist(2,i)
        buf_int (3*i  ) = stacklist(3,i)
        buf_real(3*i-2) = epsilon  (  i)
        buf_real(3*i-1) = sigma    (  i)
        buf_real(3*i  ) = theta    (  i)
      end do

    else if (direction == 2) then

      do i = 1, enefunc%num_stack_domain
        stacklist(1,i) = buf_int (3*i-2) 
        stacklist(2,i) = buf_int (3*i-1) 
        stacklist(3,i) = buf_int (3*i  ) 
        epsilon  (  i) = buf_real(3*i-2) 
        sigma    (  i) = buf_real(3*i-1) 
        theta    (  i) = buf_real(3*i  ) 
      end do

    end if

    return

  end subroutine copy_stack_information

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    copy_contact_information
  !> @brief        save/load contact data using temporary data 
  !! @authors      JJ
  !! @param[in]    direction   : 1 to temporary, 2 from temporary
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : enefunc information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine copy_contact_information(direction, domain, enefunc)

    ! formaal arguments
    integer,                 intent(in   )  :: direction
    type(s_domain),  target, intent(inout)  :: domain
    type(s_enefunc), target, intent(inout)  :: enefunc

    integer                  :: i
    real(wp),        pointer :: lj12(:), lj10(:), lj6(:)
    real(wip),       pointer :: buf_real(:)
    integer,         pointer :: contactlist(:,:), contactfunc(:)
    integer,         pointer :: buf_int(:)

    contactlist    => enefunc%contact_list
    contactfunc    => enefunc%contact_func
    lj12           => enefunc%contact_lj12
    lj10           => enefunc%contact_lj10
    lj6            => enefunc%contact_lj6

    buf_real      => domain%buf_var0_stay_real
    buf_int       => domain%buf_var0_stay_int

    if (direction == 1) then

      do i = 1, enefunc%num_contact_domain + enefunc%num_contact_boundary
        buf_int (3*i-2) = contactlist(1,i)
        buf_int (3*i-1) = contactlist(2,i)
        buf_int (3*i  ) = contactfunc(  i)
        buf_real(3*i-2) = lj12       (  i)
        buf_real(3*i-1) = lj10       (  i)
        buf_real(3*i  ) = lj6        (  i)
      end do

    else if (direction == 2) then

      do i = 1, enefunc%num_contact_domain + enefunc%num_contact_boundary
        contactlist(1,i) = buf_int (3*i-2)
        contactlist(2,i) = buf_int (3*i-1)
        contactfunc(  i) = buf_int (3*i  )
        lj12       (  i) = buf_real(3*i-2)
        lj10       (  i) = buf_real(3*i-1)
        lj6        (  i) = buf_real(3*i  )
      end do

    end if

    return

  end subroutine copy_contact_information

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    copy_pwmcos_information
  !> @brief        save/load pwmcos data using temporary data 
  !! @authors      JJ
  !! @param[in]    direction   : 1 to temporary, 2 from temporary
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : enefunc information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine copy_pwmcos_information(direction, domain, enefunc)

    ! formaal arguments
    integer,                 intent(in   )  :: direction
    type(s_domain),  target, intent(inout)  :: domain
    type(s_enefunc), target, intent(inout)  :: enefunc

    integer                  :: i, i2, k
    real(wp),        pointer :: r0(:,:)
    real(wp),        pointer :: theta1(:,:), theta2(:,:), theta3(:,:)
    real(wp),        pointer :: ene_A(:,:), ene_C(:,:), ene_G(:,:)
    real(wp),        pointer :: ene_T(:,:), gamma(:,:), eps(:,:)
    real(wip),       pointer :: buf_real(:)
    integer,         pointer :: count(:)
    integer,         pointer :: id(:), id_N(:), id_C(:)
    integer,         pointer :: buf_int(:)

    count         => enefunc%pwmcos_count
    id            => enefunc%pwmcos_protein_id
    id_N          => enefunc%pwmcos_protein_id_N
    id_C          => enefunc%pwmcos_protein_id_C
    r0            => enefunc%pwmcos_r0
    theta1        => enefunc%pwmcos_theta1
    theta2        => enefunc%pwmcos_theta2
    theta3        => enefunc%pwmcos_theta3
    ene_A         => enefunc%pwmcos_ene_A
    ene_C         => enefunc%pwmcos_ene_C
    ene_G         => enefunc%pwmcos_ene_G
    ene_T         => enefunc%pwmcos_ene_T
    gamma         => enefunc%pwmcos_gamma
    eps           => enefunc%pwmcos_eps

    buf_real      => domain%buf_var0_stay_real
    buf_int       => domain%buf_var0_stay_int

    if (direction == 1) then

      do i = 1, enefunc%num_pwmcos_domain
        buf_int (4*i-3) = id   (i)
        buf_int (4*i-2) = id_N (i)
        buf_int (4*i-1) = id_C (i)
        buf_int (4*i  ) = count(i)
        do i2 = 1, count(i)
          k = 60*(i-1) + 10*(i2-1)
          buf_real(k+ 1) = r0    (i2,i) 
          buf_real(k+ 2) = theta1(i2,i) 
          buf_real(k+ 3) = theta2(i2,i) 
          buf_real(k+ 4) = theta3(i2,i) 
          buf_real(k+ 5) = ene_A (i2,i) 
          buf_real(k+ 6) = ene_C (i2,i) 
          buf_real(k+ 7) = ene_G (i2,i) 
          buf_real(k+ 8) = ene_T (i2,i) 
          buf_real(k+ 9) = gamma (i2,i) 
          buf_real(k+10) = eps   (i2,i) 
        end do
      end do

    else if (direction == 2) then

      do i = 1, enefunc%num_pwmcos_domain
        id   (i) = buf_int(4*i-3)
        id_N (i) = buf_int(4*i-2)
        id_C (i) = buf_int(4*i-1)
        count(i) = buf_int(4*i  )
        do i2 = 1, count(i)
          k = 60*(i-1) + 10*(i2-1)
          r0    (i2,i) = buf_real(k+ 1)
          theta1(i2,i) = buf_real(k+ 2)
          theta2(i2,i) = buf_real(k+ 3)
          theta3(i2,i) = buf_real(k+ 4)
          ene_A (i2,i) = buf_real(k+ 5)
          ene_C (i2,i) = buf_real(k+ 6)
          ene_G (i2,i) = buf_real(k+ 7)
          ene_T (i2,i) = buf_real(k+ 8)
          gamma (i2,i) = buf_real(k+ 9)
          eps   (i2,i) = buf_real(k+10)
        end do
      end do

    end if

    return

  end subroutine copy_pwmcos_information

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    copy_pwmcosns_information
  !> @brief        save/load pwmcosns data using temporary data 
  !! @authors      JJ
  !! @param[in]    direction   : 1 to temporary, 2 from temporary
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : enefunc information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine copy_pwmcosns_information(direction, domain, enefunc)

    ! formaal arguments
    integer,                 intent(in   )  :: direction
    type(s_domain),  target, intent(inout)  :: domain
    type(s_enefunc), target, intent(inout)  :: enefunc

    integer                  :: i, i2, k1, k2
    real(wp),        pointer :: r0(:,:)
    real(wp),        pointer :: theta1(:,:), theta2(:,:)
    real(wp),        pointer :: ene(:,:)
    real(wip),       pointer :: buf_real(:)
    integer,         pointer :: specificity(:,:)
    integer,         pointer :: count(:)
    integer,         pointer :: id(:), id_N(:), id_C(:)
    integer,         pointer :: buf_int(:)

    count         => enefunc%pwmcosns_count
    id            => enefunc%pwmcosns_protein_id
    id_N          => enefunc%pwmcosns_protein_id_N
    id_C          => enefunc%pwmcosns_protein_id_C
    r0            => enefunc%pwmcosns_r0
    theta1        => enefunc%pwmcosns_theta1
    theta2        => enefunc%pwmcosns_theta2
    ene           => enefunc%pwmcosns_ene
    specificity   => enefunc%pwmcosns_specificity

    buf_real      => domain%buf_var0_stay_real
    buf_int       => domain%buf_var0_stay_int

    if (direction == 1) then

      do i = 1, enefunc%num_pwmcosns_domain
        buf_int (10*i-9) = id   (i)
        buf_int (10*i-8) = id_N (i)
        buf_int (10*i-7) = id_C (i)
        buf_int (10*i-6) = count(i)
        do i2 = 1, count(i)
          k1 = 24*(i-1) + 4*(i2-1)
          k2 = 10*i - 6 + i2
          buf_real(k2  ) = specificity(i2,i) 
          buf_real(k1+1) = r0         (i2,i) 
          buf_real(k1+2) = theta1     (i2,i) 
          buf_real(k1+3) = theta2     (i2,i) 
          buf_real(k1+4) = ene        (i2,i) 
        end do
      end do

    else if (direction == 2) then

      do i = 1, enefunc%num_pwmcosns_domain
        id   (i) = buf_int(10*i-9)
        id_N (i) = buf_int(10*i-8)
        id_C (i) = buf_int(10*i-7)
        count(i) = buf_int(10*i-6)
        do i2 = 1, count(i)
          k1 = 24*(i-1) + 4*(i2-1)
          k2 = 10*i - 6 + i2
          specificity(i2,i) = buf_real(k2  )
          r0         (i2,i) = buf_real(k1+1)
          theta1     (i2,i) = buf_real(k1+2)
          theta2     (i2,i) = buf_real(k1+3)
          ene        (i2,i) = buf_real(k1+4)
        end do
      end do

    end if

    return

  end subroutine copy_pwmcosns_information

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    copy_rest_information
  !> @brief        save/load rest data using temporary data 
  !! @authors      JJ
  !! @param[in]    direction   : 1 to temporary, 2 from temporary
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : enefunc information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine copy_rest_information(direction, domain, enefunc)

    ! formaal arguments
    integer,                 intent(in   )  :: direction
    type(s_domain),  target, intent(inout)  :: domain
    type(s_enefunc), target, intent(inout)  :: enefunc

    integer                  :: i
    real(wp),        pointer :: restraint_coord(:,:), restraint_force(:,:)
    real(wip),       pointer :: buf_real(:)
    integer,         pointer :: restraint_atom(:)
    integer,         pointer :: buf_int(:)

    restraint_atom   => enefunc%restraint_atom
    restraint_coord  => enefunc%restraint_coord
    restraint_force  => enefunc%restraint_force

    buf_real         => domain%buf_var0_stay_real
    buf_int          => domain%buf_var0_stay_int

    if (direction == 1) then

      do i = 1, enefunc%num_rest_domain
        buf_int (    i) = restraint_atom (  i)
        buf_real(7*i-6) = restraint_coord(1,i)
        buf_real(7*i-5) = restraint_coord(2,i)
        buf_real(7*i-4) = restraint_coord(3,i)
        buf_real(7*i-3) = restraint_force(1,i)
        buf_real(7*i-2) = restraint_force(2,i)
        buf_real(7*i-1) = restraint_force(3,i)
        buf_real(7*i  ) = restraint_force(3,i)
      end do

    else if (direction == 2) then

      do i = 1, enefunc%num_rest_domain 
        restraint_atom (  i) = buf_int (    i) 
        restraint_coord(1,i) = buf_real(7*i-6)
        restraint_coord(2,i) = buf_real(7*i-5)
        restraint_coord(3,i) = buf_real(7*i-4)
        restraint_force(1,i) = buf_real(7*i-3)
        restraint_force(2,i) = buf_real(7*i-2)
        restraint_force(3,i) = buf_real(7*i-1)
        restraint_force(4,i) = buf_real(7*i  )
      end do

    end if

    return

  end subroutine copy_rest_information

end module cg_enefunc_mod
