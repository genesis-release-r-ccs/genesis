!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_energy_reac_coord_mod
!> @brief   calculate energy and force from reaction coordinates
!! @authors Junichiro Koga (JK)
!
!  (c) Copyright 2020 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_energy_reac_coord_mod
  use at_reac_coord_str_mod
  use at_enefunc_str_mod
  use at_restraints_str_mod
  use at_experiments_mod
  use at_boundary_str_mod
  use dihedral_libs_mod
  use mpi_parallel_mod
  use timers_mod
  use messages_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: compute_energy_reac_coord
  public  :: compute_reac_coord
  private :: compute_energy_reac_coord_dist
  private :: get_reac_coord_coordinates
  private :: get_reac_coord_forces

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_reac_coord
  !> @brief        calculate reaction coordinate and it is forces
  !! @authors      JK
  !! @param[in]    domain     : domain information
  !! @param[inout] enefunc    : potential energy functions information
  !! @param[in]    coord      : coordinates of target systems
  !! @param[inout] force      : forces of target systems
  !! @param[inout] virial     : virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_reac_coord(get_coord, calc_force, boundary,&
                                       enefunc, coord, force, virial)

    ! formal arguments
    logical,                 intent(in)    :: get_coord
    logical,                 intent(in)    :: calc_force
    type(s_boundary),        intent(in)    :: boundary
    type(s_enefunc), target, intent(inout) :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(dp),                intent(inout) :: virial(3,3)

    ! local variables
    real(dp)                 :: val, eangle, edihed, etilt, eexp, erepul, efb
    real(wp),        pointer :: bonds_coord(:,:), bonds_force(:,:)
    integer,         pointer :: bonds_to_atom(:)
    integer                  :: i, j, inum, jnum, dockid
    real(dp),    allocatable :: viri(:,:)

    allocate(viri(3,3))
    bonds_coord   => enefunc%reac_coord%bonds_coord
    bonds_force   => enefunc%reac_coord%bonds_force

    if (enefunc%bias_use) dockid = my_country_no + 1

    do i = 1, enefunc%abmd%nreac_coord
      inum = i
      if (enefunc%bias_use) then
        jnum = enefunc%abmd%dockid2repsetid(dockid)
        inum = enefunc%abmd%ifunc(enefunc%abmd%repsetid2paramid(jnum,1), i)
      end if
      bonds_to_atom => enefunc%abmd%bondslist_to_atomlist(:, i)

      if (get_coord) then
        call get_reac_coord_coordinates(inum, enefunc, coord, &
                                        bonds_to_atom, bonds_coord)
      end if
      bonds_force(1:3,1:enefunc%reac_coord%num_atoms_bonds_reac_coord(inum)) = 0.0_wp

      val    = 0.0_dp
      viri   = 0.0_dp

      select case (enefunc%reac_coord%kind(inum))

      case (ReacCoordsFuncPOSI, ReacCoordsFuncPOSICOM)
        call compute_energy_reac_coord_posi(inum, calc_force, enefunc,     &
                               bonds_coord, val, bonds_force, viri)
      case (ReacCoordsFuncDIST, ReacCoordsFuncDISTCOM)
        call compute_energy_reac_coord_dist(inum, calc_force, enefunc,     &
                               bonds_coord, val, bonds_force, viri)
      case (ReacCoordsFuncANGLE, ReacCoordsFuncANGLECOM)
        call compute_energy_reac_coord_angle(inum, calc_force, enefunc,    &
                               bonds_coord, val, bonds_force, viri)
      case (ReacCoordsFuncDIHED, ReacCoordsFuncDIHEDCOM)
        call compute_energy_reac_coord_dihed(inum, calc_force, enefunc,    &
                               bonds_coord, val, bonds_force, viri)

      end select
      enefunc%reac_coord%xi(inum) = val

      do j = 1, enefunc%reac_coord%num_atoms_bonds_reac_coord(inum)
        enefunc%abmd%force_xi(i,1:3,bonds_to_atom(j)) = 0.0_wp
      end do

      if (calc_force) then
        call get_reac_coord_forces(inum, enefunc, bonds_to_atom,   &
             bonds_force, enefunc%abmd%force_xi(i,:,:))
        if(main_rank .or. replica_main_rank) then
          enefunc%abmd%virial_xi(i,:,:) = viri(:,:)
        endif
      end if
    enddo

    return

  end subroutine compute_energy_reac_coord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_reac_coord
  !> @brief        calculate reaction coordinate
  !! @authors      SI
  !! @param[inout] enefunc    : potential energy functions information
  !! @param[in]    coord      : coordinates of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_reac_coord(enefunc, coord)

    ! formal arguments
    type(s_enefunc), target, intent(inout) :: enefunc
    real(wp),                intent(in)    :: coord(:,:)

    ! local variables
    integer                  :: i
    real(dp)                 :: val
    real(wp),        pointer :: bonds_coord(:,:)
    integer,         pointer :: bonds_to_atom(:)

    bonds_coord   => enefunc%reac_coord%bonds_coord

    do i = 1, enefunc%reac_coord%nreac_coord
      bonds_to_atom => enefunc%reac_coord%bondslist_to_atomlist(:, i)

      call get_reac_coord_coordinates(i, enefunc, coord, &
                                      bonds_to_atom, bonds_coord)
      
      val    = 0.0_dp

      select case (enefunc%reac_coord%kind(i))

      case (ReacCoordsFuncPOSI, ReacCoordsFuncPOSICOM)
        call compute_energy_reac_coord_posi(i, .false., enefunc,     &
                               bonds_coord, val)
      case (ReacCoordsFuncDIST, ReacCoordsFuncDISTCOM)
        call compute_energy_reac_coord_dist(i, .false., enefunc,     &
                               bonds_coord, val)
      case (ReacCoordsFuncANGLE, ReacCoordsFuncANGLECOM)
        call compute_energy_reac_coord_angle(i, .false., enefunc,    &
                               bonds_coord, val)
      case (ReacCoordsFuncDIHED, ReacCoordsFuncDIHEDCOM)
        call compute_energy_reac_coord_dihed(i, .false., enefunc,    &
                               bonds_coord, val)

      end select
      enefunc%reac_coord%xi(i) = val
    enddo

    return

  end subroutine compute_reac_coord
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_reac_coord_posi
  !> @brief        calculate position reaction coordinate energy
  !! @authors      JK
  !! @param[in]    inum       : ID for the target reaction coordinate
  !! @param[in]    enefunc    : potential energy functions information
  !! @param[in]    pos_coord  : reaction coordinate bonds coordinates
  !! @param[inout] pos_force  : reaction coordinate bonds forces
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] dist       : the target distance
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_reac_coord_posi(inum,calc_force, enefunc, &
                                            pos_coord, pos, pos_force, virial)

    ! formal arguments
    integer,                 intent(in)    :: inum
    logical,                 intent(in)    :: calc_force
    type(s_enefunc), target, intent(inout) :: enefunc
    real(wp),                intent(in)    :: pos_coord(:,:)
    real(wp),      optional, intent(inout) :: pos_force(:,:)
    real(dp),      optional, intent(inout) :: virial(:,:)
    real(dp),                intent(inout) :: pos

    ! local variable
    real(wp)                 :: r_sum, coef, coefi, grad_coef
    real(wp)                 :: tmp_wp
    integer                  :: numgrp
    integer                  :: i, j, k, jx, index, id, id1, id2, ind
    integer                  :: iatm
    integer                  :: max_atoms
    integer                  :: i_rpath_dim
    integer                  :: alloc_stat, dealloc_stat

    real(wp)                 :: wtmp
    real(wp)                 :: d(3)
    real(wp),        pointer :: const(:,:)
    real(wp),        pointer :: ref(:,:)
    real(wp),        pointer :: wtind(:,:)
    real(wp),        pointer :: masscoef(:,:)
    real(wp),        pointer :: com(:,:), dij(:,:), rij(:), coefgrp(:,:), unitv(:)
    integer,         pointer :: num_funcs, direction, numatoms(:)
    integer,         pointer :: bondslist(:,:)
    integer,         pointer :: kind(:)
    integer,         pointer :: grplist(:,:)
    integer,         pointer :: funcgrp(:)

    num_funcs  => enefunc%reac_coord%nreac_coord
    funcgrp    => enefunc%reac_coord%funcgrp
    grplist    => enefunc%reac_coord%grouplist
    numatoms   => enefunc%reac_coord%numatoms
    bondslist  => enefunc%reac_coord%bondslist(:,:,inum)
    wtind      => enefunc%reac_coord%weight_dist
    kind       => enefunc%reac_coord%kind
    masscoef   => enefunc%reac_coord%masscoef
    direction  => enefunc%reac_coord%direction(inum)
    unitv      => enefunc%reac_coord%unitv(:,inum)

    ! initialization
    !

    numgrp = funcgrp(inum)

    max_atoms      = enefunc%reac_coord%max_atoms
    allocate(coefgrp(max_atoms, numgrp), &
             com(3, numgrp),             &
             stat = alloc_stat)

    if (kind(inum) == ReacCoordsFuncPOSI) then
      do j = 1, numgrp
        coefgrp(1:numatoms(grplist(j,inum)),j) = &
                1.0_wp/real(numatoms(grplist(j,inum)),wp)
      end do
    else
      do j = 1, numgrp
        coefgrp(1:numatoms(grplist(j,inum)),j) = &
                masscoef(1:numatoms(grplist(j,inum)),grplist(j,inum))
      end do
    end if

!    d = 1.0_wp
!    if(direction == 2) then !X only
!      d(2) = 0.0_wp
!      d(3) = 0.0_wp
!    else if (direction == 3) then !Y only
!      d(1) = 0.0_wp
!      d(3) = 0.0_wp
!    else if (direction == 4) then ! Z only
!      d(1) = 0.0_wp
!      d(2) = 0.0_wp
!    endif
    d(:) = unitv(:)

    ! check atom group range and calculate center of mass of each group
    !
    do j = 1, numgrp
      com(1:3,j) = 0.0_wp
      do jx = 1, numatoms(grplist(j,inum))
        index = bondslist(jx,grplist(j,inum))
        com(1:3,j) = com(1:3,j) + pos_coord(1:3,index) * coefgrp(jx,j) * d(1:3)
      end do
    end do

    ! calculation of posi energy
    !
    r_sum = 0.0_wp
    do j = 1, numgrp
      r_sum = r_sum + sum(com(1:3,j))
    enddo

    pos = r_sum

    ! calculation of force
    !
    if (calc_force) then
      do j = 1, numgrp
        do jx = 1, numatoms(grplist(j,inum))
          index = bondslist(jx,grplist(j,inum))
          pos_force(1:3, index) = pos_force(1:3, index) - coefgrp(jx, j) * d(1:3)
        end do
      enddo

      ! virial
      !
      if(main_rank .or. replica_main_rank) then
        do j = 1, numgrp
          do jx = 1, numatoms(grplist(j,inum))
            index = bondslist(jx,grplist(j,inum))
            do i = 1, 3
              do k = j+1, 3
                wtmp = -pos_coord(k,index)*coefgrp(jx, j) * d(k) * d(i)
                virial(k,i) = virial(k,i) + wtmp
                virial(i,k) = virial(i,k) + wtmp
              end do
              wtmp = -pos_coord(j,index)*coefgrp(jx, j)
              virial(i,i) = virial(i,i) + wtmp
            end do
          enddo
        enddo
      endif
    end if

    deallocate(coefgrp, com, stat = dealloc_stat)

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc    

    return

  end subroutine compute_energy_reac_coord_posi

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_reac_coord_dist
  !> @brief        calculate distance reaction coordinate energy
  !! @authors      JK
  !! @param[in]    inum       : ID for the target reaction coordinate
  !! @param[in]    enefunc    : potential energy functions information
  !! @param[in]    dist_coord : reaction coordinate bonds coordinates
  !! @param[inout] dist_force : reaction coordinate bonds forces
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] dist       : the target distance
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_reac_coord_dist(inum,calc_force, enefunc, &
                                            dist_coord, dist, dist_force, virial)

    ! formal arguments
    integer,                 intent(in)    :: inum
    logical,                 intent(in)    :: calc_force
    type(s_enefunc), target, intent(inout) :: enefunc
    real(wp),                intent(in)    :: dist_coord(:,:)
    real(wp),      optional, intent(inout) :: dist_force(:,:)
    real(dp),      optional, intent(inout) :: virial(:,:)
    real(dp),                intent(inout) :: dist

    ! local variable
    real(wp)                 :: sig(1:2)
    real(wp)                 :: work(1:3)
    real(wp)                 :: r_sum, dtmp, rtmp, r_dif, coef, coefi, grad_coef
    real(wp)                 :: tmp_wp
    integer                  :: numgrp, iexpo, idist, ndist
    integer                  :: i, j, jx, index, id, id1, id2, ind
    integer                  :: iatm
    integer                  :: max_atoms
    integer                  :: i_rpath_dim
    integer                  :: alloc_stat, dealloc_stat

    real(wp),        pointer :: const(:,:)
    real(wp),        pointer :: ref(:,:)
    real(wp),        pointer :: wtind(:,:)
    real(wp),        pointer :: masscoef(:,:)
    real(wp),        pointer :: com(:,:), dij(:,:), rij(:), coefgrp(:,:)
    integer,         pointer :: num_funcs, numatoms(:)
    integer,         pointer :: bondslist(:,:)
    integer,         pointer :: direction(:)
    integer,         pointer :: kind(:)
    integer,         pointer :: grplist(:,:)
    integer,         pointer :: funcgrp(:)
    integer,         pointer :: expoind(:,:)
    
    num_funcs  => enefunc%reac_coord%nreac_coord
    funcgrp    => enefunc%reac_coord%funcgrp
    grplist    => enefunc%reac_coord%grouplist
    numatoms   => enefunc%reac_coord%numatoms
    bondslist  => enefunc%reac_coord%bondslist(:,:,inum)
    expoind    => enefunc%reac_coord%exponent_dist
    wtind      => enefunc%reac_coord%weight_dist
    kind       => enefunc%reac_coord%kind
    masscoef   => enefunc%reac_coord%masscoef
    const      => enefunc%reac_coord%const
    direction  => enefunc%reac_coord%direction

    ! initialization
    !
    numgrp = funcgrp(inum)
    ndist  = int(numgrp/2)

    max_atoms      = enefunc%reac_coord%max_atoms
    allocate(coefgrp(max_atoms, numgrp), &
             com(3, numgrp),             &
             dij(3, ndist),              &
             rij(ndist),                 &
             stat = alloc_stat)

    if (kind(inum) == ReacCoordsFuncDIST) then
      do j = 1, numgrp
        coefgrp(1:numatoms(grplist(j,inum)),j) = &
                1.0_wp/real(numatoms(grplist(j,inum)),wp)
      end do
    else
      do j = 1, numgrp
        coefgrp(1:numatoms(grplist(j,inum)),j) = &
                masscoef(1:numatoms(grplist(j,inum)),grplist(j,inum))
      end do
    end if

    ! check atom group range and calculate center of mass of each group
    !
    do j = 1, numgrp
      com(1:3,j) = 0.0_wp
      do jx = 1, numatoms(grplist(j,inum))
        index = bondslist(jx,grplist(j,inum))
        com(1:3,j) = com(1:3,j) + dist_coord(1:3,index) * coefgrp(jx,j)
      end do
    end do

    ! calculation of distance energy
    !
    r_sum = 0.0_wp

    do idist = 1, ndist

      id1 = idist*2 - 1
      id2 = idist*2
      dtmp = 0.0_wp
      dij(1:3,idist) = com(1:3,id1) - com(1:3,id2)
      dij(1:3,idist) = dij(1:3,idist) * const(1:3, inum)

      if (direction(inum) == 5 .or. direction(inum) == 6 .or. direction(inum) == 7) then
        rij(idist)     = dij(direction(inum)-4,idist)

      else
        rij(idist)     = sqrt(dij(1,idist)*dij(1,idist) +  &
                              dij(2,idist)*dij(2,idist) +  &
                              dij(3,idist)*dij(3,idist))
      end if 

      if (expoind(idist,inum) == 1) then
        r_sum = r_sum + wtind(idist,inum)*rij(idist)
      else
        r_sum = r_sum + wtind(idist,inum)*rij(idist)**(expoind(idist,inum))
      end if

    end do

    dist = r_sum

    sig(1) = -1.0_wp
    sig(2) = 1.0_wp
    
    ! calculation of force
    !
    if (calc_force) then
      do idist = 1,ndist

        if (rij(idist) < 1e-10_wp) then
          coefi = 0.0_wp
        else
!          coefi = wtind(idist,inum) / rij(idist)
           if (expoind(idist,inum) == 1) then
             coefi =  wtind(idist,inum) / rij(idist)
           else
             coefi = expoind(idist,inum) &
                    * rij(idist)**(expoind(idist,inum)-2)
           end if
        end if

        do i = 1, 2
          id = idist*2-2+i
          do j = 1,numatoms(grplist(id,inum))
            iatm = bondslist(j,grplist(id,inum))
            work(1:3) = sig(i) * coefi  * coefgrp(j,id) * dij(1:3,idist)
            dist_force(1:3,iatm) = dist_force(1:3,iatm) + work(1:3)
          end do
        end do

        ! calculation of virial
        !
        if (main_rank .or. replica_main_rank) then
          do i = 1, 3
            work(i)= coefi * dij(i,idist)
            do j = 1, 3
              virial(j,i) = virial(j,i) - work(i) * dij(j,idist)
            end do
          end do
        end if
      end do
    end if
    
    deallocate(coefgrp, com, dij, rij, stat = dealloc_stat)

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc    

    return

  end subroutine compute_energy_reac_coord_dist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_reac_coord_angle
  !> @brief        calculate angle restraint energy
  !! @authors      JK
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    angle_coord : reaction coordinate bonds coordinates
  !! @param[inout] angle_force : reaction coordinate bonds forces
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] angle       : value of the angle
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_reac_coord_angle(inum, calc_force, enefunc, angle_coord, &
                                             angle, angle_force, virial)

    ! formal arguments
    integer,                 intent(in)    :: inum
    logical,                 intent(in)    :: calc_force
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: angle_coord(:,:)
    real(wp),      optional, intent(inout) :: angle_force(:,:)
    real(dp),      optional, intent(inout) :: virial(:,:)
    real(dp),                intent(inout) :: angle

    ! local variable
    real(wp)                 :: dij(1:3), dkj(1:3)
    real(wp)                 :: inv_rij2, inv_rkj2
    real(wp)                 :: rij2, rkj2, rijrkj,inv_rijrkj, vtmp
    real(wp)                 :: cos_t, t_dif, coef, theta, sin_t, grad_coef
    real(wp)                 :: work(1:6)
    real(wp)                 :: radref(1:2)
    integer                  :: numgrp, iexpo
    integer                  :: j, k, jx, index, ind
    integer                  :: iatm
    logical                  :: anglecheck

    real(wp),        pointer :: const(:,:)
    real(wp),        pointer :: ref(:,:)
    real(wp),        pointer :: masscoef(:,:)
    real(wp),        pointer :: com(:,:), coefgrp(:,:)
    integer,         pointer :: num_funcs, numatoms(:)
    integer,         pointer :: bondslist(:,:)
    integer,         pointer :: kind(:)
    integer,         pointer :: grplist(:,:)
    integer,         pointer :: funcgrp(:)
    integer,         pointer :: expo(:)
    integer                  :: max_atoms
    integer                  :: alloc_stat, dealloc_stat


    funcgrp    => enefunc%reac_coord%funcgrp
    grplist    => enefunc%reac_coord%grouplist
    numatoms   => enefunc%reac_coord%numatoms
    bondslist  => enefunc%reac_coord%bondslist(:,:,inum)
    kind       => enefunc%reac_coord%kind
    masscoef   => enefunc%reac_coord%masscoef

    ! initialization
    !
    numgrp = funcgrp(inum)

    max_atoms      = enefunc%reac_coord%max_atoms
    allocate(coefgrp(max_atoms, numgrp), &
             com(3, numgrp), stat = alloc_stat)

    if (kind(inum) == ReacCoordsFuncANGLE) then
      do j = 1, numgrp
        coefgrp(1:numatoms(grplist(j,inum)),j) = &
                1.0_wp/real(numatoms(grplist(j,inum)),wp)
      end do
    else
      do j = 1, numgrp
        coefgrp(1:numatoms(grplist(j,inum)),j) = &
                masscoef(1:numatoms(grplist(j,inum)),grplist(j,inum))
      end do
    end if

    ! check atom group range and calculate center of mass of each group
    !
    do j = 1, numgrp
      com(1:3,j) = 0.0_wp
      do jx = 1, numatoms(grplist(j,inum))
        index = bondslist(jx,grplist(j,inum))
        com(1:3,j) = com(1:3,j) + angle_coord(1:3,index) * coefgrp(jx,j)
      end do
    end do

    ! calculation of angle energy
    !
    dij(1:3) = com(1:3,1) - com(1:3,2)
    dkj(1:3) = com(1:3,3) - com(1:3,2)
    rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
    rkj2  = dkj(1)*dkj(1) + dkj(2)*dkj(2) + dkj(3)*dkj(3)
    rijrkj = sqrt(rij2*rkj2)

    inv_rijrkj = 1.0_wp /rijrkj
    inv_rij2 = 1.0_wp /rij2
    inv_rkj2 = 1.0_wp /rkj2

    cos_t = (dij(1)*dkj(1) + dij(2)*dkj(2) + dij(3)*dkj(3) )*inv_rijrkj
    cos_t  = min(  1.0_wp, cos_t )
    cos_t  = max( -1.0_wp, cos_t )
    theta = acos(cos_t)
    anglecheck = .true.
    do while(anglecheck)
      if (theta > pi)  theta = theta - 2.0_wp * pi
      if (theta < -pi) theta = theta + 2.0_wp * pi
      if (theta >= -pi  .and. theta <= pi)  anglecheck = .false.
    end do

    angle      = 180.0_wp * theta/pi
    grad_coef  = 180.0_wp/pi
    
    sin_t = (1.0_wp-cos_t*cos_t)
    sin_t = sqrt(sin_t)
    sin_t  = max( EPS, sin_t )
    grad_coef = -grad_coef / sin_t

    work(1:3) = grad_coef * ( dkj(1:3) * inv_rijrkj -  &
                              dij(1:3) * cos_t * inv_rij2 )

    work(4:6) = grad_coef * ( dij(1:3) * inv_rijrkj -  &
                              dkj(1:3) * cos_t * inv_rkj2 )
 
    ! compute force
    !
    if (calc_force)  then    
      do j = 1,numatoms(grplist(1,inum))
        iatm = bondslist(j,grplist(1,inum))
        angle_force(1:3,iatm) = angle_force(1:3,iatm)                     &
                               - work(1:3)*coefgrp(j,1)
      end do
 
      do j = 1,numatoms(grplist(2,inum))
        iatm = bondslist(j,grplist(2,inum))
        angle_force(1:3,iatm) = angle_force(1:3,iatm) +                   &
                              (work(1:3) + work(4:6))*coefgrp(j,2)
      end do
 
      do j = 1,numatoms(grplist(3,inum))
        iatm = bondslist(j,grplist(3,inum))
        angle_force(1:3,iatm) = angle_force(1:3,iatm)                     &
                               - work(4:6)*coefgrp(j,3)
      end do

      ! calculation of virial
      !
      if (main_rank .or. replica_main_rank) then
        do j = 1, 3
     
          do k = j+1, 3
            vtmp = -(dij(k)*work(j) + dkj(k)*work(j+3))
            virial(k,j) = virial(k,j) + vtmp
            virial(j,k) = virial(j,k) + vtmp
          end do
     
          vtmp = -(dij(j)*work(j) + dkj(j)*work(j+3))
          virial(j,j) = virial(j,j) + vtmp
        end do
      end if
    end if

    deallocate(coefgrp, com, stat = dealloc_stat)

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc    

    return

  end subroutine compute_energy_reac_coord_angle


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_reac_coord_dihed
  !> @brief        calculate dihedral restraint energy
  !! @authors      JK
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    dihed_coord : reaction coordinate bonds coordinates
  !! @param[inout] dihed_force : reaction coordinate bonds forces
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] dihed       : dihedral angle
  !! @note         Blondel and Karplus, J. Comput. Chem., 17, 1132-1141 (1996)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_reac_coord_dihed(inum, calc_force, enefunc, dihed_coord, &
                                             dihed, dihed_force, virial)

    ! formal arguments
    integer,                 intent(in)    :: inum
    logical,                 intent(in)    :: calc_force
    type(s_enefunc), target, intent(inout) :: enefunc
    real(wp),                intent(in)    :: dihed_coord(:,:)
    real(wp),      optional, intent(inout) :: dihed_force(:,:)
    real(dp),      optional, intent(inout) :: virial(:,:)
    real(dp),                intent(inout) :: dihed

    ! local variables
    real(wp)                 :: cos_dih, sin_dih, coef
    real(wp)                 :: cospha, sinpha, grad_coef, diffphi
    real(wp)                 :: cosdif, sindif, theta, vtmp
    real(wp)                 :: work(1:9)
    real(wp)                 :: radref(1:2)
    real(wp)                 :: tmp_wp
    integer                  :: aindex(1:4)
    integer                  :: numgrp, iexpo
    integer                  :: i, j, k, iatm, ind
    integer                  :: max_atoms
    integer                  :: alloc_stat, dealloc_stat
    integer                  :: i_rpath_dim
    real(wp)                 :: grad(1:9), v(1:3,1:3)
    real(wp)                 :: erestdihed

    real(wp),        pointer :: const(:,:)
    real(wp),        pointer :: ref(:,:)
    real(wp),        pointer :: masscoef(:,:)
    real(wp),        pointer :: com(:,:)
    real(wp),        pointer :: coefgrp(:,:)
    integer,         pointer :: num_funcs, numatoms(:)
    integer,         pointer :: bondslist(:,:)
    integer,         pointer :: kind(:)
    integer,         pointer :: grplist(:,:)
    integer,         pointer :: funcgrp(:)

    funcgrp    => enefunc%reac_coord%funcgrp 
    grplist    => enefunc%reac_coord%grouplist
    numatoms   => enefunc%reac_coord%numatoms
    bondslist  => enefunc%reac_coord%bondslist(:,:,inum)
    kind       => enefunc%reac_coord%kind
    masscoef   => enefunc%reac_coord%masscoef

    numgrp  = funcgrp(inum)
    max_atoms      = enefunc%reac_coord%max_atoms

    allocate(coefgrp(max_atoms, numgrp), &
             com(3, numgrp), stat = alloc_stat)

    if (kind(inum) == ReacCoordsFuncDIHED) then
      do j = 1, numgrp
        coefgrp(1:numatoms(grplist(j,inum)),j) = &
            1.0_wp/real(numatoms(grplist(j,inum)),wp)
      end do
    else
      do j = 1, numgrp
        coefgrp(1:numatoms(grplist(j,inum)),j) = &
                masscoef(1:numatoms(grplist(j,inum)),grplist(j,inum))
      end do
    end if
 
    ! calculation of com
    !
    
    com(1:3,1:4)=0.0_wp
 
    do i = 1,4 
      do j = 1,numatoms(grplist(i,inum))
        iatm = bondslist(j,grplist(i,inum))
        com(1:3,i) = com(1:3,i) + dihed_coord(1:3,iatm) * coefgrp(j,i)
      end do
      aindex(i) = i
    end do

    call calculate_dihedral(aindex, com, cos_dih, sin_dih, grad, v)


    if (cos_dih > 1.0E-1_wp) then
      theta = asin(sin_dih)
    else
      cos_dih = max(cos_dih, -1.0_wp)
      theta = sign(1.0_wp,sin_dih)*acos(cos_dih)
    endif

    cospha = cos(0.0_wp)
    sinpha = sin(0.0_wp)
      
    cosdif =  cos_dih
    sindif = -sin_dih
      
    if (cosdif > 1.0E-1_wp) then
      diffphi = asin(sindif)
    else
      cosdif  = max(cosdif, -1.0_wp)
      diffphi = sign(1.0_wp,sindif)*acos(cosdif)
    endif
    dihed     = -diffphi*180.0_wp/pi
    grad_coef = -180.0_wp/pi

    if (calc_force) then       
      ! compute virial
      !
      if (main_rank .or. replica_main_rank) then
        do j = 1, 3
          do k = j+1, 3
            vtmp = grad_coef*v(k,j)
         
            virial(k,j) = virial(k,j) + vtmp
            virial(j,k) = virial(j,k) + vtmp
          end do
     
          vtmp = grad_coef*v(j,j)
          virial(j,j) = virial(j,j) + vtmp
        end do
      endif
      if (abs(theta) < 1e-10_wp) then
        grad_coef = 0.0_wp
      end if
      work(1:9) = grad_coef*grad(1:9)

      dihed_force = 0.0_wp
      ! compute force
      !
      do j = 1,numatoms(grplist(1,inum))
        iatm = bondslist(j,grplist(1,inum))
        dihed_force(1:3,iatm) = dihed_force(1:3,iatm)             &
                               - work(1:3)*coefgrp(j,1)
      end do
   
      do j = 1,numatoms(grplist(2,inum))
        iatm = bondslist(j,grplist(2,inum))
        dihed_force(1:3,iatm) = dihed_force(1:3,iatm) +           &
                              (work(1:3) - work(4:6))*coefgrp(j,2)
      end do
   
      do j = 1,numatoms(grplist(3,inum))
        iatm = bondslist(j,grplist(3,inum))
        dihed_force(1:3,iatm) = dihed_force(1:3,iatm) +           &
                              (work(4:6) + work(7:9))*coefgrp(j,3)
      end do
   
      do j = 1,numatoms(grplist(4,inum))
        iatm = bondslist(j,grplist(4,inum))
        dihed_force(1:3,iatm) = dihed_force(1:3,iatm)             &
                               - work(7:9)*coefgrp(j,4)
      end do
    end if

    deallocate(coefgrp, com, stat = dealloc_stat)

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc 
    
    return

  end subroutine compute_energy_reac_coord_dihed

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_reac_coord_coordinates
  !> @brief        get coordinates for reaction coordinates calculation
  !! @authors      JK, SI
  !! @param[in]    enefunc       : potential energy functions information
  !! @param[in]    coord         : coordinates of target systems
  !! @param[in]    bonds_to_atom : bonds list to atom list for restrains
  !! @param[inout] bonds_coord   : reaction coordinate bonds coordinates
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine get_reac_coord_coordinates(inum, enefunc, coord, bonds_to_atom, &
                                        bonds_coord)

    ! formal arguments
    integer,                 intent(in)    :: inum
    type(s_enefunc),         intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    integer,                 intent(in)    :: bonds_to_atom(:)
    real(wp),                intent(inout) :: bonds_coord(:,:)

    ! local variable
    integer                  :: i


    bonds_coord(1:3,1:enefunc%reac_coord%num_atoms_bonds_reac_coord(inum)) = 0.0_wp

    do i = 1, enefunc%reac_coord%num_atoms_bonds_reac_coord(inum)
      bonds_coord(1:3,i) = coord(1:3,bonds_to_atom(i))
    end do

#ifdef HAVE_MPI_GENESIS
    !call mpi_allreduce(mpi_in_place, bonds_coord, &
    !                   3*enefunc%reac_coord%num_atoms_bonds_reac_coord(inum), &
    !                   mpi_wp_real, mpi_sum, mpi_comm_country, ierror)
#endif

    return

  end subroutine get_reac_coord_coordinates

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_reac_coord_forces
  !> @brief        get reaction coordinate forces
  !! @authors      JK, SI
  !! @param[in]    enefunc       : potential energy functions information
  !! @param[in]    dist_to_atom  : coordinates of target systems
  !! @param[in]    dist_force    : bonds list to atom list for restrains
  !! @param[inout] force         : reaction coordinate bonds coordinates
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine get_reac_coord_forces(inum, enefunc, dist_to_atom, &
                                   dist_force, force)

    ! formal arguments
    integer,                 intent(in)    :: inum
    type(s_enefunc),         intent(in)    :: enefunc
    integer,                 intent(in)    :: dist_to_atom(:)
    real(wp),                intent(inout) :: dist_force(:,:)
    real(wp),                intent(inout) :: force(:,:)

    ! local variable
    integer                  :: i


    do i = 1, enefunc%reac_coord%num_atoms_bonds_reac_coord(inum)
      force(1:3,dist_to_atom(i)) = force(1:3,dist_to_atom(i)) + dist_force(1:3,i)
    end do

    return

  end subroutine get_reac_coord_forces

end module at_energy_reac_coord_mod

