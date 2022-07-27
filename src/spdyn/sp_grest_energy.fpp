!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_grest_energy_mod
!> @brief   compute energy with grest parameters
!! @authors Jaewoon Jung(JJ), Yuji Sugita (YS), Takaharu Mori (TM),
!!          Chigusa Kobayashi (CK), Norio Takase (NT)
!  
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_grest_energy_mod

  use sp_energy_pme_mod
  use sp_energy_mod          
  use sp_boundary_str_mod
  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_energy_str_mod
  use sp_ensemble_str_mod
  use sp_remd_str_mod
  use sp_ensemble_str_mod
  use sp_domain_str_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use math_libs_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: compute_grest_energy_output

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_grest_energy_output
  !> @brief        compute energy output of each replica with all parameters
  !! @authors      JJ
  !! @param[in]    int_scheme  : integrator scheme
  !! @param[in]    alloc_ref   : flag for allocate reference or not
  !! @param[in]    dealloc_ref : flag for deallocate reference or not
  !! @param[in]    boundary    : boundary information
  !! @param[in]    pairlist    : pairlist information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : potential energy function information
  !! @param[inout] remd        : remd information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_grest_energy_output(int_scheme, alloc_ref, dealloc_ref,   &
                                         boundary, pairlist, ensemble, domain, &
                                         enefunc, remd)

    ! formal arguments
    integer,                 intent(in)    :: int_scheme
    logical,                 intent(in)    :: alloc_ref, dealloc_ref
    type(s_boundary),        intent(in)    :: boundary
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_ensemble),        intent(in)    :: ensemble 
    type(s_domain),  target, intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_remd),            intent(inout) :: remd

    ! local variables
    real(wip),       pointer :: coord(:,:,:), coord_ref(:,:,:), coord_old(:,:,:)
    integer,         pointer :: natom(:)
    integer                  :: ncell
    integer                  :: i, ix, j, k, ifunc
    integer                  :: replica_id, param_id
    real(wp)                 :: rest_param


    coord     => domain%coord
    coord_ref => domain%coord_ref
    coord_old => domain%coord_old
    natom     => domain%num_atom
    ncell     =  domain%num_cell_local+domain%num_cell_boundary

    if (int_scheme == 1) then

      ! we use coord_ref as coordinate
      !
      do i = 1, ncell
        do ix = 1, natom(i)
          coord_old(1:3,ix,i) = coord    (1:3,ix,i)
          coord    (1:3,ix,i) = coord_ref(1:3,ix,i)
        end do
      end do

    end if

    ! allocate reference force constant
    !
    if (alloc_ref) &
      call memory_alloc_grest_ref(domain, enefunc, remd)
    
    ! save the current force field to the reference
    !
    call copy_grest_to_ref(domain, enefunc, remd)
    do i = 1, remd%dimension
      if (remd%types(i) == RemdSoluteTempering) then
        remd%solute_tempering(i)%rest_param_full_ref = &
                      remd%solute_tempering(i)%rest_param_full
        remd%solute_tempering(i)%rest_param_half_ref = &
                      remd%solute_tempering(i)%rest_param_half
      end if
    end do

    ! energy calculation by changing grest parameters
    ! 
    replica_id = my_country_no + 1
    param_id   = remd%repid2parmsetid(replica_id)

    do i = 1, remd%dimension

      if (remd%types(i) == RemdSoluteTempering) then

        call compute_energy(domain, enefunc, pairlist, boundary, coord, &
                            .false., .true., .true., .true., .false.,   &
                            remd%grest_energy%energy,                   &
                            domain%atmcls_pbc,                          &
                            domain%translated,                          &
                            remd%grest_energy%force,                    &
                            remd%grest_energy%force_long,               &
                            domain%force_omp,                           &
                            domain%force_pbc,                           &
                            domain%virial_cellpair,                     &
                            remd%grest_energy%virial,                   &
                            remd%grest_energy%virial_long,              &
                            remd%grest_energy%virial_extern)

        remd%grest_energy%potential_energy(remd%parmidsets(param_id,i)) = &
                                           remd%grest_energy%energy%total

        do k = 1, remd%nreplicas(i)
 
          if (k == remd%parmidsets(param_id,i)) cycle
          rest_param = remd%dparameters(i,k)
          call assign_condition_solute_tempering (remd%solute_tempering(i), &
                                                  rest_param, domain,       &
                                                  ensemble, enefunc)

          ! energy calculation from the changed parameter
          !
          call compute_energy(domain, enefunc, pairlist, boundary, coord, &
                              .false., .true., .true., .true., .false.,   &
                              remd%grest_energy%energy,                   &
                              domain%atmcls_pbc,                          &
                              domain%translated,                          &
                              remd%grest_energy%force,                    &
                              remd%grest_energy%force_long,               &
                              domain%force_omp,                           &
                              domain%force_pbc,                           &
                              domain%virial_cellpair,                     &
                              remd%grest_energy%virial,                   &
                              remd%grest_energy%virial_long,              &
                              remd%grest_energy%virial_extern)             

          remd%grest_energy%potential_energy(k) = &
                           remd%grest_energy%energy%total
 
        end do

      end if
    end do 

    ! parameters are back
    !
    call copy_grest_from_ref(domain, remd, enefunc)
    do i = 1, remd%dimension
      if (remd%types(i) == RemdSoluteTempering) then
        remd%solute_tempering(i)%rest_param_full = &
                      remd%solute_tempering(i)%rest_param_full_ref
        remd%solute_tempering(i)%rest_param_half = &
                      remd%solute_tempering(i)%rest_param_half_ref
      end if
    end do

    if (int_scheme == 1) then

      ! original coordinate value
      !
      do i = 1, ncell
        do ix = 1, natom(i)
          coord(1:3,ix,i) = coord_old(1:3,ix,i) 
        end do
      end do

    end if

    if (dealloc_ref) &
      call memory_dealloc_grest_ref(domain, enefunc, remd)

    return

  end subroutine compute_grest_energy_output

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    memory_alloc_grest_ref
  !> @brief        allocation of arrays for saving the current force constant
  !! @authors      JJ
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy function information
  !! @param[inout] remd        : remd information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine memory_alloc_grest_ref(domain, enefunc, remd)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    type(s_remd),    target, intent(inout) :: remd

    ! local variables
    integer                          :: alloc_stat
    integer                          :: i, ncell_local, ncell
    integer                          :: ncmap_type, ncls
    type(s_soltemp),         pointer :: soltemp(:)
 

    soltemp     => remd%solute_tempering
 
    ncell_local = domain%num_cell_local
    ncell       = ncell_local + domain%num_cell_boundary

    do i = 1, remd%dimension

      if (remd%types(i) == RemdSoluteTempering) then
    
        ! basic array for coordinate and force
        !
        if (.not. allocated(remd%grest_energy%force) ) then
          allocate(remd%grest_energy%force      (3,MaxAtom,ncell),   &   
                   remd%grest_energy%force_long (3,MaxAtom,ncell),   &   
                   stat = alloc_stat)
          if (alloc_stat /= 0)   call error_msg_alloc
        end if

        ! allocation for charge array
        ! 
        if (soltemp(i)%sw_charge) then
          if (.not. allocated(remd%grest_energy%charge) ) then
            allocate(remd%grest_energy%charge(MaxAtom, ncell), &
                     stat = alloc_stat)
            if (alloc_stat /= 0)   call error_msg_alloc
          end if
        end if 

        ! allocation for bond force constant
        !
        if (soltemp(i)%sw_bonds) then
          if (.not. allocated(remd%grest_energy%bond_force_const) ) then
            allocate(remd%grest_energy%bond_force_const(MaxBond, ncell_local), &
                     stat = alloc_stat)
            if (alloc_stat /= 0)   call error_msg_alloc
          end if
        end if 

        ! allocation for angle/urey force constant
        !
        if (soltemp(i)%sw_angles) then
          if (.not. allocated(remd%grest_energy%angle_force_const) ) then
            allocate(remd%grest_energy%angle_force_const(MaxAngle, ncell_local), &
                     stat = alloc_stat)
            if (alloc_stat /= 0)   call error_msg_alloc
          end if
        end if

        if (soltemp(i)%sw_ureys) then
          if (.not. allocated(remd%grest_energy%urey_force_const) ) then
            allocate(remd%grest_energy%urey_force_const(MaxAngle, ncell_local), &
                     stat = alloc_stat)
            if (alloc_stat /= 0)   call error_msg_alloc
          end if
        end if

        ! allocation for dihedral force constant
        !
        if (soltemp(i)%sw_dihedrals) then
          if (.not. allocated(remd%grest_energy%dihe_force_const) ) then
            allocate(remd%grest_energy%dihe_force_const(MaxDihe, ncell_local), &
                     stat = alloc_stat)
            if (alloc_stat /= 0)   call error_msg_alloc
          end if
        end if

        ! allocation for rb_dihedral
        !
        if (soltemp(i)%sw_rb_dihedrals) then
          if (.not. allocated(remd%grest_energy%rb_dihe_c) ) then
            allocate(remd%grest_energy%rb_dihe_c(6, MaxDihe, ncell_local), &
                     stat = alloc_stat)
            if (alloc_stat /= 0)   call error_msg_alloc
          end if
        end if

        ! allocatation for improper dihedral force constant
        !
        if (soltemp(i)%sw_impropers) then
          if (.not. allocated(remd%grest_energy%impr_force_const) ) then
            allocate(remd%grest_energy%impr_force_const(MaxImpr, ncell_local), &
                     stat = alloc_stat)
            if (alloc_stat /= 0)   call error_msg_alloc
          end if
        end if

        ! allocation for cmap coefficients
        !
        if (soltemp(i)%num_cmap_type > 0 .and. soltemp(i)%sw_cmaps) then
          if (.not. allocated(remd%grest_energy%cmap_coef) ) then
            ncmap_type = size(enefunc%cmap_coef(1,1,1,1,:))
            allocate(remd%grest_energy%cmap_coef(4,4,24,24,ncmap_type), &
                     stat = alloc_stat)
            if (alloc_stat /= 0)   call error_msg_alloc
          end if
        end if

        ! allocation for lj
        !
        if (soltemp(i)%sw_lj) then
          ncls = enefunc%num_atom_cls
          if (.not. allocated(remd%grest_energy%nb14_lj6) .and. &
              allocated(enefunc%nb14_lj6) ) then
            allocate(remd%grest_energy%nb14_lj6 (ncls,ncls), &
                     stat = alloc_stat)
            if (alloc_stat /= 0)   call error_msg_alloc
          end if
          if (.not. allocated(remd%grest_energy%nb14_lj12) .and. &
              allocated(enefunc%nb14_lj12) ) then
            allocate(remd%grest_energy%nb14_lj12(ncls,ncls), &
                     stat = alloc_stat)
            if (alloc_stat /= 0)   call error_msg_alloc
          end if
          if (.not. allocated(remd%grest_energy%nonb_lj6) .and. &
              allocated(enefunc%nonb_lj6) ) then
            allocate(remd%grest_energy%nonb_lj6 (ncls,ncls), &
                     stat = alloc_stat)
            if (alloc_stat /= 0)   call error_msg_alloc
          end if
          if (.not. allocated(remd%grest_energy%nonb_lj12) .and. &
              allocated(enefunc%nonb_lj12) ) then
            allocate(remd%grest_energy%nonb_lj12(ncls,ncls), &
                     stat = alloc_stat)
            if (alloc_stat /= 0)   call error_msg_alloc
          end if
        end if
     
      end if
    end do

    return

  end subroutine memory_alloc_grest_ref
     
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    memory_dealloc_grest_ref
  !> @brief        deallocation of arrays for saving the current force constant
  !! @authors      JJ
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy function information
  !! @param[inout] remd        : remd information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine memory_dealloc_grest_ref(domain, enefunc, remd)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    type(s_remd),    target, intent(inout) :: remd

    ! local variables
    integer                          :: dealloc_stat
    integer                          :: i, ncell_local, ncell
    integer                          :: ncmap_type, ncls


    if (allocated(remd%grest_energy%force)) then
      deallocate(remd%grest_energy%force,            &
                 remd%grest_energy%force_long,       &
                 stat = dealloc_stat)
      if (dealloc_stat /= 0)   call error_msg_dealloc
    end if

    if ( allocated(remd%grest_energy%charge) ) then
      deallocate(remd%grest_energy%charge, stat = dealloc_stat)
      if (dealloc_stat /= 0)   call error_msg_dealloc
    end if 

    if (allocated(remd%grest_energy%bond_force_const) ) then
      deallocate(remd%grest_energy%bond_force_const, stat = dealloc_stat)
      if (dealloc_stat /= 0)   call error_msg_dealloc
    end if

    if (allocated(remd%grest_energy%angle_force_const) ) then
      deallocate(remd%grest_energy%angle_force_const, stat = dealloc_stat)
      if (dealloc_stat /= 0)   call error_msg_dealloc
    end if

    if (allocated(remd%grest_energy%urey_force_const) ) then
      deallocate(remd%grest_energy%urey_force_const, stat = dealloc_stat)
      if (dealloc_stat /= 0)   call error_msg_dealloc
    end if

    if (allocated(remd%grest_energy%dihe_force_const) ) then
      deallocate(remd%grest_energy%dihe_force_const, stat = dealloc_stat)
      if (dealloc_stat /= 0)   call error_msg_dealloc
    end if

    if (allocated(remd%grest_energy%dihe_force_const) ) then
      deallocate(remd%grest_energy%dihe_force_const, stat = dealloc_stat)
      if (dealloc_stat /= 0)   call error_msg_dealloc
    end if

    if (allocated(remd%grest_energy%impr_force_const) ) then
      deallocate(remd%grest_energy%impr_force_const, stat = dealloc_stat)
      if (dealloc_stat /= 0)   call error_msg_dealloc
    end if

    if (allocated(remd%grest_energy%cmap_coef) ) then
      deallocate(remd%grest_energy%cmap_coef, stat = dealloc_stat)
      if (dealloc_stat /= 0)   call error_msg_dealloc
    end if

    if (allocated(remd%grest_energy%nb14_lj6) ) then
      deallocate(remd%grest_energy%nb14_lj6, stat = dealloc_stat)
      if (dealloc_stat /= 0)   call error_msg_dealloc
    end if

    if (allocated(remd%grest_energy%nb14_lj12) ) then
      deallocate(remd%grest_energy%nb14_lj12, stat = dealloc_stat)
      if (dealloc_stat /= 0)   call error_msg_dealloc
    end if

    if (allocated(remd%grest_energy%nonb_lj6) ) then
      deallocate(remd%grest_energy%nonb_lj6, stat = dealloc_stat)
      if (dealloc_stat /= 0)   call error_msg_dealloc
    end if

    if (allocated(remd%grest_energy%nonb_lj12) ) then
      deallocate(remd%grest_energy%nonb_lj12, stat = dealloc_stat)
      if (dealloc_stat /= 0)   call error_msg_dealloc
    end if

    return

  end subroutine memory_dealloc_grest_ref
     
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    copy_grest_to_ref
  !> @brief        copy current force field to reference array
  !! @authors      JJ
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy function information
  !! @param[inout] remd        : remd information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine copy_grest_to_ref(domain, enefunc, remd)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    type(s_remd),    target, intent(inout) :: remd

    ! local variables
    integer                          :: i, icel, ix, jx, ig
    integer                          :: ncell_local, ncell
    integer                          :: ncmap_type, ncls
    integer                          :: list(4)
    type(s_soltemp),         pointer :: soltemp(:)
    integer,                 pointer :: natom(:), nsolute(:), id_l2g_sol(:,:)
    integer,                 pointer :: nbond(:), bond_list(:,:,:)
    integer,                 pointer :: nangl(:), angl_list(:,:,:)
    integer,                 pointer :: ndihe(:), dihe_list(:,:,:)
    integer,                 pointer :: rb_dihe(:), rb_dihe_list(:,:,:)
    integer,                 pointer :: nimpr(:), impr_list(:,:,:)
    integer,                 pointer :: ncontact(:), contact_list(:,:,:)
    real(wp),                pointer :: charge_ref(:,:), charge(:,:)
    real(wp),                pointer :: force_bond_ref(:,:), force_bond(:,:)
    real(wp),                pointer :: force_angl_ref(:,:), force_angl(:,:)
    real(wp),                pointer :: force_urey_ref(:,:), force_urey(:,:)
    real(wp),                pointer :: force_dihe_ref(:,:), force_dihe(:,:)
    real(wp),                pointer :: rb_dihe_c_ref(:,:,:), rb_dihe_c(:,:,:)
    real(wp),                pointer :: force_impr_ref(:,:), force_impr(:,:)
    real(wp),                pointer :: contact_lj6(:,:), contact_lj12(:,:)
    real(wp),                pointer :: contact_lj6_ref(:,:)
    real(wp),                pointer :: contact_lj12_ref(:,:)
    real(wp),                pointer :: nb14_lj6_ref(:,:), nb14_lj6(:,:)
    real(wp),                pointer :: nb14_lj12_ref(:,:), nb14_lj12(:,:)
    real(wp),                pointer :: nonb_lj6_ref(:,:), nonb_lj6(:,:)
    real(wp),                pointer :: nonb_lj12_ref(:,:), nonb_lj12(:,:)
    real(wp),                pointer :: cmap_coef_ref(:,:,:,:,:)
    real(wp),                pointer :: cmap_coef(:,:,:,:,:)
 
    soltemp          => remd%solute_tempering
    natom            => domain%num_atom
    nsolute          => domain%num_solute
    id_l2g_sol       => domain%id_l2g_solute
    charge           => domain%charge
    nbond            => enefunc%num_bond
    bond_list        => enefunc%bond_list
    nangl            => enefunc%num_angle
    angl_list        => enefunc%angle_list
    ndihe            => enefunc%num_dihedral 
    dihe_list        => enefunc%dihe_list
    rb_dihe          => enefunc%num_rb_dihedral
    rb_dihe_list     => enefunc%rb_dihe_list
    nimpr            => enefunc%num_improper 
    impr_list        => enefunc%impr_list
    ncontact         => enefunc%num_contact
    contact_list     => enefunc%contact_list
    force_bond       => enefunc%bond_force_const
    force_angl       => enefunc%angle_force_const
    force_urey       => enefunc%urey_force_const
    force_dihe       => enefunc%dihe_force_const
    force_impr       => enefunc%impr_force_const
    rb_dihe_c        => enefunc%rb_dihe_c
    contact_lj6      => enefunc%contact_lj6
    contact_lj12     => enefunc%contact_lj12
    nb14_lj6         => enefunc%nb14_lj6
    nb14_lj12        => enefunc%nb14_lj12
    nonb_lj6         => enefunc%nonb_lj6
    nonb_lj12        => enefunc%nonb_lj12
    cmap_coef        => enefunc%cmap_coef
    charge_ref       => remd%grest_energy%charge
    force_bond_ref   => remd%grest_energy%bond_force_const
    force_angl_ref   => remd%grest_energy%angle_force_const
    force_urey_ref   => remd%grest_energy%urey_force_const
    force_dihe_ref   => remd%grest_energy%dihe_force_const
    rb_dihe_c_ref    => remd%grest_energy%rb_dihe_c
    force_impr_ref   => remd%grest_energy%impr_force_const
    contact_lj6_ref  => remd%grest_energy%contact_lj6
    contact_lj12_ref => remd%grest_energy%contact_lj12
    nb14_lj6_ref     => remd%grest_energy%nb14_lj6
    nb14_lj12_ref    => remd%grest_energy%nb14_lj12
    nonb_lj6_ref     => remd%grest_energy%nonb_lj6
    nonb_lj12_ref    => remd%grest_energy%nonb_lj12
    cmap_coef_ref    => remd%grest_energy%cmap_coef
 
    ncell_local = domain%num_cell_local
    ncell       = ncell_local + domain%num_cell_boundary
 
    do i = 1, remd%dimension

      if (remd%types(i) == RemdSoluteTempering) then
       
        ! copy charge
        ! 
        if (soltemp(i)%sw_charge) then
          do icel = 1, ncell
            do ix = 1, nsolute(icel)
              ig = id_l2g_sol(ix,icel)
              if (soltemp(i)%is_solute(ig) > 0) then
                charge_ref(ix,icel) = charge(ix,icel)
              end if
            end do
          end do
        end if

        ! copy bond force constant
        !
        if (soltemp(i)%sw_bonds) then
          call copy_solute_tempering_internal(soltemp(i), ncell_local, 2, &
                                              nbond, bond_list,           &
                                              force_bond, force_bond_ref)
        end if 
        
        ! copy angle force constant
        !
        if (soltemp(i)%sw_angles) then
          call copy_solute_tempering_internal(soltemp(i), ncell_local, 3, &
                                              nangl, angl_list,           &
                                              force_angl, force_angl_ref)
        end if

        ! copy urey force constant
        !
        if (soltemp(i)%sw_ureys) then
          do icel = 1, ncell_local
            do ix = 1, nangl(icel)
              list(1) = angl_list(1,ix,icel)
              list(3) = angl_list(3,ix,icel)
              if ((soltemp(i)%is_solute(list(1)) > 0   .or.  &
                  (soltemp(i)%is_solute(list(3)) > 0)) .and. &
                  force_urey(ix,icel) > EPS) then
                force_urey_ref(ix,icel) = force_urey(ix,icel)
              end if
            end do
          end do
        end if
             
        ! copy dihedral force constant
        !
        if (soltemp(i)%sw_dihedrals) then
          call copy_solute_tempering_internal(soltemp(i), ncell_local, 4, &
                                              ndihe, dihe_list,           &
                                              force_dihe, force_dihe_ref)
        end if

        ! copy rb dihedral coefficients
        !
        if (soltemp(i)%sw_rb_dihedrals) then
          do icel = 1, ncell_local
            do ix = 1, rb_dihe(icel)
              list(1:4) = rb_dihe_list(1:4,ix,icel)
              if (soltemp(i)%is_solute(list(1)) > 0 .or. &
                  soltemp(i)%is_solute(list(2)) > 0 .or. &
                  soltemp(i)%is_solute(list(3)) > 0 .or. &
                  soltemp(i)%is_solute(list(4)) > 0) then
                rb_dihe_c_ref(1:6,ix,icel) = rb_dihe_c(1:6,ix,icel)
              end if
            end do
          end do
        end if

        ! copy improper dihedral force constant
        !
        if (soltemp(i)%sw_impropers) then
          call copy_solute_tempering_internal(soltemp(i), ncell_local, 4, &
                                              nimpr, impr_list,           &
                                              force_impr, force_impr_ref)
        end if

        ! copy contact
        !
        if (soltemp(i)%sw_contacts) then
          do icel = 1, ncell_local
            do ix = 1, ncontact(icel)
              list(1:2) = contact_list(1:2,ix,icel)      
              if (soltemp(i)%is_solute(list(1)) > 0 .or. &
                  soltemp(i)%is_solute(list(2)) > 0) then
                contact_lj6_ref(ix,icel)  = contact_lj6(ix,icel)
                contact_lj12_ref(ix,icel) = contact_lj12(ix,icel)
              end if
            end do
          end do
        end if

        ! copy lj coeffecieints
        !
        if (soltemp(i)%sw_lj) then
          ncls = enefunc%num_atom_cls
          if (allocated(enefunc%nb14_lj6)) then
            do ix = 1, ncls
              do jx = 1, ncls
                nb14_lj6_ref(jx,ix) = nb14_lj6(jx,ix)
              end do
            end do
          end if
          if (allocated(enefunc%nb14_lj12)) then
            do ix = 1, ncls
              do jx = 1, ncls
                nb14_lj12_ref(jx,ix) = nb14_lj12(jx,ix)
              end do
            end do
          end if
          if (allocated(enefunc%nonb_lj6)) then
            do ix = 1, ncls
              do jx = 1, ncls
                nonb_lj6_ref(jx,ix) = nonb_lj6(jx,ix)
              end do
            end do
          end if
          if (allocated(enefunc%nonb_lj12)) then
            do ix = 1, ncls
              do jx = 1, ncls
                nonb_lj12_ref(jx,ix) = nonb_lj12(jx,ix)
              end do
            end do
          end if
        end if

        ! copy cmap coeffcients
        !
        if (soltemp(i)%num_cmap_type > 0 .and. soltemp(i)%sw_cmaps) then
          ncmap_type = size(enefunc%cmap_coef(1,1,1,1,:))
          cmap_coef_ref(1:4,1:4,1:24,1:24,1:ncmap_type) = &
                              cmap_coef(1:4,1:4,1:24,1:24,1:ncmap_type) 
        end if
           
      end if
    end do

    return

  end subroutine copy_grest_to_ref
     
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    copy_grest_from_ref
  !> @brief        copy current force field from reference array
  !! @authors      JJ
  !! @param[inout] domain      : domain information
  !! @param[in]    remd        : remd information
  !! @param[inout] enefunc     : potential energy function information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine copy_grest_from_ref(domain, remd, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(inout) :: domain
    type(s_remd),    target, intent(in)    :: remd
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                          :: i, icel, ix, jx, ig, n
    integer                          :: ncell_local, ncell
    integer                          :: ncmap_type, ncls
    integer                          :: list(4)
    real(wp)                         :: el_fact, alpha
    type(s_soltemp),         pointer :: soltemp(:)
    integer,                 pointer :: natom(:), nsolute(:), id_l2g_sol(:,:)
    integer,                 pointer :: nbond(:), bond_list(:,:,:)
    integer,                 pointer :: nangl(:), angl_list(:,:,:)
    integer,                 pointer :: ndihe(:), dihe_list(:,:,:)
    integer,                 pointer :: rb_dihe(:), rb_dihe_list(:,:,:)
    integer,                 pointer :: nimpr(:), impr_list(:,:,:)
    integer,                 pointer :: ncontact(:), contact_list(:,:,:)
    real(wp),                pointer :: charge_ref(:,:), charge(:,:)
    real(wp),                pointer :: force_bond_ref(:,:), force_bond(:,:)
    real(wp),                pointer :: force_angl_ref(:,:), force_angl(:,:)
    real(wp),                pointer :: force_urey_ref(:,:), force_urey(:,:)
    real(wp),                pointer :: force_dihe_ref(:,:), force_dihe(:,:)
    real(wp),                pointer :: rb_dihe_c_ref(:,:,:), rb_dihe_c(:,:,:)
    real(wp),                pointer :: force_impr_ref(:,:), force_impr(:,:)
    real(wp),                pointer :: contact_lj6(:,:), contact_lj12(:,:)
    real(wp),                pointer :: contact_lj6_ref(:,:)
    real(wp),                pointer :: contact_lj12_ref(:,:)
    real(wp),                pointer :: nb14_lj6_ref(:,:), nb14_lj6(:,:)
    real(wp),                pointer :: nb14_lj12_ref(:,:), nb14_lj12(:,:)
    real(wp),                pointer :: nonb_lj6_ref(:,:), nonb_lj6(:,:)
    real(wp),                pointer :: nonb_lj12_ref(:,:), nonb_lj12(:,:)
    real(wp),                pointer :: cmap_coef_ref(:,:,:,:,:)
    real(wp),                pointer :: cmap_coef(:,:,:,:,:)
 

    soltemp          => remd%solute_tempering
    natom            => domain%num_atom
    nsolute          => domain%num_solute
    id_l2g_sol       => domain%id_l2g_solute
    charge           => domain%charge
    nbond            => enefunc%num_bond
    bond_list        => enefunc%bond_list
    nangl            => enefunc%num_angle
    angl_list        => enefunc%angle_list
    ndihe            => enefunc%num_dihedral 
    dihe_list        => enefunc%dihe_list
    rb_dihe          => enefunc%num_rb_dihedral
    rb_dihe_list     => enefunc%rb_dihe_list
    nimpr            => enefunc%num_improper 
    impr_list        => enefunc%impr_list
    ncontact         => enefunc%num_contact
    contact_list     => enefunc%contact_list
    force_bond       => enefunc%bond_force_const
    force_angl       => enefunc%angle_force_const
    force_urey       => enefunc%urey_force_const
    force_dihe       => enefunc%dihe_force_const
    force_impr       => enefunc%impr_force_const
    rb_dihe_c        => enefunc%rb_dihe_c
    contact_lj6      => enefunc%contact_lj6
    contact_lj12     => enefunc%contact_lj12
    nb14_lj6         => enefunc%nb14_lj6
    nb14_lj12        => enefunc%nb14_lj12
    nonb_lj6         => enefunc%nonb_lj6
    nonb_lj12        => enefunc%nonb_lj12
    cmap_coef        => enefunc%cmap_coef
    charge_ref       => remd%grest_energy%charge
    force_bond_ref   => remd%grest_energy%bond_force_const
    force_angl_ref   => remd%grest_energy%angle_force_const
    force_urey_ref   => remd%grest_energy%urey_force_const
    force_dihe_ref   => remd%grest_energy%dihe_force_const
    rb_dihe_c_ref    => remd%grest_energy%rb_dihe_c
    force_impr_ref   => remd%grest_energy%impr_force_const
    contact_lj6_ref  => remd%grest_energy%contact_lj6
    contact_lj12_ref => remd%grest_energy%contact_lj12
    nb14_lj6_ref     => remd%grest_energy%nb14_lj6
    nb14_lj12_ref    => remd%grest_energy%nb14_lj12
    nonb_lj6_ref     => remd%grest_energy%nonb_lj6
    nonb_lj12_ref    => remd%grest_energy%nonb_lj12
    cmap_coef_ref    => remd%grest_energy%cmap_coef
 
    ncell_local = domain%num_cell_local
    ncell       = ncell_local + domain%num_cell_boundary
 
    do i = 1, remd%dimension

      if (remd%types(i) == RemdSoluteTempering) then
       
        ! copy charge
        ! 
        if (soltemp(i)%sw_charge) then
          do icel = 1, ncell
            do ix = 1, nsolute(icel)
              ig = id_l2g_sol(ix,icel)
              if (soltemp(i)%is_solute(ig) > 0) then
                charge(ix,icel) = charge_ref(ix,icel)
              end if
            end do
          end do
#ifdef USE_GPU
          if (domain%nonbond_kernel == NBK_GPU) then
            n = domain%num_atom_domain
            do icel = 1, ncell
              ig = domain%start_atom(icel)
              do ix = 1, domain%num_atom(icel)
                domain%translated(    ig+ix,1,1) = domain%coord(1,ix,icel) + domain%trans_vec(1,ix,icel)
                domain%translated(  n+ig+ix,1,1) = domain%coord(2,ix,icel) + domain%trans_vec(2,ix,icel)
                domain%translated(2*n+ig+ix,1,1) = domain%coord(3,ix,icel) + domain%trans_vec(3,ix,icel)
                domain%translated(3*n+ig+ix,1,1) = charge(ix,icel)
              end do
            end do
            call gpu_upload_charge(domain%translated);
          end if
#endif

          ! PME self energy from the reference charge
          !
          u_self = 0.0_dp
          el_fact = ELECOEF / enefunc%dielec_const
          alpha = enefunc%pme_alpha
          do icel = 1, ncell_local
            do ix = 1, natom(icel)
              u_self = u_self + charge(ix,icel)*charge(ix,icel)
            end do
          end do
          u_self = - u_self * el_fact * alpha / sqrt(PI)

        end if

        ! copy bond force constant
        !
        if (soltemp(i)%sw_bonds) then
          call copy_solute_tempering_internal(soltemp(i), ncell_local, 2, &
                                              nbond, bond_list,           &
                                              force_bond_ref, force_bond)
        end if 
        
        ! copy angle force constant
        !
        if (soltemp(i)%sw_angles) then
          call copy_solute_tempering_internal(soltemp(i), ncell_local, 3, &
                                              nangl, angl_list,           &
                                              force_angl_ref, force_angl)
        end if

        ! copy urey force constant
        !
        if (soltemp(i)%sw_ureys) then
          do icel = 1, ncell_local
            do ix = 1, nangl(icel)
              list(1) = angl_list(1,ix,icel)
              list(3) = angl_list(3,ix,icel)
              if ((soltemp(i)%is_solute(list(1)) > 0   .or.  &
                  (soltemp(i)%is_solute(list(3)) > 0)) .and. &
                  force_urey(ix,icel) > EPS) then
                force_urey(ix,icel) = force_urey_ref(ix,icel)
              end if
            end do
          end do
        end if
             
        ! copy dihedral force constant
        !
        if (soltemp(i)%sw_dihedrals) then
          call copy_solute_tempering_internal(soltemp(i), ncell_local, 4, &
                                              ndihe, dihe_list,           &
                                              force_dihe_ref, force_dihe)
        end if

        ! copy rb dihedral coefficients
        !
        if (soltemp(i)%sw_rb_dihedrals) then
          do icel = 1, ncell_local
            do ix = 1, rb_dihe(icel)
              list(1:4) = rb_dihe_list(1:4,ix,icel)
              if (soltemp(i)%is_solute(list(1)) > 0 .or. &
                  soltemp(i)%is_solute(list(2)) > 0 .or. &
                  soltemp(i)%is_solute(list(3)) > 0 .or. &
                  soltemp(i)%is_solute(list(4)) > 0) then
                rb_dihe_c(1:6,ix,icel) = rb_dihe_c_ref(1:6,ix,icel)
              end if
            end do
          end do
        end if

        ! copy improper dihedral force constant
        !
        if (soltemp(i)%sw_impropers) then
          call copy_solute_tempering_internal(soltemp(i), ncell_local, 4, &
                                              nimpr, impr_list,           &
                                              force_impr_ref, force_impr)
        end if

        ! copy contact
        !
        if (soltemp(i)%sw_contacts) then
          do icel = 1, ncell_local
            do ix = 1, ncontact(icel)
              list(1:2) = contact_list(1:2,ix,icel)      
              if (soltemp(i)%is_solute(list(1)) > 0 .or. &
                  soltemp(i)%is_solute(list(2)) > 0) then
                contact_lj6(ix,icel)  = contact_lj6_ref(ix,icel)
                contact_lj12(ix,icel) = contact_lj12_ref(ix,icel)
              end if
            end do
          end do
        end if

        ! copy lj coeffecieints
        !
        if (soltemp(i)%sw_lj) then
          ncls = enefunc%num_atom_cls
          if (allocated(enefunc%nb14_lj6)) then
            do ix = 1, ncls
              do jx = 1, ncls
                nb14_lj6(jx,ix) = nb14_lj6_ref(jx,ix)
              end do
            end do
          end if
          if (allocated(enefunc%nb14_lj12)) then
            do ix = 1, ncls
              do jx = 1, ncls
                nb14_lj12(jx,ix) = nb14_lj12_ref(jx,ix)
              end do
            end do
          end if
          if (allocated(enefunc%nonb_lj6)) then
            do ix = 1, ncls
              do jx = 1, ncls
                nonb_lj6(jx,ix) = nonb_lj6_ref(jx,ix)
              end do
            end do
          end if
          if (allocated(enefunc%nonb_lj12)) then
            do ix = 1, ncls
              do jx = 1, ncls
                nonb_lj12(jx,ix) = nonb_lj12_ref(jx,ix)
              end do
            end do
          end if
#ifdef USE_GPU
          if (allocated(enefunc%nonb_lj12) .or. &
              allocated(enefunc%nonb_lj6)) then
            call gpu_upload_lj_coeffs(ncls, nonb_lj12, nonb_lj6);
          end if
#endif /* USE_GPU */

        end if

        ! copy cmap coeffcients
        !
        if (soltemp(i)%num_cmap_type > 0 .and. soltemp(i)%sw_cmaps) then
          ncmap_type = size(enefunc%cmap_coef(1,1,1,1,:))
          cmap_coef(1:4,1:4,1:24,1:24,1:ncmap_type) = &
                            cmap_coef_ref(1:4,1:4,1:24,1:24,1:ncmap_type) 
        end if
           
      end if
    end do

    return

  end subroutine copy_grest_from_ref

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_condition_solte_tempering
  !> @brief        reassign force constant involving REST
  !! @authors      MK
  !! @param[inout] soltemp     : REST information
  !! @param[in]    rest_param  : target temperature
  !! @param[inout] domain      : domain information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_condition_solute_tempering(soltemp, rest_param, &
                                                domain, ensemble, enefunc)

    ! formal arguments
    type(s_soltemp),  target, intent(inout) :: soltemp
    real(wp),                 intent(in)    :: rest_param
    type(s_domain),   target, intent(inout) :: domain
    type(s_ensemble),         intent(in)    :: ensemble
    type(s_enefunc),  target, intent(inout) :: enefunc

    integer,          pointer :: natom(:), nsolute(:), id_l2g_sol(:,:)
    integer(kind=1),  pointer :: is_sol(:)
    integer,          pointer :: nbond(:), bond_list(:,:,:)
    integer,          pointer :: nangl(:), angl_list(:,:,:)
    integer,          pointer :: ndihe(:), dihe_list(:,:,:)
    integer,          pointer :: rb_dihe(:), rb_dihe_list(:,:,:)
    integer,          pointer :: nimpr(:), impr_list(:,:,:)
    integer,          pointer :: ncontact(:), contact_list(:,:,:)
    real(wp),         pointer :: charge(:,:)
    real(wp),         pointer :: force_bond(:,:)
    real(wp),         pointer :: force_angl(:,:), force_urey(:,:)
    real(wp),         pointer :: force_dihe(:,:)
    real(wp),         pointer :: rb_dihe_c(:,:,:)
    real(wp),         pointer :: force_impr(:,:)
    real(wp),         pointer :: contact_lj6(:,:), contact_lj12(:,:)
    real(wp),         pointer :: nb14_lj12(:,:), nb14_lj6(:,:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: cmap_coef(:,:,:,:,:)
    integer                   :: alist(8), i, j, ix, num, tgt, n, ig
    integer                   :: org, ncell, ncell_local, counter
    real(wp)                  :: coeff_full, coeff_half, wgt, tmp1, tmp2
    real(wp)                  :: el_fact, alpha
    real(dp)                  :: u_self_org


    id_l2g_sol   => domain%id_l2g_solute
    natom        => domain%num_atom
    nsolute      => domain%num_solute
    charge       => domain%charge
    nbond        => enefunc%num_bond
    bond_list    => enefunc%bond_list
    force_bond   => enefunc%bond_force_const
    nangl        => enefunc%num_angle
    angl_list    => enefunc%angle_list
    force_angl   => enefunc%angle_force_const
    force_urey   => enefunc%urey_force_const
    ndihe        => enefunc%num_dihedral
    dihe_list    => enefunc%dihe_list
    force_dihe   => enefunc%dihe_force_const
    rb_dihe      => enefunc%num_rb_dihedral
    rb_dihe_list => enefunc%rb_dihe_list
    rb_dihe_c    => enefunc%rb_dihe_c
    nimpr        => enefunc%num_improper
    impr_list    => enefunc%impr_list
    force_impr   => enefunc%impr_force_const
    ncontact     => enefunc%num_contact
    contact_list => enefunc%contact_list
    contact_lj12 => enefunc%contact_lj12
    contact_lj6  => enefunc%contact_lj6
    nb14_lj12    => enefunc%nb14_lj12 
    nb14_lj6     => enefunc%nb14_lj6  
    nonb_lj12    => enefunc%nonb_lj12 
    nonb_lj6     => enefunc%nonb_lj6  
    cmap_coef    => enefunc%cmap_coef
    is_sol       => soltemp%is_solute
    ncell_local  =  domain%num_cell_local
    ncell        =  ncell_local + domain%num_cell_boundary

    coeff_full = ensemble%temperature / rest_param
    coeff_half = sqrt( coeff_full )

    if ( soltemp%done_setup ) then
      tmp1 = coeff_full
      tmp2 = coeff_half

      coeff_full = coeff_full / soltemp%rest_param_full
      coeff_half = coeff_half / soltemp%rest_param_half

      ! remember coeffs for next exchange
      !
      soltemp%rest_param_full = tmp1
      soltemp%rest_param_half = tmp2
    else
      soltemp%rest_param_full = coeff_full
      soltemp%rest_param_half = coeff_half
      soltemp%done_setup = .true.
    end if

    if (soltemp%sw_charge) then
      do i = 1, ncell
        ! charge
        do ix = 1, nsolute(i)
          ig = id_l2g_sol(ix,i)
          if ( is_sol(ig) > 0 ) then
            charge(ix,i) = coeff_half * charge(ix,i)
          end if
        end do
      end do
#ifdef USE_GPU
      if (domain%nonbond_kernel == NBK_GPU) then
        n = domain%num_atom_domain
        do i = 1, ncell
          ig = domain%start_atom(i)
          do ix = 1, domain%num_atom(i)
            domain%translated(    ig+ix,1,1) = domain%coord(1,ix,i) + domain%trans_vec(1,ix,i)
            domain%translated(  n+ig+ix,1,1) = domain%coord(2,ix,i) + domain%trans_vec(2,ix,i)
            domain%translated(2*n+ig+ix,1,1) = domain%coord(3,ix,i) + domain%trans_vec(3,ix,i)
            domain%translated(3*n+ig+ix,1,1) = charge(ix,i)
          end do
        end do
        call gpu_upload_charge(domain%translated);
      end if
#endif

      ! need to update PME self energy
      !
      u_self = 0.0_dp
      el_fact = ELECOEF / enefunc%dielec_const
      alpha = enefunc%pme_alpha
      do i = 1, ncell_local
        do ix = 1, natom(i)
          u_self = u_self + charge(ix,i)*charge(ix,i)
        end do
      end do
      u_self = - u_self * el_fact * alpha / sqrt(PI)
    end if

    do i = 1, ncell_local

      ! bond
      if (soltemp%sw_bonds) then
        call assign_condition_solute_tempering_internal( &
                  soltemp, 2, coeff_full,                &
                  nbond(i), bond_list(:,:,i), force_bond(:,i))
      end if

      ! angle and urey
      do ix = 1, nangl(i)
        alist(1:3) = angl_list(1:3,ix,i)
        num = 0
        do j = 1, 3
          if (is_sol(alist(j)) > 0) num = num + 1
        end do
        if (num > 0 .and. soltemp%sw_angles) then
          wgt = real(num,wp) / 3.0_wp
          wgt = coeff_full ** wgt
          force_angl(ix,i) = wgt * force_angl(ix,i)
        end if
        if (force_urey(ix,i) > EPS .and. soltemp%sw_ureys) then
          num = 0
          if (is_sol(alist(1)) > 0) num = num + 1
          if (is_sol(alist(3)) > 0) num = num + 1
          if (num > 0) then
            wgt = real(num,wp) / 2.0_wp
            wgt = coeff_full ** wgt
            force_urey(ix,i) = wgt * force_urey(ix,i)
          end if
        end if
      end do

      ! dihedral
      if (soltemp%sw_dihedrals) then
        call assign_condition_solute_tempering_internal( &
                  soltemp, 4, coeff_full,                &
                  ndihe(i), dihe_list(:,:,i), force_dihe(:,i))
      end if

      ! rb_dihedral
      if (soltemp%sw_rb_dihedrals) then
        do ix = 1, rb_dihe(i)
          alist(1:4) = rb_dihe_list(1:4,ix,i)
          num = 0
          do j = 1, 4
            if (is_sol(alist(j)) > 0) num = num + 1
          end do
          if (num > 0) then
            wgt = real(num,wp) / real(n,wp)
            wgt = coeff_full ** wgt
            rb_dihe_c(1:6,ix,i) = wgt * rb_dihe_c(1:6,ix,i)
          end if
        end do
      end if

      ! impropers
      if (soltemp%sw_impropers) then
        call assign_condition_solute_tempering_internal( &
                  soltemp, 4, coeff_full,                &
                  nimpr(i), impr_list(:,:,i), force_impr(:,i))
      end if

      ! contact
      if (soltemp%sw_contacts) then
        do ix = 1, ncontact(i)
          alist(1:2) = contact_list(1:2,ix,i)
          num = 0
          do j = 1, 2
            if (is_sol(alist(j)) > 0) num = num + 1
          end do
          if (num > 0) then
            wgt = real(num,wp) / 2.0_wp
            wgt = coeff_full ** wgt
            contact_lj6(ix,i)  = wgt * contact_lj6(ix,i)
            contact_lj12(ix,i) = wgt * contact_lj12(ix,i)
          end if
        end do
      end if

    end do

    ! lj
    if (soltemp%sw_lj) then
      n = enefunc%num_atom_cls
      if (allocated(enefunc%nb14_lj6)) then
        call assign_condition_solute_tempering_lj( &
                 n, soltemp%rest_param_half, soltemp%rest_param_full, &
                 soltemp, nb14_lj6)
      end if
      if (allocated(enefunc%nb14_lj12)) then
        call assign_condition_solute_tempering_lj( &
                 n, soltemp%rest_param_half, soltemp%rest_param_full, &
                 soltemp, nb14_lj12)
      end if
      if (allocated(enefunc%nonb_lj6)) then
        call assign_condition_solute_tempering_lj( &
                 n, soltemp%rest_param_half, soltemp%rest_param_full, &
                 soltemp, nonb_lj6)
      end if
      if (allocated(enefunc%nonb_lj12)) then
        call assign_condition_solute_tempering_lj( &
                 n, soltemp%rest_param_half, soltemp%rest_param_full, &
                 soltemp, nonb_lj12)
      end if
#ifdef USE_GPU
      if (allocated(enefunc%nonb_lj12) .or. allocated(enefunc%nonb_lj6)) then
        call gpu_upload_lj_coeffs(n, enefunc%nonb_lj12, enefunc%nonb_lj6);
      end if
#endif /* USE_GPU */
    end if

    ! cmap
    if (soltemp%num_cmap_type > 0 .and. soltemp%sw_cmaps) then
      do i = 1, soltemp%num_cmap_type
        tgt = i + soltemp%istart_cmap_type - 1
        org = soltemp%cmap_type_org(i)
        wgt = soltemp%rest_param_full ** soltemp%cmap_weight(i)
        cmap_coef(1:4,1:4,1:24,1:24,tgt) = wgt*cmap_coef(1:4,1:4,1:24,1:24,org)
      end do
    end if

    return

  end subroutine assign_condition_solute_tempering

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_condition_solute_tempering_lj
  !> @brief        reassign force constant involving REST for LJ
  !! @authors      MK
  !! @param[in]    n           : array size of interaction matrix
  !! @param[in]    coeff_half  : coefficient for solute-solvent
  !! @param[in]    coeff_full  : coefficient for solute-solute
  !! @param[in]    soltemp     : REST information
  !! @param[inout] nbcoeff     : LJ interaction matrix
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_condition_solute_tempering_lj(n, coeff_half, coeff_full, &
                                                  soltemp, nbcoeff)

    ! formal arguments
    integer,                 intent(in)    :: n
    real(wp),                intent(in)    :: coeff_half
    real(wp),                intent(in)    :: coeff_full
    type(s_soltemp),         intent(in)    :: soltemp
    real(wp),                intent(inout) :: nbcoeff(n,n)

    ! local variables
    integer :: i, j, oldcount, newcount, org, orgj


    oldcount = soltemp%istart_atom_cls - 1
    newcount = oldcount + soltemp%num_atom_cls

    do i = 1, soltemp%num_atom_cls
      org = soltemp%atom_cls_no_org(i)
      do j = 1, oldcount
        nbcoeff(j,i+oldcount) = coeff_half * nbcoeff(j,org)
        nbcoeff(i+oldcount,j) = coeff_half * nbcoeff(j,org)
      end do
      do j = oldcount + 1, newcount
        orgj = soltemp%atom_cls_no_org(j-oldcount)
        nbcoeff(j,i+oldcount) = coeff_full * nbcoeff(orgj,org)
        nbcoeff(i+oldcount,j) = coeff_full * nbcoeff(orgj,org)
      end do
    end do

    return

  end subroutine assign_condition_solute_tempering_lj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_condition_solte_tempering_internal
  !> @brief        reassign force constant involving REST for internal
  !! @authors      MK
  !! @param[in]    soltemp     : solute tempering information
  !! @param[in]    n           : number of indexes involved
  !! @param[in]    coeff_full  : coefficient
  !! @param[in]    n_internal  : number of terms
  !! @param[in]    aindex      : atom list
  !! @param[inout] fc          : force constant
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_condition_solute_tempering_internal(soltemp, n, &
                                  coeff_full, n_internal, aindex, fc)

    ! formal arguments
    type(s_soltemp),         intent(in)    :: soltemp
    integer,                 intent(in)    :: n
    real(wp),                intent(in)    :: coeff_full
    integer,                 intent(in)    :: n_internal
    integer,                 intent(in)    :: aindex(:,:)
    real(wp),                intent(inout) :: fc(:)

    ! local variables
    integer  :: alist(1:n), ix, j, num
    real(wp) :: wgt


    do ix = 1, n_internal
      alist(1:n) = aindex(1:n,ix)
      num = 0
      do j = 1, n
        if (soltemp%is_solute(alist(j)) > 0) num = num + 1
      end do
      if (num > 0) then
        wgt = real(num,wp) / real(n,wp)
        wgt = coeff_full ** wgt
        fc(ix) = wgt * fc(ix)
      end if
    end do

    return

  end subroutine assign_condition_solute_tempering_internal

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    copy_solute_tempering_internal
  !> @brief        copy the internal force constant
  !! @authors      JJ
  !! @param[in]    soltemp       : solute tempering information
  !! @param[in]    ncell         : number of cells
  !! @param[in]    num           : number of list
  !! @param[in]    ninternal     : number of terms
  !! @param[in]    internal_list : internal list
  !! @param[in]    force_from    : target force
  !! @param[inout] force_to      : reference force
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine copy_solute_tempering_internal(soltemp, ncell, num,       &
                                            ninternal, internal_list,  &
                                            force_from, force_to)

    ! formal arguments
    type(s_soltemp), target, intent(in)    :: soltemp
    integer,                 intent(in)    :: ncell
    integer,                 intent(in)    :: num
    integer,                 intent(in)    :: ninternal(:)
    integer,                 intent(in)    :: internal_list(:,:,:)
    real(wp),                intent(in)    :: force_from(:,:)
    real(wp),                intent(inout) :: force_to(:,:)

    ! local variables
    integer                          :: j, icel, ix, k
    integer                          :: list(num)


    do icel = 1, ncell
 
      do ix = 1, ninternal(icel)

        list(1:num) = internal_list(1:num,ix,icel)
        k = 0
        do j = 1, num
          if (soltemp%is_solute(list(j)) > 0) k = k + 1
        end do
        if (k > 0) force_to(ix,icel) = force_from(ix,icel)

      end do
    end do

    return 
  
  end subroutine copy_solute_tempering_internal

end module sp_grest_energy_mod
