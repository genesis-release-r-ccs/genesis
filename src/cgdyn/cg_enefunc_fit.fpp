!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_enefunc_fit_mod
!> @brief   Provides utilities for fitting.
!! @authors Koichi Tamura (KT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8
  
#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_enefunc_fit_mod

  use cg_restraints_str_mod
  use cg_enefunc_str_mod
  use cg_domain_str_mod
  use molecules_str_mod
  use fitting_mod
  use fitting_str_mod
  use string_mod
  use select_mod
  use select_atoms_mod
  use select_atoms_str_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public   :: setup_fitting_cgdyn
  public   :: fitting_sel
  public   :: fitting_sel_2d
  private  :: setup_enefunc_fit_domain
  private  :: setup_enefunc_fit_refcoord
  private  :: reduce_com

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_fitting_cgdyn
  !> @brief        Sets up fitting for CGDYN.
  !! @authors      CK
  !! @param[in]    fit_info   : fitting information
  !! @param[in]    sel_info   : selection information
  !! @param[in]    domain     : domain   information
  !! @param[in]    molecule   : molecule information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine setup_fitting_cgdyn(fit_info, sel_info,  &
                                 domain, molecule, enefunc)

    ! formal arguments
    type(s_fit_info),        intent(in)    :: fit_info
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_domain),          intent(in)    :: domain
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                                :: i
    integer                                :: fitting_atom_idx
    logical                                :: fit_check
    type(s_selatoms)                       :: selatoms
    
    enefunc%fitting_method  = fit_info%fitting_method

    if (enefunc%pressure_rmsd) then
      if (enefunc%fitting_method .ne. FittingMethodTR .and.      &
          enefunc%fitting_method .ne. FittingMethodTR_ROT .and.  &
          enefunc%fitting_method .ne. FittingMethodTR_ZROT) then
          call error_msg('Error: Setup_Fitting_cgdyn - the pressure_rmsd '//&
                         'option is only allowed when Translation fitting '//&
                         'is used.')
      endif
    endif

    if (enefunc%fitting_method == FittingMethodNO) then
      if (main_rank) then
        fit_check = .false.
        do i = 1,enefunc%num_restraintfuncs
          if (enefunc%restraint_kind(i) == RestraintsFuncRMSD .or. &
            enefunc%restraint_kind(i) == RestraintsFuncRMSDCOM) then
            fit_check = .true.
            exit
          endif
        end do
        if (fit_check .and. .not. fit_info%force_no_fitting) then
          call error_msg('Error: Setup_Fitting_cgdyn - Fitting is required '//&
                         'for RMSD restraints. Please specify a fitting '//&
                         'method.')
        endif
        if (fit_check .and. fit_info%force_no_fitting) then
          write(MsgOut,'(A)') "Setup_Fitting_cgdyn> RMSD restraint without fitting."
          write(MsgOut,*) 
        endif
      endif
      do i = 1,enefunc%num_restraintfuncs
        if (enefunc%restraint_kind(i) .eq. RestraintsFuncPC .or. &
            enefunc%restraint_kind(i) .eq. RestraintsFuncPCCOM) then
           call error_msg('Error: Setup_Fitting_cgdyn - Fitting is required '//&
                          'for PC/PCCOM. Please specify a fitting method.')
        endif
      end do
      return
    endif

    fit_check       = .false.
    do i = 1,enefunc%num_restraintfuncs
      if (enefunc%restraint_kind(i) == RestraintsFuncRMSD .or. &
          enefunc%restraint_kind(i) == RestraintsFuncRMSDCOM .or. &
          enefunc%restraint_kind(i) == RestraintsFuncPC .or. &
          enefunc%restraint_kind(i) == RestraintsFuncPCCOM) then
          fit_check = .true.
          exit
      endif
    end do

    if (.not. fit_check) then
       if (main_rank) then
         write(MsgOut,'(A)') "Setup_Fitting_cgdyn> NO fitting is applied, skip"
         write(MsgOut,'(A)') "  fitting method  =  NO"
         write(MsgOut,*)
       endif
       enefunc%fitting_method=FittingMethodNO
       return
    endif

    if (enefunc%fitting_method /= FittingMethodTR_ROT .and. &
        enefunc%fitting_method /= FittingMethodXYTR_ZROT) &
      call error_msg('Error: Setup_Fitting_cgdyn - Only NO/TR+ROT/XYTR+ZROT '//&
             'is allowed.')

    fitting_atom_idx      = fit_info%fitting_atom

    call select_atom(molecule, sel_info%groups(fitting_atom_idx), selatoms)
    enefunc%mass_weight   = fit_info%mass_weight

    enefunc%fitting_file  = FittingFileREF
    enefunc%fitting_move  = FittingMoveREF

    do i = 1,enefunc%num_restraintfuncs
      if (enefunc%restraint_kind(i) == RestraintsFuncRMSDCOM .and. &
        (.not. enefunc%mass_weight) ) &
        call error_msg('Error: Setup_Fitting_cgdyn - RESTRAINTS and FITTING '//&
	               'are inconsistent. Please check the [FITTING] section.')
      if (enefunc%restraint_kind(i) .eq. RestraintsFuncPCCOM .and. &
        (.not. enefunc%mass_weight) ) &
        call error_msg('Error: Setup_Fitting_cgdyn - RESTRAINTS and FITTING '//&
	               'are inconsistent. Please check the [FITTING] section.')
      if (enefunc%restraint_kind(i) == RestraintsFuncRMSD .and. &
        (enefunc%mass_weight) ) &
        call error_msg('Error: Setup_Fitting_cgdyn - RESTRAINTS and FITTING '//&
	               'are inconsistent. Please check the [FITTING] section.')
      if (enefunc%restraint_kind(i) .eq. RestraintsFuncPC .and. &
        (enefunc%mass_weight) ) &
        call error_msg('Error: Setup_Fitting_cgdyn - RESTRAINTS and FITTING '//&
	               'are inconsistent. Please check the [FITTING] section.')      
      if (enefunc%restraint_kind(i) == RestraintsFuncRMSD .or.    &
          enefunc%restraint_kind(i) == RestraintsFuncRMSDCOM .or. &
          enefunc%restraint_kind(i) == RestraintsFuncPC .or.      &
          enefunc%restraint_kind(i) == RestraintsFuncPCCOM) then
        if (fitting_atom_idx /= enefunc%restraint_grouplist(1,i)) then
          call error_msg('Error: Setup_Fitting_cgdyn - Fitting group should '//&
		         'be identical to the group in restraint.')
        endif
      endif
    end do

    call setup_enefunc_fit_refcoord(enefunc%fitting_file, molecule, enefunc)

    call setup_enefunc_fit_domain(domain, enefunc, selatoms%idx)

    return

  end subroutine setup_fitting_cgdyn

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_fit_domain
  !> @brief        Sets up fitting information about the domain.
  !! @authors      KT
  !! @param[in]    domain        : domain information
  !! @param[inout] enefunc       : energy potential functions information
  !! @param[in]    fitting_index : atom indices for fitting
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_fit_domain(domain, enefunc, fitting_index)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc
    integer,                 intent(in)    :: fitting_index(:)

    ! local variable
    integer                  :: i, iatm, icel, step, nfit

    integer,         pointer :: atom_2_cell(:), id_g2l(:), ncell
    real(wp),        pointer :: fit_refcoord(:,:), fit_coord(:,:)


    ncell           => domain%num_cell_local
    id_g2l          => domain%id_g2l
    atom_2_cell     => domain%atom_2_cell

    fit_refcoord    => enefunc%fit_refcoord

    call alloc_enefunc(enefunc, EneFuncFitd, ncell)

    fit_coord       => enefunc%fit_coord

    do step = 1, 2

      nfit = 0

      do i = 1, size(fitting_index)
        iatm = fitting_index(i)
        icel = atom_2_cell(id_g2l(iatm))
        if (icel > 0 .and. icel <= ncell) then
          nfit = nfit + 1
          if (step == 2) then
            enefunc%fitting_atom(nfit) = iatm
            enefunc%fit_coord(1:3,nfit) = fit_refcoord(1:3,iatm)
          end if
        end if
      end do

      if (step == 1) then
        enefunc%num_fit_domain = nfit
        call mpi_allreduce(nfit, MaxFit, 1, mpi_integer, &
                           mpi_max, mpi_comm_country, ierror)
        MaxFit = MaxFit * 2
        call alloc_enefunc(enefunc, EneFuncFitd, MaxFit)
      end if

    end do

    return

  end subroutine setup_enefunc_fit_domain

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_fit_refcoord
  !> @brief        Sets up energy function fitting reference coordinates.
  !! @authors      KT, CK
  !! @param[in]    fitting_file : fitfile or reffile
  !! @param[in]    molecule     : molecule information
  !! @param[inout] enefunc      : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_fit_refcoord(fitting_file, molecule, enefunc)

    ! formal arguments
    integer,                 intent(in)    :: fitting_file
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: j


    if (fitting_file == FittingFileFIT) then
      if (size(molecule%atom_fitcoord(1,:)) /= molecule%num_atoms) &
      call error_msg('Error: Setup_Enefunc_Fit_Refcoord - Invalid fitfile '//&
                     'in [INPUT].')

    else if (fitting_file == FittingFileREF) then
      if (size(enefunc%restraint_refcoord(1,:)) /= molecule%num_atoms) &
      call error_msg('Error: Setup_Enefunc_Fit_Refcoord - Invalid reffile '//&
                     'in [INPUT].')
     
    endif

    enefunc%do_fitting = .true.

    call alloc_enefunc(enefunc, EneFuncFitc, molecule%num_atoms)

    if (fitting_file == FittingFileFIT) then

      do j = 1, molecule%num_atoms
        enefunc%fit_refcoord(1,j) = molecule%atom_fitcoord(1,j)
        enefunc%fit_refcoord(2,j) = molecule%atom_fitcoord(2,j)
        enefunc%fit_refcoord(3,j) = molecule%atom_fitcoord(3,j)
      end do

    else

      do j = 1, molecule%num_atoms
        enefunc%fit_refcoord(1,j) = enefunc%restraint_refcoord(1,j)
        enefunc%fit_refcoord(2,j) = enefunc%restraint_refcoord(2,j)
        enefunc%fit_refcoord(3,j) = enefunc%restraint_refcoord(3,j)
      end do

    endif

    return

  end subroutine setup_enefunc_fit_refcoord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    fitting_sel
  !> @brief        Performs fitting.
  !! @authors      KT
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine fitting_sel(domain, enefunc, coord, rotated_coord)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wip),               intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: rotated_coord(:,:)

    ! local variables
    integer                  :: i, i1
    integer,         pointer :: id_g2l(:)
    integer,         pointer :: restraint_atom(:)
    integer,         pointer :: fitting_atom(:)
  
    integer                  :: ig, ierr
    integer                  :: ikind
    real(dp)                 :: m
    real(wip),       pointer :: mass(:)
    real(wp),        pointer :: restraint_coord(:,:), restraint_refcoord(:,:)
    real(wp),        pointer :: fit_coord(:,:), fit_refcoord(:,:)
    real(dp)                 :: com_ref(1:3), com_mov(1:3), tot_mass
    real(dp)                 :: sym_matrix(1:4,1:4)
    real(dp)                 :: rot_matrix(1:3,1:3)
    real(dp)                 :: eval(1:4), evec(1:4), work(1:11)
    real(dp)                 :: dref(1:3), dmov(1:3), dsub(1:3), dadd(1:3)

  
    id_g2l          => domain%id_g2l
    mass            => domain%mass

    restraint_atom  => enefunc%restraint_atom
    restraint_coord => enefunc%restraint_coord
    restraint_refcoord => enefunc%restraint_refcoord
    fit_refcoord    => enefunc%fit_refcoord
    
    fitting_atom    => enefunc%fitting_atom
    fit_coord       => enefunc%fit_coord

    ikind           = enefunc%fitting_move
    
    com_ref(1:3) = 0.0_dp
    com_mov(1:3) = 0.0_dp
    tot_mass     = 0.0_dp

    do i = 1, enefunc%num_fit_domain
      ig = fitting_atom(i)
      i1 = id_g2l(ig) 
      !m = mass(i1,i)
      m = 1.0_dp
      if (enefunc%mass_weight) m = mass(i1)
      select case (ikind)
      case (FittingMoveSYS)
        com_mov(1:3) = com_mov(1:3) + m * coord(i1,1:3)
        com_ref(1:3) = com_ref(1:3) + m * fit_refcoord(1:3,ig)
      case (FittingMoveREF) 
        com_mov(1:3) = com_mov(1:3) + m * restraint_refcoord(1:3,ig)
        com_ref(1:3) = com_ref(1:3) + m * coord(i1,1:3)
      end select
      tot_mass       = tot_mass     + m
    end do
    call reduce_com(com_ref, com_mov, tot_mass)
    com_ref(1:3) = com_ref(1:3) / tot_mass
    com_mov(1:3) = com_mov(1:3) / tot_mass

    ! symmetric matrix
    !
    sym_matrix(1:4,1:4) = 0.0_dp 

    do i = 1, enefunc%num_fit_domain
      ig = fitting_atom(i)
      i1 = id_g2l(ig)
      m = 1.0_dp
      if (enefunc%mass_weight) m = mass(i1)
      select case (ikind)
      case (FittingMoveSYS) 
        dmov(1:3) = coord(i1,1:3)                - com_mov(1:3)
        dref(1:3) = fit_refcoord(1:3,ig) - com_ref(1:3)
      case (FittingMoveREF) 
        dmov(1:3) = restraint_refcoord(1:3,ig)   - com_mov(1:3)
        dref(1:3) = coord(i1,1:3)                - com_ref(1:3)
      end select
      dsub(1:3) = m * (dref(1:3) - dmov(1:3))
      dadd(1:3) = m * (dref(1:3) + dmov(1:3))
      sym_matrix(1,1) = sym_matrix(1,1) + dsub(1)*dsub(1) &
                                        + dsub(2)*dsub(2) &
                                        + dsub(3)*dsub(3)
      sym_matrix(1,2) = sym_matrix(1,2) + dadd(2)*dsub(3) - dsub(2)*dadd(3)
      sym_matrix(1,3) = sym_matrix(1,3) + dsub(1)*dadd(3) - dadd(1)*dsub(3)
      sym_matrix(1,4) = sym_matrix(1,4) + dadd(1)*dsub(2) - dsub(1)*dadd(2)
      sym_matrix(2,2) = sym_matrix(2,2) + dsub(1)*dsub(1) &
                                        + dadd(2)*dadd(2) &
                                        + dadd(3)*dadd(3)
      sym_matrix(2,3) = sym_matrix(2,3) + dsub(1)*dsub(2) - dadd(1)*dadd(2)
      sym_matrix(2,4) = sym_matrix(2,4) + dsub(1)*dsub(3) - dadd(1)*dadd(3)
      sym_matrix(3,3) = sym_matrix(3,3) + dadd(1)*dadd(1) &
                                        + dsub(2)*dsub(2) &
                                        + dadd(3)*dadd(3)
      sym_matrix(3,4) = sym_matrix(3,4) + dsub(2)*dsub(3) - dadd(2)*dadd(3)
      sym_matrix(4,4) = sym_matrix(4,4) + dadd(1)*dadd(1) &
                                        + dadd(2)*dadd(2) &
                                        + dsub(3)*dsub(3)
    end do

    sym_matrix(2,1) = sym_matrix(1,2)
    sym_matrix(3,1) = sym_matrix(1,3)
    sym_matrix(3,2) = sym_matrix(2,3)
    sym_matrix(4,1) = sym_matrix(1,4)
    sym_matrix(4,2) = sym_matrix(2,4)
    sym_matrix(4,3) = sym_matrix(3,4)
    call mpi_allreduce(mpi_in_place, sym_matrix, 16, mpi_real8, mpi_sum, &
                       mpi_comm_country, ierror)

    ! eigenvalue and eigenvector
    !
#ifdef LAPACK
    call dsyev('V', 'U', 4, sym_matrix, 4, eval, work, 11, ierr)
#else
    call error_msg('Error: Fit_Trrot - Rotation for RMSD requires LAPACK.')
#endif

    if (ierr /= 0) then
      ierr = -2
      return
    end if

    ! rotation matrix
    !
    evec(1:4) = sym_matrix(1:4,1)
    
    rot_matrix(1,1) =           evec(1)*evec(1) + evec(2)*evec(2)  &
                               -evec(3)*evec(3) - evec(4)*evec(4)
    rot_matrix(1,2) = 2.0_wp * (evec(2)*evec(3) + evec(1)*evec(4))
    rot_matrix(1,3) = 2.0_wp * (evec(2)*evec(4) - evec(1)*evec(3))
    rot_matrix(2,1) = 2.0_wp * (evec(2)*evec(3) - evec(1)*evec(4))
    rot_matrix(2,2) =           evec(1)*evec(1) - evec(2)*evec(2)  &
                               +evec(3)*evec(3) - evec(4)*evec(4)
    rot_matrix(2,3) = 2.0_wp * (evec(3)*evec(4) + evec(1)*evec(2))
    rot_matrix(3,1) = 2.0_wp * (evec(2)*evec(4) + evec(1)*evec(3))
    rot_matrix(3,2) = 2.0_wp * (evec(3)*evec(4) - evec(1)*evec(2))
    rot_matrix(3,3) =           evec(1)*evec(1) - evec(2)*evec(2)  &
                               -evec(3)*evec(3) + evec(4)*evec(4)

    do i = 1, enefunc%num_rest_domain
      ig = restraint_atom(i)
      i1 = id_g2l(ig)
      select case (ikind)
      case (FittingMoveSYS) 
        dmov(1:3) = coord(i1,1:3) - com_mov(1:3)
      case (FittingMoveREF) 
        dmov(1:3) = restraint_refcoord(1:3,ig) - com_mov(1:3)
      end select
      rotated_coord(1:3,i) = real(rot_matrix(1:3,1)*dmov(1),wp) + &
                                  real(rot_matrix(1:3,2)*dmov(2),wp) + &
                                  real(rot_matrix(1:3,3)*dmov(3),wp)
      rotated_coord(1:3,i) = rotated_coord(1:3,i) + real(com_ref(1:3),wp)
    end do

    return

  end subroutine fitting_sel

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    fitting_sel_2d
  !> @brief        Performs fitting in XY space.
  !! @authors      KT
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine fitting_sel_2d(domain, enefunc, coord, rotated_coord)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wip),               intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: rotated_coord(:,:)

    ! local variables
    integer                  :: i, ix, i1, start_i
    integer,         pointer :: id_g2l(:)
    integer,         pointer :: restraint_atom(:)
    integer,         pointer :: fitting_atom(:)
  
    integer                  :: ig, ierr
    integer                  :: ikind
    real(dp)                 :: m
    real(wip),       pointer :: mass(:)
    real(wp),        pointer :: restraint_coord(:,:), restraint_refcoord(:,:)
    real(wp),        pointer :: fit_coord(:,:)
    real(dp)                 :: com_ref(1:3), com_mov(1:3), tot_mass
    real(dp)                 :: detu, detvt, sig
    real(dp)                 :: cor_matrix(1:2,1:2)
    real(dp)                 :: rot_matrix(1:2,1:2)
    real(dp)                 :: s(1:2), u(1:2,1:2), vt(1:2,1:2), work(1:10)
    real(dp)                 :: dref(1:3), dmov(1:3)
  
    id_g2l          => domain%id_g2l
    mass            => domain%mass

    restraint_atom  => enefunc%restraint_atom
    restraint_coord => enefunc%restraint_coord
    restraint_refcoord => enefunc%restraint_refcoord
    
    fitting_atom    => enefunc%fitting_atom
    fit_coord       => enefunc%fit_coord

    ikind           = enefunc%fitting_move
    
    com_ref(1:3) = 0.0_dp
    com_mov(1:3) = 0.0_dp
    tot_mass     = 0.0_dp
    do i = 1, enefunc%num_fit_domain
      ig = fitting_atom(i)
      i1 = id_g2l(ig)
      m = 1.0_dp
      if (enefunc%mass_weight) m = mass(i1)
      select case (ikind)
      case (FittingMoveSYS)
        com_mov(1:3) = com_mov(1:3) + m * coord(i1,1:3)
        com_ref(1:3) = com_ref(1:3) + m * enefunc%fit_refcoord(1:3,ig)
      case (FittingMoveREF) 
        com_mov(1:3) = com_mov(1:3) + m * restraint_refcoord(1:3,ig)
        com_ref(1:3) = com_ref(1:3) + m * coord(i1,1:3)
      end select
      tot_mass       = tot_mass     + m
    end do
    call reduce_com(com_ref, com_mov, tot_mass)
    com_ref(1:3) = com_ref(1:3) / tot_mass
    com_mov(1:3) = com_mov(1:3) / tot_mass

    ! correlation matrix
    !
    cor_matrix(1:2,1:2) = 0.0_dp 

    do i = 1, enefunc%num_fit_domain
      ig = fitting_atom(i)
      i1 = id_g2l(ig)
      m = 1.0_dp
      select case (ikind)
      case (FittingMoveSYS)
        dmov(1:2) = coord(i1,1:2)                - com_mov(1:2)
        dref(1:2) = enefunc%fit_refcoord(1:2,ig) - com_ref(1:2)
      case (FittingMoveREF)
        dmov(1:2) = restraint_refcoord(1:2,ig)   - com_mov(1:2)
        dref(1:2) = coord(i1,1:2)                - com_ref(1:2)
      end select
      cor_matrix(1,1) = cor_matrix(1,1) + dmov(1) * dref(1)
      cor_matrix(1,2) = cor_matrix(1,2) + dmov(2) * dref(1)
      cor_matrix(2,1) = cor_matrix(2,1) + dmov(1) * dref(2)
      cor_matrix(2,2) = cor_matrix(2,2) + dmov(2) * dref(2)
    end do

    call mpi_allreduce(mpi_in_place, cor_matrix, 4, mpi_real8, mpi_sum, &
                       mpi_comm_country, ierror)

    ! eigenvalue and eigenvector
    !
#ifdef LAPACK
    call dgesvd('A', 'A', 2, 2, cor_matrix, 2, s, u, 2, vt, 2, work, 10, ierr)
#else
    call error_msg('Error: Fit_Trrot - Rotation for RMSD requires LAPACK.')
#endif

    if (ierr /= 0) then
      ierr = -2
      return
    end if

    ! check reflection
    !
    detu  =  u(1,1) *  u(2,2) - u(1,2) *  u(2,1)
    detvt = vt(1,1) * vt(2,2) -vt(1,2) * vt(2,1)
    sig   = detu * detvt
    if (sig < 0.0_dp) then
      s(2)   = -s(2)
      u(:,2) = -u(:,2)
    end if


    ! rotation matrix
    !
    rot_matrix = matmul(u,vt)
    
    do i = 1, enefunc%num_rest_domain
      ig = restraint_atom(i)
      i1 = id_g2l(ig)
      select case (ikind)
      case (FittingMoveSYS)
        dmov(1:2) = coord(i1,1:2) - com_mov(1:2)
        rotated_coord(3,i)   = real(coord(i1,3),wp)
      case (FittingMoveREF)
        dmov(1:2) = restraint_refcoord(1:2,ig) - com_mov(1:2)
        rotated_coord(3,i) = restraint_refcoord(3,ig)
      end select
      rotated_coord(1:2,i) = real(rot_matrix(1:2,1)*dmov(1),wp) &
                           + real(rot_matrix(1:2,2)*dmov(2),wp)
      rotated_coord(1:2,i) = rotated_coord(1:2,start_i+ix) &
                           + real(com_ref(1:2),wp)
    end do

    return

  end subroutine fitting_sel_2d

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    reduce_com
  !> @brief        Reduces COM coordinates, velocities, and total mass over MPI.
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine reduce_com(val1, val2, val3)

  ! formal arguments
    real(dp),                intent(inout) :: val1(:), val2(:), val3

    ! local variables
    real(dp)                      :: before_reduce(7), after_reduce(7)


    before_reduce(1:3) = val1(1:3)
    before_reduce(4:6) = val2(1:3)
    before_reduce(7)   = val3

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(before_reduce, after_reduce, 7, mpi_real8, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    after_reduce(1:7) = before_reduce(1:7)
#endif

    val1(1:3)    = after_reduce(1:3)
    val2(1:3)    = after_reduce(4:6)
    val3         = after_reduce(7)

    return

  end subroutine reduce_com

end module cg_enefunc_fit_mod
