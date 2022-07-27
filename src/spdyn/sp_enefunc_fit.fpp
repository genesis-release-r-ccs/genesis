!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_enefunc_fit_mod
!> @brief   for the use of fit
!! @authors Koichi Tamura (KT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8
  
#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_enefunc_fit_mod

  use sp_restraints_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
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
#ifdef HAVE_MPI_GENESIS
#ifdef MSMPI
!GCC$ ATTRIBUTES DLLIMPORT :: MPI_BOTTOM, MPI_IN_PLACE
#endif
#endif
  private

  ! subroutines
  public  :: setup_fitting_spdyn
  public  :: fitting_sel
  public  :: fitting_sel_2d
  private :: setup_enefunc_fit_domain
  private :: setup_enefunc_fit_refcoord
  private :: reduce_com

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_fitting_spdyn
  !> @brief        setup fitting for spdyn
  !! @authors      CK
  !! @param[in]    is_rpath   : flag if rpath or not
  !! @param[in]    fit_info   : fitting information
  !! @param[in]    sel_info   : selection information
  !! @param[in]    domain     : domain   information
  !! @param[in]    molecule   : molecule information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_fitting_spdyn(is_rpath, fit_info, sel_info,  &
                                 domain, molecule, enefunc)

    ! formal arguments
    logical,                 intent(in)    :: is_rpath
    type(s_fit_info),        intent(in)    :: fit_info
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_domain),          intent(in)    :: domain
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                                :: i
    integer                                :: fitting_atom_idx
    integer                                :: fitting_move
    logical                                :: fit_check, fit_check_rpath
    type(s_selatoms)                       :: selatoms
    

    enefunc%fitting_method  = fit_info%fitting_method

    if (enefunc%pressure_rmsd .and. main_rank) then
      if (enefunc%fitting_method .ne. FittingMethodTR .and.      &
          enefunc%fitting_method .ne. FittingMethodTR_ROT .and.  &
          enefunc%fitting_method .ne. FittingMethodTR_ZROT) then
          write(MsgOut,'(A)') "Setup_Fitting_Spdyn> pressure_rmsd option is set without translate fitting"
      end if
    end if

    if (enefunc%fitting_method == FittingMethodNO) then
      if (main_rank) then
        fit_check = .false.
        do i = 1,enefunc%num_restraintfuncs
          if (enefunc%restraint_kind(i) == RestraintsFuncRMSD .or. &
            enefunc%restraint_kind(i) == RestraintsFuncRMSDCOM) then
            fit_check = .true.
            exit
          end if
        end do
        if (fit_check .and. .not. fit_info%force_no_fitting) then
          call error_msg('Setup_Fitting_Spdyn> No fit is not allowed '//&
                      'in RMSD restraint')
        endif
        if (fit_check .and. fit_info%force_no_fitting) then
          write(MsgOut,'(A)') "Setup_Fitting_Spdyn> RMSD restraint without FITTING"
          write(MsgOut,*) 
        end if
      end if
      return
    end if

    fit_check       = .false.
    do i = 1,enefunc%num_restraintfuncs
      if (enefunc%restraint_kind(i) == RestraintsFuncRMSD .or. &
          enefunc%restraint_kind(i) == RestraintsFuncRMSDCOM .or. &
          enefunc%restraint_kind(i) == RestraintsFuncPC .or. &
          enefunc%restraint_kind(i) == RestraintsFuncPCCOM) then
          fit_check = .true.
          exit
      end if
    end do

    fit_check_rpath = .false.
    do i = 1,enefunc%num_restraintfuncs
      if (enefunc%restraint_kind(i) == RestraintsFuncPOSI .and. &
          is_rpath) then
        fit_check_rpath = .true.
        exit
      end if
    end do

    if (.not. fit_check .and. .not. fit_check_rpath) then
       if (main_rank) then
         write(MsgOut,'(A)') "Setup_Fitting_Spdyn> NO fitting is applied, skip"
         write(MsgOut,'(A)') "  fitting method  =  NO"
         write(MsgOut,*) 
       end if
       enefunc%fitting_method=FittingMethodNO
       return
    end if

    if (enefunc%fitting_method /= FittingMethodTR_ROT .and. &
        enefunc%fitting_method /= FittingMethodXYTR_ZROT) &
      call error_msg('Setup_Fitting_Spdyn> NO/TR+ROT/XYTR+ZROT is allowed')

    fitting_atom_idx      = fit_info%fitting_atom
    enefunc%mass_weight   = fit_info%mass_weight

    if (is_rpath) then
      enefunc%fitting_file  = FittingFileFIT
      enefunc%fitting_move  = FittingMoveSYS
    else
      enefunc%fitting_file  = FittingFileREF
      enefunc%fitting_move  = FittingMoveREF
    end if

    do i = 1,enefunc%num_restraintfuncs
      if (enefunc%restraint_kind(i) == RestraintsFuncRMSDCOM .and. &
        (.not. enefunc%mass_weight) ) &
        call error_msg('Setup_Fitting_Spdyn> RESTRAINTS and FITTING '//&
                       'is inconsistent.  Please check [FITTING] section')
      if (enefunc%restraint_kind(i) == RestraintsFuncRMSD .and. &
        (enefunc%mass_weight) ) &
        call error_msg('Setup_Fitting_Spdyn> RESTRAINTS and FITTING'//&
                       'is inconsistent.  Please check [FITTING] section')
      if (enefunc%restraint_kind(i) == RestraintsFuncRMSD .or. &
          enefunc%restraint_kind(i) == RestraintsFuncRMSDCOM) then
        if (fitting_atom_idx /= enefunc%restraint_grouplist(1,i)) then
          call error_msg('Setup_Fitting_Spdyn> fitting group should be'//&
                       'same with restraint group in RMSD restraint')
        end if
      end if
    end do

    call setup_enefunc_fit_refcoord(enefunc%fitting_file, molecule, enefunc)

    call select_atom(molecule, sel_info%groups(fitting_atom_idx), selatoms)
    call setup_enefunc_fit_domain(domain, enefunc, selatoms%idx)

    return

  end subroutine setup_fitting_spdyn

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_fit_domain
  !> @brief        setup fit information about domain
  !! @authors      KT
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_fit_domain(domain, enefunc, fitting_index)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc
    integer,                 intent(in)    :: fitting_index(:)

    ! local variable
    integer                  :: i, j, k, ix, iatm, icel

    integer,         pointer :: ncell
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: nfitting(:)
    real(wp),        pointer :: fit_refcoord(:,:), fit_coord(:,:,:)


    ncell           => domain%num_cell_local
    id_g2l          => domain%id_g2l

    nfitting        => enefunc%nfitting
    fit_refcoord    => enefunc%fit_refcoord

    call alloc_enefunc(enefunc, EneFuncFitd, ncell)

    fit_coord       => enefunc%fit_coord

    do i = 1, size(fitting_index)

      iatm = fitting_index(i)
      icel = id_g2l(1,iatm)

      if (icel > 0 .and. icel <= ncell) then
        nfitting(icel) = nfitting(icel) + 1
        enefunc%fitting_atom(nfitting(icel),icel) = iatm
        fit_coord(1:3,nfitting(icel),icel) = fit_refcoord(1:3,iatm)
        
        !if (my_country_no + 1 == 7) then 
        !  write(MsgOut,*) " fit "  
        !  write(MsgOut,*) my_country_no + 1, icel, nfitting(icel)
        !  write(MsgOut,*) iatm, fit_refcoord(1:3,iatm)
        !  write(MsgOut,*) " "
        !end if
      end if

    end do

    return

  end subroutine setup_enefunc_fit_domain

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_fit_refcoord
  !> @brief        setup enefunc fit reference coords
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
      call error_msg('Setup_Enefunc_Fit_Refcoord> bad fitfile in [INPUT]')

    else if (fitting_file == FittingFileREF) then
      if (size(enefunc%restraint_refcoord(1,:)) /= molecule%num_atoms) &
      call error_msg('Setup_Enefunc_Fit_Refcoord> bad reffile in [INPUT]')
     
    end if

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

    end if

    return

  end subroutine setup_enefunc_fit_refcoord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    fitting_sel
  !> @brief        do fitting
  !! @authors      KT
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine fitting_sel(domain, enefunc, coord, rotated_coord)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wip),               intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: rotated_coord(:,:,:)

    ! local variables
    integer                  :: i, ix, i1, j, k
    integer                  :: ncell_local
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: nrestraint(:), restraint_atom(:,:)
    integer,         pointer :: nfitting(:), fitting_atom(:,:)
  
    integer                  :: ig, ierr, icel
    integer                  :: ikind
    real(dp)                 :: m, rmsd(1:2)
    real(wip),       pointer :: mass(:,:)
    real(wp),        pointer :: restraint_coord(:,:,:), restraint_refcoord(:,:)
    real(wp),        pointer :: fit_coord(:,:,:)
    real(dp)                 :: com_ref(1:3), com_mov(1:3), tot_mass
    real(dp)                 :: sym_matrix(1:4,1:4)
    real(dp)                 :: rot_matrix(1:3,1:3)
    real(dp)                 :: eval(1:4), evec(1:4), work(1:11)
    real(dp)                 :: dref(1:3), dmov(1:3), dsub(1:3), dadd(1:3)

  
    id_g2l          => domain%id_g2l
    mass            => domain%mass

    nrestraint      => enefunc%num_restraint
    restraint_atom  => enefunc%restraint_atom
    restraint_coord => enefunc%restraint_coord
    restraint_refcoord => enefunc%restraint_refcoord
    
    nfitting        => enefunc%nfitting
    fitting_atom    => enefunc%fitting_atom
    fit_coord       => enefunc%fit_coord

    ncell_local = domain%num_cell_local

    ikind           = enefunc%fitting_move
    
    com_ref(1:3) = 0.0_dp 
    com_mov(1:3) = 0.0_dp 
    tot_mass     = 0.0_dp 
    do i = 1, ncell_local
      do ix = 1, nfitting(i)
        ig = fitting_atom(ix,i)
        i1 = id_g2l(2,ig)
        icel = id_g2l(1,ig)
        !m = mass(i1,i)
        m = 1.0_dp
        if (enefunc%mass_weight) m = mass(i1, i)
        select case (ikind)
        case (FittingMoveSYS)
          com_mov(1:3) = com_mov(1:3) + m * coord(1:3,i1,i)
          com_ref(1:3) = com_ref(1:3) + m * enefunc%fit_refcoord(1:3,ig)
        case (FittingMoveREF) 
          com_mov(1:3) = com_mov(1:3) + m * restraint_refcoord(1:3,ig)
          com_ref(1:3) = com_ref(1:3) + m * coord(1:3,i1,i)
        end select
        tot_mass       = tot_mass     + m
      end do
    end do
    call reduce_com(com_ref, com_mov, tot_mass)
    com_ref(1:3) = com_ref(1:3) / tot_mass
    com_mov(1:3) = com_mov(1:3) / tot_mass

 
    ! symmetric matrix
    !
    sym_matrix(1:4,1:4) = 0.0_dp  

    do i = 1, ncell_local
      do ix = 1, nfitting(i)
        ig = fitting_atom(ix,i)
        i1 = id_g2l(2,ig)
        m = 1.0_dp 
        if (enefunc%mass_weight) m = real(mass(i1, i),dp)
        select case (ikind)
        case (FittingMoveSYS) 
          dmov(1:3) = coord(1:3,i1,i)              - com_mov(1:3)
          dref(1:3) = enefunc%fit_refcoord(1:3,ig) - com_ref(1:3)
        case (FittingMoveREF) 
          dmov(1:3) = restraint_refcoord(1:3,ig)   - com_mov(1:3)
          dref(1:3) = coord(1:3,i1,i)              - com_ref(1:3)
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
    end do

    sym_matrix(2,1) = sym_matrix(1,2)
    sym_matrix(3,1) = sym_matrix(1,3)
    sym_matrix(3,2) = sym_matrix(2,3)
    sym_matrix(4,1) = sym_matrix(1,4)
    sym_matrix(4,2) = sym_matrix(2,4)
    sym_matrix(4,3) = sym_matrix(3,4)
#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, sym_matrix, 16, mpi_real8, mpi_sum, &
                       mpi_comm_country, ierror)
#endif

    !if (my_country_no + 1 == 1) then
      !write(MsgOut,*) "sym(1,1)", sym_matrix(1,1)
      !write(MsgOut,*) "sym(1,2)", sym_matrix(1,2)
      !write(MsgOut,*) "sym(1,3)", sym_matrix(1,3)
      !write(MsgOut,*) "sym(1,4)", sym_matrix(1,4)
      !write(MsgOut,*) "sym(2,2)", sym_matrix(2,2)
      !write(MsgOut,*) "sym(2,3)", sym_matrix(2,3)
      !write(MsgOut,*) "sym(2,4)", sym_matrix(2,4)
      !write(MsgOut,*) "sym(3,3)", sym_matrix(3,3)
      !write(MsgOut,*) "sym(3,4)", sym_matrix(3,4)
      !write(MsgOut,*) "sym(4,4)", sym_matrix(4,4)
    !end if


    ! eigenvalue and eigenvector
    !
#ifdef LAPACK
    call dsyev('V', 'U', 4, sym_matrix, 4, eval, work, 11, ierr)
#else
    call error_msg('Fit_Trrot> ERROR: Rotation for RMSD needs LAPACK.')
#endif

    if (ierr /= 0) then
      ierr = -2
      return
    end if

    ! rotation matrix
    !
    evec(1:4) = sym_matrix(1:4,1)
    
    !if (my_country_no + 1 == 1) then
    !  write(MsgOut,*) "evec(1:4)", evec(1:4)
    !end if

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

    !if (my_country_no + 1 == 1) then
      !write(MsgOut,*) "rot_matrix(1,1)", rot_matrix(1,1)
      !write(MsgOut,*) "rot_matrix(1,2)", rot_matrix(1,2)
      !write(MsgOut,*) "rot_matrix(1,3)", rot_matrix(1,3)
      !write(MsgOut,*) "rot_matrix(2,1)", rot_matrix(2,1)
      !write(MsgOut,*) "rot_matrix(2,2)", rot_matrix(2,2)
      !write(MsgOut,*) "rot_matrix(2,3)", rot_matrix(2,3)
      !write(MsgOut,*) "rot_matrix(3,1)", rot_matrix(3,1)
      !write(MsgOut,*) "rot_matrix(3,2)", rot_matrix(3,2)
      !write(MsgOut,*) "rot_matrix(3,3)", rot_matrix(3,3)
    !end if
  

    do i = 1, ncell_local
      do ix = 1, nrestraint(i)
        ig = restraint_atom(ix,i)
        i1 = id_g2l(2,ig)
        select case (ikind)
        case (FittingMoveSYS) 
          dmov(1:3) = coord(1:3,i1,i) - com_mov(1:3)
        case (FittingMoveREF) 
          dmov(1:3) = restraint_refcoord(1:3,ig) - com_mov(1:3)
        end select
        rotated_coord(1:3,ix,i) = rot_matrix(1:3,1)*dmov(1) + &
                                  rot_matrix(1:3,2)*dmov(2) + &
                                  rot_matrix(1:3,3)*dmov(3)
        rotated_coord(1:3,ix,i) = rotated_coord(1:3,ix,i) + com_ref(1:3)
      end do
    end do

!    rmsd(1:2) = 0.0_wip
!
!    do i = 1, ncell_local
!      do ix = 1, nfitting(i)
!        ig = fitting_atom(ix,i)
!        i1 = id_g2l(2,ig)
!        m = 1.0_wip
!        if (enefunc%mass_weight) m = mass(i1, i)
!        if (ikind == FittingMoveSYS) then
!          dref(1:3) = enefunc%fit_refcoord(1:3,ig)
!          dadd(1:3) = rotated_coord(1:3,i1,i)
!        else if (ikind == FittingMoveREF) then
!          dref(1:3) = restraint_refcoord(1:3,ig)
!          dadd(1:3) = rotated_coord(1:3,ix,i)
!        end if
!        dsub(1:3) = dref(1:3) - dadd(1:3)
!        rmsd(1) = rmsd(1) + m * &
!                        (dsub(1)*dsub(1)+dsub(2)*dsub(2)+dsub(3)*dsub(3))
!        if (ikind == FittingMoveSYS) then
!          dadd(1:3) = coord(1:3,i1,i)
!        else if (ikind == FittingMoveREF) then
!          dadd(1:3) = restraint_refcoord(1:3,ig)
!        end if
!        dsub(1:2) = dref(1:2) - dadd(1:2)
!        rmsd(2) = rmsd(2) + m * &
!                        (dsub(1)*dsub(1)+dsub(2)*dsub(2)+dsub(3)*dsub(3))
!      end do
!    end do
!    call mpi_allreduce(mpi_in_place, rmsd, 2, mpi_wip_real, mpi_sum, &
!                         mpi_comm_country, ierror)
!    rmsd = sqrt(rmsd/tot_mass)
!!    if (my_country_no + 1 == 7) then
!    if (replica_main_rank) then
!!    if (rmsd > 1.0_wip) then
!     write(MsgOut,*) "aftres:",my_country_no + 1, rmsd(1:2)
!     write(MsgOut,*) " "
!    end if

    return

  end subroutine fitting_sel

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    fitting_sel_2d
  !> @brief        do fitting in XY space
  !! @authors      KT
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine fitting_sel_2d(domain, enefunc, coord, rotated_coord)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wip),               intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: rotated_coord(:,:,:)

    ! local variables
    integer                  :: i, ix, i1, j, k
    integer                  :: ncell_local
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: nrestraint(:), restraint_atom(:,:)
    integer,         pointer :: nfitting(:), fitting_atom(:,:)
  
    integer                  :: ig, ierr, icel
    integer                  :: ikind
    real(dp)                 :: m, rmsd(1:2)
    real(wip),       pointer :: mass(:,:)
    real(wp),        pointer :: restraint_coord(:,:,:), restraint_refcoord(:,:)
    real(wp),        pointer :: fit_coord(:,:,:)
    real(dp)                 :: com_ref(1:3), com_mov(1:3), tot_mass
    real(dp)                 :: detu, detvt, sig
    real(dp)                 :: cor_matrix(1:2,1:2)
    real(dp)                 :: rot_matrix(1:2,1:2)
    real(dp)                 :: s(1:2), u(1:2,1:2), vt(1:2,1:2), work(1:10)
    real(dp)                 :: dref(1:3), dmov(1:3), dsub(1:3), dadd(1:3)

  
    id_g2l          => domain%id_g2l
    mass            => domain%mass

    nrestraint      => enefunc%num_restraint
    restraint_atom  => enefunc%restraint_atom
    restraint_coord => enefunc%restraint_coord
    restraint_refcoord => enefunc%restraint_refcoord
    
    nfitting        => enefunc%nfitting
    fitting_atom    => enefunc%fitting_atom
    fit_coord       => enefunc%fit_coord

    ncell_local = domain%num_cell_local

    ikind           = enefunc%fitting_move
    
    com_ref(1:3) = 0.0_dp 
    com_mov(1:3) = 0.0_dp 
    tot_mass     = 0.0_dp 
    do i = 1, ncell_local
      do ix = 1, nfitting(i)
        ig = fitting_atom(ix,i)
        i1 = id_g2l(2,ig)
        icel = id_g2l(1,ig)
        m = 1.0_dp 
        if (enefunc%mass_weight) m = mass(i1, i)
        select case (ikind)
        case (FittingMoveSYS)
          com_mov(1:3) = com_mov(1:3) + m * coord(1:3,i1,i)
          com_ref(1:3) = com_ref(1:3) + m * enefunc%fit_refcoord(1:3,ig)
        case (FittingMoveREF) 
          com_mov(1:3) = com_mov(1:3) + m * restraint_refcoord(1:3,ig)
          com_ref(1:3) = com_ref(1:3) + m * coord(1:3,i1,i)
        end select
        tot_mass       = tot_mass     + m
      end do
    end do
    call reduce_com(com_ref, com_mov, tot_mass)
    com_ref(1:3) = com_ref(1:3) / tot_mass
    com_mov(1:3) = com_mov(1:3) / tot_mass

 
    ! correlation matrix
    !
    cor_matrix(1:2,1:2) = 0.0_dp 

    do i = 1, ncell_local
      do ix = 1, nfitting(i)
        ig = fitting_atom(ix,i)
        i1 = id_g2l(2,ig)
        m = 1.0_wip
        select case (ikind)
        case (FittingMoveSYS)
          dmov(1:2) = coord(1:2,i1,i)              - com_mov(1:2)
          dref(1:2) = enefunc%fit_refcoord(1:2,ig) - com_ref(1:2)
        case (FittingMoveREF)
          dmov(1:2) = restraint_refcoord(1:2,ig)   - com_mov(1:2)
          dref(1:2) = coord(1:2,i1,i)              - com_ref(1:2)
        end select
        cor_matrix(1,1) = cor_matrix(1,1) + dmov(1) * dref(1)
        cor_matrix(1,2) = cor_matrix(1,2) + dmov(2) * dref(1)
        cor_matrix(2,1) = cor_matrix(2,1) + dmov(1) * dref(2)
        cor_matrix(2,2) = cor_matrix(2,2) + dmov(2) * dref(2)
      end do
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, cor_matrix, 4, mpi_real8, mpi_sum, &
                       mpi_comm_country, ierror)
#endif

    ! eigenvalue and eigenvector
    !
#ifdef LAPACK
    call dgesvd('A', 'A', 2, 2, cor_matrix, 2, s, u, 2, vt, 2, work, 10, ierr)
#else
    call error_msg('Fit_Trrot> ERROR: Rotation for RMSD needs LAPACK.')
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
    if (sig < 0.0_wip) then
      s(2)   = -s(2)
      u(:,2) = -u(:,2)
    end if


    ! rotation matrix
    !
    rot_matrix = matmul(u,vt)
    
    do i = 1, ncell_local
      do ix = 1, nrestraint(i)
        ig = restraint_atom(ix,i)
        i1 = id_g2l(2,ig)
        select case (ikind)
        case (FittingMoveSYS)
          dmov(1:2) = coord(1:2,i1,i) - com_mov(1:2)
          rotated_coord(3,ix,i)   = coord(3,i1,i) 
        case (FittingMoveREF)
          dmov(1:2) = restraint_refcoord(1:2,ig) - com_mov(1:2)
          rotated_coord(3,ix,i)   = restraint_refcoord(3,ig)
        end select
        rotated_coord(1:2,ix,i) = rot_matrix(1:2,1)*dmov(1) + &
                                  rot_matrix(1:2,2)*dmov(2)
        rotated_coord(1:2,ix,i) = rotated_coord(1:2,ix,i) + com_ref(1:2)
      end do
    end do


!    rmsd(1:2) = 0.0_wip
!
!    do i = 1, ncell_local
!      do ix = 1, nfitting(i)
!        ig = fitting_atom(ix,i)
!        i1 = id_g2l(2,ig)
!        dref(1:2) = enefunc%fit_refcoord(1:2,ig)
!        if (ikind == FittingMoveSYS) then
!          dref(1:2) = enefunc%fit_refcoord(1:2,ig)
!          dadd(1:2) = rotated_coord(1:2,i1,i)
!        else if (ikind == FittingMoveREF) then
!          dref(1:2) = restraint_refcoord(1:2,ig)
!          dadd(1:2) = rotated_coord(1:2,ix,i)
!        end if
!        m = 1.0_wip
!        if (enefunc%mass_weight) m = mass(i1, i)
!        dsub(1:2) = dref(1:2) - dadd(1:2)
!        rmsd(1) = rmsd(1) + m * &
!                        (dsub(1)*dsub(1)+dsub(2)*dsub(2))
!        if (ikind == FittingMoveSYS) then
!          dadd(1:2) = coord(1:2,i1,i)
!        else if (ikind == FittingMoveREF) then
!          dadd(1:2) = restraint_refcoord(1:2,ig)
!        end if
!        dsub(1:2) = dref(1:2) - dadd(1:2)
!        rmsd(2) = rmsd(2) + m * &
!                        (dsub(1)*dsub(1)+dsub(2)*dsub(2))
!      end do
!    end do
!    call mpi_allreduce(mpi_in_place, rmsd, 2, mpi_wip_real, mpi_sum, &
!                         mpi_comm_country, ierror)
!    rmsd = sqrt(rmsd/tot_mass)
!!    if (my_country_no + 1 == 7) then
!    if (replica_main_rank) then
!!    if (rmsd > 1.0_wip) then
!     write(MsgOut,*) "aftres:",my_country_no + 1, rmsd(1:2)
!     write(MsgOut,*) " "
!    end if

    return

  end subroutine fitting_sel_2d

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    reduce_com
  !> @brief
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine reduce_com(val1, val2, val3)

  ! formal arguments
    real(dp),                intent(inout) :: val1(:), val2(:), val3

    ! local variables
    real(dp)                 :: before_reduce(7), after_reduce(7)


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

end module sp_enefunc_fit_mod
